module clump_finder_module
  use clump_merger_module
contains
subroutine m_clump_finder(pst,create_output,keep_alive)
  use output_clump_module
  use output_part_module, only: r_output_part
  use output_amr_module, only: output_header
  use amr_parameters, only: flen
  use mdl_module, only: mdl_mkdir, mdl_wtime
  use ramses_commons, only: pst_t
#ifdef GRAV
  use rho_fine_module, only: m_rho_fine
#endif
#ifdef _CUDA
  use gpu_manager, only: r_transfer_grid_host
#endif
  implicit none
  type(pst_t)::pst
  logical::create_output,keep_alive
  !-----------------------------------------------------------------------
  ! This is the master routine for the RAMSES clump finder.
  !-----------------------------------------------------------------------
  character(LEN=5)::nchar
  character(LEN=flen)::filename,filedir
  integer,dimension(1:flen/4)::input_array
  double precision::ttend, ttstart=0.0
  integer::dummy(1)
  integer::no_halo

#if NDIM==3 && defined(GRAV)

  associate(r=>pst%s%r,g=>pst%s%g,mdl=>pst%s%mdl,p=>pst%s%p,star=>pst%s%star,sink=>pst%s%sink,tree=>pst%s%tree)

  write(*,*)'Entering clump finder'
  ttstart = mdl_wtime(mdl)

  !--------------------------------------------------------------
  ! Compute rho from dark matter or stars or sinks or gas density
  !--------------------------------------------------------------
  call m_rho_fine(pst,r%levelmin,r%rho_type_clump)

  !-------------------------------------
  ! Find relevant peak patches and halos
  !-------------------------------------
  call r_clump_finder(pst)

  !---------------------------------
  ! Output clumps properties to file
  !---------------------------------
  if(create_output.and.pst%s%c%npeak_tot>0)then
     call m_timer('clump - output','start')
     call title(g%ifout,nchar)
     filedir='output_'//TRIM(nchar)//'/'
     if(r%clump_only)then
        call title(g%ifout-1,nchar)
        filedir='catalog_output_'//TRIM(nchar)//'/'
     endif
     call mdl_mkdir(mdl,filedir)
     input_array=transfer(filedir,input_array)
     if(r%output_clump)then
        write(*,*)'Writing clump files'
     endif
     if(r%output_peak_grid)then
        write(*,*)'Writing peak grid files'
        filename=TRIM(filedir)//'peak_grid_header.txt'
        call file_descriptor_clump(r,filename)
     endif
     if(r%output_peak_part.and.r%part)then
        write(*,*)'Writing peak part files'
        filename=TRIM(filedir)//'peak_part_header.txt'
        call output_peak_header(r,g,p   ,filename)
     endif
     if(r%output_peak_star.and.r%star)then
        write(*,*)'Writing peak star files'
        filename=TRIM(filedir)//'peak_star_header.txt'
        call output_peak_header(r,g,star,filename)
     endif
     if(r%output_peak_sink.and.r%sink)then
        write(*,*)'Writing peak sink files'
        filename=TRIM(filedir)//'peak_sink_header.txt'
        call output_peak_header(r,g,sink,filename)
     endif
     if(r%output_peak_tree.and.r%tree)then
        write(*,*)'Writing peak tree files'
        filename=TRIM(filedir)//'peak_tree_header.txt'
        call output_peak_header(r,g,tree,filename)
     endif
     if(r%output_peak_core.and.r%core)then
        write(*,*)'Writing peak core files'
        filename=TRIM(filedir)//'peak_core_header.txt'
        call output_peak_header(r,g,pst%s%core,filename)
     endif
     call r_output_clump(pst,input_array,flen/4,dummy,0)

     ! In clump_only mode, also write the particles in the order matching the peak files
     if(r%clump_only.and.r%pic)then
#ifdef _CUDA
        call r_transfer_grid_host(pst)
#endif
        if(r%part)then
           filename=TRIM(filedir)//'part_header.txt'
           call output_header(r,g,p   ,filename)
        endif
        if(r%star)then
           filename=TRIM(filedir)//'star_header.txt'
           call output_header(r,g,star,filename)
        endif
        if(r%sink)then
           filename=TRIM(filedir)//'sink_header.txt'
           call output_header(r,g,sink,filename)
        endif
        if(r%tree)then
           filename=TRIM(filedir)//'tree_header.txt'
           call output_header(r,g,tree,filename)
        endif
        if(r%trac)then
           filename=TRIM(filedir)//'trac_header.txt'
           call output_header(r,g,pst%s%trac,filename)
        endif
        if(r%dust)then
           filename=TRIM(filedir)//'dust_header.txt'
           call output_header(r,g,pst%s%dust,filename)
        endif
        if(r%core)then
           filename=TRIM(filedir)//'core_header.txt'
           call output_header(r,g,pst%s%core,filename)
        endif
        call r_output_part(pst,input_array,flen/4,dummy,0)
     endif
     call m_timer('clump','start')
  endif

  !------------------------------------------
  ! Deallocate all peak arrays if needed
  !------------------------------------------
  if(.not. keep_alive)then
     call r_deallocate_clump(pst)
  endif

  ttend = mdl_wtime(mdl)
  print '(A,F0.7)',' Time elapsed in finding clumps: ',ttend-ttstart

  end associate

#endif
end subroutine m_clump_finder
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_clump_finder(pst)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
#ifdef _CUDA
  use clump_device, only: gpu_clump_finder
#endif
  implicit none
  type(pst_t)::pst

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_CLUMP_FINDER,pst%iUpper+1)
     call r_clump_finder(pst%pLower)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
#ifdef _CUDA
     call gpu_clump_finder(pst%s)
#else
     call clump_finder(pst%s)
#endif
  endif

end subroutine r_clump_finder
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine clump_finder(s)
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t)::s

#if NDIM==3 && defined(GRAV)

  !-------------------------------------------------
  ! Store clump finder parameters in clump object.
  !-------------------------------------------------
  s%c%relevance_threshold = s%r%relevance_threshold
  s%c%density_threshold = s%r%density_threshold
  s%c%saddle_threshold = s%r%saddle_threshold
  s%c%mass_threshold = s%r%mass_threshold
  s%c%purity_threshold = s%r%purity_threshold
  s%c%fraction_threshold = s%r%fraction_threshold
  !----------------------------------------------------------------------
  ! Count and collect all cells above the prescribed density threshold.
  ! We call these cell test particles for the watershed algorithm.
  !----------------------------------------------------------------------
  call collect_test(s)
  if(s%c%ntest_tot==0)return
  !----------------------------------------------------------------------
  ! Count and collect all density peaks.
  ! We also compute for each test particle the coordinates of its
  ! densest neighbor.
  !----------------------------------------------------------------------
  call collect_peak(s)
  if(s%c%npeak_tot==0)return
  !----------------------------------------------------------------------
  ! Perform a segmentation of the density field using the watershed
  ! algorithm. We get well defined peak-patches around each peak.
  ! As a result, each pair of neighboring peak-patches are separated
  ! by their saddle surface.
  !----------------------------------------------------------------------
  call collect_patch(s)
  !----------------------------------------------------------------------
  ! Allocate all peak-patch based arrays.
  !----------------------------------------------------------------------
  call allocate_peak_patch_arrays(s)
  !----------------------------------------------------------------------
  ! Compute the maximum density saddle point and its corresponding
  ! unique neighboring peak.
  !----------------------------------------------------------------------
  call collect_saddle(s)
  !----------------------------------------------------------------------
  ! Merge peak-patches based on a relevance criterion.
  ! Peak-patches that are due to random noise fluctuations are merged
  ! into relevant peak-patches.
  !----------------------------------------------------------------------
  call merge_clumps(s,'relevance')
  !----------------------------------------------------------------------
  ! Compute peak-patches properties such as mass, number of cells...
  ! Global peak-patch ids are stored into flag2 for possible later use.
  !----------------------------------------------------------------------
  call compute_clump_properties(s,s%r%rho_type_clump)
  !----------------------------------------------------------------------
  ! Compute particle peak ids and store them in array pid.
  ! Compute particle-based peak-patch properties.
  !----------------------------------------------------------------------
  if(s%r%part)then
     call particle_peak_id(s,s%p)
     if(s%r%rho_type_clump.eq.1)then
        call particle_clump_properties(s,s%p)
     endif
  endif
  if(s%r%star)then
     call particle_peak_id(s,s%star)
     if(s%r%rho_type_clump.eq.2)then
        call particle_clump_properties(s,s%star)
     endif
  endif
  if(s%r%sink)then
     call particle_peak_id(s,s%sink)
  endif
  if(s%r%core)then
     call particle_peak_id(s,s%core)
  endif
  if(s%r%tree)then
     call particle_peak_id(s,s%tree)
  endif
  !----------------------------------------------------------------------
  ! Merge all neighboring peak-patches above the prescribed density
  ! threshold into halo-patches, if their saddle point density is larger
  ! that the prescribed saddle density threshold.
  ! Global halo-patch ids are stored into flag1 for possible later use.
  !----------------------------------------------------------------------
  if(s%c%saddle_threshold>0)then
     call merge_clumps(s,'saddleden')
  endif
  !----------------------------------------------------------------------
  ! Compute gravitational self-potential of all clumps hierarchically.
  ! Particle peak ids need to be reset after that.
  !----------------------------------------------------------------------
  if(s%c%saddle_threshold>0)then
     if(s%r%part.and.s%r%rho_type_clump.eq.1)then
        call particle_potential(s,s%p)
        call particle_peak_id(s,s%p)
     endif
  endif
  !----------------------------------------------------------------------
  ! Unbind particles from clumps hierarchically.
  ! The particle peak ids are updated accordingly.
  ! The halo patch is used as garbage collector.
  !----------------------------------------------------------------------
  if(s%c%saddle_threshold>0)then
     if(s%r%part.and.s%r%rho_type_clump.eq.1)then
        call particle_unbind(s,s%p)
     endif
  endif
  !----------------------------------------------------------------------
  ! Compute final clump mass profiles hierarchically.
  !----------------------------------------------------------------------
  if(s%c%saddle_threshold>0)then
     if(s%r%part.and.s%r%rho_type_clump.eq.1)then
        call mass_profile(s,s%p)
     endif
  endif
  !----------------------------------------------------------------------
  ! Remove all peak-patches (resp. halo-patches) that are below
  ! the relevance threshold or the mass threshold by setting their
  ! flag2 (resp. flag1) field values to zero.
  !----------------------------------------------------------------------
  call trim_clumps(s)
  !----------------------------------------------------------------------
  ! Compute particle halo ids and store them in array hid.
  !----------------------------------------------------------------------
  if(s%c%saddle_threshold>0)then
     if(s%r%part)then
        call particle_halo_id(s,s%p)
     endif
     if(s%r%star)then
        call particle_halo_id(s,s%star)
     endif
     if(s%r%sink)then
        call particle_halo_id(s,s%sink)
     endif
     if(s%r%core)then
        call particle_halo_id(s,s%core)
     endif
     if(s%r%tree)then
        call particle_halo_id(s,s%tree)
     endif
  else
     ! Without halo merging, each clump is its own halo
     if(s%r%part)s%p%hid=s%p%pid
     if(s%r%star)s%star%hid=s%star%pid
     if(s%r%sink)s%sink%hid=s%sink%pid
     if(s%r%core)s%core%hid=s%core%pid
     if(s%r%tree)s%tree%hid=s%tree%pid
  endif
#endif
end subroutine clump_finder
!################################################################
!################################################################
!################################################################
!################################################################
#if NDIM==3 && defined(GRAV)
subroutine collect_test(s)
  use amr_parameters, only: twotondim, ndim
  use ramses_commons, only: ramses_t
  use multigrid_fine_coarse, only: pack_fetch_rho, unpack_fetch_rho
  use nbors_utils
  use cache_commons
  use cache
  use mdl_module, only: mdl_abort
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(ramses_t)::s
  !------------------------------------------------------------------
  ! This is the clump finder routine for collecting test particles
  ! also known as all cells above the prescribed density threshold.
  ! Count number of test particles and share info across processors
  ! Written by Ziyong Wu (mini-ramses version December 2023).
  !------------------------------------------------------------------
#ifndef WITHOUTMPI
  integer::info
  integer,dimension(1:s%g%ncpu)::ntest_cpu_all
#endif
  integer::ind,igrid,idim,icpu,ngrid,nleaf,now_level,next_level
  integer::istep,ipart,jpart,ip,itest
  integer::ilevel,i,levelmin_part
  integer,dimension(1:s%g%ncpu)::ntest_cpu
  integer(kind=8),dimension(0:s%g%ncpu)::ntest_cum
  integer(kind=8)::ntest8
  logical::verbose_all=.false.
  integer::action,icellp,igridp
  logical::ok
  real(kind=8)::dx,vol
  real(kind=8)::d,dx_loc
  real(kind=8),allocatable,dimension(:)::dens
  integer,allocatable,dimension(:)::isort
  integer,allocatable,dimension(:)::iswap
  integer(kind=8),dimension(0:ndim)::hash_key
  type(msg_small_realdp)::dummy_small_realdp

#ifdef GRAV

  associate(r=>s%r,g=>s%g,m=>s%m,c=>s%c,mdl=>s%mdl)

  if(g%myid==1.and.r%verbose)write(*,*)'Entering collect_test'

  !---------------------------------------------------------
  ! Make density field smooth, positive and monotonous.
  !---------------------------------------------------------
  do igrid=m%head(r%levelmin),m%tail(r%levelmin) ! Loop over grids
     do ind=1,twotondim ! Loop over cells
        m%nref(ind,igrid) = m%rho(ind,igrid) ! Save for true mass later
     end do
  end do
  do ilevel=r%levelmin+1,r%nlevelmax ! Loop over levels
     call open_cache(mdl, m, pack_size=storage_size(dummy_small_realdp)/32, &
          pack=pack_fetch_rho, unpack=unpack_fetch_rho)
     hash_key(0)=ilevel
     do igrid=m%head(ilevel),m%tail(ilevel) ! Loop over grids
        hash_key(1:ndim)=m%grid(igrid)%ckey(1:ndim) ! Get parent cell and grid index
        call get_parent_cell(s,hash_key,igridp,icellp,flush_cache=.false.,fetch_cache=.true.)
        do ind=1,twotondim ! Loop over cells
           m%nref(ind,igrid) = m%rho(ind,igrid) ! Save for true mass later
           m%rho(ind,igrid) = 0.5d0*m%rho(ind,igrid) + 0.5d0*m%rho(icellp,igridp)
        end do
     end do
     call close_cache(mdl)
  end do

  !---------------------------------------------------------
  ! Count cells above threshold. We name them test particles
  !---------------------------------------------------------
  ntest8 = 0
  do ilevel=r%levelmin,r%nlevelmax ! Loop over levels
     do igrid=m%head(ilevel),m%tail(ilevel) ! Loop over grids
        do ind=1,twotondim ! Loop over cells
           ok = .not. m%grid(igrid)%refined(ind) ! Select leaf cells
           d = m%rho(ind,igrid)
           ok = ok .and. d > c%density_threshold
           if(ok)then
              ntest8=ntest8+1
           endif
        end do
     end do
  end do
  if(ntest8 > m%ngridmax)then
     write(*,*)'collect_test: cells above density threshold (',ntest8,') exceed ngridmax (',m%ngridmax,')'
     write(*,*)'Increase density_threshold (clump finder, merger tree) or sink_density_threshold (sink formation)'
     call mdl_abort(mdl)
  endif
  c%ntest = int(ntest8, kind=4)

  !-------------------------------------------------
  ! Compute number of test particles across all CPUs
  !-------------------------------------------------
  ntest_cpu=0
  ntest_cpu(g%myid)=c%ntest

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(ntest_cpu,ntest_cpu_all,g%ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  ntest_cpu=ntest_cpu_all
#endif
  ntest_cum=0
  do icpu=1,g%ncpu
     ntest_cum(icpu)=ntest_cum(icpu-1)+int(ntest_cpu(icpu),kind=8)
  end do
  c%ntest_tot=ntest_cum(g%ncpu)
  if(g%myid==1)then
     write(*,'(" Found",I12,"  cells above threshold")')c%ntest_tot
  endif

  if (c%ntest>0) then

     !------------------------------------------
     ! Allocate arrays for cells above threshold
     !------------------------------------------
     allocate(c%grid(c%ntest),c%cell(c%ntest),c%level(c%ntest),c%hash(1:c%ntest,0:ndim))
     c%grid=0;  c%cell=0; c%level=0; c%hash=0

     !---------------------------------------------------
     ! Compute and store arrays for cells above threshold
     !---------------------------------------------------
     allocate(dens(c%ntest))
     itest=0
     do ilevel=r%levelmin,r%nlevelmax ! Loop over levels
        do igrid=m%head(ilevel),m%tail(ilevel) ! Loop over grids
           do ind=1,twotondim ! Loop over cells
              ok=.not.m%grid(igrid)%refined(ind) ! Select leaf cells
              d=m%rho(ind,igrid)
              ok=ok.and.d>c%density_threshold
              if(ok)then
                 itest=itest+1
                 dens(itest)=d
                 c%grid(itest)=igrid
                 c%cell(itest)=ind
                 c%level(itest)=ilevel
              endif
           end do
        end do
     end do

     !-----------------------------------------------------------------------
     ! Sort cells above threshold according to their density
     !-----------------------------------------------------------------------
     allocate(isort(c%ntest))
     do i=1,c%ntest
        dens(i)=-dens(i)
        isort(i)=i
     end do
     call quick_sort_dp(dens,isort,c%ntest)
     deallocate(dens)
     ! Swap arrays to sort them permanently
     allocate(iswap(1:c%ntest))
     do i=1,c%ntest
        itest=isort(i)
        iswap(i)=c%grid(itest)
     end do
     c%grid=iswap
     do i=1,c%ntest
        itest=isort(i)
        iswap(i)=c%cell(itest)
     end do
     c%cell=iswap
     do i=1,c%ntest
        itest=isort(i)
        iswap(i)=c%level(itest)
     end do
     c%level=iswap
     deallocate(iswap,isort)

  endif

  end associate
#endif

end subroutine collect_test
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine collect_peak(s)
  use amr_parameters, only: twotondim, ndim
  use ramses_commons, only: ramses_t
  use cache_commons
  use cache
  use nbors_utils
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(ramses_t)::s
  !------------------------------------------------------------------
  ! This is the clump finder routine for collecting densest
  ! neighbors. Only scan cells above the density threshold.
  ! Density peaks are cells without densest neighbors.
  ! Count number of density peaks and share info across processors.
  ! Store the hash key of the densest neighbor for later usage.
  ! Written by Ziyong Wu (mini-ramses version December 2023).
  !------------------------------------------------------------------
#ifndef WITHOUTMPI
  integer::info
  integer,dimension(1:s%g%ncpu)::npeak_cpu_all
#endif
  type(msg_int4_small_realdp)::dummy_int4_small_realdp
  integer::igridn
  integer::ilevel
  integer::icpu,next_level,now_level,icelln,idim,igrid,ind,itest,j,k
  integer::ipart,jpart,ip,i,icellp,icellpm,ipeak
  integer,dimension(1:s%g%ncpu)::npeak_cpu
  integer,dimension(1:ndim)::ckey,ckey_nbor
  integer(kind=8),dimension(0:ndim)::hash_cell,hash_nbor
  real(kind=8),dimension(1:ndim)::xcen,xnei
  integer,parameter::nSnei=56
  integer,dimension(1:nSnei)::igrid_nbor
  integer(kind=8),dimension(1:nSnei)::icell_nbor,level_nbor
  real(kind=8),dimension(1:ndim,1:nSnei)::xSnei
  real(kind=8)::dens_nbor,density_max,x,y,z
  logical::ok,ok_peak,inside

  associate(r=>s%r,g=>s%g,m=>s%m,c=>s%c,mdl=>s%mdl)

  if(g%myid==1.and.r%verbose)write(*,*)'Entering collect_peak'

  !--------------------------------------------------------
  ! Arrays to define neighbors (center=[0,0,0])
  ! normalized to dx = 1 = size of the central leaf cell
  ! from -0.75 to 0.75
  !--------------------------------------------------------
  ind=0
  do k=1,4
     do j=1,4
        do i=1,4
           ok=.true.
           if((i==2.or.i==3).and.(j==2.or.j==3).and.(k==2.or.k==3)) ok=.false. ! centre
           if(ok)then
              ind = ind+1
              x = (i-1)+0.5d0 - 2
              y = (j-1)+0.5d0 - 2
              z = (k-1)+0.5d0 - 2
              xSnei(1,ind) = x/2d0
              xSnei(2,ind) = y/2d0
              xSnei(3,ind) = z/2d0
           endif
        enddo
     enddo
  enddo

  !----------------------------------------
  ! Compute hash key of densest neighbor
  !----------------------------------------
  call open_cache(mdl, m, pack_size=storage_size(dummy_int4_small_realdp)/32, &
       pack=pack_fetch_saddle, unpack=unpack_fetch_saddle)

  c%npeak=0
  do itest=1,c%ntest
     ilevel=c%level(itest)
     igrid=c%grid(itest)
     ind=c%cell(itest)

     ! Set grid index to null
     icelln=0
     igridn=0
     do j=1,nSnei
        igrid_nbor(j)=0
     end do

     xcen(1)=2*m%grid(igrid)%ckey(1)+MOD((ind-1)  ,2)+0.5
#if NDIM>1
     xcen(2)=2*m%grid(igrid)%ckey(2)+MOD((ind-1)/2,2)+0.5
#endif
#if NDIM>2
     xcen(3)=2*m%grid(igrid)%ckey(3)+MOD((ind-1)/4,2)+0.5
#endif
     ! Collect all neighboring cell from hash table
     do j=1,nSnei

        ! Compute neighboring cell coordinates
        xnei(1:ndim)=xcen(1:ndim)+xSnei(1:ndim,j)
        ! Periodic boundary conditions
        inside=.true.
        do idim=1,ndim
           if(r%periodic(idim))then
              if(xnei(idim)< m%box_ckey_min(idim,ilevel+1))xnei(idim)=xnei(idim)-m%box_ckey_min(idim,ilevel+1)+m%box_ckey_max(idim,ilevel+1)
              if(xnei(idim)>=m%box_ckey_max(idim,ilevel+1))xnei(idim)=xnei(idim)+m%box_ckey_min(idim,ilevel+1)-m%box_ckey_max(idim,ilevel+1)
           endif
           inside=inside.and.xnei(idim)>=m%box_ckey_min(idim,ilevel+1).and.xnei(idim)<m%box_ckey_max(idim,ilevel+1)
        end do
        ! Skip neighbors outside a non-periodic box
        if(.not.inside)cycle

        ! Get neighboring cell at ilevel
        ckey_nbor(1:ndim)=int(xnei(1:ndim))
        hash_nbor(0)=ilevel+1
        hash_nbor(1:ndim)=ckey_nbor(1:ndim)
        call get_parent_cell(s,hash_nbor,igridn,icelln,flush_cache=.false.,fetch_cache=.true.)

        ! If missing, get neighboring cell at ilevel-1
        if(igridn==0)then
           ckey_nbor(1:ndim)=int(xnei(1:ndim)/2.0)
           hash_nbor(0)=ilevel
           hash_nbor(1:ndim)=ckey_nbor(1:ndim)
           call get_parent_cell(s,hash_nbor,igridn,icelln,flush_cache=.false.,fetch_cache=.true.)

           ! If refined, get neighboring cell at ilevel+1
        else if (m%grid(igridn)%refined(icelln))then
           ckey_nbor(1:ndim)=int(xnei(1:ndim)*2.0)
           hash_nbor(0)=ilevel+2
           hash_nbor(1:ndim)=ckey_nbor(1:ndim)
           call get_parent_cell(s,hash_nbor,igridn,icelln,flush_cache=.false.,fetch_cache=.true.)
        endif

        ! Lock grid in cache
        call lock_cache(m,igridn)

        igrid_nbor(j) = igridn
        icell_nbor(j) = icelln
        level_nbor(j) = hash_nbor(0)

     end do

     density_max=1.0001d0*m%rho(ind,igrid)
     ok_peak=.true.

     do j=1,nSnei
        igridn = igrid_nbor(j) ! Gather neighboring grid
        icelln = icell_nbor(j)
        if(igridn>0)then
           dens_nbor = m%rho(icelln,igridn)
           if(dens_nbor > density_max)then
              ok_peak=.false.
              density_max=dens_nbor
              ! Store hash key of densest neighbor
              c%hash(itest,0)=level_nbor(j)
              c%hash(itest,1)=2*m%grid(igridn)%ckey(1)+MOD((icelln-1)  ,2)
#if NDIM>1
              c%hash(itest,2)=2*m%grid(igridn)%ckey(2)+MOD((icelln-1)/2,2)
#endif
#if NDIM>2
              c%hash(itest,3)=2*m%grid(igridn)%ckey(3)+MOD((icelln-1)/4,2)
#endif
           endif
        endif
     end do

     if(ok_peak)then
        c%npeak=c%npeak+1
        c%hash(itest,0:ndim)=0
     endif

     ! Unlock neighboring grids
     do j=1,nSnei
        igridn = igrid_nbor(j)
        call unlock_cache(m,igridn)
     end do

  end do

  call close_cache(mdl)

  !------------------------------------------------
  ! Compute total number of peaks across all CPUs
  ! Determine offset for global peak IDs
  !------------------------------------------------
  allocate(c%npeak_cum(0:g%ncpu))
  npeak_cpu=0
  npeak_cpu(g%myid)=c%npeak
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(npeak_cpu,npeak_cpu_all,g%ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  npeak_cpu=npeak_cpu_all
#endif
  c%npeak_cum=0
  do icpu=1,g%ncpu
     c%npeak_cum(icpu)=c%npeak_cum(icpu-1)+int(npeak_cpu(icpu),kind=8)
  end do
#ifdef WITHOUTMPI
  c%npeak_tot=int(c%npeak,kind=8)
#else
  c%npeak_tot=c%npeak_cum(g%ncpu)
#endif
  if (g%myid==1)then
     write(*,'(" Found",I12,"  density peaks")')int(c%npeak_tot,kind=4)
  endif
  end associate

end subroutine collect_peak
!################################################################
!################################################################
!################################################################
!################################################################
subroutine collect_patch(s)
  use amr_parameters, only: twotondim,ndim
  use ramses_commons, only: ramses_t
  use cache_commons
  use cache
  use nbors_utils
  use boundaries, only: init_bound_flag
  use marshal, only: pack_fetch_flag, unpack_fetch_flag
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(ramses_t)::s
  !------------------------------------------------------------------
  ! This is the clump finder routine for segmenting the density
  ! field into peak patches around each density peak.
  ! - loop over cells in descending density order
  ! - propagate peak id from densest neighbor
  ! - nmove is the number of peak id's passed along
  ! - done when nmove_tot=0
  ! For single core, only one sweep is necessary.
  ! Written by Ziyong Wu (mini-ramses version December 2023).
  !------------------------------------------------------------------
#ifndef WITHOUTMPI
  integer::info
  integer(kind=8)::nmove_all,nzero_all
#endif
  integer(kind=8),dimension(0:ndim)::hash_nbor
  type(msg_int4)::dummy_int4
  integer::igridn,icelln,igrid,ind,ipeak,istep,itest,nmove,nzero
  integer(kind=8)::nmove_tot,nzero_tot

  associate(r=>s%r,g=>s%g,m=>s%m,c=>s%c,mdl=>s%mdl)

  if(g%myid==1.and.r%verbose)write(*,*)'Entering collect_patch'

  !---------------------------------
  ! Allocate first peak based arrays
  !---------------------------------
  allocate(c%peak_cell(c%npeak+c%ncachemax))
  allocate(c%peak_grid(c%npeak+c%ncachemax))
  allocate(c%peak_level(c%npeak+c%ncachemax))
  allocate(c%max_dens(c%npeak+c%ncachemax))
  c%max_dens=0d0; c%peak_cell=0; c%peak_grid=0; c%peak_level=0

  !-------------------------------------------------
  ! Flag peaks with global peak id using flag1 array
  !-------------------------------------------------
  do igrid=1,m%ngridmax
     m%flag1(1:twotondim,igrid)=0
     m%flag2(1:twotondim,igrid)=0
  end do
  ipeak = 0
  do itest=1,c%ntest
     if(c%hash(itest,0)==0)then
        ipeak=ipeak+1
        igrid=c%grid(itest)
        ind=c%cell(itest)
        m%flag1(ind,igrid)=ipeak+c%npeak_cum(g%myid-1)
        c%peak_grid(ipeak)=igrid
        c%peak_cell(ipeak)=ind
        c%peak_level(ipeak)=c%level(itest)
        c%max_dens(ipeak)=m%rho(ind,igrid)
     endif
  end do

  !----------------------------------------
  ! Determine peak-patches around each peak
  !----------------------------------------
  nmove_tot=1
  istep=0
  do while (nmove_tot.gt.0)

     call open_cache(mdl, m, pack_size=storage_size(dummy_int4)/32, &
          pack=pack_fetch_flag, unpack=unpack_fetch_flag, bound=init_bound_flag)

     nmove=0
     nzero=0
     do itest=1,c%ntest
        hash_nbor=c%hash(itest,0:ndim)
        if(hash_nbor(0)>0)then
           igrid=c%grid(itest)
           ind=c%cell(itest)
           call get_parent_cell(s,hash_nbor,igridn,icelln,flush_cache=.false.,fetch_cache=.true.)
           if(m%flag1(ind,igrid).ne.m%flag1(icelln,igridn))nmove=nmove+1
           m%flag1(ind,igrid)=m%flag1(icelln,igridn)
           if(m%flag1(ind,igrid).eq.0)nzero=nzero+1
        endif
     end do

     call close_cache(mdl)

     istep=istep+1
     nmove_tot=int(nmove,kind=8)
     nzero_tot=int(nzero,kind=8)
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(nmove_tot,nmove_all,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(nzero_tot,nzero_all,1,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
     nmove_tot=nmove_all
     nzero_tot=nzero_all
#endif
     if(c%ntest_tot>0.and.r%verbose.and.g%myid==1)write(*,*)"istep=",istep,"nmove=",nmove_tot
  end do

  end associate

end subroutine collect_patch
!################################################################
!################################################################
!################################################################
!################################################################
#endif
end module clump_finder_module
