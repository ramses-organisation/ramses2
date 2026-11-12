module tree_formation_module
#ifdef _CUDA
  use clump_device, only: gpu_tree_clump, gpu_tree_formation
#endif

contains
!###########################################################
!###########################################################
!###########################################################
!###########################################################
recursive subroutine m_tree_formation(pst)
  use mdl_module
  use ramses_commons, only: pst_t
  use clump_merger_module, only: r_deallocate_clump
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  double precision::ttend, ttstart

  write(*,*)'Entering merger tree particle formation'
  ttstart = mdl_wtime(pst%s%mdl)

  !----------------------------
  ! Find tree formation sites
  !----------------------------
  call m_tree_formation_site(pst)

  !----------------------------
  ! Create tree particles
  !----------------------------
  if(pst%s%c%npeak_tot>0)then
     call r_tree_formation(pst)
  endif

  !------------------------------
  ! Deallocate all peak arrays
  !------------------------------
  call r_deallocate_clump(pst)

  ttend = mdl_wtime(pst%s%mdl)
  print '(A,F0.7)',' Time elapsed in creating trees: ',ttend-ttstart

end subroutine m_tree_formation
!###########################################################
!###########################################################
!###########################################################
!###########################################################
recursive subroutine r_tree_formation(pst)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_TREE_FORMATION,pst%iUpper+1)
     call r_tree_formation(pst%pLower)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
#ifdef _CUDA
     call gpu_tree_formation(pst%s)
#else
     call tree_formation(pst%s%r,pst%s%g,pst%s%m,pst%s%tree,pst%s%c)
#endif
  endif

end subroutine r_tree_formation
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine tree_formation(r,g,m,p,c)
  use rng
  use constants
  use hydro_parameters, only: nvar
  use amr_parameters, only: ndim, twotondim
  use amr_commons, only: run_t, global_t, mesh_t
  use pm_commons, only: part_t
  use clfind_commons, only:clump_t
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  type(clump_t)::c
  !-------------------------------------------------------------------
  ! Spawn tree particles from clumps using various formation criteria.
  ! We use the RAMSES clump finder PHEW for the clumps detection.
  !-------------------------------------------------------------------
#ifndef WITHOUTMPI
  integer::info
  integer,dimension(1:g%ncpu)::nsite_cpu_tot,ntree_cpu_tot
#endif
  integer(kind=8),dimension(0:g%ncpu)::nsite_cum,ntree_cum
  integer,dimension(1:g%ncpu)::nsite_cpu,ntree_cpu
  integer::i,j,icpu,nsite,ntree,ntree_loc,peak_nr,idim
  real(kind=8)::purity
  logical::ok

#if NDIM>2
  !---------------------------
  ! Count tree formation sites
  !---------------------------
  nsite=0
  c%form_tree=0
  ! Loop over peaks
  do j=1,c%npeak
     ok=.true.
     !-------------------------------------
     ! Add here all tree formation criteria
     !-------------------------------------
     if(c%saddle_threshold>0)then
        if(c%ind_halo(j).NE.j+c%npeak_cum(g%myid-1))ok=.false.
        if(c%halo_mass(j)<=10d0*c%mass_threshold)ok=.false.
     else
        if(c%clump_mass(j)<=10d0*c%mass_threshold)ok=.false.
     endif
     if(c%relevance(j)<=c%relevance_threshold)ok=.false.
     if(c%npart(j)==0)ok=.false.
     if(c%ntree(j)>0)ok=.false.
     if(c%particle_mass(j)>0)then
        purity=c%npart(j)*g%mp_min/c%particle_mass(j)
        if(purity<=c%purity_threshold)ok=.false.
     endif
     ! Set tree formation flag
     if(ok)c%form_tree(j)=1
     if(ok)nsite=nsite+1
  end do

  !---------------------------------------------------------
  ! Compute number of tree formation sites across all CPUs.
  !---------------------------------------------------------
  nsite_cpu=0
  nsite_cpu(g%myid)=nsite
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(nsite_cpu,nsite_cpu_tot,g%ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  nsite_cpu=nsite_cpu_tot
#endif
  nsite_cum=0
  do icpu=1,g%ncpu
     nsite_cum(icpu)=nsite_cum(icpu-1)+int(nsite_cpu(icpu),kind=8)
  end do

  !--------------------------
  ! Create new tree particles
  !--------------------------
  ntree_loc=0
  ! Loop over peaks
  do j=1,c%npeak
     if(c%form_tree(j).eq.1)then
        ntree_loc=ntree_loc+1
        p%npart=p%npart+1
        if(p%npart>r%ntreemax)then
           write(*,*)'Not enough memory for tree particle'
           write(*,*)'Increase ntreemax in the namelist'
           stop
        endif
        if(allocated(p%phip))p%phip(p%npart)=0.0
        ! Compute tree particle position from peak position
        p%xp(p%npart,1)=c%peak_com(j,1)
        p%xp(p%npart,2)=c%peak_com(j,2)
        p%xp(p%npart,3)=c%peak_com(j,3)
        do idim=1,ndim
           if(r%periodic(idim))then
              if(p%xp(p%npart,idim)< 0.0d0           )p%xp(p%npart,idim)=p%xp(p%npart,idim)+r%box_size(idim)
              if(p%xp(p%npart,idim)>=r%box_size(idim))p%xp(p%npart,idim)=p%xp(p%npart,idim)-r%box_size(idim)
           endif
        end do
        ! Compute tree particle velocity from peak velocity
        p%vp(p%npart,1)=c%peak_vel(j,1)-c%peak_acc(j,1)*0.5d0*g%dtnew(c%peak_level(j))
        p%vp(p%npart,2)=c%peak_vel(j,2)-c%peak_acc(j,2)*0.5d0*g%dtnew(c%peak_level(j))
        p%vp(p%npart,3)=c%peak_vel(j,3)-c%peak_acc(j,3)*0.5d0*g%dtnew(c%peak_level(j))
        ! Set mass to zero
        p%mp(p%npart)=0d0
        ! Set tracking number to parent clump global id
        p%idt(p%npart)=j+c%npeak_cum(g%myid-1)
        ! Compute tree particle birth time using proper time
        p%tp(p%npart)=g%texp
        ! Set merging time to -1000
        p%tm(p%npart)=-1000d0
        ! Set merging id to 0
        p%idm(p%npart)=0
        ! Compute level
        p%levelp(p%npart)=c%peak_level(j)
     endif
  end do
  p%tailp(r%nlevelmax)=p%tailp(r%nlevelmax)+ntree_loc

  !---------------------------------------------------------
  ! Compute number of new trees across all CPUs.
  !---------------------------------------------------------
  ntree_cpu=0
  ntree_cpu(g%myid)=ntree_loc
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(ntree_cpu,ntree_cpu_tot,g%ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  ntree_cpu=ntree_cpu_tot
#endif
  ntree_cum=0
  do icpu=1,g%ncpu
     ntree_cum(icpu)=ntree_cum(icpu-1)+int(ntree_cpu(icpu),kind=8)
  end do

  !--------------------------------------
  ! Compute new tree particle index
  !--------------------------------------
  do i=1,ntree_loc
     p%idp(p%npart-ntree_loc+i)=p%npart_tot+ntree_cum(g%myid-1)+i
  end do
  p%npart_tot=p%npart_tot+ntree_cum(g%ncpu)

  if(g%myid==1)write(*,'(" Found",I12,"  new merger tree particles for a total of ",I0)')int(ntree_cum(g%ncpu),kind=4),int(p%npart_tot,kind=4)

#endif

end subroutine tree_formation
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine m_tree_formation_site(pst)
  use amr_parameters, only: flen
  use mdl_module, only: mdl_wtime
  use ramses_commons, only: pst_t
#ifdef GRAV
  use rho_fine_module, only: m_rho_fine
#endif
  implicit none
  type(pst_t)::pst
  logical::keep_alive
  !----------------------------------------------------------------------
  ! This is the master routine for the RAMSES tree formation sites finder
  !----------------------------------------------------------------------

#if NDIM==3 && defined(GRAV)
  !--------------------------------------------------------------
  ! Compute rho from gas density or dark matter or star particles
  !--------------------------------------------------------------
  call m_rho_fine(pst,pst%s%r%levelmin,pst%s%r%rho_type_clump)

  !----------------------------------------------
  ! Find relevant peak patches as formation sites
  !----------------------------------------------
  call r_tree_clump(pst)
#endif

end subroutine m_tree_formation_site
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_tree_clump(pst)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_TREE_CLUMP,pst%iUpper+1)
     call r_tree_clump(pst%pLower)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
#ifdef _CUDA
     call gpu_tree_clump(pst%s)
#else
     call tree_clump(pst%s)
#endif
  endif

end subroutine r_tree_clump
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine tree_clump(s)
  use ramses_commons, only: ramses_t
  use clump_finder_module
  use clump_merger_module
  implicit none
  type(ramses_t)::s

#if NDIM==3 && defined(GRAV)

  !-----------------------------------------------
  ! Load clump finder parameters in clump object.
  !-----------------------------------------------
  s%c%relevance_threshold = s%r%relevance_threshold
  s%c%density_threshold = s%r%density_threshold
  s%c%saddle_threshold = s%r%saddle_threshold
  s%c%mass_threshold = s%r%mass_threshold
  s%c%purity_threshold = s%r%purity_threshold
  s%c%fraction_threshold = 2d0
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
  ! algorithm. We get well defined peak patches around each peak.
  ! As a result, each pair of neighboring peak patches are separated
  ! by their saddle surface.
  !----------------------------------------------------------------------
  call collect_patch(s)
  !----------------------------------------------------------------------
  ! Allocate all peak patch based arrays
  !----------------------------------------------------------------------
  call allocate_peak_patch_arrays(s)
  !----------------------------------------------------------------------
  ! We compute the densest saddle point and its corresponding
  ! neighboring peak.
  !----------------------------------------------------------------------
  call collect_saddle(s)
  !----------------------------------------------------------------------
  ! Merge peaks based on a relevance criterion.
  ! Peaks that are due to random noise fluctuations or peaks that
  ! have similar peak density values are merged into relevant peaks
  !----------------------------------------------------------------------
  call merge_clumps(s,'relevance')
  !----------------------------------------------------------------------
  ! Compute relevant peak properties such as mass and number of cells
  !----------------------------------------------------------------------
  call compute_clump_properties(s,s%r%rho_type_clump)
  !----------------------------------------------------------------------
  ! Compute additional particle-based clump properties.
  !----------------------------------------------------------------------
  if(s%r%part)then
     call particle_peak_id(s,s%p)
     if(s%r%rho_type_clump.eq.1)then
        call particle_clump_properties(s,s%p)
     endif
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
  !---------------------------------------------
  ! Reset tree position to peak position
  ! Count trees in each clump hierarchically.
  !---------------------------------------------
  call particle_peak_id(s,s%tree)
  call tree_in_peak(s,.true.,.true.)

#endif
end subroutine tree_clump
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine tree_in_peak(s,reset_tree_pos,count_tree)
  use amr_parameters, only: ndim, i8b
  use ramses_commons, only: ramses_t
  use clump_merger_module
  use cache_commons, only: msg_tree_clump,msg_tree_minid
  use cache
  implicit none
  type(ramses_t)::s
  logical::count_tree,reset_tree_pos

  integer::i,ipart,idim,peak_nr,halo_nr,ilev,ipeak,jpeak
  integer(kind=8)::global_peak_id,merge_to
  type(msg_tree_clump)::dummy_tree_clump
  type(msg_tree_minid)::dummy_tree_minid

  logical::bound
  real(kind=8)::pi,grav,r2,v2,rad2,vel2,radius
  real(kind=8),dimension(1:ndim)::xpart

  associate(r=>s%r,g=>s%g,m=>s%m,c=>s%c,p=>s%tree,mdl=>s%mdl)

  ! Compute constants
  pi=ACOS(-1.0D0)
  grav=1d0
  if(s%r%cosmo)grav=3d0/8d0/pi*s%g%omega_m*s%g%aexp

  !-------------------------------------------------------
  ! Sort tree particles according to their global clump id
  !-------------------------------------------------------
  do ipart=1,p%npart
     p%sortp(ipart)=ipart
     p%workp(ipart)=p%pid(ipart)
  end do
  call quick_sort_int_int(p%workp(1),p%sortp(1),p%npart)

  !-----------------------------------------------
  ! Reset tree particle positions to peak position
  !-----------------------------------------------
  if(reset_tree_pos)then
     c%min_tree_id=huge(0_i8b)
     call open_cache_clump(mdl, c, pack_size=storage_size(dummy_tree_clump)/32, &
          pack=pack_fetch_tree, unpack=unpack_fetch_tree, &
          init=init_flush_minid, flush=pack_flush_minid, combine=unpack_flush_minid)
     do i=1+p%norphan_peak,p%npart
        ipart=p%sortp(i)
        global_peak_id=p%workp(i)
        call get_peak(s,global_peak_id,peak_nr,fetch_cache=.true.,flush_cache=.true.)
        ! Compute peak's central core properties
        radius = 2d0*r%boxlen/2**c%peak_level(peak_nr)
        rad2 = radius**2
        ! Compute relative velocity
        v2 =     (p%vp(ipart,1) - c%peak_vel(peak_nr,1))**2 &
             & + (p%vp(ipart,2) - c%peak_vel(peak_nr,2))**2 &
             & + (p%vp(ipart,3) - c%peak_vel(peak_nr,3))**2
        ! Compute relative distance
        xpart(1:ndim)=p%xp(ipart,1:ndim)-c%peak_com(peak_nr,1:ndim)
        r2=0d0
        do idim=1,ndim
           if(r%periodic(idim))then
              if(xpart(idim)> r%box_size(idim)*0.5d0)xpart(idim)=xpart(idim)-r%box_size(idim)
              if(xpart(idim)<-r%box_size(idim)*0.5d0)xpart(idim)=xpart(idim)+r%box_size(idim)
           endif
           r2=r2+xpart(idim)**2
        end do
        ! Compute boundness criteria (a peak without particle mass binds nothing)
        bound = .false.
        if(c%particle_mass(peak_nr)>0d0)then
           vel2 = grav*c%particle_mass(peak_nr)/radius
           bound = ( v2/vel2 + 2d0*sqrt(r2/rad2) < 20d0 )
        endif

        if(bound)then
           ! If not merged yet then mark as merger candidate
           if(p%idm(ipart)==0)p%idm(ipart)=-2
           ! Compute tree particle position from peak position
           p%xp(ipart,1)=c%peak_com(peak_nr,1)
           p%xp(ipart,2)=c%peak_com(peak_nr,2)
           p%xp(ipart,3)=c%peak_com(peak_nr,3)
           do idim=1,ndim
              if(r%periodic(idim))then
                 if(p%xp(ipart,idim)< 0.0d0           )p%xp(ipart,idim)=p%xp(ipart,idim)+r%box_size(idim)
                 if(p%xp(ipart,idim)>=r%box_size(idim))p%xp(ipart,idim)=p%xp(ipart,idim)-r%box_size(idim)
              endif
           end do
           ! Compute tree particle velocity from peak velocity
           p%vp(ipart,1)=c%peak_vel(peak_nr,1)-c%peak_acc(peak_nr,1)*0.5d0*g%dtnew(c%peak_level(peak_nr))
           p%vp(ipart,2)=c%peak_vel(peak_nr,2)-c%peak_acc(peak_nr,2)*0.5d0*g%dtnew(c%peak_level(peak_nr))
           p%vp(ipart,3)=c%peak_vel(peak_nr,3)-c%peak_acc(peak_nr,3)*0.5d0*g%dtnew(c%peak_level(peak_nr))
           ! Update minimum tree id in each clump
           c%min_tree_id(peak_nr)=min(c%min_tree_id(peak_nr),p%idp(ipart))
        endif
     end do
     call close_cache(mdl)
  endif

  !----------------------------------------------------
  ! Merge all tree particles that sit in the same clump
  !----------------------------------------------------
  call open_cache_clump(mdl, c, pack_size=storage_size(dummy_tree_minid)/32, &
       pack=pack_fetch_minid, unpack=unpack_fetch_minid)
  ! Set tracking id to zero for orphan merger tree tracer particles
  do i=1,p%norphan_peak
     ipart=p%sortp(i)
     p%idt(ipart)=0
  end do
  do i=1+p%norphan_peak,p%npart
     ipart=p%sortp(i)
     global_peak_id=p%workp(i)
     call get_peak(s,global_peak_id,peak_nr,fetch_cache=.true.,flush_cache=.false.)
     ! If tree particle is not bound and not merged, set as orphan with zero tracking id
     if(p%idm(ipart).EQ.0)then
        p%idt(ipart)=0
     endif
     ! If tree particle is bound and just merged, update merging age, merge-to clump id and clump tracking id
     if(p%idm(ipart).EQ.-2.AND.p%idp(ipart).NE.c%min_tree_id(peak_nr))then
        p%idm(ipart)=c%min_tree_id(peak_nr)
        p%tm(ipart)=g%texp
        p%idt(ipart)=global_peak_id
     endif
     ! If tree particle is bound and not merged, update its clump tracking id
     if(p%idm(ipart).EQ.-2.AND.p%idp(ipart).EQ.c%min_tree_id(peak_nr))then
        p%idm(ipart)=0
        p%idt(ipart)=global_peak_id
     endif
  end do
  call close_cache(mdl)

  !------------------------------------
  ! Count tree particles in each halo
  !------------------------------------
  if(count_tree)then
     ! Count trees in each peak
     c%ntree=0
     call open_cache_clump(mdl, c, pack_size=storage_size(dummy_tree_clump)/32, &
          init=init_flush_tree, flush=pack_flush_tree, combine=unpack_flush_tree)
     do i=1+p%norphan_peak,p%npart
        global_peak_id=p%workp(i)
        call get_peak(s,global_peak_id,peak_nr,fetch_cache=.false.,flush_cache=.true.)
        c%ntree(peak_nr)=c%ntree(peak_nr)+1
     end do
     call close_cache(mdl)
     ! Count trees hierarchically in each halo
     if(c%saddle_threshold>0)then
        do ilev=0,c%merge_levelmax
           call open_cache_clump(mdl, c, pack_size=storage_size(dummy_tree_clump)/32, &
                init=init_flush_tree, flush=pack_flush_tree, combine=unpack_flush_tree)
           do ipeak=1,c%npeak
              if(c%lev_peak(ipeak)==ilev)then
                 merge_to=c%new_peak(ipeak)
                 if(merge_to.NE.ipeak+c%npeak_cum(g%myid-1))then
                    call get_peak(s,merge_to,jpeak,flush_cache=.true.,fetch_cache=.false.)
                    c%ntree(jpeak)=c%ntree(jpeak)+c%ntree(ipeak)
                 endif
              endif
           end do
           call close_cache(mdl)
        end do
     endif
  endif

  end associate

end subroutine tree_in_peak
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_fetch_tree(c,local_peak_id,msg_size,msg_array)
  use amr_parameters, only: ndim
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_tree_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_tree_clump)::msg

  msg%lev = c%peak_level(local_peak_id)
  msg%rad = c%clump_rad(local_peak_id)
  msg%mass = c%particle_mass(local_peak_id)
  msg%pos(1:ndim) = c%peak_com(local_peak_id,1:ndim)
  msg%vel(1:ndim) = c%peak_vel(local_peak_id,1:ndim)
  msg%acc(1:ndim) = c%peak_acc(local_peak_id,1:ndim)

  msg_array = transfer(msg,msg_array)

end subroutine pack_fetch_tree
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_fetch_tree(c,local_peak_id,msg_size,msg_array)
  use amr_parameters, only: ndim
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_tree_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_tree_clump)::msg

  msg = transfer(msg_array,msg)

  c%peak_level(local_peak_id) = msg%lev
  c%clump_rad(local_peak_id) = msg%rad
  c%particle_mass(local_peak_id) = msg%mass
  c%peak_com(local_peak_id,1:ndim) = msg%pos(1:ndim)
  c%peak_vel(local_peak_id,1:ndim) = msg%vel(1:ndim)
  c%peak_acc(local_peak_id,1:ndim) = msg%acc(1:ndim)

end subroutine unpack_fetch_tree
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_flush_tree(c,local_peak_id)
  use clfind_commons, only: clump_t
  type(clump_t)::c
  integer::local_peak_id

  c%ntree(local_peak_id)=0

end subroutine init_flush_tree
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_flush_tree(c,local_peak_id,msg_size,msg_array)
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_tree_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_tree_clump)::msg

  msg%lev=c%ntree(local_peak_id)

  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_tree
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_flush_tree(c,local_peak_id,msg_size,msg_array)
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_tree_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_tree_clump)::msg

  msg=transfer(msg_array,msg)

  c%ntree(local_peak_id)=c%ntree(local_peak_id)+msg%lev

end subroutine unpack_flush_tree
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_fetch_minid(c,local_peak_id,msg_size,msg_array)
  use amr_parameters, only: ndim
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_tree_minid
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_tree_minid)::msg

  msg%id = c%min_tree_id(local_peak_id)

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_minid
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_fetch_minid(c,local_peak_id,msg_size,msg_array)
  use amr_parameters, only: ndim
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_tree_minid
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_tree_minid)::msg

  msg=transfer(msg_array,msg)

  c%min_tree_id(local_peak_id) = msg%id

end subroutine unpack_fetch_minid
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_flush_minid(c,local_peak_id)
  use amr_parameters, only: i8b
  use clfind_commons, only: clump_t
  type(clump_t)::c
  integer::local_peak_id

  c%min_tree_id(local_peak_id)=huge(0_i8b)

end subroutine init_flush_minid
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_flush_minid(c,local_peak_id,msg_size,msg_array)
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_tree_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_tree_clump)::msg

  msg%id=c%min_tree_id(local_peak_id)

  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_minid
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_flush_minid(c,local_peak_id,msg_size,msg_array)
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_tree_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_tree_clump)::msg

  msg=transfer(msg_array,msg)

  c%min_tree_id(local_peak_id)=min(c%min_tree_id(local_peak_id),msg%id)

end subroutine unpack_flush_minid
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
end module tree_formation_module
