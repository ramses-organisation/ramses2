module sink_formation_module
#ifdef _CUDA
  use clump_device, only: gpu_sink_formation, gpu_sink_clump
#endif

  type :: out_sink_formation_t
     real(kind=8)::mass
  end type out_sink_formation_t

contains
!###########################################################
!###########################################################
!###########################################################
!###########################################################
recursive subroutine m_sink_formation(pst)
  use mdl_module
  use ramses_commons, only: pst_t
  use clump_merger_module, only: r_deallocate_clump
  use mdl_parameters
  implicit none
  type(pst_t)::pst

  type(out_sink_formation_t)::output_sink
  double precision::ttend, ttstart

  write(*,*)'Entering sink formation'
  ttstart = mdl_wtime(pst%s%mdl)

  !----------------------------
  ! Find sink formation sites
  !----------------------------
  call m_formation_site(pst)

  !----------------------------
  ! Create sink particles
  !----------------------------
  if(pst%s%c%npeak_tot>0)then
     call r_sink_formation(pst,pst%s%r%levelmin,1,output_sink,2)
     if(output_sink%mass>0)then
        pst%s%g%mass_sink_tot=pst%s%g%mass_sink_tot+output_sink%mass
     endif
  endif

  !----------------------------
  ! Output sink properties
  !----------------------------
!  if(pst%s%r%sink_dump)call dump_sink_particles(pst)

  !----------------------------
  ! Deallocate all peak arrays
  !----------------------------
  call r_deallocate_clump(pst)

  ttend = mdl_wtime(pst%s%mdl)
  print '(A,F0.7)',' Time elapsed in creating sinks: ',ttend-ttstart

end subroutine m_sink_formation
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine dump_sink_particles(pst)
    use amr_parameters, only: ndim,flen
    use ramses_commons, only: pst_t
    use output_part_module, only: r_output_part
    use output_clump_module, only: r_output_clump
    use mdl_module, only: mdl_mkdir
    use output_amr_module, only: r_output_amr,output_info
    use output_poisson_module, only: r_output_poisson,in_output_poisson_t
    implicit none
    type(pst_t)::pst
    type(in_output_poisson_t)::in_output_poisson
    ! Local variables
    integer::i,dummy(1)
    character(LEN=flen)::filename,filedir,filecmd
    integer,dimension(1:flen/4)::input_array
    character(len=20) :: str

    filedir='dump/'
    call mdl_mkdir(pst%s%mdl,filedir)
    write(str,'(I0)') pst%s%g%nstep_coarse
    filedir='dump/'//TRIM(str)//'_'

    filename=TRIM(filedir) ! Note that suffix will be added later
    input_array=transfer(filename,input_array)
    !in_output_poisson%filename=TRIM(filedir)//'grav.'
    if(pst%s%r%verbose)write(*,*)'Writing particle files'
    if(pst%s%c%npeak_tot>0)then
       call r_output_clump(pst,input_array,flen/4,dummy,0)
    endif

    !call r_output_poisson(pst,in_output_poisson,storage_size(in_output_poisson)/32)
    !call r_output_part(pst,input_array,flen/4,dummy,0)
    !filename=TRIM(filedir)//'amr.'
    !input_array=transfer(filename,input_array)
    !call r_output_amr(pst,input_array,flen/4,dummy,0)
    !filename=TRIM(filedir)//'info.txt'
    !call output_info(pst%s%r,pst%s%g,filename)

end subroutine dump_sink_particles
!###########################################################
!###########################################################
!###########################################################
!###########################################################
recursive subroutine r_sink_formation(pst,ilevel,input_size,output,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  type(out_sink_formation_t)::output,next_output

  integer::ilevel
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_SINK_FORMATION,pst%iUpper+1,input_size,output_size,ilevel)
     call r_sink_formation(pst%pLower,ilevel,input_size,output,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_output)
     output%mass=output%mass+next_output%mass
  else
#ifdef _CUDA
     call gpu_sink_formation(pst%s,output%mass)
#else
     call sink_formation(pst%s%r,pst%s%g,pst%s%m,pst%s%sink,pst%s%c,output%mass)
#endif
  endif

end subroutine r_sink_formation
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine sink_formation(r,g,m,p,c,msink_loc)
  use rng
  use constants
  use hydro_parameters, only: nvar
  use amr_parameters, only: ndim, twotondim
  use amr_commons, only: run_t, global_t, mesh_t
  use pm_commons, only: part_t
  use clfind_commons, only: clump_t
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  type(clump_t)::c
  real(kind=8)::msink_loc
  !-------------------------------------------------------------------
  ! Spawn sink particles from clumps using various formation criteria.
  ! We use the RAMSES clump finder PHEW for the clumps detection.
  !-------------------------------------------------------------------
#ifndef WITHOUTMPI
  integer::info
  integer,dimension(1:g%ncpu)::nsink_cpu_tot
#endif
  real(kind=8)::dx,vol,scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m
  integer(kind=8),dimension(0:g%ncpu)::nsink_cum
  integer,dimension(1:g%ncpu)::nsink_cpu
  integer::i,j,icpu,ilevel,nsink,nsink_loc,peak_nr,igrid,ind
  real(kind=8)::purity,mgas,mseed,dgas,dseed,mass_frac
  logical::ok

#if NDIM>2
  ! Conversion factor from user units to cgs units
  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=scale_d*scale_l**3
  mseed=r%sink_mseed*M_sun/scale_m
  !---------------------------
  ! Flag sink formation sites
  !---------------------------
  c%form_sink=0
  ! Loop over peaks
  do j=1,c%npeak
     ilevel=c%peak_level(j)
     igrid=c%peak_grid(j)
     ind=c%peak_cell(j)
     dx=r%boxlen/2**ilevel
     vol=dx**ndim
     ok=.true.
     !-------------------------------------
     ! Add here all sink formation criteria
     !-------------------------------------
     if(c%saddle_threshold>0)then
        if(c%ind_halo(j).NE.j+c%npeak_cum(g%myid-1))ok=.false.
     endif
     if(c%relevance(j)<=c%relevance_threshold)ok=.false.
     if(c%clump_mass(j)<=c%mass_threshold)ok=.false.
     if(c%nsink(j)>0)ok=.false.
     if(r%rho_type_sink==1)then
        if(c%particle_mass(j)>0)then
           purity=c%npart(j)*g%mp_min/c%particle_mass(j)
           if(purity<=c%purity_threshold)ok=.false.
        endif
        if(c%npart(j)==0)ok=.false.
     endif
#ifdef HYDRO
     dgas=m%uold(ind,1,igrid)*scale_nH
     if(dgas<=r%sink_nstar_frac*r%n_star)ok=.false.
     mgas=m%uold(ind,1,igrid)*vol
     if(mseed>0.75*mgas)ok=.false.
#endif
     ! Set sink formation flag
     if(ok)c%form_sink(j)=1
  end do

  !--------------------------
  ! Create new sink particles
  !--------------------------
  nsink_loc=0
  msink_loc=0.0d0
  ! Loop over peaks
  do j=1,c%npeak
     if(c%form_sink(j).eq.1)then
        ilevel=c%peak_level(j)
        igrid=c%peak_grid(j)
        ind=c%peak_cell(j)
        dx=r%boxlen/2**ilevel
        vol=dx**ndim
        nsink_loc=nsink_loc+1
        p%npart=p%npart+1
        if(p%npart>r%nsinkmax)then
           write(*,*)'Not enough memory for sink particle'
           write(*,*)'Increase nsinkmax in the namelist'
#ifndef WITHOUTMPI
           call MPI_ABORT(MPI_COMM_WORLD,1,info)
#endif
           stop
        endif
        ! Compute sink particle position from peak position
        p%xp(p%npart,1)=c%peak_com(j,1)
        p%xp(p%npart,2)=c%peak_com(j,2)
        p%xp(p%npart,3)=c%peak_com(j,3)
        ! Compute sink particle velocity from peak velocity
        p%vp(p%npart,1)=c%peak_vel(j,1)-c%peak_acc(j,1)*0.5d0*g%dtnew(c%peak_level(j))
        p%vp(p%npart,2)=c%peak_vel(j,2)-c%peak_acc(j,2)*0.5d0*g%dtnew(c%peak_level(j))
        p%vp(p%npart,3)=c%peak_vel(j,3)-c%peak_acc(j,3)*0.5d0*g%dtnew(c%peak_level(j))
        ! Set initial angular momentum to almost zero
        p%jp(p%npart,1)=0.0d0
        p%jp(p%npart,2)=0.0d0
        p%jp(p%npart,3)=1.0d-20
        ! Compute sink particle old force from peak acceleration
        p%fp(p%npart,1)=c%peak_acc(j,1)
        p%fp(p%npart,2)=c%peak_acc(j,2)
        p%fp(p%npart,3)=c%peak_acc(j,3)
        ! Compute sink particle mass
        p%mp(p%npart)=0
#ifdef HYDRO
        p%mp(p%npart)=mseed
        dseed=mseed/vol
        mass_frac=1.0d0-dseed/m%uold(ind,1,igrid)
        m%uold(ind,1:nvar,igrid)=m%uold(ind,1:nvar,igrid)*mass_frac
        msink_loc=msink_loc+mseed
#endif
        ! Compute sink particle birth time using proper time
        p%tp(p%npart)=g%texp
        ! Set merging time to -1000
        p%tm(p%npart)=-1000
        ! Set merging id to 0
        p%idm(p%npart)=0
        ! Compute level
        p%levelp(p%npart)=ilevel
     endif
  end do
  p%tailp(r%nlevelmax)=p%tailp(r%nlevelmax)+nsink_loc

  !---------------------------------------------------------
  ! Compute number of new sinks across all CPUs.
  !---------------------------------------------------------
  nsink_cpu=0
  nsink_cpu(g%myid)=nsink_loc
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(nsink_cpu,nsink_cpu_tot,g%ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  nsink_cpu=nsink_cpu_tot
#endif
  nsink_cum=0
  do icpu=1,g%ncpu
     nsink_cum(icpu)=nsink_cum(icpu-1)+int(nsink_cpu(icpu),kind=8)
  end do

  !--------------------------------------
  ! Compute new sink particle index
  !--------------------------------------
  do i=1,nsink_loc
     p%idp(p%npart-nsink_loc+i)=p%npart_tot+nsink_cum(g%myid-1)+i
  end do
  p%npart_tot=p%npart_tot+nsink_cum(g%ncpu)

  if(g%myid==1)write(*,*)'Formed',int(nsink_cum(g%ncpu),kind=4),' new sinks for a total of',int(p%npart_tot,kind=4)

#endif

end subroutine sink_formation
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine m_formation_site(pst)
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
  ! This is the master routine for the RAMSES sink formation sites finder
  !----------------------------------------------------------------------

#if NDIM==3 && defined(GRAV)

  associate(r=>pst%s%r,g=>pst%s%g,mdl=>pst%s%mdl,p=>pst%s%p,star=>pst%s%star)

  !--------------------------------------------------------------
  ! Compute rho from gas density or dark matter or star particles
  !--------------------------------------------------------------
  call m_rho_fine(pst,r%levelmin,r%rho_type_sink)

  !----------------------------------------------
  ! Find relevant peak patches as formation sites
  !----------------------------------------------
  call r_sink_clump(pst,r%levelmin,1)

  end associate
#endif

end subroutine m_formation_site
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_sink_clump(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::ilevel

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_SINK_CLUMP,pst%iUpper+1,input_size,0,ilevel)
     call r_sink_clump(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
#ifdef _CUDA
     call gpu_sink_clump(pst%s)
#else
     call sink_clump(pst%s)
#endif
  endif

end subroutine r_sink_clump
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine sink_clump(s)
  use ramses_commons, only: ramses_t
  use clump_finder_module
  use clump_merger_module
  implicit none
  type(ramses_t)::s

#if NDIM==3 && defined(GRAV)

  !-----------------------------------------------
  ! Load clump finder parameters in clump object.
  !-----------------------------------------------
  s%c%relevance_threshold = s%r%sink_relevance_threshold
  s%c%density_threshold = s%r%sink_density_threshold
  s%c%saddle_threshold = s%r%sink_saddle_threshold
  s%c%mass_threshold = s%r%sink_mass_threshold
  s%c%purity_threshold = s%r%sink_purity_threshold
  s%c%fraction_threshold = s%r%sink_fraction_threshold
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
  call compute_clump_properties(s,s%r%rho_type_sink)
  !----------------------------------------------------------------------
  ! Compute additional particle-based clump properties.
  !----------------------------------------------------------------------
  if(s%r%part)then
     call particle_peak_id(s,s%p)
     if(s%r%rho_type_sink.eq.1)then
        call particle_clump_properties(s,s%p)
     endif
  endif
  if(s%r%star)then
     call particle_peak_id(s,s%star)
     if(s%r%rho_type_sink.eq.2)then
        call particle_clump_properties(s,s%star)
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
  ! Reset sink position to peak position
  ! Count sinks in each clump hierarchically.
  !---------------------------------------------
  call particle_peak_id(s,s%sink)
  call sink_in_peak(s,s%r%sink_recenter_on_peak,.true.)

#endif
end subroutine sink_clump
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine sink_in_peak(s,reset_sink_pos,count_sink)
  use ramses_commons, only: ramses_t
  use clump_merger_module
  use cache_commons, only: msg_sink_clump,msg_saddle_clump
  use cache
  implicit none
  type(ramses_t)::s
  logical::count_sink,reset_sink_pos

  integer::i,ipart,peak_nr,halo_nr,ilev,ipeak,jpeak
  integer(kind=8)::global_peak_id,merge_to
  type(msg_sink_clump)::dummy_sink_clump

  associate(r=>s%r,g=>s%g,m=>s%m,c=>s%c,p=>s%sink,mdl=>s%mdl)

  !--------------------------------------------
  ! Sort particles according to global clump id
  !--------------------------------------------
  do ipart=1,p%npart
     p%sortp(ipart)=ipart
     p%workp(ipart)=p%pid(ipart)
  end do
  call quick_sort_int_int(p%workp(1),p%sortp(1),p%npart)

  !-------------------------------------
  ! Reset sink position to peak position
  !-------------------------------------
  if(reset_sink_pos)then
     call open_cache_clump(mdl, c, pack_size=storage_size(dummy_sink_clump)/32, &
          pack=pack_fetch_sink, unpack=unpack_fetch_sink)
     do i=1+p%norphan_peak,p%npart
        ipart=p%sortp(i)
        global_peak_id=p%workp(i)
        call get_peak(s,global_peak_id,peak_nr,fetch_cache=.true.,flush_cache=.false.)
        ! Compute sink particle position from peak position
        p%xp(ipart,1)=c%peak_com(peak_nr,1)
        p%xp(ipart,2)=c%peak_com(peak_nr,2)
        p%xp(ipart,3)=c%peak_com(peak_nr,3)
        ! Compute sink particle velocity from peak velocity
        p%vp(ipart,1)=c%peak_vel(peak_nr,1)-c%peak_acc(peak_nr,1)*0.5d0*g%dtnew(c%peak_level(peak_nr))
        p%vp(ipart,2)=c%peak_vel(peak_nr,2)-c%peak_acc(peak_nr,2)*0.5d0*g%dtnew(c%peak_level(peak_nr))
        p%vp(ipart,3)=c%peak_vel(peak_nr,3)-c%peak_acc(peak_nr,3)*0.5d0*g%dtnew(c%peak_level(peak_nr))
        ! Compute sink particle old force from peak acceleration
        p%fp(ipart,1)=c%peak_acc(peak_nr,1)
        p%fp(ipart,2)=c%peak_acc(peak_nr,2)
        p%fp(ipart,3)=c%peak_acc(peak_nr,3)
     end do
     call close_cache(mdl)
  endif

  !------------------------------------
  ! Count sink particles in each halo
  !------------------------------------
  if(count_sink)then
     ! Count sinks in each peak
     c%nsink=0
     call open_cache_clump(mdl, c, pack_size=storage_size(dummy_sink_clump)/32, &
          init=init_flush_sink, flush=pack_flush_sink, combine=unpack_flush_sink)
     do i=1+p%norphan_peak,p%npart
        global_peak_id=p%workp(i)
        call get_peak(s,global_peak_id,peak_nr,fetch_cache=.false.,flush_cache=.true.)
        c%nsink(peak_nr)=c%nsink(peak_nr)+1
     end do
     call close_cache(mdl)
     ! Count sinks hierarchically in each halo
     if(c%saddle_threshold>0)then
        do ilev=0,c%merge_levelmax
           call open_cache_clump(mdl, c, pack_size=storage_size(dummy_sink_clump)/32, &
                init=init_flush_sink, flush=pack_flush_sink, combine=unpack_flush_sink)
           do ipeak=1,c%npeak
              if(c%lev_peak(ipeak)==ilev)then
                 merge_to=c%new_peak(ipeak)
                 if(merge_to.NE.ipeak+c%npeak_cum(g%myid-1))then
                    call get_peak(s,merge_to,jpeak,flush_cache=.true.,fetch_cache=.false.)
                    c%nsink(jpeak)=c%nsink(jpeak)+c%nsink(ipeak)
                 endif
              endif
           end do
           call close_cache(mdl)
        end do
     endif
  endif

  end associate

end subroutine sink_in_peak
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_fetch_sink(c,local_peak_id,msg_size,msg_array)
  use amr_parameters, only: ndim
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_sink_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_sink_clump)::msg

  msg%lev = c%peak_level(local_peak_id)
  msg%pos(1:ndim) = c%peak_com(local_peak_id,1:ndim)
  msg%vel(1:ndim) = c%peak_vel(local_peak_id,1:ndim)
  msg%acc(1:ndim) = c%peak_acc(local_peak_id,1:ndim)

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_sink
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_fetch_sink(c,local_peak_id,msg_size,msg_array)
  use amr_parameters, only: ndim
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_sink_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_sink_clump)::msg

  msg=transfer(msg_array,msg)

  c%peak_level(local_peak_id)=msg%lev
  c%peak_com(local_peak_id,1:ndim)=msg%pos(1:ndim)
  c%peak_vel(local_peak_id,1:ndim)=msg%vel(1:ndim)
  c%peak_acc(local_peak_id,1:ndim)=msg%acc(1:ndim)

end subroutine unpack_fetch_sink
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_flush_sink(c,local_peak_id)
  use clfind_commons, only: clump_t
  type(clump_t)::c
  integer::local_peak_id

  c%nsink(local_peak_id)=0

end subroutine init_flush_sink
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pack_flush_sink(c,local_peak_id,msg_size,msg_array)
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_sink_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_sink_clump)::msg

  msg%lev=c%nsink(local_peak_id)

  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_sink
!################################################################
!################################################################
!################################################################
!################################################################
subroutine unpack_flush_sink(c,local_peak_id,msg_size,msg_array)
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_sink_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_sink_clump)::msg

  msg=transfer(msg_array,msg)

  c%nsink(local_peak_id)=c%nsink(local_peak_id)+msg%lev

end subroutine unpack_flush_sink
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
end module sink_formation_module
