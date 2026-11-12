module core_formation_module
  use amr_parameters, only: ndim
#ifdef _CUDA
  use phew_device, only: gpu_core_clump, gpu_core_formation
#endif

  integer,parameter::CORE_CLUMP_FIND=1
  integer,parameter::CORE_CLUMP_RESTORE=2

  type :: out_core_formation_t
     real(kind=8)::mass
  end type out_core_formation_t

  ! Cache message of the core criteria: clump centre in, virial terms out
  type :: msg_core_clump
     real(kind=8)::mass4,kin,therm,mag,grav,rad,idot
     real(kind=8),dimension(1:ndim)::pos,com,vel
     real(kind=8),dimension(1:ndim,1:ndim)::inertia,idot_tensor
  end type msg_core_clump

  ! flag1 and flag2 of the refinement map, kept while the clump finder uses them
  integer,allocatable,dimension(:,:)::flag1_keep,flag2_keep

  ! Gas mass, centre of mass, velocity, dI/dt and inertia tensors of each peak
  real(kind=8),allocatable,dimension(:)::gas_mass,gas_idot
  real(kind=8),allocatable,dimension(:,:)::gas_com,gas_vel
  real(kind=8),allocatable,dimension(:,:,:)::gas_inertia,gas_idot_tensor

contains
!###########################################################
!###########################################################
!###########################################################
!###########################################################
recursive subroutine m_core_formation(pst)
  use ramses_commons, only: pst_t
  use clump_merger_module, only: r_deallocate_clump
  implicit none
  type(pst_t)::pst

  type(out_core_formation_t)::output_core

  !-----------------------------------------------------
  ! Create core particles in the clumps found this step
  !-----------------------------------------------------
  if(pst%s%c%npeak_tot>0)then
     call r_core_formation(pst,pst%s%g%nstep_coarse,1,output_core,2)
     if(output_core%mass>0d0)then
        pst%s%g%mass_sink_tot=pst%s%g%mass_sink_tot+output_core%mass
     endif
  endif

  !------------------------------
  ! Deallocate all peak arrays
  !------------------------------
  call r_deallocate_clump(pst)

end subroutine m_core_formation
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine m_core_formation_site(pst)
  use ramses_commons, only: pst_t
#ifdef GRAV
  use rho_fine_module, only: m_rho_fine
  use core_clump_output_module, only: output_core_clump_catalogue
#endif
  implicit none
  type(pst_t)::pst
  !----------------------------------------------------------------------
  ! This is the master routine for the core formation sites finder.
  !----------------------------------------------------------------------

#ifdef GRAV
  !--------------------------------------------------------------
  ! Compute rho from the core source
  ! (1: DM, 2: stars, 3: BH sinks, 4: gas)
  !--------------------------------------------------------------
  call m_rho_fine(pst,pst%s%r%levelmin,pst%s%r%rho_type_core)

  !----------------------------------------------
  ! Find relevant peak patches as formation sites
  !----------------------------------------------
  call r_core_clump(pst,CORE_CLUMP_FIND,1)

  !------------------------------------------------------------
  ! Write the core clump catalogue while flag1 and flag2 still
  ! hold the halo and peak ids, then restore the refinement map
  !------------------------------------------------------------
  if(pst%s%c%npeak_tot>0)call output_core_clump_catalogue(pst)
  call r_core_clump(pst,CORE_CLUMP_RESTORE,1)
#endif

end subroutine m_core_formation_site
!################################################################
!################################################################
!################################################################
!################################################################
recursive subroutine r_core_clump(pst,action,input_size)
  use mdl_module
  use amr_parameters, only: twotondim
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::action

  integer::rID,n
#ifdef _CUDA
  real(kind=8)::threshold
#endif

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_CORE_CLUMP,pst%iUpper+1,input_size,0,action)
     call r_core_clump(pst%pLower,action,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     select case(action)
     case(CORE_CLUMP_FIND)
        ! Keep flag1 and flag2, which hold the refinement map
        n=pst%s%m%noct_used
        allocate(flag1_keep(1:twotondim,1:n),flag2_keep(1:twotondim,1:n))
        flag1_keep=pst%s%m%flag1(:,1:n)
        flag2_keep=pst%s%m%flag2(:,1:n)
#ifdef _CUDA
        call core_clump_collection_threshold(pst%s%r,pst%s%g,threshold)
        call gpu_core_clump(pst%s,threshold)
#else
        call core_clump(pst%s)
#endif
     case(CORE_CLUMP_RESTORE)
        n=size(flag1_keep,2)
        pst%s%m%flag1(:,1:n)=flag1_keep
        pst%s%m%flag2(:,1:n)=flag2_keep
        deallocate(flag1_keep,flag2_keep)
     end select
  endif

end subroutine r_core_clump
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine core_clump(s)
  use ramses_commons, only: ramses_t
  use clump_finder_module
  use clump_merger_module
  implicit none
  type(ramses_t)::s

#if NDIM==3 && defined(GRAV)

  !-----------------------------------------------
  ! Load clump finder parameters in clump object.
  !-----------------------------------------------
  s%c%relevance_threshold = s%r%core_relevance_threshold
  call core_clump_collection_threshold(s%r,s%g,s%c%density_threshold)
  s%c%saddle_threshold = s%r%core_saddle_threshold
  s%c%mass_threshold = s%r%core_mass_threshold
  s%c%fraction_threshold = s%r%core_fraction_threshold
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
  ! Compute relevant peak properties such as mass and number of cells,
  ! and the core energy and contraction criteria
  !----------------------------------------------------------------------
  call compute_clump_properties(s,s%r%rho_type_core)
  if(s%r%core_use_clump_prop_criteria.or.s%r%core_use_contracting_all_axes.or. &
       s%r%output_core_clump)then
     call compute_core_clump_criteria(s,s%r%rho_type_core)
  endif
  !----------------------------------------------------------------------
  ! Compute additional particle-based clump properties.
  !----------------------------------------------------------------------
  if(s%r%part)then
     call particle_peak_id(s,s%p)
     if(s%r%rho_type_core.eq.1)then
        call particle_clump_properties(s,s%p)
     endif
  endif
  if(s%r%star)then
     call particle_peak_id(s,s%star)
     if(s%r%rho_type_core.eq.2)then
        call particle_clump_properties(s,s%star)
     endif
  endif
  if(s%r%sink)then
     call particle_peak_id(s,s%sink)
     if(s%r%rho_type_core.eq.3)then
        call particle_clump_properties(s,s%sink)
     endif
  endif
  !----------------------------------------------------------------------
  ! Merge all neighboring peak-patches above the prescribed density
  ! threshold into halo-patches, if their saddle point density is larger
  ! than the prescribed saddle density threshold.
  !----------------------------------------------------------------------
  if(s%c%saddle_threshold>0)then
     call merge_clumps(s,'saddleden')
  endif
  !---------------------------------------------
  ! Count cores in each clump hierarchically.
  !---------------------------------------------
  call particle_peak_id(s,s%core)
  call count_existing_cores_in_peak(s)
  !----------------------------------------------------------------------
  ! Compute particle peak and halo ids for the core clump catalogue.
  !----------------------------------------------------------------------
  if(s%g%output_done)then
     if(s%r%tree)then
        call particle_peak_id(s,s%tree)
     endif
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
        call particle_halo_id(s,s%core)
        if(s%r%tree)then
           call particle_halo_id(s,s%tree)
        endif
     else
        ! Without halo merging, each clump is its own halo
        if(s%r%part)s%p%hid=s%p%pid
        if(s%r%star)s%star%hid=s%star%pid
        if(s%r%sink)s%sink%hid=s%sink%pid
        s%core%hid=s%core%pid
        if(s%r%tree)s%tree%hid=s%tree%pid
     endif
  endif

#endif
end subroutine core_clump
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine core_clump_collection_threshold(r,g,density_threshold)
  use amr_commons, only: run_t, global_t
  use constants, only: pi
  implicit none
  type(run_t)::r
  type(global_t)::g
  real(kind=8),intent(out)::density_threshold
  !------------------------------------------------------------------
  ! Clump finder collection threshold for cores, in code density.
  !------------------------------------------------------------------
  real(kind=8)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(kind=8)::clump_threshold_nH,aexp_ref
  real(kind=8)::dx_lp,c_s2_lp,factG_lp,d_lp

  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  if(r%core_clump_density_threshold_LP_frac>0d0)then
     ! Fraction of the LP density at nlevelmax for the given T and mu
     dx_lp=r%boxlen/2d0**r%nlevelmax
     c_s2_lp=r%gamma*(r%core_clump_density_threshold_LP_T &
          & /r%core_clump_density_threshold_LP_mu)/scale_T2
     factG_lp=1d0
     if(r%cosmo)factG_lp=3d0/8d0/pi*g%omega_m*g%aexp
     d_lp=8.86d0*c_s2_lp/pi/factG_lp/dx_lp**2
     density_threshold = r%core_clump_density_threshold_LP_frac*d_lp
  else if(r%core_clump_density_threshold_nH>0d0)then
     ! Fixed H/cc threshold, scaled as (aexp_ref/aexp)^2 in cosmology
     clump_threshold_nH=r%core_clump_density_threshold_nH
     if(r%cosmo)then
        aexp_ref=1d0/(1d0+r%core_clump_density_threshold_zref)
        clump_threshold_nH=clump_threshold_nH*(aexp_ref/g%aexp)**2
     endif
     density_threshold = clump_threshold_nH/scale_nH
  else if(r%core_clump_density_threshold>0d0)then
     density_threshold = r%core_clump_density_threshold/scale_nH
  else
     ! Default: a tenth of the core formation density threshold
     density_threshold = r%core_density_threshold/scale_nH/1.0d1
  endif
  if(r%verbose.and.g%myid==1)write(*,'(A,ES16.8,A,ES16.8)') &
       & ' Core clump collection threshold [H/cc]=', &
       & density_threshold*scale_nH,'  code=',density_threshold

end subroutine core_clump_collection_threshold
!###########################################################
!###########################################################
!###########################################################
!###########################################################
recursive subroutine r_core_formation(pst,nstep_coarse,input_size,output,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  type(out_core_formation_t)::output,next_output

  integer::nstep_coarse
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_CORE_FORMATION,pst%iUpper+1,input_size,output_size,nstep_coarse)
     call r_core_formation(pst%pLower,nstep_coarse,input_size,output,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_output)
     output%mass=output%mass+next_output%mass
  else
#ifdef _CUDA
     call gpu_core_formation(pst%s,nstep_coarse,output%mass)
#else
     call core_formation(pst%s%mdl,pst%s%r,pst%s%g,pst%s%m,pst%s%core,pst%s%c,nstep_coarse,output%mass)
#endif
  endif

end subroutine r_core_formation
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine core_formation(mdl,r,g,m,p,c,nstep_coarse,mcore_loc)
  use mdl_module, only: mdl_t, mdl_abort
  use amr_parameters, only: ndim
  use amr_commons, only: run_t, global_t, mesh_t
  use pm_commons, only: part_t
  use clfind_commons, only: clump_t
  use constants, only: M_sun
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(mdl_t)::mdl
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::p
  type(clump_t)::c
  integer::nstep_coarse
  real(kind=8)::mcore_loc
  !-------------------------------------------------------------------
  ! Spawn core particles from clumps.
  !-------------------------------------------------------------------
#ifndef WITHOUTMPI
  integer::info
  integer,dimension(1:g%ncpu)::ncore_cpu_tot
#endif
  integer(kind=8),dimension(0:g%ncpu)::ncore_cum
  integer,dimension(1:g%ncpu)::ncore_cpu
  integer::i,j,icpu,ncore_loc,ilevel,igrid,ind,idim
  real(kind=8)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(kind=8)::mseed_core

  mcore_loc=0d0

  ! Seed mass in code units
  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  mseed_core = r%core_seed_mass*M_sun/(scale_d*scale_l**ndim)

  !---------------------------
  ! Flag core formation sites
  !---------------------------
  call select_core_candidates(r,g,m,c,mseed_core)

  !--------------------------
  ! Create new core particles
  !--------------------------
  ncore_loc=0
  do j=1,c%npeak
     if(c%form_sink(j).eq.1)then
        ilevel=c%peak_level(j)
        igrid=c%peak_grid(j)
        ind=c%peak_cell(j)
        ncore_loc=ncore_loc+1
        p%npart=p%npart+1
        if(p%npart>r%ncoremax)then
           write(*,*)'Error: not enough memory for core particles'
           write(*,*)'Increase ncoremax in the namelist'
           call mdl_abort(mdl)
        endif
        if(allocated(p%phip))p%phip(p%npart)=0.0
        ! Compute core particle position from peak position
        p%xp(p%npart,1)=c%peak_com(j,1)
        p%xp(p%npart,2)=c%peak_com(j,2)
        p%xp(p%npart,3)=c%peak_com(j,3)
        do idim=1,ndim
           if(r%periodic(idim))then
              if(p%xp(p%npart,idim)< 0.0d0           )p%xp(p%npart,idim)=p%xp(p%npart,idim)+r%box_size(idim)
              if(p%xp(p%npart,idim)>=r%box_size(idim))p%xp(p%npart,idim)=p%xp(p%npart,idim)-r%box_size(idim)
           endif
        end do
        ! Compute core particle velocity from peak velocity
        p%vp(p%npart,1)=c%peak_vel(j,1)-c%peak_acc(j,1)*0.5d0*g%dtnew(c%peak_level(j))
        p%vp(p%npart,2)=c%peak_vel(j,2)-c%peak_acc(j,2)*0.5d0*g%dtnew(c%peak_level(j))
        p%vp(p%npart,3)=c%peak_vel(j,3)-c%peak_acc(j,3)*0.5d0*g%dtnew(c%peak_level(j))
        ! Compute core particle old force from peak acceleration
        p%fp(p%npart,1)=c%peak_acc(j,1)
        p%fp(p%npart,2)=c%peak_acc(j,2)
        p%fp(p%npart,3)=c%peak_acc(j,3)
        ! Compute core particle mass, removed from the host cell
        p%mp(p%npart)=mseed_core
        call core_debit_host_cell(r,m,ilevel,igrid,ind,mseed_core, &
             p%zp(p%npart),p%up(p%npart))
        p%jp(p%npart,1:ndim)=0d0
        ! Compute core particle birth time using proper time
        p%tp(p%npart)=g%texp
        ! Compute level
        p%levelp(p%npart)=ilevel
        ! Count steps without accretion from the formation step
        p%nacc(p%npart)=nstep_coarse
     endif
  end do
  p%tailp(r%nlevelmax)=p%tailp(r%nlevelmax)+ncore_loc

  !---------------------------------------------
  ! Compute number of new cores across all CPUs
  !---------------------------------------------
  ncore_cpu=0
  ncore_cpu(g%myid)=ncore_loc
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(ncore_cpu,ncore_cpu_tot,g%ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  ncore_cpu=ncore_cpu_tot
#endif
  ncore_cum=0
  do icpu=1,g%ncpu
     ncore_cum(icpu)=ncore_cum(icpu-1)+int(ncore_cpu(icpu),kind=8)
  end do

  !---------------------------------
  ! Compute new core particle index
  !---------------------------------
  do i=1,ncore_loc
     p%idp(p%npart-ncore_loc+i)=p%npart_tot+ncore_cum(g%myid-1)+i
  end do
  p%npart_tot=p%npart_tot+ncore_cum(g%ncpu)
  mcore_loc=dble(ncore_loc)*mseed_core

end subroutine core_formation
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine core_debit_host_cell(r,m,ilevel,igrid,ind,mseed_core,z_init,u_init)
  use amr_parameters, only: ndim, dp
  use hydro_parameters, only: nvar, nener
  use amr_commons, only: run_t, mesh_t
  implicit none
  type(run_t)::r
  type(mesh_t)::m
  integer,intent(in)::ilevel,igrid,ind
  real(kind=8),intent(in)::mseed_core
  real(dp),intent(out)::z_init,u_init
  !------------------------------------------------------------------
  ! Remove the seed mass from the host cell, keeping velocity and specific
  ! thermal energy. Return the cell metallicity and thermal energy.
  !------------------------------------------------------------------
  real(kind=8)::dx,vol_loc,d,dnew,delta_d,e,ethermal,ekin
  real(kind=8),dimension(1:ndim)::vv
#ifdef MHD
  real(kind=8)::bx,by,bz,emag
#endif
#if NENER>0
  integer::irad
#endif
#if NVAR>5+NENER
  integer::ivar
  real(kind=8),dimension(1:nvar)::qpass
#endif

#ifdef HYDRO
  dx=r%boxlen/2d0**ilevel
  vol_loc=dx**ndim
  d=max(real(m%uold(ind,1,igrid),kind=8),r%smallr)
  delta_d=mseed_core/vol_loc
  dnew=d-delta_d
  vv(1:ndim)=m%uold(ind,2:ndim+1,igrid)/d
  e=m%uold(ind,5,igrid)/d
#ifdef MHD
  bx=0.5d0*(m%bold(ind,1,igrid)+m%bold(ind,4,igrid))
  by=0.5d0*(m%bold(ind,2,igrid)+m%bold(ind,5,igrid))
  bz=0.5d0*(m%bold(ind,3,igrid)+m%bold(ind,6,igrid))
  emag=0.5d0*(bx**2+by**2+bz**2)
  e=e-emag/d
#endif
#if NENER>0
  do irad=1,nener
     e=e-m%uold(ind,5+irad,igrid)/d
  end do
#endif
  ekin=0.5d0*sum(vv(1:ndim)**2)
  ethermal=max(e-ekin,r%smallc**2/(r%gamma-1d0))
  u_init=ethermal
  z_init=0d0
  if(r%metal)z_init=m%uold(ind,r%imetal,igrid)/d
#if NVAR>5+NENER
  do ivar=6+nener,nvar
     qpass(ivar)=m%uold(ind,ivar,igrid)/d
  end do
#endif

  ! Magnetic and non-thermal energies stay on the mesh
  m%uold(ind,1,igrid)=dnew
  m%uold(ind,2:ndim+1,igrid)=dnew*vv(1:ndim)
  m%uold(ind,5,igrid)=dnew*(ethermal+ekin)
#ifdef MHD
  m%uold(ind,5,igrid)=m%uold(ind,5,igrid)+emag
#endif
#if NENER>0
  do irad=1,nener
     m%uold(ind,5,igrid)=m%uold(ind,5,igrid)+m%uold(ind,5+irad,igrid)
  end do
#endif
#if NVAR>5+NENER
  do ivar=6+nener,nvar
     m%uold(ind,ivar,igrid)=dnew*qpass(ivar)
  end do
#endif
#endif

end subroutine core_debit_host_cell
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine select_core_candidates(r,g,m,c,mseed)
  use amr_parameters, only: ndim
  use amr_commons, only: run_t, global_t, mesh_t
  use clfind_commons, only: clump_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(clump_t)::c
  real(kind=8)::mseed
  integer::j,ilevel,igrid,ind
  logical::ok
  real(kind=8)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(kind=8)::dx,vol,mgas

  !-------------------------------------------------------------------
  ! Flag the peaks that form a core.
  !-------------------------------------------------------------------
  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  c%form_sink=0

  do j=1,c%npeak
     ilevel=c%peak_level(j)
     igrid=c%peak_grid(j)
     ind=c%peak_cell(j)
     dx=r%boxlen/2**ilevel
     vol=dx**ndim
     ok=.true.
     ! Without halo merging, each peak is its own halo
     if(c%saddle_threshold>0d0)then
        if(c%ind_halo(j).NE.j+c%npeak_cum(g%myid-1))ok=.false.
     endif
     if(c%relevance(j)<=c%relevance_threshold)ok=.false.
     if(c%ncore(j)>0)ok=.false.
     if(.not.core_candidate_valid_now(r,g,m,c,j))ok=.false.
#ifdef HYDRO
     mgas=m%uold(ind,1,igrid)*vol
     if(mseed>0.75*mgas)ok=.false.
#endif
     ! With LP the cell test replaces the fixed threshold
     if(.not.r%core_use_LP_criteria.and. &
          r%core_density_threshold>0d0)then
        if(c%max_dens(j)<=r%core_density_threshold/scale_nH)ok=.false.
     endif
     if(r%core_use_clump_prop_criteria.or.r%core_use_contracting_all_axes)then
        if(r%core_use_mass4_criteria)then
           if(c%clump_mass4(j)<=r%core_mass4_threshold)ok=.false.
        endif
        if(r%core_use_energy_criteria)then
           if(.not.c%bound(j))ok=.false.
        endif
        if(r%core_use_contracting_criteria.or.r%core_use_contracting_all_axes)then
           if(.not.c%contracting(j))ok=.false.
        endif
     endif
     if(ok)c%form_sink(j)=1
  end do

end subroutine select_core_candidates
!###########################################################
!###########################################################
!###########################################################
!###########################################################
logical function core_candidate_valid_now(r,g,m,c,ipeak) result(ok)
  use amr_parameters, only: ndim
#if NENER>0
  use hydro_parameters, only: nener
#endif
  use amr_commons, only: run_t, global_t, mesh_t
  use clfind_commons, only: clump_t
  use constants, only: pi
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(clump_t)::c
  integer::ipeak
  !------------------------------------------------------------------
  ! Density test of the peak cell at formation (LP or fixed threshold).
  !------------------------------------------------------------------
  integer::ilevel,igrid,ind
  real(kind=8)::dx,dens
  real(kind=8)::e,ekin,ethermal,c_s2,d_core_LP,factG
  real(kind=8)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(kind=8),dimension(1:ndim)::vv
#ifdef MHD
  real(kind=8)::bx,by,bz,emag
#endif
#if NENER>0
  integer::irad
#endif

  ok=.true.
#ifdef HYDRO
  if(r%rho_type_core.ne.4)return
  ilevel=c%peak_level(ipeak)
  igrid=c%peak_grid(ipeak)
  ind=c%peak_cell(ipeak)
  dx=r%boxlen/2d0**ilevel
  dens=m%uold(ind,1,igrid)
  if(r%core_use_LP_criteria)then
     ! Larson-Penston density for this cell size and sound speed
     vv(1:ndim)=m%uold(ind,2:ndim+1,igrid)/max(dens,r%smallr)
     e=m%uold(ind,5,igrid)/max(dens,r%smallr)
#ifdef MHD
     bx=0.5d0*(m%bold(ind,1,igrid)+m%bold(ind,4,igrid))
     by=0.5d0*(m%bold(ind,2,igrid)+m%bold(ind,5,igrid))
     bz=0.5d0*(m%bold(ind,3,igrid)+m%bold(ind,6,igrid))
     emag=0.5d0*(bx**2+by**2+bz**2)
     e=e-emag/max(dens,r%smallr)
#endif
#if NENER>0
     do irad=1,nener
        e=e-m%uold(ind,5+irad,igrid)/max(dens,r%smallr)
     end do
#endif
     ekin=0.5d0*sum(vv(1:ndim)**2)
     ethermal=max(e-ekin,r%smallc**2/(r%gamma-1d0))
     c_s2=r%gamma*(r%gamma-1d0)*ethermal
     factG=1d0
     if(r%cosmo)factG=3d0/8d0/pi*g%omega_m*g%aexp
     d_core_LP=8.86d0*c_s2/pi/factG/dx**2
     ok=dens>d_core_LP
  else if(r%core_density_threshold>0d0)then
     call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
     ok=dens>r%core_density_threshold/scale_nH
  endif
#endif

end function core_candidate_valid_now
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine count_existing_cores_in_peak(s)
  use ramses_commons, only: ramses_t
  use clump_merger_module, only: get_peak
  use cache_commons, only: msg_sink_clump
  use cache
  implicit none
  type(ramses_t)::s

  integer::i,ipart,peak_nr,ilev,ipeak,jpeak
  integer(kind=8)::global_peak_id,merge_to
  type(msg_sink_clump)::dummy_sink_clump

  associate(c=>s%c,g=>s%g,p=>s%core,mdl=>s%mdl)

  ! Count cores in each peak
  c%ncore=0
  do ipart=1,p%npart
     p%sortp(ipart)=ipart
     p%workp(ipart)=p%pid(ipart)
  end do
  call quick_sort_int_int(p%workp(1),p%sortp(1),p%npart)

  call open_cache_clump(mdl, c, pack_size=storage_size(dummy_sink_clump)/32, &
       init=init_flush_core, flush=pack_flush_core, combine=unpack_flush_core)
  do i=1+p%norphan_peak,p%npart
     if(p%mp(p%sortp(i))<=0d0)cycle
     global_peak_id=p%workp(i)
     call get_peak(s,global_peak_id,peak_nr,fetch_cache=.false.,flush_cache=.true.)
     c%ncore(peak_nr)=c%ncore(peak_nr)+1
  end do
  call close_cache(mdl)

  ! Count cores hierarchically in each halo
  if(c%saddle_threshold>0)then
     do ilev=0,c%merge_levelmax
        call open_cache_clump(mdl, c, pack_size=storage_size(dummy_sink_clump)/32, &
             init=init_flush_core, flush=pack_flush_core, combine=unpack_flush_core)
        do ipeak=1,c%npeak
           if(c%lev_peak(ipeak)==ilev)then
              merge_to=c%new_peak(ipeak)
              if(merge_to.NE.ipeak+c%npeak_cum(g%myid-1))then
                 call get_peak(s,merge_to,jpeak,flush_cache=.true.,fetch_cache=.false.)
                 c%ncore(jpeak)=c%ncore(jpeak)+c%ncore(ipeak)
              endif
           endif
        end do
        call close_cache(mdl)
     end do
  endif

  end associate

end subroutine count_existing_cores_in_peak
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine init_flush_core(c,local_peak_id)
  use clfind_commons, only: clump_t
  type(clump_t)::c
  integer::local_peak_id

  c%ncore(local_peak_id)=0

end subroutine init_flush_core
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine pack_flush_core(c,local_peak_id,msg_size,msg_array)
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_sink_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_sink_clump)::msg

  msg%lev=c%ncore(local_peak_id)
  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_core
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine unpack_flush_core(c,local_peak_id,msg_size,msg_array)
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_sink_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_sink_clump)::msg

  msg=transfer(msg_array,msg)
  c%ncore(local_peak_id)=c%ncore(local_peak_id)+msg%lev

end subroutine unpack_flush_core
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine compute_core_clump_criteria(s,rtype)
  use amr_parameters, only: ndim
  use hydro_parameters, only: nener
  use ramses_commons, only: ramses_t
  use clump_merger_module, only: get_peak, pack_fetch_part, unpack_fetch_part
  use cache_commons, only: msg_prop_clump
  use cache
#ifdef DO_RT
  use rt_parameters, only: nrtgrp
  use constants, only: eV2erg, c_cgs
#endif
  implicit none
  type(ramses_t)::s
  integer::rtype
  !------------------------------------------------------------------
  ! Virial and contraction criteria of each gas clump, from its test cells.
  !------------------------------------------------------------------
  type(msg_prop_clump)::dummy_prop_clump
  type(msg_core_clump)::dummy_core_clump
  integer(kind=8)::global_peak_id
  integer::itest,ipeak,peak_nr,ilevel,igrid,ind,idim,jdim,iaxis,irad
  real(kind=8)::dx_loc,dx_min,vol,d,etot,ekin,emag,erad,pres,pmag,support
  real(kind=8),dimension(1:ndim)::xcell,vel,vrel,rrel,fgrav,xpeak
  real(kind=8),dimension(1:3)::axis_rate
  real(kind=8),dimension(1:3,1:3)::inertia_mat,idot_mat,eigvec
  logical::all_axes_contracting
#ifdef DO_RT
  integer::ig,iN
  real(kind=8),dimension(1:ndim)::frad
  real(kind=8)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(kind=8)::scale_Np,scale_Fp,scale_kappa,c_code,kappa,Np2Ep_flux
#endif
#ifdef MHD
  real(kind=8)::bx,by,bz
#endif

#ifdef HYDRO
  associate(r=>s%r,g=>s%g,m=>s%m,c=>s%c,mdl=>s%mdl)

  if(rtype.ne.4)return

  allocate(gas_mass(1:c%npeak_max),gas_idot(1:c%npeak_max))
  allocate(gas_com(1:c%npeak_max,1:ndim),gas_vel(1:c%npeak_max,1:ndim))
  allocate(gas_inertia(1:c%npeak_max,1:ndim,1:ndim),gas_idot_tensor(1:c%npeak_max,1:ndim,1:ndim))

  !-----------------------------------------------------------------
  ! Gas mass, centre of mass and velocity of each peak patch, with
  ! the peak position fetched for the periodic image
  !-----------------------------------------------------------------
  gas_mass=0d0; gas_com=0d0; gas_vel=0d0
  call open_cache_clump(mdl, c, storage_size(dummy_prop_clump)/32, &
       pack=pack_fetch_part, unpack=unpack_fetch_part, &
       init=init_flush_core_mass, flush=pack_flush_core_mass, combine=unpack_flush_core_mass)
  do itest=1,c%ntest
     ilevel=c%level(itest); igrid=c%grid(itest); ind=c%cell(itest)
     global_peak_id=m%flag2(ind,igrid)
     if(global_peak_id==0)cycle
     call get_peak(s,global_peak_id,peak_nr,flush_cache=.true.,fetch_cache=.true.)
     dx_loc=r%boxlen/2d0**ilevel
     vol=dx_loc**ndim
     xpeak(1:ndim)=c%peak_pos(peak_nr,1:ndim)
     call core_clump_cell_position(r,m,igrid,ind,dx_loc,xpeak,xcell)
     d=max(real(m%uold(ind,1,igrid),kind=8),r%smallr)
     vel(1:ndim)=m%uold(ind,2:ndim+1,igrid)/d
     gas_mass(peak_nr)=gas_mass(peak_nr)+d*vol
     gas_com(peak_nr,1:ndim)=gas_com(peak_nr,1:ndim)+d*vol*xcell(1:ndim)
     gas_vel(peak_nr,1:ndim)=gas_vel(peak_nr,1:ndim)+d*vol*vel(1:ndim)
  end do
  call close_cache(mdl)

  do ipeak=1,c%npeak
     if(gas_mass(ipeak)>0d0)then
        gas_com(ipeak,1:ndim)=gas_com(ipeak,1:ndim)/gas_mass(ipeak)
        gas_vel(ipeak,1:ndim)=gas_vel(ipeak,1:ndim)/gas_mass(ipeak)
     endif
  end do

  !-----------------------------------------------------------------
  ! Virial terms, dI/dt and inertia tensors about the centre of mass,
  ! with the peak position and centre of mass fetched
  !-----------------------------------------------------------------
  c%clump_mass4=0d0; c%kinetic_support=0d0; c%thermal_support=0d0
  c%magnetic_support=0d0; c%grav_term=0d0; c%rad_term=0d0
  gas_idot=0d0; gas_inertia=0d0; gas_idot_tensor=0d0
  dx_min=r%boxlen/2d0**r%nlevelmax
#ifdef DO_RT
  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  call rt_units(r,g,scale_Np,scale_Fp)
  scale_kappa=(scale_d*scale_l)**(-1d0)
  c_code=c_cgs/scale_v
#endif

  call open_cache_clump(mdl, c, storage_size(dummy_core_clump)/32, &
       pack=pack_fetch_core_com, unpack=unpack_fetch_core_com, &
       init=init_flush_core_virial, flush=pack_flush_core_virial, combine=unpack_flush_core_virial)
  do itest=1,c%ntest
     ilevel=c%level(itest); igrid=c%grid(itest); ind=c%cell(itest)
     global_peak_id=m%flag2(ind,igrid)
     if(global_peak_id==0)cycle
     call get_peak(s,global_peak_id,peak_nr,flush_cache=.true.,fetch_cache=.true.)
     dx_loc=r%boxlen/2d0**ilevel
     vol=dx_loc**ndim
     xpeak(1:ndim)=c%peak_pos(peak_nr,1:ndim)
     call core_clump_cell_position(r,m,igrid,ind,dx_loc,xpeak,xcell)
     d=max(real(m%uold(ind,1,igrid),kind=8),r%smallr)
     vel(1:ndim)=m%uold(ind,2:ndim+1,igrid)/d
     etot=m%uold(ind,5,igrid)
     ekin=0d0
     do idim=1,ndim
        ekin=ekin+0.5d0*d*vel(idim)**2
     end do
     emag=0d0; pmag=0d0
#ifdef MHD
     bx=0.5d0*(m%bold(ind,1,igrid)+m%bold(ind,4,igrid))
     by=0.5d0*(m%bold(ind,2,igrid)+m%bold(ind,5,igrid))
     bz=0.5d0*(m%bold(ind,3,igrid)+m%bold(ind,6,igrid))
     emag=0.5d0*(bx*bx+by*by+bz*bz)
     pmag=emag
#endif
     erad=0d0
     do irad=1,nener
        erad=erad+m%uold(ind,5+irad,igrid)
     end do
     pres=(etot-ekin-emag-erad)*(r%gamma-1d0)+erad/3d0
     rrel(1:ndim)=xcell(1:ndim)-gas_com(peak_nr,1:ndim)
     vrel(1:ndim)=vel(1:ndim)-gas_vel(peak_nr,1:ndim)
     if(sum((xcell(1:ndim)-xpeak(1:ndim))**2) <= (4d0*dx_min)**2)then
        c%clump_mass4(peak_nr)=c%clump_mass4(peak_nr)+d*vol
     endif
     c%kinetic_support(peak_nr)=c%kinetic_support(peak_nr)+sum(vrel(1:ndim)**2)*vol*d
     c%thermal_support(peak_nr)=c%thermal_support(peak_nr)+3d0*vol*pres
     gas_idot(peak_nr)=gas_idot(peak_nr)+2d0*sum(rrel(1:ndim)*vrel(1:ndim))*vol*d
     do idim=1,ndim
        do jdim=1,ndim
           gas_inertia(peak_nr,idim,jdim)=gas_inertia(peak_nr,idim,jdim)+rrel(idim)*rrel(jdim)*vol*d
           gas_idot_tensor(peak_nr,idim,jdim)=gas_idot_tensor(peak_nr,idim,jdim)+ &
                (vrel(jdim)*rrel(idim)+vrel(idim)*rrel(jdim))*vol*d
        end do
     end do
#ifdef GRAV
     fgrav(1:ndim)=m%f(ind,1:ndim,igrid)
     c%grav_term(peak_nr)=c%grav_term(peak_nr)+sum(fgrav(1:ndim)*rrel(1:ndim))*vol*d
#endif
#ifdef DO_RT
     frad=0d0
     do ig=1,nrtgrp
        iN=1+(ig-1)*(ndim+1)
        Np2Ep_flux=(scale_Fp*r%group_egy(ig)*eV2erg)/(scale_d*scale_v**3)
        kappa=r%kappaSc(ig)/scale_kappa
        frad(1:ndim)=frad(1:ndim)+m%rtuold(ind,iN+1:iN+ndim,igrid)*Np2Ep_flux*kappa/c_code
     end do
     c%rad_term(peak_nr)=c%rad_term(peak_nr)+sum(frad(1:ndim)*rrel(1:ndim))*vol*d
#endif
     c%magnetic_support(peak_nr)=c%magnetic_support(peak_nr)+3d0*vol*pmag
  end do
  call close_cache(mdl)

  !-----------------------------------------------------------------
  ! Virial ratio, bound and contracting flags of each local peak
  !-----------------------------------------------------------------
  do ipeak=1,c%npeak
     support=c%kinetic_support(ipeak)+c%thermal_support(ipeak)
#ifdef DO_RT
     if(r%core_include_radiation_support)support=support+c%rad_term(ipeak)
#endif
#ifdef MHD
     if(r%core_include_magnetic_support)support=support+c%magnetic_support(ipeak)
#endif
     if(c%grav_term(ipeak)<0d0)then
        c%virial_ratio(ipeak)=support/abs(c%grav_term(ipeak))
        ! Bound if kinetic and thermal terms are each below abs(W)
        c%bound(ipeak)=c%kinetic_support(ipeak)<abs(c%grav_term(ipeak)).and. &
             c%thermal_support(ipeak)<abs(c%grav_term(ipeak))
     else
        c%virial_ratio(ipeak)=huge(0d0)
        c%bound(ipeak)=.false.
     endif
     if(r%core_use_contracting_all_axes)then
        inertia_mat=0d0; idot_mat=0d0; axis_rate=0d0
        inertia_mat(1:ndim,1:ndim)=gas_inertia(ipeak,1:ndim,1:ndim)
        idot_mat(1:ndim,1:ndim)=gas_idot_tensor(ipeak,1:ndim,1:ndim)
        call jacobi_sym3(inertia_mat,eigvec)
        do iaxis=1,ndim
           do idim=1,ndim
              do jdim=1,ndim
                 axis_rate(iaxis)=axis_rate(iaxis)+idot_mat(idim,jdim)*eigvec(idim,iaxis)*eigvec(jdim,iaxis)
              end do
           end do
        end do
        all_axes_contracting=c%n_cells(ipeak)>1
        do iaxis=1,ndim
           all_axes_contracting=all_axes_contracting.and.axis_rate(iaxis)<0d0
        end do
        c%contracting(ipeak)=all_axes_contracting
     else
        c%contracting(ipeak)=gas_idot(ipeak)<0d0
     endif
  end do

  deallocate(gas_mass,gas_idot,gas_com,gas_vel,gas_inertia,gas_idot_tensor)

  end associate
#endif

contains
  subroutine jacobi_sym3(a,v)
    implicit none
    real(kind=8),dimension(1:3,1:3)::a,v
    integer::iter,pidx,qidx,j
    real(kind=8)::app,aqq,apq,tau,t,c,s,tmp,offmax

    v=0d0
    do j=1,3
       v(j,j)=1d0
    end do
    do iter=1,50
       pidx=1; qidx=2; offmax=abs(a(1,2))
       if(abs(a(1,3))>offmax)then
          pidx=1; qidx=3; offmax=abs(a(1,3))
       endif
       if(abs(a(2,3))>offmax)then
          pidx=2; qidx=3; offmax=abs(a(2,3))
       endif
       if(offmax==0d0)exit
       app=a(pidx,pidx); aqq=a(qidx,qidx); apq=a(pidx,qidx)
       tau=(aqq-app)/(2d0*apq)
       if(tau>=0d0)then
          t=1d0/(tau+sqrt(1d0+tau*tau))
       else
          t=-1d0/(-tau+sqrt(1d0+tau*tau))
       endif
       c=1d0/sqrt(1d0+t*t); s=t*c
       do j=1,3
          if(j/=pidx.and.j/=qidx)then
             tmp=a(j,pidx)
             a(j,pidx)=c*tmp-s*a(j,qidx)
             a(pidx,j)=a(j,pidx)
             a(j,qidx)=s*tmp+c*a(j,qidx)
             a(qidx,j)=a(j,qidx)
          endif
       end do
       a(pidx,pidx)=c*c*app-2d0*s*c*apq+s*s*aqq
       a(qidx,qidx)=s*s*app+2d0*s*c*apq+c*c*aqq
       a(pidx,qidx)=0d0; a(qidx,pidx)=0d0
       do j=1,3
          tmp=v(j,pidx)
          v(j,pidx)=c*tmp-s*v(j,qidx)
          v(j,qidx)=s*tmp+c*v(j,qidx)
       end do
    end do
  end subroutine jacobi_sym3

  subroutine core_clump_cell_position(r,m,igrid,ind,dx_loc,xpeak,xcell)
    use amr_parameters, only: ndim
    use amr_commons, only: run_t, mesh_t
    implicit none
    type(run_t)::r
    type(mesh_t)::m
    integer::igrid,ind,nstride,idim
    real(kind=8)::dx_loc
    real(kind=8),dimension(1:ndim)::xpeak,xcell
    ! Cell centre, in the periodic image closest to the peak
    do idim=1,ndim
       nstride=2**(idim-1)
       xcell(idim)=(2*m%grid(igrid)%ckey(idim)+MOD((ind-1)/nstride,2)+0.5d0)*dx_loc-m%skip(idim)
       if(r%periodic(idim))then
          if(xcell(idim)-xpeak(idim)> r%box_size(idim)*0.5)xcell(idim)=xcell(idim)-r%box_size(idim)
          if(xcell(idim)-xpeak(idim)<-r%box_size(idim)*0.5)xcell(idim)=xcell(idim)+r%box_size(idim)
       endif
    end do
  end subroutine core_clump_cell_position

end subroutine compute_core_clump_criteria
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine init_flush_core_mass(c,local_peak_id)
  use clfind_commons, only: clump_t
  type(clump_t)::c
  integer::local_peak_id

  gas_mass(local_peak_id)=0d0
  gas_com(local_peak_id,1:ndim)=0d0
  gas_vel(local_peak_id,1:ndim)=0d0

end subroutine init_flush_core_mass
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine pack_flush_core_mass(c,local_peak_id,msg_size,msg_array)
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_prop_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_prop_clump)::msg

  msg%mass=gas_mass(local_peak_id)
  msg%pos(1:ndim)=gas_com(local_peak_id,1:ndim)
  msg%vel(1:ndim)=gas_vel(local_peak_id,1:ndim)

  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_core_mass
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine unpack_flush_core_mass(c,local_peak_id,msg_size,msg_array)
  use clfind_commons, only: clump_t
  use cache_commons, only: msg_prop_clump
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_prop_clump)::msg

  msg=transfer(msg_array,msg)

  gas_mass(local_peak_id)=gas_mass(local_peak_id)+msg%mass
  gas_com(local_peak_id,1:ndim)=gas_com(local_peak_id,1:ndim)+msg%pos(1:ndim)
  gas_vel(local_peak_id,1:ndim)=gas_vel(local_peak_id,1:ndim)+msg%vel(1:ndim)

end subroutine unpack_flush_core_mass
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine pack_fetch_core_com(c,local_peak_id,msg_size,msg_array)
  use clfind_commons, only: clump_t
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_core_clump)::msg

  msg%pos(1:ndim)=c%peak_pos(local_peak_id,1:ndim)
  msg%com(1:ndim)=gas_com(local_peak_id,1:ndim)
  msg%vel(1:ndim)=gas_vel(local_peak_id,1:ndim)

  msg_array=transfer(msg,msg_array)

end subroutine pack_fetch_core_com
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine unpack_fetch_core_com(c,local_peak_id,msg_size,msg_array)
  use clfind_commons, only: clump_t
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_core_clump)::msg

  msg=transfer(msg_array,msg)

  c%peak_pos(local_peak_id,1:ndim)=msg%pos(1:ndim)
  gas_com(local_peak_id,1:ndim)=msg%com(1:ndim)
  gas_vel(local_peak_id,1:ndim)=msg%vel(1:ndim)

end subroutine unpack_fetch_core_com
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine init_flush_core_virial(c,local_peak_id)
  use clfind_commons, only: clump_t
  type(clump_t)::c
  integer::local_peak_id

  c%clump_mass4(local_peak_id)=0d0
  c%kinetic_support(local_peak_id)=0d0
  c%thermal_support(local_peak_id)=0d0
  c%magnetic_support(local_peak_id)=0d0
  c%grav_term(local_peak_id)=0d0
  c%rad_term(local_peak_id)=0d0
  gas_idot(local_peak_id)=0d0
  gas_inertia(local_peak_id,1:ndim,1:ndim)=0d0
  gas_idot_tensor(local_peak_id,1:ndim,1:ndim)=0d0

end subroutine init_flush_core_virial
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine pack_flush_core_virial(c,local_peak_id,msg_size,msg_array)
  use clfind_commons, only: clump_t
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_core_clump)::msg

  msg%mass4=c%clump_mass4(local_peak_id)
  msg%kin=c%kinetic_support(local_peak_id)
  msg%therm=c%thermal_support(local_peak_id)
  msg%mag=c%magnetic_support(local_peak_id)
  msg%grav=c%grav_term(local_peak_id)
  msg%rad=c%rad_term(local_peak_id)
  msg%idot=gas_idot(local_peak_id)
  msg%inertia(1:ndim,1:ndim)=gas_inertia(local_peak_id,1:ndim,1:ndim)
  msg%idot_tensor(1:ndim,1:ndim)=gas_idot_tensor(local_peak_id,1:ndim,1:ndim)

  msg_array=transfer(msg,msg_array)

end subroutine pack_flush_core_virial
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine unpack_flush_core_virial(c,local_peak_id,msg_size,msg_array)
  use clfind_commons, only: clump_t
  type(clump_t)::c
  integer::local_peak_id
  integer::msg_size
  integer,dimension(1:msg_size),optional::msg_array

  type(msg_core_clump)::msg

  msg=transfer(msg_array,msg)

  c%clump_mass4(local_peak_id)=c%clump_mass4(local_peak_id)+msg%mass4
  c%kinetic_support(local_peak_id)=c%kinetic_support(local_peak_id)+msg%kin
  c%thermal_support(local_peak_id)=c%thermal_support(local_peak_id)+msg%therm
  c%magnetic_support(local_peak_id)=c%magnetic_support(local_peak_id)+msg%mag
  c%grav_term(local_peak_id)=c%grav_term(local_peak_id)+msg%grav
  c%rad_term(local_peak_id)=c%rad_term(local_peak_id)+msg%rad
  gas_idot(local_peak_id)=gas_idot(local_peak_id)+msg%idot
  gas_inertia(local_peak_id,1:ndim,1:ndim)=gas_inertia(local_peak_id,1:ndim,1:ndim)+msg%inertia(1:ndim,1:ndim)
  gas_idot_tensor(local_peak_id,1:ndim,1:ndim)=gas_idot_tensor(local_peak_id,1:ndim,1:ndim)+msg%idot_tensor(1:ndim,1:ndim)

end subroutine unpack_flush_core_virial

end module core_formation_module
