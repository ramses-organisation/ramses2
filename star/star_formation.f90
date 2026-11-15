module star_formation_module

  type :: out_star_formation_t
     real(kind=8)::mass
  end type out_star_formation_t

contains
!###########################################################
!###########################################################
!###########################################################
!###########################################################
recursive subroutine r_star_formation(pst,ilevel,input_size,output,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  use SED_module,only: update_SED_group_props
#ifdef _CUDA
  use gpu_runner, only: gpu_star_formation
#endif
#if defined(_CUDA) && defined(DO_RT)
  use gpu_manager, only: gpu_star_to_host
#endif
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  type(out_star_formation_t)::output,next_output

  integer::ilevel
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_STAR_FORMATION,pst%iUpper+1,input_size,output_size,ilevel)
     call r_star_formation(pst%pLower,ilevel,input_size,output,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_output)
     output%mass=output%mass+next_output%mass
  else
#ifdef _CUDA
     call gpu_star_formation(pst%s,ilevel,output%mass)
#else
     call star_formation(pst%s%r,pst%s%g,pst%s%m,pst%s%star,ilevel,output%mass)
#endif

     ! Check if radiadion advection should be turned on
     if(pst%s%r%rt .and. pst%s%r%rt_star .and. .not. pst%s%r%rt_advect .and. pst%s%star%npart_tot .gt. 0) then
        if(pst%s%g%myid==1) then
          write(*,*) '*****************************************'
          if(pst%s%r%cosmo) then
            write(*,*) 'Stellar RT turned on at a=',pst%s%g%aexp
          else
            write(*,*) 'Stellar RT turned on at t=',pst%s%g%t
          endif
          write(*,*) '*****************************************'
        endif
        pst%s%r%rt_advect=.true.
        ! Update cross sections based on star properties
        if(pst%s%r%sedprops_update .gt. 0) then
#if defined(_CUDA) && defined(DO_RT)
          call gpu_star_to_host(pst%s)
#endif
          call update_SED_group_props(pst%s%r, pst%s%g, pst%s%sed, pst%s%star)
        endif
     endif

     ! Check if CR advection should be turned on
     if(pst%s%r%cr .and. .not. pst%s%r%cr_advect .and. pst%s%star%npart_tot .gt. 0) then
        if(pst%s%g%myid==1) then
          write(*,*) '*****************************************'
          if(pst%s%r%cosmo) then
            write(*,*) 'CR advection turned on at a=',pst%s%g%aexp
          else
            write(*,*) 'CR advection turned on at t=',pst%s%g%t
          endif
          write(*,*) '*****************************************'
        endif
        pst%s%r%cr_advect=.true.
     endif

  endif

end subroutine r_star_formation
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine star_formation(r,g,m,s,ilevel,mstar_loc)
  use rng
  use constants
  use hydro_parameters, only: nvar
  use amr_parameters, only: ndim, twotondim
  use amr_commons, only: run_t, global_t, mesh_t
  use pm_commons, only: part_t
  use input_hydro_condinit_module, only: cons_from_prim, prim_from_cons
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(part_t)::s
  integer::ilevel
  real(kind=8)::mstar_loc
  !-------------------------------------------------------------------
  ! Spawn star particles according to various star formation models.
  ! We use a random Poisson process.
  !-------------------------------------------------------------------
#ifndef WITHOUTMPI
  integer::info
  integer,dimension(1:g%ncpu)::nsite_cpu_tot,nstar_cpu_tot
#endif
  integer(kind=8),dimension(0:g%ncpu)::nsite_cum,nstar_cum
  integer,dimension(1:g%ncpu)::nsite_cpu,nstar_cpu
  integer::i,ind,igrid,idim,icpu,ngrid,nleaf,nsite,nstar,nstar_loc
  integer(kind=8)::nsite8
  real(kind=8)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,erfc
  real(kind=8)::dx,vol,factG,n_star,nCOM,d,d0,mstar,dstar,t_ff,mcell,mgas,mask,PoissMean,Rand
  real(kind=8)::sfr_ff,t_dyn,p,cs2,sigma2,alpha_vir,sigs,scrit,Mach2,b_turb
#if NENER>0
  integer::irad
#endif
  type(RngStream)::RngStream_CreateStream, gg
  real(kind=8)::RngStream_RandUni
  logical::ok

#if NDIM>2
#ifdef HYDRO
  ! Conversion factor from user units to cgs units
  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Set some constants
  dx=r%boxlen/2**ilevel
  vol=dx**ndim
  factG=1d0
  if(r%cosmo)factG=3d0/4d0/twopi*g%omega_m*g%aexp

  ! Density threshold for star formation
  n_star = r%n_star

  ! This is used only for cosmological runs
  if(r%cosmo)then
     nCOM = 200.0*g%omega_b*rhoc*(g%h0/100)**2/g%aexp**3/mH
     n_star = MAX(nCOM,n_star)
  endif
  d0 = n_star/scale_nH
  mstar = r%m_star*r%mass_sph
  dstar = mstar/vol

  !---------------------------------------------------------
  ! Convert conservative variables to primitive variables
  !---------------------------------------------------------
  call prim_from_cons(r,g,m,ilevel)

  !---------------------------------------------------------
  ! Count potential star formation sites.
  ! This is important for the pseudo-random number generator
  !---------------------------------------------------------
  nsite8=0
  ! Loop over octs with vector sweeps
  do igrid=m%head(ilevel),m%tail(ilevel)
     ! Loop over cells
     do ind=1,twotondim
        ! Select leaf cells
        ok = .not. m%grid(igrid)%refined(ind)
        ! Select dense enough cells
        d = m%uold(ind,1,igrid)
        ok = ok .and. d > d0
        ! Select cells in zoom region
        if(r%ivar_refine>0)then
           mask = m%uold(ind,r%ivar_refine,igrid)
           ok = ok .and. mask > r%var_cut_refine
        endif
        ! Count number of random numbers
        if(ok)then
           nsite8=nsite8+1
        endif
     end do
  end do
  if(nsite8 > m%ngridmax)then
     write(*,*)'star_formation: cells above star formation threshold (',nsite8,') exceed ngridmax (',m%ngridmax,')'
     write(*,*)'Increase n_star in the namelist'
     stop
  endif
  nsite = int(nsite8, kind=4)

  !---------------------------------------------------------
  ! Compute number of star formation sites across all CPUs.
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

  !------------------------------------------------
  ! Initialize RNG stream based on the global seed
  !------------------------------------------------
  call RngStream_SetPackageSeed(r%seed)
  ! Get current rng stream state
  gg = RngStream_CreateStream('ramses_stream')
  ! Advance to the state corresponding to the current cpu
  if(g%myid>1)then
     if(nsite_cum(g%myid-1)>0)then
        call RngStream_AdvanceState(gg,0_8,2*nsite_cum(g%myid-1))
     endif
  endif

  !------------------------------------------------
  ! Create new stars based on local SF efficiency
  !------------------------------------------------
  nstar_loc=0
  mstar_loc=0.0d0
  ! Loop over octs
  do igrid=m%head(ilevel),m%tail(ilevel)
     ! Loop over cells
     do ind=1,twotondim
        ! Select leaf cells
        ok = .not. m%grid(igrid)%refined(ind)
        ! Select dense enough cells
        d = m%uold(ind,1,igrid)
        ok = ok .and. d > d0
        ! Select cells in zoom region
        if(r%ivar_refine>0)then
           mask = m%uold(ind,r%ivar_refine,igrid)
           ok = ok .and. mask > r%var_cut_refine
        endif
        ! Compute sound speed squared
        p = m%uold(ind,5,igrid)
        cs2 = max(dble(r%gamma)*p/d,dble(r%smallc)**2)
        ! Turbulence 1D velocity dispersion
        sigma2 = 0d0
        if(r%sgs_turb)then
           sigma2 = m%uold(ind,r%iturb,igrid)*2d0/3d0
        endif

        ! Draw Poisson process based on local SF rate
        if(ok)then

           ! Compute local SF efficiency
           select case (r%sf_model)

           case(1)
              ! Constant efficiency
              sfr_ff = r%eps_star

           case(2)
              ! Padoan, Haugbolle & Nordlund 2012
              t_dyn = dx/(2d0*sqrt(3d0*(sigma2+cs2)))
              t_ff = 0.5427d0/sqrt(factG*d)
              sfr_ff = r%eps_star*exp(-1.6d0*t_ff/t_dyn)

           case(3)
              ! Multi-freefall model a la Krumholz & McKee
              alpha_vir = (15d0*(sigma2+cs2))/(pi*factG*d*dx**2)
              Mach2 = max(sigma2/cs2,dble(r%smallr))
              b_turb = 1.0 ! Turbulent forcing parameter (Federrath 2008 & 2010)
              sigs = log(1d0+b_turb**2*Mach2)
              scrit = log(alpha_vir*(1d0+(2d0*Mach2**2/(1d0+Mach2))))
              sfr_ff = (r%eps_star/2d0)*exp(3d0/8d0*sigs)*(2d0-erfc((sigs-scrit)/sqrt(2d0*sigs)))

           end select
           ! Compute stellar mass formed
           mcell=d*vol
           t_ff=0.5427d0/sqrt(factG*d)
           mgas=g%dtnew(ilevel)*(sfr_ff/t_ff)*mcell
           ! Poisson mean
           PoissMean=mgas/mstar
           call poissdev(RngStream_RandUni(gg),PoissMean,nstar)
           ! Draw a second random number for the formation age
           Rand=RngStream_RandUni(gg)
           ! Check that we don't empty the gas cell
           mgas=nstar*mstar
           if(mgas>0.9*mcell)then
              nstar=int(0.9*mcell/mstar)
           endif
           ! Create a new star particle
           if(nstar>0)then
              nstar_loc=nstar_loc+1
              mstar_loc=mstar_loc+nstar*mstar
              dstar=nstar*mstar/vol
              s%npart=s%npart+1
              if(s%npart>r%nstarmax)then
                 write(*,*)'Not enough memory for star particle'
                 write(*,*)'Increase nstarmax in the namelist'
                 stop
              endif
              if(allocated(s%phip))s%phip(s%npart)=0.0
              ! Compute star particle coordinate from cell centers
              s%xp(s%npart,1)=(2*m%grid(igrid)%ckey(1)+MOD((ind-1)  ,2)+0.5)*dx-m%skip(1)
              s%xp(s%npart,2)=(2*m%grid(igrid)%ckey(2)+MOD((ind-1)/2,2)+0.5)*dx-m%skip(2)
              s%xp(s%npart,3)=(2*m%grid(igrid)%ckey(3)+MOD((ind-1)/4,2)+0.5)*dx-m%skip(3)
              ! Compute star particle velocity from gas velocity
              s%vp(s%npart,1)=m%uold(ind,2,igrid)
              s%vp(s%npart,2)=m%uold(ind,3,igrid)
              s%vp(s%npart,3)=m%uold(ind,4,igrid)
#ifdef GRAV
              ! Remove half a kick (will be added later)
              s%vp(s%npart,1)=s%vp(s%npart,1)-m%f(ind,1,igrid)*0.5d0*g%dtnew(ilevel)
              s%vp(s%npart,2)=s%vp(s%npart,2)-m%f(ind,2,igrid)*0.5d0*g%dtnew(ilevel)
              s%vp(s%npart,3)=s%vp(s%npart,3)-m%f(ind,3,igrid)*0.5d0*g%dtnew(ilevel)
#endif
              ! Compute star particle mass
              s%mp(s%npart)=nstar*mstar
              ! Compute star particle birth time using proper time
              s%tp(s%npart)=g%texp+Rand*g%dtnew(ilevel)*g%aexp**2
              ! Compute star particle metallicity
              if(r%metal)then
                 s%zp(s%npart)=m%uold(ind,r%imetal,igrid)
              else
                 s%zp(s%npart)=r%z_ave*0.02
              endif
              ! Compute level
              s%levelp(s%npart)=ilevel
              ! Remove corresponding mass from gas mass density
              m%uold(ind,1,igrid)=d-dstar
              ! Star formation is isothermal, adjust pressure
              m%uold(ind,5,igrid)=m%uold(ind,1,igrid)*p/d
              ! If dual energy, adjust entropy
              if(r%entropy.and.r%dual_energy.GE.0)then
                 m%uold(ind,r%ientropy,igrid)=m%uold(ind,5,igrid)/m%uold(ind,1,igrid)**r%gamma
              endif
           endif
        endif
     end do
  end do
  s%tailp(r%nlevelmax)=s%tailp(r%nlevelmax)+nstar_loc

  !---------------------------------------------------------
  ! Convert primitive variables to conservative variables
  !---------------------------------------------------------
  call cons_from_prim(r,g,m,ilevel)

  !----------------------------------------------
  ! Get the new global seed to the final state
  !----------------------------------------------
  call RngStream_SetPackageSeed(r%seed)
  gg = RngStream_CreateStream('ramses_stream')
  ! Advance to the state corresponding to the final global state
  if(nsite_cum(g%ncpu)>0)then
     call RngStream_AdvanceState(gg,0_8,2*nsite_cum(g%ncpu))
  endif
  ! Get the new seed and store it
  call RngStream_GetState(gg,r%seed)

  !---------------------------------------------------------
  ! Compute number of new stars across all CPUs.
  !---------------------------------------------------------
  nstar_cpu=0
  nstar_cpu(g%myid)=nstar_loc
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(nstar_cpu,nstar_cpu_tot,g%ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  nstar_cpu=nstar_cpu_tot
#endif
  nstar_cum=0
  do icpu=1,g%ncpu
     nstar_cum(icpu)=nstar_cum(icpu-1)+int(nstar_cpu(icpu),kind=8)
  end do

  !--------------------------------------
  ! Compute new star particle index
  !--------------------------------------
  do i=1,nstar_loc
     s%idp(s%npart-nstar_loc+i)=s%npart_tot+nstar_cum(g%myid-1)+i
  end do
  s%npart_tot=s%npart_tot+nstar_cum(g%ncpu)

!  if(g%myid==1)write(*,*)'Formed',int(nstar_cum(g%ncpu),kind=4),' new stars for a total of',int(s%npart_tot,kind=4)

#endif
#endif

end subroutine star_formation
!###########################################################
!###########################################################
!###########################################################
!###########################################################
function erfc(x)
  ! complementary error function
  implicit none
  real(kind=8) :: erfc
  real(kind=8) :: x, y
  real(kind=8) :: pv, ph
  real(kind=8) :: q0, q1, q2, q3, q4, q5, q6, q7
  real(kind=8) :: p0, p1, p2, p3, p4, p5, p6, p7
  parameter(pv= 1.26974899965115684d+01, ph= 6.10399733098688199d+00)
  parameter(p0= 2.96316885199227378d-01, p1= 1.81581125134637070d-01)
  parameter(p2= 6.81866451424939493d-02, p3= 1.56907543161966709d-02)
  parameter(p4= 2.21290116681517573d-03, p5= 1.91395813098742864d-04)
  parameter(p6= 9.71013284010551623d-06, p7= 1.66642447174307753d-07)
  parameter(q0= 6.12158644495538758d-02, q1= 5.50942780056002085d-01)
  parameter(q2= 1.53039662058770397d+00, q3= 2.99957952311300634d+00)
  parameter(q4= 4.95867777128246701d+00, q5= 7.41471251099335407d+00)
  parameter(q6= 1.04765104356545238d+01, q7= 1.48455557345597957d+01)
  y = x*x
  y = exp(-y)*x*(p7/(y+q7)+p6/(y+q6) + p5/(y+q5)+p4/(y+q4)+p3/(y+q3) &
       &       + p2/(y+q2)+p1/(y+q1)+p0/(y+q0))
  if (x < ph) y = y+2d0/(exp(pv*x)+1.0)
  erfc = y
  return
end function erfc
!###########################################################
!###########################################################
!###########################################################
!###########################################################
end module star_formation_module
