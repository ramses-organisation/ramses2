module courant_fine_module
#ifdef _CUDA
  use gpu_runner, only: gpu_cmpdt, gpu_cmpdt_2
#endif
#ifdef _METAL
  use metal_runner, only: metal_cmpdt
#endif

  type :: out_courant_fine_t
     real(kind=8)::mass,ekin,eint,emag,dt
  end type out_courant_fine_t

contains
!###########################################################
!###########################################################
!###########################################################
!###########################################################
recursive subroutine r_courant_fine(pst,ilevel,input_size,output,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  type(out_courant_fine_t)::output, next_output

  integer::ilevel
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_COURANT_FINE,pst%iUpper+1,input_size,output_size,ilevel)
     call r_courant_fine(pst%pLower,ilevel,input_size,output,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size,next_output)
     output%mass=output%mass+next_output%mass
     output%ekin=output%ekin+next_output%ekin
     output%eint=output%eint+next_output%eint
     output%emag=output%emag+next_output%emag
     output%dt=MIN(output%dt,next_output%dt)
  else
#ifdef _CUDA
     call gpu_cmpdt_2(pst%s,ilevel,output%mass,output%ekin,output%eint,output%emag,output%dt)
!     call gpu_cmpdt(pst%s,ilevel,output%mass,output%ekin,output%eint,output%emag,output%dt)
#elif defined(_METAL)
     call metal_cmpdt(pst%s,ilevel,output%mass,output%ekin,output%eint,output%emag,output%dt)
#else
     call courant_fine(pst%s%r,pst%s%g,pst%s%m,ilevel,output%mass,output%ekin,output%eint,output%emag,output%dt)
#endif
  endif

end subroutine r_courant_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine courant_fine(r,g,m,ilevel,mass,ekin,eint,emag,dt)
  use amr_parameters, only: nvector, ndim, twotondim
  use hydro_parameters, only: nvar, nener
  use amr_commons, only: run_t,global_t,mesh_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer::ilevel
  real(kind=8)::mass,ekin,eint,emag,dt
  !----------------------------------------------------------------------
  ! Using the Courant-Friedrich-Levy stability condition,               !
  ! this routine computes the maximum allowed time-step.                !
  !----------------------------------------------------------------------
  integer::ivar,idim,ind,igrid
  real(kind=8)::dt_lev,dx,vol
  real(kind=8),dimension(1:nvar)::uu
  real(kind=8),dimension(1:ndim)::gg
  real(kind=8),dimension(1:6)::bb

#ifdef HYDRO

  ! Mesh spacing at that level
  dx=r%boxlen/2**ilevel
  vol=dx**ndim

  mass=0.0D0
  ekin=0.0D0
  eint=0.0D0
  emag=0.0D0
  dt=dx/r%smallc

  if(r%induction)call reset_init(r,g,m,ilevel)

  ! Loop over active grids by vector sweeps
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(igrid, ind, ivar, idim, uu, bb, gg, dt_lev) REDUCTION(+:mass, ekin, eint, emag) REDUCTION(MIN:dt)
  do igrid=m%head(ilevel),m%tail(ilevel)

     ! Loop over cells
     do ind=1,twotondim

        ! Gather leaf cells
        if(.NOT. m%grid(igrid)%refined(ind))then

           ! Gather hydro variables
           do ivar=1,nvar
              uu(ivar)=m%uold(ind,ivar,igrid)
           end do

           ! Gather MHD variables
           bb=0.0d0
#ifdef MHD
           do ivar=1,6
              bb(ivar)=m%bold(ind,ivar,igrid)
           end do
#endif
           ! Gather gravitational acceleration
#ifdef GRAV
           do idim=1,ndim
              gg(idim)=m%f(ind,idim,igrid)
           end do
#else
           do idim=1,ndim
              gg(idim)=r%constant_gravity(idim)
           end do
#endif
           ! Gather turbulent driving
#ifdef TURB
           if(r%turb)then
              do idim=1,ndim
                 gg(idim)=gg(idim)+m%fturb(ind,idim,igrid)
              end do
           endif
#endif
           ! Compute total mass
           mass=mass+uu(1)*vol

           ! Compute total energy
           ekin=ekin+uu(5)*vol

           ! Compute total internal energy
           eint=eint+uu(5)*vol
           do idim=1,3
              eint=eint-0.5D0*uu(1+idim)**2/max(uu(1),r%smallr)*vol
           end do
#if NENER>0
           do ivar=1,nener
              eint=eint-uu(5+ivar)*vol
           end do
#endif
#ifdef MHD
           ! Compute total magnetic energy
           do idim=1,3
              emag=emag+0.125D0*(bb(idim)+bb(idim+3))**2*vol
              eint=eint-0.125D0*(bb(idim)+bb(idim+3))**2*vol
           end do
#endif
           ! Compute CFL time-step
           call cmpdt(r,uu,bb,gg,dx,dt_lev)
           dt=min(dt,dt_lev)

        endif

     end do
     ! End loop over cells
  end do
!$OMP END PARALLEL DO
  ! End loop over grids

#endif

end subroutine courant_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine reset_init(r,g,m,ilevel)
  use amr_parameters, only: nvector, ndim, twotondim
  use hydro_parameters, only: nvar, nener
  use amr_commons, only: run_t, global_t, mesh_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer::ilevel
  !----------------------------------------------------------------------
  ! This routine overwrites all Euler variables with initial conditions.
  !----------------------------------------------------------------------
  ! Local variables
  integer::igrid,ngrid,ind,idim,nstride,i,ivar,irad
  integer::l,nfine,ii,jj,kk
  real(kind=8),dimension(1:nvector,1:ndim)::xx
#ifdef MHD
  real(kind=8),dimension(1:nvector,1:nvar+3-ndim)::qq
#else
  real(kind=8),dimension(1:nvector,1:nvar)::qq
#endif
  real(kind=8)::dx,dxmin
  real(kind=8)::rr,vx,vy,vz,pp
  real(kind=8)::bx,by,bz
  real(kind=8)::eint,ekin,emag,erad

#ifdef HYDRO

  ! Mesh size at level ilevel in code units
  dx=r%boxlen/2**ilevel
  ! Loop over grids by vector sweeps
  do igrid=m%head(ilevel),m%tail(ilevel),nvector
     ngrid=MIN(nvector,m%tail(ilevel)-igrid+1)
     ! Loop over cells
     do ind=1,twotondim
        ! Compute cell centre position in code units
        do idim=1,ndim
           nstride=2**(idim-1)
           do i=1,ngrid
              xx(i,idim)=(2*m%grid(igrid+i-1)%ckey(idim)+MOD((ind-1)/nstride,2)+0.5)*dx-m%skip(idim)
           end do
        end do
        ! Call initial condition routine
        call condinit(r,g,xx,qq,dx,ngrid)
        ! Comvert to conservative variables
        do i=1,ngrid
           rr=qq(i,1)
           vx=qq(i,2)
           vy=qq(i,3)
           vz=qq(i,4)
           pp=qq(i,5)
           ekin=0.5d0*rr*(vx**2+vy**2+vz**2)
           eint=pp/(r%gamma-1.0)
           emag=0.0d0
#ifdef MHD
           ! Compute magnetic energy for all cells
           bx=0.5d0*(m%bold(ind,1,igrid+i-1)+m%bold(ind,4,igrid+i-1))
           by=0.5d0*(m%bold(ind,2,igrid+i-1)+m%bold(ind,5,igrid+i-1))
           bz=0.5d0*(m%bold(ind,3,igrid+i-1)+m%bold(ind,6,igrid+i-1))
           emag=0.5d0*(bx**2+by**2+bz**2)
#endif
           erad=0.0d0
#if NENER>0
           ! Compute non-thermal energy densities
           do irad=1,nener
              m%uold(ind,5+irad,igrid+i-1)=qq(i,5+irad)/(r%gamma_rad(irad)-1.0d0)
              erad=erad+m%uold(ind,5+irad,igrid+i-1)
           end do
#endif
           ! Compute total fluid energy density
           m%uold(ind,5,igrid+i-1)=eint+ekin+erad+emag
           ! Compute momentum density
           do idim=1,3
              m%uold(ind,idim+1,igrid+i-1)=rr*qq(i,idim+1)
           end do
           ! Compute mass density
           m%uold(ind,1,igrid+i-1)=rr
#if NVAR>5+NENER
           ! Compute passive scalar density
           do ivar=6+nener,nvar
              m%uold(ind,ivar,igrid+i-1)=rr*qq(i,ivar)
           enddo
#endif
        end do
     end do
  end do
#endif

end subroutine reset_init
!###########################################################
!###########################################################
!###########################################################
!###########################################################
end module courant_fine_module
