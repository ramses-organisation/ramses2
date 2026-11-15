MODULE cr_flux_module

  use amr_commons, only: dp, global_t, run_t
  use amr_parameters, only: ndim
  use hydro_parameters, only: nvar
  use cr_parameters, only: ncrgrp, smallecr
  use cr_commons

  implicit none

  private   ! default
  integer,parameter::ifcr1=0                                                      ! 0
  integer,parameter::jfcr1=1-ndim/2                                               ! 0 or 1
  integer,parameter::kfcr1=1-ndim/3                                               ! 0 or 1

#ifdef DO_CR
  public cr_unsplit, invrotatevec, rotatevec

CONTAINS
!************************************************************************
SUBROUTINE cmp_cr_flux_tensors(r, kcr, iGrp, cr_c)

  ! Compute central fluxes for a CR group, for each cell in a vector
  ! of grids.
  ! The flux tensor is a three by four tensor (2*3 and 1*2 in 1D and 2D,
  ! respectively) where the first row is CR flux (x,y,z) and
  ! the other three rows compose the Eddington tensor (see Yiang&Peng 2017)
  !------------------------------------------------------------------------
  type(run_t) :: r
  type(cr_kernel_t)::kcr
  integer::iGrp
  real(kind=8),dimension(1:ndim)::Fcr
  real(kind=8)::Ecr, cr_c, nedge
  integer::i, j, k, idim, jdim
  real(kind=8)::mu1,mu2,chi,b_norm2,Fcr_norm
  real(kind=8),dimension(1:ndim)::bcell
  real(kind=8)::cr_c2,Ecr2,aniso_term
  !------------------------------------------------------------------------
  associate(iu1=>kcr%iu1, ju1=>kcr%ju1, ku1=>kcr%ku1 &
           ,iu2=>kcr%iu2, ju2=>kcr%ju2, ku2=>kcr%ku2 &
           ,iecr=>r%iecr, inener=>r%inener)

  cr_c2=cr_c**2

  ! Loop 6X6X6 cells in grid, from -1 to 4.
  do k = ku1, ku2
  do j = ju1, ju2
  do i = iu1, iu2
     nedge = 0              ! Check if we're at a corner and if so, cycle
     if(i.lt.1 .or. i.gt.2) nedge=nedge+1
     if(ndim.gt.1 .and. (j.lt.1 .or. j.gt.2)) nedge=nedge+1
     if(ndim.gt.2 .and. (k.lt.1 .or. k.gt.2)) nedge=nedge+1
     if(nedge.ge.2) cycle

     Ecr = kcr%uloc(i,j,k,iEcr+iGrp-1)                          ! CR density in cell
     Fcr = kcr%cruloc(i,j,k,1+(iGrp-1)*ndim:ndim+(iGrp-1)*ndim) ! CR flux vector
     if(Ecr .lt. 0d0) then
        write(*,*)'negative CR density in cmp_flux_tensors. -EXITING-'
        stop
     endif
     kcr%cflx(i,j,k,1,1:ndim)= Fcr  !   First row is CR flux
     ! Rest is Eddington tensor
     if(r%cr_isotropic_pressure) then
        kcr%cflx(i,j,k,2:ndim+1,1:ndim) = 0d0
        do idim = 1, ndim
           kcr%cflx(i,j,k,idim+1,idim) = Ecr*cr_c2*(r%gamma_rad(iEcr-inener+igrp)-1.0)
        enddo
     else
        ! M1 closure
        Ecr2=Ecr**2
          ! Magnetic field, needed for M1
        do idim=1,ndim
          bcell(idim) = 0.5*(kcr%bloc(i,j,k,idim)+kcr%bloc(i,j,k,3+idim))
        end do
        b_norm2    =SUM(bcell(1:ndim)**2)
        Fcr_norm=SUM(Fcr**2)
        mu1=MIN(Fcr_norm/(cr_c2*Ecr2),1.0d0)
        mu2=(3d0+4d0*mu1)/(5d0+2d0*sqrt(4d0-3d0*mu1))
        chi=0.5d0*(1d0-mu2)
        aniso_term=(1d0-3.0d0*chi)/b_norm2
        do idim = 1, ndim
           do jdim = 1, ndim
              kcr%cflx(i,j,k,idim+1,jdim) = bcell(idim)*bcell(jdim)
           end do
        end do
        kcr%cflx(i,j,k,2:ndim+1,1:ndim) = kcr%cflx(i,j,k,2:ndim+1,1:ndim)*aniso_term
        do idim = 1, ndim
           kcr%cflx(i,j,k,idim+1,idim) = kcr%cflx(i,j,k,idim+1,idim) + chi
        end do
        kcr%cflx(i,j,k,2:ndim+1,1:ndim) = kcr%cflx(i,j,k,2:ndim+1,1:ndim)*Ecr*cr_c2
     endif
  enddo
  enddo
  enddo
  end associate

end subroutine cmp_cr_flux_tensors

!************************************************************************
SUBROUTINE cmp_cr_wavespeeds(r, kcr, iGrp, cr_c, dx, dt)

  !  Compute CR wavespeeds for given vector of sub-grids.
  !  Updates lmax in the cr kernel
  !
  !  inputs/outputs
  !  r         => global vars
  !  kcr      <=> cr kernel
  !  iGrp      => CR group number
  !  cr_c      => reduced light speed in code units
  !  dx        => cell width at current level, code units
  !  dt        => current CR timestep length, code units
  !------------------------------------------------------------------------
  type(run_t) :: r
  type(cr_kernel_t)::kcr
  integer,intent(in)::iGrp
  real(kind=8),intent(in)::cr_c, dt
  real(kind=8)::dx, Ecr, va, gmone_div_twodx
  integer::i, j, k, idim, nedge, iEgrp
  real(kind=8),dimension(1:3)::bcell, gradpcr, Dcr_vec
  real(kind=8)::norm,bdotgradp,cosp,sinp,cost,sint,bxby,Dcr_dir
  !------------------------------------------------------------------------
  associate(if2=>kcr%if2, jf2=>kcr%jf2, kf2=>kcr%kf2)
  gmone_div_twodx=(r%gamma_rad(r%iEcr-r%inener+igrp)-1.0)/(2d0*dx)
  iEgrp = r%iEcr+(iGrp-1)

  ! Loop (N+2)X(N+2)X(N+2) cells in grid, where N=2*nsubgrid = 2 by
  ! default. All dimension indices go from 0 to N+1.
  ! We only need to calculate wavespeeds for those cells which have faces to
  ! the NXNXN center cells, so by skipping the 'corners' we are reduced
  ! to fewer cells to calculate (by half for the default N=2).
  do k = kfcr1, kf2
  do j = jfcr1, jf2
  do i = ifcr1, if2

     nedge=0                 ! Check if we're at a corner and if so, cycle
     if(mod(i,if2).eq.0) nedge=nedge+1
     if(ndim.gt.1 .and. mod(j,jf2).eq.0) nedge=nedge+1
     if(ndim.gt.2 .and. mod(k,kf2).eq.0) nedge=nedge+1
     if(nedge.ge.2) cycle

     Ecr  = kcr%uloc(i, j, k, iEgrp)

     ! Magnetic field, needed to rotate Dcr and for bdotgradpcr
     norm=0.
     do idim=1,3
        bcell(idim) = 0.5*(kcr%bloc(i,j,k,idim)+kcr%bloc(i,j,k,3+idim))
        norm = norm + bcell(idim)**2
     end do
     norm = max(sqrt(norm), 1d-30)

     bxby = sqrt(bcell(1)**2+bcell(2)**2)
     if(norm.gt.1e-10) then
        sint = bxby/norm
        cost = bcell(3)/norm
     else
        sint = 1d0
        cost = 0d0
     endif
     if(bxby.gt.1e-10) then
        sinp = bcell(2)/bxby
        cosp = bcell(1)/bxby
     else
        sinp = 0d0
        cosp = 1d0
     endif
     bcell = bcell/norm

     va=0.
     if(r%cr_streaming_diffusion) va=norm/sqrt(kcr%uloc(i,j,k,1))
     if(r%cr_streaming_diffusion .and. r%cr_v_alfven.gt.0.0) va = r%cr_v_alfven

     ! Calculate grad Pcr
     gradpcr(1) = (kcr%uloc(i+1,j  ,k  ,iEgrp) - kcr%uloc(i-1,j  ,k  ,iEgrp)) * gmone_div_twodx
#if NDIM>1
     gradpcr(2) = (kcr%uloc(i  ,j+1,k  ,iEgrp) - kcr%uloc(i  ,j-1,k  ,iEgrp)) * gmone_div_twodx
#endif
#if NDIM>2
     gradpcr(3) = (kcr%uloc(i  ,j  ,k+1,iEgrp) - kcr%uloc(i  ,j  ,k-1,iEgrp)) * gmone_div_twodx
#endif

     ! Calculate B dot grad Pcr
     bdotgradp = 0.
     do idim=1,ndim
        bdotgradp = bdotgradp + bcell(idim) * gradpcr(idim)
     enddo

     ! Diffusion, eq 10 in JO17
     Dcr_vec = (/ r%cr_d_code(igrp),                           &
                  r%cr_d_code(igrp)*r%cr_d_perp_factors(iGrp), &
                  r%cr_d_code(igrp)*r%cr_d_perp_factors(iGrp) /)
     if(r%cr_streaming_diffusion) &
          Dcr_vec(1) = Dcr_vec(1) + &
          min(r%cr_dmax_code, 3./max(abs(bdotgradp),1d-50) * va * r%gamma_rad(r%iEcr-r%inener+iGrp) * max(Ecr,smallecr))

     ! Rotate Dcr_vec so it is parallel with B, hence
     ! describing Dcr in the simulation coordinate system
     ! (instead of the B coordinate system)
     call invrotatevec(sint, cost, sinp, cosp, Dcr_vec(1), Dcr_vec(2), Dcr_vec(3))

     ! Calculate wavespeeds
     Dcr_dir = abs(Dcr_vec(1)) ! x component of rotated Dcr
     kcr%lmax(i,j,k,1) = &
          cmp_cr_lmax(r, dx, Dcr_dir, r%gamma_rad(r%iEcr-r%inener+igrp), cr_c, dt)

#if NDIM>1
     Dcr_dir = abs(Dcr_vec(2)) ! y component of rotated Dcr
     kcr%lmax(i,j,k,2) = &
          cmp_cr_lmax(r, dx, Dcr_dir, r%gamma_rad(r%iEcr-r%inener+igrp), cr_c, dt)
#endif

#if NDIM>2
     Dcr_dir = abs(Dcr_vec(3)) ! z component of rotated Dcr
     kcr%lmax(i,j,k,3) = &
          cmp_cr_lmax(r, dx, Dcr_dir, r%gamma_rad(r%iEcr-r%inener+igrp), cr_c, dt)
#endif

  end do
  end do
  end do
  end associate

END SUBROUTINE cmp_cr_wavespeeds

!************************************************************************
FUNCTION cmp_cr_lmax(r, dx, dcoeff, gamma_cr, cr_c, dt)

! Compute maximum local wavespeed
!------------------------------------------------------------------------
  type(run_t) :: r
  real(kind=8)::dx, cmp_cr_lmax, gamma_cr, cr_c, dt
  real(kind=8)::tau, dcoeff, r_factor
!------------------------------------------------------------------------
  tau = 0.5d0 * dx**2 / dcoeff / dt
  if(tau.lt.1e-3) then
    r_factor = sqrt((1.0 - 0.5*tau**2))
  else
    r_factor = sqrt((1.-exp(-min(tau,10.)**2))/min(tau,1e8)**2) ! Capital R on p 6 in YP17
  endif
  if(r%cr_isotropic_pressure) then
     cmp_cr_lmax = r_factor * cr_c * sqrt(gamma_cr-1.0)
  else
     cmp_cr_lmax = r_factor * cr_c
  endif

END FUNCTION cmp_cr_lmax

!************************************************************************
FUNCTION cmp_cr_face(fdn, fup, udn, uup, lminus, lplus)

! Compute HLLE intercell fluxes for all (four) CR variables.
! fdn    => flux function in the cell downwards from the border
! fup    => flux function in the cell upwards from the border
! udn    => state (CR density and flux downwards from the border)
! uup    => state (CR density and flux upwards from the border)
! lminus => minimum intercell wavespeed
! lplus  => maximum intercell wavespeed
! returns      flux vector for the given state variables, i.e. line nr dim
!              in the 3*4 flux function tensor
!------------------------------------------------------------------------
  real(kind=8),dimension(nDim+1)::fdn, fup, udn, uup, cmp_cr_face
  real(kind=8)::lminus, lplus, llmax
!------------------------------------------------------------------------
  llmax = max(abs(lplus), abs(lminus))
  cmp_cr_face = 0.5D0 * (fdn + fup - llmax * (uup - udn)) !LLF flux
END FUNCTION cmp_cr_face

!************************************************************************
SUBROUTINE cr_unsplit(r,kcr,cr_c,dx,dt)
!  Compute intercell fluxes for one CR group in all dimensions,
!  using the Eddington tensor with the Yiang+Peng'17 closure relation.
!  The intercell fluxes are the right-hand sides of the equation:
!      dq/dt = - nabla dot f   (eq 12 in YP17)
!  where q=[Ecr, Fx/ccr^2, Fy/ccr^2, Fz/ccr^2], ccr the reduced wavespeed
!  and f the Eddington pressure tensor. A flux at index i,j,k represents
!  flux across the lower faces of that cell, i.e. at i-1/2 etc.
!
!  inputs/outputs
!  r           => access to 'global' variables
!  kcr         => CR kernel
!  dx          => cell width
!  dt          => time step
!  iGrp        => CR group number
!
!  vars in kcr:
!  iu1,iu2     |First and last index of input array,
!  ju1,ju2     |cell centered,
!  ku1,ku2     |including buffer cells (6x6x6).
!  if1,if2     |First and last index of output array,
!  jf1,jf2     |edge centered, for active
!  kf1,kf2     |cells only (3x3x3).
!------------------------------------------------------------------------
  type(run_t) :: r
  type(cr_kernel_t):: kcr
  real(kind=8):: cr_c, dx, dt
  ! Upwards and downwards fluxes and states of the group
  real(kind=8),dimension(nDim+1),save:: f1,f2,f3,f4,u1,u2,u3,u4
  real(kind=8):: a2,a3
  real(kind=8):: lminus, lplus                     ! Intercell wavespeeds
  real(kind=8):: meandiffv, dtdx, prod(ndim+1)
  integer ::i, j, k, iGrp
  real(kind=8),dimension(ndim+1),save:: slopeLM,slopeRM,slopeM
  real(kind=8),dimension(ndim+1),save:: slopeLL,slopeL
!------------------------------------------------------------------------
  associate(if1=>kcr%if1, if2=>kcr%if2, jf1=>kcr%jf1                    &
           ,jf2=>kcr%jf2, kf1=>kcr%kf1, kf2=>kcr%kf2, crin=>kcr%cruloc  &
           ,cFlx=>kcr%cFlx, uin=>kcr%uloc, lmax=>kcr%lmax, iEcr=>r%iEcr)

  do iGrp = 1, ncrgrp

  ! Compute flux tensors for all the cells with correction
  call cmp_cr_flux_tensors(r, kcr, iGrp, cr_c)
  ! Wavespeeds in each cell (kcr%lmax)
  call cmp_cr_wavespeeds(r, kcr, iGrp, cr_c, dx, dt)

  ! Solve for 1D flux in X direction
  !----------------------------------------------------------------------
  dtdx=dt/dx
  do i=if1,if2                                 !
  do j=jf1,jf2                                 !  each cell in grid (1-3)
  do k=kf1,kf2                                 !
     if(ndim.gt.1 .and. j.eq.jf2) cycle
     if(ndim.gt.2 .and. k.eq.kf2) cycle
     f1 = cFlx(i-2, j, k, :, 1)               !
     f2 = cFlx(i-1, j, k, :, 1)               ! upwards and downwards
     f3 = cFlx(i,   j, k, :, 1)               ! conditions
     f4 = cFlx(i+1, j, k, :, 1)               !
     u1(1) = uin(i-2, j, k, iEcr+iGrp-1)
     u2(1) = uin(i-1, j, k, iEcr+iGrp-1)
     u3(1) = uin(i,   j, k, iEcr+iGrp-1)
     u4(1) = uin(i+1, j, k, iEcr+iGrp-1)
     u1(2:ndim+1) = crin(i-2, j, k, 1+(iGrp-1)*ndim:iGrp*ndim)
     u2(2:ndim+1) = crin(i-1, j, k, 1+(iGrp-1)*ndim:iGrp*ndim)
     u3(2:ndim+1) = crin(i,   j, k, 1+(iGrp-1)*ndim:iGrp*ndim)
     u4(2:ndim+1) = crin(i+1, j, k, 1+(iGrp-1)*ndim:iGrp*ndim)

     ! Second-order interpolation with using Van-Leer slope Limiter
     ! interpolation of U
     slopeLM = (f3-f2)/dx
     slopeRM = (f4-f3)/dx
     prod = slopeLM*slopeRM
     slopeM=0.
     where(prod.gt.0.) slopeM=2.*prod/(slopeLM+slopeRM)
     slopeLL = (f2-f1)/dx
     prod = slopeLL*slopeLM
     slopeL=0.
     where(prod.gt.0) slopeL=2.*prod/(slopeLL+slopeLM)
     f2 = f2+slopeL*0.5d0*dx
     f3 = f3-slopeM*0.5d0*dx

     ! interpolation of F
     slopeLM = (u3-u2)/dx
     slopeRM = (u4-u3)/dx
     prod = slopeLM*slopeRM
     slopeM=0.
     where(prod.gt.0) slopeM=2.*prod/(slopeLM+slopeRM)
     slopeLL = (u2-u1)/dx
     prod = slopeLL*slopeLM
     slopeL=0.
     where(prod.gt.0.) slopeL=2.*prod/(slopeLL+slopeLM)
     u2 = u2+slopeL*0.5d0*dx
     u3 = u3-slopeM*0.5d0*dx

     meandiffv = 0.5*( lmax(i-1,j,k,1) + lmax(i,j,k,1) )
     a2 = min(-meandiffv, -lmax(i-1,j,k,1))
     a2 = max(a2, -cr_c*sqrt(r%gamma_rad(iEcr-r%inener+igrp)-1.0))
     a3 = max(meandiffv, lmax(i,j,k,1))
     a3 = min(a3, cr_c*sqrt(r%gamma_rad(iEcr-r%inener+igrp)-1.0))

     lminus = min(a2,0.)
     lplus = max(a3,0.)

     kcr%crflux( i, j, k, 1+(iGrp-1)*(ndim+1):iGrp*(ndim+1), 1)=&
          cmp_cr_face( f2, f3, u2, u3, lminus, lplus)*dtdx

  end do
  end do
  end do

  ! Solve for 1D flux in Y direction
  !----------------------------------------------------------------------
#if NDIM>1
  dtdx=dt/dx
  do i=if1,if2-1
  do j=jf1,jf2
  do k=kf1,kf2
     if(ndim.gt.2 .and. k.eq.kf2) cycle
     f1 = cFlx(i, j-2, k, :, 2)
     f2 = cFlx(i, j-1, k, :, 2)
     f3 = cFlx(i, j,   k, :, 2)
     f4 = cFlx(i, j+1, k, :, 2)
     u1(1) = uin(i, j-2, k, iEcr+iGrp-1)
     u2(1) = uin(i, j-1, k, iEcr+iGrp-1)
     u3(1) = uin(i, j,   k, iEcr+iGrp-1)
     u4(1) = uin(i, j+1, k, iEcr+iGrp-1)
     u1(2:ndim+1) = crin(i, j-2, k, 1+(iGrp-1)*ndim:iGrp*ndim)
     u2(2:ndim+1) = crin(i, j-1, k, 1+(iGrp-1)*ndim:iGrp*ndim)
     u3(2:ndim+1) = crin(i, j,   k, 1+(iGrp-1)*ndim:iGrp*ndim)
     u4(2:ndim+1) = crin(i, j+1, k, 1+(iGrp-1)*ndim:iGrp*ndim)

     slopeLM = (f3-f2)/dx
     slopeRM = (f4-f3)/dx
     prod = slopeLM*slopeRM
     slopeM=0.
     where(prod.gt.0.) slopeM=2.*prod/(slopeLM+slopeRM)
     slopeLL = (f2-f1)/dx
     prod = slopeLL*slopeLM
     slopeL=0.
     where(prod.gt.0) slopeL=2.*prod/(slopeLL+slopeLM)
     f2 = f2+slopeL*0.5d0*dx
     f3 = f3-slopeM*0.5d0*dx

     slopeLM = (u3-u2)/dx
     slopeRM = (u4-u3)/dx
     prod = slopeLM*slopeRM
     slopeM=0.
     where(prod.gt.0) slopeM=2.*prod/(slopeLM+slopeRM)
     slopeLL = (u2-u1)/dx
     prod = slopeLL*slopeLM
     slopeL=0.
     where(prod.gt.0.) slopeL=2.*prod/(slopeLL+slopeLM)
     u2 = u2+slopeL*0.5d0*dx
     u3 = u3-slopeM*0.5d0*dx

     meandiffv = 0.5*( lmax(i,j-1,k,2) + lmax(i,j,k,2) )
     a2 = min(-meandiffv, -lmax(i,j-1,k,2))
     a2 = max(a2, -cr_c*sqrt(r%gamma_rad(iEcr-r%inener+igrp)-1.0))
     a3 = max(meandiffv, lmax(i,j,k,2))
     a3 = min(a3, cr_c*sqrt(r%gamma_rad(iEcr-r%inener+igrp)-1.0))
     lminus = min(a2,0.)
     lplus = max(a3,0.)

     kcr%crflux(i, j, k, 1+(iGrp-1)*(ndim+1):iGrp*(ndim+1), 2)=&
          cmp_cr_face( f2, f3, u2, u3, lminus, lplus)*dtdx

  end do
  end do
  end do
#endif

  ! Solve for 1D flux in Z direction
  !----------------------------------------------------------------------
#if NDIM>2
  dtdx=dt/dx
  do i=if1,if2-1
  do j=jf1,jf2-1
  do k=kf1,kf2
     f1 = cFlx(i, j, k-2, :, 3)
     f2 = cFlx(i, j, k-1, :, 3)
     f3 = cFlx(i, j, k,   :, 3)
     f4 = cFlx(i, j, k+1, :, 3)
     u1(1) = uin(i, j, k-2, iEcr+iGrp-1)
     u2(1) = uin(i, j, k-1, iEcr+iGrp-1)
     u3(1) = uin(i, j, k,   iEcr+iGrp-1)
     u4(1) = uin(i, j, k+1, iEcr+iGrp-1)
     u1(2:ndim+1) = crin(i, j, k-2, 1+(iGrp-1)*ndim:iGrp*ndim)
     u2(2:ndim+1) = crin(i, j, k-1, 1+(iGrp-1)*ndim:iGrp*ndim)
     u3(2:ndim+1) = crin(i, j, k,   1+(iGrp-1)*ndim:iGrp*ndim)
     u4(2:ndim+1) = crin(i, j, k+1, 1+(iGrp-1)*ndim:iGrp*ndim)

     slopeLM = (f3-f2)/dx
     slopeRM = (f4-f3)/dx
     prod = slopeLM*slopeRM
     slopeM=0.
     where(prod.gt.0.) slopeM=2.*prod/(slopeLM+slopeRM)
     slopeLL = (f2-f1)/dx
     prod = slopeLL*slopeLM
     slopeL=0.
     where(prod.gt.0) slopeL=2.*prod/(slopeLL+slopeLM)
     f2 = f2+slopeL*0.5d0*dx
     f3 = f3-slopeM*0.5d0*dx

     slopeLM = (u3-u2)/dx
     slopeRM = (u4-u3)/dx
     prod = slopeLM*slopeRM
     slopeM=0.
     where(prod.gt.0) slopeM=2.*prod/(slopeLM+slopeRM)
     slopeLL = (u2-u1)/dx
     prod = slopeLL*slopeLM
     slopeL=0.
     where(prod.gt.0.) slopeL=2.*prod/(slopeLL+slopeLM)
     u2 = u2+slopeL*0.5d0*dx
     u3 = u3-slopeM*0.5d0*dx

     meandiffv = 0.5*( lmax(i,j,k-1,3) + lmax(i,j,k,3) )
     a2 = min(-meandiffv, -lmax(i,j,k-1,3))
     a2 = max(a2, -cr_c*sqrt(r%gamma_rad(iEcr-r%inener+igrp)-1.0))
     a3 = max(meandiffv, lmax(i,j,k,3))
     a3 = min(a3, cr_c*sqrt(r%gamma_rad(iEcr-r%inener+igrp)-1.0))
     lminus = min(a2,0.)
     lplus  = max(a2,0.)

     kcr%crflux( i, j, k, 1+(iGrp-1)*(ndim+1):iGrp*(ndim+1), 3)=&
          cmp_cr_face( f2, f3, u2, u3, lminus, lplus)*dtdx

  end do
  end do
  end do
#endif

  end do ! end loop over CR groups

  END ASSOCIATE

END SUBROUTINE cr_unsplit

!************************************************************************
SUBROUTINE rotatevec(sint, cost, sinp, cosp, v1, v2, v3)
  !  Rotate vector v by t=theta and p=phi
  !  i.e. rotate to the local coordinate system from theta, phi.
  !  Hence the x-component of the result is the component of v parallel
  !  with the theta,phi vector.
  !------------------------------------------------------------------------
    implicit none
    real(dp),intent(in):: sint, cost, sinp, cosp
    real(dp),intent(inout)::v1,v2,v3
    real(dp)::newv1, newv3
  !------------------------------------------------------------------------
    ! First apply R1, then apply R2
    newv1 =  cosp * v1 + sinp * v2
    v2 = -sinp * v1 + cosp * v2

    ! Now apply R2
    v1 =  sint * newv1 + cost * v3
    newv3 = -cost * newv1 + sint * v3
    v3 = newv3
END SUBROUTINE rotatevec

!************************************************************************
SUBROUTINE invrotatevec(sint, cost, sinp, cosp, v1, v2, v3)
  !  Inverse-rotate vector v by t=theta and p=phi,
  !  i.e. rotate v onto theta, pi
  !
  !------------------------------------------------------------------------
    implicit none
    real(dp),intent(in):: sint, cost, sinp, cosp
    real(dp),intent(inout)::v1,v2,v3
    real(dp)::newv1, newv2
  !------------------------------------------------------------------------
    ! First apply R2^-1, then apply R1^-1
    newv1 = sint * v1 - cost * v3
    v3 = cost * v1 + sint * v3

    ! now apply R1^-1
    v1 = cosp * newv1 - sinp * v2
    newv2 = sinp * newv1 + cosp * v2
    v2 = newv2
END SUBROUTINE invrotatevec

#endif

END MODULE cr_flux_module
