!###########################################################
!###########################################################
!###########################################################
!###########################################################
#ifdef MHD
subroutine hydro_refine(r,ug,um,ud,bg,bm,bd,ok)
#else
subroutine hydro_refine(r,ug,um,ud,ok)
#endif
  use amr_parameters, only: ndim
  use amr_commons, only: run_t
  use hydro_parameters, only: nvar,nener
  use const
  implicit none
  ! routine arguments
  type(run_t)::r
  real(kind=8)::ug(1:nvar)
  real(kind=8)::um(1:nvar)
  real(kind=8)::ud(1:nvar)
#ifdef MHD
  real(kind=8)::bg(1:6)
  real(kind=8)::bm(1:6)
  real(kind=8)::bd(1:6)
#endif
  logical ::ok
  ! local variables
  integer::idim,ii
#if NENER>0
  integer::irad
#endif
#if NVAR>5+NENER
  integer::n
#endif
  real(kind=8)::eking,ekinm,ekind,eradg,eradm,eradd
  real(kind=8)::emagg,emagm,emagd
  real(kind=8)::dg,dm,dd,pg,pm,pd,vg,vm,vd,cg,cm,cd,error

  ! Density
  ug(1) = max(ug(1),r%smallr)
  um(1) = max(um(1),r%smallr)
  ud(1) = max(ud(1),r%smallr)

  ! Velocity
  do idim = 1,3
     ug(idim+1) = ug(idim+1)/ug(1)
     um(idim+1) = um(idim+1)/um(1)
     ud(idim+1) = ud(idim+1)/ud(1)
  end do

  ! Cell-centered magnetic field and energy
  emagg = zero
  emagm = zero
  emagd = zero
#ifdef MHD
  do idim = 1,3
     bg(idim) = half*(bg(idim)+bg(idim+3))
     bm(idim) = half*(bm(idim)+bm(idim+3))
     bd(idim) = half*(bd(idim)+bd(idim+3))
     emagg = emagg + half*bg(idim)**2
     emagm = emagm + half*bm(idim)**2
     emagd = emagd + half*bd(idim)**2
  end do
#endif

  ! Kinetic energy
  eking = zero
  ekinm = zero
  ekind = zero
  do idim = 1,3
     eking = eking + half*ug(1)*ug(idim+1)**2
     ekinm = ekinm + half*um(1)*um(idim+1)**2
     ekind = ekind + half*ud(1)*ud(idim+1)**2
  end do

  ! Non-thermal energies
  eradg = zero
  eradm = zero
  eradd = zero
#if NENER>0
  do irad = 1,nener
     eradg = eradg + ug(5+irad)
     eradm = eradm + um(5+irad)
     eradd = eradd + ud(5+irad)
  end do
#endif

  ! Thermal pressure
  pg = (r%gamma-one)*(ug(5)-eking-emagg-eradg)
  pm = (r%gamma-one)*(um(5)-ekinm-emagm-eradm)
  pd = (r%gamma-one)*(ud(5)-ekind-emagd-eradd)

  ! Passive scalars
#if NVAR>5+NENER
  do n = 6+nener,nvar
     ug(n) = ug(n)/ug(1)
     um(n) = um(n)/um(1)
     ud(n) = ud(n)/ud(1)
  end do
#endif

  ! Compute errors
  if(r%err_grad_d >= 0.)then
     dg=ug(1); dm=um(1); dd=ud(1)
     error=2.0d0*MAX( &
          & ABS((dd-dm)/(dd+dm+r%floor_d)) , &
          & ABS((dm-dg)/(dm+dg+r%floor_d)) )
     ok = ok .or. error > r%err_grad_d
  end if

  if(r%err_grad_p >= 0.)then
     error=2.0d0*MAX( &
          & ABS((pd-pm)/(pd+pm+r%floor_p)), &
          & ABS((pm-pg)/(pm+pg+r%floor_p)) )
     ok = ok .or. error > r%err_grad_p
  end if

  if(r%err_grad_u >= 0.)then
     if(r%induction)then
        do idim=1,3
           vg=ug(idim+1); vm=um(idim+1); vd=ud(idim+1)
           error=2.0d0*MAX( &
                & ABS((vd-vm)/(ABS(vd)+ABS(vm)+r%floor_u)) , &
                & ABS((vm-vg)/(ABS(vm)+ABS(vg)+r%floor_u)) )
           ok = ok .or. error > r%err_grad_u
        end do
     else
        do idim=1,3
           vg=ug(idim+1); vm=um(idim+1); vd=ud(idim+1)
           cg=sqrt(max(r%gamma*pg/ug(1),r%floor_u**2))
           cm=sqrt(max(r%gamma*pm/um(1),r%floor_u**2))
           cd=sqrt(max(r%gamma*pd/ud(1),r%floor_u**2))
           error=2.0d0*MAX( &
                & ABS((vd-vm)/(cd+cm+ABS(vd)+ABS(vm)+r%floor_u)) , &
                & ABS((vm-vg)/(cm+cg+ABS(vm)+ABS(vg)+r%floor_u)) )
           ok = ok .or. error > r%err_grad_u
        end do
     end if
  end if

#if NENER>0
  do irad = 1,nener
    if(r%err_grad_prad(irad) >= 0.)then
      pg=ug(5+irad); pm=um(5+irad); pd=ud(5+irad)
      error=2.0d0*MAX( &
        & ABS((pd-pm)/(pd+pm+r%floor_prad(irad))), &
        & ABS((pm-pg)/(pm+pg+r%floor_prad(irad))) )
      ok = ok .or. error > r%err_grad_prad(irad)
    end if
  end do
#endif


  !TODO(code): fix in the case of RTZ
  if(r%neq_chem .and. r%err_grad_xHII >= 0.)then
     ii=r%iIons
     error=2.0d0*MAX( &
          & ABS((ud(ii)-um(ii))/(ud(ii)+um(ii)+r%floor_xHII)), &
          & ABS((um(ii)-ug(ii))/(um(ii)+ug(ii)+r%floor_xHII)) )
     ok = ok .or. error > r%err_grad_xHII
  end if

  if(r%neq_chem .and. r%err_grad_xHI >= 0.)then
     ii=r%iIons
     error=2.0d0*MAX( &
          & ABS((um(ii)-ud(ii))/(2.-ud(ii)-um(ii)+r%floor_xHI)), &
          & ABS((ug(ii)-um(ii))/(2.-um(ii)-ug(ii)+r%floor_xHI)) )
     ok = ok .or. error > r%err_grad_xHI
  end if

#ifdef MHD

  if(r%err_grad_b2 >= 0.)then
     pg=emagg; pm=emagm; pd=emagd
     error=2.0d0*MAX( &
          & ABS((pd-pm)/(pd+pm+r%floor_b2)), &
          & ABS((pm-pg)/(pm+pg+r%floor_b2)) )
     ok = ok .or. error > r%err_grad_b2
  end if

  if(r%err_grad_A >= 0.)then
     idim=1
     vg=bg(idim)
     vm=bm(idim)
     vd=bd(idim)
     cg=sqrt(emagg)
     cm=sqrt(emagm)
     cd=sqrt(emagd)
     error=2.0d0*MAX( &
          & ABS((vd-vm)/(cd+cm+r%floor_A)) , &
          & ABS((vm-vg)/(cm+cg+r%floor_A)) )
     ok = ok .or. error > r%err_grad_A
  end if

  if(r%err_grad_B >= 0.)then
     idim=2
     vg=bg(idim)
     vm=bm(idim)
     vd=bd(idim)
     cg=sqrt(emagg)
     cm=sqrt(emagm)
     cd=sqrt(emagd)
     error=2.0d0*MAX( &
          & ABS((vd-vm)/(cd+cm+r%floor_B)) , &
          & ABS((vm-vg)/(cm+cg+r%floor_B)) )
     ok = ok .or. error > r%err_grad_B
  end if

  if(r%err_grad_C >= 0.)then
     idim=3
     vg=bg(idim)
     vm=bm(idim)
     vd=bd(idim)
     cg=sqrt(emagg)
     cm=sqrt(emagm)
     cd=sqrt(emagd)
     error=2.0d0*MAX( &
          & ABS((vd-vm)/(cd+cm+r%floor_C)) , &
          & ABS((vm-vg)/(cm+cg+r%floor_C)) )
     ok = ok .or. error > r%err_grad_C
  end if
#endif

end subroutine hydro_refine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmpdt(r,uu,bb,gg,dx,dt)
  use amr_parameters, only: ndim
  use hydro_parameters, only: nvar, nener
  use amr_commons, only: run_t
  use const
  implicit none
  type(run_t)::r
  real(kind=8)::dx,dt
  real(kind=8),dimension(1:nvar)::uu
  real(kind=8),dimension(1:ndim)::gg
  real(kind=8),dimension(1:6)::bb

  real(kind=8)::dtcell,smallp
  real(kind=8)::b2,a2,c2,cfast2,ctot,grav
  integer::idim
#if NENER>0
  integer::irad
#endif

  smallp = r%smallc**2/r%gamma

  ! Density
  uu(1) = max(uu(1),r%smallr)

  ! Velocity
  do idim = 1,3
     uu(idim+1) = uu(idim+1)/uu(1)
  end do

  ! Magnetic field
  b2 = zero
#ifdef MHD
  do idim = 1,3
     bb(idim) = half*(bb(idim)+bb(idim+3))
     uu(5) = uu(5) - half*bb(idim)**2
     b2 = b2 + bb(idim)**2
  end do
#endif

  ! Internal energy
  do idim = 1,3
     uu(5) = uu(5)-half*uu(1)*uu(idim+1)**2
  end do
#if NENER>0
  do irad = 1,nener
     uu(5) = uu(5)-uu(5+irad)
  end do
#endif

  ! Debug
  if(r%debug)then
     if(uu(5).le.0.0D0.or.uu(1).le.r%smallr)then
        write(*,*)'stop in cmpdt'
        write(*,*)'dx   =',dx
        write(*,*)'rho  =',uu(1)
        write(*,*)'P    =',uu(5)
        write(*,*)'vel  =',uu(2:4)
        stop
     end if
  end if

  ! Compute pressures
  uu(5) = max((r%gamma-one)*uu(5),uu(1)*smallp)
#if NENER>0
  do irad = 1,nener
     uu(5+irad) = (r%gamma_rad(irad)-one)*uu(5+irad)
  end do
#endif

  ! Compute sound speed squared
  a2 = r%gamma*uu(5)
#if NENER>0
  do irad = 1,nener
     a2 = a2 + r%gamma_rad(irad)*uu(5+irad)
  end do
#endif
  a2=a2/uu(1)

  ! Compute wave speed (note that we use ndim here, not 3)
  ctot = zero
  if(r%induction)then
     do idim = 1,ndim
        ctot = ctot+abs(uu(idim+1))+r%smallc
     end do
  else
     do idim = 1,ndim
        c2 = half*(b2/uu(1)+a2)
        cfast2 = c2+sqrt(c2**2-a2*bb(idim)**2/uu(1))
        ctot = ctot+abs(uu(idim+1))+sqrt(cfast2)
     end do
  endif

  ! Compute gravity strength ratio
  grav = zero
  do idim = 1,ndim
     grav = grav+abs(gg(idim))
  end do
  grav = grav*dx/ctot**2
  grav = MAX(grav,0.0001d0)

  ! Compute maximum time step for each authorized cell
  dt = r%courant_factor*dx/r%smallc
  dtcell = two*r%courant_factor*dx/(ctot*(one+sqrt(one+two*r%courant_factor*grav)))
  dt = min(dt,dtcell)

end subroutine cmpdt
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine riemann_llf(qleft,qright,fgdnv,params)
  use amr_parameters, only: ndim
  use hydro_parameters, only: nprim, nener
  use amr_commons, only: hydro_params_t
  use const
  implicit none
  ! dummy arguments
  type(hydro_params_t)::params
  real(kind=8),dimension(1:nprim)::qleft,qright
  real(kind=8),dimension(1:nprim+1)::fgdnv

  real(kind=8),dimension(1:nprim+1)::fleft,fright
  real(kind=8),dimension(1:nprim+1)::uleft,uright
  real(kind=8)::cmax
  integer::n
  real(kind=8)::gamma,smallr,smallc,smallp, entho
  real(kind=8),dimension(1:nener+1)::gamma_rad
  real(kind=8)::rl,ul,pl,cl
  real(kind=8)::rr,ur,pr,cr

  ! Constants
  gamma = params%gamma
  gamma_rad = params%gamma_rad
  smallr = params%smallr
  smallc = params%smallc
  smallp = smallc**2/gamma
  entho = one/(gamma-one)

  !===========================
  ! Compute maximum wave speed
  !===========================
  ! Left states
  rl = max(qleft(1),smallr)
  ul =     qleft(2)
  pl = max(qleft(5),rl*smallp)
  cl = gamma*pl
#if NENER>0
  do n = 1,nener
     cl = cl + gamma_rad(n)*qleft(5+n)
  end do
#endif
  cl = sqrt(cl/rl)
  ! Right states
  rr = max(qright(1),smallr)
  ur =     qright(2)
  pr = max(qright(5),rr*smallp)
  cr = gamma*pr
#if NENER>0
  do n = 1,nener
     cr = cr + gamma_rad(n)*qright(5+n)
  end do
#endif
  cr = sqrt(cr/rr)
  ! Local max. wave speed
  cmax = max(abs(ul)+cl,abs(ur)+cr)

  !===============================
  ! Compute conservative variables
  !===============================
  ! Mass density
  uleft (1) = qleft (1)
  uright(1) = qright(1)
  ! Normal momentum
  uleft (2) = qleft (1)*qleft (2)
  uright(2) = qright(1)*qright(2)
  ! Transverse velocities
  do n = 3, 4
     uleft (n) = qleft (1)*qleft (n)
     uright(n) = qright(1)*qright(n)
  end do
  ! Total energy
  uleft (5) = qleft (5)*entho + half*qleft (1)*qleft (2)**2
  uright(5) = qright(5)*entho + half*qright(1)*qright(2)**2
  uleft (5) = uleft (5)       + half*qleft (1)*qleft (3)**2
  uright(5) = uright(5)       + half*qright(1)*qright(3)**2
  uleft (5) = uleft (5)       + half*qleft (1)*qleft (4)**2
  uright(5) = uright(5)       + half*qright(1)*qright(4)**2
#if NENER>0
  do n = 1,nener
     uleft (5) = uleft (5) + qleft (5+n)/(gamma_rad(n)-one)
     uright(5) = uright(5) + qright(5+n)/(gamma_rad(n)-one)
  end do
#endif
  ! Non-thermal energies
#if NENER>0
  do n = 1, nener
     uleft (5+n) = qleft (5+n)/(gamma_rad(n)-one)
     uright(5+n) = qright(5+n)/(gamma_rad(n)-one)
  end do
#endif
  ! Other passively advected quantities
#if NVAR>5+NENER
  do n = 6+nener, nprim
     uleft (n) = qleft (1)*qleft (n)
     uright(n) = qright(1)*qright(n)
  end do
#endif
  ! Thermal energy
  uleft (nprim+1) = qleft (5)*entho
  uright(nprim+1) = qright(5)*entho

  !==============================
  ! Compute left and right fluxes
  !==============================
  ! Mass density
  fleft (1) = qleft (2)*uleft (1)
  fright(1) = qright(2)*uright(1)
  ! Normal momentum
  fleft (2) = qleft (2)*uleft (2) + qleft (5)
  fright(2) = qright(2)*uright(2) + qright(5)
#if NENER>0
  do n = 1,nener
     fleft (2) = fleft (2) + qleft (5+n)
     fright(2) = fright(2) + qright(5+n)
  end do
#endif
  ! Transverse momentum
  do n = 3,4
     fleft (n) = qleft (2)*uleft (n)
     fright(n) = qright(2)*uright(n)
  end do
  ! Total energy
  fleft (5) = qleft (2)*(uleft (5)+qleft (5))
  fright(5) = qright(2)*(uright(5)+qright(5))
#if NENER>0
  do n = 1,nener
     fleft (5) = fleft (5) + qleft (2)*qleft (5+n)
     fright(5) = fright(5) + qright(2)*qright(5+n)
  end do
#endif
  ! Other passively advected quantities
  do n = 6, nprim+1
     fleft (n) = qleft (2)*uleft (n)
     fright(n) = qright(2)*uright(n)
  end do

  !=============================
  ! Compute Lax-Friedrich fluxes
  !=============================
  do n = 1, nprim+1
     fgdnv(n) = half*(fleft(n)+fright(n)-cmax*(uright(n)-uleft(n)))
  end do

end subroutine riemann_llf
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine riemann_hll(qleft,qright,fgdnv,params)
  use amr_parameters, only: ndim
  use hydro_parameters, only: nprim, nener
  use amr_commons, only: hydro_params_t
  use const
  implicit none
  ! 1D HLL Riemann solver

  ! dummy arguments
  type(hydro_params_t)::params
  real(kind=8),dimension(1:nprim)::qleft,qright
  real(kind=8),dimension(1:nprim+1)::fgdnv

  real(kind=8),dimension(1:nprim+1)::fleft,fright
  real(kind=8),dimension(1:nprim+1)::uleft,uright
  real(kind=8)::SL,SR
  integer::n
  real(kind=8)::gamma,smallr,smallc,smallp, entho
  real(kind=8),dimension(1:nener+1)::gamma_rad
  real(kind=8)::rl,ul,pl,cl
  real(kind=8)::rr,ur,pr,cr

  ! Constants
  gamma = params%gamma
  gamma_rad = params%gamma_rad
  smallr = params%smallr
  smallc = params%smallc
  smallp = smallc**2/gamma
  entho = one/(gamma-one)

  !===========================
  ! Compute maximum wave speed
  !===========================
  ! Left states
  rl = max(qleft (1),smallr)
  ul =     qleft (2)
  pl = max(qleft (5),rl*smallp)
  cl = gamma*pl
#if NENER>0
  do n = 1,nener
     cl = cl + gamma_rad(n)*qleft(5+n)
  end do
#endif
  cl = sqrt(cl/rl)
  ! Right states
  rr = max(qright(1),smallr)
  ur =     qright(2)
  pr = max(qright(5),rr*smallp)
  cr = gamma*pr
#if NENER>0
  do n = 1,nener
     cr = cr + gamma_rad(n)*qright(5+n)
  end do
#endif
  cr = sqrt(cr/rr)
  ! Left and right max. wave speed
  SL=min(min(ul,ur)-max(cl,cr),zero)
  SR=max(max(ul,ur)+max(cl,cr),zero)

  !===============================
  ! Compute conservative variables
  !===============================
  ! Mass density
  uleft (1) = qleft (1)
  uright(1) = qright(1)
  ! Normal momentum
  uleft (2) = qleft (1)*qleft (2)
  uright(2) = qright(1)*qright(2)
  ! Transverse momentum
  do n = 3, 4
     uleft (n) = qleft (1)*qleft (n)
     uright(n) = qright(1)*qright(n)
  end do
  ! Total energy
  uleft (5) = qleft (5)*entho + half*qleft (1)*qleft (2)**2
  uright(5) = qright(5)*entho + half*qright(1)*qright(2)**2
  uleft (5) = uleft (5)       + half*qleft (1)*qleft (3)**2
  uright(5) = uright(5)       + half*qright(1)*qright(3)**2
  uleft (5) = uleft (5)       + half*qleft (1)*qleft (4)**2
  uright(5) = uright(5)       + half*qright(1)*qright(4)**2
#if NENER>0
  do n = 1,nener
     uleft (5) = uleft (5) + qleft (5+n)/(gamma_rad(n)-one)
     uright(5) = uright(5) + qright(5+n)/(gamma_rad(n)-one)
  end do
#endif
  ! Non-thermal energies
#if NENER>0
  do n = 1, nener
     uleft (5+n) = qleft (5+n)/(gamma_rad(n)-one)
     uright(5+n) = qright(5+n)/(gamma_rad(n)-one)
  end do
#endif
  ! Other passively advected quantities
#if NVAR>5+NENER
  do n = 6+nener, nprim
     uleft (n) = qleft (1)*qleft (n)
     uright(n) = qright(1)*qright(n)
  end do
#endif
  ! Thermal energy
  uleft (nprim+1) = qleft (5)*entho
  uright(nprim+1) = qright(5)*entho

  !==============================
  ! Compute left and right fluxes
  !==============================
  ! Mass density
  fleft (1) = uleft (2)
  fright(1) = uright(2)
  ! Normal momentum
  fleft (2) = qleft (5)+uleft (2)*qleft (2)
  fright(2) = qright(5)+uright(2)*qright(2)
#if NENER>0
  do n = 1,nener
     fleft (2) = fleft (2) + qleft (5+n)
     fright(2) = fright(2) + qright(5+n)
  end do
#endif
  ! Transverse momentum
  do n = 3, 4
     fleft (n) = qleft (2)*uleft (n)
     fright(n) = qright(2)*uright(n)
  end do
  ! Total energy
  fleft (5) = qleft (2)*(uleft (5)+qleft (5))
  fright(5) = qright(2)*(uright(5)+qright(5))
#if NENER>0
  do n = 1,nener
     fleft (5) = fleft (5) + qleft (2)*qleft (5+n)
     fright(5) = fright(5) + qright(2)*qright(5+n)
  end do
#endif
  ! Other advected quantities
  do n = 6, nprim+1
     fleft (n) = qleft (2)*uleft (n)
     fright(n) = qright(2)*uright(n)
  end do

  !===================
  ! Compute HLL fluxes
  !===================
  do n = 1, nprim+1
     fgdnv(n) = (SR*fleft(n)-SL*fright(n) &
          & + SR*SL*(uright(n)-uleft(n)))/(SR-SL)
  end do

end subroutine riemann_hll
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine riemann_hllc(qleft,qright,fgdnv,params)
  use amr_parameters, only: ndim
  use hydro_parameters, only: nprim, nener
  use amr_commons, only: hydro_params_t
  use const
  implicit none
  ! HLLC Riemann solver (Toro)

  ! dummy arguments
  type(hydro_params_t)::params
  real(kind=8),dimension(1:nprim)::qleft,qright
  real(kind=8),dimension(1:nprim+1)::fgdnv

  real(kind=8)::SL,SR
  real(kind=8)::entho
  real(kind=8)::rl,pl,ul,vl,wl,ekinl,etotl,el,ptotl
  real(kind=8)::rr,pr,ur,vr,wr,ekinr,etotr,er,ptotr
  real(kind=8)::cfastl,rcl,rstarl,estarl
  real(kind=8)::cfastr,rcr,rstarr,estarr
  real(kind=8)::etotstarl,etotstarr
  real(kind=8)::ustar,ptotstar
  real(kind=8)::ro,uo,vo,wo,ptoto,etoto,eo
  real(kind=8)::gamma,smallr,smallc,smallp
  real(kind=8),dimension(1:nener+1)::gamma_rad
#if NENER>0
  INTEGER::irad
  real(kind=8),dimension(1:nener)::eradl,eradr,erado
  real(kind=8),dimension(1:nener)::eradstarl,eradstarr
#endif
#if NVAR>5+NENER
  INTEGER::ivar
#endif

  ! constants
  gamma = params%gamma
  gamma_rad = params%gamma_rad
  smallr = params%smallr
  smallc = params%smallc
  smallp = smallc**2/gamma
  entho = one/(gamma-one)

  ! Left variables
  rl=max(qleft(1),smallr)
  ul=qleft(2); vl=qleft(3); wl=qleft(4)
  Pl=max(qleft(5),rl*smallp)

  el=Pl*entho
  ekinl=half*rl*(ul*ul+vl*vl+wl*wl)
  etotl=el+ekinl
#if NENER>0
  do irad=1,nener
     eradl(irad)=qleft(5+irad)/(gamma_rad(irad)-one)
     etotl=etotl+eradl(irad)
  end do
#endif
  Ptotl = Pl
#if NENER>0
  do irad=1,nener
     Ptotl=Ptotl+qleft(5+irad)
  end do
#endif

  ! Right variables
  rr=max(qright(1),smallr)
  ur=qright(2); vr=qright(3); wr=qright(4)
  Pr=max(qright(5),rr*smallp)

  er=Pr*entho
  ekinr=half*rr*(ur*ur+vr*vr+wr*wr)
  etotr=er+ekinr
#if NENER>0
  do irad=1,nener
     eradr(irad)=qright(5+irad)/(gamma_rad(irad)-one)
     etotr=etotr+eradr(irad)
  end do
#endif
  Ptotr = Pr
#if NENER>0
  do irad=1,nener
     Ptotr=Ptotr+qright(5+irad)
  end do
#endif

  ! Find the largest eigenvalues in the normal direction to the interface
  cfastl=gamma*Pl
#if NENER>0
  do irad = 1,nener
     cfastl = cfastl + gamma_rad(irad)*qleft(5+irad)
  end do
#endif
  cfastl=sqrt(max(cfastl/rl,smallc**2))

  cfastr=gamma*Pr
#if NENER>0
  do irad = 1,nener
     cfastr = cfastr + gamma_rad(irad)*qright(5+irad)
  end do
#endif
  cfastr=sqrt(max(cfastr/rr,smallc**2))

  ! Compute HLL wave speed
  SL=min(ul,ur)-max(cfastl,cfastr)
  SR=max(ul,ur)+max(cfastl,cfastr)

  ! Compute lagrangian sound speed
  rcl=rl*(ul-SL)
  rcr=rr*(SR-ur)

  ! Compute acoustic star state
  ustar   =(rcr*ur   +rcl*ul   +  (Ptotl-Ptotr))/(rcr+rcl)
  Ptotstar=(rcr*Ptotl+rcl*Ptotr+rcl*rcr*(ul-ur))/(rcr+rcl)

  ! Left star region variables
  rstarl=rl*(SL-ul)/(SL-ustar)
  etotstarl=((SL-ul)*etotl-Ptotl*ul+Ptotstar*ustar)/(SL-ustar)
  estarl=el*(SL-ul)/(SL-ustar)
#if NENER>0
  do irad = 1,nener
     eradstarl(irad)=eradl(irad)*(SL-ul)/(SL-ustar)
  end do
#endif

  ! Right star region variables
  rstarr=rr*(SR-ur)/(SR-ustar)
  etotstarr=((SR-ur)*etotr-Ptotr*ur+Ptotstar*ustar)/(SR-ustar)
  estarr=er*(SR-ur)/(SR-ustar)
#if NENER>0
  do irad = 1,nener
     eradstarr(irad)=eradr(irad)*(SR-ur)/(SR-ustar)
  end do
#endif

  ! Sample the solution at x/t=0
  if(SL>0d0)then
     ro=rl
     uo=ul
     vo=vl
     wo=wl
     Ptoto=Ptotl
     etoto=etotl
     eo=el
#if NENER>0
     do irad = 1,nener
        erado(irad)=eradl(irad)
     end do
#endif
  else if(ustar>0d0)then
     ro=rstarl
     uo=ustar
     vo=vl
     wo=wl
     Ptoto=Ptotstar
     etoto=etotstarl
     eo=estarl
#if NENER>0
     do irad = 1,nener
        erado(irad)=eradstarl(irad)
     end do
#endif
  else if (SR>0d0)then
     ro=rstarr
     uo=ustar
     vo=vr
     wo=wr
     Ptoto=Ptotstar
     etoto=etotstarr
     eo=estarr
#if NENER>0
     do irad = 1,nener
        erado(irad)=eradstarr(irad)
     end do
#endif
  else
     ro=rr
     uo=ur
     vo=vr
     wo=wr
     Ptoto=Ptotr
     etoto=etotr
     eo=er
#if NENER>0
     do irad = 1,nener
        erado(irad)=eradr(irad)
     end do
#endif
  end if

  !=========================
  ! Compute the Godunov flux
  !=========================
  fgdnv(1) = ro*uo
  fgdnv(2) = ro*uo*uo+Ptoto
  fgdnv(3) = ro*uo*vo
  fgdnv(4) = ro*uo*wo
  fgdnv(5) = (etoto+Ptoto)*uo
  ! Non-thermal energies
#if NENER>0
  do irad = 1,nener
     fgdnv(5+irad) = uo*erado(irad)
  end do
#endif
  ! Other passively advected quantities
#if NVAR>5+NENER
  do ivar = 6+nener,nprim
     if(ustar>0)then
        fgdnv(ivar) = ro*uo*qleft (ivar)
     else
        fgdnv(ivar) = ro*uo*qright(ivar)
     endif
  end do
#endif
  ! Thermal energy
  fgdnv(nprim+1) = uo*eo

end subroutine riemann_hllc
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#ifdef MHD
SUBROUTINE riemann_upwind_mhd(qleft,qright,fgdnv,zero_flux,params)
  use amr_parameters, only: ndim
  use hydro_parameters, only: nprim, nener
  use amr_commons, only: hydro_params_t
  use const
  ! 1D Upwind Riemann solver
  IMPLICIT NONE
  type(hydro_params_t)::params
  real(kind=8),DIMENSION(1:nprim)::qleft,qright
  real(kind=8),DIMENSION(1:nprim+1)::fgdnv
  real(kind=8)::zero_flux

  real(kind=8),DIMENSION(1:nprim+1)::fleft,fright,fmean
  real(kind=8),DIMENSION(1:nprim+1)::uleft,uright,udiff
  real(kind=8):: vleft,bx_mean

  ! Enforce continuity of normal component
  bx_mean=half*(qleft(6)+qright(6))
  qleft (6)=bx_mean
  qright(6)=bx_mean

  CALL find_mhd_flux(qleft ,uleft ,fleft ,params)
  CALL find_mhd_flux(qright,uright,fright,params)

  ! find the mean flux
  fmean =  half * ( fright + fleft ) * zero_flux

  ! find the mean normal velocity
  vleft = half * ( qleft(2) + qright(2) )

  ! difference between the 2 states
  udiff = half * ( uright - uleft )

  ! the Upwind flux
  fgdnv = fmean - ABS(vleft) * udiff

END SUBROUTINE riemann_upwind_mhd
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE riemann_llf_mhd(qleft,qright,fgdnv,zero_flux,params)
  use amr_parameters, only: ndim
  use hydro_parameters, only: nprim, nener
  use amr_commons, only: hydro_params_t
  use const
  IMPLICIT NONE
  ! 1D local Lax-Friedrich Riemann solver
  type(hydro_params_t)::params
  real(kind=8),DIMENSION(1:nprim)::qleft,qright
  real(kind=8),DIMENSION(1:nprim+1)::fgdnv
  real(kind=8)::zero_flux

  real(kind=8),DIMENSION(1:nprim+1)::fleft,fright,fmean
  real(kind=8),DIMENSION(1:nprim+1)::uleft,uright,udiff
  real(kind=8)::vleft,vright,bx_mean

  ! Enforce continuity of normal component
  bx_mean=half*(qleft(6)+qright(6))
  qleft (6)=bx_mean
  qright(6)=bx_mean

  CALL find_mhd_flux(qleft ,uleft ,fleft ,params)
  CALL find_mhd_flux(qright,uright,fright,params)

  ! find the mean flux
  fmean =  half * ( fright + fleft ) * zero_flux

  ! find the largest eigenvalue in the normal direction to the interface
  CALL find_speed_info(qleft ,vleft ,params)
  CALL find_speed_info(qright,vright,params)

  ! difference between the 2 states
  udiff  = half * ( uright - uleft )

  ! the local Lax-Friedrich flux
  fgdnv = fmean - MAX(vleft,vright) * udiff

END SUBROUTINE riemann_llf_mhd
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE riemann_hll_mhd(qleft,qright,fgdnv,params)
  use amr_parameters, only: ndim
  use hydro_parameters, only: nprim, nener
  use amr_commons, only: hydro_params_t
  use const
  IMPLICIT NONE
  ! 1D HLL Riemann solver
  type(hydro_params_t)::params
  real(kind=8),DIMENSION(1:nprim)::qleft,qright
  real(kind=8),DIMENSION(1:nprim+1)::fgdnv

  real(kind=8),DIMENSION(1:nprim+1)::fleft,fright
  real(kind=8),DIMENSION(1:nprim+1)::uleft,uright
  real(kind=8)::vleft,vright,bx_mean,cfleft,cfright,SL,SR

  ! Enforce continuity of normal component
  bx_mean=half*(qleft(6)+qright(6))
  qleft (6)=bx_mean
  qright(6)=bx_mean

  CALL find_mhd_flux(qleft ,uleft ,fleft ,params)
  CALL find_mhd_flux(qright,uright,fright,params)

  ! find the largest eigenvalue in the normal direction to the interface
  CALL find_speed_fast(qleft ,cfleft ,params)
  CALL find_speed_fast(qright,cfright,params)
  vleft =qleft (2)
  vright=qright(2)
  SL=min(min(vleft,vright)-max(cfleft,cfright),zero)
  SR=max(max(vleft,vright)+max(cfleft,cfright),zero)

  ! the HLL flux
  fgdnv = (SR*fleft-SL*fright+SR*SL*(uright-uleft))/(SR-SL)

END SUBROUTINE riemann_hll_mhd
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE riemann_hlld(qleft,qright,fgdnv,params)
  use amr_parameters, only: ndim
  use hydro_parameters, only: nprim, nener
  use amr_commons, only: hydro_params_t
  use const
  IMPLICIT NONE
  ! HLLD Riemann solver (Miyoshi & Kusano, 2005, JCP, 208, 315)
  type(hydro_params_t)::params
  real(kind=8),DIMENSION(1:nprim)::qleft,qright
  real(kind=8),DIMENSION(1:nprim+1)::fgdnv

  real(kind=8)::SL,SR,SAL,SAR
  real(kind=8)::gamma,smallc,entho,A,sgnm
  real(kind=8),dimension(1:nener+1)::gamma_rad
  real(kind=8)::rl,pl,ul,vl,wl,cl,ekinl,emagl,etotl,ptotl,vdotbl,bl,el
  real(kind=8)::rr,pr,ur,vr,wr,cr,ekinr,emagr,etotr,ptotr,vdotbr,br,er
  real(kind=8)::cfastl,calfvenl,rcl,rstarl,vstarl,wstarl,bstarl,cstarl,vdotbstarl
  real(kind=8)::cfastr,calfvenr,rcr,rstarr,vstarr,wstarr,bstarr,cstarr,vdotbstarr
  real(kind=8)::sqrrstarl,etotstarl,etotstarstarl
  real(kind=8)::sqrrstarr,etotstarr,etotstarstarr
  real(kind=8)::ustar,ptotstar,estar,vstarstar,wstarstar,bstarstar,cstarstar,vdotbstarstar
  real(kind=8)::ro,uo,vo,wo,bo,co,ptoto,etoto,vdotbo
  real(kind=8)::einto,eintl,eintr,eintstarr,eintstarl

#if NVAR>5+NENER
  INTEGER ::ivar
#endif
#if NENER>0
  INTEGER ::irad
  real(kind=8),dimension(1:nener)::erado,eradl,eradr
  real(kind=8),dimension(1:nener)::eradstarl,eradstarr
#endif

  gamma = params%gamma
  gamma_rad = params%gamma_rad
  smallc = params%smallc
  entho = one/(gamma-one)

  ! Enforce continuity of normal component
  A=half*(qleft(6)+qright(6))
  sgnm=sign(one,A)
  qleft(6)=A; qright(6)=A

  ! Left variables
  rl=qleft(1); ul=qleft(2); vl=qleft(3); wl=qleft(4);
  Pl=qleft(5); Bl=qleft(7); Cl=qleft(8)
  ekinl = half*(ul*ul+vl*vl+wl*wl)*rl
  emagl = half*(A*A+Bl*Bl+Cl*Cl)
  etotl = Pl*entho+ekinl+emagl
  Ptotl = Pl + emagl
  vdotBl= ul*A+vl*Bl+wl*cl
#if NENER>0
  do irad = 1,nener
     eradl(irad) = qleft(8+irad)/(gamma_rad(irad)-1.0d0)
     etotl = etotl + eradl(irad)
     Ptotl = Ptotl + qleft(8+irad)
  end do
#endif
  eintl=Pl*entho

  ! Right variables
  rr=qright(1); ur=qright(2); vr=qright(3); wr=qright(4)
  Pr=qright(5); Br=qright(7); Cr=qright(8)
  ekinr = half*(ur*ur+vr*vr+wr*wr)*rr
  emagr = half*(A*A+Br*Br+Cr*Cr)
  etotr = Pr*entho+ekinr+emagr
  Ptotr = Pr + emagr
  vdotBr= ur*A+vr*Br+wr*Cr
#if NENER>0
  do irad = 1,nener
     eradr(irad) = qright(8+irad)/(gamma_rad(irad)-1.0d0)
     etotr = etotr + eradr(irad)
     Ptotr = Ptotr + qright(8+irad)
  end do
#endif
  eintr=Pr*entho

  ! Find the largest eigenvalues in the normal direction to the interface
  CALL find_speed_fast(qleft ,cfastl,params)
  CALL find_speed_fast(qright,cfastr,params)

  ! Compute HLL wave speed
  SL=min(ul,ur)-max(cfastl,cfastr)
  SR=max(ul,ur)+max(cfastl,cfastr)
!  SL=ul-cfastl
!  SR=ur+cfastr

  ! Compute lagrangian sound speed
  rcl=rl*(ul-SL)
  rcr=rr*(SR-ur)

  ! Compute acoustic star state
  ustar   =(rcr*ur   +rcl*ul   +  (Ptotl-Ptotr))/(rcr+rcl)
  Ptotstar=(rcr*Ptotl+rcl*Ptotr+rcl*rcr*(ul-ur))/(rcr+rcl)

  ! Left star region variables
  rstarl=rl*(SL-ul)/(SL-ustar)
  estar =rl*(SL-ul)*(SL-ustar)-A**2
  el    =rl*(SL-ul)*(SL-ul   )-A**2
#if NENER>0
  do irad = 1,nener
     eradstarl(irad)=eradl(irad)*(SL-ul)/(SL-ustar)
  end do
#endif
  eintstarl=eintl*(SL-ul)/(SL-ustar)
  if(abs(estar)<1e-4*A**2)then
     vstarl=vl
     Bstarl=Bl
     wstarl=wl
     Cstarl=Cl
  else
     vstarl=vl-A*Bl*(ustar-ul)/estar
     Bstarl=Bl*el/estar
     wstarl=wl-A*Cl*(ustar-ul)/estar
     Cstarl=Cl*el/estar
  endif
  vdotBstarl=ustar*A+vstarl*Bstarl+wstarl*Cstarl
  etotstarl=((SL-ul)*etotl-Ptotl*ul+Ptotstar*ustar+A*(vdotBl-vdotBstarl))/(SL-ustar)
  sqrrstarl=sqrt(rstarl)
  calfvenl=abs(A)/sqrrstarl
  SAL=ustar-calfvenl

  ! Right star region variables
  rstarr=rr*(SR-ur)/(SR-ustar)
  estar =rr*(SR-ur)*(SR-ustar)-A**2
  er    =rr*(SR-ur)*(SR-ur   )-A**2
#if NENER>0
  do irad = 1,nener
     eradstarr(irad)=eradr(irad)*(SR-ur)/(SR-ustar)
  end do
#endif
  eintstarr=eintr*(SR-ur)/(SR-ustar)
  if(abs(estar)<1e-4*A**2)then
     vstarr=vr
     Bstarr=Br
     wstarr=wr
     Cstarr=Cr
  else
     vstarr=vr-A*Br*(ustar-ur)/estar
     Bstarr=Br*er/estar
     wstarr=wr-A*Cr*(ustar-ur)/estar
     Cstarr=Cr*er/estar
  endif
  vdotBstarr=ustar*A+vstarr*Bstarr+wstarr*Cstarr
  etotstarr=((SR-ur)*etotr-Ptotr*ur+Ptotstar*ustar+A*(vdotBr-vdotBstarr))/(SR-ustar)
  sqrrstarr=sqrt(rstarr)
  calfvenr=abs(A)/sqrrstarr
  SAR=ustar+calfvenr

  ! Double star region variables
  vstarstar=(sqrrstarl*vstarl+sqrrstarr*vstarr+sgnm*(Bstarr-Bstarl))/(sqrrstarl+sqrrstarr)
  wstarstar=(sqrrstarl*wstarl+sqrrstarr*wstarr+sgnm*(Cstarr-Cstarl))/(sqrrstarl+sqrrstarr)
  Bstarstar=(sqrrstarl*Bstarr+sqrrstarr*Bstarl+sgnm*sqrrstarl*sqrrstarr*(vstarr-vstarl))/(sqrrstarl+sqrrstarr)
  Cstarstar=(sqrrstarl*Cstarr+sqrrstarr*Cstarl+sgnm*sqrrstarl*sqrrstarr*(wstarr-wstarl))/(sqrrstarl+sqrrstarr)
  vdotBstarstar=ustar*A+vstarstar*Bstarstar+wstarstar*Cstarstar
  etotstarstarl=etotstarl-sgnm*sqrrstarl*(vdotBstarl-vdotBstarstar)
  etotstarstarr=etotstarr+sgnm*sqrrstarr*(vdotBstarr-vdotBstarstar)

  ! Sample the solution at x/t=0
  if(SL>0d0)then
     ro=rl
     uo=ul
     vo=vl
     wo=wl
     Bo=Bl
     Co=Cl
     Ptoto=Ptotl
     etoto=etotl
     vdotBo=vdotBl
#if NENER>0
     do irad = 1,nener
        erado(irad)=eradl(irad)
     end do
#endif
     einto=eintl
  else if(SAL>0d0)then
     ro=rstarl
     uo=ustar
     vo=vstarl
     wo=wstarl
     Bo=Bstarl
     Co=Cstarl
     Ptoto=Ptotstar
     etoto=etotstarl
     vdotBo=vdotBstarl
#if NENER>0
     do irad = 1,nener
        erado(irad)=eradstarl(irad)
     end do
#endif
     einto=eintstarl
  else if(ustar>0d0)then
     ro=rstarl
     uo=ustar
     vo=vstarstar
     wo=wstarstar
     Bo=Bstarstar
     Co=Cstarstar
     Ptoto=Ptotstar
     etoto=etotstarstarl
     vdotBo=vdotBstarstar
#if NENER>0
     do irad = 1,nener
        erado(irad)=eradstarl(irad)
     end do
#endif
     einto=eintstarl
  else if(SAR>0d0)then
     ro=rstarr
     uo=ustar
     vo=vstarstar
     wo=wstarstar
     Bo=Bstarstar
     Co=Cstarstar
     Ptoto=Ptotstar
     etoto=etotstarstarr
     vdotBo=vdotBstarstar
#if NENER>0
     do irad = 1,nener
        erado(irad)=eradstarr(irad)
     end do
#endif
     einto=eintstarr
  else if (SR>0d0)then
     ro=rstarr
     uo=ustar
     vo=vstarr
     wo=wstarr
     Bo=Bstarr
     Co=Cstarr
     Ptoto=Ptotstar
     etoto=etotstarr
     vdotBo=vdotBstarr
#if NENER>0
     do irad = 1,nener
        erado(irad)=eradstarr(irad)
     end do
#endif
     einto=eintstarr
  else
     ro=rr
     uo=ur
     vo=vr
     wo=wr
     Bo=Br
     Co=Cr
     Ptoto=Ptotr
     etoto=etotr
     vdotBo=vdotBr
#if NENER>0
     do irad = 1,nener
        erado(irad)=eradr(irad)
     end do
#endif
     einto=eintr
  end if

  ! Compute the Godunov flux
  fgdnv(1) = ro*uo
  fgdnv(2) = ro*uo*uo+Ptoto-A*A
  fgdnv(3) = ro*uo*vo-A*Bo
  fgdnv(4) = ro*uo*wo-A*Co
  fgdnv(5) = (etoto+Ptoto)*uo-A*vdotBo
  fgdnv(6) = zero
  fgdnv(7) = Bo*uo-A*vo
  fgdnv(8) = Co*uo-A*wo
#if NENER>0
  do irad = 1,nener
     fgdnv(8+irad) = uo*erado(irad)
  end do
#endif
#if NVAR>5+NENER
  do ivar = 9+nener,nprim
     if(fgdnv(1)>0)then
        fgdnv(ivar) = fgdnv(1)*qleft (ivar)
     else
        fgdnv(ivar) = fgdnv(1)*qright(ivar)
     endif
  end do
#endif
  !Thermal energy
  fgdnv(nprim+1) = uo*einto
END SUBROUTINE riemann_hlld
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE find_mhd_flux(qvar,cvar,ff,params)
  use amr_parameters, only: ndim
  use hydro_parameters, only: nprim, nener
  use amr_commons, only: hydro_params_t
  use const
  IMPLICIT NONE
  !! compute the 1D MHD fluxes from the conservative variables
  !! the structure of qvar is : rho, Pressure, Vnormal, Bnormal,
  !! Vtransverse1, Btransverse1, Vtransverse2, Btransverse2
  type(hydro_params_t)::params
  real(kind=8),DIMENSION(1:nprim)::qvar
  real(kind=8),DIMENSION(1:nprim+1)::cvar,ff

  real(kind=8)::ekin,emag,etot,d,u,v,w,A,B,C,P,Ptot,entho
  real(kind=8)::gamma
  real(kind=8),dimension(1:nener+1)::gamma_rad
#if NVAR>5+NENER
  INTEGER::ivar
#endif
#if NENER>0
  INTEGER::irad
#endif

  ! Constants
  gamma = params%gamma
  gamma_rad = params%gamma_rad
  entho = one/(gamma-one)

  ! Local variables
  d=qvar(1); u=qvar(2); v=qvar(3); w=qvar(4)
  P=qvar(5); A=qvar(6); B=qvar(7); C=qvar(8)
  ekin = half*(u*u+v*v+w*w)*d
  emag = half*(A*A+B*B+C*C)
  etot = P*entho + ekin + emag
  Ptot = P + emag
#if NENER>0
  do irad = 1,nener
     etot = etot + qvar(8+irad)/(gamma_rad(irad)-one)
     Ptot = Ptot + qvar(8+irad)
  end do
#endif

  ! Compute conservative variables
  cvar(1) = d
  cvar(2) = d*u
  cvar(3) = d*v
  cvar(4) = d*w
  cvar(5) = etot
  cvar(6) = A
  cvar(7) = B
  cvar(8) = C
#if NENER>0
  do irad = 1,nener
     cvar(8+irad) = qvar(8+irad)/(gamma_rad(irad)-one)
  end do
#endif
#if NVAR>5+NENER
  do ivar = 9+nener,nprim
     cvar(ivar) = d*qvar(ivar)
  end do
#endif
  ! Thermal energy
  cvar(nprim+1)=P*entho

  ! Compute fluxes
  ff(1) = d*u
  ff(2) = d*u*u+Ptot-A*A
  ff(3) = d*u*v-A*B
  ff(4) = d*u*w-A*C
  ff(5) = (etot+Ptot)*u-A*(A*u+B*v+C*w)
  ff(6) = zero
  ff(7) = B*u-A*v
  ff(8) = C*u-A*w
#if NENER>0
  do irad = 1,nener
     ff(8+irad) = u*cvar(8+irad)
  end do
#endif
#if NVAR>5+NENER
  do ivar = 9+nener,nprim
     ff(ivar) = d*u*qvar(ivar)
  end do
#endif
  ! Thermal energy
  ff(nprim+1)=P*entho*u

END SUBROUTINE find_mhd_flux
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE find_speed_info(qvar,vel_info,params)
  use amr_parameters, only: ndim
  use hydro_parameters, only: nprim, nener
  use amr_commons, only: hydro_params_t
  use const
  IMPLICIT NONE
  !! calculate the fastest velocity at which information is exchanged
  !! at the interface: we mean c_fast + abs(v).
  !! the structure of qvar is : rho, Pressure, Vnormal, Bnormal,
  !! Vtransverse1,Btransverse1,Vtransverse2,Btransverse2
  type(hydro_params_t)::params
  real(kind=8),DIMENSION(1:nprim)::qvar
  real(kind=8)::vel_info

  real(kind=8)::gamma
  real(kind=8),dimension(1:nener+1)::gamma_rad
  real(kind=8)::d,P,u,A,B,C,B2,c2,d2,cf
#if NENER>0
  INTEGER::irad
#endif

  ! Constants
  gamma = params%gamma
  gamma_rad = params%gamma_rad

  d=qvar(1); u=qvar(2); P=qvar(5)
  A=qvar(6); B=qvar(7); C=qvar(8)
  B2 = A*A+B*B+C*C
  c2 = gamma*P/d
#if NENER>0
  do irad = 1,nener
     c2 = c2 + gamma_rad(irad)*qvar(8+irad)/d
  end do
#endif

  d2 = half*(B2/d+c2)
  cf = sqrt( d2 + sqrt(d2**2-c2*A*A/d) )
  vel_info = cf+abs(u)

END SUBROUTINE find_speed_info
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE find_speed_fast(qvar,vel_info,params)
  use amr_parameters, only: ndim
  use hydro_parameters, only: nprim, nener
  use amr_commons, only: hydro_params_t
  use const
  IMPLICIT NONE
  !! calculate the fast magnetosonic velocity
  !! the structure of qvar is : rho, Pressure, Vnormal, Bnormal,
  !! Vtransverse1,Btransverse1,Vtransverse2,Btransverse2
  type(hydro_params_t)::params
  real(kind=8),DIMENSION(1:nprim)::qvar
  real(kind=8)::vel_info

  real(kind=8)::d,P,A,B,C,B2,c2,d2,cf
  real(kind=8)::gamma
  real(kind=8),dimension(1:nener+1)::gamma_rad
#if NENER>0
  INTEGER::irad
#endif

  ! Constants
  gamma = params%gamma
  gamma_rad = params%gamma_rad

  d=qvar(1); P=qvar(5); A=qvar(6); B=qvar(7); C=qvar(8)
  B2 = A*A+B*B+C*C
  c2 = gamma*P/d
#if NENER>0
  do irad = 1,nener
     c2 = c2 + gamma_rad(irad)*qvar(8+irad)/d
  end do
#endif
  d2 = half*(B2/d+c2)
  cf = sqrt( d2 + sqrt(d2**2-c2*A*A/d) )
  vel_info = cf

END SUBROUTINE find_speed_fast
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE find_speed_alfven(qvar,vel_info)
  USE amr_parameters
  USE const
  USE hydro_parameters
  !! calculate the alfven velocity
  !! the structure of qvar is : rho, Pressure, Vnormal, Bnormal,
  !! Vtransverse1,Btransverse1,Vtransverse2,Btransverse2
  IMPLICIT NONE

  real(kind=8),DIMENSION(1:nprim):: qvar
  real(kind=8)::vel_info

  real(kind=8)::d,A

  d=qvar(1); A=qvar(6)
  vel_info = sqrt(A*A/d)

END SUBROUTINE find_speed_alfven
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE riemann_roe_mhd(qleft,qright,fmean,zero_flux,params)
  use amr_parameters, only: ndim
  use hydro_parameters, only: nprim, nener
  use amr_commons, only: hydro_params_t
  use const
  IMPLICIT NONE

  type(hydro_params_t)::params
  real(kind=8),DIMENSION(1:nprim)::qleft,qright
  real(kind=8),DIMENSION(1:nprim+1)::fmean
  real(kind=8)::zero_flux

  real(kind=8),DIMENSION(1:nprim+1)::fleft,fright
  real(kind=8),DIMENSION(1:nprim+1)::uleft,uright,udiff

  real(kind=8), DIMENSION(7,7)::lem, rem
  real(kind=8), DIMENSION(7)  ::lambda, lambdal, lambdar, a

  real(kind=8)::droe, vxroe, vyroe, vzroe, sqrtdl, sqrtdr
  real(kind=8)::hroe, byroe, bzroe, Xfactor, Yfactor
  real(kind=8)::pbr, pbl

  real(kind=8)::fluxd, fluxmx, fluxmy, fluxmz, fluxe
  real(kind=8)::fluxby, fluxbz, byl, byr, bzl, bzr, bx
  real(kind=8)::dl, vxl, vyl, vzl, pl, hl
  real(kind=8)::dr, vxr, vyr, vzr, pr, hr
  real(kind=8)::mxl, myl, mzl, el
  real(kind=8)::mxr, myr, mzr, er
  real(kind=8)::coef, bx_mean
  real(kind=8)::vleft,vright
  real(kind=8)::dim,eim,mxm,mym,mzm,bym,bzm,etm,l1,l2

  INTEGER :: n
  LOGICAL :: llf

  ! Enforce continuity of normal component
  bx_mean=0.5d0*(qleft(6)+qright(6))
  qleft (6)=bx_mean
  qright(6)=bx_mean

  ! compute the fluxes and the conserved variables
  ! remember the convention : rho, E, rhovx, bx, rhovy, by, vz, rhobz
  CALL find_mhd_flux(qleft ,uleft ,fleft ,params)
  CALL find_mhd_flux(qright,uright,fright,params)

  ! define the primitive quantities explicitly
  dl  = qleft (1)
  dr  = qright(1)

  vxl = qleft (2)
  vxr = qright(2)

  vyl = qleft (3)
  vyr = qright(3)

  vzl = qleft (4)
  vzr = qright(4)

  pl  = qleft (5)
  pr  = qright(5)

  byl = qleft (7)
  byr = qright(7)

  bzl = qleft (8)
  bzr = qright(8)

  ! finally attribute bx value (left and right are identical)
  bx  = 0.5d0*(qleft(6)+qright(6))

  ! define explicitly the conserved quantities
  mxl = uleft (2)
  mxr = uright(2)

  myl = uleft (3)
  myr = uright(3)

  mzl = uleft (4)
  mzr = uright(4)

  el  = uleft (5)
  er  = uright(5)

  ! the magnetic pressure
  pbl = half * (bx*bx + byl*byl + bzl*bzl)
  pbr = half * (bx*bx + byr*byr + bzr*bzr)

  ! the total (specific) enthalpy
  hl = (el+pl+pbl)/dl
  hr = (er+pr+pbr)/dr
  !
  ! Step 2 : Compute Roe-averaged data from left and right states
  !
  sqrtdl  = sqrt(dl)
  sqrtdr  = sqrt(dr)
  droe    = sqrtdl*sqrtdr
  vxroe   = (sqrtdl*vxl + sqrtdr*vxr)/(sqrtdl+sqrtdr)
  vyroe   = (sqrtdl*vyl + sqrtdr*vyr)/(sqrtdl+sqrtdr)
  vzroe   = (sqrtdl*vzl + sqrtdr*vzr)/(sqrtdl+sqrtdr)
  byroe   = (sqrtdr*byl + sqrtdl*byr)/(sqrtdl+sqrtdr)
  bzroe   = (sqrtdr*bzl + sqrtdl*bzr)/(sqrtdl+sqrtdr)
  hroe    = (sqrtdl*hl  + sqrtdr*hr )/(sqrtdl+sqrtdr)
  Xfactor = ((byroe*byroe-byl*byr)+(bzroe*bzroe-bzl*bzr))/(2*droe)
  Yfactor = (dl+dr)/(2*droe)
  !
  ! Step 3 : Compute eigenvalues and eigenmatrices from Roe-averaged values
  !
  call eigen_cons(droe,vxroe,vyroe,vzroe,hroe,bx,byroe,bzroe,Xfactor,Yfactor,lambda,rem,lem,params)
  !
  ! Step 4: Compute eigenvalues from left and right states
  !
  call eigenvalues(dl,vxl,vyl,vzl,pl,bx,byl,bzl,lambdal,params)
  call eigenvalues(dr,vxr,vyr,vzr,pr,bx,byr,bzr,lambdar,params)
  !
  ! Step 5 : Create intermediate states from eigenmatrices
  !
  do n = 1, 7
     a(n)  = 0.0
     a(n)  = a(n)  + (dr -dl ) * lem(1,n)
     a(n)  = a(n)  + (mxr-mxl) * lem(2,n)
     a(n)  = a(n)  + (myr-myl) * lem(3,n)
     a(n)  = a(n)  + (mzr-mzl) * lem(4,n)
     a(n)  = a(n)  + (er - el) * lem(5,n)
     a(n)  = a(n)  + (byr-byl) * lem(6,n)
     a(n)  = a(n)  + (bzr-bzl) * lem(7,n)
  end do

  llf = .false.
  dim = dl
  mxm = mxl
  mym = myl
  mzm = mzl
  eim = el
  bym = byl
  bzm = bzl
  do n = 1, 7
     dim = dim + a(n) * rem(n,1)
     mxm = mxm + a(n) * rem(n,2)
     mym = mym + a(n) * rem(n,3)
     mzm = mzm + a(n) * rem(n,4)
     eim = eim + a(n) * rem(n,5)
     bym = bym + a(n) * rem(n,6)
     bzm = bzm + a(n) * rem(n,7)
     etm = eim-0.5*(mxm*mxm+mym*mym+mzm*mzm)/dim-0.5*(bx*bx+bym*bym+bzm*bzm)
     if(dim.le.zero.or.etm.le.zero)then
        llf=.true.
     endif
  end do

  IF( llf ) THEN
     fmean = half * ( fright + fleft ) * zero_flux
     CALL find_speed_info(qleft ,vleft ,params)
     CALL find_speed_info(qright,vright,params)
     udiff = half * ( uright - uleft )
     fmean = fmean - MAX(vleft,vright) * udiff
     RETURN
  END IF
  !
  ! Step 6 : Entropy fix for genuinely non linear waves
  !
  do n = 1, 7, 2
     l1 = min(lambdal(n),lambda(n))
     l2 = max(lambdar(n),lambda(n))
     if(l1.lt.zero.and.l2.gt.zero)then
        lambda(n)=(lambda(n)*(l2+l1)-two*l2*l1)/(l2-l1)
     endif
  end do
  !
  ! Step 6 : Compute fluxes at interface using  Roe  solver
  !
  ! add the left and right fluxes
  ! remember the convention : rho, E, rhovx, bx, rhovy, by, vz, rhobz
  fleft  = fleft  * zero_flux
  fright = fright * zero_flux
  fluxd  = fleft(1) + fright(1)
  fluxmx = fleft(2) + fright(2)
  fluxmy = fleft(3) + fright(3)
  fluxmz = fleft(4) + fright(4)
  fluxe  = fleft(5) + fright(5)
  fluxby = fleft(7) + fright(7)
  fluxbz = fleft(8) + fright(8)

  ! now compute the Roe fluxes
  ! remember that the convention of athena's eigenvalues is :
  ! rho,rhovx,rhovy,rhovz,E,by,bz
  DO n = 1, 7
      coef = ABS(lambda(n))*a(n)
      fluxd   = fluxd  - coef*rem(n,1)
      fluxmx  = fluxmx - coef*rem(n,2)
      fluxmy  = fluxmy - coef*rem(n,3)
      fluxmz  = fluxmz - coef*rem(n,4)
      fluxe   = fluxe  - coef*rem(n,5)
      fluxby  = fluxby - coef*rem(n,6)
      fluxbz  = fluxbz - coef*rem(n,7)
  ENDDO
  ! take half and put into the fmean variables
  fmean(1) = half * fluxd
  fmean(2) = half * fluxmx
  fmean(3) = half * fluxmy
  fmean(4) = half * fluxmz
  fmean(5) = half * fluxe
  fmean(6) = zero
  fmean(7) = half * fluxby
  fmean(8) = half * fluxbz
#if NVAR>5
  DO n = 9, nprim
     if(fmean(1)>0)then
        fmean(n)=qleft (n)*fmean(1)
     else
        fmean(n)=qright(n)*fmean(1)
     endif
  END DO
#endif

END SUBROUTINE riemann_roe_mhd
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine eigenvalues(d,vx,vy,vz,p,bx,by,bz,lambda,params)
  use amr_parameters, only: ndim
  use hydro_parameters, only: nprim, nener
  use amr_commons, only: hydro_params_t
  use const
  IMPLICIT NONE
!
! MHD adiabatic eigenvalues
!
! Input Arguments:
!   Bx    = magnetic field in sweep direction
!   d     = density
!   vx    = X velocity
!   vy    = Y velocity
!   vz    = Z velocity
!   by    = Y magnetic field
!   bz    = Z magnetic field
!   p     = thermal pressure
!
! Output Arguments:
!
!   lambda  = eigenvalues
!
  type(hydro_params_t)::params
  real(kind=8),intent(IN)::d, vx, vy, vz, p
  real(kind=8),intent(IN)::bx, by, bz
  real(kind=8),dimension(1:7),intent(OUT)::lambda
  ! local variables
  real(kind=8)::btsq, bt, vsq, vax,  vaxsq
  real(kind=8)::asq, astarsq, cfsq, cfast, cssq, cslow
  real(kind=8)::gamma, smallc

  gamma = params%gamma
  smallc = params%smallc

  vsq = vx**2+vy**2+vz**2
  btsq = by**2+bz**2
  bt = sqrt(btsq)
  vaxsq = Bx**2/d
  vax = sqrt(vaxsq)
  asq = gamma*p/d
  asq = MAX(asq,smallc**2)
  astarsq = asq+vaxsq+btsq/d

  cfsq = .5*(astarsq + sqrt(astarsq**2-4.0*asq*vaxsq))
  cfast = sqrt(cfsq)

  cssq = .5*(astarsq - sqrt(astarsq**2-4.0*asq*vaxsq))
  if (cssq .le. 0.) cssq = 0.
  cslow = sqrt(cssq)

  lambda(1) = vx - cfast
  lambda(2) = vx - vax
  lambda(3) = vx - cslow
  lambda(4) = vx
  lambda(5) = vx + cslow
  lambda(6) = vx + vax
  lambda(7) = vx + cfast

end subroutine eigenvalues
!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE eigen_cons(d,vx,vy,vz,h,Bx,by,bz,Xfac,Yfac,lambda,rem,lem,params)
  use amr_parameters, only: ndim
  use hydro_parameters, only: nprim, nener
  use amr_commons, only: hydro_params_t
  use const
  IMPLICIT NONE
!
! Input Arguments:
!   d     = Roe density
!   vx    = Roe X velocity
!   vy    = Roe Y velocity
!   vz    = Roe Z velocity
!   Bx    = magnetic field in sweep direction
!   by    = Roe Y magnetic field
!   bz    = Roe Z magnetic field
!   h     = Roe enthalpy
!   Xfac  = ((by^2-byl*byr)+bz^2-bzl*bzr))/(2*rho)
!   Yfac  = (rho_l+rho_r)/(2*rho)
!
! Output Arguments:
!
!   lambda = eigenvalues
!   rem     = right eigenmatrix
!   lem     = left  eigenmatrix
!
  type(hydro_params_t)::params
  real(kind=8)::d, vx, vy, vz, h
  real(kind=8)::Bx, by, bz, Xfac, Yfac

  real(kind=8),DIMENSION(7,7)::lem, rem
  real(kind=8),DIMENSION(7)  ::lambda

  real(kind=8)::gamma, smallc
  real(kind=8)::btsq, bt_starsq, bt, bt_star
  real(kind=8)::vsq, vax,  vaxsq, hp, twid_asq, q_starsq
  real(kind=8)::cfsq, cfast, cssq, cslow
  real(kind=8)::beta_y, beta_z, beta_ystar, beta_zstar, beta_starsq, vbeta
  real(kind=8)::alpha_f, alpha_s, droot, s, twid_a
  real(kind=8)::Qfast, Qslow, af_prime, as_prime, Afpbb, Aspbb, na
  real(kind=8)::cff, css, af, as, Afpb, Aspb, vqstr, norm
  real(kind=8)::Q_ystar, Q_zstar

  gamma = params%gamma
  smallc = params%smallc

  vsq = vx*vx+vy*vy+vz*vz
  btsq = by*by+bz*bz
  bt_starsq = (gamma-1. - (gamma-2.)*Yfac)*btsq
  bt=sqrt(btsq)
  bt_star = sqrt(bt_starsq)

  vaxsq = Bx*Bx/d
  vax = sqrt(vaxsq)

  hp = h - (vaxsq + btsq/d)
  twid_asq = ((gamma-1.)*(hp-.5*vsq)-(gamma-2.)*Xfac)
  twid_asq = MAX(twid_asq,smallc**2)
  q_starsq = twid_asq+(vaxsq+bt_starsq/d)

  cfsq = .5*(q_starsq + sqrt(q_starsq*q_starsq-4.0*twid_asq*vaxsq))
  cfast = sqrt(cfsq)

  cssq = .5*(q_starsq - sqrt(q_starsq*q_starsq-4.0*twid_asq*vaxsq))
  if (cssq .le. 0.) cssq = 0.
  cslow = sqrt(cssq)

  if (bt .eq. 0) then
     beta_y = .5*sqrt(2.)
     beta_z = .5*sqrt(2.)
     beta_ystar = .5*sqrt(2.)
     beta_zstar = .5*sqrt(2.)
  else
     beta_y = by/bt
     beta_z = bz/bt
     beta_ystar = by/bt_star
     beta_zstar = bz/bt_star
  endif
  beta_starsq = beta_ystar*beta_ystar + beta_zstar*beta_zstar
  vbeta = vy*beta_ystar + vz*beta_zstar

  if ( (cfsq - cssq) .eq. 0.) then
     alpha_f = 1.0
     alpha_s = 0.0
  else if ( (twid_asq - cssq) .le. 0.) then
     alpha_f = 0.0
     alpha_s = 1.0
  else if ( (cfsq - twid_asq) .le. 0.) then
     alpha_f = 1.0
     alpha_s = 0.0
  else
     alpha_f = sqrt((twid_asq-cssq)/(cfsq-cssq))
     alpha_s = sqrt((cfsq-twid_asq)/(cfsq-cssq))
  endif
  !
  ! compute Qs and As for eigenmatrices
  !
  droot = sqrt(d)
  s  = SIGN(one,bx)
  twid_a = sqrt(twid_asq)
  Qfast = s*cfast*alpha_f
  Qslow = s*cslow*alpha_s
  af_prime = twid_a*alpha_f/droot
  as_prime = twid_a*alpha_s/droot
  Afpbb = af_prime*bt_star*beta_starsq
  Aspbb = as_prime*bt_star*beta_starsq
  !
  ! eigenvalues
  !
  lambda(1) = vx - cfast
  lambda(2) = vx - vax
  lambda(3) = vx - cslow
  lambda(4) = vx
  lambda(5) = vx + cslow
  lambda(6) = vx + vax
  lambda(7) = vx + cfast
  !
  ! eigenmatrix
  !
  rem(1,1) = alpha_f
  rem(1,2) = alpha_f*(vx-cfast)
  rem(1,3) = alpha_f*vy + Qslow*beta_ystar
  rem(1,4) = alpha_f*vz + Qslow*beta_zstar
  rem(1,5) = alpha_f*(hp-vx*cfast) + Qslow*vbeta + Aspbb
  rem(1,6) = as_prime*beta_ystar
  rem(1,7) = as_prime*beta_zstar

  rem(2,1) = 0.
  rem(2,2) = 0.
  rem(2,3) = -beta_z
  rem(2,4) =  beta_y
  rem(2,5) = -(vy*beta_z - vz*beta_y)
  rem(2,6) = -s*beta_z/droot
  rem(2,7) =  s*beta_y/droot

  rem(3,1) = alpha_s
  rem(3,2) = alpha_s*(vx-cslow)
  rem(3,3) = alpha_s*vy - Qfast*beta_ystar
  rem(3,4) = alpha_s*vz - Qfast*beta_zstar
  rem(3,5) = alpha_s*(hp - vx*cslow) - Qfast*vbeta - Afpbb
  rem(3,6) = -af_prime*beta_ystar
  rem(3,7) = -af_prime*beta_zstar

  rem(4,1) = 1.0
  rem(4,2) = vx
  rem(4,3) = vy
  rem(4,4) = vz
  rem(4,5) = 0.5*vsq + (gamma-2.)*Xfac/(gamma-1.)
  rem(4,6) = 0.
  rem(4,7) = 0.

  rem(5,1) =  alpha_s
  rem(5,2) =  alpha_s*(vx+cslow)
  rem(5,3) =  alpha_s*vy + Qfast*beta_ystar
  rem(5,4) =  alpha_s*vz + Qfast*beta_zstar
  rem(5,5) =  alpha_s*(hp+vx*cslow) + Qfast*vbeta - Afpbb
  rem(5,6) =  rem(3,6)
  rem(5,7) =  rem(3,7)

  rem(6,1) = 0.
  rem(6,2) = 0.
  rem(6,3) =  beta_z
  rem(6,4) = -beta_y
  rem(6,5) = -rem(2,5)
  rem(6,6) =  rem(2,6)
  rem(6,7) =  rem(2,7)

  rem(7,1) =  alpha_f
  rem(7,2) = alpha_f*(vx+cfast)
  rem(7,3) = alpha_f*vy - Qslow*beta_ystar
  rem(7,4) = alpha_f*vz - Qslow*beta_zstar
  rem(7,5) = alpha_f*(hp + vx*cfast) - Qslow*vbeta + Aspbb
  rem(7,6) =  rem(1,6)
  rem(7,7) =  rem(1,7)
  !
  ! Left eignematrix
  !
  ! normalize some of the quantities by 1/(2a^2)
  ! some by (gamma-1)/2a^2
  !
  na  = 0.5/twid_asq
  cff = na*alpha_f*cfast
  css = na*alpha_s*cslow
  Qfast = Qfast*na
  Qslow = Qslow*na
  af = na*af_prime*d
  as = na*as_prime*d
  Afpb = na*af_prime*bt_star
  Aspb = na*as_prime*bt_star
  !
  alpha_f = (gamma-1.)*na*alpha_f
  alpha_s = (gamma-1.)*na*alpha_s
  Q_ystar = beta_ystar/beta_starsq
  Q_zstar = beta_zstar/beta_starsq
  vqstr   = (vy*Q_ystar+vz*Q_zstar)
  norm = (gamma-1.)*2.*na

  lem(1,1) = alpha_f*(vsq-hp) + cff*(cfast+vx) - Qslow*vqstr - Aspb
  lem(2,1) = -alpha_f*vx - cff
  lem(3,1) = -alpha_f*vy + Qslow*Q_ystar
  lem(4,1) = -alpha_f*vz + Qslow*Q_zstar
  lem(5,1) =  alpha_f
  lem(6,1) = as*Q_ystar - alpha_f*by
  lem(7,1) = as*Q_zstar - alpha_f*bz

  lem(1,2) =  0.5*(vy*beta_z - vz*beta_y)
  lem(2,2) =  0.
  lem(3,2) = -0.5*beta_z
  lem(4,2) =  0.5*beta_y
  lem(5,2) =  0.
  lem(6,2) = -0.5*droot*beta_z*s
  lem(7,2) =  0.5*droot*beta_y*s

  lem(1,3) =  alpha_s*(vsq-hp) + css*(cslow+vx) + Qfast*vqstr + Afpb
  lem(2,3) = -alpha_s*vx - css
  lem(3,3) = -alpha_s*vy - Qfast*Q_ystar
  lem(4,3) = -alpha_s*vz - Qfast*Q_zstar
  lem(5,3) =  alpha_s
  lem(6,3) = -af*Q_ystar - alpha_s*by
  lem(7,3) = -af*Q_zstar - alpha_s*bz

  ! Attention, il y a une difference de signe avec ramses_mhd...
  lem(1,4) =  1. - norm*(.5*vsq - (gamma-2.)*Xfac/(gamma-1.))
  ! This is the old version...
  !  lem(1,4) =  1. - norm*(.5*vsq + (gamma-2.)*Xfac/(gamma-1.))
  lem(2,4) =  norm*vx
  lem(3,4) =  norm*vy
  lem(4,4) =  norm*vz
  lem(5,4) = -norm
  lem(6,4) =  norm*by
  lem(7,4) =  norm*bz

  lem(1,5) =  alpha_s*(vsq-hp) + css*(cslow-vx) - Qfast*vqstr + Afpb
  lem(2,5) = -alpha_s*vx + css
  lem(3,5) = -alpha_s*vy + Qfast*Q_ystar
  lem(4,5) = -alpha_s*vz + Qfast*Q_zstar
  lem(5,5) =  alpha_s
  lem(6,5) =  lem(6,3)
  lem(7,5) =  lem(7,3)

  lem(1,6) = -lem(1,2)
  lem(2,6) =  0.
  lem(3,6) = -lem(3,2)
  lem(4,6) = -lem(4,2)
  lem(5,6) =  0.
  lem(6,6) =  lem(6,2)
  lem(7,6) =  lem(7,2)

  lem(1,7) =  alpha_f*(vsq-hp) + cff*(cfast-vx) + Qslow*vqstr -Aspb
  lem(2,7) = -alpha_f*vx + cff
  lem(3,7) = -alpha_f*vy - Qslow*Q_ystar
  lem(4,7) = -alpha_f*vz - Qslow*Q_zstar
  lem(5,7) =  alpha_f
  lem(6,7) =  lem(6,1)
  lem(7,7) =  lem(7,1)

END SUBROUTINE eigen_cons
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine riemann2d_hlld(qLL,qLR,qRL,qRR,E,params)
  use amr_parameters, only: ndim
  use hydro_parameters, only: nprim, nener
  use amr_commons, only: hydro_params_t
  use const
  ! HLLD 2D Riemann solver (Miyoshi & Kusano, 2005, JCP, 208, 315)
  IMPLICIT NONE
  type(hydro_params_t)::params
  real(kind=8),DIMENSION(1:nprim)::qLL,qLR,qRL,qRR
  real(kind=8)::E

  real(kind=8),DIMENSION(1:nprim)::qtmp
  real(kind=8)::ELL,ERL,ELR,ERR,SL,SR,SB,ST,SAL,SAR,SAT,SAB
  real(kind=8)::cfastLLx,cfastRLx,cfastLRx,cfastRRx,cfastLLy,cfastRLy,cfastLRy,cfastRRy
  real(kind=8)::calfvenR,calfvenL,calfvenT,calfvenB
  real(kind=8)::rLL,rLR,rRL,rRR,pLL,pLR,pRL,pRR,uLL,uLR,uRL,uRR,vLL,vLR,vRL,vRR
  real(kind=8)::ALL,ALR,ARL,ARR,BLL,BLR,BRL,BRR,CLL,CLR,CRL,CRR
  real(kind=8)::PtotLL,PtotLR,PtotRL,PtotRR,rcLLx,rcLRx,rcRLx,rcRRx,rcLLy,rcLRy,rcRLy,rcRRy
  real(kind=8)::ustar,vstar,rstarLLx,rstarLRx,rstarRLx,rstarRRx,rstarLLy,rstarLRy,rstarRLy,rstarRRy
  real(kind=8)::rstarLL,rstarLR,rstarRL,rstarRR,AstarLL,AstarLR,AstarRL,AstarRR,BstarLL,BstarLR,BstarRL,BstarRR
  real(kind=8)::EstarLLx,EstarLRx,EstarRLx,EstarRRx
  real(kind=8)::EstarLLy,EstarLRy,EstarRLy,EstarRRy
  real(kind=8)::EstarLL,EstarLR,EstarRL,EstarRR
  real(kind=8)::AstarT,AstarB,BstarR,BstarL
  real(kind=8)::rmin,pmin,Smax
  real(kind=8)::smallc,switch_llf_dmin,switch_llf_pmin
#if NENER>0
  INTEGER::irad
#endif
  logical::switch_to_llf

  switch_llf_dmin = params%switch_llf_dmin
  switch_llf_pmin = params%switch_llf_pmin
  smallc = params%smallc

  rLL=qLL(1); pLL=qLL(5); uLL=qLL(2); vLL=qLL(3); ALL=qLL(6); BLL=qLL(7) ; CLL=qLL(8)
  rLR=qLR(1); pLR=qLR(5); uLR=qLR(2); vLR=qLR(3); ALR=qLR(6); BLR=qLR(7) ; CLR=qLR(8)
  rRL=qRL(1); pRL=qRL(5); uRL=qRL(2); vRL=qRL(3); ARL=qRL(6); BRL=qRL(7) ; CRL=qRL(8)
  rRR=qRR(1); pRR=qRR(5); uRR=qRR(2); vRR=qRR(3); ARR=qRR(6); BRR=qRR(7) ; CRR=qRR(8)
  rmin=MIN(rLL, rLR, rRL, rRR)
  pmin=MIN(pLL, pLR, pRL, pRR)

  ! Compute 4 fast magnetosonic velocity relative to x direction
  qtmp(1)=rLL; qtmp(5)=pLL; qtmp(8)=CLL
  qtmp(6)=ALL; qtmp(7)=BLL
#if NENER>0
  do irad = 1,nener
     qtmp(8+irad) = qLL(8+irad)
  end do
#endif
  call find_speed_fast(qtmp,cfastLLx,params)
  qtmp(1)=rLR; qtmp(5)=pLR; qtmp(8)=CLR
  qtmp(6)=ALR; qtmp(7)=BLR
#if NENER>0
  do irad = 1,nener
     qtmp(8+irad) = qLR(8+irad)
  end do
#endif
  call find_speed_fast(qtmp,cfastLRx,params)
  qtmp(1)=rRL; qtmp(5)=pRL; qtmp(8)=CRL
  qtmp(6)=ARL; qtmp(7)=BRL
#if NENER>0
  do irad = 1,nener
     qtmp(8+irad) = qRL(8+irad)
  end do
#endif
  call find_speed_fast(qtmp,cfastRLx,params)
  qtmp(1)=rRR; qtmp(5)=pRR; qtmp(8)=CRR
  qtmp(6)=ARR; qtmp(7)=BRR
#if NENER>0
  do irad = 1,nener
     qtmp(8+irad) = qRR(8+irad)
  end do
#endif
  call find_speed_fast(qtmp,cfastRRx,params)

  ! Compute 4 fast magnetosonic velocity relative to y direction
  qtmp(1)=rLL; qtmp(5)=pLL; qtmp(8)=CLL
  qtmp(6)=BLL; qtmp(7)=ALL
#if NENER>0
  do irad = 1,nener
     qtmp(8+irad) = qLL(8+irad)
  end do
#endif
  call find_speed_fast(qtmp,cfastLLy,params)
  qtmp(1)=rLR; qtmp(5)=pLR; qtmp(8)=CLR
  qtmp(6)=BLR; qtmp(7)=ALR
#if NENER>0
  do irad = 1,nener
     qtmp(8+irad) = qLR(8+irad)
  end do
#endif
  call find_speed_fast(qtmp,cfastLRy,params)
  qtmp(1)=rRL; qtmp(5)=pRL; qtmp(8)=CRL
  qtmp(6)=BRL; qtmp(7)=ARL
#if NENER>0
  do irad = 1,nener
     qtmp(8+irad) = qRL(8+irad)
  end do
#endif
  call find_speed_fast(qtmp,cfastRLy,params)
  qtmp(1)=rRR; qtmp(5)=pRR; qtmp(8)=CRR
  qtmp(6)=BRR; qtmp(7)=ARR
#if NENER>0
  do irad = 1,nener
     qtmp(8+irad) = qRR(8+irad)
  end do
#endif
  call find_speed_fast(qtmp,cfastRRy,params)

  SL=min(uLL,uLR,uRL,uRR)-max(cfastLLx,cfastLRx,cfastRLx,cfastRRx)
  SR=max(uLL,uLR,uRL,uRR)+max(cfastLLx,cfastLRx,cfastRLx,cfastRRx)
  SB=min(vLL,vLR,vRL,vRR)-max(cfastLLy,cfastLRy,cfastRLy,cfastRRy)
  ST=max(vLL,vLR,vRL,vRR)+max(cfastLLy,cfastLRy,cfastRLy,cfastRRy)
  Smax = max(abs(SR), abs(ST), abs(SL), abs(SB))

  ELL=uLL*BLL-vLL*ALL
  ELR=uLR*BLR-vLR*ALR
  ERL=uRL*BRL-vRL*ARL
  ERR=uRR*BRR-vRR*ARR

  ! Switch to llf and exit
  switch_to_llf=.false.
  if(switch_llf_dmin.gt.0)switch_to_llf=rmin.lt.switch_llf_dmin
  if(switch_llf_pmin.gt.0.and.pmin.lt.switch_llf_pmin)switch_to_llf=.true.
  if(switch_to_llf)then
     E = forth*(ERR+ERL+ELR+ELL)+half*Smax*(qRR(6)-qLL(6))-half*Smax*(qRR(7)-qLL(7))
     return
  endif

#if NENER>0
 do irad = 1,nener
     pLL = pLL + qLL(8+irad)
     pLR = pLR + qLR(8+irad)
     pRL = pRL + qRL(8+irad)
     pRR = pRR + qRR(8+irad)
  end do
#endif
  PtotLL=pLL+half*(ALL*ALL+BLL*BLL+CLL*CLL)
  PtotLR=pLR+half*(ALR*ALR+BLR*BLR+CLR*CLR)
  PtotRL=pRL+half*(ARL*ARL+BRL*BRL+CRL*CRL)
  PtotRR=pRR+half*(ARR*ARR+BRR*BRR+CRR*CRR)

  rcLLx=rLL*(uLL-SL); rcRLx=rRL*(SR-uRL)
  rcLRx=rLR*(uLR-SL); rcRRx=rRR*(SR-uRR)
  rcLLy=rLL*(vLL-SB); rcLRy=rLR*(ST-vLR)
  rcRLy=rRL*(vRL-SB); rcRRy=rRR*(ST-vRR)

  ustar=(rcLLx*uLL+rcLRx*uLR+rcRLx*uRL+rcRRx*uRR+(PtotLL-PtotRL+PtotLR-PtotRR))/(rcLLx+rcLRx+rcRLx+rcRRx)
  vstar=(rcLLy*vLL+rcLRy*vLR+rcRLy*vRL+rcRRy*vRR+(PtotLL-PtotLR+PtotRL-PtotRR))/(rcLLy+rcLRy+rcRLy+rcRRy)

  rstarLLx=rLL*(SL-uLL)/(SL-ustar); BstarLL=BLL*(SL-uLL)/(SL-ustar)
  rstarLLy=rLL*(SB-vLL)/(SB-vstar); AstarLL=ALL*(SB-vLL)/(SB-vstar)
  rstarLL =rLL*(SL-uLL)/(SL-ustar)*(SB-vLL)/(SB-vstar)
  EstarLLx=ustar*BstarLL-vLL  *ALL
  EstarLLy=uLL  *BLL    -vstar*AstarLL
  EstarLL =ustar*BstarLL-vstar*AstarLL

  rstarLRx=rLR*(SL-uLR)/(SL-ustar); BstarLR=BLR*(SL-uLR)/(SL-ustar)
  rstarLRy=rLR*(ST-vLR)/(ST-vstar); AstarLR=ALR*(ST-vLR)/(ST-vstar)
  rstarLR =rLR*(SL-uLR)/(SL-ustar)*(ST-vLR)/(ST-vstar)
  EstarLRx=ustar*BstarLR-vLR  *ALR
  EstarLRy=uLR  *BLR    -vstar*AstarLR
  EstarLR =ustar*BstarLR-vstar*AstarLR

  rstarRLx=rRL*(SR-uRL)/(SR-ustar); BstarRL=BRL*(SR-uRL)/(SR-ustar)
  rstarRLy=rRL*(SB-vRL)/(SB-vstar); AstarRL=ARL*(SB-vRL)/(SB-vstar)
  rstarRL =rRL*(SR-uRL)/(SR-ustar)*(SB-vRL)/(SB-vstar)
  EstarRLx=ustar*BstarRL-vRL  *ARL
  EstarRLy=uRL  *BRL    -vstar*AstarRL
  EstarRL =ustar*BstarRL-vstar*AstarRL

  rstarRRx=rRR*(SR-uRR)/(SR-ustar); BstarRR=BRR*(SR-uRR)/(SR-ustar)
  rstarRRy=rRR*(ST-vRR)/(ST-vstar); AstarRR=ARR*(ST-vRR)/(ST-vstar)
  rstarRR =rRR*(SR-uRR)/(SR-ustar)*(ST-vRR)/(ST-vstar)
  EstarRRx=ustar*BstarRR-vRR  *ARR
  EstarRRy=uRR  *BRR    -vstar*AstarRR
  EstarRR =ustar*BstarRR-vstar*AstarRR

  calfvenL=max(abs(ALR)/sqrt(rstarLRx),abs(AstarLR)/sqrt(rstarLR), &
       &       abs(ALL)/sqrt(rstarLLx),abs(AstarLL)/sqrt(rstarLL),smallc)
  calfvenR=max(abs(ARR)/sqrt(rstarRRx),abs(AstarRR)/sqrt(rstarRR), &
       &       abs(ARL)/sqrt(rstarRLx),abs(AstarRL)/sqrt(rstarRL),smallc)
  calfvenB=max(abs(BLL)/sqrt(rstarLLy),abs(BstarLL)/sqrt(rstarLL), &
       &       abs(BRL)/sqrt(rstarRLy),abs(BstarRL)/sqrt(rstarRL),smallc)
  calfvenT=max(abs(BLR)/sqrt(rstarLRy),abs(BstarLR)/sqrt(rstarLR), &
       &       abs(BRR)/sqrt(rstarRRy),abs(BstarRR)/sqrt(rstarRR),smallc)
  SAL=min(ustar-calfvenL,zero); SAR=max(ustar+calfvenR,zero)
  SAB=min(vstar-calfvenB,zero); SAT=max(vstar+calfvenT,zero)
  AstarT=(SAR*AstarRR-SAL*AstarLR)/(SAR-SAL); AstarB=(SAR*AstarRL-SAL*AstarLL)/(SAR-SAL)
  BstarR=(SAT*BstarRR-SAB*BstarRL)/(SAT-SAB); BstarL=(SAT*BstarLR-SAB*BstarLL)/(SAT-SAB)

  if(SB>0d0)then
     if(SL>0d0)then
        E=ELL
     else if(SR<0d0)then
        E=ERL
     else
        E=(SAR*EstarLLx-SAL*EstarRLx+SAR*SAL*(BRL-BLL))/(SAR-SAL)
     endif
  else if (ST<0d0)then
     if(SL>0d0)then
        E=ELR
     else if(SR<0d0)then
        E=ERR
     else
        E=(SAR*EstarLRx-SAL*EstarRRx+SAR*SAL*(BRR-BLR))/(SAR-SAL)
     endif
  else if(SL>0d0)then
     E=(SAT*EstarLLy-SAB*EstarLRy-SAT*SAB*(ALR-ALL))/(SAT-SAB)
  else if(SR<0d0)then
     E=(SAT*EstarRLy-SAB*EstarRRy-SAT*SAB*(ARR-ARL))/(SAT-SAB)
  else
     E=(SAL*SAB*EstarRR-SAL*SAT*EstarRL-SAR*SAB*EstarLR+SAR*SAT*EstarLL)/(SAR-SAL)/(SAT-SAB) &
          & -SAT*SAB/(SAT-SAB)*(AstarT-AstarB)+SAR*SAL/(SAR-SAL)*(BstarR-BstarL)
  endif

end subroutine riemann2d_hlld
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine riemann2d_hll_fast(qLL,qLR,qRL,qRR,E,params)
  use amr_parameters, only: ndim
  use hydro_parameters, only: nprim, nener
  use amr_commons, only: hydro_params_t
  use const
  ! HLL 2D Riemann solver using the fast magnetosonic speed
  IMPLICIT NONE
  type(hydro_params_t)::params
  real(kind=8),DIMENSION(1:nprim)::qLL,qLR,qRL,qRR
  real(kind=8)::E

  real(kind=8),DIMENSION(1:nprim)::qtmp
  real(kind=8)::ELL,ERL,ELR,ERR,SL,SR,SB,ST
  real(kind=8)::cLLx,cRLx,cLRx,cRRx,cLLy,cRLy,cLRy,cRRy
  real(kind=8)::uLL,uRL,uLR,uRR,vLL,vRL,vLR,vRR
#if NENER>0
  INTEGER::irad
#endif

  ! vx*by - vy*bx at the four edge centers
  ELL = qLL(2)*qLL(7) - qLL(3)*qLL(6)
  ERL = qRL(2)*qRL(7) - qRL(3)*qRL(6)
  ELR = qLR(2)*qLR(7) - qLR(3)*qLR(6)
  ERR = qRR(2)*qRR(7) - qRR(3)*qRR(6)

  ! Compute 4 fast magnetosonic velocity relative to x direction
  qtmp(1)=qLL(1); qtmp(5)=qLL(5); qtmp(8)=qLL(8)
  qtmp(6)=qLL(6); qtmp(7)=qLL(7)
#if NENER>0
  do irad = 1,nener
     qtmp(8+irad) = qLL(8+irad)
  end do
#endif
  uLL=qtmp(2); call find_speed_fast(qtmp,cLLx,params)
  qtmp(1)=qLR(1); qtmp(5)=qLR(5); qtmp(8)=qLR(8)
  qtmp(6)=qLR(6); qtmp(7)=qLR(7)
#if NENER>0
  do irad = 1,nener
     qtmp(8+irad) = qLR(8+irad)
  end do
#endif
  uLR=qtmp(2); call find_speed_fast(qtmp,cLRx,params)
  qtmp(1)=qRL(1); qtmp(5)=qRL(5); qtmp(8)=qRL(8)
  qtmp(6)=qRL(6); qtmp(7)=qRL(7)
#if NENER>0
  do irad = 1,nener
     qtmp(8+irad) = qRL(8+irad)
  end do
#endif
  uRL=qtmp(2); call find_speed_fast(qtmp,cRLx,params)
  qtmp(1)=qRR(1); qtmp(5)=qRR(5); qtmp(8)=qRR(8)
  qtmp(6)=qRR(6); qtmp(7)=qRR(7)
#if NENER>0
  do irad = 1,nener
     qtmp(8+irad) = qRR(8+irad)
  end do
#endif
  uRR=qtmp(2); call find_speed_fast(qtmp,cRRx,params)

  ! Compute 4 fast magnetosonic velocity relative to y direction
  qtmp(1)=qLL(1); qtmp(5)=qLL(5); qtmp(8)=qLL(8)
  qtmp(6)=qLL(7); qtmp(7)=qLL(6)
#if NENER>0
  do irad = 1,nener
     qtmp(8+irad) = qLL(8+irad)
  end do
#endif
  vLL=qtmp(3); call find_speed_fast(qtmp,cLLy,params)
  qtmp(1)=qLR(1); qtmp(5)=qLR(5); qtmp(8)=qLR(8)
  qtmp(6)=qLR(7); qtmp(7)=qLR(6)
#if NENER>0
  do irad = 1,nener
     qtmp(8+irad) = qLR(8+irad)
  end do
#endif
  vLR=qtmp(3); call find_speed_fast(qtmp,cLRy,params)
  qtmp(1)=qRL(1); qtmp(5)=qRL(5); qtmp(8)=qRL(8)
  qtmp(6)=qRL(7); qtmp(7)=qRL(6)
#if NENER>0
  do irad = 1,nener
     qtmp(8+irad) = qRL(8+irad)
  end do
#endif
  vRL=qtmp(3); call find_speed_fast(qtmp,cRLy,params)
  qtmp(1)=qRR(1); qtmp(5)=qRR(5); qtmp(8)=qRR(8)
  qtmp(6)=qRR(7); qtmp(7)=qRR(6)
#if NENER>0
  do irad = 1,nener
     qtmp(8+irad) = qRR(8+irad)
  end do
#endif
  vRR=qtmp(3); call find_speed_fast(qtmp,cRRy,params)

  SL=min(min(uLL,uLR,uRL,uRR)-max(cLLx,cLRx,cRLx,cRRx),zero)
  SR=max(max(uLL,uLR,uRL,uRR)+max(cLLx,cLRx,cRLx,cRRx),zero)
  SB=min(min(vLL,vLR,vRL,vRR)-max(cLLy,cLRy,cRLy,cRRy),zero)
  ST=max(max(vLL,vLR,vRL,vRR)+max(cLLy,cLRy,cRLy,cRRy),zero)

  E = (SL*SB*ERR-SL*ST*ERL-SR*SB*ELR+SR*ST*ELL)/(SR-SL)/(ST-SB) &
       & -ST*SB/(ST-SB)*(qRR(6)-qLL(6)) &
       & +SR*SL/(SR-SL)*(qRR(7)-qLL(7))

end subroutine riemann2d_hll_fast
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine riemann2d_hll_alfven(qLL,qLR,qRL,qRR,E)
  use amr_parameters, only: ndim
  use hydro_parameters, only: nprim, nener
  use const
  ! HLL 2D Riemann solver using the Alfven speed
  IMPLICIT NONE
  real(kind=8),DIMENSION(1:nprim)::qLL,qLR,qRL,qRR
  real(kind=8)::E

  real(kind=8),DIMENSION(1:nprim)::qtmp
  real(kind=8)::ELL,ERL,ELR,ERR,SL,SR,SB,ST
  real(kind=8)::cLLx,cRLx,cLRx,cRRx,cLLy,cRLy,cLRy,cRRy
  real(kind=8)::uLL,uRL,uLR,uRR,vLL,vRL,vLR,vRR

  ! vx*by - vy*bx at the four edge centers
  ELL = qLL(2)*qLL(7) - qLL(3)*qLL(6)
  ERL = qRL(2)*qRL(7) - qRL(3)*qRL(6)
  ELR = qLR(2)*qLR(7) - qLR(3)*qLR(6)
  ERR = qRR(2)*qRR(7) - qRR(3)*qRR(6)

  ! Compute 4 Alfven velocity relative to x direction
  qtmp(1)=qLL(1); qtmp(6)=qLL(6)
  uLL=qtmp(2); call find_speed_alfven(qtmp,cLLx)
  qtmp(1)=qLR(1); qtmp(6)=qLR(6)
  uLR=qtmp(2); call find_speed_alfven(qtmp,cLRx)
  qtmp(1)=qRL(1); qtmp(6)=qRL(6)
  uRL=qtmp(2); call find_speed_alfven(qtmp,cRLx)
  qtmp(1)=qRR(1); qtmp(6)=qRR(6)
  uRR=qtmp(2); call find_speed_alfven(qtmp,cRRx)

  ! Compute 4 Alfven relative to y direction
  qtmp(1)=qLL(1); qtmp(6)=qLL(7)
  vLL=qtmp(3); call find_speed_alfven(qtmp,cLLy)
  qtmp(1)=qLR(1); qtmp(6)=qLR(7)
  vLR=qtmp(3); call find_speed_alfven(qtmp,cLRy)
  qtmp(1)=qRL(1); qtmp(6)=qRL(7)
  vRL=qtmp(3); call find_speed_alfven(qtmp,cRLy)
  qtmp(1)=qRR(1); qtmp(6)=qRR(7)
  vRR=qtmp(3); call find_speed_alfven(qtmp,cRRy)

  SL=min(min(uLL,uLR,uRL,uRR)-max(cLLx,cLRx,cRLx,cRRx),zero)
  SR=max(max(uLL,uLR,uRL,uRR)+max(cLLx,cLRx,cRLx,cRRx),zero)
  SB=min(min(vLL,vLR,VRL,vRR)-max(cLLy,cLRy,cRLy,cRRy),zero)
  ST=max(max(vLL,vLR,VRL,vRR)+max(cLLy,cLRy,cRLy,cRRy),zero)

  E = (SL*SB*ERR-SL*ST*ERL-SR*SB*ELR+SR*ST*ELL)/(SR-SL)/(ST-SB) &
       & -ST*SB/(ST-SB)*(qRR(6)-qLL(6)) &
       & +SR*SL/(SR-SL)*(qRR(7)-qLL(7))

end subroutine riemann2d_hll_alfven
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine riemann2d_simple(qLL,qLR,qRL,qRR,E,params)
  use amr_parameters, only: ndim
  use amr_commons, only: hydro_params_t
  use hydro_parameters
  use const
  ! HLL 2D Riemann solver using the Alfven speed
  IMPLICIT NONE
  type(hydro_params_t)::params
  real(kind=8),DIMENSION(1:nprim)::qLL,qLR,qRL,qRR
  real(kind=8)::E

  integer::riemann2d
  real(kind=8),DIMENSION(1:nprim)::qleft,qright
  real(kind=8),DIMENSION(1:nprim+1)::fmean_x,fmean_y
  real(kind=8)::ELL,ERL,ELR,ERR,zero_flux
#if NENER>0
  INTEGER::irad
#endif

  riemann2d=params%riemann2d

  ! vx*by - vy*bx at the four edge centers
  ELL = qLL(2)*qLL(7) - qLL(3)*qLL(6)
  ERL = qRL(2)*qRL(7) - qRL(3)*qRL(6)
  ELR = qLR(2)*qLR(7) - qLR(3)*qLR(6)
  ERR = qRR(2)*qRR(7) - qRR(3)*qRR(6)

  ! find the average value of E
  E = forth*(ELL+ERL+ELR+ERR)

  ! call the first solver in the x direction
  ! density
  qleft (1) = half*(qLL(1)+qLR(1))
  qright(1) = half*(qRR(1)+qRL(1))

  ! vt1 becomes normal velocity
  qleft (2) = half*(qLL(2)+qLR(2))
  qright(2) = half*(qRR(2)+qRL(2))

  ! vt2 becomes transverse velocity field
  qleft (3) = half*(qLL(3)+qLR(3))
  qright(3) = half*(qRR(3)+qRL(3))

  ! velocity component perp. to the plane is now transverse
  qleft (4) = half*(qLL(4)+qLR(4))
  qright(4) = half*(qRR(4)+qRL(4))

  ! pressure
  qleft (5) = half*(qLL(5)+qLR(5))
  qright(5) = half*(qRR(5)+qRL(5))

  ! bt1 becomes normal magnetic field
  qleft (6) = half*(qLL(6)+qLR(6))
  qright(6) = half*(qRR(6)+qRL(6))

  ! bt2 becomes transverse magnetic field
  qleft (7) = half*(qLL(7)+qLR(7))
  qright(7) = half*(qRR(7)+qRL(7))

  ! magnetic field component perp. to the plane is now transverse
  qleft (8) = half*(qLL(8)+qLR(8))
  qright(8) = half*(qRR(8)+qRL(8))

  ! non-thermal energies
#if NENER>0
  do irad = 1,nener
     qleft (8+irad) = half*(qLL(8+irad)+qLR(8+irad))
     qright(8+irad) = half*(qRR(8+irad)+qRL(8+irad))
  end do
#endif

  zero_flux = 0.
  SELECT CASE (riemann2d)
  CASE (solver2d_roe)
     CALL riemann_roe_mhd(qleft,qright,fmean_x,zero_flux,params)
  CASE (solver2d_llf)
     CALL riemann_llf_mhd(qleft,qright,fmean_x,zero_flux,params)
  CASE (solver2d_upwind)
     CALL riemann_upwind_mhd(qleft,qright,fmean_x,zero_flux,params)
  CASE DEFAULT
     write(*,*)'unknown 2D riemann solver'
     stop
  END SELECT

  ! call the second solver in the y direction
  ! density
  qleft (1) = half*(qLL(1)+qRL(1))
  qright(1) = half*(qRR(1)+qLR(1))

  ! vt2 becomes normal velocity
  qleft (2) = half*(qLL(3)+qRL(3))
  qright(2) = half*(qRR(3)+qLR(3))

  ! vt1 becomes transverse velocity field
  qleft (3) = half*(qLL(2)+qRL(2))
  qright(3) = half*(qRR(2)+qLR(2))

  ! velocity component perp. to the plane is now transverse
  qleft (4) = half*(qLL(4)+qRL(4))
  qright(4) = half*(qRR(4)+qLR(4))

  ! pressure
  qleft (5) = half*(qLL(5)+qRL(5))
  qright(5) = half*(qRR(5)+qLR(5))

  ! bt2 becomes normal magnetic field
  qleft (6) = half*(qLL(7)+qRL(7))
  qright(6) = half*(qRR(7)+qLR(7))

  ! bt1 becomes transverse magnetic field
  qleft (7) = half*(qLL(6)+qRL(6))
  qright(7) = half*(qRR(6)+qLR(6))

  ! magnetic field component perp. to the plane is now transverse
  qleft (8) = half*(qLL(8)+qRL(8))
  qright(8) = half*(qRR(8)+qLR(8))

  ! non-thermal energies
#if NENER>0
  do irad = 1,nener
     qleft (8+irad) = half*(qLL(8+irad)+qRL(8+irad))
     qright(8+irad) = half*(qRR(8+irad)+qLR(8+irad))
  end do
#endif

  zero_flux = 0.
  SELECT CASE (riemann2d)
  CASE (solver2d_roe)
     CALL riemann_roe_mhd(qleft,qright,fmean_y,zero_flux,params)
  CASE (solver2d_llf)
     CALL riemann_llf_mhd(qleft,qright,fmean_y,zero_flux,params)
  CASE (solver2d_upwind)
     CALL riemann_upwind_mhd(qleft,qright,fmean_y,zero_flux,params)
  CASE DEFAULT
     write(*,*)'unknown 2D riemann solver'
     stop
  END SELECT

  ! compute the final value of E including the 2D diffusive
  ! terms that ensure stability
  E = E + (fmean_x(7) - fmean_y(7))

END SUBROUTINE riemann2d_simple
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#endif
