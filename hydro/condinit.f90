!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(r,g,x,q,dx,nn)
  use amr_parameters, only: ndim, nvector
  use hydro_parameters, only: nvar, nener
  use amr_commons, only: run_t, global_t
  use input_hydro_condinit_module, only: region_condinit
  use constants, only: kB, mH, M_sun, factG_in_cgs
  implicit none
  type(run_t)::r
  type(global_t)::g
  integer ::nn                            ! Number of cells
  real(kind=8)::dx                            ! Cell size
#ifdef MHD
  real(kind=8),dimension(1:nvector,1:nvar+3-ndim)::q ! Primitive variables
#else
  real(kind=8),dimension(1:nvector,1:nvar)::q ! Primitive variables
#endif
  real(kind=8),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user (aka code) units:
  ! x(i,1:ndim) are in [0,box_size]**ndim.
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:4):u,v,w and Q(i,5): P.
  ! If nvar >= 6, remaining variables are treated as passive
  ! scalars or non-thermal energies in the hydro solver.
  ! For 1D MHD, Q(i,nvar+1) is By and Q(i,nvar+2) is Bz.
  ! For 2D MHD, Q(i,nvar+1) is Bz.
  ! Q(:,:) are in user (aka code) units.
  !================================================================
#define COEUR 1
#define INSTA 2
#define DOUBLEMACH 3
#define OT 4
#define PONO 5
#define ABC 6
#define CURRENTSHEET 7
#define RTZEQM 8
#define PANCAKE 9
#define ALFVENWAVE 10
#define COLLAPSE 11

  integer::i
  real(kind=8)::xx,yy,zz,rr,theta,pi,xcenter,ttmin,ttmax
#if INIT==COEUR
  real(kind=8)::r2,rx,ry,rz,d,p,vx,vy,vz,r_trunc,r2_trunc,c2
  real(kind=8)::omega_code,AU,Msol,pi,M,sigma,r_min,r2_min,omega_const,r_vortex,invr2_vortex
  real(kind=8)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m
#elif INIT==INSTA
  integer::id,iu,iv,iw,ip,ix,iy
  real(kind=8)::x0,lambday,ky,lambdaz,kz,rho1,rho2,p0,v0,v1,v2
#elif INIT==DOUBLEMACH
  integer::id,iu,iv,iw,ip
  real(kind=8)::pi,xp
#elif INIT==OT
  real(kind=8)::pi,xc,yc
#elif INIT==PONO
  real(kind=8)::vx,vy,vz,tt,omega,R0,twopi
#elif INIT==ABC
  real(kind=8)::vx,vy,vz,A0,twopi
#elif INIT==CURRENTSHEET
  real(kind=8)::pi,xc,yc,beta,v0
#elif INIT==RTZEQM
  real(kind=8)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m
#elif INIT==PANCAKE
  real(kind=8)::pi,del_ini
  real(kind=8)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m
#elif INIT==ALFVENWAVE
  real(kind=8)::pi,del_ini
  real(kind=8)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m
#elif INIT==COLLAPSE
  real(kind=8)::x0,y0,z0,xx,yy,zz,rc,rs,phi
  real(kind=8)::r0,d0,p0,omega0,B0,mass_c_cu,scale_m
  real(kind=8)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(kind=8),parameter::pi=3.14159265358979323846d0
  real(kind=8),parameter::delta_rho=0.1d0             ! m=2 density perturbation amplitude
  real(kind=8),parameter::alpha_dense_core=0.1d0      ! thermal-to-gravitational energy ratio
  real(kind=8),parameter::beta_dense_core=0.01d0      ! rotational-to-gravitational energy ratio
  real(kind=8),parameter::crit_dense_core=0.08d0      ! 1/mu for Bfield strength
  real(kind=8),parameter::theta_mag=0.0d0             ! angle in degrees for rotation misalignment between Bfield and rotation
  real(kind=8),parameter::mass_c=1.0d0                ! mass of the collapsing core in solar masses
  real(kind=8),parameter::Mach=0.0d0
  real(kind=8),parameter::T_eos=10.0d0
  real(kind=8),parameter::mu_gas=2.31d0
#else
  ! Call built-in initial condition generator
  call region_condinit(r,g,x,q,dx,nn)
#endif

  ! Add here, if you wish, some user-defined initial conditions
  ! ........

#if INIT==COEUR
  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=scale_d*scale_l**3
  ! constants
  AU=1.49598d13
  Msol= 1.98892d33
  pi=3.14159
  ! mass, radius, and ratio of rotational to gravitational energy
  r_trunc=25*4000.*AU/scale_l
  r_min=10.*AU/scale_l
  r_vortex=4000.*AU/scale_l
  M=100.*Msol/scale_m
  sigma=M/(4*pi*r_trunc)
  r2_trunc=r_trunc**2
  r2_min=r_min**2
  invr2_vortex=1./r_vortex**2
  omega_const=0.1*sqrt(1./r2_trunc+invr2_vortex)*sqrt(M/r_trunc)
  c2=(18939.2/(scale_l/scale_t))**2
  do i=1,nn
     rx=x(i,1)-r%box_size(1)/2.
     ry=x(i,2)-r%box_size(2)/2.
     rz=x(i,3)-r%box_size(3)/2.
     !density
     r2=rx**2+ry**2+rz**2
     d=sigma/(r2+r2_min)
     omega_code=omega_const/sqrt(1.+invr2_vortex*r2)
     if (r2>=r2_trunc)then
        d=d*1.e-4
        omega_code=omega_code/sqrt(r2)*exp(10.*(r2_trunc-r2))
     end if
     !pressure
     p=d*c2
     !velocity
     vx=-omega_code*ry
     vy=omega_code*rx
     vz=0.
     ! primitive variables
     q(i,1)=d
     q(i,2)=vx
     q(i,3)=vy
     q(i,4)=vz
     q(i,5)=p
  end do
#endif

#if INIT==INSTA
  id=1; iu=2; iv=3; iw=4; ip=5
  x0=r%x_center(1)
  if(r%constant_gravity(2) .ne. 0)then
     ix=2
     iy=1
     iu=3
     iv=2
  else
     ix=1
     iy=2
     iu=2
     iv=3
  endif

  lambday=0.25
  ky=2.*acos(-1.0d0)/lambday
  lambdaz=0.25
  kz=2.*acos(-1.0d0)/lambdaz
  rho1=r%d_region(1)
  rho2=r%d_region(2)
  v1=r%v_region(1)
  v2=r%v_region(2)
  v0=0.1
  p0=10.
  do i=1,nn
     if(x(i,ix) < x0)then
        q(i,id)=rho1
        q(i,iu)=0.0
        q(i,iu)=v0*cos(ky*(x(i,iy)-lambday/2.))*exp(+ky*(x(i,ix)-x0))
        q(i,iv)=v1
        q(i,iw)=0.0D0
        q(i,ip)=p0+rho1*r%constant_gravity(ix)*x(i,ix)
     else
        q(i,id)=rho2
        q(i,iu)=0.0
        q(i,iu)=v0*cos(ky*(x(i,iy)-lambday/2.))*exp(-ky*(x(i,ix)-x0))
        q(i,iv)=v2
        q(i,iw)=0.0D0
        q(i,ip)=p0+rho1*r%constant_gravity(ix)*x0+rho2*r%constant_gravity(ix)*(x(i,ix)-x0)
     endif
  end do
#endif

#if INIT==DOUBLEMACH
  id=1; iu=2; iv=3; iw=4; ip=5
  pi=acos(-1.0d0)
  do i=1,nn
     xp=x(i,1)-x(i,2)/tan(pi/3.0)-10./sin(pi/3.0)*g%t
     if(xp<1./6.)then
        q(i,id)=8.
        q(i,iu)=7.145
        q(i,iv)=-4.125
        q(i,iw)=0
        q(i,ip)=116.5
     else
        q(i,id)=r%gamma
        q(i,iu)=0.0
        q(i,iv)=0.0
        q(i,iw)=0.0
        q(i,ip)=1.0
     endif
  end do
#endif

#if INIT==OT
  pi=acos(-1.0d0)
  do i=1,nn
     xc=x(i,1)
     yc=x(i,2)
     q(i,1)=25.0/(36.0*pi)
     q(i,2)=-sin(2.0*pi*yc)
     q(i,3)=+sin(2.0*pi*xc)
     q(i,4)=0.0
     q(i,5)=5.0/(12.0*pi)
     q(i,nvar+1)=0.0 ! Bz
  end do
#endif

#if INIT==PONO
  R0=1.0
  twopi=2d0*ACOS(-1d0)
  do i=1,nn
     q(i,1)=1.0
     q(i,5)=1.0*(r%gamma-1.0)
     xx=x(i,1)-r%box_size(1)/2.
     yy=x(i,2)-r%box_size(2)/2.
     rr = SQRT(xx**2+yy**2)
     if(rr < 1.0)then
        omega=0.609711
        vz=0.792624
     else
        omega=0.0
        vz=0.0
     endif
     if(rr > 0.0)then
        if(yy > 0.0)then
           tt=acos(xx/rr)
        else
           tt=-acos(xx/rr)+twopi
        endif
        vx=-sin(tt)*rr*omega
        vy=+cos(tt)*rr*omega
     else
        vx=0.0
        vy=0.0
     endif
     q(i,2)=vx
     q(i,3)=vy
     q(i,4)=vz
  end do
#endif

#if INIT==ABC
  A0=1.0
  twopi=2d0*ACOS(-1d0)
  do i=1,nn
     q(i,1)=1.0
     q(i,5)=1.0*(r%gamma-1.0)
     xx=x(i,1)-r%box_size(1)/2.
     yy=x(i,2)-r%box_size(2)/2.
     zz=x(i,3)-r%box_size(3)/2.
     vx=A0*(cos(twopi*yy)+sin(twopi*zz))
     vy=A0*(sin(twopi*xx)+cos(twopi*zz))
     vz=A0*(cos(twopi*xx)+sin(twopi*yy))
     q(i,2)=vx
     q(i,3)=vy
     q(i,4)=vz
  end do
#endif

#if INIT==CURRENTSHEET
  pi = acos(-1.0d0)
  beta = 0.1
  v0 = 0.1
  do i = 1,nn
     xc = x(i,1)
     yc = x(i,2)
     q(i,1) = 1.0
     q(i,2) = v0*sin(pi*yc)
     q(i,3) = 0.0
     q(i,4) = 0.0
     q(i,5) = 0.5*beta
     q(i,nvar+1) = 0.0 ! Bz
  end do
#endif

#if INIT==RTZEQM
  ! get the units
  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  ! Smoothly interpolate gas density between
  ! 1e-3 and 1e5, fix T to 10^4, and convert everything to code units
  ! note that this assumes a boxsize of 1 and Nx = Ny and unigrid
  do i = 1,nn
     q(i,1) = (10.d0**(x(i,1) * 8.d0 - 3.d0)) / scale_nH
     q(i,2) = 0.0 ! Vx
     q(i,3) = 0.0 ! Vy
     q(i,4) = 0.0 ! Vz
     q(i,5) = 1.d4 / scale_T2 ! Temperature is 10^4 K
  end do
#endif

#if INIT==PANCAKE
  pi = acos(-1.0d0)
  del_ini = 0.1
  ! get cgs units
  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  do i = 1,nn
     q(i,1) = g%omega_b/g%omega_m/(1+del_ini*COS(2.0d0*pi*x(i,1)))
     q(i,2) = del_ini*g%vfact(1)*SIN(2.0d0*pi*x(i,1))/(2.0d0*pi)
     q(i,3) = 0.0 ! Vy
     q(i,4) = 0.0 ! Vz
     q(i,5) = 100./scale_T2 ! Temperature is 10^2 K
  end do
#endif

#if INIT==ALFVENWAVE
  pi = acos(-1.0d0)
  del_ini = 0.1
  ! get cgs units
  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  do i = 1,nn
     q(i,1) = 1 ! rho
     q(i,2) = 0 ! Vx
     q(i,3) = 0.1*SIN(2.0d0*pi*x(i,1)) ! Vy
     q(i,4) = 0.1*COS(2.0d0*pi*x(i,1)) ! Vz
     q(i,5) = 0.1 ! Pressure
     q(i,6) = 0.1*SIN(2.0d0*pi*x(i,1)) ! By
     q(i,7) = 0.1*COS(2.0d0*pi*x(i,1)) ! Bz
  end do
#endif

#if INIT==COLLAPSE
  if(abs(theta_mag)>0.0d0)then
     write(*,*)'COLLAPSE condinit currently supports theta_mag=0 only'
     stop
  endif
  if(abs(Mach)>0.0d0)then
     write(*,*)'COLLAPSE condinit currently supports Mach=0 only'
     stop
  endif

  x0=0.5d0*r%box_size(1)
  y0=0.5d0*r%box_size(2)
  z0=0.5d0*r%box_size(3)

  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=scale_d*scale_l**ndim
  mass_c_cu=mass_c*(M_sun/scale_m)

  r0=alpha_dense_core*2.0d0*factG_in_cgs*mass_c_cu*scale_m*mu_gas*mH &
       & /(5.0d0*kB*T_eos)/scale_l
  d0=3.0d0*mass_c_cu/(4.0d0*pi*r0**3)
  omega0=sqrt(beta_dense_core*4.0d0*pi*d0)
  p0=alpha_dense_core*d0*d0*r0*r0*8.0d0*pi/15.0d0
  B0=sqrt(4.0d0*pi/5.0d0)/0.53d0*crit_dense_core*d0*r0

#ifdef MHD
#if NDIM==3
  r%A_ave=0.0d0
  r%B_ave=0.0d0
  r%C_ave=B0
#endif
#endif

  do i=1,nn
     xx=x(i,1)-x0
#if NDIM>1
     yy=x(i,2)-y0
#else
     yy=0.0d0
#endif
#if NDIM>2
     zz=x(i,3)-z0
#else
     zz=0.0d0
#endif
     rc=sqrt(xx**2+yy**2)
     rs=sqrt(xx**2+yy**2+zz**2)

     if(rc>0.0d0)then
        phi=atan2(yy,xx)
     else
        phi=0.0d0
     endif

     if(rs<=r0)then
        q(i,1)=d0*(1.0d0+delta_rho*cos(2.0d0*phi))
        q(i,2)= omega0*yy
        q(i,3)=-omega0*xx
        q(i,4)=0.0d0
        q(i,5)=p0
     else
        q(i,1)=d0/100.0d0
        q(i,2)=0.0d0
        q(i,3)=0.0d0
        q(i,4)=0.0d0
        q(i,5)=p0/100.0d0
     endif

#ifdef MHD
#if NDIM==1
     q(i,nvar+1)=0.0d0
     q(i,nvar+2)=B0
#endif
#if NDIM==2
     q(i,nvar+1)=B0
#endif
#endif
  end do
#endif

  ! Compute entropy if needed
  if(r%entropy)then
     q(1:nn,r%ientropy)=q(1:nn,5)/q(1:nn,1)**r%gamma
  endif

  ! Compute metallicity if needed
  if(r%metal)then
     q(1:nn,r%imetal)=r%z_ave*0.02
  endif

#ifdef DO_CR
  if(r%cr_test_setup=='streaming_triangle') then
     q(1:nn,6)=(2d0-1d0*sqrt((x(1:nn,1)-r%box_size(1)*0.5d0)**2))/3.
  else if(r%cr_test_setup=='diffusion') then
    q(1:nn,6)=exp(-40d0*(x(1:nn,1)-r%box_size(1)*0.5d0)**2)/3.
  else if(r%cr_test_setup=='1d_cr_cloud') then
     q(1:nn,1)=0.1d0+(10d0-0.1d0)*(1d0+tanh((x(1:nn,1)-200d0)/25d0)) &
          &                       *(1d0+tanh((200d0-x(1:nn,1))/25d0))
  else if(r%cr_test_setup=='circular_diffusion') then
     ! CR energy: enhanced on one arc of the loop (around the z-axis).
     ! atan2 replaces atan(yy/xx) to avoid a divide-by-zero FPE at xx=0;
     ! it is identical in the xx>0 region that the arc condition selects.
     pi=acos(-1d0)
     ttmin=-pi/12d0
     ttmax= pi/12d0
     xcenter=r%box_size(1)*0.5d0
     do i=1,nn
        xx=x(i,1)-xcenter
#if NDIM>1
        yy=x(i,2)-xcenter
#endif
        rr=sqrt(xx**2+yy**2)
        theta=atan2(yy,xx)
        if(rr>0.25d0*r%box_size(1) .and. rr<0.35d0*r%box_size(1) .and. theta>ttmin .and. &
            & theta<ttmax .and. xx>0d0 .and. zz>-0.1d0*r%box_size(1) .and. zz<0.1d0*r%box_size(1))then
          q(i,6)=1.2d1
        else
          q(i,6)=1.0d1
        endif
     end do
  endif
#endif


end subroutine condinit
