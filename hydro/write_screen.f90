subroutine write_screen(r,m)
  use amr_parameters, only: ndim
  use amr_commons, only: mesh_t, run_t
  use hydro_parameters, only: nener
  use color_module, only: c_blue, c_reset
  implicit none
  type(run_t)::r
  type(mesh_t)::m

  ! Local variables
  integer::ilevel,nleaf,i
  integer,dimension(:),allocatable::ll,ii
  real(kind=8),dimension(:),allocatable::xx,dd,uu,vv,ww,pp
#if NION>0 && !defined(DO_RTZ)
  ! Ionization fractions of non-equilibrium chemistry (xHII, xHeII, xHeIII)
  integer::nx
  integer,dimension(1:3)::ix
  real(kind=8),dimension(:,:),allocatable::xi
  character(len=33)::xhead
#endif
#ifdef MHD
  real(kind=8),dimension(:),allocatable::AA,BB,CC
#endif
#if NENER>0
  integer::irad
  real(kind=8),dimension(:,:),allocatable::prad
#endif
  integer::igrid,ind
  logical::leaf
  real(kind=8)::ekin,emag,erad,dx

#ifdef HYDRO

  if(ndim>1)return

  ! Compute number of cells
  nleaf = 0
  do ilevel=r%levelmin,r%nlevelmax
     if(m%noct_tot(ilevel)>0)then
        do igrid=m%head(ilevel),m%tail(ilevel)
           do ind=1,2
              leaf = .not. m%grid(igrid)%refined(ind)
              if(leaf)nleaf=nleaf+1
           end do
        end do
     endif
  end do

  if(nleaf>0)then
     allocate(ll(1:nleaf))
     allocate(xx(1:nleaf))
     allocate(dd(1:nleaf))
     allocate(uu(1:nleaf))
     allocate(vv(1:nleaf))
     allocate(ww(1:nleaf))
#if NENER>0
     allocate(prad(1:nener,1:nleaf))
#endif
#ifdef MHD
     allocate(AA(1:nleaf))
     allocate(BB(1:nleaf))
     allocate(CC(1:nleaf))
#endif
     allocate(pp(1:nleaf))
     allocate(ii(1:nleaf))
#if NION>0 && !defined(DO_RTZ)
     nx=0
     xhead=''
     if(r%neq_chem)then
        nx=1
        ix(1)=r%iIons-1+r%ixHII
        xhead='       xHII'
        if(r%isHe)then
           nx=3
           ix(2)=r%iIons-1+r%ixHeII
           ix(3)=r%iIons-1+r%ixHeIII
           xhead='       xHII      xHeII     xHeIII'
        endif
     endif
     allocate(xi(1:3,1:nleaf))
#endif

     nleaf=0
     do ilevel=r%levelmin,r%nlevelmax
        if(m%noct_tot(ilevel)>0)then
           dx=r%boxlen/2**ilevel
           do igrid=m%head(ilevel),m%tail(ilevel)
              do ind=1,2
                 leaf=.not.m%grid(igrid)%refined(ind)
                 if(leaf)then
                    nleaf=nleaf+1
                    ll(nleaf)=m%grid(igrid)%lev
                    xx(nleaf)=(2*m%grid(igrid)%ckey(1)+ind-0.5)*dx-m%skip(1)
                    dd(nleaf)=m%uold(ind,1,igrid)
                    uu(nleaf)=m%uold(ind,2,igrid)/m%uold(ind,1,igrid)
                    vv(nleaf)=m%uold(ind,3,igrid)/m%uold(ind,1,igrid)
                    ww(nleaf)=m%uold(ind,4,igrid)/m%uold(ind,1,igrid)
                    ekin=0.5*dd(nleaf)*(uu(nleaf)**2+vv(nleaf)**2+ww(nleaf)**2)
                    emag=0.0
#ifdef MHD
                    AA(nleaf)=0.5*(m%bold(ind,1,igrid)+m%bold(ind,4,igrid))
                    BB(nleaf)=0.5*(m%bold(ind,2,igrid)+m%bold(ind,5,igrid))
                    CC(nleaf)=0.5*(m%bold(ind,3,igrid)+m%bold(ind,6,igrid))
                    emag=0.5*(AA(nleaf)**2+BB(nleaf)**2+CC(nleaf)**2)
#endif
                    erad=0.0d0
#if NENER>0
                    do irad=1,nener
                       erad=erad+m%uold(ind,5+irad,igrid)
                       prad(irad,nleaf)=(r%gamma_rad(irad)-1)*m%uold(ind,5+irad,igrid)
                    end do
#endif
                    pp(nleaf)=m%uold(ind,5,igrid)-ekin-emag-erad
                    pp(nleaf)=(r%gamma-1)*pp(nleaf)
#if NION>0 && !defined(DO_RTZ)
                    do i=1,nx
                       xi(i,nleaf)=m%uold(ind,ix(i),igrid)/m%uold(ind,1,igrid)
                       if(ABS(xi(i,nleaf))<1d-99)xi(i,nleaf)=0.0D0
                    end do
#endif

                    ! Values below 1d-99 would need a 3-digit exponent in the E10.3 format
                    if(ABS(uu(nleaf))<1d-99)uu(nleaf)=0.0D0
                    if(ABS(vv(nleaf))<1d-99)vv(nleaf)=0.0D0
                    if(ABS(ww(nleaf))<1d-99)ww(nleaf)=0.0D0
#ifdef MHD
                    if(ABS(AA(nleaf))<1d-99)AA(nleaf)=0.0D0
                    if(ABS(BB(nleaf))<1d-99)BB(nleaf)=0.0D0
                    if(ABS(CC(nleaf))<1d-99)CC(nleaf)=0.0D0
#endif
                 endif
              end do
           end do
        endif
     end do
     call quick_sort_dp(xx,ii,nleaf)

#ifdef MHD
     write(*,*)'================================================'
     write(*,114)c_blue(),nleaf,c_reset()
     write(*,*)'lev      x           d          u          v          w          P          A          B          C'
     do i=1,nleaf
        write(*,115)      &
             & ll(ii(i)), &
             & xx(i),     &
             & dd(ii(i)), &
             & uu(ii(i)), &
             & vv(ii(i)), &
             & ww(ii(i)), &
             & pp(ii(i)), &
             & AA(ii(i)), &
             & BB(ii(i)), &
             & CC(ii(i))
     end do
     write(*,*)'================================================'
#elif NENER>0
     write(*,*)'================================================'
     write(*,114)c_blue(),nleaf,c_reset()
     write(*,*)'lev      x           d          u          P        P_rad'
     do i=1,nleaf
        write(*,116)      &
             & ll(ii(i)), &
             & xx(i),     &
             & dd(ii(i)), &
             & uu(ii(i)), &
             & pp(ii(i)), &
             & prad(1,ii(i))
     end do
     write(*,*)'================================================'
#elif NION>0 && !defined(DO_RTZ)
     write(*,*)'================================================'
     write(*,114)c_blue(),nleaf,c_reset()
     write(*,*)'lev      x           d          u          P'//trim(xhead)
     do i=1,nleaf
        write(*,117)      &
             & ll(ii(i)), &
             & xx(i),     &
             & dd(ii(i)), &
             & uu(ii(i)), &
             & pp(ii(i)), &
             & xi(1:nx,ii(i))
     end do
     write(*,*)'================================================'
#else
     write(*,*)'================================================'
     write(*,114)c_blue(),nleaf,c_reset()
     write(*,*)'lev      x           d          u          P'
     do i=1,nleaf
        write(*,113)      &
             & ll(ii(i)), &
             & xx(i),     &
             & dd(ii(i)), &
             & uu(ii(i)), &
             & pp(ii(i))
     end do
     write(*,*)'================================================'
#endif
     deallocate(ll)
     deallocate(xx)
     deallocate(dd)
     deallocate(uu)
     deallocate(vv)
     deallocate(ww)
     deallocate(pp)
     deallocate(ii)
#if NION>0 && !defined(DO_RTZ)
     deallocate(xi)
#endif
#if NENER>0
     deallocate(prad)
#endif
#ifdef MHD
     deallocate(AA)
     deallocate(BB)
     deallocate(CC)
#endif

  endif

113 format(i3,1x,1pe12.5,1x,3(1pe10.3,1x))
114 format(' ',a,'Output ',i5,' cells',a)
115 format(i3,1x,1pe12.5,1x,8(1pe10.3,1x))
116 format(i3,1x,1pe12.5,1x,4(1pe10.3,1x))
117 format(i3,1x,1pe12.5,1x,6(1pe10.3,1x))

#endif

end subroutine write_screen
