module boundaries
contains
subroutine get_bound(s,ix,ilevel,ibound)
  use amr_parameters, only: ndim,twotondim
  use ramses_commons, only: ramses_t
  implicit none
  type(ramses_t)::s
  integer(kind=8),dimension(1:ndim)::ix
  integer::ibound,ilevel

  logical::in_bound
  integer::i,idim

  associate(r=>s%r,g=>s%g,m=>s%m,mdl=>s%mdl)

    ibound = 0
    do i = 1, r%nbound
       in_bound = .true.
       in_bound = in_bound .and. ix(1) .ge. m%bound_ckey_min(1,i,ilevel) .and. ix(1) .lt. m%bound_ckey_max(1,i,ilevel)
#if NDIM>1
       in_bound = in_bound .and. ix(2) .ge. m%bound_ckey_min(2,i,ilevel) .and. ix(2) .lt. m%bound_ckey_max(2,i,ilevel)
#endif
#if NDIM>2
       in_bound = in_bound .and. ix(3) .ge. m%bound_ckey_min(3,i,ilevel) .and. ix(3) .lt. m%bound_ckey_max(3,i,ilevel)
#endif
       if (in_bound) then
          ibound = i
          exit
       endif
    end do

  end associate

end subroutine get_bound
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_bound_refine(r,g,m,igrid,igrid_ref,ibound)
  use amr_parameters, only: ndim, twotondim, nvector
  use hydro_parameters, only: nvar, nener
  use rt_parameters, only: nrtvar, nrtgrp
  use cr_parameters, only: ncruvar, ncrgrp
  use amr_commons, only: run_t, global_t, mesh_t, oct
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer::igrid
  integer::igrid_ref
  integer::ibound
  real(kind=8)::scale_nH, scale_T2, scale_l, scale_d, scale_t, scale_v
  real(kind=8)::scale_np, scale_fp, dx_cgs, dt_cgs

  integer,dimension(1:8,1:3)::ref_right=reshape(&
       & (/ 0,1,0,3,0,5,0,7,&
       &    0,0,1,2,0,0,5,6,&
       &    0,0,0,0,1,2,3,4/),(/8,3/))
  integer,dimension(1:8,1:3)::ref_left=reshape(&
       & (/ 2,0,4,0,6,0,8,0,&
       &    3,4,0,0,7,8,0,0,&
       &    5,6,7,8,0,0,0,0/),(/8,3/))
  integer,dimension(1:8,1:3)::ind1_right=reshape(&
       & (/ 2,1,4,3,6,5,8,7,&
       &    3,4,1,2,7,8,5,6,&
       &    5,6,7,8,1,2,3,4/),(/8,3/))
  integer,dimension(1:8,1:3)::ind1_left=reshape(&
       & (/ 2,1,4,3,6,5,8,7,&
       &    3,4,1,2,7,8,5,6,&
       &    5,6,7,8,1,2,3,4/),(/8,3/))
  integer,dimension(1:8,1:3)::ind2_right=reshape(&
       & (/ 1,1,3,3,5,5,7,7,&
       &    1,2,1,2,5,6,5,6,&
       &    1,2,3,4,1,2,3,4/),(/8,3/))
  integer,dimension(1:8,1:3)::ind2_left=reshape(&
       & (/ 2,2,4,4,6,6,8,8,&
       &    3,4,3,4,7,8,7,8,&
       &    5,6,7,8,5,6,7,8/),(/8,3/))

  integer::idim, ind, ivar, irad
  integer::type, dir, shift, nstride
  real(kind=8)::reverse, ek_bound
  real(kind=8),dimension(1:nvector,1:ndim)::xx
  real(kind=8),dimension(1:nvector,1:nvar)::uu
  real(kind=8),dimension(1:nvector,1:nvar)::qq
  real(kind=8),dimension(1:nvector,1:ndim)::ff
  real(kind=8),dimension(1:nvector)::phi
  real(kind=8)::dx,rr,vx,vy,vz,pp,eint,ekin,emag,erad
#if defined(MHD) && NDIM>1
  integer::alt
  real(kind=8)::bface,dbn
#endif
#ifdef DO_CR
  real(kind=8),dimension(1:nvector,1:ncruvar)::cruu
#endif

  type = r%bound_type(ibound)
  dir = r%bound_dir(ibound)
  shift = r%bound_shift(ibound)

  ! Mesh size at level ilevel in code units
  dx=r%boxlen/2**m%grid(igrid)%lev

  ! Set refinement map
  if(shift==+1)then
     do ind=1,twotondim
        if(ref_right(ind,dir)>0)then
           m%grid(igrid)%refined(ind)=m%grid(igrid_ref)%refined(ref_right(ind,dir))
        else
           m%grid(igrid)%refined(ind)=.false.
        endif
     end do
  end if
  if(shift==-1)then
     do ind=1,twotondim
        if(ref_left(ind,dir)>0)then
           m%grid(igrid)%refined(ind)=m%grid(igrid_ref)%refined(ref_left(ind,dir))
        else
           m%grid(igrid)%refined(ind)=.false.
        endif
     end do
  end if

#ifdef HYDRO

  ! Reflexive BC
  if(type == 1)then

     if(shift==+1)then
        do ivar=1,nvar
           reverse=1
           if(ivar==1+dir)reverse=-1
           do ind=1,twotondim
              m%uold(ind,ivar,igrid)=m%uold(ind1_right(ind,dir),ivar,igrid_ref)*reverse
           end do
        end do
#ifdef MHD
#if NDIM==1
        do ind=1,twotondim
           m%bold(ind,1,igrid)=r%A_ave
           m%bold(ind,4,igrid)=r%A_ave
           m%bold(ind,2,igrid)=m%bold(ind1_right(ind,dir),2,igrid_ref)
           m%bold(ind,5,igrid)=m%bold(ind1_right(ind,dir),5,igrid_ref)
           m%bold(ind,3,igrid)=m%bold(ind1_right(ind,dir),3,igrid_ref)
           m%bold(ind,6,igrid)=m%bold(ind1_right(ind,dir),6,igrid_ref)
        end do
#endif
#if NDIM>1
        ! Tangential field copied, normal field mirrored about its boundary face value,
        ! and the total energy given the ghost cell's own magnetic energy
        do ind=1,twotondim
           emag=0.125d0*sum((m%bold(ind1_right(ind,dir),1:3,igrid_ref)+m%bold(ind1_right(ind,dir),4:6,igrid_ref))**2)
           m%bold(ind,1:6,igrid)=m%bold(ind1_right(ind,dir),1:6,igrid_ref)
           bface=m%bold(ind2_right(ind,dir),dir,igrid_ref)
           m%bold(ind,dir,igrid)=2d0*bface-m%bold(ind1_right(ind,dir),dir+3,igrid_ref)
           m%bold(ind,dir+3,igrid)=2d0*bface-m%bold(ind1_right(ind,dir),dir,igrid_ref)
           m%uold(ind,5,igrid)=m%uold(ind,5,igrid)-emag+0.125d0*sum((m%bold(ind,1:3,igrid)+m%bold(ind,4:6,igrid))**2)
        end do
#endif
#endif
#ifdef DO_RT
        do ivar=1,nrtvar
           reverse=1
           if(mod(ivar-1, ndim+1) == dir) reverse=-1
           do ind=1,twotondim
              m%rtuold(ind,ivar,igrid)=m%rtuold(ind1_right(ind,dir),ivar,igrid_ref)*reverse
           end do
        end do
#endif
#ifdef DO_CR
        do ivar=1,ncruvar
           reverse=1
           if((mod(ivar-1, ndim)+1) == dir) reverse=-1
           do ind=1,twotondim
              m%cruold(ind,ivar,igrid)=m%cruold(ind1_right(ind,dir),ivar,igrid_ref)*reverse
           end do
        end do
#endif
     endif

     if(shift==-1)then
        do ivar=1,nvar
           reverse=1
           if(ivar==1+dir)reverse=-1
           do ind=1,twotondim
              m%uold(ind,ivar,igrid)=m%uold(ind1_left(ind,dir),ivar,igrid_ref)*reverse
           end do
        end do
#ifdef MHD
#if NDIM==1
        do ind=1,twotondim
           m%bold(ind,1,igrid)=r%A_ave
           m%bold(ind,4,igrid)=r%A_ave
           m%bold(ind,2,igrid)=m%bold(ind1_left(ind,dir),2,igrid_ref)
           m%bold(ind,5,igrid)=m%bold(ind1_left(ind,dir),5,igrid_ref)
           m%bold(ind,3,igrid)=m%bold(ind1_left(ind,dir),3,igrid_ref)
           m%bold(ind,6,igrid)=m%bold(ind1_left(ind,dir),6,igrid_ref)
        end do
#endif
#if NDIM>1
        ! Tangential field copied, normal field mirrored about its boundary face value,
        ! and the total energy given the ghost cell's own magnetic energy
        do ind=1,twotondim
           emag=0.125d0*sum((m%bold(ind1_left(ind,dir),1:3,igrid_ref)+m%bold(ind1_left(ind,dir),4:6,igrid_ref))**2)
           m%bold(ind,1:6,igrid)=m%bold(ind1_left(ind,dir),1:6,igrid_ref)
           bface=m%bold(ind2_left(ind,dir),dir+3,igrid_ref)
           m%bold(ind,dir,igrid)=2d0*bface-m%bold(ind1_left(ind,dir),dir+3,igrid_ref)
           m%bold(ind,dir+3,igrid)=2d0*bface-m%bold(ind1_left(ind,dir),dir,igrid_ref)
           m%uold(ind,5,igrid)=m%uold(ind,5,igrid)-emag+0.125d0*sum((m%bold(ind,1:3,igrid)+m%bold(ind,4:6,igrid))**2)
        end do
#endif
#endif
#ifdef DO_RT
        do ivar=1,nrtvar
           reverse=1
           if(mod(ivar-1, ndim+1) == dir) reverse=-1
           do ind=1,twotondim
              m%rtuold(ind,ivar,igrid)=m%rtuold(ind1_left(ind,dir),ivar,igrid_ref)*reverse
           end do
        end do
#endif
#ifdef DO_CR
        do ivar=1,ncruvar
           reverse=1
           if((mod(ivar-1, ndim)+1) == dir) reverse=-1
           do ind=1,twotondim
              m%cruold(ind,ivar,igrid)=m%cruold(ind1_left(ind,dir),ivar,igrid_ref)*reverse
           end do
        end do
#endif
     endif

  endif

  ! Zero gradient BC
  if(type == 2)then

     if(shift==+1)then
        do ivar=1,nvar
           do ind=1,twotondim
              m%uold(ind,ivar,igrid)=m%uold(ind2_right(ind,dir),ivar,igrid_ref)
           end do
        end do
#ifdef MHD
#if NDIM==1
        do ind=1,twotondim
           m%bold(ind,1,igrid)=r%A_ave
           m%bold(ind,4,igrid)=r%A_ave
           m%bold(ind,2,igrid)=m%bold(ind2_right(ind,dir),2,igrid_ref)
           m%bold(ind,5,igrid)=m%bold(ind2_right(ind,dir),5,igrid_ref)
           m%bold(ind,3,igrid)=m%bold(ind2_right(ind,dir),3,igrid_ref)
           m%bold(ind,6,igrid)=m%bold(ind2_right(ind,dir),6,igrid_ref)
        end do
#endif
#if NDIM>1
        ! Tangential field copied, normal field extrapolated linearly from the boundary cell,
        ! and the total energy given the ghost cell's own magnetic energy
        do ind=1,twotondim
           emag=0.125d0*sum((m%bold(ind2_right(ind,dir),1:3,igrid_ref)+m%bold(ind2_right(ind,dir),4:6,igrid_ref))**2)
           m%bold(ind,1:6,igrid)=m%bold(ind2_right(ind,dir),1:6,igrid_ref)
           alt=mod((ind-1)/2**(dir-1),2)-2
           dbn=m%bold(ind2_right(ind,dir),dir+3,igrid_ref)-m%bold(ind2_right(ind,dir),dir,igrid_ref)
           m%bold(ind,dir,igrid)=m%bold(ind,dir,igrid)+alt*dbn
           m%bold(ind,dir+3,igrid)=m%bold(ind,dir+3,igrid)+alt*dbn
           m%uold(ind,5,igrid)=m%uold(ind,5,igrid)-emag+0.125d0*sum((m%bold(ind,1:3,igrid)+m%bold(ind,4:6,igrid))**2)
        end do
#endif
#endif
#ifdef DO_RT
        do ivar=1,nrtvar
           do ind=1,twotondim
              m%rtuold(ind,ivar,igrid)=m%rtuold(ind2_right(ind,dir),ivar,igrid_ref)
           end do
        end do
#endif
#ifdef DO_CR
        do ivar=1,ncruvar
           do ind=1,twotondim
              m%cruold(ind,ivar,igrid)=m%cruold(ind2_right(ind,dir),ivar,igrid_ref)
           end do
        end do
#endif
     endif

     if(shift==-1)then
        do ivar=1,nvar
           do ind=1,twotondim
              m%uold(ind,ivar,igrid)=m%uold(ind2_left(ind,dir),ivar,igrid_ref)
           end do
        end do
#ifdef MHD
#if NDIM==1
        do ind=1,twotondim
           m%bold(ind,1,igrid)=r%A_ave
           m%bold(ind,4,igrid)=r%A_ave
           m%bold(ind,2,igrid)=m%bold(ind2_left(ind,dir),2,igrid_ref)
           m%bold(ind,5,igrid)=m%bold(ind2_left(ind,dir),5,igrid_ref)
           m%bold(ind,3,igrid)=m%bold(ind2_left(ind,dir),3,igrid_ref)
           m%bold(ind,6,igrid)=m%bold(ind2_left(ind,dir),6,igrid_ref)
        end do
#endif
#if NDIM>1
        ! Tangential field copied, normal field extrapolated linearly from the boundary cell,
        ! and the total energy given the ghost cell's own magnetic energy
        do ind=1,twotondim
           emag=0.125d0*sum((m%bold(ind2_left(ind,dir),1:3,igrid_ref)+m%bold(ind2_left(ind,dir),4:6,igrid_ref))**2)
           m%bold(ind,1:6,igrid)=m%bold(ind2_left(ind,dir),1:6,igrid_ref)
           alt=mod((ind-1)/2**(dir-1),2)+1
           dbn=m%bold(ind2_left(ind,dir),dir+3,igrid_ref)-m%bold(ind2_left(ind,dir),dir,igrid_ref)
           m%bold(ind,dir,igrid)=m%bold(ind,dir,igrid)+alt*dbn
           m%bold(ind,dir+3,igrid)=m%bold(ind,dir+3,igrid)+alt*dbn
           m%uold(ind,5,igrid)=m%uold(ind,5,igrid)-emag+0.125d0*sum((m%bold(ind,1:3,igrid)+m%bold(ind,4:6,igrid))**2)
        end do
#endif
#endif
#ifdef DO_RT
        do ivar=1,nrtvar
           do ind=1,twotondim
              m%rtuold(ind,ivar,igrid)=m%rtuold(ind2_left(ind,dir),ivar,igrid_ref)
           end do
        end do
#endif
#ifdef DO_CR
        do ivar=1,ncruvar
           do ind=1,twotondim
              m%cruold(ind,ivar,igrid)=m%cruold(ind2_left(ind,dir),ivar,igrid_ref)
           end do
        end do
#endif
     endif

  endif

  ! Imposed BC from namelist constant values
  if(type == 3)then

     do ind=1,twotondim

        m%uold(ind,1,igrid)=r%d_bound(ibound)
        m%uold(ind,2,igrid)=r%d_bound(ibound)*r%u_bound(ibound)
        ek_bound=0.5d0*r%d_bound(ibound)*r%u_bound(ibound)**2
        m%uold(ind,3,igrid)=r%d_bound(ibound)*r%v_bound(ibound)
        ek_bound=ek_bound+0.5d0*r%d_bound(ibound)*r%v_bound(ibound)**2
        m%uold(ind,4,igrid)=r%d_bound(ibound)*r%w_bound(ibound)
        ek_bound=ek_bound+0.5d0*r%d_bound(ibound)*r%w_bound(ibound)**2
#if NENER>0
        do ivar=1,nener
           m%uold(ind,5+ivar,igrid)=r%prad_bound(ibound,ivar)/(r%gamma_rad(ivar)-1.0d0)
           ek_bound=ek_bound+r%prad_bound(ibound,ivar)/(r%gamma_rad(ivar)-1.0d0)
        enddo
#endif
        m%uold(ind,5,igrid)=r%p_bound(ibound)/(r%gamma-1.0d0)+ek_bound
#if NVAR>5+NENER
        do ivar=6+nener,nvar
           m%uold(ind,ivar,igrid)=r%d_bound(ibound)*r%var_bound(ibound,ivar-5-nener)
        end do
#endif
#ifdef MHD
#if NDIM==1
        m%bold(ind,1,igrid)=r%A_ave
        m%bold(ind,4,igrid)=r%A_ave
        m%bold(ind,2,igrid)=r%B_bound(ibound)
        m%bold(ind,5,igrid)=r%B_bound(ibound)
        m%bold(ind,3,igrid)=r%C_bound(ibound)
        m%bold(ind,6,igrid)=r%C_bound(ibound)
#endif
#if NDIM==2
        m%bold(ind,3,igrid)=r%C_bound(ibound)
        m%bold(ind,6,igrid)=r%C_bound(ibound)
#endif
#endif
     end do
#ifdef DO_RT
     call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
     call rt_units(r,g,scale_np,scale_fp)
     do ind=1,twotondim
        do ivar=1,nrtgrp
            m%rtuold(ind,1+(ivar-1)*(ndim+1),igrid)= &
              r%rt_n_bound(ibound,ivar)/g%rt_c(m%grid(igrid)%lev)/scale_fp
            m%rtuold(ind,2+(ivar-1)*(ndim+1),igrid)= &
              r%rt_n_bound(ibound,ivar)*r%rt_u_bound(ibound,ivar)/scale_fp
            if(ndim>1) m%rtuold(ind,3+(ivar-1)*(ndim+1),igrid)= &
                r%rt_n_bound(ibound,ivar)*r%rt_v_bound(ibound,ivar)/scale_fp
            if(ndim>2) m%rtuold(ind,4+(ivar-1)*(ndim+1),igrid)= &
                r%rt_n_bound(ibound,ivar)*r%rt_w_bound(ibound,ivar)/scale_fp
        end do
     end do
#endif
#ifdef DO_CR
     do ind=1,twotondim
        do ivar=1,ncrgrp
            m%cruold(ind,1+(ivar-1)*ndim,igrid)= r%cr_fx_bound(ibound,ivar)
            if(ndim>1) m%cruold(ind,2+(ivar-1)*ndim,igrid)= r%cr_fy_bound(ibound,ivar)
            if(ndim>2) m%cruold(ind,3+(ivar-1)*ndim,igrid)= r%cr_fz_bound(ibound,ivar)
        end do
     end do
#endif

  endif

  ! Imposed BC from condinit
  if(type == 4)then

     do ind=1,twotondim
        do idim=1,ndim
           nstride=2**(idim-1)
           xx(1,idim)=(2*m%grid(igrid)%ckey(idim)+MOD((ind-1)/nstride,2)+0.5)*dx-m%skip(idim)
        end do
        ! Call initial condition routine
        call condinit(r,g,xx,qq,dx,1)
        ! Compute kinetic and internal energy densities
        rr=qq(1,1)
        vx=qq(1,2)
        vy=qq(1,3)
        vz=qq(1,4)
        pp=qq(1,5)
        ekin=0.5d0*rr*(vx**2+vy**2+vz**2)
        eint=pp/(r%gamma-1.0)
        emag=0.0d0
        erad=0.0d0
#if NENER>0
        ! Compute non-thermal energy densities
        do irad=1,nener
           uu(1,5+irad)=qq(1,5+irad)/(r%gamma_rad(irad)-1.0d0)
           erad=erad+uu(1,5+irad)
        end do
#endif
#if NVAR>5+NENER
        ! Compute passive scalar density
        do ivar=6+nener,nvar
           uu(1,ivar)=rr*qq(1,ivar)
        enddo
#endif
        ! Compute density
        uu(1,1)=qq(1,1)
        ! Compute momentum density
        do idim=1,3
           uu(1,idim+1)=rr*qq(1,idim+1)
        end do
        ! Compute total fluid energy density
        uu(1,5)=eint+ekin+erad+emag
        ! Convert to conservative variables
        do ivar=1,nvar
           m%uold(ind,ivar,igrid)=uu(1,ivar)
        end do
     end do

  endif

  ! Imposed BC from boundana
  if(type == 5)then

     do ind=1,twotondim
        do idim=1,ndim
           nstride=2**(idim-1)
           xx(1,idim)=(2*m%grid(igrid)%ckey(idim)+MOD((ind-1)/nstride,2)+0.5)*dx-m%skip(idim)
        end do
        ! Call initial condition routine
        call boundana(r,g,xx,uu,dx,ibound,1)
        ! Scatter variables to main memory
        do ivar=1,nvar
           m%uold(ind,ivar,igrid)=uu(1,ivar)
        end do
     end do

#ifdef DO_CR
     do ind=1,twotondim
        do idim=1,ndim
           nstride=2**(idim-1)
           xx(1,idim)=(2*m%grid(igrid)%ckey(idim)+MOD((ind-1)/nstride,2)+0.5)*dx-m%skip(idim)
        end do
        ! Call initial condition routine
        call cr_boundana(r,g,xx,cruu,dx,ibound,1)
        ! Scatter variables to main memory
        do ivar=1,ncruvar
           m%cruold(ind,ivar,igrid)=cruu(1,ivar)
        end do
     end do
#endif

  endif

#endif

#ifdef GRAV

  ! Without the Poisson solver, boundary grids have no force
  if(.not.r%poisson)then
     m%f(1:twotondim,1:ndim,igrid)=0d0
     return
  endif

  do ind=1,twotondim
     do idim=1,ndim
        nstride=2**(idim-1)
        xx(1,idim)=(2*m%grid(igrid)%ckey(idim)+MOD((ind-1)/nstride,2)+0.5)*dx-m%skip(idim)
     end do
     ! Call analytical acceleration routine
     call gravana(r,g,xx,ff,dx,1)
     do idim=1,ndim
        m%f(ind,idim,igrid)=ff(1,idim)
     end do
     ! Call analytical potential routine
     call phiana(r,g,xx,phi,dx,1)
     m%phi(ind,igrid)=phi(1)
     m%phi_old(ind,igrid)=phi(1)
  end do

#endif

end subroutine init_bound_refine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_bound_phi(r,g,m,igrid,igrid_ref,ibound)
  use amr_parameters, only: ndim, twotondim, nvector
  use amr_commons, only: run_t, global_t, mesh_t, oct
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer::igrid
  integer::igrid_ref
  integer::ibound

  integer::idim, ind, nstride
  real(kind=8)::dx
  real(kind=8),dimension(1:nvector,1:ndim)::xx
  real(kind=8),dimension(1:nvector)::pp

#ifdef GRAV

  ! Mesh size at level ilevel in code units
  dx=r%boxlen/2**m%grid(igrid)%lev

  do ind=1,twotondim
     do idim=1,ndim
        nstride=2**(idim-1)
        xx(1,idim)=(2*m%grid(igrid)%ckey(idim)+MOD((ind-1)/nstride,2)+0.5)*dx-m%skip(idim)
     end do
     ! Call analytical potential routine
     call phiana(r,g,xx,pp,dx,1)
     m%phi(ind,igrid)=pp(1)
     m%phi_old(ind,igrid)=pp(1)
     ! Set mask to -1
     m%f(ind,3,igrid)=-1
  end do

#endif

end subroutine init_bound_phi
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_bound_mg(r,g,m,igrid,igrid_ref,ibound)
  use amr_parameters, only: ndim, twotondim, nvector
  use amr_commons, only: run_t, global_t, mesh_t, oct
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer::igrid
  integer::igrid_ref
  integer::ibound

  integer::ind

#ifdef GRAV
  do ind=1,twotondim
     m%phi(ind,igrid)=0.0
     m%f(ind,3,igrid)=-1
  end do
#endif

end subroutine init_bound_mg
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_bound_grav(r,g,m,igrid,igrid_ref,ibound)
  use amr_parameters, only: ndim, twotondim, nvector
  use amr_commons, only: run_t, global_t, mesh_t, oct
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer::igrid
  integer::igrid_ref
  integer::ibound

  integer::idim, ind, nstride
  real(kind=8)::dx
  real(kind=8),dimension(1:nvector,1:ndim)::xx
  real(kind=8),dimension(1:nvector,1:ndim)::ff

#ifdef GRAV

  ! Without the Poisson solver, boundary grids have no force
  if(.not.r%poisson)then
     m%f(1:twotondim,1:ndim,igrid)=0d0
     return
  endif

  ! Mesh size at level ilevel in code units
  dx=r%boxlen/2**m%grid(igrid)%lev

  do ind=1,twotondim
     do idim=1,ndim
        nstride=2**(idim-1)
        xx(1,idim)=(2*m%grid(igrid)%ckey(idim)+MOD((ind-1)/nstride,2)+0.5)*dx-m%skip(idim)
     end do
     ! Call analytical acceleration routine
     call gravana(r,g,xx,ff,dx,1)
     do idim=1,ndim
        m%f(ind,idim,igrid)=ff(1,idim)
     end do
  end do

#endif

end subroutine init_bound_grav
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_bound_flag(r,g,m,igrid,igrid_ref,ibound)
  use amr_parameters, only: ndim, twotondim
  use hydro_parameters, only: nvar, nener
  use amr_commons, only: run_t, global_t, mesh_t, oct
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer::igrid
  integer::igrid_ref
  integer::ibound

  integer,dimension(1:8,1:3)::ref_right=reshape(&
       & (/ 0,1,0,3,0,5,0,7,&
       &    0,0,1,2,0,0,5,6,&
       &    0,0,0,0,1,2,3,4/),(/8,3/))
  integer,dimension(1:8,1:3)::ref_left=reshape(&
       & (/ 2,0,4,0,6,0,8,0,&
       &    3,4,0,0,7,8,0,0,&
       &    5,6,7,8,0,0,0,0/),(/8,3/))

  integer::idim, ind
  integer::type, dir, shift

  type = r%bound_type(ibound)
  dir = r%bound_dir(ibound)
  shift = r%bound_shift(ibound)

  ! Set refinement map
  if(shift==+1)then
     do ind=1,twotondim
        if(ref_right(ind,dir)>0)then
           m%flag1(ind,igrid)=m%flag1(ref_right(ind,dir),igrid_ref)
        else
           m%flag1(ind,igrid)=0
        endif
     end do
  end if
  if(shift==-1)then
     do ind=1,twotondim
        if(ref_left(ind,dir)>0)then
           m%flag1(ind,igrid)=m%flag1(ref_left(ind,dir),igrid_ref)
        else
           m%flag1(ind,igrid)=0
        endif
     end do
  end if

end subroutine init_bound_flag
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_bound_flag2(r,g,m,igrid,igrid_ref,ibound)
  use amr_parameters, only: ndim, twotondim
  use hydro_parameters, only: nvar, nener
  use amr_commons, only: run_t, global_t, mesh_t, oct
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer::igrid
  integer::igrid_ref
  integer::ibound

  integer,dimension(1:8,1:3)::ref_right=reshape(&
       & (/ 0,1,0,3,0,5,0,7,&
       &    0,0,1,2,0,0,5,6,&
       &    0,0,0,0,1,2,3,4/),(/8,3/))
  integer,dimension(1:8,1:3)::ref_left=reshape(&
       & (/ 2,0,4,0,6,0,8,0,&
       &    3,4,0,0,7,8,0,0,&
       &    5,6,7,8,0,0,0,0/),(/8,3/))

  integer::idim, ind
  integer::type, dir, shift

  type = r%bound_type(ibound)
  dir = r%bound_dir(ibound)
  shift = r%bound_shift(ibound)

  ! Set refinement map
  if(shift==+1)then
     do ind=1,twotondim
        if(ref_right(ind,dir)>0)then
           m%flag2(ind,igrid)=m%flag2(ref_right(ind,dir),igrid_ref)
        else
          m%flag2(ind,igrid)=0
        endif
     end do
  end if
  if(shift==-1)then
     do ind=1,twotondim
        if(ref_left(ind,dir)>0)then
           m%flag2(ind,igrid)=m%flag2(ref_left(ind,dir),igrid_ref)
        else
           m%flag2(ind,igrid)=0
        endif
     end do
  end if

end subroutine init_bound_flag2
!################################################################
!################################################################
!################################################################
!################################################################
end module boundaries
