module input_hydro_condinit_module
contains
!###############################################
!###############################################
!###############################################
!###############################################
recursive subroutine r_input_hydro_condinit(pst,ilevel,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size

  integer::ilevel
  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_INPUT_HYDRO_CONDINIT,pst%iUpper+1,input_size,0,ilevel)
     call r_input_hydro_condinit(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call input_hydro_condinit(pst%s%r,pst%s%g,pst%s%m,ilevel)
  endif

end subroutine r_input_hydro_condinit
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine input_hydro_condinit(r,g,m,ilevel)
  use amr_parameters, only: ndim, twotondim, nvector
  use hydro_parameters, only: nvar, nener
  use amr_commons, only: run_t,global_t,mesh_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer::ilevel
  !------------------------------------------------------------
  ! Compute hydro primitive variables from subroutine condinit.
  ! Input magnetic field if MHD is activated.
  ! Finally convert primitive to conservative variables.
  !------------------------------------------------------------
  integer::igrid,ngrid,ind,idim,nstride,i,ivar,irad
  real(kind=8),dimension(1:nvector,1:ndim)::xx
#ifdef MHD
  real(kind=8),dimension(1:nvector,1:nvar+3-ndim)::qq
#else
  real(kind=8),dimension(1:nvector,1:nvar)::qq
#endif
  real(kind=8)::dx

  if(m%noct(ilevel)==0)return

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
        ! Scatter primitive variables to main memory
        do ivar=1,nvar
           do i=1,ngrid
              m%uold(ind,ivar,igrid+i-1)=qq(i,ivar)
           end do
        end do
#ifdef MHD
#if NDIM==1
        do i=1,ngrid
           m%bold(ind,1,igrid+i-1)=r%A_ave
           m%bold(ind,4,igrid+i-1)=r%A_ave
           m%bold(ind,2,igrid+i-1)=qq(i,nvar+1)
           m%bold(ind,5,igrid+i-1)=qq(i,nvar+1)
           m%bold(ind,3,igrid+i-1)=qq(i,nvar+2)
           m%bold(ind,6,igrid+i-1)=qq(i,nvar+2)
        end do
#endif
#if NDIM==2
        do i=1,ngrid
           m%bold(ind,3,igrid+i-1)=qq(i,nvar+1)
           m%bold(ind,6,igrid+i-1)=qq(i,nvar+1)
        end do
#endif
#endif
     end do
     ! End loop over cells
  end do
  ! End loop over grids
#endif

  ! Compute initial magnetic field
  call input_hydro_vecpot(r,g,m,ilevel)

  ! Convert primitive to conservative
  call cons_from_prim(r,g,m,ilevel)

end subroutine input_hydro_condinit
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine input_hydro_vecpot(r,g,m,ilevel)
  use amr_parameters, only: ndim, twotondim, nvector
  use amr_commons, only: run_t, global_t, mesh_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer::ilevel
  !------------------------------------------------------------------------------
  ! Compute magnetic field from vector potential from subroutine vecpotentialinit
  !------------------------------------------------------------------------------
  integer::igrid,ngrid,ind,i
  integer::l,nfine,ii,jj,kk
  real(kind=8),dimension(1:nvector,1:ndim)::xll,xrl,xlr,xrr
  real(kind=8),dimension(1:nvector)::ax,ay,az
  real(kind=8),dimension(1:nvector)::axll,axrl,axlr,axrr
  real(kind=8),dimension(1:nvector)::ayll,ayrl,aylr,ayrr
  real(kind=8),dimension(1:nvector)::azll,azrl,azlr,azrr
  real(kind=8)::dx,dxmin

  if(m%noct(ilevel)==0)return

#ifdef MHD
#if NDIM>1
  ! Mesh size at level ilevel in code units
  dx=r%boxlen/2**ilevel
  dxmin=r%boxlen/2**r%nlevelmax
  nfine=2**(r%nlevelmax-ilevel)
  ! Loop over grids by vector sweeps
  do igrid=m%head(ilevel),m%tail(ilevel),nvector
     ngrid=MIN(nvector,m%tail(ilevel)-igrid+1)
#if NDIM==2
     do ii=0,1 ! Loop over cells
        do jj=0,1
           ind=1+ii+2*jj
           do i=1,ngrid ! Compute 4 edges positions in code units
              xll(i,1)=(2*m%grid(igrid+i-1)%ckey(1)+ii  )*dx-m%skip(1)
              xll(i,2)=(2*m%grid(igrid+i-1)%ckey(2)+jj  )*dx-m%skip(2)
              xrl(i,1)=(2*m%grid(igrid+i-1)%ckey(1)+ii+1)*dx-m%skip(1)
              xrl(i,2)=(2*m%grid(igrid+i-1)%ckey(2)+jj  )*dx-m%skip(2)
              xlr(i,1)=(2*m%grid(igrid+i-1)%ckey(1)+ii  )*dx-m%skip(1)
              xlr(i,2)=(2*m%grid(igrid+i-1)%ckey(2)+jj+1)*dx-m%skip(2)
              xrr(i,1)=(2*m%grid(igrid+i-1)%ckey(1)+ii+1)*dx-m%skip(1)
              xrr(i,2)=(2*m%grid(igrid+i-1)%ckey(2)+jj+1)*dx-m%skip(2)
           end do
           call vecpotentialinit(r,g,xll,azll,3,ngrid)
           call vecpotentialinit(r,g,xrl,azrl,3,ngrid)
           call vecpotentialinit(r,g,xlr,azlr,3,ngrid)
           call vecpotentialinit(r,g,xrr,azrr,3,ngrid)
           do i=1,ngrid
              ! bx = d Az / dy
              m%bold(ind,1,igrid+i-1)=r%A_ave+(azlr(i)-azll(i))/dx
              m%bold(ind,4,igrid+i-1)=r%A_ave+(azrr(i)-azrl(i))/dx
              ! by = - d Az / dx
              m%bold(ind,2,igrid+i-1)=r%B_ave-(azrl(i)-azll(i))/dx
              m%bold(ind,5,igrid+i-1)=r%B_ave-(azrr(i)-azlr(i))/dx
           end do
        end do
     end do ! End loop over cells
#endif
#if NDIM==3
     do ii=0,1 ! Loop over cells
        do jj=0,1
           do kk=0,1
              ind=1+ii+2*jj+4*kk
              axll(1:ngrid)=0.0d0; axrl(1:ngrid)=0.0d0; axlr(1:ngrid)=0.0d0; axrr(1:ngrid)=0.0d0
              ayll(1:ngrid)=0.0d0; ayrl(1:ngrid)=0.0d0; aylr(1:ngrid)=0.0d0; ayrr(1:ngrid)=0.0d0
              azll(1:ngrid)=0.0d0; azrl(1:ngrid)=0.0d0; azlr(1:ngrid)=0.0d0; azrr(1:ngrid)=0.0d0
              ! Ax
              do i=1,ngrid ! Compute edge position in code units
                 xll(i,1)=(2*m%grid(igrid+i-1)%ckey(1)+ii)*dx-m%skip(1)+dxmin/2
                 xrl(i,1)=(2*m%grid(igrid+i-1)%ckey(1)+ii)*dx-m%skip(1)+dxmin/2
                 xlr(i,1)=(2*m%grid(igrid+i-1)%ckey(1)+ii)*dx-m%skip(1)+dxmin/2
                 xrr(i,1)=(2*m%grid(igrid+i-1)%ckey(1)+ii)*dx-m%skip(1)+dxmin/2
                 xll(i,2)=(2*m%grid(igrid+i-1)%ckey(2)+jj  )*dx-m%skip(2)
                 xll(i,3)=(2*m%grid(igrid+i-1)%ckey(3)+kk  )*dx-m%skip(3)
                 xrl(i,2)=(2*m%grid(igrid+i-1)%ckey(2)+jj+1)*dx-m%skip(2)
                 xrl(i,3)=(2*m%grid(igrid+i-1)%ckey(3)+kk  )*dx-m%skip(3)
                 xlr(i,2)=(2*m%grid(igrid+i-1)%ckey(2)+jj  )*dx-m%skip(2)
                 xlr(i,3)=(2*m%grid(igrid+i-1)%ckey(3)+kk+1)*dx-m%skip(3)
                 xrr(i,2)=(2*m%grid(igrid+i-1)%ckey(2)+jj+1)*dx-m%skip(2)
                 xrr(i,3)=(2*m%grid(igrid+i-1)%ckey(3)+kk+1)*dx-m%skip(3)
              end do
              do l=1,nfine
                 call vecpotentialinit(r,g,xll,ax,1,ngrid)
                 axll(1:ngrid)=axll(1:ngrid)+ax(1:ngrid)/dble(nfine)
                 call vecpotentialinit(r,g,xrl,ax,1,ngrid)
                 axrl(1:ngrid)=axrl(1:ngrid)+ax(1:ngrid)/dble(nfine)
                 call vecpotentialinit(r,g,xlr,ax,1,ngrid)
                 axlr(1:ngrid)=axlr(1:ngrid)+ax(1:ngrid)/dble(nfine)
                 call vecpotentialinit(r,g,xrr,ax,1,ngrid)
                 axrr(1:ngrid)=axrr(1:ngrid)+ax(1:ngrid)/dble(nfine)
                 xll(1:ngrid,1)=xll(1:ngrid,1)+dxmin
                 xrl(1:ngrid,1)=xrl(1:ngrid,1)+dxmin
                 xlr(1:ngrid,1)=xlr(1:ngrid,1)+dxmin
                 xrr(1:ngrid,1)=xrr(1:ngrid,1)+dxmin
              end do
              ! Ay
              do i=1,ngrid ! Compute edge position in code units
                 xll(i,2)=(2*m%grid(igrid+i-1)%ckey(2)+jj)*dx-m%skip(2)+dxmin/2
                 xrl(i,2)=(2*m%grid(igrid+i-1)%ckey(2)+jj)*dx-m%skip(2)+dxmin/2
                 xlr(i,2)=(2*m%grid(igrid+i-1)%ckey(2)+jj)*dx-m%skip(2)+dxmin/2
                 xrr(i,2)=(2*m%grid(igrid+i-1)%ckey(2)+jj)*dx-m%skip(2)+dxmin/2
                 xll(i,1)=(2*m%grid(igrid+i-1)%ckey(1)+ii  )*dx-m%skip(1)
                 xll(i,3)=(2*m%grid(igrid+i-1)%ckey(3)+kk  )*dx-m%skip(3)
                 xrl(i,1)=(2*m%grid(igrid+i-1)%ckey(1)+ii+1)*dx-m%skip(1)
                 xrl(i,3)=(2*m%grid(igrid+i-1)%ckey(3)+kk  )*dx-m%skip(3)
                 xlr(i,1)=(2*m%grid(igrid+i-1)%ckey(1)+ii  )*dx-m%skip(1)
                 xlr(i,3)=(2*m%grid(igrid+i-1)%ckey(3)+kk+1)*dx-m%skip(3)
                 xrr(i,1)=(2*m%grid(igrid+i-1)%ckey(1)+ii+1)*dx-m%skip(1)
                 xrr(i,3)=(2*m%grid(igrid+i-1)%ckey(3)+kk+1)*dx-m%skip(3)
              end do
              do l=1,nfine
                 call vecpotentialinit(r,g,xll,ay,2,ngrid)
                 ayll(1:ngrid)=ayll(1:ngrid)+ay(1:ngrid)/dble(nfine)
                 call vecpotentialinit(r,g,xrl,ay,2,ngrid)
                 ayrl(1:ngrid)=ayrl(1:ngrid)+ay(1:ngrid)/dble(nfine)
                 call vecpotentialinit(r,g,xlr,ay,2,ngrid)
                 aylr(1:ngrid)=aylr(1:ngrid)+ay(1:ngrid)/dble(nfine)
                 call vecpotentialinit(r,g,xrr,ay,2,ngrid)
                 ayrr(1:ngrid)=ayrr(1:ngrid)+ay(1:ngrid)/dble(nfine)
                 xll(1:ngrid,2)=xll(1:ngrid,2)+dxmin
                 xrl(1:ngrid,2)=xrl(1:ngrid,2)+dxmin
                 xlr(1:ngrid,2)=xlr(1:ngrid,2)+dxmin
                 xrr(1:ngrid,2)=xrr(1:ngrid,2)+dxmin
              end do
              ! Az
              do i=1,ngrid ! Compute edge position in code units
                 xll(i,3)=(2*m%grid(igrid+i-1)%ckey(3)+kk)*dx-m%skip(3)+dxmin/2
                 xrl(i,3)=(2*m%grid(igrid+i-1)%ckey(3)+kk)*dx-m%skip(3)+dxmin/2
                 xlr(i,3)=(2*m%grid(igrid+i-1)%ckey(3)+kk)*dx-m%skip(3)+dxmin/2
                 xrr(i,3)=(2*m%grid(igrid+i-1)%ckey(3)+kk)*dx-m%skip(3)+dxmin/2
                 xll(i,1)=(2*m%grid(igrid+i-1)%ckey(1)+ii  )*dx-m%skip(1)
                 xll(i,2)=(2*m%grid(igrid+i-1)%ckey(2)+jj  )*dx-m%skip(2)
                 xrl(i,1)=(2*m%grid(igrid+i-1)%ckey(1)+ii+1)*dx-m%skip(1)
                 xrl(i,2)=(2*m%grid(igrid+i-1)%ckey(2)+jj  )*dx-m%skip(2)
                 xlr(i,1)=(2*m%grid(igrid+i-1)%ckey(1)+ii  )*dx-m%skip(1)
                 xlr(i,2)=(2*m%grid(igrid+i-1)%ckey(2)+jj+1)*dx-m%skip(2)
                 xrr(i,1)=(2*m%grid(igrid+i-1)%ckey(1)+ii+1)*dx-m%skip(1)
                 xrr(i,2)=(2*m%grid(igrid+i-1)%ckey(2)+jj+1)*dx-m%skip(2)
              end do
              do l=1,nfine
                 call vecpotentialinit(r,g,xll,az,3,ngrid)
                 azll(1:ngrid)=azll(1:ngrid)+az(1:ngrid)/dble(nfine)
                 call vecpotentialinit(r,g,xrl,az,3,ngrid)
                 azrl(1:ngrid)=azrl(1:ngrid)+az(1:ngrid)/dble(nfine)
                 call vecpotentialinit(r,g,xlr,az,3,ngrid)
                 azlr(1:ngrid)=azlr(1:ngrid)+az(1:ngrid)/dble(nfine)
                 call vecpotentialinit(r,g,xrr,az,3,ngrid)
                 azrr(1:ngrid)=azrr(1:ngrid)+az(1:ngrid)/dble(nfine)
                 xll(1:ngrid,3)=xll(1:ngrid,3)+dxmin
                 xrl(1:ngrid,3)=xrl(1:ngrid,3)+dxmin
                 xlr(1:ngrid,3)=xlr(1:ngrid,3)+dxmin
                 xrr(1:ngrid,3)=xrr(1:ngrid,3)+dxmin
              end do
              do i=1,ngrid
                 ! bx = d Az / dy - d Ay /dz
                 m%bold(ind,1,igrid+i-1)=r%A_ave+(azlr(i)-azll(i)-(aylr(i)-ayll(i)))/dx
                 m%bold(ind,4,igrid+i-1)=r%A_ave+(azrr(i)-azrl(i)-(ayrr(i)-ayrl(i)))/dx
                 ! by = d Ax / dz - d Az /dx
                 m%bold(ind,2,igrid+i-1)=r%B_ave+(axlr(i)-axll(i)-(azrl(i)-azll(i)))/dx
                 m%bold(ind,5,igrid+i-1)=r%B_ave+(axrr(i)-axrl(i)-(azrr(i)-azlr(i)))/dx
                 ! bz = d Ay / dx - d Ax /dy
                 m%bold(ind,3,igrid+i-1)=r%C_ave+(ayrl(i)-ayll(i)-(axrl(i)-axll(i)))/dx
                 m%bold(ind,6,igrid+i-1)=r%C_ave+(ayrr(i)-aylr(i)-(axrr(i)-axlr(i)))/dx
              end do
           end do
        end do
     end do ! End loop over cells
#endif
  end do
  ! End loop over grids
#endif
#endif
end subroutine input_hydro_vecpot
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cons_from_prim(r,g,m,ilevel)
  use amr_parameters, only: ndim, twotondim, nvector
  use hydro_parameters, only: nvar, nener
  use amr_commons, only: run_t,global_t,mesh_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer::ilevel
  !--------------------------------------------
  ! Convert primitive to conservative variables
  !--------------------------------------------
  integer::igrid,ind,idim,i,ivar,irad
  real(kind=8)::rr,vx,vy,vz,pp
  real(kind=8)::bx,by,bz
  real(kind=8)::eint,ekin,emag,erad

  if(m%noct(ilevel)==0)return

#ifdef HYDRO
  ! Loop over grids
  do igrid=m%head(ilevel),m%tail(ilevel)
     ! Loop over cells
     do ind=1,twotondim

        ! Compute kinetic and internal energy densities
        rr=m%uold(ind,1,igrid)
        vx=m%uold(ind,2,igrid)
        vy=m%uold(ind,3,igrid)
        vz=m%uold(ind,4,igrid)
        pp=m%uold(ind,5,igrid)
        ekin=0.5d0*rr*(vx**2+vy**2+vz**2)
        eint=pp/(r%gamma-1.0)
        emag=0.0d0
#ifdef MHD
        ! Compute magnetic energy for all cells
        bx=0.5d0*(m%bold(ind,1,igrid)+m%bold(ind,4,igrid))
        by=0.5d0*(m%bold(ind,2,igrid)+m%bold(ind,5,igrid))
        bz=0.5d0*(m%bold(ind,3,igrid)+m%bold(ind,6,igrid))
        emag=0.5d0*(bx**2+by**2+bz**2)
#endif
        erad=0.0d0
#if NENER>0
        ! Compute non-thermal energy densities
        do irad=1,nener
           m%uold(ind,5+irad,igrid)=m%uold(ind,5+irad,igrid)/(r%gamma_rad(irad)-1.0d0)
           erad=erad+m%uold(ind,5+irad,igrid)
        end do
#endif
        ! Compute total fluid energy density
        m%uold(ind,5,igrid)=eint+ekin+erad+emag
        ! Compute momentum density
        do idim=1,3
           m%uold(ind,idim+1,igrid)=rr*m%uold(ind,idim+1,igrid)
        end do
#if NVAR>5+NENER
        ! Compute passive scalar density
        do ivar=6+nener,nvar
           m%uold(ind,ivar,igrid)=rr*m%uold(ind,ivar,igrid)
        enddo
#endif
     end do
     ! End loop over cells
  end do
  ! End loop over grids
#endif

end subroutine cons_from_prim
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine prim_from_cons(r,g,m,ilevel)
  use amr_parameters, only: ndim, twotondim, nvector
  use hydro_parameters, only: nvar, nener
  use amr_commons, only: run_t, global_t, mesh_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  integer::ilevel
  !--------------------------------------------
  ! Convert primitive to conservative variables
  !--------------------------------------------
  integer::igrid,ind,idim,i,ivar,irad
  real(kind=8)::rr,vx,vy,vz,pp
  real(kind=8)::bx,by,bz
  real(kind=8)::eint,ekin,emag,erad,etot

  if(m%noct(ilevel)==0)return

#ifdef HYDRO
  ! Loop over grids
  do igrid=m%head(ilevel),m%tail(ilevel)
     ! Loop over cells
     do ind=1,twotondim
        ! Compute velocities and kinetic energy density
        rr=m%uold(ind,1,igrid)
        vx=m%uold(ind,2,igrid)/rr
        vy=m%uold(ind,3,igrid)/rr
        vz=m%uold(ind,4,igrid)/rr
        ekin=0.5d0*rr*(vx**2+vy**2+vz**2)
        m%uold(ind,2,igrid)=vx
        m%uold(ind,3,igrid)=vy
        m%uold(ind,4,igrid)=vz
        emag=0.0d0
#ifdef MHD
        ! Compute magnetic energy for all cells
        bx=0.5d0*(m%bold(ind,1,igrid)+m%bold(ind,4,igrid))
        by=0.5d0*(m%bold(ind,2,igrid)+m%bold(ind,5,igrid))
        bz=0.5d0*(m%bold(ind,3,igrid)+m%bold(ind,6,igrid))
        emag=0.5d0*(bx**2+by**2+bz**2)
#endif
        erad=0.0d0
#if NENER>0
        ! Compute non-thermal pressures
        do irad=1,nener
           erad=erad+m%uold(ind,5+irad,igrid)
           m%uold(ind,5+irad,igrid)=m%uold(ind,5+irad,igrid)*(r%gamma_rad(irad)-1.0d0)
        end do
#endif
        ! Compute internal energy
        etot=m%uold(ind,5,igrid)
        eint=etot-ekin-emag-erad
        ! Compute thermal pressure
        pp=(r%gamma-1.0)*eint
        m%uold(ind,5,igrid)=pp
#if NVAR>5+NENER
        ! Compute passive scalar mass fraction
        do ivar=6+nener,nvar
           m%uold(ind,ivar,igrid)=m%uold(ind,ivar,igrid)/rr
        enddo
#endif
     end do
     ! End loop over cells
  end do
  ! End loop over grids
#endif

end subroutine prim_from_cons
!################################################################
!################################################################
!################################################################
!################################################################
subroutine region_condinit(r,g,x,q,dx,nn)
  use amr_parameters, only: nvector, ndim
  use hydro_parameters, only: nvar, nener
  use amr_commons, only: run_t, global_t
#ifdef DO_RTZ
  use rtz_module, only: elements, n_elements
#endif
  implicit none
  type(run_t)::r
  type(global_t)::g
  integer ::nn
  real(kind=8)::dx
#ifdef MHD
  real(kind=8),dimension(1:nvector,1:nvar+3-ndim)::q
#else
  real(kind=8),dimension(1:nvector,1:nvar)::q
#endif
  real(kind=8),dimension(1:nvector,1:ndim)::x
  !----------------------------------------------------
  ! This routine sets simple pre-defined initial
  ! conditions, like points, squares, etc.
  !----------------------------------------------------
#if NVAR>5
  integer::ivar
#endif
  integer::i,k
  real(kind=8)::vol,rad,weight,xn,yn,zn,en
#ifdef DO_RTZ
  integer::counter
  real(kind=8)::total_mass_density
#endif
  ! Set some (tiny) default values in case n_region=0
  q(1:nn,1)=r%smallr
  q(1:nn,2)=0.0d0
  q(1:nn,3)=0.0d0
  q(1:nn,4)=0.0d0
  q(1:nn,5)=r%smallr*r%smallc**2/r%gamma
#if NVAR>5
  do ivar=6,nvar
     q(1:nn,ivar)=0.0d0
  enddo
#endif
#ifdef MHD
#if NDIM<3
  q(1:nn,nvar+1)=0.0d0
#if NDIM==1
  q(1:nn,nvar+2)=0.0d0
#endif
#endif
#endif

  ! Loop over initial conditions regions
  do k=1,r%nregion

     ! For "square" regions only:
     if(r%region_type(k) .eq. 'square')then
        ! Exponent of choosen norm
        en=r%exp_region(k)
        do i=1,nn
           ! Compute position in normalized coordinates
           xn=0.0d0; yn=0.0d0; zn=0.0d0
           xn=2.0d0*abs(x(i,1)-r%x_center(k))/r%length_x(k)
#if NDIM>1
           yn=2.0d0*abs(x(i,2)-r%y_center(k))/r%length_y(k)
#endif
#if NDIM>2
           zn=2.0d0*abs(x(i,3)-r%z_center(k))/r%length_z(k)
#endif
           ! Compute cell "radius" relative to region center
           if(r%exp_region(k)<10)then
              rad=(xn**en+yn**en+zn**en)**(1.0/en)
           else
              rad=max(xn,yn,zn)
           end if
           ! If cell lies within region,
           ! REPLACE primitive variables by region values
           if(rad<1.0)then
              q(i,1)=r%d_region(k)
              q(i,2)=r%u_region(k)
              q(i,3)=r%v_region(k)
              q(i,4)=r%w_region(k)
              q(i,5)=r%p_region(k)
#if NENER>0
              do ivar=6,5+nener
                 q(i,ivar)=r%prad_region(k,ivar-5)
              enddo
#endif
#if NVAR>5+NENER
              do ivar=6+nener,nvar
                 q(i,ivar)=r%var_region(k,ivar-5-nener)
              end do
#ifdef DO_RTZ
              counter = 0
              total_mass_density = 0.d0
              do ivar=1,n_elements ! loop over elements
                 if (elements(ivar)%atomic_number.gt.0) then
                  !   write(*,*) g%myid,elements(ivar)%atomic_number,elements(ivar)%z_solar,r%z_ave
                    ! This gives us the total mass density
                    q(i,r%ichem+counter) = elements(ivar)%z_solar * r%z_ave * elements(ivar)%atomic_mass
                    total_mass_density = total_mass_density + q(i,r%ichem+counter)
                    counter = counter + 1
                 end if
              end do ! end loop over elements

              ! Loop again and renormalize
              counter = 0
              do ivar=1,n_elements ! loop over elements
                 if (elements(ivar)%atomic_number.gt.0) then
                    ! This gives us the total mass density
                    q(i,r%ichem+counter) = (1.d0 / total_mass_density) * q(i,r%ichem+counter) ! HK note will be multiplied by rho later
                    counter = counter + 1
                 end if
              end do ! end loop over elements
#endif
#endif
#ifdef MHD
#if NDIM==1
              q(i,nvar+1)=r%B_region(k)
              q(i,nvar+2)=r%C_region(k)
#endif
#if NDIM==2
              q(i,nvar+1)=r%C_region(k)
#endif
#endif
           end if
        end do
     end if

     ! For "point" regions only:
     if(r%region_type(k) .eq. 'point')then
        ! Volume elements
        vol=dx**ndim
        ! Compute CIC weights relative to region center
        do i=1,nn
           xn=1.0; yn=1.0; zn=1.0
           xn=max(1.0-abs(x(i,1)-r%x_center(k))/dx,0.0d0)
#if NDIM>1
           yn=max(1.0-abs(x(i,2)-r%y_center(k))/dx,0.0d0)
#endif
#if NDIM>2
           zn=max(1.0-abs(x(i,3)-r%z_center(k))/dx,0.0d0)
#endif
           weight=xn*yn*zn
           ! If cell lies within CIC cloud,
           ! ADD to primitive variables the region values
           q(i,1)=q(i,1)+r%d_region(k)*weight/vol
           q(i,2)=q(i,2)+r%u_region(k)*weight
           q(i,3)=q(i,3)+r%v_region(k)*weight
           q(i,4)=q(i,4)+r%w_region(k)*weight
           q(i,5)=q(i,5)+r%p_region(k)*weight/vol
#if NENER>0
           do ivar=6,5+nener
              q(i,ivar)=q(i,ivar)+r%prad_region(k,ivar-5)*weight/vol
           enddo
#endif
#if NVAR>5+NENER
           do ivar=6+nener,nvar
              q(i,ivar)=r%var_region(k,ivar-5-nener)
           end do
#endif
        end do
     end if
  end do

  return
end subroutine region_condinit
!################################################################
!################################################################
!################################################################
!################################################################
end module input_hydro_condinit_module
