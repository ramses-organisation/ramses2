module init_xion_module
contains
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_init_xion(pst)

! Initialize hydrogen ionization state in all cells at all levels from
! density and temperature in the cells, assuming chemical equilibrium.
!-------------------------------------------------------------------------
  use ramses_commons, only: pst_t
  use upload_module, only: m_upload_fine
  implicit none
  type(pst_t)::pst
  integer::ilevel

  if(pst%s%g%myid==1)then
     if(pst%s%r%cosmo .and. pst%s%r%nrestart==0)then
        write(*,*) 'Initialising to residual (freeze-out) ionisation fractions'
     else
        write(*,*) 'Initialising to CIE ionisation fractions'
     endif
  endif
  do ilevel=pst%s%r%nlevelmax,pst%s%r%levelmin,-1
      call r_init_xion(pst,ilevel,1)
      call m_upload_fine(pst,ilevel)
  end do

end subroutine m_init_xion

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
recursive subroutine r_init_xion(pst,ilevel,input_size)

! Initialize hydrogen ionization state in all cells at given level from
! density and temperature in the cells, assuming chemical equilibrium.
!-------------------------------------------------------------------------
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::ilevel

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request( &
             pst%s%mdl,MDL_INIT_XION,pst%iUpper+1,input_size,0,ilevel)
     call r_init_xion(pst%pLower,ilevel,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     call init_xion(pst%s%r,pst%s%g,pst%s%m,pst%s%tables,ilevel)
  endif

end subroutine r_init_xion
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine init_xion(r,g,m,tables,ilevel)
  use constants
  use amr_parameters, only:ndim,nvector,twotondim
  use hydro_parameters, only: nener, nion
  use amr_commons, only:run_t,global_t,mesh_t
  use cooling_module, only:cooling_t
#ifdef DO_RTZ
  use rtz_module, only: elements, n_elements
#else
  use neq_cooling_module, only: cmp_equilibrium_abundances
#endif
  use coolrates_module, only: neq_cooling_t
  implicit none
  type(run_t)::r
  type(global_t)::g
  type(mesh_t)::m
  type(neq_cooling_t)::tables
  integer::ilevel
  integer::i,ind,igrid,idim,ngrid,nleaf
  real(kind=8)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(kind=8)::mu,x
  integer,dimension(1:nvector)::ind_leaf
  real(kind=8),dimension(1:nvector)::nH,T2,ekk,err,emag,Zsolar
  real(kind=8),dimension(nion)::phI_rates       ! Photoionization rates [s-1]
  real(kind=8),dimension(1:nvector, 7)::nSpec    !          Species abundances
#if NENER>0
  integer::irad
#endif
#ifdef DO_RTZ
  integer::counter, iE, iI
#endif

  if(r%verbose.and.g%myid==1)write(*,'("   Entering init_xion for level ",I2)')ilevel

  associate(ixHI=>r%ixHi, ixHII=>r%ixHII, ixHeII=>r%ixHeII, ixHeIII=>r%ixHeIII)

  ! Conversion factor from user units to cgs units
  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  pHI_rates(:)=0.0                   ! No UV background for the time being

#ifdef HYDRO
  ! Loop over cells
  do ind=1,twotondim
     ! Loop over octs with vector sweeps
     do igrid=m%head(ilevel),m%tail(ilevel),nvector

        ! Collect a vector of leaf cells
        ngrid=MIN(nvector,m%tail(ilevel)-igrid+1)
        nleaf=0
        do i=1,ngrid
           if(.NOT. m%grid(igrid+i-1)%refined(ind))then
              nleaf=nleaf+1
              ind_leaf(nleaf)=igrid+i-1
           endif
        end do
        if(nleaf.eq.0)cycle

        ! Compute rho
        do i=1,nleaf
           nH(i)=MAX(dble(m%uold(ind,1,ind_leaf(i))),r%smallr)
        end do

        ! Compute metallicity in solar units
        if(r%metal)then
           do i=1,nleaf
              Zsolar(i)=m%uold(ind,r%imetal,ind_leaf(i))/nH(i)/0.02d0
           end do
        else
           do i=1,nleaf
              Zsolar(i)=r%z_ave
           end do
        endif

        do i=1,nleaf ! Total energy
           T2(i)=m%uold(ind,5,ind_leaf(i))
        end do

        do i=1,nleaf ! Kinetic energy
           ekk(i)=0.0d0
        end do
        do idim=1,3
           do i=1,nleaf
              ekk(i)=ekk(i)+0.5d0*m%uold(ind,idim+1,ind_leaf(i))**2/nH(i)
           end do
        end do

        do i=1,nleaf ! Non-thermal energies
           err(i)=0.0d0
        end do
#if NENER>0
        do irad=1,nener
           do i=1,nleaf
              err(i)=err(i)+m%uold(ind,5+irad,ind_leaf(i))
           end do
        end do
#endif
        do i=1,nleaf ! Magnetic energy
           emag(i)=0.0d0
        end do
#ifdef MHD
        do idim=1,3
           do i=1,nleaf
              emag(i)=emag(i)+0.125d0*(m%bold(ind,idim  ,ind_leaf(i)) &
                   &                  +m%bold(ind,idim+3,ind_leaf(i)))**2
           end do
        end do
#endif
        do i=1,nleaf ! Gas thermal pressure
           T2(i)=(r%gamma-1.0d0)*(T2(i)-ekk(i)-err(i)-emag(i))
        end do

        do i=1,nleaf ! Compute T2=T/mu in Kelvin
           T2(i)=T2(i)/nH(i)*scale_T2
        end do

        do i=1,nleaf ! Compute nH in H/cc
           nH(i)=nH(i)*scale_nH
        end do

        ! Do the main computation of equilibrium abundances
#ifdef DO_RTZ
        ! In the case of RTZ, we set everything to neutral to start,
        ! except for the residual HII fraction in cosmological runs
        do i=1,nleaf ! loop over cells
          counter = 1
          do iE=1,n_elements ! loop over elements
             if (elements(iE)%atomic_number.gt.0) then
                do iI=1,elements(iE)%n_ions ! loop over ions
                   x = 0.d0
                   if (iI.eq.1) then
                      x = 1.d0
                   end if
                   if (r%cosmo .and. r%nrestart==0 .and. iE.eq.1) then
                      if (iI.eq.1) x = 1.d0 - g%xHII_start
                      if (iI.eq.2) x = g%xHII_start
                   end if
                   m%uold(ind,r%iIons-1+counter,ind_leaf(i)) = x*m%uold(ind,1,ind_leaf(i))
                   counter = counter + 1
                end do  ! end loop over ions
             end if
          end do ! end loop over elements
          ! Molecules come after all ions (same layout as read_params)
          if (elements(1)%atomic_number.gt.0 .and. r%isH2_rtz) then
             m%uold(ind,r%iIons-1+counter,ind_leaf(i)) = 0.d0
             counter = counter + 1
          end if
        end do ! end loop over cells
#else
        if(r%cosmo .and. r%nrestart==0)then
           ! Use the abundances from the single-cell evolution in neq_set_model
           ! (CIE at T2_start would wipe out the residual ionization)
           do i=1,nleaf
              nSpec(i,:)=0d0
              if(r%isH2) then
                 nSpec(i,2)=0.5d0*(1d0-g%xion_start(ixHI)-g%xion_start(ixHII))
                 nSpec(i,3)=g%xion_start(ixHI)
              else
                 nSpec(i,3)=1d0-g%xion_start(ixHII)
              endif
              nSpec(i,4)=g%xion_start(ixHII)
              if(r%isHe)then
                 nSpec(i,5)=1d0-g%xion_start(ixHeII)-g%xion_start(ixHeIII)
                 nSpec(i,6)=g%xion_start(ixHeII)
                 nSpec(i,7)=g%xion_start(ixHeIII)
              endif
              nSpec(i,:)=nSpec(i,:)*nH(i)
           end do
        else
           do i=1,nleaf
              call cmp_equilibrium_abundances(r,tables,T2(i),nH(i),pHI_rates(:),mu,nSpec(i,:),Zsolar(i))
           end do
        endif

        ! nspec: 1=ne, 2=nH2, 3=nHI, 4=nHII, 5=nHeI, 6=nHeII, 7=n_HEIII
        do i=1,nleaf ! Update ionization states in uold

          if(r%isH2) then
              x = nSpec(i,3)/(2.*nSpec(i,2)+nSpec(i,3)+nspec(i,4))   ! HI fraction
              m%uold(ind,r%iIons-1+ixHI,ind_leaf(i)) = x*m%uold(ind,1,ind_leaf(i))
          endif

          x = nSpec(i,4)/(2.*nSpec(i,2)+nSpec(i,3)+nspec(i,4))        ! HII fraction
          m%uold(ind,r%iIons-1+ixHII,ind_leaf(i)) = x*m%uold(ind,1,ind_leaf(i))

          if(r%Y_He .gt. 0d0 .and. r%isHe) then
              x = nSpec(i,6)/(nSpec(i,5)+nSpec(i,6)+nSpec(i,7)) !  HeII fraction
              m%uold(ind,r%iIons-1+ixHeII,ind_leaf(i)) = x*m%uold(ind,1,ind_leaf(i))
              x = nSpec(i,7)/(nSpec(i,5)+nSpec(i,6)+nSpec(i,7)) ! HeIII fraction
              m%uold(ind,r%iIons-1+ixHeIII,ind_leaf(i)) = x*m%uold(ind,1,ind_leaf(i))
          endif
        end do
#endif

     end do
     ! End loop over grid
  end do
  ! End loop over cells
#endif
  end associate

end subroutine init_xion
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
#ifndef DO_RTZ
subroutine calc_equilibrium_xion(s, uold, &
#ifdef MHD
     & bold, &
#endif
#ifdef DO_RT
     & rtuold, &
#endif
     & ilevel, xion)

! Calculate and return photoionization equilibrium abundance states for
! a cell
! gridp     => Pointer to (type oct) grid containing the cell in question
! icell     => Index for cell in grid
! ilevel    => Cell refinement level
! xion      <= Returned equilibrium ionization fractions of cell
!-------------------------------------------------------------------------
  use amr_commons, only: oct
  use amr_parameters, only: ndim
  use hydro_parameters, only: nvar, nener, nion
  use neq_cooling_module, only: cmp_equilibrium_abundances
  use ramses_commons, only: ramses_t
  use rt_parameters, only: nrtvar, nrtgrp
  implicit none
  type(ramses_t)::s
  real(kind=8),dimension(1:nvar)::uold
#ifdef MHD
  real(kind=8),dimension(1:6)::bold
#endif
#ifdef DO_RT
  real(kind=8),dimension(1:nrtvar)::rtuold
#endif
  integer::ilevel
  real(kind=8),dimension(nion)::xion
!-------------------------------------------------------------------------
  integer::ip, iI, idim, iNp
  real(kind=8)::scale_nH, scale_T2, scale_l, scale_d, scale_t, scale_v
  real(kind=8)::scale_Np,scale_Fp,nH,T2,ekk,err,emag,mu,Zsolar,ss_factor
  real(kind=8),dimension(nion)::phI_rates    ! Photoionization rates [s-1]
  real(kind=8),dimension(7)::ns              !          Species abundances
#if NENER>0
  integer::irad
#endif
!-------------------------------------------------------------------------
#ifdef HYDRO
  ! Conversion factor from user units to cgs units
  call units(s%r,s%g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
#ifdef DO_RT
  call rt_units(s%r,s%g,scale_Np,scale_Fp)
#endif
  ! Calculate photoionization rates:
  phI_rates(:)=0.0
  do ip=1, nrtgrp
     iNp=1+(ip-1)*(ndim+1)
#ifdef DO_RT
     do iI=1,nion
        phI_rates(iI) = phI_rates(iI) &
                      + rtuold(iNp) &
                      * scale_Np*s%tables%signc(ip,iI,ilevel)
     end do
#endif
  end do

  nH = MAX(dble(uold(1)),s%r%smallr)      !   Nb density of gas [UU]

  Zsolar = s%r%z_ave
  if(s%r%metal) Zsolar=uold(s%r%imetal)/nH/0.02  ! Z (Solar U)

  ! Compute temperature from energy density
  T2 = uold(5)               ! Energy density (kin+heat) [UU]
  ekk = 0.0d0                !            Kinetic energy [UU]
  do idim=1,3
     ekk=ekk+0.5d0*uold(1+idim)**2/nH
  end do
  err = 0.0d0
#if NENER>0
  do irad=1,nener
     err=err+uold(5+irad)
  end do
#endif
  emag=0.0d0
#ifdef MHD
  do idim=1,3
     emag=emag+0.125d0*(bold(idim)+bold(idim+3))**2
  end do
#endif
  ! Gas thermal pressure
  T2=(s%r%gamma-1.0d0)*(T2-ekk-err-emag)
  T2 = T2/nH*scale_T2                       !                T/mu [Kelvin]
  nH = nH*scale_nH                          !        Number density [H/cc]

  ! UV background photoionization
  ss_factor = 1d0
  if(s%r%self_shielding) ss_factor = exp(-nH/1d-2)
  if(s%r%haardt_madau) &
     phI_rates = phI_rates + s%tables%UVrates(:,1) * ss_factor

  call cmp_Equilibrium_Abundances(s%r, s%tables, T2, nH, pHI_rates, mu   &
                                 ,ns, Zsolar)

  if(s%r%isH2) xion(s%r%ixHI)=ns(3)/(2.*ns(2)+ns(3)+ns(4)) !   HI fraction
  xion(s%r%ixHII)=ns(4)/(2.*ns(2)+ns(3)+ns(4))             !  HII fraction
  if(s%r%Y_He .gt. 0d0 .and. s%r%isHe) then
     xion(s%r%ixHeII) = ns(6)/(ns(5)+ns(6)+ns(7))          ! HeII fraction
     xion(s%r%ixHeIII) = ns(7)/(ns(5)+ns(6)+ns(7))         !HeIII fraction
  endif
#endif
end subroutine calc_equilibrium_xion
#endif
end module init_xion_module
