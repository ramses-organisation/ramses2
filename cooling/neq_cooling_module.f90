! Non-equlibrium (in H2, HI, HII, HeI, HeII, HeIII)
! cooling module for radiation-hydrodynamics.
! For details, see Rosdahl et al. (2013), Rosdahl & Teyssier (2015),
! and Nickerson, Teyssier, & Rosdahl (2018).
! Joki Rosdahl, Sarah Nickerson, Andreas Bleuler, and Romain Teyssier.
! NOTE: T2=T/mu, Np = photon density, Fp = photon flux,
module neq_cooling_module
  use amr_parameters, only: ndim, nvector
  use amr_commons, only: run_t
  use hydro_parameters, only: nion
  use rt_parameters
  use coolrates_module
  use constants
  implicit none

  private   ! default
  public neq_set_model, neq_solve_cooling, updateRTGroups_CoolConstants &
       ,update_metal_cooling, cmp_equilibrium_abundances, getmu, HsurH0

  real(kind=8),parameter::T2_min_fix=1d-2 ! Min temperature [K]

  ! cosmic ray ionisation rates, primary and secondary
  ! see Nickerson, Teyssier, & Rosdahl (2018)
  real(kind=8),parameter::cosray_H2 = 7.525d-16 ! Indriolo 2012, Gong 2017[s-1]
  real(kind=8),parameter::cosray_HI = 4.45d-16  ! Indriolo 2015, Gong 2017[s-1]
  ! for HeI ionisation, use Glover 2010 cosray_HeI = 1.1 * cosray_HI

  real(kind=8),parameter::T_min=0.1, T_frac=0.1
  real(kind=8),parameter::x_min=1d-20, x_fm=1d-6, x_frac=0.1
  real(kind=8),parameter::Np_min=1d-13, Np_frac=0.2
  real(kind=8),parameter::Fp_frac=0.5

  ! IR group index
  integer,parameter::iIR=1

CONTAINS

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE neq_set_model(r, tables, h, omegab, omega0, omegaL, astart_sim, T2_sim, &
     & xion_sim)
  ! Initialize cooling. All these parameters are unused at the moment and
  ! are only there for the original cooling-module.
  ! h (dble)            => H0/100
  ! omegab (dble)       => Omega Baryons
  ! omega0 (dble)       => Omega Matter total
  ! omegaL (dble)       => Omega Lambda
  ! astart_sim (dble)   => Redshift at which we start the simulation
  ! T2_sim (dble)      <=  Starting temperature in simulation?
  ! xion_sim (dble)    <=  Starting ionization fractions in simulation
  !-------------------------------------------------------------------------
!  use UV_module
  type(run_t) :: r
  type(neq_cooling_t) :: tables
  real(kind=8) :: h, omegab, omega0, omegaL, astart_sim, T2_sim
  real(kind=8),dimension(1:nion) :: xion_sim
  real(kind=8) :: astart=0.0001, aend, dasura, T2end=T2_min_fix, mu=1., ne, Ha
  !-------------------------------------------------------------------------

  ! Calculate initial temperature
  if (astart_sim < astart) then
     write(*,*) 'ERROR in set_model : astart_sim is too small.'
     write(*,*) 'astart     =',astart
     write(*,*) 'astart_sim =',astart_sim
     STOP
  endif
  aend=astart_sim
  dasura=0.02d0
  Ha = 1d7 * h / Mpc2cm                                          &
       * HsurH0(1d0/astart-1d0, omega0, omegaL, 1d0-omega0-omegaL)

!  call init_UV_background
  call init_coolrates_tables(r, tables, astart, Ha)
  call updateRTGroups_CoolConstants(r, tables)

  xion_sim=0d0
  if(r%nrestart==0 .and. r%cosmo)                                        &
      call neq_evol_single_cell(r, tables, astart, aend, dasura,         &
       &        h, omegab, omega0, omegaL, T2end, mu, ne, xion_sim, .false.)

  ! Set the a-dependent rates to astart_sim
  Ha = 1d7 * h / Mpc2cm                                          &
       * HsurH0(1d0/astart_sim-1d0, omega0, omegaL, 1d0-omega0-omegaL)
  call update_coolrates_tables(r, tables, astart_sim, Ha)
  T2_sim=T2end

END SUBROUTINE neq_set_model

!!$!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!!$SUBROUTINE update_UVrates(aexp)
!!$! Set the UV ionization and heating rates according to the given aexp.
!!$!-------------------------------------------------------------------------
!!$  use UV_module
!!$  use amr_parameters,only:haardt_madau
!!$  real(kind=8)::aexp
!!$!------------------------------------------------------------------------
!!$  UVrates=0.
!!$  if(.not. haardt_madau) RETURN
!!$
!!$  call inp_UV_rates_table(1./aexp - 1., UVrates, .true.)
!!$
!!$END SUBROUTINE update_UVrates

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE neq_solve_cooling(r, tables, T2, xion, nH, Zsolar, &
#ifdef DO_RT
     & Np, Fp, p_gas, dNpdt, dFpdt, ilevel, &
#endif
     & dt, nCell)
  ! Semi-implicitly solve for new temperature, ionization states,
  ! photon density/flux, and gas velocity in a number of cells.
  ! Parameters:
  ! tables  => object containing all non equilibrium chemistry tables
  ! T2     <=> T/mu [K]
  ! xion   <=> nion ionization fractions
  ! Np     <=> nrtgrp photon number densities [cm-3]
  ! Fp     <=> nrtgrp * ndim photon number fluxes [cm-2 s-1]
  ! p_gas  <=> ndim gas momentum densities [cm s-1 g cm-3]
  ! dNpdt   =>  Op split increment in photon densities during dt
  ! dFpdt   =>  Op split increment in photon flux magnitudes during dt
  ! nH      =>  Hydrogen number densities [cm-3]
  ! c_switch=>  Cooling switch (1 for cool/heat, 0 for no cool/heat) (OFF)
  ! Zsolar  =>  Cell metallicities [solar fraction]
  ! dt      =>  Timestep size             [s]
  ! nCell   =>  Number of cells (length of all the above vectors)
  ! ilevel  =>  Refinement levels
  !
  ! We use a slightly modified method of Anninos et al. (1997).
  !-------------------------------------------------------------------------
  implicit none
  type(run_t):: r
  type(neq_cooling_t):: tables
  real(kind=8),dimension(1:nvector):: T2
  real(kind=8),dimension(1:nion, 1:nvector):: xion
  real(kind=8),dimension(1:nvector):: nH, Zsolar
#ifdef DO_RT
  real(kind=8),dimension(1:ndim, 1:nvector):: p_gas
  real(kind=8),dimension(1:nrtgrp, 1:nvector):: Np, dNpdt
  real(kind=8),dimension(1:ndim, 1:nrtgrp, 1:nvector):: Fp, dFpdt
  integer::ilevel
#endif
!  logical,dimension(1:nvector):: c_switch
  real(kind=8)::dt
  integer::ncell
  !--------------------------------------------------------
  real(kind=8),dimension(1:nvector):: tLeft, ddt
  logical:: dt_ok
  real(kind=8):: dt_rec
  real(kind=8):: dT2
  real(kind=8),dimension(nion):: dXion
  integer::i, ia, nAct, nAct_next, loopcnt, code
  integer,dimension(1:nvector):: indAct              ! Active cell indexes
  real(kind=8):: one_over_x_FRAC, one_over_T_FRAC
#ifdef DO_RT
  integer::ig
  real(kind=8):: one_over_rt_c_cgs, one_over_egy_IR_erg
  real(kind=8):: one_over_Np_FRAC, one_over_Fp_FRAC
  real(kind=8),dimension(1:ndim):: dp_gas
  real(kind=8),dimension(nrtgrp):: dNp
  real(kind=8),dimension(1:ndim, 1:nrtgrp):: dFp
  real(kind=8),dimension(1:nrtgrp):: group_egy_ratio, group_egy_erg
#endif
  integer*8,dimension(20)::loopCodes=0

  associate(ixHI=>r%ixHi, ixHII=>r%ixHII, ixHeII=>r%ixHeII, ixHeIII=>r%ixHeIII)

  ! Store some temporary variables reduce computations
  one_over_T_FRAC = 1d0 / T_FRAC
  one_over_x_FRAC = 1d0 / x_FRAC
#ifdef DO_RT
  one_over_Np_FRAC = 1d0 / Np_FRAC
  one_over_Fp_FRAC = 1d0 / Fp_FRAC
  one_over_rt_c_cgs = 1d0 / tables%rt_c_cgs(ilevel)
  group_egy_erg(1:nrtgrp) = r%group_egy(1:nrtgrp) * eV2erg
  if(r%rt_isIR) then
     group_egy_ratio(1:nrtgrp) = r%group_egy(1:nrtgrp) / r%group_egy(iIR)
     one_over_egy_IR_erg = 1d0 / group_egy_erg(iIR)
  endif
#endif
  !-----------------------------------------------------------------------
  tleft(1:ncell) = dt                !       Time left in dt for each cell
  ddt(1:ncell) = dt                  ! First guess at sub-timestep lengths

  do i=1,ncell
     indact(i) = i                   !      Set up indexes of active cells
     ! Ensure all state vars are legal:
     T2(i) = MAX(T2(i), T2_min_fix)
     xion(1:nion,i) = MIN(MAX(xion(1:nion,i), x_MIN),1d0)
     if(r%isH2) then
        ! Ensure the total hydrogen fraction is 1:
        if(xion(ixHI,i)+xion(ixHII,i) .gt. 1d0) then
           if(xion(ixHI,i) .gt. xion(ixHII,i)) then
              xion(ixHI,i)=1d0-xion(ixHII,i)
           else
              xion(ixHII,i)=1d0-xion(ixHI,i)
           endif
        endif ! total hydrogen fraction
     endif ! isH2
     if(r%isHe) then                                        ! Helium species
        ! Ensure the total helium fraction is 1:
        if(xion(ixHeII,i)+xion(ixHeIII,i) .gt. 1d0) then
           if(xion(ixHeII,i) .gt. xion(ixHeIII,i)) then
              xion(ixHeII,i)=1d0-xion(ixHeIII,i)
           else
              xion(ixHeIII,i)=1d0-xion(ixHeII,i)
           endif
        endif
     endif ! isHe
#ifdef DO_RT
     do ig=1,nrtgrp
        Np(ig,i) = MAX(smallNp, Np(ig,i))
        call reduce_flux(Fp(:,ig,i),Np(ig,i)*tables%rt_c_cgs(ilevel))
     end do
#endif
  end do

  ! Loop until all cells have tleft=0
  ! **********************************************
  nAct=nCell                                      ! Currently active cells
  loopcnt=0 !; n_cool_cells=n_cool_cells+nCell     !             Statistics
  do while (nAct .gt. 0)      ! Iterate while there are still active cells
     loopcnt=loopcnt+1 !  ;   tot_cool_loopcnt=tot_cool_loopcnt+nAct
     nAct_next=0                     ! Active cells for the next iteration
     do ia=1,nAct                             ! Loop over the active cells
        i = indAct(ia)                        !                 Cell index
        call cool_step(i)
        if(loopcnt .gt. 100000) then
           call display_coolinfo(.true., loopcnt, i, dt-tleft(i), dt, ddt(i), nH(i), &
#ifdef DO_RT
                &                Np(:,i), Fp(:,:,i), p_gas(:,i), dNp, dFp, dp_gas, ilevel, &
#endif
                &                T2(i), xion(:,i), dT2, dXion, code)
        endif
        if(.not. dt_ok) then
           ddt(i)=ddt(i)/2.                    ! Try again with smaller dt
           nAct_next=nAct_next+1 ; indAct(nAct_next) = i
           loopCodes(code) = loopCodes(code)+1
           cycle
        endif
        ! Update the cell state (advance the time by ddt):
        T2(i) = T2(i) + dT2
        xion(:,i) = xion(:,i) + dXion(:)
#ifdef DO_RT
        Np(:,i) = Np(:,i) + dNp(:)
        Fp(:,:,i) = Fp(:,:,i) + dFp(:,:)
        p_gas(:,i) = p_gas(:,i) + dp_gas(:)
#endif
        tleft(i)=tleft(i)-ddt(i)
        if(tleft(i) .gt. 0.) then           ! Not finished with this cell
           nAct_next=nAct_next+1 ; indAct(nAct_next) = i
        else if(tleft(i) .lt. 0.) then        ! Overshot by abs(tleft(i))
           print*,'In neq_solve_cooling: tleft < 0  !!'
           stop
        endif
        ddt(i)=min(dt_rec,tleft(i))    ! Use recommended dt from cool_step
     end do ! end loop over active cells
     nAct=nAct_next
  end do ! end iterative loop
  ! loop statistics
!  max_cool_loopcnt=max(max_cool_loopcnt,loopcnt)

  end associate

contains

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  SUBROUTINE cool_step(icell)
    ! Compute change in cell state in timestep ddt(icell), or set in dt_rec
    ! a recommendation for new timestep if ddt(icell) proves too large.
    ! T2      => T/mu [K]                               -- dT2 is new value
    ! xion    => nion ionization fractions              --     dXion is new
    ! Np      => nrtgrp photon number densities [cm-3]  -- dNp is new value
    ! Fp      => nrtgrp * ndim photon fluxes [cm-2 s-1] -- dFp is new value
    ! p_gas   => ndim gas momenta [cm s-1 g cm-3]       --    dp_gas is new
    ! dNpdt   =>  Op split increment in photon densities during dt
    ! dFpdt   =>  Op split increment in photon flux magnitudes during dt
    ! nH      =>  Hydrogen number densities [cm-3]
    ! c_switch=>  Cooling switch (1 for cool/heat, 0 for no cool/heat) (OFF)
    ! Zsolar  =>  Cell metallicities [solar fraction]
    ! dt      =>  Timestep size [s]
    ! dt_ok   <=  .f. if timestep constraints were broken, .t. otherwise
    ! dt_rec  <=  Recommended timesteps for next iteration
    ! code    <= Error code in cool step, if dt_ok=.f.
    !
    ! The original values, T2, xion etc, must stay unchanged, while dT2,
    ! dxion etc contain the new values (the difference at the end of the
    ! routine).
    !-----------------------------------------------------------------------
    use amr_commons
    use const
    implicit none
    integer, intent(in):: icell
    !-----------------------------------------------------------------------
    real(kind=8),dimension(nion):: alpha, beta, nN, nI
    real(kind=8):: dUU, fracMax, x_tot
    real(kind=8):: mu, TK, nHe, ne, neInit, Hrate
    real(kind=8):: xHI,dxHI, xH2=0d0,dXH2=0d0, xHeI,dxHeI
    real(kind=8):: Crate, dCdT2, X_nHkb, rate, dRate, cr, de=0d0
    real(kind=8):: photoRate, metal_tot, metal_prime, ss_factor, f_dust
#ifdef DO_RT
    integer:: iion,igroup,idim
    real(kind=8),dimension(ndim):: dmom
    real(kind=8),dimension(nrtgrp):: recRad, phAbs, phSc, dustAbs
    real(kind=8),dimension(nrtgrp):: dustSc, kAbs_loc, kSc_loc
    real(kind=8),dimension(nrtgrp,nion)::signc
    real(kind=8):: rt_c_fraction, rt_c_cgs
    real(kind=8):: TR, one_over_C_v, E_rad, dE_T
    real(kind=8):: G0, eff_peh, cdex
    real(kind=8):: fluxMag, mom_fact
#endif
    real(kind=8):: rho, ncr, Cpeebles, fHm
    logical:: newAtomicCons=.true.
    !-----------------------------------------------------------------------
    associate(ixHI=>r%ixHi,ixHII=>r%ixHII,ixHeII=>r%ixHeII,ixHeIII=>r%ixHeIII)

#ifdef DO_RT
    signc=tables%signc(:,:,ilevel)
    rt_c_fraction = r%rt_c_fraction(ilevel)
    rt_c_cgs = tables%rt_c_cgs(ilevel)
#endif
    dt_ok=.false.
    nHe=0.25*nH(icell)*r%Y_He/r%X_H       ! Helium number density
    ! U contains the original values, dU the updated ones
    dT2 = T2(icell) ; dXion(:) = xion(:,icell)
#ifdef DO_RT
    dNp(:) = Np(:,icell) ; dFp(:,:) = Fp(:,:,icell)
    dp_gas(:) = p_gas(:,icell)
#endif
    ! nN='neutral' species (pre-ionized), nI=their ionized counterparts
    ! nN(1) == nN(ixHI)    == nH2         !      nI(1) == nI(ixHI)    == nHI
    ! nN(2) == nN(ixHII)   == nHI         !     nI(2) == nI(ixHII)   == nHII
    ! nN(3) == nN(ixHeII)  == nHeI        !    nI(3) == nI(ixHeII)  == nHeII
    ! nN(4) == nN(ixHeIII) == nHeII       !   nI(4) == nI(ixHeIII) == nHeIII
    ! Hydrogen chemistry
    xHI = MAX(1d0-dxion(ixHII),x_min)       ! need in case of .not. isH2
    if(r%isH2) xHI = MAX(dxion(ixHI),x_min)
    if(r%isH2) xH2 = MAX((1.-dxion(ixHI)-dxion(ixHII))/2.,x_min)
    ! Helium chemistry
    if(r%isHe) xHeI = MAX(1.-dxion(ixHeII)-dxion(ixHeIII),x_min)
    ! nN='neutral' species (pre-ionized)
    nN=0d0
    if(r%isH2) nN(ixHI) = nH(icell) * xH2                        !      nH2
    nN(ixHII) = nH(icell) * xHI                                  !      nHI
    if(r%isHe) nN(ixHeII)  = nHe*xHeI                            !     nHeI
    if(r%isHe) nN(ixHeIII) = nHe*dxion(ixHeII)                   !    nHeII
    ! nI=ionized counterparts of the neutral species
    nI=0d0
    if(r%isH2) nI(ixHI)  = nN(ixHII)                             !      nHI
    nI(ixHII) = nH(icell) * dxion(ixHII)                         !     nHII
    if(r%isHe) nI(ixHeII)  = nN(ixHeIII)                         !    nHeII
    if(r%isHe) nI(ixHeIII) = nHe*dxion(ixHeIII)                  !   nHeIII
    f_dust = (1.-dxion(ixHII))                     ! No dust in ionised gas

    mu = getMu(r, dxion, dT2)
    TK = dT2 * mu                                        !      Temperature
    if(r%neq_isTconst) TK=r%neq_Tconst                   ! Force constant T
    ne = nH(icell)*dxion(ixHII)
    if(r%isHe) ne=ne+nHe*(dxion(ixHeII)+2.*dxion(ixHeIII)) !  Elec. density
    neInit = ne
    fracMax = 0d0 ! Max fractional update, to check if dt can be increased
    ss_factor = 1d0                  ! UV background self_shielding factor
    if(r%self_shielding) ss_factor = exp(-nH(icell)/1d-2)
    rho = nH(icell) / r%X_H * mH
#ifdef DO_RT
    ! Set dust opacities--------------------------------------------------
    kAbs_loc = r%kappaAbs
    kSc_loc = r%kappaSc
    if(r%is_kIR_T) then                            ! k_IR depends on T_rad
          ! For the radiation temperature,  weigh the energy in each group
            ! by its opacity over IR opacity (derived from IR temperature)
       E_rad = group_egy_erg(iIR) * dNp(iIR)
       TR = max(0d0,(E_rad*rt_c_fraction/a_r)**0.25)      ! IR temperature
       kAbs_loc(iIR) = r%kappaAbs(iIR) * (TR/10d0)**2
       do iGroup=1,nrtgrp
          if(iGroup .ne. iIR)                                            &
               E_rad = E_rad + kAbs_loc(iGroup) / kAbs_loc(iIR)          &
               * group_egy_erg(iGroup) * dNp(iGroup)
       end do
       TR = max(0d0,(E_rad*rt_c_fraction/a_r)**0.25)  ! Rad. temperature
       if(r%rt_T_rad) then      ! Use radiation temperature for everything
          dT2 = TR/mu ;   TK = TR
       endif
                 ! Set the IR opacities according to the rad. temperature:
       kAbs_loc(iIR) = r%kappaAbs(iIR) * (TR/10d0)**2 * exp(-TR/1d3)
       kSc_loc(iIR)  = r%kappaSc(iIR)  * (TR/10d0)**2 * exp(-TR/1d3)
    endif ! if(is_kIR_T)
                         ! Set dust absorption and scattering rates [s-1]:
    dustAbs(:) =kAbs_loc(:) *rho*Zsolar(icell)*f_dust*rt_c_cgs
    dustSc(iIR)=kSc_loc(iIR)*rho*Zsolar(icell)*f_dust*rt_c_cgs

    ! UPDATE PHOTON DENSITY AND FLUX *************************************
    if(r%rt_advect) then
       recRad(1:nrtgrp)=0. ; phAbs(1:nrtgrp)=0.
        ! Scattering rate; reduce the photon flux, but not photon density:
       phSc(1:nrtgrp)=0.

       ! EMISSION FROM GAS
       if(.not. r%rt_otsa .and. r%rt_advect) then  ! ------ Rec. radiation
          if(r%isH2) alpha(ixHI) = 0d0        ! H2 emits no rec. radiation
          alpha(ixHII) = inp_coolrates_table(tables,tables%tbl_alphaA_HII, TK,.false.) &
                       - inp_coolrates_table(tables,tables%tbl_alphaB_HII, TK,.false.)
          if(r%isHe) then
                  ! alpha(2) A-B becomes negative around 1K, hence the max
             alpha(ixHeII) = &
                  MAX(0d0,inp_coolrates_table(tables,tables%tbl_alphaA_HeII,TK,.false.) &
                          -inp_coolrates_table(tables,tables%tbl_alphaB_HeII, TK,.false.))
             alpha(ixHeIII) = inp_coolrates_table(tables,tables%tbl_alphaA_HeIII, TK,.false.) &
                            - inp_coolrates_table(tables,tables%tbl_alphaB_HeIII, TK,.false.)
          endif
          do iion=1,nion
             if(r%spec2group(iion) .gt. 0) &   ! Contribution ion -> group
                  recRad(r%spec2group(iion)) = &
                  recRad(r%spec2group(iion)) + alpha(iion) * nI(iion) * ne
          enddo
       endif

       ! ABSORPTION/SCATTERING OF PHOTONS BY GAS
       do igroup=1,nrtgrp       ! ----------------Ionization absorbtion
          phAbs(igroup) = SUM(nN(:)*signc(igroup,:)*r%ssh2(igroup)) ! s-1
       end do
       ! IR, optical and UV depletion by dust absorption: ----------------
       ! IR scattering/abs on dust (abs after T update)
       if(r%rt_isIR) phSc(iIR)  = phSc(iIR) + dustSc(iIR)
       do igroup=1,nrtgrp      ! Deplete photons, since they go into IR
          if( .not. (r%rt_isIR .and. igroup.eq.iIR) ) & ! IR done elsewhere
               phAbs(igroup) = phAbs(igroup) + dustAbs(igroup)
       end do

       if(r%iPEH_group .gt. 0) then
          ! Photoelectric absorption: the effective PEH cross section
          ! is photoelectric heating rate / Habing flux
          ! Note: as this absorption is done separately, kappaAbs
          !       should not include PEH absorption when PEH is included.
          ! from Bakes and Tielens 1994 and Wolfire 2003
          G0 = group_egy_erg(r%iPEH_group)                               &
             * dNp(r%iPEH_group) * rt_c_cgs / 1.6d-3
          eff_peh = 4.87d-2                                              &
                  / (1d0 + 4d-3 * (G0*sqrt(TK)/ne*2.)**0.73)             &
                  + 3.65d-2 * (TK/1d4)**0.7                              &
                  / (1d0 + 2d-4 * (G0 * sqrt(TK) / ne*2. ))
          phAbs(r%iPEH_group) = phAbs(r%iPEH_group)                      &
                            + 8.125d-22 * eff_peh * rt_c_cgs *nH(icell)  &
                            * Zsolar(icell) * f_dust
       endif

       dmom(1:ndim)=0d0
       do igroup=1,nrtgrp     ! ----------------- Do the update of N and F
          dNp(igroup)= MAX(smallNp,                                      &
                        (ddt(icell)*(recRad(igroup)+dNpdt(igroup,icell)) &
                                    +dNp(igroup))                        &
                        / (1d0+ddt(icell)*phAbs(igroup)))

          dUU = ABS(dNp(igroup)-Np(igroup,icell))                        &
                /(Np(igroup,icell)+Np_MIN) * one_over_Np_FRAC
          if(dUU .gt. 1d0) then
             code=1 ;   RETURN                        ! ddt(icell) too big
          endif
          fracMax=MAX(fracMax,dUU)      ! To check if ddt can be increased

          do idim=1,ndim
             dFp(idim,igroup) = &
                  (ddt(icell)*dFpdt(idim,igroup,icell)+dFp(idim,igroup)) &
                  /(1d0+ddt(icell)*(phAbs(igroup)+phSc(igroup)))
          end do
          call reduce_flux(dFp(:,igroup),dNp(igroup)*rt_c_cgs)

          do idim=1,ndim
             dUU = ABS(dFp(idim,igroup)-Fp(idim,igroup,icell))           &
                  / (ABS(Fp(idim,igroup,icell))+Np_MIN*rt_c_cgs)         &
                  * one_over_Fp_FRAC
             if(dUU .gt. 1d0) then
                code=2 ;   RETURN                     ! ddt(icell) too big
             endif
             fracMax=MAX(fracMax,dUU)   ! To check if ddt can be increased
          end do

       end do

       do igroup=1,nrtgrp ! -------Momentum transfer from photons to gas:
          mom_fact = ddt(icell) * (phAbs(igroup) + phSc(igroup)) &
               * group_egy_erg(igroup) * one_over_clight

          if(r%rt_isoPress .and. .not. (r%rt_isIR .and. igroup==iIR)) then
             ! rt_isoPress: assume f=1, where f is reduced flux.
             fluxMag=sqrt(sum((dFp(:,igroup))**2))
             if(fluxMag .gt. 0d0) then
                mom_fact = mom_fact * dNp(igroup) / fluxMag
             else
                mom_fact = 0d0
             endif
          else
             mom_fact = mom_fact * one_over_rt_c_cgs
          end if

          do idim = 1, ndim
             dmom(idim) = dmom(idim) + dFp(idim,igroup) * mom_fact
          end do
       end do
       dp_gas = dp_gas + dmom * r%rt_pressBoost      ! update gas momentum

       ! Add absorbed UV/optical energy to IR:----------------------------
       if(r%rt_isIR) then
          do igroup=iIR+1,nrtgrp
             dNp(iIR) = dNp(iIR) + dustAbs(igroup) * ddt(icell)          &
                  * dNp(igroup) * group_egy_ratio(igroup)
          end do
       endif
       ! -----------------------------------------------------------------
    endif !if(rt)
#endif

    ! UPDATE TEMPERATURE *************************************************
    !if(c_switch(icell) .and. .not. rt_isTconst .and. .not. r%rt_T_rad) then
    if(.not. r%neq_isTconst .and. .not. r%rt_T_rad) then
       Hrate = 0.                           !  Heating rate [erg cm-3 s-1]
       if(r%haardt_madau) Hrate = Hrate + SUM(nN(:)*tables%UVrates(:,2)) * ss_factor
#ifdef DO_RT
       if(r%rt_advect) then
          do igroup=1,nrtgrp                            !  Photoheating
             Hrate = Hrate + dNp(igroup) * SUM(nN(:)                     &
                   * tables%PHrate(igroup,:,ilevel))
          end do
       endif
       if(r%iPEH_group .gt. 0 .and. r%rt_advect) then
                              ! Photoelectric heating Bakes & Tielens 1994
                                        ! and Wolfire 2003, [erg cm-3 s-1]
          Hrate = Hrate + 1.3d-24 * eff_peh * G0 * nH(icell)             &
                * Zsolar(icell) * f_dust
       endif
#endif
       if(r%isH2) then
#ifdef DO_RT
                                              ! UV pumping, Baczynski 2015
          cdex  = 1d-12 * (1.4 * exp(-18100. / (TK + 1200.)) * xH2       &
                + exp(-1000. / TK) * dxion(ixHI))                        &
                * sqrt(TK) * nH(icell)                             ! [s-1]
          Hrate = Hrate + 6.94 * SUM(dNp(:) * signc(:,ixHI)       &
                * r%isLW(:)) * 2. * eV2erg * cdex / (cdex + 2d-7)        &
                * nH(icell) * xH2
#endif
                                       ! H2 formation heating, Omukai 2000
                           ! and Hollenbach and McKee 1976, [erg cm-3 s-1]

          ! Fraction of H- surviving CMB photodissociation
          fHm = 1.0d0 / (1.0d0                                           &
                + inp_coolrates_table(tables,tables%tbl_beta_cmb_Hm,TK,.false.) &
                / nN(ixHII))

          ncr   = 1d6 * TK**(-0.5)                                       &
                / (1.6 * dxion(ixHI) * exp(-(400. / TK)**2)              &
                + 1.4 * xH2 * exp(-12000. / (TK + 1200.)))        ! [cm-3]
          Hrate = Hrate + eV2erg                                         &
                * ((0.2 + 4.2 / (1. + ncr / nH(icell)))                  &
                * inp_coolrates_table(tables,tables%tbl_AlphaZ_H2,TK,.false.) &
                * Zsolar(icell) * f_dust                                 &
                * nH(icell)**2 * dxion(ixHI)                             &
                + 3.53 / (1. + ncr/nH(icell))                            &
                * inp_coolrates_table(tables,tables%tbl_AlphaGP_H2,TK,.false.) &
                * fHm * nH(icell) * dxion(ixHI) * ne                     &
                + 4.48 / (1.+ncr / nH(icell))                            &
                * inp_coolrates_table(tables,tables%tbl_Beta_H3B,TK,.false.) &
                * nH(icell)**3 * dxion(ixHI)**2 * (dxion(ixHI) + xH2/8.))
       endif
       if (r%cosmic_rays) then                 ! CR heating [erg cm-3 s-1]
                                       ! Glassgold 2012, ~10 ev/ionisation
          Hrate = Hrate + 10. * eV2erg                                   &
                * nH(icell) * xHI * cosray_HI
          if (r%isH2) Hrate = Hrate + 10. * eV2erg                       &
                          * nH(icell) * xH2 * cosray_H2
          if (r%isHe) Hrate = Hrate + 10. * eV2erg                       &
                          * nHe * xHeI * 1.1 * cosray_HI
       endif
       Crate = compCoolrate(r,tables,TK,ne,nN,nI,dCdT2)          ! Cooling
       dCdT2 = dCdT2 * mu                            ! dC/dT2 = mu * dC/dT
       metal_tot=0d0 ; metal_prime=0d0                     ! Metal cooling
       if(Zsolar(icell) .gt. 0d0) &
            call neq_cmp_metals(tables, T2(icell), nH(icell), mu,        &
            &                  metal_tot, metal_prime)
       X_nHkb = r%X_H/(1.5 * nH(icell) * kB)        ! Multiplication factor
       rate  = X_nHkb*(Hrate - Crate - Zsolar(icell)*metal_tot)
       dRate = -X_nHkb*(dCdT2 + Zsolar(icell)*metal_prime)     ! dRate/dT2
                                                     ! 1st order dt constr
       dUU   = ABS(MAX(T2_min_fix, T2(icell)+rate*ddt(icell))-T2(icell))
                                                            ! New T2 value
       dT2   = MAX(T2_min_fix &
                  ,T2(icell)+rate*ddt(icell)/(1.-dRate*ddt(icell)))
       dUU   = MAX(dUU, ABS(dT2-T2(icell))) / (T2(icell)+T_MIN) &
                        *one_over_T_FRAC
       if(dUU .gt. 1.) then                                     ! 10% rule
          code=3 ; RETURN
       endif
       fracMax=MAX(fracMax,dUU)
       TK=dT2*mu
    endif

#ifdef DO_RT
    if(r%rt_isIR) then
       if(kAbs_loc(iIR) .gt. 0d0 .and. .not. r%rt_T_rad) then
          ! Evolve IR-Dust equilibrium temperature------------------------
          ! Delta (Cv T)= ( c_red/lambda E - c/lambda a T^4)
          !           / ( 1/Delta t + 4 c/lambda/C_v a T^3 + c_red/lambda)
          one_over_C_v = mH*mu*(r%gamma-1d0) / (rho*kB)
          E_rad = group_egy_erg(iIR) * dNp(iIR)
          dE_T = (rt_c_cgs * E_rad - c_cgs*a_r*TK**4)                    &
               /(1d0/(kAbs_loc(iIR) * Zsolar(icell) * rho * ddt(icell))  &
               +4d0*c_cgs * one_over_C_v *a_r*TK**3+rt_c_cgs)
          dT2 = dT2 + 1d0/mu * one_over_C_v * dE_T
          dNp(iIR) = dNp(iIR) - dE_T * one_over_egy_IR_erg

          dT2 = max(T2_min_fix,dT2)
          dNp(iIR) = max(dNp(iIR), smallNp)
          ! 10% rule for photon density:
          dUU = ABS(dNp(iIR)-Np(iIR,icell)) / (Np(iIR,icell)+Np_MIN)     &
                                            * one_over_Np_FRAC
          if(dUU .gt. 1.) then
             code=4 ;   RETURN
          endif
          fracMax=MAX(fracMax,dUU)

          dUU   = ABS(dT2-T2(icell)) / (T2(icell)+T_MIN) * one_over_T_FRAC
          if(dUU .gt. 1.) then                           ! 10% rule for T2
             code=5 ; RETURN
          endif
          fracMax=MAX(fracMax,dUU)
          TK=dT2*mu
          call reduce_flux(dFp(:,iIR),dNp(iIR)*rt_c_cgs)
       endif
    endif
#endif

    ! HYDROGEN UPDATE*****************************************************
    ! Update xH2**********************************************************
    dxH2=xH2
    if(r%isH2) then
       alpha(ixHI) = inp_coolrates_table(tables,tables%tbl_AlphaZ_H2, TK,.false.)      &
                   * Zsolar(icell) * f_dust * nH(icell)                  &
                   + inp_coolrates_table(tables,tables%tbl_AlphaGP_H2,TK,.false.)      &
                   * fHm * ne                                            &
                   + inp_coolrates_table(tables,tables%tbl_Beta_H3B,TK,.false.)        &
                   * nH(icell)**2 * dxion(ixHI) * (dxion(ixHI) + xH2/ 8.)
       beta(ixHI)  = inp_coolrates_table(tables,tables%tbl_Beta_H2HI, TK,.false.)      &
                   * dxion(ixHI)                                         &
                   + inp_coolrates_table(tables,tables%tbl_Beta_H2H2, TK,.false.) * xH2
       cr = alpha(ixHI) * dxion(ixHI)                        ! H2 Creation
       photoRate=0.
#ifdef DO_RT
       photoRate = SUM(signc(:,ixHI)*dNp)
#endif
       if(r%haardt_madau) photoRate = photoRate + tables%UVrates(ixHI,1)*ss_factor
       de = beta(ixHI) * nH(icell) + photoRate            ! H2 Destruction
       if(r%cosmic_rays) de = de + cosray_H2
       dxH2 = (cr*ddt(icell)+xH2)/(1.+de*ddt(icell))
       dxH2 = MIN(MAX(dxH2, x_min), 0.5d0)
    endif !if(isH2)
    ! Update xHI (also if .not. isH2, for stability)**********************
    if(r%rt_otsa .or. .not. r%rt_advect) then     !    Recombination rates
       alpha(ixHII) = inp_coolrates_table(tables,tables%tbl_AlphaB_HII, TK,.false.)
    else
       alpha(ixHII) = inp_coolrates_table(tables,tables%tbl_AlphaA_HII, TK,.false.)
    endif
    ! Peebles C factor
    Cpeebles = 1d0 / (1d0                                                &
        + inp_coolrates_table(tables,tables%tbl_beta_cmb_HI, TK,.false.) &
        / (8.2206d0 + 1d0 / (nN(ixHII)                                   &
        * inp_coolrates_table(tables,tables%tbl_beta_ly_red, TK,.false.))))

    alpha(ixHII) = alpha(ixHII) * Cpeebles
    beta(ixHII) = inp_coolrates_table(tables,tables%tbl_Beta_HI, TK,.false.) ! Coll.ion.rate
    cr = alpha(ixHII) * ne * dxion(ixHII) + 2. * de * dxH2    !HI creation
    photoRate=0.
#ifdef DO_RT
    photoRate = SUM(signc(:,ixHII)*dNp)           !                  [s-1]
#endif
    photoRate = photoRate + Cpeebles                                    &
        * inp_coolrates_table(tables,tables%tbl_gamma_cmb_HI, TK,.false.)
    if(r%haardt_madau) photoRate = photoRate + tables%UVrates(ixHII,1)*ss_factor
    de = beta(ixHII) * ne + photoRate             !         HI destruction
    if(r%cosmic_rays) de = de + cosray_HI
    if(r%isH2) de = de + 2. * alpha(ixHI)
    dxHI = (cr*ddt(icell)+xHI)/(1.+de*ddt(icell))
    dxHI = MIN(MAX(dxHI, x_min),1d0)
    if(r%isH2) dxion(ixHI)=dxHI
    ! Update xHII*********************************************************
    cr = (beta(ixHII)*ne+photoRate)*dxHI            !             Creation
    if(r%cosmic_rays) cr = cr + cosray_HI * dxHI
    de = alpha(ixHII) * ne                          !          Destruction
    dxion(ixHII) = (cr*ddt(icell)+dxion(ixHII))/(1.+de*ddt(icell))
    dxion(ixHII) = MIN(MAX(dxion(ixHII), x_MIN),1d0)
    ! Atomic conservation of H *******************************************
    if(newAtomicCons) then
       x_tot = 2.*dxH2 + dxHI + dxion(ixHII)
       dxH2 = dxH2/x_tot
       dxHI = dxHI/x_tot
       if(r%isH2) dxion(ixHI) = dxHI
       dxion(ixHII)=dxion(ixHII)/x_tot
    else
       if(r%isH2) then
          if(dxH2.ge.dxion(ixHII)) then      !   H2 or HI is most abundant
             if(dxH2.le.dxion(ixHI)) &       !    -> HI
                  dxion(ixHI)=1.-2.*dxH2-dxion(ixHII)
          else                               !  HI or HII is most abundant
             if(dxion(ixHI).le.dxion(ixHII)) then
                dxion(ixHII) = 1.-2.*dxH2-dxion(ixHI)             ! -> HII
             else
                dxion(ixHI) = 1.-2.*dxH2-dxion(ixHII)             ! ->  HI
             endif
          endif
       else
          if(dxHI.le.dxion(ixHII)) dxion(ixHII) = 1.-dxHI         ! -> HII
       endif ! if(isH2)
    endif ! if(newAtomicCons)
    ! Timestep ok for hydrogen update? ***********************************
    dUU = 0d0
    if(r%isH2) dUU = MAX(dUU,ABS((dxH2-xH2)/(xH2+x_FM)))
    dUU = MAX(dUU,ABS((dxHI-xHI)/(xHI+x_FM)))
    dUU = MAX(dUU &
         ,ABS(dXion(ixHII)-xion(ixHII,icell))/(xion(ixHII,icell)+x_FM))  &
        * one_over_x_FRAC
    if(dUU .gt. 1.) then
       code=6 ; RETURN
    endif
    fracMax=MAX(fracMax,dUU)

    ! HELIUM UPDATE*******************************************************
    if(r%isHe) then
       ! Update ne because of changed hydrogen ionisation:
       ne= nH(icell)*dXion(ixHII)+nHE*(dXion(ixHeII)+2.*dXion(ixHeIII))
       mu = getMu(r, dXion, dT2)
       if(.not. r%neq_isTconst) TK=dT2*mu! Update TK because of changed  mu
       ! Update xHeI *****************************************************
       if(r%rt_otsa .or. .not. r%rt_advect) then
          alpha(ixHeII)  = inp_coolrates_table(tables,tables%tbl_alphaB_HeII, TK,.false.)
          alpha(ixHeIII) = inp_coolrates_table(tables,tables%tbl_alphaB_HeIII, TK,.false.)
       else
          alpha(ixHeII)  = inp_coolrates_table(tables,tables%tbl_alphaA_HeII, TK,.false.)
          alpha(ixHeIII) = inp_coolrates_table(tables,tables%tbl_alphaA_HeIII, TK,.false.)
       endif
       beta(ixHeII)  =  inp_coolrates_table(tables,tables%tbl_beta_HeI, TK,.false.)
       beta(ixHeIII) = inp_coolrates_table(tables,tables%tbl_beta_HeII, TK,.false.)
       ! Creation = recombination of HeII and electrons
       cr = alpha(ixHeII) * ne * dXion(ixHeII)
       ! Destruction = collisional ionization+photoionization of HeI
       de = beta(ixHeII) * ne
       if(r%cosmic_rays) de = de + 1.1 * cosray_HI
#ifdef DO_RT
       de = de + SUM(signc(:,ixHeII)*dNp)
#endif
       if(r%haardt_madau) de = de + tables%UVrates(ixHeII,1) * ss_factor
       dxHeI = (cr*ddt(icell)+xHeI)/(1.+de*ddt(icell))         ! The update
       dxHeI = MIN(MAX(dxHeI, x_min),1d0)
       ! Update xHeII ****************************************************
       ! Creation = coll.- and photo-ionization of HI + rec. of HeIII
       cr = de * xHeI + alpha(ixHeIII) * ne * dXion(ixHeIII)
       ! Destruction = rec. of HeII + coll.- and photo-ionization of HeII
       photoRate = 0.
#ifdef DO_RT
       photoRate = SUM(signc(:,ixHeIII)*dNp)
#endif
       if(r%haardt_madau) photoRate = photoRate + tables%UVrates(ixHeIII,1) * ss_factor
       de = (alpha(ixHeII) + beta(ixHeIII)) * ne + photoRate
       dXion(ixHeII) = (cr*ddt(icell)+dXion(ixHeII))/(1.+de*ddt(icell))
       dXion(ixHeII) = MIN(MAX(dXion(ixHeII), x_MIN),1d0)
       ! Update xHeIII ***************************************************
       ! Creation = coll.- and photo-ionization of HeII
       cr = (beta(ixHeIII) * ne + photoRate) * dXion(ixHeII)   ! new xHeII
       ! Destruction = rec. of HeIII and e
       de = alpha(ixHeIII) * ne
       dXion(ixHeIII) = (cr*ddt(icell)+dXion(ixHeIII))/(1.+de*ddt(icell))
       dXion(ixHeIII) = MIN(MAX(dXion(ixHeIII), x_MIN),1d0)
       ! Atomic conservation of He ***************************************
       if(newAtomicCons) then
          x_tot = dxHeI + dXion(ixHeII) + dxion(ixHeIII)
          dxHeI          = dxHeI / x_tot
          dxion(ixHeII)  = dxion(ixHeII)/x_tot
          dxion(ixHeIII) = dxion(ixHeIII)/x_tot
       else
          if(dxHeI .ge. dXion(ixHeIII)) then! Either HeI or HeII most abundant
             if(dxHeI.le.dXion(ixHeII)) dXion(ixHeII) = 1.-dxHeI-dXion(ixHeIII)
          else                        ! Either HeII or HeIII is most abundant
             if(dxion(ixHeII) .le. dxion(ixHeIII)) then
                dxion(ixHeIII) = 1. - dxHeI-dxion(ixHeII)
             else
                dxion(ixHeII) = 1. - dxHeI-dxion(ixHeIII)              !  HeII
             endif
          endif
       endif
       ! Timestep ok for helium update?
       dUU = ABS(dxHeI-xHeI)/(xHeI+x_FM)
       dUU = MAX(dUU, ABS(dxion(ixHeII)-xion(ixHeII,icell)) &
                      /(xion(ixHeII,icell)+x_FM))
       dUU = MAX(dUU, ABS(dxion(ixHeIII)-xion(ixHeIII,icell)) &
                      /(xion(ixHeIII,icell)+x_FM)) &
           * one_over_x_FRAC
       if(dUU .gt. 1.) then
          code=7 ; RETURN
       endif
       fracMax=MAX(fracMax,dUU)
    endif !if(isHe) ! END HELIUM UPDATE***********************************

    ! CHECK FRACTIONAL CHANGE IN ELECTRON ABUNDANCE **********************
    ne = nH(icell)*dxion(ixHII)
    if(r%isHe) ne = ne + nHe*(dxion(ixHeII)+2.*dxion(ixHeIII))
    dUU=ABS((ne-neInit)) / (neInit+x_FM) * one_over_x_FRAC
    if(dUU .gt. 1.) then
       code=8 ; RETURN
    endif
    fracMax=MAX(fracMax,dUU)

    if(r%neq_isTconst)then
       mu = getMu(r, dXion, dT2)
       dT2 = r%neq_Tconst/mu
    endif

    ! CLEAN UP AND RETURN ************************************************
    dT2 = dT2-T2(icell) ; dXion(:) = dXion(:)-xion(:,icell)
#ifdef DO_RT
    dNp(:) = dNp(:)-Np(:,icell) ; dFp(:,:) = dFp(:,:)-Fp(:,:,icell)
    dp_gas(:)= dp_gas(:)-p_gas(:,icell)
#endif
    ! Now the dUs are really changes, not new values
    ! Check if we are safe to use a bigger timestep in next iteration:
    if(fracMax .lt. 0.5) then
       dt_rec=ddt(icell)*2.
    else
       dt_rec=ddt(icell)
    endif
    dt_ok=.true.
    code=0

    end associate

  END SUBROUTINE cool_step

  !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  SUBROUTINE display_coolinfo(stopRun, loopcnt, i, dtDone, dt, ddt, nH,    &
#ifdef DO_RT
       &                      Np,  Fp,  p_gas, dNp, dFp, dp_gas,ilevel,    &
#endif
       &                      T2,  xion, dT2, dXion, code)
    ! Print cooling information to standard output, and maybe stop execution.
    !------------------------------------------------------------------------
    real(kind=8),dimension(nion):: xion, dXion
#ifdef DO_RT
    real(kind=8),dimension(nrtgrp):: Np, dNp
    real(kind=8),dimension(ndim, nrtgrp):: Fp, dFp
    real(kind=8),dimension(ndim):: p_gas, dp_gas
    integer::ilevel
#endif
    real(kind=8)::T2, dT2, dtDone, dt, ddt, nH
    logical::stopRun
    integer::loopcnt,i, code
    !------------------------------------------------------------------------
    if(stopRun) write(*, 111) loopcnt
    if(.true.) then
#ifdef DO_RT
       write(*,900) loopcnt, code, i, dtDone, dt, ddt, tables%rt_c_cgs(ilevel), nH
       write(*,901) T2,      xion,      Np,      Fp,      p_gas
       write(*,902) dT2,     dXion,     dNp,     dFp,     dp_gas
       write(*,903) dT2/ddt, dXion/ddt, dNp/ddt, dFp/ddt, dp_gas/ddt
       write(*,904) abs(dT2)/(T2+T_MIN), abs(dxion)/(xion+x_FM),          &
            abs(dNp)/(Np+Np_MIN), abs(dFp)/(Fp+Np_MIN*tables%rt_c_cgs(ilevel))
#else
       write(*,900) loopcnt, code, i, dtDone, dt, ddt, nH
       write(*,901) T2,      xion
       write(*,902) dT2,     dXion
       write(*,903) dT2/ddt, dXion/ddt
       write(*,904) abs(dT2)/(T2+T_MIN), abs(dxion)/(xion+x_FM)
#endif
    endif
    print*,loopcodes
    if(stopRun) then
       print *,'XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
       STOP
    endif

111 format(' Stopping because of large number of timestesps in', &
         ' neq_solve_cooling (', I6, ')')
900 format (I3, ' code=', I2, ' i=', I5, ' t=', 1pe12.3,xs&
         '/', 1pe12.3, ' ddt=', 1pe12.3, ' c=', 1pe12.3, &
         ' nH=', 1pe12.3)
901 format ('  U      =', 20(1pe12.3))
902 format ('  dU     =', 20(1pe12.3))
903 format ('  dU/dt  =', 20(1pe12.3))
904 format ('  dU/U % =', 20(1pe12.3))
  END SUBROUTINE display_coolinfo

END SUBROUTINE neq_solve_cooling

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE cmp_equilibrium_abundances(r, tables, &
                                   & T2, nH, phI_rates, mu, nSpec, Zsolar)

  ! Compute chemical equilibrium abundances of e,H2,HI,HII,HeI,HeII,HeIII.
  ! tables     => Object containing all non-equilibrium cooling tables
  ! T2         => T/mu in Kelvin
  ! nH         => Hydrogen density in cm^-3
  ! phI_rates  => Photoionization rates [s-1] for H2, HI, HeI, HeII
  ! Zsolar     => Metallicity in Solar units
  ! nSpec      <= Resulting species number densities
  ! mu         <= Resulting average particle mass in units of proton mass
  !-------------------------------------------------------------------------
  implicit none
  type(run_t)::r
  type(neq_cooling_t)::tables
  real(kind=8)::T2,nH
  real(kind=8),dimension(nion)::phI_rates
  real(kind=8)::mu,Zsolar
  real(kind=8),dimension(1:7)::nSpec
  !-------------------------------------------------------------------------
  real(kind=8)::mu_mid, err_mu, mu_left, mu_right, T
  integer::niter
  !-------------------------------------------------------------------------
  ! Iteration to find mu                     ! n_E     = n_spec(1) ! e
                                             ! n_H2    = n_spec(2) ! H2
  err_mu=1.                                  ! n_HI    = n_spec(3) ! H
  mu_left=0.5                                ! n_HII   = n_spec(4) ! H+
  mu_right=2.3                               ! n_HEI   = n_spec(5) ! He
  niter=0                                    ! n_HEII  = n_spec(6) ! He+
                                             ! n_HEIII = n_spec(7) ! He++
  do while (err_mu > 1d-4 .and. niter <= 50)
     mu_mid=0.5*(mu_left+mu_right)
     T = T2*mu_mid
     call cmp_chem_eq(r, tables, T, nH, phI_rates, nSpec, mu, Zsolar)
     if(mu - mu_mid > 0.)then
        mu_left = mu_mid
     else
        mu_right= mu_mid
     end if
     err_mu = ABS(mu_right-mu_left)/mu_mid
     niter=niter+1
  end do
  if (niter > 50) then
     write(*,*) 'ERROR in cmp_equilibrium_abundances : too many iterations.'
     STOP
  endif

END SUBROUTINE cmp_equilibrium_abundances

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE cmp_chem_eq(r, tables, TK, nH, t_rad_spec, nSpec, mu, Zsol)

  ! Compute chemical equilibrium abundances of e,H2,HI,HII,HeI,HeII,HeIII.
  ! tables     => Object containing all non-equilibrium cooling tables
  ! TK         => Temperature in Kelvin
  ! nH         => Hydrogen density in cm^-3
  ! t_rad_spec => Photoionization rates [s-1] for H2, HI, HeI, HeII
  ! Zsol       => Metallicity in Solar units
  ! nSpec      <= Resulting species number densities
  ! mu         <= Resulting average particle mass in units of proton mass
  !------------------------------------------------------------------------
  implicit none
  type(run_t)::r
  type(neq_cooling_t)::tables
  real(kind=8),intent(in)::TK, nH, Zsol
  real(kind=8),intent(out)::mu
  real(kind=8),dimension(nion),intent(in)::t_rad_spec
  real(kind=8),dimension(1:7),intent(out)::nSpec
  !------------------------------------------------------------------------
  real(kind=8)::nHe, nTot
  real(kind=8)::n_H2, n_HI, n_HII, n_HEI, n_HEII, n_HEIII, n_E, n_E_min
  real(kind=8)::g_H2=0,   g_HI=0,    g_HEI=0, g_HEII=0,  g_cmb_HI=0   ! Photoion/dissoc
  real(kind=8)::aZ_H2=0,  aGP_H2=0,  a_HI=0,  a_HEI=0,   a_HEII=0  ! Formation
  real(kind=8)::b_H2HI=0, b_H2H2=0,  b_H3B,   b_HI=0,    b_HEI=0, b_HEII=0!Col
  real(kind=8)::b_Hm=0,   fHm=0                                  ! CMB
  real(kind=8)::C_HII=0,  C_H2=0,    D_H2=0,  f_HII=0,   f_H2=0  ! Cre & destr
  real(kind=8)::D_HEI=0,  C_HEIII=0, f_HeI=0, f_HeIII=0, f_dust=0! Cre & destr
  real(kind=8)::err_nE, err_nH2, n_H2_old, n_H2_new, beta_H2
  !-------------------------------------------------------------------------
  associate(ixHI=>r%ixHi, ixHII=>r%ixHII, ixHeII=>r%ixHeII, ixHeIII=>r%ixHeIII)

  beta_H2 = 0.1                            ! Stabilising H2 damping factor
  g_HI   = t_rad_spec(ixHII)               !         Photoionization [s-1]
  if(r%isH2) then
     g_H2   = t_rad_spec(ixHI)                !    Photodissociation [s-1]
     aZ_H2  = inp_coolrates_table(tables,tables%tbl_AlphaZ_H2, TK,.false.)  ! Dust form   [cm3 s-1]
     aGP_H2 = inp_coolrates_table(tables,tables%tbl_AlphaGP_H2, TK,.false.) ! Gas phase   [cm3 s-1]
     b_H3B  = inp_coolrates_table(tables,tables%tbl_Beta_H3B, TK,.false.)   ! 3 body H2   [cm3 s-1]
     b_H2HI = inp_coolrates_table(tables,tables%tbl_Beta_H2HI, TK,.false.)  !     Cdiss   [cm3 s-1]
     b_H2H2 = inp_coolrates_table(tables,tables%tbl_Beta_H2H2, TK,.false.)  !     Cdiss   [cm3 s-1]
     b_Hm   = inp_coolrates_table(tables,tables%tbl_Beta_cmb_Hm, TK,.false.)! CMB H- diss [cm3 s-1]
  endif
  if(r%rt_otsa) then                           !   Recombination [cm3 s-1]
     a_HI   = inp_coolrates_table(tables,tables%tbl_AlphaB_HII, TK,.false.)
     a_HEI  = inp_coolrates_table(tables,tables%tbl_AlphaB_HeII, TK,.false.)
     a_HEII = inp_coolrates_table(tables,tables%tbl_AlphaB_HeIII, TK,.false.)
  else
     a_HI   = inp_coolrates_table(tables,tables%tbl_AlphaA_HII, TK,.false.)
     a_HEI  = inp_coolrates_table(tables,tables%tbl_AlphaA_HeII, TK,.false.)
     a_HEII = inp_coolrates_table(tables,tables%tbl_AlphaA_HeIII, TK,.false.)
  endif
  b_HI     = inp_coolrates_table(tables,tables%tbl_Beta_HI, TK,.false.)      !  Cion [cm3 s-1]
  g_cmb_HI = inp_coolrates_table(tables,tables%tbl_gamma_cmb_HI, TK,.false.) !  CMBion [cm3 s-1]
  if(r%isHe) then
     nHe = r%Y_He/(1.-r%Y_He)/4.*nH
     g_HEI  = t_rad_spec(ixHeII)
     g_HEII = t_rad_spec(ixHeIII)
     b_HEI  = inp_coolrates_table(tables,tables%tbl_Beta_HeI, TK,.false.)
     b_HEII = inp_coolrates_table(tables,tables%tbl_Beta_HeII, TK,.false.)
  endif

  n_E = nH     ; n_H2 = 0d0    ; n_H2_old = nH/2d0
  n_HeI = 0d0  ; n_HeII=0d0    ; n_HeIII=0d0
  err_nE = 1d0 ; err_nH2=0d0   ! err_nH2 initialisation in case of no H2
  n_HI = 0.0d0 ; n_HII=nH

  do while(err_nE > 1d-8 .or. err_nH2 > 1d-8)
     n_E_min = MAX(n_E,1e-15*nH)
     C_HII = b_HI * n_E_min + g_HI + g_cmb_HI !  HII creation (s-1)
     if(r%cosmic_rays) C_HII = C_HII + cosray_HI
     f_HII = max(C_HII / a_HI / n_E_min, 1d-20)     ! Cre/Destr [unitless]
     f_H2 = 0d0
     if(r%isH2) then
        f_dust = 1d0-n_HII/nH
        fHm = 1.0d0 / (1.0d0 + b_Hm / n_HI)
        C_H2   = aZ_H2 * nH * Zsol * f_dust                              &
               + fHm * aGP_H2 * n_E_min                                  &
               + b_H3B * n_HI * (n_HI + n_H2/ 8.)
        D_H2   = b_H2HI * n_HI + b_H2H2 * n_H2 + g_H2    ! H2 destr. (s-1)
        if(r%cosmic_rays) D_H2 = D_H2 + cosray_H2
        f_H2   = C_H2 / max(D_H2,1d-50)             ! Cre/Destr [unitless]
        n_H2_new = nH / (2d0 + 1d0/f_H2 + f_HII/f_H2)
        n_H2     = (1d0-beta_H2) * n_H2 + beta_H2 * n_H2_new
     endif ! if(isH2)
     n_HI  = nH / (1d0 + f_HII + 2d0*f_H2)
     n_HII = nH / (1d0 + 1d0/f_HII + 2d0*f_H2/f_HII)
     if(r%isHe) then
        D_HeI   = b_HEI*n_E_min  + g_HEI               !  HeI destr. (s-1)
        if(r%cosmic_rays) D_HeI = D_HeI + 1.1 * cosray_HI
        C_HeIII = max(b_HEII*n_E_min + g_HEII,1d-99)   !  HeIII cre. (s-1)
        f_HeI   = D_HeI / a_HeI / n_E_min          !  Destr/Cre [unitless]
        f_HeIII = a_HeII * n_E_min / C_HeIII       !  Destr/Cre [unitless]

        n_HEI   = nHe / (1d0 + f_HeI + f_HeI/max(f_HeIII,1d-99))
        n_HEII  = nHe / (1d0 + 1d0/max(f_HeI,1d-99) + 1d0/max(f_HeIII,1d-99))
        n_HEIII = nHe / (1d0 + f_HeIII + f_HeIII/max(f_HeI,1d-99))
     endif ! if(isHe)

     err_nE = ABS((n_E - (n_HII + n_HEII + 2.*n_HEIII))/nH)
     n_E = 0.5*n_E+0.5*(n_HII + n_HEII + 2.*n_HEIII)

     if(r%isH2) then
        err_nH2 = ABS((n_H2_old - n_H2)/nH)
        n_H2_old = n_H2
     endif

  end do

  nTOT     = n_E+n_H2+n_HI+n_HII+n_HEI+n_HEII+n_HEIII
  mu       = nH/(1.-r%Y_He)/nTOT
  nSpec(1) = n_E
  nSpec(2) = n_H2
  nSpec(3) = n_HI
  nSpec(4) = n_HII
  nSpec(5) = n_HEI
  nSpec(6) = n_HEII
  nSpec(7) = n_HEIII

  end associate

END SUBROUTINE cmp_chem_eq

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE neq_evol_single_cell(r, tables, astart, aend, dasura, &
     & h, omegab, omega0, omegaL, T2end, mu, ne, xend, if_write_result)
  !-------------------------------------------------------------------------
  ! Used for initialization of thermal state in cosmological simulations.
  !
  ! tables : object containing all non-equilibrium cooling tables
  ! astart : valeur du facteur d'expansion au debut du calcul
  ! aend   : valeur du facteur d'expansion a la fin du calcul
  ! dasura : la valeur de da/a entre 2 pas de temps
  ! h      : la valeur de H0/100
  ! omegab : la valeur de Omega baryons
  ! omega0 : la valeur de Omega matiere (total)
  ! omegaL : la valeur de Omega Lambda
  ! T2end  : Le T/mu en output
  ! mu     : le poids moleculaire en output
  ! ne     : le ne en output
  ! xend   : les fractions d'ionisation en output
  ! if_write_result : .true. pour ecrire l'evolution de la temperature
  !          et de n_e sur l'ecran.
  !-------------------------------------------------------------------------
!  use UV_module
  implicit none
  type(run_t)::r
  type(neq_cooling_t)::tables
  real(kind=8)::astart,aend,T2end,h,omegab,omega0,omegaL,ne,dasura
  real(kind=8),dimension(1:nion)::xend
  logical::if_write_result
  !-------------------------------------------------------------------------
  real(kind=8)::aexp, daexp, dt_cool, T2_com, nH_com
  real(kind=8),dimension(nion)::pHI_rates
  real(kind=8)::mu
  real(kind=8)::n_spec(1:7)
  real(kind=8),dimension(1:nvector):: T2
  real(kind=8),dimension(1:nion, 1:nvector):: xion
#ifdef DO_RT
  real(kind=8),dimension(1:nrtgrp, 1:nvector):: Np, dNpdt
  real(kind=8),dimension(1:ndim, 1:nrtgrp, 1:nvector):: Fp, dFpdt
  real(kind=8),dimension(1:ndim, 1:nvector):: p_gas
#endif
  real(kind=8),dimension(1:nvector)::nH, Zsolar
!  logical,dimension(1:nvector)::c_switch
  !-------------------------------------------------------------------------
  associate(ixHI=>r%ixHi, ixHII=>r%ixHII, ixHeII=>r%ixHeII, ixHeIII=>r%ixHeIII)

  aexp = astart
  T2_com = 2.726d0 / aexp * aexp**2                    ! Temperature T [K]
  nH_com = omegab*rhoc*h**2*r%X_H/mH
  pHI_rates = 0.                              ! Initially no UV background

  call cmp_chem_eq(r, tables, T2_com/aexp**2, nH_com/aexp**3             &
        ,phI_rates, n_Spec, mu, 0d0)

  T2_com = T2_com / mu                              ! Temperature T/mu [K]

  ! Initialize cell state
  T2(1)=T2_com                                          !      Temperature
  if(r%isH2) xion(ixHI,1)=n_Spec(3)/(nH_com/aexp**3)    !   HI    fraction
  xion(ixHII,1)=n_Spec(4)/(nH_com/aexp**3)              !   HII   fraction
  if(r%isHe) xion(ixHeII,1)=n_Spec(6)/(nH_com/aexp**3)  !   HeII  fraction
  if(r%isHe) xion(ixHeIII,1)=n_Spec(7)/(nH_com/aexp**3) !   HeIII fraction
#ifdef DO_RT
  Np(:,1)=0. ; Fp(:,:,1)=0.                  ! Photon densities and fluxes
  dNpdt(:,1)=0. ; dFpdt(:,:,1)=0.              ! are set to zero initially
  p_gas(:,1)=0.                                 ! Gas momemtum set to zero
#endif
  Zsolar(1)=0.                                   ! Metallicity set to zero
!  c_switch(1)=.true.
  do while (aexp < aend)
!     call update_UVrates(aexp)
     call update_coolrates_tables(r, tables, aexp, 1d7 * h / Mpc2cm      &
         * HsurH0(1d0/dble(aexp)-1d0, omega0, omegaL, 1.0d0-omega0-omegaL))! Update Compton heating

     daexp = dasura*aexp
     dt_cool = daexp                                                     &
             / (aexp * 1d7 * h / Mpc2cm)                               &
             / HsurH0(1.0/dble(aexp)-1.,omega0,omegaL,1.-omega0-omegaL)

     nH(1) = nH_com/aexp**3
     T2(1) = T2(1)/aexp**2
     call neq_solve_cooling(r, tables, T2, xion, nH, Zsolar, &
#ifdef DO_RT
          & Np, Fp, p_gas, dNpdt, dFpdt, r%levelmin, &
#endif
          & dt_cool, 1)
     T2(1) = T2(1)*aexp**2
     aexp = aexp + daexp
     if (if_write_result) write(*,'(4(1pe10.3))')aexp,nH(1),T2_com*mu/aexp**2,n_spec(1)/nH(1)
  end do
  T2end = T2(1)/(aexp-daexp)**2
  xend(1:nion) = xion(1:nion,1)
  ne = (n_spec(3)+(n_spec(5)+2.*n_spec(6))*0.25*r%Y_He/r%X_H)

  end associate

end subroutine neq_evol_single_cell

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
FUNCTION HsurH0(z,omega0,omegaL,OmegaR)
  !-------------------------------------------------------------------------
  implicit none
  real(kind=8) :: HsurH0,z,omega0,omegaL,omegaR
  !-------------------------------------------------------------------------
  HsurH0=sqrt(Omega0*(1d0+z)**3+OmegaR*(1d0+z)**2+OmegaL)
END FUNCTION HsurH0

!=========================================================================
subroutine update_metal_cooling(tables, aexp)
  ! Compute the UV background effect on metal cooling
  ! as calibrated on Cloudy
  !-------------------------------------------------------------------------
  implicit none
  type(neq_cooling_t)::tables
  real(kind=8)::aexp
  !-------------------------------------------------------------------------
  real(kind=8),dimension(1:50),parameter::z_courty=(/                         &
       & 0.00000,0.04912,0.10060,0.15470,0.21140,0.27090,0.33330,0.39880, &
       & 0.46750,0.53960,0.61520,0.69450,0.77780,0.86510,0.95670,1.05300, &
       & 1.15400,1.25900,1.37000,1.48700,1.60900,1.73700,1.87100,2.01300, &
       & 2.16000,2.31600,2.47900,2.64900,2.82900,3.01700,3.21400,3.42100, &
       & 3.63800,3.86600,4.10500,4.35600,4.61900,4.89500,5.18400,5.48800, &
       & 5.80700,6.14100,6.49200,6.85900,7.24600,7.65000,8.07500,8.52100, &
       & 8.98900,9.50000 /)
  real(kind=8),dimension(1:50),parameter::phi_courty=(/                             &
       & 0.0499886,0.0582622,0.0678333,0.0788739,0.0915889,0.1061913,0.1229119, &
       & 0.1419961,0.1637082,0.1883230,0.2161014,0.2473183,0.2822266,0.3210551, &
       & 0.3639784,0.4111301,0.4623273,0.5172858,0.5752659,0.6351540,0.6950232, &
       & 0.7529284,0.8063160,0.8520859,0.8920522,0.9305764,0.9682031,1.0058810, &
       & 1.0444020,1.0848160,1.1282190,1.1745120,1.2226670,1.2723200,1.3231350, &
       & 1.3743020,1.4247480,1.4730590,1.5174060,1.5552610,1.5833640,1.5976390, &
       & 1.5925270,1.5613110,1.4949610,1.3813710,1.2041510,0.9403100,0.5555344, &
       & 0.0000000 /)
  real(kind=8)::ZZ,deltaZ
  integer::iZ
  !-------------------------------------------------------------------------
  ! This is a simple model to take into account the ionization background
  ! on metal cooling (calibrated using CLOUDY).
  ZZ=1d0/aexp-1d0
  iZ=1
  do while(ZZ>z_courty(iZ+1).and.iZ<49)
     iZ=iZ+1
  end do
  deltaZ=z_courty(iZ+1)-z_courty(iZ)
  ZZ=min(ZZ,z_courty(50))
  tables%phi = (phi_courty(iZ+1)*(ZZ-z_courty(iZ))/deltaZ &
       & + phi_courty(iZ)*(z_courty(iZ+1)-ZZ)/deltaZ )

end subroutine update_metal_cooling

!=========================================================================
subroutine neq_cmp_metals(tables, T2, nH, mu, metal_tot, metal_prime)
  ! Taken from the equilibrium cooling_module of RAMSES
  ! Compute cooling enhancement due to metals
  ! T2           => Temperature in Kelvin, divided by mu
  ! nH           => Hydrogen number density (H/cc)
  ! mu           => Average mass per particle in terms of mH
  ! metal_tot   <=  Metal cooling contribution to de/dt [erg s-1 cm-3]
  ! metal_prime <=  d(metal_tot)/dT2 [erg s-1 cm-3 K-1]
  !=========================================================================
  implicit none
  type(neq_cooling_t)::tables
  real(kind=8) ::T2, nH, mu, metal_tot, metal_prime
  ! Cloudy at solar metalicity
  real(kind=8),dimension(1:91),parameter :: temperature_cc07 = (/ &
       & 3.9684,4.0187,4.0690,4.1194,4.1697,4.2200,4.2703, &
       & 4.3206,4.3709,4.4212,4.4716,4.5219,4.5722,4.6225, &
       & 4.6728,4.7231,4.7734,4.8238,4.8741,4.9244,4.9747, &
       & 5.0250,5.0753,5.1256,5.1760,5.2263,5.2766,5.3269, &
       & 5.3772,5.4275,5.4778,5.5282,5.5785,5.6288,5.6791, &
       & 5.7294,5.7797,5.8300,5.8804,5.9307,5.9810,6.0313, &
       & 6.0816,6.1319,6.1822,6.2326,6.2829,6.3332,6.3835, &
       & 6.4338,6.4841,6.5345,6.5848,6.6351,6.6854,6.7357, &
       & 6.7860,6.8363,6.8867,6.9370,6.9873,7.0376,7.0879, &
       & 7.1382,7.1885,7.2388,7.2892,7.3395,7.3898,7.4401, &
       & 7.4904,7.5407,7.5911,7.6414,7.6917,7.7420,7.7923, &
       & 7.8426,7.8929,7.9433,7.9936,8.0439,8.0942,8.1445, &
       & 8.1948,8.2451,8.2955,8.3458,8.3961,8.4464,8.4967 /)
  ! Cooling from metals only (without the contribution of H and He)
  ! log cooling rate in [erg s-1 cm3]
  ! S. Ploeckinger 06/2015
  real(kind=8),dimension(1:91) :: excess_cooling_cc07 = (/ &
       &  -24.9082, -24.9082, -24.5503, -24.0898, -23.5328, -23.0696, -22.7758, &
       &  -22.6175, -22.5266, -22.4379, -22.3371, -22.2289, -22.1181, -22.0078, &
       &  -21.8992, -21.7937, -21.6921, -21.5961, -21.5089, -21.4343, -21.3765, &
       &  -21.3431, -21.3274, -21.3205, -21.3142, -21.3040, -21.2900, -21.2773, &
       &  -21.2791, -21.3181, -21.4006, -21.5045, -21.6059, -21.6676, -21.6877, &
       &  -21.6934, -21.7089, -21.7307, -21.7511, -21.7618, -21.7572, -21.7532, &
       &  -21.7668, -21.7860, -21.8129, -21.8497, -21.9035, -21.9697, -22.0497, &
       &  -22.1327, -22.2220, -22.3057, -22.3850, -22.4467, -22.4939, -22.5205, &
       &  -22.5358, -22.5391, -22.5408, -22.5408, -22.5475, -22.5589, -22.5813, &
       &  -22.6122, -22.6576, -22.7137, -22.7838, -22.8583, -22.9348, -23.0006, &
       &  -23.0547, -23.0886, -23.1101, -23.1139, -23.1147, -23.1048, -23.1017, &
       &  -23.0928, -23.0969, -23.0968, -23.1105, -23.1191, -23.1388, -23.1517, &
       &  -23.1717, -23.1837, -23.1986, -23.2058, -23.2134, -23.2139, -23.2107 /)
  real(kind=8),dimension(1:91),parameter :: excess_prime_cc07 = (/           &
       &   2.0037,  4.7267, 12.2283, 13.5820,  9.8755,  4.8379,  1.8046, &
       &   1.4574,  1.8086,  2.0685,  2.2012,  2.2250,  2.2060,  2.1605, &
       &   2.1121,  2.0335,  1.9254,  1.7861,  1.5357,  1.1784,  0.7628, &
       &   0.1500, -0.1401,  0.1272,  0.3884,  0.2761,  0.1707,  0.2279, &
       &  -0.2417, -1.7802, -3.0381, -2.3511, -0.9864, -0.0989,  0.1854, &
       &  -0.1282, -0.8028, -0.7363, -0.0093,  0.3132,  0.1894, -0.1526, &
       &  -0.3663, -0.3873, -0.3993, -0.6790, -1.0615, -1.4633, -1.5687, &
       &  -1.7183, -1.7313, -1.8324, -1.5909, -1.3199, -0.8634, -0.5542, &
       &  -0.1961, -0.0552,  0.0646, -0.0109, -0.0662, -0.2539, -0.3869, &
       &  -0.6379, -0.8404, -1.1662, -1.3930, -1.6136, -1.5706, -1.4266, &
       &  -1.0460, -0.7244, -0.3006, -0.1300,  0.1491,  0.0972,  0.2463, &
       &   0.0252,  0.1079, -0.1893, -0.1033, -0.3547, -0.2393, -0.4280, &
       &  -0.2735, -0.3670, -0.2033, -0.2261, -0.0821, -0.0754,  0.0634 /)
  real(kind=8)::TT,lTT,deltaT,lcool1,lcool2,lcool1_prime,lcool2_prime
  real(kind=8)::c1=0.4,c2=10.0,TT0=1d5,TTC=1d6,alpha1=0.15
  real(kind=8)::ux,g_courty,f_courty=1d0,g_courty_prime,f_courty_prime
  integer::iT
  !-------------------------------------------------------------------------
  TT=T2*mu
  lTT=log10(TT)
  ux=1d-4*tables%phi/nH
  g_courty=c1*(TT/TT0)**alpha1+c2*exp(-TTC/TT)
  g_courty_prime=(c1*alpha1*(TT/TT0)**alpha1+c2*exp(-TTC/TT)*TTC/TT)/TT
  f_courty=1d0/(1d0+ux/g_courty)
  f_courty_prime=ux/g_courty/(1d0+ux/g_courty)**2*g_courty_prime/g_courty

  if(lTT.ge.temperature_cc07(91))then
     metal_tot=0d0 !1d-100
     metal_prime=0d0
  else if(lTT.ge.1.0)then
     lcool1=-100d0
     lcool1_prime=0d0
      if(lTT.ge.temperature_cc07(1))then
        iT=1+int((lTT-temperature_cc07(1)) /                             &
             (temperature_cc07(91)-temperature_cc07(1))*90.0)
        iT=min(iT,90)
        iT=max(iT,1)
        deltaT = temperature_cc07(iT+1) - temperature_cc07(iT)
        lcool1 = &
             excess_cooling_cc07(iT+1)*(lTT-temperature_cc07(iT))/deltaT &
           + excess_cooling_cc07(iT)*(temperature_cc07(iT+1)-lTT)/deltaT
        lcool1_prime  =                                                  &
             excess_prime_cc07(iT+1)*(lTT-temperature_cc07(iT))/deltaT   &
           + excess_prime_cc07(iT)*(temperature_cc07(iT+1)-lTT)/deltaT
     endif
     ! Fine structure cooling from infrared lines
     lcool2=-31.522879+2.0*lTT-20.0/TT-TT*4.342944d-5
     lcool2_prime=2d0+(20d0/TT-TT*4.342944d-5)*log(10d0)
     ! Total metal cooling and temperature derivative
     metal_tot=10d0**lcool1+10d0**lcool2
     metal_prime=(10d0**lcool1*lcool1_prime+10d0**lcool2*lcool2_prime)/metal_tot
     metal_prime=metal_prime*f_courty+metal_tot*f_courty_prime
     metal_tot=metal_tot*f_courty
  else
     metal_tot=0d0 !1d-100
     metal_prime=0d0
  endif

  metal_tot = metal_tot * nH**2
  ! Convert from DlogLambda/DlogT to DLambda/DT
  metal_prime = metal_prime * metal_tot/TT * mu

end subroutine neq_cmp_metals

!*************************************************************************
FUNCTION getMu(r, xion, Tmu)
  ! Returns the mean particle mass, in units of the proton mass.
  ! xion => Hydrogen and helium ionisation fractions
  ! Tmu => T/mu in Kelvin
  !-------------------------------------------------------------------------
  implicit none
  type(run_t),intent(in)::r
  real(kind=8),intent(in)::Tmu
  real(kind=8),intent(in),dimension(nion)::xion
  real(kind=8)::getMu
  !-------------------------------------------------------------------------
  real(kind=8)::xHI, xHII, xHeII, xHeIII
  !-------------------------------------------------------------------------
  xHII=0d0 ; xHeII=0d0 ; xHeIII=0d0
  if(r%isH2) then
     xHI=xion(r%ixHI)
  else
     xHI=1.-xion(r%ixHII)
  endif
  xHII=xion(r%ixHII)
  if(r%isHe) xHeII=xion(r%ixHeII)
  if(r%isHe) xHeIII=xion(r%ixHeIII)
  getMu = 1./(r%X_H*(0.5+0.5*xHI+1.5*xHII) + 0.25*r%Y_He*(1.+xHeII+2.*xHeIII))
  if(r%is_mu_H2) getMu = getMu + exp(-1d0*(Tmu/r%Tmu_dissoc)**2) * (2.33-getMu)
END FUNCTION getMu
!************************************************************************
SUBROUTINE updateRTGroups_CoolConstants(r,tables)
  ! Update photon group cooling and heating constants, to reflect an update
  ! in rt_c_cgs and in the cross-sections and energies in the groups.
  !------------------------------------------------------------------------
  implicit none
  type(run_t)::r
  type(neq_cooling_t)::tables
#ifdef DO_RT
  !------------------------------------------------------------------------
  integer::iP, iI, i
  !------------------------------------------------------------------------
  do i=r%nlevelmax,r%levelmin,-1
    tables%signc(:,:,i) = r%group_csn*tables%rt_c_cgs(i)        ! [cm3 s-1]
    tables%sigec(:,:,i) = r%group_cse*tables%rt_c_cgs(i)        ! [cm3 s-1]
    do iP = 1,nrtgrp
      do iI = 1,nion               ! Photoheating rates for photons on ions
        tables%PHrate(iP,iI,i) =  eV2erg * &    ! See eq (19) in Aubert(08)
             (tables%sigec(iP,iI,i) * r%group_egy(iP)  &
             -tables%signc(iP,iI,i)*r%ionEvs(iI))
        tables%PHrate(iP,iI,i) = max(tables%PHrate(iP,iI,i),0d0)!Heating>0
      end do
    end do
  end do
#endif
END SUBROUTINE updateRTGroups_CoolConstants
!************************************************************************
#ifdef DO_RT
SUBROUTINE reduce_flux(Fp, cNp)
  ! Make sure the reduced photon flux is less than one
  !------------------------------------------------------------------------
  implicit none
  real(kind=8),dimension(ndim)::Fp
  real(kind=8)::cNp
  !------------------------------------------------------------------------
  real(kind=8)::fred
  fred = sqrt(sum(Fp**2))/cNp
  if(fred .gt. 1d0) Fp = Fp/fred
END SUBROUTINE reduce_flux
#endif
END MODULE neq_cooling_module
