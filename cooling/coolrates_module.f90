MODULE coolrates_module

  ! Module for returning cooling and thermochemistry interaction rates.
  ! The temperature dependence is tabulated in rate (nonlog) versus the
  ! log of T, because it is expensive to calculate on the fly.
  ! The rates are interpolated using cubic splines, and extrapolated in
  ! log-log space if temperature is above table boundaries.
  ! Joki Rosdahl and Andreas Bleuler, September 2015.
  ! Nickerson, Teyssier, and Rosdahl, March 2018.

  use constants
  use amr_parameters, only: MAXLEVEL
  use hydro_parameters, only: nion
  use rt_parameters, only: nrtgrp
  use amr_commons, only: run_t
  implicit none

  private ! default
  public init_coolrates_tables, update_coolrates_tables, inp_coolrates_table, &
       & compCoolrate, neq_cooling_t, update_rt_c

  ! Default cooling rates table parameters
  integer, parameter :: nbinT = 1001
  real(kind=8), parameter :: Tmin = 1d-2
  real(kind=8), parameter :: Tmax = 1d+9

  type coolrates_table
     ! Cooling and interaction rates (log):
     real(kind=8),dimension(nbinT) :: rates = 0d0
     ! Temperature derivatives of those rates (drate/dlog(T)):
     real(kind=8),dimension(nbinT) :: primes = 0d0
  end type coolrates_table

  type neq_cooling_t

     ! Reduced speed of light in cgs
     real(kind=8),dimension(1:MAXLEVEL)::rt_c_cgs

     ! Contribution of UV background to metal cooling
     real(kind=8) :: phi

     ! UV background heating and ionization rates
     real(kind=8),dimension(nion,2) :: UVrates
#ifdef DO_RT
     ! Contribution of each group to photo-ionization rates
#ifdef DO_RTZ
     real(kind=8),dimension(nrtgrp,1:27,1:27,MAXLEVEL) :: signc

     ! Contribution of each group to photo-heating rates
     real(kind=8),dimension(nrtgrp,1:27,1:27,MAXLEVEL) :: sigec, PHrate
#else
     real(kind=8),dimension(nrtgrp,nion,MAXLEVEL) :: signc

     ! Contribution of each group to photo-heating rates
     real(kind=8),dimension(nrtgrp,nion,MAXLEVEL) :: sigec, PHrate
#endif
#endif
     real(kind=8) :: dlogTinv ! Inverse of the bin space (in K)
     real(kind=8) :: hTable, h2Table, h3Table   ! Interpol constants
     real(kind=8) :: one_over_lnTen, one_over_hTable, one_over_h2Table
     real(kind=8) :: three_over_h2Table, two_over_h3Table

     ! Lookup temperature in log K
     real(kind=8),dimension(nbinT) :: T_lookup = 0d0

     ! Interpolation factors saved for later use
     integer :: iT = 1
     logical :: extrap
     real(kind=8) :: facT, yy, yy2, yy3, Tlast=-1

     type(coolrates_table) :: tbl_alphaZ_H2  ! H2 dust grain formation
     type(coolrates_table) :: tbl_alphaGP_H2 ! H2 gas-phase formation
     type(coolrates_table) :: tbl_alphaA_HII ! Case A rec. coefficients
     type(coolrates_table) :: tbl_alphaA_HeII
     type(coolrates_table) :: tbl_alphaA_HeIII
     type(coolrates_table) :: tbl_alphaB_HII ! Case B rec. coefficients
     type(coolrates_table) :: tbl_alphaB_HeII
     type(coolrates_table) :: tbl_alphaB_HeIII
     type(coolrates_table) :: tbl_beta_HI ! Collisional ionisation rates
     type(coolrates_table) :: tbl_beta_HeI
     type(coolrates_table) :: tbl_beta_HeII
     type(coolrates_table) :: tbl_cr_ci_HI ! Coll. ionisation cooling
     type(coolrates_table) :: tbl_cr_ci_HeI
     type(coolrates_table) :: tbl_cr_ci_HeII
     type(coolrates_table) :: tbl_cr_ce_HI ! Coll. excitation cooling
     type(coolrates_table) :: tbl_cr_ce_HeI
     type(coolrates_table) :: tbl_cr_ce_HeII
     type(coolrates_table) :: tbl_cr_r_HII ! Recombination cooling
     type(coolrates_table) :: tbl_cr_r_HeII
     type(coolrates_table) :: tbl_cr_r_HeIII
     type(coolrates_table) :: tbl_cr_bre ! Bremsstrahlung cooling rates
     type(coolrates_table) :: tbl_cr_com ! Compton cooling rates
     type(coolrates_table) :: tbl_cr_die ! Dielectronic cooling rates
     type(coolrates_table) :: tbl_beta_H2HI ! Coll. dissociation rates
     type(coolrates_table) :: tbl_beta_H2H2 ! Coll. dissociation rates
     type(coolrates_table) :: tbl_beta_H3B ! Coll. dissociation rates
     type(coolrates_table) :: tbl_cr_H2HI ! Collisional diss. cooling
     type(coolrates_table) :: tbl_cr_H2H2 ! Collisional diss. cooling
     type(coolrates_table) :: tbl_beta_ly_red ! Lyman Alpha redshifting
     type(coolrates_table) :: tbl_beta_cmb_HI ! CMB ion. from excited states
     type(coolrates_table) :: tbl_gamma_cmb_HI ! Total CMB ion. rate
     type(coolrates_table) :: tbl_beta_cmb_Hm ! H- destr. by CMB over H2 cr.

  end type neq_cooling_t

CONTAINS

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE init_coolrates_tables(r, tables, aexp, Ha)

  ! Initialise the cooling rates tables.
  !-------------------------------------------------------------------------
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(run_t)::r
  type(neq_cooling_t)::tables
  real(kind=8)::aexp, Ha
  !-------------------------------------------------------------------------
#ifndef WITHOUTMPI
  integer::ierr
#endif
  integer::myid, ncpu, iT
  !-------------------------------------------------------------------------
#ifndef WITHOUTMPI
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu,ierr)
#endif
#ifdef WITHOUTMPI
  myid=0
  ncpu=1
#endif

  ! Initialise the table lookup temperatures -----------------------------
  do iT=1, nbinT
     tables%T_lookup(iT) =  log10(Tmin) + (dble(iT)-1d0) / (dble(nbinT)-1d0)     &
                         * (log10(Tmax) - log10(Tmin))
  end do
  tables%dlogTinv = dble(nbinT-1)/(tables%T_lookup(nbinT)-tables%T_lookup(1)) ! (space)^-1
  tables%hTable = 1d0/tables%dlogTinv                                         !
  tables%h2Table = tables%hTable*tables%hTable                                ! Constants for table
  tables%h3Table = tables%h2Table*tables%hTable                               ! interpolation

  tables%one_over_lnTen = 1d0/log(10d0)
  tables%one_over_hTable = 1d0/tables%hTable
  tables%one_over_h2Table = 1d0/tables%h2Table
  tables%three_over_h2Table = 3d0/tables%h2Table
  tables%two_over_h3Table = 2d0/tables%h3Table

  do iT = myid+1, nbinT, ncpu ! Loop over TK and assign rates
     call comp_table_rates(r, tables, iT)
     ! Compton rates at aexp (comp_table_rates leaves them out)
     call update_table_rates(r, tables, iT, aexp, Ha)
  end do ! end TK loop

  ! Distribute the complete table between cpus ---------------------------
#ifndef WITHOUTMPI
  call mpi_distribute_coolrates_table(tables%tbl_alphaZ_H2)
  call mpi_distribute_coolrates_table(tables%tbl_alphaGP_H2)
  call mpi_distribute_coolrates_table(tables%tbl_alphaA_HII)
  call mpi_distribute_coolrates_table(tables%tbl_alphaA_HeII)
  call mpi_distribute_coolrates_table(tables%tbl_alphaA_HeIII)

  call mpi_distribute_coolrates_table(tables%tbl_alphaB_HII)
  call mpi_distribute_coolrates_table(tables%tbl_alphaB_HeII)
  call mpi_distribute_coolrates_table(tables%tbl_alphaB_HeIII)

  call mpi_distribute_coolrates_table(tables%tbl_beta_HI)
  call mpi_distribute_coolrates_table(tables%tbl_beta_HeI)
  call mpi_distribute_coolrates_table(tables%tbl_beta_HeII)

  call mpi_distribute_coolrates_table(tables%tbl_cr_ci_HI)
  call mpi_distribute_coolrates_table(tables%tbl_cr_ci_HeI)
  call mpi_distribute_coolrates_table(tables%tbl_cr_ci_HeII)

  call mpi_distribute_coolrates_table(tables%tbl_cr_ce_HI)
  call mpi_distribute_coolrates_table(tables%tbl_cr_ce_HeI)
  call mpi_distribute_coolrates_table(tables%tbl_cr_ce_HeII)

  call mpi_distribute_coolrates_table(tables%tbl_cr_r_HII)
  call mpi_distribute_coolrates_table(tables%tbl_cr_r_HeII)
  call mpi_distribute_coolrates_table(tables%tbl_cr_r_HeIII)

  call mpi_distribute_coolrates_table(tables%tbl_cr_bre)
  call mpi_distribute_coolrates_table(tables%tbl_cr_com)
  call mpi_distribute_coolrates_table(tables%tbl_cr_die)

  call mpi_distribute_coolrates_table(tables%tbl_beta_H2HI)
  call mpi_distribute_coolrates_table(tables%tbl_beta_H2H2)
  call mpi_distribute_coolrates_table(tables%tbl_beta_H3B)
  call mpi_distribute_coolrates_table(tables%tbl_cr_H2HI)
  call mpi_distribute_coolrates_table(tables%tbl_cr_H2H2)

  call mpi_distribute_coolrates_table(tables%tbl_beta_ly_red)
  call mpi_distribute_coolrates_table(tables%tbl_beta_cmb_HI)
  call mpi_distribute_coolrates_table(tables%tbl_gamma_cmb_HI)
  call mpi_distribute_coolrates_table(tables%tbl_beta_cmb_Hm)

#endif

  if(myid==0) print*,'Coolrates tables initialised '

END SUBROUTINE init_coolrates_tables

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE update_coolrates_tables(r, tables, aexp, Ha)
  ! Update cooling rates lookup tables which depend on aexp
  !-------------------------------------------------------------------------
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(run_t)::r
  type(neq_cooling_t)::tables
  real(kind=8)::aexp, Ha
  !-------------------------------------------------------------------------
#ifndef WITHOUTMPI
  integer::ierr
#endif
  integer::myid, ncpu, iT
  !-------------------------------------------------------------------------
#ifndef WITHOUTMPI
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu,ierr)
#endif
#ifdef WITHOUTMPI
  myid=0
  ncpu=1
#endif
  tables%tbl_cr_com%rates  = 0d0 ; tables%tbl_cr_com%primes = 0d0
  tables%tbl_beta_ly_red%rates  = 0d0 ; tables%tbl_beta_ly_red%primes = 0d0
  tables%tbl_beta_cmb_HI%rates  = 0d0 ; tables%tbl_beta_cmb_HI%primes = 0d0
  tables%tbl_gamma_cmb_HI%rates  = 0d0 ; tables%tbl_gamma_cmb_HI%primes = 0d0
  tables%tbl_beta_cmb_Hm%rates  = 0d0 ; tables%tbl_beta_cmb_Hm%primes = 0d0
  do iT = myid+1, nbinT, ncpu ! Loop over TK and assign rates
     call update_table_rates(r, tables, iT, aexp, Ha)
  end do ! end TK loop

  ! Distribute the complete table between cpus ---------------------------
#ifndef WITHOUTMPI
  call mpi_distribute_coolrates_table(tables%tbl_cr_com)
  call mpi_distribute_coolrates_table(tables%tbl_beta_ly_red)
  call mpi_distribute_coolrates_table(tables%tbl_beta_cmb_HI)
  call mpi_distribute_coolrates_table(tables%tbl_gamma_cmb_HI)
  call mpi_distribute_coolrates_table(tables%tbl_beta_cmb_Hm)
#endif

END SUBROUTINE update_coolrates_tables

#ifndef WITHOUTMPI
!PRIVATEXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE mpi_distribute_coolrates_table(table)
  ! Distribute table between all cpus, assuming table contains only partial
  ! entries on each cpu, but the whole table is acquired by summing those
  ! partial tables
  !-------------------------------------------------------------------------
  use mpi
  implicit none
  type(coolrates_table) :: table
  real(kind=8),dimension(:),allocatable :: table_mpi_sum
  integer::ierr
  !-------------------------------------------------------------------------
  allocate(table_mpi_sum(nbinT))
  call MPI_ALLREDUCE(table%rates,table_mpi_sum                           &
                  ,nbinT,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  table%rates = table_mpi_sum
  call MPI_ALLREDUCE(table%primes,table_mpi_sum                          &
                  ,nbinT,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  table%primes = table_mpi_sum
  deallocate(table_mpi_sum)

END SUBROUTINE mpi_distribute_coolrates_table
#endif

!PRIVATEXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE comp_table_rates(r, tables, iT)
  ! Fill in index iTK in all rates tables.
  ! Note that Compton cooling/heating is done elsewhere.
  !-------------------------------------------------------------------------
  implicit none
  type(run_t)::r
  type(neq_cooling_t)::tables
  integer::iT
  !-------------------------------------------------------------------------
  real(kind=8)::T, T2, T5, lambda, f, hf, laHII, laHeII, laHeIII
  real(kind=8)::lowrleft,lowrright,lowr_hi,lowr_h2
  real(kind=8)::lowvleft_hi,lowvright_hi,lowv_hi,lowv_h2
  real(kind=8)::TT
  !-------------------------------------------------------------------------
  ! Rates are stored in non-log, while temperature derivatives (primes)
  ! are stored in dRate/dlogT (= dRate/dT * T * ln(10)).

  ! The log-log primes are just the normal primes times T/rate,
  ! i.e. dlogL/dlogT = T/L dL/dT

  T = 10d0**tables%T_lookup(iT)

  ! Creation rate of H2 on dust [cm3 s-1] (functio of T from--------------
  ! Hollenback & McKee 1979, average rate between Habart 2004 dense PDRs
  ! and Jura 1974 diffuse ISM)
  ! When using, must multiply by Z*f_dust
  T2 = T/1d2
  f = (1.0+0.4*sqrt(T2)+0.2*T2+0.08*T2**2)
  tables%tbl_alphaZ_H2%rates(iT)   = 9.0d-17 * sqrt(T2) / f
  tables%tbl_alphaZ_H2%primes(iT)  = (-0.12*T2**2 - 0.1*T2 + 0.5) / f           &
                            * log(10d0) * tables%tbl_alphaZ_H2%rates(iT)

  ! Gas phase formation rate rate for H2 on H- [cm3 s-1] assuming---------
  ! equilibrium abundances for H-, as explained in the Appendix of
  ! McKee and Krumholz (2010)
  tables%tbl_alphaGP_H2%rates(iT)  = 1.83d-18*T**0.88
  tables%tbl_alphaGP_H2%primes(iT) = 0.88                                       &
                            * log(10d0) * tables%tbl_alphaGP_H2%rates(iT)
  ! Case A rec. coefficient [cm3 s-1] for HII (Hui&Gnedin'97)-------------
  lambda = 315614./T                                ! 2d0 * 157807d0 / T
  f = 1d0+(lambda/0.522)**0.47
  tables%tbl_alphaA_HII%rates(iT)  =  1.269d-13 * lambda**1.503 / f**1.923
  tables%tbl_alphaA_HII%primes(iT) = ( 0.90381*(f-1.)/f - 1.503 )               &
                            * log(10d0) * tables%tbl_alphaA_HII%rates(iT)

  ! Case A rec. coefficient [cm3 s-1] for HeII (Hui&Gnedin'97)------------
  lambda = 570670./T
  tables%tbl_alphaA_HeII%rates(iT)  = 3d-14 * lambda**0.654
  tables%tbl_alphaA_HeII%primes(iT) = -0.654                                    &
                             * log(10d0) * tables%tbl_alphaA_HeII%rates(iT)

  ! Case A rec. coefficient [cm3 s-1] for HeIII (Hui&Gnedin'97)-----------
  lambda =  1263030./T
  f= 1d0+(lambda/0.522)**0.47
  tables%tbl_alphaA_HeIII%rates(iT)  =  2.538d-13 * lambda**1.503 / f**1.923
  tables%tbl_alphaA_HeIII%primes(iT) =  ( 0.90381*(f-1.)/f - 1.503 )            &
                              * log(10d0) * tables%tbl_alphaA_HeIII%rates(iT)

  ! Case B rec. coefficient [cm3 s-1] for HII (Hui&Gnedin'97)-------------
  lambda = 315614./T
  f= 1d0+(lambda/2.74)**0.407
  tables%tbl_alphaB_HII%rates(iT)  = 2.753d-14 * lambda**1.5 / f**2.242
  tables%tbl_alphaB_HII%primes(iT) = ( 0.912494*(f-1.)/f - 1.5 )                &
                             * log(10d0) * tables%tbl_alphaB_HII%rates(iT)

  ! Case B rec. coefficient [cm3 s-1] for HeII (Hui&Gnedin'97)------------
  lambda = 570670./T
  tables%tbl_alphaB_HeII%rates(iT)  = 1.26d-14 * lambda**0.75
  tables%tbl_alphaB_HeII%primes(iT) =  -0.75                                    &
                             * log(10d0) * tables%tbl_alphaB_HeII%rates(iT)

  ! Case B rec. coefficient [cm3 s-1] for HeIII (Hui&Gnedin'97)-----------
  lambda = 1263030./T
  f= 1d0+(lambda/2.74)**0.407
  tables%tbl_alphaB_HeIII%rates(iT)  = 5.506d-14 * lambda**1.5 / f**2.242
  tables%tbl_alphaB_HeIII%primes(iT) = ( 0.912494*(f-1.)/f - 1.5 )              &
                              * log(10d0) * tables%tbl_alphaB_HeIII%rates(iT)

  ! Collisional ionization rate [cm3 s-1] of HI (Maselli&'03)-------------
  T5 = T/1d5
  f = 1d0+sqrt(T5) ; hf=0.5d0/f
  tables%tbl_beta_HI%rates(iT)  = 5.85d-11 * sqrt(T) / f * exp(-157809.1d0/T)
  tables%tbl_beta_HI%primes(iT) = (hf+157809.1d0/T)                             &
                         * log(10d0) * tables%tbl_beta_HI%rates(iT)

  ! Collisional ionization rate [cm3 s-1] of HeI (Maselli&'03)------------
  tables%tbl_beta_HeI%rates(iT)  = 2.38d-11 * sqrt(T) / f * exp(-285335.4d0/T)
  tables%tbl_beta_HeI%primes(iT) = (hf+285335.4d0/T)                            &
                          * log(10d0) * tables%tbl_beta_HeI%rates(iT)

  ! Collisional ionization rate [cm3 s-1] of HeII (Maselli&'03)-----------
  tables%tbl_beta_HeII%rates(iT)  = 5.68d-12 * sqrt(T) / f * exp(-631515d0/T)
  tables%tbl_beta_HeII%primes(iT) = (hf+631515d0/T)                            &
                           * log(10d0) * tables%tbl_beta_HeII%rates(iT)

  ! BEGIN COOLING RATES-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  T5 = T/1d5
  f = 1d0+sqrt(T5) ; hf=0.5d0/f

  ! Coll. Ionization Cooling from Cen 1992 (via Maselli et al 2003)
  tables%tbl_cr_ci_HI%rates(iT)  = 1.27d-21 * sqrt(T) / f * exp(-157809.1/T)
  tables%tbl_cr_ci_HI%primes(iT) = (hf+157809.1/T)                              &
                          * log(10d0) * tables%tbl_cr_ci_HI%rates(iT)

  tables%tbl_cr_ci_HeI%rates(iT)  = 9.38d-22 * sqrt(T) / f * exp(-285335.4/T)
  tables%tbl_cr_ci_HeI%primes(iT) = (hf+285335.4/T)                             &
                           * log(10d0) * tables%tbl_cr_ci_HeI%rates(iT)

  tables%tbl_cr_ci_HeII%rates(iT)  = 4.95d-22 * sqrt(T) / f * exp(-631515. /T)
  tables%tbl_cr_ci_HeII%primes(iT) = (hf+631515.0/T)                            &
                            * log(10d0) * tables%tbl_cr_ci_HeII%rates(iT)

  ! Collisional excitation cooling from Cen'92
  tables%tbl_cr_ce_HI%rates(iT)  = 7.5d-19 / f * exp(-118348./T)
  tables%tbl_cr_ce_HI%primes(iT) = (118348./T - 0.5d0 * sqrt(T5) / f )          &
                          * log(10d0) * tables%tbl_cr_ce_HI%rates(iT)

  tables%tbl_cr_ce_HeI%rates(iT)  = 9.10d-27 * T**(-0.1687) / f * exp(-13179./T)
  tables%tbl_cr_ce_HeI%primes(iT) = (13179./T - 0.1687 - 0.5d0 * sqrt(T5) / f)  &
                           * log(10d0) * tables%tbl_cr_ce_HeI%rates(iT)

  tables%tbl_cr_ce_HeII%rates(iT) = 5.54d-17 * T**(-0.397)  / f * exp(-473638./T)
  tables%tbl_cr_ce_HeII%primes(iT) = (473638./T - 0.397 - 0.5d0 * sqrt(T5) / f) &
                            * log(10d0) * tables%tbl_cr_ce_HeII%rates(iT)

  ! Recombination Cooling (Hui&Gnedin'97)
  laHII    = 315614./T
  laHeII   = 570670./T
  laHeIII  = 1263030./T
  if(.not. r%rt_otsa) then ! Case A
     f = 1d0+(laHII/0.541)**0.502
     tables%tbl_cr_r_HII%rates(iT)    = 1.778d-29 * laHII**1.965 / f**2.697 * T
     tables%tbl_cr_r_HII%primes(iT)   = (-0.965 + 1.35389*(f-1.)/f)             &
                               * log(10d0) * tables%tbl_cr_r_HII%rates(iT)

     tables%tbl_cr_r_HeII%rates(iT)   = 3d-14 * laHeII**0.654 * kB * T
     tables%tbl_cr_r_HeII%primes(iT)  = 0.346                                   &
                               * log(10d0) * tables%tbl_cr_r_HeII%rates(iT)

     f = 1d0+(laHeIII/0.541)**0.502
     tables%tbl_cr_r_HeIII%rates(iT) = 14.224d-29 * laHeIII**1.965 / f**2.697 * T
     tables%tbl_cr_r_HeIII%primes(iT)= (-0.965 + 1.35389*(f-1.)/f)              &
                              * log(10d0) * tables%tbl_cr_r_HeIII%rates(iT)
  else ! Case B
     f = 1d0+(laHII/2.25)**0.376
     tables%tbl_cr_r_HII%rates(iT)    = 3.435d-30 * laHII**1.97 / f**3.72 * T
     tables%tbl_cr_r_HII%primes(iT)   = (-0.97 + 1.39827*(f-1.)/f)              &
                               * log(10d0) * tables%tbl_cr_r_HII%rates(iT)

     tables%tbl_cr_r_HeII%rates(iT)   = 1.26d-14 * laHeII**0.75 * kB * T
     tables%tbl_cr_r_HeII%primes(iT)  = 0.25                                    &
                               * log(10d0) * tables%tbl_cr_r_HeII%rates(iT)

     f = 1d0+(laHeIII/2.25)**0.376
     tables%tbl_cr_r_HeIII%rates(iT)  = 27.48d-30 * laHeIII**1.97 / f**3.72 * T
     tables%tbl_cr_r_HeIII%primes(iT) = (-0.97 + 1.39827*(f-1.)/f)              &
                               * log(10d0) * tables%tbl_cr_r_HeIII%rates(iT)
  endif

  ! Bremsstrahlung from Osterbrock & Ferland 2006
  tables%tbl_cr_bre%rates(iT)  = 1.42d-27 * 1.5 * sqrt(T)
  tables%tbl_cr_bre%primes(iT) = 0.5                                            &
                        * log(10d0) * tables%tbl_cr_bre%rates(iT)

  ! Dielectronic recombination cooling, from Black 1981
  f = 1.24d-13*T**(-1.5d0)*exp(-470000d0/T)
  tables%tbl_cr_die%rates(iT) = f*(1d0+0.3d0*exp(-94000d0/T))
  tables%tbl_cr_die%primes(iT)=0d0
  if(tables%tbl_cr_die%rates(iT) .gt. 0d0) then ! Can simplify w algebra
     tables%tbl_cr_die%primes(iT) = (tables%tbl_cr_die%rates(iT)*(564000.-1.5*T)       &
                              - f*94000.) /T**2 * T * log(10d0)
  endif

  ! BEGIN H2 STUFF --------------------------------------------------------

  ! Collisional dissociation ground state (Dove&Mandy 1986) [cm3 s-1]
  tables%tbl_Beta_H2HI%rates(iT) =                                              &
       7.073d-19*(T**2.012)*exp(-5.179d4/T)/(1d0+2.130d-5*T)**3.512
  ! Col dissociation ground state (Martin&Keogh&Mandy 1998)  [cm3 s-1]
  tables%tbl_Beta_H2H2%rates(iT) = &
       5.996d-30*(T**4.1881)*exp(-5.466d4/T)/(1d0+6.761d-6*T)**5.6881
  ! Three body H2 formation 3H -> H2 + H (Forrey 2013) [cm6 s-1]
  ! For H2+2H -> 2H2 use Beta_H3B/8 (Palla 1983)
  tables%tbl_Beta_H3B%rates(iT) = &
       6d-32*T**(-0.25)+2d-31*T**(-0.5)
  ! Prefer to use numerical prime, rather than analytic:
  TT=T*1.0001
  tables%tbl_Beta_H2HI%primes(iT) = &
       7.073d-19*(TT**2.012)*exp(-5.179d4/TT)/(1d0+2.130d-5*TT)**3.512
  tables%tbl_Beta_H2HI%primes(iT) = &
       (tables%tbl_Beta_H2HI%primes(iT) - tables%tbl_Beta_H2HI%rates(iT)) &
       / (0.0001*T) * T * log(10d0) ! last two terms for log-log
  tables%tbl_Beta_H2H2%primes(iT) = &
       5.996d-30*(TT**4.1881)*exp(-5.466d4/TT)/(1d0+6.761d-6*TT)**5.6881
  tables%tbl_Beta_H2H2%primes(iT) = &
       (tables%tbl_Beta_H2H2%primes(iT) - tables%tbl_Beta_H2H2%rates(iT)) &
       / (0.0001*T) * T * log(10d0) ! last two terms for log-log
  tables%tbl_Beta_H3B%primes(iT) = &
       6d-32*TT**(-0.25)+2d-31*TT**(-0.5)
  tables%tbl_Beta_H3B%primes(iT) = &
       (tables%tbl_Beta_H3B%primes(iT) - tables%tbl_Beta_H3B%rates(iT)) &
       / (0.0001*T) * T * log(10d0) ! last two terms for log-log

  ! Collisional H2 cooling************************************************
  ! Hallenbach McKee (1979) + Halle Combes (2012) in cgs
  lowrleft=1.25*exp(-1.70d2/T)*2.35d-14
  lowrright=1.75*exp(-5.05d2/T)*6.97d-14
  lowr_hi=lowrleft*gamma_hi(T,2d0)+lowrright*gamma_hi(T,3d0)
  lowr_h2=lowrleft*gamma_h2(T,2d0)+lowrright*gamma_h2(T,3d0)
  lowvleft_hi=exp(-5860.0/T)*8.09d-13*1.0d-12*sqrt(T)*exp(-1.0d3/T)
  lowvright_hi=exp(-11720.0/T)*1.6d-12*1.6d-12*sqrt(T)*exp(-1.0*(4.0d2/T)**2)
  lowv_hi=lowvleft_hi+lowvright_hi
  lowv_h2=exp(-5860.0/T)*8.09d-13*1.4d-12*sqrt(T)*exp(-1.2d4/(T+1.2d3))

  tables%tbl_cr_H2HI%rates(iT)=lowr_hi+lowv_hi
  tables%tbl_cr_H2H2%rates(iT)=lowr_h2+lowv_h2

  ! Collisional H2 cooling derivatives (calculated numerically)
  TT=T*1.0001
  lowrleft=1.25*exp(-1.70d2/TT)*2.35d-14
  lowrright=1.75*exp(-5.05d2/TT)*6.97d-14
  lowr_hi=lowrleft*gamma_hi(TT,2d0)+lowrright*gamma_hi(TT,3d0)
  lowr_h2=lowrleft*gamma_h2(TT,2d0)+lowrright*gamma_h2(TT,3d0)
  lowvleft_hi=exp(-5860.0/TT)*8.09d-13*1.0d-12*sqrt(TT)*exp(-1.0d3/TT)
  lowvright_hi=exp(-11720.0/TT)*1.6d-12*1.6d-12*sqrt(TT)*exp(-1.0*(4.0d2/TT)**2)
  lowv_hi=lowvleft_hi+lowvright_hi
  lowv_h2=exp(-5860.0/TT)*8.09d-13*1.4d-12*sqrt(TT)*exp(-1.2d4/(TT+1.2d3))

  tables%tbl_cr_H2HI%primes(iT)=lowr_hi+lowv_hi
  tables%tbl_cr_H2HI%primes(iT) = &
       (tables%tbl_cr_H2HI%primes(iT) - tables%tbl_cr_H2HI%rates(iT)) &
       / (0.0001*T) * T * log(10d0) ! last two terms for log-log

  tables%tbl_cr_H2H2%primes(iT)=lowr_h2+lowv_h2
  tables%tbl_cr_H2H2%primes(iT) = &
       (tables%tbl_cr_H2H2%primes(iT) - tables%tbl_cr_H2H2%rates(iT)) &
       / (0.0001*T) * T * log(10d0) ! last two terms for log-log

  ! END H2 STUFF ----------------------------------------------------------

END SUBROUTINE comp_table_rates

!PRIVATEXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE update_rt_c(r, g, tables)
  ! Update reduced speed of light in bith code units and cgs units.
  !-------------------------------------------------------------------------
  use amr_commons, only: run_t, global_t
  implicit none
  type(run_t) :: r
  type(global_t) :: g
  type(neq_cooling_t) :: tables
  real(kind=8)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  integer::i
  !-------------------------------------------------------------------------
  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  do i=r%nlevelmax,r%levelmin,-1
    g%rt_c_cgs(i) = c_cgs * r%rt_c_fraction(i)
    g%rt_c(i) = g%rt_c_cgs(i) / scale_v
    tables%rt_c_cgs(i) = g%rt_c_cgs(i)
  enddo

END SUBROUTINE update_rt_c

!PRIVATEXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
SUBROUTINE update_table_rates(r, tables, iT, aexp, Ha)
  ! Update index iTK in all tables due to change in aexp
  !-------------------------------------------------------------------------
  implicit none
  type(run_t)::r
  type(neq_cooling_t)::tables
  integer::iT
  real(kind=8)::aexp, Ha
  !-------------------------------------------------------------------------
  real(kind=8)::T, Ta, saha_factor, lyman_exp
  !-------------------------------------------------------------------------
  ! Rates are stored in log, while temperature derivatives (primes) are
  ! stored in non-log

  ! Gas temperature
  T  = 10d0**tables%T_lookup(iT)
  ! CMB temperature
  Ta = 2.725d0/aexp

  ! Compton Cooling from Haimann et al. 96, via Maselli et al.
  tables%tbl_cr_com%rates(iT)   = 1.017d-37 * Ta**4 * (T-Ta)
  tables%tbl_cr_com%primes(iT)  = T / (T-Ta)                                    &
                         * log(10d0) * tables%tbl_cr_com%rates(iT)

  ! Recombination physics from Peebles (1993)

  ! Lyman Alpha photons redshifted out of resonance (eq. 6.107)
  tables%tbl_beta_ly_red%rates(iT) = (1.21567d-5)**3 / (8d0*pi) / Ha

  ! Ionization of excited states (eq. 6.110) with fudge factor (1.14d0) from Seager et al. (1999)
  saha_factor = 1.14d0 * 2.4147d15*Ta**1.5d0 * exp(-39451.75d0/Ta)
  tables%tbl_beta_cmb_HI%rates(iT) = tables%tbl_AlphaB_HII%rates(iT) * saha_factor
  tables%tbl_beta_cmb_HI%primes(iT) = tables%tbl_AlphaB_HII%primes(iT) * saha_factor

  ! CMB contribution to ionization (eq. 6.112)
  lyman_exp = exp(-118355.25d0/Ta)
  tables%tbl_gamma_cmb_HI%rates(iT) = tables%tbl_beta_cmb_HI%rates(iT) * lyman_exp
  tables%tbl_gamma_cmb_HI%primes(iT) = tables%tbl_beta_cmb_HI%primes(iT) * lyman_exp

  ! H- photodissociation rate by the CMB from Galli and Palla (1998) over
  ! the H2 formation rate by Schmeltekopf et al. (1967)
  tables%tbl_beta_cmb_Hm%rates(iT) = 0.11d0*Ta**2.13d0*exp(-8823d0/Ta) / 1.3d-9

END SUBROUTINE update_table_rates

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
FUNCTION inp_coolrates_table(tables, interp, T, allow_negative, retPrime)
  ! Returns TABULATED rate value from given table
  ! tables       => Object containing all the tables and constants
  ! interp       => Rates (and primes) table to interpolate
  ! T            => Temperature [K]
  ! retPrime     <= Temperature derivative of rate at T (optional).
  !-------------------------------------------------------------------------
  implicit none
  type(neq_cooling_t) :: tables
  type(coolrates_table) :: interp
  real(kind=8),intent(in) :: T
  logical,intent(in) :: allow_negative
  real(kind=8),optional :: retPrime
  real(kind=8) :: inp_coolrates_table
  !-------------------------------------------------------------------------
  real(kind=8) :: fa, fb, fprimea, fprimeb
  real(kind=8) :: alpha, beta, gamma
  !-------------------------------------------------------------------------
  if (.not. (T .eq. tables%Tlast)) then    ! Reuse index if same T from last call
     ! Log of T, snapped to table at the lower boundary, but allowed
     ! to go above upper boundary, in which case we use extrapolation:
     tables%facT = MAX( log10(T), tables%T_lookup(1) )
     tables%extrap=.false.
     if(tables%facT .gt. tables%T_lookup(nbinT)) tables%extrap=.true. ! Above upper limit
     ! Lower closest index in table:
     tables%iT = MIN(MAX(int((tables%facT-tables%T_lookup(1))*tables%dlogTinv)+1, 1), nbinT-1)
     tables%yy=tables%facT-tables%T_lookup(tables%iT)  ! Dist., in log(T), from T to lower table index
     tables%yy2=tables%yy*tables%yy             ! That distance squared
     tables%yy3=tables%yy2*tables%yy            ! ...and cubed
     tables%Tlast = T
  endif

  if(tables%extrap) then ! TK above upper table limit, so extrapolate in log-log:
     alpha = (log10(interp%rates(nbinT))             &
              -log10(interp%rates(nbinT-1)) )        &
           / (tables%T_lookup(nbinT)-tables%T_lookup(nbinT-1))
     inp_coolrates_table = 10d0**(log10(interp%rates(nbinT))        &
                                  + alpha * (tables%facT - tables%T_lookup(nbinT)))
     if( present(retPrime) )                                             &
          retPrime = alpha * inp_coolrates_table / T
     return
  endif

  fa = interp%rates(tables%iT)           !
  fb = interp%rates(tables%iT+1)         !  Values at neighbouring table
  fprimea = interp%primes(tables%iT)     !  indexes
  fprimeb = interp%primes(tables%iT+1)   !

  ! Spline interpolation:
  alpha = fprimea
  beta = (fb-fa) * tables%three_over_h2Table - (2d0*fprimea+fprimeb) * tables%one_over_hTable
  gamma = (fprimea+fprimeb) * tables%one_over_h2Table - (fb-fa) * tables%two_over_h3Table
  inp_coolrates_table = fa + alpha * tables%yy + beta * tables%yy2 + gamma * tables%yy3
  if(.not. allow_negative) then
     ! Only allow positive rates (spline can go negative):
     inp_coolrates_table = max(0d0,inp_coolrates_table)
  endif
  if( present(retPrime) )                                                &
       retPrime = (alpha+2d0*beta*tables%yy+3d0*gamma*tables%yy2) / T * tables%one_over_lnTen
END FUNCTION inp_coolrates_table

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
FUNCTION compCoolrate(r, tables, T, ne, nN, nI, dcooldT)

  ! Compute cooling rate in a cell, using interpolation from cooling rate
  ! tables
  ! T        => Cell emperature [K]
  ! nN       => Neutral and molecular abundances for each ionisation species
  ! nI       => Ionized and dissociated abundances for each species
  ! dcooldT  <= Temperature derivative of the rate
  ! dcooldx  <= Ionized fraction derivative of the rate
  ! returns:    Cooling rate [erg s-1 cm-3]
  !-------------------------------------------------------------------------
  implicit none
  type(run_t)::r
  type(neq_cooling_t)::tables
  real(kind=8)::T, ne
  real(kind=8),dimension(nion)::nN, nI
  real(kind=8)::compCoolrate, dcooldT
  !-------------------------------------------------------------------------
  real(kind=8)::ci_HI, ci_HeI, ci_HeII
  real(kind=8)::ci_HI_prime, ci_HeI_prime, ci_HeII_prime
  real(kind=8)::cr_H2HI, cr_H2H2, cr_H2HI_prime, cr_H2H2_prime
  real(kind=8)::ce_HI, ce_HeI, ce_HeII
  real(kind=8)::ce_HI_prime, ce_HeI_prime, ce_HeII_prime
  real(kind=8)::r_HII, r_HeII, r_HeIII
  real(kind=8)::r_HII_prime, r_HeII_prime, r_HeIII_prime
  real(kind=8)::bre, brefac, bre_prime, com, com_prime, die, die_prime
!-------------------------------------------------------------------------
  associate(isHe=>r%isHe, isH2=>r%isH2, ixHI=>r%ixHi, ixHII=>r%ixHII, ixHeII=>r%ixHeII, ixHeIII=>r%ixHeIII)

  ! Coll. Ionization Cooling
  ci_HI       = inp_coolrates_table(tables,tables%tbl_cr_ci_HI, T, .false., ci_HI_prime) * ne * nN(ixHII) ! ne * nHI
  ci_HI_prime = ci_HI_prime * ne * nN(ixHII)
  ! Collisional excitation cooling
  ce_HI          = inp_coolrates_table(tables,tables%tbl_cr_ce_HI, T, .false., ce_HI_prime) * ne * nN(ixHII) ! ne * nHI
  ce_HI_prime    = ce_HI_prime * ne * nN(ixHII)
  ! Recombination Cooling
  r_HII         = inp_coolrates_table(tables,tables%tbl_cr_r_HII, T, .false., r_HII_prime) * ne * nI(ixHII) ! ne * nHII
  r_HII_prime   = r_HII_prime * ne * nI(ixHII)

  if(isHe) then
     ! Coll. Ionization Cooling
     ci_HeI       = inp_coolrates_table(tables,tables%tbl_cr_ci_HeI, T, .false. , ci_HeI_prime)  * ne * nN(ixHeII)  ! ne * nHeI
     ci_HeI_prime = ci_HeI_prime  * ne * nN(ixHeII)
     ci_HeII       = inp_coolrates_table(tables,tables%tbl_cr_ci_HeII, T, .false. , ci_HeII_prime) * ne * nN(ixHeIII) ! ne * nHeII
     ci_HeII_prime = ci_HeII_prime * ne * nN(ixHeIII)
     ! Collisional excitation cooling
     ce_HeI       = inp_coolrates_table(tables,tables%tbl_cr_ce_HeI, T, .false.  , ce_HeI_prime) * ne * nN(ixHeII) ! ne * nHeI
     ce_HeI_prime = ce_HeI_prime *ne*nN(ixHeII)
     ce_HeII       = inp_coolrates_table(tables,tables%tbl_cr_ce_HeII, T, .false. , ce_HeII_prime) * ne * nN(ixHeIII) ! ne * nHeII
     ce_HeII_prime = ce_HeII_prime * ne * nN(ixHeIII)
     ! Recombination Cooling
     r_HeII       = inp_coolrates_table(tables,tables%tbl_cr_r_HeII, T, .false. , r_HeII_prime) * ne * nI(ixHeII) ! ne * nHeII
     r_HeII_prime = r_HeII_prime * ne * nI(ixHeII)
     r_HeIII       = inp_coolrates_table(tables,tables%tbl_cr_r_HeIII, T, .false. , r_HeIII_prime) * ne * nI(ixHeIII) ! ne * nHeIII
     r_HeIII_prime = r_HeIII_prime * ne * nI(ixHeIII)
  endif

  ! Bremsstrahlung
  bre    = inp_coolrates_table(tables,tables%tbl_cr_bre, T, .false., bre_prime)
  brefac = ne * nI(ixHII)
  if(isHe) then
     brefac = brefac + ne * ( nI(ixHeII) + 4. * nI(ixHeIII) )
  endif
  bre       = bre * brefac
  bre_prime = bre_prime * brefac

  ! Compton Cooling
  com       = inp_coolrates_table(tables,tables%tbl_cr_com, T, .true., com_prime) * ne
  com_prime = com_prime * ne

  ! Dielectronic recombination cooling
  if(isHe) then
     die       = inp_coolrates_table(tables,tables%tbl_cr_die, T, .false., die_prime) * ne * nN(ixHeIII) ! ne * nHeII
     die_prime = die_prime * ne * nN(ixHeIII)
  endif

  if(isH2) then
     ! Collisional dissociation (H2) cooling
     cr_H2HI       = inp_coolrates_table(tables,tables%tbl_cr_H2HI, T, .false. , cr_H2HI_prime) * nN(ixHI) * nN(ixHII) ! nH2 * nHI
     cr_H2HI_prime = cr_H2HI_prime * nN(ixHI)*nN(ixHII)

     cr_H2H2       = inp_coolrates_table(tables,tables%tbl_cr_H2H2, T, .false. , cr_H2H2_prime) * nN(ixHI)**2 ! nH2 * nH2
     cr_H2H2_prime = cr_H2H2_prime * nN(ixHI)**2
  endif

  ! Overall Cooling
  compCoolrate          = ci_HI    + r_HII    + ce_HI    + com + bre
  if(isHe) compCoolrate = compCoolrate                                   &
                        + ci_HeI   + r_HeII   + ce_HeI                   &
                        + ci_HeII  + r_HeIII  + ce_HeII  + die
  if(isH2) compCoolrate = compCoolrate                                   &
                        + cr_H2HI  + cr_H2H2

  dCooldT          = ci_HI_prime + r_HII_prime + ce_HI_prime + com_prime + bre_prime
  if(isHe) dCooldT = dCooldT                                             &
                   + ci_HeI_prime  + r_HeII_prime  + ce_HeI_prime        &
                   + ci_HeII_prime + r_HeIII_prime + ce_HeII_prime       &
                   + die_prime
  if(isH2) dCooldT = dCooldT                                             &
                   + cr_H2HI_prime + cr_H2H2_prime
  end associate

END FUNCTION compCoolrate

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION gamma_hi(T,J)
  implicit none
  real(kind=8),intent(in)::T,J
  !-------------------------------------------------------------------------
  real(kind=8)::gamma_hi, T3, jconsts
  !-------------------------------------------------------------------------
  T3 = T/1d3
  jconsts=0.33+0.9*exp(-1.0d0*((J-3.5d0)/0.9d0)**2)
  gamma_hi = &
       jconsts*(1.0d-11*sqrt(T3)/(1.0d0+60.0d0*T3**(-4)) + 1.0d-12*T3)

END FUNCTION gamma_hi
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
ELEMENTAL FUNCTION gamma_h2(T,J)
  implicit none
  real(kind=8),intent(in)::T,J
  !-------------------------------------------------------------------------
  real(kind=8)::gamma_h2, T3, jconsts
  !-------------------------------------------------------------------------
  T3 = T/1d3
  jconsts=(0.276*J**2)*exp(-1.0d0*(J/3.18d0)**1.7)
  gamma_h2 = jconsts*(3.3d-12 +6.6d-12*T3)

END FUNCTION gamma_h2
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

END MODULE coolrates_module
