! These are modules for reading, integrating and interpolating RT-relevant
! values from spectral tables, specifically SED (spectral energy
! distrbution) tables for stellar particle sources, as functions of age
! and metallicity, and UV-spectrum tables, as functions of redshift.
! The modules are:
! spectrum_integrator_module
!    For dealing with the integration of spectra in general
! SED_module
!    For dealing with the SED data
! UV_module
!    For dealing with the UV data
!_________________________________________________________________________

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! Module for integrating a wavelength-dependent spectrum
!_________________________________________________________________________
!
MODULE spectrum_integrator_module
!_________________________________________________________________________
  use amr_commons, only: run_t
  use constants, only: c_cgs, eV2erg, hplanck
  use cross_sections_module, only: getCrosssection_rtz
  implicit none

  PUBLIC integrateSpectrum, f1, fLambda, fdivLambda, fSig, fSigLambda,   &
         fSigdivLambda, trapz1

  PRIVATE   ! default

CONTAINS

!*************************************************************************
FUNCTION integrateSpectrum(run, X, Y, N, e0, e1, species, ion, func)

  ! Integrate spectral weighted function in energy interval [e0,e1]
  ! X      => Wavelengths [angstrom]
  ! Y      => Spectral luminosity per angstrom at wavelenghts [XX A-1]
  ! N      => Length of X and Y
  ! e0,e1  => Integrated interval [eV]
  ! species=> ion species, used as an argument in fx
  ! ionization=> ionization state of the species --> only needed for RTZ
  ! func   => Function which is integrated (of X, photon energy, Y, species)
  !-------------------------------------------------------------------------
  type(run_t) :: run
  real(kind=8) :: integrateSpectrum, X(N), Y(N), e0, e1
  integer :: N, species, ion
  interface
     real(kind=8) function func(run, wavelength, energy, intensity, species, ion)
       use amr_commons, only: run_t
       type(run_t) :: run
       real(kind=8) :: wavelength,energy,intensity
       integer :: species, ion
     end function func
  end interface
  !-------------------------------------------------------------------------
  real(kind=8),dimension(:),allocatable:: xx, yy, f
  real(kind=8):: la0, la1, emin, emax
  integer :: i
  !-------------------------------------------------------------------------
  integrateSpectrum=0.
  if(N .le. 2) RETURN
  ! Convert energy interval to wavelength interval
  la0 = X(1) ; la1 = X(N)
  if(e1.gt.0) la0 = max(la0, 1d8 * hplanck * c_cgs / e1 / eV2erg)
  if(e0.gt.0) la1 = min(la1, 1d8 * hplanck * c_cgs / e0 / eV2erg)
  if(la0 .ge. la1) then
     print*,'The energy limits do not overlap with SED range, so stopping'
     stop
  endif
  ! Photon energies of the integrand are held strictly inside (e0,e1)
  emin = 0d0 ; emax = huge(1d0)
  if(e0.gt.0) emin = nearest(e0, 1d0)
  if(e1.gt.0) emax = nearest(e1,-1d0)
  ! If we get here, the [la0, la1] inverval is completely within X
  allocate(xx(N)) ; allocate(yy(N)) ; allocate(f(N))
  xx = la0; yy = 0.; f = 0.
  i = 2
  do while ( i.lt.N .and. X(i).le.la0 )
     i = i+1                      !              Below wavelength interval
  enddo                           !   X(i) is now the first entry .gt. la0
  ! Interpolate to value at la0
  yy(i-1) = Y(i-1) + (xx(i-1)-X(i-1))*(Y(i)-Y(i-1))/(X(i)-X(i-1))
  f(i-1) = func(run, xx(i-1), photonEnergy(xx(i-1),emin,emax), yy(i-1), species, ion)
  do while ( i.lt.N .and. X(i).le.la1 )              ! Now within interval
     xx(i) = X(i) ; yy(i) = Y(i)
     f(i) = func(run, xx(i), photonEnergy(xx(i),emin,emax), yy(i), species, ion)
     i = i+1
  enddo                          ! i=N or X(i) is the first entry .gt. la1
  xx(i:) = la1                   !             Interpolate to value at la1
  yy(i) = Y(i-1) + (xx(i)-X(i-1))*(Y(i)-Y(i-1))/(X(i)-X(i-1))
  f(i) = func(run, xx(i), photonEnergy(xx(i),emin,emax), yy(i), species, ion)
  integrateSpectrum = trapz1(xx,f,i)
  deallocate(xx) ; deallocate(yy) ; deallocate(f)
END FUNCTION integrateSpectrum

!*************************************************************************
FUNCTION photonEnergy(lambda, emin, emax)

  ! Photon energy in eV at wavelength lambda [angstrom], within [emin,emax]
  !-------------------------------------------------------------------------
  real(kind=8) :: photonEnergy, lambda, emin, emax
  !-------------------------------------------------------------------------
  photonEnergy = min(max(hplanck * c_cgs/(lambda*1d-8) / eV2erg, emin), emax)
END FUNCTION photonEnergy

!*************************************************************************
! FUNCTIONS FOR USE WITH integrateSpectrum:
! lambda  => wavelengths in Angstrom
! energy  => photon energy at lambda in eV
! f       => function of wavelength (a spectrum in some units)
! species => 1=HI, 2=HeI or 3=HeII
!_________________________________________________________________________
FUNCTION f1(run, lambda, energy, f, species, ion)
  type(run_t) :: run
  real(kind=8) :: f1, lambda, energy, f
  integer :: species, ion
  f1 = f
END FUNCTION f1

FUNCTION fLambda(run, lambda, energy, f, species, ion)
  type(run_t) :: run
  real(kind=8):: fLambda, lambda, energy, f
  integer :: species, ion
  fLambda = f * lambda
END FUNCTION fLambda

FUNCTION fdivLambda(run, lambda, energy, f, species, ion)
  type(run_t) :: run
  real(kind=8):: fdivlambda, lambda, energy, f
  integer :: species, ion
  fdivLambda = f / lambda
END FUNCTION fdivLambda

FUNCTION fSig(run, lambda, energy, f, species, ion)
  type(run_t) :: run
  real(kind=8):: fSig, lambda, energy, f
  integer :: species, ion
#ifdef DO_RTZ
  fSig = f * getCrosssection_rtz(lambda, species, ion)
#else
  fSig = f * getCrosssection(run, energy, species)
#endif
END FUNCTION fSig

FUNCTION fSigLambda(run, lambda, energy, f, species, ion)
  type(run_t) :: run
  real(kind=8):: fSigLambda, lambda, energy, f
  integer :: species, ion
#ifdef DO_RTZ
  fSigLambda = f * lambda * getCrosssection_rtz(lambda, species, ion)
#else
  fSigLambda = f * lambda * getCrosssection(run, energy, species)
#endif
END FUNCTION fSigLambda

FUNCTION fSigdivLambda(run, lambda, energy, f, species, ion)
  type(run_t) :: run
  real(kind=8):: fSigdivLambda, lambda, energy, f
  integer :: species, ion
#ifdef DO_RTZ
  fSigdivLambda = f / lambda * getCrosssection_rtz(lambda, species, ion)
#else
  fSigdivLambda = f / lambda * getCrosssection(run, energy, species)
#endif
END FUNCTION fSigdivLambda
!_________________________________________________________________________

!*************************************************************************
FUNCTION trapz1(X, Y, N, cum)

  ! Integrates function Y(X) along the whole interval 1..N, using a very
  ! simple staircase method and returns the result.
  ! Optionally, the cumulative integral is returned in the cum argument.
  !-------------------------------------------------------------------------
  integer :: N,i
  real(kind=8) :: trapz1
  real(kind=8) :: X(N),Y(N)
  real(kind=8), optional :: cum(N)
  real(kind=8), allocatable :: cumInt(:)
  !-------------------------------------------------------------------------
  trapz1=0.
  if (N.le.1) RETURN
  allocate(cumInt(N))
  cumInt(:)=0d0
  do i=2,N
     cumInt(i)= cumInt(i-1) + abs(X(i)-X(i-1)) * (Y(i)+Y(i-1)) / 2d0
  end do
  trapz1 = cumInt(N)
  if(present(cum)) cum=cumInt
  deallocate(cumInt)
END FUNCTION trapz1

!*************************************************************************
#ifndef DO_RTZ
FUNCTION getCrosssection(run, E, species)

  ! Gives an atom-photon cross-section of given species at given photon energy,
  ! as given by Hui and Gnedin (1997) for HI and He, and Abel97 for H2.
  ! E       => Photon energy in eV
  ! species => ixHI (H2), ixHII (HI), ixHeII (HeI) or ixHeIII (HeII)
  ! returns :  photoionization or photodissociation cross-section in cm^2
  !------------------------------------------------------------------------
  type(run_t)  :: run
  real(kind=8) :: E, getCrosssection
  integer      :: species
  real(kind=8) :: E0=1., cs0=0., P=1., ya=1., yw=0., y0=0., y1=1.
  real(kind=8) :: x, y
  !------------------------------------------------------------------------
  if ( E .lt. run%ionEvs(species) ) then            ! below ionization energy
     getCrosssection=0.
     RETURN
  endif
  if(species .eq. run%ixHI) then ! H2 ionization cs, Abel1997 eqn A24
    getCrosssection = 0.
    if(E .ge. 11.2d0  .and. E .lt. 13.6d0) getCrosssection=2.1d-19 ! Sternberg2014
    if(E .ge. 15.42d0 .and. E .lt. 16.5d0) getCrosssection=6.2e-18*E-9.4e-17
    if(E .ge. 16.5d0  .and. E .lt. 17.7d0) getCrosssection=1.4e-18*E-1.48e-17
    if(E .ge. 17.7d0) getCrosssection=2.5e-14*E**(-2.71)
    RETURN
  endif
  if(species .eq. run%ixHII) then ! HI
     E0 = 4.298d-1 ; cs0 = 5.475d-14  ; P  = 2.963
     ya = 32.88    ; yw  = 0          ; y0 = 0         ; y1 = 0
  endif
  if(species .eq. run%ixHeII) then ! HeI
     E0 = 1.361d1  ; cs0 = 9.492d-16  ; P  = 3.188
     ya = 1.469    ; yw  = 2.039      ; y0 = 0.4434    ; y1 = 2.136
  endif
  if(species .eq. run%ixHeIII) then ! HeII
     E0 = 1.720    ; cs0 = 1.369d-14  ; P  = 2.963
     ya = 32.88    ; yw  = 0          ; y0 = 0         ; y1 = 0
  endif
  x = E/E0 - y0
  y = sqrt(x**2+y1**2)

  getCrosssection = cs0 * ((x-1.)**2 + yw**2) * y**(0.5*P-5.5)/(1.+sqrt(y/ya))**P

END FUNCTION getCrosssection
#endif

END MODULE spectrum_integrator_module

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! Module for Stellar Energy Distribution table.
!_________________________________________________________________________
!
MODULE SED_module
  !_________________________________________________________________________
  use amr_commons, only: run_t, global_t
  use constants, only: L_sun, c_cgs, eV2erg, hplanck, Gyr2sec
  use hydro_parameters, only: nion
  use rt_parameters, only: nrtgrp
  implicit none

  PUBLIC sed_table_t, init_SED_table, inp_SED_table &
        ,update_SED_group_props, getNPhotonsEmitted
#ifdef DO_RTZ
  PUBLIC initialize_cross_sections_from_blackbody &
        ,initialize_group_energies_from_blackbody
#endif

  PRIVATE   ! default

  type sed_table_t

     !-------------------------------------------------------------------------
     logical :: is_SED_single_Z
     !-------------------------------------------------------------------------
     ! Light properties for different spectral energy distributions------------
     integer :: nA, nZ=8              ! Number of age bins and Z bins
     !-------------------------------------------------------------------------
     ! Age and z logarithmic intervals and lowest values:
     real(kind=8) :: dlgA=0.02d0
     real(kind=8) :: dlgZ
     real(kind=8) :: lgA0, lgZ0
     real(kind=8), allocatable, dimension(:) :: ages, zeds ! [Gyr], [m_met/m_gas]
     !-------------------------------------------------------------------------
     ! SED_table: iAges, imetallicities, igroups, properties
     !                                         (Lum, Lum-acc, egy, csn, cse)
     ! Lum is photons per sec per solar mass (eV per sec per solar mass).
     ! Lum-acc is accumulated lum.
     real(kind=8), allocatable, dimension(:,:,:,:) :: table
     !-------------------------------------------------------------------------
     ! ia, iz: lower indexes: 0<ia<sed_nA etc.
     ! da0, da1, dz0, dz1: proportional distances from edges:
     ! 0<=da0<=1, 0<=da1<=1 etc.
     integer :: ia, iz
     real(kind=8) :: da, da0, da1, dz, dz0, dz1
     !-------------------------------------------------------------------------

  end type sed_table_t

CONTAINS

!*************************************************************************
SUBROUTINE init_SED_table(r, g, SED)

  ! Initiate SED properties table, which gives photon luminosities,
  ! integrated luminosities, average photon cross sections and energies of
  ! each photon group as a function of stellar population age and
  ! metallicity.  The SED is read from a directory specified by sed_dir.
  !-------------------------------------------------------------------------
  use spectrum_integrator_module
#ifdef DO_RTZ
  use rtz_module
#endif
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(run_t) :: r
  type(global_t) :: g
  type(sed_table_t) :: SED
#ifndef WITHOUTMPI
  real(kind=8),allocatable::tbl2(:,:,:)
  integer::ierr
#endif
  ! Temporary SSP/SED parameters (read from SED files):
  integer:: nAges, nzs, nLs              ! # of bins of age, z, wavelength
  real(kind=8),allocatable::ages(:), Zs(:), Ls(:), rebAges(:)
  real(kind=8),allocatable::SEDs(:,:,:)           ! SEDs f(lambda,age,met)
  real(kind=8),allocatable::tbl(:,:,:), reb_tbl(:,:,:)
  integer::i,ia,iz,ip,ii,jj,dum,counter
  character(len=128)::fZs, fAges, fSEDs                        ! Filenames
  logical::ok,okAge,okZ
  real(kind=8)::dlgA, pL0, pL1, tmp
  integer::locid,ncpu2
  integer::nv
  integer,parameter::tag=1132

! set nv to the correct value depending on whether we are using rtz or not
#ifdef DO_RTZ
  nv=3+2
  do i=1,n_elements
     if (elements(i)%atomic_number.gt.0) then
        nv = nv + (2 * elements(i)%n_ions)
     end if
     if (r%isH2_rtz) then
        nv = nv + 2
     end if
  end do
#else
  nv=3+2*nion ! # vars in SED table: L,Lacc,egy,nions*(csn,egy)
#endif
  !-------------------------------------------------------------------------
  if(g%myid==1)write(*,*) 'Stars are photon emitting, so initializing SED table'
  if(r%sed_dir.eq.'')call get_environment_variable('RAMSES_SED_DIR',r%sed_dir)
  inquire(FILE=TRIM(r%sed_dir)//'/all_seds.dat', exist=ok)
  if(.not. ok)then
     if(g%myid.eq.1) then
        write(*,*)'Cannot access SED directory ',TRIM(r%sed_dir)
        write(*,*)'Directory '//TRIM(r%sed_dir)//' not found'
        write(*,*)'You need to set the RAMSES_SED_DIR envvar' // &
                  ' to the correct path, or use the namelist.'
     endif
     stop
  end if
  write(fZs,'(a,a)')   trim(r%sed_dir),"/metallicity_bins.dat"
  write(fAges,'(a,a)') trim(r%sed_dir),"/age_bins.dat"
  write(fSEDs,'(a,a)') trim(r%sed_dir),"/all_seds.dat"
  inquire(file=fZs, exist=okZ)
  inquire(file=fAges, exist=okAge)
  inquire(file=fSEDs, exist=ok)
  if(.not. ok .or. .not. okAge .or. .not. okZ) then
     if(g%myid.eq.1) then
        write(*,*) 'Cannot read SED files...'
        write(*,*) 'Check if SED-directory contains the files ',  &
                   'metallicity_bins.dat, age_bins.dat, all_seds.dat'
     endif
     stop
  end if

  ! READ METALLICITY BINS-------------------------------------------------
  open(unit=10,file=fZs,status='old',form='formatted')
  read(10,'(i8)') nzs
  allocate(zs(nzs))
  do i = 1, nzs
     read(10,'(e14.6)') zs(i)
  end do
  close(10)
  if(nzs.eq.1)SED%is_SED_single_Z=.true.

  ! READ AGE BINS---------------------------------------------------------
  open(unit=10,file=fAges,status='old',form='formatted')
  read(10,'(i8)') nAges
  allocate(ages(nAges))
  do i = 1, nAges
     read(10,'(e14.6)') ages(i)
  end do
  close(10)
  if(nAges.lt.2)then
      if(g%myid==1) print*,'WARNING! Only one age bin found - check if interpolated values make sense'
  endif
  ages = ages*1.e-9 ! Convert from yr to Gyr
  if(ages(1) .ne. 0.) ages(1) = 0.

  ! READ SEDS-------------------------------------------------------------
  open(unit=10,file=fSEDs,status='old',form='unformatted')
  read(10) nLs, dum
  allocate(Ls(nLs))
  read(10) Ls(:)
  allocate(SEDs(nLs,nAges,nzs))
  do iz = 1, nzs
     do ia = 1, nAges
        read(10) SEDs(:,ia,iz)
     end do
  end do
  close(10)

  ! Do not interpolate and update SEDs if single metallicity and age
  if(nZs.eq.1 .and. nAges<3)then
     SED%nZ=1
  end if

  ! If MPI then share the SED integration between the cpus:
#ifndef WITHOUTMPI
  call MPI_COMM_RANK(MPI_COMM_WORLD,locid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu2,ierr)
#endif
#ifdef WITHOUTMPI
  locid=0
  ncpu2=1
#endif

  ! Perform SED integration of luminosity, csn and egy per (age,Z) bin----
  allocate(tbl(nAges,nZs,nv))
  do ip = 1, nrtgrp                                ! Loop photon groups
     tbl = 0.
     pL0 = r%group_L0(ip) ; pL1 = r%group_L1(ip)! eV interval of photon group ip
     do iz = 1, nzs                                     ! Loop metallicity
     do ia = locid+1,nAges,ncpu2                                ! Loop age
        tbl(ia,iz,1) = getSEDLuminosity(r,Ls,SEDs(:,ia,iz),nLs,pL0,pL1)
        tbl(ia,iz,3) = getSEDEgy(r,Ls,SEDs(:,ia,iz),nLs,pL0,pL1)
#ifdef DO_RTZ
        counter = 1
        do ii=1, n_elements ! Loop over elements
           ! Cross sections for atomic species
           if (elements(ii)%atomic_number.gt.0) then
              do jj=1,elements(ii)%n_ions-1 !loop over ionization states
                 tbl(ia,iz,3+counter) = getSEDcsn(r,Ls,SEDs(:,ia,iz),nLs,pL0,pL1,ii,jj)
                 counter = counter + 1
                 tbl(ia,iz,3+counter) = getSEDcse(r,Ls,SEDs(:,ia,iz),nLs,pL0,pL1,ii,jj)
                 counter = counter + 1
              end do
           end if
        end do

        ! Deal with molecules separately
        if (elements(1)%atomic_number.gt.0 .and. r%isH2_rtz) then
           tbl(ia,iz,2+counter) = getSEDcsn(r,Ls,SEDs(:,ia,iz),nLs,pL0,pL1,1,3)
           counter = counter + 1
           tbl(ia,iz,2+counter) = getSEDcse(r,Ls,SEDs(:,ia,iz),nLs,pL0,pL1,1,3)
           counter = counter + 1
        end if
#else
        do ii = 1, nion                                     ! Loop species --> put dummy ion
           tbl(ia,iz,2+ii*2) = getSEDcsn(r,Ls,SEDs(:,ia,iz),nLs,pL0,pL1,ii,1)
           tbl(ia,iz,3+ii*2) = getSEDcse(r,Ls,SEDs(:,ia,iz),nLs,pL0,pL1,ii,1)
        end do ! End species loop
#endif
     end do ! End age loop
     end do ! End Z loop

#ifndef WITHOUTMPI
     allocate(tbl2(nAges,nzs,nv))
     call MPI_ALLREDUCE(tbl,tbl2,nAges*nzs*nv,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
     tbl = tbl2
     deallocate(tbl2)
#endif

     ! Now the SED properties are in tbl...just need to rebin it, to get !
     ! even log-intervals between bins for fast interpolation.           !
     dlgA = SED%dlgA ; SED%dlgZ = -SED%nz
     call rebin_log(dlgA, SED%dlgZ                                       &
          , tbl(2:nAges,:,:), nAges-1, nZs, ages(2:nAges), zs, nv        &
          , reb_tbl, SED%nA, SED%nZ, rebAges, SED%Zeds)
     SED%nA = SED%nA + 1                              ! Make room for zero age
     if(ip .eq. 1) allocate(SED%table(SED%nA, SED%nZ, nrtgrp, nv))
     SED%table(1, :,ip,:) = reb_tbl(1,:,:)            ! Zero age properties
     SED%table(1, :,ip,2) = 0.                        !  Lacc=0 at zero age
     SED%table(2:,:,ip,:) = reb_tbl
     deallocate(reb_tbl)

     if(ip .eq. 1) then
        SED%lgZ0 = log10(SED%Zeds(1))                  ! Interpolation intervals
        SED%lgA0 = log10(rebAges(1))
        allocate(SED%ages(SED%nA))
        SED%ages(1)=0d0 ; SED%ages(2:)=rebAges ;    ! Must have zero initial age
     end if

     ! Integrate the cumulative luminosities:
     SED%table(:,:,ip,2)=0d0
     do iz = 1, SED%nZ ! Loop metallicity
        tmp = trapz1( SED%ages, SED%table(:,iz,ip,1), SED%nA, SED%table(:,iz,ip,2) )
        SED%table(:,iz,ip,2) = SED%table(:,iz,ip,2) * Gyr2sec
     end do

  end do ! End photon group loop

  deallocate(SEDs) ; deallocate(tbl)
  deallocate(ages) ; deallocate(rebAges)
  deallocate(zs)
  deallocate(Ls)

  if (g%myid==1) call write_SED_table(SED)

END SUBROUTINE init_SED_table

!*************************************************************************
SUBROUTINE update_SED_group_props(r, g, SED, p)

  ! Compute and assign to all SED photon groups an average of the
  ! quantities group_csn and group_egy from all the star particles in the
  ! simulation, weighted by the luminosity of the particles.
  ! If there are no stars, we assign the SED properties of zero-age,
  ! zero-metallicity stellar populations.
  !-------------------------------------------------------------------------
#ifndef WITHOUTMPI
  use mpi
#endif
  use pm_commons, only: part_t
  use amr_commons, only: global_t
  use hydro_parameters, only: nion
  use constants, only: Gyr2sec
#ifdef DO_RTZ
  use rtz_module, only: elements, n_elements
#endif
  type(run_t) :: r
  type(part_t) :: p
  type(global_t) :: g
  type(sed_table_t) :: SED
  !-------------------------------------------------------------------------
  integer :: i, ip, ii, jj, counter
#ifndef WITHOUTMPI
  integer::info
#endif
  real(kind=8), dimension(1:nrtgrp) :: L_star
  real(kind=8), dimension(1:nrtgrp) :: egy_star
  real(kind=8), dimension(1:nrtgrp) :: sum_L_cpu, sum_L_all
  real(kind=8), dimension(1:nrtgrp) :: sum_egy_cpu, sum_egy_all
#ifdef DO_RTZ
  real(kind=8), dimension(1:nrtgrp,1:27,1:27) :: csn_star, cse_star
  real(kind=8), dimension(1:nrtgrp,1:27,1:27) :: sum_csn_cpu, sum_csn_all
  real(kind=8), dimension(1:nrtgrp,1:27,1:27) :: sum_cse_cpu, sum_cse_all
#else
  real(kind=8), dimension(1:nrtgrp,1:nion) :: csn_star, cse_star
  real(kind=8), dimension(1:nrtgrp,1:nion) :: sum_csn_cpu, sum_csn_all
  real(kind=8), dimension(1:nrtgrp,1:nion) :: sum_cse_cpu, sum_cse_all
#endif
  real(kind=8) :: mass, age, Z, t_SN_gyr, scale_t_Gyr
  real(kind=8) :: scale_nH, scale_T2, scale_l, scale_d, scale_t, scale_v
  !-------------------------------------------------------------------------

  ! Conversion factor from user units to cgs units
  call units(r,g,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  scale_t_Gyr = (scale_t/g%aexp**2)/Gyr2sec

  ! Supernovae progenitors life time from Myr to proper time in code units
  !t_SN = r%t_SNII*1d6*(365.*24.*3600.)/(scale_t/g%aexp**2)
  t_SN_gyr = r%t_SNII/1d3

  sum_L_cpu   = 0d0 ! Accumulated luminosity, avg cross sections and
  sum_egy_cpu = 0d0 ! photon energies for all stars belonging to
  sum_csn_cpu = 0d0 ! 'this' cpu
  sum_cse_cpu = 0d0

  ! Loop over all star particles
  do i = 1, p%npart
     mass = p%mp(i)
     age = (g%texp-p%tp(i)) * scale_t_Gyr ! stellar age in Gyrs
     ! Account for stellar mass loss - SED uses initial population mass
     if(age.gt.t_SN_gyr)then
        mass = mass / (1d0-r%eta_SNII)
     endif
     if(r%metal) then
        Z = max(dble(p%zp(i)), 1d-5)                               ! [m_metals/m_tot]
     else
        Z = max(dble(r%z_ave)*0.02d0, 1d-5)                          ! [m_metals/m_tot]
     endif
     call inp_SED_table(SED, age, Z, 1, .false., L_star)     !  [# s-1 M_sun-1]
     call inp_SED_table(SED, age, Z, 3, .true., egy_star)    !             [eV]

#ifdef DO_RTZ
     counter = 1
     do ii=1, n_elements ! Loop over elements
     ! Cross sections for atomic species
        if (elements(ii)%atomic_number.gt.0) then
           do jj=1,elements(ii)%n_ions-1 !loop over ionization states
              call inp_SED_table(SED, age, Z, 3+counter, .true., csn_star(1:nrtgrp,ii,jj))! [cm^2]
              counter = counter + 1
              call inp_SED_table(SED, age, Z, 3+counter, .true., cse_star(1:nrtgrp,ii,jj))! [cm^2]
              counter = counter + 1
           end do
        end if
     end do

     ! Deal with molecules separately
     if (elements(1)%atomic_number.gt.0 .and. r%isH2_rtz) then
        call inp_SED_table(SED, age, Z, 3+counter, .true., csn_star(1:nrtgrp,1,3))! [cm^2]
        counter = counter + 1
        call inp_SED_table(SED, age, Z, 3+counter, .true., cse_star(1:nrtgrp,1,3))! [cm^2]
        counter = counter + 1
     end if
#else
     do ii = 1, nion
        call inp_SED_table(SED, age, Z, 2+2*ii, .true., csn_star(1:nrtgrp,ii))! [cm^2]
        call inp_SED_table(SED, age, Z, 3+2*ii, .true., cse_star(1:nrtgrp,ii))! [cm^2]
     end do
#endif

     do ip = 1, nrtgrp
        L_star(ip) = L_star(ip) * mass                  !       [# photons s-1]
        sum_L_cpu(ip) = sum_L_cpu(ip) + L_star(ip)
        sum_egy_cpu(ip) = sum_egy_cpu(ip) + L_star(ip) * egy_star(ip)
#ifdef DO_RTZ
        sum_csn_cpu(ip,1:27,1:27) = sum_csn_cpu(ip,1:27,1:27) + L_star(ip) * csn_star(ip,1:27,1:27)
        sum_cse_cpu(ip,1:27,1:27) = sum_cse_cpu(ip,1:27,1:27) + L_star(ip) * cse_star(ip,1:27,1:27)
#else
        sum_csn_cpu(ip,1:nion) = sum_csn_cpu(ip,1:nion) + L_star(ip) * csn_star(ip,1:nion)
        sum_cse_cpu(ip,1:nion) = sum_cse_cpu(ip,1:nion) + L_star(ip) * cse_star(ip,1:nion)
#endif
     end do

  end do

  ! Sum up for all cpus
#ifdef WITHOUTMPI
  sum_L_all = sum_L_cpu
  sum_egy_all = sum_egy_cpu
  sum_csn_all = sum_csn_cpu
  sum_cse_all = sum_cse_cpu
#else
  call MPI_ALLREDUCE(sum_L_cpu,   sum_L_all,   nrtgrp,               &
                     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(sum_egy_cpu, sum_egy_all, nrtgrp,               &
                     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
#ifdef DO_RTZ
  call MPI_ALLREDUCE(sum_csn_cpu, sum_csn_all, nrtgrp*27*27,         &
                     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(sum_cse_cpu, sum_cse_all, nrtgrp*27*27,         &
                     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
#else
  call MPI_ALLREDUCE(sum_csn_cpu, sum_csn_all, nrtgrp*nion,         &
                     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
  call MPI_ALLREDUCE(sum_cse_cpu, sum_cse_all, nrtgrp*nion,         &
                     MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
#endif
#endif

  ! ...and take averages weighted by luminosities
  do ip = 1, nrtgrp
     ! No update for non-SED groups (L0>L1):
     if(r%group_L0(ip).ne. 0d0 .and. r%group_L1(ip) .ne. 0d0 .and. &
          &  (r%group_L0(ip) .ge. r%group_L1(ip)) ) cycle
     ! We have star particles already
     if(sum_L_all(ip) .gt. 0.) then
        r%group_egy(ip) = sum_egy_all(ip) / sum_L_all(ip)
#ifdef DO_RTZ
        r%group_csn(ip,1:27,1:27) = sum_csn_all(ip,1:27,1:27) / sum_L_all(ip)
        r%group_cse(ip,1:27,1:27) = sum_cse_all(ip,1:27,1:27) / sum_L_all(ip)
#else
        r%group_csn(ip,1:nion) = sum_csn_all(ip,1:nion) / sum_L_all(ip)
        r%group_cse(ip,1:nion) = sum_cse_all(ip,1:nion) / sum_L_all(ip)
#endif
     else ! no stars -> assign zero-age zero-metallicity props
        r%group_egy(ip) = SED%table(1,1,ip,3)
#ifdef DO_RTZ
     counter = 1
     do ii=1, n_elements ! Loop over elements
     ! Cross sections for atomic species
        if (elements(ii)%atomic_number.gt.0) then
           do jj=1,elements(ii)%n_ions-1 !loop over ionization states
              r%group_csn(ip,ii,jj) = SED%table(1,1,ip,3+counter)
              counter = counter + 1
              r%group_cse(ip,ii,jj) = SED%table(1,1,ip,3+counter)
              counter = counter + 1
           end do
        end if
     end do

     ! Deal with molecules separately
     if (elements(1)%atomic_number.gt.0 .and. r%isH2_rtz) then
        r%group_csn(ip,1,3) = SED%table(1,1,ip,3+counter)
        counter = counter + 1
        r%group_cse(ip,1,3) = SED%table(1,1,ip,3+counter)
        counter = counter + 1
     end if
#else
        do ii = 1, nion
           r%group_csn(ip,ii) = SED%table(1,1,ip,2+2*ii)
           r%group_cse(ip,ii) = SED%table(1,1,ip,3+2*ii)
        enddo
#endif
     endif
  end do

  if(g%myid==1) write(*,*)'SED Photon groups updated through stellar polling'

END SUBROUTINE update_SED_group_props

!*************************************************************************
FUNCTION getSEDLuminosity(run, X, Y, N, e0, e1)

  ! Compute and return luminosity in energy interval (e0,e1) [eV]
  ! in SED Y(X). Assumes X is in Angstroms and Y in Lo/Angstroms/Msun.
  ! (Lo=[Lo_sun], Lo_sun=[erg s-1]. total solar luminosity is
  ! Lo_sun=10^33.58 erg/s)
  ! returns: Photon luminosity in, [# s-1 Msun-1]
  !-------------------------------------------------------------------------
  use spectrum_integrator_module
  type(run_t) :: run
  integer :: N
  real(kind=8) :: getSEDLuminosity, X(N), Y(N), e0, e1
  !-------------------------------------------------------------------------
  integer :: species = 1                 ! irrelevant but must be included
  integer :: ion = 1                 ! irrelevant but must be included
  real(kind=8), parameter :: const=1.0e-8/hplanck/c_cgs
  ! const is a div by ph energy => ph count.  1e-8 is a conversion into
  ! cgs, since wly=[angstrom] h=[erg s-1], c=[cm s-1]
  !-------------------------------------------------------------------------
  getSEDLuminosity = const * integrateSpectrum(run, X, Y, N, e0, e1, species, ion, fLambda)
  getSEDLuminosity = getSEDLuminosity * L_sun  ! Scale by solar luminosity

END FUNCTION getSEDLuminosity

!*************************************************************************
FUNCTION getSEDEgy(run, X, Y, N, e0, e1)

  ! Compute average energy, in eV, in energy interval (e0,e1) [eV] in SED
  ! Y(X). Assumes X is in Angstroms and Y is energy weight per angstrom
  ! (not photon count).
  !-------------------------------------------------------------------------
  use spectrum_integrator_module
  type(run_t) :: run
  integer :: N
  real(kind=8) :: getSEDEgy, X(N), Y(N), e0, e1
  !-------------------------------------------------------------------------
  integer :: species = 1
  integer :: ion = 1
  real(kind=8) :: norm
  real(kind=8), parameter :: const=1d8*hplanck*c_cgs/eV2erg ! energy conversion
  !-------------------------------------------------------------------------
  norm = integrateSpectrum(run, X, Y, N, e0, e1, species, ion, fLambda)
  getSEDEgy = const * integrateSpectrum(run, X, Y, N, e0, e1, species, ion, f1) / norm

END FUNCTION getSEDEgy

!*************************************************************************
FUNCTION getSEDcsn(run, X, Y, N, e0, e1, species, ion)

  ! Compute and return average photoionization
  ! cross-section, in cm^2, for a given energy interval (e0,e1) [eV] in
  ! SED Y(X). Assumes X is in Angstroms and that Y is energy weight per
  ! angstrom (not photon #).
  ! Species is a code for the ion in question: 1=HI, 2=HeI, 3=HeIII
  !-------------------------------------------------------------------------
  use spectrum_integrator_module
  type(run_t) :: run
  integer :: N, species, ion
  real(kind=8) :: getSEDcsn, X(N), Y(N), e0, e1
  !-------------------------------------------------------------------------
  real(kind=8) :: norm
  !-------------------------------------------------------------------------
#ifdef DO_RTZ
  if(e1 .gt. 0. .and. e1 .le. run%ionEvs(species, ion)) then
#else
  if(e1 .gt. 0. .and. e1 .le. run%ionEvs(species)) then
#endif
     getSEDcsn=0. ; RETURN    ! [e0,e1] below ionization energy of species
  endif
  norm = integrateSpectrum(run, X, Y, N, e0, e1, species, ion, fLambda)
  getSEDcsn = integrateSpectrum(run, X, Y, N, e0, e1, species, ion, fSigLambda) / norm

END FUNCTION getSEDcsn

!************************************************************************
FUNCTION getSEDcse(run, X, Y, N, e0, e1, species, ion)

  ! Compute and return average energy weighted photoionization
  ! cross-section, in cm^2, for a given energy interval (e0,e1) [eV] in
  ! SED Y(X). Assumes X is in Angstroms and that Y is energy weight per
  ! angstrom (not photon #).
  ! Species is a code for the ion in question: 1=HI, 2=HeI, 3=HeIII
  !-------------------------------------------------------------------------
  use spectrum_integrator_module
  type(run_t) :: run
  integer :: N, species, ion
  real(kind=8) :: getSEDcse, X(N), Y(N), e0, e1
  !-------------------------------------------------------------------------
  real(kind=8) :: norm
  !-------------------------------------------------------------------------
#ifdef DO_RTZ
  if(e1 .gt. 0. .and. e1 .le. run%ionEvs(species, ion)) then
#else
  if(e1 .gt. 0. .and. e1 .le. run%ionEvs(species)) then
#endif
     getSEDcse = 0. ; RETURN    ! [e0,e1] below ionization energy of species
  endif

  norm = integrateSpectrum(run, X, Y, N, e0, e1, species, ion, f1)
  getSEDcse = integrateSpectrum(run, X, Y, N, e0, e1, species, ion, fSig) / norm

END FUNCTION getSEDcse

!*************************************************************************
SUBROUTINE rebin_log(xint_log, yint_log,                         &
     &               data, nx, ny, x, y, nz,                     &
     &               new_data, new_nx, new_ny, new_x, new_y      )

  ! Rebin the given 2d data into constant logarithmic intervals, using
  ! linear interpolation.
  ! xint_log,  => x and y intervals in the rebinned data. If negative,
  ! yint_log      these values represent the new number of bins.
  ! data       => The 2d data to be rebinned
  ! nx,ny      => Number of points in x and y in the original data
  ! nz         => Number of values in the data
  ! x,y        => x and y values for the data
  ! new_data  <=  The rebinned 2d data
  ! new_nx,ny <=  Number of points in x and y in the rebinned data
  ! new_x,y   <=  x and y point values for the rebinned data
  !-------------------------------------------------------------------------
  integer, intent(in) :: nx, ny, nz
  integer :: new_nx, new_ny
  real(kind=8) :: xint_log, yint_log
  real(kind=8), intent(in) :: x(nx),y(ny)
  real(kind=8), intent(in) :: data(nx,ny,nz)
  real(kind=8), dimension(:), allocatable :: new_x, new_y
  real(kind=8), dimension(:,:,:), allocatable :: new_data
  !-------------------------------------------------------------------------
  real(kind=8), dimension(:), allocatable :: new_lgx, new_lgy
  real(kind=8) :: dx0, dx1, dy0, dy1, x_step, y_step
  real(kind=8) :: x0lg, x1lg, y0lg, y1lg
  integer :: i, j, ix, iy, ix1, iy1
  !-------------------------------------------------------------------------
  if(allocated(new_x)) deallocate(new_x)
  if(allocated(new_y)) deallocate(new_y)
  if(allocated(new_data)) deallocate(new_data)

  ! Find dimensions of the new_data and initialize it
  x0lg = log10(x(1));   x1lg = log10(x(nx))
  y0lg = log10(y(1));   y1lg = log10(y(ny))

  if(xint_log .lt. 0 .and. nx .gt. 1) then
     new_nx=int(-xint_log)                        ! xint represents wanted
     xint_log = (x1lg-x0lg)/(new_nx-1)            !     number of new bins
  else
     new_nx =int((x1lg-x0lg)/xint_log) + 1
  endif
  allocate(new_x(new_nx)) ; allocate(new_lgx(new_nx))
  do i = 0, new_nx-1                              !  initialize the x-axis
     new_lgx(i+1) = x0lg + i*xint_log
  end do
  new_x=10d0**new_lgx

  if(yint_log .lt. 0 .and. ny .gt. 1) then        ! yint represents wanted
     new_ny=int(-yint_log)                        !     number of new bins
     yint_log = (y1lg-y0lg)/(new_ny-1)
  else
     new_ny = int((y1lg-y0lg)/yint_log) + 1
  endif
  allocate(new_y(new_ny)) ; allocate(new_lgy(new_ny))
  do j = 0, new_ny-1                              !      ...and the y-axis
     new_lgy(j+1) = y0lg + j*yint_log
  end do
  new_y=10d0**new_lgy

  ! Initialize new_data and find values for each point in it
  allocate(new_data(new_nx, new_ny, nz))
  do j = 1, new_ny
     call locate(y, ny, new_y(j), iy)
     ! y(iy) <= new_y(j) <= y(iy+1)
     ! iy is lower bound, so it can be zero but not larger than ny
     if(iy < 1) iy=1
     if (iy < ny) then
        iy1  = iy + 1
        y_step = y(iy1) - y(iy)
        dy0  = max(new_y(j) - y(iy),    0.0d0)  / y_step
        dy1  = min(y(iy1)   - new_y(j), y_step) / y_step
     else
        iy1  = iy
        dy0  = 0.0d0 ;  dy1  = 1.0d0
     end if

     do i = 1, new_nx
        call locate(x, nx, new_x(i), ix)
        if(ix < 1) ix=1
        if (ix < nx) then
           ix1  = ix+1
           x_step = x(ix1)-x(ix)
           dx0  = max(new_x(i) - x(ix),    0.0d0)  / x_step
           dx1  = min(x(ix1)   - new_x(i), x_step) / x_step
        else
           ix1  = ix
           dx0  = 0.0d0 ;  dx1  = 1.0d0
        end if

        if (abs(dx0+dx1-1.0d0) .gt. 1.0d-6 .or.                          &
                                         abs(dy0+dy1-1.0d0) > 1.0d-6) then
           write(*,*) 'Screwed up the rebin interpolation ... '
           write(*,*) dx0+dx1,dy0+dy1
           stop
        end if

        new_data(i,j,:) =                                                &
             dx0 * dy0 * data(ix1,iy1,:) + dx1 * dy0 * data(ix, iy1,:) + &
             dx0 * dy1 * data(ix1,iy, :) + dx1 * dy1 * data(ix, iy, :)
     end do
  end do

  deallocate(new_lgx)
  deallocate(new_lgy)

END SUBROUTINE rebin_log

!*************************************************************************
SUBROUTINE write_SED_table(SED)
  ! Write the SED properties to a file (this is just in debugging, to check
  ! if the SEDs are being read correctly).
  ! Photon group properties: age [Gyr], metal mass fraction,
  ! luminosity [photons s-1 Msun-1], cumulative luminosity [photons Msun-1],
  ! group energy [ergs], csn_HX [cm2], cse_HX [cm2], where HX=H2, HI, HeI,
  ! and HeII; H2 and He are optional
  !-------------------------------------------------------------------------
  use hydro_parameters, only: nion
#ifdef DO_RTZ
  use rtz_module, only: n_elements, elements
#endif
  type(sed_table_t) :: SED
  !-------------------------------------------------------------------------
  character(len=128) :: filename
  integer :: ip, i, j, k, nv
  !-------------------------------------------------------------------------

#ifdef DO_RTZ
  do i=1,n_elements
     if (elements(i)%atomic_number.gt.0) then
        nv = nv + (2 * elements(i)%n_ions)
     end if
#if N_H2 > 0
     nv = nv + 2
#endif
  end do
#endif

  do ip = 1, nrtgrp
     write(filename,'(A, I1, A)') 'SEDtable', ip, '.list'
     open(10, file=filename, status='unknown')
     write(10,*) SED%nA, SED%nZ

     do j = 1,SED%nz
        do i = 1,SED%nA
           write(10,900,advance='no')                                    &
                SED%ages(i)        ,    SED%zeds(j)        ,            &
                SED%table(i,j,ip,1),    SED%table(i,j,ip,2),            &
                SED%table(i,j,ip,3)
#ifdef DO_RTZ
           do k=1,(nv/2)-1
              write(10,901,advance='no') SED%table(i,j,ip,2+2*k), SED%table(i,j,ip,3+2*k)
           end do
           write(10,901) SED%table(i,j,ip,2+2*(nv/2)), SED%table(i,j,ip,3+2*(nv/2))
#else
           if(nion .gt. 1) then
              do k = 1,nion-1
                 write(10,901,advance='no')                              &
                      SED%table(i,j,ip,2+2*k), SED%table(i,j,ip,3+2*k)
              enddo
           endif
           write(10,901)                                                 &
                SED%table(i,j,ip,2+2*nion), SED%table(i,j,ip,3+2*nion)
#endif
        end do
     end do
     close(10)
  end do
900 format (ES15.4, ES15.4, ES15.4, ES15.4, f15.4)
901 format (ES15.4, ES15.4)

END SUBROUTINE write_SED_table

!*************************************************************************
SUBROUTINE inp_SED_table(SED, age, Z, nProp, same, ret)

  ! Compute SED property by interpolation from table.
  ! input/output:
  ! age   => Star population age [Gyrs]
  ! Z     => Star population metallicity [m_metals/m_tot]
  ! nprop => Number of property to fetch
  !          1=log(photon # intensity [# Msun-1 s-1]),
  !          2=log(cumulative photon # intensity [# Msun-1]),
  !          3=avg_egy, 2+2*iIon=avg_csn, 3+2*iIon=avg_cse
  ! same  => If true then assume same age and Z as used in last call.
  !          In this case the interpolation indexes can be recycled.
  ! ret   => The interpolated values of the sed property for every photon
  !          group
  !-------------------------------------------------------------------------
  type(sed_table_t) :: SED
  real(kind=8), intent(in) :: age, Z
  integer :: nProp
  logical :: same
  real(kind=8), dimension(1:nrtgrp) :: ret
  !-------------------------------------------------------------------------
  real(kind=8) :: lgAge, lgZ, da, dz
  !-------------------------------------------------------------------------
  if(.not. same) then
     if(age.le.0d0) then
        lgAge = -4d0
     else
        lgAge = log10(age)
     endif
     SED%ia = min(max(floor((lgAge-SED%lgA0)/SED%dlgA ) + 2, 1  ),  SED%nA-1 )
     da = SED%ages(SED%ia+1)-SED%ages(SED%ia)
     SED%da0 = min( max(   (age-SED%ages(SED%ia)) /da,       0d0 ), 1d0          )
     SED%da1 = min( max(  (SED%ages(SED%ia+1)-age)/da,       0d0 ), 1d0          )

     if(SED%is_SED_single_Z)then
        SED%iz = 1
        SED%dz0 = 0.0d0
        SED%dz1 = 1.0d0
     else
        lgZ = log10(Z)
        SED%iz = min(max(floor((lgZ-SED%lgZ0)/SED%dlgZ ) + 1,   1  ),  SED%nZ-1 )
        dz = SED%Zeds(SED%iz+1)-SED%Zeds(SED%iz)
        SED%dz0 = min( max(   (Z-SED%Zeds(SED%iz)) /dz,         0d0 ),  1d0         )
        SED%dz1 = min( max(  (SED%Zeds(SED%iz+1)-Z)/dz,         0d0 ),  1d0         )
     endif
     if (abs(SED%da0+SED%da1-1.0d0) > 1.0d-5 .or. abs(SED%dz0+SED%dz1-1.0d0) > 1.0d-5) then
        write(*,*) 'Screwed up the sed interpolation ... '
        write(*,*) SED%da0+SED%da1,SED%dz0+SED%dz1
        stop
     end if
  endif

  ret = SED%da0 * SED%dz0 * SED%table(SED%ia+1, SED%iz+1, :, nProp) + &
        SED%da1 * SED%dz0 * SED%table(SED%ia,   SED%iz+1, :, nProp) + &
        SED%da0 * SED%dz1 * SED%table(SED%ia+1, SED%iz,   :, nProp) + &
        SED%da1 * SED%dz1 * SED%table(SED%ia,   SED%iz,   :, nProp)

END SUBROUTINE inp_SED_table

!*************************************************************************
SUBROUTINE getNPhotonsEmitted(run, SED, age1_Gyr, dt_Gyr, Z, ret)

  ! Compute number of photons emitted by a stellar particle per solar mass
  ! over a timestep.
  ! input/output:
  ! age1_Gyr => Star population age [Gyrs] at timestep end
  ! dt_gyr   => Timestep length in Gyr
  ! Z        => Star population metallicity [m_metals/m_tot]
  ! ret      => # of photons emitted per solar mass over timestep
  !-------------------------------------------------------------------------
  type(run_t) :: run
  type(sed_table_t) :: SED
  real(kind=8),intent(in) :: age1_Gyr, dt_Gyr, Z
  real(kind=8),dimension(nrtgrp) :: ret
  !-------------------------------------------------------------------------
  real(kind=8),dimension(nrtgrp) :: Lc0, Lc1
  !-------------------------------------------------------------------------
  ! Lc0 = cumulative emitted photons at the start of the timestep
  call inp_SED_table(SED, age1_Gyr-dt_Gyr, Z, 2, .false., Lc0)
  ! Lc1 = cumulative emitted photons at the end of the timestep
  call inp_SED_table(SED, age1_Gyr, Z, 2, .false., Lc1)
  ret = max(Lc1-Lc0,0d0)

END SUBROUTINE getNPhotonsEmitted

#ifdef DO_RTZ
FUNCTION blackbody(T, lambda) result(B_lam)
  ! Blackbody function B_lam
  ! Harley addition so that it is easier to make default
  ! cross sections
  ! T --> temeprature [K]
  ! lambda --> wavelengths [A]
  use constants, only: c_cgs, hplanck, kB
  implicit none
  real(kind=8), intent(in):: T, lambda
  real(kind=8):: B_lam
  real(kind=8):: lambda_cm

  ! convert lambda in A to cm
  lambda_cm = 1.d-8 * lambda

  ! now compute B_lam
  B_lam = 2.d0 * hplanck * c_cgs * c_cgs / (lambda_cm**5.d0)
  B_lam = B_lam * (1.d0 / (exp(hplanck * c_cgs / (lambda_cm * kB * T)) - 1.d0))

END FUNCTION blackbody

SUBROUTINE initialize_cross_sections_from_blackbody(r, T, group_L0, group_L1, group_csn, group_cse, isH2_rtz)
  ! This subroutine initializes the cross sections of
  ! each species to be consistent with a blackbody of
  ! a given temperature
  use constants, only: c_cgs, hplanck, eV2erg
  use rt_parameters, only: nrtgrp
  use rtz_module, only: elements, n_elements
  implicit none
  type(run_t) :: r
  real(kind=8), intent(in):: T
  logical, intent(in):: isH2_rtz
  real(kind=8), intent(inout):: group_csn(nrtgrp,1:27,1:27), group_cse(nrtgrp,1:27,1:27)
  real(kind=8), intent(in):: group_L0(nrtgrp), group_L1(nrtgrp)
  real(kind=8):: lambda_min, lambda_max, delta_lambda, tmp
  real(kind=8):: X(1000), Y(1000)
  integer:: ip, ii, jj

  ! Loop over groups
  do ip = 1, nrtgrp
     ! No update for non-SED groups (L0>L1):
     if(group_L0(ip).ne. 0d0 .and. group_L1(ip) .ne. 0d0 .and. &
          &  (group_L0(ip) .ge. group_L1(ip)) ) cycle

     ! Fill out the X and Y arrays for integration
     lambda_max = (hplanck * c_cgs / (group_L0(ip)*eV2erg)) * 1d8 ! [A]
     lambda_min = (hplanck * c_cgs / (group_L1(ip)*eV2erg)) * 1d8 ! [A]
     delta_lambda = (lambda_max - lambda_min) / 999.d0

     ! Initialize X and Y arrays and fill them out
     X = 0.d0
     Y = 0.d0
     do ii=1, 1000
        X(ii) = lambda_min + (delta_lambda * (real(ii,kind=8)-1.d0))
        Y(ii) = blackbody(T, X(ii))
     end do

     ! Loop over elements
     do ii=1, n_elements
        ! Check if we actually use the element
        if (elements(ii)%atomic_number.gt.0) then
           ! Loop over ionization states
           do jj=1,elements(ii)%n_ions-1 !loop over ionization states
              group_csn(ip,ii,jj) = getSEDcsn(r, X, Y, 1000, group_L0(ip), group_L1(ip), ii, jj)
              group_cse(ip,ii,jj) = getSEDcse(r, X, Y, 1000, group_L0(ip), group_L1(ip), ii, jj)
           end do ! End loop over ionization states
        end if
     end do ! End loop over elements

     ! Deal with molecules separately
     if (elements(1)%atomic_number.gt.0 .and. isH2_rtz) then
        group_csn(ip,1,3) = getSEDcsn(r, X, Y, 1000, group_L0(ip), group_L1(ip), 1, 3)
        group_cse(ip,1,3) = getSEDcse(r, X, Y, 1000, group_L0(ip), group_L1(ip), 1, 3)
     end if
  end do ! End loop over groups

END SUBROUTINE initialize_cross_sections_from_blackbody

SUBROUTINE initialize_group_energies_from_blackbody(r, T, group_L0, group_L1, group_egy)
  ! This subroutine initializes the cross sections of
  ! each species to be consistent with a blackbody of
  ! a given temperature
  use constants, only: c_cgs, hplanck, eV2erg
  use rt_parameters, only: nrtgrp
  use rtz_module, only: elements, n_elements
  implicit none
  type(run_t) :: r
  real(kind=8), intent(in):: T
  real(kind=8), intent(in):: group_L0(nrtgrp), group_L1(nrtgrp)
  real(kind=8), intent(inout):: group_egy(nrtgrp)
  real(kind=8):: lambda_min, lambda_max, delta_lambda
  real(kind=8):: X(1000), Y(1000)
  integer:: ip, ii

  ! Loop over groups
  do ip = 1, nrtgrp
     ! No update for non-SED groups (L0>L1):
     if(group_L0(ip).ne. 0d0 .and. group_L1(ip) .ne. 0d0 .and. &
          &  (group_L0(ip) .ge. group_L1(ip)) ) cycle

     ! Fill out the X and Y arrays for integration
     lambda_max = (hplanck * c_cgs / (group_L0(ip)*eV2erg)) * 1d8 ! [A]
     lambda_min = (hplanck * c_cgs / (group_L1(ip)*eV2erg)) * 1d8 ! [A]
     delta_lambda = (lambda_max - lambda_min) / 1000.0

     ! Initialize X and Y arrays and fill them out
     X = 0.d0
     Y = 0.d0
     do ii=1, 1000
        X(ii) = lambda_min + (delta_lambda * (real(ii,kind=8)-1.d0))
        Y(ii) = blackbody(T, X(ii))
     end do

     group_egy(ip) = getSEDEgy(r, X, Y, 1000, group_L0(ip), group_L1(ip))

  end do

END SUBROUTINE initialize_group_energies_from_blackbody

#endif

END MODULE SED_module
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

!*************************************************************************
SUBROUTINE locate(xx,n,x,j)
  ! Locates position j of a value x in an ordered array xx of n elements
  ! After: xx(j) <= x <= xx(j+1) (assuming increasing order)
  ! j is lower bound, so it can be zero but not larger than n
  !-------------------------------------------------------------------------
  integer ::  n,j
  real(kind=8) :: xx(n),x
  !-------------------------------------------------------------------------
  integer :: jl,ju,jm
  !-------------------------------------------------------------------------
  jl = 0
  ju = n+1
  do while (ju-jl > 1)
     jm = (ju+jl)/2
     if ((xx(n) > xx(1)) .eqv. (x > xx(jm))) then
        jl = jm
     else
        ju = jm
     endif
  enddo
  j = jl
END SUBROUTINE locate

!*************************************************************************
SUBROUTINE inp_1d(xax,nx,x,ix0,ix1,dx0,dx1)
  ! Compute variables by interpolation from table with non-equal intervals.
  ! xax      => Axis of x-values in table
  ! nx       => Length of x-axis
  ! x        => x-value to interpolate to
  ! ix0,ix1 <=  Lower and upper boundaries of x in xax
  ! dx0,dx1 <=  Weights of ix0 and ix1 indexes
  !-------------------------------------------------------------------------
  integer:: nx, ix0, ix1
  real(kind=8), intent(in)::xax(nx), x
  real(kind=8):: dx0, dx1
  !-------------------------------------------------------------------------
  real(kind=8):: x_step
  !-------------------------------------------------------------------------
  call locate(xax, nx, x, ix0)
  if(ix0 < 1) ix0=1
  if (ix0 < nx) then
     ix1 = ix0+1
     x_step = xax(ix1) - xax(ix0)
     dx0 = max(           x - xax(ix0), 0.0d0 ) / x_step
     dx1 = min(xax(ix1) - x           , x_step) / x_step
  else
     ix1 = ix0
     dx0 = 0.0d0 ;  dx1 = 1.0d0
  end if

  if (abs(dx0+dx1-1.0d0) .gt. 1.0d-5) then
     write(*,*) 'Screwed up the 1d interpolation ... '
     write(*,*) dx0+dx1
     stop
  end if

END SUBROUTINE inp_1d
