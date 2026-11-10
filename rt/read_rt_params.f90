module rt_params_module
contains
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine m_read_rt_params(pst)
  use amr_parameters
  use constants,only:c_cgs
  use hydro_parameters
  use rt_parameters
  use ramses_commons, only: pst_t
  use mdl_module
  use movie_module, only: set_movie_vars
#ifdef DO_RTZ
  use SED_module, only: initialize_cross_sections_from_blackbody, initialize_group_energies_from_blackbody
  use cross_sections_module, only: initialize_cross_sections
#endif
  implicit none
  type(pst_t)::pst

#ifdef DO_RT
  !--------------------------------------------------
  ! Local variables
  !--------------------------------------------------
  character(LEN=flen)::infile
  logical::nml_ok, rt_vsla=.false.

  !--------------------------------------------------
  ! RT namelist variables
  !--------------------------------------------------

  ! RT_PARAMS namelist
  logical::rt_advect=.false.           ! Advection of photons?                           !
  logical::rt_smooth=.true.            ! Smooth the discrete RT update of op. splitting  !
  logical::rt_star=.true.              ! Activate radiation from star particles          !
  logical::rt_sink=.false.             ! Activate radiation from sink particles          !
  logical::rt_emission_stats=.false.   ! Print stellar emission info                     !
  real(kind=8)::rt_esc_frac=1d0            ! Photon escape fraction from stellar particles   !
  !character(LEN=10)::rt_flux_scheme='glf'                                               !
  !logical::rt_use_hll=.false.         ! Use hll flux (or the default glf)               !
  !logical::rt_is_outflow_bound=.false.! Make all boundaries=outflow for RT              !
  real(kind=8)::rt_courant_factor=0.8d0    ! Courant factor for RT timesteps                 !
  real(kind=8),dimension(1:MAXLEVEL)::rt_c_fraction=1d0 ! Lightspeed fraction for RT         !
  !logical::rt_vsla=.false.            ! Are we using level variable light speed?        !
  integer::rt_nsubcycle=1              ! Maximum number of RT subcycles per hydro step   !
  logical::rt_otsa=.true.              ! Use on-the-spot approximation                   !
  logical::rt_isIR=.false.             ! Use IR photon groups                            !
  logical::is_kIR_T=.false.            ! Kappa IR depends on T_rad                       !
  logical::rt_T_rad=.false.            ! Use radiation temperature for everything        !
  logical::rt_isoPress=.false.         ! Use cE, not F, for rad. pressure                !
  real(kind=8)::rt_pressBoost=1d0          ! Boost on RT pressure                            !
  integer::iPEH_group=-1               ! Radiation group used for photo-electric heating !
  !logical::rt_isDiffuseUVsrc=.false.  ! UV emission from low-density cells              !
  !real(kind=8)::rt_UVsrc_nHmax=-1d0       ! Density threshold for UV emission               !
  !logical::upload_equilibrium_x=.true.! Enforce equilibrium xion when uploading         !
  !integer::heat_unresolved_HII=0      ! Subgrid model heating unresolved HII regions    !
  !integer::iHIIheat=6                 ! Var index for HII heating                       !
  !logical::cosmic_rays=.false.        ! Include cosmic ray ionisation                   !

  !character(LEN=128)::hll_evals_file=''! File HLL eigenvalues                           !
  character(LEN=128)::sed_dir=''       ! Dir containing stellar energy distributions     !
  !character(LEN=128)::uv_file=''      ! File containing stellar energy distributions    !

  ! RT source regions parameters----------------------------------------------------------
  integer                           ::rt_nsource=0
  character(LEN=10),dimension(1:MAXREGION)::rt_source_type='square'
  real(kind=8),dimension(1:MAXREGION)   ::rt_src_x_center=0.
  real(kind=8),dimension(1:MAXREGION)   ::rt_src_y_center=0.
  real(kind=8),dimension(1:MAXREGION)   ::rt_src_z_center=0.
  real(kind=8),dimension(1:MAXREGION)   ::rt_src_length_x=1.E10
  real(kind=8),dimension(1:MAXREGION)   ::rt_src_length_y=1.E10
  real(kind=8),dimension(1:MAXREGION)   ::rt_src_length_z=1.E10
  real(kind=8),dimension(1:MAXREGION)   ::rt_exp_source=2.0
  integer, dimension(1:MAXREGION)   ::rt_src_group=1
  !integer, dimension(1:MAXREGION)   ::rt_src_trace_group=1
  real(kind=8),dimension(1:MAXREGION)   ::rt_n_source=0.                      ! Photon density
  real(kind=8),dimension(1:MAXREGION)   ::rt_u_source=0.                      !    Photon flux
  real(kind=8),dimension(1:MAXREGION)   ::rt_v_source=0.                      !    Photon flux
  real(kind=8),dimension(1:MAXREGION)   ::rt_w_source=0.                      !    Photon flux

  ! RT_GROUPS namelist---------------------------------------------------------------------
  integer::sedprops_update=-1           ! Update sedprops from stellar populations        !
  ! negative: never update, 0:update on init, pos x: update every x coarse steps
  ! logical::SED_isEgy=.false. ! Integrate energy out of SEDs rather than photon count

  ! Group props: avg and energy weigthed photoionization c-section (cm2), avg. energy (ev)
  ! Indices nrtgrp, nion stand for photon group vs species (e.g. 1=H, 2=He).
#ifdef DO_RTZ
  real(kind=8),dimension(nrtgrp,1:27,1:27)::group_csn=0                   !    Cross sections (cm2)
  real(kind=8),dimension(nrtgrp,1:27,1:27)::group_cse=0                   !    Cross sections (cm2)
#else
  real(kind=8),dimension(nrtgrp,nion)::group_csn=0                   !    Cross sections (cm2)
  real(kind=8),dimension(nrtgrp,nion)::group_cse=0                   !    Cross sections (cm2)
#endif
  real(kind=8),dimension(nrtgrp)::group_egy=0                        !  Avg photon energy (ev)
  real(kind=8),dimension(nrtgrp)::group_L0=13.60                     ! Wavelength lower limits
  real(kind=8),dimension(nrtgrp)::group_L1=0                         ! Wavelength upper limits
  real(kind=8),dimension(nrtgrp)::kappaAbs=0                         ! Dust absorption opacity
  real(kind=8),dimension(nrtgrp)::kappaSc=0                          ! Dust scattering opacity
  integer,dimension(nion)::spec2group=0                 ! Ion -> group # in recombinations

  integer::i, igroup_HI=0, igroup_HII=0, igroup_HeII=0, igroup_HeIII=0

  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------
  namelist/rt_params/rt_advect, rt_otsa, rt_c_fraction, rt_nsubcycle     &
       & ,rt_smooth, rt_isIR, is_kIR_T, rt_T_rad, rt_isoPress            &
       & ,rt_pressBoost, iPEH_group
  namelist/rt_sources/rt_star, rt_sink, rt_esc_frac                      &
       & ,rt_emission_stats ,rt_nsource, rt_source_type                  &
       & ,rt_src_x_center, rt_src_y_center, rt_src_z_center              &
       & ,rt_src_length_x, rt_src_length_y, rt_src_length_z              &
       & ,rt_exp_source, rt_src_group                                    &
       & ,rt_n_source, rt_u_source, rt_v_source, rt_w_source
  namelist/rt_groups/group_csn, group_cse, group_egy, spec2group         &
       & ,group_L0, group_L1, kappaAbs, kappaSc, sed_dir, sedprops_update

  associate(s=>pst%s, ionEvs=>pst%s%r%ionEvs, isH2=>pst%s%r%isH2         &
    ,isHe=>pst%s%r%isHe, ixHI=>pst%s%r%ixHI, ixHII=>pst%s%r%ixHII        &
    ,ixHeII=>pst%s%r%ixHeII, ixHeIII=>pst%s%r%ixHeIII)

  write(*,'(" Working with ",I2," photon groups ")') nrtgrp

  if(nrtgrp .le. 0) then
     s%r%rt = .false.
     write(*,'(" Turning off RT, since no RT groups ")')
     return
  endif

  !--------------------------------------------------
  ! Set defaults for radiation groups
  !--------------------------------------------------
#if NRTGRP>0
#ifndef DO_RTZ
  !  Use H2, HI, HeI, HeII ionization energies  as default group intervals
  group_L0(1:min(nrtgrp,nion))=ionEvs(1:min(nrtgrp,nion))   ! Lower bounds
  group_L1(1:min(nrtgrp,nion-1))=ionEvs(2:min(nrtgrp+1,nion)) ! Upper bnds
  group_L1(min(nrtgrp,nion))=0.                     ! Upper bound=infinity

  i=0
  if(isH2) then ! Set index for H2 dissociating group
     i=i+1 ; igroup_HI=i
  endif
  if(i .lt. nrtgrp) then ! Set index for HI ionizing group
     i=i+1 ; igroup_HII=i
  endif
  if(i .lt. nrtgrp .and. isHe) then ! Set index for HeI ionizing group
     i=i+1 ; igroup_HeII=i
  endif
  if(i .lt. nrtgrp .and. isHe) then ! Set index for HeII ionizing group
     i=i+1 ; igroup_HeIII=i
  endif

  ! Default groups are all blackbodies at 10^5 Kelvin:
  group_csn=0d0 ; group_cse=0d0 ; group_egy=0d0         ! Default all zero
  if(igroup_HI .gt. 0) then
     if(ixHI .gt. 0) then                                ! H2 dissociation
        group_csn(igroup_HI,ixHI)=2.1d-19
        group_cse(igroup_HI,ixHI)=2.1d-19
     endif
     group_egy(igroup_HI)=12.44
  endif
  if(igroup_HII .gt. 0) then
     if(ixHI .gt. 0) then                   ! H2 ionization by HI photons
        group_csn(igroup_HII,ixHI)=5.0d-18
        group_cse(igroup_HII,ixHI)=5.3d-18
     endif
     group_csn(igroup_HII,ixHII)=3.007d-18                ! HI ionization
     group_cse(igroup_HII,ixHII)=2.781d-18
     group_egy(igroup_HII)=18.85
  endif
  if(igroup_HeII .gt. 0) then
     if(ixHI .gt. 0) then                  ! H2 ionization by HeI photons
        group_csn(igroup_HeII,ixHI)=2.9d-18
        group_cse(igroup_HeII,ixHI)=2.8d-18
     endif
     group_csn(igroup_HeII,ixHII)=5.687d-19! HI ionization by HeI photons
     group_cse(igroup_HeII,ixHII)=5.042d-19
     if(ixHeII .gt. 0) then                              ! HeI ionization
        group_csn(igroup_HeII,ixHeII)=4.478d-18
        group_cse(igroup_HeII,ixHeII)=4.130d-18
     endif
     group_egy(igroup_HeII)=35.079
  endif
  if(igroup_HeIII .gt. 0) then
     if(ixHI .gt. 0) then                 ! H2 ionization by HeII photons
        group_csn(igroup_HeIII,ixHI)=4.1d-19
        group_cse(igroup_HeIII,ixHI)=4.1d-19
     endif
     group_csn(igroup_HeIII,ixHII)=7.889d-20  ! HI ioniz. by HeII photons
     group_cse(igroup_HeIII,ixHII)=7.456d-20
     if(ixHeII .gt. 0) then                  ! HeI ioniz. by HeII photons
        group_csn(igroup_HeIII,ixHeII)=1.197d-18
        group_cse(igroup_HeIII,ixHeII)=1.142d-18
     endif
     if(ixHeIII .gt. 0) then                            ! HeII ionization
        group_csn(igroup_HeIII,ixHeIII)=1.055d-18
        group_cse(igroup_HeIII,ixHeIII)=1.001d-18
     endif
     group_egy(igroup_HeIII)=65.666
  endif
#endif
#endif

  do i=1,min(nion,nrtgrp)
     spec2group(i)=i                   ! Species contributions to groups
  end do

  !-------------------------------------------------
  ! Read the namelist file
  !-------------------------------------------------
  CALL getarg(1,infile)
  namelist_file=TRIM(infile)
  INQUIRE(file=infile,exist=nml_ok)
  if(.not. nml_ok)then
     write(*,*)'File '//TRIM(infile)//' does not exist'
     call mdl_abort(s%mdl)
  end if
  open(1,file=infile)
  rewind(1)
  read(1,NML=rt_params,END=113)
113 continue
  rewind(1)
  read(1,NML=rt_groups,END=114)
114 continue
  rewind(1)
  read(1,NML=rt_sources,END=115)
115 continue
  close(1)

  !--------------------------------------------------
  ! Check the IR and photoelectric heating options
  !--------------------------------------------------
  if(iPEH_group>nrtgrp)then
     write(*,*)'Error: iPEH_group must be at most nrtgrp, not ',iPEH_group
     call mdl_abort(s%mdl)
  endif
  if(is_kIR_T.and.nrtgrp>1.and.kappaAbs(1)<=0d0)then
     write(*,*)'Error: is_kIR_T=.true. with nrtgrp>1 needs kappaAbs(1)>0'
     call mdl_abort(s%mdl)
  endif

#ifdef _CUDA
  !--------------------------------------------------
  ! Check the RT options supported on GPU
  !--------------------------------------------------
  ! The GPU chemistry has no IR dust coupling
  if(rt_isIR)then
     write(*,*)'Error: rt_isIR is not supported on GPU'
     call mdl_abort(s%mdl)
  endif
  ! The GPU chemistry has no IR opacity that depends on the radiation temperature
  if(is_kIR_T.and.(rt_advect.or.rt_T_rad.or.(s%r%star.and.rt_star)))then
     write(*,*)'Error: is_kIR_T with rt_advect=.true., rt_T_rad=.true. or star=.true. with rt_star=.true. is not supported on GPU'
     call mdl_abort(s%mdl)
  endif
  ! Without rt_advect the GPU cooling runs, and it updates the photon variables only with neq_chem
  if(s%r%isothermal.and..not.s%r%neq_chem.and..not.rt_advect)then
     write(*,*)'Error: rt with isothermal and rt_advect=.false. but without neq_chem is not supported on GPU'
     call mdl_abort(s%mdl)
  endif
#endif

  !--------------------------------------------------
  ! Reduced light speed. First check if only one light speed fraction set
  ! and if so, set that same light speed fraction at all levels. If more
  ! than one fraction is set, we are using a variable speed of light.
  !--------------------------------------------------
  if(rt_c_fraction(1).ne.rt_c_fraction(2) .and. rt_c_fraction(2).eq.1.0  &
                                .and. all((rt_c_fraction(2:)).eq.1.0)) then
      rt_c_fraction(2:) = rt_c_fraction(1)
  endif

  ! Shift lightspeed fractions from levels [1:] to [levelmin:]
  do i=s%r%nlevelmax,s%r%levelmin,-1
     rt_c_fraction(i)=rt_c_fraction(i-s%r%levelmin+1)
  end do
  do i=1,s%r%levelmin-1     ! Just a dummy lightspeed for non-leaf levels
     rt_c_fraction(i)=1.0
  end do
  !do i=s%r%nlevelmax,s%r%levelmin,-1 !Set the light speed(s) according to f_c
  !   s%g%rt_c_cgs(i) = c_cgs * rt_c_fraction(i)
  !end do

  ! Print a message if using level-variable speed of light
  do i=s%r%levelmin,s%r%nlevelmax-1
     if(rt_c_fraction(i) .ne. rt_c_fraction(s%r%nlevelmax)) then
        write(*,212) rt_c_fraction(s%r%levelmin:s%r%nlevelmax)
        rt_vsla=.true.
        exit
     endif
  end do
  if(.not. rt_vsla) write(*,213) rt_c_fraction(s%r%levelmin)

#if NRTGRP>0
#ifdef DO_RTZ
  ! in the case of RTZ, perform initialization after reading in group
  ! energies

  ! Frist initialize the cross sections data
  call initialize_cross_sections()

  ! Initialize cross sections to be a blackbody at 1e5 K
  call initialize_cross_sections_from_blackbody(s%r, 1.d5, group_L0, group_L1, group_csn, group_cse, .true.)

  ! Initialize group energies for the same black body
  call initialize_group_energies_from_blackbody(s%r, 1.d5, group_L0, group_L1, group_egy)

#endif
#endif

  ! Fill in all run parameters in corresponding structure
  s%r%rt_otsa=rt_otsa
  s%r%rt_advect=rt_advect
  s%r%rt_smooth=rt_smooth
  s%r%rt_c_fraction=rt_c_fraction
  s%r%rt_nsubcycle=rt_nsubcycle
  s%r%rt_courant_factor=rt_courant_factor
  s%r%rt_isIR=rt_isIR
  s%r%is_kIR_T=is_kIR_T
  s%r%rt_T_rad=rt_T_rad
  s%r%rt_isoPress=rt_isoPress
  s%r%rt_pressBoost=rt_pressBoost
  s%r%iPEH_group=iPEH_group

  s%r%rt_star=rt_star
  s%r%rt_sink=rt_sink
  s%r%rt_esc_frac=rt_esc_frac
  s%r%sed_dir=sed_dir
  s%r%rt_emission_stats=rt_emission_stats
  s%r%rt_nsource=rt_nsource
  s%r%rt_source_type=rt_source_type
  s%r%rt_src_x_center=rt_src_x_center
  s%r%rt_src_y_center=rt_src_y_center
  s%r%rt_src_z_center=rt_src_z_center
  s%r%rt_src_length_x=rt_src_length_x
  s%r%rt_src_length_y=rt_src_length_y
  s%r%rt_src_length_z=rt_src_length_z
  s%r%rt_exp_source=rt_exp_source
  s%r%rt_src_group=rt_src_group
  s%r%rt_n_source=rt_n_source
  s%r%rt_u_source=rt_u_source
  s%r%rt_v_source=rt_v_source
  s%r%rt_w_source=rt_w_source

  s%r%sed_dir=sed_dir
  s%r%sedprops_update=sedprops_update
  s%r%group_csn=group_csn
  s%r%group_cse=group_cse
  s%r%group_egy=group_egy
  s%r%spec2group=spec2group
  s%r%group_L0=group_L0
  s%r%group_L1=group_L1
  s%r%kappaAbs=kappaAbs
  s%r%kappaSc=kappaSc

  if(minval(s%r%group_egy) .le. 0d0) then
     write(*,*) '========================================================='
     write(*,*) 'WARNING! Some photon groups have zero or negative energy!'
     write(*,*) 'This could have unwanted effects, so be careful!!!'
     write(*,*) '========================================================='
  endif

  if(s%r%isH2.or.s%r%isH2_rtz) then
     do i=1,nrtgrp
        if(    (s%r%group_L0(i) .ge. 11.2) .and. (s%r%group_L1(i) .le. 13.6) .and. &
             & (s%r%group_L0(i) .le. 13.6) .and. (s%r%group_L1(i) .ge. 11.2) )then
           s%r%ssh2(i) = 4d2 ! H2 self-shielding factor
           s%r%isLW(i) = 1d0 ! Index for LW groups
        endif
     enddo
  endif

  end associate

212 format (' Using a level-variable speed of light, with f_c= '20(1pe12.3))
213 format (' Using a uniform reduced speed of light fraction of f_c='1pe10.3)
#endif
end subroutine m_read_rt_params
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
end module rt_params_module
