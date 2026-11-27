module cr_params_module
contains
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine m_read_cr_params(pst)
  use amr_parameters
  use constants,only:c_cgs
  use hydro_parameters
  use cr_parameters
  use ramses_commons, only: pst_t
  use mdl_module
  use movie_module, only: set_movie_vars
  implicit none
  type(pst_t)::pst

#ifdef DO_CR
  !--------------------------------------------------
  ! Local variables
  !--------------------------------------------------
  character(LEN=flen)::infile
  logical::nml_ok

  !--------------------------------------------------
  ! CR namelist variables
  !--------------------------------------------------

  logical::cr_advect=.false.              ! Advection of cosmic rays?                       !
  logical::cr_streaming_diffusion=.false. ! Streaming diffusion of cosmic rays?             !
  logical::cr_streaming_heating=.false.   ! Streaming heating of cosmic rays?               !
  logical::cr_cooling=.false.             ! CR cooling?                                     !
  logical::cr_isotropic_pressure=.true.   ! Isotropic CR pressure?                          !
  logical::cr_reduced_flux_correction=.false.  ! Make sure F<c*E always?                    !
  real(kind=8)::cr_c_fraction=1.0
  real(kind=8)::cr_dmax=1d30              ! Max CR streaming diffusion coefficient in cgs   !
  integer::cr_nsubcycle=1                 ! Maximum number of CR subcycles per hydro step   !
  real(kind=8)::cr_courant_factor=0.8d0   ! Courant factor for CR timesteps                 !
  character(LEN=100)::cr_test_setup='none'! Setup for standard CR tests                     !
  real(kind=8),dimension(1:ncrgrp)::cr_d=1.0d29
  real(kind=8),dimension(1:ncrgrp)::cr_d_perp_factors=1d-6 ! perp CR diffusion suppression  !
  real(kind=8),dimension(1:ncrgrp)::fecr=0d0               ! SN fraction of CR energy       !
  real(kind=8)::cr_v_alfven=0.0           ! For idealised tests

  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------
  namelist/cr_params/ &
       &  cr_advect, cr_streaming_diffusion, cr_streaming_heating        &
       & ,cr_cooling, cr_isotropic_pressure, cr_reduced_flux_correction  &
       & ,cr_c_fraction, cr_dmax, cr_nsubcycle, cr_courant_factor        &
       & ,cr_test_setup, cr_v_alfven

  namelist/cr_groups/cr_d, cr_d_perp_factors, fecr

  associate(s=>pst%s)

  write(*,'(" Working with ",I2," CR groups ")') ncrgrp

  if(ncrgrp .le. 0) then
     s%r%cr = .false.
     write(*,'(" Turning off CRs, since no CR groups ")')
     return
  endif

  !--------------------------------------------------
  ! Set defaults for CR groups
  !--------------------------------------------------

  ! ...

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
  read(1,NML=cr_params,END=113)
113 continue
  rewind(1)
  read(1,NML=cr_groups,END=114)
114 continue
  close(1)

  ! Fill in all run parameters in corresponding structure
  s%r%cr_advect=cr_advect
  s%r%cr_streaming_diffusion=cr_streaming_diffusion
  s%r%cr_streaming_heating=cr_streaming_heating
  s%r%cr_cooling=cr_cooling
  s%r%cr_isotropic_pressure=cr_isotropic_pressure
  s%r%cr_reduced_flux_correction=cr_reduced_flux_correction
  s%r%cr_c_fraction=cr_c_fraction
  s%r%cr_dmax=cr_dmax
  s%r%cr_nsubcycle=cr_nsubcycle
  s%r%cr_courant_factor=cr_courant_factor
  s%r%cr_test_setup=cr_test_setup

  s%r%cr_d=cr_d
  s%r%cr_d_perp_factors=cr_d_perp_factors
  s%r%fecr=fecr
  s%r%cr_v_alfven=cr_v_alfven

  s%r%iecr=s%r%inener
  print*,'The index for cosmic ray energies is ', s%r%iecr
  if(ncrgrp>nener) then
     write(*,*)'There are not enough NENER variables to contain cosmic ray groups ', ncrgrp, nener
     call mdl_abort(s%mdl)
  endif
#ifdef _CUDA
  ! The GPU CR step runs without streaming, CR cooling and CR feedback, in 3D
  if(s%r%cr .and. (cr_streaming_diffusion .or. cr_streaming_heating))then
     write(*,*)'Error: cr_streaming_diffusion and cr_streaming_heating are not supported on GPU'
     call mdl_abort(s%mdl)
  endif
  if(s%r%cr .and. cr_cooling)then
     write(*,*)'Error: cr_cooling is not supported on GPU'
     call mdl_abort(s%mdl)
  endif
  if(s%r%cr .and. any(fecr/=0d0))then
     write(*,*)'Error: fecr is not supported on GPU'
     call mdl_abort(s%mdl)
  endif
  if(s%r%cr .and. ndim/=3)then
     write(*,*)'Error: cr is supported in 3D only on GPU'
     call mdl_abort(s%mdl)
  endif
#endif

  end associate
#endif

end subroutine m_read_cr_params
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
end module cr_params_module
