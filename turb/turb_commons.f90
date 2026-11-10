module turb_commons
  use amr_parameters, only: ndim, nvector, twotondim
  use amr_commons, only: run_t

  integer, parameter :: ILP = selected_int_kind(r=15) ! integer long precision

  ! Turbulence variables
  integer, parameter :: TURB_GS=64                    ! Turbulent grid size
  integer, parameter :: TGRID_X=TURB_GS-1             ! Limit of grid, x dimension
#if NDIM>1
  integer, parameter :: TGRID_Y=TURB_GS-1             ! Limit of grid, y dimension
#else
  integer, parameter :: TGRID_Y=0                     ! Limit of grid, y dimension
#endif
#if NDIM>2
  integer, parameter :: TGRID_Z=TURB_GS-1             ! Limit of grid, z dimension
#else
  integer, parameter :: TGRID_Z=0                     ! Limit of grid, z dimension
#endif
  character(len=16), parameter :: precision_str='DOUBLE_PRECISION'
  real(kind=8), parameter :: turb_gs_real=real(TURB_GS,kind=8)

  ! Turbulence main object
  type turb_t

     complex(kind=8), allocatable :: turb_last(:,:,:,:)   ! Turbulent spectrum at time = t
     complex(kind=8), allocatable :: turb_next(:,:,:,:)   ! Turbulent spectrum at time = t + dt
     real(kind=8), allocatable :: power_spec(:,:,:)    ! Power spectrum of turbulence

     real(kind=8), allocatable :: afield_last(:,:,:,:) ! Forcing field at time = t
     real(kind=8), allocatable :: afield_next(:,:,:,:) ! Forcing field at time = t + dt
     real(kind=8), allocatable :: afield_now(:,:,:,:)  ! Forcing field now

     real(kind=8) :: comp_amp         ! Longitudinal (compressive) mode amplitude
     real(kind=8) :: sol_amp          ! Amplitude of each transverse (solenoidal) mode
     real(kind=8) :: turb_last_time   ! Time of old turbulent field
     real(kind=8) :: turb_next_time   ! Time of next turbulent field
     real(kind=8) :: turb_dt          ! Turbulent velocity evolution timestep
     real(kind=8) :: turb_decay_frac  ! Decay fraction per dt
     real(kind=8) :: turb_space(1:3)  ! Grid spacing
     real(kind=8) :: turb_norm        ! Normalizing constant from combination
                                      ! of Ornstein-Uhlenbeck process, initial
                                      ! power spectrum and projection

     character(len=256) :: turb_file_last    ! filename for 'last' field
     character(len=256) :: turb_file_next    ! filename for 'last' field
     character(len=256) :: turb_file_header  ! filename for header file

     integer(ILP) :: kiss64_state(1:4) ! State variables for KISS64 PRNG (Marsaglia)

  end type turb_t

contains
#if NDIM==3
#define HERMITIAN_FIELD 1
#endif
#define PI 3.141592653589793
  !=====================================================================================
  !=====================================================================================
  !=====================================================================================
  !=====================================================================================
  subroutine read_turb_fields(run, turb)
    implicit none
    type(run_t)        :: run
    type(turb_t)       :: turb

    integer            :: ilun, nturbtemp_tmp
    character(len=256) :: filedir
    logical            :: ok
    character(len=5)   :: nchar
    character(len=16)  :: precision_str_tmp
    real(kind=8)       :: aux
    character(len=17)  :: turb_last_time_z ! turb_last_time in hex
    character(len=17)  :: turb_next_time_z ! turb_next_time in hex

    ilun = 100
    call title(run%nrestart,nchar)
    filedir = 'backup_'//TRIM(nchar)//'/'
    turb%turb_file_header = trim(filedir)//'turb_fields.dat'
    inquire(file=turb%turb_file_header, exist=ok)
    if (.not. ok)then
       write(*,*)'Restart failed:'
       write(*,*)'File '//TRIM(turb%turb_file_header)//' not found'
       stop
    end if

    ! Read header file of important information
    open(ilun,file=turb%turb_file_header,status="old",form="formatted")
    read(ilun,*) nturbtemp_tmp ! Not used, SEREN compatibility
    read(ilun,*) precision_str_tmp
    read(ilun,*) turb%turb_file_last
    read(ilun,*) aux, turb_last_time_z
    read(ilun,*) turb%turb_file_next
    read(ilun,*) aux, turb_next_time_z
    read(ilun,*); read(ilun,*); read(ilun,*); read(ilun,*); read(ilun,*)
    read(ilun,*) turb%kiss64_state
    close(ilun)
    ! We don't care about some of the stuff in the file
    ! as it should be set in the params file

    read(turb_last_time_z, '(Z16)') turb%turb_last_time
    read(turb_next_time_z, '(Z16)') turb%turb_next_time

    ! Read turbulent fields
    open(ilun,file=trim(filedir)//trim(turb%turb_file_last),status="unknown",access="stream",form="unformatted")
    read(ilun) turb%turb_last
    close(ilun)

    open(ilun,file=trim(filedir)//trim(turb%turb_file_next),status="unknown",access="stream",form="unformatted")
    read(ilun) turb%turb_next
    close(ilun)

  end subroutine read_turb_fields
  !=====================================================================================
  !=====================================================================================
  !=====================================================================================
  !=====================================================================================
  subroutine write_turb_fields(run, turb, output_dir)
    implicit none
    type(run_t)                  :: run
    type(turb_t)                 :: turb
    character(len=*), intent(in) :: output_dir  ! Output directory

    integer                      :: ilun        ! File I/O unit
    character(len=50000)         :: file_buffer ! Buffer for header file
    character(len=1)             :: c           ! Mostly pointless variable
    character(len=17)            :: turb_last_time_z ! turb_last_time in hex
    character(len=17)            :: turb_next_time_z ! turb_next_time in hex

    ilun = 10

    turb%turb_file_last = 'turb_last.dat'
    turb%turb_file_next = 'turb_next.dat'
    turb%turb_file_header = trim(output_dir)//'turb_fields.dat'

    write(turb_last_time_z, '(X,Z16)') turb%turb_last_time
    write(turb_next_time_z, '(X,Z16)') turb%turb_next_time

    ! Write turbulent fields
    open(ilun,file=trim(output_dir)//trim(turb%turb_file_last),status="unknown",access="stream",form="unformatted")
    write(ilun) turb%turb_last
    close(ilun)

    open(ilun,file=trim(output_dir)//trim(turb%turb_file_next),status="unknown",access="stream",form="unformatted")
    write(ilun) turb%turb_next
    close(ilun)

    ! Write header file of important information
    write(file_buffer,*) 0, NEW_LINE(c),&
         & precision_str, NEW_LINE(c),&
         & trim(turb%turb_file_last), NEW_LINE(c),&
         & turb%turb_last_time, turb_last_time_z, NEW_LINE(c),&
         & trim(turb%turb_file_next), NEW_LINE(c),&
         & turb%turb_next_time, turb_next_time_z, NEW_LINE(c),&
         & TURB_GS, TURB_GS, TURB_GS, NEW_LINE(c),&
         & (/0.d0, 0.d0, 0.d0/), NEW_LINE(c),&      ! turb_min
         & (/1.d0, 1.d0, 1.d0/), NEW_LINE(c),&      ! turb_max
         & run%turb_T, turb%turb_dt, NEW_LINE(c),&
         & run%comp_frac, 1.d0 - run%comp_frac, NEW_LINE(c),&  ! compressive, solenoidal power
         & turb%kiss64_state, NEW_LINE(c)                ! random number state

    ! Now write atomically (safest)
    open(ilun,file=turb%turb_file_header,status="unknown",access='stream',form="formatted")
    write(ilun,*) trim(file_buffer)
    close(ilun)

  end subroutine write_turb_fields
  !=====================================================================================
  !=====================================================================================
  !=====================================================================================
  !=====================================================================================
  subroutine turb_next_field(run, turb)
    implicit none
    type(run_t) :: run
    type(turb_t) :: turb

    ! Set turbulent field times
    turb%turb_last_time = turb%turb_next_time
    turb%turb_next_time = turb%turb_last_time + turb%turb_dt

    ! Copy current field, and add more turbulence
    turb%turb_last = turb%turb_next
    turb%afield_last = turb%afield_next
    call add_turbulence(turb, turb%turb_next, turb%turb_dt)

    ! Subtract decay of turbulence
    call decay_turbulence(turb%turb_decay_frac, turb%turb_last, turb%turb_next)

    ! Fourier transform
#if NDIM==1
    call FFT_1D(turb%turb_next(1,:,0,0), turb%afield_next(1,:,0,0))
#elif NDIM==2
    call FFT_2D(turb%turb_next(:,:,:,0), turb%afield_next(:,:,:,0))
#else
    call FFT_3D(turb%turb_next, turb%afield_next)
#endif
    turb%afield_next = turb%afield_next * turb%turb_norm * run%turb_rms

  end subroutine turb_next_field
  !=====================================================================================
  !=====================================================================================
  !=====================================================================================
  !=====================================================================================
subroutine find_conj_pair(i,j,k,ii,jj,kk)
    implicit none
    integer, intent(in)  :: i,j,k         ! Grid coordinates
    integer, intent(out) :: ii,jj,kk      ! Hermitian pair

    ii = TURB_GS - i
    if (ii >= TURB_GS) ii = ii - TURB_GS
    jj = TURB_GS - j
    if (jj >= TURB_GS) jj = jj - TURB_GS
    kk = TURB_GS - k
    if (kk >= TURB_GS) kk = kk - TURB_GS

  end subroutine find_conj_pair
  !=====================================================================================
  !=====================================================================================
  !=====================================================================================
  !=====================================================================================
  subroutine find_unitk(i,j,k,limit,unitk)
    implicit none
    integer, intent(in)        :: i,j,k,limit
    integer                    :: ii,jj,kk
    real (kind=8), intent(out) :: unitk(1:NDIM)

    ii = i
    if (i > limit) ii = i - TURB_GS
#if NDIM>1
    jj = j
    if (j > limit) jj = j - TURB_GS
#endif
#if NDIM>2
    kk = k
    if (k > limit) kk = k - TURB_GS
#endif

#if NDIM==1
    unitk = (/real(ii,kind=8)/)
#elif NDIM==2
    unitk = (/real(ii,kind=8),real(jj,kind=8)/)
#else
    unitk = (/real(ii,kind=8),real(jj,kind=8),real(kk,kind=8)/)
#endif
#if NDIM>1
    unitk = unitk / sqrt(sum(unitk**2))
#endif

  end subroutine find_unitk
  !=====================================================================================
  !=====================================================================================
  !=====================================================================================
  !=====================================================================================
  subroutine calc_power_spectrum(run, k, power_spectrum)
    implicit none
    type(run_t)               :: run
    integer, intent(in)       :: k(1:3)         ! Wavevector
    real(kind=8), intent(out) :: power_spectrum ! Power value
    real(kind=8)              :: k_mag          ! Wavevector magnitude

    ! Remark that the components of k are between -TURB_GS and TURB_GS
    ! with k=1 corresponding to the box size
    select case(run%forcing_power_spectrum)

    case('power_law')
       ! alpha^-2 power spectrum
       if (all(k==0)) then
          power_spectrum = 0
          return
       end if

       k_mag = sqrt(real(sum(k**2),kind=8))
       if (k_mag > (TURB_GS/2)) then
          power_spectrum = 0
          return
       end if
       power_spectrum = k_mag**(-2)

    case('parabolic')
       ! 'parabola' large-scale modes power spectrum
       power_spectrum = 0.
       k_mag = sqrt(real(sum(k**2),kind=8))
       if ((k_mag > 1.0) .AND. (k_mag < 3.0)) then
          power_spectrum = 1.0 - (k_mag-2.0)**2
       end if

    case('konstandin')
       ! forcing between k=1 (max) and k=1 (zero) as in Konstandin 2015
       power_spectrum = 0.
       k_mag = sqrt(real(sum(k**2),kind=8))
       if ((k_mag >= 0.999999999999999) .AND. (k_mag < 2.0)) then
          power_spectrum = 2.0 - (k_mag)
       end if

    case('test')
       power_spectrum = 0.
       if (k(1)==1 .AND. k(2)==0 .AND. k(3)==0) then
          power_spectrum = 1.0
       end if

    case default
       write (6,*) "Unknown forcing_power_spectrum!"
       write (6,*) "Use 'power_law', 'parabolic', 'konstandin' or 'test'"
       stop
    end select

  end subroutine calc_power_spectrum
  !=====================================================================================
  !=====================================================================================
  !=====================================================================================
  !=====================================================================================
  subroutine gaussian_cmplx(turb, G)
    implicit none
    ! Generates 3 complex random variates, where each variate has a
    ! complex value drawn from a Gaussian distribution EITHER
    ! by assigning a Gaussian magnitude and uniformly distributed argument OR
    ! by assigning Gaussian real and imaginary components
    type(turb_t) :: turb
    complex(kind=8), intent(out) :: G(1:ndim) ! Random gaussian values
    ! Local variables
    integer      :: d           ! Dimension counter
    real(kind=8) :: Rnd(1:3)    ! Random numbers
    real(kind=8) :: w           ! Aux. variable for random Gaussian
    real(kind=8) :: mag(1:ndim), arg(1:ndim) ! Magnitude and argument
    ! Create Gaussian distributed random numbers
    ! Marsaglia version of the Box-Muller algorithm
    ! This creates a Gaussian with a spread of 1
    ! and centred on 0 (hopefully)
    ! You get two independent variates each time

    do d=1,ndim
       do
          call kiss64_double(2, turb%kiss64_state, Rnd(1:2))
          Rnd(1) = 2.0 * Rnd(1) - 1.0
          Rnd(2) = 2.0 * Rnd(2) - 1.0
          w = Rnd(1)**2 + Rnd(2)**2
          if (w<1.0) exit
       end do
       w = sqrt( (-2.0 * log( w ) ) / w )

       mag(d) = Rnd(1) * w
       ! Throwing away second random variate (something of a waste)
    end do

    call kiss64_double(ndim, turb%kiss64_state, Rnd(1:ndim))
    arg = (Rnd(1:ndim) * 2.0 * PI) - 1.0

    G = cmplx(mag * cos(arg), mag * sin(arg), kind=8)

  end subroutine gaussian_cmplx
  !=====================================================================================
  !=====================================================================================
  !=====================================================================================
  !=====================================================================================
  subroutine add_turbulence(turb, turb_field, dt)
    implicit none
    ! Take a complex Fourier field of turbulence, and add a random component
    ! from a Wiener process (in this case Gaussian perturbations) multiplied
    ! by an initial distribution of power and then projected by Helmholtz
    ! decomposition.
    ! Optionally, this field can be Hermitian i.e. purely real after the
    ! transform, but I'm not sure if this produces a purely even field...
    type(turb_t)                   :: turb
    complex(kind=8), intent(inout) :: turb_field(1:NDIM,0:TGRID_X,0:TGRID_Y,0:TGRID_Z)
                                           ! Complex field to add to
    real(kind=8), intent(in)       :: dt   ! Width of Gaussian which drives
                                           ! Wiener process
    integer              :: i, j, k         ! Loop variables
    integer              :: halfway         ! Half of grid size

#if defined(HERMITIAN_FIELD)
    integer              :: ii,jj,kk        ! Hermitian pair
    logical              :: hermitian_pair  ! Do we need to take conjugate?
    logical              :: own_conjg       ! Are we are own conjugate pair?
#endif

    complex(kind=8)    :: complex_vec(1:NDIM) ! Complex Fourier vector
#if NDIM>1
    real(kind=8)       :: unitk(1:NDIM)       ! Unit vector parallel to k
    complex(kind=8)    :: unitk_cmplx(1:NDIM) ! Unit vector parallel to k
    complex(kind=8)    :: comp_cmplx(1:NDIM)  ! Compressive component
    complex(kind=8)    :: sol_cmplx(1:NDIM)   ! Solenoidal components
#endif

    halfway = TURB_GS / 2 ! integer division

    do k = 0, TGRID_Z
       do j = 0, TGRID_Y
          do i = 0, TGRID_X
#if defined(HERMITIAN_FIELD)
             hermitian_pair = .FALSE.
#endif
             ! Check there is any power in this mode, else cycle
             if (turb%power_spec(i,j,k) == 0.0) cycle

#if defined(HERMITIAN_FIELD)
             ! Test if we are own conjugate or need to use another conjugate
             own_conjg = .FALSE.
             call find_conj_pair(i,j,k,ii,jj,kk)
             if (i==ii .AND. j==jj .AND. k==kk) own_conjg = .TRUE.
             if (k > halfway) then
                hermitian_pair = .TRUE.
             else if (k == 0 .OR. 2*k == TURB_GS) then
                if (j > halfway) then
                   hermitian_pair = .TRUE.
                else if (j == 0 .OR. 2*j == TURB_GS) then
                   if (i > halfway) hermitian_pair = .TRUE.
                end if
             end if
             if (hermitian_pair) then
                turb_field(1:3,i,j,k) = &
                     & conjg(turb_field(1:3,ii,jj,kk))
                cycle
             end if
#endif

             ! Ornstein-Uhlenbeck process
             ! dF(k,t) = F_0(k) P(k) dW - F(k,t) dt / T
             ! where F(k,t) is the vector fourier amplitude,
             ! F_0 is the power spectrum distribution of amplitudes,
             ! P(k) is the projection (Helmholtz decomposition),
             ! dt is the timestep of integration,
             ! and dW is the Wiener process, described by a random variate
             ! from a normal distribution with a mean of zero and a
             ! variance of dt i.e. a standard deviation of sqrt(dt)

             ! In this we only calculate the first term - the random
             ! growth term.

             ! Random power/phase (complex Gaussian) multiplied by power
             ! spectrum and multiplied by width dt
             call gaussian_cmplx(turb, complex_vec)
             complex_vec = complex_vec * sqrt(dt) * turb%power_spec(i,j,k)

#if defined(HERMITIAN_FIELD)
             if (own_conjg) then
                ! To be own conjugate, phase must be zero
                ! (in order to be real)
                complex_vec = abs(complex_vec)
             end if
#endif

#if NDIM>1
             ! Helmholtz decomposition!
             call find_unitk(i, j, k, halfway, unitk)
             unitk_cmplx = cmplx(unitk, kind=8)
             comp_cmplx = unitk_cmplx *  dot_product(unitk_cmplx,complex_vec)
             sol_cmplx = complex_vec - comp_cmplx
             ! comp_amp and sol_amp make comp_frac the compressive share
             ! of the forcing power (see init_turb)
             complex_vec = comp_cmplx*turb%comp_amp + sol_cmplx*turb%sol_amp
             turb_field(1:NDIM,i,j,k) = turb_field(1:NDIM,i,j,k) + &
                  & complex_vec
#else
             turb_field(1:NDIM,i,j,k) = turb_field(1:NDIM,i,j,k) + &
                  & complex_vec
#endif

          end do
       end do
    end do

  end subroutine add_turbulence
  !=====================================================================================
  !=====================================================================================
  !=====================================================================================
  !=====================================================================================
  subroutine decay_turbulence(turb_decay_frac, old_turb_field, new_turb_field)
    implicit none
    ! Take a old complex Fourier field of turbulence, and find
    ! decay. Remove this decay from new turbulent field.
    real(kind=8), intent(in)       :: turb_decay_frac
    complex(kind=8), intent(in)    :: old_turb_field(1:NDIM,0:TGRID_X,0:TGRID_Y,0:TGRID_Z)
                                           ! Old field for reference
    complex(kind=8), intent(inout) :: new_turb_field(1:NDIM,0:TGRID_X,0:TGRID_Y,0:TGRID_Z)
                                           ! Complex field to subtract from
    ! Ornstein-Uhlenbeck process
    ! dF(k,t) = F_0(k) P(k) dW - F(k,t) dt / T
    ! where F(k,t) is the vector fourier amplitude,
    ! F_0 is the power spectrum distribution of amplitudes,
    ! P(k) is the projection (Helmholtz decomposition),
    ! dt is the timestep of integration,
    ! and dW is the Wiener process, described by a random variate
    ! from a normal distribution with a mean of zero and a
    ! variance of dt i.e. a standard deviation of sqrt(dt)

    ! In this we only calculate the second term - the exponential
    ! decay term.

    new_turb_field = new_turb_field - (old_turb_field * turb_decay_frac)

  end subroutine decay_turbulence
  !=====================================================================================
  !=====================================================================================
  !=====================================================================================
  !=====================================================================================
  subroutine FFT_1D(complex_field, real_field)
    implicit none
#ifdef TURB
#include "fftw3.f"
#endif
    ! Transform complex field into purely real field for 1D vector field

    complex(kind=8), intent(in)  :: complex_field(0:TGRID_X)
                                           ! Complex field to transform
    real(kind=8), intent(out)    :: real_field(0:TGRID_X)
                                           ! Result of transforms
#ifdef TURB
    integer (kind=ILP)   :: plan           ! FFTW plan
    complex(kind=8), allocatable :: fftfield(:) ! Memory for FFT

    ! Allocate storage for performing FFTs
    allocate(fftfield(0:TGRID_X))

    call dfftw_plan_dft_1d(plan, TURB_GS, fftfield, fftfield, FFTW_BACKWARD, FFTW_ESTIMATE)
    fftfield = complex_field(:)

    call dfftw_execute_dft(plan, fftfield, fftfield)
    real_field(:) = real(fftfield, kind=8) / (turb_gs_real)

    call dfftw_destroy_plan(plan)

    deallocate(fftfield)
#else
    real_field(:)=0
#endif
  end subroutine FFT_1D
  !=====================================================================================
  !=====================================================================================
  !=====================================================================================
  !=====================================================================================
  subroutine FFT_2D(complex_field, real_field)
    implicit none
#ifdef TURB
#include "fftw3.f"
#endif
    ! Transform complex field into purely real field for 2D vector field

    complex(kind=8), intent(in)  :: complex_field(1:2,0:TGRID_X,0:TGRID_Y)
                                           ! Complex field to transform
    real(kind=8), intent(out)    :: real_field(1:2,0:TGRID_X,0:TGRID_Y)
                                           ! Result of transforms
#ifdef TURB
    integer              :: d               ! Dimension counter
    integer (kind=ILP)   :: plan            ! FFTW plan
    complex(kind=8), allocatable :: fftfield(:,:) ! Memory for FFT

    ! Allocate storage for performing FFTs
    allocate(fftfield(0:TGRID_X,0:TGRID_Y))

    call dfftw_plan_dft_2d(plan, TURB_GS, TURB_GS, fftfield, fftfield, FFTW_BACKWARD, FFTW_ESTIMATE)

    do d=1,2
       fftfield = complex_field(d,:,:)
       call dfftw_execute_dft(plan, fftfield, fftfield)
       real_field(d,:,:) = real(fftfield, kind=8) / (turb_gs_real**2)
    end do

    call dfftw_destroy_plan(plan)

    deallocate(fftfield)
#else
    real_field=0
#endif
  end subroutine FFT_2D
  !=====================================================================================
  !=====================================================================================
  !=====================================================================================
  !=====================================================================================
  subroutine FFT_3D(complex_field, real_field)
    implicit none
#ifdef TURB
#include "fftw3.f"
#endif
    ! Transform complex field into purely real field for 3D vector field

    complex(kind=8), intent(in)  :: complex_field(1:3,0:TGRID_X,0:TGRID_Y,0:TGRID_Z)
                                           ! Complex field to transform
    real(kind=8), intent(out)    :: real_field(1:3,0:TGRID_X,0:TGRID_Y,0:TGRID_Z)
                                           ! Result of transforms
#ifdef TURB
    integer              :: d               ! Dimension counter
    integer (kind=ILP)   :: plan            ! FFTW plan
    complex(kind=8), allocatable :: fftfield(:,:,:) ! Memory for FFT

    ! Allocate storage for performing FFTs
    allocate(fftfield(0:TGRID_X,0:TGRID_Y,0:TGRID_Z))

    call dfftw_plan_dft_3d(plan, TURB_GS, TURB_GS, TURB_GS, fftfield, fftfield, FFTW_BACKWARD, FFTW_ESTIMATE)

    do d=1,3
       fftfield = complex_field(d,:,:,:)
       call dfftw_execute_dft(plan, fftfield, fftfield)
       real_field(d,:,:,:) = real(fftfield, kind=8) / (turb_gs_real**3)
    end do

    call dfftw_destroy_plan(plan)

    deallocate(fftfield)
#else
    real_field=0
#endif
  end subroutine FFT_3D
  !=====================================================================================
  !=====================================================================================
  !=====================================================================================
  !=====================================================================================
  subroutine proj_rms_norm(turb, P)
    implicit none
    ! rms of the projected random vector relative to the unprojected one.
    ! Each of the NDIM components has the same variance, one along k and
    ! NDIM-1 across it, so the ratio is exact and the same for every comp_frac.
    type(turb_t), intent(in)      :: turb
    real(kind=8), intent(out)     :: P           ! Normalization constant

#if NDIM>1
    P = sqrt((turb%comp_amp**2 + (NDIM - 1)*turb%sol_amp**2) / real(NDIM, kind=8))
#else
    P = 1.0d0
#endif

  end subroutine proj_rms_norm
  !=====================================================================================
  !=====================================================================================
  !=====================================================================================
  !=====================================================================================
  subroutine power_rms_norm(power_in, P)
    implicit none
    ! Calculate the normalizations for FFT of initial power spectrum
    real(kind=8), intent(in)  :: power_in(0:TGRID_X,0:TGRID_Y,0:TGRID_Z)
    real(kind=8), intent(out) :: P  ! Normalization constant

    complex(kind=8  )         :: complex_field(1:NDIM,0:TGRID_X,0:TGRID_Y, 0:TGRID_Z)
                                                 ! Complex field to transform
    real(kind=8)              :: real_field(1:NDIM,0:TGRID_X,0:TGRID_Y, 0:TGRID_Z)
                                                 ! Result of transforms
    integer                   :: d                ! Dimension counter

    do d=1,NDIM
       complex_field(d,:,:,:) = cmplx(power_in, kind=8)
    end do

#if NDIM==1
    call FFT_1D(complex_field, real_field)
#elif NDIM==2
    call FFT_2D(complex_field, real_field)
#else
    call FFT_3D(complex_field, real_field)
#endif

    P = sqrt(sum(real_field**2)/size(real_field))

  end subroutine power_rms_norm
  !=====================================================================================
  !=====================================================================================
  !=====================================================================================
  !=====================================================================================
  subroutine turb_force_calc(run, turb, ncache, x_cell, rho, aturb)
    implicit none
    type(run_t)               :: run
    type(turb_t)              :: turb
    integer, intent(in)       :: ncache
    real(kind=8), intent(in)  :: x_cell(1:nvector, 1:ndim) ! Positions
    real(kind=8), intent(in)  :: rho(1:nvector)            ! Densities
    real(kind=8), intent(out) :: aturb(1:nvector, 1:ndim)  ! Turbulent forcing

    integer                   :: nok                       ! no. of OK cells
    integer                   :: ok_cell(1:nvector)        ! 'ok' cells
    integer                   :: i                         ! cell counter
    real(kind=8)              :: r(1:ndim,1:nvector)       ! Position in turb grid
    real(kind=8)              :: dr1(1:ndim,1:nvector)     ! Position in cell
    real(kind=8)              :: dr2(1:ndim,1:nvector)     ! 1 - position in cell
    integer                   :: bmin(1:ndim,1:nvector)    ! 'top-left' corner of box
    integer                   :: bmax(1:ndim,1:nvector)    ! 'bottom-right' corner of box
    real(kind=8)              :: cube_vals(1:ndim,1:nvector,1:twotondim)
                                              ! Values of turbulent forcing
                                              ! from top left to bottom right
    real(kind=8)              :: interp(1:ndim,1:nvector,1:twotondim)
                                              ! Interpolation weights
    aturb = 0.0

    nok = 0
    do i = 1, ncache
       ! Position of particle in 'grid' space
       if (rho(i) < run%turb_min_rho) then
          continue ! Less than minimum density for adding turbulence
       else
          nok = nok + 1
          ok_cell(nok) = i
          r(:,nok) = x_cell(i,:)/run%boxlen*turb_gs_real
       end if
    end do

    ! Find ids of top left and bottom right corner of encompassing cube
    ! These can be the same if a particle is right on the boundary, but won't
    ! affect the result.
    bmin(:,1:nok) = floor(r(:,1:nok))
    bmax(:,1:nok) = ceiling(r(:,1:nok))

    ! Right-hand edge is equal to left-hand edge (periodic),
    ! which gives TURB_GS + 1 interpolation points
    where (bmin(:,1:nok)==TURB_GS) bmin(:,1:nok)=0
                                    ! this only happens for r==TURB_GS
    where (bmax(:,1:nok)==TURB_GS) bmax(:,1:nok)=0

    dr1(:,1:nok) = r(:,1:nok) - real(floor(r(:,1:nok)),kind=8)
    dr2(:,1:nok) = 1.0 - dr1(:,1:nok)

    ! Find cube values
#if NDIM==1
    do i=1,nok
       cube_vals(:,i,1) = turb%afield_now(:, bmin(1,i), 0, 0)
       cube_vals(:,i,2) = turb%afield_now(:, bmax(1,i), 0, 0)
    end do

    do i=1,nok
       interp(:,i,1) = dr2(1,i)
       interp(:,i,2) = dr1(1,i)
    end do
#elif NDIM==2
    do i=1,nok
       cube_vals(:,i,1) = turb%afield_now(:, bmin(1,i), bmin(2,i), 0)
       cube_vals(:,i,2) = turb%afield_now(:, bmax(1,i), bmin(2,i), 0)
       cube_vals(:,i,3) = turb%afield_now(:, bmin(1,i), bmax(2,i), 0)
       cube_vals(:,i,4) = turb%afield_now(:, bmax(1,i), bmax(2,i), 0)
    end do

    do i=1,nok
       interp(:,i,1) = dr2(1,i) * dr2(2,i)
       interp(:,i,2) = dr1(1,i) * dr2(2,i)
       interp(:,i,3) = dr2(1,i) * dr1(2,i)
       interp(:,i,4) = dr1(1,i) * dr1(2,i)
    end do
#else
    do i=1,nok
       cube_vals(:,i,1) = turb%afield_now(:, bmin(1,i), bmin(2,i), bmin(3,i))
       cube_vals(:,i,2) = turb%afield_now(:, bmax(1,i), bmin(2,i), bmin(3,i))
       cube_vals(:,i,3) = turb%afield_now(:, bmin(1,i), bmax(2,i), bmin(3,i))
       cube_vals(:,i,4) = turb%afield_now(:, bmax(1,i), bmax(2,i), bmin(3,i))
       cube_vals(:,i,5) = turb%afield_now(:, bmin(1,i), bmin(2,i), bmax(3,i))
       cube_vals(:,i,6) = turb%afield_now(:, bmax(1,i), bmin(2,i), bmax(3,i))
       cube_vals(:,i,7) = turb%afield_now(:, bmin(1,i), bmax(2,i), bmax(3,i))
       cube_vals(:,i,8) = turb%afield_now(:, bmax(1,i), bmax(2,i), bmax(3,i))
    end do

    ! Find interpolation values
    do i=1,nok
       interp(:,i,1) = dr2(1,i) * dr2(2,i) * dr2(3,i)
       interp(:,i,2) = dr1(1,i) * dr2(2,i) * dr2(3,i)
       interp(:,i,3) = dr2(1,i) * dr1(2,i) * dr2(3,i)
       interp(:,i,4) = dr1(1,i) * dr1(2,i) * dr2(3,i)
       interp(:,i,5) = dr2(1,i) * dr2(2,i) * dr1(3,i)
       interp(:,i,6) = dr1(1,i) * dr2(2,i) * dr1(3,i)
       interp(:,i,7) = dr2(1,i) * dr1(2,i) * dr1(3,i)
       interp(:,i,8) = dr1(1,i) * dr1(2,i) * dr1(3,i)
    end do
#endif

    ! Find interpolated value
    dr1(1:ndim,1:nok) = sum(interp(:,1:nok,:) * cube_vals(:,1:nok,:), dim=3)
    do i=1,nok
       aturb(ok_cell(i),:) = dr1(:,i)
    end do

  end subroutine turb_force_calc
  !=====================================================================================
  !=====================================================================================
  !=====================================================================================
  !=====================================================================================
  subroutine current_turb_rms(turb, rms_val)
    implicit none
    type(turb_t)              :: turb
    real (kind=8), intent(out) :: rms_val

    integer                    :: i, j, k
    real (kind=8)              :: asqd

    rms_val = 0.0
    do k = 0, TGRID_Z
       do j = 0, TGRID_Y
          do i = 0, TGRID_X
             asqd = sum(turb%afield_now(1:ndim, i, j, k)**2)
             rms_val = rms_val + asqd
          end do
       end do
    end do

    rms_val = sqrt(rms_val / (turb_gs_real**NDIM))

  end subroutine current_turb_rms
  !=====================================================================================
  !=====================================================================================
  !=====================================================================================
  !=====================================================================================
  ! PRNG of Marsaglia
  subroutine spin_up(s)
    integer(ILP), intent(inout)   :: s(4)
    integer                       :: i
    integer(ILP)                  :: r
    do i=1,10000
       call kiss64_core(s, r)
    end do
  end subroutine spin_up
  !=====================================================================================
  !=====================================================================================
  !=====================================================================================
  !=====================================================================================
  subroutine kiss64_core(s, r)
    integer(ILP), intent(inout) :: s(4)
    integer(ILP), intent(out)   :: r
    integer(ILP) :: t
    t = ishft(s(1), 58) + s(4)
    if (ishft(s(1), -63) .eq. ishft(t, -63)) then
       s(4) = ishft(s(1), -6) + ishft(s(1), -63)
    else
       s(4) = ishft(s(1), -6) + 1 - ishft(s(1) + t, -63)
    endif
    s(1) = t + s(1)
    s(2) = ieor(s(2), ishft(s(2),13))
    s(2) = ieor(s(2), ishft(s(2),-17))
    s(2) = ieor(s(2), ishft(s(2),43))
    s(3)= 6906969069_ILP * s(3) + 1234567_ILP
    r = s(1) + s(2) + s(3)
  end subroutine kiss64_core
  !=====================================================================================
  !=====================================================================================
  !=====================================================================================
  !=====================================================================================
  subroutine kiss64_double(N, s, out_array)
    integer, intent(in)           :: N
    integer(ILP), intent(inout)   :: s(4)
    real(kind=8), intent(out)         :: out_array(1:N)
    integer                       :: i
    integer(ILP)                  :: int_array(1:N)
    integer(ILP), parameter       :: randmax = 9223372036854775807_ILP
    double precision              :: dbl_array(1:N)

    do i=1,N
       call kiss64_core(s, int_array(i))
    end do

    int_array = iand(int_array, z'FFFFFFFFFFFFF')
    int_array = ieor(int_array, z'3FF0000000000000')
    out_array = transfer(int_array, out_array) - 1.0

  end subroutine kiss64_double
  !=====================================================================================
  !=====================================================================================
  !=====================================================================================
  !=====================================================================================
end module turb_commons
