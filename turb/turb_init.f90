module turb_init_module
  use ramses_commons, only: pst_t
  use amr_commons, only: run_t, global_t
  use turb_commons
contains
!################################################################
!################################################################
!################################################################
!################################################################
  recursive subroutine r_init_turb(pst)
    use mdl_module
    use mdl_parameters
#ifdef _CUDA
    use gpu_manager, only: gpu_turb_init_fields
#endif
    implicit none
    type(pst_t)::pst

    integer::rID

    if(pst%nLower>0)then
       rID = mdl_send_request(pst%s%mdl,MDL_INIT_TURB,pst%iUpper+1)
       call r_init_turb(pst%pLower)
       call mdl_get_reply(pst%s%mdl,rID,0)
    else
       call init_turb(pst%s%r,pst%s%g, pst%s%turb)
#ifdef _CUDA
       call gpu_turb_init_fields(pst)
#endif
    endif

  end subroutine r_init_turb
!################################################################
!################################################################
!################################################################
!################################################################
  subroutine init_turb(run, global, turb)
    implicit none
    type(run_t)  :: run
    type(global_t) :: global
    type(turb_t) :: turb
    !--------------------------------------------------
    ! Local variables
    !--------------------------------------------------
    integer      :: i, j, k         ! Loop variables
    integer      :: k_vec(1:3)      ! Wavevector
    integer      :: all_stat(1:3)   ! Allocation statuses

    integer      :: n_seed=4        ! Length of random seed, 4 for KISS64
    integer      :: clock           ! Integer clock time

    real(kind=8) :: power_norm      ! Normalization from power spectrum
    real(kind=8) :: proj_norm       ! Normalization from projection
    real(kind=8) :: OU_norm         ! Normalization for OU process

    real(kind=8) :: turb_last_tfrac ! Time fraction since last
    real(kind=8) :: turb_next_tfrac ! Time fraction until next

    real(kind=8)  :: cur_percent     ! Percentage currently tracking

    all_stat = 0

    ! Allocate grids
    allocate(turb%afield_last(1:NDIM,0:TGRID_X,0:TGRID_Y,0:TGRID_Z), stat=all_stat(1))
    allocate(turb%afield_next(1:NDIM,0:TGRID_X,0:TGRID_Y,0:TGRID_Z), stat=all_stat(2))
    allocate(turb%afield_now (1:NDIM,0:TGRID_X,0:TGRID_Y,0:TGRID_Z), stat=all_stat(3))

    if (any(all_stat /= 0)) stop 'Out of memory in init_turb!'

    ! Set grid spacing
    turb%turb_space = (/1.d0, 1.d0, 1.d0/) / turb_gs_real

    ! Set turbulence update time from autocorrelation time and number of substeps
    turb%turb_dt = run%turb_T / real(run%turb_Ndt,kind=8)

    if (run%nrestart == 0) then
       ! Set up random seed
       if (run%turb_seed == -1) then
          turb%kiss64_state = 1234
       else
          turb%kiss64_state = run%turb_seed
       end if
       call spin_up(turb%kiss64_state)
    end if

    ! Allocate grids
    allocate(turb%turb_last(1:NDIM,0:TGRID_X,0:TGRID_Y,0:TGRID_Z), stat=all_stat(1))
    allocate(turb%turb_next(1:NDIM,0:TGRID_X,0:TGRID_Y,0:TGRID_Z), stat=all_stat(2))
    allocate(turb%power_spec      (0:TGRID_X,0:TGRID_Y,0:TGRID_Z), stat=all_stat(3))

    if (any(all_stat /= 0)) stop 'Out of memory in init_turb!'

    ! Set decay fraction per timestep dt
    turb%turb_decay_frac = turb%turb_dt / run%turb_T

    ! comp_frac is the compressive (longitudinal) share of the forcing power.
    ! A mode has one longitudinal and NDIM-1 transverse degrees of freedom,
    ! so equal weights give comp_frac = 1/NDIM, the natural mixture
    ! (Federrath et al. 2010).
    if (run%comp_frac < 0.0d0 .or. run%comp_frac > 1.0d0) stop 'init_turb: comp_frac must lie in [0,1]'
    turb%comp_amp = sqrt(run%comp_frac)
#if NDIM>1
    turb%sol_amp = sqrt((1.0d0 - run%comp_frac) / real(NDIM - 1, kind=8))
#else
    turb%sol_amp = 0.0d0
#endif

    ! Set turbulent field time
    turb%turb_next_time = global%t

    ! Set initial power distribution (should be parameterised in some fashion)
    do k = 0, TGRID_Z
       if (k > TURB_GS / 2) then
          k_vec(3) = k - TURB_GS
       else
          k_vec(3) = k
       end if
       do j = 0, TGRID_Y
          if (j > TURB_GS / 2) then
             k_vec(2) = j - TURB_GS
          else
             k_vec(2) = j
          end if
          do i = 0, TGRID_X
             if (i > TURB_GS / 2) then
                k_vec(1) = i - TURB_GS
             else
                k_vec(1) = i
             end if
             call calc_power_spectrum(run, k_vec, turb%power_spec(i,j,k))
          end do
       end do
    end do

    ! Calculate turbulent normalization
    ! Power normalization comes from FFT of initial power spectrum
    call power_rms_norm(turb%power_spec, power_norm)
    ! Projection normalization of the random vectors
    call proj_rms_norm(turb, proj_norm)
    ! OU norm comes from standard deviation of OU process
    OU_norm = sqrt(run%turb_T / 2.0)
    ! Combination of all factored (reciprocal for easy multiplication)
    turb%turb_norm = 1.0 / (power_norm * proj_norm * OU_norm)

    ! Restart
    if (run%nrestart > 0) then
       ! Load turbulent fields from files and perform FFTs
       call read_turb_fields(run, turb)
#if NDIM==1
       call FFT_1D(turb%turb_last(1,:,0,0), turb%afield_last(1,:,0,0))
       call FFT_1D(turb%turb_next(1,:,0,0), turb%afield_next(1,:,0,0))
#elif NDIM==2
       call FFT_2D(turb%turb_last(:,:,:,0), turb%afield_last(:,:,:,0))
       call FFT_2D(turb%turb_next(:,:,:,0), turb%afield_next(:,:,:,0))
#else
       call FFT_3D(turb%turb_last, turb%afield_last)
       call FFT_3D(turb%turb_next, turb%afield_next)
#endif
       turb%afield_last = turb%afield_last * turb%turb_norm * run%turb_rms
       turb%afield_next = turb%afield_next * turb%turb_norm * run%turb_rms

    ! Not a restart
    else
       ! Set up initial field and perform FFT
       turb%turb_next = cmplx(0, 0, kind=8)
       call add_turbulence(turb, turb%turb_next, turb%turb_dt)
#if NDIM==1
       call FFT_1D(turb%turb_next(1,:,0,0), turb%afield_next(1,:,0,0))
#elif NDIM==2
       call FFT_2D(turb%turb_next(:,:,:,0), turb%afield_next(:,:,:,0))
#else
       call FFT_3D(turb%turb_next, turb%afield_next)
#endif
       turb%afield_next = turb%afield_next * turb%turb_norm * run%turb_rms

       ! Call turb_next_field to create second field
       call turb_next_field(run, turb)
    end if

    ! Set up afield_now
    turb_last_tfrac = real((global%t - turb%turb_last_time) / turb%turb_dt, kind=8)
    turb_next_tfrac = 1.0 - turb_last_tfrac

    turb%afield_now = turb_next_tfrac*turb%afield_last + turb_last_tfrac*turb%afield_next

  end subroutine init_turb
!################################################################
!################################################################
!################################################################
!################################################################
end module turb_init_module
