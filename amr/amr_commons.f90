module amr_commons
  use amr_parameters
  use hydro_parameters
  use rt_parameters
  use cr_parameters
  use oct_commons
  use hydro_commons
  use rt_commons
  use cr_commons
  use hash
  use domain_m

  type multipole_t
    real(kind=8),dimension(1:ndim+1)::q
  end type multipole_t

  type run_t

     ! Run control
     logical::cosmo   =.false.   ! Cosmology activated
     logical::pic     =.false.   ! Particle In Cell activated
     logical::poisson =.false.   ! Poisson solver activated
     logical::hydro   =.false.   ! Hydro activated
     logical::rt      =.false.   ! RT activated
     logical::cr      =.false.   ! Cosmic rays activated
     logical::part    =.false.   ! Dark matter particles activated
     logical::star    =.false.   ! Stars and star formation activated
     logical::sink    =.false.   ! Sinks and sink formation activated
     logical::tree    =.false.   ! Merger tree particles activated
     logical::trac    =.false.   ! Tracer particles activated
     logical::dust    =.false.   ! Dust particles activated
     logical::orphan  =.false.   ! Orphan particles activated
     logical::verbose =.false.   ! Write everything
     logical::debug   =.false.   ! Debug mode activated
     logical::clump_only =.false. ! Only clump finding
     integer::nrestart=0         ! New run or backup file number
     integer::ncontrol=1         ! Write control variables
     integer::nstepmax=1000000   ! Maximum number of time steps
     integer,dimension(1:MAXLEVEL)::nsubcycle=2 ! Subcycling at each level
     integer::nremap=10          ! Load balancing frequency (0: never)
     logical::static_mesh=.false.! Static mesh refinement mode activated
     logical::static_gas=.false. ! Turn off gas dynamics
     integer::geom    =1         ! 1: cartesian, 2: cylindrical, 3: spherical
     integer::overload=1         ! MPI domain overloading
     integer::nsuperoct=0        ! Number of superoct levels

     ! Output parameters
     integer::noutput=1          ! Total number of outputs
     integer::foutput=1000000    ! Frequency of outputs
     real(kind=8),dimension(1:MAXOUT)::aout=1.1 ! Output expansion factors
     real(kind=8),dimension(1:MAXOUT)::tout=0.0 ! Output times
     integer::output_mode=0      ! Output mode (for hires runs)
     logical::gadget_output=.false. ! Output in gadget format
     real(kind=8)::bkp_time_hrs=2   ! Backup file frequency in hours
     real(kind=8)::run_time_hrs=0   ! Estimated run time in hrs
     real(kind=8)::bkp_last_min=10  ! Backup file before the end of run in min
     integer::bkp_modulo=0       ! Use modulo for backup file count
     integer::nfile=1            ! Number of file used per snapshot. Use -1 for nfile=ncpu
     logical::output_part=.true. ! Output particle data in dumps
     logical::output_grav=.true. ! Output gravity data in dumps
     logical::output_hydro=.true.! Output hydro data in dumps
     logical::output_amr=.true.  ! Output AMR data in dumps

     ! Trajectory output parameters (per-step particle traces)
     integer::ntrajectories=0
     integer,dimension(1:MAXOUT)::trajectories=0

     ! Mesh parameters
     integer::levelmin=1         ! Full refinement up to levelmin
     integer::nlevelmax=1        ! Maximum number of level
     integer::ngridmax=0         ! Maximum number of grids
     integer::ncachemax=10000    ! Maximum number of cache lines
     integer::npartmax=0         ! Maximum number of DM particles
     integer::nstarmax=0         ! Maximum number of star particles
     integer::nsinkmax=0         ! Maximum number of sink particles
     integer::ntreemax=0         ! Maximum number of tree particles
     integer::ntracmax=0         ! Maximum number of tracer particles
     integer::ndustmax=0         ! Maximum number of dust particles
     integer::ntrac_per_cell=1   ! Number of tracer particles per cell in ICs
     integer::ndust_per_cell=1   ! Number of dust particles per cell in ICs
     logical::part_subcell_positions=.true.  ! Use subcell positions for particles
     character(LEN=32)::tracer_kick_pdf='piecewise_skew_uniform' ! Tracer kick PDF type
     real(kind=8)::dust_to_gas_mass_ratio=0.0d0   ! Dust to gas mass ratio
     real(kind=8)::grain_size_parameter=0.0d0   ! Dust particle size parameter (dimensionless)
     real(kind=8)::grain_charge_parameter=0.0d0  ! Dust particle charge-to-mass ratio (dimensionless)
     real(kind=8)::dust_gyro_factor=0.1d0  ! Fraction of gyro period allowed per step
     logical :: analytic_dust_force = .false. ! If true use analytic drag+Lorentz
     integer,dimension(1:MAXLEVEL)::nexpand=1 ! Number of mesh expansion
     real(kind=8)::boxlen=1.0        ! Cell size at level 0 (total box size)

     ! Poisson solver parameters
     logical :: gravity_test=.false.  ! Use file rho_ana.f90 to test the Poisson solvers
     real(kind=8)::epsilon=1.0D-4     ! Convergence criterion for Poisson solvers
     integer :: nvcycle = -1          ! Desired number of V-cycles
     real(kind=8),dimension(1:10)::gravity_params=0.0 ! Gravity parameters
     integer :: gravity_type=0     ! Type of force computation
     integer :: cic_levelmax=0     ! Maximum level for CIC dark matter interpolation
     integer :: cg_levelmin=999    ! Min level for CG solver
     logical :: fast_solver = .false. ! Fast solver with MPI pre-fetch (memory intensive)
     integer :: part_mass_deposition_scheme=1     ! part mass deposition schemes (CIC 1, TSC 2, PCS 3)
     integer :: part_dep_algo=2                   ! part GPU CIC deposition algorithm
     ! 1: large 27-offset
     ! 2: medium shifted 8-offset (fastest, default)
     ! 3: small shifted 8-offset prefix-sum)
     integer :: part_force_interpolation_scheme=1 ! part force interpolation schemes (CIC 1, TSC 2, PCS 3)
     integer :: star_mass_deposition_scheme=1     ! star mass deposition schemes
     integer :: star_dep_algo=2                   ! star GPU CIC deposition algorithm
     integer :: star_force_interpolation_scheme=1 ! star force interpolation schemes
     integer :: sink_mass_deposition_scheme=1     ! sink mass deposition schemes
     integer :: sink_dep_algo=2                   ! sink GPU CIC deposition algorithm
     integer :: sink_force_interpolation_scheme=1 ! sink force interpolation schemes
     integer :: tree_mass_deposition_scheme=1     ! tree mass deposition schemes
     integer :: tree_force_interpolation_scheme=1 ! tree force interpolation schemes
     integer :: trac_interpolation_scheme=1 ! tracer force interpolation schemes
     integer :: dust_mass_deposition_scheme=1 ! dust mass deposition schemes
     integer :: dust_force_interpolation_scheme=1 ! dust force interpolation schemes

     ! Movie parameters
     integer::levelmax_frame=0
     integer::nw_frame=512 ! width of frame in pixels
     integer::nh_frame=512 ! height of frame in pixels
     integer::ivar_frame=1
     real(kind=8),dimension(1:20)::xcentre_frame=0d0
     real(kind=8),dimension(1:20)::ycentre_frame=0d0
     real(kind=8),dimension(1:20)::zcentre_frame=0d0
     real(kind=8),dimension(1:10)::deltax_frame=0d0
     real(kind=8),dimension(1:10)::deltay_frame=0d0
     real(kind=8),dimension(1:10)::deltaz_frame=0d0
     logical::movie=.false.
     logical::zoom_only=.false.
     integer::imovout=0    ! Increment for output times
     integer::imov=1       ! Initialize
     real(kind=8)::tstartmov=0d0
     real(kind=8)::astartmov=0d0
     real(kind=8)::tendmov=0.
     real(kind=8)::aendmov=0.
     character(LEN=5)::proj_axis='z' ! x->x, y->y, projection along z
     integer,dimension(0:NVAR+2+nrtgrp+ncrgrp)::movie_vars=0
     character(len=5),dimension(0:NVAR+2+nrtgrp+ncrgrp)::movie_vars_txt=''
     ! Movie camera and rendering options (per-projection, NMOV=5)
     real(kind=8),dimension(1:5)::theta_camera=0d0
     real(kind=8),dimension(1:5)::phi_camera=0d0
     real(kind=8),dimension(1:5)::dtheta_camera=0d0
     real(kind=8),dimension(1:5)::dphi_camera=0d0
     real(kind=8),dimension(1:5)::tstart_theta_camera=0d0
     real(kind=8),dimension(1:5)::tstart_phi_camera=0d0
     real(kind=8),dimension(1:5)::tend_theta_camera=0d0
     real(kind=8),dimension(1:5)::tend_phi_camera=0d0
     real(kind=8),dimension(1:5)::focal_camera=0d0
     real(kind=8),dimension(1:5)::dist_camera=0d0
     real(kind=8),dimension(1:5)::ddist_camera=0d0
     real(kind=8),dimension(1:5)::smooth_frame=1d0
     real(kind=8),dimension(1:5)::varmin_frame=-1d60
     real(kind=8),dimension(1:5)::varmax_frame=1d60
     logical,dimension(1:5)::perspective_camera=.false.
     logical,dimension(1:5)::zoom_only_frame=.false.
     character(LEN=6),dimension(1:5)::shader_frame='square'
     character(LEN=10),dimension(1:5)::method_frame='mean_mass'

     ! Hydro solver parameters
     real(kind=8)::gamma=1.4d0
     real(kind=8)::courant_factor=0.5d0
     real(kind=8)::smallc=1.d-10
     real(kind=8)::smallr=1.d-10
     integer ::niter_riemann=10
     integer ::slope_type=1
     integer ::slope_mag_type=1
     real(kind=8)::difmag=0.0d0
     real(kind=8)::etamag=0.0d0
     real(kind=8),dimension(1:nener+1)::gamma_rad=1.33333333334d0
     logical ::induction=.false.
     logical ::entropy=.false.
     logical ::sgs_turb=.false.
     logical ::equilibrium_sgs=.false.       ! Equilibrium SGS turbulence
     real(kind=8)::smagorinsky_lilly_constant=1.0d0 ! Smagorinsky-Lilly constant for SGS turbulence
     real(kind=8)::dual_energy=-1
     real(kind=8)::T2_fix=0d0
     character(LEN=10)::scheme='muscl'
     integer::riemann=0
     integer::riemann2d=0
     real(kind=8)::switch_llf_dmin=-1
     real(kind=8)::switch_llf_pmin=-1
     real(kind=8),dimension(1:3)::constant_gravity
     integer::inener,ientropy,imetal,iturb,ichem,iecr

     ! Physics parameters
     real(kind=8)::units_density=1.0 ! [g/cm^3]
     real(kind=8)::units_time=1.0    ! [seconds]
     real(kind=8)::units_length=1.0  ! [cm]
     real(kind=8)::units_np=1.0      ! [#photon/cm^3]

     ! Cosmological parameters (others are read from file)
     real(kind=8)::boxlen_ini     ! Box size in h-1 Mpc
     real(kind=8)::omega_b=0.0D0  ! Omega Baryon
     real(kind=8)::omega_m=1.0D0  ! Omega Matter
     real(kind=8)::omega_l=0.0D0  ! Omega Lambda
     real(kind=8)::h0=1.0D0       ! Hubble constant in km/s/Mpc
     real(kind=8)::aexp_ini=10.   ! Starting expansion factor

     ! Refinement parameters for each level
     real(kind=8),dimension(1:MAXLEVEL)::m_refine = -1.0 ! Lagrangian threshold
     real(kind=8),dimension(1:MAXLEVEL)::r_refine = -1.0 ! Radius of refinement region
     real(kind=8),dimension(1:MAXLEVEL)::x_refine = 0.0 ! Center of refinement region
     real(kind=8),dimension(1:MAXLEVEL)::y_refine = 0.0 ! Center of refinement region
     real(kind=8),dimension(1:MAXLEVEL)::z_refine = 0.0 ! Center of refinement region
     real(kind=8),dimension(1:MAXLEVEL)::exp_refine = 2.0 ! Exponent for distance
     real(kind=8),dimension(1:MAXLEVEL)::a_refine = 1.0 ! Ellipticity (Y/X)
     real(kind=8),dimension(1:MAXLEVEL)::b_refine = 1.0 ! Ellipticity (Z/X)
     real(kind=8),dimension(1:MAXLEVEL)::jeans_refine=-1.0 ! Number of cells per Jeans length
     real(kind=8)::var_cut_refine=-1.0 ! Threshold for variable-based refinement
     real(kind=8)::mass_cut_refine=-1.0 ! Mass threshold for particle-based refinement
     integer::ivar_refine=-1 ! Variable index for refinement
     real(kind=8)::aexp_lock_refine=-1.0
     logical::pic_lock_refine=.false.

     ! Refinement parameters for hydro
     integer::interpol_var=0   ! Interpolated variables
     integer::interpol_type=1  ! Interpolation scheme
     real(kind=8)::err_grad_d=-1.0  ! Density gradient
     real(kind=8)::err_grad_u=-1.0  ! Velocity gradient
     real(kind=8)::err_grad_p=-1.0  ! Pressure gradient
     real(kind=8)::err_grad_xHI=-1  ! xHI gradient
     real(kind=8)::err_grad_xHII=-1 ! xHII gradient
     real(kind=8)::floor_d=1.d-10   ! Density floor
     real(kind=8)::floor_u=1.d-10   ! Velocity floor
     real(kind=8)::floor_p=1.d-10   ! Pressure floor
     real(kind=8)::floor_xHI=1d-10  ! xHI floor
     real(kind=8)::floor_xHII=1d-10 ! xHII floor
     real(kind=8)::mass_sph=0.0D0   ! mass_sph
#ifdef MHD
     real(kind=8)::err_grad_b2=-1.0
     real(kind=8)::err_grad_A=-1.0
     real(kind=8)::err_grad_B=-1.0
     real(kind=8)::err_grad_C=-1.0
     real(kind=8)::floor_b2=1.d-10
     real(kind=8)::floor_A=1.d-10
     real(kind=8)::floor_B=1.d-10
     real(kind=8)::floor_C=1.d-10
#endif
#if NENER>0
     real(kind=8),dimension(1:NENER)::err_grad_prad=-1.0
     real(kind=8),dimension(1:NENER)::floor_prad=1d-10
#endif
#if NVAR>5+NENER
     real(kind=8),dimension(1:NVAR-5-NENER)::err_grad_var=-1.0
#endif

     ! Initial condition regions parameters
     integer::nregion=0
     character(LEN=10),dimension(1:MAXREGION)::region_type='square'
     real(kind=8),dimension(1:MAXREGION)::x_center=0.
     real(kind=8),dimension(1:MAXREGION)::y_center=0.
     real(kind=8),dimension(1:MAXREGION)::z_center=0.
     real(kind=8),dimension(1:MAXREGION)::length_x=1.d10
     real(kind=8),dimension(1:MAXREGION)::length_y=1.d10
     real(kind=8),dimension(1:MAXREGION)::length_z=1.d10
     real(kind=8),dimension(1:MAXREGION)::exp_region=2.0

     ! Initial conditions hydro variables
     real(kind=8),dimension(1:MAXREGION)::d_region=0.
     real(kind=8),dimension(1:MAXREGION)::u_region=0.
     real(kind=8),dimension(1:MAXREGION)::v_region=0.
     real(kind=8),dimension(1:MAXREGION)::w_region=0.
     real(kind=8),dimension(1:MAXREGION)::p_region=0.
#if NENER>0
     real(kind=8),dimension(1:MAXREGION,1:NENER)::prad_region=0.0
#endif
#if NVAR>5+NENER
     real(kind=8),dimension(1:MAXREGION,1:NVAR-5-NENER)::var_region=0.0
#endif
#ifdef MHD
     real(kind=8),dimension(1:MAXREGION)::B_region=0.
     real(kind=8),dimension(1:MAXREGION)::C_region=0.
     real(kind=8)::A_ave=0.,B_ave=0.,C_ave=0.
#endif
     ! Initial condition files for each level
     character(LEN=20)::filetype='ascii'
     logical::multiple=.false.
     character(LEN=80),dimension(1:MAXLEVEL)::initfile=' '
     real(kind=8)::ic_scale_l=1.0d0
     real(kind=8)::ic_scale_v=1.0d0
     real(kind=8)::ic_scale_m=1.0d0

     ! Boundary conditions parameters
     logical,dimension(1:3)::periodic=.true.
     integer::nbound=0
     logical::no_inflow=.false.
     integer::bound_levelmin=1   ! AMR level to define boundary geometry
     real(kind=8),dimension(1:3)::box_size=0.0  ! Box length of active domain along each direction
     integer::box_xmin,box_xmax  ! Min and max Cartesian keys at levelmin
     integer::box_ymin,box_ymax  ! Min and max Cartesian keys at levelmin
     integer::box_zmin,box_zmax  ! Min and max Cartesian keys at levelmin
     integer,dimension(1:MAXBOUND)::bound_type=0
     integer,dimension(1:MAXBOUND)::bound_dir=0
     integer,dimension(1:MAXBOUND)::bound_shift=0
     integer,dimension(1:MAXBOUND)::bound_xmin=0
     integer,dimension(1:MAXBOUND)::bound_xmax=0
     integer,dimension(1:MAXBOUND)::bound_ymin=0
     integer,dimension(1:MAXBOUND)::bound_ymax=0
     integer,dimension(1:MAXBOUND)::bound_zmin=0
     integer,dimension(1:MAXBOUND)::bound_zmax=0
     real(kind=8),dimension(1:MAXBOUND)::d_bound=0
     real(kind=8),dimension(1:MAXBOUND)::u_bound=0
     real(kind=8),dimension(1:MAXBOUND)::v_bound=0
     real(kind=8),dimension(1:MAXBOUND)::w_bound=0
     real(kind=8),dimension(1:MAXBOUND)::p_bound=0
#if NENER>0
     real(kind=8),dimension(1:MAXBOUND,1:NENER)::prad_bound=0
#endif
#if NVAR>5+NENER
     real(kind=8),dimension(1:MAXBOUND,1:NVAR-5-NENER)::var_bound=0
#endif
#ifdef MHD
     real(kind=8),dimension(1:MAXBOUND)::B_bound=0
     real(kind=8),dimension(1:MAXBOUND)::C_bound=0
#endif
     ! Cooling parameters
     logical::neq_chem=.false.
     logical::cooling=.false.
     logical::cooling_ism=.false.
     logical::metal=.false.
     logical::isothermal=.false.
     logical::haardt_madau=.false.
     logical::self_shielding=.false.
     real(kind=8)::J21=0d0,a_spec=1d0,z_ave=0d0,z_reion=8.5d0,cooling_uvb_delta=0.05d0
     integer::eos_type=1 ! 1=isothermal, 2=polytrope, 3=isothermal+polytrope
     real(kind=8)::eos_nH=1d50,eos_index=1d0,eos_T2=10d0
     real(kind=8)::T2max
     real(kind=8)::mu_mol       = 1.2195d0           ! Mean molecular weight (std ISM value)
     real(kind=8)::X_H          = 0.7600d0           !                Hydrogen mass fraction
     real(kind=8)::Y_He         = 0.2400d0           !                  Helium mass fraction
     logical::is_init_xion=.false.          ! Initialize ionization from T profile (neq only)
     logical::isHe=.true.                             !      He ionization fractions tracked?
     logical::isH2=.false.                            !                           H2 tracked?
     integer::iIons, ixHI, ixHII, ixHeII, ixHeIII     !       Indexes of ionization fractions
#ifdef DO_RTZ
     integer::element_first_idx(1:27)       !  Start idx of elements in uold
     integer::molecules_first_idx(1:3)     !  Start idx of molecules in uold
     real(kind=8),dimension(1:27,1:27)::ionEvs        !                   Ionization energies
#else
     real(kind=8),dimension(nIon)::ionEvs             !                   Ionization energies
#endif
     real(kind=8)::neq_Tconst=-1           ! If positive use this value for all T-dependent rates
     logical::neq_isTconst=.false.                    !             Constant rates activated?
     logical::upload_equilibrium_x=.false.          ! Enforce equilibrium xion when uploading

     ! RTZ cooling parameters
     logical::rtz_cooling=.false.
     integer::rtz_equilibrium_test=-1
     logical::rtz_include_collisional_ionization=.true.
     logical::rtz_include_photoionization=.true.
     logical::rtz_include_cosmic_ray_ionization=.true.
     logical::rtz_include_charge_exchange=.true.
     logical::rtz_include_dust_recombination=.true.
     logical::rtz_include_HM12_UVB=.true.
     logical::isH2_rtz=.false.
     real(kind=8)::rtz_UV_background_G0=0.d0
     real(kind=8)::rtz_primary_cosmic_ray_ionization_rate=0.d0
     real(kind=8)::rtz_max_cool_timestep=1.d11
     integer::rtz_eqm_min_its=3000

     ! Star formation parameters
     integer::sf_model=1
     real(kind=8)::T2_star=2e4
     real(kind=8)::n_star=0.1
     real(kind=8)::eps_star=0.01
     integer(kind=8),dimension(1:6)::seed
     real(kind=8)::m_star=1

     ! Supernovae feedback parameters
     real(kind=8)::M_SNII=10.
     real(kind=8)::E_SNII=1d51
     real(kind=8)::t_SNII=20.
     real(kind=8)::eta_SNII=0.1
     real(kind=8)::yield_SNII=0.1
     logical::thermal_feedback=.false.
     logical::mechanical_feedback=.false.

     ! Clump finder parameters
     logical::clump_finder=.false.
     logical::clump_info=.false.
     logical::output_clump=.false.
     logical::output_peak_grid=.false.
     logical::output_peak_part=.false.
     logical::output_peak_star=.false.
     logical::output_peak_sink=.false.
     logical::output_peak_tree=.false.
     logical::output_peak_trac=.false.
     logical::output_peak_dust=.false.
     integer::rho_type_clump=1
     integer::nsteps_per_tree=1
     real(kind=8)::relevance_threshold=2
     real(kind=8)::density_threshold=-1
     real(kind=8)::saddle_threshold=-1
     real(kind=8)::mass_threshold=0
     real(kind=8)::purity_threshold=-1
     real(kind=8)::fraction_threshold=2d0

     ! Lightcone parameters
     logical::lightcone = .false.   ! Lightcone activated
     real(kind=8)::cone_z_min = 0.0 ! Minimum redshift
     real(kind=8)::cone_z_max = 1.0 ! Maximum redshift
     real(kind=8)::cone_opening_angle_y = 0.0 ! Opening angle in y direction in degrees
     real(kind=8)::cone_opening_angle_z = 0.0 ! Opening angle in z direction in degrees
     real(kind=8)::cone_theta = 0.0 ! Rotation of the cone's x-axis around the box's y-axis in degrees
     real(kind=8)::cone_phi = 0.0 ! Rotation of the cone's x-axis around the box's z-axis in degrees
     real(kind=8),dimension(1:3) :: cone_observer = (/0.0, 0.0, 0.0/) ! Observer position in code units

     ! Sink parameters
     integer::rho_type_sink=1
     logical::sink_descent=.false.
     real(kind=8)::fudge_descent=0.5d0
     real(kind=8)::sink_relevance_threshold=2
     real(kind=8)::sink_density_threshold=-1
     real(kind=8)::sink_saddle_threshold=-1
     real(kind=8)::sink_mass_threshold=0
     real(kind=8)::sink_purity_threshold=-1
     real(kind=8)::sink_fraction_threshold=2d0
     real(kind=8)::sink_radius=-1
     real(kind=8)::sink_delta_tout=0
     real(kind=8)::sink_nstar_frac=-1
     real(kind=8)::sink_mseed
     logical::sink_form=.false.
     logical::sink_merge=.false.
     logical::sink_refine=.false.
     logical::sink_dump=.false.
     logical::static_sink=.false.
     logical::fix_sink_mass=.false.
     logical::drag_sink=.false.

     ! Black hole parameters
     integer::accretion_type = 0 ! 0: None, 1: Bondi
     real(kind=8)::acc_sink_boost = 1.0d0 ! Boost for bondi accretion
     logical::bondi_use_vrel = .true. ! Whether to use the relative sink velocity for BHL accretion
     real(kind=8)::eddington_cap = -1 ! Factor of Eddington rate to cap accretion at
     integer::sink_b_spline_order = 4 ! Order of B-spline interpolation used for sink accretion and dynamics
     logical::verbose_sink = .false. ! Whether to print verbose statements for sink particles
     logical::bondi_use_gas_mass = .true. ! Whether to include the local gas mass in the Bondi calculation
     logical::use_local_bondi_rate = .false. ! Switch to average after (true) or before (false) computing the Bondi rate
     logical::use_rho_inf = .true. ! Whether to use bondi_alpha(x) to extrapolate density at infinity from Bondi solution
     real(kind=8)::t_start_black_hole = -1 ! Time after which to start using sink particle/black hole routines
     logical::use_bondi_lambda = .true.
     logical::mass_weighting = .true.
     logical::momentum_conserving = .false.
     real(kind=8)::eddington_floor = -1 ! Accretion rate floor below which nothing happens
     real(kind=8)::manual_accretion_rate = -1 ! Manual accretion rate (fraction of Eddington)

     ! AGN Feedback parameters
     logical::agn = .false. ! Whether to activate AGN feedback around black hole/sink particles
     integer::agn_feedback_radius = 4 ! Radius (in dx_min) of feedback region (should be geq sink_b_spline_order/2)
     integer::agn_weighting_scheme = 1 ! Which AGN weighting scheme (psy_function) to use
     real(kind=8)::epsilon_rad = 0.1d0 ! Radiative efficiency
     real(kind=8)::epsilon_radio = 1.0d0 ! Efficiency of momentum feedback for jet
     real(kind=8)::epsilon_quasar = 0.15d0 ! Efficiency of thermal feedback for quasar
     real(kind=8)::momentum_boost = 10.0d0 ! Momentum boost in units of L/c for the jet
     real(kind=8)::agn_fbk_mode_switch_threshold = 0.01d0 ! Threshold accretion rate to switch from jet to quasar mode
     real(kind=8)::agn_jet_opening_angle = 60.0d0 !  Outflow cone opening angle; in deg
     logical::agn_use_mass_weighting = .false. ! Whether to use a mass-weighted feedback scheme

     ! Gadget initial conditions parameters
     character(len=flen)::ic_file, ic_format
     integer,dimension(1:6)::ic_skip_type=-1
     character(len=4)::ic_head_name = 'HEAD'
     character(len=4)::ic_pos_name = 'POS '
     character(len=4)::ic_vel_name = 'VEL '
     character(len=4)::ic_id_name = 'ID  '
     character(len=4)::ic_mass_name = 'MASS'
     character(len=4)::ic_u_name = 'U   '
     character(len=4)::ic_metal_name = 'Z   '
     character(len=4)::ic_age_name = 'AGE '
     real(kind=8)::gadget_scale_l = 3.085677581282D21 ! Gadget units in cgs
     real(kind=8)::gadget_scale_v = 1.0D5
     real(kind=8)::gadget_scale_m = 1.9891D43
     real(kind=8)::gadget_scale_t = 1.0D6*365*24*3600
     real(kind=8)::IG_rho = 1.0D-6
     real(kind=8)::IG_T2 = 1.0D7
     real(kind=8)::IG_metal = 0.01

     ! RT parameters. Some parameters are not (yet) used
     logical::rt_advect=.false.            ! Advection of photons?                           !
     logical::rt_smooth=.true.             ! Smooth the discrete RT update of op. splitting  !
     logical::rt_star=.false.              ! Activate radiation from star particles          !
     logical::rt_sink=.false.              ! Activate radiation from sink particles          !
     real(kind=8)::rt_esc_frac=1d0             ! Photon escape fraction from stellar particles   !
     logical::rt_is_outflow_bound=.false.  ! Make all boundaries=outflow for RT              !
     real(kind=8)::rt_courant_factor=0.8d0     ! Courant factor for RT timesteps                 !
     real(kind=8)::rt_err_grad_cn(nrtgrp)=-1   ! Photon flux gradient for refinement             !
     real(kind=8)::rt_floor_cn(nrtgrp)=1d-10   ! Photon flux floor for refinement                !
     real(kind=8)::rt_refine_aexp=-1           ! Expansion factor for rt refinements             !
     real(kind=8),dimension(1:MAXLEVEL)::rt_c_fraction=1d0 ! Reduced light speed on each level   !
     integer::rt_nsubcycle=1               ! Maximum number of RT-steps during one hydro/    !
                                           ! gravity/etc timestep                            !
     logical::rt_otsa=.true.               ! Use on-the-spot approximation                   !
     logical::rt_isIR=.false.              ! Use IR photon groups                            !
     logical::is_kIR_T=.false.             ! Kappa IT depends on T_rad                       !
     logical::rt_T_rad=.false.             ! Use radiation temperature for everything        !
     logical::rt_isoPress=.false.          ! Use cE, not F, for rad. pressure                !
     real(kind=8)::rt_pressBoost=1d0           ! Boost on RT pressure                            !
     logical::is_mu_H2=.false.
     real(kind=8)::Tmu_dissoc=1d3              ! Dissociation temperature [K]                    !
     integer::iPEH_group=-1                ! Radiation group used for photo-electric heating !
     logical::cosmic_rays=.false.          ! Include cosmic ray ionisation                   !
     character(LEN=128)::sed_dir=''        ! Dir containing stellar energy distributions     !
     !character(LEN=128)::uv_file=''       ! File containing UV background                   !
     ! SED statistics: Radiation emitted, total, last coarse step [#photons/10^50]-----------
     logical::rt_emission_stats=.false.    ! Print info about stellar emission in log        !

     ! Initial condition RT regions parameters----------------------------------------------
     integer                           ::rt_nregion=0
     character(LEN=10),dimension(1:MAXREGION)::rt_region_type='square'
     real(kind=8),dimension(1:MAXREGION)   ::rt_reg_x_center=0.
     real(kind=8),dimension(1:MAXREGION)   ::rt_reg_y_center=0.
     real(kind=8),dimension(1:MAXREGION)   ::rt_reg_z_center=0.
     real(kind=8),dimension(1:MAXREGION)   ::rt_reg_length_x=1.E10
     real(kind=8),dimension(1:MAXREGION)   ::rt_reg_length_y=1.E10
     real(kind=8),dimension(1:MAXREGION)   ::rt_reg_length_z=1.E10
     real(kind=8),dimension(1:MAXREGION)   ::rt_exp_region=2.0
     integer,dimension(1:MAXREGION)    ::rt_reg_group=1
     real(kind=8),dimension(1:MAXREGION)   ::rt_n_region=0.                     ! Photon density
     real(kind=8),dimension(1:MAXREGION)   ::rt_u_region=0.                     ! Photon flux
     real(kind=8),dimension(1:MAXREGION)   ::rt_v_region=0.                     ! Photon flux
     real(kind=8),dimension(1:MAXREGION)   ::rt_w_region=0.                     ! Photon flux

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
     integer, dimension(1:MAXREGION)   ::rt_src_trace_group=1
     real(kind=8),dimension(1:MAXREGION)   ::rt_n_source=0.                     ! Photon density
     real(kind=8),dimension(1:MAXREGION)   ::rt_u_source=0.                     ! Photon flux
     real(kind=8),dimension(1:MAXREGION)   ::rt_v_source=0.                     ! Photon flux
     real(kind=8),dimension(1:MAXREGION)   ::rt_w_source=0.                     ! Photon flux

     ! RT boundary condition parameters-------------------------------------------------------
     real(kind=8),dimension(1:MAXBOUND,1:nrtgrp)::rt_n_bound=0.0d0
     real(kind=8),dimension(1:MAXBOUND,1:nrtgrp)::rt_u_bound=0.0d0
     real(kind=8),dimension(1:MAXBOUND,1:nrtgrp)::rt_v_bound=0.0d0
     real(kind=8),dimension(1:MAXBOUND,1:nrtgrp)::rt_w_bound=0.0d0

     ! RT groups parameters-------------------------------------------------------------------
     integer::sedprops_update=-1           ! Update sedprops from stellar populations
     ! negative: never update, 0:update on init, pos x: update every x coarse steps
     ! logical::SED_isEgy=.false. ! Integrate energy out of SEDs rather than photon count
     ! Group props: avg and energy weigthed photoionization c-section (cm2), avg. energy (ev)
     ! Indexes nrtgrp, nIon stand for photon group vs species (e.g. 1=H, 2=He)
     ! Note that this is not true in the case of RTZ
#ifdef DO_RTZ
     real(kind=8),dimension(nrtgrp,1:27,1:27)::group_csn=0, group_cse=0         ! Cross sections (cm2)
#else
     real(kind=8),dimension(nrtgrp,nIon)::group_csn=0, group_cse=0         ! Cross sections (cm2)
#endif
     real(kind=8),dimension(nrtgrp)::group_egy=0                        !  Avg photon energy (ev)
     real(kind=8),dimension(nrtgrp)::group_L0=13.60                     ! Wavelength lower limits
     real(kind=8),dimension(nrtgrp)::group_L1=0                         ! Wavelength upper limits
     real(kind=8),dimension(nrtgrp)::kappaAbs=0                         ! Dust absorption opacity
     real(kind=8),dimension(nrtgrp)::kappaSc=0                          ! Dust scattering opacity
     real(kind=8),dimension(nrtgrp)::isLW=0d0                          ! Use to find the LW group
     real(kind=8),dimension(nrtgrp)::ssh2=1d0                      ! Self-shielding factor for H2
     ! HK note --> OTSA required for RTZ
     integer,dimension(nIon)::spec2group=0                 ! Ion -> group # in recombinations

     ! CR boundary condition parameters-------------------------------------------------------
     real(kind=8),dimension(1:MAXBOUND,1:ncrgrp)::cr_fx_bound=0.0d0
     real(kind=8),dimension(1:MAXBOUND,1:ncrgrp)::cr_fy_bound=0.0d0
     real(kind=8),dimension(1:MAXBOUND,1:ncrgrp)::cr_fz_bound=0.0d0

     ! CR namelist parameters.
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
     ! CR group parameters---------------------------------------------------------------------
     real(kind=8),dimension(1:ncrgrp)::cr_d=1.0d29 !Classical value cm^2/s (e.g., Jockipii 1999)
     real(kind=8),dimension(1:ncrgrp)::cr_d_perp_factors=1d-6 ! perp diffusion CR suppression  !
     real(kind=8),dimension(1:ncrgrp)::fecr=0d0     ! SN fraction of CR energy
     real(kind=8)::cr_v_alfven=0.0                  ! For idealised tests
     ! Initial condition CR regions parameters-------------------------------------------------
     integer                           ::cr_nregion=0
     character(LEN=10),dimension(1:MAXREGION)::cr_region_type='square'
     real(kind=8),dimension(1:MAXREGION)   ::cr_reg_x_center=0.
     real(kind=8),dimension(1:MAXREGION)   ::cr_reg_y_center=0.
     real(kind=8),dimension(1:MAXREGION)   ::cr_reg_z_center=0.
     real(kind=8),dimension(1:MAXREGION)   ::cr_reg_length_x=1.E10
     real(kind=8),dimension(1:MAXREGION)   ::cr_reg_length_y=1.E10
     real(kind=8),dimension(1:MAXREGION)   ::cr_reg_length_z=1.E10
     real(kind=8),dimension(1:MAXREGION)   ::cr_exp_region=2.0
     integer,dimension(1:MAXREGION)        ::cr_reg_group=1
     real(kind=8),dimension(1:MAXREGION)   ::cr_fx_region=0.                    ! CR flux
     real(kind=8),dimension(1:MAXREGION)   ::cr_fy_region=0.                    ! CR flux
     real(kind=8),dimension(1:MAXREGION)   ::cr_fz_region=0.                    ! CR flux
     ! CR derived parameters.
     real(kind=8),dimension(1:ncrgrp)::cr_d_code ! CR diffusion coefficients in code units    !
     real(kind=8)::cr_dmax_code                  ! Max diffusion coefficient in code units    !

     ! Turbulence driving parameters
     logical  :: turb=.false.            ! Use turbulence?
     integer  :: turb_seed=-1            ! Turbulent seed (-1=random)
     character (LEN=100) :: forcing_power_spectrum='parabolic'
                                         ! Power spectrum type of turbulent forcing
     real(kind=8) :: comp_frac=0.3333    ! Compressive fraction
     real(kind=8) :: turb_T=1.0          ! Turbulent velocity autocorrelation time
     integer      :: turb_Ndt=100        ! Number of timesteps per autocorr. time
     real(kind=8) :: turb_rms=1.0        ! rms turbulent forcing acceleration
     real(kind=8) :: turb_min_rho=1d-50  ! Minimum density for turbulence

  end type run_t

  type hydro_params_t
     integer::slope_type,slope_mag_type,riemann,riemann2d
     real(kind=8),dimension(1:nener+1)::gamma_rad
     real(kind=8)::gamma,smallr,smallc
     real(kind=8)::difmag
     real(kind=8)::switch_llf_dmin,switch_llf_pmin
#ifdef MHD
     real(kind=8)::etamag
     logical::induction
#endif
  end type hydro_params_t

  type global_t

     ! MPI variables
     integer::ncpu, myid

     integer::iout=1             ! Increment for output times
     integer::ifout=1            ! Increment for output files
     integer::ifbkp=1            ! Increment for backup files

     logical::output_done=.false.                  ! Output just performed
     logical::init=.false.                         ! Set up or run
     integer::nstep=0                              ! Time step
     integer::nstep_coarse=0                       ! Coarse step
     integer::nstep_coarse_old=0                   ! Old coarse step
     integer::nflag,ncreate,nkill                  ! Refinements
     logical::first_coarse_restart =.false.        ! Flag that indicates first course time step after restart

     real(kind=8)::ekin_tot=0.0D0                      ! Total kinetic energy
     real(kind=8)::eint_tot=0.0D0                      ! Total internal energy
     real(kind=8)::emag_tot=0.0D0                      ! Total magnetic energy
     real(kind=8)::epot_tot=0.0D0                      ! Total potential energy
     real(kind=8)::epot_tot_old=0.0D0                  ! Old potential energy
     real(kind=8)::epot_tot_int=0.0D0                  ! Time integrated potential
     real(kind=8)::const=0.0D0                         ! Energy conservation
     real(kind=8)::aexp_old=1.0D0                      ! Old expansion factor
     real(kind=8)::rho_tot=0.0D0                       ! Mean density in the box
     real(kind=8)::t=0.0D0                             ! Time variable
     real(kind=8)::mass_tot=0.0D0                      ! Total mass in gas
     real(kind=8)::mass_star_tot=0.0D0                 ! Total mass in new stars
     real(kind=8)::mass_sink_tot=0.0D0                 ! Total mass in new sinks
     real(kind=8)::tot_nPhot=0.0D0, step_nPhot=0.0D0   ! RT bookkeeping
     real(kind=8)::step_nStar=0.0D0, step_mStar=0.0D0  ! RT bookkeeping

     real(kind=8)::mass_tot_0=0.0D0                    ! Initial total mass (gas+star+sink)

     ! Level related arrays
     real(kind=8),dimension(1:MAXLEVEL)::dtold,dtnew   ! Time step at each level
     real(kind=8),dimension(1:MAXLEVEL)::rho_max       ! Maximum density at each level
     integer,dimension(1:MAXLEVEL)::isubcycle      ! Current subcycling step at each level

     ! Only one process can write at a time in an I/O group
     integer::IOGROUPSIZE=0           ! Main snapshot
     integer::IOGROUPSIZECONE=0       ! Lightcone
     integer::IOGROUPSIZEREP=0        ! Subfolder size
     logical::withoutmkdir=.false.    ! If true mkdir should be done before the run
     logical::print_when_io=.false.   ! If true print when IO
     logical::synchro_when_io=.false. ! If true synchronize when IO

     ! Cosmology parameters
     real(kind=8)::boxlen_ini     ! Box size in h-1 Mpc
     real(kind=8)::omega_b=0.0D0  ! Omega Baryon
     real(kind=8)::omega_m=1.0D0  ! Omega Matter
     real(kind=8)::omega_l=0.0D0  ! Omega Lambda
     real(kind=8)::omega_k=0.0D0  ! Omega Curvature
     real(kind=8)::h0=1.0D0       ! Hubble constant in km/s/Mpc
     real(kind=8)::aexp=1.0D0     ! Current expansion factor
     real(kind=8)::hexp=0.0D0     ! Current Hubble parameter
     real(kind=8)::texp=0.0D0     ! Current proper time
     logical ::use_proper_time=.false.

     ! Executable identification
     CHARACTER(LEN=300)::builddate,buildcommand,patchdir
     CHARACTER(LEN=300)::gitrepo,gitbranch,githash

     ! Save namelist filename
     CHARACTER(LEN=300)::namelist_file

     ! Friedman model variables
     real(kind=8),dimension(0:n_frw)::aexp_frw,hexp_frw,tau_frw,t_frw

     ! Initial conditions parameters from grafic
     integer::nlevelmax_part
     real(kind=8)::aexp_ini=10.
     real(kind=8)::T2_start          ! Starting gas temperature
     real(kind=8),dimension(1:MAXLEVEL)::dfact=1.0d0,astart
     real(kind=8),dimension(1:MAXLEVEL)::vfact=1.0d0
     real(kind=8),dimension(1:MAXLEVEL)::xoff1,xoff2,xoff3,dxini
     integer ,dimension(1:MAXLEVEL)::n1,n2,n3

     ! Minimum particle mass
     real(kind=8) :: mp_min=-1.0

     ! Minimum MG level
     integer :: levelmin_mg

     ! Multigrid safety switch
     logical, dimension(1:MAXLEVEL)::safe_mode=.false.

     ! Multipole coefficients
     type(multipole_t)::multipole

     ! RT global variables
     real(kind=8),dimension(1:MAXLEVEL)::rt_c=1d0            ! Reduced lightspeed in code units
     real(kind=8),dimension(1:MAXLEVEL)::rt_c_cgs            ! Reduced lightspeed in [cm s-1]

     ! CR global variables
     real(kind=8),dimension(1:MAXLEVEL)::cr_c=1d0            ! Reduced lightspeed in code units
     !real(kind=8),dimension(1:MAXLEVEL)::cr_c_cgs            ! Reduced lightspeed in [cm s-1]

  end type global_t

  type mesh_t
     ! For GPU, is data on device
     logical::data_on_device=.false.

     ! Level related arrays
     integer(kind=4),allocatable,dimension(:)::head      ! Starting index for each level
     integer(kind=4),allocatable,dimension(:)::tail      ! Final index for each level
     integer(kind=4),allocatable,dimension(:)::noct      ! Number of octs for each level
     integer::ifree                                      ! Index of first oct in free memory

     ! Cache related arrays
     integer(kind=4),allocatable,dimension(:)::head_cache
     integer(kind=4),allocatable,dimension(:)::tail_cache
     integer(kind=4),allocatable,dimension(:)::noct_cache
     integer::ifree_cache

     integer(kind=4),allocatable,dimension(:)::noct_min  ! Min. number of octs across cpus
     integer(kind=4),allocatable,dimension(:)::noct_max  ! Max. number of octs across cpus
     integer(kind=8),allocatable,dimension(:)::noct_tot  ! Total number of octs across cpus

     integer(kind=4),allocatable,dimension(:)::ckey_max        ! Max. Cartesian key per level
     integer(kind=8),allocatable,dimension(:)::key_off         ! Key offset for GPU nly
     integer(kind=8),allocatable,dimension(:,:)::hkey_max      ! Max. Hilbert key per level
     integer(kind=4),allocatable,dimension(:,:)::box_ckey_min  ! Min. Cartesian key per level for the box
     integer(kind=4),allocatable,dimension(:,:)::box_ckey_max  ! Max. Cartesian key per level for the box
     integer(kind=4),allocatable,dimension(:,:,:)::bound_ckey_min  ! Min. Cartesian key per level for the boundaries
     integer(kind=4),allocatable,dimension(:,:,:)::bound_ckey_max  ! Max. Cartesian key per level for the boundaries

     integer(kind=4)::noct_used,noct_used_max,noct_used_tot ! Total used octs in local memory

     integer(kind=4)::nx,ny,nz                   ! Size of mesh at levelmin
     real(kind=8),allocatable,dimension(:)::skip ! Coordinates of lower left corner of the box

     ! Size of the oct array
     integer(kind=4)::ngridmax

     ! Persistent array for the grid
     type(oct),allocatable,dimension(:)::grid

     ! Grid hash table
     type(hash_table)::grid_dict
     integer(kind=4)::hash_size
     integer(kind=4)::hash_used

     ! Grid variables
     integer,allocatable,dimension(:,:)::flag1
     integer,allocatable,dimension(:,:)::flag2
#ifdef HYDRO
     real(dp),allocatable,dimension(:,:,:)::uold
     real(dp),allocatable,dimension(:,:,:)::unew
#ifdef TRCFLX
     real(dp),allocatable,dimension(:,:,:)::mflux      ! Time-integrated mass flux for tracers
#endif
#endif
#ifdef MHD
     real(dp),allocatable,dimension(:,:,:)::bold
     real(dp),allocatable,dimension(:,:,:)::bnew
#endif
#ifdef DO_RT
     real(dp),allocatable,dimension(:,:,:)::rtuold
     real(dp),allocatable,dimension(:,:,:)::rtunew
     real(dp),allocatable,dimension(:,:,:)::emissivity
#endif
#ifdef DO_CR
     real(dp),allocatable,dimension(:,:,:)::cruold
     real(dp),allocatable,dimension(:,:,:)::crunew
#endif
#ifdef TURB
     real(dp),allocatable,dimension(:,:,:)::fturb
#endif
#ifdef GRAV
     real(dp),allocatable,dimension(:,:,:)::f
     real(dp),allocatable,dimension(:,:)::rho
     real(dp),allocatable,dimension(:,:)::phi
     real(dp),allocatable,dimension(:,:)::phi_old
     real(dp),allocatable,dimension(:,:)::nref
#endif

     ! Clean/dirty octs for first neighbors
     integer(kind=4),allocatable,dimension(:)::indx_clean, head_clean, tail_clean, noct_clean
     integer(kind=4),allocatable,dimension(:)::indx_dirty, head_dirty, tail_dirty, noct_dirty

     ! Software cache array for the AMR grid
     logical,allocatable,dimension(:)::dirty
     logical,allocatable,dimension(:)::occupied
     logical,allocatable,dimension(:)::locked
     integer,allocatable,dimension(:)::parent_cpu
     integer,allocatable,dimension(:)::ghost_parent_grid
     integer,allocatable,dimension(:)::ghost_parent_cell
     integer::free_cache,ncache,ncachemax,nlocked,nlocked_max

     ! Software cache array for failed requests
     logical,allocatable,dimension(:)::occupied_null
     integer,allocatable,dimension(:)::lev_null
     integer,allocatable,dimension(:,:)::ckey_null
     integer::free_null,nnull

     ! Peano-Hilbert key boundaries for cpu domains
     type(domain_t),allocatable,dimension(:)::domain

     ! Hydro kernel workspace
     type(hydro_workspace_t)::hydro_w

     ! RT kernel workspace
     type(rt_workspace_t)::rt_w

     ! CR kernel workspace
     type(cr_workspace_t)::cr_w

  end type mesh_t

contains

  subroutine print_run_parameters(run_p)
    class(run_t)::run_p
    type(run_t)::run_params
    namelist /run_parameters/ run_params
    run_params = run_p
    write(*,NML=run_parameters)
  end subroutine print_run_parameters

  subroutine set_hydro_parameters(r,h_params)
    type(run_t)::r
    type(hydro_params_t)::h_params
    ! Transfer hydro parameters
    h_params%gamma=r%gamma
    h_params%gamma_rad=r%gamma_rad
    h_params%smallr=r%smallr
    h_params%smallc=r%smallc
    h_params%slope_type=r%slope_type
    h_params%slope_mag_type=r%slope_mag_type
    h_params%riemann=r%riemann
    h_params%riemann2d=r%riemann2d
    h_params%difmag=r%difmag
    h_params%switch_llf_dmin=r%switch_llf_dmin
    h_params%switch_llf_pmin=r%switch_llf_pmin
#ifdef MHD
    h_params%etamag=r%etamag
    h_params%induction=r%induction
#endif
  end subroutine set_hydro_parameters

end module amr_commons
