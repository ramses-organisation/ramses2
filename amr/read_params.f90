module params_module

contains

subroutine m_read_params(pst)
  use amr_parameters
  use hydro_parameters
  use rt_parameters, only: nrtgrp
  use cr_parameters, only: ncrgrp
  use ramses_commons, only: pst_t
  use mdl_module
  use movie_module, only: set_movie_vars
  use rt_params_module
  use cr_params_module
  use constants
#ifdef DO_RTZ
  use rtz_module, only: elements, n_elements !, initialize_elements
#endif
  implicit none
  type(pst_t)::pst

  !--------------------------------------------------
  ! Local variables
  !--------------------------------------------------
  integer::i,narg,levelmax
  character(LEN=80)::infile
  character(LEN=80)::cmdarg
  integer(kind=8)::ngridtot=0
  integer(kind=8)::nparttot=0
  integer(kind=8)::nstartot=0
  integer(kind=8)::nsinktot=0
  integer(kind=8)::ntreetot=0
  integer(kind=8)::ntractot=0
  integer(kind=8)::ndusttot=0
  real(kind=8)::delta_tout=0,tend=0
  real(kind=8)::delta_aout=0,aend=0
  logical::nml_ok

  !--------------------------------------------------
  ! Namelist variables
  !--------------------------------------------------

  ! Run control
  logical::cosmo   =.false.    ! Cosmology activated
  logical::pic     =.false.    ! Particle In Cell activated
  logical::poisson =.false.    ! Poisson solver activated
  logical::hydro   =.false.    ! Hydro activated
  logical::rt      =.false.    ! RT activated
  logical::cr      =.false.    ! Cosmic rays activated
  logical::star    =.false.    ! Stars and star formation activated
  logical::sink    =.false.    ! Sinks and sink formation activated
  logical::part    =.false.    ! Dark matter particles activated
  logical::trac    =.false.    ! Tracer particles activated
  logical::dust    =.false.    ! Dust particles activated
  logical::merger_tree=.false. ! Merger tree particles activated
  logical::orphan  =.false.    ! Orphan particles activated
  logical::verbose =.false.    ! Write everything
  logical::debug   =.false.    ! Debug mode activated
  logical::static_mesh=.false. ! Static mesh refinement activated
  logical::static_gas=.false.  ! Hydro is turned off
  logical::clump_only=.false.  ! Only clump finding

  ! Step parameters
  integer::nrestart=0         ! New run or backup file number
  integer::nstepmax=1000000   ! Maximum number of time steps
  integer::ncontrol=1         ! Write control variables
  integer::nremap=10          ! Particle load balancing frequency (0: never)

  ! Maximum number of allocatable particles
  integer::npartmax=0
  integer::nstarmax=0
  integer::nsinkmax=0
  integer::ntreemax=0
  integer::ntracmax=0
  integer::ndustmax=0

  ! IC subcell multiplicity
  integer::ntrac_per_cell=1
  integer::ndust_per_cell=1
  logical :: part_subcell_positions=.true.

  ! Dust parameters
  real(kind=8)::dust_to_gas_mass_ratio=0.0d0
  real(kind=8)::grain_size_parameter=0.0d0
  real(kind=8)::grain_charge_parameter=0.0d0
  real(kind=8)::dust_gyro_factor=0.1d0 ! At least 10 steps per gyro-period.
  logical :: analytic_dust_force = .false.

  ! Tracer parameters
  character(LEN=32)::tracer_kick_pdf='piecewise_skew_uniform' ! Tracer kick PDF

  ! Number of superoct levels
  integer::nsuperoct=0

  ! MPI domain overloading
  integer::overload=1

  ! Mesh parameters
  integer::geom=1             ! 1: cartesian, 2: cylindrical, 3: spherical
  integer::levelmin=1         ! Full refinement up to levelmin
  integer::nlevelmax=1        ! Maximum number of level
  integer::ngridmax=0         ! Maximum number of grids
  integer::ncachemax=10000    ! Maximum number of cache lines
  real(kind=8)::boxlen=1.0D0      ! Cell sixe at level 0

  ! Output parameters
  integer::noutput=1          ! Total number of outputs
  integer::foutput=1000000    ! Frequency of outputs
  integer::output_mode=0      ! Output mode (for hires runs)
  logical::gadget_output=.false. ! Output in gadget format
  real(kind=8)::bkp_time_hrs=2   ! Backup file frequency in hours
  real(kind=8)::run_time_hrs=0   ! Estimated run time in hrs
  real(kind=8)::bkp_last_min=10  ! Backup file before the end of run in min
  integer::bkp_modulo=0          ! Use modulo for backup file count
  integer::nfile=1               ! Number of file per snapshot. Use -1 for nfile=ncpu
  logical::output_part=.true. ! Output particle data in regular dumps (default: true)
  logical::output_grav=.true. ! Output gravity data in regular dumps (default: true)
  logical::output_hydro=.true.! Output hydro data in regular dumps (default: true)
  logical::output_amr=.true.  ! Output AMR data in regular dumps (default: true)

  ! Output times
  real(kind=8),dimension(1:MAXOUT)::aout=1.1  ! Output expansion factors
  real(kind=8),dimension(1:MAXOUT)::tout=0.0  ! Output times

  ! Trajectory output parameters
  integer::ntrajectories=0
  integer,dimension(1:MAXOUT)::trajectories=0

  ! Movie
  integer::imovout=0     ! Increment for output times
  integer::imov=1        ! Initialize
  real(kind=8)::tstartmov=0d0,astartmov=0d0
  real(kind=8)::tendmov=0.,aendmov=0.
  logical::movie=.false.
  logical::zoom_only=.false.
  integer::nw_frame=512 ! width of frame in pixels
  integer::nh_frame=512 ! height of frame in pixels
  integer::levelmax_frame=0
  integer::ivar_frame=1
  real(kind=8),dimension(1:20)::xcentre_frame=0d0
  real(kind=8),dimension(1:20)::ycentre_frame=0d0
  real(kind=8),dimension(1:20)::zcentre_frame=0d0
  real(kind=8),dimension(1:10)::deltax_frame=0d0
  real(kind=8),dimension(1:10)::deltay_frame=0d0
  real(kind=8),dimension(1:10)::deltaz_frame=0d0
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

  ! Refinement parameters for each level
  integer ,dimension(1:MAXLEVEL)::nexpand = 1 ! Number of mesh expansion
  integer ,dimension(1:MAXLEVEL)::nsubcycle = 2 ! Subcycling at each level
  real(kind=8),dimension(1:MAXLEVEL)::m_refine =-1.0 ! Lagrangian threshold
  real(kind=8),dimension(1:MAXLEVEL)::r_refine =-1.0 ! Radius of refinement region
  real(kind=8),dimension(1:MAXLEVEL)::x_refine = 0.0 ! Center of refinement region
  real(kind=8),dimension(1:MAXLEVEL)::y_refine = 0.0 ! Center of refinement region
  real(kind=8),dimension(1:MAXLEVEL)::z_refine = 0.0 ! Center of refinement region
  real(kind=8),dimension(1:MAXLEVEL)::exp_refine = 2.0 ! Exponent for distance
  real(kind=8),dimension(1:MAXLEVEL)::a_refine = 1.0 ! Ellipticity (Y/X)
  real(kind=8),dimension(1:MAXLEVEL)::b_refine = 1.0 ! Ellipticity (Z/X)
  real(kind=8)::var_cut_refine=-1.0 ! Threshold for variable-based refinement
  real(kind=8)::mass_cut_refine=-1.0 ! Mass threshold for particle-based refinement
  integer::ivar_refine=-1 ! Variable index for refinement
  real(kind=8)::aexp_lock_refine=-1.0
  logical::pic_lock_refine=.false.

  ! Default units
  real(kind=8)::units_density=1.0 ! [g/cm^3]
  real(kind=8)::units_time=1.0    ! [seconds]
  real(kind=8)::units_length=1.0  ! [cm]
  real(kind=8)::units_np=1.0      ! [#photon/cm^3]

  ! Initial conditions parameters from grafic
  real(kind=8)::aexp_ini=10.
  real(kind=8)::boxlen_ini=1.0
  real(kind=8)::omega_b=0.045
  real(kind=8)::omega_m=1.0
  real(kind=8)::omega_l=0.0
  real(kind=8)::h0=1.0

  ! Initial condition regions parameters
  integer::nregion=0
  character(LEN=10),dimension(1:MAXREGION)::region_type='square'
  real(kind=8),dimension(1:MAXREGION)::x_center=0.
  real(kind=8),dimension(1:MAXREGION)::y_center=0.
  real(kind=8),dimension(1:MAXREGION)::z_center=0.
  real(kind=8),dimension(1:MAXREGION)::length_x=1.E10
  real(kind=8),dimension(1:MAXREGION)::length_y=1.E10
  real(kind=8),dimension(1:MAXREGION)::length_z=1.E10
  real(kind=8),dimension(1:MAXREGION)::exp_region=2.0

  ! Initial condition files for each level
  logical::multiple=.false.
  character(LEN=20)::filetype='ascii'
  character(LEN=80),dimension(1:MAXLEVEL)::initfile=' '
  real(kind=8)::ic_scale_l=1.0d0
  real(kind=8)::ic_scale_v=1.0d0
  real(kind=8)::ic_scale_m=1.0d0

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

  ! Initial condition rt variables
#ifdef DO_RT
  integer::rt_nregion=0
  character(LEN=10),dimension(1:MAXREGION)::rt_region_type='square'
  real(kind=8),dimension(1:MAXREGION)::rt_reg_x_center=0.
  real(kind=8),dimension(1:MAXREGION)::rt_reg_y_center=0.
  real(kind=8),dimension(1:MAXREGION)::rt_reg_z_center=0.
  real(kind=8),dimension(1:MAXREGION)::rt_reg_length_x=1.E10
  real(kind=8),dimension(1:MAXREGION)::rt_reg_length_y=1.E10
  real(kind=8),dimension(1:MAXREGION)::rt_reg_length_z=1.E10
  real(kind=8),dimension(1:MAXREGION)::rt_exp_region=2.0
  integer ,dimension(1:MAXREGION)::rt_reg_group=1
  real(kind=8),dimension(1:MAXREGION)::rt_n_region=0.0 ! Photon density
  real(kind=8),dimension(1:MAXREGION)::rt_u_region=0.0 !    Photon flux
  real(kind=8),dimension(1:MAXREGION)::rt_v_region=0.0 !    Photon flux
  real(kind=8),dimension(1:MAXREGION)::rt_w_region=0.0 !    Photon flux
#endif

  ! Initial condition CR variables
#ifdef DO_CR
  integer::cr_nregion=0
  character(LEN=10),dimension(1:MAXREGION)::cr_region_type='square'
  real(kind=8),dimension(1:MAXREGION)::cr_reg_x_center=0.
  real(kind=8),dimension(1:MAXREGION)::cr_reg_y_center=0.
  real(kind=8),dimension(1:MAXREGION)::cr_reg_z_center=0.
  real(kind=8),dimension(1:MAXREGION)::cr_reg_length_x=1.E10
  real(kind=8),dimension(1:MAXREGION)::cr_reg_length_y=1.E10
  real(kind=8),dimension(1:MAXREGION)::cr_reg_length_z=1.E10
  real(kind=8),dimension(1:MAXREGION)::cr_exp_region=2.0
  integer ,dimension(1:MAXREGION)::cr_reg_group=1
  real(kind=8),dimension(1:MAXREGION)::cr_fx_region=0.0 ! CR flux
  real(kind=8),dimension(1:MAXREGION)::cr_fy_region=0.0 ! CR flux
  real(kind=8),dimension(1:MAXREGION)::cr_fz_region=0.0 ! CR flux
#endif

  ! Refinement parameters for hydro
  real(kind=8)::err_grad_d=-1.0  ! Density gradient
  real(kind=8)::err_grad_u=-1.0  ! Velocity gradient
  real(kind=8)::err_grad_p=-1.0  ! Pressure gradient
  real(kind=8)::err_grad_xHI=-1.0! xHI gradient
  real(kind=8)::err_grad_xHII=-1.0 ! xHII gradient
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
  real(kind=8),dimension(1:MAXLEVEL)::jeans_refine=-1.0

  ! Refinement parameters for rt
#ifdef DO_RT
  real(kind=8)::rt_err_grad_cn(nrtgrp)=-1 ! Photon flux gradient for refinement
  real(kind=8)::rt_floor_cn(nrtgrp)=1d-10 ! Photon flux floor for refinement
  real(kind=8)::rt_refine_aexp=-1.0      ! Start expansion factor for RT refinements
#endif

  ! Hydro solver parameters
  integer ::niter_riemann=10
  integer ::slope_type=1
  integer ::slope_mag_type=-1
  real(kind=8)::gamma=1.4d0
  real(kind=8),dimension(1:512)::gamma_rad=1.33333333334d0
  real(kind=8)::courant_factor=0.5d0
  real(kind=8)::difmag=0.0d0
  real(kind=8)::etamag=0.0d0
  real(kind=8)::smallc=1.d-10
  real(kind=8)::smallr=1.d-10
  character(LEN=10)::scheme='muscl'
  character(LEN=10)::riemann='llf'
  character(LEN=10)::riemann2d='none'
  logical ::induction=.false.
  logical ::entropy=.false.
  logical ::sgs_turb=.false.
  logical ::equilibrium_sgs=.false.
  real(kind=8)::dual_energy=-1
  real(kind=8)::T2_fix=0d0
  real(kind=8),dimension(1:3)::constant_gravity=0.0d0
  real(kind=8)::switch_llf_dmin=-1
  real(kind=8)::switch_llf_pmin=-1
  real(kind=8)::smagorinsky_lilly_constant=1.0d0

  ! Non-thernal energies and passive scalars index
  integer ::inener,ientropy,imetal,iturb,ichem

  ! Interpolation parameters
  integer ::interpol_var=0
  integer ::interpol_type=1

  ! Poisson solver parameters
  logical :: gravity_test=.false. ! Use file rho_ana.f90 to test the Poisson solers.
  real(kind=8)::epsilon=1.0D-4    ! Convergence criterion
  integer :: nvcycle = -1         ! Desired number of V-cycles
  real(kind=8),dimension(1:10)::gravity_params=0.0 ! Gravity parameters
  integer :: gravity_type=0 ! Type of gravity calculations (see user guide)
  integer :: cic_levelmax=0 ! Maximum level for CIC dark matter interpolation
  integer :: cg_levelmin=999   ! Min level for CG solver
  ! level < cg_levelmin uses fine multigrid
  ! level >=cg_levelmin uses conjugate gradient
  logical :: fast_solver=.false.   ! Fast solver with MPI pre-fetch (memory intensive)
  integer :: part_mass_deposition_scheme=1     ! part mass deposition schemes (CIC 1, TSC 2, PCS 3)
  integer :: part_dep_algo=2                   ! part GPU CIC deposition algorithm
  integer :: part_force_interpolation_scheme=1 ! part force interpolation schemes (CIC 1, TSC 2, PCS 3)
  integer :: star_mass_deposition_scheme=1     ! star mass deposition schemes
  integer :: star_dep_algo=2                   ! star GPU CIC deposition algorithm
  integer :: star_force_interpolation_scheme=1 ! star force interpolation schemes
  integer :: sink_mass_deposition_scheme=1     ! sink mass deposition schemes
  integer :: sink_dep_algo=2                   ! sink GPU CIC deposition algorithm
  integer :: sink_force_interpolation_scheme=1 ! sink force interpolation schemes
  integer :: tree_mass_deposition_scheme=1     ! tree mass deposition schemes
  integer :: tree_force_interpolation_scheme=1 ! tree force interpolation schemes
  integer :: trac_interpolation_scheme=1 ! tracer interpolation/numerical schemes
  integer :: dust_mass_deposition_scheme=1 ! dust mass deposition schemes
  integer :: dust_force_interpolation_scheme=1 ! dust force interpolation/numerical schemes

  ! Boundary conditions parameters
  integer::nbound=0
  logical::no_inflow=.false.
  logical,dimension(1:3)::periodic=.true.
  integer::bound_levelmin=1 ! AMR level for boundary geometry
  real(kind=8),dimension(1:3)::box_size=0.0D0 ! Box length in active domain along each direction
  integer::box_xmin=0 ! Min. Cartesian key for the box at levelmin in x direction
  integer::box_xmax=0 ! Max. Cartesian key for the box at levelmin in x direction
  integer::box_ymin=0 ! Min. Cartesian key for the box at levelmin in y direction
  integer::box_ymax=0 ! Max. Cartesian key for the box at levelmin in y direction
  integer::box_zmin=0 ! Min. Cartesian key for the box at levelmin in z direction
  integer::box_zmax=0 ! Max. Cartesian key for the box at levelmin in z direction
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
  real(kind=8),dimension(1:MAXBOUND)::p_bound=0
  real(kind=8),dimension(1:MAXBOUND)::u_bound=0
  real(kind=8),dimension(1:MAXBOUND)::v_bound=0
  real(kind=8),dimension(1:MAXBOUND)::w_bound=0
#if NENER>0
  real(kind=8),dimension(1:MAXBOUND,1:NENER)::prad_bound=0
#endif
#if NVAR>5+NENER
  real(kind=8),dimension(1:MAXBOUND,1:NVAR-5-NENER)::var_bound=0
#endif
#ifdef DO_RT
  real(kind=8),dimension(1:MAXBOUND,1:nrtgrp)::rt_n_bound=0.0d0
  real(kind=8),dimension(1:MAXBOUND,1:nrtgrp)::rt_u_bound=0.0d0
  real(kind=8),dimension(1:MAXBOUND,1:nrtgrp)::rt_v_bound=0.0d0
  real(kind=8),dimension(1:MAXBOUND,1:nrtgrp)::rt_w_bound=0.0d0
#endif
#ifdef DO_CR
  real(kind=8),dimension(1:MAXBOUND,1:ncrgrp)::cr_fx_bound=0.0d0
  real(kind=8),dimension(1:MAXBOUND,1:ncrgrp)::cr_fy_bound=0.0d0
  real(kind=8),dimension(1:MAXBOUND,1:ncrgrp)::cr_fz_bound=0.0d0
#endif

  ! Cooling parameters
  logical::cooling=.false.
  logical::cooling_ism=.false.
  logical::metal=.false.
  logical::isothermal=.false. ! Force temperature to eos value
  logical::haardt_madau=.false.
  logical::self_shielding=.false.
  real(kind=8)::J21=0,a_spec=1,z_ave=0,z_reion=8.5,cooling_uvb_delta=0.05d0
  integer::eos_type=1 ! 1=isothermal, 2=polytrope, 3=isothermal+polytrope
  real(kind=8)::eos_nH=huge(1d0),eos_index=1,eos_T2=10
  real(kind=8)::T2max=huge(1d0)
  logical::neq_chem=.false. ! Non-equilibrium cooling -------------------
  logical::is_init_xion=.false.   ! Initialize ionization from T profile?
  logical::upload_equilibrium_x=.false.! Enforce equ. xion when uploading
  logical::isHe=.true.            !      He ionization fractions tracked?
  logical::isH2=.false.           !                 H2 tracked (via xHI)?
  real(kind=8)::neq_Tconst=-1         !             Const T in neq chemistry?
  real(kind=8) ::mu_mol=1.2195d0  ! Mean molecular weight (std ISM value)
  real(kind=8) ::X_H=0.7600d0     !                Hydrogen mass fraction
  real(kind=8) ::Y_He=0.2400d0    !                  Helium mass fraction
  integer::iIons,ixHI=0,ixHII=0,ixHeII=0,ixHeIII=0 !   Ionization indices
#ifdef DO_RTZ
  integer::element_first_idx(1:27)       !  Start idx of elements in uold
  integer::molecules_first_idx(1:3)     !  Start idx of molecules in uold
  real(kind=8),dimension(1:27,1:27)::ionEvs=0.     !  Ionization energies
#else
  real(kind=8),dimension(nion)::ionEvs=0.          !  Ionization energies
#endif
  integer::icount

  ! RTZ cooling parameters
  logical::rtz_cooling=.false. ! RTZ Non-equilibrium cooling ------------
  integer::rtz_equilibrium_test=-1 ! RTZ EQM test ------------------
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
  integer::rtz_eqm_min_its

  ! Star formation parameters
  integer::sf_model=1
  real(kind=8)::T2_star=2e4
  real(kind=8)::n_star=0.1
  real(kind=8)::eps_star=0.01
  integer(kind=8),dimension(1:6)::seed=(/123,456,789,1,1,1/)
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
  integer::rho_type_clump=1 ! 1: DM, 2: stars, 3: sinks, 4: gas
  integer::nsteps_per_tree=1 ! Call clump finder for tree formation every n coarse steps
  real(kind=8)::relevance_threshold=2
  real(kind=8)::density_threshold=-1
  real(kind=8)::saddle_threshold=-1
  real(kind=8)::mass_threshold=0
  real(kind=8)::purity_threshold=-1
  real(kind=8)::fraction_threshold=2d0

  ! Lightcone parameters
  logical :: lightcone = .false.   ! Lightcone activated
  real(kind=8) :: cone_z_min = 0.0 ! Minimum redshift
  real(kind=8) :: cone_z_max = 1.0 ! Maximum redshift
  real(kind=8) :: cone_opening_angle_y = 0.0 ! Opening angle in y direction in degrees
  real(kind=8) :: cone_opening_angle_z = 0.0 ! Opening angle in z direction in degrees
  real(kind=8) :: cone_theta = 0.0 ! Rotation of the cone's x-axis around the box's y-axis in degrees
  real(kind=8) :: cone_phi = 0.0 ! Rotation of the cone's x-axis around the box's z-axis in degrees
  real(kind=8), dimension(1:3) :: cone_observer = (/0.0, 0.0, 0.0/) ! Observer position in code units

  ! Sink formation/dynamical parameters
  integer::rho_type_sink=1
  logical::sink_descent=.false.
  real(kind=8)::fudge_descent=0.5d0
  real(kind=8)::sink_relevance_threshold=2
  real(kind=8)::sink_density_threshold=-1
  real(kind=8)::sink_saddle_threshold=-1
  real(kind=8)::sink_mass_threshold=0
  real(kind=8)::sink_purity_threshold=-1
  real(kind=8)::sink_fraction_threshold=2d0
  real(kind=8)::sink_nstar_frac=-1 ! Gas density threshold in units of SF density
  real(kind=8)::sink_delta_tout=0 ! Time interval in code units between each sink high frequency dump
  real(kind=8)::sink_mseed=1e6 ! Sink seed mass in solar masses
  logical::sink_form=.false.
  logical::sink_merge=.false.
  logical::sink_refine=.false.
  logical::sink_dump=.false.
  logical::static_sink=.false.
  logical::fix_sink_mass=.false.
  logical::drag_sink=.false. ! Whether to use dynamical friction for black hole dynamics

  ! Sink accretion parameters
  integer::accretion_type = 0 ! 0: None, 1: Bondi, 2: Flux
  real(kind=8)::acc_sink_boost = 1.0d0 ! Boost for bondi accretion
  logical::bondi_use_vrel = .true. ! Whether to use the relative sink velocity for BHL accretion
  real(kind=8)::eddington_cap = -1 ! Factor of Eddington rate to cap accretion at
  integer::sink_b_spline_order = 4 ! Order of B-spline interpolation used for sink accretion and dynamics
  logical::verbose_sink = .false. ! Whether to print verbose statements for sink particles
  logical::bondi_use_gas_mass = .false. ! Whether to include the local gas mass in the Bondi calculation
  logical::use_local_bondi_rate = .false. ! Switch to average after (true) or before (false) computing the Bondi rate
  logical::use_rho_inf = .true. ! Whether to use bondi_alpha(x) to extrapolate density at infinity from Bondi solution
  real(kind=8)::t_start_black_hole = -1 ! Time after which to start using sink particle/black hole routines (code units)
  logical::use_bondi_lambda = .true.
  logical::mass_weighting = .true.
  logical::momentum_conserving = .false.
  real(kind=8)::eddington_floor = -1 ! Accretion rate floor below which nothing happens
  real(kind=8)::manual_accretion_rate = -1 ! Manual accretion rate (fraction of Eddington)

  ! Sink feedback parameters
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

#ifdef DO_RTZ
  integer::i_Element, i_Iion
#endif

  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------
  ! Global run parameter
  namelist/run_params/cosmo,pic,poisson,hydro,rt,cr &
       & ,verbose,debug,nrestart,ncontrol,nstepmax,nsubcycle,nremap &
       & ,static_mesh,static_gas,geom,overload,nsuperoct &
       & ,clump_only
  ! Output parameters
  namelist/output_params/foutput,aout,tout,output_mode &
       & ,tend,delta_tout,aend,delta_aout,gadget_output &
       & ,run_time_hrs,bkp_time_hrs,bkp_last_min,bkp_modulo,nfile &
       & ,output_part,output_grav,output_hydro,output_amr
  ! Trajectory output parameters
  namelist/traj_params/ntrajectories,trajectories
  ! AMR grid basic parameters
  namelist/amr_params/levelmin,levelmax,ngridmax,ncachemax,ngridtot &
       & ,npartmax,nparttot,nexpand,boxlen
  ! Poisson solver parameters
  namelist/poisson_params/epsilon,nvcycle,gravity_type,gravity_params &
       & ,cg_levelmin,cic_levelmax,fast_solver,gravity_test &
       & ,part_mass_deposition_scheme,part_dep_algo,part_force_interpolation_scheme &
       & ,star_mass_deposition_scheme,star_dep_algo,star_force_interpolation_scheme &
       & ,sink_mass_deposition_scheme,sink_dep_algo,sink_force_interpolation_scheme &
       & ,tree_mass_deposition_scheme,tree_force_interpolation_scheme
  ! Movies parameters
  namelist/movie_params/levelmax_frame,nw_frame,nh_frame,ivar_frame &
       & ,xcentre_frame,ycentre_frame,zcentre_frame &
       & ,deltax_frame,deltay_frame,deltaz_frame,movie,zoom_only,zoom_only_frame &
       & ,imovout,imov,tstartmov,astartmov,tendmov,aendmov,proj_axis &
       & ,movie_vars,movie_vars_txt &
       & ,theta_camera,phi_camera,dtheta_camera,dphi_camera &
       & ,tstart_theta_camera,tstart_phi_camera,tend_theta_camera,tend_phi_camera &
       & ,focal_camera,dist_camera,ddist_camera &
       & ,perspective_camera,smooth_frame,shader_frame,method_frame &
       & ,varmin_frame,varmax_frame
  ! Initial conditions parameters
  namelist/init_params/filetype,initfile,multiple,nregion,region_type &
       & ,x_center,y_center,z_center,aexp_ini,omega_b,omega_m,omega_l,h0 &
       & ,length_x,length_y,length_z,exp_region,boxlen_ini &
       & ,ic_scale_l,ic_scale_v,ic_scale_m &
#if NENER>0
       & ,prad_region &
#endif
#if NVAR>5+NENER
       & ,var_region &
#endif
#ifdef MHD
       & ,B_region,C_region,A_ave,B_ave,C_ave &
#endif
#ifdef DO_RT
       & ,rt_nregion, rt_region_type                           &
       & ,rt_reg_x_center, rt_reg_y_center, rt_reg_z_center    &
       & ,rt_reg_length_x, rt_reg_length_y, rt_reg_length_z    &
       & ,rt_exp_region, rt_reg_group                          &
       & ,rt_n_region, rt_u_region, rt_v_region, rt_w_region   &
#endif
#ifdef DO_CR
       & ,cr_nregion, cr_region_type                             &
       & ,cr_reg_x_center, cr_reg_y_center, cr_reg_z_center      &
       & ,cr_reg_length_x, cr_reg_length_y, cr_reg_length_z      &
       & ,cr_exp_region, cr_reg_group                            &
       & ,cr_fx_region, cr_fy_region, cr_fz_region  &
#endif
       & ,d_region,u_region,v_region,w_region,p_region
  ! Hydro solver parameters
  namelist/hydro_params/gamma,courant_factor,smallr,smallc &
       & ,slope_type,slope_mag_type,difmag,etamag,gamma_rad &
       & ,dual_energy,T2_fix,induction,entropy,sgs_turb,equilibrium_sgs,riemann,riemann2d,constant_gravity &
       & ,niter_riemann,scheme,switch_llf_dmin,switch_llf_pmin,smagorinsky_lilly_constant
  ! Grid refinement parameters
  namelist/refine_params/x_refine,y_refine,z_refine,r_refine &
       & ,a_refine,b_refine,exp_refine,jeans_refine,mass_cut_refine &
#ifdef MHD
       & ,err_grad_b2,err_grad_A,err_grad_B,err_grad_C &
       & ,floor_b2,floor_A,floor_B,floor_C &
#endif
#if NENER>0
       & ,err_grad_prad, floor_prad &
#endif
#if NVAR>5+NENER
       & ,err_grad_var &
#endif
#ifdef DO_RT
       & ,rt_err_grad_cn, rt_floor_cn, rt_refine_aexp &
#endif
       & ,err_grad_xHI, err_grad_xHII, floor_xHI, floor_xHII &
       & ,m_refine,mass_sph,err_grad_d,err_grad_p,err_grad_u &
       & ,floor_d,floor_u,floor_p,ivar_refine,var_cut_refine &
       & ,interpol_var,interpol_type &
       & ,aexp_lock_refine,pic_lock_refine,sink_refine
  ! Units parameters
  namelist/units_params/units_density,units_time,units_length,units_np
  ! Boundary conditions parameters
  namelist/boundary_params/periodic,nbound,bound_type,bound_dir,bound_shift &
       & ,bound_xmin,bound_xmax,bound_ymin,bound_ymax,bound_zmin,bound_zmax &
       & ,bound_levelmin,box_size,box_xmin,box_xmax,box_ymin,box_ymax,box_zmin,box_zmax &
#if NENER>0
       & ,prad_bound &
#endif
#if NVAR>5+NENER
       & ,var_bound &
#endif
#ifdef DO_RT
       & ,rt_n_bound,rt_u_bound,rt_v_bound,rt_w_bound &
#endif
#ifdef DO_CR
       & ,cr_fx_bound,cr_fy_bound,cr_fz_bound &
#endif
       & ,d_bound,u_bound,v_bound,w_bound,p_bound
  ! Cooling / basic chemistry parameters
  namelist/cooling_params/neq_chem,cooling,metal,isothermal,haardt_madau,J21 &
       & ,eos_type,eos_nH,eos_index,eos_T2, mu_mol, X_H, Y_He &
       & ,a_spec,self_shielding,z_ave,z_reion,T2max,cooling_ism &
       & ,cooling_uvb_delta &
       & ,isHe, isH2, is_init_xion, neq_Tconst, upload_equilibrium_x &
       & ,rtz_cooling, rtz_equilibrium_test, rtz_include_collisional_ionization &
       & ,rtz_include_photoionization, rtz_include_cosmic_ray_ionization &
       & ,rtz_include_charge_exchange, rtz_include_dust_recombination, rtz_UV_background_G0 &
       & ,rtz_primary_cosmic_ray_ionization_rate, rtz_include_HM12_UVB, isH2_rtz &
       & ,rtz_max_cool_timestep, rtz_eqm_min_its
  ! Tracer particles parameters
  namelist/trac_params/trac,ntracmax,ntractot,ntrac_per_cell,trac_interpolation_scheme,&
       & part_subcell_positions,tracer_kick_pdf
  namelist/dust_params/dust,ndustmax,ndusttot,ndust_per_cell,dust_to_gas_mass_ratio,&
       & grain_size_parameter,grain_charge_parameter,&
       & dust_mass_deposition_scheme,dust_force_interpolation_scheme,dust_gyro_factor,&
       & analytic_dust_force
  ! Star particles and star formation recipe
  namelist/star_params/star,nstarmax,nstartot,T2_star,n_star,eps_star,seed,m_star,sf_model
  ! Sink particles and black hole parameters
  namelist/sink_params/sink,nsinkmax,nsinktot,rho_type_sink,sink_descent,fudge_descent &
       & ,sink_relevance_threshold,sink_density_threshold,sink_saddle_threshold &
       & ,sink_mass_threshold,sink_purity_threshold,sink_fraction_threshold &
       & ,sink_nstar_frac,sink_mseed,sink_delta_tout &
       & ,sink_form,sink_merge,verbose_sink,sink_dump,drag_sink
  ! Black Hole accretion parameters
  namelist/sink_accretion_params/accretion_type,acc_sink_boost,bondi_use_vrel,use_rho_inf &
       & ,eddington_cap,sink_b_spline_order,bondi_use_gas_mass,use_bondi_lambda &
       & ,t_start_black_hole,use_local_bondi_rate,static_sink &
       & ,fix_sink_mass,eddington_floor,mass_weighting,momentum_conserving &
       & ,manual_accretion_rate
  ! AGN Feedback parameters
  namelist/sink_feedback_params/agn,agn_feedback_radius,agn_weighting_scheme,epsilon_rad &
       & ,epsilon_radio,epsilon_quasar,momentum_boost,agn_fbk_mode_switch_threshold &
       & ,agn_jet_opening_angle,agn_use_mass_weighting
  ! Supernovae feedback parameters
  namelist/feedback_params/M_SNII,E_SNII,t_SNII,eta_SNII,yield_SNII,thermal_feedback,mechanical_feedback
  ! Clump finder parameters
  namelist/clump_params/clump_finder,clump_info &
       & ,output_clump,output_peak_grid,output_peak_part,output_peak_star,output_peak_sink,output_peak_tree &
       & ,output_peak_trac,output_peak_dust &
       & ,relevance_threshold,density_threshold,saddle_threshold &
       & ,mass_threshold,purity_threshold,fraction_threshold &
       & ,merger_tree,orphan,ntreemax,ntreetot,rho_type_clump,nsteps_per_tree
  ! Lightcone parameters
  namelist/lightcone_params/lightcone,cone_z_min,cone_z_max,cone_opening_angle_y,cone_opening_angle_z &
       & ,cone_theta,cone_phi,cone_observer
  ! Gadget initial conditions parameters
  namelist/gadget_params/ic_file,ic_format,IG_rho,IG_T2,IG_metal &
       & ,ic_head_name,ic_pos_name,ic_vel_name,ic_id_name,ic_mass_name &
       & ,ic_u_name,ic_metal_name,ic_age_name &
       & ,gadget_scale_l, gadget_scale_v, gadget_scale_m ,gadget_scale_t &
       & ,ic_skip_type
  ! Turbulence driving parameters
  namelist/turb_params/turb, turb_seed, comp_frac,&
       & forcing_power_spectrum, turb_T, turb_Ndt, turb_rms, turb_min_rho

  associate(s=>pst%s)

  !--------------------------------------------------
  ! Advertise RAMSES
  !--------------------------------------------------
  write(*,*)'_/_/_/       _/_/     _/    _/    _/_/_/   _/_/_/_/    _/_/_/  '
  write(*,*)'_/    _/    _/  _/    _/_/_/_/   _/    _/  _/         _/    _/ '
  write(*,*)'_/    _/   _/    _/   _/ _/ _/   _/        _/         _/       '
  write(*,*)'_/_/_/     _/_/_/_/   _/    _/     _/_/    _/_/_/       _/_/   '
  write(*,*)'_/    _/   _/    _/   _/    _/         _/  _/               _/ '
  write(*,*)'_/    _/   _/    _/   _/    _/   _/    _/  _/         _/    _/ '
  write(*,*)'_/    _/   _/    _/   _/    _/    _/_/_/   _/_/_/_/    _/_/_/  '
  write(*,*)'                        Version 3.0                            '
  write(*,*)'       written by Romain Teyssier (Princeton University)       '
  write(*,*)'        (c) CEA 1999-2007, UZH 2008-2021, PU 2022-2026         '
  write(*,*)' '

  write(*,'(" Working with ndim = ",I0)')ndim
#ifdef GRAV
  write(*,'(" Using gravity solver")')
#endif
#ifdef HYDRO
  write(*,'(" Using hydro solver with nvar = ",I0)')nvar
  ! Check nvar is not too small
  if(nvar<5)then
     write(*,*)'You should have: nvar>=5'
     write(*,'(" Please recompile with -DNVAR=5")')
     call mdl_abort(s%mdl)
  endif
#endif
#ifdef DO_RT
  write(*,'(" Using radiation solver with nrtgrp = ",I0)')nrtgrp
#endif
#ifdef DO_CR
  write(*,'(" Using 2-moment cosmic rays solver with ncrgrp = ",I2)')ncrgrp
#endif

  ! Write information about git version
  call write_gitinfo

  ! Read namelist filename from command line argument
  narg = command_argument_count()
  IF(narg .LT. 1)THEN
     write(*,*)'You should type: ramses3d input.nml [nrestart]'
     write(*,*)'File input.nml should contain a parameter namelist'
     write(*,*)'nrestart is optional'
     call mdl_abort(s%mdl)
  END IF
  CALL getarg(1,infile)

  !-------------------------------------------------
  ! Read the namelist
  !-------------------------------------------------
  namelist_file=TRIM(infile)
  INQUIRE(file=infile,exist=nml_ok)
  if(.not. nml_ok)then
     write(*,*)'File '//TRIM(infile)//' does not exist'
     call mdl_abort(s%mdl)
  end if

  open(1,file=infile)
  rewind(1)
  read(1,NML=run_params)
  rewind(1)
  read(1,NML=output_params)
  rewind(1)
  read(1,NML=traj_params,END=83)
83 continue
  rewind(1)
  read(1,NML=amr_params)
  rewind(1)
  read(1,NML=movie_params,END=82)
82 continue
  rewind(1)
  read(1,NML=poisson_params,END=81)
81 continue

  !-------------------------------------------------
  ! Read optional nrestart command-line argument
  !-------------------------------------------------
  if (narg==2) then
    CALL getarg(2,cmdarg)
    read(cmdarg,*) nrestart
  endif

  !-------------------------------------------------
  ! Compute time step for outputs
  !-------------------------------------------------
  if(tend>0)then
     if(delta_tout==0)delta_tout=tend
     noutput=MIN(int(tend/delta_tout),MAXOUT)
     do i=1,noutput
        tout(i)=dble(i)*delta_tout
     end do
  else if(aend>0)then
     if(delta_aout==0)delta_aout=aend
     noutput=MIN(int(aend/delta_aout),MAXOUT)
     do i=1,noutput
        aout(i)=dble(i)*delta_aout
     end do
  endif
  noutput=max(count(aout>0.and.aout<=1),count(tout>0))
  noutput=MIN(noutput,MAXOUT)
  if(imovout>0) then
     if(tendmov==0.and.aendmov==0)movie=.false.
  endif
  if(nfile==-1)nfile=s%g%ncpu
  nfile=max(nfile,1)
  nfile=min(nfile,s%g%ncpu)

  !-------------------------------------------------
  ! Set default initialisation of ionisation states:
  ! -Off if restarting, but can set to true (for postprocessing)
  ! -On otherwise, but can set to false
  !-------------------------------------------------
  if(nrestart .gt. 0) then
     is_init_xion=.false.
  else
     is_init_xion=.true.
  endif
  !--------------------------------------------------
  ! Check for errors in the namelist so far
  !--------------------------------------------------
  levelmin=MAX(levelmin,1)
  nlevelmax=levelmax
  nsuperoct=MIN(nsuperoct,5)
  nml_ok=.true.
  if(levelmin<1)then
     write(*,*)'Error in the namelist:'
     write(*,*)'levelmin should not be lower than 1 !!!'
     nml_ok=.false.
  end if
  if(nlevelmax<levelmin)then
     write(*,*)'Error in the namelist:'
     write(*,*)'levelmax should not be lower than levelmin'
     nml_ok=.false.
  end if
  if(part_dep_algo<1 .or. part_dep_algo>3)then
     write(*,*)'Error in the namelist:'
     write(*,*)'part_dep_algo must be 1 (large), 2 (medium) or 3 (small)'
     nml_ok=.false.
  end if
  if(star_dep_algo<1 .or. star_dep_algo>3)then
     write(*,*)'Error in the namelist:'
     write(*,*)'star_dep_algo must be 1 (large), 2 (medium) or 3 (small)'
     nml_ok=.false.
  end if
  if(ngridmax==0)then
     if(ngridtot==0)then
        write(*,*)'Error in the namelist:'
        write(*,*)'Allocate some space for refinements !!!'
        nml_ok=.false.
     else
        ngridmax=int(ngridtot/int(s%g%ncpu,kind=8),kind=4)
     endif
  end if
  !--------------------------------------------------
  ! Compute maximum number of particles:
  ! dm, stars, sinks, trees, tracers, and dust
  !--------------------------------------------------

  if(npartmax==0)then
     npartmax=int(nparttot/int(s%g%ncpu,kind=8),kind=4)
  endif
  if(pic.and.npartmax>0)part=.true.
  if(nstarmax==0)then
     nstarmax=int(nstartot/int(s%g%ncpu,kind=8),kind=4)
     if(nstarmax==0)star=.false.
  endif
  if(nsinkmax==0)then
     nsinkmax=int(nsinktot/int(s%g%ncpu,kind=8),kind=4)
     if(nsinkmax==0)sink=.false.
  endif
  if(ntreemax==0)then
     ntreemax=int(ntreetot/int(s%g%ncpu,kind=8),kind=4)
     if(ntreemax==0)then
        ntreemax=int(npartmax/100)
     endif
     if(ntreemax==0)merger_tree=.false.
  endif
#ifdef HYDRO
  if(.not. hydro)then
     write(*,*)'You are not using the hydro solver but'
     write(*,*)'the code was compiled with HYDRO=1'
     write(*,*)'This might not be optimal but I am still running.'
  endif
#else
  if(hydro)then
     write(*,*)'You are using the hydro solver but'
     write(*,*)'the code was compiled with HYDRO=0'
     write(*,*)'Please recompile with HYDRO=1'
     call mdl_abort(s%mdl)
  endif
#endif
#ifdef GRAV
  if(.not. poisson)then
     write(*,*)'You are not using the poisson solver but'
     write(*,*)'the code was compiled with GRAV=1'
     write(*,*)'This might not be optimal but I am still running.'
  endif
#else
  if(poisson)then
     write(*,*)'You are using the poisson solver but'
     write(*,*)'the code was compiled with GRAV=0'
     write(*,*)'Please recompile with GRAV=1'
     call mdl_abort(s%mdl)
  endif
#endif
#ifdef DO_RT
  if(.not. rt)then
     write(*,*)'You are not using the rt solver but'
     write(*,*)'the code was compiled with RT=1'
     write(*,*)'This is just a warning and RAMSES will continue'
  endif
#else
  if(rt)then
     write(*,*)'You are trying to use the rt solver but'
     write(*,*)'the code was compiled with RT=0'
     write(*,*)'Please recompile with RT=1'
     call mdl_abort(s%mdl)
  endif
#endif
#ifdef DO_CR
  if(.not. cr)then
     write(*,*)'You are not using the cr solver but'
     write(*,*)'the code was compiled with CR=1'
     write(*,*)'This is just a warning and RAMSES will continue'
  endif
#else
  if(cr)then
     write(*,*)'You are trying to use the cr solver but'
     write(*,*)'the code was compiled with CR=0'
     write(*,*)'Please recompile with CR=1'
     call mdl_abort(s%mdl)
  endif
#endif
#ifndef MHD
  if(cr)then
     write(*,*)'You are trying to use the cr solver but'
     write(*,*)'the code was compiled with MHD=0'
     write(*,*)'Please recompile with MHD=1'
     call mdl_abort(s%mdl)
  endif
#endif

  !----------------------------
  ! Read hydro parameters
  !----------------------------
  rewind(1)
  read(1,NML=init_params,END=101)
  goto 102
101 write(*,*)' You need to set up namelist &INIT_PARAMS in parameter file'
  call mdl_abort(s%mdl)
102 rewind(1)
  if(nlevelmax>levelmin)read(1,NML=refine_params)
  rewind(1)
  if(hydro)read(1,NML=hydro_params)
  rewind(1)
  read(1,NML=units_params,END=105)
105 continue
  rewind(1)
  read(1,NML=boundary_params,END=106)
106 continue
  rewind(1)
  read(1,NML=cooling_params,END=107)
107 continue
  rewind(1)
  read(1,NML=star_params,END=108)
108 continue
  rewind(1)
  read(1,NML=feedback_params,END=109)
109 continue
  rewind(1)
  read(1,NML=clump_params,END=110)
110 continue
  rewind(1)
  read(1,NML=gadget_params,END=111)
111 continue
  rewind(1)
  read(1,NML=sink_params,END=112)
112 continue
  rewind(1)
  read(1,NML=sink_accretion_params,END=113)
113 continue
  rewind(1)
  read(1,NML=sink_feedback_params,END=114)
114 continue
  rewind(1)
  read(1, NML=lightcone_params, END=115)
115 continue
  rewind(1)
  read(1, NML=turb_params, END=116)
116 continue
  rewind(1)
  read(1,NML=trac_params,END=117)
117 continue
  rewind(1)
  read(1,NML=dust_params,END=118)
118 continue
  close(1)

  ! Compute maximum number of tracer particles
  if(ntracmax==0)then
     ntracmax=int(ntractot/int(s%g%ncpu,kind=8),kind=4)
     if(ntracmax==0)trac=.false.
  endif
  if(ndustmax==0)then
     ndustmax=int(ndusttot/int(s%g%ncpu,kind=8),kind=4)
     if(ndustmax==0)dust=.false.
  endif

  !-----------------
  ! Max size checks
  !-----------------
  if(nlevelmax>MAXLEVEL)then
     write(*,*) 'Error: nlevelmax>MAXLEVEL'
     call mdl_abort(s%mdl)
  end if
  if(nregion>MAXREGION)then
     write(*,*) 'Error: nregion>MAXREGION'
     call mdl_abort(s%mdl)
  end if

  !-----------------------------------
  ! Rearrange level dependent arrays
  !-----------------------------------
  do i=nlevelmax,levelmin,-1
     nexpand   (i)=nexpand   (i-levelmin+1)
     nsubcycle (i)=nsubcycle (i-levelmin+1)
     r_refine  (i)=r_refine  (i-levelmin+1)
     a_refine  (i)=a_refine  (i-levelmin+1)
     b_refine  (i)=b_refine  (i-levelmin+1)
     x_refine  (i)=x_refine  (i-levelmin+1)
     y_refine  (i)=y_refine  (i-levelmin+1)
     z_refine  (i)=z_refine  (i-levelmin+1)
     m_refine  (i)=m_refine  (i-levelmin+1)
     exp_refine(i)=exp_refine(i-levelmin+1)
     initfile  (i)=initfile  (i-levelmin+1)
     jeans_refine(i)=jeans_refine(i-levelmin+1)
  end do
  do i=1,levelmin-1
     nexpand   (i)= 1
     nsubcycle (i)= 1
     r_refine  (i)=-1.0
     a_refine  (i)= 1.0
     b_refine  (i)= 1.0
     x_refine  (i)= 0.0
     y_refine  (i)= 0.0
     z_refine  (i)= 0.0
     m_refine  (i)=-1.0
     exp_refine(i)= 2.0
     initfile  (i)= ' '
     jeans_refine(i)=-1.0
  end do

  !--------------------------------------------------
  ! Check for non-thermal energies
  !--------------------------------------------------
#if NENER>0
  if(nvar<5+nener)then
     write(*,*)'Error: non-thermal energy need nvar >= 5+nener'
     write(*,*)'Modify NENER and recompile'
     nml_ok=.false.
  endif
#endif

  !--------------------------------------------------
  ! Check for metals
  !--------------------------------------------------
  if(metal.and.nvar<6)then
     write(*,*)'Error: metal=.true. need nvar >= 6'
     write(*,*)'Modify NVAR and recompile'
     nml_ok=.false.
  endif

  !--------------------------------------------------
  ! Check for entropy
  !--------------------------------------------------
  if(entropy.and.nvar<6)then
     write(*,*)'Error: entropy=.true. need nvar >= 6'
     write(*,*)'Modify NVAR and recompile'
     nml_ok=.false.
  endif

  !--------------------------------------------------
  ! Compute indices for passive scalars
  ! and non-thermal energies
  !--------------------------------------------------
  inener=6
  ientropy=inener+nener
  imetal=ientropy
  if(entropy)imetal=ientropy+1
  iturb=imetal
  if(metal)then
     iturb=imetal+1
  endif
  ichem=iturb
  if(sgs_turb)then
     ichem=iturb+1
  endif

  !--------------------------------------------------
  ! Check for sgs_turb
  !--------------------------------------------------
  if(sgs_turb.and.iturb>nvar)then
     write(*,*)'Error: sgs_turb=.true. needs nvar >= ',iturb
     write(*,*)'Currently nvar=',nvar,' but iturb=',iturb
     write(*,*)'Modify NVAR and recompile'
     nml_ok=.false.
  endif
  if(hydro.and.(nvar>5)) then
     write(*,'(A50)')"__________________________________________________"
     write(*,*) 'Hydro var extra indices:'
#if NENER>0
                 write(*,*) '   inener   = ',inener
#endif
     if(entropy) write(*,*) '   ientropy = ',ientropy
     if(metal)   write(*,*) '   imetal   = ',imetal
     if(sgs_turb)write(*,*) '   iturb    = ',iturb
     if(ichem.LE.nvar)then
                 write(*,*) '   ichem    = ',ichem
     endif
     write(*,'(A50)')"__________________________________________________"
  endif

  ! Force non-equilibrium chemistry if using rt
  if(rt .and. cooling) then
    print*,'Using radiation hydrodynamics, so switching to non-equilibrum cooling'
    cooling=.false.
    neq_chem=.true.
  endif
  !----------------------------------------------------------------
  ! Set number of used ionisation fractions, indices of ionization
  ! fractions, and ionization energies, and check if we have enough
  ! ionization variables (NION)
  !----------------------------------------------------------------
  if(neq_chem) then
     iCount=0
#ifdef DO_RTZ
     do i_Element=1,n_elements ! loop over elements
        if (elements(i_Element)%atomic_number.gt.0) then
           do i_Iion=1,elements(i_Element)%n_ions ! loop over ions
              icount = icount + 1
              if (i_Iion.eq.1) then
                 element_first_idx(i_Element) = icount
              end if
           end do ! end loop over ions
        end if
     end do ! end loop over elements

     ! Deal with molecules separately
     if (elements(1)%atomic_number.gt.0 .and. isH2_rtz) then
        icount = icount + 1
        molecules_first_idx(1) = icount
     end if
#else
     ! HI fraction and ionization energy
     if(isH2) then
        iCount=iCount+1 ; ixHI=iCount ; ionEvs(ixHI)=ionEv_HI
     endif
     ! HII fraction and ionization energy
     iCount=iCount+1 ; ixHII=iCount ; ionEvs(ixHII)=ionEv_HII
     ! HeII and HeIII fractions and ionization energies
     if(isHe) then
        iCount=iCount+1 ; ixHeII=iCount ; ionEvs(ixHeII)=ionEv_HeII
        iCount=iCount+1 ; ixHeIII=iCount ; ionEvs(ixHeIII)=ionEv_HeIII
     endif
#endif
     ! Check that we have enough chemical species
     if(iCount .gt. NION) then
        write(*,*) 'Not enough variables for ionization fractions'
        write(*,*) 'Have NION=',NION
        write(*,*) 'STOPPING!'
        call mdl_abort(s%mdl)
     endif
     if(iCount .lt. NION) then
        write(*,*) 'Too many variables for ionization fractions'
        write(*,*) 'Have NION=',NION
        write(*,*) 'Need NION=',iCount
        write(*,*) 'Probably no harm, so still continuing...'
     endif
     ! Starting index of ionized species
     iIons = ichem
     ichem = ichem + icount
     ! Check we have enough passive scalar
     if(iIons+iCount-1 .gt. nvar) then
        write(*,*) 'Something wrong with NVAR.'
        write(*,*) 'Have NVAR=',nvar
        write(*,*) 'Should have NVAR=',iIons+iCount-1
        write(*,*) 'STOPPING!'
        call mdl_abort(s%mdl)
     endif
     ! Output indices for the user to check
     write(*,'(A39, I4)') 'The number of ionization fractions is:',iCount
     write(*,*) 'Their indices in U are:'
#ifdef DO_RTZ
     do i_Element=1,n_elements ! Loop over elements
        if (elements(i_Element)%atomic_number.gt.0) then
           write(*,'(A10, I4)') 'i' // trim(elements(i_Element)%element_name) // '=', iIons-1+element_first_idx(i_Element)
        end if
     end do

     ! Deal with molecules separately
     if (elements(1)%atomic_number.gt.0 .and. isH2_rtz) then
        write(*,'(A10, I4)') 'iH2 =', iIons-1+molecules_first_idx(1)
     end if
#else
     if(isH2) write(*,'(A10, I2)') '  iHI =    ', iIons-1+ixHI
     write(*,'(A10, I2)')          '  iHII =   ', iIons-1+ixHII
     if(isHe) write(*,'(A10, I2)') '  iHeII =  ', iIons-1+ixHeII
     if(isHe) write(*,'(A10, I2)') '  iHeIII = ', iIons-1+ixHeIII
     !if(isHe) print '(I3, A9)', iIons-1+ixHeIII, 'iHeIII'
#endif
  endif

  ! Check that lightcone parameters are consistent
  if (lightcone) then
     if (ndim /= 3) then
        write(*, *) 'Error: lightcone is only supported in 3D'
        nml_ok = .false.
     endif
     if (cone_z_min >= cone_z_max) then
        write(*, *) 'Error: cone_z_min >= cone_z_max'
        nml_ok = .false.
     endif
     if (cone_z_min < 0.0) then
        write(*, *) 'Error: cone_z_min must be greater than 0'
        nml_ok = .false.
     endif
     if (cone_opening_angle_y <= 0.0 .or. cone_opening_angle_z <= 0.0) then
        write(*, *) 'Error: cone opening angles must be positive'
        nml_ok = .false.
     endif
     if (cone_opening_angle_y >= 90.0 .or. cone_opening_angle_z >= 90.0) then
        write(*, *) 'Error: cone opening angles must be < 90 degrees'
        nml_ok = .false.
     endif
  endif

  if(.not. nml_ok)then
     write(*,*)'Too many errors in the namelist'
     write(*,*)'Aborting...'
     call mdl_abort(s%mdl)
  end if

  ! Fill in all run parameters in corresponding structure

  s%r%cosmo=cosmo
  s%r%pic=pic
  s%r%poisson=poisson
  s%r%hydro=hydro
  s%r%rt=rt
  s%r%cr=cr
  s%r%part=part
  s%r%star=star
  s%r%sink=sink
  s%r%trac=trac
  s%r%dust=dust
  s%r%tree=merger_tree
  s%r%orphan=orphan
  s%r%verbose=verbose
  s%r%debug=debug
  s%r%nrestart=nrestart
  s%r%ncontrol=ncontrol
  s%r%nstepmax=nstepmax
  s%r%nsubcycle=nsubcycle
  s%r%nremap=nremap
  s%r%static_mesh=static_mesh
  s%r%static_gas=static_gas
  s%r%geom=geom
  s%r%overload=overload
  s%r%nsuperoct=nsuperoct
  s%r%clump_only=clump_only

  s%r%noutput=noutput
  s%r%foutput=foutput
  s%r%aout=aout
  s%r%tout=tout
  s%r%output_mode=output_mode
  s%r%gadget_output=gadget_output
  s%r%run_time_hrs=run_time_hrs
  s%r%bkp_time_hrs=bkp_time_hrs
  s%r%bkp_last_min=bkp_last_min
  s%r%bkp_modulo=bkp_modulo
  s%r%nfile=nfile
  s%r%output_part=output_part
  s%r%output_grav=output_grav
  s%r%output_hydro=output_hydro
  s%r%output_amr=output_amr

  s%r%levelmin=levelmin
  s%r%nlevelmax=nlevelmax
  s%r%ngridmax=ngridmax
  s%r%ncachemax=ncachemax
  s%r%npartmax=npartmax
  s%r%nstarmax=nstarmax
  s%r%nsinkmax=nsinkmax
  s%r%ntreemax=ntreemax
  s%r%ntracmax=ntracmax
  s%r%ndustmax=ndustmax
  s%r%ntrac_per_cell=ntrac_per_cell
  s%r%ndust_per_cell=ndust_per_cell
  s%r%part_subcell_positions=part_subcell_positions
  s%r%tracer_kick_pdf=tracer_kick_pdf
  s%r%dust_to_gas_mass_ratio=dust_to_gas_mass_ratio
  s%r%grain_size_parameter=grain_size_parameter
  s%r%grain_charge_parameter=grain_charge_parameter
  s%r%dust_gyro_factor=dust_gyro_factor
  s%r%analytic_dust_force=analytic_dust_force
  s%r%nexpand=nexpand
  s%r%boxlen=boxlen

  s%r%gravity_test=gravity_test
  s%r%epsilon=epsilon
  s%r%nvcycle=nvcycle
  s%r%gravity_type=gravity_type
  s%r%gravity_params=gravity_params
  s%r%cic_levelmax=cic_levelmax
  s%r%cg_levelmin=cg_levelmin
  s%r%fast_solver=fast_solver
  s%r%part_mass_deposition_scheme=part_mass_deposition_scheme
  s%r%part_dep_algo=part_dep_algo
  s%r%part_force_interpolation_scheme=part_force_interpolation_scheme
  s%r%star_mass_deposition_scheme=star_mass_deposition_scheme
  s%r%star_dep_algo=star_dep_algo
  s%r%star_force_interpolation_scheme=star_force_interpolation_scheme
  s%r%sink_mass_deposition_scheme=sink_mass_deposition_scheme
  s%r%sink_dep_algo=sink_dep_algo
  s%r%sink_force_interpolation_scheme=sink_force_interpolation_scheme
  s%r%tree_mass_deposition_scheme=tree_mass_deposition_scheme
  s%r%tree_force_interpolation_scheme=tree_force_interpolation_scheme
  s%r%trac_interpolation_scheme=trac_interpolation_scheme
  s%r%dust_mass_deposition_scheme=dust_mass_deposition_scheme
  s%r%dust_force_interpolation_scheme=dust_force_interpolation_scheme

  s%r%nw_frame=nw_frame
  s%r%nh_frame=nh_frame
  s%r%levelmax_frame=levelmax_frame
  s%r%ivar_frame=ivar_frame
  s%r%xcentre_frame=xcentre_frame
  s%r%ycentre_frame=ycentre_frame
  s%r%zcentre_frame=zcentre_frame
  s%r%deltax_frame=deltax_frame
  s%r%deltay_frame=deltay_frame
  s%r%deltaz_frame=deltaz_frame
  s%r%movie=movie
  s%r%zoom_only=zoom_only
  s%r%zoom_only_frame=zoom_only_frame
  s%r%imovout=imovout
  s%r%imov=imov
  s%r%tstartmov=tstartmov
  s%r%astartmov=astartmov
  s%r%tendmov=tendmov
  s%r%aendmov=aendmov
  s%r%proj_axis=proj_axis
  s%r%movie_vars_txt=movie_vars_txt
  s%r%theta_camera=theta_camera
  s%r%phi_camera=phi_camera
  s%r%dtheta_camera=dtheta_camera
  s%r%dphi_camera=dphi_camera
  s%r%tstart_theta_camera=tstart_theta_camera
  s%r%tstart_phi_camera=tstart_phi_camera
  s%r%tend_theta_camera=tend_theta_camera
  s%r%tend_phi_camera=tend_phi_camera
  s%r%focal_camera=focal_camera
  s%r%dist_camera=dist_camera
  s%r%ddist_camera=ddist_camera
  s%r%perspective_camera=perspective_camera
  s%r%smooth_frame=smooth_frame
  s%r%shader_frame=shader_frame
  s%r%method_frame=method_frame
  s%r%varmin_frame=varmin_frame
  s%r%varmax_frame=varmax_frame
  if(s%r%movie)call set_movie_vars(s%r)

  ! Trajectory output params
  s%r%ntrajectories=ntrajectories
  s%r%trajectories=trajectories

  s%r%gamma=gamma
  s%r%courant_factor=courant_factor
  s%r%smallc=smallc
  s%r%smallr=smallr
  s%r%niter_riemann=niter_riemann
  s%r%slope_type=slope_type
  if(slope_mag_type<0)then
     s%r%slope_mag_type=slope_type
  else
     s%r%slope_mag_type=slope_mag_type
  endif
  s%r%difmag=difmag
  s%r%etamag=etamag
  s%r%gamma_rad=gamma_rad(1:nener+1)
  s%r%dual_energy=dual_energy
  s%r%T2_fix=T2_fix
  s%r%induction=induction
  s%r%entropy=entropy
  s%r%sgs_turb=sgs_turb
  s%r%equilibrium_sgs=equilibrium_sgs
  s%r%inener=inener
  s%r%ientropy=ientropy
  s%r%imetal=imetal
  s%r%ichem=ichem
  s%r%iturb=iturb
  s%r%scheme=scheme
  s%r%constant_gravity=constant_gravity
#ifndef MHD
  if(riemann=='llf')s%r%riemann=solver_llf
  if(riemann=='hll')s%r%riemann=solver_hll
  if(riemann=='hllc')s%r%riemann=solver_hllc
#endif
#ifdef MHD
  if(riemann=='llf')s%r%riemann=solver_llf
  if(riemann=='hll')s%r%riemann=solver_hll
  if(riemann=='hlld')s%r%riemann=solver_hlld
  if(riemann=='roe')s%r%riemann=solver_roe
  if(riemann=='upwind')s%r%riemann=solver_upwind
#ifdef _METAL
  if(riemann=='uct-hlld')s%r%riemann=solver_uct_hlld
#endif
#endif
#ifdef MHD
  if(riemann2d=='none')riemann2d=riemann
  if(riemann2d=='llf')s%r%riemann2d=solver2d_llf
  if(riemann2d=='hll')s%r%riemann2d=solver2d_hllf
  if(riemann2d=='hllf')s%r%riemann2d=solver2d_hllf
  if(riemann2d=='hlla')s%r%riemann2d=solver2d_hlla
  if(riemann2d=='hlld')s%r%riemann2d=solver2d_hlld
  if(riemann2d=='roe')s%r%riemann2d=solver2d_roe
  if(riemann2d=='upwind')s%r%riemann2d=solver2d_upwind
#endif
  s%r%switch_llf_dmin=switch_llf_dmin
  s%r%switch_llf_pmin=switch_llf_pmin
  s%r%smagorinsky_lilly_constant=smagorinsky_lilly_constant

  s%r%units_density=units_density
  s%r%units_time=units_time
  s%r%units_length=units_length
  s%r%units_np=units_np

  s%r%m_refine=m_refine
  s%r%r_refine=r_refine
  s%r%x_refine=x_refine
  s%r%y_refine=y_refine
  s%r%z_refine=z_refine
  s%r%exp_refine=exp_refine
  s%r%a_refine=a_refine
  s%r%b_refine=b_refine
  s%r%jeans_refine=jeans_refine
  s%r%var_cut_refine=var_cut_refine
  s%r%mass_cut_refine=mass_cut_refine
  s%r%ivar_refine=ivar_refine
  s%r%aexp_lock_refine=aexp_lock_refine
  s%r%pic_lock_refine=pic_lock_refine

  s%r%interpol_var=interpol_var
  s%r%interpol_type=interpol_type
  s%r%err_grad_d=err_grad_d
  s%r%err_grad_u=err_grad_u
  s%r%err_grad_p=err_grad_p
  s%r%floor_d=floor_d
  s%r%floor_u=floor_u
  s%r%floor_p=floor_p
  s%r%mass_sph=mass_sph
#ifdef MHD
  s%r%err_grad_b2=err_grad_b2
  s%r%err_grad_A=err_grad_A
  s%r%err_grad_B=err_grad_B
  s%r%err_grad_C=err_grad_C
  s%r%floor_b2=floor_b2
  s%r%floor_A=floor_A
  s%r%floor_B=floor_B
  s%r%floor_C=floor_C
#endif
#if NENER>0
  s%r%err_grad_prad=err_grad_prad
  s%r%floor_prad=floor_prad
#endif
#if NVAR>5+NENER
  s%r%err_grad_var=err_grad_var
#endif
  s%r%err_grad_xHI=err_grad_xHI
  s%r%err_grad_xHII=err_grad_xHII
  s%r%floor_xHI=floor_xHI
  s%r%floor_xHII=floor_xHII
#ifdef DO_RT
  s%r%rt_err_grad_cn=rt_err_grad_cn
  s%r%rt_floor_cn=rt_floor_cn
  s%r%rt_refine_aexp=rt_refine_aexp
#endif

  if(nrestart>0)filetype='restart'
  s%r%filetype=filetype
  s%r%initfile=initfile
  s%r%multiple=multiple
  s%r%aexp_ini=aexp_ini
  s%r%omega_b=omega_b
  s%r%boxlen_ini=boxlen_ini
  s%r%omega_m=omega_m
  s%r%omega_l=omega_l
  s%r%h0=h0
  s%r%ic_scale_m=ic_scale_m

  s%r%nregion=nregion
  s%r%region_type=region_type
  s%r%x_center=x_center
  s%r%y_center=y_center
  s%r%z_center=z_center
  s%r%length_x=length_x
  s%r%length_y=length_y
  s%r%length_z=length_z
  s%r%exp_region=exp_region
  s%r%d_region=d_region
  s%r%u_region=u_region
  s%r%v_region=v_region
  s%r%w_region=w_region
  s%r%p_region=p_region
#if NENER>0
  s%r%prad_region=prad_region
#endif
#if NVAR>5+NENER
  s%r%var_region=var_region
#endif
#ifdef MHD
  s%r%B_region=B_region
  s%r%C_region=C_region
  s%r%A_ave=A_ave
  s%r%B_ave=B_ave
  s%r%C_ave=C_ave
#endif
#ifdef DO_RT
  s%r%rt_nregion=rt_nregion
  s%r%rt_region_type=rt_region_type
  s%r%rt_reg_x_center=rt_reg_x_center
  s%r%rt_reg_y_center=rt_reg_y_center
  s%r%rt_reg_z_center=rt_reg_z_center
  s%r%rt_reg_length_x=rt_reg_length_x
  s%r%rt_reg_length_y=rt_reg_length_y
  s%r%rt_reg_length_z=rt_reg_length_z
  s%r%rt_exp_region=rt_exp_region
  s%r%rt_reg_group=rt_reg_group
  s%r%rt_n_region=rt_n_region
  s%r%rt_u_region=rt_u_region
  s%r%rt_v_region=rt_v_region
  s%r%rt_w_region=rt_w_region
#endif
#ifdef DO_CR
  s%r%cr_nregion=cr_nregion
  s%r%cr_region_type=cr_region_type
  s%r%cr_reg_x_center=cr_reg_x_center
  s%r%cr_reg_y_center=cr_reg_y_center
  s%r%cr_reg_z_center=cr_reg_z_center
  s%r%cr_reg_length_x=cr_reg_length_x
  s%r%cr_reg_length_y=cr_reg_length_y
  s%r%cr_reg_length_z=cr_reg_length_z
  s%r%cr_exp_region=cr_exp_region
  s%r%cr_reg_group=cr_reg_group
  s%r%cr_fx_region=cr_fx_region
  s%r%cr_fy_region=cr_fy_region
  s%r%cr_fz_region=cr_fz_region
#endif

  s%r%periodic=periodic
  s%r%nbound=nbound
  s%r%no_inflow=no_inflow
  s%r%bound_levelmin=bound_levelmin
  s%r%box_size=box_size
  s%r%box_xmin=box_xmin
  s%r%box_xmax=box_xmax
  s%r%box_ymin=box_ymin
  s%r%box_ymax=box_ymax
  s%r%box_zmin=box_zmin
  s%r%box_zmax=box_zmax
  s%r%bound_dir=bound_dir
  s%r%bound_type=bound_type
  s%r%bound_shift=bound_shift
  s%r%bound_xmin=bound_xmin
  s%r%bound_xmax=bound_xmax
  s%r%bound_ymin=bound_ymin
  s%r%bound_ymax=bound_ymax
  s%r%bound_zmin=bound_zmin
  s%r%bound_zmax=bound_zmax
  s%r%d_bound=d_bound
  s%r%u_bound=u_bound
  s%r%v_bound=v_bound
  s%r%w_bound=w_bound
  s%r%p_bound=p_bound
#if NENER>0
  s%r%prad_bound=prad_bound
#endif
#if NVAR>5+NENER
  s%r%var_bound=var_bound
#endif
#ifdef DO_RT
  s%r%rt_n_bound=rt_n_bound
  s%r%rt_u_bound=rt_u_bound
  s%r%rt_v_bound=rt_v_bound
  s%r%rt_w_bound=rt_w_bound
#endif
#ifdef DO_CR
  s%r%cr_fx_bound=cr_fx_bound
  s%r%cr_fy_bound=cr_fy_bound
  s%r%cr_fz_bound=cr_fz_bound
#endif

  s%r%cooling=cooling
  s%r%cooling_ism=cooling_ism
  s%r%metal=metal
  s%r%isothermal=isothermal
  s%r%haardt_madau=haardt_madau
  s%r%self_shielding=self_shielding
  s%r%J21=J21
  s%r%a_spec=a_spec
  s%r%z_ave=z_ave
  s%r%z_reion=z_reion
  s%r%cooling_uvb_delta=cooling_uvb_delta
  s%r%eos_type=eos_type
  s%r%eos_nH=eos_nH
  s%r%eos_index=eos_index
  s%r%eos_T2=eos_T2
  s%r%T2max=T2max
  s%r%mu_mol=mu_mol
  s%r%X_H=X_H
  s%r%Y_He=Y_He
  s%r%neq_chem=neq_chem
  s%r%is_init_xion=is_init_xion
  s%r%isHe=isHe
  s%r%isH2=isH2
  s%r%neq_Tconst=neq_Tconst
  if(neq_Tconst .ge. 0d0) s%r%neq_isTconst=.true.
  s%r%upload_equilibrium_x = upload_equilibrium_x

  s%r%rtz_cooling=rtz_cooling
  s%r%rtz_equilibrium_test=rtz_equilibrium_test
  s%r%rtz_include_collisional_ionization=rtz_include_collisional_ionization
  s%r%rtz_include_photoionization=rtz_include_photoionization
  s%r%rtz_include_cosmic_ray_ionization=rtz_include_cosmic_ray_ionization
  s%r%rtz_include_charge_exchange=rtz_include_charge_exchange
  s%r%rtz_include_dust_recombination=rtz_include_dust_recombination
  s%r%rtz_include_HM12_UVB=rtz_include_HM12_UVB
  s%r%isH2_rtz=isH2_rtz
  s%r%rtz_UV_background_G0=rtz_UV_background_G0
  s%r%rtz_primary_cosmic_ray_ionization_rate=rtz_primary_cosmic_ray_ionization_rate
  s%r%rtz_max_cool_timestep=rtz_max_cool_timestep
  s%r%rtz_eqm_min_its=rtz_eqm_min_its
#ifdef DO_RTZ
  s%r%element_first_idx=element_first_idx
  s%r%molecules_first_idx=molecules_first_idx
#endif

  s%r%iIons=iIons
  s%r%ixHI=ixHI
  s%r%ixHII=ixHII
  s%r%ixHeII=ixHeII
  s%r%ixHeIII=ixHeIII
  s%r%ionEvs=ionEvs

  s%r%T2_star=T2_star
  s%r%n_star=n_star
  s%r%eps_star=eps_star
  s%r%seed=seed
  s%r%m_star=m_star
  s%r%sf_model=sf_model

  s%r%M_SNII=M_SNII
  s%r%E_SNII=E_SNII
  s%r%t_SNII=t_SNII
  s%r%eta_SNII=eta_SNII
  s%r%yield_SNII=yield_SNII
  s%r%thermal_feedback=thermal_feedback
  s%r%mechanical_feedback=mechanical_feedback

  s%r%clump_finder=clump_finder
  s%r%clump_info=clump_info
  s%r%output_clump=output_clump
  s%r%output_peak_grid=output_peak_grid
  s%r%output_peak_part=output_peak_part
  s%r%output_peak_star=output_peak_star
  s%r%output_peak_sink=output_peak_sink
  s%r%output_peak_tree=output_peak_tree
  s%r%output_peak_trac=output_peak_trac
  s%r%output_peak_dust=output_peak_dust
  s%r%relevance_threshold=relevance_threshold
  s%r%density_threshold=density_threshold
  s%r%saddle_threshold=saddle_threshold
  s%r%mass_threshold=mass_threshold
  s%r%purity_threshold=purity_threshold
  s%r%fraction_threshold=fraction_threshold
  s%r%rho_type_clump=rho_type_clump
  s%r%nsteps_per_tree=nsteps_per_tree

  s%r%lightcone = lightcone
  s%r%cone_z_min = cone_z_min
  s%r%cone_z_max = cone_z_max
  s%r%cone_opening_angle_y = cone_opening_angle_y
  s%r%cone_opening_angle_z = cone_opening_angle_z
  s%r%cone_theta = cone_theta
  s%r%cone_phi = cone_phi
  s%r%cone_observer = cone_observer

  s%r%rho_type_sink=rho_type_sink
  s%r%sink_descent=sink_descent
  s%r%fudge_descent=fudge_descent
  s%r%sink_relevance_threshold=sink_relevance_threshold
  s%r%sink_density_threshold=sink_density_threshold
  s%r%sink_saddle_threshold=sink_saddle_threshold
  s%r%sink_mass_threshold=sink_mass_threshold
  s%r%sink_purity_threshold=sink_purity_threshold
  s%r%sink_fraction_threshold=sink_fraction_threshold
  s%r%static_sink=static_sink
  s%r%sink_delta_tout=sink_delta_tout
  s%r%sink_nstar_frac=sink_nstar_frac
  s%r%sink_mseed=sink_mseed
  s%r%fix_sink_mass=fix_sink_mass
  s%r%drag_sink=drag_sink

  s%r%accretion_type = accretion_type
  s%r%acc_sink_boost = acc_sink_boost
  s%r%bondi_use_vrel = bondi_use_vrel
  s%r%eddington_cap = eddington_cap
  s%r%sink_form = sink_form
  s%r%sink_merge = sink_merge
  s%r%sink_refine = sink_refine
  s%r%sink_dump = sink_dump
  s%r%sink_b_spline_order = sink_b_spline_order
  s%r%verbose_sink = verbose_sink
  s%r%bondi_use_gas_mass = bondi_use_gas_mass
  s%r%use_local_bondi_rate = use_local_bondi_rate
  s%r%use_rho_inf = use_rho_inf
  s%r%t_start_black_hole = t_start_black_hole
  s%r%use_bondi_lambda = use_bondi_lambda
  s%r%mass_weighting = mass_weighting
  s%r%momentum_conserving = momentum_conserving
  s%r%eddington_floor = eddington_floor
  s%r%manual_accretion_rate = manual_accretion_rate

  s%r%agn = agn
  s%r%agn_feedback_radius = agn_feedback_radius
  s%r%agn_weighting_scheme = agn_weighting_scheme
  s%r%epsilon_rad = epsilon_rad
  s%r%epsilon_radio = epsilon_radio
  s%r%epsilon_quasar = epsilon_quasar
  s%r%momentum_boost = momentum_boost
  s%r%agn_fbk_mode_switch_threshold = agn_fbk_mode_switch_threshold
  s%r%agn_jet_opening_angle = agn_jet_opening_angle
  s%r%agn_use_mass_weighting = agn_use_mass_weighting

  s%r%ic_file=ic_file
  s%r%ic_format=ic_format
  s%r%IG_rho=IG_rho
  s%r%IG_T2=IG_T2
  s%r%IG_metal=IG_metal
  s%r%ic_head_name=ic_head_name
  s%r%ic_pos_name=ic_pos_name
  s%r%ic_vel_name=ic_vel_name
  s%r%ic_id_name=ic_id_name
  s%r%ic_mass_name=ic_mass_name
  s%r%ic_u_name=ic_u_name
  s%r%ic_metal_name=ic_metal_name
  s%r%ic_age_name=ic_age_name
  s%r%gadget_scale_l=gadget_scale_l
  s%r%gadget_scale_v=gadget_scale_v
  s%r%gadget_scale_m=gadget_scale_m
  s%r%gadget_scale_t=gadget_scale_t
  s%r%ic_skip_type=ic_skip_type

  s%r%turb=turb
  s%r%turb_seed=turb_seed
  s%r%forcing_power_spectrum=forcing_power_spectrum
  s%r%comp_frac=comp_frac
  s%r%turb_T=turb_T
  s%r%turb_Ndt=turb_Ndt
  s%r%turb_rms=turb_rms
  s%r%turb_min_rho=turb_min_rho

  ! Read RT parameters from namelist
  if(rt)call m_read_rt_params(pst)

  ! Read CR parameters from namelist
  if(cr)call m_read_cr_params(pst)

  ! Broadcast parameters to all CPUs.
  call m_broadcast_params(pst)

  end associate

end subroutine m_read_params
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_broadcast_params(pst)
  use ramses_commons, only: pst_t
  implicit none
  type(pst_t)::pst
  !--------------------------------------------------------------------
  ! This routine is the master procedure to broadcast the run
  ! parameters to all the CPUs.
  !--------------------------------------------------------------------

  ! Broadcast parameters to all CPUs.
  call r_broadcast_params(pst,pst%s%r,storage_size(pst%s%r)/32)

end subroutine m_broadcast_params
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
recursive subroutine r_broadcast_params(pst,input,input_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use amr_commons, only: run_t, set_hydro_parameters
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  type(run_t)::input

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_BCAST_PARAMS,pst%iUpper+1,input_size,0,input)
     call r_broadcast_params(pst%pLower,input,input_size)
     call mdl_get_reply(pst%s%mdl,rID,0)
  else
     pst%s%r=input
     call set_hydro_parameters(pst%s%r,pst%s%h_params)
  endif

end subroutine r_broadcast_params
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine m_broadcast_global(pst)
  use ramses_commons, only: pst_t
  implicit none
  type(pst_t)::pst
  !--------------------------------------------------------------------
  ! This routine is the master procedure to broadcast the run
  ! parameters to all the CPUs.
  !--------------------------------------------------------------------
  integer::input_size
  integer,dimension(1:1)::dummy
  integer,dimension(:),allocatable::input_array

  ! Broadcast parameters to all CPUs.
  input_size=storage_size(pst%s%g)/32
  allocate(input_array(1:input_size))
  input_array=transfer(pst%s%g,input_array)
  call r_broadcast_global(pst,input_array,input_size,dummy,0)
  deallocate(input_array)

end subroutine m_broadcast_global
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
recursive subroutine r_broadcast_global(pst,input_array,input_size,output_array,output_size)
  use mdl_module
  use ramses_commons, only: pst_t
  use mdl_parameters
  implicit none
  type(pst_t)::pst
  integer,VALUE::input_size
  integer::output_size
  integer,dimension(1:input_size)::input_array
  integer,dimension(1:output_size)::output_array

  integer::rID

  if(pst%nLower>0)then
     rID = mdl_send_request(pst%s%mdl,MDL_BCAST_GLOBAL,pst%iUpper+1,input_size,output_size,input_array)
     call r_broadcast_global(pst%pLower,input_array,input_size,output_array,output_size)
     call mdl_get_reply(pst%s%mdl,rID,output_size)
  else
     pst%s%g=transfer(input_array,pst%s%g)
     pst%s%g%myid=mdl_self(pst%s%mdl)
  endif

end subroutine r_broadcast_global
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
end module params_module
