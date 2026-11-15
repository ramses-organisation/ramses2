!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine mdl_init
  use call_back
  use clean_stop_module, only: r_clean_stop
  use amr_parameters, only: flen
  use mdl_module
  use ramses_commons, only: pst_t, ramses_t
#ifdef DO_RTZ
  use rtz_module, only: initialize_elements
#endif
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_FUNLOC, C_SIZEOF
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  type(mdl_t),pointer::mdl
  type(pst_t),allocatable::pst

  allocate(mdl)
  call mdl_initialize(mdl)

  pst = worker_init(mdl)

#ifndef WITHOUTMPI
  ! Allocate input and output buffer sizes
  allocate(mdl%mpi_input_buffer(1:32+mdl%MDL_INPUT_MAXSIZE))

  ! Initialize software cache
  call init_cache(mdl)
#endif

#ifdef DO_RTZ
  call initialize_elements()
#endif

  ! For slave workers, go into waiting loop
  if(mdl_self(mdl)>1)then
     call mdl_wait(pst)
  else
     call master(pst)
  endif

#ifndef WITHOUTMPI
  call MPI_FINALIZE(info)
#endif

contains
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine master(pst)
  use init_amr_module, only: r_set_add
  implicit none
  type(pst_t)::pst
  call r_set_add(pst,mdl_threads(pst%s%mdl),1)
  call adaptive_loop(pst)
  call r_clean_stop(pst)
end subroutine master
!##############################################################
!##############################################################
!##############################################################
!##############################################################
function worker_init(mdl) result(pst)
  use ramses_commons, only: pst_t, ramses_t
  use call_back, only: ramses_function
  use mdl_parameters
  use flag_utils, only: r_init_flag, r_ensure_ref_rules, r_ensure_subgrid, r_user_flag
  use init_amr_module, only: r_init_amr, r_set_add
  use params_module, only: r_broadcast_params,r_broadcast_global
  use init_time_module, only: r_init_time
  use init_hydro_module, only: r_init_hydro
  use init_part_module, only: r_init_part, r_deallocate_gas
  use input_part_grafic_module, only: r_input_part_grafic, r_input_trac_grafic, r_input_dust_grafic
  use input_part_zoom_module, only: r_input_part_zoom
  use input_part_ascii_module, only: r_input_part_ascii, r_input_star_ascii, r_input_sink_ascii
  use input_part_restart_module, only: r_input_part_restart
  use input_part_ramses_module, only: r_input_part_ramses
  use input_part_gadget_module, only: r_input_part_gadget
  use input_part_module, only: r_npart_max, r_mass_min_part, r_broadcast_mp_min, r_check_part_emission
  use update_time_module, only: r_broadcast_aexp, r_hash_stats
  use init_refine_basegrid_module, only:r_init_refine_basegrid,r_collect_noct,r_noct_tot,r_noct_min,r_noct_max,r_noct_used_max
  use init_refine_restart_module, only: r_init_refine_restart
  use init_refine_ramses_module, only: r_init_refine_ramses
  use load_balance_module, only: r_load_balance,r_balance_part,r_broadcast_bound_key,r_collect_bound_key
  use refine_utils, only: r_refine_fine, r_clean_dirty
  use smooth_module, only: r_smooth_fine
  use input_hydro_condinit_module, only: r_input_hydro_condinit
  use input_hydro_grafic_module, only: r_input_hydro_grafic, r_input_refmap_grafic
  use input_hydro_gadget_module, only: r_input_hydro_gadget
  use upload_module, only: r_upload_fine
  use move_fine_module, only: r_kick_drift_part
  use output_amr_module, only: r_output_amr
  use output_hydro_module, only: r_output_hydro
  use output_poisson_module, only: r_output_poisson
  use output_part_module, only: r_output_part
  use synchro_hydro_fine_module, only: r_synchro_hydro_fine, r_gravity_hydro_fine
  use source_hydro_fine_module, only: r_source_hydro_fine
  use interpol_phi_module, only: r_save_phi_old
  use courant_fine_module, only: r_courant_fine
  use godunov_fine_module, only: r_godunov_fine,r_set_unew,r_set_uold,r_nimhd_fine
  use cooling_fine_module, only: r_cooling_fine
  use init_xion_module, only: r_init_xion
  use star_formation_module, only: r_star_formation
  use core_to_star_formation_module, only: r_core_to_star_formation
  use sink_formation_module, only: r_sink_formation,r_sink_clump
  use core_formation_module, only: r_core_formation,r_core_clump
  use input_part_ascii_module, only: r_input_core_ascii
  use tree_formation_module, only: r_tree_formation,r_tree_clump
  use feedback_module, only: r_thermal_feedback, r_mechanical_feedback
  use sink_evolution_module, only: r_sink_evolution
  use sink_mechanical_module, only: r_sink_mechanical
  use core_evolution_module, only: r_core_evolution,in_core_acc_t,out_core_acc_t
  use core_retirement_module, only: r_core_retirement, &
       in_core_retire_t,out_core_retire_t
  use sink_merger_module, only: r_sink_merger
  use core_merger_module, only: r_core_merger
  use newdt_fine_module, only: r_newdt_part,r_broadcast_dt,r_max_B_and_Q
  use rho_fine_module, only: r_split_part,r_sort_part
#ifdef GRAV
  use force_fine_module, only: r_force_analytic,r_compute_epot,r_compute_rhomax,r_gradient_phi
  use rho_fine_module, only: r_multipole_leaf_cells,r_multipole_split_cells,r_broadcast_multipole,r_collect_multipole,&
                            r_cic_multipole,r_cic_part,r_reset_rho
  use phi_fine_cg_module, only: r_cmp_pAp_cg,r_cmp_r2_cg,r_cmp_residual_cg,r_cmp_rhs_norm,&
                                r_make_initial_phi,r_recurrence_on_p,r_recurrence_x_and_r
  use multigrid_fine_commons, only: r_init_mg,r_build_mg,r_cleanup_mg,r_make_mask,r_make_bc_rhs,&
                                r_cmp_rhs_norm2
  use multigrid_fine_coarse, only: r_restrict_mask,r_cmp_residual_mg,r_cmp_residual_norm2,r_restrict_residual,&
                                r_reset_correction,r_set_scan_flag,r_gauss_seidel_mg,r_interpolate_and_correct
#endif
  use clump_finder_module, only: r_clump_finder
  use clump_merger_module, only: r_deallocate_clump
  use output_clump_module, only: r_output_clump
  use core_clump_output_module, only: r_output_core_clump
  use lightcone_module, only: r_output_lightcone
  use movie_module, only: r_output_frame
  use amr_parameters, only: nhilbert
  use init_rt_module, only: r_init_rt
  use rt_input_condinit_module, only: r_rt_input_condinit, r_rt_input_source_regions
  use rt_upload_module, only: r_rt_upload_fine
  use output_rt_module, only: r_output_rt
  use rt_godunov_fine_module, only: r_rt_godunov_fine,r_set_rtunew,r_set_rtuold,r_set_emissivity
  use update_rt_c_module, only: r_rt_neq_updates
  use rt_star_feedback, only: r_star_rt_feedback
  use init_cr_module, only: r_init_cr
  use cr_godunov_fine_module, only: r_cr_godunov_fine,r_set_crunew,r_set_cruold,r_conserve_cr_flux
  use cr_source_terms_module, only: r_cr_source_terms
  use update_cr_c_module, only: r_cr_updates, r_broadcast_cr_c
  use cr_input_condinit_module, only: r_cr_input_condinit
  use cr_upload_module, only: r_cr_upload_fine
  use output_cr_module, only: r_output_cr
#if defined(_CUDA)
  use gpu_manager, only: r_set_grid_device, r_transfer_grid_host
#elif defined(_METAL)
  use metal_runner, only: r_set_grid_device, r_transfer_grid_host
#endif
  use turb_driving, only: r_drive_turb
  use turb_hydro_module, only: r_turb_hydro
  use turb_init_module, only: r_init_turb
  use turb_update_module, only: r_update_turb

  implicit none

  integer::dummy
  integer(kind=8)::dummy8
  type(in_core_acc_t)::dummy_core_acc_in
  type(out_core_acc_t)::dummy_core_acc_out
  type(in_core_retire_t)::dummy_core_retire_in
  type(out_core_retire_t)::dummy_core_retire_out
  integer::ncpu
  type(mdl_t),target::mdl
  type(pst_t),allocatable::pst

  allocate(pst)
  allocate(pst%s)
  pst%s%mdl => mdl
  ncpu = mdl_threads(pst%s%mdl)

  ! Store cpu info as a global variable
  pst%s%g%myid=mdl_self(pst%s%mdl)
  pst%s%g%ncpu=mdl_threads(pst%s%mdl)

  ! Input and output sizes are in units of integer (single precision)
  call mdl_add_service(pst%s%mdl,MDL_CLEAN_STOP,             pst,C_FUNLOC(r_clean_stop),0,0,"clean_stop")
  call mdl_add_service(pst%s%mdl,MDL_SET_ADD,                pst,C_FUNLOC(r_set_add),storage_size(dummy)/32,0, "set_add")
  call mdl_add_service(pst%s%mdl,MDL_BCAST_PARAMS,           pst,C_FUNLOC(r_broadcast_params),storage_size(pst%s%r)/32,0,"broadcast_params")
  call mdl_add_service(pst%s%mdl,MDL_HASH_STATS,             pst,C_FUNLOC(r_hash_stats),0,0,"hash_stats")
  call mdl_add_service(pst%s%mdl,MDL_BCAST_GLOBAL,           pst,C_FUNLOC(r_broadcast_global),storage_size(pst%s%g)/32,0,"broadcast_global")
  call mdl_add_service(pst%s%mdl,MDL_INIT_AMR,               pst,C_FUNLOC(r_init_amr),0,0,"init_amr")
  call mdl_add_service(pst%s%mdl,MDL_INIT_TIME,              pst,C_FUNLOC(r_init_time),0,0,"init_time")
  call mdl_add_service(pst%s%mdl,MDL_INIT_HYDRO,             pst,C_FUNLOC(r_init_hydro),0,0,"init_hydro")
  call mdl_add_service(pst%s%mdl,MDL_INIT_PART,              pst,C_FUNLOC(r_init_part),0,0,"init_part")
  call mdl_add_service(pst%s%mdl,MDL_INPUT_PART_GRAFIC,      pst,C_FUNLOC(r_input_part_grafic),storage_size(pst%s%p%npart_tot)/32,0,"input_part_grafic")
  call mdl_add_service(pst%s%mdl,MDL_INPUT_TRAC_GRAFIC,      pst,C_FUNLOC(r_input_trac_grafic),storage_size(pst%s%p%npart_tot)/32,0,"input_trac_grafic")
  call mdl_add_service(pst%s%mdl,MDL_INPUT_DUST_GRAFIC,      pst,C_FUNLOC(r_input_dust_grafic),storage_size(pst%s%p%npart_tot)/32,0,"input_dust_grafic")
  call mdl_add_service(pst%s%mdl,MDL_INPUT_PART_ZOOM,        pst,C_FUNLOC(r_input_part_zoom),1,3,"input_part_zoom")
  call mdl_add_service(pst%s%mdl,MDL_INPUT_PART_ASCII,       pst,C_FUNLOC(r_input_part_ascii),storage_size(pst%s%p%npart_tot)/32,0,"input_part_ascii")
  call mdl_add_service(pst%s%mdl,MDL_INPUT_STAR_ASCII,       pst,C_FUNLOC(r_input_star_ascii),storage_size(pst%s%p%npart_tot)/32,0,"input_star_ascii")
  call mdl_add_service(pst%s%mdl,MDL_INPUT_SINK_ASCII,       pst,C_FUNLOC(r_input_sink_ascii),storage_size(pst%s%p%npart_tot)/32,0,"input_sink_ascii")
  call mdl_add_service(pst%s%mdl,MDL_INPUT_CORE_ASCII,       pst,C_FUNLOC(r_input_core_ascii),storage_size(pst%s%p%npart_tot)/32,0,"input_core_ascii")
  call mdl_add_service(pst%s%mdl,MDL_INPUT_PART_RESTART,     pst,C_FUNLOC(r_input_part_restart),(MDL_MAX_CPU+1),0,"input_part_restart")
  call mdl_add_service(pst%s%mdl,MDL_INPUT_PART_RAMSES,      pst,C_FUNLOC(r_input_part_ramses),(MDL_MAX_CPU+1),0,"input_part_ramses")
  call mdl_add_service(pst%s%mdl,MDL_INPUT_PART_GADGET,      pst,C_FUNLOC(r_input_part_gadget),MDL_MAX_CPU,6,"input_part_gadget")
  call mdl_add_service(pst%s%mdl,MDL_DEALLOCATE_GAS,         pst,C_FUNLOC(r_deallocate_gas),MDL_MAX_CPU,0,"input_part_gadget")
  call mdl_add_service(pst%s%mdl,MDL_NPART_MAX,              pst,C_FUNLOC(r_npart_max),0,1,"npart_max")
  call mdl_add_service(pst%s%mdl,MDL_INIT_FLAG,              pst,C_FUNLOC(r_init_flag),1,1,"init_flag")
  call mdl_add_service(pst%s%mdl,MDL_USER_FLAG,              pst,C_FUNLOC(r_user_flag),1,1,"user_flag")
  call mdl_add_service(pst%s%mdl,MDL_ENSURE_SUBGRID,         pst,C_FUNLOC(r_ensure_subgrid),1,0,"ensure_subgrid")
  call mdl_add_service(pst%s%mdl,MDL_ENSURE_REF_RULES,       pst,C_FUNLOC(r_ensure_ref_rules),1,0,"ensure_ref_rules")
  call mdl_add_service(pst%s%mdl,MDL_COLLECT_NOCT,           pst,C_FUNLOC(r_collect_noct),1,ncpu*storage_size(dummy)/32,"collect_noct")
  call mdl_add_service(pst%s%mdl,MDL_NOCT_TOT,               pst,C_FUNLOC(r_noct_tot),1,2,"noct_tot")
  call mdl_add_service(pst%s%mdl,MDL_NOCT_MIN,               pst,C_FUNLOC(r_noct_min),1,1,"noct_min")
  call mdl_add_service(pst%s%mdl,MDL_NOCT_MAX,               pst,C_FUNLOC(r_noct_max),1,1,"noct_max")
  call mdl_add_service(pst%s%mdl,MDL_NOCT_USED_MAX,          pst,C_FUNLOC(r_noct_used_max),1,1,"noct_used_max")
  call mdl_add_service(pst%s%mdl,MDL_INIT_REFINE_BASEGRID,   pst,C_FUNLOC(r_init_refine_basegrid),1,0,"init_refine_basegrid")
  call mdl_add_service(pst%s%mdl,MDL_INIT_REFINE_RESTART,    pst,C_FUNLOC(r_init_refine_restart),0,2*nhilbert*(pst%s%g%ncpu+1),"init_refine_restart")
  call mdl_add_service(pst%s%mdl,MDL_INIT_REFINE_RAMSES,     pst,C_FUNLOC(r_init_refine_ramses),0,2*nhilbert*(pst%s%g%ncpu+1),"init_refine_ramses")
  call mdl_add_service(pst%s%mdl,MDL_COLLECT_BOUND_KEY,      pst,C_FUNLOC(r_collect_bound_key),(MDL_MAX_CPU+1),nhilbert*(ncpu+1)*storage_size(dummy8)/32,"collect_bound_key")
  call mdl_add_service(pst%s%mdl,MDL_BROADCAST_BOUND_KEY,    pst,C_FUNLOC(r_broadcast_bound_key),nhilbert*(ncpu+1)*storage_size(dummy8)/32 + 1,0,"broadcast_bound_key")
  call mdl_add_service(pst%s%mdl,MDL_LOAD_BALANCE,           pst,C_FUNLOC(r_load_balance),1,0,"load_balance")
  call mdl_add_service(pst%s%mdl,MDL_BALANCE_PART,           pst,C_FUNLOC(r_balance_part),0,0,"balance_part")
  call mdl_add_service(pst%s%mdl,MDL_REFINE_FINE,            pst,C_FUNLOC(r_refine_fine),1,2,"refine_fine")
  call mdl_add_service(pst%s%mdl,MDL_CLEAN_DIRTY,            pst,C_FUNLOC(r_clean_dirty),1,0,"clean_dirty")
  call mdl_add_service(pst%s%mdl,MDL_SMOOTH_FINE,            pst,C_FUNLOC(r_smooth_fine),1,1,"smooth_fine")
  call mdl_add_service(pst%s%mdl,MDL_INPUT_HYDRO_CONDINIT,   pst,C_FUNLOC(r_input_hydro_condinit),1,0,"input_hydro_condinit")
  call mdl_add_service(pst%s%mdl,MDL_INPUT_HYDRO_GRAFIC,     pst,C_FUNLOC(r_input_hydro_grafic),0,0,"input_hydro_grafic")
  call mdl_add_service(pst%s%mdl,MDL_INPUT_REFMAP_GRAFIC,    pst,C_FUNLOC(r_input_refmap_grafic),0,0,"input_refmap_grafic")
  call mdl_add_service(pst%s%mdl,MDL_INPUT_HYDRO_GADGET,     pst,C_FUNLOC(r_input_hydro_gadget),0,0,"input_hydro_gadget")
  call mdl_add_service(pst%s%mdl,MDL_UPLOAD_FINE,            pst,C_FUNLOC(r_upload_fine),1,0,"upload_fine")
  call mdl_add_service(pst%s%mdl,MDL_KICK_DRIFT_PART,        pst,C_FUNLOC(r_kick_drift_part),0,0,"kick_drift_part")
  call mdl_add_service(pst%s%mdl,MDL_MASS_MIN_PART,          pst,C_FUNLOC(r_mass_min_part),0,2,"mass_min_part")
  call mdl_add_service(pst%s%mdl,MDL_BROADCAST_MP_MIN,       pst,C_FUNLOC(r_broadcast_mp_min),0,0,"broadcast_mp_min")
  call mdl_add_service(pst%s%mdl,MDL_OUTPUT_AMR,             pst,C_FUNLOC(r_output_amr),flen,0,"output_amr")
  call mdl_add_service(pst%s%mdl,MDL_OUTPUT_HYDRO,           pst,C_FUNLOC(r_output_hydro),flen,0,"output_hydro")
  call mdl_add_service(pst%s%mdl,MDL_OUTPUT_CLUMP,           pst,C_FUNLOC(r_output_clump),flen,0,"output_clump")
  call mdl_add_service(pst%s%mdl,MDL_OUTPUT_CORE_CLUMP, &
       pst,C_FUNLOC(r_output_core_clump),flen,0,"output_core_clump")
  call mdl_add_service(pst%s%mdl,MDL_OUTPUT_FRAME,           pst,C_FUNLOC(r_output_frame),2,2*pst%s%r%nw_frame*pst%s%r%nh_frame,"output_frame")
  call mdl_add_service(pst%s%mdl,MDL_OUTPUT_POISSON,         pst,C_FUNLOC(r_output_poisson),flen/4,0,"output_poisson")
  call mdl_add_service(pst%s%mdl,MDL_OUTPUT_PART,            pst,C_FUNLOC(r_output_part),flen/4,0,"output_part")
  call mdl_add_service(pst%s%mdl,MDL_BROADCAST_AEXP,         pst,C_FUNLOC(r_broadcast_aexp),4,0,"broadcast_aexp")
  call mdl_add_service(pst%s%mdl,MDL_COURANT_FINE,           pst,C_FUNLOC(r_courant_fine),1,8,"courant_fine")
  call mdl_add_service(pst%s%mdl,MDL_GODUNOV_FINE,           pst,C_FUNLOC(r_godunov_fine),1,0,"godunov_fine")
  call mdl_add_service(pst%s%mdl,MDL_SET_UNEW,               pst,C_FUNLOC(r_set_unew),1,0,"set_unew")
  call mdl_add_service(pst%s%mdl,MDL_SET_UOLD,               pst,C_FUNLOC(r_set_uold),1,0,"set_uold")
  call mdl_add_service(pst%s%mdl,MDL_NIMHD_FINE,             pst,C_FUNLOC(r_nimhd_fine),2,0,"nimhd_fine")
  call mdl_add_service(pst%s%mdl,MDL_COOLING_FINE,           pst,C_FUNLOC(r_cooling_fine),1,0,"cooling_fine")
  call mdl_add_service(pst%s%mdl,MDL_INIT_XION,              pst,C_FUNLOC(r_init_xion),1,0,"init_xion")
  call mdl_add_service(pst%s%mdl,MDL_STAR_FORMATION,         pst,C_FUNLOC(r_star_formation),1,2,"star_formation")
  call mdl_add_service(pst%s%mdl,MDL_CORE_TO_STAR_FORMATION, &
       pst,C_FUNLOC(r_core_to_star_formation),1,2,"core_to_star_formation")
  call mdl_add_service(pst%s%mdl,MDL_SINK_FORMATION,         pst,C_FUNLOC(r_sink_formation),1,2,"sink_formation")
  call mdl_add_service(pst%s%mdl,MDL_SINK_CLUMP,             pst,C_FUNLOC(r_sink_clump),1,2,"sink_clump")
  call mdl_add_service(pst%s%mdl,MDL_CORE_FORMATION, pst,C_FUNLOC(r_core_formation),1,2,"core_formation")
  call mdl_add_service(pst%s%mdl,MDL_CORE_CLUMP,     pst,C_FUNLOC(r_core_clump),1,2,"core_clump")
  call mdl_add_service(pst%s%mdl,MDL_TREE_FORMATION,         pst,C_FUNLOC(r_tree_formation),1,2,"tree_formation")
  call mdl_add_service(pst%s%mdl,MDL_TREE_CLUMP,             pst,C_FUNLOC(r_tree_clump),1,2,"tree_clump")
  call mdl_add_service(pst%s%mdl,MDL_THERMAL_FEEDBACK,       pst,C_FUNLOC(r_thermal_feedback),4,2,"thermal_feedback")
  call mdl_add_service(pst%s%mdl,MDL_MECHANICAL_FEEDBACK,    pst,C_FUNLOC(r_mechanical_feedback),4,2,"mechanical_feedback")
  call mdl_add_service(pst%s%mdl,MDL_SINK_EVOLUTION,         pst,C_FUNLOC(r_sink_evolution),1,2,"sink_evolution")
  call mdl_add_service(pst%s%mdl,MDL_SINK_MECHANICAL,        pst,C_FUNLOC(r_sink_mechanical),6,2,"sink_mechanical")
  call mdl_add_service(pst%s%mdl,MDL_CORE_EVOLUTION, pst,C_FUNLOC(r_core_evolution), &
       storage_size(dummy_core_acc_in)/32,storage_size(dummy_core_acc_out)/32,"core_evolution")
  call mdl_add_service(pst%s%mdl,MDL_CORE_RETIREMENT, &
       pst,C_FUNLOC(r_core_retirement), &
       storage_size(dummy_core_retire_in)/32, &
       storage_size(dummy_core_retire_out)/32,"core_retirement")
  call mdl_add_service(pst%s%mdl,MDL_SINK_MERGER,            pst,C_FUNLOC(r_sink_merger),1,0,"sink_merger")
  call mdl_add_service(pst%s%mdl,MDL_CORE_MERGER,    pst,C_FUNLOC(r_core_merger),1,0,"core_merger")
  call mdl_add_service(pst%s%mdl,MDL_NEWDT_PART,             pst,C_FUNLOC(r_newdt_part),0,0,"newdt_part")
  call mdl_add_service(pst%s%mdl,MDL_MAX_B_AND_Q,            pst,C_FUNLOC(r_max_B_and_Q),1,4,"max_B_and_Q")
  call mdl_add_service(pst%s%mdl,MDL_BROADCAST_DT,           pst,C_FUNLOC(r_broadcast_dt),24,0,"broadcast_dt")
  call mdl_add_service(pst%s%mdl,MDL_SYNCHRO_HYDRO_FINE,     pst,C_FUNLOC(r_synchro_hydro_fine),3,0,"synchro_hydro_fine")
  call mdl_add_service(pst%s%mdl,MDL_GRAVITY_HYDRO_FINE,     pst,C_FUNLOC(r_gravity_hydro_fine),1,0,"gravity_hydro_fine")
  call mdl_add_service(pst%s%mdl,MDL_SOURCE_HYDRO_FINE,      pst,C_FUNLOC(r_source_hydro_fine),1,0,"source_hydro_fine")
  call mdl_add_service(pst%s%mdl,MDL_CLUMP_FINDER,           pst,C_FUNLOC(r_clump_finder),1,0,"clump_finder")
  call mdl_add_service(pst%s%mdl,MDL_OUTPUT_LIGHTCONE,       pst,C_FUNLOC(r_output_lightcone),1,0,"output_lightcone")
  call mdl_add_service(pst%s%mdl,MDL_CLUMP_DEALLOC,          pst,C_FUNLOC(r_deallocate_clump),1,0,"deallocate_clump")
  call mdl_add_service(pst%s%mdl,MDL_SPLIT_PART,             pst,C_FUNLOC(r_split_part),0,0,"split_part")
  call mdl_add_service(pst%s%mdl,MDL_SORT_PART,              pst,C_FUNLOC(r_sort_part),0,0,"sort_part")
#ifdef GRAV
  call mdl_add_service(pst%s%mdl,MDL_MULTIPOLE_LEAF_CELLS,   pst,C_FUNLOC(r_multipole_leaf_cells),1,0,"multipole_leaf_cells")
  call mdl_add_service(pst%s%mdl,MDL_MULTIPOLE_SPLIT_CELLS,  pst,C_FUNLOC(r_multipole_split_cells),1,0,"multipole_split_cells")
  call mdl_add_service(pst%s%mdl,MDL_RESET_RHO,              pst,C_FUNLOC(r_reset_rho),1,0,"reset_rho")
  call mdl_add_service(pst%s%mdl,MDL_CIC_MULTIPOLE,          pst,C_FUNLOC(r_cic_multipole),1,0,"cic_multipole")
  call mdl_add_service(pst%s%mdl,MDL_CIC_PART,               pst,C_FUNLOC(r_cic_part),0,0,"cic_part")
  call mdl_add_service(pst%s%mdl,MDL_COLLECT_MULTIPOLE,      pst,C_FUNLOC(r_collect_multipole),1,storage_size(pst%s%g%multipole)/32,"collect_multipole")
  call mdl_add_service(pst%s%mdl,MDL_COMPUTE_RHOMAX,         pst,C_FUNLOC(r_compute_rhomax),1,2,"compute_rhomax")
  call mdl_add_service(pst%s%mdl,MDL_BROADCAST_MULTIPOLE,    pst,C_FUNLOC(r_broadcast_multipole),storage_size(pst%s%g%multipole)/32,0,"broadcast_multipole")
  call mdl_add_service(pst%s%mdl,MDL_SAVE_PHI_OLD,           pst,C_FUNLOC(r_save_phi_old),1,0,"save_phi_old")
  call mdl_add_service(pst%s%mdl,MDL_FORCE_ANALYTIC,         pst,C_FUNLOC(r_force_analytic),1,0,"force_analytic")
  call mdl_add_service(pst%s%mdl,MDL_GRADIENT_PHI,           pst,C_FUNLOC(r_gradient_phi),2,0,"gradient_phi")
  call mdl_add_service(pst%s%mdl,MDL_COMPUTE_EPOT,           pst,C_FUNLOC(r_compute_epot),1,2,"compute_epot")
  call mdl_add_service(pst%s%mdl,MDL_MAKE_INITIAL_PHI,       pst,C_FUNLOC(r_make_initial_phi),2,0,"make_initial_phi")
  call mdl_add_service(pst%s%mdl,MDL_RECURRENCE_ON_P,        pst,C_FUNLOC(r_recurrence_on_p),0,0,"recurrence_on_p")
  call mdl_add_service(pst%s%mdl,MDL_RECURRENCE_X_AND_R,     pst,C_FUNLOC(r_recurrence_x_and_r),0,0,"recurrence_x_and_r")
  call mdl_add_service(pst%s%mdl,MDL_CMP_RESIDUAL_CG,        pst,C_FUNLOC(r_cmp_residual_cg),0,0,"cmp_residual_cg")
  call mdl_add_service(pst%s%mdl,MDL_CMP_R2_CG,              pst,C_FUNLOC(r_cmp_r2_cg),0,2,"cmp_r2_cg")
  call mdl_add_service(pst%s%mdl,MDL_CMP_PAP_CG,             pst,C_FUNLOC(r_cmp_pAp_cg),0,2,"cmp_pAp_cg")
  call mdl_add_service(pst%s%mdl,MDL_CMP_RHS_NORM,           pst,C_FUNLOC(r_cmp_rhs_norm),0,2,"cmp_rhs_norm")
  call mdl_add_service(pst%s%mdl,MDL_INIT_MG,                pst,C_FUNLOC(r_init_mg),1,0,"init_mg")
  call mdl_add_service(pst%s%mdl,MDL_BUILD_MG,               pst,C_FUNLOC(r_build_mg),1,0,"build_mg")
  call mdl_add_service(pst%s%mdl,MDL_CLEANUP_MG,             pst,C_FUNLOC(r_cleanup_mg),1,0,"cleanup_mg")
  call mdl_add_service(pst%s%mdl,MDL_MAKE_MASK,              pst,C_FUNLOC(r_make_mask),4,0,"make_mask")
  call mdl_add_service(pst%s%mdl,MDL_MAKE_BC_RHS,            pst,C_FUNLOC(r_make_bc_rhs),2,0,"make_bc_rhs")
  call mdl_add_service(pst%s%mdl,MDL_RESTRICT_MASK,          pst,C_FUNLOC(r_restrict_mask),1,1,"restrict_mask")
  call mdl_add_service(pst%s%mdl,MDL_CMP_RESIDUAL_MG,        pst,C_FUNLOC(r_cmp_residual_mg),2,0,"cmp_residual_mg")
  call mdl_add_service(pst%s%mdl,MDL_GAUSS_SEIDEL_MG,        pst,C_FUNLOC(r_gauss_seidel_mg),4,0,"gauss_seidel_mg")
  call mdl_add_service(pst%s%mdl,MDL_RESET_CORRECTION,       pst,C_FUNLOC(r_reset_correction),1,0,"reset_correction")
  call mdl_add_service(pst%s%mdl,MDL_RESTRICT_RESIDUAL,      pst,C_FUNLOC(r_restrict_residual),1,0,"restrict_residual")
  call mdl_add_service(pst%s%mdl,MDL_INTERPOLATE_AND_CORRECT,pst,C_FUNLOC(r_interpolate_and_correct),1,0,"interpolate_and_correct")
  call mdl_add_service(pst%s%mdl,MDL_SET_SCAN_FLAG,          pst,C_FUNLOC(r_set_scan_flag),2,0,"set_scan_flag")
  call mdl_add_service(pst%s%mdl,MDL_CMP_RESIDUAL_NORM2,     pst,C_FUNLOC(r_cmp_residual_norm2),1,2,"cmp_residual_norm2")
  call mdl_add_service(pst%s%mdl,MDL_CMP_RHS_NORM2,          pst,C_FUNLOC(r_cmp_rhs_norm2),1,2,"cmp_rhs_norm2")
#endif
  call mdl_add_service(pst%s%mdl,MDL_INIT_RT,                pst,C_FUNLOC(r_init_rt),0,0,"init_rt")
  call mdl_add_service(pst%s%mdl,MDL_RT_UPLOAD_FINE,         pst,C_FUNLOC(r_rt_upload_fine),1,0,"rt_upload_fine")
  call mdl_add_service(pst%s%mdl,MDL_RT_INPUT_CONDINIT,      pst,C_FUNLOC(r_rt_input_condinit),1,0,"rt_input_condinit")
  call mdl_add_service(pst%s%mdl,MDL_RT_INPUT_SOURCE_REGIONS,pst,C_FUNLOC(r_rt_input_source_regions),1,0,"rt_input_source_regions")
  call mdl_add_service(pst%s%mdl,MDL_OUTPUT_RT,              pst,C_FUNLOC(r_output_rt),flen,0,"output_rt")
  call mdl_add_service(pst%s%mdl,MDL_RT_GODUNOV_FINE,        pst,C_FUNLOC(r_rt_godunov_fine),1,0,"rt_godunov_fine")
  call mdl_add_service(pst%s%mdl,MDL_SET_EMISSIVITY,         pst,C_FUNLOC(r_set_emissivity),1,0,"set_emissivity")
  call mdl_add_service(pst%s%mdl,MDL_STAR_RT_FEEDBACK,       pst,C_FUNLOC(r_star_rt_feedback),1,0,"star_rt_feedback")
  call mdl_add_service(pst%s%mdl,MDL_SET_RTUNEW,             pst,C_FUNLOC(r_set_rtunew),1,0,"set_rtunew")
  call mdl_add_service(pst%s%mdl,MDL_SET_RTUOLD,             pst,C_FUNLOC(r_set_rtuold),1,0,"set_rtuold")
  call mdl_add_service(pst%s%mdl,MDL_RT_NEQ_UPDATES,         pst,C_FUNLOC(r_rt_neq_updates),1,0,"rt_neq_updates")
  call mdl_add_service(pst%s%mdl,MDL_CHECK_PART_EMISSION,    pst,C_FUNLOC(r_check_part_emission),0,0,"check_part_emission")
  call mdl_add_service(pst%s%mdl,MDL_INIT_CR,                pst,C_FUNLOC(r_init_cr),0,0,"init_cr")
  call mdl_add_service(pst%s%mdl,MDL_CR_UPLOAD_FINE,         pst,C_FUNLOC(r_cr_upload_fine),1,0,"cr_upload_fine")
  call mdl_add_service(pst%s%mdl,MDL_CR_INPUT_CONDINIT,      pst,C_FUNLOC(r_cr_input_condinit),1,0,"cr_input_condinit")
  call mdl_add_service(pst%s%mdl,MDL_OUTPUT_CR,              pst,C_FUNLOC(r_output_cr),flen,0,"output_cr")
  call mdl_add_service(pst%s%mdl,MDL_CR_GODUNOV_FINE,        pst,C_FUNLOC(r_cr_godunov_fine),1,0,"cr_godunov_fine")
  call mdl_add_service(pst%s%mdl,MDL_CR_SOURCE_TERMS,        pst,C_FUNLOC(r_cr_source_terms),1,0,"cr_source_terms")
  call mdl_add_service(pst%s%mdl,MDL_SET_CRUNEW,             pst,C_FUNLOC(r_set_crunew),1,0,"set_crunew")
  call mdl_add_service(pst%s%mdl,MDL_SET_CRUOLD,             pst,C_FUNLOC(r_set_cruold),1,0,"set_cruold")
  call mdl_add_service(pst%s%mdl,MDL_CONSERVE_CR_FLUX,       pst,C_FUNLOC(r_conserve_cr_flux),1,0,"conserve_cr_flux")
  call mdl_add_service(pst%s%mdl,MDL_CR_UPDATES,             pst,C_FUNLOC(r_cr_updates),1,0,"cr_updates")
  call mdl_add_service(pst%s%mdl,MDL_BROADCAST_CR_C,         pst,C_FUNLOC(r_broadcast_cr_c),24,0,"broadcast_cr_c")
#if defined(_CUDA) || defined(_METAL)
  call mdl_add_service(pst%s%mdl,MDL_SET_GRID_DEVICE,        pst,C_FUNLOC(r_set_grid_device),0,0,"set_grid_device")
  call mdl_add_service(pst%s%mdl,MDL_TRANSFER_GRID_HOST,     pst,C_FUNLOC(r_transfer_grid_host),0,0,"transfer_grid_host")
#endif
  call mdl_add_service(pst%s%mdl,MDL_INIT_TURB,              pst,C_FUNLOC(r_init_turb),1,0,"init_turb")
  call mdl_add_service(pst%s%mdl,MDL_UPDATE_TURB,            pst,C_FUNLOC(r_update_turb),1,0,"update_turb")
  call mdl_add_service(pst%s%mdl,MDL_DRIVE_TURB,             pst,C_FUNLOC(r_drive_turb),1,0,"drive_turb")
  call mdl_add_service(pst%s%mdl,MDL_TURB_HYDRO,             pst,C_FUNLOC(r_turb_hydro),1,0,"turb_hydro")
end function worker_init
!##############################################################
!##############################################################
!##############################################################
!##############################################################
end subroutine mdl_init
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine mdl_wait(pst)
  use ramses_commons, only: pst_t
  use call_back, only: call_back_f, ramses_function
  use mdl_module
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_F_PROCPOINTER, C_PTR, C_LOC, C_F_POINTER
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  type(pst_t)::pst

  type generic
  end type generic

  procedure(ramses_function),pointer::mdl_function
  type(generic),pointer::ipar,opar

#ifndef WITHOUTMPI

  logical::stop_order_received,order_received
  integer::order_id,order_tag=101,output_tag=203,output_id
  integer,dimension(MPI_STATUS_SIZE)::order_status,output_status
  integer::input_size
  integer::output_size,source_cpu,function_id
  integer,dimension(:),allocatable,target::input_array,output_array
  integer,dimension(1)::dummy
  integer,dimension(1:32)::header

  associate(s=>pst%s,mdl=>pst%s%mdl)

  ! Post the first RECV for a launch order
  call MPI_IRECV(mdl%mpi_input_buffer,mdl%MDL_INPUT_MAXSIZE+32,MPI_INTEGER,MPI_ANY_SOURCE,order_tag,MPI_COMM_WORLD,order_id,info)

  stop_order_received=.false.
  do while(.NOT. stop_order_received)

     call MPI_Test(order_id,order_received,order_status,info)

     if(order_received)then

        ! Execute call-back functions
        header=mdl%mpi_input_buffer(1:32)
        function_id=header(1)
        if(function_id==0)stop_order_received=.true.
        ! Get source cpu
        source_cpu=order_status(MPI_SOURCE)

        ! Allocate input and output arrays
        input_size=header(2)
        output_size=header(3)

        allocate(input_array(1:input_size))
        input_array(1:input_size)=mdl%mpi_input_buffer(33:32+input_size)

        allocate(output_array(1:output_size))
        output_array=0

        ! Launch the corresponding call-back function
        call c_f_pointer(c_loc(input_array),ipar)
        call c_f_pointer(c_loc(output_array),opar)
        CALL C_F_PROCPOINTER (mdl%callback(function_id), mdl_function)
        call mdl_function(pst,ipar,input_size,opar,output_size)

        ! Deallocate input array
        deallocate(input_array)

        ! Always send the output back to the source cpu (even if not allocated = handshake)
        if(output_size>0)then
           call MPI_ISEND(output_array,output_size,MPI_INTEGER,source_cpu,output_tag,MPI_COMM_WORLD,output_id,info)
        else
           dummy(1)=1
           call MPI_ISEND(dummy,1,MPI_INTEGER,source_cpu,output_tag,MPI_COMM_WORLD,output_id,info)
        endif
        call MPI_WAIT(output_id,output_status,info)

        ! Deallocate output array
        deallocate(output_array)

        ! Post a new RECV for the next launch order
        if(.NOT. stop_order_received)then
           call MPI_IRECV(mdl%mpi_input_buffer,mdl%MDL_INPUT_MAXSIZE+32,MPI_INTEGER,MPI_ANY_SOURCE,order_tag,MPI_COMM_WORLD,order_id,info)
        endif

     endif
  end do

  end associate

#endif

end subroutine mdl_wait
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine init_cache(mdl)
  use mdl_module
  use cache_commons
  implicit none
  type(mdl_t)::mdl

  integer::ncpu,ibuf
  type(msg_large_realdp)::dummy_large_realdp

  ncpu=mdl_threads(mdl)

#ifndef WITHOUTMPI

  ! Allocate all communication and cache-related variables
  allocate(mdl%reply_id(1:ncpu))

  ! Compute largest possible message size
  mdl%size_request_array=1+ndim
  mdl%size_msg_array=storage_size(dummy_large_realdp)/32 ! Largest message size
  mdl%size_flush_array=1+(1+ndim+mdl%size_msg_array)*nflushmax
  mdl%size_fetch_array=2+(1+ndim+mdl%size_msg_array)*ntilemax

  ! Allocate large enough message communication buffers
  allocate(mdl%recv_request_array(1:mdl%size_request_array))
  allocate(mdl%send_request_array(1:mdl%size_request_array))
  allocate(mdl%recv_fetch_array(1:mdl%size_fetch_array))
  allocate(mdl%recv_flush_array(1:mdl%size_flush_array))

  allocate(mdl%cpu2buf_fetch(1:ncpu))
  allocate(mdl%cpu2buf_flush(1:ncpu))
  allocate(mdl%send_fetch(1:ncpu))
  allocate(mdl%send_flush(1:ncpu))

  mdl%nbuffer_fetch=MIN(ncpu,nbuffermax)
  mdl%nbuffer_flush=MIN(ncpu,nbuffermax)
  do ibuf=1,mdl%nbuffer_fetch
     allocate(mdl%send_fetch(ibuf)%array(mdl%size_fetch_array))
  end do
  do ibuf=1,mdl%nbuffer_flush
     allocate(mdl%send_flush(ibuf)%array(mdl%size_flush_array))
  end do

  mdl%ibuffer_fetch=0
  mdl%ibuffer_flush=0
  mdl%cpu2buf_fetch=0
  mdl%cpu2buf_flush=0

#endif

end subroutine init_cache
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine kill_cache(mdl)
  use mdl_module
  use cache_commons
  implicit none
  type(mdl_t)::mdl

  integer::ibuf

#ifndef WITHOUTMPI

  do ibuf=1,mdl%nbuffer_fetch
     deallocate(mdl%send_fetch(ibuf)%array)
  end do
  do ibuf=1,mdl%nbuffer_flush
     deallocate(mdl%send_flush(ibuf)%array)
  end do

  deallocate(mdl%cpu2buf_fetch)
  deallocate(mdl%cpu2buf_flush)
  deallocate(mdl%send_fetch)
  deallocate(mdl%send_flush)

  deallocate(mdl%recv_request_array)
  deallocate(mdl%send_request_array)
  deallocate(mdl%recv_fetch_array)
  deallocate(mdl%recv_flush_array)

  deallocate(mdl%reply_id)

#endif

end subroutine kill_cache
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine init_cache_clump(mdl)
  use mdl_module
  use cache_commons
  implicit none
  type(mdl_t)::mdl

  integer::ncpu,ibuf
  type(msg_unbind_clump)::dummy_large_clump

  ncpu=mdl_threads(mdl)

#ifndef WITHOUTMPI

  ! Allocate all communication and cache-related variables
  allocate(mdl%reply_id_clump(1:ncpu))

  ! Compute largest possible message size
  mdl%size_request_array_clump=2
  mdl%size_msg_array_clump=storage_size(dummy_large_clump)/32 ! Largest message size
  mdl%size_flush_array_clump=1+(2+mdl%size_msg_array_clump)*nflushmax
  mdl%size_fetch_array_clump=2+(2+mdl%size_msg_array_clump)*ntilemax

  ! Allocate large enough message communication buffers
  allocate(mdl%recv_request_array_clump(1:mdl%size_request_array_clump))
  allocate(mdl%send_request_array_clump(1:mdl%size_request_array_clump))
  allocate(mdl%recv_fetch_array_clump(1:mdl%size_fetch_array_clump))
  allocate(mdl%recv_flush_array_clump(1:mdl%size_flush_array_clump))

  allocate(mdl%cpu2buf_fetch_clump(1:ncpu))
  allocate(mdl%cpu2buf_flush_clump(1:ncpu))
  allocate(mdl%send_fetch_clump(1:ncpu))
  allocate(mdl%send_flush_clump(1:ncpu))

  mdl%nbuffer_fetch_clump=MIN(ncpu,nbuffermax)
  mdl%nbuffer_flush_clump=MIN(ncpu,nbuffermax)
  do ibuf=1,mdl%nbuffer_fetch_clump
     allocate(mdl%send_fetch_clump(ibuf)%array(mdl%size_fetch_array_clump))
  end do
  do ibuf=1,mdl%nbuffer_flush_clump
     allocate(mdl%send_flush_clump(ibuf)%array(mdl%size_flush_array_clump))
  end do

  mdl%ibuffer_fetch_clump=0
  mdl%ibuffer_flush_clump=0
  mdl%cpu2buf_fetch_clump=0
  mdl%cpu2buf_flush_clump=0

#endif

end subroutine init_cache_clump
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine kill_cache_clump(mdl)
  use mdl_module
  use cache_commons
  implicit none
  type(mdl_t)::mdl

  integer::ibuf

#ifndef WITHOUTMPI

  do ibuf=1,mdl%nbuffer_fetch_clump
     deallocate(mdl%send_fetch_clump(ibuf)%array)
  end do
  do ibuf=1,mdl%nbuffer_flush_clump
     deallocate(mdl%send_flush_clump(ibuf)%array)
  end do

  deallocate(mdl%cpu2buf_fetch_clump)
  deallocate(mdl%cpu2buf_flush_clump)
  deallocate(mdl%send_fetch_clump)
  deallocate(mdl%send_flush_clump)

  deallocate(mdl%recv_request_array_clump)
  deallocate(mdl%send_request_array_clump)
  deallocate(mdl%recv_fetch_array_clump)
  deallocate(mdl%recv_flush_array_clump)

  deallocate(mdl%reply_id_clump)

#endif

end subroutine kill_cache_clump
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine check_mail(mdl,comm_id)
  use mdl_module
  use amr_parameters, only: ndim, nhilbert, twotondim
  use amr_commons, only: mesh_t
  use clfind_commons, only: clump_t
  use cache, only: cache_oct, get_tile
  use cache_commons
  use hilbert
  use hash
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(mdl_t)::mdl
  integer::comm_id
  !
  ! This routine checks for incoming messages.
  ! It can be a cache request, and the routine
  ! assembles a fetch message and sends it back.
  ! It can be a flush message, and the routine
  ! unpacks it and combine it in the local memory.
  !
  integer::i,info,iskip
  integer::grid_cpu,peak_cpu,ilevel,itile,ntile_reply,nflush,ibuf
  integer::buffer_size_fetch_array,buffer_size_msg_array
  integer::buffer_size_fetch_array_clump,buffer_size_msg_array_clump
  logical::comm_completed,request_received,flush_received=.false.
#ifndef WITHOUTMPI
  integer,dimension(MPI_STATUS_SIZE)::reply_status,request_status,flush_status,comm_status
  integer,dimension(MPI_STATUS_SIZE)::reply_status_clump,request_status_clump,flush_status_clump
#endif
  integer(kind=8), dimension(1:nhilbert)::hk
  integer(kind=8), dimension(1:ndim)::ix
  integer(kind=8), dimension(0:ndim)::hash_key,hash_child
  integer(kind=8)::global_peak_id
  integer::local_peak_id
  integer::child
  type(cache_oct),dimension(1:ntilemax)::tile
  type(msg_large_realdp)::dummy_large_realdp
  type(msg_unbind_clump)::dummy_large_clump
  logical::child_exist

#ifndef WITHOUTMPI

  associate(m=>mdl%m,c=>mdl%c)

  comm_completed=.false.
  do while (.not. comm_completed)

     ! USE A NESTED DO WHILE LOOP FOR THE FOLLOWING TESTS

     if(mdl%cache_opened)then
        !===========================
        ! Check for incoming request
        !===========================
        request_received=.true.
        do while(request_received)
           call MPI_TEST(mdl%request_id,request_received,request_status,info)
           if(request_received)then

              ! Assemble a reply and send it back
              ilevel=mdl%recv_request_array(1)
              hash_key(0)=ilevel
              hash_key(1:ndim)=mdl%recv_request_array(2:ndim+1)
              child=hash_getp(m%grid_dict,hash_key)
              grid_cpu=request_status(MPI_SOURCE)+1

              ! Identify the corresponding communication buffer
              ibuf=mdl%cpu2buf_fetch(grid_cpu)

              ! If none is available, get the next free one
              if(ibuf==0)then
                 mdl%ibuffer_fetch=mdl%ibuffer_fetch+1
                 ! If there is no communication buffer left, allocate a new one
                 if(mdl%ibuffer_fetch>mdl%nbuffer_fetch)then
                    mdl%nbuffer_fetch=mdl%nbuffer_fetch+1
                    buffer_size_msg_array=storage_size(dummy_large_realdp)/32 ! Largest message size
                    buffer_size_fetch_array=2+(1+ndim+buffer_size_msg_array)*ntilemax
                    allocate(mdl%send_fetch(mdl%nbuffer_fetch)%array(1:buffer_size_fetch_array))
                 endif
                 mdl%cpu2buf_fetch(grid_cpu)=mdl%ibuffer_fetch
                 ibuf=mdl%ibuffer_fetch
              endif

              child_exist=.true.
              child_exist=child_exist.and.(child>0)
              child_exist=child_exist.and.(child<=m%ngridmax)

              ! If grid does not exist, send a null reply
              if(.not.child_exist)then

                 ! Store type corresponding to a null reply
                 iskip=1
                 mdl%send_fetch(ibuf)%array(iskip)=-1

                 ! Otherwise, assemble a proper reply with a complete tile
              else
                 iskip=1
                 ! Record the location of where we want the keys to be stored, then ask for a tile
                 call get_tile(m,child,tile,ntile_reply)

                 ! Store type of reply and number of entries
                 mdl%send_fetch(ibuf)%array(iskip)=1
                 iskip=iskip+1
                 mdl%send_fetch(ibuf)%array(iskip)=ntile_reply
                 iskip=iskip+1

                 ! Store data, depending on reply type
                 do i=1,ntile_reply
                    mdl%send_fetch(ibuf)%array(iskip:iskip+ndim)=int(tile(i)%key(0:ndim),kind=4)
                    iskip=iskip+1+ndim ! Skip over the key already present
                    call pack_fetch%proc(m,tile(i)%igrid,mdl%size_msg_array,mdl%send_fetch(ibuf)%array(iskip:iskip+mdl%size_msg_array-1))
                    iskip=iskip+mdl%size_msg_array
                 end do
              endif

              ! Wait for the old SEND to free memory in corresponding MPI buffer
              call MPI_WAIT(mdl%reply_id(grid_cpu),reply_status,info)

              ! Send back the reply
              iskip=1
              call MPI_ISEND(mdl%send_fetch(ibuf)%array(iskip),mdl%size_fetch_array,MPI_INTEGER,grid_cpu-1,msg_tag,MPI_COMM_WORLD,mdl%reply_id(grid_cpu),info)

              !=================================
              ! Post a new RECV for request
              !=================================
              call MPI_IRECV(mdl%recv_request_array,mdl%size_request_array,MPI_INTEGER,MPI_ANY_SOURCE,request_tag,MPI_COMM_WORLD,mdl%request_id,info)

           endif
        end do

        !=========================
        ! Check for incoming flush
        !=========================
        flush_received=.true.
        do while(flush_received)
           call MPI_TEST(mdl%flush_id,flush_received,flush_status,info)
           if(flush_received)then

              ! Combine received data to local memory only if grid exists
              if(mdl%combiner_rule.eq.COMBINER_EXIST)then
                 iskip=1
                 nflush=mdl%recv_flush_array(iskip)
                 iskip=iskip+1

                 do i=1,nflush
                    ilevel=mdl%recv_flush_array(iskip)
                    hash_child(0)=ilevel
                    hash_child(1:ndim)=mdl%recv_flush_array(iskip+1:iskip+ndim)
                    iskip=iskip+ndim+1
                    ! Get grid from hash table
                    child=hash_getp(m%grid_dict,hash_child)
                    if(child>0)then
                       call unpack_flush%proc(m,child,mdl%size_msg_array,mdl%recv_flush_array(iskip:iskip+mdl%size_msg_array-1),hash_child)
                    endif

                    iskip=iskip+mdl%size_msg_array

                 end do

              endif

              ! Combine received data to local memory only if grid does not exist
              if(mdl%combiner_rule.eq.COMBINER_CREATE)then
                 iskip=1
                 nflush=mdl%recv_flush_array(iskip)
                 iskip=iskip+1

                 do i=1,nflush
                    ilevel=mdl%recv_flush_array(iskip)
                    hash_child(0)=ilevel
                    hash_child(1:ndim)=mdl%recv_flush_array(iskip+1:iskip+ndim)
                    iskip=iskip+ndim+1

                    ! Create new grid if grid does not exist
                    if(hash_getp(m%grid_dict,hash_child)<=0)then

                       ! Compute Hilbert keys of new octs
                       ix(1:ndim)=hash_child(1:ndim)
                       hk(1:nhilbert)=hilbert_key(ix,ilevel-1)

                       ! Set grid index to a virtual grid in local main memory
                       child=m%ifree

                       ! Go to next main memory free line
                       m%ifree=m%ifree+1
                       if(m%ifree.GT.m%ngridmax)then
                          write(*,*)'No more free memory'
                          write(*,*)'Increase ngridmax'
                          call mdl_abort(mdl)
                       endif
                       m%grid(child)%hkey(1:nhilbert)=hk(1:nhilbert)
                       m%flag1(1:twotondim,child)=0
                       m%flag2(1:twotondim,child)=0

                       ! Insert new grid in hash table
                       call hash_setp(m%grid_dict,hash_child,child)

                       ! Unpack message content
                       call unpack_flush%proc(m,child,mdl%size_msg_array,mdl%recv_flush_array(iskip:iskip+mdl%size_msg_array-1),hash_child)

                    endif

                    iskip=iskip+mdl%size_msg_array

                 end do

              endif

              !=================================
              ! Post a new RECV for flush
              !=================================
              call MPI_IRECV(mdl%recv_flush_array,mdl%size_flush_array,MPI_INTEGER,MPI_ANY_SOURCE,flush_tag,MPI_COMM_WORLD,mdl%flush_id,info)

           endif
        end do
     endif

     if(mdl%cache_opened_clump)then
        !=================================
        ! Check for incoming clump request
        !=================================
        request_received=.true.
        do while(request_received)
           call MPI_TEST(mdl%request_id_clump,request_received,request_status_clump,info)
           if(request_received)then

              ! Assemble a reply and send it back
              global_peak_id=transfer(mdl%recv_request_array_clump(1:2),global_peak_id)
              local_peak_id=int(global_peak_id-c%npeak_cum(mdl_self(mdl)-1),kind=4)
              peak_cpu=request_status_clump(MPI_SOURCE)+1

              ! Identify the corresponding communication buffer
              ibuf=mdl%cpu2buf_fetch_clump(peak_cpu)

              ! If none is available, get the next free one or allocate a new one
              if(ibuf==0)then
                 mdl%ibuffer_fetch_clump=mdl%ibuffer_fetch_clump+1
                 if(mdl%ibuffer_fetch_clump>mdl%nbuffer_fetch_clump)then
                    mdl%nbuffer_fetch_clump=mdl%nbuffer_fetch_clump+1
                    buffer_size_msg_array_clump=storage_size(dummy_large_clump)/32 ! Largest message size
                    buffer_size_fetch_array_clump=2+(2+buffer_size_msg_array_clump)*ntilemax
                    allocate(mdl%send_fetch_clump(mdl%nbuffer_fetch_clump)%array(1:buffer_size_fetch_array_clump))
                 endif
                 mdl%cpu2buf_fetch_clump(peak_cpu)=mdl%ibuffer_fetch_clump
                 ibuf=mdl%ibuffer_fetch_clump
              endif

              ! Assemble a proper reply with a complete tile
              iskip=1
              itile=(local_peak_id-1)/ntilemax
              ntile_reply=MIN(c%npeak-itile*ntilemax,ntilemax)

              ! Store type of reply and number of entries
              mdl%send_fetch_clump(ibuf)%array(iskip)=1
              iskip=iskip+1
              mdl%send_fetch_clump(ibuf)%array(iskip)=ntile_reply
              iskip=iskip+1

              ! Store data, depending on reply type
              do i=1,ntile_reply
                 local_peak_id=itile*ntilemax+i
                 global_peak_id=c%npeak_cum(mdl_self(mdl)-1)+local_peak_id
                 mdl%send_fetch_clump(ibuf)%array(iskip:iskip+1)=transfer(global_peak_id,local_peak_id,2)
                 iskip=iskip+2
                 call pack_fetch_clump%proc(c,local_peak_id,mdl%size_msg_array_clump,mdl%send_fetch_clump(ibuf)%array(iskip:iskip+mdl%size_msg_array_clump-1))
                 iskip=iskip+mdl%size_msg_array_clump
              end do

              ! Wait for the old SEND to free memory in corresponding MPI buffer
              call MPI_WAIT(mdl%reply_id_clump(peak_cpu),reply_status_clump,info)

              ! Send back the reply
              iskip=1
              call MPI_ISEND(mdl%send_fetch_clump(ibuf)%array(iskip),mdl%size_fetch_array_clump,MPI_INTEGER,peak_cpu-1,msg_tag_clump,MPI_COMM_WORLD,mdl%reply_id_clump(peak_cpu),info)

              !=================================
              ! Post a new RECV for request
              !=================================
              call MPI_IRECV(mdl%recv_request_array_clump,mdl%size_request_array_clump,MPI_INTEGER,MPI_ANY_SOURCE,request_tag_clump,MPI_COMM_WORLD,mdl%request_id_clump,info)

           endif
        end do

        !===============================
        ! Check for incoming clump flush
        !===============================
        flush_received=.true.
        do while(flush_received)
           call MPI_TEST(mdl%flush_id_clump,flush_received,flush_status_clump,info)
           if(flush_received)then

              ! Combine received data to local memory
              iskip=1
              nflush=mdl%recv_flush_array_clump(iskip)
              iskip=iskip+1
              do i=1,nflush
                 global_peak_id=transfer(mdl%recv_flush_array_clump(iskip:iskip+1),global_peak_id)
                 iskip=iskip+2
                 local_peak_id=int(global_peak_id-c%npeak_cum(mdl_self(mdl)-1),kind=4)
                 call unpack_flush_clump%proc(c,local_peak_id,mdl%size_msg_array_clump,mdl%recv_flush_array_clump(iskip:iskip+mdl%size_msg_array_clump-1))
                 iskip=iskip+mdl%size_msg_array_clump
              end do

              !=================================
              ! Post a new RECV for flush
              !=================================
              call MPI_IRECV(mdl%recv_flush_array_clump,mdl%size_flush_array_clump,MPI_INTEGER,MPI_ANY_SOURCE,flush_tag_clump,MPI_COMM_WORLD,mdl%flush_id_clump,info)

           endif
        end do
     endif

     !=================================
     ! Check for input comm. completion
     !=================================
     if(comm_id==MPI_REQUEST_NULL)then
        comm_completed=.true.
     else
        call MPI_TEST(comm_id,comm_completed,comm_status,info)
     endif
  end do

  end associate

#endif
end subroutine check_mail
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine destage(mdl,igrid)
  use amr_parameters, only: ndim
  use mdl_module
  use cache_commons
  use hash
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(mdl_t)::mdl
  integer::igrid
  !
  ! This routine frees the cache memory.
  ! It assembles flush messages, and when the message
  ! buffer is full, it sends it to the target CPU.
  !
  integer::info,icache,iflush,grid_cpu,ibuf
  integer::send_flush_id,iskip,nflush
  integer::buffer_size_flush_array,buffer_size_msg_array
  integer(kind=8),dimension(0:ndim)::hash_key
  type(msg_large_realdp)::dummy_large_realdp

  associate(m=>mdl%m)

  hash_key(0)=m%grid(igrid)%lev
  hash_key(1:ndim)=m%grid(igrid)%ckey(1:ndim)

  if(hash_getp(m%grid_dict,hash_key)<=0)then
     write(*,*)'PE ',mdl_self(mdl),' trying to free non existing grid'
     stop
  endif

  call hash_free(m%grid_dict,hash_key)

#ifndef WITHOUTMPI

  ! Check if the destage requires a flush
  icache=igrid-m%ngridmax

  if(m%dirty(icache))then
     grid_cpu=m%parent_cpu(icache)
     m%dirty(icache)=.false.

     ! Identify the corresponding communication buffer
     ibuf=mdl%cpu2buf_flush(grid_cpu)

     ! If none is available, get the next free one
     if(ibuf==0)then
        mdl%ibuffer_flush=mdl%ibuffer_flush+1
        ! If there is no communication buffer left, allocate a new one
        if(mdl%ibuffer_flush>mdl%nbuffer_flush)then
           mdl%nbuffer_flush=mdl%nbuffer_flush+1
           buffer_size_msg_array=storage_size(dummy_large_realdp)/32 ! Largest message size
           buffer_size_flush_array=1+(1+ndim+buffer_size_msg_array)*nflushmax
           allocate(mdl%send_flush(mdl%nbuffer_flush)%array(1:buffer_size_flush_array))
           mdl%send_flush(mdl%nbuffer_flush)%array(1)=0
        endif
        mdl%cpu2buf_flush(grid_cpu)=mdl%ibuffer_flush
        ibuf=mdl%ibuffer_flush
     endif

     ! Filling the flush buffer
     iskip=1
     nflush=mdl%send_flush(ibuf)%array(iskip)

     ! If buffer full, send it to remote CPU.
     if(nflush==nflushmax)then
        ! Post send
        call MPI_ISSEND(mdl%send_flush(ibuf)%array(iskip),mdl%size_flush_array,MPI_INTEGER,grid_cpu-1,flush_tag,MPI_COMM_WORLD,send_flush_id,info)
        ! While waiting for completion, check on incoming messages and perform actions
        call check_mail(mdl,send_flush_id)
        ! Reset counter
        mdl%send_flush(ibuf)%array(iskip)=0
     endif

     ! Increment counter
     mdl%send_flush(ibuf)%array(iskip)=mdl%send_flush(ibuf)%array(iskip)+1

     ! Skip to the last available position
     iflush=mdl%send_flush(ibuf)%array(iskip)
     iskip=iskip+(1+ndim+mdl%size_msg_array)*(iflush-1)+1

     ! Pack message header
     mdl%send_flush(ibuf)%array(iskip)=m%grid(igrid)%lev
     mdl%send_flush(ibuf)%array(iskip+1:iskip+ndim)=m%grid(igrid)%ckey(1:ndim)
     iskip=iskip+ndim+1

     ! Pack message content
     call pack_flush%proc(m,igrid,mdl%size_msg_array,mdl%send_flush(ibuf)%array(iskip:iskip+mdl%size_msg_array-1))

  endif
#endif

  end associate

end subroutine destage
!##############################################################
!##############################################################
!##############################################################
!##############################################################
subroutine destage_clump(mdl,local_peak_id)
  use amr_parameters, only: ndim
  use mdl_module
  use cache_commons
  use hash
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
  type(mdl_t)::mdl
  integer::local_peak_id
  !
  ! This routine frees the clump cache memory.
  ! It assembles flush messages, and when the message
  ! buffer is full, it sends it to the target CPU.
  !
  integer::info,icache,iflush,peak_cpu,ibuf
  integer::send_flush_id_clump,iskip,nflush
  integer::buffer_size_flush_array_clump,buffer_size_msg_array_clump
  integer(kind=8)::global_peak_id
  type(msg_unbind_clump)::dummy_large_clump

  associate(c=>mdl%c)

  ! Free hash table at global peak position
  icache = local_peak_id-c%npeak
  global_peak_id = c%gid(icache)

  if(hash_getp_simple(c%peak_dict,global_peak_id)==0)then
     write(*,*)'PE ',mdl_self(mdl),' trying to free non existing peak'
     stop
  endif

  call hash_free_simple(c%peak_dict,global_peak_id)

#ifndef WITHOUTMPI

  ! Check if the destage requires a flush
  if(c%dirty(icache))then
     peak_cpu=c%parent_cpu(icache)
     c%dirty(icache)=.false.

     ! Identify the corresponding communication buffer
     ibuf=mdl%cpu2buf_flush_clump(peak_cpu)

     ! If none is available, get the next free one
     if(ibuf==0)then
        mdl%ibuffer_flush_clump=mdl%ibuffer_flush_clump+1
        ! If there is no communication buffer left, allocate a new one
        if(mdl%ibuffer_flush_clump>mdl%nbuffer_flush_clump)then
           mdl%nbuffer_flush_clump=mdl%nbuffer_flush_clump+1
           buffer_size_msg_array_clump=storage_size(dummy_large_clump)/32 ! Largest message size
           buffer_size_flush_array_clump=1+(2+buffer_size_msg_array_clump)*nflushmax
           allocate(mdl%send_flush_clump(mdl%nbuffer_flush_clump)%array(1:buffer_size_flush_array_clump))
           mdl%send_flush_clump(mdl%nbuffer_flush_clump)%array(1)=0
        endif
        mdl%cpu2buf_flush_clump(peak_cpu)=mdl%ibuffer_flush_clump
        ibuf=mdl%ibuffer_flush_clump
     endif

     ! Filling the flush buffer
     iskip=1
     nflush=mdl%send_flush_clump(ibuf)%array(iskip)

     ! If buffer full, send it to remote CPU.
     if(nflush==nflushmax)then
        ! Post send
        call MPI_ISSEND(mdl%send_flush_clump(ibuf)%array(iskip),mdl%size_flush_array_clump,MPI_INTEGER,peak_cpu-1,flush_tag_clump,MPI_COMM_WORLD,send_flush_id_clump,info)
        ! While waiting for completion, check on incoming messages and perform actions
        call check_mail(mdl,send_flush_id_clump)
        ! Reset counter
        mdl%send_flush_clump(ibuf)%array(iskip)=0
     endif

     ! Increment counter
     mdl%send_flush_clump(ibuf)%array(iskip)=mdl%send_flush_clump(ibuf)%array(iskip)+1

     ! Skip to the last available position
     iflush=mdl%send_flush_clump(ibuf)%array(iskip)
     iskip=iskip+(2+mdl%size_msg_array_clump)*(iflush-1)+1

     ! Pack message header
     mdl%send_flush_clump(ibuf)%array(iskip:iskip+1)=transfer(global_peak_id,local_peak_id,2)
     iskip=iskip+2

     ! Pack message content
     call pack_flush_clump%proc(c,local_peak_id,mdl%size_msg_array_clump,mdl%send_flush_clump(ibuf)%array(iskip:iskip+mdl%size_msg_array_clump-1))

  endif
#endif

  end associate

end subroutine destage_clump
!##############################################################
!##############################################################
!##############################################################
!##############################################################
