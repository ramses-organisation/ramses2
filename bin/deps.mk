# Auto-generated Fortran module dependency rules for parallel make.
# Regenerate: python3 utils/py/gen_deps.py > bin/deps.mk

adaptive_loop.o: amr_step.o clump_finder.o init_amr.o init_cr.o init_hydro.o init_part.o init_refine_adaptive.o init_refine_basegrid.o init_refine_ramses.o init_refine_restart.o init_rt.o init_time.o init_xion.o input_part.o load_balance.o mdl.o ramses_commons.o read_params.o refine_utils.o turb_init.o update_time.o
amr_commons.o: amr_parameters.o cr_commons.o cr_parameters.o domain_hilbert.o hash.o hydro_commons.o hydro_parameters.o oct_commons.o rt_commons.o rt_parameters.o
amr_step.o: clump_finder.o cooling_fine.o core_step.o core_to_star_formation.o cr_godunov_fine.o cr_step.o feedback.o flag_utils.o godunov_fine.o interpol_phi.o lightcone.o move_fine.o movie.o newdt_fine.o output_amr.o pm_parameters.o ramses_commons.o refine_utils.o rho_fine.o rt_godunov_fine.o rt_step.o sink_evolution.o sink_formation.o sink_mechanical.o sink_merger.o source_hydro_fine.o star_formation.o synchro_hydro_fine.o tree_formation.o turb_driving.o turb_hydro.o update_time.o upload.o
boundana.o: amr_commons.o amr_parameters.o hydro_parameters.o
boundaries.o: amr_commons.o amr_parameters.o cr_parameters.o hydro_parameters.o ramses_commons.o rt_parameters.o
cache.o: amr_commons.o amr_parameters.o cache_commons.o clfind_commons.o hash.o hilbert.o mdl.o
cache_commons.o: amr_parameters.o call_back.o cr_parameters.o hydro_parameters.o rt_parameters.o
call_back.o: amr_commons.o amr_parameters.o clfind_commons.o ramses_commons.o
clfind_commons.o: amr_parameters.o hash.o
clump_finder.o: amr_parameters.o boundaries.o cache.o cache_commons.o clump_merger.o marshal.o mdl.o mdl_commons.o multigrid_fine_coarse.o nbors_utils.o output_amr.o output_clump.o output_part.o ramses_commons.o
clump_merger.o: amr_commons.o amr_parameters.o boundaries.o cache.o cache_commons.o clfind_commons.o hash.o hilbert.o marshal.o mdl.o mdl_commons.o nbors_utils.o pm_commons.o ramses_commons.o
condinit.o: amr_commons.o amr_parameters.o constants.o hydro_parameters.o input_hydro_condinit.o
cooling_fine.o: amr_commons.o amr_parameters.o constants.o cooling_module.o coolrates_module.o hydro_parameters.o mdl.o mdl_commons.o ramses_commons.o rt_parameters.o
cooling_module.o: constants.o
cooling_module_ism.o: constants.o
coolrates_module.o: amr_commons.o amr_parameters.o constants.o hydro_parameters.o rt_parameters.o
core_clump_output.o: amr_parameters.o mdl.o mdl_commons.o output_clump.o ramses_commons.o
core_evolution.o: amr_commons.o amr_parameters.o boundaries.o cache.o cache_commons.o constants.o hydro_parameters.o marshal.o mdl.o mdl_commons.o nbors_utils.o pm_commons.o ramses_commons.o sink_evolution.o
core_formation.o: amr_commons.o amr_parameters.o cache.o cache_commons.o clfind_commons.o clump_finder.o clump_merger.o constants.o hydro_parameters.o mdl.o mdl_commons.o pm_commons.o ramses_commons.o
core_merger.o: amr_commons.o amr_parameters.o cache.o constants.o mdl.o mdl_commons.o nbors_utils.o pm_commons.o ramses_commons.o
core_params.o: ramses_commons.o
core_retirement.o: amr_commons.o amr_parameters.o mdl.o mdl_commons.o pm_commons.o ramses_commons.o star_formation_common.o
core_sf_support.o: amr_commons.o
core_step.o: core_evolution.o core_formation.o core_merger.o core_retirement.o ramses_commons.o rho_fine.o
core_to_star_formation.o: amr_commons.o amr_parameters.o constants.o core_sf_support.o mdl.o mdl_commons.o pm_commons.o ramses_commons.o star_formation_common.o
courant_fine.o: amr_commons.o amr_parameters.o hydro_parameters.o mdl.o mdl_commons.o ramses_commons.o
cr_boundana.o: amr_commons.o amr_parameters.o cr_parameters.o
cr_commons.o: amr_parameters.o cr_parameters.o hydro_parameters.o oct_commons.o
cr_condinit.o: amr_commons.o amr_parameters.o cr_input_condinit.o cr_parameters.o
cr_flux_module.o: amr_commons.o amr_parameters.o cr_commons.o cr_parameters.o hydro_parameters.o
cr_godunov_fine.o: amr_commons.o amr_parameters.o boundaries.o cache.o cache_commons.o cr_commons.o cr_parameters.o hash.o hydro_parameters.o marshal.o mdl.o mdl_commons.o nbors_utils.o ramses_commons.o
cr_godunov_utils.o: amr_commons.o
cr_init_flow_fine.o: amr_commons.o amr_parameters.o constants.o cr_input_condinit.o cr_parameters.o mdl.o mdl_commons.o ramses_commons.o
cr_input_condinit.o: amr_commons.o amr_parameters.o cr_parameters.o mdl.o mdl_commons.o ramses_commons.o
cr_parameters.o: amr_parameters.o
cr_source_terms.o: amr_commons.o amr_parameters.o boundaries.o cache.o cache_commons.o cr_parameters.o hydro_parameters.o marshal.o mdl.o mdl_commons.o nbors_utils.o ramses_commons.o
cr_step.o: amr_parameters.o constants.o cr_godunov_fine.o cr_source_terms.o cr_upload.o mdl.o newdt_fine.o ramses_commons.o update_time.o
cr_upload.o: amr_commons.o amr_parameters.o cache.o cache_commons.o cr_parameters.o hydro_parameters.o mdl.o mdl_commons.o nbors_utils.o ramses_commons.o
cross_sections_module.o: constants.o
domain_hilbert.o: amr_parameters.o hilbert.o
feedback.o: amr_parameters.o boundaries.o cache.o cache_commons.o cr_parameters.o godunov_fine.o hilbert.o hydro_parameters.o marshal.o mdl.o mdl_commons.o nbors_utils.o pm_commons.o ramses_commons.o
flag_utils.o: amr_commons.o amr_parameters.o boundaries.o cache.o cache_commons.o hilbert.o hydro_flag.o marshal.o mdl.o mdl_commons.o nbors_utils.o pm_commons.o pm_parameters.o ramses_commons.o rt_flag.o smooth.o
force_fine.o: multigrid_fine_coarse.o
gadgetreadfile.o: amr_parameters.o
godunov_fine.o: amr_commons.o amr_parameters.o boundaries.o cache.o cache_commons.o hydro_commons.o hydro_parameters.o marshal.o mdl.o mdl_commons.o nbors_utils.o ramses_commons.o
godunov_utils.o: amr_commons.o amr_parameters.o hydro_parameters.o
grav_ana.o: amr_commons.o amr_parameters.o
hash.o: amr_parameters.o
hilbert.o: amr_parameters.o
hydro_commons.o: amr_parameters.o hydro_parameters.o rt_parameters.o
hydro_flag.o: amr_commons.o amr_parameters.o boundaries.o cache.o cache_commons.o hydro_parameters.o nbors_utils.o ramses_commons.o
hydro_parameters.o: amr_parameters.o
init_amr.o: amr_commons.o amr_parameters.o cr_parameters.o hash.o hilbert.o hydro_parameters.o mdl.o mdl_commons.o output_amr.o ramses_commons.o rt_parameters.o
init_cooling.o: amr_commons.o cooling_module.o
init_cr.o: amr_commons.o mdl.o mdl_commons.o ramses_commons.o
init_flow_fine.o: amr_parameters.o input_hydro_condinit.o input_hydro_gadget.o input_hydro_grafic.o ramses_commons.o
init_hydro.o: amr_commons.o mdl.o mdl_commons.o ramses_commons.o
init_neq_chem.o: amr_commons.o coolrates_module.o
init_part.o: amr_commons.o amr_parameters.o mdl.o mdl_commons.o pm_commons.o pm_parameters.o ramses_commons.o
init_refine_adaptive.o: cr_upload.o flag_utils.o init_part.o input_hydro_grafic.o input_part.o input_part_zoom.o ramses_commons.o refine_utils.o rt_upload.o upload.o
init_refine_basegrid.o: amr_parameters.o flag_utils.o hash.o hilbert.o input_hydro_grafic.o mdl.o mdl_commons.o ramses_commons.o
init_refine_ramses.o: amr_parameters.o hash.o hilbert.o hydro_parameters.o init_refine_basegrid.o input_hydro_condinit.o load_balance.o mdl.o mdl_commons.o output_amr.o ramses_commons.o read_params.o
init_refine_restart.o: amr_parameters.o cr_parameters.o hash.o hilbert.o hydro_parameters.o init_refine_basegrid.o load_balance.o mdl.o mdl_commons.o output_amr.o output_rt.o ramses_commons.o read_params.o rt_parameters.o
init_rt.o: amr_commons.o mdl.o mdl_commons.o ramses_commons.o
init_time.o: amr_commons.o amr_parameters.o cooling_module.o coolrates_module.o cr_init_flow_fine.o gadgetreadfile.o init_cooling.o init_neq_chem.o mdl.o mdl_commons.o ramses_commons.o rt_spectra.o
init_xion.o: amr_commons.o amr_parameters.o constants.o cooling_module.o coolrates_module.o hydro_parameters.o mdl.o mdl_commons.o ramses_commons.o upload.o
input_hydro_condinit.o: amr_commons.o amr_parameters.o hydro_parameters.o mdl.o mdl_commons.o ramses_commons.o
input_hydro_gadget.o: amr_parameters.o boundaries.o cache.o cache_commons.o godunov_fine.o hydro_parameters.o input_hydro_condinit.o marshal.o mdl.o mdl_commons.o nbors_utils.o ramses_commons.o
input_hydro_grafic.o: amr_commons.o amr_parameters.o hydro_parameters.o mdl.o mdl_commons.o ramses_commons.o
input_part.o: input_part_ascii.o input_part_gadget.o input_part_grafic.o input_part_ramses.o input_part_restart.o mdl.o mdl_commons.o ramses_commons.o rt_spectra.o
input_part_ascii.o: amr_commons.o amr_parameters.o mdl.o mdl_commons.o pm_commons.o ramses_commons.o
input_part_gadget.o: amr_commons.o amr_parameters.o gadgetreadfile.o init_part.o mdl.o mdl_commons.o output_amr.o pm_commons.o ramses_commons.o
input_part_grafic.o: amr_commons.o amr_parameters.o mdl.o mdl_commons.o pm_commons.o ramses_commons.o
input_part_ramses.o: amr_commons.o amr_parameters.o mdl.o mdl_commons.o output_amr.o pm_commons.o pm_parameters.o ramses_commons.o
input_part_restart.o: amr_commons.o amr_parameters.o mdl.o mdl_commons.o output_amr.o pm_commons.o pm_parameters.o ramses_commons.o
input_part_zoom.o: amr_commons.o amr_parameters.o mdl.o mdl_commons.o pm_commons.o ramses_commons.o
interpol_cr.o: amr_parameters.o cr_parameters.o
interpol_hydro.o: amr_parameters.o hydro_parameters.o
interpol_phi.o: amr_commons.o amr_parameters.o mdl.o mdl_commons.o ramses_commons.o
interpol_rt.o: amr_parameters.o rt_parameters.o
lightcone.o: amr_parameters.o lightcone_buffer.o lightcone_io.o lightcone_utils.o mdl.o mdl_commons.o pm_commons.o ramses_commons.o
lightcone_buffer.o: amr_parameters.o
lightcone_io.o: amr_parameters.o lightcone_buffer.o ramses_commons.o
lightcone_test.o: lightcone_io.o lightcone_utils.o ramses_commons.o
lightcone_utils.o: amr_commons.o
load_balance.o: amr_commons.o amr_parameters.o cache.o cache_commons.o cr_parameters.o domain_hilbert.o hash.o hilbert.o hydro_parameters.o init_refine_basegrid.o marshal.o mdl.o mdl_commons.o nbors_utils.o oct_commons.o pm_commons.o pm_parameters.o ramses_commons.o rho_fine.o rt_parameters.o
marshal.o: amr_commons.o amr_parameters.o cache_commons.o cr_parameters.o hydro_parameters.o rt_parameters.o
mdl.o: amr_commons.o amr_parameters.o clfind_commons.o mdl_commons.o
move_fine.o: amr_commons.o amr_parameters.o boundaries.o cache.o cache_commons.o cr_parameters.o hydro_parameters.o mdl.o mdl_commons.o nbors_utils.o pm_commons.o pm_parameters.o ramses_commons.o rho_fine.o rngstream.o rt_parameters.o
movie.o: amr_commons.o amr_parameters.o hilbert.o hydro_parameters.o mdl.o mdl_commons.o output_amr.o pm_commons.o ramses_commons.o rt_parameters.o
nbors_utils.o: amr_commons.o amr_parameters.o boundaries.o cache_commons.o hash.o hilbert.o hydro_parameters.o mdl.o ramses_commons.o
neq_cooling_module.o: amr_commons.o amr_parameters.o constants.o coolrates_module.o hydro_parameters.o rt_parameters.o
newdt_fine.o: amr_commons.o amr_parameters.o constants.o courant_fine.o cr_init_flow_fine.o mdl.o mdl_commons.o pm_commons.o pm_parameters.o ramses_commons.o
oct_commons.o: amr_parameters.o
output_amr.o: amr_commons.o amr_parameters.o cooling_module.o gadgetreadfile.o hydro_commons.o mdl.o mdl_commons.o output_cr.o output_hydro.o output_part.o output_poisson.o output_rt.o pm_commons.o posix_io.o ramses_commons.o turb_commons.o
output_clump.o: amr_commons.o amr_parameters.o clfind_commons.o clump_merger.o mdl.o mdl_commons.o pm_commons.o posix_io.o ramses_commons.o
output_cr.o: amr_commons.o amr_parameters.o constants.o cr_parameters.o mdl.o mdl_commons.o posix_io.o ramses_commons.o
output_hydro.o: amr_commons.o amr_parameters.o hydro_parameters.o mdl.o mdl_commons.o posix_io.o ramses_commons.o
output_part.o: amr_commons.o amr_parameters.o mdl.o mdl_commons.o pm_commons.o posix_io.o ramses_commons.o
output_poisson.o: amr_commons.o amr_parameters.o hydro_parameters.o mdl.o mdl_commons.o posix_io.o ramses_commons.o
output_rt.o: amr_commons.o amr_parameters.o constants.o hydro_parameters.o mdl.o mdl_commons.o posix_io.o ramses_commons.o rt_parameters.o
phi_fine_cg.o: multigrid_fine_coarse.o
pm_commons.o: amr_parameters.o
poisson_flag.o: amr_commons.o amr_parameters.o hash.o hydro_parameters.o ramses_commons.o
posix_io.o: amr_parameters.o
ramses_commons.o: amr_commons.o amr_parameters.o clfind_commons.o cooling_module.o coolrates_module.o cr_parameters.o hydro_parameters.o mdl.o pm_commons.o rt_parameters.o rt_spectra.o turb_commons.o
read_cr_params.o: amr_parameters.o constants.o cr_parameters.o hydro_parameters.o mdl.o movie.o ramses_commons.o
read_params.o: amr_commons.o amr_parameters.o constants.o core_params.o cr_parameters.o hilbert.o hydro_parameters.o mdl.o mdl_commons.o movie.o ramses_commons.o read_cr_params.o read_rt_params.o rt_parameters.o
read_rt_params.o: amr_parameters.o constants.o hydro_parameters.o mdl.o movie.o ramses_commons.o rt_parameters.o
refine_utils.o: amr_commons.o amr_parameters.o boundaries.o cache.o cache_commons.o call_back.o cr_parameters.o hash.o hilbert.o hydro_parameters.o init_refine_basegrid.o load_balance.o marshal.o mdl.o mdl_commons.o nbors_utils.o oct_commons.o ramses_commons.o rt_parameters.o
rho_ana.o: amr_parameters.o
rho_fine.o: amr_commons.o amr_parameters.o cache.o cache_commons.o hilbert.o mdl.o mdl_commons.o nbors_utils.o pm_commons.o pm_parameters.o ramses_commons.o
rt_commons.o: amr_parameters.o oct_commons.o rt_parameters.o
rt_condinit.o: amr_commons.o amr_parameters.o rt_input_condinit.o rt_parameters.o
rt_flag.o: amr_commons.o amr_parameters.o boundaries.o cache.o cache_commons.o nbors_utils.o ramses_commons.o rt_parameters.o
rt_flux_module.o: amr_parameters.o rt_parameters.o
rt_godunov_fine.o: amr_commons.o amr_parameters.o boundaries.o cache.o cache_commons.o hash.o marshal.o mdl.o mdl_commons.o nbors_utils.o ramses_commons.o rt_commons.o rt_flux_module.o rt_parameters.o
rt_godunov_utils.o: amr_commons.o amr_parameters.o rt_parameters.o
rt_init_flow_fine.o: amr_commons.o amr_parameters.o coolrates_module.o mdl.o mdl_commons.o output_rt.o ramses_commons.o rt_input_condinit.o rt_parameters.o rt_spectra.o
rt_input_condinit.o: amr_commons.o amr_parameters.o mdl.o mdl_commons.o ramses_commons.o rt_parameters.o
rt_parameters.o: amr_parameters.o
rt_spectra.o: amr_commons.o constants.o cross_sections_module.o hydro_parameters.o pm_commons.o rt_parameters.o
rt_star_feedback.o: amr_commons.o amr_parameters.o cache.o cache_commons.o constants.o hilbert.o hydro_parameters.o mdl.o mdl_commons.o nbors_utils.o pm_commons.o ramses_commons.o rt_parameters.o rt_spectra.o
rt_step.o: amr_parameters.o cooling_fine.o mdl.o newdt_fine.o ramses_commons.o rt_godunov_fine.o rt_input_condinit.o rt_star_feedback.o rt_upload.o update_time.o
rt_units.o: amr_commons.o amr_parameters.o
rt_upload.o: amr_commons.o amr_parameters.o cache.o cache_commons.o mdl.o mdl_commons.o nbors_utils.o ramses_commons.o rt_flag.o rt_parameters.o
rtz_cooling_module.o: amr_commons.o amr_parameters.o charge_exchange_module.o collisional_ionization_module.o constants.o coolrates_module.o cosmic_ray_ionization_module.o dust_recombination_module.o hydro_parameters.o molecules_module.o photoionization_UVB_module.o recombination_module.o rt_parameters.o rtz_coolrates_module.o rtz_module.o
rtz_coolrates_module.o: amr_commons.o charge_exchange_module.o coolrates_module.o cosmic_ray_ionization_module.o molecules_module.o photoionization_UVB_module.o rt_parameters.o rtz_module.o
sink_evolution.o: amr_commons.o amr_parameters.o boundaries.o cache.o cache_commons.o constants.o godunov_fine.o hilbert.o hydro_parameters.o marshal.o mdl.o mdl_commons.o nbors_utils.o pm_commons.o ramses_commons.o read_params.o rho_fine.o
sink_formation.o: amr_commons.o amr_parameters.o cache.o cache_commons.o clfind_commons.o clump_finder.o clump_merger.o constants.o hydro_parameters.o mdl.o mdl_commons.o output_amr.o output_clump.o output_part.o output_poisson.o pm_commons.o ramses_commons.o rngstream.o
sink_mechanical.o: amr_parameters.o boundaries.o cache.o cache_commons.o constants.o godunov_fine.o hilbert.o hydro_parameters.o marshal.o mdl.o mdl_commons.o nbors_utils.o pm_commons.o ramses_commons.o read_params.o sink_evolution.o
sink_merger.o: amr_commons.o amr_parameters.o cache.o cache_commons.o constants.o marshal.o mdl.o mdl_commons.o nbors_utils.o pm_commons.o ramses_commons.o sink_evolution.o
smooth.o: amr_parameters.o boundaries.o cache.o cache_commons.o marshal.o mdl.o mdl_commons.o nbors_utils.o ramses_commons.o
source_hydro_fine.o: amr_parameters.o boundaries.o cache.o cache_commons.o hydro_flag.o hydro_parameters.o mdl.o mdl_commons.o nbors_utils.o ramses_commons.o
star_formation.o: amr_commons.o amr_parameters.o constants.o hydro_parameters.o input_hydro_condinit.o mdl.o mdl_commons.o pm_commons.o ramses_commons.o rngstream.o rt_spectra.o
star_formation_common.o: ramses_commons.o rt_spectra.o
synchro_hydro_fine.o: amr_commons.o amr_parameters.o mdl.o mdl_commons.o ramses_commons.o
task_manager.o: amr_commons.o amr_parameters.o cache.o cache_commons.o call_back.o clfind_commons.o clump_finder.o clump_merger.o cooling_fine.o core_clump_output.o core_evolution.o core_formation.o core_merger.o core_retirement.o core_to_star_formation.o courant_fine.o cr_godunov_fine.o cr_init_flow_fine.o cr_input_condinit.o cr_source_terms.o cr_upload.o feedback.o flag_utils.o godunov_fine.o hash.o hilbert.o init_amr.o init_cr.o init_hydro.o init_part.o init_refine_basegrid.o init_refine_ramses.o init_refine_restart.o init_rt.o init_time.o init_xion.o input_hydro_condinit.o input_hydro_gadget.o input_hydro_grafic.o input_part.o input_part_ascii.o input_part_gadget.o input_part_grafic.o input_part_ramses.o input_part_restart.o input_part_zoom.o interpol_phi.o lightcone.o load_balance.o mdl.o mdl_commons.o move_fine.o movie.o newdt_fine.o output_amr.o output_clump.o output_cr.o output_hydro.o output_part.o output_poisson.o output_rt.o ramses_commons.o read_params.o refine_utils.o rho_fine.o rt_godunov_fine.o rt_init_flow_fine.o rt_input_condinit.o rt_star_feedback.o rt_upload.o sink_evolution.o sink_formation.o sink_mechanical.o sink_merger.o smooth.o source_hydro_fine.o star_formation.o synchro_hydro_fine.o tree_formation.o turb_driving.o turb_hydro.o turb_init.o turb_update.o update_time.o upload.o
test_cooling.o: amr_commons.o amr_parameters.o constants.o cooling_module.o coolrates_module.o hydro_parameters.o init_neq_chem.o neq_cooling_module.o rt_parameters.o
timer.o: amr_parameters.o mdl.o
tree_formation.o: amr_commons.o amr_parameters.o cache.o cache_commons.o clfind_commons.o clump_finder.o clump_merger.o constants.o hydro_parameters.o mdl.o mdl_commons.o pm_commons.o ramses_commons.o rngstream.o
turb_commons.o: amr_commons.o amr_parameters.o
turb_driving.o: amr_commons.o amr_parameters.o mdl.o mdl_commons.o ramses_commons.o turb_commons.o
turb_hydro.o: amr_commons.o amr_parameters.o mdl.o mdl_commons.o ramses_commons.o
turb_init.o: amr_commons.o mdl.o mdl_commons.o ramses_commons.o turb_commons.o
turb_update.o: amr_commons.o mdl.o mdl_commons.o ramses_commons.o turb_commons.o
umuscl.o: amr_commons.o amr_parameters.o hydro_parameters.o
units.o: amr_commons.o amr_parameters.o constants.o
update_time.o: amr_parameters.o cr_init_flow_fine.o hash.o mdl.o mdl_commons.o ramses_commons.o rt_init_flow_fine.o turb_update.o
upload.o: amr_commons.o amr_parameters.o cache.o cache_commons.o hydro_parameters.o mdl.o mdl_commons.o nbors_utils.o ramses_commons.o
vecpotentialinit.o: amr_commons.o amr_parameters.o
write_gitinfo.o: amr_parameters.o
write_screen.o: amr_commons.o amr_parameters.o hydro_parameters.o

ifeq ($(COMPILER),METAL)
adaptive_loop.o: metal_runner.o
amr_step.o: metal_runner.o
cooling_fine.o: metal_runner.o
courant_fine.o: metal_runner.o
flag_utils.o: metal_runner.o
force_fine.o: metal_runner.o
godunov_fine.o: metal_runner.o
init_amr.o: metal_runner.o
init_part.o: metal_runner.o
interpol_phi.o: metal_runner.o
mdl.o: metal_interface.o
metal_runner.o: amr_parameters.o constants.o cooling_module.o hydro_parameters.o mdl.o mdl_commons.o metal_interface.o ramses_commons.o
move_fine.o: metal_runner.o
multigrid_fine_coarse.o: metal_runner.o
multigrid_fine_commons.o: metal_runner.o
newdt_fine.o: metal_runner.o
phi_fine_cg.o: metal_runner.o
refine_utils.o: metal_runner.o
rho_fine.o: metal_runner.o
smooth.o: metal_runner.o
synchro_hydro_fine.o: metal_runner.o
task_manager.o: metal_runner.o
timer.o: metal_interface.o
upload.o: metal_runner.o
endif

ifeq ($(COMPILER),NVHPC)
adaptive_loop.o: gpu_manager.o
amr_step.o: gpu_manager.o
clump_finder.o: gpu_clump.o gpu_manager.o
clump_merger.o: gpu_clump.o
cooling_fine.o: gpu_cooling.o gpu_neq_cooling.o gpu_runner.o
core_evolution.o: gpu_core.o
core_formation.o: gpu_phew.o
core_merger.o: gpu_core_merger.o
core_retirement.o: gpu_core_star.o
core_to_star_formation.o: gpu_core_star.o
courant_fine.o: gpu_runner.o
feedback.o: gpu_runner.o
flag_utils.o: gpu_runner.o gpu_utils.o
force_fine.o: gpu_runner.o
godunov_fine.o: gpu_runner.o
gpu_clump.o: amr_parameters.o constants.o gpu_part.o gpu_runner.o gpu_utils.o hydro_parameters.o oct_commons.o pm_commons.o ramses_commons.o
gpu_cooling.o: amr_parameters.o cooling_module.o hydro_parameters.o oct_commons.o
gpu_core.o: amr_parameters.o constants.o cub_module_radix_sort_f.o gpu_feedback.o gpu_runner.o gpu_sink.o gpu_utils.o hydro_parameters.o ramses_commons.o
gpu_core_merger.o: amr_parameters.o constants.o cub_inclusive_scan_f.o cub_module_radix_sort_f.o gpu_feedback.o gpu_runner.o gpu_utils.o ramses_commons.o
gpu_core_star.o: amr_parameters.o constants.o core_sf_support.o gpu_runner.o gpu_utils.o ramses_commons.o star_formation_common.o
gpu_feedback.o: amr_parameters.o gpu_utils.o hydro_parameters.o oct_commons.o
gpu_flag.o: amr_parameters.o gpu_feedback.o gpu_utils.o hydro_parameters.o oct_commons.o
gpu_hilbert.o: amr_parameters.o hilbert.o
gpu_hydro.o: amr_parameters.o gpu_utils.o hydro_parameters.o oct_commons.o
gpu_manager.o: amr_parameters.o gpu_refine.o gpu_runner.o gpu_utils.o hydro_parameters.o mdl.o mdl_commons.o ramses_commons.o rt_parameters.o turb_commons.o
gpu_mg.o: amr_parameters.o gpu_utils.o hydro_parameters.o oct_commons.o
gpu_mpi.o: amr_parameters.o cache.o cache_commons.o gpu_hydro.o gpu_utils.o hash.o hydro_commons.o hydro_parameters.o marshal.o mdl.o nbors_utils.o oct_commons.o ramses_commons.o
gpu_nbor.o: amr_parameters.o gpu_utils.o oct_commons.o
gpu_neq_cooling.o: amr_commons.o amr_parameters.o constants.o coolrates_module.o gpu_runner.o gpu_utils.o hydro_parameters.o oct_commons.o ramses_commons.o rt_parameters.o
gpu_part.o: amr_parameters.o cub_module_radix_sort_f.o gpu_refine.o gpu_runner.o gpu_sink.o gpu_utils.o oct_commons.o ramses_commons.o
gpu_phew.o: amr_parameters.o constants.o cub_module_radix_sort_f.o gpu_clump.o gpu_part.o gpu_runner.o gpu_utils.o hydro_parameters.o oct_commons.o ramses_commons.o
gpu_refine.o: amr_parameters.o gpu_utils.o hydro_parameters.o oct_commons.o rt_parameters.o
gpu_rho.o: amr_parameters.o gpu_utils.o hydro_parameters.o oct_commons.o
gpu_rt.o: amr_parameters.o gpu_utils.o hydro_parameters.o oct_commons.o rt_parameters.o
gpu_runner.o: amr_parameters.o cache.o cache_commons.o constants.o cub_inclusive_scan_f.o cub_module_radix_sort_f.o gpu_cooling.o gpu_feedback.o gpu_flag.o gpu_hydro.o gpu_mg.o gpu_refine.o gpu_rho.o gpu_rt.o gpu_star.o gpu_turb.o gpu_utils.o hydro_parameters.o oct_commons.o ramses_commons.o rt_parameters.o turb_commons.o
gpu_scan.o: amr_parameters.o
gpu_sink.o: amr_parameters.o constants.o gpu_feedback.o gpu_runner.o gpu_utils.o hydro_parameters.o oct_commons.o ramses_commons.o
gpu_star.o: amr_parameters.o gpu_utils.o hydro_parameters.o oct_commons.o
gpu_turb.o: amr_parameters.o hydro_parameters.o oct_commons.o turb_commons.o
gpu_utils.o: amr_parameters.o oct_commons.o
init_amr.o: gpu_clump.o gpu_manager.o
init_part.o: gpu_manager.o gpu_phew.o
interpol_phi.o: gpu_runner.o
move_fine.o: gpu_part.o
multigrid_fine_coarse.o: gpu_runner.o
multigrid_fine_commons.o: gpu_runner.o
newdt_fine.o: gpu_part.o
output_amr.o: gpu_manager.o
output_clump.o: gpu_clump.o
phi_fine_cg.o: gpu_runner.o
refine_utils.o: gpu_runner.o
rho_fine.o: gpu_part.o gpu_runner.o
rt_godunov_fine.o: gpu_runner.o
rt_init_flow_fine.o: gpu_neq_cooling.o
rt_input_condinit.o: gpu_runner.o
rt_upload.o: gpu_runner.o
sink_evolution.o: gpu_sink.o
sink_formation.o: gpu_clump.o
sink_mechanical.o: gpu_sink.o
sink_merger.o: gpu_sink.o
smooth.o: gpu_runner.o
source_hydro_fine.o: gpu_runner.o
star_formation.o: gpu_runner.o
synchro_hydro_fine.o: gpu_runner.o
task_manager.o: gpu_manager.o
tree_formation.o: gpu_clump.o
turb_driving.o: gpu_runner.o
turb_hydro.o: gpu_runner.o
turb_init.o: gpu_manager.o
turb_update.o: gpu_manager.o
upload.o: gpu_runner.o
endif

ifeq ($(DO_CR),1)
cr_godunov_fine.o: cr_flux_module.o
cr_source_terms.o: cr_flux_module.o
refine_utils.o: cr_parameters.o
endif

ifeq ($(DO_RT),1)
core_formation.o: constants.o rt_parameters.o
endif

ifeq ($(DO_RTZ),0)
cooling_fine.o: neq_cooling_module.o
init_neq_chem.o: neq_cooling_module.o
init_xion.o: amr_commons.o amr_parameters.o hydro_parameters.o neq_cooling_module.o ramses_commons.o rt_parameters.o
refine_utils.o: init_xion.o
rt_init_flow_fine.o: neq_cooling_module.o
endif

ifeq ($(DO_RTZ),1)
cooling_fine.o: rtz_cooling_module.o rtz_module.o
init_neq_chem.o: charge_exchange_module.o cosmic_ray_ionization_module.o photoionization_UVB_module.o rtz_cooling_module.o rtz_coolrates_module.o
init_xion.o: rtz_module.o
input_hydro_condinit.o: rtz_module.o
output_rt.o: rtz_module.o
read_params.o: rtz_module.o
read_rt_params.o: cross_sections_module.o rt_spectra.o
rt_init_flow_fine.o: rtz_cooling_module.o
rt_spectra.o: constants.o rt_parameters.o rtz_module.o
task_manager.o: rtz_module.o
endif

ifeq ($(GRAV),1)
amr_step.o: force_fine.o multigrid_fine_commons.o phi_fine_cg.o
clump_finder.o: rho_fine.o
core_formation.o: core_clump_output.o rho_fine.o
force_fine.o: amr_commons.o amr_parameters.o boundaries.o cache.o cache_commons.o interpol_phi.o mdl.o mdl_commons.o nbors_utils.o phi_fine_cg.o ramses_commons.o
init_refine_adaptive.o: rho_fine.o
init_refine_basegrid.o: rho_fine.o
init_refine_ramses.o: rho_fine.o
init_refine_restart.o: rho_fine.o
multigrid_fine_coarse.o: amr_commons.o amr_parameters.o boundaries.o cache.o cache_commons.o hash.o hilbert.o mdl.o mdl_commons.o nbors_utils.o ramses_commons.o
multigrid_fine_commons.o: amr_commons.o amr_parameters.o boundaries.o cache.o cache_commons.o hash.o hilbert.o interpol_phi.o mdl.o mdl_commons.o multigrid_fine_coarse.o nbors_utils.o phi_fine_cg.o poisson_parameters.o ramses_commons.o
phi_fine_cg.o: amr_commons.o amr_parameters.o boundaries.o cache.o cache_commons.o call_back.o interpol_phi.o mdl.o mdl_commons.o nbors_utils.o ramses_commons.o
rho_fine.o: amr_commons.o amr_parameters.o cache.o cache_commons.o hilbert.o mdl.o mdl_commons.o multigrid_fine_coarse.o nbors_utils.o pm_commons.o pm_parameters.o ramses_commons.o
sink_formation.o: rho_fine.o
task_manager.o: force_fine.o multigrid_fine_coarse.o multigrid_fine_commons.o phi_fine_cg.o rho_fine.o
tree_formation.o: rho_fine.o
endif

ifeq ($(MHD),1)
godunov_utils.o: amr_commons.o amr_parameters.o hydro_parameters.o
umuscl.o: amr_commons.o amr_parameters.o hydro_parameters.o
endif
