cmake_minimum_required(VERSION 3.20)
cmake_policy(VERSION 3.20...3.28)

include(CMakePrintHelpers)

# =============================================================================
# Compiler selection
# Set COMPILER=GNU|INTEL|NVHPC before the first cmake run (or via FC env var).
# Example: cmake -DCOMPILER=NVHPC -DCMAKE_Fortran_COMPILER=nvfortran ..
#      or: FC=nvfortran cmake -DCOMPILER=NVHPC ..
# Changing COMPILER after the initial configure requires a fresh build dir.
# =============================================================================
set(COMPILER "GNU" CACHE STRING "Compiler toolchain: GNU, INTEL, NVHPC, METAL")

# Map COMPILER to the right executable on a fresh configure only.
if(NOT CMAKE_Fortran_COMPILER)
    if(COMPILER STREQUAL "NVHPC")
        set(CMAKE_Fortran_COMPILER nvfortran)
    elseif(COMPILER STREQUAL "INTEL")
        set(CMAKE_Fortran_COMPILER ifort)
    endif()
endif()

# Register .cuf as a Fortran extension before project() so that CMake's
# dependency scanner tracks MODULE / USE relationships inside .cuf files.
# Without this, parallel builds (make -j / ninja) fail because gpu_runner.mod
# is not yet present when poisson/*.f90 files that USE it start compiling.
# Must come before project() / enable_language(Fortran); it replaces (not
# augments) the built-in list, so all default extensions must be included.
if(COMPILER STREQUAL "NVHPC")
    set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS
        f F fpp FPP f77 F77 f90 F90 for For FOR f95 F95 f03 F03 f08 F08
        cuf CUF)
endif()

project(ramses LANGUAGES Fortran)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "No build type selected, defaulting to RelWithDebInfo")
    set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()

# =============================================================================
# Simulation parameters
# =============================================================================
set(NDIM     3    CACHE STRING "Number of dimensions (1, 2, or 3)")
set(NVECTOR  32   CACHE STRING "SIMD vector length")
set(NHILBERT 1    CACHE STRING "Hilbert key length exponent")
set(NSUBGRID 1    CACHE STRING "CUDA subgrid size in octs per dimension (NVHPC: 1, 2 or 4; MHD: 1 or 2)")
set(NPRE     8    CACHE STRING "Floating-point precision in bytes")
if(COMPILER STREQUAL "METAL")
    # Apple GPU has no float64; force single precision regardless of user setting.
    set(NPRE 4 CACHE STRING "Floating-point precision in bytes" FORCE)
endif()
set(NENER    0    CACHE STRING "Number of non-thermal energy variables")
set(NMETAL   0    CACHE STRING "Number of metal species")
set(NPSCAL   0    CACHE STRING "Number of passive scalar tracers")
set(NION     0    CACHE STRING "Number of ionisation species (auto-computed when RTZ=ON)")
set(NRTGRP   1    CACHE STRING "Number of RT photon groups (used when DO_RT=ON)")
set(NCRGRP   1    CACHE STRING "Number of CR groups (used when DO_CR=ON)")
set(UNITS    ""   CACHE STRING "Unit system (e.g. COSMO, COEUR; empty = default)")
set(INIT     ""   CACHE STRING "Initial condition module (e.g. LOOP, COEUR; empty = default)")
set(PATCH    ""   CACHE STRING "Patch directories (colon-separated, relative to the source tree) whose files replace the same-name sources")

# RTZ per-species ion counts (only relevant when DO_RTZ=ON)
set(N_HYDROGEN_IONS  0  CACHE STRING "RTZ: number of hydrogen  ion states")
set(N_H2             0  CACHE STRING "RTZ: number of H2 states")
set(N_HELIUM_IONS    0  CACHE STRING "RTZ: number of helium    ion states")
set(N_CARBON_IONS    0  CACHE STRING "RTZ: number of carbon    ion states")
set(N_NITROGEN_IONS  0  CACHE STRING "RTZ: number of nitrogen  ion states")
set(N_OXYGEN_IONS    0  CACHE STRING "RTZ: number of oxygen    ion states")
set(N_NEON_IONS      0  CACHE STRING "RTZ: number of neon      ion states")
set(N_MAGNESIUM_IONS 0  CACHE STRING "RTZ: number of magnesium ion states")
set(N_SILICON_IONS   0  CACHE STRING "RTZ: number of silicon   ion states")
set(N_SULFUR_IONS    0  CACHE STRING "RTZ: number of sulfur    ion states")
set(N_IRON_IONS      0  CACHE STRING "RTZ: number of iron      ion states")

# Physics toggles
option(HYDRO  "Activate hydrodynamics solver"          ON)
option(GRAV   "Activate self-gravity / Poisson solver" OFF)
option(MHD    "Activate MHD solver"                    OFF)
option(DO_RT  "Activate radiative transfer"            OFF)
option(DO_RTZ "Activate RTZ (metal-line RT) solver"    OFF)
option(DO_CR  "Activate cosmic rays"                   OFF)
option(TURB   "Activate turbulence driving"            OFF)
option(TRCFLX "Activate tracer-flux advection"         OFF)

# Build options
option(DEBUG        "Enable debug compilation flags"               OFF)
option(MPI          "Enable MPI parallelism"                       OFF)
option(HDF5         "Enable HDF5 output"                           OFF)
option(POSIXIO      "Write bulk output data with POSIX pwrite"     ON)
option(CUDA_PROFILE "Enable NVTX profiling instrumentation"        OFF)
option(PAPER        "Enable PAPER timing mode (NVHPC only)"        OFF)
option(CUB_SORT_PART   "Use CUB radix sort for particles (NVHPC)"  ON)
option(CUB_SORT_REFINE "Use CUB radix sort for refine (NVHPC)"     ON)
option(CUB_SCAN_REFINE "Use CUB inclusive scan for refine (NVHPC)" ON)

set(CUDA_ARCH "sm_80" CACHE STRING "CUDA target architecture (e.g. sm_80, sm_90)")
set(FFTW_DIR  ""      CACHE PATH   "FFTW3 install prefix (used when TURB=ON; empty = system search)")

set(EXEC ramses${NDIM}d)

if(COMPILER STREQUAL "METAL" AND MPI)
    message(FATAL_ERROR "COMPILER=METAL does not support MPI: set MPI=OFF")
endif()

# =============================================================================
# Derived variables
# =============================================================================

# Compute NION from per-species counts when DO_RTZ is on
if(DO_RTZ)
    math(EXPR NION "${N_HYDROGEN_IONS} + ${N_H2} + ${N_HELIUM_IONS} + ${N_CARBON_IONS} + ${N_NITROGEN_IONS} + ${N_OXYGEN_IONS} + ${N_NEON_IONS} + ${N_MAGNESIUM_IONS} + ${N_SILICON_IONS} + ${N_SULFUR_IONS} + ${N_IRON_IONS}")
endif()

# CR adds NCRGRP energy variables and forces MHD
if(DO_CR)
    math(EXPR NENER "${NENER} + ${NCRGRP}")
endif()

# NVAR = 5 (hydro) + NENER + NPSCAL + NMETAL + NION
math(EXPR NVAR "5 + ${NENER} + ${NPSCAL} + ${NMETAL} + ${NION}")

cmake_print_variables(COMPILER NDIM NVAR NENER NION NRTGRP NCRGRP UNITS INIT HYDRO GRAV MHD DO_RT DO_RTZ DO_CR TURB)

# =============================================================================
# Preprocessor defines
# =============================================================================
set(defines
    NDIM=${NDIM}
    NVAR=${NVAR}
    NVECTOR=${NVECTOR}
    NPRE=${NPRE}
    NHILBERT=${NHILBERT}
)

if(NENER GREATER 0)
    list(APPEND defines NENER=${NENER})
endif()
if(NION GREATER 0)
    list(APPEND defines NION=${NION})
endif()
if(NOT NSUBGRID EQUAL 1)
    list(APPEND defines NSUBGRID=${NSUBGRID})
endif()
if(CUDA_ARCH MATCHES "^sm_([0-9]+)[af]?$")
    list(APPEND defines GPU_SM=${CMAKE_MATCH_1})
endif()
if(NOT UNITS STREQUAL "")
    list(APPEND defines UNITS=${UNITS})
endif()
if(NOT INIT STREQUAL "")
    list(APPEND defines INIT=${INIT})
endif()
if(HYDRO)
    list(APPEND defines HYDRO)
endif()
if(GRAV)
    list(APPEND defines GRAV)
endif()
if(MHD)
    list(APPEND defines MHD)
endif()
if(DO_RT)
    list(APPEND defines DO_RT NRTGRP=${NRTGRP})
endif()
if(DO_RTZ)
    list(APPEND defines DO_RTZ
        N_HYDROGEN_IONS=${N_HYDROGEN_IONS}
        N_H2=${N_H2}
        N_HELIUM_IONS=${N_HELIUM_IONS}
        N_CARBON_IONS=${N_CARBON_IONS}
        N_NITROGEN_IONS=${N_NITROGEN_IONS}
        N_OXYGEN_IONS=${N_OXYGEN_IONS}
        N_NEON_IONS=${N_NEON_IONS}
        N_MAGNESIUM_IONS=${N_MAGNESIUM_IONS}
        N_SILICON_IONS=${N_SILICON_IONS}
        N_SULFUR_IONS=${N_SULFUR_IONS}
        N_IRON_IONS=${N_IRON_IONS}
    )
endif()
if(DO_CR)
    list(APPEND defines MHD DO_CR NCRGRP=${NCRGRP})
endif()
if(TURB)
    list(APPEND defines TURB)
endif()
if(TRCFLX)
    list(APPEND defines TRCFLX)
endif()
if(NOT MPI)
    list(APPEND defines WITHOUTMPI)
endif()
if(POSIXIO)
    list(APPEND defines POSIXIO)
endif()
if(COMPILER STREQUAL "METAL")
    list(APPEND defines _METAL)
endif()

# =============================================================================
# Executable
# =============================================================================
add_executable(${EXEC} "")

# =============================================================================
# MPI
# =============================================================================
if(MPI)
    find_package(MPI REQUIRED Fortran)
    target_include_directories(${EXEC} PRIVATE ${MPI_Fortran_INCLUDE_PATH})
    target_link_libraries(${EXEC} PRIVATE ${MPI_Fortran_LIBRARIES})
    separate_arguments(mpi_compile_flags NATIVE_COMMAND "${MPI_Fortran_COMPILE_FLAGS}")
    target_compile_options(${EXEC} PRIVATE
        $<$<COMPILE_LANGUAGE:Fortran>:${mpi_compile_flags}>)
endif()

# =============================================================================
# FFTW (turbulence driving)
# =============================================================================
if(TURB)
    find_path(FFTW_INCLUDE_DIR fftw3.f HINTS ${FFTW_DIR} PATH_SUFFIXES include REQUIRED)
    find_library(FFTW_LIBRARY fftw3 HINTS ${FFTW_DIR} PATH_SUFFIXES lib lib64 REQUIRED)
    target_include_directories(${EXEC} PRIVATE ${FFTW_INCLUDE_DIR})
    target_link_libraries(${EXEC} PRIVATE ${FFTW_LIBRARY})
endif()

# =============================================================================
# Compiler flags
# =============================================================================
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(fort_flags -cpp -ffree-line-length-none -fimplicit-none -frecord-marker=4)
    if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
        list(APPEND fort_flags -fallow-argument-mismatch)
    endif()
    if(DEBUG AND COMPILER STREQUAL "METAL")
        # Same as bin/Makefile: no FPE traps or NaN init for Metal
        list(APPEND fort_flags -O0 -fcheck=bounds -g -fbacktrace -Wall)
    elseif(DEBUG)
        list(APPEND fort_flags
            -O0 -g -fcheck=bounds -fbacktrace
            -ffpe-trap=zero,overflow,invalid -finit-real=nan
            -Wuninitialized -Wall)
    else()
        list(APPEND fort_flags -O3)
    endif()
    if(APPLE)
        target_link_options(${EXEC} PRIVATE -Wl,-no_compact_unwind)
    endif()
    if(COMPILER STREQUAL "METAL")
        # -Wl,-framework,X keeps each flag as one token through gfortran to ld.
        target_link_options(${EXEC} PRIVATE
            -Wl,-framework,Metal
            -Wl,-framework,Foundation
            -lstdc++)
    endif()

elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(fort_flags -cpp -recursive -fp-model source)
    if(DEBUG)
        list(APPEND fort_flags
            -warn all -O0 -g -traceback -check bounds
            -fpe0 -ftrapuv -init=zero -init=snan -init=arrays)
    else()
        list(APPEND fort_flags -O3 -diag-disable=10448)
    endif()

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    list(APPEND defines CUB_SORT_PART CUB_SORT_REFINE CUB_SCAN_REFINE)
    if(NOT CUB_SORT_PART)
        list(REMOVE_ITEM defines CUB_SORT_PART)
    endif()
    if(NOT CUB_SORT_REFINE)
        list(REMOVE_ITEM defines CUB_SORT_REFINE)
    endif()
    if(NOT CUB_SCAN_REFINE)
        list(REMOVE_ITEM defines CUB_SCAN_REFINE)
    endif()
    if(PAPER)
        list(APPEND defines PAPER)
    endif()
    set(nvhpc_cudalibs nvtx,curand)
    if(TURB)
        set(nvhpc_cudalibs nvtx,curand,cufft)
    endif()
    set(fort_flags -cpp -gpu=${CUDA_ARCH} -cuda -cudalib=${nvhpc_cudalibs} -Mfreeform)
    if(DEBUG)
        list(APPEND fort_flags
            -O0 -traceback -g
            -gpu=debug,lineinfo -Mprof=dwarf
            -ffpe-trap=zero,overflow,invalid -Minit-real=snan)
    else()
        list(APPEND fort_flags -O4)
    endif()
    if(CUDA_PROFILE)
        list(APPEND fort_flags -gpu=lineinfo,ptxinfo -Minstrument -traceback -gopt)
    endif()
    # Link flags mirror LIBS in the Makefile
    set(nvhpc_link_flags -gpu=${CUDA_ARCH} -cuda -cudalib=${nvhpc_cudalibs} -lstdc++)
    if(CUDA_PROFILE)
        list(APPEND nvhpc_link_flags -gpu=lineinfo,ptxinfo -Minstrument -traceback -gopt)
    endif()
    target_link_options(${EXEC} PRIVATE ${nvhpc_link_flags})

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
    set(fort_flags -e ZR)
endif()

target_compile_definitions(${EXEC} PRIVATE ${defines})
target_compile_options(${EXEC} PRIVATE
    $<$<COMPILE_LANGUAGE:Fortran>:${fort_flags}>)
set_target_properties(${EXEC} PROPERTIES
    Fortran_PREPROCESS ON
    Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules"
)

# =============================================================================
# Generated source files (git info, build command, patch)
# =============================================================================
# Regenerated at every build, like write_gitinfo.o: FORCE in bin/Makefile
list(JOIN defines " -D" _build_defines)
add_custom_target(gitinfo ALL
    COMMAND ${CMAKE_SOURCE_DIR}/utils/scripts/cr_write_gitinfo.sh
            ${CMAKE_BINARY_DIR}/write_gitinfo.f90
            "${PATCH}" "cmake -D${_build_defines}"
    BYPRODUCTS ${CMAKE_BINARY_DIR}/write_gitinfo.f90
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    VERBATIM
)
add_dependencies(${EXEC} gitinfo)
add_custom_command(
    OUTPUT  ${CMAKE_BINARY_DIR}/write_makefile.f90
    COMMAND ${CMAKE_SOURCE_DIR}/utils/scripts/cr_write_makefile.sh
            ${CMAKE_SOURCE_DIR}/CMakeLists.txt
            ${CMAKE_BINARY_DIR}/write_makefile.f90
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    VERBATIM
)
add_custom_command(
    OUTPUT  ${CMAKE_BINARY_DIR}/write_patch.f90
    COMMAND ${CMAKE_SOURCE_DIR}/utils/scripts/cr_write_patch.sh
            "" ${CMAKE_BINARY_DIR}/write_patch.f90
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    VERBATIM
)

# =============================================================================
# Fortran source files
# =============================================================================
set(fortran_files
  # AMR core
  amr/amr_parameters.f90
  amr/amr_commons.f90
  amr/oct_commons.f90
  amr/mdl_commons.f90
  amr/ramses_commons.f90
  amr/cache_commons.f90
  amr/hilbert.f90
  amr/domain_hilbert.f90
  amr/hash.f90
  amr/rngstream.f90
  amr/call_back.f90
  amr/boundaries.f90
  amr/nbors_utils.f90
  amr/flag_utils.f90
  amr/refine_utils.f90
  amr/smooth.f90
  amr/sort.f90
  amr/marshal.f90
  amr/cache.f90
  amr/timer.f90
  amr/title.f90
  amr/units.f90
  amr/posix_io.f90
  amr/output_amr.f90
  amr/movie.f90
  amr/load_balance.f90
  amr/read_params.f90
  amr/init_amr.f90
  amr/init_time.f90
  amr/init_refine_basegrid.f90
  amr/init_refine_adaptive.f90
  amr/init_refine_restart.f90
  amr/init_refine_ramses.f90
  amr/update_time.f90
  amr/amr_step.f90
  amr/adaptive_loop.f90
  amr/task_manager.f90
  amr/ramses.f90

  # Hydro
  hydro/hydro_parameters.f90
  hydro/hydro_commons.f90
  hydro/hydro_flag.f90
  hydro/boundana.f90
  hydro/godunov_utils.f90
  hydro/write_screen.f90
  hydro/output_hydro.f90
  hydro/interpol_hydro.f90
  hydro/init_hydro.f90
  hydro/init_flow_fine.f90
  hydro/input_hydro_condinit.f90
  hydro/input_hydro_grafic.f90
  hydro/input_hydro_gadget.f90
  hydro/courant_fine.f90
  hydro/umuscl.f90
  hydro/godunov_fine.f90
  hydro/synchro_hydro_fine.f90
  hydro/source_hydro_fine.f90
  hydro/upload.f90
  hydro/vecpotentialinit.f90
  hydro/condinit.f90

  # Particle mesh
  pm/pm_parameters.f90
  pm/pm_commons.f90
  pm/gadgetreadfile.f90
  pm/output_part.f90
  pm/init_part.f90
  pm/input_part.f90
  pm/input_part_restart.f90
  pm/input_part_ascii.f90
  pm/input_part_grafic.f90
  pm/input_part_gadget.f90
  pm/input_part_ramses.f90
  pm/input_part_zoom.f90
  pm/newdt_fine.f90
  pm/rho_fine.f90
  pm/move_fine.f90

  # Poisson / gravity
  poisson/poisson_parameters.f90
  poisson/multigrid_fine_commons.f90
  poisson/output_poisson.f90
  poisson/poisson_flag.f90
  poisson/rho_ana.f90
  poisson/grav_ana.f90
  poisson/interpol_phi.f90
  poisson/phi_fine_cg.f90
  poisson/force_fine.f90
  poisson/multigrid_fine_coarse.f90

  # Cooling
  cooling/constants.f90
  cooling/cooling_module.f90
  cooling/cooling_module_ism.f90
  cooling/coolrates_module.f90
  cooling/collisional_ionization_module.f90
  cooling/recombination_module.f90
  cooling/cosmic_ray_ionization_module.f90
  cooling/photoionization_UVB_module.f90
  cooling/dust_recombination_module.f90
  cooling/molecules_module.f90
  cooling/charge_exchange_module.f90
  cooling/rtz_module.f90
  cooling/init_cooling.f90
  cooling/init_xion.f90
  cooling/init_neq_chem.f90
  cooling/cooling_fine.f90

  # Cosmic rays
  cr/cr_parameters.f90
  cr/cr_commons.f90
  cr/cr_flux_module.f90
  cr/cr_godunov_utils.f90
  cr/cr_godunov_fine.f90
  cr/cr_source_terms.f90
  cr/cr_input_condinit.f90
  cr/cr_init_flow_fine.f90
  cr/cr_condinit.f90
  cr/cr_boundana.f90
  cr/cr_upload.f90
  cr/interpol_cr.f90
  cr/read_cr_params.f90
  cr/output_cr.f90
  cr/cr_step.f90
  cr/init_cr.f90

  # Radiative transfer
  rt/rt_parameters.f90
  rt/rt_commons.f90
  rt/cross_sections_module.f90
  rt/rt_spectra.f90
  rt/rt_flux_module.f90
  rt/rt_godunov_utils.f90
  rt/rt_flag.f90
  rt/rt_units.f90
  rt/output_rt.f90
  rt/interpol_rt.f90
  rt/rt_init_flow_fine.f90
  rt/rt_upload.f90
  rt/rt_condinit.f90
  rt/rt_input_condinit.f90
  rt/read_rt_params.f90
  rt/rt_godunov_fine.f90
  rt/rt_step.f90
  rt/init_rt.f90

  # Turbulence driving
  turb/turb_commons.f90
  turb/turb_hydro.f90
  turb/turb_driving.f90
  turb/turb_init.f90
  turb/turb_update.f90

  # Clump finder
  clump/clfind_commons.f90
  clump/output_clump.f90
  clump/tree_formation.f90
  clump/clump_merger.f90
  clump/clump_finder.f90

  # Stars, sinks, feedback
  star/star_formation.f90
  sink/sink_formation.f90
  sink/sink_evolution.f90
  sink/sink_mechanical.f90
  sink/sink_merger.f90
  feedback/feedback.f90
  feedback/rt_star_feedback.f90

  # Lightcone
  lightcone/lightcone_utils.f90
  lightcone/lightcone_buffer.f90
  lightcone/lightcone_io.f90
  lightcone/lightcone_test.f90
  lightcone/lightcone.f90

  # MDL
  mdl1/mdl.f90

  # Generated build-info stubs
  ${CMAKE_BINARY_DIR}/write_gitinfo.f90
  ${CMAKE_BINARY_DIR}/write_makefile.f90
  ${CMAKE_BINARY_DIR}/write_patch.f90
)

# DO_RTZ-conditional cooling modules (neq vs rtz variants)
if(DO_RTZ)
    list(APPEND fortran_files
        cooling/rtz_coolrates_module.f90
        cooling/rtz_cooling_module.f90
    )
else()
    list(APPEND fortran_files cooling/neq_cooling_module.f90)
endif()

# Core particles
list(APPEND fortran_files
    core/core_params.f90
    core/core_sf_support.f90
    core/star_formation_common.f90
    core/core_to_star_formation.f90
    core/core_clump_output.f90
    core/core_formation.f90
    core/core_evolution.f90
    core/core_retirement.f90
    core/core_merger.f90
    core/core_step.f90
)

# =============================================================================
# NVHPC / GPU sources
# =============================================================================
if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")

    # Regular CUDA Fortran sources – compiled with nvfortran + -cuda -gpu=...
    # Order matches bin/Makefile MODOBJ: gpu_runner must precede gpu_part.
    set(gpu_cuf_files
        gpu/gpu_utils.cuf
        gpu/gpu_hydro.cuf
        gpu/gpu_cooling.cuf
        gpu/gpu_star.cuf
        gpu/gpu_feedback.cuf
        gpu/gpu_rho.cuf
        gpu/gpu_mg.cuf
        gpu/gpu_turb.cuf

        gpu/gpu_refine.cuf
        gpu/gpu_flag.cuf
        gpu/gpu_rt.cuf
        gpu/gpu_runner.cuf
        gpu/gpu_neq_cooling.cuf
        gpu/gpu_clump.cuf
        gpu/gpu_sink.cuf
        gpu/gpu_part.cuf
        gpu/gpu_core.cuf
        gpu/gpu_core_merger.cuf
        gpu/gpu_core_star.cuf
        gpu/gpu_phew.cuf
        gpu/gpu_manager.cuf
    )
    list(APPEND fortran_files ${gpu_cuf_files})
    set_source_files_properties(${gpu_cuf_files} PROPERTIES LANGUAGE Fortran)

    # -------------------------------------------------------------------------
    # CUB wrappers: .cu files compiled with nvcc, then .cuf wrappers compiled
    # by nvfortran receiving the .o as input (for GPU device-code merging).
    # This mirrors the hand-written rules in bin/Makefile.
    # -------------------------------------------------------------------------
    set(cub_c_radix ${CMAKE_BINARY_DIR}/cub_module_radix_sort_c.o)
    set(cub_c_scan  ${CMAKE_BINARY_DIR}/cub_inclusive_scan_c.o)
    set(cub_f_radix ${CMAKE_BINARY_DIR}/cub_module_radix_sort_f.o)
    set(cub_f_scan  ${CMAKE_BINARY_DIR}/cub_inclusive_scan_f.o)


    add_custom_command(
        OUTPUT  ${cub_c_radix}
        COMMAND nvcc -arch=${CUDA_ARCH} -O3 -lineinfo
                -c ${CMAKE_SOURCE_DIR}/gpu/cub_module_radix_sort.cu
                -o ${cub_c_radix}
        DEPENDS ${CMAKE_SOURCE_DIR}/gpu/cub_module_radix_sort.cu
                ${CMAKE_SOURCE_DIR}/gpu/cub_module_radix_sort.h
        COMMENT "Compiling cub_module_radix_sort.cu (nvcc)"
        VERBATIM
    )
    add_custom_command(
        OUTPUT  ${cub_c_scan}
        COMMAND nvcc -arch=${CUDA_ARCH} -O3 -lineinfo
                -c ${CMAKE_SOURCE_DIR}/gpu/cub_inclusive_scan.cu
                -o ${cub_c_scan}
        DEPENDS ${CMAKE_SOURCE_DIR}/gpu/cub_inclusive_scan.cu
                ${CMAKE_SOURCE_DIR}/gpu/cub_inclusive_scan.h
        COMMENT "Compiling cub_inclusive_scan.cu (nvcc)"
        VERBATIM
    )
    add_custom_command(
        OUTPUT  ${cub_f_radix}
        COMMAND nvfortran ${fort_flags} -c
                ${cub_c_radix}
                ${CMAKE_SOURCE_DIR}/gpu/cub_module_radix_sort.cuf
                -o ${cub_f_radix}
        DEPENDS ${cub_c_radix}
                ${CMAKE_SOURCE_DIR}/gpu/cub_module_radix_sort.cuf
        BYPRODUCTS ${CMAKE_BINARY_DIR}/cub_module.mod
        COMMENT "Compiling cub_module_radix_sort.cuf (nvfortran)"
        VERBATIM
    )
    add_custom_command(
        OUTPUT  ${cub_f_scan}
        COMMAND nvfortran ${fort_flags} -c
                ${cub_c_scan}
                ${CMAKE_SOURCE_DIR}/gpu/cub_inclusive_scan.cuf
                -o ${cub_f_scan}
        DEPENDS ${cub_c_scan}
                ${CMAKE_SOURCE_DIR}/gpu/cub_inclusive_scan.cuf
        BYPRODUCTS ${CMAKE_BINARY_DIR}/cub_inclusive_scan_mod.mod
        COMMENT "Compiling cub_inclusive_scan.cuf (nvfortran)"
        VERBATIM
    )

    # Mark as pre-built external objects so CMake passes them straight to the linker
    set(cub_objects ${cub_c_radix} ${cub_c_scan} ${cub_f_radix} ${cub_f_scan})
    set_source_files_properties(${cub_objects}
        PROPERTIES EXTERNAL_OBJECT TRUE GENERATED TRUE)
    target_sources(${EXEC} PRIVATE ${cub_objects})

endif()

# =============================================================================
# Metal sources (Apple Silicon, NPRE=4, no MPI)
# =============================================================================
if(COMPILER STREQUAL "METAL")
    list(APPEND fortran_files
        metal/metal_interface.f90
        metal/metal_runner.f90)

    set(_metal_bridge_flags)
    if(MHD)
        list(APPEND _metal_bridge_flags -DMHD)
    endif()
    if(GRAV)
        list(APPEND _metal_bridge_flags -DGRAV)
    endif()
    if(DEBUG)
        list(APPEND _metal_bridge_flags -g -O0)
    endif()

    # ObjC++ bridge: compiled explicitly with xcrun clang++, same pattern as
    # the CUB wrappers in the NVHPC block — avoids relying on CMake's CXX
    # compiler detection (which may pick g++ instead of clang++ on macOS).
    set(_metal_bridge_o ${CMAKE_BINARY_DIR}/metal_bridge.o)
    add_custom_command(
        OUTPUT  ${_metal_bridge_o}
        COMMAND xcrun clang++ -std=c++17 -ObjC++ -fobjc-arc -O2
                -I${CMAKE_SOURCE_DIR}/metal
                ${_metal_bridge_flags}
                -c ${CMAKE_SOURCE_DIR}/metal/metal_bridge.mm
                -o ${_metal_bridge_o}
        DEPENDS ${CMAKE_SOURCE_DIR}/metal/metal_bridge.mm
                ${CMAKE_SOURCE_DIR}/metal/metal_types.h
                ${CMAKE_SOURCE_DIR}/metal/metal_config.h
        COMMENT "Compiling metal_bridge.mm (xcrun clang++)"
        VERBATIM
    )
    set_source_files_properties(${_metal_bridge_o}
        PROPERTIES EXTERNAL_OBJECT TRUE GENERATED TRUE)
    target_sources(${EXEC} PRIVATE ${_metal_bridge_o})

    # Metal kernel pipeline: .metal -> .air -> .metallib
    # NVAR is already a plain integer from math(EXPR) above.
    set(_metalc     xcrun -sdk macosx metal)
    set(_metallib   xcrun -sdk macosx metallib)
    set(_metalflags -std=metal3.1 -O2
        -DNDIM=${NDIM} -DNVAR=${NVAR}
        -I${CMAKE_SOURCE_DIR}/metal
        -I${CMAKE_SOURCE_DIR}/metal/kernels)
    if(HYDRO)
        list(APPEND _metalflags -DHYDRO)
    endif()
    if(GRAV)
        list(APPEND _metalflags -DGRAV)
    endif()
    if(MHD)
        list(APPEND _metalflags -DMHD)
    endif()
    if(DEBUG)
        list(APPEND _metalflags -gline-tables-only -frecord-sources)
    endif()

    set(_metal_kernels hydro refine flag scan rho mg cooling part)
    set(_metal_air)
    foreach(_kernel IN LISTS _metal_kernels)
        list(APPEND _metal_air ${CMAKE_BINARY_DIR}/${_kernel}.air)
        set(_kernel_depends ${CMAKE_SOURCE_DIR}/metal/kernels/${_kernel}.metal)
        if(MHD AND _kernel STREQUAL "hydro")
            list(APPEND _kernel_depends ${CMAKE_SOURCE_DIR}/metal/kernels/mhd_uct.metal)
        endif()
        add_custom_command(
            OUTPUT  ${CMAKE_BINARY_DIR}/${_kernel}.air
            COMMAND ${_metalc} ${_metalflags}
                    -c ${CMAKE_SOURCE_DIR}/metal/kernels/${_kernel}.metal
                    -o ${CMAKE_BINARY_DIR}/${_kernel}.air
            DEPENDS ${_kernel_depends}
                    ${CMAKE_SOURCE_DIR}/metal/metal_types.h
                    ${CMAKE_SOURCE_DIR}/metal/metal_config.h
                    ${CMAKE_SOURCE_DIR}/metal/kernels/metal_utils.h
            COMMENT "Compiling ${_kernel}.metal"
            VERBATIM
        )
    endforeach()

    add_custom_command(
        OUTPUT  ${CMAKE_BINARY_DIR}/ramses_kernels.metallib
        COMMAND ${_metallib} ${_metal_air}
                -o ${CMAKE_BINARY_DIR}/ramses_kernels.metallib
        DEPENDS ${_metal_air}
        COMMENT "Linking ramses_kernels.metallib"
        VERBATIM
    )

    add_custom_target(metal_kernels ALL
        DEPENDS ${CMAKE_BINARY_DIR}/ramses_kernels.metallib)
    add_dependencies(${EXEC} metal_kernels)
endif()

# =============================================================================
# Patch directories: a file in a patch directory replaces the source file of
# the same name, the first directory listed winning (like VPATH in bin/Makefile)
# =============================================================================
if(NOT PATCH STREQUAL "")
    string(REPLACE ":" ";" _patch_dirs "${PATCH}")
    set(_patched_files)
    foreach(_src IN LISTS fortran_files)
        set(_new "${_src}")
        if(NOT IS_ABSOLUTE "${_src}")
            get_filename_component(_name "${_src}" NAME)
            foreach(_dir IN LISTS _patch_dirs)
                get_filename_component(_absdir "${_dir}" ABSOLUTE BASE_DIR ${CMAKE_SOURCE_DIR})
                if(EXISTS "${_absdir}/${_name}")
                    file(RELATIVE_PATH _new ${CMAKE_SOURCE_DIR} "${_absdir}/${_name}")
                    message(STATUS "Patch: ${_src} -> ${_new}")
                    break()
                endif()
            endforeach()
        endif()
        list(APPEND _patched_files "${_new}")
    endforeach()
    set(fortran_files ${_patched_files})
endif()

# =============================================================================
# NVHPC: enforce GPU module compilation order via explicit OBJECT_DEPENDS.
#  2. Scan every non-GPU source file at configure time for "USE gpu_*".
#     Any file that USEs a GPU module is made to depend on the last GPU
#     object in the chain, which transitively ensures all GPU modules exist.
# =============================================================================
if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    set(_gpu_obj_dir "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${EXEC}.dir/gpu")

    # Step 1: GPU module dependencies (using object files)
    set_source_files_properties(gpu/gpu_hydro.cuf gpu/gpu_cooling.cuf gpu/gpu_star.cuf gpu/gpu_feedback.cuf gpu/gpu_rho.cuf
                                gpu/gpu_mg.cuf gpu/gpu_refine.cuf gpu/gpu_flag.cuf gpu/gpu_rt.cuf
        PROPERTIES OBJECT_DEPENDS "${_gpu_obj_dir}/gpu_utils.cuf.o")

    set_source_files_properties(gpu/gpu_runner.cuf
        PROPERTIES OBJECT_DEPENDS "${_gpu_obj_dir}/gpu_utils.cuf.o;${_gpu_obj_dir}/gpu_hydro.cuf.o;${_gpu_obj_dir}/gpu_cooling.cuf.o;${_gpu_obj_dir}/gpu_star.cuf.o;${_gpu_obj_dir}/gpu_feedback.cuf.o;${_gpu_obj_dir}/gpu_rho.cuf.o;${_gpu_obj_dir}/gpu_mg.cuf.o;${_gpu_obj_dir}/gpu_turb.cuf.o;${_gpu_obj_dir}/gpu_refine.cuf.o;${_gpu_obj_dir}/gpu_flag.cuf.o;${_gpu_obj_dir}/gpu_rt.cuf.o")

    set_source_files_properties(gpu/gpu_neq_cooling.cuf
        PROPERTIES OBJECT_DEPENDS "${_gpu_obj_dir}/gpu_utils.cuf.o;${_gpu_obj_dir}/gpu_runner.cuf.o")

    set_source_files_properties(gpu/gpu_clump.cuf
        PROPERTIES OBJECT_DEPENDS "${_gpu_obj_dir}/gpu_utils.cuf.o;${_gpu_obj_dir}/gpu_runner.cuf.o")

    set_source_files_properties(gpu/gpu_sink.cuf
        PROPERTIES OBJECT_DEPENDS "${_gpu_obj_dir}/gpu_utils.cuf.o;${_gpu_obj_dir}/gpu_feedback.cuf.o;${_gpu_obj_dir}/gpu_runner.cuf.o")

    set_source_files_properties(gpu/gpu_part.cuf
        PROPERTIES OBJECT_DEPENDS "${_gpu_obj_dir}/gpu_utils.cuf.o;${_gpu_obj_dir}/gpu_runner.cuf.o")

    set_source_files_properties(gpu/gpu_manager.cuf
        PROPERTIES OBJECT_DEPENDS "${_gpu_obj_dir}/gpu_runner.cuf.o")

    set_source_files_properties(gpu/gpu_core.cuf
        PROPERTIES OBJECT_DEPENDS "${_gpu_obj_dir}/gpu_utils.cuf.o;${_gpu_obj_dir}/gpu_feedback.cuf.o;${_gpu_obj_dir}/gpu_runner.cuf.o;${_gpu_obj_dir}/gpu_sink.cuf.o")

    set_source_files_properties(gpu/gpu_core_merger.cuf
        PROPERTIES OBJECT_DEPENDS "${_gpu_obj_dir}/gpu_utils.cuf.o;${_gpu_obj_dir}/gpu_feedback.cuf.o;${_gpu_obj_dir}/gpu_runner.cuf.o")

    set_source_files_properties(gpu/gpu_core_star.cuf
        PROPERTIES OBJECT_DEPENDS "${_gpu_obj_dir}/gpu_utils.cuf.o;${_gpu_obj_dir}/gpu_runner.cuf.o")

    set_source_files_properties(gpu/gpu_phew.cuf
        PROPERTIES OBJECT_DEPENDS "${_gpu_obj_dir}/gpu_utils.cuf.o;${_gpu_obj_dir}/gpu_runner.cuf.o;${_gpu_obj_dir}/gpu_clump.cuf.o;${_gpu_obj_dir}/gpu_part.cuf.o")

    # Step 2 – scan non-GPU sources for GPU module USE statements
    foreach(_src IN LISTS fortran_files)
        if(IS_ABSOLUTE "${_src}")
            continue()  # skip generated files (absolute paths)
        endif()
        if("${_src}" MATCHES "^gpu/")
            continue()  # skip GPU sources themselves
        endif()
        set(_full "${CMAKE_SOURCE_DIR}/${_src}")
        if(NOT EXISTS "${_full}")
            continue()
        endif()

        set(_deps "")
        file(STRINGS "${_full}" _use_runner REGEX "[Uu][Ss][Ee][ \t]+gpu_runner([ \t,;]|$)")
        if(_use_runner)
            list(APPEND _deps "${_gpu_obj_dir}/gpu_runner.cuf.o")
        endif()

        file(STRINGS "${_full}" _use_part REGEX "[Uu][Ss][Ee][ \t]+part_device([ \t,;]|$)")
        if(_use_part)
            list(APPEND _deps "${_gpu_obj_dir}/gpu_part.cuf.o")
        endif()

        file(STRINGS "${_full}" _use_sink REGEX "[Uu][Ss][Ee][ \t]+sink_device([ \t,;]|$)")
        if(_use_sink)
            list(APPEND _deps "${_gpu_obj_dir}/gpu_sink.cuf.o")
        endif()

        file(STRINGS "${_full}" _use_utils REGEX "[Uu][Ss][Ee][ \t]+gpu_utils([ \t,;]|$)")
        if(_use_utils)
            list(APPEND _deps "${_gpu_obj_dir}/gpu_utils.cuf.o")
        endif()

        file(STRINGS "${_full}" _use_manager REGEX "[Uu][Ss][Ee][ \t]+gpu_manager([ \t,;]|$)")
        if(_use_manager)
            list(APPEND _deps "${_gpu_obj_dir}/gpu_manager.cuf.o")
        endif()

        file(STRINGS "${_full}" _use_clump REGEX "[Uu][Ss][Ee][ \t]+clump_device([ \t,;]|$)")
        if(_use_clump)
            list(APPEND _deps "${_gpu_obj_dir}/gpu_clump.cuf.o")
        endif()

        file(STRINGS "${_full}" _use_cooling REGEX "[Uu][Ss][Ee][ \t]+cooling_device([ \t,;]|$)")
        if(_use_cooling)
            list(APPEND _deps "${_gpu_obj_dir}/gpu_cooling.cuf.o")
        endif()

        file(STRINGS "${_full}" _use_neq_cooling REGEX "[Uu][Ss][Ee][ \t]+neq_cooling_device([ \t,;]|$)")
        if(_use_neq_cooling)
            list(APPEND _deps "${_gpu_obj_dir}/gpu_neq_cooling.cuf.o")
        endif()

        file(STRINGS "${_full}" _use_core REGEX "[Uu][Ss][Ee][ \t]+core_device([ \t,;]|$)")
        if(_use_core)
            list(APPEND _deps "${_gpu_obj_dir}/gpu_core.cuf.o")
        endif()

        file(STRINGS "${_full}" _use_core_merger REGEX "[Uu][Ss][Ee][ \t]+core_merger_device([ \t,;]|$)")
        if(_use_core_merger)
            list(APPEND _deps "${_gpu_obj_dir}/gpu_core_merger.cuf.o")
        endif()

        file(STRINGS "${_full}" _use_core_star REGEX "[Uu][Ss][Ee][ \t]+core_star_device([ \t,;]|$)")
        if(_use_core_star)
            list(APPEND _deps "${_gpu_obj_dir}/gpu_core_star.cuf.o")
        endif()

        file(STRINGS "${_full}" _use_phew REGEX "[Uu][Ss][Ee][ \t]+phew_device([ \t,;]|$)")
        if(_use_phew)
            list(APPEND _deps "${_gpu_obj_dir}/gpu_phew.cuf.o")
        endif()

        if(_deps)
            set_source_files_properties("${_src}"
                PROPERTIES OBJECT_DEPENDS "${_deps}")
        endif()
    endforeach()
endif()

# =============================================================================
# Metal: enforce module compilation order via OBJECT_DEPENDS.
# metal_runner uses metal_interface, so metal_interface.f90.o must exist first.
# Files in amr/ and hydro/ that USE metal_runner/metal_interface are handled
# the same way the NVHPC block above handles gpu_runner dependencies.
# =============================================================================
if(COMPILER STREQUAL "METAL")
    set(_metal_obj_dir "${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/${EXEC}.dir/metal")

    set_source_files_properties(metal/metal_runner.f90
        PROPERTIES OBJECT_DEPENDS "${_metal_obj_dir}/metal_interface.f90.o")

    foreach(_src IN LISTS fortran_files)
        if(IS_ABSOLUTE "${_src}")
            continue()
        endif()
        if("${_src}" MATCHES "^metal/")
            continue()
        endif()
        set(_full "${CMAKE_SOURCE_DIR}/${_src}")
        if(NOT EXISTS "${_full}")
            continue()
        endif()
        set(_deps "")
        file(STRINGS "${_full}" _use_runner REGEX "[Uu][Ss][Ee][ \t]+metal_runner([ \t,;]|$)")
        file(STRINGS "${_full}" _use_iface  REGEX "[Uu][Ss][Ee][ \t]+metal_interface([ \t,;]|$)")
        if(_use_runner)
            list(APPEND _deps "${_metal_obj_dir}/metal_runner.f90.o")
        endif()
        if(_use_iface)
            list(APPEND _deps "${_metal_obj_dir}/metal_interface.f90.o")
        endif()
        if(_deps)
            set_source_files_properties("${_src}"
                PROPERTIES OBJECT_DEPENDS "${_deps}")
        endif()
    endforeach()
endif()

# =============================================================================
# Attach all Fortran sources
# =============================================================================
target_sources(${EXEC} PRIVATE ${fortran_files})
