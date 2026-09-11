
## Documentation ##

This directory contains the basic documentation for RAMSES II.

We describe the meaning and default values of the namelist variables in each of the following namelist blocks:

[Compilation parameters](./compilation_params.md): parameters used to compile the code using the `make` command.

[Global run parameters](./run_params.md): Main global parameters of the code controlling its execution.

[Output parameters](./output_params.md): Parameters controlling the output strategy.

[AMR grid parameters](./amr_params.md): Parameters controlling mesh generation and memory.

[Poisson solver parameters](./poisson_params.md): Parameters controlling the Poisson solver for self-gravity.

[Initial conditions parameters](./init_params.md): Parameters used to setup the initial conditions.

[Hydro and MHD solver parameters](./hydro_params.md): Parameters controlling the MUSCL-Hancock scheme for solving the Euler and ideal MHD equations.

[Mesh refinement parameters](./refine_params.md): Parameters controlling the refinement strategy of the AMR grid.

[Unit system parameters](./units_params.md): Parameters controlling the conversion between code units and cgs units.

[Boundary conditions parameters](./boundary_params.md): Parameters controlling the boundary conditions of the hydro solver.

[Cooling parameters](./cooling_params.md): Parameters controlling the cooling function used in the code.

[Radiative transfer parameters](./rt_params.md): Parameters controlling the M1 radiation solver.

[Radiation sources parameters](./rt_sources.md): Parameters controlling the sources of radiation for the M1 solver.

[Radiation groups parameters](./rt_groups.md): Parameters controlling the radiation energy groups for the M1 solver.

[Cosmic ray parameters](./cr_params.md): Parameters controlling the cosmic ray two-moment (M1) solver.

[Cosmic ray groups parameters](./cr_groups.md): Parameters controlling the per-group properties for cosmic ray multigroup transport.

[Star formation parameters](./star_params.md): Parameters controlling the subgrid star formation model used in the code.

[Supernovae feedback parameters](./feedback_params.md): Parameters controlling the supernova feedback model.

[Turbulence driving parameters](./turb_params.md): Parameters controlling turbulent driving and forcing.

[Gadget initial conditions parameters](./gadget_params.md): Parameters controlling the reading and scaling of Gadget IC files.

[Clump finder parameters](./clump_params.md): Parameters controlling the clump finder.

[Sink particle parameters](./sink_params.md): Parameters controlling sink particle formation, accretion and AGN feedback.

[Lightcone parameters](./lightcone_params.md): Parameters controlling the light cone output in cosmological simulations.

[Movie parameters](./movie_params.md): Parameters controlling the generation of movie frames during simulation runs.
