## Run parameters

This block, called `&RUN_PARAMS`, contains the run global control
parameters.

| Variable name, syntax, default value | Fortran type  | Description               |
|:---------------------------- |:------------- |:------------------------- |
| `cosmo=.false.`              |  `logical`    | Activate cosmological supercomoving cooordinates and computes the expansion factor |
| `pic=.false.`                |  `logical`    | Activate Particle-In-Cell solver |
| `poisson=.false.`            |  `logical`    | Activate Poisson solver for self-gravity |
| `hydro=.false.`              |  `logical`    | Activate hydro or MHD solver. |
| `rt=.false.`                 |  `logical`    | Activate radiative transfer solver. |
| `cr=.false.`                 |  `logical`    | Activate cosmic rays solver. |
| `clump_only=.false.`         |  `logical`    | Run only the clump finder on the initial conditions and stops the simulation |
| `verbose=.false.`            |  `logical`    | Activate verbose mode |
| `debug=.false.`              |  `logical`    | Activate debug mode |
| `static_mesh=.false.`        |  `logical`    | Activate static mesh refinement |
| `static_gas=.false.`         |  `logical`    | Turn off hydro variables update |
| `nrestart=0`                 |  `integer`    | Backup file number from which the code loads checkpoint/restart data and resumes the simulation, The default value, zero, is for a fresh start from the beginning (time=0).   |
| `nstepmax=1000000`           |  `integer`    | Maximum number of coarse time step. This can be used to terminate the simulation after a fixed amount of main steps. |
| `ncontrol=1`                 |  `integer`    | Frequency of screen output for control lines (to stdout), usually redirected into a log file. |
| `nremap=10`                  |  `integer`    | Frequency of call, in units of coarse time steps, for the particle load balancing routine, for MPI runs only. The default value of 10 means every 10 coarse step. |
| `nsubcycle=2,2,2,2,2,2,`     |  `integer array`    | Number of fine level sub-cycling steps within one coarser level time step. Each value in the array corresponds to a given level of refinement, starting from the coarse grid at `levelmin` up to the finest level at `levelmax`. `nsubcycle(1)=1` means that `levelmin` and `levelmin+1` are synchronized. To enforce single time stepping across the whole AMR hierarchy, you need to set `nsubcycle=1,1,1,1,1,1,1,`|
| `nsuperoct=0`                |  `integer`    | Collect `2**nsuperoct` neighboring octs per dimensions to assemble a block of octs when calling the hydro solver (if possible). The maximum number is 5 and corresponds to a block of `32**ndim` octs. |
