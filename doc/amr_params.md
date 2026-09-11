This sets of parameters, contained in the namelist block `&AMR_PARAMS`, controls the AMR grid global properties. Parameters specifying the refinement strategy are described in the namelist block `&REFINE_PARAMS` and are used only if `levelmax>levelmin`.

| Variable name, syntax, default value | Fortran type  | Description       |
|:---------------------------- |:------------- |:------------------------- |
| `levelmin=1`                 |  `integer`    | Minimum level of refinement. This parameter sets the size of the coarse (or base) grid to `nx=2**levelmin`.|
| `levelmax=1`                 |  `integer`    | Maximum level of refinement. If `levelmax=levelmin`, the simulation will be executed on a standad Cartesian grid of linear size `nx=2**levelmin`|
| `ngridmax=0`                 |  `integer`    | Maximum number of grids (or octs) that can be allocated during the run within each MPI process. |
| `npartmax=0`                 |  `integer`    | Maximum number of dark matter particles that can be allocated during the run within each MPI process. |
| `ncachemax=10000`            |  `integer`    | Maximum number of cache lines per MPI rank allocated for the software cache. |
| `nexpand=1`                  |  `integer`    | Number of times the mesh expansion is applied to the refinement map (see mesh smoothing).|
| `boxlen=1.0`                 |  `real`       | Logical box size in code units. It corresponds to the square box in which Cartesian indices as well as Hilbert indices are defined.|
