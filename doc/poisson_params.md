The namelist block `&POISSON_PARAMS` is used to specify runtime parameters for the gravitational acceleration. The `&RUN_PARAMS` parameter `poisson=.true.` will activate self-gravity. If `poisson=.false.` then gravity is considered an external force with constant acceleration specified in the `&HYDRO_PARAMS`namelist.

Two different Poisson solvers are available in RAMSES: conjugate gradient (CG) and multigrid (MG). Unlike the CG solver, MG has an initialization overhead cost (at every call of the solver), but is much more efficient on very big levels with few "holes".  MG is always used for `levelmin`. MG can also be used on refined levels in conjuction with CG. The parameter `cg_levelmin` selects the Poisson solver as follows:

* The coarse level at `levelmin` is solved with MG
* Refined levels with `level<cg_levelmin` are solved with MG
* Refined levels with `level>=cg_levelmin` are solved with CG

| Variable name | Fortran type | Default value  | Description      |
|:------------------- |:-------|:----- |:------------------------- |
| `gravity_type`      | `integer`  | 0  | Type of gravity solver. `gravity_type=0` stands for self-gravity. `gravity_type>0` implements an external acceleration set in file `grava_ana.f90`. Finally, `gravity_type<0` combines self-gravity and an external acceleration. |
| `gravity_test`      | `logical`  | .false.  | Add an external density field set in file `rho_ana.f90` to test the Poisson solvers. |
| `gravity_params`    | `real array` | 0.0   | Runtime parameters used in the analytical expression of the external acceleration set in `grav_ana.f90` or the analytical density set in `rho_ana.f90`. |
| `epsilon`           | `real`  | 1e-4  | Stopping criterion for the iterative Poisson solver: residual 2-norm should be lower than `epsilon` times the right hand side 2-norm. |
| `nvcycle`           | `integer` | -1  | Desired number of Multigrid V-cycles (if set to `-1`, iterative solver cycles until residual satisfies `epsilon`). |
| `mg_safe_mode_reset` | `integer` | 0  | Number of solves a level stays in Multigrid safe mode once it has switched to it, counting the solve that triggered the switch (`0` = never reset, `1` = only that solve). The state is saved in the restart files. |
| `ngridmax_mg`       | `integer` | 0  | Maximum number of octs per MPI rank for the multigrid solver. 0 means `ngridmax`/7. |
| `cg_levelmin`       | `integer`  | 999 | Minimum level from which the Conjugate Gradient solver is used in place of the Multigrid solver. |

## Particle mass deposition and force interpolation

Each particle family has its own scheme for depositing mass onto the grid (which builds the right hand side of the Poisson equation) and for interpolating the gravitational force back to the particles. The families are dark matter (`part_*`), stars (`star_*`), sinks (`sink_*`) and tree particles (`tree_*`). The schemes are 1 = Cloud-In-Cell (CIC), 2 = Triangular-Shaped Cloud (TSC) and 3 = Piecewise Cubic Spline (PCS). Other values switch off the deposition or the interpolation of that family.

| Variable name | Fortran type | Default value  | Description      |
|:------------------- |:-------|:----- |:------------------------- |
| `part_mass_deposition_scheme`     | `integer` | 1 | Mass deposition scheme for dark matter particles (1 = CIC, 2 = TSC, 3 = PCS). |
| `part_force_interpolation_scheme` | `integer` | 1 | Force interpolation scheme for dark matter particles (1 = CIC, 2 = TSC, 3 = PCS). |
| `part_dep_algo`                   | `integer` | 2 | GPU CIC deposit kernel for dark matter particles: 1 = large (unshifted 3^ndim stencil), 2 = medium (half-cell shifted 2^ndim stencil), 3 = small (one segmented scan per stencil offset over runs of particles that share a source cell). Used with CIC only: TSC and PCS use a direct per-particle deposit on the GPU. |
| `star_mass_deposition_scheme`     | `integer` | 1 | Mass deposition scheme for star particles (1 = CIC, 2 = TSC, 3 = PCS). |
| `star_force_interpolation_scheme` | `integer` | 1 | Force interpolation scheme for star particles (1 = CIC, 2 = TSC, 3 = PCS). |
| `star_dep_algo`                   | `integer` | 2 | GPU CIC deposit kernel for star particles (same values as `part_dep_algo`). |
| `sink_mass_deposition_scheme`     | `integer` | 1 | Mass deposition scheme for sink particles (1 = CIC, 2 = TSC, 3 = PCS). |
| `sink_force_interpolation_scheme` | `integer` | 1 | Force interpolation scheme for sink particles (1 = CIC, 2 = TSC, 3 = PCS). |
| `sink_dep_algo`                   | `integer` | 2 | GPU CIC deposit kernel for sink particles (same values as `part_dep_algo`). |
| `tree_mass_deposition_scheme`     | `integer` | 1 | Mass deposition scheme for tree particles. Read but currently unused: tree particles do not deposit mass. |
| `tree_force_interpolation_scheme` | `integer` | 1 | Force interpolation scheme for tree particles (1 = CIC, 2 = TSC, 3 = PCS). |
