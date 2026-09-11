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
| `cg_levelmin`       | `integer`  | 999 | Minimum level from which the Conjugate Gradient solver is used in place of the Multigrid solver. |
