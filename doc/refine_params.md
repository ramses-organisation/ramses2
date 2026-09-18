The block named `&REFINE_PARAMS` contains the parameters related to grid refinement.

| Variable name | Fortran type | Default value  | Description               |
|:------------------- |:-------|:----- |:------------------------- |
| `interpol_var`      | `int`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  | 0     | Variables used to perform interpolation (prolongation) and averaging (restriction). `interpol_type=0`: conservatives; `interpol_type=1`: primitives |
| `interpol_type`     | `int`  | 1     | Type of slope limiter used in the interpolation scheme for newly refined cells. `interpol_type=0`: Straight injection (1st order), `interpol_type=1`: MinMod limiter, `interpol_type=2`: MonCen limiter, `interpol_type=3`: unlimited central slope. |
| `x_refine`          | `real array` | 0.0   | Geometry-based strategy: center of the refined region at each level of the AMR grid. |
| `y_refine`          | `real array` | 0.0   | Geometry-based strategy: center of the refined region at each level of the AMR grid. |
| `z_refine`          | `real array` | 0.0   | Geometry-based strategy: center of the refined region at each level of the AMR grid. |
| `r_refine`          | `real array` | 1e10  | Geometry-based strategy: radius of the refined region at each level. |
| `a_refine`          | `real array` | 1.0   | Geometry-based strategy: ratio Y/X of the refined region at each level. |
| `b_refine`          | `real array` | 1.0   | Geometry-based strategy: ratio Z/X of the refined region at each level. |
| `exp_refine`        | `real array` | 2.0   | Geometry-based strategy: exponent of the norm. |
| `err_grad_d`        | `real` | -1.0  | Discontinuity-based strategy: density gradient relative variations above which a cell is refined |
| `err_grad_u`        | `real` | -1.0  | Discontinuity-based strategy: velocity gradient relative variations above which a cell is refined |
| `err_grad_p`        | `real` | -1.0  | Discontinuity-based strategy: pressure gradient relative variations above which a cell is refined |
| `err_grad_b2`       | `real` | -1.0  | Discontinuity-based strategy: magnetic energy gradient relative variations above which a cell is refined |
| `err_grad_A`        | `real` | -1.0  | Discontinuity-based strategy: vector potential / magnetic field component A gradient threshold (`MHD=1`) |
| `err_grad_B`        | `real` | -1.0  | Discontinuity-based strategy: vector potential / magnetic field component B gradient threshold (`MHD=1`) |
| `err_grad_C`        | `real` | -1.0  | Discontinuity-based strategy: vector potential / magnetic field component C gradient threshold (`MHD=1`) |
| `err_grad_prad`     | `real array` | -1.0  | Discontinuity-based strategy: non-thermal energy gradient relative variations above which a cell is refined |
| `err_grad_var `     | `real array` | -1.0  | Discontinuity-based strategy: passive scalar gradient relative variations above which a cell is refined |
| `rt_err_grad_cn`     | `real array` | -1.0  | Discontinuity-based strategy: photon flux gradient relative variations above which a cell is refined |
| `err_grad_xHI`   | `real` | -1.0  | Discontinuity-based strategy: HI fraction gradient relative variations above which a cell is refined |
| `err_grad_xHII`  | `real` | -1.0  | Discontinuity-based strategy: HII fraction gradient relative variations above which a cell is refined |
| `floor_d`           | `real` | 1e-10 | Discontinuity-based strategy: density floor below which gradients are ignored |
| `floor_u`           | `real` | 1e-10 | Discontinuity-based strategy: velocity floor below which gradients are ignored |
| `floor_p`           | `real` | 1e-10 | Discontinuity-based strategy: pressure floor below which gradients are ignored |
| `floor_b2`          | `real` | 1e-10 | Discontinuity-based strategy: magnetic energy floor below which gradients are ignored |
| `floor_A`           | `real` | 1e-10 | Discontinuity-based strategy: magnetic field component A floor below which gradients are ignored (`MHD=1`) |
| `floor_B`           | `real` | 1e-10 | Discontinuity-based strategy: magnetic field component B floor below which gradients are ignored (`MHD=1`) |
| `floor_C`           | `real` | 1e-10 | Discontinuity-based strategy: magnetic field component C floor below which gradients are ignored (`MHD=1`) |
| `floor_prad`        | `real array` | 1e-10 | Discontinuity-based strategy: non-thermal pressure floor below which gradients are ignored |
| `rt_floor_cn`        | `real array` | -1.0  | Discontinuity-based strategy: photon group flux floor below which gradients are ignored |
| `rt_refine_aexp`    | `real` | -1.0  | Starting expansion factor for RT-based refinements in cosmological simulations |
| `floor_xHI`      | `real` | 1e-10 | Discontinuity-based strategy: HI fraction floor below which gradients are ignored |
| `floor_xHII`     | `real` | 1e-10 | Discontinuity-based strategy: HII fraction floor below which gradients are ignored |
| `jeans_refine`      | `real array` | -1.0   | Jeans refinement strategy: each level is refined if the cell size exceeds the local Jeans length divided by jeans_refine(ilevel). |
| `m_refine`          | `real array` | -1.0   | Quasi-Lagrangian strategy: each level is refined if the baryons mass in a cell exceeds `m_refine(ilevel)*mass_sph`, or if the number of dark matter particles exceeds `m_refine(ilevel)`, whatever the mass is. |
| `mass_sph`          | `real` | 0.0   | Quasi-Lagrangian strategy: `mass_sph` is used to set a typical baryonic mass scale. For cosmo runs, its value is set automatically to the initial mass resolution. |
| `mass_cut_refine`   | `real` | -1.0   | Quasi-Lagrangian strategy: mass threshold in code units for particle-based refinement. Particles more massive than this value are ignored. |
| `sink_refine`       | `logical` | `.false.` | Force maximum refinement level around sink particles |
| `ivar_refine`       | `int`  | -1    | Refinement map strategy: variable index (usually a passive scalar) used to define the refinement map. |
| `var_cut_refine`    | `real` | -1.0  | Refinement map strategy: threshold used on the refinement map to allow refinements. |
| `aexp_lock_refine`  | `real` | -1.0  | Activate progressive unlocking of levels when expansion factor doubles. This is used in cosmological simulations to enforce a quasi-constant resolution in physical units although length scales are in comoving units. The maximum level set by `levelmax` is activated after `aexp=aexp_lock_refine`. |
| `pic_lock_refine`   | `logical` | `.false.` | Lock particle-mesh refinement until `aexp >= aexp_lock_refine` |
