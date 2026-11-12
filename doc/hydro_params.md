This namelist is called &HYDRO_PARAMS, and is used to specify runtime parameters for the hydro and MHD solver. These parameters are based on rather standard concepts in computational fluid dynamics. We briefly describe them now.

| Variable name, syntax, default value | Fortran type  | Description               |
|:---------------------------- |:------------- |:------------------------- |
| `gamma=1.4`&nbsp;&nbsp;&nbsp;&nbsp;           |  `Real`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;    | Adiabatic exponent for the perfect gas EOS |
| `gamma_rad=1.333`    |  `Real array` | Adiabatic exponent for each non-thermal pressure EOS (used only if `NENER>1`) |
| `courant_factor=0.5` |  `Real`    | CFL number for time step control (less than 1) |
| `constant_gravity=0.0`|  `Real array` | 3 components of the constant gravitational acceleration in code units in case `poisson=.false.` |
| `smallr=1d-10 `      |  `Real`    | Minimum density to prevent floating exceptions |
| `smallc=1d-10 `      |  `Real`    | Minimum sound speed to prevent floating exceptions |
| `riemann=’llf’`      |  `Character LEN=20`| Name of the 1D Riemann solver. For the hydro solver (`MHD=0`), possible choices are `llf`, `hll` or `hllc`. For the MHD solver (`MHD=1`), possible choices are `llf`, `hll`, `roe`, `hlld`, `upwind`, and `uct-hlld` on GPUs (one MPI rank; on CUDA it takes about 2 KB of GPU memory per oct of `ngridmax`+`ncachemax` in double precision and 1 KB in single, more with passive scalars). |
| `riemann2d=’none’`   |  `Characher LEN=20`| Name of the 2D Riemann solver for the induction equation (`MHD=1` only). Possible choices are `llf`, `hll`, `roe`,`hlla`, `hllf`, `hlld` and `upwind`. If not set or set to `none`, the code will use the same as the 1D Riemann solver. Not used with `uct-hlld`. |
| `switch_llf_dmin=-1` |  `Real`| Minimum density in code units below which the adopted Riemann solvers (1D and 2D) will switch to the more diffusive LLF Riemann solver. |
| `switch_llf_pmin=-1` |  `Real`| Minimum pressure in code units below which the adopted Riemann solvers (1D and 2D) will switch to the more diffusive LLF Riemann solver. |
| `slope_type=1` |  `Integer`    | Type of slope limiter used for the piecewise linear reconstruction of volume-averaged quantities: `slope_type=0`: first order scheme, `slope_type=1`: MinMod limiter, `slope_type=2`: MonCen limiter, `slope_type=3`: Multi-dimensional MonCen limiter. In 1D runs only, it is also possible to choose: `slope_type=4`: Superbee limiter, `slope_type=5`: Ultrabee limiter. |
| `slope_mag_type=-1`  |  `Integer`    | Type of slope limiter used for the piecewise linear reconstruction of face-averaged quantities (`MHD=1` only): `slope_type=0`: first order scheme, `slope_type=1`: MinMod limiter, `slope_type=2`: MonCen limiter. If not set or set to `-1`, then the code will use `slope_mag_type=slope_type`.|
| `turb=.false.`       |  `Logical`    | Solve the turbulent kinetic energy equation using a standard LES model using an additional passive scalar. You must compile the code with `NVAR>5`. |
| `induction=.false.`  |  `Logical`    | Limit the MHD solver to the induction equation. The velocity field is set in the initial conditions but none of the hydro variables are updated. |
| `entropy=.false.`    |  `Logical`    | Solve for the conservation of entropy using an additional passive scalar. You must compile the code with `NVAR>5`. It can also be used in conjunction with the `dual_energy` parameter. |
| `dual_energy=-1`     |  `Real`    | Activate dual energy scheme for high-Mach flows if `dual_energy>=0`. Useful to prevent negative temperatures. The chosen value is used to set the fraction of the energy truncation error. Recommended values are between `0.0` and `0.5`. It must be used in conjunction with `entropy=.true.`. |
| `difmag=0d0` |  `Real`    | Add explicit diffusion for all volume-averaged conservative variable. |
| `etamag=0d0` |  `Real`    | Add explicit magnetic diffusivity (Ohm's law). `eta_ohm` is a synonym. |
| `eta_ad=0d0` |  `Real`    | Add ambipolar diffusion, with EMF `eta_ad*[(J.B)B-B^2 J]`. Needs `nimhd_solver='explicit'` or `'sts'`. |
| `nimhd_solver='unsplit'` |  `Character LEN=10`| Non-ideal MHD solver. `unsplit` adds `etamag` to the EMF of the Godunov step. `explicit` (substeps) and `sts` (super-time-stepping) solve Ohmic and ambipolar diffusion in a separate split step; with AMR, use `nsubcycle=1`. GPU, 3D and single rank only. The coefficients are set in `nimhd_coeff` (`gpu/gpu_hydro.cuf`). |
| `nimhd_split='strang'` |  `Character LEN=10`| Splitting of `explicit` and `sts`: `strang` diffuses half a time step before and half after the Godunov step (second order), `godunov` a full time step after it (first order). |
| `nimhd_courant=0.8` |  `Real`    | Safety factor of the explicit diffusion time step used to set the number of substeps or stages. |
| `nimhd_alpha=0.5` |  `Real`    | Gegenbauer parameter of `sts` (Mattia et al. 2026), at least 0.5. `0.5` gives the Runge-Kutta-Legendre scheme (Meyer et al. 2014). Larger values are more robust for anisotropic diffusion. |
