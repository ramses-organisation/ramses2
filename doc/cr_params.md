The namelist block `&CR_PARAMS` is used to specify parameters controlling the cosmic ray two-moment (M1) solver. Cosmic rays require the code to be compiled with `DO_CR=1`, which also turns on MHD, and activated at runtime with `cr=.true.` in `&RUN_PARAMS`.

| Variable name | Fortran type | Default value  | Description      |
|:------------------- |:-------|:----- |:------------------------- |
| `cr_advect`         | `logical`   | `.false.`  | Turn on and off cosmic ray advection. It is off by default and should be turned on when cosmic ray sources are present in the simulation. |
| `cr_streaming_diffusion` | `logical` | `.false.` | Turn on and off cosmic ray streaming diffusion. When enabled, CRs stream along magnetic field lines at a velocity related to the Alfvén speed. |
| `cr_streaming_heating`   | `logical` | `.false.` | Turn on and off the heating of the thermal gas by cosmic ray streaming. This represents the energy transfer from CRs to the gas via streaming instability damping. |
| `cr_cooling`        | `logical`   | `.false.`  | Turn on and off cosmic ray energy losses (cooling). |
| `cr_isotropic_pressure`  | `logical` | `.true.`  | If `.true.`, CR pressure is treated as isotropic. If `.false.`, CR pressure is anisotropic and only acts along the magnetic field direction. |
| `cr_reduced_flux_correction` | `logical` | `.false.` | If `.true.`, enforce that the CR flux magnitude never exceeds the free-streaming limit F < c × E, where c is the (reduced) speed of light and E the CR energy density. |
| `cr_c_fraction`     | `real`      | 1.0        | Fraction of the true speed of light used to define the adopted reduced speed of light for cosmic rays. Similar to the RT reduced speed of light approximation. |
| `cr_dmax`           | `real`      | 1.0d30     | Maximum CR streaming diffusion coefficient in CGS units (cm²/s). Used to cap the streaming diffusion coefficient. |
| `cr_nsubcycle`      | `integer`   | 1          | Maximum number of cosmic ray subcycles allowed per hydro step. Increasing this allows CR transport to take smaller timesteps than the hydro, which can improve stability for high diffusion coefficients. |
| `cr_courant_factor` | `real`      | 0.8        | Courant factor used to compute the CR transport timestep. A value smaller than 1 improves stability at the cost of smaller timesteps. |
| `cr_test_setup`     | `character(LEN=100)` | `'none'`  | Name of a predefined setup for standard CR test problems. When set to `'none'`, no special test configuration is applied. |
| `cr_v_alfven`       | `real`      | 0.0        | Manually prescribed Alfvén velocity for idealised test problems (in code units). When set to 0.0, the Alfvén speed is computed self-consistently from the local magnetic field and density. |

The thermal energy is the total energy minus the kinetic, magnetic and CR energies. Where the CR energy dominates and the gas expands, it can become negative. Running with `entropy=.true.` and `dual_energy` (see `&HYDRO_PARAMS`) prevents this.
