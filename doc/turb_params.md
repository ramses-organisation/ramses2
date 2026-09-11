The namelist block `&TURB_PARAMS` is used to specify parameters controlling the continuous turbulent driving (forcing) module in RAMSES II.

| Variable name | Fortran type | Default value | Description |
|:---|:---|:---|:---|
| `turb` | `logical` | `.false.` | Activate turbulent driving when set to `.true.`. |
| `turb_seed` | `integer` | `-1` | Random number generator seed for generating the stochastic turbulent acceleration field (`-1` uses a random seed). |
| `forcing_power_spectrum` | `character(100)` | `'parabolic'` | Power spectrum shape of the turbulent forcing. Available options are `'parabolic'` (parabolic decay from k=1 to k=3), `'power_law'` (slope -2), `'konstandin'` (linear decay from k=1 to k=2), or `'test'`. |
| `comp_frac` | `real` | `0.3333` | Compressive fraction for the Helmholtz decomposition of the forcing field (0.0 for purely solenoidal / divergence-free, 1.0 for purely compressive / curl-free). |
| `turb_T` | `real` | `1.0` | Autocorrelation time of the Ornstein-Uhlenbeck stochastic driving process in code time units. |
| `turb_Ndt` | `integer` | `100` | Number of sub-timesteps per autocorrelation time `turb_T`. |
| `turb_rms` | `real` | `1.0` | RMS amplitude of the turbulent forcing acceleration in code units. |
| `turb_min_rho` | `real` | `1.0d-50` | Minimum gas density threshold below which turbulent forcing is not applied. |
