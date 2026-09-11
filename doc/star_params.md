The namelist block `&STAR_PARAMS` is used to specify parameters controlling the subgrid star formation model and star particle properties.

| Variable name | Fortran type | Default value  | Description      |
|:------------------- |:-------|:----- |:------------------------- |
| `star   `           | `logical`   | `.false.`  | Turn the star formation model on or off. |
| `nstarmax`          | `integer`   | 0          | Maximum number of star particles that can be allocated during the run within each MPI process. |
| `n_star`            | `real`      | 0.1        | Density threshold in [H/cc] above which star formation is triggered. |
| `T2_star`           | `real`      | 2e4        | Temperature ceiling in [K] below which star formation is allowed. |
| `m_star`            | `real`      | 1          | Star particle mass in units of `mass_sph`. |
| `eps_star`          | `real`      | 0.01       | Star formation efficiency in units of the local gas freefall time. |
| `sf_model`          | `integer`   | 1          | Subgrid star formation model recipe (1: standard Schmidt law with constant threshold `n_star`). |
| `seed`              | `int array` | (123,456,789,1,1,1) | Random number generator seed. |
