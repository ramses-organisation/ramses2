The namelist block `&CLUMP_PARAMS` is used to specify parameters controlling the clump finder in the code. Clump finding is performed before each snapshot is outputted to disk. If you set `clump_only=.true.` in the `&RUN_PARAMS` namelist block then the code will only perform clump finding and stop right after outputting the clump catalog to disk.

| Variable name      | Fortran type | Default value       | Description                                                                                  |
|:-------------------|:-------------|:--------------------|:---------------------------------------------------------------------------------------------|
| `clump_finder`     | `logical`    | `.false.`           | Activate the clump finder if set to .true.                                                   |
| `merger_tree`      | `logical`    | `.false.`           | Activate merger tree particles if set to .true.                                              |
| `clump_info`       | `logical`    | `.false.`           | Turn clump finder diagnostics on or off.                                                     |
| `output_clump`     | `logical`    | `.false.`           | Output clump information to disk as an ascii file.                                           |
| `output_peak_grid` | `logical`    | `.false.`           | Output AMR grid density, peak id and halo id fields.                                         |
| `output_peak_part` | `logical`    | `.false.`           | Output dark matter particle peak id and halo id.                                             |
| `output_peak_tree` | `logical`    | `.false.`           | Output merger tree particle peak id and halo id.                                             |
| `output_peak_star` | `logical`    | `.false.`           | Output star particle peak id and halo id.                                                    |
| `output_peak_sink` | `logical`    | `.false.`           | Output sink particle peak id and halo id.                                                    |
| `output_peak_core` | `logical`    | `.false.`           | Output core particle peak id and halo id (needs `core=.true.`).                              |
| `rho_type_clump`   | `integer`    | 1                   | Set the fluid used to compute the density field (0: all, 1: DM, 2: stars, 4: gas).         |
| `nstep_per_tree`   | `integer`    | 1                   | Call the clump finder for merger-tree particle formation every `nstep_per_tree` coarse steps. |
| `ntestmax`         | `integer`    | 0                   | GPU only: maximum number of leaf cells above `density_threshold`, allocated at startup (36 bytes each). 0 means 1.5×`ngridmax`. The run stops if exceeded. |
| `nclumpmax`        | `integer`    | 0                   | GPU only: maximum number of density peaks, allocated at startup (~240 bytes each, plus 16×`nbin` bytes with halo mass profiles). 0 means `ngridmax`/50. The run stops if exceeded. |
| `relevance_threshold` | `real`    | 2                   | Set the relevance (prominence) threshold to remove noisy peaks.                              |
| `density_threshold`| `real`       | -1                  | Set the density threshold in code units above which peaks are detected. Must be > 0 when `clump_finder` or `merger_tree` is set. |
| `saddle_threshold` | `real`       | -1                  | Set the saddle point density in code units above which peaks are merged into haloes.         |
| `mass_threshold`   | `real`       | 0                   | Set the mass threshold in code units below which small clumps are discarded.                 |
| `purity_threshold` | `real`       | -1                  | Set the required fraction of mass in high-res particles (for zooms). Usually set to 0.98.    |
