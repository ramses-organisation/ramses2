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
| `rho_type_clump`   | `integer`    | 1                   | Set the fluid used to compute the density field (1: DM, 2: stars, 3: sinks, 4: gas).         |
| `nsteps_per_tree`  | `integer`    | 1                   | Call the clump finder for merger-tree particle formation every `nsteps_per_tree` coarse steps. |
| `relevance_threshold` | `real`    | 2                   | Set the relevance (prominence) threshold to remove noisy peaks.                              |
| `density_threshold`| `real`       | -1                  | Set the density threshold in code units above which peaks are detected.                      |
| `saddle_threshold` | `real`       | -1                  | Set the saddle point density in code units above which peaks are merged into haloes.         |
| `mass_threshold`   | `real`       | 0                   | Set the mass threshold in code units below which small clumps are discarded.                 |
| `purity_threshold` | `real`       | -1                  | Set the required fraction of mass in high-res particles (for zooms). Usually set to 0.98.    |
