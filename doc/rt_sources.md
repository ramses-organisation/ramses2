The namelist block `&RT_SOURCES` is used to specify the sources of radiation.

| Variable name | Fortran type | Default value  | Description      |
|:------------------- |:-------|:----- |:------------------------- |
| `rt_star`           | `logical`   | `.true.`   | Activate the emission of radiation from stellar particles. |
| `rt_sink`           | `logical`   | `.false.`  | Activate the emission of radiation from sink particles (this is for now not working). |
| `rt_esc_frac`       | `real`      | 1.0        | Escape fraction of photons from stellar particles, essentially a multiplication factor for the particle emission. |
| `rt_emission_stats` | `logical`   | `.false.`  | Print information about photon emission to standard output |
| `rt_nsource`        | `integer`   | 0          | Number of independent source (photon emission) regions in the computational box. |
| `rt_source_type`    | `character(len=128) array` |  `` | Geometry defining each source region. ’square’ defines a generalized ellipsoidal shape with photons injected everywhere inside, ’shell’ defines a finite width spherical shell into which photons are injected, and ’point’ defines a point source. |
| `rt_src_x_center`   | `real array` | 0.0        | X coordinates (0 to boxlen) of the center of each source region. |
| `rt_src_y_center`   | `real array` | 0.0        | Y coordinates (0 to boxlen) of the center of each source region. |
| `rt_src_z_center`   | `real array` | 0.0        | Z coordinates (0 to boxlen) of the center of each source region. |
| `rt_src_length_x`   | `real array` | 0.0        | X size of each source region. |
| `rt_src_length_y`   | `real array` | 0.0        | Y size of each source region. |
| `rt_src_length_z`   | `real array` | 0.0        | Z size of each source region. |
| `rt_exp_source`     | `real array` | 0.0        | Exponents defining the norm used to compute distances for the generalized ellipsoid. 2 corresponds to a spheroid, 1 to a diamond shape, 10 to a perfect square. |
| `rt_src_group`      | `integer`    | 1          | Photon groups into which photons are emitted in each source region (1 to M, where M is the number of groups). |
| `rt_n_source`       | `real array` | 0.0        | Photon injection for each source. For point region, this is a luminosity, photons per time, in code units. For square region, it is a flux in photons per time and area, also in code units.  |
| `rt_u_source`       | `real array` | 0.0        | Photon injection x-flux for each source. Betwen zero (for no flux in x-direction) and 1 (for full flux in x-diretion) |
| `rt_v_source`       | `real array` | 0.0        | Photon injection y-flux for each source. Betwen zero (for no flux in y-direction) and 1 (for full flux in y-diretion) |
| `rt_w_source`       | `real array` | 0.0        | Photon injection z-flux for each source. Betwen zero (for no flux in z-direction) and 1 (for full flux in z-diretion) |
