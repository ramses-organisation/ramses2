This sets of parameters, contained in the namelist block `&INIT_PARAMS`. This is used to set up the initial conditions.

| Variable name, syntax, default value | Fortran type | Description |
|:---------------------------- |:------------- |:------------------------- |
| `initfile=' ' ` | `80*char` | Directory where IC files are stored.
| `filetype='ascii'` | `20*char` | Type of initial conditions file for particles and gas cells. Possible choices are `ascii`, `gadget`, `grafic` or `grafic_zoom`. Default value is `ascii`. In this case, only particles are read in the ascii file called `ic_part`, `ic_star`or `ic_sink`. For gas cells, initial conditions are set using file `condinit.f90` via compilation time parameters or via a patch. If nothing is done to set user defined initial conditions, then the default setup strategy is basd on uniform regions as explained now. |
| `aexp_ini=10.0` | `real` | This parameter sets the starting expansion factor for cosmology runs only. Default value is read in the IC file. |
| `boxlen_ini=1.0` | `real` | Box size in comoving Mpc/h used in grafic initial conditions. |
| `omega_b=0.045` | `real` | This parameter sets the baryonic density parameter for cosmology runs.  |
| `omega_m=1.0`   | `real` | Matter density parameter at z=0 for cosmology runs. |
| `omega_l=0.0`   | `real` | Cosmological constant density parameter at z=0 for cosmology runs. |
| `h0=1.0`        | `real` | Hubble parameter in units of H_0 = 100 km/s/Mpc for cosmology runs. |
| `ic_scale_l=1.0` | `real` | Length rescaling factor for initial conditions. |
| `ic_scale_v=1.0` | `real` | Velocity rescaling factor for initial conditions. |
| `ic_scale_m=1.0` | `real` | Rescale the dark matter particle masses only for ascii input files. |
| `A_ave=0.0`     | `real` | Average magnetic field in the box (x component). |
| `B_ave=0.0`     | `real` | Average magnetic field in the box (y component). |
| `C_ave=0.0`     | `real` | Average magnetic field in the box (z component). |
| `nregion=1`     | `integer` | Number of independent regions in the computational box used to set up initial flow variables. |
| `region_type='square'`  | `10*char` | Geometry defining each region. `square` defines a generalized ellipsoidal shape, while `point` defines a delta function. |
| `x_center=0.0`  | `real arrays` | X coordinate of the center of each region. |
| `y_center=0.0`  | `real arrays` | Y coordinate of the center of each region. |
| `z_center=0.0`  | `real arrays` | Z coordinate of the center of each region. |
| `length_x=0.0`  | `real arrays` | Size in X direction of each region. Used only for `square` regions.|
| `length_y=0.0`  | `real arrays` | Size in Y direction of each region. Used only for `square` regions.|
| `length_z=0.0`  | `real arrays` | Size in Z direction of each region. Used only for `square` regions.|
| `exp_region=2.0`  | `real arrays` | Exponent defining the norm used to compute distances for the generalized ellipsoid. `exp_region=2` corresponds to a spheroid, `exp_region=1` to a diamond shape, `exp_region>=10` to a perfect square. |
| `d_region=0.0`  | `real arrays` | Density. For `point` regions this is used to define mass. |
| `u_region=0.0`  | `real arrays` | X velocity. For `point` regions this is used to define momentum. |
| `v_region=0.0`  | `real arrays` | Y velocity. For `point` regions this is used to define momentum. |
| `w_region=0.0`  | `real arrays` | Z velocity. For `point` regions this is used to define momentum. |
| `B_region=0.0`  | `real arrays` | Y component of the magnetic field (only in 1D). |
| `C_region=0.0`  | `real arrays` | Z component of the magnetic field (only in 1D or 2D). |
| `p_region=0.0`  | `real arrays` | Pressure. For `point` regions this is used to define energy. |
| `prad_region=0.0` | `real arrays` | Non-thermal pressure. For `point` regions this is used to define energy. This is used only if compilation option `NENER>0`. |
| `var_region=0.0` | `real arrays` | Passive scalars. This is used only if compilation option `NVAR>5+NENER`. |
| `rt_nregion=1`  | `integer` | Number of independent regions in the computational box used to set up initial flow variables. |
| `rt_region_type='square'`  | `10*char` | Geometry defining each region. `square` defines a generalized ellipsoidal shape, while `point` defines a delta function. |
| `rt_reg_x_center=0.0`  | `real arrays` | X coordinate of the center of each region. |
| `rt_reg_y_center=0.0`  | `real arrays` | Y coordinate of the center of each region. |
| `rt_reg_z_center=0.0`  | `real arrays` | Z coordinate of the center of each region. |
| `rt_reg_length_x=0.0`  | `real arrays` | Size in X direction of each region. Used only for `square` regions.|
| `rt_reg_length_y=0.0`  | `real arrays` | Size in Y direction of each region. Used only for `square` regions.|
| `rt_reg_length_z=0.0`  | `real arrays` | Size in Z direction of each region. Used only for `square` regions.|
| `rt_exp_region=2.0`  | `real arrays` | Exponent defining the norm used to compute distances for the generalized ellipsoid. `exp_region=2` corresponds to a spheroid, `exp_region=1` to a diamond shape, `exp_region>=10` to a perfect square. |
| `rt_reg_group=1`  | `integer` | Index of the radiation group in which the rt region emits photons. |
| `rt_n_region=0.0`  | `real arrays` | Photon number density. For `point` regions this is used to define luminosity in photon per sec. |
| `rt_u_region=0.0`  | `real arrays` | Photon reduced flux direction (x coordinate). |
| `rt_v_region=0.0`  | `real arrays` | Photon reduced flux direction (y coordinate). |
| `rt_w_region=0.0`  | `real arrays` | Photon reduced flux direction (z coordinate). |
| `cr_nregion=0`  | `integer` | Number of independent regions in the computational box used to set up initial cosmic ray flux variables. Only available if compiled with `CR=1`. |
| `cr_region_type='square'`  | `10*char` | Geometry defining each CR region. `square` defines a generalized ellipsoidal shape. |
| `cr_reg_x_center=0.0`  | `real arrays` | X coordinate of the center of each CR region. |
| `cr_reg_y_center=0.0`  | `real arrays` | Y coordinate of the center of each CR region. |
| `cr_reg_z_center=0.0`  | `real arrays` | Z coordinate of the center of each CR region. |
| `cr_reg_length_x=1.E10`  | `real arrays` | Size in X direction of each CR region. Used only for `square` regions. |
| `cr_reg_length_y=1.E10`  | `real arrays` | Size in Y direction of each CR region. Used only for `square` regions. |
| `cr_reg_length_z=1.E10`  | `real arrays` | Size in Z direction of each CR region. Used only for `square` regions. |
| `cr_exp_region=2.0`  | `real arrays` | Exponent defining the norm used to compute distances for the generalized ellipsoid. `cr_exp_region=2` corresponds to a spheroid, `cr_exp_region=1` to a diamond shape, `cr_exp_region>=10` to a perfect square. |
| `cr_reg_group=1`  | `integer` | Index of the cosmic ray group into which the CR region injects flux. |
| `cr_fx_region=0.0`  | `real arrays` | X component of the cosmic ray flux in the region. |
| `cr_fy_region=0.0`  | `real arrays` | Y component of the cosmic ray flux in the region. |
| `cr_fz_region=0.0`  | `real arrays` | Z component of the cosmic ray flux in the region. |
