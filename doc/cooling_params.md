The namelist block `&COOLING_PARAMS` is used to specify parameters controlling gas cooling.

| Variable name | Fortran type | Default value  | Description      |
|:------------------- |:-------|:----- |:------------------------- |
| `cooling`           | `logical`  | `.false.`  | Turn galaxy formation equilibrium cooling on or off. |
| `neq_chem`          | `logical`  | `.false.`  | Turn galaxy formation non-equilibrium cooling and chemistry on or off. |
| `cooling_ism`       | `logical`  | `.false.`  | Turn ISM equilibrium cooling on or off. |
| `metal`             | `logical`  | `.false.`  | Turn metal advection on or off. This adds metallicity as an additional passive scalar. You must set `NVAR>5` to account for this additional hydro variable. |
| `isothermal`        | `logical`  | `.false.`  | Turn on or off a density-dependent equation of state (EoS). This forces the gas temperature to follow a polytropic EoS whose parameters can be set using the variables defined below. |
| `eos_type`          | `int`      | 1  | Type of EoS used in case `isothermal=.true.`. `eos_type=1` correspnds to a strict isothermal EoS. `eos_type=2` corresponds to a polytropic EoS. `eos_type=3` corresponds to the sum of an isothermal and a polytropic EoS, `eos_type=4` corresponds to the max of an isothermal and a polytropic EoS. |
| `eos_nH`            | `real`     | 1d50    | Characteristic density (in units of [H/cc]) used in the polytropic EoS. |
| `eos_T2`            | `real`     | 10      | Temperature (T/mu in units of [K]) of the isothermal EoS or temperature of the polytropic EoS at the characteristic density. |
| `eos_index`         | `real`     | 1       | Power law density dependance of the polytropic EoS. `eos_index=1` is equivalent to the isothermal case `eos_type=1`. |
| `mu_mol`            | `real`     | 1.2195  | Mean molecular weight of the gas. Used only if `cooling_ism=.true.`. |
| `X_H`               | `real`     | 0.76    | Hydrogen mass fraction. |
| `Y_He`              | `real`     | 0.24    | Helium mass fraction. Try and use `X_H + Y_He = 1`. |
| `haardt_madau`      | `logical`  | .false. | Turn on or off the default UV background model in the code. |
| `self_shielding`    | `logical`  | .false. | Turn on or off a simple self-shielding model to reduce the UV radiation using the following attenuation factor: `exp(-nH/0.01)`. |
| `J21`               | `real`     | 0       | Nornalization of the UV background Spectral Energy Distribution at 13.6 eV in units of `10^21 erg`. |
| `a_spec`            | `real`     | -1      | Power law of the UV background Spectral Energy Distribution as `J(nu)=J21*(nu/13.6)^(a_spec)`. |
| `z_reion`           | `real`     | 8.5     | Reionization redshift that marks the starting eoph of the adopted UV background. |
| `z_ave`             | `real`     | 0       | Average metallicity (in solar units) used in the cooling function in case `metal=.false.`. If `metal=.true.`, it is used as the default value in the initial conditions for the corresponding passive scalar. Note that a metallicity of 1 in solar units corresponds to the passive scalar value `Z=0.02`. |
| `T2max`             | `real`     | 1d50    | Maximum temperature (T/mu in units of [K]) allowed in the cooling routine. If for some reason `T2>T2max`, then the code resets `T2=T2max` in the corresponding cell. |
| `neq_Tconst`        | `real`     | -1      | If positive, set the value of the constant gas temperature for non-equilibrium chemistry. |
| `is_init_xion`      | `logical`  |`.false.`| Turn on or off equilibrium abundances for non-equilibrium chemistry initial conditions. |
| `upload_equilibrium_x` | `logical`  |`.true.`| Set photoionisation equilibrium ionisation fractions in the parent cell when de-refining. This is to avoid unnatural emission rates which can occur when the children ionization fractions are averaged (relevant mostly for mock observations) |
| `cooling_uvb_delta` | `real`     | 0.05    | Redshift spacing $\Delta z$ for the UV background interpolation table. |
| `isHe`              | `logical`  |`.false.`| Turn on or off non-equilibrium chemistry for Helium. |
| `isH2`              | `logical`  |`.false.`| Turn on or off non-equilibrium chemistry for H2 molecule. |
| `rtz_cooling`       | `logical`  |`.false.`| Turn on or off the RTZ non-equilibrium multi-element cooling and chemistry solver (`DO_RTZ`). |
| `rtz_equilibrium_test` | `integer` | -1    | Equilibrium test mode for RTZ chemistry. |
| `rtz_include_collisional_ionization` | `logical` | `.true.` | Include collisional ionization in RTZ chemistry. |
| `rtz_include_photoionization` | `logical` | `.true.` | Include photoionization processes in RTZ chemistry. |
| `rtz_include_cosmic_ray_ionization` | `logical` | `.true.` | Include cosmic ray ionization in RTZ chemistry. |
| `rtz_include_charge_exchange` | `logical` | `.true.` | Include charge exchange reactions between ions in RTZ chemistry. |
| `rtz_include_dust_recombination` | `logical` | `.true.` | Include grain-assisted recombination in RTZ chemistry. |
| `rtz_include_HM12_UVB` | `logical` | `.true.` | Include Haardt & Madau (2012) UV background in RTZ chemistry. |
| `isH2_rtz`          | `logical`  |`.false.`| Track molecular hydrogen ($H_2$) in RTZ chemistry. |
| `rtz_UV_background_G0` | `real`  | 0.0     | Far-UV background field intensity in Habing ($G_0$) units. |
| `rtz_primary_cosmic_ray_ionization_rate` | `real` | 0.0 | Primary cosmic ray ionization rate $\zeta_{\mathrm{CR}}$ in [$\text{s}^{-1}$]. |
| `rtz_max_cool_timestep` | `real` | 1.0d11  | Maximum allowable cooling timestep for RTZ chemistry in seconds. |
| `rtz_eqm_min_its`   | `integer`  | 0       | Minimum number of iterations in the RTZ equilibrium solver. |
