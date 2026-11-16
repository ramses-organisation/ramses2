The namelist block `&FEEDBACK_PARAMS` is used to specify parameters controlling the subgrid supernova (Type II) feedback model and stellar mass/metal return into the interstellar medium.

| Variable name | Fortran type | Default value | Description |
|:---|:---|:---|:---|
| `thermal_feedback` | `logical` | `.false.` | Turn on or off thermal supernova feedback (injects thermal energy into the host cell). |
| `mechanical_feedback` | `logical` | `.false.` | Turn on or off mechanical supernova feedback (injects momentum and kinetic/thermal energy using an outflow bubble model). |
| `mechanical_feedback_kimm` | `logical` | `.false.` | Turn on or off the cell-grouped mechanical supernova feedback of Taysun Kimm (one deposit per host cell, neighbor energy set from the old state plus all deposits, thermal dump). Cannot be combined with `mechanical_feedback`. CPU only. |
| `mechanical_r_thermal_dump` | `real` | `0.1` | For `mechanical_feedback_kimm`: when the mean ratio of the swept-up mass to the transition mass over the stencil is below this value, the supernova energy is injected as heat. Set to 0 to turn the thermal dump off. |
| `M_SNII` | `real` | `10.0` | Typical progenitor mass per Type II supernova in solar masses. |
| `E_SNII` | `real` | `1.0d51` | Energy released per Type II supernova in ergs. |
| `t_SNII` | `real` | `20.0` | Lifetime of supernova progenitor stars in Myr. Supernova explosions occur when star particles reach an age of `t_SNII`. |
| `eta_SNII` | `real` | `0.1` | Mass fraction of newly formed star particles returned to the gas as supernova ejecta (corresponding to a certain stellar initial mass function). |
| `yield_SNII` | `real` | `0.1` | Metal yield: fraction of the supernova ejecta mass released as new metals. |
