This namelist block, called `&OUTPUT_PARAMS`, is used to set up the output strategy of data to disk.

| Variable name, syntax, default value | Fortran type | Description |
|:---------------------------- |:------------- |:------------------------- |
| `nfile=1     ` | `integer` | Number of files per family (amr, hydro, grav, part, star...) per snapshot. These files are stored in directories called `output_00001`, `output_00002`, etc. Default value is `nfile=1`, even if one uses a large number of processors. `nfile=-1` will set `nfile=ncpu`, which is the maximum allowed value. |
| `tend=10     ` | `real`    | End time of the simulation, in time code units. |
| `delta_tout=0` | `real`    | Frequency of outputs in time code units. |
| `aend=0`       | `real`    | End time of the simulation, in expansion factor (for cosmology runs only). |
| `delta_aout=0` | `real`    | Frequency of outputs in expansion factor (for cosmology runs only). |
| `tout=0.0,0.0,0.0,` | `real array` | Value of specified output times with `tout>0`. |
| `nlanding=0` | `integer` | When positive, the coarse time step is reduced over the last coarse steps before each output time `tout` (or output expansion factor `aout` for cosmology runs), so that the snapshot is written at `tout` (or `aout`) instead of overshooting it: within `nlanding` coarse steps of the output, the remaining time is split into equal steps, so each step is reduced by at most a fraction `1/nlanding`. Typical values are 5 to 10. Default value of 0 means this is not used. With frequent outputs, keep `nlanding` below the number of coarse steps between two outputs. |
| `aout=1.1,1.1,1.1,` | `real array` | Value of specified output expansion factor (for cosmology runs only) with `0<aout<=1` where `aout=1` means "present epoch" or "redshift zero". |
| `bkp_time_hrs=2` | `real` | Wallclock time in hours between two consecutive backup files. Backup files are used to perform a checkpoint/restart of a simulation that ended prematurely.|
| `bkp_modulo=0`   | `integer` | Backup files are written using a numbering scheme that cycles modulo the prescribed number. Default value of 0 means never. A value of 1 corresponds to only 1 backup file at all time, with the risk of loosing it if something goes wrong during the write. Preferred values are 2 or 3. |
| `run_time_hrs=0` | `real` | Esimated wallclock time of the current job in hours. This is used to dump a last backup file just before the end of the simulation. Default value of 0 means this is not used.|
| `bkp_last_min=10` | `real` | Dump the last backup file just before the end of the simulation using the prescribed time in minutes. |
| `foutput=N` | `integer` | Frequency of outputs in units of main time steps. If N=0 then no output files are generated. |
| `output_part=.true.` | `logical` | When set to `.false.`, prevents writing individual dark matter particle data to regular output dumps. All other types of particle data (stars, sinks, tree, tracers) will still be written to standard output files if they are requested elsewhere in the namelist file. All particle data will still be written to backup files. Default is `.true.` |
| `output_grav=.true.` | `logical` | When set to `.false.`, prevents writing gravity/potential data to regular output dumps. Gravity data will still be written to backup files. Default is `.true.` |
| `output_hydro=.true.` | `logical` | When set to `.false.`, prevents writing hydro data to regular output dumps. Hydro data will still be written to backup files. Default is `.true.` |
| `output_amr=.true.` | `logical` | When set to `.false.`, prevents writing AMR grid structure data to regular output dumps. AMR data will still be written to backup files. Default is `.true.` |
