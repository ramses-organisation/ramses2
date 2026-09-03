We discuss here the various compilation time parameters. These parameters are used with the Makefile available in directory `bin/`.
You can compile the code by setting these parameters to your preffered value using:
`make PARAM=value`. You can see many examples in the directory `namelist/`. For example, to compile the code for a galaxy formation cosmological zoom simulation, go in directory `bin/` and type the commande `make NDIM=3 MPI=1 UNITS=COSMO HYDRO=1 GRAV=1 NPSCAL=3`.

| Variable, name, syntax, default value | Type | Description |
| -------- | ---- | ----------------- |
| `COMPILER=GNU`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;| `GNU`&nbsp;,&nbsp;`INTEL`,&nbsp;`NVHPC`&nbsp;or&nbsp;`METAL` | This sets the Fortran compiler you would like to use. The `GNU` option uses `gfortran` from the Gnu Compiler Collection, `INTEL` uses `ifort` from Intel, Inc., `NVHPC` uses the `nvfortran` compiler from the NVIDIA HPC SDK, and `METAL` uses `gfortran` with Apple Metal GPU kernels on macOS. |
| `NDIM = 3`       | `1`, `2` or `3`  | This sets the dimensionality of the simulation. |
| `HYDRO = 1`      | `0` or `1`  | This turns on or off the hydro solver. |
| `MHD = 0`        | `0` or `1`  | This turns on or off the MHD solver. |
| `GRAV = 0`       | `0` or `1`  | This turns on or off the self-gravity solver. |
| `DO_RT = 0`      | `0` or `1`  | This turns on or off the radiative transfer solver. |
| `DO_CR = 0`      | `0` or `1`  | This turns on or off the cosmic ray solver. Requires `MHD=1`. |
| `DO_RTZ = 0`     | `0` or `1`  | This turns on or off the radiative transfer solver with complex chemistry. |
| `MPI = 0`        | `0` or `1`  | This turns on or off parallel computing using the Message Passing Interface (MPI) library. MPI must be properly installed on your machine. In particular, you have to check that the MPI library was configured for your Fortran compiler. |
| `NPSCAL = 0`     | `integer`  | Number of passive scalars advected by the hydro solver. |
| `NENER = 0`      | `integer`  | Number of non-thermal energies used in the hydro solver. |
| `NMETAL = 0`     | `integer`  | Number of metal species advected by the hydro solver. |
| `NION = 0`       | `integer`  | Number of atomic ionized species advected by the hydro solver. |
| `NRTGRP = 1`     | `integer`  | Number of radiative transfer groups transported by the RT solver. |
| `NCRGRP = 1`     | `integer`  | Number of cosmic ray energy groups transported by the CR solver. |
| `INIT = `        | `string`  | This sets the adopted initial condition for your simulation. This parameter can be an empty string (default) wich uses the default initial condition generator in the code, as explained in the corresponding namelist block documentation. Possible values are `COEUR`, `INSTA`, `DOUBLEMACH`, `OT`, `PONO` and more. These initial conditions are set in the file `hydro/condinit.f90`. You can add your own here and share them with the community later via a pull request. |
| `UNITS = `       | `string`  | This sets the adopted unit system to convert variables from code units to cgs units. This parameter can be an empty string (default) wich uses the default unit system that can be set via the namelist, as explained in the corresponding namelist block documentation. You can also use predefined unit systems. Possible values are `COEUR`, `COSMO`, `MERGER` and more. These different unit systems are set in the file `amr/units.f90`. You can add your own here and share them with the community later via a pull request. |
| `MPIF90 = mpif90` | `string`  | This sets the name of the MPI Fortran compiler on your system. |
| `EXEC=ramses`    | `string` | This sets the name of the executable for the compiled code. The number of dimensions used at compilation via parameter `NDIM` is appended at the end. |
| `NHILBERT = 1`   | `1`, `2` or `3`  | This sets the number of long integers uswed to code the Hilbert key. This choice sets also the maximum level of refinement. In 3D, `NHILBERT=1` corresponds to `levelmax<21`, `NHILBERT=2` to `levelmax<42` and `NHILBERT=3` to `levelmax<63`.|
| `NVECTOR = 32`   | `integer`  | This sets the size of the vector sweeps used in many subroutines of the code. This is used for optimization purposes. This is highly problem and architecture dependant. |
| `NPRE = 8`       | `4` or `8`  | This sets the number of bytes used to code floating point numbers. `4` means single precision (not recommanded) `8` means double precision. |
| `DEBUG = 0`      | `0` or `1`  | This turns on or off the debug mode. |
