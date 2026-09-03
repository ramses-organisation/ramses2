
[1]: https://github.com/ramses-organisation/ramses
[2]: https://github.com/ramses-organisation/ramses2

## RAMSES II ##

The ramses2 repository contains a major evolution of the original [RAMSES code][1] called [RAMSES II][2]. The main new features are perfect level-by-level load-balancing of grid cells and particles. Linked lists and octrees are not used anymore. They have been replaced by flat arrays and hash tables. The code works on Nvidia GPU using Nvidia Fortran. The CPU only code handles MPI communication using a low memory-footprint software cache.

The main development branch is the `develop` branch (set as default). This branch is still work in progress and the code is rapidly evolving.

You can download the code by cloning the git repository using
```
$ git clone git@github.com:ramses-organisation/ramses2.git
```

If you want to contribute to ramses2, you can fork this repository. To bring changes back into the `develop` branch of ramses2, simply issue a pull request.

To compile and execute the standard test cases, please follow these steps.

1- Shock tube test in 1D:

```
$ cd bin
$ make clean
$ make NDIM=1 HYDRO=1
$ cd ..
$ bin/ramses1d namelist/tube1d.nml
```

2- Sedov explosion in 2D:

```
$ cd bin
$ make clean
$ make NDIM=2 HYDRO=1
$ cd ..
$ bin/ramses2d namelist/sedov2d.nml
```

3- Magnetic loop advection in 2D:

```
$ cd bin
$ make clean
$ make NDIM=2 HYDRO=1 MHD=1 INIT=LOOP
$ cd ..
$ bin/ramses2d namelist/loop.nml
```

4- Sedov explosion test in 3D:

```
$ cd bin
$ make clean
$ make NDIM=3 HYDRO=1
$ cd ..
$ bin/ramses3d namelist/sedov3d.nml
```

5- Cosmological N body simulation in 3D

```
$ cd bin
$ make clean
$ make NDIM=3 HYDRO=0 GRAV=1 UNITS=COSMO
$ cd ..
$ utils/scripts/load_cosmo_ic.sh
$ bin/ramses3d namelist/dmo.nml
```

6- Molecular core test in 3D:

```
$ cd bin
$ make clean
$ make NDIM=3 HYDRO=1 GRAV=1 UNITS=COEUR INIT=COEUR
$ cd ..
$ bin/ramses3d namelist/coeur.nml
```

You get the picture now ;-)

To visualize the 2D and 3D results, compile the map making executable in the utils/f90 directory.

```
$ cd utils/f90
$ gfortran amr2map.f90 -o amr2map
$ cd ../..
$ utils/f90/amr2map -inp output_00002 -out dens.map -typ 1
$ utils/py/map2img.py dens.map --log
```

In the molecular cloud collapse case, you can also explore the movie1 directory and use the python function directly on any of the maps in there.
If you have the MPI library properly installed on your system, you can repeat all the tests above using the MPI=1 compilation option.
