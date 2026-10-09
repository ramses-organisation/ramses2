! Launch configuration of hydro_integrator_kernel in pure-hydro builds.
! Chosen by timing several thread counts and register limits for each NSUBGRID and precision
! (HLLC, uniform periodic 3D Sedov blast at 512^3, H200) and keeping the fastest; the same
! settings were then checked on an A100. launch_bounds(T,B) limits registers to 65536/(T*B).
! Double precision is register-limited, so a bound helps; single precision is fastest unbounded.
#ifndef HYDRO_LAUNCH_H
#define HYDRO_LAUNCH_H

#ifndef NSUBGRID
#define NSUBGRID 1
#endif

#ifndef MHD
#if defined(NPRE) && NPRE == 4
#define HYDRO_LB
#if NSUBGRID == 1
#define HYDRO_THREADS 64
#elif NSUBGRID == 2
#define HYDRO_THREADS 128
#else
#define HYDRO_THREADS 512
#endif
#else
#if NSUBGRID == 1
#define HYDRO_THREADS 64
#define HYDRO_LB launch_bounds(64,16)
#elif NSUBGRID == 2
#define HYDRO_THREADS 256
#define HYDRO_LB launch_bounds(256,4)
#else
#define HYDRO_THREADS 640
#define HYDRO_LB launch_bounds(640,1)
#endif
#endif
#endif

#endif
