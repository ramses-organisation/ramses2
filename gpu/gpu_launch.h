! Launch configuration of the tuned CUDA kernels, included by gpu_hydro.cuf and gpu_runner.cuf.
!
! Each block sets, for one build, the threads per block (_THREADS) and the minimum resident blocks
! per SM (_MINB) of every tuned kernel. A kernel is compiled with launch_bounds(_THREADS,_MINB), or
! without a bound of its own when _MINB is 0, and launched with the _THREADS that hydro_device
! exports as parameters with the names of the macros in lower case. gpu_print_launch (gpu_runner.cuf)
! writes the values of the kernels a run uses at startup, with debug=.true.
!
! GPU_LAUNCH_TARGET is the GPU whose measured blocks a build uses (0: generic blocks). The B200
! blocks cover 3D builds without passive scalars, the H200 and A100 blocks cover MHD. Measured on
! Orszag-Tang at 512^3 (A100: 256^3) for MHD and on Sedov at 512^3 for hydro. To tune another GPU,
! add a target here and a copy of the generic blocks for it.
#if defined(GPU_SM) && GPU_SM == 100 && NDIM == 3 && NVAR == 5 && !defined(PAPER) \
    && (NSUBGRID <= 2 || (!defined(MHD) && NSUBGRID == 4))
#define GPU_LAUNCH_TARGET 100
#elif defined(GPU_SM) && GPU_SM == 90 && defined(MHD)
#define GPU_LAUNCH_TARGET 90
#elif defined(GPU_SM) && GPU_SM == 80 && defined(MHD)
#define GPU_LAUNCH_TARGET 80
#else
#define GPU_LAUNCH_TARGET 0
#endif

! Pure hydro without passive scalars traces pairs of cells at slope_type=3 in its own kernel, where
! it is faster: NSUBGRID=4, and NSUBGRID=2 in single precision.
#if !defined(MHD) && !(NVAR>5+NENER) && defined(NSUBGRID) && (NSUBGRID == 4 || (NSUBGRID == 2 && defined(NPRE) && NPRE == 4))
#define HYDRO_PAIR
#endif

! Kernels: HYDRO (hydro_integrator_kernel), HYDRO_PAIR (hydro_integrator_pair_kernel, hydro at
! slope_type=3), UCT_FACE, UCT_GHOST_FACE, UCT_UPDATE and UCT_UPDATE_AMR (uct_*_kernel, uct- solvers),
! NIMHD_DT and NIMHD_STAGE (nimhd_*_kernel, nimhd_solver sts and explicit). The non-ideal MHD values
! are the same on every GPU; NIMHD_DT_MINB stays 0, because nimhd_dt_kernel is defined before
! kernels without a bound (see below), and NIMHD_STAGE_MINB must be positive.
! B200 (sm_100), MHD
#if GPU_LAUNCH_TARGET == 100 && defined(MHD) && NPRE == 4 && NSUBGRID <= 1 && defined(GPU_FASTMATH)
#define HYDRO_THREADS 64
#define HYDRO_MINB 14
#define UCT_FACE_THREADS 32
#define UCT_FACE_MINB 20
#define UCT_GHOST_FACE_THREADS 32
#define UCT_GHOST_FACE_MINB 6
#define UCT_UPDATE_THREADS 32
#define UCT_UPDATE_MINB 20
#define UCT_UPDATE_AMR_THREADS 32
#define UCT_UPDATE_AMR_MINB 6
#define NIMHD_DT_THREADS 64
#define NIMHD_DT_MINB 0
#define NIMHD_STAGE_THREADS 64
#define NIMHD_STAGE_MINB 5
#elif GPU_LAUNCH_TARGET == 100 && defined(MHD) && NPRE == 4 && NSUBGRID <= 1 && !defined(GPU_FASTMATH)
#define HYDRO_THREADS 64
#define HYDRO_MINB 12
#define UCT_FACE_THREADS 32
#define UCT_FACE_MINB 20
#define UCT_GHOST_FACE_THREADS 32
#define UCT_GHOST_FACE_MINB 6
#define UCT_UPDATE_THREADS 32
#define UCT_UPDATE_MINB 20
#define UCT_UPDATE_AMR_THREADS 32
#define UCT_UPDATE_AMR_MINB 6
#define NIMHD_DT_THREADS 64
#define NIMHD_DT_MINB 0
#define NIMHD_STAGE_THREADS 64
#define NIMHD_STAGE_MINB 5
#elif GPU_LAUNCH_TARGET == 100 && defined(MHD) && NPRE == 4 && NSUBGRID == 2 && defined(GPU_FASTMATH)
#define HYDRO_THREADS 256
#define HYDRO_MINB 3
#define UCT_FACE_THREADS 128
#define UCT_FACE_MINB 7
#define UCT_GHOST_FACE_THREADS 256
#define UCT_GHOST_FACE_MINB 3
#define UCT_UPDATE_THREADS 96
#define UCT_UPDATE_MINB 12
#define UCT_UPDATE_AMR_THREADS 128
#define UCT_UPDATE_AMR_MINB 7
#define NIMHD_DT_THREADS 256
#define NIMHD_DT_MINB 0
#define NIMHD_STAGE_THREADS 256
#define NIMHD_STAGE_MINB 5
#elif GPU_LAUNCH_TARGET == 100 && defined(MHD) && NPRE == 4 && NSUBGRID == 2 && !defined(GPU_FASTMATH)
#define HYDRO_THREADS 256
#define HYDRO_MINB 3
#define UCT_FACE_THREADS 128
#define UCT_FACE_MINB 7
#define UCT_GHOST_FACE_THREADS 256
#define UCT_GHOST_FACE_MINB 3
#define UCT_UPDATE_THREADS 128
#define UCT_UPDATE_MINB 9
#define UCT_UPDATE_AMR_THREADS 128
#define UCT_UPDATE_AMR_MINB 7
#define NIMHD_DT_THREADS 256
#define NIMHD_DT_MINB 0
#define NIMHD_STAGE_THREADS 256
#define NIMHD_STAGE_MINB 5
#elif GPU_LAUNCH_TARGET == 100 && defined(MHD) && NPRE != 4 && NSUBGRID <= 1
#define HYDRO_THREADS 64
#define HYDRO_MINB 5
#define UCT_FACE_THREADS 32
#define UCT_FACE_MINB 12
#define UCT_GHOST_FACE_THREADS 32
#define UCT_GHOST_FACE_MINB 6
#define UCT_UPDATE_THREADS 32
#define UCT_UPDATE_MINB 20
#define UCT_UPDATE_AMR_THREADS 32
#define UCT_UPDATE_AMR_MINB 16
#define NIMHD_DT_THREADS 64
#define NIMHD_DT_MINB 0
#define NIMHD_STAGE_THREADS 64
#define NIMHD_STAGE_MINB 5
#elif GPU_LAUNCH_TARGET == 100 && defined(MHD) && NPRE != 4 && NSUBGRID == 2
#define HYDRO_THREADS 320
#define HYDRO_MINB 1
#define UCT_FACE_THREADS 128
#define UCT_FACE_MINB 0
#define UCT_GHOST_FACE_THREADS 256
#define UCT_GHOST_FACE_MINB 2
#define UCT_UPDATE_THREADS 128
#define UCT_UPDATE_MINB 6
#define UCT_UPDATE_AMR_THREADS 128
#define UCT_UPDATE_AMR_MINB 6
#define NIMHD_DT_THREADS 256
#define NIMHD_DT_MINB 0
#define NIMHD_STAGE_THREADS 256
#define NIMHD_STAGE_MINB 5
! B200 (sm_100), hydro
#elif GPU_LAUNCH_TARGET == 100 && !defined(MHD) && NPRE == 4 && NSUBGRID <= 1 && defined(GPU_FASTMATH)
#define HYDRO_THREADS 32
#define HYDRO_MINB 32
#elif GPU_LAUNCH_TARGET == 100 && !defined(MHD) && NPRE == 4 && NSUBGRID <= 1 && !defined(GPU_FASTMATH)
#define HYDRO_THREADS 32
#define HYDRO_MINB 0
#elif GPU_LAUNCH_TARGET == 100 && !defined(MHD) && NPRE == 4 && NSUBGRID == 2
#define HYDRO_THREADS 96
#define HYDRO_MINB 0
#define HYDRO_PAIR_THREADS 96
#define HYDRO_PAIR_MINB 8
#elif GPU_LAUNCH_TARGET == 100 && !defined(MHD) && NPRE == 4 && NSUBGRID == 4 && defined(GPU_FASTMATH)
#define HYDRO_THREADS 608
#define HYDRO_MINB 2
#define HYDRO_PAIR_THREADS 448
#define HYDRO_PAIR_MINB 2
#elif GPU_LAUNCH_TARGET == 100 && !defined(MHD) && NPRE == 4 && NSUBGRID == 4 && !defined(GPU_FASTMATH)
#define HYDRO_THREADS 512
#define HYDRO_MINB 0
#define HYDRO_PAIR_THREADS 448
#define HYDRO_PAIR_MINB 2
#elif GPU_LAUNCH_TARGET == 100 && !defined(MHD) && NPRE != 4 && NSUBGRID <= 1
#define HYDRO_THREADS 32
#define HYDRO_MINB 16
#elif GPU_LAUNCH_TARGET == 100 && !defined(MHD) && NPRE != 4 && NSUBGRID == 2
#define HYDRO_THREADS 256
#define HYDRO_MINB 4
#elif GPU_LAUNCH_TARGET == 100 && !defined(MHD) && NPRE != 4 && NSUBGRID == 4
#define HYDRO_THREADS 640
#define HYDRO_MINB 1
#define HYDRO_PAIR_THREADS 512
#define HYDRO_PAIR_MINB 1
! H200 (sm_90), MHD
#elif GPU_LAUNCH_TARGET == 90 && defined(MHD) && NPRE == 4 && NSUBGRID <= 1
#define HYDRO_THREADS 64
#define HYDRO_MINB 12
#define UCT_FACE_THREADS 64
#define UCT_FACE_MINB 14
#define UCT_GHOST_FACE_THREADS 64
#define UCT_GHOST_FACE_MINB 4
#define UCT_UPDATE_THREADS 64
#define UCT_UPDATE_MINB 20
#define UCT_UPDATE_AMR_THREADS 64
#define UCT_UPDATE_AMR_MINB 4
#define NIMHD_DT_THREADS 64
#define NIMHD_DT_MINB 0
#define NIMHD_STAGE_THREADS 64
#define NIMHD_STAGE_MINB 5
#elif GPU_LAUNCH_TARGET == 90 && defined(MHD) && NPRE == 4 && NSUBGRID == 2
#define HYDRO_THREADS 256
#define HYDRO_MINB 3
#define UCT_FACE_THREADS 128
#define UCT_FACE_MINB 7
#define UCT_GHOST_FACE_THREADS 256
#define UCT_GHOST_FACE_MINB 3
#define UCT_UPDATE_THREADS 128
#define UCT_UPDATE_MINB 9
#define UCT_UPDATE_AMR_THREADS 128
#define UCT_UPDATE_AMR_MINB 7
#define NIMHD_DT_THREADS 256
#define NIMHD_DT_MINB 0
#define NIMHD_STAGE_THREADS 256
#define NIMHD_STAGE_MINB 5
#elif GPU_LAUNCH_TARGET == 90 && defined(MHD) && NPRE != 4 && NSUBGRID <= 1
#define HYDRO_THREADS 64
#define HYDRO_MINB 5
#define UCT_FACE_THREADS 64
#define UCT_FACE_MINB 8
#define UCT_GHOST_FACE_THREADS 64
#define UCT_GHOST_FACE_MINB 4
#define UCT_UPDATE_THREADS 64
#define UCT_UPDATE_MINB 16
#define UCT_UPDATE_AMR_THREADS 64
#define UCT_UPDATE_AMR_MINB 12
#define NIMHD_DT_THREADS 64
#define NIMHD_DT_MINB 0
#define NIMHD_STAGE_THREADS 64
#define NIMHD_STAGE_MINB 5
#elif GPU_LAUNCH_TARGET == 90 && defined(MHD) && NPRE != 4 && NSUBGRID == 2
#define HYDRO_THREADS 256
#define HYDRO_MINB 0
#define UCT_FACE_THREADS 128
#define UCT_FACE_MINB 0
#define UCT_GHOST_FACE_THREADS 256
#define UCT_GHOST_FACE_MINB 2
#define UCT_UPDATE_THREADS 128
#define UCT_UPDATE_MINB 6
#define UCT_UPDATE_AMR_THREADS 128
#define UCT_UPDATE_AMR_MINB 6
#define NIMHD_DT_THREADS 256
#define NIMHD_DT_MINB 0
#define NIMHD_STAGE_THREADS 256
#define NIMHD_STAGE_MINB 5
! A100 (sm_80), MHD
#elif GPU_LAUNCH_TARGET == 80 && defined(MHD) && NPRE == 4 && NSUBGRID <= 1
#define HYDRO_THREADS 64
#define HYDRO_MINB 0
#define UCT_FACE_THREADS 64
#define UCT_FACE_MINB 12
#define UCT_GHOST_FACE_THREADS 64
#define UCT_GHOST_FACE_MINB 4
#define UCT_UPDATE_THREADS 64
#define UCT_UPDATE_MINB 24
#define UCT_UPDATE_AMR_THREADS 64
#define UCT_UPDATE_AMR_MINB 4
#define NIMHD_DT_THREADS 64
#define NIMHD_DT_MINB 0
#define NIMHD_STAGE_THREADS 64
#define NIMHD_STAGE_MINB 5
#elif GPU_LAUNCH_TARGET == 80 && defined(MHD) && NPRE == 4 && NSUBGRID == 2
#define HYDRO_THREADS 256
#define HYDRO_MINB 0
#define UCT_FACE_THREADS 192
#define UCT_FACE_MINB 4
#define UCT_GHOST_FACE_THREADS 256
#define UCT_GHOST_FACE_MINB 3
#define UCT_UPDATE_THREADS 192
#define UCT_UPDATE_MINB 8
#define UCT_UPDATE_AMR_THREADS 192
#define UCT_UPDATE_AMR_MINB 7
#define NIMHD_DT_THREADS 256
#define NIMHD_DT_MINB 0
#define NIMHD_STAGE_THREADS 256
#define NIMHD_STAGE_MINB 5
#elif GPU_LAUNCH_TARGET == 80 && defined(MHD) && NPRE != 4 && NSUBGRID <= 1
#define HYDRO_THREADS 64
#define HYDRO_MINB 6
#define UCT_FACE_THREADS 64
#define UCT_FACE_MINB 8
#define UCT_GHOST_FACE_THREADS 64
#define UCT_GHOST_FACE_MINB 4
#define UCT_UPDATE_THREADS 64
#define UCT_UPDATE_MINB 16
#define UCT_UPDATE_AMR_THREADS 64
#define UCT_UPDATE_AMR_MINB 12
#define NIMHD_DT_THREADS 64
#define NIMHD_DT_MINB 0
#define NIMHD_STAGE_THREADS 64
#define NIMHD_STAGE_MINB 5
#elif GPU_LAUNCH_TARGET == 80 && defined(MHD) && NPRE != 4 && NSUBGRID == 2
#define HYDRO_THREADS 256
#define HYDRO_MINB 1
#define UCT_FACE_THREADS 192
#define UCT_FACE_MINB 0
#define UCT_GHOST_FACE_THREADS 256
#define UCT_GHOST_FACE_MINB 2
#define UCT_UPDATE_THREADS 192
#define UCT_UPDATE_MINB 4
#define UCT_UPDATE_AMR_THREADS 192
#define UCT_UPDATE_AMR_MINB 4
#define NIMHD_DT_THREADS 256
#define NIMHD_DT_MINB 0
#define NIMHD_STAGE_THREADS 256
#define NIMHD_STAGE_MINB 5
! Generic, MHD
#elif GPU_LAUNCH_TARGET == 0 && defined(MHD) && NPRE == 4 && NSUBGRID <= 1
#define HYDRO_THREADS 64
#define HYDRO_MINB 0
#define UCT_FACE_THREADS 64
#define UCT_FACE_MINB 14
#define UCT_GHOST_FACE_THREADS 64
#define UCT_GHOST_FACE_MINB 4
#define UCT_UPDATE_THREADS 64
#define UCT_UPDATE_MINB 20
#define UCT_UPDATE_AMR_THREADS 64
#define UCT_UPDATE_AMR_MINB 4
#define NIMHD_DT_THREADS 64
#define NIMHD_DT_MINB 0
#define NIMHD_STAGE_THREADS 64
#define NIMHD_STAGE_MINB 5
#elif GPU_LAUNCH_TARGET == 0 && defined(MHD) && NPRE == 4 && NSUBGRID == 2
#define HYDRO_THREADS 256
#define HYDRO_MINB 0
#define UCT_FACE_THREADS 128
#define UCT_FACE_MINB 7
#define UCT_GHOST_FACE_THREADS 256
#define UCT_GHOST_FACE_MINB 3
#define UCT_UPDATE_THREADS 128
#define UCT_UPDATE_MINB 9
#define UCT_UPDATE_AMR_THREADS 128
#define UCT_UPDATE_AMR_MINB 7
#define NIMHD_DT_THREADS 256
#define NIMHD_DT_MINB 0
#define NIMHD_STAGE_THREADS 256
#define NIMHD_STAGE_MINB 5
#elif GPU_LAUNCH_TARGET == 0 && defined(MHD) && NPRE != 4 && NSUBGRID <= 1
#define HYDRO_THREADS 64
#define HYDRO_MINB 0
#define UCT_FACE_THREADS 64
#define UCT_FACE_MINB 8
#define UCT_GHOST_FACE_THREADS 64
#define UCT_GHOST_FACE_MINB 4
#define UCT_UPDATE_THREADS 64
#define UCT_UPDATE_MINB 16
#define UCT_UPDATE_AMR_THREADS 64
#define UCT_UPDATE_AMR_MINB 12
#define NIMHD_DT_THREADS 64
#define NIMHD_DT_MINB 0
#define NIMHD_STAGE_THREADS 64
#define NIMHD_STAGE_MINB 5
#elif GPU_LAUNCH_TARGET == 0 && defined(MHD) && NPRE != 4 && NSUBGRID == 2
#define HYDRO_THREADS 256
#define HYDRO_MINB 0
#define UCT_FACE_THREADS 128
#define UCT_FACE_MINB 0
#define UCT_GHOST_FACE_THREADS 256
#define UCT_GHOST_FACE_MINB 2
#define UCT_UPDATE_THREADS 128
#define UCT_UPDATE_MINB 6
#define UCT_UPDATE_AMR_THREADS 128
#define UCT_UPDATE_AMR_MINB 6
#define NIMHD_DT_THREADS 256
#define NIMHD_DT_MINB 0
#define NIMHD_STAGE_THREADS 256
#define NIMHD_STAGE_MINB 5
! Generic, hydro
#elif GPU_LAUNCH_TARGET == 0 && !defined(MHD) && NPRE == 4 && NSUBGRID <= 1
#define HYDRO_THREADS 64
#define HYDRO_MINB 0
#elif GPU_LAUNCH_TARGET == 0 && !defined(MHD) && NPRE == 4 && NSUBGRID == 2
#define HYDRO_THREADS 128
#define HYDRO_MINB 0
#define HYDRO_PAIR_THREADS 128
#define HYDRO_PAIR_MINB 8
#elif GPU_LAUNCH_TARGET == 0 && !defined(MHD) && NPRE == 4 && NSUBGRID == 4
#define HYDRO_THREADS 512
#define HYDRO_MINB 0
#define HYDRO_PAIR_THREADS 512
#define HYDRO_PAIR_MINB 2
#elif GPU_LAUNCH_TARGET == 0 && !defined(MHD) && NPRE != 4 && NSUBGRID <= 1
#define HYDRO_THREADS 64
#define HYDRO_MINB 16
#elif GPU_LAUNCH_TARGET == 0 && !defined(MHD) && NPRE != 4 && NSUBGRID == 2
#define HYDRO_THREADS 256
#define HYDRO_MINB 4
#elif GPU_LAUNCH_TARGET == 0 && !defined(MHD) && NPRE != 4 && NSUBGRID == 4
#define HYDRO_THREADS 640
#define HYDRO_MINB 1
#define HYDRO_PAIR_THREADS 640
#define HYDRO_PAIR_MINB 1
#else
#error "gpu_launch.h: CUDA builds support NSUBGRID=1, 2 or 4, and 1 or 2 with MHD"
#endif

! NVHPC applies a kernel's launch_bounds to every kernel that follows it in the module, so a kernel
! with _MINB 0 keeps the bound of the bounded kernel before it, if any: the double precision UCT face
! kernel at NSUBGRID=2 that of the fused integrator, which must admit its launch size.
#if HYDRO_MINB > 0
#define HYDRO_LB launch_bounds(HYDRO_THREADS,HYDRO_MINB)
#else
#define HYDRO_LB
#endif
#ifdef HYDRO_PAIR
#if HYDRO_PAIR_MINB > 0
#define HYDRO_PAIR_LB launch_bounds(HYDRO_PAIR_THREADS,HYDRO_PAIR_MINB)
#else
#define HYDRO_PAIR_LB
#endif
#endif
#ifdef MHD
#if UCT_FACE_MINB > 0
#define UCT_FACE_LB launch_bounds(UCT_FACE_THREADS,UCT_FACE_MINB)
#else
#define UCT_FACE_LB
#endif
#if UCT_GHOST_FACE_MINB > 0
#define UCT_GHOST_FACE_LB launch_bounds(UCT_GHOST_FACE_THREADS,UCT_GHOST_FACE_MINB)
#else
#define UCT_GHOST_FACE_LB
#endif
#if UCT_UPDATE_MINB > 0
#define UCT_UPDATE_LB launch_bounds(UCT_UPDATE_THREADS,UCT_UPDATE_MINB)
#else
#define UCT_UPDATE_LB
#endif
#if UCT_UPDATE_AMR_MINB > 0
#define UCT_UPDATE_AMR_LB launch_bounds(UCT_UPDATE_AMR_THREADS,UCT_UPDATE_AMR_MINB)
#else
#define UCT_UPDATE_AMR_LB
#endif
#if NIMHD_DT_MINB > 0
#define NIMHD_DT_LB launch_bounds(NIMHD_DT_THREADS,NIMHD_DT_MINB)
#else
#define NIMHD_DT_LB
#endif
! The stage kernel keeps the bound of its NSUBGRID=2 launch, 256 threads, at both tiles.
#define NIMHD_STAGE_LB launch_bounds(256,NIMHD_STAGE_MINB)
#endif
