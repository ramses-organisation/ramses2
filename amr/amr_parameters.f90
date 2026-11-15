module amr_parameters

  ! Define real types
  integer,parameter::sp=kind(1.0E0)
#ifndef NPRE
  integer,parameter::dp=kind(1.0D0) ! default
#else
#if NPRE==4
  integer,parameter::dp=kind(1.0E0) ! real*4
#else
  integer,parameter::dp=kind(1.0D0) ! real*8
#endif
#endif

  integer,parameter::MAXOUT=1000
  integer,parameter::MAXLEVEL=65
  integer,parameter::MAXREGION=100
  integer,parameter::MAXBOUND=100

  ! Define integer types (for particle IDs only)
  integer,parameter::i4b=4
#ifndef LONGINT
  integer,parameter::i8b=4  ! default particle IDs are short int
#else
  integer,parameter::i8b=8  ! longint particle IDs are long int
#endif
  ! Length of a file path, and room for the name a reader appends to a directory
  integer,parameter::flen=512
  integer,parameter::flen_tail=64

  ! Number of dimensions
#ifndef NDIM
  integer,parameter::ndim=1
#else
  integer,parameter::ndim=NDIM
#endif
  integer,parameter::twotondim=2**ndim
  integer,parameter::threetondim=3**ndim
  integer,parameter::fourtondim=4**ndim
  integer,parameter::fivetondim=5**ndim

  ! Number of octs per dimension in a block (CPU and GPU kernels)
#ifndef NSUBGRID
  integer,parameter::nsubgrid=1
#else
#if NSUBGRID!=1 && NSUBGRID!=2 && NSUBGRID!=4 && NSUBGRID!=8 && NSUBGRID!=16 && NSUBGRID!=32
#error "NSUBGRID must be 1, 2, 4, 8, 16 or 32"
#endif
  integer,parameter::nsubgrid=NSUBGRID
#endif
  integer,parameter::twondim=2*ndim

  ! Number of 64-bit integers needed to store one Hilbert key
#ifndef NHILBERT
  integer, parameter :: nhilbert = 1
#else
  integer, parameter :: nhilbert = NHILBERT
#endif

  ! Vectorization parameter
#ifndef NVECTOR
  integer,parameter::nvector=32  ! Size of vector sweeps
#else
  integer,parameter::nvector=NVECTOR
#endif

  ! Size in bytes of the buffers used for grid output
  integer,parameter::nbuf_io=4194304

  ! Expansion factor look-up table size
  integer,parameter::n_frw=1000

  ! Number of bins for halo mass profiles
  integer, parameter::nbin=100

  ! Executable identification
  CHARACTER(LEN=300)::builddate,buildcommand,patchdir
  CHARACTER(LEN=300)::gitrepo,gitbranch,githash

  ! Save namelist filename
  CHARACTER(LEN=300)::namelist_file

end module amr_parameters
