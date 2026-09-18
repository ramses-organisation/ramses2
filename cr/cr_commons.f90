module cr_commons
  use amr_parameters
  use oct_commons

  type cr_kernel_t
     integer::iu1,iu2,ju1,ju2,ku1,ku2
     integer::if1,if2,jf1,jf2,kf1,kf2
     integer::io1,io2,jo1,jo2,ko1,ko2
     logical,dimension(:,:,:),allocatable::inkernel
     logical,dimension(:,:,:),allocatable::okloc
     integer,dimension(:,:,:),allocatable::cellloc
     integer,dimension(:,:,:),allocatable::childloc
     integer,dimension(:,:,:),allocatable::gridloc
     integer,dimension(:,:,:,:),allocatable::nborloc
     real(kind=8),dimension(:,:,:,:),allocatable::uloc
#ifdef MHD
     real(kind=8),dimension(:,:,:,:),allocatable::bloc
#endif
#ifdef DO_CR
     real(kind=8),dimension(:,:,:,:),allocatable::cruloc
     real(kind=8),dimension(:,:,:,:,:),allocatable::crflux
     real(kind=8),dimension(:,:,:,:,:),allocatable::cFlx
     real(kind=8),dimension(:,:,:,:),allocatable::lmax
#endif
   contains
     procedure :: init => init_cr_kernel
     procedure :: size => size_cr_kernel
  end type cr_kernel_t

  type cr_workspace_t
     type(cr_kernel_t)::kernel_1,kernel_2,kernel_4,kernel_8,kernel_16,kernel_32
  end type cr_workspace_t

contains

  subroutine init_cr_kernel(h,nn)
    use amr_parameters, only: ndim
    use hydro_parameters, only: nvar
    use cr_parameters, only: ncruvar, ncrgrp

    integer::nn
    class(cr_kernel_t)::h

    h%iu1=-1; h%iu2=nn+2
    h%ju1=(1-ndim/2)-1*(ndim/2); h%ju2=(1-ndim/2)+(nn+2)*(ndim/2)
    h%ku1=(1-ndim/3)-1*(ndim/3); h%ku2=(1-ndim/3)+(nn+2)*(ndim/3)
    h%if1=1; h%if2=nn+1
    h%jf1=1; h%jf2=(1-ndim/2)+(nn+1)*(ndim/2)
    h%kf1=1; h%kf2=(1-ndim/3)+(nn+1)*(ndim/3)
    h%io1=0; h%io2=nn/2+1
    h%jo1=(1-ndim/2); h%jo2=(1-ndim/2)+(nn/2+1)*(ndim/2)
    h%ko1=(1-ndim/3); h%ko2=(1-ndim/3)+(nn/2+1)*(ndim/3)

    allocate(h%uloc (h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2,1:nvar))
    allocate(h%okloc(h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2))
    allocate(h%childloc(h%io1:h%io2,h%jo1:h%jo2,h%ko1:h%ko2))
    allocate(h%gridloc (h%io1:h%io2,h%jo1:h%jo2,h%ko1:h%ko2))
    allocate(h%cellloc (h%io1:h%io2,h%jo1:h%jo2,h%ko1:h%ko2))
    allocate(h%nborloc (h%io1:h%io2,h%jo1:h%jo2,h%ko1:h%ko2,1:twondim))
#ifdef MHD
    allocate(h%bloc(h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2,1:6))
#endif
#ifdef DO_CR
    allocate(h%cruloc (h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2,1:ncruvar))
    allocate(h%crflux(h%if1:h%if2,h%jf1:h%jf2,h%kf1:h%kf2,1:ncruvar+ncrgrp,1:ndim))
    allocate(h%cFlx(h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2,1:ndim+1,1:ndim))
    allocate(h%lmax(h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2,1:ndim))
#endif
    allocate(h%inkernel(h%io1:h%io2,h%jo1:h%jo2,h%ko1:h%ko2))
  end subroutine init_cr_kernel

  function size_cr_kernel(h)
    integer::size_cr_kernel
    class(cr_kernel_t)::h

    integer::nint

    nint=0
    nint=nint+size(transfer(h%okloc,(/1/)))

    nint=nint+size(transfer(h%inkernel,(/1/)))
    nint=nint+size(transfer(h%childloc,(/1/)))
    nint=nint+size(transfer(h%gridloc ,(/1/)))
    nint=nint+size(transfer(h%cellloc ,(/1/)))
    nint=nint+size(transfer(h%nborloc ,(/1/)))

    nint=nint+size(transfer(h%uloc,(/1/)))
#ifdef MHD
    nint=nint+size(transfer(h%bloc,(/1/)))
#endif
#ifdef DO_CR
    nint=nint+size(transfer(h%cruloc,(/1/)))
    nint=nint+size(transfer(h%crflux,(/1/)))
    nint=nint+size(transfer(h%cFlx,(/1/)))
    nint=nint+size(transfer(h%lmax,(/1/)))
#endif
    size_cr_kernel = nint

  end function size_cr_kernel

end module cr_commons
