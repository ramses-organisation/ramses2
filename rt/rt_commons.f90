module rt_commons
  use amr_parameters
  use oct_commons

  type rt_kernel_t
     integer::iu1,iu2,ju1,ju2,ku1,ku2
     integer::if1,if2,jf1,jf2,kf1,kf2
     integer::io1,io2,jo1,jo2,ko1,ko2
     logical,dimension(:,:,:),allocatable::inkernel
     logical,dimension(:,:,:),allocatable::okloc
     integer,dimension(:,:,:),allocatable::cellloc
     integer,dimension(:,:,:),allocatable::childloc
     integer,dimension(:,:,:),allocatable::gridloc
     integer,dimension(:,:,:,:),allocatable::nborloc
#ifdef DO_RT
     real(kind=8),dimension(:,:,:,:),allocatable::rtuloc
     real(kind=8),dimension(:,:,:,:,:),allocatable::rtflux
     real(kind=8),dimension(:,:,:,:,:),allocatable::cFlx
#endif
   contains
     procedure :: init => init_rt_kernel
     procedure :: size => size_rt_kernel
  end type rt_kernel_t

  type rt_workspace_t
     type(rt_kernel_t)::kernel
  end type rt_workspace_t

contains

  subroutine init_rt_kernel(h,nn)
    use amr_parameters, only: ndim
    use rt_parameters, only: nrtvar

    integer::nn
    class(rt_kernel_t)::h

    h%iu1=-1; h%iu2=nn+2
    h%ju1=(1-ndim/2)-1*(ndim/2); h%ju2=(1-ndim/2)+(nn+2)*(ndim/2)
    h%ku1=(1-ndim/3)-1*(ndim/3); h%ku2=(1-ndim/3)+(nn+2)*(ndim/3)
    h%if1=1; h%if2=nn+1
    h%jf1=1; h%jf2=(1-ndim/2)+(nn+1)*(ndim/2)
    h%kf1=1; h%kf2=(1-ndim/3)+(nn+1)*(ndim/3)
    h%io1=0; h%io2=nn/2+1
    h%jo1=(1-ndim/2); h%jo2=(1-ndim/2)+(nn/2+1)*(ndim/2)
    h%ko1=(1-ndim/3); h%ko2=(1-ndim/3)+(nn/2+1)*(ndim/3)

    allocate(h%okloc(h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2))

    allocate(h%inkernel(h%io1:h%io2,h%jo1:h%jo2,h%ko1:h%ko2))
    allocate(h%childloc(h%io1:h%io2,h%jo1:h%jo2,h%ko1:h%ko2))
    allocate(h%gridloc (h%io1:h%io2,h%jo1:h%jo2,h%ko1:h%ko2))
    allocate(h%cellloc (h%io1:h%io2,h%jo1:h%jo2,h%ko1:h%ko2))
    allocate(h%nborloc (h%io1:h%io2,h%jo1:h%jo2,h%ko1:h%ko2,1:twondim))

#ifdef DO_RT
    allocate(h%rtuloc (h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2,1:nrtvar))
    allocate(h%rtflux(h%if1:h%if2,h%jf1:h%jf2,h%kf1:h%kf2,1:nrtvar,1:ndim))
    allocate(h%cFlx(h%iu1:h%iu2,h%ju1:h%ju2,h%ku1:h%ku2,1:ndim+1,1:ndim))
#endif
  end subroutine init_rt_kernel

  function size_rt_kernel(h)
    integer::size_rt_kernel
    class(rt_kernel_t)::h

    integer::nint

    nint=0
    nint=nint+size(transfer(h%okloc,(/1/)))

    nint=nint+size(transfer(h%inkernel,(/1/)))
    nint=nint+size(transfer(h%childloc,(/1/)))
    nint=nint+size(transfer(h%gridloc ,(/1/)))
    nint=nint+size(transfer(h%cellloc ,(/1/)))
    nint=nint+size(transfer(h%nborloc ,(/1/)))

#ifdef DO_RT
    nint=nint+size(transfer(h%rtuloc,(/1/)))
    nint=nint+size(transfer(h%rtflux,(/1/)))
    nint=nint+size(transfer(h%cFlx,(/1/)))
#endif
    size_rt_kernel = nint

  end function size_rt_kernel

end module rt_commons
