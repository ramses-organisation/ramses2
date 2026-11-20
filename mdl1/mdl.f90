module mdl_module

  use mdl_parameters
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_FUNPTR, C_PTR
  use amr_commons, only: mesh_t
  use clfind_commons, only: clump_t

  type ::comm_buff
     integer(kind=4),dimension(:),allocatable::array
  end type comm_buff

  type :: mdl_t

     integer,PRIVATE::myid
     integer,PRIVATE::ncpu
#ifdef _CUDA
     integer,PRIVATE::mydev
     integer,PRIVATE::ngpu
#endif
     integer::MDL_INPUT_MAXSIZE=1
     integer,dimension(:),allocatable::mpi_input_buffer

     ! Cache opened or closed
     logical::cache_opened=.false.
     logical::cache_opened_clump=.false.

     ! Communication-related objects
     integer::mail_counter
     integer::request_id,flush_id
     integer::request_id_clump,flush_id_clump
     integer,dimension(:),allocatable::reply_id
     integer,dimension(:),allocatable::reply_id_clump

     ! Adopted combiner rule
     integer::combiner_rule

     ! Message sizes
     integer::size_msg_array
     integer::size_request_array
     integer::size_flush_array
     integer::size_fetch_array

     ! Message sizes clump
     integer::size_msg_array_clump
     integer::size_request_array_clump
     integer::size_flush_array_clump
     integer::size_fetch_array_clump

     ! Message buffers
     integer(kind=4),dimension(:),allocatable::recv_request_array
     integer(kind=4),dimension(:),allocatable::send_request_array
     integer(kind=4),dimension(:),allocatable::recv_fetch_array
     integer(kind=4),dimension(:),allocatable::recv_flush_array

     ! Message buffers clump
     integer(kind=4),dimension(:),allocatable::recv_request_array_clump
     integer(kind=4),dimension(:),allocatable::send_request_array_clump
     integer(kind=4),dimension(:),allocatable::recv_fetch_array_clump
     integer(kind=4),dimension(:),allocatable::recv_flush_array_clump

     ! Send buffers
     integer(kind=4)::nbuffer_fetch
     integer(kind=4)::nbuffer_flush
     integer(kind=4)::ibuffer_fetch
     integer(kind=4)::ibuffer_flush
     integer(kind=4),dimension(:),allocatable::cpu2buf_fetch
     integer(kind=4),dimension(:),allocatable::cpu2buf_flush
     type(comm_buff),dimension(:),allocatable::send_fetch
     type(comm_buff),dimension(:),allocatable::send_flush

     ! Send buffers clump
     integer(kind=4)::nbuffer_fetch_clump
     integer(kind=4)::nbuffer_flush_clump
     integer(kind=4)::ibuffer_fetch_clump
     integer(kind=4)::ibuffer_flush_clump
     integer(kind=4),dimension(:),allocatable::cpu2buf_fetch_clump
     integer(kind=4),dimension(:),allocatable::cpu2buf_flush_clump
     type(comm_buff),dimension(:),allocatable::send_fetch_clump
     type(comm_buff),dimension(:),allocatable::send_flush_clump

     ! Callback functions
     type(c_funptr),dimension(0:200)::callback
     type(c_ptr),dimension(0:200)::p1opaque
     integer(kind=4),dimension(0:200)::input_size, output_size

     ! Pointers to data structurs
     type(mesh_t),pointer::m
     type(clump_t),pointer::c

  end type mdl_t

  interface mdl_send_request
    module procedure mdl_send_request_array, mdl_send_request_scalar
  end interface mdl_send_request
  PRIVATE :: mdl_send_request_array, mdl_send_request_scalar

  interface mdl_get_reply
    module procedure mdl_get_reply_array, mdl_get_reply_scalar
  end interface mdl_get_reply
  PRIVATE :: mdl_get_reply_array, mdl_get_reply_scalar

contains
  !##############################################################
  !##############################################################
  !##############################################################
  !##############################################################
  subroutine mdl_abort(mdl)
    type(mdl_t)::mdl
#ifndef WITHOUTMPI
    include 'mpif.h'
    integer::info
    ! Abort with error code 1
    call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
    stop 1
#endif
  end subroutine mdl_abort
  !##############################################################
  !##############################################################
  !##############################################################
  !##############################################################
  subroutine mdl_add_service(mdl,sid,p1,service,input_size,output_size,name)
    USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_FUNPTR, C_LOC
    type(mdl_t)::mdl
    integer::sid
    type(*),target::p1
    type(c_funptr), intent(in), value :: service
    integer :: input_size, output_size ! NOTE: size in units of integer
    character(len=*), intent(in),optional :: name
    mdl%callback(sid) = service
    mdl%p1opaque(sid) = c_loc(p1)
    mdl%input_size(sid) = input_size
    mdl%output_size(sid) = output_size
    mdl%MDL_INPUT_MAXSIZE = MAX(mdl%MDL_INPUT_MAXSIZE,input_size)

  end subroutine mdl_add_service
  !##############################################################
  !##############################################################
  !##############################################################
  !##############################################################
  integer function mdl_send_request_array(mdl,mdl_function_id,target_cpu,input_size,output_size,input_array)
    implicit none
    type(mdl_t)::mdl
    integer,intent(in)::mdl_function_id
    integer,intent(in)::target_cpu,input_size,output_size
    integer,intent(in),dimension(1:input_size)::input_array

#ifndef WITHOUTMPI
    include 'mpif.h'
    integer::info
    integer::launch_id,launch_tag=101
    integer,dimension(MPI_STATUS_SIZE)::launch_status
    integer,dimension(1:32)::header=0

!      write(*,*)'ARRAY:',mdl_function_id,input_size
!      if (mdl%input_size(mdl_function_id) .lt. input_size) then
!        write(*,*) 'BROKEN PROMISE: function_id=',mdl_function_id,' INPUT SIZE:',input_size,'MAX:',mdl%input_size(mdl_function_id)
!      end if
!      if (mdl%output_size(mdl_function_id) .lt. output_size) then
!        write(*,*) 'BROKEN PROMISE: function_id=',mdl_function_id,' OUTPUT SIZE:',output_size,'MAX:',mdl%output_size(mdl_function_id)
!      end if
!      if (mod(input_size,4).ne.0) then
!        write(*,*) 'SUSPICIOUS SIZE: function_id=',mdl_function_id,'SIZE:',input_size
!      end if

    ! Assemble MPI message
    header(1)=mdl_function_id
    header(2)=input_size
    header(3)=output_size
    mdl%mpi_input_buffer(1:32)=header
    if(input_size>0)then
       mdl%mpi_input_buffer(33:32+input_size)=input_array
    endif

    ! Send input array to the target cpu
    call MPI_ISEND(mdl%mpi_input_buffer,input_size+32,MPI_INTEGER,target_cpu-1,launch_tag,MPI_COMM_WORLD,launch_id,info)

    ! Wait for ISEND completion to free memory in corresponding MPI buffer
    call MPI_WAIT(launch_id,launch_status,info)
#endif
    mdl_send_request_array = target_cpu
  end function mdl_send_request_array
  !##############################################################
  !##############################################################
  !##############################################################
  !##############################################################
  integer function mdl_send_request_scalar(mdl,mdl_function_id,target_cpu,input_size,output_size,input)
    USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_LOC, C_F_POINTER
    implicit none
    type(mdl_t)::mdl
    integer,intent(in)::mdl_function_id
    integer,intent(in)::target_cpu
    integer,intent(in),optional::input_size,output_size
    type(*),intent(in),optional,target::input

    integer::isize
#ifndef WITHOUTMPI
    include 'mpif.h'
    integer::info
    integer::launch_id,launch_tag=101
    integer,dimension(MPI_STATUS_SIZE)::launch_status
    integer,dimension(1:32)::header=0
    integer,dimension(:),pointer::dummy

    if (present(input_size)) then
       isize = input_size
       if (mdl%input_size(mdl_function_id) .lt. input_size) then
!          write(*,*) 'BROKEN PROMISE: function_id=',mdl_function_id,'INPUT SIZE:',input_size,'MAX:',mdl%input_size(mdl_function_id)
       end if
!        if (mod(input_size,4).ne.0) then
!          write(*,*) 'SUSPICIOUS SIZE: function_id=',mdl_function_id,'SIZE:',input_size
!        end if
    else
       isize = 0
    end if
    if (present(output_size)) then
       if (mdl%output_size(mdl_function_id) .lt. output_size) then
!          write(*,*) 'BROKEN PROMISE: function_id=',mdl_function_id,' OUTPUT SIZE:',output_size,'MAX:',mdl%output_size(mdl_function_id)
       end if
    end if
!      write(*,*)'SCALAR:',mdl_function_id,isize

    ! Assemble MPI message
    header(1)=mdl_function_id
    header(2)=isize
    if (present(output_size)) then
       header(3)=output_size
    else
       header(3)=0
    end if
    mdl%mpi_input_buffer(1:32)=header
    if (isize .gt. 0) then
       call c_f_pointer(c_loc(input),dummy,[input_size])
       mdl%mpi_input_buffer(33:32+input_size) = transfer(dummy,mdl%mpi_input_buffer(33:32+input_size))
    endif

    ! Send input array to the target cpu
    call MPI_ISEND(mdl%mpi_input_buffer,isize+32,MPI_INTEGER,target_cpu-1,launch_tag,MPI_COMM_WORLD,launch_id,info)

    ! Wait for ISEND completion to free memory in corresponding MPI buffer
    call MPI_WAIT(launch_id,launch_status,info)
#endif
    mdl_send_request_scalar = target_cpu
  end function mdl_send_request_scalar
  !##############################################################
  !##############################################################
  !##############################################################
  !##############################################################
  subroutine mdl_get_reply_array(mdl,target_cpu,output_size,output_array)
    implicit none
    type(mdl_t)::mdl
    integer::target_cpu
    integer::output_size
    integer,dimension(1:output_size)::output_array

#ifndef WITHOUTMPI
    include 'mpif.h'
    integer::info
    integer::output_tag=203,output_id
    integer,dimension(1)::dummy
    integer,dimension(MPI_STATUS_SIZE)::output_status

    ! Post a RECV for the output back from target_cpu
    if(output_size>0)then
       call MPI_IRECV(output_array,output_size,MPI_INTEGER,target_cpu-1,output_tag,MPI_COMM_WORLD,output_id,info)
    else
       call MPI_IRECV(dummy,1,MPI_INTEGER,target_cpu-1,output_tag,MPI_COMM_WORLD,output_id,info)
    endif

    ! Wait for ISEND completion to free memory in corresponding MPI buffer
    call MPI_WAIT(output_id,output_status,info)
#endif
  end subroutine mdl_get_reply_array
  !##############################################################
  !##############################################################
  !##############################################################
  !##############################################################
  subroutine mdl_get_reply_scalar(mdl,target_cpu,output_length,output)
    USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_LOC, C_F_POINTER
    implicit none
    type(mdl_t)::mdl
    integer,value,intent(in)::target_cpu
    integer,optional::output_length
    type(*),optional,target::output
    integer,dimension(:),pointer::buffer

#ifndef WITHOUTMPI
    include 'mpif.h'
    integer::info
    integer::output_tag=203,output_id
    integer,dimension(1)::dummy
    integer,dimension(MPI_STATUS_SIZE)::output_status

    ! Post a RECV for the output back from target_cpu
    if(present(output)) then
       call c_f_pointer(c_loc(output),buffer,[mdl%MDL_INPUT_MAXSIZE/4])
       call MPI_IRECV(buffer,output_length,MPI_INTEGER,target_cpu-1,output_tag,MPI_COMM_WORLD,output_id,info)
    else
       call MPI_IRECV(dummy,1,MPI_INTEGER,target_cpu-1,output_tag,MPI_COMM_WORLD,output_id,info)
    endif

    ! Wait for ISEND completion to free memory in corresponding MPI buffer
    call MPI_WAIT(output_id,output_status,info)
#endif
  end subroutine mdl_get_reply_scalar
  !##############################################################
  !##############################################################
  !##############################################################
  !##############################################################
  subroutine mdl_initialize(mdl)
#ifdef _CUDA
    use cudafor
#endif
#ifdef _METAL
    use metal_interface, only: mtl_init
#endif
    type(mdl_t)::mdl
#ifndef WITHOUTMPI
    include 'mpif.h'
    integer::info
#endif
#ifdef _CUDA
    integer::status, local_rank
    type(cudaDeviceProp) :: prop
#ifndef WITHOUTMPI
    integer::node_rank, node_size, node_comm
    character(len=16)::value

    ! Prefer the MPI launcher's local rank when it runs inside a Slurm task.
    local_rank=-1
    call get_environment_variable('OMPI_COMM_WORLD_LOCAL_RANK',value)
    if(value==' ')call get_environment_variable('MPI_LOCALRANKID',value)
    if(value==' ')call get_environment_variable('SLURM_LOCALID',value)
    if(value/=' ')then
       read(value,*,iostat=info)local_rank
       if(info/=0)local_rank=-1
    endif
#else
    local_rank=0
#endif
    ! CUDA-aware MPI can select a transport during MPI_INIT.
    status=cudaGetDeviceCount(mdl%ngpu)
    mdl%mydev=-1
    if(status==cudaSuccess)then
       if(mdl%ngpu==0)error stop 'No CUDA device visible'
       if(mdl%ngpu==1)local_rank=0
       if(local_rank>=0)then
          mdl%mydev=mod(local_rank,mdl%ngpu)
          status=cudaSetDevice(mdl%mydev)
       endif
    endif
#endif
#ifndef WITHOUTMPI
    call MPI_INIT(info)
    call MPI_COMM_RANK(MPI_COMM_WORLD,mdl%myid,info)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,mdl%ncpu,info)
    mdl%myid=mdl%myid+1 ! Careful with this...
    if(mdl_self(mdl)==1)then
       write(*,'(" Launching MPI with ntask = ",I0)')mdl%ncpu
    endif
#else
    write(*,'(" Serial execution (no MPI).")')
    mdl%ncpu=1
    mdl%myid=1
#endif
#ifdef _CUDA
    call check_cuda(status,'initial CUDA device selection')
#ifndef WITHOUTMPI
    node_rank=0
    if(mdl%ncpu>1)then
       call MPI_COMM_SPLIT_TYPE(MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,0,MPI_INFO_NULL,node_comm,info)
       call MPI_COMM_RANK(node_comm,node_rank,info)
       call MPI_COMM_SIZE(node_comm,node_size,info)
       call MPI_COMM_FREE(node_comm,info)
       if(node_rank==0.and.node_size>mdl%ngpu)write(*,'(" Rank ",I0,": ",I0," ranks share ",I0," GPUs on this node")') &
            & mdl%myid,node_size,mdl%ngpu
    endif
    if(mdl%mydev/=mod(node_rank,mdl%ngpu))then
       mdl%mydev=mod(node_rank,mdl%ngpu)
       call check_cuda(cudaSetDevice(mdl%mydev),'cudaSetDevice')
    endif
#endif
    if(mdl%myid==1)write(*,'(" Launching CUDA with ndevice/task = ",I0)')mdl%ngpu
    call check_cuda(cudaGetDeviceProperties(prop,mdl%mydev),'cudaGetDeviceProperties')
    write(*,'(" Rank ",I0," uses device ",I0,": ",A)')mdl%myid,mdl%mydev,trim(prop%name)
    write(*,'(" Compute Capability: ",i0,".",i0)')prop%major,prop%minor
    write(*,'(" Shared Memory per Block: ",i0)')prop%sharedMemPerBlock
    write(*,'(" Maximum Shared Memory per Block: ",i0)')prop%sharedMemPerBlockOptIn
#endif
#ifdef _METAL
    call mtl_init()
#endif
#ifdef _CUDA
  contains
    subroutine check_cuda(status,operation)
      integer,intent(in)::status
      character(len=*),intent(in)::operation
      if(status/=cudaSuccess)then
         write(*,*)'Error: ',operation,' failed on rank',mdl%myid,': ',trim(cudaGetErrorString(status))
         call mdl_abort(mdl)
      endif
    end subroutine check_cuda
#endif
  end subroutine mdl_initialize
  !##############################################################
  !##############################################################
  !##############################################################
  !##############################################################
  double precision function mdl_wtime(mdl)
    type(mdl_t)::mdl
    real::tt
    call cpu_time(tt)
    mdl_wtime=tt
  end function mdl_wtime
  !##############################################################
  !##############################################################
  !##############################################################
  !##############################################################
  integer function mdl_threads(mdl)
    type(mdl_t)::mdl
    mdl_threads = mdl%ncpu
  end function mdl_threads
  !##############################################################
  !##############################################################
  !##############################################################
  !##############################################################
  integer function mdl_self(mdl)
    type(mdl_t)::mdl
    mdl_self = mdl%myid
  end function mdl_self
  !##############################################################
  !##############################################################
  !##############################################################
  !##############################################################
  integer function mdl_core(mdl)
    type(mdl_t)::mdl
    mdl_core = 1
  end function mdl_core
  !##############################################################
  !##############################################################
  !##############################################################
  !##############################################################
  integer function mdl_cores(mdl)
    type(mdl_t)::mdl
    mdl_cores = 1
  end function mdl_cores
  !##############################################################
  !##############################################################
  !##############################################################
  !##############################################################
  subroutine mdl_mkdir(mdl,filedir)
    use amr_parameters, only: flen
    type(mdl_t)::mdl
    character(len=*), intent(in) :: filedir
    character(LEN=flen)::filecmd
    filecmd='mkdir -p '//TRIM(filedir)
#ifdef NOSYSTEM
    call PXFMKDIR(TRIM(filedir),LEN(TRIM(filedir)),O'755',ierr)
#else
    call system(filecmd)
#endif
  end subroutine mdl_mkdir
  !##############################################################
  !##############################################################
  !##############################################################
  !##############################################################
end module mdl_module
