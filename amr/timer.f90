!################################################################
!################################################################
!################################################################
!################################################################
!=======================================================================
real(kind=8) function wallclock()
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none
#ifdef WITHOUTMPI
  integer(kind=8), save :: tstart
  integer(kind=8)       :: tcur
  integer(kind=8)       :: count_rate
#else
  real(kind=8), save :: tstart
  real(kind=8)       :: tcur
#endif
  logical,      save :: first_call=.true.
  real(kind=8), save :: norm
  !---------------------------------------------------------------------
  if (first_call) then
#ifdef WITHOUTMPI
     call system_clock(count=tstart, count_rate=count_rate)
     norm=1d0/count_rate
#else
     norm = 1d0
     tstart = MPI_Wtime()
#endif
     first_call=.false.
  end if
#ifdef WITHOUTMPI
  call system_clock(count=tcur)
#else
  tcur = MPI_Wtime()
#endif
  wallclock = (tcur-tstart)*norm
end function wallclock
!################################################################
!################################################################
!################################################################
!################################################################
module timer_module
  implicit none
  integer,            parameter         :: mtimer=200       ! max nr of timers
  real(kind=8),       dimension(mtimer) :: start, time
  integer                               :: ntimer=0, itimer=0
  character(len=72), dimension(mtimer)  :: labels
contains
  subroutine findit (label)
    implicit none
    character(len=*) label
    do itimer=1,ntimer
       if (trim(label) == trim(labels(itimer))) return
    end do
    ntimer = ntimer+1
    itimer = ntimer
    labels(itimer) = label
    time(itimer) = 0.
  end subroutine findit
end module timer_module
!################################################################
!################################################################
!################################################################
!################################################################
subroutine m_timer(label,cmd)
  use mdl_module
  use timer_module
#ifdef _CUDA
  use cudafor
#elif defined(_METAL)
  use metal_interface
#endif
  implicit none
  character(len=*) label, cmd
  real(kind=8) wallclock, current
#ifdef _CUDA
  integer :: cuda_ierr
  cuda_ierr = cudaDeviceSynchronize()
#elif defined(_METAL)
  call mtl_device_sync()
#endif
  current = wallclock()                                                 ! current time
  if (itimer > 0) then                                                  ! if timer is active ..
     time(itimer) = time(itimer) + current - start(itimer)              ! add to it
  end if
  call findit (label)                                                   ! locate timer slot
  if (cmd == 'start') then                                              ! start command
     start(itimer) = current                                            ! register start time
  else if (cmd == 'stop') then                                          ! stop command
     itimer = 0                                                         ! turn off timer
  end if
end subroutine m_timer
!################################################################
!################################################################
!################################################################
!################################################################
subroutine m_output_timer(write_file,filename)
  use amr_parameters, only: flen
  use mdl_module
  use timer_module
  use color_module, only: c_bold, c_reset
  implicit none
  real(kind=8)::total
  character(len=*)::filename
  character(len=flen)::fileloc
  logical::write_file
  integer::ilun,jtimer
  character(len=:),allocatable::cb,cr
  if(write_file)then
     ilun=11
     fileloc=TRIM(filename) ! Open file for timing info
     open(unit=ilun,file=fileloc,form='formatted')
  else
     ilun=6
  endif
  if(write_file)then
     cb=''; cr=''
  else
     cb=c_bold(); cr=c_reset()
  endif
  write (ilun,'(/a)') cb//'     seconds         %    STEP'//cr
  total = 1e-9
  do jtimer = 1,ntimer
     total = total + time(jtimer)
  end do
  do jtimer = 1,ntimer
     if (time(jtimer)/total >= 0.001) write (ilun,'(f12.3,4x,f6.1,4x,a24)') &
          time(jtimer), 100.*time(jtimer)/total,labels(jtimer)
  end do
  write (ilun,'(a,f12.3,4x,f6.1,4x,a)') cb, total, 100., 'TOTAL'//cr
  write(ilun,*)
  if(write_file)close(ilun)
end subroutine m_output_timer
!################################################################
!################################################################
!################################################################
!################################################################
