!################################################################
!################################################################
!################################################################
!################################################################
module posix_io
  !--------------------------------------------------------------
  ! Bulk binary output at explicit byte offsets.
  ! The file is created and opened as a Fortran stream unit by the
  ! caller (headers are written with Fortran I/O). With -DPOSIXIO
  ! the bulk data is then written with POSIX pwrite, bypassing the
  ! Fortran runtime (NVHPC splits every write into 8 KiB calls).
  ! Without it, the same bytes are written with Fortran stream I/O.
  !--------------------------------------------------------------
  use iso_c_binding
  use amr_parameters, only: flen
  implicit none
  private
  public::pio_t,pio_open,pio_write,pio_write_array,pio_close

  type pio_t
     integer::ilun=-1
     integer(c_int)::fd=-1
  end type pio_t

  ! Write a whole 1D array and advance the offset
  interface pio_write_array
     module procedure pio_write_r4,pio_write_r8,pio_write_i4,pio_write_i8
  end interface pio_write_array

#ifdef POSIXIO
  integer(c_int),parameter::O_WRONLY=1

  interface
     integer(c_int) function c_open(path,flags) bind(C,name='open')
       import::c_int,c_char
       character(kind=c_char),dimension(*)::path
       integer(c_int),value::flags
     end function c_open
     integer(c_long) function c_pwrite(fd,buf,nbyte,offset) bind(C,name='pwrite')
       import::c_int,c_long,c_size_t,c_ptr
       integer(c_int),value::fd
       type(c_ptr),value::buf
       integer(c_size_t),value::nbyte
       integer(c_long),value::offset
     end function c_pwrite
     integer(c_int) function c_close(fd) bind(C,name='close')
       import::c_int
       integer(c_int),value::fd
     end function c_close
     ! Print the message, then the system error text of errno
     subroutine c_perror(s) bind(C,name='perror')
       import::c_char
       character(kind=c_char),dimension(*)::s
     end subroutine c_perror
  end interface
#endif

contains
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pio_open(f,ilun)
  !--------------------------------------------------------------
  ! Attach to the file already opened as Fortran stream unit ilun
  !--------------------------------------------------------------
  type(pio_t)::f
  integer::ilun
#ifdef POSIXIO
  character(LEN=flen)::fileloc
#endif

  f%ilun=ilun
  ! Push buffered header bytes to the file first
  flush(ilun)
#ifdef POSIXIO
  inquire(unit=ilun,name=fileloc)
  f%fd=c_open(TRIM(fileloc)//c_null_char,O_WRONLY)
  if(f%fd<0)then
     call c_perror(' pio_open'//c_null_char)
     write(*,*)'pio_open: cannot open ',TRIM(fileloc)
     stop
  endif
#endif

end subroutine pio_open
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pio_write(f,buf,nbyte,offset)
  !--------------------------------------------------------------
  ! Write nbyte bytes starting at address buf to the file,
  ! at the 0-based byte offset (stream POS=offset+1)
  !--------------------------------------------------------------
  type(pio_t)::f
  type(c_ptr)::buf
  integer(kind=8)::nbyte,offset

  integer(kind=1),pointer,dimension(:)::bytes
#ifdef POSIXIO
  integer(c_long)::nw
  integer(kind=8)::ndone
  character(LEN=flen)::fileloc
#endif

  if(nbyte<=0)return
  call c_f_pointer(buf,bytes,[nbyte])
#ifdef POSIXIO
  ! pwrite may write fewer bytes than requested
  ndone=0
  do while(ndone<nbyte)
     nw=c_pwrite(f%fd,c_loc(bytes(ndone+1)),int(nbyte-ndone,c_size_t),int(offset+ndone,c_long))
     if(nw<=0)then
        ! perror first, before any Fortran I/O can change errno
        if(nw<0)call c_perror(' pio_write'//c_null_char)
        inquire(unit=f%ilun,name=fileloc)
        write(*,*)'pio_write: pwrite failed in ',TRIM(fileloc)
        stop
     endif
     ndone=ndone+nw
  end do
#else
  write(f%ilun,POS=offset+1)bytes
#endif

end subroutine pio_write
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pio_write_r4(f,x,offset)
  type(pio_t)::f
  real(kind=4),dimension(:),contiguous,target::x
  integer(kind=8)::offset
  integer(kind=8)::nbyte

  if(size(x)==0)return
  nbyte=4*int(size(x),kind=8)
  call pio_write(f,c_loc(x),nbyte,offset)
  offset=offset+nbyte

end subroutine pio_write_r4
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pio_write_r8(f,x,offset)
  type(pio_t)::f
  real(kind=8),dimension(:),contiguous,target::x
  integer(kind=8)::offset
  integer(kind=8)::nbyte

  if(size(x)==0)return
  nbyte=8*int(size(x),kind=8)
  call pio_write(f,c_loc(x),nbyte,offset)
  offset=offset+nbyte

end subroutine pio_write_r8
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pio_write_i4(f,x,offset)
  type(pio_t)::f
  integer(kind=4),dimension(:),contiguous,target::x
  integer(kind=8)::offset
  integer(kind=8)::nbyte

  if(size(x)==0)return
  nbyte=4*int(size(x),kind=8)
  call pio_write(f,c_loc(x),nbyte,offset)
  offset=offset+nbyte

end subroutine pio_write_i4
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pio_write_i8(f,x,offset)
  type(pio_t)::f
  integer(kind=8),dimension(:),contiguous,target::x
  integer(kind=8)::offset
  integer(kind=8)::nbyte

  if(size(x)==0)return
  nbyte=8*int(size(x),kind=8)
  call pio_write(f,c_loc(x),nbyte,offset)
  offset=offset+nbyte

end subroutine pio_write_i8
!################################################################
!################################################################
!################################################################
!################################################################
subroutine pio_close(f)
  type(pio_t)::f
#ifdef POSIXIO
  integer(c_int)::ierr

  ierr=c_close(f%fd)
  f%fd=-1
#endif
  f%ilun=-1

end subroutine pio_close

end module posix_io
