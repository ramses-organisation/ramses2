module lightcone_io_module
  use ramses_commons, only: pst_t
  use amr_parameters, only: flen, i8b
#ifndef WITHOUTMPI
  use mpi
#endif
  implicit none

contains

  subroutine open_lightcone_file(ilun, filename)
    integer :: ilun
    character(LEN=flen), intent(in) :: filename
#ifndef WITHOUTMPI
    integer :: ierr
#endif

#ifndef WITHOUTMPI
    call MPI_FILE_OPEN(MPI_COMM_WORLD, TRIM(filename), IOR(MPI_MODE_WRONLY,MPI_MODE_CREATE), MPI_INFO_NULL, ilun, ierr)
    call MPI_FILE_SET_SIZE(ilun, 0_MPI_OFFSET_KIND, ierr)
#else
    open(ilun, file=TRIM(filename), form='unformatted', access='stream', status='replace')
#endif
  end subroutine open_lightcone_file

  subroutine close_lightcone_file(ilun)
    integer :: ilun
#ifndef WITHOUTMPI
    integer :: ierr
#endif

#ifndef WITHOUTMPI
    call MPI_FILE_CLOSE(ilun, ierr)
#else
    close(ilun)
#endif
  end subroutine close_lightcone_file

  subroutine write_buffer(ilun, buffer, nbefore, ntotal, nthbuffer)
    use lightcone_buffer_module

    type(lightcone_buffer), intent(in) :: buffer
    integer, intent(in) :: ilun, nthbuffer
    integer(kind=8), intent(in) :: nbefore, ntotal
    integer :: idim
    integer(kind=8) :: offset
#ifndef WITHOUTMPI
    integer :: ierr, status(MPI_STATUS_SIZE)
#endif

    if (.not. buffer_is_empty(buffer)) then

       ! Write particle IDs (property 1)
       offset = calculate_write_offset(nbefore, ntotal, 1, nthbuffer, buffer%nstride)
#ifndef WITHOUTMPI
       call MPI_FILE_WRITE_AT(ilun, int(offset, kind=MPI_OFFSET_KIND), buffer%idp(1:buffer%ncurrent), buffer%ncurrent, MPI_INTEGER4, status, ierr)
#else
       write(unit=ilun, pos=offset+1) buffer%idp(1:buffer%ncurrent) ! pos is 1-based
#endif

       ! Write positions (properties 2, 3, 4)
       do idim = 1, 3
          offset = calculate_write_offset(nbefore, ntotal, idim+1, nthbuffer, buffer%nstride)
#ifndef WITHOUTMPI
          call MPI_FILE_WRITE_AT(ilun, int(offset, kind=MPI_OFFSET_KIND), buffer%xp(1:buffer%ncurrent, idim), buffer%ncurrent, MPI_REAL, status, ierr)
#else
          write(unit=ilun, pos=offset+1) buffer%xp(1:buffer%ncurrent, idim) ! pos is 1-based
#endif
       end do

       ! Write velocities (properties 5, 6, 7)
       do idim = 1, 3
          offset = calculate_write_offset(nbefore, ntotal, idim+4, nthbuffer, buffer%nstride)
#ifndef WITHOUTMPI
          call MPI_FILE_WRITE_AT(ilun, int(offset, kind=MPI_OFFSET_KIND), buffer%vp(1:buffer%ncurrent, idim), buffer%ncurrent, MPI_REAL, status, ierr)
#else
          write(unit=ilun, pos=offset+1) buffer%vp(1:buffer%ncurrent, idim) ! pos is 1-based
#endif
       end do

       ! Write masses (property 8)
       offset = calculate_write_offset(nbefore, ntotal, 8, nthbuffer, buffer%nstride)
#ifndef WITHOUTMPI
       call MPI_FILE_WRITE_AT(ilun, int(offset, kind=MPI_OFFSET_KIND), buffer%mp(1:buffer%ncurrent), buffer%ncurrent, MPI_REAL, status, ierr)
#else
       write(unit=ilun, pos=offset+1) buffer%mp(1:buffer%ncurrent) ! pos is 1-based
#endif
    end if
  end subroutine write_buffer

  subroutine write_buffer_grav(ilun, buffer, nbefore, ntotal, nthbuffer)
    use lightcone_buffer_module

    type(lightcone_buffer_grav), intent(in) :: buffer
    integer, intent(in) :: ilun, nthbuffer
    integer(kind=8), intent(in) :: nbefore, ntotal
    integer :: idim
    integer(kind=8) :: offset
#ifndef WITHOUTMPI
    integer :: ierr, status(MPI_STATUS_SIZE)
#endif

    if (.not. buffer_grav_is_empty(buffer)) then

       ! Write positions (properties 1, 2, 3)
       do idim = 1, 3
          offset = calculate_write_offset(nbefore, ntotal, idim, nthbuffer, buffer%nstride)
#ifndef WITHOUTMPI
          call MPI_FILE_WRITE_AT(ilun, int(offset, kind=MPI_OFFSET_KIND), buffer%xc(1:buffer%ncurrent, idim), buffer%ncurrent, MPI_REAL, status, ierr)
#else
          write(unit=ilun, pos=offset+1) buffer%xc(1:buffer%ncurrent, idim) ! pos is 1-based
#endif
       end do

       ! Write rho (property 4)
       offset = calculate_write_offset(nbefore, ntotal, 4, nthbuffer, buffer%nstride)
#ifndef WITHOUTMPI
       call MPI_FILE_WRITE_AT(ilun, int(offset, kind=MPI_OFFSET_KIND), buffer%rho(1:buffer%ncurrent), buffer%ncurrent, MPI_REAL, status, ierr)
#else
       write(unit=ilun, pos=offset+1) buffer%rho(1:buffer%ncurrent) ! pos in 1-based
#endif

       ! Write phi (property 5)
       offset = calculate_write_offset(nbefore, ntotal, 5, nthbuffer, buffer%nstride)
#ifndef WITHOUTMPI
       call MPI_FILE_WRITE_AT(ilun, int(offset, kind=MPI_OFFSET_KIND), buffer%phi(1:buffer%ncurrent), buffer%ncurrent, MPI_REAL, status, ierr)
#else
       write(unit=ilun, pos=offset+1) buffer%phi(1:buffer%ncurrent) ! pos is 1-based
#endif

       ! Write acceleration (properties 6, 7, 8)
       do idim = 1, 3
          offset = calculate_write_offset(nbefore, ntotal, idim+5, nthbuffer, buffer%nstride)
#ifndef WITHOUTMPI
          call MPI_FILE_WRITE_AT(ilun, int(offset, kind=MPI_OFFSET_KIND), buffer%accel(1:buffer%ncurrent, idim), buffer%ncurrent, MPI_REAL, status, ierr)
#else
          write(unit=ilun, pos=offset+1) buffer%accel(1:buffer%ncurrent, idim) ! pos is 1-based
#endif
       end do

       ! Write dphidt (property 9)
       offset = calculate_write_offset(nbefore, ntotal, 9, nthbuffer, buffer%nstride)
#ifndef WITHOUTMPI
       call MPI_FILE_WRITE_AT(ilun, int(offset, kind=MPI_OFFSET_KIND), buffer%dphidt(1:buffer%ncurrent), buffer%ncurrent, MPI_REAL, status, ierr)
#else
       write(unit=ilun, pos=offset+1) buffer%dphidt(1:buffer%ncurrent) ! pos is 1-based
#endif
    end if
  end subroutine write_buffer_grav

  subroutine write_lightcone_txt_file(filename, ntotal, aexp_old, aexp)
    character(LEN=flen), intent(in) :: filename
    integer(kind=8), intent(in) :: ntotal
    real(kind=8), intent(in) :: aexp_old, aexp
    integer :: ilun

    ilun = 1171 ! random number

    open(ilun, file=TRIM(filename)//".txt", form='formatted')
    rewind(ilun)
    write(ilun, *) ntotal
    write(ilun, *) aexp_old
    write(ilun, *) aexp
    close(ilun)
  end subroutine write_lightcone_txt_file

  subroutine output_lightcone_parameters(pst)
    type(pst_t), intent(in) :: pst

    write(*, *) ""
    write(*, *) "Lightcone parameters"
    write(*, *) "  opening_angle_y: ", pst%s%r%cone_opening_angle_y
    write(*, *) "  opening_angle_z: ", pst%s%r%cone_opening_angle_z
    write(*, *) "  cone_x_theta: ", pst%s%r%cone_theta
    write(*, *) "  cone_x_phi: ", pst%s%r%cone_phi
    write(*, *) "  z_max: ", pst%s%r%cone_z_max
    write(*, *) "  z_min: ", pst%s%r%cone_z_min
    write(*, *) "  observer position: ", pst%s%r%cone_observer
    write(*, *) ""
  end subroutine output_lightcone_parameters

  integer(kind=8) function calculate_write_offset(nbefore, ntotal, nthproperty, nthbuffer, nstride)
    ! Calculate the offset at which to write the buffer
    ! nbefore is the number of particles in the preceding processes
    ! ntotal is the total number of particles across all processes
    ! nthproperty is the property number (idp=1, x=2, y=3, z=4, vx=5, vy=6, vz=7, mass=8)
    ! nthbuffer is the buffer number (i.e if current process has already written 2 buffers, nthbuffer=3)

    ! All arithmetic is done in 64-bit: files larger than 2 GiB overflow 32-bit offsets
    integer(kind=8), intent(in) :: nbefore, ntotal
    integer, intent(in) :: nthproperty, nthbuffer, nstride
    integer(kind=8), parameter :: bpp = 4 ! bytes per property

    calculate_write_offset = bpp * ((nthproperty - 1) * ntotal + nbefore + int(nthbuffer - 1, kind=8) * nstride)
  end function calculate_write_offset

end module lightcone_io_module
