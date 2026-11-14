module lightcone_test_module
  use ramses_commons, only: pst_t
  use lightcone_utils
  use lightcone_io_module
  implicit none
contains

  subroutine lightcone_test(pst)
    ! This subroutine is used only during development to test the lightcone utility functions
    ! It is best called from adaptive_loop.f90 with a stop command right after
    type(pst_t), intent(in) :: pst
    integer :: i

    write(*,*) "==== STARTING LIGHTCONE TEST ===="

    ! Test that the namelist parameters are read correctly
    call output_lightcone_parameters(pst)

    ! Test the rotation matrices
    call test_rotation_matrices(pst)

    ! Test the comoving_distance function
    call test_comoving_distance(pst) ! Results can be compared against Wright cosmology calculator online

    ! Test the comoving2code function
    call test_comoving2code(pst)

    ! Test the calculate_write_offset function
    call test_calculate_write_offset(pst)

    write(*,*) "==== ENDING LIGHTCONE TEST ===="
  end subroutine lightcone_test

  subroutine test_rotation_matrices(pst)
    type(pst_t), intent(in) :: pst
    real(kind=8) :: x(3) ! Coordinates of a point in the some coordinate system
    real(kind=8) :: x_prime(3) ! Coordinates of the same point in a different coordinate system
    real(kind=8) :: x_back(3) ! Coordinates the same point in the original coordinate system
    real(kind=8) :: cone_to_box_rotation(3,3), box_to_cone_rotation(3,3)

    cone_to_box_rotation = rotation_matrix(deg2rad(pst%s%r%cone_theta), deg2rad(pst%s%r%cone_phi))
    box_to_cone_rotation = transpose(cone_to_box_rotation)

    write(*,*) ""
    write(*,*) "Testing rotation matrices"
    write(*,*) "  R (cone to box):"
    call print_dp_matrix(cone_to_box_rotation)
    write(*,*) ""
    write(*,*) "  R^-1 (box to cone):"
    call print_dp_matrix(box_to_cone_rotation)
    write(*,*) ""
    write(*,*) "  R * R^-1:"
    call print_dp_matrix(matmul(cone_to_box_rotation, box_to_cone_rotation))
    write(*,*) ""

    write(*,*) ""
    write(*,*) "Testing cone to box coordinates (cone -> box -> cone)"
    ! Test case 0: Origin
    x = (/0.0, 0.0, 0.0/)
    x_prime = cone_to_box_coordinates(cone_to_box_rotation, pst%s%r%cone_observer, x)
    x_back = box_to_cone_coordinates(box_to_cone_rotation, pst%s%r%cone_observer, x_prime)
    write(*, '(3F12.2, A, 3F12.2, A, 3F12.2)') x, '  -> ', x_prime, '  -> ', x_back

    ! Test case 1: Along x-axis
    x = (/1.0, 0.0, 0.0/)
    x_prime = cone_to_box_coordinates(cone_to_box_rotation, pst%s%r%cone_observer, x)
    x_back = box_to_cone_coordinates(box_to_cone_rotation, pst%s%r%cone_observer, x_prime)
    write(*, '(3F12.2, A, 3F12.2, A, 3F12.2)') x, '  -> ', x_prime, '  -> ', x_back

    ! Test case 2: Along y-axis
    x = (/0.0, 1.0, 0.0/)
    x_prime = cone_to_box_coordinates(cone_to_box_rotation, pst%s%r%cone_observer, x)
    x_back = box_to_cone_coordinates(box_to_cone_rotation, pst%s%r%cone_observer, x_prime)
    write(*, '(3F12.2, A, 3F12.2, A, 3F12.2)') x, '  -> ', x_prime, '  -> ', x_back

    ! Test case 3: Along z-axis
    x = (/0.0, 0.0, 1.0/)
    x_prime = cone_to_box_coordinates(cone_to_box_rotation, pst%s%r%cone_observer, x)
    x_back = box_to_cone_coordinates(box_to_cone_rotation, pst%s%r%cone_observer, x_prime)
    write(*, '(3F12.2, A, 3F12.2, A, 3F12.2)') x, '  -> ', x_prime, '  -> ', x_back
    write(*,*) ""

  end subroutine test_rotation_matrices

  subroutine test_comoving_distance(pst)
    type(pst_t), intent(in) :: pst
    real(kind=8) :: z, dist, Omega0, OmegaL, OmegaR, h0, coverH0

    Omega0 = pst%s%g%omega_m
    OmegaL = pst%s%g%omega_l
    OmegaR = 1.0d0 - Omega0 - OmegaL
    h0 = pst%s%g%h0
    coverH0 = 2.9979246d+5/h0

    write(*,*) ""
    write(*,*) "Testing comoving_distance with the following cosmology:"
    write(*,*) "  omega_m = ", Omega0
    write(*,*) "  omega_l = ", OmegaL
    write(*,*) "  omega_r = ", OmegaR
    write(*,*) "  h0 = ", h0
    write(*,*) ""

    z = 0.5d0
    dist = comoving_distance(z, Omega0, OmegaL, OmegaR, coverH0)
    write(*,*) "    Distance at z = ", z, " is ", dist, " Mpc"

    z = 1.0d0
    dist = comoving_distance(z, Omega0, OmegaL, OmegaR, coverH0)
    write(*,*) "    Distance at z = ", z, " is ", dist, " Mpc"

    z = 2.0d0
    dist = comoving_distance(z, Omega0, OmegaL, OmegaR, coverH0)
    write(*,*) "    Distance at z = ", z, " is ", dist, " Mpc"

    z = 10.0d0
    dist = comoving_distance(z, Omega0, OmegaL, OmegaR, coverH0)
    write(*,*) "    Distance at z = ", z, " is ", dist, " Mpc"

    z = 100.0d0
    dist = comoving_distance(z, Omega0, OmegaL, OmegaR, coverH0)
    write(*,*) "    Distance at z = ", z, " is ", dist, " Mpc"

    z = 1000.0d0
    dist = comoving_distance(z, Omega0, OmegaL, OmegaR, coverH0)
    write(*,*) "    Distance at z = ", z, " is ", dist, " Mpc"
  end subroutine test_comoving_distance

  subroutine test_comoving2code(pst)
    type(pst_t), intent(in) :: pst
    real(kind=8) :: l, code_l

    write(*,*) ""
    write(*,*) "Testing conversion of comoving distance to code units"
    l = 1.0d0
    code_l = comoving2code(pst%s%g, l) ! This tells us how many code units per Mpc, the inverse of which is how many Mpc per code unit
    write(*,*) "    Box comoving size is ", 1.0d0 / code_l, " Mpc"
    write(*,*) ""
  end subroutine test_comoving2code

  subroutine print_dp_matrix(matrix)
    real(kind=8), dimension(:,:), intent(in) :: matrix
    integer :: i, n, m
    character(len=32) :: fmt

    ! Get matrix dimensions
    n = size(matrix, 1)
    m = size(matrix, 2)

    ! Create format string for m columns
    write(fmt,'(A,I0,A)') '(',m,'F12.2)'

    ! Loop over rows
    do i = 1, n
       write(*,fmt) matrix(i, :)
    end do

  end subroutine print_dp_matrix

  subroutine test_calculate_write_offset(pst)
    type(pst_t), intent(in) :: pst
    integer(kind=8) :: nbefore, ntotal
    integer :: nthproperty, nthbuffer, nstride
    integer(kind=8) :: offset

    nbefore = 80000
    ntotal = 100000
    nthproperty = 3
    nthbuffer = 2
    nstride = 20000
    offset = calculate_write_offset(nbefore, ntotal, nthproperty, nthbuffer, nstride)

    write(*,*) ""
    write(*,*) "Testing calculate_write_offset"
    write(*,*) "    nbefore = ", nbefore
    write(*,*) "    ntotal = ", ntotal
    write(*,*) "    nthproperty = ", nthproperty
    write(*,*) "    nthbuffer = ", nthbuffer
    write(*,*) "    nstride = ", nstride
    write(*,*) "    Offset is ", offset
    write(*,*) ""
  end subroutine test_calculate_write_offset

end module lightcone_test_module
