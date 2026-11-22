#include <qmckl_f.F90>

program main
  use iso_c_binding
  use qmckl
  implicit none

  ! Declare variables
  integer :: argc
  character(128) :: trexio_filename, err_msg
  integer(qmckl_exit_code) :: rc
  integer(qmckl_context) :: context
  integer*8 :: nucl_num, mo_num, point_num
  double precision, allocatable :: nucl_coord(:,:)
  integer*8  :: nx(3)
  double precision, dimension(3) :: shift, step, rmin, rmax
  double precision, allocatable :: mo_value(:,:), overlap(:,:), point(:), linspace(:,:)
  double precision :: before, after, dr
  integer*8 :: i,j,k,m

  ! Initialize variables
  err_msg = ' ' 
  argc = command_argument_count()
  if (argc /= 1) then
    print *, "Usage: ./program <TREXIO filename>"
    stop -1
  end if
  call get_command_argument(1, trexio_filename)
  rc = QMCKL_SUCCESS

  ! Create a QMCkl context
  context = qmckl_context_create()

  ! Read the TREXIO file into the context
  rc = qmckl_trexio_read(context, trim(trexio_filename), len(trexio_filename)*1_8)
  if (rc /= QMCKL_SUCCESS) then
    call qmckl_string_of_error(rc, err_msg)
    write(*,*) "Error reading TREXIO file:", err_msg
    stop -1
  end if

  ! Retrieve the number of nuclei
  rc = qmckl_get_nucleus_num(context, nucl_num)
  if (rc /= QMCKL_SUCCESS) then
    call qmckl_string_of_error(rc, err_msg)
    write(*,*) "Error getting nucl_num:", err_msg
    stop -1
  end if

  ! Retrieve the nuclear coordinates
  allocate(nucl_coord(3, nucl_num))
  rc = qmckl_get_nucleus_coord(context, 'N', nucl_coord, nucl_num * 3_8)
  if (rc /= QMCKL_SUCCESS) then
    call qmckl_string_of_error(rc, err_msg)
    write(*,*) "Error getting nucl_coord:", err_msg
    stop -1
  end if

  ! Retrieve the number of MOs
  rc = qmckl_get_mo_basis_mo_num(context, mo_num)
  if (rc /= QMCKL_SUCCESS) then
    call qmckl_string_of_error(rc, err_msg)
    write(*,*) "Error getting mo_num:", err_msg
    stop -1
  end if

  ! Initialize grid points for the calculation
  nx = (/ 120, 120, 120 /)
  shift = (/ 5.d0, 5.d0, 5.d0 /)
  point_num = nx(1) * nx(2) * nx(3)

  ! Initialize rmin and rmax
  rmin = nucl_coord(:,1)
  rmax = nucl_coord(:,1)

  ! Update rmin and rmax based on nucl_coord
  do i = 1, 3
    do j = 1, nucl_num
      rmin(i) = min(nucl_coord(i,j), rmin(i))
      rmax(i) = max(nucl_coord(i,j), rmax(i))
    end do
  end do

  ! Apply shift
  rmin = rmin - shift
  rmax = rmax + shift

  ! Initialize linspace and step
  allocate(linspace(3, maxval(nx)))

  do i = 1, 3
    step(i) = (rmax(i) - rmin(i)) / real(nx(i) - 1, 8)
    do j = 1, nx(i)
      linspace(i, j) = rmin(i) + (j - 1) * step(i)
    end do
  end do

  ! Initialize point array
  allocate(point(3 * point_num))
  m = 1
  do i = 1, nx(1)
    do j = 1, nx(2)
      do k = 1, nx(3)
        point(m) = linspace(1, i); m = m + 1
        point(m) = linspace(2, j); m = m + 1
        point(m) = linspace(3, k); m = m + 1
      end do
    end do
  end do

  deallocate(linspace)


  ! Set points in QMCKL context
  rc = qmckl_set_point(context, 'N', point_num, point, point_num * 3)
  if (rc /= QMCKL_SUCCESS) then
    call qmckl_string_of_error(rc, err_msg)
    write(*,*) "Error setting point:", err_msg
    stop -1
  end if




  ! Perform the actual calculation and measure the time taken
  call cpu_time(before)
  allocate(mo_value(point_num, mo_num))
  rc = qmckl_get_mo_basis_mo_value(context, mo_value, point_num * mo_num)
  if (rc /= QMCKL_SUCCESS) then
    call qmckl_string_of_error(rc, err_msg)
    write(*,*) "Error getting mo_value:", err_msg
    stop
  end if
  call cpu_time(after)

  write(*,*) "Number of MOs:", mo_num
  write(*,*) "Number of grid points:", point_num
  write(*,*) "Execution time:", (after - before), "seconds"

  ! Compute the overlap matrix
  dr = step(1) * step(2) * step(3)

  allocate(overlap(mo_num, mo_num))
  rc = qmckl_dgemm(context, 'N', 'T', mo_num, mo_num, point_num, dr, &
                   mo_value, mo_num, mo_value, mo_num, 0.d0, overlap, mo_num)

  ! Print the overlap matrix
  do i = 1, mo_num
    write(*,'(i4)', advance='no') i
    do j = 1, mo_num
      write(*,'(f8.4)', advance='no') overlap(i, j)
    end do
    write(*,*)
  end do

  ! Clean-up and exit
  deallocate(mo_value, overlap)
  rc = qmckl_context_destroy(context)
  if (rc /= QMCKL_SUCCESS) then
    call qmckl_string_of_error(rc, err_msg)
    write(*,*) "Error destroying context:", err_msg
    stop -1
  end if

end program main
