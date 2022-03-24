subroutine qmckl_check_error(rc, message)
  use qmckl
  implicit none
  integer(qmckl_exit_code), intent(in) :: rc
  character(len=*)        , intent(in) :: message
  character(len=128)                   :: str_buffer
  if (rc /= QMCKL_SUCCESS) then
     print *, message
     call qmckl_string_of_error(rc, str_buffer)
     print *, str_buffer
     call exit(rc)
  end if
end subroutine qmckl_check_error

program ao_grid
  use qmckl
  implicit none

  integer(qmckl_context)    :: qmckl_ctx  ! QMCkl context
  integer(qmckl_exit_code)  :: rc         ! Exit code of QMCkl functions

  character(len=128)            :: trexio_filename
  character(len=128)            :: str_buffer
  integer                       :: ao_id
  integer                       :: point_num_x

  integer(c_int64_t)            :: nucl_num
  double precision, allocatable :: nucl_coord(:,:)

  integer(c_int64_t)            :: point_num
  integer(c_int64_t)            :: ao_num
  integer(c_int64_t)            :: ipoint, i, j, k
  double precision              :: x, y, z, dr(3)
  double precision              :: rmin(3), rmax(3)
  double precision, allocatable :: points(:,:)
  double precision, allocatable :: ao_vgl(:,:,:)

if (iargc() /= 3) then
   print *, 'Syntax: ao_grid <trexio_file> <AO_id> <point_num>'
   call exit(-1)
end if
call getarg(1, trexio_filename)
call getarg(2, str_buffer)
read(str_buffer, *) ao_id
call getarg(3, str_buffer)
read(str_buffer, *) point_num_x

if (point_num_x < 0 .or. point_num_x > 300) then
   print *, 'Error: 0 < point_num < 300'
   call exit(-1)
end if

qmckl_ctx = qmckl_context_create()
rc  = qmckl_trexio_read(qmckl_ctx, trexio_filename, 1_8*len(trim(trexio_filename)))
call qmckl_check_error(rc, 'Read TREXIO')

rc = qmckl_get_ao_basis_ao_num(qmckl_ctx, ao_num)
call qmckl_check_error(rc, 'Getting ao_num')

if (ao_id < 0 .or. ao_id > ao_num) then
   print *, 'Error: 0 < ao_id < ', ao_num
   call exit(-1)
end if

rc = qmckl_get_nucleus_num(qmckl_ctx, nucl_num)
call qmckl_check_error(rc, 'Get nucleus num')

allocate( nucl_coord(3, nucl_num) )
rc = qmckl_get_nucleus_coord(qmckl_ctx, 'N', nucl_coord, 3_8*nucl_num)
call qmckl_check_error(rc, 'Get nucleus coord')

rmin(1) = minval( nucl_coord(1,:) ) - 5.d0
rmin(2) = minval( nucl_coord(2,:) ) - 5.d0
rmin(3) = minval( nucl_coord(3,:) ) - 5.d0

rmax(1) = maxval( nucl_coord(1,:) ) + 5.d0
rmax(2) = maxval( nucl_coord(2,:) ) + 5.d0
rmax(3) = maxval( nucl_coord(3,:) ) + 5.d0

dr(1:3) = (rmax(1:3) - rmin(1:3)) / dble(point_num_x-1)

point_num = point_num_x**3
allocate( points(point_num, 3) )
ipoint=0
z = rmin(3)
do k=1,point_num_x
   y = rmin(2)
   do j=1,point_num_x
      x = rmin(1)
      do i=1,point_num_x
         ipoint = ipoint+1
         points(ipoint,1) = x
         points(ipoint,2) = y
         points(ipoint,3) = z
         x = x + dr(1)
      end do
      y = y + dr(2)
   end do
   z = z + dr(3)
end do

rc = qmckl_set_point(qmckl_ctx, 'T', points, point_num)
call qmckl_check_error(rc, 'Setting points')

allocate( ao_vgl(ao_num, 5, point_num) )
rc = qmckl_get_ao_basis_ao_vgl(qmckl_ctx, ao_vgl, ao_num*5_8*point_num)
call qmckl_check_error(rc, 'Setting points')

do ipoint=1, point_num
   print '(3(F16.10,X),E20.10)', points(ipoint, 1:3), ao_vgl(ao_id,1,ipoint)
end do

deallocate( nucl_coord, points, ao_vgl )
end program ao_grid
