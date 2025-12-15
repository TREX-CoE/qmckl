! Compute

!     #+NAME: qmckl_nn_distance_args
!    | ~qmckl_context~ | ~context~                         | in  | Global state                   |
!    | ~int64_t~       | ~nucl_num~                        | in  | Number of nuclei               |
!    | ~double~        | ~coord[3][nucl_num]~              | in  | Nuclear coordinates (au)       |
!    | ~double~        | ~nn_distance[nucl_num][nucl_num]~ | out | Nucleus-nucleus distances (au) |


integer function qmckl_compute_nn_distance_f(context, nucl_num, coord, nn_distance) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in)  :: context
  integer*8             , intent(in)  :: nucl_num
  double precision      , intent(in)  :: coord(nucl_num,3)
  double precision      , intent(out) :: nn_distance(nucl_num,nucl_num)

  integer*8 :: k

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (nucl_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  info = qmckl_distance(context, 'T', 'T', nucl_num, nucl_num, &
          coord, nucl_num, &
          coord, nucl_num, &
          nn_distance, nucl_num)

end function qmckl_compute_nn_distance_f




! #+CALL: generate_c_interface(table=qmckl_nn_distance_args,rettyp="qmckl_exit_code",fname="qmckl_compute_nn_distance")

! #+RESULTS:

integer(c_int32_t) function qmckl_compute_nn_distance &
    (context, nucl_num, coord, nn_distance) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  real    (c_double ) , intent(in)          :: coord(nucl_num,3)
  real    (c_double ) , intent(out)         :: nn_distance(nucl_num,nucl_num)

  integer(c_int32_t), external :: qmckl_compute_nn_distance_f
  info = qmckl_compute_nn_distance_f &
         (context, nucl_num, coord, nn_distance)

end function qmckl_compute_nn_distance

! Compute

!     #+NAME: qmckl_nucleus_repulsion_args
!    | ~qmckl_context~ | ~context~                         | in  | Global state                   |
!    | ~int64_t~       | ~nucl_num~                        | in  | Number of nuclei               |
!    | ~double~        | ~charge[nucl_num]~                | in  | Nuclear charges (au)           |
!    | ~double~        | ~nn_distance[nucl_num][nucl_num]~ | in  | Nucleus-nucleus distances (au) |
!    | ~double~        | ~energy~                          | out | Nuclear repulsion energy       |


integer function qmckl_compute_nucleus_repulsion_f(context, nucl_num, charge, nn_distance, energy) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in)  :: context
  integer*8             , intent(in)  :: nucl_num
  double precision      , intent(in)  :: charge(nucl_num)
  double precision      , intent(in)  :: nn_distance(nucl_num,nucl_num)
  double precision      , intent(out) :: energy

  integer*8 :: i, j

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (nucl_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  energy = 0.d0
  do j=2, nucl_num
     do i=1, j-1
        if (dabs(nn_distance(i,j)) > 1e-5) then
          energy = energy + charge(i) * charge(j) / nn_distance(i,j)
        endif
     end do
  end do

end function qmckl_compute_nucleus_repulsion_f



! #+CALL: generate_c_interface(table=qmckl_nucleus_repulsion_args,rettyp="qmckl_exit_code",fname="qmckl_compute_nucleus_repulsion")

! #+RESULTS:

integer(c_int32_t) function qmckl_compute_nucleus_repulsion &
    (context, nucl_num, charge, nn_distance, energy) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  real    (c_double ) , intent(in)          :: charge(nucl_num)
  real    (c_double ) , intent(in)          :: nn_distance(nucl_num,nucl_num)
  real    (c_double ) , intent(out)         :: energy

  integer(c_int32_t), external :: qmckl_compute_nucleus_repulsion_f
  info = qmckl_compute_nucleus_repulsion_f &
         (context, nucl_num, charge, nn_distance, energy)

end function qmckl_compute_nucleus_repulsion
