! Compute
!    :PROPERTIES:
!    :Name:     qmckl_compute_ee_distance
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!     #+NAME: qmckl_ee_distance_args
!     | Variable      | Type                                   | In/Out | Description                 |
!     |---------------+----------------------------------------+--------+-----------------------------|
!     | ~context~     | ~qmckl_context~                        | in     | Global state                |
!     | ~elec_num~    | ~int64_t~                              | in     | Number of electrons         |
!     | ~walk_num~    | ~int64_t~                              | in     | Number of walkers           |
!     | ~coord~       | ~double[3][walk_num][elec_num]~        | in     | Electron coordinates        |
!     | ~ee_distance~ | ~double[walk_num][elec_num][elec_num]~ | out    | Electron-electron distances |


function qmckl_compute_ee_distance(context, elec_num, walk_num, coord, ee_distance) &
     result(info) bind(C)
  use qmckl_constants
  use qmckl, only : qmckl_distance
  implicit none
  integer (qmckl_context) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: elec_num
  integer (c_int64_t) , intent(in)  , value :: walk_num
  real    (c_double ) , intent(in)          :: coord(elec_num,3,walk_num)
  real    (c_double ) , intent(out)         :: ee_distance(elec_num,elec_num,walk_num)
  integer (qmckl_exit_code) :: info

  integer*8 :: k, i, j
  double precision :: x, y, z

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (elec_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  do k=1,walk_num
     info = qmckl_distance(context, 'T', 'T', elec_num, elec_num, &
          coord(1,k,1), elec_num * walk_num, &
          coord(1,k,1), elec_num * walk_num, &
          ee_distance(1,1,k), elec_num)
     if (info /= QMCKL_SUCCESS) then
        exit
     endif
  end do

end function qmckl_compute_ee_distance

! Compute
!    :PROPERTIES:
!    :Name:     qmckl_compute_ee_potential
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!     #+NAME: qmckl_ee_potential_args
!     | Variable       | Type                                   | In/Out | Description                 |
!     |----------------+----------------------------------------+--------+-----------------------------|
!     | ~context~      | ~qmckl_context~                        | in     | Global state                |
!     | ~elec_num~     | ~int64_t~                              | in     | Number of electrons         |
!     | ~walk_num~     | ~int64_t~                              | in     | Number of walkers           |
!     | ~ee_distance~  | ~double[walk_num][elec_num][elec_num]~ | in     | Electron-electron distances |
!     | ~ee_potential~ | ~double[walk_num]~                     | out    | Electron-electron potential |


function qmckl_compute_ee_potential(context, elec_num, walk_num, &
     ee_distance, ee_potential) result(info) bind(C)
  use qmckl_constants
  implicit none
  integer(qmckl_context), intent(in), value :: context
  integer (c_int64_t) , intent(in)  , value :: elec_num
  integer (c_int64_t) , intent(in)  , value :: walk_num
  real    (c_double ) , intent(in)          :: ee_distance(elec_num,elec_num,walk_num)
  real    (c_double ) , intent(out)         :: ee_potential(walk_num)
  integer (qmckl_exit_code) :: info

  integer*8 :: nw, i, j

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (elec_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  ee_potential = 0.0d0
  do nw=1,walk_num
    do j=2,elec_num
      do i=1,j-1
        if (dabs(ee_distance(i,j,nw)) > 1e-5) then
          ee_potential(nw) = ee_potential(nw) + 1.0d0/(ee_distance(i,j,nw))
        endif
      end do
    end do
  end do
  
end function qmckl_compute_ee_potential

! Compute
!    :PROPERTIES:
!    :Name:     qmckl_compute_en_distance
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!     #+NAME: qmckl_en_distance_args
!     | Variable      | Type                          | In/Out | Description                |
!     |---------------+-------------------------------+--------+----------------------------|
!     | ~context~     | ~qmckl_context~               | in     | Global state               |
!     | ~point_num~   | ~int64_t~                     | in     | Number of points           |
!     | ~nucl_num~    | ~int64_t~                     | in     | Number of nuclei           |
!     | ~elec_coord~  | ~double[3][point_num]~        | in     | Electron coordinates       |
!     | ~nucl_coord~  | ~double[3][nucl_num]~         | in     | Nuclear coordinates        |
!     | ~en_distance~ | ~double[point_num][nucl_num]~ | out    | Electron-nucleus distances |


function qmckl_compute_en_distance(context, &
     point_num, nucl_num, elec_coord, nucl_coord, en_distance) &
     result(info) bind(C)
  use qmckl_constants
  use qmckl, only : qmckl_distance
  implicit none
  integer (qmckl_context), intent(in), value :: context
  integer (c_int64_t)    , intent(in), value :: point_num
  integer (c_int64_t)    , intent(in), value :: nucl_num
  real    (c_double )    , intent(in)        :: elec_coord(point_num,3)
  real    (c_double )    , intent(in)        :: nucl_coord(nucl_num,3)
  real    (c_double )    , intent(out)       :: en_distance(nucl_num,point_num)
  integer(qmckl_exit_code) :: info
  
  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (point_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (nucl_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  info = qmckl_distance(context, 'T', 'T', nucl_num, point_num, &
          nucl_coord, nucl_num, &
          elec_coord, point_num, &
          en_distance, nucl_num)

end function qmckl_compute_en_distance

! Compute
!    :PROPERTIES:
!    :Name:     qmckl_compute_en_potential
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!     #+NAME: qmckl_en_potential_args
!     | Variable       | Type                                   | In/Out | Description                          |
!     |----------------+----------------------------------------+--------+--------------------------------------|
!     | ~context~      | ~qmckl_context~                        | in     | Global state                         |
!     | ~elec_num~     | ~int64_t~                              | in     | Number of electrons                  |
!     | ~nucl_num~     | ~int64_t~                              | in     | Number of nuclei                     |
!     | ~walk_num~     | ~int64_t~                              | in     | Number of walkers                    |
!     | ~charge~       | ~double[nucl_num]~                     | in     | charge of nucleus                    |
!     | ~en_distance~  | ~double[walk_num][elec_num][nucl_num]~ | in     | Electron-electron distances          |
!     | ~en_potential~ | ~double[walk_num]~                     | out    | Electron-electron potential          |


function qmckl_compute_en_potential(context, elec_num, nucl_num, walk_num, &
     charge, en_distance, en_potential) &
     result(info) bind(C)
  use qmckl
  implicit none
  integer (qmckl_context), intent(in), value :: context
  integer (c_int64_t) , intent(in)   , value :: elec_num
  integer (c_int64_t) , intent(in)   , value :: nucl_num
  integer (c_int64_t) , intent(in)   , value :: walk_num
  real    (c_double ) , intent(in)           :: charge(nucl_num)
  real    (c_double ) , intent(in)           :: en_distance(nucl_num,elec_num,walk_num)
  real    (c_double ) , intent(out)          :: en_potential(walk_num)

  integer(qmckl_exit_code) :: info
  integer*8 :: nw, i, j

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (elec_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  en_potential = 0.0d0
  do nw=1,walk_num
    do i=1,elec_num
      do j=1,nucl_num
        if (dabs(en_distance(j,i,nw)) > 1.d-6) then
          en_potential(nw) = en_potential(nw) - charge(j)/(en_distance(j,i,nw))
        endif
      end do
    end do
  end do

end function qmckl_compute_en_potential
