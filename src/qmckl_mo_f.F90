! Compute
!    :PROPERTIES:
!    :Name:     qmckl_compute_mo_basis_mo_value
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!     #+NAME: qmckl_mo_basis_mo_value_args
!     | Variable        | Type                        | In/Out | Description                                     |
!     |-----------------+-----------------------------+--------+-------------------------------------------------|
!     | ~context~       | ~qmckl_context~             | in     | Global state                                    |
!     | ~ao_num~        | ~int64_t~                   | in     | Number of AOs                                   |
!     | ~mo_num~        | ~int64_t~                   | in     | Number of MOs                                   |
!     | ~point_num~     | ~int64_t~                   | in     | Number of points                                |
!     | ~coefficient_t~ | ~double[mo_num][ao_num]~    | in     | Transpose of the AO to MO transformation matrix |
!     | ~ao_value~      | ~double[point_num][ao_num]~ | in     | Value of the AOs                                |
!     | ~mo_value~      | ~double[point_num][mo_num]~ | out    | Value of the MOs                                |


!     The matrix of AO values is very sparse, so we use a sparse-dense
!     matrix multiplication instead of a dgemm, as exposed in
!     https://dx.doi.org/10.1007/978-3-642-38718-0_14.




integer function qmckl_compute_mo_basis_mo_value_doc_f(context, &
     ao_num, mo_num, point_num, &
     coefficient_t, ao_value, mo_value) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in)  :: context
  integer*8             , intent(in)  :: ao_num, mo_num
  integer*8             , intent(in)  :: point_num
  double precision      , intent(in)  :: ao_value(ao_num,point_num)
  double precision      , intent(in)  :: coefficient_t(mo_num,ao_num)
  double precision      , intent(out) :: mo_value(mo_num,point_num)
  integer*8 :: j,k

  info = QMCKL_SUCCESS

  do j=1,point_num
     mo_value(:,j) = 0.d0
     do k=1,ao_num
        if (ao_value(k,j) == 0.d0) cycle
        mo_value(:,j) = mo_value(:,j) + coefficient_t(:,k) * ao_value(k,j)
     end do
  end do

end function qmckl_compute_mo_basis_mo_value_doc_f



! #+CALL: generate_c_interface(table=qmckl_mo_basis_mo_value_args,rettyp=get_value("CRetType"),fname="qmckl_compute_mo_basis_mo_value_doc"))

!  #+RESULTS:

integer(c_int32_t) function qmckl_compute_mo_basis_mo_value_doc &
    (context, ao_num, mo_num, point_num, coefficient_t, ao_value, mo_value) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: ao_num
  integer (c_int64_t) , intent(in)  , value :: mo_num
  integer (c_int64_t) , intent(in)  , value :: point_num
  real    (c_double ) , intent(in)          :: coefficient_t(ao_num,mo_num)
  real    (c_double ) , intent(in)          :: ao_value(ao_num,point_num)
  real    (c_double ) , intent(out)         :: mo_value(mo_num,point_num)

  integer(c_int32_t), external :: qmckl_compute_mo_basis_mo_value_doc_f
  info = qmckl_compute_mo_basis_mo_value_doc_f &
         (context, ao_num, mo_num, point_num, coefficient_t, ao_value, mo_value)

end function qmckl_compute_mo_basis_mo_value_doc

! Compute
!    :PROPERTIES:
!    :Name:     qmckl_compute_mo_basis_mo_vgl
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!     #+NAME: qmckl_mo_basis_mo_vgl_args
!     | Variable            | Type                           | In/Out | Description                                     |
!     |---------------------+--------------------------------+--------+-------------------------------------------------|
!     | ~context~           | ~qmckl_context~                | in     | Global state                                    |
!     | ~ao_num~            | ~int64_t~                      | in     | Number of AOs                                   |
!     | ~mo_num~            | ~int64_t~                      | in     | Number of MOs                                   |
!     | ~point_num~         | ~int64_t~                      | in     | Number of points                                |
!     | ~coefficient_t~     | ~double[mo_num][ao_num]~       | in     | Transpose of the AO to MO transformation matrix |
!     | ~ao_vgl~            | ~double[point_num][5][ao_num]~ | in     | Value, gradients and Laplacian of the AOs       |
!     | ~mo_vgl~            | ~double[point_num][5][mo_num]~ | out    | Value, gradients and Laplacian of the MOs       |


!     The matrix of AO values is very sparse, so we use a sparse-dense
!     matrix multiplication instead of a dgemm, as exposed in
!     https://dx.doi.org/10.1007/978-3-642-38718-0_14.




integer function qmckl_compute_mo_basis_mo_vgl_doc_f(context, &
     ao_num, mo_num, point_num, &
     coefficient_t, ao_vgl, mo_vgl) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in)  :: context
  integer*8             , intent(in)  :: ao_num, mo_num
  integer*8             , intent(in)  :: point_num
  double precision      , intent(in)  :: ao_vgl(ao_num,5,point_num)
  double precision      , intent(in)  :: coefficient_t(mo_num,ao_num)
  double precision      , intent(out) :: mo_vgl(mo_num,5,point_num)
  integer*8 :: i,j,k
  double precision :: c1, c2, c3, c4, c5

  info = QMCKL_SUCCESS

  do j=1,point_num
     mo_vgl(:,:,j) = 0.d0
     do k=1,ao_num
        if (ao_vgl(k,1,j) /= 0.d0) then
           c1 = ao_vgl(k,1,j)
           c2 = ao_vgl(k,2,j)
           c3 = ao_vgl(k,3,j)
           c4 = ao_vgl(k,4,j)
           c5 = ao_vgl(k,5,j)
           do i=1,mo_num
              mo_vgl(i,1,j) = mo_vgl(i,1,j) + coefficient_t(i,k) * c1
              mo_vgl(i,2,j) = mo_vgl(i,2,j) + coefficient_t(i,k) * c2
              mo_vgl(i,3,j) = mo_vgl(i,3,j) + coefficient_t(i,k) * c3
              mo_vgl(i,4,j) = mo_vgl(i,4,j) + coefficient_t(i,k) * c4
              mo_vgl(i,5,j) = mo_vgl(i,5,j) + coefficient_t(i,k) * c5
           end do
        end if
     end do
  end do

end function qmckl_compute_mo_basis_mo_vgl_doc_f



! #+CALL: generate_c_interface(table=qmckl_mo_basis_mo_vgl_args,rettyp=get_value("CRetType"),fname="qmckl_compute_mo_basis_mo_vgl_doc"))

!  #+RESULTS:

integer(c_int32_t) function qmckl_compute_mo_basis_mo_vgl_doc &
    (context, ao_num, mo_num, point_num, coefficient_t, ao_vgl, mo_vgl) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: ao_num
  integer (c_int64_t) , intent(in)  , value :: mo_num
  integer (c_int64_t) , intent(in)  , value :: point_num
  real    (c_double ) , intent(in)          :: coefficient_t(ao_num,mo_num)
  real    (c_double ) , intent(in)          :: ao_vgl(ao_num,5,point_num)
  real    (c_double ) , intent(out)         :: mo_vgl(mo_num,5,point_num)

  integer(c_int32_t), external :: qmckl_compute_mo_basis_mo_vgl_doc_f
  info = qmckl_compute_mo_basis_mo_vgl_doc_f &
         (context, ao_num, mo_num, point_num, coefficient_t, ao_vgl, mo_vgl)

end function qmckl_compute_mo_basis_mo_vgl_doc

! Compute
!    :PROPERTIES:
!    :Name:     qmckl_compute_mo_basis_mo_value_cusp
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!     #+NAME: qmckl_mo_basis_mo_value_cusp_args
!     | Variable        | Type                          | In/Out | Description                                     |
!     |-----------------+-------------------------------+--------+-------------------------------------------------|
!     | ~context~       | ~qmckl_context~               | in     | Global state                                    |
!     | ~nucl_num~      | ~int64_t~                     | in     | Number of nuclei                                |
!     | ~ao_num~        | ~int64_t~                     | in     | Number of AOs                                   |
!     | ~mo_num~        | ~int64_t~                     | in     | Number of MOs                                   |
!     | ~point_num~     | ~int64_t~                     | in     | Number of points                                |
!     | ~ao_nucl~       | ~int64_t[ao_num]~             | in     | Nucleus on which the AO is centered             |
!     | ~ao_ang_mom~    | ~int32_t[ao_num]~             | in     | Angular momentum of the shell                   |
!     | ~en_distance~   | ~double[point_num][nucl_num]~ | in     | Electron-nucleus distances                      |
!     | ~r_cusp~        | ~double[nucl_num]~            | in     | Cusp-adjustment radius                          |
!     | ~cusp_param~    | ~double[nucl_num][4][mo_num]~ | in     | Cusp-adjustment parameters                      |
!     | ~coefficient_t~ | ~double[mo_num][ao_num]~      | in     | Transpose of the AO to MO transformation matrix |
!     | ~ao_value~      | ~double[point_num][ao_num]~   | in     | Value of the AOs                                |
!     | ~mo_value~      | ~double[point_num][mo_num]~   | out    | Cusp correction for the values of the MOs       |





integer function qmckl_compute_mo_basis_mo_value_cusp_doc_f(context, &
     nucl_num, ao_num, mo_num, point_num, ao_nucl, ao_ang_mom, en_distance, &
     r_cusp, cusp_param, coefficient_t, ao_value, mo_value) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in)  :: context
  integer*8             , intent(in)  :: nucl_num, ao_num, mo_num, point_num
  integer*8             , intent(in)  :: ao_nucl(ao_num)
  integer*4             , intent(in)  :: ao_ang_mom(ao_num)
  double precision      , intent(in)  :: en_distance(nucl_num, point_num)
  double precision      , intent(in)  :: r_cusp(nucl_num)
  double precision      , intent(in)  :: cusp_param(mo_num, 4, nucl_num)
  double precision      , intent(in)  :: coefficient_t(mo_num, ao_num)
  double precision      , intent(in)  :: ao_value(ao_num, point_num)
  double precision      , intent(out) :: mo_value(mo_num, point_num)

  integer*8 :: i, j, k, inucl
  double precision :: r

  info = QMCKL_SUCCESS

  do i=1,point_num
     mo_value(:,i) = 0.d0
     do k=1,ao_num
        if (ao_value(k,i) == 0.d0) cycle
        inucl = ao_nucl(k)+1
        if ( (en_distance(inucl,i) < r_cusp(inucl)) .and. (ao_ang_mom(k) == 0) ) cycle
        mo_value(:,i) = mo_value(:,i) + coefficient_t(:,k) * ao_value(k,i)
     end do ! k

     do inucl=1,nucl_num
        r = en_distance(inucl,i)
        if (r > r_cusp(inucl)) cycle

        do j=1,mo_num
           mo_value(j,i) = mo_value(j,i) + &
                cusp_param(j,1,inucl) + r*(cusp_param(j,2,inucl) + r*(  &
                cusp_param(j,3,inucl) + r* cusp_param(j,4,inucl)     ))
        enddo
     enddo ! inucl
  enddo ! i

end function qmckl_compute_mo_basis_mo_value_cusp_doc_f



! #+CALL: generate_c_interface(table=qmckl_mo_basis_mo_value_cusp_args,rettyp=get_value("CRetType"),fname="qmckl_compute_mo_basis_mo_value_cusp_doc"))

!  #+RESULTS:

integer(c_int32_t) function qmckl_compute_mo_basis_mo_value_cusp_doc &
    (context, &
     nucl_num, &
     ao_num, &
     mo_num, &
     point_num, &
     ao_nucl, &
     ao_ang_mom, &
     en_distance, &
     r_cusp, &
     cusp_param, &
     coefficient_t, &
     ao_value, &
     mo_value) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  integer (c_int64_t) , intent(in)  , value :: ao_num
  integer (c_int64_t) , intent(in)  , value :: mo_num
  integer (c_int64_t) , intent(in)  , value :: point_num
  integer (c_int64_t) , intent(in)          :: ao_nucl(ao_num)
  integer (c_int32_t) , intent(in)          :: ao_ang_mom(ao_num)
  real    (c_double ) , intent(in)          :: en_distance(nucl_num,point_num)
  real    (c_double ) , intent(in)          :: r_cusp(nucl_num)
  real    (c_double ) , intent(in)          :: cusp_param(mo_num,4,nucl_num)
  real    (c_double ) , intent(in)          :: coefficient_t(ao_num,mo_num)
  real    (c_double ) , intent(in)          :: ao_value(ao_num,point_num)
  real    (c_double ) , intent(out)         :: mo_value(mo_num,point_num)

  integer(c_int32_t), external :: qmckl_compute_mo_basis_mo_value_cusp_doc_f
  info = qmckl_compute_mo_basis_mo_value_cusp_doc_f &
         (context, &
     nucl_num, &
     ao_num, &
     mo_num, &
     point_num, &
     ao_nucl, &
     ao_ang_mom, &
     en_distance, &
     r_cusp, &
     cusp_param, &
     coefficient_t, &
     ao_value, &
     mo_value)

end function qmckl_compute_mo_basis_mo_value_cusp_doc

! Compute
!    :PROPERTIES:
!    :Name:     qmckl_compute_mo_basis_mo_vgl
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!     #+NAME: qmckl_mo_basis_mo_vgl_cusp_args
!     | Variable        | Type                           | In/Out | Description                                     |
!     |-----------------+--------------------------------+--------+-------------------------------------------------|
!     | ~context~       | ~qmckl_context~                | in     | Global state                                    |
!     | ~nucl_num~      | ~int64_t~                      | in     | Number of nuclei                                |
!     | ~ao_num~        | ~int64_t~                      | in     | Number of AOs                                   |
!     | ~mo_num~        | ~int64_t~                      | in     | Number of MOs                                   |
!     | ~point_num~     | ~int64_t~                      | in     | Number of points                                |
!     | ~ao_nucl~       | ~int64_t[ao_num]~              | in     | Nucleus on which the AO is centered             |
!     | ~ao_ang_mom~    | ~int32_t[ao_num]~              | in     | Angular momentum of the shell                   |
!     | ~en_distance~   | ~double[point_num][nucl_num]~  | in     | Electron-nucleus distances                      |
!     | ~nucl_coord~    | ~double[3][nucl_num]~          | in     | Nuclear coordinates                             |
!     | ~point_coord~   | ~double[3][point_num]~         | in     | Electron coordinates                            |
!     | ~r_cusp~        | ~double[nucl_num]~             | in     | Cusp-adjustment radius                          |
!     | ~cusp_param~    | ~double[nucl_num][4][mo_num]~  | in     | Cusp-adjustment parameters                      |
!     | ~coefficient_t~ | ~double[mo_num][ao_num]~       | in     | Transpose of the AO to MO transformation matrix |
!     | ~ao_vgl~        | ~double[point_num][5][ao_num]~ | in     | Value, gradients and Laplacian of the AOs       |
!     | ~mo_vgl~        | ~double[point_num][5][mo_num]~ | out    | Value, gradients and Laplacian of the MOs       |




integer function qmckl_compute_mo_basis_mo_vgl_cusp_doc_f(context, &
     nucl_num, ao_num, mo_num, point_num, ao_nucl, ao_ang_mom, en_distance, &
     nucl_coord, point_coord, r_cusp, cusp_param, coefficient_t, ao_vgl, mo_vgl) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in)  :: context
  integer*8             , intent(in)  :: nucl_num, ao_num, mo_num, point_num
  integer*8             , intent(in)  :: ao_nucl(ao_num)
  integer*4             , intent(in)  :: ao_ang_mom(ao_num)
  double precision      , intent(in)  :: en_distance(nucl_num, point_num)
  double precision      , intent(in)  :: nucl_coord(nucl_num,3)
  double precision      , intent(in)  :: point_coord(point_num,3)
  double precision      , intent(in)  :: r_cusp(nucl_num)
  double precision      , intent(in)  :: cusp_param(mo_num,4,nucl_num)
  double precision      , intent(in)  :: coefficient_t(mo_num,ao_num)
  double precision      , intent(in)  :: ao_vgl(ao_num,5,point_num)
  double precision      , intent(out) :: mo_vgl(mo_num,5,point_num)
  integer*8 :: i,j,k, inucl
  double precision :: c1, c2, c3, c4, c5
  double precision :: r, r_inv, r_vec(3)

  do j=1,point_num

     ! Initial contribution of the MO
     mo_vgl(:,:,j) = 0.d0
     do k=1,ao_num
        if (ao_vgl(k,1,j) /= 0.d0) then
          inucl = ao_nucl(k)+1
          if ( (en_distance(inucl,j) > r_cusp(inucl)) .or. (ao_ang_mom(k) > 0) ) then
             c1 = ao_vgl(k,1,j)
             c2 = ao_vgl(k,2,j)
             c3 = ao_vgl(k,3,j)
             c4 = ao_vgl(k,4,j)
             c5 = ao_vgl(k,5,j)
             do i=1,mo_num
                mo_vgl(i,1,j) = mo_vgl(i,1,j) + coefficient_t(i,k) * c1
                mo_vgl(i,2,j) = mo_vgl(i,2,j) + coefficient_t(i,k) * c2
                mo_vgl(i,3,j) = mo_vgl(i,3,j) + coefficient_t(i,k) * c3
                mo_vgl(i,4,j) = mo_vgl(i,4,j) + coefficient_t(i,k) * c4
                mo_vgl(i,5,j) = mo_vgl(i,5,j) + coefficient_t(i,k) * c5
             end do
          end if
       end if
     end do

    ! Cusp adjustment
    do inucl=1,nucl_num
       r = en_distance(inucl,j)
       if (r < r_cusp(inucl)) then
          
          r_vec(1:3) = point_coord(j,1:3) - nucl_coord(inucl,1:3)
          r_inv = 1.d0/r
          
          do i=1,mo_num
             mo_vgl(i,1,j) = mo_vgl(i,1,j) +  &
                  cusp_param(i,1,inucl) + r*( &
                  cusp_param(i,2,inucl) + r*( &
                  cusp_param(i,3,inucl) + r*( &
                  cusp_param(i,4,inucl)    )))
             
             c1 = r_inv * cusp_param(i,2,inucl) + 2.d0*cusp_param(i,3,inucl) +  &
                  r * 3.d0 * cusp_param(i,4,inucl)
             
             mo_vgl(i,2,j) = mo_vgl(i,2,j) + r_vec(1) * c1
             mo_vgl(i,3,j) = mo_vgl(i,3,j) + r_vec(2) * c1
             mo_vgl(i,4,j) = mo_vgl(i,4,j) + r_vec(3) * c1
             
             mo_vgl(i,5,j) = mo_vgl(i,5,j) +         &
                  2.d0*cusp_param(i,2,inucl)*r_inv + &
                  6.d0*cusp_param(i,3,inucl) +       &
                  12.d0*cusp_param(i,4,inucl)*r

          end do
       end if
    end do ! inucl
  end do
  info = QMCKL_SUCCESS

end function qmckl_compute_mo_basis_mo_vgl_cusp_doc_f



! #+CALL: generate_c_interface(table=qmckl_mo_basis_mo_vgl_cusp_args,rettyp=get_value("CRetType"),fname="qmckl_compute_mo_basis_mo_vgl_cusp_doc"))

!  #+RESULTS:

integer(c_int32_t) function qmckl_compute_mo_basis_mo_vgl_cusp_doc &
    (context, &
     nucl_num, &
     ao_num, &
     mo_num, &
     point_num, &
     ao_nucl, &
     ao_ang_mom, &
     en_distance, &
     nucl_coord, &
     point_coord, &
     r_cusp, &
     cusp_param, &
     coefficient_t, &
     ao_vgl, &
     mo_vgl) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  integer (c_int64_t) , intent(in)  , value :: ao_num
  integer (c_int64_t) , intent(in)  , value :: mo_num
  integer (c_int64_t) , intent(in)  , value :: point_num
  integer (c_int64_t) , intent(in)          :: ao_nucl(ao_num)
  integer (c_int32_t) , intent(in)          :: ao_ang_mom(ao_num)
  real    (c_double ) , intent(in)          :: en_distance(nucl_num,point_num)
  real    (c_double ) , intent(in)          :: nucl_coord(nucl_num,3)
  real    (c_double ) , intent(in)          :: point_coord(point_num,3)
  real    (c_double ) , intent(in)          :: r_cusp(nucl_num)
  real    (c_double ) , intent(in)          :: cusp_param(mo_num,4,nucl_num)
  real    (c_double ) , intent(in)          :: coefficient_t(ao_num,mo_num)
  real    (c_double ) , intent(in)          :: ao_vgl(ao_num,5,point_num)
  real    (c_double ) , intent(out)         :: mo_vgl(mo_num,5,point_num)

  integer(c_int32_t), external :: qmckl_compute_mo_basis_mo_vgl_cusp_doc_f
  info = qmckl_compute_mo_basis_mo_vgl_cusp_doc_f &
         (context, &
     nucl_num, &
     ao_num, &
     mo_num, &
     point_num, &
     ao_nucl, &
     ao_ang_mom, &
     en_distance, &
     nucl_coord, &
     point_coord, &
     r_cusp, &
     cusp_param, &
     coefficient_t, &
     ao_vgl, &
     mo_vgl)

end function qmckl_compute_mo_basis_mo_vgl_cusp_doc
