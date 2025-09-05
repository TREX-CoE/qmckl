! Compute alpha
!    :PROPERTIES:
!    :Name:     qmckl_compute_det_vgl_alpha
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!     #+NAME: qmckl_compute_det_vgl_alpha_args
!     | Variable         | Type                                                       | In/Out | Description                               |
!     |------------------+------------------------------------------------------------+--------+-------------------------------------------|
!     | ~context~        | ~qmckl_context~                                            | in     | Global state                              |
!     | ~det_num_alpha~  | ~int64_t~                                                  | in     | Number of determinants                    |
!     | ~walk_num~       | ~int64_t~                                                  | in     | Number of walkers                         |
!     | ~alpha_num~      | ~int64_t~                                                  | in     | Number of electrons                       |
!     | ~beta_num~       | ~int64_t~                                                  | in     | Number of electrons                       |
!     | ~elec_num~       | ~int64_t~                                                  | in     | Number of electrons                       |
!     | ~mo_index_alpha~ | ~int64_t[det_num_alpha][walk_num][alpha_num]~              | in     | MO indices for electrons                  |
!     | ~mo_num~         | ~int64_t~                                                  | in     | Number of MOs                             |
!     | ~mo_vgl~         | ~double[5][elec_num][mo_num]~                              | in     | Value, gradients and Laplacian of the MOs |
!     | ~det_vgl_alpha~  | ~double[det_num_alpha][walk_num][5][alpha_num][alpha_num]~ | out    | Value, gradients and Laplacian of the Det |


integer function qmckl_compute_det_vgl_alpha_f(context, &
     det_num_alpha, walk_num, alpha_num, beta_num, elec_num, &
     mo_index_alpha, mo_num, mo_vgl, det_vgl_alpha) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context)  , intent(in)  :: context
  integer*8, intent(in)             :: det_num_alpha
  integer*8, intent(in)             :: walk_num
  integer*8, intent(in)             :: alpha_num
  integer*8, intent(in)             :: beta_num
  integer*8, intent(in)             :: elec_num
  integer*8, intent(in)             :: mo_num
  integer*8, intent(in)             :: mo_index_alpha(alpha_num, walk_num, det_num_alpha)
  double precision, intent(in)      :: mo_vgl(mo_num, elec_num, 5)
  double precision, intent(inout)   :: det_vgl_alpha(alpha_num, alpha_num, 5, walk_num, det_num_alpha)
  integer*8 :: idet, iwalk, ielec, mo_id, imo

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (alpha_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  do idet  = 1, det_num_alpha
  do iwalk = 1, walk_num
    do ielec = 1, alpha_num
      do imo = 1, alpha_num
        mo_id = mo_index_alpha(imo,iwalk,idet)
        ! Value
        det_vgl_alpha(imo, ielec, 1, iwalk, idet) = mo_vgl(mo_id, ielec, 1)

        ! Grad_x
        det_vgl_alpha(imo, ielec, 2, iwalk, idet) = mo_vgl(mo_id, ielec, 2)

        ! Grad_y
        det_vgl_alpha(imo, ielec, 3, iwalk, idet) = mo_vgl(mo_id, ielec, 3)

        ! Grad_z
        det_vgl_alpha(imo, ielec, 4, iwalk, idet) = mo_vgl(mo_id, ielec, 4)

        ! Lap
        det_vgl_alpha(imo, ielec, 5, iwalk, idet) = mo_vgl(mo_id, ielec, 5)
      end do
    end do
  end do
  end do

end function qmckl_compute_det_vgl_alpha_f



! #+CALL: generate_c_interface(table=qmckl_compute_det_vgl_alpha_args,rettyp=get_value("CRetType"),fname="qmckl_compute_det_vgl_alpha"))

!  #+RESULTS:

integer(c_int32_t) function qmckl_compute_det_vgl_alpha &
    (context, &
     det_num_alpha, &
     walk_num, &
     alpha_num, &
     beta_num, &
     elec_num, &
     mo_index_alpha, &
     mo_num, &
     mo_vgl, &
     det_vgl_alpha) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: det_num_alpha
  integer (c_int64_t) , intent(in)  , value :: walk_num
  integer (c_int64_t) , intent(in)  , value :: alpha_num
  integer (c_int64_t) , intent(in)  , value :: beta_num
  integer (c_int64_t) , intent(in)  , value :: elec_num
  integer (c_int64_t) , intent(in)          :: mo_index_alpha(alpha_num,walk_num,det_num_alpha)
  integer (c_int64_t) , intent(in)  , value :: mo_num
  real    (c_double ) , intent(in)          :: mo_vgl(mo_num,elec_num,5)
  real    (c_double ) , intent(out)         :: det_vgl_alpha(alpha_num,alpha_num,5,walk_num,det_num_alpha)

  integer(c_int32_t), external :: qmckl_compute_det_vgl_alpha_f
  info = qmckl_compute_det_vgl_alpha_f &
         (context, &
     det_num_alpha, &
     walk_num, &
     alpha_num, &
     beta_num, &
     elec_num, &
     mo_index_alpha, &
     mo_num, &
     mo_vgl, &
     det_vgl_alpha)

end function qmckl_compute_det_vgl_alpha

! Compute beta
!    :PROPERTIES:
!    :Name:     qmckl_compute_det_vgl_beta
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!     #+NAME: qmckl_compute_det_vgl_beta_args
!     | Variable        | Type                                                    | In/Out | Description                               |
!     |-----------------+---------------------------------------------------------+--------+-------------------------------------------|
!     | ~context~       | ~qmckl_context~                                         | in     | Global state                              |
!     | ~det_num_beta~  | ~int64_t~                                               | in     | Number of determinants                    |
!     | ~walk_num~      | ~int64_t~                                               | in     | Number of walkers                         |
!     | ~alpha_num~     | ~int64_t~                                               | in     | Number of electrons                       |
!     | ~beta_num~      | ~int64_t~                                               | in     | Number of electrons                       |
!     | ~elec_num~      | ~int64_t~                                               | in     | Number of electrons                       |
!     | ~mo_index_beta~ | ~int64_t[det_num_beta][walk_num][beta_num]~             | in     | Number of electrons                       |
!     | ~mo_num~        | ~int64_t~                                               | in     | Number of MOs                             |
!     | ~mo_vgl~        | ~double[5][elec_num][mo_num]~                           | in     | Value, gradients and Laplacian of the MOs |
!     | ~det_vgl_beta~  | ~double[det_num_beta][walk_num][5][beta_num][beta_num]~ | out    | Value, gradients and Laplacian of the Det |


integer function qmckl_compute_det_vgl_beta_f(context, &
     det_num_beta, walk_num, alpha_num, beta_num, elec_num, &
     mo_index_beta, mo_num, mo_vgl, det_vgl_beta) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context)  , intent(in)  :: context
  integer*8, intent(in)             :: det_num_beta
  integer*8, intent(in)             :: walk_num
  integer*8, intent(in)             :: alpha_num
  integer*8, intent(in)             :: beta_num
  integer*8, intent(in)             :: elec_num
  integer*8, intent(in)             :: mo_num
  integer*8, intent(in)             :: mo_index_beta(beta_num,walk_num,det_num_beta)
  double precision, intent(in)      :: mo_vgl(mo_num, elec_num, 5)
  double precision, intent(inout)   :: det_vgl_beta(beta_num, beta_num, 5, walk_num, det_num_beta)
  integer*8 :: idet, iwalk, ielec, mo_id, imo

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (beta_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  do idet = 1, det_num_beta
  do iwalk = 1, walk_num
    do ielec = 1, beta_num
      do imo = 1, beta_num
        mo_id = mo_index_beta(imo, iwalk, idet)
        ! Value
        det_vgl_beta(imo, ielec, 1, iwalk, idet) = mo_vgl(mo_id, alpha_num + ielec, 1)

        ! Grad_x
        det_vgl_beta(imo, ielec, 2, iwalk, idet) = mo_vgl(mo_id, alpha_num + ielec, 2)

        ! Grad_y
        det_vgl_beta(imo, ielec, 3, iwalk, idet) = mo_vgl(mo_id, alpha_num + ielec, 3)

        ! Grad_z
        det_vgl_beta(imo, ielec, 4, iwalk, idet) = mo_vgl(mo_id, alpha_num + ielec, 4)

        ! Lap
        det_vgl_beta(imo, ielec, 5, iwalk, idet) = mo_vgl(mo_id, alpha_num + ielec, 5)
      end do
    end do
  end do
  end do

end function qmckl_compute_det_vgl_beta_f



! #+CALL: generate_c_interface(table=qmckl_compute_det_vgl_beta_args,rettyp=get_value("CRetType"),fname="qmckl_compute_det_vgl_beta"))

!  #+RESULTS:

integer(c_int32_t) function qmckl_compute_det_vgl_beta &
    (context, &
     det_num_beta, &
     walk_num, &
     alpha_num, &
     beta_num, &
     elec_num, &
     mo_index_beta, &
     mo_num, &
     mo_vgl, &
     det_vgl_beta) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: det_num_beta
  integer (c_int64_t) , intent(in)  , value :: walk_num
  integer (c_int64_t) , intent(in)  , value :: alpha_num
  integer (c_int64_t) , intent(in)  , value :: beta_num
  integer (c_int64_t) , intent(in)  , value :: elec_num
  integer (c_int64_t) , intent(in)          :: mo_index_beta(beta_num,walk_num,det_num_beta)
  integer (c_int64_t) , intent(in)  , value :: mo_num
  real    (c_double ) , intent(in)          :: mo_vgl(mo_num,elec_num,5)
  real    (c_double ) , intent(out)         :: det_vgl_beta(beta_num,beta_num,5,walk_num,det_num_beta)

  integer(c_int32_t), external :: qmckl_compute_det_vgl_beta_f
  info = qmckl_compute_det_vgl_beta_f &
         (context, &
     det_num_beta, &
     walk_num, &
     alpha_num, &
     beta_num, &
     elec_num, &
     mo_index_beta, &
     mo_num, &
     mo_vgl, &
     det_vgl_beta)

end function qmckl_compute_det_vgl_beta

! Compute alpha
!    :PROPERTIES:
!    :Name:     qmckl_compute_det_inv_matrix_alpha
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!    #+NAME: qmckl_det_inv_matrix_alpha_args
!     | Variable               | Type                                                       | In/Out | Description                                                  |
!     |------------------------+------------------------------------------------------------+--------+--------------------------------------------------------------|
!     | ~context~              | ~qmckl_context~                                            | in     | Global state                                                 |
!     | ~det_num_alpha~        | ~int64_t~                                                  | in     | Number of determinants                                       |
!     | ~walk_num~             | ~int64_t~                                                  | in     | Number of walkers                                            |
!     | ~alpha_num~            | ~int64_t~                                                  | in     | Number of electrons                                          |
!     | ~det_vgl_alpha~        | ~double[det_num_alpha][walk_num][5][alpha_num][alpha_num]~ | in     | determinant matrix Value, gradients and Laplacian of the MOs |
!     | ~det_value_alpha~      | ~double[det_num_alpha][walk_num]~                          | out    | value of determinant matrix                                  |
!     | ~det_adj_matrix_alpha~ | ~double[det_num_alpha][walk_num][alpha_num][alpha_num]~    | out    | adjoint of determinant matrix                                |
!     | ~det_inv_matrix_alpha~ | ~double[det_num_alpha][walk_num][alpha_num][alpha_num]~    | out    | inverse of determinant matrix                                |


integer function qmckl_compute_det_inv_matrix_alpha_f(context, &
     det_num_alpha, walk_num, alpha_num, det_vgl_alpha, det_value_alpha, det_adj_matrix_alpha, det_inv_matrix_alpha) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context)  , intent(in)  :: context
  integer*8, intent(in)             :: det_num_alpha
  integer*8, intent(in)             :: walk_num
  integer*8, intent(in)             :: alpha_num
  double precision, intent(in)      :: det_vgl_alpha(alpha_num, alpha_num, 5, walk_num, det_num_alpha)
  double precision, intent(inout)   :: det_value_alpha(walk_num, det_num_alpha)
  double precision, intent(inout)   :: det_adj_matrix_alpha(alpha_num, alpha_num, walk_num, det_num_alpha)
  double precision, intent(inout)   :: det_inv_matrix_alpha(alpha_num, alpha_num, walk_num, det_num_alpha)
  double precision,dimension(:,:),allocatable  :: matA
  double precision                  :: det_l
  integer*8 :: idet, iwalk, ielec, mo_id, imo, LDA, res, i, j

  allocate(matA(alpha_num, alpha_num))

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (det_num_alpha <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  if (alpha_num <= 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  LDA = alpha_num
  do idet = 1, det_num_alpha
     do iwalk = 1, walk_num
        ! Value
        matA(1:alpha_num,1:alpha_num) = &
             det_vgl_alpha(1:alpha_num, 1:alpha_num, 1, iwalk, idet)

        res = qmckl_adjugate(context,                 &
             alpha_num, matA, LDA,                    &
             det_adj_matrix_alpha(1:alpha_num, 1:alpha_num, iwalk, idet), &
             int(size(det_adj_matrix_alpha,1),8),     &
             det_l)

        det_inv_matrix_alpha(1:alpha_num, 1:alpha_num, iwalk, idet) = &
             (1.d0/det_l) * &
             det_adj_matrix_alpha(1:alpha_num, 1:alpha_num, iwalk, idet)

        det_value_alpha(iwalk, idet) = det_l
     end do
  end do

  deallocate(matA)
end function qmckl_compute_det_inv_matrix_alpha_f



! #+CALL: generate_c_interface(table=qmckl_det_inv_matrix_alpha_args,rettyp=get_value("CRetType"),fname="qmckl_compute_det_inv_matrix_alpha"))

! #+RESULTS:

integer(c_int32_t) function qmckl_compute_det_inv_matrix_alpha &
    (context, &
     det_num_alpha, &
     walk_num, &
     alpha_num, &
     det_vgl_alpha, &
     det_value_alpha, &
     det_adj_matrix_alpha, &
     det_inv_matrix_alpha) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: det_num_alpha
  integer (c_int64_t) , intent(in)  , value :: walk_num
  integer (c_int64_t) , intent(in)  , value :: alpha_num
  real    (c_double ) , intent(in)          :: det_vgl_alpha(alpha_num,alpha_num,5,walk_num,det_num_alpha)
  real    (c_double ) , intent(out)         :: det_value_alpha(walk_num,det_num_alpha)
  real    (c_double ) , intent(out)         :: det_adj_matrix_alpha(alpha_num,alpha_num,walk_num,det_num_alpha)
  real    (c_double ) , intent(out)         :: det_inv_matrix_alpha(alpha_num,alpha_num,walk_num,det_num_alpha)

  integer(c_int32_t), external :: qmckl_compute_det_inv_matrix_alpha_f
  info = qmckl_compute_det_inv_matrix_alpha_f &
         (context, &
     det_num_alpha, &
     walk_num, &
     alpha_num, &
     det_vgl_alpha, &
     det_value_alpha, &
     det_adj_matrix_alpha, &
     det_inv_matrix_alpha)

end function qmckl_compute_det_inv_matrix_alpha

! Compute beta
!    :PROPERTIES:
!    :Name:     qmckl_compute_det_inv_matrix_beta
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!     #+NAME: qmckl_det_inv_matrix_beta_args
!     | Variable              | Type                                                    | In/Out | Description                                                  |
!     |-----------------------+---------------------------------------------------------+--------+--------------------------------------------------------------|
!     | ~context~             | ~qmckl_context~                                         | in     | Global state                                                 |
!     | ~det_num_beta~        | ~int64_t~                                               | in     | Number of determinants                                       |
!     | ~walk_num~            | ~int64_t~                                               | in     | Number of walkers                                            |
!     | ~beta_num~            | ~int64_t~                                               | in     | Number of electrons                                          |
!     | ~det_vgl_beta~        | ~double[det_num_beta][walk_num][5][beta_num][beta_num]~ | in     | determinant matrix Value, gradients and Laplacian of the MOs |
!     | ~det_value_beta~      | ~double[det_num_beta][walk_num]~                        | out    | value of determinant matrix                                  |
!     | ~det_adj_matrix_beta~ | ~double[det_num_beta][walk_num][beta_num][beta_num]~    | out    | adjoint of determinant matrix                                |
!     | ~det_inv_matrix_beta~ | ~double[det_num_beta][walk_num][beta_num][beta_num]~    | out    | inverse of determinant matrix                                |


integer function qmckl_compute_det_inv_matrix_beta_f(context, &
     det_num_beta, walk_num, beta_num, det_vgl_beta, det_value_beta, det_adj_matrix_beta, det_inv_matrix_beta) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context)  , intent(in)  :: context
  integer*8, intent(in)             :: det_num_beta
  integer*8, intent(in)             :: walk_num
  integer*8, intent(in)             :: beta_num
  double precision, intent(in)      :: det_vgl_beta(beta_num, beta_num, 5, walk_num, det_num_beta)
  double precision, intent(inout)   :: det_value_beta(walk_num, det_num_beta)
  double precision, intent(inout)   :: det_adj_matrix_beta(beta_num, beta_num, walk_num, det_num_beta)
  double precision, intent(inout)   :: det_inv_matrix_beta(beta_num, beta_num, walk_num, det_num_beta)
  double precision,dimension(:,:),allocatable  :: matA
  double precision                  :: det_l
  integer*8 :: idet, iwalk, ielec, mo_id, imo, LDA, res

  allocate(matA(beta_num, beta_num))

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (det_num_beta <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  if (beta_num <= 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  LDA = beta_num
  do idet = 1, det_num_beta
     do iwalk = 1, walk_num
        ! Value
        matA(1:beta_num,1:beta_num) = &
             det_vgl_beta(1:beta_num, 1:beta_num, 1, iwalk, idet)

        res = qmckl_adjugate(context,                 &
             beta_num, matA, LDA,                    &
             det_adj_matrix_beta(1, 1, iwalk, idet), &
             int(size(det_adj_matrix_beta,1),8),     &
             det_l)

        det_inv_matrix_beta(1:beta_num, 1:beta_num, iwalk, idet) = &
             (1.d0/det_l) * &
             det_adj_matrix_beta(1:beta_num, 1:beta_num, iwalk, idet)

        det_value_beta(iwalk, idet) = det_l
     end do
  end do


  deallocate(matA)
end function qmckl_compute_det_inv_matrix_beta_f



! #+CALL: generate_c_interface(table=qmckl_det_inv_matrix_beta_args,rettyp=get_value("CRetType"),fname="qmckl_compute_det_inv_matrix_beta"))

! #+RESULTS:

integer(c_int32_t) function qmckl_compute_det_inv_matrix_beta &
    (context, &
     det_num_beta, &
     walk_num, &
     beta_num, &
     det_vgl_beta, &
     det_value_beta, &
     det_adj_matrix_beta, &
     det_inv_matrix_beta) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: det_num_beta
  integer (c_int64_t) , intent(in)  , value :: walk_num
  integer (c_int64_t) , intent(in)  , value :: beta_num
  real    (c_double ) , intent(in)          :: det_vgl_beta(beta_num,beta_num,5,walk_num,det_num_beta)
  real    (c_double ) , intent(out)         :: det_value_beta(walk_num,det_num_beta)
  real    (c_double ) , intent(out)         :: det_adj_matrix_beta(beta_num,beta_num,walk_num,det_num_beta)
  real    (c_double ) , intent(out)         :: det_inv_matrix_beta(beta_num,beta_num,walk_num,det_num_beta)

  integer(c_int32_t), external :: qmckl_compute_det_inv_matrix_beta_f
  info = qmckl_compute_det_inv_matrix_beta_f &
         (context, &
     det_num_beta, &
     walk_num, &
     beta_num, &
     det_vgl_beta, &
     det_value_beta, &
     det_adj_matrix_beta, &
     det_inv_matrix_beta)

end function qmckl_compute_det_inv_matrix_beta
