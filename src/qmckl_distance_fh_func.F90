! Performance

!     This function is more efficient when ~A~ and ~B~ are
!     transposed.

!    #+CALL: generate_f_interface(table=qmckl_distance_sq_args,fname=get_value("Name"))

!    #+RESULTS:

interface
  integer(qmckl_exit_code) function qmckl_distance_sq &
      (context, transa, transb, m, n, A, lda, B, ldb, C, ldc) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (qmckl_context), intent(in)  , value :: context
    character(c_char  ) , intent(in)  , value :: transa
    character(c_char  ) , intent(in)  , value :: transb
    integer (c_int64_t) , intent(in)  , value :: m
    integer (c_int64_t) , intent(in)  , value :: n
    integer (c_int64_t) , intent(in)  , value :: lda
    integer (c_int64_t) , intent(in)  , value :: ldb
    integer (c_int64_t) , intent(in)  , value :: ldc
    real    (c_double ) , intent(in)          :: A(lda,*)
    real    (c_double ) , intent(in)          :: B(ldb,*)
    real    (c_double ) , intent(out)         :: C(ldc,n)

  end function qmckl_distance_sq
end interface



! #+CALL: generate_f_interface(table=qmckl_distance_args,fname="qmckl_distance")

! #+RESULTS:

interface
  integer(qmckl_exit_code) function qmckl_distance &
      (context, transa, transb, m, n, A, lda, B, ldb, C, ldc) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (qmckl_context), intent(in)  , value :: context
    character(c_char  ) , intent(in)  , value :: transa
    character(c_char  ) , intent(in)  , value :: transb
    integer (c_int64_t) , intent(in)  , value :: m
    integer (c_int64_t) , intent(in)  , value :: n
    integer (c_int64_t) , intent(in)  , value :: lda
    integer (c_int64_t) , intent(in)  , value :: ldb
    integer (c_int64_t) , intent(in)  , value :: ldc
    real    (c_double ) , intent(in)          :: A(lda,*)
    real    (c_double ) , intent(in)          :: B(ldb,*)
    real    (c_double ) , intent(out)         :: C(ldc,n)

  end function qmckl_distance
end interface

! C interface                                                     :noexport:

!    #+CALL: generate_f_interface(table=qmckl_distance_rescaled_args,fname="qmckl_distance_rescaled")

!    #+RESULTS:

interface
  integer(qmckl_exit_code) function qmckl_distance_rescaled &
      (context, transa, transb, m, n, A, lda, B, ldb, C, ldc, rescale_factor_kappa) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (qmckl_context), intent(in)  , value :: context
    character(c_char  ) , intent(in)  , value :: transa
    character(c_char  ) , intent(in)  , value :: transb
    integer (c_int64_t) , intent(in)  , value :: m
    integer (c_int64_t) , intent(in)  , value :: n
    integer (c_int64_t) , intent(in)  , value :: lda
    integer (c_int64_t) , intent(in)  , value :: ldb
    integer (c_int64_t) , intent(in)  , value :: ldc
    real    (c_double ) , intent(in)  , value :: rescale_factor_kappa
    real    (c_double ) , intent(in)          :: A(lda,*)
    real    (c_double ) , intent(in)          :: B(ldb,*)
    real    (c_double ) , intent(out)         :: C(ldc,n)

  end function qmckl_distance_rescaled
end interface



!  This function is more efficient when ~A~ and ~B~ are transposed.


! #+CALL: generate_f_interface(table=qmckl_distance_rescaled_gl_args,rettyp=get_value("FRetType"),fname=get_value("Name"))

! #+RESULTS:

interface
  integer(qmckl_exit_code) function qmckl_distance_rescaled_gl &
      (context, transa, transb, m, n, A, lda, B, ldb, C, ldc, rescale_factor_kappa) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (qmckl_context), intent(in)  , value :: context
    character(c_char  ) , intent(in)  , value :: transa
    character(c_char  ) , intent(in)  , value :: transb
    integer (c_int64_t) , intent(in)  , value :: m
    integer (c_int64_t) , intent(in)  , value :: n
    integer (c_int64_t) , intent(in)  , value :: lda
    integer (c_int64_t) , intent(in)  , value :: ldb
    integer (c_int64_t) , intent(in)  , value :: ldc
    real    (c_double ) , intent(in)  , value :: rescale_factor_kappa
    real    (c_double ) , intent(in)          :: A(lda,*)
    real    (c_double ) , intent(in)          :: B(ldb,*)
    real    (c_double ) , intent(out)         :: C(4,ldc,n)

  end function qmckl_distance_rescaled_gl
end interface
