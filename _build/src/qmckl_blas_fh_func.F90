! C interface                                                    :noexport:

!     #+CALL: generate_f_interface(table=qmckl_dgemm_args,rettyp="qmckl_exit_code",fname="qmckl_dgemm")

!     #+RESULTS:

interface
  integer(qmckl_exit_code) function qmckl_dgemm &
      (context, TransA, TransB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (qmckl_context), intent(in)  , value :: context
    character(c_char  ) , intent(in)  , value :: TransA
    character(c_char  ) , intent(in)  , value :: TransB
    integer (c_int64_t) , intent(in)  , value :: m
    integer (c_int64_t) , intent(in)  , value :: n
    integer (c_int64_t) , intent(in)  , value :: k
    real    (c_double ) , intent(in)  , value :: alpha
    integer (c_int64_t) , intent(in)  , value :: lda
    integer (c_int64_t) , intent(in)  , value :: ldb
    real    (c_double ) , intent(in)  , value :: beta
    integer (c_int64_t) , intent(in)  , value :: ldc
    real    (c_double ) , intent(in)          :: A(lda,*)
    real    (c_double ) , intent(in)          :: B(ldb,*)
    real    (c_double ) , intent(out)         :: C(ldc,*)

  end function qmckl_dgemm
end interface

! C interface                                                    :noexport:

!     #+CALL: generate_f_interface(table=qmckl_dgemm_safe_args,rettyp="qmckl_exit_code",fname="qmckl_dgemm_safe")

!     #+RESULTS:

interface
  integer(qmckl_exit_code) function qmckl_dgemm_safe &
      (context, TransA, TransB, m, n, k, alpha, A, size_max_A, lda, B, size_max_B, ldb, beta, C, size_max_C, ldc) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (qmckl_context), intent(in)  , value :: context
    character(c_char  ) , intent(in)  , value :: TransA
    character(c_char  ) , intent(in)  , value :: TransB
    integer (c_int64_t) , intent(in)  , value :: m
    integer (c_int64_t) , intent(in)  , value :: n
    integer (c_int64_t) , intent(in)  , value :: k
    real    (c_double ) , intent(in)  , value :: alpha
    integer (c_int64_t) , intent(in)  , value :: size_max_A
    integer (c_int64_t) , intent(in)  , value :: lda
    integer (c_int64_t) , intent(in)  , value :: size_max_B
    integer (c_int64_t) , intent(in)  , value :: ldb
    real    (c_double ) , intent(in)  , value :: beta
    integer (c_int64_t) , intent(in)  , value :: size_max_C
    integer (c_int64_t) , intent(in)  , value :: ldc
    real    (c_double ) , intent(in)          :: A(lda,*)
    real    (c_double ) , intent(in)          :: B(ldb,*)
    real    (c_double ) , intent(out)         :: C(ldc,*)

  end function qmckl_dgemm_safe
end interface



! #+CALL: generate_f_interface(table=qmckl_adjugate_args,rettyp="qmckl_exit_code",fname="qmckl_adjugate")

! #+RESULTS:

interface
  integer(qmckl_exit_code) function qmckl_adjugate &
      (context, n, A, lda, B, ldb, det_l) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (qmckl_context), intent(in)  , value :: context
    integer (c_int64_t) , intent(in)  , value :: n
    integer (c_int64_t) , intent(in)  , value :: lda
    integer (c_int64_t) , intent(in)  , value :: ldb
    real    (c_double ) , intent(inout)        :: det_l
    real    (c_double ) , intent(in)          :: A(lda,*)
    real    (c_double ) , intent(out)         :: B(ldb,*)

  end function qmckl_adjugate
end interface

! C interface

!     #+CALL: generate_f_interface(table=qmckl_adjugate_safe_args,rettyp="qmckl_exit_code",fname="qmckl_adjugate_safe")

!     #+RESULTS:

interface
  integer(qmckl_exit_code) function qmckl_adjugate_safe &
      (context, n, A, size_max_A, lda, B, size_max_B, ldb, det_l) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (qmckl_context), intent(in)  , value :: context
    integer (c_int64_t) , intent(in)  , value :: n
    integer (c_int64_t) , intent(in)  , value :: size_max_A
    integer (c_int64_t) , intent(in)  , value :: lda
    integer (c_int64_t) , intent(in)  , value :: size_max_B
    integer (c_int64_t) , intent(in)  , value :: ldb
    real    (c_double ) , intent(inout)        :: det_l
    real    (c_double ) , intent(in)          :: A(lda,*)
    real    (c_double ) , intent(out)         :: B(ldb,*)

  end function qmckl_adjugate_safe
end interface

! C interface                                                    :noexport:

!     #+CALL: generate_f_interface(table=qmckl_dgemv_args,rettyp="qmckl_exit_code",fname="qmckl_dgemv")

!     #+RESULTS:

interface
  integer(qmckl_exit_code) function qmckl_dgemv &
      (context, Trans, m, n, alpha, A, lda, X, incx, beta, Y, incy) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (qmckl_context), intent(in)  , value :: context
    character(c_char  ) , intent(in)  , value :: Trans
    integer (c_int64_t) , intent(in)  , value :: m
    integer (c_int64_t) , intent(in)  , value :: n
    real    (c_double ) , intent(in)  , value :: alpha
    integer (c_int64_t) , intent(in)  , value :: lda
    integer (c_int64_t) , intent(in)  , value :: incx
    real    (c_double ) , intent(in)  , value :: beta
    integer (c_int64_t) , intent(in)  , value :: incy
    real    (c_double ) , intent(in)          :: A(lda,*)
    real    (c_double ) , intent(in)          :: X(*)
    real    (c_double ) , intent(out)         :: Y(*)

  end function qmckl_dgemv
end interface

! C interface                                                    :noexport:

!     #+CALL: generate_f_interface(table=qmckl_dger_args,rettyp="qmckl_exit_code",fname="qmckl_dger")

!     #+RESULTS:

interface
  integer(qmckl_exit_code) function qmckl_dger &
      (context, m, n, alpha, X, incx, Y, incy, A, lda) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (qmckl_context), intent(in)  , value :: context
    integer (c_int64_t) , intent(in)  , value :: m
    integer (c_int64_t) , intent(in)  , value :: n
    real    (c_double ) , intent(in)  , value :: alpha
    integer (c_int64_t) , intent(in)  , value :: lda
    integer (c_int64_t) , intent(in)  , value :: incx
    integer (c_int64_t) , intent(in)  , value :: incy
    real    (c_double ) , intent(in)          :: X(*)
    real    (c_double ) , intent(in)          :: Y(*)
    real    (c_double ) , intent(inout)       :: A(lda,*)

  end function qmckl_dger
end interface
