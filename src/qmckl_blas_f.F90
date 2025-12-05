function qmckl_dgemm(context, TransA, TransB, &
     m, n, k, alpha, A, LDA, B, LDB, beta, C, LDC) &
     bind(C) result(info)
  use qmckl_constants
#ifdef HAVE_LIBQMCKLDGEMM
  use qmckl_dgemm_tiled_module
#endif
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
  integer (c_int64_t) , intent(in)  , value :: ldc
  real    (c_double ) , intent(in)  , value :: beta
  real    (c_double ) , intent(in)          :: A(lda,*)
  real    (c_double ) , intent(in)          :: B(ldb,*)
  real    (c_double ) , intent(out)         :: C(ldc,*)
  integer(qmckl_exit_code) :: info

#ifdef HAVE_LIBQMCKLDGEMM
  double precision,allocatable,dimension(:,:) :: A1
  double precision,allocatable,dimension(:,:) :: B1
  double precision,allocatable,dimension(:,:) :: C1
#endif

  integer*8                           :: i, j, LDA1, LDB1, LDC1

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (m <= 0_8) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  if (n <= 0_8) then
     info = QMCKL_INVALID_ARG_5
     return
  endif

  if (k <= 0_8) then
     info = QMCKL_INVALID_ARG_6
     return
  endif

  if (LDA <= 0) then
     info = QMCKL_INVALID_ARG_9
     return
  endif

  if (LDB <= 0) then
     info = QMCKL_INVALID_ARG_11
     return
  endif

  if (LDC <= 0) then
     info = QMCKL_INVALID_ARG_14
     return
  endif

#ifdef HAVE_LIBQMCKLDGEMM
  ! Copy A to A1
  allocate(A1(k,m))
  do j=1,m
     do i=1,k
        A1(i,j) = A(j,i)
     end do
  end do

  ! Copy B to B1
  allocate(B1(n,k))
  do j=1,k
     do i=1,n
        B1(i,j) = B(j,i)
     end do
  end do

  ! Copy C to C1
  allocate(C1(n,m))
  do j=1,m
     do i=1,n
        C1(i,j) = C(j,i)
     end do
  end do

  LDA1 = size(A1,1)
  LDB1 = size(B1,1)
  LDC1 = size(C1,1)

  info = qmckl_dgemm_tiled_avx2(int(m,8), int(n,8), int(k,8), &
       A1, int(LDA1,8), B1, int(LDB1,8), C1, int(LDC1,8))

  do j=1,n
     do i=1,m
        transpose         C(i,j) = alpha*C1(j,i) + beta*C(i,j)
     end do
  end do

  deallocate(A1,B1,C1)
#else
  call dgemm(transA, transB, int(m,4), int(n,4), int(k,4), &
       alpha, A, int(LDA,4), B, int(LDB,4), beta, C, int(LDC,4))
#endif


end function qmckl_dgemm

function qmckl_dgemm_safe(context, TransA, TransB, &
     m, n, k, alpha, A, size_A, LDA, B, size_B, LDB, beta, C, size_C, LDC) &
     result(info) bind(C)
  use qmckl_constants
  implicit none
  integer (qmckl_context), intent(in)  , value :: context
  character(c_char  ) , intent(in)  , value :: TransA
  character(c_char  ) , intent(in)  , value :: TransB
  integer (c_int64_t) , intent(in)  , value :: m
  integer (c_int64_t) , intent(in)  , value :: n
  integer (c_int64_t) , intent(in)  , value :: k
  real    (c_double ) , intent(in)  , value :: alpha
  integer (c_int64_t) , intent(in)  , value :: size_A
  integer (c_int64_t) , intent(in)  , value :: lda
  integer (c_int64_t) , intent(in)  , value :: size_B
  integer (c_int64_t) , intent(in)  , value :: ldb
  real    (c_double ) , intent(in)  , value :: beta
  integer (c_int64_t) , intent(in)  , value :: size_C
  integer (c_int64_t) , intent(in)  , value :: ldc
  real    (c_double ) , intent(in)          :: A(lda,*)
  real    (c_double ) , intent(in)          :: B(ldb,*)
  real    (c_double ) , intent(out)         :: C(ldc,*)

  integer(qmckl_exit_code) :: info
  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (m <= 0_8) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  if (n <= 0_8) then
     info = QMCKL_INVALID_ARG_5
     return
  endif

  if (k <= 0_8) then
     info = QMCKL_INVALID_ARG_6
     return
  endif

  if (LDA <= 0) then
     info = QMCKL_INVALID_ARG_10
     return
  endif

  if (LDB <= 0) then
     info = QMCKL_INVALID_ARG_13
     return
  endif

  if (LDC <= 0) then
     info = QMCKL_INVALID_ARG_17
     return
  endif

  if (size_A <= 0) then
     info = QMCKL_INVALID_ARG_9
     return
  endif

  if (size_B <= 0) then
     info = QMCKL_INVALID_ARG_12
     return
  endif

  if (size_C <= 0) then
     info = QMCKL_INVALID_ARG_16
     return
  endif

  call dgemm(transA, transB, int(m,4), int(n,4), int(k,4), &
       alpha, A, int(LDA,4), B, int(LDB,4), beta, C, int(LDC,4))

end function qmckl_dgemm_safe

function qmckl_adjugate(context, n, A, LDA, B, ldb, det_l) &
     result(info) bind(C)
  use qmckl_constants
  implicit none

  integer (qmckl_context), intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: n
  integer (c_int64_t) , intent(in)  , value :: lda
  integer (c_int64_t) , intent(in)  , value :: ldb
  real    (c_double ) , intent(inout)        :: det_l
  real    (c_double ) , intent(in)          :: A(lda,*)
  real    (c_double ) , intent(out)         :: B(ldb,*)
  integer(qmckl_exit_code) :: info

  info = QMCKL_SUCCESS

  select case (n)
  case(1)
     det_l = a(1,1)
     b(1,1) = 1.d0
  case (2)
     det_l  =  a(1,1)*a(2,2) - a(1,2)*a(2,1)
     b(1,1) =  a(2,2)
     b(2,1) = -a(2,1)
     b(1,2) = -a(1,2)
     b(2,2) =  a(1,1)
  case (3)
     det_l = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) &
            -a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1)) &
            +a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))

     b(1,1) =  a(2,2)*a(3,3) - a(2,3)*a(3,2)
     b(2,1) =  a(2,3)*a(3,1) - a(2,1)*a(3,3)
     b(3,1) =  a(2,1)*a(3,2) - a(2,2)*a(3,1)

     b(1,2) =  a(1,3)*a(3,2) - a(1,2)*a(3,3)
     b(2,2) =  a(1,1)*a(3,3) - a(1,3)*a(3,1)
     b(3,2) =  a(1,2)*a(3,1) - a(1,1)*a(3,2)

     b(1,3) =  a(1,2)*a(2,3) - a(1,3)*a(2,2)
     b(2,3) =  a(1,3)*a(2,1) - a(1,1)*a(2,3)
     b(3,3) =  a(1,1)*a(2,2) - a(1,2)*a(2,1)

  case (4)
     call adjugate4(A,LDA,B,LDB,n,det_l)
  case (5)
     call adjugate5(A,LDA,B,LDB,n,det_l)
  case default
     if (context == QMCKL_NULL_CONTEXT) then
        info = QMCKL_INVALID_CONTEXT
        return
     endif

     if (n <= 0_8) then
        info = QMCKL_INVALID_ARG_2
        return
     endif
    
     if (LDA < n) then
        info = QMCKL_INVALID_ARG_4
        return
     endif

     call adjugate_general(context, n, A, LDA, B, LDB, det_l)
  end select

end function qmckl_adjugate

subroutine cofactor2(A,LDA,B,LDB,na,det_l)
  implicit none
  integer*8, intent(in)           :: LDA, LDB
  integer*8, intent(in)           :: na
  double precision, intent(in)    :: A (LDA,na)
  double precision, intent(out)   :: B (LDB,na)
  double precision, intent(inout) :: det_l

  double precision :: C(2,2)

  call adjugate2(A,LDA,C,2_8,na,det_l)

  B(1,1) = C(1,1)
  B(2,1) = C(1,2)
  B(1,2) = C(2,1)
  B(2,2) = C(2,2)

end subroutine cofactor2

subroutine cofactor3(a,LDA,B,LDB,na,det_l)
  implicit none
  integer*8, intent(in)           :: LDA, LDB
  integer*8, intent(in)           :: na
  double precision, intent(in)    :: A (LDA,na)
  double precision, intent(out)   :: B (LDB,na)
  double precision, intent(inout) :: det_l

  double precision :: C(4,3)

  call adjugate3(A,LDA,C,4_8,na,det_l)

  B(1,1) = C(1,1)
  B(1,2) = C(2,1)
  B(1,3) = C(3,1)
  B(2,1) = C(1,2)
  B(2,2) = C(2,2)
  B(2,3) = C(3,2)
  B(3,1) = C(1,3)
  B(3,2) = C(2,3)
  B(3,3) = C(3,3)

end subroutine cofactor3

subroutine cofactor4(a,LDA,B,LDB,na,det_l)
  implicit none
  integer*8, intent(in)           :: LDA, LDB
  integer*8, intent(in)           :: na
  double precision, intent(in)    :: A (LDA,na)
  double precision, intent(out)   :: B (LDB,na)
  double precision, intent(inout) :: det_l

  double precision :: C(4,4)

  call adjugate4(A,LDA,C,4_8,na,det_l)
  B(1,1) = C(1,1)
  B(1,2) = C(2,1)
  B(1,3) = C(3,1)
  B(1,4) = C(4,1)
  B(2,1) = C(1,2)
  B(2,2) = C(2,2)
  B(2,3) = C(3,2)
  B(2,4) = C(4,2)
  B(3,1) = C(1,3)
  B(3,2) = C(2,3)
  B(3,3) = C(3,3)
  B(3,4) = C(4,3)
  B(4,1) = C(1,4)
  B(4,2) = C(2,4)
  B(4,3) = C(3,4)
  B(4,4) = C(4,4)

end subroutine cofactor4

subroutine cofactor5(A,LDA,B,LDB,na,det_l)
  implicit none
  integer*8, intent(in)           :: LDA, LDB
  integer*8, intent(in)           :: na
  double precision, intent(in)    :: A (LDA,na)
  double precision, intent(out)   :: B (LDB,na)
  double precision, intent(inout) :: det_l

  double precision  :: C(8,5)

  call adjugate5(A,LDA,C,8_8,na,det_l)

  B(1,1) = C(1,1)
  B(1,2) = C(2,1)
  B(1,3) = C(3,1)
  B(1,4) = C(4,1)
  B(1,5) = C(5,1)
  B(2,1) = C(1,2)
  B(2,2) = C(2,2)
  B(2,3) = C(3,2)
  B(2,4) = C(4,2)
  B(2,5) = C(5,2)
  B(3,1) = C(1,3)
  B(3,2) = C(2,3)
  B(3,3) = C(3,3)
  B(3,4) = C(4,3)
  B(3,5) = C(5,3)
  B(4,1) = C(1,4)
  B(4,2) = C(2,4)
  B(4,3) = C(3,4)
  B(4,4) = C(4,4)
  B(4,5) = C(5,4)
  B(5,1) = C(1,5)
  B(5,2) = C(2,5)
  B(5,3) = C(3,5)
  B(5,4) = C(4,5)
  B(5,5) = C(5,5)

end subroutine cofactor5

subroutine adjugate2(a,LDA,b,LDB,na,det_l)
  implicit none
  integer*8, intent(in)           :: LDA, LDB
  integer*8        :: na
  double precision, intent(in)    :: A (LDA,na)
  double precision, intent(out)   :: B (LDB,na)
  double precision :: det_l

  det_l  =  a(1,1)*a(2,2) - a(1,2)*a(2,1)
  b(1,1) =  a(2,2)
  b(2,1) = -a(2,1)
  b(1,2) = -a(1,2)
  b(2,2) =  a(1,1)
end subroutine adjugate2

subroutine adjugate3(a,LDA,b,LDB,na,det_l)
  implicit none
  integer*8, intent(in)           :: LDA, LDB
  integer*8, intent(in)           :: na
  double precision, intent(in)    :: A (LDA,na)
  double precision, intent(out)   :: B (LDB,na)
  double precision, intent(inout) :: det_l
  integer :: i

  det_l = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) &
         -a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1)) &
         +a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))

  b(1,1) =  a(2,2)*a(3,3) - a(2,3)*a(3,2)
  b(2,1) =  a(2,3)*a(3,1) - a(2,1)*a(3,3)
  b(3,1) =  a(2,1)*a(3,2) - a(2,2)*a(3,1)

  b(1,2) =  a(1,3)*a(3,2) - a(1,2)*a(3,3)
  b(2,2) =  a(1,1)*a(3,3) - a(1,3)*a(3,1)
  b(3,2) =  a(1,2)*a(3,1) - a(1,1)*a(3,2)

  b(1,3) =  a(1,2)*a(2,3) - a(1,3)*a(2,2)
  b(2,3) =  a(1,3)*a(2,1) - a(1,1)*a(2,3)
  b(3,3) =  a(1,1)*a(2,2) - a(1,2)*a(2,1)

end subroutine adjugate3

subroutine adjugate4(a,LDA,b,LDB,na,det_l)
  implicit none
  integer*8, intent(in)           :: LDA, LDB
  integer*8, intent(in)           :: na
  double precision, intent(in)    :: A (LDA,na)
  double precision, intent(out)   :: B (LDB,na)
  double precision, intent(inout) :: det_l
  integer :: i,j
  det_l =  a(1,1)*(a(2,2)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))  &
                  -a(2,3)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))  &
                  +a(2,4)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))) &
          -a(1,2)*(a(2,1)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))  &
                  -a(2,3)*(a(3,1)*a(4,4)-a(3,4)*a(4,1))  &
                  +a(2,4)*(a(3,1)*a(4,3)-a(3,3)*a(4,1))) &
          +a(1,3)*(a(2,1)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))  &
                  -a(2,2)*(a(3,1)*a(4,4)-a(3,4)*a(4,1))  &
                  +a(2,4)*(a(3,1)*a(4,2)-a(3,2)*a(4,1))) &
          -a(1,4)*(a(2,1)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))  &
                  -a(2,2)*(a(3,1)*a(4,3)-a(3,3)*a(4,1))  &
                  +a(2,3)*(a(3,1)*a(4,2)-a(3,2)*a(4,1)))

  b(1,1) =  a(2,2)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))-a(2,3)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))+a(2,4)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))
  b(2,1) = -a(2,1)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))+a(2,3)*(a(3,1)*a(4,4)-a(3,4)*a(4,1))-a(2,4)*(a(3,1)*a(4,3)-a(3,3)*a(4,1))
  b(3,1) =  a(2,1)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))-a(2,2)*(a(3,1)*a(4,4)-a(3,4)*a(4,1))+a(2,4)*(a(3,1)*a(4,2)-a(3,2)*a(4,1))
  b(4,1) = -a(2,1)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))+a(2,2)*(a(3,1)*a(4,3)-a(3,3)*a(4,1))-a(2,3)*(a(3,1)*a(4,2)-a(3,2)*a(4,1))

  b(1,2) = -a(1,2)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))+a(1,3)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))-a(1,4)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))
  b(2,2) =  a(1,1)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))-a(1,3)*(a(3,1)*a(4,4)-a(3,4)*a(4,1))+a(1,4)*(a(3,1)*a(4,3)-a(3,3)*a(4,1))
  b(3,2) = -a(1,1)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))+a(1,2)*(a(3,1)*a(4,4)-a(3,4)*a(4,1))-a(1,4)*(a(3,1)*a(4,2)-a(3,2)*a(4,1))
  b(4,2) =  a(1,1)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))-a(1,2)*(a(3,1)*a(4,3)-a(3,3)*a(4,1))+a(1,3)*(a(3,1)*a(4,2)-a(3,2)*a(4,1))

  b(1,3) =  a(1,2)*(a(2,3)*a(4,4)-a(2,4)*a(4,3))-a(1,3)*(a(2,2)*a(4,4)-a(2,4)*a(4,2))+a(1,4)*(a(2,2)*a(4,3)-a(2,3)*a(4,2))
  b(2,3) = -a(1,1)*(a(2,3)*a(4,4)-a(2,4)*a(4,3))+a(1,3)*(a(2,1)*a(4,4)-a(2,4)*a(4,1))-a(1,4)*(a(2,1)*a(4,3)-a(2,3)*a(4,1))
  b(3,3) =  a(1,1)*(a(2,2)*a(4,4)-a(2,4)*a(4,2))-a(1,2)*(a(2,1)*a(4,4)-a(2,4)*a(4,1))+a(1,4)*(a(2,1)*a(4,2)-a(2,2)*a(4,1))
  b(4,3) = -a(1,1)*(a(2,2)*a(4,3)-a(2,3)*a(4,2))+a(1,2)*(a(2,1)*a(4,3)-a(2,3)*a(4,1))-a(1,3)*(a(2,1)*a(4,2)-a(2,2)*a(4,1))

  b(1,4) = -a(1,2)*(a(2,3)*a(3,4)-a(2,4)*a(3,3))+a(1,3)*(a(2,2)*a(3,4)-a(2,4)*a(3,2))-a(1,4)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))
  b(2,4) =  a(1,1)*(a(2,3)*a(3,4)-a(2,4)*a(3,3))-a(1,3)*(a(2,1)*a(3,4)-a(2,4)*a(3,1))+a(1,4)*(a(2,1)*a(3,3)-a(2,3)*a(3,1))
  b(3,4) = -a(1,1)*(a(2,2)*a(3,4)-a(2,4)*a(3,2))+a(1,2)*(a(2,1)*a(3,4)-a(2,4)*a(3,1))-a(1,4)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
  b(4,4) =  a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))-a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1))+a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))

end subroutine adjugate4

subroutine adjugate5(A,LDA,B,LDB,na,det_l)
  implicit none
  integer*8, intent(in)           :: LDA, LDB
  integer*8, intent(in)           :: na
  double precision, intent(in)    :: A (LDA,na)
  double precision, intent(out)   :: B (LDB,na)
  double precision, intent(inout) :: det_l
  integer :: i,j

 det_l = a(1,1)*(a(2,2)*(a(3,3)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(3,4)*( &
 a(4,3)*a(5,5)-a(4,5)*a(5,3))+a(3,5)*(a(4,3)*a(5,4)-a(4,4)*a(5,3)))- &
 a(2,3)*(a(3,2)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(3,4)*(a(4,2)*a(5,5)- &
 a(4,5)*a(5,2))+a(3,5)*(a(4,2)*a(5,4)-a(4,4)*a(5,2)))+a(2,4)*(a(3,2)*( &
 a(4,3)*a(5,5)-a(4,5)*a(5,3))-a(3,3)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))+ &
 a(3,5)*(a(4,2)*a(5,3)-a(4,3)*a(5,2)))-a(2,5)*(a(3,2)*(a(4,3)*a(5,4)- &
 a(4,4)*a(5,3))-a(3,3)*(a(4,2)*a(5,4)-a(4,4)*a(5,2))+a(3,4)*(a(4,2)* &
 a(5,3)-a(4,3)*a(5,2))))-a(1,2)*(a(2,1)*(a(3,3)*(a(4,4)*a(5,5)-a(4,5)* &
 a(5,4))-a(3,4)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))+a(3,5)*(a(4,3)*a(5,4)- &
 a(4,4)*a(5,3)))-a(2,3)*(a(3,1)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(3,4)*( &
 a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,4)-a(4,4)*a(5,1)))+ &
 a(2,4)*(a(3,1)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))-a(3,3)*(a(4,1)*a(5,5)- &
 a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,3)-a(4,3)*a(5,1)))-a(2,5)*(a(3,1)*( &
 a(4,3)*a(5,4)-a(4,4)*a(5,3))-a(3,3)*(a(4,1)*a(5,4)-a(4,4)*a(5,1))+ &
 a(3,4)*(a(4,1)*a(5,3)-a(4,3)*a(5,1))))+a(1,3)*(a(2,1)*(a(3,2)*(a(4,4)* &
 a(5,5)-a(4,5)*a(5,4))-a(3,4)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))+a(3,5)*( &
 a(4,2)*a(5,4)-a(4,4)*a(5,2)))-a(2,2)*(a(3,1)*(a(4,4)*a(5,5)-a(4,5)* &
 a(5,4))-a(3,4)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,4)- &
 a(4,4)*a(5,1)))+a(2,4)*(a(3,1)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))-a(3,2)*( &
 a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,2)-a(4,2)*a(5,1)))- &
 a(2,5)*(a(3,1)*(a(4,2)*a(5,4)-a(4,4)*a(5,2))-a(3,2)*(a(4,1)*a(5,4)- &
 a(4,4)*a(5,1))+a(3,4)*(a(4,1)*a(5,2)-a(4,2)*a(5,1))))-a(1,4)*(a(2,1)*( &
 a(3,2)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))-a(3,3)*(a(4,2)*a(5,5)-a(4,5)* &
 a(5,2))+a(3,5)*(a(4,2)*a(5,3)-a(4,3)*a(5,2)))-a(2,2)*(a(3,1)*(a(4,3)* &
 a(5,5)-a(4,5)*a(5,3))-a(3,3)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*( &
 a(4,1)*a(5,3)-a(4,3)*a(5,1)))+a(2,3)*(a(3,1)*(a(4,2)*a(5,5)-a(4,5)* &
 a(5,2))-a(3,2)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,2)- &
 a(4,2)*a(5,1)))-a(2,5)*(a(3,1)*(a(4,2)*a(5,3)-a(4,3)*a(5,2))-a(3,2)*( &
 a(4,1)*a(5,3)-a(4,3)*a(5,1))+a(3,3)*(a(4,1)*a(5,2)-a(4,2)*a(5,1))))+ &
 a(1,5)*(a(2,1)*(a(3,2)*(a(4,3)*a(5,4)-a(4,4)*a(5,3))-a(3,3)*(a(4,2)* &
 a(5,4)-a(4,4)*a(5,2))+a(3,4)*(a(4,2)*a(5,3)-a(4,3)*a(5,2)))-a(2,2)*( &
 a(3,1)*(a(4,3)*a(5,4)-a(4,4)*a(5,3))-a(3,3)*(a(4,1)*a(5,4)-a(4,4)* &
 a(5,1))+a(3,4)*(a(4,1)*a(5,3)-a(4,3)*a(5,1)))+a(2,3)*(a(3,1)*(a(4,2)* &
 a(5,4)-a(4,4)*a(5,2))-a(3,2)*(a(4,1)*a(5,4)-a(4,4)*a(5,1))+a(3,4)*( &
 a(4,1)*a(5,2)-a(4,2)*a(5,1)))-a(2,4)*(a(3,1)*(a(4,2)*a(5,3)-a(4,3)* &
 a(5,2))-a(3,2)*(a(4,1)*a(5,3)-a(4,3)*a(5,1))+a(3,3)*(a(4,1)*a(5,2)- &
 a(4,2)*a(5,1))))

 b(1,1) = &
 (a(2,2)*(a(3,3)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(3,4)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))+a(3,5)*(a(4,3)*a(5,4)-a(4,4)*a(5,3)))-a(2,3)* &
 (a(3,2)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(3,4)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))+a(3,5)*(a(4,2)*a(5,4)-a(4,4)*a(5,2)))+a(2,4)* &
 (a(3,2)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))-a(3,3)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))+a(3,5)*(a(4,2)*a(5,3)-a(4,3)*a(5,2)))-a(2,5)* &
 (a(3,2)*(a(4,3)*a(5,4)-a(4,4)*a(5,3))-a(3,3)*(a(4,2)*a(5,4)-a(4,4)*a(5,2))+a(3,4)*(a(4,2)*a(5,3)-a(4,3)*a(5,2))))
 b(2,1) = &
 (-a(2,1)*(a(3,3)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(3,4)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))+a(3,5)*(a(4,3)*a(5,4)-a(4,4)*a(5,3)))+a(2,3)* &
 (a(3,1)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(3,4)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,4)-a(4,4)*a(5,1)))-a(2,4)* &
 (a(3,1)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))-a(3,3)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,3)-a(4,3)*a(5,1)))+a(2,5)* &
 (a(3,1)*(a(4,3)*a(5,4)-a(4,4)*a(5,3))-a(3,3)*(a(4,1)*a(5,4)-a(4,4)*a(5,1))+a(3,4)*(a(4,1)*a(5,3)-a(4,3)*a(5,1))))
 b(3,1) = &
 (a(2,1)*(a(3,2)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(3,4)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))+a(3,5)*(a(4,2)*a(5,4)-a(4,4)*a(5,2)))-a(2,2)* &
 (a(3,1)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(3,4)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,4)-a(4,4)*a(5,1)))+a(2,4)* &
 (a(3,1)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))-a(3,2)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,2)-a(4,2)*a(5,1)))-a(2,5)* &
 (a(3,1)*(a(4,2)*a(5,4)-a(4,4)*a(5,2))-a(3,2)*(a(4,1)*a(5,4)-a(4,4)*a(5,1))+a(3,4)*(a(4,1)*a(5,2)-a(4,2)*a(5,1))))
 b(4,1) = &
 (-a(2,1)*(a(3,2)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))-a(3,3)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))+a(3,5)*(a(4,2)*a(5,3)-a(4,3)*a(5,2)))+a(2,2)* &
 (a(3,1)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))-a(3,3)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,3)-a(4,3)*a(5,1)))-a(2,3)* &
 (a(3,1)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))-a(3,2)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,2)-a(4,2)*a(5,1)))+a(2,5)* &
 (a(3,1)*(a(4,2)*a(5,3)-a(4,3)*a(5,2))-a(3,2)*(a(4,1)*a(5,3)-a(4,3)*a(5,1))+a(3,3)*(a(4,1)*a(5,2)-a(4,2)*a(5,1))))
 b(5,1) = &
 (a(2,1)*(a(3,2)*(a(4,3)*a(5,4)-a(4,4)*a(5,3))-a(3,3)*(a(4,2)*a(5,4)-a(4,4)*a(5,2))+a(3,4)*(a(4,2)*a(5,3)-a(4,3)*a(5,2)))-a(2,2)* &
 (a(3,1)*(a(4,3)*a(5,4)-a(4,4)*a(5,3))-a(3,3)*(a(4,1)*a(5,4)-a(4,4)*a(5,1))+a(3,4)*(a(4,1)*a(5,3)-a(4,3)*a(5,1)))+a(2,3)* &
 (a(3,1)*(a(4,2)*a(5,4)-a(4,4)*a(5,2))-a(3,2)*(a(4,1)*a(5,4)-a(4,4)*a(5,1))+a(3,4)*(a(4,1)*a(5,2)-a(4,2)*a(5,1)))-a(2,4)* &
 (a(3,1)*(a(4,2)*a(5,3)-a(4,3)*a(5,2))-a(3,2)*(a(4,1)*a(5,3)-a(4,3)*a(5,1))+a(3,3)*(a(4,1)*a(5,2)-a(4,2)*a(5,1))))

 b(1,2) = &
 (-a(1,2)*(a(3,3)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(3,4)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))+a(3,5)*(a(4,3)*a(5,4)-a(4,4)*a(5,3)))+a(1,3)* &
 (a(3,2)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(3,4)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))+a(3,5)*(a(4,2)*a(5,4)-a(4,4)*a(5,2)))-a(1,4)* &
 (a(3,2)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))-a(3,3)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))+a(3,5)*(a(4,2)*a(5,3)-a(4,3)*a(5,2)))+a(1,5)* &
 (a(3,2)*(a(4,3)*a(5,4)-a(4,4)*a(5,3))-a(3,3)*(a(4,2)*a(5,4)-a(4,4)*a(5,2))+a(3,4)*(a(4,2)*a(5,3)-a(4,3)*a(5,2))))
 b(2,2) = &
 (a(1,1)*(a(3,3)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(3,4)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))+a(3,5)*(a(4,3)*a(5,4)-a(4,4)*a(5,3)))-a(1,3)* &
 (a(3,1)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(3,4)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,4)-a(4,4)*a(5,1)))+a(1,4)* &
 (a(3,1)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))-a(3,3)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,3)-a(4,3)*a(5,1)))-a(1,5)* &
 (a(3,1)*(a(4,3)*a(5,4)-a(4,4)*a(5,3))-a(3,3)*(a(4,1)*a(5,4)-a(4,4)*a(5,1))+a(3,4)*(a(4,1)*a(5,3)-a(4,3)*a(5,1))))
 b(3,2) = &
 (-a(1,1)*(a(3,2)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(3,4)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))+a(3,5)*(a(4,2)*a(5,4)-a(4,4)*a(5,2)))+a(1,2)* &
 (a(3,1)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(3,4)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,4)-a(4,4)*a(5,1)))-a(1,4)* &
 (a(3,1)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))-a(3,2)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,2)-a(4,2)*a(5,1)))+a(1,5)* &
 (a(3,1)*(a(4,2)*a(5,4)-a(4,4)*a(5,2))-a(3,2)*(a(4,1)*a(5,4)-a(4,4)*a(5,1))+a(3,4)*(a(4,1)*a(5,2)-a(4,2)*a(5,1))))
 b(4,2) = &
 (a(1,1)*(a(3,2)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))-a(3,3)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))+a(3,5)*(a(4,2)*a(5,3)-a(4,3)*a(5,2)))-a(1,2)* &
 (a(3,1)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))-a(3,3)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,3)-a(4,3)*a(5,1)))+a(1,3)* &
 (a(3,1)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))-a(3,2)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(3,5)*(a(4,1)*a(5,2)-a(4,2)*a(5,1)))-a(1,5)* &
 (a(3,1)*(a(4,2)*a(5,3)-a(4,3)*a(5,2))-a(3,2)*(a(4,1)*a(5,3)-a(4,3)*a(5,1))+a(3,3)*(a(4,1)*a(5,2)-a(4,2)*a(5,1))))
 b(5,2) = &
 (-a(1,1)*(a(3,2)*(a(4,3)*a(5,4)-a(4,4)*a(5,3))-a(3,3)*(a(4,2)*a(5,4)-a(4,4)*a(5,2))+a(3,4)*(a(4,2)*a(5,3)-a(4,3)*a(5,2)))+a(1,2)* &
 (a(3,1)*(a(4,3)*a(5,4)-a(4,4)*a(5,3))-a(3,3)*(a(4,1)*a(5,4)-a(4,4)*a(5,1))+a(3,4)*(a(4,1)*a(5,3)-a(4,3)*a(5,1)))-a(1,3)* &
 (a(3,1)*(a(4,2)*a(5,4)-a(4,4)*a(5,2))-a(3,2)*(a(4,1)*a(5,4)-a(4,4)*a(5,1))+a(3,4)*(a(4,1)*a(5,2)-a(4,2)*a(5,1)))+a(1,4)* &
 (a(3,1)*(a(4,2)*a(5,3)-a(4,3)*a(5,2))-a(3,2)*(a(4,1)*a(5,3)-a(4,3)*a(5,1))+a(3,3)*(a(4,1)*a(5,2)-a(4,2)*a(5,1))))

 b(1,3) = &
 (a(1,2)*(a(2,3)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(2,4)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))+a(2,5)*(a(4,3)*a(5,4)-a(4,4)*a(5,3)))-a(1,3)* &
 (a(2,2)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(2,4)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))+a(2,5)*(a(4,2)*a(5,4)-a(4,4)*a(5,2)))+a(1,4)* &
 (a(2,2)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))-a(2,3)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))+a(2,5)*(a(4,2)*a(5,3)-a(4,3)*a(5,2)))-a(1,5)* &
 (a(2,2)*(a(4,3)*a(5,4)-a(4,4)*a(5,3))-a(2,3)*(a(4,2)*a(5,4)-a(4,4)*a(5,2))+a(2,4)*(a(4,2)*a(5,3)-a(4,3)*a(5,2))))
 b(2,3) = &
 (-a(1,1)*(a(2,3)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(2,4)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))+a(2,5)*(a(4,3)*a(5,4)-a(4,4)*a(5,3)))+a(1,3)* &
 (a(2,1)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(2,4)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(2,5)*(a(4,1)*a(5,4)-a(4,4)*a(5,1)))-a(1,4)* &
 (a(2,1)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))-a(2,3)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(2,5)*(a(4,1)*a(5,3)-a(4,3)*a(5,1)))+a(1,5)* &
 (a(2,1)*(a(4,3)*a(5,4)-a(4,4)*a(5,3))-a(2,3)*(a(4,1)*a(5,4)-a(4,4)*a(5,1))+a(2,4)*(a(4,1)*a(5,3)-a(4,3)*a(5,1))))
 b(3,3) = &
 (a(1,1)*(a(2,2)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(2,4)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))+a(2,5)*(a(4,2)*a(5,4)-a(4,4)*a(5,2)))-a(1,2)* &
 (a(2,1)*(a(4,4)*a(5,5)-a(4,5)*a(5,4))-a(2,4)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(2,5)*(a(4,1)*a(5,4)-a(4,4)*a(5,1)))+a(1,4)* &
 (a(2,1)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))-a(2,2)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(2,5)*(a(4,1)*a(5,2)-a(4,2)*a(5,1)))-a(1,5)* &
 (a(2,1)*(a(4,2)*a(5,4)-a(4,4)*a(5,2))-a(2,2)*(a(4,1)*a(5,4)-a(4,4)*a(5,1))+a(2,4)*(a(4,1)*a(5,2)-a(4,2)*a(5,1))))
 b(4,3) = &
 (-a(1,1)*(a(2,2)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))-a(2,3)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))+a(2,5)*(a(4,2)*a(5,3)-a(4,3)*a(5,2)))+a(1,2)* &
 (a(2,1)*(a(4,3)*a(5,5)-a(4,5)*a(5,3))-a(2,3)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(2,5)*(a(4,1)*a(5,3)-a(4,3)*a(5,1)))-a(1,3)* &
 (a(2,1)*(a(4,2)*a(5,5)-a(4,5)*a(5,2))-a(2,2)*(a(4,1)*a(5,5)-a(4,5)*a(5,1))+a(2,5)*(a(4,1)*a(5,2)-a(4,2)*a(5,1)))+a(1,5)* &
 (a(2,1)*(a(4,2)*a(5,3)-a(4,3)*a(5,2))-a(2,2)*(a(4,1)*a(5,3)-a(4,3)*a(5,1))+a(2,3)*(a(4,1)*a(5,2)-a(4,2)*a(5,1))))
 b(5,3) = &
 (a(1,1)*(a(2,2)*(a(4,3)*a(5,4)-a(4,4)*a(5,3))-a(2,3)*(a(4,2)*a(5,4)-a(4,4)*a(5,2))+a(2,4)*(a(4,2)*a(5,3)-a(4,3)*a(5,2)))-a(1,2)* &
 (a(2,1)*(a(4,3)*a(5,4)-a(4,4)*a(5,3))-a(2,3)*(a(4,1)*a(5,4)-a(4,4)*a(5,1))+a(2,4)*(a(4,1)*a(5,3)-a(4,3)*a(5,1)))+a(1,3)* &
 (a(2,1)*(a(4,2)*a(5,4)-a(4,4)*a(5,2))-a(2,2)*(a(4,1)*a(5,4)-a(4,4)*a(5,1))+a(2,4)*(a(4,1)*a(5,2)-a(4,2)*a(5,1)))-a(1,4)* &
 (a(2,1)*(a(4,2)*a(5,3)-a(4,3)*a(5,2))-a(2,2)*(a(4,1)*a(5,3)-a(4,3)*a(5,1))+a(2,3)*(a(4,1)*a(5,2)-a(4,2)*a(5,1))))

 b(1,4) = &
 (-a(1,2)*(a(2,3)*(a(3,4)*a(5,5)-a(3,5)*a(5,4))-a(2,4)*(a(3,3)*a(5,5)-a(3,5)*a(5,3))+a(2,5)*(a(3,3)*a(5,4)-a(3,4)*a(5,3)))+a(1,3)* &
 (a(2,2)*(a(3,4)*a(5,5)-a(3,5)*a(5,4))-a(2,4)*(a(3,2)*a(5,5)-a(3,5)*a(5,2))+a(2,5)*(a(3,2)*a(5,4)-a(3,4)*a(5,2)))-a(1,4)* &
 (a(2,2)*(a(3,3)*a(5,5)-a(3,5)*a(5,3))-a(2,3)*(a(3,2)*a(5,5)-a(3,5)*a(5,2))+a(2,5)*(a(3,2)*a(5,3)-a(3,3)*a(5,2)))+a(1,5)* &
 (a(2,2)*(a(3,3)*a(5,4)-a(3,4)*a(5,3))-a(2,3)*(a(3,2)*a(5,4)-a(3,4)*a(5,2))+a(2,4)*(a(3,2)*a(5,3)-a(3,3)*a(5,2))))
 b(2,4) = &
 (a(1,1)*(a(2,3)*(a(3,4)*a(5,5)-a(3,5)*a(5,4))-a(2,4)*(a(3,3)*a(5,5)-a(3,5)*a(5,3))+a(2,5)*(a(3,3)*a(5,4)-a(3,4)*a(5,3)))-a(1,3)* &
 (a(2,1)*(a(3,4)*a(5,5)-a(3,5)*a(5,4))-a(2,4)*(a(3,1)*a(5,5)-a(3,5)*a(5,1))+a(2,5)*(a(3,1)*a(5,4)-a(3,4)*a(5,1)))+a(1,4)* &
 (a(2,1)*(a(3,3)*a(5,5)-a(3,5)*a(5,3))-a(2,3)*(a(3,1)*a(5,5)-a(3,5)*a(5,1))+a(2,5)*(a(3,1)*a(5,3)-a(3,3)*a(5,1)))-a(1,5)* &
 (a(2,1)*(a(3,3)*a(5,4)-a(3,4)*a(5,3))-a(2,3)*(a(3,1)*a(5,4)-a(3,4)*a(5,1))+a(2,4)*(a(3,1)*a(5,3)-a(3,3)*a(5,1))))
 b(3,4) = &
 (-a(1,1)*(a(2,2)*(a(3,4)*a(5,5)-a(3,5)*a(5,4))-a(2,4)*(a(3,2)*a(5,5)-a(3,5)*a(5,2))+a(2,5)*(a(3,2)*a(5,4)-a(3,4)*a(5,2)))+a(1,2)* &
 (a(2,1)*(a(3,4)*a(5,5)-a(3,5)*a(5,4))-a(2,4)*(a(3,1)*a(5,5)-a(3,5)*a(5,1))+a(2,5)*(a(3,1)*a(5,4)-a(3,4)*a(5,1)))-a(1,4)* &
 (a(2,1)*(a(3,2)*a(5,5)-a(3,5)*a(5,2))-a(2,2)*(a(3,1)*a(5,5)-a(3,5)*a(5,1))+a(2,5)*(a(3,1)*a(5,2)-a(3,2)*a(5,1)))+a(1,5)* &
 (a(2,1)*(a(3,2)*a(5,4)-a(3,4)*a(5,2))-a(2,2)*(a(3,1)*a(5,4)-a(3,4)*a(5,1))+a(2,4)*(a(3,1)*a(5,2)-a(3,2)*a(5,1))))
 b(4,4) = &
 (a(1,1)*(a(2,2)*(a(3,3)*a(5,5)-a(3,5)*a(5,3))-a(2,3)*(a(3,2)*a(5,5)-a(3,5)*a(5,2))+a(2,5)*(a(3,2)*a(5,3)-a(3,3)*a(5,2)))-a(1,2)* &
 (a(2,1)*(a(3,3)*a(5,5)-a(3,5)*a(5,3))-a(2,3)*(a(3,1)*a(5,5)-a(3,5)*a(5,1))+a(2,5)*(a(3,1)*a(5,3)-a(3,3)*a(5,1)))+a(1,3)* &
 (a(2,1)*(a(3,2)*a(5,5)-a(3,5)*a(5,2))-a(2,2)*(a(3,1)*a(5,5)-a(3,5)*a(5,1))+a(2,5)*(a(3,1)*a(5,2)-a(3,2)*a(5,1)))-a(1,5)* &
 (a(2,1)*(a(3,2)*a(5,3)-a(3,3)*a(5,2))-a(2,2)*(a(3,1)*a(5,3)-a(3,3)*a(5,1))+a(2,3)*(a(3,1)*a(5,2)-a(3,2)*a(5,1))))
 b(5,4) = &
 (-a(1,1)*(a(2,2)*(a(3,3)*a(5,4)-a(3,4)*a(5,3))-a(2,3)*(a(3,2)*a(5,4)-a(3,4)*a(5,2))+a(2,4)*(a(3,2)*a(5,3)-a(3,3)*a(5,2)))+a(1,2)* &
 (a(2,1)*(a(3,3)*a(5,4)-a(3,4)*a(5,3))-a(2,3)*(a(3,1)*a(5,4)-a(3,4)*a(5,1))+a(2,4)*(a(3,1)*a(5,3)-a(3,3)*a(5,1)))-a(1,3)* &
 (a(2,1)*(a(3,2)*a(5,4)-a(3,4)*a(5,2))-a(2,2)*(a(3,1)*a(5,4)-a(3,4)*a(5,1))+a(2,4)*(a(3,1)*a(5,2)-a(3,2)*a(5,1)))+a(1,4)* &
 (a(2,1)*(a(3,2)*a(5,3)-a(3,3)*a(5,2))-a(2,2)*(a(3,1)*a(5,3)-a(3,3)*a(5,1))+a(2,3)*(a(3,1)*a(5,2)-a(3,2)*a(5,1))))

 b(1,5) = &
 (a(1,2)*(a(2,3)*(a(3,4)*a(4,5)-a(3,5)*a(4,4))-a(2,4)*(a(3,3)*a(4,5)-a(3,5)*a(4,3))+a(2,5)*(a(3,3)*a(4,4)-a(3,4)*a(4,3)))-a(1,3)* &
 (a(2,2)*(a(3,4)*a(4,5)-a(3,5)*a(4,4))-a(2,4)*(a(3,2)*a(4,5)-a(3,5)*a(4,2))+a(2,5)*(a(3,2)*a(4,4)-a(3,4)*a(4,2)))+a(1,4)* &
 (a(2,2)*(a(3,3)*a(4,5)-a(3,5)*a(4,3))-a(2,3)*(a(3,2)*a(4,5)-a(3,5)*a(4,2))+a(2,5)*(a(3,2)*a(4,3)-a(3,3)*a(4,2)))-a(1,5)* &
 (a(2,2)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))-a(2,3)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))+a(2,4)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))))
 b(2,5) = &
 (-a(1,1)*(a(2,3)*(a(3,4)*a(4,5)-a(3,5)*a(4,4))-a(2,4)*(a(3,3)*a(4,5)-a(3,5)*a(4,3))+a(2,5)*(a(3,3)*a(4,4)-a(3,4)*a(4,3)))+a(1,3)* &
 (a(2,1)*(a(3,4)*a(4,5)-a(3,5)*a(4,4))-a(2,4)*(a(3,1)*a(4,5)-a(3,5)*a(4,1))+a(2,5)*(a(3,1)*a(4,4)-a(3,4)*a(4,1)))-a(1,4)* &
 (a(2,1)*(a(3,3)*a(4,5)-a(3,5)*a(4,3))-a(2,3)*(a(3,1)*a(4,5)-a(3,5)*a(4,1))+a(2,5)*(a(3,1)*a(4,3)-a(3,3)*a(4,1)))+a(1,5)* &
 (a(2,1)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))-a(2,3)*(a(3,1)*a(4,4)-a(3,4)*a(4,1))+a(2,4)*(a(3,1)*a(4,3)-a(3,3)*a(4,1))))
 b(3,5) = &
 (a(1,1)*(a(2,2)*(a(3,4)*a(4,5)-a(3,5)*a(4,4))-a(2,4)*(a(3,2)*a(4,5)-a(3,5)*a(4,2))+a(2,5)*(a(3,2)*a(4,4)-a(3,4)*a(4,2)))-a(1,2)* &
 (a(2,1)*(a(3,4)*a(4,5)-a(3,5)*a(4,4))-a(2,4)*(a(3,1)*a(4,5)-a(3,5)*a(4,1))+a(2,5)*(a(3,1)*a(4,4)-a(3,4)*a(4,1)))+a(1,4)* &
 (a(2,1)*(a(3,2)*a(4,5)-a(3,5)*a(4,2))-a(2,2)*(a(3,1)*a(4,5)-a(3,5)*a(4,1))+a(2,5)*(a(3,1)*a(4,2)-a(3,2)*a(4,1)))-a(1,5)* &
 (a(2,1)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))-a(2,2)*(a(3,1)*a(4,4)-a(3,4)*a(4,1))+a(2,4)*(a(3,1)*a(4,2)-a(3,2)*a(4,1))))
 b(4,5) = &
 (-a(1,1)*(a(2,2)*(a(3,3)*a(4,5)-a(3,5)*a(4,3))-a(2,3)*(a(3,2)*a(4,5)-a(3,5)*a(4,2))+a(2,5)*(a(3,2)*a(4,3)-a(3,3)*a(4,2)))+a(1,2)* &
 (a(2,1)*(a(3,3)*a(4,5)-a(3,5)*a(4,3))-a(2,3)*(a(3,1)*a(4,5)-a(3,5)*a(4,1))+a(2,5)*(a(3,1)*a(4,3)-a(3,3)*a(4,1)))-a(1,3)* &
 (a(2,1)*(a(3,2)*a(4,5)-a(3,5)*a(4,2))-a(2,2)*(a(3,1)*a(4,5)-a(3,5)*a(4,1))+a(2,5)*(a(3,1)*a(4,2)-a(3,2)*a(4,1)))+a(1,5)* &
 (a(2,1)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))-a(2,2)*(a(3,1)*a(4,3)-a(3,3)*a(4,1))+a(2,3)*(a(3,1)*a(4,2)-a(3,2)*a(4,1))))
 b(5,5) = &
 (a(1,1)*(a(2,2)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))-a(2,3)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))+a(2,4)*(a(3,2)*a(4,3)-a(3,3)*a(4,2)))-a(1,2)* &
 (a(2,1)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))-a(2,3)*(a(3,1)*a(4,4)-a(3,4)*a(4,1))+a(2,4)*(a(3,1)*a(4,3)-a(3,3)*a(4,1)))+a(1,3)* &
 (a(2,1)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))-a(2,2)*(a(3,1)*a(4,4)-a(3,4)*a(4,1))+a(2,4)*(a(3,1)*a(4,2)-a(3,2)*a(4,1)))-a(1,4)* &
 (a(2,1)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))-a(2,2)*(a(3,1)*a(4,3)-a(3,3)*a(4,1))+a(2,3)*(a(3,1)*a(4,2)-a(3,2)*a(4,1))))

end

subroutine adjugate_general(context, na, A, LDA, B, LDB, det_l)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in) :: context
  integer*8, intent(in)              :: LDA
  integer*8, intent(in)              :: LDB
  integer*8, intent(in)              :: na
  double precision, intent(in)       :: A (LDA,na)
  double precision, intent(out)      :: B (LDB,na)
  double precision, intent(inout)    :: det_l

  double precision :: work(LDA*max(na,64))
  integer          :: inf
  integer          :: ipiv(LDA)
  integer          :: lwork
  integer(8)       :: i, j

B(1:na,1:na) = A(1:na,1:na)

call dgetrf(na, na, B, LDB, ipiv, inf )

det_l = 1.d0
j=0_8
do i=1,na
 j = j+min(abs(ipiv(i)-i),1)
 det_l = det_l*B(i,i)
enddo

if (iand(j,1_8) /= 0_8)  then
  det_l = -det_l
endif

lwork = SIZE(work)
call dgetri(na, B, LDB, ipiv, work, lwork, inf )

B(1:na,1:na) = B(1:na,1:na)*det_l

end subroutine adjugate_general

function qmckl_adjugate_safe(context, &
     na, A, size_A, LDA, B, size_B, LDB, det_l) &
     result(info) bind(C)
  use qmckl_constants
  use qmckl, only: qmckl_adjugate

  implicit none

  integer (qmckl_context), intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: na
  integer (c_int64_t) , intent(in)  , value :: size_A
  integer (c_int64_t) , intent(in)  , value :: lda
  integer (c_int64_t) , intent(in)  , value :: size_B
  integer (c_int64_t) , intent(in)  , value :: ldb
  real    (c_double ) , intent(inout)        :: det_l
  real    (c_double ) , intent(in)          :: A(lda,*)
  real    (c_double ) , intent(out)         :: B(ldb,*)

  integer(qmckl_exit_code) :: info

  info = QMCKL_SUCCESS

  if (size_A < na) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  if (size_B <= 0_8) then
     info = QMCKL_INVALID_ARG_7
     return
  endif

  info = qmckl_adjugate(context, na, A, LDA, B, LDB, det_l)

  if (info == QMCKL_INVALID_ARG_4) then
     info = QMCKL_INVALID_ARG_5
     return
  endif

  if (info == QMCKL_INVALID_ARG_6) then
     info = QMCKL_INVALID_ARG_8
     return
  endif

end function qmckl_adjugate_safe

function qmckl_dgemv(context, Trans, &
     m, n, alpha, A, LDA, X, incx, beta, Y, incy) &
     bind(C) result(info)
  use qmckl_constants
  implicit none
  integer (qmckl_context), intent(in)  , value :: context
  character(c_char  ) , intent(in)  , value :: Trans
  integer (c_int64_t) , intent(in)  , value :: m
  integer (c_int64_t) , intent(in)  , value :: n
  real    (c_double ) , intent(in)  , value :: alpha
  integer (c_int64_t) , intent(in)  , value :: lda
  integer (c_int64_t) , intent(in)  , value :: incx
  integer (c_int64_t) , intent(in)  , value :: incy
  real    (c_double ) , intent(in)  , value :: beta
  real    (c_double ) , intent(in)          :: A(lda,*)
  real    (c_double ) , intent(in)          :: X(*)
  real    (c_double ) , intent(out)         :: Y(*)
  integer(qmckl_exit_code) :: info

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (m <= 0_8) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  if (n <= 0_8) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  if (LDA <= 0) then
     info = QMCKL_INVALID_ARG_7
     return
  endif

  if (incx <= 0) then
     info = QMCKL_INVALID_ARG_9
     return
  endif

  if (incy <= 0) then
     info = QMCKL_INVALID_ARG_12
     return
  endif

  call dgemv(trans, int(m,4), int(n,4), &
       alpha, A, int(LDA,4), X, int(incx,4), beta, Y, int(incy,4))


end function qmckl_dgemv

function qmckl_dger(context, &
     m, n, alpha, X, incx, Y, incy, A, LDA) &
     bind(C) result(info)
  use qmckl_constants
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
  integer(qmckl_exit_code) :: info

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (m <= 0_8) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (n <= 0_8) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  if (LDA <= 0) then
     info = QMCKL_INVALID_ARG_10
     return
  endif

  if (incx <= 0) then
     info = QMCKL_INVALID_ARG_6
     return
  endif

  if (incy <= 0) then
     info = QMCKL_INVALID_ARG_8
     return
  endif

  call dger(int(m,4), int(n,4), alpha, X, int(incx,4), Y, int(incy,4), A, int(LDA,4))


end function qmckl_dger
