integer(qmckl_exit_code) function test_qmckl_dgemm(context) bind(C)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in), value :: context

  double precision, allocatable :: A(:,:), B(:,:), C(:,:), D(:,:)
  integer*8                     :: m, n, k, LDA, LDB, LDC
  integer*8                     :: i,j,l
  character                     :: TransA, TransB
  double precision              :: x, alpha, beta

  TransA = 'N'
  TransB = 'N'
  m = 1_8
  k = 4_8
  n = 6_8
  LDA = m
  LDB = k
  LDC = m

  allocate( A(LDA,k), B(LDB,n) , C(LDC,n), D(LDC,n))

  A = 0.d0
  B = 0.d0
  C = 0.d0
  D = 0.d0
  alpha = 1.0d0
  beta  = 0.0d0
  do j=1,k
     do i=1,m
        A(i,j) = -10.d0 + dble(i+j)
     end do
  end do

  do j=1,n
     do i=1,k
        B(i,j) = -10.d0 + dble(i+j)
     end do
  end do

  test_qmckl_dgemm = qmckl_dgemm(context, TransA, TransB, m, n, k, &
       alpha, A, LDA, B, LDB, beta, C, LDC)

  if (test_qmckl_dgemm /= QMCKL_SUCCESS) return

  test_qmckl_dgemm = QMCKL_FAILURE

  x = 0.d0
  do j=1,n
     do i=1,m
        do l=1,k
           D(i,j) = D(i,j) + A(i,l)*B(l,j)
        end do
        x = x + (D(i,j) - C(i,j))**2
     end do
  end do

  if (dabs(x) <= 1.d-12) then
     test_qmckl_dgemm = QMCKL_SUCCESS
  endif

  deallocate(A,B,C,D)

end function test_qmckl_dgemm

integer(qmckl_exit_code) function test_qmckl_adjugate(context) bind(C)
  use qmckl_constants
  implicit none
  integer(qmckl_context), intent(in), value :: context

  double precision, allocatable :: A(:,:), B(:,:)
  integer*8                     :: m, n, k, LDA, LDB
  integer*8                     :: i,j,l
  double precision              :: x, det_l, det_l_ref

  LDA = 6_8
  LDB = 6_8

  allocate( A(LDA,6), B(LDB,6))

  A = 0.1d0
  A(1,1) = 1.0d0;
  A(2,2) = 2.0d0;
  A(3,3) = 3.0d0;
  A(4,4) = 4.0d0;
  A(5,5) = 5.0d0;
  A(6,6) = 6.0d0;

  test_qmckl_adjugate = QMCKL_SUCCESS

deallocate(A,B)

end function test_qmckl_adjugate
