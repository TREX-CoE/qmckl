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
  use qmckl
  implicit none
  integer(qmckl_context), intent(in), value :: context

  double precision, allocatable :: A(:,:), B(:,:)
  integer*8                     :: m, n, k, LDA, LDB
  integer*8                     :: i,j,l
  double precision              :: x, det_l, det_l_ref

  LDA = 6_8
  LDB = 6_8

  allocate( A(LDA,6), B(LDB,6))

  do i=1,6
     do j=1,i-1
        A(i,j) = 0.1d0*i-j
     enddo
     A(i,i) = 1.d0+i
     do j=i+1,6
        A(i,j) = 0.1d0*i+j
     enddo
  enddo
  
  test_qmckl_adjugate = QMCKL_SUCCESS

! N =  1
  test_qmckl_adjugate = qmckl_adjugate(context, 1_8, A, LDA, B, LDB, det_l)
  if (test_qmckl_adjugate /= QMCKL_SUCCESS) return
  if (dabs((det_l - (2.0d0))/det_l) > 1.d-13) then
      print *, 'N = 1: det = ', det_l, 2.0d0
      test_qmckl_adjugate = 1
      return
  end if

  x = 0.d0
  x = x + (B(1,1) - (1.0) )**2
  if (dabs(x / det_l) > 1.d-10) then
      print *, 'N = 1: x = ', x
      print *, '         det = ', det_l
      test_qmckl_adjugate = 1
      return
  end if

! N =  2
  test_qmckl_adjugate = qmckl_adjugate(context, 2_8, A, LDA, B, LDB, det_l)
  if (test_qmckl_adjugate /= QMCKL_SUCCESS) return
  if (dabs((det_l - (7.6800000000000015d0))/det_l) > 1.d-13) then
      print *, 'N = 2: det = ', det_l, 7.6800000000000015d0
      test_qmckl_adjugate = 2
      return
  end if

  x = 0.d0
  x = x + (B(1,1) - (3.0000000000000004) )**2
  x = x + (B(1,2) - (-2.100000000000001) )**2
  x = x + (B(2,1) - (0.8000000000000003) )**2
  x = x + (B(2,2) - (2.0000000000000004) )**2
  if (dabs(x / det_l) > 1.d-10) then
      print *, 'N = 2: x = ', x
      print *, '         det = ', det_l
      test_qmckl_adjugate = 2
      return
  end if

! N =  3
  test_qmckl_adjugate = qmckl_adjugate(context, 3_8, A, LDA, B, LDB, det_l)
  if (test_qmckl_adjugate /= QMCKL_SUCCESS) return
  if (dabs((det_l - (47.62200000000001d0))/det_l) > 1.d-13) then
      print *, 'N = 3: det = ', det_l, 47.62200000000001d0
      test_qmckl_adjugate = 3
      return
  end if

  x = 0.d0
  x = x + (B(1,1) - (17.44) )**2
  x = x + (B(1,2) - (-13.670000000000003) )**2
  x = x + (B(1,3) - (-2.579999999999999) )**2
  x = x + (B(2,1) - (0.96) )**2
  x = x + (B(2,2) - (10.170000000000002) )**2
  x = x + (B(2,3) - (-8.880000000000003) )**2
  x = x + (B(3,1) - (3.4600000000000004) )**2
  x = x + (B(3,2) - (1.9300000000000006) )**2
  x = x + (B(3,3) - (7.6800000000000015) )**2
  if (dabs(x / det_l) > 1.d-10) then
      print *, 'N = 3: x = ', x
      print *, '         det = ', det_l
      test_qmckl_adjugate = 3
      return
  end if

! N =  4
  test_qmckl_adjugate = qmckl_adjugate(context, 4_8, A, LDA, B, LDB, det_l)
  if (test_qmckl_adjugate /= QMCKL_SUCCESS) return
  if (dabs((det_l - (397.2747999999998d0))/det_l) > 1.d-13) then
      print *, 'N = 4: det = ', det_l, 397.2747999999998d0
      test_qmckl_adjugate = 4
      return
  end if

  x = 0.d0
  x = x + (B(1,1) - (144.16799999999992) )**2
  x = x + (B(1,2) - (-114.86199999999995) )**2
  x = x + (B(1,3) - (-21.787999999999975) )**2
  x = x + (B(1,4) - (-2.996000000000003) )**2
  x = x + (B(2,1) - (4.275999999999998) )**2
  x = x + (B(2,2) - (82.51399999999997) )**2
  x = x + (B(2,3) - (-74.82799999999996) )**2
  x = x + (B(2,4) - (-8.465999999999992) )**2
  x = x + (B(3,1) - (4.475999999999997) )**2
  x = x + (B(3,2) - (0.8980000000000001) )**2
  x = x + (B(3,3) - (59.175999999999966) )**2
  x = x + (B(3,4) - (-55.31599999999997) )**2
  x = x + (B(4,1) - (20.99599999999999) )**2
  x = x + (B(4,2) - (13.087999999999996) )**2
  x = x + (B(4,3) - (4.211999999999997) )**2
  x = x + (B(4,4) - (47.62199999999998) )**2
  if (dabs(x / det_l) > 1.d-10) then
      print *, 'N = 4: x = ', x
      print *, '         det = ', det_l
      test_qmckl_adjugate = 4
      return
  end if

! N =  5
  test_qmckl_adjugate = qmckl_adjugate(context, 5_8, A, LDA, B, LDB, det_l)
  if (test_qmckl_adjugate /= QMCKL_SUCCESS) return
  if (dabs((det_l - (4096.899599999997d0))/det_l) > 1.d-13) then
      print *, 'N = 5: det = ', det_l, 4096.899599999997d0
      test_qmckl_adjugate = 5
      return
  end if

  x = 0.d0
  x = x + (B(1,1) - (1484.1379999999992) )**2
  x = x + (B(1,2) - (-1186.3349999999991) )**2
  x = x + (B(1,3) - (-225.31799999999956) )**2
  x = x + (B(1,4) - (-31.121999999999847) )**2
  x = x + (B(1,5) - (-6.319600000000006) )**2
  x = x + (B(2,1) - (40.574) )**2
  x = x + (B(2,2) - (848.4569999999994) )**2
  x = x + (B(2,3) - (-772.5179999999996) )**2
  x = x + (B(2,4) - (-87.61199999999984) )**2
  x = x + (B(2,5) - (-8.575599999999962) )**2
  x = x + (B(3,1) - (28.733999999999984) )**2
  x = x + (B(3,2) - (-2.955000000000003) )**2
  x = x + (B(3,3) - (606.0299999999996) )**2
  x = x + (B(3,4) - (-571.9619999999996) )**2
  x = x + (B(3,5) - (-42.42359999999996) )**2
  x = x + (B(4,1) - (29.79399999999999) )**2
  x = x + (B(4,2) - (4.065000000000015) )**2
  x = x + (B(4,3) - (-1.8180000000000054) )**2
  x = x + (B(4,4) - (474.8639999999997) )**2
  x = x + (B(4,5) - (-454.61959999999965) )**2
  x = x + (B(5,1) - (163.17399999999986) )**2
  x = x + (B(5,2) - (114.39299999999992) )**2
  x = x + (B(5,3) - (39.54599999999997) )**2
  x = x + (B(5,4) - (14.190000000000001) )**2
  x = x + (B(5,5) - (397.2747999999997) )**2
  if (dabs(x / det_l) > 1.d-10) then
      print *, 'N = 5: x = ', x
      print *, '         det = ', det_l
      test_qmckl_adjugate = 5
      return
  end if

! N =  6
  test_qmckl_adjugate = qmckl_adjugate(context, 6_8, A, LDA, B, LDB, det_l)
  if (test_qmckl_adjugate /= QMCKL_SUCCESS) return
  if (dabs((det_l - (50129.481119999975d0))/det_l) > 1.d-13) then
      print *, 'N = 6: det = ', det_l, 50129.481119999975d0
      test_qmckl_adjugate = 6
      return
  end if

  x = 0.d0
  x = x + (B(1,1) - (18145.49839999999) )**2
  x = x + (B(1,2) - (-14527.349999999993) )**2
  x = x + (B(1,3) - (-2761.183199999996) )**2
  x = x + (B(1,4) - (-382.3847999999964) )**2
  x = x + (B(1,5) - (-78.1278399999993) )**2
  x = x + (B(1,6) - (-38.203200000000116) )**2
  x = x + (B(2,1) - (487.2176) )**2
  x = x + (B(2,2) - (10374.332399999996) )**2
  x = x + (B(2,3) - (-9455.203199999998) )**2
  x = x + (B(2,4) - (-1073.0327999999997) )**2
  x = x + (B(2,5) - (-105.44704000000002) )**2
  x = x + (B(2,6) - (-24.613199999999587) )**2
  x = x + (B(3,1) - (337.07760000000064) )**2
  x = x + (B(3,2) - (-47.693999999999484) )**2
  x = x + (B(3,3) - (7411.106399999997) )**2
  x = x + (B(3,4) - (-7000.0967999999975) )**2
  x = x + (B(3,5) - (-519.9038399999989) )**2
  x = x + (B(3,6) - (-38.63519999999899) )**2
  x = x + (B(4,1) - (259.5495999999998) )**2
  x = x + (B(4,2) - (-33.74999999999984) )**2
  x = x + (B(4,3) - (-52.99920000000007) )**2
  x = x + (B(4,4) - (5798.870399999997) )**2
  x = x + (B(4,5) - (-5568.571039999997) )**2
  x = x + (B(4,6) - (-279.5952000000003) )**2
  x = x + (B(5,1) - (258.8727999999996) )**2
  x = x + (B(5,2) - (18.099599999999644) )**2
  x = x + (B(5,3) - (-25.048799999999975) )**2
  x = x + (B(5,4) - (-17.411999999999882) )**2
  x = x + (B(5,5) - (4763.966559999997) )**2
  x = x + (B(5,6) - (-4626.839999999998) )**2
  x = x + (B(6,1) - (1538.685599999999) )**2
  x = x + (B(6,2) - (1223.3639999999996) )**2
  x = x + (B(6,3) - (450.6407999999997) )**2
  x = x + (B(6,4) - (169.15919999999994) )**2
  x = x + (B(6,5) - (85.95215999999988) )**2
  x = x + (B(6,6) - (4096.899599999998) )**2
  if (dabs(x / det_l) > 1.d-10) then
      print *, 'N = 6: x = ', x
      print *, '         det = ', det_l
      test_qmckl_adjugate = 6
      return
  end if

deallocate(A,B)

end function test_qmckl_adjugate
