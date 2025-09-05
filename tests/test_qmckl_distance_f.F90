integer(qmckl_exit_code) function test_qmckl_distance_sq(context) bind(C)

  use qmckl
  use qmckl_verificarlo_f

  implicit none

  integer(qmckl_context), intent(in), value :: context
  logical(C_BOOL) :: vfc_err

  double precision, allocatable :: A(:,:), B(:,:), C(:,:)
  integer*8                     :: m, n, LDA, LDB, LDC
  double precision              :: x
  integer*8                     :: i,j

  m = 5
  n = 6
  LDA = m
  LDB = n
  LDC = 5

  allocate( A(LDA,m), B(LDB,n), C(LDC,n) )
  do j=1,m
     do i=1,m
        A(i,j) = -10.d0 + dble(i+j)
     end do
  end do
  do j=1,n
     do i=1,n
        B(i,j) = -1.d0 + dble(i*j)
     end do
  end do

  test_qmckl_distance_sq = &
       qmckl_distance_sq(context, 'X', 't', m, n, A, LDA, B, LDB, C, LDC)

  vfc_err = qmckl_probe("distance", "distance_sq_Xt_2_2", C(2,2))

  if (test_qmckl_distance_sq == 0) return


  test_qmckl_distance_sq = &
       qmckl_distance_sq(context, 't', 'X', m, n, A, LDA, B, LDB, C, LDC)

  vfc_err = qmckl_probe("distance", "distance_sq_tX_2_2", C(2,2))

  if (test_qmckl_distance_sq == 0) return


  test_qmckl_distance_sq = &
       qmckl_distance_sq(context, 'T', 't', m, n, A, LDA, B, LDB, C, LDC)

  vfc_err = qmckl_probe_check("distance", "distance_sq_Tt_2_2", C(2,2), 0.d0, 1.d-14)

  if (test_qmckl_distance_sq == 0) return


  test_qmckl_distance_sq = QMCKL_FAILURE

  do j=1,n
     do i=1,m
        x =  (A(i,1)-B(j,1))**2 + &
             (A(i,2)-B(j,2))**2 + &
             (A(i,3)-B(j,3))**2
#ifndef VFC_CI
        if ( dabs(1.d0 - C(i,j)/x) > 1.d-14) return
#endif
     end do
  end do

  test_qmckl_distance_sq = &
       qmckl_distance_sq(context, 'n', 'T', m, n, A, LDA, B, LDB, C, LDC)

  vfc_err = qmckl_probe_check("distance", "distance_sq_nT_2_2", C(2,2), 0.d0, 1.d-14)


  test_qmckl_distance_sq = QMCKL_FAILURE

  do j=1,n
     do i=1,m
        x =  (A(1,i)-B(j,1))**2 + &
             (A(2,i)-B(j,2))**2 + &
             (A(3,i)-B(j,3))**2
#ifndef VFC_CI
        if ( dabs(1.d0 - C(i,j)/x) > 1.d-14) return
#endif
     end do
  end do

  test_qmckl_distance_sq = &
       qmckl_distance_sq(context, 'T', 'n', m, n, A, LDA, B, LDB, C, LDC)

  vfc_err =  qmckl_probe_check("distance", "distance_sq_Tn_2_2", C(2,2), 0.d0, 1.d-14)

  if (test_qmckl_distance_sq == 0) return

  test_qmckl_distance_sq = QMCKL_FAILURE

  do j=1,n
     do i=1,m
        x =  (A(i,1)-B(1,j))**2 + &
             (A(i,2)-B(2,j))**2 + &
             (A(i,3)-B(3,j))**2
#ifndef VFC_CI
        if ( dabs(1.d0 - C(i,j)/x) > 1.d-14) return
#endif
     end do
  end do

  test_qmckl_distance_sq = &
       qmckl_distance_sq(context, 'n', 'N', m, n, A, LDA, B, LDB, C, LDC)

  vfc_err = qmckl_probe_check("distance", "distance_sq_nN_2_2", C(2,2), 0.d0, 1.d-14)

  test_qmckl_distance_sq = QMCKL_FAILURE

  do j=1,n
     do i=1,m
        x =  (A(1,i)-B(1,j))**2 + &
             (A(2,i)-B(2,j))**2 + &
             (A(3,i)-B(3,j))**2
        if ( dabs(1.d0 - C(i,j)/x) > 1.d-14 ) return
     end do
  end do

  test_qmckl_distance_sq = QMCKL_SUCCESS

  deallocate(A,B,C)
end function test_qmckl_distance_sq

integer(qmckl_exit_code) function test_qmckl_dist(context) bind(C)

  use qmckl
  use qmckl_verificarlo_f

  implicit none

  integer(qmckl_context), intent(in), value :: context
  logical(C_BOOL) :: vfc_err

  double precision, allocatable :: A(:,:), B(:,:), C(:,:)
  integer*8                     :: m, n, LDA, LDB, LDC
  double precision              :: x
  integer*8                     :: i,j

  m = 5
  n = 6
  LDA = m
  LDB = n
  LDC = 5

  allocate( A(LDA,m), B(LDB,n), C(LDC,n) )

  do j=1,m
     do i=1,m
        A(i,j) = -10.d0 + dble(i+j)
     end do
  end do
  do j=1,n
     do i=1,n
        B(i,j) = -1.d0 + dble(i*j)
     end do
  end do

  test_qmckl_dist = &
       qmckl_distance(context, 'X', 't', m, n, A, LDA, B, LDB, C, LDC)

  vfc_err = qmckl_probe("distance", "distance_Xt_2_2", C(2,2))

  if (test_qmckl_dist == 0) return

  test_qmckl_dist = &
       qmckl_distance(context, 't', 'X', m, n, A, LDA, B, LDB, C, LDC)

  vfc_err = qmckl_probe("distance", "distance_tX_2_2", C(2,2))

  if (test_qmckl_dist == 0) return

  test_qmckl_dist = &
       qmckl_distance(context, 'T', 't', m, n, A, LDA, B, LDB, C, LDC)

  vfc_err = qmckl_probe_check("distance", "distance_Tt_2_2", C(2,2), 0.d0, 1.d-14)

  if (test_qmckl_dist == 0) return


  test_qmckl_dist = QMCKL_FAILURE

  do j=1,n
     do i=1,m
        x =  dsqrt((A(i,1)-B(j,1))**2 + &
                   (A(i,2)-B(j,2))**2 + &
                   (A(i,3)-B(j,3))**2)
#ifndef VFC_CI
        if ( dabs(1.d0 - C(i,j)/x) > 1.d-14) return
#endif
     end do
  end do

  test_qmckl_dist = &
       qmckl_distance(context, 'n', 'T', m, n, A, LDA, B, LDB, C, LDC)

  vfc_err = qmckl_probe_check("distance", "distance_nT_2_2", C(2,2), 0.d0, 1.d-14)

  if (test_qmckl_dist == 0) return


  test_qmckl_dist = QMCKL_FAILURE

  do j=1,n
     do i=1,m
        x = dsqrt((A(1,i)-B(j,1))**2 + &
                  (A(2,i)-B(j,2))**2 + &
                  (A(3,i)-B(j,3))**2)
#ifndef VFC_CI
        if ( dabs(1.d0 - C(i,j)/x) > 1.d-14) return
#endif
     end do
  end do

  test_qmckl_dist = &
       qmckl_distance(context, 'T', 'n', m, n, A, LDA, B, LDB, C, LDC)

  vfc_err = qmckl_probe_check("distance", "distance_Tn_2_2", C(2,2), 0.d0, 1.d-14)

  if (test_qmckl_dist == 0) return


  test_qmckl_dist = QMCKL_FAILURE

  do j=1,n
     do i=1,m
        x =  dsqrt((A(i,1)-B(1,j))**2 + &
                   (A(i,2)-B(2,j))**2 + &
                   (A(i,3)-B(3,j))**2)
#ifndef VFC_CI
        if ( dabs(1.d0 - C(i,j)/x) > 1.d-14) return
#endif
     end do
  end do

  test_qmckl_dist = &
       qmckl_distance(context, 'n', 'N', m, n, A, LDA, B, LDB, C, LDC)

  vfc_err = qmckl_probe_check("distance", "distance_nN_2_2", C(2,2), 0.d0, 1.d-14)

  if (test_qmckl_dist == 0) return


  test_qmckl_dist = QMCKL_FAILURE

  do j=1,n
     do i=1,m
        x = dsqrt((A(1,i)-B(1,j))**2 + &
                  (A(2,i)-B(2,j))**2 + &
                  (A(3,i)-B(3,j))**2)
#ifndef VFC_CI
        if ( dabs(1.d0 - C(i,j)/x) > 1.d-14) return
#endif
     end do
  end do

  test_qmckl_dist = QMCKL_SUCCESS

  deallocate(A,B,C)
end function test_qmckl_dist

integer(qmckl_exit_code) function test_qmckl_dist_rescaled(context) bind(C)

  use qmckl

  implicit none

  integer(qmckl_context), intent(in), value :: context
  integer*8                     :: m, n, LDA, LDB, LDC
  double precision              :: x
  integer*8                     :: i,j

  double precision, parameter :: kappa = 0.6d0
  double precision, parameter :: kappa_inv = 1.d0/kappa

  integer*8, parameter :: elec_num = 10_8
  integer*8, parameter :: nucl_num = 2_8

  double precision :: nucl_coord(nucl_num,3) = reshape( (/ &
        0.0d0,  0.0d0 ,  &
        0.0d0,  0.0d0 ,  &
        0.0d0,  2.059801d0 /), shape(nucl_coord) )

  double precision :: elec_coord(3,elec_num) = reshape( (/ &
        -0.250655104764153d0    ,  0.503070975550133d0    ,  -0.166554344502303d0 , &
        -0.587812193472177d0    , -0.128751981129274d0    ,   0.187773606533075d0 , &
         1.61335569047166d0     , -0.615556732874863d0    ,  -1.43165470979934d0  , &
        -4.901239896295210d-003 , -1.120440036458986d-002 ,   1.99761909330422d0  , &
         0.766647499681200d0    , -0.293515395797937d0    ,   3.66454589201239d0  , &
        -0.127732483187947d0    , -0.138975497694196d0    ,  -8.669850480215846d-002 , &
        -0.232271834949124d0    , -1.059321673434182d-002 ,  -0.504862241464867d0 , &
         1.09360863531826d0     , -2.036103063808752d-003 ,  -2.702796910818986d-002 , &
        -0.108090166832043d0    ,  0.189161729653261d0    ,   2.15398313919894d0 , &
         0.397978144318712d0    , -0.254277292595981d0    ,   2.54553335476344d0  /), &
        shape(elec_coord))

  double precision :: ref_ee(elec_num,elec_num) = reshape( (/ &
       0.d0, 0.63475074d0, 1.29816415d0, 1.23147027d0, 1.51933127d0, &
       0.54402406d0, 0.51452479d0, 0.96538731d0, 1.25878564d0, 1.3722995d0 , &
       0.63475074d0, 0.d0, 1.35148664d0, 1.13524156d0, 1.48940503d0, &
       0.4582292d0, 0.62758076d0, 1.06560856d0, 1.179133d0, 1.30763703d0 , &
       1.29816415d0, 1.35148664d0, 0.d0, 1.50021375d0, 1.59200788d0, &
       1.23488312d0, 1.20844259d0, 1.0355537d0, 1.52064535d0, 1.53049239d0 , &
       1.23147027d0, 1.13524156d0, 1.50021375d0, 0.d0, 1.12016142d0, &
       1.19158954d0, 1.29762585d0, 1.24824277d0, 0.25292267d0, 0.58609336d0 , &
       1.51933127d0, 1.48940503d0, 1.59200788d0, 1.12016142d0, 0.d0, &
       1.50217017d0, 1.54012828d0, 1.48753895d0, 1.10441805d0, 0.84504381d0 , &
       0.54402406d0, 0.4582292d0, 1.23488312d0, 1.19158954d0, 1.50217017d0, &
       0.d0, 0.39417354d0, 0.87009603d0, 1.23838502d0, 1.33419121d0 , &
       0.51452479d0, 0.62758076d0, 1.20844259d0, 1.29762585d0, 1.54012828d0, &
       0.39417354d0, 0.d0, 0.95118809d0, 1.33068934d0, 1.41097406d0 , &
       0.96538731d0, 1.06560856d0, 1.0355537d0, 1.24824277d0, 1.48753895d0, &
       0.87009603d0, 0.95118809d0, 0.d0, 1.29422213d0, 1.33222493d0 , &
       1.25878564d0, 1.179133d0, 1.52064535d0, 0.25292267d0, 1.10441805d0, &
       1.23838502d0, 1.33068934d0, 1.29422213d0, 0.d0, 0.62196802d0 , &
       1.3722995d0, 1.30763703d0, 1.53049239d0, 0.58609336d0, 0.84504381d0, &
       1.33419121d0, 1.41097406d0, 1.33222493d0, 0.62196802d0, 0.d0 /), shape(ref_ee) )

  double precision :: ref_en(elec_num, nucl_num) = reshape( (/ &
       0.49421587d0, 0.52486545d0, 1.23280503d0, 1.16396998d0, 1.49156627d0, &
       0.1952946d0, 0.4726453d0, 0.80211227d0, 1.21198526d0, 1.31411513d0, &
       1.24641375d0, 1.15444238d0, 1.50565145d0, 0.06218339d0, 1.10153184d0, &
       1.20919677d0, 1.3111856d0, 1.26122875d0, 0.22122563d0, 0.55669168d0 /), shape(ref_en) )

  double precision, allocatable :: distance_ee(:,:), distance_en(:,:)

  allocate( distance_ee(elec_num,elec_num), distance_en(elec_num,nucl_num) )

  print *, 'ee'
  test_qmckl_dist_rescaled = &
       qmckl_distance_rescaled(context, 'N', 'N', elec_num, elec_num, elec_coord, &
       size(elec_coord,1)*1_8, elec_coord, size(elec_coord,1)*1_8, &
       distance_ee, size(distance_ee,1)*1_8, kappa)

  if (test_qmckl_dist_rescaled /= QMCKL_SUCCESS) return

  test_qmckl_dist_rescaled = QMCKL_FAILURE

  do j=1,elec_num
     do i=1,elec_num
        print *, i,j,real(distance_ee(i,j)), real(ref_ee(i,j))
        if (dabs(distance_ee(i,j) - ref_ee(i,j)) > 1.d-7) then
           return
        endif
     end do
  end do

  print *, 'en'
  test_qmckl_dist_rescaled = &
       qmckl_distance_rescaled(context, 'N', 'T', elec_num, nucl_num, elec_coord, &
       size(elec_coord,1)*1_8, nucl_coord, size(nucl_coord,1)*1_8, &
       distance_en, size(distance_en,1)*1_8, kappa)

  if (test_qmckl_dist_rescaled /= QMCKL_SUCCESS) return

  test_qmckl_dist_rescaled = QMCKL_FAILURE

  do j=1,nucl_num
     do i=1,elec_num
        print *, i,j,real(distance_en(i,j)), real(ref_en(i,j))
        if (dabs(distance_en(i,j) - ref_en(i,j)) > 1.d-7) then
           return
        endif
     end do
  end do

  test_qmckl_dist_rescaled = QMCKL_SUCCESS
end function test_qmckl_dist_rescaled
