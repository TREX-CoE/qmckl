function qmckl_distance_sq(context, transa, transb, m, n, &
     A, LDA, B, LDB, C, LDC) &
     bind(C) result(info)

  use qmckl_constants
  implicit none

  integer (qmckl_context) , intent(in)  , value :: context
  character(c_char)   , intent(in)  , value :: transa
  character(c_char)   , intent(in)  , value :: transb
  integer (c_int64_t) , intent(in)  , value :: m
  integer (c_int64_t) , intent(in)  , value :: n
  integer (c_int64_t) , intent(in)  , value :: lda
  integer (c_int64_t) , intent(in)  , value :: ldb
  integer (c_int64_t) , intent(in)  , value :: ldc
  real    (c_double ) , intent(in)          :: A(lda,*)
  real    (c_double ) , intent(in)          :: B(ldb,*)
  real    (c_double ) , intent(out)         :: C(ldc,n)

  integer(qmckl_exit_code) :: info

  integer*8 :: i,j
  real*8    :: x, y, z
  integer   :: transab

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

  if (transa == 'N' .or. transa == 'n') then
     transab = 0
  else if (transa == 'T' .or. transa == 't') then
     transab = 1
  else
     transab = -100
  endif

  if (transb == 'N' .or. transb == 'n') then
     continue
  else if (transb == 'T' .or. transb == 't') then
     transab = transab + 2
  else
     transab = -100
  endif

  if (transab < 0) then
     info = QMCKL_INVALID_ARG_1
     return
  endif

  if (iand(transab,1) == 0 .and. LDA < 3) then
     info = QMCKL_INVALID_ARG_7
     return
  endif

  if (iand(transab,1) == 1 .and. LDA < m) then
     info = QMCKL_INVALID_ARG_7
     return
  endif

  if (iand(transab,2) == 0 .and. LDB < 3) then
     info = QMCKL_INVALID_ARG_7
     return
  endif

  if (iand(transab,2) == 2 .and. LDB < n) then
     info = QMCKL_INVALID_ARG_7
     return
  endif


  select case (transab)

  case(0)

     do j=1,n
        do i=1,m
           x = A(1,i) - B(1,j)
           y = A(2,i) - B(2,j)
           z = A(3,i) - B(3,j)
           C(i,j) = x*x + y*y + z*z
        end do
     end do

  case(1)

     do j=1,n
        do i=1,m
           x = A(i,1) - B(1,j)
           y = A(i,2) - B(2,j)
           z = A(i,3) - B(3,j)
           C(i,j) = x*x + y*y + z*z
        end do
     end do

  case(2)

     do j=1,n
        do i=1,m
           x = A(1,i) - B(j,1)
           y = A(2,i) - B(j,2)
           z = A(3,i) - B(j,3)
           C(i,j) = x*x + y*y + z*z
        end do
     end do

  case(3)

     do j=1,n
        do i=1,m
           x = A(i,1) - B(j,1)
           y = A(i,2) - B(j,2)
           z = A(i,3) - B(j,3)
           C(i,j) = x*x + y*y + z*z
        end do
     end do

  end select

end function qmckl_distance_sq

function qmckl_distance(context, transa, transb, m, n, &
     A, LDA, B, LDB, C, LDC) &
     bind(C) result(info)
  use qmckl_constants
  implicit none
  integer(qmckl_context), intent(in), value  :: context
  character(c_char)   , intent(in)  , value :: transa
  character(c_char)   , intent(in)  , value :: transb
  integer (c_int64_t) , intent(in)  , value :: m
  integer (c_int64_t) , intent(in)  , value :: n
  integer (c_int64_t) , intent(in)  , value :: lda
  integer (c_int64_t) , intent(in)  , value :: ldb
  integer (c_int64_t) , intent(in)  , value :: ldc
  real    (c_double ) , intent(in)          :: A(lda,*)
  real    (c_double ) , intent(in)          :: B(ldb,*)
  real    (c_double ) , intent(out)         :: C(ldc,n)
  integer (qmckl_exit_code) :: info

  integer*8 :: i,j
  real*8    :: x, y, z
  integer   :: transab

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

  if (transa == 'N' .or. transa == 'n') then
     transab = 0
  else if (transa == 'T' .or. transa == 't') then
     transab = 1
  else
     transab = -100
  endif

  if (transb == 'N' .or. transb == 'n') then
     continue
  else if (transb == 'T' .or. transb == 't') then
     transab = transab + 2
  else
     transab = -100
  endif

  if (transab < 0) then
     info = QMCKL_INVALID_ARG_1
     return
  endif

  ! check for LDA
  if (iand(transab,1) == 0 .and. LDA < 3) then
     info = QMCKL_INVALID_ARG_7
     return
  endif

  if (iand(transab,1) == 1 .and. LDA < m) then
     info = QMCKL_INVALID_ARG_7
     return
  endif

  ! check for LDB
  if (iand(transab,1) == 0 .and. LDB < 3) then
     info = QMCKL_INVALID_ARG_9
     return
  endif

  if (iand(transab,1) == 1 .and. LDB < n) then
     info = QMCKL_INVALID_ARG_9
     return
  endif

  ! check for LDC
  if (LDC < m) then
     info = QMCKL_INVALID_ARG_11
     return
  endif


  select case (transab)

  case(0)

     do j=1,n
        do i=1,m
           x = A(1,i) - B(1,j)
           y = A(2,i) - B(2,j)
           z = A(3,i) - B(3,j)
           C(i,j) = x*x + y*y + z*z
        end do
        C(:,j) = dsqrt(C(:,j))
     end do

  case(1)

     do j=1,n
        do i=1,m
           x = A(i,1) - B(1,j)
           y = A(i,2) - B(2,j)
           z = A(i,3) - B(3,j)
           C(i,j) = x*x + y*y + z*z
        end do
        C(:,j) = dsqrt(C(:,j))
     end do

  case(2)

     do j=1,n
        do i=1,m
           x = A(1,i) - B(j,1)
           y = A(2,i) - B(j,2)
           z = A(3,i) - B(j,3)
           C(i,j) = x*x + y*y + z*z
        end do
        C(:,j) = dsqrt(C(:,j))
     end do

  case(3)

     do j=1,n
        do i=1,m
           x = A(i,1) - B(j,1)
           y = A(i,2) - B(j,2)
           z = A(i,3) - B(j,3)
           C(i,j) = x*x + y*y + z*z
        end do
        C(:,j) = dsqrt(C(:,j))
     end do

  end select

end function qmckl_distance

function qmckl_distance_rescaled(context, transa, transb, m, n, &
     A, LDA, B, LDB, C, LDC, rescale_factor_kappa) &
     bind(C) result(info)
  use qmckl_constants
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
  integer(qmckl_exit_code) :: info

  integer*8 :: i,j
  real*8    :: x, y, z, dist, rescale_factor_kappa_inv
  integer   :: transab

  rescale_factor_kappa_inv = 1.0d0/rescale_factor_kappa;

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

  if (transa == 'N' .or. transa == 'n') then
     transab = 0
  else if (transa == 'T' .or. transa == 't') then
     transab = 1
  else
     transab = -100
  endif

  if (transb == 'N' .or. transb == 'n') then
     continue
  else if (transb == 'T' .or. transb == 't') then
     transab = transab + 2
  else
     transab = -100
  endif

  ! check for LDA
  if (transab < 0) then
     info = QMCKL_INVALID_ARG_1
     return
  endif

  if (iand(transab,1) == 0 .and. LDA < 3) then
     info = QMCKL_INVALID_ARG_7
     return
  endif

  if (iand(transab,1) == 1 .and. LDA < m) then
     info = QMCKL_INVALID_ARG_7
     return
  endif

  ! check for LDB
  if (iand(transab,2) == 0 .and. LDB < 3) then
     info = QMCKL_INVALID_ARG_9
     return
  endif

  if (iand(transab,2) == 2 .and. LDB < n) then
     info = QMCKL_INVALID_ARG_9
     return
  endif

  ! check for LDC
  if (LDC < m) then
     info = QMCKL_INVALID_ARG_11
     return
  endif


  select case (transab)

  case(0)

     do j=1,n
        do i=1,m
           x = A(1,i) - B(1,j)
           y = A(2,i) - B(2,j)
           z = A(3,i) - B(3,j)
           dist = dsqrt(x*x + y*y + z*z)
           C(i,j) = (1.0d0 - dexp(-rescale_factor_kappa * dist)) * rescale_factor_kappa_inv
        end do
     end do

  case(1)

     do j=1,n
        do i=1,m
           x = A(i,1) - B(1,j)
           y = A(i,2) - B(2,j)
           z = A(i,3) - B(3,j)
           dist = dsqrt(x*x + y*y + z*z)
           C(i,j) = (1.0d0 - dexp(-rescale_factor_kappa * dist)) * rescale_factor_kappa_inv
        end do
     end do

  case(2)

     do j=1,n
        do i=1,m
           x = A(1,i) - B(j,1)
           y = A(2,i) - B(j,2)
           z = A(3,i) - B(j,3)
           dist = dsqrt(x*x + y*y + z*z)
           C(i,j) = (1.0d0 - dexp(-rescale_factor_kappa * dist)) * rescale_factor_kappa_inv
        end do
     end do

  case(3)

     do j=1,n
        do i=1,m
           x = A(i,1) - B(j,1)
           y = A(i,2) - B(j,2)
           z = A(i,3) - B(j,3)
           dist = dsqrt(x*x + y*y + z*z)
           C(i,j) = (1.0d0 - dexp(-rescale_factor_kappa * dist)) * rescale_factor_kappa_inv
        end do
     end do

  end select

end function qmckl_distance_rescaled

function qmckl_distance_rescaled_gl(context, transa, transb, m, n, &
     A, LDA, B, LDB, C, LDC, rescale_factor_kappa) &
     bind(C) result(info)
  use qmckl_constants
  implicit none

  integer(qmckl_exit_code) :: info
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

  integer*8 :: i,j
  real*8    :: x, y, z, dist, dist_inv
  real*8    :: rij
  integer   :: transab

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

  if (transa == 'N' .or. transa == 'n') then
     transab = 0
  else if (transa == 'T' .or. transa == 't') then
     transab = 1
  else
     transab = -100
  endif

  if (transb == 'N' .or. transb == 'n') then
     continue
  else if (transb == 'T' .or. transb == 't') then
     transab = transab + 2
  else
     transab = -100
  endif

  ! check for LDA
  if (transab < 0) then
     info = QMCKL_INVALID_ARG_1
     return
  endif

  if (iand(transab,1) == 0 .and. LDA < 3) then
     info = QMCKL_INVALID_ARG_7
     return
  endif

  if (iand(transab,1) == 1 .and. LDA < m) then
     info = QMCKL_INVALID_ARG_7
     return
  endif

  ! check for LDB
  if (iand(transab,2) == 0 .and. LDB < 3) then
     info = QMCKL_INVALID_ARG_9
     return
  endif

  if (iand(transab,2) == 2 .and. LDB < n) then
     info = QMCKL_INVALID_ARG_9
     return
  endif

  ! check for LDC
  if (LDC < m) then
     info = QMCKL_INVALID_ARG_11
     return
  endif


  select case (transab)

  case(0)

     do j=1,n
        do i=1,m
           x = A(1,i) - B(1,j)
           y = A(2,i) - B(2,j)
           z = A(3,i) - B(3,j)
           dist = max(1.d-20, dsqrt(x*x + y*y + z*z))
           dist_inv = 1.0d0/dist
           rij = dexp(-rescale_factor_kappa * dist)
           C(1,i,j) = x * dist_inv * rij
           C(2,i,j) = y * dist_inv * rij
           C(3,i,j) = z * dist_inv * rij
           C(4,i,j) = (2.0d0 * dist_inv - rescale_factor_kappa) * rij
        end do
     end do

  case(1)

     do j=1,n
        do i=1,m
           x = A(i,1) - B(1,j)
           y = A(i,2) - B(2,j)
           z = A(i,3) - B(3,j)
           dist = max(1.d-20, dsqrt(x*x + y*y + z*z))
           dist_inv = 1.0d0/dist
           rij = dexp(-rescale_factor_kappa * dist)
           C(1,i,j) = x * dist_inv * rij
           C(2,i,j) = y * dist_inv * rij
           C(3,i,j) = z * dist_inv * rij
           C(4,i,j) = (2.0d0 * dist_inv - rescale_factor_kappa) * rij
        end do
     end do

  case(2)

     do j=1,n
        do i=1,m
           x = A(1,i) - B(j,1)
           y = A(2,i) - B(j,2)
           z = A(3,i) - B(j,3)
           dist = max(1.d-20, dsqrt(x*x + y*y + z*z))
           dist_inv = 1.0d0/dist
           rij = dexp(-rescale_factor_kappa * dist)
           C(1,i,j) = x * dist_inv * rij
           C(2,i,j) = y * dist_inv * rij
           C(3,i,j) = z * dist_inv * rij
           C(4,i,j) = (2.0d0 * dist_inv - rescale_factor_kappa) * rij
        end do
     end do

  case(3)

     do j=1,n
        do i=1,m
           x = A(i,1) - B(j,1)
           y = A(i,2) - B(j,2)
           z = A(i,3) - B(j,3)
           dist = max(1.d-20, dsqrt(x*x + y*y + z*z))
           dist_inv = 1.0d0/dist
           rij = dexp(-rescale_factor_kappa * dist)
           C(1,i,j) = x * dist_inv * rij
           C(2,i,j) = y * dist_inv * rij
           C(3,i,j) = z * dist_inv * rij
           C(4,i,j) = (2.0d0 * dist_inv - rescale_factor_kappa) * rij
        end do
     end do

  end select

end function qmckl_distance_rescaled_gl
