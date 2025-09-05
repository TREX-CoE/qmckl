function test_qmckl_ao_gaussian_vgl(context) bind(C)
  use qmckl
  implicit none

  integer(qmckl_context), intent(in), value :: context
  integer(qmckl_exit_code) :: test_qmckl_ao_gaussian_vgl

  integer*8                     :: n, ldv, j, i
  double precision              :: X(3), R(3), Y(3), r2, z
  double precision, allocatable :: VGL(:,:), A(:)
  double precision              :: epsilon

  epsilon = 3.d0 * qmckl_get_numprec_epsilon(context)

  X = (/ 1.1 , 2.2 ,  3.3 /)
  R = (/ 0.1 , 1.2 , -2.3 /)
  Y(:) = X(:) - R(:)
  r2 = Y(1)**2 + Y(2)**2 + Y(3)**2

  n = 10;
  ldv = 100;

  allocate (A(n), VGL(ldv,5))
  do i=1,n
     A(i) = 0.0013 * dble(ishft(1,i))
  end do


  test_qmckl_ao_gaussian_vgl = &
       qmckl_ao_gaussian_vgl(context, X, R, n, A, VGL, ldv)
  if (test_qmckl_ao_gaussian_vgl /= 0) return

  test_qmckl_ao_gaussian_vgl = -1

  do i=1,n
     test_qmckl_ao_gaussian_vgl = -11
     z = dabs(1.d0 - VGL(i,1) / (dexp(-A(i) * r2)) )
     if ( z > epsilon ) then
        print *, z, epsilon
        return
     end if

     test_qmckl_ao_gaussian_vgl = -12
     z = dabs(1.d0 - VGL(i,2) / (&
          -2.d0 * A(i) * Y(1) * dexp(-A(i) * r2) ))
     if ( z > epsilon ) then
        print *, z, epsilon
        return
     end if

     test_qmckl_ao_gaussian_vgl = -13
     z = dabs(1.d0 - VGL(i,3) / (&
          -2.d0 * A(i) * Y(2) * dexp(-A(i) * r2) ))
     if ( z > epsilon ) then
        print *, z, epsilon
        return
     end if

     test_qmckl_ao_gaussian_vgl = -14
     z = dabs(1.d0 - VGL(i,4) / (&
          -2.d0 * A(i) * Y(3) * dexp(-A(i) * r2) ))
     if ( z > epsilon ) then
        print *, z, epsilon
        return
     end if

     test_qmckl_ao_gaussian_vgl = -15
     z = dabs(1.d0 - VGL(i,5) / (&
          A(i) * (4.d0*r2*A(i) - 6.d0) * dexp(-A(i) * r2) ))
     if ( z > epsilon ) then
        print *, z, epsilon
        return
     end if
  end do

  test_qmckl_ao_gaussian_vgl = 0

  deallocate(VGL)
end function test_qmckl_ao_gaussian_vgl

function test_qmckl_ao_slater_vgl(context) bind(C)
  use qmckl
  implicit none

  integer(qmckl_context), intent(in), value :: context
  integer(qmckl_exit_code) :: test_qmckl_ao_slater_vgl

  integer*8                     :: num_slater, ldv, j, i
  double precision              :: X(3), R(3), Y_vec(3), radius, z
  double precision, allocatable :: VGL(:,:), A(:)
  integer*8, allocatable        :: N(:)
  double precision              :: epsilon
  double precision              :: expected_val, expected_grad, expected_lapl
  double precision              :: alpha, n_val, rn, exp_alpha_r
  double precision              :: grad_coef, lapl_coef, r_inv, r_inv_2

  epsilon = 3.d0 * qmckl_get_numprec_epsilon(context)

  X = (/ 1.1d0, 2.2d0, 3.3d0 /)
  R = (/ 0.1d0, 1.2d0, -2.3d0 /)
  Y_vec(:) = X(:) - R(:)
  radius = dsqrt(Y_vec(1)**2 + Y_vec(2)**2 + Y_vec(3)**2)

  num_slater = 5
  ldv = 100

  allocate (A(num_slater), N(num_slater), VGL(ldv,5))
  do i=1,num_slater
     A(i) = 0.5d0 * dble(i)
     N(i) = i
  end do

  test_qmckl_ao_slater_vgl = &
       qmckl_ao_slater_vgl(context, X, R, num_slater, N, A, VGL, ldv)
  if (test_qmckl_ao_slater_vgl /= 0) return

  test_qmckl_ao_slater_vgl = -1

  r_inv = 1.d0 / radius
  r_inv_2 = r_inv * r_inv

  ! Test values and analytical derivatives
  do i=1,num_slater
     alpha = A(i)
     n_val = dble(N(i))
     rn = radius**N(i)
     exp_alpha_r = dexp(-alpha * radius)
     
     ! Test value: r^n * exp(-alpha * r)
     expected_val = rn * exp_alpha_r
     z = dabs(1.d0 - VGL(i,1) / expected_val)
     if (z > epsilon) then
        print *, 'Value test failed for i=', i, ', z=', z, ', epsilon=', epsilon
        return
     end if
     
     ! Test gradient components
     grad_coef = (n_val * r_inv - alpha) * expected_val
     expected_grad = grad_coef * Y_vec(1) * r_inv
     z = dabs(VGL(i,2) - expected_grad) / (dabs(expected_grad) + epsilon)
     if (z > epsilon) then
        print *, 'X-gradient test failed for i=', i, ', z=', z
        return
     end if
     
     expected_grad = grad_coef * Y_vec(2) * r_inv
     z = dabs(VGL(i,3) - expected_grad) / (dabs(expected_grad) + epsilon)
     if (z > epsilon) then
        print *, 'Y-gradient test failed for i=', i, ', z=', z
        return
     end if
     
     expected_grad = grad_coef * Y_vec(3) * r_inv
     z = dabs(VGL(i,4) - expected_grad) / (dabs(expected_grad) + epsilon)
     if (z > epsilon) then
        print *, 'Z-gradient test failed for i=', i, ', z=', z
        return
     end if
     
     ! Test Laplacian
     lapl_coef = n_val * (n_val - 1.d0) * r_inv_2 - 2.d0 * n_val * alpha * r_inv + alpha * alpha
     expected_lapl = lapl_coef * expected_val
     z = dabs(VGL(i,5) - expected_lapl) / (dabs(expected_lapl) + epsilon)
     if (z > epsilon) then
        print *, 'Laplacian test failed for i=', i, ', z=', z
        return
     end if
  end do

  test_qmckl_ao_slater_vgl = 0

  deallocate(VGL, A, N)
end function test_qmckl_ao_slater_vgl

function test_qmckl_ao_power(context) bind(C)
  use qmckl
  implicit none

  integer(qmckl_context), intent(in), value :: context
  integer(qmckl_exit_code) :: test_qmckl_ao_power

  integer*8                     :: n, LDP
  integer, allocatable          :: LMAX(:)
  double precision, allocatable :: X(:), P(:,:)
  integer*8                     :: i,j
  double precision              :: epsilon

  epsilon = qmckl_get_numprec_precision(context)

  n = 100;
  LDP = 10;

  allocate(X(n), P(LDP,n), LMAX(n))

  do j=1,n
     X(j) = -5.d0 + 0.1d0 * dble(j)
     LMAX(j) = 1 + int(mod(j, 5),4)
  end do

  test_qmckl_ao_power = qmckl_ao_power(context, n, X, LMAX, P, LDP)
  if (test_qmckl_ao_power /= QMCKL_SUCCESS) return

  test_qmckl_ao_power = QMCKL_FAILURE

  do j=1,n
     do i=1,LMAX(j)
        if ( X(j)**i == 0.d0 ) then
           if ( P(i,j) /= 0.d0) return
        else
           if ( dabs(1.d0 - P(i,j) / (X(j)**i)) > epsilon ) return
        end if
     end do
  end do

  test_qmckl_ao_power = QMCKL_SUCCESS
  deallocate(X,P,LMAX)
end function test_qmckl_ao_power

function test_qmckl_ao_polynomial_vgl(context) bind(C)
  use qmckl
  implicit none

  integer(qmckl_context), intent(in), value :: context

  integer                       :: lmax, d, i
  integer, allocatable          :: L(:,:)
  integer*8                     :: n, ldl, ldv, j
  double precision              :: X(3), R(3), Y(3)
  double precision, allocatable :: VGL(:,:)
  double precision              :: w
  double precision              :: epsilon
  integer(qmckl_exit_code) ::  test_qmckl_ao_polynomial_vgl

  epsilon = qmckl_get_numprec_precision(context)

  X = (/ 1.1 , 2.2 ,  3.3 /)
  R = (/ 0.1 , 1.2 , -2.3 /)
  Y(:) = X(:) - R(:)

  lmax = 4;
  ldl = 3;
  ldv = 100;

  d = (lmax+1)*(lmax+2)*(lmax+3)/6

  allocate (L(ldl,d), VGL(ldv,d))

  test_qmckl_ao_polynomial_vgl = &
       qmckl_ao_polynomial_vgl(context, X, R, lmax, n, L, ldl, VGL, ldv)

  if (test_qmckl_ao_polynomial_vgl /= QMCKL_SUCCESS) return
  if (n /= d) return

  do j=1,n
     test_qmckl_ao_polynomial_vgl = QMCKL_FAILURE
     do i=1,3
        if (L(i,j) < 0) return
     end do
     test_qmckl_ao_polynomial_vgl = QMCKL_FAILURE
     if (dabs(1.d0 - VGL(1,j) / (&
          Y(1)**L(1,j) * Y(2)**L(2,j) * Y(3)**L(3,j)  &
          )) > epsilon ) return

     test_qmckl_ao_polynomial_vgl = QMCKL_FAILURE
     if (L(1,j) < 1) then
        if (VGL(2,j) /= 0.d0) return
     else
        if (dabs(1.d0 - VGL(2,j) / (&
             L(1,j) * Y(1)**(L(1,j)-1) * Y(2)**L(2,j) * Y(3)**L(3,j) &
             )) > epsilon ) return
     end if

     test_qmckl_ao_polynomial_vgl = QMCKL_FAILURE
     if (L(2,j) < 1) then
        if (VGL(3,j) /= 0.d0) return
     else
        if (dabs(1.d0 - VGL(3,j) / (&
             L(2,j) * Y(1)**L(1,j) * Y(2)**(L(2,j)-1) * Y(3)**L(3,j) &
             )) > epsilon ) return
     end if

     test_qmckl_ao_polynomial_vgl = QMCKL_FAILURE
     if (L(3,j) < 1) then
        if (VGL(4,j) /= 0.d0) return
     else
        if (dabs(1.d0 - VGL(4,j) / (&
             L(3,j) * Y(1)**L(1,j) * Y(2)**L(2,j) * Y(3)**(L(3,j)-1) &
             )) > epsilon ) return
     end if

     test_qmckl_ao_polynomial_vgl = QMCKL_FAILURE
     w = 0.d0
     if (L(1,j) > 1) then
        w = w + L(1,j) * (L(1,j)-1) * Y(1)**(L(1,j)-2) * Y(2)**L(2,j) * Y(3)**L(3,j)
     end if
     if (L(2,j) > 1) then
        w = w + L(2,j) * (L(2,j)-1) * Y(1)**L(1,j) * Y(2)**(L(2,j)-2) * Y(3)**L(3,j)
     end if
     if (L(3,j) > 1) then
        w = w + L(3,j) * (L(3,j)-1) * Y(1)**L(1,j) * Y(2)**L(2,j) * Y(3)**(L(3,j)-2)
     end if
     if (w /= 0.d0) then
       if (dabs(1.d0 - VGL(5,j) / w) > epsilon ) return
     endif
  end do

  test_qmckl_ao_polynomial_vgl = QMCKL_SUCCESS

  deallocate(L,VGL)
end function test_qmckl_ao_polynomial_vgl
