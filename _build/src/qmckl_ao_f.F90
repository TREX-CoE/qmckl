#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

function qmckl_ao_gaussian_vgl(context, X, R, n, A, VGL, ldv) &
  bind(C) result(info)
  use qmckl_constants
  implicit none
  integer (qmckl_context) , intent(in) , value :: context
  real    (c_double)  , intent(in)         :: X(3), R(3)
  integer (c_int64_t) , intent(in) , value :: n
  integer (c_int64_t) , intent(in) , value :: ldv
  real    (c_double)  , intent(in)         :: A(n)
  real    (c_double)  , intent(out)        :: VGL(ldv,5)
  integer (qmckl_exit_code) :: info

  integer*8         :: i,j
  double precision  :: Y(3), r2, t, u, v

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (n <= 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  if (ldv < n) then
     info = QMCKL_INVALID_ARG_7
     return
  endif


  do i=1,3
     Y(i) = X(i) - R(i)
  end do
  r2 = Y(1)*Y(1) + Y(2)*Y(2) + Y(3)*Y(3)

  do i=1,n
     VGL(i,1) = dexp(-A(i) * r2)
  end do

  do i=1,n
     VGL(i,5) = A(i) * VGL(i,1)
  end do

  t = -2.d0 * ( X(1) - R(1) )
  u = -2.d0 * ( X(2) - R(2) )
  v = -2.d0 * ( X(3) - R(3) )

  do i=1,n
     VGL(i,2) = t * VGL(i,5)
     VGL(i,3) = u * VGL(i,5)
     VGL(i,4) = v * VGL(i,5)
  end do

  t = 4.d0 * r2
  do i=1,n
     VGL(i,5) = (t * A(i) - 6.d0) *  VGL(i,5)
  end do

end function qmckl_ao_gaussian_vgl

function qmckl_ao_slater_vgl(context, X, R, num_slater, N, A, VGL, ldv) &
  bind(C) result(info)
  use qmckl_constants
  implicit none
  integer (qmckl_context) , intent(in) , value :: context
  real    (c_double)  , intent(in)         :: X(3), R(3)
  integer (c_int64_t) , intent(in) , value :: num_slater
  integer (c_int64_t) , intent(in) , value :: ldv
  integer (c_int64_t) , intent(in)         :: N(num_slater)
  real    (c_double)  , intent(in)         :: A(num_slater)
  real    (c_double)  , intent(out)        :: VGL(ldv,5)
  integer (qmckl_exit_code) :: info

  integer*8         :: i
  double precision  :: Y_vec(3), radius, r_inv, r_inv_2, x_over_r, y_over_r, z_over_r
  double precision  :: alpha, n_val, rn, exp_alpha_r
  double precision  :: grad_coef, lapl_coef

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (num_slater <= 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  if (ldv < num_slater) then
     info = QMCKL_INVALID_ARG_8
     return
  endif

  ! Compute distance vector and radius
  do i=1,3
     Y_vec(i) = X(i) - R(i)
  end do
  radius = dsqrt(Y_vec(1)*Y_vec(1) + Y_vec(2)*Y_vec(2) + Y_vec(3)*Y_vec(3))

  ! Handle r = 0 case
  if (radius < 1.d-15) then
     do i=1,num_slater
        if (N(i) == 1) then
           VGL(i,1) = 1.d0
        else
           VGL(i,1) = 0.d0
        end if
        VGL(i,2) = 0.d0
        VGL(i,3) = 0.d0
        VGL(i,4) = 0.d0
        VGL(i,5) = 0.d0
     end do
     return
  end if

  r_inv = 1.d0 / radius
  r_inv_2 = r_inv * r_inv
  x_over_r = Y_vec(1) * r_inv
  y_over_r = Y_vec(2) * r_inv
  z_over_r = Y_vec(3) * r_inv

  ! Compute values, gradients and Laplacians for each Slater orbital
  do i=1,num_slater
     alpha = A(i)
     n_val = dble(N(i))

     ! Compute r^n * exp(-alpha * r)
     if (N(i) == 1) then
        rn = radius
     else
        rn = radius**N(i)
     end if
     exp_alpha_r = dexp(-alpha * radius)

     ! Value: r^n * exp(-alpha * r)
     VGL(i,1) = rn * exp_alpha_r

     ! Gradient coefficients
     grad_coef = (n_val * r_inv - alpha) * VGL(i,1)
     VGL(i,2) = grad_coef * x_over_r
     VGL(i,3) = grad_coef * y_over_r
     VGL(i,4) = grad_coef * z_over_r

     ! Laplacian: 1/r² [ n² - 2nar + ar(ar-2) + n] * v
     lapl_coef = r_inv*r_inv * (n_val*n_val - 2.d0 * n_val * alpha*radius + &
          alpha*radius * (alpha*radius - 2.d0) + n_val)
     VGL(i,5) = lapl_coef * VGL(i,1)
  end do

end function qmckl_ao_slater_vgl

function qmckl_compute_ao_basis_primitive_gaussian_vgl &
     (context, prim_num, point_num, nucl_num, nucleus_prim_index, coord, nucl_coord, expo, primitive_vgl) &
     bind(C) result(info)

  use qmckl_constants

  use qmckl, only: qmckl_get_numprec_epsilon
  implicit none
  integer (qmckl_context), intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: prim_num
  integer (c_int64_t) , intent(in)  , value :: point_num
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  integer (c_int64_t) , intent(in)          :: nucleus_prim_index(nucl_num+1)
  real    (c_double ) , intent(in)          :: coord(point_num,3)
  real    (c_double ) , intent(in)          :: nucl_coord(nucl_num,3)
  real    (c_double ) , intent(in)          :: expo(prim_num)
  real    (c_double ) , intent(out)         :: primitive_vgl(prim_num,5,point_num)
  integer(qmckl_exit_code) :: info

  integer*8 :: inucl, iprim, ipoint
  double precision :: x, y, z, two_a, ar2, r2, v, cutoff

  info = QMCKL_SUCCESS

  ! Don't compute exponentials when the result will be almost zero.
  cutoff = -dlog(qmckl_get_numprec_epsilon(context))

  do inucl=1,nucl_num
     ! C is zero-based, so shift bounds by one
     do iprim = nucleus_prim_index(inucl)+1, nucleus_prim_index(inucl+1)
        do ipoint = 1, point_num
           x = coord(ipoint,1) - nucl_coord(inucl,1)
           y = coord(ipoint,2) - nucl_coord(inucl,2)
           z = coord(ipoint,3) - nucl_coord(inucl,3)

           r2 = x*x + y*y + z*z
           ar2 = expo(iprim)*r2
           if (ar2 > cutoff) cycle

           v = dexp(-ar2)
           two_a = -2.d0 * expo(iprim) * v

           primitive_vgl(iprim, 1, ipoint) = v
           primitive_vgl(iprim, 2, ipoint) = two_a * x
           primitive_vgl(iprim, 3, ipoint) = two_a * y
           primitive_vgl(iprim, 4, ipoint) = two_a * z
           primitive_vgl(iprim, 5, ipoint) = two_a * (3.d0 - 2.d0*ar2)

        end do
     end do
  end do

end function qmckl_compute_ao_basis_primitive_gaussian_vgl

function qmckl_compute_ao_basis_shell_gaussian_vgl( &
     context, prim_num, shell_num, point_num, nucl_num,       &
     nucleus_shell_num, nucleus_index, nucleus_range,         &
     shell_prim_index, shell_prim_num, coord, nucl_coord,     &
     expo, coef_normalized, shell_vgl)                        &
     bind(C) result(info)

  use qmckl_constants
  use qmckl, only: qmckl_get_numprec_epsilon

  implicit none
  integer (qmckl_context), intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: prim_num
  integer (c_int64_t) , intent(in)  , value :: shell_num
  integer (c_int64_t) , intent(in)  , value :: point_num
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  integer (c_int64_t) , intent(in)          :: nucleus_shell_num(nucl_num)
  integer (c_int64_t) , intent(in)          :: nucleus_index(nucl_num)
  real    (c_double ) , intent(in)          :: nucleus_range(nucl_num)
  integer (c_int64_t) , intent(in)          :: shell_prim_index(shell_num)
  integer (c_int64_t) , intent(in)          :: shell_prim_num(shell_num)
  real    (c_double ) , intent(in)          :: coord(point_num,3)
  real    (c_double ) , intent(in)          :: nucl_coord(nucl_num,3)
  real    (c_double ) , intent(in)          :: expo(prim_num)
  real    (c_double ) , intent(in)          :: coef_normalized(prim_num)
  real    (c_double ) , intent(out)         :: shell_vgl(shell_num,5,point_num)
  integer(qmckl_exit_code) :: info

  integer*8 :: inucl, iprim, ipoint, ishell
  integer*8 :: ishell_start, ishell_end
  integer*8 :: iprim_start , iprim_end
  double precision :: x, y, z, two_a, ar2, r2, v, cutoff

  info = QMCKL_SUCCESS

  ! Don't compute exponentials when the result will be almost zero.
  cutoff = -dlog(qmckl_get_numprec_epsilon(context))

  do ipoint = 1, point_num

     do inucl=1,nucl_num

        x = coord(ipoint,1) - nucl_coord(inucl,1)
        y = coord(ipoint,2) - nucl_coord(inucl,2)
        z = coord(ipoint,3) - nucl_coord(inucl,3)

        r2 = x*x + y*y + z*z

        if (r2 > cutoff*nucleus_range(inucl)) then
           cycle
        end if

        ! C is zero-based, so shift bounds by one
        ishell_start = nucleus_index(inucl) + 1
        ishell_end   = nucleus_index(inucl) + nucleus_shell_num(inucl)

        do ishell=ishell_start, ishell_end

           shell_vgl(ishell, 1, ipoint) = 0.d0
           shell_vgl(ishell, 2, ipoint) = 0.d0
           shell_vgl(ishell, 3, ipoint) = 0.d0
           shell_vgl(ishell, 4, ipoint) = 0.d0
           shell_vgl(ishell, 5, ipoint) = 0.d0

           iprim_start = shell_prim_index(ishell) + 1
           iprim_end   = shell_prim_index(ishell) + shell_prim_num(ishell)

           do iprim = iprim_start, iprim_end

              ar2 = expo(iprim)*r2
              if (ar2 > cutoff) then
                 cycle
              end if

              v = coef_normalized(iprim) * dexp(-ar2)
              two_a = -2.d0 * expo(iprim) * v

              shell_vgl(ishell, 1, ipoint) = &
                   shell_vgl(ishell, 1, ipoint) + v

              shell_vgl(ishell, 2, ipoint) = &
                   shell_vgl(ishell, 2, ipoint) + two_a * x

              shell_vgl(ishell, 3, ipoint) = &
                   shell_vgl(ishell, 3, ipoint) + two_a * y

              shell_vgl(ishell, 4, ipoint) = &
                   shell_vgl(ishell, 4, ipoint) + two_a * z

              shell_vgl(ishell, 5, ipoint) = &
                   shell_vgl(ishell, 5, ipoint) + two_a * (3.d0 - 2.d0*ar2)

           end do

        end do
     end do

  end do

end function qmckl_compute_ao_basis_shell_gaussian_vgl

function qmckl_compute_ao_basis_shell_slater_vgl( &
     context, prim_num, shell_num, point_num, nucl_num,       &
     nucleus_shell_num, nucleus_index, nucleus_range,         &
     shell_prim_index, shell_prim_num, r_power,               &
     coord, nucl_coord, expo, coef_normalized, shell_vgl)     &
     bind(C) result(info)

  use qmckl_constants
  use qmckl, only: qmckl_get_numprec_epsilon

  implicit none
  integer (qmckl_context), intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: prim_num
  integer (c_int64_t) , intent(in)  , value :: shell_num
  integer (c_int64_t) , intent(in)  , value :: point_num
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  integer (c_int64_t) , intent(in)          :: nucleus_shell_num(nucl_num)
  integer (c_int64_t) , intent(in)          :: nucleus_index(nucl_num)
  real    (c_double ) , intent(in)          :: nucleus_range(nucl_num)
  integer (c_int64_t) , intent(in)          :: shell_prim_index(shell_num)
  integer (c_int64_t) , intent(in)          :: shell_prim_num(shell_num)
  integer (c_int32_t) , intent(in)          :: r_power(shell_num)
  real    (c_double ) , intent(in)          :: coord(point_num,3)
  real    (c_double ) , intent(in)          :: nucl_coord(nucl_num,3)
  real    (c_double ) , intent(in)          :: expo(prim_num)
  real    (c_double ) , intent(in)          :: coef_normalized(prim_num)
  real    (c_double ) , intent(out)         :: shell_vgl(shell_num,5,point_num)
  integer(qmckl_exit_code) :: info

  integer*8 :: inucl, iprim, ipoint, ishell
  integer*8 :: ishell_start, ishell_end
  integer*8 :: iprim_start , iprim_end
  double precision :: x, y, z, r, r_inv, ar, n_val, v, cutoff
  double precision :: grad_coef, lapl_coef

  info = QMCKL_SUCCESS

  ! Don't compute exponentials when the result will be almost zero.
  cutoff = -dlog(qmckl_get_numprec_epsilon(context))

  do ipoint = 1, point_num

     do inucl=1,nucl_num

        x = coord(ipoint,1) - nucl_coord(inucl,1)
        y = coord(ipoint,2) - nucl_coord(inucl,2)
        z = coord(ipoint,3) - nucl_coord(inucl,3)

        r = dsqrt(x*x + y*y + z*z)

        if (r > cutoff*nucleus_range(inucl)) then
           cycle
        end if

        ! C is zero-based, so shift bounds by one
        ishell_start = nucleus_index(inucl) + 1
        ishell_end   = nucleus_index(inucl) + nucleus_shell_num(inucl)

        do ishell=ishell_start, ishell_end

           shell_vgl(ishell, 1, ipoint) = 0.d0
           shell_vgl(ishell, 2, ipoint) = 0.d0
           shell_vgl(ishell, 3, ipoint) = 0.d0
           shell_vgl(ishell, 4, ipoint) = 0.d0
           shell_vgl(ishell, 5, ipoint) = 0.d0

           iprim_start = shell_prim_index(ishell) + 1
           iprim_end   = shell_prim_index(ishell) + shell_prim_num(ishell)

           ! Handle r = 0 case for Slater orbitals
           if (r < 1.d-15) then
              do iprim = iprim_start, iprim_end
                 n_val = dble(r_power(ishell) + 1)  ! Principal quantum number
                 if (n_val == 1.d0) then
                    shell_vgl(ishell, 1, ipoint) = shell_vgl(ishell, 1, ipoint) + coef_normalized(iprim)
                 end if
                 ! Gradients and Laplacian are zero at r=0 for all n>1
              end do
              cycle
           end if

           r_inv = 1.d0 / r

           do iprim = iprim_start, iprim_end

              ar = expo(iprim) * r
              if (ar > cutoff) then
                 cycle
              end if

              n_val = dble(r_power(ishell) + 1)  ! Principal quantum number for Slater orbitals

              ! Compute r^n * exp(-a*r)
              v = coef_normalized(iprim) * (r ** int(n_val)) * dexp(-ar)

              ! Value
              shell_vgl(ishell, 1, ipoint) = shell_vgl(ishell, 1, ipoint) + v

              ! Gradient coefficient: (n/r - a)
              grad_coef = (n_val * r_inv - expo(iprim)) * v

              ! Gradients: gradient_coefficient * (x/r, y/r, z/r)
              shell_vgl(ishell, 2, ipoint) = shell_vgl(ishell, 2, ipoint) + grad_coef * x * r_inv
              shell_vgl(ishell, 3, ipoint) = shell_vgl(ishell, 3, ipoint) + grad_coef * y * r_inv
              shell_vgl(ishell, 4, ipoint) = shell_vgl(ishell, 4, ipoint) + grad_coef * z * r_inv

              ! Laplacian: 1/r² [ n² - 2nar + ar(ar-2) + n] * v
              lapl_coef = r_inv*r_inv * ( &
                   n_val*n_val - 2.d0 * n_val * ar + &
                   ar * (ar - 2.d0) + n_val)

              shell_vgl(ishell, 5, ipoint) = shell_vgl(ishell, 5, ipoint) + lapl_coef * v

           end do

        end do
     end do

  end do

end function qmckl_compute_ao_basis_shell_slater_vgl

function qmckl_compute_ao_basis_shell_gaussian_hessian( &
     context, prim_num, shell_num, point_num, nucl_num,       &
     nucleus_shell_num, nucleus_index, nucleus_range,         &
     shell_prim_index, shell_prim_num, coord, nucl_coord,     &
     expo, coef_normalized, shell_hessian)                        &
     bind(C) result(info)

  use qmckl_constants
  use qmckl, only: qmckl_get_numprec_epsilon

  implicit none
  integer (qmckl_context), intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: prim_num
  integer (c_int64_t) , intent(in)  , value :: shell_num
  integer (c_int64_t) , intent(in)  , value :: point_num
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  integer (c_int64_t) , intent(in)          :: nucleus_shell_num(nucl_num)
  integer (c_int64_t) , intent(in)          :: nucleus_index(nucl_num)
  real    (c_double ) , intent(in)          :: nucleus_range(nucl_num)
  integer (c_int64_t) , intent(in)          :: shell_prim_index(shell_num)
  integer (c_int64_t) , intent(in)          :: shell_prim_num(shell_num)
  real    (c_double ) , intent(in)          :: coord(point_num,3)
  real    (c_double ) , intent(in)          :: nucl_coord(nucl_num,3)
  real    (c_double ) , intent(in)          :: expo(prim_num)
  real    (c_double ) , intent(in)          :: coef_normalized(prim_num)
  real    (c_double ) , intent(out)         :: shell_hessian(shell_num,4,3,point_num)
  integer(qmckl_exit_code) :: info

  double precision :: xyz(3)

  integer*8 :: inucl, iprim, ipoint, ishell
  integer*8 :: ishell_start, ishell_end
  integer*8 :: iprim_start , iprim_end, i, j
  double precision :: x, y, z, two_a, ar2, r2, v, cutoff

  info = QMCKL_SUCCESS

  ! Don't compute exponentials when the result will be almost zero.
  cutoff = -dlog(qmckl_get_numprec_epsilon(context))

  do ipoint = 1, point_num

     do inucl=1,nucl_num

        xyz(1) = coord(ipoint,1) - nucl_coord(inucl,1)
        xyz(2) = coord(ipoint,2) - nucl_coord(inucl,2)
        xyz(3) = coord(ipoint,3) - nucl_coord(inucl,3)

        r2 = xyz(1)*xyz(1) + xyz(2)*xyz(2) + xyz(3)*xyz(3)

        if (r2 > cutoff*nucleus_range(inucl)) then
           cycle
        end if

        ! C is zero-based, so shift bounds by one
        ishell_start = nucleus_index(inucl) + 1
        ishell_end   = nucleus_index(inucl) + nucleus_shell_num(inucl)

        do ishell=ishell_start, ishell_end

           do i = 1, 4
              do j = 1, 3
                 shell_hessian(ishell, i, j, ipoint) = 0.d0
              end do
           end do

           iprim_start = shell_prim_index(ishell) + 1
           iprim_end   = shell_prim_index(ishell) + shell_prim_num(ishell)

           do iprim = iprim_start, iprim_end

              ar2 = expo(iprim)*r2
              if (ar2 > cutoff) then
                 cycle
              end if

              v = coef_normalized(iprim) * dexp(-ar2)
              two_a = 2.d0 * expo(iprim)

              do i = 1, 3
                 do j = 1, 3
                    if (i == j) then
                        shell_hessian(ishell, i, j, ipoint) = &
                              shell_hessian(ishell, i, j, ipoint) - two_a * v
                    end if
                    shell_hessian(ishell, i, j, ipoint) = &
                         shell_hessian(ishell, i, j, ipoint) + two_a * two_a * xyz(i) * xyz(j) * v

                 end do
                 shell_hessian(ishell,4,i,ipoint) = shell_hessian(ishell,4,i,ipoint) &
                      + (5.d0 * two_a * two_a * xyz(i) &
                      - two_a * two_a * two_a * r2 * xyz(i)) * v
              end do

           end do

        end do
     end do

  end do

end function qmckl_compute_ao_basis_shell_gaussian_hessian

function qmckl_compute_ao_basis_shell_slater_hessian( &
     context, prim_num, shell_num, point_num, nucl_num,       &
     nucleus_shell_num, nucleus_index, nucleus_range,         &
     shell_prim_index, shell_prim_num, r_power,               &
     coord, nucl_coord, expo, coef_normalized, shell_hessian) &
     bind(C) result(info)

  use qmckl_constants
  use qmckl, only: qmckl_get_numprec_epsilon

  implicit none
  integer (qmckl_context), intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: prim_num
  integer (c_int64_t) , intent(in)  , value :: shell_num
  integer (c_int64_t) , intent(in)  , value :: point_num
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  integer (c_int64_t) , intent(in)          :: nucleus_shell_num(nucl_num)
  integer (c_int64_t) , intent(in)          :: nucleus_index(nucl_num)
  real    (c_double ) , intent(in)          :: nucleus_range(nucl_num)
  integer (c_int64_t) , intent(in)          :: shell_prim_index(shell_num)
  integer (c_int64_t) , intent(in)          :: shell_prim_num(shell_num)
  integer (c_int32_t) , intent(in)          :: r_power(shell_num)
  real    (c_double ) , intent(in)          :: coord(point_num,3)
  real    (c_double ) , intent(in)          :: nucl_coord(nucl_num,3)
  real    (c_double ) , intent(in)          :: expo(prim_num)
  real    (c_double ) , intent(in)          :: coef_normalized(prim_num)
  real    (c_double ) , intent(out)         :: shell_hessian(shell_num,4,3,point_num)
  integer(qmckl_exit_code) :: info

  integer*8 :: inucl, iprim, ipoint, ishell
  integer*8 :: ishell_start, ishell_end
  integer*8 :: iprim_start , iprim_end
  double precision :: x, y, z, r, r_inv, r_inv2, ar, n_val, v, cutoff
  double precision :: alpha, grad_coef, hess_diag, hess_off

  info = QMCKL_SUCCESS

  ! Don't compute exponentials when the result will be almost zero.
  cutoff = -dlog(qmckl_get_numprec_epsilon(context))

  do ipoint = 1, point_num

     do inucl=1,nucl_num

        x = coord(ipoint,1) - nucl_coord(inucl,1)
        y = coord(ipoint,2) - nucl_coord(inucl,2)
        z = coord(ipoint,3) - nucl_coord(inucl,3)

        r = dsqrt(x*x + y*y + z*z)

        if (r > nucleus_range(inucl)) then
           cycle
        end if

        ! C is zero-based, so shift bounds by one
        ishell_start = nucleus_index(inucl) + 1
        ishell_end   = nucleus_index(inucl) + nucleus_shell_num(inucl)

        do ishell=ishell_start, ishell_end

           ! Initialize Hessian components to zero
           shell_hessian(ishell, 1, 1, ipoint) = 0.d0  ! d²/dx²
           shell_hessian(ishell, 1, 2, ipoint) = 0.d0  ! d²/dy²
           shell_hessian(ishell, 1, 3, ipoint) = 0.d0  ! d²/dz²
           shell_hessian(ishell, 2, 1, ipoint) = 0.d0  ! d²/dxdy
           shell_hessian(ishell, 2, 2, ipoint) = 0.d0  ! d²/dxdz
           shell_hessian(ishell, 2, 3, ipoint) = 0.d0  ! d²/dydz
           shell_hessian(ishell, 3, 1, ipoint) = 0.d0
           shell_hessian(ishell, 3, 2, ipoint) = 0.d0
           shell_hessian(ishell, 3, 3, ipoint) = 0.d0
           shell_hessian(ishell, 4, 1, ipoint) = 0.d0
           shell_hessian(ishell, 4, 2, ipoint) = 0.d0
           shell_hessian(ishell, 4, 3, ipoint) = 0.d0

           iprim_start = shell_prim_index(ishell) + 1
           iprim_end   = shell_prim_index(ishell) + shell_prim_num(ishell)

           ! Handle r = 0 case for Slater orbitals
           if (r < 1.d-15) then
              ! Hessian is zero at r=0 for Slater orbitals
              cycle
           end if

           r_inv = 1.d0 / r
           r_inv2 = r_inv * r_inv

           do iprim = iprim_start, iprim_end

              alpha = expo(iprim)
              ar = alpha * r
              if (ar > cutoff) then
                 cycle
              end if

              n_val = dble(r_power(ishell) + 1)

              ! Compute r^n * exp(-a*r)
              v = coef_normalized(iprim) * (r ** int(n_val)) * dexp(-ar)

              ! For Slater orbitals: f = r^n * exp(-a*r)
              ! Hessian components involve second derivatives

              ! Diagonal Hessian terms: d²f/dx²
              ! v * (n^2 x^2+n (r^2-2 \[Alpha] r x^2-2 x^2)+\[Alpha] r (-r^2+\[Alpha] r x^2+x^2))/r^4

              shell_hessian(ishell, 1, 1, ipoint) = shell_hessian(ishell, 1, 1, ipoint) + &
                   v * r_inv2 * r_inv2 * ( &
                   n_val*n_val*x*x + n_val*(r*r-2.d0*alpha*r*x*x-2.d0*x*x) + &
                   alpha*r*(-r*r+alpha*r*x*x+x*x))

              shell_hessian(ishell, 1, 2, ipoint) = shell_hessian(ishell, 1, 2, ipoint) + &
                   v * r_inv2 * r_inv2 * ( &
                   n_val*n_val*y*y + n_val*(r*r-2.d0*alpha*r*y*y-2.d0*y*y) + &
                   alpha*r*(-r*r+alpha*r*y*y+y*y))

              shell_hessian(ishell, 1, 3, ipoint) = shell_hessian(ishell, 1, 3, ipoint) + &
                   v * r_inv2 * r_inv2 * ( &
                   n_val*n_val*z*z + n_val*(r*r-2.d0*alpha*r*z*z-2.d0*z*z) + &
                   alpha*r*(-r*r+alpha*r*z*z+z*z))

              ! Off-diagonal Hessian terms: d²f/dxdy, etc.
              hess_off = -v * r_inv2 * r_inv2 * (n_val*n_val - 2.d0*n_val * (alpha*r+1.d0) + &
                   alpha*r*(alpha*r+1.d0))

              shell_hessian(ishell, 2, 1, ipoint) = shell_hessian(ishell, 2, 1, ipoint) + &
                hess_off * x * y

              shell_hessian(ishell, 2, 2, ipoint) = shell_hessian(ishell, 2, 2, ipoint) + &
                hess_off * x * z

              shell_hessian(ishell, 2, 3, ipoint) = shell_hessian(ishell, 2, 3, ipoint) + &
                hess_off * y * z

           end do

        end do
     end do

  end do

end function qmckl_compute_ao_basis_shell_slater_hessian

function qmckl_ao_power(context, n, X, LMAX, P, ldp) &
     bind(C) result(info)
  use qmckl_constants
  implicit none

  integer (qmckl_context), intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: n
  integer (c_int64_t) , intent(in)  , value :: ldp
  real    (c_double ) , intent(in)          :: X(n)
  integer (c_int32_t) , intent(in)          :: LMAX(n)
  real    (c_double ) , intent(out)         :: P(ldp,n)

  integer(qmckl_exit_code) :: info
  integer(c_int64_t)  :: i,k

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (n <= ldp) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  k = MAXVAL(LMAX)
  if (LDP < k) then
     info = QMCKL_INVALID_ARG_6
     return
  endif

  if (k <= 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  do i=1,n
     P(1,i) = X(i)
     do k=2,LMAX(i)
        P(k,i) = P(k-1,i) * X(i)
     end do
  end do

end function qmckl_ao_power

function qmckl_ao_polynomial_vgl_doc (context,    &
     X, R, lmax, n, L, ldl, VGL, ldv) &
     bind(C) result(info)
  use qmckl_constants
  implicit none

  integer (qmckl_context), intent(in)  , value :: context
  real    (c_double ) , intent(in)          :: X(3)
  real    (c_double ) , intent(in)          :: R(3)
  integer (c_int32_t) , intent(in)  , value :: lmax
  integer (c_int64_t) , intent(inout)        :: n
  integer (c_int64_t) , intent(in)  , value :: ldl
  integer (c_int64_t) , intent(in)  , value :: ldv
  integer (c_int32_t) , intent(out)         :: L(ldl,(lmax+1)*(lmax+2)*(lmax+3)/6)
  real    (c_double ) , intent(out)         :: VGL(ldv,(lmax+1)*(lmax+2)*(lmax+3)/6)

  integer(qmckl_exit_code) :: info

  integer*8         :: i,j
  integer           :: a,b,c,d
  double precision  :: Y(3)
  double precision  :: pows(-2:lmax,3)
  double precision  :: xy, yz, xz
  double precision  :: da, db, dc, dd

  info = 0

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (lmax < 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  if (ldl < 3) then
     info = QMCKL_INVALID_ARG_7
     return
  endif

  if (ldv < 5) then
     info = QMCKL_INVALID_ARG_9
     return
  endif


  ! The shift below is such that polynomials will not make the AO equal to zero at the nodes of the orbitals
  do i=1,3
     Y(i) = (X(i) - R(i)) + 1.d-20
  end do

  if (lmax == 0) then
     VGL(1,1) = 1.d0
     VGL(2,1) = 0.d0
     VGL(3,1) = 0.d0
     VGL(4,1) = 0.d0
     VGL(5,1) = 0.d0
     l(1,1) = 0
     l(2,1) = 0
     l(3,1) = 0
     n=1
  else if (lmax > 0) then
     pows(-2:0,1:3) = 1.d0
     do i=1,lmax
        pows(i,1) = pows(i-1,1) * Y(1)
        pows(i,2) = pows(i-1,2) * Y(2)
        pows(i,3) = pows(i-1,3) * Y(3)
     end do

     VGL(1:5,1:4) = 0.d0
     VGL(1,1) = 1.d0
     VGL(1,2) = pows(1,1)
     VGL(2,2) = 1.d0
     VGL(1,3) = pows(1,2)
     VGL(3,3) = 1.d0
     VGL(1,4) = pows(1,3)
     VGL(4,4) = 1.d0

     l  (1:3,1:4) = 0
     l  (1,2) = 1
     l  (2,3) = 1
     l  (3,4) = 1

     n=4
  endif

  ! l>=2
  dd = 2.d0
  do d=2,lmax
     da = dd
     do a=d,0,-1
        db = dd-da
        do b=d-a,0,-1
           c  = d  - a  - b
           dc = dd - da - db
           n = n+1

           l(1,n) = a
           l(2,n) = b
           l(3,n) = c

           xy = pows(a,1) * pows(b,2)
           yz = pows(b,2) * pows(c,3)
           xz = pows(a,1) * pows(c,3)

           VGL(1,n) = xy * pows(c,3)

           xy = dc * xy
           xz = db * xz
           yz = da * yz

           VGL(2,n) = pows(a-1,1) * yz
           VGL(3,n) = pows(b-1,2) * xz
           VGL(4,n) = pows(c-1,3) * xy

           VGL(5,n) = &
                (da-1.d0) * pows(a-2,1) * yz + &
                (db-1.d0) * pows(b-2,2) * xz + &
                (dc-1.d0) * pows(c-2,3) * xy

           db = db - 1.d0
        end do
        da = da - 1.d0
     end do
     dd = dd + 1.d0
  end do

  info = QMCKL_SUCCESS

end function qmckl_ao_polynomial_vgl_doc

function qmckl_ao_polynomial_transp_vgl_doc (context,    &
     X, R, lmax, n, L, ldl, VGL, ldv) &
     bind(C) result(info)

  use qmckl_constants
  implicit none

  integer (qmckl_context), intent(in)  , value :: context
  real    (c_double ) , intent(in)          :: X(3)
  real    (c_double ) , intent(in)          :: R(3)
  integer (c_int32_t) , intent(in)  , value :: lmax
  integer (c_int64_t) , intent(inout)        :: n
  integer (c_int64_t) , intent(in)  , value :: ldl
  integer (c_int64_t) , intent(in)  , value :: ldv
  integer (c_int32_t) , intent(out)         :: L(ldl,(lmax+1)*(lmax+2)*(lmax+3)/6)
  real    (c_double ) , intent(out)         :: VGL(ldv,5)

  integer(qmckl_exit_code) :: info

  integer*8         :: i,j
  integer           :: a,b,c,d
  real*8            :: Y(3)
  real*8            :: pows(-2:21,3) ! lmax < 22
  double precision  :: xy, yz, xz
  double precision  :: da, db, dc, dd

  info = 0

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (lmax < 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  if (ldl < 3) then
     info = QMCKL_INVALID_ARG_7
     return
  endif

  if (ldv < (lmax+1)*(lmax+2)*(lmax+3)/6) then
     info = QMCKL_INVALID_ARG_9
     return
  endif


  if (lmax > 0) then

     do i=1,3
        Y(i) = X(i) - R(i)
     end do
     pows(-2:0,1:3) = 1.d0
     do i=1,lmax
        pows(i,1) = pows(i-1,1) * Y(1)
        pows(i,2) = pows(i-1,2) * Y(2)
        pows(i,3) = pows(i-1,3) * Y(3)
     end do

     l  (1:3,1:4) = 0
     VGL(1:4,1:5) = 0.d0

     VGL(1  ,1  ) = 1.d0

     l  (1,2) = 1
     VGL(2,1) = Y(1)
     VGL(2,2) = 1.d0

     l  (2,3) = 1
     VGL(3,1) = Y(2)
     VGL(3,3) = 1.d0

     l  (3,4) = 1
     VGL(4,1) = Y(3)
     VGL(4,4) = 1.d0

     n=4
  else
     VGL(1,1) = 1.d0
     VGL(1,2:5) = 0.d0
     l(1:3,1) = 0
     n=1
     return
  endif

  ! l>=2
  dd = 2.d0
  do d=2,lmax
     da = dd
     do a=d,0,-1
        db = dd-da
        do b=d-a,0,-1
           c  = d  - a  - b
           dc = dd - da - db
           n = n+1

           xy = pows(a,1) * pows(b,2)
           yz = pows(b,2) * pows(c,3)
           xz = pows(a,1) * pows(c,3)

           l(1,n) = a
           l(2,n) = b
           l(3,n) = c

           VGL(n,1) = xy * pows(c,3)

           xy = dc * xy
           xz = db * xz
           yz = da * yz

           VGL(n,2) = pows(a-1,1) * yz
           VGL(n,3) = pows(b-1,2) * xz
           VGL(n,4) = pows(c-1,3) * xy

           VGL(n,5) = &
                (da-1.d0) * pows(a-2,1) * yz + &
                (db-1.d0) * pows(b-2,2) * xz + &
                (dc-1.d0) * pows(c-2,3) * xy

           db = db - 1.d0
        end do
        da = da - 1.d0
     end do
     dd = dd + 1.d0
  end do

  info = QMCKL_SUCCESS

end function qmckl_ao_polynomial_transp_vgl_doc

function qmckl_compute_ao_polynomial_hessian_doc (context,    &
     X, R, lmax, n, hessian) &
     bind(C) result(info)
  use qmckl_constants
  implicit none

  integer (qmckl_context), intent(in)  , value :: context
  real    (c_double ) , intent(in)          :: X(3)
  real    (c_double ) , intent(in)          :: R(3)
  integer (c_int32_t) , intent(in)  , value :: lmax
  integer (c_int64_t) , intent(inout)        :: n
  real    (c_double ) , intent(out)         :: hessian(3,4,(lmax+1)*(lmax+2)*(lmax+3)/6)

  integer(qmckl_exit_code) :: info

  integer           :: i,j,m,k
  integer           :: a,b,c,d
  double precision  :: Y(3)
  double precision  :: pows(-2:lmax,3)
  double precision  :: xy, yz, xz
  double precision  :: da, db, dc, dd

  info = 0

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (lmax < 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif


  ! The shift below is such that polynomials will not make the AO equal to zero at the nodes of the orbitals
  do i=1,3
     Y(i) = (X(i) - R(i)) + 1.d-20
  end do

  if (lmax == 0) then
     do k = 1, 3
        do m = 1, 3
           hessian(k,m,1) = 0.d0
        end do
      end do
     n=1
  else if (lmax > 0) then
     pows(-2:0,1:3) = 1.d0
     do i=1,lmax
        pows(i,1) = pows(i-1,1) * Y(1)
        pows(i,2) = pows(i-1,2) * Y(2)
        pows(i,3) = pows(i-1,3) * Y(3)
     end do

     hessian(1:3,1:4,1:4) = 0.d0

     n=4
  endif

  ! l>=2
  dd = 2.d0
  do d=2,lmax
     da = dd
     do a=d,0,-1
        db = dd-da
        do b=d-a,0,-1
           c  = d  - a  - b
           dc = dd - da - db
           n = n+1

           xy = pows(a,1) * pows(b,2)
           yz = pows(b,2) * pows(c,3)
           xz = pows(a,1) * pows(c,3)

           xy = dc * xy
           xz = db * xz
           yz = da * yz

          hessian(1,1,n) = da * (da-1.d0) * pows(a-2,1) * pows(b,2) * pows(c,3)
          hessian(2,2,n) = db * (db-1.d0) * pows(a,1) * pows(b-2,2) * pows(c,3)
          hessian(3,3,n) = dc * (dc-1.d0) * pows(a,1) * pows(b,2) * pows(c-2,3)

          hessian(1,2,n) = da * db * pows(a-1,1) * pows(b-1,2) * pows(c,3)
          hessian(2,1,n) = da * db * pows(a-1,1) * pows(b-1,2) * pows(c,3)

          hessian(1,3,n) = da * dc * pows(a-1,1) * pows(b,2) * pows(c-1,3)
          hessian(3,1,n) = da * dc * pows(a-1,1) * pows(b,2) * pows(c-1,3)

          hessian(2,3,n) = db * dc * pows(a,1) * pows(b-1,2) * pows(c-1,3)
          hessian(3,2,n) = db * dc * pows(a,1) * pows(b-1,2) * pows(c-1,3)

          hessian(1,4,n) = (da * db * (db-1.d0) * pows(a-1,1) * pows(b-2,2) * pows(c,3)) + &
                           (da * dc * (dc-1.d0) * pows(a-1,1) * pows(b,2) * pows(c-2,3))

          hessian(2,4,n) = (db * da * (da-1.d0) * pows(a-2,1) * pows(b-1,2) * pows(c,3)) + &
                           (db * dc * (dc-1.d0) * pows(a,1) * pows(b-1,2) * pows(c-2,3))

          hessian(3,4,n) = (dc * da * (da-1.d0) * pows(a-2,1) * pows(b,2) * pows(c-1,3)) + &
                           (dc * db * (db-1.d0) * pows(a,1) * pows(b-2,2) * pows(c-1,3))

          if (da > 2) then
            hessian(1,4,n) = hessian(1,4,n) + (da * (da-1.d0) * (da-2.d0) * &
                 pows(a-3,1) * pows(b,2) * pows(c,3))
          end if
          if (db > 2) then
            hessian(2,4,n) = hessian(2,4,n) + (db * (db-1.d0) * (db-2.d0) * &
                 pows(a,1) * pows(b-3,2) * pows(c,3))
          end if
          if (dc > 2) then
            hessian(3,4,n) = hessian(3,4,n) + (dc * (dc-1.d0) * (dc-2.d0) * &
                 pows(a,1) * pows(b,2) * pows(c-3,3))
          end if

           db = db - 1.d0
        end do
        da = da - 1.d0
     end do
     dd = dd + 1.d0
  end do

  info = QMCKL_SUCCESS

end function qmckl_compute_ao_polynomial_hessian_doc

function qmckl_compute_nucleus_range_gaussian(context, &
     ao_num, shell_num, prim_num, nucl_num, &
     nucl_coord, nucleus_index, nucleus_shell_num, &
     nucleus_max_ang_mom, shell_prim_index, shell_prim_num, r_power, &
     ao_factor, expo, coef_normalized, nucleus_range) &
     bind(C) result(info)
  use qmckl_constants
  use qmckl, only: qmckl_ao_polynomial_vgl, qmckl_get_numprec_epsilon
  implicit none
  integer (qmckl_context), intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: ao_num
  integer (c_int64_t) , intent(in)  , value :: shell_num
  integer (c_int64_t) , intent(in)  , value :: prim_num
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  real    (c_double ) , intent(in)          :: nucl_coord(nucl_num,3)
  integer (c_int64_t) , intent(in)          :: nucleus_index(nucl_num)
  integer (c_int64_t) , intent(in)          :: nucleus_shell_num(nucl_num)
  integer (c_int32_t) , intent(in)          :: nucleus_max_ang_mom(nucl_num)
  integer (c_int64_t) , intent(in)          :: shell_prim_index(shell_num)
  integer (c_int64_t) , intent(in)          :: shell_prim_num(shell_num)
  integer (c_int32_t) , intent(in)          :: r_power(shell_num)
  real    (c_double ) , intent(in)          :: ao_factor(ao_num)
  real    (c_double ) , intent(in)          :: expo(prim_num)
  real    (c_double ) , intent(in)          :: coef_normalized(prim_num)
  real    (c_double ) , intent(out)         :: nucleus_range(nucl_num,53)

  integer(qmckl_exit_code) :: info

  double precision  :: e_coord(3), n_coord(3)
  integer*8         :: n_poly, iprim, iprim_start, iprim_end
  integer           :: l, il, k, iprecision
  integer*8         :: ipoint, inucl, ishell
  integer*8         :: ishell_start, ishell_end
  integer           :: lstart(0:20)
  double precision  :: shell_vgl(3)
  double precision  :: x, y, z, r2, ar2, two_a
  double precision  :: vmax, cutoff, v

  double precision, allocatable  :: poly_vgl(:,:), ao_vgl(:,:)
  integer         , allocatable  :: powers(:,:), ao_index(:)

  allocate(poly_vgl(5,ao_num), powers(3,ao_num), ao_index(ao_num))

  ! Pre-computed data
  do l=0,20
     lstart(l) = l*(l+1)*(l+2)/6 +1
  end do

  k=1
  do inucl=1,nucl_num
     ishell_start = nucleus_index(inucl) + 1
     ishell_end   = nucleus_index(inucl) + nucleus_shell_num(inucl)
     do ishell = ishell_start, ishell_end
        l = r_power(ishell)
        ao_index(ishell) = k
        k = k + lstart(l+1) - lstart(l)
     end do
  end do
  info = QMCKL_SUCCESS

  n_coord(1) = 0.d0
  n_coord(2) = 0.d0
  n_coord(3) = 0.d0
  e_coord(2) = 0.d0
  e_coord(3) = 0.d0

  nucleus_range = 50.d0

  do inucl=1,nucl_num

     x = 50.d0
     do iprecision = 53,2,-1

        cutoff = 2.d0** (1-iprecision)

        vmax = 0.d0
        do while ( (vmax < cutoff) .and. (x > 0.d0) )
           x = x - .1d0
           vmax = 0.d0
           e_coord(1) = x
           r2 = x*x

           ! Compute polynomials
           info = qmckl_ao_polynomial_vgl(context, e_coord, n_coord, &
                nucleus_max_ang_mom(inucl), n_poly, powers, 3_8, &
                poly_vgl, 5_8)

           ! Loop over shells
           ishell_start = nucleus_index(inucl) + 1
           ishell_end   = nucleus_index(inucl) + nucleus_shell_num(inucl)
           do ishell = ishell_start, ishell_end
              shell_vgl(1) = 0.d0
              shell_vgl(2) = 0.d0
              shell_vgl(3) = 0.d0

              iprim_start = shell_prim_index(ishell) + 1
              iprim_end   = shell_prim_index(ishell) + shell_prim_num(ishell)

              do iprim = iprim_start, iprim_end

                 ar2 = expo(iprim)*r2

                 v = coef_normalized(iprim) * dexp(-ar2)
                 two_a = -2.d0 * expo(iprim) * v

                 shell_vgl(1) = shell_vgl(1) + v
                 shell_vgl(2) = shell_vgl(2) + two_a * x
                 shell_vgl(3) = shell_vgl(3) + two_a * (3.d0 - 2.d0*ar2)

              end do

              k = ao_index(ishell)
              l = r_power(ishell)
              do il = lstart(l), lstart(l+1)-1
                 vmax = max(vmax, poly_vgl(1,il) * shell_vgl(1) * ao_factor(k))
                 vmax = max(vmax, ( poly_vgl(2,il) * shell_vgl(1) + poly_vgl(1,il) * shell_vgl(2) ) * ao_factor(k))
                 vmax = max(vmax, ( poly_vgl(5,il) * shell_vgl(3) + poly_vgl(1,il) * shell_vgl(3) + &
                      2.d0 * (poly_vgl(2,il) * shell_vgl(2) ) ) * ao_factor(k) )

                 k = k+1
              end do
           end do
        end do
        nucleus_range(inucl,iprecision) = x
     end do
  end do



  deallocate(poly_vgl, powers)
end function qmckl_compute_nucleus_range_gaussian

function qmckl_compute_nucleus_range_slater(context, &
     ao_num, shell_num, prim_num, nucl_num, &
     nucl_coord, nucleus_index, nucleus_shell_num, &
     nucleus_max_ang_mom, shell_prim_index, shell_prim_num, r_power, &
     ao_factor, expo, coef_normalized, nucleus_range) &
     bind(C) result(info)
  use qmckl_constants
  use qmckl, only: qmckl_ao_polynomial_vgl
  implicit none
  integer (qmckl_context), intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: ao_num
  integer (c_int64_t) , intent(in)  , value :: shell_num
  integer (c_int64_t) , intent(in)  , value :: prim_num
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  real    (c_double ) , intent(in)          :: nucl_coord(nucl_num,3)
  integer (c_int64_t) , intent(in)          :: nucleus_index(nucl_num)
  integer (c_int64_t) , intent(in)          :: nucleus_shell_num(nucl_num)
  integer (c_int32_t) , intent(in)          :: nucleus_max_ang_mom(nucl_num)
  integer (c_int64_t) , intent(in)          :: shell_prim_index(shell_num)
  integer (c_int64_t) , intent(in)          :: shell_prim_num(shell_num)
  integer (c_int32_t) , intent(in)          :: r_power(shell_num)
  real    (c_double ) , intent(in)          :: ao_factor(ao_num)
  real    (c_double ) , intent(in)          :: expo(prim_num)
  real    (c_double ) , intent(in)          :: coef_normalized(prim_num)
  real    (c_double ) , intent(out)         :: nucleus_range(nucl_num,53)

  integer(qmckl_exit_code) :: info

  double precision  :: e_coord(3), n_coord(3)
  integer*8         :: n_poly, iprim, iprim_start, iprim_end
  integer           :: l, il, k, iprecision
  integer*8         :: ipoint, inucl, ishell
  integer*8         :: ishell_start, ishell_end
  integer           :: lstart(0:20)
  double precision  :: shell_vgl(3)
  double precision  :: x, y, z, r, ar, two_a
  double precision  :: vmax, cutoff, v

  double precision, allocatable  :: poly_vgl(:,:), ao_vgl(:,:)
  integer         , allocatable  :: powers(:,:), ao_index(:)
  double precision  :: r_inv, r_inv_2, x_over_r, alpha, rn, grad_coef, lapl_coef, n_val

  allocate(poly_vgl(5,ao_num), powers(3,ao_num), ao_index(ao_num))

  ! Pre-computed data
  do l=0,20
     lstart(l) = l*(l+1)*(l+2)/6 +1
  end do

  k=1
  do inucl=1,nucl_num
     ishell_start = nucleus_index(inucl) + 1
     ishell_end   = nucleus_index(inucl) + nucleus_shell_num(inucl)
     do ishell = ishell_start, ishell_end
        l = r_power(ishell)
        ao_index(ishell) = k
        k = k + lstart(l+1) - lstart(l)
     end do
  end do
  info = QMCKL_SUCCESS

  n_coord(1) = 0.d0
  n_coord(2) = 0.d0
  n_coord(3) = 0.d0
  e_coord(2) = 0.d0
  e_coord(3) = 0.d0

  nucleus_range = 50.d0

  do inucl=1,nucl_num

     x = 50.d0
     do iprecision = 53,2,-1

        cutoff = 2.d0** (1-iprecision)

        vmax = 0.d0
        do while ( (vmax < cutoff) .and. (x > 0.d0) )
           x = x - .1d0
           vmax = 0.d0
           e_coord(1) = x
           r = x

           r_inv = 1.d0 / r
           r_inv_2 = r_inv * r_inv
           x_over_r = e_coord(1) * r_inv

           ! Compute polynomials
           info = qmckl_ao_polynomial_vgl(context, e_coord, n_coord, &
                nucleus_max_ang_mom(inucl), n_poly, powers, 3_8, &
                poly_vgl, 5_8)

           ! Loop over shells
           ishell_start = nucleus_index(inucl) + 1
           ishell_end   = nucleus_index(inucl) + nucleus_shell_num(inucl)
           do ishell = ishell_start, ishell_end
              shell_vgl(1) = 0.d0
              shell_vgl(2) = 0.d0
              shell_vgl(3) = 0.d0

              iprim_start = shell_prim_index(ishell) + 1
              iprim_end   = shell_prim_index(ishell) + shell_prim_num(ishell)

              do iprim = iprim_start, iprim_end

                 alpha = expo(iprim)
                 ar = alpha*r
                 n_val = dble(r_power(ishell))
                 rn = r**(r_power(ishell))

                 v = coef_normalized(iprim) * rn*dexp(-ar)

                 grad_coef = (n_val * r_inv - alpha) * v
                 lapl_coef = n_val * (n_val - 1.d0) * r_inv_2 - 2.d0 * n_val * alpha * r_inv + alpha * alpha
                 shell_vgl(1) = shell_vgl(1) + v
                 shell_vgl(2) = shell_vgl(2) + grad_coef * x_over_r
                 shell_vgl(3) = shell_vgl(3) + lapl_coef * v

              end do

              k = ao_index(ishell)
              l = r_power(ishell)
              do il = lstart(l), lstart(l+1)-1
                 vmax = max(vmax, poly_vgl(1,il) * shell_vgl(1) * ao_factor(k))
                 vmax = max(vmax, ( poly_vgl(2,il) * shell_vgl(1) + poly_vgl(1,il) * shell_vgl(2) ) * ao_factor(k))
                 vmax = max(vmax, ( poly_vgl(5,il) * shell_vgl(3) + poly_vgl(1,il) * shell_vgl(3) + &
                      2.d0 * (poly_vgl(2,il) * shell_vgl(2) ) ) * ao_factor(k) )

                 k = k+1
              end do
           end do
        end do
        nucleus_range(inucl,iprecision) = x
     end do
  end do



  deallocate(poly_vgl, powers)
end function qmckl_compute_nucleus_range_slater

! Unoptimized version
!      #+NAME: qmckl_ao_value_args_doc
!     | Variable              | Type                              | In/Out | Description                                              |
!     |-----------------------+-----------------------------------+--------+----------------------------------------------------------|
!     | ~context~             | ~qmckl_context~                   | in     | Global state                                             |
!     | ~ao_num~              | ~int64_t~                         | in     | Number of AOs                                            |
!     | ~shell_num~           | ~int64_t~                         | in     | Number of shells                                         |
!     | ~point_num~           | ~int64_t~                         | in     | Number of points                                         |
!     | ~nucl_num~            | ~int64_t~                         | in     | Number of nuclei                                         |
!     | ~coord~               | ~double[3][point_num]~            | in     | Coordinates                                              |
!     | ~nucl_coord~          | ~double[3][nucl_num]~             | in     | Nuclear  coordinates                                     |
!     | ~nucleus_index~       | ~int64_t[nucl_num]~               | in     | Index of the 1st shell of each nucleus                   |
!     | ~nucleus_shell_num~   | ~int64_t[nucl_num]~               | in     | Number of shells per nucleus                             |
!     | ~nucleus_range~       | ~double[nucl_num]~                | in     | Range beyond which all is zero                           |
!     | ~nucleus_max_ang_mom~ | ~int32_t[nucl_num]~               | in     | Maximum angular momentum per nucleus                     |
!     | ~shell_ang_mom~       | ~int32_t[shell_num]~              | in     | Angular momentum of each shell                           |
!     | ~ao_factor~           | ~double[ao_num]~                  | in     | Normalization factor of the AOs                          |
!     | ~shell_vgl~           | ~double[point_num][5][shell_num]~ | in     | Value, gradients and Laplacian of the shells             |
!     | ~ao_value~            | ~double[point_num][ao_num]~       | out    | Values of the AOs                                        |


function qmckl_compute_ao_value_doc(context, &
     ao_num, shell_num, point_num, nucl_num, &
     coord, nucl_coord, nucleus_index, nucleus_shell_num, &
     nucleus_range, nucleus_max_ang_mom, shell_ang_mom, &
     ao_factor, shell_vgl, ao_value) &
     bind(C) result(info)
  use qmckl_constants
  use qmckl, only : qmckl_ao_polynomial_vgl, qmckl_get_numprec_epsilon
  implicit none
  integer (qmckl_context), intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: ao_num
  integer (c_int64_t) , intent(in)  , value :: shell_num
  integer (c_int64_t) , intent(in)  , value :: point_num
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  real    (c_double ) , intent(in)          :: coord(point_num,3)
  real    (c_double ) , intent(in)          :: nucl_coord(nucl_num,3)
  integer (c_int64_t) , intent(in)          :: nucleus_index(nucl_num)
  integer (c_int64_t) , intent(in)          :: nucleus_shell_num(nucl_num)
  real    (c_double ) , intent(in)          :: nucleus_range(nucl_num)
  integer (c_int32_t) , intent(in)          :: nucleus_max_ang_mom(nucl_num)
  integer (c_int32_t) , intent(in)          :: shell_ang_mom(shell_num)
  real    (c_double ) , intent(in)          :: ao_factor(ao_num)
  real    (c_double ) , intent(in)          :: shell_vgl(shell_num,5,point_num)
  real    (c_double ) , intent(out)         :: ao_value(ao_num,point_num)

  integer(qmckl_exit_code) :: info

  double precision  :: e_coord(3), n_coord(3)
  integer*8         :: n_poly
  integer           :: l, il, k
  integer*8         :: ipoint, inucl, ishell
  integer*8         :: ishell_start, ishell_end
  integer           :: lstart(0:20)
  double precision  :: x, y, z, r2
  double precision  :: cutoff

  double precision, allocatable  :: poly_vgl(:,:)
  integer         , allocatable  :: powers(:,:), ao_index(:)

  allocate(poly_vgl(5,ao_num), powers(3,ao_num), ao_index(ao_num))

  ! Pre-computed data
  do l=0,20
     lstart(l) = l*(l+1)*(l+2)/6 +1
  end do

  k=1
  do inucl=1,nucl_num
     ishell_start = nucleus_index(inucl) + 1
     ishell_end   = nucleus_index(inucl) + nucleus_shell_num(inucl)
     do ishell = ishell_start, ishell_end
        l = shell_ang_mom(ishell)
        ao_index(ishell) = k
        k = k + lstart(l+1) - lstart(l)
     end do
  end do
  info = QMCKL_SUCCESS

  ! Don't compute polynomials when the radial part is zero.
  cutoff = -dlog(qmckl_get_numprec_epsilon(context))

  do ipoint = 1, point_num
     e_coord(1) = coord(ipoint,1)
     e_coord(2) = coord(ipoint,2)
     e_coord(3) = coord(ipoint,3)
     ao_value(:,ipoint) = 0.d0
     do inucl=1,nucl_num
        n_coord(1) = nucl_coord(inucl,1)
        n_coord(2) = nucl_coord(inucl,2)
        n_coord(3) = nucl_coord(inucl,3)

        ! Test if the point is in the range of the nucleus
        x = e_coord(1) - n_coord(1)
        y = e_coord(2) - n_coord(2)
        z = e_coord(3) - n_coord(3)

        r2 = x*x + y*y + z*z

        if (r2 > cutoff*nucleus_range(inucl)) then
           cycle
        end if

        ! Compute polynomials
        info = qmckl_ao_polynomial_vgl(context, e_coord, n_coord, &
             nucleus_max_ang_mom(inucl), n_poly, powers, 3_8, &
             poly_vgl, 5_8)

        ! Loop over shells
        ishell_start = nucleus_index(inucl) + 1
        ishell_end   = nucleus_index(inucl) + nucleus_shell_num(inucl)
        do ishell = ishell_start, ishell_end
           k = ao_index(ishell)
           l = shell_ang_mom(ishell)
           do il = lstart(l), lstart(l+1)-1
              ! Value
              ao_value(k,ipoint) = &
                   poly_vgl(1,il) * shell_vgl(ishell,1,ipoint) * ao_factor(k)
              k = k+1
           end do
        end do
     end do
  end do

  deallocate(poly_vgl, powers)
end function qmckl_compute_ao_value_doc

! Reference version
!      #+NAME: qmckl_ao_vgl_args_doc
!     | Variable              | Type                              | In/Out | Description                                  |
!     |-----------------------+-----------------------------------+--------+----------------------------------------------|
!     | ~context~             | ~qmckl_context~                   | in     | Global state                                 |
!     | ~ao_num~              | ~int64_t~                         | in     | Number of AOs                                |
!     | ~shell_num~           | ~int64_t~                         | in     | Number of shells                             |
!     | ~point_num~           | ~int64_t~                         | in     | Number of points                             |
!     | ~nucl_num~            | ~int64_t~                         | in     | Number of nuclei                             |
!     | ~coord~               | ~double[3][point_num]~            | in     | Coordinates                                  |
!     | ~nucl_coord~          | ~double[3][nucl_num]~             | in     | Nuclear  coordinates                         |
!     | ~nucleus_index~       | ~int64_t[nucl_num]~               | in     | Index of the 1st shell of each nucleus       |
!     | ~nucleus_shell_num~   | ~int64_t[nucl_num]~               | in     | Number of shells per nucleus                 |
!     | ~nucleus_range~       | ~double[nucl_num]~                | in     | Range beyond which all is zero               |
!     | ~nucleus_max_ang_mom~ | ~int32_t[nucl_num]~               | in     | Maximum angular momentum per nucleus         |
!     | ~shell_ang_mom~       | ~int32_t[shell_num]~              | in     | Angular momentum of each shell               |
!     | ~ao_factor~           | ~double[ao_num]~                  | in     | Normalization factor of the AOs              |
!     | ~shell_vgl~           | ~double[point_num][5][shell_num]~ | in     | Value, gradients and Laplacian of the shells |
!     | ~ao_vgl~              | ~double[point_num][5][ao_num]~    | out    | Value, gradients and Laplacian of the AOs    |


function qmckl_compute_ao_vgl_doc(context, &
     ao_num, shell_num, point_num, nucl_num, &
     coord, nucl_coord, nucleus_index, nucleus_shell_num, &
     nucleus_range, nucleus_max_ang_mom, shell_ang_mom, &
     ao_factor, shell_vgl, ao_vgl) &
     bind(C) result(info)
  use qmckl_constants
  use qmckl, only : qmckl_ao_polynomial_vgl, qmckl_get_numprec_epsilon
  implicit none
  integer (qmckl_context), intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: ao_num
  integer (c_int64_t) , intent(in)  , value :: shell_num
  integer (c_int64_t) , intent(in)  , value :: point_num
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  real    (c_double ) , intent(in)          :: coord(point_num,3)
  real    (c_double ) , intent(in)          :: nucl_coord(nucl_num,3)
  integer (c_int64_t) , intent(in)          :: nucleus_index(nucl_num)
  integer (c_int64_t) , intent(in)          :: nucleus_shell_num(nucl_num)
  real    (c_double ) , intent(in)          :: nucleus_range(nucl_num)
  integer (c_int32_t) , intent(in)          :: nucleus_max_ang_mom(nucl_num)
  integer (c_int32_t) , intent(in)          :: shell_ang_mom(shell_num)
  real    (c_double ) , intent(in)          :: ao_factor(ao_num)
  real    (c_double ) , intent(in)          :: shell_vgl(shell_num,5,point_num)
  real    (c_double ) , intent(out)         :: ao_vgl(ao_num,5,point_num)
  integer(qmckl_exit_code) :: info

  double precision  :: e_coord(3), n_coord(3)
  integer*8         :: n_poly
  integer           :: l, il, k
  integer*8         :: ipoint, inucl, ishell
  integer*8         :: ishell_start, ishell_end
  integer           :: lstart(0:20)
  double precision  :: x, y, z, r2
  double precision  :: cutoff

  double precision, allocatable  :: poly_vgl(:,:)
  integer         , allocatable  :: powers(:,:), ao_index(:)

  allocate(poly_vgl(5,ao_num), powers(3,ao_num), ao_index(ao_num))

  ! Pre-computed data
  do l=0,20
     lstart(l) = l*(l+1)*(l+2)/6 +1
  end do

  k=1
  do inucl=1,nucl_num
     ishell_start = nucleus_index(inucl) + 1
     ishell_end   = nucleus_index(inucl) + nucleus_shell_num(inucl)
     do ishell = ishell_start, ishell_end
        l = shell_ang_mom(ishell)
        ao_index(ishell) = k
        k = k + lstart(l+1) - lstart(l)
     end do
  end do
  info = QMCKL_SUCCESS

  ! Don't compute polynomials when the radial part is zero.
  cutoff = -dlog(qmckl_get_numprec_epsilon(context))

  do ipoint = 1, point_num
     e_coord(1) = coord(ipoint,1)
     e_coord(2) = coord(ipoint,2)
     e_coord(3) = coord(ipoint,3)
     ao_vgl(:,:,ipoint) = 0.d0
     do inucl=1,nucl_num
        n_coord(1) = nucl_coord(inucl,1)
        n_coord(2) = nucl_coord(inucl,2)
        n_coord(3) = nucl_coord(inucl,3)

        ! Test if the point is in the range of the nucleus
        x = e_coord(1) - n_coord(1)
        y = e_coord(2) - n_coord(2)
        z = e_coord(3) - n_coord(3)

        r2 = x*x + y*y + z*z

        if (r2 > cutoff*nucleus_range(inucl)) then
           cycle
        end if

        ! Compute polynomials
        info = qmckl_ao_polynomial_vgl(context, e_coord, n_coord, &
             nucleus_max_ang_mom(inucl), n_poly, powers, 3_8, &
             poly_vgl, 5_8)

        ! Loop over shells
        ishell_start = nucleus_index(inucl) + 1
        ishell_end   = nucleus_index(inucl) + nucleus_shell_num(inucl)
        do ishell = ishell_start, ishell_end
           k = ao_index(ishell)
           l = shell_ang_mom(ishell)
           do il = lstart(l), lstart(l+1)-1
              ! Value
              ao_vgl(k,1,ipoint) = &
                   poly_vgl(1,il) * shell_vgl(ishell,1,ipoint) * ao_factor(k)

              ! Grad_x
              ao_vgl(k,2,ipoint) = ( &
                   poly_vgl(2,il) * shell_vgl(ishell,1,ipoint) + &
                   poly_vgl(1,il) * shell_vgl(ishell,2,ipoint) &
                   ) * ao_factor(k)

              ! Grad_y
              ao_vgl(k,3,ipoint) = ( &
                   poly_vgl(3,il) * shell_vgl(ishell,1,ipoint) + &
                   poly_vgl(1,il) * shell_vgl(ishell,3,ipoint) &
                   ) * ao_factor(k)

              ! Grad_z
              ao_vgl(k,4,ipoint) = ( &
                   poly_vgl(4,il) * shell_vgl(ishell,1,ipoint) + &
                   poly_vgl(1,il) * shell_vgl(ishell,4,ipoint) &
                   ) * ao_factor(k)

              ! Lapl_z
              ao_vgl(k,5,ipoint) = ( &
                   poly_vgl(5,il) * shell_vgl(ishell,1,ipoint) + &
                   poly_vgl(1,il) * shell_vgl(ishell,5,ipoint) + &
                   2.d0 * ( &
                   poly_vgl(2,il) * shell_vgl(ishell,2,ipoint) + &
                   poly_vgl(3,il) * shell_vgl(ishell,3,ipoint) + &
                   poly_vgl(4,il) * shell_vgl(ishell,4,ipoint) ) &
                   ) * ao_factor(k)

              k = k+1
           end do
        end do
     end do
  end do

  deallocate(poly_vgl, powers)
end function qmckl_compute_ao_vgl_doc

! Compute

!    :PROPERTIES:
!    :Name:     qmckl_compute_ao_hessian
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!      #+NAME: qmckl_ao_hessian_args_doc
!     | Variable              | Type                                        | In/Out | Description                                  |
!     |-----------------------+---------------------------------------------+--------+----------------------------------------------|
!     | ~context~             | ~qmckl_context~                             | in     | Global state                                 |
!     | ~ao_num~              | ~int64_t~                                   | in     | Number of AOs                                |
!     | ~shell_num~           | ~int64_t~                                   | in     | Number of shells                             |
!     | ~point_num~           | ~int64_t~                                   | in     | Number of points                             |
!     | ~nucl_num~            | ~int64_t~                                   | in     | Number of nuclei                             |
!     | ~coord~               | ~double[3][point_num]~                      | in     | Coordinates                                  |
!     | ~nucl_coord~          | ~double[3][nucl_num]~                       | in     | Nuclear  coordinates                         |
!     | ~nucleus_index~       | ~int64_t[nucl_num]~                         | in     | Index of the 1st shell of each nucleus       |
!     | ~nucleus_shell_num~   | ~int64_t[nucl_num]~                         | in     | Number of shells per nucleus                 |
!     | ~nucleus_range~       | ~double[nucl_num]~                          | in     | Range beyond which all is zero               |
!     | ~nucleus_max_ang_mom~ | ~int32_t[nucl_num]~                         | in     | Maximum angular momentum per nucleus         |
!     | ~shell_ang_mom~       | ~int32_t[shell_num]~                        | in     | Angular momentum of each shell               |
!     | ~ao_factor~           | ~double[ao_num]~                            | in     | Normalization factor of the AOs              |
!     | ~shell_vgl~           | ~double[point_num][5][shell_num]~           | in     | Value, gradients and Laplacian of the shells |
!     | ~shell_hessian~       | ~double[point_num][3][4][shell_num]~        | in     | Hessian of the shells                        |
!     | ~ao_hessian~          | ~double[3][point_num][4][ao_num]~           | out    | Hessian of the AOs                           |
!     |-----------------------+---------------------------------------------+--------+----------------------------------------------|


function qmckl_compute_ao_hessian_doc(context, &
     ao_num, shell_num, point_num, nucl_num, &
     coord, nucl_coord, nucleus_index, nucleus_shell_num, &
     nucleus_range, nucleus_max_ang_mom, shell_ang_mom, &
     ao_factor, shell_vgl, shell_hessian, ao_hessian) &
     bind(C) result(info)
  use qmckl_constants
  use qmckl, only : qmckl_ao_polynomial_hessian, qmckl_ao_polynomial_vgl, qmckl_get_numprec_epsilon
  implicit none
  integer (qmckl_context), intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: ao_num
  integer (c_int64_t) , intent(in)  , value :: shell_num
  integer (c_int64_t) , intent(in)  , value :: point_num
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  real    (c_double ) , intent(in)          :: coord(point_num,3)
  real    (c_double ) , intent(in)          :: nucl_coord(nucl_num,3)
  integer (c_int64_t) , intent(in)          :: nucleus_index(nucl_num)
  integer (c_int64_t) , intent(in)          :: nucleus_shell_num(nucl_num)
  real    (c_double ) , intent(in)          :: nucleus_range(nucl_num)
  integer (c_int32_t) , intent(in)          :: nucleus_max_ang_mom(nucl_num)
  integer (c_int32_t) , intent(in)          :: shell_ang_mom(shell_num)
  real    (c_double ) , intent(in)          :: ao_factor(ao_num)
  real    (c_double ) , intent(in)          :: shell_vgl(shell_num,5,point_num)
  real    (c_double ) , intent(in)          :: shell_hessian(shell_num,4,3,point_num)
  real    (c_double ) , intent(out)         :: ao_hessian(ao_num,4,point_num,3)
  integer(qmckl_exit_code) :: info

  double precision  :: e_coord(3), n_coord(3)
  integer*8         :: n_poly
  integer           :: l, il, k, i, j
  integer*8         :: ipoint, inucl, ishell
  integer*8         :: ishell_start, ishell_end
  integer           :: lstart(0:20)
  double precision  :: x, y, z, r2
  double precision  :: cutoff

  double precision, allocatable  :: poly_hessian(:,:,:), poly_vgl(:,:)
  integer         , allocatable  :: powers(:,:), ao_index(:)

  allocate(poly_vgl(5,ao_num), poly_hessian(3,4,ao_num), powers(3,ao_num), ao_index(ao_num))

  ! Pre-computed data
  do l=0,20
     lstart(l) = l*(l+1)*(l+2)/6 +1
  end do

  k=1
  do inucl=1,nucl_num
     ishell_start = nucleus_index(inucl) + 1
     ishell_end   = nucleus_index(inucl) + nucleus_shell_num(inucl)
     do ishell = ishell_start, ishell_end
        l = shell_ang_mom(ishell)
        ao_index(ishell) = k
        k = k + lstart(l+1) - lstart(l)
     end do
  end do
  info = QMCKL_SUCCESS

  ! Don't compute polynomials when the radial part is zero.
  cutoff = -dlog(qmckl_get_numprec_epsilon(context))

  !ao_hessian = 0.d0

  do ipoint = 1, point_num
     e_coord(1) = coord(ipoint,1)
     e_coord(2) = coord(ipoint,2)
     e_coord(3) = coord(ipoint,3)
     ao_hessian(:,:,ipoint,:) = 0.d0
     do inucl=1,nucl_num
        n_coord(1) = nucl_coord(inucl,1)
        n_coord(2) = nucl_coord(inucl,2)
        n_coord(3) = nucl_coord(inucl,3)

        ! Test if the point is in the range of the nucleus
        x = e_coord(1) - n_coord(1)
        y = e_coord(2) - n_coord(2)
        z = e_coord(3) - n_coord(3)

        r2 = x*x + y*y + z*z
        if (r2 > cutoff*nucleus_range(inucl)) then
           !ao_hessian(:,:,ipoint,:) = 0.d0
           cycle
        end if

        ! Compute polynomials
        info = qmckl_ao_polynomial_hessian(context, e_coord, n_coord, &
             nucleus_max_ang_mom(inucl), n_poly, poly_hessian)

        info = qmckl_ao_polynomial_vgl(context, e_coord, n_coord, &
             nucleus_max_ang_mom(inucl), n_poly, powers, 3_8, &
             poly_vgl, 5_8)

        ! Loop over shells
        ishell_start = nucleus_index(inucl) + 1
        ishell_end   = nucleus_index(inucl) + nucleus_shell_num(inucl)
        do ishell = ishell_start, ishell_end
           k = ao_index(ishell)
           l = shell_ang_mom(ishell)
           do il = lstart(l), lstart(l+1)-1

              do i = 1, 3
                 do j = 1, 3
                    ao_hessian(k,j,ipoint,i) = (&
                         poly_hessian(j,i,il) * shell_vgl(ishell,1,ipoint) + &
                         poly_vgl(1,il) * shell_hessian(ishell,j,i,ipoint) + &
                         poly_vgl(j+1,il) * shell_vgl(ishell,i+1,ipoint) + &
                         poly_vgl(i+1,il) * shell_vgl(ishell,j+1,ipoint) &
                         ) * ao_factor(k)
                 end do

                 ao_hessian(k,4,ipoint,i) = (&
                      poly_hessian(i,4,il) * shell_vgl(ishell,1,ipoint) + &
                      poly_vgl(1,il) * shell_hessian(ishell,4,i,ipoint) + &
                      poly_vgl(5,il) * shell_vgl(ishell,i+1,ipoint) + &
                      poly_vgl(i+1,il) * shell_vgl(ishell,5,ipoint) + &
                      2.d0 * (&
                      poly_hessian(1,i,il) * shell_vgl(ishell,2,ipoint) + &
                      poly_vgl(2,il) * shell_hessian(ishell,1,i,ipoint) + &
                      poly_hessian(2,i,il) * shell_vgl(ishell,3,ipoint) + &
                      poly_vgl(3,il) * shell_hessian(ishell,2,i,ipoint) + &
                      poly_hessian(3,i,il) * shell_vgl(ishell,4,ipoint) + &
                      poly_vgl(4,il) * shell_hessian(ishell,3,i,ipoint) &
                      )) * ao_factor(k)

              end do

              k = k+1
           end do
        end do
     end do
  end do

  deallocate(poly_vgl, powers, poly_hessian, ao_index)
end function qmckl_compute_ao_hessian_doc
