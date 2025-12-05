! Fortran interface


interface
  integer(c_int32_t) function qmckl_set_ao_basis_type (context, &
       basis_type) bind(C)
    use qmckl_constants
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    character(c_char)   , intent(in)  , value :: basis_type
  end function qmckl_set_ao_basis_type
end interface

interface
  integer(c_int32_t) function qmckl_set_ao_basis_shell_num(context, &
       num) bind(C)
    use qmckl_constants
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(in)  , value :: num
  end function qmckl_set_ao_basis_shell_num
end interface

interface
  integer(c_int32_t) function qmckl_set_ao_basis_prim_num(context, &
       num) bind(C)
    use qmckl_constants
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(in)  , value :: num
  end function qmckl_set_ao_basis_prim_num
end interface

interface
  integer(c_int32_t) function qmckl_set_ao_basis_nucleus_index(context, &
       idx, size_max) bind(C)
    use qmckl_constants
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(in)          :: idx(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_set_ao_basis_nucleus_index
end interface

interface
  integer(c_int32_t) function qmckl_set_ao_basis_nucleus_shell_num(context, &
       shell_num, size_max) bind(C)
    use qmckl_constants
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(in)          :: shell_num(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_set_ao_basis_nucleus_shell_num
end interface

interface
  integer(c_int32_t) function qmckl_set_ao_basis_shell_ang_mom(context, &
       shell_ang_mom, size_max) bind(C)
    use qmckl_constants
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int32_t) , intent(in)          :: shell_ang_mom(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_set_ao_basis_shell_ang_mom
end interface

interface
  integer(c_int32_t) function qmckl_set_ao_basis_shell_prim_num(context, &
       shell_prim_num, size_max) bind(C)
    use qmckl_constants
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(in)          :: shell_prim_num(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_set_ao_basis_shell_prim_num
end interface

interface
  integer(c_int32_t) function qmckl_set_ao_basis_shell_prim_index(context, &
       shell_prim_index, size_max) bind(C)
    use qmckl_constants
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(in)          :: shell_prim_index(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_set_ao_basis_shell_prim_index
end interface

interface
  integer(c_int32_t) function qmckl_set_ao_basis_shell_factor(context, &
       shell_factor, size_max) bind(C)
    use qmckl_constants
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double)  , intent(in)          :: shell_factor(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_set_ao_basis_shell_factor
end interface

interface
  integer(c_int32_t) function qmckl_set_ao_basis_exponent(context, &
       exponent, size_max) bind(C)
    use qmckl_constants
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double)  , intent(in)          :: exponent(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_set_ao_basis_exponent
end interface

interface
  integer(c_int32_t) function qmckl_set_ao_basis_coefficient(context, &
       coefficient, size_max)  bind(C)
    use qmckl_constants
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double)  , intent(in)          :: coefficient(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_set_ao_basis_coefficient
end interface

interface
  integer(c_int32_t) function qmckl_set_ao_basis_prim_factor(context, &
       prim_factor, size_max) bind(C)
    use qmckl_constants
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double)  , intent(in)          :: prim_factor(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_set_ao_basis_prim_factor
end interface

interface
  integer(c_int32_t) function qmckl_set_ao_basis_ao_num(context, &
       num) bind(C)
    use qmckl_constants
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(in)  , value :: num
  end function qmckl_set_ao_basis_ao_num
end interface

interface
  integer(c_int32_t) function qmckl_set_ao_basis_cartesian(context, &
       cartesian) bind(C)
    use qmckl_constants
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    logical (c_bool)    , intent(in)  , value :: cartesian
  end function qmckl_set_ao_basis_cartesian
end interface

interface
  integer(c_int32_t) function qmckl_set_ao_basis_ao_factor(context, &
       ao_factor, size_max) bind(C)
    use qmckl_constants
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double)  , intent(in)          :: ao_factor(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_set_ao_basis_ao_factor
end interface

! Fortran interface


interface
  integer(c_int32_t) function qmckl_get_ao_basis_type (context, &
       basis_type) bind(C)
    use qmckl_constants
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    character(c_char)   , intent(out)         :: basis_type
  end function qmckl_get_ao_basis_type
end interface

interface
  integer(c_int32_t) function qmckl_get_ao_basis_shell_num(context, &
       num) bind(C)
    use qmckl_constants
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(out)         :: num
  end function qmckl_get_ao_basis_shell_num
end interface

interface
  integer(c_int32_t) function qmckl_get_ao_basis_prim_num(context, &
       num) bind(C)
    use qmckl_constants
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(out)         :: num
  end function qmckl_get_ao_basis_prim_num
end interface

interface
  integer(c_int32_t) function qmckl_get_ao_basis_nucleus_shell_num(context, &
       shell_num, size_max) bind(C)
    use qmckl_constants
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(out)         :: shell_num(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_get_ao_basis_nucleus_shell_num
end interface

interface
  integer(c_int32_t) function qmckl_get_ao_basis_nucleus_index(context, &
       idx, size_max) bind(C)
    use qmckl_constants
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(out)         :: idx(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_get_ao_basis_nucleus_index
end interface

interface
  integer(c_int32_t) function qmckl_get_ao_basis_shell_ang_mom(context, &
       shell_ang_mom, size_max) bind(C)
    use qmckl_constants
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int32_t) , intent(out)         :: shell_ang_mom(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_get_ao_basis_shell_ang_mom
end interface

interface
  integer(c_int32_t) function qmckl_get_ao_basis_shell_prim_num(context, &
       shell_prim_num, size_max) bind(C)
    use qmckl_constants
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(out)         :: shell_prim_num(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_get_ao_basis_shell_prim_num
end interface

interface
  integer(c_int32_t) function qmckl_get_ao_basis_shell_prim_index(context, &
       shell_prim_index, size_max) bind(C)
    use qmckl_constants
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(out)         :: shell_prim_index(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_get_ao_basis_shell_prim_index
end interface

interface
  integer(c_int32_t) function qmckl_get_ao_basis_shell_factor(context, &
       shell_factor, size_max) bind(C)
    use qmckl_constants
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double)  , intent(out)         :: shell_factor(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_get_ao_basis_shell_factor
end interface

interface
  integer(c_int32_t) function qmckl_get_ao_basis_exponent(context, &
       exponent, size_max) bind(C)
    use qmckl_constants
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double)  , intent(out)         :: exponent(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_get_ao_basis_exponent
end interface

interface
  integer(c_int32_t) function qmckl_get_ao_basis_coefficient(context, &
       coefficient, size_max) bind(C)
    use qmckl_constants
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double)  , intent(out)         :: coefficient(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_get_ao_basis_coefficient
end interface

interface
  integer(c_int32_t) function qmckl_get_ao_basis_prim_factor(context, &
       prim_factor, size_max) bind(C)
    use qmckl_constants
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double)  , intent(out)         :: prim_factor(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_get_ao_basis_prim_factor
end interface

interface
  integer(c_int32_t) function qmckl_get_ao_basis_ao_num(context, &
       num) bind(C)
    use qmckl_constants
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(out)         :: num
  end function qmckl_get_ao_basis_ao_num
end interface

interface
  integer(c_int32_t) function qmckl_get_ao_basis_cartesian(context, &
       cartesian) bind(C)
    use qmckl_constants
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    logical (c_bool)    , intent(out)         :: cartesian
  end function qmckl_get_ao_basis_cartesian
end interface

interface
  integer(c_int32_t) function qmckl_get_ao_basis_ao_factor(context, &
       ao_factor, size_max) bind(C)
    use qmckl_constants
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double)  , intent(out)         :: ao_factor(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_get_ao_basis_ao_factor
end interface

interface
  integer(c_int32_t) function qmckl_get_ao_basis_primitive_vgl &
       (context, primitive_vgl, size_max) &
       bind(C)
    use qmckl_constants
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    real(c_double),       intent(out)         :: primitive_vgl(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function
end interface

interface
  integer(c_int32_t) function qmckl_get_ao_basis_shell_vgl &
       (context, shell_vgl, size_max) &
       bind(C)
    use qmckl_constants
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    real(c_double),       intent(out)         :: shell_vgl(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function
end interface

interface
   integer(c_int32_t) function qmckl_get_ao_basis_ao_vgl (context, &
        ao_vgl, size_max) bind(C)
     use qmckl_constants
     import
     implicit none
     integer (c_int64_t) , intent(in)  , value :: context
     real(c_double),       intent(out)         :: ao_vgl(*)
     integer (c_int64_t) , intent(in)  , value :: size_max
   end function qmckl_get_ao_basis_ao_vgl
end interface

interface
   integer(c_int32_t) function qmckl_get_ao_basis_ao_vgl_inplace (context, &
        ao_vgl, size_max) bind(C)
     use qmckl_constants
     import
     implicit none
     integer (c_int64_t) , intent(in)  , value :: context
     real(c_double),       intent(out)         :: ao_vgl(*)
     integer (c_int64_t) , intent(in)  , value :: size_max
   end function qmckl_get_ao_basis_ao_vgl_inplace
end interface

interface
   integer(c_int32_t) function qmckl_get_ao_basis_ao_value (context, &
        ao_value, size_max) bind(C)
     use qmckl_constants
     import
     implicit none
     integer (c_int64_t) , intent(in)  , value :: context
     real(c_double)      , intent(out)         :: ao_value(*)
     integer (c_int64_t) , intent(in)  , value :: size_max
   end function qmckl_get_ao_basis_ao_value
end interface

interface
   integer(c_int32_t) function qmckl_get_ao_basis_ao_value_inplace (context, &
        ao_value, size_max) bind(C)
     use qmckl_constants
     import
     implicit none
     integer (c_int64_t) , intent(in)  , value :: context
     real(c_double)      , intent(out)         :: ao_value(*)
     integer (c_int64_t) , intent(in)  , value :: size_max
   end function qmckl_get_ao_basis_ao_value_inplace
end interface

interface
   function qmckl_ao_gaussian_vgl(context, &
        X, R, n, A, VGL, ldv) bind(C) result(info)
     use qmckl_constants
     integer (qmckl_context) , intent(in) , value :: context
     integer (c_int64_t) , intent(in) , value :: ldv
     integer (c_int64_t) , intent(in) , value :: n
     real    (c_double)  , intent(in)         :: X(3), R(3), A(n)
     real    (c_double)  , intent(out)        :: VGL(ldv,5)
     integer(qmckl_exit_code) :: info
   end function qmckl_ao_gaussian_vgl
end interface

interface
  integer(c_int32_t) function qmckl_get_ao_basis_shell_hessian &
       (context, shell_hessian, size_max) &
       bind(C)
    use qmckl_constants
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    real(c_double),       intent(out)         :: shell_hessian(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function
end interface



! #+CALL: generate_f_interface(table=qmckl_ao_power_args,rettyp=get_value("CRetType"),fname="qmckl_ao_power")

! #+RESULTS:

interface
  integer(qmckl_exit_code) function qmckl_ao_power &
      (context, n, X, LMAX, P, ldp) &
      bind(C)
    use :: qmckl_constants
    import
    implicit none

    integer (qmckl_context), intent(in)  , value :: context
    integer (c_int64_t) , intent(in)  , value :: n
    integer (c_int64_t) , intent(in)  , value :: ldp
    real    (c_double ) , intent(in)          :: X(n)
    integer (c_int32_t) , intent(in)          :: LMAX(n)
    real    (c_double ) , intent(out)         :: P(ldp,n)

  end function qmckl_ao_power
end interface



! #+CALL: generate_f_interface(table=qmckl_ao_polynomial_vgl_args,rettyp=get_value("FRetType"),fname="qmckl_ao_polynomial_vgl_doc" )

! #+RESULTS:

interface
  integer(qmckl_exit_code) function qmckl_ao_polynomial_vgl_doc &
      (context, X, R, lmax, n, L, ldl, VGL, ldv) &
      bind(C)
    use qmckl_constants
    import
    implicit none

    integer (qmckl_context), intent(in)  , value :: context
    real    (c_double ) , intent(in)          :: X(3)
    real    (c_double ) , intent(in)          :: R(3)
    integer (c_int32_t) , intent(in)  , value :: lmax
    integer (c_int64_t) , intent(inout)        :: n
    integer (c_int64_t) , intent(in)  , value :: ldl
    integer (c_int64_t) , intent(in)  , value :: ldv
    integer (c_int32_t) , intent(out)         :: L(ldl,n)
    real    (c_double ) , intent(out)         :: VGL(ldv,n)

  end function qmckl_ao_polynomial_vgl_doc
end interface



! #+CALL: generate_f_interface(table=qmckl_ao_polynomial_vgl_args,rettyp=get_value("FRetType"),fname="qmckl_ao_polynomial_vgl" )

! #+RESULTS:

interface
  integer(qmckl_exit_code) function qmckl_ao_polynomial_vgl &
      (context, X, R, lmax, n, L, ldl, VGL, ldv) &
      bind(C)
    use qmckl_constants
    import
    implicit none

    integer (qmckl_context), intent(in)  , value :: context
    real    (c_double ) , intent(in)          :: X(3)
    real    (c_double ) , intent(in)          :: R(3)
    integer (c_int32_t) , intent(in)  , value :: lmax
    integer (c_int64_t) , intent(inout)       :: n
    integer (c_int64_t) , intent(in)  , value :: ldl
    integer (c_int64_t) , intent(in)  , value :: ldv
    integer (c_int32_t) , intent(out)         :: L(ldl,n)
    real    (c_double ) , intent(out)         :: VGL(ldv,n)

  end function qmckl_ao_polynomial_vgl
end interface

interface
  integer(qmckl_exit_code) function qmckl_ao_polynomial_transp_vgl &
      (context, X, R, lmax, n, L, ldl, VGL, ldv) &
      bind(C)
    use qmckl_constants
    import
    implicit none

    integer (qmckl_context), intent(in)  , value :: context
    real    (c_double ) , intent(in)          :: X(3)
    real    (c_double ) , intent(in)          :: R(3)
    integer (c_int32_t) , intent(in)  , value :: lmax
    integer (c_int64_t) , intent(inout)        :: n
    integer (c_int64_t) , intent(in)  , value :: ldl
    integer (c_int64_t) , intent(in)  , value :: ldv
    integer (c_int32_t) , intent(out)         :: L(ldl,n)
    real    (c_double ) , intent(out)         :: VGL(ldv,5)

  end function qmckl_ao_polynomial_transp_vgl
end interface

interface
  integer(qmckl_exit_code) function qmckl_ao_polynomial_hessian &
      (context, X, R, lmax, n, hessian) &
      bind(C)
    use qmckl_constants
    import
    implicit none

    integer (qmckl_context), intent(in)  , value :: context
    real    (c_double ) , intent(in)          :: X(3)
    real    (c_double ) , intent(in)          :: R(3)
    integer (c_int32_t) , intent(in)  , value :: lmax
    integer (c_int64_t) , intent(inout)       :: n
    real    (c_double ) , intent(out)         :: hessian(3,4,(lmax+1)*(lmax+2)*(lmax+3)/6)

  end function qmckl_ao_polynomial_hessian
end interface

interface
  integer(qmckl_exit_code) function qmckl_get_ao_basis_ao_hessian(context, ao_hessian, size_max) &
    bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double ) , intent(out)         :: ao_hessian(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function
end interface
