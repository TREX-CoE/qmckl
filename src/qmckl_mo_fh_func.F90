interface
  integer(qmckl_exit_code) function qmckl_set_mo_basis_coefficient (context, &
       coefficient, size_max) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real(c_double)      , intent(in)          :: coefficient(*)
    integer (c_int64_t) , intent(in), value   :: size_max
  end function qmckl_set_mo_basis_coefficient
end interface

! Fortran interfaces


interface
  integer(qmckl_exit_code) function qmckl_get_mo_basis_mo_num (context, &
       mo_num) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(out)         :: mo_num
  end function qmckl_get_mo_basis_mo_num
end interface

interface
  integer(qmckl_exit_code) function qmckl_get_mo_basis_coefficient(context, &
       coefficient, size_max) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in), value   :: context
    real(c_double)      , intent(out)         :: coefficient(*)
    integer (c_int64_t) , intent(in), value   :: size_max
  end function qmckl_get_mo_basis_coefficient
end interface

interface
  integer(qmckl_exit_code) function qmckl_set_mo_basis_r_cusp(context, &
       r_cusp, size_max) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in), value   :: context
    real(c_double)      , intent(in)          :: r_cusp(*)
    integer (c_int64_t) , intent(in), value   :: size_max
  end function qmckl_set_mo_basis_r_cusp
end interface

! Fortran interface


interface
  integer(qmckl_exit_code) function qmckl_mo_basis_select_mo (context, &
       keep, size_max) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in), value :: context
    integer (c_int32_t) , intent(in)        :: keep(*)
    integer (c_int64_t) , intent(in), value :: size_max
  end function qmckl_mo_basis_select_mo
end interface

interface
   integer(qmckl_exit_code) function qmckl_get_mo_basis_mo_value (context, &
        mo_value, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none

     integer (c_int64_t) , intent(in)  , value :: context
     real(c_double),       intent(out)         :: mo_value(*)
     integer (c_int64_t) , intent(in)  , value :: size_max
   end function qmckl_get_mo_basis_mo_value
end interface

interface
   integer(qmckl_exit_code) function qmckl_get_mo_basis_mo_value_inplace (context, &
        mo_value, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (c_int64_t) , intent(in)  , value :: context
     real(c_double),       intent(out)         :: mo_value(*)
     integer (c_int64_t) , intent(in)  , value :: size_max
   end function qmckl_get_mo_basis_mo_value_inplace
end interface

interface
   integer(qmckl_exit_code) function qmckl_get_mo_basis_mo_vgl (context, &
        mo_vgl, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none

     integer (c_int64_t) , intent(in)  , value :: context
     real(c_double),       intent(out)         :: mo_vgl(*)
     integer (c_int64_t) , intent(in)  , value :: size_max
   end function qmckl_get_mo_basis_mo_vgl
end interface

interface
   integer(qmckl_exit_code) function qmckl_get_mo_basis_mo_vgl_inplace (context, &
        mo_vgl, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (c_int64_t) , intent(in)  , value :: context
     real(c_double),       intent(out)         :: mo_vgl(*)
     integer (c_int64_t) , intent(in)  , value :: size_max
   end function qmckl_get_mo_basis_mo_vgl_inplace
end interface

! Fortran interface


interface
  integer(qmckl_exit_code) function qmckl_mo_basis_rescale (context, &
       scaling_factor) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in), value :: context
    real    (c_double)  , intent(in), value :: scaling_factor
  end function qmckl_mo_basis_rescale
end interface
