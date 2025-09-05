interface
  integer(qmckl_exit_code) function qmckl_set_single_point(context, &
       transp, num, coord, size_max) bind(C, name="qmckl_set_single_point_f")
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    character(c_char)   , intent(in)  , value :: transp
    integer (c_int64_t) , intent(in)  , value :: num
    real    (c_double ) , intent(in)          :: coord(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function
end interface

interface
  integer(qmckl_exit_code) function qmckl_get_single_electron_ee_distance(context, distance, size_max) &
    bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double ) , intent(out)         :: distance(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function
end interface

interface
  integer(qmckl_exit_code) function qmckl_get_single_electron_en_distance(context, distance, size_max) &
    bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double ) , intent(out)         :: distance(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function
end interface

interface
   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_single_een (context, &
        delta_een, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
     integer(c_int64_t), intent(in), value       :: size_max
     real(c_double),   intent(out)               :: delta_een(size_max)
   end function
end interface

interface
   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_single_een_gl (context, &
        delta_een_gl, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
     integer(c_int64_t), intent(in), value       :: size_max
     real(c_double),   intent(out)               :: delta_een_gl(size_max)
   end function
end interface

interface
   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_single_een_g (context, &
        delta_een_g, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
     integer(c_int64_t), intent(in), value       :: size_max
     real(c_double),   intent(out)               :: delta_een_g(size_max)
   end function
end interface

interface
   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_single_ee (context, &
        delta_ee, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
     integer(c_int64_t), intent(in), value       :: size_max
     real(c_double),   intent(out)               :: delta_ee(size_max)
   end function qmckl_get_jastrow_champ_single_ee
end interface

interface
   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_single_ee_gl (context, &
        delta_ee_gl, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
     integer(c_int64_t), intent(in), value       :: size_max
     real(c_double),   intent(out)               :: delta_ee_gl(size_max)
   end function qmckl_get_jastrow_champ_single_ee_gl
end interface

interface
   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_single_en (context, &
        delta_en, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
     integer(c_int64_t), intent(in), value       :: size_max
     real(c_double),   intent(out)               :: delta_en(size_max)
   end function qmckl_get_jastrow_champ_single_en
end interface

interface
   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_single_en_gl (context, &
        delta_en_gl, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
     integer(c_int64_t), intent(in), value       :: size_max
     real(c_double),   intent(out)               :: delta_en_gl(size_max)
   end function qmckl_get_jastrow_champ_single_en_gl
end interface

interface
   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_single_accept (context) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
   end function qmckl_get_jastrow_champ_single_accept
end interface
