interface
   integer(qmckl_exit_code) function qmckl_get_forces_jastrow_en (context, &
        forces_jastrow_en, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
     integer(c_int64_t), intent(in), value       :: size_max
     real(c_double),   intent(out)               :: forces_jastrow_en(size_max)
   end function qmckl_get_forces_jastrow_en
end interface

interface
   integer(qmckl_exit_code) function qmckl_get_forces_jastrow_en_g (context, &
        forces_jastrow_en_g, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
     integer(c_int64_t), intent(in), value       :: size_max
     real(c_double),   intent(out)               :: forces_jastrow_en_g(size_max)
   end function qmckl_get_forces_jastrow_en_g
end interface

interface
   integer(qmckl_exit_code) function qmckl_get_forces_jastrow_en_l (context, &
        forces_jastrow_en_l, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
     integer(c_int64_t), intent(in), value       :: size_max
     real(c_double),   intent(out)               :: forces_jastrow_en_l(size_max)
   end function qmckl_get_forces_jastrow_en_l
end interface

interface
   integer(qmckl_exit_code) function qmckl_get_forces_jastrow_single_en (context, &
        forces_jastrow_single_en, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none

     integer (c_int64_t) , intent(in)  , value :: context
     integer (c_int64_t) , intent(in)  , value :: size_max
     real(c_double),       intent(out)         :: forces_jastrow_single_en(size_max)
   end function qmckl_get_forces_jastrow_single_en
end interface

interface
   integer(qmckl_exit_code) function qmckl_get_forces_jastrow_een (context, &
        forces_jastrow_een, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
     integer(c_int64_t), intent(in), value       :: size_max
     real(c_double),   intent(out)               :: forces_jastrow_een(size_max)
   end function qmckl_get_forces_jastrow_een
end interface

interface
   integer(qmckl_exit_code) function qmckl_get_forces_jastrow_een_g (context, &
        forces_jastrow_een_g, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
     integer(c_int64_t), intent(in), value       :: size_max
     real(c_double),   intent(out)               :: forces_jastrow_een_g(size_max)
   end function qmckl_get_forces_jastrow_een_g
end interface

interface
   integer(qmckl_exit_code) function qmckl_get_forces_jastrow_een_l (context, &
        forces_jastrow_een_l, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
     integer(c_int64_t), intent(in), value       :: size_max
     real(c_double),   intent(out)               :: forces_jastrow_een_l(size_max)
   end function qmckl_get_forces_jastrow_een_l
end interface

interface
   integer(qmckl_exit_code) function qmckl_get_forces_jastrow_single_een (context, &
        forces_jastrow_single_een, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none

     integer (c_int64_t) , intent(in)  , value :: context
     integer (c_int64_t) , intent(in)  , value :: size_max
     real(c_double),       intent(out)         :: forces_jastrow_single_een(size_max)
   end function qmckl_get_forces_jastrow_single_een
end interface

interface
   integer(qmckl_exit_code) function qmckl_get_forces_mo_value (context, &
        forces_mo_value, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none

     integer (c_int64_t) , intent(in)  , value :: context
     real(c_double),       intent(out)         :: forces_mo_value(*)
     integer (c_int64_t) , intent(in)  , value :: size_max
   end function qmckl_get_forces_mo_value
end interface

interface
   integer(qmckl_exit_code) function qmckl_get_forces_mo_value_inplace (context, &
        forces_mo_value, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (c_int64_t) , intent(in)  , value :: context
     real(c_double),       intent(out)         :: forces_mo_value(*)
     integer (c_int64_t) , intent(in)  , value :: size_max
   end function qmckl_get_forces_mo_value_inplace
end interface

interface
   integer(qmckl_exit_code) function qmckl_get_forces_mo_g (context, &
        forces_mo_g, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none

     integer (c_int64_t) , intent(in)  , value :: context
     real(c_double),       intent(out)         :: forces_mo_g(*)
     integer (c_int64_t) , intent(in)  , value :: size_max
   end function qmckl_get_forces_mo_g
end interface

interface
   integer(qmckl_exit_code) function qmckl_get_forces_mo_g_inplace (context, &
        forces_mo_g, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none

     integer (c_int64_t) , intent(in)  , value :: context
     real(c_double),       intent(out)         :: forces_mo_g(*)
     integer (c_int64_t) , intent(in)  , value :: size_max
   end function qmckl_get_forces_mo_g_inplace
end interface

interface
   integer(qmckl_exit_code) function qmckl_get_forces_mo_l (context, &
        forces_mo_l, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none

     integer (c_int64_t) , intent(in)  , value :: context
     real(c_double),       intent(out)         :: forces_mo_l(*)
     integer (c_int64_t) , intent(in)  , value :: size_max
   end function qmckl_get_forces_mo_l
end interface
