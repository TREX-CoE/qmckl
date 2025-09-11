integer, parameter :: QMCKL_DEFAULT_PRECISION        = 53
integer, parameter :: QMCKL_DEFAULT_RANGE            = 11

interface
   integer (qmckl_exit_code) function qmckl_set_numprec_precision(context, precision) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer (qmckl_context), intent(in), value :: context
     integer (c_int32_t), intent(in), value :: precision
   end function qmckl_set_numprec_precision
end interface

interface
   integer (qmckl_exit_code) function qmckl_get_numprec_precision(context) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer (qmckl_context), intent(in), value :: context
   end function qmckl_get_numprec_precision
end interface

interface
   integer (qmckl_exit_code) function qmckl_set_numprec_range(context, range) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer (qmckl_context), intent(in), value :: context
     integer (c_int32_t), intent(in), value :: range
   end function qmckl_set_numprec_range
end interface

interface
   integer (qmckl_exit_code) function qmckl_get_numprec_range(context) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer (qmckl_context), intent(in), value :: context
   end function qmckl_get_numprec_range
end interface

interface
   real (c_double) function qmckl_get_numprec_epsilon(context) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer (qmckl_context), intent(in), value :: context
   end function qmckl_get_numprec_epsilon
end interface

interface
   integer (c_int) function qmckl_test_precision_32(a,b) bind(C)
     use, intrinsic :: iso_c_binding
     import
     real (c_float), intent(in), value :: a, b
   end function qmckl_test_precision_32
end interface

interface
   integer (c_int) function qmckl_test_precision_64(a,b) bind(C)
     use, intrinsic :: iso_c_binding
     import
     real (c_double), intent(in), value :: a, b
   end function qmckl_test_precision_64
end interface
