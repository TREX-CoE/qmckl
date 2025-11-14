! Fortran interface


interface
   integer(qmckl_exit_code) function qmckl_set_jastrow_champ_rescale_factor_ee (context, &
        kappa_ee) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in)  , value :: context
     real(c_double),   intent(in), value  :: kappa_ee
   end function qmckl_set_jastrow_champ_rescale_factor_ee

   integer(qmckl_exit_code) function qmckl_set_jastrow_champ_rescale_factor_en (context, &
        kappa_en, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in)  , value :: context
     integer(c_int64_t), intent(in), value  :: size_max
     real(c_double),   intent(in) :: kappa_en(size_max)
   end function qmckl_set_jastrow_champ_rescale_factor_en

   integer(qmckl_exit_code) function qmckl_set_jastrow_champ_aord_num (context, &
        aord_num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in)  , value :: context
     integer(c_int64_t), intent(in), value  :: aord_num
   end function qmckl_set_jastrow_champ_aord_num

   integer(qmckl_exit_code) function qmckl_set_jastrow_champ_bord_num (context, &
        bord_num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in)  , value :: context
     integer(c_int64_t), intent(in), value  :: bord_num
   end function qmckl_set_jastrow_champ_bord_num

   integer(qmckl_exit_code) function qmckl_set_jastrow_champ_cord_num (context, &
        cord_num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in)  , value :: context
     integer(c_int64_t), intent(in), value  :: cord_num
   end function qmckl_set_jastrow_champ_cord_num

   integer(qmckl_exit_code) function qmckl_set_jastrow_champ_type_nucl_num (context, &
        type_nucl_num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in)  , value :: context
     integer(c_int64_t), intent(in), value  :: type_nucl_num
   end function qmckl_set_jastrow_champ_type_nucl_num

   integer(qmckl_exit_code) function qmckl_set_jastrow_champ_type_nucl_vector (context, &
        type_nucl_vector, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in)  , value :: context
     integer(c_int64_t), intent(in), value  :: size_max
     integer(c_int64_t), intent(in) :: type_nucl_vector(size_max)
   end function qmckl_set_jastrow_champ_type_nucl_vector

   integer(qmckl_exit_code) function qmckl_set_jastrow_champ_a_vector(context, &
        a_vector, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in)  , value :: context
     integer(c_int64_t), intent(in), value  :: size_max
     real(c_double),   intent(in) :: a_vector(size_max)
   end function qmckl_set_jastrow_champ_a_vector

   integer(qmckl_exit_code) function qmckl_set_jastrow_champ_b_vector(context, &
        b_vector, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in)  , value :: context
     integer(c_int64_t), intent(in), value  :: size_max
     real(c_double),   intent(in) :: b_vector(size_max)
   end function qmckl_set_jastrow_champ_b_vector

   integer(qmckl_exit_code) function qmckl_set_jastrow_champ_c_vector(context, &
        c_vector, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in)  , value :: context
     integer(c_int64_t), intent(in), value  :: size_max
     real(c_double),   intent(in) :: c_vector(size_max)
   end function qmckl_set_jastrow_champ_c_vector

   integer(qmckl_exit_code) function qmckl_set_jastrow_champ_spin_independent(context, &
        spin_independent) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer(qmckl_context) , intent(in)  , value :: context
     integer(c_int32_t),    intent(in), value  :: spin_independent
   end function qmckl_set_jastrow_champ_spin_independent

end interface

! Fortran interface


interface
   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_rescale_factor_ee (context, &
        kappa_ee) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in)  , value :: context
     real(c_double),   intent(out) :: kappa_ee
   end function qmckl_get_jastrow_champ_rescale_factor_ee

   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_rescale_factor_en (context, &
        kappa_en, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
     integer(c_int64_t), intent(in), value       :: size_max
     real(c_double),   intent(out)               :: kappa_en(size_max)
   end function qmckl_get_jastrow_champ_rescale_factor_en

   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_aord_num (context, &
        aord_num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
     integer(c_int64_t), intent(out)             :: aord_num
   end function qmckl_get_jastrow_champ_aord_num

   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_bord_num (context, &
        bord_num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
     integer(c_int64_t), intent(out)             :: bord_num
   end function qmckl_get_jastrow_champ_bord_num

   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_cord_num (context, &
        cord_num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
     integer(c_int64_t), intent(out)             :: cord_num
   end function qmckl_get_jastrow_champ_cord_num

   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_type_nucl_num (context, &
        type_nucl_num) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
     integer(c_int64_t), intent(out)              :: type_nucl_num
   end function qmckl_get_jastrow_champ_type_nucl_num

   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_type_nucl_vector (context, &
        type_nucl_vector, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context), intent(in), value :: context
     integer(c_int64_t), intent(in), value      :: size_max
     integer(c_int64_t), intent(out)            :: type_nucl_vector(size_max)
   end function qmckl_get_jastrow_champ_type_nucl_vector

   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_a_vector(context, &
        a_vector, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
     integer(c_int64_t), intent(in), value       :: size_max
     real(c_double),   intent(out)               :: a_vector(size_max)
   end function qmckl_get_jastrow_champ_a_vector

   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_b_vector(context, &
        b_vector, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in)  , value :: context
     integer(c_int64_t), intent(in), value         :: size_max
     real(c_double),   intent(out)                 :: b_vector(size_max)
   end function qmckl_get_jastrow_champ_b_vector

   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_c_vector(context, &
        c_vector, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in)  , value :: context
     integer(c_int64_t), intent(in), value         :: size_max
     real(c_double),   intent(out)                 :: c_vector(size_max)
   end function qmckl_get_jastrow_champ_c_vector

   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_spin_independent(context, &
        spin_independent) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer(qmckl_context) , intent(in)  , value :: context
     integer(c_int32_t),   intent(out)                :: spin_independent
   end function qmckl_get_jastrow_champ_spin_independent

end interface

! Fortran interface


interface
   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_asymp_jasb(context, &
        asymp_jasb, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
     integer(c_int64_t), intent(in), value       :: size_max
     real(c_double),   intent(out)               :: asymp_jasb(size_max)
   end function qmckl_get_jastrow_champ_asymp_jasb
end interface

! Fortran interface


interface
   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_factor_ee (context, &
        factor_ee, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
     integer(c_int64_t), intent(in), value       :: size_max
     real(c_double),   intent(out)               :: factor_ee(size_max)
   end function qmckl_get_jastrow_champ_factor_ee
end interface

! Fortran interface


interface
   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_factor_ee_gl (context, &
        factor_ee_gl, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
     integer(c_int64_t), intent(in), value       :: size_max
     real(c_double),   intent(out)               :: factor_ee_gl(size_max)
   end function qmckl_get_jastrow_champ_factor_ee_gl
end interface

! Fortran interface


interface
   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_asymp_jasa(context, &
        asymp_jasa, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
     integer(c_int64_t), intent(in), value       :: size_max
     real(c_double),   intent(out)               :: asymp_jasa(size_max)
   end function qmckl_get_jastrow_champ_asymp_jasa
end interface

! Fortran interface


interface
   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_factor_en (context, &
        factor_en, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
     integer(c_int64_t), intent(in), value       :: size_max
     real(c_double),   intent(out)               :: factor_en(size_max)
   end function qmckl_get_jastrow_champ_factor_en
end interface

! Fortran interface


interface
   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_factor_en_gl (context, &
        factor_en_gl, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
     integer(c_int64_t), intent(in), value       :: size_max
     real(c_double),   intent(out)               :: factor_en_gl(size_max)
   end function qmckl_get_jastrow_champ_factor_en_gl
end interface

! Fortran interface


interface
   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_factor_een (context, &
        factor_een, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
     integer(c_int64_t), intent(in), value       :: size_max
     real(c_double),   intent(out)               :: factor_een(size_max)
   end function qmckl_get_jastrow_champ_factor_een
end interface

! Fortran interface


interface
   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_factor_een_gl (context, &
        factor_een_gl, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
     integer(c_int64_t), intent(in), value       :: size_max
     real(c_double),   intent(out)               :: factor_een_gl(size_max)
   end function qmckl_get_jastrow_champ_factor_een_gl
end interface

interface
   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_factor_een_grad (context, &
        factor_een_grad, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
     integer(c_int64_t), intent(in), value       :: size_max
     real(c_double),   intent(out)               :: factor_een_grad(size_max)
   end function qmckl_get_jastrow_champ_factor_een_grad
end interface

! Fortran interface


interface
   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_value (context, &
        value, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
     integer(c_int64_t), intent(in), value       :: size_max
     real(c_double),   intent(out)               :: value(size_max)
   end function qmckl_get_jastrow_champ_value
end interface

! Fortran interface


interface
   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_gl (context, &
        gl, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
     integer(c_int64_t), intent(in), value       :: size_max
     real(c_double),   intent(out)               :: gl(size_max)
   end function qmckl_get_jastrow_champ_gl
end interface
