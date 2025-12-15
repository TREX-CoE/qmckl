!
!    ------------------------------------------
!     QMCkl - Quantum Monte Carlo kernel library
!     ------------------------------------------
!    
!     Documentation : https://trex-coe.github.io/qmckl
!     Issues        : https://github.com/trex-coe/qmckl/issues
!    
!     BSD 3-Clause License
!     
!     Copyright (c) 2020, TREX Center of Excellence
!     All rights reserved.
!     
!     Redistribution and use in source and binary forms, with or without
!     modification, are permitted provided that the following conditions are met:
!     
!     1. Redistributions of source code must retain the above copyright notice, this
!        list of conditions and the following disclaimer.
!     
!     2. Redistributions in binary form must reproduce the above copyright notice,
!        this list of conditions and the following disclaimer in the documentation
!        and/or other materials provided with the distribution.
!     
!     3. Neither the name of the copyright holder nor the names of its
!        contributors may be used to endorse or promote products derived from
!        this software without specific prior written permission.
!     
!     THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!     AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!     IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!     DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
!     FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
!     DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
!     SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!     CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
!     OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!     OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!     
!     
!    
!    
!
module qmckl_constants
  use, intrinsic :: iso_c_binding
integer  , parameter :: qmckl_context = c_int64_t
integer*8, parameter :: QMCKL_NULL_CONTEXT = 0
integer  , parameter :: qmckl_exit_code = c_int32_t

integer(qmckl_exit_code), parameter :: QMCKL_SUCCESS                  = 0
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_1            = 1
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_2            = 2
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_3            = 3
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_4            = 4
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_5            = 5
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_6            = 6
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_7            = 7
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_8            = 8
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_9            = 9
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_10           = 10
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_11           = 11
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_12           = 12
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_13           = 13
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_14           = 14
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_15           = 15
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_16           = 16
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_17           = 17
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_18           = 18
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_19           = 19
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_ARG_20           = 20
integer(qmckl_exit_code), parameter :: QMCKL_FAILURE                  = 101
integer(qmckl_exit_code), parameter :: QMCKL_ERRNO                    = 102
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_CONTEXT          = 103
integer(qmckl_exit_code), parameter :: QMCKL_ALLOCATION_FAILED        = 104
integer(qmckl_exit_code), parameter :: QMCKL_DEALLOCATION_FAILED      = 105
integer(qmckl_exit_code), parameter :: QMCKL_NOT_PROVIDED             = 106
integer(qmckl_exit_code), parameter :: QMCKL_OUT_OF_BOUNDS            = 107
integer(qmckl_exit_code), parameter :: QMCKL_ALREADY_SET              = 108
integer(qmckl_exit_code), parameter :: QMCKL_INVALID_EXIT_CODE        = 109
end module qmckl_constants

module qmckl
  use qmckl_constants
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
  integer(c_int32_t) function qmckl_set_ao_basis_r_power(context, &
       r_power, size_max) bind(C)
    use qmckl_constants
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int32_t) , intent(in)          :: r_power(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_set_ao_basis_r_power
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
  integer(c_int32_t) function qmckl_get_ao_basis_r_power(context, &
       r_power, size_max) bind(C)
    use qmckl_constants
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int32_t) , intent(out)         :: r_power(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_get_ao_basis_r_power
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
   integer(c_int32_t) function qmckl_ao_slater_vgl(context, &
        X, R, num_slater, N, A, VGL, ldv) bind(C)
     use qmckl_constants
     import
     implicit none
     integer (c_int64_t) , intent(in)  , value :: context
     real    (c_double)  , intent(in)          :: X(3), R(3)
     integer (c_int64_t) , intent(in)  , value :: num_slater
     integer (c_int64_t) , intent(in)  , value :: ldv
     integer (c_int64_t) , intent(in)          :: N(num_slater)
     real    (c_double)  , intent(in)          :: A(num_slater)
     real    (c_double)  , intent(out)         :: VGL(ldv,5)
   end function qmckl_ao_slater_vgl
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
! C interface                                                    :noexport:

!     #+CALL: generate_f_interface(table=qmckl_dgemm_args,rettyp="qmckl_exit_code",fname="qmckl_dgemm")

!     #+RESULTS:

interface
  integer(qmckl_exit_code) function qmckl_dgemm &
      (context, TransA, TransB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (qmckl_context), intent(in)  , value :: context
    character(c_char  ) , intent(in)  , value :: TransA
    character(c_char  ) , intent(in)  , value :: TransB
    integer (c_int64_t) , intent(in)  , value :: m
    integer (c_int64_t) , intent(in)  , value :: n
    integer (c_int64_t) , intent(in)  , value :: k
    real    (c_double ) , intent(in)  , value :: alpha
    integer (c_int64_t) , intent(in)  , value :: lda
    integer (c_int64_t) , intent(in)  , value :: ldb
    real    (c_double ) , intent(in)  , value :: beta
    integer (c_int64_t) , intent(in)  , value :: ldc
    real    (c_double ) , intent(in)          :: A(lda,*)
    real    (c_double ) , intent(in)          :: B(ldb,*)
    real    (c_double ) , intent(out)         :: C(ldc,*)

  end function qmckl_dgemm
end interface

! C interface                                                    :noexport:

!     #+CALL: generate_f_interface(table=qmckl_dgemm_safe_args,rettyp="qmckl_exit_code",fname="qmckl_dgemm_safe")

!     #+RESULTS:

interface
  integer(qmckl_exit_code) function qmckl_dgemm_safe &
      (context, TransA, TransB, m, n, k, alpha, A, size_max_A, lda, B, size_max_B, ldb, beta, C, size_max_C, ldc) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (qmckl_context), intent(in)  , value :: context
    character(c_char  ) , intent(in)  , value :: TransA
    character(c_char  ) , intent(in)  , value :: TransB
    integer (c_int64_t) , intent(in)  , value :: m
    integer (c_int64_t) , intent(in)  , value :: n
    integer (c_int64_t) , intent(in)  , value :: k
    real    (c_double ) , intent(in)  , value :: alpha
    integer (c_int64_t) , intent(in)  , value :: size_max_A
    integer (c_int64_t) , intent(in)  , value :: lda
    integer (c_int64_t) , intent(in)  , value :: size_max_B
    integer (c_int64_t) , intent(in)  , value :: ldb
    real    (c_double ) , intent(in)  , value :: beta
    integer (c_int64_t) , intent(in)  , value :: size_max_C
    integer (c_int64_t) , intent(in)  , value :: ldc
    real    (c_double ) , intent(in)          :: A(lda,*)
    real    (c_double ) , intent(in)          :: B(ldb,*)
    real    (c_double ) , intent(out)         :: C(ldc,*)

  end function qmckl_dgemm_safe
end interface



! #+CALL: generate_f_interface(table=qmckl_adjugate_args,rettyp="qmckl_exit_code",fname="qmckl_adjugate")

! #+RESULTS:

interface
  integer(qmckl_exit_code) function qmckl_adjugate &
      (context, n, A, lda, B, ldb, det_l) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (qmckl_context), intent(in)  , value :: context
    integer (c_int64_t) , intent(in)  , value :: n
    integer (c_int64_t) , intent(in)  , value :: lda
    integer (c_int64_t) , intent(in)  , value :: ldb
    real    (c_double ) , intent(inout)        :: det_l
    real    (c_double ) , intent(in)          :: A(lda,*)
    real    (c_double ) , intent(out)         :: B(ldb,*)

  end function qmckl_adjugate
end interface

! C interface

!     #+CALL: generate_f_interface(table=qmckl_adjugate_safe_args,rettyp="qmckl_exit_code",fname="qmckl_adjugate_safe")

!     #+RESULTS:

interface
  integer(qmckl_exit_code) function qmckl_adjugate_safe &
      (context, n, A, size_max_A, lda, B, size_max_B, ldb, det_l) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (qmckl_context), intent(in)  , value :: context
    integer (c_int64_t) , intent(in)  , value :: n
    integer (c_int64_t) , intent(in)  , value :: size_max_A
    integer (c_int64_t) , intent(in)  , value :: lda
    integer (c_int64_t) , intent(in)  , value :: size_max_B
    integer (c_int64_t) , intent(in)  , value :: ldb
    real    (c_double ) , intent(inout)        :: det_l
    real    (c_double ) , intent(in)          :: A(lda,*)
    real    (c_double ) , intent(out)         :: B(ldb,*)

  end function qmckl_adjugate_safe
end interface

! C interface                                                    :noexport:

!     #+CALL: generate_f_interface(table=qmckl_dgemv_args,rettyp="qmckl_exit_code",fname="qmckl_dgemv")

!     #+RESULTS:

interface
  integer(qmckl_exit_code) function qmckl_dgemv &
      (context, Trans, m, n, alpha, A, lda, X, incx, beta, Y, incy) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (qmckl_context), intent(in)  , value :: context
    character(c_char  ) , intent(in)  , value :: Trans
    integer (c_int64_t) , intent(in)  , value :: m
    integer (c_int64_t) , intent(in)  , value :: n
    real    (c_double ) , intent(in)  , value :: alpha
    integer (c_int64_t) , intent(in)  , value :: lda
    integer (c_int64_t) , intent(in)  , value :: incx
    real    (c_double ) , intent(in)  , value :: beta
    integer (c_int64_t) , intent(in)  , value :: incy
    real    (c_double ) , intent(in)          :: A(lda,*)
    real    (c_double ) , intent(in)          :: X(*)
    real    (c_double ) , intent(out)         :: Y(*)

  end function qmckl_dgemv
end interface

! C interface                                                    :noexport:

!     #+CALL: generate_f_interface(table=qmckl_dger_args,rettyp="qmckl_exit_code",fname="qmckl_dger")

!     #+RESULTS:

interface
  integer(qmckl_exit_code) function qmckl_dger &
      (context, m, n, alpha, X, incx, Y, incy, A, lda) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (qmckl_context), intent(in)  , value :: context
    integer (c_int64_t) , intent(in)  , value :: m
    integer (c_int64_t) , intent(in)  , value :: n
    real    (c_double ) , intent(in)  , value :: alpha
    integer (c_int64_t) , intent(in)  , value :: lda
    integer (c_int64_t) , intent(in)  , value :: incx
    integer (c_int64_t) , intent(in)  , value :: incy
    real    (c_double ) , intent(in)          :: X(*)
    real    (c_double ) , intent(in)          :: Y(*)
    real    (c_double ) , intent(inout)       :: A(lda,*)

  end function qmckl_dger
end interface
interface
  integer (qmckl_context) function qmckl_context_touch(context) bind(C)
    use, intrinsic :: iso_c_binding
    import
    integer (qmckl_context), intent(in), value :: context
  end function qmckl_context_touch
end interface

interface
   integer (qmckl_context) function qmckl_context_create() bind(C)
     use, intrinsic :: iso_c_binding
     import
   end function qmckl_context_create
end interface

interface
   integer (qmckl_context) function qmckl_context_copy(context) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer (qmckl_context), intent(in), value :: context
   end function qmckl_context_copy
end interface

interface
   integer (qmckl_exit_code) function qmckl_context_destroy(context) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer (qmckl_context), intent(in), value :: context
   end function qmckl_context_destroy
end interface
! Performance

!     This function is more efficient when ~A~ and ~B~ are
!     transposed.

!    #+CALL: generate_f_interface(table=qmckl_distance_sq_args,fname=get_value("Name"))

!    #+RESULTS:

interface
  integer(qmckl_exit_code) function qmckl_distance_sq &
      (context, transa, transb, m, n, A, lda, B, ldb, C, ldc) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (qmckl_context), intent(in)  , value :: context
    character(c_char  ) , intent(in)  , value :: transa
    character(c_char  ) , intent(in)  , value :: transb
    integer (c_int64_t) , intent(in)  , value :: m
    integer (c_int64_t) , intent(in)  , value :: n
    integer (c_int64_t) , intent(in)  , value :: lda
    integer (c_int64_t) , intent(in)  , value :: ldb
    integer (c_int64_t) , intent(in)  , value :: ldc
    real    (c_double ) , intent(in)          :: A(lda,*)
    real    (c_double ) , intent(in)          :: B(ldb,*)
    real    (c_double ) , intent(out)         :: C(ldc,n)

  end function qmckl_distance_sq
end interface



! #+CALL: generate_f_interface(table=qmckl_distance_args,fname="qmckl_distance")

! #+RESULTS:

interface
  integer(qmckl_exit_code) function qmckl_distance &
      (context, transa, transb, m, n, A, lda, B, ldb, C, ldc) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (qmckl_context), intent(in)  , value :: context
    character(c_char  ) , intent(in)  , value :: transa
    character(c_char  ) , intent(in)  , value :: transb
    integer (c_int64_t) , intent(in)  , value :: m
    integer (c_int64_t) , intent(in)  , value :: n
    integer (c_int64_t) , intent(in)  , value :: lda
    integer (c_int64_t) , intent(in)  , value :: ldb
    integer (c_int64_t) , intent(in)  , value :: ldc
    real    (c_double ) , intent(in)          :: A(lda,*)
    real    (c_double ) , intent(in)          :: B(ldb,*)
    real    (c_double ) , intent(out)         :: C(ldc,n)

  end function qmckl_distance
end interface

! C interface                                                     :noexport:

!    #+CALL: generate_f_interface(table=qmckl_distance_rescaled_args,fname="qmckl_distance_rescaled")

!    #+RESULTS:

interface
  integer(qmckl_exit_code) function qmckl_distance_rescaled &
      (context, transa, transb, m, n, A, lda, B, ldb, C, ldc, rescale_factor_kappa) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
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

  end function qmckl_distance_rescaled
end interface



!  This function is more efficient when ~A~ and ~B~ are transposed.


! #+CALL: generate_f_interface(table=qmckl_distance_rescaled_gl_args,rettyp=get_value("FRetType"),fname=get_value("Name"))

! #+RESULTS:

interface
  integer(qmckl_exit_code) function qmckl_distance_rescaled_gl &
      (context, transa, transb, m, n, A, lda, B, ldb, C, ldc, rescale_factor_kappa) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
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
    real    (c_double ) , intent(out)         :: C(4,ldc,n)

  end function qmckl_distance_rescaled_gl
end interface
interface
  integer(qmckl_exit_code) function qmckl_set_electron_num(context, alpha, beta) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(in)  , value :: alpha
    integer (c_int64_t) , intent(in)  , value :: beta
  end function
end interface

interface
  integer(qmckl_exit_code) function qmckl_set_electron_coord(context, transp, walk_num, coord, size_max) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    character(c_char)   , intent(in)  , value :: transp
    integer (c_int64_t) , intent(in)  , value :: walk_num
    real(c_double)      , intent(in)          :: coord(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function
end interface

interface
  integer(qmckl_exit_code) function qmckl_get_electron_coord(context, transp, coord, size_max) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    character(c_char)   , intent(in)  , value :: transp
    real(c_double)      , intent(in)          :: coord(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function
end interface

interface
  integer(qmckl_exit_code) function qmckl_get_electron_ee_distance(context, distance, size_max) &
    bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double ) , intent(out)         :: distance(*)
    integer (c_int64_t) , intent(in)          :: size_max
  end function
end interface

interface
  integer(c_int32_t) function qmckl_get_electron_en_distance(context, distance, size_max) &
    bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double ) , intent(out)         :: distance(*)
    integer (c_int64_t) , intent(in)          :: size_max
  end function
end interface
interface
   subroutine qmckl_string_of_error (error, string) bind(C, name='qmckl_string_of_error_f')
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_exit_code), intent(in), value :: error
     character(c_char), intent(out) :: string(128)
   end subroutine qmckl_string_of_error
end interface

interface
   subroutine qmckl_last_error (context, string) bind(C, name='qmckl_last_error')
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (c_int64_t) , intent(in), value :: context
     character(c_char), intent(out) :: string(*)
   end subroutine qmckl_last_error
end interface

interface
   function qmckl_check (context, rc) bind(C, name='qmckl_check')
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer(qmckl_exit_code) :: qmckl_check
     integer (c_int64_t) , intent(in), value :: context
     integer(qmckl_exit_code), intent(in), value :: rc
   end function qmckl_check
end interface
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
   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_asymp_jasb_pderiv(context, &
        asymp_jasb_pderiv, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
     integer(c_int64_t), intent(in), value       :: size_max
     real(c_double),   intent(out)               :: asymp_jasb_pderiv(size_max)
   end function qmckl_get_jastrow_champ_asymp_jasb_pderiv
end interface

interface
   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_factor_ee_pderiv (context, &
        factor_ee_pderiv, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
     integer(c_int64_t), intent(in), value       :: size_max
     real(c_double),   intent(out)               :: factor_ee_pderiv(size_max)
   end function qmckl_get_jastrow_champ_factor_ee_pderiv
end interface

! Fortran interface


interface
   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_factor_ee_gl_pderiv (context, &
        factor_ee_gl_pderiv, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
     integer(c_int64_t), intent(in), value       :: size_max
     real(c_double),   intent(out)               :: factor_ee_gl_pderiv(size_max)
   end function qmckl_get_jastrow_champ_factor_ee_gl_pderiv
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
   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_asymp_jasa_pderiv(context, &
        asymp_jasa_pderiv, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
     integer(c_int64_t), intent(in), value       :: size_max
     real(c_double),   intent(out)               :: asymp_jasa_pderiv(size_max)
   end function qmckl_get_jastrow_champ_asymp_jasa_pderiv
end interface

interface
   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_factor_en_pderiv (context, &
        factor_en_pderiv, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
     integer(c_int64_t), intent(in), value       :: size_max
     real(c_double),   intent(out)               :: factor_en_pderiv(size_max)
   end function qmckl_get_jastrow_champ_factor_en_pderiv
end interface

! Fortran interface


interface
   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_factor_en_gl_pderiv (context, &
        factor_en_gl_pderiv, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
     integer(c_int64_t), intent(in), value       :: size_max
     real(c_double),   intent(out)               :: factor_en_gl_pderiv(size_max)
   end function qmckl_get_jastrow_champ_factor_en_gl_pderiv
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
   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_factor_een_pderiv (context, &
        factor_een_pderiv, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
     integer(c_int64_t), intent(in), value       :: size_max
     real(c_double),   intent(out)               :: factor_een_pderiv(size_max)
   end function qmckl_get_jastrow_champ_factor_een_pderiv
end interface

! Fortran interface


interface
   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_factor_een_gl_pderiv (context, &
        factor_een_gl_pderiv, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
     integer(c_int64_t), intent(in), value       :: size_max
     real(c_double),   intent(out)               :: factor_een_gl_pderiv(size_max)
   end function qmckl_get_jastrow_champ_factor_een_gl_pderiv
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
   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_single_tmp (context, &
        tmp, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
     integer(c_int64_t), intent(in), value       :: size_max
     real(c_double),   intent(out)               :: tmp(size_max)
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
   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_single_een_pderiv (context, &
        delta_een_pderiv, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
     integer(c_int64_t), intent(in), value       :: size_max
     real(c_double),   intent(out)               :: delta_een_pderiv(size_max)
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
   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_single_ee_pderiv (context, &
        delta_ee_pderiv, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
     integer(c_int64_t), intent(in), value       :: size_max
     real(c_double),   intent(out)               :: delta_ee_pderiv(size_max)
   end function qmckl_get_jastrow_champ_single_ee_pderiv
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
   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_single_en_pderiv (context, &
        delta_en_pderiv, size_max) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
     integer(c_int64_t), intent(in), value       :: size_max
     real(c_double),   intent(out)               :: delta_en_pderiv(size_max)
   end function qmckl_get_jastrow_champ_single_en_pderiv
end interface

interface
   integer(qmckl_exit_code) function qmckl_get_jastrow_champ_single_accept (context) bind(C)
     use, intrinsic :: iso_c_binding
     import
     implicit none
     integer (qmckl_context) , intent(in), value :: context
   end function qmckl_get_jastrow_champ_single_accept
end interface
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
interface
  integer(c_int32_t) function qmckl_get_nucleus_num(context, num) &
    bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(out)         :: num
  end function qmckl_get_nucleus_num
end interface

interface
  integer(c_int32_t) function qmckl_get_nucleus_charge(context, charge, size_max) &
    bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double)  , intent(out)         :: charge(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_get_nucleus_charge
end interface

interface
  integer(c_int32_t) function qmckl_get_nucleus_coord(context, transp, coord, size_max) &
    bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    character(c_char)   , intent(in)  , value :: transp
    real    (c_double)  , intent(out)         :: coord(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_get_nucleus_coord
end interface

interface
  integer(c_int32_t) function qmckl_set_nucleus_num(context, num) &
    bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(in)  , value :: num
  end function qmckl_set_nucleus_num
end interface

interface
  integer(c_int32_t) function qmckl_set_nucleus_charge(context, charge, size_max) &
    bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double)  , intent(in)          :: charge(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function
end interface

interface
  integer(c_int32_t) function qmckl_set_nucleus_coord(context, transp, coord, size_max) &
    bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    character(c_char)   , intent(in)  , value :: transp
    real    (c_double)  , intent(in)          :: coord(*)
    integer (c_int64_t) , intent(in)  , value :: size_max
  end function qmckl_set_nucleus_coord
end interface

interface
  integer(c_int32_t) function qmckl_get_nucleus_nn_distance(context, distance, size_max) &
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
  integer(c_int32_t) function qmckl_get_nucleus_repulsion(context, energy) &
    bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none
    integer (c_int64_t) , intent(in)  , value :: context
    real    (c_double ) , intent(out)         :: energy
  end function
end interface
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
interface
  integer(c_int32_t) function qmckl_get_point_num(context, num) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(out)         :: num
  end function
end interface

interface
  integer(c_int32_t) function qmckl_get_point(context, transp, coord, size_max) bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    character(c_char)   , intent(in)  , value :: transp
    real    (c_double ) , intent(out)         :: coord(*)
    integer (c_int64_t) , intent(in)          :: size_max
  end function
end interface

interface
  integer(c_int32_t) function qmckl_set_point(context, &
       transp, num, coord, size_max) bind(C)
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
! Fortran interfaces (exposed in qmckl_f.F90)
!    :PROPERTIES:
!    :Name:     qmckl_sm_naive
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

! #+CALL: generate_f_interface(table=qmckl_sm_naive_args,rettyp=get_value("FRetType"),fname=get_value("Name"))

! #+RESULTS:

interface
  integer(c_int32_t) function qmckl_sm_naive &
      (context, LDS, Dim, N_updates, Updates, Updates_index, breakdown, Slater_inv, determinant) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(in)  , value :: LDS
    integer (c_int64_t) , intent(in)  , value :: Dim
    integer (c_int64_t) , intent(in)  , value :: N_updates
    real    (c_double ) , intent(in)          :: Updates(N_updates*LDS)
    integer (c_int64_t) , intent(in)          :: Updates_index(N_updates)
    real    (c_double ) , intent(in)  , value :: breakdown
    real    (c_double ) , intent(inout)        :: Slater_inv(Dim*LDS)
    real    (c_double ) , intent(inout)        :: determinant

  end function qmckl_sm_naive
end interface



! #+CALL: generate_f_interface(table=qmckl_sm_naive_args,rettyp=get_value("FRetType"),fname="qmckl_sm_naive_hpc")

! #+RESULTS:

interface
  integer(c_int32_t) function qmckl_sm_naive_hpc &
      (context, LDS, Dim, N_updates, Updates, Updates_index, breakdown, Slater_inv, determinant) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(in)  , value :: LDS
    integer (c_int64_t) , intent(in)  , value :: Dim
    integer (c_int64_t) , intent(in)  , value :: N_updates
    real    (c_double ) , intent(in)          :: Updates(N_updates*LDS)
    integer (c_int64_t) , intent(in)          :: Updates_index(N_updates)
    real    (c_double ) , intent(in)  , value :: breakdown
    real    (c_double ) , intent(inout)        :: Slater_inv(Dim*LDS)
    real    (c_double ) , intent(inout)        :: determinant

  end function qmckl_sm_naive_hpc
end interface



! #+CALL: generate_f_interface(table=qmckl_sm_naive_args,rettyp=get_value("FRetType"),fname="qmckl_sm_naive_doc")

! #+RESULTS:

interface
  integer(c_int32_t) function qmckl_sm_naive_doc &
      (context, LDS, Dim, N_updates, Updates, Updates_index, breakdown, Slater_inv, determinant) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(in)  , value :: LDS
    integer (c_int64_t) , intent(in)  , value :: Dim
    integer (c_int64_t) , intent(in)  , value :: N_updates
    real    (c_double ) , intent(in)          :: Updates(N_updates*LDS)
    integer (c_int64_t) , intent(in)          :: Updates_index(N_updates)
    real    (c_double ) , intent(in)  , value :: breakdown
    real    (c_double ) , intent(inout)        :: Slater_inv(Dim*LDS)
    real    (c_double ) , intent(inout)        :: determinant

  end function qmckl_sm_naive_doc
end interface

! Fortran interfaces (exposed in qmckl_f.F90)
! :PROPERTIES:
! :Name:     qmckl_sm_splitting_core
! :CRetType: qmckl_exit_code
! :FRetType: qmckl_exit_code
! :END:

! #+CALL: generate_f_interface(table=qmckl_sm_splitting_core_args,rettyp=get_value("FRetType"),fname="qmckl_sm_splitting_core_hpc")

! #+RESULTS:

interface
  integer(c_int32_t) function qmckl_sm_splitting_core_hpc &
      (context, LDS, Dim, N_updates, Updates, Updates_index, breakdown, &
      Slater_inv, later_updates, later_index, later, determinant) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(in)  , value :: LDS
    integer (c_int64_t) , intent(in)  , value :: Dim
    integer (c_int64_t) , intent(in)  , value :: N_updates
    real    (c_double ) , intent(in)          :: Updates(N_updates*LDS)
    integer (c_int64_t) , intent(in)          :: Updates_index(N_updates)
    real    (c_double ) , intent(in)  , value :: breakdown
    real    (c_double ) , intent(inout)        :: Slater_inv(Dim*LDS)
    real    (c_double ) , intent(inout)        :: later_updates(N_updates*LDS)
    integer (c_int64_t) , intent(inout)        :: later_index(N_updates)
    integer (c_int64_t) , intent(inout)        :: later
    real    (c_double ) , intent(inout)        :: determinant

  end function qmckl_sm_splitting_core_hpc
end interface



! #+CALL: generate_f_interface(table=qmckl_sm_splitting_core_args,rettyp=get_value("FRetType"),fname="qmckl_sm_splitting_core_doc")

! #+RESULTS:

interface
  integer(c_int32_t) function qmckl_sm_splitting_core_doc &
      (context, LDS, Dim, N_updates, Updates, Updates_index, breakdown, &
      Slater_inv, later_updates, later_index, later, determinant) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(in)  , value :: LDS
    integer (c_int64_t) , intent(in)  , value :: Dim
    integer (c_int64_t) , intent(in)  , value :: N_updates
    real    (c_double ) , intent(in)          :: Updates(N_updates*LDS)
    integer (c_int64_t) , intent(in)          :: Updates_index(N_updates)
    real    (c_double ) , intent(in)  , value :: breakdown
    real    (c_double ) , intent(inout)        :: Slater_inv(Dim*LDS)
    real    (c_double ) , intent(inout)        :: later_updates(N_updates*LDS)
    integer (c_int64_t) , intent(inout)        :: later_index(N_updates)
    integer (c_int64_t) , intent(inout)        :: later
    real    (c_double ) , intent(inout)        :: determinant

  end function qmckl_sm_splitting_core_doc
end interface



! #+CALL: generate_f_interface(table=qmckl_sm_splitting_core_args,rettyp=get_value("FRetType"),fname=get_value("Name"))

! #+RESULTS:

interface
  integer(c_int32_t) function qmckl_sm_splitting_core &
      (context, LDS, Dim, N_updates, Updates, Updates_index, breakdown, &
      Slater_inv, later_updates, later_index, later, determinant) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(in)  , value :: LDS
    integer (c_int64_t) , intent(in)  , value :: Dim
    integer (c_int64_t) , intent(in)  , value :: N_updates
    real    (c_double ) , intent(in)          :: Updates(N_updates*LDS)
    integer (c_int64_t) , intent(in)          :: Updates_index(N_updates)
    real    (c_double ) , intent(in)  , value :: breakdown
    real    (c_double ) , intent(inout)        :: Slater_inv(Dim*LDS)
    real    (c_double ) , intent(inout)        :: later_updates(N_updates*LDS)
    integer (c_int64_t) , intent(inout)        :: later_index(N_updates)
    integer (c_int64_t) , intent(inout)        :: later
    real    (c_double ) , intent(inout)        :: determinant

  end function qmckl_sm_splitting_core
end interface

! Fortran interfaces (exposed in qmckl_f.F90)
! :PROPERTIES:
! :Name:     qmckl_woodbury_2x2
! :CRetType: qmckl_exit_code
! :FRetType: qmckl_exit_code
! :END:

! #+CALL: generate_f_interface(table=qmckl_woodbury_2x2_args,rettyp=get_value("FRetType"),fname=get_value("Name"))

! #+RESULTS:

interface
  integer(c_int32_t) function qmckl_woodbury_2x2 &
      (context, LDS, Dim, Updates, Updates_index, breakdown, Slater_inv, determinant) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(in)  , value :: LDS
    integer (c_int64_t) , intent(in)  , value :: Dim
    real    (c_double ) , intent(in)          :: Updates(2*Dim)
    integer (c_int64_t) , intent(in)          :: Updates_index(2)
    real    (c_double ) , intent(in)  , value :: breakdown
    real    (c_double ) , intent(inout)        :: Slater_inv(LDS*Dim)
    real    (c_double ) , intent(inout)        :: determinant

  end function qmckl_woodbury_2x2
end interface



! #+CALL: generate_f_interface(table=qmckl_woodbury_2x2_args,rettyp=get_value("FRetType"),fname="qmckl_woodbury_2x2_hpc")

! #+RESULTS:

interface
  integer(c_int32_t) function qmckl_woodbury_2x2_hpc &
      (context, LDS, Dim, Updates, Updates_index, breakdown, Slater_inv, determinant) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(in)  , value :: LDS
    integer (c_int64_t) , intent(in)  , value :: Dim
    real    (c_double ) , intent(in)          :: Updates(2*Dim)
    integer (c_int64_t) , intent(in)          :: Updates_index(2)
    real    (c_double ) , intent(in)  , value :: breakdown
    real    (c_double ) , intent(inout)        :: Slater_inv(LDS*Dim)
    real    (c_double ) , intent(inout)        :: determinant

  end function qmckl_woodbury_2x2_hpc
end interface



! #+CALL: generate_f_interface(table=qmckl_woodbury_2x2_args,rettyp=get_value("FRetType"),fname="qmckl_woodbury_2x2_doc")

! #+RESULTS:

interface
  integer(c_int32_t) function qmckl_woodbury_2x2_doc &
      (context, LDS, Dim, Updates, Updates_index, breakdown, Slater_inv, determinant) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(in)  , value :: LDS
    integer (c_int64_t) , intent(in)  , value :: Dim
    real    (c_double ) , intent(in)          :: Updates(2*Dim)
    integer (c_int64_t) , intent(in)          :: Updates_index(2)
    real    (c_double ) , intent(in)  , value :: breakdown
    real    (c_double ) , intent(inout)        :: Slater_inv(LDS*Dim)
    real    (c_double ) , intent(inout)        :: determinant

  end function qmckl_woodbury_2x2_doc
end interface

! Fortran interfaces (exposed in qmckl_f.F90)
! :PROPERTIES:
! :Name:     qmckl_woodbury_3x3
! :CRetType: qmckl_exit_code
! :FRetType: qmckl_exit_code
! :END:

! #+CALL: generate_f_interface(table=qmckl_woodbury_3x3_args,rettyp=get_value("FRetType"),fname=get_value("Name"))

! #+RESULTS:

interface
  integer(c_int32_t) function qmckl_woodbury_3x3 &
      (context, LDS, Dim, Updates, Updates_index, breakdown, Slater_inv, determinant) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(in)  , value :: LDS
    integer (c_int64_t) , intent(in)  , value :: Dim
    real    (c_double ) , intent(in)          :: Updates(3*Dim)
    integer (c_int64_t) , intent(in)          :: Updates_index(3)
    real    (c_double ) , intent(in)  , value :: breakdown
    real    (c_double ) , intent(inout)        :: Slater_inv(LDS*Dim)
    real    (c_double ) , intent(inout)        :: determinant

  end function qmckl_woodbury_3x3
end interface



! #+CALL: generate_f_interface(table=qmckl_woodbury_3x3_args,rettyp=get_value("FRetType"),fname="qmckl_woodbury_3x3_hpc")

! #+RESULTS:

interface
  integer(c_int32_t) function qmckl_woodbury_3x3_hpc &
      (context, LDS, Dim, Updates, Updates_index, breakdown, Slater_inv, determinant) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(in)  , value :: LDS
    integer (c_int64_t) , intent(in)  , value :: Dim
    real    (c_double ) , intent(in)          :: Updates(3*Dim)
    integer (c_int64_t) , intent(in)          :: Updates_index(3)
    real    (c_double ) , intent(in)  , value :: breakdown
    real    (c_double ) , intent(inout)        :: Slater_inv(LDS*Dim)
    real    (c_double ) , intent(inout)        :: determinant

  end function qmckl_woodbury_3x3_hpc
end interface



! #+CALL: generate_f_interface(table=qmckl_woodbury_3x3_args,rettyp=get_value("FRetType"),fname="qmckl_woodbury_3x3_doc")

! #+RESULTS:

interface
  integer(c_int32_t) function qmckl_woodbury_3x3_doc &
      (context, LDS, Dim, Updates, Updates_index, breakdown, Slater_inv, determinant) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(in)  , value :: LDS
    integer (c_int64_t) , intent(in)  , value :: Dim
    real    (c_double ) , intent(in)          :: Updates(3*Dim)
    integer (c_int64_t) , intent(in)          :: Updates_index(3)
    real    (c_double ) , intent(in)  , value :: breakdown
    real    (c_double ) , intent(inout)        :: Slater_inv(LDS*Dim)
    real    (c_double ) , intent(inout)        :: determinant

  end function qmckl_woodbury_3x3_doc
end interface

! Fortran interfaces (exposed in qmckl_f.F90)
!    :PROPERTIES:
!    :Name:     qmckl_sm_naive
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

! #+CALL: generate_f_interface(table=qmckl_sm_splitting_args,rettyp=get_value("FRetType"),fname="qmckl_sm_splitting_hpc")

! #+RESULTS:

interface
  integer(c_int32_t) function qmckl_sm_splitting_hpc &
      (context, LDS, Dim, N_updates, Updates, Updates_index, breakdown, Slater_inv, determinant) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(in)  , value :: LDS
    integer (c_int64_t) , intent(in)  , value :: Dim
    integer (c_int64_t) , intent(in)  , value :: N_updates
    real    (c_double ) , intent(in)          :: Updates(N_updates*LDS)
    integer (c_int64_t) , intent(in)          :: Updates_index(N_updates)
    real    (c_double ) , intent(in)  , value :: breakdown
    real    (c_double ) , intent(inout)        :: Slater_inv(Dim*LDS)
    real    (c_double ) , intent(inout)        :: determinant

  end function qmckl_sm_splitting_hpc
end interface



! #+CALL: generate_f_interface(table=qmckl_sm_splitting_args,rettyp=get_value("FRetType"),fname="qmckl_sm_splitting_doc")

! #+RESULTS:

interface
  integer(c_int32_t) function qmckl_sm_splitting_doc &
      (context, LDS, Dim, N_updates, Updates, Updates_index, breakdown, Slater_inv, determinant) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(in)  , value :: LDS
    integer (c_int64_t) , intent(in)  , value :: Dim
    integer (c_int64_t) , intent(in)  , value :: N_updates
    real    (c_double ) , intent(in)          :: Updates(N_updates*LDS)
    integer (c_int64_t) , intent(in)          :: Updates_index(N_updates)
    real    (c_double ) , intent(in)  , value :: breakdown
    real    (c_double ) , intent(inout)        :: Slater_inv(Dim*LDS)
    real    (c_double ) , intent(inout)        :: determinant

  end function qmckl_sm_splitting_doc
end interface



! #+CALL: generate_f_interface(table=qmckl_sm_splitting_args,rettyp=get_value("FRetType"),fname=get_value("Name"))

! #+RESULTS:

interface
  integer(c_int32_t) function qmckl_sm_splitting &
      (context, LDS, Dim, N_updates, Updates, Updates_index, breakdown, Slater_inv, determinant) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer (c_int64_t) , intent(in)  , value :: context
    integer (c_int64_t) , intent(in)  , value :: LDS
    integer (c_int64_t) , intent(in)  , value :: Dim
    integer (c_int64_t) , intent(in)  , value :: N_updates
    real    (c_double ) , intent(in)          :: Updates(N_updates*LDS)
    integer (c_int64_t) , intent(in)          :: Updates_index(N_updates)
    real    (c_double ) , intent(in)  , value :: breakdown
    real    (c_double ) , intent(inout)        :: Slater_inv(Dim*LDS)
    real    (c_double ) , intent(inout)        :: determinant

  end function qmckl_sm_splitting
end interface
interface
  integer(c_int32_t) function qmckl_trexio_read &
      (context, file_name, size_max) &
      bind(C)
    use, intrinsic :: iso_c_binding
    import
    implicit none

    integer  (c_int64_t) , intent(in)  , value :: context
    integer  (c_int64_t) , intent(in)  , value :: size_max
    character(c_char   ) , intent(in)          :: file_name(size_max)

  end function qmckl_trexio_read
end interface
end module qmckl
