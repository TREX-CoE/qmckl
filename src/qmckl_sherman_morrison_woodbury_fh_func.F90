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
