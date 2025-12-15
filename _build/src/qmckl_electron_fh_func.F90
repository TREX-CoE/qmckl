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
