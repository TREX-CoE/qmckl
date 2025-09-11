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
