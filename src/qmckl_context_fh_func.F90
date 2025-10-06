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

! interface
!    integer (qmckl_context) function qmckl_context_copy(context) bind(C)
!      use, intrinsic :: iso_c_binding
!      import
!      integer (qmckl_context), intent(in), value :: context
!    end function qmckl_context_copy
! end interface

interface
   integer (qmckl_exit_code) function qmckl_context_destroy(context) bind(C)
     use, intrinsic :: iso_c_binding
     import
     integer (qmckl_context), intent(in), value :: context
   end function qmckl_context_destroy
end interface
