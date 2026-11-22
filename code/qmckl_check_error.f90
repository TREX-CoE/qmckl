subroutine qmckl_check_error(rc, message)
  use qmckl
  implicit none
  integer(qmckl_exit_code), intent(in) :: rc
  character(len=*)        , intent(in) :: message
  character(len=128)                   :: str_buffer
  if (rc /= QMCKL_SUCCESS) then
     print *, message
     call qmckl_string_of_error(rc, str_buffer)
     print *, str_buffer
     call exit(rc)
  end if
end subroutine qmckl_check_error
