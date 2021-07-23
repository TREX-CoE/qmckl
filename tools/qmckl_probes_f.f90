module qmckl_probes_f
    interface
        logical(c_bool) function qmckl_probe &
             (testName, varName, val) &
             bind(C, name="qmckl_probe_f")

            use, intrinsic :: iso_c_binding
            import
            implicit none

            character(C_CHAR), dimension(*) :: testName
            character(C_CHAR), dimension(*) :: varName

            real(C_DOUBLE) :: val
        end function qmckl_probe

        logical(c_bool) function qmckl_probe_check &
             (testName, varName, val, expectedValue, accuracyTarget) &
             bind(C, name="qmckl_probe_check_f")

            use, intrinsic :: iso_c_binding
            import
            implicit none

            character(C_CHAR), dimension(*) :: testName
            character(C_CHAR), dimension(*) :: varName

            real(C_DOUBLE) :: val
            real(C_DOUBLE) :: expectedValue
            real(C_DOUBLE) :: accuracyTarget
        end function qmckl_probe_check

        logical(c_bool) function qmckl_probe_check_relative &
             (testName, varName, val, expectedValue, accuracyTarget) &
             bind(C, name="qmckl_probe_check_relative_f")

            use, intrinsic :: iso_c_binding
            import
            implicit none

            character(C_CHAR), dimension(*) :: testName
            character(C_CHAR), dimension(*) :: varName

            real(C_DOUBLE) :: val
            real(C_DOUBLE) :: expectedValue
            real(C_DOUBLE) :: accuracyTarget
        end function qmckl_probe_check_relative
    end interface
end module qmckl_probes_f
