module qmckl_probes_f
    interface
        logical(c_bool) function qmckl_probe_c &
             (testName, varName, val) &
             bind(C, name="qmckl_probe")

            use, intrinsic :: iso_c_binding
            import
            implicit none

            character(C_CHAR), dimension(*) :: testName
            character(C_CHAR), dimension(*) :: varName

            real(C_DOUBLE), value :: val
        end function qmckl_probe_c

        logical(c_bool) function qmckl_probe_check_c &
             (testName, varName, val, expectedValue, accuracyTarget) &
             bind(C, name="qmckl_probe_check")

            use, intrinsic :: iso_c_binding
            import
            implicit none

            character(C_CHAR), dimension(*) :: testName
            character(C_CHAR), dimension(*) :: varName

            real(C_DOUBLE), value :: val
            real(C_DOUBLE), value :: expectedValue
            real(C_DOUBLE), value :: accuracyTarget
        end function qmckl_probe_check_c

        logical(c_bool) function qmckl_probe_check_relative_c &
             (testName, varName, val, expectedValue, accuracyTarget) &
             bind(C, name="qmckl_probe_check_relative")

            use, intrinsic :: iso_c_binding
            import
            implicit none

            character(C_CHAR), dimension(*) :: testName
            character(C_CHAR), dimension(*) :: varName

            real(C_DOUBLE), value :: val
            real(C_DOUBLE), value :: expectedValue
            real(C_DOUBLE), value :: accuracyTarget
        end function qmckl_probe_check_relative_c
    end interface

    logical function qmckl_probe (testName, varName, val)
        implicit none
        character, dimension(*) :: testName
        character, dimension(*) :: varName
        double precision        :: val
        return qmckl_probe_c(testName//C_NULL_CHAR, varName//C_NULL_CHAR, val)
    end function qmckl_probe

    logical function qmckl_probe_check
         (testName, varName, val, expectedValue, accuracyTarget)
        implicit none
        character, dimension(*) :: testName
        character, dimension(*) :: varName
        double precision        :: val
        double precision        :: expectedValue
        double precision        :: accuracyTarget
        return qmckl_probe_check_c(testName//C_NULL_CHAR, varName//C_NULL_CHAR, &
            val, expectedValue, accuracyTarget)
    end function qmckl_probe_check

    logical function qmckl_probe_check_relative &
         (testName, varName, val, expectedValue, accuracyTarget)
        implicit none
        character, dimension(*) :: testName
        character, dimension(*) :: varName
        double precision        :: val
        double precision        :: expectedValue
        double precision        :: accuracyTarget
        return qmckl_probe_check_relative_c(testName//C_NULL_CHAR, varName//C_NULL_CHAR, &
            val, expectedValue, accuracyTarget)
    end function qmckl_probe_check_relative
end module qmckl_probes_f
