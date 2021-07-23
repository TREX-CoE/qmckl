#include <stdbool.h>
#include <stdio.h>

#ifdef VFC_CI
#include <vfc_probes.h>
vfc_probes probes;
#endif


// Wrappers to Verificarlo functions

#ifdef VFC_CI
void __attribute__((constructor)) qmckl_init_probes(){
	probes = vfc_init_probes();
}
#endif


bool qmckl_probe(
    char * testName,
    char * varName,
    double value
) {
#ifdef VFC_CI
    return vfc_probe(&probes, testName, varName, value);
#else
	return false;
#endif
}


bool qmckl_probe_check(
    char * testName,
    char * varName,
    double value,
    double expectedValue,
    double accuracyTarget
) {
#ifdef VFC_CI
    return vfc_probe_check(&probes, testName, varName, value, accuracyTarget);
#else
    return !(value == expectedValue);
#endif
}


bool qmckl_probe_check_relative (
    char * testName,
    char * varName,
    double value,
    double expectedValue,
    double accuracyTarget
) {
#ifdef VFC_CI
    return vfc_probe_check_relative(&probes, testName, varName, value, accuracyTarget);
#else
    return !(value <= expectedValue + accuracyTarget || value >= expectedValue - accuracyTarget);
#endif
}


void __attribute__((destructor)) qmckl_dump_probes(){
#ifdef VFC_CI
    vfc_dump_probes(&probes);
#endif
}


// Fortran wrappers

bool qmckl_probe_f(
    char * testName,
    char * varName,
    double * value
) {
    return qmckl_probe(testName, varName, *value);
}


bool qmckl_probe_check_f(
    char * testName,
    char * varName,
    double * value,
    double * expectedValue,
    double * accuracyTarget
) {
    return qmckl_probe_check(
        testName, varName,
        *value, *expectedValue, *accuracyTarget
    );
}


bool qmckl_probe_check_relative_f(
    char * testName,
    char * varName,
    double * value,
    double * expectedValue,
    double * accuracyTarget
) {
    return qmckl_probe_check_relative(
        testName, varName,
        *value, *expectedValue, *accuracyTarget
    );
}
