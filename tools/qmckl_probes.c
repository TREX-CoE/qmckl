#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>

#ifdef VFC_CI
#include <vfc_probes.h>
vfc_probes probes;
#endif


// QMCkl is a wrapper to Verificarlo's vfc_probes system. The goal of QMCkl 
// probes isto simplify the use of vfc_probes, and to provide functions that 
// can be called either wit or without vfc_ci support by using #ifndef 
// statements : 
//
// - when vfc_ci is disabled, qmckl_probes functions will either return false 
// (no error) or perform a check based on a reference value
// - when vfc_ci is enabled, qmckl_probe functions will simply encapsulate 
// calls to vfc_probe 
//
// Moreover, one does not have to worry about the life cycle of the probes 
// structure, as it is automatically created, dumped and freed by this wrapper.
//
// vfc_ci support can be enabled by using the following configure command : 
// QMCKL_DEVEL=1 ./configure --prefix=$PWD/_install --enable-silent-rules 
// --enable-maintainer-mode CC=verificarlo-f FC=verificarlo-f --host=x86_64
// 
// Finally, this wrapper also comes with a Fortran interface (in its dedicated 
// file).
//
// To learn more about Verificarlo CI : 
// https://github.com/verificarlo/verificarlo/blob/master/doc/06-Postprocessing.md#verificarlo-ci


// Automatically initialize the vfc_probe object if VFC_CI is defined
#ifdef VFC_CI
void __attribute__((constructor)) qmckl_init_probes(){
	probes = vfc_init_probes();
}
#endif


// Standard probe, without check
// - if VFC_CI is defined, place a standard probe
// - if VFC_CI is undefined, return false (no error)
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


// Probe with absolute check
// - if VFC_CI is defined, place a probe with an absolute check
// - if VFC_CI is undefined, perform an absolute check based on target value
// and accuracy
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
    return !(abs(value - expectedValue) < accuracyTarget);
#endif
}


// Probe with relative check
// - if VFC_CI is defined, place a probe with a relative check
// - if VFC_CI is undefined, perform a relative check based on target value
// and accuracy
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
    return !(abs(value - expectedValue) / abs(expectedValue) < accuracyTarget);
#endif
}


// Automatically delete and dump the vfc_probe object if VFC_CI is defined
#ifdef VFC_CI
void __attribute__((destructor)) qmckl_dump_probes(){
    vfc_dump_probes(&probes);
}
#endif

