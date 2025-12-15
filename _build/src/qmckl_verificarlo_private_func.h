#ifndef QMCKL_VERIFICARLO_HPT
#define QMCKL_VERIFICARLO_HPT

#include <stdbool.h>


#ifdef VFC_CI
#include <vfc_probes.h>
extern vfc_probes * probes;
#endif

void qmckl_init_probes();

bool qmckl_probe(
    char * testName,
    char * varName,
    double value
);

bool qmckl_probe_check(
    char * testName,
    char * varName,
    double value,
    double expectedValue,
    double accuracyTarget
);

bool qmckl_probe_check_relative(
    char * testName,
    char * varName,
    double value,
    double expectedValue,
    double accuracyTarget
);

void qmckl_dump_probes();

bool qmckl_probe_f(
    char * testName,
    char * varName,
    double * value
);

bool qmckl_probe_check_f(
    char * testName,
    char * varName,
    double * value,
    double * expectedValue,
    double * accuracyTarget
);


bool qmckl_probe_check_relative_f(
    char * testName,
    char * varName,
    double * value,
    double * expectedValue,
    double * accuracyTarget
);

/* [[file:../../org/qmckl_verificarlo.org::*End of files][End of files:2]] */
#endif
/* End of files:2 ends here */
