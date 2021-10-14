#include <stdbool.h>


#ifdef VFC_CI
#include <vfc_probes.h>
extern vfc_probes * probes;
#endif

// Wrappers to Verificarlo functions

#ifdef VFC_CI
void qmckl_init_probes() __attribute__((constructor));
#endif

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

#ifdef VFC_CI
void qmckl_dump_probes() __attribute__((destructor));
#endif


// Fortran wrappers

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
