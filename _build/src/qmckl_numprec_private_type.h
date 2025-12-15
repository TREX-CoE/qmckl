#ifndef QMCKL_NUMPREC_HPT
#define QMCKL_NUMPREC_HPT

#ifdef HAVE_STDINT_H
#include <stdint.h>
#elif HAVE_INTTYPES_H
#include <inttypes.h>
#endif


/* :end: */


typedef struct qmckl_numprec_struct {
  uint32_t  precision;
  uint32_t  range;
} qmckl_numprec_struct;

/* [[file:../../org/qmckl_numprec.org::*End of files][End of files:1]] */
#endif
/* End of files:1 ends here */
