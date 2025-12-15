#ifndef QMCKL_MO_HPT
#define QMCKL_MO_HPT

#include <stdbool.h>
#include "qmckl_blas_private_type.h"

/* Data structure */


typedef struct qmckl_mo_basis_struct {
  int64_t  mo_num;
  double * restrict coefficient;
  double * restrict coefficient_t;
  double * restrict r_cusp;

  double * restrict mo_vgl;
  double * restrict mo_value;
  qmckl_tensor cusp_param;

  uint64_t  mo_vgl_date;
  uint64_t  mo_value_date;

  int32_t   uninitialized;
  bool      provided;
} qmckl_mo_basis_struct;

#endif
