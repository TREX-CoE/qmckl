#ifndef QMCKL_AO_HPT
#define QMCKL_AO_HPT

#include <stdbool.h>
#include <stdio.h>
#include "qmckl_blas_private_type.h"

/* Data structure                                                 :noexport: */


typedef struct qmckl_ao_basis_struct {
  int64_t            shell_num;
  int64_t            prim_num;
  int64_t            ao_num;
  int64_t * restrict nucleus_index;
  int64_t * restrict nucleus_shell_num;
  int32_t * restrict r_power;
  int32_t * restrict shell_ang_mom;
  int64_t * restrict shell_prim_num;
  int64_t * restrict shell_prim_index;
  double  * restrict shell_factor;
  double  * restrict exponent;
  double  * restrict coefficient;
  double  * restrict prim_factor;
  double  * restrict ao_factor;
  int64_t * restrict ao_nucl;
  int32_t * restrict ao_ang_mom;

  int64_t * restrict nucleus_prim_index;
  double  * restrict coefficient_normalized;
  int32_t * restrict nucleus_max_ang_mom;
  double  * restrict nucleus_range;
  double  * restrict primitive_vgl;
  uint64_t           primitive_vgl_date;
  double  * restrict shell_vgl;
  uint64_t           shell_vgl_date;
  uint64_t           shell_vgl_maxsize;
  double  * restrict ao_vgl;
  uint64_t           ao_vgl_date;
  uint64_t           ao_vgl_maxsize;
  double  * restrict ao_value;
  uint64_t           ao_value_date;
  uint64_t           ao_value_maxsize;

  double * restrict  shell_hessian;
  uint64_t           shell_hessian_date;

  double * restrict  ao_hessian;
  uint64_t           ao_hessian_date;

  int32_t            uninitialized;
  bool               provided;
  bool               ao_cartesian;
  char               type;
#ifdef HAVE_HPC
  /* HPC specific data structures */
  int32_t*  restrict   prim_num_per_nucleus;
  qmckl_tensor         coef_per_nucleus;
  qmckl_matrix         expo_per_nucleus;
#endif
} qmckl_ao_basis_struct;

#endif
