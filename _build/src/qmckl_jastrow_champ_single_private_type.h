#ifndef QMCKL_JASTROW_CHAMP_SINGLE_HPT
#define QMCKL_JASTROW_CHAMP_SINGLE_HPT
#include <stdbool.h>
#include "qmckl_blas_private_type.h"

/* Data structure */


typedef struct qmckl_jastrow_champ_single_struct{
  int64_t      num;
  uint64_t     date;
  qmckl_matrix coord;
  double  *    een_rescaled_single_e;
  uint64_t     een_rescaled_single_e_date;
  uint64_t     een_rescaled_single_e_maxsize;
  double  *    een_rescaled_single_n;
  uint64_t     een_rescaled_single_n_date;
  uint64_t     een_rescaled_single_n_maxsize;
  double*      single_ee_distance;
  uint64_t     single_ee_distance_date;
  uint64_t     single_ee_distance_maxsize;
  double*      single_en_distance;
  uint64_t     single_en_distance_date;
  uint64_t     single_en_distance_maxsize;
  double*      delta_een;
  void*        delta_een_hpc_function;
  uint64_t     delta_een_date;
  uint64_t     delta_een_maxsize;
  double*      delta_een_pderiv;
  uint64_t     delta_een_pderiv_date;
  uint64_t     delta_een_pderiv_maxsize;
  double*      ee_rescaled_single;
  uint64_t     ee_rescaled_single_date;
  uint64_t     ee_rescaled_single_maxsize;
  double*      en_rescaled_single;
  uint64_t     en_rescaled_single_date;
  uint64_t     en_rescaled_single_maxsize;
  double*      delta_en;
  uint64_t     delta_en_date;
  uint64_t     delta_en_maxsize;
  double*      delta_en_pderiv;
  uint64_t     delta_en_pderiv_date;
  uint64_t     delta_en_pderiv_maxsize;
  double*      delta_ee;
  uint64_t     delta_ee_date;
  uint64_t     delta_ee_maxsize;
  double*      delta_ee_pderiv;
  uint64_t     delta_ee_pderiv_date;
  uint64_t     delta_ee_pderiv_maxsize;
  double  *    een_rescaled_single_e_gl;
  uint64_t     een_rescaled_single_e_gl_date;
  uint64_t     een_rescaled_single_e_gl_maxsize;
  double  *    een_rescaled_single_n_gl;
  uint64_t     een_rescaled_single_n_gl_date;
  uint64_t     een_rescaled_single_n_gl_maxsize;
  double*      delta_een_gl;
  uint64_t     delta_een_gl_date;
  uint64_t     delta_een_gl_maxsize;
  double*      delta_een_g;
  uint64_t     delta_een_g_date;
  uint64_t     delta_een_g_maxsize;
  double*      ee_rescaled_single_gl;
  uint64_t     ee_rescaled_single_gl_date;
  uint64_t     ee_rescaled_single_gl_maxsize;
  double*      en_rescaled_single_gl;
  uint64_t     en_rescaled_single_gl_date;
  uint64_t     en_rescaled_single_gl_maxsize;
  double*      delta_en_gl;
  uint64_t     delta_en_gl_date;
  uint64_t     delta_en_gl_maxsize;
  double*      delta_ee_gl;
  uint64_t     delta_ee_gl_date;
  uint64_t     delta_ee_gl_maxsize;
  double*      tmp;
  uint64_t     tmp_date;
  uint64_t     tmp_maxsize;

} qmckl_jastrow_champ_single_struct;

#endif
