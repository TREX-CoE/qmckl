#ifndef QMCKL_JASTROW_CHAMP_HPT
#define QMCKL_JASTROW_CHAMP_HPT
#include <stdbool.h>

/* Data structure */


typedef struct qmckl_jastrow_champ_struct{
  int64_t * restrict lkpm_combined_index;
  int64_t * restrict type_nucl_vector;
  double  * restrict asymp_jasa;
  double  * restrict asymp_jasa_pderiv;
  double             asymp_jasb[2];
  double  * restrict asymp_jasb_pderiv;
  double  * restrict a_vector;
  double  * restrict b_vector;
  double  * restrict c_vector;
  double  * restrict c_vector_full;
  double  * restrict dtmp_c;
  double  * restrict ee_distance_rescaled;
  double  * restrict ee_distance_rescaled_gl;
  double  * restrict een_rescaled_e;
  double  * restrict een_rescaled_e_gl;
  double  * restrict een_rescaled_n;
  double  * restrict een_rescaled_n_gl;
  double  * restrict en_distance_rescaled;
  double  * restrict en_distance_rescaled_gl;
  double  * restrict factor_ee;
  double  * restrict factor_ee_gl;
  double  * restrict factor_ee_pderiv;
  double  * restrict factor_ee_gl_pderiv;
  double  * restrict factor_een;
  double  * restrict factor_een_gl;
  double  * restrict factor_een_grad;
  double  * restrict factor_een_pderiv;
  double  * restrict factor_een_gl_pderiv;
  double  * restrict factor_en;
  double  * restrict factor_en_gl;
  double  * restrict factor_en_pderiv;
  double  * restrict factor_en_gl_pderiv;
  double  * restrict rescale_factor_en;
  double  * restrict tmp_c;
  double  * restrict value;
  double  * restrict gl;
  double  * restrict grad;
  int64_t   aord_num;
  int64_t   bord_num;
  int64_t   cord_num;
  int64_t   dim_c_vector;
  int64_t   type_nucl_num;
  uint64_t  asymp_jasa_date;
  uint64_t  asymp_jasa_pderiv_date;
  uint64_t  asymp_jasb_date;
  uint64_t  asymp_jasb_pderiv_date;
  uint64_t  c_vector_full_date;
  uint64_t  dim_c_vector_date;
  uint64_t  dtmp_c_date;
  uint64_t  ee_distance_rescaled_date;
  uint64_t  ee_distance_rescaled_gl_date;
  uint64_t  een_rescaled_e_date;
  uint64_t  een_rescaled_e_gl_date;
  uint64_t  een_rescaled_n_date;
  uint64_t  een_rescaled_n_gl_date;
  uint64_t  en_distance_rescaled_date;
  uint64_t  en_distance_rescaled_gl_date;
  uint64_t  factor_ee_date;
  uint64_t  factor_ee_gl_date;
  uint64_t  factor_ee_pderiv_date;
  uint64_t  factor_ee_gl_pderiv_date;
  uint64_t  factor_een_date;
  uint64_t  factor_een_gl_date;
  uint64_t  factor_een_grad_date;
  uint64_t  factor_een_pderiv_date;
  uint64_t  factor_een_gl_pderiv_date;
  uint64_t  factor_en_date;
  uint64_t  factor_en_gl_date;
  uint64_t  factor_en_pderiv_date;
  uint64_t  factor_en_gl_pderiv_date;
  uint64_t  lkpm_combined_index_date;
  uint64_t  tmp_c_date;
  uint64_t  value_date;
  uint64_t  gl_date;
  uint64_t  grad_date;
  double    rescale_factor_ee;
  int32_t   uninitialized;
  int32_t   spin_independent;
  bool      provided;

} qmckl_jastrow_champ_struct;

#endif
