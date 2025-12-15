#ifndef QMCKL_DETERMINANT_HPT
#define QMCKL_DETERMINANT_HPT

#include <stdbool.h>

/* Data structure */


typedef struct qmckl_determinant_struct {
  char     type;
  int64_t  det_num_alpha;
  int64_t  det_num_beta ;
  int64_t  up_num;
  int64_t  down_num;
  int64_t* mo_index_alpha;
  int64_t* mo_index_beta;

  double  * det_value_alpha;
  double  * det_value_beta;
  double  * det_vgl_alpha;
  double  * det_adj_matrix_alpha;
  double  * det_inv_matrix_alpha;
  double  * det_vgl_beta;
  double  * det_adj_matrix_beta;
  double  * det_inv_matrix_beta;
  uint64_t  det_value_alpha_date;
  uint64_t  det_vgl_alpha_date;
  uint64_t  det_adj_matrix_alpha_date;
  uint64_t  det_inv_matrix_alpha_date;
  uint64_t  det_value_beta_date;
  uint64_t  det_vgl_beta_date;
  uint64_t  det_adj_matrix_beta_date;
  uint64_t  det_inv_matrix_beta_date;

  int32_t   uninitialized;
  bool      provided;
} qmckl_determinant_struct;

#endif
