#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_STDINT_H
#include <stdint.h>
#elif HAVE_INTTYPES_H
#include <inttypes.h>
#endif

#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>

#include "qmckl.h"
#include "qmckl_context_private_type.h"
#include "qmckl_memory_private_type.h"
#include "qmckl_memory_private_func.h"
#include "qmckl_ao_private_type.h"
#include "qmckl_ao_private_func.h"
#include "qmckl_mo_private_type.h"
#include "qmckl_mo_private_func.h"
#include "qmckl_determinant_private_type.h"
#include "qmckl_determinant_private_func.h"

qmckl_exit_code qmckl_init_determinant(qmckl_context context) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return false;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  ctx->det.uninitialized = (1 << 5) - 1;

  return QMCKL_SUCCESS;
}

bool qmckl_determinant_provided(const qmckl_context context) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return false;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  return ctx->det.provided;
}

char qmckl_get_determinant_type (const qmckl_context context) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return (char) 0;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1;

  if ( (ctx->det.uninitialized & mask) != 0) {
    return (char) 0;
  }

  assert (ctx->det.type != (char) 0);
  return ctx->det.type;
}

int64_t qmckl_get_determinant_det_num_alpha (const qmckl_context context) {
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return (int64_t) 0;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 1;

  if ( (ctx->det.uninitialized & mask) != 0) {
    return (int64_t) 0;
  }

  assert (ctx->det.det_num_alpha > (int64_t) 0);
  return ctx->det.det_num_alpha;
}

int64_t qmckl_get_determinant_det_num_beta (const qmckl_context context) {
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return (int64_t) 0;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 2;

  if ( (ctx->det.uninitialized & mask) != 0) {
    return (int64_t) 0;
  }

  assert (ctx->det.det_num_beta > (int64_t) 0);
  return ctx->det.det_num_beta;
}

int64_t* qmckl_get_determinant_mo_index_alpha (const qmckl_context context) {
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return NULL;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 3;

  if ( (ctx->det.uninitialized & mask) != 0) {
    return NULL;
  }

  assert (ctx->det.mo_index_alpha != NULL);
  return ctx->det.mo_index_alpha;
}

int64_t* qmckl_get_determinant_mo_index_beta (const qmckl_context context) {
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return NULL;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 4;

  if ( (ctx->det.uninitialized & mask) != 0) {
    return NULL;
  }

  assert (ctx->det.mo_index_beta != NULL);
  return ctx->det.mo_index_beta;
}

qmckl_exit_code qmckl_set_determinant_type(qmckl_context context, const char t) {
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
   }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;

  if (t != 'G' && t != 'S') {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_set_determinant_type",
                           NULL);
  }

  int32_t mask = 1;
  ctx->det.type = t;

  ctx->det.uninitialized &= ~mask;
  ctx->det.provided = (ctx->det.uninitialized == 0);
  if (ctx->det.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_determinant(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
   }
  
  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_set_determinant_det_num_alpha(qmckl_context context, const int64_t det_num_alpha) {
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
   }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;

  if (det_num_alpha <= 0) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_set_determinant_det_num_alpha",
                           "det_num_alpha <= 0");
  }

  int32_t mask = 1 << 1;
  ctx->det.det_num_alpha = det_num_alpha;

  ctx->det.uninitialized &= ~mask;
  ctx->det.provided = (ctx->det.uninitialized == 0);
  if (ctx->det.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_determinant(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
   }
  
  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_set_determinant_det_num_beta(qmckl_context context, const int64_t det_num_beta) {
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
   }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;

  if (det_num_beta <= 0) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_set_determinant_det_num_beta",
                           "det_num_beta <= 0");
  }

  int32_t mask = 1 << 2;
  ctx->det.det_num_beta = det_num_beta;

  ctx->det.uninitialized &= ~mask;
  ctx->det.provided = (ctx->det.uninitialized == 0);
  if (ctx->det.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_determinant(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
   }
  
  return QMCKL_SUCCESS;
}

qmckl_exit_code  qmckl_set_determinant_mo_index_alpha(qmckl_context context, const int64_t* mo_index_alpha, const int64_t size_max) {
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
   }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;

  int32_t mask = 1 << 3;

  if (ctx->det.mo_index_alpha != NULL) {
    qmckl_exit_code rc = qmckl_free(context, ctx->det.mo_index_alpha);
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context, rc,
                             "qmckl_set_determinant_mo_index_alpha",
                             NULL);
    }
  }

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = ctx->electron.walker.num * ctx->det.det_num_alpha *
                  ctx->electron.up_num * sizeof(int64_t);
  int64_t* new_array = (int64_t*) qmckl_malloc(context, mem_info);
  if (new_array == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_set_determinant_mo_index_alpha",
                           NULL);
  }

  if (size_max * sizeof(int64_t) < mem_info.size) {
      return qmckl_failwith( context,
                             QMCKL_INVALID_ARG_3,
                             "qmckl_set_determinant_mo_index_alpha",
                             "input array too small");
  }

  memcpy(new_array, mo_index_alpha, mem_info.size);

  ctx->det.mo_index_alpha = new_array;

  ctx->det.uninitialized &= ~mask;
  ctx->det.provided = (ctx->det.uninitialized == 0);
  if (ctx->det.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_determinant(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
   }
  
  return QMCKL_SUCCESS;
}

qmckl_exit_code  qmckl_set_determinant_mo_index_beta(qmckl_context context, const int64_t* mo_index_beta, const int64_t size_max) {
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
   }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;

  int32_t mask = 1 << 4;

  if (ctx->det.mo_index_beta != NULL) {
    qmckl_exit_code rc = qmckl_free(context, ctx->det.mo_index_beta);
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context, rc,
                             "qmckl_set_determinant_mo_index_beta",
                             NULL);
    }
  }

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = ctx->electron.walker.num * ctx->det.det_num_beta *
                  ctx->electron.down_num * sizeof(int64_t);
  int64_t* new_array = (int64_t*) qmckl_malloc(context, mem_info);
  if (new_array == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_set_determinant_mo_index_beta",
                           NULL);
  }

  if (size_max * sizeof(int64_t) < mem_info.size) {
      return qmckl_failwith( context,
                             QMCKL_INVALID_ARG_3,
                             "qmckl_set_determinant_mo_index_beta",
                             "input array too small");
  }

  memcpy(new_array, mo_index_beta, mem_info.size);

  ctx->det.mo_index_beta = new_array;

  ctx->det.uninitialized &= ~mask;
  ctx->det.provided = (ctx->det.uninitialized == 0);
  if (ctx->det.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_determinant(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
   }
  
  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_finalize_determinant(qmckl_context context) {
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_finalize_determinant",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  qmckl_exit_code rc;
  rc = qmckl_provide_det_vgl_alpha(context);
  if (rc != QMCKL_SUCCESS) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_finalize_determinant",
                            NULL);
  }
  rc = qmckl_provide_det_vgl_beta(context);
  if (rc != QMCKL_SUCCESS) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_finalize_determinant",
                            NULL);
  }
  rc = qmckl_provide_det_inv_matrix_alpha(context);
  if (rc != QMCKL_SUCCESS) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_finalize_determinant",
                            NULL);
  }
  rc = qmckl_provide_det_inv_matrix_beta(context);
  if (rc != QMCKL_SUCCESS) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_finalize_determinant",
                            NULL);
  }
  return rc;
}

qmckl_exit_code qmckl_get_det_vgl_alpha(qmckl_context context, double * const det_vgl_alpha) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_ao_basis_ao_vgl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_mo_basis_mo_vgl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_det_vgl_alpha(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  size_t sze = 5 * ctx->det.det_num_alpha * ctx->electron.walker.num *
               ctx->electron.up_num * ctx->electron.up_num * sizeof(double);
  memcpy(det_vgl_alpha, ctx->det.det_vgl_alpha, sze);

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_get_det_vgl_beta(qmckl_context context, double * const det_vgl_beta) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_ao_basis_ao_vgl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_mo_basis_mo_vgl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_det_vgl_beta(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  size_t sze = 5 * ctx->det.det_num_beta * ctx->electron.walker.num *
               ctx->electron.down_num * ctx->electron.down_num * sizeof(double);
  memcpy(det_vgl_beta, ctx->det.det_vgl_beta, sze);

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_det_vgl_alpha(qmckl_context context) {

  qmckl_exit_code rc;
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if(!(ctx->nucleus.provided)) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_electron",
                           NULL);
  }

  if(!(ctx->electron.provided)) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_electron",
                           NULL);
  }

  if (!ctx->ao_basis.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_ao_basis",
                           NULL);
  }

  if (!ctx->mo_basis.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_mo_basis",
                           NULL);
  }
  rc = qmckl_provide_mo_basis_mo_vgl(context);
  if (rc != QMCKL_SUCCESS) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_mo_basis",
                           NULL);
  }


  if (!ctx->det.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_determinant",
                           NULL);
  }

  /* Compute if necessary */
  if (ctx->electron.walker.point.date > ctx->det.det_vgl_alpha_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      free(ctx->det.det_vgl_alpha);
      ctx->det.det_vgl_alpha = NULL;
    }
    
    /* Allocate array */
    if (ctx->det.det_vgl_alpha == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = 5 * ctx->electron.walker.num * ctx->det.det_num_alpha *
                      ctx->electron.up_num * ctx->electron.up_num * sizeof(double);
      double* det_vgl_alpha = (double*) qmckl_malloc(context, mem_info);

      if (det_vgl_alpha == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_det_vgl_alpha",
                               NULL);
      }
      ctx->det.det_vgl_alpha = det_vgl_alpha;
    }

    if (ctx->det.type == 'G') {
      rc = qmckl_compute_det_vgl_alpha(context,
                                 ctx->det.det_num_alpha,
                                 ctx->electron.walker.num,
                                 ctx->electron.up_num,
                                 ctx->electron.down_num,
                                 ctx->electron.num,
                                 ctx->det.mo_index_alpha,
                                 ctx->mo_basis.mo_num,
                                 ctx->mo_basis.mo_vgl,
                                 ctx->det.det_vgl_alpha);
    } else {
      return qmckl_failwith( context,
                             QMCKL_FAILURE,
                             "compute_det_vgl_alpha",
                             "Not yet implemented");
    }
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->det.det_vgl_alpha_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_det_vgl_beta(qmckl_context context) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if(!(ctx->nucleus.provided)) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_electron",
                           NULL);
  }

  if(!(ctx->electron.provided)) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_electron",
                           NULL);
  }

  if (!ctx->ao_basis.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_ao_basis",
                           NULL);
  }

  if (!ctx->mo_basis.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_mo_basis",
                           NULL);
  }

  if (!ctx->det.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_mo_basis",
                           NULL);
  }

  /* Compute if necessary */
  if (ctx->electron.walker.point.date > ctx->det.det_vgl_beta_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      free(ctx->det.det_vgl_beta);
      ctx->det.det_vgl_beta = NULL;
    }
    
    /* Allocate array */
    if (ctx->det.det_vgl_beta == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = 5 * ctx->electron.walker.num * ctx->det.det_num_beta *
                      ctx->electron.down_num * ctx->electron.down_num * sizeof(double);
      double* det_vgl_beta = (double*) qmckl_malloc(context, mem_info);

      if (det_vgl_beta == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_det_vgl_beta",
                               NULL);
      }
      ctx->det.det_vgl_beta = det_vgl_beta;
    }

    qmckl_exit_code rc;
    if (ctx->det.type == 'G') {
      rc = qmckl_compute_det_vgl_beta(context,
                                 ctx->det.det_num_beta,
                                 ctx->electron.walker.num,
                                 ctx->electron.up_num,
                                 ctx->electron.down_num,
                                 ctx->electron.num,
                                 ctx->det.mo_index_beta,
                                 ctx->mo_basis.mo_num,
                                 ctx->mo_basis.mo_vgl,
                                 ctx->det.det_vgl_beta);
    } else {
      return qmckl_failwith( context,
                             QMCKL_FAILURE,
                             "compute_det_vgl_beta",
                             "Not yet implemented");
    }
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->det.det_vgl_beta_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_get_det_inv_matrix_alpha(qmckl_context context, double * const det_inv_matrix_alpha) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_ao_basis_ao_vgl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_mo_basis_mo_vgl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_det_vgl_alpha(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_det_inv_matrix_alpha(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  size_t sze = ctx->det.det_num_alpha * ctx->electron.walker.num * ctx->electron.up_num * ctx->electron.up_num;
  memcpy(det_inv_matrix_alpha, ctx->det.det_inv_matrix_alpha, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_get_det_inv_matrix_beta(qmckl_context context, double * const det_inv_matrix_beta) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_ao_basis_ao_vgl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_mo_basis_mo_vgl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_det_vgl_beta(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_det_inv_matrix_beta(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  size_t sze = ctx->det.det_num_alpha * ctx->electron.walker.num * ctx->electron.down_num * ctx->electron.down_num;
  memcpy(det_inv_matrix_beta, ctx->det.det_inv_matrix_beta, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_get_det_adj_matrix_alpha(qmckl_context context, double * const det_adj_matrix_alpha) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_ao_basis_ao_vgl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_mo_basis_mo_vgl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_det_vgl_alpha(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_det_inv_matrix_alpha(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  size_t sze = ctx->det.det_num_alpha * ctx->electron.walker.num * ctx->electron.up_num * ctx->electron.up_num;
  memcpy(det_adj_matrix_alpha, ctx->det.det_adj_matrix_alpha, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_get_det_adj_matrix_beta(qmckl_context context, double * const det_adj_matrix_beta) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_ao_basis_ao_vgl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_mo_basis_mo_vgl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_det_vgl_beta(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_det_inv_matrix_beta(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  size_t sze = ctx->det.det_num_alpha * ctx->electron.walker.num * ctx->electron.down_num * ctx->electron.down_num;
  memcpy(det_adj_matrix_beta, ctx->det.det_adj_matrix_beta, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_get_det_alpha(qmckl_context context, double * const det_value_alpha) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_ao_basis_ao_vgl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_mo_basis_mo_vgl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_det_vgl_alpha(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_det_inv_matrix_alpha(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  size_t sze = ctx->det.det_num_alpha * ctx->electron.walker.num;
  memcpy(det_value_alpha, ctx->det.det_value_alpha, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_get_det_beta(qmckl_context context, double * const det_value_beta) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_ao_basis_ao_vgl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_mo_basis_mo_vgl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_det_vgl_beta(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_det_inv_matrix_beta(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  size_t sze = ctx->det.det_num_alpha * ctx->electron.walker.num;
  memcpy(det_value_beta, ctx->det.det_value_beta, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_det_inv_matrix_alpha(qmckl_context context) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if(!(ctx->nucleus.provided)) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_electron",
                           NULL);
  }

  if(!(ctx->electron.provided)) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_electron",
                           NULL);
  }

  if (!ctx->ao_basis.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_ao_basis",
                           NULL);
  }

  if (!ctx->mo_basis.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_mo_basis",
                           NULL);
  }

  if (!ctx->det.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_mo_basis",
                           NULL);
  }

  /* Compute if necessary */
  if (ctx->electron.walker.point.date > ctx->det.det_inv_matrix_alpha_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      free(ctx->det.det_inv_matrix_alpha);
      ctx->det.det_inv_matrix_alpha = NULL;
    }
    
    /* Allocate array */
    if (ctx->det.det_inv_matrix_alpha == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.walker.num * ctx->det.det_num_alpha *
                      ctx->electron.up_num * ctx->electron.up_num * sizeof(double);
      double* det_inv_matrix_alpha = (double*) qmckl_malloc(context, mem_info);

      if (det_inv_matrix_alpha == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_det_inv_matrix_alpha",
                               NULL);
      }
      ctx->det.det_inv_matrix_alpha = det_inv_matrix_alpha;
    }

    if (ctx->det.det_adj_matrix_alpha == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.walker.num * ctx->det.det_num_alpha *
                      ctx->electron.up_num * ctx->electron.up_num * sizeof(double);
      double* det_adj_matrix_alpha = (double*) qmckl_malloc(context, mem_info);

      if (det_adj_matrix_alpha == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_det_adj_matrix_alpha",
                               NULL);
      }
      ctx->det.det_adj_matrix_alpha = det_adj_matrix_alpha;
    }

    if (ctx->det.det_value_alpha == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.walker.num * ctx->det.det_num_alpha * sizeof(double);
      double* det_value_alpha = (double*) qmckl_malloc(context, mem_info);

      if (det_value_alpha == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_det_value_alpha",
                               NULL);
      }
      ctx->det.det_value_alpha = det_value_alpha;
    }

    qmckl_exit_code rc;
    if (ctx->det.type == 'G') {
      rc = qmckl_compute_det_inv_matrix_alpha(context,
                                        ctx->det.det_num_alpha,
                                        ctx->electron.walker.num,
                                        ctx->electron.up_num,
                                        ctx->det.det_vgl_alpha,
                                        ctx->det.det_value_alpha,
                                        ctx->det.det_adj_matrix_alpha,
                                        ctx->det.det_inv_matrix_alpha);
    } else {
      return qmckl_failwith( context,
                             QMCKL_FAILURE,
                             "compute_det_inv_matrix_alpha",
                             "Not yet implemented");
    }
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->det.det_value_alpha_date = ctx->date;
    ctx->det.det_adj_matrix_alpha_date = ctx->date;
    ctx->det.det_inv_matrix_alpha_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_det_inv_matrix_beta(qmckl_context context) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if(!(ctx->nucleus.provided)) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_electron",
                           NULL);
  }

  if(!(ctx->electron.provided)) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_electron",
                           NULL);
  }

  if (!ctx->ao_basis.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_ao_basis",
                           NULL);
  }

  if (!ctx->mo_basis.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_mo_basis",
                           NULL);
  }

  if (!ctx->det.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_mo_basis",
                           NULL);
  }

  /* Compute if necessary */
  if (ctx->electron.walker.point.date > ctx->det.det_inv_matrix_beta_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      free(ctx->det.det_inv_matrix_beta);
      ctx->det.det_inv_matrix_beta = NULL;
    }
    
    /* Allocate array */
    if (ctx->det.det_inv_matrix_beta == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.walker.num * ctx->det.det_num_beta *
                      ctx->electron.down_num * ctx->electron.down_num * sizeof(double);
      double* det_inv_matrix_beta = (double*) qmckl_malloc(context, mem_info);

      if (det_inv_matrix_beta == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_det_inv_matrix_beta",
                               NULL);
      }
      ctx->det.det_inv_matrix_beta = det_inv_matrix_beta;
    }

    if (ctx->det.det_adj_matrix_beta == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.walker.num * ctx->det.det_num_beta *
                      ctx->electron.down_num * ctx->electron.down_num * sizeof(double);
      double* det_adj_matrix_beta = (double*) qmckl_malloc(context, mem_info);

      if (det_adj_matrix_beta == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_det_adj_matrix_beta",
                               NULL);
      }
      ctx->det.det_adj_matrix_beta = det_adj_matrix_beta;
    }

    if (ctx->det.det_value_beta == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.walker.num * ctx->det.det_num_beta * sizeof(double);
      double* det_value_beta = (double*) qmckl_malloc(context, mem_info);

      if (det_value_beta == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_det_value_beta",
                               NULL);
      }
      ctx->det.det_value_beta = det_value_beta;
    }

    qmckl_exit_code rc;
    if (ctx->det.type == 'G') {
      rc = qmckl_compute_det_inv_matrix_beta(context,
                                        ctx->det.det_num_beta,
                                        ctx->electron.walker.num,
                                        ctx->electron.down_num,
                                        ctx->det.det_vgl_beta,
                                        ctx->det.det_value_beta,
                                        ctx->det.det_adj_matrix_beta,
                                        ctx->det.det_inv_matrix_beta);
    } else {
      return qmckl_failwith( context,
                             QMCKL_FAILURE,
                             "compute_det_inv_matrix_beta",
                             "Not yet implemented");
    }
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->det.det_value_beta_date = ctx->date;
    ctx->det.det_adj_matrix_beta_date = ctx->date;
    ctx->det.det_inv_matrix_beta_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}
