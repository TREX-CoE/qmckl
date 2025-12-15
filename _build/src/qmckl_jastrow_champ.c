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
#include <math.h>


#include <stdio.h>

#include "qmckl.h"
#include "qmckl_context_private_type.h"
#include "qmckl_memory_private_type.h"
#include "qmckl_memory_private_func.h"
#include "qmckl_jastrow_champ_private_type.h"
#include "qmckl_jastrow_champ_private_func.h"

qmckl_exit_code qmckl_init_jastrow_champ(qmckl_context context) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return false;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  ctx->jastrow_champ.uninitialized = (1 << 11) - 1;

  /* Default values */
  ctx->jastrow_champ.aord_num = -1;
  ctx->jastrow_champ.bord_num = -1;
  ctx->jastrow_champ.cord_num = -1;
  ctx->jastrow_champ.dim_c_vector = -1;
  ctx->jastrow_champ.type_nucl_num = -1;
  ctx->jastrow_champ.spin_independent = -1;

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_set_jastrow_champ_aord_num(qmckl_context context, const int64_t aord_num)
{

  int32_t mask = 1 << 0;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
   }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  
  if (mask != 0 && !(ctx->jastrow_champ.uninitialized & mask)) {
      return qmckl_failwith( context,
                             QMCKL_ALREADY_SET,
                             "qmckl_set_jastrow_champ_*",
                             NULL);
   }

      if (aord_num < 0) {
        return qmckl_failwith( context,
                               QMCKL_INVALID_ARG_2,
                               "qmckl_set_jastrow_champ_aord_num",
                               "aord_num < 0");
      }

  ctx->jastrow_champ.aord_num = aord_num;
  ctx->jastrow_champ.uninitialized |= (1 << 5);

  ctx->jastrow_champ.uninitialized &= ~mask;
  ctx->jastrow_champ.provided = (ctx->jastrow_champ.uninitialized == 0);
  if (ctx->jastrow_champ.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_jastrow_champ(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
   }
  
  return QMCKL_SUCCESS;
      }

qmckl_exit_code
qmckl_set_jastrow_champ_bord_num(qmckl_context context, const int64_t bord_num)
{

  int32_t mask = 1 << 1;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
   }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  
  if (mask != 0 && !(ctx->jastrow_champ.uninitialized & mask)) {
      return qmckl_failwith( context,
                             QMCKL_ALREADY_SET,
                             "qmckl_set_jastrow_champ_*",
                             NULL);
   }

      if (bord_num < 0) {
        return qmckl_failwith( context,
                               QMCKL_INVALID_ARG_2,
                               "qmckl_set_jastrow_champ_bord_num",
                               "bord_num < 0");
      }

  ctx->jastrow_champ.bord_num = bord_num;
  ctx->jastrow_champ.uninitialized |= (1 << 6);

  ctx->jastrow_champ.uninitialized &= ~mask;
  ctx->jastrow_champ.provided = (ctx->jastrow_champ.uninitialized == 0);
  if (ctx->jastrow_champ.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_jastrow_champ(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
   }
  
  return QMCKL_SUCCESS;
      }

qmckl_exit_code
qmckl_set_jastrow_champ_cord_num(qmckl_context context, const int64_t cord_num)
{

  int32_t mask = 1 << 2;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
   }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  
  if (mask != 0 && !(ctx->jastrow_champ.uninitialized & mask)) {
      return qmckl_failwith( context,
                             QMCKL_ALREADY_SET,
                             "qmckl_set_jastrow_champ_*",
                             NULL);
   }

      if (cord_num < 0) {
        return qmckl_failwith( context,
                               QMCKL_INVALID_ARG_2,
                               "qmckl_set_jastrow_champ_cord_num",
                               "cord_num < 0");
      }

  int64_t dim_c_vector = -1;
  qmckl_exit_code rc = qmckl_compute_dim_c_vector(context, cord_num, &dim_c_vector);
  assert (rc == QMCKL_SUCCESS);

  ctx->jastrow_champ.cord_num = cord_num;
  ctx->jastrow_champ.dim_c_vector = dim_c_vector;

  // If cord_num == 0, a_vector can't be set
  if (cord_num > 0) {
    ctx->jastrow_champ.uninitialized |= (1 << 7);
  } else {
    ctx->jastrow_champ.uninitialized &= ~(1 << 7);
  }

  ctx->jastrow_champ.uninitialized &= ~mask;
  ctx->jastrow_champ.provided = (ctx->jastrow_champ.uninitialized == 0);
  if (ctx->jastrow_champ.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_jastrow_champ(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
   }
  
  return QMCKL_SUCCESS;
           }


qmckl_exit_code
qmckl_set_jastrow_champ_type_nucl_num(qmckl_context context, const int64_t type_nucl_num)
{
  int32_t mask = 1 << 3;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
   }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  
  if (mask != 0 && !(ctx->jastrow_champ.uninitialized & mask)) {
      return qmckl_failwith( context,
                             QMCKL_ALREADY_SET,
                             "qmckl_set_jastrow_champ_*",
                             NULL);
   }

      if (type_nucl_num <= 0) {
        return qmckl_failwith( context,
                               QMCKL_INVALID_ARG_2,
                               "qmckl_set_jastrow_champ_type_nucl_num",
                               "type_nucl_num < 0");
      }

  ctx->jastrow_champ.type_nucl_num = type_nucl_num;

  ctx->jastrow_champ.uninitialized &= ~mask;
  ctx->jastrow_champ.provided = (ctx->jastrow_champ.uninitialized == 0);
  if (ctx->jastrow_champ.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_jastrow_champ(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
   }
  
  return QMCKL_SUCCESS;
      }


qmckl_exit_code
qmckl_set_jastrow_champ_type_nucl_vector(qmckl_context context,
                                         int64_t const * type_nucl_vector,
                                         const int64_t nucl_num)
{

  int32_t mask = 1 << 4;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
   }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  
  if (mask != 0 && !(ctx->jastrow_champ.uninitialized & mask)) {
      return qmckl_failwith( context,
                             QMCKL_ALREADY_SET,
                             "qmckl_set_jastrow_champ_*",
                             NULL);
   }

      int64_t type_nucl_num = ctx->jastrow_champ.type_nucl_num;

  if (type_nucl_num <= 0) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_set_jastrow_champ_type_nucl_vector",
                           "type_nucl_num not initialized");
  }

  if (type_nucl_vector == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_set_jastrow_champ_type_nucl_vector",
                           "type_nucl_vector = NULL");
  }

  for (int i=0 ; i<nucl_num ; ++i) {
    if (type_nucl_vector[i] < 0) {
      return qmckl_failwith( context, QMCKL_INVALID_ARG_2,
                             "qmckl_set_jastrow_champ_type_nucl_vector",
                             "Inconsistent values of type_nucl_vector (<0)" );
    }
    if (type_nucl_vector[i] >= type_nucl_num) {
      return qmckl_failwith( context, QMCKL_INVALID_ARG_2,
                             "qmckl_set_jastrow_champ_type_nucl_vector",
                             "Inconsistent values of type_nucl_vector (>=nucl_num). Values should use 0-based indexing as in C." );
    }
  }

  if (ctx->jastrow_champ.type_nucl_vector != NULL) {
    qmckl_exit_code rc = qmckl_free(context, ctx->jastrow_champ.type_nucl_vector);
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context, rc,
                             "qmckl_set_jastrow_champ_type_nucl_vector",
                             "Unable to free ctx->jastrow_champ.type_nucl_vector");
    }
  }

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = nucl_num * sizeof(int64_t);
  int64_t* new_array = (int64_t*) qmckl_malloc(context, mem_info);

  if(new_array == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_set_jastrow_champ_type_nucl_vector",
                           NULL);
  }

  memcpy(new_array, type_nucl_vector, mem_info.size);

  ctx->jastrow_champ.type_nucl_vector = new_array;

  ctx->jastrow_champ.uninitialized &= ~mask;
  ctx->jastrow_champ.provided = (ctx->jastrow_champ.uninitialized == 0);
  if (ctx->jastrow_champ.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_jastrow_champ(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
   }
  
  return QMCKL_SUCCESS;
      }


qmckl_exit_code
qmckl_set_jastrow_champ_a_vector(qmckl_context context,
                                 double const * a_vector,
                                 const int64_t size_max)
{
  int32_t mask = 1 << 5;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
   }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  
  if (mask != 0 && !(ctx->jastrow_champ.uninitialized & mask)) {
      return qmckl_failwith( context,
                             QMCKL_ALREADY_SET,
                             "qmckl_set_jastrow_champ_*",
                             NULL);
   }

      int64_t aord_num = ctx->jastrow_champ.aord_num;
  if (aord_num < 0) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_set_jastrow_champ_a_vector",
                           "aord_num not initialized");
  }

  int64_t type_nucl_num = ctx->jastrow_champ.type_nucl_num;

  if (type_nucl_num <= 0) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_set_jastrow_champ_a_vector",
                           "type_nucl_num not initialized");
  }

  if (a_vector == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_set_jastrow_champ_a_vector",
                           "a_vector = NULL");
  }

  if (ctx->jastrow_champ.a_vector != NULL) {
    qmckl_exit_code rc = qmckl_free(context, ctx->jastrow_champ.a_vector);
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context, rc,
                             "qmckl_set_jastrow_champ_a_vector",
                             "Unable to free ctx->jastrow_champ.a_vector");
    }
  }

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = (aord_num + 1) * type_nucl_num * sizeof(double);

  if (size_max < (aord_num+1)*type_nucl_num ) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_set_jastrow_champ_a_vector",
                           "Array too small. Expected (aord_num+1)*type_nucl_num");
  }

  double* new_array = (double*) qmckl_malloc(context, mem_info);

  if(new_array == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_set_jastrow_champ_a_vector",
                           NULL);
  }

  memcpy(new_array, a_vector, mem_info.size);

  ctx->jastrow_champ.a_vector = new_array;

  ctx->jastrow_champ.uninitialized &= ~mask;
  ctx->jastrow_champ.provided = (ctx->jastrow_champ.uninitialized == 0);
  if (ctx->jastrow_champ.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_jastrow_champ(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
   }
  
  return QMCKL_SUCCESS;
      }


qmckl_exit_code
qmckl_set_jastrow_champ_b_vector(qmckl_context context,
                                 double const * b_vector,
                                 const int64_t size_max)
{
  int32_t mask = 1 << 6;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
   }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  
  if (mask != 0 && !(ctx->jastrow_champ.uninitialized & mask)) {
      return qmckl_failwith( context,
                             QMCKL_ALREADY_SET,
                             "qmckl_set_jastrow_champ_*",
                             NULL);
   }

      int64_t bord_num = ctx->jastrow_champ.bord_num;
  if (bord_num < 0) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_set_jastrow_champ_b_vector",
                           "bord_num not initialized");
  }

  if (b_vector == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_set_jastrow_champ_b_vector",
                           "b_vector = NULL");
  }

  if (ctx->jastrow_champ.b_vector != NULL) {
    qmckl_exit_code rc = qmckl_free(context, ctx->jastrow_champ.b_vector);
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context, rc,
                             "qmckl_set_jastrow_champ_b_vector",
                             "Unable to free ctx->jastrow_champ.b_vector");
    }
  }

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = (bord_num + 1) * sizeof(double);

  if (size_max < (bord_num+1)) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_set_jastrow_champ_b_vector",
                           "Array too small. Expected (bord_num+1)");
  }

  double* new_array = (double*) qmckl_malloc(context, mem_info);

  if(new_array == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_set_jastrow_champ_b_vector",
                           NULL);
  }

  memcpy(new_array, b_vector, mem_info.size);

  ctx->jastrow_champ.b_vector = new_array;

  ctx->jastrow_champ.uninitialized &= ~mask;
  ctx->jastrow_champ.provided = (ctx->jastrow_champ.uninitialized == 0);
  if (ctx->jastrow_champ.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_jastrow_champ(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
   }
  
  return QMCKL_SUCCESS;
      }


qmckl_exit_code
qmckl_set_jastrow_champ_c_vector(qmckl_context context,
                                 double const * c_vector,
                                 const int64_t size_max)
{
  int32_t mask = 1 << 7;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
   }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  
  if (mask != 0 && !(ctx->jastrow_champ.uninitialized & mask)) {
      return qmckl_failwith( context,
                             QMCKL_ALREADY_SET,
                             "qmckl_set_jastrow_champ_*",
                             NULL);
   }

      int64_t type_nucl_num = ctx->jastrow_champ.type_nucl_num;
  if (type_nucl_num <= 0) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_set_jastrow_champ_c_vector",
                           "type_nucl_num not initialized");
  }

  int64_t dim_c_vector = ctx->jastrow_champ.dim_c_vector;
  if (dim_c_vector < 0) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_set_jastrow_champ_c_vector",
                           "cord_num not initialized");
  }

  if (c_vector == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_set_jastrow_champ_c_vector",
                           "c_vector = NULL");
  }

  if (ctx->jastrow_champ.c_vector != NULL) {
    qmckl_exit_code rc = qmckl_free(context, ctx->jastrow_champ.c_vector);
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context, rc,
                             "qmckl_set_jastrow_champ_c_vector",
                             "Unable to free ctx->jastrow_champ.c_vector");
    }
  }

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = dim_c_vector*type_nucl_num * sizeof(double);

  if (size_max < dim_c_vector*type_nucl_num) {
    char msg[256];
    sprintf(msg, "Array too small. Expected dim_c_vector*type_nucl_num = %ld", (long)
            (dim_c_vector*type_nucl_num) );
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_set_jastrow_champ_c_vector",
                           msg);
  }

  double* new_array = (double*) qmckl_malloc(context, mem_info);

  if(new_array == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_set_jastrow_champ_c_vector",
                           NULL);
  }

  memcpy(new_array, c_vector, mem_info.size);

  ctx->jastrow_champ.c_vector = new_array;

  ctx->jastrow_champ.uninitialized &= ~mask;
  ctx->jastrow_champ.provided = (ctx->jastrow_champ.uninitialized == 0);
  if (ctx->jastrow_champ.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_jastrow_champ(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
   }
  
  return QMCKL_SUCCESS;
      }

qmckl_exit_code
qmckl_set_jastrow_champ_rescale_factor_ee(qmckl_context context,
                                          const double rescale_factor_ee) {

  int32_t mask = 1 << 8;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
   }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  
  if (mask != 0 && !(ctx->jastrow_champ.uninitialized & mask)) {
      return qmckl_failwith( context,
                             QMCKL_ALREADY_SET,
                             "qmckl_set_jastrow_champ_*",
                             NULL);
   }

      if (rescale_factor_ee <= 0.0) {
        return qmckl_failwith( context,
                               QMCKL_INVALID_ARG_2,
                               "qmckl_set_jastrow_champ_rescale_factor_ee",
                               "rescale_factor_ee <= 0.0");
      }

  ctx->jastrow_champ.rescale_factor_ee = rescale_factor_ee;

  ctx->jastrow_champ.uninitialized &= ~mask;
  ctx->jastrow_champ.provided = (ctx->jastrow_champ.uninitialized == 0);
  if (ctx->jastrow_champ.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_jastrow_champ(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
   }
  
  return QMCKL_SUCCESS;
      }


qmckl_exit_code
qmckl_set_jastrow_champ_rescale_factor_en(qmckl_context context,
                                          const double* rescale_factor_en,
                                          const int64_t size_max) {

  int32_t mask = 1 << 9;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
   }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  
  if (mask != 0 && !(ctx->jastrow_champ.uninitialized & mask)) {
      return qmckl_failwith( context,
                             QMCKL_ALREADY_SET,
                             "qmckl_set_jastrow_champ_*",
                             NULL);
   }

      if (ctx->jastrow_champ.type_nucl_num <= 0) {
        return qmckl_failwith( context,
                               QMCKL_NOT_PROVIDED,
                               "qmckl_set_jastrow_champ_rescale_factor_en",
                               "type_nucl_num not set");
      }


  if (rescale_factor_en == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_set_jastrow_champ_rescale_factor_en",
                           "Null pointer");
  }

  if (size_max < ctx->jastrow_champ.type_nucl_num) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_set_jastrow_champ_rescale_factor_en",
                           "Array too small. Expected type_nucl_num.");
  }


  if (ctx->jastrow_champ.rescale_factor_en != NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_set_jastrow_champ_rescale_factor_en",
                           "Already set");
  }

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = ctx->jastrow_champ.type_nucl_num * sizeof(double);
  ctx->jastrow_champ.rescale_factor_en = (double*) qmckl_malloc(context, mem_info);

  for (int64_t i=0 ; i<ctx->jastrow_champ.type_nucl_num ; ++i) {
    if (rescale_factor_en[i] <= 0.0) {
      return qmckl_failwith( context,
                             QMCKL_INVALID_ARG_2,
                             "qmckl_set_jastrow_champ_rescale_factor_en",
                             "rescale_factor_en <= 0.0");
    }
    ctx->jastrow_champ.rescale_factor_en[i] = rescale_factor_en[i];
  }

  ctx->jastrow_champ.uninitialized &= ~mask;
  ctx->jastrow_champ.provided = (ctx->jastrow_champ.uninitialized == 0);
  if (ctx->jastrow_champ.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_jastrow_champ(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
   }
  
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_set_jastrow_champ_spin_independent(qmckl_context context, const int32_t spin_independent)
{
  int32_t mask = 1 << 10;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
   }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  
  if (mask != 0 && !(ctx->jastrow_champ.uninitialized & mask)) {
      return qmckl_failwith( context,
                             QMCKL_ALREADY_SET,
                             "qmckl_set_jastrow_champ_*",
                             NULL);
   }

  ctx->jastrow_champ.spin_independent = spin_independent;

  ctx->jastrow_champ.uninitialized &= ~mask;
  ctx->jastrow_champ.provided = (ctx->jastrow_champ.uninitialized == 0);
  if (ctx->jastrow_champ.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_jastrow_champ(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
   }
  
  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_finalize_jastrow_champ(qmckl_context context) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_finalize_jastrow_champ",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  /* ----------------------------------- */
  /* Check for the necessary information */
  /* ----------------------------------- */

  if (!(ctx->electron.provided)) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_electron",
                           NULL);
  }

  if (!(ctx->nucleus.provided)) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_nucleus",
                           NULL);
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_jastrow_champ_asymp_jasa(context);
  assert(rc == QMCKL_SUCCESS);

  rc = qmckl_provide_jastrow_champ_asymp_jasb(context);
  assert(rc == QMCKL_SUCCESS);

  rc = qmckl_provide_jastrow_champ_asymp_jasa_pderiv(context);
  assert(rc == QMCKL_SUCCESS);

  rc = qmckl_provide_jastrow_champ_asymp_jasb_pderiv(context);
  assert(rc == QMCKL_SUCCESS);

  rc = qmckl_context_touch(context);
  return rc;


}



/* Performs a deep copy of the Jastrow Champ structure from ~src~ to ~dest~. */
/* Memory allocations are done using the provided context. */


qmckl_exit_code qmckl_copy_jastrow_champ(qmckl_context context, const qmckl_jastrow_champ_struct* src, qmckl_jastrow_champ_struct* dest) {
  
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  }

  if (src == NULL || dest == NULL) {
    return qmckl_failwith(context,
                          QMCKL_INVALID_ARG_2,
                          "qmckl_copy_jastrow_champ",
                          "NULL pointer");
  }

  /* If source is not provided/initialized, just copy the struct */
  if (!src->provided) {
    *dest = *src;
    return QMCKL_SUCCESS;
  }

  /* Copy scalar fields */
  dest->aord_num = src->aord_num;
  dest->bord_num = src->bord_num;
  dest->cord_num = src->cord_num;
  dest->dim_c_vector = src->dim_c_vector;
  dest->type_nucl_num = src->type_nucl_num;
  dest->rescale_factor_ee = src->rescale_factor_ee;
  dest->uninitialized = src->uninitialized;
  dest->spin_independent = src->spin_independent;
  dest->provided = src->provided;

  /* Copy asymp_jasb array (fixed size of 2) */
  dest->asymp_jasb[0] = src->asymp_jasb[0];
  dest->asymp_jasb[1] = src->asymp_jasb[1];

  /* Copy all date fields */
  dest->asymp_jasa_date = src->asymp_jasa_date;
  dest->asymp_jasa_pderiv_date = src->asymp_jasa_pderiv_date;
  dest->asymp_jasb_date = src->asymp_jasb_date;
  dest->asymp_jasb_pderiv_date = src->asymp_jasb_pderiv_date;
  dest->c_vector_full_date = src->c_vector_full_date;
  dest->dim_c_vector_date = src->dim_c_vector_date;
  dest->dtmp_c_date = src->dtmp_c_date;
  dest->ee_distance_rescaled_date = src->ee_distance_rescaled_date;
  dest->ee_distance_rescaled_gl_date = src->ee_distance_rescaled_gl_date;
  dest->een_rescaled_e_date = src->een_rescaled_e_date;
  dest->een_rescaled_e_gl_date = src->een_rescaled_e_gl_date;
  dest->een_rescaled_n_date = src->een_rescaled_n_date;
  dest->een_rescaled_n_gl_date = src->een_rescaled_n_gl_date;
  dest->en_distance_rescaled_date = src->en_distance_rescaled_date;
  dest->en_distance_rescaled_gl_date = src->en_distance_rescaled_gl_date;
  dest->factor_ee_date = src->factor_ee_date;
  dest->factor_ee_gl_date = src->factor_ee_gl_date;
  dest->factor_ee_pderiv_date = src->factor_ee_pderiv_date;
  dest->factor_ee_gl_pderiv_date = src->factor_ee_gl_pderiv_date;
  dest->factor_een_date = src->factor_een_date;
  dest->factor_een_gl_date = src->factor_een_gl_date;
  dest->factor_een_grad_date = src->factor_een_grad_date;
  dest->factor_een_pderiv_date = src->factor_een_pderiv_date;
  dest->factor_een_gl_pderiv_date = src->factor_een_gl_pderiv_date;
  dest->factor_en_date = src->factor_en_date;
  dest->factor_en_gl_date = src->factor_en_gl_date;
  dest->factor_en_pderiv_date = src->factor_en_pderiv_date;
  dest->factor_en_gl_pderiv_date = src->factor_en_gl_pderiv_date;
  dest->lkpm_combined_index_date = src->lkpm_combined_index_date;
  dest->tmp_c_date = src->tmp_c_date;
  dest->value_date = src->value_date;
  dest->gl_date = src->gl_date;
  dest->grad_date = src->grad_date;

  /* Get nucl_num from context for array sizes */
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  int64_t nucl_num = ctx->nucleus.num;

  /* Deep copy type_nucl_vector array */
  if (src->type_nucl_vector != NULL && nucl_num > 0) {
    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = nucl_num * sizeof(int64_t);
    dest->type_nucl_vector = (int64_t*) qmckl_malloc(context, mem_info);
    if (dest->type_nucl_vector == NULL) return QMCKL_ALLOCATION_FAILED;
    memcpy(dest->type_nucl_vector, src->type_nucl_vector, mem_info.size);
  } else {
    dest->type_nucl_vector = NULL;
  }

  /* Deep copy a_vector array */
  if (src->a_vector != NULL && src->aord_num >= 0 && src->type_nucl_num > 0) {
    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = (src->aord_num + 1) * src->type_nucl_num * sizeof(double);
    dest->a_vector = (double*) qmckl_malloc(context, mem_info);
    if (dest->a_vector == NULL) return QMCKL_ALLOCATION_FAILED;
    memcpy(dest->a_vector, src->a_vector, mem_info.size);
  } else {
    dest->a_vector = NULL;
  }

  /* Deep copy b_vector array */
  if (src->b_vector != NULL && src->bord_num >= 0) {
    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = (src->bord_num + 1) * sizeof(double);
    dest->b_vector = (double*) qmckl_malloc(context, mem_info);
    if (dest->b_vector == NULL) return QMCKL_ALLOCATION_FAILED;
    memcpy(dest->b_vector, src->b_vector, mem_info.size);
  } else {
    dest->b_vector = NULL;
  }

  /* Deep copy c_vector array */
  if (src->c_vector != NULL && src->dim_c_vector >= 0 && src->type_nucl_num > 0) {
    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = src->dim_c_vector * src->type_nucl_num * sizeof(double);
    dest->c_vector = (double*) qmckl_malloc(context, mem_info);
    if (dest->c_vector == NULL) return QMCKL_ALLOCATION_FAILED;
    memcpy(dest->c_vector, src->c_vector, mem_info.size);
  } else {
    dest->c_vector = NULL;
  }

  /* Deep copy rescale_factor_en array */
  if (src->rescale_factor_en != NULL && src->type_nucl_num > 0) {
    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = src->type_nucl_num * sizeof(double);
    dest->rescale_factor_en = (double*) qmckl_malloc(context, mem_info);
    if (dest->rescale_factor_en == NULL) return QMCKL_ALLOCATION_FAILED;
    memcpy(dest->rescale_factor_en, src->rescale_factor_en, mem_info.size);
  } else {
    dest->rescale_factor_en = NULL;
  }

  /* All remaining arrays are computed/cached values.
     Set them to NULL and reset their dates to 0 so they will be recomputed when needed.
     This is more efficient than copying potentially large cached arrays. */
  
  dest->lkpm_combined_index = NULL;
  dest->lkpm_combined_index_date = 0;
  
  dest->asymp_jasa = NULL;
  dest->asymp_jasa_date = 0;
  
  dest->asymp_jasa_pderiv = NULL;
  dest->asymp_jasa_pderiv_date = 0;
  
  dest->asymp_jasb_pderiv = NULL;
  dest->asymp_jasb_pderiv_date = 0;
  
  dest->c_vector_full = NULL;
  dest->c_vector_full_date = 0;
  
  dest->dtmp_c = NULL;
  dest->dtmp_c_date = 0;
  
  dest->ee_distance_rescaled = NULL;
  dest->ee_distance_rescaled_date = 0;
  
  dest->ee_distance_rescaled_gl = NULL;
  dest->ee_distance_rescaled_gl_date = 0;
  
  dest->een_rescaled_e = NULL;
  dest->een_rescaled_e_date = 0;
  
  dest->een_rescaled_e_gl = NULL;
  dest->een_rescaled_e_gl_date = 0;
  
  dest->een_rescaled_n = NULL;
  dest->een_rescaled_n_date = 0;
  
  dest->een_rescaled_n_gl = NULL;
  dest->een_rescaled_n_gl_date = 0;
  
  dest->en_distance_rescaled = NULL;
  dest->en_distance_rescaled_date = 0;
  
  dest->en_distance_rescaled_gl = NULL;
  dest->en_distance_rescaled_gl_date = 0;
  
  dest->factor_ee = NULL;
  dest->factor_ee_date = 0;
  
  dest->factor_ee_gl = NULL;
  dest->factor_ee_gl_date = 0;
  
  dest->factor_ee_pderiv = NULL;
  dest->factor_ee_pderiv_date = 0;
  
  dest->factor_ee_gl_pderiv = NULL;
  dest->factor_ee_gl_pderiv_date = 0;
  
  dest->factor_een = NULL;
  dest->factor_een_date = 0;
  
  dest->factor_een_gl = NULL;
  dest->factor_een_gl_date = 0;
  
  dest->factor_een_grad = NULL;
  dest->factor_een_grad_date = 0;
  
  dest->factor_een_pderiv = NULL;
  dest->factor_een_pderiv_date = 0;
  
  dest->factor_een_gl_pderiv = NULL;
  dest->factor_een_gl_pderiv_date = 0;
  
  dest->factor_en = NULL;
  dest->factor_en_date = 0;
  
  dest->factor_en_gl = NULL;
  dest->factor_en_gl_date = 0;
  
  dest->factor_en_pderiv = NULL;
  dest->factor_en_pderiv_date = 0;
  
  dest->factor_en_gl_pderiv = NULL;
  dest->factor_en_gl_pderiv_date = 0;
  
  dest->tmp_c = NULL;
  dest->tmp_c_date = 0;
  
  dest->value = NULL;
  dest->value_date = 0;
  
  dest->gl = NULL;
  dest->gl_date = 0;
  
  dest->grad = NULL;
  dest->grad_date = 0;

  return QMCKL_SUCCESS;
}

bool qmckl_jastrow_champ_provided(const qmckl_context context) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return false;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  return ctx->jastrow_champ.provided;
}

qmckl_exit_code qmckl_get_jastrow_champ_aord_num (const qmckl_context context, int64_t* const aord_num) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_jastrow_champ_aord_num",
                           NULL);
  }

  if (aord_num == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_champ_aord_num",
                           "aord_num is a null pointer");
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 0;

  if ( (ctx->jastrow_champ.uninitialized & mask) != 0) {
    return QMCKL_NOT_PROVIDED;
  }

  assert (ctx->jastrow_champ.aord_num > 0);
  *aord_num = ctx->jastrow_champ.aord_num;
  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_get_jastrow_champ_bord_num (const qmckl_context context, int64_t* const bord_num) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_jastrow_champ_bord_num",
                           NULL);
  }

  if (bord_num == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_champ_bord_num",
                           "aord_num is a null pointer");
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 1;

  if ( (ctx->jastrow_champ.uninitialized & mask) != 0) {
    return QMCKL_NOT_PROVIDED;
  }

  assert (ctx->jastrow_champ.bord_num > 0);
  *bord_num = ctx->jastrow_champ.bord_num;
  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_get_jastrow_champ_cord_num (const qmckl_context context, int64_t* const cord_num) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_jastrow_champ_cord_num",
                           NULL);
  }

  if (cord_num == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_champ_cord_num",
                           "aord_num is a null pointer");
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 2;

  if ( (ctx->jastrow_champ.uninitialized & mask) != 0) {
    return QMCKL_NOT_PROVIDED;
  }

  assert (ctx->jastrow_champ.cord_num > 0);
  *cord_num = ctx->jastrow_champ.cord_num;
  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_get_jastrow_champ_type_nucl_num (const qmckl_context context, int64_t* const type_nucl_num) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_jastrow_champ_type_nucl_num",
                           NULL);
  }

  if (type_nucl_num == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_champ_type_nucl_num",
                           "type_nucl_num is a null pointer");
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 3;

  if ( (ctx->jastrow_champ.uninitialized & mask) != 0) {
    return QMCKL_NOT_PROVIDED;
  }

  assert (ctx->jastrow_champ.type_nucl_num > 0);
  *type_nucl_num = ctx->jastrow_champ.type_nucl_num;
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_jastrow_champ_type_nucl_vector (const qmckl_context context,
                                    int64_t* const type_nucl_vector,
                                    const int64_t size_max)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_jastrow_champ_type_nucl_vector",
                           NULL);
  }

  if (type_nucl_vector == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_champ_type_nucl_vector",
                           "type_nucl_vector is a null pointer");
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 4;

  if ( (ctx->jastrow_champ.uninitialized & mask) != 0) {
    return QMCKL_NOT_PROVIDED;
  }

  assert (ctx->jastrow_champ.type_nucl_vector != NULL);
  if (size_max < ctx->jastrow_champ.type_nucl_num) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_type_nucl_vector",
                           "Array too small. Expected jastrow_champ.type_nucl_num");
  }

  memcpy(type_nucl_vector, ctx->jastrow_champ.type_nucl_vector, ctx->jastrow_champ.type_nucl_num*sizeof(int64_t));
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_jastrow_champ_a_vector (const qmckl_context context,
                               double * const a_vector,
                               const int64_t size_max) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_jastrow_champ_a_vector",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 5;

  if ( (ctx->jastrow_champ.uninitialized & mask) != 0) {
    return QMCKL_NOT_PROVIDED;
  }

  if (a_vector == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_champ_a_vector",
                           "Null pointer");
  }

  assert (ctx->jastrow_champ.a_vector != NULL);

  const int64_t sze = (ctx->jastrow_champ.aord_num + 1)*ctx->jastrow_champ.type_nucl_num;

  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_a_vector",
                           "Array too small. Expected (aord_num + 1)*type_nucl_num");
  }

  memcpy(a_vector, ctx->jastrow_champ.a_vector, sze*sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_jastrow_champ_b_vector (const qmckl_context context,
                               double * const b_vector,
                               const int64_t size_max) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_jastrow_champ_b_vector",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 6;

  if ( (ctx->jastrow_champ.uninitialized & mask) != 0) {
    return QMCKL_NOT_PROVIDED;
  }

  assert (ctx->jastrow_champ.b_vector != NULL);

  if (b_vector == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_champ_b_vector",
                           "Null pointer");
  }

  const int64_t sze=ctx->jastrow_champ.bord_num +1;

  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_b_vector",
                           "Array too small. Expected bord_num + 1");
  }

  memcpy(b_vector, ctx->jastrow_champ.b_vector, sze*sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_jastrow_champ_c_vector (const qmckl_context context,
                               double * const c_vector,
                               const int64_t size_max) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_jastrow_champ_c_vector",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 7;

  if ( (ctx->jastrow_champ.uninitialized & mask) != 0) {
    return QMCKL_NOT_PROVIDED;
  }

  assert (ctx->jastrow_champ.c_vector != NULL);

  int64_t dim_c_vector;
  qmckl_exit_code rc = qmckl_get_jastrow_champ_dim_c_vector(context, &dim_c_vector);
  if (rc != QMCKL_SUCCESS) return rc;

  if (c_vector == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_champ_c_vector",
                           "c_vector is a null pointer");
  }

  const int64_t sze=dim_c_vector * ctx->jastrow_champ.type_nucl_num;

  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_c_vector",
                           "Array too small. Expected dim_c_vector*type_nucl_num");
  }

  memcpy(c_vector, ctx->jastrow_champ.c_vector, sze*sizeof(double));
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_jastrow_champ_rescale_factor_ee (const qmckl_context context,
                                           double* const rescale_factor_ee) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_jastrow_champ_rescale_factor_ee",
                           NULL);
  }

  if (rescale_factor_ee == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_champ_rescale_factor_ee",
                           "rescale_factor_ee is a null pointer");
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 8;

  if ( (ctx->jastrow_champ.uninitialized & mask) != 0) {
    return QMCKL_NOT_PROVIDED;
  }
  assert (ctx->jastrow_champ.rescale_factor_ee > 0.0);

  *rescale_factor_ee = ctx->jastrow_champ.rescale_factor_ee;
  return QMCKL_SUCCESS;
}


qmckl_exit_code
qmckl_get_jastrow_champ_rescale_factor_en (const qmckl_context context,
                                      double* const rescale_factor_en,
                                      const int64_t size_max) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_jastrow_champ_rescale_factor_en",
                           NULL);
  }

  if (rescale_factor_en == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_champ_rescale_factor_en",
                           "rescale_factor_en is a null pointer");
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 9;

  if ( (ctx->jastrow_champ.uninitialized & mask) != 0) {
    return QMCKL_NOT_PROVIDED;
  }

  if (size_max < ctx->jastrow_champ.type_nucl_num) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_rescale_factor_en",
                           "Array to small. Expected type_nucl_num.");
  }

  assert(ctx->jastrow_champ.rescale_factor_en != NULL);
  for (int64_t i=0 ; i<ctx->jastrow_champ.type_nucl_num ; ++i) {
    rescale_factor_en[i] = ctx->jastrow_champ.rescale_factor_en[i];
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_get_jastrow_champ_dim_c_vector(qmckl_context context, int64_t* const dim_c_vector)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_jastrow_champ_dim_c_vector",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  *dim_c_vector = ctx->jastrow_champ.dim_c_vector;

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_get_jastrow_champ_spin_independent(const qmckl_context context, int32_t* const spin_independent) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_jastrow_champ_spin_independent",
                           NULL);
  }

  if (spin_independent == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_champ_spin_independent",
                           "spin_independent is a null pointer");
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 10;

  if ( (ctx->jastrow_champ.uninitialized & mask) != 0) {
    return QMCKL_NOT_PROVIDED;
  }

  *spin_independent = ctx->jastrow_champ.spin_independent ;
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_jastrow_champ_asymp_jasb(qmckl_context context,
                             double* const asymp_jasb,
                             const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_jastrow_champ_asymp_jasb",
                           NULL);
  }


  /* Provided in finalize_jastrow */
  /*
  qmckl_exit_code rc;
  rc = qmckl_provide_jastrow_champ_asymp_jasb(context);
  if(rc != QMCKL_SUCCESS) return rc;
  */

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  const int64_t sze = 2;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_asymp_jasb",
                           "Array too small. Expected 2");
  }
  memcpy(asymp_jasb, ctx->jastrow_champ.asymp_jasb, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_jastrow_champ_asymp_jasb(qmckl_context context)
{

  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_provide_jastrow_champ_asymp_jasb",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->jastrow_champ.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_provide_jastrow_champ_asymp_jasb",
                           NULL);
  }

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow_champ.asymp_jasb_date) {

    rc = qmckl_compute_jastrow_champ_asymp_jasb(context,
                                  ctx->jastrow_champ.bord_num,
                                  ctx->jastrow_champ.b_vector,
                                  ctx->jastrow_champ.rescale_factor_ee,
                                  ctx->jastrow_champ.spin_independent,
                                  ctx->jastrow_champ.asymp_jasb);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow_champ.asymp_jasb_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_compute_jastrow_champ_asymp_jasb_hpc (const qmckl_context context,
                                            const int64_t bord_num,
                                            const double* b_vector,
                                            const double rescale_factor_ee,
                                            const int32_t spin_independent,
                                            double* const asymp_jasb )
{

  if (context == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  }

  if (bord_num < 0) {
    return QMCKL_INVALID_ARG_2;
  }

  const double kappa_inv = 1.0 / rescale_factor_ee;
  const double asym_one = b_vector[0] * kappa_inv / (1.0 + b_vector[1] * kappa_inv);

  double f = 0.;
  double x = kappa_inv;
  for (int k = 2; k <= bord_num; ++k) {
    x *= kappa_inv;
    f = f + b_vector[k]*x;
  }

  asymp_jasb[0] = spin_independent == 1 ? asym_one + f : 0.5 * asym_one + f;
  asymp_jasb[1] = asym_one + f;

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_compute_jastrow_champ_asymp_jasb (const qmckl_context context,
                                                        const int64_t bord_num,
                                                        const double* b_vector,
                                                        const double rescale_factor_ee,
                                                        const int32_t spin_independent,
                                                        double* const asymp_jasb )
{
#ifdef HAVE_HPC
  return qmckl_compute_jastrow_champ_asymp_jasb_hpc
#else
  return qmckl_compute_jastrow_champ_asymp_jasb_doc
#endif
    (context, bord_num, b_vector, rescale_factor_ee, spin_independent, asymp_jasb);
}

qmckl_exit_code qmckl_get_jastrow_champ_ee_distance_rescaled(qmckl_context context,
                                                             double* const distance_rescaled,
                                                             int64_t const size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_ee_distance_rescaled(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (distance_rescaled == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_champ_ee_distance_rescaled",
                           "Null pointer");
  }

  const int64_t sze = ctx->electron.num * ctx->electron.num * ctx->electron.walker.num;

  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_champ_ee_distance_rescaled",
                           "Array too small. Expected elec_num*elec_num*walk_num.");
  }
  memcpy(distance_rescaled, ctx->jastrow_champ.ee_distance_rescaled, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_ee_distance_rescaled(qmckl_context context)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);


  /* Compute if necessary */
  if (ctx->electron.walker.point.date > ctx->jastrow_champ.ee_distance_rescaled_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      if (ctx->jastrow_champ.ee_distance_rescaled != NULL) {
        qmckl_exit_code rc = qmckl_free(context, ctx->jastrow_champ.ee_distance_rescaled);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_ee_distance_rescaled",
                                 "Unable to free ctx->jastrow_champ.ee_distance_rescaled");
        }
        ctx->jastrow_champ.ee_distance_rescaled = NULL;
      }
    }

    /* Allocate array */
    if (ctx->jastrow_champ.ee_distance_rescaled == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.num * ctx->electron.num *
        ctx->electron.walker.num * sizeof(double);
      double* ee_distance_rescaled = (double*) qmckl_malloc(context, mem_info);

      if (ee_distance_rescaled == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_ee_distance_rescaled",
                               NULL);
      }
      ctx->jastrow_champ.ee_distance_rescaled = ee_distance_rescaled;
    }

    qmckl_exit_code rc =
      qmckl_compute_ee_distance_rescaled(context,
                                ctx->electron.num,
                                ctx->jastrow_champ.rescale_factor_ee,
                                ctx->electron.walker.num,
                                ctx->electron.walker.point.coord.data,
                                ctx->jastrow_champ.ee_distance_rescaled);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow_champ.ee_distance_rescaled_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_compute_ee_distance_rescaled_hpc (
          const qmckl_context context,
          const int64_t elec_num,
          const double rescale_factor_ee,
          const int64_t walk_num,
          const double* coord,
          double* const ee_distance_rescaled )
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  }

  if (elec_num <= 0) {
    return QMCKL_INVALID_ARG_2;
  }

  if (walk_num <= 0) {
    return QMCKL_INVALID_ARG_4;
  }

  if (coord == NULL) {
    return QMCKL_INVALID_ARG_5;
  }

  if (ee_distance_rescaled == NULL) {
    return QMCKL_INVALID_ARG_6;
  }


  const int64_t sze = elec_num*walk_num;
  const int64_t elec_num2= elec_num*elec_num;


  qmckl_exit_code result = QMCKL_SUCCESS;
#ifdef HAVE_OPENMP
  #pragma omp parallel
  {
#endif
    qmckl_exit_code rc = QMCKL_SUCCESS;
#ifdef HAVE_OPENMP
    #pragma omp for
#endif
    for (int64_t k=0 ; k<walk_num ; ++k)
    {
      rc |= qmckl_distance_rescaled(context, 'T', 'T', elec_num, elec_num,
                                    &(coord[k*elec_num]), sze, &(coord[k*elec_num]), sze,
                                    &(ee_distance_rescaled[k*elec_num2]), elec_num, rescale_factor_ee);
    }
#ifdef HAVE_OPENMP
    #pragma omp critical
#endif
    result |= rc;
#ifdef HAVE_OPENMP
  }
#endif
  return result;
}

qmckl_exit_code qmckl_compute_ee_distance_rescaled (
          const qmckl_context context,
          const int64_t elec_num,
          const double rescale_factor_ee,
          const int64_t walk_num,
          const double* coord,
          double* const ee_distance_rescaled )
{
#ifdef HAVE_HPC
  return qmckl_compute_ee_distance_rescaled_hpc
#else
  return qmckl_compute_ee_distance_rescaled_doc
#endif
 (context, elec_num, rescale_factor_ee, walk_num, coord, ee_distance_rescaled);
}

qmckl_exit_code qmckl_get_jastrow_champ_ee_distance_rescaled_gl(qmckl_context context,
                                                                double* const distance_rescaled_gl,
                                                                const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_ee_distance_rescaled_gl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (distance_rescaled_gl == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_champ_ee_distance_rescaled_gl",
                           "Null pointer.");
  }

  const int64_t sze = 4 * ctx->electron.num * ctx->electron.num * ctx->electron.walker.num;

  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_champ_ee_distance_rescaled_gl",
                           "Array too small. Expected 4*elec_num*elec_num*walk_num");
  }

  memcpy(distance_rescaled_gl, ctx->jastrow_champ.ee_distance_rescaled_gl, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_ee_distance_rescaled_gl(qmckl_context context)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);


  /* Compute if necessary */
  if (ctx->electron.walker.point.date > ctx->jastrow_champ.ee_distance_rescaled_gl_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      if (ctx->jastrow_champ.ee_distance_rescaled_gl != NULL) {
        qmckl_exit_code rc = qmckl_free(context, ctx->jastrow_champ.ee_distance_rescaled_gl);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_ee_distance_rescaled_gl",
                                 "Unable to free ctx->jastrow_champ.ee_distance_rescaled_gl");
        }
        ctx->jastrow_champ.ee_distance_rescaled_gl = NULL;
      }
    }

    /* Allocate array */
    if (ctx->jastrow_champ.ee_distance_rescaled_gl == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = 4 * ctx->electron.num * ctx->electron.num *
        ctx->electron.walker.num * sizeof(double);
      double* ee_distance_rescaled_gl = (double*) qmckl_malloc(context, mem_info);

      if (ee_distance_rescaled_gl == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_ee_distance_rescaled_gl",
                               NULL);
      }
      ctx->jastrow_champ.ee_distance_rescaled_gl = ee_distance_rescaled_gl;
    }

    qmckl_exit_code rc =
      qmckl_compute_ee_distance_rescaled_gl(context,
                                ctx->electron.num,
                                ctx->jastrow_champ.rescale_factor_ee,
                                ctx->electron.walker.num,
                                ctx->electron.walker.point.coord.data,
                                ctx->jastrow_champ.ee_distance_rescaled_gl);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow_champ.ee_distance_rescaled_gl_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_compute_ee_distance_rescaled_gl_hpc (
          const qmckl_context context,
          const int64_t elec_num,
          const double rescale_factor_ee,
          const int64_t walk_num,
          const double* coord,
          double* const ee_distance_rescaled_gl )
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) return QMCKL_INVALID_CONTEXT;
  if (elec_num <= 0) return QMCKL_INVALID_ARG_2;
  if (walk_num <= 0) return QMCKL_INVALID_ARG_4;
  if (coord == NULL) return QMCKL_INVALID_ARG_5;
  if (ee_distance_rescaled_gl == NULL) return QMCKL_INVALID_ARG_6;

  const int64_t sze = elec_num*walk_num;
  const int64_t elec_num2= elec_num*elec_num*4;

  qmckl_exit_code result = QMCKL_SUCCESS;
#ifdef HAVE_OPENMP
  #pragma omp parallel
#endif
  {
    qmckl_exit_code rc = QMCKL_SUCCESS;
#ifdef HAVE_OPENMP
    #pragma omp for
#endif
    for (int64_t k=0 ; k<walk_num ; ++k)
    {
      rc |= qmckl_distance_rescaled_gl(context, 'T', 'T', elec_num, elec_num,
                                       &(coord[k*elec_num]), sze,
                                       &(coord[k*elec_num]), sze,
                                       &(ee_distance_rescaled_gl[k*elec_num2]), elec_num,
                                       rescale_factor_ee);
    }
#ifdef HAVE_OPENMP
    #pragma omp critical
#endif
    result |= rc;
  }
  return result;
}

qmckl_exit_code qmckl_compute_ee_distance_rescaled_gl (
          const qmckl_context context,
          const int64_t elec_num,
          const double rescale_factor_ee,
          const int64_t walk_num,
          const double* coord,
          double* const ee_distance_rescaled_gl )
{
#ifdef HAVE_HPC
  return qmckl_compute_ee_distance_rescaled_gl_hpc
#else
  return qmckl_compute_ee_distance_rescaled_gl_doc
#endif
   (context, elec_num, rescale_factor_ee, walk_num, coord,
   ee_distance_rescaled_gl);
}

qmckl_exit_code
qmckl_get_jastrow_champ_factor_ee(qmckl_context context,
                            double* const factor_ee,
                            const int64_t size_max)
{
  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_jastrow_champ_factor_ee",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  rc = qmckl_provide_jastrow_champ_factor_ee(context);
  if (rc != QMCKL_SUCCESS) return rc;

  if (factor_ee == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_champ_factor_ee",
                           "Null pointer");
  }

  const int64_t sze = ctx->electron.walker.num;

  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_factor_ee",
                           "Array too small. Expected walk_num");
  }
  memcpy(factor_ee, ctx->jastrow_champ.factor_ee, sze*sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_jastrow_champ_factor_ee(qmckl_context context)
{

  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_provide_jastrow_champ_factor_ee",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->jastrow_champ.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_provide_jastrow_champ_factor_ee",
                           NULL);
  }

  rc = qmckl_provide_ee_distance_rescaled(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Provided in finalize_jastrow */
  /*
  rc = qmckl_provide_jastrow_champ_asymp_jasb(context);
  if(rc != QMCKL_SUCCESS) return rc;
  */

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow_champ.factor_ee_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      if (ctx->jastrow_champ.factor_ee != NULL) {
        rc = qmckl_free(context, ctx->jastrow_champ.factor_ee);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_jastrow_champ_factor_ee",
                                 "Unable to free ctx->jastrow_champ.factor_ee");
        }
        ctx->jastrow_champ.factor_ee = NULL;
      }
    }

    /* Allocate array */
    if (ctx->jastrow_champ.factor_ee == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.walker.num * sizeof(double);
      double* factor_ee = (double*) qmckl_malloc(context, mem_info);

      if (factor_ee == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_jastrow_champ_factor_ee",
                               NULL);
      }
      ctx->jastrow_champ.factor_ee = factor_ee;
    }

    rc = qmckl_compute_jastrow_champ_factor_ee(context,
                                 ctx->electron.walker.num,
                                 ctx->electron.num,
                                 ctx->electron.up_num,
                                 ctx->jastrow_champ.bord_num,
                                 ctx->jastrow_champ.b_vector,
                                 ctx->jastrow_champ.ee_distance_rescaled,
                                 ctx->jastrow_champ.asymp_jasb,
                                 ctx->jastrow_champ.spin_independent,
                                 ctx->jastrow_champ.factor_ee);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow_champ.factor_ee_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_ee_hpc (const qmckl_context context,
                                           const int64_t walk_num,
                                           const int64_t elec_num,
                                           const int64_t up_num,
                                           const int64_t bord_num,
                                           const double* restrict b_vector,
                                           const double* restrict ee_distance_rescaled,
                                           const double* restrict asymp_jasb,
                                           const int32_t spin_independent,
                                           double* restrict const factor_ee )
{

  if (context == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  }

  if (walk_num <= 0) {
    return QMCKL_INVALID_ARG_2;
  }

  if (elec_num <= 0) {
    return QMCKL_INVALID_ARG_3;
  }

  if (bord_num < 0) {
    return QMCKL_INVALID_ARG_4;
  }

  const int64_t dn_num = elec_num - up_num;
  const double fshift = 0.5 * (double) ((dn_num-1)*dn_num + (up_num-1)*up_num) * asymp_jasb[0] +
    (float) (up_num*dn_num) * asymp_jasb[1];

#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
  for (int nw = 0; nw < walk_num; ++nw) {
    double result = 0.;

    size_t ishift = nw * elec_num * elec_num;

    if (spin_independent == 1) {

      for (int j = 0; j < elec_num; ++j ) {
        const double* restrict xj = &(ee_distance_rescaled[j * elec_num + ishift]);
        for (int i = 0; i < j ; ++i) {
          result = result + b_vector[0]*xj[i] / (1. + b_vector[1]*xj[i]);
        }
      }

    } else {

      for (int j = 0; j < up_num; ++j ) {
        const double* restrict xj = &(ee_distance_rescaled[j * elec_num + ishift]);
        for (int i = 0; i < j ; ++i) {
          result = result + 0.5 * b_vector[0]*xj[i] / (1. + b_vector[1]*xj[i]);
        }
      }

      for (int j = up_num ; j < elec_num; ++j ) {
        const double* restrict xj = &(ee_distance_rescaled[j * elec_num + ishift]);
        for (int i = 0; i < up_num; ++i) {
          result = result + b_vector[0]*xj[i] / (1. + b_vector[1]*xj[i]);
        }
        for (int i = up_num ; i < j ; ++i) {
          result = result + 0.5 * b_vector[0]*xj[i] / (1. + b_vector[1]*xj[i]);
        }

      }
    }

    result = result - fshift;

    switch (bord_num) {
    case 2:
      {
        const double b2 = b_vector[2];
        for (int j=0; j < elec_num; ++j ) {
          const double* restrict xj = &(ee_distance_rescaled[j * elec_num + ishift]);
          for (int i=0; i < j ; ++i) {
            const double x = xj[i];
            double x2 = x*x;
            result = result + b2 * x2;
          }
        }
        break;
      }
    case 3:
      {
        const double b2 = b_vector[2];
        const double b3 = b_vector[3];
        for (int j=0; j < elec_num; ++j ) {
          const double* restrict xj = &(ee_distance_rescaled[j * elec_num + ishift]);
          for (int i=0; i < j ; ++i) {
            const double x = xj[i];
            double x2 = x*x;
            double x3 = x2*x;
            result = result + b2 * x2 + b3 * x3;
          
          }
        }
        break;
      }
    case 4:
      {
        const double b2 = b_vector[2];
        const double b3 = b_vector[3];
        const double b4 = b_vector[4];
        for (int j=0; j < elec_num; ++j ) {
          const double* restrict xj = &(ee_distance_rescaled[j * elec_num + ishift]);
          for (int i=0; i < j ; ++i) {
            const double x = xj[i];
            double x2 = x*x;
            double x3 = x2*x;
            double x4 = x2*x2;
            result = result + b2*x2 + b3*x3 + b4*x4;
          }
        }
        break;
      }
    case 5:
      {
        const double b2 = b_vector[2];
        const double b3 = b_vector[3];
        const double b4 = b_vector[4];
        const double b5 = b_vector[5];
        for (int j=0; j < elec_num; ++j ) {
          const double* restrict xj = &(ee_distance_rescaled[j * elec_num + ishift]);
          for (int i=0; i < j ; ++i) {
            const double x = xj[i];
            double x2 = x*x;
            double x3 = x2*x;
            double x4 = x2*x2;
            double x5 = x3*x2;
            result = result + b2*x2 + b3*x3 + b4*x4 + b5*x5;
          }
        }
        break;
      }
    default:
      {
        for (int j=0; j < elec_num; ++j ) {
          const double* restrict xj = &(ee_distance_rescaled[j * elec_num + ishift]);
          for (int i=0; i < j ; ++i) {
            const double x = xj[i];
            double xk = x;
            for (int k = 2; k <= bord_num; ++k) {
              xk *= x;
              result = result + b_vector[k] * xk;
            }
          
          }
        }
        break;
      }
    }
    factor_ee[nw] = result;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_ee (const qmckl_context context,
                                       const int64_t walk_num,
                                       const int64_t elec_num,
                                       const int64_t up_num,
                                       const int64_t bord_num,
                                       const double* b_vector,
                                       const double* ee_distance_rescaled,
                                       const double* asymp_jasb,
                                       const int32_t spin_independent,
                                       double* const factor_ee )
{

#ifdef HAVE_HPC
  return qmckl_compute_jastrow_champ_factor_ee_hpc
#else
  return qmckl_compute_jastrow_champ_factor_ee_doc
#endif
    (context, walk_num, elec_num, up_num, bord_num, b_vector,
     ee_distance_rescaled, asymp_jasb, spin_independent, factor_ee);
}

qmckl_exit_code
qmckl_get_jastrow_champ_factor_ee_gl(qmckl_context context,
                                    double* const factor_ee_gl,
                                    const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_jastrow_champ_factor_ee_gl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (factor_ee_gl == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_champ_factor_ee_gl",
                           "Null pointer");
  }

  const int64_t sze = ctx->electron.walker.num * 4 * ctx->electron.num;

  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_factor_ee_gl",
                           "Array too small. Expected 4*walk_num*elec_num");
  }

  memcpy(factor_ee_gl, ctx->jastrow_champ.factor_ee_gl, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_jastrow_champ_factor_ee_gl(qmckl_context context)
{

  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_provide_jastrow_champ_factor_ee_gl",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->jastrow_champ.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_provide_jastrow_champ_factor_ee_gl",
                           NULL);
  }

  /* Check if ee rescaled distance is provided */
  rc = qmckl_provide_ee_distance_rescaled(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Check if ee rescaled distance deriv e is provided */
  rc = qmckl_provide_ee_distance_rescaled_gl(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow_champ.factor_ee_gl_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      if (ctx->jastrow_champ.factor_ee_gl != NULL) {
        rc = qmckl_free(context, ctx->jastrow_champ.factor_ee_gl);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_jastrow_champ_factor_ee_gl",
                                 "Unable to free ctx->jastrow_champ.factor_ee_gl");
        }
        ctx->jastrow_champ.factor_ee_gl = NULL;
      }
    }

    /* Allocate array */
    if (ctx->jastrow_champ.factor_ee_gl == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.walker.num * 4 * ctx->electron.num * sizeof(double);
      double* factor_ee_gl = (double*) qmckl_malloc(context, mem_info);

      if (factor_ee_gl == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_jastrow_champ_factor_ee_gl",
                               NULL);
      }
      ctx->jastrow_champ.factor_ee_gl = factor_ee_gl;
    }

    rc = qmckl_compute_jastrow_champ_factor_ee_gl(context,
                                                  ctx->electron.walker.num,
                                                  ctx->electron.num,
                                                  ctx->electron.up_num,
                                                  ctx->jastrow_champ.bord_num,
                                                  ctx->jastrow_champ.b_vector,
                                                  ctx->jastrow_champ.ee_distance_rescaled,
                                                  ctx->jastrow_champ.ee_distance_rescaled_gl,
                                                  ctx->jastrow_champ.spin_independent,
                                                  ctx->jastrow_champ.factor_ee_gl);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow_champ.factor_ee_gl_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_ee_gl_hpc(const qmckl_context context,
                                             const int64_t walk_num,
                                             const int64_t elec_num,
                                             const int64_t up_num,
                                             const int64_t bord_num,
                                             const double* b_vector,
                                             const double* ee_distance_rescaled,
                                             const double* ee_distance_rescaled_gl,
                                             const int32_t spin_independent,
                                             double* const factor_ee_gl )
{

  if (context == QMCKL_NULL_CONTEXT) return QMCKL_INVALID_CONTEXT;
  if (walk_num <= 0) return QMCKL_INVALID_ARG_2;
  if (elec_num <= 0) return QMCKL_INVALID_ARG_3;
  if (up_num <= 0) return QMCKL_INVALID_ARG_4;
  if (bord_num < 0) return QMCKL_INVALID_ARG_5;
  if (b_vector == NULL) return QMCKL_INVALID_ARG_6;
  if (ee_distance_rescaled == NULL) return QMCKL_INVALID_ARG_7;
  if (ee_distance_rescaled_gl == NULL) return QMCKL_INVALID_ARG_8;
  if (spin_independent & (int32_t) (-2)) return QMCKL_INVALID_ARG_8;
  if (factor_ee_gl == NULL) return QMCKL_INVALID_ARG_9;

  double kf[bord_num+1];
  for (int k=0 ; k<=bord_num ; ++k) {
    kf[k] = (double) k;
  }

#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
  for (int nw = 0; nw < walk_num; ++nw) {
    double xk[bord_num+1];
    bool touched = false;

    for (int j = 0; j < elec_num; ++j) {
      const double* dxj = &ee_distance_rescaled_gl[4*elec_num*(j+nw*elec_num)];
      const double*  xj = &ee_distance_rescaled   [  elec_num*(j+nw*elec_num)];

      double * restrict factor_ee_gl_0 = &(factor_ee_gl[nw*elec_num*4]);
      double * restrict factor_ee_gl_1 = factor_ee_gl_0 + elec_num;
      double * restrict factor_ee_gl_2 = factor_ee_gl_1 + elec_num;
      double * restrict factor_ee_gl_3 = factor_ee_gl_2 + elec_num;

      for (int i = 0; i < elec_num; ++i) {
        if (j == i) {
          if (j == 0) {
            factor_ee_gl_0[i] = 0.0;
            factor_ee_gl_1[i] = 0.0;
            factor_ee_gl_2[i] = 0.0;
            factor_ee_gl_3[i] = 0.0;
          }
          continue;
        }

        double x = xj[i];

        const double denom      = 1.0 + b_vector[1]*x;
        const double invdenom   = 1.0  / denom;
        const double invdenom2  = invdenom * invdenom;

        const double* restrict dx = dxj + 4*i;

        const double grad_c2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];

        double f = b_vector[0] * invdenom2;
        if ((spin_independent == 0) && (
            ((i <  up_num) && (j <  up_num)) ||
            ((i >= up_num) && (j >= up_num))) ) {
          f *= 0.5;
        }

        if (touched) {
          factor_ee_gl_0[i] = factor_ee_gl_0[i] + f*dx[0];
          factor_ee_gl_1[i] = factor_ee_gl_1[i] + f*dx[1];
          factor_ee_gl_2[i] = factor_ee_gl_2[i] + f*dx[2];
          factor_ee_gl_3[i] = factor_ee_gl_3[i] + f*dx[3];
        } else {
          factor_ee_gl_0[i] = f*dx[0];
          factor_ee_gl_1[i] = f*dx[1];
          factor_ee_gl_2[i] = f*dx[2];
          factor_ee_gl_3[i] = f*dx[3];
        }
        factor_ee_gl_3[i] = factor_ee_gl_3[i] - f*grad_c2*invdenom*2.0 * b_vector[1];

        xk[0] = 1.0;
        for (int k=1 ; k<= bord_num ; ++k) {
          xk[k] = xk[k-1]*x;
        }

        for (int k=2 ; k<= bord_num ; ++k) {
          const double f1 = b_vector[k] * kf[k] * xk[k-2];
          const double f2 = f1*xk[1];
          factor_ee_gl_0[i] = factor_ee_gl_0[i] + f2*dx[0];
          factor_ee_gl_1[i] = factor_ee_gl_1[i] + f2*dx[1];
          factor_ee_gl_2[i] = factor_ee_gl_2[i] + f2*dx[2];
          factor_ee_gl_3[i] = factor_ee_gl_3[i] + f2*dx[3];
          factor_ee_gl_3[i] = factor_ee_gl_3[i] + f1*kf[k-1]*grad_c2;
        }
      }
      touched = true;
    }
    if (!touched) {
      memset(&(factor_ee_gl[nw*4*elec_num]), 0, elec_num*4*sizeof(double));
    }
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_ee_gl (const qmckl_context context,
                                          const int64_t walk_num,
                                          const int64_t elec_num,
                                          const int64_t up_num,
                                          const int64_t bord_num,
                                          const double* b_vector,
                                          const double* ee_distance_rescaled,
                                          const double* ee_distance_rescaled_gl,
                                          const int32_t spin_independent,
                                          double* const factor_ee_gl )
{
#ifdef HAVE_HPC
  return qmckl_compute_jastrow_champ_factor_ee_gl_hpc
#else
    return qmckl_compute_jastrow_champ_factor_ee_gl_doc
#endif
    (context, walk_num, elec_num, up_num, bord_num, b_vector,
     ee_distance_rescaled, ee_distance_rescaled_gl, spin_independent, factor_ee_gl );
}

qmckl_exit_code
qmckl_get_jastrow_champ_asymp_jasb_pderiv(qmckl_context context,
                                          double* const asymp_jasb_pderiv,
                                          const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_jastrow_champ_asymp_jasb_pderiv",
                           NULL);
  }


  /* Provided in finalize_jastrow */
  /*
  qmckl_exit_code rc;
  rc = qmckl_provide_jastrow_champ_asymp_jasb_pderiv(context);
  if(rc != QMCKL_SUCCESS) return rc;
  */

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  const int64_t sze = (ctx->jastrow_champ.bord_num + 1) * 2;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_asymp_jasb_pderiv",
                           "Array too small. Expected 2");
  }
  memcpy(asymp_jasb_pderiv, ctx->jastrow_champ.asymp_jasb_pderiv, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_jastrow_champ_asymp_jasb_pderiv(qmckl_context context)
{

  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_provide_jastrow_champ_asymp_jasb_pderiv",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->jastrow_champ.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_provide_jastrow_champ_asymp_jasb_pderiv",
                           NULL);
  }

  /* Compute if necessary */
  //if (ctx->date > ctx->jastrow_champ.asymp_jasb_pderiv_date) {
  /* Allocate array */
  if (ctx->jastrow_champ.asymp_jasb_pderiv == NULL) {

    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = 2 * (ctx->jastrow_champ.bord_num + 1) * sizeof(double);
    double* asymp_jasb_pderiv = (double*) qmckl_malloc(context, mem_info);

    if (asymp_jasb_pderiv == NULL) {
      return qmckl_failwith( context,
                              QMCKL_ALLOCATION_FAILED,
                              "qmckl_asymp_jasb_pderiv",
                              NULL);
    }
    ctx->jastrow_champ.asymp_jasb_pderiv = asymp_jasb_pderiv;
  }

    rc = qmckl_compute_jastrow_champ_asymp_jasb_pderiv_doc(context,
                                  ctx->jastrow_champ.bord_num,
                                  ctx->jastrow_champ.b_vector,
                                  ctx->jastrow_champ.rescale_factor_ee,
                                  ctx->jastrow_champ.spin_independent,
                                  ctx->jastrow_champ.asymp_jasb_pderiv);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow_champ.asymp_jasb_pderiv_date = ctx->date;
 // }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_jastrow_champ_factor_ee_pderiv(qmckl_context context,
                                         double* const factor_ee_pderiv,
                                         const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_jastrow_champ_factor_ee_pderiv(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (factor_ee_pderiv == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_champ_factor_ee_pderiv",
                           "Null pointer");
  }

  const int64_t sze = ctx->jastrow_champ.bord_num + 1;

  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_factor_ee_pderiv",
                           "Array too small. Expected bord_num+1");
  }

  memcpy(factor_ee_pderiv, ctx->jastrow_champ.factor_ee_pderiv, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_jastrow_champ_factor_ee_pderiv(qmckl_context context)
{

  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_provide_jastrow_champ_factor_ee_pderiv",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->jastrow_champ.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_provide_jastrow_champ_factor_ee_pderiv",
                           NULL);
  }

  /* Check if ee rescaled distance is provided */
  rc = qmckl_provide_ee_distance_rescaled(context);
  if(rc != QMCKL_SUCCESS) return rc;


  /* Provided in finalize_jastrow */
  /*
  rc = qmckl_provide_jastrow_champ_asymp_jasb_pderiv(context);
  if(rc != QMCKL_SUCCESS) return rc;
  */

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow_champ.factor_ee_pderiv_date) {

    /* Allocate array */
    if (ctx->jastrow_champ.factor_ee_pderiv == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = (ctx->jastrow_champ.bord_num + 1) * sizeof(double);
      double* factor_ee_pderiv = (double*) qmckl_malloc(context, mem_info);

      if (factor_ee_pderiv == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_jastrow_champ_factor_ee_pderiv",
                               NULL);
      }
      ctx->jastrow_champ.factor_ee_pderiv = factor_ee_pderiv;
    }

    rc = qmckl_compute_jastrow_champ_factor_ee_pderiv(context,
                                                      ctx->electron.walker.num,
                                                      ctx->electron.num,
                                                      ctx->electron.up_num,
                                                      ctx->jastrow_champ.bord_num,
                                                      ctx->jastrow_champ.b_vector,
                                                      ctx->jastrow_champ.ee_distance_rescaled,
                                                      ctx->jastrow_champ.asymp_jasb_pderiv,
                                                      ctx->jastrow_champ.spin_independent,
                                                      ctx->jastrow_champ.factor_ee_pderiv);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow_champ.factor_ee_pderiv_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_ee_pderiv (const qmckl_context context,
                                              const int64_t walk_num,
                                              const int64_t elec_num,
                                              const int64_t up_num,
                                              const int64_t bord_num,
                                              const double* b_vector,
                                              const double* ee_distance_rescaled,
                                              const double* asymp_jasb_pderiv,
                                              const int32_t spin_independent,
                                              double* const factor_ee_pderiv )
{
    return qmckl_compute_jastrow_champ_factor_ee_pderiv_doc
    (context, walk_num, elec_num, up_num, bord_num, b_vector,
     ee_distance_rescaled, asymp_jasb_pderiv, spin_independent, factor_ee_pderiv );
}

qmckl_exit_code
qmckl_get_jastrow_champ_factor_ee_gl_pderiv(qmckl_context context,
                                    double* const factor_ee_gl_pderiv,
                                    const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_jastrow_champ_factor_ee_gl_pderiv(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (factor_ee_gl_pderiv == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_champ_factor_ee_gl_pderiv",
                           "Null pointer");
  }

  const int64_t sze = 4 * ctx->electron.num * (ctx->jastrow_champ.bord_num+1);

  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_factor_ee_gl_pderiv",
                           "Array too small. Expected 4*elec_num*(bord_num + 1)");
  }

  memcpy(factor_ee_gl_pderiv, ctx->jastrow_champ.factor_ee_gl_pderiv, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_jastrow_champ_factor_ee_gl_pderiv(qmckl_context context)
{

  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_provide_jastrow_champ_factor_ee_gl_pderiv",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->jastrow_champ.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_provide_jastrow_champ_factor_ee_gl_pderiv",
                           NULL);
  }

  /* Check if ee rescaled distance is provided */
  rc = qmckl_provide_ee_distance_rescaled(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Check if ee rescaled distance deriv e is provided */
  rc = qmckl_provide_ee_distance_rescaled_gl(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow_champ.factor_ee_gl_pderiv_date) {

    /* Allocate array */
    if (ctx->jastrow_champ.factor_ee_gl_pderiv == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = 4 * ctx->electron.num * (ctx->jastrow_champ.bord_num + 1) * sizeof(double);
      double* factor_ee_gl_pderiv = (double*) qmckl_malloc(context, mem_info);

      if (factor_ee_gl_pderiv == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_jastrow_champ_factor_ee_gl_pderiv",
                               NULL);
      }
      ctx->jastrow_champ.factor_ee_gl_pderiv = factor_ee_gl_pderiv;
    }

    rc = qmckl_compute_jastrow_champ_factor_ee_gl_pderiv(context,
                                                  ctx->electron.walker.num,
                                                  ctx->electron.num,
                                                  ctx->electron.up_num,
                                                  ctx->jastrow_champ.bord_num,
                                                  ctx->jastrow_champ.b_vector,
                                                  ctx->jastrow_champ.ee_distance_rescaled,
                                                  ctx->jastrow_champ.ee_distance_rescaled_gl,
                                                  ctx->jastrow_champ.spin_independent,
                                                  ctx->jastrow_champ.factor_ee_gl_pderiv);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow_champ.factor_ee_gl_pderiv_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_ee_gl_pderiv (const qmckl_context context,
                                          const int64_t walk_num,
                                          const int64_t elec_num,
                                          const int64_t up_num,
                                          const int64_t bord_num,
                                          const double* b_vector,
                                          const double* ee_distance_rescaled,
                                          const double* ee_distance_rescaled_gl,
                                          const int32_t spin_independent,
                                          double* const factor_ee_gl_pderiv )
{
    return qmckl_compute_jastrow_champ_factor_ee_gl_pderiv_doc
    (context, walk_num, elec_num, up_num, bord_num, b_vector,
     ee_distance_rescaled, ee_distance_rescaled_gl, spin_independent, factor_ee_gl_pderiv );
}

qmckl_exit_code
qmckl_get_jastrow_champ_asymp_jasa(qmckl_context context,
                                   double* const asymp_jasa,
                                   const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_jastrow_champ_asymp_jasa",
                           NULL);
  }


  /* Provided in finalize_jastrow */
  /*
  qmckl_exit_code rc;
  rc = qmckl_provide_jastrow_champ_asymp_jasa(context);
  if(rc != QMCKL_SUCCESS) return rc;
  */

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (asymp_jasa == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_champ_asymp_jasa",
                           "Null pointer");
  }

  const int64_t sze = ctx->jastrow_champ.type_nucl_num;

  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_asymp_jasa",
                           "Array too small. Expected nucl_num");
  }

  memcpy(asymp_jasa, ctx->jastrow_champ.asymp_jasa, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_jastrow_champ_asymp_jasa(qmckl_context context)
{

  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_provide_jastrow_champ_asymp_jasa",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->jastrow_champ.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_provide_jastrow_champ_asymp_jasa",
                           NULL);
  }

//  /* Compute if necessary */
//  if (ctx->date > ctx->jastrow_champ.asymp_jasa_date) {

    /* Allocate array */
    if (ctx->jastrow_champ.asymp_jasa == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->jastrow_champ.type_nucl_num * sizeof(double);
      double* asymp_jasa = (double*) qmckl_malloc(context, mem_info);

      if (asymp_jasa == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_asymp_jasa",
                               NULL);
      }
      ctx->jastrow_champ.asymp_jasa = asymp_jasa;
    }

    rc = qmckl_compute_jastrow_champ_asymp_jasa(context,
                                  ctx->jastrow_champ.aord_num,
                                  ctx->jastrow_champ.type_nucl_num,
                                  ctx->jastrow_champ.a_vector,
                                  ctx->jastrow_champ.rescale_factor_en,
                                  ctx->jastrow_champ.asymp_jasa);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow_champ.asymp_jasa_date = ctx->date;
//  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_jastrow_champ_en_distance_rescaled(qmckl_context context,
                                             double* const distance_rescaled,
                                             const int64_t size_max)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_en_distance_rescaled(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (distance_rescaled == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_champ_en_distance_rescaled",
                           "Null pointer");
  }

  const int64_t sze = ctx->electron.num * ctx->nucleus.num * ctx->electron.walker.num;

  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_en_distance_rescaled",
                           "Array too small. Expected walk_num*nucl_num*elec_num");
  }

  memcpy(distance_rescaled, ctx->jastrow_champ.en_distance_rescaled, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_en_distance_rescaled(qmckl_context context)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!(ctx->nucleus.provided)) {
    return QMCKL_NOT_PROVIDED;
  }

  /* Compute if necessary */
  if (ctx->electron.walker.point.date > ctx->jastrow_champ.en_distance_rescaled_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      if (ctx->jastrow_champ.en_distance_rescaled != NULL) {
        qmckl_exit_code rc = qmckl_free(context, ctx->jastrow_champ.en_distance_rescaled);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_en_distance_rescaled",
                                 "Unable to free ctx->jastrow_champ.en_distance_rescaled");
        }
        ctx->jastrow_champ.en_distance_rescaled = NULL;
      }
    }

    /* Allocate array */
    if (ctx->jastrow_champ.en_distance_rescaled == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.num * ctx->nucleus.num *
        ctx->electron.walker.num * sizeof(double);
      double* en_distance_rescaled = (double*) qmckl_malloc(context, mem_info);

      if (en_distance_rescaled == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_en_distance_rescaled",
                               NULL);
      }
      ctx->jastrow_champ.en_distance_rescaled = en_distance_rescaled;
    }

    qmckl_exit_code rc =
      qmckl_compute_en_distance_rescaled(context,
                                ctx->electron.num,
                                ctx->nucleus.num,
                                ctx->jastrow_champ.type_nucl_num,
                                ctx->jastrow_champ.type_nucl_vector,
                                ctx->jastrow_champ.rescale_factor_en,
                                ctx->electron.walker.num,
                                ctx->electron.walker.point.coord.data,
                                ctx->nucleus.coord.data,
                                ctx->jastrow_champ.en_distance_rescaled);

    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow_champ.en_distance_rescaled_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_compute_en_distance_rescaled_hpc (
          const qmckl_context context,
          const int64_t elec_num,
          const int64_t nucl_num,
          const int64_t type_nucl_num,
          const int64_t* type_nucl_vector,
          const double* rescale_factor_en,
          const int64_t walk_num,
          const double* elec_coord,
          const double* nucl_coord,
          double* const en_distance_rescaled )
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) return QMCKL_INVALID_CONTEXT;
  if (elec_num <= 0) return QMCKL_INVALID_ARG_2;
  if (nucl_num <= 0) return QMCKL_INVALID_ARG_3;
  if (type_nucl_num <= 0) return QMCKL_INVALID_ARG_4;
  if (type_nucl_vector == NULL) return QMCKL_INVALID_ARG_5;
  if (rescale_factor_en == NULL) return QMCKL_INVALID_ARG_6;
  if (walk_num <= 0) return QMCKL_INVALID_ARG_7;
  if (elec_coord == NULL) return QMCKL_INVALID_ARG_8;
  if (nucl_coord == NULL) return QMCKL_INVALID_ARG_9;
  if (en_distance_rescaled == NULL) return QMCKL_INVALID_ARG_10;

  const int64_t sze = elec_num*walk_num;

  qmckl_exit_code result = QMCKL_SUCCESS;
#ifdef HAVE_OPENMP
  #pragma omp parallel
#endif
  {
    qmckl_exit_code rc = QMCKL_SUCCESS;
#ifdef HAVE_OPENMP
    #pragma omp for
#endif
    for (int64_t k=0 ; k<walk_num ; ++k)
    {
      for (int64_t a=0 ; a<nucl_num ; ++a) {
        const double coord[3] = { nucl_coord[a], nucl_coord[a+nucl_num], nucl_coord[a+2*nucl_num] };
        rc |= qmckl_distance_rescaled(context, 'T', 'N', elec_num, 1,
                                      &(elec_coord[k*elec_num]), sze,
                                      coord, 3,
                                      &(en_distance_rescaled[elec_num*(a+nucl_num*k)]), elec_num,
                                      rescale_factor_en[type_nucl_vector[a]]);
      }
    }
#ifdef HAVE_OPENMP
    #pragma omp critical
#endif
    result |= rc;
  }
  return result;
}

qmckl_exit_code qmckl_compute_en_distance_rescaled (
          const qmckl_context context,
          const int64_t elec_num,
          const int64_t nucl_num,
          const int64_t type_nucl_num,
          const int64_t* type_nucl_vector,
          const double* rescale_factor_en,
          const int64_t walk_num,
          const double* elec_coord,
          const double* nucl_coord,
          double* const en_distance_rescaled )
{
#ifdef HAVE_HPC
  return qmckl_compute_en_distance_rescaled_hpc
#else
  return qmckl_compute_en_distance_rescaled_doc
#endif
  (context, elec_num, nucl_num, type_nucl_num, type_nucl_vector,
  rescale_factor_en, walk_num, elec_coord, nucl_coord, en_distance_rescaled );
}

qmckl_exit_code
qmckl_get_jastrow_champ_en_distance_rescaled_gl(qmckl_context context,
                                                double* const distance_rescaled_gl,
                                                const int64_t size_max)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_en_distance_rescaled_gl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (distance_rescaled_gl == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_champ_en_distance_rescaled_gl",
                           "Null pointer");
  }

  const int64_t sze = 4 * ctx->electron.num * ctx->nucleus.num * ctx->electron.walker.num;

  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_en_distance_rescaled_gl",
                           "Array too small. Expected walk_num*elec_num*4");
  }

  memcpy(distance_rescaled_gl, ctx->jastrow_champ.en_distance_rescaled_gl, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_en_distance_rescaled_gl(qmckl_context context)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!(ctx->nucleus.provided)) {
    return QMCKL_NOT_PROVIDED;
  }

  /* Compute if necessary */
  if (ctx->electron.walker.point.date > ctx->jastrow_champ.en_distance_rescaled_gl_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      if (ctx->jastrow_champ.en_distance_rescaled_gl != NULL) {
        qmckl_exit_code rc = qmckl_free(context, ctx->jastrow_champ.en_distance_rescaled_gl);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_en_distance_rescaled_gl",
                                 "Unable to free ctx->jastrow_champ.en_distance_rescaled_gl");
        }
        ctx->jastrow_champ.en_distance_rescaled_gl = NULL;
      }
    }

    /* Allocate array */
    if (ctx->jastrow_champ.en_distance_rescaled_gl == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = 4 * ctx->electron.num * ctx->nucleus.num *
        ctx->electron.walker.num * sizeof(double);
      double* en_distance_rescaled_gl = (double*) qmckl_malloc(context, mem_info);

      if (en_distance_rescaled_gl == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_en_distance_rescaled_gl",
                               NULL);
      }
      ctx->jastrow_champ.en_distance_rescaled_gl = en_distance_rescaled_gl;
    }

    qmckl_exit_code rc =
      qmckl_compute_en_distance_rescaled_gl(context,
                                ctx->electron.num,
                                ctx->nucleus.num,
                                ctx->jastrow_champ.type_nucl_num,
                                ctx->jastrow_champ.type_nucl_vector,
                                ctx->jastrow_champ.rescale_factor_en,
                                ctx->electron.walker.num,
                                ctx->electron.walker.point.coord.data,
                                ctx->nucleus.coord.data,
                                ctx->jastrow_champ.en_distance_rescaled_gl);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow_champ.en_distance_rescaled_gl_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_compute_en_distance_rescaled_gl_hpc (
          const qmckl_context context,
          const int64_t elec_num,
          const int64_t nucl_num,
          const int64_t type_nucl_num,
          int64_t* const type_nucl_vector,
          const double*  rescale_factor_en,
          const int64_t walk_num,
          const double* elec_coord,
          const double* nucl_coord,
          double* const en_distance_rescaled_gl )
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) return QMCKL_INVALID_CONTEXT;
  if (elec_num <= 0) return QMCKL_INVALID_ARG_2;
  if (nucl_num <= 0) return QMCKL_INVALID_ARG_3;
  if (type_nucl_num <= 0) return QMCKL_INVALID_ARG_4;
  if (type_nucl_vector == NULL) return QMCKL_INVALID_ARG_5;
  if (rescale_factor_en == NULL) return QMCKL_INVALID_ARG_6;
  if (walk_num <= 0) return QMCKL_INVALID_ARG_7;
  if (elec_coord == NULL) return QMCKL_INVALID_ARG_8;
  if (nucl_coord == NULL) return QMCKL_INVALID_ARG_9;
  if (en_distance_rescaled_gl == NULL) return QMCKL_INVALID_ARG_10;

  const int64_t sze = elec_num*walk_num;

  qmckl_exit_code result = QMCKL_SUCCESS;
#ifdef HAVE_OPENMP
  #pragma omp parallel
#endif
  {
    qmckl_exit_code rc = QMCKL_SUCCESS;
#ifdef HAVE_OPENMP
    #pragma omp for
#endif
    for (int64_t k=0 ; k<walk_num ; ++k)
    {
      for (int64_t a=0 ; a<nucl_num ; ++a) {
        const double coord[3] = { nucl_coord[a], nucl_coord[a+nucl_num], nucl_coord[a+2*nucl_num] };
        rc |= qmckl_distance_rescaled_gl(context, 'T', 'T', elec_num, 1,
                                         &(elec_coord[k*elec_num]), sze,
                                         coord, 1,
                                         &(en_distance_rescaled_gl[4*elec_num*(a+nucl_num*k)]), elec_num,
                                         rescale_factor_en[type_nucl_vector[a]]);
      }
    }
#ifdef HAVE_OPENMP
    #pragma omp critical
#endif
    result |= rc;
  }
  return result;
}

qmckl_exit_code qmckl_compute_en_distance_rescaled_gl (
          const qmckl_context context,
          const int64_t elec_num,
          const int64_t nucl_num,
          const int64_t type_nucl_num,
          int64_t* const type_nucl_vector,
          const double*  rescale_factor_en,
          const int64_t walk_num,
          const double* elec_coord,
          const double* nucl_coord,
          double* const en_distance_rescaled_gl )
{
#ifdef HAVE_HPC
  return qmckl_compute_en_distance_rescaled_gl_hpc
#else
  return qmckl_compute_en_distance_rescaled_gl_doc
#endif
    (context, elec_num, nucl_num, type_nucl_num, type_nucl_vector, rescale_factor_en,
     walk_num, elec_coord, nucl_coord, en_distance_rescaled_gl );
}

qmckl_exit_code
qmckl_get_jastrow_champ_factor_en(qmckl_context context,
                            double* const factor_en,
                            const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_jastrow_champ_factor_en",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  qmckl_exit_code rc;

  rc = qmckl_provide_jastrow_champ_factor_en(context);
  if (rc != QMCKL_SUCCESS) return rc;

  if (factor_en == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_champ_factor_en",
                           "Null pointer");
  }

  const int64_t sze=ctx->electron.walker.num;

  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_factor_en",
                           "Array too small. Expected walk_num");
  }

  memcpy(factor_en, ctx->jastrow_champ.factor_en, sze*sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_jastrow_champ_factor_en(qmckl_context context)
{

  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_provide_jastrow_champ_factor_en",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->jastrow_champ.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_provide_jastrow_champ_factor_en",
                           NULL);
  }

  /* Check if en rescaled distance is provided */
  rc = qmckl_provide_en_distance_rescaled(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Provided in finalize_jastrow */
  /*
  rc = qmckl_provide_jastrow_champ_asymp_jasa(context);
  if(rc != QMCKL_SUCCESS) return rc;
  */

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow_champ.factor_en_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      if (ctx->jastrow_champ.factor_en != NULL) {
        rc = qmckl_free(context, ctx->jastrow_champ.factor_en);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_jastrow_champ_factor_en",
                                 "Unable to free ctx->jastrow_champ.factor_en");
        }
        ctx->jastrow_champ.factor_en = NULL;
      }
    }

    /* Allocate array */
    if (ctx->jastrow_champ.factor_en == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.walker.num * sizeof(double);
      double* factor_en = (double*) qmckl_malloc(context, mem_info);

      if (factor_en == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_jastrow_champ_factor_en",
                               NULL);
      }
      ctx->jastrow_champ.factor_en = factor_en;
    }

    rc = qmckl_compute_jastrow_champ_factor_en(context,
                                 ctx->electron.walker.num,
                                 ctx->electron.num,
                                 ctx->nucleus.num,
                                 ctx->jastrow_champ.type_nucl_num,
                                 ctx->jastrow_champ.type_nucl_vector,
                                 ctx->jastrow_champ.aord_num,
                                 ctx->jastrow_champ.a_vector,
                                 ctx->jastrow_champ.en_distance_rescaled,
                                 ctx->jastrow_champ.asymp_jasa,
                                 ctx->jastrow_champ.factor_en);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow_champ.factor_en_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_compute_jastrow_champ_factor_en (
         const qmckl_context context,
         const int64_t walk_num,
         const int64_t elec_num,
         const int64_t nucl_num,
         const int64_t type_nucl_num,
         const int64_t* type_nucl_vector,
         const int64_t aord_num,
         const double* a_vector,
         const double* en_distance_rescaled,
         const double* asymp_jasa,
         double* const factor_en )
{
#ifdef HAVE_HPC
  return qmckl_compute_jastrow_champ_factor_en_hpc
#else
  return qmckl_compute_jastrow_champ_factor_en_doc
#endif
    (context, walk_num, elec_num, nucl_num, type_nucl_num,
     type_nucl_vector, aord_num, a_vector, en_distance_rescaled,
     asymp_jasa, factor_en );
}

qmckl_exit_code qmckl_compute_jastrow_champ_factor_en_hpc (
         const qmckl_context context,
         const int64_t walk_num,
         const int64_t elec_num,
         const int64_t nucl_num,
         const int64_t type_nucl_num,
         const int64_t* type_nucl_vector,
         const int64_t aord_num,
         const double* a_vector,
         const double* en_distance_rescaled,
         const double* asymp_jasa,
         double* const factor_en )
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) return QMCKL_NULL_CONTEXT;
  if (walk_num <= 0) return QMCKL_INVALID_ARG_3;
  if (elec_num <= 0) return QMCKL_INVALID_ARG_3;
  if (nucl_num <= 0) return QMCKL_INVALID_ARG_4;
  if (type_nucl_num <= 0) return QMCKL_INVALID_ARG_5;
  if (type_nucl_vector == NULL) return QMCKL_INVALID_ARG_6;
  if (aord_num < 0) return QMCKL_INVALID_ARG_7;
  if (a_vector == NULL) return QMCKL_INVALID_ARG_8;
  if (en_distance_rescaled == NULL) return QMCKL_INVALID_ARG_9;
  if (asymp_jasa == NULL) return QMCKL_INVALID_ARG_10;
  if (factor_en == NULL) return QMCKL_INVALID_ARG_11;

  const double de = (double) elec_num;
#ifdef HAVE_OPENMP
  #pragma omp parallel for
#endif
  for (int64_t nw=0 ; nw<walk_num ; ++nw) {
    factor_en[nw] = 0.;
    const double* en_distance_rescaled_ = &(en_distance_rescaled[nw*elec_num*nucl_num]);

    for (int64_t a=0 ; a<nucl_num ; ++a) {
      const double* en_distance_rescaled__ = &(en_distance_rescaled_[a*elec_num]);
      const double* a_vec = &(a_vector[(aord_num+1)*type_nucl_vector[a]]);

      factor_en[nw] = factor_en[nw] - asymp_jasa[type_nucl_vector[a]]*de;

      for (int64_t i=0 ; i<elec_num ; ++i) {
        double x = en_distance_rescaled__[i];
        factor_en[nw] = factor_en[nw] + a_vec[0]*x / (1.0 + a_vec[1]*x);

        for (int64_t p=2 ; p <= aord_num ; ++p) {
          x *= en_distance_rescaled__[i];
          factor_en[nw] = factor_en[nw] + a_vec[p]*x;
        }
      }
    }

  }
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_jastrow_champ_factor_en_gl(qmckl_context context,
                                    double* const factor_en_gl,
                                    const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_jastrow_champ_factor_en_gl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (factor_en_gl == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_champ_factor_en_gl",
                           "Null pointer");
  }

  const int64_t sze = ctx->electron.walker.num * 4 * ctx->electron.num;

  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_factor_en_gl",
                           "Array too small. Expected 4*walk_num*elec_num");
  }

  memcpy(factor_en_gl, ctx->jastrow_champ.factor_en_gl, sze*sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_jastrow_champ_factor_en_gl(qmckl_context context)
{

  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_provide_jastrow_champ_factor_en_gl",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->jastrow_champ.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_provide_jastrow_champ_factor_en_gl",
                           NULL);
  }

  /* Check if en rescaled distance is provided */
  rc = qmckl_provide_en_distance_rescaled(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Check if en rescaled distance derivatives is provided */
  rc = qmckl_provide_en_distance_rescaled_gl(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow_champ.factor_en_gl_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      if (ctx->jastrow_champ.factor_en_gl != NULL) {
        rc = qmckl_free(context, ctx->jastrow_champ.factor_en_gl);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_jastrow_champ_factor_en_gl",
                                 "Unable to free ctx->jastrow_champ.factor_en_gl");
        }
        ctx->jastrow_champ.factor_en_gl = NULL;
      }
    }

    /* Allocate array */
    if (ctx->jastrow_champ.factor_en_gl == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.walker.num * 4 * ctx->electron.num * sizeof(double);
      double* factor_en_gl = (double*) qmckl_malloc(context, mem_info);

      if (factor_en_gl == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_jastrow_champ_factor_en_gl",
                               NULL);
      }
      ctx->jastrow_champ.factor_en_gl = factor_en_gl;
    }

    rc = qmckl_compute_jastrow_champ_factor_en_gl(context,
                                         ctx->electron.walker.num,
                                         ctx->electron.num,
                                         ctx->nucleus.num,
                                         ctx->jastrow_champ.type_nucl_num,
                                         ctx->jastrow_champ.type_nucl_vector,
                                         ctx->jastrow_champ.aord_num,
                                         ctx->jastrow_champ.a_vector,
                                         ctx->jastrow_champ.en_distance_rescaled,
                                         ctx->jastrow_champ.en_distance_rescaled_gl,
                                         ctx->jastrow_champ.factor_en_gl);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow_champ.factor_en_gl_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_en_gl_hpc (const qmckl_context context,
                                              const int64_t walk_num,
                                              const int64_t elec_num,
                                              const int64_t nucl_num,
                                              const int64_t type_nucl_num,
                                              const int64_t* type_nucl_vector,
                                              const int64_t aord_num,
                                              const double* a_vector,
                                              const double* en_distance_rescaled,
                                              const double* en_distance_rescaled_gl,
                                              double* const factor_en_gl )
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) return QMCKL_NULL_CONTEXT;
  if (walk_num <= 0) return QMCKL_INVALID_ARG_3;
  if (elec_num <= 0) return QMCKL_INVALID_ARG_3;
  if (nucl_num <= 0) return QMCKL_INVALID_ARG_4;
  if (type_nucl_num <= 0) return QMCKL_INVALID_ARG_5;
  if (type_nucl_vector == NULL) return QMCKL_INVALID_ARG_6;
  if (aord_num < 0) return QMCKL_INVALID_ARG_7;
  if (a_vector == NULL) return QMCKL_INVALID_ARG_8;
  if (en_distance_rescaled == NULL) return QMCKL_INVALID_ARG_9;
  if (en_distance_rescaled_gl == NULL) return QMCKL_INVALID_ARG_10;
  if (factor_en_gl == NULL) return QMCKL_INVALID_ARG_11;

  double kf[aord_num+1];
  for (int k=0 ; k<=aord_num ; ++k) {
    kf[k] = (double) k;
  }

#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
  for (int64_t nw = 0; nw < walk_num; ++nw) {
    bool touched = false;

    for (int64_t a = 0; a < nucl_num; ++a) {
      const double* dxj = &en_distance_rescaled_gl[4*elec_num*(a+nw*nucl_num)];
      const double*  xj = &en_distance_rescaled   [  elec_num*(a+nw*nucl_num)];
      const double* a_vec = &( a_vector[(aord_num+1)*type_nucl_vector[a]] );

      double * restrict factor_en_gl_0 = &(factor_en_gl[nw*elec_num*4]);
      double * restrict factor_en_gl_1 = factor_en_gl_0 + elec_num;
      double * restrict factor_en_gl_2 = factor_en_gl_1 + elec_num;
      double * restrict factor_en_gl_3 = factor_en_gl_2 + elec_num;

      for (int64_t i = 0; i < elec_num; ++i) {

        double x = xj[i];
        if (x < 1.e-12) continue;

        const double denom      = 1.0 + a_vec[1]*x;
        const double invdenom   = 1.0  / denom;
        const double invdenom2  = invdenom * invdenom;

        const double* restrict dx = dxj + 4*i;

        const double grad_c2 = dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2];

        double f = a_vec[0] * invdenom2;

        if (touched) {
          factor_en_gl_0[i] = factor_en_gl_0[i] + f*dx[0];
          factor_en_gl_1[i] = factor_en_gl_1[i] + f*dx[1];
          factor_en_gl_2[i] = factor_en_gl_2[i] + f*dx[2];
          factor_en_gl_3[i] = factor_en_gl_3[i] + f*dx[3];
        } else {
          factor_en_gl_0[i] = f*dx[0];
          factor_en_gl_1[i] = f*dx[1];
          factor_en_gl_2[i] = f*dx[2];
          factor_en_gl_3[i] = f*dx[3];
        }

        factor_en_gl_3[i] = factor_en_gl_3[i] - f*grad_c2*invdenom*2.0 * a_vec[1];

        double xk[aord_num+1];
        xk[0] = 1.0;
        for (int k=1 ; k<= aord_num ; ++k) {
          xk[k] = xk[k-1]*x;
        }

        for (int k=2 ; k<= aord_num ; ++k) {
          const double f1 = a_vec[k] * kf[k] * xk[k-2];
          const double f2 = f1*xk[1];
          factor_en_gl_0[i] = factor_en_gl_0[i] + f2*dx[0];
          factor_en_gl_1[i] = factor_en_gl_1[i] + f2*dx[1];
          factor_en_gl_2[i] = factor_en_gl_2[i] + f2*dx[2];
          factor_en_gl_3[i] = factor_en_gl_3[i] + f2*dx[3];
          factor_en_gl_3[i] = factor_en_gl_3[i] + f1*kf[k-1]*grad_c2;
        }
      }
      touched = true;
    }
    if (!touched) {
      memset(&(factor_en_gl[nw*4*elec_num]), 0, elec_num*4*sizeof(double));
    }
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_en_gl (const qmckl_context context,
                                          const int64_t walk_num,
                                          const int64_t elec_num,
                                          const int64_t nucl_num,
                                          const int64_t type_nucl_num,
                                          const int64_t* type_nucl_vector,
                                          const int64_t aord_num,
                                          const double* a_vector,
                                          const double* en_distance_rescaled,
                                          const double* en_distance_rescaled_gl,
                                          double* const factor_en_gl )
{
#ifdef HAVE_HPC
  return qmckl_compute_jastrow_champ_factor_en_gl_hpc
#else
  return qmckl_compute_jastrow_champ_factor_en_gl_doc
#endif
    (context, walk_num, elec_num, nucl_num, type_nucl_num, type_nucl_vector, aord_num,
     a_vector, en_distance_rescaled, en_distance_rescaled_gl, factor_en_gl );
}

qmckl_exit_code
qmckl_get_jastrow_champ_asymp_jasa_pderiv(qmckl_context context,
                                          double* const asymp_jasa_pderiv,
                                          const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_jastrow_champ_asymp_jasa_pderiv",
                           NULL);
  }


  /* Provided in finalize_jastrow */
  /*
  qmckl_exit_code rc;
  rc = qmckl_provide_jastrow_champ_asymp_jasa_pderiv(context);
  if(rc != QMCKL_SUCCESS) return rc;
  */

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (asymp_jasa_pderiv == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_champ_asymp_jasa_pderiv",
                           "Null pointer");
  }

  const int64_t sze = ctx->jastrow_champ.type_nucl_num * (ctx->jastrow_champ.aord_num + 1);

  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_asymp_jasa_pderiv",
                           "Array too small. Expected nucl_num");
  }

  memcpy(asymp_jasa_pderiv, ctx->jastrow_champ.asymp_jasa_pderiv, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_jastrow_champ_asymp_jasa_pderiv(qmckl_context context)
{

  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_provide_jastrow_champ_asymp_jasa_pderiv",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->jastrow_champ.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_provide_jastrow_champ_asymp_jasa_pderiv",
                           NULL);
  }

//  /* Compute if necessary */
//  if (ctx->date > ctx->jastrow_champ.asymp_jasa_pderiv_date) {

    /* Allocate array */
    if (ctx->jastrow_champ.asymp_jasa_pderiv == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->jastrow_champ.type_nucl_num * (ctx->jastrow_champ.aord_num + 1) * sizeof(double);
      double* asymp_jasa_pderiv = (double*) qmckl_malloc(context, mem_info);

      if (asymp_jasa_pderiv == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_asymp_jasa_pderiv",
                               NULL);
      }
      ctx->jastrow_champ.asymp_jasa_pderiv = asymp_jasa_pderiv;
    }

    rc = qmckl_compute_jastrow_champ_asymp_jasa_pderiv(context,
                                  ctx->jastrow_champ.aord_num,
                                  ctx->jastrow_champ.type_nucl_num,
                                  ctx->jastrow_champ.a_vector,
                                  ctx->jastrow_champ.rescale_factor_en,
                                  ctx->jastrow_champ.asymp_jasa_pderiv);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow_champ.asymp_jasa_pderiv_date = ctx->date;
//  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_jastrow_champ_factor_en_pderiv(qmckl_context context,
                                         double* const factor_en_pderiv,
                                         const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_jastrow_champ_factor_en_pderiv",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  qmckl_exit_code rc;

  rc = qmckl_provide_jastrow_champ_factor_en_pderiv(context);
  if (rc != QMCKL_SUCCESS) return rc;

  if (factor_en_pderiv == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_champ_factor_en_pderiv",
                           "Null pointer");
  }

  const int64_t sze=ctx->jastrow_champ.type_nucl_num * (ctx->jastrow_champ.aord_num + 1);

  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_factor_en_pderiv",
                           "Array too small. Expected walk_num");
  }

  memcpy(factor_en_pderiv, ctx->jastrow_champ.factor_en_pderiv, sze*sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_jastrow_champ_factor_en_pderiv(qmckl_context context)
{

  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_provide_jastrow_champ_factor_en_pderiv",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->jastrow_champ.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_provide_jastrow_champ_factor_en_pderiv",
                           NULL);
  }

  /* Check if en rescaled distance is provided */
  rc = qmckl_provide_en_distance_rescaled(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Provided in finalize_jastrow */
  /*
  rc = qmckl_provide_jastrow_champ_asymp_jasa_pderiv(context);
  if(rc != QMCKL_SUCCESS) return rc;
  */

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow_champ.factor_en_pderiv_date) {

    /* Allocate array */
    if (ctx->jastrow_champ.factor_en_pderiv == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->jastrow_champ.type_nucl_num * (ctx->jastrow_champ.aord_num + 1) * sizeof(double);
      double* factor_en_pderiv = (double*) qmckl_malloc(context, mem_info);

      if (factor_en_pderiv == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_jastrow_champ_factor_en_pderiv",
                               NULL);
      }
      ctx->jastrow_champ.factor_en_pderiv = factor_en_pderiv;
    }

    rc = qmckl_compute_jastrow_champ_factor_en_pderiv(context,
                                 ctx->electron.walker.num,
                                 ctx->electron.num,
                                 ctx->nucleus.num,
                                 ctx->jastrow_champ.type_nucl_num,
                                 ctx->jastrow_champ.type_nucl_vector,
                                 ctx->jastrow_champ.aord_num,
                                 ctx->jastrow_champ.a_vector,
                                 ctx->jastrow_champ.en_distance_rescaled,
                                 ctx->jastrow_champ.asymp_jasa_pderiv,
                                 ctx->jastrow_champ.factor_en_pderiv);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow_champ.factor_en_pderiv_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_compute_jastrow_champ_factor_en_pderiv (
         const qmckl_context context,
         const int64_t walk_num,
         const int64_t elec_num,
         const int64_t nucl_num,
         const int64_t type_nucl_num,
         const int64_t* type_nucl_vector,
         const int64_t aord_num,
         const double* a_vector,
         const double* en_distance_rescaled,
         const double* asymp_jasa_pderiv,
         double* const factor_en_pderiv )
{
  return qmckl_compute_jastrow_champ_factor_en_pderiv_doc
    (context, walk_num, elec_num, nucl_num, type_nucl_num,
     type_nucl_vector, aord_num, a_vector, en_distance_rescaled,
     asymp_jasa_pderiv, factor_en_pderiv );
}

qmckl_exit_code
qmckl_get_jastrow_champ_factor_en_gl_pderiv(qmckl_context context,
                                    double* const factor_en_gl_pderiv,
                                    const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_jastrow_champ_factor_en_gl_pderiv(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (factor_en_gl_pderiv == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_champ_factor_en_gl_pderiv",
                           "Null pointer");
  }

  const int64_t sze = 4 * ctx->electron.num * (ctx->jastrow_champ.aord_num + 1) * ctx->jastrow_champ.type_nucl_num;

  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_factor_en_gl_pderiv",
                           "Array too small. Expected 4*walk_num*elec_num");
  }

  memcpy(factor_en_gl_pderiv, ctx->jastrow_champ.factor_en_gl_pderiv, sze*sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_jastrow_champ_factor_en_gl_pderiv(qmckl_context context)
{

  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_provide_jastrow_champ_factor_en_gl_pderiv",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->jastrow_champ.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_provide_jastrow_champ_factor_en_gl_pderiv",
                           NULL);
  }

  /* Check if en rescaled distance is provided */
  rc = qmckl_provide_en_distance_rescaled(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Check if en rescaled distance derivatives is provided */
  rc = qmckl_provide_en_distance_rescaled_gl(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow_champ.factor_en_gl_pderiv_date) {

    /* Allocate array */
    if (ctx->jastrow_champ.factor_en_gl_pderiv == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = 4 * ctx->electron.num * (ctx->jastrow_champ.aord_num + 1) * ctx->jastrow_champ.type_nucl_num * sizeof(double);
      double* factor_en_gl_pderiv = (double*) qmckl_malloc(context, mem_info);

      if (factor_en_gl_pderiv == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_jastrow_champ_factor_en_gl_pderiv",
                               NULL);
      }
      ctx->jastrow_champ.factor_en_gl_pderiv = factor_en_gl_pderiv;
    }

    rc = qmckl_compute_jastrow_champ_factor_en_gl_pderiv(context,
                                         ctx->electron.walker.num,
                                         ctx->electron.num,
                                         ctx->nucleus.num,
                                         ctx->jastrow_champ.type_nucl_num,
                                         ctx->jastrow_champ.type_nucl_vector,
                                         ctx->jastrow_champ.aord_num,
                                         ctx->jastrow_champ.a_vector,
                                         ctx->jastrow_champ.en_distance_rescaled,
                                         ctx->jastrow_champ.en_distance_rescaled_gl,
                                         ctx->jastrow_champ.factor_en_gl_pderiv);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow_champ.factor_en_gl_pderiv_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_en_gl_pderiv (const qmckl_context context,
                                          const int64_t walk_num,
                                          const int64_t elec_num,
                                          const int64_t nucl_num,
                                          const int64_t type_nucl_num,
                                          const int64_t* type_nucl_vector,
                                          const int64_t aord_num,
                                          const double* a_vector,
                                          const double* en_distance_rescaled,
                                          const double* en_distance_rescaled_gl,
                                          double* const factor_en_gl_pderiv )
{
  return qmckl_compute_jastrow_champ_factor_en_gl_pderiv_doc
    (context, walk_num, elec_num, nucl_num, type_nucl_num, type_nucl_vector, aord_num,
     a_vector, en_distance_rescaled, en_distance_rescaled_gl, factor_en_gl_pderiv );
}

qmckl_exit_code
qmckl_get_jastrow_champ_een_rescaled_e(qmckl_context context,
                                                double* const een_rescaled_e,
                                                const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_een_rescaled_e(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (een_rescaled_e == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_champ_een_rescaled_e",
                           "Null pointer");
  }

  const int64_t sze = ctx->electron.num * ctx->electron.num * ctx->electron.walker.num * (ctx->jastrow_champ.cord_num + 1);

  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_een_rescaled_e",
                           "Array too small. Expected elec_num*elec_num*walk_num*(cord_num + 1)");
  }

  memcpy(een_rescaled_e, ctx->jastrow_champ.een_rescaled_e, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_een_rescaled_e(qmckl_context context)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  /* Check if ee distance is provided */
  qmckl_exit_code rc = qmckl_provide_ee_distance(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow_champ.een_rescaled_e_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      if (ctx->jastrow_champ.een_rescaled_e != NULL) {
        rc = qmckl_free(context, ctx->jastrow_champ.een_rescaled_e);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_een_rescaled_e",
                                 "Unable to free ctx->jastrow_champ.een_rescaled_e");
        }
        ctx->jastrow_champ.een_rescaled_e = NULL;
      }
    }

    /* Allocate array */
    if (ctx->jastrow_champ.een_rescaled_e == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.num * ctx->electron.num *
        ctx->electron.walker.num * (ctx->jastrow_champ.cord_num + 1) * sizeof(double);
      double* een_rescaled_e = (double*) qmckl_malloc(context, mem_info);

      if (een_rescaled_e == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_een_rescaled_e",
                               NULL);
      }
      ctx->jastrow_champ.een_rescaled_e = een_rescaled_e;
    }

    rc = qmckl_compute_een_rescaled_e(context,
                                      ctx->electron.walker.num,
                                      ctx->electron.num,
                                      ctx->jastrow_champ.cord_num,
                                      ctx->jastrow_champ.rescale_factor_ee,
                                      ctx->electron.ee_distance,
                                      ctx->jastrow_champ.een_rescaled_e);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow_champ.een_rescaled_e_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_compute_een_rescaled_e_hpc (const qmckl_context context,
                                                  const int64_t walk_num,
                                                  const int64_t elec_num,
                                                  const int64_t cord_num,
                                                  const double rescale_factor_ee,
                                                  const double* ee_distance,
                                                  double* const een_rescaled_e ) {

  if (context == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  }

  if (walk_num <= 0) {
    return QMCKL_INVALID_ARG_2;
  }

  if (elec_num <= 0) {
    return QMCKL_INVALID_ARG_3;
  }

  if (cord_num < 0) {
    return QMCKL_INVALID_ARG_4;
  }

#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
  for (int64_t nw=0 ; nw<walk_num ; ++nw) {
    double* const restrict een_rescaled_e_ =
      &een_rescaled_e[elec_num*elec_num*((cord_num+1)*nw)];
    
#ifdef HAVE_OPENMP
#pragma omp simd
#endif
    for (int64_t j=0 ; j<elec_num ; ++j) {
      for (int64_t i=0 ; i<elec_num ; ++i) {
        een_rescaled_e_[i + elec_num*j] = 1.0;
      }
    }
    for (int64_t j=0 ; j<elec_num ; ++j) {
      een_rescaled_e_[j + elec_num*j] = 0.0;
    }
    if (cord_num > 0) {
      const double kappa = -rescale_factor_ee;
      for (int64_t j=0 ; j<elec_num ; ++j) {
#ifdef HAVE_OPENMP
#pragma omp simd
#endif
        for (int64_t i=0 ; i<elec_num ; ++i) {
          een_rescaled_e_[i + elec_num*(j + elec_num)] =
            exp(kappa * ee_distance[i + elec_num*(j + elec_num*nw)]);
        }
      }
      for (int64_t j=0 ; j<elec_num ; ++j) {
        een_rescaled_e_[j + elec_num*(j + elec_num)] = 0.0;
      }
    }
    for (int64_t l=2 ; l<=cord_num ; ++l) {
      for (int64_t j=0 ; j<elec_num ; ++j) {
#ifdef HAVE_OPENMP
#pragma omp simd
#endif
        for (int64_t i=0 ; i<elec_num ; ++i) {
          een_rescaled_e_[i + elec_num*(j + elec_num*l)] =
            een_rescaled_e_[i + elec_num*(j + elec_num*(l-1))] *
            een_rescaled_e_[i + elec_num*(j + elec_num)];
        }
      }
      for (int64_t j=0 ; j<elec_num ; ++j) {
        een_rescaled_e_[j + elec_num*(j + elec_num*l)] = 0.0;
      }
    }
  } // OpenMP
  
  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_compute_een_rescaled_e (
	  const qmckl_context context,
	  const int64_t walk_num,
	  const int64_t elec_num,
	  const int64_t cord_num,
	  const double rescale_factor_ee,
	  const double* ee_distance,
	  double* const een_rescaled_e ) {

#ifdef HAVE_HPC
    return qmckl_compute_een_rescaled_e_hpc
#else
    return qmckl_compute_een_rescaled_e_doc
#endif
      (context, walk_num, elec_num, cord_num, rescale_factor_ee, ee_distance, een_rescaled_e);
    }

qmckl_exit_code
qmckl_get_jastrow_champ_een_rescaled_e_gl(qmckl_context context,
                                          double* const een_rescaled_e_gl,
                                          const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_een_rescaled_e_gl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (een_rescaled_e_gl == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_champ_een_rescaled_e_gl",
                           "Null pointer");
  }

  const int64_t sze = ctx->electron.num * 4 * ctx->electron.num *
    ctx->electron.walker.num * (ctx->jastrow_champ.cord_num + 1);

  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_een_rescaled_e_gl",
                           "Array too small. Expected elec_num*4*elec_num*walk_num*(cord_num + 1)");
  }

  memcpy(een_rescaled_e_gl, ctx->jastrow_champ.een_rescaled_e_gl, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_een_rescaled_e_gl(qmckl_context context)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  /* Check if ee distance is provided */
  qmckl_exit_code rc = qmckl_provide_ee_distance(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Check if een rescaled distance is provided */
  rc = qmckl_provide_een_rescaled_e(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow_champ.een_rescaled_e_gl_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      if (ctx->jastrow_champ.een_rescaled_e_gl != NULL) {
        rc = qmckl_free(context, ctx->jastrow_champ.een_rescaled_e_gl);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_een_rescaled_e_gl",
                                 "Unable to free ctx->jastrow_champ.een_rescaled_e_gl");
        }
        ctx->jastrow_champ.een_rescaled_e_gl = NULL;
      }
    }

    /* Allocate array */
    if (ctx->jastrow_champ.een_rescaled_e_gl == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.num * 4 * ctx->electron.num *
        ctx->electron.walker.num * (ctx->jastrow_champ.cord_num + 1) * sizeof(double);

      double* een_rescaled_e_gl = (double*) qmckl_malloc(context, mem_info);

      if (een_rescaled_e_gl == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_een_rescaled_e_gl",
                               NULL);
      }
      ctx->jastrow_champ.een_rescaled_e_gl = een_rescaled_e_gl;
    }

    rc = qmckl_compute_jastrow_champ_factor_een_rescaled_e_gl(context,
                                                     ctx->electron.walker.num,
                                                     ctx->electron.num,
                                                     ctx->jastrow_champ.cord_num,
                                                     ctx->jastrow_champ.rescale_factor_ee,
                                                     ctx->electron.walker.point.coord.data,
                                                     ctx->electron.ee_distance,
                                                     ctx->jastrow_champ.een_rescaled_e,
                                                     ctx->jastrow_champ.een_rescaled_e_gl);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow_champ.een_rescaled_e_gl_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_een_rescaled_e_gl_hpc (const qmckl_context context,
                                                          const int64_t walk_num,
                                                          const int64_t elec_num,
                                                          const int64_t cord_num,
                                                          const double rescale_factor_ee,
                                                          const double* restrict coord_ee,
                                                          const double* restrict ee_distance,
                                                          const double* restrict een_rescaled_e,
                                                          double* restrict const een_rescaled_e_gl )
{
  if (context == QMCKL_NULL_CONTEXT) return QMCKL_INVALID_CONTEXT;
  if (walk_num <= 0) return QMCKL_INVALID_ARG_2;
  if (elec_num <= 0) return QMCKL_INVALID_ARG_3;
  if (cord_num < 0) return QMCKL_INVALID_ARG_4;

#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
  {

    double* restrict elec_dist_gl0 = (double*) calloc(elec_num * elec_num, sizeof(double));
    double* restrict elec_dist_gl1 = (double*) calloc(elec_num * elec_num, sizeof(double));
    double* restrict elec_dist_gl2 = (double*) calloc(elec_num * elec_num, sizeof(double));
    double* restrict elec_dist_gl3 = (double*) calloc(elec_num * elec_num, sizeof(double));
    assert (elec_dist_gl0 != NULL);
    assert (elec_dist_gl1 != NULL);
    assert (elec_dist_gl2 != NULL);
    assert (elec_dist_gl3 != NULL);

#ifdef HAVE_OPENMP
#pragma omp for
#endif
    for (int64_t nw = 0; nw < walk_num; ++nw) {

      const double* restrict ee_ = &(ee_distance[nw*elec_num*elec_num]);
      const double* restrict eex_ = &(coord_ee[elec_num*nw]);
      const double* restrict eey_ = &(coord_ee[elec_num*(nw+walk_num)]);
      const double* restrict eez_ = &(coord_ee[elec_num*(nw+2*walk_num)]);
      for (int64_t j=0; j<elec_num; ++j) {

        double rij_inv[elec_num];
        const double* restrict ee = &(ee_[j*elec_num]);
        double* restrict elec_dist_gl0_ = &(elec_dist_gl0[j*elec_num]);
        double* restrict elec_dist_gl1_ = &(elec_dist_gl1[j*elec_num]);
        double* restrict elec_dist_gl2_ = &(elec_dist_gl2[j*elec_num]);
        double* restrict elec_dist_gl3_ = &(elec_dist_gl3[j*elec_num]);

#ifdef HAVE_OPENMP
#pragma omp simd
#endif
        for (int64_t i = 0; i < elec_num ; ++i) {
          rij_inv[i] = ee[i];
        }

        rij_inv[j] = 1.0;
#ifdef HAVE_OPENMP
#pragma omp simd
#endif
        for (int64_t i = 0; i < elec_num ; ++i) {
          rij_inv[i] = 1.0/rij_inv[i];
        }
        rij_inv[j] = 0.0;

        const double xj = eex_[j];
        const double yj = eey_[j];
        const double zj = eez_[j];

#ifdef HAVE_OPENMP
#pragma omp simd
#endif
        for (int64_t i = 0; i < elec_num ; ++i) {

          const double xi = eex_[i];
          const double yi = eey_[i];
          const double zi = eez_[i];

          elec_dist_gl0_[i] = rij_inv[i] * (xi-xj);
          elec_dist_gl1_[i] = rij_inv[i] * (yi-yj);
          elec_dist_gl2_[i] = rij_inv[i] * (zi-zj);
          elec_dist_gl3_[i] = rij_inv[i] + rij_inv[i];
        }
      }

      double* restrict eegl = &een_rescaled_e_gl[elec_num*4*elec_num*(cord_num+1)*nw];
      memset(eegl, 0, 4*elec_num*elec_num*sizeof(double));

      for (int64_t l=1; l<=cord_num; ++l) {

        const double kappa_l  = -rescale_factor_ee * (double) l;

        for (int64_t j=0; j<elec_num; ++j) {

          double* restrict eegl =
            &een_rescaled_e_gl[elec_num*4*(j+elec_num*(l+(cord_num+1)*nw))];

          const double* restrict ee =
            &een_rescaled_e[elec_num*(j+elec_num*(l+(cord_num+1)*nw))];

          double* restrict elec_dist_gl0_ = &(elec_dist_gl0[j*elec_num]);
          double* restrict elec_dist_gl1_ = &(elec_dist_gl1[j*elec_num]);
          double* restrict elec_dist_gl2_ = &(elec_dist_gl2[j*elec_num]);
          double* restrict elec_dist_gl3_ = &(elec_dist_gl3[j*elec_num]);

          double kee[elec_num];
#ifdef HAVE_OPENMP
#pragma omp simd
#endif
          for (int64_t i=0; i<elec_num; ++i) {
            kee[i] = kappa_l * ee[i];
          }
          
#ifdef HAVE_OPENMP
#pragma omp simd
#endif
          for (int64_t i=0; i<elec_num; ++i) {
            eegl[i           ] = kee[i] *  elec_dist_gl0_[i];
            eegl[i+elec_num  ] = kee[i] *  elec_dist_gl1_[i];
            eegl[i+elec_num*2] = kee[i] *  elec_dist_gl2_[i];
            eegl[i+elec_num*3] = kee[i] * (elec_dist_gl3_[i] + kappa_l);
          }
          eegl[j           ] = 0.0;
          eegl[j+elec_num  ] = 0.0;
          eegl[j+elec_num*2] = 0.0;
          eegl[j+elec_num*3] = 0.0;
        }
      }
    }

    free(elec_dist_gl0);
    free(elec_dist_gl1);
    free(elec_dist_gl2);
    free(elec_dist_gl3);
  }
  
  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_compute_jastrow_champ_factor_een_rescaled_e_gl (
         const qmckl_context context,
         const int64_t walk_num,
         const int64_t elec_num,
         const int64_t cord_num,
         const double rescale_factor_ee,
         const double* coord_ee,
         const double* ee_distance,
         const double* een_rescaled_e,
         double* const een_rescaled_e_gl )
{
#ifdef HAVE_HPC
  return qmckl_compute_jastrow_champ_factor_een_rescaled_e_gl_hpc
#else
  return qmckl_compute_jastrow_champ_factor_een_rescaled_e_gl_doc
#endif
    (context, walk_num, elec_num, cord_num, rescale_factor_ee,
    coord_ee, ee_distance, een_rescaled_e, een_rescaled_e_gl );
}

qmckl_exit_code
qmckl_get_jastrow_champ_een_rescaled_n(qmckl_context context,
                                 double* const een_rescaled_n,
                                 const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_een_rescaled_n(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (een_rescaled_n == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_champ_een_rescaled_n",
                           "Null pointer");
  }

  const int64_t sze = ctx->electron.num * ctx->nucleus.num * ctx->electron.walker.num * (ctx->jastrow_champ.cord_num + 1);

  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_een_rescaled_n",
                           "Array too small. Expected elec_num*nucl_num*walk_num*(cord_num + 1)");
  }

  memcpy(een_rescaled_n, ctx->jastrow_champ.een_rescaled_n, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_een_rescaled_n(qmckl_context context)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  /* Check if ee distance is provided */
  qmckl_exit_code rc = qmckl_provide_en_distance(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow_champ.een_rescaled_n_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      if (ctx->jastrow_champ.een_rescaled_n != NULL) {
        rc = qmckl_free(context, ctx->jastrow_champ.een_rescaled_n);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_een_rescaled_n",
                                 "Unable to free ctx->jastrow_champ.een_rescaled_n");
        }
        ctx->jastrow_champ.een_rescaled_n = NULL;
      }
    }

    /* Allocate array */
    if (ctx->jastrow_champ.een_rescaled_n == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.num * ctx->nucleus.num *
        ctx->electron.walker.num * (ctx->jastrow_champ.cord_num + 1) * sizeof(double);
      double* een_rescaled_n = (double*) qmckl_malloc(context, mem_info);

      if (een_rescaled_n == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_een_rescaled_n",
                               NULL);
      }
      ctx->jastrow_champ.een_rescaled_n = een_rescaled_n;
    }

    rc = qmckl_compute_een_rescaled_n(context,
                                      ctx->electron.walker.num,
                                      ctx->electron.num,
                                      ctx->nucleus.num,
                                      ctx->jastrow_champ.type_nucl_num,
                                      ctx->jastrow_champ.type_nucl_vector,
                                      ctx->jastrow_champ.cord_num,
                                      ctx->jastrow_champ.rescale_factor_en,
                                      ctx->electron.en_distance,
                                      ctx->jastrow_champ.een_rescaled_n);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow_champ.een_rescaled_n_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

/*
qmckl_exit_code qmckl_compute_een_rescaled_n_hpc (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t type_nucl_num,
      int64_t* const type_nucl_vector,
      const int64_t cord_num,
      const double* rescale_factor_en,
      const double* en_distance,
      double* const een_rescaled_n ) {


  if (context == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  }

  if (walk_num <= 0) {
    return QMCKL_INVALID_ARG_2;
  }

  if (elec_num <= 0) {
    return QMCKL_INVALID_ARG_3;
  }

  if (nucl_num <= 0) {
    return QMCKL_INVALID_ARG_4;
  }

  if (cord_num < 0) {
    return QMCKL_INVALID_ARG_5;
  }

  // Prepare table of exponentiated distances raised to appropriate power
  for (int i = 0; i < (walk_num*(cord_num+1)*nucl_num*elec_num); ++i) {
    een_rescaled_n[i] = 1.0;
  }

  for (int nw = 0; nw < walk_num; ++nw) {
    for (int a = 0; a < nucl_num; ++a) {
      for (int i = 0; i < elec_num; ++i) {
        een_rescaled_n[i + a*elec_num + nw * elec_num*nucl_num*(cord_num+1)] = 1.0;
        een_rescaled_n[i + a*elec_num + elec_num*nucl_num + nw*elec_num*nucl_num*(cord_num+1)] =
          exp(-rescale_factor_en[type_nucl_vector[a]] * en_distance[a + i*nucl_num + nw*elec_num*nucl_num]);
      }
    }

    for (int l = 2; l < (cord_num+1); ++l){
      for (int a = 0; a < nucl_num; ++a) {
        for (int i = 0; i < elec_num; ++i) {
          een_rescaled_n[i + a*elec_num + l*elec_num*nucl_num + nw*elec_num*nucl_num*(cord_num+1)] =
            een_rescaled_n[i + a*elec_num + (l-1)*elec_num*nucl_num + nw*elec_num*nucl_num*(cord_num+1)] *
            een_rescaled_n[i + a*elec_num +       elec_num*nucl_num + nw*elec_num*nucl_num*(cord_num+1)];
        }
      }
    }

  }

  return QMCKL_SUCCESS;
}
*/

qmckl_exit_code
qmckl_get_jastrow_champ_een_rescaled_n_gl(qmckl_context context,
                                         double* const een_rescaled_n_gl,
                                         const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_een_rescaled_n_gl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (een_rescaled_n_gl == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_champ_een_rescaled_n_gl",
                           "Null pointer");
  }

  const int64_t sze = ctx->electron.num * 4 * ctx->nucleus.num * ctx->electron.walker.num * (ctx->jastrow_champ.cord_num + 1);

  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_een_rescaled_n_gl",
                           "Array too small. Expected ctx->electron.num * 4 * ctx->nucleus.num * ctx->electron.walker.num * (ctx->jastrow_champ.cord_num + 1)");
  }

  memcpy(een_rescaled_n_gl, ctx->jastrow_champ.een_rescaled_n_gl, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_een_rescaled_n_gl(qmckl_context context)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  /* Check if ee distance is provided */
  qmckl_exit_code rc = qmckl_provide_en_distance(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Check if ee distance is provided */
  rc = qmckl_provide_een_rescaled_n(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow_champ.een_rescaled_n_gl_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      if (ctx->jastrow_champ.een_rescaled_n_gl != NULL) {
        rc = qmckl_free(context, ctx->jastrow_champ.een_rescaled_n_gl);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_een_rescaled_n_gl",
                                 "Unable to free ctx->jastrow_champ.een_rescaled_n_gl");
        }
        ctx->jastrow_champ.een_rescaled_n_gl = NULL;
      }
    }

    /* Allocate array */
    if (ctx->jastrow_champ.een_rescaled_n_gl == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.num * 4 * ctx->nucleus.num *
        ctx->electron.walker.num * (ctx->jastrow_champ.cord_num + 1) * sizeof(double);
      double* een_rescaled_n_gl = (double*) qmckl_malloc(context, mem_info);

      if (een_rescaled_n_gl == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_een_rescaled_n_gl",
                               NULL);
      }
      ctx->jastrow_champ.een_rescaled_n_gl = een_rescaled_n_gl;
    }

    rc = qmckl_compute_jastrow_champ_factor_een_rescaled_n_gl(context,
                                                     ctx->electron.walker.num,
                                                     ctx->electron.num,
                                                     ctx->nucleus.num,
                                                     ctx->jastrow_champ.type_nucl_num,
                                                     ctx->jastrow_champ.type_nucl_vector,
                                                     ctx->jastrow_champ.cord_num,
                                                     ctx->jastrow_champ.rescale_factor_en,
                                                     ctx->electron.walker.point.coord.data,
                                                     ctx->nucleus.coord.data,
                                                     ctx->electron.en_distance,
                                                     ctx->jastrow_champ.een_rescaled_n,
                                                     ctx->jastrow_champ.een_rescaled_n_gl);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow_champ.een_rescaled_n_gl_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_compute_dim_c_vector_hpc (
      const qmckl_context context,
      const int64_t cord_num,
      int64_t* const dim_c_vector){

  int  lmax;


  if (context == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  }

  if (cord_num < 0) {
    return QMCKL_INVALID_ARG_2;
  }

  *dim_c_vector = 0;

  for (int p=2; p <= cord_num; ++p){
    for (int k=p-1; k >= 0; --k) {
      if (k != 0) {
        lmax = p - k;
      } else {
        lmax = p - k - 2;
      }
      for (int l = lmax; l >= 0; --l) {
        if ( ((p - k - l) & 1)==1) continue;
        *dim_c_vector=*dim_c_vector+1;
      }
    }
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_compute_dim_c_vector (
      const qmckl_context context,
      const int64_t cord_num,
      int64_t* const dim_c_vector)
{
#ifdef QMCKL_HPC
  return qmckl_compute_dim_c_vector_hpc
#else
  return qmckl_compute_dim_c_vector_doc
#endif
    (context, cord_num, dim_c_vector);
}

qmckl_exit_code
qmckl_get_jastrow_champ_tmp_c(qmckl_context context,
                              double* const tmp_c,
                              const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_jastrow_champ_c_vector_full(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_tmp_c(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (tmp_c == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_champ_tmp_c",
                           "Null pointer");
  }

  const int64_t sze = (ctx->jastrow_champ.cord_num) * (ctx->jastrow_champ.cord_num + 1) *
    ctx->electron.num * ctx->nucleus.num * ctx->electron.walker.num;

  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_tmp_c",
                           "Array too small. Expected cord_num*(cord_num+1)*walk_num*elec_num*nucl_num");
  }

  memcpy(tmp_c, ctx->jastrow_champ.tmp_c, sze * sizeof(double));

  return QMCKL_SUCCESS;
}


qmckl_exit_code
qmckl_get_jastrow_champ_dtmp_c(qmckl_context context,
                               double* const dtmp_c,
                               const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_jastrow_champ_c_vector_full(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_dtmp_c(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (dtmp_c == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_champ_dtmp_c",
                           "Null pointer");
  }

  const int64_t sze = (ctx->jastrow_champ.cord_num) * (ctx->jastrow_champ.cord_num + 1)*
    4* ctx->electron.num * ctx->nucleus.num * ctx->electron.walker.num;

  if (size_max < sze) {
    return qmckl_failwith(context,
                          QMCKL_INVALID_ARG_3,
                          "qmckl_get_jastrow_champ_dtmp_c",
                          "Array too small. Expected 4*cord_num*(cord_num+1)*walk_num*elec_num*nucl_num");
  }

  memcpy(dtmp_c, ctx->jastrow_champ.dtmp_c, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_jastrow_champ_c_vector_full(qmckl_context context)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  qmckl_exit_code rc = QMCKL_SUCCESS;

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow_champ.c_vector_full_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      if (ctx->jastrow_champ.c_vector_full != NULL) {
        rc = qmckl_free(context, ctx->jastrow_champ.c_vector_full);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_jastrow_champ_c_vector_full",
                                 "Unable to free ctx->jastrow_champ.c_vector_full");
        }
        ctx->jastrow_champ.c_vector_full = NULL;
      }
    }

    /* Allocate array */
    if (ctx->jastrow_champ.c_vector_full == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->jastrow_champ.dim_c_vector * ctx->nucleus.num * sizeof(double);
      double* c_vector_full = (double*) qmckl_malloc(context, mem_info);

      if (c_vector_full == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_jastrow_champ_c_vector_full",
                               NULL);
      }
      ctx->jastrow_champ.c_vector_full = c_vector_full;
    }

    rc = qmckl_compute_c_vector_full(context,
                                      ctx->nucleus.num,
                                      ctx->jastrow_champ.dim_c_vector,
                                      ctx->jastrow_champ.type_nucl_num,
                                      ctx->jastrow_champ.type_nucl_vector,
                                      ctx->jastrow_champ.c_vector,
                                      ctx->jastrow_champ.c_vector_full);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow_champ.c_vector_full_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_lkpm_combined_index(qmckl_context context)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  qmckl_exit_code rc = QMCKL_SUCCESS;

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow_champ.lkpm_combined_index_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      if (ctx->jastrow_champ.lkpm_combined_index != NULL) {
        rc = qmckl_free(context, ctx->jastrow_champ.lkpm_combined_index);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_jastrow_champ_factor_ee",
                                 "Unable to free ctx->jastrow_champ.lkpm_combined_index");
        }
        ctx->jastrow_champ.lkpm_combined_index = NULL;
      }
    }

    /* Allocate array */
    if (ctx->jastrow_champ.lkpm_combined_index == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = 4 * ctx->jastrow_champ.dim_c_vector * sizeof(int64_t);
      int64_t* lkpm_combined_index = (int64_t*) qmckl_malloc(context, mem_info);

      if (lkpm_combined_index == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_lkpm_combined_index",
                               NULL);
      }
      ctx->jastrow_champ.lkpm_combined_index = lkpm_combined_index;
    }

    rc = qmckl_compute_lkpm_combined_index(context,
                                           ctx->jastrow_champ.cord_num,
                                           ctx->jastrow_champ.dim_c_vector,
                                           ctx->jastrow_champ.lkpm_combined_index);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow_champ.lkpm_combined_index_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_tmp_c(qmckl_context context)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  qmckl_exit_code rc = QMCKL_SUCCESS;

  rc = qmckl_provide_een_rescaled_e(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_een_rescaled_n(context);
  if (rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow_champ.tmp_c_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      if (ctx->jastrow_champ.tmp_c != NULL) {
        rc = qmckl_free(context, ctx->jastrow_champ.tmp_c);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_tmp_c",
                                 "Unable to free ctx->jastrow_champ.tmp_c");
        }
        ctx->jastrow_champ.tmp_c = NULL;
      }
    }

    /* Allocate array */
    if (ctx->jastrow_champ.tmp_c == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = (ctx->jastrow_champ.cord_num) * (ctx->jastrow_champ.cord_num + 1)
                      * ctx->electron.num * ctx->nucleus.num * ctx->electron.walker.num * sizeof(double);
      double* tmp_c = (double*) qmckl_malloc(context, mem_info);

      if (tmp_c == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_tmp_c",
                               NULL);
      }
      ctx->jastrow_champ.tmp_c = tmp_c;
    }

    rc = qmckl_compute_tmp_c(context,
                               ctx->jastrow_champ.cord_num,
                               ctx->electron.num,
                               ctx->nucleus.num,
                               ctx->electron.walker.num,
                               ctx->jastrow_champ.een_rescaled_e,
                               ctx->jastrow_champ.een_rescaled_n,
                               ctx->jastrow_champ.tmp_c);

    ctx->jastrow_champ.tmp_c_date = ctx->date;
  }

  return rc;
}

qmckl_exit_code qmckl_provide_dtmp_c(qmckl_context context)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  rc = qmckl_provide_een_rescaled_e_gl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_een_rescaled_n(context);
  if (rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow_champ.dtmp_c_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      if (ctx->jastrow_champ.dtmp_c != NULL) {
        rc = qmckl_free(context, ctx->jastrow_champ.dtmp_c);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_dtmp_c",
                                 "Unable to free ctx->jastrow_champ.dtmp_c");
        }
        ctx->jastrow_champ.dtmp_c = NULL;
      }
    }

    /* Allocate array */
    if (ctx->jastrow_champ.dtmp_c == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = (ctx->jastrow_champ.cord_num) * (ctx->jastrow_champ.cord_num + 1)
                      * 4 * ctx->electron.num * ctx->nucleus.num * ctx->electron.walker.num * sizeof(double);
      double* dtmp_c = (double*) qmckl_malloc(context, mem_info);

      if (dtmp_c == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_dtmp_c",
                               NULL);
      }
      ctx->jastrow_champ.dtmp_c = dtmp_c;
    }


    rc = qmckl_compute_dtmp_c(context,
                              ctx->jastrow_champ.cord_num,
                              ctx->electron.num,
                              ctx->nucleus.num,
                              ctx->electron.walker.num,
                              ctx->jastrow_champ.een_rescaled_e_gl,
                              ctx->jastrow_champ.een_rescaled_n,
                              ctx->jastrow_champ.dtmp_c);

    if (rc != QMCKL_SUCCESS) {
      return rc;
    }


    ctx->jastrow_champ.dtmp_c_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_compute_c_vector_full_hpc (
	  const qmckl_context context,
	  const int64_t nucl_num,
	  const int64_t dim_c_vector,
	  const int64_t type_nucl_num,
	  const int64_t* type_nucl_vector,
	  const double* c_vector,
	  double* const c_vector_full ) {

  if (context == QMCKL_NULL_CONTEXT) return QMCKL_INVALID_CONTEXT;
  if (nucl_num <= 0)                 return QMCKL_INVALID_ARG_2;
  if (dim_c_vector < 0)              return QMCKL_INVALID_ARG_3;
  if (type_nucl_num <= 0)            return QMCKL_INVALID_ARG_4;
  if (type_nucl_vector == NULL)      return QMCKL_INVALID_ARG_5;
  if (c_vector == NULL)              return QMCKL_INVALID_ARG_6;
  if (c_vector_full == NULL)         return QMCKL_INVALID_ARG_7;

  for (int i=0; i < dim_c_vector; ++i) {
    for (int a=0; a < nucl_num; ++a){
      c_vector_full[a + i*nucl_num] = c_vector[i + type_nucl_vector[a]*dim_c_vector];
    }
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_compute_c_vector_full (
	  const qmckl_context context,
	  const int64_t nucl_num,
	  const int64_t dim_c_vector,
	  const int64_t type_nucl_num,
	  const int64_t* type_nucl_vector,
	  const double* c_vector,
	  double* const c_vector_full ) {

#ifdef HAVE_HPC
      return qmckl_compute_c_vector_full_hpc
#else
      return qmckl_compute_c_vector_full_doc
#endif
(context, nucl_num, dim_c_vector, type_nucl_num, type_nucl_vector, c_vector, c_vector_full);
    }

qmckl_exit_code qmckl_compute_lkpm_combined_index_hpc (
      const qmckl_context context,
      const int64_t cord_num,
      const int64_t dim_c_vector,
      int64_t* const lkpm_combined_index ) {

  int kk, lmax, m;

  if (context == QMCKL_NULL_CONTEXT) return QMCKL_INVALID_CONTEXT;
  if (cord_num < 0) return QMCKL_INVALID_ARG_2;
  if (dim_c_vector < 0) return QMCKL_INVALID_ARG_3;

  kk = 0;
  for (int p = 2; p <= cord_num; ++p) {
    for (int k=(p-1); k >= 0; --k) {
      if (k != 0) {
        lmax = p - k;
      } else {
        lmax = p - k - 2;
      }
      for (int l=lmax; l >= 0; --l) {
        if (((p - k - l) & 1) == 1) continue;
        m = (p - k - l)/2;
        lkpm_combined_index[kk                 ] = l;
        lkpm_combined_index[kk +   dim_c_vector] = k;
        lkpm_combined_index[kk + 2*dim_c_vector] = p;
        lkpm_combined_index[kk + 3*dim_c_vector] = m;
        kk = kk + 1;
      }
    }
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_compute_lkpm_combined_index (
      const qmckl_context context,
      const int64_t cord_num,
      const int64_t dim_c_vector,
      int64_t* const lkpm_combined_index ) {

    #ifdef HAVE_HPC
      return qmckl_compute_lkpm_combined_index_hpc
    #else
      return qmckl_compute_lkpm_combined_index_doc
    #endif
        (context, cord_num, dim_c_vector, lkpm_combined_index);
}

/* Compute tmp_c */
/*      :PROPERTIES: */
/*      :Name:     qmckl_compute_tmp_c */
/*      :CRetType: qmckl_exit_code */
/*      :FRetType: qmckl_exit_code */
/*      :END: */

/*      #+NAME: qmckl_factor_tmp_c_args */
/*      |------------------+------------------------------------------------------------------+--------+-----------------------------------| */
/*      | Variable         | Type                                                             | In/Out | Description                       | */
/*      |------------------+------------------------------------------------------------------+--------+-----------------------------------| */
/*      | ~context~        | ~qmckl_context~                                                  | in     | Global state                      | */
/*      | ~cord_num~       | ~int64_t~                                                        | in     | Order of polynomials              | */
/*      | ~elec_num~       | ~int64_t~                                                        | in     | Number of electrons               | */
/*      | ~nucl_num~       | ~int64_t~                                                        | in     | Number of nuclei                 | */
/*      | ~walk_num~       | ~int64_t~                                                        | in     | Number of walkers                 | */
/*      | ~een_rescaled_e~ | ~double[walk_num][0:cord_num][elec_num][elec_num]~               | in     | Electron-electron rescaled factor | */
/*      | ~een_rescaled_n~ | ~double[walk_num][0:cord_num][nucl_num][elec_num]~               | in     | Electron-nucleus rescaled factor  | */
/*      | ~tmp_c~          | ~double[walk_num][0:cord_num-1][0:cord_num][nucl_num][elec_num]~ | out    | vector of non-zero coefficients   | */
/*      |------------------+------------------------------------------------------------------+--------+-----------------------------------| */


qmckl_exit_code qmckl_compute_tmp_c (const qmckl_context context,
                                     const int64_t cord_num,
                                     const int64_t elec_num,
                                     const int64_t nucl_num,
                                     const int64_t walk_num,
                                     const double* een_rescaled_e,
                                     const double* een_rescaled_n,
                                     double* const tmp_c )
{
#ifdef HAVE_HPC
  return qmckl_compute_tmp_c_hpc
#else
  return qmckl_compute_tmp_c_doc
#endif
    (context, cord_num, elec_num, nucl_num, walk_num,
     een_rescaled_e, een_rescaled_n, tmp_c);
}

/* HPC                                                          :noexport: */


qmckl_exit_code qmckl_compute_tmp_c_hpc (
      const qmckl_context context,
      const int64_t cord_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t walk_num,
      const double* een_rescaled_e,
      const double* een_rescaled_n,
      double* const tmp_c ) {

  if (context == QMCKL_NULL_CONTEXT) return QMCKL_INVALID_CONTEXT;
  if (cord_num <  0)                 return QMCKL_INVALID_ARG_2;
  if (elec_num <= 0)                 return QMCKL_INVALID_ARG_3;
  if (nucl_num <= 0)                 return QMCKL_INVALID_ARG_4;
  if (walk_num <= 0)                 return QMCKL_INVALID_ARG_5;

  const double alpha = 1.0;
  const double beta  = 0.0;

  const int64_t M = elec_num;
  const int64_t K = elec_num;
  const int64_t N = nucl_num;

  const int64_t LDA = M;
  const int64_t LDB = K;
  const int64_t LDC = M;

  const int64_t af = LDA*elec_num;
  const int64_t bf = LDB*nucl_num;
  const int64_t cf = LDC*nucl_num;
  
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
  for (int64_t nw=0; nw < walk_num; ++nw) {
    // (i,0)
    { int64_t j=0 ;
      for (int64_t i=1; i<cord_num; ++i){
        double* rowsum = calloc(M, sizeof(double));
        for (int64_t k=0 ; k<K ; ++k) {
          for (int64_t m=0 ; m<M ; ++m) {
            rowsum[m] += een_rescaled_e[m + k*LDA + af*(i+nw*(cord_num+1))];
          }
        }
        for (int64_t n=0 ; n<N ; ++n) {
          for (int64_t m=0 ; m<M ; ++m) {
            tmp_c[m + n*LDC + cf*(j+(cord_num+1)*(i+nw*cord_num))] = rowsum[m];
          }
        }
        free(rowsum);
      }
    }

    // (0,j)
    { int64_t i=0;
      for (int64_t j=0; j<=cord_num; ++j){
        for (int64_t l=0 ; l<N ; ++l) {
          double colsum = 0.0;
          for (int64_t k=0 ; k<K ; ++k) {
            colsum += een_rescaled_n[k + l*elec_num + bf*(j+nw*(cord_num+1))];
          }
          for (int64_t k=0 ; k<M ; ++k) {
            tmp_c[k + l*elec_num + cf*(j+(cord_num+1)*(i+nw*cord_num))] = colsum - 
              een_rescaled_n[k + l*elec_num + bf*(j+nw*(cord_num+1))];
          }
        }
      }
    }
  }

 double eps = qmckl_get_numprec_epsilon(context);

#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
  for (int64_t nw=0; nw < walk_num; ++nw) {

    CSCMatrix *eN[cord_num];
    for (int64_t j=1; j<cord_num; ++j){
      eN[j] = qmckl_dense_to_csc(K, N, &(een_rescaled_n[bf*(j+nw*(cord_num+1))]),
                                 LDB, eps);
    }
    
    for (int64_t i=1; i<cord_num; ++i){
      for (int64_t j=1; j<=cord_num-i; ++j){
        qmckl_dgemm_sparse_b_nn(context, M, N, K, alpha,
                                &(een_rescaled_e[af*(i+nw*(cord_num+1))]), LDA,
                                eN[j], beta,
                                &(tmp_c[cf*(j+(cord_num+1)*(i+nw*cord_num))]), LDC);
      }
    }
    for (int64_t j=1; j<cord_num; ++j){
      free_csc(eN[j]);
    }
  }
    
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_compute_dtmp_c (const qmckl_context context,
                      const int64_t cord_num,
                      const int64_t elec_num,
                      const int64_t nucl_num,
                      const int64_t walk_num,
                      const double* een_rescaled_e_gl,
                      const double* een_rescaled_n,
                      double* const dtmp_c )
{
#ifdef HAVE_HPC
  return qmckl_compute_dtmp_c_hpc
#else
  return qmckl_compute_dtmp_c_doc
#endif
    (context, cord_num, elec_num, nucl_num, walk_num, een_rescaled_e_gl,
     een_rescaled_n, dtmp_c );
}

/* CPU                                                          :noexport: */


qmckl_exit_code
qmckl_compute_dtmp_c_hpc (const qmckl_context context,
                          const int64_t cord_num,
                          const int64_t elec_num,
                          const int64_t nucl_num,
                          const int64_t walk_num,
                          const double* een_rescaled_e_gl,
                          const double* een_rescaled_n,
                          double* const dtmp_c )
{

  if (context == QMCKL_NULL_CONTEXT) return QMCKL_INVALID_CONTEXT;
  if (cord_num <  0)                 return QMCKL_INVALID_ARG_2;
  if (elec_num <= 0)                 return QMCKL_INVALID_ARG_3;
  if (nucl_num <= 0)                 return QMCKL_INVALID_ARG_4;
  if (walk_num <= 0)                 return QMCKL_INVALID_ARG_5;

  qmckl_exit_code  info = QMCKL_SUCCESS;

  const double alpha = 1.0;
  const double beta  = 0.0;

  const int64_t M = 4*elec_num;
  const int64_t N = nucl_num;
  const int64_t K = elec_num;

  const int64_t LDA = M;
  const int64_t LDB = K;
  const int64_t LDC = M;

  const int64_t af = LDA*elec_num;
  const int64_t bf = LDB*nucl_num;
  const int64_t cf = LDC*nucl_num;

#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
  for (int64_t nw=0; nw < walk_num; ++nw) {
    // (i,0)
    { int64_t j=0 ;
      for (int64_t i=1; i<cord_num; ++i){
        double rowsum[M];
        for (int64_t m=0 ; m<M ; ++m) {
          rowsum[m] = 0.0;
        }
        for (int64_t k=0 ; k<K ; ++k) {
          for (int64_t m=0 ; m<M ; ++m) {
            rowsum[m] += een_rescaled_e_gl[m + k*LDA + af*(i+nw*(cord_num+1))];
          }
        }
        for (int64_t n=0 ; n<N ; ++n) {
          for (int64_t m=0 ; m<M ; ++m) {
            dtmp_c[m + n*LDC + cf*(j+(cord_num+1)*(i+nw*cord_num))] = rowsum[m];
          }
        }
      }
    }

    // (0,j)
    { int64_t i=0;
      for (int64_t j=0; j<=cord_num; ++j){
        for (int64_t k=0 ; k<N*M ; ++k) {
          dtmp_c[k + cf*(j+(cord_num+1)*(i+nw*cord_num))] = 0.0;
        }
      }
    }
  }

 double eps = qmckl_get_numprec_epsilon(context);

#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
  for (int64_t nw=0; nw < walk_num; ++nw) {

    CSCMatrix *eN[cord_num];
    for (int64_t j=1; j<cord_num; ++j){
      eN[j] = qmckl_dense_to_csc(K, N, &(een_rescaled_n[bf*(j+nw*(cord_num+1))]),
                                 LDB, eps);
    }
    
    for (int64_t i=1; i<cord_num; ++i){
      for (int64_t j=1; j<=cord_num-i; ++j){
        qmckl_dgemm_sparse_b_nn(context, M, N, K, alpha,
                                &(een_rescaled_e_gl[af*(i+nw*(cord_num+1))]), LDA,
                                eN[j], beta,
                                &(dtmp_c[cf*(j+(cord_num+1)*(i+nw*cord_num))]), LDC);
      }
    }
    for (int64_t j=1; j<cord_num; ++j){
      free_csc(eN[j]);
    }
  }

  return info;
}

qmckl_exit_code
qmckl_get_jastrow_champ_factor_een(qmckl_context context,
                             double* const factor_een,
                             const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_jastrow_champ_factor_een(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (factor_een == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_champ_factor_een",
                           "Null pointer");
  }

  const int64_t sze = ctx->electron.walker.num;

  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_factor_een",
                           "Array too small. Expected walk_num");
  }

  memcpy(factor_een, ctx->jastrow_champ.factor_een, sze*sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_jastrow_champ_factor_een(qmckl_context context)
{

  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (ctx->jastrow_champ.cord_num > 0) {

    /* Check if en rescaled distance is provided */
    rc = qmckl_provide_een_rescaled_e(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if en rescaled distance derivatives is provided */
    rc = qmckl_provide_een_rescaled_n(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if en rescaled distance derivatives is provided */
    rc = qmckl_provide_jastrow_champ_c_vector_full(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if en rescaled distance derivatives is provided */
    rc = qmckl_provide_lkpm_combined_index(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if tmp_c is provided */
    rc = qmckl_provide_tmp_c(context);
    if(rc != QMCKL_SUCCESS) return rc;

  }

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow_champ.factor_een_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      if (ctx->jastrow_champ.factor_een != NULL) {
        rc = qmckl_free(context, ctx->jastrow_champ.factor_een);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_jastrow_champ_factor_een",
                                 "Unable to free ctx->jastrow_champ.factor_een");
        }
        ctx->jastrow_champ.factor_een = NULL;
      }
    }

    /* Allocate array */
    if (ctx->jastrow_champ.factor_een == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.walker.num * sizeof(double);
      double* factor_een = (double*) qmckl_malloc(context, mem_info);

      if (factor_een == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_jastrow_champ_factor_een",
                               NULL);
      }
      ctx->jastrow_champ.factor_een = factor_een;
    }

    rc = qmckl_compute_jastrow_champ_factor_een(context,
                                                ctx->electron.walker.num,
                                                ctx->electron.num,
                                                ctx->nucleus.num,
                                                ctx->jastrow_champ.cord_num,
                                                ctx->jastrow_champ.dim_c_vector,
                                                ctx->jastrow_champ.c_vector_full,
                                                ctx->jastrow_champ.lkpm_combined_index,
                                                ctx->jastrow_champ.tmp_c,
                                                ctx->jastrow_champ.een_rescaled_n,
                                                ctx->jastrow_champ.factor_een);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow_champ.factor_een_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_een_hpc(const qmckl_context context,
                                           const int64_t walk_num,
                                           const int64_t elec_num,
                                           const int64_t nucl_num,
                                           const int64_t cord_num,
                                           const int64_t dim_c_vector,
                                           const double* restrict c_vector_full,
                                           const int64_t* restrict lkpm_combined_index,
                                           const double* restrict tmp_c,
                                           const double* restrict een_rescaled_n,
                                           double* restrict const factor_een)
{

  int64_t info = QMCKL_SUCCESS;

  if (context == QMCKL_NULL_CONTEXT) return QMCKL_INVALID_CONTEXT;
  if (walk_num <= 0)                 return QMCKL_INVALID_ARG_2;
  if (elec_num <= 0)                 return QMCKL_INVALID_ARG_3;
  if (nucl_num <= 0)                 return QMCKL_INVALID_ARG_4;
  if (cord_num <  0)                 return QMCKL_INVALID_ARG_5;


  memset(factor_een, 0, walk_num*sizeof(double));

  if (cord_num == 0) {
    return QMCKL_SUCCESS;
  }

#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
  for (size_t nw = 0; nw < (size_t) walk_num; ++nw) {
    for (size_t n = 0; n < (size_t) dim_c_vector; ++n) {
      const size_t l = lkpm_combined_index[n];
      const size_t k = lkpm_combined_index[n+  dim_c_vector];
      const size_t m = lkpm_combined_index[n+3*dim_c_vector];

      const size_t en = elec_num*nucl_num;
      const size_t len = l*en;
      const size_t c1 = cord_num+1;
      const size_t cn = cord_num*nw;
      const size_t addr0 = en*(m+c1*(k+cn));
      const size_t addr1 = en*(m+c1*nw);

      const double* restrict tmp_c_mlkn = &(tmp_c[addr0]) + len;
      const double* restrict een_rescaled_n_mnw = &(een_rescaled_n[addr1]);

      const double* restrict cv = &(c_vector_full[n*nucl_num]);
      for (size_t a = 0; a < (size_t) nucl_num; a++) {
        double cn = cv[a];
        if (cn == 0.0) continue;

        const size_t ishift  = elec_num*a;

        const double* restrict tmp_c_amlkn = tmp_c_mlkn + ishift;
        const double* restrict een_rescaled_n_amnw  = een_rescaled_n_mnw + ishift;

        double accu = 0.;
#ifdef HAVE_OPENMP
#pragma omp simd
#endif
        for (size_t j = 0; j < (size_t) elec_num; ++j) {
          accu = accu + tmp_c_amlkn[j] * een_rescaled_n_amnw [j];
        }
        factor_een[nw] += cn * accu;

      }
    }
  }
  return info;
}

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_een (const qmckl_context context,
                          const int64_t walk_num,
                          const int64_t elec_num,
                          const int64_t nucl_num,
                          const int64_t cord_num,
                          const int64_t dim_c_vector,
                          const double* c_vector_full,
                          const int64_t* lkpm_combined_index,
                          const double* tmp_c,
                          const double* een_rescaled_n,
                          double* const factor_een )
{
#ifdef HAVE_HPC
  return qmckl_compute_jastrow_champ_factor_een_hpc
#else
  return qmckl_compute_jastrow_champ_factor_een_doc
#endif
    (context, walk_num, elec_num, nucl_num, cord_num, dim_c_vector,
     c_vector_full, lkpm_combined_index, tmp_c, een_rescaled_n,
     factor_een );
}

qmckl_exit_code
qmckl_get_jastrow_champ_factor_een_gl(qmckl_context context,
                                      double* const factor_een_gl,
                                      const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_jastrow_champ_factor_een_gl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (factor_een_gl == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_champ_factor_een_gl",
                           "Null pointer");
  }

  const int64_t sze = ctx->electron.walker.num * 4 * ctx->electron.num;

  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_factor_een_gl",
                           "Array too small. Expected 4*walk_num*elec_num");
  }

  memcpy(factor_een_gl, ctx->jastrow_champ.factor_een_gl, sze*sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_jastrow_champ_factor_een_grad(qmckl_context context,
                                        double* const factor_een_grad,
                                        const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_jastrow_champ_factor_een_grad(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t sze = ctx->electron.walker.num * 3 * ctx->electron.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_factor_een_grad",
                           "Array too small. Expected 3*walk_num*elec_num");
  }
  memcpy(factor_een_grad, ctx->jastrow_champ.factor_een_grad, sze*sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_jastrow_champ_factor_een_gl(qmckl_context context)
{

  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (ctx->jastrow_champ.cord_num > 0) {

    /* Check if en rescaled distance is provided */
    rc = qmckl_provide_een_rescaled_e(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if en rescaled distance derivatives is provided */
    rc = qmckl_provide_een_rescaled_n(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if en rescaled distance is provided */
    rc = qmckl_provide_een_rescaled_e_gl(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if en rescaled distance derivatives is provided */
    rc = qmckl_provide_een_rescaled_n_gl(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if en rescaled distance derivatives is provided */
    rc = qmckl_provide_jastrow_champ_c_vector_full(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if en rescaled distance derivatives is provided */
    rc = qmckl_provide_lkpm_combined_index(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if tmp_c is provided */
    rc = qmckl_provide_tmp_c(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if dtmp_c is provided */
    rc = qmckl_provide_dtmp_c(context);
    if(rc != QMCKL_SUCCESS) return rc;

  }

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow_champ.factor_een_gl_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      if (ctx->jastrow_champ.factor_een_gl != NULL) {
        rc = qmckl_free(context, ctx->jastrow_champ.factor_een_gl);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_jastrow_champ_factor_een_gl",
                                 "Unable to free ctx->jastrow_champ.factor_een_gl");
        }
        ctx->jastrow_champ.factor_een_gl = NULL;
      }
    }

    /* Allocate array */
    if (ctx->jastrow_champ.factor_een_gl == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = 4 * ctx->electron.num * ctx->electron.walker.num * sizeof(double);
      double* factor_een_gl = (double*) qmckl_malloc(context, mem_info);

      if (factor_een_gl == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_jastrow_champ_factor_een_gl",
                               NULL);
      }
      ctx->jastrow_champ.factor_een_gl = factor_een_gl;
    }

    rc = qmckl_compute_jastrow_champ_factor_een_gl(context,
                                                   ctx->electron.walker.num,
                                                   ctx->electron.num,
                                                   ctx->nucleus.num,
                                                   ctx->jastrow_champ.cord_num,
                                                   ctx->jastrow_champ.dim_c_vector,
                                                   ctx->jastrow_champ.c_vector_full,
                                                   ctx->jastrow_champ.lkpm_combined_index,
                                                   ctx->jastrow_champ.tmp_c,
                                                   ctx->jastrow_champ.dtmp_c,
                                                   ctx->jastrow_champ.een_rescaled_n,
                                                   ctx->jastrow_champ.een_rescaled_n_gl,
                                                   ctx->jastrow_champ.factor_een_gl);

    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow_champ.factor_een_gl_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}


qmckl_exit_code qmckl_provide_jastrow_champ_factor_een_grad(qmckl_context context)
{

  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (ctx->jastrow_champ.cord_num > 0) {

    /* Check if en rescaled distance is provided */
    rc = qmckl_provide_een_rescaled_e(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if en rescaled distance derivatives is provided */
    rc = qmckl_provide_een_rescaled_n(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if en rescaled distance is provided */
    rc = qmckl_provide_een_rescaled_e_gl(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if en rescaled distance derivatives is provided */
    rc = qmckl_provide_een_rescaled_n_gl(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if en rescaled distance derivatives is provided */
    rc = qmckl_provide_jastrow_champ_c_vector_full(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if en rescaled distance derivatives is provided */
    rc = qmckl_provide_lkpm_combined_index(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if tmp_c is provided */
    rc = qmckl_provide_tmp_c(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if dtmp_c is provided */
    rc = qmckl_provide_dtmp_c(context);
    if(rc != QMCKL_SUCCESS) return rc;

  }

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow_champ.factor_een_grad_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      if (ctx->jastrow_champ.factor_een_grad != NULL) {
        rc = qmckl_free(context, ctx->jastrow_champ.factor_een_grad);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_jastrow_champ_factor_een_grad",
                                 "Unable to free ctx->jastrow_champ.factor_een_grad");
        }
        ctx->jastrow_champ.factor_een_grad = NULL;
      }
    }

    /* Allocate array */
    if (ctx->jastrow_champ.factor_een_grad == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = 4 * ctx->electron.num * ctx->electron.walker.num * sizeof(double);
      double* factor_een_grad = (double*) qmckl_malloc(context, mem_info);

      if (factor_een_grad == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_jastrow_champ_factor_een_grad",
                               NULL);
      }
      ctx->jastrow_champ.factor_een_grad = factor_een_grad;
    }

    rc = qmckl_compute_jastrow_champ_factor_een_grad(context,
                                                     ctx->electron.walker.num,
                                                     ctx->electron.num,
                                                     ctx->nucleus.num,
                                                     ctx->jastrow_champ.cord_num,
                                                     ctx->jastrow_champ.dim_c_vector,
                                                     ctx->jastrow_champ.c_vector_full,
                                                     ctx->jastrow_champ.lkpm_combined_index,
                                                     ctx->jastrow_champ.tmp_c,
                                                     ctx->jastrow_champ.dtmp_c,
                                                     ctx->jastrow_champ.een_rescaled_n,
                                                     ctx->jastrow_champ.een_rescaled_n_gl,
                                                     ctx->jastrow_champ.factor_een_grad);

    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow_champ.factor_een_grad_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_een_gl(const qmckl_context context,
                                               const int64_t walk_num,
                                               const int64_t elec_num,
                                               const int64_t nucl_num,
                                               const int64_t cord_num,
                                               const int64_t dim_c_vector,
                                               const double *c_vector_full,
                                               const int64_t *lkpm_combined_index,
                                               const double *tmp_c,
                                               const double *dtmp_c,
                                               const double *een_rescaled_n,
                                               const double *een_rescaled_n_gl,
                                               double* const factor_een_gl)
{
#ifdef HAVE_HPC
  return qmckl_compute_jastrow_champ_factor_een_gl_hpc
#else
  return qmckl_compute_jastrow_champ_factor_een_gl_doc
#endif
    (context, walk_num, elec_num, nucl_num,
     cord_num, dim_c_vector, c_vector_full,
     lkpm_combined_index, tmp_c, dtmp_c,
     een_rescaled_n, een_rescaled_n_gl,
     factor_een_gl);
}

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_een_gl_hpc(const qmckl_context context,
                                                   const int64_t walk_num,
                                                   const int64_t elec_num,
                                                   const int64_t nucl_num,
                                                   const int64_t cord_num,
                                                   const int64_t dim_c_vector,
                                                   const double* restrict c_vector_full,
                                                   const int64_t* restrict lkpm_combined_index,
                                                   const double* restrict tmp_c,
                                                   const double* restrict dtmp_c,
                                                   const double* restrict een_rescaled_n,
                                                   const double* restrict een_rescaled_n_gl,
                                                   double* restrict const factor_een_gl)
{

  int64_t info = QMCKL_SUCCESS;

  if (context == QMCKL_NULL_CONTEXT) return QMCKL_INVALID_CONTEXT;
  if (walk_num <= 0)                 return QMCKL_INVALID_ARG_2;
  if (elec_num <= 0)                 return QMCKL_INVALID_ARG_3;
  if (nucl_num <= 0)                 return QMCKL_INVALID_ARG_4;
  if (cord_num <  0)                 return QMCKL_INVALID_ARG_5;


  if (cord_num == 0) {
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for (size_t nw = 0; nw < (size_t) walk_num; ++nw) {
      double* const restrict factor_een_gl_0nw = &(factor_een_gl[elec_num*4*nw]);
      memset(factor_een_gl_0nw, 0, elec_num*4*sizeof(double));
    }
    return QMCKL_SUCCESS;
  }

  const size_t elec_num2 = elec_num + elec_num;
  const size_t elec_num3 = elec_num2 + elec_num;

#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
  for (size_t nw = 0; nw < (size_t) walk_num; ++nw) {
    double* const restrict factor_een_gl_0nw = &(factor_een_gl[elec_num*4*nw]);
    memset(factor_een_gl_0nw, 0, elec_num*4*sizeof(double));
    for (size_t n = 0; n < (size_t) dim_c_vector; ++n) {
      const size_t l = lkpm_combined_index[n];
      const size_t k = lkpm_combined_index[n+  dim_c_vector];
      const size_t m = lkpm_combined_index[n+3*dim_c_vector];

      const size_t en = elec_num*nucl_num;
      const size_t len = l*en;
      const size_t len4 = len*4;
      const size_t c1 = cord_num+1;
      const size_t cn = cord_num*nw;
      const size_t addr0 = en*(m+c1*(k+cn));
      const size_t addr1 = en*(m+c1*nw);

      const double* restrict tmp_c_mkn = &(tmp_c[addr0]);
      const double* restrict tmp_c_mlkn = tmp_c_mkn + len;
      const double* restrict een_rescaled_n_mnw = &(een_rescaled_n[addr1]);
      const double* restrict een_rescaled_n_mlnw = een_rescaled_n_mnw + len;
      const double* restrict dtmp_c_mknw = &(dtmp_c[addr0*4]);
      const double* restrict dtmp_c_mlknw = dtmp_c_mknw + len4;
      const double* restrict een_rescaled_n_gl_mnw = &(een_rescaled_n_gl[addr1*4]);   // ?
      const double* restrict een_rescaled_n_gl_mlnw = een_rescaled_n_gl_mnw + len4;

      for (size_t a = 0; a < (size_t) nucl_num; a++) {
        double cn = c_vector_full[a+n*nucl_num];
        if (cn == 0.0) continue;

        const size_t ishift  = elec_num*a;
        const size_t ishift4 = ishift*4;

        const double* restrict tmp_c_amkn  = tmp_c_mkn  + ishift;
        const double* restrict tmp_c_amlkn = tmp_c_mlkn + ishift;
        const double* restrict een_rescaled_n_amnw  = een_rescaled_n_mnw  + ishift;
        const double* restrict een_rescaled_n_amlnw = een_rescaled_n_mlnw + ishift;
        const double* restrict dtmp_c_0amknw  = dtmp_c_mknw  + ishift4;
        const double* restrict dtmp_c_0amlknw = dtmp_c_mlknw + ishift4;
        const double* restrict een_rescaled_n_gl_0amnw  = een_rescaled_n_gl_mnw  + ishift4;
        const double* restrict een_rescaled_n_gl_0amlnw = een_rescaled_n_gl_mlnw + ishift4;

        const double* restrict dtmp_c_1amknw  = dtmp_c_0amknw  + elec_num;
        const double* restrict dtmp_c_1amlknw = dtmp_c_0amlknw + elec_num;
        const double* restrict dtmp_c_2amknw  = dtmp_c_0amknw  + elec_num2;
        const double* restrict dtmp_c_2amlknw = dtmp_c_0amlknw + elec_num2;
        const double* restrict dtmp_c_3amknw  = dtmp_c_0amknw  + elec_num3;
        const double* restrict dtmp_c_3amlknw = dtmp_c_0amlknw + elec_num3;

        const double* restrict een_rescaled_n_gl_1amnw  = een_rescaled_n_gl_0amnw  + elec_num;
        const double* restrict een_rescaled_n_gl_1amlnw = een_rescaled_n_gl_0amlnw + elec_num;
        const double* restrict een_rescaled_n_gl_2amnw  = een_rescaled_n_gl_0amnw  + elec_num2;
        const double* restrict een_rescaled_n_gl_2amlnw = een_rescaled_n_gl_0amlnw + elec_num2;
        const double* restrict een_rescaled_n_gl_3amnw  = een_rescaled_n_gl_0amnw  + elec_num3;
        const double* restrict een_rescaled_n_gl_3amlnw = een_rescaled_n_gl_0amlnw + elec_num3;

        double* const restrict factor_een_gl_1nw = factor_een_gl_0nw + elec_num;
        double* const restrict factor_een_gl_2nw = factor_een_gl_0nw + elec_num2;
        double* const restrict factor_een_gl_3nw = factor_een_gl_0nw + elec_num3;

        double tmp3[elec_num];

#ifdef HAVE_OPENMP
#pragma omp simd
#endif
        for (size_t j = 0; j < (size_t) elec_num; ++j) {

            factor_een_gl_0nw[j] = factor_een_gl_0nw[j] + cn * (
               dtmp_c_0amknw [j] * een_rescaled_n_amlnw[j]  +
               dtmp_c_0amlknw[j] * een_rescaled_n_amnw [j]  +
               tmp_c_amkn [j] * een_rescaled_n_gl_0amlnw[j] +
               tmp_c_amlkn[j] * een_rescaled_n_gl_0amnw [j] );

            tmp3[j] =
              dtmp_c_0amknw [j]  * een_rescaled_n_gl_0amlnw[j] +
              dtmp_c_0amlknw[j]  * een_rescaled_n_gl_0amnw [j];
          }

#ifdef HAVE_OPENMP
#pragma omp simd
#endif
        for (size_t j = 0; j < (size_t) elec_num; ++j) {

            factor_een_gl_1nw[j] = factor_een_gl_1nw[j] + cn * (
               dtmp_c_1amknw [j] * een_rescaled_n_amlnw[j]  +
               dtmp_c_1amlknw[j] * een_rescaled_n_amnw [j]  +
               tmp_c_amkn [j] * een_rescaled_n_gl_1amlnw[j] +
               tmp_c_amlkn[j] * een_rescaled_n_gl_1amnw [j]);

            tmp3[j] = tmp3[j] +
              dtmp_c_1amknw [j] * een_rescaled_n_gl_1amlnw[j] +
              dtmp_c_1amlknw[j] * een_rescaled_n_gl_1amnw [j];
        }

#ifdef HAVE_OPENMP
#pragma omp simd
#endif
        for (size_t j = 0; j < (size_t) elec_num; ++j) {

            factor_een_gl_2nw[j] = factor_een_gl_2nw[j] + cn * (
               dtmp_c_2amknw [j] * een_rescaled_n_amlnw[j]  +
               dtmp_c_2amlknw[j] * een_rescaled_n_amnw [j]  +
               tmp_c_amkn [j] * een_rescaled_n_gl_2amlnw[j] +
               tmp_c_amlkn[j] * een_rescaled_n_gl_2amnw [j]);

            tmp3[j] = tmp3[j] +
              dtmp_c_2amknw [j] * een_rescaled_n_gl_2amlnw[j] +
              dtmp_c_2amlknw[j] * een_rescaled_n_gl_2amnw [j];
        }

#ifdef HAVE_OPENMP
#pragma omp simd
#endif
        for (size_t j = 0; j < (size_t) elec_num; ++j) {
            factor_een_gl_3nw[j] = factor_een_gl_3nw[j] + cn * (
               dtmp_c_3amknw [j] * een_rescaled_n_amlnw[j] +
               dtmp_c_3amlknw[j] * een_rescaled_n_amnw [j] +
               tmp_c_amkn [j] * een_rescaled_n_gl_3amlnw[j] +
               tmp_c_amlkn[j] * een_rescaled_n_gl_3amnw [j] +
               tmp3[j]*2.0);
          }

      }
    }
  }
  return info;
}

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_een_grad(const qmckl_context context,
                                               const int64_t walk_num,
                                               const int64_t elec_num,
                                               const int64_t nucl_num,
                                               const int64_t cord_num,
                                               const int64_t dim_c_vector,
                                               const double *c_vector_full,
                                               const int64_t *lkpm_combined_index,
                                               const double *tmp_c,
                                               const double *dtmp_c,
                                               const double *een_rescaled_n,
                                               const double *een_rescaled_n_gl,
                                               double* const factor_een_grad)
{
#ifdef HAVE_HPC
  return qmckl_compute_jastrow_champ_factor_een_grad_hpc
#else
  return qmckl_compute_jastrow_champ_factor_een_grad_doc
#endif
    (context, walk_num, elec_num, nucl_num,
     cord_num, dim_c_vector, c_vector_full,
     lkpm_combined_index, tmp_c, dtmp_c,
     een_rescaled_n, een_rescaled_n_gl,
     factor_een_grad);
}

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_een_grad_hpc(const qmckl_context context,
                                                const int64_t walk_num,
                                                const int64_t elec_num,
                                                const int64_t nucl_num,
                                                const int64_t cord_num,
                                                const int64_t dim_c_vector,
                                                const double *c_vector_full,
                                                const int64_t *lkpm_combined_index,
                                                const double *tmp_c,
                                                const double *dtmp_c,
                                                const double *een_rescaled_n,
                                                const double *een_rescaled_n_gl,
                                                double* const factor_een_grad)
{

  int64_t info = QMCKL_SUCCESS;

  if (context == QMCKL_NULL_CONTEXT) return QMCKL_INVALID_CONTEXT;
  if (walk_num <= 0)                 return QMCKL_INVALID_ARG_2;
  if (elec_num <= 0)                 return QMCKL_INVALID_ARG_3;
  if (nucl_num <= 0)                 return QMCKL_INVALID_ARG_4;
  if (cord_num <  0)                 return QMCKL_INVALID_ARG_5;


  if (cord_num == 0) {
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for (size_t nw = 0; nw < (size_t) walk_num; ++nw) {
      double* const restrict factor_een_grad_0nw = &(factor_een_grad[elec_num*3*nw]);
      memset(factor_een_grad_0nw, 0, elec_num*3*sizeof(double));
    }
    return QMCKL_SUCCESS;
  }

  const size_t elec_num2 = elec_num + elec_num;

#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
  for (size_t nw = 0; nw < (size_t) walk_num; ++nw) {
    double* const restrict factor_een_grad_0nw = &(factor_een_grad[elec_num*3*nw]);
    memset(factor_een_grad_0nw, 0, elec_num*3*sizeof(double));
    for (size_t n = 0; n < (size_t) dim_c_vector; ++n) {
      const size_t l = lkpm_combined_index[n];
      const size_t k = lkpm_combined_index[n+  dim_c_vector];
      const size_t m = lkpm_combined_index[n+3*dim_c_vector];

      const size_t en = elec_num*nucl_num;
      const size_t len = l*en;
      const size_t len4 = len*4;
      const size_t cn = cord_num*nw;
      const size_t c1 = cord_num+1;
      const size_t addr0 = en*(m+c1*(k+cn));
      const size_t addr1 = en*(m+c1*nw);

      const double* restrict tmp_c_mkn = &(tmp_c[addr0]);
      const double* restrict tmp_c_mlkn = tmp_c_mkn + len;
      const double* restrict een_rescaled_n_mnw = &(een_rescaled_n[addr1]);
      const double* restrict een_rescaled_n_mlnw = een_rescaled_n_mnw + len;
      const double* restrict dtmp_c_mknw = &(dtmp_c[addr0*4]);
      const double* restrict dtmp_c_mlknw = dtmp_c_mknw + len4;
      const double* restrict een_rescaled_n_gl_mnw = &(een_rescaled_n_gl[addr1*4]);
      const double* restrict een_rescaled_n_gl_mlnw = een_rescaled_n_gl_mnw + len4;

      for (size_t a = 0; a < (size_t) nucl_num; a++) {
        double cn = c_vector_full[a+n*nucl_num];
        if (cn == 0.0) continue;

        const size_t ishift  = elec_num*a;
        const size_t ishift4 = ishift*4;

        const double* restrict tmp_c_amlkn = tmp_c_mlkn + ishift;
        const double* restrict tmp_c_amkn = tmp_c_mkn + ishift;
        const double* restrict een_rescaled_n_amnw  = een_rescaled_n_mnw  + ishift;
        const double* restrict een_rescaled_n_amlnw = een_rescaled_n_mlnw + ishift;
        const double* restrict dtmp_c_0amknw  = dtmp_c_mknw  + ishift4;
        const double* restrict dtmp_c_0amlknw = dtmp_c_mlknw + ishift4;
        const double* restrict een_rescaled_n_gl_0amnw  = een_rescaled_n_gl_mnw  + ishift4;
        const double* restrict een_rescaled_n_gl_0amlnw = een_rescaled_n_gl_mlnw + ishift4;

        const double* restrict dtmp_c_1amknw  = dtmp_c_0amknw  + elec_num;
        const double* restrict dtmp_c_1amlknw = dtmp_c_0amlknw + elec_num;
        const double* restrict dtmp_c_2amknw  = dtmp_c_0amknw  + elec_num2;
        const double* restrict dtmp_c_2amlknw = dtmp_c_0amlknw + elec_num2;

        const double* restrict een_rescaled_n_gl_1amnw  = een_rescaled_n_gl_0amnw  + elec_num;
        const double* restrict een_rescaled_n_gl_1amlnw = een_rescaled_n_gl_0amlnw + elec_num;
        const double* restrict een_rescaled_n_gl_2amnw  = een_rescaled_n_gl_0amnw  + elec_num2;
        const double* restrict een_rescaled_n_gl_2amlnw = een_rescaled_n_gl_0amlnw + elec_num2;

        double* const restrict factor_een_grad_1nw = factor_een_grad_0nw + elec_num;
        double* const restrict factor_een_grad_2nw = factor_een_grad_0nw + elec_num2;

#ifdef HAVE_OPENMP
#pragma omp simd
#endif
        for (size_t j = 0; j < (size_t) elec_num; ++j) {
          factor_een_grad_0nw[j] = factor_een_grad_0nw[j] + cn * (
             dtmp_c_0amknw [j] * een_rescaled_n_amlnw[j]  +
             dtmp_c_0amlknw[j] * een_rescaled_n_amnw [j]  +
             tmp_c_amkn [j] * een_rescaled_n_gl_0amlnw[j] +
             tmp_c_amlkn[j] * een_rescaled_n_gl_0amnw [j] );

        }

#ifdef HAVE_OPENMP
#pragma omp simd
#endif
        for (size_t j = 0; j < (size_t) elec_num; ++j) {
          factor_een_grad_1nw[j] = factor_een_grad_1nw[j] + cn * (
             dtmp_c_1amknw [j] * een_rescaled_n_amlnw[j]  +
             dtmp_c_1amlknw[j] * een_rescaled_n_amnw [j]  +
             tmp_c_amkn [j] * een_rescaled_n_gl_1amlnw[j] +
             tmp_c_amlkn[j] * een_rescaled_n_gl_1amnw [j]);
        }

#ifdef HAVE_OPENMP
#pragma omp simd
#endif
        for (size_t j = 0; j < (size_t) elec_num; ++j) {
          factor_een_grad_2nw[j] = factor_een_grad_2nw[j] + cn * (
             dtmp_c_2amknw [j] * een_rescaled_n_amlnw[j]  +
             dtmp_c_2amlknw[j] * een_rescaled_n_amnw [j]  +
             tmp_c_amkn [j] * een_rescaled_n_gl_2amlnw[j] +
             tmp_c_amlkn[j] * een_rescaled_n_gl_2amnw [j]);
        }

      }
    }
  }
  return info;
}

qmckl_exit_code
qmckl_get_jastrow_champ_factor_een_pderiv(qmckl_context context,
                             double* const factor_een_pderiv,
                             const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_jastrow_champ_factor_een_pderiv(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (factor_een_pderiv == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_champ_factor_een_pderiv",
                           "Null pointer");
  }

  const int64_t sze = ctx->jastrow_champ.dim_c_vector * ctx->jastrow_champ.type_nucl_num;

  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_factor_een_pderiv",
                           "Array too small. Expected dim_c_vector * type_nucl_num");
  }

  memcpy(factor_een_pderiv, ctx->jastrow_champ.factor_een_pderiv, sze*sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_jastrow_champ_factor_een_pderiv(qmckl_context context)
{

  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (ctx->jastrow_champ.cord_num == 0) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_provide_jastrow_champ_factor_een_gl_pderiv",
                           "Expected cord_num > 0");
  }

  if (ctx->jastrow_champ.cord_num > 0) {

    /* Check if en rescaled distance is provided */
    rc = qmckl_provide_een_rescaled_e(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if en rescaled distance derivatives is provided */
    rc = qmckl_provide_een_rescaled_n(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if en rescaled distance derivatives is provided */
    rc = qmckl_provide_jastrow_champ_c_vector_full(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if en rescaled distance derivatives is provided */
    rc = qmckl_provide_lkpm_combined_index(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if tmp_c is provided */
    rc = qmckl_provide_tmp_c(context);
    if(rc != QMCKL_SUCCESS) return rc;
  }


  /* Compute if necessary */
  if (ctx->date > ctx->jastrow_champ.factor_een_pderiv_date) {

    /* Allocate array */
    if (ctx->jastrow_champ.factor_een_pderiv == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->jastrow_champ.dim_c_vector * ctx->jastrow_champ.type_nucl_num * sizeof(double);
      double* factor_een_pderiv = (double*) qmckl_malloc(context, mem_info);

      if (factor_een_pderiv == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_jastrow_champ_factor_een_pderiv",
                               NULL);
      }
      ctx->jastrow_champ.factor_een_pderiv = factor_een_pderiv;
    }

    rc = qmckl_compute_jastrow_champ_factor_een_pderiv(context,
                                                       ctx->electron.walker.num,
                                                       ctx->electron.num,
                                                       ctx->nucleus.num,
                                                       ctx->jastrow_champ.type_nucl_num,
                                                       ctx->jastrow_champ.type_nucl_vector,
                                                       ctx->jastrow_champ.cord_num,
                                                       ctx->jastrow_champ.dim_c_vector,
                                                       ctx->jastrow_champ.c_vector_full,
                                                       ctx->jastrow_champ.lkpm_combined_index,
                                                       ctx->jastrow_champ.tmp_c,
                                                       ctx->jastrow_champ.een_rescaled_n,
                                                       ctx->jastrow_champ.factor_een_pderiv);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow_champ.factor_een_pderiv_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_een_pderiv (const qmckl_context context,
                          const int64_t walk_num,
                          const int64_t elec_num,
                          const int64_t nucl_num,
                          const int64_t type_nucl_num,
                          const int64_t* type_nucl_vector,
                          const int64_t cord_num,
                          const int64_t dim_c_vector,
                          const double* c_vector_full,
                          const int64_t* lkpm_combined_index,
                          const double* tmp_c,
                          const double* een_rescaled_n,
                          double* const factor_een_pderiv )
{
  return qmckl_compute_jastrow_champ_factor_een_pderiv_doc
    (context, walk_num, elec_num, nucl_num, type_nucl_num, type_nucl_vector, 
    cord_num, dim_c_vector, c_vector_full, lkpm_combined_index, tmp_c, 
    een_rescaled_n, factor_een_pderiv );
}

qmckl_exit_code
qmckl_get_jastrow_champ_factor_een_gl_pderiv(qmckl_context context,
                                      double* const factor_een_gl_pderiv,
                                      const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_jastrow_champ_factor_een_gl_pderiv(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (factor_een_gl_pderiv == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_champ_factor_een_gl_pderiv",
                           "Null pointer");
  }

  const int64_t sze = 4 * ctx->electron.num * ctx->jastrow_champ.dim_c_vector * ctx->jastrow_champ.type_nucl_num;

  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_factor_een_gl_pderiv",
                           "Array too small. Expected 4*elec_num*dim_c_vector*type_nucl_num");
  }

  memcpy(factor_een_gl_pderiv, ctx->jastrow_champ.factor_een_gl_pderiv, sze*sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_jastrow_champ_factor_een_gl_pderiv(qmckl_context context)
{

  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (ctx->jastrow_champ.cord_num == 0) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_provide_jastrow_champ_factor_een_gl_pderiv",
                           "Expected cord_num > 0");
  }

  if (ctx->jastrow_champ.cord_num > 0) {

    /* Check if en rescaled distance is provided */
    rc = qmckl_provide_een_rescaled_e(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if en rescaled distance derivatives is provided */
    rc = qmckl_provide_een_rescaled_n(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if en rescaled distance is provided */
    rc = qmckl_provide_een_rescaled_e_gl(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if en rescaled distance derivatives is provided */
    rc = qmckl_provide_een_rescaled_n_gl(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if en rescaled distance derivatives is provided */
    rc = qmckl_provide_jastrow_champ_c_vector_full(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if en rescaled distance derivatives is provided */
    rc = qmckl_provide_lkpm_combined_index(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if tmp_c is provided */
    rc = qmckl_provide_tmp_c(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if dtmp_c is provided */
    rc = qmckl_provide_dtmp_c(context);
    if(rc != QMCKL_SUCCESS) return rc;
  }


    /* Compute if necessary */
  if (ctx->date > ctx->jastrow_champ.factor_een_gl_pderiv_date) {

    /* Allocate array */
    if (ctx->jastrow_champ.factor_een_gl_pderiv == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = 4 * ctx->electron.num * ctx->jastrow_champ.dim_c_vector * ctx->jastrow_champ.type_nucl_num * sizeof(double);
      double* factor_een_gl_pderiv = (double*) qmckl_malloc(context, mem_info);

      if (factor_een_gl_pderiv == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_jastrow_champ_factor_een_gl_pderiv",
                               NULL);
      }
      ctx->jastrow_champ.factor_een_gl_pderiv = factor_een_gl_pderiv;
    }

    rc = qmckl_compute_jastrow_champ_factor_een_gl_pderiv(context,
                                                   ctx->electron.walker.num,
                                                   ctx->electron.num,
                                                   ctx->nucleus.num,
                                                   ctx->jastrow_champ.type_nucl_num,
                                                   ctx->jastrow_champ.type_nucl_vector,
                                                   ctx->jastrow_champ.cord_num,
                                                   ctx->jastrow_champ.dim_c_vector,
                                                   ctx->jastrow_champ.c_vector_full,
                                                   ctx->jastrow_champ.lkpm_combined_index,
                                                   ctx->jastrow_champ.tmp_c,
                                                   ctx->jastrow_champ.dtmp_c,
                                                   ctx->jastrow_champ.een_rescaled_n,
                                                   ctx->jastrow_champ.een_rescaled_n_gl,
                                                   ctx->jastrow_champ.factor_een_gl_pderiv);

    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->jastrow_champ.factor_een_gl_pderiv_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_een_gl_pderiv(const qmckl_context context,
                                               const int64_t walk_num,
                                               const int64_t elec_num,
                                               const int64_t nucl_num,
                                               const int64_t type_nucl_num,
                                               const int64_t* type_nucl_vector,
                                               const int64_t cord_num,
                                               const int64_t dim_c_vector,
                                               const double *c_vector_full,
                                               const int64_t *lkpm_combined_index,
                                               const double *tmp_c,
                                               const double *dtmp_c,
                                               const double *een_rescaled_n,
                                               const double *een_rescaled_n_gl,
                                               double* const factor_een_gl_pderiv)
{
  return qmckl_compute_jastrow_champ_factor_een_gl_pderiv_doc
    (context, walk_num, elec_num, nucl_num,
     type_nucl_num, type_nucl_vector,
     cord_num, dim_c_vector, c_vector_full,
     lkpm_combined_index, tmp_c, dtmp_c,
     een_rescaled_n, een_rescaled_n_gl,
     factor_een_gl_pderiv);
}

qmckl_exit_code
qmckl_get_jastrow_champ_value(qmckl_context context,
                            double* const value,
                            const int64_t size_max)
{
  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_jastrow_champ_value",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  rc = qmckl_provide_jastrow_champ_value(context);
  if (rc != QMCKL_SUCCESS) return rc;

  if (value == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_champ_value",
                           "Null pointer");
  }

  const int64_t sze = ctx->electron.walker.num;

  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_value",
                           "Array too small. Expected walk_num");
  }

  memcpy(value, ctx->jastrow_champ.value, sze*sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_jastrow_champ_value(qmckl_context context)
{

  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_provide_jastrow_champ_value",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->jastrow_champ.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_provide_jastrow_champ_value",
                           NULL);
  }


  rc = qmckl_provide_jastrow_champ_factor_ee(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_jastrow_champ_factor_en(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_jastrow_champ_factor_een(context);
  if (rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow_champ.value_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      if (ctx->jastrow_champ.value != NULL) {
        rc = qmckl_free(context, ctx->jastrow_champ.value);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_jastrow_champ_value",
                                 "Unable to free ctx->jastrow_champ.value");
        }
        ctx->jastrow_champ.value = NULL;
      }
    }

    /* Allocate array */
    if (ctx->jastrow_champ.value == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.walker.num * sizeof(double);
      double* value = (double*) qmckl_malloc(context, mem_info);

      if (value == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_jastrow_champ_value",
                               NULL);
      }
      ctx->jastrow_champ.value = value;
    }

    rc = qmckl_compute_jastrow_champ_value_doc(context,
                                               ctx->electron.walker.num,
                                               ctx->jastrow_champ.factor_ee,
                                               ctx->jastrow_champ.factor_en,
                                               ctx->jastrow_champ.factor_een,
                                               ctx->jastrow_champ.value);

    ctx->jastrow_champ.value_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

inline qmckl_exit_code
qmckl_compute_jastrow_champ_value_hpc (
          const qmckl_context context,
          const int64_t walk_num,
          const double* factor_ee,
          const double* factor_en,
          const double* factor_een,
          double* const value)
{

  if (context    == QMCKL_NULL_CONTEXT) return QMCKL_INVALID_CONTEXT;
  if (walk_num   <= 0                 ) return QMCKL_INVALID_ARG_2;
  if (factor_ee  == NULL              ) return QMCKL_INVALID_ARG_3;
  if (factor_en  == NULL              ) return QMCKL_INVALID_ARG_4;
  if (factor_een == NULL              ) return QMCKL_INVALID_ARG_5;
  if (value      == NULL              ) return QMCKL_INVALID_ARG_6;

  for (int64_t i = 0; i < walk_num; ++i) {
    value[i] = exp(factor_ee[i] + factor_en[i] + factor_een[i]);
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_compute_jastrow_champ_value (
          const qmckl_context context,
          const int64_t walk_num,
          const double* factor_ee,
          const double* factor_en,
          const double* factor_een,
          double* const value)
{

#ifdef HAVE_HPC
  return qmckl_compute_jastrow_champ_value_hpc
#else
  return qmckl_compute_jastrow_champ_value_doc
#endif
    (context, walk_num, factor_ee, factor_en, factor_een, value);
}

qmckl_exit_code
qmckl_get_jastrow_champ_gl(qmckl_context context,
                            double* const gl,
                            const int64_t size_max)
{
  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_jastrow_champ_gl",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  rc = qmckl_provide_jastrow_champ_gl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  if (gl == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_champ_gl",
                           "Null pointer");
  }

  const int64_t sze = 4 * ctx->electron.walker.num * ctx->electron.num;

  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_gl",
                           "Array too small. Expected walk_num*elec_num*4");
  }

  memcpy(gl, ctx->jastrow_champ.gl, sze*sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_jastrow_champ_grad(qmckl_context context,
                            double* const grad,
                            const int64_t size_max)
{
  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_jastrow_champ_grad",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  rc = qmckl_provide_jastrow_champ_grad(context);
  if (rc != QMCKL_SUCCESS) return rc;

  int64_t sze = 3 * ctx->electron.walker.num * ctx->electron.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_grad",
                           "Array too small. Expected walker.num * electron.num * 3");
  }
  memcpy(grad, ctx->jastrow_champ.grad, sze*sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_jastrow_champ_gl(qmckl_context context)
{

  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_provide_jastrow_champ_gl",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->jastrow_champ.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_provide_jastrow_champ_gl",
                           NULL);
  }


  rc = qmckl_provide_jastrow_champ_value(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_jastrow_champ_factor_ee_gl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_jastrow_champ_factor_en_gl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_jastrow_champ_factor_een_gl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow_champ.gl_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      if (ctx->jastrow_champ.gl != NULL) {
        rc = qmckl_free(context, ctx->jastrow_champ.gl);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_jastrow_champ_gl",
                                 "Unable to free ctx->jastrow_champ.gl");
        }
        ctx->jastrow_champ.gl = NULL;
      }
    }

    /* Allocate array */
    if (ctx->jastrow_champ.gl == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.walker.num * ctx->electron.num * 4 * sizeof(double);
      double* gl = (double*) qmckl_malloc(context, mem_info);

      if (gl == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_jastrow_champ_gl",
                               NULL);
      }
      ctx->jastrow_champ.gl = gl;
    }

    rc = qmckl_compute_jastrow_champ_gl_doc(context,
                                            ctx->electron.walker.num,
                                            ctx->electron.num,
                                            ctx->jastrow_champ.value,
                                            ctx->jastrow_champ.factor_ee_gl,
                                            ctx->jastrow_champ.factor_en_gl,
                                            ctx->jastrow_champ.factor_een_gl,
                                            ctx->jastrow_champ.gl);

    ctx->jastrow_champ.gl_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_jastrow_champ_grad(qmckl_context context)
{

  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_provide_jastrow_champ_grad",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->jastrow_champ.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_provide_jastrow_champ_grad",
                           NULL);
  }


  rc = qmckl_provide_jastrow_champ_value(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_jastrow_champ_factor_ee_gl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_jastrow_champ_factor_en_gl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_jastrow_champ_factor_een_grad(context);
  if (rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->date > ctx->jastrow_champ.grad_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      if (ctx->jastrow_champ.grad != NULL) {
        rc = qmckl_free(context, ctx->jastrow_champ.grad);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_jastrow_champ_grad",
                                 "Unable to free ctx->jastrow_champ.grad");
        }
        ctx->jastrow_champ.grad = NULL;
      }
    }

    /* Allocate array */
    if (ctx->jastrow_champ.grad == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.walker.num * ctx->electron.num * 3 * sizeof(double);
      double* grad = (double*) qmckl_malloc(context, mem_info);

      if (grad == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_jastrow_champ_grad",
                               NULL);
      }
      ctx->jastrow_champ.grad = grad;
    }

    rc = qmckl_compute_jastrow_champ_grad_doc(context,
                                              ctx->electron.walker.num,
                                              ctx->electron.num,
                                              ctx->jastrow_champ.value,
                                              ctx->jastrow_champ.factor_ee_gl,
                                              ctx->jastrow_champ.factor_en_gl,
                                              ctx->jastrow_champ.factor_een_grad,
                                              ctx->jastrow_champ.grad);

    ctx->jastrow_champ.grad_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

inline qmckl_exit_code
qmckl_compute_jastrow_champ_gl_hpc (const qmckl_context context,
                                    const int64_t walk_num,
                                    const int64_t elec_num,
                                    const double* value,
                                    const double* gl_ee,
                                    const double* gl_en,
                                    const double* gl_een,
                                    double* const gl)
{

  if (context    == QMCKL_NULL_CONTEXT) return QMCKL_INVALID_CONTEXT;
  if (walk_num   <= 0                 ) return QMCKL_INVALID_ARG_2;
  if (elec_num   <= 0                 ) return QMCKL_INVALID_ARG_3;
  if (value      == NULL              ) return QMCKL_INVALID_ARG_4;
  if (gl_ee      == NULL              ) return QMCKL_INVALID_ARG_5;
  if (gl_en      == NULL              ) return QMCKL_INVALID_ARG_6;
  if (gl_een     == NULL              ) return QMCKL_INVALID_ARG_7;
  if (gl         == NULL              ) return QMCKL_INVALID_ARG_8;

  for (int64_t k = 0; k < walk_num; ++k) {
    for (int64_t j = 0; j < 4; ++j) {
      for (int64_t i = 0; i < elec_num; ++i) {
        gl[i + elec_num*(j + k*4)] = gl_ee[i + elec_num*(j + k*4)] +
          gl_en[i + elec_num*(j + k*4)] + gl_een[i + elec_num*(j + k*4)];
      }
    }
    for (int64_t i = 0; i < elec_num; ++i) {
      gl[i + elec_num*(3 + walk_num*4)] +=
        gl_ee[i + elec_num*(0 + k*4)] * gl_ee[i + elec_num*(0 + k*4)] +
        gl_ee[i + elec_num*(1 + k*4)] * gl_ee[i + elec_num*(1 + k*4)] +
        gl_ee[i + elec_num*(2 + k*4)] * gl_ee[i + elec_num*(2 + k*4)];
    }
    for (int64_t j = 0; j < 4; ++j) {
      for (int64_t i = 0; i < elec_num; ++i) {
        gl[i + elec_num*(j + k*4)] *= value[k];
      }
    }
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_compute_jastrow_champ_gl (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const double* value,
      const double* gl_ee,
      const double* gl_en,
      const double* gl_een,
      double* const gl)
{

#ifdef HAVE_HPC
  return qmckl_compute_jastrow_champ_gl_hpc
#else
  return qmckl_compute_jastrow_champ_gl_doc
#endif
    (context, walk_num, elec_num, value, gl_ee, gl_en, gl_een, gl);
}

inline qmckl_exit_code
qmckl_compute_jastrow_champ_grad_hpc (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const double* value,
      const double* gl_ee,
      const double* gl_en,
      const double* grad_een,
      double* const grad)
{

  if (context    == QMCKL_NULL_CONTEXT) return QMCKL_INVALID_CONTEXT;
  if (walk_num   <= 0                 ) return QMCKL_INVALID_ARG_2;
  if (elec_num   <= 0                 ) return QMCKL_INVALID_ARG_3;
  if (value      == NULL              ) return QMCKL_INVALID_ARG_4;
  if (gl_ee      == NULL              ) return QMCKL_INVALID_ARG_5;
  if (gl_en      == NULL              ) return QMCKL_INVALID_ARG_6;
  if (grad_een   == NULL              ) return QMCKL_INVALID_ARG_7;
  if (grad       == NULL              ) return QMCKL_INVALID_ARG_8;

  for (int64_t k = 0; k < walk_num; ++k) {
    for (int64_t j = 0; j < 3; ++j) {
      for (int64_t i = 0; i < elec_num; ++i) {
        grad[i + elec_num*(j + k*3)] = ( gl_ee[i + elec_num*(j + k*4)] +
          gl_en[i + elec_num*(j + k*4)] + grad_een[i + elec_num*(j + k*3)] )* value[k];
      }
    }
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_compute_jastrow_champ_grad (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const double* value,
      const double* gl_ee,
      const double* gl_en,
      const double* grad_een,
      double* const grad)
{

#ifdef HAVE_HPC
  return qmckl_compute_jastrow_champ_grad_hpc
#else
  return qmckl_compute_jastrow_champ_grad_doc
#endif
    (context, walk_num, elec_num, value, gl_ee, gl_en, grad_een, grad);
}
