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
#include <stdio.h>
#include <math.h>

#include "qmckl.h"
#include "qmckl_context_private_type.h"
#include "qmckl_memory_private_type.h"
#include "qmckl_blas_private_type.h"
#include "qmckl_memory_private_func.h"
#include "qmckl_ao_private_type.h"
#include "qmckl_ao_private_func.h"

qmckl_exit_code qmckl_init_ao_basis(qmckl_context context) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_init_ao_basis",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  ctx->ao_basis.uninitialized = (1 << 14) - 1;

  /* Default values */
  ctx->ao_basis.ao_cartesian = true;

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_set_ao_basis_type(qmckl_context context,
                        const char basis_type)
{
  int32_t mask = 1;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
      return qmckl_failwith( context,
                             QMCKL_NULL_CONTEXT,
                             "qmckl_set_ao_*",
                             NULL);
  }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  
  if (mask != 0 && !(ctx->ao_basis.uninitialized & mask)) {
      return qmckl_failwith( context,
                             QMCKL_ALREADY_SET,
                             "qmckl_set_ao_*",
                             NULL);
   }

  if (basis_type != 'G' && basis_type != 'S') {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_set_ao_basis_type",
                           NULL);
  }

  ctx->ao_basis.type = basis_type;

  ctx->ao_basis.uninitialized &= ~mask;
  ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
  if (ctx->ao_basis.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_basis(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
   }
  
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_set_ao_basis_shell_num (qmckl_context context,
                              const int64_t shell_num)
{
  int32_t mask = 1 << 1;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
      return qmckl_failwith( context,
                             QMCKL_NULL_CONTEXT,
                             "qmckl_set_ao_*",
                             NULL);
  }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  
  if (mask != 0 && !(ctx->ao_basis.uninitialized & mask)) {
      return qmckl_failwith( context,
                             QMCKL_ALREADY_SET,
                             "qmckl_set_ao_*",
                             NULL);
   }

  if (shell_num <= 0) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_set_ao_basis_shell_num",
                           "shell_num <= 0");
  }

  const int64_t prim_num = ctx->ao_basis.prim_num;

  if (0L < prim_num && prim_num < shell_num) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_set_ao_basis_shell_num",
                           "shell_num > prim_num");
  }

  ctx->ao_basis.shell_num = shell_num;

  ctx->ao_basis.uninitialized &= ~mask;
  ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
  if (ctx->ao_basis.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_basis(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
   }
  
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_set_ao_basis_prim_num (qmckl_context context,
                             const int64_t prim_num)
{
  int32_t mask = 1 << 2;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
      return qmckl_failwith( context,
                             QMCKL_NULL_CONTEXT,
                             "qmckl_set_ao_*",
                             NULL);
  }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  
  if (mask != 0 && !(ctx->ao_basis.uninitialized & mask)) {
      return qmckl_failwith( context,
                             QMCKL_ALREADY_SET,
                             "qmckl_set_ao_*",
                             NULL);
   }

  if (prim_num <= 0) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_set_ao_basis_shell_num",
                           "prim_num must be positive");
  }

  const int64_t shell_num = ctx->ao_basis.shell_num;

  if (shell_num <= 0L) {
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "qmckl_set_ao_basis_shell_num",
                           "shell_num is not set");
  }

  if (prim_num < shell_num) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_set_ao_basis_shell_num",
                           "prim_num < shell_num");
  }

  ctx->ao_basis.prim_num = prim_num;

  ctx->ao_basis.uninitialized &= ~mask;
  ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
  if (ctx->ao_basis.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_basis(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
   }
  
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_set_ao_basis_nucleus_shell_num (qmckl_context context,
                                      const int64_t* nucleus_shell_num,
                                      const int64_t size_max)
{
  int32_t mask = 1 << 3;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
      return qmckl_failwith( context,
                             QMCKL_NULL_CONTEXT,
                             "qmckl_set_ao_*",
                             NULL);
  }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  
  if (mask != 0 && !(ctx->ao_basis.uninitialized & mask)) {
      return qmckl_failwith( context,
                             QMCKL_ALREADY_SET,
                             "qmckl_set_ao_*",
                             NULL);
   }

  const int64_t nucl_num = ctx->nucleus.num;

  if (nucl_num <= 0L) {
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "qmckl_set_ao_basis_nucleus_shell_num",
                           "shell_num is not set");
  }

  if (size_max < nucl_num) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_set_ao_basis_nucleus_shell_num",
                           "input array too small");
  }

  if (ctx->ao_basis.nucleus_shell_num != NULL) {
    qmckl_exit_code rc = qmckl_free(context, ctx->ao_basis.nucleus_shell_num);
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context, rc,
                             "qmckl_set_ao_basis_nucleus_shell_num",
                             NULL);
    }
  }

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = nucl_num * sizeof(int64_t);
  int64_t* new_array = (int64_t*) qmckl_malloc(context, mem_info);

  if (new_array == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_set_ao_basis_nucleus_shell_num",
                           NULL);
  }

  memcpy(new_array, nucleus_shell_num, mem_info.size);

  ctx->ao_basis.nucleus_shell_num = new_array;

  ctx->ao_basis.uninitialized &= ~mask;
  ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
  if (ctx->ao_basis.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_basis(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
   }
  
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_set_ao_basis_nucleus_index (qmckl_context context,
                                  const int64_t* nucleus_index,
                                  const int64_t size_max)
{
  int32_t mask = 1 << 4;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
      return qmckl_failwith( context,
                             QMCKL_NULL_CONTEXT,
                             "qmckl_set_ao_*",
                             NULL);
  }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  
  if (mask != 0 && !(ctx->ao_basis.uninitialized & mask)) {
      return qmckl_failwith( context,
                             QMCKL_ALREADY_SET,
                             "qmckl_set_ao_*",
                             NULL);
   }

  const int64_t nucl_num = ctx->nucleus.num;

  if (nucl_num <= 0L) {
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "qmckl_set_ao_basis_nucleus_index",
                           "nucl_num is not set");
  }

  if (size_max < nucl_num) {
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "qmckl_set_ao_basis_nucleus_index",
                           "input array too small");
  }

  if (ctx->ao_basis.nucleus_index != NULL) {
    qmckl_exit_code rc = qmckl_free(context, ctx->ao_basis.nucleus_index);
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context, rc,
                             "qmckl_set_ao_basis_nucleus_index",
                             NULL);
    }
  }

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = nucl_num * sizeof(int64_t);
  int64_t* new_array = (int64_t*) qmckl_malloc(context, mem_info);

  if (new_array == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_set_ao_basis_nucleus_index",
                           NULL);
  }

  memcpy(new_array, nucleus_index, mem_info.size);

  ctx->ao_basis.nucleus_index = new_array;

  ctx->ao_basis.uninitialized &= ~mask;
  ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
  if (ctx->ao_basis.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_basis(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
   }
  
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_set_ao_basis_shell_ang_mom (qmckl_context context,
                                  const int32_t* shell_ang_mom,
                                  const int64_t size_max)
{
  int32_t mask = 1 << 5;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
      return qmckl_failwith( context,
                             QMCKL_NULL_CONTEXT,
                             "qmckl_set_ao_*",
                             NULL);
  }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  
  if (mask != 0 && !(ctx->ao_basis.uninitialized & mask)) {
      return qmckl_failwith( context,
                             QMCKL_ALREADY_SET,
                             "qmckl_set_ao_*",
                             NULL);
   }

  const int64_t shell_num = ctx->ao_basis.shell_num;

  if (shell_num == 0L) {
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "qmckl_set_ao_basis_shell_ang_mom",
                           "shell_num is not set");
  }

  if (size_max < shell_num) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_set_ao_basis_shell_ang_mom",
                           "input array too small");
  }

  if (ctx->ao_basis.shell_ang_mom != NULL) {
    qmckl_exit_code rc = qmckl_free(context, ctx->ao_basis.shell_ang_mom);
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context, rc,
                             "qmckl_set_ao_basis_shell_ang_mom",
                             NULL);
    }
  }


  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = shell_num * sizeof(int32_t);
  int32_t * new_array = (int32_t*) qmckl_malloc(context, mem_info);

  if (new_array == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_set_ao_basis_shell_ang_mom",
                           NULL);
  }

  memcpy(new_array, shell_ang_mom, mem_info.size);

  ctx->ao_basis.shell_ang_mom = new_array;

  ctx->ao_basis.uninitialized &= ~mask;
  ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
  if (ctx->ao_basis.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_basis(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
   }
  
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_set_ao_basis_shell_prim_num (qmckl_context context,
                                   const int64_t* shell_prim_num,
                                   const int64_t size_max)
{
  int32_t mask = 1 << 6;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
      return qmckl_failwith( context,
                             QMCKL_NULL_CONTEXT,
                             "qmckl_set_ao_*",
                             NULL);
  }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  
  if (mask != 0 && !(ctx->ao_basis.uninitialized & mask)) {
      return qmckl_failwith( context,
                             QMCKL_ALREADY_SET,
                             "qmckl_set_ao_*",
                             NULL);
   }

  const int64_t shell_num = ctx->ao_basis.shell_num;

  if (shell_num <= 0L) {
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "qmckl_set_ao_basis_shell_prim_num",
                           "shell_num is not set");
  }

  if (size_max < shell_num) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_set_ao_basis_shell_prim_num",
                           "input array too small");
  }

  if (ctx->ao_basis.shell_prim_num != NULL) {
    qmckl_exit_code rc = qmckl_free(context, ctx->ao_basis.shell_prim_num);
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context, rc,
                             "qmckl_set_ao_basis_shell_prim_num",
                             NULL);
    }
  }


  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = shell_num * sizeof(int64_t);
  int64_t* new_array = (int64_t*) qmckl_malloc(context, mem_info);

  if (new_array == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_set_ao_basis_shell_prim_num",
                           NULL);
  }

  memcpy(new_array, shell_prim_num, mem_info.size);

  ctx->ao_basis.shell_prim_num = new_array;

  ctx->ao_basis.uninitialized &= ~mask;
  ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
  if (ctx->ao_basis.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_basis(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
   }
  
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_set_ao_basis_shell_prim_index (qmckl_context context,
                                     const int64_t* shell_prim_index,
                                     const int64_t size_max)
{
  int32_t mask = 1 << 7;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
      return qmckl_failwith( context,
                             QMCKL_NULL_CONTEXT,
                             "qmckl_set_ao_*",
                             NULL);
  }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  
  if (mask != 0 && !(ctx->ao_basis.uninitialized & mask)) {
      return qmckl_failwith( context,
                             QMCKL_ALREADY_SET,
                             "qmckl_set_ao_*",
                             NULL);
   }

  const int64_t shell_num = ctx->ao_basis.shell_num;

  if (shell_num <= 0L) {
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "qmckl_set_ao_basis_shell_prim_index",
                           "shell_num is not set");
  }

  if (size_max < shell_num) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_set_ao_basis_shell_prim_index",
                           "input array too small");
  }

  if (ctx->ao_basis.shell_prim_index != NULL) {
    qmckl_exit_code rc = qmckl_free(context, ctx->ao_basis.shell_prim_index);
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context, rc,
                             "qmckl_set_ao_basis_shell_prim_index",
                             NULL);
    }
  }

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = shell_num * sizeof(int64_t);
  int64_t* new_array = (int64_t*) qmckl_malloc(context, mem_info);

  if (new_array == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_set_ao_basis_shell_prim_index",
                           NULL);
  }

  memcpy(new_array, shell_prim_index, mem_info.size);

  ctx->ao_basis.shell_prim_index = new_array;

  ctx->ao_basis.uninitialized &= ~mask;
  ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
  if (ctx->ao_basis.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_basis(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
   }
  
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_set_ao_basis_shell_factor (qmckl_context context,
                                 const double* shell_factor,
                                 const int64_t size_max)
{
  int32_t mask = 1 << 8;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
      return qmckl_failwith( context,
                             QMCKL_NULL_CONTEXT,
                             "qmckl_set_ao_*",
                             NULL);
  }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  
  if (mask != 0 && !(ctx->ao_basis.uninitialized & mask)) {
      return qmckl_failwith( context,
                             QMCKL_ALREADY_SET,
                             "qmckl_set_ao_*",
                             NULL);
   }

  const int64_t shell_num = ctx->ao_basis.shell_num;

  if (shell_num <= 0L) {
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "qmckl_set_ao_basis_shell_factor",
                           "shell_num is not set");
  }

  if (size_max < shell_num) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_set_ao_basis_shell_factor",
                           "input array too small");
  }

  if (ctx->ao_basis.shell_factor != NULL) {
    qmckl_exit_code rc = qmckl_free(context, ctx->ao_basis.shell_factor);
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context, rc,
                             "qmckl_set_ao_basis_shell_factor",
                             NULL);
    }
  }

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = shell_num * sizeof(double);
  double* new_array = (double*) qmckl_malloc(context, mem_info);

  if (new_array == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_set_ao_basis_shell_factor",
                           NULL);
  }

  memcpy(new_array, shell_factor, mem_info.size);

  ctx->ao_basis.shell_factor = new_array;

  ctx->ao_basis.uninitialized &= ~mask;
  ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
  if (ctx->ao_basis.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_basis(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
   }
  
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_set_ao_basis_exponent (qmckl_context context,
                             const double* exponent,
                             const int64_t size_max)
{
  int32_t mask = 1 << 9;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
      return qmckl_failwith( context,
                             QMCKL_NULL_CONTEXT,
                             "qmckl_set_ao_*",
                             NULL);
  }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  
  if (mask != 0 && !(ctx->ao_basis.uninitialized & mask)) {
      return qmckl_failwith( context,
                             QMCKL_ALREADY_SET,
                             "qmckl_set_ao_*",
                             NULL);
   }

  const int64_t prim_num = ctx->ao_basis.prim_num;

  if (prim_num <= 0L) {
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "qmckl_set_ao_basis_exponent",
                           "prim_num is not set");
  }

  if (size_max < prim_num) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_set_ao_basis_exponent",
                           "input array too small");
  }

  if (ctx->ao_basis.exponent != NULL) {
    qmckl_exit_code rc = qmckl_free(context, ctx->ao_basis.exponent);
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context, rc,
                             "qmckl_set_ao_basis_exponent",
                             NULL);
    }
  }

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = prim_num * sizeof(double);
  double* new_array = (double*) qmckl_malloc(context, mem_info);

  if (new_array == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_set_ao_basis_exponent",
                           NULL);
  }

  memcpy(new_array, exponent, mem_info.size);

  ctx->ao_basis.exponent = new_array;

  ctx->ao_basis.uninitialized &= ~mask;
  ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
  if (ctx->ao_basis.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_basis(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
   }
  
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_set_ao_basis_coefficient (qmckl_context context,
                                const double* coefficient,
                                const int64_t size_max)
{
  int32_t mask = 1 << 10;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
      return qmckl_failwith( context,
                             QMCKL_NULL_CONTEXT,
                             "qmckl_set_ao_*",
                             NULL);
  }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  
  if (mask != 0 && !(ctx->ao_basis.uninitialized & mask)) {
      return qmckl_failwith( context,
                             QMCKL_ALREADY_SET,
                             "qmckl_set_ao_*",
                             NULL);
   }

  const int64_t prim_num = ctx->ao_basis.prim_num;

  if (prim_num <= 0L) {
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "qmckl_set_ao_basis_coefficient",
                           "prim_num is not set");
  }

  if (size_max < prim_num) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_set_ao_basis_coefficient",
                           "input array too small");
  }

  if (ctx->ao_basis.coefficient != NULL) {
    qmckl_exit_code rc = qmckl_free(context, ctx->ao_basis.coefficient);
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context, rc,
                             "qmckl_set_ao_basis_coefficient",
                             NULL);
    }
  }

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = prim_num * sizeof(double);
  double* new_array = (double*) qmckl_malloc(context, mem_info);

  if (new_array == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_set_ao_basis_coefficient",
                           NULL);
  }

  memcpy(new_array, coefficient, mem_info.size);

  ctx->ao_basis.coefficient = new_array;

  ctx->ao_basis.uninitialized &= ~mask;
  ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
  if (ctx->ao_basis.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_basis(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
   }
  
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_set_ao_basis_prim_factor (qmckl_context context,
                                const double* prim_factor,
                                const int64_t size_max)
{
  int32_t mask = 1 << 11;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
      return qmckl_failwith( context,
                             QMCKL_NULL_CONTEXT,
                             "qmckl_set_ao_*",
                             NULL);
  }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  
  if (mask != 0 && !(ctx->ao_basis.uninitialized & mask)) {
      return qmckl_failwith( context,
                             QMCKL_ALREADY_SET,
                             "qmckl_set_ao_*",
                             NULL);
   }

  const int64_t prim_num = ctx->ao_basis.prim_num;

  if (prim_num <= 0L) {
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "qmckl_set_ao_basis_prim_factor",
                           "prim_num is not set");
  }

  if (size_max < prim_num) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_set_ao_basis_prim_factor",
                           "input array too small");
  }

  if (ctx->ao_basis.prim_factor != NULL) {
    qmckl_exit_code rc = qmckl_free(context, ctx->ao_basis.prim_factor);
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context, rc,
                             "qmckl_set_ao_basis_prim_factor",
                             NULL);
    }
  }

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = prim_num * sizeof(double);
  double* new_array = (double*) qmckl_malloc(context, mem_info);

  if (new_array == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_set_ao_basis_prim_factor",
                           NULL);
  }

  memcpy(new_array, prim_factor, mem_info.size);

  ctx->ao_basis.prim_factor = new_array;

  ctx->ao_basis.uninitialized &= ~mask;
  ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
  if (ctx->ao_basis.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_basis(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
   }
  
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_set_ao_basis_ao_num (qmckl_context context,
                           const int64_t ao_num)
{
  int32_t mask = 1 << 12;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
      return qmckl_failwith( context,
                             QMCKL_NULL_CONTEXT,
                             "qmckl_set_ao_*",
                             NULL);
  }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  
  if (mask != 0 && !(ctx->ao_basis.uninitialized & mask)) {
      return qmckl_failwith( context,
                             QMCKL_ALREADY_SET,
                             "qmckl_set_ao_*",
                             NULL);
   }

  if (ao_num <= 0) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_set_ao_basis_shell_num",
                           "ao_num must be positive");
  }

  const int64_t shell_num = ctx->ao_basis.shell_num;
  if (shell_num <= 0L) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_set_ao_basis_shell_num",
                           "shell_num is not set");
  }

  if (ao_num < shell_num) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_set_ao_basis_shell_num",
                           "ao_num < shell_num");
  }

  ctx->ao_basis.ao_num = ao_num;

  ctx->ao_basis.uninitialized &= ~mask;
  ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
  if (ctx->ao_basis.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_basis(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
   }
  
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_set_ao_basis_ao_factor (qmckl_context context,
                              const double* ao_factor,
                              const int64_t size_max)
{
  int32_t mask = 1 << 13;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
      return qmckl_failwith( context,
                             QMCKL_NULL_CONTEXT,
                             "qmckl_set_ao_*",
                             NULL);
  }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  
  if (mask != 0 && !(ctx->ao_basis.uninitialized & mask)) {
      return qmckl_failwith( context,
                             QMCKL_ALREADY_SET,
                             "qmckl_set_ao_*",
                             NULL);
   }

  const int64_t ao_num = ctx->ao_basis.ao_num;

  if (ao_num <= 0L) {
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "qmckl_set_ao_basis_ao_factor",
                           "ao_num is not set");
  }

  if (size_max < ao_num) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_set_ao_basis_ao_factor",
                           "input array too small");
  }

  if (ctx->ao_basis.ao_factor != NULL) {
    qmckl_exit_code rc = qmckl_free(context, ctx->ao_basis.ao_factor);
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context, rc,
                             "qmckl_set_ao_basis_ao_factor",
                             NULL);
    }
  }

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = ao_num * sizeof(double);
  double* new_array = (double*) qmckl_malloc(context, mem_info);

  if (new_array == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_set_ao_basis_ao_factor",
                           NULL);
  }

  memcpy(new_array, ao_factor, mem_info.size);

  ctx->ao_basis.ao_factor = new_array;

  ctx->ao_basis.uninitialized &= ~mask;
  ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
  if (ctx->ao_basis.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_basis(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
   }
  
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_set_ao_basis_r_power (qmckl_context context,
                            const int32_t* r_power,
                            const int64_t size_max)
{
//  Optional, so we don't need the mask here
  int32_t mask = 0;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
      return qmckl_failwith( context,
                             QMCKL_NULL_CONTEXT,
                             "qmckl_set_ao_*",
                             NULL);
  }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  
  if (mask != 0 && !(ctx->ao_basis.uninitialized & mask)) {
      return qmckl_failwith( context,
                             QMCKL_ALREADY_SET,
                             "qmckl_set_ao_*",
                             NULL);
   }

  const int64_t shell_num = ctx->ao_basis.shell_num;

  if (shell_num == 0L) {
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "qmckl_set_ao_basis_r_power",
                           "shell_num is not set");
  }

  if (size_max < shell_num) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_set_ao_basis_r_power",
                           "input array too small");
  }

  if (ctx->ao_basis.r_power != NULL) {
    qmckl_exit_code rc = qmckl_free(context, ctx->ao_basis.r_power);
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context, rc,
                             "qmckl_set_ao_basis_r_power",
                             NULL);
    }
  }


  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = shell_num * sizeof(int32_t);
  int32_t * new_array = (int32_t*) qmckl_malloc(context, mem_info);

  if (new_array == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_set_ao_basis_r_power",
                           NULL);
  }

  memcpy(new_array, r_power, mem_info.size);

  ctx->ao_basis.r_power = new_array;

  ctx->ao_basis.uninitialized &= ~mask;
  ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
  if (ctx->ao_basis.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_basis(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
   }
  
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_set_ao_basis_cartesian (qmckl_context context,
                              const bool cartesian)
{
  int32_t mask = 1;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
      return qmckl_failwith( context,
                             QMCKL_NULL_CONTEXT,
                             "qmckl_set_ao_*",
                             NULL);
  }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  
  if (mask != 0 && !(ctx->ao_basis.uninitialized & mask)) {
      return qmckl_failwith( context,
                             QMCKL_ALREADY_SET,
                             "qmckl_set_ao_*",
                             NULL);
   }

  ctx->ao_basis.ao_cartesian = cartesian;

  ctx->ao_basis.uninitialized &= ~mask;
  ctx->ao_basis.provided = (ctx->ao_basis.uninitialized == 0);
  if (ctx->ao_basis.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_basis(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
   }
  
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_ao_basis_type (const qmckl_context context,
                         char* const basis_type)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_ao_basis_type",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1;

  if ( (ctx->ao_basis.uninitialized & mask) != 0) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_get_ao_basis_type",
                           NULL);

  }

  assert (ctx->ao_basis.type != (char) 0);

  basis_type[0] = ctx->ao_basis.type;
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_ao_basis_shell_num (const qmckl_context context,
                              int64_t* const shell_num)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_ao_basis_shell_factor",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 1;

  if ( (ctx->ao_basis.uninitialized & mask) != 0) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_get_ao_basis_shell_num",
                           NULL);
  }

  assert (ctx->ao_basis.shell_num > (int64_t) 0);
  *shell_num = ctx->ao_basis.shell_num;
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_ao_basis_prim_num (const qmckl_context context,
                             int64_t* const prim_num)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_ao_basis_prim_num",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 2;

  if ( (ctx->ao_basis.uninitialized & mask) != 0) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_get_ao_basis_prim_num",
                           NULL);
  }

  assert (ctx->ao_basis.prim_num > (int64_t) 0);

  *prim_num = ctx->ao_basis.prim_num;
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_ao_basis_nucleus_shell_num (const qmckl_context context,
                                      int64_t* const nucleus_shell_num,
                                      const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_ao_basis_nucleus_shell_num",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 3;

  if ( (ctx->ao_basis.uninitialized & mask) != 0) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_get_ao_basis_nucleus_shell_num",
                           NULL);
  }

  if (nucleus_shell_num == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_ao_basis_nucleus_shell_num",
                           "NULL pointer");
  }

  if (size_max < ctx->nucleus.num) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_ao_basis_nucleus_shell_num",
                           "Array too small. Expected nucl_num");
  }

  assert (ctx->ao_basis.nucleus_shell_num != NULL);
  memcpy(nucleus_shell_num, ctx->ao_basis.nucleus_shell_num,
         (size_t) ctx->nucleus.num * sizeof(int64_t));
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_ao_basis_nucleus_index (const qmckl_context context,
                                  int64_t* const nucleus_index,
                                  const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_ao_basis_nucleus_index",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 4;

  if ( (ctx->ao_basis.uninitialized & mask) != 0) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_get_ao_basis_nucleus_index",
                           NULL);
  }

  if (nucleus_index == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_ao_basis_nucleus_index",
                           "NULL pointer");
  }

  if (size_max < ctx->nucleus.num) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_ao_basis_nucleus_index",
                           "Array too small. Expected shell_num");
  }

  assert (ctx->ao_basis.nucleus_index != NULL);
  memcpy(nucleus_index, ctx->ao_basis.nucleus_index,
         (size_t) ctx->nucleus.num * sizeof(int64_t));
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_ao_basis_r_power (const qmckl_context context,
                            int32_t* const r_power,
                            const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_ao_basis_r_power",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if ( (ctx->ao_basis.uninitialized) != 0) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_get_ao_basis_r_power",
                           NULL);
  }

  if (r_power == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_ao_basis_r_power",
                           "NULL pointer");
  }

  if (size_max < ctx->ao_basis.shell_num) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_ao_basis_r_power",
                           "Array too small. Expected shell_num");
  }

  assert (ctx->ao_basis.r_power != NULL);
  memcpy(r_power, ctx->ao_basis.r_power,
         (size_t) ctx->ao_basis.shell_num * sizeof(int32_t));
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_ao_basis_shell_ang_mom (const qmckl_context context,
                                  int32_t* const shell_ang_mom,
                                  const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_ao_basis_shell_ang_mom",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 5;

  if ( (ctx->ao_basis.uninitialized & mask) != 0) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_get_ao_basis_shell_ang_mom",
                           NULL);
  }

  if (shell_ang_mom == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_ao_basis_shell_ang_mom",
                           "NULL pointer");
  }

  if (size_max < ctx->ao_basis.shell_num) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_ao_basis_shell_ang_mom",
                           "Array too small. Expected shell_num");
  }

  assert (ctx->ao_basis.shell_ang_mom != NULL);
  memcpy(shell_ang_mom, ctx->ao_basis.shell_ang_mom,
         (size_t) ctx->ao_basis.shell_num * sizeof(int32_t));
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_ao_basis_shell_prim_num (const qmckl_context context,
                                   int64_t* const shell_prim_num,
                                   const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_ao_basis_shell_prim_num",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 6;

  if ( (ctx->ao_basis.uninitialized & mask) != 0) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_get_ao_basis_shell_prim_num",
                           NULL);
  }

  if (shell_prim_num == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_ao_basis_shell_prim_num",
                           "NULL pointer");
  }

  if (size_max < ctx->ao_basis.shell_num) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_ao_basis_shell_prim_num",
                           "Array too small. Expected shell_num");
  }

  assert (ctx->ao_basis.shell_prim_num != NULL);
  memcpy(shell_prim_num, ctx->ao_basis.shell_prim_num,
         (size_t) ctx->ao_basis.shell_num * sizeof(int64_t));
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_ao_basis_shell_prim_index (const qmckl_context context,
                                     int64_t* const shell_prim_index,
                                     const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_ao_basis_shell_prim_index",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 7;

  if ( (ctx->ao_basis.uninitialized & mask) != 0) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_get_ao_basis_shell_prim_index",
                           NULL);
  }

  if (shell_prim_index == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_ao_basis_shell_prim_index",
                           "NULL pointer");
  }

  if (size_max < ctx->ao_basis.shell_num) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_ao_basis_shell_prim_index",
                           "Array too small. Expected shell_num");
  }

  assert (ctx->ao_basis.shell_prim_index != NULL);
  memcpy(shell_prim_index, ctx->ao_basis.shell_prim_index,
         (size_t) ctx->ao_basis.shell_num * sizeof(int64_t));
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_ao_basis_shell_factor (const qmckl_context context,
                                 double* const shell_factor,
                                 const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_ao_basis_shell_factor",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 8;

  if ( (ctx->ao_basis.uninitialized & mask) != 0) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_get_ao_basis_shell_factor",
                           NULL);
  }

  if (shell_factor == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_ao_basis_shell_factor",
                           "NULL pointer");
  }

  if (size_max < ctx->ao_basis.shell_num) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_ao_basis_shell_factor",
                           "Array too small. Expected shell_num");
  }

  assert (ctx->ao_basis.shell_factor != NULL);
  memcpy(shell_factor, ctx->ao_basis.shell_factor,
         (size_t) ctx->ao_basis.shell_num * sizeof(double));
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_ao_basis_exponent (const qmckl_context context,
                             double* const exponent,
                             const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_ao_basis_exponent",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 9;

  if ( (ctx->ao_basis.uninitialized & mask) != 0) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_get_ao_basis_exponent",
                           NULL);
  }

  if (exponent == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_ao_basis_exponent",
                           "NULL pointer");
  }

  if (size_max < ctx->ao_basis.prim_num) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_ao_basis_exponent",
                           "Array too small. Expected prim_num");
  }

  assert (ctx->ao_basis.exponent != NULL);
  memcpy(exponent, ctx->ao_basis.exponent,
         (size_t) ctx->ao_basis.prim_num * sizeof(double));
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_ao_basis_coefficient (const qmckl_context context,
                                double* const coefficient,
                                const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_ao_basis_coefficient",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 10;

  if ( (ctx->ao_basis.uninitialized & mask) != 0) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_get_ao_basis_coefficient",
                           NULL);
  }

  if (coefficient == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_ao_basis_coefficient",
                           "NULL pointer");
  }

  if (size_max < ctx->ao_basis.prim_num) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_ao_basis_coefficient",
                           "Array too small. Expected prim_num");
  }
  assert (ctx->ao_basis.coefficient != NULL);
  memcpy(coefficient, ctx->ao_basis.coefficient,
         (size_t) ctx->ao_basis.prim_num * sizeof(double));
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_ao_basis_prim_factor (const qmckl_context context,
                                double* const prim_factor,
                                const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_ao_basis_prim_factor",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 11;

  if ( (ctx->ao_basis.uninitialized & mask) != 0) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_get_ao_basis_prim_factor",
                           NULL);
  }

  if (prim_factor == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_ao_basis_prim_factor",
                           "NULL pointer");
  }

  if (size_max < ctx->ao_basis.prim_num) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_ao_basis_prim_factor",
                           "Array too small. Expected prim_num");
  }

  assert (ctx->ao_basis.prim_factor != NULL);
  memcpy(prim_factor, ctx->ao_basis.prim_factor,
         (size_t) ctx->ao_basis.prim_num * sizeof(double));
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_ao_basis_ao_num (const qmckl_context context,
                           int64_t* const ao_num)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_ao_basis_ao_num",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 12;

  if ( (ctx->ao_basis.uninitialized & mask) != 0) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_get_ao_basis_ao_num",
                           NULL);
  }

  assert (ctx->ao_basis.ao_num > (int64_t) 0);

  *ao_num = ctx->ao_basis.ao_num;
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_ao_basis_ao_factor (const qmckl_context context,
                              double* const ao_factor,
                              const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_ao_basis_ao_factor",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 13;

  if ( (ctx->ao_basis.uninitialized & mask) != 0) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_get_ao_basis_ao_factor",
                           NULL);
  }

  if (ao_factor == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_ao_basis_ao_factor",
                           "NULL pointer");
  }

  if (size_max < ctx->ao_basis.ao_num) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_ao_basis_ao_factor",
                           "Array too small. Expected ao_num");
  }

  assert (ctx->ao_basis.ao_factor != NULL);
  memcpy(ao_factor, ctx->ao_basis.ao_factor, ctx->ao_basis.ao_num * sizeof(double));
  return QMCKL_SUCCESS;
}

bool qmckl_ao_basis_provided(const qmckl_context context) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return false;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  return ctx->ao_basis.provided;
}



/* Performs a deep copy of the AO basis structure from ~src~ to ~dest~. */
/* Memory allocations are done using the provided context. */


qmckl_exit_code qmckl_copy_ao_basis(qmckl_context context, const qmckl_ao_basis_struct* src, qmckl_ao_basis_struct* dest) {
  
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  }

  if (src == NULL || dest == NULL) {
    return qmckl_failwith(context,
                          QMCKL_INVALID_ARG_2,
                          "qmckl_copy_ao_basis",
                          "NULL pointer");
  }

  /* If source is not provided/initialized, just copy the struct */
  if (!src->provided) {
    *dest = *src;
    return QMCKL_SUCCESS;
  }

  /* Use memcpy to copy all scalar fields and pointer values at once */
  memcpy(dest, src, sizeof(qmckl_ao_basis_struct));

  /* Now update pointer fields with deep copies */
  
  /* Deep copy nucleus_index array */
  if (src->nucleus_index != NULL && src->shell_num > 0) {
    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = src->shell_num * sizeof(int64_t);
    dest->nucleus_index = (int64_t*) qmckl_malloc(context, mem_info);
    if (dest->nucleus_index == NULL) return QMCKL_ALLOCATION_FAILED;
    memcpy(dest->nucleus_index, src->nucleus_index, mem_info.size);
  }

  /* Deep copy nucleus_shell_num array */
  if (src->nucleus_shell_num != NULL && src->shell_num > 0) {
    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = src->shell_num * sizeof(int64_t);
    dest->nucleus_shell_num = (int64_t*) qmckl_malloc(context, mem_info);
    if (dest->nucleus_shell_num == NULL) return QMCKL_ALLOCATION_FAILED;
    memcpy(dest->nucleus_shell_num, src->nucleus_shell_num, mem_info.size);
  }

  /* Deep copy r_power array */
  if (src->r_power != NULL && src->shell_num > 0) {
    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = src->shell_num * sizeof(int32_t);
    dest->r_power = (int32_t*) qmckl_malloc(context, mem_info);
    if (dest->r_power == NULL) return QMCKL_ALLOCATION_FAILED;
    memcpy(dest->r_power, src->r_power, mem_info.size);
  }

  /* Deep copy shell_ang_mom array */
  if (src->shell_ang_mom != NULL && src->shell_num > 0) {
    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = src->shell_num * sizeof(int32_t);
    dest->shell_ang_mom = (int32_t*) qmckl_malloc(context, mem_info);
    if (dest->shell_ang_mom == NULL) return QMCKL_ALLOCATION_FAILED;
    memcpy(dest->shell_ang_mom, src->shell_ang_mom, mem_info.size);
  }

  /* Deep copy shell_prim_num array */
  if (src->shell_prim_num != NULL && src->shell_num > 0) {
    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = src->shell_num * sizeof(int64_t);
    dest->shell_prim_num = (int64_t*) qmckl_malloc(context, mem_info);
    if (dest->shell_prim_num == NULL) return QMCKL_ALLOCATION_FAILED;
    memcpy(dest->shell_prim_num, src->shell_prim_num, mem_info.size);
  }

  /* Deep copy shell_prim_index array */
  if (src->shell_prim_index != NULL && src->shell_num > 0) {
    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = src->shell_num * sizeof(int64_t);
    dest->shell_prim_index = (int64_t*) qmckl_malloc(context, mem_info);
    if (dest->shell_prim_index == NULL) return QMCKL_ALLOCATION_FAILED;
    memcpy(dest->shell_prim_index, src->shell_prim_index, mem_info.size);
  }

  /* Deep copy shell_factor array */
  if (src->shell_factor != NULL && src->shell_num > 0) {
    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = src->shell_num * sizeof(double);
    dest->shell_factor = (double*) qmckl_malloc(context, mem_info);
    if (dest->shell_factor == NULL) return QMCKL_ALLOCATION_FAILED;
    memcpy(dest->shell_factor, src->shell_factor, mem_info.size);
  }

  /* Deep copy exponent array */
  if (src->exponent != NULL && src->prim_num > 0) {
    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = src->prim_num * sizeof(double);
    dest->exponent = (double*) qmckl_malloc(context, mem_info);
    if (dest->exponent == NULL) return QMCKL_ALLOCATION_FAILED;
    memcpy(dest->exponent, src->exponent, mem_info.size);
  }

  /* Deep copy coefficient array */
  if (src->coefficient != NULL && src->prim_num > 0) {
    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = src->prim_num * sizeof(double);
    dest->coefficient = (double*) qmckl_malloc(context, mem_info);
    if (dest->coefficient == NULL) return QMCKL_ALLOCATION_FAILED;
    memcpy(dest->coefficient, src->coefficient, mem_info.size);
  }

  /* Deep copy prim_factor array */
  if (src->prim_factor != NULL && src->prim_num > 0) {
    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = src->prim_num * sizeof(double);
    dest->prim_factor = (double*) qmckl_malloc(context, mem_info);
    if (dest->prim_factor == NULL) return QMCKL_ALLOCATION_FAILED;
    memcpy(dest->prim_factor, src->prim_factor, mem_info.size);
  }

  /* Deep copy ao_factor array */
  if (src->ao_factor != NULL && src->ao_num > 0) {
    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = src->ao_num * sizeof(double);
    dest->ao_factor = (double*) qmckl_malloc(context, mem_info);
    if (dest->ao_factor == NULL) return QMCKL_ALLOCATION_FAILED;
    memcpy(dest->ao_factor, src->ao_factor, mem_info.size);
  }

  /* Deep copy ao_nucl array */
  if (src->ao_nucl != NULL && src->ao_num > 0) {
    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = src->ao_num * sizeof(int64_t);
    dest->ao_nucl = (int64_t*) qmckl_malloc(context, mem_info);
    if (dest->ao_nucl == NULL) return QMCKL_ALLOCATION_FAILED;
    memcpy(dest->ao_nucl, src->ao_nucl, mem_info.size);
  }

  /* Deep copy ao_ang_mom array */
  if (src->ao_ang_mom != NULL && src->ao_num > 0) {
    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = src->ao_num * sizeof(int32_t);
    dest->ao_ang_mom = (int32_t*) qmckl_malloc(context, mem_info);
    if (dest->ao_ang_mom == NULL) return QMCKL_ALLOCATION_FAILED;
    memcpy(dest->ao_ang_mom, src->ao_ang_mom, mem_info.size);
  }

  /* Deep copy nucleus_prim_index array (derived data) */
  if (src->nucleus_prim_index != NULL && src->shell_num > 0) {
    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = src->shell_num * sizeof(int64_t);
    dest->nucleus_prim_index = (int64_t*) qmckl_malloc(context, mem_info);
    if (dest->nucleus_prim_index == NULL) return QMCKL_ALLOCATION_FAILED;
    memcpy(dest->nucleus_prim_index, src->nucleus_prim_index, mem_info.size);
  }

  /* Deep copy coefficient_normalized array (derived data) */
  if (src->coefficient_normalized != NULL && src->prim_num > 0) {
    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = src->prim_num * sizeof(double);
    dest->coefficient_normalized = (double*) qmckl_malloc(context, mem_info);
    if (dest->coefficient_normalized == NULL) return QMCKL_ALLOCATION_FAILED;
    memcpy(dest->coefficient_normalized, src->coefficient_normalized, mem_info.size);
  }

  /* Deep copy nucleus_max_ang_mom array (derived data) */
  if (src->nucleus_max_ang_mom != NULL && src->shell_num > 0) {
    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = src->shell_num * sizeof(int32_t);
    dest->nucleus_max_ang_mom = (int32_t*) qmckl_malloc(context, mem_info);
    if (dest->nucleus_max_ang_mom == NULL) return QMCKL_ALLOCATION_FAILED;
    memcpy(dest->nucleus_max_ang_mom, src->nucleus_max_ang_mom, mem_info.size);
  }

  /* Deep copy nucleus_range array (derived data) */
  if (src->nucleus_range != NULL && src->shell_num > 0) {
    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = src->shell_num * sizeof(double);
    dest->nucleus_range = (double*) qmckl_malloc(context, mem_info);
    if (dest->nucleus_range == NULL) return QMCKL_ALLOCATION_FAILED;
    memcpy(dest->nucleus_range, src->nucleus_range, mem_info.size);
  }

  /* For computed arrays (primitive_vgl, shell_vgl, ao_vgl, ao_value, shell_hessian, ao_hessian),
     we set them to NULL as they can be recomputed from the basis data */
  /* Set dates to zero so they will be recomputed */
  dest->primitive_vgl = NULL;
  dest->primitive_vgl_date = 0;
  dest->shell_vgl = NULL;
  dest->shell_vgl_date = 0;
  dest->ao_vgl = NULL;
  dest->ao_vgl_date = 0;
  dest->ao_value = NULL;
  dest->ao_value_date = 0;
  dest->shell_hessian = NULL;
  dest->shell_hessian_date = 0;
  dest->ao_hessian = NULL;
  dest->ao_hessian_date = 0;

#ifdef HAVE_HPC
  /* Deep copy HPC-specific data structures */
  
  /* Deep copy prim_num_per_nucleus array */
  if (src->prim_num_per_nucleus != NULL && src->shell_num > 0) {
    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = src->shell_num * sizeof(int32_t);
    dest->prim_num_per_nucleus = (int32_t*) qmckl_malloc(context, mem_info);
    if (dest->prim_num_per_nucleus == NULL) return QMCKL_ALLOCATION_FAILED;
    memcpy(dest->prim_num_per_nucleus, src->prim_num_per_nucleus, mem_info.size);
  } else {
    dest->prim_num_per_nucleus = NULL;
  }

  /* Deep copy expo_per_nucleus matrix */
  if (src->expo_per_nucleus.data != NULL && src->expo_per_nucleus.size[0] > 0) {
    size_t size = src->expo_per_nucleus.size[0] * src->expo_per_nucleus.size[1] * sizeof(double);
    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = size;
    dest->expo_per_nucleus.data = (double*) qmckl_malloc(context, mem_info);
    if (dest->expo_per_nucleus.data == NULL) return QMCKL_ALLOCATION_FAILED;
    memcpy(dest->expo_per_nucleus.data, src->expo_per_nucleus.data, size);
    dest->expo_per_nucleus.size[0] = src->expo_per_nucleus.size[0];
    dest->expo_per_nucleus.size[1] = src->expo_per_nucleus.size[1];
  } else {
    dest->expo_per_nucleus.data = NULL;
    dest->expo_per_nucleus.size[0] = 0;
    dest->expo_per_nucleus.size[1] = 0;
  }

  /* Deep copy coef_per_nucleus tensor */
  if (src->coef_per_nucleus.data != NULL) {
    size_t total_size = 1;
    for (int i = 0; i < src->coef_per_nucleus.order; i++) {
      total_size *= src->coef_per_nucleus.size[i];
    }
    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = total_size * sizeof(double);
    dest->coef_per_nucleus.data = (double*) qmckl_malloc(context, mem_info);
    if (dest->coef_per_nucleus.data == NULL) return QMCKL_ALLOCATION_FAILED;
    memcpy(dest->coef_per_nucleus.data, src->coef_per_nucleus.data, mem_info.size);
    dest->coef_per_nucleus.order = src->coef_per_nucleus.order;
    for (int i = 0; i < src->coef_per_nucleus.order; i++) {
      dest->coef_per_nucleus.size[i] = src->coef_per_nucleus.size[i];
    }
  } else {
    dest->coef_per_nucleus.data = NULL;
    dest->coef_per_nucleus.order = 0;
  }
#endif

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_finalize_basis(qmckl_context context) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_finalize_basis",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t nucl_num = 0;

  qmckl_exit_code rc = qmckl_get_nucleus_num(context, &nucl_num);
  if (rc != QMCKL_SUCCESS) return rc;

  /* r_power */
  {
    if (ctx->ao_basis.r_power == NULL) {

      int32_t* r_power = (int32_t*) malloc( (ctx->ao_basis.shell_num) * sizeof(int32_t));

      if (r_power == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "ao_basis.r_power",
                               NULL);
      }

      for (int64_t i=0 ; i< ctx->ao_basis.shell_num ; ++i) {
        r_power[i] = 0;
      }

      rc = qmckl_set_ao_basis_r_power(context, r_power, ctx->ao_basis.shell_num);
      free(r_power);

      if (rc != QMCKL_SUCCESS) {
        return qmckl_failwith(context,
                             QMCKL_FAILURE,
                             "ao_basis.r_power",
                             NULL);
      }

    }
  }
  /* nucleus_prim_index */
  {
    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = (ctx->nucleus.num + (int64_t) 1) * sizeof(int64_t);

    ctx->ao_basis.nucleus_prim_index = (int64_t*) qmckl_malloc(context, mem_info);

    if (ctx->ao_basis.nucleus_prim_index == NULL) {
      return qmckl_failwith( context,
                             QMCKL_ALLOCATION_FAILED,
                             "ao_basis.nucleus_prim_index",
                             NULL);
    }

    for (int64_t i=0 ; i<nucl_num ; ++i) {
      int64_t shell_idx = ctx->ao_basis.nucleus_index[i];
      ctx->ao_basis.nucleus_prim_index[i] = ctx->ao_basis.shell_prim_index[shell_idx];
    }
    ctx->ao_basis.nucleus_prim_index[nucl_num] = ctx->ao_basis.prim_num;
  }

  /* ao_ang_mom */
  {
    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = ctx->ao_basis.ao_num * sizeof(int32_t);

    ctx->ao_basis.ao_ang_mom = (int32_t*) qmckl_malloc(context, mem_info);

    if (ctx->ao_basis.ao_ang_mom == NULL) {
      return qmckl_failwith( context,
                             QMCKL_ALLOCATION_FAILED,
                             "ao_basis.ao_ang_mom",
                             NULL);
    }

    mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = ctx->ao_basis.ao_num * sizeof(int64_t);
    ctx->ao_basis.ao_nucl = (int64_t*) qmckl_malloc(context, mem_info);

    if (ctx->ao_basis.ao_nucl == NULL) {
      return qmckl_failwith( context,
                             QMCKL_ALLOCATION_FAILED,
                             "ao_basis.ao_nucl",
                             NULL);
    }

    int64_t lstart[32];
    for (int32_t l=0 ; l<32 ; ++l) {
      lstart[l] = l*(l+1)*(l+2)/6;
    }

    int64_t ao_idx = 0;
    for (int64_t inucl=0 ; inucl<nucl_num ; ++inucl) {
      const int64_t ishell_start = ctx->ao_basis.nucleus_index[inucl];
      const int64_t ishell_end   = ctx->ao_basis.nucleus_index[inucl] + ctx->ao_basis.nucleus_shell_num[inucl];
      for (int64_t ishell = ishell_start ; ishell < ishell_end ; ++ishell) {
        const int l = ctx->ao_basis.shell_ang_mom[ishell];
        assert (l<32);
          for (int m=lstart[l] ; m < lstart[l+1]; m++) {
            ctx->ao_basis.ao_ang_mom[ao_idx] = l;
            ctx->ao_basis.ao_nucl[ao_idx] = inucl;
            ++ao_idx;
          }
      }
    }
    assert( ao_idx == ctx->ao_basis.ao_num );
  }

  /* Normalize coefficients */
  {
    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = ctx->ao_basis.prim_num * sizeof(double);

    ctx->ao_basis.coefficient_normalized = (double *) qmckl_malloc(context, mem_info);

    if (ctx->ao_basis.coefficient_normalized == NULL) {
      return qmckl_failwith( context,
                             QMCKL_ALLOCATION_FAILED,
                             "ao_basis.coefficient_normalized",
                             NULL);
    }

    for (int64_t ishell=0 ; ishell < ctx->ao_basis.shell_num ; ++ishell) {
      for (int64_t iprim=ctx->ao_basis.shell_prim_index[ishell] ;
           iprim < ctx->ao_basis.shell_prim_index[ishell]+ctx->ao_basis.shell_prim_num[ishell] ;
           ++iprim) {
        ctx->ao_basis.coefficient_normalized[iprim] =
          ctx->ao_basis.coefficient[iprim] * ctx->ao_basis.prim_factor[iprim] *
          ctx->ao_basis.shell_factor[ishell];
      }
    }
  }


  /* Find max angular momentum on each nucleus */
  {
    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = ctx->nucleus.num * sizeof(int32_t);

    ctx->ao_basis.nucleus_max_ang_mom = (int32_t *) qmckl_malloc(context, mem_info);

    if (ctx->ao_basis.nucleus_max_ang_mom == NULL) {
      return qmckl_failwith( context,
                             QMCKL_ALLOCATION_FAILED,
                             "ao_basis.nucleus_max_ang_mom",
                             NULL);
    }

    for (int64_t inucl=0 ; inucl < nucl_num ; ++inucl) {
      ctx->ao_basis.nucleus_max_ang_mom[inucl] = 0;
      for (int64_t ishell=ctx->ao_basis.nucleus_index[inucl]  ;
           ishell < ctx->ao_basis.nucleus_index[inucl] + ctx->ao_basis.nucleus_shell_num[inucl] ;
           ++ishell) {
        ctx->ao_basis.nucleus_max_ang_mom[inucl] =
          ctx->ao_basis.nucleus_max_ang_mom[inucl] > ctx->ao_basis.shell_ang_mom[ishell] ?
          ctx->ao_basis.nucleus_max_ang_mom[inucl] : ctx->ao_basis.shell_ang_mom[ishell] ;
      }
    }
  }

  /*
  // Find distance beyond which all AOs are zero using computed Gaussians.
  {
    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = ctx->nucleus.num * 53 * sizeof(double);

    ctx->ao_basis.nucleus_range = (double *) qmckl_malloc(context, mem_info);
    if (ctx->ao_basis.nucleus_range == NULL) {
      return qmckl_failwith( context,
                             QMCKL_ALLOCATION_FAILED,
                             "ao_basis.nucleus_range",
                             NULL);
    }
    for (int64_t i=0 ; i<ctx->nucleus.num * 53 ; ++i) {
      ctx->ao_basis.nucleus_range[i] = 50.;
    }
    if (ctx->ao_basis.type == 'G') {
      rc = qmckl_compute_nucleus_range_gaussian(context,
                                                ctx->ao_basis.ao_num,
                                                ctx->ao_basis.shell_num,
                                                ctx->ao_basis.prim_num,
                                                nucl_num,
                                                ctx->nucleus.coord.data,
                                                ctx->ao_basis.nucleus_index,
                                                ctx->ao_basis.nucleus_shell_num,
                                                ctx->ao_basis.nucleus_max_ang_mom,
                                                ctx->ao_basis.shell_prim_index,
                                                ctx->ao_basis.shell_prim_num,
                                                ctx->ao_basis.shell_ang_mom,
                                                ctx->ao_basis.ao_factor,
                                                ctx->ao_basis.exponent,
                                                ctx->ao_basis.coefficient_normalized,
                                                ctx->ao_basis.nucleus_range);
    } else if (ctx->ao_basis.type == 'S') {
      rc = qmckl_compute_nucleus_range_slater(context,
                                              ctx->ao_basis.ao_num,
                                              ctx->ao_basis.shell_num,
                                              ctx->ao_basis.prim_num,
                                              nucl_num,
                                              ctx->nucleus.coord.data,
                                              ctx->ao_basis.nucleus_index,
                                              ctx->ao_basis.nucleus_shell_num,
                                              ctx->ao_basis.nucleus_max_ang_mom,
                                              ctx->ao_basis.shell_prim_index,
                                              ctx->ao_basis.shell_prim_num,
                                              ctx->ao_basis.shell_ang_mom,
                                              ctx->ao_basis.ao_factor,
                                              ctx->ao_basis.exponent,
                                              ctx->ao_basis.coefficient_normalized,
                                              ctx->ao_basis.nucleus_range);
    }
  }
*/

  /* Find distance beyond which all AOs are zero.
     The distance is obtained by sqrt(-log(epsilon)*range) */
  {
    if (ctx->ao_basis.type == 'G' || ctx->ao_basis.type == 'S') {
      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->nucleus.num * sizeof(double);

      ctx->ao_basis.nucleus_range = (double *) qmckl_malloc(context, mem_info);

      if (ctx->ao_basis.nucleus_range == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "ao_basis.nucleus_range",
                               NULL);
      }

      for (int64_t inucl=0 ; inucl < ctx->nucleus.num ; ++inucl) {
        ctx->ao_basis.nucleus_range[inucl] = 0.;
        for (int64_t ishell=ctx->ao_basis.nucleus_index[inucl]  ;
             ishell < ctx->ao_basis.nucleus_index[inucl] + ctx->ao_basis.nucleus_shell_num[inucl] ;
             ++ishell) {
          for (int64_t iprim=ctx->ao_basis.shell_prim_index[ishell] ;
               iprim < ctx->ao_basis.shell_prim_index[ishell] + ctx->ao_basis.shell_prim_num[ishell] ;
               ++iprim) {
            double n_val = (double)(ctx->ao_basis.r_power[ishell] + 1);
            double range = n_val/ctx->ao_basis.exponent[iprim];
            ctx->ao_basis.nucleus_range[inucl] =
              ctx->ao_basis.nucleus_range[inucl] > range ?
              ctx->ao_basis.nucleus_range[inucl] : range;
          }
        }
      }
    }
  }

#ifdef HAVE_HPC
  rc = qmckl_finalize_basis_hpc(context);
  if (rc != QMCKL_SUCCESS) return rc;
#endif

  return qmckl_context_touch(context);
}

/* Data structure for sorting */
struct combined {
  double expo;
  int64_t index;
};

/* Comparison function */
int compare_basis( const void * l, const void * r )
{
  const struct combined * restrict _l= (const struct combined *)l;
  const struct combined * restrict _r= (const struct combined *)r;
  if( _l->expo > _r->expo ) return 1;
  if( _l->expo < _r->expo ) return -1;
  return 0;
}

#ifdef HAVE_HPC
qmckl_exit_code qmckl_finalize_basis_hpc (qmckl_context context)
{
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

  mem_info.size = ctx->nucleus.num * sizeof(int32_t);
  ctx->ao_basis.prim_num_per_nucleus = (int32_t*) qmckl_malloc(context, mem_info);

  /* Find max number of primitives per nucleus */
  int64_t shell_max = 0;
  int64_t prim_max  = 0;
  int64_t nucl_num = ctx->nucleus.num;
  for (int inucl=0 ; inucl < nucl_num ; ++inucl) {
    shell_max = ctx->ao_basis.nucleus_shell_num[inucl] > shell_max ?
      ctx->ao_basis.nucleus_shell_num[inucl] : shell_max;

    int64_t prim_num = 0;
    const int64_t ishell_start = ctx->ao_basis.nucleus_index[inucl];
    const int64_t ishell_end   = ctx->ao_basis.nucleus_index[inucl] + ctx->ao_basis.nucleus_shell_num[inucl];
    for (int64_t ishell = ishell_start ; ishell < ishell_end ; ++ishell) {
      prim_num += ctx->ao_basis.shell_prim_num[ishell];
    }
    prim_max = prim_num > prim_max ?
      prim_num : prim_max;
    ctx->ao_basis.prim_num_per_nucleus[inucl] = prim_num;
  }


  int64_t size[3] = { prim_max, shell_max, nucl_num };
  ctx->ao_basis.coef_per_nucleus = qmckl_tensor_alloc( context, 3, size );
  ctx->ao_basis.coef_per_nucleus = qmckl_tensor_set(ctx->ao_basis.coef_per_nucleus, 0.);

  ctx->ao_basis.expo_per_nucleus = qmckl_matrix_alloc( context, prim_max, nucl_num );
  ctx->ao_basis.expo_per_nucleus = qmckl_matrix_set(ctx->ao_basis.expo_per_nucleus, 0.);

  struct combined expo[prim_max];
  double coef[shell_max][prim_max];
  double newcoef[prim_max];
  for (int64_t inucl=0 ; inucl < nucl_num ; ++inucl) {
    memset(expo, 0, sizeof(expo));
    memset(coef, 0, sizeof(expo));

    int64_t idx = 0;
    const int64_t ishell_start = ctx->ao_basis.nucleus_index[inucl];
    const int64_t ishell_end   = ctx->ao_basis.nucleus_index[inucl] + ctx->ao_basis.nucleus_shell_num[inucl];
    for (int64_t ishell = ishell_start ; ishell < ishell_end ; ++ishell) {

      const int64_t iprim_start = ctx->ao_basis.shell_prim_index[ishell];
      const int64_t iprim_end   = ctx->ao_basis.shell_prim_index[ishell] + ctx->ao_basis.shell_prim_num[ishell];
      for (int64_t iprim = iprim_start ; iprim < iprim_end ; ++iprim) {
        expo[idx].expo = ctx->ao_basis.exponent[iprim];
        expo[idx].index = idx;
        idx += 1;
      }
    }

    /* Sort exponents */
    qsort( expo, (size_t) idx, sizeof(struct combined), compare_basis );

    idx = 0;
    int64_t idx2 = 0;
    for (int64_t ishell = ishell_start ; ishell < ishell_end ; ++ishell) {

      memset(newcoef, 0, sizeof(newcoef));
      const int64_t iprim_start = ctx->ao_basis.shell_prim_index[ishell];
      const int64_t iprim_end   = ctx->ao_basis.shell_prim_index[ishell] + ctx->ao_basis.shell_prim_num[ishell];
      for (int64_t iprim = iprim_start ; iprim < iprim_end ; ++iprim) {
        newcoef[idx] = ctx->ao_basis.coefficient_normalized[iprim];
        idx += 1;
      }
      for (int32_t i=0 ; i<ctx->ao_basis.prim_num_per_nucleus[inucl] ; ++i) {
        idx2 = expo[i].index;
        coef[ishell - ishell_start][i] = newcoef[idx2];
      }
    }

 /* Apply ordering to coefficients */

 /* Remove duplicates */
    int64_t newidx[prim_max];
    int64_t idxmax = 0;
    idx = 0;
    newidx[0] = 0;
    for (int32_t i=1 ; i<ctx->ao_basis.prim_num_per_nucleus[inucl] ; ++i) {
      if (expo[i].expo != expo[i-1].expo) {
        idx += 1;
      }
      newidx[i] = idx;
    }
    idxmax = idx;

    for (int32_t j=0 ; j<ishell_end-ishell_start ; ++j) {
      memset(newcoef, 0, sizeof(newcoef));
      for (int32_t i=0 ; i<ctx->ao_basis.prim_num_per_nucleus[inucl] ; ++i) {
        newcoef[newidx[i]] += coef[j][i];
      }
      for (int32_t i=0 ; i<ctx->ao_basis.prim_num_per_nucleus[inucl] ; ++i) {
        coef[j][i] = newcoef[i];
      }
    }

    for (int32_t i=0 ; i<ctx->ao_basis.prim_num_per_nucleus[inucl] ; ++i) {
        expo[newidx[i]].expo = expo[i].expo;
    }
    ctx->ao_basis.prim_num_per_nucleus[inucl] = idxmax+1;

    for (int32_t i=0 ; i<ctx->ao_basis.prim_num_per_nucleus[inucl] ; ++i) {
      qmckl_mat( ctx->ao_basis.expo_per_nucleus, i, inucl ) = expo[i].expo;
    }

    for (int32_t j=0 ; j<ishell_end-ishell_start ; ++j) {
      for (int32_t i=0 ; i<ctx->ao_basis.prim_num_per_nucleus[inucl] ; ++i) {
        qmckl_ten3( ctx->ao_basis.coef_per_nucleus, i, j, inucl ) = coef[j][i];
      }
    }

/*
    for (int32_t i=0 ; i<ctx->ao_basis.prim_num_per_nucleus[inucl] ; ++i) {
      printf("%4ld %4ld %15.5e | ", inucl, i, qmckl_mat( ctx->ao_basis.expo_per_nucleus, i, inucl ));
      for (int64_t j=0 ; j<ishell_end-ishell_start ; ++j) {
        printf("%8.5f ", qmckl_ten3( ctx->ao_basis.coef_per_nucleus, i, j, inucl ));
      }
      printf("\n");
    }
    printf("\n");
*/
  }


  return QMCKL_SUCCESS;
}
#endif



/* Returns the array of values, gradients and Laplacian of primitive */
/* basis functions evaluated at the current coordinates. */
/* See section [[Computation of primitives]]. */


qmckl_exit_code
qmckl_get_ao_basis_primitive_vgl (qmckl_context context,
                                  double* const primitive_vgl,
                                  const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_ao_basis_primitive_vgl",
                           NULL);
  }

  if (size_max <= 0) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_ao_basis_primitive_vgl",
                           "size_max <= 0");
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_ao_basis_primitive_vgl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t sze = ctx->ao_basis.prim_num * 5 * ctx->point.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_ao_basis_primitive_vgl",
                           "input array too small");
  }
  memcpy(primitive_vgl, ctx->ao_basis.primitive_vgl, (size_t) sze * sizeof(double));

  return QMCKL_SUCCESS;
}



/* Returns the array of values, gradients and Laplacian of contracted shells */
/* evaluated at the current coordinates. See section [[Computation of shells]]. */


qmckl_exit_code
qmckl_get_ao_basis_shell_vgl (qmckl_context context,
                              double* const shell_vgl,
                              const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_ao_basis_shell_vgl",
                           NULL);
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_ao_basis_shell_vgl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t sze = ctx->ao_basis.shell_num * 5 * ctx->point.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_ao_basis_shell_vgl",
                           "input array too small");
  }
  memcpy(shell_vgl, ctx->ao_basis.shell_vgl, (size_t)sze * sizeof(double));

  return QMCKL_SUCCESS;
}



/* Returns the array of values, gradients and Laplacian of the atomic orbitals */
/* evaluated at the current coordinates. */
/* See section [[Combining radial and polynomial parts]]. */


qmckl_exit_code
qmckl_get_ao_basis_ao_vgl (qmckl_context context,
                           double* const ao_vgl,
                           const int64_t size_max)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_ao_basis_ao_vgl",
                           NULL);
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_ao_basis_ao_vgl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t sze = ctx->ao_basis.ao_num * 5 * ctx->point.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_ao_basis_ao_vgl",
                           "input array too small");
  }
  memcpy(ao_vgl, ctx->ao_basis.ao_vgl, (size_t) sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_ao_basis_ao_vgl_inplace (qmckl_context context,
                                   double* const ao_vgl,
                                   const int64_t size_max)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_ao_basis_ao_vgl",
                           NULL);
  }

  qmckl_exit_code rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t sze = ctx->ao_basis.ao_num * 5 * ctx->point.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_ao_basis_ao_vgl",
                           "input array too small");
  }

  rc = qmckl_context_touch(context);
  if (rc != QMCKL_SUCCESS) return rc;

  double* old_array = ctx->ao_basis.ao_vgl;

  ctx->ao_basis.ao_vgl = ao_vgl;

  rc = qmckl_provide_ao_basis_ao_vgl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  ctx->ao_basis.ao_vgl = old_array;

  return QMCKL_SUCCESS;
}



/* Returns the array of values of the atomic orbitals evaluated at */
/* the current coordinates. See section [[Combining radial and polynomial parts]]. */


qmckl_exit_code
qmckl_get_ao_basis_ao_value (qmckl_context context,
                             double* const ao_value,
                             const int64_t size_max)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_ao_basis_ao_value",
                           NULL);
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_ao_basis_ao_value(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t sze = ctx->ao_basis.ao_num * ctx->point.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_ao_basis_ao_value",
                           "input array too small");
  }
  memcpy(ao_value, ctx->ao_basis.ao_value, (size_t) sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_ao_basis_ao_value_inplace (qmckl_context context,
                                     double* const ao_value,
                                     const int64_t size_max)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_ao_basis_ao_value",
                           NULL);
  }

  qmckl_exit_code rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t sze = ctx->ao_basis.ao_num * ctx->point.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_ao_basis_ao_value",
                           "input array too small");
  }

  rc = qmckl_context_touch(context);
  if (rc != QMCKL_SUCCESS) return rc;

  double* old_array = ctx->ao_basis.ao_value;

  ctx->ao_basis.ao_value = ao_value;

  rc = qmckl_provide_ao_basis_ao_value(context);
  if (rc != QMCKL_SUCCESS) return rc;

  ctx->ao_basis.ao_value = old_array;

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_ao_slater_vgl_hpc(const qmckl_context context,
                         const double *X,
                         const double *R,
                         const int64_t *num_slater,
                         const int64_t *N,
                         const double *A,
                         double *VGL,
                         const int64_t ldv)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  }

  if (*num_slater <= 0) {
    return QMCKL_INVALID_ARG_4;
  }

  if (ldv < *num_slater) {
    return QMCKL_INVALID_ARG_8;
  }

  const double Y[3] = {X[0] - R[0], X[1] - R[1], X[2] - R[2]};
  const double r = sqrt(Y[0]*Y[0] + Y[1]*Y[1] + Y[2]*Y[2]);

  // Handle r = 0 case
  if (r < 1.e-15) {
    for (int64_t i = 0; i < *num_slater; ++i) {
      const int64_t idx = i;
      VGL[idx] = (N[i] == 1) ? 1.0 : 0.0;
      VGL[ldv + idx] = 0.0;
      VGL[2*ldv + idx] = 0.0;
      VGL[3*ldv + idx] = 0.0;
      VGL[4*ldv + idx] = 0.0;
    }
    return QMCKL_SUCCESS;
  }

  const double r_inv = 1.0 / r;
  const double r_inv_2 = r_inv * r_inv;
  const double x_over_r = Y[0] * r_inv;
  const double y_over_r = Y[1] * r_inv;
  const double z_over_r = Y[2] * r_inv;

  // Vectorized computation for HPC performance
  #pragma omp simd
  for (int64_t i = 0; i < *num_slater; ++i) {
    const double alpha = A[i];
    const double n_val = (double) N[i];

    // Compute r^n * exp(-alpha * r) efficiently
    const double rn = pow(r, N[i]);
    const double exp_alpha_r = exp(-alpha * r);

    // Value: r^n * exp(-alpha * r)
    const double val = rn * exp_alpha_r;
    VGL[i] = val;

    // Gradient coefficients
    const double grad_coef = (n_val * r_inv - alpha) * val;
    VGL[ldv + i] = grad_coef * x_over_r;
    VGL[2*ldv + i] = grad_coef * y_over_r;
    VGL[3*ldv + i] = grad_coef * z_over_r;

    // Laplacian coefficient
    const double lapl_coef = n_val * (n_val - 1.0) * r_inv_2 - 2.0 * n_val * alpha * r_inv + alpha * alpha;
    VGL[4*ldv + i] = lapl_coef * val;
  }

  return QMCKL_SUCCESS;
}



/* #+CALL: write_provider_pre( group="ao_basis", data="primitive_vgl", dimension="ctx->ao_basis.prim_num * 5 * ctx->point.num") */

/* #+RESULTS: */

qmckl_exit_code qmckl_provide_ao_basis_primitive_vgl(qmckl_context context)
{

  qmckl_exit_code rc = QMCKL_SUCCESS;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_provide_ao_basis_primitive_vgl",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->ao_basis.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_provide_ao_basis_primitive_vgl",
                           NULL);
  }

  /* Compute if necessary */
  if (ctx->point.date > ctx->ao_basis.primitive_vgl_date) {

    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = ctx->ao_basis.prim_num * 5 * ctx->point.num * sizeof(double);

    if (ctx->ao_basis.primitive_vgl != NULL) {
      qmckl_memory_info_struct mem_info_test = qmckl_memory_info_struct_zero;
      rc = qmckl_get_malloc_info(context, ctx->ao_basis.primitive_vgl, &mem_info_test);

      /* if rc != QMCKL_SUCCESS, we are maybe in an _inplace function because the
         memory was not allocated with qmckl_malloc */

      if ((rc == QMCKL_SUCCESS) && (mem_info_test.size != mem_info.size)) {
        rc = qmckl_free(context, ctx->ao_basis.primitive_vgl);
        assert (rc == QMCKL_SUCCESS);
        ctx->ao_basis.primitive_vgl = NULL;
      }
    }

    /* Allocate array */
    if (ctx->ao_basis.primitive_vgl == NULL) {

      double* primitive_vgl = (double*) qmckl_malloc(context, mem_info);

      if (primitive_vgl == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_ao_basis_primitive_vgl",
                               NULL);
      }
      ctx->ao_basis.primitive_vgl = primitive_vgl;
    }

if (ctx->ao_basis.type == 'G') {
  rc = qmckl_compute_ao_basis_primitive_gaussian_vgl(context,
                                                     ctx->ao_basis.prim_num,
                                                     ctx->point.num,
                                                     ctx->nucleus.num,
                                                     ctx->ao_basis.nucleus_prim_index,
                                                     ctx->point.coord.data,
                                                     ctx->nucleus.coord.data,
                                                     ctx->ao_basis.exponent,
                                                     ctx->ao_basis.primitive_vgl);
} else {
  return qmckl_failwith( context,
                         QMCKL_FAILURE,
                         "compute_ao_basis_primitive_vgl",
                         "Not yet implemented");
}



/* #+CALL: write_provider_post( group="ao_basis", data="shell_vgl" ) */

/* #+RESULTS: */

if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->ao_basis.shell_vgl_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}



/* #+CALL: write_provider_pre( group="ao_basis", data="shell_vgl", dimension="ctx->ao_basis.shell_num * 5 * ctx->point.num") */

/* #+RESULTS: */

qmckl_exit_code qmckl_provide_ao_basis_shell_vgl(qmckl_context context)
{

  qmckl_exit_code rc = QMCKL_SUCCESS;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_provide_ao_basis_shell_vgl",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->ao_basis.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_provide_ao_basis_shell_vgl",
                           NULL);
  }

  /* Compute if necessary */
  if (ctx->point.date > ctx->ao_basis.shell_vgl_date) {

    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = ctx->ao_basis.shell_num * 5 * ctx->point.num * sizeof(double);

    if (mem_info.size > ctx->ao_basis.shell_vgl_maxsize) {
      if (ctx->ao_basis.shell_vgl != NULL) {
          rc = qmckl_free(context, ctx->ao_basis.shell_vgl);
          assert (rc == QMCKL_SUCCESS);
          ctx->ao_basis.shell_vgl = NULL;
      }
    }

    /* Allocate array */
    if (ctx->ao_basis.shell_vgl == NULL) {

      double* shell_vgl = (double*) qmckl_malloc(context, mem_info);
      ctx->ao_basis.shell_vgl_maxsize = mem_info.size;
      if (shell_vgl == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_ao_basis_shell_vgl",
                               NULL);
      }
      ctx->ao_basis.shell_vgl = shell_vgl;
    }

if (ctx->ao_basis.type == 'G') {
  rc = qmckl_compute_ao_basis_shell_gaussian_vgl(context,
                                                 ctx->ao_basis.prim_num,
                                                 ctx->ao_basis.shell_num,
                                                 ctx->point.num,
                                                 ctx->nucleus.num,
                                                 ctx->ao_basis.nucleus_shell_num,
                                                 ctx->ao_basis.nucleus_index,
                                                 ctx->ao_basis.nucleus_range,
                                                 ctx->ao_basis.shell_prim_index,
                                                 ctx->ao_basis.shell_prim_num,
                                                 ctx->point.coord.data,
                                                 ctx->nucleus.coord.data,
                                                 ctx->ao_basis.exponent,
                                                 ctx->ao_basis.coefficient_normalized,
                                                 ctx->ao_basis.shell_vgl);
} else if (ctx->ao_basis.type == 'S') {
  rc = qmckl_compute_ao_basis_shell_slater_vgl(context,
                                               ctx->ao_basis.prim_num,
                                               ctx->ao_basis.shell_num,
                                               ctx->point.num,
                                               ctx->nucleus.num,
                                               ctx->ao_basis.nucleus_shell_num,
                                               ctx->ao_basis.nucleus_index,
                                               ctx->ao_basis.nucleus_range,
                                               ctx->ao_basis.shell_prim_index,
                                               ctx->ao_basis.shell_prim_num,
                                               ctx->ao_basis.r_power,
                                               ctx->point.coord.data,
                                               ctx->nucleus.coord.data,
                                               ctx->ao_basis.exponent,
                                               ctx->ao_basis.coefficient_normalized,
                                               ctx->ao_basis.shell_vgl);
} else {
  return qmckl_failwith( context,
                         QMCKL_FAILURE,
                         "compute_ao_basis_shell_vgl",
                         "Not yet implemented");
}


/* #+CALL: write_provider_post( group="ao_basis", data="shell_vgl" ) */

/* #+RESULTS: */

if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->ao_basis.shell_vgl_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_ao_basis_shell_hessian (qmckl_context context,
                              double* const shell_hessian,
                              const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_ao_basis_shell_hessian",
                           NULL);
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_ao_basis_shell_hessian(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t sze = ctx->ao_basis.shell_num * 4 * 3 * ctx->point.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_ao_basis_shell_hessian",
                           "input array too small");
  }
  memcpy(shell_hessian, ctx->ao_basis.shell_hessian, (size_t)sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_ao_basis_shell_hessian(qmckl_context context)
{

  qmckl_exit_code rc = QMCKL_SUCCESS;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_provide_ao_basis_shell_hessian",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->ao_basis.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_provide_ao_basis_shell_hessian",
                           NULL);
  }

  /* Compute if necessary */
  if (ctx->point.date > ctx->ao_basis.shell_hessian_date) {

    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = ctx->ao_basis.shell_num * 4 * 3 * ctx->point.num * sizeof(double);

    if (ctx->ao_basis.shell_hessian != NULL) {
      qmckl_memory_info_struct mem_info_test = qmckl_memory_info_struct_zero;
      rc = qmckl_get_malloc_info(context, ctx->ao_basis.shell_hessian, &mem_info_test);

      if ((rc == QMCKL_SUCCESS) && (mem_info_test.size != mem_info.size)) {
        rc = qmckl_free(context, ctx->ao_basis.shell_hessian);
        assert (rc == QMCKL_SUCCESS);
        ctx->ao_basis.shell_hessian = NULL;
      }
    }

    /* Allocate array */
    if (ctx->ao_basis.shell_hessian == NULL) {

      double* shell_hessian = (double*) qmckl_malloc(context, mem_info);

      if (shell_hessian == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_ao_basis_shell_hessian",
                               NULL);
      }
      ctx->ao_basis.shell_hessian = shell_hessian;
    }

    if (ctx->ao_basis.type == 'G') {
      rc = qmckl_compute_ao_basis_shell_gaussian_hessian(context,
                                                     ctx->ao_basis.prim_num,
                                                     ctx->ao_basis.shell_num,
                                                     ctx->point.num,
                                                     ctx->nucleus.num,
                                                     ctx->ao_basis.nucleus_shell_num,
                                                     ctx->ao_basis.nucleus_index,
                                                     ctx->ao_basis.nucleus_range,
                                                     ctx->ao_basis.shell_prim_index,
                                                     ctx->ao_basis.shell_prim_num,
                                                     ctx->point.coord.data,
                                                     ctx->nucleus.coord.data,
                                                     ctx->ao_basis.exponent,
                                                     ctx->ao_basis.coefficient_normalized,
                                                     ctx->ao_basis.shell_hessian);

    } else if (ctx->ao_basis.type == 'S') {
      rc = qmckl_compute_ao_basis_shell_slater_hessian(context,
                                                       ctx->ao_basis.prim_num,
                                                       ctx->ao_basis.shell_num,
                                                       ctx->point.num,
                                                       ctx->nucleus.num,
                                                       ctx->ao_basis.nucleus_shell_num,
                                                       ctx->ao_basis.nucleus_index,
                                                       ctx->ao_basis.nucleus_range,
                                                       ctx->ao_basis.shell_prim_index,
                                                       ctx->ao_basis.shell_prim_num,
                                                       ctx->ao_basis.shell_ang_mom,
                                                       ctx->point.coord.data,
                                                       ctx->nucleus.coord.data,
                                                       ctx->ao_basis.exponent,
                                                       ctx->ao_basis.coefficient_normalized,
                                                       ctx->ao_basis.shell_hessian);
    } else {
      return qmckl_failwith( context,
                             QMCKL_FAILURE,
                             "compute_ao_basis_shell_hessian",
                             "Not yet implemented");
    }

    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->ao_basis.shell_hessian_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_ao_polynomial_vgl (const qmckl_context context,
                         const double* X,
                         const double* R,
                         const int32_t lmax,
                         int64_t* n,
                         int32_t* const L,
                         const int64_t ldl,
                         double* const VGL,
                         const int64_t ldv )
{
#ifdef HAVE_HPC
  //return qmckl_ao_polynomial_vgl_hpc (context, X, R, lmax, n, L, ldl, VGL, ldv);
  return qmckl_ao_polynomial_vgl_doc
#else
  return qmckl_ao_polynomial_vgl_doc
#endif
    (context, X, R, lmax, n, L, ldl, VGL, ldv);
}

qmckl_exit_code
qmckl_ao_polynomial_transp_vgl (const qmckl_context context,
                                const double* X,
                                const double* R,
                                const int32_t lmax,
                                int64_t* n,
                                int32_t* const L,
                                const int64_t ldl,
                                double* const VGL,
                                const int64_t ldv )
{
#ifdef HAVE_HPC
  return qmckl_ao_polynomial_transp_vgl_hpc
#else
  return qmckl_ao_polynomial_transp_vgl_doc
#endif
    (context, X, R, lmax, n, L, ldl, VGL, ldv);
}

static inline qmckl_exit_code
qmckl_ao_polynomial_transp_vgl_hpc_inline (const qmckl_context context,
                                    const double* restrict X,
                                    const double* restrict R,
                                    const int32_t lmax,
                                    int64_t* n,
                                    int32_t* restrict const L,
                                    const int64_t ldl,
                                    double* restrict const VGL,
                                    const int64_t ldv )
{
  const qmckl_context_struct* ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL && X != NULL && R != NULL && n != NULL && L != NULL && VGL != NULL);
  if (lmax < 0) return QMCKL_INVALID_ARG_4;
  if (ldl < 3) return QMCKL_INVALID_ARG_7;

  int32_t m;
  const double shift=1.e-20;

  switch (lmax) {
  case 0:
    {
      if (ldv < 1) return QMCKL_INVALID_ARG_9;
      L[0] = 0; L[1] = 0; L[2] = 0;

      VGL[0         ] = 1.0;
      VGL[ldv       ] = 0.0;
      VGL[ldv<<1    ] = 0.0;
      VGL[(ldv<<1)+ldv] = 0.0;
      VGL[ldv<<2    ] = 0.0;
      m=1;
      break;
    }
  case 1:
    {
      if (ldv < 4) return QMCKL_INVALID_ARG_9;

      if (ldl == 3) {

        const int32_t ltmp[12] = {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1};
        for (int i=0 ; i<12 ; ++i)
          L[i] = ltmp[i];

      } else {

        int32_t* restrict const l[4] = {L, L+ldl, L+(ldl<<1), L+ldl+(ldl<<1)};
        l[0][0] = 0; l[0][1] = 0; l[0][2] = 0;
        l[1][0] = 1; l[1][1] = 0; l[1][2] = 0;
        l[2][0] = 0; l[2][1] = 1; l[2][2] = 0;
        l[3][0] = 0; l[3][1] = 0; l[3][2] = 1;
      }


      if (ldv == 4) {

        const double tmp[20] = {1.0, (X[0]-R[0])+shift, (X[1]-R[1])+shift, (X[2]-R[2])+shift,
                                0.0, 1.0, 0.0, 0.0,
                                0.0, 0.0, 1.0, 0.0,
                                0.0, 0.0, 0.0, 1.0,
                                0.0, 0.0, 0.0, 0.0};

        for (int i=0 ; i<20 ; ++i)
          VGL[i] = tmp[i];

      } else {

        double* restrict const vgl1 = VGL;
        double* restrict const vgl2 = VGL + ldv;
        double* restrict const vgl3 = VGL + (ldv << 1);
        double* restrict const vgl4 = VGL + ldv + (ldv << 1);
        double* restrict const vgl5 = VGL + (ldv << 2);

        for (int32_t k=0 ; k<4 ; ++k) {
          vgl2[k] = 0.0;
          vgl3[k] = 0.0;
          vgl4[k] = 0.0;
          vgl5[k] = 0.0;
        }
        vgl1[0] = 1.0;
        vgl1[1] = (X[0]-R[0])+shift;
        vgl1[2] = (X[1]-R[1])+shift;
        vgl1[3] = (X[2]-R[2])+shift;
        vgl2[1] = 1.0;
        vgl3[2] = 1.0;
        vgl4[3] = 1.0;
      }
      m=4;
      break;
    }
  case 2:
    {
      if (ldv < 10) return QMCKL_INVALID_ARG_9;
      if (ldl == 3) {
        const int32_t ltmp[30] = {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1,
                                  2, 0, 0, 1, 1, 0, 1, 0, 1, 0, 2, 0,
                                  0, 1, 1, 0, 0, 2};
        for (int i=0 ; i<30 ; ++i)
          L[i] = ltmp[i];

      } else {
        int32_t* restrict l[10];
        for (int32_t i=0 ; i<10 ; ++i) {
          l[i] = L + i*ldl;
        }

        l[0][0] = 0; l[0][1] = 0; l[0][2] = 0;
        l[1][0] = 1; l[1][1] = 0; l[1][2] = 0;
        l[2][0] = 0; l[2][1] = 1; l[2][2] = 0;
        l[3][0] = 0; l[3][1] = 0; l[3][2] = 1;
        l[4][0] = 2; l[4][1] = 0; l[4][2] = 0;
        l[5][0] = 1; l[5][1] = 1; l[5][2] = 0;
        l[6][0] = 1; l[6][1] = 0; l[6][2] = 1;
        l[7][0] = 0; l[7][1] = 2; l[7][2] = 0;
        l[8][0] = 0; l[8][1] = 1; l[8][2] = 1;
        l[9][0] = 0; l[9][1] = 0; l[9][2] = 2;
      }

      const double Y[3] = { X[0]-R[0],
                            X[1]-R[1],
                            X[2]-R[2] };

      if (ldv == 50) {
         const double tmp[50] = {
            1.0, Y[0]+shift, Y[1]+shift, Y[2]+shift,
            Y[0] * Y[0]+shift, Y[0] * Y[1]+shift, Y[0] * Y[2]+shift,
            Y[1] * Y[1]+shift, Y[1] * Y[2]+shift, Y[2] * Y[2]+shift,
            0.0, 1.0, 0.0, 0.0,
            Y[0] + Y[0], Y[1], Y[2], 0.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0, 0.0,
            Y[0], 0.0, Y[1] + Y[1], Y[2], 0.0,
            0.0, 0.0, 0.0, 1.0, 0.0,
            0.0, Y[0], 0.0, Y[1], Y[2] + Y[2],
            0.0, 0.0, 0.0, 0.0, 2.0,
            0.0, 0.0, 2.0, 0., 2.0 };
         for (int i=0 ; i<50 ; ++i)
           VGL[i] = tmp[i];
      } else {
        double* restrict const vgl1 = VGL;
        double* restrict const vgl2 = VGL + ldv;
        double* restrict const vgl3 = VGL + (ldv << 1);
        double* restrict const vgl4 = VGL + 3*ldv;
        double* restrict const vgl5 = VGL + (ldv << 2);

        vgl1[0] = 1.0            ; vgl1[1] = Y[0]+shift     ; vgl1[2] = Y[1]+shift;
        vgl1[3] = Y[2]+shift     ; vgl1[4] = Y[0]*Y[0]+shift; vgl1[5] = Y[0]*Y[1]+shift;
        vgl1[6] = Y[0]*Y[2]+shift; vgl1[7] = Y[1]*Y[1]+shift; vgl1[8] = Y[1]*Y[2]+shift;
        vgl1[9] = Y[2]*Y[2]+shift;

        vgl2[0] = 0.0 ; vgl2[1] = 1.0      ; vgl2[2] = 0.0 ;
        vgl2[3] = 0.0 ; vgl2[4] = Y[0]+Y[0]; vgl2[5] = Y[1];
        vgl2[6] = Y[2]; vgl2[7] = 0.0      ; vgl2[8] = 0.0 ;
        vgl2[9] = 0.0 ;

        vgl3[0] = 0.0; vgl3[1] = 0.0      ; vgl3[2] = 1.0 ;
        vgl3[3] = 0.0; vgl3[4] = 0.0      ; vgl3[5] = Y[0];
        vgl3[6] = 0.0; vgl3[7] = Y[1]+Y[1]; vgl3[8] = Y[2];
        vgl3[9] = 0.0;

        vgl4[0] = 0.0      ; vgl4[1] = 0.0; vgl4[2] = 0.0 ;
        vgl4[3] = 1.0      ; vgl4[4] = 0.0; vgl4[5] = 0.0 ;
        vgl4[6] = Y[0]     ; vgl4[7] = 0.0; vgl4[8] = Y[1];
        vgl4[9] = Y[2]+Y[2];

        vgl5[0] = 0.0; vgl5[1] = 0.0; vgl5[2] = 0.0;
        vgl5[3] = 0.0; vgl5[4] = 2.0; vgl5[5] = 0.0;
        vgl5[6] = 0.0; vgl5[7] = 2.0; vgl5[8] = 0.0;
        vgl5[9] = 2.0;
      }
      m=10;
      break;
    }
  default:
    {
      const int32_t size_max = (lmax+1)*(lmax+2)*(lmax+3)/6;
      if (ldv < size_max) return QMCKL_INVALID_ARG_9;

      double* restrict const vgl1 = VGL;
      double* restrict const vgl2 = VGL + ldv;
      double* restrict const vgl3 = VGL + (ldv<<1);
      double* restrict const vgl4 = VGL + ldv + (ldv<<1);
      double* restrict const vgl5 = VGL + (ldv<<2);

      const double Y[3] = { X[0]-R[0],
                            X[1]-R[1],
                            X[2]-R[2]};

      assert(size_max > lmax+3);
      double pows[3][size_max];

      for (int32_t i=0 ; i<3 ; ++i) {
        pows[0][i] = 1.0;
        pows[1][i] = 1.0;
        pows[2][i] = 1.0;
      }

      for (int32_t i=3 ; i<=lmax+2 ; ++i) {
        pows[0][i] = pows[0][i-1] * Y[0];
        pows[1][i] = pows[1][i-1] * Y[1];
        pows[2][i] = pows[2][i-1] * Y[2];
      }

      int32_t* l[size_max];
      for (int32_t i=0 ; i<size_max ; ++i) {
        l[i] = &(L[i*ldl]);
      }

      for (int32_t i=0 ; i<4 ; ++i) {
        l[i][0] = 0;
        l[i][1] = 0;
        l[i][2] = 0;
      }
      l[1][0] = 1;
      l[2][1] = 1;
      l[3][2] = 1;

      for (int32_t k=0 ; k<4 ; ++k) {
        vgl2[k] = 0.0;
        vgl3[k] = 0.0;
        vgl4[k] = 0.0;
        vgl5[k] = 0.0;
      }
      vgl1[0] = 1.0;
      vgl1[1] = Y[0]+shift;
      vgl1[2] = Y[1]+shift;
      vgl1[3] = Y[2]+shift;

      vgl2[1] = 1.0;
      vgl3[2] = 1.0;
      vgl4[3] = 1.0;
      m=4;

      double dd = 2.0;

      for (int32_t d=2 ; d<= lmax ; ++d) {
        double da = dd;

        for (int32_t a=d ; a>=0 ; --a) {
          double db = dd-da;

          for (int32_t b=d-a ; b>=0 ; --b) {
            const int32_t c  = d - a - b;
            const double dc = dd - da - db;

            double xy = pows[0][a+2] * pows[1][b+2];
            double yz = pows[1][b+2] * pows[2][c+2];
            double xz = pows[0][a+2] * pows[2][c+2];

            l[m][0] = a;
            l[m][1] = b;
            l[m][2] = c;

            vgl1[m] = xy * pows[2][c+2]+shift;

            xy *= dc;
            xz *= db;
            yz *= da;

            vgl2[m] = pows[0][a+1] * yz;
            vgl3[m] = pows[1][b+1] * xz;
            vgl4[m] = pows[2][c+1] * xy;

            vgl5[m] = (da-1.) * pows[0][a] * yz +
              (db-1.) * pows[1][b] * xz +
              (dc-1.) * pows[2][c] * xy;

            db -= 1.0;
            ++m;
          }
          da -= 1.0;
        }
        dd += 1.0;
      }
    }
  }
  *n = m;
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_ao_polynomial_transp_vgl_hpc (const qmckl_context context,
                                    const double* restrict X,
                                    const double* restrict R,
                                    const int32_t lmax,
                                    int64_t* n,
                                    int32_t* restrict const L,
                                    const int64_t ldl,
                                    double* restrict const VGL,
                                    const int64_t ldv )
{
  return qmckl_ao_polynomial_transp_vgl_hpc_inline (context,
                                    X, R, lmax, n, L, ldl, VGL, ldv );
  }

qmckl_exit_code
qmckl_ao_polynomial_hessian (const qmckl_context context,
                         const double* X,
                         const double* R,
                         const int32_t lmax,
                         int64_t* n,
                         double* const hessian)
{
#ifdef HAVE_HPC
  return qmckl_compute_ao_polynomial_hessian_doc
#else
  return qmckl_compute_ao_polynomial_hessian_doc
#endif
    (context, X, R, lmax, n, hessian);
}

#ifdef HAVE_HPC
qmckl_exit_code
qmckl_compute_ao_value_hpc_gaussian (const qmckl_context context,
                                     const int64_t ao_num,
                                     const int64_t shell_num,
                                     const int32_t* restrict prim_num_per_nucleus,
                                     const int64_t point_num,
                                     const int64_t nucl_num,
                                     const double* restrict coord,
                                     const double* restrict nucl_coord,
                                     const int64_t* restrict nucleus_index,
                                     const int64_t* restrict nucleus_shell_num,
                                     const double* nucleus_range,
                                     const int32_t* restrict nucleus_max_ang_mom,
                                     const int32_t* restrict shell_ang_mom,
                                     const double* restrict ao_factor,
                                     const qmckl_matrix expo_per_nucleus,
                                     const qmckl_tensor coef_per_nucleus,
                                     double* restrict const ao_value )
{
  int32_t lstart[32] __attribute__((aligned(64)));
  for (int32_t l=0 ; l<32 ; ++l) {
    lstart[l] = l*(l+1)*(l+2)/6;
  }
  
  int64_t ao_index[shell_num+1] __attribute__((aligned(64)));
  int64_t size_max = 0;
  int64_t prim_max = 0;
  int64_t shell_max = 0;
  {
    int64_t k=0;
    for (int inucl=0 ; inucl < nucl_num ; ++inucl) {
      prim_max = prim_num_per_nucleus[inucl] > prim_max ?
        prim_num_per_nucleus[inucl] : prim_max;
      shell_max = nucleus_shell_num[inucl] > shell_max ?
        nucleus_shell_num[inucl] : shell_max;
      const int64_t ishell_start = nucleus_index[inucl];
      const int64_t ishell_end   = nucleus_index[inucl] + nucleus_shell_num[inucl];
      for (int64_t ishell = ishell_start ; ishell < ishell_end ; ++ishell) {
        const int l = shell_ang_mom[ishell];
        ao_index[ishell] = k;
        k += lstart[l+1] - lstart[l];
        size_max = size_max < lstart[l+1] ? lstart[l+1] : size_max;
      }
    }
    ao_index[shell_num] = ao_num+1;
  }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);
  
  /* Don't compute polynomials when the radial part is zero. */
  const int precision = ctx->numprec.precision;
  const double cutoff = log( (double) ( ((uint64_t) 1) << (precision-1) ) );
  const bool double_precision = precision > 23;
  

#ifdef HAVE_OPENMP
#pragma omp parallel if (point_num > 16)
#endif
  {
    qmckl_exit_code rc;
    double ar2[prim_max] __attribute__((aligned(64)));
    float  ar2_sp[prim_max] __attribute__((aligned(64)));
    int32_t powers[3*size_max] __attribute__((aligned(64)));
    double  poly_vgl[5*size_max] __attribute__((aligned(64)));

    double exp_mat[prim_max] __attribute__((aligned(64)));
    float  exp_mat_sp[prim_max] __attribute__((aligned(64)));
    double ce_mat[shell_max] __attribute__((aligned(64)));

    int32_t (*coef_mat_sparse_idx)[shell_max][prim_max+1] __attribute__((aligned(64))) =
      malloc(nucl_num * sizeof (*coef_mat_sparse_idx));

    double  (*coef_mat_sparse)[shell_max][prim_max+1] __attribute__((aligned(64))) =
      malloc(nucl_num * sizeof (*coef_mat_sparse));

    for (int i=0 ; i<nucl_num ; ++i) {
      for (int j=0 ; j<shell_max; ++j) {
        int l=1;
        for (int k=0 ; k<prim_max; ++k) {
          const double c = qmckl_ten3(coef_per_nucleus,k, j, i);
          if (c != 0.) {
            coef_mat_sparse_idx[i][j][l] = k;
            coef_mat_sparse[i][j][l] = c;
            l+=1;
          }
        }
        coef_mat_sparse_idx[i][j][0] = l;
      }
    }

#ifdef HAVE_OPENMP
#pragma omp for schedule(guided)
#endif
    for (int64_t ipoint=0 ; ipoint < point_num ; ++ipoint) {

/*
      double* const ptr = &ao_value[ipoint*ao_num];
      for (size_t k=0 ; k<ao_num ; ++k) {
        ptr[k] = 0.;
      }
*/

      memset(&ao_value[ipoint*ao_num], 0, ao_num*sizeof(double));

      const double e_coord[3]  __attribute__((aligned(64))) =
        { coord[ipoint],
          coord[ipoint + point_num],
          coord[ipoint + 2*point_num] };

      for (int64_t inucl=0 ; inucl < nucl_num ; ++inucl) {
        const double n_coord[3] __attribute__((aligned(64))) =
          { nucl_coord[inucl],
            nucl_coord[inucl + nucl_num],
            nucl_coord[inucl + 2*nucl_num] };

        /* Shift to avoid having exact zeros at the nodes of the AOs */
        const double shift = 1.e-20;

        const double x = e_coord[0] - n_coord[0];
        const double y = e_coord[1] - n_coord[1];
        const double z = e_coord[2] - n_coord[2];

        const double x2 = x*x;
        const double y2 = y*y;
        const double z2 = z*z;
        const double r2 = x2 + y2 + z2;

        /* Test if the point is in the range of the nucleus */
        if (r2 > cutoff * nucleus_range[inucl]) {
          continue;
        }

        int64_t n_poly;
        switch (nucleus_max_ang_mom[inucl]) {
        case 0:
          break;

        case 1:
          poly_vgl[0] = 1.;
          poly_vgl[1] = x;
          poly_vgl[2] = y;
          poly_vgl[3] = z;
          break;

        case 2:
          poly_vgl[0] = 1.;
          poly_vgl[1] = x+shift;
          poly_vgl[2] = y+shift;
          poly_vgl[3] = z+shift;
          poly_vgl[4] = x2;
          poly_vgl[5] = x*y;
          poly_vgl[6] = x*z;
          poly_vgl[7] = y2;
          poly_vgl[8] = y*z;
          poly_vgl[9] = z2;
          break;

        default:
          rc = qmckl_ao_polynomial_transp_vgl_hpc_inline(context, e_coord, n_coord,
                                                         nucleus_max_ang_mom[inucl],
                                                         &n_poly, powers, (int64_t) 3,
                                                         poly_vgl, size_max);
          assert (rc == QMCKL_SUCCESS);
          break;
        }

        /* Compute all exponents */

        int64_t nidx = 0;
        if (double_precision) {

          for (int64_t iprim = 0 ; iprim < prim_num_per_nucleus[inucl] ; ++iprim) {
            const double v = qmckl_mat(expo_per_nucleus, iprim, inucl) * r2;
            if (v <= cutoff) {
              ar2[iprim] = v;
              ++nidx;
            } else {
              break;
            }
          }

#ifdef HAVE_OPENMP
#pragma omp simd
#endif
          for (int64_t iprim = 0 ; iprim < nidx ; ++iprim) {
            exp_mat[iprim] = exp(-ar2[iprim]);
          }

        } else {

          for (int64_t iprim = 0 ; iprim < prim_num_per_nucleus[inucl] ; ++iprim) {
            const double v = qmckl_mat(expo_per_nucleus, iprim, inucl) * r2;
            if (v <= cutoff) {
              ar2_sp[iprim] = (float) v;
              ++nidx;
            } else {
              break;
            }
          }

#ifdef HAVE_OPENMP
#pragma omp simd
#endif
          for (int64_t iprim = 0 ; iprim < nidx ; ++iprim) {
            exp_mat_sp[iprim] = expf(-ar2_sp[iprim]);
          }

#ifdef HAVE_OPENMP
#pragma omp simd
#endif
          for (int64_t iprim = 0 ; iprim < nidx ; ++iprim) {
            exp_mat[iprim] = exp_mat_sp[iprim];
          }

        }

        for (int i=0 ; i<nucleus_shell_num[inucl] ; ++i) {
          double ce = 0.;
          const int32_t* restrict idx = &(coef_mat_sparse_idx[inucl][i][0]);
          const double*  restrict v = &(coef_mat_sparse[inucl][i][0]);
          for (int l=1 ; l<coef_mat_sparse_idx[inucl][i][0]; ++l) {
            if (l+1 < idx[0] && idx[l+1] < nidx)
              __builtin_prefetch(&(exp_mat[ idx[l+1] ]), 0, 1);
            const int k = idx[l];
            if (k >= nidx) break;
            ce += v[l] * exp_mat[k];
          }
          ce_mat[i] = ce;
        }

        const int64_t ishell_start = nucleus_index[inucl];
        const int64_t ishell_end   = nucleus_index[inucl] + nucleus_shell_num[inucl];

        const int64_t ipoint_ao_num = ipoint*ao_num;
        for (int64_t ishell = ishell_start ; ishell < ishell_end ; ++ishell) {

          const double s1 = ce_mat[ishell-ishell_start];

          if (s1 == 0.0) {
            continue;
          }

          const int64_t k = ao_index[ishell];
          double* restrict const ao_value_1 = ao_value + ipoint_ao_num + k;

          const int32_t l = shell_ang_mom[ishell];
          const int32_t n = lstart[l+1]-lstart[l];

          double* restrict poly_vgl_1 = NULL;
          if (nidx > 0) {
            const double* restrict f = ao_factor + k;
            const int64_t idx = lstart[l];

            poly_vgl_1 = &(poly_vgl[idx]);

            switch (n) {
            case 1:
              ao_value_1[0] = s1 * f[0];
              break;
            case 3:
#ifdef HAVE_OPENMP
#pragma omp simd
#endif
              for (int il=0 ; il<3 ; ++il) {
                ao_value_1[il] =  poly_vgl_1[il] * s1 * f[il];
              }
              break;
            case(6):
#ifdef HAVE_OPENMP
#pragma omp simd
#endif
              for (int il=0 ; il<6 ; ++il) {
                ao_value_1[il] =  poly_vgl_1[il] * s1 * f[il];
              }
              break;
            default:
#ifdef HAVE_OPENMP
#pragma omp simd
#endif
              for (int il=0 ; il<n ; ++il) {
                ao_value_1[il] =  poly_vgl_1[il] * s1 * f[il];
              }
              break;
            }
          }
        }
      }
    }
    free(coef_mat_sparse_idx);
    free(coef_mat_sparse);
  }

  return QMCKL_SUCCESS;
}
#endif



/* #+CALL: write_provider_pre( group="ao_basis", data="ao_value", dimension="ctx->ao_basis.ao_num * ctx->point.num") */

/* #+RESULTS: */

qmckl_exit_code qmckl_provide_ao_basis_ao_value(qmckl_context context)
{

  qmckl_exit_code rc = QMCKL_SUCCESS;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_provide_ao_basis_ao_value",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->ao_basis.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_provide_ao_basis_ao_value",
                           NULL);
  }

  /* Compute if necessary */
  if (ctx->point.date > ctx->ao_basis.ao_value_date) {

    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = ctx->ao_basis.ao_num * ctx->point.num * sizeof(double);

    if (mem_info.size > ctx->ao_basis.ao_value_maxsize) {
      if (ctx->ao_basis.ao_value != NULL) {
          rc = qmckl_free(context, ctx->ao_basis.ao_value);
          assert (rc == QMCKL_SUCCESS);
          ctx->ao_basis.ao_value = NULL;
        }
    }

    /* Allocate array */
    if (ctx->ao_basis.ao_value == NULL) {

      double* ao_value = (double*) qmckl_malloc(context, mem_info);
      ctx->ao_basis.ao_value_maxsize = mem_info.size;

      if (ao_value == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_ao_basis_ao_value",
                               NULL);
      }
      ctx->ao_basis.ao_value = ao_value;
    }

if (ctx->ao_basis.ao_vgl_date == ctx->point.date) {

      // ao_vgl has been computed at this step: Just copy the data.

      double * v = &(ctx->ao_basis.ao_value[0]);
      double * vgl = &(ctx->ao_basis.ao_vgl[0]);
      for (int i=0 ; i<ctx->point.num ; ++i) {
        for (int k=0 ; k<ctx->ao_basis.ao_num ; ++k) {
          v[k] = vgl[k];
        }
        v   += ctx->ao_basis.ao_num;
        vgl += ctx->ao_basis.ao_num * 5;
      }

    } else {

#ifdef HAVE_HPC
      if (ctx->ao_basis.type == 'G') {
        rc = qmckl_compute_ao_value_hpc_gaussian(context,
                                                 ctx->ao_basis.ao_num,
                                                 ctx->ao_basis.shell_num,
                                                 ctx->ao_basis.prim_num_per_nucleus,
                                                 ctx->point.num,
                                                 ctx->nucleus.num,
                                                 ctx->point.coord.data,
                                                 ctx->nucleus.coord.data,
                                                 ctx->ao_basis.nucleus_index,
                                                 ctx->ao_basis.nucleus_shell_num,
                                                 ctx->ao_basis.nucleus_range,
                                                 ctx->ao_basis.nucleus_max_ang_mom,
                                                 ctx->ao_basis.shell_ang_mom,
                                                 ctx->ao_basis.ao_factor,
                                                 ctx->ao_basis.expo_per_nucleus,
                                                 ctx->ao_basis.coef_per_nucleus,
                                                 ctx->ao_basis.ao_value);
        /*
          } else if (ctx->ao_basis.type == 'S') {
          rc = qmck_compute_ao_value_hpc_slater(context,
          ctx->ao_basis.ao_num,
          ctx->ao_basis.shell_num,
          ctx->ao_basis.prim_num,
          ctx->point.num,
          ctx->nucleus.num,
          ctx->point.coord.data,
          ctx->nucleus.coord.data,
          ctx->ao_basis.nucleus_index,
          ctx->ao_basis.nucleus_shell_num,
          ctx->ao_basis.nucleus_range,
          ctx->ao_basis.nucleus_max_ang_mom,
          ctx->ao_basis.shell_ang_mom,
          ctx->ao_basis.shell_prim_index,
          ctx->ao_basis.shell_prim_num,
          ctx->ao_basis.ao_factor,
          ctx->ao_basis.exponent,
          ctx->ao_basis.coefficient_normalized,
          ctx->ao_basis.ao_value);
        */
      } else {
        /* Provide required data */
        rc = qmckl_provide_ao_basis_shell_vgl(context);
        if (rc != QMCKL_SUCCESS) {
            return qmckl_failwith( context, rc, "qmckl_provide_ao_basis_shell_vgl", NULL);
        }

        rc = qmckl_compute_ao_value_doc(context,
                                        ctx->ao_basis.ao_num,
                                        ctx->ao_basis.shell_num,
                                        ctx->point.num,
                                        ctx->nucleus.num,
                                        ctx->point.coord.data,
                                        ctx->nucleus.coord.data,
                                        ctx->ao_basis.nucleus_index,
                                        ctx->ao_basis.nucleus_shell_num,
                                        ctx->ao_basis.nucleus_range,
                                        ctx->ao_basis.nucleus_max_ang_mom,
                                        ctx->ao_basis.shell_ang_mom,
                                        ctx->ao_basis.ao_factor,
                                        ctx->ao_basis.shell_vgl,
                                        ctx->ao_basis.ao_value);
      }
#else
      /* Provide required data */
      rc = qmckl_provide_ao_basis_shell_vgl(context);
      if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc, "qmckl_provide_ao_basis_shell_vgl", NULL);
      }

      rc = qmckl_compute_ao_value_doc(context,
                                      ctx->ao_basis.ao_num,
                                      ctx->ao_basis.shell_num,
                                      ctx->point.num,
                                      ctx->nucleus.num,
                                      ctx->point.coord.data,
                                      ctx->nucleus.coord.data,
                                      ctx->ao_basis.nucleus_index,
                                      ctx->ao_basis.nucleus_shell_num,
                                      ctx->ao_basis.nucleus_range,
                                      ctx->ao_basis.nucleus_max_ang_mom,
                                      ctx->ao_basis.shell_ang_mom,
                                      ctx->ao_basis.ao_factor,
                                      ctx->ao_basis.shell_vgl,
                                      ctx->ao_basis.ao_value);
#endif
    }



/* #+CALL: write_provider_post( group="ao_basis", data="ao_value" ) */

/* #+RESULTS: */

if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->ao_basis.ao_value_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

/* HPC version                                                    :noexport: */
/*      #+NAME: qmckl_ao_vgl_args_hpc_gaussian */
/*     | Variable              | Type                           | In/Out | Description                               | */
/*     |-----------------------+--------------------------------+--------+-------------------------------------------| */
/*     | ~context~             | ~qmckl_context~                | in     | Global state                              | */
/*     | ~ao_num~              | ~int64_t~                      | in     | Number of AOs                             | */
/*     | ~shell_num~           | ~int64_t~                      | in     | Number of shells                          | */
/*     | ~prim_num~            | ~int64_t~                      | in     | Number of primitives                      | */
/*     | ~point_num~           | ~int64_t~                      | in     | Number of points                          | */
/*     | ~nucl_num~            | ~int64_t~                      | in     | Number of nuclei                          | */
/*     | ~coord~               | ~double[3][point_num]~         | in     | Coordinates                               | */
/*     | ~nucl_coord~          | ~double[3][nucl_num]~          | in     | Nuclear  coordinates                      | */
/*     | ~nucleus_index~       | ~int64_t[nucl_num]~            | in     | Index of the 1st shell of each nucleus    | */
/*     | ~nucleus_shell_num~   | ~int64_t[nucl_num]~            | in     | Number of shells per nucleus              | */
/*     | ~nucleus_range~       | ~double[nucl_num]~             | in     | Range beyond which all is zero            | */
/*     | ~nucleus_max_ang_mom~ | ~int32_t[nucl_num]~            | in     | Maximum angular momentum per nucleus      | */
/*     | ~shell_ang_mom~       | ~int32_t[shell_num]~           | in     | Angular momentum of each shell            | */
/*     | ~shell_prim_index~    | ~int64_t[shell_num]~           | in     | Index of the 1st primitive of each shell  | */
/*     | ~shell_prim_num~      | ~int64_t[shell_num]~           | in     | Number of primitives per shell            | */
/*     | ~ao_factor~           | ~double[ao_num]~               | in     | Normalization factor of the AOs           | */
/*     | ~ao_expo~             | ~double[prim_num]~             | in     | Exponents of the primitives               | */
/*     | ~coef_normalized~     | ~double[prim_num]~             | in     | Normalized coefficients of the primitives | */
/*     | ~ao_vgl~              | ~double[point_num][5][ao_num]~ | out    | Value, gradients and Laplacian of the AOs | */



#ifdef HAVE_HPC
qmckl_exit_code
qmckl_compute_ao_vgl_hpc_gaussian (
                                   const qmckl_context context,
                                   const int64_t ao_num,
                                   const int64_t shell_num,
                                   const int32_t* restrict prim_num_per_nucleus,
                                   const int64_t point_num,
                                   const int64_t nucl_num,
                                   const double* restrict coord,
                                   const double* restrict nucl_coord,
                                   const int64_t* restrict nucleus_index,
                                   const int64_t* restrict nucleus_shell_num,
                                   const double* nucleus_range,
                                   const int32_t* restrict nucleus_max_ang_mom,
                                   const int32_t* restrict shell_ang_mom,
                                   const double* restrict ao_factor,
                                   const qmckl_matrix expo_per_nucleus,
                                   const qmckl_tensor coef_per_nucleus,
                                   double* restrict const ao_vgl )
{
  int32_t lstart[32] __attribute__((aligned(64)));
  for (int32_t l=0 ; l<32 ; ++l) {
    lstart[l] = l*(l+1)*(l+2)/6;
  }
  
  int64_t ao_index[shell_num+1] __attribute__((aligned(64)));
  int64_t size_max = 0;
  int64_t prim_max = 0;
  int64_t shell_max = 0;
  {
    int64_t k=0;
    for (int inucl=0 ; inucl < nucl_num ; ++inucl) {
      prim_max = prim_num_per_nucleus[inucl] > prim_max ?
        prim_num_per_nucleus[inucl] : prim_max;
      shell_max = nucleus_shell_num[inucl] > shell_max ?
        nucleus_shell_num[inucl] : shell_max;
      const int64_t ishell_start = nucleus_index[inucl];
      const int64_t ishell_end   = nucleus_index[inucl] + nucleus_shell_num[inucl];
      for (int64_t ishell = ishell_start ; ishell < ishell_end ; ++ishell) {
        const int l = shell_ang_mom[ishell];
        ao_index[ishell] = k;
        k += lstart[l+1] - lstart[l];
        size_max = size_max < lstart[l+1] ? lstart[l+1] : size_max;
      }
    }
    ao_index[shell_num] = ao_num+1;
  }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);
  
  /* Don't compute polynomials when the radial part is zero. */
  const int precision = ctx->numprec.precision;
  const double cutoff = log( (double) ( ((uint64_t) 1) << (precision-1) ) );
  const bool double_precision = precision > 23;
  

#ifdef HAVE_OPENMP
#pragma omp parallel if (point_num > 16)
#endif
  {
    qmckl_exit_code rc;
    double ar2[prim_max]  __attribute__((aligned(64)));
    float  ar2_sp[prim_max] __attribute__((aligned(64)));
    int32_t powers[3*size_max]  __attribute__((aligned(64)));

    // 0: value 1-3: grad 4: Laplacian. We use 8 for vectorization/padding
    double exp_mat[prim_max][8] __attribute__((aligned(64))) ;
    float  exp_mat_sp[prim_max] __attribute__((aligned(64)));
    double ce_mat[shell_max][8] __attribute__((aligned(64))) ;

    int32_t (*coef_mat_sparse_idx)[shell_max][prim_max+1] __attribute__((aligned(64))) =
      malloc(nucl_num * sizeof (*coef_mat_sparse_idx));

    double  (*coef_mat_sparse)[shell_max][prim_max+1] __attribute__((aligned(64))) =
      malloc(nucl_num * sizeof (*coef_mat_sparse));

    for (int i=0 ; i<nucl_num ; ++i) {
      for (int j=0 ; j<shell_max; ++j) {
        int l=1;
        for (int k=0 ; k<prim_max; ++k) {
          const double c = qmckl_ten3(coef_per_nucleus,k, j, i);
          if (c != 0.) {
            coef_mat_sparse_idx[i][j][l] = k;
            coef_mat_sparse[i][j][l] = c;
            l+=1;
          }
        }
        coef_mat_sparse_idx[i][j][0] = l;
      }
    }

    double  poly_vgl_l1[4][4] __attribute__((aligned(64))) =
      {{1.0, 0.0, 0.0, 0.0},
       {0.0, 1.0, 0.0, 0.0},
       {0.0, 0.0, 1.0, 0.0},
       {0.0, 0.0, 0.0, 1.0}};
    double  poly_vgl_l2[5][10]__attribute__((aligned(64))) =
      {{1., 0., 0., 0., 0., 0., 0., 0., 0., 0.},
       {0., 1., 0., 0., 0., 0., 0., 0., 0., 0.},
       {0., 0., 1., 0., 0., 0., 0., 0., 0., 0.},
       {0., 0., 0., 1., 0., 0., 0., 0., 0., 0.},
       {0., 0., 0., 0., 2., 0., 0., 2., 0., 2.}};
    double  poly_vgl[5][size_max]  __attribute__((aligned(64)));

#ifdef HAVE_OPENMP
#pragma omp for
#endif
    for (int64_t ipoint=0 ; ipoint < point_num ; ++ipoint) {

/*
      double* const ptr = &ao_vgl[ipoint*5*ao_num];
      for (size_t k=0 ; k<5*ao_num ; ++k) {
        ptr[k] = 0.;
      }
*/

      memset(&ao_vgl[ipoint*ao_num*5], 0, 5*ao_num*sizeof(double));

      const double e_coord[3]  __attribute__((aligned(64))) =
        { coord[ipoint],
          coord[ipoint + point_num],
          coord[ipoint + 2*point_num] };

      for (int64_t inucl=0 ; inucl < nucl_num ; ++inucl) {
        const double n_coord[3] __attribute__((aligned(64))) =
          { nucl_coord[inucl],
            nucl_coord[inucl + nucl_num],
            nucl_coord[inucl + 2*nucl_num] };

        const double shift = 1.e-20;

        const double x = e_coord[0] - n_coord[0];
        const double y = e_coord[1] - n_coord[1];
        const double z = e_coord[2] - n_coord[2];

        const double x2 = x*x;
        const double y2 = y*y;
        const double z2 = z*z;
        const double r2 = x2 + y2 + z2;

        /* Test if the point is in the range of the nucleus */
        if (r2 > cutoff * nucleus_range[inucl]) {
          continue;
        }

        int64_t n_poly;
        switch (nucleus_max_ang_mom[inucl]) {
        case 0:
          break;

        case 1:
          poly_vgl_l1[0][1] = x+shift;
          poly_vgl_l1[0][2] = y+shift;
          poly_vgl_l1[0][3] = z+shift;
          break;

        case 2:
          poly_vgl_l2[0][1] = x+shift;
          poly_vgl_l2[0][2] = y+shift;
          poly_vgl_l2[0][3] = z+shift;
          poly_vgl_l2[0][4] = x2+shift;
          poly_vgl_l2[0][5] = x*y+shift;
          poly_vgl_l2[0][6] = x*z+shift;
          poly_vgl_l2[0][7] = y2+shift;
          poly_vgl_l2[0][8] = y*z+shift;
          poly_vgl_l2[0][9] = z2+shift;
          poly_vgl_l2[1][4] = x+x;
          poly_vgl_l2[1][5] = y;
          poly_vgl_l2[1][6] = z;
          poly_vgl_l2[2][5] = x;
          poly_vgl_l2[2][7] = y+y;
          poly_vgl_l2[2][8] = z;
          poly_vgl_l2[3][6] = x;
          poly_vgl_l2[3][8] = y;
          poly_vgl_l2[3][9] = z+z;
          break;

        default:
          rc = qmckl_ao_polynomial_transp_vgl_hpc_inline(context, e_coord, n_coord,
                                                         nucleus_max_ang_mom[inucl],
                                                         &n_poly, powers, (int64_t) 3,
                                                         &(poly_vgl[0][0]), size_max);
          assert (rc == QMCKL_SUCCESS);
          break;
        }

        /* Compute all exponents */

        int64_t nidx = 0;

        if (double_precision) {

          for (int64_t iprim = 0 ; iprim < prim_num_per_nucleus[inucl] ; ++iprim) {
            const double v = qmckl_mat(expo_per_nucleus, iprim, inucl) * r2;
            if (v <= cutoff) {
              ar2[iprim] = v;
              ++nidx;
            } else {
              break;
            }
          }

#ifdef HAVE_OPENMP
#pragma omp simd
#endif
          for (int64_t iprim = 0 ; iprim < nidx ; ++iprim) {
            exp_mat[iprim][0] = exp(-ar2[iprim]);
          }

        } else {

          for (int64_t iprim = 0 ; iprim < prim_num_per_nucleus[inucl] ; ++iprim) {
            const double v = qmckl_mat(expo_per_nucleus, iprim, inucl) * r2;
            if (v <= cutoff) {
              ar2[iprim] = v;
              ar2_sp[iprim] = (float) v;
              ++nidx;
            } else {
              break;
            }
          }

#ifdef HAVE_OPENMP
#pragma omp simd
#endif
          for (int64_t iprim = 0 ; iprim < nidx ; ++iprim) {
            exp_mat_sp[iprim] = expf(-ar2_sp[iprim]);
          }

#ifdef HAVE_OPENMP
#pragma omp simd
#endif
          for (int64_t iprim = 0 ; iprim < nidx ; ++iprim) {
            exp_mat[iprim][0] = exp_mat_sp[iprim];
          }

        }

#ifdef HAVE_OPENMP
#pragma omp simd
#endif
        for (int64_t iprim = 0 ; iprim < nidx ; ++iprim) {
          double f = qmckl_mat(expo_per_nucleus, iprim, inucl) * exp_mat[iprim][0];
          f = -f-f;
          exp_mat[iprim][1] = f * x;
          exp_mat[iprim][2] = f * y;
          exp_mat[iprim][3] = f * z;
          exp_mat[iprim][4] = f * (3.0 - 2.0 * ar2[iprim]);
        }


        /* --- */

        //        if (do_sparse) {
        for (int i=0 ; i<nucleus_shell_num[inucl] ; ++i) {
          double ce[8] = {0.};
          const int32_t* restrict idx = &(coef_mat_sparse_idx[inucl][i][0]);
          const double* restrict v = &(coef_mat_sparse[inucl][i][0]);
          for (int l=1 ; l<idx[0]; ++l) {
            if (l+1 < idx[0] && idx[l+1] < nidx)
              __builtin_prefetch(&(exp_mat[ idx[l+1] ][0]), 0, 1);
            if (idx[l] >= nidx) break;
            double alpha = v[l];
            const double* restrict e = exp_mat[ idx[l] ];
            for (int j=0 ; j<8 ; ++j) {
              ce[j] += alpha * e[j];
            }
          }
          for (int j=0 ; j<8 ; ++j) {
            ce_mat[i][j] = ce[j];
          }
        }

        const int64_t ishell_start = nucleus_index[inucl];
        const int64_t ishell_end   = nucleus_index[inucl] + nucleus_shell_num[inucl];

        const int64_t ipoint_5_ao_num = ipoint*5*ao_num;

        for (int64_t ishell = ishell_start ; ishell < ishell_end ; ++ishell) {

          const double s1 = ce_mat[ishell-ishell_start][0];

          const int64_t k = ao_index[ishell];
          double* restrict const ao_vgl_1 = ao_vgl + ipoint_5_ao_num + k;

          const int32_t l = shell_ang_mom[ishell];
          const int32_t n = lstart[l+1]-lstart[l];

          double* restrict const ao_vgl_2 = ao_vgl_1 + ao_num;
          double* restrict const ao_vgl_3 = ao_vgl_1 + (ao_num<<1);
          double* restrict const ao_vgl_4 = ao_vgl_1 + (ao_num<<1) + ao_num;
          double* restrict const ao_vgl_5 = ao_vgl_1 + (ao_num<<2);

          if (s1 == 0.0) {
/*
            for (int64_t il=0 ; il<n ; ++il) {
              ao_vgl_1[il] = 0.0;
              ao_vgl_2[il] = 0.0;
              ao_vgl_3[il] = 0.0;
              ao_vgl_4[il] = 0.0;
              ao_vgl_5[il] = 0.0;
            }
*/
            continue;
          }

          const double s2 = ce_mat[ishell-ishell_start][1];
          const double s3 = ce_mat[ishell-ishell_start][2];
          const double s4 = ce_mat[ishell-ishell_start][3];
          const double s5 = ce_mat[ishell-ishell_start][4];

          double* restrict poly_vgl_1 = NULL;
          double* restrict poly_vgl_2 = NULL;
          double* restrict poly_vgl_3 = NULL;
          double* restrict poly_vgl_4 = NULL;
          double* restrict poly_vgl_5 = NULL;
          if (nidx > 0) {
            const double* restrict f = ao_factor + k;
            const int64_t idx = lstart[l];

            switch (nucleus_max_ang_mom[inucl]) {
            case 0:
              break;
            case 1:
              poly_vgl_1 = &(poly_vgl_l1[0][idx]);
              poly_vgl_2 = &(poly_vgl_l1[1][idx]);
              poly_vgl_3 = &(poly_vgl_l1[2][idx]);
              poly_vgl_4 = &(poly_vgl_l1[3][idx]);
              break;
            case 2:
              poly_vgl_1 = &(poly_vgl_l2[0][idx]);
              poly_vgl_2 = &(poly_vgl_l2[1][idx]);
              poly_vgl_3 = &(poly_vgl_l2[2][idx]);
              poly_vgl_4 = &(poly_vgl_l2[3][idx]);
              poly_vgl_5 = &(poly_vgl_l2[4][idx]);
              break;
            default:
              poly_vgl_1 = &(poly_vgl[0][idx]);
              poly_vgl_2 = &(poly_vgl[1][idx]);
              poly_vgl_3 = &(poly_vgl[2][idx]);
              poly_vgl_4 = &(poly_vgl[3][idx]);
              poly_vgl_5 = &(poly_vgl[4][idx]);
            }
            switch (n) {
            case 1:
              ao_vgl_1[0] = s1 * f[0];
              ao_vgl_2[0] = s2 * f[0];
              ao_vgl_3[0] = s3 * f[0];
              ao_vgl_4[0] = s4 * f[0];
              ao_vgl_5[0] = s5 * f[0];
              break;
            case 3:
#ifdef HAVE_OPENMP
#pragma omp simd
#endif
              for (int il=0 ; il<3 ; ++il) {
                ao_vgl_1[il] =  poly_vgl_1[il] * s1 * f[il];
                ao_vgl_2[il] = (poly_vgl_2[il] * s1 + poly_vgl_1[il] * s2) * f[il];
                ao_vgl_3[il] = (poly_vgl_3[il] * s1 + poly_vgl_1[il] * s3) * f[il];
                ao_vgl_4[il] = (poly_vgl_4[il] * s1 + poly_vgl_1[il] * s4) * f[il];
                ao_vgl_5[il] = (poly_vgl_1[il] * s5 +
                                2.0*(poly_vgl_2[il] * s2 +
                                     poly_vgl_3[il] * s3 +
                                     poly_vgl_4[il] * s4 )) * f[il];
              }
              break;
            case 6:
#ifdef HAVE_OPENMP
#pragma omp simd
#endif
              for (int il=0 ; il<6 ; ++il) {
                ao_vgl_1[il] =  poly_vgl_1[il] * s1 * f[il];
                ao_vgl_2[il] = (poly_vgl_2[il] * s1 + poly_vgl_1[il] * s2) * f[il];
                ao_vgl_3[il] = (poly_vgl_3[il] * s1 + poly_vgl_1[il] * s3) * f[il];
                ao_vgl_4[il] = (poly_vgl_4[il] * s1 + poly_vgl_1[il] * s4) * f[il];
                ao_vgl_5[il] = (poly_vgl_5[il] * s1 + poly_vgl_1[il] * s5 +
                                2.0*(poly_vgl_2[il] * s2 +
                                     poly_vgl_3[il] * s3 +
                                     poly_vgl_4[il] * s4 )) * f[il];
              }
              break;
            default:
#ifdef HAVE_OPENMP
#pragma omp simd
#endif
              for (int il=0 ; il<n ; ++il) {
                ao_vgl_1[il] =  poly_vgl_1[il] * s1 * f[il];
                ao_vgl_2[il] = (poly_vgl_2[il] * s1 + poly_vgl_1[il] * s2) * f[il];
                ao_vgl_3[il] = (poly_vgl_3[il] * s1 + poly_vgl_1[il] * s3) * f[il];
                ao_vgl_4[il] = (poly_vgl_4[il] * s1 + poly_vgl_1[il] * s4) * f[il];
                ao_vgl_5[il] = (poly_vgl_5[il] * s1 + poly_vgl_1[il] * s5 +
                                2.0*(poly_vgl_2[il] * s2 +
                                     poly_vgl_3[il] * s3 +
                                     poly_vgl_4[il] * s4 )) * f[il];
              }
              break;
            }
/*
          } else {
            for (int64_t il=0 ; il<n ; ++il) {
              ao_vgl_1[il] = 0.0;
              ao_vgl_2[il] = 0.0;
              ao_vgl_3[il] = 0.0;
              ao_vgl_4[il] = 0.0;
              ao_vgl_5[il] = 0.0;
            }
*/
          }
        }
      }
    }
    free(coef_mat_sparse_idx);
    free(coef_mat_sparse);
  }

  return QMCKL_SUCCESS;
}
#endif



/* #+CALL: write_provider_pre( group="ao_basis", data="ao_vgl", dimension="ctx->ao_basis.ao_num * 5 * ctx->point.num") */

/* #+RESULTS: */

qmckl_exit_code qmckl_provide_ao_basis_ao_vgl(qmckl_context context)
{

  qmckl_exit_code rc = QMCKL_SUCCESS;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_provide_ao_basis_ao_vgl",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->ao_basis.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_provide_ao_basis_ao_vgl",
                           NULL);
  }

  /* Compute if necessary */
  if (ctx->point.date > ctx->ao_basis.ao_vgl_date) {

    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = ctx->ao_basis.ao_num * 5 * ctx->point.num * sizeof(double);


    if (mem_info.size > ctx->ao_basis.ao_vgl_maxsize) {
      if (ctx->ao_basis.ao_vgl != NULL) {
          rc = qmckl_free(context, ctx->ao_basis.ao_vgl);
          assert (rc == QMCKL_SUCCESS);
          ctx->ao_basis.ao_vgl = NULL;
        }
    }

    /* Allocate array */
    if (ctx->ao_basis.ao_vgl == NULL) {

      double* ao_vgl = (double*) qmckl_malloc(context, mem_info);
      ctx->ao_basis.ao_vgl_maxsize = mem_info.size;
      if (ao_vgl == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_ao_basis_ao_vgl",
                               NULL);
      }
      ctx->ao_basis.ao_vgl = ao_vgl;
    }

#ifdef HAVE_HPC
if (ctx->ao_basis.type == 'G') {
  rc = qmckl_compute_ao_vgl_hpc_gaussian(context,
                                         ctx->ao_basis.ao_num,
                                         ctx->ao_basis.shell_num,
                                         ctx->ao_basis.prim_num_per_nucleus,
                                         ctx->point.num,
                                         ctx->nucleus.num,
                                         ctx->point.coord.data,
                                         ctx->nucleus.coord.data,
                                         ctx->ao_basis.nucleus_index,
                                         ctx->ao_basis.nucleus_shell_num,
                                         ctx->ao_basis.nucleus_range,
                                         ctx->ao_basis.nucleus_max_ang_mom,
                                         ctx->ao_basis.shell_ang_mom,
                                         ctx->ao_basis.ao_factor,
                                         ctx->ao_basis.expo_per_nucleus,
                                         ctx->ao_basis.coef_per_nucleus,
                                         ctx->ao_basis.ao_vgl);

  /* DEBUG
     rc = qmckl_provide_ao_basis_shell_vgl(context);
     if (rc != QMCKL_SUCCESS) {
     return qmckl_failwith( context, rc, "qmckl_provide_ao_basis_shell_vgl", NULL);
     }
     int64_t K= ctx->ao_basis.ao_num * 5 * ctx->point.num;
     double* check = malloc(K*sizeof(double));
     rc = qmckl_compute_ao_vgl_doc(context,
     ctx->ao_basis.ao_num,
     ctx->ao_basis.shell_num,
     ctx->point.num,
     ctx->nucleus.num,
     ctx->point.coord.data,
     ctx->nucleus.coord.data,
     ctx->ao_basis.nucleus_index,
     ctx->ao_basis.nucleus_shell_num,
     ctx->ao_basis.nucleus_range,
     ctx->ao_basis.nucleus_max_ang_mom,
     ctx->ao_basis.shell_ang_mom,
     ctx->ao_basis.ao_factor,
     ctx->ao_basis.shell_vgl,
     check);
     for (int64_t i=0 ; i<K ; ++i) {
     if (fabs(check[i] - ctx->ao_basis.ao_vgl[i]) > 1.e-10) {
     int a, b, c;
     a = i/(ctx->ao_basis.ao_num*5);
     b = (i-a*ctx->ao_basis.ao_num*5)/ctx->ao_basis.ao_num;
     c = (i-a*ctx->ao_basis.ao_num*5 -b*ctx->ao_basis.ao_num);
     printf("%d: %d, %d, %d, %e  %e\n", i, a, b, c, check[i], ctx->ao_basis.ao_vgl[i]);
     }
     }
  */

  /*
    } else if (ctx->ao_basis.type == 'S') {
    rc = qmck_compute_ao_vgl_hpc_slater(context,
    ctx->ao_basis.ao_num,
    ctx->ao_basis.shell_num,
    ctx->ao_basis.prim_num,
    ctx->point.num,
    ctx->nucleus.num,
    ctx->point.coord.data,
    ctx->nucleus.coord.data,
    ctx->ao_basis.nucleus_index,
    ctx->ao_basis.nucleus_shell_num,
    ctx->ao_basis.nucleus_range,
    ctx->ao_basis.nucleus_max_ang_mom,
    ctx->ao_basis.shell_ang_mom,
    ctx->ao_basis.shell_prim_index,
    ctx->ao_basis.shell_prim_num,
    ctx->ao_basis.ao_factor,
    ctx->ao_basis.exponent,
    ctx->ao_basis.coefficient_normalized,
    ctx->ao_basis.ao_vgl);
  */
 } else {
  /* Provide required data */

  rc = qmckl_provide_ao_basis_shell_vgl(context);
  if (rc != QMCKL_SUCCESS) {
    return qmckl_failwith( context, rc, "qmckl_provide_ao_basis_shell_vgl", NULL);
  }

  rc = qmckl_compute_ao_vgl_doc(context,
                                ctx->ao_basis.ao_num,
                                ctx->ao_basis.shell_num,
                                ctx->point.num,
                                ctx->nucleus.num,
                                ctx->point.coord.data,
                                ctx->nucleus.coord.data,
                                ctx->ao_basis.nucleus_index,
                                ctx->ao_basis.nucleus_shell_num,
                                ctx->ao_basis.nucleus_range,
                                ctx->ao_basis.nucleus_max_ang_mom,
                                ctx->ao_basis.shell_ang_mom,
                                ctx->ao_basis.ao_factor,
                                ctx->ao_basis.shell_vgl,
                                ctx->ao_basis.ao_vgl);
 }
#else
rc = qmckl_provide_ao_basis_shell_vgl(context);
if (rc != QMCKL_SUCCESS) {
  return qmckl_failwith( context, rc, "qmckl_provide_ao_basis_shell_vgl", NULL);
 }

rc = qmckl_compute_ao_vgl_doc(context,
                              ctx->ao_basis.ao_num,
                              ctx->ao_basis.shell_num,
                              ctx->point.num,
                              ctx->nucleus.num,
                              ctx->point.coord.data,
                              ctx->nucleus.coord.data,
                              ctx->ao_basis.nucleus_index,
                              ctx->ao_basis.nucleus_shell_num,
                              ctx->ao_basis.nucleus_range,
                              ctx->ao_basis.nucleus_max_ang_mom,
                              ctx->ao_basis.shell_ang_mom,
                              ctx->ao_basis.ao_factor,
                              ctx->ao_basis.shell_vgl,
                              ctx->ao_basis.ao_vgl);
#endif



/* #+CALL: write_provider_post( group="ao_basis", data="ao_vgl" ) */

/* #+RESULTS: */

if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->ao_basis.ao_vgl_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_ao_basis_ao_hessian(qmckl_context context,
                          double* const ao_hessian,
                          const int64_t size_max)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_ao_basis_ao_hessian(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  const int64_t sze = ctx->ao_basis.ao_num * 4 * 3 * ctx->point.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_ao_basis_ao_hessian",
                           "input array too small");
  }

  memcpy(ao_hessian, ctx->ao_basis.ao_hessian, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_ao_basis_ao_hessian(qmckl_context context)
{

  qmckl_exit_code rc = QMCKL_SUCCESS;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_provide_ao_basis_ao_hessian",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->ao_basis.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_provide_ao_basis_ao_hessian",
                           NULL);
  }

  /* Compute if necessary */
  if (ctx->point.date > ctx->ao_basis.ao_hessian_date) {

    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = ctx->ao_basis.ao_num * 4 * 3 * ctx->point.num * sizeof(double);

    if (ctx->ao_basis.ao_hessian != NULL) {
      qmckl_memory_info_struct mem_info_test = qmckl_memory_info_struct_zero;
      rc = qmckl_get_malloc_info(context, ctx->ao_basis.ao_hessian, &mem_info_test);

      /* if rc != QMCKL_SUCCESS, we are maybe in an _inplace function because the
         memory was not allocated with qmckl_malloc */

      if ((rc == QMCKL_SUCCESS) && (mem_info_test.size != mem_info.size)) {
        rc = qmckl_free(context, ctx->ao_basis.ao_hessian);
        assert (rc == QMCKL_SUCCESS);
        ctx->ao_basis.ao_hessian = NULL;
      }
    }

    /* Allocate array */
    if (ctx->ao_basis.ao_hessian == NULL) {

      double* ao_hessian = (double*) qmckl_malloc(context, mem_info);

      if (ao_hessian == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_ao_basis_ao_hessian",
                               NULL);
      }
      ctx->ao_basis.ao_hessian = ao_hessian;
    }


    rc = qmckl_provide_ao_basis_shell_vgl(context);
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context, rc, "qmckl_provide_ao_basis_shell_vgl", NULL);
    }

    rc = qmckl_provide_ao_basis_shell_hessian(context);
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context, rc, "qmckl_provide_ao_basis_shell_hessian", NULL);
    }

    //memset(ctx->ao_basis.ao_hessian, 0, mem_info.size);

    rc = qmckl_compute_ao_hessian_doc(context,
                                  ctx->ao_basis.ao_num,
                                  ctx->ao_basis.shell_num,
                                  ctx->point.num,
                                  ctx->nucleus.num,
                                  ctx->point.coord.data,
                                  ctx->nucleus.coord.data,
                                  ctx->ao_basis.nucleus_index,
                                  ctx->ao_basis.nucleus_shell_num,
                                  ctx->ao_basis.nucleus_range,
                                  ctx->ao_basis.nucleus_max_ang_mom,
                                  ctx->ao_basis.shell_ang_mom,
                                  ctx->ao_basis.ao_factor,
                                  ctx->ao_basis.shell_vgl,
                                  ctx->ao_basis.shell_hessian,
                                  ctx->ao_basis.ao_hessian);

    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->ao_basis.ao_hessian_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}
