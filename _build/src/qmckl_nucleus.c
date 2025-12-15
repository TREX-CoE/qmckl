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
#include "qmckl_nucleus_private_func.h"

qmckl_exit_code qmckl_init_nucleus(qmckl_context context) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return false;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  ctx->nucleus.uninitialized = (1 << 3) - 1;

  /* Default values */

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_nucleus_num (const qmckl_context context, int64_t* const num) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  }

  if (num == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_nucleus_num",
                           "num is a null pointer");
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 0;

  if ( (ctx->nucleus.uninitialized & mask) != 0) {
    *num = (int64_t) 0;
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_get_nucleus_num",
                           "nucleus data is not provided");
  }

  assert (ctx->nucleus.num >= (int64_t) 0);
  *num = ctx->nucleus.num;

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_nucleus_charge (const qmckl_context context,
                          double* const charge,
                          const int64_t size_max)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  }

  if (charge == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_nucleus_charge",
                           "charge is a null pointer");
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 1;

  if ( (ctx->nucleus.uninitialized & mask) != 0) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_get_nucleus_charge",
                           "nucleus data is not provided");
  }

  assert (ctx->nucleus.charge.data != NULL);

  if (ctx->nucleus.num > size_max) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_nucleus_charge",
                           "Array too small");
  }

  qmckl_exit_code rc;
  rc = qmckl_double_of_vector(context, ctx->nucleus.charge, charge, size_max);

  return rc;
}

qmckl_exit_code
qmckl_get_nucleus_coord (const qmckl_context context,
                         const char transp,
                         double* const coord,
                         const int64_t size_max)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  }

  if (transp != 'N' && transp != 'T') {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_nucleus_coord",
                           "transp should be 'N' or 'T'");
  }

  if (coord == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_nucleus_coord",
                           "coord is a null pointer");
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 2;

  if ( (ctx->nucleus.uninitialized & mask) != 0) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_get_nucleus_coord",
                           "nucleus data is not provided");
  }

  assert (ctx->nucleus.coord.data != NULL);

  qmckl_exit_code rc;

  if (transp == 'N') {
    qmckl_matrix At = qmckl_matrix_alloc(context, 3, ctx->nucleus.coord.size[0]);
    rc = qmckl_transpose(context, ctx->nucleus.coord, At);
    if (rc != QMCKL_SUCCESS) return rc;
    rc = qmckl_double_of_matrix(context, At, coord, size_max);
    qmckl_matrix_free(context, &At);
  } else {
    rc = qmckl_double_of_matrix(context, ctx->nucleus.coord, coord, size_max);
  }

  return rc;
}

bool qmckl_nucleus_provided(const qmckl_context context) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return false;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  return ctx->nucleus.provided;
}



/* Sets the number of nuclei. */


qmckl_exit_code
qmckl_set_nucleus_num(qmckl_context context,
                      const int64_t num)
{
  int32_t mask = 1 << 0;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
      return qmckl_failwith( context,
                             QMCKL_NULL_CONTEXT,
                             "qmckl_set_nucleus_*",
                             NULL);
  }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  
  if (mask != 0 && !(ctx->nucleus.uninitialized & mask)) {
      return qmckl_failwith( context,
                             QMCKL_ALREADY_SET,
                             "qmckl_set_nucleus_*",
                             NULL);
  }

  if (num <= 0) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_set_nucleus_num",
                           "num <= 0");
  }

  ctx->nucleus.num = num;

  ctx->nucleus.uninitialized &= ~mask;
  ctx->nucleus.provided = (ctx->nucleus.uninitialized == 0);
  
  return QMCKL_SUCCESS;
}



/* Sets the nuclear charges of all the atoms. */


qmckl_exit_code
qmckl_set_nucleus_charge(qmckl_context context,
                         const double* charge,
                         const int64_t size_max)
{
  int32_t mask = 1 << 1;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
      return qmckl_failwith( context,
                             QMCKL_NULL_CONTEXT,
                             "qmckl_set_nucleus_*",
                             NULL);
  }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  
  if (mask != 0 && !(ctx->nucleus.uninitialized & mask)) {
      return qmckl_failwith( context,
                             QMCKL_ALREADY_SET,
                             "qmckl_set_nucleus_*",
                             NULL);
  }

  if (charge == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_set_nucleus_charge",
                           "charge is a null pointer");
  }

  int64_t num;
  qmckl_exit_code rc;

  rc = qmckl_get_nucleus_num(context, &num);
  if (rc != QMCKL_SUCCESS) return rc;

  if (num > size_max) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_set_nucleus_charge",
                           "Array too small");
  }

  ctx->nucleus.charge = qmckl_vector_alloc(context, num);
  rc = qmckl_vector_of_double(context, charge, num, &(ctx->nucleus.charge));
  if (rc != QMCKL_SUCCESS) {
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "qmckl_set_nucleus_charge",
                           "Error in vector->double* conversion");
  }

  ctx->nucleus.uninitialized &= ~mask;
  ctx->nucleus.provided = (ctx->nucleus.uninitialized == 0);
  
  return QMCKL_SUCCESS;
}



/* Sets the nuclear coordinates of all the atoms. The coordinates */
/* are be given in atomic units. */


qmckl_exit_code
qmckl_set_nucleus_coord(qmckl_context context,
                        const char transp,
                        const double* coord,
                        const int64_t size_max)
{
  int32_t mask = 1 << 2;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
      return qmckl_failwith( context,
                            QMCKL_NULL_CONTEXT,
                            "qmckl_set_nucleus_*",
                            NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;

  qmckl_exit_code rc;

  const int64_t nucl_num = (int64_t) ctx->nucleus.num;

  if (ctx->nucleus.coord.data != NULL) {
    rc = qmckl_matrix_free(context, &(ctx->nucleus.coord));
    if (rc != QMCKL_SUCCESS) return rc;
  }

  ctx->nucleus.coord = qmckl_matrix_alloc(context, nucl_num, 3);

  if (ctx->nucleus.coord.data == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_set_nucleus_coord",
                           NULL);
  }

  if (size_max < 3*nucl_num) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_4,
                           "qmckl_set_nucleus_coord",
                           "Array too small");
  }

  if (transp == 'N') {
    qmckl_matrix At;
    At = qmckl_matrix_alloc(context, 3, nucl_num);
    rc = qmckl_matrix_of_double(context, coord, 3*nucl_num, &At);
    if (rc != QMCKL_SUCCESS) return rc;
    rc = qmckl_transpose(context, At, ctx->nucleus.coord);
  } else {
    rc = qmckl_matrix_of_double(context, coord, nucl_num*3, &(ctx->nucleus.coord));
  }
  if (rc != QMCKL_SUCCESS) return rc;

  ctx->nucleus.uninitialized &= ~mask;
  ctx->nucleus.provided = (ctx->nucleus.uninitialized == 0);
  
  return QMCKL_SUCCESS;
}



/* Performs a deep copy of the nucleus structure from ~src~ to ~dest~. */
/* Memory allocations are done using the provided context. */


qmckl_exit_code qmckl_copy_nucleus(qmckl_context context, const qmckl_nucleus_struct* src, qmckl_nucleus_struct* dest) {
  
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  }

  if (src == NULL || dest == NULL) {
    return qmckl_failwith(context,
                          QMCKL_INVALID_ARG_2,
                          "qmckl_copy_nucleus",
                          "NULL pointer");
  }

  /* Copy scalar fields */
  dest->num = src->num;
  dest->repulsion_date = src->repulsion_date;
  dest->nn_distance_date = src->nn_distance_date;
  dest->coord_date = src->coord_date;
  dest->repulsion = src->repulsion;
  dest->uninitialized = src->uninitialized;
  dest->provided = src->provided;

  /* Deep copy charge vector */
  if (src->charge.data != NULL && src->charge.size > 0) {
    size_t size = src->charge.size * sizeof(double);
    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = size;
    dest->charge.data = (double*) qmckl_malloc(context, mem_info);
    if (dest->charge.data == NULL) {
      return QMCKL_ALLOCATION_FAILED;
    }
    memcpy(dest->charge.data, src->charge.data, size);
    dest->charge.size = src->charge.size;
  } else {
    dest->charge.data = NULL;
    dest->charge.size = 0;
  }

  /* Deep copy coord matrix */
  if (src->coord.data != NULL && src->coord.size[0] > 0) {
    size_t size = src->coord.size[0] * src->coord.size[1] * sizeof(double);
    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = size;
    dest->coord.data = (double*) qmckl_malloc(context, mem_info);
    if (dest->coord.data == NULL) {
      return QMCKL_ALLOCATION_FAILED;
    }
    memcpy(dest->coord.data, src->coord.data, size);
    dest->coord.size[0] = src->coord.size[0];
    dest->coord.size[1] = src->coord.size[1];
  } else {
    dest->coord.data = NULL;
    dest->coord.size[0] = 0;
    dest->coord.size[1] = 0;
  }

  /* Deep copy nn_distance matrix */
  if (src->nn_distance.data != NULL && src->nn_distance.size[0] > 0) {
    size_t size = src->nn_distance.size[0] * src->nn_distance.size[1] * sizeof(double);
    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = size;
    dest->nn_distance.data = (double*) qmckl_malloc(context, mem_info);
    if (dest->nn_distance.data == NULL) {
      return QMCKL_ALLOCATION_FAILED;
    }
    memcpy(dest->nn_distance.data, src->nn_distance.data, size);
    dest->nn_distance.size[0] = src->nn_distance.size[0];
    dest->nn_distance.size[1] = src->nn_distance.size[1];
  } else {
    dest->nn_distance.data = NULL;
    dest->nn_distance.size[0] = 0;
    dest->nn_distance.size[1] = 0;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_nucleus_nn_distance(qmckl_context context,
                              double* distance,
                              const int64_t size_max)
{
  /* Check input parameters */
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return (char) 0;
  }

  qmckl_exit_code rc = qmckl_provide_nn_distance(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  const int64_t sze = ctx->nucleus.num * ctx->nucleus.num;
  if (sze > size_max) {
    return qmckl_failwith(context,
                          QMCKL_INVALID_ARG_3,
                          "qmckl_get_nucleus_nn_distance",
                          "Array too small");
  }
  rc = qmckl_double_of_matrix(context, ctx->nucleus.nn_distance, distance, size_max);

  return rc;
}

qmckl_exit_code qmckl_provide_nn_distance(qmckl_context context)
{
  /* Check input parameters */
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return (char) 0;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->nucleus.provided) return QMCKL_NOT_PROVIDED;

  /* Allocate array */
  if (ctx->nucleus.nn_distance.data == NULL) {
    ctx->nucleus.nn_distance =
      qmckl_matrix_alloc(context, ctx->nucleus.num, ctx->nucleus.num);

    if (ctx->nucleus.nn_distance.data == NULL) {
      return qmckl_failwith( context,
                             QMCKL_ALLOCATION_FAILED,
                             "qmckl_nn_distance",
                             NULL);
    }
  }

  qmckl_exit_code rc =
    qmckl_compute_nn_distance(context,
                              ctx->nucleus.num,
                              ctx->nucleus.coord.data,
                              ctx->nucleus.nn_distance.data);
  if (rc != QMCKL_SUCCESS) {
    return rc;
  }

  ctx->nucleus.nn_distance_date = ctx->date;

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_get_nucleus_repulsion(qmckl_context context, double* const energy)
{
  /* Check input parameters */
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return (char) 0;
  }

  qmckl_exit_code rc = qmckl_provide_nucleus_repulsion(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  *energy = ctx->nucleus.repulsion;

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_nucleus_repulsion(qmckl_context context)
{
  /* Check input parameters */
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return (char) 0;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  qmckl_exit_code rc;

  if (!ctx->nucleus.provided) return QMCKL_NOT_PROVIDED;

  rc = qmckl_provide_nn_distance(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_compute_nucleus_repulsion(context,
                                ctx->nucleus.num,
                                ctx->nucleus.charge.data,
                                ctx->nucleus.nn_distance.data,
                                &(ctx->nucleus.repulsion));
  if (rc != QMCKL_SUCCESS) {
    return rc;
    }

  ctx->nucleus.repulsion_date = ctx->date;

  return QMCKL_SUCCESS;
}
