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
#include "qmckl_blas_private_type.h"

#include "qmckl_point_private_func.h"
#include "qmckl_memory_private_func.h"
#include "qmckl_blas_private_func.h"

qmckl_exit_code qmckl_init_point(qmckl_context context) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return false;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  memset(&(ctx->point), 0, sizeof(qmckl_point_struct));

  return QMCKL_SUCCESS;
}



/* Returns the number of points stored in the context. */


qmckl_exit_code
qmckl_get_point_num (const qmckl_context context, int64_t* const num) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  }

  if (num == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_point_num",
                           "num is a null pointer");
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  assert (ctx->point.num >= (int64_t) 0);
  *num = ctx->point.num;
  return QMCKL_SUCCESS;
}



/* Returns the point coordinates as sequences of (x,y,z). */
/* The pointer is assumed to point on a memory block of size */
/* ~size_max~ \ge ~3 * point_num~. */


qmckl_exit_code
qmckl_get_point(const qmckl_context context,
                const char transp,
                double* const coord,
                const int64_t size_max)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  }

  if (coord == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_point_coord",
                           "coord is a null pointer");
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t point_num = ctx->point.num;
  if (point_num == 0) return QMCKL_NOT_PROVIDED;

  assert (ctx->point.coord.data != NULL);

  if (size_max < 3*point_num) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_point_coord",
                           "size_max too small");
  }

  qmckl_exit_code rc;
  if (transp == 'N') {
    qmckl_matrix At = qmckl_matrix_alloc( context, 3, point_num);
    rc = qmckl_transpose( context, ctx->point.coord, At );
    if (rc != QMCKL_SUCCESS) return rc;
    rc = qmckl_double_of_matrix( context, At, coord, size_max);
    if (rc != QMCKL_SUCCESS) return rc;
    rc = qmckl_matrix_free(context, &At);
  } else {
    rc = qmckl_double_of_matrix( context, ctx->point.coord, coord, size_max);
  }
  if (rc != QMCKL_SUCCESS) return rc;

  return rc;
}



/* It copies a sequence of ~num~ points $(x,y,z)$ into the context. */


qmckl_exit_code
qmckl_set_point (qmckl_context context,
                 const char transp,
                 const int64_t num,
                 const double* coord,
                 const int64_t size_max)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  if (num <= 0) {
      return qmckl_failwith( context,
                             QMCKL_INVALID_ARG_3,
                             "qmckl_set_point",
                             "Number of points should be >0.");
  }

  if (size_max < 3*num) {
      return qmckl_failwith( context,
                             QMCKL_INVALID_ARG_4,
                             "qmckl_set_point",
                             "Array too small");
  }

  if (transp != 'N' && transp != 'T') {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_set_point",
                           "transp should be 'N' or 'T'");
  }

  if (coord == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_set_point",
                           "coord is a NULL pointer");
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  qmckl_exit_code rc;
  if (num != ctx->point.num) {

    if (ctx->point.coord.data != NULL) {
        rc = qmckl_matrix_free(context, &(ctx->point.coord));
        assert (rc == QMCKL_SUCCESS);
    }

    ctx->point.coord = qmckl_matrix_alloc(context, num, 3);
    if (ctx->point.coord.data == NULL) {
      return qmckl_failwith( context,
                             QMCKL_ALLOCATION_FAILED,
                             "qmckl_set_point",
                             NULL);
    }
  };

  ctx->point.num = num;

  if (transp == 'T') {
    double *a = ctx->point.coord.data;
#ifdef HAVE_OPENMP
    #pragma omp for
#endif
    for (int64_t i=0 ; i<3*num ; ++i) {
      a[i] = coord[i];
    }
  } else {
#ifdef HAVE_OPENMP
    #pragma omp for
#endif
    for (int64_t i=0 ; i<num ; ++i) {
      qmckl_mat(ctx->point.coord, i, 0) = coord[3*i  ];
      qmckl_mat(ctx->point.coord, i, 1) = coord[3*i+1];
      qmckl_mat(ctx->point.coord, i, 2) = coord[3*i+2];
    }
  }

  /* Increment the date of the context */
  rc = qmckl_context_touch(context);
  assert (rc == QMCKL_SUCCESS);

  return QMCKL_SUCCESS;

}



/* Performs a deep copy of the point structure from ~src~ to ~dest~. */
/* Memory allocations are done using the provided context. */


qmckl_exit_code qmckl_copy_point(qmckl_context context, const qmckl_point_struct* src, qmckl_point_struct* dest) {
  
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  }

  if (src == NULL || dest == NULL) {
    return qmckl_failwith(context,
                          QMCKL_INVALID_ARG_2,
                          "qmckl_copy_point",
                          "NULL pointer");
  }

  /* Copy scalar fields */
  dest->num = src->num;
  dest->date = src->date;

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

  return QMCKL_SUCCESS;
}
