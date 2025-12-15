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
#include "qmckl_electron_private_func.h"

qmckl_exit_code qmckl_init_electron(qmckl_context context) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return false;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  ctx->electron.uninitialized = (1 << 1) - 1;

  return QMCKL_SUCCESS;
}

bool qmckl_electron_provided(const qmckl_context context) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return false;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  return ctx->electron.provided;
}



/* To set the number of electrons, we give the number of up-spin and */
/* down-spin electrons to the context and we set the number of walkers. */


qmckl_exit_code
qmckl_set_electron_num(qmckl_context context,
                       const int64_t up_num,
                       const int64_t down_num) {
  int32_t mask = 1 << 0;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
   }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  
  if (mask != 0 && !(ctx->electron.uninitialized & mask)) {
      return qmckl_failwith( context,
                             QMCKL_ALREADY_SET,
                             "qmckl_set_electron_*",
                             NULL);
  }

  if (up_num <= 0) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_set_electron_num",
                           "up_num <= 0");
  }

  if (down_num < 0) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_set_electron_num",
                           "down_num < 0");
  }

  ctx->electron.up_num = up_num;
  ctx->electron.down_num = down_num;
  ctx->electron.num = up_num + down_num;

  ctx->electron.uninitialized &= ~mask;
  ctx->electron.provided = (ctx->electron.uninitialized == 0);
  
  return QMCKL_SUCCESS;
}




/* The following function sets the electron coordinates of all the */
/* walkers. When this is done, the pointers to the old and new sets */
/* of coordinates are swapped, and the new coordinates are */
/* overwritten. This can be done only when the data relative to */
/* electrons have been set. */

/* ~size_max~ should be equal equal or geater than ~elec_num * */
/* walker.num * 3~, to be symmetric with ~qmckl_get_electron_coord~. */

/* *Important*: changing the electron coordinates increments the date */
/* in the context. */


qmckl_exit_code
qmckl_set_electron_coord(qmckl_context context,
                         const char transp,
                         const int64_t walk_num,
                         const double* coord,
                         const int64_t size_max)
{

  int32_t mask = 0;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
   }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  
  if (mask != 0 && !(ctx->electron.uninitialized & mask)) {
      return qmckl_failwith( context,
                             QMCKL_ALREADY_SET,
                             "qmckl_set_electron_*",
                             NULL);
  }

  if (transp != 'N' && transp != 'T') {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_set_electron_coord",
                           "transp should be 'N' or 'T'");
  }

  if (walk_num <= 0) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_set_electron_coord",
                           "walk_num <= 0");
  }

  if (coord == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_4,
                           "qmckl_set_electron_coord",
                           "coord is a null pointer");
  }

  const int64_t elec_num = ctx->electron.num;

  if (elec_num == 0L) {
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "qmckl_set_electron_coord",
                           "elec_num is not set");
  }

  /* Swap pointers */
  qmckl_walker tmp         = ctx->electron.walker_old;
  ctx->electron.walker_old = ctx->electron.walker;
  ctx->electron.walker     = tmp;

  memcpy(&(ctx->point), &(ctx->electron.walker.point), sizeof(qmckl_point_struct));

  qmckl_exit_code rc;
  rc = qmckl_set_point(context, transp, walk_num*elec_num, coord, size_max);
  if (rc != QMCKL_SUCCESS) return rc;

  ctx->electron.walker.num = walk_num;
  memcpy(&(ctx->electron.walker.point), &(ctx->point), sizeof(qmckl_point_struct));

  // If it is the first time we set the electrons, we set also walkers_old.
  if (ctx->electron.walker_old.num == 0) {
    qmckl_set_electron_coord(context, transp, walk_num, coord, size_max);
  }
  return QMCKL_SUCCESS;

}

qmckl_exit_code
qmckl_get_electron_num (const qmckl_context context, int64_t* const num) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  }

  if (num == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_electron_num",
                           "num is a null pointer");
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 0;

  if ( (ctx->electron.uninitialized & mask) != 0) {
    return QMCKL_NOT_PROVIDED;
  }

  assert (ctx->electron.num > (int64_t) 0);
  *num = ctx->electron.num;
  return QMCKL_SUCCESS;
}


qmckl_exit_code
qmckl_get_electron_up_num (const qmckl_context context, int64_t* const up_num) {
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  }

  if (up_num == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_electron_up_num",
                           "up_num is a null pointer");
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 0;

  if ( (ctx->electron.uninitialized & mask) != 0) {
    return QMCKL_NOT_PROVIDED;
  }

  assert (ctx->electron.up_num > (int64_t) 0);
  *up_num = ctx->electron.up_num;
  return QMCKL_SUCCESS;
}


qmckl_exit_code
qmckl_get_electron_down_num (const qmckl_context context, int64_t* const down_num) {
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  }

  if (down_num == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_electron_down_num",
                           "down_num is a null pointer");
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 0;

  if ( (ctx->electron.uninitialized & mask) != 0) {
    return QMCKL_NOT_PROVIDED;
  }

  assert (ctx->electron.down_num >= (int64_t) 0);
  *down_num = ctx->electron.down_num;
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_electron_walk_num (const qmckl_context context, int64_t* const walk_num) {
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  }

  if (walk_num == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_electron_walk_num",
                           "walk_num is a null pointer");
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  *walk_num = ctx->electron.walker.num;
  return QMCKL_SUCCESS;
}



/* As the ~walker~ attribute is equal to ~points~, returning the */
/* current electron coordinates is equivalent to returning the */
/* current points. */


qmckl_exit_code
qmckl_get_electron_coord (const qmckl_context context,
                          const char transp,
                          double* const coord,
                          const int64_t size_max)
{
  if (transp != 'N' && transp != 'T') {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_electron_coord",
                           "transp should be 'N' or 'T'");
  }

  if (coord == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_electron_coord",
                           "coord is a null pointer");
  }

  if (size_max <= 0) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_4,
                           "qmckl_get_electron_coord",
                           "size_max should be > 0");
  }

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->electron.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_get_electron_coord",
                           NULL);
  }

  assert (ctx->point.num  == ctx->electron.walker.point.num);
  assert (ctx->point.coord.data == ctx->electron.walker.point.coord.data);

  return qmckl_get_point(context, transp, coord, size_max);
}



/* Performs a deep copy of the electron structure from ~src~ to ~dest~. */
/* Memory allocations are done using the provided context. */


qmckl_exit_code qmckl_copy_electron(qmckl_context context, const qmckl_electron_struct* src, qmckl_electron_struct* dest) {
  qmckl_exit_code rc;
  
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  }

  if (src == NULL || dest == NULL) {
    return qmckl_failwith(context,
                          QMCKL_INVALID_ARG_2,
                          "qmckl_copy_electron",
                          "NULL pointer");
  }

  /* Copy scalar fields */
  dest->num = src->num;
  dest->up_num = src->up_num;
  dest->down_num = src->down_num;
  dest->ee_distance_date = src->ee_distance_date;
  dest->en_distance_date = src->en_distance_date;
  dest->ee_potential_date = src->ee_potential_date;
  dest->en_potential_date = src->en_potential_date;
  dest->uninitialized = src->uninitialized;
  dest->provided = src->provided;

  /* Deep copy walker structures */
  dest->walker.num = src->walker.num;
  rc = qmckl_copy_point(context, &(src->walker.point), &(dest->walker.point));
  if (rc != QMCKL_SUCCESS) return rc;

  dest->walker_old.num = src->walker_old.num;
  rc = qmckl_copy_point(context, &(src->walker_old.point), &(dest->walker_old.point));
  if (rc != QMCKL_SUCCESS) return rc;

  /* Deep copy distance and potential arrays */
  /* Note: These are computed arrays and may not need copying in all cases */
  /* For now, we set them to NULL as they can be recomputed */
  /* Set dates to zero so they will be recomputed */
  dest->ee_distance = NULL;
  dest->ee_distance_date = 0;
  dest->en_distance = NULL;
  dest->en_distance_date = 0;
  dest->ee_potential = NULL;
  dest->ee_potential_date = 0;
  dest->en_potential = NULL;
  dest->en_potential_date = 0;

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_electron_ee_distance(qmckl_context context,
                               double* const distance,
                               const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_ee_distance(context);
  if (rc != QMCKL_SUCCESS) return rc;

  if (distance == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_electron_ee_distance",
                           "distance is a null pointer");
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  const int64_t sze = ctx->electron.num * ctx->electron.num * ctx->electron.walker.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_electron_ee_distance",
                           "size_max < num*num*walk_num");
  }
  memcpy(distance, ctx->electron.ee_distance, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_ee_distance(qmckl_context context)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);


  /* Compute if necessary */
  if (ctx->point.date > ctx->electron.ee_distance_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      free(ctx->electron.ee_distance);
      ctx->electron.ee_distance = NULL;
    }

    /* Allocate array */
    if (ctx->electron.ee_distance == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.num * ctx->electron.num *
        ctx->electron.walker.num * sizeof(double);
      double* ee_distance = (double*) qmckl_malloc(context, mem_info);

      if (ee_distance == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_ee_distance",
                               NULL);
      }
      ctx->electron.ee_distance = ee_distance;
    }

    qmckl_exit_code rc =
      qmckl_compute_ee_distance(context,
                                ctx->electron.num,
                                ctx->electron.walker.num,
                                ctx->electron.walker.point.coord.data,
                                ctx->electron.ee_distance);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->electron.ee_distance_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_get_electron_ee_potential(qmckl_context context, double* const ee_potential)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_ee_potential(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  size_t sze = ctx->electron.walker.num * sizeof(double);
  memcpy(ee_potential, ctx->electron.ee_potential, sze);

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_ee_potential(qmckl_context context)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->electron.provided) return QMCKL_NOT_PROVIDED;

  qmckl_exit_code rc = qmckl_provide_ee_distance(context);
  if (rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->point.date > ctx->electron.ee_potential_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      free(ctx->electron.ee_potential);
    }

    /* Allocate array */
    if (ctx->electron.ee_potential == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.walker.num * sizeof(double);
      double* ee_potential = (double*) qmckl_malloc(context, mem_info);

      if (ee_potential == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_ee_potential",
                               NULL);
      }
      ctx->electron.ee_potential = ee_potential;
    }

    rc = qmckl_compute_ee_potential(context,
                                ctx->electron.num,
                                ctx->electron.walker.num,
                                ctx->electron.ee_distance,
                                ctx->electron.ee_potential);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->electron.ee_potential_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_electron_en_distance(qmckl_context context,
                               double* const distance,
                               const int64_t size_max)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_en_distance(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (distance == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_electron_en_distance",
                           "distance is a null pointer");
  }

  int64_t sze = ctx->point.num * ctx->nucleus.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_electron_en_distance",
                           "size_max < num*nucl_num");
  }
  memcpy(distance, ctx->electron.en_distance, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_en_distance(qmckl_context context)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!(ctx->nucleus.provided)) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_provide_en_distance",
                           NULL);
  }

  /* Compute if necessary */
  if (ctx->point.date > ctx->electron.en_distance_date) {

    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = ctx->point.num * ctx->nucleus.num * sizeof(double);

    if (ctx->electron.en_distance != NULL) {
      qmckl_memory_info_struct mem_info_test = qmckl_memory_info_struct_zero;
      qmckl_exit_code rc = qmckl_get_malloc_info(context, ctx->electron.en_distance, &mem_info_test);

      /* if rc != QMCKL_SUCCESS, we are maybe in an _inplace function because the
         memory was not allocated with qmckl_malloc */

      if ((rc == QMCKL_SUCCESS) && (mem_info_test.size != mem_info.size)) {
        rc = qmckl_free(context, ctx->electron.en_distance);
        assert (rc == QMCKL_SUCCESS);
        ctx->electron.en_distance = NULL;
      }
    }

    /* Allocate array */
    if (ctx->electron.en_distance == NULL) {

      double* en_distance = (double*) qmckl_malloc(context, mem_info);

      if (en_distance == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_en_distance",
                               NULL);
      }
      ctx->electron.en_distance = en_distance;
    }

    qmckl_exit_code rc =
      qmckl_compute_en_distance(context,
                                ctx->point.num,
                                ctx->nucleus.num,
                                ctx->point.coord.data,
                                ctx->nucleus.coord.data,
                                ctx->electron.en_distance);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->electron.en_distance_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_get_electron_en_potential(qmckl_context context, double* const en_potential)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_en_potential(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  size_t sze = ctx->electron.walker.num * sizeof(double);
  memcpy(en_potential, ctx->electron.en_potential, sze);

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_en_potential(qmckl_context context)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->electron.provided) return QMCKL_NOT_PROVIDED;
  if (!ctx->nucleus.provided) return QMCKL_NOT_PROVIDED;

  qmckl_exit_code rc = qmckl_provide_en_distance(context);
  if (rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->point.date > ctx->electron.en_potential_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      free(ctx->electron.en_potential);
      ctx->electron.en_potential = NULL;
    }

    /* Allocate array */
    if (ctx->electron.en_potential == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.walker.num * sizeof(double);
      double* en_potential = (double*) qmckl_malloc(context, mem_info);

      if (en_potential == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_en_potential",
                               NULL);
      }
      ctx->electron.en_potential = en_potential;
    }

    rc = qmckl_compute_en_potential(context,
                                ctx->electron.num,
                                ctx->nucleus.num,
                                ctx->electron.walker.num,
                                ctx->nucleus.charge.data,
                                ctx->electron.en_distance,
                                ctx->electron.en_potential);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->electron.en_potential_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}
