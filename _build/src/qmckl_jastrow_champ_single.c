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
#include <sys/time.h>


#include <stdio.h>

#include "qmckl.h"
#include "qmckl_context_private_type.h"
#include "qmckl_memory_private_type.h"
#include "qmckl_memory_private_func.h"
#include "qmckl_jastrow_champ_private_type.h"
#include "qmckl_jastrow_champ_private_func.h"
#include "qmckl_jastrow_champ_single_private_type.h"
#include "qmckl_jastrow_champ_single_private_func.h"

qmckl_exit_code
qmckl_copy_jastrow_champ_single(qmckl_context context,
                                 const qmckl_jastrow_champ_single_struct* src,
                                 qmckl_jastrow_champ_single_struct* dest)
{
  /* Use memcpy to copy the entire struct (including scalars and pointer values) */
  memcpy(dest, src, sizeof(qmckl_jastrow_champ_single_struct));

  /* Deep copy the coord matrix data */
  if (src->coord.data != NULL && src->coord.size[0] > 0 && src->coord.size[1] > 0) {
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

  /* For computed arrays, set them to NULL and their dates to 0 so they can be recomputed */
  dest->een_rescaled_single_e = NULL;
  dest->een_rescaled_single_e_date = 0;
  dest->een_rescaled_single_n = NULL;
  dest->een_rescaled_single_n_date = 0;
  dest->single_ee_distance = NULL;
  dest->single_ee_distance_date = 0;
  dest->single_en_distance = NULL;
  dest->single_en_distance_date = 0;
  dest->delta_een = NULL;
  dest->delta_een_hpc_function = NULL;
  dest->delta_een_date = 0;
  dest->delta_een_pderiv = NULL;
  dest->delta_een_pderiv_date = 0;
  dest->ee_rescaled_single = NULL;
  dest->ee_rescaled_single_date = 0;
  dest->en_rescaled_single = NULL;
  dest->en_rescaled_single_date = 0;
  dest->delta_en = NULL;
  dest->delta_en_date = 0;
  dest->delta_en_pderiv = NULL;
  dest->delta_en_pderiv_date = 0;
  dest->delta_ee = NULL;
  dest->delta_ee_date = 0;
  dest->delta_ee_pderiv = NULL;
  dest->delta_ee_pderiv_date = 0;
  dest->een_rescaled_single_e_gl = NULL;
  dest->een_rescaled_single_e_gl_date = 0;
  dest->een_rescaled_single_n_gl = NULL;
  dest->een_rescaled_single_n_gl_date = 0;
  dest->delta_een_gl = NULL;
  dest->delta_een_gl_date = 0;
  dest->delta_een_g = NULL;
  dest->delta_een_g_date = 0;
  dest->ee_rescaled_single_gl = NULL;
  dest->ee_rescaled_single_gl_date = 0;
  dest->en_rescaled_single_gl = NULL;
  dest->en_rescaled_single_gl_date = 0;
  dest->delta_en_gl = NULL;
  dest->delta_en_gl_date = 0;
  dest->delta_ee_gl = NULL;
  dest->delta_ee_gl_date = 0;
  dest->tmp = NULL;
  dest->tmp_date = 0;

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_set_single_point (qmckl_context context,
                 const char transp,
                 const int64_t num,
                 const double* coord,
                 const int64_t size_max)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  if (num < 0) {
      return qmckl_failwith( context,
                             QMCKL_INVALID_ARG_3,
                             "qmckl_set_single_point",
                             "Incorrect point number");
  }

  if (transp != 'N' && transp != 'T') {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_set_single_point",
                           "transp should be 'N' or 'T'");
  }

  if (coord == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_set_single_point",
                           "coord is a NULL pointer");
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t walk_num = ctx->electron.walker.num;

  if (size_max < 3*walk_num) {
      return qmckl_failwith( context,
                             QMCKL_INVALID_ARG_4,
                             "qmckl_set_single_point",
                             "Array too small");
  }

  //qmckl_exit_code rc;

  //if (ctx->single_point.coord.data != NULL) {
  //  rc = qmckl_matrix_free(context, &(ctx->single_point.coord));
  //  assert (rc == QMCKL_SUCCESS);
  //}

  if (ctx->single_point.coord.data == NULL) {
    ctx->single_point.coord = qmckl_matrix_alloc(context, walk_num, 3);
    if (ctx->single_point.coord.data == NULL) {
      return qmckl_failwith( context,
                            QMCKL_ALLOCATION_FAILED,
                            "qmckl_set_single_point",
                            NULL);
    }
  }

  ctx->single_point.num = num;

  if (transp == 'N') {
    double *a = ctx->single_point.coord.data;
    for (int64_t i=0 ; i<3*walk_num ; ++i) {
      a[i] = coord[i];
    }
  } else {
    for (int64_t i=0 ; i<walk_num ; ++i) {
      qmckl_mat(ctx->single_point.coord, i, 0) = coord[i*walk_num + 0];
      qmckl_mat(ctx->single_point.coord, i, 1) = coord[i*walk_num + 1];
      qmckl_mat(ctx->single_point.coord, i, 2) = coord[i*walk_num + 2];
    }
  }

  /* Increment the date of the single point */
  ctx->single_point.date += 1UL;

  return QMCKL_SUCCESS;

}

qmckl_exit_code
qmckl_set_single_point_f (qmckl_context context,
                 const char transp,
                 const int64_t num,
                 const double* coord,
                 const int64_t size_max)
{
  return qmckl_set_single_point(context, transp, num-1, coord, size_max);
}

qmckl_exit_code
qmckl_single_touch(const qmckl_context context)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_single_touch",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;

  ctx->date += 1UL;
  ctx->point.date = ctx-> date;
  ctx->electron.walker.point.date = ctx-> date;
  ctx->single_point.date = ctx-> date;
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_single_electron_ee_distance(qmckl_context context,
                                      double* const distance,
                                      const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_single_ee_distance(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (distance == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_single_electron_ee_distance",
                           "distance is a NULL pointer");
  }

  int64_t sze = ctx->electron.num * ctx->electron.walker.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_single_electron_ee_distance",
                           "Array too small. Expected ctx->electron.num * ctx->electron.walker.num");
  }
  memcpy(distance, ctx->single_point.single_ee_distance, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_single_ee_distance(qmckl_context context)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  /* Compute if necessary */
  if (ctx->single_point.date > ctx->single_point.single_ee_distance_date) {

    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size =  ctx->electron.num * ctx->electron.walker.num * sizeof(double);

    if (mem_info.size > ctx->single_point.single_ee_distance_maxsize) {
      if (ctx->single_point.single_ee_distance != NULL) {
        rc = qmckl_free(context, ctx->single_point.single_ee_distance);
        assert(rc == QMCKL_SUCCESS);
        ctx->single_point.single_ee_distance = NULL;
      }
    }

    /* Allocate array */
    if (ctx->single_point.single_ee_distance == NULL) {
      double* single_ee_distance = (double*) qmckl_malloc(context, mem_info);
      ctx->single_point.single_ee_distance_maxsize = mem_info.size;
      if (single_ee_distance == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_single_ee_distance",
                               NULL);
      }
      ctx->single_point.single_ee_distance = single_ee_distance;
    }

    rc =
      qmckl_compute_single_ee_distance(context,
                                ctx->single_point.num,
                                ctx->electron.num,
                                ctx->electron.walker.num,
                                ctx->electron.walker.point.coord.data,
                                ctx->single_point.coord.data,
                                ctx->single_point.single_ee_distance);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->single_point.single_ee_distance_date = ctx->single_point.date;
  }

  //printf("single_ee_distance_date %u\n", ctx->single_point.single_ee_distance_date);

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_single_electron_en_distance(qmckl_context context,
                                      double* distance,
                                      const int64_t size_max)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_single_en_distance(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (distance == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_single_electron_en_distance",
                           "distance is a NULL pointer");
  }

  int64_t sze = ctx->nucleus.num * ctx->electron.walker.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_single_electron_en_distance",
                           "Array too small. Expected ctx->nucleus.num * ctx->electron.walker.num");
  }
  memcpy(distance, ctx->single_point.single_en_distance, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_single_en_distance(qmckl_context context)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!(ctx->nucleus.provided)) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_provide_single_en_distance",
                           NULL);
  }

  /* Compute if necessary */
  if (ctx->single_point.date > ctx->single_point.single_en_distance_date) {

    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = ctx->nucleus.num * ctx->electron.walker.num * sizeof(double);

    if (mem_info.size > ctx->single_point.single_en_distance_maxsize) {
      if (ctx->single_point.single_en_distance != NULL) {
          rc = qmckl_free(context, ctx->single_point.single_en_distance);
          assert (rc == QMCKL_SUCCESS);
          ctx->single_point.single_en_distance = NULL;
        }
    }

    /* Allocate array */
    if (ctx->single_point.single_en_distance == NULL) {

      double* single_en_distance = (double*) qmckl_malloc(context, mem_info);
      ctx->single_point.single_en_distance_maxsize = mem_info.size;

      if (single_en_distance == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_single_en_distance",
                               NULL);
      }
      ctx->single_point.single_en_distance = single_en_distance;
    }

    qmckl_exit_code rc =
      qmckl_compute_single_en_distance(context,
                                ctx->nucleus.num,
                                ctx->electron.walker.num,
                                ctx->single_point.coord.data,
                                ctx->nucleus.coord.data,
                                ctx->single_point.single_en_distance);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->single_point.single_en_distance_date = ctx->single_point.date;
  }

  // printf("single_en_distance_date %u\n", ctx->single_point.single_en_distance_date);

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_een_rescaled_single_e(qmckl_context context,
                                 double* const distance_rescaled,
                                 const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_een_rescaled_single_e(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t sze = ctx->electron.num * ctx->electron.walker.num * (ctx->jastrow_champ.cord_num + 1);
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_een_rescaled_single_e",
                           "Array too small. Expected ctx->electron.num * ctx->electron.num * ctx->electron.walker.num * (ctx->jastrow_champ.cord_num + 1)");
  }
  memcpy(distance_rescaled, ctx->single_point.een_rescaled_single_e, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_een_rescaled_single_e(qmckl_context context)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  /* Check if ee distance is provided */
  qmckl_exit_code rc = qmckl_provide_single_ee_distance(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Check if rescaled ee distance is provided */
  rc = qmckl_provide_een_rescaled_e(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->single_point.date > ctx->single_point.een_rescaled_single_e_date) {

    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = ctx->electron.num * ctx->electron.walker.num *
        (ctx->jastrow_champ.cord_num + 1) * sizeof(double);

    if (mem_info.size > ctx->single_point.een_rescaled_single_e_maxsize) {
      if (ctx->single_point.een_rescaled_single_e!= NULL) {
        rc = qmckl_free(context, ctx->single_point.een_rescaled_single_e);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_een_rescaled_single_e",
                                 "Unable to free ctx->single_point.een_rescaled_single_e");
        }
        ctx->single_point.een_rescaled_single_e = NULL;
      }
    }

    /* Allocate array */
    if (ctx->single_point.een_rescaled_single_e == NULL) {


      double* een_rescaled_single_e = (double*) qmckl_malloc(context, mem_info);
      ctx->single_point.een_rescaled_single_e_maxsize = mem_info.size;

      if (een_rescaled_single_e == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_een_rescaled_single_e",
                               NULL);
      }

      ctx->single_point.een_rescaled_single_e = een_rescaled_single_e;
    }

    rc = qmckl_compute_een_rescaled_single_e(context,
                                      ctx->single_point.num,
                                      ctx->electron.walker.num,
                                      ctx->electron.num,
                                      ctx->jastrow_champ.cord_num,
                                      ctx->jastrow_champ.rescale_factor_ee,
                                      ctx->single_point.single_ee_distance,
                                      ctx->jastrow_champ.een_rescaled_e,
                                      ctx->single_point.een_rescaled_single_e);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->single_point.een_rescaled_single_e_date = ctx->single_point.date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_compute_een_rescaled_single_e (const qmckl_context context,
                                     const int64_t num,
                                     const int64_t walk_num,
                                     const int64_t elec_num,
                                     const int64_t cord_num,
                                     const double rescale_factor_ee,
                                     const double* single_ee_distance,
                                     const double* een_rescaled_e,
                                     double* const een_rescaled_single_e )
{

#ifdef HAVE_HPC
  return qmckl_compute_een_rescaled_single_e_doc
#else
  return qmckl_compute_een_rescaled_single_e_doc
#endif
    (context, num, walk_num, elec_num, cord_num, rescale_factor_ee, single_ee_distance, een_rescaled_e, een_rescaled_single_e);
}

qmckl_exit_code
qmckl_get_een_rescaled_single_n(qmckl_context context,
                                 double* const distance_rescaled,
                                 const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_een_rescaled_single_n(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t sze = ctx->nucleus.num * ctx->electron.walker.num * (ctx->jastrow_champ.cord_num + 1);
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "mckl_get_een_rescaled_single_n",
                           "Array too small. Expected ctx->nucleus.num * ctx->electron.walker.num * (ctx->jastrow_champ.cord_num + 1)");
  }
  memcpy(distance_rescaled, ctx->single_point.een_rescaled_single_n, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_een_rescaled_single_n(qmckl_context context)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  /* Check if ee distance is provided */
  qmckl_exit_code rc = qmckl_provide_single_en_distance(context);
  if(rc != QMCKL_SUCCESS) return rc;

    /* Check if een rescaled distance is provided */
  rc = qmckl_provide_een_rescaled_n(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->single_point.date > ctx->single_point.een_rescaled_single_n_date) {

    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size =  ctx->nucleus.num *
        ctx->electron.walker.num * (ctx->jastrow_champ.cord_num + 1) * sizeof(double);

     if (mem_info.size > ctx->single_point.een_rescaled_single_n_maxsize) {
      if (ctx->single_point.een_rescaled_single_n != NULL) {
        rc = qmckl_free(context, ctx->single_point.een_rescaled_single_n);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_een_rescaled_single_n",
                                 "Unable to free ctx->single_point.een_rescaled_single_n");
        }
        ctx->single_point.een_rescaled_single_n = NULL;
      }
    }

    /* Allocate array */
    if (ctx->single_point.een_rescaled_single_n == NULL) {

      double* een_rescaled_single_n = (double*) qmckl_malloc(context, mem_info);
      ctx->single_point.een_rescaled_single_n_maxsize = mem_info.size;

      if (een_rescaled_single_n == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_een_rescaled_single_n",
                               NULL);
      }
      ctx->single_point.een_rescaled_single_n = een_rescaled_single_n;
    }

    rc = qmckl_compute_een_rescaled_single_n(context,
                                      ctx->single_point.num,
                                      ctx->electron.walker.num,
                                      ctx->electron.num,
                                      ctx->nucleus.num,
                                      ctx->jastrow_champ.type_nucl_num,
                                      ctx->jastrow_champ.type_nucl_vector,
                                      ctx->jastrow_champ.cord_num,
                                      ctx->jastrow_champ.rescale_factor_en,
                                      ctx->single_point.single_en_distance,
                                      ctx->jastrow_champ.een_rescaled_n,
                                      ctx->single_point.een_rescaled_single_n);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->single_point.een_rescaled_single_n_date = ctx->single_point.date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_jastrow_champ_single_een(qmckl_context context,
                                 double* const delta_een,
                                 const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_jastrow_champ_single_een(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t sze = ctx->electron.walker.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_single_een",
                           "Array too small. Expected ctx->electron.walker.num");
  }
  memcpy(delta_een, ctx->single_point.delta_een, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_jastrow_champ_single_een(qmckl_context context)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  qmckl_exit_code rc;

  if (ctx->jastrow_champ.cord_num > 0) {

    rc = qmckl_provide_een_rescaled_n(context);
    if(rc != QMCKL_SUCCESS) return rc;

    rc = qmckl_provide_een_rescaled_e(context);
    if(rc != QMCKL_SUCCESS) return rc;

    rc = qmckl_provide_een_rescaled_single_n(context);
    if(rc != QMCKL_SUCCESS) return rc;

    rc = qmckl_provide_een_rescaled_single_e(context);
    if(rc != QMCKL_SUCCESS) return rc;

  }


  /* Compute if necessary */
  if (ctx->single_point.date > ctx->single_point.delta_een_date) {

    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = ctx->electron.walker.num * sizeof(double);

    if (mem_info.size > ctx->single_point.delta_een_maxsize) {
      if (ctx->single_point.delta_een != NULL) {
        rc = qmckl_free(context, ctx->single_point.delta_een);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_jastrow_champ_single_een",
                                 "Unable to free ctx->single_point.delta_een");
        }
        ctx->single_point.delta_een = NULL;
      }
    }

    /* Allocate array */
    if (ctx->single_point.delta_een == NULL) {

      double* delta_een = (double*) qmckl_malloc(context, mem_info);
      ctx->single_point.delta_een_maxsize = mem_info.size;

      if (delta_een == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_jastrow_champ_single_een",
                               NULL);
      }
      ctx->single_point.delta_een = delta_een;
    }

    rc = qmckl_compute_jastrow_champ_factor_single_een(context,
                                                ctx->single_point.num,
                                                ctx->electron.walker.num,
                                                ctx->electron.num,
                                                ctx->nucleus.num,
                                                ctx->jastrow_champ.cord_num,
                                                ctx->jastrow_champ.dim_c_vector,
                                                ctx->jastrow_champ.c_vector_full,
                                                ctx->jastrow_champ.lkpm_combined_index,
                                                ctx->jastrow_champ.een_rescaled_n,
                                                ctx->jastrow_champ.een_rescaled_e,
                                                ctx->single_point.een_rescaled_single_n,
                                                ctx->single_point.een_rescaled_single_e,
                                                ctx->single_point.delta_een);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->single_point.delta_een_date = ctx->single_point.date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_single_een_hpc_c(
                                                  const qmckl_context context,
                                                  const int64_t num,
                                                  const int64_t walk_num,
                                                  const int64_t elec_num,
                                                  const int64_t nucl_num,
                                                  const int64_t cord_num,      /* cord_num in Fortran; Fortran used 0:cord_num -> cord_dim = cord_num+1 */
                                                  const int64_t dim_c_vector,
                                                  const double * restrict c_vector_full, /* shape (nucl_num, dim_c_vector), column-major */
                                                  const int64_t * restrict lkpm_combined_index, /* shape (dim_c_vector,4), column-major */
                                                  /* een_rescaled_n(elec_num, nucl_num, 0:cord_num, walk_num) */
                                                  const double * restrict een_rescaled_n,
                                                  /* een_rescaled_e(elec_num, elec_num, 0:cord_num, walk_num) */
                                                  const double * restrict een_rescaled_e,
                                                  /* een_rescaled_single_n(nucl_num, 0:cord_num, walk_num) */
                                                  const double * restrict een_rescaled_single_n,
                                                  /* een_rescaled_single_e(elec_num, 0:cord_num, walk_num) */
                                                  const double * restrict een_rescaled_single_e,
                                                  double * restrict delta_een            /* length walk_num, output */
                                                  )
{
  if (!qmckl_context_check(context))
    return qmckl_failwith(context, QMCKL_INVALID_CONTEXT, "qmckl_compute_jastrow_champ_factor_single_een_hpc_c", "NULL context");

  if (walk_num <= 0) return qmckl_failwith(context, QMCKL_INVALID_ARG_3, "qmckl_compute_jastrow_champ_factor_single_een_hpc_c", "walk_num<=0");
  if (elec_num <= 0) return qmckl_failwith(context, QMCKL_INVALID_ARG_4, "qmckl_compute_jastrow_champ_factor_single_een_hpc_c", "elec_num<=0");
  if (nucl_num <= 0) return qmckl_failwith(context, QMCKL_INVALID_ARG_5, "qmckl_compute_jastrow_champ_factor_single_een_hpc_c", "nucl_num<=0");
  if (cord_num < 0)  return qmckl_failwith(context, QMCKL_INVALID_ARG_6, "qmckl_compute_jastrow_champ_factor_single_een_hpc_c", "cord_num<0");

  if (cord_num == 0) {                            /* trivial case */
    for (int64_t nw = 0; nw < walk_num; ++nw) delta_een[nw] = 0.0;
    return QMCKL_SUCCESS;
  }

  if (!lkpm_combined_index || !c_vector_full || !een_rescaled_n || !een_rescaled_e ||
      !een_rescaled_single_n || !een_rescaled_single_e || !delta_een)
    return qmckl_failwith(context, QMCKL_FAILURE,
                          "qmckl_compute_jastrow_champ_factor_single_een_hpc_c",
                          "NULL pointer argument");

  const int64_t cord_dim = cord_num + 1;
  const int64_t stride_e = elec_num + 1;
  const int64_t t_count = 2;
  const size_t total_nonzeros_elems = stride_e * nucl_num * t_count;
  int64_t * restrict non_zeros = (int64_t *) malloc(total_nonzeros_elems * sizeof(int64_t));
  if (!non_zeros) return qmckl_failwith(context, -20, "qmckl_compute_jastrow_champ_factor_single_een_hpc_c", "alloc non_zeros");
  memset(non_zeros, 0, total_nonzeros_elems * sizeof(int64_t));

  /* helper macros for indexing Fortran column-major arrays (0-based in C):
     Fortran indexing order: A(i,j,k,...) in memory => linear index:
     i + j*dim_i + k*(dim_i*dim_j) + ...
     We provide macros for each array shape to avoid mistakes. */

#define IDX_CVECTOR(a,n)  ((a) + (n) * nucl_num)
  /* lkpm_combined_index: shape (dim_c_vector, 4) */
#define IDX_LKPM(n,col)   ((n) + (col) * dim_c_vector)
  /* een_rescaled_n: (elec_num, nucl_num, cord_dim, walk_num) */
#define IDX_EEN_N(j,a,l,nw) ((j) + (a) * elec_num +                     \
                             (l) * elec_num * nucl_num +                \
                             (nw) * elec_num * nucl_num * cord_dim )
  /* een_rescaled_e: (elec_num, elec_num, cord_dim, walk_num) */
#define IDX_EEN_E(j,col,l,nw) (                                         \
                                                                        (j) + (col) * elec_num + \
                                                                        (l) * elec_num * elec_num + \
                                                                        (nw) * elec_num * elec_num * cord_dim )
  /* een_rescaled_single_n: (nucl_num, cord_dim, walk_num) */
#define IDX_EEN_SING_N(a,l,nw) (                                        \
                                                                        (a) + (l) * nucl_num + \
                                                                        (nw) * nucl_num * cord_dim )
  /* een_rescaled_single_e: (elec_num, cord_dim, walk_num) */
#define IDX_EEN_SING_E(j,l,nw) (                                        \
                                                                        (j) + (l) * elec_num + \
                                                                        (nw) * elec_num * cord_dim )

  const double eps = qmckl_get_numprec_epsilon(context);

  /* Allocate per-walker temporary non_zeros array on heap (thread-private). */
  int64_t * restrict __attribute__((aligned(64))) nz = (int64_t *) malloc(total_nonzeros_elems * sizeof(int64_t));
  assert (nz != NULL);
  /* zero only the needed prefix (we will fill counts and entries) */
  memset(nz, 0, total_nonzeros_elems * sizeof(int64_t));

  for (int64_t nw = 0; nw < walk_num; ++nw) {

    /* Build the non_zeros arrays: for each nucleus a, scan electrons j */
    for (int64_t a = 0; a < nucl_num; ++a) {
      int64_t k = 0;
      int64_t p = 0;
      const double * restrict n2_ = &(een_rescaled_n[IDX_EEN_N(0, a, 2, nw)]);
      const double * restrict n1_ = &(een_rescaled_n[IDX_EEN_N(0, a, 1, nw)]);
      int64_t * const restrict nzk_ = &(nz[a*stride_e]);
      int64_t * const restrict nzp_ = &(nz[stride_e*(a+nucl_num)]);
      for (int64_t j = 0; j < elec_num; ++j) {
        const double v_m = n2_[j];
        if (v_m > eps) {
          k++; p++;
          nzk_[k] = j;
          nzp_[p] = j;
        } else {
          const double v_lm = n1_[j];
          if (v_lm > eps) {
            k++;
            nzk_[k] = j;
          }
        }
      }
      nzk_[0] = k;
      nzp_[0] = p;
    }

    double delta_val = 0.0;

    /* Main accumulation loop: n = 1 .. dim_c_vector in Fortran -> 0..dim_c_vector-1 in C */
    for (int64_t n = 0; n < dim_c_vector; ++n) {
      /* In Fortran: l = lkpm_combined_index(n,1); k = lkpm_combined_index(n,2); m = lkpm_combined_index(n,4) */
      /* In C: columns are 0..3; we read values exactly (these are offsets used with 0:cord_num) */
      const int64_t l = lkpm_combined_index[ IDX_LKPM(n, 0) ]; /* column 1 in Fortran is index 0 here */
      const int64_t k = lkpm_combined_index[ IDX_LKPM(n, 1) ];
      const int64_t m = lkpm_combined_index[ IDX_LKPM(n, 3) ];

      /* c_vector_full(a,n) in Fortran => index a + n*nucl_num */
      for (int64_t a = 0; a < nucl_num; ++a) {
        const double cn = c_vector_full[ IDX_CVECTOR(a, n) ];
        if (cn == 0.0) continue;

        double accu = 0.0;


        const double s_n_lpm = een_rescaled_single_n[ IDX_EEN_SING_N(a, l + m, nw) ];
        const double n_num_lpm = een_rescaled_n[ IDX_EEN_N( num, a, l + m, nw) ];
        const double * restrict es_ = &(een_rescaled_single_e[ IDX_EEN_SING_E(0, k, nw) ]);
        const double * restrict e_  = &(een_rescaled_e[ IDX_EEN_E(0, num, k, nw) ]);

        {
          const int64_t idx = (m < 2) ? m-1 : 1;
          const int64_t * restrict nz_ = &(nz[ a * stride_e + idx * stride_e * nucl_num ]);
          const double * restrict n_  = &(een_rescaled_n[ IDX_EEN_N(0, a, m, nw) ]);

          if ( s_n_lpm > eps && n_num_lpm > eps ) {
            if (idx >= 0) {
              for (int64_t p = 1; p <= nz_[0] ; ++p) {
                int64_t j = nz_[p];
                accu += (s_n_lpm * es_[j] - n_num_lpm * e_[j]) *n_[j];
              }
            } else {
              for (int64_t j = 0; j < elec_num; ++j) {
                double term = s_n_lpm * es_[j] - n_num_lpm * e_[j];
                accu += term;
              }
            }
          } else if ( s_n_lpm > eps ) {
            if (idx >= 0) {
              for (int64_t p = 1; p <= nz_[0]; ++p) {
                int64_t j = nz_[p];
                accu += s_n_lpm * es_[j] * n_[j];
              }
            } else {
              for (int64_t j = 0; j < elec_num; ++j) {
                accu += s_n_lpm * es_[j];
              }
            }
          } else if (n_num_lpm > eps ) {
            if (idx >= 0) {
              for (int64_t p = 1; p <= nz_[0]; ++p) {
                int64_t j = nz_[p];
                accu -= n_num_lpm * e_[j] * n_[j];
              }
            } else {
#ifdef HAVE_OPENMP
#pragma omp simd reduction(+:accu)
#endif
              for (int64_t j = 0; j < elec_num; ++j) {
                accu -= n_num_lpm * e_[j];
              }
            }
          }
        }

        {
          int64_t idx = ((l + m) < 2) ? (l + m)-1 : 1;
          const double s_n_m = een_rescaled_single_n[ IDX_EEN_SING_N(a, m, nw) ];
          const double n_num_m = een_rescaled_n[ IDX_EEN_N( num, a, m, nw) ];
          const int64_t * restrict nz_ = &(nz[ a * stride_e + idx * stride_e * nucl_num ]);
          const double * restrict n_  = &(een_rescaled_n[ IDX_EEN_N(0, a, l+m, nw) ]);

          if ( s_n_m > eps && n_num_m > eps ) {

            double accu1 = 0.;
            double accu2 = 0.;
            for (int64_t p = 1; p <= nz_[0]; ++p) {
              int64_t j = nz_[p];
              accu1 += es_[j] * n_[j];
              accu2 += e_[j] * n_[j];
//              accu += (s_n_m * es_[j] - n_num_m * e_[j]) * n_[j];
            }
            accu += s_n_m * accu1 - n_num_m*accu2;
          } else if ( s_n_m > eps ) {
            for (int64_t p = 1; p <= nz_[0]; ++p) {
              int64_t j = nz_[p];
              accu += s_n_m * es_[j] * n_[j];
            }
          } else if ( n_num_m > eps ) {
            for (int64_t p = 1; p <= nz_[0]; ++p) {
              int64_t j = nz_[p];
              accu -= n_num_m * e_[j] * n_[j];
            }
          }

          delta_val += accu * cn;
        } /* end loop a (nuclei) */
      }
    } /* end loop n (dim_c_vector) */

    delta_een[nw] = delta_val;
  } /* end parallel over nw */
  free(nz);

  free(non_zeros);
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_single_een (qmckl_context context,
                          const int64_t num,
                          const int64_t walk_num,
                          const int64_t elec_num,
                          const int64_t nucl_num,
                          const int64_t cord_num,
                          const int64_t dim_c_vector,
                          const double* c_vector_full,
                          const int64_t* lkpm_combined_index,
                          const double* een_rescaled_n,
                          const double* een_rescaled_e,
                          const double* een_rescaled_single_n,
                          const double* een_rescaled_single_e,
                          double* const delta_een )
{
if (cord_num == 0) {
    return qmckl_compute_jastrow_champ_factor_single_een_doc(context, num, walk_num, elec_num, nucl_num, cord_num, dim_c_vector,
     c_vector_full, lkpm_combined_index, een_rescaled_n, een_rescaled_e, een_rescaled_single_n, een_rescaled_single_e, delta_een );
 }
#ifdef HAVE_HPC

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  qmckl_exit_code (*fptr) (qmckl_context context,
                          const int64_t num,
                          const int64_t walk_num,
                          const int64_t elec_num,
                          const int64_t nucl_num,
                          const int64_t cord_num,
                          const int64_t dim_c_vector,
                          const double* c_vector_full,
                          const int64_t* lkpm_combined_index,
                          const double* een_rescaled_n,
                          const double* een_rescaled_e,
                          const double* een_rescaled_single_n,
                          const double* een_rescaled_single_e,
                          double* const delta_een ) = ctx->single_point.delta_een_hpc_function ;
  if (fptr == NULL) {
    long t1, t2, t3;

    long start, end;

    {
      start = clock();
      for (int i=0 ; i<1000 ; ++i) {
        qmckl_compute_jastrow_champ_factor_single_een_doc(context, num, walk_num, elec_num, nucl_num, cord_num, dim_c_vector,
                                                          c_vector_full, lkpm_combined_index, een_rescaled_n, een_rescaled_e,
                                                          een_rescaled_single_n, een_rescaled_single_e, delta_een );
      }
      end = clock();
      t1 = end-start;
    }
    {
      start = clock();
      for (int i=0 ; i<1000 ; ++i) {
        qmckl_compute_jastrow_champ_factor_single_een_hpc_c(context, num, walk_num, elec_num, nucl_num, cord_num, dim_c_vector,
                                                          c_vector_full, lkpm_combined_index, een_rescaled_n, een_rescaled_e,
                                                          een_rescaled_single_n, een_rescaled_single_e, delta_een );
      }
      end = clock();
      t2 = end-start;
    }
    {
      start = clock();
      for (int i=0 ; i<1000 ; ++i) {
        qmckl_compute_jastrow_champ_factor_single_een_hpc_f(context, num, walk_num, elec_num, nucl_num, cord_num, dim_c_vector,
                                                          c_vector_full, lkpm_combined_index, een_rescaled_n, een_rescaled_e,
                                                          een_rescaled_single_n, een_rescaled_single_e, delta_een );
      }
      end = clock();
      t3 = end-start;
    }

    if (t1 < t2 && t1 < t3) {
      fptr = qmckl_compute_jastrow_champ_factor_single_een_doc;
    } else if (t2 < t1 && t2 < t3) {
      fptr = qmckl_compute_jastrow_champ_factor_single_een_hpc_c;
    } else {
      fptr = qmckl_compute_jastrow_champ_factor_single_een_hpc_f;
    }

    /*
    printf("doc: %ld\nhpc_c: %ld\n hpc_f: %ld\n", t1, t2, t3);
    if (t1 < t2 && t1 < t3) {
      printf("Chose doc\n");
    } else if (t2 < t1 && t2 < t3) {
      printf("Chose hpc_c\n");
    } else if (t3 < t1 && t3 < t2) {
      printf("Chose hpc_f\n");
    }
    */

    ctx->single_point.delta_een_hpc_function = (void*) fptr;
  }
  return fptr
#else
  return qmckl_compute_jastrow_champ_factor_single_een_doc
#endif
    (context, num, walk_num, elec_num, nucl_num, cord_num, dim_c_vector,
     c_vector_full, lkpm_combined_index, een_rescaled_n, een_rescaled_e, een_rescaled_single_n, een_rescaled_single_e, delta_een );

}

qmckl_exit_code
qmckl_get_een_rescaled_single_n_gl(qmckl_context context,
                                         double* const distance_rescaled,
                                         const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_een_rescaled_single_n_gl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t sze = 4 * ctx->nucleus.num * ctx->electron.walker.num * (ctx->jastrow_champ.cord_num + 1);
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_single_een_gl",
                           "Array too small. Expected 4 * ctx->nucleus.num * ctx->electron.walker.num * (ctx->jastrow_champ.cord_num + 1)");
  }
  memcpy(distance_rescaled, ctx->single_point.een_rescaled_single_n_gl, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_een_rescaled_single_n_gl(qmckl_context context)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  /* Check if en distance is provided */
  qmckl_exit_code rc = qmckl_provide_single_en_distance(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Check if ee distance is provided */
  rc = qmckl_provide_een_rescaled_single_n(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->single_point.date > ctx->single_point.een_rescaled_single_n_gl_date) {

    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size =  4 * ctx->nucleus.num * ctx->electron.walker.num * (ctx->jastrow_champ.cord_num + 1) * sizeof(double);

    if (mem_info.size > ctx->single_point.een_rescaled_single_n_gl_maxsize) {
      if (ctx->single_point.een_rescaled_single_n_gl != NULL) {
        rc = qmckl_free(context, ctx->single_point.een_rescaled_single_n_gl);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_een_rescaled_single_n_gl",
                                 "Unable to free ctx->single_pont.een_rescaled_single_n_gl");
        }
        ctx->single_point.een_rescaled_single_n_gl = NULL;
      }
    }

    /* Allocate array */
    if (ctx->single_point.een_rescaled_single_n_gl == NULL) {

      double* een_rescaled_single_n_gl = (double*) qmckl_malloc(context, mem_info);
      ctx->single_point.een_rescaled_single_n_gl_maxsize = mem_info.size;

      if (een_rescaled_single_n_gl == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_een_rescaled_single_n_gl",
                               NULL);
      }
      ctx->single_point.een_rescaled_single_n_gl = een_rescaled_single_n_gl;
    }

    rc = qmckl_compute_een_rescaled_single_n_gl(context,
                                                     ctx->electron.walker.num,
                                                     ctx->nucleus.num,
                                                     ctx->jastrow_champ.type_nucl_num,
                                                     ctx->jastrow_champ.type_nucl_vector,
                                                     ctx->jastrow_champ.cord_num,
                                                     ctx->jastrow_champ.rescale_factor_en,
                                                     ctx->single_point.coord.data,
                                                     ctx->nucleus.coord.data,
                                                     ctx->single_point.single_en_distance,
                                                     ctx->single_point.een_rescaled_single_n,
                                                     ctx->single_point.een_rescaled_single_n_gl);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->single_point.een_rescaled_single_n_gl_date = ctx->single_point.date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_een_rescaled_single_e_gl(qmckl_context context,
                                         double* const distance_rescaled,
                                         const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_een_rescaled_single_e_gl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t sze = 4 * ctx->electron.num * ctx->electron.walker.num * (ctx->jastrow_champ.cord_num + 1);
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_factor_een_gl",
                           "Array too small. Expected 4 * ctx->electron.num * ctx->electron.walker.num * (ctx->jastrow_champ.cord_num + 1)");
  }
  memcpy(distance_rescaled, ctx->single_point.een_rescaled_single_e_gl, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_een_rescaled_single_e_gl(qmckl_context context)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  /* Check if rescaled een-ee distance is provided */
  qmckl_exit_code rc = qmckl_provide_een_rescaled_single_e(context);
  if(rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_single_ee_distance(context);
  if(rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_single_en_distance(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->single_point.date > ctx->single_point.een_rescaled_single_e_gl_date) {

    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size =  4 * ctx->electron.num * ctx->electron.walker.num * (ctx->jastrow_champ.cord_num + 1) * sizeof(double);

    if (mem_info.size > ctx->single_point.een_rescaled_single_e_gl_maxsize) {
      if (ctx->single_point.een_rescaled_single_e_gl != NULL) {
        rc = qmckl_free(context, ctx->single_point.een_rescaled_single_e_gl);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_een_rescaled_e_gl",
                                 "Unable to free ctx->single_point.een_rescaled_single_e_gl");
        }
        ctx->single_point.een_rescaled_single_e_gl = NULL;
      }
    }

    /* Allocate array */
    if (ctx->single_point.een_rescaled_single_e_gl == NULL) {

      double* een_rescaled_single_e_gl = (double*) qmckl_malloc(context, mem_info);
      ctx->single_point.een_rescaled_single_e_gl_maxsize = mem_info.size;

      if (een_rescaled_single_e_gl == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_een_rescaled_single_e_gl",
                               NULL);
      }
      ctx->single_point.een_rescaled_single_e_gl = een_rescaled_single_e_gl;
    }

    rc = qmckl_compute_een_rescaled_single_e_gl(context,
                                                     ctx->single_point.num,
                                                     ctx->electron.walker.num,
                                                     ctx->electron.num,
                                                     ctx->jastrow_champ.cord_num,
                                                     ctx->jastrow_champ.rescale_factor_ee,
                                                     ctx->single_point.coord.data,
                                                     ctx->electron.walker.point.coord.data,
                                                     ctx->single_point.single_ee_distance,
                                                     ctx->single_point.een_rescaled_single_e,
                                                     ctx->single_point.een_rescaled_single_e_gl);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->single_point.een_rescaled_single_e_gl_date = ctx->single_point.date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_compute_een_rescaled_single_e_gl (
         const qmckl_context context,
          const int64_t num,
         const int64_t walk_num,
         const int64_t elec_num,
         const int64_t cord_num,
         const double rescale_factor_ee,
         const double* coord,
         const double* coord_ee,
         const double* single_ee_distance,
         const double* een_rescaled_single_e,
         double* const een_rescaled_single_e_gl )
{
#ifdef HAVE_HPC
  return qmckl_compute_een_rescaled_single_e_gl_doc
#else
  return qmckl_compute_een_rescaled_single_e_gl_doc
#endif
    (context, num, walk_num, elec_num, cord_num, rescale_factor_ee, coord,
    coord_ee, single_ee_distance, een_rescaled_single_e, een_rescaled_single_e_gl );
}

qmckl_exit_code
qmckl_get_jastrow_champ_single_tmp(qmckl_context context,
                                 double* const tmp,
                                 const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_jastrow_champ_single_tmp(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (tmp == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_champ_single_tmp",
                           "Array is NULL.");
  }

  int64_t sze = 2 * ctx->electron.walker.num * ctx->nucleus.num * (ctx->jastrow_champ.cord_num + 1) * (ctx->jastrow_champ.cord_num + 1);
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_single_tmp",
                           "Array too small. Expected 2 * ctx->electron.walker.num * ctx->nucleus.num * (ctx->jastrow_champ.cord_num + 1) * (ctx->jastrow_champ.cord_num + 1)");
  }
  memcpy(tmp, ctx->single_point.tmp, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_jastrow_champ_single_tmp(qmckl_context context)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  qmckl_exit_code rc;

  if (ctx->jastrow_champ.cord_num > 0) {

    rc = qmckl_provide_een_rescaled_n(context);
    if(rc != QMCKL_SUCCESS) return rc;

    rc = qmckl_provide_een_rescaled_e(context);
    if(rc != QMCKL_SUCCESS) return rc;

    rc = qmckl_provide_een_rescaled_single_e(context);
    if(rc != QMCKL_SUCCESS) return rc;

  }


  /* Compute if necessary */
  if (ctx->single_point.date > ctx->single_point.tmp_date) {

    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = 2 * ctx->electron.walker.num * ctx->nucleus.num * (ctx->jastrow_champ.cord_num + 1) * (ctx->jastrow_champ.cord_num + 1) * sizeof(double);

    if (mem_info.size > ctx->single_point.tmp_maxsize) {
      if (ctx->single_point.tmp != NULL) {
        rc = qmckl_free(context, ctx->single_point.tmp);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_jastrow_champ_single_tmp",
                                 "Unable to free ctx->single_point.tmp");
        }
        ctx->single_point.tmp = NULL;
      }
    }

    /* Allocate array */
    if (ctx->single_point.tmp == NULL) {

      double* tmp = (double*) qmckl_malloc(context, mem_info);
      ctx->single_point.tmp_maxsize = mem_info.size;

      if (tmp == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_jastrow_champ_single_tmp",
                               NULL);
      }
      ctx->single_point.tmp = tmp;
    }

    rc = qmckl_compute_jastrow_champ_factor_single_tmp_doc(context,
                                                ctx->single_point.num,
                                                ctx->electron.walker.num,
                                                ctx->electron.num,
                                                ctx->nucleus.num,
                                                ctx->jastrow_champ.cord_num,
                                                ctx->jastrow_champ.dim_c_vector,
                                                ctx->jastrow_champ.c_vector_full,
                                                ctx->jastrow_champ.lkpm_combined_index,
                                                ctx->jastrow_champ.een_rescaled_n,
                                                ctx->jastrow_champ.een_rescaled_e,
                                                ctx->single_point.een_rescaled_single_e,
                                                ctx->single_point.tmp);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->single_point.tmp_date = ctx->single_point.date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_jastrow_champ_single_een_gl(qmckl_context context,
                                 double* const delta_een_gl,
                                 const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_jastrow_champ_single_een_gl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (delta_een_gl == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_champ_single_een_gl",
                           "Array is NULL.");
  }

  int64_t sze = 4 * ctx->electron.num * ctx->electron.walker.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_single_een_gl",
                           "Array too small. Expected 4 * ctx->electron.num * ctx->electron.walker.num");
  }
  memcpy(delta_een_gl, ctx->single_point.delta_een_gl, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_jastrow_champ_single_een_gl(qmckl_context context)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  qmckl_exit_code rc;

  if (ctx->jastrow_champ.cord_num > 0) {

    rc = qmckl_provide_een_rescaled_n(context);
    if(rc != QMCKL_SUCCESS) return rc;

    rc = qmckl_provide_een_rescaled_single_n(context);
    if(rc != QMCKL_SUCCESS) return rc;

    rc = qmckl_provide_een_rescaled_n_gl(context);
    if(rc != QMCKL_SUCCESS) return rc;

    rc = qmckl_provide_een_rescaled_single_n_gl(context);
    if(rc != QMCKL_SUCCESS) return rc;

    rc = qmckl_provide_een_rescaled_e(context);
    if(rc != QMCKL_SUCCESS) return rc;

    rc = qmckl_provide_een_rescaled_single_e(context);
    if(rc != QMCKL_SUCCESS) return rc;

      rc = qmckl_provide_een_rescaled_e_gl(context);
    if(rc != QMCKL_SUCCESS) return rc;

    rc = qmckl_provide_een_rescaled_single_e_gl(context);
    if(rc != QMCKL_SUCCESS) return rc;

  }


  /* Compute if necessary */
  if (ctx->single_point.date > ctx->single_point.delta_een_gl_date) {

    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = 4 * ctx->electron.num * ctx->electron.walker.num * sizeof(double);

    if (mem_info.size > ctx->single_point.delta_een_gl_maxsize) {
      if (ctx->single_point.delta_een_gl != NULL) {
        rc = qmckl_free(context, ctx->single_point.delta_een_gl);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_jastrow_champ_single_een_gl",
                                 "Unable to free ctx->single_point.delta_een_gl");
        }
        ctx->single_point.delta_een_gl = NULL;
      }
    }

    /* Allocate array */
    if (ctx->single_point.delta_een_gl == NULL) {

      double* delta_een_gl = (double*) qmckl_malloc(context, mem_info);
      ctx->single_point.delta_een_gl_maxsize = mem_info.size;

      if (delta_een_gl == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_jastrow_champ_single_een_gl",
                               NULL);
      }
      ctx->single_point.delta_een_gl = delta_een_gl;
    }

    rc = qmckl_compute_jastrow_champ_factor_single_een_gl(context,
                                                ctx->single_point.num,
                                                ctx->electron.walker.num,
                                                ctx->electron.num,
                                                ctx->nucleus.num,
                                                ctx->jastrow_champ.cord_num,
                                                ctx->jastrow_champ.dim_c_vector,
                                                ctx->jastrow_champ.c_vector_full,
                                                ctx->jastrow_champ.lkpm_combined_index,
                                                ctx->jastrow_champ.een_rescaled_n,
                                                ctx->single_point.een_rescaled_single_n,
                                                ctx->jastrow_champ.een_rescaled_n_gl,
                                                ctx->single_point.een_rescaled_single_n_gl,
                                                ctx->jastrow_champ.een_rescaled_e,
                                                ctx->single_point.een_rescaled_single_e,
                                                ctx->jastrow_champ.een_rescaled_e_gl,
                                                ctx->single_point.een_rescaled_single_e_gl,
                                                ctx->single_point.delta_een_gl);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->single_point.delta_een_gl_date = ctx->single_point.date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_single_een_gl (const qmckl_context context,
                              const int64_t num,
                              const int64_t walk_num,
                              const int64_t elec_num,
                              const int64_t nucl_num,
                              const int64_t cord_num,
                              const int64_t dim_c_vector,
                              const double* c_vector_full,
                              const int64_t* lkpm_combined_index,
                              const double* een_rescaled_n,
                              const double* een_rescaled_single_n,
                              const double* een_rescaled_n_gl,
                              const double* een_rescaled_single_n_gl,
                              const double* een_rescaled_e,
                              const double* een_rescaled_single_e,
                              const double* een_rescaled_e_gl,
                              const double* een_rescaled_single_e_gl,
                          double* const delta_een_gl )
{
#ifdef HAVE_HPC
  return qmckl_compute_jastrow_champ_factor_single_een_gl_hpc
#else
  return qmckl_compute_jastrow_champ_factor_single_een_gl_doc
#endif
    (context, num, walk_num, elec_num, nucl_num, cord_num, dim_c_vector,
     c_vector_full, lkpm_combined_index, een_rescaled_n, een_rescaled_single_n, een_rescaled_n_gl, een_rescaled_single_n_gl, een_rescaled_e, een_rescaled_single_e, een_rescaled_e_gl, een_rescaled_single_e_gl, delta_een_gl );

}

qmckl_exit_code
qmckl_get_jastrow_champ_single_een_g(qmckl_context context,
                                 double* const delta_een_g,
                                 const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_jastrow_champ_single_een_g(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (delta_een_g == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_jastrow_champ_single_een_g",
                           "Array is NULL.");
  }

  int64_t sze = 4 * ctx->electron.num * ctx->electron.walker.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_single_een_g",
                           "Array too small. Expected 4 * ctx->electron.num * ctx->electron.walker.num");
  }
  memcpy(delta_een_g, ctx->single_point.delta_een_g, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_jastrow_champ_single_een_g(qmckl_context context)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  qmckl_exit_code rc;

  if (ctx->jastrow_champ.cord_num > 0) {

    rc = qmckl_provide_een_rescaled_n(context);
    if(rc != QMCKL_SUCCESS) return rc;

    rc = qmckl_provide_een_rescaled_single_n(context);
    if(rc != QMCKL_SUCCESS) return rc;

    rc = qmckl_provide_een_rescaled_n_gl(context);
    if(rc != QMCKL_SUCCESS) return rc;

    rc = qmckl_provide_een_rescaled_single_n_gl(context);
    if(rc != QMCKL_SUCCESS) return rc;

    rc = qmckl_provide_een_rescaled_e(context);
    if(rc != QMCKL_SUCCESS) return rc;

    rc = qmckl_provide_een_rescaled_single_e(context);
    if(rc != QMCKL_SUCCESS) return rc;

    rc = qmckl_provide_een_rescaled_e_gl(context);
    if(rc != QMCKL_SUCCESS) return rc;

    rc = qmckl_provide_een_rescaled_single_e_gl(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if single tmp matrix is provided */
    rc = qmckl_provide_jastrow_champ_single_tmp(context);
    if(rc != QMCKL_SUCCESS) return rc;

  }


  /* Compute if necessary */
  if (ctx->single_point.date > ctx->single_point.delta_een_g_date) {

    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = 4 * ctx->electron.num * ctx->electron.walker.num * sizeof(double);

    if (mem_info.size > ctx->single_point.delta_een_g_maxsize) {
      if (ctx->single_point.delta_een_g != NULL) {
        rc = qmckl_free(context, ctx->single_point.delta_een_g);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_jastrow_champ_single_een_g",
                                 "Unable to free ctx->single_point.delta_een_g");
        }
        ctx->single_point.delta_een_g = NULL;
      }
    }

    /* Allocate array */
    if (ctx->single_point.delta_een_g == NULL) {

      double* delta_een_g = (double*) qmckl_malloc(context, mem_info);
      ctx->single_point.delta_een_g_maxsize = mem_info.size;

      if (delta_een_g == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_jastrow_champ_single_een_g",
                               NULL);
      }
      ctx->single_point.delta_een_g = delta_een_g;
    }

    rc = qmckl_compute_jastrow_champ_factor_single_een_g_doc(context,
                                                ctx->single_point.num,
                                                ctx->electron.walker.num,
                                                ctx->electron.num,
                                                ctx->nucleus.num,
                                                ctx->jastrow_champ.cord_num,
                                                ctx->jastrow_champ.dim_c_vector,
                                                ctx->jastrow_champ.c_vector_full,
                                                ctx->jastrow_champ.lkpm_combined_index,
                                                ctx->jastrow_champ.een_rescaled_n,
                                                ctx->single_point.een_rescaled_single_n,
                                                ctx->jastrow_champ.een_rescaled_n_gl,
                                                ctx->single_point.een_rescaled_single_n_gl,
                                                ctx->jastrow_champ.een_rescaled_e,
                                                ctx->single_point.een_rescaled_single_e,
                                                ctx->jastrow_champ.een_rescaled_e_gl,
                                                ctx->single_point.een_rescaled_single_e_gl,
                                                ctx->single_point.tmp,
                                                ctx->single_point.delta_een_g);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    //ctx->single_point.delta_een_g_date = ctx->single_point.date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_single_een_g (const qmckl_context context,
                              const int64_t num,
                              const int64_t walk_num,
                              const int64_t elec_num,
                              const int64_t nucl_num,
                              const int64_t cord_num,
                              const int64_t dim_c_vector,
                              const double* c_vector_full,
                              const int64_t* lkpm_combined_index,
                              const double* een_rescaled_n,
                              const double* een_rescaled_single_n,
                              const double* een_rescaled_n_gl,
                              const double* een_rescaled_single_n_gl,
                              const double* een_rescaled_e,
                              const double* een_rescaled_single_e,
                              const double* een_rescaled_e_gl,
                              const double* een_rescaled_single_e_gl,
                              const double* tmp,
                          double* const delta_een_g )
{
#ifdef HAVE_HPC
  return qmckl_compute_jastrow_champ_factor_single_een_g_doc
#else
  return qmckl_compute_jastrow_champ_factor_single_een_g_doc
#endif
    (context, num, walk_num, elec_num, nucl_num, cord_num, dim_c_vector,
     c_vector_full, lkpm_combined_index, een_rescaled_n, een_rescaled_single_n, een_rescaled_n_gl,
     een_rescaled_single_n_gl, een_rescaled_e, een_rescaled_single_e, een_rescaled_e_gl, een_rescaled_single_e_gl, tmp, delta_een_g );

}

qmckl_exit_code
qmckl_get_jastrow_champ_single_een_pderiv(qmckl_context context,
                                 double* const delta_een_pderiv,
                                 const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_jastrow_champ_single_een_pderiv(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t sze = ctx->jastrow_champ.dim_c_vector * ctx->jastrow_champ.type_nucl_num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_single_een_pderiv",
                           "Array too small. Expected ctx->jastrow_champ.dim_c_vector * ctx->jastrow_champ.type_nucl_num");
  }
  memcpy(delta_een_pderiv, ctx->single_point.delta_een_pderiv, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_jastrow_champ_single_een_pderiv(qmckl_context context)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  qmckl_exit_code rc;

  if (ctx->jastrow_champ.cord_num == 0) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_provide_jastrow_champ_factor_een_pderiv",
                           "Expected cord_num > 0");
  }
  if (ctx->jastrow_champ.cord_num > 0) {

    rc = qmckl_provide_een_rescaled_n(context);
    if(rc != QMCKL_SUCCESS) return rc;

    rc = qmckl_provide_een_rescaled_single_n(context);
    if(rc != QMCKL_SUCCESS) return rc;

    rc = qmckl_provide_een_rescaled_e(context);
    if(rc != QMCKL_SUCCESS) return rc;

    rc = qmckl_provide_een_rescaled_single_e(context);
    if(rc != QMCKL_SUCCESS) return rc;
  }


  /* Compute if necessary */
  if (ctx->single_point.date > ctx->single_point.delta_een_pderiv_date) {

    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = ctx->jastrow_champ.type_nucl_num * ctx->jastrow_champ.dim_c_vector * sizeof(double);


    /* Allocate array */
    if (ctx->single_point.delta_een_pderiv == NULL) {

      double* delta_een_pderiv = (double*) qmckl_malloc(context, mem_info);
      ctx->single_point.delta_een_pderiv_maxsize = mem_info.size;

      if (delta_een_pderiv == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_jastrow_champ_single_een_pderiv",
                               NULL);
      }
      ctx->single_point.delta_een_pderiv = delta_een_pderiv;
    }

    rc = qmckl_compute_jastrow_champ_factor_single_een_pderiv_doc(context,
                                                ctx->single_point.num,
                                                ctx->electron.walker.num,
                                                ctx->electron.num,
                                                ctx->nucleus.num,
                                                ctx->jastrow_champ.type_nucl_num,
                                                ctx->jastrow_champ.type_nucl_vector,
                                                ctx->jastrow_champ.cord_num,
                                                ctx->jastrow_champ.dim_c_vector,
                                                ctx->jastrow_champ.c_vector_full,
                                                ctx->jastrow_champ.lkpm_combined_index,
                                                ctx->jastrow_champ.een_rescaled_n,
                                                ctx->jastrow_champ.een_rescaled_e,
                                                ctx->single_point.een_rescaled_single_n,
                                                ctx->single_point.een_rescaled_single_e,
                                                ctx->single_point.delta_een_pderiv);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->single_point.delta_een_pderiv_date = ctx->single_point.date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_single_een_pderiv (const qmckl_context context,
                          const int64_t num,
                          const int64_t walk_num,
                          const int64_t elec_num,
                          const int64_t nucl_num,
                          const int64_t type_nucl_num,
                          const int64_t* type_nucl_vector,
                          const int64_t cord_num,
                          const int64_t dim_c_vector,
                          const double* c_vector_full,
                          const int64_t* lkpm_combined_index,
                          const double* een_rescaled_n,
                          const double* een_rescaled_e,
                          const double* een_rescaled_single_n,
                          const double* een_rescaled_single_e,
                          double* const delta_een_pderiv )
{
#ifdef HAVE_HPC
  return qmckl_compute_jastrow_champ_factor_single_een_pderiv_doc
#else
  return qmckl_compute_jastrow_champ_factor_single_een_pderiv_doc
#endif
    (context, num, walk_num, elec_num, nucl_num, type_nucl_num, type_nucl_vector, cord_num, dim_c_vector,
     c_vector_full, lkpm_combined_index, een_rescaled_n, een_rescaled_e, een_rescaled_single_n, een_rescaled_single_e, delta_een_pderiv );

}

qmckl_exit_code
qmckl_get_ee_rescaled_single(qmckl_context context,
                                 double* const distance_rescaled,
                                 const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_ee_rescaled_single(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t sze = ctx->electron.num * ctx->electron.walker.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "todo",
                           "Array too small. Expected ctx->electron.num * ctx->electron.walker.num ");
  }
  memcpy(distance_rescaled, ctx->single_point.ee_rescaled_single, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_ee_rescaled_single(qmckl_context context)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);


  /* Check if ee distance is provided */
  qmckl_exit_code rc = qmckl_provide_single_ee_distance(context);
  if(rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_ee_distance(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->single_point.date > ctx->single_point.ee_rescaled_single_date) {

    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = ctx->electron.num * ctx->electron.walker.num  * sizeof(double);

    if (mem_info.size > ctx->single_point.ee_rescaled_single_maxsize) {
      if (ctx->single_point.ee_rescaled_single!= NULL) {
        rc = qmckl_free(context, ctx->single_point.ee_rescaled_single);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_ee_rescaled_single",
                                 "Unable to free ctx->single_point.ee_rescaled_single");
        }
        ctx->single_point.ee_rescaled_single = NULL;
      }
    }

    /* Allocate array */
    if (ctx->single_point.ee_rescaled_single == NULL) {

      double* ee_rescaled_single = (double*) qmckl_malloc(context, mem_info);
      ctx->single_point.ee_rescaled_single_maxsize = mem_info.size;

      if (ee_rescaled_single == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_ee_rescaled_single",
                               NULL);
      }

      ctx->single_point.ee_rescaled_single = ee_rescaled_single;
    }

    rc = qmckl_compute_ee_rescaled_single(context,
                                      ctx->electron.num,
                                      ctx->jastrow_champ.rescale_factor_ee,
                                      ctx->electron.walker.num,
                                      ctx->single_point.single_ee_distance,
                                      ctx->single_point.ee_rescaled_single);

    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->single_point.ee_rescaled_single_date = ctx->single_point.date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_compute_ee_rescaled_single (
          const qmckl_context context,
          const int64_t elec_num,
          const double rescale_factor_ee,
          const int64_t walk_num,
          const double* single_ee_distance,
          double* const ee_rescaled_single )
{
#ifdef HAVE_HPC
  return qmckl_compute_ee_rescaled_single_doc
#else
  return qmckl_compute_ee_rescaled_single_doc
#endif
 (context, elec_num, rescale_factor_ee, walk_num, single_ee_distance, ee_rescaled_single);
}

qmckl_exit_code
qmckl_get_jastrow_champ_single_ee(qmckl_context context,
                            double* const delta_ee,
                            const int64_t size_max)
{
  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_jastrow_champ_single_ee",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  rc = qmckl_provide_jastrow_champ_single_ee(context);
  if (rc != QMCKL_SUCCESS) return rc;

  int64_t sze=ctx->electron.walker.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_single_ee",
                           "Array too small. Expected walker.num");
  }
  memcpy(delta_ee, ctx->single_point.delta_ee, sze*sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_jastrow_champ_single_ee(qmckl_context context)
{

  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_provide_jastrow_champ_single_ee",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->jastrow_champ.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_provide_jastrow_champ_single_ee",
                           NULL);
  }

  rc = qmckl_provide_ee_distance_rescaled(context);
  if(rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_ee_rescaled_single(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->single_point.date > ctx->single_point.delta_ee_date) {

    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = ctx->electron.walker.num * sizeof(double);

    if (mem_info.size > ctx->single_point.delta_ee_maxsize) {
      if (ctx->single_point.delta_ee != NULL) {
        rc = qmckl_free(context, ctx->single_point.delta_ee);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_jastrow_champ_single_ee",
                                 "Unable to free ctx->single_point.delta_ee");
        }
        ctx->single_point.delta_ee = NULL;
      }
    }

    /* Allocate array */
    if (ctx->single_point.delta_ee == NULL) {

      double* delta_ee = (double*) qmckl_malloc(context, mem_info);
      ctx->single_point.delta_ee_maxsize = mem_info.size;

      if (delta_ee == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_jastrow_champ_single_ee",
                               NULL);
      }
      ctx->single_point.delta_ee = delta_ee;
    }

    rc = qmckl_compute_jastrow_champ_single_ee(context,
                                 ctx->single_point.num,
                                 ctx->electron.walker.num,
                                 ctx->electron.num,
                                 ctx->electron.up_num,
                                 ctx->jastrow_champ.bord_num,
                                 ctx->jastrow_champ.b_vector,
                                 ctx->jastrow_champ.ee_distance_rescaled,
                                 ctx->single_point.ee_rescaled_single,
                                 ctx->jastrow_champ.spin_independent,
                                 ctx->single_point.delta_ee);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->single_point.delta_ee_date = ctx->single_point.date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_compute_jastrow_champ_single_ee (const qmckl_context context,
                                        const int64_t num,
                                       const int64_t walk_num,
                                       const int64_t elec_num,
                                       const int64_t up_num,
                                       const int64_t bord_num,
                                       const double* b_vector,
                                       const double* ee_distance_rescaled,
                                       const double* ee_rescaled_single,
                                       const int32_t spin_independent,
                                       double* const delta_ee )
{

#ifdef HAVE_HPC
  return qmckl_compute_jastrow_champ_single_ee_doc
#else
  return qmckl_compute_jastrow_champ_single_ee_doc
#endif
    (context, num, walk_num, elec_num, up_num, bord_num, b_vector,
     ee_distance_rescaled, ee_rescaled_single, spin_independent, delta_ee);
}

qmckl_exit_code qmckl_get_ee_rescaled_single_gl(qmckl_context context,
                                                double* const distance_rescaled_gl,
                                                const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_ee_rescaled_single_gl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (distance_rescaled_gl == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_ee_rescaled_single_gl",
                           "Array is NULL");
  }

  int64_t sze = 4 * ctx->electron.num  * ctx->electron.walker.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_ee_rescaled_single_gl",
                           "Array too small. Expected 4 * ctx->electron.num * ctx->electron.walker.num");
  }
  memcpy(distance_rescaled_gl, ctx->single_point.ee_rescaled_single_gl, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_ee_rescaled_single_gl(qmckl_context context)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  qmckl_exit_code rc = qmckl_provide_single_ee_distance(context);
  if(rc != QMCKL_SUCCESS) return rc;


  /* Compute if necessary */
  if (ctx->single_point.date > ctx->single_point.ee_rescaled_single_gl_date) {

    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = 4  * ctx->electron.num * ctx->electron.walker.num * sizeof(double);

    if (mem_info.size > ctx->single_point.ee_rescaled_single_gl_maxsize) {
      if (ctx->single_point.ee_rescaled_single_gl != NULL) {
        rc = qmckl_free(context, ctx->single_point.ee_rescaled_single_gl);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_ee_rescaled_single_gl",
                                 "Unable to free ctx->single_point.ee_rescaled_single_gl");
        }
        ctx->single_point.ee_rescaled_single_gl = NULL;
      }
    }

    /* Allocate array */
    if (ctx->single_point.ee_rescaled_single_gl == NULL) {

      double* ee_rescaled_single_gl = (double*) qmckl_malloc(context, mem_info);
      ctx->single_point.ee_rescaled_single_gl_maxsize = mem_info.size;

      if (ee_rescaled_single_gl == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_ee_rescaled_single_gl",
                               NULL);
      }
      ctx->single_point.ee_rescaled_single_gl = ee_rescaled_single_gl;
    }

    qmckl_exit_code rc =
      qmckl_compute_ee_rescaled_single_gl(context,
                                ctx->single_point.num,
                                ctx->electron.num,
                                ctx->jastrow_champ.rescale_factor_ee,
                                ctx->electron.walker.num,
                                ctx->single_point.single_ee_distance,
                                ctx->electron.walker.point.coord.data,
                                ctx->single_point.coord.data,
                                ctx->single_point.ee_rescaled_single_gl);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->single_point.ee_rescaled_single_gl_date = ctx->single_point.date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_compute_ee_rescaled_single_gl (
          const qmckl_context context,
          const int64_t num,
          const int64_t elec_num,
          const double rescale_factor_ee,
          const int64_t walk_num,
          const double* single_ee_distance,
          const double* elec_coord,
          const double* coord,
          double* const ee_rescaled_single_gl )
{
#ifdef HAVE_HPC
  return qmckl_compute_ee_rescaled_single_gl_doc
#else
  return qmckl_compute_ee_rescaled_single_gl_doc
#endif
   (context, num, elec_num, rescale_factor_ee, walk_num, single_ee_distance, elec_coord, coord,
   ee_rescaled_single_gl);
}

qmckl_exit_code
qmckl_get_jastrow_champ_single_ee_gl(qmckl_context context,
                                    double* const delta_ee_gl,
                                    const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_jastrow_champ_single_ee_gl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t sze = ctx->electron.walker.num * 4 * ctx->electron.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_single_ee_gl",
                           "Array too small. Expected 4*walk_num*elec_num");
  }

  memcpy(delta_ee_gl, ctx->single_point.delta_ee_gl, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_jastrow_champ_single_ee_gl(qmckl_context context)
{

  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_provide_jastrow_champ_single_ee_gl",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->jastrow_champ.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_provide_jastrow_champ_single_ee_gl",
                           NULL);
  }

  /* Check if ee rescaled distance is provided */
  rc = qmckl_provide_ee_distance_rescaled(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Check if ee rescaled distance deriv e is provided */
  rc = qmckl_provide_ee_distance_rescaled_gl(context);
  if(rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_ee_rescaled_single(context);
  if(rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_ee_rescaled_single_gl(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->single_point.date > ctx->single_point.delta_ee_gl_date) {

    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = ctx->electron.walker.num * 4 * ctx->electron.num * sizeof(double);

    if (mem_info.size > ctx->single_point.delta_ee_gl_maxsize) {
      if (ctx->single_point.delta_ee_gl != NULL) {
        rc = qmckl_free(context, ctx->single_point.delta_ee_gl);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_jastrow_champ_single_ee_gl",
                                 "Unable to free ctx->single_point.delta_ee_gl");
        }
        ctx->single_point.delta_ee_gl = NULL;
      }
    }

    /* Allocate array */
    if (ctx->single_point.delta_ee_gl == NULL) {

      double* delta_ee_gl = (double*) qmckl_malloc(context, mem_info);
      ctx->single_point.delta_ee_gl_maxsize = mem_info.size;

      if (delta_ee_gl == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_jastrow_champ_single_ee_gl",
                               NULL);
      }
      ctx->single_point.delta_ee_gl = delta_ee_gl;
    }

    rc = qmckl_compute_jastrow_champ_single_ee_gl(context,
                                                  ctx->single_point.num,
                                                  ctx->electron.walker.num,
                                                  ctx->electron.num,
                                                  ctx->electron.up_num,
                                                  ctx->jastrow_champ.bord_num,
                                                  ctx->jastrow_champ.b_vector,
                                                  ctx->jastrow_champ.ee_distance_rescaled,
                                                  ctx->jastrow_champ.ee_distance_rescaled_gl,
                                                  ctx->single_point.ee_rescaled_single,
                                                  ctx->single_point.ee_rescaled_single_gl,
                                                  ctx->jastrow_champ.spin_independent,
                                                  ctx->single_point.delta_ee_gl);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->single_point.delta_ee_gl_date = ctx->single_point.date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_compute_jastrow_champ_single_ee_gl (const qmckl_context context,
                                          const int64_t num,
                                          const int64_t walk_num,
                                          const int64_t elec_num,
                                          const int64_t up_num,
                                          const int64_t bord_num,
                                          const double* b_vector,
                                          const double* ee_distance_rescaled,
                                          const double* ee_distance_rescaled_gl,
                                          const double* ee_rescaled_single,
                                          const double* ee_rescaled_single_gl,
                                          const int32_t spin_independent,
                                          double* const delta_ee_gl )
{
#ifdef HAVE_HPC
  return qmckl_compute_jastrow_champ_single_ee_gl_doc
#else
    return qmckl_compute_jastrow_champ_single_ee_gl_doc
#endif
    (context, num, walk_num, elec_num, up_num, bord_num, b_vector,
     ee_distance_rescaled, ee_distance_rescaled_gl, ee_rescaled_single, ee_rescaled_single_gl, spin_independent, delta_ee_gl );
}

qmckl_exit_code
qmckl_get_jastrow_champ_single_ee_pderiv(qmckl_context context,
                            double* const delta_ee_pderiv,
                            const int64_t size_max)
{
  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_jastrow_champ_single_ee_pderiv",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  rc = qmckl_provide_jastrow_champ_single_ee_pderiv(context);
  if (rc != QMCKL_SUCCESS) return rc;

  int64_t sze = ctx->jastrow_champ.bord_num + 1;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_single_ee_pderiv",
                           "Array too small. Expected bord_num + 1");
  }
  memcpy(delta_ee_pderiv, ctx->single_point.delta_ee_pderiv, sze*sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_jastrow_champ_single_ee_pderiv(qmckl_context context)
{

  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_provide_jastrow_champ_single_ee_pderiv",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->jastrow_champ.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_provide_jastrow_champ_single_ee_pderiv",
                           NULL);
  }

  rc = qmckl_provide_ee_distance_rescaled(context);
  if(rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_ee_rescaled_single(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->single_point.date > ctx->single_point.delta_ee_pderiv_date) {

    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = (ctx->jastrow_champ.bord_num+1) * sizeof(double);

    /* Allocate array */
    if (ctx->single_point.delta_ee_pderiv == NULL) {

      double* delta_ee_pderiv = (double*) qmckl_malloc(context, mem_info);
      ctx->single_point.delta_ee_pderiv_maxsize = mem_info.size;

      if (delta_ee_pderiv == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_jastrow_champ_single_ee_pderiv",
                               NULL);
      }
      ctx->single_point.delta_ee_pderiv = delta_ee_pderiv;
    }

    rc = qmckl_compute_jastrow_champ_single_ee_pderiv(context,
                                 ctx->single_point.num,
                                 ctx->electron.walker.num,
                                 ctx->electron.num,
                                 ctx->electron.up_num,
                                 ctx->jastrow_champ.bord_num,
                                 ctx->jastrow_champ.b_vector,
                                 ctx->jastrow_champ.ee_distance_rescaled,
                                 ctx->single_point.ee_rescaled_single,
                                 ctx->jastrow_champ.spin_independent,
                                 ctx->single_point.delta_ee_pderiv);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->single_point.delta_ee_pderiv_date = ctx->single_point.date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_compute_jastrow_champ_single_ee_pderiv (const qmckl_context context,
                                        const int64_t num,
                                       const int64_t walk_num,
                                       const int64_t elec_num,
                                       const int64_t up_num,
                                       const int64_t bord_num,
                                       const double* b_vector,
                                       const double* ee_distance_rescaled,
                                       const double* ee_rescaled_single,
                                       const int32_t spin_independent,
                                       double* const delta_ee_pderiv )
{

#ifdef HAVE_HPC
  return qmckl_compute_jastrow_champ_single_ee_pderiv_doc
#else
  return qmckl_compute_jastrow_champ_single_ee_pderiv_doc
#endif
    (context, num, walk_num, elec_num, up_num, bord_num, b_vector,
     ee_distance_rescaled, ee_rescaled_single, spin_independent, delta_ee_pderiv);
}

qmckl_exit_code
qmckl_get_en_rescaled_single(qmckl_context context,
                                 double* const distance_rescaled,
                                 const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_en_rescaled_single(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t sze = ctx->nucleus.num * ctx->electron.walker.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "todo",
                           "Array too small. Expected ctx->nucleus.num * ctx->electron.walker.num ");
  }
  memcpy(distance_rescaled, ctx->single_point.en_rescaled_single, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_en_rescaled_single(qmckl_context context)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  /* Check if ee distance is provided */
  qmckl_exit_code rc = qmckl_provide_single_en_distance(context);
  if(rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_en_distance(context);
  if(rc != QMCKL_SUCCESS) return rc;


  /* Compute if necessary */
  if (ctx->single_point.date > ctx->single_point.en_rescaled_single_date) {

    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = ctx->nucleus.num * ctx->electron.walker.num  * sizeof(double);

    if (mem_info.size > ctx->single_point.en_rescaled_single_maxsize) {
      if (ctx->single_point.en_rescaled_single!= NULL) {
        rc = qmckl_free(context, ctx->single_point.en_rescaled_single);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_en_rescaled_single",
                                 "Unable to free ctx->single_point.en_rescaled_single");
        }
        ctx->single_point.en_rescaled_single = NULL;
      }
    }

    /* Allocate array */
    if (ctx->single_point.en_rescaled_single == NULL) {

      double* en_rescaled_single = (double*) qmckl_malloc(context, mem_info);
      ctx->single_point.en_rescaled_single_maxsize = mem_info.size;

      if (en_rescaled_single == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_en_rescaled_single",
                               NULL);
      }

      ctx->single_point.en_rescaled_single = en_rescaled_single;
    }

    rc = qmckl_compute_en_rescaled_single(context,
                                ctx->nucleus.num,
                                ctx->jastrow_champ.type_nucl_num,
                                ctx->jastrow_champ.type_nucl_vector,
                                ctx->jastrow_champ.rescale_factor_en,
                                ctx->electron.walker.num,
                                ctx->single_point.single_en_distance,
                                ctx->single_point.en_rescaled_single);

    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->single_point.en_rescaled_single_date = ctx->single_point.date;

  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_compute_en_rescaled_single(
          const qmckl_context context,
          const int64_t nucl_num,
          const int64_t type_nucl_num,
          const int64_t* type_nucl_vector,
          const double* rescale_factor_en,
          const int64_t walk_num,
          const double* single_en_distance,
          double* const en_rescaled_single )
{
#ifdef HAVE_HPC
  return qmckl_compute_en_rescaled_single_doc
#else
  return qmckl_compute_en_rescaled_single_doc
#endif
  (context, nucl_num, type_nucl_num, type_nucl_vector,
  rescale_factor_en, walk_num, single_en_distance, en_rescaled_single );
}

qmckl_exit_code
qmckl_get_jastrow_champ_single_en(qmckl_context context,
                            double* const delta_en,
                            const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_jastrow_champ_single_en",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  qmckl_exit_code rc;

  rc = qmckl_provide_jastrow_champ_single_en(context);
  if (rc != QMCKL_SUCCESS) return rc;

  int64_t sze=ctx->electron.walker.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_single_en",
                           "Array too small. Expected walker.num");
  }
  memcpy(delta_en, ctx->single_point.delta_en, sze*sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_jastrow_champ_single_en(qmckl_context context)
{

  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_provide_jastrow_champ_single_en",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->jastrow_champ.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_provide_jastrow_champ_single_en",
                           NULL);
  }

  /* Check if en rescaled distance is provided */
  rc = qmckl_provide_en_distance_rescaled(context);
  if(rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_en_rescaled_single(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->single_point.date > ctx->single_point.delta_en_date) {

    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = ctx->electron.walker.num * sizeof(double);

    if (mem_info.size > ctx->single_point.delta_en_maxsize) {
      if (ctx->single_point.delta_en != NULL) {
        rc = qmckl_free(context, ctx->single_point.delta_en);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_jastrow_champ_single_en",
                                 "Unable to free ctx->single_point.delta_en");
        }
        ctx->single_point.delta_en = NULL;
      }
    }

    /* Allocate array */
    if (ctx->single_point.delta_en == NULL) {

      double* delta_en = (double*) qmckl_malloc(context, mem_info);
      ctx->single_point.delta_en_maxsize = mem_info.size;

      if (delta_en == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_jastrow_champ_single_en",
                               NULL);
      }
      ctx->single_point.delta_en = delta_en;
    }

    rc = qmckl_compute_jastrow_champ_single_en(context,
                                 ctx->single_point.num,
                                 ctx->electron.walker.num,
                                 ctx->electron.num,
                                 ctx->nucleus.num,
                                 ctx->jastrow_champ.type_nucl_num,
                                 ctx->jastrow_champ.type_nucl_vector,
                                 ctx->jastrow_champ.aord_num,
                                 ctx->jastrow_champ.a_vector,
                                 ctx->jastrow_champ.en_distance_rescaled,
                                 ctx->single_point.en_rescaled_single,
                                 ctx->single_point.delta_en);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->single_point.delta_en_date = ctx->single_point.date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_compute_jastrow_champ_single_en (
         const qmckl_context context,
         const int64_t num,
         const int64_t walk_num,
         const int64_t elec_num,
         const int64_t nucl_num,
         const int64_t type_nucl_num,
         const int64_t* type_nucl_vector,
         const int64_t aord_num,
         const double* a_vector,
         const double* en_distance_rescaled,
         const double* en_rescaled_single,
         double* const delta_en )
{
#ifdef HAVE_HPC
  return qmckl_compute_jastrow_champ_single_en_doc
#else
  return qmckl_compute_jastrow_champ_single_en_doc
#endif
    (context, num, walk_num, elec_num, nucl_num, type_nucl_num,
     type_nucl_vector, aord_num, a_vector, en_distance_rescaled,
     en_rescaled_single, delta_en );
}

qmckl_exit_code qmckl_get_en_rescaled_single_gl(qmckl_context context,
                                                double* distance_rescaled_gl,
                                                const size_t size_max)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_en_rescaled_single_gl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  size_t sze = 4 * ctx->nucleus.num * ctx->electron.walker.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_en_rescaled_single_gl",
                           "Array too small. Expected 4 * ctx->nucleus.num * ctx->electron.walker.num");
  }
  memcpy(distance_rescaled_gl, ctx->single_point.en_rescaled_single_gl, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_en_rescaled_single_gl(qmckl_context context)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!(ctx->nucleus.provided)) {
    return QMCKL_NOT_PROVIDED;
  }

  qmckl_exit_code rc = qmckl_provide_single_en_distance(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->single_point.date > ctx->single_point.en_rescaled_single_gl_date) {

    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = 4 * ctx->nucleus.num * ctx->electron.walker.num * sizeof(double);

    if (mem_info.size > ctx->single_point.en_rescaled_single_gl_maxsize) {
      if (ctx->single_point.en_rescaled_single_gl != NULL) {
        rc = qmckl_free(context, ctx->single_point.en_rescaled_single_gl);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_en_rescaled_single_gl",
                                 "Unable to free ctx->single_point.en_rescaled_single_gl");
        }
        ctx->single_point.en_rescaled_single_gl = NULL;
      }
    }

    /* Allocate array */
    if (ctx->single_point.en_rescaled_single_gl == NULL) {

      double* en_rescaled_single_gl = (double*) qmckl_malloc(context, mem_info);
      ctx->single_point.en_rescaled_single_gl_maxsize = mem_info.size;

      if (en_rescaled_single_gl == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_en_rescaled_single_gl",
                               NULL);
      }
      ctx->single_point.en_rescaled_single_gl = en_rescaled_single_gl;
    }

    qmckl_exit_code rc =
      qmckl_compute_en_rescaled_single_gl(context,
                                ctx->nucleus.num,
                                ctx->jastrow_champ.type_nucl_num,
                                ctx->jastrow_champ.type_nucl_vector,
                                ctx->jastrow_champ.rescale_factor_en,
                                ctx->electron.walker.num,
                                ctx->single_point.single_en_distance,
                                ctx->single_point.coord.data,
                                ctx->nucleus.coord.data,
                                ctx->single_point.en_rescaled_single_gl);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->single_point.en_rescaled_single_gl_date = ctx->single_point.date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_compute_en_rescaled_single_gl (
          const qmckl_context context,
          const int64_t nucl_num,
          const int64_t type_nucl_num,
          int64_t* const type_nucl_vector,
          const double*  rescale_factor_en,
          const int64_t walk_num,
          const double* single_en_distance,
          const double* coord,
          const double* nucl_coord,
          double* const en_rescaled_single_gl )
{
#ifdef HAVE_HPC
  return qmckl_compute_en_rescaled_single_gl_doc
#else
  return qmckl_compute_en_rescaled_single_gl_doc
#endif
    (context, nucl_num, type_nucl_num, type_nucl_vector, rescale_factor_en,
     walk_num, single_en_distance, coord, nucl_coord, en_rescaled_single_gl );
}

qmckl_exit_code
qmckl_get_jastrow_champ_single_en_gl(qmckl_context context,
                                    double* const delta_en_gl,
                                    const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_jastrow_champ_single_en_gl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t sze = ctx->electron.walker.num * 4 * ctx->electron.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_single_en_gl",
                           "Array too small. Expected 4*walker.num*elec_num");
  }
  memcpy(delta_en_gl, ctx->single_point.delta_en_gl, sze*sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_jastrow_champ_single_en_gl(qmckl_context context)
{

  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_provide_jastrow_champ_single_en_gl",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->jastrow_champ.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_provide_jastrow_champ_single_en_gl",
                           NULL);
  }

  /* Check if en rescaled distance is provided */
  rc = qmckl_provide_en_distance_rescaled(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Check if en rescaled distance derivatives is provided */
  rc = qmckl_provide_en_distance_rescaled_gl(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Check if en rescaled distance is provided */
  rc = qmckl_provide_en_rescaled_single(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Check if en rescaled distance derivatives is provided */
  rc = qmckl_provide_en_rescaled_single_gl(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->single_point.date > ctx->single_point.delta_en_gl_date) {

    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = ctx->electron.walker.num * 4 * ctx->electron.num * sizeof(double);

    if (mem_info.size > ctx->single_point.delta_en_gl_maxsize) {
      if (ctx->single_point.delta_en_gl != NULL) {
        rc = qmckl_free(context, ctx->single_point.delta_en_gl);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_jastrow_champ_single_en_gl",
                                 "Unable to free ctx->single_point.delta_en_gl");
        }
        ctx->single_point.delta_en_gl = NULL;
      }
    }

    /* Allocate array */
    if (ctx->single_point.delta_en_gl == NULL) {

      double* delta_en_gl = (double*) qmckl_malloc(context, mem_info);
      ctx->single_point.delta_en_gl_maxsize = mem_info.size;

      if (delta_en_gl == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_jastrow_champ_single_en_gl",
                               NULL);
      }
      ctx->single_point.delta_en_gl = delta_en_gl;
    }

    rc = qmckl_compute_jastrow_champ_single_en_gl(context,
                                         ctx->single_point.num,
                                         ctx->electron.walker.num,
                                         ctx->electron.num,
                                         ctx->nucleus.num,
                                         ctx->jastrow_champ.type_nucl_num,
                                         ctx->jastrow_champ.type_nucl_vector,
                                         ctx->jastrow_champ.aord_num,
                                         ctx->jastrow_champ.a_vector,
                                         ctx->jastrow_champ.en_distance_rescaled,
                                         ctx->jastrow_champ.en_distance_rescaled_gl,
                                         ctx->single_point.en_rescaled_single,
                                         ctx->single_point.en_rescaled_single_gl,
                                         ctx->single_point.delta_en_gl);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->single_point.delta_en_gl_date = ctx->single_point.date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_compute_jastrow_champ_single_en_gl (const qmckl_context context,
                                          const int64_t num,
                                          const int64_t walk_num,
                                          const int64_t elec_num,
                                          const int64_t nucl_num,
                                          const int64_t type_nucl_num,
                                          const int64_t* type_nucl_vector,
                                          const int64_t aord_num,
                                          const double* a_vector,
                                          const double* en_distance_rescaled,
                                          const double* en_distance_rescaled_gl,
                                          const double* en_rescaled_single,
                                          const double* en_rescaled_single_gl,
                                          double* const delta_en_gl )
{
#ifdef HAVE_HPC
  return qmckl_compute_jastrow_champ_single_en_gl_doc
#else
  return qmckl_compute_jastrow_champ_single_en_gl_doc
#endif
    (context, num, walk_num, elec_num, nucl_num, type_nucl_num, type_nucl_vector, aord_num,
     a_vector, en_distance_rescaled, en_distance_rescaled_gl, en_rescaled_single, en_rescaled_single_gl, delta_en_gl );
}

qmckl_exit_code
qmckl_get_jastrow_champ_single_en_pderiv(qmckl_context context,
                            double* const delta_en_pderiv,
                            const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_jastrow_champ_single_en_pderiv",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  qmckl_exit_code rc;

  rc = qmckl_provide_jastrow_champ_single_en_pderiv(context);
  if (rc != QMCKL_SUCCESS) return rc;

  int64_t sze=(ctx->jastrow_champ.aord_num+1) * ctx->jastrow_champ.type_nucl_num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_jastrow_champ_single_en_pderiv",
                           "Array too small. Expected aord_num + 1");
  }
  memcpy(delta_en_pderiv, ctx->single_point.delta_en_pderiv, sze*sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_jastrow_champ_single_en_pderiv(qmckl_context context)
{

  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_provide_jastrow_champ_single_en_pderiv",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->jastrow_champ.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_provide_jastrow_champ_single_en_pderiv",
                           NULL);
  }

  /* Check if en rescaled distance is provided */
  rc = qmckl_provide_en_distance_rescaled(context);
  if(rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_en_rescaled_single(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->single_point.date > ctx->single_point.delta_en_pderiv_date) {

    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = (ctx->jastrow_champ.aord_num + 1) * ctx->jastrow_champ.type_nucl_num * sizeof(double);

    /* Allocate array */
    if (ctx->single_point.delta_en_pderiv == NULL) {

      double* delta_en_pderiv = (double*) qmckl_malloc(context, mem_info);
      // ctx->single_point.delta_en_pderiv_maxsize = mem_info.size;

      if (delta_en_pderiv == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_jastrow_champ_single_en_pderiv",
                               NULL);
      }
      ctx->single_point.delta_en_pderiv = delta_en_pderiv;
    }

    rc = qmckl_compute_jastrow_champ_single_en_pderiv(context,
                                 ctx->single_point.num,
                                 ctx->electron.walker.num,
                                 ctx->electron.num,
                                 ctx->nucleus.num,
                                 ctx->jastrow_champ.type_nucl_num,
                                 ctx->jastrow_champ.type_nucl_vector,
                                 ctx->jastrow_champ.aord_num,
                                 ctx->jastrow_champ.a_vector,
                                 ctx->jastrow_champ.en_distance_rescaled,
                                 ctx->single_point.en_rescaled_single,
                                 ctx->single_point.delta_en_pderiv);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->single_point.delta_en_pderiv_date = ctx->single_point.date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_compute_jastrow_champ_single_en_pderiv (
         const qmckl_context context,
         const int64_t num,
         const int64_t walk_num,
         const int64_t elec_num,
         const int64_t nucl_num,
         const int64_t type_nucl_num,
         const int64_t* type_nucl_vector,
         const int64_t aord_num,
         const double* a_vector,
         const double* en_distance_rescaled,
         const double* en_rescaled_single,
         double* const delta_en_pderiv )
{
#ifdef HAVE_HPC
  return qmckl_compute_jastrow_champ_single_en_pderiv_doc
#else
  return qmckl_compute_jastrow_champ_single_en_pderiv_doc
#endif
    (context, num, walk_num, elec_num, nucl_num, type_nucl_num,
     type_nucl_vector, aord_num, a_vector, en_distance_rescaled,
     en_rescaled_single, delta_en_pderiv );
}

qmckl_exit_code
qmckl_get_jastrow_champ_single_accept(qmckl_context context)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  qmckl_exit_code rc;

  rc = qmckl_provide_jastrow_champ_single_en_gl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_jastrow_champ_single_ee_gl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  if (ctx->jastrow_champ.cord_num > 0) {
    rc = qmckl_provide_jastrow_champ_single_een_gl(context);
    if (rc != QMCKL_SUCCESS) return rc;
  }

  double metric[4] = {-1.0, -1.0, -1.0, 1.0};

  int shift1, shift2, shift3, shift4, shift5, shift6, shift7;

  int shift8, shift9, shift10, shift11, shift12, shift13, shift14;

  shift2 = ctx->electron.num*ctx->electron.num;
  shift5 = ctx->electron.num*4*ctx->electron.num;
  shift6 = ctx->electron.num*4;
  shift11 = ctx->nucleus.num*4*ctx->electron.num;
  shift13 = ctx->nucleus.num*4;
  shift14 = ctx->nucleus.num*ctx->electron.num;


  if (ctx->jastrow_champ.cord_num > 0) {

    shift1 = (ctx->jastrow_champ.cord_num+1)*ctx->electron.num*ctx->electron.num;

    shift3 = (ctx->jastrow_champ.cord_num+1)*ctx->electron.num;

    shift4 = (ctx->jastrow_champ.cord_num+1)*ctx->electron.num*4*ctx->electron.num;

    shift7 = (ctx->jastrow_champ.cord_num+1)*ctx->electron.num*4;

    shift8 = (ctx->jastrow_champ.cord_num+1)*ctx->nucleus.num*ctx->electron.num;
    shift9 = (ctx->jastrow_champ.cord_num+1)*ctx->nucleus.num;

    shift10 = (ctx->jastrow_champ.cord_num+1)*ctx->nucleus.num*4*ctx->electron.num;

    shift12 = (ctx->jastrow_champ.cord_num+1)*ctx->nucleus.num*4;

    for (int nw = 0; nw < ctx->electron.walker.num; nw++) {
      ctx->jastrow_champ.factor_een[nw] = ctx->jastrow_champ.factor_een[nw] + ctx->single_point.delta_een[nw];
      for (int l = 0; l <= ctx->jastrow_champ.cord_num; l++){
        for (int i = 0; i < ctx->electron.num; i++) {
          ctx->jastrow_champ.een_rescaled_e[nw*shift1
            + l*shift2
            + i*ctx->electron.num
            + ctx->single_point.num] =
          ctx->single_point.een_rescaled_single_e[nw*shift3
            + l*ctx->electron.num
            + i];

          ctx->jastrow_champ.een_rescaled_e[nw*shift1
            + l*shift2
            + ctx->single_point.num*ctx->electron.num
            + i] =
          ctx->single_point.een_rescaled_single_e[nw*shift3
            + l*ctx->electron.num
            + i];
          for (int k = 0; k < 4; k++){
            ctx->jastrow_champ.een_rescaled_e_gl[nw*shift4
              + l*shift5
              + i*shift6
              + k*ctx->electron.num
              + ctx->single_point.num] =
            ctx->single_point.een_rescaled_single_e_gl[nw*shift7
              + l*shift6
              + k*ctx->electron.num
              + i];
            ctx->jastrow_champ.een_rescaled_e_gl[nw*shift4
              + l*shift5
              + ctx->single_point.num*shift6
              + k*ctx->electron.num
              + i] =
            metric[k] * ctx->single_point.een_rescaled_single_e_gl[nw*shift7
              + l*shift6
              + k*ctx->electron.num
              + i];
          }
        }
        for (int a = 0; a < ctx->nucleus.num; a++){
          ctx->jastrow_champ.een_rescaled_n[nw*shift8
            + l*shift14
            + a*ctx->electron.num
            + ctx->single_point.num] =
          ctx->single_point.een_rescaled_single_n[nw*shift9
            + l*ctx->nucleus.num
            + a];
          for (int k = 0; k < 4; k++){
            ctx->jastrow_champ.een_rescaled_n_gl[nw*shift10
              + l*shift11
              + a*shift6
              + k*ctx->electron.num
              + ctx->single_point.num] =
            ctx->single_point.een_rescaled_single_n_gl[nw*shift12
              + l*shift13
              + a*4
              + k];
          }
        }
      }
    }
    for (int i = 0; i < ctx->electron.walker.num * 4 * ctx->electron.num; i++) {
      ctx->jastrow_champ.factor_een_gl[i] = ctx->jastrow_champ.factor_een_gl[i] + ctx->single_point.delta_een_gl[i];
    }
  }



  for (int nw = 0; nw < ctx->electron.walker.num; nw++) {
    ctx->jastrow_champ.factor_en[nw] = ctx->jastrow_champ.factor_en[nw] + ctx->single_point.delta_en[nw];
    ctx->jastrow_champ.factor_ee[nw] = ctx->jastrow_champ.factor_ee[nw] + ctx->single_point.delta_ee[nw];
    for (int a = 0; a < ctx->nucleus.num; a++) {
      ctx->electron.en_distance[nw*shift14
        + ctx->single_point.num*ctx->nucleus.num
        + a] =
      ctx->single_point.single_en_distance[nw*ctx->nucleus.num
        + a];

      ctx->jastrow_champ.en_distance_rescaled[nw*shift14
        + a*ctx->electron.num
        + ctx->single_point.num] =
      ctx->single_point.en_rescaled_single[nw*ctx->nucleus.num
        + a];
      for (int k = 0; k < 4; k++){
        ctx->jastrow_champ.en_distance_rescaled_gl[nw*shift11
          + a*shift6
          + ctx->single_point.num*4
          + k] =
        ctx->single_point.en_rescaled_single_gl[nw*shift13
          + a*4
          + k];
      }
    }
    for (int i = 0; i < ctx->electron.num; i++) {
      ctx->jastrow_champ.ee_distance_rescaled[nw*shift2
        + i*ctx->electron.num
        + ctx->single_point.num] =
      ctx->single_point.ee_rescaled_single[nw*ctx->electron.num
        + i];

      ctx->jastrow_champ.ee_distance_rescaled[nw*shift2
        + ctx->single_point.num*ctx->electron.num
        + i] =
      ctx->single_point.ee_rescaled_single[nw*ctx->electron.num
        + i];

      ctx->electron.ee_distance[nw*shift2
        + i*ctx->electron.num
        + ctx->single_point.num] =
      ctx->single_point.single_ee_distance[nw*ctx->electron.num
        + i];

      ctx->electron.ee_distance[nw*shift2
        + ctx->single_point.num*ctx->electron.num
        + i] =
      ctx->single_point.single_ee_distance[nw*ctx->electron.num
        + i];

      for (int k = 0; k < 4; k++){
        ctx->jastrow_champ.ee_distance_rescaled_gl[nw*shift5
          + i*shift6
          + ctx->single_point.num*4
          + k] =
        metric[k] * ctx->single_point.ee_rescaled_single_gl[nw*shift6
          + i*4
          + k];
        ctx->jastrow_champ.ee_distance_rescaled_gl[nw*shift5
          + ctx->single_point.num*shift6
          + i*4
          + k] =
        ctx->single_point.ee_rescaled_single_gl[nw*shift6
          + i*4
          + k];
      }
    }
    for (int k = 0; k < 4; k++){
      for (int i = 0; i < ctx->electron.num; i++) {

        ctx->jastrow_champ.factor_ee_gl[nw*shift6
          + k*ctx->electron.num
          + i] =
        ctx->jastrow_champ.factor_ee_gl[nw*shift6
          + k*ctx->electron.num
          + i] +
        ctx->single_point.delta_ee_gl[nw*shift6
          + i*4
          + k];

        ctx->jastrow_champ.factor_en_gl[nw*shift6
          + k*ctx->electron.num
          + i] =
        ctx->jastrow_champ.factor_en_gl[nw*shift6
          + k*ctx->electron.num
          + i] +
        ctx->single_point.delta_en_gl[nw*shift6
          + i*4
          + k];
      }
    }
  }

  for (int nw = 0; nw < ctx->electron.walker.num; nw++) {
    for (int k = 0; k < 3; k++) {
      ctx->point.coord.data[k*ctx->electron.walker.num*ctx->electron.num + nw*ctx->electron.num + ctx->single_point.num] = ctx->single_point.coord.data[nw*3 + k];
    }
  }

  rc = qmckl_context_touch(context);
  if (rc != QMCKL_SUCCESS) return rc;


  if (ctx->jastrow_champ.cord_num > 0){
    ctx->jastrow_champ.een_rescaled_e_date = ctx->date;
    ctx->jastrow_champ.een_rescaled_e_gl_date = ctx->date;
    ctx->jastrow_champ.een_rescaled_n_date = ctx->date;
    ctx->jastrow_champ.een_rescaled_n_gl_date = ctx->date;
    ctx->jastrow_champ.factor_een_date = ctx->date;
    ctx->jastrow_champ.factor_een_gl_date = ctx->date;

  }
  ctx->jastrow_champ.ee_distance_rescaled_date = ctx->date;
  ctx->jastrow_champ.ee_distance_rescaled_gl_date = ctx->date;
  ctx->jastrow_champ.en_distance_rescaled_date = ctx->date;
  ctx->jastrow_champ.en_distance_rescaled_gl_date = ctx->date;
  ctx->jastrow_champ.factor_ee_date = ctx->date;
  ctx->jastrow_champ.factor_ee_gl_date = ctx->date;
  ctx->jastrow_champ.factor_en_date = ctx->date;
  ctx->jastrow_champ.factor_en_gl_date = ctx->date;

  ctx->electron.ee_distance_date = ctx->date;
  ctx->electron.en_distance_date = ctx->date;

  ctx->single_point.date = ctx->date;

  if (ctx->jastrow_champ.cord_num > 0){
    ctx->single_point.een_rescaled_single_e_date = ctx->single_point.date;
    ctx->single_point.een_rescaled_single_n_date = ctx->single_point.date;
    ctx->single_point.delta_een_date = ctx->single_point.date;
    ctx->single_point.een_rescaled_single_e_gl_date = ctx->single_point.date;
    ctx->single_point.een_rescaled_single_n_gl_date = ctx->single_point.date;
    ctx->single_point.delta_een_gl_date = ctx->single_point.date;
    ctx->single_point.delta_een_g_date = ctx->single_point.date;
    ctx->forces.forces_jastrow_single_een_date = ctx->single_point.date;
    ctx->single_point.tmp_date = ctx->single_point.date;
  }
  ctx->single_point.single_ee_distance_date = ctx->single_point.date;
  ctx->single_point.single_en_distance_date = ctx->single_point.date;
  ctx->single_point.ee_rescaled_single_date = ctx->single_point.date;
  ctx->single_point.en_rescaled_single_date = ctx->single_point.date;
  ctx->single_point.delta_en_date = ctx->single_point.date;
  ctx->single_point.delta_ee_date = ctx->single_point.date;
  ctx->single_point.ee_rescaled_single_gl_date = ctx->single_point.date;
  ctx->single_point.en_rescaled_single_gl_date = ctx->single_point.date;
  ctx->single_point.delta_en_gl_date = ctx->single_point.date;
  ctx->single_point.delta_ee_gl_date = ctx->single_point.date;
  ctx->forces.forces_jastrow_single_en_date = ctx->single_point.date;


  return QMCKL_SUCCESS;
}
