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

#include "qmckl.h"
#include "qmckl_context_private_type.h"
#include "qmckl_memory_private_type.h"
#include "qmckl_memory_private_func.h"
#include "qmckl_ao_private_type.h"
#include "qmckl_ao_private_func.h"
#include "qmckl_mo_private_type.h"
#include "qmckl_mo_private_func.h"

qmckl_exit_code qmckl_init_mo_basis(qmckl_context context) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return false;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  ctx->mo_basis.r_cusp = NULL;

  ctx->mo_basis.uninitialized = (1 << 2) - 1;

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_set_mo_basis_mo_num(qmckl_context context, const int64_t mo_num) {

  int32_t mask = 1 ;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
   }
  
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  
  if (mask != 0 && !(ctx->mo_basis.uninitialized & mask)) {
      return qmckl_failwith( context,
                             QMCKL_ALREADY_SET,
                             "qmckl_set_mo_*",
                             NULL);
   }

  if (mo_num <= 0) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_set_mo_basis_mo_num",
                           "mo_num <= 0");
  }

  ctx->mo_basis.mo_num = mo_num;

  ctx->mo_basis.uninitialized &= ~mask;
  ctx->mo_basis.provided = (ctx->mo_basis.uninitialized == 0);
  if (ctx->mo_basis.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_mo_basis(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
  }
  return QMCKL_SUCCESS;

}

qmckl_exit_code  qmckl_set_mo_basis_coefficient(qmckl_context context,
                                                const double* coefficient,
                                                const int64_t size_max)
{
  int32_t mask = 1 << 1;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
     return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;

  if (ctx->mo_basis.coefficient != NULL) {
    qmckl_exit_code rc = qmckl_free(context, ctx->mo_basis.coefficient);
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context, rc,
                             "qmckl_set_mo_basis_coefficient",
                             NULL);
    }
  }

  if (size_max < ctx->ao_basis.ao_num * ctx->mo_basis.mo_num) {
    return qmckl_failwith( context, QMCKL_INVALID_ARG_3,
                             "qmckl_set_mo_basis_coefficient",
                             "Array too small");
  }

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = ctx->ao_basis.ao_num * ctx->mo_basis.mo_num * sizeof(double);
  double* new_array = (double*) qmckl_malloc(context, mem_info);
  if (new_array == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_set_mo_basis_coefficient",
                           NULL);
  }

  memcpy(new_array, coefficient, mem_info.size);

  ctx->mo_basis.coefficient = new_array;

  ctx->mo_basis.uninitialized &= ~mask;
  ctx->mo_basis.provided = (ctx->mo_basis.uninitialized == 0);
  if (ctx->mo_basis.provided) {
    qmckl_exit_code rc_ = qmckl_finalize_mo_basis(context);
    if (rc_ != QMCKL_SUCCESS) return rc_;
  }
  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_finalize_mo_basis(qmckl_context context) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_finalize_mo_basis",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (ctx->mo_basis.coefficient_t != NULL) {
    qmckl_exit_code rc = qmckl_free(context, ctx->mo_basis.coefficient_t);
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context, rc,
                             "qmckl_finalize_mo_basis",
                             NULL);
    }
  }

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = ctx->ao_basis.ao_num * ctx->mo_basis.mo_num * sizeof(double);
  double* new_array = (double*) qmckl_malloc(context, mem_info);
  if (new_array == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_finalize_mo_basis",
                           NULL);
  }

  assert (ctx->mo_basis.coefficient != NULL);

  for (int64_t i=0 ; i<ctx->ao_basis.ao_num ; ++i) {
    for (int64_t j=0 ; j<ctx->mo_basis.mo_num ; ++j) {
      new_array[i*ctx->mo_basis.mo_num + j] = ctx->mo_basis.coefficient[j*ctx->ao_basis.ao_num + i];
    }
  }

  ctx->mo_basis.coefficient_t = new_array;

  qmckl_exit_code rc;
  if (ctx->mo_basis.mo_vgl != NULL) {
    rc = qmckl_free(context, ctx->mo_basis.mo_vgl);
    if (rc != QMCKL_SUCCESS) return rc;
    ctx->mo_basis.mo_vgl = NULL;
    ctx->mo_basis.mo_vgl_date = 0;
  }

  if (ctx->mo_basis.mo_value != NULL) {
    rc = qmckl_free(context, ctx->mo_basis.mo_value);
    if (rc != QMCKL_SUCCESS) return rc;
    ctx->mo_basis.mo_value = NULL;
    ctx->mo_basis.mo_value_date = 0;
  }

  return qmckl_context_touch(context);
}



/* Performs a deep copy of the MO basis structure from ~src~ to ~dest~. */
/* Memory allocations are done using the provided context. */


qmckl_exit_code qmckl_copy_mo_basis(qmckl_context context, const qmckl_mo_basis_struct* src, qmckl_mo_basis_struct* dest) {
  
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  }

  if (src == NULL || dest == NULL) {
    return qmckl_failwith(context,
                          QMCKL_INVALID_ARG_2,
                          "qmckl_copy_mo_basis",
                          "NULL pointer");
  }

  /* If source is not provided/initialized, just copy the struct */
  if (!src->provided) {
    *dest = *src;
    return QMCKL_SUCCESS;
  }

  /* Copy scalar fields */
  dest->mo_num = src->mo_num;
  dest->mo_vgl_date = src->mo_vgl_date;
  dest->mo_value_date = src->mo_value_date;
  dest->uninitialized = src->uninitialized;
  dest->provided = src->provided;

  /* Get ao_num from the context for array size calculations */
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  int64_t ao_num = ctx->ao_basis.ao_num;
  int64_t mo_num = src->mo_num;

  /* Deep copy coefficient array */
  if (src->coefficient != NULL && ao_num > 0 && mo_num > 0) {
    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = ao_num * mo_num * sizeof(double);
    dest->coefficient = (double*) qmckl_malloc(context, mem_info);
    if (dest->coefficient == NULL) return QMCKL_ALLOCATION_FAILED;
    memcpy(dest->coefficient, src->coefficient, mem_info.size);
  } else {
    dest->coefficient = NULL;
  }

  /* Deep copy coefficient_t array */
  if (src->coefficient_t != NULL && ao_num > 0 && mo_num > 0) {
    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = ao_num * mo_num * sizeof(double);
    dest->coefficient_t = (double*) qmckl_malloc(context, mem_info);
    if (dest->coefficient_t == NULL) return QMCKL_ALLOCATION_FAILED;
    memcpy(dest->coefficient_t, src->coefficient_t, mem_info.size);
  } else {
    dest->coefficient_t = NULL;
  }

  /* Deep copy r_cusp array */
  if (src->r_cusp != NULL && ctx->nucleus.num > 0) {
    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = ctx->nucleus.num * sizeof(double);
    dest->r_cusp = (double*) qmckl_malloc(context, mem_info);
    if (dest->r_cusp == NULL) return QMCKL_ALLOCATION_FAILED;
    memcpy(dest->r_cusp, src->r_cusp, mem_info.size);
  } else {
    dest->r_cusp = NULL;
  }

  /* Deep copy mo_vgl array - this is a computed array, set to NULL to recompute */
  dest->mo_vgl = NULL;
  dest->mo_vgl_date = 0;

  /* Deep copy mo_value array - this is a computed array, set to NULL to recompute */
  dest->mo_value = NULL;
  dest->mo_value_date = 0;

  /* Deep copy cusp_param tensor */
  if (src->cusp_param.data != NULL && src->cusp_param.size[0] > 0) {
    dest->cusp_param = qmckl_tensor_alloc(context, src->cusp_param.order, src->cusp_param.size);
    if (dest->cusp_param.data == NULL) return QMCKL_ALLOCATION_FAILED;
    
    /* Calculate total size of tensor data */
    int64_t tensor_size = 1;
    for (int32_t i = 0; i < src->cusp_param.order; ++i) {
      tensor_size *= src->cusp_param.size[i];
    }
    memcpy(dest->cusp_param.data, src->cusp_param.data, tensor_size * sizeof(double));
  } else {
    memset(&(dest->cusp_param), 0, sizeof(qmckl_tensor));
  }

  return QMCKL_SUCCESS;
}

/* Cusp adjsutment functions */

/*    To activate the cusp adjustment, the user must enter the radius of */
/*    the fitting function for each atom. */

/*    This function requires the computation of the value and gradients */
/*    of the $s$ AOs at the distance equal to the radius, and the values */
/*    of the non-$s$ AOs at the center. */


qmckl_exit_code
qmckl_set_mo_basis_r_cusp(qmckl_context context,
                          const double* r_cusp,
                          const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;

  if (r_cusp == NULL) {
      return qmckl_failwith( context, QMCKL_INVALID_ARG_2,
                             "qmckl_set_mo_basis_r_cusp",
                             "r_cusp: Null pointer");
  }

  if (size_max < ctx->nucleus.num) {
    return qmckl_failwith( context, QMCKL_INVALID_ARG_3,
                             "qmckl_set_mo_basis_r_cusp",
                             "Array too small");
  }


  // Nullify r_cusp
  if (ctx->mo_basis.r_cusp != NULL) {
    rc = qmckl_free(context, ctx->mo_basis.r_cusp);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }
    ctx->mo_basis.r_cusp = NULL;
  }


  // Save old points
  int64_t old_point_num = ctx->point.num;
  double* old_coord = NULL;
  if (old_point_num > 0) {
    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = old_point_num * 3 * sizeof(double);
    old_coord = (double*) qmckl_malloc(context, mem_info);

    rc = qmckl_get_point(context, 'T', old_coord,  (old_point_num * 3));
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context,
                             QMCKL_FAILURE,
                             "qmckl_set_mo_basis_r_cusp",
                             "Unable to get old coordinates");
    }
  }

  double* coord = (double*) malloc(ctx->nucleus.num * 3 * sizeof(double));

  // Set r_cusp
  {
    assert (ctx->mo_basis.r_cusp == NULL);
    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = ctx->nucleus.num * sizeof(double);
    ctx->mo_basis.r_cusp = (double*) qmckl_malloc(context, mem_info);
    if (ctx->mo_basis.r_cusp == NULL) {
      return qmckl_failwith( context,
                             QMCKL_ALLOCATION_FAILED,
                             "qmckl_set_mo_basis_r_cusp",
                             NULL);
    }
    memcpy(ctx->mo_basis.r_cusp, r_cusp, mem_info.size);
  }



  // Allocate cusp parameters and set them to zero
  {
    if (ctx->mo_basis.cusp_param.size[0] != 0) {
      rc = qmckl_tensor_free(context, &(ctx->mo_basis.cusp_param));
      if (rc != QMCKL_SUCCESS) return rc;
    }

    int64_t sze[3] = { ctx->mo_basis.mo_num, 4, ctx->nucleus.num };
    ctx->mo_basis.cusp_param = qmckl_tensor_alloc(context, 3, &(sze[0]));
    ctx->mo_basis.cusp_param = qmckl_tensor_set(ctx->mo_basis.cusp_param, 0.);
  }


  // Evaluate MO value at nucleus without s components
  qmckl_matrix mo_value_at_nucl_no_s;
  {
    mo_value_at_nucl_no_s = qmckl_matrix_alloc(context, ctx->mo_basis.mo_num, ctx->nucleus.num);

    rc = qmckl_double_of_matrix(context, ctx->nucleus.coord, coord, ctx->nucleus.num * 3);
    if (rc != QMCKL_SUCCESS) return rc;

    rc = qmckl_set_point(context, 'T', ctx->nucleus.num, coord, (ctx->nucleus.num * 3));
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context,
                             QMCKL_FAILURE,
                             "qmckl_set_mo_basis_r_cusp",
                             "Unable to set coordinates at the nuclei");
    }

    rc = qmckl_get_mo_basis_mo_value(context,
                                     &(qmckl_mat(mo_value_at_nucl_no_s,0,0)),
                                     ctx->mo_basis.mo_num * ctx->nucleus.num);
    if (rc != QMCKL_SUCCESS) return rc;
  }


  // Evaluate MO vgl at r_cusp without s components
  qmckl_tensor mo_vgl_at_r_cusp_s;
  {
    int64_t sze[3] = { ctx->mo_basis.mo_num, 3, ctx->nucleus.num };
    mo_vgl_at_r_cusp_s = qmckl_tensor_alloc(context, 3, &(sze[0]));
  }

  {
    qmckl_tensor ao_vgl_at_r_cusp_s;
    int64_t sze[3] = { ctx->ao_basis.ao_num, 5, ctx->nucleus.num };
    ao_vgl_at_r_cusp_s = qmckl_tensor_alloc(context, 3, &(sze[0]));

    rc = qmckl_double_of_matrix(context, ctx->nucleus.coord, coord, ctx->nucleus.num * 3);
    if (rc != QMCKL_SUCCESS) return rc;

    for (int64_t i=0 ; i<ctx->nucleus.num ; ++i) {
      coord[2*ctx->nucleus.num + i] += r_cusp[i];
    }

    rc = qmckl_set_point(context, 'T', ctx->nucleus.num, coord, (ctx->nucleus.num * 3));
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context,
                             QMCKL_FAILURE,
                             "qmckl_set_mo_basis_r_cusp",
                             "Unable to set coordinates at r_cusp");
    }

    rc = qmckl_get_ao_basis_ao_vgl(context,
                                   &(qmckl_ten3(ao_vgl_at_r_cusp_s,0,0,0)),
                                   ctx->ao_basis.ao_num * 5 * ctx->point.num);

    for (int64_t inucl=0 ; inucl<ctx->nucleus.num ; ++inucl) {
      for (int64_t i=0 ; i<ctx->mo_basis.mo_num ; ++i) {
        qmckl_ten3(mo_vgl_at_r_cusp_s,i,0,inucl) = 0.;
        qmckl_ten3(mo_vgl_at_r_cusp_s,i,1,inucl) = 0.;
        qmckl_ten3(mo_vgl_at_r_cusp_s,i,2,inucl) = 0.;
        for (int64_t k=0 ; k<ctx->ao_basis.ao_num ; ++k) {
          if ( ctx->ao_basis.ao_nucl[k] == inucl && ctx->ao_basis.ao_ang_mom[k] == 0) {
            const double ck = ctx->mo_basis.coefficient[k + i*ctx->ao_basis.ao_num];
            qmckl_ten3(mo_vgl_at_r_cusp_s,i,0,inucl) = qmckl_ten3(mo_vgl_at_r_cusp_s,i,0,inucl) + ck * qmckl_ten3(ao_vgl_at_r_cusp_s,k,0,inucl);
            qmckl_ten3(mo_vgl_at_r_cusp_s,i,1,inucl) = qmckl_ten3(mo_vgl_at_r_cusp_s,i,1,inucl) + ck * qmckl_ten3(ao_vgl_at_r_cusp_s,k,3,inucl);
            qmckl_ten3(mo_vgl_at_r_cusp_s,i,2,inucl) = qmckl_ten3(mo_vgl_at_r_cusp_s,i,2,inucl) + ck * qmckl_ten3(ao_vgl_at_r_cusp_s,k,4,inucl);
          }
        }
      }
    }
    rc = qmckl_tensor_free(context,&ao_vgl_at_r_cusp_s);

    if (rc != QMCKL_SUCCESS) return rc;
  }

  // Compute parameters
  {
     for (int64_t inucl=0 ; inucl < ctx->nucleus.num ; ++inucl) {
       const double Z = qmckl_vec(ctx->nucleus.charge,inucl);
       if (Z < 0.1) continue;  // Avoid dummy atoms
       const double R = r_cusp[inucl];
       assert (R != 0.);
       assert (fabs(2.*R*Z-6.) > 1.e-6);
       assert (fabs(R*R*Z-3.*R) > 1.e-6);
       assert (fabs(R*R*R*Z-3.*R*R) > 1.e-6);
       for (int64_t i=0 ; i<ctx->mo_basis.mo_num ; ++i) {
         const double phi        = qmckl_ten3(mo_vgl_at_r_cusp_s,i,0,inucl);
         const double grad_phi   = qmckl_ten3(mo_vgl_at_r_cusp_s,i,1,inucl);
         const double lap_phi    = qmckl_ten3(mo_vgl_at_r_cusp_s,i,2,inucl);
         const double eta        = qmckl_mat(mo_value_at_nucl_no_s,i,inucl);

         qmckl_ten3(ctx->mo_basis.cusp_param,i,0,inucl) =
           -(R*(2.*eta*Z-6.*grad_phi)+lap_phi*R*R+6.*phi)/(2.*R*Z-6.);

         qmckl_ten3(ctx->mo_basis.cusp_param,i,1,inucl) =
           (lap_phi*R*R*Z-6.*grad_phi*R*Z+6.*phi*Z+6.*eta*Z)/(2.*R*Z-6.);

         qmckl_ten3(ctx->mo_basis.cusp_param,i,2,inucl) =
           -(R*(-5.*grad_phi*Z-1.5*lap_phi)+lap_phi*R*R*Z+3.*phi*Z+3.*eta*Z+6.*grad_phi)/(R*R*Z-3.*R);

         qmckl_ten3(ctx->mo_basis.cusp_param,i,3,inucl) =
           (R*(-2.*grad_phi*Z-lap_phi)+0.5*lap_phi*R*R*Z+phi*Z+eta*Z+3.*grad_phi)/(R*R*R*Z-3.*R*R);

       }
     }
  }

  free(coord);
  qmckl_matrix_free(context, &mo_value_at_nucl_no_s);
  qmckl_tensor_free(context, &mo_vgl_at_r_cusp_s);

  // Restore old points
  if (old_point_num > 0) {
    rc = qmckl_set_point(context, 'T', old_point_num, old_coord, (old_point_num * 3));
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context,
                             QMCKL_FAILURE,
                             "qmckl_set_mo_basis_r_cusp",
                             "Unable to set old coordinates");
    }
    rc = qmckl_free(context, old_coord);
    if (rc != QMCKL_SUCCESS) return rc;
    old_coord = NULL;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_mo_basis_mo_num (const qmckl_context context,
                           int64_t* mo_num)
{
   if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_mo_basis_mo_num",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1;

  if ( (ctx->mo_basis.uninitialized & mask) != 0) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_get_mo_basis_mo_num",
                           NULL);
  }

  assert (ctx->mo_basis.mo_num > (int64_t) 0);
  *mo_num = ctx->mo_basis.mo_num;
  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_mo_basis_coefficient (const qmckl_context context,
                                double* const coefficient,
                                const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_mo_basis_coefficient",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int32_t mask = 1 << 1;

  if ( (ctx->mo_basis.uninitialized & mask) != 0) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_get_mo_basis_coefficient",
                           NULL);
  }

  if (coefficient == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_mo_basis_coefficient",
                           "NULL pointer");
  }

  if (size_max < ctx->ao_basis.ao_num * ctx->mo_basis.mo_num) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_mo_basis_coefficient",
                           "Array too small: expected mo_num * ao_num.");
  }

  assert (ctx->mo_basis.coefficient != NULL);
  memcpy(coefficient, ctx->mo_basis.coefficient,
         ctx->ao_basis.ao_num * ctx->mo_basis.mo_num * sizeof(double));

  return QMCKL_SUCCESS;
}

bool qmckl_mo_basis_provided(const qmckl_context context) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return false;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  return ctx->mo_basis.provided;
}

qmckl_exit_code
qmckl_mo_basis_select_mo (const qmckl_context context,
                          const int32_t* keep,
                          const int64_t size_max)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_mo_basis_select_mo",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if ( !(qmckl_mo_basis_provided(context)) ) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_get_mo_basis_select_mo",
                           NULL);
  }

  if (keep == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_get_mo_basis_select_mo",
                           "NULL pointer");
  }

  const int64_t mo_num = ctx->mo_basis.mo_num;
  const int64_t ao_num = ctx->ao_basis.ao_num;

  if (size_max < mo_num) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_mo_basis_select_mo",
                           "Array too small: expected mo_num.");
  }

/*
  if (ctx->mo_basis.r_cusp != NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALREADY_SET,
                           "qmckl_get_mo_basis_select_mo",
                           "r_cusp is already set. Please set it only after selecting MOs.");
  }
*/

  int64_t mo_num_new = 0;
  for (int64_t i=0 ; i<mo_num ; ++i) {
    if (keep[i] != 0) ++mo_num_new;
  }

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = ao_num * mo_num_new * sizeof(double);
  double* restrict coefficient   = (double*) qmckl_malloc(context, mem_info);

  int64_t k = 0;
  for (int64_t i=0 ; i<mo_num ; ++i) {
    if (keep[i] != 0) {
      memcpy(&(coefficient[k*ao_num]), &(ctx->mo_basis.coefficient[i*ao_num]), ao_num*sizeof(double));
      ++k;
    }
  }

  qmckl_exit_code rc = qmckl_free(context, ctx->mo_basis.coefficient);
  if (rc != QMCKL_SUCCESS) return rc;

  ctx->mo_basis.coefficient = coefficient;
  ctx->mo_basis.mo_num = mo_num_new;

  rc = qmckl_finalize_mo_basis(context);

  if (ctx->mo_basis.r_cusp != NULL) {
    double * tmp = (double*) calloc (ctx->nucleus.num, sizeof(double));
    assert(tmp != NULL);
    memcpy(tmp, ctx->mo_basis.r_cusp, ctx->nucleus.num * sizeof(double));
    rc = qmckl_set_mo_basis_r_cusp(context, tmp,
        mo_num_new * 4 * ctx->nucleus.num);
    free(tmp);
  }

  return rc;
}

qmckl_exit_code
qmckl_get_mo_basis_mo_value(qmckl_context context,
                            double* const mo_value,
                            const int64_t size_max)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_mo_basis_mo_value(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  const int64_t sze = ctx->point.num * ctx->mo_basis.mo_num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_mo_basis_mo_value",
                           "input array too small");
  }
  memcpy(mo_value, ctx->mo_basis.mo_value, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_mo_basis_mo_value_inplace (qmckl_context context,
                                     double* const mo_value,
                                     const int64_t size_max)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_mo_basis_mo_value",
                           NULL);
  }

  qmckl_exit_code rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  const int64_t sze = ctx->mo_basis.mo_num * ctx->point.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_mo_basis_mo_value",
                           "input array too small");
  }

  ctx->mo_basis.mo_value_date = ctx->point.date - 1UL;

  double* old_array = ctx->mo_basis.mo_value;

  ctx->mo_basis.mo_value = mo_value;

  rc = qmckl_provide_mo_basis_mo_value(context);
  if (rc != QMCKL_SUCCESS) return rc;

  ctx->mo_basis.mo_value = old_array;

  ctx->mo_basis.mo_value_date = ctx->point.date - 1UL;

  return QMCKL_SUCCESS;
}



/* #+CALL: write_provider_pre( group="mo_basis", data="mo_value", dimension="ctx->mo_basis.mo_num * ctx->point.num") */

/* #+RESULTS: */

qmckl_exit_code qmckl_provide_mo_basis_mo_value(qmckl_context context)
{

  qmckl_exit_code rc = QMCKL_SUCCESS;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_provide_mo_basis_mo_value",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->mo_basis.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_provide_mo_basis_mo_value",
                           NULL);
  }

  /* Compute if necessary */
  if (ctx->point.date > ctx->mo_basis.mo_value_date) {

    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = ctx->mo_basis.mo_num * ctx->point.num * sizeof(double);

    if (ctx->mo_basis.mo_value != NULL) {
      qmckl_memory_info_struct mem_info_test = qmckl_memory_info_struct_zero;
      rc = qmckl_get_malloc_info(context, ctx->mo_basis.mo_value, &mem_info_test);

      /* if rc != QMCKL_SUCCESS, we are maybe in an _inplace function because the
         memory was not allocated with qmckl_malloc */

      if ((rc == QMCKL_SUCCESS) && (mem_info_test.size != mem_info.size)) {
        rc = qmckl_free(context, ctx->mo_basis.mo_value);
        assert (rc == QMCKL_SUCCESS);
        ctx->mo_basis.mo_value = NULL;
      }
    }

    /* Allocate array */
    if (ctx->mo_basis.mo_value == NULL) {

      double* mo_value = (double*) qmckl_malloc(context, mem_info);

      if (mo_value == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_mo_basis_mo_value",
                               NULL);
      }
      ctx->mo_basis.mo_value = mo_value;
    }

if (ctx->mo_basis.mo_vgl_date == ctx->point.date) {

  // mo_vgl has been computed at this step: Just copy the data.

  double * v = &(ctx->mo_basis.mo_value[0]);
  double * vgl = &(ctx->mo_basis.mo_vgl[0]);
  for (int i=0 ; i<ctx->point.num ; ++i) {
    for (int k=0 ; k<ctx->mo_basis.mo_num ; ++k) {
      v[k] = vgl[k];
    }
    v   += ctx->mo_basis.mo_num;
    vgl += ctx->mo_basis.mo_num * 5;
  }

} else {

  rc = qmckl_provide_ao_basis_ao_value(context);
  if (rc != QMCKL_SUCCESS) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_ao_value",
                           NULL);
  }

  if (ctx->mo_basis.r_cusp == NULL) {
    /* No cusp correction */
    rc = qmckl_compute_mo_basis_mo_value(context,
                                         ctx->ao_basis.ao_num,
                                         ctx->mo_basis.mo_num,
                                         ctx->point.num,
                                         ctx->mo_basis.coefficient_t,
                                         ctx->ao_basis.ao_value,
                                         ctx->mo_basis.mo_value);
  } else {
    rc = qmckl_provide_en_distance(context);
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context,
                             QMCKL_NOT_PROVIDED,
                             "qmckl_provide_mo_basis_mo_value",
                             "en_distance");
    }

    rc = qmckl_compute_mo_basis_mo_value_cusp(context,
                                              ctx->nucleus.num,
                                              ctx->ao_basis.ao_num,
                                              ctx->mo_basis.mo_num,
                                              ctx->point.num,
                                              ctx->ao_basis.ao_nucl,
                                              ctx->ao_basis.ao_ang_mom,
                                              ctx->electron.en_distance,
                                              ctx->mo_basis.r_cusp,
                                              ctx->mo_basis.cusp_param,
                                              ctx->mo_basis.coefficient_t,
                                              ctx->ao_basis.ao_value,
                                              ctx->mo_basis.mo_value);
  }
}





/* #+CALL: write_provider_post( group="mo_basis", data="mo_value" ) */

/* #+RESULTS: */

if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->mo_basis.mo_value_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_compute_mo_basis_mo_value (const qmckl_context context,
                                 const int64_t ao_num,
                                 const int64_t mo_num,
                                 const int64_t point_num,
                                 const double* coefficient_t,
                                 const double* ao_value,
                                 double* const mo_value )
{
#ifdef HAVE_HPC
  return qmckl_compute_mo_basis_mo_value_hpc (context, ao_num, mo_num, point_num, coefficient_t, ao_value, mo_value);
#else
  return qmckl_compute_mo_basis_mo_value_doc (context, ao_num, mo_num, point_num, coefficient_t, ao_value, mo_value);
#endif
}

/* Single-precision */


#ifdef HAVE_HPC
qmckl_exit_code
qmckl_compute_mo_basis_mo_value_hpc_sp (const qmckl_context context,
                                        const int64_t ao_num,
                                        const int64_t mo_num,
                                        const int64_t point_num,
                                        const double* restrict coefficient_t,
                                        const double* restrict ao_value,
                                        double* restrict const mo_value )
{
  assert (context != QMCKL_NULL_CONTEXT);

  float* __attribute__((aligned(64))) coefficient_t_sp = calloc(ao_num*mo_num, sizeof(float));
  if (coefficient_t_sp == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_compute_mo_basis_mo_value_hpc_sp",
                           "coefficient_t_sp");
  };

#ifdef HAVE_OPENMP
#pragma omp simd
#endif
    for (int64_t i=0 ; i<mo_num*ao_num ; ++i) {
      coefficient_t_sp[i] = (float) coefficient_t[i];
    }

#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
    {
      int64_t* __attribute__((aligned(64))) idx = calloc((size_t) ao_num, sizeof(int64_t));
      float*   __attribute__((aligned(64))) av1 = calloc((size_t) ao_num, sizeof(float));
      assert (idx != NULL);
      assert (av1 != NULL);

#ifdef HAVE_OPENMP
#pragma omp for
#endif
    for (int64_t ipoint=0 ; ipoint < point_num ; ++ipoint) {
      double* restrict const vgl1  = &(mo_value[ipoint*mo_num]);
      const double* restrict avgl1 = &(ao_value[ipoint*ao_num]);

      float*   __attribute__((aligned(64))) vgl_sp = calloc((size_t) mo_num, sizeof(float));
      assert (vgl_sp != NULL);

      int64_t nidx=0;
      for (int64_t k=0 ; k<ao_num ; ++k) {
        if (avgl1[k] != 0.) {
          idx[nidx] = k;
          av1[nidx] = (float) avgl1[k];
          ++nidx;
        }
      }

      int64_t n=0;

      for (n=0 ; n < nidx-8 ; n+=8) {
        const float* restrict ck1 = coefficient_t_sp + idx[n  ]*mo_num;
        const float* restrict ck2 = coefficient_t_sp + idx[n+1]*mo_num;
        const float* restrict ck3 = coefficient_t_sp + idx[n+2]*mo_num;
        const float* restrict ck4 = coefficient_t_sp + idx[n+3]*mo_num;
        const float* restrict ck5 = coefficient_t_sp + idx[n+4]*mo_num;
        const float* restrict ck6 = coefficient_t_sp + idx[n+5]*mo_num;
        const float* restrict ck7 = coefficient_t_sp + idx[n+6]*mo_num;
        const float* restrict ck8 = coefficient_t_sp + idx[n+7]*mo_num;

        const float a11 = av1[n  ];
        const float a21 = av1[n+1];
        const float a31 = av1[n+2];
        const float a41 = av1[n+3];
        const float a51 = av1[n+4];
        const float a61 = av1[n+5];
        const float a71 = av1[n+6];
        const float a81 = av1[n+7];

#ifdef HAVE_OPENMP
#pragma omp simd
#endif
          for (int64_t i=0 ; i<mo_num ; ++i) {
            vgl_sp[i] = vgl_sp[i] +
                         ck1[i] * a11 + ck2[i] * a21 + ck3[i] * a31 + ck4[i] * a41 +
                         ck5[i] * a51 + ck6[i] * a61 + ck7[i] * a71 + ck8[i] * a81;
          }
      }

      for (int64_t m=n ; m < nidx ; m+=1) {
        const float* restrict ck = coefficient_t_sp + idx[m]*mo_num;
        const float a1 = av1[m];

#ifdef HAVE_OPENMP
#pragma omp simd
#endif
          for (int64_t i=0 ; i<mo_num ; ++i) {
            vgl_sp[i] = vgl_sp[i] + ck[i] * a1;
          }
      }

#ifdef HAVE_OPENMP
#pragma omp simd
#endif
      for (int64_t i=0 ; i<mo_num ; ++i) {
        vgl1[i] = (double) vgl_sp[i];
      }
      free(vgl_sp);
    }
    free(av1);
    free(idx);
  }
  free(coefficient_t_sp);
  return QMCKL_SUCCESS;
}
#endif

/* Double-precision */


#ifdef HAVE_HPC
qmckl_exit_code
qmckl_compute_mo_basis_mo_value_hpc (const qmckl_context context,
                                     const int64_t ao_num,
                                     const int64_t mo_num,
                                     const int64_t point_num,
                                     const double* restrict coefficient_t,
                                     const double* restrict ao_value,
                                     double* restrict const mo_value )
{
  assert (context != QMCKL_NULL_CONTEXT);
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  /* Don't compute polynomials when the radial part is zero. */
  const int precision = ctx->numprec.precision;
  const bool single_precision = precision <= 24;

  if (single_precision) {
    return qmckl_compute_mo_basis_mo_value_hpc_sp (context,
                                                   ao_num,
                                                   mo_num,
                                                   point_num,
                                                   coefficient_t,
                                                   ao_value,
                                                   mo_value );
  }

#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
  {
    int64_t* __attribute__((aligned(64))) idx = calloc(ao_num, sizeof(int64_t));
    double*  __attribute__((aligned(64))) av1 = calloc(ao_num, sizeof(double));
    assert (idx != NULL);
    assert (av1 != NULL);

#ifdef HAVE_OPENMP
#pragma omp for
#endif
    for (int64_t ipoint=0 ; ipoint < point_num ; ++ipoint) {
      double* restrict const vgl1  = &(mo_value[ipoint*mo_num]);
      const double* restrict avgl1 = &(ao_value[ipoint*ao_num]);

      memset(vgl1, 0, mo_num*sizeof(double));

      int64_t nidx=0;
      int64_t idx[ao_num];
      double  av1[ao_num];

      for (int64_t k=0 ; k<ao_num ; ++k) {
        if (avgl1[k] != 0.) {
          idx[nidx] = k;
          av1[nidx] = avgl1[k];
          ++nidx;
        }
      }

      int64_t n=0;

      for (n=0 ; n < nidx-8 ; n+=8) {
        const double* restrict ck1 = coefficient_t + idx[n  ]*mo_num;
        const double* restrict ck2 = coefficient_t + idx[n+1]*mo_num;
        const double* restrict ck3 = coefficient_t + idx[n+2]*mo_num;
        const double* restrict ck4 = coefficient_t + idx[n+3]*mo_num;
        const double* restrict ck5 = coefficient_t + idx[n+4]*mo_num;
        const double* restrict ck6 = coefficient_t + idx[n+5]*mo_num;
        const double* restrict ck7 = coefficient_t + idx[n+6]*mo_num;
        const double* restrict ck8 = coefficient_t + idx[n+7]*mo_num;

        const double a11 = av1[n  ];
        const double a21 = av1[n+1];
        const double a31 = av1[n+2];
        const double a41 = av1[n+3];
        const double a51 = av1[n+4];
        const double a61 = av1[n+5];
        const double a71 = av1[n+6];
        const double a81 = av1[n+7];

#ifdef HAVE_OPENMP
#pragma omp simd
#endif
          for (int64_t i=0 ; i<mo_num ; ++i) {
            vgl1[i] = vgl1[i] +
                       ck1[i] * a11 + ck2[i] * a21 + ck3[i] * a31 + ck4[i] * a41 +
                       ck5[i] * a51 + ck6[i] * a61 + ck7[i] * a71 + ck8[i] * a81;
          }
      }

      for (int64_t m=n ; m < nidx ; m+=1) {
        const double* restrict ck = coefficient_t + idx[m]*mo_num;
        const double a1 = av1[m];

#ifdef HAVE_OPENMP
#pragma omp simd
#endif
          for (int64_t i=0 ; i<mo_num ; ++i) {
            vgl1[i] = vgl1[i] + ck[i] * a1;
          }
      }
    }
    free(av1);
    free(idx);
  }
  return QMCKL_SUCCESS;
}
#endif

qmckl_exit_code
qmckl_get_mo_basis_mo_vgl(qmckl_context context,
                          double* const mo_vgl,
                          const int64_t size_max)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_mo_basis_mo_vgl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  const int64_t sze = ctx->point.num * 5 * ctx->mo_basis.mo_num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_mo_basis_mo_vgl",
                           "input array too small");
  }
  memcpy(mo_vgl, ctx->mo_basis.mo_vgl, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_mo_basis_mo_vgl_inplace (qmckl_context context,
                                   double* const mo_vgl,
                                   const int64_t size_max)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_mo_basis_mo_vgl",
                           NULL);
  }

  qmckl_exit_code rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  const int64_t sze = ctx->mo_basis.mo_num * 5 * ctx->point.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_mo_basis_mo_vgl",
                           "input array too small");
  }

  ctx->mo_basis.mo_vgl_date = ctx->point.date - 1UL;

  double* old_array = ctx->mo_basis.mo_vgl;

  ctx->mo_basis.mo_vgl = mo_vgl;

  rc = qmckl_provide_mo_basis_mo_vgl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  ctx->mo_basis.mo_vgl = old_array;

  ctx->mo_basis.mo_vgl_date = ctx->point.date - 1UL;

  return QMCKL_SUCCESS;
}



/* #+CALL: write_provider_pre( group="mo_basis", data="mo_vgl", dimension="5 * ctx->mo_basis.mo_num * ctx->point.num") */

/* #+RESULTS: */

qmckl_exit_code qmckl_provide_mo_basis_mo_vgl(qmckl_context context)
{

  qmckl_exit_code rc = QMCKL_SUCCESS;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_provide_mo_basis_mo_vgl",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->mo_basis.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_provide_mo_basis_mo_vgl",
                           NULL);
  }

  /* Compute if necessary */
  if (ctx->point.date > ctx->mo_basis.mo_vgl_date) {

    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = 5 * ctx->mo_basis.mo_num * ctx->point.num * sizeof(double);

    if (ctx->mo_basis.mo_vgl != NULL) {
      qmckl_memory_info_struct mem_info_test = qmckl_memory_info_struct_zero;
      rc = qmckl_get_malloc_info(context, ctx->mo_basis.mo_vgl, &mem_info_test);

      /* if rc != QMCKL_SUCCESS, we are maybe in an _inplace function because the
         memory was not allocated with qmckl_malloc */

      if ((rc == QMCKL_SUCCESS) && (mem_info_test.size != mem_info.size)) {
        rc = qmckl_free(context, ctx->mo_basis.mo_vgl);
        assert (rc == QMCKL_SUCCESS);
        ctx->mo_basis.mo_vgl = NULL;
      }
    }

    /* Allocate array */
    if (ctx->mo_basis.mo_vgl == NULL) {

      double* mo_vgl = (double*) qmckl_malloc(context, mem_info);

      if (mo_vgl == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_mo_basis_mo_vgl",
                               NULL);
      }
      ctx->mo_basis.mo_vgl = mo_vgl;
    }

rc = qmckl_provide_ao_basis_ao_vgl(context);
if (rc != QMCKL_SUCCESS) {
  return qmckl_failwith( context,
                         QMCKL_NOT_PROVIDED,
                         "qmckl_ao_basis",
                         NULL);
}

if (ctx->mo_basis.r_cusp == NULL) {
  /* No cusp correction */
  rc = qmckl_compute_mo_basis_mo_vgl(context,
                                     ctx->ao_basis.ao_num,
                                     ctx->mo_basis.mo_num,
                                     ctx->point.num,
                                     ctx->mo_basis.coefficient_t,
                                     ctx->ao_basis.ao_vgl,
                                     ctx->mo_basis.mo_vgl);
} else {
  rc = qmckl_provide_en_distance(context);
  if (rc != QMCKL_SUCCESS) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_electron_en_distance",
                           NULL);
  }
  rc = qmckl_compute_mo_basis_mo_vgl_cusp(context,
                                          ctx->nucleus.num,
                                          ctx->ao_basis.ao_num,
                                          ctx->mo_basis.mo_num,
                                          ctx->point.num,
                                          ctx->ao_basis.ao_nucl,
                                          ctx->ao_basis.ao_ang_mom,
                                          ctx->electron.en_distance,
                                          ctx->nucleus.coord,
                                          ctx->point.coord,
                                          ctx->mo_basis.r_cusp,
                                          ctx->mo_basis.cusp_param,
                                          ctx->mo_basis.coefficient_t,
                                          ctx->ao_basis.ao_vgl,
                                          ctx->mo_basis.mo_vgl);
}



/* #+CALL: write_provider_post( group="mo_basis", data="mo_vgl" ) */

/* #+RESULTS: */

if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->mo_basis.mo_vgl_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_compute_mo_basis_mo_vgl (const qmckl_context context,
                            const int64_t ao_num,
                            const int64_t mo_num,
                            const int64_t point_num,
                            const double* coefficient_t,
                            const double* ao_vgl,
                            double* const mo_vgl )
{
#ifdef HAVE_HPC
  return qmckl_compute_mo_basis_mo_vgl_hpc (context, ao_num, mo_num, point_num, coefficient_t, ao_vgl, mo_vgl);
#else
  return qmckl_compute_mo_basis_mo_vgl_doc (context, ao_num, mo_num, point_num, coefficient_t, ao_vgl, mo_vgl);
#endif
}

#ifdef HAVE_HPC
qmckl_exit_code
qmckl_compute_mo_basis_mo_vgl_hpc (const qmckl_context context,
                                   const int64_t ao_num,
                                   const int64_t mo_num,
                                   const int64_t point_num,
                                   const double* restrict coefficient_t,
                                   const double* restrict ao_vgl,
                                   double* restrict const mo_vgl )
{
  assert (context != QMCKL_NULL_CONTEXT);
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  /* Don't compute polynomials when the radial part is zero. */
  const int precision = ctx->numprec.precision;
  const bool single_precision = precision <= 26;

  if (single_precision) {
    return qmckl_compute_mo_basis_mo_vgl_hpc_sp (context,
                                                 ao_num,
                                                 mo_num,
                                                 point_num,
                                                 coefficient_t,
                                                 ao_vgl,
                                                 mo_vgl );
  }

#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
  {
    int64_t* __attribute__((aligned(64))) idx = calloc(ao_num, sizeof(int64_t));
    double*  __attribute__((aligned(64))) av1 = calloc(ao_num, sizeof(double));
    double*  __attribute__((aligned(64))) av2 = calloc(ao_num, sizeof(double));
    double*  __attribute__((aligned(64))) av3 = calloc(ao_num, sizeof(double));
    double*  __attribute__((aligned(64))) av4 = calloc(ao_num, sizeof(double));
    double*  __attribute__((aligned(64))) av5 = calloc(ao_num, sizeof(double));
    assert (idx != NULL);
    assert (av1 != NULL);
    assert (av2 != NULL);
    assert (av3 != NULL);
    assert (av4 != NULL);
    assert (av5 != NULL);

#ifdef HAVE_OPENMP
#pragma omp for
#endif
  for (int64_t ipoint=0 ; ipoint < point_num ; ++ipoint) {
    double* restrict const vgl1 = &(mo_vgl[ipoint*5*mo_num]);
    double* restrict const vgl2 =  vgl1 + mo_num;
    double* restrict const vgl3 =  vgl1 + (mo_num << 1);
    double* restrict const vgl4 =  vgl1 + (mo_num << 1) + mo_num;
    double* restrict const vgl5 =  vgl1 + (mo_num << 2);

    const double* restrict avgl1 = &(ao_vgl[ipoint*5*ao_num]);
    const double* restrict avgl2 = avgl1 + ao_num;
    const double* restrict avgl3 = avgl1 + (ao_num << 1);
    const double* restrict avgl4 = avgl1 + (ao_num << 1) + ao_num;
    const double* restrict avgl5 = avgl1 + (ao_num << 2);

    memset(vgl1,0,5*mo_num*sizeof(double));

    int64_t nidx=0;
    for (int64_t k=0 ; k<ao_num ; ++k) {
      if (avgl1[k] != 0.) {
        idx[nidx] = k;
        av1[nidx] = avgl1[k];
        av2[nidx] = avgl2[k];
        av3[nidx] = avgl3[k];
        av4[nidx] = avgl4[k];
        av5[nidx] = avgl5[k];
        ++nidx;
      }
    }

    int64_t n=0;

    for (n=0 ; n < nidx-4 ; n+=4) {
      const double* restrict ck1 = coefficient_t + idx[n  ]*mo_num;
      const double* restrict ck2 = coefficient_t + idx[n+1]*mo_num;
      const double* restrict ck3 = coefficient_t + idx[n+2]*mo_num;
      const double* restrict ck4 = coefficient_t + idx[n+3]*mo_num;

      const double a11 = av1[n  ];
      const double a21 = av1[n+1];
      const double a31 = av1[n+2];
      const double a41 = av1[n+3];

      const double a12 = av2[n  ];
      const double a22 = av2[n+1];
      const double a32 = av2[n+2];
      const double a42 = av2[n+3];

      const double a13 = av3[n  ];
      const double a23 = av3[n+1];
      const double a33 = av3[n+2];
      const double a43 = av3[n+3];

      const double a14 = av4[n  ];
      const double a24 = av4[n+1];
      const double a34 = av4[n+2];
      const double a44 = av4[n+3];

      const double a15 = av5[n  ];
      const double a25 = av5[n+1];
      const double a35 = av5[n+2];
      const double a45 = av5[n+3];

#ifdef HAVE_OPENMP
#pragma omp simd
#endif
      for (int64_t i=0 ; i<mo_num ; ++i) {
        vgl1[i] = vgl1[i] + ck1[i] * a11 + ck2[i] * a21 + ck3[i] * a31 + ck4[i] * a41;
        vgl2[i] = vgl2[i] + ck1[i] * a12 + ck2[i] * a22 + ck3[i] * a32 + ck4[i] * a42;
        vgl3[i] = vgl3[i] + ck1[i] * a13 + ck2[i] * a23 + ck3[i] * a33 + ck4[i] * a43;
        vgl4[i] = vgl4[i] + ck1[i] * a14 + ck2[i] * a24 + ck3[i] * a34 + ck4[i] * a44;
        vgl5[i] = vgl5[i] + ck1[i] * a15 + ck2[i] * a25 + ck3[i] * a35 + ck4[i] * a45;
      }
    }

    for (int64_t m=n ; m < nidx ; m+=1) {
      const double* restrict ck = coefficient_t + idx[m]*mo_num;
      const double a1 = av1[m];
      const double a2 = av2[m];
      const double a3 = av3[m];
      const double a4 = av4[m];
      const double a5 = av5[m];

#ifdef HAVE_OPENMP
  #pragma omp simd
#endif
      for (int64_t i=0 ; i<mo_num ; ++i) {
        vgl1[i] = vgl1[i] + ck[i] * a1;
        vgl2[i] = vgl2[i] + ck[i] * a2;
        vgl3[i] = vgl3[i] + ck[i] * a3;
        vgl4[i] = vgl4[i] + ck[i] * a4;
        vgl5[i] = vgl5[i] + ck[i] * a5;
      }
    }
  }
  free(idx);
  free(av1);
  free(av2);
  free(av3);
  free(av4);
  free(av5);
  }
  return QMCKL_SUCCESS;
}
#endif

#ifdef HAVE_HPC
qmckl_exit_code
qmckl_compute_mo_basis_mo_vgl_hpc_sp (const qmckl_context context,
                                      const int64_t ao_num,
                                      const int64_t mo_num,
                                      const int64_t point_num,
                                      const double* restrict coefficient_t,
                                      const double* restrict ao_vgl,
                                      double* restrict const mo_vgl )
{
  assert (context != QMCKL_NULL_CONTEXT);

  float* __attribute__((aligned(64))) coefficient_t_sp = calloc(ao_num*mo_num, sizeof(float));
  if (coefficient_t_sp == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_compute_mo_basis_mo_vgl_hpc_sp",
                           "coefficient_t_sp");
  };

#ifdef HAVE_OPENMP
#pragma omp simd
#endif
    for (int64_t i=0 ; i<mo_num*ao_num ; ++i) {
      coefficient_t_sp[i] = (float) coefficient_t[i];
    }


#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
  {
    int64_t* __attribute__((aligned(64))) idx = calloc(ao_num, sizeof(int64_t));
    float*  __attribute__((aligned(64))) av1 = calloc(ao_num, sizeof(double));
    float*  __attribute__((aligned(64))) av2 = calloc(ao_num, sizeof(double));
    float*  __attribute__((aligned(64))) av3 = calloc(ao_num, sizeof(double));
    float*  __attribute__((aligned(64))) av4 = calloc(ao_num, sizeof(double));
    float*  __attribute__((aligned(64))) av5 = calloc(ao_num, sizeof(double));
    assert (idx != NULL);
    assert (av1 != NULL);
    assert (av2 != NULL);
    assert (av3 != NULL);
    assert (av4 != NULL);
    assert (av5 != NULL);

#ifdef HAVE_OPENMP
#pragma omp for
#endif
    for (int64_t ipoint=0 ; ipoint < point_num ; ++ipoint) {
      float* __attribute__((aligned(64))) vgl_sp1 = calloc(mo_num, sizeof(float));
      float* __attribute__((aligned(64))) vgl_sp2 = calloc(mo_num, sizeof(float));
      float* __attribute__((aligned(64))) vgl_sp3 = calloc(mo_num, sizeof(float));
      float* __attribute__((aligned(64))) vgl_sp4 = calloc(mo_num, sizeof(float));
      float* __attribute__((aligned(64))) vgl_sp5 = calloc(mo_num, sizeof(float));
      assert (vgl_sp1 != NULL);
      assert (vgl_sp2 != NULL);
      assert (vgl_sp3 != NULL);
      assert (vgl_sp4 != NULL);
      assert (vgl_sp5 != NULL);

      double* restrict const vgl1 = &(mo_vgl[ipoint*5*mo_num]);
      double* restrict const vgl2 =  vgl1 + mo_num;
      double* restrict const vgl3 =  vgl1 + (mo_num << 1);
      double* restrict const vgl4 =  vgl1 + (mo_num << 1) + mo_num;
      double* restrict const vgl5 =  vgl1 + (mo_num << 2);

      const double* restrict avgl1 = &(ao_vgl[ipoint*5*ao_num]);
      const double* restrict avgl2 = avgl1 + ao_num;
      const double* restrict avgl3 = avgl1 + (ao_num << 1);
      const double* restrict avgl4 = avgl1 + (ao_num << 1) + ao_num;
      const double* restrict avgl5 = avgl1 + (ao_num << 2);

      int64_t nidx=0;
      for (int64_t k=0 ; k<ao_num ; ++k) {
        if (avgl1[k] != 0.) {
          idx[nidx] = k;
          av1[nidx] = (float) avgl1[k];
          av2[nidx] = (float) avgl2[k];
          av3[nidx] = (float) avgl3[k];
          av4[nidx] = (float) avgl4[k];
          av5[nidx] = (float) avgl5[k];
          ++nidx;
        }
      }

      int64_t n=0;

      for (n=0 ; n < nidx-4 ; n+=4) {
        const float* restrict ck1 = coefficient_t_sp + idx[n  ]*mo_num;
        const float* restrict ck2 = coefficient_t_sp + idx[n+1]*mo_num;
        const float* restrict ck3 = coefficient_t_sp + idx[n+2]*mo_num;
        const float* restrict ck4 = coefficient_t_sp + idx[n+3]*mo_num;

        const float a11 = av1[n  ];
        const float a21 = av1[n+1];
        const float a31 = av1[n+2];
        const float a41 = av1[n+3];

        const float a12 = av2[n  ];
        const float a22 = av2[n+1];
        const float a32 = av2[n+2];
        const float a42 = av2[n+3];

        const float a13 = av3[n  ];
        const float a23 = av3[n+1];
        const float a33 = av3[n+2];
        const float a43 = av3[n+3];

        const float a14 = av4[n  ];
        const float a24 = av4[n+1];
        const float a34 = av4[n+2];
        const float a44 = av4[n+3];

        const float a15 = av5[n  ];
        const float a25 = av5[n+1];
        const float a35 = av5[n+2];
        const float a45 = av5[n+3];

#ifdef HAVE_OPENMP
#pragma omp simd
#endif
          for (int64_t i=0 ; i<mo_num ; ++i) {
            vgl_sp1[i] = vgl_sp1[i] + ck1[i] * a11 + ck2[i] * a21 + ck3[i] * a31 + ck4[i] * a41;
            vgl_sp2[i] = vgl_sp2[i] + ck1[i] * a12 + ck2[i] * a22 + ck3[i] * a32 + ck4[i] * a42;
            vgl_sp3[i] = vgl_sp3[i] + ck1[i] * a13 + ck2[i] * a23 + ck3[i] * a33 + ck4[i] * a43;
            vgl_sp4[i] = vgl_sp4[i] + ck1[i] * a14 + ck2[i] * a24 + ck3[i] * a34 + ck4[i] * a44;
            vgl_sp5[i] = vgl_sp5[i] + ck1[i] * a15 + ck2[i] * a25 + ck3[i] * a35 + ck4[i] * a45;
          }
      }

      for (int64_t m=n ; m < nidx ; m+=1) {
        const float* restrict ck = coefficient_t_sp + idx[m]*mo_num;
        const float a1 = av1[m];
        const float a2 = av2[m];
        const float a3 = av3[m];
        const float a4 = av4[m];
        const float a5 = av5[m];

#ifdef HAVE_OPENMP
#pragma omp simd
#endif
          for (int64_t i=0 ; i<mo_num ; ++i) {
            vgl_sp1[i] = vgl_sp1[i] + ck[i] * a1;
            vgl_sp2[i] = vgl_sp2[i] + ck[i] * a2;
            vgl_sp3[i] = vgl_sp3[i] + ck[i] * a3;
            vgl_sp4[i] = vgl_sp4[i] + ck[i] * a4;
            vgl_sp5[i] = vgl_sp5[i] + ck[i] * a5;
          }
      }
#ifdef HAVE_OPENMP
#pragma omp simd
#endif
      for (int64_t i=0 ; i<mo_num ; ++i) {
        vgl1[i] = (double) vgl_sp1[i];
        vgl2[i] = (double) vgl_sp2[i];
        vgl3[i] = (double) vgl_sp3[i];
        vgl4[i] = (double) vgl_sp4[i];
        vgl5[i] = (double) vgl_sp5[i];
      }
      free(vgl_sp1);
      free(vgl_sp2);
      free(vgl_sp3);
      free(vgl_sp4);
      free(vgl_sp5);
    }
    free(idx);
    free(av1);
    free(av2);
    free(av3);
    free(av4);
    free(av5);
  }
  free(coefficient_t_sp);
  return QMCKL_SUCCESS;
}
#endif

qmckl_exit_code
qmckl_compute_mo_basis_mo_value_cusp (const qmckl_context context,
                                      const int64_t nucl_num,
                                      const int64_t ao_num,
                                      const int64_t mo_num,
                                      const int64_t point_num,
                                      const int64_t* ao_nucl,
                                      const int32_t* ao_ang_mom,
                                      const double* en_distance,
                                      const double* r_cusp,
                                      const qmckl_tensor cusp_param_tensor,
                                      const double* coefficient_t,
                                      const double* ao_value,
                                      double* const mo_value )
{
  qmckl_exit_code rc;

#ifdef HAVE_HPC
  rc = qmckl_compute_mo_basis_mo_value_cusp_hpc (context, nucl_num, ao_num, mo_num, point_num,
                                                 ao_nucl, ao_ang_mom, en_distance, r_cusp,
                                                 cusp_param_tensor, coefficient_t, ao_value, mo_value );
#else
  double * cusp_param  = qmckl_alloc_double_of_tensor(context, cusp_param_tensor);

  rc = qmckl_compute_mo_basis_mo_value_cusp_doc (context, nucl_num, ao_num, mo_num, point_num,
                                                 ao_nucl, ao_ang_mom, en_distance, r_cusp,
                                                 cusp_param, coefficient_t, ao_value, mo_value );

  qmckl_free(context, cusp_param);
#endif
  return rc;
}

#ifdef HAVE_HPC
qmckl_exit_code
qmckl_compute_mo_basis_mo_value_cusp_hpc (const qmckl_context context,
                                          const int64_t nucl_num,
                                          const int64_t ao_num,
                                          const int64_t mo_num,
                                          const int64_t point_num,
                                          const int64_t* ao_nucl,
                                          const int32_t* ao_ang_mom,
                                          const double* en_distance,
                                          const double* r_cusp,
                                          const qmckl_tensor cusp_param,
                                          const double* coefficient_t,
                                          const double* ao_value,
                                          double* const mo_value)
{
  assert (context != QMCKL_NULL_CONTEXT);

#ifdef HAVE_OPENMP
  #pragma omp parallel for
#endif
  for (int64_t ipoint=0 ; ipoint < point_num ; ++ipoint) {
    double* restrict const vgl1  = &(mo_value[ipoint*mo_num]);
    const double* restrict avgl1 = &(ao_value[ipoint*ao_num]);
    const double* restrict ria   = &(en_distance[ipoint*nucl_num]);

    for (int64_t i=0 ; i<mo_num ; ++i) {
      vgl1[i] = 0.;
    }

    int64_t nidx=0;
    int64_t idx[ao_num];
    double  av1[ao_num];
    for (int64_t k=0 ; k<ao_num ; ++k) {
      if (avgl1[k] != 0.) {
        const int64_t inucl = ao_nucl[k];
        if (ria[inucl] > r_cusp[inucl] || ao_ang_mom[k] > 0) {
          idx[nidx] = k;
          av1[nidx] = avgl1[k];
          ++nidx;
        }
      }
    }

    int64_t n=0;

    for (n=0 ; n < nidx-4 ; n+=4) {
      const double* restrict ck1 = coefficient_t + idx[n  ]*mo_num;
      const double* restrict ck2 = coefficient_t + idx[n+1]*mo_num;
      const double* restrict ck3 = coefficient_t + idx[n+2]*mo_num;
      const double* restrict ck4 = coefficient_t + idx[n+3]*mo_num;

      const double a11 = av1[n  ];
      const double a21 = av1[n+1];
      const double a31 = av1[n+2];
      const double a41 = av1[n+3];

#ifdef HAVE_OPENMP
#pragma omp simd
#endif
      for (int64_t i=0 ; i<mo_num ; ++i) {
        vgl1[i] = vgl1[i] + ck1[i] * a11 + ck2[i] * a21 + ck3[i] * a31 + ck4[i] * a41;
      }
    }

    for (int64_t m=n ; m < nidx ; m+=1) {
      const double* restrict ck = coefficient_t + idx[m]*mo_num;
      const double a1 = av1[m];

#ifdef HAVE_OPENMP
#pragma omp simd
#endif
      for (int64_t i=0 ; i<mo_num ; ++i) {
        vgl1[i] = vgl1[i] + ck[i] * a1;
      }
    }

    for (int64_t inucl=0 ; inucl<nucl_num ; ++inucl) {
      if (ria[inucl] < r_cusp[inucl]) {
        const double r = ria[inucl];
IVDEP
        for (int64_t i=0 ; i<mo_num ; ++i) {
          vgl1[i] = vgl1[i] + qmckl_ten3(cusp_param,i,0,inucl) + r*(
                     qmckl_ten3(cusp_param,i,1,inucl) + r*(
                     qmckl_ten3(cusp_param,i,2,inucl) + r*(
                     qmckl_ten3(cusp_param,i,3,inucl) )));
        }
      }
    }

  }
  return QMCKL_SUCCESS;
}
#endif

qmckl_exit_code
qmckl_compute_mo_basis_mo_vgl_cusp (const qmckl_context context,
                                    const int64_t nucl_num,
                                    const int64_t ao_num,
                                    const int64_t mo_num,
                                    const int64_t point_num,
                                    const int64_t* ao_nucl,
                                    const int32_t* ao_ang_mom,
                                    const double* en_distance,
                                    const qmckl_matrix nucl_coord_matrix,
                                    const qmckl_matrix point_coord_matrix,
                                    const double* r_cusp,
                                    const qmckl_tensor cusp_param_tensor,
                                    const double* coefficient_t,
                                    const double* ao_vgl,
                                    double* const mo_vgl )
{
  qmckl_exit_code rc;

#ifdef HAVE_HPC
  rc = qmckl_compute_mo_basis_mo_vgl_cusp_hpc (context, nucl_num, ao_num, mo_num, point_num,
                                               ao_nucl, ao_ang_mom, en_distance, nucl_coord_matrix,
                                               point_coord_matrix, r_cusp, cusp_param_tensor,
                                               coefficient_t, ao_vgl, mo_vgl );
#else
  double * nucl_coord  = qmckl_alloc_double_of_matrix(context, nucl_coord_matrix);
  double * point_coord = qmckl_alloc_double_of_matrix(context, point_coord_matrix);
  double * cusp_param  = qmckl_alloc_double_of_tensor(context, cusp_param_tensor);

  rc = qmckl_compute_mo_basis_mo_vgl_cusp_doc (context, nucl_num, ao_num, mo_num, point_num,
                                               ao_nucl, ao_ang_mom, en_distance, nucl_coord,
                                               point_coord, r_cusp, cusp_param, coefficient_t,
                                               ao_vgl, mo_vgl );

  qmckl_free(context, nucl_coord);
  qmckl_free(context, point_coord);
  qmckl_free(context, cusp_param);
#endif
  return rc;
}

#ifdef HAVE_HPC
qmckl_exit_code
qmckl_compute_mo_basis_mo_vgl_cusp_hpc (const qmckl_context context,
                                        const int64_t nucl_num,
                                        const int64_t ao_num,
                                        const int64_t mo_num,
                                        const int64_t point_num,
                                        const int64_t* ao_nucl,
                                        const int32_t* ao_ang_mom,
                                        const double* en_distance,
                                        const qmckl_matrix nucl_coord,
                                        const qmckl_matrix point_coord,
                                        const double* r_cusp,
                                        const qmckl_tensor cusp_param,
                                        const double* coefficient_t,
                                        const double* ao_vgl,
                                        double* const mo_vgl )
{
  assert (context != QMCKL_NULL_CONTEXT);

#ifdef HAVE_OPENMP
  #pragma omp parallel for
#endif
  for (int64_t ipoint=0 ; ipoint < point_num ; ++ipoint) {
    double* restrict const vgl1 = &(mo_vgl[ipoint*5*mo_num]);
    double* restrict const vgl2 =  vgl1 + mo_num;
    double* restrict const vgl3 =  vgl1 + (mo_num << 1);
    double* restrict const vgl4 =  vgl1 + (mo_num << 1) + mo_num;
    double* restrict const vgl5 =  vgl1 + (mo_num << 2);

    const double* restrict avgl1 = &(ao_vgl[ipoint*5*ao_num]);
    const double* restrict avgl2 = avgl1 + ao_num;
    const double* restrict avgl3 = avgl1 + (ao_num << 1);
    const double* restrict avgl4 = avgl1 + (ao_num << 1) + ao_num;
    const double* restrict avgl5 = avgl1 + (ao_num << 2);

    for (int64_t i=0 ; i<mo_num ; ++i) {
      vgl1[i] = 0.;
      vgl2[i] = 0.;
      vgl3[i] = 0.;
      vgl4[i] = 0.;
      vgl5[i] = 0.;
    }

    const double* restrict ria   = &(en_distance[ipoint*nucl_num]);

    int64_t nidx=0;
    int64_t idx[ao_num];
    double  av1[ao_num];
    double  av2[ao_num];
    double  av3[ao_num];
    double  av4[ao_num];
    double  av5[ao_num];
    for (int64_t k=0 ; k<ao_num ; ++k) {
      if (avgl1[k] != 0.) {
        const int64_t inucl = ao_nucl[k];
        if (ria[inucl] > r_cusp[inucl] || ao_ang_mom[k] > 0) {
          idx[nidx] = k;
          av1[nidx] = avgl1[k];
          av2[nidx] = avgl2[k];
          av3[nidx] = avgl3[k];
          av4[nidx] = avgl4[k];
          av5[nidx] = avgl5[k];
          ++nidx;
        }
      }
    }

    int64_t n=0;

    for (n=0 ; n < nidx-4 ; n+=4) {
      const double* restrict ck1 = coefficient_t + idx[n  ]*mo_num;
      const double* restrict ck2 = coefficient_t + idx[n+1]*mo_num;
      const double* restrict ck3 = coefficient_t + idx[n+2]*mo_num;
      const double* restrict ck4 = coefficient_t + idx[n+3]*mo_num;

      const double a11 = av1[n  ];
      const double a21 = av1[n+1];
      const double a31 = av1[n+2];
      const double a41 = av1[n+3];

      const double a12 = av2[n  ];
      const double a22 = av2[n+1];
      const double a32 = av2[n+2];
      const double a42 = av2[n+3];

      const double a13 = av3[n  ];
      const double a23 = av3[n+1];
      const double a33 = av3[n+2];
      const double a43 = av3[n+3];

      const double a14 = av4[n  ];
      const double a24 = av4[n+1];
      const double a34 = av4[n+2];
      const double a44 = av4[n+3];

      const double a15 = av5[n  ];
      const double a25 = av5[n+1];
      const double a35 = av5[n+2];
      const double a45 = av5[n+3];

#ifdef HAVE_OPENMP
#pragma omp simd
#endif
      for (int64_t i=0 ; i<mo_num ; ++i) {
        vgl1[i] = vgl1[i] + ck1[i] * a11 + ck2[i] * a21 + ck3[i] * a31 + ck4[i] * a41;
        vgl2[i] = vgl2[i] + ck1[i] * a12 + ck2[i] * a22 + ck3[i] * a32 + ck4[i] * a42;
        vgl3[i] = vgl3[i] + ck1[i] * a13 + ck2[i] * a23 + ck3[i] * a33 + ck4[i] * a43;
        vgl4[i] = vgl4[i] + ck1[i] * a14 + ck2[i] * a24 + ck3[i] * a34 + ck4[i] * a44;
        vgl5[i] = vgl5[i] + ck1[i] * a15 + ck2[i] * a25 + ck3[i] * a35 + ck4[i] * a45;
      }
    }

    for (int64_t m=n ; m < nidx ; m+=1) {
      const double* restrict ck = coefficient_t + idx[m]*mo_num;
      const double a1 = av1[m];
      const double a2 = av2[m];
      const double a3 = av3[m];
      const double a4 = av4[m];
      const double a5 = av5[m];

#ifdef HAVE_OPENMP
  #pragma omp simd
#endif
      for (int64_t i=0 ; i<mo_num ; ++i) {
        vgl1[i] = vgl1[i] + ck[i] * a1;
        vgl2[i] = vgl2[i] + ck[i] * a2;
        vgl3[i] = vgl3[i] + ck[i] * a3;
        vgl4[i] = vgl4[i] + ck[i] * a4;
        vgl5[i] = vgl5[i] + ck[i] * a5;
      }
    }

    // TODO
    for (int64_t inucl=0 ; inucl<nucl_num ; ++inucl) {
      const double r = ria[inucl];
      if (r < r_cusp[inucl]) {
        const double r_vec[3] = {
          qmckl_mat(point_coord,ipoint,0) - qmckl_mat(nucl_coord,inucl,0),
          qmckl_mat(point_coord,ipoint,1) - qmckl_mat(nucl_coord,inucl,1),
          qmckl_mat(point_coord,ipoint,2) - qmckl_mat(nucl_coord,inucl,2) };
        const double r_inv = 1./r;

        for (int64_t i=0 ; i<mo_num ; ++i) {
          vgl1[i] = vgl1[i] +
                     qmckl_ten3(cusp_param,i,0,inucl) + r*(
                     qmckl_ten3(cusp_param,i,1,inucl) + r*(
                     qmckl_ten3(cusp_param,i,2,inucl) + r*(
                     qmckl_ten3(cusp_param,i,3,inucl) )));

          const double c1 = r_inv * qmckl_ten3(cusp_param,i,1,inucl) +
            2.0*qmckl_ten3(cusp_param,i,2,inucl) +
            r * 3.0 * qmckl_ten3(cusp_param,i,3,inucl);

           vgl2[i] = vgl2[i] + r_vec[0] * c1;
           vgl3[i] = vgl3[i] + r_vec[1] * c1;
           vgl4[i] = vgl4[i] + r_vec[2] * c1;

           vgl5[i] = vgl5[i] +  2.0*qmckl_ten3(cusp_param,i,1,inucl)*r_inv +
                       6.0*qmckl_ten3(cusp_param,i,2,inucl) +
                      12.0*qmckl_ten3(cusp_param,i,3,inucl)*r;
        }
      }
    }
  }
  return QMCKL_SUCCESS;
}
#endif

qmckl_exit_code
qmckl_mo_basis_rescale(qmckl_context context,
                          const double scaling_factor)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
      return qmckl_failwith( context,
                             QMCKL_NULL_CONTEXT,
                             "qmckl_mo_basis_rescale",
                             NULL);
  }

  if (scaling_factor == 0.) {
      return qmckl_failwith( context,
                             QMCKL_INVALID_ARG_2,
                             "qmckl_mo_basis_rescale",
                             "scaling factor can't be zero");
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_mo_basis_mo_vgl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->mo_basis.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_mo_basis_rescale",
                           NULL);
  }

  for (int64_t i=0 ; i<ctx->ao_basis.ao_num * ctx->mo_basis.mo_num ; ++i) {
    ctx->mo_basis.coefficient[i] *= scaling_factor;
    ctx->mo_basis.coefficient_t[i] *= scaling_factor;
  }
  rc = qmckl_context_touch(context);


  return rc;
}
