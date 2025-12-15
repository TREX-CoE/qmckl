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
#include "qmckl_jastrow_champ_single_private_type.h"
#include "qmckl_jastrow_champ_single_private_func.h"
#include "qmckl_forces_private_type.h"
#include "qmckl_forces_private_func.h"



/* Performs a deep copy of the forces structure from ~src~ to ~dest~. */
/* Memory allocations are done using the provided context. */
/* Note: Most force arrays are computed arrays and will be set to NULL to allow recomputation. */


qmckl_exit_code qmckl_copy_forces(qmckl_context context, const qmckl_forces_struct* src, qmckl_forces_struct* dest) {
  
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  }

  if (src == NULL || dest == NULL) {
    return qmckl_failwith(context,
                          QMCKL_INVALID_ARG_2,
                          "qmckl_copy_forces",
                          "NULL pointer");
  }

  /* Copy all scalar fields (dates and maxsize fields) */
  dest->forces_jastrow_en_date = src->forces_jastrow_en_date;
  dest->forces_jastrow_en_g_date = src->forces_jastrow_en_g_date;
  dest->forces_jastrow_en_l_date = src->forces_jastrow_en_l_date;
  dest->forces_tmp_c_date = src->forces_tmp_c_date;
  dest->forces_dtmp_c_date = src->forces_dtmp_c_date;
  dest->forces_een_n_date = src->forces_een_n_date;
  dest->forces_jastrow_een_date = src->forces_jastrow_een_date;
  dest->forces_jastrow_een_g_date = src->forces_jastrow_een_g_date;
  dest->forces_jastrow_een_l_date = src->forces_jastrow_een_l_date;
  dest->forces_ao_value_date = src->forces_ao_value_date;
  dest->forces_ao_value_maxsize = src->forces_ao_value_maxsize;
  dest->forces_mo_value_date = src->forces_mo_value_date;
  dest->forces_mo_value_maxsize = src->forces_mo_value_maxsize;
  dest->forces_mo_g_date = src->forces_mo_g_date;
  dest->forces_mo_l_date = src->forces_mo_l_date;
  dest->forces_jastrow_single_en_date = src->forces_jastrow_single_en_date;
  dest->forces_jastrow_single_een_date = src->forces_jastrow_single_een_date;

  /* All force arrays are computed/cached values.
     For a deep copy, we set them to NULL and reset their dates to 0
     so they will be recomputed when needed.
     This is more efficient than copying potentially large cached arrays. */
  
  dest->forces_jastrow_en = NULL;
  dest->forces_jastrow_en_date = 0;
  
  dest->forces_jastrow_en_g = NULL;
  dest->forces_jastrow_en_g_date = 0;
  
  dest->forces_jastrow_en_l = NULL;
  dest->forces_jastrow_en_l_date = 0;
  
  dest->forces_tmp_c = NULL;
  dest->forces_tmp_c_date = 0;
  
  dest->forces_dtmp_c = NULL;
  dest->forces_dtmp_c_date = 0;
  
  dest->forces_een_n = NULL;
  dest->forces_een_n_date = 0;
  
  dest->forces_jastrow_een = NULL;
  dest->forces_jastrow_een_date = 0;
  
  dest->forces_jastrow_een_g = NULL;
  dest->forces_jastrow_een_g_date = 0;
  
  dest->forces_jastrow_een_l = NULL;
  dest->forces_jastrow_een_l_date = 0;
  
  dest->forces_ao_value = NULL;
  dest->forces_ao_value_date = 0;
  
  dest->forces_mo_value = NULL;
  dest->forces_mo_value_date = 0;
  
  dest->forces_mo_g = NULL;
  dest->forces_mo_g_date = 0;
  
  dest->forces_mo_l = NULL;
  dest->forces_mo_l_date = 0;
  
  dest->forces_jastrow_single_en = NULL;
  dest->forces_jastrow_single_en_date = 0;
  
  dest->forces_jastrow_single_een = NULL;
  dest->forces_jastrow_single_een_date = 0;

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_finite_difference_deriv_n(
    qmckl_context context,
    const double delta_x,
    function_callback get_function,
    double* const derivative_output,
    int64_t const size)
{
    // Finite difference coefficients for a 9-point stencil
    double coef[9] = { 1.0/280.0, -4.0/105.0, 1.0/5.0, -4.0/5.0, 0.0, 4.0/5.0, -1.0/5.0, 4.0/105.0, -1.0/280.0 };

    qmckl_exit_code rc;

    int64_t walk_num;
    rc = qmckl_get_electron_walk_num(context, &walk_num);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    int64_t nucl_num;
    rc = qmckl_get_nucleus_num(context, &nucl_num);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }


    double* nucleus_coord = (double*) malloc(3 * nucl_num * sizeof(double));
    if (nucleus_coord == NULL) {
      return QMCKL_ALLOCATION_FAILED;
    }
    rc = qmckl_get_nucleus_coord (context, 'N', nucleus_coord, 3*nucl_num);


    double* temp_coord = (double*) malloc(3 * nucl_num * sizeof(double));
    if (temp_coord == NULL) {
      free(nucleus_coord);
      return QMCKL_ALLOCATION_FAILED;
    }

    double* function_values = (double*) malloc(walk_num*size * sizeof(double));
    if (function_values == NULL) {
        free(nucleus_coord);
        free(temp_coord);
        return QMCKL_ALLOCATION_FAILED;
    }

    memset(derivative_output, 0, nucl_num*3*walk_num*size*sizeof(double));

    // Copy original coordinates
    for (int i = 0; i < 3 * nucl_num; i++) {
      temp_coord[i] = nucleus_coord[i];
    }

    for (int64_t a = 0; a < nucl_num; a++) {
      for (int64_t k = 0; k < 3; k++) {
        for (int64_t m = -4; m <= 4; m++) {

          // Apply finite difference displacement
          temp_coord[k+a*3] = nucleus_coord[k+3*a] + (double) m * delta_x;

          // Update coordinates in the context
          rc = qmckl_set_nucleus_coord(context, 'N', temp_coord, 3*nucl_num);
          assert(rc == QMCKL_SUCCESS);

          rc = qmckl_context_touch(context);
          assert(rc == QMCKL_SUCCESS);

          rc = qmckl_single_touch(context);
          assert(rc == QMCKL_SUCCESS);

          // Call the provided function
          rc = get_function(context, function_values, walk_num*size);
          assert(rc == QMCKL_SUCCESS);

          // Accumulate derivative using finite-difference coefficients
          for (int64_t nw=0 ; nw<walk_num ; nw++) {
            int64_t shift = nucl_num*3*size*nw + size*(k + 3*a);
            for (int64_t i = 0; i < size; i++) {
              derivative_output[i+shift] += coef[m + 4] * function_values[nw*size+i];
            }
          }
        }
        temp_coord[k+a*3] = nucleus_coord[k+3*a];
      }
    }

    // Reset coordinates in the context
    rc = qmckl_set_nucleus_coord(context, 'N', temp_coord, 3*nucl_num);
    assert(rc == QMCKL_SUCCESS);

    rc = qmckl_context_touch(context);
    assert(rc == QMCKL_SUCCESS);

    // Normalize by the step size
    for (int64_t i = 0; i < size*3*nucl_num*walk_num ; i++) {
      derivative_output[i] /= delta_x;
    }

    free(nucleus_coord);
    free(temp_coord);
    free(function_values);
    return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_forces_jastrow_en(qmckl_context context,
                                    double* const forces_jastrow_en,
                                    const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_forces_jastrow_en(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t sze = 3 * ctx->nucleus.num * ctx->electron.walker.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_forces_jastrow_en",
                           "Array too small. Expected 3*nucl_num*walk_num");
  }

  memcpy(forces_jastrow_en, ctx->forces.forces_jastrow_en, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_forces_jastrow_en(qmckl_context context)
{

  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_provide_forces_jastrow_en",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->jastrow_champ.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_provide_forces_jastrow_en",
                           NULL);
  }

  /* Check if en rescaled distance is provided */
  rc = qmckl_provide_en_distance_rescaled(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Check if en rescaled distance derivatives is provided */
  rc = qmckl_provide_en_distance_rescaled_gl(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->date > ctx->forces.forces_jastrow_en_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      if (ctx->forces.forces_jastrow_en != NULL) {
        rc = qmckl_free(context, ctx->forces.forces_jastrow_en);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_forces_jastrow_en",
                                 "Unable to free ctx->forces.forces_jastrow_en");
        }
        ctx->forces.forces_jastrow_en = NULL;
      }
    }

    /* Allocate array */
    if (ctx->forces.forces_jastrow_en == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.walker.num * 3 * ctx->nucleus.num * sizeof(double);
      double* forces_jastrow_en = (double*) qmckl_malloc(context, mem_info);

      if (forces_jastrow_en == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_forces_jastrow_en",
                               NULL);
      }
      ctx->forces.forces_jastrow_en = forces_jastrow_en;
    }

    rc = qmckl_compute_forces_jastrow_en(context,
                                         ctx->electron.walker.num,
                                         ctx->electron.num,
                                         ctx->nucleus.num,
                                         ctx->jastrow_champ.type_nucl_num,
                                         ctx->jastrow_champ.type_nucl_vector,
                                         ctx->jastrow_champ.aord_num,
                                         ctx->jastrow_champ.a_vector,
                                         ctx->jastrow_champ.en_distance_rescaled,
                                         ctx->jastrow_champ.en_distance_rescaled_gl,
                                         ctx->forces.forces_jastrow_en);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->forces.forces_jastrow_en_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_compute_forces_jastrow_en (const qmckl_context context,
                                          const int64_t walk_num,
                                          const int64_t elec_num,
                                          const int64_t nucl_num,
                                          const int64_t type_nucl_num,
                                          const int64_t* type_nucl_vector,
                                          const int64_t aord_num,
                                          const double* a_vector,
                                          const double* en_distance_rescaled,
                                          const double* en_distance_rescaled_gl,
                                          double* const forces_jastrow_en )
{
#ifdef HAVE_HPC
  return qmckl_compute_forces_jastrow_en_doc
#else
  return qmckl_compute_forces_jastrow_en_doc
#endif
    (context, walk_num, elec_num, nucl_num, type_nucl_num, type_nucl_vector, aord_num,
     a_vector, en_distance_rescaled, en_distance_rescaled_gl, forces_jastrow_en);
}

qmckl_exit_code
qmckl_get_forces_jastrow_en_g(qmckl_context context,
                                    double* const forces_jastrow_en_g,
                                    const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_forces_jastrow_en_g(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t sze = 3*3 * ctx->nucleus.num * ctx->electron.walker.num * ctx->electron.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_forces_jastrow_en_g",
                           "Array too small. Expected 3*3*nucl_num*walk_num_elec_num");
  }

  memcpy(forces_jastrow_en_g, ctx->forces.forces_jastrow_en_g, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_forces_jastrow_en_g(qmckl_context context)
{

  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_provide_forces_jastrow_en_g",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->jastrow_champ.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_provide_forces_jastrow_en_g",
                           NULL);
  }


  /* Check if en rescaled distance is provided */
  rc = qmckl_provide_en_distance_rescaled(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Check if en rescaled distance derivatives is provided */
  rc = qmckl_provide_en_distance_rescaled_gl(context);
  if(rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_en_distance(context);
  if(rc != QMCKL_SUCCESS) return rc;


  /* Compute if necessary */
  if (ctx->date > ctx->forces.forces_jastrow_en_g_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      if (ctx->forces.forces_jastrow_en_g != NULL) {
        rc = qmckl_free(context, ctx->forces.forces_jastrow_en_g);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_forces_jastrow_en",
                                 "Unable to free ctx->forces.forces_jastrow_en_g");
        }
        ctx->forces.forces_jastrow_en_g = NULL;
      }
    }

    /* Allocate array */
    if (ctx->forces.forces_jastrow_en_g == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.walker.num * 3 * 3 * ctx->nucleus.num * ctx->electron.num * sizeof(double);
      double* forces_jastrow_en_g = (double*) qmckl_malloc(context, mem_info);

      if (forces_jastrow_en_g == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_forces_jastrow_en_g",
                               NULL);
      }
      ctx->forces.forces_jastrow_en_g = forces_jastrow_en_g;
    }

    rc = qmckl_compute_forces_jastrow_en_g(context,
                                         ctx->electron.walker.num,
                                         ctx->electron.num,
                                         ctx->nucleus.num,
                                         ctx->jastrow_champ.type_nucl_num,
                                         ctx->jastrow_champ.type_nucl_vector,
                                         ctx->jastrow_champ.aord_num,
                                         ctx->jastrow_champ.a_vector,
                                         ctx->jastrow_champ.rescale_factor_en,
                                         ctx->electron.en_distance,
                                         ctx->jastrow_champ.en_distance_rescaled,
                                         ctx->jastrow_champ.en_distance_rescaled_gl,
                                         ctx->forces.forces_jastrow_en_g);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->forces.forces_jastrow_en_g_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_compute_forces_jastrow_en_g (const qmckl_context context,
                                          const int64_t walk_num,
                                          const int64_t elec_num,
                                          const int64_t nucl_num,
                                          const int64_t type_nucl_num,
                                          const int64_t* type_nucl_vector,
                                          const int64_t aord_num,
                                          const double* a_vector,
                                          const double* rescale_factor_en,
                                          const double* en_distance,
                                          const double* en_distance_rescaled,
                                          const double* en_distance_rescaled_gl,
                                          double* const forces_jastrow_en_g )
{
#ifdef HAVE_HPC
  return qmckl_compute_forces_jastrow_en_g_doc
#else
  return qmckl_compute_forces_jastrow_en_g_doc
#endif
    (context, walk_num, elec_num, nucl_num, type_nucl_num, type_nucl_vector, aord_num,
     a_vector, rescale_factor_en, en_distance, en_distance_rescaled, en_distance_rescaled_gl, forces_jastrow_en_g);
}

qmckl_exit_code
qmckl_get_forces_jastrow_en_l(qmckl_context context,
                                    double* const forces_jastrow_en_l,
                                    const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_forces_jastrow_en_l(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t sze = 3 * ctx->nucleus.num * ctx->electron.walker.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_forces_jastrow_en_l",
                           "Array too small. Expected 3*nucl_num*walk_num");
  }

  memcpy(forces_jastrow_en_l, ctx->forces.forces_jastrow_en_l, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_forces_jastrow_en_l(qmckl_context context)
{

  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_provide_forces_jastrow_en_l",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->jastrow_champ.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_provide_forces_jastrow_en_l",
                           NULL);
  }


  /* Check if en rescaled distance is provided */
  rc = qmckl_provide_en_distance_rescaled(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Check if en rescaled distance derivatives is provided */
  rc = qmckl_provide_en_distance_rescaled_gl(context);
  if(rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_en_distance(context);
  if(rc != QMCKL_SUCCESS) return rc;


  /* Compute if necessary */
  if (ctx->date > ctx->forces.forces_jastrow_en_l_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      if (ctx->forces.forces_jastrow_en_l != NULL) {
        rc = qmckl_free(context, ctx->forces.forces_jastrow_en_l);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_forces_jastrow_en_l",
                                 "Unable to free ctx->forces.forces_jastrow_en_l");
        }
        ctx->forces.forces_jastrow_en_l = NULL;
      }
    }

    /* Allocate array */
    if (ctx->forces.forces_jastrow_en_l == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.walker.num * 3 * ctx->nucleus.num  * sizeof(double);
      double* forces_jastrow_en_l = (double*) qmckl_malloc(context, mem_info);

      if (forces_jastrow_en_l == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_forces_jastrow_en_l",
                               NULL);
      }
      ctx->forces.forces_jastrow_en_l = forces_jastrow_en_l;
    }

    rc = qmckl_compute_forces_jastrow_en_l(context,
                                         ctx->electron.walker.num,
                                         ctx->electron.num,
                                         ctx->nucleus.num,
                                         ctx->jastrow_champ.type_nucl_num,
                                         ctx->jastrow_champ.type_nucl_vector,
                                         ctx->jastrow_champ.aord_num,
                                         ctx->jastrow_champ.a_vector,
                                         ctx->jastrow_champ.rescale_factor_en,
                                         ctx->electron.en_distance,
                                         ctx->jastrow_champ.en_distance_rescaled,
                                         ctx->jastrow_champ.en_distance_rescaled_gl,
                                         ctx->forces.forces_jastrow_en_l);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->forces.forces_jastrow_en_l_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_compute_forces_jastrow_en_l (const qmckl_context context,
                                          const int64_t walk_num,
                                          const int64_t elec_num,
                                          const int64_t nucl_num,
                                          const int64_t type_nucl_num,
                                          const int64_t* type_nucl_vector,
                                          const int64_t aord_num,
                                          const double* a_vector,
                                          const double* rescale_factor_en,
                                          const double* en_distance,
                                          const double* en_distance_rescaled,
                                          const double* en_distance_rescaled_gl,
                                          double* const forces_jastrow_en_l )
{
#ifdef HAVE_HPC
  return qmckl_compute_forces_jastrow_en_l_doc
#else
  return qmckl_compute_forces_jastrow_en_l_doc
#endif
    (context, walk_num, elec_num, nucl_num, type_nucl_num, type_nucl_vector, aord_num,
     a_vector, rescale_factor_en, en_distance, en_distance_rescaled, en_distance_rescaled_gl, forces_jastrow_en_l);
}

qmckl_exit_code
qmckl_get_forces_jastrow_single_en(qmckl_context context,
                          double* const forces_jastrow_single_en,
                          const int64_t size_max)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_forces_jastrow_single_en(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  const int64_t sze = ctx->electron.walker.num * 3 * ctx->nucleus.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_forces_jastrow_single_en",
                           "input array too small");
  }

  memcpy(forces_jastrow_single_en, ctx->forces.forces_jastrow_single_en, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_forces_jastrow_single_en(qmckl_context context)
{

  qmckl_exit_code rc = QMCKL_SUCCESS;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_provide_forces_jastrow_single_en",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->jastrow_champ.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_provide_forces_jastrow_single_en",
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
  if (ctx->single_point.date > ctx->forces.forces_jastrow_single_en_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      if (ctx->forces.forces_jastrow_single_en != NULL) {
        rc = qmckl_free(context, ctx->forces.forces_jastrow_single_en);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_forces_jastrow_single_en",
                                 "Unable to free ctx->forces.forces_jastrow_single_en");
        }
        ctx->forces.forces_jastrow_single_en = NULL;
      }
    }


    /* Allocate array */
    if (ctx->forces.forces_jastrow_single_en == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.walker.num * 3 * ctx->nucleus.num * sizeof(double);
      double* forces_jastrow_single_en = (double*) qmckl_malloc(context, mem_info);

      if (forces_jastrow_single_en == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_forces_jastrow_single_en",
                               NULL);
      }
      ctx->forces.forces_jastrow_single_en = forces_jastrow_single_en;
    }

    rc = qmckl_compute_forces_jastrow_single_en(context,
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
                                         ctx->forces.forces_jastrow_single_en);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->forces.forces_jastrow_single_en_date = ctx->single_point.date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_compute_forces_jastrow_single_en (const qmckl_context context,
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
                                          double* const forces_jastrow_single_en )
{
#ifdef HAVE_HPC
  return qmckl_compute_forces_jastrow_single_en_doc
#else
  return qmckl_compute_forces_jastrow_single_en_doc
#endif
    (context, num, walk_num, elec_num, nucl_num, type_nucl_num, type_nucl_vector, aord_num,
     a_vector, en_distance_rescaled, en_distance_rescaled_gl, en_rescaled_single, en_rescaled_single_gl, forces_jastrow_single_en );
}

qmckl_exit_code
qmckl_get_forces_tmp_c(qmckl_context context,
                                    double* const forces_tmp_c,
                                    const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_forces_tmp_c(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t sze = 4 * ctx->electron.walker.num * ctx->electron.num * ctx->nucleus.num *
       ctx->jastrow_champ.cord_num * (ctx->jastrow_champ.cord_num+1);


  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_forces_tmp_c",
                           "Array too small. Expected 4*walk_num*elec_num*nucl_num*cord_num*(cord_num+1)");
  }

  memcpy(forces_tmp_c, ctx->forces.forces_tmp_c, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_forces_tmp_c(qmckl_context context)
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

    /* Check if en rescaled distance gl derivatives is provided */
    rc = qmckl_provide_een_rescaled_n_gl(context);
    if(rc != QMCKL_SUCCESS) return rc;

  }

  /* Compute if necessary */
  if (ctx->date > ctx->forces.forces_tmp_c_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      if (ctx->forces.forces_tmp_c != NULL) {
        rc = qmckl_free(context, ctx->forces.forces_tmp_c);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_forces_tmp_c",
                                 "Unable to free ctx->forces.forces_tmp_c");
        }
        ctx->forces.forces_tmp_c = NULL;
      }
    }

    /* Allocate array */
    if (ctx->forces.forces_tmp_c == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = 4 * ctx->electron.num * ctx->jastrow_champ.cord_num *
      (ctx->jastrow_champ.cord_num+1) * ctx->nucleus.num * ctx->electron.walker.num * sizeof(double);
      double* forces_tmp_c = (double*) qmckl_malloc(context, mem_info);

      if (forces_tmp_c == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_forces_tmp_c",
                               NULL);
      }
      ctx->forces.forces_tmp_c = forces_tmp_c;
    }


    rc = qmckl_compute_forces_tmp_c(context,
                                          ctx->electron.walker.num,
                                          ctx->electron.num,
                                          ctx->nucleus.num,
                                          ctx->jastrow_champ.cord_num,
                                          ctx->jastrow_champ.een_rescaled_e,
                                          ctx->jastrow_champ.een_rescaled_n_gl,
                                          ctx->forces.forces_tmp_c);

    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->forces.forces_tmp_c_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_forces_dtmp_c(qmckl_context context,
                                    double* const forces_dtmp_c,
                                    const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_forces_dtmp_c(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t sze = 4 * 4 * ctx->electron.walker.num * ctx->electron.num * ctx->nucleus.num *
       ctx->jastrow_champ.cord_num * (ctx->jastrow_champ.cord_num+1);


  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_forces_dtmp_c",
                           "Array too small. Expected 4*4*walk_num*elec_num*nucl_num*cord_num*(cord_num+1)");
  }

  memcpy(forces_dtmp_c, ctx->forces.forces_dtmp_c, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_forces_dtmp_c(qmckl_context context)
{

  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (ctx->jastrow_champ.cord_num > 0) {

    /* Check if en rescaled distance is provided */
    rc = qmckl_provide_een_rescaled_e_gl(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if en rescaled distance gl derivatives is provided */
    rc = qmckl_provide_een_rescaled_n_gl(context);
    if(rc != QMCKL_SUCCESS) return rc;

  }

  /* Compute if necessary */
  if (ctx->date > ctx->forces.forces_dtmp_c_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      if (ctx->forces.forces_dtmp_c != NULL) {
        rc = qmckl_free(context, ctx->forces.forces_dtmp_c);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_forces_dtmp_c",
                                 "Unable to free ctx->forces.forces_dtmp_c");
        }
        ctx->forces.forces_dtmp_c = NULL;
      }
    }

    /* Allocate array */
    if (ctx->forces.forces_dtmp_c == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = 4* 4 * ctx->electron.num * ctx->jastrow_champ.cord_num *
      (ctx->jastrow_champ.cord_num+1) * ctx->nucleus.num * ctx->electron.walker.num * sizeof(double);
      double* forces_dtmp_c = (double*) qmckl_malloc(context, mem_info);

      if (forces_dtmp_c == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_forces_dtmp_c",
                               NULL);
      }
      ctx->forces.forces_dtmp_c = forces_dtmp_c;
    }


    rc = qmckl_compute_forces_dtmp_c(context,
                                          ctx->electron.walker.num,
                                          ctx->electron.num,
                                          ctx->nucleus.num,
                                          ctx->jastrow_champ.cord_num,
                                          ctx->jastrow_champ.een_rescaled_e_gl,
                                          ctx->jastrow_champ.een_rescaled_n_gl,
                                          ctx->forces.forces_dtmp_c);

    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->forces.forces_dtmp_c_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_compute_forces_dtmp_c (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t cord_num,
      const double* een_rescaled_e_gl,
      const double* een_rescaled_n_gl,
      double* const forces_dtmp_c )
{
#ifdef HAVE_HPC
  return qmckl_compute_forces_dtmp_c_hpc
#else
  return qmckl_compute_forces_dtmp_c_doc
#endif
    (context, walk_num, elec_num, nucl_num, cord_num, een_rescaled_e_gl,
     een_rescaled_n_gl, forces_dtmp_c );
}

qmckl_exit_code
qmckl_get_forces_jastrow_een(qmckl_context context,
                                    double* const forces_jastrow_een,
                                    const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_forces_jastrow_een(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t sze = 3 * ctx->nucleus.num * ctx->electron.walker.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_forces_jastrow_een",
                           "Array too small. Expected 3*nucl_num*walk_num");
  }

  memcpy(forces_jastrow_een, ctx->forces.forces_jastrow_een, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_forces_jastrow_een(qmckl_context context)
{

  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (ctx->jastrow_champ.cord_num > 0) {

    /* Check if en rescaled distance derivatives is provided */
    rc = qmckl_provide_een_rescaled_n(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if en rescaled distance gl derivatives is provided */
    rc = qmckl_provide_een_rescaled_n_gl(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if tmp_c is provided */
    rc = qmckl_provide_tmp_c(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if forces_tmp_c is provided */
    rc = qmckl_provide_forces_tmp_c(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if en rescaled distance derivatives is provided */
    rc = qmckl_provide_jastrow_champ_c_vector_full(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if en rescaled distance derivatives is provided */
    rc = qmckl_provide_lkpm_combined_index(context);
    if(rc != QMCKL_SUCCESS) return rc;



  }

  /* Compute if necessary */
  if (ctx->date > ctx->forces.forces_jastrow_een_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      if (ctx->forces.forces_jastrow_een != NULL) {
        rc = qmckl_free(context, ctx->forces.forces_jastrow_een);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_forces_jastrow_een",
                                 "Unable to free ctx->forces.forces_jastrow_een");
        }
        ctx->forces.forces_jastrow_een = NULL;
      }
    }

    /* Allocate array */
    if (ctx->forces.forces_jastrow_een == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = 3 * ctx->nucleus.num * ctx->electron.walker.num * sizeof(double);
      double* forces_jastrow_een = (double*) qmckl_malloc(context, mem_info);

      if (forces_jastrow_een == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_forces_jastrow_een",
                               NULL);
      }
      ctx->forces.forces_jastrow_een = forces_jastrow_een;
    }


    rc = qmckl_compute_forces_jastrow_een(context,
                                                ctx->electron.walker.num,
                                                ctx->electron.num,
                                                ctx->nucleus.num,
                                                ctx->jastrow_champ.cord_num,
                                                ctx->jastrow_champ.dim_c_vector,
                                                ctx->jastrow_champ.c_vector_full,
                                                ctx->jastrow_champ.lkpm_combined_index,
                                                ctx->jastrow_champ.een_rescaled_n,
                                                ctx->jastrow_champ.een_rescaled_n_gl,
                                                ctx->jastrow_champ.tmp_c,
                                                ctx->forces.forces_tmp_c,
                                                ctx->forces.forces_jastrow_een);

    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->forces.forces_jastrow_een_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_forces_een_rescaled_n_gl(qmckl_context context,
                                         double* const forces_een_n,
                                         const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_forces_een_rescaled_n_gl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t sze = 3* ctx->electron.num * 4 * ctx->nucleus.num * ctx->electron.walker.num * (ctx->jastrow_champ.cord_num + 1);
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_forces_een_rescaled_n_gl",
                           "Array too small. Expected ctx->electron.num * 3 * 4 * ctx->nucleus.num * ctx->electron.walker.num * (ctx->jastrow_champ.cord_num + 1)");
  }
  memcpy(forces_een_n, ctx->forces.forces_een_n, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_forces_een_rescaled_n_gl(qmckl_context context)
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

  /* Check if ee gl distance is provided */
  rc = qmckl_provide_een_rescaled_n_gl(context);
  if(rc != QMCKL_SUCCESS) return rc;

  /* Compute if necessary */
  if (ctx->date > ctx->forces.forces_een_n_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      if (ctx->forces.forces_een_n != NULL) {
        rc = qmckl_free(context, ctx->forces.forces_een_n);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_forces_een_rescaled_n_gl",
                                 "Unable to free ctx->forces.forces_een_n");
        }
        ctx->forces.forces_een_n = NULL;
      }
    }

    /* Allocate array */
    if (ctx->forces.forces_een_n == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.num * 3 * 4 * ctx->nucleus.num *
        ctx->electron.walker.num * (ctx->jastrow_champ.cord_num + 1) * sizeof(double);
      double* forces_een_n = (double*) qmckl_malloc(context, mem_info);

      if (forces_een_n == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_forces_een_rescaled_n_gl",
                               NULL);
      }
      ctx->forces.forces_een_n = forces_een_n;
    }

    rc = qmckl_compute_forces_een_rescaled_n_gl(context,
                                                     ctx->electron.walker.num,
                                                     ctx->electron.num,
                                                     ctx->nucleus.num,
                                                     ctx->jastrow_champ.type_nucl_num,
                                                     ctx->jastrow_champ.type_nucl_vector,
                                                     ctx->jastrow_champ.cord_num,
                                                     ctx->jastrow_champ.rescale_factor_en,
                                                     ctx->electron.en_distance,
                                                     ctx->jastrow_champ.een_rescaled_n,
                                                     ctx->jastrow_champ.een_rescaled_n_gl,
                                                     ctx->forces.forces_een_n);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->forces.forces_een_n_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_forces_jastrow_een_g(qmckl_context context,
                                     double* const forces_jastrow_een_g,
                                     const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_forces_jastrow_een_g(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t sze = ctx->electron.walker.num * 3 * ctx->electron.num * 3 * ctx->nucleus.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_forces_jastrow_een_g",
                           "Array too small. Expected 3*3*walk_num*elec_num*nucl_num");
  }
  memcpy(forces_jastrow_een_g, ctx->forces.forces_jastrow_een_g, sze*sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_forces_jastrow_een_g(qmckl_context context)
{

  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (ctx->jastrow_champ.cord_num > 0) {

    /* Check if en rescaled distance derivatives is provided */
    rc = qmckl_provide_een_rescaled_n(context);
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

    /* Check if forces_tmp_c is provided */
    rc = qmckl_provide_forces_tmp_c(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if forces_dtmp_c is provided */
    rc = qmckl_provide_forces_dtmp_c(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if ne distance is provided */
    rc = qmckl_provide_en_distance(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if en rescaled distance derivatives is provided */
    rc = qmckl_provide_forces_een_rescaled_n_gl(context);
    if(rc != QMCKL_SUCCESS) return rc;
  }

  /* Compute if necessary */
  if (ctx->date > ctx->forces.forces_jastrow_een_g_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      if (ctx->forces.forces_jastrow_een_g != NULL) {
        rc = qmckl_free(context, ctx->forces.forces_jastrow_een_g);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_forces_jastrow_een_g",
                                 "Unable to free ctx->forces.forces_jastrow_een_g");
        }
        ctx->forces.forces_jastrow_een_g = NULL;
      }
    }

    /* Allocate array */
    if (ctx->forces.forces_jastrow_een_g == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = 3 * 3 * ctx->electron.num * ctx->electron.walker.num * ctx->nucleus.num * sizeof(double);
      double* forces_jastrow_een_g = (double*) qmckl_malloc(context, mem_info);

      if (forces_jastrow_een_g == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_forces_jastrow_een_g",
                               NULL);
      }
      ctx->forces.forces_jastrow_een_g = forces_jastrow_een_g;
    }


    rc = qmckl_compute_forces_jastrow_een_g(context,
                                                ctx->electron.walker.num,
                                                ctx->electron.num,
                                                ctx->nucleus.num,
                                                ctx->jastrow_champ.cord_num,
                                                ctx->jastrow_champ.dim_c_vector,
                                                ctx->jastrow_champ.c_vector_full,
                                                ctx->jastrow_champ.lkpm_combined_index,
                                                ctx->electron.en_distance,
                                                ctx->jastrow_champ.tmp_c,
                                                ctx->jastrow_champ.dtmp_c,
                                                ctx->forces.forces_tmp_c,
                                                ctx->forces.forces_dtmp_c,
                                                ctx->forces.forces_een_n,
                                                ctx->jastrow_champ.een_rescaled_n,
                                                ctx->jastrow_champ.een_rescaled_n_gl,
                                                ctx->forces.forces_jastrow_een_g);

    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->forces.forces_jastrow_een_g_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_forces_jastrow_een_l(qmckl_context context,
                                     double* const forces_jastrow_een_l,
                                     const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_forces_jastrow_een_l(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t sze = ctx->electron.walker.num  * 3 * ctx->nucleus.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_forces_jastrow_een_l",
                           "Array too small. Expected 3*walk_num*nucl_num");
  }
  memcpy(forces_jastrow_een_l, ctx->forces.forces_jastrow_een_l, sze*sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_forces_jastrow_een_l(qmckl_context context)
{

  qmckl_exit_code rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (ctx->jastrow_champ.cord_num > 0) {

    /* Check if en rescaled distance derivatives is provided */
    rc = qmckl_provide_een_rescaled_n(context);
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

    /* Check if forces_tmp_c is provided */
    rc = qmckl_provide_forces_tmp_c(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if forces_dtmp_c is provided */
    rc = qmckl_provide_forces_dtmp_c(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if ne distance is provided */
    rc = qmckl_provide_en_distance(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if en rescaled distance derivatives is provided */
    rc = qmckl_provide_forces_een_rescaled_n_gl(context);
    if(rc != QMCKL_SUCCESS) return rc;
  }

  /* Compute if necessary */
  if (ctx->date > ctx->forces.forces_jastrow_een_l_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      if (ctx->forces.forces_jastrow_een_l != NULL) {
        rc = qmckl_free(context, ctx->forces.forces_jastrow_een_l);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_forces_jastrow_een_l",
                                 "Unable to free ctx->forces.forces_jastrow_een_l");
        }
        ctx->forces.forces_jastrow_een_l = NULL;
      }
    }

    /* Allocate array */
    if (ctx->forces.forces_jastrow_een_l == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = 3 * ctx->electron.walker.num * ctx->nucleus.num * sizeof(double);
      double* forces_jastrow_een_l = (double*) qmckl_malloc(context, mem_info);

      if (forces_jastrow_een_l == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_forces_jastrow_een_l",
                               NULL);
      }
      ctx->forces.forces_jastrow_een_l = forces_jastrow_een_l;
    }


    rc = qmckl_compute_forces_jastrow_een_l(context,
                                                ctx->electron.walker.num,
                                                ctx->electron.num,
                                                ctx->nucleus.num,
                                                ctx->jastrow_champ.cord_num,
                                                ctx->jastrow_champ.dim_c_vector,
                                                ctx->jastrow_champ.c_vector_full,
                                                ctx->jastrow_champ.lkpm_combined_index,
                                                ctx->electron.en_distance,
                                                ctx->jastrow_champ.tmp_c,
                                                ctx->jastrow_champ.dtmp_c,
                                                ctx->forces.forces_tmp_c,
                                                ctx->forces.forces_dtmp_c,
                                                ctx->forces.forces_een_n,
                                                ctx->jastrow_champ.een_rescaled_n,
                                                ctx->jastrow_champ.een_rescaled_n_gl,
                                                ctx->forces.forces_jastrow_een_l);

    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->forces.forces_jastrow_een_l_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_forces_jastrow_single_een(qmckl_context context,
                          double* const forces_jastrow_single_een,
                          const int64_t size_max)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_forces_jastrow_single_een(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  const int64_t sze = ctx->electron.walker.num * 3 * ctx->nucleus.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_forces_jastrow_single_een",
                           "input array too small");
  }

  memcpy(forces_jastrow_single_een, ctx->forces.forces_jastrow_single_een, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_forces_jastrow_single_een(qmckl_context context)
{

  qmckl_exit_code rc = QMCKL_SUCCESS;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_provide_forces_jastrow_single_een",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->jastrow_champ.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_provide_forces_jastrow_single_een",
                           NULL);
  }

  if (ctx->jastrow_champ.cord_num > 0) {

    /* Check if en rescaled distance is provided */
    rc = qmckl_provide_een_rescaled_n(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if en rescaled distance derivatives is provided */
    rc = qmckl_provide_een_rescaled_n_gl(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if en rescaled distance is provided */
    rc = qmckl_provide_een_rescaled_single_n(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if en rescaled distance derivatives is provided */
    rc = qmckl_provide_een_rescaled_single_n_gl(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if ee rescaled distance is provided */
    rc = qmckl_provide_een_rescaled_e(context);
    if(rc != QMCKL_SUCCESS) return rc;

    /* Check if single eerescaled distance is provided */
    rc = qmckl_provide_een_rescaled_single_e(context);
    if(rc != QMCKL_SUCCESS) return rc;

      /* Check if single tmp matrix is provided */
    rc = qmckl_provide_jastrow_champ_single_tmp(context);
    if(rc != QMCKL_SUCCESS) return rc;

  }


  /* Compute if necessary */
  if (ctx->single_point.date > ctx->forces.forces_jastrow_single_een_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      if (ctx->forces.forces_jastrow_single_een != NULL) {
        rc = qmckl_free(context, ctx->forces.forces_jastrow_single_een);
        if (rc != QMCKL_SUCCESS) {
          return qmckl_failwith( context, rc,
                                 "qmckl_provide_forces_jastrow_single_een",
                                 "Unable to free ctx->forces.forces_jastrow_single_een");
        }
        ctx->forces.forces_jastrow_single_een = NULL;
      }
    }


    /* Allocate array */
    if (ctx->forces.forces_jastrow_single_een == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.walker.num * 3 * ctx->nucleus.num * sizeof(double);
      double* forces_jastrow_single_een = (double*) qmckl_malloc(context, mem_info);

      if (forces_jastrow_single_een == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_forces_jastrow_single_een",
                               NULL);
      }
      ctx->forces.forces_jastrow_single_een = forces_jastrow_single_een;
    }

    rc = qmckl_compute_forces_jastrow_single_een(context,
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
                                         ctx->single_point.tmp,
                                         ctx->forces.forces_jastrow_single_een);
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->forces.forces_jastrow_single_een_date = ctx->single_point.date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_compute_forces_jastrow_single_een (const qmckl_context context,
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
                                          const double* tmp,
                                          double* const forces_jastrow_single_een )
{
#ifdef HAVE_HPC
  return qmckl_compute_forces_jastrow_single_een_doc
#else
  return qmckl_compute_forces_jastrow_single_een_doc
#endif
    (context, num, walk_num, elec_num, nucl_num, cord_num, dim_c_vector, c_vector_full, lkpm_combined_index,
     een_rescaled_n, een_rescaled_single_n,  een_rescaled_n_gl, een_rescaled_single_n_gl, 
     een_rescaled_e, een_rescaled_single_e, tmp, forces_jastrow_single_een );
}

qmckl_exit_code
qmckl_get_forces_ao_value(qmckl_context context,
                                     double* const forces_ao_value,
                                     const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_forces_ao_value(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  int64_t sze = ctx->ao_basis.ao_num * ctx->nucleus.num * 3 * ctx->point.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_forces_ao_value",
                           "Array too small. Expected walk_num*nucl_num*point_num*3");
  }
  memcpy(forces_ao_value, ctx->forces.forces_ao_value, sze*sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_forces_ao_value(qmckl_context context)
{
  qmckl_exit_code rc = QMCKL_SUCCESS;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_provide_forces_ao_value",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->ao_basis.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_provide_forces_ao_value",
                           NULL);
  }

  /* Compute if necessary */
  if (ctx->point.date > ctx->forces.forces_ao_value_date) {

    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = ctx->ao_basis.ao_num * 3 * ctx->nucleus.num * ctx->point.num * sizeof(double);

    if (mem_info.size > ctx->forces.forces_ao_value_maxsize) {
      if (ctx->forces.forces_ao_value != NULL) {
          rc = qmckl_free(context, ctx->forces.forces_ao_value);
          if (rc != QMCKL_SUCCESS) {
            return qmckl_failwith( context, rc,
                                  "qmckl_provide_forces_ao_value",
                                  "Unable to free ctx->forces.forces_ao_value");
          }
          ctx->forces.forces_ao_value = NULL;

      }
    }
    /* Allocate array */
    if (ctx->forces.forces_ao_value == NULL) {

      double* forces_ao_value = (double*) qmckl_malloc(context, mem_info);
      ctx->forces.forces_ao_value_maxsize = mem_info.size;

      if (forces_ao_value == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_forces_ao_value",
                               NULL);
      }
      ctx->forces.forces_ao_value = forces_ao_value;
    }

    rc = qmckl_provide_ao_basis_ao_vgl(context);
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context, rc, "qmckl_provide_ao_basis_ao_vgl", NULL);
    }

    //memset(ctx->forces.forces_ao_value, 0, mem_info.size);

    rc = qmckl_compute_forces_ao_value_doc(context,
                                  ctx->ao_basis.ao_num,
                                  ctx->ao_basis.shell_num,
                                  ctx->point.num,
                                  ctx->nucleus.num,
                                  ctx->ao_basis.nucleus_index,
                                  ctx->ao_basis.nucleus_shell_num,
                                  ctx->ao_basis.shell_ang_mom,
                                  ctx->ao_basis.ao_factor,
                                  ctx->ao_basis.ao_vgl,
                                  ctx->forces.forces_ao_value);

    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->forces.forces_ao_value_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_forces_mo_value(qmckl_context context,
                          double* const forces_mo_value,
                          const int64_t size_max)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_forces_mo_value(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  const int64_t sze = ctx->point.num * 3 * ctx->mo_basis.mo_num * ctx->nucleus.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_forces_mo_value",
                           "input array too small");
  }
  memcpy(forces_mo_value, ctx->forces.forces_mo_value, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_forces_mo_value_inplace (qmckl_context context,
                                     double* const forces_mo_value,
                                     const int64_t size_max)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_forces_mo_value",
                           NULL);
  }

  qmckl_exit_code rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  const int64_t sze = ctx->point.num * 3 * ctx->mo_basis.mo_num * ctx->nucleus.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_mo_basis_mo_value",
                           "input array too small");
  }

  ctx->forces.forces_mo_value_date = ctx->point.date - 1UL;
  double* old_array = ctx->forces.forces_mo_value;

  ctx->forces.forces_mo_value = forces_mo_value;

  rc = qmckl_provide_forces_mo_value(context);
  if (rc != QMCKL_SUCCESS) return rc;

  ctx->forces.forces_mo_value = old_array;

  ctx->forces.forces_mo_value_date = ctx->point.date - 1UL;

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_forces_mo_value(qmckl_context context)
{

  qmckl_exit_code rc = QMCKL_SUCCESS;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_provide_forces_mo_value",
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
  if (ctx->point.date > ctx->forces.forces_mo_value_date) {

    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = 3 * ctx->mo_basis.mo_num * ctx->point.num * ctx->nucleus.num * sizeof(double);
    
    if (ctx->forces.forces_mo_value != NULL) {
      qmckl_memory_info_struct mem_info_test = qmckl_memory_info_struct_zero;
      rc = qmckl_get_malloc_info(context, ctx->forces.forces_mo_value, &mem_info_test);

      /* if rc != QMCKL_SUCCESS, we are maybe in an _inplace function because the
         memory was not allocated with qmckl_malloc */

      if ((rc == QMCKL_SUCCESS) && (mem_info_test.size != mem_info.size)) {
        rc = qmckl_free(context, ctx->forces.forces_mo_value);
        assert (rc == QMCKL_SUCCESS);
        ctx->forces.forces_mo_value = NULL;
      }
    }

    /* Allocate array */
    if (ctx->forces.forces_mo_value == NULL) {

      double* forces_mo_value = (double*) qmckl_malloc(context, mem_info);
      ctx->forces.forces_mo_value_maxsize = mem_info.size;

      if (forces_mo_value == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_forces_mo_value",
                               NULL);
      }
      ctx->forces.forces_mo_value = forces_mo_value;
    }

    rc = qmckl_provide_ao_basis_ao_vgl(context);
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context,
                             QMCKL_NOT_PROVIDED,
                             "qmckl_forces_ao_value",
                             NULL);
    }

    rc = qmckl_compute_forces_mo_value_doc(context,
                                    ctx->nucleus.num,
                                    ctx->ao_basis.ao_num,
                                    ctx->mo_basis.mo_num,
                                    ctx->point.num,
                                    ctx->ao_basis.shell_num,
                                    ctx->ao_basis.nucleus_index,
                                    ctx->ao_basis.nucleus_shell_num,
                                    ctx->ao_basis.shell_ang_mom,
                                    ctx->mo_basis.coefficient_t,
                                    ctx->ao_basis.ao_vgl,
                                    ctx->forces.forces_mo_value);

    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->forces.forces_mo_value_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_forces_mo_g(qmckl_context context,
                          double* const forces_mo_g,
                          const int64_t size_max)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_forces_mo_g",
                           NULL);
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_forces_mo_g(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  const int64_t sze = ctx->point.num * 3 * 3 * ctx->mo_basis.mo_num * ctx->nucleus.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_forces_mo_g",
                           "input array too small");
  }

  memcpy(forces_mo_g, ctx->forces.forces_mo_g, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_forces_mo_g_inplace (qmckl_context context,
                                     double* const forces_mo_g,
                                     const int64_t size_max)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_get_forces_mo_g",
                           NULL);
  }

  qmckl_exit_code rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  const int64_t sze = ctx->point.num * 3 * 3 * ctx->mo_basis.mo_num * ctx->nucleus.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_forces_mo_g",
                           "input array too small");
  }

  ctx->forces.forces_mo_g_date = ctx->point.date - 1UL;

  double* old_array =  ctx->forces.forces_mo_g;

  ctx->forces.forces_mo_g = forces_mo_g;

  rc = qmckl_provide_forces_mo_g(context);
  if (rc != QMCKL_SUCCESS) return rc;

  ctx->forces.forces_mo_g = old_array;

  ctx->forces.forces_mo_g_date = ctx->point.date - 1UL;

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_forces_mo_g(qmckl_context context)
{

  qmckl_exit_code rc = QMCKL_SUCCESS;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_provide_forces_mo_g",
                           NULL);
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  if (!ctx->mo_basis.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_provide_forces_mo_g",
                           NULL);
  }

  /* Compute if necessary */
  if (ctx->point.date > ctx->forces.forces_mo_g_date) {

    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = 3 * 3 * ctx->mo_basis.mo_num * ctx->point.num * ctx->nucleus.num * sizeof(double);

    if (ctx->forces.forces_mo_g != NULL) {
      qmckl_memory_info_struct mem_info_test = qmckl_memory_info_struct_zero;
      rc = qmckl_get_malloc_info(context, ctx->forces.forces_mo_g, &mem_info_test);

      /* if rc != QMCKL_SUCCESS, we are maybe in an _inplace function because the
         memory was not allocated with qmckl_malloc */

      if ((rc == QMCKL_SUCCESS) && (mem_info_test.size != mem_info.size)) {
        rc = qmckl_free(context, ctx->forces.forces_mo_g);
        assert (rc == QMCKL_SUCCESS);
        ctx->forces.forces_mo_g = NULL;
      }
    }

    /* Allocate array */
    if (ctx->forces.forces_mo_g == NULL) {

      double* forces_mo_g = (double*) qmckl_malloc(context, mem_info);

      if (forces_mo_g == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_forces_mo_g",
                               NULL);
      }
      ctx->forces.forces_mo_g = forces_mo_g;
    }

    rc = qmckl_provide_ao_basis_ao_hessian(context);
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context,
                             QMCKL_NOT_PROVIDED,
                             "qmckl_ao_basis_ao_hessian",
                             NULL);
    }

    rc = qmckl_compute_forces_mo_g(context,
                                   ctx->ao_basis.ao_num,
                                   ctx->mo_basis.mo_num,
                                   ctx->point.num,
                                   ctx->nucleus.num,
                                   ctx->ao_basis.shell_num,
                                   ctx->ao_basis.nucleus_index,
                                   ctx->ao_basis.nucleus_shell_num,
                                   ctx->ao_basis.shell_ang_mom,
                                   ctx->mo_basis.coefficient_t,
                                   ctx->ao_basis.ao_hessian,
                                   ctx->forces.forces_mo_g);

    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->forces.forces_mo_g_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_compute_forces_mo_g (
         const qmckl_context context,
         const int64_t ao_num,
         const int64_t mo_num,
         const int64_t point_num,
         const int64_t nucl_num,
         const int64_t shell_num,
         const int64_t* nucleus_index,
         const int64_t* nucleus_shell_num,
         const int32_t* shell_ang_mom,
         const double* coefficient_t,
         const double* ao_hessian,
         double* const forces_mo_g )
{
#ifdef HAVE_HPC
  return qmckl_compute_forces_mo_g_hpc
#else
  return qmckl_compute_forces_mo_g_hpc
#endif
   (context, ao_num, mo_num, point_num, nucl_num, shell_num, nucleus_index,
   nucleus_shell_num, shell_ang_mom, coefficient_t, ao_hessian, forces_mo_g );
}

qmckl_exit_code
qmckl_get_forces_mo_l(qmckl_context context,
                          double* const forces_mo_l,
                          const int64_t size_max)
{

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  rc = qmckl_provide_forces_mo_l(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  const int64_t sze = ctx->point.num * 3 * ctx->mo_basis.mo_num * ctx->nucleus.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_forces_mo_l",
                           "input array too small");
  }

  memcpy(forces_mo_l, ctx->forces.forces_mo_l, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_forces_mo_l(qmckl_context context)
{

  qmckl_exit_code rc = QMCKL_SUCCESS;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_CONTEXT,
                           "qmckl_provide_forces_mo_l",
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
  if (ctx->point.date > ctx->forces.forces_mo_l_date) {

    qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
    mem_info.size = 3  * ctx->mo_basis.mo_num * ctx->point.num * ctx->nucleus.num * sizeof(double);

    if (ctx->forces.forces_mo_l != NULL) {
      qmckl_memory_info_struct mem_info_test = qmckl_memory_info_struct_zero;
      rc = qmckl_get_malloc_info(context, ctx->forces.forces_mo_l, &mem_info_test);

      /* if rc != QMCKL_SUCCESS, we are maybe in an _inplace function because the
         memory was not allocated with qmckl_malloc */

      if ((rc == QMCKL_SUCCESS) && (mem_info_test.size != mem_info.size)) {
        rc = qmckl_free(context, ctx->forces.forces_mo_l);
        assert (rc == QMCKL_SUCCESS);
        ctx->forces.forces_mo_l= NULL;
      }
    }

    /* Allocate array */
    if (ctx->forces.forces_mo_l == NULL) {

      double* forces_mo_l = (double*) qmckl_malloc(context, mem_info);

      if (forces_mo_l == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_forces_mo_l",
                               NULL);
      }
      ctx->forces.forces_mo_l = forces_mo_l;
    }

    rc = qmckl_provide_ao_basis_ao_hessian(context);
    if (rc != QMCKL_SUCCESS) {
      return qmckl_failwith( context,
                             QMCKL_NOT_PROVIDED,
                             "qmckl_ao_basis_ao_hessian",
                             NULL);
    }


    rc = qmckl_compute_forces_mo_l_doc(context,
                                         ctx->ao_basis.ao_num,
                                         ctx->mo_basis.mo_num,
                                         ctx->point.num,
                                         ctx->nucleus.num,
                                         ctx->ao_basis.shell_num,
                                         ctx->ao_basis.nucleus_index,
                                         ctx->ao_basis.nucleus_shell_num,
                                         ctx->ao_basis.shell_ang_mom,
                                         ctx->mo_basis.coefficient_t,
                                         ctx->ao_basis.ao_hessian,
                                         ctx->forces.forces_mo_l);

    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->forces.forces_mo_l_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}
