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

bool qmckl_local_energy_provided(const qmckl_context context) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return false;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  qmckl_exit_code rc;

  if(!qmckl_electron_provided(context)) return false;

  if(!qmckl_nucleus_provided(context)) return false;

  rc = qmckl_provide_ao_basis_ao_vgl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_mo_basis_mo_vgl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  ctx->local_energy.provided = (ctx->local_energy.uninitialized == 0);
  return ctx->local_energy.provided;
}

qmckl_exit_code qmckl_get_kinetic_energy(qmckl_context context, double * const kinetic_energy) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  if(!qmckl_electron_provided(context)) return QMCKL_NOT_PROVIDED;

  if(!qmckl_nucleus_provided(context)) return QMCKL_NOT_PROVIDED;

  rc = qmckl_provide_ao_basis_ao_vgl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_mo_basis_mo_vgl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_kinetic_energy(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  size_t sze = ctx->electron.walker.num * sizeof(double);
  memcpy(kinetic_energy, ctx->local_energy.e_kin, sze);

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_kinetic_energy(qmckl_context context) {

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

  if (!ctx->det.provided) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_mo_basis",
                           NULL);
  }
  rc = qmckl_provide_det_inv_matrix_alpha(context);
  if (rc != QMCKL_SUCCESS) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_det_inv_matrix_alpha",
                           NULL);
  }

  rc = qmckl_provide_det_inv_matrix_beta(context);
  if (rc != QMCKL_SUCCESS) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_det_inv_matrix_beta",
                           NULL);
  }

  /* Compute if necessary */
  if (ctx->electron.walker.point.date > ctx->local_energy.e_kin_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      free(ctx->local_energy.e_kin);
      ctx->local_energy.e_kin = NULL;
    }

    /* Allocate array */
    if (ctx->local_energy.e_kin == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.walker.num * sizeof(double);
      double* e_kin = (double*) qmckl_malloc(context, mem_info);

      if (e_kin == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_e_kin",
                               NULL);
      }
      ctx->local_energy.e_kin = e_kin;
    }

    if (ctx->det.type == 'G') {
      rc = qmckl_compute_kinetic_energy(context,
                                 ctx->electron.walker.num,
                                 ctx->det.det_num_alpha,
                                 ctx->det.det_num_beta,
                                 ctx->electron.up_num,
                                 ctx->electron.down_num,
                                 ctx->electron.num,
                                 ctx->det.mo_index_alpha,
                                 ctx->det.mo_index_beta,
                                 ctx->mo_basis.mo_num,
                                 ctx->mo_basis.mo_vgl,
                                 ctx->det.det_value_alpha,
                                 ctx->det.det_value_beta,
                                 ctx->det.det_inv_matrix_alpha,
                                 ctx->det.det_inv_matrix_beta,
                                 ctx->local_energy.e_kin);
    } else {
      return qmckl_failwith( context,
                             QMCKL_FAILURE,
                             "compute_kinetic_energy",
                             "Not yet implemented");
    }
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->local_energy.e_kin_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_get_potential_energy(qmckl_context context, double * const potential_energy) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  if(!qmckl_electron_provided(context)) return QMCKL_NOT_PROVIDED;

  if(!qmckl_nucleus_provided(context)) return QMCKL_NOT_PROVIDED;

  rc = qmckl_provide_ao_basis_ao_vgl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_mo_basis_mo_vgl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_potential_energy(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  size_t sze = ctx->electron.walker.num * sizeof(double);
  memcpy(potential_energy, ctx->local_energy.e_pot, sze);

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_potential_energy(qmckl_context context) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  qmckl_exit_code rc;
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

  rc = qmckl_provide_nucleus_repulsion(context);
  if (rc != QMCKL_SUCCESS) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_nucleus_repulsion",
                           NULL);
  }

  rc = qmckl_provide_ee_potential(context);
  if (rc != QMCKL_SUCCESS) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_ee_potential",
                           NULL);
  }

  rc = qmckl_provide_en_potential(context);
  if (rc != QMCKL_SUCCESS) {
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_en_potential",
                           NULL);
  }

  /* Compute if necessary */
  if (ctx->electron.walker.point.date > ctx->local_energy.e_pot_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      free(ctx->local_energy.e_pot);
      ctx->local_energy.e_pot = NULL;
    }

    /* Allocate array */
    if (ctx->local_energy.e_pot == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.walker.num * sizeof(double);
      double* e_pot = (double*) qmckl_malloc(context, mem_info);

      if (e_pot == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_provide_potential_energy",
                               NULL);
      }
      ctx->local_energy.e_pot = e_pot;
    }

    if (ctx->det.type == 'G') {
      rc = qmckl_compute_potential_energy(context,
                                 ctx->electron.walker.num,
                                 ctx->electron.num,
                                 ctx->nucleus.num,
                                 ctx->electron.ee_potential,
                                 ctx->electron.en_potential,
                                 ctx->nucleus.repulsion,
                                 ctx->local_energy.e_pot);
    } else {
      return qmckl_failwith( context,
                             QMCKL_FAILURE,
                             "compute_potential_energy",
                             "Not yet implemented");
    }
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->local_energy.e_pot_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_get_local_energy(qmckl_context context, double * const local_energy, const int64_t size_max) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  if(!qmckl_electron_provided(context)) return QMCKL_NOT_PROVIDED;

  if(!qmckl_nucleus_provided(context)) return QMCKL_NOT_PROVIDED;

  rc = qmckl_provide_ao_basis_ao_vgl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_mo_basis_mo_vgl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_local_energy(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  const int64_t sze = ctx->electron.walker.num;
  if (size_max < sze) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_get_local_energy",
                           "input array too small");
  }
  memcpy(local_energy, ctx->local_energy.e_local, sze * sizeof(double));

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_local_energy(qmckl_context context) {

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

  qmckl_exit_code rc;
  rc = qmckl_provide_kinetic_energy(context);
  if(rc != QMCKL_SUCCESS){
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_kinetic_energy",
                           NULL);
  }

  rc = qmckl_provide_potential_energy(context);
  if(rc != QMCKL_SUCCESS){
    return qmckl_failwith( context,
                           QMCKL_NOT_PROVIDED,
                           "qmckl_potential_energy",
                           NULL);
  }

  /* Compute if necessary */
  if (ctx->electron.walker.point.date > ctx->local_energy.e_local_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      free(ctx->local_energy.e_local);
      ctx->local_energy.e_local = NULL;
    }

    /* Allocate array */
    if (ctx->local_energy.e_local == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.walker.num * sizeof(double);
      double* local_energy = (double*) qmckl_malloc(context, mem_info);

      if (local_energy == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_local_energy",
                               NULL);
      }
      ctx->local_energy.e_local = local_energy;
    }

    if (ctx->det.type == 'G') {
      rc = qmckl_compute_local_energy(context,
                                 ctx->electron.walker.num,
                                 ctx->local_energy.e_kin,
                                 ctx->local_energy.e_pot,
                                 ctx->local_energy.e_local);
    } else {
      return qmckl_failwith( context,
                             QMCKL_FAILURE,
                             "compute_local_energy",
                             "Not yet implemented");
    }
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->local_energy.e_local_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_get_drift_vector(qmckl_context context, double * const drift_vector) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_NULL_CONTEXT;
  }

  qmckl_exit_code rc;

  if(!qmckl_electron_provided(context)) return QMCKL_NOT_PROVIDED;

  if(!qmckl_nucleus_provided(context)) return QMCKL_NOT_PROVIDED;

  rc = qmckl_provide_ao_basis_ao_vgl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_mo_basis_mo_vgl(context);
  if (rc != QMCKL_SUCCESS) return rc;

  rc = qmckl_provide_drift_vector(context);
  if (rc != QMCKL_SUCCESS) return rc;

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  size_t sze = ctx->electron.walker.num * ctx->electron.num * 3 * sizeof(double);
  memcpy(drift_vector, ctx->local_energy.r_drift, sze);

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_provide_drift_vector(qmckl_context context) {

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
  if (ctx->electron.walker.point.date > ctx->local_energy.r_drift_date) {

    if (ctx->electron.walker.num > ctx->electron.walker_old.num) {
      free(ctx->local_energy.r_drift);
      ctx->local_energy.r_drift = NULL;
    }

    /* Allocate array */
    if (ctx->local_energy.r_drift == NULL) {

      qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
      mem_info.size = ctx->electron.walker.num * ctx->electron.num * 3 * sizeof(double);
      double* r_drift = (double*) qmckl_malloc(context, mem_info);

      if (r_drift == NULL) {
        return qmckl_failwith( context,
                               QMCKL_ALLOCATION_FAILED,
                               "qmckl_r_drift",
                               NULL);
      }
      ctx->local_energy.r_drift = r_drift;
    }

    qmckl_exit_code rc;
    if (ctx->det.type == 'G') {
      rc = qmckl_compute_drift_vector(context,
                                 ctx->electron.walker.num,
                                 ctx->det.det_num_alpha,
                                 ctx->det.det_num_beta,
                                 ctx->electron.up_num,
                                 ctx->electron.down_num,
                                 ctx->electron.num,
                                 ctx->det.mo_index_alpha,
                                 ctx->det.mo_index_beta,
                                 ctx->mo_basis.mo_num,
                                 ctx->mo_basis.mo_vgl,
                                 ctx->det.det_inv_matrix_alpha,
                                 ctx->det.det_inv_matrix_beta,
                                 ctx->local_energy.r_drift);
    } else {
      return qmckl_failwith( context,
                             QMCKL_FAILURE,
                             "compute_drift_vector",
                             "Not yet implemented");
    }
    if (rc != QMCKL_SUCCESS) {
      return rc;
    }

    ctx->local_energy.r_drift_date = ctx->date;
  }

  return QMCKL_SUCCESS;
}
