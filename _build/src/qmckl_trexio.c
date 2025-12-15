#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_STDINT_H
#include <stdint.h>
#elif HAVE_INTTYPES_H
#include <inttypes.h>
#endif

#ifdef HAVE_TREXIO
#include <trexio.h>
#endif

#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>

#include <stdio.h>

#include "qmckl.h"
#include "qmckl_memory_private_type.h"
#include "qmckl_memory_private_func.h"

#ifdef HAVE_TREXIO
trexio_t* qmckl_trexio_open_X(const char* file_name, qmckl_exit_code* rc)
{
  *rc = QMCKL_SUCCESS;
  trexio_t* file = NULL;

  file = trexio_open(file_name, 'r', TREXIO_TEXT, rc);
  if (file != NULL) return file;

  file = trexio_open(file_name, 'r', TREXIO_HDF5, rc);
  if (file != NULL) return file;

  *rc = QMCKL_FAILURE;
  /* TODO
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "trexio_read_electron_up_num",
                           trexio_string_of_error(rcio));
                           */
  return NULL;
}
#endif

#ifdef HAVE_TREXIO
qmckl_exit_code
qmckl_trexio_read_electron_X(qmckl_context context, trexio_t* const file)
{
  assert (context != (qmckl_context) 0);
  assert (file != NULL);

  int rcio = 0;

  int64_t up_num = 0L;
  int64_t dn_num = 0L;

  rcio = trexio_read_electron_up_num_64(file, &up_num);
  if (rcio != TREXIO_SUCCESS) {
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "trexio_read_electron_up_num",
                           trexio_string_of_error(rcio));
  }

  assert (up_num >= 0L);

  rcio = trexio_read_electron_dn_num_64(file, &dn_num);
  if (rcio != TREXIO_SUCCESS) {
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "trexio_read_electron_dn_num",
                           trexio_string_of_error(rcio));
  }

  assert (dn_num >= 0L);


  qmckl_exit_code rc;
  rc = qmckl_set_electron_num(context, up_num, dn_num);
  return rc;
}
#endif

#ifdef HAVE_TREXIO
qmckl_exit_code
qmckl_trexio_read_nucleus_X(qmckl_context context, trexio_t* const file)
{
  assert (context != (qmckl_context) 0);
  assert (file != NULL);

  qmckl_exit_code rc;
  int rcio = 0;

int64_t nucleus_num = 0L;

rcio = trexio_read_nucleus_num_64(file, &nucleus_num);
if (rcio != TREXIO_SUCCESS) {
  return qmckl_failwith( context,
                         QMCKL_FAILURE,
                         "trexio_read_nucleus_num",
                         trexio_string_of_error(rcio));
}

assert (nucleus_num > 0);
rc = qmckl_set_nucleus_num(context, nucleus_num);

if (rc != QMCKL_SUCCESS)
  return rc;

{
  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = nucleus_num * sizeof(double);

  double* nucl_charge = (double*) qmckl_malloc(context, mem_info);

  if (nucl_charge == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_trexio_read_nucleus_X",
                           NULL);
  }

  assert (nucl_charge != NULL);

  rcio = trexio_read_safe_nucleus_charge_64(file, nucl_charge, nucleus_num);
  if (rcio != TREXIO_SUCCESS) {
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "trexio_read_nucleus_charge",
                           trexio_string_of_error(rcio));
  }

  rc = qmckl_set_nucleus_charge(context, nucl_charge, nucleus_num);

  qmckl_free(context, nucl_charge);
  nucl_charge = NULL;

  if (rc != QMCKL_SUCCESS)
    return rc;

}

{
  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = nucleus_num * 3 * sizeof(double);

  double* nucl_coord = (double*) qmckl_malloc(context, mem_info);

  if (nucl_coord == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_trexio_read_nucleus_X",
                           NULL);
  }

  assert (nucl_coord != NULL);

  rcio = trexio_read_safe_nucleus_coord_64(file, nucl_coord, 3*nucleus_num);
  if (rcio != TREXIO_SUCCESS) {
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "trexio_read_nucleus_charge",
                           trexio_string_of_error(rcio));
  }

  rc = qmckl_set_nucleus_coord(context, 'N', nucl_coord, 3*nucleus_num);

  qmckl_free(context, nucl_coord);
  nucl_coord = NULL;

  if (rc != QMCKL_SUCCESS) {
    return rc;
  }
}

assert ( qmckl_nucleus_provided(context) );
  return QMCKL_SUCCESS;
}
#endif

#ifdef HAVE_TREXIO
qmckl_exit_code
qmckl_trexio_read_ao_X(qmckl_context context, trexio_t* const file)
{
  assert (context != (qmckl_context) 0);
  assert (file != NULL);

  qmckl_exit_code rc;
  int rcio = 0;
  int64_t nucleus_num = 0L;

  rc = qmckl_get_nucleus_num(context, &nucleus_num);
  if (rc != QMCKL_SUCCESS)
    return rc;

#define MAX_STR_LEN 1024
  char basis_type[MAX_STR_LEN];

  rcio = trexio_read_basis_type(file, basis_type, MAX_STR_LEN);
  if (rcio != TREXIO_SUCCESS) {
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "trexio_read_basis_type",
                           trexio_string_of_error(rcio));
  }

  if (basis_type[0] == 'G' ||
      basis_type[0] == 'S' ||
      basis_type[0] == 'N' ) {
    rc = qmckl_set_ao_basis_type(context, basis_type[0]);
  } else {
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "trexio_read_basis_type",
                           "Invalid basis type");
  }

  if (rc != QMCKL_SUCCESS)
    return rc;

int64_t shell_num = 0L;

rcio = trexio_read_basis_shell_num_64(file, &shell_num);
if (rcio != TREXIO_SUCCESS) {
  return qmckl_failwith( context,
                         QMCKL_FAILURE,
                         "trexio_read_basis_shell_num",
                         trexio_string_of_error(rcio));
}

assert (shell_num > 0);
rc = qmckl_set_ao_basis_shell_num(context, shell_num);

if (rc != QMCKL_SUCCESS)
  return rc;

int64_t prim_num = 0L;

rcio = trexio_read_basis_prim_num_64(file, &prim_num);
if (rcio != TREXIO_SUCCESS) {
  return qmckl_failwith( context,
                         QMCKL_FAILURE,
                         "trexio_read_basis_prim_num",
                         trexio_string_of_error(rcio));
}

assert (prim_num > 0);
rc = qmckl_set_ao_basis_prim_num(context, prim_num);

if (rc != QMCKL_SUCCESS)
  return rc;

int64_t ao_num = 0;

rcio = trexio_read_ao_num_64(file, &ao_num);
if (rcio != TREXIO_SUCCESS) {
  return qmckl_failwith( context,
                         QMCKL_FAILURE,
                         "trexio_read_ao_num",
                         trexio_string_of_error(rcio));
}

assert (ao_num > 0);
rc = qmckl_set_ao_basis_ao_num(context, ao_num);

if (rc != QMCKL_SUCCESS)
  return rc;

{
  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

  /* Allocate array for data */
  mem_info.size = nucleus_num * sizeof(int64_t);
  int64_t* nucleus_index = (int64_t*) qmckl_malloc(context, mem_info);

  if (nucleus_index == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_trexio_read_basis_nucleus_index_X",
                           NULL);
  }

  assert (nucleus_index != NULL);

  /* Allocate temporary array */
  mem_info.size = shell_num * sizeof(int64_t);
  int64_t* tmp_array = (int64_t*) qmckl_malloc(context, mem_info);

  if (tmp_array == NULL) {
    qmckl_free(context, nucleus_index);
    nucleus_index = NULL;
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_trexio_read_basis_nucleus_index_X",
                           NULL);
  }

  assert (tmp_array != NULL);

  /* Read in the temporary array */
  rcio = trexio_read_safe_basis_nucleus_index_64(file, tmp_array, shell_num);
  if (rcio != TREXIO_SUCCESS) {
    qmckl_free(context, tmp_array);
    tmp_array = NULL;
    qmckl_free(context, nucleus_index);
    nucleus_index = NULL;
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "trexio_read_basis_nucleus_index",
                           trexio_string_of_error(rcio));
  }

  /* Reformat data */
  for (int i=shell_num-1 ; i>=0 ; --i) {
    const int k = tmp_array[i];
    if (k < 0 || k >= nucleus_num) {
      qmckl_free(context, tmp_array);
      tmp_array = NULL;
      qmckl_free(context, nucleus_index);
      nucleus_index = NULL;
      return qmckl_failwith( context,
                              QMCKL_FAILURE,
                              "trexio_read_basis_nucleus_index",
                              "Irrelevant data in TREXIO file");
    }
    nucleus_index[k] = i;
  }

  qmckl_free(context, tmp_array);
  tmp_array = NULL;

  /* Store data */
  rc = qmckl_set_ao_basis_nucleus_index(context, nucleus_index, shell_num);

  qmckl_free(context, nucleus_index);
  nucleus_index = NULL;

  if (rc != QMCKL_SUCCESS)
    return rc;
}

{
  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

  /* Allocate array for data */
  mem_info.size = nucleus_num * sizeof(int64_t);
  int64_t* nucleus_shell_num = (int64_t*) qmckl_malloc(context, mem_info);

  if (nucleus_shell_num == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_trexio_read_basis_nucleus_shell_num_X",
                           NULL);
  }

  assert (nucleus_shell_num != NULL);

  /* Allocate temporary array */
  mem_info.size = shell_num * sizeof(int64_t);
  int64_t* tmp_array = (int64_t*) qmckl_malloc(context, mem_info);

  if (tmp_array == NULL) {
    qmckl_free(context, nucleus_shell_num);
    nucleus_shell_num = NULL;
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_trexio_read_basis_nucleus_shell_num_X",
                           NULL);
  }

  assert (tmp_array != NULL);


  /* Read in the temporary array */
  rcio = trexio_read_safe_basis_nucleus_index_64(file, tmp_array, shell_num);
  if (rcio != TREXIO_SUCCESS) {
    qmckl_free(context, tmp_array);
    tmp_array = NULL;
    qmckl_free(context, nucleus_shell_num);
    nucleus_shell_num = NULL;
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "trexio_read_basis_nucleus_shell_num",
                           trexio_string_of_error(rcio));
  }

  /* Reformat data */
  for (int i=0 ; i<nucleus_num ; ++i) {
    nucleus_shell_num[i] = 0;
  }

  for (int i=0 ; i<shell_num ; ++i) {
    const int k = tmp_array[i];
    if (k < 0 || k >= nucleus_num) {
      qmckl_free(context, tmp_array);
      tmp_array = NULL;
      qmckl_free(context, nucleus_shell_num);
      nucleus_shell_num = NULL;
      return qmckl_failwith( context,
                             QMCKL_FAILURE,
                             "trexio_read_basis_nucleus_shell_num",
                             "Irrelevant data in TREXIO file");
    }
    nucleus_shell_num[k] += 1;
  }

  qmckl_free(context, tmp_array);
  tmp_array = NULL;

  /* Store data */
  rc = qmckl_set_ao_basis_nucleus_shell_num(context, nucleus_shell_num, shell_num);

  qmckl_free(context, nucleus_shell_num);
  nucleus_shell_num = NULL;

  if (rc != QMCKL_SUCCESS)
    return rc;
}

{
  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

  /* Allocate array for data */
  mem_info.size = shell_num * sizeof(int32_t);

  int32_t* r_power = (int32_t*) qmckl_malloc(context, mem_info);

  if (r_power == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_trexio_read_basis_r_power",
                           NULL);
  }

  assert (r_power != NULL);

  /* Read data */
  rcio = trexio_read_safe_basis_r_power_32(file, r_power, shell_num);
  if (rcio != TREXIO_SUCCESS) {
    for (int i=0 ; i<shell_num ; ++i) {
      r_power[i] = 0;
    }
  }

  /* Store data */
  rc = qmckl_set_ao_basis_r_power(context, r_power, shell_num);

  qmckl_free(context, r_power);
  r_power = NULL;

  if (rc != QMCKL_SUCCESS)
    return rc;
}

{
  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

  /* Allocate array for data */
  mem_info.size = shell_num * sizeof(int32_t);

  int32_t* shell_ang_mom = (int32_t*) qmckl_malloc(context, mem_info);

  if (shell_ang_mom == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_trexio_read_basis_shell_ang_mom_X",
                           NULL);
  }

  assert (shell_ang_mom != NULL);

  /* Read data */
  rcio = trexio_read_safe_basis_shell_ang_mom_32(file, shell_ang_mom, shell_num);
  if (rcio != TREXIO_SUCCESS) {
    qmckl_free(context, shell_ang_mom);
    shell_ang_mom = NULL;
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "trexio_read_basis_shell_ang_mom",
                           trexio_string_of_error(rcio));
  }

  /* Store data */
  rc = qmckl_set_ao_basis_shell_ang_mom(context, shell_ang_mom, shell_num);

  qmckl_free(context, shell_ang_mom);
  shell_ang_mom = NULL;

  if (rc != QMCKL_SUCCESS)
    return rc;
}

{
  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

  /* Allocate array for data */
  mem_info.size = shell_num * sizeof(int64_t);

  int64_t* shell_prim_num = (int64_t*) qmckl_malloc(context, mem_info);

  if (shell_prim_num == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_trexio_read_basis_shell_index",
                           NULL);
  }

  assert (shell_prim_num != NULL);

  /* Allocate temporary array */
  mem_info.size = prim_num * sizeof(int64_t);

  int64_t* tmp_array = (int64_t*) qmckl_malloc(context, mem_info);

  if (tmp_array == NULL) {
    qmckl_free(context, shell_prim_num);
    shell_prim_num = NULL;
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_trexio_read_basis_shell_index",
                           NULL);
  }

  assert (tmp_array != NULL);

  /* Read data */
  rcio = trexio_read_safe_basis_shell_index_64 (file, tmp_array, prim_num);
  if (rcio != TREXIO_SUCCESS) {
    qmckl_free(context, shell_prim_num);
    shell_prim_num = NULL;
    qmckl_free(context, tmp_array);
    tmp_array = NULL;
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "qmckl_trexio_read_basis_shell_index",
                           trexio_string_of_error(rcio));
  }

  /* Reformat data */
  for (int i=0 ; i<shell_num ; ++i) {
    shell_prim_num[i] = 0;
  }

  for (int i=0 ; i<prim_num ; ++i) {
    const int k = tmp_array[i];
    if (k < 0 || k >= shell_num) {
      qmckl_free(context, tmp_array);
      qmckl_free(context, shell_prim_num);
      char msg[128];
      sprintf(&msg[0], "Irrelevant data in TREXIO file: k = %d", k);
      return qmckl_failwith( context,
                              QMCKL_FAILURE,
                              "qmckl_trexio_read_basis_shell_index",
                              &msg[0]);
    }
    shell_prim_num[k] += 1;
  }

  qmckl_free(context, tmp_array);
  tmp_array = NULL;

  /* Store data */
  rc = qmckl_set_ao_basis_shell_prim_num(context, shell_prim_num, shell_num);

  qmckl_free(context, shell_prim_num);
  shell_prim_num = NULL;

  if (rc != QMCKL_SUCCESS)
    return rc;
}

{
  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

  /* Allocate array for data */
  mem_info.size = shell_num * sizeof(int64_t);

  int64_t* shell_prim_index = (int64_t*) qmckl_malloc(context, mem_info);

  if (shell_prim_index == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_trexio_read_basis_shell_prim_index_X",
                           NULL);
  }

  assert (shell_prim_index != NULL);

  /* Allocate temporary array */
  mem_info.size = prim_num * sizeof(int64_t);

  int64_t* tmp_array = (int64_t*) qmckl_malloc(context, mem_info);

  if (tmp_array == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_trexio_read_basis_shell_prim_index_X",
                           NULL);
  }

  assert (tmp_array != NULL);

  /* Read data */
  rcio = trexio_read_safe_basis_shell_index_64(file, tmp_array, prim_num);
  if (rcio != TREXIO_SUCCESS) {
    qmckl_free(context, shell_prim_index);
    shell_prim_index = NULL;
    qmckl_free(context, tmp_array);
    tmp_array = NULL;
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "trexio_read_basis_shell_prim_index",
                           trexio_string_of_error(rcio));
  }

  /* Reformat data */
  for (int i=prim_num-1 ; i>=0 ; --i) {
    const int k = tmp_array[i];
    if (k < 0 || k >= shell_num) {
      qmckl_free(context, tmp_array);
      tmp_array = NULL;
      qmckl_free(context, shell_prim_index);
      shell_prim_index = NULL;
      return qmckl_failwith( context,
                             QMCKL_FAILURE,
                             "trexio_read_basis_shell_prim_index",
                             "Irrelevant data in TREXIO file");
    }
    shell_prim_index[k] = i;
  }

  qmckl_free(context, tmp_array);
  tmp_array = NULL;

  /* Store data */
  rc = qmckl_set_ao_basis_shell_prim_index(context, shell_prim_index, shell_num);

  qmckl_free(context, shell_prim_index);
  shell_prim_index = NULL;

  if (rc != QMCKL_SUCCESS)
    return rc;
}

{
  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

  /* Allocate array for data */
  mem_info.size = shell_num * sizeof(double);

  double* shell_factor = (double*) qmckl_malloc(context, mem_info);

  if (shell_factor == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_trexio_read_basis_shell_factor_X",
                           NULL);
  }

  assert (shell_factor != NULL);

  /* Read data */
  rcio = trexio_read_safe_basis_shell_factor_64(file, shell_factor, shell_num);
  if (rcio != TREXIO_SUCCESS) {
    qmckl_free(context, shell_factor);
    shell_factor = NULL;
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "trexio_read_basis_shell_factor",
                           trexio_string_of_error(rcio));
  }

  /* Store data */
  rc = qmckl_set_ao_basis_shell_factor(context, shell_factor, shell_num);

  qmckl_free(context, shell_factor);
  shell_factor = NULL;

  if (rc != QMCKL_SUCCESS)
    return rc;
}

{
  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

  /* Allocate array for data */
  mem_info.size = prim_num * sizeof(double);

  double* exponent = (double*) qmckl_malloc(context, mem_info);

  if (exponent == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_trexio_read_basis_exponent_X",
                           NULL);
  }

  assert (exponent != NULL);

  /* Read data */
  rcio = trexio_read_safe_basis_exponent_64(file, exponent, prim_num);
  if (rcio != TREXIO_SUCCESS) {
    qmckl_free(context, exponent);
    exponent = NULL;
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "trexio_read_basis_exponent",
                           trexio_string_of_error(rcio));
  }

  /* Store data */
  rc = qmckl_set_ao_basis_exponent(context, exponent, prim_num);

  qmckl_free(context, exponent);
  exponent = NULL;

  if (rc != QMCKL_SUCCESS)
    return rc;
}

{
  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

  /* Allocate array for data */
  mem_info.size = prim_num * sizeof(double);

  double* coefficient = (double*) qmckl_malloc(context, mem_info);

  if (coefficient == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_trexio_read_basis_coefficient_X",
                           NULL);
  }

  assert (coefficient != NULL);

  /* Read data */
  rcio = trexio_read_safe_basis_coefficient_64(file, coefficient, prim_num);
  if (rcio != TREXIO_SUCCESS) {
    qmckl_free(context, coefficient);
    coefficient = NULL;
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "trexio_read_basis_coefficient",
                           trexio_string_of_error(rcio));
  }

  /* Store data */
  rc = qmckl_set_ao_basis_coefficient(context, coefficient, prim_num);

  qmckl_free(context, coefficient);
  coefficient = NULL;

  if (rc != QMCKL_SUCCESS)
    return rc;
}

{
  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

  /* Allocate array for data */
  mem_info.size = prim_num * sizeof(double);

  double* prim_factor = (double*) qmckl_malloc(context, mem_info);

  if (prim_factor == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_trexio_read_basis_prim_factor_X",
                           NULL);
  }

  assert (prim_factor != NULL);

  /* Read data */
  rcio = trexio_read_safe_basis_prim_factor_64(file, prim_factor, prim_num);
  if (rcio != TREXIO_SUCCESS) {
    qmckl_free(context, prim_factor);
    prim_factor = NULL;
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "trexio_read_basis_prim_factor",
                           trexio_string_of_error(rcio));
  }

  /* Read data */
  rc = qmckl_set_ao_basis_prim_factor(context, prim_factor, prim_num);

  qmckl_free(context, prim_factor);
  prim_factor = NULL;

  if (rc != QMCKL_SUCCESS)
    return rc;
}

{
  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;

  /* Allocate array for data */
  mem_info.size = ao_num * sizeof(double);

  double* ao_normalization = (double*) qmckl_malloc(context, mem_info);

  if (ao_normalization == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_trexio_read_ao_normalization_X",
                           NULL);
  }

  assert (ao_normalization != NULL);

  /* Read data */
  rcio = trexio_read_safe_ao_normalization_64(file, ao_normalization, ao_num);
  if (rcio != TREXIO_SUCCESS) {
    qmckl_free(context, ao_normalization);
    ao_normalization = NULL;
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "trexio_read_ao_normalization",
                           trexio_string_of_error(rcio));
  }

  /* Store data */
  rc = qmckl_set_ao_basis_ao_factor(context, ao_normalization, ao_num);

  qmckl_free(context, ao_normalization);
  ao_normalization = NULL;

  if (rc != QMCKL_SUCCESS)
    return rc;
}

return QMCKL_SUCCESS;
}
#endif

#ifdef HAVE_TREXIO
qmckl_exit_code
qmckl_trexio_read_mo_X(qmckl_context context, trexio_t* const file)
{
  assert (context != (qmckl_context) 0);
  assert (file != NULL);

  qmckl_exit_code rc;
  int rcio = 0;
  int64_t ao_num = 0L;

  rc = qmckl_get_ao_basis_ao_num(context, &ao_num);
  if (rc != QMCKL_SUCCESS)
    return rc;

int64_t mo_num = 0L;

rcio = trexio_read_mo_num_64(file, &mo_num);
if (rcio != TREXIO_SUCCESS) {
  return qmckl_failwith( context,
                         QMCKL_FAILURE,
                         "trexio_read_mo_num",
                         trexio_string_of_error(rcio));
}

assert (mo_num > 0);
rc = qmckl_set_mo_basis_mo_num(context, mo_num);

if (rc != QMCKL_SUCCESS)
  return rc;

{
  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = ao_num * mo_num * sizeof(double);

  double* mo_coef = (double*) qmckl_malloc(context, mem_info);

  if (mo_coef == NULL) {
    return qmckl_failwith( context,
                           QMCKL_ALLOCATION_FAILED,
                           "qmckl_trexio_read_mo_X",
                           NULL);
  }

  assert (mo_coef != NULL);

  rcio = trexio_read_mo_coefficient_64(file, mo_coef);
  if (rcio != TREXIO_SUCCESS) {
    return qmckl_failwith( context,
                           QMCKL_FAILURE,
                           "trexio_read_mo_coefficient",
                           trexio_string_of_error(rcio));
  }

  rc = qmckl_set_mo_basis_coefficient(context, mo_coef, ao_num*mo_num);

  qmckl_free(context, mo_coef);
  mo_coef = NULL;

  if (rc != QMCKL_SUCCESS)
    return rc;

}

return QMCKL_SUCCESS;
}
#endif

qmckl_exit_code
qmckl_trexio_read(const qmckl_context context, const char* file_name, const int64_t size_max)
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return false;
  }

  qmckl_exit_code rc;
  char file_name_new[size_max+1];
  strncpy(file_name_new, file_name, size_max);
  file_name_new[size_max] = '\0';

#ifdef HAVE_TREXIO
  trexio_t* file = qmckl_trexio_open_X(file_name_new, &rc);
  if (file == NULL) {
    trexio_close(file);
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_trexio_read",
                           trexio_string_of_error(rc));
  }

  assert (file != NULL);

  rc = qmckl_trexio_read_electron_X(context, file);
  if (rc != QMCKL_SUCCESS) {
    trexio_close(file);
    return rc;
  }

  rc = qmckl_trexio_read_nucleus_X(context, file);
  if (rc != QMCKL_SUCCESS) {
    trexio_close(file);
    return rc;
  }

  rc = qmckl_trexio_read_ao_X(context, file);
  if (rc != QMCKL_SUCCESS) {
    trexio_close(file);
    return rc;
  }

  rc = qmckl_trexio_read_mo_X(context, file);
  if (rc != QMCKL_SUCCESS) {
    trexio_close(file);
    return rc;
  }

  trexio_close(file);
  file = NULL;
#else

  rc = qmckl_failwith( context,
                       QMCKL_FAILURE,
                       "qmckl_trexio_read",
                       "QMCkl was compiled without TREXIO");
#endif
  return rc;
}
