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
#include "qmckl_blas_private_type.h"

#include "qmckl_memory_private_func.h"
#include "qmckl_blas_private_func.h"



/* Allocates a new vector of the specified size. The function returns a */
/* ~qmckl_vector~ structure containing the allocated memory and size information. */
/* If the allocation fails, the returned vector will have a size of zero, */
/* which can be checked by the caller to detect allocation failures. */


qmckl_vector
qmckl_vector_alloc( qmckl_context context,
                    const int64_t size)
{
  /* Should always be true by contruction */
  assert (size > (int64_t) 0);

  qmckl_vector result;
  result.size = size;

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = size * sizeof(double);
  result.data = (double*) qmckl_malloc (context, mem_info);

  if (result.data == NULL) {
    result.size = (int64_t) 0;
  }

  return result;
}

qmckl_exit_code
qmckl_vector_free( qmckl_context context,
                   qmckl_vector* vector)
{
  if (vector == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_vector_free",
                           "Null pointer");
  }

  /* Always true */
  assert (vector->data != NULL);

  qmckl_exit_code rc;

  rc = qmckl_free(context, vector->data);
  if (rc != QMCKL_SUCCESS) {
    return rc;
  }

  vector->size = (int64_t) 0;
  vector->data = NULL;
  return QMCKL_SUCCESS;
}



/* Allocates a new matrix. If the allocation failed the sizes are zero. */


qmckl_matrix
qmckl_matrix_alloc( qmckl_context context,
                    const int64_t size1,
                    const int64_t size2)
{
  /* Should always be true by contruction */
  assert (size1 * size2 > (int64_t) 0);

  qmckl_matrix result;

  result.size[0] = size1;
  result.size[1] = size2;

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = size1 * size2 * sizeof(double);
  result.data = (double*) qmckl_malloc (context, mem_info);

  if (result.data == NULL) {
    result.size[0] = (int64_t) 0;
    result.size[1] = (int64_t) 0;
  }

  return result;
}

qmckl_exit_code
qmckl_matrix_free( qmckl_context context,
                   qmckl_matrix* matrix)
{
  if (matrix == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_matrix_free",
                           "Null pointer");
  }

  /* Always true */
  assert (matrix->data != NULL);

  qmckl_exit_code rc;

  rc = qmckl_free(context, matrix->data);
  if (rc != QMCKL_SUCCESS) {
    return rc;
  }
  matrix->data = NULL;
  matrix->size[0] = (int64_t) 0;
  matrix->size[1] = (int64_t) 0;

  return QMCKL_SUCCESS;
}



/* Allocates memory for a tensor. If the allocation failed, the size */
/* is zero. */


qmckl_tensor
qmckl_tensor_alloc( qmckl_context context,
                    const int32_t  order,
                    const int64_t* size)
{
  /* Should always be true by contruction */
  assert (order > 0);
  assert (order <= QMCKL_TENSOR_ORDER_MAX);
  assert (size  != NULL);

  qmckl_tensor result;
  memset(&result, 0, sizeof(qmckl_tensor));
  result.order = order;

  int64_t prod_size = (int64_t) 1;
  for (int32_t i=0 ; i<order ; ++i) {
    assert (size[i] > (int64_t) 0);
    result.size[i] = size[i];
    prod_size *= size[i];
  }

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = prod_size * sizeof(double);

  result.data = (double*) qmckl_malloc (context, mem_info);

  if (result.data == NULL) {
    memset(&result, 0, sizeof(qmckl_tensor));
  }

  return result;
}

qmckl_exit_code
qmckl_tensor_free( qmckl_context context,
                   qmckl_tensor* tensor)
{
  if (tensor == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_tensor_free",
                           "Null pointer");
  }

  /* Always true */
  assert (tensor->data != NULL);

  qmckl_exit_code rc;

  rc = qmckl_free(context, tensor->data);
  if (rc != QMCKL_SUCCESS) {
    return rc;
  }

  memset(tensor, 0, sizeof(qmckl_tensor));

  return QMCKL_SUCCESS;
}



/* Reshapes a vector into a matrix. */


qmckl_matrix
qmckl_matrix_of_vector(const qmckl_vector vector,
                       const int64_t size1,
                       const int64_t size2)
{
  /* Always true */
  assert (size1 * size2 == vector.size);

  qmckl_matrix result;

  result.size[0] = size1;
  result.size[1] = size2;
  result.data    = vector.data;

  return result;
}



/* Reshapes a vector into a tensor. */


qmckl_tensor
qmckl_tensor_of_vector(const qmckl_vector vector,
                       const int32_t order,
                       const int64_t* size)
{
  qmckl_tensor result;

  int64_t prod_size = 1;
  for (int32_t i=0 ; i<order ; ++i) {
    result.size[i] = size[i];
    prod_size *= size[i];
  }
  assert (prod_size == vector.size);

  result.order = order;
  result.data = vector.data;

  return result;
}



/* Reshapes a matrix into a vector. */


qmckl_vector
qmckl_vector_of_matrix(const qmckl_matrix matrix)
{
  qmckl_vector result;

  result.size = matrix.size[0] * matrix.size[1];
  result.data = matrix.data;

  return result;
}



/* Reshapes a matrix into a tensor. */


qmckl_tensor
qmckl_tensor_of_matrix(const qmckl_matrix matrix,
                       const int32_t order,
                       const int64_t* size)
{
  qmckl_tensor result;

  int64_t prod_size = 1;
  for (int32_t i=0 ; i<order ; ++i) {
    result.size[i] = size[i];
    prod_size *= size[i];
  }
  assert (prod_size == matrix.size[0] * matrix.size[1]);

  result.data = matrix.data;

  return result;
}



/* Reshapes a tensor into a vector. */


qmckl_vector
qmckl_vector_of_tensor(const qmckl_tensor tensor)
{
  int64_t prod_size = (int64_t) tensor.size[0];
  for (int32_t i=1 ; i<tensor.order ; i++) {
    prod_size *= tensor.size[i];
  }

  qmckl_vector result;

  result.size = prod_size;
  result.data = tensor.data;

  return result;
}



/* Reshapes a tensor into a vector. */


qmckl_matrix
qmckl_matrix_of_tensor(const qmckl_tensor tensor,
                       const int64_t size1,
                       const int64_t size2)
{
  /* Always true */
  int64_t prod_size = (int64_t) 1;
  for (int32_t i=0 ; i<tensor.order ; i++) {
    prod_size *= tensor.size[i];
  }
  assert (prod_size == size1 * size2);

  qmckl_matrix result;

  result.size[0] = size1;
  result.size[1] = size2;
  result.data = tensor.data;

  return result;
}

qmckl_vector
qmckl_vector_set(qmckl_vector vector, double value)
{
  for (int64_t i=0 ; i<vector.size ; ++i) {
    qmckl_vec(vector, i) = value;
  }
  return vector;
}

qmckl_matrix
qmckl_matrix_set(qmckl_matrix matrix, double value)
{
  qmckl_vector vector = qmckl_vector_of_matrix(matrix);
  for (int64_t i=0 ; i<vector.size ; ++i) {
    qmckl_vec(vector, i) = value;
  }
  return qmckl_matrix_of_vector(vector, matrix.size[0], matrix.size[1]);
}

qmckl_tensor
qmckl_tensor_set(qmckl_tensor tensor, double value)
{
  qmckl_vector vector = qmckl_vector_of_tensor(tensor);
  for (int64_t i=0 ; i<vector.size ; ++i) {
    qmckl_vec(vector, i) = value;
  }
  return qmckl_tensor_of_vector(vector, tensor.order, tensor.size);
}



/* Converts a vector to a ~double*~. */


qmckl_exit_code
qmckl_double_of_vector(const qmckl_context context,
                       const qmckl_vector vector,
                       double* const target,
                       const int64_t size_max)
{
  /* Always true by construction */
  assert (qmckl_context_check(context) != QMCKL_NULL_CONTEXT);
  assert (vector.size > (int64_t) 0);
  assert (target != NULL);
  assert (size_max > (int64_t) 0);
  assert (size_max >= vector.size);
  memcpy(target, vector.data, vector.size*sizeof(double));
  return QMCKL_SUCCESS;

}



/* Converts a matrix to a ~double*~. */


qmckl_exit_code
qmckl_double_of_matrix(const qmckl_context context,
                       const qmckl_matrix matrix,
                       double* const target,
                       const int64_t size_max)
{
  qmckl_vector vector = qmckl_vector_of_matrix(matrix);
  return qmckl_double_of_vector(context, vector, target, size_max);
}



/* Converts a tensor to a ~double*~. */


qmckl_exit_code
qmckl_double_of_tensor(const qmckl_context context,
                       const qmckl_tensor tensor,
                       double* const target,
                       const int64_t size_max)
{
  qmckl_vector vector = qmckl_vector_of_tensor(tensor);
  return qmckl_double_of_vector(context, vector, target, size_max);
}



/* Converts a ~double*~ to a vector. */


qmckl_exit_code
qmckl_vector_of_double(const qmckl_context context,
                       const double* target,
                       const int64_t size_max,
                       qmckl_vector* vector_out)
{
  qmckl_vector vector = *vector_out;
  /* Always true by construction */
  assert (qmckl_context_check(context) != QMCKL_NULL_CONTEXT);

  if (vector.size == 0) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_4,
                           "qmckl_vector_of_double",
                           "Vector not allocated");
  }

  if (vector.size != size_max) {
      return qmckl_failwith( context,
                             QMCKL_INVALID_ARG_4,
                             "qmckl_vector_of_double",
                             "Wrong vector size");
  }

  for (int64_t i=0 ; i<vector.size ; ++i) {
    vector.data[i] = target[i];
  }

  *vector_out = vector;
  return QMCKL_SUCCESS;

}



/* Converts a ~double*~ to a matrix. */


qmckl_exit_code
qmckl_matrix_of_double(const qmckl_context context,
                       const double* target,
                       const int64_t size_max,
                       qmckl_matrix* matrix)
{
  qmckl_vector vector = qmckl_vector_of_matrix(*matrix);
  qmckl_exit_code rc =
    qmckl_vector_of_double(context, target, size_max, &vector);
  *matrix = qmckl_matrix_of_vector(vector, matrix->size[0], matrix->size[1]);
  return rc;
}



/* Converts a ~double*~ to a tensor. */


qmckl_exit_code
qmckl_tensor_of_double(const qmckl_context context,
                       const double* target,
                       const int64_t size_max,
                       qmckl_tensor* tensor)
{
  qmckl_vector vector = qmckl_vector_of_tensor(*tensor);
  qmckl_exit_code rc =
    qmckl_vector_of_double(context, target, size_max, &vector);
  *tensor = qmckl_tensor_of_vector(vector, tensor->order, tensor->size);
  return rc;
}

double* qmckl_alloc_double_of_vector(const qmckl_context context,
                                     const qmckl_vector vector)
{
  /* Always true by construction */
  assert (qmckl_context_check(context) != QMCKL_NULL_CONTEXT);
  assert (vector.size > (int64_t) 0);

  qmckl_memory_info_struct mem_info = qmckl_memory_info_struct_zero;
  mem_info.size = vector.size * sizeof(double);

  double* target = (double*) qmckl_malloc(context, mem_info);
  if (target == NULL) {
     return NULL;
  }

  qmckl_exit_code rc;
  rc = qmckl_double_of_vector(context, vector, target, vector.size);
  assert (rc == QMCKL_SUCCESS);
  if (rc != QMCKL_SUCCESS) {
    rc = qmckl_free(context, target);
    target = NULL;
  }

  return target;
}

double* qmckl_alloc_double_of_matrix(const qmckl_context context,
                               const qmckl_matrix matrix)
{
  qmckl_vector vector = qmckl_vector_of_matrix(matrix);
  return qmckl_alloc_double_of_vector(context, vector);
}

double* qmckl_alloc_double_of_tensor(const qmckl_context context,
                                     const qmckl_tensor tensor)
{
  qmckl_vector vector = qmckl_vector_of_tensor(tensor);
  return qmckl_alloc_double_of_vector(context, vector);
}

qmckl_exit_code
qmckl_matmul (const qmckl_context context,
              const char TransA,
              const char TransB,
              const double alpha,
              const qmckl_matrix A,
              const qmckl_matrix B,
              const double beta,
              qmckl_matrix* const C )
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  }

  qmckl_exit_code rc = QMCKL_SUCCESS;

  if (TransA != 'N' && TransA != 'T') {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_matmul",
                           "TransA should be 'N' or 'T'");
  }

  if (TransB != 'N' && TransB != 'T') {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_matmul",
                           "TransB should be 'N' or 'T'");
  }

  if (A.size[0] < 1) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_5,
                           "qmckl_matmul",
                           "Invalid size for A");
  }

  if (B.size[0] < 1) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_6,
                           "qmckl_matmul",
                           "Invalid size for B");
  }

  if (C == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_8,
                           "qmckl_matmul",
                           "Null pointer");
  }

  int t = 0;
  if (TransA == 'T') t +=1;
  if (TransB == 'T') t +=2;
  /*
    | t | TransA | TransB |
    +---+--------+--------+
    | 0 | N      | N      |
    | 1 | T      | N      |
    | 2 | N      | T      |
    | 3 | T      | T      |
  */

  switch (t) {
  case 0:
    if (A.size[1] != B.size[0]) {
      return qmckl_failwith( context,
                             QMCKL_INVALID_ARG_2,
                             "qmckl_matmul",
                             "A and B have incompatible dimensions");
    }
    C->size[0] = A.size[0];
    C->size[1] = B.size[1];
    rc = qmckl_dgemm (context, 'N', 'N',
                      C->size[0], C->size[1], A.size[1],
                      alpha,
                      A.data, A.size[0],
                      B.data, B.size[0],
                      beta,
                      C->data, C->size[0]);
    break;
  case 1:
    if (A.size[0] != B.size[0]) {
      return qmckl_failwith( context,
                             QMCKL_INVALID_ARG_2,
                             "qmckl_matmul",
                             "A and B have incompatible dimensions");
    }
    C->size[0] = A.size[1];
    C->size[1] = B.size[1];
    rc = qmckl_dgemm (context, 'T', 'N',
                      C->size[0], C->size[1], A.size[0],
                      alpha,
                      A.data, A.size[0],
                      B.data, B.size[0],
                      beta,
                      C->data, C->size[0]);
    break;
  case 2:
    if (A.size[1] != B.size[1]) {
      return qmckl_failwith( context,
                             QMCKL_INVALID_ARG_2,
                             "qmckl_matmul",
                             "A and B have incompatible dimensions");
    }
    C->size[0] = A.size[0];
    C->size[1] = B.size[0];
    rc = qmckl_dgemm (context, 'N', 'T',
                      C->size[0], C->size[1], A.size[1],
                      alpha,
                      A.data, A.size[0],
                      B.data, B.size[0],
                      beta,
                      C->data, C->size[0]);
    break;
  case 3:
    if (A.size[0] != B.size[1]) {
      return qmckl_failwith( context,
                             QMCKL_INVALID_ARG_2,
                             "qmckl_matmul",
                             "A and B have incompatible dimensions");
    }
    C->size[0] = A.size[1];
    C->size[1] = B.size[0];
    rc = qmckl_dgemm (context, 'T', 'T',
                      C->size[0], C->size[1], A.size[0],
                      alpha,
                      A.data, A.size[0],
                      B.data, B.size[0],
                      beta,
                      C->data, C->size[0]);
    break;
  }
  return rc;
}

qmckl_exit_code
qmckl_transpose (qmckl_context context,
                 const qmckl_matrix A,
                 qmckl_matrix At )
{
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return QMCKL_INVALID_CONTEXT;
  }

  if (A.size[0] < 1) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_2,
                           "qmckl_transpose",
                           "Invalid size for A");
  }

  if (At.data == NULL) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_transpose",
                           "Output matrix not allocated");
  }

  if (At.size[0] != A.size[1] || At.size[1] != A.size[0]) {
    return qmckl_failwith( context,
                           QMCKL_INVALID_ARG_3,
                           "qmckl_transpose",
                           "Invalid size for At");
  }

  for (int64_t j=0 ; j<At.size[1] ; ++j)
    for (int64_t i=0 ; i<At.size[0] ; ++i)
      qmckl_mat(At, i, j) = qmckl_mat(A, j, i);

  return QMCKL_SUCCESS;
}

/* Convert dense (column-major) matrix A (m x n) to CSR, dropping |x| <= eps */
CSRMatrix* qmckl_dense_to_csr(const int64_t m, const int64_t n,
                              const double *A, const int64_t lda,
                              const double eps)
{
  assert(m >= 0 && n >= 0);
  int64_t *row_ptr = malloc((m + 1) * sizeof(int64_t));
  if (!row_ptr) return NULL;

  /* worst-case nnz = m * n (add +1 to avoid zero-size malloc) */
  int64_t nnz_est = (m * n) + 1;
  int64_t *col_idx = malloc(nnz_est * sizeof(int64_t));
  double  *val     = malloc(nnz_est * sizeof(double));
  if (!col_idx || !val) {
    free(row_ptr);
    free(col_idx);
    free(val);
    return NULL;
  }

  int64_t nnz = 0;
  row_ptr[0] = 0;
  for (int64_t i = 0; i < m; ++i) {
    for (int64_t j = 0; j < n; ++j) {
      double x = A[i + j * lda]; /* column-major layout */
      if (fabs(x) > eps) {
        col_idx[nnz] = j;
        val[nnz]     = x;
        ++nnz;
      }
    }
    row_ptr[i + 1] = nnz;
  }

  /*
  // shrink to actual size
  col_idx = realloc(col_idx, (nnz > 0) ? nnz * sizeof(int64_t) : 1);
  val     = realloc(val,     (nnz > 0) ? nnz * sizeof(double) : 1);
  */

  CSRMatrix *csr = malloc(sizeof(CSRMatrix));
  if (!csr) {
    free(row_ptr); free(col_idx); free(val);
    return NULL;
  }
  csr->nrows = m; csr->ncols = n; csr->nnz = nnz;
  csr->row_ptr = row_ptr; csr->col_idx = col_idx; csr->val = val;
  csr->A = A; csr->lda = lda;
  csr->density = (double) nnz / ( m * (double) n);
  return csr;
}

void free_csr(CSRMatrix *A) {
  if (!A) return;
  free(A->row_ptr);
  free(A->col_idx);
  free(A->val);
  free(A);
}

/* Convert dense (column-major) matrix A (m x n) to CSR, dropping |x| <= eps */
CSCMatrix* qmckl_dense_to_csc(const int64_t m, const int64_t n,
                              const double *A, const int64_t lda,
                              const double eps)
{
  assert(m >= 0 && n >= 0);
  int64_t *col_ptr = malloc((n + 1) * sizeof(int64_t));
  if (!col_ptr) return NULL;

  /* worst-case nnz = m * n (add +1 to avoid zero-size malloc) */
  int64_t nnz_est = (m * n) + 1;
  int64_t *row_idx = malloc(nnz_est * sizeof(int64_t));
  double  *val     = malloc(nnz_est * sizeof(double));
  if (!row_idx || !val) {
    free(col_ptr);
    free(row_idx);
    free(val);
    return NULL;
  }

  int64_t nnz = 0;
  col_ptr[0] = 0;
  for (int64_t j = 0; j < n; ++j) {
    for (int64_t i = 0; i < m; ++i) {
      double x = A[i + j * lda]; /* column-major layout */
      if (fabs(x) > eps) {
        row_idx[nnz] = i;
        val[nnz]     = x;
        ++nnz;
      }
    }
    col_ptr[j + 1] = nnz;
  }

  /*
  // shrink to actual size
  row_idx = realloc(row_idx, (nnz > 0) ? nnz * sizeof(int64_t) : 1);
  val     = realloc(val,     (nnz > 0) ? nnz * sizeof(double) : 1);
  */

  CSCMatrix *csc = malloc(sizeof(CSCMatrix));
  if (!csc) {
    free(col_ptr); free(row_idx); free(val);
    return NULL;
  }
  csc->nrows = m; csc->ncols = n; csc->nnz = nnz;
  csc->col_ptr = col_ptr; csc->row_idx = row_idx; csc->val = val;
  csc->A = A; csc->lda = lda;
  csc->density = (double) nnz / ( m * (double) n);
  return csc;
}

void free_csc(CSCMatrix *A) {
  if (!A) return;
  free(A->col_ptr);
  free(A->row_idx);
  free(A->val);
  free(A);
}

qmckl_exit_code
qmckl_dgemm_sparse_nn(const qmckl_context context,
                      const int64_t m, const int64_t n, const int64_t k,
                      const double a,
                      const CSRMatrix *Acsr, const CSRMatrix *Bcsr,
                      const double b,
                      double *C, const int64_t ldc)
{
  /* Basic context / arg checks (replace error codes/messages as desired) */
  if (qmckl_context_check(context) == 0) {
    return qmckl_failwith(context, QMCKL_INVALID_CONTEXT, "qmckl_dgemm_sparse_nn", "NULL context");
  }
  if (m <= 0) return qmckl_failwith(context, QMCKL_INVALID_ARG_2, "qmckl_dgemm_sparse_nn", "m<=0");
  if (n <= 0) return qmckl_failwith(context, QMCKL_INVALID_ARG_3, "qmckl_dgemm_sparse_nn", "n<=0");
  if (k <= 0) return qmckl_failwith(context, QMCKL_INVALID_ARG_4, "qmckl_dgemm_sparse_nn", "k<=0");
  if (Acsr == NULL) return qmckl_failwith(context, QMCKL_INVALID_ARG_6, "qmckl_dgemm_sparse_nn", "A NULL");
  if (Bcsr == NULL) return qmckl_failwith(context, QMCKL_INVALID_ARG_7, "qmckl_dgemm_sparse_nn", "B NULL");
  if (C == NULL) return qmckl_failwith(context, QMCKL_INVALID_ARG_9, "qmckl_dgemm_sparse_nn", "C NULL");

  /* Scale C by beta (column-major) */
  if (b == 0.0) {
    for (int64_t j = 0; j < n; ++j)
      for (int64_t i = 0; i < m; ++i)
        C[i + j * ldc] = 0.0;
  } else if (b != 1.0) {
    for (int64_t j = 0; j < n; ++j)
      for (int64_t i = 0; i < m; ++i)
        C[i + j * ldc] *= b;
  }

  // If either matrix is "dense enough", use tuned BLAS
  const double threshold = 0.20; // 20% nonzero
  if (Acsr->density < threshold || Bcsr->density < threshold) {

    /* Multiply: for each row i of A, iterate nonzeros (l, valA),
       then iterate nonzeros in row l of B (col j, valB) to update C[i,j]. */
    for (int64_t i = 0; i < m; ++i) {
      int64_t row_start = Acsr->row_ptr[i];
      int64_t row_end   = Acsr->row_ptr[i + 1];

      for (int64_t idxA = row_start; idxA < row_end; ++idxA) {
        const int64_t l = Acsr->col_idx[idxA];      /* column index in A (0..k-1) */
        const double valA = a * Acsr->val[idxA];

        /* traverse row l of B (B is k x n stored as CSR) */
        int64_t brow_start = Bcsr->row_ptr[l];
        int64_t brow_end   = Bcsr->row_ptr[l + 1];
        for (int64_t idxB = brow_start; idxB < brow_end; ++idxB) {
          int64_t j = Bcsr->col_idx[idxB];    /* column index in B (0..n-1) */
          double valB = Bcsr->val[idxB];

          /* update dense C (column-major) */
          C[i + j * ldc] += valA * valB;
        }
      }
    }

  } else {

    qmckl_dgemm(context, 'N', 'N', m, n, k, a, Acsr->A, Acsr->lda, Bcsr->A, Bcsr->lda, b, C, ldc);

  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_dgemm_sparse_b_nn(const qmckl_context context,
                        const int64_t m, const int64_t n, const int64_t k,
                        const double a,
                        const double * restrict A, const int64_t lda,
                        const CSCMatrix * restrict Bcsc, const double b,
                        double * restrict C, const int64_t ldc)
{
  /* Basic context / arg checks (replace error codes/messages as desired) */
  if (qmckl_context_check(context) == 0) {
    return qmckl_failwith(context, QMCKL_INVALID_CONTEXT, "qmckl_dgemm_sparse_nn", "NULL context");
  }
  if (m <= 0) return qmckl_failwith(context, QMCKL_INVALID_ARG_2, "qmckl_dgemm_sparse_nn", "m<=0");
  if (n <= 0) return qmckl_failwith(context, QMCKL_INVALID_ARG_3, "qmckl_dgemm_sparse_nn", "n<=0");
  if (k <= 0) return qmckl_failwith(context, QMCKL_INVALID_ARG_4, "qmckl_dgemm_sparse_nn", "k<=0");
  if (A == NULL) return qmckl_failwith(context, QMCKL_INVALID_ARG_6, "qmckl_dgemm_sparse_nn", "A NULL");
  if (Bcsc == NULL) return qmckl_failwith(context, QMCKL_INVALID_ARG_7, "qmckl_dgemm_sparse_nn", "B NULL");
  if (C == NULL) return qmckl_failwith(context, QMCKL_INVALID_ARG_9, "qmckl_dgemm_sparse_nn", "C NULL");

  const double threshold = 0.25; // 25% nonzero
  if (Bcsc->density < threshold) {

    /* Scale C by beta (column-major) */
    if (b == 0.0) {
      for (int64_t j = 0; j < n; ++j)
        for (int64_t i = 0; i < m; ++i)
          C[i + j * ldc] = 0.0;
    } else if (b != 1.0) {
      for (int64_t j = 0; j < n; ++j)
        for (int64_t i = 0; i < m; ++i)
          C[i + j * ldc] *= b;
    }

    for (int64_t j=0 ; j<n ; ++j) {
      const int64_t bcol_start = Bcsc->col_ptr[j];
      const int64_t bcol_end   = Bcsc->col_ptr[j+1];
      const int64_t bmax = bcol_start + ((bcol_end-bcol_start) & (~((int64_t)0x3)));
      double* restrict const C_ = &(C[j*ldc]);

      for (int64_t idxB = bcol_start ; idxB < bmax; idxB+=4) {
        const int64_t kk1 = Bcsc->row_idx[idxB+0];
        const int64_t kk2 = Bcsc->row_idx[idxB+1];
        const int64_t kk3 = Bcsc->row_idx[idxB+2];
        const int64_t kk4 = Bcsc->row_idx[idxB+3];
        const double valB1 = a*Bcsc->val[idxB+0];
        const double valB2 = a*Bcsc->val[idxB+1];
        const double valB3 = a*Bcsc->val[idxB+2];
        const double valB4 = a*Bcsc->val[idxB+3];
        const double* restrict A1 = &(A[kk1*lda]);
        const double* restrict A2 = &(A[kk2*lda]);
        const double* restrict A3 = &(A[kk3*lda]);
        const double* restrict A4 = &(A[kk4*lda]);

#ifdef HAVE_OPENMP
#pragma omp simd
#endif
        for (int64_t i=0 ; i<m ; ++i) {
          C_[i] += A1[i]*valB1 + A2[i]*valB2 + A3[i]*valB3 + A4[i]*valB4 ;
        }
      }

      for (int64_t idxB = bmax; idxB < bcol_end ; idxB+=1) {
        const int64_t kk = Bcsc->row_idx[idxB];
        const double* restrict A1 = &(A[kk*lda]);
        const double valB = a*Bcsc->val[idxB];
#ifdef HAVE_OPENMP
#pragma omp simd
#endif
        for (int64_t i=0 ; i<m ; ++i) {
          C_[i] += A1[i] * valB;
        }
      }
    }

  } else {

    qmckl_dgemm(context, 'N', 'N', m, n, k, a, A, lda, Bcsc->A, Bcsc->lda, b, C, ldc);

  }

  return QMCKL_SUCCESS;
}

int mkl_serv_intel_cpu_true() {
  return 1;
}
