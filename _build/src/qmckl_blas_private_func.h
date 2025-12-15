#ifndef QMCKL_BLAS_HPF
#define QMCKL_BLAS_HPF
#include "qmckl_blas_private_type.h"

/* Allocation and deallocation */


qmckl_vector
qmckl_vector_alloc( qmckl_context context,
                    const int64_t size);

qmckl_exit_code
qmckl_vector_free( qmckl_context context,
                   qmckl_vector* vector);

qmckl_matrix
qmckl_matrix_alloc( qmckl_context context,
                    const int64_t size1,
                    const int64_t size2);

qmckl_exit_code
qmckl_matrix_free( qmckl_context context,
                   qmckl_matrix* matrix);

qmckl_tensor
qmckl_tensor_alloc( qmckl_context context,
                    const int32_t order,
                    const int64_t* size);

qmckl_exit_code
qmckl_tensor_free (qmckl_context context,
                   qmckl_tensor* tensor);

/* Vector -> Matrix */


qmckl_matrix
qmckl_matrix_of_vector(const qmckl_vector vector,
                       const int64_t size1,
                       const int64_t size2);

/* Vector -> Tensor */


qmckl_tensor
qmckl_tensor_of_vector(const qmckl_vector vector,
                       const int32_t order,
                       const int64_t* size);

/* Matrix -> Vector */


qmckl_vector
qmckl_vector_of_matrix(const qmckl_matrix matrix);

/* Matrix -> Tensor */


qmckl_tensor
qmckl_tensor_of_matrix(const qmckl_matrix matrix,
                       const int32_t order,
                       const int64_t* size);

/* Tensor -> Vector */


qmckl_vector
qmckl_vector_of_tensor(const qmckl_tensor tensor);

/* Tensor -> Matrix */


qmckl_matrix
qmckl_matrix_of_tensor(const qmckl_tensor tensor,
                       const int64_t size1,
                       const int64_t size2);

#define qmckl_vec(v, i) v.data[i]
#define qmckl_mat(m, i, j) m.data[(i) + (j)*m.size[0]]

#define qmckl_ten3(t, i, j, k) t.data[(i) + t.size[0]*((j) + t.size[1]*(k))]
#define qmckl_ten4(t, i, j, k, l) t.data[(i) + t.size[0]*((j) + t.size[1]*((k) + t.size[2]*(l)))]
#define qmckl_ten5(t, i, j, k, l, m) t.data[(i) + t.size[0]*((j) + t.size[1]*((k) + t.size[2]*((l) + t.size[3]*(m))))]
#define qmckl_ten6(t, i, j, k, l, m, n) t.data[(i) + t.size[0]*((j) + t.size[1]*((k) + t.size[2]*((l) + t.size[3]*((m) + t.size[4]*(n)))))]

/* Vector */


qmckl_vector
qmckl_vector_set(qmckl_vector vector, double value);

/* Matrix */


qmckl_matrix
qmckl_matrix_set(qmckl_matrix matrix, double value);

/* Tensor */


qmckl_tensor
qmckl_tensor_set(qmckl_tensor tensor, double value);

/* Copy to/from to ~double*~ */


qmckl_exit_code
qmckl_double_of_vector(const qmckl_context context,
                       const qmckl_vector vector,
                       double* const target,
                       const int64_t size_max);

qmckl_exit_code
qmckl_double_of_matrix(const qmckl_context context,
                       const qmckl_matrix matrix,
                       double* const target,
                       const int64_t size_max);

qmckl_exit_code
qmckl_double_of_tensor(const qmckl_context context,
                       const qmckl_tensor tensor,
                       double* const target,
                       const int64_t size_max);

qmckl_exit_code
qmckl_vector_of_double(const qmckl_context context,
                       const double* target,
                       const int64_t size_max,
                       qmckl_vector* vector);

qmckl_exit_code
qmckl_matrix_of_double(const qmckl_context context,
                       const double* target,
                       const int64_t size_max,
                       qmckl_matrix* matrix);

qmckl_exit_code
qmckl_tensor_of_double(const qmckl_context context,
                       const double* target,
                       const int64_t size_max,
                       qmckl_tensor* tensor);

/* Allocate and copy to ~double*~ */


double* qmckl_alloc_double_of_vector(const qmckl_context context,
                                     const qmckl_vector vector);

double* qmckl_alloc_double_of_matrix(const qmckl_context context,
                                     const qmckl_matrix matrix);

double* qmckl_alloc_double_of_tensor(const qmckl_context context,
                                     const qmckl_tensor tensor);

/* ~qmckl_matmul~ */

/*    Matrix multiplication using the =qmckl_matrix= data type: */

/*    \[ */
/*    C_{ij} = \beta C_{ij} + \alpha \sum_{k} A_{ik} \cdot B_{kj} */
/*    \] */

/* #  TODO: Add description about the external library dependence. */

/*    #+NAME: qmckl_matmul_args */
/*    | Variable  | Type            | In/Out | Description       | */
/*    |-----------+-----------------+--------+-------------------| */
/*    | ~context~ | ~qmckl_context~ | in     | Global state      | */
/*    | ~TransA~  | ~char~          | in     | 'T' is transposed | */
/*    | ~TransB~  | ~char~          | in     | 'T' is transposed | */
/*    | ~alpha~   | ~double~        | in     | \alpha            | */
/*    | ~A~       | ~qmckl_matrix~  | in     | Matrix $A$        | */
/*    | ~B~       | ~qmckl_matrix~  | in     | Matrix $B$        | */
/*    | ~beta~    | ~double~        | in     | \beta             | */
/*    | ~C~       | ~qmckl_matrix~  | out    | Matrix $C$        | */

/*    #+CALL: generate_c_header(table=qmckl_matmul_args,rettyp="qmckl_exit_code",fname="qmckl_matmul") */

/*     #+RESULTS: */

qmckl_exit_code
qmckl_matmul (const qmckl_context context,
              const char TransA,
              const char TransB,
              const double alpha,
              const qmckl_matrix A,
              const qmckl_matrix B,
              const double beta,
              qmckl_matrix* const C );

/* ~qmckl_transpose~ */

/*    Transposes a matrix: $A^\dagger_{ji} = A_{ij}$. */

/*    | Variable  | Type            | In/Out | Description       | */
/*    |-----------+-----------------+--------+-------------------| */
/*    | ~context~ | ~qmckl_context~ | in     | Global state      | */
/*    | ~A~       | ~qmckl_matrix~  | in     | Input matrix      | */
/*    | ~At~      | ~qmckl_matrix~  | out    | Transposed matrix | */


qmckl_exit_code
qmckl_transpose (qmckl_context context,
                 const qmckl_matrix A,
                 qmckl_matrix At );

CSRMatrix* qmckl_dense_to_csr(const int64_t m, const int64_t n, const double *A, const int64_t lda,
                              const double eps);
CSCMatrix* qmckl_dense_to_csc(const int64_t m, const int64_t n, const double *A, const int64_t lda,
                              const double eps);

void free_csr(CSRMatrix* Acsr);
void free_csc(CSCMatrix* Acsc);

qmckl_exit_code qmckl_dgemm_sparse_nn (const qmckl_context context,
                                       const int64_t m,
                                       const int64_t n,
                                       const int64_t k,
                                       const double alpha,
                                       const CSRMatrix* A,
                                       const CSRMatrix* B,
                                       const double beta,
                                       double* C,
                                       const int64_t ldc);

/* ~qmckl_dgemm_sparse_b_nn~ */

/*    Sparse matrix multiplication where $A$ is dense and $B$ is in CSC (Compressed Sparse Column) format: */

/*    \[ */
/*    C_{ij} = \beta C_{ij} + \alpha \sum_{k} A_{ik} \cdot B_{kj} */
/*    \] */

/*    This function performs matrix multiplication where $A$ is a dense matrix and $B$ is stored */
/*    in CSC (Compressed Sparse Column) format. The function automatically switches between */
/*    sparse and dense BLAS implementations based on the matrix $B$ density threshold. */

/*    #+NAME: qmckl_dgemm_sparse_b_nn_args */
/*    | Variable  | Type             | In/Out | Description                                    | */
/*    |-----------+------------------+--------+------------------------------------------------| */
/*    | ~context~ | ~qmckl_context~  | in     | Global state                                   | */
/*    | ~m~       | ~int64_t~        | in     | Number of rows of matrix $A$ and $C$           | */
/*    | ~n~       | ~int64_t~        | in     | Number of columns of matrix $B$ and $C$        | */
/*    | ~k~       | ~int64_t~        | in     | Number of columns of $A$ and rows of $B$       | */
/*    | ~alpha~   | ~double~         | in     | Scalar $\alpha$                                | */
/*    | ~A~       | ~double*~        | in     | Dense matrix $A$ ($m \times k$)                | */
/*    | ~LDA~     | ~int64_t~        | in     | Leading dimension of array ~A~                 | */
/*    | ~B~       | ~CSCMatrix*~     | in     | CSC format sparse matrix $B$ ($k \times n$)    | */
/*    | ~beta~    | ~double~         | in     | Scalar $\beta$                                 | */
/*    | ~C~       | ~double*~        | out    | Dense output matrix $C$ ($m \times n$)         | */
/*    | ~ldc~     | ~int64_t~        | in     | Leading dimension of array ~C~                 | */

/*    Requirements: */

/*     - ~context~ is not ~QMCKL_NULL_CONTEXT~ */
/*     - ~m > 0~ */
/*     - ~n > 0~ */
/*     - ~k > 0~ */
/*     - ~A~ is not ~NULL~ and allocated with at least $m \times k$ elements */
/*     - ~LDA >= m~ */
/*     - ~B~ is not ~NULL~ and properly initialized as CSC matrix */
/*     - ~C~ is allocated with at least $m \times n$ elements */
/*     - ~ldc >= m~ */


qmckl_exit_code qmckl_dgemm_sparse_b_nn (const qmckl_context context,
                                         const int64_t m,
                                         const int64_t n,
                                         const int64_t k,
                                         const double alpha,
                                         const double* A,
                                         const int64_t LDA,
                                         const CSCMatrix* B,
                                         const double beta,
                                         double* C,
                                         const int64_t ldc);

#endif
