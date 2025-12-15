#ifndef QMCKL_BLAS_HPT
#define QMCKL_BLAS_HPT

/* Vector */

/*   A vector represents a one-dimensional array of double-precision floating-point */
/*   numbers. The structure includes both the data pointer and the size, ensuring */
/*   that operations on vectors can verify dimensional compatibility. */

/*   | Variable | Type      | Description             | */
/*   |----------+-----------+-------------------------| */
/*   | ~size~   | ~int64_t~ | Dimension of the vector | */
/*   | ~data~   | ~double*~ | Elements                | */


typedef struct qmckl_vector {
  double* restrict data;
  int64_t size;
} qmckl_vector;

/* Matrix */

/*   | Variable | Type         | Description                 | */
/*   |----------+--------------+-----------------------------| */
/*   | ~size~   | ~int64_t[2]~ | Dimension of each component | */
/*   | ~data~   | ~double*~    | Elements                    | */

/*   The dimensions use Fortran ordering: two elements differing by one */
/*   in the first dimension are consecutive in memory. */


typedef struct qmckl_matrix {
  double* restrict data;
  int64_t size[2];
} qmckl_matrix;

/* Tensor */

/*   | Variable | Type                              | Description                 | */
/*   |----------+-----------------------------------+-----------------------------| */
/*   | ~order~  | ~int32_t~                         | Order of the tensor         | */
/*   | ~size~   | ~int64_t[QMCKL_TENSOR_ORDER_MAX]~ | Dimension of each component | */
/*   | ~data~   | ~double*~                         | Elements                    | */

/*   The dimensions use Fortran ordering: two elements differing by one */
/*   in the first dimension are consecutive in memory. */


#define QMCKL_TENSOR_ORDER_MAX 16

typedef struct qmckl_tensor {
  double* restrict data;
  int32_t order;
  int32_t __space__;
  int64_t size[QMCKL_TENSOR_ORDER_MAX];
} qmckl_tensor;

/* ~qmckl_dgemm_sparse_nn~ */

/*    Sparse matrix multiplication using CSR (Compressed Sparse Row) format for both matrices: */

/*    \[ */
/*    C_{ij} = \beta C_{ij} + \alpha \sum_{k} A_{ik} \cdot B_{kj} */
/*    \] */

/*    This function performs matrix multiplication where both $A$ and $B$ are stored in */
/*    CSR (Compressed Sparse Row) format. The function automatically switches between */
/*    sparse and dense BLAS implementations based on the matrix density threshold. */

/*    #+NAME: qmckl_dgemm_sparse_nn_args */
/*    | Variable  | Type            | In/Out | Description                                    | */
/*    |-----------+-----------------+--------+------------------------------------------------| */
/*    | ~context~ | ~qmckl_context~ | in     | Global state                                   | */
/*    | ~m~       | ~int64_t~       | in     | Number of rows of matrix $A$ and $C$           | */
/*    | ~n~       | ~int64_t~       | in     | Number of columns of matrix $B$ and $C$        | */
/*    | ~k~       | ~int64_t~       | in     | Number of columns of $A$ and rows of $B$       | */
/*    | ~alpha~   | ~double~        | in     | Scalar $\alpha$                                | */
/*    | ~A~       | ~CSRMatrix*~    | in     | CSR format sparse matrix $A$ ($m \times k$)    | */
/*    | ~B~       | ~CSRMatrix*~    | in     | CSR format sparse matrix $B$ ($k \times n$)    | */
/*    | ~beta~    | ~double~        | in     | Scalar $\beta$                                 | */
/*    | ~C~       | ~double*~       | out    | Dense output matrix $C$ ($m \times n$)         | */
/*    | ~ldc~     | ~int64_t~       | in     | Leading dimension of array ~C~                 | */

/*    Requirements: */

/*     - ~context~ is not ~QMCKL_NULL_CONTEXT~ */
/*     - ~m > 0~ */
/*     - ~n > 0~ */
/*     - ~k > 0~ */
/*     - ~A~ is not ~NULL~ and properly initialized as CSR matrix */
/*     - ~B~ is not ~NULL~ and properly initialized as CSR matrix */
/*     - ~C~ is allocated with at least $m \times n$ elements */
/*     - ~ldc >= m~ */


typedef struct {
  int64_t nrows, ncols, nnz, lda;
  int64_t *row_ptr;   // size nrows+1
  int64_t *col_idx;   // size nnz
  double  *val;       // size nnz
  const double *A;
  double density;
} CSRMatrix;


typedef struct {
  int64_t nrows, ncols, nnz, lda;
  int64_t *col_ptr;   // size ncols+1
  int64_t *row_idx;   // size nnz
  double  *val;       // size nnz
  const double *A;
  double density;
} CSCMatrix;

#endif
