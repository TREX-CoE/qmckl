/* ~qmckl_dgemm~ */

/*    Matrix multiplication with a BLAS interface: */

/*    \[ */
/*    C_{ij} = \beta C_{ij} + \alpha \sum_{k} A_{ik} \cdot B_{kj} */
/*    \] */

/* #  TODO: Add description about the external library dependence. */

/*    #+NAME: qmckl_dgemm_args */
/*    | Variable  | Type            | In/Out | Description                           | */
/*    |-----------+-----------------+--------+---------------------------------------| */
/*    | ~context~ | ~qmckl_context~ | in     | Global state                          | */
/*    | ~TransA~  | ~char~          | in     | 'T' is transposed                     | */
/*    | ~TransB~  | ~char~          | in     | 'T' is transposed                     | */
/*    | ~m~       | ~int64_t~       | in     | Number of rows of the input matrix    | */
/*    | ~n~       | ~int64_t~       | in     | Number of columns of the input matrix | */
/*    | ~k~       | ~int64_t~       | in     | Number of columns of the input matrix | */
/*    | ~alpha~   | ~double~        | in     | \alpha                                | */
/*    | ~A~       | ~double[][lda]~ | in     | Array containing the matrix $A$       | */
/*    | ~lda~     | ~int64_t~       | in     | Leading dimension of array ~A~        | */
/*    | ~B~       | ~double[][ldb]~ | in     | Array containing the matrix $B$       | */
/*    | ~ldb~     | ~int64_t~       | in     | Leading dimension of array ~B~        | */
/*    | ~beta~    | ~double~        | in     | \beta                                 | */
/*    | ~C~       | ~double[][ldc]~ | out    | Array containing the matrix $C$       | */
/*    | ~ldc~     | ~int64_t~       | in     | Leading dimension of array ~C~        | */

/*    Requirements: */

/*     - ~context~ is not ~QMCKL_NULL_CONTEXT~ */
/*     - ~m > 0~ */
/*     - ~n > 0~ */
/*     - ~k > 0~ */
/*     - ~lda >= m~ */
/*     - ~ldb >= n~ */
/*     - ~ldc >= n~ */
/*     - ~A~ is allocated with at least $m \times k \times 8$ bytes */
/*     - ~B~ is allocated with at least $k \times n \times 8$ bytes */
/*     - ~C~ is allocated with at least $m \times n \times 8$ bytes */

/*     #+CALL: generate_c_header(table=qmckl_dgemm_args,rettyp="qmckl_exit_code",fname="qmckl_dgemm") */

/*     #+RESULTS: */

qmckl_exit_code qmckl_dgemm (
      const qmckl_context context,
      const char TransA,
      const char TransB,
      const int64_t m,
      const int64_t n,
      const int64_t k,
      const double alpha,
      const double* A,
      const int64_t lda,
      const double* B,
      const int64_t ldb,
      const double beta,
      double* const C,
      const int64_t ldc );

/* ~qmckl_dgemm_safe~ */

/*    "Size-safe" proxy function with the same functionality as ~qmckl_dgemm~ */
/*    but with 3 additional arguments. These arguments ~size_max_M~ (where ~M~ is a matix) */
/*    are required primarily for the Python API, where compatibility with */
/*    NumPy arrays implies that sizes of the input and output arrays are provided. */

/*    #+NAME: qmckl_dgemm_safe_args */
/*    | Variable     | Type            | In/Out | Description                           | */
/*    |--------------+-----------------+--------+---------------------------------------| */
/*    | ~context~    | ~qmckl_context~ | in     | Global state                          | */
/*    | ~TransA~     | ~char~          | in     | 'T' is transposed                     | */
/*    | ~TransB~     | ~char~          | in     | 'T' is transposed                     | */
/*    | ~m~          | ~int64_t~       | in     | Number of rows of the input matrix    | */
/*    | ~n~          | ~int64_t~       | in     | Number of columns of the input matrix | */
/*    | ~k~          | ~int64_t~       | in     | Number of columns of the input matrix | */
/*    | ~alpha~      | ~double~        | in     | \alpha                                | */
/*    | ~A~          | ~double[][lda]~ | in     | Array containing the matrix $A$       | */
/*    | ~size_max_A~ | ~int64_t~       | in     | Size of the matrix $A$                | */
/*    | ~lda~        | ~int64_t~       | in     | Leading dimension of array ~A~        | */
/*    | ~B~          | ~double[][ldb]~ | in     | Array containing the matrix $B$       | */
/*    | ~size_max_B~ | ~int64_t~       | in     | Size of the matrix $B$                | */
/*    | ~ldb~        | ~int64_t~       | in     | Leading dimension of array ~B~        | */
/*    | ~beta~       | ~double~        | in     | \beta                                 | */
/*    | ~C~          | ~double[][ldc]~ | out    | Array containing the matrix $C$       | */
/*    | ~size_max_C~ | ~int64_t~       | in     | Size of the matrix $C$                | */
/*    | ~ldc~        | ~int64_t~       | in     | Leading dimension of array ~C~        | */

/*    Requirements: */

/*     - ~context~ is not ~QMCKL_NULL_CONTEXT~ */
/*     - ~m > 0~ */
/*     - ~n > 0~ */
/*     - ~k > 0~ */
/*     - ~lda >= m~ */
/*     - ~ldb >= n~ */
/*     - ~ldc >= n~ */
/*     - ~A~ is allocated with at least $m \times k \times 8$ bytes */
/*     - ~B~ is allocated with at least $k \times n \times 8$ bytes */
/*     - ~C~ is allocated with at least $m \times n \times 8$ bytes */
/*     - ~size_max_A >= m * k~ */
/*     - ~size_max_B >= k * n~ */
/*     - ~size_max_C >= m * n~ */

/*     #+CALL: generate_c_header(table=qmckl_dgemm_safe_args,rettyp="qmckl_exit_code",fname="qmckl_dgemm_safe") */

/*     #+RESULTS: */

qmckl_exit_code qmckl_dgemm_safe (
      const qmckl_context context,
      const char TransA,
      const char TransB,
      const int64_t m,
      const int64_t n,
      const int64_t k,
      const double alpha,
      const double* A,
      const int64_t size_max_A,
      const int64_t lda,
      const double* B,
      const int64_t size_max_B,
      const int64_t ldb,
      const double beta,
      double* const C,
      const int64_t size_max_C,
      const int64_t ldc );

/* ~qmckl_adjugate~ */

/*    Given a matrix $\mathbf{A}$, the adjugate matrix */
/*    $\text{adj}(\mathbf{A})$ is the transpose of the cofactors matrix */
/*    of $\mathbf{A}$. */

/*    \[ */
/*    \mathbf{B} = \text{adj}(\mathbf{A}) = \text{det}(\mathbf{A}) \, \mathbf{A}^{-1} */
/*    \] */

/*    See also: https://en.wikipedia.org/wiki/Adjugate_matrix */

/*    #+NAME: qmckl_adjugate_args */
/*    | Variable  | Type            | In/Out | Description                                    | */
/*    |-----------+-----------------+--------+------------------------------------------------| */
/*    | ~context~ | ~qmckl_context~ | in     | Global state                                   | */
/*    | ~n~       | ~int64_t~       | in     | Number of rows and columns of the input matrix | */
/*    | ~A~       | ~double[][lda]~ | in     | Array containing the $n \times n$ matrix $A$   | */
/*    | ~lda~     | ~int64_t~       | in     | Leading dimension of array ~A~                 | */
/*    | ~B~       | ~double[][ldb]~ | out    | Adjugate of $A$                                | */
/*    | ~ldb~     | ~int64_t~       | in     | Leading dimension of array ~B~                 | */
/*    | ~det_l~   | ~double~        | inout  | determinant of $A$                             | */

/*    Requirements: */

/*     - ~context~ is not ~QMCKL_NULL_CONTEXT~ */
/*     - ~n > 0~ */
/*     - ~lda >= m~ */
/*     - ~A~ is allocated with at least $m \times m \times 8$ bytes */
/*     - ~ldb >= m~ */
/*     - ~B~ is allocated with at least $m \times m \times 8$ bytes */

/*     #+CALL: generate_c_header(table=qmckl_adjugate_args,rettyp="qmckl_exit_code",fname="qmckl_adjugate") */

/*     #+RESULTS: */

qmckl_exit_code qmckl_adjugate (
      const qmckl_context context,
      const int64_t n,
      const double* A,
      const int64_t lda,
      double* const B,
      const int64_t ldb,
      double* det_l );

/* ~qmckl_adjugate_safe~ */

/*    "Size-safe" proxy function with the same functionality as ~qmckl_adjugate~ */
/*    but with 2 additional arguments. These arguments ~size_max_M~ (where ~M~ is a matix) */
/*    are required primarily for the Python API, where compatibility with */
/*    NumPy arrays implies that sizes of the input and output arrays are provided. */


/*    #+NAME: qmckl_adjugate_safe_args */
/*    | Variable     | Type            | In/Out | Description                                    | */
/*    |--------------+-----------------+--------+------------------------------------------------| */
/*    | ~context~    | ~qmckl_context~ | in     | Global state                                   | */
/*    | ~n~          | ~int64_t~       | in     | Number of rows and columns of the input matrix | */
/*    | ~A~          | ~double[][lda]~ | in     | Array containing the $n \times n$ matrix $A$   | */
/*    | ~size_max_A~ | ~int64_t~       | in     | Size of the matrix $A$                         | */
/*    | ~lda~        | ~int64_t~       | in     | Leading dimension of array ~A~                 | */
/*    | ~B~          | ~double[][ldb]~ | out    | Adjugate of $A$                                | */
/*    | ~size_max_B~ | ~int64_t~       | in     | Size of the matrix $B$                         | */
/*    | ~ldb~        | ~int64_t~       | in     | Leading dimension of array ~B~                 | */
/*    | ~det_l~      | ~double~        | inout  | determinant of $A$                             | */

/*    Requirements: */

/*     - ~context~ is not ~QMCKL_NULL_CONTEXT~ */
/*     - ~n > 0~ */
/*     - ~lda >= m~ */
/*     - ~A~ is allocated with at least $m \times m \times 8$ bytes */
/*     - ~ldb >= m~ */
/*     - ~B~ is allocated with at least $m \times m \times 8$ bytes */
/*     - ~size_max_A >= m * m~ */
/*     - ~size_max_B >= m * m~ */

/*     #+CALL: generate_c_header(table=qmckl_adjugate_safe_args,rettyp="qmckl_exit_code",fname="qmckl_adjugate_safe") */

/*     #+RESULTS: */

qmckl_exit_code qmckl_adjugate_safe (
      const qmckl_context context,
      const int64_t n,
      const double* A,
      const int64_t size_max_A,
      const int64_t lda,
      double* const B,
      const int64_t size_max_B,
      const int64_t ldb,
      double* det_l );

/* ~qmckl_dgemv~ */

/*    Matrix-vector multiplication with a BLAS interface: */

/*    \[ */
/*    Y_{i} = \sum_{j} \alpha A_{ij} X_j + \beta Y_i */
/*    \] */

/*    #+NAME: qmckl_dgemv_args */
/*    | Variable  | Type            | In/Out | Description                           | */
/*    |-----------+-----------------+--------+---------------------------------------| */
/*    | ~context~ | ~qmckl_context~ | in     | Global state                          | */
/*    | ~Trans~   | ~char~          | in     | 'T' is transposed                     | */
/*    | ~m~       | ~int64_t~       | in     | Number of rows of the input matrix    | */
/*    | ~n~       | ~int64_t~       | in     | Number of columns of the input matrix | */
/*    | ~alpha~   | ~double~        | in     | $\alpha$                              | */
/*    | ~A~       | ~double[][lda]~ | in     | Array containing the matrix $A$       | */
/*    | ~lda~     | ~int64_t~       | in     | Leading dimension of array ~A~        | */
/*    | ~X~       | ~double[]~      | in     | Array containing the vector $x$       | */
/*    | ~incx~    | ~int64_t~       | in     | Increment of the elements of ~x~      | */
/*    | ~beta~    | ~double~        | in     | $\beta$                               | */
/*    | ~y~       | ~double[]~      | in     | Array containing the vector $y$       | */
/*    | ~incy~    | ~int64_t~       | in     | Increment of the elements of ~y~      | */
/*     #+CALL: generate_c_header(table=qmckl_dgemv_args,rettyp="qmckl_exit_code",fname="qmckl_dgemv") */

/*     #+RESULTS: */

qmckl_exit_code qmckl_dgemv (
      const qmckl_context context,
      const char Trans,
      const int64_t m,
      const int64_t n,
      const double alpha,
      const double* A,
      const int64_t lda,
      const double* X,
      const int64_t incx,
      const double beta,
      double* const Y,
      const int64_t incy);

/* ~qmckl_dger~ */

/*    Rank-1 update */

/*    \[ */
/*    A = \alpha x y^{\dagger} + A */
/*    \] */

/*    #+NAME: qmckl_dger_args */
/*    | Variable  | Type            | In/Out | Description                           | */
/*    |-----------+-----------------+--------+---------------------------------------| */
/*    | ~context~ | ~qmckl_context~ | in     | Global state                          | */
/*    | ~m~       | ~int64_t~       | in     | Number of rows of the input matrix    | */
/*    | ~n~       | ~int64_t~       | in     | Number of columns of the input matrix | */
/*    | ~alpha~   | ~double~        | in     | $\alpha$                              | */
/*    | ~X~       | ~double[]~      | in     | Array containing the vector $x$       | */
/*    | ~incx~    | ~int64_t~       | in     | Increment of the elements of ~x~      | */
/*    | ~y~       | ~double[]~      | in     | Array containing the vector $y$       | */
/*    | ~incy~    | ~int64_t~       | in     | Increment of the elements of ~y~      | */
/*    | ~A~       | ~double[][lda]~ | inout  | Array containing the matrix $A$       | */
/*    | ~lda~     | ~int64_t~       | in     | Leading dimension of array ~A~        | */

/*     #+CALL: generate_c_header(table=qmckl_dger_args,rettyp="qmckl_exit_code",fname="qmckl_dger") */
/*     #+RESULTS: */

qmckl_exit_code qmckl_dger (
      const qmckl_context context,
      const int64_t m,
      const int64_t n,
      const double alpha,
      const double* X,
      const int64_t incx,
      const double* Y,
      const int64_t incy,
      double* const A,
      const int64_t lda);
