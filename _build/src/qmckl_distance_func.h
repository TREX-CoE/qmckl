/* ~qmckl_distance_sq~ */
/*    :PROPERTIES: */
/*    :Name:     qmckl_distance_sq */
/*    :CRetType: qmckl_exit_code */
/*    :FRetType: qmckl_exit_code */
/*    :END: */

/*    The ~qmckl_distance_sq~ function computes the matrix of squared Euclidean */
/*    distances between all pairs of points from two sets. For each point in set A */
/*    and each point in set B, it calculates the squared distance, which avoids the */
/*    computational cost of square root operations that would be needed for actual */
/*    distances. In many quantum chemistry applications, squared distances are */
/*    directly needed (e.g., in Gaussian functions), making this an efficient choice. */
   
/*    The function computes: */

/*    \[ */
/*    C_{ij} = \sum_{k=1}^3 (A_{k,i}-B_{k,j})^2 */
/*    \] */
   
/*    where $A$ and $B$ are sets of 3D points and $C$ is the resulting distance matrix. */
/*    The function supports both normal and transposed input matrices through the */
/*    ~transa~ and ~transb~ parameters, allowing efficient handling of data in */
/*    different memory layouts. */

/*    #+NAME: qmckl_distance_sq_args */
/*    | Variable  | Type             | In/Out | Description                                   | */
/*    |-----------+------------------+--------+-----------------------------------------------| */
/*    | ~context~ | ~qmckl_context~  | in     | Global state                                  | */
/*    | ~transa~  | ~char~           | in     | Array ~A~ is ~'N'~: Normal, ~'T'~: Transposed | */
/*    | ~transb~  | ~char~           | in     | Array ~B~ is ~'N'~: Normal, ~'T'~: Transposed | */
/*    | ~m~       | ~int64_t~        | in     | Number of points in the first set             | */
/*    | ~n~       | ~int64_t~        | in     | Number of points in the second set            | */
/*    | ~A~       | ~double[][lda]~  | in     | Array containing the $m \times 3$ matrix $A$  | */
/*    | ~lda~     | ~int64_t~        | in     | Leading dimension of array ~A~                | */
/*    | ~B~       | ~double[][ldb]~  | in     | Array containing the $n \times 3$ matrix $B$  | */
/*    | ~ldb~     | ~int64_t~        | in     | Leading dimension of array ~B~                | */
/*    | ~C~       | ~double[n][ldc]~ | out    | Array containing the $m \times n$ matrix $C$  | */
/*    | ~ldc~     | ~int64_t~        | in     | Leading dimension of array ~C~                | */

/*    Requirements: */

/*     - ~context~ is not ~QMCKL_NULL_CONTEXT~ */
/*     - ~m > 0~ */
/*     - ~n > 0~ */
/*     - ~lda >= 3~ if ~transa == 'N'~ */
/*     - ~lda >= m~ if ~transa == 'T'~ */
/*     - ~ldb >= 3~ if ~transb == 'N'~ */
/*     - ~ldb >= n~ if ~transb == 'T'~ */
/*     - ~ldc >= m~ */
/*     - ~A~ is allocated with at least $3 \times m \times 8$ bytes */
/*     - ~B~ is allocated with at least $3 \times n \times 8$ bytes */
/*     - ~C~ is allocated with at least $m \times n \times 8$ bytes */

/*     #+CALL: generate_c_header(table=qmckl_distance_sq_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */

/*     #+RESULTS: */

qmckl_exit_code qmckl_distance_sq (
      const qmckl_context context,
      const char transa,
      const char transb,
      const int64_t m,
      const int64_t n,
      const double* A,
      const int64_t lda,
      const double* B,
      const int64_t ldb,
      double* const C,
      const int64_t ldc );

/* C header */

/*     #+CALL: generate_c_header(table=qmckl_distance_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */

/*     #+RESULTS: */

qmckl_exit_code qmckl_distance (
      const qmckl_context context,
      const char transa,
      const char transb,
      const int64_t m,
      const int64_t n,
      const double* A,
      const int64_t lda,
      const double* B,
      const int64_t ldb,
      double* const C,
      const int64_t ldc );

/* C header */

/*     #+CALL: generate_c_header(table=qmckl_distance_rescaled_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */

/*     #+RESULTS: */

qmckl_exit_code qmckl_distance_rescaled (
      const qmckl_context context,
      const char transa,
      const char transb,
      const int64_t m,
      const int64_t n,
      const double* A,
      const int64_t lda,
      const double* B,
      const int64_t ldb,
      double* const C,
      const int64_t ldc,
      const double rescale_factor_kappa );

/* ~qmckl_distance_rescaled_gl~ */
/*    :PROPERTIES: */
/*    :Name:     qmckl_distance_rescaled_gl */
/*    :CRetType: qmckl_exit_code */
/*    :FRetType: qmckl_exit_code */
/*    :END: */

/*    ~qmckl_distance_rescaled_gl~ computes the matrix of the gradient and Laplacian of the */
/*    rescaled distance with respect to the electron coordinates. The derivative is a rank 3 tensor. */
/*    The first dimension has a dimension of 4 of which the first three coordinates */
/*    contains the gradient vector and the last index is the Laplacian. */


/*    \[ */
/*    C(r_{ij}) = \left( 1 - \exp(-\kappa\, r_{ij})\right)/\kappa */
/*    \] */

/*    Here the gradient is defined as follows: */

/*    \[ */
/*    \nabla_i C(r_{ij}) = \left(\frac{\partial C(r_{ij})}{\partial x_i},\frac{\partial C(r_{ij})}{\partial y_i},\frac{\partial C(r_{ij})}{\partial z_i} \right) */
/*    \] */
/*    and the Laplacian is defined as follows: */

/*    \[ */
/*    \Delta_i C(r_{ij}) = \frac{\partial^2}{\partial x_i^2} + \frac{\partial^2}{\partial y_i^2} + \frac{\partial^2}{\partial z_i^2} */
/*    \] */

/*    Using the above three formulas, the expression for the gradient and Laplacian is: */

/*    \[ */
/*    \frac{\partial  C(r_{ij})}{\partial x_i} = \frac{|(x_i - */
/*    x_j)|}{r_{ij}} \exp (- \kappa \, r_{ij}) */
/*    \] */

/*    \[ */
/*    \Delta C_{ij}(r_{ij}) = \left[ \frac{2}{r_{ij}} - \kappa  \right] \exp (- \kappa \, r_{ij}) */
/*    \] */

/*    If the input array is normal (~'N'~), the xyz coordinates are in */
/*    the leading dimension: ~[n][3]~ in C and ~(3,n)~ in Fortran. */

/*    #+NAME: qmckl_distance_rescaled_gl_args */
/*    | Variable               | Type                | In/Out | Description                                           | */
/*    |------------------------+---------------------+--------+-------------------------------------------------------| */
/*    | ~context~              | ~qmckl_context~     | in     | Global state                                          | */
/*    | ~transa~               | ~char~              | in     | Array ~A~ is ~'N'~: Normal, ~'T'~: Transposed         | */
/*    | ~transb~               | ~char~              | in     | Array ~B~ is ~'N'~: Normal, ~'T'~: Transposed         | */
/*    | ~m~                    | ~int64_t~           | in     | Number of points in the first set                     | */
/*    | ~n~                    | ~int64_t~           | in     | Number of points in the second set                    | */
/*    | ~A~                    | ~double[][lda]~     | in     | Array containing the $m \times 3$ matrix $A$          | */
/*    | ~lda~                  | ~int64_t~           | in     | Leading dimension of array ~A~                        | */
/*    | ~B~                    | ~double[][ldb]~     | in     | Array containing the $n \times 3$ matrix $B$          | */
/*    | ~ldb~                  | ~int64_t~           | in     | Leading dimension of array ~B~                        | */
/*    | ~C~                    | ~double[n][ldc][4]~ | out    | Array containing the $4 \times m \times n$ matrix $C$ | */
/*    | ~ldc~                  | ~int64_t~           | in     | Leading dimension of array ~C~                        | */
/*    | ~rescale_factor_kappa~ | ~double~            | in     | Factor for calculating rescaled distances derivatives | */

/*    Requirements: */

/*     - ~context~ is not ~QMCKL_NULL_CONTEXT~ */
/*     - ~m > 0~ */
/*     - ~n > 0~ */
/*     - ~lda >= 3~ if ~transa == 'N'~ */
/*     - ~lda >= m~ if ~transa == 'T'~ */
/*     - ~ldb >= 3~ if ~transb == 'N'~ */
/*     - ~ldb >= n~ if ~transb == 'T'~ */
/*     - ~ldc >= m~ */
/*     - ~A~ is allocated with at least $3 \times m \times 8$ bytes */
/*     - ~B~ is allocated with at least $3 \times n \times 8$ bytes */
/*     - ~C~ is allocated with at least $4 \times m \times n \times 8$ bytes */

/*     #+CALL: generate_c_header(table=qmckl_distance_rescaled_gl_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */

/*     #+RESULTS: */

qmckl_exit_code qmckl_distance_rescaled_gl (
      const qmckl_context context,
      const char transa,
      const char transb,
      const int64_t m,
      const int64_t n,
      const double* A,
      const int64_t lda,
      const double* B,
      const int64_t ldb,
      double* const C,
      const int64_t ldc,
      const double rescale_factor_kappa );
