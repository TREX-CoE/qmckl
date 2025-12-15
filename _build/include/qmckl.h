/*
 *    ------------------------------------------
 *     QMCkl - Quantum Monte Carlo kernel library
 *     ------------------------------------------
 *
 *     Documentation : https://trex-coe.github.io/qmckl
 *     Issues        : https://github.com/trex-coe/qmckl/issues
 *
 *     BSD 3-Clause License
 *
 *     Copyright (c) 2020, TREX Center of Excellence
 *     All rights reserved.
 *
 *     Redistribution and use in source and binary forms, with or without
 *     modification, are permitted provided that the following conditions are met:
 *
 *     1. Redistributions of source code must retain the above copyright notice, this
 *        list of conditions and the following disclaimer.
 *
 *     2. Redistributions in binary form must reproduce the above copyright notice,
 *        this list of conditions and the following disclaimer in the documentation
 *        and/or other materials provided with the distribution.
 *
 *     3. Neither the name of the copyright holder nor the names of its
 *        contributors may be used to endorse or promote products derived from
 *        this software without specific prior written permission.
 *
 *     THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *     AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *     IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *     DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 *     FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 *     DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 *     SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *     CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 *     OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *     OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 *
 *
 */

#ifndef __QMCKL_H__
#define __QMCKL_H__

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>
/* /home/runner/work/qmckl/qmckl/_build/src/qmckl_context_type.h */
/* Context handling */

/*   The context is the central data structure in QMCkl, serving as a handle */
/*   for the complete state of the library. All QMCkl functions require a context */
/*   as their first argument, and all computed data is stored within the context. */
  
/*   The context variable is a handle for the state of the library, and is stored  */
/*   in a data structure which can't be seen outside of the library. This  */
/*   encapsulation provides a clean API boundary and allows internal implementation  */
/*   details to change without affecting user code. */
  
/*   To simplify compatibility with other languages (particularly Fortran), the */
/*   pointer to the internal data structure is converted into a 64-bit */
/*   signed integer, defined in the ~qmckl_context~ type. This approach avoids */
/*   issues with language interoperability related to opaque pointer types. */
/*   A value of ~QMCKL_NULL_CONTEXT~ for the context is equivalent to a */
/*   ~NULL~ pointer and represents an invalid or uninitialized context. */

/*   #+NAME: qmckl_context */

typedef int64_t qmckl_context ;
#define QMCKL_NULL_CONTEXT (qmckl_context) 0
/* /home/runner/work/qmckl/qmckl/_build/src/qmckl_error_type.h */
/* Exit codes */

/*    All the functions in the QMCkl library return with an exit code, defined as */
/*    #+NAME: type-exit-code */

typedef int32_t qmckl_exit_code;



/* #+RESULTS: */
/* :results: */

#define  QMCKL_SUCCESS                  ((qmckl_exit_code) 0)
#define  QMCKL_INVALID_ARG_1            ((qmckl_exit_code) 1)
#define  QMCKL_INVALID_ARG_2            ((qmckl_exit_code) 2)
#define  QMCKL_INVALID_ARG_3            ((qmckl_exit_code) 3)
#define  QMCKL_INVALID_ARG_4            ((qmckl_exit_code) 4)
#define  QMCKL_INVALID_ARG_5            ((qmckl_exit_code) 5)
#define  QMCKL_INVALID_ARG_6            ((qmckl_exit_code) 6)
#define  QMCKL_INVALID_ARG_7            ((qmckl_exit_code) 7)
#define  QMCKL_INVALID_ARG_8            ((qmckl_exit_code) 8)
#define  QMCKL_INVALID_ARG_9            ((qmckl_exit_code) 9)
#define  QMCKL_INVALID_ARG_10           ((qmckl_exit_code) 10)
#define  QMCKL_INVALID_ARG_11           ((qmckl_exit_code) 11)
#define  QMCKL_INVALID_ARG_12           ((qmckl_exit_code) 12)
#define  QMCKL_INVALID_ARG_13           ((qmckl_exit_code) 13)
#define  QMCKL_INVALID_ARG_14           ((qmckl_exit_code) 14)
#define  QMCKL_INVALID_ARG_15           ((qmckl_exit_code) 15)
#define  QMCKL_INVALID_ARG_16           ((qmckl_exit_code) 16)
#define  QMCKL_INVALID_ARG_17           ((qmckl_exit_code) 17)
#define  QMCKL_INVALID_ARG_18           ((qmckl_exit_code) 18)
#define  QMCKL_INVALID_ARG_19           ((qmckl_exit_code) 19)
#define  QMCKL_INVALID_ARG_20           ((qmckl_exit_code) 20)
#define  QMCKL_FAILURE                  ((qmckl_exit_code) 101)
#define  QMCKL_ERRNO                    ((qmckl_exit_code) 102)
#define  QMCKL_INVALID_CONTEXT          ((qmckl_exit_code) 103)
#define  QMCKL_ALLOCATION_FAILED        ((qmckl_exit_code) 104)
#define  QMCKL_DEALLOCATION_FAILED      ((qmckl_exit_code) 105)
#define  QMCKL_NOT_PROVIDED             ((qmckl_exit_code) 106)
#define  QMCKL_OUT_OF_BOUNDS            ((qmckl_exit_code) 107)
#define  QMCKL_ALREADY_SET              ((qmckl_exit_code) 108)
#define  QMCKL_INVALID_EXIT_CODE        ((qmckl_exit_code) 109)
/* /home/runner/work/qmckl/qmckl/_build/src/qmckl_numprec_type.h */


/* #+RESULTS: */
/* :results: */

#define  QMCKL_DEFAULT_PRECISION        53
#define  QMCKL_DEFAULT_RANGE            11
/* /home/runner/work/qmckl/qmckl/_build/src/qmckl_context_func.h */


/* The ~qmckl_context_check~ function checks if the domain pointed to by */
/* the pointer is a valid context by verifying the magic tag value. It returns  */
/* the input ~qmckl_context~ if the context is valid (tag matches ~VALID_TAG~),  */
/* or ~QMCKL_NULL_CONTEXT~ if the context is invalid. This function should be  */
/* called at the beginning of every public API function to ensure the context */
/* parameter is valid before attempting to use it. */


qmckl_context
qmckl_context_check (const qmckl_context context) ;



/* The context keeps a /date/ that allows to check which data needs */
/* to be recomputed. The date is incremented when the context is touched. */

/* When a new element is added to the context, the functions */
/* [[Creation][=qmckl_context_create=]] [[Destroy][=qmckl_context_destroy=]] and [[Copy][=qmckl_context_copy=]] */
/* should be updated in order to make deep copies. */

/* When the electron coordinates have changed, the context is touched */
/* using the following function. */


qmckl_exit_code
qmckl_context_touch (const qmckl_context context);

/* Creation */

/*    To create a new context, ~qmckl_context_create()~ should be used. */
/*    - Upon success, it returns a pointer to a new context with the ~qmckl_context~ type */
/*    - It returns ~QMCKL_NULL_CONTEXT~ upon failure to allocate the internal data structure */
/*    - A new context always has all its members initialized with a NULL value */

/*    # Header */

qmckl_context qmckl_context_create();

/* Locking */

/*    For thread safety, the context may be locked/unlocked. The lock is */
/*    initialized with the ~PTHREAD_MUTEX_RECURSIVE~ attribute, and the */
/*    number of times the thread has locked it is saved in the */
/*    ~lock_count~ attribute. */

/*    # Header */

void qmckl_lock  (qmckl_context context);
void qmckl_unlock(qmckl_context context);

/* Copy */

/*    ~qmckl_context_copy~ makes a deep copy of a context. It returns */
/*    ~QMCKL_NULL_CONTEXT~ upon failure. */

/*    # Header */

qmckl_context qmckl_context_copy(const qmckl_context context);

/* Destroy */

/*    The context is destroyed with ~qmckl_context_destroy~, leaving the ancestors untouched. */
/*    It frees the context, and returns the previous context. */

/*    # Header */

qmckl_exit_code
qmckl_context_destroy (const qmckl_context context);
/* /home/runner/work/qmckl/qmckl/_build/src/qmckl_error_func.h */
/* Decoding errors */

/*    To facilitate debugging and error reporting, QMCkl provides the  */
/*    ~qmckl_string_of_error~ function which converts an error code into  */
/*    a descriptive string. This allows applications to present meaningful  */
/*    error messages to users or write detailed logs for troubleshooting. */
   
/*    The function takes an error code as input and returns a constant string */
/*    describing the error condition. Both C and Fortran interfaces are provided */
/*    for maximum compatibility. */


const char*
qmckl_string_of_error (const qmckl_exit_code error);

/* Updating errors in the context */

/*    The error state in the context is updated using the ~qmckl_set_error~ */
/*    function. This function provides a centralized mechanism for recording */
/*    errors that occur during library operations. */
   
/*    When an error is set in the context, it is mandatory to specify */
/*    from which function the error is triggered, and a message */
/*    explaining the error. The exit code can't be ~QMCKL_SUCCESS~. */
   
/*    This detailed error information enables precise debugging and helps */
/*    users understand exactly what went wrong and where, making it easier */
/*    to diagnose and fix issues in code using the QMCkl library. */

/*    # Header */

qmckl_exit_code
qmckl_set_error(qmckl_context context,
                const qmckl_exit_code exit_code,
                const char* function_name,
                const char* message);

/* Get the error */

/*   Upon error, the calling program can retrieve detailed error information from the */
/*   context using ~qmckl_get_error~. This function provides access to the error  */
/*   code, the name of the function where the error occurred, and a descriptive  */
/*   message explaining the error condition. */
  
/*   The error message and function name are returned in the variables provided  */
/*   by the caller. Therefore, passing valid pointers for the function name and  */
/*   message is mandatory. The caller must ensure that the provided buffers are */
/*   large enough to hold the error information. */
  
/*   This retrieval mechanism allows applications to implement sophisticated error */
/*   handling strategies, such as retry logic, fallback mechanisms, or detailed */
/*   logging for post-mortem analysis. */

/*   # Header */

qmckl_exit_code
qmckl_get_error(qmckl_context context,
                qmckl_exit_code *exit_code,
                char* function_name,
                char* message);

/* Failing */

/*    To make a function fail, the ~qmckl_failwith~ function should be */
/*    called, such that information about the failure is stored in */
/*    the context. The desired exit code is given as an argument, as */
/*    well as the name of the function and an error message. If the */
/*    message is ~NULL~, then the default message obtained by */
/*    ~qmckl_string_of_error~ is used. The return code of the function is */
/*    the desired return code. */
/*    Upon failure, a ~QMCKL_NULL_CONTEXT~ is returned. */


qmckl_exit_code
qmckl_failwith(qmckl_context context,
               const qmckl_exit_code exit_code,
               const char* function,
               const char* message) ;

/* Last error */

/*   Returns a string describing the last error, using ~qmckl_get_error~. */

/*   # Header */

qmckl_exit_code
qmckl_last_error(qmckl_context context, char* buffer);

/* Helper functions for debugging */

/*   The following function prints to ~stderr~ an error message is the return code is */
/*   not ~QMCKL_SUCCESS~. */

/*   # Header */

qmckl_exit_code
qmckl_check(qmckl_context context, qmckl_exit_code rc);
/* /home/runner/work/qmckl/qmckl/_build/src/qmckl_blas_func.h */
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
/* /home/runner/work/qmckl/qmckl/_build/src/qmckl_numprec_func.h */
/* Precision */
/*   ~qmckl_context_set_numprec_precision~ modifies the parameter for the */
/*   numerical precision in the context. */

/*   # Header */

qmckl_exit_code qmckl_set_numprec_precision(const qmckl_context context, const int precision);



/* ~qmckl_get_numprec_precision~ returns the value of the numerical precision in the context. */


int32_t qmckl_get_numprec_precision(const qmckl_context context);

/* Range */

/*    ~qmckl_set_numprec_range~ modifies the parameter for the numerical */
/*    range in a given context. */

/*    # Header */

qmckl_exit_code qmckl_set_numprec_range(const qmckl_context context, const int range);



/* ~qmckl_get_numprec_range~ returns the value of the numerical range in the context. */


int32_t qmckl_get_numprec_range(const qmckl_context context);

/* Epsilon */

/*    ~qmckl_get_numprec_epsilon~ returns $\epsilon = 2^{1-n}$ where ~n~ is the precision. */
/*    We need to remove the sign bit from the precision. */


double qmckl_get_numprec_epsilon(const qmckl_context context);

int32_t qmckl_test_precision_64(double a, double b);
int32_t qmckl_test_precision_32(float a, float b);
/* /home/runner/work/qmckl/qmckl/_build/src/qmckl_point_func.h */
/* Number of points */


qmckl_exit_code qmckl_get_point_num (const qmckl_context context, int64_t* const num);

/* Point coordinates */


qmckl_exit_code qmckl_get_point(const qmckl_context context,
                                const char transp,
                                double* const coord,
                                const int64_t size_max);

/* Initialization functions */

/*    When the data is set in the context, if the arrays are large */
/*    enough, we overwrite the data contained in them. */

/*    To set the data relative to the points in the context, the */
/*    following function need to be called. Here, ~num~ is the number of */
/*    points to set. */


qmckl_exit_code qmckl_set_point (qmckl_context context,
                                 const char transp,
                                 const int64_t num,
                                 const double* coord,
                                 const int64_t size_max);
/* /home/runner/work/qmckl/qmckl/_build/src/qmckl_nucleus_func.h */
/* Access functions */


qmckl_exit_code
qmckl_get_nucleus_num(const qmckl_context context,
                      int64_t* const num);

qmckl_exit_code
qmckl_get_nucleus_charge(const qmckl_context context,
                         double* const charge,
                         const int64_t size_max);

qmckl_exit_code
qmckl_get_nucleus_coord(const qmckl_context context,
                        const char transp,
                        double* const coord,
                        const int64_t size_max);



/* When all the data relative to nuclei have been set, the following */
/* function returns ~true~. */


bool qmckl_nucleus_provided (const qmckl_context context);



/* To set the data relative to the nuclei in the context, the */
/* following functions need to be called. */


qmckl_exit_code
qmckl_set_nucleus_num(qmckl_context context,
                      const int64_t num);

qmckl_exit_code
qmckl_set_nucleus_charge(qmckl_context context,
                         const double* charge,
                         const int64_t size_max);

qmckl_exit_code
qmckl_set_nucleus_coord(qmckl_context context,
                        const char transp,
                        const double* coord,
                        const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_nucleus_nn_distance(qmckl_context context,
                              double* distance,
                              const int64_t size_max);

/* Get */


qmckl_exit_code qmckl_get_nucleus_repulsion(qmckl_context context, double* const energy);
/* /home/runner/work/qmckl/qmckl/_build/src/qmckl_electron_func.h */
bool qmckl_electron_provided (const qmckl_context context);

/* Initialization functions */

/*    To set the data relative to the electrons in the context, the */
/*    following functions need to be called. When the data structure is */
/*    initialized, the internal ~coord_new~ and ~coord_old~ arrays are */
/*    both not allocated. */


qmckl_exit_code qmckl_set_electron_num      (qmckl_context context, const int64_t up_num, const int64_t down_num);
qmckl_exit_code qmckl_set_electron_coord    (qmckl_context context, const char transp, const int64_t walk_num, const double* coord, const int64_t size_max);

/* Number of electrons */


qmckl_exit_code qmckl_get_electron_num        (const qmckl_context context, int64_t* const num);
qmckl_exit_code qmckl_get_electron_up_num     (const qmckl_context context, int64_t* const up_num);
qmckl_exit_code qmckl_get_electron_down_num   (const qmckl_context context, int64_t* const down_num);

/* Number of walkers */

/*     A walker is a set of electron coordinates that are arguments of */
/*     the wave function. ~walk_num~ is the number of walkers. */


qmckl_exit_code qmckl_get_electron_walk_num   (const qmckl_context context, int64_t* const walk_num);

/* Electron coordinates */

/*     Returns the current electron coordinates. The pointer is assumed */
/*     to point on a memory block of size ~size_max~ \ge ~3 * elec_num * walker.num~. */
/*     The order of the indices is: */

/*     |         | Normal                     | Transposed                 | */
/*     |---------+----------------------------+----------------------------| */
/*     | C       | ~[walker.num*elec_num][3]~ | ~[3][walker.num*elec_num]~ | */
/*     | Fortran | ~(3,walker.num*elec_num)~  | ~(walker.num*elec_num, 3)~ | */



qmckl_exit_code
qmckl_get_electron_coord (const qmckl_context context,
                          const char transp,
                          double* const coord,
                          const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_electron_ee_distance(qmckl_context context,
                               double* const distance,
                               const int64_t size_max);

/* Get */


qmckl_exit_code qmckl_get_electron_ee_potential(qmckl_context context, double* const ee_potential);

/* Get */


qmckl_exit_code
qmckl_get_electron_en_distance(qmckl_context context,
                               double* const distance,
                               const int64_t size_max);

/* Get */


qmckl_exit_code qmckl_get_electron_en_potential(qmckl_context context, double* const en_potential);

qmckl_exit_code qmckl_compute_en_potential (
      const qmckl_context context,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t walk_num,
      const double* charge,
      const double* en_distance,
      double* const en_potential );
/* /home/runner/work/qmckl/qmckl/_build/src/qmckl_distance_func.h */
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
/* /home/runner/work/qmckl/qmckl/_build/src/qmckl_ao_func.h */


/* To set the basis set, all the following functions need to be */
/* called. */


qmckl_exit_code
qmckl_set_ao_basis_type (qmckl_context context,
                         const char basis_type);

qmckl_exit_code
qmckl_set_ao_basis_shell_num (qmckl_context context,
                              const int64_t shell_num);

qmckl_exit_code
qmckl_set_ao_basis_prim_num (qmckl_context context,
                             const int64_t prim_num);

qmckl_exit_code
qmckl_set_ao_basis_nucleus_shell_num (qmckl_context context,
                                      const int64_t* nucleus_shell_num,
                                      const int64_t size_max);

qmckl_exit_code
qmckl_set_ao_basis_nucleus_index (qmckl_context context,
                                  const int64_t* nucleus_index,
                                  const int64_t size_max);

qmckl_exit_code
qmckl_set_ao_basis_shell_ang_mom (qmckl_context context,
                                  const int32_t* shell_ang_mom,
                                  const int64_t size_max);

qmckl_exit_code
qmckl_set_ao_basis_shell_prim_num (qmckl_context context,
                                   const int64_t* shell_prim_num,
                                   const int64_t size_max);

qmckl_exit_code
qmckl_set_ao_basis_shell_prim_index (qmckl_context context,
                                     const int64_t* shell_prim_index,
                                     const int64_t size_max);

qmckl_exit_code
qmckl_set_ao_basis_shell_factor (qmckl_context context,
                                 const double* shell_factor,
                                 const int64_t size_max);

qmckl_exit_code
qmckl_set_ao_basis_exponent (qmckl_context context,
                             const double* exponent,
                             const int64_t size_max);

qmckl_exit_code
qmckl_set_ao_basis_coefficient (qmckl_context context,
                                const double* coefficient,
                                const int64_t size_max);

qmckl_exit_code
qmckl_set_ao_basis_prim_factor (qmckl_context context,
                                const double* prim_factor,
                                const int64_t size_max);

qmckl_exit_code
qmckl_set_ao_basis_ao_num (qmckl_context context,
                           const int64_t ao_num);

qmckl_exit_code
qmckl_set_ao_basis_ao_factor (qmckl_context context,
                              const double* ao_factor,
                              const int64_t size_max);

qmckl_exit_code
qmckl_set_ao_basis_r_power (qmckl_context context,
                            const int32_t* r_power,
                            const int64_t size_max);

qmckl_exit_code
qmckl_set_ao_basis_cartesian (qmckl_context context,
                              const bool cartesian);

/* C interface */


qmckl_exit_code
qmckl_get_ao_basis_type (const qmckl_context context,
                         char* const basis_type);

qmckl_exit_code
qmckl_get_ao_basis_shell_num (const qmckl_context context,
                              int64_t* const shell_num);

qmckl_exit_code
qmckl_get_ao_basis_prim_num (const qmckl_context context,
                             int64_t* const prim_num);

qmckl_exit_code
qmckl_get_ao_basis_nucleus_shell_num (const qmckl_context context,
                                      int64_t* const nucleus_shell_num,
                                      const int64_t size_max);

qmckl_exit_code
qmckl_get_ao_basis_nucleus_index (const qmckl_context context,
                                  int64_t* const nucleus_index,
                                  const int64_t size_max);

qmckl_exit_code
qmckl_get_ao_basis_r_power (const qmckl_context context,
                            int32_t* const r_power,
                            const int64_t size_max);

qmckl_exit_code
qmckl_get_ao_basis_shell_ang_mom (const qmckl_context context,
                                  int32_t* const shell_ang_mom,
                                  const int64_t size_max);

qmckl_exit_code
qmckl_get_ao_basis_shell_prim_num (const qmckl_context context,
                                   int64_t* const shell_prim_num,
                                   const int64_t size_max);

qmckl_exit_code
qmckl_get_ao_basis_shell_prim_index (const qmckl_context context,
                                     int64_t* const shell_prim_index,
                                     const int64_t size_max);

qmckl_exit_code
qmckl_get_ao_basis_shell_factor (const qmckl_context context,
                                 double*  const shell_factor,
                                 const int64_t size_max);

qmckl_exit_code
qmckl_get_ao_basis_exponent (const qmckl_context context,
                             double*  const exponent,
                             const int64_t size_max);

qmckl_exit_code
qmckl_get_ao_basis_coefficient (const qmckl_context context,
                                double*  const coefficient,
                                const int64_t size_max);

qmckl_exit_code
qmckl_get_ao_basis_prim_factor (const qmckl_context context,
                                double*  const prim_factor,
                                const int64_t size_max);

qmckl_exit_code
qmckl_get_ao_basis_ao_num (const qmckl_context context,
                           int64_t* const ao_num);

qmckl_exit_code
qmckl_get_ao_basis_ao_factor (const qmckl_context context,
                              double* const ao_factor,
                              const int64_t size_max);




/* When all the data for the AOs have been provided, the following */
/* function returns ~true~. */


bool qmckl_ao_basis_provided (const qmckl_context context);

/* Access functions */


qmckl_exit_code
qmckl_get_ao_basis_primitive_vgl (qmckl_context context,
                                  double* const primitive_vgl,
                                  const int64_t size_max);

qmckl_exit_code
qmckl_get_ao_basis_shell_vgl (qmckl_context context,
                              double* const shell_vgl,
                              const int64_t size_max);

qmckl_exit_code
qmckl_get_ao_basis_ao_vgl (qmckl_context context,
                           double* const ao_vgl,
                           const int64_t size_max);



/* Uses the given array to compute the VGL. */


qmckl_exit_code
qmckl_get_ao_basis_ao_vgl_inplace (qmckl_context context,
                                   double* const ao_vgl,
                                   const int64_t size_max);

qmckl_exit_code
qmckl_get_ao_basis_ao_value (qmckl_context context,
                             double* const ao_value,
                             const int64_t size_max);



/* Uses the given array to compute the value. */


qmckl_exit_code
qmckl_get_ao_basis_ao_value_inplace (qmckl_context context,
                                     double* const ao_value,
                                     const int64_t size_max);

qmckl_exit_code
qmckl_ao_gaussian_vgl(const qmckl_context context,
                      const double *X,
                      const double *R,
                      const int64_t *n,
                      const int64_t *A,
                      const double *VGL,
                      const int64_t ldv);

qmckl_exit_code
qmckl_ao_slater_vgl(const qmckl_context context,
                    const double *X,
                    const double *R,
                    const int64_t *num_slater,
                    const int64_t *N,
                    const double *A,
                    const double *VGL,
                    const int64_t ldv);

qmckl_exit_code
qmckl_ao_slater_vgl_hpc(const qmckl_context context,
                         const double *X,
                         const double *R,
                         const int64_t *num_slater,
                         const int64_t *N,
                         const double *A,
                         double *VGL,
                         const int64_t ldv);

/* Computation of primitives */
/*    :PROPERTIES: */
/*    :Name:     qmckl_compute_ao_basis_primitive_gaussian_vgl */
/*    :CRetType: qmckl_exit_code */
/*    :FRetType: qmckl_exit_code */
/*    :END: */

/*    #+NAME: qmckl_ao_basis_primitive_gaussian_vgl_args */
/*    | Variable             | Type                             | In/Out | Description                                      | */
/*    |----------------------+----------------------------------+--------+--------------------------------------------------| */
/*    | ~context~            | ~qmckl_context~                  | in     | Global state                                     | */
/*    | ~prim_num~           | ~int64_t~                        | in     | Number of primitives                             | */
/*    | ~point_num~          | ~int64_t~                        | in     | Number of points considered                      | */
/*    | ~nucl_num~           | ~int64_t~                        | in     | Number of nuclei                                 | */
/*    | ~nucleus_prim_index~ | ~int64_t[nucl_num+1]~            | in     | Index of the 1st primitive of each nucleus       | */
/*    | ~coord~              | ~double[3][point_num]~           | in     | Coordinates                                      | */
/*    | ~nucl_coord~         | ~double[3][nucl_num]~            | in     | Nuclear coordinates                              | */
/*    | ~expo~               | ~double[prim_num]~               | in     | Exponents of the primitives                      | */
/*    | ~primitive_vgl~      | ~double[point_num][5][prim_num]~ | out    | Value, gradients and Laplacian of the primitives | */

/*    #+CALL: generate_private_c_header(table=qmckl_ao_basis_primitive_gaussian_vgl_args,rettyp=get_value("CRetType"),fname="qmckl_compute_ao_basis_primitive_gaussian_vgl")) */

/*    #+RESULTS: */

qmckl_exit_code qmckl_compute_ao_basis_primitive_gaussian_vgl (
      const qmckl_context context,
      const int64_t prim_num,
      const int64_t point_num,
      const int64_t nucl_num,
      const int64_t* nucleus_prim_index,
      const double* coord,
      const double* nucl_coord,
      const double* expo,
      double* const primitive_vgl );

/* Gaussians */
/*    :PROPERTIES: */
/*    :Name:     qmckl_compute_ao_basis_shell_gaussian_vgl */
/*    :CRetType: qmckl_exit_code */
/*    :FRetType: qmckl_exit_code */
/*    :END: */

/*    #+NAME: qmckl_ao_basis_shell_gaussian_vgl_args */
/*    | Variable            | Type                              | In/Out | Description                                  | */
/*    |---------------------+-----------------------------------+--------+----------------------------------------------| */
/*    | ~context~           | ~qmckl_context~                   | in     | Global state                                 | */
/*    | ~prim_num~          | ~int64_t~                         | in     | Number of primitives                         | */
/*    | ~shell_num~         | ~int64_t~                         | in     | Number of shells                             | */
/*    | ~point_num~         | ~int64_t~                         | in     | Number of points                             | */
/*    | ~nucl_num~          | ~int64_t~                         | in     | Number of nuclei                             | */
/*    | ~nucleus_shell_num~ | ~int64_t[nucl_num]~               | in     | Number of shells for each nucleus            | */
/*    | ~nucleus_index~     | ~int64_t[nucl_num]~               | in     | Index of the 1st shell of each nucleus       | */
/*    | ~nucleus_range~     | ~double[nucl_num]~                | in     | Range of the nucleus                         | */
/*    | ~shell_prim_index~  | ~int64_t[shell_num]~              | in     | Index of the 1st primitive of each shell     | */
/*    | ~shell_prim_num~    | ~int64_t[shell_num]~              | in     | Number of primitives per shell               | */
/*    | ~coord~             | ~double[3][point_num]~            | in     | Coordinates                                  | */
/*    | ~nucl_coord~        | ~double[3][nucl_num]~             | in     | Nuclear coordinates                          | */
/*    | ~expo~              | ~double[prim_num]~                | in     | Exponents of the primitives                  | */
/*    | ~coef_normalized~   | ~double[prim_num]~                | in     | Coefficients of the primitives               | */
/*    | ~shell_vgl~         | ~double[point_num][5][shell_num]~ | out    | Value, gradients and Laplacian of the shells | */

/*    #+CALL: generate_private_c_header(table=qmckl_ao_basis_shell_gaussian_vgl_args,rettyp=get_value("CRetType"),fname="qmckl_compute_ao_basis_shell_gaussian_vgl")) */

/*    #+RESULTS: */

qmckl_exit_code qmckl_compute_ao_basis_shell_gaussian_vgl (
      const qmckl_context context,
      const int64_t prim_num,
      const int64_t shell_num,
      const int64_t point_num,
      const int64_t nucl_num,
      const int64_t* nucleus_shell_num,
      const int64_t* nucleus_index,
      const double* nucleus_range,
      const int64_t* shell_prim_index,
      const int64_t* shell_prim_num,
      const double* coord,
      const double* nucl_coord,
      const double* expo,
      const double* coef_normalized,
      double* const shell_vgl );

/* Slater */
/*    :PROPERTIES: */
/*    :Name:     qmckl_compute_ao_basis_shell_slater_vgl */
/*    :CRetType: qmckl_exit_code */
/*    :FRetType: qmckl_exit_code */
/*    :END: */

/*    #+NAME: qmckl_ao_basis_shell_slater_vgl_args */
/*    | Variable            | Type                              | In/Out | Description                                  | */
/*    |---------------------+-----------------------------------+--------+----------------------------------------------| */
/*    | ~context~           | ~qmckl_context~                   | in     | Global state                                 | */
/*    | ~prim_num~          | ~int64_t~                         | in     | Number of primitives                         | */
/*    | ~shell_num~         | ~int64_t~                         | in     | Number of shells                             | */
/*    | ~point_num~         | ~int64_t~                         | in     | Number of points                             | */
/*    | ~nucl_num~          | ~int64_t~                         | in     | Number of nuclei                             | */
/*    | ~nucleus_shell_num~ | ~int64_t[nucl_num]~               | in     | Number of shells per nucleus                 | */
/*    | ~nucleus_index~     | ~int64_t[nucl_num]~               | in     | Index of the first shell of each nucleus     | */
/*    | ~nucleus_range~     | ~double[nucl_num]~                | in     | Range beyond which all shells are zero       | */
/*    | ~shell_prim_index~  | ~int64_t[shell_num]~              | in     | Index of the first primitive of each shell   | */
/*    | ~shell_prim_num~    | ~int64_t[shell_num]~              | in     | Number of primitives per shell               | */
/*    | ~r_power~           | ~int32_t[shell_num]~              | in     | Power of r prefix                            | */
/*    | ~coord~             | ~double[3][point_num]~            | in     | Coordinates of the points                    | */
/*    | ~nucl_coord~        | ~double[3][nucl_num]~             | in     | Nuclear coordinates                          | */
/*    | ~expo~              | ~double[prim_num]~                | in     | Exponents of the primitives                  | */
/*    | ~coef_normalized~   | ~double[prim_num]~                | in     | Coefficients of the primitives               | */
/*    | ~shell_vgl~         | ~double[point_num][5][shell_num]~ | out    | Value, gradients and Laplacian of the shells | */

/*    #+CALL: generate_private_c_header(table=qmckl_ao_basis_shell_slater_vgl_args,rettyp=get_value("CRetType"),fname="qmckl_compute_ao_basis_shell_slater_vgl")) */

/*    #+RESULTS: */

qmckl_exit_code qmckl_compute_ao_basis_shell_slater_vgl (
      const qmckl_context context,
      const int64_t prim_num,
      const int64_t shell_num,
      const int64_t point_num,
      const int64_t nucl_num,
      const int64_t* nucleus_shell_num,
      const int64_t* nucleus_index,
      const double* nucleus_range,
      const int64_t* shell_prim_index,
      const int64_t* shell_prim_num,
      const int32_t* r_power,
      const double* coord,
      const double* nucl_coord,
      const double* expo,
      const double* coef_normalized,
      double* const shell_vgl );

/* Get */


qmckl_exit_code
qmckl_get_ao_basis_shell_hessian (qmckl_context context,
                              double* const shell_hessian,
                              const int64_t size_max);

/* Computation of shell Hessians with Gaussian functions */
/*    :PROPERTIES: */
/*    :Name:     qmckl_compute_ao_basis_shell_gaussian_hessian */
/*    :CRetType: qmckl_exit_code */
/*    :FRetType: qmckl_exit_code */
/*    :END: */

/*    #+NAME: qmckl_ao_basis_shell_gaussian_hessian_args */
/*    | Variable            | Type                                 | In/Out | Description                              | */
/*    |---------------------+--------------------------------------+--------+------------------------------------------| */
/*    | ~context~           | ~qmckl_context~                      | in     | Global state                             | */
/*    | ~prim_num~          | ~int64_t~                            | in     | Number of primitives                     | */
/*    | ~shell_num~         | ~int64_t~                            | in     | Number of shells                         | */
/*    | ~point_num~         | ~int64_t~                            | in     | Number of points                         | */
/*    | ~nucl_num~          | ~int64_t~                            | in     | Number of nuclei                         | */
/*    | ~nucleus_shell_num~ | ~int64_t[nucl_num]~                  | in     | Number of shells for each nucleus        | */
/*    | ~nucleus_index~     | ~int64_t[nucl_num]~                  | in     | Index of the 1st shell of each nucleus   | */
/*    | ~nucleus_range~     | ~double[nucl_num]~                   | in     | Range of the nucleus                     | */
/*    | ~shell_prim_index~  | ~int64_t[shell_num]~                 | in     | Index of the 1st primitive of each shell | */
/*    | ~shell_prim_num~    | ~int64_t[shell_num]~                 | in     | Number of primitives per shell           | */
/*    | ~coord~             | ~double[3][point_num]~               | in     | Coordinates                              | */
/*    | ~nucl_coord~        | ~double[3][nucl_num]~                | in     | Nuclear coordinates                      | */
/*    | ~expo~              | ~double[prim_num]~                   | in     | Exponents of the primitives              | */
/*    | ~coef_normalized~   | ~double[prim_num]~                   | in     | Coefficients of the primitives           | */
/*    | ~shell_hessian~     | ~double[point_num][3][4][shell_num]~ | out    | Hessian of the shells                    | */
/*    |---------------------+--------------------------------------+--------+------------------------------------------| */


qmckl_exit_code qmckl_compute_ao_basis_shell_gaussian_hessian (
      const qmckl_context context,
      const int64_t prim_num,
      const int64_t shell_num,
      const int64_t point_num,
      const int64_t nucl_num,
      const int64_t* nucleus_shell_num,
      const int64_t* nucleus_index,
      const double* nucleus_range,
      const int64_t* shell_prim_index,
      const int64_t* shell_prim_num,
      const double* coord,
      const double* nucl_coord,
      const double* expo,
      const double* coef_normalized,
      double* const shell_hessian );

/* TODO Computation of shell Hessians with Slater functions */
/*    :PROPERTIES: */
/*    :Name:     qmckl_compute_ao_basis_shell_slater_hessian */
/*    :CRetType: qmckl_exit_code */
/*    :FRetType: qmckl_exit_code */
/*    :END: */

/*    #+NAME: qmckl_ao_basis_shell_slater_hessian_args */
/*    | Variable            | Type                                 | In/Out | Description                                | */
/*    |---------------------+--------------------------------------+--------+--------------------------------------------| */
/*    | ~context~           | ~qmckl_context~                      | in     | Global state                               | */
/*    | ~prim_num~          | ~int64_t~                            | in     | Number of primitives                       | */
/*    | ~shell_num~         | ~int64_t~                            | in     | Number of shells                           | */
/*    | ~point_num~         | ~int64_t~                            | in     | Number of points                           | */
/*    | ~nucl_num~          | ~int64_t~                            | in     | Number of nuclei                           | */
/*    | ~nucleus_shell_num~ | ~int64_t[nucl_num]~                  | in     | Number of shells per nucleus               | */
/*    | ~nucleus_index~     | ~int64_t[nucl_num]~                  | in     | Index of the first shell of each nucleus   | */
/*    | ~nucleus_range~     | ~double[nucl_num]~                   | in     | Range beyond which all shells are zero     | */
/*    | ~shell_prim_index~  | ~int64_t[shell_num]~                 | in     | Index of the first primitive of each shell | */
/*    | ~shell_prim_num~    | ~int64_t[shell_num]~                 | in     | Number of primitives per shell             | */
/*    | ~r_power~           | ~int32_t[shell_num]~                 | in     | Power of the r prefix                      | */
/*    | ~coord~             | ~double[3][point_num]~               | in     | Coordinates of the points                  | */
/*    | ~nucl_coord~        | ~double[3][nucl_num]~                | in     | Nuclear coordinates                        | */
/*    | ~expo~              | ~double[prim_num]~                   | in     | Exponents of the primitives                | */
/*    | ~coef_normalized~   | ~double[prim_num]~                   | in     | Coefficients of the primitives             | */
/*    | ~shell_hessian~     | ~double[point_num][3][4][shell_num]~ | out    | Hessian of the shells                      | */
/*    |---------------------+--------------------------------------+--------+--------------------------------------------| */


qmckl_exit_code qmckl_compute_ao_basis_shell_slater_hessian (
      const qmckl_context context,
      const int64_t prim_num,
      const int64_t shell_num,
      const int64_t point_num,
      const int64_t nucl_num,
      const int64_t* nucleus_shell_num,
      const int64_t* nucleus_index,
      const double* nucleus_range,
      const int64_t* shell_prim_index,
      const int64_t* shell_prim_num,
      const int32_t* r_power,
      const double* coord,
      const double* nucl_coord,
      const double* expo,
      const double* coef_normalized,
      double* const shell_hessian );

/* General functions for Powers of $x-X_i$ */
/*    :PROPERTIES: */
/*    :Name:     qmckl_ao_power */
/*    :CRetType: qmckl_exit_code */
/*    :FRetType: qmckl_exit_code */
/*    :END: */

/*    The ~qmckl_ao_power~ function computes all the powers of the ~n~ */
/*    input data up to the given maximum value given in input for each of */
/*    the $n$ points: */

/*    \[ P_{ik} = X_i^k \] */

/*    #+NAME: qmckl_ao_power_args */
/*    | Variable  | Type            | In/Out | Description                                       | */
/*    |-----------+-----------------+--------+---------------------------------------------------| */
/*    | ~context~ | ~qmckl_context~ | in     | Global state                                      | */
/*    | ~n~       | int64_t         | in     | Number of values                                  | */
/*    | ~X~       | double[n]       | in     | Array containing the input values                 | */
/*    | ~LMAX~    | int32_t[n]      | in     | Array containing the maximum power for each value | */
/*    | ~P~       | double[n][ldp]  | out    | Array containing all the powers of ~X~            | */
/*    | ~ldp~     | int64_t         | in     | Leading dimension of array ~P~                    | */

/*    Requirements: */

/*    - ~context~ is not ~QMCKL_NULL_CONTEXT~ */
/*    - ~n~ > 0 */
/*    - ~X~ is allocated with at least $n \times 8$ bytes */
/*    - ~LMAX~ is allocated with at least $n \times 4$ bytes */
/*    - ~P~ is allocated with at least $n \times \max_i \text{LMAX}_i \times 8$ bytes */
/*    - ~LDP~ >= $\max_i$ ~LMAX[i]~ */

/*    #+CALL: generate_c_header(table=qmckl_ao_power_args,rettyp=get_value("CRetType"),fname="qmckl_ao_power") */

/*    #+RESULTS: */

qmckl_exit_code qmckl_ao_power (
      const qmckl_context context,
      const int64_t n,
      const double* X,
      const int32_t* LMAX,
      double* const P,
      const int64_t ldp );

/* General functions for Value, Gradient and Laplacian of a polynomial */
/*    :PROPERTIES: */
/*    :Name:     qmckl_ao_polynomial_vgl */
/*    :CRetType: qmckl_exit_code */
/*    :FRetType: qmckl_exit_code */
/*    :END: */

/*    A polynomial is centered on a nucleus $\mathbf{R}_i$ */

/*    \[ */
/*    P_l(\mathbf{r},\mathbf{R}_i)  =   (x-X_i)^a (y-Y_i)^b (z-Z_i)^c */
/*    \] */

/*    The gradients with respect to electron coordinates are */

/*    \begin{eqnarray*} */
/*    \frac{\partial }{\partial x} P_l\left(\mathbf{r},\mathbf{R}_i \right) & */
/*                   = & a (x-X_i)^{a-1} (y-Y_i)^b (z-Z_i)^c \\ */
/*    \frac{\partial }{\partial y} P_l\left(\mathbf{r},\mathbf{R}_i \right) & */
/*                   = & b (x-X_i)^a (y-Y_i)^{b-1} (z-Z_i)^c \\ */
/*    \frac{\partial }{\partial z} P_l\left(\mathbf{r},\mathbf{R}_i \right) & */
/*                   = & c (x-X_i)^a (y-Y_i)^b (z-Z_i)^{c-1} \\ */
/*    \end{eqnarray*} */

/*    and the Laplacian is */

/*    \begin{eqnarray*} */
/*    \left( \frac{\partial }{\partial x^2} + */
/*               \frac{\partial }{\partial y^2} + */
/*               \frac{\partial }{\partial z^2} \right) P_l */
/*               \left(\mathbf{r},\mathbf{R}_i \right) &  = & */
/*             a(a-1) (x-X_i)^{a-2} (y-Y_i)^b (z-Z_i)^c + \\ */
/*          && b(b-1) (x-X_i)^a (y-Y_i)^{b-2} (z-Z_i)^c + \\ */
/*          && c(c-1) (x-X_i)^a (y-Y_i)^b (z-Z_i)^{c-2}. */
/*    \end{eqnarray*} */

/*    ~qmckl_ao_polynomial_vgl~ computes the values, gradients and */
/*    Laplacians at a given point in space, of all polynomials with an */
/*    angular momentum up to ~lmax~. */

/*    #+NAME: qmckl_ao_polynomial_vgl_args */
/*    | Variable  | Type              | In/Out | Description                                          | */
/*    |-----------+-------------------+--------+------------------------------------------------------| */
/*    | ~context~ | ~qmckl_context~   | in     | Global state                                         | */
/*    | ~X~       | ~double[3]~       | in     | Array containing the coordinates of the points       | */
/*    | ~R~       | ~double[3]~       | in     | Array containing the x,y,z coordinates of the center | */
/*    | ~lmax~    | ~int32_t~         | in     | Maximum angular momentum                             | */
/*    | ~n~       | ~int64_t~         | inout  | Number of computed polynomials                       | */
/*    | ~L~       | ~int32_t[n][ldl]~ | out    | Contains a,b,c for all ~n~ results                   | */
/*    | ~ldl~     | ~int64_t~         | in     | Leading dimension of ~L~                             | */
/*    | ~VGL~     | ~double[n][ldv]~  | out    | Value, gradients and Laplacian of the polynomials    | */
/*    | ~ldv~     | ~int64_t~         | in     | Leading dimension of array ~VGL~                     | */
/*    |-----------+-------------------+--------+------------------------------------------------------| */

/*    Requirements: */

/*     - ~context~ \ne ~QMCKL_NULL_CONTEXT~ */
/*     - ~n~ > 0 */
/*     - ~lmax~ >= 0 */
/*     - ~ldl~ >= 3 */
/*     - ~ldv~ >= 5 */
/*     - ~X~ is allocated with at least $3 \times 8$ bytes */
/*     - ~R~ is allocated with at least $3 \times 8$ bytes */
/*     - ~n~ >= ~(lmax+1)(lmax+2)(lmax+3)/6~ */
/*     - ~L~ is allocated with at least $3 \times n \times 4$ bytes */
/*     - ~VGL~ is allocated with at least $5 \times n \times 8$ bytes */
/*     - On output, ~n~ should be equal to ~(lmax+1)(lmax+2)(lmax+3)/6~ */
/*     - On output, the powers are given in the following order (l=a+b+c): */
/*       - Increasing values of ~l~ */
/*       - Within a given value of ~l~, alphabetical order of the */
/*         string made by a*"x" + b*"y" + c*"z" (in Python notation). */
/*         For example, with a=0, b=2 and c=1 the string is "yyz" */

/*     #+CALL: generate_c_header(table=qmckl_ao_polynomial_vgl_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */

/*     #+RESULTS: */

qmckl_exit_code qmckl_ao_polynomial_vgl (
      const qmckl_context context,
      const double* X,
      const double* R,
      const int32_t lmax,
      int64_t* n,
      int32_t* const L,
      const int64_t ldl,
      double* const VGL,
      const int64_t ldv );



/* #+CALL: generate_c_header(table=qmckl_ao_polynomial_vgl_args,rettyp=get_value("CRetType"),fname="qmckl_ao_polynomial_vgl_doc") */

/* #+RESULTS: */

qmckl_exit_code qmckl_ao_polynomial_vgl_doc (
      const qmckl_context context,
      const double* X,
      const double* R,
      const int32_t lmax,
      int64_t* n,
      int32_t* const L,
      const int64_t ldl,
      double* const VGL,
      const int64_t ldv );



/* #+CALL: generate_c_header(table=qmckl_ao_polynomial_vgl_args,rettyp=get_value("CRetType"),fname="qmckl_ao_polynomial_transp_vgl") */

/* #+RESULTS: */

qmckl_exit_code qmckl_ao_polynomial_transp_vgl (
      const qmckl_context context,
      const double* X,
      const double* R,
      const int32_t lmax,
      int64_t* n,
      int32_t* const L,
      const int64_t ldl,
      double* const VGL,
      const int64_t ldv );



/* #+CALL: generate_c_header(table=qmckl_ao_polynomial_vgl_args,rettyp=get_value("CRetType"),fname="qmckl_ao_polynomial_transp_vgl_doc") */

/* #+RESULTS: */

qmckl_exit_code qmckl_ao_polynomial_transp_vgl_doc (
      const qmckl_context context,
      const double* X,
      const double* R,
      const int32_t lmax,
      int64_t* n,
      int32_t* const L,
      const int64_t ldl,
      double* const VGL,
      const int64_t ldv );



/* #   #+CALL: generate_c_header(table=qmckl_ao_polynomial_vgl_args,rettyp=get_value("CRetType"),fname="qmckl_ao_polynomial_transp_vgl_hpc") */

/*     #+RESULTS: */

qmckl_exit_code qmckl_ao_polynomial_transp_vgl_hpc (
      const qmckl_context context,
      const double* X,
      const double* R,
      const int32_t lmax,
      int64_t* n,
      int32_t* const L,
      const int64_t ldl,
      double* const VGL,
      const int64_t ldv );

/* Hessian */

/* Compute the Hessian of the polynomial part of the atomic orbitals. Similarly to the shells, the ~hessian[:][3][:]~ component is reserved for the derivative of the Laplacian. */



/*   :PROPERTIES: */
/*    :Name:     qmckl_compute_ao_polynomial_hessian */
/*    :CRetType: qmckl_exit_code */
/*    :FRetType: qmckl_exit_code */
/*    :END: */

/*    #+NAME: qmckl_compute_ao_polynomial_hessian_args */
/*    | Variable  | Type                | In/Out | Description                                          | */
/*    |-----------+---------------------+--------+------------------------------------------------------| */
/*    | ~context~ | ~qmckl_context~     | in     | Global state                                         | */
/*    | ~X~       | ~double[3]~         | in     | Array containing the coordinates of the points       | */
/*    | ~R~       | ~double[3]~         | in     | Array containing the x,y,z coordinates of the center | */
/*    | ~lmax~    | ~int32_t~           | in     | Maximum angular momentum                             | */
/*    | ~n~       | ~int64_t~           | inout  | Number of computed polynomials                       | */
/*    | ~hessian~ | ~double[ldv][4][3]~ | out    | Hessian of the polynomials                           | */
/*    |-----------+---------------------+--------+------------------------------------------------------| */


qmckl_exit_code qmckl_compute_ao_polynomial_hessian_doc (
      const qmckl_context context,
      const double* X,
      const double* R,
      const int32_t lmax,
      int64_t* n,
      double* const hessian);

qmckl_exit_code qmckl_ao_polynomial_hessian (
      const qmckl_context context,
      const double* X,
      const double* R,
      const int32_t lmax,
      int64_t* n,
      double* const hessian );

qmckl_exit_code qmckl_compute_ao_polynomial_hessian_doc (
      const qmckl_context context,
      const double* X,
      const double* R,
      const int32_t lmax,
      int64_t* n,
      double* const hessian);

/* Get */


qmckl_exit_code
qmckl_get_ao_basis_ao_hessian(qmckl_context context,
                          double* const ao_hessian,
                          const int64_t size_max);
/* /home/runner/work/qmckl/qmckl/_build/src/qmckl_mo_func.h */
/* Initialization functions */

/*    To set the basis set, all the following functions need to be */
/*    called. */


qmckl_exit_code  qmckl_set_mo_basis_mo_num           (qmckl_context context, const int64_t   mo_num);
qmckl_exit_code  qmckl_set_mo_basis_coefficient      (qmckl_context context, const double  * coefficient, const int64_t size_max);
qmckl_exit_code  qmckl_set_mo_basis_r_cusp           (qmckl_context context, const double  * r_cusp, const int64_t size_max);

/* Access functions */


qmckl_exit_code
qmckl_get_mo_basis_mo_num (const qmckl_context context,
                           int64_t* mo_num);

qmckl_exit_code
qmckl_get_mo_basis_coefficient (const qmckl_context context,
                                double* const coefficient,
                                const int64_t size_max);



/* When all the data for the AOs have been provided, the following */
/* function returns ~true~. */


bool qmckl_mo_basis_provided (const qmckl_context context);

/* Update */

/*    It may be desirable to remove certain molecular orbitals (MOs) that */
/*    do not significantly contribute to the wave function.  In */
/*    particular, in a single determinant calculation, the virtual MOs */
/*    can be removed as they do not participate in the ground state */
/*    configuration. */

/*    To select a subset of MOs that will be kept, an array of integers of */
/*    size ~mo_num~ can be created. If the integer corresponding to an MO is */
/*    zero, that MO is dropped and will not be included in the */
/*    calculation. If the integer is non-zero, the MO will be kept. */



qmckl_exit_code
qmckl_mo_basis_select_mo (const qmckl_context context,
                          const int32_t* keep,
                          const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_mo_basis_mo_value(qmckl_context context,
                            double* const mo_value,
                            const int64_t size_max);



/* Uses the given array to compute the values. */


qmckl_exit_code
qmckl_get_mo_basis_mo_value_inplace (qmckl_context context,
                                     double* const mo_value,
                                     const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_mo_basis_mo_vgl(qmckl_context context,
                          double* const mo_vgl,
                          const int64_t size_max);



/* Uses the given array to compute the VGL. */


qmckl_exit_code
qmckl_get_mo_basis_mo_vgl_inplace (qmckl_context context,
                                   double* const mo_vgl,
                                   const int64_t size_max);



/* #+CALL: generate_private_c_header(table=qmckl_mo_basis_mo_value_cusp_args,rettyp=get_value("CRetType"),fname="qmckl_compute_mo_basis_mo_value_cusp_doc")) */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_mo_basis_mo_value_cusp_doc (
      const qmckl_context context,
      const int64_t nucl_num,
      const int64_t ao_num,
      const int64_t mo_num,
      const int64_t point_num,
      const int64_t* ao_nucl,
      const int32_t* ao_ang_mom,
      const double* en_distance,
      const double* r_cusp,
      const double* cusp_param,
      const double* coefficient_t,
      const double* ao_value,
      double* const mo_value );



/* #+CALL: generate_private_c_header(table=qmckl_mo_basis_mo_vgl_cusp_args,rettyp=get_value("CRetType"),fname="qmckl_compute_mo_basis_mo_vgl_cusp_doc")) */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_mo_basis_mo_vgl_cusp_doc (
      const qmckl_context context,
      const int64_t nucl_num,
      const int64_t ao_num,
      const int64_t mo_num,
      const int64_t point_num,
      const int64_t* ao_nucl,
      const int32_t* ao_ang_mom,
      const double* en_distance,
      const double* nucl_coord,
      const double* point_coord,
      const double* r_cusp,
      const double* cusp_param,
      const double* coefficient_t,
      const double* ao_vgl,
      double* const mo_vgl );

/* Rescaling of MO coefficients */

/*    When evaluating Slater determinants, the value of the determinants */
/*    may get out of the range of double precision. A simple fix is to */
/*    rescale the MO coefficients to put back the determinants in the */
/*    correct range. */


qmckl_exit_code
qmckl_mo_basis_rescale(qmckl_context context,
                          const double scaling_factor);
/* /home/runner/work/qmckl/qmckl/_build/src/qmckl_determinant_func.h */


/* When all the data for the slater determinants have been provided, the following */
/* function returns ~true~. */


bool      qmckl_determinant_provided           (const qmckl_context context);

/* Initialization functions */

/*    To set the basis set, all the following functions need to be */
/*    called. */


qmckl_exit_code  qmckl_set_determinant_type             (const qmckl_context context, const char t);
qmckl_exit_code  qmckl_set_determinant_det_num_alpha    (const qmckl_context context, const int64_t det_num_alpha);
qmckl_exit_code  qmckl_set_determinant_det_num_beta     (const qmckl_context context, const int64_t det_num_beta);
qmckl_exit_code  qmckl_set_determinant_mo_index_alpha   (const qmckl_context context, const int64_t* mo_index_alpha, const int64_t size_max);
qmckl_exit_code  qmckl_set_determinant_mo_index_beta    (const qmckl_context context, const int64_t* mo_index_beta, const int64_t size_max);

/* Get */


qmckl_exit_code qmckl_get_det_vgl_alpha(qmckl_context context, double* const det_vgl_alpha);
qmckl_exit_code qmckl_get_det_vgl_beta(qmckl_context context, double* const det_vgl_beta);



/*  #+CALL: generate_private_c_header(table=qmckl_compute_det_vgl_alpha_args,rettyp=get_value("CRetType"),fname="qmckl_compute_det_vgl_alpha")) */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_det_vgl_alpha (
      const qmckl_context context,
      const int64_t det_num_alpha,
      const int64_t walk_num,
      const int64_t alpha_num,
      const int64_t beta_num,
      const int64_t elec_num,
      const int64_t* mo_index_alpha,
      const int64_t mo_num,
      const double* mo_vgl,
      double* const det_vgl_alpha );



/*  #+CALL: generate_private_c_header(table=qmckl_compute_det_vgl_beta_args,rettyp=get_value("CRetType"),fname="qmckl_compute_det_vgl_beta")) */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_det_vgl_beta (
      const qmckl_context context,
      const int64_t det_num_beta,
      const int64_t walk_num,
      const int64_t alpha_num,
      const int64_t beta_num,
      const int64_t elec_num,
      const int64_t* mo_index_beta,
      const int64_t mo_num,
      const double* mo_vgl,
      double* const det_vgl_beta );

/* Get */


qmckl_exit_code qmckl_get_det_inv_matrix_alpha(qmckl_context context, double* const det_inv_matrix_alpha);
qmckl_exit_code qmckl_get_det_inv_matrix_beta(qmckl_context context, double* const det_inv_matrix_beta);
qmckl_exit_code qmckl_get_det_adj_matrix_alpha(qmckl_context context, double* const det_adj_matrix_alpha);
qmckl_exit_code qmckl_get_det_adj_matrix_beta(qmckl_context context, double* const det_adj_matrix_beta);
qmckl_exit_code qmckl_get_det_alpha(qmckl_context context, double* const det_adj_matrix_alpha);
qmckl_exit_code qmckl_get_det_beta(qmckl_context context, double* const det_adj_matrix_beta);



/* #+CALL: generate_private_c_header(table=qmckl_det_inv_matrix_alpha_args,rettyp=get_value("CRetType"),fname="qmckl_compute_det_inv_matrix_alpha")) */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_det_inv_matrix_alpha (
      const qmckl_context context,
      const int64_t det_num_alpha,
      const int64_t walk_num,
      const int64_t alpha_num,
      const double* det_vgl_alpha,
      double* const det_value_alpha,
      double* const det_adj_matrix_alpha,
      double* const det_inv_matrix_alpha );



/*  #+CALL: generate_private_c_header(table=qmckl_det_inv_matrix_beta_args,rettyp=get_value("CRetType"),fname="qmckl_compute_det_inv_matrix_beta")) */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_det_inv_matrix_beta (
      const qmckl_context context,
      const int64_t det_num_beta,
      const int64_t walk_num,
      const int64_t beta_num,
      const double* det_vgl_beta,
      double* const det_value_beta,
      double* const det_adj_matrix_beta,
      double* const det_inv_matrix_beta );
/* /home/runner/work/qmckl/qmckl/_build/src/qmckl_sherman_morrison_woodbury_func.h */
qmckl_exit_code qmckl_sm_naive (
      const qmckl_context context,
      const uint64_t LDS,
      const uint64_t Dim,
      const uint64_t N_updates,
      const double* Updates,
      const uint64_t* Updates_index,
      const double breakdown,
      double* Slater_inv,
      double* determinant );



/* #+CALL: generate_c_header(table=qmckl_sm_naive_args,rettyp=get_value("CRetType"),fname="qmckl_sm_naive_doc") */

/* #+RESULTS: */

qmckl_exit_code qmckl_sm_naive_doc (
      const qmckl_context context,
      const uint64_t LDS,
      const uint64_t Dim,
      const uint64_t N_updates,
      const double* Updates,
      const uint64_t* Updates_index,
      const double breakdown,
      double* Slater_inv,
      double* determinant );

/* C headers (exposed in qmckl.h) */
/* #+CALL: generate_c_header(table=qmckl_sm_splitting_core_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */

/* #+RESULTS: */

qmckl_exit_code qmckl_sm_splitting_core (
      const qmckl_context context,
      const uint64_t LDS,
      const uint64_t Dim,
      const uint64_t N_updates,
      const double* Updates,
      const uint64_t* Updates_index,
      const double breakdown,
      double* Slater_inv,
      double* later_updates,
      uint64_t* later_index,
      uint64_t* later,
      double* determinant );

qmckl_exit_code qmckl_sm_splitting_core_doc (
      const qmckl_context context,
      const uint64_t LDS,
      const uint64_t Dim,
      const uint64_t N_updates,
      const double* Updates,
      const uint64_t* Updates_index,
      const double breakdown,
      double* Slater_inv,
      double* later_updates,
      uint64_t* later_index,
      uint64_t* later,
      double* determinant );

/* C headers (exposed in qmckl.h) */
/* #+CALL: generate_c_header(table=qmckl_woodbury_2x2_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */

/* #+RESULTS: */

qmckl_exit_code qmckl_woodbury_2x2 (
      const qmckl_context context,
      const uint64_t LDS,
      const uint64_t Dim,
      const double* Updates,
      const uint64_t* Updates_index,
      const double breakdown,
      double* Slater_inv,
      double* determinant );



/* #+CALL: generate_c_header(table=qmckl_woodbury_2x2_args,rettyp=get_value("CRetType"),fname="qmckl_woodbury_2x2_hpc") */

/* #+RESULTS: */

qmckl_exit_code qmckl_woodbury_2x2_hpc (
      const qmckl_context context,
      const uint64_t LDS,
      const uint64_t Dim,
      const double* Updates,
      const uint64_t* Updates_index,
      const double breakdown,
      double* Slater_inv,
      double* determinant );



/* #+CALL: generate_c_header(table=qmckl_woodbury_2x2_args,rettyp=get_value("CRetType"),fname="qmckl_woodbury_2x2_doc") */

/* #+RESULTS: */

qmckl_exit_code qmckl_woodbury_2x2_doc (
      const qmckl_context context,
      const uint64_t LDS,
      const uint64_t Dim,
      const double* Updates,
      const uint64_t* Updates_index,
      const double breakdown,
      double* Slater_inv,
      double* determinant );

/* C headers (exposed in qmckl.h) */
/* #+CALL: generate_c_header(table=qmckl_woodbury_3x3_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */

/* #+RESULTS: */

qmckl_exit_code qmckl_woodbury_3x3 (
      const qmckl_context context,
      const uint64_t LDS,
      const uint64_t Dim,
      const double* Updates,
      const uint64_t* Updates_index,
      const double breakdown,
      double* Slater_inv,
      double* determinant );



/* #+CALL: generate_c_header(table=qmckl_woodbury_3x3_args,rettyp=get_value("CRetType"),fname="qmckl_woodbury_3x3_hpc") */

/* #+RESULTS: */

qmckl_exit_code qmckl_woodbury_3x3_hpc (
      const qmckl_context context,
      const uint64_t LDS,
      const uint64_t Dim,
      const double* Updates,
      const uint64_t* Updates_index,
      const double breakdown,
      double* Slater_inv,
      double* determinant );



/* #+CALL: generate_c_header(table=qmckl_woodbury_3x3_args,rettyp=get_value("CRetType"),fname="qmckl_woodbury_3x3_doc") */

/* #+RESULTS: */

qmckl_exit_code qmckl_woodbury_3x3_doc (
      const qmckl_context context,
      const uint64_t LDS,
      const uint64_t Dim,
      const double* Updates,
      const uint64_t* Updates_index,
      const double breakdown,
      double* Slater_inv,
      double* determinant );

/* C headers (exposed in qmckl.h) */

/* #+CALL: generate_c_header(table=qmckl_sm_splitting_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */

/* #+RESULTS: */

qmckl_exit_code qmckl_sm_splitting (
      const qmckl_context context,
      const uint64_t LDS,
      const uint64_t Dim,
      const uint64_t N_updates,
      const double* Updates,
      const uint64_t* Updates_index,
      const double breakdown,
      double* Slater_inv,
      double* determinant );



/* #+CALL: generate_c_header(table=qmckl_sm_splitting_args,rettyp=get_value("CRetType"),fname="qmckl_sm_splitting_doc") */

/* #+RESULTS: */

qmckl_exit_code qmckl_sm_splitting_doc (
      const qmckl_context context,
      const uint64_t LDS,
      const uint64_t Dim,
      const uint64_t N_updates,
      const double* Updates,
      const uint64_t* Updates_index,
      const double breakdown,
      double* Slater_inv,
      double* determinant );
/* /home/runner/work/qmckl/qmckl/_build/src/qmckl_jastrow_champ_func.h */
/* Initialization functions */

/*    To prepare for the Jastrow and its derivative, all the following functions need to be */
/*    called. */


qmckl_exit_code  qmckl_set_jastrow_champ_rescale_factor_ee (qmckl_context context, const double  kappa_ee);
qmckl_exit_code  qmckl_set_jastrow_champ_rescale_factor_en (qmckl_context context, const double* kappa_en, const int64_t size_max);
qmckl_exit_code  qmckl_set_jastrow_champ_aord_num          (qmckl_context context, const int64_t aord_num);
qmckl_exit_code  qmckl_set_jastrow_champ_bord_num          (qmckl_context context, const int64_t bord_num);
qmckl_exit_code  qmckl_set_jastrow_champ_cord_num          (qmckl_context context, const int64_t cord_num);
qmckl_exit_code  qmckl_set_jastrow_champ_type_nucl_num     (qmckl_context context, const int64_t type_nucl_num);
qmckl_exit_code  qmckl_set_jastrow_champ_type_nucl_vector  (qmckl_context context, const int64_t* type_nucl_vector, const int64_t size_max);
qmckl_exit_code  qmckl_set_jastrow_champ_a_vector          (qmckl_context context, const double * a_vector, const int64_t size_max);
qmckl_exit_code  qmckl_set_jastrow_champ_b_vector          (qmckl_context context, const double * b_vector, const int64_t size_max);
qmckl_exit_code  qmckl_set_jastrow_champ_c_vector          (qmckl_context context, const double * c_vector, const int64_t size_max);
qmckl_exit_code  qmckl_set_jastrow_champ_spin_independent  (qmckl_context context, const int32_t spin_independent);

/* Access functions */


qmckl_exit_code  qmckl_get_jastrow_champ_aord_num          (qmckl_context context, int64_t* const aord_num);
qmckl_exit_code  qmckl_get_jastrow_champ_bord_num          (qmckl_context context, int64_t* const bord_num);
qmckl_exit_code  qmckl_get_jastrow_champ_cord_num          (qmckl_context context, int64_t* const bord_num);
qmckl_exit_code  qmckl_get_jastrow_champ_type_nucl_num     (qmckl_context context, int64_t* const type_nucl_num);
qmckl_exit_code  qmckl_get_jastrow_champ_type_nucl_vector  (qmckl_context context, int64_t* const type_nucl_num, const int64_t size_max);
qmckl_exit_code  qmckl_get_jastrow_champ_a_vector          (qmckl_context context, double * const a_vector, const int64_t size_max);
qmckl_exit_code  qmckl_get_jastrow_champ_b_vector          (qmckl_context context, double * const b_vector, const int64_t size_max);
qmckl_exit_code  qmckl_get_jastrow_champ_c_vector          (qmckl_context context, double * const c_vector, const int64_t size_max);
qmckl_exit_code  qmckl_get_jastrow_champ_rescale_factor_ee (const qmckl_context context, double* const rescale_factor_ee);
qmckl_exit_code  qmckl_get_jastrow_champ_rescale_factor_en (const qmckl_context context, double* const rescale_factor_en, const int64_t size_max);
qmckl_exit_code  qmckl_get_jastrow_champ_dim_c_vector      (qmckl_context context, int64_t* const dim_c_vector);
qmckl_exit_code  qmckl_get_jastrow_champ_spin_independent  (qmckl_context context, int32_t* const spin_independent);




/* Along with these core functions, calculation of the jastrow factor */
/* requires the following additional information to be set: */


/* When all the data for the AOs have been provided, the following */
/* function returns ~true~. */


bool      qmckl_jastrow_champ_provided           (const qmckl_context context);

/* Get */

qmckl_exit_code
qmckl_get_jastrow_champ_asymp_jasb(qmckl_context context,
                             double* const asymp_jasb,
                             const int64_t size_max);

/* Get */


qmckl_exit_code qmckl_get_jastrow_champ_ee_distance_rescaled(qmckl_context context,
                                                             double* const distance_rescaled,
                                                             int64_t const size_max);

/* Get */


qmckl_exit_code qmckl_get_jastrow_champ_ee_distance_rescaled_gl(qmckl_context context,
                                                                double* const distance_rescaled_gl,
                                                                const int64_t size_max);

/* Get */

qmckl_exit_code
qmckl_get_jastrow_champ_factor_ee(qmckl_context context,
                            double* const factor_ee,
                            const int64_t size_max);

/* Get */

qmckl_exit_code
qmckl_get_jastrow_champ_factor_ee_gl(qmckl_context context,
                                    double* const factor_ee_gl,
                                    const int64_t size_max);

/* Get */

qmckl_exit_code
qmckl_get_jastrow_champ_asymp_jasb_pderiv(qmckl_context context,
                                          double* const asymp_jasb_pderiv,
                                          const int64_t size_max);

/* Get */

qmckl_exit_code
qmckl_get_jastrow_champ_factor_ee_pderiv(qmckl_context context,
                                         double* const factor_ee_pderiv,
                                         const int64_t size_max);

/* Get */

qmckl_exit_code
qmckl_get_jastrow_champ_factor_ee_gl_pderiv(qmckl_context context,
                                            double* const factor_ee_gl_pderiv,
                                            const int64_t size_max);

/* Get */

qmckl_exit_code
qmckl_get_jastrow_champ_asymp_jasa(qmckl_context context,
                                   double* const asymp_jasa,
                                   const int64_t size_max);



/* #+CALL: generate_private_c_header(table=qmckl_asymp_jasa_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_jastrow_champ_asymp_jasa (
      const qmckl_context context,
      const int64_t aord_num,
      const int64_t type_nucl_num,
      const double* a_vector,
      const double* rescale_factor_en,
      double* const asymp_jasa );

/* Get */


qmckl_exit_code
qmckl_get_jastrow_champ_en_distance_rescaled(qmckl_context context,
                                             double* const distance_rescaled,
                                             const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_jastrow_champ_en_distance_rescaled_gl(qmckl_context context,
                                                double* const distance_rescaled_gl,
                                                const int64_t size_max);

/* Get */

qmckl_exit_code
qmckl_get_jastrow_champ_factor_en(qmckl_context context,
                            double* const factor_en,
                            const int64_t size_max);



/* #+CALL: generate_private_c_header(table=qmckl_factor_en_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_jastrow_champ_factor_en (
         const qmckl_context context,
         const int64_t walk_num,
         const int64_t elec_num,
         const int64_t nucl_num,
         const int64_t type_nucl_num,
         const int64_t* type_nucl_vector,
         const int64_t aord_num,
         const double* a_vector,
         const double* en_distance_rescaled,
         const double* asymp_jasa,
         double* const factor_en );

qmckl_exit_code qmckl_compute_jastrow_champ_factor_en_doc (
         const qmckl_context context,
         const int64_t walk_num,
         const int64_t elec_num,
         const int64_t nucl_num,
         const int64_t type_nucl_num,
         const int64_t* type_nucl_vector,
         const int64_t aord_num,
         const double* a_vector,
         const double* en_distance_rescaled,
         const double* asymp_jasa,
         double* const factor_en );

qmckl_exit_code qmckl_compute_jastrow_champ_factor_en_hpc (
         const qmckl_context context,
         const int64_t walk_num,
         const int64_t elec_num,
         const int64_t nucl_num,
         const int64_t type_nucl_num,
         const int64_t* type_nucl_vector,
         const int64_t aord_num,
         const double* a_vector,
         const double* en_distance_rescaled,
         const double* asymp_jasa,
         double* const factor_en );

/* Get */

qmckl_exit_code
qmckl_get_jastrow_champ_factor_en_gl(qmckl_context context,
                                    double* const factor_en_gl,
                                    const int64_t size_max);

/* Get */

qmckl_exit_code
qmckl_get_jastrow_champ_asymp_jasa_pderiv(qmckl_context context,
                                          double* const asymp_jasa_pderiv,
                                          const int64_t size_max);

qmckl_exit_code qmckl_compute_jastrow_champ_asymp_jasa_pderiv (
      const qmckl_context context,
      const int64_t aord_num,
      const int64_t type_nucl_num,
      const double* a_vector,
      const double* rescale_factor_en,
      double* const asymp_jasa_pderiv );

/* Get */

qmckl_exit_code
qmckl_get_jastrow_champ_factor_en_pderiv(qmckl_context context,
                                         double* const factor_en_pderiv,
                                         const int64_t size_max);

qmckl_exit_code qmckl_compute_jastrow_champ_factor_en_pderiv (
         const qmckl_context context,
         const int64_t walk_num,
         const int64_t elec_num,
         const int64_t nucl_num,
         const int64_t type_nucl_num,
         const int64_t* type_nucl_vector,
         const int64_t aord_num,
         const double* a_vector,
         const double* en_distance_rescaled,
         const double* asymp_jasa_pderiv,
         double* const factor_en_pderiv );

qmckl_exit_code qmckl_compute_jastrow_champ_factor_en_pderiv_doc (
         const qmckl_context context,
         const int64_t walk_num,
         const int64_t elec_num,
         const int64_t nucl_num,
         const int64_t type_nucl_num,
         const int64_t* type_nucl_vector,
         const int64_t aord_num,
         const double* a_vector,
         const double* en_distance_rescaled,
         const double* asymp_jasa_pderiv,
         double* const factor_en_pderiv );

/* Get */

qmckl_exit_code
qmckl_get_jastrow_champ_factor_en_gl_pderiv(qmckl_context context,
                                    double* const factor_en_gl_pderiv,
                                    const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_jastrow_champ_een_rescaled_e(qmckl_context context,
                                                double* const een_rescaled_e,
                                                const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_jastrow_champ_een_rescaled_e_gl(qmckl_context context,
                                          double* const een_rescaled_e_gl,
                                          const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_jastrow_champ_een_rescaled_n(qmckl_context context,
                                 double* const een_rescaled_n,
                                 const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_jastrow_champ_een_rescaled_n_gl(qmckl_context context,
                                         double* const een_rescaled_n_gl,
                                         const int64_t size_max);


/* #+CALL: generate_private_c_header(table=lkpm_combined_index_args,rettyp=get_value("CRetType"),fname="qmckl_compute_lkpm_combined_index_doc") */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_lkpm_combined_index_doc (
      const qmckl_context context,
      const int64_t cord_num,
      const int64_t dim_c_vector,
      int64_t* const lkpm_combined_index );



/* #+CALL: generate_private_c_header(table=lkpm_combined_index_args,rettyp=get_value("CRetType"),fname="qmckl_compute_lkpm_combined_index_hpc") */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_lkpm_combined_index_hpc (
      const qmckl_context context,
      const int64_t cord_num,
      const int64_t dim_c_vector,
      int64_t* const lkpm_combined_index );





/* #+CALL: generate_private_c_header(table=qmckl_factor_tmp_c_args,rettyp=get_value("CRetType"),fname="qmckl_compute_tmp_c") */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_tmp_c (
          const qmckl_context context,
          const int64_t cord_num,
          const int64_t elec_num,
          const int64_t nucl_num,
          const int64_t walk_num,
          const double* een_rescaled_e,
          const double* een_rescaled_n,
          double* const tmp_c );

/* Get */

qmckl_exit_code
qmckl_get_jastrow_champ_factor_een(qmckl_context context,
                             double* const factor_een,
                             const int64_t size_max);

/* Get */

qmckl_exit_code
qmckl_get_jastrow_champ_factor_een_gl(qmckl_context context,
                                     double* const factor_een_gl,
                                     const int64_t size_max);

qmckl_exit_code
qmckl_get_jastrow_champ_factor_een_grad(qmckl_context context,
                                        double* const factor_een_grad,
                                        const int64_t size_max);

/* Get */

qmckl_exit_code
qmckl_get_jastrow_champ_factor_een_pderiv(qmckl_context context,
                                          double* const factor_een_pderiv,
                                          const int64_t size_max);

/* Get */

qmckl_exit_code
qmckl_get_jastrow_champ_factor_een_gl_pderiv(qmckl_context context,
                                     double* const factor_een_gl_pderiv,
                                     const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_jastrow_champ_value(qmckl_context context,
                            double* const value,
                            const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_jastrow_champ_gl(qmckl_context context,
                            double* const gl,
                            const int64_t size_max);

qmckl_exit_code
qmckl_get_jastrow_champ_grad(qmckl_context context,
                            double* const grad,
                            const int64_t size_max);

qmckl_exit_code
qmckl_get_jastrow_champ_grad(qmckl_context context,
                            double* const grad,
                            const int64_t size_max);
/* /home/runner/work/qmckl/qmckl/_build/src/qmckl_jastrow_champ_single_func.h */
/* Set */

/*    We set the coordinates of the ~num~-th electron for all walkers, where ~num~ is the electron which has to be moved. */
/*    The dimension of ~coord~ is */
/*    - [walk_num][3] if ~transp~ is ~'N'~ */
/*    - [3][walk_num] if ~transp~ is ~'T'~ */

/*    Internally, the coordinates are stored in 'N' format as opposed to elec_coord. */
/*    This function has to be called before any other functions from this module. */


qmckl_exit_code qmckl_set_single_point (qmckl_context context,
                                 const char transp,
                                 const int64_t num,
                                 const double* coord,
                                 const int64_t size_max);

/* Touch */


qmckl_exit_code
qmckl_single_touch (const qmckl_context context);

/* Get */


qmckl_exit_code qmckl_get_single_electron_ee_distance(qmckl_context context,
                                                      double* const distance,
                                                      const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_single_electron_en_distance(qmckl_context context,
                                      double* distance,
                                      const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_een_rescaled_single_e(qmckl_context context,
                                 double* const distance_rescaled,
                                 const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_een_rescaled_single_n(qmckl_context context,
                                 double* const distance_rescaled,
                                 const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_jastrow_champ_single_een(qmckl_context context,
                                 double* const delta_een,
                                 const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_een_rescaled_single_n_gl(qmckl_context context,
                                         double* const distance_rescaled,
                                         const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_een_rescaled_single_e_gl(qmckl_context context,
                                         double* const distance_rescaled,
                                         const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_jastrow_champ_single_tmp(qmckl_context context,
                                 double* const tmp,
                                 const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_jastrow_champ_single_een_gl(qmckl_context context,
                                 double* const delta_een_gl,
                                 const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_jastrow_champ_single_een_g(qmckl_context context,
                                 double* const delta_een_g,
                                 const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_jastrow_champ_single_een_pderiv(qmckl_context context,
                                 double* const delta_een_pderiv,
                                 const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_ee_rescaled_single(qmckl_context context,
                                 double* const distance_rescaled,
                                 const int64_t size_max);

/* Get */

qmckl_exit_code
qmckl_get_jastrow_champ_single_ee(qmckl_context context,
                            double* const delta_ee,
                            const int64_t size_max);

/* Get */


qmckl_exit_code qmckl_get_ee_rescaled_single_gl(qmckl_context context,
                                                double* const distance_rescaled_gl,
                                                const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_jastrow_champ_single_ee_gl(qmckl_context context,
                                    double* const delta_ee_gl,
                                    const int64_t size_max);

/* Get */

qmckl_exit_code
qmckl_get_jastrow_champ_single_ee_pderiv(qmckl_context context,
                            double* const delta_ee_pderiv,
                            const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_en_rescaled_single(qmckl_context context,
                                 double* const distance_rescaled,
                                 const int64_t size_max);

/* Get */

qmckl_exit_code
qmckl_get_jastrow_champ_single_en(qmckl_context context,
                            double* const delta_en,
                            const int64_t size_max);

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
         double* const delta_en );

qmckl_exit_code qmckl_compute_jastrow_champ_single_en_doc (
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
         double* const delta_en );

/* Get */


qmckl_exit_code qmckl_get_en_rescaled_single_gl(qmckl_context context,
                                                double* distance_rescaled_gl,
                                                const size_t size_max);

/* Get */

qmckl_exit_code
qmckl_get_jastrow_champ_single_en_gl(qmckl_context context,
                                    double* const delta_en_gl,
                                    const int64_t size_max);

/* Get */

qmckl_exit_code
qmckl_get_jastrow_champ_single_en_pderiv(qmckl_context context,
                            double* const delta_en_pderiv,
                            const int64_t size_max);

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
         double* const delta_en_pderiv );

qmckl_exit_code qmckl_compute_jastrow_champ_single_en_pderiv_doc (
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
         double* const delta_en_pderiv );

/* Code */


qmckl_exit_code
qmckl_get_jastrow_champ_single_accept(qmckl_context context);
/* /home/runner/work/qmckl/qmckl/_build/src/qmckl_local_energy_func.h */
/* Access functions */

/*    When all the data for the local energy have been provided, the following */
/*    function returns ~true~. */


bool      qmckl_local_energy_provided           (const qmckl_context context);

/* Get */


qmckl_exit_code qmckl_get_kinetic_energy(qmckl_context context, double* const kinetic_energy);



/* #+CALL: generate_private_c_header(table=qmckl_compute_kinetic_energy_args,rettyp=get_value("CRetType"),fname="qmckl_compute_kinetic_energy")) */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_kinetic_energy (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t det_num_alpha,
      const int64_t det_num_beta,
      const int64_t alpha_num,
      const int64_t beta_num,
      const int64_t elec_num,
      const int64_t* mo_index_alpha,
      const int64_t* mo_index_beta,
      const int64_t mo_num,
      const double* mo_vgl,
      const double* det_value_alpha,
      const double* det_value_beta,
      const double* det_inv_matrix_alpha,
      const double* det_inv_matrix_beta,
      double* const e_kin );

/* Get */


qmckl_exit_code qmckl_get_potential_energy(qmckl_context context, double* const potential_energy);



/*  #+CALL: generate_private_c_header(table=qmckl_compute_potential_energy_args,rettyp=get_value("CRetType"),fname="qmckl_compute_potential_energy")) */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_potential_energy (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const double* ee_potential,
      const double* en_potential,
      const double repulsion,
      double* const e_pot );

/* Get */


qmckl_exit_code qmckl_get_local_energy(qmckl_context context, double* const local_energy, const int64_t size_max);



/*  #+CALL: generate_private_c_header(table=qmckl_compute_local_energy_args,rettyp=get_value("CRetType"),fname="qmckl_compute_local_energy")) */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_local_energy (
      const qmckl_context context,
      const int64_t walk_num,
      const double* e_kin,
      const double* e_pot,
      double* const e_local );

/* Get */


qmckl_exit_code qmckl_get_drift_vector(qmckl_context context, double* const drift_vector);



/*  #+CALL: generate_private_c_header(table=qmckl_compute_drift_vector_args,rettyp=get_value("CRetType"),fname="qmckl_compute_drift_vector")) */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_drift_vector (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t det_num_alpha,
      const int64_t det_num_beta,
      const int64_t alpha_num,
      const int64_t beta_num,
      const int64_t elec_num,
      const int64_t* mo_index_alpha,
      const int64_t* mo_index_beta,
      const int64_t mo_num,
      const double* mo_vgl,
      const double* det_inv_matrix_alpha,
      const double* det_inv_matrix_beta,
      double* const r_drift );
/* /home/runner/work/qmckl/qmckl/_build/src/qmckl_forces_func.h */
/* Get */



qmckl_exit_code
qmckl_get_forces_jastrow_en(qmckl_context context,
                                    double* const forces_jastrow_en,
                                    const int64_t size_max);

/* Get */



qmckl_exit_code
qmckl_get_forces_jastrow_en_g(qmckl_context context,
                                    double* const forces_jastrow_en_g,
                                    const int64_t size_max);

/* Get */



qmckl_exit_code
qmckl_get_forces_jastrow_en_l(qmckl_context context,
                                    double* const forces_jastrow_en_l,
                                    const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_forces_jastrow_single_en(qmckl_context context,
                          double* const forces_jastrow_single_en,
                          const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_forces_tmp_c(qmckl_context context,
                                    double* const forces_tmp_c,
                                    const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_forces_dtmp_c(qmckl_context context,
                                    double* const forces_dtmp_c,
                                    const int64_t size_max);

/* Get */



qmckl_exit_code
qmckl_get_forces_jastrow_een(qmckl_context context,
                                    double* const forces_jastrow_een,
                                    const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_forces_een_rescaled_n_gl(qmckl_context context,
                                         double* const forces_een_n,
                                         const int64_t size_max);

/* Get */

qmckl_exit_code
qmckl_get_forces_jastrow_een_g(qmckl_context context,
                                     double* const forces_jastrow_een_g,
                                     const int64_t size_max);

/* Get */

qmckl_exit_code
qmckl_get_forces_jastrow_een_l(qmckl_context context,
                                     double* const forces_jastrow_een_l,
                                     const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_forces_jastrow_single_een(qmckl_context context,
                          double* const forces_jastrow_single_een,
                          const int64_t size_max);

/* Get */

qmckl_exit_code
qmckl_get_forces_ao_value(qmckl_context context,
                                     double* const forces_ao_value,
                                     const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_forces_mo_value(qmckl_context context,
                          double* const forces_mo_value,
                          const int64_t size_max);

qmckl_exit_code
qmckl_get_forces_mo_value_inplace (qmckl_context context,
                                     double* const forces_mo_value,
                                     const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_forces_mo_g(qmckl_context context,
                          double* const forces_mo_g,
                          const int64_t size_max);

qmckl_exit_code
qmckl_get_forces_mo_g_inplace(qmckl_context context,
                              double* const forces_mo_g,
                              const int64_t size_max);

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
         double* const forces_mo_g );

qmckl_exit_code qmckl_compute_forces_mo_g_doc (
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
         double* const forces_mo_g );

qmckl_exit_code qmckl_compute_forces_mo_g_hpc (
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
         double* const forces_mo_g );

/* Get */


qmckl_exit_code
qmckl_get_forces_mo_l(qmckl_context context,
                          double* const forces_mo_l,
                          const int64_t size_max);



/* #+RESULTS: */

qmckl_exit_code qmckl_compute_forces_mo_l_doc (
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
      double* const forces_mo_l );
/* /home/runner/work/qmckl/qmckl/_build/src/qmckl_trexio_func.h */
qmckl_exit_code
qmckl_trexio_read(const qmckl_context context,
                  const char* file_name,
                  const int64_t size_max);
#ifdef __cplusplus
}
#endif
#endif
