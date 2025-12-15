/* [[file:../../org/qmckl_blas.org::*Headers][Headers:5]] */
#include "qmckl.h"
#include <stdio.h>
#include <assert.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "qmckl_memory_private_type.h"
#include "qmckl_blas_private_type.h"

#include "qmckl_memory_private_func.h"
#include "qmckl_blas_private_func.h"

int main() {
  qmckl_context context;
  context = qmckl_context_create();
/* Headers:5 ends here */

/* [[file:../../org/qmckl_blas.org::*Tests][Tests:1]] */
{
  int64_t m = 3;
  int64_t n = 4;
  int64_t p = m*n;
  qmckl_vector vec = qmckl_vector_alloc(context, p);

  for (int64_t i=0 ; i<p ; ++i)
    qmckl_vec(vec, i) = (double) i;

  for (int64_t i=0 ; i<p ; ++i)
    assert( vec.data[i] == (double) i );

  printf("qmckl_vector ok\n");

  qmckl_matrix mat = qmckl_matrix_of_vector(vec, m, n);
  assert (mat.size[0] == m);
  assert (mat.size[1] == n);
  assert (mat.data == vec.data);

  for (int64_t j=0 ; j<n ; ++j)
    for (int64_t i=0 ; i<m ; ++i)
      assert ( qmckl_mat(mat, i, j) == qmckl_vec(vec, i+j*m)) ;

  printf("qmckl_matrix_of_vector ok\n");

  qmckl_vector vec2 = qmckl_vector_of_matrix(mat);
  assert (vec2.size == p);
  assert (vec2.data == vec.data);
  for (int64_t i=0 ; i<p ; ++i)
      assert ( qmckl_vec(vec2, i) == qmckl_vec(vec, i) ) ;

  printf("qmckl_vector_of_matrix ok\n");

  double* dbl = qmckl_alloc_double_of_matrix(context, mat);
  for (int64_t i=0 ; i<p ; ++i)
      assert ( dbl[i] == qmckl_vec(vec, i) ) ;

  printf("qmckl_double_of_matrix ok\n");

  qmckl_exit_code rc = qmckl_free(context, dbl);
  assert (rc == QMCKL_SUCCESS);
  printf("qmckl_free ok\n");

  qmckl_vector_free(context, &vec);
  printf("qmckl_vector_free ok\n");

}
/* Tests:1 ends here */

/* [[file:../../org/qmckl_blas.org::*Test][Test:2]] */
qmckl_exit_code test_qmckl_dgemm(qmckl_context context);
assert(QMCKL_SUCCESS == test_qmckl_dgemm(context));
printf("qmckl_dgemm ok\n");
/* Test:2 ends here */

/* [[file:../../org/qmckl_blas.org::*Test][Test:2]] */
{
  double a[12] = { 1.,  5.,  9.,
                   2.,  6., 10.,
                   3.,  7., 11.,
                   4.,  8., 12. };

  double b[20] = { 1.,  5.,  9., 10.,
                   -2., -6., 10., 11.,
                   3.,  7., 11., 12.,
                   4.,  8., 12., 15.,
                   5.,  9., 13., 14. };

  double c[15] = { 39.,  89., 139.,
                   30.,  56.,  82.,
                   49., 115., 181.,
                   58., 136., 214.,
                   59., 141., 223. };

  qmckl_exit_code rc;
  qmckl_matrix A = qmckl_matrix_alloc(context, 3, 4);
  rc = qmckl_matrix_of_double(context, a, 12, &A);
  assert(rc == QMCKL_SUCCESS);
  printf("A ok\n");

  qmckl_matrix B = qmckl_matrix_alloc(context, 4, 5);
  rc = qmckl_matrix_of_double(context, b, 20, &B);
  assert(rc == QMCKL_SUCCESS);
  printf("B ok\n");

  qmckl_matrix C = qmckl_matrix_alloc(context, 3, 5);
  rc = qmckl_matmul(context, 'N', 'N', 0.5, A, B, 0., &C);
  printf("C ok\n");
  assert(rc == QMCKL_SUCCESS);

  double cnew[15];
  rc = qmckl_double_of_matrix(context, C, &(cnew[0]), 15);
  assert(rc == QMCKL_SUCCESS);
  printf("cnew ok\n");

  for (int i=0 ; i<15 ; ++i) {
    printf("%f %f\n", cnew[i], c[i]);
    assert (c[i] == cnew[i]);
  }
  printf("qmckl_matmul ok\n");
}
/* Test:2 ends here */

/* [[file:../../org/qmckl_blas.org::*Test][Test:5]] */
qmckl_exit_code test_qmckl_adjugate(qmckl_context context);
assert(QMCKL_SUCCESS == test_qmckl_adjugate(context));
printf("qmckl_adjugate ok\n");
/* Test:5 ends here */

/* [[file:../../org/qmckl_blas.org::*Test][Test:1]] */
{
  qmckl_matrix A;
  qmckl_matrix At;
  A  = qmckl_matrix_alloc(context, 2, 3);
  At = qmckl_matrix_alloc(context, 3, 2);
  for (int j=0 ; j<3 ; ++j)
    for (int i=0 ; i<2 ; ++i)
      qmckl_mat(A, i, j) = (double) 10*i+j;

  qmckl_exit_code rc = qmckl_transpose(context, A, At);
  assert(rc == QMCKL_SUCCESS);

  assert(A.size[0] == At.size[1]);
  assert(A.size[1] == At.size[0]);
  for (int j=0 ; j<3 ; ++j)
    for (int i=0 ; i<2 ; ++i)
      assert (qmckl_mat(A, i, j) == qmckl_mat(At, j, i));

printf("qmckl_transpose ok\n");
  qmckl_matrix_free(context, &A);
  qmckl_matrix_free(context, &At);
}
/* Test:1 ends here */

/* [[file:../../org/qmckl_blas.org::*End of files][End of files:3]] */
assert (qmckl_context_destroy(context) == QMCKL_SUCCESS);
  return 0;
}
/* End of files:3 ends here */
