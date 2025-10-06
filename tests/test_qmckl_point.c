#include "qmckl.h"
#include <assert.h>
#include <math.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "chbrclf.h"
#include "qmckl_blas_private_type.h"
#include "qmckl_blas_private_func.h"

int main() {
  qmckl_context context;
  context = qmckl_context_create();

/* Reference input data */
int64_t point_num = chbrclf_elec_num;
const double* coord     = &(chbrclf_elec_coord[0][0][0]);

/* --- */

qmckl_exit_code rc;
double coord2[point_num*3];
double coord3[point_num*3];


rc = qmckl_get_point (context, 'N', coord2, (point_num*3));
assert(rc == QMCKL_NOT_PROVIDED);

rc = qmckl_set_point (context, 'N', point_num, coord, (point_num*3));
assert(rc == QMCKL_SUCCESS);

int64_t n;
rc = qmckl_get_point_num (context, &n);
assert(rc == QMCKL_SUCCESS);
assert(n == point_num);

rc = qmckl_get_point (context, 'N', coord2, (point_num*3));
assert(rc == QMCKL_SUCCESS);

for (int64_t i=0 ; i<3*point_num ; ++i) {
  assert( coord[i] == coord2[i] );
}

rc = qmckl_get_point (context, 'T', coord2, (point_num*3));
assert(rc == QMCKL_SUCCESS);

for (int64_t i=0 ; i<point_num ; ++i) {
  assert( coord[3*i+0] == coord2[i] );
  assert( coord[3*i+1] == coord2[i+point_num] );
  assert( coord[3*i+2] == coord2[i+point_num*2] );
}

rc = qmckl_set_point (context, 'T', point_num, coord2, (point_num*3));
assert(rc == QMCKL_SUCCESS);

rc = qmckl_get_point (context, 'N', coord3, (point_num*3));
assert(rc == QMCKL_SUCCESS);

for (int64_t i=0 ; i<3*point_num ; ++i) {
  assert( coord[i] == coord3[i] );
}

if (qmckl_context_destroy(context) != QMCKL_SUCCESS)
    return QMCKL_FAILURE;
  return 0;
}
