/* [[file:../../org/qmckl_distance.org::*Headers][Headers:2]] */
#include "qmckl.h"
#include "assert.h"
#include <stdio.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
int main() {
  qmckl_context context;
  context = qmckl_context_create();

#ifdef VFC_CI
  qmckl_init_probes();
#endif
/* Headers:2 ends here */

/* [[file:../../org/qmckl_distance.org::*Test][Test:2]] */
qmckl_exit_code test_qmckl_distance_sq(qmckl_context context);
assert(test_qmckl_distance_sq(context) == QMCKL_SUCCESS);
/* Test:2 ends here */

/* [[file:../../org/qmckl_distance.org::*Test][Test:2]] */
qmckl_exit_code test_qmckl_dist(qmckl_context context);
assert(test_qmckl_dist(context) == QMCKL_SUCCESS);
/* Test:2 ends here */

/* [[file:../../org/qmckl_distance.org::*Test][Test:3]] */
qmckl_exit_code test_qmckl_dist_rescaled(qmckl_context context);
assert(test_qmckl_dist_rescaled(context) == QMCKL_SUCCESS);
/* Test:3 ends here */

/* [[file:../../org/qmckl_distance.org::*End of files][End of files:1]] */
assert (qmckl_context_destroy(context) == QMCKL_SUCCESS);

#ifdef VFC_CI
  qmckl_dump_probes();
#endif
  return 0;
}
/* End of files:1 ends here */
