#include "qmckl.h"
#include "assert.h"
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
int main() {

/* [[file:../org/qmckl_context.org::*Creation][Creation:4]] */
assert( qmckl_context_check(QMCKL_NULL_CONTEXT) ==  QMCKL_NULL_CONTEXT);

qmckl_context context = qmckl_context_create();
assert( context != QMCKL_NULL_CONTEXT );
assert( qmckl_context_check(context) ==  context );
/* Creation:4 ends here */

/* [[file:../org/qmckl_context.org::*Copy][Copy:4]] */
qmckl_context new_context = qmckl_context_copy(context);
assert(new_context != QMCKL_NULL_CONTEXT);
assert(new_context != context);
assert(qmckl_context_check(new_context) == new_context);
qmckl_context_destroy(new_context);
/* Copy:4 ends here */

/* Destroy valid context */
assert(qmckl_context_check(context) == context);
assert(qmckl_context_destroy(context) == QMCKL_SUCCESS);

/* Check that context is destroyed  */
#ifndef DEBUG 
assert(qmckl_context_check(context) != context);
assert(qmckl_context_check(context) == QMCKL_NULL_CONTEXT);

/* Destroy invalid context */
assert(qmckl_context_destroy(QMCKL_NULL_CONTEXT) == QMCKL_INVALID_CONTEXT);
#endif

/* [[file:../org/qmckl_context.org::*Test][Test:1]] */
return 0;
}
/* Test:1 ends here */
