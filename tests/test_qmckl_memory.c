#include "qmckl.h"
#include "assert.h"
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "qmckl_context_private_type.h"
#include "qmckl_memory_private_func.h"
int main() {

/* Create a context */
qmckl_context context = qmckl_context_create();

qmckl_memory_info_struct info = qmckl_memory_info_struct_zero;
info.size = (size_t) 3*sizeof(int);

/* Allocate an array of ints */
int *a = (int*) qmckl_malloc(context, info);

/* Check that array of ints is OK */
assert(a != NULL);
a[0] = 1;  assert(a[0] == 1);
a[1] = 2;  assert(a[1] == 2);
a[2] = 3;  assert(a[2] == 3);

/* Allocate another array of ints */
int *b = (int*) qmckl_malloc(context, info);

/* Check that array of ints is OK */
assert(b != NULL);
b[0] = 1;  assert(b[0] == 1);
b[1] = 2;  assert(b[1] == 2);
b[2] = 3;  assert(b[2] == 3);

qmckl_exit_code rc;
/* Assert that both arrays are allocated */
assert(a != NULL);
assert(b != NULL);

/* Free in NULL context */
rc = qmckl_free(QMCKL_NULL_CONTEXT, a);
assert(rc == QMCKL_INVALID_CONTEXT);

/* Free NULL pointer */
rc = qmckl_free(context, NULL);
assert(rc == QMCKL_INVALID_ARG_2);

/* Free for the first time */
rc = qmckl_free(context, a);
assert(rc == QMCKL_SUCCESS);

/* Free again */
rc = qmckl_free(context, a);
assert(rc == QMCKL_INVALID_ARG_2);

/* Clean up */
rc = qmckl_context_destroy(context);
assert(rc == QMCKL_SUCCESS);

/* Create a context */
context = qmckl_context_create();

info = qmckl_memory_info_struct_zero;
info.size = (size_t) 3*sizeof(int);

/* Allocate an array of ints */
a = (int*) qmckl_malloc(context, info);

/* Check that the size of a is 3*sizeof(int) */
info = qmckl_memory_info_struct_zero;
rc = qmckl_get_malloc_info(context, NULL, &info);
assert (rc == QMCKL_INVALID_ARG_2);
rc = qmckl_get_malloc_info(context, &rc, &info);
assert (rc == QMCKL_INVALID_ARG_2);
rc = qmckl_get_malloc_info(context, a, &info);
assert (rc == QMCKL_SUCCESS);
assert (info.size == 3*sizeof(int));
rc = qmckl_context_destroy(context);

/* Test */

return 0;
}
