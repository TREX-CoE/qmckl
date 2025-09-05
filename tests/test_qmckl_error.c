#include <string.h>
#include <stdio.h>
#include "qmckl.h"
#include "assert.h"
#include "qmckl_error_private_type.h"
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
int main() {
  qmckl_context context;
  context = qmckl_context_create();

/* [[file:../org/qmckl_error.org::*End of files][End of files:2]] */
/* Initialize the variables */
  char function_name[QMCKL_MAX_FUN_LEN]="";
  char message[QMCKL_MAX_MSG_LEN]="";

  /* Set the error code to be different from Success */
  qmckl_exit_code exit_code;
  exit_code = 1;

  assert (qmckl_set_error(context, exit_code, "my_function", "Message") == QMCKL_SUCCESS);

  assert (qmckl_get_error(context, &exit_code, function_name, message) == QMCKL_SUCCESS);
  assert (exit_code == 1);

  assert (strcmp(function_name,"my_function") == 0);
  assert (strcmp(message,"Message") == 0);

  exit_code = qmckl_context_destroy(context);
  assert(exit_code == QMCKL_SUCCESS);

  return 0;
}
/* End of files:2 ends here */
