#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_STDINT_H
#include <stdint.h>
#elif HAVE_INTTYPES_H
#include <inttypes.h>
#endif

#include <string.h>
#include <assert.h>
#include <pthread.h>
#include <errno.h>

#include "qmckl.h"
#include "qmckl_context_private_type.h"



/* #+RESULTS: cases */
/* #+begin_example */
/* case QMCKL_SUCCESS: */
/*       return "Success"; */

/* case QMCKL_INVALID_ARG_1: */
/*       return "Invalid argument 1"; */

/* case QMCKL_INVALID_ARG_2: */
/*       return "Invalid argument 2"; */

/* case QMCKL_INVALID_ARG_3: */
/*       return "Invalid argument 3"; */

/* case QMCKL_INVALID_ARG_4: */
/*       return "Invalid argument 4"; */

/* case QMCKL_INVALID_ARG_5: */
/*       return "Invalid argument 5"; */

/* case QMCKL_INVALID_ARG_6: */
/*       return "Invalid argument 6"; */

/* case QMCKL_INVALID_ARG_7: */
/*       return "Invalid argument 7"; */

/* case QMCKL_INVALID_ARG_8: */
/*       return "Invalid argument 8"; */

/* case QMCKL_INVALID_ARG_9: */
/*       return "Invalid argument 9"; */

/* case QMCKL_INVALID_ARG_10: */
/*       return "Invalid argument 10"; */

/* case QMCKL_INVALID_ARG_11: */
/*       return "Invalid argument 11"; */

/* case QMCKL_INVALID_ARG_12: */
/*       return "Invalid argument 12"; */

/* case QMCKL_INVALID_ARG_13: */
/*       return "Invalid argument 13"; */

/* case QMCKL_INVALID_ARG_14: */
/*       return "Invalid argument 14"; */

/* case QMCKL_INVALID_ARG_15: */
/*       return "Invalid argument 15"; */

/* case QMCKL_INVALID_ARG_16: */
/*       return "Invalid argument 16"; */

/* case QMCKL_INVALID_ARG_17: */
/*       return "Invalid argument 17"; */

/* case QMCKL_INVALID_ARG_18: */
/*       return "Invalid argument 18"; */

/* case QMCKL_INVALID_ARG_19: */
/*       return "Invalid argument 19"; */

/* case QMCKL_INVALID_ARG_20: */
/*       return "Invalid argument 20"; */

/* case QMCKL_FAILURE: */
/*       return "Failure"; */

/* case QMCKL_ERRNO: */
/*       return strerror(errno); */

/* case QMCKL_INVALID_CONTEXT: */
/*       return "Invalid context"; */

/* case QMCKL_ALLOCATION_FAILED: */
/*       return "Allocation failed"; */

/* case QMCKL_DEALLOCATION_FAILED: */
/*       return "De-allocation failed"; */

/* case QMCKL_NOT_PROVIDED: */
/*       return "Not provided"; */

/* case QMCKL_OUT_OF_BOUNDS: */
/*       return "Index out of bounds"; */

/* case QMCKL_ALREADY_SET: */
/*       return "Already set"; */

/* case QMCKL_INVALID_EXIT_CODE: */
/*       return "Invalid exit code"; */
/* #+end_example */

/* # Source */

const char* qmckl_string_of_error(const qmckl_exit_code error) {
  switch (error) {
  case QMCKL_SUCCESS:
        return "Success";
      
  case QMCKL_INVALID_ARG_1:
        return "Invalid argument 1";
      
  case QMCKL_INVALID_ARG_2:
        return "Invalid argument 2";
      
  case QMCKL_INVALID_ARG_3:
        return "Invalid argument 3";
      
  case QMCKL_INVALID_ARG_4:
        return "Invalid argument 4";
      
  case QMCKL_INVALID_ARG_5:
        return "Invalid argument 5";
      
  case QMCKL_INVALID_ARG_6:
        return "Invalid argument 6";
      
  case QMCKL_INVALID_ARG_7:
        return "Invalid argument 7";
      
  case QMCKL_INVALID_ARG_8:
        return "Invalid argument 8";
      
  case QMCKL_INVALID_ARG_9:
        return "Invalid argument 9";
      
  case QMCKL_INVALID_ARG_10:
        return "Invalid argument 10";
      
  case QMCKL_INVALID_ARG_11:
        return "Invalid argument 11";
      
  case QMCKL_INVALID_ARG_12:
        return "Invalid argument 12";
      
  case QMCKL_INVALID_ARG_13:
        return "Invalid argument 13";
      
  case QMCKL_INVALID_ARG_14:
        return "Invalid argument 14";
      
  case QMCKL_INVALID_ARG_15:
        return "Invalid argument 15";
      
  case QMCKL_INVALID_ARG_16:
        return "Invalid argument 16";
      
  case QMCKL_INVALID_ARG_17:
        return "Invalid argument 17";
      
  case QMCKL_INVALID_ARG_18:
        return "Invalid argument 18";
      
  case QMCKL_INVALID_ARG_19:
        return "Invalid argument 19";
      
  case QMCKL_INVALID_ARG_20:
        return "Invalid argument 20";
      
  case QMCKL_FAILURE:
        return "Failure";
      
  case QMCKL_ERRNO:
        return strerror(errno);
      
  case QMCKL_INVALID_CONTEXT:
        return "Invalid context";
      
  case QMCKL_ALLOCATION_FAILED:
        return "Allocation failed";
      
  case QMCKL_DEALLOCATION_FAILED:
        return "De-allocation failed";
      
  case QMCKL_NOT_PROVIDED:
        return "Not provided";
      
  case QMCKL_OUT_OF_BOUNDS:
        return "Index out of bounds";
      
  case QMCKL_ALREADY_SET:
        return "Already set";
      
  case QMCKL_INVALID_EXIT_CODE:
        return "Invalid exit code";
  }
  return "Unknown error";
}

void qmckl_string_of_error_f(const qmckl_exit_code error, char result[128]) {
  strncpy(result, qmckl_string_of_error(error), 128-1);
}

qmckl_exit_code
qmckl_set_error(qmckl_context context,
                const qmckl_exit_code exit_code,
                const char* function_name,
                const char* message)
{
  /* Passing a function name and a message is mandatory. */
  assert (function_name != NULL);
  assert (message != NULL);

  /* Exit codes are assumed valid. */
  assert (exit_code >= 0);
  assert (exit_code != QMCKL_SUCCESS);
  assert (exit_code < QMCKL_INVALID_EXIT_CODE);

  /* The context is assumed to exist. */
  assert (qmckl_context_check(context) != QMCKL_NULL_CONTEXT);

  qmckl_lock(context);
  {
    qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
    assert (ctx != NULL); /* Impossible because the context is valid. */

    ctx->error.exit_code = exit_code;
    strncpy(ctx->error.function, function_name, QMCKL_MAX_FUN_LEN-1);
    strncpy(ctx->error.message, message, QMCKL_MAX_MSG_LEN-1);
  }
  qmckl_unlock(context);

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_error(qmckl_context context,
                qmckl_exit_code *exit_code,
                char* function_name,
                char* message)
{
  /* Passing a function name and a message is mandatory. */
  assert (function_name != NULL);
  assert (message != NULL);

  /* The context is assumed to exist. */
  assert (qmckl_context_check(context) != QMCKL_NULL_CONTEXT);

  qmckl_lock(context);
  {
    qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
    assert (ctx != NULL); /* Impossible because the context is valid. */

    size_t sizeCp;

    sizeCp = strlen(ctx->error.function);
    sizeCp = sizeCp > QMCKL_MAX_FUN_LEN ? QMCKL_MAX_FUN_LEN : sizeCp;
    memcpy(function_name, ctx->error.function, sizeCp);

    sizeCp = strlen(ctx->error.message);
    sizeCp = sizeCp > QMCKL_MAX_MSG_LEN ? QMCKL_MAX_MSG_LEN : sizeCp;
    memcpy(message, ctx->error.message, sizeCp);

    (*exit_code) = ctx->error.exit_code;
  }
  qmckl_unlock(context);

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_failwith(qmckl_context context,
               const qmckl_exit_code exit_code,
               const char* function,
               const char* message)
{
  assert (exit_code > 0);
  assert (exit_code < QMCKL_INVALID_EXIT_CODE);
  assert (function != NULL);
  assert (strlen(function) < QMCKL_MAX_FUN_LEN);
  if (message != NULL) {
    assert (strlen(message)  < QMCKL_MAX_MSG_LEN);
  }

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT)
    return QMCKL_INVALID_CONTEXT;

  if (message == NULL) {
    qmckl_exit_code rc =
      qmckl_set_error(context, exit_code, function, qmckl_string_of_error(exit_code));
    assert (rc == QMCKL_SUCCESS);
  } else {
    qmckl_exit_code rc =
      qmckl_set_error(context, exit_code, function, message);
    assert (rc == QMCKL_SUCCESS);
  }

  return exit_code;
}

qmckl_exit_code
qmckl_last_error(qmckl_context context, char* buffer) {

  char function_name[QMCKL_MAX_FUN_LEN];
  char message[QMCKL_MAX_MSG_LEN];

  qmckl_exit_code rc, last_rc;

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    strncpy(buffer, "Null context", 13);
    return QMCKL_FAILURE;
  }

  rc = qmckl_get_error(context, &last_rc, function_name, message);
  if (rc != QMCKL_SUCCESS) {
    return rc;
  }


  sprintf(buffer, "Error -- %s -- in %s\n%s",
          qmckl_string_of_error(last_rc),
          function_name, message);

  return QMCKL_SUCCESS;
}

#include <stdio.h>

qmckl_exit_code
qmckl_check(qmckl_context context, qmckl_exit_code rc)
{

  if (rc != QMCKL_SUCCESS) {
    char fname[QMCKL_MAX_FUN_LEN];
    char message[QMCKL_MAX_MSG_LEN];

    fprintf(stderr, "===========\nQMCKL ERROR\n%s\n", qmckl_string_of_error(rc));
    qmckl_get_error(context, &rc, fname, message);
    fprintf(stderr, "Function: %s\nMessage: %s\n===========\n", fname, message);
  }

  return rc;
}
