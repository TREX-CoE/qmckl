#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_STDINT_H
#include <stdint.h>
#elif HAVE_INTTYPES_H
#include <inttypes.h>
#endif

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef HAVE_FPE
#define _GNU_SOURCE
#include <fenv.h>
#include <signal.h>
#include <stdio.h>
#include <execinfo.h>


#define MAX_BACKTRACE_SIZE 100

void floatingPointExceptionHandler(int signal) {
    void* backtraceArray[MAX_BACKTRACE_SIZE];
    int backtraceSize = backtrace(backtraceArray, MAX_BACKTRACE_SIZE);
    char** backtraceSymbols = backtrace_symbols(backtraceArray, backtraceSize);

    // Print the backtrace
    for (int i = 0; i < backtraceSize; ++i) {
            printf("[%d] %s\n", i, backtraceSymbols[i]);
    }

    // Clean up the memory used by backtrace_symbols
    free(backtraceSymbols);

    exit(EXIT_FAILURE);
}

static void __attribute__ ((constructor))
trapfpe ()
{
  /* Enable some exceptions.  At startup all exceptions are masked.  */

  feenableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
  signal(SIGFPE, floatingPointExceptionHandler);
}
#endif

#include "qmckl.h"
#include "qmckl_context_private_type.h"

qmckl_exit_code qmckl_set_numprec_precision(const qmckl_context context, const int precision) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT)
    return QMCKL_INVALID_CONTEXT;

  if (precision <  2) {
    return qmckl_failwith(context,
                          QMCKL_INVALID_ARG_2,
                          "qmckl_update_numprec_precision",
                          "precision < 2");
  }

  if (precision > 53) {
    return qmckl_failwith(context,
                          QMCKL_INVALID_ARG_2,
                          "qmckl_update_numprec_precision",
                          "precision > 53");
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;

  /* This should be always true because the context is valid */
  assert (ctx != NULL);

  qmckl_lock(context);
  {
    ctx->numprec.precision = (uint32_t) precision;
  }
  qmckl_unlock(context);

  return QMCKL_SUCCESS;
}

int qmckl_get_numprec_precision(const qmckl_context context) {
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
      return qmckl_failwith(context,
                      QMCKL_INVALID_CONTEXT,
                      "qmckl_get_numprec_precision",
                      "");
  }

  const qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  return ctx->numprec.precision;
}

qmckl_exit_code qmckl_set_numprec_range(const qmckl_context context, const int range) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT)
    return QMCKL_INVALID_CONTEXT;

  if (range <  2) {
    return qmckl_failwith(context,
                    QMCKL_INVALID_ARG_2,
                    "qmckl_set_numprec_range",
                    "range < 2");
  }

  if (range > 11) {
    return qmckl_failwith(context,
                    QMCKL_INVALID_ARG_2,
                    "qmckl_set_numprec_range",
                    "range > 11");
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;

  /* This should be always true because the context is valid */
  assert (ctx != NULL);

  qmckl_lock(context);
  {
    ctx->numprec.range = (uint32_t) range;
  }
  qmckl_unlock(context);

  return QMCKL_SUCCESS;
}

int qmckl_get_numprec_range(const qmckl_context context) {
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
      return qmckl_failwith(context,
                      QMCKL_INVALID_CONTEXT,
                      "qmckl_get_numprec_range",
                      "");
  }

  const qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  return ctx->numprec.range;
}

double qmckl_get_numprec_epsilon(const qmckl_context context) {
  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT)
    return QMCKL_INVALID_CONTEXT;
  const qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  const int precision = ctx->numprec.precision;
  return 1. /  (double) ( ((uint64_t) 1) << (precision-1));
}

int64_t countUlpDifference_64(double a, double b) {

  union int_or_float {
    int64_t i;
    double f;
  } x, y;

  x.f = a;
  y.f = b;

  // Handle sign bit discontinuity: if the signs are different and either value is not zero
  if ((x.i < 0) != (y.i < 0) && (x.f != 0.0) && (y.f != 0.0)) {
    // Use the absolute values and add the distance to zero for both numbers
    int64_t distanceToZeroForX = x.i < 0 ? INT64_MAX + x.i : INT64_MAX - x.i;
    int64_t distanceToZeroForY = y.i < 0 ? INT64_MAX + y.i : INT64_MAX - y.i;
    return distanceToZeroForX + distanceToZeroForY;
  }

  // Calculate the difference in their binary representations
  int64_t result = x.i - y.i;
  result = result > 0 ? result : -result;
  return result;
}

int32_t qmckl_test_precision_64(double a, double b) {

  int64_t diff = countUlpDifference_64(a,b);

  if (diff == 0) return 53;

  int32_t result = 53;

  for (int i=0 ; i<53 && diff != 0 ; ++i) {
    diff >>= 1;
    result--;
  }

  return result;
}

int32_t qmckl_test_precision_32(float a, float b) {
  return qmckl_test_precision_64( (double) a, (double) b );
}

float fastExpf(float x)
{
  const float a = 12102203.0;
  const float b = 1064986816.0;
  x = a * x + b;

  const float c = 8388608.0;
  const float d = 2139095040.0;
  if (x < c || x > d)
    x = (x < c) ? 0.0f : d;

  uint32_t n = (uint32_t) x;
  memcpy(&x, &n, 4);
  return x;
}


double fastExp(double x)
{
  const double a = 6497320848556798.0;
  const double b = 4606985713057410560.0;
  x = a * x + b;

  const double c = 4503599627370496.0;
  const double d = 9218868437227405312.0;
  if (x < c || x > d)
    x = (x < c) ? 0.0 : d;

  uint64_t n = (uint64_t) x;
  memcpy(&x, &n, 8);
  return x;
}
