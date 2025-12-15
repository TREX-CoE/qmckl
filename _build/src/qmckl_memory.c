#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_STDINT_H
#include <stdint.h>
#elif HAVE_INTTYPES_H
#include <inttypes.h>
#endif

#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "qmckl.h"
#include "qmckl_memory_private_type.h"
#include "qmckl_context_private_type.h"
#include "qmckl_memory_private_func.h"

void* qmckl_malloc(qmckl_context context, const qmckl_memory_info_struct info) {

  assert (qmckl_context_check(context) != QMCKL_NULL_CONTEXT);

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;

  /* Allocate memory and zero it */
  void * pointer = NULL;
#if defined(HAVE_HPC) && defined(HAVE_POSIX_MEMALIGN)
  if (posix_memalign(&pointer, 64, info.size) != 0) pointer = NULL;
#else
  pointer = malloc(info.size);
#endif
  if (pointer == NULL) {
    return NULL;
  }
  memset(pointer, 0, info.size);

  qmckl_lock(context);
  {
    /* If qmckl_memory_struct is full, reallocate a larger one */
    if (ctx->memory.n_allocated == ctx->memory.array_size) {
      const size_t old_size = ctx->memory.array_size;
      qmckl_memory_info_struct * new_array = realloc(ctx->memory.element,
                                                          2L * old_size *
                                                          sizeof(qmckl_memory_info_struct));
      if (new_array == NULL) {
        qmckl_unlock(context);
        free(pointer);
        return NULL;
      }

      memset( &(new_array[old_size]), 0, old_size * sizeof(qmckl_memory_info_struct) );
      ctx->memory.element = new_array;
      ctx->memory.array_size = 2L * old_size;
    }

    /* Find first NULL entry */
    size_t pos = (size_t) 0;
    while ( pos < ctx->memory.array_size && ctx->memory.element[pos].size > (size_t) 0) {
      pos += (size_t) 1;
    }
    assert (ctx->memory.element[pos].size == (size_t) 0);

    /* Copy info at the new location */
    memcpy(&(ctx->memory.element[pos]), &info, sizeof(qmckl_memory_info_struct));
    ctx->memory.element[pos].pointer = pointer;
    ctx->memory.n_allocated += (size_t) 1;
//printf("MALLOC: %5ld  %p\n", ctx->memory.n_allocated, ctx->memory.element[pos].pointer);
  }
  qmckl_unlock(context);

  return pointer;
}

qmckl_exit_code qmckl_free(qmckl_context context, void * const ptr) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
      return qmckl_failwith(context,
                            QMCKL_INVALID_CONTEXT,
                            "qmckl_free",
                            NULL);
  }

  if (ptr == NULL) {
    return qmckl_failwith(context,
                          QMCKL_INVALID_ARG_2,
                          "qmckl_free",
                          "NULL pointer");
  }

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;

  qmckl_lock(context);
  {
    /* Find pointer in array of saved pointers */
    size_t pos = (size_t) 0;
    while ( pos < ctx->memory.array_size && ctx->memory.element[pos].pointer != ptr) {
      pos += (size_t) 1;
    }

    if (pos >= ctx->memory.array_size) {
      /* Not found */
      qmckl_unlock(context);
      return qmckl_failwith(context,
                            QMCKL_INVALID_ARG_2,
                            "qmckl_free",
                            "Pointer not found in context");
    }

    /* Found */

    free(ptr);

    ctx->memory.n_allocated -= (size_t) 1;
//printf("FREE  : %5ld  %p\n", ctx->memory.n_allocated, ctx->memory.element[pos].pointer);
    ctx->memory.element[pos] = qmckl_memory_info_struct_zero;
  }
  qmckl_unlock(context);

  return QMCKL_SUCCESS;
}

qmckl_exit_code
qmckl_get_malloc_info(qmckl_context context,
                      const void* ptr, 
                      qmckl_memory_info_struct* info) 
{

  assert (qmckl_context_check(context) != QMCKL_NULL_CONTEXT);

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;

  if (ptr == NULL) {
    return qmckl_failwith(context,
                          QMCKL_INVALID_ARG_2,
                          "qmckl_get_malloc_info",
                          "Null pointer");
  }

  if (info == NULL) {
    return qmckl_failwith(context,
                          QMCKL_INVALID_ARG_3,
                          "qmckl_get_malloc_info",
                          "Null pointer");
  }

  qmckl_lock(context);
  {
    /* Find the pointer entry */
    size_t pos = (size_t) 0;
    while ( pos < ctx->memory.array_size && ctx->memory.element[pos].pointer != ptr) {
      pos += (size_t) 1;
    }

    if (pos >= ctx->memory.array_size) {
      /* Not found */
      qmckl_unlock(context);
      return qmckl_failwith(context,
                            QMCKL_INVALID_ARG_2,
                            "qmckl_get_malloc_info",
                            "Pointer not found in context");
    }

    /* Copy info */
    memcpy(info, &(ctx->memory.element[pos]), sizeof(qmckl_memory_info_struct));
  }
  qmckl_unlock(context);

  return QMCKL_SUCCESS;
}
