#ifndef QMCKL_MEMORY_HPT
#define QMCKL_MEMORY_HPT

#include <stdint.h>
#include <stdlib.h>

typedef struct qmckl_memory_info_struct {
  size_t size;
  void*  pointer;
} qmckl_memory_info_struct;

static const qmckl_memory_info_struct qmckl_memory_info_struct_zero =
  {
   .size = (size_t) 0,
   .pointer = NULL
  };

typedef struct qmckl_memory_struct {
  size_t                    n_allocated;
  size_t                    array_size;
  qmckl_memory_info_struct* element;
} qmckl_memory_struct;

#endif
