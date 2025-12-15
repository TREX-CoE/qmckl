#ifndef QMCKL_MEMORY_HPF
#define QMCKL_MEMORY_HPF

void* qmckl_malloc(qmckl_context context,
                   const qmckl_memory_info_struct info);

qmckl_exit_code qmckl_free(qmckl_context context,
                           void * const ptr);

qmckl_exit_code
qmckl_get_malloc_info(qmckl_context context,
                      const void* pointer, 
                      qmckl_memory_info_struct* info);

/* End of files                                                     :noexport: */


#endif
