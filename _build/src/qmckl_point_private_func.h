#ifndef QMCKL_POINT_HPF
#define QMCKL_POINT_HPF
#include "qmckl_point_private_type.h"
#include "qmckl_blas_private_type.h"
#include "qmckl_blas_private_func.h"

qmckl_exit_code qmckl_init_point(qmckl_context context);

/* Deep copy */


qmckl_exit_code qmckl_copy_point(qmckl_context context, const qmckl_point_struct* src, qmckl_point_struct* dest);

#endif
