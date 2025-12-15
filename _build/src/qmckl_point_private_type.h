#ifndef QMCKL_POINT_HPT
#define QMCKL_POINT_HPT
#include <stdbool.h>
#include "qmckl_blas_private_type.h"

/* Data structure */

/*   The following data stored in the context: */

/*   | Variable     | Type           | Description                               | */
/*   |--------------+----------------+-------------------------------------------| */
/*   | ~num~        | ~int64_t~      | Total number of points                    | */
/*   | ~date~       | ~uint64_t~     | Last modification date of the coordinates | */
/*   | ~coord~      | ~qmckl_matrix~ | ~num~ \times 3 matrix                     | */

/*   We consider that the matrix is stored 'transposed' and 'normal' */
/*   corresponds to the 3 \times ~num~ matrix. */



typedef struct qmckl_point_struct {
  int64_t      num;
  uint64_t     date;
  qmckl_matrix coord;
} qmckl_point_struct;

#endif
