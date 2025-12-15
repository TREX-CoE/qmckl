#ifndef QMCKL_NUCLEUS_HPT
#define QMCKL_NUCLEUS_HPT
#include <stdbool.h>

/* Data structure */

/*    The nucleus data structure consolidates all nuclear information in a single */
/*    C structure. The date fields enable automatic recomputation detection: when */
/*    coordinates change, the ~coord_date~ is updated, triggering recomputation */
/*    of dependent quantities like distances and repulsion energy when they are */
/*    next requested. */


typedef struct qmckl_nucleus_struct {
  int64_t      num;
  int64_t      repulsion_date;
  int64_t      nn_distance_date;
  int64_t      coord_date;
  qmckl_vector charge;
  qmckl_matrix coord;
  qmckl_matrix nn_distance;
  double       repulsion;
  int32_t      uninitialized;
  bool         provided;
} qmckl_nucleus_struct;

#endif
