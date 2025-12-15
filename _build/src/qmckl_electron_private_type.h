#ifndef QMCKL_ELECTRON_HPT
#define QMCKL_ELECTRON_HPT
#include <stdbool.h>
#include "qmckl_point_private_type.h"

/* Data structures */

/*    The electron data is organized using two levels of structures. The ~qmckl_walker~ */
/*    structure represents a collection of electronic configurations (walkers) being */
/*    sampled in the Monte Carlo process. Each walker contains multiple electrons at */
/*    specific positions. */
   
/*    The ~qmckl_electron_struct~ contains all electron-related data, including the */
/*    current and previous walker states (for computing acceptance ratios), and */
/*    cached computed quantities like inter-particle distances and potential energies. */
/*    The date fields enable automatic invalidation of cached results when electron */
/*    positions change. */


typedef struct qmckl_walker_struct {
  int64_t num;
  qmckl_point_struct point;
} qmckl_walker;

typedef struct qmckl_electron_struct {
  int64_t        num;
  int64_t        up_num;
  int64_t        down_num;
  qmckl_walker   walker;
  qmckl_walker   walker_old;
  uint64_t       ee_distance_date;
  uint64_t       en_distance_date;
  uint64_t       ee_potential_date;
  uint64_t       en_potential_date;
  double*        ee_distance;
  double*        en_distance;
  double*        ee_potential;
  double*        en_potential;
  int32_t        uninitialized;
  bool           provided;
} qmckl_electron_struct;

#endif
