#ifndef QMCKL_LOCAL_ENERGY_HPT
#define QMCKL_LOCAL_ENERGY_HPT

#include <stdbool.h>

/* Data structure */


typedef struct qmckl_local_energy_struct {
  double  * e_kin;
  double  * e_pot;
  double  * e_local;
  double  * accep_prob;
  double  * r_drift;
  double  * y_move;
  uint64_t   e_kin_date;
  uint64_t   e_pot_date;
  uint64_t   e_local_date;
  uint64_t   accep_prob_date;
  uint64_t   r_drift_date;
  uint64_t   y_move_date;

  int32_t   uninitialized;
  bool      provided;
} qmckl_local_energy_struct;

#endif
