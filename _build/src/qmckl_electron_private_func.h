#ifndef QMCKL_ELECTRON_HPF
#define QMCKL_ELECTRON_HPF
#include "qmckl_electron_private_type.h"
#include "qmckl_point_private_type.h"
#include "qmckl_point_private_func.h"



/* The ~uninitialized~ integer contains one bit set to one for each */
/* initialization function which has not been called. It becomes equal */
/* to zero after all initialization functions have been called. The */
/* struct is then initialized and ~provided == true~. */
/* Some values are initialized by default, and are not concerned by */
/* this mechanism. */


qmckl_exit_code qmckl_init_electron(qmckl_context context);

/* Deep copy */


qmckl_exit_code qmckl_copy_electron(qmckl_context context, const qmckl_electron_struct* src, qmckl_electron_struct* dest);

/* Provide                                                        :noexport: */


qmckl_exit_code qmckl_provide_ee_distance(qmckl_context context);

qmckl_exit_code qmckl_compute_ee_distance (
          const qmckl_context context,
          const int64_t elec_num,
          const int64_t walk_num,
          const double* coord,
          double* const ee_distance );

/* Provide                                                        :noexport: */


qmckl_exit_code qmckl_provide_ee_potential(qmckl_context context);

qmckl_exit_code qmckl_compute_ee_potential (const qmckl_context context,
                                            const int64_t elec_num,
                                            const int64_t walk_num,
                                            const double* ee_distance,
                                            double* const ee_potential );

/* Provide                                                        :noexport: */


qmckl_exit_code qmckl_provide_en_distance(qmckl_context context);

qmckl_exit_code qmckl_compute_en_distance (
          const qmckl_context context,
          const int64_t point_num,
          const int64_t nucl_num,
          const double* elec_coord,
          const double* nucl_coord,
          double* const en_distance );

/* Provide                                                        :noexport: */


qmckl_exit_code qmckl_provide_en_potential(qmckl_context context);

#endif
