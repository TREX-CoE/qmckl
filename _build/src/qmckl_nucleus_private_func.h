#ifndef QMCKL_NUCLEUS_HPF
#define QMCKL_NUCLEUS_HPF
#include "qmckl_nucleus_private_type.h"

/* Initialization tracking */

/*    The ~uninitialized~ integer serves as a bit field where each bit represents */
/*    whether a specific initialization function has been called. Initially, certain */
/*    bits are set to one, and each initialization function clears its corresponding */
/*    bit. When all required initialization functions have been called, the */
/*    ~uninitialized~ field becomes zero, the struct is fully initialized, and */
/*    ~provided~ is set to ~true~. */
   
/*    This mechanism ensures that the nucleus data is complete before being used */
/*    in calculations. Some values may be initialized by default and are not tracked */
/*    by this mechanism. */


qmckl_exit_code qmckl_init_nucleus(qmckl_context context);

/* Deep copy */


qmckl_exit_code qmckl_copy_nucleus(qmckl_context context, const qmckl_nucleus_struct* src, qmckl_nucleus_struct* dest);

/* Provide                                                        :noexport: */


qmckl_exit_code qmckl_provide_nn_distance(qmckl_context context);

qmckl_exit_code qmckl_compute_nn_distance (
          const qmckl_context context,
          const int64_t nucl_num,
          const double* coord,
          double* const nn_distance );

/* Provide                                                        :noexport: */


qmckl_exit_code qmckl_provide_nucleus_repulsion(qmckl_context context);

qmckl_exit_code qmckl_compute_nucleus_repulsion (
     const qmckl_context context,
     const int64_t nucl_num,
     const double* charge,
     const double* nn_distance,
     double* energy
  );

#endif
