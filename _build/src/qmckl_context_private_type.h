#ifndef QMCKL_CONTEXT_HPT
#define QMCKL_CONTEXT_HPT

#ifdef HAVE_STDINT_H
#include <stdint.h>
#elif HAVE_INTTYPES_H
#include <inttypes.h>
#endif

#include <pthread.h>

#include "qmckl_error_private_type.h"
#include "qmckl_memory_private_type.h"
#include "qmckl_numprec_private_type.h"
#include "qmckl_point_private_type.h"
#include "qmckl_nucleus_private_type.h"
#include "qmckl_electron_private_type.h"
#include "qmckl_ao_private_type.h"
#include "qmckl_mo_private_type.h"
#include "qmckl_jastrow_champ_private_type.h"
#include "qmckl_jastrow_champ_single_private_type.h"
#include "qmckl_forces_private_type.h"
#include "qmckl_determinant_private_type.h"
#include "qmckl_local_energy_private_type.h"
#include "qmckl_point_private_func.h"
#include "qmckl_nucleus_private_func.h"
#include "qmckl_electron_private_func.h"
#include "qmckl_ao_private_func.h"
#include "qmckl_mo_private_func.h"
#include "qmckl_jastrow_champ_private_func.h"
#include "qmckl_jastrow_champ_single_private_func.h"
#include "qmckl_forces_private_func.h"
#include "qmckl_determinant_private_func.h"
#include "qmckl_local_energy_private_func.h"

/* Data structure */

/*    The context data structure contains all the state information needed by the */
/*    library, organized into logical sections. This includes the molecular system */
/*    description (nuclei, electrons, basis sets), computational parameters */
/*    (numerical precision), computed quantities (atomic orbitals, molecular */
/*    orbitals, Jastrow factors, etc.), and internal state management (error */
/*    handling, memory allocation tracking, thread synchronization). */


typedef struct qmckl_context_struct {
  /* -- State of the library -- */

  /* Validity checking */
  uint64_t                 tag;

  /* Numerical precision */
  qmckl_numprec_struct     numprec;

  /* Thread lock */
  int                      lock_count;
  pthread_mutex_t          mutex;

  /* Error handling */
  qmckl_error_struct       error;

  /* Memory allocation */
  qmckl_memory_struct      memory;

  /* Current date */
  uint64_t                 date;

  /* Points */
  qmckl_point_struct         point;
  qmckl_jastrow_champ_single_struct single_point;

  /* -- Molecular system -- */
  qmckl_nucleus_struct        nucleus;
  qmckl_electron_struct       electron;
  qmckl_ao_basis_struct       ao_basis;
  qmckl_mo_basis_struct       mo_basis;
  qmckl_jastrow_champ_struct  jastrow_champ;
  qmckl_forces_struct         forces;
  qmckl_determinant_struct    det;
  qmckl_local_energy_struct   local_energy;

  /* To be implemented:
  */

  /* Pointer to implementation-specific data */

  void* qmckl_extra;

} qmckl_context_struct;

/* Validity checking */

/*    A tag is used internally to check if the memory domain pointed to */
/*    by a pointer is a valid context. This provides a sanity check that allows */
/*    verification that even if the pointer associated with a context is non-null,  */
/*    it actually points to the expected data structure and not to arbitrary memory. */
/*    This is particularly important for catching user errors like passing an */
/*    uninitialized context or a context that has already been destroyed. */


#define VALID_TAG   0xBEEFFACE
#define INVALID_TAG 0xDEADBEEF

/* [[file:../../org/qmckl_context.org::*End of files][End of files:1]] */
#endif
/* End of files:1 ends here */
