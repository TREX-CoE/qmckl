#ifndef QMCKL_AO_HPF
#define QMCKL_AO_HPF
#include "qmckl_ao_private_type.h"
#include "qmckl_blas_private_type.h"



/* The ~uninitialized~ integer contains one bit set to one for each */
/* initialization function which has not been called. It becomes equal */
/* to zero after all initialization functions have been called. The */
/* struct is then initialized and ~provided == true~. */
/* Some values are initialized by default, and are not concerned by */
/* this mechanism. */


qmckl_exit_code qmckl_init_ao_basis(qmckl_context context);

/* Deep copy */


qmckl_exit_code qmckl_copy_ao_basis(qmckl_context context, const qmckl_ao_basis_struct* src, qmckl_ao_basis_struct* dest);

/* After initialization */

/*     When the basis set is completely entered, extra data structures may be */
/*     computed to accelerate the calculations. The primitives within each */
/*     contraction are sorted in ascending order of their exponents, such */
/*     that as soon as a primitive is zero all the following functions */
/*     vanish. Also, it is possible to compute a nuclear radius beyond which */
/*     all the primitives are zero up to the numerical accuracy defined in */
/*     the context. */


qmckl_exit_code qmckl_finalize_basis (qmckl_context context);

#ifdef HAVE_HPC
qmckl_exit_code qmckl_finalize_basis_hpc (qmckl_context context);
#endif

/* Provide                                                        :noexport: */

/* #+CALL: write_provider_header( group="ao_basis", data="primitive_vgl" ) */

/* #+RESULTS: */

qmckl_exit_code qmckl_provide_ao_basis_primitive_vgl(qmckl_context context);

/* Provide                                                        :noexport: */

/* #+CALL: write_provider_header( group="ao_basis", data="shell_vgl" ) */

/* #+RESULTS: */

qmckl_exit_code qmckl_provide_ao_basis_shell_vgl(qmckl_context context);

/* Provide                                                        :noexport: */


qmckl_exit_code qmckl_provide_ao_basis_shell_hessian(qmckl_context context);

/* Determination of nucleus ranges */

/*    The range of a nucleus is defined as the distance beyond which all */
/*    the AOs are zero, up to a given precision. For all possible numbers */
/*    of bits of precision (1-53), a range is given. */


qmckl_exit_code qmckl_compute_nucleus_range_gaussian(
          const qmckl_context context,
          const int64_t ao_num,
          const int64_t shell_num,
          const int64_t prim_num,
          const int64_t nucl_num,
          const double* nucl_coord,
          const int64_t* nucleus_index,
          const int64_t* nucleus_shell_num,
          const int32_t* nucleus_max_ang_mom,
          const int64_t* shell_prim_index,
          const int64_t* shell_prim_num,
          const int32_t* r_power,
          const double* ao_factor,
          const double* ao_expo,
          const double* coefficient_normalized,
          double* const nucleus_range);

qmckl_exit_code qmckl_compute_nucleus_range_slater(
          const qmckl_context context,
          const int64_t ao_num,
          const int64_t shell_num,
          const int64_t prim_num,
          const int64_t nucl_num,
          const double* nucl_coord,
          const int64_t* nucleus_index,
          const int64_t* nucleus_shell_num,
          const int32_t* nucleus_max_ang_mom,
          const int64_t* shell_prim_index,
          const int64_t* shell_prim_num,
          const int32_t* r_power,
          const double* ao_factor,
          const double* expo,
          const double* coefficient_normalized,
          double* const nucleus_range);

/* Interfaces */
/*  #  #+CALL: generate_private_c_header(table=qmckl_ao_value_args_doc,rettyp=get_value("CRetType"),fname="qmckl_compute_ao_value")) */
/*  #  (Commented because the header needs to go into h_private_func) */


qmckl_exit_code qmckl_compute_ao_value_doc (
      const qmckl_context context,
      const int64_t ao_num,
      const int64_t shell_num,
      const int64_t point_num,
      const int64_t nucl_num,
      const double* coord,
      const double* nucl_coord,
      const int64_t* nucleus_index,
      const int64_t* nucleus_shell_num,
      const double* nucleus_range,
      const int32_t* nucleus_max_ang_mom,
      const int32_t* shell_ang_mom,
      const double* ao_factor,
      const double* shell_vgl,
      double* const ao_value );

#ifdef HAVE_HPC
    qmckl_exit_code qmckl_compute_ao_value_hpc_gaussian (
          const qmckl_context context,
          const int64_t ao_num,
          const int64_t shell_num,
          const int32_t* prim_num_per_nucleus,
          const int64_t point_num,
          const int64_t nucl_num,
          const double* coord,
          const double* nucl_coord,
          const int64_t* nucleus_index,
          const int64_t* nucleus_shell_num,
          const double* nucleus_range,
          const int32_t* nucleus_max_ang_mom,
          const int32_t* shell_ang_mom,
          const double* ao_factor,
          const qmckl_matrix expo_per_nucleus,
          const qmckl_tensor coef_per_nucleus,
          double* const ao_value );
#endif

/* Provide                                                       :noexport: */

/* #+CALL: write_provider_header( group="ao_basis", data="ao_value" ) */

/* #+RESULTS: */

qmckl_exit_code qmckl_provide_ao_basis_ao_value(qmckl_context context);

/* Interfaces */
/*     #+CALL: generate_private_c_header(table=qmckl_ao_vgl_args_doc,rettyp=get_value("CRetType"),fname="qmckl_compute_ao_vgl") */


qmckl_exit_code qmckl_compute_ao_vgl_doc (
      const qmckl_context context,
      const int64_t ao_num,
      const int64_t shell_num,
      const int64_t point_num,
      const int64_t nucl_num,
      const double* coord,
      const double* nucl_coord,
      const int64_t* nucleus_index,
      const int64_t* nucleus_shell_num,
      const double* nucleus_range,
      const int32_t* nucleus_max_ang_mom,
      const int32_t* shell_ang_mom,
      const double* ao_factor,
      const double* shell_vgl,
      double* const ao_vgl );

#ifdef HAVE_HPC
    qmckl_exit_code qmckl_compute_ao_vgl_hpc_gaussian (
          const qmckl_context context,
          const int64_t ao_num,
          const int64_t shell_num,
          const int32_t* prim_num_per_nucleus,
          const int64_t point_num,
          const int64_t nucl_num,
          const double* coord,
          const double* nucl_coord,
          const int64_t* nucleus_index,
          const int64_t* nucleus_shell_num,
          const double* nucleus_range,
          const int32_t* nucleus_max_ang_mom,
          const int32_t* shell_ang_mom,
          const double* ao_factor,
          const qmckl_matrix expo_per_nucleus,
          const qmckl_tensor coef_per_nucleus,
          double* const ao_vgl );
#endif

/* Provide                                                       :noexport: */

/* #+CALL: write_provider_header( group="ao_basis", data="ao_vgl" ) */

/* #+RESULTS: */

qmckl_exit_code qmckl_provide_ao_basis_ao_vgl(qmckl_context context);

/* Provide :noexport: */


qmckl_exit_code qmckl_provide_ao_basis_ao_hessian(qmckl_context context);

qmckl_exit_code qmckl_compute_ao_hessian_doc (
      const qmckl_context context,
      const int64_t ao_num,
      const int64_t shell_num,
      const int64_t point_num,
      const int64_t nucl_num,
      const double* coord,
      const double* nucl_coord,
      const int64_t* nucleus_index,
      const int64_t* nucleus_shell_num,
      const double* nucleus_range,
      const int32_t* nucleus_max_ang_mom,
      const int32_t* shell_ang_mom,
      const double* ao_factor,
      const double* shell_vgl,
      const double* shell_hessian,
      double* const ao_hessian );

#endif
