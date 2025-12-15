#ifndef QMCKL_FORCES_HPF
#define QMCKL_FORCES_HPF
#include "qmckl_forces_private_type.h"

/* Deep copy */


qmckl_exit_code qmckl_copy_forces(qmckl_context context, const qmckl_forces_struct* src, qmckl_forces_struct* dest);

/* Finite-difference function */

/* We introduce here a general function to compute the derivatives of any quantity with respect to nuclear coordinates. */
/* using finite-differences. */


typedef qmckl_exit_code (*function_callback)(qmckl_context context, double* const output, const int64_t size);

qmckl_exit_code qmckl_finite_difference_deriv_n(
    qmckl_context context,
    const double delta_x,                 // Step size for finite difference
    function_callback get_function,       // Function to compute values
    double* const derivative_output,      // Output derivative array: nucl_num*3*size
    int64_t const size);                  // Size of the object to differentiate

/* Provide                                                       :noexport: */

qmckl_exit_code qmckl_provide_forces_jastrow_en(qmckl_context context);

qmckl_exit_code qmckl_compute_forces_jastrow_en_doc (
          const qmckl_context context,
          const int64_t walk_num,
          const int64_t elec_num,
          const int64_t nucl_num,
          const int64_t type_nucl_num,
          const int64_t* type_nucl_vector,
          const int64_t aord_num,
          const double* a_vector,
          const double* en_distance_rescaled,
          const double* en_distance_rescaled_gl,
          double* const forces_jastrow_en );

qmckl_exit_code qmckl_compute_forces_jastrow_en (
          const qmckl_context context,
          const int64_t walk_num,
          const int64_t elec_num,
          const int64_t nucl_num,
          const int64_t type_nucl_num,
          const int64_t* type_nucl_vector,
          const int64_t aord_num,
          const double* a_vector,
          const double* en_distance_rescaled,
          const double* en_distance_rescaled_gl,
          double* const forces_jastrow_en );

qmckl_exit_code qmckl_compute_forces_jastrow_en (
          const qmckl_context context,
          const int64_t walk_num,
          const int64_t elec_num,
          const int64_t nucl_num,
          const int64_t type_nucl_num,
          const int64_t* type_nucl_vector,
          const int64_t aord_num,
          const double* a_vector,
          const double* en_distance_rescaled,
          const double* en_distance_rescaled_gl,
          double* const forces_jastrow_en );

/* Provide                                                       :noexport: */


qmckl_exit_code qmckl_provide_forces_jastrow_en_g(qmckl_context context);

qmckl_exit_code qmckl_compute_forces_jastrow_en_g_doc (
          const qmckl_context context,
          const int64_t walk_num,
          const int64_t elec_num,
          const int64_t nucl_num,
          const int64_t type_nucl_num,
          const int64_t* type_nucl_vector,
          const int64_t aord_num,
          const double* a_vector,
          const double* rescale_factor_en,
          const double* en_distance,
          const double* en_distance_rescaled,
          const double* en_distance_rescaled_gl,
          double* const forces_jastrow_en_g );

qmckl_exit_code qmckl_compute_forces_jastrow_en_g (
          const qmckl_context context,
          const int64_t walk_num,
          const int64_t elec_num,
          const int64_t nucl_num,
          const int64_t type_nucl_num,
          const int64_t* type_nucl_vector,
          const int64_t aord_num,
          const double* a_vector,
          const double* rescale_factor_en,
          const double* en_distance,
          const double* en_distance_rescaled,
          const double* en_distance_rescaled_gl,
          double* const forces_jastrow_en_g );

qmckl_exit_code qmckl_compute_forces_jastrow_en_g (
          const qmckl_context context,
          const int64_t walk_num,
          const int64_t elec_num,
          const int64_t nucl_num,
          const int64_t type_nucl_num,
          const int64_t* type_nucl_vector,
          const int64_t aord_num,
          const double* a_vector,
          const double* rescale_factor_en,
          const double* en_distance,
          const double* en_distance_rescaled,
          const double* en_distance_rescaled_gl,
          double* const forces_jastrow_en_g );

/* Provide                                                       :noexport: */

qmckl_exit_code qmckl_provide_forces_jastrow_en_l(qmckl_context context);

qmckl_exit_code qmckl_compute_forces_jastrow_en_l_doc (
          const qmckl_context context,
          const int64_t walk_num,
          const int64_t elec_num,
          const int64_t nucl_num,
          const int64_t type_nucl_num,
          const int64_t* type_nucl_vector,
          const int64_t aord_num,
          const double* a_vector,
          const double* rescale_factor_en,
          const double* en_distance,
          const double* en_distance_rescaled,
          const double* en_distance_rescaled_gl,
          double* const forces_jastrow_en_l );

qmckl_exit_code qmckl_compute_forces_jastrow_en_l (
          const qmckl_context context,
          const int64_t walk_num,
          const int64_t elec_num,
          const int64_t nucl_num,
          const int64_t type_nucl_num,
          const int64_t* type_nucl_vector,
          const int64_t aord_num,
          const double* a_vector,
          const double* rescale_factor_en,
          const double* en_distance,
          const double* en_distance_rescaled,
          const double* en_distance_rescaled_gl,
          double* const forces_jastrow_en_l );

qmckl_exit_code qmckl_compute_forces_jastrow_en_l (
          const qmckl_context context,
          const int64_t walk_num,
          const int64_t elec_num,
          const int64_t nucl_num,
          const int64_t type_nucl_num,
          const int64_t* type_nucl_vector,
          const int64_t aord_num,
          const double* a_vector,
          const double* rescale_factor_en,
          const double* en_distance,
          const double* en_distance_rescaled,
          const double* en_distance_rescaled_gl,
          double* const forces_jastrow_en_l );

/* Provide                                                         :noexport: */


qmckl_exit_code qmckl_provide_forces_jastrow_single_en(qmckl_context context);

qmckl_exit_code qmckl_compute_forces_jastrow_single_en_doc (
          const qmckl_context context,
          const int64_t num,
          const int64_t walk_num,
          const int64_t elec_num,
          const int64_t nucl_num,
          const int64_t type_nucl_num,
          const int64_t* type_nucl_vector,
          const int64_t aord_num,
          const double* a_vector,
          const double* en_distance_rescaled,
          const double* en_distance_rescaled_gl,
          const double* en_rescaled_single,
          const double* en_rescaled_single_gl,
          double* const forces_jastrow_single_en );

qmckl_exit_code qmckl_compute_forces_jastrow_single_en (
          const qmckl_context context,
          const int64_t num,
          const int64_t walk_num,
          const int64_t elec_num,
          const int64_t nucl_num,
          const int64_t type_nucl_num,
          const int64_t* type_nucl_vector,
          const int64_t aord_num,
          const double* a_vector,
          const double* en_distance_rescaled,
          const double* en_distance_rescaled_gl,
          const double* en_rescaled_single,
          const double* en_rescaled_single_gl,
          double* const forces_jastrow_single_en );

/* Provide                                                       :noexport: */

qmckl_exit_code qmckl_provide_forces_tmp_c(qmckl_context context);

qmckl_exit_code qmckl_compute_forces_tmp_c (
  const qmckl_context context,
  const int64_t walk_num,
  const int64_t elec_num,
  const int64_t nucl_num,
  const int64_t cord_num,
  const double* een_rescaled_e,
  const double* een_rescaled_n_gl,
  double* const forces_tmp_c );

/* Provide                                                       :noexport: */

qmckl_exit_code qmckl_provide_forces_dtmp_c(qmckl_context context);

qmckl_exit_code qmckl_compute_forces_dtmp_c_doc (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t cord_num,
      const double* een_rescaled_e_gl,
      const double* een_rescaled_n_gl,
      double* const forces_dtmp_c );

qmckl_exit_code qmckl_compute_forces_dtmp_c_hpc (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t cord_num,
      const double* een_rescaled_e_gl,
      const double* een_rescaled_n_gl,
      double* const forces_dtmp_c );

qmckl_exit_code qmckl_compute_forces_dtmp_c (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t cord_num,
      const double* een_rescaled_e_gl,
      const double* een_rescaled_n_gl,
      double* const forces_dtmp_c );

/* Provide                                                       :noexport: */


qmckl_exit_code qmckl_provide_forces_jastrow_een(qmckl_context context);

qmckl_exit_code qmckl_compute_forces_jastrow_een (
  const qmckl_context context,
  const int64_t walk_num,
  const int64_t elec_num,
  const int64_t nucl_num,
  const int64_t cord_num,
  const int64_t dim_c_vector,
  const double* c_vector_full,
  const int64_t* lkpm_combined_index,
  const double* een_rescaled_n,
  const double* een_rescaled_n_gl,
  const double* tmp_c,
  const double* forces_tmp_c,
  double* const forces_jastrow_een );

/* Provide                                                       :noexport: */


qmckl_exit_code qmckl_provide_forces_een_rescaled_n_gl(qmckl_context context);

qmckl_exit_code qmckl_compute_forces_een_rescaled_n_gl (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t type_nucl_num,
      int64_t* const type_nucl_vector,
      const int64_t cord_num,
      const double* rescale_factor_en,
      const double* en_distance,
      const double* een_rescaled_n,
      const double* een_rescaled_n_gl,
      double* const forces_een_n );

/* Provide                                                       :noexport: */

qmckl_exit_code qmckl_provide_forces_jastrow_een_g(qmckl_context context);

qmckl_exit_code qmckl_compute_forces_jastrow_een_g(
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t cord_num,
      const int64_t dim_c_vector,
      const double* c_vector_full,
      const int64_t* lkpm_combined_index,
      const double* en_distance,
      const double* tmp_c,
      const double* dtmp_c,
      const double* forces_tmp_c,
      const double* forces_dtmp_c,
      const double* forces_een_n,
      const double* een_rescaled_n,
      const double* een_rescaled_n_gl,
      double* const forces_jastrow_een_g );

/* Provide                                                       :noexport: */

qmckl_exit_code qmckl_provide_forces_jastrow_een_l(qmckl_context context);

qmckl_exit_code qmckl_compute_forces_jastrow_een_l(
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t cord_num,
      const int64_t dim_c_vector,
      const double* c_vector_full,
      const int64_t* lkpm_combined_index,
      const double* en_distance,
      const double* tmp_c,
      const double* dtmp_c,
      const double* forces_tmp_c,
      const double* forces_dtmp_c,
      const double* forces_een_n,
      const double* een_rescaled_n,
      const double* een_rescaled_n_gl,
      double* const forces_jastrow_een_l );

/* Provide                                                         :noexport: */


qmckl_exit_code qmckl_provide_forces_jastrow_single_een(qmckl_context context);

qmckl_exit_code qmckl_compute_forces_jastrow_single_een_doc (
          const qmckl_context context,
          const int64_t num,
          const int64_t walk_num,
          const int64_t elec_num,
          const int64_t nucl_num,
          const int64_t cord_num,
          const int64_t dim_c_vector,
          const double* c_vector_full,
          const int64_t* lkpm_combined_index,
          const double* een_rescaled_n,
          const double* een_rescaled_single_n,
          const double* een_rescaled_n_gl,
          const double* een_rescaled_single_n_gl,
          const double* een_rescaled_e,
          const double* een_rescaled_single_e,
          const double* tmp,
          double* const forces_jastrow_single_een );

qmckl_exit_code qmckl_compute_forces_jastrow_single_een (
          const qmckl_context context,
          const int64_t num,
          const int64_t walk_num,
          const int64_t elec_num,
          const int64_t nucl_num,
          const int64_t cord_num,
          const int64_t dim_c_vector,
          const double* c_vector_full,
          const int64_t* lkpm_combined_index,
          const double* een_rescaled_n,
          const double* een_rescaled_single_n,
          const double* een_rescaled_n_gl,
          const double* een_rescaled_single_n_gl,
          const double* een_rescaled_e,
          const double* een_rescaled_single_e,
          const double* tmp,
          double* const forces_jastrow_single_een );

/* Provide                                                       :noexport: */


qmckl_exit_code qmckl_provide_forces_ao_value(qmckl_context context);

qmckl_exit_code qmckl_compute_forces_ao_value_doc(
      const qmckl_context context,
      const int64_t ao_num,
      const int64_t shell_num,
      const int64_t point_num,
      const int64_t nucl_num,
      const int64_t* nucleus_index,
      const int64_t* nucleus_shell_num,
      const int32_t* shell_ang_mom,
      const double* ao_factor,
      const double* ao_vgl,
      double* const forces_ao_value );

/* Provide :noexport: */


qmckl_exit_code qmckl_provide_forces_mo_value(qmckl_context context);

qmckl_exit_code qmckl_compute_forces_mo_value_doc (
      const qmckl_context context,
      const int64_t nucl_num,
      const int64_t ao_num,
      const int64_t mo_num,
      const int64_t point_num,
      const int64_t shell_num,
      const int64_t* nucleus_index,
      const int64_t* nucleus_shell_num,
      const int32_t* shell_ang_mom,
      const double* coefficient_t,
      const double* ao_vgl,
      double* const forces_mo_value );

/* Provide :noexport: */


qmckl_exit_code qmckl_provide_forces_mo_g(qmckl_context context);

/* Provide :noexport: */


qmckl_exit_code qmckl_provide_forces_mo_l(qmckl_context context);

#endif
