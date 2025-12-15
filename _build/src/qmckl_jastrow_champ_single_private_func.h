#ifndef QMCKL_JASTROW_CHAMP_SINGLE_HPF
#define QMCKL_JASTROW_CHAMP_SINGLE_HPF

/* Deep copy */


#include "qmckl_jastrow_champ_single_private_type.h"
qmckl_exit_code qmckl_copy_jastrow_champ_single(qmckl_context context, const qmckl_jastrow_champ_single_struct* src, qmckl_jastrow_champ_single_struct* dest);



/* The Fortran function shifts the ~num~ by 1 because of 1-based */
/* indexing. */


qmckl_exit_code qmckl_set_single_point_f (qmckl_context context,
                                 const char transp,
                                 const int64_t num,
                                 const double* coord,
                                 const int64_t size_max);

/* Provide                                                        :noexport: */


qmckl_exit_code qmckl_provide_single_ee_distance(qmckl_context context);

qmckl_exit_code qmckl_compute_single_ee_distance (
          const qmckl_context context,
          const int64_t num,
          const int64_t elec_num,
          const int64_t walk_num,
          const double* coord,
          const double* single_coord,
          double* const single_ee_distance );

/* Provide                                                        :noexport: */


qmckl_exit_code qmckl_provide_single_en_distance(qmckl_context context);

qmckl_exit_code qmckl_compute_single_en_distance (
          const qmckl_context context,
          const int64_t nucl_num,
          const int64_t walk_num,
          const double* elec_coord,
          const double* nucl_coord,
          double* const single_en_distance );

/* Provide                                                       :noexport: */


qmckl_exit_code qmckl_provide_een_rescaled_single_e(qmckl_context context);

qmckl_exit_code qmckl_compute_een_rescaled_single_e (
      const qmckl_context context,
      const int64_t num,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t cord_num,
      const double rescale_factor_ee,
      const double* single_ee_distance,
      const double* een_rescaled_e,
      double* const een_rescaled_single_e );

qmckl_exit_code qmckl_compute_een_rescaled_single_e_doc (
      const qmckl_context context,
      const int64_t num,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t cord_num,
      const double rescale_factor_ee,
      const double* single_ee_distance,
      const double* een_rescaled_e,
      double* const een_rescaled_single_e );

/* Provide                                                       :noexport: */


qmckl_exit_code qmckl_provide_een_rescaled_single_n(qmckl_context context);

qmckl_exit_code qmckl_compute_een_rescaled_single_n (
      const qmckl_context context,
      const int64_t num,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t type_nucl_num,
      int64_t* const type_nucl_vector,
      const int64_t cord_num,
      const double* rescale_factor_en,
      const double* single_en_distance,
      const double* een_rescaled_n,
      double* const een_rescaled_single_n );

/* Provide                                                       :noexport: */


qmckl_exit_code qmckl_provide_jastrow_champ_single_een(qmckl_context context);

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_single_een_doc (const qmckl_context context,
                              const int64_t num,
                              const int64_t walk_num,
                              const int64_t elec_num,
                              const int64_t nucl_num,
                              const int64_t cord_num,
                              const int64_t dim_c_vector,
                              const double* c_vector_full,
                              const int64_t* lkpm_combined_index,
                              const double* een_rescaled_n,
                              const double* een_rescaled_e,
                              const double* een_rescaled_single_n,
                              const double* een_rescaled_single_e,
                              double* const delta_een );

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_single_een_hpc_c (const qmckl_context context,
                              const int64_t num,
                              const int64_t walk_num,
                              const int64_t elec_num,
                              const int64_t nucl_num,
                              const int64_t cord_num,
                              const int64_t dim_c_vector,
                              const double* c_vector_full,
                              const int64_t* lkpm_combined_index,
                              const double* een_rescaled_n,
                              const double* een_rescaled_e,
                              const double* een_rescaled_single_n,
                              const double* een_rescaled_single_e,
                              double* const delta_een );
qmckl_exit_code
qmckl_compute_jastrow_champ_factor_single_een_hpc_f (const qmckl_context context,
                              const int64_t num,
                              const int64_t walk_num,
                              const int64_t elec_num,
                              const int64_t nucl_num,
                              const int64_t cord_num,
                              const int64_t dim_c_vector,
                              const double* c_vector_full,
                              const int64_t* lkpm_combined_index,
                              const double* een_rescaled_n,
                              const double* een_rescaled_e,
                              const double* een_rescaled_single_n,
                              const double* een_rescaled_single_e,
                              double* const delta_een );

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_single_een (const qmckl_context context,
                          const int64_t num,
                          const int64_t walk_num,
                          const int64_t elec_num,
                          const int64_t nucl_num,
                          const int64_t cord_num,
                          const int64_t dim_c_vector,
                          const double* c_vector_full,
                          const int64_t* lkpm_combined_index,
                          const double* een_rescaled_n,
                          const double* een_rescaled_e,
                          const double* een_rescaled_single_n,
                          const double* een_rescaled_single_e,
                          double* const delta_een );

/* Provide                                                       :noexport: */


qmckl_exit_code qmckl_provide_een_rescaled_single_n_gl(qmckl_context context);

qmckl_exit_code qmckl_compute_een_rescaled_single_n_gl (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t nucl_num,
      const int64_t type_nucl_num,
      int64_t* const type_nucl_vector,
      const int64_t cord_num,
      const double* rescale_factor_en,
      const double* coord_ee,
      const double* coord_n,
      const double* single_en_distance,
      const double* een_rescaled_single_n,
      double* const een_rescaled_single_n_gl );

/* Provide                                                       :noexport: */


qmckl_exit_code qmckl_provide_een_rescaled_single_e_gl(qmckl_context context);

qmckl_exit_code qmckl_compute_een_rescaled_single_e_gl (
      const qmckl_context context,
      const int64_t num,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t cord_num,
      const double rescale_factor_ee,
      const double* coord,
      const double* coord_ee,
      const double* single_ee_distance,
      const double* een_rescaled_single_e,
      double* const een_rescaled_single_e_gl );

qmckl_exit_code qmckl_compute_een_rescaled_single_e_gl_doc (
      const qmckl_context context,
       const int64_t num,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t cord_num,
      const double rescale_factor_ee,
      const double* coord,
      const double* coord_ee,
      const double* single_ee_distance,
      const double* een_rescaled_single_e,
      double* const een_rescaled_single_e_gl );

/* Provide                                                       :noexport: */


qmckl_exit_code qmckl_provide_jastrow_champ_single_tmp(qmckl_context context);

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_single_tmp_doc (const qmckl_context context,
                              const int64_t num,
                              const int64_t walk_num,
                              const int64_t elec_num,
                              const int64_t nucl_num,
                              const int64_t cord_num,
                              const int64_t dim_c_vector,
                              const double* c_vector_full,
                              const int64_t* lkpm_combined_index,
                              const double* een_rescaled_n,
                              const double* een_rescaled_e,
                              const double* een_rescaled_single_e,
                              double* const tmp);

/* Provide                                                       :noexport: */


qmckl_exit_code qmckl_provide_jastrow_champ_single_een_gl(qmckl_context context);

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_single_een_gl_doc (const qmckl_context context,
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
                              const double* een_rescaled_e_gl,
                              const double* een_rescaled_single_e_gl,
                              double* const delta_een_gl );

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_single_een_gl_hpc (const qmckl_context context,
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
                              const double* een_rescaled_e_gl,
                              const double* een_rescaled_single_e_gl,
                              double* const delta_een_gl );


qmckl_exit_code
qmckl_compute_jastrow_champ_factor_single_een_gl (const qmckl_context context,
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
                              const double* een_rescaled_e_gl,
                              const double* een_rescaled_single_e_gl,
                          double* const delta_een_gl );

/* Provide                                                       :noexport: */


qmckl_exit_code qmckl_provide_jastrow_champ_single_een_g(qmckl_context context);

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_single_een_g_doc (const qmckl_context context,
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
                              const double* een_rescaled_e_gl,
                              const double* een_rescaled_single_e_gl,
                              const double* tmp,
                              double* const delta_een_g );

qmckl_exit_code

qmckl_compute_jastrow_champ_factor_single_een_g (const qmckl_context context,
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
                              const double* een_rescaled_e_gl,
                              const double* een_rescaled_single_e_gl,
                              const double* tmp,
                          double* const delta_een_g );

/* Provide                                                       :noexport: */


qmckl_exit_code qmckl_provide_jastrow_champ_single_een_pderiv(qmckl_context context);

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_single_een_pderiv_doc (const qmckl_context context,
                              const int64_t num,
                              const int64_t walk_num,
                              const int64_t elec_num,
                              const int64_t nucl_num,
                              const int64_t type_nucl_num,
                              const int64_t* type_nucl_vector,
                              const int64_t cord_num,
                              const int64_t dim_c_vector,
                              const double* c_vector_full,
                              const int64_t* lkpm_combined_index,
                              const double* een_rescaled_n,
                              const double* een_rescaled_e,
                              const double* een_rescaled_single_n,
                              const double* een_rescaled_single_e,
                              double* const delta_een_pderiv );

qmckl_exit_code

qmckl_compute_jastrow_champ_factor_single_een_pderiv (const qmckl_context context,
                          const int64_t num,
                          const int64_t walk_num,
                          const int64_t elec_num,
                          const int64_t nucl_num,
                          const int64_t type_nucl_num,
                          const int64_t* type_nucl_vector,
                          const int64_t cord_num,
                          const int64_t dim_c_vector,
                          const double* c_vector_full,
                          const int64_t* lkpm_combined_index,
                          const double* een_rescaled_n,
                          const double* een_rescaled_e,
                          const double* een_rescaled_single_n,
                          const double* een_rescaled_single_e,
                          double* const delta_een_pderiv );

/* Provide                                                       :noexport: */


qmckl_exit_code qmckl_provide_ee_rescaled_single(qmckl_context context);

qmckl_exit_code qmckl_compute_ee_rescaled_single_doc (
          const qmckl_context context,
          const int64_t elec_num,
          const double rescale_factor_ee,
          const int64_t walk_num,
          const double* single_ee_distance,
          double* const ee_rescaled_single );

qmckl_exit_code qmckl_compute_ee_rescaled_single(
          const qmckl_context context,
          const int64_t elec_num,
          const double rescale_factor_ee,
          const int64_t walk_num,
          const double* single_ee_distance,
          double* const ee_rescaled_single );

/* Provide                                                       :noexport: */

qmckl_exit_code qmckl_provide_jastrow_champ_single_ee(qmckl_context context);

qmckl_exit_code
qmckl_compute_jastrow_champ_single_ee_doc (const qmckl_context context,
                                           const int64_t num,
                                           const int64_t walk_num,
                                           const int64_t elec_num,
                                           const int64_t up_num,
                                           const int64_t bord_num,
                                           const double* b_vector,
                                           const double* ee_distance_rescaled,
                                           const double* ee_rescaled_singe,
                                           const int32_t spin_independent,
                                           double* const delta_ee );

qmckl_exit_code
qmckl_compute_jastrow_champ_single_ee (const qmckl_context context,
                                        const int64_t num,
                                       const int64_t walk_num,
                                       const int64_t elec_num,
                                       const int64_t up_num,
                                       const int64_t bord_num,
                                       const double* b_vector,
                                       const double* ee_distance_rescaled,
                                       const double* ee_rescaled_single,
                                       const int32_t spin_independent,
                                       double* const delta_ee );

/* Provide                                                       :noexport: */


qmckl_exit_code qmckl_provide_ee_rescaled_single_gl(qmckl_context context);

qmckl_exit_code qmckl_compute_ee_rescaled_single_gl_doc (
          const qmckl_context context,
          const int64_t num,
          const int64_t elec_num,
          const double rescale_factor_ee,
          const int64_t walk_num,
          const double* single_ee_distance,
          const double* elec_coord,
          const double* coord,
          double* const ee_rescaled_single_gl );


qmckl_exit_code qmckl_compute_ee_rescaled_single_gl (
          const qmckl_context context,
          const int64_t num,
          const int64_t elec_num,
          const double rescale_factor_ee,
          const int64_t walk_num,
          const double* single_ee_distance,
          const double* elec_coord,
          const double* coord,
          double* const ee_rescaled_single_gl );

/* Provide                                                       :noexport: */

qmckl_exit_code qmckl_provide_jastrow_champ_single_ee_gl(qmckl_context context);

qmckl_exit_code
qmckl_compute_jastrow_champ_single_ee_gl (const qmckl_context context,
                                          const int64_t num,
                                          const int64_t walk_num,
                                          const int64_t elec_num,
                                          const int64_t up_num,
                                          const int64_t bord_num,
                                          const double* b_vector,
                                          const double* ee_distance_rescaled,
                                          const double* ee_distance_rescaled_gl,
                                          const double* ee_rescaled_single,
                                          const double* ee_rescaled_single_gl,
                                          const int32_t spin_independent,
                                          double* const delta_ee_gl );

qmckl_exit_code
qmckl_compute_jastrow_champ_single_ee_gl_doc (const qmckl_context context,
                                              const int64_t num,
                                              const int64_t walk_num,
                                              const int64_t elec_num,
                                              const int64_t up_num,
                                              const int64_t bord_num,
                                              const double* b_vector,
                                              const double* ee_distance_rescaled,
                                              const double* ee_distance_rescaled_gl,
                                              const double* ee_rescaled_single,
                                              const double* ee_rescaled_single_gl,
                                              const int32_t spin_independent,
                                              double* const delta_ee_gl );

/* Provide                                                     :noexport: */

qmckl_exit_code qmckl_provide_jastrow_champ_single_ee_pderiv(qmckl_context context);

qmckl_exit_code
qmckl_compute_jastrow_champ_single_ee_pderiv_doc (const qmckl_context context,
                                           const int64_t num,
                                           const int64_t walk_num,
                                           const int64_t elec_num,
                                           const int64_t up_num,
                                           const int64_t bord_num,
                                           const double* b_vector,
                                           const double* ee_distance_rescaled,
                                           const double* ee_rescaled_singe,
                                           const int32_t spin_independent,
                                           double* const delta_ee_pderiv );

qmckl_exit_code
qmckl_compute_jastrow_champ_single_ee_pderiv (const qmckl_context context,
                                       const int64_t num,
                                       const int64_t walk_num,
                                       const int64_t elec_num,
                                       const int64_t up_num,
                                       const int64_t bord_num,
                                       const double* b_vector,
                                       const double* ee_distance_rescaled,
                                       const double* ee_rescaled_single,
                                       const int32_t spin_independent,
                                       double* const delta_ee_pderiv );

/* Provide                                                       :noexport: */


qmckl_exit_code qmckl_provide_en_rescaled_single(qmckl_context context);

qmckl_exit_code qmckl_compute_en_rescaled_single_doc (
          const qmckl_context context,
          const int64_t nucl_num,
          const int64_t type_nucl_num,
          const int64_t* type_nucl_vector,
          const double*  rescale_factor_en,
          const int64_t walk_num,
          const double* single_en_distance,
          double* const en_rescaled_single );

qmckl_exit_code qmckl_compute_en_rescaled_single (
          const qmckl_context context,
          const int64_t nucl_num,
          const int64_t type_nucl_num,
          const int64_t* type_nucl_vector,
          const double*  rescale_factor_en,
          const int64_t walk_num,
          const double* single_en_distance,
          double* const en_rescaled_single );

/* Provide                                                       :noexport: */

qmckl_exit_code qmckl_provide_jastrow_champ_single_en(qmckl_context context);

/* Provide                                       :noexport: */


qmckl_exit_code qmckl_provide_en_rescaled_single_gl(qmckl_context context);

qmckl_exit_code qmckl_compute_en_rescaled_single_gl_doc (
          const qmckl_context context,
          const int64_t nucl_num,
          const int64_t type_nucl_num,
          int64_t* const type_nucl_vector,
          const double*  rescale_factor_en,
          const int64_t walk_num,
          const double* single_en_distance,
          const double* coord,
          const double* nucl_coord,
          double* const en_rescaled_single_gl );

qmckl_exit_code qmckl_compute_en_rescaled_single_gl (
          const qmckl_context context,
          const int64_t nucl_num,
          const int64_t type_nucl_num,
          int64_t* const type_nucl_vector,
          const double*  rescale_factor_en,
          const int64_t walk_num,
          const double* single_en_distance,
          const double* coord,
          const double* nucl_coord,
          double* const en_rescaled_single_gl );

/* Provide                                                       :noexport: */

qmckl_exit_code qmckl_provide_jastrow_champ_single_en_gl(qmckl_context context);

qmckl_exit_code qmckl_compute_jastrow_champ_single_en_gl_doc (
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
          double* const delta_en_gl );

qmckl_exit_code qmckl_compute_jastrow_champ_single_en_gl (
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
          double* const delta_en_gl );

/* Provide                                                  :noexport: */

qmckl_exit_code qmckl_provide_jastrow_champ_single_en_pderiv(qmckl_context context);

#endif
