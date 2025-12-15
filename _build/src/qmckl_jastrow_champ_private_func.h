#ifndef QMCKL_JASTROW_CHAMP_HPF
#define QMCKL_JASTROW_CHAMP_HPF
#include "qmckl_jastrow_champ_private_type.h"



/* The ~uninitialized~ integer contains one bit set to one for each */
/* initialization function which has not been called. It becomes equal */
/* to zero after all initialization functions have been called. The */
/* struct is then initialized and ~provided == true~. */
/* Some values are initialized by default, and are not concerned by */
/* this mechanism. */


qmckl_exit_code qmckl_init_jastrow_champ(qmckl_context context);



/* When the required information is completely entered, other data structures are */
/* computed to accelerate the calculations. The intermediates factors */
/* are precontracted using BLAS LEVEL 3 operations. */


qmckl_exit_code qmckl_finalize_jastrow_champ(qmckl_context context);

/* Deep copy */


qmckl_exit_code qmckl_copy_jastrow_champ(qmckl_context context, const qmckl_jastrow_champ_struct* src, qmckl_jastrow_champ_struct* dest);

/* Provide                                                       :noexport: */

qmckl_exit_code qmckl_provide_jastrow_champ_asymp_jasb(qmckl_context context);

qmckl_exit_code qmckl_compute_jastrow_champ_asymp_jasb_doc (const qmckl_context context,
                                                            const int64_t bord_num,
                                                            const double* b_vector,
                                                            const double rescale_factor_ee,
                                                            const int32_t spin_independent,
                                                            double* const asymp_jasb);

qmckl_exit_code qmckl_compute_jastrow_champ_asymp_jasb_hpc (const qmckl_context context,
                                                            const int64_t bord_num,
                                                            const double* b_vector,
                                                            const double rescale_factor_ee,
                                                            const int32_t spin_independent,
                                                            double* const asymp_jasb );

qmckl_exit_code qmckl_compute_jastrow_champ_asymp_jasb (const qmckl_context context,
                                                        const int64_t bord_num,
                                                        const double* b_vector,
                                                        const double rescale_factor_ee,
                                                        const int32_t spin_independent,
                                                        double* const asymp_jasb );

/* Provide                                                       :noexport: */


qmckl_exit_code qmckl_provide_ee_distance_rescaled(qmckl_context context);

qmckl_exit_code qmckl_compute_ee_distance_rescaled_doc (
          const qmckl_context context,
          const int64_t elec_num,
          const double rescale_factor_ee,
          const int64_t walk_num,
          const double* coord,
          double* const ee_distance_rescaled );

qmckl_exit_code qmckl_compute_ee_distance_rescaled_hpc (
          const qmckl_context context,
          const int64_t elec_num,
          const double rescale_factor_ee,
          const int64_t walk_num,
          const double* coord,
          double* const ee_distance_rescaled );

qmckl_exit_code qmckl_compute_ee_distance_rescaled (
          const qmckl_context context,
          const int64_t elec_num,
          const double rescale_factor_ee,
          const int64_t walk_num,
          const double* coord,
          double* const ee_distance_rescaled );

/* Provide                                                       :noexport: */


qmckl_exit_code qmckl_provide_ee_distance_rescaled_gl(qmckl_context context);

qmckl_exit_code qmckl_compute_ee_distance_rescaled_gl_doc (
          const qmckl_context context,
          const int64_t elec_num,
          const double rescale_factor_ee,
          const int64_t walk_num,
          const double* coord,
          double* const ee_distance_rescaled_gl );

qmckl_exit_code qmckl_compute_ee_distance_rescaled_gl_hpc (
          const qmckl_context context,
          const int64_t elec_num,
          const double rescale_factor_ee,
          const int64_t walk_num,
          const double* coord,
          double* const ee_distance_rescaled_gl );

qmckl_exit_code qmckl_compute_ee_distance_rescaled_gl (
          const qmckl_context context,
          const int64_t elec_num,
          const double rescale_factor_ee,
          const int64_t walk_num,
          const double* coord,
          double* const ee_distance_rescaled_gl );

/* Provide                                                       :noexport: */

qmckl_exit_code qmckl_provide_jastrow_champ_factor_ee(qmckl_context context);

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_ee_doc (const qmckl_context context,
                                           const int64_t walk_num,
                                           const int64_t elec_num,
                                           const int64_t up_num,
                                           const int64_t bord_num,
                                           const double* b_vector,
                                           const double* ee_distance_rescaled,
                                           const double* asymp_jasb,
                                           const int32_t spin_independent,
                                           double* const factor_ee );

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_ee_hpc (const qmckl_context context,
                                           const int64_t walk_num,
                                           const int64_t elec_num,
                                           const int64_t up_num,
                                           const int64_t bord_num,
                                           const double* b_vector,
                                           const double* ee_distance_rescaled,
                                           const double* asymp_jasb,
                                           const int32_t spin_independent,
                                           double* const factor_ee );

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_ee (const qmckl_context context,
                                       const int64_t walk_num,
                                       const int64_t elec_num,
                                       const int64_t up_num,
                                       const int64_t bord_num,
                                       const double* b_vector,
                                       const double* ee_distance_rescaled,
                                       const double* asymp_jasb,
                                       const int32_t spin_independent,
                                       double* const factor_ee );

/* Provide                                                       :noexport: */

qmckl_exit_code qmckl_provide_jastrow_champ_factor_ee_gl(qmckl_context context);

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_ee_gl (const qmckl_context context,
                                          const int64_t walk_num,
                                          const int64_t elec_num,
                                          const int64_t up_num,
                                          const int64_t bord_num,
                                          const double* b_vector,
                                          const double* ee_distance_rescaled,
                                          const double* ee_distance_rescaled_gl,
                                          const int32_t spin_independent,
                                          double* const factor_ee_gl );

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_ee_gl_hpc (const qmckl_context context,
                                              const int64_t walk_num,
                                              const int64_t elec_num,
                                              const int64_t up_num,
                                              const int64_t bord_num,
                                              const double* b_vector,
                                              const double* ee_distance_rescaled,
                                              const double* ee_distance_rescaled_gl,
                                              const int32_t spin_independent,
                                              double* const factor_ee_gl );

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_ee_gl_doc (const qmckl_context context,
                                              const int64_t walk_num,
                                              const int64_t elec_num,
                                              const int64_t up_num,
                                              const int64_t bord_num,
                                              const double* b_vector,
                                              const double* ee_distance_rescaled,
                                              const double* ee_distance_rescaled_gl,
                                              const int32_t spin_independent,
                                              double* const factor_ee_gl );

/* Provide                                                       :noexport: */

qmckl_exit_code qmckl_provide_jastrow_champ_asymp_jasb_pderiv(qmckl_context context);

qmckl_exit_code qmckl_compute_jastrow_champ_asymp_jasb_pderiv_doc (const qmckl_context context,
                                                                   const int64_t bord_num,
                                                                   const double* b_vector,
                                                                   const double rescale_factor_ee,
                                                                   const int32_t spin_independent,
                                                                   double* const asymp_jasb_pderiv);

/* Provide                                                       :noexport: */


qmckl_exit_code qmckl_provide_jastrow_champ_factor_ee_pderiv(qmckl_context context);

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_ee_pderiv (const qmckl_context context,
                                              const int64_t walk_num,
                                              const int64_t elec_num,
                                              const int64_t up_num,
                                              const int64_t bord_num,
                                              const double* b_vector,
                                              const double* ee_distance_rescaled,
                                              const double* asymp_jasb_pderiv,
                                              const int32_t spin_independent,
                                              double* const factor_ee_pderiv );

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_ee_pderiv_doc (const qmckl_context context,
                                                  const int64_t walk_num,
                                                  const int64_t elec_num,
                                                  const int64_t up_num,
                                                  const int64_t bord_num,
                                                  const double* b_vector,
                                                  const double* ee_distance_rescaled,
                                                  const double* asymp_jasb_pderiv,
                                                  const int32_t spin_independent,
                                                  double* const factor_ee_pderiv );

/* Provide                                                       :noexport: */


qmckl_exit_code qmckl_provide_jastrow_champ_factor_ee_gl_pderiv(qmckl_context context);

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_ee_gl_pderiv (const qmckl_context context,
                                          const int64_t walk_num,
                                          const int64_t elec_num,
                                          const int64_t up_num,
                                          const int64_t bord_num,
                                          const double* b_vector,
                                          const double* ee_distance_rescaled,
                                          const double* ee_distance_rescaled_gl,
                                          const int32_t spin_independent,
                                          double* const factor_ee_gl_pderiv );

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_ee_gl_pderiv_doc (const qmckl_context context,
                                              const int64_t walk_num,
                                              const int64_t elec_num,
                                              const int64_t up_num,
                                              const int64_t bord_num,
                                              const double* b_vector,
                                              const double* ee_distance_rescaled,
                                              const double* ee_distance_rescaled_gl,
                                              const int32_t spin_independent,
                                              double* const factor_ee_gl_pderiv );

/* Provide                                                       :noexport: */

qmckl_exit_code qmckl_provide_jastrow_champ_asymp_jasa(qmckl_context context);

/* Provide                                                       :noexport: */


qmckl_exit_code qmckl_provide_en_distance_rescaled(qmckl_context context);

qmckl_exit_code qmckl_compute_en_distance_rescaled_doc (
          const qmckl_context context,
          const int64_t elec_num,
          const int64_t nucl_num,
          const int64_t type_nucl_num,
          const int64_t* type_nucl_vector,
          const double*  rescale_factor_en,
          const int64_t walk_num,
          const double* elec_coord,
          const double* nucl_coord,
          double* const en_distance_rescaled );

qmckl_exit_code qmckl_compute_en_distance_rescaled_hpc (
          const qmckl_context context,
          const int64_t elec_num,
          const int64_t nucl_num,
          const int64_t type_nucl_num,
          const int64_t* type_nucl_vector,
          const double* rescale_factor_en,
          const int64_t walk_num,
          const double* elec_coord,
          const double* nucl_coord,
          double* const en_distance_rescaled );

qmckl_exit_code qmckl_compute_en_distance_rescaled (
          const qmckl_context context,
          const int64_t elec_num,
          const int64_t nucl_num,
          const int64_t type_nucl_num,
          const int64_t* type_nucl_vector,
          const double*  rescale_factor_en,
          const int64_t walk_num,
          const double* elec_coord,
          const double* nucl_coord,
          double* const en_distance_rescaled );

/* Provide                                                       :noexport: */


qmckl_exit_code qmckl_provide_en_distance_rescaled_gl(qmckl_context context);

qmckl_exit_code qmckl_compute_en_distance_rescaled_gl_doc (
          const qmckl_context context,
          const int64_t elec_num,
          const int64_t nucl_num,
          const int64_t type_nucl_num,
          int64_t* const type_nucl_vector,
          const double*  rescale_factor_en,
          const int64_t walk_num,
          const double* elec_coord,
          const double* nucl_coord,
          double* const en_distance_rescaled_gl );

qmckl_exit_code qmckl_compute_en_distance_rescaled_gl_hpc (
          const qmckl_context context,
          const int64_t elec_num,
          const int64_t nucl_num,
          const int64_t type_nucl_num,
          int64_t* const type_nucl_vector,
          const double*  rescale_factor_en,
          const int64_t walk_num,
          const double* elec_coord,
          const double* nucl_coord,
          double* const en_distance_rescaled_gl );

qmckl_exit_code qmckl_compute_en_distance_rescaled_gl (
          const qmckl_context context,
          const int64_t elec_num,
          const int64_t nucl_num,
          const int64_t type_nucl_num,
          int64_t* const type_nucl_vector,
          const double*  rescale_factor_en,
          const int64_t walk_num,
          const double* elec_coord,
          const double* nucl_coord,
          double* const en_distance_rescaled_gl );

/* Provide                                                       :noexport: */

qmckl_exit_code qmckl_provide_jastrow_champ_factor_en(qmckl_context context);

/* Provide                                                       :noexport: */

qmckl_exit_code qmckl_provide_jastrow_champ_factor_en_gl(qmckl_context context);



/* #   #+CALL: generate_private_c_header(table=qmckl_factor_en_gl_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */


qmckl_exit_code qmckl_compute_jastrow_champ_factor_en_gl_doc (
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
          double* const factor_en_gl );

qmckl_exit_code qmckl_compute_jastrow_champ_factor_en_gl_hpc (
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
          double* const factor_en_gl );

qmckl_exit_code qmckl_compute_jastrow_champ_factor_en_gl (
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
          double* const factor_en_gl );

/* Provide                                                       :noexport: */

qmckl_exit_code qmckl_provide_jastrow_champ_asymp_jasa_pderiv(qmckl_context context);

/* Provide                                                       :noexport: */

qmckl_exit_code qmckl_provide_jastrow_champ_factor_en_pderiv(qmckl_context context);

/* Provide                                                       :noexport: */

qmckl_exit_code qmckl_provide_jastrow_champ_factor_en_gl_pderiv(qmckl_context context);

qmckl_exit_code qmckl_compute_jastrow_champ_factor_en_gl_pderiv_doc (
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
          double* const factor_en_gl_pderiv );


qmckl_exit_code qmckl_compute_jastrow_champ_factor_en_gl_pderiv (
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
          double* const factor_en_gl_pderiv );

/* Provide                                                       :noexport: */


qmckl_exit_code qmckl_provide_een_rescaled_e(qmckl_context context);



/* #   #+CALL: generate_private_c_header(table=qmckl_factor_een_rescaled_e_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */


qmckl_exit_code qmckl_compute_een_rescaled_e (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t cord_num,
      const double rescale_factor_ee,
      const double* ee_distance,
      double* const een_rescaled_e );



/* #   #+CALL: generate_private_c_header(table=qmckl_factor_een_rescaled_e_args,rettyp=get_value("CRetType"),fname="qmckl_compute_een_rescaled_e_doc") */


qmckl_exit_code qmckl_compute_een_rescaled_e (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t cord_num,
      const double rescale_factor_ee,
      const double* ee_distance,
      double* const een_rescaled_e );

qmckl_exit_code qmckl_compute_een_rescaled_e_doc (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t cord_num,
      const double rescale_factor_ee,
      const double* ee_distance,
      double* const een_rescaled_e );

qmckl_exit_code qmckl_compute_een_rescaled_e_hpc (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t cord_num,
      const double rescale_factor_ee,
      const double* ee_distance,
      double* const een_rescaled_e );

/* Provide                                                       :noexport: */


qmckl_exit_code qmckl_provide_een_rescaled_e_gl(qmckl_context context);



/* #  #+CALL: generate_private_c_header(table=qmckl_factor_een_rescaled_e_gl_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */


qmckl_exit_code qmckl_compute_jastrow_champ_factor_een_rescaled_e_gl (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t cord_num,
      const double rescale_factor_ee,
      const double* coord_ee,
      const double* ee_distance,
      const double* een_rescaled_e,
      double* const een_rescaled_e_gl );

qmckl_exit_code qmckl_compute_jastrow_champ_factor_een_rescaled_e_gl_doc (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t cord_num,
      const double rescale_factor_ee,
      const double* coord_ee,
      const double* ee_distance,
      const double* een_rescaled_e,
      double* const een_rescaled_e_gl );

qmckl_exit_code qmckl_compute_jastrow_champ_factor_een_rescaled_e_gl_hpc (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t cord_num,
      const double rescale_factor_ee,
      const double* coord_ee,
      const double* ee_distance,
      const double* een_rescaled_e,
      double* const een_rescaled_e_gl );

/* Provide                                                       :noexport: */


qmckl_exit_code qmckl_provide_een_rescaled_n(qmckl_context context);



/* #   #+CALL: generate_private_c_header(table=qmckl_factor_een_rescaled_n_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */


qmckl_exit_code qmckl_compute_een_rescaled_n (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t type_nucl_num,
      int64_t* const type_nucl_vector,
      const int64_t cord_num,
      const double* rescale_factor_en,
      const double* en_distance,
      double* const een_rescaled_n );

/* Provide                                                       :noexport: */


qmckl_exit_code qmckl_provide_een_rescaled_n_gl(qmckl_context context);



/* #   #+CALL: generate_private_c_header(table=qmckl_compute_jastrow_champ_factor_een_rescaled_n_gl_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */


qmckl_exit_code qmckl_compute_jastrow_champ_factor_een_rescaled_n_gl (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t type_nucl_num,
      int64_t* const type_nucl_vector,
      const int64_t cord_num,
      const double* rescale_factor_en,
      const double* coord_ee,
      const double* coord_n,
      const double* en_distance,
      const double* een_rescaled_n,
      double* const een_rescaled_n_gl );




/* #   #+CALL: generate_private_c_header(table=qmckl_factor_dim_c_vector_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */


qmckl_exit_code qmckl_compute_dim_c_vector(
      const qmckl_context context,
      const int64_t cord_num,
      int64_t* const dim_c_vector );

qmckl_exit_code qmckl_compute_dim_c_vector_doc(
      const qmckl_context context,
      const int64_t cord_num,
      int64_t* const dim_c_vector );

qmckl_exit_code qmckl_compute_dim_c_vector_hpc(
      const qmckl_context context,
      const int64_t cord_num,
      int64_t* const dim_c_vector );

/* Get */


qmckl_exit_code
qmckl_get_jastrow_champ_tmp_c(qmckl_context context,
                              double* const tmp_c,
                              const int64_t size_max);

qmckl_exit_code
qmckl_get_jastrow_champ_dtmp_c(qmckl_context context,
                               double* const dtmp_c,
                               const int64_t size_max);

/* Provide                                                       :noexport: */


qmckl_exit_code qmckl_provide_jastrow_champ_c_vector_full(qmckl_context context);
qmckl_exit_code qmckl_provide_lkpm_combined_index(qmckl_context context);
qmckl_exit_code qmckl_provide_tmp_c(qmckl_context context);
qmckl_exit_code qmckl_provide_dtmp_c(qmckl_context context);




/* #   #+CALL: generate_private_c_header(table=qmckl_factor_c_vector_full_args,rettyp=get_value("CRetType"),fname="qmckl_compute_c_vector_full_doc") */


qmckl_exit_code qmckl_compute_c_vector_full (
      const qmckl_context context,
      const int64_t nucl_num,
      const int64_t dim_c_vector,
      const int64_t type_nucl_num,
      const int64_t* type_nucl_vector,
      const double* c_vector,
      double* const c_vector_full );

qmckl_exit_code qmckl_compute_c_vector_full_doc (
      const qmckl_context context,
      const int64_t nucl_num,
      const int64_t dim_c_vector,
      const int64_t type_nucl_num,
      const int64_t* type_nucl_vector,
      const double* c_vector,
      double* const c_vector_full );

qmckl_exit_code qmckl_compute_c_vector_full_hpc (
      const qmckl_context context,
      const int64_t nucl_num,
      const int64_t dim_c_vector,
      const int64_t type_nucl_num,
      const int64_t* type_nucl_vector,
      const double* c_vector,
      double* const c_vector_full );




/* #+CALL: generate_private_c_header(table=lkpm_combined_index_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */


qmckl_exit_code qmckl_compute_lkpm_combined_index (
      const qmckl_context context,
      const int64_t cord_num,
      const int64_t dim_c_vector,
      int64_t* const lkpm_combined_index );

qmckl_exit_code qmckl_compute_lkpm_combined_index (
      const qmckl_context context,
      const int64_t cord_num,
      const int64_t dim_c_vector,
      int64_t* const lkpm_combined_index );


/* #+CALL: generate_private_c_header(table=lkpm_combined_index_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */


qmckl_exit_code qmckl_compute_lkpm_combined_index (
      const qmckl_context context,
      const int64_t cord_num,
      const int64_t dim_c_vector,
      int64_t* const lkpm_combined_index );



/* #   #+CALL: generate_private_c_header(table=qmckl_factor_tmp_c_args,rettyp=get_value("CRetType"),fname="qmckl_compute_tmp_c") */


qmckl_exit_code qmckl_compute_tmp_c (
      const qmckl_context context,
      const int64_t cord_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t walk_num,
      const double* een_rescaled_e,
      const double* een_rescaled_n,
      double* const tmp_c );

qmckl_exit_code qmckl_compute_tmp_c_doc (
          const qmckl_context context,
          const int64_t cord_num,
          const int64_t elec_num,
          const int64_t nucl_num,
          const int64_t walk_num,
          const double* een_rescaled_e,
          const double* een_rescaled_n,
          double* const tmp_c );



/* #   #+CALL: generate_private_c_header(table=qmckl_factor_tmp_c_args,rettyp=get_value("CRetType"),fname="qmckl_compute_tmp_c_doc") */

/*     #+RESULTS: */

qmckl_exit_code qmckl_compute_tmp_c_doc (
          const qmckl_context context,
          const int64_t cord_num,
          const int64_t elec_num,
          const int64_t nucl_num,
          const int64_t walk_num,
          const double* een_rescaled_e,
          const double* een_rescaled_n,
          double* const tmp_c );



/* #   #+CALL: generate_private_c_header(table=qmckl_factor_tmp_c_args,rettyp=get_value("CRetType"),fname="qmckl_compute_tmp_c_hpc") */

/*     #+RESULTS: */


qmckl_exit_code qmckl_compute_tmp_c_hpc (const qmckl_context context,
                                         const int64_t cord_num,
                                         const int64_t elec_num,
                                         const int64_t nucl_num,
                                         const int64_t walk_num,
                                         const double* een_rescaled_e,
                                         const double* een_rescaled_n,
                                         double* const tmp_c );

/* Compute dtmp_c */
/*      :PROPERTIES: */
/*      :Name:     qmckl_compute_dtmp_c */
/*      :CRetType: qmckl_exit_code */
/*      :FRetType: qmckl_exit_code */
/*      :END: */

/*      #+NAME: qmckl_factor_dtmp_c_args */
/*      |---------------------+---------------------------------------------------------------------+--------+-----------------------------------------------| */
/*      | Variable            | Type                                                                | In/Out | Description                                   | */
/*      |---------------------+---------------------------------------------------------------------+--------+-----------------------------------------------| */
/*      | ~context~           | ~qmckl_context~                                                     | in     | Global state                                  | */
/*      | ~cord_num~          | ~int64_t~                                                           | in     | Order of polynomials                          | */
/*      | ~elec_num~          | ~int64_t~                                                           | in     | Number of electrons                           | */
/*      | ~nucl_num~          | ~int64_t~                                                           | in     | Number of nuclei                              | */
/*      | ~walk_num~          | ~int64_t~                                                           | in     | Number of walkers                             | */
/*      | ~een_rescaled_e_gl~ | ~double[walk_num][0:cord_num][elec_num][4][elec_num]~               | in     | Electron-electron rescaled factor derivatives | */
/*      | ~een_rescaled_n~    | ~double[walk_num][0:cord_num][nucl_num][elec_num]~                  | in     | Electron-nucleus rescaled factor              | */
/*      | ~dtmp_c~            | ~double[walk_num][0:cord_num-1][0:cord_num][nucl_num][4][elec_num]~ | out    | vector of non-zero coefficients               | */
/*      |---------------------+---------------------------------------------------------------------+--------+-----------------------------------------------| */



qmckl_exit_code
qmckl_compute_dtmp_c (const qmckl_context context,
                      const int64_t cord_num,
                      const int64_t elec_num,
                      const int64_t nucl_num,
                      const int64_t walk_num,
                      const double* een_rescaled_e_gl,
                      const double* een_rescaled_n,
                      double* const dtmp_c );

qmckl_exit_code qmckl_compute_dtmp_c_doc (
      const qmckl_context context,
      const int64_t cord_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t walk_num,
      const double* een_rescaled_e_gl,
      const double* een_rescaled_n,
      double* const dtmp_c );

qmckl_exit_code qmckl_compute_dtmp_c_hpc (
      const qmckl_context context,
      const int64_t cord_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t walk_num,
      const double* een_rescaled_e_gl,
      const double* een_rescaled_n,
      double* const dtmp_c );

/* Provide                                                       :noexport: */

qmckl_exit_code qmckl_provide_jastrow_champ_factor_een(qmckl_context context);



/* #   #+CALL: generate_private_c_header(table=qmckl_factor_een_naive_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */


qmckl_exit_code qmckl_compute_jastrow_champ_factor_een_naive (
  const qmckl_context context,
  const int64_t walk_num,
  const int64_t elec_num,
  const int64_t nucl_num,
  const int64_t cord_num,
  const int64_t dim_c_vector,
  const double* c_vector_full,
  const int64_t* lkpm_combined_index,
  const double* een_rescaled_e,
  const double* een_rescaled_n,
  double* const factor_een );


/* # #+CALL: generate_private_c_header(table=qmckl_factor_een_args,rettyp=qmckl_exit_code),fname=get_value("Name")) */


qmckl_exit_code
qmckl_compute_jastrow_champ_factor_een_doc (const qmckl_context context,
                              const int64_t walk_num,
                              const int64_t elec_num,
                              const int64_t nucl_num,
                              const int64_t cord_num,
                              const int64_t dim_c_vector,
                              const double* c_vector_full,
                              const int64_t* lkpm_combined_index,
                              const double* tmp_c,
                              const double* een_rescaled_n,
                              double* const factor_een );

qmckl_exit_code

qmckl_compute_jastrow_champ_factor_een (const qmckl_context context,
                          const int64_t walk_num,
                          const int64_t elec_num,
                          const int64_t nucl_num,
                          const int64_t cord_num,
                          const int64_t dim_c_vector,
                          const double* c_vector_full,
                          const int64_t* lkpm_combined_index,
                          const double* tmp_c,
                          const double* een_rescaled_n,
                          double* const factor_een );

/* Provide                                                       :noexport: */

qmckl_exit_code qmckl_provide_jastrow_champ_factor_een_gl(qmckl_context context);
qmckl_exit_code qmckl_provide_jastrow_champ_factor_een_grad(qmckl_context context);



/* #   #+CALL: generate_private_c_header(table=qmckl_factor_een_gl_naive_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */


qmckl_exit_code qmckl_compute_jastrow_champ_factor_een_gl_naive (
  const qmckl_context context,
  const int64_t walk_num,
  const int64_t elec_num,
  const int64_t nucl_num,
  const int64_t cord_num,
  const int64_t dim_c_vector,
  const double* c_vector_full,
  const int64_t* lkpm_combined_index,
  const double* een_rescaled_e,
  const double* een_rescaled_n,
  const double* een_rescaled_e_gl,
  const double* een_rescaled_n_gl,
  double* const factor_een_gl );



/* #+CALL: generate_private_c_header(table=qmckl_factor_een_gl_args,rettyp=get_value("CRetType"),fname="qmckl_compute_jastrow_champ_factor_een_gl_doc" ) */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_jastrow_champ_factor_een_gl_doc (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t cord_num,
      const int64_t dim_c_vector,
      const double* c_vector_full,
      const int64_t* lkpm_combined_index,
      const double* tmp_c,
      const double* dtmp_c,
      const double* een_rescaled_n,
      const double* een_rescaled_n_gl,
      double* const factor_een_gl );



/* #+CALL: generate_private_c_header(table=qmckl_factor_een_gl_args,rettyp=get_value("CRetType"),fname="qmckl_compute_jastrow_champ_factor_een_gl" ) */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_jastrow_champ_factor_een_gl (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t cord_num,
      const int64_t dim_c_vector,
      const double* c_vector_full,
      const int64_t* lkpm_combined_index,
      const double* tmp_c,
      const double* dtmp_c,
      const double* een_rescaled_n,
      const double* een_rescaled_n_gl,
      double* const factor_een_gl );

/* HPC implementation                                           :noexport: */
/*    #+CALL: generate_private_c_header(table=qmckl_factor_een_gl_args,rettyp=get_value("CRetType"),fname="qmckl_compute_jastrow_champ_factor_een_gl_hpc" ) */

/*    #+RESULTS: */

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_een_gl_hpc (
         const qmckl_context context,
         const int64_t walk_num,
         const int64_t elec_num,
         const int64_t nucl_num,
         const int64_t cord_num,
         const int64_t dim_c_vector,
         const double* c_vector_full,
         const int64_t* lkpm_combined_index,
         const double* tmp_c,
         const double* dtmp_c,
         const double* een_rescaled_n,
         const double* een_rescaled_n_gl,
         double* const factor_een_gl );



/* #+CALL: generate_private_c_header(table=qmckl_factor_een_grad_args,rettyp=get_value("CRetType"),fname="qmckl_compute_jastrow_champ_factor_een_grad_doc" ) */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_jastrow_champ_factor_een_grad_doc (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t cord_num,
      const int64_t dim_c_vector,
      const double* c_vector_full,
      const int64_t* lkpm_combined_index,
      const double* tmp_c,
      const double* dtmp_c,
      const double* een_rescaled_n,
      const double* een_rescaled_n_gl,
      double* const factor_een_grad );



/* #+CALL: generate_private_c_header(table=qmckl_factor_een_grad_args,rettyp=get_value("CRetType"),fname="qmckl_compute_jastrow_champ_factor_een_grad" ) */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_jastrow_champ_factor_een_grad (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t cord_num,
      const int64_t dim_c_vector,
      const double* c_vector_full,
      const int64_t* lkpm_combined_index,
      const double* tmp_c,
      const double* dtmp_c,
      const double* een_rescaled_n,
      const double* een_rescaled_n_gl,
      double* const factor_een_grad );

/* HPC implementation                                           :noexport: */
/*    #+CALL: generate_private_c_header(table=qmckl_factor_een_grad_args,rettyp=get_value("CRetType"),fname="qmckl_compute_jastrow_champ_factor_een_grad_hpc" ) */

/*    #+RESULTS: */

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_een_grad_hpc (
         const qmckl_context context,
         const int64_t walk_num,
         const int64_t elec_num,
         const int64_t nucl_num,
         const int64_t cord_num,
         const int64_t dim_c_vector,
         const double* c_vector_full,
         const int64_t* lkpm_combined_index,
         const double* tmp_c,
         const double* dtmp_c,
         const double* een_rescaled_n,
         const double* een_rescaled_n_gl,
         double* const factor_een_grad );

/* Provide                                                       :noexport: */

qmckl_exit_code qmckl_provide_jastrow_champ_factor_een_pderiv(qmckl_context context);

qmckl_exit_code
qmckl_compute_jastrow_champ_factor_een_pderiv_doc (const qmckl_context context,
                              const int64_t walk_num,
                              const int64_t elec_num,
                              const int64_t nucl_num,
                              const int64_t type_nucl_num,
                              const int64_t* type_nucl_vector,
                              const int64_t cord_num,
                              const int64_t dim_c_vector,
                              const double* c_vector_full,
                              const int64_t* lkpm_combined_index,
                              const double* tmp_c,
                              const double* een_rescaled_n,
                              double* const factor_een_pderiv );

qmckl_exit_code

qmckl_compute_jastrow_champ_factor_een_pderiv (const qmckl_context context,
                          const int64_t walk_num,
                          const int64_t elec_num,
                          const int64_t nucl_num,
                          const int64_t type_nucl_num,
                          const int64_t* type_nucl_vector,
                          const int64_t cord_num,
                          const int64_t dim_c_vector,
                          const double* c_vector_full,
                          const int64_t* lkpm_combined_index,
                          const double* tmp_c,
                          const double* een_rescaled_n,
                          double* const factor_een_pderiv );

/* Provide                                                      :noexport: */

qmckl_exit_code qmckl_provide_jastrow_champ_factor_een_gl_pderiv(qmckl_context context);

qmckl_exit_code qmckl_compute_jastrow_champ_factor_een_gl_pderiv_doc (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t type_nucl_num,
      const int64_t* type_nucl_vector,
      const int64_t cord_num,
      const int64_t dim_c_vector,
      const double* c_vector_full,
      const int64_t* lkpm_combined_index,
      const double* tmp_c,
      const double* dtmp_c,
      const double* een_rescaled_n,
      const double* een_rescaled_n_gl,
      double* const factor_een_gl_pderiv );

qmckl_exit_code qmckl_compute_jastrow_champ_factor_een_gl_pderiv (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t type_nucl_num,
      const int64_t* type_nucl_vector,
      const int64_t cord_num,
      const int64_t dim_c_vector,
      const double* c_vector_full,
      const int64_t* lkpm_combined_index,
      const double* tmp_c,
      const double* dtmp_c,
      const double* een_rescaled_n,
      const double* een_rescaled_n_gl,
      double* const factor_een_gl_pderiv );

/* Provide                                                       :noexport: */

qmckl_exit_code qmckl_provide_jastrow_champ_value(qmckl_context context);



/* #+CALL: generate_private_c_header(table=qmckl_jastrow_champ_value_args,rettyp=get_value("CRetType"),fname="qmckl_compute_jastrow_champ_value") */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_jastrow_champ_value (
      const qmckl_context context,
      const int64_t walk_num,
      const double* f_ee,
      const double* f_en,
      const double* f_een,
      double* const value );



/* #+CALL: generate_private_c_header(table=qmckl_jastrow_champ_value_args,rettyp=get_value("CRetType"),fname="qmckl_compute_jastrow_champ_value_doc") */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_jastrow_champ_value_doc (
      const qmckl_context context,
      const int64_t walk_num,
      const double* f_ee,
      const double* f_en,
      const double* f_een,
      double* const value );



/* #+CALL: generate_private_c_header(table=qmckl_jastrow_champ_value_args,rettyp=get_value("CRetType"),fname="qmckl_compute_jastrow_champ_value_hpc") */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_jastrow_champ_value_hpc (
      const qmckl_context context,
      const int64_t walk_num,
      const double* f_ee,
      const double* f_en,
      const double* f_een,
      double* const value );

/* Provide                                                       :noexport: */

qmckl_exit_code qmckl_provide_jastrow_champ_gl(qmckl_context context);
qmckl_exit_code qmckl_provide_jastrow_champ_grad(qmckl_context context);



/* #+CALL: generate_private_c_header(table=qmckl_jastrow_champ_gl_args,rettyp=get_value("CRetType"),fname="qmckl_compute_jastrow_champ_gl") */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_jastrow_champ_gl (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const double* value,
      const double* gl_ee,
      const double* gl_en,
      const double* gl_een,
      double* const gl );



/* #+CALL: generate_private_c_header(table=qmckl_jastrow_champ_gl_args,rettyp=get_value("CRetType"),fname="qmckl_compute_jastrow_champ_gl_doc") */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_jastrow_champ_gl_doc (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const double* value,
      const double* gl_ee,
      const double* gl_en,
      const double* gl_een,
      double* const gl );



/* #+CALL: generate_private_c_header(table=qmckl_jastrow_champ_gl_args,rettyp=get_value("CRetType"),fname="qmckl_compute_jastrow_champ_gl_hpc") */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_jastrow_champ_gl_hpc (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const double* value,
      const double* gl_ee,
      const double* gl_en,
      const double* gl_een,
      double* const gl );



/* #+CALL: generate_private_c_header(table=qmckl_jastrow_champ_grad_args,rettyp=get_value("CRetType"),fname="qmckl_compute_jastrow_champ_grad") */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_jastrow_champ_grad (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const double* value,
      const double* gl_ee,
      const double* gl_en,
      const double* grad_een,
      double* const grad );



/* #+CALL: generate_private_c_header(table=qmckl_jastrow_champ_grad_args,rettyp=get_value("CRetType"),fname="qmckl_compute_jastrow_champ_grad_doc") */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_jastrow_champ_grad_doc (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const double* value,
      const double* gl_ee,
      const double* gl_en,
      const double* grad_een,
      double* const grad );



/* #+CALL: generate_private_c_header(table=qmckl_jastrow_champ_grad_args,rettyp=get_value("CRetType"),fname="qmckl_compute_jastrow_champ_grad_hpc") */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_jastrow_champ_grad_hpc (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const double* value,
      const double* gl_ee,
      const double* gl_en,
      const double* grad_een,
      double* const grad );

#endif
