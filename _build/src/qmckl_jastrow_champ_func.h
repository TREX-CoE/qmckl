/* Initialization functions */

/*    To prepare for the Jastrow and its derivative, all the following functions need to be */
/*    called. */


qmckl_exit_code  qmckl_set_jastrow_champ_rescale_factor_ee (qmckl_context context, const double  kappa_ee);
qmckl_exit_code  qmckl_set_jastrow_champ_rescale_factor_en (qmckl_context context, const double* kappa_en, const int64_t size_max);
qmckl_exit_code  qmckl_set_jastrow_champ_aord_num          (qmckl_context context, const int64_t aord_num);
qmckl_exit_code  qmckl_set_jastrow_champ_bord_num          (qmckl_context context, const int64_t bord_num);
qmckl_exit_code  qmckl_set_jastrow_champ_cord_num          (qmckl_context context, const int64_t cord_num);
qmckl_exit_code  qmckl_set_jastrow_champ_type_nucl_num     (qmckl_context context, const int64_t type_nucl_num);
qmckl_exit_code  qmckl_set_jastrow_champ_type_nucl_vector  (qmckl_context context, const int64_t* type_nucl_vector, const int64_t size_max);
qmckl_exit_code  qmckl_set_jastrow_champ_a_vector          (qmckl_context context, const double * a_vector, const int64_t size_max);
qmckl_exit_code  qmckl_set_jastrow_champ_b_vector          (qmckl_context context, const double * b_vector, const int64_t size_max);
qmckl_exit_code  qmckl_set_jastrow_champ_c_vector          (qmckl_context context, const double * c_vector, const int64_t size_max);
qmckl_exit_code  qmckl_set_jastrow_champ_spin_independent  (qmckl_context context, const int32_t spin_independent);

/* Access functions */


qmckl_exit_code  qmckl_get_jastrow_champ_aord_num          (qmckl_context context, int64_t* const aord_num);
qmckl_exit_code  qmckl_get_jastrow_champ_bord_num          (qmckl_context context, int64_t* const bord_num);
qmckl_exit_code  qmckl_get_jastrow_champ_cord_num          (qmckl_context context, int64_t* const bord_num);
qmckl_exit_code  qmckl_get_jastrow_champ_type_nucl_num     (qmckl_context context, int64_t* const type_nucl_num);
qmckl_exit_code  qmckl_get_jastrow_champ_type_nucl_vector  (qmckl_context context, int64_t* const type_nucl_num, const int64_t size_max);
qmckl_exit_code  qmckl_get_jastrow_champ_a_vector          (qmckl_context context, double * const a_vector, const int64_t size_max);
qmckl_exit_code  qmckl_get_jastrow_champ_b_vector          (qmckl_context context, double * const b_vector, const int64_t size_max);
qmckl_exit_code  qmckl_get_jastrow_champ_c_vector          (qmckl_context context, double * const c_vector, const int64_t size_max);
qmckl_exit_code  qmckl_get_jastrow_champ_rescale_factor_ee (const qmckl_context context, double* const rescale_factor_ee);
qmckl_exit_code  qmckl_get_jastrow_champ_rescale_factor_en (const qmckl_context context, double* const rescale_factor_en, const int64_t size_max);
qmckl_exit_code  qmckl_get_jastrow_champ_dim_c_vector      (qmckl_context context, int64_t* const dim_c_vector);
qmckl_exit_code  qmckl_get_jastrow_champ_spin_independent  (qmckl_context context, int32_t* const spin_independent);




/* Along with these core functions, calculation of the jastrow factor */
/* requires the following additional information to be set: */


/* When all the data for the AOs have been provided, the following */
/* function returns ~true~. */


bool      qmckl_jastrow_champ_provided           (const qmckl_context context);

/* Get */

qmckl_exit_code
qmckl_get_jastrow_champ_asymp_jasb(qmckl_context context,
                             double* const asymp_jasb,
                             const int64_t size_max);

/* Get */


qmckl_exit_code qmckl_get_jastrow_champ_ee_distance_rescaled(qmckl_context context,
                                                             double* const distance_rescaled,
                                                             int64_t const size_max);

/* Get */


qmckl_exit_code qmckl_get_jastrow_champ_ee_distance_rescaled_gl(qmckl_context context,
                                                                double* const distance_rescaled_gl,
                                                                const int64_t size_max);

/* Get */

qmckl_exit_code
qmckl_get_jastrow_champ_factor_ee(qmckl_context context,
                            double* const factor_ee,
                            const int64_t size_max);

/* Get */

qmckl_exit_code
qmckl_get_jastrow_champ_factor_ee_gl(qmckl_context context,
                                    double* const factor_ee_gl,
                                    const int64_t size_max);

/* Get */

qmckl_exit_code
qmckl_get_jastrow_champ_asymp_jasb_pderiv(qmckl_context context,
                                          double* const asymp_jasb_pderiv,
                                          const int64_t size_max);

/* Get */

qmckl_exit_code
qmckl_get_jastrow_champ_factor_ee_pderiv(qmckl_context context,
                                         double* const factor_ee_pderiv,
                                         const int64_t size_max);

/* Get */

qmckl_exit_code
qmckl_get_jastrow_champ_factor_ee_gl_pderiv(qmckl_context context,
                                            double* const factor_ee_gl_pderiv,
                                            const int64_t size_max);

/* Get */

qmckl_exit_code
qmckl_get_jastrow_champ_asymp_jasa(qmckl_context context,
                                   double* const asymp_jasa,
                                   const int64_t size_max);



/* #+CALL: generate_private_c_header(table=qmckl_asymp_jasa_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_jastrow_champ_asymp_jasa (
      const qmckl_context context,
      const int64_t aord_num,
      const int64_t type_nucl_num,
      const double* a_vector,
      const double* rescale_factor_en,
      double* const asymp_jasa );

/* Get */


qmckl_exit_code
qmckl_get_jastrow_champ_en_distance_rescaled(qmckl_context context,
                                             double* const distance_rescaled,
                                             const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_jastrow_champ_en_distance_rescaled_gl(qmckl_context context,
                                                double* const distance_rescaled_gl,
                                                const int64_t size_max);

/* Get */

qmckl_exit_code
qmckl_get_jastrow_champ_factor_en(qmckl_context context,
                            double* const factor_en,
                            const int64_t size_max);



/* #+CALL: generate_private_c_header(table=qmckl_factor_en_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_jastrow_champ_factor_en (
         const qmckl_context context,
         const int64_t walk_num,
         const int64_t elec_num,
         const int64_t nucl_num,
         const int64_t type_nucl_num,
         const int64_t* type_nucl_vector,
         const int64_t aord_num,
         const double* a_vector,
         const double* en_distance_rescaled,
         const double* asymp_jasa,
         double* const factor_en );

qmckl_exit_code qmckl_compute_jastrow_champ_factor_en_doc (
         const qmckl_context context,
         const int64_t walk_num,
         const int64_t elec_num,
         const int64_t nucl_num,
         const int64_t type_nucl_num,
         const int64_t* type_nucl_vector,
         const int64_t aord_num,
         const double* a_vector,
         const double* en_distance_rescaled,
         const double* asymp_jasa,
         double* const factor_en );

qmckl_exit_code qmckl_compute_jastrow_champ_factor_en_hpc (
         const qmckl_context context,
         const int64_t walk_num,
         const int64_t elec_num,
         const int64_t nucl_num,
         const int64_t type_nucl_num,
         const int64_t* type_nucl_vector,
         const int64_t aord_num,
         const double* a_vector,
         const double* en_distance_rescaled,
         const double* asymp_jasa,
         double* const factor_en );

/* Get */

qmckl_exit_code
qmckl_get_jastrow_champ_factor_en_gl(qmckl_context context,
                                    double* const factor_en_gl,
                                    const int64_t size_max);

/* Get */

qmckl_exit_code
qmckl_get_jastrow_champ_asymp_jasa_pderiv(qmckl_context context,
                                          double* const asymp_jasa_pderiv,
                                          const int64_t size_max);

qmckl_exit_code qmckl_compute_jastrow_champ_asymp_jasa_pderiv (
      const qmckl_context context,
      const int64_t aord_num,
      const int64_t type_nucl_num,
      const double* a_vector,
      const double* rescale_factor_en,
      double* const asymp_jasa_pderiv );

/* Get */

qmckl_exit_code
qmckl_get_jastrow_champ_factor_en_pderiv(qmckl_context context,
                                         double* const factor_en_pderiv,
                                         const int64_t size_max);

qmckl_exit_code qmckl_compute_jastrow_champ_factor_en_pderiv (
         const qmckl_context context,
         const int64_t walk_num,
         const int64_t elec_num,
         const int64_t nucl_num,
         const int64_t type_nucl_num,
         const int64_t* type_nucl_vector,
         const int64_t aord_num,
         const double* a_vector,
         const double* en_distance_rescaled,
         const double* asymp_jasa_pderiv,
         double* const factor_en_pderiv );

qmckl_exit_code qmckl_compute_jastrow_champ_factor_en_pderiv_doc (
         const qmckl_context context,
         const int64_t walk_num,
         const int64_t elec_num,
         const int64_t nucl_num,
         const int64_t type_nucl_num,
         const int64_t* type_nucl_vector,
         const int64_t aord_num,
         const double* a_vector,
         const double* en_distance_rescaled,
         const double* asymp_jasa_pderiv,
         double* const factor_en_pderiv );

/* Get */

qmckl_exit_code
qmckl_get_jastrow_champ_factor_en_gl_pderiv(qmckl_context context,
                                    double* const factor_en_gl_pderiv,
                                    const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_jastrow_champ_een_rescaled_e(qmckl_context context,
                                                double* const een_rescaled_e,
                                                const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_jastrow_champ_een_rescaled_e_gl(qmckl_context context,
                                          double* const een_rescaled_e_gl,
                                          const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_jastrow_champ_een_rescaled_n(qmckl_context context,
                                 double* const een_rescaled_n,
                                 const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_jastrow_champ_een_rescaled_n_gl(qmckl_context context,
                                         double* const een_rescaled_n_gl,
                                         const int64_t size_max);


/* #+CALL: generate_private_c_header(table=lkpm_combined_index_args,rettyp=get_value("CRetType"),fname="qmckl_compute_lkpm_combined_index_doc") */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_lkpm_combined_index_doc (
      const qmckl_context context,
      const int64_t cord_num,
      const int64_t dim_c_vector,
      int64_t* const lkpm_combined_index );



/* #+CALL: generate_private_c_header(table=lkpm_combined_index_args,rettyp=get_value("CRetType"),fname="qmckl_compute_lkpm_combined_index_hpc") */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_lkpm_combined_index_hpc (
      const qmckl_context context,
      const int64_t cord_num,
      const int64_t dim_c_vector,
      int64_t* const lkpm_combined_index );





/* #+CALL: generate_private_c_header(table=qmckl_factor_tmp_c_args,rettyp=get_value("CRetType"),fname="qmckl_compute_tmp_c") */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_tmp_c (
          const qmckl_context context,
          const int64_t cord_num,
          const int64_t elec_num,
          const int64_t nucl_num,
          const int64_t walk_num,
          const double* een_rescaled_e,
          const double* een_rescaled_n,
          double* const tmp_c );

/* Get */

qmckl_exit_code
qmckl_get_jastrow_champ_factor_een(qmckl_context context,
                             double* const factor_een,
                             const int64_t size_max);

/* Get */

qmckl_exit_code
qmckl_get_jastrow_champ_factor_een_gl(qmckl_context context,
                                     double* const factor_een_gl,
                                     const int64_t size_max);

qmckl_exit_code
qmckl_get_jastrow_champ_factor_een_grad(qmckl_context context,
                                        double* const factor_een_grad,
                                        const int64_t size_max);

/* Get */

qmckl_exit_code
qmckl_get_jastrow_champ_factor_een_pderiv(qmckl_context context,
                                          double* const factor_een_pderiv,
                                          const int64_t size_max);

/* Get */

qmckl_exit_code
qmckl_get_jastrow_champ_factor_een_gl_pderiv(qmckl_context context,
                                     double* const factor_een_gl_pderiv,
                                     const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_jastrow_champ_value(qmckl_context context,
                            double* const value,
                            const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_jastrow_champ_gl(qmckl_context context,
                            double* const gl,
                            const int64_t size_max);

qmckl_exit_code
qmckl_get_jastrow_champ_grad(qmckl_context context,
                            double* const grad,
                            const int64_t size_max);

qmckl_exit_code
qmckl_get_jastrow_champ_grad(qmckl_context context,
                            double* const grad,
                            const int64_t size_max);
