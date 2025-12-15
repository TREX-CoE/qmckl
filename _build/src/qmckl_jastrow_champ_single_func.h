/* Set */

/*    We set the coordinates of the ~num~-th electron for all walkers, where ~num~ is the electron which has to be moved. */
/*    The dimension of ~coord~ is */
/*    - [walk_num][3] if ~transp~ is ~'N'~ */
/*    - [3][walk_num] if ~transp~ is ~'T'~ */

/*    Internally, the coordinates are stored in 'N' format as opposed to elec_coord. */
/*    This function has to be called before any other functions from this module. */


qmckl_exit_code qmckl_set_single_point (qmckl_context context,
                                 const char transp,
                                 const int64_t num,
                                 const double* coord,
                                 const int64_t size_max);

/* Touch */


qmckl_exit_code
qmckl_single_touch (const qmckl_context context);

/* Get */


qmckl_exit_code qmckl_get_single_electron_ee_distance(qmckl_context context,
                                                      double* const distance,
                                                      const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_single_electron_en_distance(qmckl_context context,
                                      double* distance,
                                      const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_een_rescaled_single_e(qmckl_context context,
                                 double* const distance_rescaled,
                                 const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_een_rescaled_single_n(qmckl_context context,
                                 double* const distance_rescaled,
                                 const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_jastrow_champ_single_een(qmckl_context context,
                                 double* const delta_een,
                                 const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_een_rescaled_single_n_gl(qmckl_context context,
                                         double* const distance_rescaled,
                                         const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_een_rescaled_single_e_gl(qmckl_context context,
                                         double* const distance_rescaled,
                                         const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_jastrow_champ_single_tmp(qmckl_context context,
                                 double* const tmp,
                                 const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_jastrow_champ_single_een_gl(qmckl_context context,
                                 double* const delta_een_gl,
                                 const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_jastrow_champ_single_een_g(qmckl_context context,
                                 double* const delta_een_g,
                                 const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_jastrow_champ_single_een_pderiv(qmckl_context context,
                                 double* const delta_een_pderiv,
                                 const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_ee_rescaled_single(qmckl_context context,
                                 double* const distance_rescaled,
                                 const int64_t size_max);

/* Get */

qmckl_exit_code
qmckl_get_jastrow_champ_single_ee(qmckl_context context,
                            double* const delta_ee,
                            const int64_t size_max);

/* Get */


qmckl_exit_code qmckl_get_ee_rescaled_single_gl(qmckl_context context,
                                                double* const distance_rescaled_gl,
                                                const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_jastrow_champ_single_ee_gl(qmckl_context context,
                                    double* const delta_ee_gl,
                                    const int64_t size_max);

/* Get */

qmckl_exit_code
qmckl_get_jastrow_champ_single_ee_pderiv(qmckl_context context,
                            double* const delta_ee_pderiv,
                            const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_en_rescaled_single(qmckl_context context,
                                 double* const distance_rescaled,
                                 const int64_t size_max);

/* Get */

qmckl_exit_code
qmckl_get_jastrow_champ_single_en(qmckl_context context,
                            double* const delta_en,
                            const int64_t size_max);

qmckl_exit_code qmckl_compute_jastrow_champ_single_en (
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
         const double* en_rescaled_single,
         double* const delta_en );

qmckl_exit_code qmckl_compute_jastrow_champ_single_en_doc (
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
         const double* en_rescaled_single,
         double* const delta_en );

/* Get */


qmckl_exit_code qmckl_get_en_rescaled_single_gl(qmckl_context context,
                                                double* distance_rescaled_gl,
                                                const size_t size_max);

/* Get */

qmckl_exit_code
qmckl_get_jastrow_champ_single_en_gl(qmckl_context context,
                                    double* const delta_en_gl,
                                    const int64_t size_max);

/* Get */

qmckl_exit_code
qmckl_get_jastrow_champ_single_en_pderiv(qmckl_context context,
                            double* const delta_en_pderiv,
                            const int64_t size_max);

qmckl_exit_code qmckl_compute_jastrow_champ_single_en_pderiv (
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
         const double* en_rescaled_single,
         double* const delta_en_pderiv );

qmckl_exit_code qmckl_compute_jastrow_champ_single_en_pderiv_doc (
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
         const double* en_rescaled_single,
         double* const delta_en_pderiv );

/* Code */


qmckl_exit_code
qmckl_get_jastrow_champ_single_accept(qmckl_context context);
