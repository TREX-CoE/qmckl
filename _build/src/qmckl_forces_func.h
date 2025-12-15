/* Get */



qmckl_exit_code
qmckl_get_forces_jastrow_en(qmckl_context context,
                                    double* const forces_jastrow_en,
                                    const int64_t size_max);

/* Get */



qmckl_exit_code
qmckl_get_forces_jastrow_en_g(qmckl_context context,
                                    double* const forces_jastrow_en_g,
                                    const int64_t size_max);

/* Get */



qmckl_exit_code
qmckl_get_forces_jastrow_en_l(qmckl_context context,
                                    double* const forces_jastrow_en_l,
                                    const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_forces_jastrow_single_en(qmckl_context context,
                          double* const forces_jastrow_single_en,
                          const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_forces_tmp_c(qmckl_context context,
                                    double* const forces_tmp_c,
                                    const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_forces_dtmp_c(qmckl_context context,
                                    double* const forces_dtmp_c,
                                    const int64_t size_max);

/* Get */



qmckl_exit_code
qmckl_get_forces_jastrow_een(qmckl_context context,
                                    double* const forces_jastrow_een,
                                    const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_forces_een_rescaled_n_gl(qmckl_context context,
                                         double* const forces_een_n,
                                         const int64_t size_max);

/* Get */

qmckl_exit_code
qmckl_get_forces_jastrow_een_g(qmckl_context context,
                                     double* const forces_jastrow_een_g,
                                     const int64_t size_max);

/* Get */

qmckl_exit_code
qmckl_get_forces_jastrow_een_l(qmckl_context context,
                                     double* const forces_jastrow_een_l,
                                     const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_forces_jastrow_single_een(qmckl_context context,
                          double* const forces_jastrow_single_een,
                          const int64_t size_max);

/* Get */

qmckl_exit_code
qmckl_get_forces_ao_value(qmckl_context context,
                                     double* const forces_ao_value,
                                     const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_forces_mo_value(qmckl_context context,
                          double* const forces_mo_value,
                          const int64_t size_max);

qmckl_exit_code
qmckl_get_forces_mo_value_inplace (qmckl_context context,
                                     double* const forces_mo_value,
                                     const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_forces_mo_g(qmckl_context context,
                          double* const forces_mo_g,
                          const int64_t size_max);

qmckl_exit_code
qmckl_get_forces_mo_g_inplace(qmckl_context context,
                              double* const forces_mo_g,
                              const int64_t size_max);

qmckl_exit_code qmckl_compute_forces_mo_g (
         const qmckl_context context,
         const int64_t ao_num,
         const int64_t mo_num,
         const int64_t point_num,
         const int64_t nucl_num,
         const int64_t shell_num,
         const int64_t* nucleus_index,
         const int64_t* nucleus_shell_num,
         const int32_t* shell_ang_mom,
         const double* coefficient_t,
         const double* ao_hessian,
         double* const forces_mo_g );

qmckl_exit_code qmckl_compute_forces_mo_g_doc (
         const qmckl_context context,
         const int64_t ao_num,
         const int64_t mo_num,
         const int64_t point_num,
         const int64_t nucl_num,
         const int64_t shell_num,
         const int64_t* nucleus_index,
         const int64_t* nucleus_shell_num,
         const int32_t* shell_ang_mom,
         const double* coefficient_t,
         const double* ao_hessian,
         double* const forces_mo_g );

qmckl_exit_code qmckl_compute_forces_mo_g_hpc (
         const qmckl_context context,
         const int64_t ao_num,
         const int64_t mo_num,
         const int64_t point_num,
         const int64_t nucl_num,
         const int64_t shell_num,
         const int64_t* nucleus_index,
         const int64_t* nucleus_shell_num,
         const int32_t* shell_ang_mom,
         const double* coefficient_t,
         const double* ao_hessian,
         double* const forces_mo_g );

/* Get */


qmckl_exit_code
qmckl_get_forces_mo_l(qmckl_context context,
                          double* const forces_mo_l,
                          const int64_t size_max);



/* #+RESULTS: */

qmckl_exit_code qmckl_compute_forces_mo_l_doc (
      const qmckl_context context,
      const int64_t ao_num,
      const int64_t mo_num,
      const int64_t point_num,
      const int64_t nucl_num,
      const int64_t shell_num,
      const int64_t* nucleus_index,
      const int64_t* nucleus_shell_num,
      const int32_t* shell_ang_mom,
      const double* coefficient_t,
      const double* ao_hessian,
      double* const forces_mo_l );
