

/* When all the data for the slater determinants have been provided, the following */
/* function returns ~true~. */


bool      qmckl_determinant_provided           (const qmckl_context context);

/* Initialization functions */

/*    To set the basis set, all the following functions need to be */
/*    called. */


qmckl_exit_code  qmckl_set_determinant_type             (const qmckl_context context, const char t);
qmckl_exit_code  qmckl_set_determinant_det_num_alpha    (const qmckl_context context, const int64_t det_num_alpha);
qmckl_exit_code  qmckl_set_determinant_det_num_beta     (const qmckl_context context, const int64_t det_num_beta);
qmckl_exit_code  qmckl_set_determinant_mo_index_alpha   (const qmckl_context context, const int64_t* mo_index_alpha, const int64_t size_max);
qmckl_exit_code  qmckl_set_determinant_mo_index_beta    (const qmckl_context context, const int64_t* mo_index_beta, const int64_t size_max);

/* Get */


qmckl_exit_code qmckl_get_det_vgl_alpha(qmckl_context context, double* const det_vgl_alpha);
qmckl_exit_code qmckl_get_det_vgl_beta(qmckl_context context, double* const det_vgl_beta);



/*  #+CALL: generate_private_c_header(table=qmckl_compute_det_vgl_alpha_args,rettyp=get_value("CRetType"),fname="qmckl_compute_det_vgl_alpha")) */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_det_vgl_alpha (
      const qmckl_context context,
      const int64_t det_num_alpha,
      const int64_t walk_num,
      const int64_t alpha_num,
      const int64_t beta_num,
      const int64_t elec_num,
      const int64_t* mo_index_alpha,
      const int64_t mo_num,
      const double* mo_vgl,
      double* const det_vgl_alpha );



/*  #+CALL: generate_private_c_header(table=qmckl_compute_det_vgl_beta_args,rettyp=get_value("CRetType"),fname="qmckl_compute_det_vgl_beta")) */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_det_vgl_beta (
      const qmckl_context context,
      const int64_t det_num_beta,
      const int64_t walk_num,
      const int64_t alpha_num,
      const int64_t beta_num,
      const int64_t elec_num,
      const int64_t* mo_index_beta,
      const int64_t mo_num,
      const double* mo_vgl,
      double* const det_vgl_beta );

/* Get */


qmckl_exit_code qmckl_get_det_inv_matrix_alpha(qmckl_context context, double* const det_inv_matrix_alpha);
qmckl_exit_code qmckl_get_det_inv_matrix_beta(qmckl_context context, double* const det_inv_matrix_beta);
qmckl_exit_code qmckl_get_det_adj_matrix_alpha(qmckl_context context, double* const det_adj_matrix_alpha);
qmckl_exit_code qmckl_get_det_adj_matrix_beta(qmckl_context context, double* const det_adj_matrix_beta);
qmckl_exit_code qmckl_get_det_alpha(qmckl_context context, double* const det_adj_matrix_alpha);
qmckl_exit_code qmckl_get_det_beta(qmckl_context context, double* const det_adj_matrix_beta);



/* #+CALL: generate_private_c_header(table=qmckl_det_inv_matrix_alpha_args,rettyp=get_value("CRetType"),fname="qmckl_compute_det_inv_matrix_alpha")) */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_det_inv_matrix_alpha (
      const qmckl_context context,
      const int64_t det_num_alpha,
      const int64_t walk_num,
      const int64_t alpha_num,
      const double* det_vgl_alpha,
      double* const det_value_alpha,
      double* const det_adj_matrix_alpha,
      double* const det_inv_matrix_alpha );



/*  #+CALL: generate_private_c_header(table=qmckl_det_inv_matrix_beta_args,rettyp=get_value("CRetType"),fname="qmckl_compute_det_inv_matrix_beta")) */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_det_inv_matrix_beta (
      const qmckl_context context,
      const int64_t det_num_beta,
      const int64_t walk_num,
      const int64_t beta_num,
      const double* det_vgl_beta,
      double* const det_value_beta,
      double* const det_adj_matrix_beta,
      double* const det_inv_matrix_beta );
