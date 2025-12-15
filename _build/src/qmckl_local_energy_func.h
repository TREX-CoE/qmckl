/* Access functions */

/*    When all the data for the local energy have been provided, the following */
/*    function returns ~true~. */


bool      qmckl_local_energy_provided           (const qmckl_context context);

/* Get */


qmckl_exit_code qmckl_get_kinetic_energy(qmckl_context context, double* const kinetic_energy);



/* #+CALL: generate_private_c_header(table=qmckl_compute_kinetic_energy_args,rettyp=get_value("CRetType"),fname="qmckl_compute_kinetic_energy")) */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_kinetic_energy (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t det_num_alpha,
      const int64_t det_num_beta,
      const int64_t alpha_num,
      const int64_t beta_num,
      const int64_t elec_num,
      const int64_t* mo_index_alpha,
      const int64_t* mo_index_beta,
      const int64_t mo_num,
      const double* mo_vgl,
      const double* det_value_alpha,
      const double* det_value_beta,
      const double* det_inv_matrix_alpha,
      const double* det_inv_matrix_beta,
      double* const e_kin );

/* Get */


qmckl_exit_code qmckl_get_potential_energy(qmckl_context context, double* const potential_energy);



/*  #+CALL: generate_private_c_header(table=qmckl_compute_potential_energy_args,rettyp=get_value("CRetType"),fname="qmckl_compute_potential_energy")) */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_potential_energy (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t elec_num,
      const int64_t nucl_num,
      const double* ee_potential,
      const double* en_potential,
      const double repulsion,
      double* const e_pot );

/* Get */


qmckl_exit_code qmckl_get_local_energy(qmckl_context context, double* const local_energy, const int64_t size_max);



/*  #+CALL: generate_private_c_header(table=qmckl_compute_local_energy_args,rettyp=get_value("CRetType"),fname="qmckl_compute_local_energy")) */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_local_energy (
      const qmckl_context context,
      const int64_t walk_num,
      const double* e_kin,
      const double* e_pot,
      double* const e_local );

/* Get */


qmckl_exit_code qmckl_get_drift_vector(qmckl_context context, double* const drift_vector);



/*  #+CALL: generate_private_c_header(table=qmckl_compute_drift_vector_args,rettyp=get_value("CRetType"),fname="qmckl_compute_drift_vector")) */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_drift_vector (
      const qmckl_context context,
      const int64_t walk_num,
      const int64_t det_num_alpha,
      const int64_t det_num_beta,
      const int64_t alpha_num,
      const int64_t beta_num,
      const int64_t elec_num,
      const int64_t* mo_index_alpha,
      const int64_t* mo_index_beta,
      const int64_t mo_num,
      const double* mo_vgl,
      const double* det_inv_matrix_alpha,
      const double* det_inv_matrix_beta,
      double* const r_drift );
