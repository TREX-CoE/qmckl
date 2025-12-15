

/* #+CALL: generate_private_c_header(table=qmckl_sm_naive_args,rettyp=get_value("CRetType"),fname="qmckl_sm_naive_hpc") */

/* #+RESULTS: */

qmckl_exit_code qmckl_sm_naive_hpc (
      const qmckl_context context,
      const uint64_t LDS,
      const uint64_t Dim,
      const uint64_t N_updates,
      const double* Updates,
      const uint64_t* Updates_index,
      const double breakdown,
      double* Slater_inv,
      double* determinant );



/* #+CALL: generate_c_header(table=qmckl_sm_splitting_core_args,rettyp=get_value("CRetType"),fname="qmckl_sm_splitting_core_hpc") */

/* #+RESULTS: */

qmckl_exit_code qmckl_sm_splitting_core_hpc (
      const qmckl_context context,
      const uint64_t LDS,
      const uint64_t Dim,
      const uint64_t N_updates,
      const double* Updates,
      const uint64_t* Updates_index,
      const double breakdown,
      double* Slater_inv,
      double* later_updates,
      uint64_t* later_index,
      uint64_t* later,
      double* determinant );



/* #+CALL: generate_private_c_header(table=qmckl_sm_splitting_args,rettyp=get_value("CRetType"),fname="qmckl_sm_splitting_hpc") */

/* #+RESULTS: */

qmckl_exit_code qmckl_sm_splitting_hpc (
      const qmckl_context context,
      const uint64_t LDS,
      const uint64_t Dim,
      const uint64_t N_updates,
      const double* Updates,
      const uint64_t* Updates_index,
      const double breakdown,
      double* Slater_inv,
      double* determinant );
