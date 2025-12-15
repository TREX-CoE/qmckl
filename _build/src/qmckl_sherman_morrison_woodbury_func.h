qmckl_exit_code qmckl_sm_naive (
      const qmckl_context context,
      const uint64_t LDS,
      const uint64_t Dim,
      const uint64_t N_updates,
      const double* Updates,
      const uint64_t* Updates_index,
      const double breakdown,
      double* Slater_inv,
      double* determinant );



/* #+CALL: generate_c_header(table=qmckl_sm_naive_args,rettyp=get_value("CRetType"),fname="qmckl_sm_naive_doc") */

/* #+RESULTS: */

qmckl_exit_code qmckl_sm_naive_doc (
      const qmckl_context context,
      const uint64_t LDS,
      const uint64_t Dim,
      const uint64_t N_updates,
      const double* Updates,
      const uint64_t* Updates_index,
      const double breakdown,
      double* Slater_inv,
      double* determinant );

/* C headers (exposed in qmckl.h) */
/* #+CALL: generate_c_header(table=qmckl_sm_splitting_core_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */

/* #+RESULTS: */

qmckl_exit_code qmckl_sm_splitting_core (
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

qmckl_exit_code qmckl_sm_splitting_core_doc (
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

/* C headers (exposed in qmckl.h) */
/* #+CALL: generate_c_header(table=qmckl_woodbury_2x2_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */

/* #+RESULTS: */

qmckl_exit_code qmckl_woodbury_2x2 (
      const qmckl_context context,
      const uint64_t LDS,
      const uint64_t Dim,
      const double* Updates,
      const uint64_t* Updates_index,
      const double breakdown,
      double* Slater_inv,
      double* determinant );



/* #+CALL: generate_c_header(table=qmckl_woodbury_2x2_args,rettyp=get_value("CRetType"),fname="qmckl_woodbury_2x2_hpc") */

/* #+RESULTS: */

qmckl_exit_code qmckl_woodbury_2x2_hpc (
      const qmckl_context context,
      const uint64_t LDS,
      const uint64_t Dim,
      const double* Updates,
      const uint64_t* Updates_index,
      const double breakdown,
      double* Slater_inv,
      double* determinant );



/* #+CALL: generate_c_header(table=qmckl_woodbury_2x2_args,rettyp=get_value("CRetType"),fname="qmckl_woodbury_2x2_doc") */

/* #+RESULTS: */

qmckl_exit_code qmckl_woodbury_2x2_doc (
      const qmckl_context context,
      const uint64_t LDS,
      const uint64_t Dim,
      const double* Updates,
      const uint64_t* Updates_index,
      const double breakdown,
      double* Slater_inv,
      double* determinant );

/* C headers (exposed in qmckl.h) */
/* #+CALL: generate_c_header(table=qmckl_woodbury_3x3_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */

/* #+RESULTS: */

qmckl_exit_code qmckl_woodbury_3x3 (
      const qmckl_context context,
      const uint64_t LDS,
      const uint64_t Dim,
      const double* Updates,
      const uint64_t* Updates_index,
      const double breakdown,
      double* Slater_inv,
      double* determinant );



/* #+CALL: generate_c_header(table=qmckl_woodbury_3x3_args,rettyp=get_value("CRetType"),fname="qmckl_woodbury_3x3_hpc") */

/* #+RESULTS: */

qmckl_exit_code qmckl_woodbury_3x3_hpc (
      const qmckl_context context,
      const uint64_t LDS,
      const uint64_t Dim,
      const double* Updates,
      const uint64_t* Updates_index,
      const double breakdown,
      double* Slater_inv,
      double* determinant );



/* #+CALL: generate_c_header(table=qmckl_woodbury_3x3_args,rettyp=get_value("CRetType"),fname="qmckl_woodbury_3x3_doc") */

/* #+RESULTS: */

qmckl_exit_code qmckl_woodbury_3x3_doc (
      const qmckl_context context,
      const uint64_t LDS,
      const uint64_t Dim,
      const double* Updates,
      const uint64_t* Updates_index,
      const double breakdown,
      double* Slater_inv,
      double* determinant );

/* C headers (exposed in qmckl.h) */

/* #+CALL: generate_c_header(table=qmckl_sm_splitting_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */

/* #+RESULTS: */

qmckl_exit_code qmckl_sm_splitting (
      const qmckl_context context,
      const uint64_t LDS,
      const uint64_t Dim,
      const uint64_t N_updates,
      const double* Updates,
      const uint64_t* Updates_index,
      const double breakdown,
      double* Slater_inv,
      double* determinant );



/* #+CALL: generate_c_header(table=qmckl_sm_splitting_args,rettyp=get_value("CRetType"),fname="qmckl_sm_splitting_doc") */

/* #+RESULTS: */

qmckl_exit_code qmckl_sm_splitting_doc (
      const qmckl_context context,
      const uint64_t LDS,
      const uint64_t Dim,
      const uint64_t N_updates,
      const double* Updates,
      const uint64_t* Updates_index,
      const double breakdown,
      double* Slater_inv,
      double* determinant );
