/* Initialization functions */

/*    To set the basis set, all the following functions need to be */
/*    called. */


qmckl_exit_code  qmckl_set_mo_basis_mo_num           (qmckl_context context, const int64_t   mo_num);
qmckl_exit_code  qmckl_set_mo_basis_coefficient      (qmckl_context context, const double  * coefficient, const int64_t size_max);
qmckl_exit_code  qmckl_set_mo_basis_r_cusp           (qmckl_context context, const double  * r_cusp, const int64_t size_max);

/* Access functions */


qmckl_exit_code
qmckl_get_mo_basis_mo_num (const qmckl_context context,
                           int64_t* mo_num);

qmckl_exit_code
qmckl_get_mo_basis_coefficient (const qmckl_context context,
                                double* const coefficient,
                                const int64_t size_max);



/* When all the data for the AOs have been provided, the following */
/* function returns ~true~. */


bool qmckl_mo_basis_provided (const qmckl_context context);

/* Update */

/*    It may be desirable to remove certain molecular orbitals (MOs) that */
/*    do not significantly contribute to the wave function.  In */
/*    particular, in a single determinant calculation, the virtual MOs */
/*    can be removed as they do not participate in the ground state */
/*    configuration. */

/*    To select a subset of MOs that will be kept, an array of integers of */
/*    size ~mo_num~ can be created. If the integer corresponding to an MO is */
/*    zero, that MO is dropped and will not be included in the */
/*    calculation. If the integer is non-zero, the MO will be kept. */



qmckl_exit_code
qmckl_mo_basis_select_mo (const qmckl_context context,
                          const int32_t* keep,
                          const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_mo_basis_mo_value(qmckl_context context,
                            double* const mo_value,
                            const int64_t size_max);



/* Uses the given array to compute the values. */


qmckl_exit_code
qmckl_get_mo_basis_mo_value_inplace (qmckl_context context,
                                     double* const mo_value,
                                     const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_mo_basis_mo_vgl(qmckl_context context,
                          double* const mo_vgl,
                          const int64_t size_max);



/* Uses the given array to compute the VGL. */


qmckl_exit_code
qmckl_get_mo_basis_mo_vgl_inplace (qmckl_context context,
                                   double* const mo_vgl,
                                   const int64_t size_max);



/* #+CALL: generate_private_c_header(table=qmckl_mo_basis_mo_value_cusp_args,rettyp=get_value("CRetType"),fname="qmckl_compute_mo_basis_mo_value_cusp_doc")) */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_mo_basis_mo_value_cusp_doc (
      const qmckl_context context,
      const int64_t nucl_num,
      const int64_t ao_num,
      const int64_t mo_num,
      const int64_t point_num,
      const int64_t* ao_nucl,
      const int32_t* ao_ang_mom,
      const double* en_distance,
      const double* r_cusp,
      const double* cusp_param,
      const double* coefficient_t,
      const double* ao_value,
      double* const mo_value );



/* #+CALL: generate_private_c_header(table=qmckl_mo_basis_mo_vgl_cusp_args,rettyp=get_value("CRetType"),fname="qmckl_compute_mo_basis_mo_vgl_cusp_doc")) */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_mo_basis_mo_vgl_cusp_doc (
      const qmckl_context context,
      const int64_t nucl_num,
      const int64_t ao_num,
      const int64_t mo_num,
      const int64_t point_num,
      const int64_t* ao_nucl,
      const int32_t* ao_ang_mom,
      const double* en_distance,
      const double* nucl_coord,
      const double* point_coord,
      const double* r_cusp,
      const double* cusp_param,
      const double* coefficient_t,
      const double* ao_vgl,
      double* const mo_vgl );

/* Rescaling of MO coefficients */

/*    When evaluating Slater determinants, the value of the determinants */
/*    may get out of the range of double precision. A simple fix is to */
/*    rescale the MO coefficients to put back the determinants in the */
/*    correct range. */


qmckl_exit_code
qmckl_mo_basis_rescale(qmckl_context context,
                          const double scaling_factor);
