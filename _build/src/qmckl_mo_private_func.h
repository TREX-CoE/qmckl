#ifndef QMCKL_MO_HPF
#define QMCKL_MO_HPF

#include "qmckl_mo_private_type.h"
#include "qmckl_blas_private_type.h"



/* The ~uninitialized~ integer contains one bit set to one for each */
/* initialization function which has not been called. It becomes equal */
/* to zero after all initialization functions have been called. The */
/* struct is then initialized and ~provided == true~. */
/* Some values are initialized by default, and are not concerned by */
/* this mechanism. */


qmckl_exit_code qmckl_init_mo_basis(qmckl_context context);


/* When the basis set is completely entered, other data structures are */
/* computed to accelerate the calculations. */


qmckl_exit_code qmckl_finalize_mo_basis(qmckl_context context);

/* Deep copy */


qmckl_exit_code qmckl_copy_mo_basis(qmckl_context context, const qmckl_mo_basis_struct* src, qmckl_mo_basis_struct* dest);

/* Provide */

/* #+CALL: write_provider_header( group="mo_basis", data="mo_value" ) */

/* #+RESULTS: */

qmckl_exit_code qmckl_provide_mo_basis_mo_value(qmckl_context context);



/*  #+CALL: generate_private_c_header(table=qmckl_mo_basis_mo_value_args,rettyp=get_value("CRetType"),fname="qmckl_compute_mo_basis_mo_value")) */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_mo_basis_mo_value (
      const qmckl_context context,
      const int64_t ao_num,
      const int64_t mo_num,
      const int64_t point_num,
      const double* coefficient_t,
      const double* ao_value,
      double* const mo_value );



/* #+CALL: generate_private_c_header(table=qmckl_mo_basis_mo_value_args,rettyp=get_value("CRetType"),fname="qmckl_compute_mo_basis_mo_value_doc")) */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_mo_basis_mo_value_doc (
      const qmckl_context context,
      const int64_t ao_num,
      const int64_t mo_num,
      const int64_t point_num,
      const double* coefficient_t,
      const double* ao_value,
      double* const mo_value );

/* HPC version                                                    :noexport: */



#ifdef HAVE_HPC
qmckl_exit_code
qmckl_compute_mo_basis_mo_value_hpc (const qmckl_context context,
                                     const int64_t ao_num,
                                     const int64_t mo_num,
                                     const int64_t point_num,
                                     const double* coefficient_t,
                                     const double* ao_value,
                                     double* const mo_value );

qmckl_exit_code
qmckl_compute_mo_basis_mo_value_hpc_sp (const qmckl_context context,
                                        const int64_t ao_num,
                                        const int64_t mo_num,
                                        const int64_t point_num,
                                        const double* coefficient_t,
                                        const double* ao_value,
                                        double* const mo_value );
#endif

/* Provide */

/* #+CALL: write_provider_header( group="mo_basis", data="mo_vgl" ) */

/* #+RESULTS: */

qmckl_exit_code qmckl_provide_mo_basis_mo_vgl(qmckl_context context);



/*  #+CALL: generate_private_c_header(table=qmckl_mo_basis_mo_vgl_args,rettyp=get_value("CRetType"),fname="qmckl_compute_mo_basis_mo_vgl")) */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_mo_basis_mo_vgl (
      const qmckl_context context,
      const int64_t ao_num,
      const int64_t mo_num,
      const int64_t point_num,
      const double* coefficient_t,
      const double* ao_vgl,
      double* const mo_vgl );



/* #+CALL: generate_private_c_header(table=qmckl_mo_basis_mo_vgl_args,rettyp=get_value("CRetType"),fname="qmckl_compute_mo_basis_mo_vgl_doc")) */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_mo_basis_mo_vgl_doc (
      const qmckl_context context,
      const int64_t ao_num,
      const int64_t mo_num,
      const int64_t point_num,
      const double* coefficient_t,
      const double* ao_vgl,
      double* const mo_vgl );

/* Double precision */


#ifdef HAVE_HPC
qmckl_exit_code
qmckl_compute_mo_basis_mo_vgl_hpc (const qmckl_context context,
                                   const int64_t ao_num,
                                   const int64_t mo_num,
                                   const int64_t point_num,
                                   const double* coefficient_t,
                                   const double* ao_vgl,
                                   double* const mo_vgl );
#endif

/* Single precision */


#ifdef HAVE_HPC
qmckl_exit_code
qmckl_compute_mo_basis_mo_vgl_hpc_sp (const qmckl_context context,
                                      const int64_t ao_num,
                                      const int64_t mo_num,
                                      const int64_t point_num,
                                      const double* coefficient_t,
                                      const double* ao_vgl,
                                      double* const mo_vgl );
#endif



/* #+CALL: generate_private_c_header(table=qmckl_mo_basis_mo_value_cusp_args,rettyp=get_value("CRetType"),fname="qmckl_compute_mo_basis_mo_value_cusp")) */

/* #+RESULTS: */

qmckl_exit_code qmckl_compute_mo_basis_mo_value_cusp (
      const qmckl_context context,
      const int64_t nucl_num,
      const int64_t ao_num,
      const int64_t mo_num,
      const int64_t point_num,
      const int64_t* ao_nucl,
      const int32_t* ao_ang_mom,
      const double* en_distance,
      const double* r_cusp,
      const qmckl_tensor cusp_param,
      const double* coefficient_t,
      const double* ao_value,
      double* const mo_value );

/* HPC version                                                    :noexport: */


#ifdef HAVE_HPC
qmckl_exit_code
qmckl_compute_mo_basis_mo_value_cusp_hpc (const qmckl_context context,
                                          const int64_t nucl_num,
                                          const int64_t ao_num,
                                          const int64_t mo_num,
                                          const int64_t point_num,
                                          const int64_t* ao_nucl,
                                          const int32_t* ao_ang_mom,
                                          const double* en_distance,
                                          const double* r_cusp,
                                          const qmckl_tensor cusp_param,
                                          const double* coefficient_t,
                                          const double* ao_value,
                                          double* const mo_value );

qmckl_exit_code
qmckl_compute_mo_basis_mo_value_cusp_hpc_sp (const qmckl_context context,
                                             const int64_t nucl_num,
                                             const int64_t ao_num,
                                             const int64_t mo_num,
                                             const int64_t point_num,
                                             const int64_t* ao_nucl,
                                             const int32_t* ao_ang_mom,
                                             const double* en_distance,
                                             const double* r_cusp,
                                             const qmckl_tensor cusp_param,
                                             const double* coefficient_t,
                                             const double* ao_value,
                                             double* const mo_value );
#endif



/* #   #+CALL: generate_private_c_header(table=qmckl_mo_basis_mo_vgl_cusp_args,rettyp=get_value("CRetType"),fname="qmckl_compute_mo_basis_mo_vgl_cusp")) */

/*    #+RESULTS: */

qmckl_exit_code qmckl_compute_mo_basis_mo_vgl_cusp (
      const qmckl_context context,
      const int64_t nucl_num,
      const int64_t ao_num,
      const int64_t mo_num,
      const int64_t point_num,
      const int64_t* ao_nucl,
      const int32_t* ao_ang_mom,
      const double* en_distance,
      const qmckl_matrix nucl_coord,
      const qmckl_matrix point_coord,
      const double* r_cusp,
      const qmckl_tensor cusp_param,
      const double* coefficient_t,
      const double* ao_vgl,
      double* const mo_vgl );

/* HPC version                                                    :noexport: */



#ifdef HAVE_HPC
qmckl_exit_code
qmckl_compute_mo_basis_mo_vgl_cusp_hpc (const qmckl_context context,
                                        const int64_t nucl_num,
                                        const int64_t ao_num,
                                        const int64_t mo_num,
                                        const int64_t point_num,
                                        const int64_t* ao_nucl,
                                        const int32_t* ao_ang_mom,
                                        const double* en_distance,
                                        const qmckl_matrix nucl_coord,
                                        const qmckl_matrix point_coord,
                                        const double* r_cusp,
                                        const qmckl_tensor cusp_param,
                                        const double* coefficient_t,
                                        const double* ao_vgl,
                                        double* const mo_vgl );
#endif

#endif
