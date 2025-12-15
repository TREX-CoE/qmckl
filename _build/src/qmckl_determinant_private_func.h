

/* The ~uninitialized~ integer contains one bit set to one for each */
/* initialization function which has not been called. It becomes equal */
/* to zero after all initialization functions have been called. The */
/* struct is then initialized and ~provided == true~. */
/* Some values are initialized by default, and are not concerned by */
/* this mechanism. */


qmckl_exit_code qmckl_init_determinant(qmckl_context context);

/* Access functions */


char      qmckl_get_determinant_type             (const qmckl_context context);
int64_t   qmckl_get_determinant_det_num_alpha    (const qmckl_context context);
int64_t   qmckl_get_determinant_det_num_beta     (const qmckl_context context);
int64_t*   qmckl_get_determinant_mo_index_alpha  (const qmckl_context context);
int64_t*   qmckl_get_determinant_mo_index_beta   (const qmckl_context context);



/* When the basis set is completely entered, other data structures are */
/* computed to accelerate the calculations. */


qmckl_exit_code qmckl_finalize_determinant(qmckl_context context);

/* Provide */


qmckl_exit_code qmckl_provide_det_vgl_alpha(qmckl_context context);
qmckl_exit_code qmckl_provide_det_vgl_beta(qmckl_context context);

/* Provide */


qmckl_exit_code qmckl_provide_det_inv_matrix_alpha(qmckl_context context);
qmckl_exit_code qmckl_provide_det_inv_matrix_beta(qmckl_context context);
