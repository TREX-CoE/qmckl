/* Access functions */


qmckl_exit_code
qmckl_get_nucleus_num(const qmckl_context context,
                      int64_t* const num);

qmckl_exit_code
qmckl_get_nucleus_charge(const qmckl_context context,
                         double* const charge,
                         const int64_t size_max);

qmckl_exit_code
qmckl_get_nucleus_coord(const qmckl_context context,
                        const char transp,
                        double* const coord,
                        const int64_t size_max);



/* When all the data relative to nuclei have been set, the following */
/* function returns ~true~. */


bool qmckl_nucleus_provided (const qmckl_context context);



/* To set the data relative to the nuclei in the context, the */
/* following functions need to be called. */


qmckl_exit_code
qmckl_set_nucleus_num(qmckl_context context,
                      const int64_t num);

qmckl_exit_code
qmckl_set_nucleus_charge(qmckl_context context,
                         const double* charge,
                         const int64_t size_max);

qmckl_exit_code
qmckl_set_nucleus_coord(qmckl_context context,
                        const char transp,
                        const double* coord,
                        const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_nucleus_nn_distance(qmckl_context context,
                              double* distance,
                              const int64_t size_max);

/* Get */


qmckl_exit_code qmckl_get_nucleus_repulsion(qmckl_context context, double* const energy);
