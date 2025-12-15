/* Number of points */


qmckl_exit_code qmckl_get_point_num (const qmckl_context context, int64_t* const num);

/* Point coordinates */


qmckl_exit_code qmckl_get_point(const qmckl_context context,
                                const char transp,
                                double* const coord,
                                const int64_t size_max);

/* Initialization functions */

/*    When the data is set in the context, if the arrays are large */
/*    enough, we overwrite the data contained in them. */

/*    To set the data relative to the points in the context, the */
/*    following function need to be called. Here, ~num~ is the number of */
/*    points to set. */


qmckl_exit_code qmckl_set_point (qmckl_context context,
                                 const char transp,
                                 const int64_t num,
                                 const double* coord,
                                 const int64_t size_max);
