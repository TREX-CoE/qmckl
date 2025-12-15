/* Precision */
/*   ~qmckl_context_set_numprec_precision~ modifies the parameter for the */
/*   numerical precision in the context. */

/*   # Header */

qmckl_exit_code qmckl_set_numprec_precision(const qmckl_context context, const int precision);



/* ~qmckl_get_numprec_precision~ returns the value of the numerical precision in the context. */


int32_t qmckl_get_numprec_precision(const qmckl_context context);

/* Range */

/*    ~qmckl_set_numprec_range~ modifies the parameter for the numerical */
/*    range in a given context. */

/*    # Header */

qmckl_exit_code qmckl_set_numprec_range(const qmckl_context context, const int range);



/* ~qmckl_get_numprec_range~ returns the value of the numerical range in the context. */


int32_t qmckl_get_numprec_range(const qmckl_context context);

/* Epsilon */

/*    ~qmckl_get_numprec_epsilon~ returns $\epsilon = 2^{1-n}$ where ~n~ is the precision. */
/*    We need to remove the sign bit from the precision. */


double qmckl_get_numprec_epsilon(const qmckl_context context);

int32_t qmckl_test_precision_64(double a, double b);
int32_t qmckl_test_precision_32(float a, float b);
