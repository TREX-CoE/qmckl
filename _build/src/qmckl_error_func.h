/* Decoding errors */

/*    To facilitate debugging and error reporting, QMCkl provides the  */
/*    ~qmckl_string_of_error~ function which converts an error code into  */
/*    a descriptive string. This allows applications to present meaningful  */
/*    error messages to users or write detailed logs for troubleshooting. */
   
/*    The function takes an error code as input and returns a constant string */
/*    describing the error condition. Both C and Fortran interfaces are provided */
/*    for maximum compatibility. */


const char*
qmckl_string_of_error (const qmckl_exit_code error);

/* Updating errors in the context */

/*    The error state in the context is updated using the ~qmckl_set_error~ */
/*    function. This function provides a centralized mechanism for recording */
/*    errors that occur during library operations. */
   
/*    When an error is set in the context, it is mandatory to specify */
/*    from which function the error is triggered, and a message */
/*    explaining the error. The exit code can't be ~QMCKL_SUCCESS~. */
   
/*    This detailed error information enables precise debugging and helps */
/*    users understand exactly what went wrong and where, making it easier */
/*    to diagnose and fix issues in code using the QMCkl library. */

/*    # Header */

qmckl_exit_code
qmckl_set_error(qmckl_context context,
                const qmckl_exit_code exit_code,
                const char* function_name,
                const char* message);

/* Get the error */

/*   Upon error, the calling program can retrieve detailed error information from the */
/*   context using ~qmckl_get_error~. This function provides access to the error  */
/*   code, the name of the function where the error occurred, and a descriptive  */
/*   message explaining the error condition. */
  
/*   The error message and function name are returned in the variables provided  */
/*   by the caller. Therefore, passing valid pointers for the function name and  */
/*   message is mandatory. The caller must ensure that the provided buffers are */
/*   large enough to hold the error information. */
  
/*   This retrieval mechanism allows applications to implement sophisticated error */
/*   handling strategies, such as retry logic, fallback mechanisms, or detailed */
/*   logging for post-mortem analysis. */

/*   # Header */

qmckl_exit_code
qmckl_get_error(qmckl_context context,
                qmckl_exit_code *exit_code,
                char* function_name,
                char* message);

/* Failing */

/*    To make a function fail, the ~qmckl_failwith~ function should be */
/*    called, such that information about the failure is stored in */
/*    the context. The desired exit code is given as an argument, as */
/*    well as the name of the function and an error message. If the */
/*    message is ~NULL~, then the default message obtained by */
/*    ~qmckl_string_of_error~ is used. The return code of the function is */
/*    the desired return code. */
/*    Upon failure, a ~QMCKL_NULL_CONTEXT~ is returned. */


qmckl_exit_code
qmckl_failwith(qmckl_context context,
               const qmckl_exit_code exit_code,
               const char* function,
               const char* message) ;

/* Last error */

/*   Returns a string describing the last error, using ~qmckl_get_error~. */

/*   # Header */

qmckl_exit_code
qmckl_last_error(qmckl_context context, char* buffer);

/* Helper functions for debugging */

/*   The following function prints to ~stderr~ an error message is the return code is */
/*   not ~QMCKL_SUCCESS~. */

/*   # Header */

qmckl_exit_code
qmckl_check(qmckl_context context, qmckl_exit_code rc);
