

/* The ~qmckl_context_check~ function checks if the domain pointed to by */
/* the pointer is a valid context by verifying the magic tag value. It returns  */
/* the input ~qmckl_context~ if the context is valid (tag matches ~VALID_TAG~),  */
/* or ~QMCKL_NULL_CONTEXT~ if the context is invalid. This function should be  */
/* called at the beginning of every public API function to ensure the context */
/* parameter is valid before attempting to use it. */


qmckl_context
qmckl_context_check (const qmckl_context context) ;



/* The context keeps a /date/ that allows to check which data needs */
/* to be recomputed. The date is incremented when the context is touched. */

/* When a new element is added to the context, the functions */
/* [[Creation][=qmckl_context_create=]] [[Destroy][=qmckl_context_destroy=]] and [[Copy][=qmckl_context_copy=]] */
/* should be updated in order to make deep copies. */

/* When the electron coordinates have changed, the context is touched */
/* using the following function. */


qmckl_exit_code
qmckl_context_touch (const qmckl_context context);

/* Creation */

/*    To create a new context, ~qmckl_context_create()~ should be used. */
/*    - Upon success, it returns a pointer to a new context with the ~qmckl_context~ type */
/*    - It returns ~QMCKL_NULL_CONTEXT~ upon failure to allocate the internal data structure */
/*    - A new context always has all its members initialized with a NULL value */

/*    # Header */

qmckl_context qmckl_context_create();

/* Locking */

/*    For thread safety, the context may be locked/unlocked. The lock is */
/*    initialized with the ~PTHREAD_MUTEX_RECURSIVE~ attribute, and the */
/*    number of times the thread has locked it is saved in the */
/*    ~lock_count~ attribute. */

/*    # Header */

void qmckl_lock  (qmckl_context context);
void qmckl_unlock(qmckl_context context);

/* Copy */

/*    ~qmckl_context_copy~ makes a deep copy of a context. It returns */
/*    ~QMCKL_NULL_CONTEXT~ upon failure. */

/*    # Header */

qmckl_context qmckl_context_copy(const qmckl_context context);

/* Destroy */

/*    The context is destroyed with ~qmckl_context_destroy~, leaving the ancestors untouched. */
/*    It frees the context, and returns the previous context. */

/*    # Header */

qmckl_exit_code
qmckl_context_destroy (const qmckl_context context);
