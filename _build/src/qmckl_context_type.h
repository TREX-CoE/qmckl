/* Context handling */

/*   The context is the central data structure in QMCkl, serving as a handle */
/*   for the complete state of the library. All QMCkl functions require a context */
/*   as their first argument, and all computed data is stored within the context. */
  
/*   The context variable is a handle for the state of the library, and is stored  */
/*   in a data structure which can't be seen outside of the library. This  */
/*   encapsulation provides a clean API boundary and allows internal implementation  */
/*   details to change without affecting user code. */
  
/*   To simplify compatibility with other languages (particularly Fortran), the */
/*   pointer to the internal data structure is converted into a 64-bit */
/*   signed integer, defined in the ~qmckl_context~ type. This approach avoids */
/*   issues with language interoperability related to opaque pointer types. */
/*   A value of ~QMCKL_NULL_CONTEXT~ for the context is equivalent to a */
/*   ~NULL~ pointer and represents an invalid or uninitialized context. */

/*   #+NAME: qmckl_context */

typedef int64_t qmckl_context ;
#define QMCKL_NULL_CONTEXT (qmckl_context) 0
