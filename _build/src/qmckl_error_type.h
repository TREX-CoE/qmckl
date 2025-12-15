/* Exit codes */

/*    All the functions in the QMCkl library return with an exit code, defined as */
/*    #+NAME: type-exit-code */

typedef int32_t qmckl_exit_code;



/* #+RESULTS: */
/* :results: */

#define  QMCKL_SUCCESS                  ((qmckl_exit_code) 0)
#define  QMCKL_INVALID_ARG_1            ((qmckl_exit_code) 1)
#define  QMCKL_INVALID_ARG_2            ((qmckl_exit_code) 2)
#define  QMCKL_INVALID_ARG_3            ((qmckl_exit_code) 3)
#define  QMCKL_INVALID_ARG_4            ((qmckl_exit_code) 4)
#define  QMCKL_INVALID_ARG_5            ((qmckl_exit_code) 5)
#define  QMCKL_INVALID_ARG_6            ((qmckl_exit_code) 6)
#define  QMCKL_INVALID_ARG_7            ((qmckl_exit_code) 7)
#define  QMCKL_INVALID_ARG_8            ((qmckl_exit_code) 8)
#define  QMCKL_INVALID_ARG_9            ((qmckl_exit_code) 9)
#define  QMCKL_INVALID_ARG_10           ((qmckl_exit_code) 10)
#define  QMCKL_INVALID_ARG_11           ((qmckl_exit_code) 11)
#define  QMCKL_INVALID_ARG_12           ((qmckl_exit_code) 12)
#define  QMCKL_INVALID_ARG_13           ((qmckl_exit_code) 13)
#define  QMCKL_INVALID_ARG_14           ((qmckl_exit_code) 14)
#define  QMCKL_INVALID_ARG_15           ((qmckl_exit_code) 15)
#define  QMCKL_INVALID_ARG_16           ((qmckl_exit_code) 16)
#define  QMCKL_INVALID_ARG_17           ((qmckl_exit_code) 17)
#define  QMCKL_INVALID_ARG_18           ((qmckl_exit_code) 18)
#define  QMCKL_INVALID_ARG_19           ((qmckl_exit_code) 19)
#define  QMCKL_INVALID_ARG_20           ((qmckl_exit_code) 20)
#define  QMCKL_FAILURE                  ((qmckl_exit_code) 101)
#define  QMCKL_ERRNO                    ((qmckl_exit_code) 102)
#define  QMCKL_INVALID_CONTEXT          ((qmckl_exit_code) 103)
#define  QMCKL_ALLOCATION_FAILED        ((qmckl_exit_code) 104)
#define  QMCKL_DEALLOCATION_FAILED      ((qmckl_exit_code) 105)
#define  QMCKL_NOT_PROVIDED             ((qmckl_exit_code) 106)
#define  QMCKL_OUT_OF_BOUNDS            ((qmckl_exit_code) 107)
#define  QMCKL_ALREADY_SET              ((qmckl_exit_code) 108)
#define  QMCKL_INVALID_EXIT_CODE        ((qmckl_exit_code) 109)
