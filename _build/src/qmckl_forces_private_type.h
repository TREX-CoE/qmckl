#ifndef QMCKL_FORCES_HPT
#define QMCKL_FORCES_HPT
#include <stdbool.h>

/* Data structure */


typedef struct qmckl_forces_struct{
    double  * restrict forces_jastrow_en;
    uint64_t  forces_jastrow_en_date;
    double  * restrict forces_jastrow_en_g;
    uint64_t  forces_jastrow_en_g_date;
    double  * restrict forces_jastrow_en_l;
    uint64_t  forces_jastrow_en_l_date;
    double  * restrict forces_tmp_c;
    uint64_t  forces_tmp_c_date;
    double  * restrict forces_dtmp_c;
    uint64_t  forces_dtmp_c_date;
    double  * restrict forces_een_n;
    uint64_t  forces_een_n_date;
    double  * restrict forces_jastrow_een;
    uint64_t  forces_jastrow_een_date;
    double  * restrict forces_jastrow_een_g;
    uint64_t  forces_jastrow_een_g_date;
    double  * restrict forces_jastrow_een_l;
    uint64_t  forces_jastrow_een_l_date;
    double  * restrict forces_ao_value;
    uint64_t  forces_ao_value_date;
    uint64_t  forces_ao_value_maxsize;
    double  * restrict forces_mo_value;
    uint64_t  forces_mo_value_date;
    uint64_t  forces_mo_value_maxsize;
    double  * restrict forces_mo_g;
    uint64_t  forces_mo_g_date;
    double  * restrict forces_mo_l;
    uint64_t  forces_mo_l_date;
    double  * forces_jastrow_single_en;
    uint64_t  forces_jastrow_single_en_date;
    double  * forces_jastrow_single_een;
    uint64_t  forces_jastrow_single_een_date;
} qmckl_forces_struct;

#endif
