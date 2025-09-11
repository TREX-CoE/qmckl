#include "qmckl.h"
#include "assert.h"
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <math.h>
#include <string.h>

#include "chbrclf.h"
#include "qmckl_ao_private_func.h"

int main() {
    qmckl_context context;
    context = qmckl_context_create();

const int64_t   nucl_num      = chbrclf_nucl_num;
const double*   nucl_charge   = chbrclf_charge;
const double*   nucl_coord    = &(chbrclf_nucl_coord[0][0]);

qmckl_exit_code rc;
rc = qmckl_set_nucleus_num (context, nucl_num);
assert(rc == QMCKL_SUCCESS);

rc = qmckl_set_nucleus_coord (context, 'T', &(nucl_coord[0]), 3*nucl_num);
assert(rc == QMCKL_SUCCESS);

rc = qmckl_set_nucleus_charge(context, nucl_charge, nucl_num);
assert(rc == QMCKL_SUCCESS);

assert(qmckl_nucleus_provided(context));


const int64_t    shell_num         =  chbrclf_shell_num;
const int64_t    prim_num          =  chbrclf_prim_num;
const int64_t    ao_num            =  chbrclf_ao_num;
const int64_t *  nucleus_index     =  &(chbrclf_basis_nucleus_index[0]);
const int64_t *  nucleus_shell_num =  &(chbrclf_basis_nucleus_shell_num[0]);
const int32_t *  shell_ang_mom     =  &(chbrclf_basis_shell_ang_mom[0]);
const int64_t *  shell_prim_num    =  &(chbrclf_basis_shell_prim_num[0]);
const int64_t *  shell_prim_index  =  &(chbrclf_basis_shell_prim_index[0]);
const double  *  shell_factor      =  &(chbrclf_basis_shell_factor[0]);
const double  *  exponent          =  &(chbrclf_basis_exponent[0]);
const double  *  coefficient       =  &(chbrclf_basis_coefficient[0]);
const double  *  prim_factor       =  &(chbrclf_basis_prim_factor[0]);
const double  *  ao_factor         =  &(chbrclf_basis_ao_factor[0]);

const char typ = 'G';

assert(!qmckl_ao_basis_provided(context));

rc = qmckl_set_ao_basis_type (context, typ);
assert(rc == QMCKL_SUCCESS);
assert(!qmckl_ao_basis_provided(context));

rc = qmckl_set_ao_basis_shell_num (context, shell_num);
assert(rc == QMCKL_SUCCESS);
assert(!qmckl_ao_basis_provided(context));

rc = qmckl_set_ao_basis_prim_num (context, prim_num);
assert(rc == QMCKL_SUCCESS);
assert(!qmckl_ao_basis_provided(context));

rc = qmckl_set_ao_basis_nucleus_index (context, nucleus_index, nucl_num);
assert(rc == QMCKL_SUCCESS);
assert(!qmckl_ao_basis_provided(context));

rc = qmckl_set_ao_basis_nucleus_index (context, nucleus_index, nucl_num);
assert(rc == QMCKL_ALREADY_SET);

rc = qmckl_set_ao_basis_nucleus_shell_num (context, nucleus_shell_num, nucl_num);
assert(rc == QMCKL_SUCCESS);
assert(!qmckl_ao_basis_provided(context));

rc = qmckl_set_ao_basis_shell_ang_mom (context, shell_ang_mom, shell_num);
assert(rc == QMCKL_SUCCESS);
assert(!qmckl_ao_basis_provided(context));

rc = qmckl_set_ao_basis_shell_factor  (context, shell_factor, shell_num);
assert(rc == QMCKL_SUCCESS);
assert(!qmckl_ao_basis_provided(context));

rc = qmckl_set_ao_basis_shell_prim_num (context, shell_prim_num, shell_num);
assert(rc == QMCKL_SUCCESS);
assert(!qmckl_ao_basis_provided(context));

rc = qmckl_set_ao_basis_shell_prim_index (context, shell_prim_index, shell_num);
assert(rc == QMCKL_SUCCESS);
assert(!qmckl_ao_basis_provided(context));

rc = qmckl_set_ao_basis_exponent      (context, exponent, prim_num);
assert(rc == QMCKL_SUCCESS);
assert(!qmckl_ao_basis_provided(context));

rc = qmckl_set_ao_basis_coefficient   (context, coefficient, prim_num);
assert(rc == QMCKL_SUCCESS);
assert(!qmckl_ao_basis_provided(context));

rc = qmckl_set_ao_basis_prim_factor (context, prim_factor, prim_num);
assert(rc == QMCKL_SUCCESS);

rc = qmckl_set_ao_basis_ao_num(context, ao_num);
assert(rc == QMCKL_SUCCESS);

rc = qmckl_set_ao_basis_ao_factor (context, ao_factor, ao_num);
assert(rc == QMCKL_SUCCESS);

assert(qmckl_ao_basis_provided(context));

int64_t    shell_num_test        ;
int64_t    prim_num_test         ;
int64_t    ao_num_test           ;
int64_t *  nucleus_index_test    ;
int64_t *  nucleus_shell_num_test;
int32_t *  shell_ang_mom_test    ;
int64_t *  shell_prim_num_test   ;
int64_t *  shell_prim_index_test ;
double  *  shell_factor_test     ;
double  *  exponent_test         ;
double  *  coefficient_test      ;
double  *  prim_factor_test      ;
double  *  ao_factor_test        ;
char    typ_test ;


rc = qmckl_get_ao_basis_type (context, &typ_test);
assert (rc == QMCKL_SUCCESS);
assert(typ == typ_test);

rc = qmckl_get_ao_basis_shell_num (context, &shell_num_test);
assert (rc == QMCKL_SUCCESS);
assert(shell_num == shell_num_test);

rc = qmckl_get_ao_basis_prim_num (context, &prim_num_test);
assert (rc == QMCKL_SUCCESS);
assert(prim_num == prim_num_test);

nucleus_index_test = (int64_t*) malloc (nucl_num * sizeof(int64_t));
rc = qmckl_get_ao_basis_nucleus_index (context, nucleus_index_test, nucl_num);
assert (rc == QMCKL_SUCCESS);
for (int64_t i=0 ; i < nucl_num ; ++i) {
  assert(nucleus_index_test[i] == nucleus_index[i]);
 }
free(nucleus_index_test);

nucleus_shell_num_test = (int64_t*) malloc ( nucl_num * sizeof(int64_t));
rc = qmckl_get_ao_basis_nucleus_shell_num (context, nucleus_shell_num_test, nucl_num);
assert (rc == QMCKL_SUCCESS);

for (int64_t i=0 ; i < nucl_num ; ++i) {
  assert(nucleus_shell_num_test[i] == nucleus_shell_num[i]);
 }
free(nucleus_shell_num_test);

shell_ang_mom_test = (int32_t*) malloc ( shell_num * sizeof(int32_t));
rc = qmckl_get_ao_basis_shell_ang_mom (context, shell_ang_mom_test, shell_num);
assert (rc == QMCKL_SUCCESS);

for (int64_t i=0 ; i < shell_num ; ++i) {
  assert(shell_ang_mom_test[i] == shell_ang_mom[i]);
 }
free(shell_ang_mom_test);

shell_factor_test = (double*) malloc ( shell_num * sizeof(double));
rc = qmckl_get_ao_basis_shell_factor  (context, shell_factor_test, shell_num);
assert (rc == QMCKL_SUCCESS);

for (int64_t i=0 ; i < shell_num ; ++i) {
  assert(shell_factor_test[i] == shell_factor[i]);
}
free(shell_factor_test);

shell_prim_num_test = (int64_t*) malloc ( shell_num * sizeof(int64_t));
rc = qmckl_get_ao_basis_shell_prim_num (context, shell_prim_num_test, shell_num);
assert (rc == QMCKL_SUCCESS);

for (int64_t i=0 ; i < shell_num ; ++i) {
  assert(shell_prim_num_test[i] == shell_prim_num[i]);
}
free(shell_prim_num_test);

shell_prim_index_test = (int64_t*) malloc ( shell_num * sizeof(int64_t));
rc = qmckl_get_ao_basis_shell_prim_index (context, shell_prim_index_test, shell_num);
assert (rc == QMCKL_SUCCESS);

for (int64_t i=0 ; i < shell_num ; ++i) {
  assert(shell_prim_index_test[i] == shell_prim_index[i]);
}
free(shell_prim_index_test);

exponent_test = (double*) malloc ( prim_num * sizeof(double));
rc = qmckl_get_ao_basis_exponent(context, exponent_test, prim_num);
assert (rc == QMCKL_SUCCESS);

for (int64_t i=0 ; i < prim_num ; ++i) {
  assert(exponent_test[i] == exponent[i]);
}
free(exponent_test);

coefficient_test = (double*) malloc ( prim_num * sizeof(double));
rc = qmckl_get_ao_basis_coefficient(context, coefficient_test, prim_num);
assert (rc == QMCKL_SUCCESS);

for (int64_t i=0 ; i < prim_num ; ++i) {
  assert(coefficient_test[i] == coefficient[i]);
}
free(coefficient_test);

prim_factor_test = (double*) malloc ( prim_num * sizeof(double));
rc = qmckl_get_ao_basis_prim_factor (context, prim_factor_test, prim_num);
assert (rc == QMCKL_SUCCESS);

for (int64_t i=0 ; i < prim_num ; ++i) {
  assert(prim_factor_test[i] == prim_factor[i]);
}
free(prim_factor_test);

rc = qmckl_get_ao_basis_ao_num(context, &ao_num_test);
assert(ao_num == ao_num_test);

ao_factor_test = (double*) malloc ( ao_num * sizeof(double));
rc = qmckl_get_ao_basis_ao_factor (context, ao_factor_test, ao_num);
assert (rc == QMCKL_SUCCESS);

for (int64_t i=0 ; i < ao_num ; ++i) {
  assert(ao_factor_test[i] == ao_factor[i]);
}
free(ao_factor_test);

int test_qmckl_ao_gaussian_vgl(qmckl_context context);
assert(0 == test_qmckl_ao_gaussian_vgl(context));

int test_qmckl_ao_slater_vgl(qmckl_context context);
assert(0 == test_qmckl_ao_slater_vgl(context));

int test_qmckl_ao_slater_shell_vgl(qmckl_context context);
assert(0 == test_qmckl_ao_slater_shell_vgl(context));

{
#define walk_num 1 // chbrclf_walk_num
#define elec_num chbrclf_elec_num
#define prim_num chbrclf_prim_num

  int64_t elec_up_num   = chbrclf_elec_up_num;
  int64_t elec_dn_num   = chbrclf_elec_dn_num;
  double* elec_coord    = &(chbrclf_elec_coord[0][0][0]);

  rc = qmckl_set_electron_num (context, elec_up_num, elec_dn_num);
  assert (rc == QMCKL_SUCCESS);

  assert(qmckl_electron_provided(context));

  int64_t point_num = elec_num;

  rc = qmckl_set_point(context, 'N', point_num, elec_coord, point_num*3);
  assert(rc == QMCKL_SUCCESS);



  double prim_vgl[point_num][5][prim_num];

  rc = qmckl_get_ao_basis_primitive_vgl(context, &(prim_vgl[0][0][0]),
          (int64_t) 5*point_num*prim_num );
  assert (rc == QMCKL_SUCCESS);

  printf("prim_vgl[26][0][7] = %e\n",prim_vgl[26][0][7]);
  assert( fabs(prim_vgl[26][0][7] - ( 1.0501570432064878E-003)) < 1.e-14 );
  printf("prim_vgl[26][1][7] = %e\n",prim_vgl[26][1][7]);
  assert( fabs(prim_vgl[26][1][7] - (-7.5014974095310560E-004)) < 1.e-14 );
  printf("prim_vgl[26][2][7] = %e\n",prim_vgl[26][2][7]);
  assert( fabs(prim_vgl[26][2][7] - (-3.8250692897610380E-003)) < 1.e-14 );
  printf("prim_vgl[26][3][7] = %e\n",prim_vgl[26][3][7]);
  assert( fabs(prim_vgl[26][3][7] - ( 3.4950559194080275E-003)) < 1.e-14 );
  printf("prim_vgl[26][4][7] = %e\n",prim_vgl[26][4][7]);
  assert( fabs(prim_vgl[26][4][7] - ( 2.0392163767356572E-002)) < 1.e-14 );

}

{
#define shell_num chbrclf_shell_num

  double* elec_coord    = &(chbrclf_elec_coord[0][0][0]);

  assert(qmckl_electron_provided(context));

  int64_t point_num = elec_num;
  rc = qmckl_set_point(context, 'N', point_num, elec_coord, point_num*3);
  assert(rc == QMCKL_SUCCESS);


  double shell_vgl[point_num][5][shell_num];

  rc = qmckl_get_ao_basis_shell_vgl(context, &(shell_vgl[0][0][0]),
        (int64_t) 5*point_num*shell_num);
  assert (rc == QMCKL_SUCCESS);

  printf(" shell_vgl[26][0][1]  %25.15e\n", shell_vgl[26][0][1]);
  printf(" shell_vgl[26][1][1]  %25.15e\n", shell_vgl[26][1][1]);
  printf(" shell_vgl[26][2][1]  %25.15e\n", shell_vgl[26][2][1]);
  printf(" shell_vgl[26][3][1]  %25.15e\n", shell_vgl[26][3][1]);
  printf(" shell_vgl[26][4][1]  %25.15e\n", shell_vgl[26][4][1]);

  assert( fabs(shell_vgl[26][0][1] - ( 3.564393437193868e-02)) < 1.e-14 );
  assert( fabs(shell_vgl[26][1][1] - (-6.030177987072189e-03)) < 1.e-14 );
  assert( fabs(shell_vgl[26][2][1] - (-3.074832579537582e-02)) < 1.e-14 );
  assert( fabs(shell_vgl[26][3][1] - ( 2.809546963519935e-02)) < 1.e-14 );
  assert( fabs(shell_vgl[26][4][1] - ( 1.896046117183968e-02)) < 1.e-14 );

}

{
#define shell_num chbrclf_shell_num

  double* elec_coord    = &(chbrclf_elec_coord[0][0][0]);

  assert(qmckl_electron_provided(context));

  int64_t point_num = elec_num;
  rc = qmckl_set_point(context, 'N', point_num, elec_coord, point_num*3);
  assert(rc == QMCKL_SUCCESS);


  double shell_hessian[point_num][3][4][shell_num];

  rc = qmckl_get_ao_basis_shell_hessian(context, &(shell_hessian[0][0][0][0]),
        (int64_t) 4*3*point_num*shell_num);
  assert (rc == QMCKL_SUCCESS);


  printf(" shell_hessian[26][0][0][1]  %25.15e\n", shell_hessian[26][0][0][1]);
  printf(" shell_hessian[26][1][0][1]  %25.15e\n", shell_hessian[26][1][0][1]);
  printf(" shell_hessian[26][2][0][1]  %25.15e\n", shell_hessian[26][2][0][1]);
  printf(" shell_hessian[26][0][1][1]  %25.15e\n", shell_hessian[26][0][1][1]);
  printf(" shell_hessian[26][1][1][1]  %25.15e\n", shell_hessian[26][1][1][1]);
  printf(" shell_hessian[26][2][1][1]  %25.15e\n", shell_hessian[26][2][1][1]);
  printf(" shell_hessian[26][0][2][1]  %25.15e\n", shell_hessian[26][0][2][1]);
  printf(" shell_hessian[26][1][2][1]  %25.15e\n", shell_hessian[26][1][2][1]);
  printf(" shell_hessian[26][2][2][1]  %25.15e\n", shell_hessian[26][2][2][1]);

  printf(" shell_hessian[26][0][3][1]  %25.15e\n", shell_hessian[26][0][3][1]);
  printf(" shell_hessian[26][1][3][1]  %25.15e\n", shell_hessian[26][1][3][1]);
  printf(" shell_hessian[26][2][3][1]  %25.15e\n", shell_hessian[26][2][3][1]);

  assert( fabs(shell_hessian[26][0][0][1] - ( -1.396360193576081e-02)) < 1.e-14 );
  assert( fabs(shell_hessian[26][1][0][1] - (  6.788393224947506e-03)) < 1.e-14 );
  assert( fabs(shell_hessian[26][2][0][1] - ( -6.202714807711193e-03)) < 1.e-14 );
  assert( fabs(shell_hessian[26][0][1][1] - (  6.788393224947506e-03)) < 1.e-14 );
  assert( fabs(shell_hessian[26][1][1][1] - (  1.931962058731147e-02)) < 1.e-14 );
  assert( fabs(shell_hessian[26][2][1][1] - ( -3.162810386893850e-02)) < 1.e-14 );
  assert( fabs(shell_hessian[26][0][2][1] - ( -6.202714807711193e-03)) < 1.e-14 );
  assert( fabs(shell_hessian[26][1][2][1] - ( -3.162810386893850e-02)) < 1.e-14 );
  assert( fabs(shell_hessian[26][2][2][1] - (  1.360444252028902e-02)) < 1.e-14 );
  //assert( fabs(shell_hessian[26][0][3][1] - (  3.122502256758253e+01)) < 1.e-14 );
  //assert( fabs(shell_hessian[26][1][3][1] - (  6.938893903907229e+00)) < 1.e-14 );
  //assert( fabs(shell_hessian[26][2][3][1] - ( -1.734723475976807e+01)) < 1.e-14 );
}

int  test_qmckl_ao_power(qmckl_context context);
assert(0 == test_qmckl_ao_power(context));

int  test_qmckl_ao_polynomial_vgl(qmckl_context context);
assert(0 == test_qmckl_ao_polynomial_vgl(context));

double X[3] = { 1.1, 2.2, 3.3 };
double R[3] = { 0.2, 1.1, 3.0 };
int32_t ldv[8] = {1, 4, 10, 20, 35, 56, 84, 120};
for (int32_t ldl=3 ; ldl<=5 ; ++ldl) {
    int64_t n;
    int32_t L0[200][ldl];
    int32_t L1[200][ldl];
    printf("ldl=%d\n", ldl);
    for (int32_t lmax=0 ; lmax<=7 ; lmax++) {
      double VGL0[5][ldv[lmax]];
      double VGL1[5][ldv[lmax]];
      memset(&L0[0][0], 0, sizeof(L0));
      memset(&L1[0][0], 0, sizeof(L1));
      memset(&VGL0[0][0], 0, sizeof(VGL0));
      memset(&VGL1[0][0], 0, sizeof(VGL1));
      rc = qmckl_ao_polynomial_transp_vgl_doc (context, X, R, lmax, &n, &(L0[0][0]), ldl, &(VGL0[0][0]), ldv[lmax]);
      assert (rc == QMCKL_SUCCESS);
      rc = qmckl_ao_polynomial_transp_vgl_hpc (context, X, R, lmax, &n, &(L1[0][0]), ldl, &(VGL1[0][0]), ldv[lmax]);
      assert (rc == QMCKL_SUCCESS);
      printf("lmax=%d\n", lmax);
      for (int32_t l=0 ; l<n ; ++l) {
        for (int32_t k=0 ; k<3 ; ++k) {
          printf("L[%d][%d] = %d   %d\n", l, k, L0[l][k], L1[l][k]);
          assert( L0[l][k] == L1[l][k] );
        }
      }

      for (int32_t k=0 ; k<5 ; ++k) {
        for (int32_t l=0 ; l<n ; ++l) {
          printf("VGL[%d][%d] = %e   %e  %e\n", k, l, VGL0[k][l], VGL1[k][l], VGL0[k][l]-VGL1[k][l]);
          assert( fabs(1.-(VGL0[k][l]+1.e-100)/(VGL1[k][l]+1.e-100)) < 1.e-15 );
        }
      }
    }
}

{
#define shell_num chbrclf_shell_num
#define ao_num chbrclf_ao_num

double* elec_coord    = &(chbrclf_elec_coord[0][0][0]);

assert(qmckl_electron_provided(context));

int64_t point_num = elec_num;
rc = qmckl_set_point(context, 'N', point_num, elec_coord, point_num*3);
assert(rc == QMCKL_SUCCESS);


double ao_value[point_num][ao_num];

rc = qmckl_get_ao_basis_ao_value(context, &(ao_value[0][0]),
         (int64_t) point_num*ao_num);
assert (rc == QMCKL_SUCCESS);

printf("\n");
printf(" ao_value ao_value[26][219] %25.15e\n", ao_value[26][219]);
printf(" ao_value ao_value[26][220] %25.15e\n", ao_value[26][220]);
printf(" ao_value ao_value[26][221] %25.15e\n", ao_value[26][221]);
printf(" ao_value ao_value[26][222] %25.15e\n", ao_value[26][222]);
printf(" ao_value ao_value[26][223] %25.15e\n", ao_value[26][223]);
printf(" ao_value ao_value[26][224] %25.15e\n", ao_value[26][224]);
printf("\n");

printf("%e %e\n", ao_value[26][219], 1.020298798341620e-08);
assert( fabs(ao_value[26][219] - (  1.020298798341620e-08)) < 1.e-14 );
assert( fabs(ao_value[26][220] - (  1.516643537739178e-08)) < 1.e-14 );
assert( fabs(ao_value[26][221] - ( -4.686370882518819e-09)) < 1.e-14 );
assert( fabs(ao_value[26][222] - (  7.514816980753531e-09)) < 1.e-14 );
assert( fabs(ao_value[26][223] - ( -4.021908374204471e-09)) < 1.e-14 );
assert( fabs(ao_value[26][224] - (  7.175045873560788e-10)) < 1.e-14 );

}

{
#define shell_num chbrclf_shell_num
#define ao_num chbrclf_ao_num

double* elec_coord    = &(chbrclf_elec_coord[0][0][0]);

assert(qmckl_electron_provided(context));

int64_t point_num = elec_num;
rc = qmckl_set_point(context, 'N', point_num, elec_coord, point_num*3);
assert(rc == QMCKL_SUCCESS);


double ao_vgl[point_num][5][ao_num];

rc = qmckl_get_ao_basis_ao_vgl(context, &(ao_vgl[0][0][0]),
         (int64_t) 5*point_num*ao_num);
assert (rc == QMCKL_SUCCESS);

printf("\n");
printf(" ao_vgl ao_vgl[26][0][219] %25.15e\n", ao_vgl[26][0][219]);
printf(" ao_vgl ao_vgl[26][1][219] %25.15e\n", ao_vgl[26][1][219]);
printf(" ao_vgl ao_vgl[26][2][219] %25.15e\n", ao_vgl[26][2][219]);
printf(" ao_vgl ao_vgl[26][3][219] %25.15e\n", ao_vgl[26][3][219]);
printf(" ao_vgl ao_vgl[26][4][219] %25.15e\n", ao_vgl[26][4][219]);
printf(" ao_vgl ao_vgl[26][0][220] %25.15e\n", ao_vgl[26][0][220]);
printf(" ao_vgl ao_vgl[26][1][220] %25.15e\n", ao_vgl[26][1][220]);
printf(" ao_vgl ao_vgl[26][2][220] %25.15e\n", ao_vgl[26][2][220]);
printf(" ao_vgl ao_vgl[26][3][220] %25.15e\n", ao_vgl[26][3][220]);
printf(" ao_vgl ao_vgl[26][4][220] %25.15e\n", ao_vgl[26][4][220]);
printf(" ao_vgl ao_vgl[26][0][221] %25.15e\n", ao_vgl[26][0][221]);
printf(" ao_vgl ao_vgl[26][1][221] %25.15e\n", ao_vgl[26][1][221]);
printf(" ao_vgl ao_vgl[26][2][221] %25.15e\n", ao_vgl[26][2][221]);
printf(" ao_vgl ao_vgl[26][3][221] %25.15e\n", ao_vgl[26][3][221]);
printf(" ao_vgl ao_vgl[26][4][221] %25.15e\n", ao_vgl[26][4][221]);
printf(" ao_vgl ao_vgl[26][0][222] %25.15e\n", ao_vgl[26][0][222]);
printf(" ao_vgl ao_vgl[26][1][222] %25.15e\n", ao_vgl[26][1][222]);
printf(" ao_vgl ao_vgl[26][2][222] %25.15e\n", ao_vgl[26][2][222]);
printf(" ao_vgl ao_vgl[26][3][222] %25.15e\n", ao_vgl[26][3][222]);
printf(" ao_vgl ao_vgl[26][4][222] %25.15e\n", ao_vgl[26][4][222]);
printf(" ao_vgl ao_vgl[26][0][223] %25.15e\n", ao_vgl[26][0][223]);
printf(" ao_vgl ao_vgl[26][1][223] %25.15e\n", ao_vgl[26][1][223]);
printf(" ao_vgl ao_vgl[26][2][223] %25.15e\n", ao_vgl[26][2][223]);
printf(" ao_vgl ao_vgl[26][3][223] %25.15e\n", ao_vgl[26][3][223]);
printf(" ao_vgl ao_vgl[26][4][223] %25.15e\n", ao_vgl[26][4][223]);
printf(" ao_vgl ao_vgl[26][0][224] %25.15e\n", ao_vgl[26][0][224]);
printf(" ao_vgl ao_vgl[26][1][224] %25.15e\n", ao_vgl[26][1][224]);
printf(" ao_vgl ao_vgl[26][2][224] %25.15e\n", ao_vgl[26][2][224]);
printf(" ao_vgl ao_vgl[26][3][224] %25.15e\n", ao_vgl[26][3][224]);
printf(" ao_vgl ao_vgl[26][4][224] %25.15e\n", ao_vgl[26][4][224]);
printf("\n");

assert( fabs(ao_vgl[26][0][219] - (  1.020298798341620e-08)) < 1.e-14 );
assert( fabs(ao_vgl[26][1][219] - ( -4.928035238010602e-08)) < 1.e-14 );
assert( fabs(ao_vgl[26][2][219] - ( -4.691009312035986e-08)) < 1.e-14 );
assert( fabs(ao_vgl[26][3][219] - (  1.449504046436699e-08)) < 1.e-14 );
assert( fabs(ao_vgl[26][4][219] - (  4.296442111843973e-07)) < 1.e-14 );
assert( fabs(ao_vgl[26][0][220] - (  1.516643537739178e-08)) < 1.e-14 );
assert( fabs(ao_vgl[26][1][220] - ( -7.725221462603871e-08)) < 1.e-14 );
assert( fabs(ao_vgl[26][2][220] - ( -6.507140835104833e-08)) < 1.e-14 );
assert( fabs(ao_vgl[26][3][220] - (  2.154644255710413e-08)) < 1.e-14 );
assert( fabs(ao_vgl[26][4][220] - (  6.365449359656352e-07)) < 1.e-14 );
assert( fabs(ao_vgl[26][0][221] - ( -4.686370882518819e-09)) < 1.e-14 );
assert( fabs(ao_vgl[26][1][221] - (  2.387064067626827e-08)) < 1.e-14 );
assert( fabs(ao_vgl[26][2][221] - (  2.154644255710412e-08)) < 1.e-14 );
assert( fabs(ao_vgl[26][3][221] - ( -1.998731863512374e-09)) < 1.e-14 );
assert( fabs(ao_vgl[26][4][221] - ( -1.966899656441993e-07)) < 1.e-14 );
assert( fabs(ao_vgl[26][0][222] - (  7.514816980753531e-09)) < 1.e-14 );
assert( fabs(ao_vgl[26][1][222] - ( -4.025889138635182e-08)) < 1.e-14 );
assert( fabs(ao_vgl[26][2][222] - ( -2.993372555126361e-08)) < 1.e-14 );
assert( fabs(ao_vgl[26][3][222] - (  1.067604670272904e-08)) < 1.e-14 );
assert( fabs(ao_vgl[26][4][222] - (  3.168199650002648e-07)) < 1.e-14 );
assert( fabs(ao_vgl[26][0][223] - ( -4.021908374204471e-09)) < 1.e-14 );
assert( fabs(ao_vgl[26][1][223] - (  2.154644255710413e-08)) < 1.e-14 );
assert( fabs(ao_vgl[26][2][223] - (  1.725594944732276e-08)) < 1.e-14 );
assert( fabs(ao_vgl[26][3][223] - ( -1.715339357718333e-09)) < 1.e-14 );
assert( fabs(ao_vgl[26][4][223] - ( -1.688020516893476e-07)) < 1.e-14 );
assert( fabs(ao_vgl[26][0][224] - (  7.175045873560788e-10)) < 1.e-14 );
assert( fabs(ao_vgl[26][1][224] - ( -3.843864637762753e-09)) < 1.e-14 );
assert( fabs(ao_vgl[26][2][224] - ( -3.298857850451910e-09)) < 1.e-14 );
assert( fabs(ao_vgl[26][3][224] - ( -4.073047518790881e-10)) < 1.e-14 );
assert( fabs(ao_vgl[26][4][224] - (  3.153244195820293e-08)) < 1.e-14 );

}

printf("AO hessian\n");

int64_t point_num = elec_num;

double delta_x = 0.00001;
double coef[9] = { 1.0/280.0, -4.0/105.0, 1.0/5.0, -4.0/5.0, 0.0, 4.0/5.0, -1.0/5.0, 4.0/105.0, -1.0/280.0 };

rc = qmckl_set_nucleus_coord(context, 'T', &(nucl_coord[0]), 3*nucl_num);
assert(rc == QMCKL_SUCCESS);

double * ao_hessian = (double*) malloc(ao_num*4*point_num*3 *sizeof(double));
rc = qmckl_get_ao_basis_ao_hessian(context, &ao_hessian[0], 3*4*ao_num*point_num);
assert(rc == QMCKL_SUCCESS);

double * finite_difference_ao_hessian= (double*) malloc(3*3 * point_num * ao_num * sizeof(double));

double* nucleus_coord = (double*) malloc(3 * nucl_num * sizeof(double));
if (nucleus_coord == NULL) {
  return QMCKL_ALLOCATION_FAILED;
}

rc = qmckl_get_nucleus_coord(context, 'N', nucleus_coord, 3*nucl_num);

double* temp_coord = (double*) malloc(3 * nucl_num * sizeof(double));
if (temp_coord == NULL) {
  free(nucleus_coord);
  return QMCKL_ALLOCATION_FAILED;
}

double* ao_vgl= (double*) malloc(5 * point_num * ao_num * sizeof(double));
if (ao_vgl == NULL) {
  free(temp_coord);
  free(nucleus_coord);
  return QMCKL_ALLOCATION_FAILED;
}


// Copy original coordinates
for (int i = 0; i < 3 * nucl_num; i++) {
  temp_coord[i] = nucleus_coord[i];
}



for (int i = 0; i < 3*3*ao_num*point_num; i++) {
  finite_difference_ao_hessian[i] = 0.0;
}
for (int64_t a = 0; a < nucl_num; a++) {
  for (int64_t k = 0; k < 3; k++) {
    for (int64_t m = -4; m <= 4; m++) {

      // Apply finite difference displacement
      temp_coord[k+a*3] = nucleus_coord[k+3*a] + (double) m * delta_x;

      // Update coordinates in the context
      rc = qmckl_set_nucleus_coord(context, 'N', temp_coord, 3*nucl_num);
      assert(rc == QMCKL_SUCCESS);

      rc = qmckl_context_touch(context);
      assert(rc == QMCKL_SUCCESS);

      // Call the provided function
      rc = qmckl_get_ao_basis_ao_vgl(context,&ao_vgl[0], 5*point_num*ao_num);
      assert(rc == QMCKL_SUCCESS);

      // Accumulate derivative using finite-difference coefficients
      for (int i = 0; i < point_num; i++) {
        for (int n = 0; n < 3; n++){
          for (int j = 0; j < ao_num; j++) {
            finite_difference_ao_hessian[k*3*ao_num*point_num  + i*3*ao_num + n*ao_num + j] -=
              coef[m + 4] * ao_vgl[i*ao_num*5 + (n+1)*ao_num + j]/delta_x;
          }
        }
      }
    }
    temp_coord[k+a*3] = nucleus_coord[k+3*a];
  }
}

// Reset coordinates in the context
rc = qmckl_set_nucleus_coord(context, 'N', temp_coord, 3*nucl_num);
assert(rc == QMCKL_SUCCESS);

rc = qmckl_context_touch(context);
assert(rc == QMCKL_SUCCESS);

free(nucleus_coord);
free(temp_coord);
free(ao_vgl);




  for (int n = 0; n < 3; n++){
    for (int j = 0; j < point_num; j++){
      for (int m = 0; m < 3; m++){
        for (int i = 0; i < ao_num; i++){
          //printf("n=%i j=%i m=%i i=%i\n", n, j, m, i);
          //printf("%.10f\n", ao_hessian[n*3*ao_num*point_num + j*3*ao_num + m*ao_num + i]);
          //printf("%.10f\n", finite_difference_ao_hessian[n*3*ao_num*point_num + j*3*ao_num + m*ao_num + i]);
          assert(fabs(finite_difference_ao_hessian[n*3*ao_num*point_num + j*3*ao_num + m*ao_num + i] - 
                      ao_hessian[n*4*ao_num*point_num + j*4*ao_num + m*ao_num + i]) < 1.e-9);
        }
      }
    }
  }


free(ao_hessian);
free(finite_difference_ao_hessian);

printf("OK\n");

rc = qmckl_context_destroy(context);
    assert (rc == QMCKL_SUCCESS);

    return 0;
}
