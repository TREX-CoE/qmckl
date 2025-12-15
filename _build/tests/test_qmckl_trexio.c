#include "qmckl.h"
#include <assert.h>
#include <math.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "chbrclf.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main() {
  qmckl_context context;
  context = qmckl_context_create();
#ifdef HAVE_TREXIO

#ifdef HAVE_TREXIO

qmckl_exit_code rc;
char filename[256];

#ifndef QMCKL_TEST_DIR
#error "QMCKL_TEST_DIR is not defined"
#endif

strncpy(filename, QMCKL_TEST_DIR,255);
strncat(filename, "/chbrclf", 255);
printf("Test file: %s\n", filename);

rc = qmckl_trexio_read(context, filename, 255);

qmckl_check(context, rc);

printf("Electrons\n");
int64_t up_num, dn_num;
rc = qmckl_get_electron_up_num(context, &up_num);
qmckl_check(context, rc);
assert (up_num == chbrclf_elec_up_num);

rc = qmckl_get_electron_down_num(context, &dn_num);
qmckl_check(context, rc);
assert (dn_num == chbrclf_elec_dn_num);

printf("Nuclei\n");

int64_t nucl_num;
rc = qmckl_get_nucleus_num(context, &nucl_num);
qmckl_check(context, rc);
assert (nucl_num == chbrclf_nucl_num);

printf("Nuclear charges\n");
double * charge = (double*) malloc (nucl_num * sizeof(double));
rc = qmckl_get_nucleus_charge(context, charge, nucl_num);
qmckl_check(context, rc);
for (int i=0 ; i<nucl_num ; i++) {
  assert (charge[i] == chbrclf_charge[i]);
 }
free(charge);
charge = NULL;

printf("Nuclear coordinates\n");
double * coord  = (double*) malloc (nucl_num * 3 * sizeof(double));
rc = qmckl_get_nucleus_coord(context, 'T', coord, 3*nucl_num);
qmckl_check(context, rc);
int k=0;
for (int j=0 ; j<3 ; ++j) {
  for (int i=0 ; i<nucl_num ; ++i) {
//    printf("%f  %f\n", coord[k], chbrclf_nucl_coord[j][i]);
    assert (coord[k] == chbrclf_nucl_coord[j][i]);
    ++k;
  }
 }
free(coord);
coord = NULL;

printf("Atomic basis\n");

char basis_type;
rc = qmckl_get_ao_basis_type(context, &basis_type);
assert (basis_type == 'G');

int64_t shell_num;
rc = qmckl_get_ao_basis_shell_num(context, &shell_num);
assert (shell_num == chbrclf_shell_num);

int64_t prim_num;
rc = qmckl_get_ao_basis_prim_num(context, &prim_num);
assert (prim_num == chbrclf_prim_num);

int64_t ao_num;
rc = qmckl_get_ao_basis_ao_num(context, &ao_num);
assert (ao_num == chbrclf_ao_num);

int64_t* nucleus_index = (int64_t*) malloc (nucl_num * sizeof(int64_t));
rc = qmckl_get_ao_basis_nucleus_index(context, nucleus_index, nucl_num);
qmckl_check(context, rc);
for (int i=0 ; i<nucl_num ; i++) {
  assert (nucleus_index[i] == chbrclf_basis_nucleus_index[i]);
 }
free(nucleus_index);
nucleus_index = NULL;

int64_t* nucleus_shell_num  = (int64_t*) malloc (nucl_num * sizeof(int64_t));
rc = qmckl_get_ao_basis_nucleus_shell_num(context, nucleus_shell_num, nucl_num);
qmckl_check(context, rc);
for (int i=0 ; i<nucl_num ; i++) {
  assert (nucleus_shell_num[i] == chbrclf_basis_nucleus_shell_num[i]);
 }
free(nucleus_shell_num);
nucleus_shell_num = NULL;

int32_t* r_power  = (int32_t*) malloc (shell_num * sizeof(int32_t));
rc = qmckl_get_ao_basis_r_power(context, r_power, shell_num);
qmckl_check(context, rc);
for (int i=0 ; i<shell_num ; i++) {
  assert (r_power[i] == 0);
 }
free(r_power);
r_power = NULL;

int32_t* shell_ang_mom  = (int32_t*) malloc (shell_num * sizeof(int32_t));
rc = qmckl_get_ao_basis_shell_ang_mom(context, shell_ang_mom, shell_num);
qmckl_check(context, rc);
for (int i=0 ; i<shell_num ; i++) {
  assert (shell_ang_mom[i] == chbrclf_basis_shell_ang_mom[i]);
 }
free(shell_ang_mom);
shell_ang_mom = NULL;

int64_t* shell_prim_num  = (int64_t*) malloc (shell_num * sizeof(int64_t));
rc = qmckl_get_ao_basis_shell_prim_num(context, shell_prim_num, shell_num);
qmckl_check(context, rc);
for (int i=0 ; i<shell_num ; i++) {
  assert (shell_prim_num[i] == chbrclf_basis_shell_prim_num[i]);
 }
free(shell_prim_num);
shell_prim_num = NULL;

int64_t* shell_prim_index  = (int64_t*) malloc (shell_num * sizeof(int64_t));
rc = qmckl_get_ao_basis_shell_prim_index(context, shell_prim_index, shell_num);
qmckl_check(context, rc);
for (int i=0 ; i<shell_num ; i++) {
  assert (shell_prim_index[i] == chbrclf_basis_shell_prim_index[i]);
 }
free(shell_prim_index);
shell_prim_index = NULL;

double* shell_factor = (double*) malloc (shell_num * sizeof(double));
rc = qmckl_get_ao_basis_shell_factor(context, shell_factor, shell_num);
qmckl_check(context, rc);
for (int i=0 ; i<shell_num ; i++) {
  assert (fabs(shell_factor[i] - chbrclf_basis_shell_factor[i]) < 1.e-6);
 }
free(shell_factor);
shell_factor = NULL;

double* exponent = (double*) malloc (prim_num * sizeof(double));
rc = qmckl_get_ao_basis_exponent(context, exponent, prim_num);
qmckl_check(context, rc);
for (int i=0 ; i<prim_num ; i++) {
  assert (fabs((exponent[i] - chbrclf_basis_exponent[i])/chbrclf_basis_exponent[i]) < 1.e-7);
 }
free(exponent);
exponent = NULL;

double* coefficient = (double*) malloc (prim_num * sizeof(double));
rc = qmckl_get_ao_basis_coefficient(context, coefficient, prim_num);
qmckl_check(context, rc);
for (int i=0 ; i<prim_num ; i++) {
  assert (fabs((coefficient[i] - chbrclf_basis_coefficient[i])/chbrclf_basis_coefficient[i]) < 1.e-7);
 }
free(coefficient);
coefficient = NULL;

double* prim_factor = (double*) malloc (prim_num * sizeof(double));
rc = qmckl_get_ao_basis_prim_factor(context, prim_factor, prim_num);
qmckl_check(context, rc);
for (int i=0 ; i<prim_num ; i++) {
  assert (fabs((prim_factor[i] - chbrclf_basis_prim_factor[i])/chbrclf_basis_prim_factor[i]) < 1.e-7);
 }
free(prim_factor);
prim_factor = NULL;

#endif

printf("MOs\n");

int64_t mo_num;
rc = qmckl_get_mo_basis_mo_num(context, &mo_num);
qmckl_check(context, rc);

assert (mo_num == chbrclf_mo_num);

printf("MO coefs\n");
double * mo_coef = (double*) malloc (ao_num * mo_num * sizeof(double));
rc = qmckl_get_mo_basis_coefficient(context, mo_coef, mo_num*ao_num);
qmckl_check(context, rc);
for (int i=0 ; i<ao_num * mo_num ; i++) {
  printf("%d %e %e %e\n", i, mo_coef[i], chbrclf_mo_coef[i],
         ( fabs(mo_coef[i] - chbrclf_mo_coef[i])/fabs(mo_coef[i])) );
  assert ( fabs(mo_coef[i] - chbrclf_mo_coef[i])/fabs(mo_coef[i]) < 1.e-12 );
 }
free(mo_coef);
charge = NULL;

if (qmckl_context_destroy(context) != QMCKL_SUCCESS)
    return QMCKL_FAILURE;
#endif
  return 0;
}
