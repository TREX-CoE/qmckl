#include "qmckl.h"
#include "assert.h"
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <math.h>
#include "chbrclf.h"
#include "qmckl_ao_private_func.h"
#include "qmckl_mo_private_func.h"
#include "qmckl_determinant_private_func.h"

int main() {
    qmckl_context context;
    context = qmckl_context_create();

    qmckl_exit_code rc;

double* elec_coord    = &(chbrclf_elec_coord[0][0][0]);
const double*   nucl_charge   = chbrclf_charge;
const double*   nucl_coord    = &(chbrclf_nucl_coord[0][0]);

rc = qmckl_set_electron_num (context, chbrclf_elec_up_num, chbrclf_elec_dn_num);
qmckl_check(context, rc);

rc = qmckl_set_electron_coord (context, 'N', chbrclf_walk_num, elec_coord, chbrclf_walk_num*chbrclf_elec_num*3);
qmckl_check(context, rc);

assert(qmckl_electron_provided(context));

rc = qmckl_set_nucleus_num (context, chbrclf_nucl_num);
qmckl_check(context, rc);

rc = qmckl_set_nucleus_coord (context, 'T', &(nucl_coord[0]), chbrclf_nucl_num*3);
qmckl_check(context, rc);

rc = qmckl_set_nucleus_charge(context, nucl_charge, chbrclf_nucl_num);
qmckl_check(context, rc);

assert(qmckl_nucleus_provided(context));

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
qmckl_check(context, rc);
assert(!qmckl_ao_basis_provided(context));

rc = qmckl_set_ao_basis_shell_num (context, chbrclf_shell_num);
qmckl_check(context, rc);
assert(!qmckl_ao_basis_provided(context));

rc = qmckl_set_ao_basis_prim_num (context, chbrclf_prim_num);
qmckl_check(context, rc);
assert(!qmckl_ao_basis_provided(context));

rc = qmckl_set_ao_basis_nucleus_index (context, nucleus_index, chbrclf_nucl_num);
qmckl_check(context, rc);
assert(!qmckl_ao_basis_provided(context));

rc = qmckl_set_ao_basis_nucleus_shell_num (context, nucleus_shell_num, chbrclf_nucl_num);
qmckl_check(context, rc);
assert(!qmckl_ao_basis_provided(context));

rc = qmckl_set_ao_basis_shell_ang_mom (context, shell_ang_mom, chbrclf_shell_num);
qmckl_check(context, rc);
assert(!qmckl_ao_basis_provided(context));

rc = qmckl_set_ao_basis_shell_factor  (context, shell_factor, chbrclf_shell_num);
qmckl_check(context, rc);
assert(!qmckl_ao_basis_provided(context));

rc = qmckl_set_ao_basis_shell_prim_num (context, shell_prim_num, chbrclf_shell_num);
qmckl_check(context, rc);
assert(!qmckl_ao_basis_provided(context));

rc = qmckl_set_ao_basis_shell_prim_index (context, shell_prim_index, chbrclf_shell_num);
qmckl_check(context, rc);
assert(!qmckl_ao_basis_provided(context));

rc = qmckl_set_ao_basis_exponent (context, exponent, chbrclf_prim_num);
qmckl_check(context, rc);
assert(!qmckl_ao_basis_provided(context));

rc = qmckl_set_ao_basis_coefficient (context, coefficient, chbrclf_prim_num);
qmckl_check(context, rc);
assert(!qmckl_ao_basis_provided(context));

rc = qmckl_set_ao_basis_prim_factor (context, prim_factor, chbrclf_prim_num);
qmckl_check(context, rc);

rc = qmckl_set_ao_basis_ao_num(context, chbrclf_ao_num);
qmckl_check(context, rc);

rc = qmckl_set_ao_basis_ao_factor (context, ao_factor, chbrclf_ao_num);
qmckl_check(context, rc);

assert(qmckl_ao_basis_provided(context));


double ao_vgl[chbrclf_walk_num*chbrclf_elec_num][5][chbrclf_ao_num];

rc = qmckl_get_ao_basis_ao_vgl(context, &(ao_vgl[0][0][0]), (int64_t) 5*chbrclf_walk_num*chbrclf_elec_num*chbrclf_ao_num);
qmckl_check(context, rc);

/* Set up MO data */
rc = qmckl_set_mo_basis_mo_num(context, chbrclf_mo_num);
qmckl_check(context, rc);

const double  * mo_coefficient          =  &(chbrclf_mo_coef[0]);

rc = qmckl_set_mo_basis_coefficient(context, mo_coefficient, chbrclf_mo_num*chbrclf_ao_num);
qmckl_check(context, rc);

assert(qmckl_mo_basis_provided(context));

double mo_vgl[chbrclf_walk_num*chbrclf_elec_num][5][chbrclf_mo_num];
rc = qmckl_get_mo_basis_mo_vgl(context, &(mo_vgl[0][0][0]), 5*chbrclf_walk_num*chbrclf_elec_num*chbrclf_mo_num);
qmckl_check(context, rc);

/* Set up determinant data */

#define det_num_alpha 1
#define det_num_beta  1
int64_t mo_index_alpha[det_num_alpha][chbrclf_walk_num][chbrclf_elec_up_num];
int64_t mo_index_beta[det_num_alpha][chbrclf_walk_num][chbrclf_elec_dn_num];

int i, j, k;
for(k = 0; k < det_num_alpha; ++k)
  for(i = 0; i < chbrclf_walk_num; ++i)
    for(j = 0; j < chbrclf_elec_up_num; ++j)
      mo_index_alpha[k][i][j] = j + 1;
for(k = 0; k < det_num_beta; ++k)
  for(i = 0; i < chbrclf_walk_num; ++i)
    for(j = 0; j < chbrclf_elec_up_num; ++j)
      mo_index_beta[k][i][j] = j + 1;

rc = qmckl_set_determinant_type (context, typ);
qmckl_check(context, rc);

rc = qmckl_set_determinant_det_num_alpha (context, det_num_alpha);
qmckl_check(context, rc);

rc = qmckl_set_determinant_det_num_beta (context, det_num_beta);
qmckl_check(context, rc);

rc = qmckl_set_determinant_mo_index_alpha (context, &(mo_index_alpha[0][0][0]), det_num_alpha*chbrclf_walk_num*chbrclf_elec_up_num);
qmckl_check(context, rc);

rc = qmckl_set_determinant_mo_index_beta (context, &(mo_index_beta[0][0][0]), det_num_beta*chbrclf_walk_num*chbrclf_elec_dn_num);
qmckl_check(context, rc);

// Get slater-determinant

double det_vgl_alpha[det_num_alpha][chbrclf_walk_num][5][chbrclf_elec_up_num][chbrclf_elec_up_num];
double det_vgl_beta[det_num_beta][chbrclf_walk_num][5][chbrclf_elec_dn_num][chbrclf_elec_dn_num];

rc = qmckl_get_det_vgl_alpha(context, &(det_vgl_alpha[0][0][0][0][0]));
qmckl_check(context, rc);

rc = qmckl_get_det_vgl_beta(context, &(det_vgl_beta[0][0][0][0][0]));
qmckl_check(context, rc);

// Get adjoint of the slater-determinant

double det_inv_matrix_alpha[det_num_alpha][chbrclf_walk_num][chbrclf_elec_up_num][chbrclf_elec_up_num];
double det_inv_matrix_beta[det_num_beta][chbrclf_walk_num][chbrclf_elec_dn_num][chbrclf_elec_dn_num];

rc = qmckl_get_det_inv_matrix_alpha(context, &(det_inv_matrix_alpha[0][0][0][0]));
qmckl_check(context, rc);

rc = qmckl_get_det_inv_matrix_beta(context, &(det_inv_matrix_beta[0][0][0][0]));
qmckl_check(context, rc);

rc = qmckl_context_destroy(context);
    qmckl_check(context, rc);

    return 0;
}
