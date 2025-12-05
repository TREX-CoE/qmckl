#include "qmckl.h"
#include <assert.h>
#include <math.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdbool.h>
#include <stdio.h>
#include "n2.h"
#include "chbrclf.h"
#include "qmckl_context_private_type.h"
#include "qmckl_jastrow_champ_private_func.h"
#include "qmckl_jastrow_champ_single_private_func.h"
#include "qmckl_forces_private_func.h"

int main() {
  qmckl_context context;
  context = qmckl_context_create();

/* Reference input data */
int64_t walk_num      = n2_walk_num;
int64_t elec_num      = n2_elec_num;
int64_t elec_up_num   = n2_elec_up_num;
int64_t elec_dn_num   = n2_elec_dn_num;
int64_t nucl_num      = n2_nucl_num;
double  rescale_factor_ee   = 0.6;
double  rescale_factor_en[2] = { 0.6, 0.6 };
double* elec_coord    = &(n2_elec_coord[0][0][0]);

const double*   nucl_charge   = n2_charge;
double*  nucl_coord    = &(n2_nucl_coord[0][0]);

/* Provide Electron data */

qmckl_exit_code rc;

assert(!qmckl_electron_provided(context));

rc = qmckl_check(context,
                 qmckl_set_electron_num (context, elec_up_num, elec_dn_num)
                 );
assert(rc == QMCKL_SUCCESS);

assert(qmckl_electron_provided(context));

rc = qmckl_check(context,
                 qmckl_set_electron_coord (context, 'N', walk_num, elec_coord, walk_num*3*elec_num)
                 );
assert(rc == QMCKL_SUCCESS);

double elec_coord2[walk_num*3*elec_num];

rc = qmckl_check(context,
                 qmckl_get_electron_coord (context, 'N', elec_coord2, walk_num*3*elec_num)
                 );
assert(rc == QMCKL_SUCCESS);
for (int64_t i=0 ; i<3*elec_num ; ++i) {
  assert( elec_coord[i] == elec_coord2[i] );
 }


/* Provide Nucleus data */

assert(!qmckl_nucleus_provided(context));

rc = qmckl_check(context,
                 qmckl_set_nucleus_num (context, nucl_num)
                 );
assert(rc == QMCKL_SUCCESS);
assert(!qmckl_nucleus_provided(context));

double nucl_coord2[3*nucl_num];

rc = qmckl_get_nucleus_coord (context, 'T', nucl_coord2, 3*nucl_num);
assert(rc == QMCKL_NOT_PROVIDED);

rc = qmckl_check(context,
                 qmckl_set_nucleus_coord (context, 'T', &(nucl_coord[0]), 3*nucl_num)
                 );
assert(rc == QMCKL_SUCCESS);

assert(!qmckl_nucleus_provided(context));

rc = qmckl_check(context,
                 qmckl_get_nucleus_coord (context, 'N', nucl_coord2, nucl_num*3)
                 );
assert(rc == QMCKL_SUCCESS);
for (int64_t k=0 ; k<3 ; ++k) {
  for (int64_t i=0 ; i<nucl_num ; ++i) {
    assert( nucl_coord[nucl_num*k+i] == nucl_coord2[3*i+k] );
  }
}

rc = qmckl_check(context,
                 qmckl_get_nucleus_coord (context, 'T', nucl_coord2, nucl_num*3)
                 );
assert(rc == QMCKL_SUCCESS);
for (int64_t i=0 ; i<3*nucl_num ; ++i) {
  assert( nucl_coord[i] == nucl_coord2[i] );
}

double nucl_charge2[nucl_num];

rc = qmckl_get_nucleus_charge(context, nucl_charge2, nucl_num);
assert(rc == QMCKL_NOT_PROVIDED);

rc = qmckl_check(context,
                 qmckl_set_nucleus_charge(context, nucl_charge, nucl_num)
                 );
assert(rc == QMCKL_SUCCESS);

rc = qmckl_check(context,
                 qmckl_get_nucleus_charge(context, nucl_charge2, nucl_num)
                 );
assert(rc == QMCKL_SUCCESS);
for (int64_t i=0 ; i<nucl_num ; ++i) {
  assert( nucl_charge[i] == nucl_charge2[i] );
 }
assert(qmckl_nucleus_provided(context));

assert(qmckl_electron_provided(context));

int64_t type_nucl_num = n2_type_nucl_num;
int64_t* type_nucl_vector = &(n2_type_nucl_vector[0]);
int64_t aord_num = n2_aord_num;
int64_t bord_num = n2_bord_num;
int64_t cord_num = n2_cord_num;
double* a_vector = &(n2_a_vector[0][0]);
double* b_vector = &(n2_b_vector[0]);
double* c_vector = &(n2_c_vector[0][0]);
int64_t dim_c_vector=0;

assert(!qmckl_jastrow_champ_provided(context));

/* Set the data */
rc = qmckl_check(context,
                 qmckl_set_jastrow_champ_spin_independent(context, 0)
                 );

rc = qmckl_check(context,
                 qmckl_set_jastrow_champ_aord_num(context, aord_num)
                 );
rc = qmckl_check(context,
                 qmckl_set_jastrow_champ_bord_num(context, bord_num)
                 );
rc = qmckl_check(context,
                 qmckl_set_jastrow_champ_cord_num(context, cord_num)
                 );
assert(rc == QMCKL_SUCCESS);
rc = qmckl_check(context,
                 qmckl_set_jastrow_champ_type_nucl_num(context, type_nucl_num)
                 );
assert(rc == QMCKL_SUCCESS);
rc = qmckl_check(context,
                 qmckl_set_jastrow_champ_type_nucl_vector(context, type_nucl_vector, nucl_num)
                 );
assert(rc == QMCKL_SUCCESS);
rc = qmckl_check(context,
                 qmckl_set_jastrow_champ_a_vector(context, a_vector,(aord_num+1)*type_nucl_num)
                 );
assert(rc == QMCKL_SUCCESS);
rc = qmckl_check(context,
                 qmckl_set_jastrow_champ_b_vector(context, b_vector,(bord_num+1))
                 );
assert(rc == QMCKL_SUCCESS);
rc = qmckl_check(context,
                 qmckl_get_jastrow_champ_dim_c_vector(context, &dim_c_vector)
                 );
assert(rc == QMCKL_SUCCESS);
rc = qmckl_check(context,
                 qmckl_set_jastrow_champ_c_vector(context, c_vector, dim_c_vector*type_nucl_num)
                 );
assert(rc == QMCKL_SUCCESS);

double k_ee = 0.;
double k_en[2] = { 0., 0. };
rc = qmckl_check(context,
                 qmckl_set_jastrow_champ_rescale_factor_en(context, rescale_factor_en, type_nucl_num)
                 );
assert(rc == QMCKL_SUCCESS);

rc = qmckl_check(context,
                 qmckl_set_jastrow_champ_rescale_factor_ee(context, rescale_factor_ee)
                 );
assert(rc == QMCKL_SUCCESS);

rc = qmckl_check(context,
                 qmckl_get_jastrow_champ_rescale_factor_ee (context, &k_ee)
                 );
assert(rc == QMCKL_SUCCESS);
assert(k_ee == rescale_factor_ee);

rc = qmckl_check(context,
                 qmckl_get_jastrow_champ_rescale_factor_en (context, &(k_en[0]), type_nucl_num)
                 );
assert(rc == QMCKL_SUCCESS);
for (int i=0 ; i<type_nucl_num ; ++i) {
  assert(k_en[i] == rescale_factor_en[i]);
 }

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

double delta_x = 0.00001;
double coef[9] = { 1.0/280.0, -4.0/105.0, 1.0/5.0, -4.0/5.0, 0.0, 4.0/5.0, -1.0/5.0, 4.0/105.0, -1.0/280.0 };

printf("Forces Jastrow en\n");

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

rc = qmckl_set_nucleus_coord(context, 'T', &(nucl_coord[0]), 3*nucl_num);
assert(rc == QMCKL_SUCCESS);

double forces_jastrow_en[walk_num][nucl_num][3];
rc = qmckl_get_forces_jastrow_en(context, &forces_jastrow_en[0][0][0], 3*nucl_num*walk_num);
assert(rc == QMCKL_SUCCESS);

double finite_difference_force_en[walk_num][nucl_num][3];
rc = qmckl_finite_difference_deriv_n(context, delta_x, &qmckl_get_jastrow_champ_factor_en, &(finite_difference_force_en[0][0][0]), 1);

for (int nw = 0; nw < walk_num; nw++){
  for (int a = 0; a < nucl_num; a++) {
    for (int k = 0; k < 3; k++){
      //printf("%.10f\t", finite_difference_force_en[nw][a][k]);
      //printf("%.10f\n", forces_jastrow_en[nw][a][k]);
    }
  }
}
for (int nw = 0; nw < walk_num; nw++){
  for (int a = 0; a < nucl_num; a++) {
    for (int k = 0; k < 3; k++){
      assert(fabs(finite_difference_force_en[nw][a][k] - forces_jastrow_en[nw][a][k]) < 1.e-8);
    }
  }
}
printf("OK\n");

printf("Forces Jastrow en G\n");

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

rc = qmckl_set_nucleus_coord(context, 'T', &(nucl_coord[0]), 3*nucl_num);
assert(rc == QMCKL_SUCCESS);

double forces_jastrow_en_g[walk_num][nucl_num][3][elec_num][3];
rc = qmckl_get_forces_jastrow_en_g(context, &forces_jastrow_en_g[0][0][0][0][0], 3*3*nucl_num*walk_num*elec_num);
assert(rc == QMCKL_SUCCESS);

double finite_difference_force_en_g[walk_num][nucl_num][3][4][elec_num];
rc = qmckl_finite_difference_deriv_n(context, delta_x, &qmckl_get_jastrow_champ_factor_en_gl, &finite_difference_force_en_g[0][0][0][0][0], 4*elec_num);

for (int nw = 0; nw < walk_num; nw++){
  for (int a = 0; a < nucl_num; a++) {
    for (int k = 0; k < 3; k++){
      for (int i = 0; i < elec_num; i++){
        for (int l = 1; l < 3; l++){
          //printf("finite_difference_force_en_g[%i][%i][%i][%i][%i] %+3.10f \n", nw,a,k,l,i,finite_difference_force_en_g[nw][a][k][l][i]);
          //printf("forces_jastrow_en_g         [%i][%i][%i][%i][%i] %+3.10f\n", nw,a,k,i,l,forces_jastrow_en_g[nw][a][k][i][l]);
        }
      }
    }
  }
}
for (int nw = 0; nw < walk_num; nw++){
  for (int a = 0; a < nucl_num; a++) {
    for (int k = 0; k < 3; k++){
      for (int i = 0; i < elec_num; i++){
        for (int l = 1; l < 3; l++){
          assert(fabs(finite_difference_force_en_g[nw][a][k][l][i] - forces_jastrow_en_g[nw][a][k][i][l]) < 1.e-8);
        }
      }
    }
  }
}
printf("OK\n");

printf("Forces Jastrow en L\n");

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

rc = qmckl_set_nucleus_coord(context, 'T', &(nucl_coord[0]), 3*nucl_num);
assert(rc == QMCKL_SUCCESS);

double forces_jastrow_en_l[walk_num][nucl_num][3];
rc = qmckl_get_forces_jastrow_en_l(context, &forces_jastrow_en_l[0][0][0], 3*nucl_num*walk_num);
assert(rc == QMCKL_SUCCESS);

double finite_difference_force_en_l[walk_num][nucl_num][3][4][elec_num];
rc = qmckl_finite_difference_deriv_n(context, delta_x, &qmckl_get_jastrow_champ_factor_en_gl, &finite_difference_force_en_l[0][0][0][0][0], 4*elec_num);


double finite_difference_force_en_l_sum[walk_num][nucl_num][3];
for (int nw = 0; nw < walk_num; nw++){
  for (int a = 0; a < nucl_num; a++) {
    for (int k = 0; k < 3; k++){
      finite_difference_force_en_l_sum[nw][a][k] = 0;
      for (int i = 0; i < elec_num; i++){
        finite_difference_force_en_l_sum[nw][a][k] = finite_difference_force_en_l_sum[nw][a][k] + finite_difference_force_en_l[nw][a][k][3][i];
      }
    }
  }
}

for (int nw = 0; nw < walk_num; nw++){
  for (int a = 0; a < nucl_num; a++) {
    for (int k = 0; k < 3; k++){
      //printf("finite_difference_force_en_l_sum[%i][%i][%i] %+3.10f \n", nw,a,k,finite_difference_force_en_l_sum[nw][a][k]);
      //printf("forces_jastrow_en_l             [%i][%i][%i] %+3.10f\n", nw,a,k,forces_jastrow_en_l[nw][a][k]);

    }
  }
}
for (int nw = 0; nw < walk_num; nw++){
  for (int a = 0; a < nucl_num; a++) {
    for (int k = 0; k < 3; k++){
          assert(fabs(finite_difference_force_en_l_sum[nw][a][k] - forces_jastrow_en_l[nw][a][k]) < 1.e-8);
    }
  }
}
printf("OK\n");

printf("Forces Jastrow single en\n");

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

rc = qmckl_set_nucleus_coord(context, 'T', &(nucl_coord[0]), 3*nucl_num);
assert(rc == QMCKL_SUCCESS);

double new_coords[6] = {1.0,2.0,3.0,4.0,5.0,6.0};

rc = qmckl_set_single_point(context, 'N', 2, new_coords, 3*walk_num);
assert (rc == QMCKL_SUCCESS);

double forces_jastrow_single_en[walk_num][nucl_num][3];
rc = qmckl_get_forces_jastrow_single_en(context, &forces_jastrow_single_en[0][0][0], 3*nucl_num*walk_num);
assert(rc == QMCKL_SUCCESS);

double finite_difference_force_single_en[walk_num][nucl_num][3];
rc = qmckl_finite_difference_deriv_n(context, delta_x, &qmckl_get_jastrow_champ_single_en, &(finite_difference_force_single_en[0][0][0]), 1);

for (int nw = 0; nw < walk_num; nw++){
  for (int a = 0; a < nucl_num; a++) {
    for (int k = 0; k < 3; k++){
      //printf("%.10f\t", finite_difference_force_single_en[nw][a][k]);
      //printf("%.10f\n", forces_jastrow_single_en[nw][a][k]);
      assert(fabs(finite_difference_force_single_en[nw][a][k] - forces_jastrow_single_en[nw][a][k]) < 1.e-8);
    }
  }
}
printf("OK\n");

printf("Forces tmp_c\n");

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

rc = qmckl_set_nucleus_coord(context, 'T', &(nucl_coord[0]), 3*nucl_num);
assert(rc == QMCKL_SUCCESS);

double forces_tmp_c[walk_num][cord_num][cord_num+1][nucl_num][4][elec_num];
rc = qmckl_get_forces_tmp_c(context, &forces_tmp_c[0][0][0][0][0][0], 4*nucl_num*walk_num*elec_num*(cord_num+1)*cord_num);
assert(rc == QMCKL_SUCCESS);


double finite_difference_force_tmp_c[walk_num][cord_num][cord_num+1][nucl_num][3][elec_num];

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
double output[walk_num][cord_num][cord_num+1][nucl_num][elec_num];

// Copy original coordinates
for (int i = 0; i < 3 * nucl_num; i++) {
  temp_coord[i] = nucleus_coord[i];
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
      rc = qmckl_get_jastrow_champ_tmp_c(context,
           &output[0][0][0][0][0],
           4*nucl_num*walk_num*elec_num*(cord_num+1)*cord_num);
      assert(rc == QMCKL_SUCCESS);

      // Accumulate derivative using finite-difference coefficients
      for (int nw=0 ; nw<walk_num ; nw++) {
        for (int l = 0; l < cord_num; l++) {
          for (int mm = 0; mm <= cord_num; mm++) {
            for (int i = 0; i < elec_num; i++) {

              if (m == -4) {
                finite_difference_force_tmp_c[nw][l][mm][a][k][i] = 0.0;
              }
              finite_difference_force_tmp_c[nw][l][mm][a][k][i] += coef[m + 4] * output[nw][l][mm][a][i]/delta_x;
            }
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


for (int nw = 0; nw < walk_num; nw++){
  for (int l = 0; l < cord_num; l++){
    for (int m = 0; m <= cord_num; m++){
      for (int a = 0; a < nucl_num; a++) {
        for (int k = 0; k < 3; k++){
          for (int i = 0; i < elec_num; i++){
            //printf("nw=%i l=%i m=%i a=%i k=%i i=%i\n",nw,l,m,a,k,i);
            //printf("%.10f\t", finite_difference_force_tmp_c[nw][l][m][a][k][i]);
            //printf("%.10f\n", forces_tmp_c[nw][l][m][a][k][i]);

          }
        }
      }
    }
  }
}
for (int nw = 0; nw < walk_num; nw++){
  for (int l = 0; l < cord_num; l++){
    for (int m = 0; m <= cord_num; m++){
      for (int a = 0; a < nucl_num; a++) {
        for (int k = 0; k < 3; k++){
          for (int i = 0; i < elec_num; i++){
            assert(fabs(finite_difference_force_tmp_c[nw][l][m][a][k][i] - forces_tmp_c[nw][l][m][a][k][i]) < 1.e-8);
          }
        }
      }
    }
  }
}


printf("OK\n");

printf("Forces dtmp_c\n");

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

rc = qmckl_set_nucleus_coord(context, 'T', &(nucl_coord[0]), 3*nucl_num);
assert(rc == QMCKL_SUCCESS);

double forces_dtmp_c[walk_num][cord_num][cord_num+1][4][nucl_num][4][elec_num];
rc = qmckl_get_forces_dtmp_c(context, &forces_dtmp_c[0][0][0][0][0][0][0], 4*4*nucl_num*walk_num*elec_num*(cord_num+1)*cord_num);
assert(rc == QMCKL_SUCCESS);


double finite_difference_force_dtmp_c[walk_num][cord_num][cord_num+1][3][nucl_num][4][elec_num];

nucleus_coord = (double*) malloc(3 * nucl_num * sizeof(double));
if (nucleus_coord == NULL) {
  return QMCKL_ALLOCATION_FAILED;
}

rc = qmckl_get_nucleus_coord(context, 'N', nucleus_coord, 3*nucl_num);

temp_coord = (double*) malloc(3 * nucl_num * sizeof(double));
if (temp_coord == NULL) {
  free(nucleus_coord);
  return QMCKL_ALLOCATION_FAILED;
}
double doutput[walk_num][cord_num][cord_num+1][nucl_num][4][elec_num];

// Copy original coordinates
for (int i = 0; i < 3 * nucl_num; i++) {
  temp_coord[i] = nucleus_coord[i];
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
      rc = qmckl_get_jastrow_champ_dtmp_c(context,
              &doutput[0][0][0][0][0][0],
              4*4*nucl_num*walk_num*elec_num*(cord_num+1)*cord_num);
      assert(rc == QMCKL_SUCCESS);

      // Accumulate derivative using finite-difference coefficients
      for (int nw=0 ; nw<walk_num ; nw++) {
        for (int l = 0; l < cord_num; l++) {
          for (int mm = 0; mm <= cord_num; mm++) {
            for (int i = 0; i < elec_num; i++) {
              for (int j = 0; j < 4; j++) {
                if (m == -4) {
                  finite_difference_force_dtmp_c[nw][l][mm][k][a][j][i] = 0.0;
                }
                finite_difference_force_dtmp_c[nw][l][mm][k][a][j][i] += coef[m + 4] * doutput[nw][l][mm][a][j][i]/delta_x;
              }

            }
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


for (int nw = 0; nw < walk_num; nw++){
  for (int l = 0; l < cord_num; l++){
    for (int m = 0; m <= cord_num; m++){
      for (int a = 0; a < nucl_num; a++) {
        for (int k = 0; k < 3; k++){
          for (int i = 0; i < elec_num; i++){
            for (int kk = 0; kk<4; kk++){
            //printf("nw=%i l=%i m=%i k=%i a=%i kk=%i i=%i\n",nw,l,m,k,a,kk,i);
            //printf("%.10f\t", finite_difference_force_dtmp_c[nw][l][m][k][a][kk][i]);
            //printf("%.10f\n", forces_dtmp_c[nw][l][m][k][a][kk][i]);
            //assert(fabs(finite_difference_force_dtmp_c[nw][l][m][k][a][kk][i] - forces_dtmp_c[nw][l][m][k][a][kk][i]) < 1.e-8);
            }

          }
        }
      }
    }
  }
}
for (int nw = 0; nw < walk_num; nw++){
  for (int l = 0; l < cord_num; l++){
    for (int m = 0; m <= cord_num; m++){
      for (int a = 0; a < nucl_num; a++) {
        for (int k = 0; k < 3; k++){
          for (int i = 0; i < elec_num; i++){
            for (int kk = 0; kk<4; kk++){
              assert(fabs(finite_difference_force_dtmp_c[nw][l][m][k][a][kk][i] - forces_dtmp_c[nw][l][m][k][a][kk][i]) < 1.e-8);
            }
          }
        }
      }
    }
  }
}


printf("OK\n");

printf("Forces Jastrow een\n");

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

rc = qmckl_set_nucleus_coord(context, 'T', &(nucl_coord[0]), 3*nucl_num);
assert(rc == QMCKL_SUCCESS);

double forces_jastrow_een[walk_num][nucl_num][3];
rc = qmckl_get_forces_jastrow_een(context, &forces_jastrow_een[0][0][0], 3*nucl_num*walk_num);
assert(rc == QMCKL_SUCCESS);

double finite_difference_force_een[walk_num][nucl_num][3];
rc = qmckl_finite_difference_deriv_n(context, delta_x, &qmckl_get_jastrow_champ_factor_een, &(finite_difference_force_een[0][0][0]), 1);

for (int nw = 0; nw < walk_num; nw++){
  for (int a = 0; a < nucl_num; a++) {
    for (int k = 0; k < 3; k++){
      //printf("%.10f\t", finite_difference_force_een[nw][a][k]);
      //printf("%.10f\n", forces_jastrow_een[nw][a][k]);
    }
  }
}
for (int nw = 0; nw < walk_num; nw++){
  for (int a = 0; a < nucl_num; a++) {
    for (int k = 0; k < 3; k++){
      assert(fabs(finite_difference_force_een[nw][a][k] - forces_jastrow_een[nw][a][k]) < 1.e-8);
    }
  }
}
printf("OK\n");

printf("Forces of een_rescaled_n_gl\n");

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

rc = qmckl_set_nucleus_coord(context, 'T', &(nucl_coord[0]), 3*nucl_num);
assert(rc == QMCKL_SUCCESS);

double forces_een_n[walk_num][cord_num+1][3][nucl_num][4][elec_num];
rc = qmckl_get_forces_een_rescaled_n_gl(context, &forces_een_n[0][0][0][0][0][0], 3*4*nucl_num*walk_num*elec_num*(cord_num+1));
assert(rc == QMCKL_SUCCESS);


double finite_difference_force_een_n[walk_num][cord_num+1][3][nucl_num][4][elec_num];

nucleus_coord = (double*) malloc(3 * nucl_num * sizeof(double));
if (nucleus_coord == NULL) {
  return QMCKL_ALLOCATION_FAILED;
}

rc = qmckl_get_nucleus_coord(context, 'N', nucleus_coord, 3*nucl_num);

temp_coord = (double*) malloc(3 * nucl_num * sizeof(double));
if (temp_coord == NULL) {
  free(nucleus_coord);
  return QMCKL_ALLOCATION_FAILED;
}
double ddoutput[walk_num][cord_num+1][nucl_num][4][elec_num];

// Copy original coordinates
for (int i = 0; i < 3 * nucl_num; i++) {
  temp_coord[i] = nucleus_coord[i];
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
      rc = qmckl_get_jastrow_champ_een_rescaled_n_gl(context,&ddoutput[0][0][0][0][0], 4*elec_num*nucl_num*(cord_num+1)*walk_num);
      assert(rc == QMCKL_SUCCESS);

      // Accumulate derivative using finite-difference coefficients
      for (int nw=0 ; nw<walk_num ; nw++) {
        for (int l = 0; l <= cord_num; l++) {
            for (int i = 0; i < elec_num; i++) {
              for (int j = 0; j < 4; j++) {
                if (m == -4) {
                  finite_difference_force_een_n[nw][l][k][a][j][i] = 0.0;
                }
                finite_difference_force_een_n[nw][l][k][a][j][i] += coef[m + 4] * ddoutput[nw][l][a][j][i]/delta_x;
              }


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


for (int nw = 0; nw < walk_num; nw++){
  for (int l = 0; l <=cord_num; l++){
      for (int a = 0; a < nucl_num; a++) {
        for (int k = 0; k < 3; k++){
          for (int i = 0; i < elec_num; i++){
            for (int kk = 0; kk<4; kk++){
            //printf("nw=%i l=%i k=%i a=%i kk=%i i=%i\n",nw,l,k,a,kk,i);
            //printf("%.10f\t", finite_difference_force_een_n[nw][l][k][a][kk][i]);
            //printf("%.10f\n", forces_een_n[nw][l][k][a][kk][i]);
            //assert(fabs(finite_difference_force_een_n[nw][l][k][a][kk][i] - forces_een_n[nw][l][k][a][kk][i]) < 1.e-6);
            }

          }
        }
    }
  }
}
for (int nw = 0; nw < walk_num; nw++){
  for (int l = 0; l <= cord_num; l++){
      for (int a = 0; a < nucl_num; a++) {
        for (int k = 0; k < 3; k++){
          for (int i = 0; i < elec_num; i++){
            for (int kk = 0; kk<4; kk++){
              assert(fabs(finite_difference_force_een_n[nw][l][k][a][kk][i] - forces_een_n[nw][l][k][a][kk][i]) < 1.e-6);
            }
          }
      }
    }
  }
}


printf("OK\n");

printf("Forces Jastrow een G\n");

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

rc = qmckl_set_nucleus_coord(context, 'T', &(nucl_coord[0]), 3*nucl_num);
assert(rc == QMCKL_SUCCESS);

double forces_jastrow_een_g[walk_num][3][nucl_num][3][elec_num];
rc = qmckl_get_forces_jastrow_een_g(context, &forces_jastrow_een_g[0][0][0][0][0], 3*3*nucl_num*walk_num*elec_num);
assert(rc == QMCKL_SUCCESS);

double finite_difference_force_een_g[walk_num][nucl_num][3][4][elec_num];
rc = qmckl_finite_difference_deriv_n(context, delta_x, &qmckl_get_jastrow_champ_factor_een_gl, &finite_difference_force_een_g[0][0][0][0][0], 4*elec_num);

for (int nw = 0; nw < walk_num; nw++){
  for (int a = 0; a < nucl_num; a++) {
    for (int k = 0; k < 3; k++){
      for (int i = 0; i < elec_num; i++){
        for (int l = 0; l < 3; l++){
          //printf("finite_difference_force_een_g[%i][%i][%i][%i][%i] %+3.10f \n", nw,a,k,l,i,finite_difference_force_een_g[nw][a][k][l][i]);
          //printf("forces_jastrow_een_g         [%i][%i][%i][%i][%i] %+3.10f\n", nw,k,a,l,i,forces_jastrow_een_g[nw][k][a][l][i]);
          //assert(fabs(finite_difference_force_een_g[nw][a][k][l][i] - forces_jastrow_een_g[nw][k][a][l][i]) < 1.e-8);
        }
      }
    }
  }
}
for (int nw = 0; nw < walk_num; nw++){
  for (int a = 0; a < nucl_num; a++) {
    for (int k = 0; k < 3; k++){
      for (int i = 0; i < elec_num; i++){
        for (int l = 1; l < 3; l++){
          assert(fabs(finite_difference_force_een_g[nw][a][k][l][i] - forces_jastrow_een_g[nw][k][a][l][i]) < 1.e-8);
        }
      }
    }
  }
}
printf("OK\n");

printf("Forces Jastrow een l\n");

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

rc = qmckl_set_nucleus_coord(context, 'T', &(nucl_coord[0]), 3*nucl_num);
assert(rc == QMCKL_SUCCESS);

double forces_jastrow_een_l[walk_num][3][nucl_num];
rc = qmckl_get_forces_jastrow_een_l(context, &forces_jastrow_een_l[0][0][0], 3*nucl_num*walk_num);
assert(rc == QMCKL_SUCCESS);

double finite_difference_force_een_l[walk_num][nucl_num][3][4][elec_num];
rc = qmckl_finite_difference_deriv_n(context, delta_x, &qmckl_get_jastrow_champ_factor_een_gl, &finite_difference_force_een_l[0][0][0][0][0], 4*elec_num);

double finite_difference_force_een_l_sum[walk_num][nucl_num][3];
for (int nw = 0; nw < walk_num; nw++){
  for (int a = 0; a < nucl_num; a++) {
    for (int k = 0; k < 3; k++){
      finite_difference_force_een_l_sum[nw][a][k] = 0;
      for (int i = 0; i < elec_num; i++){
        finite_difference_force_een_l_sum[nw][a][k] = finite_difference_force_een_l_sum[nw][a][k] + finite_difference_force_een_l[nw][a][k][3][i];
      }
    }
  }
}

for (int nw = 0; nw < walk_num; nw++){
  for (int a = 0; a < nucl_num; a++) {
    for (int k = 0; k < 3; k++){
      //printf("finite_difference_force_een_l_sum[%i][%i][%i] %+3.10f \n", nw,a,k,finite_difference_force_een_l_sum[nw][a][k]);
      //printf("forces_jastrow_een_l         [%i][%i][%i] %+3.10f\n", nw,k,a,forces_jastrow_een_l[nw][k][a]);
      assert(fabs(finite_difference_force_een_l_sum[nw][a][k] - forces_jastrow_een_l[nw][k][a]) < 1.e-8);
    }
  }
}
for (int nw = 0; nw < walk_num; nw++){
  for (int a = 0; a < nucl_num; a++) {
    for (int k = 0; k < 3; k++){
      assert(fabs(finite_difference_force_een_l_sum[nw][a][k] - forces_jastrow_een_l[nw][k][a]) < 1.e-8);
    }
  }
}
printf("OK\n");

printf("Forces Jastrow single een\n");

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

rc = qmckl_set_nucleus_coord(context, 'T', &(nucl_coord[0]), 3*nucl_num);
assert(rc == QMCKL_SUCCESS);

rc = qmckl_set_single_point(context, 'N', 2, new_coords, 3*walk_num);
assert (rc == QMCKL_SUCCESS);

double forces_jastrow_single_een[walk_num][nucl_num][3];
rc = qmckl_get_forces_jastrow_single_een(context, &forces_jastrow_single_een[0][0][0], 3*nucl_num*walk_num);
assert(rc == QMCKL_SUCCESS);

double finite_difference_force_single_een[walk_num][nucl_num][3];
rc = qmckl_finite_difference_deriv_n(context, delta_x, &qmckl_get_jastrow_champ_single_een, &(finite_difference_force_single_een[0][0][0]), 1);

for (int nw = 0; nw < walk_num; nw++){
  for (int a = 0; a < nucl_num; a++) {
    for (int k = 0; k < 3; k++){
      //printf("%.10f\t", finite_difference_force_single_een[nw][a][k]);
      //printf("%.10f\n", forces_jastrow_single_een[nw][a][k]);
      assert(fabs(finite_difference_force_single_een[nw][a][k] - forces_jastrow_single_een[nw][a][k]) < 1.e-8);
    }
  }
}
printf("OK\n");

rc = qmckl_context_destroy(context);
  assert (rc == QMCKL_SUCCESS);
   context = qmckl_context_create();

   nucl_num      = chbrclf_nucl_num;
   nucl_charge   = chbrclf_charge;
   nucl_coord    = &(chbrclf_nucl_coord[0][0]);

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

#define walk_num 1 // chbrclf_walk_num
#define elec_num chbrclf_elec_num
#define prim_num chbrclf_prim_num

  elec_up_num   = chbrclf_elec_up_num;
  elec_dn_num   = chbrclf_elec_dn_num;
  elec_coord    = &(chbrclf_elec_coord[0][0][0]);

  rc = qmckl_set_electron_num (context, elec_up_num, elec_dn_num);
  assert (rc == QMCKL_SUCCESS);

  assert(qmckl_electron_provided(context));

  int64_t point_num = elec_num;

  rc = qmckl_set_point(context, 'N', point_num, elec_coord, point_num*3);
  assert(rc == QMCKL_SUCCESS);

  int64_t mo_num = chbrclf_mo_num;
  rc = qmckl_set_mo_basis_mo_num(context, mo_num);
  assert (rc == QMCKL_SUCCESS);

  const double  * mo_coefficient          =  &(chbrclf_mo_coef[0]);

  rc = qmckl_set_mo_basis_coefficient(context, mo_coefficient, chbrclf_mo_num*chbrclf_ao_num);
  assert (rc == QMCKL_SUCCESS);

printf("Forces AO value\n");

rc = qmckl_set_nucleus_coord(context, 'T', &(nucl_coord[0]), 3*nucl_num);
assert(rc == QMCKL_SUCCESS);

double * forces_ao_value = (double*) malloc(nucl_num * 3 *point_num * ao_num *sizeof(double));
rc = qmckl_get_forces_ao_value(context, &forces_ao_value[0], 3*nucl_num*ao_num*point_num);
assert(rc == QMCKL_SUCCESS);

double * finite_difference_force_ao_value = (double*) malloc(3 *nucl_num * point_num * ao_num *sizeof(double));

nucleus_coord = (double*) malloc(3 * nucl_num * sizeof(double));
if (nucleus_coord == NULL) {
  return QMCKL_ALLOCATION_FAILED;
}

rc = qmckl_get_nucleus_coord(context, 'N', nucleus_coord, 3*nucl_num);

temp_coord = (double*) malloc(3 * nucl_num * sizeof(double));
if (temp_coord == NULL) {
  free(nucleus_coord);
  return QMCKL_ALLOCATION_FAILED;
}
double * ao_output = (double*) malloc(point_num*ao_num*sizeof(double));

// Copy original coordinates
for (int i = 0; i < 3 * nucl_num; i++) {
  temp_coord[i] = nucleus_coord[i];
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
      rc = qmckl_get_ao_basis_ao_value(context,&ao_output[0], point_num*ao_num);
      assert(rc == QMCKL_SUCCESS);

      // Accumulate derivative using finite-difference coefficients
      for (int i = 0; i < point_num; i++) {
        for (int j = 0; j < ao_num; j++) {
          if (m == -4) {
            finite_difference_force_ao_value[k*ao_num*point_num*nucl_num + a*ao_num*point_num + i*ao_num + j] = 0.0;
          }
          finite_difference_force_ao_value[k*ao_num*point_num*nucl_num + a*ao_num*point_num + i*ao_num + j] += coef[m + 4] * ao_output[i*ao_num + j]/delta_x;
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
free(ao_output);


for (int a = 0; a < nucl_num; a++) {
  for (int i = 0; i < point_num; i++){
    for (int k = 0; k < 3; k++){
      for (int j = 0; j < ao_num; j++){
//        printf("k=%i a=%i i=%i j=%i\n", k, a, i, j);
//        printf("%.10e\t", finite_difference_force_ao_value[k*ao_num*point_num*nucl_num + a*ao_num*point_num + i*ao_num + j]);
//        printf("%.10e\n", forces_ao_value[a*3*ao_num*point_num + i*ao_num*3 + k*ao_num + j]);
        assert(fabs(finite_difference_force_ao_value[k*ao_num*point_num*nucl_num + a*ao_num*point_num + i*ao_num + j] -
                    forces_ao_value[a*3*ao_num*point_num + i*ao_num*3 + k*ao_num + j]) < 1.e-9);
        }
    }
  }
}

free(forces_ao_value);
free(finite_difference_force_ao_value);

printf("OK\n");

printf("Forces MO value\n");

rc = qmckl_set_nucleus_coord(context, 'T', &(nucl_coord[0]), 3*nucl_num);
assert(rc == QMCKL_SUCCESS);

double * forces_mo_value = (double*) malloc(nucl_num * point_num* 3 * mo_num *sizeof(double));
//double forces_mo_value[nucl_num][3][point_num][mo_num];
rc = qmckl_get_forces_mo_value(context, &forces_mo_value[0], 3*nucl_num*mo_num*point_num);
assert(rc == QMCKL_SUCCESS);

//double finite_difference_force_mo_value[3][nucl_num][point_num][mo_num];
double * finite_difference_force_mo_value = (double*) malloc(nucl_num* 3 * point_num * mo_num * sizeof(double));

nucleus_coord = (double*) malloc(3 * nucl_num * sizeof(double));
if (nucleus_coord == NULL) {
  return QMCKL_ALLOCATION_FAILED;
}

rc = qmckl_get_nucleus_coord(context, 'N', nucleus_coord, 3*nucl_num);

temp_coord = (double*) malloc(3 * nucl_num * sizeof(double));
if (temp_coord == NULL) {
  free(nucleus_coord);
  return QMCKL_ALLOCATION_FAILED;
}

double* mo_output = (double*) malloc(point_num * mo_num * sizeof(double));
if (mo_output == NULL) {
  free(temp_coord);
  free(nucleus_coord);
  return QMCKL_ALLOCATION_FAILED;
}


// Copy original coordinates
for (int i = 0; i < 3 * nucl_num; i++) {
  temp_coord[i] = nucleus_coord[i];
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
      rc = qmckl_get_mo_basis_mo_value(context,&mo_output[0], point_num*mo_num);
      assert(rc == QMCKL_SUCCESS);

      // Accumulate derivative using finite-difference coefficients
      for (int i = 0; i < point_num; i++) {
        for (int j = 0; j < mo_num; j++) {
          if (m == -4) {
            finite_difference_force_mo_value[k*mo_num*point_num*nucl_num + a*mo_num*point_num + i*mo_num + j]= 0.0;
          }
          finite_difference_force_mo_value[k*mo_num*point_num*nucl_num + a*mo_num*point_num + i*mo_num + j] += coef[m + 4] * mo_output[i*mo_num+ j]/delta_x;
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
free(mo_output);


for (int a = 0; a < nucl_num; a++) {
  for (int i = 0; i < point_num; i++){
    for (int k = 0; k < 3; k++){
      for (int j = 0; j < mo_num; j++){
        //printf("k=%i a=%i i=%i j=%i\n", k, a, i, j);
        //printf("%.10f\t", finite_difference_force_mo_value[k*mo_num*point_num*nucl_num + a*mo_num*point_num + i*mo_num + j]);
        //printf("%.10f\n", forces_mo_value[a*3*mo_num*point_num + k*mo_num*point_num + i*mo_num + j]);
        assert(fabs(finite_difference_force_mo_value[k*mo_num*point_num*nucl_num + a*mo_num*point_num + i*mo_num + j] -
                    forces_mo_value[a*3*mo_num*point_num + i*mo_num*3 + k*mo_num + j]) < 1.e-9);
      }
    }
  }
}

free(forces_mo_value);
free(finite_difference_force_mo_value);

printf("OK\n");

printf("Forces MO gradient\n");

qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
assert (ctx != NULL);

rc = qmckl_set_nucleus_coord(context, 'T', &(nucl_coord[0]), 3*nucl_num);
assert(rc == QMCKL_SUCCESS);


double * forces_mo_g_doc = (double*) malloc(nucl_num* 3 * point_num* 3 * mo_num *sizeof(double));
rc = qmckl_get_forces_mo_g(context, &forces_mo_g_doc[0], 3*3*nucl_num*mo_num*point_num);

rc = qmckl_compute_forces_mo_g_doc(context,
                                   ctx->ao_basis.ao_num,
                                   ctx->mo_basis.mo_num,
                                   ctx->point.num,
                                   ctx->nucleus.num,
                                   ctx->ao_basis.shell_num,
                                   ctx->ao_basis.nucleus_index,
                                   ctx->ao_basis.nucleus_shell_num,
                                   ctx->ao_basis.shell_ang_mom,
                                   ctx->mo_basis.coefficient_t,
                                   ctx->ao_basis.ao_hessian,
                                   &forces_mo_g_doc[0]);

assert(rc == QMCKL_SUCCESS);
double * forces_mo_g = (double*) malloc(nucl_num* 3 * point_num* 3 * mo_num *sizeof(double));

rc = qmckl_compute_forces_mo_g_hpc(context,
                                   ctx->ao_basis.ao_num,
                                   ctx->mo_basis.mo_num,
                                   ctx->point.num,
                                   ctx->nucleus.num,
                                   ctx->ao_basis.shell_num,
                                   ctx->ao_basis.nucleus_index,
                                   ctx->ao_basis.nucleus_shell_num,
                                   ctx->ao_basis.shell_ang_mom,
                                   ctx->mo_basis.coefficient_t,
                                   ctx->ao_basis.ao_hessian,
                                   &forces_mo_g[0]);

assert(rc == QMCKL_SUCCESS);

double * finite_difference_force_mo_g = (double*) malloc(nucl_num* 3 * point_num* 3 * mo_num * sizeof(double));

nucleus_coord = (double*) malloc(3 * nucl_num * sizeof(double));
if (nucleus_coord == NULL) {
  return QMCKL_ALLOCATION_FAILED;
}

rc = qmckl_get_nucleus_coord(context, 'N', nucleus_coord, 3*nucl_num);

temp_coord = (double*) malloc(3 * nucl_num * sizeof(double));
if (temp_coord == NULL) {
  free(nucleus_coord);
  return QMCKL_ALLOCATION_FAILED;
}

mo_output = (double*) malloc(5 * point_num * mo_num * sizeof(double));
if (mo_output == NULL) {
  free(temp_coord);
  free(nucleus_coord);
  return QMCKL_ALLOCATION_FAILED;
}


// Copy original coordinates
for (int i = 0; i < 3 * nucl_num; i++) {
  temp_coord[i] = nucleus_coord[i];
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
      rc = qmckl_get_mo_basis_mo_vgl(context,&mo_output[0], 5*point_num*mo_num);
      assert(rc == QMCKL_SUCCESS);

      // Accumulate derivative using finite-difference coefficients
      for (int i = 0; i < point_num; i++) {
        for (int n = 0; n < 3; n++){
          for (int j = 0; j < mo_num; j++) {
            if (m == -4) {
              finite_difference_force_mo_g[k*3*mo_num*point_num*nucl_num + a*3*mo_num*point_num + i*3*mo_num + n*mo_num + j] = 0.0;
            }
            finite_difference_force_mo_g[k*3*mo_num*point_num*nucl_num + a*3*mo_num*point_num + i*3*mo_num + n*mo_num + j] +=
              coef[m + 4] * mo_output[i*mo_num*5 + (n+1)*mo_num + j]/delta_x;
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
free(mo_output);


for (int a = 0; a < nucl_num; a++) {
  for (int n = 0; n < 3; n++){
    for (int j = 0; j < point_num; j++){
      for (int m = 0; m < 3; m++){
        for (int i = 0; i < mo_num; i++){
          //printf("a=%i n=%i j=%i m=%i i=%i\n", a, n, j, m, i);
          //printf("hpc %.10f\n", forces_mo_g[a*9*mo_num*point_num + n*3*mo_num*point_num + j*3*mo_num + m*mo_num + i]);
          //printf("doc %.10f\n", forces_mo_g_doc[a*9*mo_num*point_num + n*3*mo_num*point_num + j*3*mo_num + m*mo_num + i]);
          //printf("fd  %.10f\n", finite_difference_force_mo_g[n*3*mo_num*point_num*nucl_num + a*3*mo_num*point_num + j*3*mo_num + m*mo_num + i]);
          //assert(fabs(forces_mo_g[a*9*mo_num*point_num + n*3*mo_num*point_num + j*3*mo_num + m*mo_num + i] -
          //            forces_mo_g_doc[a*9*mo_num*point_num + n*3*mo_num*point_num + j*3*mo_num + m*mo_num + i]
          //            ) < 1.e-9);
          assert(fabs(finite_difference_force_mo_g[n*3*mo_num*point_num*nucl_num +
                                                   a*3*mo_num*point_num + j*3*mo_num +
                                                   m*mo_num + i] - 
                      forces_mo_g[a*9*mo_num*point_num + n*3*mo_num*point_num +
                                  j*3*mo_num + m*mo_num + i]) < 1.e-9);
        }
      }
    }
  }
}

free(forces_mo_g_doc);
free(forces_mo_g);
free(finite_difference_force_mo_g);

printf("OK\n");

printf("Forces MO laplacian\n");

rc = qmckl_set_nucleus_coord(context, 'T', &(nucl_coord[0]), 3*nucl_num);
assert(rc == QMCKL_SUCCESS);

double * forces_mo_l = (double*) malloc(nucl_num * point_num* 3 * mo_num *sizeof(double));
rc = qmckl_get_forces_mo_l(context, &forces_mo_l[0], 3*nucl_num*mo_num*point_num);
assert(rc == QMCKL_SUCCESS);

double * finite_difference_force_mo_l= (double*) malloc(nucl_num* 3 * point_num * mo_num * sizeof(double));

nucleus_coord = (double*) malloc(3 * nucl_num * sizeof(double));
if (nucleus_coord == NULL) {
  return QMCKL_ALLOCATION_FAILED;
}

rc = qmckl_get_nucleus_coord(context, 'N', nucleus_coord, 3*nucl_num);

temp_coord = (double*) malloc(3 * nucl_num * sizeof(double));
if (temp_coord == NULL) {
  free(nucleus_coord);
  return QMCKL_ALLOCATION_FAILED;
}

mo_output = (double*) malloc(5 * point_num * mo_num * sizeof(double));
if (mo_output == NULL) {
  free(temp_coord);
  free(nucleus_coord);
  return QMCKL_ALLOCATION_FAILED;
}


// Copy original coordinates
for (int i = 0; i < 3 * nucl_num; i++) {
  temp_coord[i] = nucleus_coord[i];
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
      rc = qmckl_get_mo_basis_mo_vgl(context,&mo_output[0], 5*point_num*mo_num);
      assert(rc == QMCKL_SUCCESS);

      // Accumulate derivative using finite-difference coefficients
      for (int i = 0; i < point_num; i++) {
        for (int j = 0; j < mo_num; j++) {
          if (m == -4) {
            finite_difference_force_mo_l[k*mo_num*point_num*nucl_num + a*mo_num*point_num + i*mo_num + j] = 0.0;
          }
          finite_difference_force_mo_l[k*mo_num*point_num*nucl_num + a*mo_num*point_num + i*mo_num + j] += coef[m + 4] * mo_output[i*mo_num*5 + 4*mo_num + j]/delta_x;
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
free(mo_output);



for (int j = 0; j < mo_num; j++){
  for (int i = 0; i < point_num; i++){
    for (int a = 0; a < nucl_num; a++) {
      for (int k = 0; k < 3; k++){
        //printf("k=%i a=%i i=%i j=%i\n", k, a, i, j);
        //printf("%.10f\t", finite_difference_force_mo_l[k*mo_num*point_num*nucl_num + a*mo_num*point_num + i*mo_num + j]);
        //printf("%.10f\n", forces_mo_l[a*3*mo_num*point_num + k*mo_num*point_num + i*mo_num + j]);
        assert(fabs(finite_difference_force_mo_l[k*mo_num*point_num*nucl_num + a*mo_num*point_num + i*mo_num + j] - forces_mo_l[a*3*mo_num*point_num + k*mo_num*point_num + i*mo_num + j]) < 1.e-8);
      }
    }
    }
}

free(forces_mo_l);
free(finite_difference_force_mo_l);

printf("OK\n");

rc = qmckl_context_destroy(context);
    assert (rc == QMCKL_SUCCESS);

    return 0;
}
