#include "qmckl.h"
#include <string.h>
#include <assert.h>
#include <math.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdbool.h>
#include <stdio.h>
#include "n2.h"
#include "qmckl_context_private_type.h"
#include "qmckl_jastrow_champ_private_func.h"

#ifdef __GNUC__
#pragma GCC push_options
#pragma GCC optimize ("O0")
#endif

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
int64_t size_max;

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

printf("asymp_jasb\n");
double asymp_jasb[2];
rc = qmckl_check(context,
                 qmckl_get_jastrow_champ_asymp_jasb(context, &(asymp_jasb[0]),2)
                 );

// calculate asymp_jasb
assert(fabs(asymp_jasb[0]-0.7115733522582638) < 1.e-12);
assert(fabs(asymp_jasb[1]-1.043287918508297 ) < 1.e-12);

printf("asymp_jasb_hpc\n");
double asymp_jasb_doc[2];
double asymp_jasb_hpc[2];
// calculate asymp_jasb
rc = qmckl_check(context,
                 qmckl_compute_jastrow_champ_asymp_jasb_doc (context,
                                            bord_num,
                                            b_vector,
                                            rescale_factor_ee,
                                            0,
                                            &(asymp_jasb_doc[0]) )
                 );
rc = qmckl_check(context,
                 qmckl_compute_jastrow_champ_asymp_jasb_hpc (context,
                                            bord_num,
                                            b_vector,
                                            rescale_factor_ee,
                                            0,
                                            &(asymp_jasb_hpc[0]) )
                 );
assert(fabs(asymp_jasb_doc[0]-asymp_jasb_hpc[0]) < 1.e-8);
assert(fabs(asymp_jasb_doc[1]-asymp_jasb_hpc[1]) < 1.e-8);

assert(qmckl_electron_provided(context));

{
  printf("ee_distance_rescaled\n");
  double ee_distance_rescaled[walk_num * elec_num * elec_num];
  rc = qmckl_get_jastrow_champ_ee_distance_rescaled(context,
                                                    ee_distance_rescaled,
                                                    walk_num*elec_num*elec_num);

  // (e1,e2,w)
  // (0,0,0) == 0.
  assert(ee_distance_rescaled[0] == 0.);

  // (1,0,0) == (0,1,0)
  assert(ee_distance_rescaled[1] == ee_distance_rescaled[elec_num]);

  // value of (1,0,0)
  assert(fabs(ee_distance_rescaled[1]-0.6347507420688708) < 1.e-12);

  // (0,0,1) == 0.
  assert(ee_distance_rescaled[5*elec_num + 5] == 0.);

  // (1,0,1) == (0,1,1)
  assert(ee_distance_rescaled[5*elec_num+6] == ee_distance_rescaled[6*elec_num+5]);

  // value of (1,0,1)
  assert(fabs(ee_distance_rescaled[5*elec_num+6]-0.3941735387855409) < 1.e-12);

  printf("ee_distance_rescaled_hpc\n");
  double ee_distance[walk_num * elec_num * elec_num];
  rc = qmckl_get_electron_ee_distance(context, &(ee_distance[0]), walk_num*elec_num*elec_num);
  assert(rc == QMCKL_SUCCESS);

  double ee_distance_rescaled_doc[walk_num * elec_num * elec_num * (cord_num+1)];
  memset(ee_distance_rescaled_doc, 0, sizeof(ee_distance_rescaled_doc));

  rc = qmckl_compute_een_rescaled_e_doc (context, walk_num,
                                         elec_num, cord_num,
                                         rescale_factor_ee,
                                         &(ee_distance[0]),
                                         &(ee_distance_rescaled_doc[0]));

  assert(rc == QMCKL_SUCCESS);

  double ee_distance_rescaled_hpc[walk_num * elec_num * elec_num * (cord_num+1)];
  memset(ee_distance_rescaled_hpc, 0, sizeof(ee_distance_rescaled_hpc));

  rc = qmckl_compute_een_rescaled_e_hpc (context, walk_num,
                                         elec_num, cord_num,
                                         rescale_factor_ee,
                                         &(ee_distance[0]),
                                         &(ee_distance_rescaled_hpc[0]));
  assert(rc == QMCKL_SUCCESS);

  for (int64_t i=0 ; i<walk_num*elec_num*elec_num*(cord_num+1) ; i++) {
    if (fabs(ee_distance_rescaled_hpc[i] - ee_distance_rescaled_doc[i]) > 1.e-10) {
      printf("i = %ld\n", i);
      printf("ee_distance_rescaled_hpc = %f\n", ee_distance_rescaled_hpc[i]);
      printf("ee_distance_rescaled_doc = %f\n", ee_distance_rescaled_doc[i]);
      fflush(stdout);
    }
    assert(fabs(ee_distance_rescaled_hpc[i] - ee_distance_rescaled_doc[i]) < 1.e-10);
  }
}

assert(qmckl_electron_provided(context));

{
  printf("ee_distance_rescaled_gl\n");
  double fd[walk_num][elec_num][elec_num][4];

  double delta_x = 0.001;

  // Finite difference coefficients for gradients
  double coef[9] = { 1.0/280.0, -4.0/105.0,  1.0/5.0, -4.0/5.0,         0.0, 4.0/5.0, -1.0/5.0, 4.0/105.0, -1.0/280.0 };

  // Finite difference coefficients for Laplacian
  double coef2[9]= {-1.0/560.0,	 8.0/315.0, -1.0/5.0,  8.0/5.0, -205.0/72.0, 8.0/5.0, -1.0/5.0,	8.0/315.0, -1.0/560.0 };

  qmckl_exit_code rc;

  int64_t elec_num;
  rc = qmckl_get_electron_num(context, &elec_num);
  if (rc != QMCKL_SUCCESS) {
    return rc;
  }

  double elec_coord[walk_num][elec_num][3];
  rc = qmckl_get_electron_coord (context, 'N', &(elec_coord[0][0][0]), 3*walk_num*elec_num);

  double temp_coord[walk_num][elec_num][3];
  memcpy(&(temp_coord[0][0][0]), &(elec_coord[0][0][0]), sizeof(temp_coord));

  double function_values[walk_num][elec_num][elec_num];

  memset(&(fd[0][0][0][0]), 0, sizeof(fd));

  for (int64_t i = 0; i < elec_num; i++) {
    for (int64_t k = 0; k < 3; k++) {
      for (int64_t m = -4; m <= 4; m++) {          // Apply finite difference displacement

        for (int64_t nw=0 ; nw<walk_num ; nw++) {
          temp_coord[nw][i][k] = elec_coord[nw][i][k] + (double) m * delta_x;
        }

        // Update coordinates in the context
        rc = qmckl_set_electron_coord (context, 'N', walk_num,
                                       &(temp_coord[0][0][0]),
                                       walk_num*3*elec_num);
        assert(rc == QMCKL_SUCCESS);

        // Call the provided function
        rc = qmckl_get_jastrow_champ_ee_distance_rescaled(context,
                                                          &(function_values[0][0][0]),
                                                          elec_num*elec_num*walk_num);
        assert(rc == QMCKL_SUCCESS);

        // Accumulate derivative using finite-difference coefficients
        for (int64_t nw=0 ; nw<walk_num ; nw++) {
          for (int64_t j = 0; j < elec_num; j++) {
            fd[nw][j][i][k] += coef [m + 4] * function_values[nw][j][i];
            fd[nw][j][i][3] += coef2[m + 4] * function_values[nw][j][i];
          }
        }
      }
      for (int64_t nw=0 ; nw<walk_num ; nw++) {
        temp_coord[nw][i][k] = elec_coord[nw][i][k];
      }
    }
  }

  // Reset coordinates in the context
  rc = qmckl_set_electron_coord (context, 'N', walk_num,
                                 &(elec_coord[0][0][0]),
                                 walk_num*3*elec_num);
  assert(rc == QMCKL_SUCCESS);

  // Normalize by the step size
  for (int64_t nw=0 ; nw<walk_num ; nw++) {
    for (int64_t i = 0; i < elec_num; i++) {
      for (int64_t k = 0; k < 4; k++) {
        for (int64_t j = 0; j < elec_num; j++) {
          fd[nw][i][j][k] /= delta_x;
        }
      }
      for (int64_t j = 0; j < elec_num; j++) {
        fd[nw][i][j][3] /= delta_x;
      }
    }
  }


  double ee_distance_rescaled_gl[walk_num][elec_num][elec_num][4];
  rc = qmckl_check(context,
                   qmckl_get_jastrow_champ_ee_distance_rescaled_gl(context,
                                                                   &(ee_distance_rescaled_gl[0][0][0][0]),
                                                                   walk_num*elec_num*4*elec_num)
                   );

  assert(rc == QMCKL_SUCCESS);

  for (int nw = 0; nw < walk_num; nw++){
    for (int i = 0; i < elec_num; i++) {
      for (int j = 0; j < elec_num; j++) {
        for (int k = 0; k < 3; k++){
          if (fabs(fd[nw][i][j][k] - ee_distance_rescaled_gl[nw][i][j][k]) > 1.e-12) {
            printf("nw=%d i=%d j=%d k=%d\n", nw, i, j, k);
            printf("fd                     =%f\n", fd[nw][i][j][k]);
            printf("ee_distance_rescaled_gl=%f\n", ee_distance_rescaled_gl[nw][i][j][k]);
            fflush(stdout);
          }
          assert(fabs(fd[nw][i][j][k] - ee_distance_rescaled_gl[nw][i][j][k]) < 1.e-8);
        }
        int k=3;
        if (i != j) {
          if (fabs(fd[nw][i][j][k] - ee_distance_rescaled_gl[nw][i][j][k]) > 1.e-12) {
            printf("nw=%d i=%d j=%d k=%d\n", nw, i, j, k);
            printf("fd                     =%f\n", fd[nw][i][j][k]);
            printf("ee_distance_rescaled_gl=%f\n", ee_distance_rescaled_gl[nw][i][j][k]);
            fflush(stdout);
          }
          assert(fabs(fd[nw][i][j][k] - ee_distance_rescaled_gl[nw][i][j][k]) < 1.e-6);
        }
      }
    }
  }
  printf("OK\n");

  printf("ee_distance_rescaled_gl_hpc\n");

  double ee_distance_rescaled_gl_doc[walk_num*elec_num*elec_num*4];
  rc = qmckl_compute_ee_distance_rescaled_gl_doc (context,
                                                  elec_num,
                                                  rescale_factor_ee,
                                                  walk_num,
                                                  &(elec_coord[0][0][0]),
                                                  &(ee_distance_rescaled_gl_doc[0]));
  assert(rc == QMCKL_SUCCESS);

  double ee_distance_rescaled_gl_hpc[walk_num*elec_num*elec_num*4];
  rc = qmckl_compute_ee_distance_rescaled_gl_hpc (context,
                                                  elec_num,
                                                  rescale_factor_ee,
                                                  walk_num,
                                                  &(elec_coord[0][0][0]),
                                                  &(ee_distance_rescaled_gl_hpc[0]));
  assert(rc == QMCKL_SUCCESS);

  for (int64_t i = 0; i < walk_num*nucl_num*elec_num*4; i++) {
    if (fabs(ee_distance_rescaled_gl_hpc[i] - ee_distance_rescaled_gl_doc[i]) > 1.e-12) {
      printf("i=%ld\n", i);
      printf("ee_distance_rescaled_gl_doc=%f\n", ee_distance_rescaled_gl_doc[i]);
      printf("ee_distance_rescaled_gl_hpc=%f\n", ee_distance_rescaled_gl_hpc[i]);
      fflush(stdout);
    }
    assert(fabs(ee_distance_rescaled_gl_doc[i] - ee_distance_rescaled_gl_hpc[i]) < 1.e-8);
  }
}

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

{
  printf("factor_ee\n");
  double factor_ee[walk_num];
  rc = qmckl_check(context,
                   qmckl_get_jastrow_champ_factor_ee(context, &(factor_ee[0]), walk_num)
                   );

  // calculate factor_ee
  printf("%20.15f\n%20.15f\n",factor_ee[0],-16.83886184243964);
  fflush(stdout);
  assert(fabs(factor_ee[0]+16.83886184243964) < 1.e-12);

  printf("factor_ee_hpc\n");
  double ee_distance_rescaled[walk_num*elec_num*elec_num];
  rc = qmckl_get_jastrow_champ_ee_distance_rescaled(context,
                                                    &(ee_distance_rescaled[0]),
                                                    walk_num*elec_num*elec_num);
  assert(rc == QMCKL_SUCCESS);

  int64_t up_num;
  rc = qmckl_get_electron_up_num(context, &up_num);
  assert(rc == QMCKL_SUCCESS);

  double factor_ee_doc[walk_num];
  rc = qmckl_compute_jastrow_champ_factor_ee_doc(context,
                                                 walk_num,
                                                 elec_num,
                                                 up_num,
                                                 bord_num,
                                                 b_vector,
                                                 &(ee_distance_rescaled[0]),
                                                 &(asymp_jasb[0]),
                                                 0,
                                                 &(factor_ee_doc[0]));
  assert (rc == QMCKL_SUCCESS);

  double factor_ee_hpc[walk_num];
  rc = qmckl_compute_jastrow_champ_factor_ee_hpc(context,
                                                 walk_num,
                                                 elec_num,
                                                 up_num,
                                                 bord_num,
                                                 b_vector,
                                                 &(ee_distance_rescaled[0]),
                                                 &(asymp_jasb[0]),
                                                 0,
                                                 &(factor_ee_hpc[0]));
  assert (rc == QMCKL_SUCCESS);

  for (int64_t i = 0; i < walk_num; i++) {
    if (fabs(factor_ee_doc[i] - factor_ee_hpc[i]) > 1.e-12) {
      printf("i=%ld\n", i);
      printf("factor_ee_doc=%f\n", factor_ee_doc[i]);
      printf("factor_ee_hpc=%f\n", factor_ee_hpc[i]);
      fflush(stdout);
    }
    assert(fabs(factor_ee_doc[i] - factor_ee_hpc[i]) < 1.e-8);
  }
}

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

//----
{
  printf("factor_ee_gl\n");
  double fd[walk_num][4][elec_num];
  double delta_x = 0.001;

  // Finite difference coefficients for gradients
  double coef[9] = { 1.0/280.0, -4.0/105.0,  1.0/5.0, -4.0/5.0,         0.0, 4.0/5.0, -1.0/5.0, 4.0/105.0, -1.0/280.0 };

  // Finite difference coefficients for Laplacian
  double coef2[9]= {-1.0/560.0,	 8.0/315.0, -1.0/5.0,  8.0/5.0, -205.0/72.0, 8.0/5.0, -1.0/5.0,	8.0/315.0, -1.0/560.0 };

  qmckl_exit_code rc;

  int64_t walk_num;
  rc = qmckl_get_electron_walk_num(context, &walk_num);
  if (rc != QMCKL_SUCCESS) {
    return rc;
  }

  int64_t elec_num;
  rc = qmckl_get_electron_num(context, &elec_num);
  if (rc != QMCKL_SUCCESS) {
    return rc;
  }

  double elec_coord[walk_num][elec_num][3];
  rc = qmckl_get_electron_coord (context, 'N', &(elec_coord[0][0][0]), 3*walk_num*elec_num);

  double temp_coord[walk_num][elec_num][3];
  memcpy(&(temp_coord[0][0][0]), &(elec_coord[0][0][0]), sizeof(temp_coord));

  double function_values[walk_num];

  memset(&(fd[0][0][0]), 0, sizeof(fd));

  for (int64_t i = 0; i < elec_num; i++) {
    for (int64_t k = 0; k < 3; k++) {
      for (int64_t m = -4; m <= 4; m++) {          // Apply finite difference displacement

        for (int64_t nw=0 ; nw<walk_num ; nw++) {
          temp_coord[nw][i][k] = elec_coord[nw][i][k] + (double) m * delta_x;
        }

        // Update coordinates in the context
        rc = qmckl_set_electron_coord (context, 'N', walk_num, &(temp_coord[0][0][0]), walk_num*3*elec_num);
        assert(rc == QMCKL_SUCCESS);

        // Call the provided function
        rc = qmckl_get_jastrow_champ_factor_ee(context, &(function_values[0]), walk_num);
        assert(rc == QMCKL_SUCCESS);

        // Accumulate derivative using finite-difference coefficients
        for (int64_t nw=0 ; nw<walk_num ; nw++) {
          fd[nw][k][i] += coef [m + 4] * function_values[nw];
          fd[nw][3][i] += coef2[m + 4] * function_values[nw];
        }
      }
      for (int64_t nw=0 ; nw<walk_num ; nw++) {
        temp_coord[nw][i][k] = elec_coord[nw][i][k];
      }
    }
  }

  // Reset coordinates in the context
  rc = qmckl_set_electron_coord (context, 'N', walk_num, &(elec_coord[0][0][0]), walk_num*3*elec_num);
  assert(rc == QMCKL_SUCCESS);

  // Normalize by the step size
  for (int64_t nw=0 ; nw<walk_num ; nw++) {
    for (int64_t k = 0; k < 4; k++) {
      for (int64_t i = 0; i < elec_num; i++) {
          fd[nw][k][i] /= delta_x;
      }
    }
    for (int64_t i = 0; i < elec_num; i++) {
      fd[nw][3][i] /= delta_x;
    }
  }


  double factor_ee_gl[walk_num][4][elec_num];
  rc = qmckl_check(context,
                   qmckl_get_jastrow_champ_factor_ee_gl(context,
                                                        &(factor_ee_gl[0][0][0]),
                                                        walk_num*4*elec_num)
                   );

  assert(rc == QMCKL_SUCCESS);

  for (int nw = 0; nw < walk_num; nw++){
    for (int i = 0; i < elec_num; i++) {
      for (int k = 0; k < 3; k++){
        if (fabs(fd[nw][k][i] - factor_ee_gl[nw][k][i]) > 1.e-12) {
          printf("nw=%d i=%d k=%d\n", nw, i, k);
          printf("fd=%f factor_ee_gl=%f\n", fd[nw][k][i], factor_ee_gl[nw][k][i]);
          fflush(stdout);
        }
        assert(fabs(fd[nw][k][i] - factor_ee_gl[nw][k][i]) < 1.e-8);
      }
      int k=3;
      if (fabs(fd[nw][k][i] - factor_ee_gl[nw][k][i]) > 1.e-12) {
        printf("nw=%d i=%d k=%d\n", nw, i, k);
        printf("fd=%f factor_ee_gl=%f\n", fd[nw][k][i], factor_ee_gl[nw][k][i]);
        fflush(stdout);
      }
      assert(fabs(fd[nw][k][i] - factor_ee_gl[nw][k][i]) < 1.e-5);
    }
  }
  printf("OK\n");
}
{
  printf("factor_ee_gl_hpc\n");
  int64_t up_num;
  rc = qmckl_get_electron_up_num(context, &up_num);
  assert(rc == QMCKL_SUCCESS);

  double ee_distance_rescaled[walk_num*elec_num*elec_num];
  rc = qmckl_get_jastrow_champ_ee_distance_rescaled(context,
                                                    &(ee_distance_rescaled[0]),
                                                    walk_num*elec_num*elec_num);
  assert(rc == QMCKL_SUCCESS);

  double ee_distance_rescaled_gl[4*walk_num*elec_num*elec_num];
  rc = qmckl_get_jastrow_champ_ee_distance_rescaled_gl(context,
                                                      &(ee_distance_rescaled_gl[0]),
                                                      4*walk_num*elec_num*elec_num);
  assert(rc == QMCKL_SUCCESS);

  double factor_ee_gl_doc[walk_num*4*elec_num];
  memset(&(factor_ee_gl_doc[0]), 0, sizeof(factor_ee_gl_doc));

  rc = qmckl_compute_jastrow_champ_factor_ee_gl_doc(context,
                                                    walk_num,
                                                    elec_num,
                                                    up_num,
                                                    bord_num,
                                                    b_vector,
                                                    &(ee_distance_rescaled[0]),
                                                    &(ee_distance_rescaled_gl[0]),
                                                    0,
                                                    &(factor_ee_gl_doc[0]));
  assert(rc == QMCKL_SUCCESS);

  double factor_ee_gl_hpc[walk_num*4*elec_num];
  memset(&(factor_ee_gl_hpc[0]), 0, sizeof(factor_ee_gl_hpc));

  rc = qmckl_compute_jastrow_champ_factor_ee_gl_hpc(context,
                                                    walk_num,
                                                    elec_num,
                                                    up_num,
                                                    bord_num,
                                                    b_vector,
                                                    &(ee_distance_rescaled[0]),
                                                    &(ee_distance_rescaled_gl[0]),
                                                    0,
                                                    &(factor_ee_gl_hpc[0]));
  assert(rc == QMCKL_SUCCESS);

  for (int64_t i = 0 ; i < walk_num*4*elec_num ; i++) {
    if (fabs(factor_ee_gl_hpc[i] - factor_ee_gl_doc[i]) > 1.e-12) {
      printf("i=%ld\nfactor_ee_gl_hpc=%f\nfactor_ee_gl_doc=%f\n", i, factor_ee_gl_hpc[i], factor_ee_gl_doc[i]);
      fflush(stdout);
    }
    assert(fabs(factor_ee_gl_hpc[i] - factor_ee_gl_doc[i]) < 1.e-8);
  }
}

{
assert(qmckl_electron_provided(context));
assert(qmckl_nucleus_provided(context));

double en_distance_rescaled[walk_num][nucl_num][elec_num];

rc = qmckl_check(context,
                 qmckl_get_jastrow_champ_en_distance_rescaled(context,
                                                              &(en_distance_rescaled[0][0][0]),
                                                              walk_num*nucl_num*elec_num)
                 );
assert (rc == QMCKL_SUCCESS);

// (e,n,w) in Fortran notation
// (1,1,1)
assert(fabs(en_distance_rescaled[0][0][0] - 0.4942158656729477) < 1.e-12);
// (1,2,1)
assert(fabs(en_distance_rescaled[0][1][0] - 1.2464137498005765) < 1.e-12);
// (2,1,1)
assert(fabs(en_distance_rescaled[0][0][1] - 0.5248654474756858) < 1.e-12);
// (1,1,2)
assert(fabs(en_distance_rescaled[0][0][5] - 0.19529459944794733) < 1.e-12);
// (1,2,2)
assert(fabs(en_distance_rescaled[0][1][5] - 1.2091967687767369) < 1.e-12);
// (2,1,2)
assert(fabs(en_distance_rescaled[0][0][6] - 0.4726452953409436) < 1.e-12);

}

{
  printf("en_distance_rescaled_hpc\n");

  double en_distance_rescaled_doc[walk_num*nucl_num*elec_num];
  memset(&(en_distance_rescaled_doc[0]), 0, walk_num*nucl_num*elec_num*sizeof(double));
  rc = qmckl_compute_en_distance_rescaled_doc(context, elec_num, nucl_num, type_nucl_num,
                                              type_nucl_vector, rescale_factor_en, walk_num,
                                              elec_coord, nucl_coord, en_distance_rescaled_doc);
  assert(rc == QMCKL_SUCCESS);

  double en_distance_rescaled_hpc[walk_num*nucl_num*elec_num];
  memset(&(en_distance_rescaled_hpc[0]), 0, walk_num*nucl_num*elec_num*sizeof(double));
  rc = qmckl_compute_en_distance_rescaled_hpc(context, elec_num, nucl_num, type_nucl_num,
                                              type_nucl_vector, rescale_factor_en, walk_num,
                                              elec_coord, nucl_coord, en_distance_rescaled_hpc);
  assert(rc == QMCKL_SUCCESS);

  for (int64_t i=0 ; i<walk_num*nucl_num*elec_num ; ++i) {
    if (fabs(en_distance_rescaled_doc[i] - en_distance_rescaled_hpc[i]) > 1.e-12) {
      printf("i = %ld, doc = %e, hpc = %e\n", i, en_distance_rescaled_doc[i], en_distance_rescaled_hpc[i]);
      fflush(stdout);
    }
    assert(fabs(en_distance_rescaled_doc[i] - en_distance_rescaled_hpc[i]) < 1.e-8);
  }
}

assert(qmckl_electron_provided(context));

{
  printf("en_distance_rescaled_gl\n");
  double fd[walk_num][nucl_num][elec_num][4];

  double delta_x = 0.001;

  // Finite difference coefficients for gradients
  double coef[9] = { 1.0/280.0, -4.0/105.0,  1.0/5.0, -4.0/5.0,         0.0, 4.0/5.0, -1.0/5.0, 4.0/105.0, -1.0/280.0 };

  // Finite difference coefficients for Laplacian
  double coef2[9]= {-1.0/560.0,	 8.0/315.0, -1.0/5.0,  8.0/5.0, -205.0/72.0, 8.0/5.0, -1.0/5.0,	8.0/315.0, -1.0/560.0 };

  qmckl_exit_code rc;

  int64_t walk_num;
  rc = qmckl_get_electron_walk_num(context, &walk_num);
  if (rc != QMCKL_SUCCESS) {
    return rc;
  }

  int64_t elec_num;
  rc = qmckl_get_electron_num(context, &elec_num);
  if (rc != QMCKL_SUCCESS) {
    return rc;
  }

  double elec_coord[walk_num][elec_num][3];
  rc = qmckl_get_electron_coord (context, 'N', &(elec_coord[0][0][0]), 3*walk_num*elec_num);

  double temp_coord[walk_num][elec_num][3];
  memcpy(&(temp_coord[0][0][0]), &(elec_coord[0][0][0]), sizeof(temp_coord));

  double function_values[walk_num][nucl_num][elec_num];

  memset(&(fd[0][0][0][0]), 0, sizeof(fd));

  for (int64_t i = 0; i < elec_num; i++) {
    for (int64_t k = 0; k < 3; k++) {
      for (int64_t m = -4; m <= 4; m++) {          // Apply finite difference displacement

        for (int64_t nw=0 ; nw<walk_num ; nw++) {
          temp_coord[nw][i][k] = elec_coord[nw][i][k] + (double) m * delta_x;
        }

        // Update coordinates in the context
        rc = qmckl_set_electron_coord (context, 'N', walk_num, &(temp_coord[0][0][0]), walk_num*3*elec_num);
        assert(rc == QMCKL_SUCCESS);

        // Call the provided function
        rc = qmckl_get_jastrow_champ_en_distance_rescaled(context, &(function_values[0][0][0]), nucl_num*elec_num*walk_num);
        assert(rc == QMCKL_SUCCESS);

        // Accumulate derivative using finite-difference coefficients
        for (int64_t nw=0 ; nw<walk_num ; nw++) {
          for (int64_t j = 0; j < nucl_num; j++) {
            fd[nw][j][i][k] += coef [m + 4] * function_values[nw][j][i];
            fd[nw][j][i][3] += coef2[m + 4] * function_values[nw][j][i];
          }
        }
      }
      for (int64_t nw=0 ; nw<walk_num ; nw++) {
        temp_coord[nw][i][k] = elec_coord[nw][i][k];
      }
    }
  }

  // Reset coordinates in the context
  rc = qmckl_set_electron_coord (context, 'N', walk_num, &(elec_coord[0][0][0]), walk_num*3*elec_num);
  assert(rc == QMCKL_SUCCESS);

  // Normalize by the step size
  for (int64_t nw=0 ; nw<walk_num ; nw++) {
    for (int64_t i = 0; i < nucl_num; i++) {
      for (int64_t k = 0; k < 4; k++) {
        for (int64_t j = 0; j < elec_num; j++) {
          fd[nw][i][j][k] /= delta_x;
        }
      }
      for (int64_t j = 0; j < elec_num; j++) {
        fd[nw][i][j][3] /= delta_x;
      }
    }
  }


  double en_distance_rescaled_gl[walk_num][nucl_num][elec_num][4];
  rc = qmckl_check(context,
                   qmckl_get_jastrow_champ_en_distance_rescaled_gl(context,
                                                                   &(en_distance_rescaled_gl[0][0][0][0]),
                                                                   walk_num*nucl_num*4*elec_num)
                   );

  assert(rc == QMCKL_SUCCESS);

  for (int nw = 0; nw < walk_num; nw++){
    for (int i = 0; i < nucl_num; i++) {
      for (int j = 0; j < elec_num; j++) {
        for (int k = 0; k < 3; k++){
          printf("%.10f\t", fd[nw][i][j][k]);
          printf("%.10f\n", en_distance_rescaled_gl[nw][i][j][k]);
          fflush(stdout);
          assert(fabs(fd[nw][i][j][k] - en_distance_rescaled_gl[nw][i][j][k]) < 1.e-8);
        }
        int k=3;
        if (i != j) {
          printf("%.10f\t", fd[nw][i][j][k]);
          printf("%.10f\n", en_distance_rescaled_gl[nw][i][j][k]);
          fflush(stdout);
          assert(fabs(fd[nw][i][j][k] - en_distance_rescaled_gl[nw][i][j][k]) < 1.e-6);
        }
      }
    }
  }
  printf("OK\n");
}

{
  printf("en_distance_rescaled_gl_hpc\n");

  double en_distance_rescaled_gl_doc[walk_num*nucl_num*elec_num*4];
  rc = qmckl_compute_en_distance_rescaled_gl_doc (context,
          elec_num, nucl_num, type_nucl_num, type_nucl_vector, rescale_factor_en,
          walk_num, elec_coord, nucl_coord,
          &(en_distance_rescaled_gl_doc[0]));
  assert(rc == QMCKL_SUCCESS);

  double en_distance_rescaled_gl_hpc[walk_num*nucl_num*elec_num*4];
  rc = qmckl_compute_en_distance_rescaled_gl_hpc (context,
          elec_num, nucl_num, type_nucl_num, type_nucl_vector, rescale_factor_en,
          walk_num, elec_coord, nucl_coord,
          &(en_distance_rescaled_gl_hpc[0]));
  assert(rc == QMCKL_SUCCESS);

  for (int i = 0; i < walk_num*nucl_num*elec_num*4; i++) {
    if (fabs(en_distance_rescaled_gl_doc[i] - en_distance_rescaled_gl_hpc[i]) > 1.e-8) {
      printf("i = %d, doc = %e, hpc = %e\n", i, en_distance_rescaled_gl_doc[i], en_distance_rescaled_gl_hpc[i]);
      fflush(stdout);
    }
    assert(fabs(en_distance_rescaled_gl_doc[i] - en_distance_rescaled_gl_hpc[i]) < 1.e-8);
  }
}

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

double factor_en[walk_num];
rc = qmckl_get_jastrow_champ_factor_en(context, factor_en,walk_num);

// calculate factor_en
printf("%f %f\n", factor_en[0], 22.781375792083587);
fflush(stdout);
assert(fabs(22.781375792083587 - factor_en[0]) < 1.e-12);

{
  printf("factor_en_hpc\n");
  double asymp_jasa[type_nucl_num];
  rc = qmckl_get_jastrow_champ_asymp_jasa(context, asymp_jasa, type_nucl_num);
  assert(rc == QMCKL_SUCCESS);

  double en_distance_rescaled[walk_num*nucl_num*elec_num];
  rc = qmckl_get_jastrow_champ_en_distance_rescaled(context,
                                                    en_distance_rescaled,
                                                    walk_num*nucl_num*elec_num);
  assert(rc == QMCKL_SUCCESS);

  double factor_en_doc[walk_num];
  memset(&(factor_en_doc[0]), 0, sizeof(factor_en_doc));
  rc = qmckl_compute_jastrow_champ_factor_en_doc (context,
          walk_num, elec_num, nucl_num, type_nucl_num, type_nucl_vector,
          aord_num, a_vector,
          en_distance_rescaled, asymp_jasa, factor_en_doc);
  assert(rc == QMCKL_SUCCESS);

  double factor_en_hpc[walk_num];
  rc = qmckl_compute_jastrow_champ_factor_en_hpc (context,
          walk_num, elec_num, nucl_num, type_nucl_num, type_nucl_vector,
          aord_num, a_vector,
          en_distance_rescaled, asymp_jasa, factor_en_hpc);
  assert(rc == QMCKL_SUCCESS);

  for (int64_t i = 0; i < walk_num; i++) {
    assert(fabs(factor_en_doc[i] - factor_en_hpc[i]) < 1.e-10);
  }
}

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

{
  printf("factor_en_gl\n");
  double fd[walk_num][4][elec_num];
  double delta_x = 0.001;

  // Finite difference coefficients for gradients
  double coef[9] = { 1.0/280.0, -4.0/105.0,  1.0/5.0, -4.0/5.0,         0.0, 4.0/5.0, -1.0/5.0, 4.0/105.0, -1.0/280.0 };

  // Finite difference coefficients for Laplacian
  double coef2[9]= {-1.0/560.0,	 8.0/315.0, -1.0/5.0,  8.0/5.0, -205.0/72.0, 8.0/5.0, -1.0/5.0,	8.0/315.0, -1.0/560.0 };

  qmckl_exit_code rc;

  int64_t walk_num;
  rc = qmckl_get_electron_walk_num(context, &walk_num);
  if (rc != QMCKL_SUCCESS) {
    return rc;
  }

  int64_t elec_num;
  rc = qmckl_get_electron_num(context, &elec_num);
  if (rc != QMCKL_SUCCESS) {
    return rc;
  }

  double elec_coord[walk_num][elec_num][3];
  rc = qmckl_get_electron_coord (context, 'N', &(elec_coord[0][0][0]), 3*walk_num*elec_num);

  double temp_coord[walk_num][elec_num][3];
  memcpy(&(temp_coord[0][0][0]), &(elec_coord[0][0][0]), sizeof(temp_coord));

  double function_values[walk_num];

  memset(&(fd[0][0][0]), 0, sizeof(fd));

  for (int64_t i = 0; i < elec_num; i++) {
    for (int64_t k = 0; k < 3; k++) {
      for (int64_t m = -4; m <= 4; m++) {          // Apply finite difference displacement

        for (int64_t nw=0 ; nw<walk_num ; nw++) {
          temp_coord[nw][i][k] = elec_coord[nw][i][k] + (double) m * delta_x;
        }

        // Update coordinates in the context
        rc = qmckl_set_electron_coord (context, 'N', walk_num, &(temp_coord[0][0][0]), walk_num*3*elec_num);
        assert(rc == QMCKL_SUCCESS);

        // Call the provided function
        rc = qmckl_get_jastrow_champ_factor_en(context, &(function_values[0]), walk_num);
        assert(rc == QMCKL_SUCCESS);

        // Accumulate derivative using finite-difference coefficients
        for (int64_t nw=0 ; nw<walk_num ; nw++) {
          fd[nw][k][i] += coef [m + 4] * function_values[nw];
          fd[nw][3][i] += coef2[m + 4] * function_values[nw];
        }
      }
      for (int64_t nw=0 ; nw<walk_num ; nw++) {
        temp_coord[nw][i][k] = elec_coord[nw][i][k];
      }
    }
  }

  // Reset coordinates in the context
  rc = qmckl_set_electron_coord (context, 'N', walk_num, &(elec_coord[0][0][0]), walk_num*3*elec_num);
  assert(rc == QMCKL_SUCCESS);

  // Normalize by the step size
  for (int64_t nw=0 ; nw<walk_num ; nw++) {
    for (int64_t k = 0; k < 4; k++) {
      for (int64_t i = 0; i < elec_num; i++) {
          fd[nw][k][i] /= delta_x;
      }
    }
    for (int64_t i = 0; i < elec_num; i++) {
      fd[nw][3][i] /= delta_x;
    }
  }


  double factor_en_gl[walk_num][4][elec_num];
  rc = qmckl_check(context,
                   qmckl_get_jastrow_champ_factor_en_gl(context,
                                                        &(factor_en_gl[0][0][0]),
                                                        walk_num*4*elec_num)
                   );

  assert(rc == QMCKL_SUCCESS);

  for (int nw = 0; nw < walk_num; nw++){
    for (int i = 0; i < elec_num; i++) {
      for (int k = 0; k < 3; k++){
        printf("%.10f\t", fd[nw][k][i]);
        printf("%.10f\n", factor_en_gl[nw][k][i]);
        fflush(stdout);
        assert(fabs(fd[nw][k][i] - factor_en_gl[nw][k][i]) < 1.e-8);
      }
      int k=3;
      if (fabs(fd[nw][k][i] - factor_en_gl[nw][k][i]) > 1.e-5) {
        printf("i=%d doc=%f hpc=%f\n", i, fd[nw][k][i], factor_en_gl[nw][k][i]);
        fflush(stdout);
      }
      assert(fabs(fd[nw][k][i] - factor_en_gl[nw][k][i]) < 1.e-5);
    }
  }
  printf("OK\n");
}

{
  printf("factor_en_gl_hpc\n");

  double en_distance_rescaled[walk_num][nucl_num][elec_num];
  rc = qmckl_get_jastrow_champ_en_distance_rescaled(context,
                                                    &(en_distance_rescaled[0][0][0]),
                                                    walk_num*nucl_num*elec_num);
  assert(rc == QMCKL_SUCCESS);

  double en_distance_rescaled_gl[walk_num][4][elec_num][nucl_num];
  rc = qmckl_get_jastrow_champ_en_distance_rescaled_gl(context,
                                                       &(en_distance_rescaled_gl[0][0][0][0]),
                                                       walk_num*4*elec_num*nucl_num);
  assert(rc == QMCKL_SUCCESS);

  double factor_en_gl_doc[walk_num*4*elec_num];
  memset(&(factor_en_gl_doc[0]), 0, sizeof(factor_en_gl_doc));
  rc = qmckl_compute_jastrow_champ_factor_en_gl_doc(context, walk_num, elec_num,
       nucl_num, type_nucl_num, type_nucl_vector, aord_num, &(a_vector[0]),
       &(en_distance_rescaled[0][0][0]), &(en_distance_rescaled_gl[0][0][0][0]), &(factor_en_gl_doc[0]));
  assert(rc == QMCKL_SUCCESS);

  double factor_en_gl_hpc[walk_num*4*elec_num];
  memset(&(factor_en_gl_hpc[0]), 0, sizeof(factor_en_gl_hpc));
  rc = qmckl_compute_jastrow_champ_factor_en_gl_hpc(context, walk_num, elec_num,
       nucl_num, type_nucl_num, type_nucl_vector, aord_num, &(a_vector[0]),
       &(en_distance_rescaled[0][0][0]), &(en_distance_rescaled_gl[0][0][0][0]), &(factor_en_gl_hpc[0]));
  assert(rc == QMCKL_SUCCESS);

  for (int64_t i = 0; i < walk_num*4*elec_num; i++) {
    if (fabs(factor_en_gl_doc[i] - factor_en_gl_hpc[i]) > 1.e-12) {
      printf("i=%ld doc=%f hpc=%f\n", i, factor_en_gl_doc[i], factor_en_gl_hpc[i]);
      fflush(stdout);
    }
    assert(fabs(factor_en_gl_doc[i] - factor_en_gl_hpc[i]) < 1.e-8);
  }
}

assert(qmckl_electron_provided(context));
{

  double een_rescaled_e[walk_num][(cord_num + 1)][elec_num][elec_num];
  rc = qmckl_get_jastrow_champ_een_rescaled_e(context, &(een_rescaled_e[0][0][0][0]),elec_num*elec_num*(cord_num+1)*walk_num);

  assert( fabs(een_rescaled_e[0][0][0][0] - (0.0)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][0][0][1] - (1.0)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][0][0][2] - (1.0)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][0][1][0] - (1.0)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][0][1][1] - (0.0)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][0][1][2] - (1.0)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][0][2][0] - (1.0)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][0][2][1] - (1.0)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][0][2][2] - (0.0)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][1][0][0] - (0.0)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][1][0][1] - (0.6191495547586775)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][1][0][2] - (0.2211015082992776)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][1][1][0] - (0.6191495547586775)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][1][1][1] - (0.0)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][1][1][2] - (0.1891080158548495)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][1][2][0] - (0.2211015082992776)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][1][2][1] - (0.1891080158548495)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][1][2][2] - (0.0)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][2][0][0] - (0.0)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][2][0][1] - (0.38334617115786856)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][2][0][2] - (0.048885876972215525)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][2][1][0] - (0.38334617115786856)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][2][1][1] - (0.0)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][2][1][2] - (0.03576184166055801)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][2][2][0] - (0.048885876972215525)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][2][2][1] - (0.03576184166055801)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][2][2][2] - (0.0)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][3][0][0] - (0.0)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][3][0][1] - (0.23734861119083814)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][3][0][2] - (0.010808741133089782)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][3][1][0] - (0.23734861119083814)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][3][1][1] - (0.0)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][3][1][2] - (0.006762850919743422)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][3][2][0] - (0.010808741133089782)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][3][2][1] - (0.006762850919743422)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][3][2][2] - (0.0)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][4][0][0] - (0.0)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][4][0][1] - (0.14695428694139787)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][4][0][2] - (0.002389828967342592)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][4][1][0] - (0.14695428694139787)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][4][1][1] - (0.0)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][4][1][2] - (0.0012789093189548221)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][4][2][0] - (0.002389828967342592)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][4][2][1] - (0.0012789093189548221)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][4][2][2] - (0.0)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][5][0][0] - (0.0)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][5][0][1] - (0.09098668132964542)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][5][0][2] - (0.0005283947892567522)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][5][1][0] - (0.09098668132964542)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][5][1][1] - (0.0)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][5][1][2] - (0.0002418520037658232)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][5][2][0] - (0.0005283947892567522)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][5][2][1] - (0.0002418520037658232)) < 1.e-10 );
  assert( fabs(een_rescaled_e[0][5][2][2] - (0.0)) < 1.e-10 );
}

{
  printf("een_rescaled_e_hpc\n");
  double ee_distance[walk_num * elec_num * elec_num];
  rc = qmckl_get_electron_ee_distance(context, &(ee_distance[0]), walk_num*elec_num*elec_num);
  assert(rc == QMCKL_SUCCESS);

  double een_rescaled_e_doc[walk_num][cord_num+1][elec_num][elec_num];
  memset(&(een_rescaled_e_doc[0][0][0][0]), 0, sizeof(een_rescaled_e_doc));
  rc = qmckl_compute_een_rescaled_e(context, walk_num, elec_num, cord_num,
                                    rescale_factor_ee, &(ee_distance[0]), &(een_rescaled_e_doc[0][0][0][0]));
  assert(rc == QMCKL_SUCCESS);

  double een_rescaled_e_hpc[walk_num][cord_num+1][elec_num][elec_num];
  memset(&(een_rescaled_e_hpc[0][0][0][0]), 0, sizeof(een_rescaled_e_hpc));
  rc = qmckl_compute_een_rescaled_e_hpc(context, walk_num, elec_num, cord_num,
                                        rescale_factor_ee, &(ee_distance[0]), &(een_rescaled_e_hpc[0][0][0][0]));
  assert(rc == QMCKL_SUCCESS);

  for (int64_t i = 0; i < walk_num; i++) {
    for (int64_t j = 0; j < cord_num+1; j++) {
      for (int64_t k = 0; k < elec_num; k++) {
        for (int64_t l = 0; l < elec_num; l++) {
          if (fabs(een_rescaled_e_doc[i][j][k][l] - een_rescaled_e_hpc[i][j][k][l]) > 1.e-12) {
            printf("i=%ld j=%ld k=%ld l=%ld doc=%f hpc=%f\n", i, j, k, l, een_rescaled_e_doc[i][j][k][l], een_rescaled_e_hpc[i][j][k][l]);
            fflush(stdout);
          }
          assert(fabs(een_rescaled_e_doc[i][j][k][l] - een_rescaled_e_hpc[i][j][k][l]) < 1.e-8);
        }
      }
    }
  }
}

assert(qmckl_electron_provided(context));

{
  double een_rescaled_e_gl[walk_num][(cord_num + 1)][elec_num][4][elec_num];
  size_max=walk_num*(cord_num + 1)*elec_num*4*elec_num;
  rc = qmckl_get_jastrow_champ_een_rescaled_e_gl(context,
            &(een_rescaled_e_gl[0][0][0][0][0]),size_max);
  printf( "%e %e\n", een_rescaled_e_gl[0][0][0][0][0], 0.0);
  assert( fabs(een_rescaled_e_gl[0][0][0][0][0] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][0][0][0][1], 0.0);
  assert( fabs(een_rescaled_e_gl[0][0][0][0][1] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][0][0][0][2], -0.0);
  assert( fabs(een_rescaled_e_gl[0][0][0][0][2] - (-0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][0][0][1][0], 0.0);
  assert( fabs(een_rescaled_e_gl[0][0][0][1][0] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][0][0][1][1], 0.0);
  assert( fabs(een_rescaled_e_gl[0][0][0][1][1] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][0][0][1][2], 0.0);
  assert( fabs(een_rescaled_e_gl[0][0][0][1][2] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][0][0][2][0], 0.0);
  assert( fabs(een_rescaled_e_gl[0][0][0][2][0] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][0][0][2][1], -0.0);
  assert( fabs(een_rescaled_e_gl[0][0][0][2][1] - (-0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][0][0][2][2], 0.0);
  assert( fabs(een_rescaled_e_gl[0][0][0][2][2] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][0][0][3][0], 0.0);
  assert( fabs(een_rescaled_e_gl[0][0][0][3][0] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][0][0][3][1], -0.0);
  assert( fabs(een_rescaled_e_gl[0][0][0][3][1] - (-0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][0][0][3][2], -0.0);
  assert( fabs(een_rescaled_e_gl[0][0][0][3][2] - (-0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][0][1][0][0], -0.0);
  assert( fabs(een_rescaled_e_gl[0][0][1][0][0] - (-0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][0][1][0][1], 0.0);
  assert( fabs(een_rescaled_e_gl[0][0][1][0][1] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][0][1][0][2], -0.0);
  assert( fabs(een_rescaled_e_gl[0][0][1][0][2] - (-0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][0][1][1][0], -0.0);
  assert( fabs(een_rescaled_e_gl[0][0][1][1][0] - (-0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][0][1][1][1], 0.0);
  assert( fabs(een_rescaled_e_gl[0][0][1][1][1] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][0][1][1][2], 0.0);
  assert( fabs(een_rescaled_e_gl[0][0][1][1][2] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][0][1][2][0], 0.0);
  assert( fabs(een_rescaled_e_gl[0][0][1][2][0] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][0][1][2][1], 0.0);
  assert( fabs(een_rescaled_e_gl[0][0][1][2][1] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][0][1][2][2], 0.0);
  assert( fabs(een_rescaled_e_gl[0][0][1][2][2] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][0][1][3][0], -0.0);
  assert( fabs(een_rescaled_e_gl[0][0][1][3][0] - (-0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][0][1][3][1], 0.0);
  assert( fabs(een_rescaled_e_gl[0][0][1][3][1] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][0][1][3][2], -0.0);
  assert( fabs(een_rescaled_e_gl[0][0][1][3][2] - (-0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][0][2][0][0], 0.0);
  assert( fabs(een_rescaled_e_gl[0][0][2][0][0] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][0][2][0][1], 0.0);
  assert( fabs(een_rescaled_e_gl[0][0][2][0][1] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][0][2][0][2], 0.0);
  assert( fabs(een_rescaled_e_gl[0][0][2][0][2] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][0][2][1][0], -0.0);
  assert( fabs(een_rescaled_e_gl[0][0][2][1][0] - (-0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][0][2][1][1], -0.0);
  assert( fabs(een_rescaled_e_gl[0][0][2][1][1] - (-0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][0][2][1][2], 0.0);
  assert( fabs(een_rescaled_e_gl[0][0][2][1][2] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][0][2][2][0], -0.0);
  assert( fabs(een_rescaled_e_gl[0][0][2][2][0] - (-0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][0][2][2][1], -0.0);
  assert( fabs(een_rescaled_e_gl[0][0][2][2][1] - (-0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][0][2][2][2], 0.0);
  assert( fabs(een_rescaled_e_gl[0][0][2][2][2] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][0][2][3][0], -0.0);
  assert( fabs(een_rescaled_e_gl[0][0][2][3][0] - (-0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][0][2][3][1], -0.0);
  assert( fabs(een_rescaled_e_gl[0][0][2][3][1] - (-0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][0][2][3][2], 0.0);
  assert( fabs(een_rescaled_e_gl[0][0][2][3][2] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][1][0][0][0], 0.0);
  assert( fabs(een_rescaled_e_gl[0][1][0][0][0] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][1][0][0][1], 0.15675618873333058);
  assert( fabs(een_rescaled_e_gl[0][1][0][0][1] - (0.15675618873333058)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][1][0][0][2], -0.09831391870751387);
  assert( fabs(een_rescaled_e_gl[0][1][0][0][2] - (-0.09831391870751387)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][1][0][1][0], 0.0);
  assert( fabs(een_rescaled_e_gl[0][1][0][1][0] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][1][0][1][1], 0.2937567144823047);
  assert( fabs(een_rescaled_e_gl[0][1][0][1][1] - (0.2937567144823047)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][1][0][1][2], 0.05900001967325228);
  assert( fabs(een_rescaled_e_gl[0][1][0][1][2] - (0.05900001967325228)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][1][0][2][0], 0.0);
  assert( fabs(een_rescaled_e_gl[0][1][0][2][0] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][1][0][2][1], -0.16473952654780463);
  assert( fabs(een_rescaled_e_gl[0][1][0][2][1] - (-0.16473952654780463)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][1][0][2][2], 0.06672545823691127);
  assert( fabs(een_rescaled_e_gl[0][1][0][2][2] - (0.06672545823691127)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][1][0][3][0], 0.0);
  assert( fabs(een_rescaled_e_gl[0][1][0][3][0] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][1][0][3][1], -0.7069765025301057);
  assert( fabs(een_rescaled_e_gl[0][1][0][3][1] - (-0.7069765025301057)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][1][0][3][2], -0.025889883334245332);
  assert( fabs(een_rescaled_e_gl[0][1][0][3][2] - (-0.025889883334245332)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][1][1][0][0], -0.15675618873333058);
  assert( fabs(een_rescaled_e_gl[0][1][1][0][0] - (-0.15675618873333058)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][1][1][0][1], 0.0);
  assert( fabs(een_rescaled_e_gl[0][1][1][0][1] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][1][1][0][2], -0.08997822461519353);
  assert( fabs(een_rescaled_e_gl[0][1][1][0][2] - (-0.08997822461519353)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][1][1][1][0], -0.2937567144823047);
  assert( fabs(een_rescaled_e_gl[0][1][1][1][0] - (-0.2937567144823047)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][1][1][1][1], 0.0);
  assert( fabs(een_rescaled_e_gl[0][1][1][1][1] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][1][1][1][2], 0.019899357798110473);
  assert( fabs(een_rescaled_e_gl[0][1][1][1][2] - (0.019899357798110473)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][1][1][2][0], 0.16473952654780463);
  assert( fabs(een_rescaled_e_gl[0][1][1][2][0] - (0.16473952654780463)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][1][1][2][1], 0.0);
  assert( fabs(een_rescaled_e_gl[0][1][1][2][1] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][1][1][2][2], 0.06619816955265034);
  assert( fabs(een_rescaled_e_gl[0][1][1][2][2] - (0.06619816955265034)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][1][1][3][0], -0.7069765025301057);
  assert( fabs(een_rescaled_e_gl[0][1][1][3][0] - (-0.7069765025301057)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][1][1][3][1], 0.0);
  assert( fabs(een_rescaled_e_gl[0][1][1][3][1] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][1][1][3][2], -0.013676100153650048);
  assert( fabs(een_rescaled_e_gl[0][1][1][3][2] - (-0.013676100153650048)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][1][2][0][0], 0.09831391870751387);
  assert( fabs(een_rescaled_e_gl[0][1][2][0][0] - (0.09831391870751387)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][1][2][0][1], 0.08997822461519353);
  assert( fabs(een_rescaled_e_gl[0][1][2][0][1] - (0.08997822461519353)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][1][2][0][2], 0.0);
  assert( fabs(een_rescaled_e_gl[0][1][2][0][2] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][1][2][1][0], -0.05900001967325228);
  assert( fabs(een_rescaled_e_gl[0][1][2][1][0] - (-0.05900001967325228)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][1][2][1][1], -0.019899357798110473);
  assert( fabs(een_rescaled_e_gl[0][1][2][1][1] - (-0.019899357798110473)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][1][2][1][2], 0.0);
  assert( fabs(een_rescaled_e_gl[0][1][2][1][2] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][1][2][2][0], -0.06672545823691127);
  assert( fabs(een_rescaled_e_gl[0][1][2][2][0] - (-0.06672545823691127)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][1][2][2][1], -0.06619816955265034);
  assert( fabs(een_rescaled_e_gl[0][1][2][2][1] - (-0.06619816955265034)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][1][2][2][2], 0.0);
  assert( fabs(een_rescaled_e_gl[0][1][2][2][2] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][1][2][3][0], -0.025889883334245332);
  assert( fabs(een_rescaled_e_gl[0][1][2][3][0] - (-0.025889883334245332)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][1][2][3][1], -0.013676100153650048);
  assert( fabs(een_rescaled_e_gl[0][1][2][3][1] - (-0.013676100153650048)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][1][2][3][2], 0.0);
  assert( fabs(een_rescaled_e_gl[0][1][2][3][2] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][2][0][0][0], 0.0);
  assert( fabs(een_rescaled_e_gl[0][2][0][0][0] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][2][0][0][1], 0.1941110489198177);
  assert( fabs(een_rescaled_e_gl[0][2][0][0][1] - (0.1941110489198177)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][2][0][0][2], -0.04347471142608776);
  assert( fabs(een_rescaled_e_gl[0][2][0][0][2] - (-0.04347471142608776)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][2][0][1][0], 0.0);
  assert( fabs(een_rescaled_e_gl[0][2][0][1][0] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][2][0][1][1], 0.36375867795818184);
  assert( fabs(een_rescaled_e_gl[0][2][0][1][1] - (0.36375867795818184)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][2][0][1][2], 0.026089986678886266);
  assert( fabs(een_rescaled_e_gl[0][2][0][1][2] - (0.026089986678886266)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][2][0][2][0], 0.0);
  assert( fabs(een_rescaled_e_gl[0][2][0][2][0] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][2][0][2][1], -0.20399680902645712);
  assert( fabs(een_rescaled_e_gl[0][2][0][2][1] - (-0.20399680902645712)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][2][0][2][2], 0.029506198916283078);
  assert( fabs(een_rescaled_e_gl[0][2][0][2][2] - (0.029506198916283078)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][2][0][3][0], 0.0);
  assert( fabs(een_rescaled_e_gl[0][2][0][3][0] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][2][0][3][1], -0.5994391302990586);
  assert( fabs(een_rescaled_e_gl[0][2][0][3][1] - (-0.5994391302990586)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][2][0][3][2], 0.023749246910207227);
  assert( fabs(een_rescaled_e_gl[0][2][0][3][2] - (0.023749246910207227)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][2][1][0][0], -0.1941110489198177);
  assert( fabs(een_rescaled_e_gl[0][2][1][0][0] - (-0.1941110489198177)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][2][1][0][1], 0.0);
  assert( fabs(een_rescaled_e_gl[0][2][1][0][1] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][2][1][0][2], -0.03403120705424245);
  assert( fabs(een_rescaled_e_gl[0][2][1][0][2] - (-0.03403120705424245)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][2][1][1][0], -0.36375867795818184);
  assert( fabs(een_rescaled_e_gl[0][2][1][1][0] - (-0.36375867795818184)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][2][1][1][1], 0.0);
  assert( fabs(een_rescaled_e_gl[0][2][1][1][1] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][2][1][1][2], 0.007526256139972797);
  assert( fabs(een_rescaled_e_gl[0][2][1][1][2] - (0.007526256139972797)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][2][1][2][0], 0.20399680902645712);
  assert( fabs(een_rescaled_e_gl[0][2][1][2][0] - (0.20399680902645712)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][2][1][2][1], 0.0);
  assert( fabs(een_rescaled_e_gl[0][2][1][2][1] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][2][1][2][2], 0.02503720899464923);
  assert( fabs(een_rescaled_e_gl[0][2][1][2][2] - (0.02503720899464923)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][2][1][3][0], -0.5994391302990586);
  assert( fabs(een_rescaled_e_gl[0][2][1][3][0] - (-0.5994391302990586)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][2][1][3][1], 0.0);
  assert( fabs(een_rescaled_e_gl[0][2][1][3][1] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][2][1][3][2], 0.020576005666223834);
  assert( fabs(een_rescaled_e_gl[0][2][1][3][2] - (0.020576005666223834)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][2][2][0][0], 0.04347471142608776);
  assert( fabs(een_rescaled_e_gl[0][2][2][0][0] - (0.04347471142608776)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][2][2][0][1], 0.03403120705424245);
  assert( fabs(een_rescaled_e_gl[0][2][2][0][1] - (0.03403120705424245)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][2][2][0][2], 0.0);
  assert( fabs(een_rescaled_e_gl[0][2][2][0][2] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][2][2][1][0], -0.026089986678886266);
  assert( fabs(een_rescaled_e_gl[0][2][2][1][0] - (-0.026089986678886266)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][2][2][1][1], -0.007526256139972797);
  assert( fabs(een_rescaled_e_gl[0][2][2][1][1] - (-0.007526256139972797)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][2][2][1][2], 0.0);
  assert( fabs(een_rescaled_e_gl[0][2][2][1][2] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][2][2][2][0], -0.029506198916283078);
  assert( fabs(een_rescaled_e_gl[0][2][2][2][0] - (-0.029506198916283078)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][2][2][2][1], -0.02503720899464923);
  assert( fabs(een_rescaled_e_gl[0][2][2][2][1] - (-0.02503720899464923)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][2][2][2][2], 0.0);
  assert( fabs(een_rescaled_e_gl[0][2][2][2][2] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][2][2][3][0], 0.023749246910207227);
  assert( fabs(een_rescaled_e_gl[0][2][2][3][0] - (0.023749246910207227)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][2][2][3][1], 0.020576005666223834);
  assert( fabs(een_rescaled_e_gl[0][2][2][3][1] - (0.020576005666223834)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][2][2][3][2], 0.0);
  assert( fabs(een_rescaled_e_gl[0][2][2][3][2] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][3][0][0][0], 0.0);
  assert( fabs(een_rescaled_e_gl[0][3][0][0][0] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][3][0][0][1], 0.18027565426866748);
  assert( fabs(een_rescaled_e_gl[0][3][0][0][1] - (0.18027565426866748)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][3][0][0][2], -0.014418486403775773);
  assert( fabs(een_rescaled_e_gl[0][3][0][0][2] - (-0.014418486403775773)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][3][0][1][0], 0.0);
  assert( fabs(een_rescaled_e_gl[0][3][0][1][0] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][3][0][1][1], 0.33783153524612014);
  assert( fabs(een_rescaled_e_gl[0][3][0][1][1] - (0.33783153524612014)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][3][0][1][2], 0.008652803109314725);
  assert( fabs(een_rescaled_e_gl[0][3][0][1][2] - (0.008652803109314725)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][3][0][2][0], 0.0);
  assert( fabs(een_rescaled_e_gl[0][3][0][2][0] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][3][0][2][1], -0.18945680022138287);
  assert( fabs(een_rescaled_e_gl[0][3][0][2][1] - (-0.18945680022138287)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][3][0][2][2], 0.009785797626853053);
  assert( fabs(een_rescaled_e_gl[0][3][0][2][2] - (0.009785797626853053)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][3][0][3][0], 0.0);
  assert( fabs(een_rescaled_e_gl[0][3][0][3][0] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][3][0][3][1], -0.3003772058582815);
  assert( fabs(een_rescaled_e_gl[0][3][0][3][1] - (-0.3003772058582815)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][3][0][3][2], 0.01954993189296513);
  assert( fabs(een_rescaled_e_gl[0][3][0][3][2] - (0.01954993189296513)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][3][1][0][0], -0.18027565426866748);
  assert( fabs(een_rescaled_e_gl[0][3][1][0][0] - (-0.18027565426866748)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][3][1][0][1], 0.0);
  assert( fabs(een_rescaled_e_gl[0][3][1][0][1] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][3][1][0][2], -0.009653361064760021);
  assert( fabs(een_rescaled_e_gl[0][3][1][0][2] - (-0.009653361064760021)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][3][1][1][0], -0.33783153524612014);
  assert( fabs(een_rescaled_e_gl[0][3][1][1][0] - (-0.33783153524612014)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][3][1][1][1], 0.0);
  assert( fabs(een_rescaled_e_gl[0][3][1][1][1] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][3][1][1][2], 0.0021349130481684514);
  assert( fabs(een_rescaled_e_gl[0][3][1][1][2] - (0.0021349130481684514)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][3][1][2][0], 0.18945680022138287);
  assert( fabs(een_rescaled_e_gl[0][3][1][2][0] - (0.18945680022138287)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][3][1][2][1], 0.0);
  assert( fabs(een_rescaled_e_gl[0][3][1][2][1] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][3][1][2][2], 0.007102105373281963);
  assert( fabs(een_rescaled_e_gl[0][3][1][2][2] - (0.007102105373281963)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][3][1][3][0], -0.3003772058582815);
  assert( fabs(een_rescaled_e_gl[0][3][1][3][0] - (-0.3003772058582815)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][3][1][3][1], 0.0);
  assert( fabs(een_rescaled_e_gl[0][3][1][3][1] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][3][1][3][2], 0.013140510401959491);
  assert( fabs(een_rescaled_e_gl[0][3][1][3][2] - (0.013140510401959491)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][3][2][0][0], 0.014418486403775773);
  assert( fabs(een_rescaled_e_gl[0][3][2][0][0] - (0.014418486403775773)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][3][2][0][1], 0.009653361064760021);
  assert( fabs(een_rescaled_e_gl[0][3][2][0][1] - (0.009653361064760021)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][3][2][0][2], 0.0);
  assert( fabs(een_rescaled_e_gl[0][3][2][0][2] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][3][2][1][0], -0.008652803109314725);
  assert( fabs(een_rescaled_e_gl[0][3][2][1][0] - (-0.008652803109314725)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][3][2][1][1], -0.0021349130481684514);
  assert( fabs(een_rescaled_e_gl[0][3][2][1][1] - (-0.0021349130481684514)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][3][2][1][2], 0.0);
  assert( fabs(een_rescaled_e_gl[0][3][2][1][2] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][3][2][2][0], -0.009785797626853053);
  assert( fabs(een_rescaled_e_gl[0][3][2][2][0] - (-0.009785797626853053)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][3][2][2][1], -0.007102105373281963);
  assert( fabs(een_rescaled_e_gl[0][3][2][2][1] - (-0.007102105373281963)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][3][2][2][2], 0.0);
  assert( fabs(een_rescaled_e_gl[0][3][2][2][2] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][3][2][3][0], 0.01954993189296513);
  assert( fabs(een_rescaled_e_gl[0][3][2][3][0] - (0.01954993189296513)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][3][2][3][1], 0.013140510401959491);
  assert( fabs(een_rescaled_e_gl[0][3][2][3][1] - (0.013140510401959491)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][3][2][3][2], 0.0);
  assert( fabs(een_rescaled_e_gl[0][3][2][3][2] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][4][0][0][0], 0.0);
  assert( fabs(een_rescaled_e_gl[0][4][0][0][0] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][4][0][0][1], 0.14882345476569966);
  assert( fabs(een_rescaled_e_gl[0][4][0][0][1] - (0.14882345476569966)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][4][0][0][2], -0.004250598788356597);
  assert( fabs(een_rescaled_e_gl[0][4][0][0][2] - (-0.004250598788356597)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][4][0][1][0], 0.0);
  assert( fabs(een_rescaled_e_gl[0][4][0][1][0] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][4][0][1][1], 0.2788909928414343);
  assert( fabs(een_rescaled_e_gl[0][4][0][1][1] - (0.2788909928414343)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][4][0][1][2], 0.0025508637579815512);
  assert( fabs(een_rescaled_e_gl[0][4][0][1][2] - (0.0025508637579815512)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][4][0][2][0], 0.0);
  assert( fabs(een_rescaled_e_gl[0][4][0][2][0] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][4][0][2][1], -0.1564027913374305);
  assert( fabs(een_rescaled_e_gl[0][4][0][2][1] - (-0.1564027913374305)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][4][0][2][2], 0.002884872820278267);
  assert( fabs(een_rescaled_e_gl[0][4][0][2][2] - (0.002884872820278267)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][4][0][3][0], 0.0);
  assert( fabs(een_rescaled_e_gl[0][4][0][3][0] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][4][0][3][1], -0.03635704449346779);
  assert( fabs(een_rescaled_e_gl[0][4][0][3][1] - (-0.03635704449346779)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][4][0][3][2], 0.009204712951216984);
  assert( fabs(een_rescaled_e_gl[0][4][0][3][2] - (0.009204712951216984)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][4][1][0][0], -0.14882345476569966);
  assert( fabs(een_rescaled_e_gl[0][4][1][0][0] - (-0.14882345476569966)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][4][1][0][1], 0.0);
  assert( fabs(een_rescaled_e_gl[0][4][1][0][1] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][4][1][0][2], -0.0024340372763829664);
  assert( fabs(een_rescaled_e_gl[0][4][1][0][2] - (-0.0024340372763829664)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][4][1][1][0], -0.2788909928414343);
  assert( fabs(een_rescaled_e_gl[0][4][1][1][0] - (-0.2788909928414343)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][4][1][1][1], 0.0);
  assert( fabs(een_rescaled_e_gl[0][4][1][1][1] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][4][1][1][2], 0.0005383055607490193);
  assert( fabs(een_rescaled_e_gl[0][4][1][1][2] - (0.0005383055607490193)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][4][1][2][0], 0.1564027913374305);
  assert( fabs(een_rescaled_e_gl[0][4][1][2][0] - (0.1564027913374305)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][4][1][2][1], 0.0);
  assert( fabs(een_rescaled_e_gl[0][4][1][2][1] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][4][1][2][2], 0.001790753407377889);
  assert( fabs(een_rescaled_e_gl[0][4][1][2][2] - (0.001790753407377889)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][4][1][3][0], -0.03635704449346779);
  assert( fabs(een_rescaled_e_gl[0][4][1][3][0] - (-0.03635704449346779)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][4][1][3][1], 0.0);
  assert( fabs(een_rescaled_e_gl[0][4][1][3][1] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][4][1][3][2], 0.005154930551874369);
  assert( fabs(een_rescaled_e_gl[0][4][1][3][2] - (0.005154930551874369)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][4][2][0][0], 0.004250598788356597);
  assert( fabs(een_rescaled_e_gl[0][4][2][0][0] - (0.004250598788356597)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][4][2][0][1], 0.0024340372763829664);
  assert( fabs(een_rescaled_e_gl[0][4][2][0][1] - (0.0024340372763829664)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][4][2][0][2], 0.0);
  assert( fabs(een_rescaled_e_gl[0][4][2][0][2] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][4][2][1][0], -0.0025508637579815512);
  assert( fabs(een_rescaled_e_gl[0][4][2][1][0] - (-0.0025508637579815512)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][4][2][1][1], -0.0005383055607490193);
  assert( fabs(een_rescaled_e_gl[0][4][2][1][1] - (-0.0005383055607490193)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][4][2][1][2], 0.0);
  assert( fabs(een_rescaled_e_gl[0][4][2][1][2] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][4][2][2][0], -0.002884872820278267);
  assert( fabs(een_rescaled_e_gl[0][4][2][2][0] - (-0.002884872820278267)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][4][2][2][1], -0.001790753407377889);
  assert( fabs(een_rescaled_e_gl[0][4][2][2][1] - (-0.001790753407377889)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][4][2][2][2], 0.0);
  assert( fabs(een_rescaled_e_gl[0][4][2][2][2] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][4][2][3][0], 0.009204712951216984);
  assert( fabs(een_rescaled_e_gl[0][4][2][3][0] - (0.009204712951216984)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][4][2][3][1], 0.005154930551874369);
  assert( fabs(een_rescaled_e_gl[0][4][2][3][1] - (0.005154930551874369)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][4][2][3][2], 0.0);
  assert( fabs(een_rescaled_e_gl[0][4][2][3][2] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][5][0][0][0], 0.0);
  assert( fabs(een_rescaled_e_gl[0][5][0][0][0] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][5][0][0][1], 0.11517996969478891);
  assert( fabs(een_rescaled_e_gl[0][5][0][0][1] - (0.11517996969478891)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][5][0][0][2], -0.0011747672541009072);
  assert( fabs(een_rescaled_e_gl[0][5][0][0][2] - (-0.0011747672541009072)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][5][0][1][0], 0.0);
  assert( fabs(een_rescaled_e_gl[0][5][0][1][0] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][5][0][1][1], 0.21584404255497447);
  assert( fabs(een_rescaled_e_gl[0][5][0][1][1] - (0.21584404255497447)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][5][0][1][2], 0.0007049997804446056);
  assert( fabs(een_rescaled_e_gl[0][5][0][1][2] - (0.0007049997804446056)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][5][0][2][0], 0.0);
  assert( fabs(een_rescaled_e_gl[0][5][0][2][0] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][5][0][2][1], -0.12104589827448056);
  assert( fabs(een_rescaled_e_gl[0][5][0][2][1] - (-0.12104589827448056)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][5][0][2][2], 0.0007973121647688946);
  assert( fabs(een_rescaled_e_gl[0][5][0][2][2] - (0.0007973121647688946)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][5][0][3][0], 0.0);
  assert( fabs(een_rescaled_e_gl[0][5][0][3][0] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][5][0][3][1], 0.13563796650527177);
  assert( fabs(een_rescaled_e_gl[0][5][0][3][1] - (0.13563796650527177)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][5][0][3][2], 0.0034950805168821176);
  assert( fabs(een_rescaled_e_gl[0][5][0][3][2] - (0.0034950805168821176)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][5][1][0][0], -0.11517996969478891);
  assert( fabs(een_rescaled_e_gl[0][5][1][0][0] - (-0.11517996969478891)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][5][1][0][1], 0.0);
  assert( fabs(een_rescaled_e_gl[0][5][1][0][1] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][5][1][0][2], -0.0005753699498169057);
  assert( fabs(een_rescaled_e_gl[0][5][1][0][2] - (-0.0005753699498169057)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][5][1][1][0], -0.21584404255497447);
  assert( fabs(een_rescaled_e_gl[0][5][1][1][0] - (-0.21584404255497447)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][5][1][1][1], 0.0);
  assert( fabs(een_rescaled_e_gl[0][5][1][1][1] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][5][1][1][2], 0.00012724737064609895);
  assert( fabs(een_rescaled_e_gl[0][5][1][1][2] - (0.00012724737064609895)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][5][1][2][0], 0.12104589827448056);
  assert( fabs(een_rescaled_e_gl[0][5][1][2][0] - (0.12104589827448056)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][5][1][2][1], 0.0);
  assert( fabs(een_rescaled_e_gl[0][5][1][2][1] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][5][1][2][2], 0.00042330727969317936);
  assert( fabs(een_rescaled_e_gl[0][5][1][2][2] - (0.00042330727969317936)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][5][1][3][0], 0.13563796650527177);
  assert( fabs(een_rescaled_e_gl[0][5][1][3][0] - (0.13563796650527177)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][5][1][3][1], 0.0);
  assert( fabs(een_rescaled_e_gl[0][5][1][3][1] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][5][1][3][2], 0.0016538819674466144);
  assert( fabs(een_rescaled_e_gl[0][5][1][3][2] - (0.0016538819674466144)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][5][2][0][0], 0.0011747672541009072);
  assert( fabs(een_rescaled_e_gl[0][5][2][0][0] - (0.0011747672541009072)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][5][2][0][1], 0.0005753699498169057);
  assert( fabs(een_rescaled_e_gl[0][5][2][0][1] - (0.0005753699498169057)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][5][2][0][2], 0.0);
  assert( fabs(een_rescaled_e_gl[0][5][2][0][2] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][5][2][1][0], -0.0007049997804446056);
  assert( fabs(een_rescaled_e_gl[0][5][2][1][0] - (-0.0007049997804446056)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][5][2][1][1], -0.00012724737064609895);
  assert( fabs(een_rescaled_e_gl[0][5][2][1][1] - (-0.00012724737064609895)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][5][2][1][2], 0.0);
  assert( fabs(een_rescaled_e_gl[0][5][2][1][2] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][5][2][2][0], -0.0007973121647688946);
  assert( fabs(een_rescaled_e_gl[0][5][2][2][0] - (-0.0007973121647688946)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][5][2][2][1], -0.00042330727969317936);
  assert( fabs(een_rescaled_e_gl[0][5][2][2][1] - (-0.00042330727969317936)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][5][2][2][2], 0.0);
  assert( fabs(een_rescaled_e_gl[0][5][2][2][2] - (0.0)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][5][2][3][0], 0.0034950805168821176);
  assert( fabs(een_rescaled_e_gl[0][5][2][3][0] - (0.0034950805168821176)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][5][2][3][1], 0.0016538819674466144);
  assert( fabs(een_rescaled_e_gl[0][5][2][3][1] - (0.0016538819674466144)) < 1.e-10 );
  printf( "%e %e\n", een_rescaled_e_gl[0][5][2][3][2], 0.0);
  assert( fabs(een_rescaled_e_gl[0][5][2][3][2] - (0.0)) < 1.e-10 );
}


{
  qmckl_context_struct* ctx = (qmckl_context_struct*) context;
  double een_rescaled_e_gl_doc[walk_num*(cord_num + 1)*elec_num*4*elec_num];
  memset(een_rescaled_e_gl_doc, 0, sizeof(een_rescaled_e_gl_doc));
  rc = qmckl_compute_jastrow_champ_factor_een_rescaled_e_gl_doc(context,
            ctx->electron.walker.num,
            ctx->electron.num,
            ctx->jastrow_champ.cord_num,
            ctx->jastrow_champ.rescale_factor_ee,
            ctx->electron.walker.point.coord.data,
            ctx->electron.ee_distance,
            ctx->jastrow_champ.een_rescaled_e,
            een_rescaled_e_gl_doc);
  assert(rc == QMCKL_SUCCESS);

  double een_rescaled_e_gl_hpc[walk_num*(cord_num + 1)*elec_num*4*elec_num];
  memset(een_rescaled_e_gl_hpc, 0, sizeof(een_rescaled_e_gl_hpc));
  rc = qmckl_compute_jastrow_champ_factor_een_rescaled_e_gl_hpc(context,
            ctx->electron.walker.num,
            ctx->electron.num,
            ctx->jastrow_champ.cord_num,
            ctx->jastrow_champ.rescale_factor_ee,
            ctx->electron.walker.point.coord.data,
            ctx->electron.ee_distance,
            ctx->jastrow_champ.een_rescaled_e,
            een_rescaled_e_gl_hpc);
  assert(rc == QMCKL_SUCCESS);

  for (int nw=0 ; nw < walk_num ; nw++) {
    for (int c=0 ; c <= cord_num ; c++) {
      for (int i=0 ; i < elec_num ; i++) {
        for (int j=0 ; j < elec_num ; j++) {
          for (int k=0 ; k < 4 ; k++) {
            if (fabs(een_rescaled_e_gl_doc[nw*(cord_num + 1)*elec_num*4*elec_num + c*elec_num*4*elec_num + i*4*elec_num + k*elec_num + j] - een_rescaled_e_gl_hpc[nw*(cord_num + 1)*elec_num*4*elec_num + c*elec_num*4*elec_num + i*4*elec_num + k*elec_num + j]) > 1.e-12) {
              printf("nw=%d c=%d i=%d k=%d j=%d doc=%e hpc=%e\n", nw, c, i, k, j, een_rescaled_e_gl_doc[nw*(cord_num + 1)*elec_num*4*elec_num + c*elec_num*4*elec_num + i*4*elec_num + k*elec_num + j], een_rescaled_e_gl_hpc[nw*(cord_num + 1)*elec_num*4*elec_num + c*elec_num*4*elec_num + i*4*elec_num + k*elec_num + j]);
              fflush(stdout);
            }
            assert(fabs(een_rescaled_e_gl_doc[nw*(cord_num + 1)*elec_num*4*elec_num + c*elec_num*4*elec_num + i*4*elec_num + k*elec_num + j] - een_rescaled_e_gl_hpc[nw*(cord_num + 1)*elec_num*4*elec_num + c*elec_num*4*elec_num + i*4*elec_num + k*elec_num + j]) < 1.e-8);
          }
        }
      }
    }
  }

}

{
  /* Finite difference test fails and I can't understand why... */

  printf("een_rescaled_e_gl\n");

  double fd[walk_num][cord_num+1][elec_num][4][elec_num];

  double delta_x = 0.001;

  // Finite difference coefficients for gradients
  double coef[9] = { 1.0/280.0, -4.0/105.0,  1.0/5.0, -4.0/5.0,         0.0, 4.0/5.0, -1.0/5.0, 4.0/105.0, -1.0/280.0 };

  // Finite difference coefficients for Laplacian
  double coef2[9]= {-1.0/560.0,	 8.0/315.0, -1.0/5.0,  8.0/5.0, -205.0/72.0, 8.0/5.0, -1.0/5.0,	8.0/315.0, -1.0/560.0 };

  qmckl_exit_code rc;

  double elec_coord[walk_num][elec_num][3];
  rc = qmckl_get_electron_coord (context, 'N', &(elec_coord[0][0][0]), 3*walk_num*elec_num);
  assert (rc == QMCKL_SUCCESS);

  double temp_coord[walk_num][elec_num][3];
  memcpy(&(temp_coord[0][0][0]), &(elec_coord[0][0][0]), sizeof(temp_coord));

  double function_values[walk_num][cord_num+1][elec_num][elec_num];
  memset(&(fd[0][0][0][0]), 0, sizeof(fd));

  for (int64_t i = 0; i < elec_num; i++) {
    for (int64_t k = 0; k < 3; k++) {
      for (int64_t m = -4; m <= 4; m++) {          // Apply finite difference displacement

        memcpy(&(temp_coord[0][0][0]), &(elec_coord[0][0][0]), sizeof(temp_coord));
        for (int64_t nw=0 ; nw<walk_num ; nw++) {
          temp_coord[nw][i][k] = elec_coord[nw][i][k] + (double) m * delta_x;
        }

        // Update coordinates in the context
        rc = qmckl_set_electron_coord (context, 'N', walk_num,
                                       &(temp_coord[0][0][0]),
                                       walk_num*3*elec_num);
        assert(rc == QMCKL_SUCCESS);

        // Call the provided function
        rc = qmckl_get_jastrow_champ_een_rescaled_e(context,
                                                    &(function_values[0][0][0][0]),
                                                    walk_num*(cord_num+1)*elec_num*elec_num);
        assert(rc == QMCKL_SUCCESS);

        // Accumulate derivative using finite-difference coefficients
        for (int64_t nw=0 ; nw<walk_num ; nw++) {
          for (int64_t c = 0; c < cord_num+1 ; c++) {
            for (int64_t j = 0; j < elec_num; j++) {
            fd[nw][c][j][k][i] += coef [m + 4] * function_values[nw][c][j][i];
            fd[nw][c][j][3][i] += coef2[m + 4] * function_values[nw][c][j][i];
            }
          }
        }

      }
    }
  }

  // Reset coordinates in the context
  rc = qmckl_set_electron_coord (context, 'N', walk_num,
                                 &(elec_coord[0][0][0]),
                                 walk_num*3*elec_num);
  assert(rc == QMCKL_SUCCESS);

  // Normalize by the step size
  for (int64_t nw=0 ; nw<walk_num ; nw++) {
    for (int64_t c = 0; c < cord_num+1 ; c++) {
      for (int64_t i = 0; i < elec_num; i++) {
        for (int64_t k = 0; k < 4; k++) {
          for (int64_t j = 0; j < elec_num; j++) {
            fd[nw][c][i][k][j] /= delta_x;
          }
        }
        for (int64_t j = 0; j < elec_num; j++) {
          fd[nw][c][i][3][j] /= delta_x;
        }
      }
    }
  }


  double een_rescaled_e_gl[walk_num][cord_num+1][elec_num][4][elec_num];
  rc = qmckl_check(context,
                   qmckl_get_jastrow_champ_een_rescaled_e_gl(context,
                                                                      &(een_rescaled_e_gl[0][0][0][0][0]),
                                                                      walk_num*(cord_num+1)*elec_num*4*elec_num)
                   );


  assert(rc == QMCKL_SUCCESS);

  for (int nw = 0; nw < walk_num; nw++){
    for (int c = 0; c < cord_num+1 ; c++) {
      for (int i = 0; i < elec_num; i++) {
        for (int j = 0; j < elec_num; j++) {
          for (int k = 0; k < 3; k++){
            if (fabs(fd[nw][c][i][k][j] - een_rescaled_e_gl[nw][c][i][k][j]) > 1.e-10) {
              printf("%2d %2d %2d %2d %2d\t", nw, c, i, k, j);
              printf("%.10e\t", fd[nw][c][i][k][j]);
              printf("%.10e\n", een_rescaled_e_gl[nw][c][i][k][j]);
              fflush(stdout);
            }
            assert(fabs(fd[nw][c][i][k][j] - een_rescaled_e_gl[nw][c][i][k][j]) < 1.e-8);
          }
          int k=3;
          if (i != j) {
            if (fabs(fd[nw][c][i][k][j] - een_rescaled_e_gl[nw][c][i][k][j]) > 1.e-8) {
              printf("%2d %2d %2d %2d %2d\t", nw, c, i, k, j);
              printf("%.10e\t", fd[nw][c][i][k][j]);
              printf("%.10e\n", een_rescaled_e_gl[nw][c][i][k][j]);
              fflush(stdout);
            }
            assert(fabs(fd[nw][c][i][k][j] - een_rescaled_e_gl[nw][c][i][k][j]) < 1.e-6);
          }
        }
      }
    }
  }
  printf("OK\n");
}

assert(qmckl_electron_provided(context));

double een_rescaled_n[walk_num][(cord_num + 1)][nucl_num][elec_num];
size_max=walk_num*(cord_num + 1)*nucl_num*elec_num;
rc = qmckl_get_jastrow_champ_een_rescaled_n(context, &(een_rescaled_n[0][0][0][0]),size_max);

// value of (0,2,1)
assert(fabs(een_rescaled_n[0][1][0][2]-0.2603169838750542 )< 1.e-12);
assert(fabs(een_rescaled_n[0][1][0][3]-0.3016180139679065 )< 1.e-12);
assert(fabs(een_rescaled_n[0][1][0][4]-0.10506023826192266)< 1.e-12);
assert(fabs(een_rescaled_n[0][2][1][3]-0.9267719759374164 )< 1.e-12);
assert(fabs(een_rescaled_n[0][2][1][4]-0.11497585238132658)< 1.e-12);
assert(fabs(een_rescaled_n[0][2][1][5]-0.07534033469115217)< 1.e-12);

assert(qmckl_electron_provided(context));

double een_rescaled_n_gl[walk_num][(cord_num + 1)][nucl_num][4][elec_num];
size_max=walk_num*(cord_num + 1)*nucl_num*4*elec_num;
rc = qmckl_get_jastrow_champ_een_rescaled_n_gl(context, &(een_rescaled_n_gl[0][0][0][0][0]),size_max);

// value of (0,2,1)
assert(fabs( -0.11234061209936878  - een_rescaled_n_gl[0][1][0][0][2])  < 1.e-12);
assert(fabs( 0.0004440109367151707 - een_rescaled_n_gl[0][1][0][0][3])  < 1.e-12);
assert(fabs( -0.012868642597346566 - een_rescaled_n_gl[0][1][0][0][4])  < 1.e-12);
assert(fabs( 0.08601122289922644   - een_rescaled_n_gl[0][2][1][0][3])  < 1.e-12);
assert(fabs( -0.058681563677207206 - een_rescaled_n_gl[0][2][1][0][4])  < 1.e-12);
assert(fabs( 0.005359281880312882  - een_rescaled_n_gl[0][2][1][0][5])  < 1.e-12);

{
  assert(qmckl_electron_provided(context));

  printf("tmp_c\n");
  double tmp_c[walk_num][cord_num][cord_num+1][nucl_num][elec_num];
  rc = qmckl_get_jastrow_champ_tmp_c(context, &(tmp_c[0][0][0][0][0]), sizeof(tmp_c)/sizeof(double));

  
}

{
  printf("dtmp_c\n");
  double dtmp_c[walk_num][cord_num][cord_num+1][nucl_num][4][elec_num];
  rc = qmckl_get_jastrow_champ_dtmp_c(context, &(dtmp_c[0][0][0][0][0][0]), sizeof(dtmp_c)/sizeof(double));

  printf("%e %e\n", dtmp_c[0][0][0][0][0][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][0][0][0][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][0][0][0][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][0][0][0][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][0][0][0][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][0][0][0][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][0][0][1][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][0][0][1][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][0][0][1][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][0][0][1][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][0][0][1][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][0][0][1][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][0][0][2][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][0][0][2][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][0][0][2][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][0][0][2][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][0][0][2][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][0][0][2][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][0][0][3][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][0][0][3][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][0][0][3][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][0][0][3][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][0][0][3][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][0][0][3][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][0][1][0][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][0][1][0][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][0][1][0][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][0][1][0][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][0][1][0][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][0][1][0][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][0][1][1][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][0][1][1][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][0][1][1][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][0][1][1][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][0][1][1][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][0][1][1][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][0][1][2][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][0][1][2][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][0][1][2][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][0][1][2][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][0][1][2][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][0][1][2][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][0][1][3][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][0][1][3][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][0][1][3][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][0][1][3][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][0][1][3][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][0][1][3][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][1][0][0][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][1][0][0][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][1][0][0][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][1][0][0][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][1][0][0][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][1][0][0][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][1][0][1][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][1][0][1][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][1][0][1][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][1][0][1][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][1][0][1][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][1][0][1][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][1][0][2][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][1][0][2][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][1][0][2][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][1][0][2][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][1][0][2][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][1][0][2][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][1][0][3][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][1][0][3][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][1][0][3][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][1][0][3][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][1][0][3][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][1][0][3][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][1][1][0][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][1][1][0][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][1][1][0][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][1][1][0][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][1][1][0][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][1][1][0][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][1][1][1][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][1][1][1][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][1][1][1][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][1][1][1][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][1][1][1][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][1][1][1][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][1][1][2][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][1][1][2][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][1][1][2][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][1][1][2][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][1][1][2][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][1][1][2][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][1][1][3][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][1][1][3][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][1][1][3][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][1][1][3][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][1][1][3][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][1][1][3][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][2][0][0][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][2][0][0][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][2][0][0][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][2][0][0][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][2][0][0][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][2][0][0][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][2][0][1][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][2][0][1][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][2][0][1][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][2][0][1][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][2][0][1][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][2][0][1][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][2][0][2][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][2][0][2][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][2][0][2][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][2][0][2][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][2][0][2][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][2][0][2][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][2][0][3][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][2][0][3][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][2][0][3][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][2][0][3][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][2][0][3][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][2][0][3][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][2][1][0][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][2][1][0][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][2][1][0][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][2][1][0][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][2][1][0][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][2][1][0][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][2][1][1][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][2][1][1][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][2][1][1][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][2][1][1][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][2][1][1][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][2][1][1][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][2][1][2][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][2][1][2][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][2][1][2][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][2][1][2][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][2][1][2][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][2][1][2][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][2][1][3][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][2][1][3][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][2][1][3][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][2][1][3][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][2][1][3][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][2][1][3][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][3][0][0][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][3][0][0][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][3][0][0][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][3][0][0][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][3][0][0][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][3][0][0][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][3][0][1][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][3][0][1][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][3][0][1][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][3][0][1][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][3][0][1][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][3][0][1][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][3][0][2][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][3][0][2][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][3][0][2][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][3][0][2][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][3][0][2][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][3][0][2][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][3][0][3][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][3][0][3][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][3][0][3][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][3][0][3][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][3][0][3][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][3][0][3][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][3][1][0][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][3][1][0][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][3][1][0][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][3][1][0][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][3][1][0][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][3][1][0][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][3][1][1][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][3][1][1][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][3][1][1][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][3][1][1][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][3][1][1][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][3][1][1][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][3][1][2][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][3][1][2][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][3][1][2][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][3][1][2][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][3][1][2][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][3][1][2][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][3][1][3][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][3][1][3][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][3][1][3][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][3][1][3][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][3][1][3][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][3][1][3][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][4][0][0][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][4][0][0][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][4][0][0][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][4][0][0][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][4][0][0][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][4][0][0][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][4][0][1][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][4][0][1][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][4][0][1][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][4][0][1][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][4][0][1][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][4][0][1][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][4][0][2][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][4][0][2][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][4][0][2][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][4][0][2][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][4][0][2][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][4][0][2][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][4][0][3][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][4][0][3][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][4][0][3][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][4][0][3][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][4][0][3][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][4][0][3][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][4][1][0][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][4][1][0][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][4][1][0][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][4][1][0][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][4][1][0][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][4][1][0][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][4][1][1][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][4][1][1][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][4][1][1][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][4][1][1][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][4][1][1][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][4][1][1][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][4][1][2][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][4][1][2][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][4][1][2][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][4][1][2][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][4][1][2][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][4][1][2][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][4][1][3][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][4][1][3][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][4][1][3][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][4][1][3][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][4][1][3][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][4][1][3][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][5][0][0][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][5][0][0][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][5][0][0][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][5][0][0][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][5][0][0][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][5][0][0][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][5][0][1][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][5][0][1][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][5][0][1][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][5][0][1][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][5][0][1][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][5][0][1][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][5][0][2][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][5][0][2][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][5][0][2][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][5][0][2][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][5][0][2][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][5][0][2][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][5][0][3][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][5][0][3][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][5][0][3][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][5][0][3][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][5][0][3][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][5][0][3][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][5][1][0][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][5][1][0][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][5][1][0][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][5][1][0][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][5][1][0][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][5][1][0][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][5][1][1][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][5][1][1][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][5][1][1][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][5][1][1][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][5][1][1][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][5][1][1][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][5][1][2][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][5][1][2][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][5][1][2][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][5][1][2][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][5][1][2][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][5][1][2][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][5][1][3][0], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][5][1][3][0] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][5][1][3][1], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][5][1][3][1] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][0][5][1][3][2], 0.0); fflush(stdout);
  assert( fabs(dtmp_c[0][0][5][1][3][2] - (0.0)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][0][0][0][0], 0.32786569627746437); fflush(stdout);
  assert( fabs(dtmp_c[0][1][0][0][0][0] - (0.32786569627746437)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][0][0][0][1], 1.175778273532559); fflush(stdout);
  assert( fabs(dtmp_c[0][1][0][0][0][1] - (1.175778273532559)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][0][0][0][2], -0.5891750295764052); fflush(stdout);
  assert( fabs(dtmp_c[0][1][0][0][0][2] - (-0.5891750295764052)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][0][0][1][0], -1.2752658765514837); fflush(stdout);
  assert( fabs(dtmp_c[0][1][0][0][1][0] - (-1.2752658765514837)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][0][0][1][1], 0.3677002476817276); fflush(stdout);
  assert( fabs(dtmp_c[0][1][0][0][1][1] - (0.3677002476817276)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][0][0][1][2], 0.27007596029075887); fflush(stdout);
  assert( fabs(dtmp_c[0][1][0][0][1][2] - (0.27007596029075887)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][0][0][2][0], 0.38986608278267915); fflush(stdout);
  assert( fabs(dtmp_c[0][1][0][0][2][0] - (0.38986608278267915)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][0][0][2][1], -0.281438284901225); fflush(stdout);
  assert( fabs(dtmp_c[0][1][0][0][2][1] - (-0.281438284901225)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][0][0][2][2], 0.6671492735630229); fflush(stdout);
  assert( fabs(dtmp_c[0][1][0][0][2][2] - (0.6671492735630229)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][0][0][3][0], -3.1027358156582245); fflush(stdout);
  assert( fabs(dtmp_c[0][1][0][0][3][0] - (-3.1027358156582245)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][0][0][3][1], -3.105256527340169); fflush(stdout);
  assert( fabs(dtmp_c[0][1][0][0][3][1] - (-3.105256527340169)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][0][0][3][2], -0.26112735005369025); fflush(stdout);
  assert( fabs(dtmp_c[0][1][0][0][3][2] - (-0.26112735005369025)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][0][1][0][0], 0.32786569627746437); fflush(stdout);
  assert( fabs(dtmp_c[0][1][0][1][0][0] - (0.32786569627746437)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][0][1][0][1], 1.175778273532559); fflush(stdout);
  assert( fabs(dtmp_c[0][1][0][1][0][1] - (1.175778273532559)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][0][1][0][2], -0.5891750295764052); fflush(stdout);
  assert( fabs(dtmp_c[0][1][0][1][0][2] - (-0.5891750295764052)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][0][1][1][0], -1.2752658765514837); fflush(stdout);
  assert( fabs(dtmp_c[0][1][0][1][1][0] - (-1.2752658765514837)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][0][1][1][1], 0.3677002476817276); fflush(stdout);
  assert( fabs(dtmp_c[0][1][0][1][1][1] - (0.3677002476817276)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][0][1][1][2], 0.27007596029075887); fflush(stdout);
  assert( fabs(dtmp_c[0][1][0][1][1][2] - (0.27007596029075887)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][0][1][2][0], 0.38986608278267915); fflush(stdout);
  assert( fabs(dtmp_c[0][1][0][1][2][0] - (0.38986608278267915)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][0][1][2][1], -0.281438284901225); fflush(stdout);
  assert( fabs(dtmp_c[0][1][0][1][2][1] - (-0.281438284901225)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][0][1][2][2], 0.6671492735630229); fflush(stdout);
  assert( fabs(dtmp_c[0][1][0][1][2][2] - (0.6671492735630229)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][0][1][3][0], -3.1027358156582245); fflush(stdout);
  assert( fabs(dtmp_c[0][1][0][1][3][0] - (-3.1027358156582245)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][0][1][3][1], -3.105256527340169); fflush(stdout);
  assert( fabs(dtmp_c[0][1][0][1][3][1] - (-3.105256527340169)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][0][1][3][2], -0.26112735005369025); fflush(stdout);
  assert( fabs(dtmp_c[0][1][0][1][3][2] - (-0.26112735005369025)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][1][0][0][0], 0.12974942861337901); fflush(stdout);
  assert( fabs(dtmp_c[0][1][1][0][0][0] - (0.12974942861337901)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][1][0][0][1], 0.737321508542927); fflush(stdout);
  assert( fabs(dtmp_c[0][1][1][0][0][1] - (0.737321508542927)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][1][0][0][2], -0.3933232128531056); fflush(stdout);
  assert( fabs(dtmp_c[0][1][1][0][0][2] - (-0.3933232128531056)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][1][0][1][0], -0.8815381591013682); fflush(stdout);
  assert( fabs(dtmp_c[0][1][1][0][1][0] - (-0.8815381591013682)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][1][0][1][1], 0.25207342867289106); fflush(stdout);
  assert( fabs(dtmp_c[0][1][1][0][1][1] - (0.25207342867289106)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][1][0][1][2], 0.16884925476964702); fflush(stdout);
  assert( fabs(dtmp_c[0][1][1][0][1][2] - (0.16884925476964702)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][1][0][2][0], 0.09968893609256244); fflush(stdout);
  assert( fabs(dtmp_c[0][1][1][0][2][0] - (0.09968893609256244)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][1][0][2][1], -0.44754252122173505); fflush(stdout);
  assert( fabs(dtmp_c[0][1][1][0][2][1] - (-0.44754252122173505)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][1][0][2][2], 0.36893090614395035); fflush(stdout);
  assert( fabs(dtmp_c[0][1][1][0][2][2] - (0.36893090614395035)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][1][0][3][0], -2.2764908066925345); fflush(stdout);
  assert( fabs(dtmp_c[0][1][1][0][3][0] - (-2.2764908066925345)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][1][0][3][1], -2.3370447999714297); fflush(stdout);
  assert( fabs(dtmp_c[0][1][1][0][3][1] - (-2.3370447999714297)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][1][0][3][2], -0.17618867383487768); fflush(stdout);
  assert( fabs(dtmp_c[0][1][1][0][3][2] - (-0.17618867383487768)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][1][1][0][0], 0.08655901222441939); fflush(stdout);
  assert( fabs(dtmp_c[0][1][1][1][0][0] - (0.08655901222441939)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][1][1][0][1], 0.3718476223462695); fflush(stdout);
  assert( fabs(dtmp_c[0][1][1][1][0][1] - (0.3718476223462695)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][1][1][0][2], -0.1880144254085473); fflush(stdout);
  assert( fabs(dtmp_c[0][1][1][1][0][2] - (-0.1880144254085473)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][1][1][1][0], -0.373219730383185); fflush(stdout);
  assert( fabs(dtmp_c[0][1][1][1][1][0] - (-0.373219730383185)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][1][1][1][1], 0.11557326827726634); fflush(stdout);
  assert( fabs(dtmp_c[0][1][1][1][1][1] - (0.11557326827726634)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][1][1][1][2], 0.08237999457321925); fflush(stdout);
  assert( fabs(dtmp_c[0][1][1][1][1][2] - (0.08237999457321925)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][1][1][2][0], 0.37003474273905146); fflush(stdout);
  assert( fabs(dtmp_c[0][1][1][1][2][0] - (0.37003474273905146)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][1][1][2][1], 0.23459398489351002); fflush(stdout);
  assert( fabs(dtmp_c[0][1][1][1][2][1] - (0.23459398489351002)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][1][1][2][2], 0.25766171588353404); fflush(stdout);
  assert( fabs(dtmp_c[0][1][1][1][2][2] - (0.25766171588353404)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][1][1][3][0], -0.8539076659161222); fflush(stdout);
  assert( fabs(dtmp_c[0][1][1][1][3][0] - (-0.8539076659161222)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][1][1][3][1], -0.8931184537969846); fflush(stdout);
  assert( fabs(dtmp_c[0][1][1][1][3][1] - (-0.8931184537969846)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][1][1][3][2], -0.05441849547515038); fflush(stdout);
  assert( fabs(dtmp_c[0][1][1][1][3][2] - (-0.05441849547515038)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][2][0][0][0], 0.06497828699933875); fflush(stdout);
  assert( fabs(dtmp_c[0][1][2][0][0][0] - (0.06497828699933875)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][2][0][0][1], 0.5299426253129281); fflush(stdout);
  assert( fabs(dtmp_c[0][1][2][0][0][1] - (0.5299426253129281)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][2][0][0][2], -0.28147113511754035); fflush(stdout);
  assert( fabs(dtmp_c[0][1][2][0][0][2] - (-0.28147113511754035)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][2][0][1][0], -0.656536894417369); fflush(stdout);
  assert( fabs(dtmp_c[0][1][2][0][1][0] - (-0.656536894417369)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][2][0][1][1], 0.17348835520190764); fflush(stdout);
  assert( fabs(dtmp_c[0][1][2][0][1][1] - (0.17348835520190764)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][2][0][1][2], 0.11300260357985709); fflush(stdout);
  assert( fabs(dtmp_c[0][1][2][0][1][2] - (0.11300260357985709)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][2][0][2][0], 0.030118259925656102); fflush(stdout);
  assert( fabs(dtmp_c[0][1][2][0][2][0] - (0.030118259925656102)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][2][0][2][1], -0.400859507927543); fflush(stdout);
  assert( fabs(dtmp_c[0][1][2][0][2][1] - (-0.400859507927543)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][2][0][2][2], 0.23669854037165033); fflush(stdout);
  assert( fabs(dtmp_c[0][1][2][0][2][2] - (0.23669854037165033)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][2][0][3][0], -1.7260664383867277); fflush(stdout);
  assert( fabs(dtmp_c[0][1][2][0][3][0] - (-1.7260664383867277)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][2][0][3][1], -1.8326386990480228); fflush(stdout);
  assert( fabs(dtmp_c[0][1][2][0][3][1] - (-1.8326386990480228)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][2][0][3][2], -0.11974756473481223); fflush(stdout);
  assert( fabs(dtmp_c[0][1][2][0][3][2] - (-0.11974756473481223)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][2][1][0][0], 0.04101895608273089); fflush(stdout);
  assert( fabs(dtmp_c[0][1][2][1][0][0] - (0.04101895608273089)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][2][1][0][1], 0.1692098996461025); fflush(stdout);
  assert( fabs(dtmp_c[0][1][2][1][0][1] - (0.1692098996461025)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][2][1][0][2], -0.08157456325334705); fflush(stdout);
  assert( fabs(dtmp_c[0][1][2][1][0][2] - (-0.08157456325334705)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][2][1][1][0], -0.14062811888255994); fflush(stdout);
  assert( fabs(dtmp_c[0][1][2][1][1][0] - (-0.14062811888255994)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][2][1][1][1], 0.04966346199008702); fflush(stdout);
  assert( fabs(dtmp_c[0][1][2][1][1][1] - (0.04966346199008702)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][2][1][1][2], 0.03397180433339092); fflush(stdout);
  assert( fabs(dtmp_c[0][1][2][1][1][2] - (0.03397180433339092)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][2][1][2][0], 0.3092486076246164); fflush(stdout);
  assert( fabs(dtmp_c[0][1][2][1][2][0] - (0.3092486076246164)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][2][1][2][1], 0.3103297743123867); fflush(stdout);
  assert( fabs(dtmp_c[0][1][2][1][2][1] - (0.3103297743123867)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][2][1][2][2], 0.14068463168890066); fflush(stdout);
  assert( fabs(dtmp_c[0][1][2][1][2][2] - (0.14068463168890066)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][2][1][3][0], -0.2771885915793018); fflush(stdout);
  assert( fabs(dtmp_c[0][1][2][1][3][0] - (-0.2771885915793018)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][2][1][3][1], -0.32750955178890273); fflush(stdout);
  assert( fabs(dtmp_c[0][1][2][1][3][1] - (-0.32750955178890273)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][2][1][3][2], -0.005414624476035229); fflush(stdout);
  assert( fabs(dtmp_c[0][1][2][1][3][2] - (-0.005414624476035229)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][3][0][0][0], 0.041518317038679226); fflush(stdout);
  assert( fabs(dtmp_c[0][1][3][0][0][0] - (0.041518317038679226)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][3][0][0][1], 0.40808685616423107); fflush(stdout);
  assert( fabs(dtmp_c[0][1][3][0][0][1] - (0.40808685616423107)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][3][0][0][2], -0.20936719974741586); fflush(stdout);
  assert( fabs(dtmp_c[0][1][3][0][0][2] - (-0.20936719974741586)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][3][0][1][0], -0.5079036503920699); fflush(stdout);
  assert( fabs(dtmp_c[0][1][3][0][1][0] - (-0.5079036503920699)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][3][0][1][1], 0.11990591460399942); fflush(stdout);
  assert( fabs(dtmp_c[0][1][3][0][1][1] - (0.11990591460399942)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][3][0][1][2], 0.07916578219778489); fflush(stdout);
  assert( fabs(dtmp_c[0][1][3][0][1][2] - (0.07916578219778489)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][3][0][2][0], 0.013175750148178606); fflush(stdout);
  assert( fabs(dtmp_c[0][1][3][0][2][0] - (0.013175750148178606)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][3][0][2][1], -0.3270385104723065); fflush(stdout);
  assert( fabs(dtmp_c[0][1][3][0][2][1] - (-0.3270385104723065)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][3][0][2][2], 0.16494307000507016); fflush(stdout);
  assert( fabs(dtmp_c[0][1][3][0][2][2] - (0.16494307000507016)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][3][0][3][0], -1.3394133357634042); fflush(stdout);
  assert( fabs(dtmp_c[0][1][3][0][3][0] - (-1.3394133357634042)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][3][0][3][1], -1.472154587343); fflush(stdout);
  assert( fabs(dtmp_c[0][1][3][0][3][1] - (-1.472154587343)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][3][0][3][2], -0.08409179172949051); fflush(stdout);
  assert( fabs(dtmp_c[0][1][3][0][3][2] - (-0.08409179172949051)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][3][1][0][0], 0.029321550458560004); fflush(stdout);
  assert( fabs(dtmp_c[0][1][3][1][0][0] - (0.029321550458560004)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][3][1][0][1], 0.10968727103181443); fflush(stdout);
  assert( fabs(dtmp_c[0][1][3][1][0][1] - (0.10968727103181443)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][3][1][0][2], -0.05056277489728557); fflush(stdout);
  assert( fabs(dtmp_c[0][1][3][1][0][2] - (-0.05056277489728557)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][3][1][1][0], -0.07490965750749531); fflush(stdout);
  assert( fabs(dtmp_c[0][1][3][1][1][0] - (-0.07490965750749531)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][3][1][1][1], 0.03161213358508026); fflush(stdout);
  assert( fabs(dtmp_c[0][1][3][1][1][1] - (0.03161213358508026)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][3][1][1][2], 0.02042682750021449); fflush(stdout);
  assert( fabs(dtmp_c[0][1][3][1][1][2] - (0.02042682750021449)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][3][1][2][0], 0.2651287007442243); fflush(stdout);
  assert( fabs(dtmp_c[0][1][3][1][2][0] - (0.2651287007442243)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][3][1][2][1], 0.29864100133392585); fflush(stdout);
  assert( fabs(dtmp_c[0][1][3][1][2][1] - (0.29864100133392585)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][3][1][2][2], 0.10131396651000343); fflush(stdout);
  assert( fabs(dtmp_c[0][1][3][1][2][2] - (0.10131396651000343)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][3][1][3][0], -0.12233799123812668); fflush(stdout);
  assert( fabs(dtmp_c[0][1][3][1][3][0] - (-0.12233799123812668)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][3][1][3][1], -0.1749982805145044); fflush(stdout);
  assert( fabs(dtmp_c[0][1][3][1][3][1] - (-0.1749982805145044)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][3][1][3][2], 0.005533912938599828); fflush(stdout);
  assert( fabs(dtmp_c[0][1][3][1][3][2] - (0.005533912938599828)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][4][0][0][0], 0.03228007764178889); fflush(stdout);
  assert( fabs(dtmp_c[0][1][4][0][0][0] - (0.03228007764178889)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][4][0][0][1], 0.32650132646618324); fflush(stdout);
  assert( fabs(dtmp_c[0][1][4][0][0][1] - (0.32650132646618324)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][4][0][0][2], -0.15983679736539755); fflush(stdout);
  assert( fabs(dtmp_c[0][1][4][0][0][2] - (-0.15983679736539755)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][4][0][1][0], -0.4023830807970329); fflush(stdout);
  assert( fabs(dtmp_c[0][1][4][0][1][0] - (-0.4023830807970329)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][4][0][1][1], 0.08299585383989319); fflush(stdout);
  assert( fabs(dtmp_c[0][1][4][0][1][1] - (0.08299585383989319)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][4][0][1][2], 0.057417997456280095); fflush(stdout);
  assert( fabs(dtmp_c[0][1][4][0][1][2] - (0.057417997456280095)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][4][0][2][0], 0.009710514963830302); fflush(stdout);
  assert( fabs(dtmp_c[0][1][4][0][2][0] - (0.009710514963830302)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][4][0][2][1], -0.2619776157592431); fflush(stdout);
  assert( fabs(dtmp_c[0][1][4][0][2][1] - (-0.2619776157592431)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][4][0][2][2], 0.12112340544952337); fflush(stdout);
  assert( fabs(dtmp_c[0][1][4][0][2][2] - (0.12112340544952337)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][4][0][3][0], -1.0586324701186223); fflush(stdout);
  assert( fabs(dtmp_c[0][1][4][0][3][0] - (-1.0586324701186223)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][4][0][3][1], -1.2022846698870553); fflush(stdout);
  assert( fabs(dtmp_c[0][1][4][0][3][1] - (-1.2022846698870553)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][4][0][3][2], -0.06127569538128786); fflush(stdout);
  assert( fabs(dtmp_c[0][1][4][0][3][2] - (-0.06127569538128786)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][4][1][0][0], 0.024568482604823683); fflush(stdout);
  assert( fabs(dtmp_c[0][1][4][1][0][0] - (0.024568482604823683)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][4][1][0][1], 0.08749273467036488); fflush(stdout);
  assert( fabs(dtmp_c[0][1][4][1][0][1] - (0.08749273467036488)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][4][1][0][2], -0.039603319111111486); fflush(stdout);
  assert( fabs(dtmp_c[0][1][4][1][0][2] - (-0.039603319111111486)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][4][1][1][0], -0.053528296761971975); fflush(stdout);
  assert( fabs(dtmp_c[0][1][4][1][1][0] - (-0.053528296761971975)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][4][1][1][1], 0.025573742786162115); fflush(stdout);
  assert( fabs(dtmp_c[0][1][4][1][1][1] - (0.025573742786162115)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][4][1][1][2], 0.01585516278079176); fflush(stdout);
  assert( fabs(dtmp_c[0][1][4][1][1][2] - (0.01585516278079176)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][4][1][2][0], 0.23388256744788105); fflush(stdout);
  assert( fabs(dtmp_c[0][1][4][1][2][0] - (0.23388256744788105)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][4][1][2][1], 0.27290602141283044); fflush(stdout);
  assert( fabs(dtmp_c[0][1][4][1][2][1] - (0.27290602141283044)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][4][1][2][2], 0.08400594794954518); fflush(stdout);
  assert( fabs(dtmp_c[0][1][4][1][2][2] - (0.08400594794954518)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][4][1][3][0], -0.07721034468420171); fflush(stdout);
  assert( fabs(dtmp_c[0][1][4][1][3][0] - (-0.07721034468420171)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][4][1][3][1], -0.12837380227095865); fflush(stdout);
  assert( fabs(dtmp_c[0][1][4][1][3][1] - (-0.12837380227095865)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][4][1][3][2], 0.0074019690684439185); fflush(stdout);
  assert( fabs(dtmp_c[0][1][4][1][3][2] - (0.0074019690684439185)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][5][0][0][0], 0.028154251977039687); fflush(stdout);
  assert( fabs(dtmp_c[0][1][5][0][0][0] - (0.028154251977039687)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][5][0][0][1], 0.2675648464894296); fflush(stdout);
  assert( fabs(dtmp_c[0][1][5][0][0][1] - (0.2675648464894296)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][5][0][0][2], -0.12452862415334812); fflush(stdout);
  assert( fabs(dtmp_c[0][1][5][0][0][2] - (-0.12452862415334812)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][5][0][1][0], -0.32448339329853426); fflush(stdout);
  assert( fabs(dtmp_c[0][1][5][0][1][0] - (-0.32448339329853426)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][5][0][1][1], 0.05739916372169026); fflush(stdout);
  assert( fabs(dtmp_c[0][1][5][0][1][1] - (0.05739916372169026)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][5][0][1][2], 0.042842019446247125); fflush(stdout);
  assert( fabs(dtmp_c[0][1][5][0][1][2] - (0.042842019446247125)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][5][0][2][0], 0.009582449802574257); fflush(stdout);
  assert( fabs(dtmp_c[0][1][5][0][2][0] - (0.009582449802574257)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][5][0][2][1], -0.2103200962586508); fflush(stdout);
  assert( fabs(dtmp_c[0][1][5][0][2][1] - (-0.2103200962586508)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][5][0][2][2], 0.09232438690101506); fflush(stdout);
  assert( fabs(dtmp_c[0][1][5][0][2][2] - (0.09232438690101506)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][5][0][3][0], -0.8499383650671706); fflush(stdout);
  assert( fabs(dtmp_c[0][1][5][0][3][0] - (-0.8499383650671706)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][5][0][3][1], -0.9945478199910702); fflush(stdout);
  assert( fabs(dtmp_c[0][1][5][0][3][1] - (-0.9945478199910702)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][5][0][3][2], -0.046223814939612445); fflush(stdout);
  assert( fabs(dtmp_c[0][1][5][0][3][2] - (-0.046223814939612445)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][5][1][0][0], 0.021674237430510352); fflush(stdout);
  assert( fabs(dtmp_c[0][1][5][1][0][0] - (0.021674237430510352)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][5][1][0][1], 0.07619490347632643); fflush(stdout);
  assert( fabs(dtmp_c[0][1][5][1][0][1] - (0.07619490347632643)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][5][1][0][2], -0.03437119170244582); fflush(stdout);
  assert( fabs(dtmp_c[0][1][5][1][0][2] - (-0.03437119170244582)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][5][1][1][0], -0.04474232067042888); fflush(stdout);
  assert( fabs(dtmp_c[0][1][5][1][1][0] - (-0.04474232067042888)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][5][1][1][1], 0.022605073588539335); fflush(stdout);
  assert( fabs(dtmp_c[0][1][5][1][1][1] - (0.022605073588539335)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][5][1][1][2], 0.013736156923521098); fflush(stdout);
  assert( fabs(dtmp_c[0][1][5][1][1][2] - (0.013736156923521098)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][5][1][2][0], 0.2102381984988093); fflush(stdout);
  assert( fabs(dtmp_c[0][1][5][1][2][0] - (0.2102381984988093)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][5][1][2][1], 0.24811662611440918); fflush(stdout);
  assert( fabs(dtmp_c[0][1][5][1][2][1] - (0.24811662611440918)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][5][1][2][2], 0.07387138323364789); fflush(stdout);
  assert( fabs(dtmp_c[0][1][5][1][2][2] - (0.07387138323364789)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][5][1][3][0], -0.061616378298223); fflush(stdout);
  assert( fabs(dtmp_c[0][1][5][1][3][0] - (-0.061616378298223)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][5][1][3][1], -0.10988792547317254); fflush(stdout);
  assert( fabs(dtmp_c[0][1][5][1][3][1] - (-0.10988792547317254)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][1][5][1][3][2], 0.007206101623458766); fflush(stdout);
  assert( fabs(dtmp_c[0][1][5][1][3][2] - (0.007206101623458766)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][0][0][0][0], 0.19018068428042578); fflush(stdout);
  assert( fabs(dtmp_c[0][2][0][0][0][0] - (0.19018068428042578)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][0][0][0][1], 1.222540367818457); fflush(stdout);
  assert( fabs(dtmp_c[0][2][0][0][0][1] - (1.222540367818457)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][0][0][0][2], -0.28452524018948755); fflush(stdout);
  assert( fabs(dtmp_c[0][2][0][0][0][2] - (-0.28452524018948755)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][0][0][1][0], -1.513804305499182); fflush(stdout);
  assert( fabs(dtmp_c[0][2][0][0][1][0] - (-1.513804305499182)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][0][0][1][1], 0.44594827720555874); fflush(stdout);
  assert( fabs(dtmp_c[0][2][0][0][1][1] - (0.44594827720555874)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][0][0][1][2], 0.14594525681695206); fflush(stdout);
  assert( fabs(dtmp_c[0][2][0][0][1][2] - (0.14594525681695206)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][0][0][2][0], 0.14000905291281396); fflush(stdout);
  assert( fabs(dtmp_c[0][2][0][0][2][0] - (0.14000905291281396)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][0][0][2][1], -0.7037590320312661); fflush(stdout);
  assert( fabs(dtmp_c[0][2][0][0][2][1] - (-0.7037590320312661)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][0][0][2][2], 0.3198894225004895); fflush(stdout);
  assert( fabs(dtmp_c[0][2][0][0][2][2] - (0.3198894225004895)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][0][0][3][0], -2.715417700311989); fflush(stdout);
  assert( fabs(dtmp_c[0][2][0][0][3][0] - (-2.715417700311989)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][0][0][3][1], -2.7241785813620303); fflush(stdout);
  assert( fabs(dtmp_c[0][2][0][0][3][1] - (-2.7241785813620303)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][0][0][3][2], 0.11019753358894852); fflush(stdout);
  assert( fabs(dtmp_c[0][2][0][0][3][2] - (0.11019753358894852)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][0][1][0][0], 0.19018068428042578); fflush(stdout);
  assert( fabs(dtmp_c[0][2][0][1][0][0] - (0.19018068428042578)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][0][1][0][1], 1.222540367818457); fflush(stdout);
  assert( fabs(dtmp_c[0][2][0][1][0][1] - (1.222540367818457)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][0][1][0][2], -0.28452524018948755); fflush(stdout);
  assert( fabs(dtmp_c[0][2][0][1][0][2] - (-0.28452524018948755)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][0][1][1][0], -1.513804305499182); fflush(stdout);
  assert( fabs(dtmp_c[0][2][0][1][1][0] - (-1.513804305499182)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][0][1][1][1], 0.44594827720555874); fflush(stdout);
  assert( fabs(dtmp_c[0][2][0][1][1][1] - (0.44594827720555874)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][0][1][1][2], 0.14594525681695206); fflush(stdout);
  assert( fabs(dtmp_c[0][2][0][1][1][2] - (0.14594525681695206)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][0][1][2][0], 0.14000905291281396); fflush(stdout);
  assert( fabs(dtmp_c[0][2][0][1][2][0] - (0.14000905291281396)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][0][1][2][1], -0.7037590320312661); fflush(stdout);
  assert( fabs(dtmp_c[0][2][0][1][2][1] - (-0.7037590320312661)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][0][1][2][2], 0.3198894225004895); fflush(stdout);
  assert( fabs(dtmp_c[0][2][0][1][2][2] - (0.3198894225004895)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][0][1][3][0], -2.715417700311989); fflush(stdout);
  assert( fabs(dtmp_c[0][2][0][1][3][0] - (-2.715417700311989)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][0][1][3][1], -2.7241785813620303); fflush(stdout);
  assert( fabs(dtmp_c[0][2][0][1][3][1] - (-2.7241785813620303)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][0][1][3][2], 0.11019753358894852); fflush(stdout);
  assert( fabs(dtmp_c[0][2][0][1][3][2] - (0.11019753358894852)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][1][0][0][0], 0.08893332816616771); fflush(stdout);
  assert( fabs(dtmp_c[0][2][1][0][0][0] - (0.08893332816616771)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][1][0][0][1], 0.8774402217805611); fflush(stdout);
  assert( fabs(dtmp_c[0][2][1][0][0][1] - (0.8774402217805611)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][1][0][0][2], -0.19640778735660136); fflush(stdout);
  assert( fabs(dtmp_c[0][2][1][0][0][2] - (-0.19640778735660136)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][1][0][1][0], -1.1166921152182026); fflush(stdout);
  assert( fabs(dtmp_c[0][2][1][0][1][0] - (-1.1166921152182026)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][1][0][1][1], 0.3054593920576255); fflush(stdout);
  assert( fabs(dtmp_c[0][2][1][0][1][1] - (0.3054593920576255)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][1][0][1][2], 0.09189401684670058); fflush(stdout);
  assert( fabs(dtmp_c[0][2][1][0][1][2] - (0.09189401684670058)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][1][0][2][0], 0.026769237816560234); fflush(stdout);
  assert( fabs(dtmp_c[0][2][1][0][2][0] - (0.026769237816560234)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][1][0][2][1], -0.6654251261931808); fflush(stdout);
  assert( fabs(dtmp_c[0][2][1][0][2][1] - (-0.6654251261931808)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][1][0][2][2], 0.19314741857168571); fflush(stdout);
  assert( fabs(dtmp_c[0][2][1][0][2][2] - (0.19314741857168571)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][1][0][3][0], -2.1311374340332185); fflush(stdout);
  assert( fabs(dtmp_c[0][2][1][0][3][0] - (-2.1311374340332185)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][1][0][3][1], -2.254921033503945); fflush(stdout);
  assert( fabs(dtmp_c[0][2][1][0][3][1] - (-2.254921033503945)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][1][0][3][2], 0.07301069545320776); fflush(stdout);
  assert( fabs(dtmp_c[0][2][1][0][3][2] - (0.07301069545320776)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][1][1][0][0], 0.04307749857261751); fflush(stdout);
  assert( fabs(dtmp_c[0][2][1][1][0][0] - (0.04307749857261751)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][1][1][0][1], 0.3561328544809898); fflush(stdout);
  assert( fabs(dtmp_c[0][2][1][1][0][1] - (0.3561328544809898)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][1][1][0][2], -0.07848527878423692); fflush(stdout);
  assert( fabs(dtmp_c[0][2][1][1][0][2] - (-0.07848527878423692)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][1][1][1][0], -0.41380706659187266); fflush(stdout);
  assert( fabs(dtmp_c[0][2][1][1][1][0] - (-0.41380706659187266)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][1][1][1][1], 0.12449267477900366); fflush(stdout);
  assert( fabs(dtmp_c[0][2][1][1][1][1] - (0.12449267477900366)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][1][1][1][2], 0.03878774252412468); fflush(stdout);
  assert( fabs(dtmp_c[0][2][1][1][1][2] - (0.03878774252412468)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][1][1][2][0], 0.17995353988724855); fflush(stdout);
  assert( fabs(dtmp_c[0][2][1][1][2][0] - (0.17995353988724855)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][1][1][2][1], 0.0006508635106195759); fflush(stdout);
  assert( fabs(dtmp_c[0][2][1][1][2][1] - (0.0006508635106195759)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][1][1][2][2], 0.09623762360751305); fflush(stdout);
  assert( fabs(dtmp_c[0][2][1][1][2][2] - (0.09623762360751305)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][1][1][3][0], -0.6561461758424563); fflush(stdout);
  assert( fabs(dtmp_c[0][2][1][1][3][0] - (-0.6561461758424563)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][1][1][3][1], -0.6639557387454029); fflush(stdout);
  assert( fabs(dtmp_c[0][2][1][1][3][1] - (-0.6639557387454029)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][1][1][3][2], 0.04092541037956077); fflush(stdout);
  assert( fabs(dtmp_c[0][2][1][1][3][2] - (0.04092541037956077)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][2][0][0][0], 0.05464801083848009); fflush(stdout);
  assert( fabs(dtmp_c[0][2][2][0][0][0] - (0.05464801083848009)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][2][0][0][1], 0.6763613229673414); fflush(stdout);
  assert( fabs(dtmp_c[0][2][2][0][0][1] - (0.6763613229673414)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][2][0][0][2], -0.14168988175736474); fflush(stdout);
  assert( fabs(dtmp_c[0][2][2][0][0][2] - (-0.14168988175736474)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][2][0][1][0], -0.8547780477843795); fflush(stdout);
  assert( fabs(dtmp_c[0][2][2][0][1][0] - (-0.8547780477843795)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][2][0][1][1], 0.2109154904123181); fflush(stdout);
  assert( fabs(dtmp_c[0][2][2][0][1][1] - (0.2109154904123181)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][2][0][1][2], 0.06071117694424413); fflush(stdout);
  assert( fabs(dtmp_c[0][2][2][0][1][2] - (0.06071117694424413)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][2][0][2][0], 0.0030694260829817206); fflush(stdout);
  assert( fabs(dtmp_c[0][2][2][0][2][0] - (0.0030694260829817206)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][2][0][2][1], -0.5500599726734945); fflush(stdout);
  assert( fabs(dtmp_c[0][2][2][0][2][1] - (-0.5500599726734945)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][2][0][2][2], 0.1260341320676603); fflush(stdout);
  assert( fabs(dtmp_c[0][2][2][0][2][2] - (0.1260341320676603)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][2][0][3][0], -1.668212060873375); fflush(stdout);
  assert( fabs(dtmp_c[0][2][2][0][3][0] - (-1.668212060873375)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][2][0][3][1], -1.8556693465648701); fflush(stdout);
  assert( fabs(dtmp_c[0][2][2][0][3][1] - (-1.8556693465648701)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][2][0][3][2], 0.053416544880742725); fflush(stdout);
  assert( fabs(dtmp_c[0][2][2][0][3][2] - (0.053416544880742725)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][2][1][0][0], 0.017830656847428992); fflush(stdout);
  assert( fabs(dtmp_c[0][2][2][1][0][0] - (0.017830656847428992)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][2][1][0][1], 0.1349488372147283); fflush(stdout);
  assert( fabs(dtmp_c[0][2][2][1][0][1] - (0.1349488372147283)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][2][1][0][2], -0.026174261381932806); fflush(stdout);
  assert( fabs(dtmp_c[0][2][2][1][0][2] - (-0.026174261381932806)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][2][1][1][0], -0.12999345215341282); fflush(stdout);
  assert( fabs(dtmp_c[0][2][2][1][1][0] - (-0.12999345215341282)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][2][1][1][1], 0.043704069259968836); fflush(stdout);
  assert( fabs(dtmp_c[0][2][2][1][1][1] - (0.043704069259968836)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][2][1][1][2], 0.012118690017597503); fflush(stdout);
  assert( fabs(dtmp_c[0][2][2][1][1][2] - (0.012118690017597503)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][2][1][2][0], 0.15424255856794733); fflush(stdout);
  assert( fabs(dtmp_c[0][2][2][1][2][0] - (0.15424255856794733)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][2][1][2][1], 0.14837264849951884); fflush(stdout);
  assert( fabs(dtmp_c[0][2][2][1][2][1] - (0.14837264849951884)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][2][1][2][2], 0.03816655598158332); fflush(stdout);
  assert( fabs(dtmp_c[0][2][2][1][2][2] - (0.03816655598158332)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][2][1][3][0], -0.1363685370803286); fflush(stdout);
  assert( fabs(dtmp_c[0][2][2][1][3][0] - (-0.1363685370803286)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][2][1][3][1], -0.1408731927292388); fflush(stdout);
  assert( fabs(dtmp_c[0][2][2][1][3][1] - (-0.1408731927292388)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][2][1][3][2], 0.02131644828832673); fflush(stdout);
  assert( fabs(dtmp_c[0][2][2][1][3][2] - (0.02131644828832673)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][3][0][0][0], 0.04262341679173909); fflush(stdout);
  assert( fabs(dtmp_c[0][2][3][0][0][0] - (0.04262341679173909)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][3][0][0][1], 0.5415799421124924); fflush(stdout);
  assert( fabs(dtmp_c[0][2][3][0][0][1] - (0.5415799421124924)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][3][0][0][2], -0.1054949091533634); fflush(stdout);
  assert( fabs(dtmp_c[0][2][3][0][0][2] - (-0.1054949091533634)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][3][0][1][0], -0.6698441905843271); fflush(stdout);
  assert( fabs(dtmp_c[0][2][3][0][1][0] - (-0.6698441905843271)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][3][0][1][1], 0.14608102999031564); fflush(stdout);
  assert( fabs(dtmp_c[0][2][3][0][1][1] - (0.14608102999031564)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][3][0][1][2], 0.0418104270030757); fflush(stdout);
  assert( fabs(dtmp_c[0][2][3][0][1][2] - (0.0418104270030757)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][3][0][2][0], 0.0013943671279129054); fflush(stdout);
  assert( fabs(dtmp_c[0][2][3][0][2][0] - (0.0013943671279129054)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][3][0][2][1], -0.44172376712831907); fflush(stdout);
  assert( fabs(dtmp_c[0][2][3][0][2][1] - (-0.44172376712831907)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][3][0][2][2], 0.0871726173290942); fflush(stdout);
  assert( fabs(dtmp_c[0][2][3][0][2][2] - (0.0871726173290942)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][3][0][3][0], -1.3168995886109034); fflush(stdout);
  assert( fabs(dtmp_c[0][2][3][0][3][0] - (-1.3168995886109034)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][3][0][3][1], -1.5356654468826247); fflush(stdout);
  assert( fabs(dtmp_c[0][2][3][0][3][1] - (-1.5356654468826247)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][3][0][3][2], 0.04072938212195742); fflush(stdout);
  assert( fabs(dtmp_c[0][2][3][0][3][2] - (0.04072938212195742)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][3][1][0][0], 0.012960729895566838); fflush(stdout);
  assert( fabs(dtmp_c[0][2][3][1][0][0] - (0.012960729895566838)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][3][1][0][1], 0.0741298097662961); fflush(stdout);
  assert( fabs(dtmp_c[0][2][3][1][0][1] - (0.0741298097662961)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][3][1][0][2], -0.012277061853367108); fflush(stdout);
  assert( fabs(dtmp_c[0][2][3][1][0][2] - (-0.012277061853367108)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][3][1][1][0], -0.05329329218896755); fflush(stdout);
  assert( fabs(dtmp_c[0][2][3][1][1][0] - (-0.05329329218896755)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][3][1][1][1], 0.022702621978253158); fflush(stdout);
  assert( fabs(dtmp_c[0][2][3][1][1][1] - (0.022702621978253158)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][3][1][1][2], 0.005258652722025287); fflush(stdout);
  assert( fabs(dtmp_c[0][2][3][1][1][2] - (0.005258652722025287)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][3][1][2][0], 0.13218601527793972); fflush(stdout);
  assert( fabs(dtmp_c[0][2][3][1][2][0] - (0.13218601527793972)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][3][1][2][1], 0.1691043488583484); fflush(stdout);
  assert( fabs(dtmp_c[0][2][3][1][2][1] - (0.1691043488583484)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][3][1][2][2], 0.021916906999824466); fflush(stdout);
  assert( fabs(dtmp_c[0][2][3][1][2][2] - (0.021916906999824466)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][3][1][3][0], -0.005512801061502169); fflush(stdout);
  assert( fabs(dtmp_c[0][2][3][1][3][0] - (-0.005512801061502169)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][3][1][3][1], -0.010076887579249677); fflush(stdout);
  assert( fabs(dtmp_c[0][2][3][1][3][1] - (-0.010076887579249677)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][3][1][3][2], 0.014930281772348983); fflush(stdout);
  assert( fabs(dtmp_c[0][2][3][1][3][2] - (0.014930281772348983)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][4][0][0][0], 0.03813233269772307); fflush(stdout);
  assert( fabs(dtmp_c[0][2][4][0][0][0] - (0.03813233269772307)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][4][0][0][1], 0.4438604831301593); fflush(stdout);
  assert( fabs(dtmp_c[0][2][4][0][0][1] - (0.4438604831301593)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][4][0][0][2], -0.0805430449279358); fflush(stdout);
  assert( fabs(dtmp_c[0][2][4][0][0][2] - (-0.0805430449279358)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][4][0][1][0], -0.5343857508081257); fflush(stdout);
  assert( fabs(dtmp_c[0][2][4][0][1][0] - (-0.5343857508081257)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][4][0][1][1], 0.10116100660240025); fflush(stdout);
  assert( fabs(dtmp_c[0][2][4][0][1][1] - (0.10116100660240025)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][4][0][1][2], 0.029871464901813947); fflush(stdout);
  assert( fabs(dtmp_c[0][2][4][0][1][2] - (0.029871464901813947)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][4][0][2][0], 0.004458479216758887); fflush(stdout);
  assert( fabs(dtmp_c[0][2][4][0][2][0] - (0.004458479216758887)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][4][0][2][1], -0.3542531456618736); fflush(stdout);
  assert( fabs(dtmp_c[0][2][4][0][2][1] - (-0.3542531456618736)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][4][0][2][2], 0.06322261252403447); fflush(stdout);
  assert( fabs(dtmp_c[0][2][4][0][2][2] - (0.06322261252403447)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][4][0][3][0], -1.0516226855570123); fflush(stdout);
  assert( fabs(dtmp_c[0][2][4][0][3][0] - (-1.0516226855570123)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][4][0][3][1], -1.281200838550817); fflush(stdout);
  assert( fabs(dtmp_c[0][2][4][0][3][1] - (-1.281200838550817)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][4][0][3][2], 0.03172716648865696); fflush(stdout);
  assert( fabs(dtmp_c[0][2][4][0][3][2] - (0.03172716648865696)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][4][1][0][0], 0.011449471373215819); fflush(stdout);
  assert( fabs(dtmp_c[0][2][4][1][0][0] - (0.011449471373215819)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][4][1][0][1], 0.05479041428299073); fflush(stdout);
  assert( fabs(dtmp_c[0][2][4][1][0][1] - (0.05479041428299073)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][4][1][0][2], -0.008194665667702015); fflush(stdout);
  assert( fabs(dtmp_c[0][2][4][1][0][2] - (-0.008194665667702015)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][4][1][1][0], -0.03105543800362735); fflush(stdout);
  assert( fabs(dtmp_c[0][2][4][1][1][0] - (-0.03105543800362735)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][4][1][1][1], 0.01652184035897789); fflush(stdout);
  assert( fabs(dtmp_c[0][2][4][1][1][1] - (0.01652184035897789)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][4][1][1][2], 0.0033403740697770428); fflush(stdout);
  assert( fabs(dtmp_c[0][2][4][1][1][2] - (0.0033403740697770428)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][4][1][2][0], 0.11681595977083911); fflush(stdout);
  assert( fabs(dtmp_c[0][2][4][1][2][0] - (0.11681595977083911)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][4][1][2][1], 0.16198600420647907); fflush(stdout);
  assert( fabs(dtmp_c[0][2][4][1][2][1] - (0.16198600420647907)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][4][1][2][2], 0.016538142923122962); fflush(stdout);
  assert( fabs(dtmp_c[0][2][4][1][2][2] - (0.016538142923122962)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][4][1][3][0], 0.025889040367487987); fflush(stdout);
  assert( fabs(dtmp_c[0][2][4][1][3][0] - (0.025889040367487987)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][4][1][3][1], 0.02066410313780114); fflush(stdout);
  assert( fabs(dtmp_c[0][2][4][1][3][1] - (0.02066410313780114)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][4][1][3][2], 0.012280066471279343); fflush(stdout);
  assert( fabs(dtmp_c[0][2][4][1][3][2] - (0.012280066471279343)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][5][0][0][0], 0.035956095375261614); fflush(stdout);
  assert( fabs(dtmp_c[0][2][5][0][0][0] - (0.035956095375261614)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][5][0][0][1], 0.3696356834752306); fflush(stdout);
  assert( fabs(dtmp_c[0][2][5][0][0][1] - (0.3696356834752306)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][5][0][0][2], -0.06280598690889702); fflush(stdout);
  assert( fabs(dtmp_c[0][2][5][0][0][2] - (-0.06280598690889702)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][5][0][1][0], -0.43276760797237757); fflush(stdout);
  assert( fabs(dtmp_c[0][2][5][0][1][0] - (-0.43276760797237757)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][5][0][1][1], 0.0698904473968921); fflush(stdout);
  assert( fabs(dtmp_c[0][2][5][0][1][1] - (0.0698904473968921)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][5][0][1][2], 0.022044004473798853); fflush(stdout);
  assert( fabs(dtmp_c[0][2][5][0][1][2] - (0.022044004473798853)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][5][0][2][0], 0.007744868429178038); fflush(stdout);
  assert( fabs(dtmp_c[0][2][5][0][2][0] - (0.007744868429178038)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][5][0][2][1], -0.28616781595041874); fflush(stdout);
  assert( fabs(dtmp_c[0][2][5][0][2][1] - (-0.28616781595041874)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][5][0][2][2], 0.04768932393175146); fflush(stdout);
  assert( fabs(dtmp_c[0][2][5][0][2][2] - (0.04768932393175146)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][5][0][3][0], -0.850003834086743); fflush(stdout);
  assert( fabs(dtmp_c[0][2][5][0][3][0] - (-0.850003834086743)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][5][0][3][1], -1.0777041011066946); fflush(stdout);
  assert( fabs(dtmp_c[0][2][5][0][3][1] - (-1.0777041011066946)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][5][0][3][2], 0.025086889680343255); fflush(stdout);
  assert( fabs(dtmp_c[0][2][5][0][3][2] - (0.025086889680343255)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][5][1][0][0], 0.0104860048592525); fflush(stdout);
  assert( fabs(dtmp_c[0][2][5][1][0][0] - (0.0104860048592525)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][5][1][0][1], 0.046763658997173364); fflush(stdout);
  assert( fabs(dtmp_c[0][2][5][1][0][1] - (0.046763658997173364)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][5][1][0][2], -0.006711548322015135); fflush(stdout);
  assert( fabs(dtmp_c[0][2][5][1][0][2] - (-0.006711548322015135)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][5][1][1][0], -0.023657571129070156); fflush(stdout);
  assert( fabs(dtmp_c[0][2][5][1][1][0] - (-0.023657571129070156)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][5][1][1][1], 0.014074513707690191); fflush(stdout);
  assert( fabs(dtmp_c[0][2][5][1][1][1] - (0.014074513707690191)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][5][1][1][2], 0.0026860087114576542); fflush(stdout);
  assert( fabs(dtmp_c[0][2][5][1][1][2] - (0.0026860087114576542)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][5][1][2][0], 0.10548456189403788); fflush(stdout);
  assert( fabs(dtmp_c[0][2][5][1][2][0] - (0.10548456189403788)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][5][1][2][1], 0.15005908607648322); fflush(stdout);
  assert( fabs(dtmp_c[0][2][5][1][2][1] - (0.15005908607648322)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][5][1][2][2], 0.014160777013293286); fflush(stdout);
  assert( fabs(dtmp_c[0][2][5][1][2][2] - (0.014160777013293286)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][5][1][3][0], 0.03174368702359653); fflush(stdout);
  assert( fabs(dtmp_c[0][2][5][1][3][0] - (0.03174368702359653)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][5][1][3][1], 0.026139955230237814); fflush(stdout);
  assert( fabs(dtmp_c[0][2][5][1][3][1] - (0.026139955230237814)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][2][5][1][3][2], 0.010807900710766975); fflush(stdout);
  assert( fabs(dtmp_c[0][2][5][1][3][2] - (0.010807900710766975)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][0][0][0][0], 0.08718192269762033); fflush(stdout);
  assert( fabs(dtmp_c[0][3][0][0][0][0] - (0.08718192269762033)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][0][0][0][1], 1.0956261439767068); fflush(stdout);
  assert( fabs(dtmp_c[0][3][0][0][0][1] - (1.0956261439767068)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][0][0][0][2], -0.11332788446698937); fflush(stdout);
  assert( fabs(dtmp_c[0][3][0][0][0][2] - (-0.11332788446698937)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][0][0][1][0], -1.4398532440806); fflush(stdout);
  assert( fabs(dtmp_c[0][3][0][0][1][0] - (-1.4398532440806)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][0][0][1][1], 0.4039631544446396); fflush(stdout);
  assert( fabs(dtmp_c[0][3][0][0][1][1] - (0.4039631544446396)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][0][0][1][2], 0.06559508352852411); fflush(stdout);
  assert( fabs(dtmp_c[0][3][0][0][1][2] - (0.06559508352852411)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][0][0][2][0], -4.476331029233324e-06); fflush(stdout);
  assert( fabs(dtmp_c[0][3][0][0][2][0] - (-4.476331029233324e-06)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][0][0][2][1], -0.8250727027485139); fflush(stdout);
  assert( fabs(dtmp_c[0][3][0][0][2][1] - (-0.8250727027485139)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][0][0][2][2], 0.140283226639315); fflush(stdout);
  assert( fabs(dtmp_c[0][3][0][0][2][2] - (0.140283226639315)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][0][0][3][0], -1.7014729854655117); fflush(stdout);
  assert( fabs(dtmp_c[0][3][0][0][3][0] - (-1.7014729854655117)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][0][0][3][1], -1.7812853021119628); fflush(stdout);
  assert( fabs(dtmp_c[0][3][0][0][3][1] - (-1.7812853021119628)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][0][0][3][2], 0.1543406487837135); fflush(stdout);
  assert( fabs(dtmp_c[0][3][0][0][3][2] - (0.1543406487837135)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][0][1][0][0], 0.08718192269762033); fflush(stdout);
  assert( fabs(dtmp_c[0][3][0][1][0][0] - (0.08718192269762033)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][0][1][0][1], 1.0956261439767068); fflush(stdout);
  assert( fabs(dtmp_c[0][3][0][1][0][1] - (1.0956261439767068)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][0][1][0][2], -0.11332788446698937); fflush(stdout);
  assert( fabs(dtmp_c[0][3][0][1][0][2] - (-0.11332788446698937)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][0][1][1][0], -1.4398532440806); fflush(stdout);
  assert( fabs(dtmp_c[0][3][0][1][1][0] - (-1.4398532440806)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][0][1][1][1], 0.4039631544446396); fflush(stdout);
  assert( fabs(dtmp_c[0][3][0][1][1][1] - (0.4039631544446396)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][0][1][1][2], 0.06559508352852411); fflush(stdout);
  assert( fabs(dtmp_c[0][3][0][1][1][2] - (0.06559508352852411)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][0][1][2][0], -4.476331029233324e-06); fflush(stdout);
  assert( fabs(dtmp_c[0][3][0][1][2][0] - (-4.476331029233324e-06)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][0][1][2][1], -0.8250727027485139); fflush(stdout);
  assert( fabs(dtmp_c[0][3][0][1][2][1] - (-0.8250727027485139)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][0][1][2][2], 0.140283226639315); fflush(stdout);
  assert( fabs(dtmp_c[0][3][0][1][2][2] - (0.140283226639315)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][0][1][3][0], -1.7014729854655117); fflush(stdout);
  assert( fabs(dtmp_c[0][3][0][1][3][0] - (-1.7014729854655117)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][0][1][3][1], -1.7812853021119628); fflush(stdout);
  assert( fabs(dtmp_c[0][3][0][1][3][1] - (-1.7812853021119628)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][0][1][3][2], 0.1543406487837135); fflush(stdout);
  assert( fabs(dtmp_c[0][3][0][1][3][2] - (0.1543406487837135)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][1][0][0][0], 0.05043422858493376); fflush(stdout);
  assert( fabs(dtmp_c[0][3][1][0][0][0] - (0.05043422858493376)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][1][0][0][1], 0.8436456791166707); fflush(stdout);
  assert( fabs(dtmp_c[0][3][1][0][0][1] - (0.8436456791166707)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][1][0][0][2], -0.07782284177832925); fflush(stdout);
  assert( fabs(dtmp_c[0][3][1][0][0][2] - (-0.07782284177832925)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][1][0][1][0], -1.0908830725680714); fflush(stdout);
  assert( fabs(dtmp_c[0][3][1][0][1][0] - (-1.0908830725680714)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][1][0][1][1], 0.27849232341448016); fflush(stdout);
  assert( fabs(dtmp_c[0][3][1][0][1][1] - (0.27849232341448016)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][1][0][1][2], 0.04032068104246521); fflush(stdout);
  assert( fabs(dtmp_c[0][3][1][0][1][2] - (0.04032068104246521)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][1][0][2][0], -0.022789922749039605); fflush(stdout);
  assert( fabs(dtmp_c[0][3][1][0][2][0] - (-0.022789922749039605)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][1][0][2][1], -0.6935584081865638); fflush(stdout);
  assert( fabs(dtmp_c[0][3][1][0][2][1] - (-0.6935584081865638)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][1][0][2][2], 0.08478779894879457); fflush(stdout);
  assert( fabs(dtmp_c[0][3][1][0][2][2] - (0.08478779894879457)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][1][0][3][0], -1.3723379299400016); fflush(stdout);
  assert( fabs(dtmp_c[0][3][1][0][3][0] - (-1.3723379299400016)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][1][0][3][1], -1.5548428450064795); fflush(stdout);
  assert( fabs(dtmp_c[0][3][1][0][3][1] - (-1.5548428450064795)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][1][0][3][2], 0.10130465952185594); fflush(stdout);
  assert( fabs(dtmp_c[0][3][1][0][3][2] - (0.10130465952185594)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][1][1][0][0], 0.014719283290584119); fflush(stdout);
  assert( fabs(dtmp_c[0][3][1][1][0][0] - (0.014719283290584119)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][1][1][0][1], 0.30163478854688897); fflush(stdout);
  assert( fabs(dtmp_c[0][3][1][1][0][1] - (0.30163478854688897)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][1][1][0][2], -0.02909671784388913); fflush(stdout);
  assert( fabs(dtmp_c[0][3][1][1][0][2] - (-0.02909671784388913)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][1][1][1][0], -0.38112709248079696); fflush(stdout);
  assert( fabs(dtmp_c[0][3][1][1][1][0] - (-0.38112709248079696)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][1][1][1][1], 0.10578546540659697); fflush(stdout);
  assert( fabs(dtmp_c[0][3][1][1][1][1] - (0.10578546540659697)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][1][1][1][2], 0.01645191921302484); fflush(stdout);
  assert( fabs(dtmp_c[0][3][1][1][1][2] - (0.01645191921302484)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][1][1][2][0], 0.06808149228834358); fflush(stdout);
  assert( fabs(dtmp_c[0][3][1][1][2][0] - (0.06808149228834358)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][1][1][2][1], -0.1267575671578415); fflush(stdout);
  assert( fabs(dtmp_c[0][3][1][1][2][1] - (-0.1267575671578415)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][1][1][2][2], 0.03699882213241174); fflush(stdout);
  assert( fabs(dtmp_c[0][3][1][1][2][2] - (0.03699882213241174)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][1][1][3][0], -0.3901176612613575); fflush(stdout);
  assert( fabs(dtmp_c[0][3][1][1][3][0] - (-0.3901176612613575)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][1][1][3][1], -0.4069815362496829); fflush(stdout);
  assert( fabs(dtmp_c[0][3][1][1][3][1] - (-0.4069815362496829)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][1][1][3][2], 0.04168307889556077); fflush(stdout);
  assert( fabs(dtmp_c[0][3][1][1][3][2] - (0.04168307889556077)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][2][0][0][0], 0.039670887881486795); fflush(stdout);
  assert( fabs(dtmp_c[0][3][2][0][0][0] - (0.039670887881486795)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][2][0][0][1], 0.6752946988413608); fflush(stdout);
  assert( fabs(dtmp_c[0][3][2][0][0][1] - (0.6752946988413608)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][2][0][0][2], -0.05556387822692927); fflush(stdout);
  assert( fabs(dtmp_c[0][3][2][0][0][2] - (-0.05556387822692927)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][2][0][1][0], -0.8455781340511063); fflush(stdout);
  assert( fabs(dtmp_c[0][3][2][0][1][0] - (-0.8455781340511063)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][2][0][1][1], 0.1929284060004117); fflush(stdout);
  assert( fabs(dtmp_c[0][3][2][0][1][1] - (0.1929284060004117)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][2][0][1][2], 0.02586200786066304); fflush(stdout);
  assert( fabs(dtmp_c[0][3][2][0][1][2] - (0.02586200786066304)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][2][0][2][0], -0.018866959924278344); fflush(stdout);
  assert( fabs(dtmp_c[0][3][2][0][2][0] - (-0.018866959924278344)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][2][0][2][1], -0.5589135055651518); fflush(stdout);
  assert( fabs(dtmp_c[0][3][2][0][2][1] - (-0.5589135055651518)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][2][0][2][2], 0.054111124870279234); fflush(stdout);
  assert( fabs(dtmp_c[0][3][2][0][2][2] - (0.054111124870279234)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][2][0][3][0], -1.0926232575426393); fflush(stdout);
  assert( fabs(dtmp_c[0][3][2][0][3][0] - (-1.0926232575426393)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][2][0][3][1], -1.3250284516498996); fflush(stdout);
  assert( fabs(dtmp_c[0][3][2][0][3][1] - (-1.3250284516498996)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][2][0][3][2], 0.07002039842410261); fflush(stdout);
  assert( fabs(dtmp_c[0][3][2][0][3][2] - (0.07002039842410261)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][2][1][0][0], 0.004534575498258681); fflush(stdout);
  assert( fabs(dtmp_c[0][3][2][1][0][0] - (0.004534575498258681)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][2][1][0][1], 0.0974854919060389); fflush(stdout);
  assert( fabs(dtmp_c[0][3][2][1][0][1] - (0.0974854919060389)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][2][1][0][2], -0.008185499871923137); fflush(stdout);
  assert( fabs(dtmp_c[0][3][2][1][0][2] - (-0.008185499871923137)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][2][1][1][0], -0.10841665521405669); fflush(stdout);
  assert( fabs(dtmp_c[0][3][2][1][1][0] - (-0.10841665521405669)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][2][1][1][1], 0.03201617001409644); fflush(stdout);
  assert( fabs(dtmp_c[0][3][2][1][1][1] - (0.03201617001409644)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][2][1][1][2], 0.00440274244080579); fflush(stdout);
  assert( fabs(dtmp_c[0][3][2][1][1][2] - (0.00440274244080579)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][2][1][2][0], 0.061337434571626805); fflush(stdout);
  assert( fabs(dtmp_c[0][3][2][1][2][0] - (0.061337434571626805)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][2][1][2][1], 0.03481005562425092); fflush(stdout);
  assert( fabs(dtmp_c[0][3][2][1][2][1] - (0.03481005562425092)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][2][1][2][2], 0.011175922738139813); fflush(stdout);
  assert( fabs(dtmp_c[0][3][2][1][2][2] - (0.011175922738139813)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][2][1][3][0], -0.06462038020428579); fflush(stdout);
  assert( fabs(dtmp_c[0][3][2][1][3][0] - (-0.06462038020428579)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][2][1][3][1], -0.05351292689086246); fflush(stdout);
  assert( fabs(dtmp_c[0][3][2][1][3][1] - (-0.05351292689086246)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][2][1][3][2], 0.013296132113609202); fflush(stdout);
  assert( fabs(dtmp_c[0][3][2][1][3][2] - (0.013296132113609202)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][3][0][0][0], 0.037059551434450194); fflush(stdout);
  assert( fabs(dtmp_c[0][3][3][0][0][0] - (0.037059551434450194)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][3][0][0][1], 0.5531176807788016); fflush(stdout);
  assert( fabs(dtmp_c[0][3][3][0][0][1] - (0.5531176807788016)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][3][0][0][2], -0.040990433687054274); fflush(stdout);
  assert( fabs(dtmp_c[0][3][3][0][0][2] - (-0.040990433687054274)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][3][0][1][0], -0.6671535149772022); fflush(stdout);
  assert( fabs(dtmp_c[0][3][3][0][1][0] - (-0.6671535149772022)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][3][0][1][1], 0.13373989270757727); fflush(stdout);
  assert( fabs(dtmp_c[0][3][3][0][1][1] - (0.13373989270757727)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][3][0][1][2], 0.01730921894958104); fflush(stdout);
  assert( fabs(dtmp_c[0][3][3][0][1][2] - (0.01730921894958104)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][3][0][2][0], -0.010287705081322483); fflush(stdout);
  assert( fabs(dtmp_c[0][3][3][0][2][0] - (-0.010287705081322483)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][3][0][2][1], -0.44790260626788075); fflush(stdout);
  assert( fabs(dtmp_c[0][3][3][0][2][1] - (-0.44790260626788075)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][3][0][2][2], 0.03639437038789638); fflush(stdout);
  assert( fabs(dtmp_c[0][3][3][0][2][2] - (0.03639437038789638)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][3][0][3][0], -0.8723378539673872); fflush(stdout);
  assert( fabs(dtmp_c[0][3][3][0][3][0] - (-0.8723378539673872)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][3][0][3][1], -1.1258292468805597); fflush(stdout);
  assert( fabs(dtmp_c[0][3][3][0][3][1] - (-1.1258292468805597)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][3][0][3][2], 0.05044211348038961); fflush(stdout);
  assert( fabs(dtmp_c[0][3][3][0][3][2] - (0.05044211348038961)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][3][1][0][0], 0.003729281832846056); fflush(stdout);
  assert( fabs(dtmp_c[0][3][3][1][0][0] - (0.003729281832846056)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][3][1][0][1], 0.043084899912543786); fflush(stdout);
  assert( fabs(dtmp_c[0][3][3][1][0][1] - (0.043084899912543786)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][3][1][0][2], -0.002869737097021148); fflush(stdout);
  assert( fabs(dtmp_c[0][3][3][1][0][2] - (-0.002869737097021148)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][3][1][1][0], -0.036047934095113754); fflush(stdout);
  assert( fabs(dtmp_c[0][3][3][1][1][0] - (-0.036047934095113754)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][3][1][1][1], 0.013377620032024761); fflush(stdout);
  assert( fabs(dtmp_c[0][3][3][1][1][1] - (0.013377620032024761)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][3][1][1][2], 0.0014069547300744648); fflush(stdout);
  assert( fabs(dtmp_c[0][3][3][1][1][2] - (0.0014069547300744648)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][3][1][2][0], 0.05135253751206031); fflush(stdout);
  assert( fabs(dtmp_c[0][3][3][1][2][0] - (0.05135253751206031)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][3][1][2][1], 0.06845441989957772); fflush(stdout);
  assert( fabs(dtmp_c[0][3][3][1][2][1] - (0.06845441989957772)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][3][1][2][2], 0.0045360749532722365); fflush(stdout);
  assert( fabs(dtmp_c[0][3][3][1][2][2] - (0.0045360749532722365)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][3][1][3][0], 0.015187044606887428); fflush(stdout);
  assert( fabs(dtmp_c[0][3][3][1][3][0] - (0.015187044606887428)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][3][1][3][1], 0.03418602160543856); fflush(stdout);
  assert( fabs(dtmp_c[0][3][3][1][3][1] - (0.03418602160543856)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][3][1][3][2], 0.00586529050396006); fflush(stdout);
  assert( fabs(dtmp_c[0][3][3][1][3][2] - (0.00586529050396006)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][4][0][0][0], 0.036494137154527005); fflush(stdout);
  assert( fabs(dtmp_c[0][3][4][0][0][0] - (0.036494137154527005)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][4][0][0][1], 0.46018838552304103); fflush(stdout);
  assert( fabs(dtmp_c[0][3][4][0][0][1] - (0.46018838552304103)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][4][0][0][2], -0.03109498371004895); fflush(stdout);
  assert( fabs(dtmp_c[0][3][4][0][0][2] - (-0.03109498371004895)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][4][0][1][0], -0.5344679494594589); fflush(stdout);
  assert( fabs(dtmp_c[0][3][4][0][1][0] - (-0.5344679494594589)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][4][0][1][1], 0.09253579359135737); fflush(stdout);
  assert( fabs(dtmp_c[0][3][4][0][1][1] - (0.09253579359135737)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][4][0][1][2], 0.012070708556363306); fflush(stdout);
  assert( fabs(dtmp_c[0][3][4][0][1][2] - (0.012070708556363306)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][4][0][2][0], -0.0025533970525859068); fflush(stdout);
  assert( fabs(dtmp_c[0][3][4][0][2][0] - (-0.0025533970525859068)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][4][0][2][1], -0.36102046492602013); fflush(stdout);
  assert( fabs(dtmp_c[0][3][4][0][2][1] - (-0.36102046492602013)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][4][0][2][2], 0.025726022899259748); fflush(stdout);
  assert( fabs(dtmp_c[0][3][4][0][2][2] - (0.025726022899259748)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][4][0][3][0], -0.7021343287791432); fflush(stdout);
  assert( fabs(dtmp_c[0][3][4][0][3][0] - (-0.7021343287791432)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][4][0][3][1], -0.9595046123161148); fflush(stdout);
  assert( fabs(dtmp_c[0][3][4][0][3][1] - (-0.9595046123161148)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][4][0][3][2], 0.037615120338098156); fflush(stdout);
  assert( fabs(dtmp_c[0][3][4][0][3][2] - (0.037615120338098156)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][4][1][0][0], 0.0038216374362923912); fflush(stdout);
  assert( fabs(dtmp_c[0][3][4][1][0][0] - (0.0038216374362923912)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][4][1][0][1], 0.02745167766942589); fflush(stdout);
  assert( fabs(dtmp_c[0][3][4][1][0][1] - (0.02745167766942589)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][4][1][0][2], -0.0014542722905223392); fflush(stdout);
  assert( fabs(dtmp_c[0][3][4][1][0][2] - (-0.0014542722905223392)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][4][1][1][0], -0.016087458287979937); fflush(stdout);
  assert( fabs(dtmp_c[0][3][4][1][1][0] - (-0.016087458287979937)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][4][1][1][1], 0.008321531332585943); fflush(stdout);
  assert( fabs(dtmp_c[0][3][4][1][1][1] - (0.008321531332585943)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][4][1][1][2], 0.0006381881097886708); fflush(stdout);
  assert( fabs(dtmp_c[0][3][4][1][1][2] - (0.0006381881097886708)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][4][1][2][0], 0.04466456765066179); fflush(stdout);
  assert( fabs(dtmp_c[0][3][4][1][2][0] - (0.04466456765066179)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][4][1][2][1], 0.07177239175469316); fflush(stdout);
  assert( fabs(dtmp_c[0][3][4][1][2][1] - (0.07177239175469316)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][4][1][2][2], 0.0027061007144890755); fflush(stdout);
  assert( fabs(dtmp_c[0][3][4][1][2][2] - (0.0027061007144890755)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][4][1][3][0], 0.033140867046369864); fflush(stdout);
  assert( fabs(dtmp_c[0][3][4][1][3][0] - (0.033140867046369864)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][4][1][3][1], 0.053011110578175744); fflush(stdout);
  assert( fabs(dtmp_c[0][3][4][1][3][1] - (0.053011110578175744)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][4][1][3][2], 0.0037420368994668983); fflush(stdout);
  assert( fabs(dtmp_c[0][3][4][1][3][2] - (0.0037420368994668983)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][5][0][0][0], 0.03593144047367517); fflush(stdout);
  assert( fabs(dtmp_c[0][3][5][0][0][0] - (0.03593144047367517)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][5][0][0][1], 0.38735355395186644); fflush(stdout);
  assert( fabs(dtmp_c[0][3][5][0][0][1] - (0.38735355395186644)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][5][0][0][2], -0.02415810793192205); fflush(stdout);
  assert( fabs(dtmp_c[0][3][5][0][0][2] - (-0.02415810793192205)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][5][0][1][0], -0.4340531088122261); fflush(stdout);
  assert( fabs(dtmp_c[0][3][5][0][1][0] - (-0.4340531088122261)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][5][0][1][1], 0.06378553712352382); fflush(stdout);
  assert( fabs(dtmp_c[0][3][5][0][1][1] - (0.06378553712352382)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][5][0][1][2], 0.008743324805611132); fflush(stdout);
  assert( fabs(dtmp_c[0][3][5][0][1][2] - (0.008743324805611132)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][5][0][2][0], 0.0032607112383044715); fflush(stdout);
  assert( fabs(dtmp_c[0][3][5][0][2][0] - (0.0032607112383044715)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][5][0][2][1], -0.2937110848433893); fflush(stdout);
  assert( fabs(dtmp_c[0][3][5][0][2][1] - (-0.2937110848433893)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][5][0][2][2], 0.019017520508765057); fflush(stdout);
  assert( fabs(dtmp_c[0][3][5][0][2][2] - (0.019017520508765057)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][5][0][3][0], -0.5707320233347631); fflush(stdout);
  assert( fabs(dtmp_c[0][3][5][0][3][0] - (-0.5707320233347631)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][5][0][3][1], -0.8215201875027386); fflush(stdout);
  assert( fabs(dtmp_c[0][3][5][0][3][1] - (-0.8215201875027386)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][5][0][3][2], 0.02887671710464598); fflush(stdout);
  assert( fabs(dtmp_c[0][3][5][0][3][2] - (0.02887671710464598)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][5][1][0][0], 0.0037695308089584303); fflush(stdout);
  assert( fabs(dtmp_c[0][3][5][1][0][0] - (0.0037695308089584303)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][5][1][0][1], 0.02211108580460943); fflush(stdout);
  assert( fabs(dtmp_c[0][3][5][1][0][1] - (0.02211108580460943)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][5][1][0][2], -0.0010341417856208133); fflush(stdout);
  assert( fabs(dtmp_c[0][3][5][1][0][2] - (-0.0010341417856208133)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][5][1][1][0], -0.010183876444822859); fflush(stdout);
  assert( fabs(dtmp_c[0][3][5][1][1][0] - (-0.010183876444822859)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][5][1][1][1], 0.006649504812089057); fflush(stdout);
  assert( fabs(dtmp_c[0][3][5][1][1][1] - (0.006649504812089057)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][5][1][1][2], 0.00042350858782990526); fflush(stdout);
  assert( fabs(dtmp_c[0][3][5][1][1][2] - (0.00042350858782990526)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][5][1][2][0], 0.04014308010355947); fflush(stdout);
  assert( fabs(dtmp_c[0][3][5][1][2][0] - (0.04014308010355947)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][5][1][2][1], 0.0683153468714047); fflush(stdout);
  assert( fabs(dtmp_c[0][3][5][1][2][1] - (0.0683153468714047)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][5][1][2][2], 0.002110124318732031); fflush(stdout);
  assert( fabs(dtmp_c[0][3][5][1][2][2] - (0.002110124318732031)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][5][1][3][0], 0.035510656097454345); fflush(stdout);
  assert( fabs(dtmp_c[0][3][5][1][3][0] - (0.035510656097454345)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][5][1][3][1], 0.05434684019178174); fflush(stdout);
  assert( fabs(dtmp_c[0][3][5][1][3][1] - (0.05434684019178174)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][3][5][1][3][2], 0.0030055325665043265); fflush(stdout);
  assert( fabs(dtmp_c[0][3][5][1][3][2] - (0.0030055325665043265)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][0][0][0][0], 0.03641894886362417); fflush(stdout);
  assert( fabs(dtmp_c[0][4][0][0][0][0] - (0.03641894886362417)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][0][0][0][1], 0.9384662053704101); fflush(stdout);
  assert( fabs(dtmp_c[0][4][0][0][0][1] - (0.9384662053704101)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][0][0][0][2], -0.04284654622509029); fflush(stdout);
  assert( fabs(dtmp_c[0][4][0][0][0][2] - (-0.04284654622509029)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][0][0][1][0], -1.2513195856973107); fflush(stdout);
  assert( fabs(dtmp_c[0][4][0][0][1][0] - (-1.2513195856973107)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][0][0][1][1], 0.32712040135830606); fflush(stdout);
  assert( fabs(dtmp_c[0][4][0][0][1][1] - (0.32712040135830606)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][0][0][1][2], 0.028013604341357287); fflush(stdout);
  assert( fabs(dtmp_c[0][4][0][0][1][2] - (0.028013604341357287)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][0][0][2][0], -0.05899001531572354); fflush(stdout);
  assert( fabs(dtmp_c[0][4][0][0][2][0] - (-0.05899001531572354)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][0][0][2][1], -0.7765388699266639); fflush(stdout);
  assert( fabs(dtmp_c[0][4][0][0][2][1] - (-0.7765388699266639)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][0][0][2][2], 0.06031496080601511); fflush(stdout);
  assert( fabs(dtmp_c[0][4][0][0][2][2] - (0.06031496080601511)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][0][0][3][0], -0.6971546735585968); fflush(stdout);
  assert( fabs(dtmp_c[0][4][0][0][3][0] - (-0.6971546735585968)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][0][0][3][1], -0.8486419020589924); fflush(stdout);
  assert( fabs(dtmp_c[0][4][0][0][3][1] - (-0.8486419020589924)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][0][0][3][2], 0.10926329133859354); fflush(stdout);
  assert( fabs(dtmp_c[0][4][0][0][3][2] - (0.10926329133859354)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][0][1][0][0], 0.03641894886362417); fflush(stdout);
  assert( fabs(dtmp_c[0][4][0][1][0][0] - (0.03641894886362417)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][0][1][0][1], 0.9384662053704101); fflush(stdout);
  assert( fabs(dtmp_c[0][4][0][1][0][1] - (0.9384662053704101)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][0][1][0][2], -0.04284654622509029); fflush(stdout);
  assert( fabs(dtmp_c[0][4][0][1][0][2] - (-0.04284654622509029)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][0][1][1][0], -1.2513195856973107); fflush(stdout);
  assert( fabs(dtmp_c[0][4][0][1][1][0] - (-1.2513195856973107)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][0][1][1][1], 0.32712040135830606); fflush(stdout);
  assert( fabs(dtmp_c[0][4][0][1][1][1] - (0.32712040135830606)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][0][1][1][2], 0.028013604341357287); fflush(stdout);
  assert( fabs(dtmp_c[0][4][0][1][1][2] - (0.028013604341357287)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][0][1][2][0], -0.05899001531572354); fflush(stdout);
  assert( fabs(dtmp_c[0][4][0][1][2][0] - (-0.05899001531572354)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][0][1][2][1], -0.7765388699266639); fflush(stdout);
  assert( fabs(dtmp_c[0][4][0][1][2][1] - (-0.7765388699266639)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][0][1][2][2], 0.06031496080601511); fflush(stdout);
  assert( fabs(dtmp_c[0][4][0][1][2][2] - (0.06031496080601511)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][0][1][3][0], -0.6971546735585968); fflush(stdout);
  assert( fabs(dtmp_c[0][4][0][1][3][0] - (-0.6971546735585968)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][0][1][3][1], -0.8486419020589924); fflush(stdout);
  assert( fabs(dtmp_c[0][4][0][1][3][1] - (-0.8486419020589924)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][0][1][3][2], 0.10926329133859354); fflush(stdout);
  assert( fabs(dtmp_c[0][4][0][1][3][2] - (0.10926329133859354)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][1][0][0][0], 0.029278165901546183); fflush(stdout);
  assert( fabs(dtmp_c[0][4][1][0][0][0] - (0.029278165901546183)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][1][0][0][1], 0.750083864178892); fflush(stdout);
  assert( fabs(dtmp_c[0][4][1][0][0][1] - (0.750083864178892)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][1][0][0][2], -0.028738531591193037); fflush(stdout);
  assert( fabs(dtmp_c[0][4][1][0][0][2] - (-0.028738531591193037)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][1][0][1][0], -0.9596070542896873); fflush(stdout);
  assert( fabs(dtmp_c[0][4][1][0][1][0] - (-0.9596070542896873)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][1][0][1][1], 0.22657669049577392); fflush(stdout);
  assert( fabs(dtmp_c[0][4][1][0][1][1] - (0.22657669049577392)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][1][0][1][2], 0.01666959839009637); fflush(stdout);
  assert( fabs(dtmp_c[0][4][1][0][1][2] - (0.01666959839009637)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][1][0][2][0], -0.046759037568902734); fflush(stdout);
  assert( fabs(dtmp_c[0][4][1][0][2][0] - (-0.046759037568902734)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][1][0][2][1], -0.6288451029825203); fflush(stdout);
  assert( fabs(dtmp_c[0][4][1][0][2][1] - (-0.6288451029825203)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][1][0][2][2], 0.03552513904985523); fflush(stdout);
  assert( fabs(dtmp_c[0][4][1][0][2][2] - (0.03552513904985523)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][1][0][3][0], -0.5850154022860471); fflush(stdout);
  assert( fabs(dtmp_c[0][4][1][0][3][0] - (-0.5850154022860471)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][1][0][3][1], -0.7963270023360557); fflush(stdout);
  assert( fabs(dtmp_c[0][4][1][0][3][1] - (-0.7963270023360557)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][1][0][3][2], 0.0689379911571399); fflush(stdout);
  assert( fabs(dtmp_c[0][4][1][0][3][2] - (0.0689379911571399)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][1][1][0][0], 0.002521820465582915); fflush(stdout);
  assert( fabs(dtmp_c[0][4][1][1][0][0] - (0.002521820465582915)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][1][1][0][1], 0.2510141485885489); fflush(stdout);
  assert( fabs(dtmp_c[0][4][1][1][0][1] - (0.2510141485885489)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][1][1][0][2], -0.010648393450435155); fflush(stdout);
  assert( fabs(dtmp_c[0][4][1][1][0][2] - (-0.010648393450435155)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][1][1][1][0], -0.3260429645891888); fflush(stdout);
  assert( fabs(dtmp_c[0][4][1][1][1][0] - (-0.3260429645891888)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][1][1][1][1], 0.08280266076739266); fflush(stdout);
  assert( fabs(dtmp_c[0][4][1][1][1][1] - (0.08280266076739266)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][1][1][1][2], 0.006876248194932467); fflush(stdout);
  assert( fabs(dtmp_c[0][4][1][1][1][2] - (0.006876248194932467)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][1][1][2][0], 0.021012063662983613); fflush(stdout);
  assert( fabs(dtmp_c[0][4][1][1][2][0] - (0.021012063662983613)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][1][1][2][1], -0.16155855486178486); fflush(stdout);
  assert( fabs(dtmp_c[0][4][1][1][2][1] - (-0.16155855486178486)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][1][1][2][2], 0.015112777954939946); fflush(stdout);
  assert( fabs(dtmp_c[0][4][1][1][2][2] - (0.015112777954939946)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][1][1][3][0], -0.14700983918956803); fflush(stdout);
  assert( fabs(dtmp_c[0][4][1][1][3][0] - (-0.14700983918956803)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][1][1][3][1], -0.190297981115932); fflush(stdout);
  assert( fabs(dtmp_c[0][4][1][1][3][1] - (-0.190297981115932)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][1][1][3][2], 0.02747728751417891); fflush(stdout);
  assert( fabs(dtmp_c[0][4][1][1][3][2] - (0.02747728751417891)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][2][0][0][0], 0.02975057359830538); fflush(stdout);
  assert( fabs(dtmp_c[0][4][2][0][0][0] - (0.02975057359830538)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][2][0][0][1], 0.6136113253237288); fflush(stdout);
  assert( fabs(dtmp_c[0][4][2][0][0][1] - (0.6136113253237288)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][2][0][0][2], -0.020077970533001966); fflush(stdout);
  assert( fabs(dtmp_c[0][4][2][0][0][2] - (-0.020077970533001966)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][2][0][1][0], -0.7487199736878327); fflush(stdout);
  assert( fabs(dtmp_c[0][4][2][0][1][0] - (-0.7487199736878327)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][2][0][1][1], 0.15718402674663834); fflush(stdout);
  assert( fabs(dtmp_c[0][4][2][0][1][1] - (0.15718402674663834)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][2][0][1][2], 0.010316581816643631); fflush(stdout);
  assert( fabs(dtmp_c[0][4][2][0][1][2] - (0.010316581816643631)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][2][0][2][0], -0.031081624443481195); fflush(stdout);
  assert( fabs(dtmp_c[0][4][2][0][2][0] - (-0.031081624443481195)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][2][0][2][1], -0.5037825434259197); fflush(stdout);
  assert( fabs(dtmp_c[0][4][2][0][2][1] - (-0.5037825434259197)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][2][0][2][2], 0.02189024676623996); fflush(stdout);
  assert( fabs(dtmp_c[0][4][2][0][2][2] - (0.02189024676623996)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][2][0][3][0], -0.47763242984655563); fflush(stdout);
  assert( fabs(dtmp_c[0][4][2][0][3][0] - (-0.47763242984655563)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][2][0][3][1], -0.7138719864278856); fflush(stdout);
  assert( fabs(dtmp_c[0][4][2][0][3][1] - (-0.7138719864278856)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][2][0][3][2], 0.04559115586167293); fflush(stdout);
  assert( fabs(dtmp_c[0][4][2][0][3][2] - (0.04559115586167293)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][2][1][0][0], -0.0004113138144010972); fflush(stdout);
  assert( fabs(dtmp_c[0][4][2][1][0][0] - (-0.0004113138144010972)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][2][1][0][1], 0.0732434988140552); fflush(stdout);
  assert( fabs(dtmp_c[0][4][2][1][0][1] - (0.0732434988140552)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][2][1][0][2], -0.0027552586768104124); fflush(stdout);
  assert( fabs(dtmp_c[0][4][2][1][0][2] - (-0.0027552586768104124)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][2][1][1][0], -0.08854462725281863); fflush(stdout);
  assert( fabs(dtmp_c[0][4][2][1][1][0] - (-0.08854462725281863)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][2][1][1][1], 0.022783301706344707); fflush(stdout);
  assert( fabs(dtmp_c[0][4][2][1][1][1] - (0.022783301706344707)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][2][1][1][2], 0.0017278506177386063); fflush(stdout);
  assert( fabs(dtmp_c[0][4][2][1][1][2] - (0.0017278506177386063)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][2][1][2][0], 0.023374437938563114); fflush(stdout);
  assert( fabs(dtmp_c[0][4][2][1][2][0] - (0.023374437938563114)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][2][1][2][1], -0.013699108092899407); fflush(stdout);
  assert( fabs(dtmp_c[0][4][2][1][2][1] - (-0.013699108092899407)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][2][1][2][2], 0.003980007037119749); fflush(stdout);
  assert( fabs(dtmp_c[0][4][2][1][2][2] - (0.003980007037119749)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][2][1][3][0], -0.016343952992919627); fflush(stdout);
  assert( fabs(dtmp_c[0][4][2][1][3][0] - (-0.016343952992919627)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][2][1][3][1], -0.014586500566322031); fflush(stdout);
  assert( fabs(dtmp_c[0][4][2][1][3][1] - (-0.014586500566322031)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][2][1][3][2], 0.007332618693239406); fflush(stdout);
  assert( fabs(dtmp_c[0][4][2][1][3][2] - (0.007332618693239406)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][3][0][0][0], 0.03152720276787019); fflush(stdout);
  assert( fabs(dtmp_c[0][4][3][0][0][0] - (0.03152720276787019)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][3][0][0][1], 0.5098331056486243); fflush(stdout);
  assert( fabs(dtmp_c[0][4][3][0][0][1] - (0.5098331056486243)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][3][0][0][2], -0.014556221738452993); fflush(stdout);
  assert( fabs(dtmp_c[0][4][3][0][0][2] - (-0.014556221738452993)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][3][0][1][0], -0.5931336291794296); fflush(stdout);
  assert( fabs(dtmp_c[0][4][3][0][1][0] - (-0.5931336291794296)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][3][0][1][1], 0.1088947034315803); fflush(stdout);
  assert( fabs(dtmp_c[0][4][3][0][1][1] - (0.1088947034315803)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][3][0][1][2], 0.006666003572591696); fflush(stdout);
  assert( fabs(dtmp_c[0][4][3][0][1][2] - (0.006666003572591696)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][3][0][2][0], -0.01783383785096143); fflush(stdout);
  assert( fabs(dtmp_c[0][4][3][0][2][0] - (-0.01783383785096143)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][3][0][2][1], -0.4052089320950355); fflush(stdout);
  assert( fabs(dtmp_c[0][4][3][0][2][1] - (-0.4052089320950355)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][3][0][2][2], 0.014182917472410976); fflush(stdout);
  assert( fabs(dtmp_c[0][4][3][0][2][2] - (0.014182917472410976)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][3][0][3][0], -0.38796584382291055); fflush(stdout);
  assert( fabs(dtmp_c[0][4][3][0][3][0] - (-0.38796584382291055)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][3][0][3][1], -0.6309135540587967); fflush(stdout);
  assert( fabs(dtmp_c[0][4][3][0][3][1] - (-0.6309135540587967)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][3][0][3][2], 0.031560200751032226); fflush(stdout);
  assert( fabs(dtmp_c[0][4][3][0][3][2] - (0.031560200751032226)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][3][1][0][0], 0.0003491661298714175); fflush(stdout);
  assert( fabs(dtmp_c[0][4][3][1][0][0] - (0.0003491661298714175)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][3][1][0][1], 0.026381431279199764); fflush(stdout);
  assert( fabs(dtmp_c[0][4][3][1][0][1] - (0.026381431279199764)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][3][1][0][2], -0.0007921202792388702); fflush(stdout);
  assert( fabs(dtmp_c[0][4][3][1][0][2] - (-0.0007921202792388702)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][3][1][1][0], -0.026098176834432208); fflush(stdout);
  assert( fabs(dtmp_c[0][4][3][1][1][0] - (-0.026098176834432208)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][3][1][1][1], 0.00784947900724862); fflush(stdout);
  assert( fabs(dtmp_c[0][4][3][1][1][1] - (0.00784947900724862)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][3][1][1][2], 0.00046510849802275064); fflush(stdout);
  assert( fabs(dtmp_c[0][4][3][1][1][2] - (0.00046510849802275064)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][3][1][2][0], 0.018791323424391004); fflush(stdout);
  assert( fabs(dtmp_c[0][4][3][1][2][0] - (0.018791323424391004)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][3][1][2][1], 0.020742693705728024); fflush(stdout);
  assert( fabs(dtmp_c[0][4][3][1][2][1] - (0.020742693705728024)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][3][1][2][2], 0.0012078291929494046); fflush(stdout);
  assert( fabs(dtmp_c[0][4][3][1][2][2] - (0.0012078291929494046)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][3][1][3][0], 0.013896101491367248); fflush(stdout);
  assert( fabs(dtmp_c[0][4][3][1][3][0] - (0.013896101491367248)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][3][1][3][1], 0.03006318719166455); fflush(stdout);
  assert( fabs(dtmp_c[0][4][3][1][3][1] - (0.03006318719166455)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][3][1][3][2], 0.0022955591393589566); fflush(stdout);
  assert( fabs(dtmp_c[0][4][3][1][3][2] - (0.0022955591393589566)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][4][0][0][0], 0.03265523388867102); fflush(stdout);
  assert( fabs(dtmp_c[0][4][4][0][0][0] - (0.03265523388867102)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][4][0][0][1], 0.42853088220057267); fflush(stdout);
  assert( fabs(dtmp_c[0][4][4][0][0][1] - (0.42853088220057267)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][4][0][0][2], -0.010902951422590604); fflush(stdout);
  assert( fabs(dtmp_c[0][4][4][0][0][2] - (-0.010902951422590604)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][4][0][1][0], -0.47647982297023356); fflush(stdout);
  assert( fabs(dtmp_c[0][4][4][0][1][0] - (-0.47647982297023356)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][4][0][1][1], 0.07518799481974753); fflush(stdout);
  assert( fabs(dtmp_c[0][4][4][0][1][1] - (0.07518799481974753)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][4][0][1][2], 0.004503917145645722); fflush(stdout);
  assert( fabs(dtmp_c[0][4][4][0][1][2] - (0.004503917145645722)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][4][0][2][0], -0.007880161843529002); fflush(stdout);
  assert( fabs(dtmp_c[0][4][4][0][2][0] - (-0.007880161843529002)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][4][0][2][1], -0.32878265875584134); fflush(stdout);
  assert( fabs(dtmp_c[0][4][4][0][2][1] - (-0.32878265875584134)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][4][0][2][2], 0.009683690060484331); fflush(stdout);
  assert( fabs(dtmp_c[0][4][4][0][2][2] - (0.009683690060484331)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][4][0][3][0], -0.3161154729336227); fflush(stdout);
  assert( fabs(dtmp_c[0][4][4][0][3][0] - (-0.3161154729336227)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][4][0][3][1], -0.5550924402920809); fflush(stdout);
  assert( fabs(dtmp_c[0][4][4][0][3][1] - (-0.5550924402920809)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][4][0][3][2], 0.022794074554597712); fflush(stdout);
  assert( fabs(dtmp_c[0][4][4][0][3][2] - (0.022794074554597712)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][4][1][0][0], 0.0009273041728807913); fflush(stdout);
  assert( fabs(dtmp_c[0][4][4][1][0][0] - (0.0009273041728807913)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][4][1][0][1], 0.013554790172679776); fflush(stdout);
  assert( fabs(dtmp_c[0][4][4][1][0][1] - (0.013554790172679776)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][4][1][0][2], -0.0002927065283570177); fflush(stdout);
  assert( fabs(dtmp_c[0][4][4][1][0][2] - (-0.0002927065283570177)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][4][1][1][0], -0.009236105302181174); fflush(stdout);
  assert( fabs(dtmp_c[0][4][4][1][1][0] - (-0.009236105302181174)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][4][1][1][1], 0.003986074954054194); fflush(stdout);
  assert( fabs(dtmp_c[0][4][4][1][1][1] - (0.003986074954054194)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][4][1][1][2], 0.0001515591661193652); fflush(stdout);
  assert( fabs(dtmp_c[0][4][4][1][1][2] - (0.0001515591661193652)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][4][1][2][0], 0.015654562520082578); fflush(stdout);
  assert( fabs(dtmp_c[0][4][4][1][2][0] - (0.015654562520082578)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][4][1][2][1], 0.027485394689259427); fflush(stdout);
  assert( fabs(dtmp_c[0][4][4][1][2][1] - (0.027485394689259427)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][4][1][2][2], 0.0005005082686599629); fflush(stdout);
  assert( fabs(dtmp_c[0][4][4][1][2][2] - (0.0005005082686599629)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][4][1][3][0], 0.01995337684220751); fflush(stdout);
  assert( fabs(dtmp_c[0][4][4][1][3][0] - (0.01995337684220751)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][4][1][3][1], 0.03948065183116187); fflush(stdout);
  assert( fabs(dtmp_c[0][4][4][1][3][1] - (0.03948065183116187)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][4][1][3][2], 0.000997506989025163); fflush(stdout);
  assert( fabs(dtmp_c[0][4][4][1][3][2] - (0.000997506989025163)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][5][0][0][0], 0.03272519348705163); fflush(stdout);
  assert( fabs(dtmp_c[0][4][5][0][0][0] - (0.03272519348705163)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][5][0][0][1], 0.3634851380697153); fflush(stdout);
  assert( fabs(dtmp_c[0][4][5][0][0][1] - (0.3634851380697153)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][5][0][0][2], -0.00839865401627028); fflush(stdout);
  assert( fabs(dtmp_c[0][4][5][0][0][2] - (-0.00839865401627028)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][5][0][1][0], -0.38772514005601383); fflush(stdout);
  assert( fabs(dtmp_c[0][4][5][0][1][0] - (-0.38772514005601383)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][5][0][1][1], 0.051648529756526944); fflush(stdout);
  assert( fabs(dtmp_c[0][4][5][0][1][1] - (0.051648529756526944)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][5][0][1][2], 0.0031782902950070464); fflush(stdout);
  assert( fabs(dtmp_c[0][4][5][0][1][2] - (0.0031782902950070464)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][5][0][2][0], -0.0008274677489341686); fflush(stdout);
  assert( fabs(dtmp_c[0][4][5][0][2][0] - (-0.0008274677489341686)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][5][0][2][1], -0.2694431499792432); fflush(stdout);
  assert( fabs(dtmp_c[0][4][5][0][2][1] - (-0.2694431499792432)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][5][0][2][2], 0.006954850054341348); fflush(stdout);
  assert( fabs(dtmp_c[0][4][5][0][2][2] - (0.006954850054341348)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][5][0][3][0], -0.25923091788381186); fflush(stdout);
  assert( fabs(dtmp_c[0][4][5][0][3][0] - (-0.25923091788381186)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][5][0][3][1], -0.48781098678008744); fflush(stdout);
  assert( fabs(dtmp_c[0][4][5][0][3][1] - (-0.48781098678008744)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][5][0][3][2], 0.0170923563068019); fflush(stdout);
  assert( fabs(dtmp_c[0][4][5][0][3][2] - (0.0170923563068019)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][5][1][0][0], 0.0011411734166630787); fflush(stdout);
  assert( fabs(dtmp_c[0][4][5][1][0][0] - (0.0011411734166630787)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][5][1][0][1], 0.009702243023321726); fflush(stdout);
  assert( fabs(dtmp_c[0][4][5][1][0][1] - (0.009702243023321726)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][5][1][0][2], -0.0001594495027813494); fflush(stdout);
  assert( fabs(dtmp_c[0][4][5][1][0][2] - (-0.0001594495027813494)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][5][1][1][0], -0.004505026302608349); fflush(stdout);
  assert( fabs(dtmp_c[0][4][5][1][1][0] - (-0.004505026302608349)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][5][1][1][1], 0.0028614620724252107); fflush(stdout);
  assert( fabs(dtmp_c[0][4][5][1][1][1] - (0.0028614620724252107)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][5][1][1][2], 7.13262143983927e-05); fflush(stdout);
  assert( fabs(dtmp_c[0][4][5][1][1][2] - (7.13262143983927e-05)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][5][1][2][0], 0.013782423535475943); fflush(stdout);
  assert( fabs(dtmp_c[0][4][5][1][2][0] - (0.013782423535475943)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][5][1][2][1], 0.027553325828820363); fflush(stdout);
  assert( fabs(dtmp_c[0][4][5][1][2][1] - (0.027553325828820363)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][5][1][2][2], 0.0003080233937426753); fflush(stdout);
  assert( fabs(dtmp_c[0][4][5][1][2][2] - (0.0003080233937426753)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][5][1][3][0], 0.020226270177046275); fflush(stdout);
  assert( fabs(dtmp_c[0][4][5][1][3][0] - (0.020226270177046275)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][5][1][3][1], 0.03967433210086987); fflush(stdout);
  assert( fabs(dtmp_c[0][4][5][1][3][1] - (0.03967433210086987)) < 1.e-10 );
  printf("%e %e\n", dtmp_c[0][4][5][1][3][2], 0.0006371438697383446); fflush(stdout);
  assert( fabs(dtmp_c[0][4][5][1][3][2] - (0.0006371438697383446)) < 1.e-10 );
}

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

double factor_een[walk_num];
rc = qmckl_get_jastrow_champ_factor_een(context, &(factor_een[0]),walk_num);

assert(fabs(factor_een[0] - (-0.38258026017432123)) < 1e-12);
      

{
  double factor_een_naive[walk_num];
  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;
  assert (ctx != NULL);

  rc = qmckl_compute_jastrow_champ_factor_een_naive(context,
                                                ctx->electron.walker.num,
                                                ctx->electron.num,
                                                ctx->nucleus.num,
                                                ctx->jastrow_champ.cord_num,
                                                ctx->jastrow_champ.dim_c_vector,
                                                ctx->jastrow_champ.c_vector_full,
                                                ctx->jastrow_champ.lkpm_combined_index,
                                                ctx->jastrow_champ.een_rescaled_e,
                                                ctx->jastrow_champ.een_rescaled_n,
                                                factor_een_naive);

  for (int64_t i = 0; i < walk_num; i++) {
    if (fabs(factor_een[i] - factor_een_naive[i]) > 1e-12) {
      printf("factor_een[%ld] = %e\n", i, factor_een[i]);
      printf("factor_een_naive[%ld] = %e\n", i, factor_een_naive[i]);
    }
    assert(fabs(factor_een[i] - factor_een_naive[i]) < 1e-8);
  }
}

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

/*
{
  double factor_een_gl[walk_num][4][elec_num];
  rc = qmckl_check(context,
                   rc = qmckl_get_jastrow_champ_factor_een_gl(context, &(factor_een_gl[0][0][0]),4*walk_num*elec_num)
                   );
  assert(rc == QMCKL_SUCCESS);

  printf("%20.15e\n", factor_een_gl[0][0][0]);
  assert(fabs(8.967809309100624e-02 - factor_een_gl[0][0][0]) < 1e-12);

  printf("%20.15e\n", factor_een_gl[0][1][1]);
  assert(fabs(3.543090132452453e-02 - factor_een_gl[0][1][1]) < 1e-12);

  printf("%20.15e\n", factor_een_gl[0][2][2]);
  assert(fabs(8.996044894431991e-04 - factor_een_gl[0][2][2]) < 1e-12);

  printf("%20.15e\n", factor_een_gl[0][3][3]);
  assert(fabs(-1.175028308456619e+00 - factor_een_gl[0][3][3]) < 1e-12);
}
*/
{
  printf("factor_een_gl_hpc\n");

  qmckl_context_struct* const ctx = (qmckl_context_struct*) context;

  assert (ctx != NULL);
  assert(walk_num == 2);
  assert(elec_num == 10);
  assert(nucl_num == 2);
  assert(cord_num == 5);
  assert(dim_c_vector == 23);

  assert(ctx->electron.walker.num == walk_num);
  assert(ctx->electron.num == elec_num);
  assert(ctx->nucleus.num == nucl_num);
  assert(ctx->jastrow_champ.cord_num == cord_num);
  assert(ctx->jastrow_champ.dim_c_vector == dim_c_vector);

  double factor_een_gl_naive[walk_num*4*elec_num];
  memset(&(factor_een_gl_naive[0]), 0, sizeof(factor_een_gl_naive));

  rc = qmckl_compute_jastrow_champ_factor_een_gl_naive(context,
                                                       ctx->electron.walker.num,
                                                       ctx->electron.num,
                                                       ctx->nucleus.num,
                                                       ctx->jastrow_champ.cord_num,
                                                       ctx->jastrow_champ.dim_c_vector,
                                                       ctx->jastrow_champ.c_vector_full,
                                                       ctx->jastrow_champ.lkpm_combined_index,
                                                       ctx->jastrow_champ.een_rescaled_e,
                                                       ctx->jastrow_champ.een_rescaled_n,
                                                       ctx->jastrow_champ.een_rescaled_e_gl,
                                                       ctx->jastrow_champ.een_rescaled_n_gl,
                                                       factor_een_gl_naive);
  assert(rc == QMCKL_SUCCESS);

  double factor_een_gl_doc[walk_num*4*elec_num];
  memset(&(factor_een_gl_doc[0]), 0, sizeof(factor_een_gl_doc));

  rc = qmckl_compute_jastrow_champ_factor_een_gl_doc(context,
                                                     ctx->electron.walker.num,
                                                     ctx->electron.num,
                                                     ctx->nucleus.num,
                                                     ctx->jastrow_champ.cord_num,
                                                     ctx->jastrow_champ.dim_c_vector,
                                                     ctx->jastrow_champ.c_vector_full,
                                                     ctx->jastrow_champ.lkpm_combined_index,
                                                     ctx->jastrow_champ.tmp_c,
                                                     ctx->jastrow_champ.dtmp_c,
                                                     ctx->jastrow_champ.een_rescaled_n,
                                                     ctx->jastrow_champ.een_rescaled_n_gl,
                                                     factor_een_gl_doc);
  assert(rc == QMCKL_SUCCESS);

  for (int64_t i = 0; i < walk_num*4*elec_num; ++i) {
      if (fabs(factor_een_gl_naive[i] - factor_een_gl_doc[i]) > 1e-12) {
        printf("i = %ld\n", i);
        printf("factor_een_gl_naive = %20.15e\n", factor_een_gl_naive[i]);
        printf("factor_een_gl_doc = %20.15e\n", factor_een_gl_doc[i]);
        fflush(stdout);
      }
      assert(fabs(factor_een_gl_naive[i] - factor_een_gl_doc[i]) < 1e-8);
  }

  double factor_een_gl_hpc[walk_num*4*elec_num];
  memset(&(factor_een_gl_hpc[0]), 0, sizeof(factor_een_gl_hpc));

  rc = qmckl_compute_jastrow_champ_factor_een_gl_hpc(context,
                                                    ctx->electron.walker.num,
                                                    ctx->electron.num,
                                                    ctx->nucleus.num,
                                                    ctx->jastrow_champ.cord_num,
                                                    ctx->jastrow_champ.dim_c_vector,
                                                    ctx->jastrow_champ.c_vector_full,
                                                    ctx->jastrow_champ.lkpm_combined_index,
                                                    ctx->jastrow_champ.tmp_c,
                                                    ctx->jastrow_champ.dtmp_c,
                                                    ctx->jastrow_champ.een_rescaled_n,
                                                    ctx->jastrow_champ.een_rescaled_n_gl,
                                                    factor_een_gl_hpc);

  for (int64_t i = 0; i < walk_num*4*elec_num; ++i) {
      if (fabs(factor_een_gl_doc[i] - factor_een_gl_hpc[i]) > 1e-12) {
        printf("i = %ld\n", i);
        printf("factor_een_gl_doc = %20.15e\n", factor_een_gl_doc[i]);
        printf("factor_een_gl_hpc = %20.15e\n", factor_een_gl_hpc[i]);
        fflush(stdout);
      }
      assert(fabs(factor_een_gl_doc[i] - factor_een_gl_hpc[i]) < 1e-8);
  }
}

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

{
  double factor_een_gl[walk_num][4][elec_num];
  rc = qmckl_get_jastrow_champ_factor_een_gl(context, &(factor_een_gl[0][0][0]),4*walk_num*elec_num);

  double factor_een_grad[walk_num][3][elec_num];
  rc = qmckl_get_jastrow_champ_factor_een_grad(context, &(factor_een_grad[0][0][0]),3*walk_num*elec_num);

  for (int nw=0 ; nw<walk_num ; nw++) {
    for (int k=0 ; k<3; k++) {
      for (int i=0 ; i<elec_num ; i++) {
        assert(fabs(factor_een_gl[nw][k][i] - factor_een_grad[nw][k][i]) < 1e-8);
      }
    }
  }
}

printf("Total Jastrow value\n");
/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

{

  double factor_ee[walk_num];
  rc = qmckl_check(context,
                   qmckl_get_jastrow_champ_factor_ee(context, &(factor_ee[0]), walk_num)
                   );
  assert(rc == QMCKL_SUCCESS);

  double factor_en[walk_num];
  rc = qmckl_check(context,
                   qmckl_get_jastrow_champ_factor_en(context, &(factor_en[0]), walk_num)
                   );
  assert(rc == QMCKL_SUCCESS);

  double factor_een[walk_num];
  rc = qmckl_check(context,
                   qmckl_get_jastrow_champ_factor_een(context, &(factor_een[0]), walk_num)
                   );
  assert(rc == QMCKL_SUCCESS);

  double total_j[walk_num];
  rc = qmckl_check(context,
                   qmckl_get_jastrow_champ_value(context, &(total_j[0]), walk_num)
                   );
  assert(rc == QMCKL_SUCCESS);


  for (int64_t i=0 ; i< walk_num ; ++i) {
    if (fabs(total_j[i] - exp(factor_ee[i] + factor_en[i] + factor_een[i])) > 1e-12) {
      printf("i = %ld\n", i);
      printf("total_j = %20.15e\n", total_j[i]);
      printf("exp(factor_ee + factor_en + factor_een) = %20.15e\n", exp(factor_ee[i] + factor_en[i] + factor_een[i]));
      fflush(stdout);
    }
    assert (fabs(total_j[i] - exp(factor_ee[i] + factor_en[i] + factor_een[i])) < 1.e-8);
  }
}

printf("Total Jastrow derivatives\n");
/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));
{

  double factor_ee_gl[walk_num][4][elec_num];
  rc = qmckl_check(context,
                   qmckl_get_jastrow_champ_factor_ee_gl(context, &(factor_ee_gl[0][0][0]), walk_num*elec_num*4)
                   );
  assert(rc == QMCKL_SUCCESS);

  double factor_en_gl[walk_num][4][elec_num];
  rc = qmckl_check(context,
                   qmckl_get_jastrow_champ_factor_en_gl(context, &(factor_en_gl[0][0][0]), walk_num*elec_num*4)
                   );
  assert(rc == QMCKL_SUCCESS);

  double factor_een_gl[walk_num][4][elec_num];
  rc = qmckl_check(context,
                   qmckl_get_jastrow_champ_factor_een_gl(context, &(factor_een_gl[0][0][0]), walk_num*elec_num*4)
                   );
  assert(rc == QMCKL_SUCCESS);

  double total_j_deriv[walk_num][4][elec_num];
  rc = qmckl_check(context,
                   qmckl_get_jastrow_champ_gl(context, &(total_j_deriv[0][0][0]), walk_num*elec_num*4)
                   );
  assert(rc == QMCKL_SUCCESS);

  double total_j[walk_num];
  rc = qmckl_check(context,
                   qmckl_get_jastrow_champ_value(context, &(total_j[0]), walk_num)
                   );
  assert(rc == QMCKL_SUCCESS);


  for (int64_t k=0 ; k< walk_num ; ++k) {
    for (int64_t m=0 ; m<4; ++m) {
      for (int64_t e=0 ; e<elec_num; ++e) {
        if (m < 3) { /* test only gradients */
          if (fabs(total_j_deriv[k][m][e]/total_j[k] - (factor_ee_gl[k][m][e] + factor_en_gl[k][m][e] + factor_een_gl[k][m][e])) > 1e-12) {
            printf("k = %ld\n", k);
            printf("m = %ld\n", m);
            printf("e = %ld\n", e);
            printf("total_j_deriv/total_j = %20.15e\n", total_j_deriv[k][m][e]/total_j[k]);
            printf("factor_ee_gl + factor_en_gl + factor_een_gl = %20.15e\n", factor_ee_gl[k][m][e] + factor_en_gl[k][m][e] + factor_een_gl[k][m][e]);
            fflush(stdout);
          }
          assert (fabs(total_j_deriv[k][m][e]/total_j[k] - (factor_ee_gl[k][m][e] + factor_en_gl[k][m][e] + factor_een_gl[k][m][e])) < 1.e-8);
        }
      }
    }
  }

printf("Total Jastrow gradient only\n");
/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

{
  double total_j_grad[walk_num][3][elec_num];
  rc = qmckl_check(context,
                  qmckl_get_jastrow_champ_grad(context, &(total_j_grad[0][0][0]), walk_num*elec_num*3)
                  );
  assert(rc == QMCKL_SUCCESS);


  double total_j_deriv[walk_num][4][elec_num];
  rc = qmckl_check(context,
                  qmckl_get_jastrow_champ_gl(context, &(total_j_deriv[0][0][0]), walk_num*elec_num*4)
                  );
  assert(rc == QMCKL_SUCCESS);



  for (int64_t k=0 ; k< walk_num ; ++k) {
    for (int64_t m=0 ; m<3; ++m) {
      for (int64_t e=0 ; e<elec_num; ++e) {
        if (fabs(total_j_grad[k][m][e] - total_j_deriv[k][m][e]) > 1e-12) {
          printf("%ld %ld %ld\n", k, m, e);
          printf("total_j_deriv = %20.15e\n", total_j_deriv[k][m][e]);
          printf("total_j_grad  = %20.15e\n", total_j_grad[k][m][e]);
          fflush(stdout);
        }
        assert (fabs(total_j_deriv[k][m][e] - total_j_grad[k][m][e]) < 1.e-8);
      }
    }
  }
}

}

rc = qmckl_context_destroy(context);
    assert (rc == QMCKL_SUCCESS);

    return 0;
}
