#include "qmckl.h"
#include <assert.h>
#include <math.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdbool.h>
#include <stdio.h>
#include "n2.h"
#include "qmckl_jastrow_champ_private_func.h"
#include "qmckl_jastrow_champ_single_private_func.h"

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

double new_coords[6] = {1.0,2.0,3.0,4.0,5.0,6.0};

double coords[walk_num][elec_num][3];

printf("Single e-e distance\n");

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

double ee_distance[walk_num][elec_num][elec_num];
double single_ee_distance[walk_num][elec_num];

for (int elec = 0; elec < elec_num; elec++){

  rc = qmckl_set_electron_coord(context, 'N', walk_num, elec_coord, walk_num*elec_num*3);
  assert(rc == QMCKL_SUCCESS);

  double ee_distance[walk_num][elec_num][elec_num];
  rc = qmckl_get_electron_ee_distance(context, &ee_distance[0][0][0], walk_num*elec_num*elec_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_electron_coord(context, 'N', &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_set_single_point(context, 'N', elec, new_coords, 3*walk_num);
  assert (rc == QMCKL_SUCCESS);

  double single_ee_distance[walk_num][elec_num];
  rc = qmckl_get_single_electron_ee_distance(context,&single_ee_distance[0][0],walk_num*elec_num);
  assert (rc == QMCKL_SUCCESS);

  coords[0][elec][0] = new_coords[0];
  coords[0][elec][1] = new_coords[1];
  coords[0][elec][2] = new_coords[2];
  coords[1][elec][0] = new_coords[3];
  coords[1][elec][1] = new_coords[4];
  coords[1][elec][2] = new_coords[5];

  rc = qmckl_set_electron_coord(context, 'N', walk_num, &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_electron_ee_distance(context, &ee_distance[0][0][0], walk_num*elec_num*elec_num);
  assert (rc == QMCKL_SUCCESS);

  for (int nw = 0; nw < walk_num; nw++){
    for (int i = 0; i < elec_num; i++) {
      if (i == elec) continue;
      assert(fabs((ee_distance[nw][elec][i]-single_ee_distance[nw][i])) < 1.e-12);
      assert(fabs((ee_distance[nw][i][elec]-single_ee_distance[nw][i])) < 1.e-12);
    }
  }
}

printf("OK\n");

printf("Single e-n distance\n");

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

double en_distance[walk_num][elec_num][nucl_num];
double single_en_distance[walk_num][nucl_num];

for (int elec = 0; elec < elec_num; elec++){

  rc = qmckl_set_electron_coord(context, 'N', walk_num, elec_coord, walk_num*elec_num*3);
  assert(rc == QMCKL_SUCCESS);

  rc = qmckl_get_electron_en_distance(context, &en_distance[0][0][0],walk_num*elec_num*nucl_num);

  rc = qmckl_get_electron_coord(context, 'N', &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_set_single_point(context, 'N', elec, new_coords, 3*walk_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_single_electron_en_distance(context, &single_en_distance[0][0],nucl_num*walk_num);
  assert (rc == QMCKL_SUCCESS);

  coords[0][elec][0] = new_coords[0];
  coords[0][elec][1] = new_coords[1];
  coords[0][elec][2] = new_coords[2];
  coords[1][elec][0] = new_coords[3];
  coords[1][elec][1] = new_coords[4];
  coords[1][elec][2] = new_coords[5];

	rc = qmckl_set_electron_coord(context, 'N', walk_num, &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_electron_en_distance(context, &en_distance[0][0][0], walk_num*elec_num*nucl_num);
  assert (rc == QMCKL_SUCCESS);

  for (int nw = 0 ; nw < walk_num ; nw++) {
    for (int a = 0; a < nucl_num; a++){
      assert(fabs((en_distance[nw][elec][a]-single_en_distance[nw][a])) < 1.e-12);
    }
  }
}

printf("OK\n");

printf("Single een rescaled e-e distance\n");

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

double rescaled_een_ee_distance[walk_num][cord_num+1][elec_num][elec_num];
double single_rescaled_een_ee_distance[walk_num][cord_num+1][elec_num];

for (int elec = 0; elec < elec_num; elec++){

  rc = qmckl_set_electron_coord(context, 'N', walk_num, elec_coord, walk_num*elec_num*3);
  assert(rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_een_rescaled_e(context,  &rescaled_een_ee_distance[0][0][0][0], walk_num*(cord_num+1)*elec_num*elec_num);

  rc = qmckl_get_electron_coord(context, 'N', &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_set_single_point(context, 'N', elec, new_coords, 3*walk_num);
  assert (rc == QMCKL_SUCCESS);


  rc = qmckl_get_een_rescaled_single_e(context, &single_rescaled_een_ee_distance[0][0][0], walk_num*(cord_num+1)*elec_num);
  assert (rc == QMCKL_SUCCESS);

  coords[0][elec][0] = new_coords[0];
  coords[0][elec][1] = new_coords[1];
  coords[0][elec][2] = new_coords[2];
  coords[1][elec][0] = new_coords[3];
  coords[1][elec][1] = new_coords[4];
  coords[1][elec][2] = new_coords[5];

  rc = qmckl_set_electron_coord(context, 'N', walk_num, &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_een_rescaled_e(context,  &rescaled_een_ee_distance[0][0][0][0], walk_num*(cord_num+1)*elec_num*elec_num);
  assert (rc == QMCKL_SUCCESS);

  for (int nw = 0; nw < walk_num; nw++){
    for (int l = 0; l <= cord_num; l++){
      for (int i = 0; i < elec_num; i++) {
        if (i == elec) continue;
        assert(fabs((rescaled_een_ee_distance[nw][l][elec][i]-single_rescaled_een_ee_distance[nw][l][i])) < 1.e-12);
        assert(fabs((rescaled_een_ee_distance[nw][l][i][elec]-single_rescaled_een_ee_distance[nw][l][i])) < 1.e-12);
      }
    }
  }

}

printf("OK\n");

printf("Single een rescaled e-n distance\n");

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

double rescaled_een_en_distance[walk_num][cord_num+1][nucl_num][elec_num];
double single_rescaled_een_en_distance[walk_num][cord_num+1][nucl_num];

for (int elec = 0; elec < elec_num; elec++){

	rc = qmckl_set_electron_coord(context, 'N', walk_num, elec_coord, walk_num*elec_num*3);
  assert(rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_een_rescaled_n(context,  &rescaled_een_en_distance[0][0][0][0], walk_num*(cord_num+1)*nucl_num*elec_num);

  rc = qmckl_get_electron_coord(context, 'N', &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_set_single_point(context, 'N', elec, new_coords, 3*walk_num);
  assert (rc == QMCKL_SUCCESS);


  rc = qmckl_get_een_rescaled_single_n(context, &single_rescaled_een_en_distance[0][0][0], walk_num*(cord_num+1)*nucl_num);
  assert (rc == QMCKL_SUCCESS);

  coords[0][elec][0] = new_coords[0];
  coords[0][elec][1] = new_coords[1];
  coords[0][elec][2] = new_coords[2];
  coords[1][elec][0] = new_coords[3];
  coords[1][elec][1] = new_coords[4];
  coords[1][elec][2] = new_coords[5];

	rc = qmckl_set_electron_coord(context, 'N', walk_num, &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_een_rescaled_n(context,  &rescaled_een_en_distance[0][0][0][0], walk_num*(cord_num+1)*nucl_num*elec_num);
  assert (rc == QMCKL_SUCCESS);

  for (int nw = 0; nw < walk_num; nw++){
    for (int l = 0; l <= cord_num; l++){
      for (int a = 0; a < nucl_num; a++) {
        assert(fabs((rescaled_een_en_distance[nw][l][a][elec]-single_rescaled_een_en_distance[nw][l][a])) < 1.e-12);
      }
    }
  }
}

printf("OK\n");

printf("Delta p\n");

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

double p_old[walk_num][cord_num][cord_num+1][nucl_num][elec_num];
double delta_p[walk_num][cord_num][cord_num+1][nucl_num][elec_num];
double p_new[walk_num][cord_num][cord_num+1][nucl_num][elec_num];

for (int elec = 0; elec < elec_num; elec++){

	rc = qmckl_set_electron_coord(context, 'N', walk_num, elec_coord, walk_num*elec_num*3);
  assert(rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_tmp_c(context,  &p_old[0][0][0][0][0], walk_num*cord_num*(cord_num+1)*nucl_num*elec_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_electron_coord(context, 'N', &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_set_single_point(context, 'N', elec, new_coords, 3*walk_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_delta_p(context, &delta_p[0][0][0][0][0], walk_num*cord_num*(cord_num+1)*nucl_num*elec_num);
  assert (rc == QMCKL_SUCCESS);

  coords[0][elec][0] = new_coords[0];
  coords[0][elec][1] = new_coords[1];
  coords[0][elec][2] = new_coords[2];
  coords[1][elec][0] = new_coords[3];
  coords[1][elec][1] = new_coords[4];
  coords[1][elec][2] = new_coords[5];

	rc = qmckl_set_electron_coord(context, 'N', walk_num, &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_tmp_c(context,  &p_new[0][0][0][0][0], walk_num*cord_num*(cord_num+1)*nucl_num*elec_num);
  assert (rc == QMCKL_SUCCESS);

  for (int nw = 0; nw < walk_num; nw++){
    for (int l = 0; l < cord_num; l++){
      for (int m = 0; m <= cord_num; m++){
        for (int a = 0; a < nucl_num; a++) {
          for (int i = 0; i < elec_num; i++){
            assert(fabs(((p_new[nw][l][m][a][i]-p_old[nw][l][m][a][i])-delta_p[nw][l][m][a][i])) < 1.e-12);
          }
        }
      }
    }
  }
}

printf("OK\n");

printf("Delta een\n");

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

double jastrow_een_old[walk_num];
double delta_een[walk_num];
double jastrow_een_new[walk_num];

for (int elec = 0; elec < elec_num; elec++){

	rc = qmckl_set_electron_coord(context, 'N', walk_num, elec_coord, walk_num*elec_num*3);
  assert(rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_factor_een(context, &jastrow_een_old[0], walk_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_electron_coord(context, 'N', &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_set_single_point(context, 'N', elec, new_coords, 3*walk_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_single_een(context, &delta_een[0], walk_num);
  assert (rc == QMCKL_SUCCESS);

  coords[0][elec][0] = new_coords[0];
  coords[0][elec][1] = new_coords[1];
  coords[0][elec][2] = new_coords[2];
  coords[1][elec][0] = new_coords[3];
  coords[1][elec][1] = new_coords[4];
  coords[1][elec][2] = new_coords[5];

	rc = qmckl_set_electron_coord(context, 'N', walk_num, &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_factor_een(context, &jastrow_een_new[0], walk_num);
  assert (rc == QMCKL_SUCCESS);

  for (int nw = 0; nw < walk_num; nw++) {
    //printf("jastrow_een_old %f\n", jastrow_een_old[nw]);
    //printf("jastrow_een_new %f\n", jastrow_een_new[nw]);
    //printf("delta_een %f\n", delta_een[nw]);
    assert(fabs((jastrow_een_new[nw]-jastrow_een_old[nw])-delta_een[nw]) < 1.e-12);
  }

}

printf("OK\n");

printf("Een rescaled single n gl\n");

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

double een_rescaled_en_gl[walk_num][cord_num+1][nucl_num][4][elec_num];
double een_rescaled_single_n_gl[walk_num][cord_num+1][nucl_num][4];

for (int elec = 0; elec < elec_num; elec++){

  rc = qmckl_set_electron_coord(context, 'N', walk_num, elec_coord, walk_num*elec_num*3);
  assert(rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_een_rescaled_n_gl(context, &een_rescaled_en_gl[0][0][0][0][0], walk_num*(cord_num+1)*nucl_num*elec_num*4);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_electron_coord(context, 'N', &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_set_single_point(context, 'N', elec, new_coords, 3*walk_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_een_rescaled_single_n_gl(context, &een_rescaled_single_n_gl[0][0][0][0], walk_num*(cord_num+1)*nucl_num*4);
  assert (rc == QMCKL_SUCCESS);

  coords[0][elec][0] = new_coords[0];
  coords[0][elec][1] = new_coords[1];
  coords[0][elec][2] = new_coords[2];
  coords[1][elec][0] = new_coords[3];
  coords[1][elec][1] = new_coords[4];
  coords[1][elec][2] = new_coords[5];

	rc = qmckl_set_electron_coord(context, 'N', walk_num, &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_een_rescaled_n_gl(context, &een_rescaled_en_gl[0][0][0][0][0], walk_num*(cord_num+1)*nucl_num*elec_num*4);
  assert (rc == QMCKL_SUCCESS);

  for (int nw = 0; nw < walk_num; nw++) {
    for (int l = 0; l < cord_num+1; l++) {
      for (int a = 0; a < nucl_num; a++) {
        for (int m = 0; m < 4; m++) {
          //printf("nw %d l %d a %d m %d\n", nw, l, a, m);
          //printf(" %f %f\n", een_rescaled_en_gl[nw][l][a][m][elec], een_rescaled_single_n_gl[nw][l][a][m]);
          assert(fabs(een_rescaled_en_gl[nw][l][a][m][elec] - een_rescaled_single_n_gl[nw][l][a][m]) < 1.e-12);
        }
      }
    }
  }
}

printf("OK\n");

printf("Een rescaled single e gl\n");

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

double een_rescaled_ee_gl[walk_num][cord_num+1][elec_num][4][elec_num];
double een_rescaled_single_e_gl[walk_num][cord_num+1][elec_num][4];

for (int elec = 0; elec < elec_num; elec++){

	rc = qmckl_set_electron_coord(context, 'N', walk_num, elec_coord, walk_num*elec_num*3);
  assert(rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_een_rescaled_e_gl(context, &een_rescaled_ee_gl[0][0][0][0][0], walk_num*(cord_num+1)*elec_num*elec_num*4);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_electron_coord(context, 'N', &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_set_single_point(context, 'N', elec, new_coords, 3*walk_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_een_rescaled_single_e_gl(context, &een_rescaled_single_e_gl[0][0][0][0], walk_num*(cord_num+1)*elec_num*4);
  assert (rc == QMCKL_SUCCESS);

  coords[0][elec][0] = new_coords[0];
  coords[0][elec][1] = new_coords[1];
  coords[0][elec][2] = new_coords[2];
  coords[1][elec][0] = new_coords[3];
  coords[1][elec][1] = new_coords[4];
  coords[1][elec][2] = new_coords[5];

	rc = qmckl_set_electron_coord(context, 'N', walk_num, &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_een_rescaled_e_gl(context, &een_rescaled_ee_gl[0][0][0][0][0], walk_num*(cord_num+1)*elec_num*elec_num*4);
  assert (rc == QMCKL_SUCCESS);

  double metric[4] = {-1.0, -1.0, -1.0, 1.0};

  for (int l = 0; l < cord_num+1; l++) {
    for (int nw = 0; nw < walk_num; nw++) {
      for (int i = 0; i < elec_num; i++) {
        for (int m = 0; m < 4; m++) {
          //printf("een_rescaled_ee_gl[nw][l][i][m][elec] %i %i %i %f \n", l, m ,i, een_rescaled_ee_gl[nw][l][i][m][elec]);
          //printf("een_rescaled_ee_gl[nw][l][elec][m][i] %i %i %i %f \n", l, m ,i, een_rescaled_ee_gl[nw][l][elec][m][i]);
          //printf("een_rescaled_single_e_gl[nw][l][i][m] %i %i %i %f\n", l, m, i,een_rescaled_single_e_gl[nw][l][i][m]);
          assert(fabs(een_rescaled_ee_gl[nw][l][i][m][elec] - een_rescaled_single_e_gl[nw][l][i][m]) < 1.e-12);
          assert(fabs(een_rescaled_ee_gl[nw][l][elec][m][i] - metric[m] * een_rescaled_single_e_gl[nw][l][i][m]) < 1.e-12);
        }
      }
    }
  }
}

printf("OK\n");

printf("Delta P gl\n");

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

double p_gl_old[walk_num][cord_num][cord_num+1][nucl_num][4][elec_num];
double delta_p_gl[walk_num][cord_num][cord_num+1][4][nucl_num][elec_num];
double p_gl_new[walk_num][cord_num][cord_num+1][nucl_num][4][elec_num];


for (int elec = 0; elec < elec_num; elec++){

	rc = qmckl_set_electron_coord(context, 'N', walk_num, elec_coord, walk_num*elec_num*3);
  assert(rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_dtmp_c(context,  &p_gl_old[0][0][0][0][0][0], walk_num*cord_num*(cord_num+1)*nucl_num*4*elec_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_electron_coord(context, 'N', &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_set_single_point(context, 'N', elec, new_coords, 3*walk_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_delta_p_gl(context, &delta_p_gl[0][0][0][0][0][0], 4*walk_num*cord_num*(cord_num+1)*nucl_num*elec_num);

  assert (rc == QMCKL_SUCCESS);

  coords[0][elec][0] = new_coords[0];
  coords[0][elec][1] = new_coords[1];
  coords[0][elec][2] = new_coords[2];
  coords[1][elec][0] = new_coords[3];
  coords[1][elec][1] = new_coords[4];
  coords[1][elec][2] = new_coords[5];

	rc = qmckl_set_electron_coord(context, 'N', walk_num, &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_dtmp_c(context, &p_gl_new[0][0][0][0][0][0], walk_num*cord_num*(cord_num+1)*nucl_num*4*elec_num);
  assert (rc == QMCKL_SUCCESS);

  for (int nw = 0; nw < walk_num; nw++){
    for (int l = 0; l < cord_num; l++){
      for (int m = 0; m <= cord_num; m++){
        for (int a = 0; a < nucl_num; a++) {
          for (int i = 0; i < elec_num; i++){
            for (int k = 0; k < 4; k++){
              if (fabs(((p_gl_new[nw][l][m][a][k][i]-p_gl_old[nw][l][m][a][k][i])-delta_p_gl[nw][l][m][k][a][i])) > 1.e-12) {
                printf("p_gl[%d][%d][%d][%d][%d][%d] = %f\n", nw, l, m, a, k, i, p_gl_new[nw][l][m][a][k][i] - p_gl_old[nw][l][m][a][k][i]);
                printf("delta_p_gl[%d][%d][%d][%d][%d][%d] = %f\n", nw, l, m, a, k, i, delta_p_gl[nw][l][m][k][a][i]);
              }
              assert(fabs(((p_gl_new[nw][l][m][a][k][i]-p_gl_old[nw][l][m][a][k][i])-delta_p_gl[nw][l][m][k][a][i])) < 1.e-12);
            }
          }
        }
      }
    }
  }
}

printf("OK\n");

printf("Delta een gl\n");

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

double een_gl_old[walk_num][4][elec_num];
double delta_een_gl[walk_num][4][elec_num];
double een_gl_new[walk_num][4][elec_num];


for (int elec = 0; elec < elec_num; elec++){

	rc = qmckl_set_electron_coord(context, 'N', walk_num, elec_coord, walk_num*elec_num*3);
  assert(rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_factor_een_gl(context, &een_gl_old[0][0][0], walk_num*4*elec_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_electron_coord(context, 'N', &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_set_single_point(context, 'N', elec, new_coords, 3*walk_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_single_een_gl(context, &delta_een_gl[0][0][0], walk_num*4*elec_num);
  assert (rc == QMCKL_SUCCESS);

  coords[0][elec][0] = new_coords[0];
  coords[0][elec][1] = new_coords[1];
  coords[0][elec][2] = new_coords[2];
  coords[1][elec][0] = new_coords[3];
  coords[1][elec][1] = new_coords[4];
  coords[1][elec][2] = new_coords[5];

	rc = qmckl_set_electron_coord(context, 'N', walk_num, &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_factor_een_gl(context, &een_gl_new[0][0][0], walk_num*4*elec_num);
  assert (rc == QMCKL_SUCCESS);

  for (int nw = 0; nw < walk_num; nw++) {
    for (int m = 0; m < 4; m++) {
      for (int i = 0; i < elec_num; i++) {
        //printf("delta_een_gl[%d][%d][%d] = %f\n", nw, i, m, delta_een_gl[nw][i][m]);
        //printf("een_gl_[%d][%d][%d] = %f\n", nw, m,i, een_gl_new[nw][m][i]-een_gl_old[nw][m][i]);

        assert(fabs((een_gl_new[nw][m][i]- een_gl_old[nw][m][i]) - delta_een_gl[nw][m][i]) < 1.e-12);

      }
    }
  }
}

printf("OK\n");

printf("Delta P g\n");

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

double delta_p_g[walk_num][cord_num][cord_num+1][4][nucl_num][elec_num];


for (int elec = 0; elec < elec_num; elec++){

	rc = qmckl_set_electron_coord(context, 'N', walk_num, elec_coord, walk_num*elec_num*3);
  assert(rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_dtmp_c(context,  &p_gl_old[0][0][0][0][0][0], walk_num*cord_num*(cord_num+1)*nucl_num*4*elec_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_electron_coord(context, 'N', &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_set_single_point(context, 'N', elec, new_coords, 3*walk_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_delta_p_g(context, &delta_p_g[0][0][0][0][0][0], 4*walk_num*cord_num*(cord_num+1)*nucl_num*elec_num);

  assert (rc == QMCKL_SUCCESS);

  coords[0][elec][0] = new_coords[0];
  coords[0][elec][1] = new_coords[1];
  coords[0][elec][2] = new_coords[2];
  coords[1][elec][0] = new_coords[3];
  coords[1][elec][1] = new_coords[4];
  coords[1][elec][2] = new_coords[5];

	rc = qmckl_set_electron_coord(context, 'N', walk_num, &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_dtmp_c(context, &p_gl_new[0][0][0][0][0][0], walk_num*cord_num*(cord_num+1)*nucl_num*4*elec_num);
  assert (rc == QMCKL_SUCCESS);

  for (int nw = 0; nw < walk_num; nw++){
    for (int l = 0; l < cord_num; l++){
      for (int m = 0; m <= cord_num; m++){
        for (int a = 0; a < nucl_num; a++) {
            for (int k = 0; k < 3; k++){
              if (fabs(((p_gl_new[nw][l][m][a][k][elec]-p_gl_old[nw][l][m][a][k][elec])-delta_p_g[nw][l][m][k][a][elec])) > 1.e-12) {
                printf("p_gl[%d][%d][%d][%d][%d][%d] = %f\n", nw, l, m, a, k, elec, p_gl_new[nw][l][m][a][k][elec] - p_gl_old[nw][l][m][a][k][elec]);
                printf("delta_p_g[%d][%d][%d][%d][%d][%d] = %f\n", nw, l, m, a, k, elec, delta_p_g[nw][l][m][k][a][elec]);
              }
              assert(fabs(((p_gl_new[nw][l][m][a][k][elec]-p_gl_old[nw][l][m][a][k][elec])-delta_p_g[nw][l][m][k][a][elec])) < 1.e-12);
            }
        }
      }
    }
  }
}

printf("OK\n");

printf("Delta een g\n");

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

double een_g_old[walk_num][4][elec_num];
double delta_een_g[walk_num][4][elec_num];
double een_g_new[walk_num][4][elec_num];


for (int elec = 0; elec < elec_num; elec++){

	rc = qmckl_set_electron_coord(context, 'N', walk_num, elec_coord, walk_num*elec_num*3);
  assert(rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_factor_een_gl(context, &een_g_old[0][0][0], walk_num*4*elec_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_electron_coord(context, 'N', &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_set_single_point(context, 'N', elec, new_coords, 3*walk_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_single_een_g(context, &delta_een_g[0][0][0], walk_num*4*elec_num);
  assert (rc == QMCKL_SUCCESS);

  coords[0][elec][0] = new_coords[0];
  coords[0][elec][1] = new_coords[1];
  coords[0][elec][2] = new_coords[2];
  coords[1][elec][0] = new_coords[3];
  coords[1][elec][1] = new_coords[4];
  coords[1][elec][2] = new_coords[5];

	rc = qmckl_set_electron_coord(context, 'N', walk_num, &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_factor_een_gl(context, &een_g_new[0][0][0], walk_num*4*elec_num);
  assert (rc == QMCKL_SUCCESS);

  for (int nw = 0; nw < walk_num; nw++) {
    for (int m = 0; m < 3; m++) {
      //for (int i = 0; i < elec_num; i++) {
        //printf("delta_een_g[%d][%d][%d] = %f\n", nw, i, m, delta_een_g[nw][i][m]);
        //printf("een_g_[%d][%d][%d] = %f\n", nw, m,i, een_g_new[nw][m][i]-een_g_old[nw][m][i]);

        assert(fabs((een_g_new[nw][m][elec]- een_g_old[nw][m][elec]) - delta_een_g[nw][m][elec]) < 1.e-12);

      //}
    }
  }
}

printf("OK\n");

printf("ee rescaled single\n");

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

double ee_rescaled[walk_num][elec_num][elec_num];
double single_ee_rescaled[walk_num][elec_num];

for (int elec = 0; elec < elec_num; elec++){

	rc = qmckl_set_electron_coord(context, 'N', walk_num, elec_coord, walk_num*elec_num*3);
  assert(rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_ee_distance_rescaled(context, &ee_rescaled[0][0][0], walk_num*elec_num*elec_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_electron_coord(context, 'N', &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_set_single_point(context, 'N', elec, new_coords, 3*walk_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_ee_rescaled_single(context, &single_ee_rescaled[0][0], walk_num*elec_num);
  assert (rc == QMCKL_SUCCESS);

  coords[0][elec][0] = new_coords[0];
  coords[0][elec][1] = new_coords[1];
  coords[0][elec][2] = new_coords[2];
  coords[1][elec][0] = new_coords[3];
  coords[1][elec][1] = new_coords[4];
  coords[1][elec][2] = new_coords[5];

	rc = qmckl_set_electron_coord(context, 'N', walk_num, &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_ee_distance_rescaled(context, &ee_rescaled[0][0][0], walk_num*elec_num*elec_num);
  assert (rc == QMCKL_SUCCESS);


  for (int nw = 0; nw < walk_num; nw++) {
    for (int i = 0; i < elec_num; i++){
      //printf("nw %d i %d %f %f\n", nw, i, ee_rescaled[nw][2][i], single_ee_rescaled[nw][i]);
      assert(fabs(ee_rescaled[nw][elec][i]-single_ee_rescaled[nw][i]) < 1.e-12);
      assert(fabs(ee_rescaled[nw][i][elec]-single_ee_rescaled[nw][i]) < 1.e-12);
    }
  }
}

printf("OK\n");

printf("Delta ee\n");

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

double jastrow_ee_old[walk_num];
double delta_ee[walk_num];
double jastrow_ee_new[walk_num];

for (int elec = 0; elec < elec_num; elec++){

	rc = qmckl_set_electron_coord(context, 'N', walk_num, elec_coord, walk_num*elec_num*3);
  assert(rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_factor_ee(context, &jastrow_ee_old[0], walk_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_electron_coord(context, 'N', &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_set_single_point(context, 'N', elec, new_coords, 3*walk_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_single_ee(context, &delta_ee[0], walk_num);
  assert (rc == QMCKL_SUCCESS);

  coords[0][elec][0] = new_coords[0];
  coords[0][elec][1] = new_coords[1];
  coords[0][elec][2] = new_coords[2];
  coords[1][elec][0] = new_coords[3];
  coords[1][elec][1] = new_coords[4];
  coords[1][elec][2] = new_coords[5];

	rc = qmckl_set_electron_coord(context, 'N', walk_num, &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_factor_ee(context, &jastrow_ee_new[0], walk_num);
  assert (rc == QMCKL_SUCCESS);

  for (int nw = 0; nw < walk_num; nw++) {
    //printf("%f %f %f %3.14f\n", jastrow_ee_new[nw], jastrow_ee_old[nw], delta_ee[nw], fabs((jastrow_ee_new[nw] - jastrow_ee_old[nw]) - delta_ee[nw]));
    assert(fabs((jastrow_ee_new[nw] - jastrow_ee_old[nw]) - delta_ee[nw]) < 1.e-12);
  }
}

printf("OK\n");

printf("ee rescaled single gl\n");     


/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

double ee_rescaled_gl[walk_num][elec_num][elec_num][4];
double single_ee_rescaled_gl[walk_num][elec_num][4];


for (int elec = 0; elec < elec_num; elec++){

	rc = qmckl_set_electron_coord(context, 'N', walk_num, elec_coord, walk_num*elec_num*3);
  assert(rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_ee_distance_rescaled_gl(context, &ee_rescaled_gl[0][0][0][0], walk_num*elec_num*elec_num*4);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_electron_coord(context, 'N', &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_set_single_point(context, 'N', elec, new_coords, 3*walk_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_ee_rescaled_single_gl(context, &single_ee_rescaled_gl[0][0][0], walk_num*elec_num*4);
  assert (rc == QMCKL_SUCCESS);

  coords[0][elec][0] = new_coords[0];
  coords[0][elec][1] = new_coords[1];
  coords[0][elec][2] = new_coords[2];
  coords[1][elec][0] = new_coords[3];
  coords[1][elec][1] = new_coords[4];
  coords[1][elec][2] = new_coords[5];

	rc = qmckl_set_electron_coord(context, 'N', walk_num, &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_ee_distance_rescaled_gl(context, &ee_rescaled_gl[0][0][0][0], walk_num*elec_num*elec_num*4);
  assert (rc == QMCKL_SUCCESS);

  double metric[4] = {-1.0, -1.0, -1.0, 1.0};

  for (int nw = 0; nw < walk_num; nw++) {
    for (int i = 0; i < elec_num; i++) {
      for (int m = 0; m < 4; m++) {
        if (i == elec) continue;
        //printf("%f\n", ee_rescaled_gl[nw][elec][i][m]);
        //printf("%f\n", single_ee_rescaled_gl[nw][i][m]);
        assert(fabs(ee_rescaled_gl[nw][elec][i][m] - single_ee_rescaled_gl[nw][i][m]) < 1.e-12);
        assert(fabs(ee_rescaled_gl[nw][i][elec][m] - metric[m] * single_ee_rescaled_gl[nw][i][m]) < 1.e-12);
      }
    }
  }
}

printf("OK\n");

printf("Delta ee gl\n");

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

double ee_gl_old[walk_num][4][elec_num];
double delta_ee_gl[walk_num][elec_num][4];
double ee_gl_new[walk_num][4][elec_num];

for (int elec = 0; elec < elec_num; elec++){

	rc = qmckl_set_electron_coord(context, 'N', walk_num, elec_coord, walk_num*elec_num*3);
  assert(rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_factor_ee_gl(context, &ee_gl_old[0][0][0], walk_num*elec_num*4);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_electron_coord(context, 'N', &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_set_single_point(context, 'N', elec, new_coords, 3*walk_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_single_ee_gl(context, &delta_ee_gl[0][0][0], walk_num*elec_num*4);
  assert (rc == QMCKL_SUCCESS);

  coords[0][elec][0] = new_coords[0];
  coords[0][elec][1] = new_coords[1];
  coords[0][elec][2] = new_coords[2];
  coords[1][elec][0] = new_coords[3];
  coords[1][elec][1] = new_coords[4];
  coords[1][elec][2] = new_coords[5];

	rc = qmckl_set_electron_coord(context, 'N', walk_num, &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_factor_ee_gl(context, &ee_gl_new[0][0][0], walk_num*elec_num*4);
  assert (rc == QMCKL_SUCCESS);

  for (int nw = 0; nw < walk_num; nw++) {
    for (int i = 0; i < elec_num; i++) {
      for (int m = 0; m < 4; m++) {
        //printf("%f\n",(ee_gl_new[nw][m][i] - ee_gl_old[nw][m][i]));
        //printf("%f\n",delta_ee_gl[nw][i][m]);
        assert(fabs((ee_gl_new[nw][m][i] - ee_gl_old[nw][m][i]) - delta_ee_gl[nw][i][m]) < 1.e-12);
      }
    }
  }
}

printf("OK\n");

printf("En rescaled single\n");


/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

double en_rescaled[walk_num][nucl_num][elec_num];
double single_en_rescaled[walk_num][nucl_num];

for (int elec = 0; elec < elec_num; elec++){

	rc = qmckl_set_electron_coord(context, 'N', walk_num, elec_coord, walk_num*elec_num*3);
  assert(rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_en_distance_rescaled(context, &en_rescaled[0][0][0], walk_num*nucl_num*elec_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_electron_coord(context, 'N', &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_set_single_point(context, 'N', elec, new_coords, 3*walk_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_en_rescaled_single(context, &single_en_rescaled[0][0], walk_num*nucl_num);
  assert (rc == QMCKL_SUCCESS);

  coords[0][elec][0] = new_coords[0];
  coords[0][elec][1] = new_coords[1];
  coords[0][elec][2] = new_coords[2];
  coords[1][elec][0] = new_coords[3];
  coords[1][elec][1] = new_coords[4];
  coords[1][elec][2] = new_coords[5];

	rc = qmckl_set_electron_coord(context, 'N', walk_num, &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_en_distance_rescaled(context, &en_rescaled[0][0][0], walk_num*nucl_num*elec_num);
  assert (rc == QMCKL_SUCCESS);

  for (int nw = 0; nw < walk_num; nw++) {
    for (int a = 0; a < nucl_num; a++){
      assert(fabs(en_rescaled[nw][a][elec]-single_en_rescaled[nw][a]) < 1.e-12);
    }
  }
}

printf("OK\n");

printf("Delta en\n");

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

double jastrow_en_old[walk_num];
double delta_en[walk_num];
double jastrow_en_new[walk_num];

for (int elec = 0; elec < elec_num; elec++){

	rc = qmckl_set_electron_coord(context, 'N', walk_num, elec_coord, walk_num*elec_num*3);
  assert(rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_factor_en(context, &jastrow_en_old[0], walk_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_electron_coord(context, 'N', &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_set_single_point(context, 'N', elec, new_coords, 3*walk_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_single_en(context, &delta_en[0], walk_num);
  assert (rc == QMCKL_SUCCESS);

  coords[0][elec][0] = new_coords[0];
  coords[0][elec][1] = new_coords[1];
  coords[0][elec][2] = new_coords[2];
  coords[1][elec][0] = new_coords[3];
  coords[1][elec][1] = new_coords[4];
  coords[1][elec][2] = new_coords[5];

  rc = qmckl_set_electron_coord(context, 'N', walk_num, &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_factor_en(context, &jastrow_en_new[0], walk_num);
  assert (rc == QMCKL_SUCCESS);

  for (int nw = 0; nw < walk_num; nw++) {
    //printf("electron %d walk %d \n", elec, nw);
    //printf("jastrow_en_new %f\n", jastrow_en_new[nw]);
    //printf("jastrow_en_old %f\n", jastrow_en_old[nw]);
    //printf("delta_en %f\n", delta_en[nw]);
    //printf("diff %f\n", jastrow_en_new[nw] - jastrow_en_old[nw] - delta_en[nw]);
    assert(fabs((jastrow_en_new[nw] - jastrow_en_old[nw]) - delta_en[nw]) < 1.e-12);
  }
}

printf("OK\n");

printf("En rescaled single gl\n");

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

double en_rescaled_gl[walk_num][nucl_num][elec_num][4];
double single_en_rescaled_gl[walk_num][nucl_num][4];

for (int elec = 0; elec < elec_num; elec++){

	rc = qmckl_set_electron_coord(context, 'N', walk_num, elec_coord, walk_num*elec_num*3);
  assert(rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_en_distance_rescaled_gl(context, &en_rescaled_gl[0][0][0][0], walk_num*nucl_num*elec_num*4);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_electron_coord(context, 'N', &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_set_single_point(context, 'N', elec, new_coords, 3*walk_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_en_rescaled_single_gl(context, &single_en_rescaled_gl[0][0][0], walk_num*nucl_num*4);
  assert (rc == QMCKL_SUCCESS);

  coords[0][elec][0] = new_coords[0];
  coords[0][elec][1] = new_coords[1];
  coords[0][elec][2] = new_coords[2];
  coords[1][elec][0] = new_coords[3];
  coords[1][elec][1] = new_coords[4];
  coords[1][elec][2] = new_coords[5];

	rc = qmckl_set_electron_coord(context, 'N', walk_num, &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_en_distance_rescaled_gl(context, &en_rescaled_gl[0][0][0][0], walk_num*nucl_num*elec_num*4);
  assert (rc == QMCKL_SUCCESS);

  for (int nw = 0; nw < walk_num; nw++) {
    for (int a = 0; a < nucl_num; a++) {
      for (int m = 0; m < 4; m++) {
        assert(fabs(en_rescaled_gl[nw][a][elec][m] - single_en_rescaled_gl[nw][a][m]) < 1.e-12);
      }
    }
  }
}

printf("OK\n");

printf("Delta en gl\n");

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

double en_gl_old[walk_num][4][elec_num];
double delta_en_gl[walk_num][elec_num][4];
double en_gl_new[walk_num][4][elec_num];

for (int elec = 0; elec < elec_num; elec++){

	rc = qmckl_set_electron_coord(context, 'N', walk_num, elec_coord, walk_num*elec_num*3);
  assert(rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_factor_en_gl(context, &en_gl_old[0][0][0], walk_num*elec_num*4);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_electron_coord(context, 'N', &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_set_single_point(context, 'N', elec, new_coords, 3*walk_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_single_en_gl(context, &delta_en_gl[0][0][0], walk_num*elec_num*4);
  assert (rc == QMCKL_SUCCESS);

  coords[0][elec][0] = new_coords[0];
  coords[0][elec][1] = new_coords[1];
  coords[0][elec][2] = new_coords[2];
  coords[1][elec][0] = new_coords[3];
  coords[1][elec][1] = new_coords[4];
  coords[1][elec][2] = new_coords[5];

	rc = qmckl_set_electron_coord(context, 'N', walk_num, &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_factor_en_gl(context, &en_gl_new[0][0][0], walk_num*elec_num*4);
  assert (rc == QMCKL_SUCCESS);

  for (int nw = 0; nw < walk_num; nw++) {
    for (int i = 0; i < elec_num; i++) {
      for (int m = 0; m < 4; m++) {
        assert(fabs((en_gl_new[nw][m][i] - en_gl_old[nw][m][i]) - delta_en_gl[nw][i][m]) < 1.e-12);
      }
    }
  }
}

printf("OK\n");

printf("Accept test 1\n");

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

for (int elec = 0; elec < elec_num; elec++){

	rc = qmckl_set_electron_coord(context, 'N', walk_num, elec_coord, walk_num*elec_num*3);
  assert(rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_factor_en(context, &jastrow_en_old[0], walk_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_factor_ee(context, &jastrow_ee_old[0], walk_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_factor_een(context, &jastrow_een_old[0], walk_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_factor_en_gl(context, &en_gl_old[0][0][0], walk_num*4*elec_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_factor_ee_gl(context, &ee_gl_old[0][0][0], walk_num*4*elec_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_factor_een_gl(context, &een_gl_old[0][0][0], walk_num*4*elec_num);
  assert (rc == QMCKL_SUCCESS);


  rc = qmckl_get_jastrow_champ_tmp_c(context, &p_old[0][0][0][0][0], walk_num*cord_num*(cord_num+1)*nucl_num*elec_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_dtmp_c(context,  &p_gl_old[0][0][0][0][0][0], 4*walk_num*cord_num*(cord_num+1)*nucl_num*elec_num);
  assert (rc == QMCKL_SUCCESS);



  // ----------------------------

  rc = qmckl_get_electron_coord(context, 'N', &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_set_single_point(context, 'N', elec, new_coords, 3*walk_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_single_en(context, &delta_en[0], walk_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_single_ee(context, &delta_ee[0], walk_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_single_een(context, &delta_een[0], walk_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_single_en_gl(context, &delta_en_gl[0][0][0], walk_num*4*elec_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_single_ee_gl(context, &delta_ee_gl[0][0][0], walk_num*4*elec_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_single_een_gl(context, &delta_een_gl[0][0][0], walk_num*4*elec_num);
  assert (rc == QMCKL_SUCCESS);



  rc = qmckl_get_een_rescaled_single_e(context, &single_rescaled_een_ee_distance[0][0][0], walk_num*(cord_num+1)*elec_num);
  assert (rc == QMCKL_SUCCESS);
  rc = qmckl_get_een_rescaled_single_n(context, &single_rescaled_een_en_distance[0][0][0], walk_num*(cord_num+1)*nucl_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_delta_p(context, &delta_p[0][0][0][0][0], walk_num*cord_num*(cord_num+1)*nucl_num*elec_num);
  assert (rc == QMCKL_SUCCESS);


  rc = qmckl_get_ee_rescaled_single(context, &single_ee_rescaled[0][0], walk_num*elec_num);
  assert (rc == QMCKL_SUCCESS);
  rc = qmckl_get_en_rescaled_single(context, &single_en_rescaled[0][0], walk_num*nucl_num);
  assert (rc == QMCKL_SUCCESS);


  rc = qmckl_get_een_rescaled_single_n_gl(context, &een_rescaled_single_n_gl[0][0][0][0], walk_num*(cord_num+1)*nucl_num*4);
  assert (rc == QMCKL_SUCCESS);
  rc = qmckl_get_een_rescaled_single_e_gl(context, &een_rescaled_single_e_gl[0][0][0][0], walk_num*(cord_num+1)*elec_num*4);
  assert (rc == QMCKL_SUCCESS);


  rc = qmckl_get_jastrow_champ_delta_p_gl(context, &delta_p_gl[0][0][0][0][0][0], 4*walk_num*cord_num*(cord_num+1)*nucl_num*elec_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_ee_rescaled_single_gl(context, &single_ee_rescaled_gl[0][0][0], walk_num*elec_num*4);
  assert (rc == QMCKL_SUCCESS);
  rc = qmckl_get_en_rescaled_single_gl(context, &single_en_rescaled_gl[0][0][0], walk_num*nucl_num*4);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_single_electron_ee_distance(context, &single_ee_distance[0][0], walk_num*elec_num);
  assert (rc == QMCKL_SUCCESS);
  rc = qmckl_get_single_electron_en_distance(context, &single_en_distance[0][0], walk_num*nucl_num);
  assert (rc == QMCKL_SUCCESS);

  // ----------------------------

  coords[0][elec][0] = new_coords[0];
  coords[0][elec][1] = new_coords[1];
  coords[0][elec][2] = new_coords[2];
  coords[1][elec][0] = new_coords[3];
  coords[1][elec][1] = new_coords[4];
  coords[1][elec][2] = new_coords[5];

  rc = qmckl_get_jastrow_champ_single_accept(context);
  assert (rc == QMCKL_SUCCESS);

  //rc = qmckl_context_touch(context);

  rc = qmckl_get_jastrow_champ_factor_en(context, &jastrow_en_new[0], walk_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_factor_ee(context, &jastrow_ee_new[0], walk_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_factor_een(context, &jastrow_een_new[0], walk_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_factor_en_gl(context, &en_gl_new[0][0][0], walk_num*4*elec_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_factor_ee_gl(context, &ee_gl_new[0][0][0], walk_num*4*elec_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_factor_een_gl(context, &een_gl_new[0][0][0], walk_num*4*elec_num);
  assert (rc == QMCKL_SUCCESS);



  rc = qmckl_get_jastrow_champ_een_rescaled_e(context, &rescaled_een_ee_distance[0][0][0][0], walk_num*(cord_num+1)*elec_num*elec_num);
  assert (rc == QMCKL_SUCCESS);
  rc = qmckl_get_jastrow_champ_een_rescaled_n(context, &rescaled_een_en_distance[0][0][0][0], walk_num*(cord_num+1)*elec_num*nucl_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_tmp_c(context,  &p_new[0][0][0][0][0], walk_num*cord_num*(cord_num+1)*nucl_num*elec_num);
  assert (rc == QMCKL_SUCCESS);


  rc = qmckl_get_jastrow_champ_ee_distance_rescaled(context, &ee_rescaled[0][0][0], walk_num*elec_num*elec_num);
  assert (rc == QMCKL_SUCCESS);
  rc = qmckl_get_jastrow_champ_en_distance_rescaled(context, &en_rescaled[0][0][0], walk_num*nucl_num*elec_num);
  assert (rc == QMCKL_SUCCESS);


  rc = qmckl_get_jastrow_champ_een_rescaled_n_gl(context, &een_rescaled_en_gl[0][0][0][0][0], walk_num*(cord_num+1)*nucl_num*elec_num*4);
  assert (rc == QMCKL_SUCCESS);
  rc = qmckl_get_jastrow_champ_een_rescaled_e_gl(context, &een_rescaled_ee_gl[0][0][0][0][0], walk_num*(cord_num+1)*elec_num*elec_num*4);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_dtmp_c(context,  &p_gl_new[0][0][0][0][0][0], 4*walk_num*cord_num*(cord_num+1)*nucl_num*elec_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_ee_distance_rescaled_gl(context, &ee_rescaled_gl[0][0][0][0], walk_num*4*elec_num*elec_num);
  assert (rc == QMCKL_SUCCESS);
  rc = qmckl_get_jastrow_champ_en_distance_rescaled_gl(context, &en_rescaled_gl[0][0][0][0], walk_num*4*nucl_num*elec_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_electron_ee_distance(context, &ee_distance[0][0][0], walk_num*elec_num*elec_num);
  assert (rc == QMCKL_SUCCESS);
  rc = qmckl_get_electron_en_distance(context, &en_distance[0][0][0], walk_num*nucl_num*elec_num);
  assert (rc == QMCKL_SUCCESS);

  double metric[4] = {-1.0, -1.0, -1.0, 1.0};

  for (int nw = 0; nw < walk_num; nw++) {
    //printf("jastrow_en_new[%d] = %f\n", nw, jastrow_en_new[nw]);
    //printf("jastrow_en_old[%d] = %f\n", nw, jastrow_en_old[nw]);
    //printf("delta_en[%d] = %f\n", nw, delta_en[nw]);
    assert(fabs((jastrow_en_new[nw] - jastrow_en_old[nw]) - delta_en[nw]) < 1.e-12);
    assert(fabs((jastrow_ee_new[nw] - jastrow_ee_old[nw]) - delta_ee[nw]) < 1.e-12);
    assert(fabs((jastrow_een_new[nw] - jastrow_een_old[nw]) - delta_een[nw]) < 1.e-12);
    for (int i = 0; i < elec_num; i++){
      for (int k = 0; k < 4; k++){
        //printf("en_gl_new[%d][%d][%d] = %f\n", nw, k, i, en_gl_new[nw][k][i]);
        //printf("en_gl_old[%d][%d][%d] = %f\n", nw, k, i, en_gl_old[nw][k][i]);
        //printf("delta_en_gl[%d][%d][%d] = %f\n", nw, i, k, delta_en_gl[nw][i][k]);
        assert(fabs((en_gl_new[nw][k][i] - en_gl_old[nw][k][i]) - delta_en_gl[nw][i][k]) < 1.e-12);
        assert(fabs((ee_gl_new[nw][k][i] - ee_gl_old[nw][k][i]) - delta_ee_gl[nw][i][k]) < 1.e-12);
        assert(fabs((een_gl_new[nw][k][i] - een_gl_old[nw][k][i]) - delta_een_gl[nw][k][i]) < 1.e-12);
      }
    }
  }

  for (int nw = 0; nw < walk_num; nw++){
    for (int l = 0; l <= cord_num; l++){
      for (int i = 0; i < elec_num; i++) {
        //printf("rescaled_een_ee_distance[%d][%d][elec][%d] = %f\n", nw, l, i, rescaled_een_ee_distance[nw][l][elec][i]);
        //printf("single_rescaled_een_ee_distance[%d][%d][%d] = %f\n", nw, l, i, single_rescaled_een_ee_distance[nw][l][i]);
        assert(fabs((rescaled_een_ee_distance[nw][l][elec][i]-single_rescaled_een_ee_distance[nw][l][i])) < 1.e-12);
        assert(fabs((rescaled_een_ee_distance[nw][l][i][elec]-single_rescaled_een_ee_distance[nw][l][i])) < 1.e-12);
      }
    }
  }

  for (int nw = 0; nw < walk_num; nw++){
    for (int l = 0; l <= cord_num; l++){
      for (int a = 0; a < nucl_num; a++) {
        //printf("rescaled_een_en_distance[%d][%d][%d][elec] = %f\n", nw, l, a, rescaled_een_en_distance[nw][l][a][elec]);
        //printf("single_rescaled_een_en_distance[%d][%d][%d] = %f\n", nw, l, a, single_rescaled_een_en_distance[nw][l][a]);
        assert(fabs((rescaled_een_en_distance[nw][l][a][elec]-single_rescaled_een_en_distance[nw][l][a])) < 1.e-12);
      }
    }
  }
  for (int nw = 0; nw < walk_num; nw++){
    for (int l = 0; l < cord_num; l++){
      for (int m = 0; m <= cord_num; m++){
        for (int a = 0; a < nucl_num; a++) {
          for (int i = 0; i < elec_num; i++){
            assert(fabs(((p_new[nw][l][m][a][i]-p_old[nw][l][m][a][i])-delta_p[nw][l][m][a][i])) < 1.e-12);
          }
        }
      }
    }
  }

  for (int nw = 0; nw < walk_num; nw++) {
    for (int i = 0; i < elec_num; i++){
      assert(fabs(ee_rescaled[nw][elec][i]-single_ee_rescaled[nw][i]) < 1.e-12);
    }
  }
  for (int nw = 0; nw < walk_num; nw++) {
    for (int a = 0; a < nucl_num; a++){
      assert(fabs(en_rescaled[nw][a][elec]-single_en_rescaled[nw][a]) < 1.e-12);
    }
  }

  for (int l = 0; l < cord_num+1; l++) {
    for (int nw = 0; nw < walk_num; nw++) {
      for (int a = 0; a < nucl_num; a++) {
        for (int m = 0; m < 4; m++) {
          assert(fabs(een_rescaled_en_gl[nw][l][a][m][elec] - een_rescaled_single_n_gl[nw][l][a][m]) < 1.e-12);
        }
      }
    }
  }

  for (int l = 0; l < cord_num+1; l++) {
    for (int nw = 0; nw < walk_num; nw++) {
      for (int i = 0; i < elec_num; i++) {
        for (int m = 0; m < 4; m++) {
          //printf("een_rescaled_ee_gl[nw][l][i][m][elec] %i %i %i %f \n", l, m ,i, een_rescaled_ee_gl[nw][l][i][m][elec]);
          //printf("een_rescaled_single_e_gl[nw][l][i][m] %i %i %i %f\n", l, m, i,een_rescaled_single_e_gl[nw][l][i][m]);
          assert(fabs(een_rescaled_ee_gl[nw][l][i][m][elec] - een_rescaled_single_e_gl[nw][l][i][m]) < 1.e-12);
          assert(fabs(een_rescaled_ee_gl[nw][l][elec][m][i] - metric[m] * een_rescaled_single_e_gl[nw][l][i][m]) < 1.e-12);
        }
      }
    }
  }

  for (int nw = 0; nw < walk_num; nw++){
    for (int l = 0; l < cord_num; l++){
      for (int m = 0; m <= cord_num; m++){
        for (int a = 0; a < nucl_num; a++) {
          for (int i = 0; i < elec_num; i++){
            for (int k = 0; k < 4; k++){
              //printf("p_gl[%d][%d][%d][%d][%d][%d] = %f\n", nw, l, m, a, k, i, p_gl_new[nw][l][m][a][k][i] - p_gl_old[nw][l][m][a][k][i]);
              //printf("delta_p_gl[%d][%d][%d][%d][%d][%d] = %f\n", nw, l, m, a, k, i, delta_p_gl[nw][l][m][k][a][i]);
              assert(fabs(((p_gl_new[nw][l][m][a][k][i]-p_gl_old[nw][l][m][a][k][i])-delta_p_gl[nw][l][m][k][a][i])) < 1.e-12);
            }
          }
        }
      }
    }
  }

  for (int nw = 0; nw < walk_num; nw++) {
    for (int i = 0; i < elec_num; i++) {
      for (int m = 0; m < 4; m++) {
        if (i == 2) continue;
        //printf("%f\n", ee_rescaled_gl[nw][elec][i][m]);
        //printf("%f\n", single_ee_rescaled_gl[nw][i][m]);
        assert(fabs(ee_rescaled_gl[nw][elec][i][m] - single_ee_rescaled_gl[nw][i][m]) < 1.e-12);
        assert(fabs(ee_rescaled_gl[nw][i][elec][m] - metric[m] * single_ee_rescaled_gl[nw][i][m]) < 1.e-12);     
      }
    }
  }

  for (int nw = 0; nw < walk_num; nw++) {
    for (int a = 0; a < nucl_num; a++) {
        for (int m = 0; m < 4; m++) {
          assert(fabs(en_rescaled_gl[nw][a][elec][m] - single_en_rescaled_gl[nw][a][m]) < 1.e-12);
        
        }
      }
    }

  for (int nw = 0; nw < walk_num; nw++){
    for (int i = 0; i < elec_num; i++) {
      if (i == 2) continue;
      assert(fabs((ee_distance[nw][elec][i]-single_ee_distance[nw][i])) < 1.e-12);
    }
  }

  for (int nw = 0; nw < walk_num; nw++){
    for (int a = 0; a < nucl_num; a++){
          assert(fabs((en_distance[nw][elec][a]-single_en_distance[nw][a])) < 1.e-12);
    }
  }
}

printf("OK\n");

printf("Accept test 2\n");

/* Check if Jastrow is properly initialized */
assert(qmckl_jastrow_champ_provided(context));

for (int elec = 0; elec < elec_num; elec++){

	rc = qmckl_set_electron_coord(context, 'N', walk_num, elec_coord, walk_num*elec_num*3);
  assert(rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_factor_en(context, &jastrow_en_old[0], walk_num);
  assert (rc == QMCKL_SUCCESS);
  rc = qmckl_get_jastrow_champ_factor_ee(context, &jastrow_ee_old[0], walk_num);
  assert (rc == QMCKL_SUCCESS);
  rc = qmckl_get_jastrow_champ_factor_een(context, &jastrow_een_old[0], walk_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_factor_en_gl(context, &en_gl_old[0][0][0], walk_num*4*elec_num);
  assert (rc == QMCKL_SUCCESS);
  rc = qmckl_get_jastrow_champ_factor_ee_gl(context, &ee_gl_old[0][0][0], walk_num*4*elec_num);
  assert (rc == QMCKL_SUCCESS);
  rc = qmckl_get_jastrow_champ_factor_een_gl(context, &een_gl_old[0][0][0], walk_num*4*elec_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_electron_ee_distance(context, &ee_distance[0][0][0], walk_num*elec_num*elec_num);
  assert (rc == QMCKL_SUCCESS);
  rc = qmckl_get_electron_en_distance(context, &en_distance[0][0][0], walk_num*nucl_num*elec_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_ee_distance_rescaled(context, &ee_rescaled[0][0][0], walk_num*elec_num*elec_num);
  assert (rc == QMCKL_SUCCESS);
  rc = qmckl_get_jastrow_champ_ee_distance_rescaled_gl(context, &ee_rescaled_gl[0][0][0][0], walk_num*4*elec_num*elec_num);
  assert (rc == QMCKL_SUCCESS);


  rc = qmckl_get_electron_coord(context, 'N', &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_set_single_point(context, 'N', 2, new_coords, 3*walk_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_single_en(context, &delta_en[0], walk_num);
  assert (rc == QMCKL_SUCCESS);
  rc = qmckl_get_jastrow_champ_single_ee(context, &delta_ee[0], walk_num);
  assert (rc == QMCKL_SUCCESS);
  rc = qmckl_get_jastrow_champ_single_een(context, &delta_een[0], walk_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_single_en_gl(context, &delta_en_gl[0][0][0], walk_num*4*elec_num);
  assert (rc == QMCKL_SUCCESS);
  rc = qmckl_get_jastrow_champ_single_ee_gl(context, &delta_ee_gl[0][0][0], walk_num*4*elec_num);
  assert (rc == QMCKL_SUCCESS);
  rc = qmckl_get_jastrow_champ_single_een_gl(context, &delta_een_gl[0][0][0], walk_num*4*elec_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_single_electron_ee_distance(context, &single_ee_distance[0][0], walk_num*elec_num);
  assert (rc == QMCKL_SUCCESS);
  rc = qmckl_get_single_electron_en_distance(context, &single_en_distance[0][0], walk_num*nucl_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_ee_rescaled_single(context, &single_ee_rescaled[0][0], walk_num*elec_num);
  assert (rc == QMCKL_SUCCESS);
  rc = qmckl_get_ee_rescaled_single_gl(context, &single_ee_rescaled_gl[0][0][0], walk_num*4*elec_num);
  assert (rc == QMCKL_SUCCESS);

  coords[0][elec][0] = new_coords[0];
  coords[0][elec][1] = new_coords[1];
  coords[0][elec][2] = new_coords[2];
  coords[1][elec][0] = new_coords[3];
  coords[1][elec][1] = new_coords[4];
  coords[1][elec][2] = new_coords[5];

  rc = qmckl_get_jastrow_champ_single_accept(context);
  assert (rc == QMCKL_SUCCESS);

  //rc = qmckl_context_touch(context);

  rc = qmckl_get_jastrow_champ_factor_en(context, &jastrow_en_old[0], walk_num);
  assert (rc == QMCKL_SUCCESS);
  rc = qmckl_get_jastrow_champ_factor_ee(context, &jastrow_ee_old[0], walk_num);
  assert (rc == QMCKL_SUCCESS);
  rc = qmckl_get_jastrow_champ_factor_een(context, &jastrow_een_old[0], walk_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_factor_en_gl(context, &en_gl_old[0][0][0], walk_num*4*elec_num);
  assert (rc == QMCKL_SUCCESS);
  rc = qmckl_get_jastrow_champ_factor_ee_gl(context, &ee_gl_old[0][0][0], walk_num*4*elec_num);
  assert (rc == QMCKL_SUCCESS);
  rc = qmckl_get_jastrow_champ_factor_een_gl(context, &een_gl_old[0][0][0], walk_num*4*elec_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_electron_ee_distance(context, &ee_distance[0][0][0], walk_num*elec_num*elec_num);
  assert (rc == QMCKL_SUCCESS);
  rc = qmckl_get_electron_en_distance(context, &en_distance[0][0][0], walk_num*nucl_num*elec_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_ee_distance_rescaled(context, &ee_rescaled[0][0][0], walk_num*elec_num*elec_num);
  assert (rc == QMCKL_SUCCESS);
  rc = qmckl_get_jastrow_champ_ee_distance_rescaled_gl(context, &ee_rescaled_gl[0][0][0][0], walk_num*4*elec_num*elec_num);
  assert (rc == QMCKL_SUCCESS);


  double new_coords2[6] = {3.0,2.0,1.0,3.0,2.0,1.0};


  rc = qmckl_set_single_point(context, 'N', elec, new_coords2, 3*walk_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_single_en(context, &delta_en[0], walk_num);
  assert (rc == QMCKL_SUCCESS);
  rc = qmckl_get_jastrow_champ_single_ee(context, &delta_ee[0], walk_num);
  assert (rc == QMCKL_SUCCESS);
  rc = qmckl_get_jastrow_champ_single_een(context, &delta_een[0], walk_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_single_en_gl(context, &delta_en_gl[0][0][0], walk_num*4*elec_num);
  assert (rc == QMCKL_SUCCESS);
  rc = qmckl_get_jastrow_champ_single_ee_gl(context, &delta_ee_gl[0][0][0], walk_num*4*elec_num);
  assert (rc == QMCKL_SUCCESS);
  rc = qmckl_get_jastrow_champ_single_een_gl(context, &delta_een_gl[0][0][0], walk_num*4*elec_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_single_electron_ee_distance(context, &single_ee_distance[0][0], walk_num*elec_num);
  assert (rc == QMCKL_SUCCESS);
  rc = qmckl_get_single_electron_en_distance(context, &single_en_distance[0][0], walk_num*nucl_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_ee_rescaled_single(context, &single_ee_rescaled[0][0], walk_num*elec_num);
  assert (rc == QMCKL_SUCCESS);
  rc = qmckl_get_ee_rescaled_single_gl(context, &single_ee_rescaled_gl[0][0][0], walk_num*4*elec_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_electron_coord(context, 'N', &coords[0][0][0], walk_num*elec_num*3);
  assert (rc == QMCKL_SUCCESS);

  coords[0][elec][0] = new_coords2[0];
  coords[0][elec][1] = new_coords2[1];
  coords[0][elec][2] = new_coords2[2];
  coords[1][elec][0] = new_coords2[3];
  coords[1][elec][1] = new_coords2[4];
  coords[1][elec][2] = new_coords2[5];

	rc = qmckl_set_electron_coord(context, 'N', walk_num, &coords[0][0][0], walk_num*elec_num*3);
  rc = qmckl_context_touch(context);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_factor_en(context, &jastrow_en_new[0], walk_num);
  assert (rc == QMCKL_SUCCESS);
  rc = qmckl_get_jastrow_champ_factor_ee(context, &jastrow_ee_new[0], walk_num);
  assert (rc == QMCKL_SUCCESS);
  rc = qmckl_get_jastrow_champ_factor_een(context, &jastrow_een_new[0], walk_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_factor_en_gl(context, &en_gl_new[0][0][0], walk_num*4*elec_num);
  assert (rc == QMCKL_SUCCESS);
  rc = qmckl_get_jastrow_champ_factor_ee_gl(context, &ee_gl_new[0][0][0], walk_num*4*elec_num);
  assert (rc == QMCKL_SUCCESS);
  rc = qmckl_get_jastrow_champ_factor_een_gl(context, &een_gl_new[0][0][0], walk_num*4*elec_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_electron_ee_distance(context, &ee_distance[0][0][0], walk_num*elec_num*elec_num);
  assert (rc == QMCKL_SUCCESS);
  rc = qmckl_get_electron_en_distance(context, &en_distance[0][0][0], walk_num*nucl_num*elec_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_jastrow_champ_ee_distance_rescaled(context, &ee_rescaled[0][0][0], walk_num*elec_num*elec_num);
  assert (rc == QMCKL_SUCCESS);
  rc = qmckl_get_jastrow_champ_ee_distance_rescaled_gl(context, &ee_rescaled_gl[0][0][0][0], walk_num*4*elec_num*elec_num);
  assert (rc == QMCKL_SUCCESS);

  for (int nw = 0; nw < walk_num; nw++) {
    //printf("jastrow_en_new[%d] = %f\t", nw, jastrow_en_new[nw]);
    //printf("jastrow_en_old[%d] = %f\t", nw, jastrow_en_old[nw]);
    //printf("delta_en[%d] = %f\n", nw, delta_en[nw]);
    //printf("jastrow_ee_new[%d] = %f\t", nw, jastrow_ee_new[nw]);
    //printf("jastrow_ee_old[%d] = %f\t", nw, jastrow_ee_old[nw]);
    //printf("delta_ee[%d] = %f\n", nw, delta_ee[nw]);
    assert(fabs((jastrow_en_new[nw] - jastrow_en_old[nw]) - delta_en[nw]) < 1.e-12);
    assert(fabs((jastrow_ee_new[nw] - jastrow_ee_old[nw]) - delta_ee[nw]) < 1.e-12);
    assert(fabs((jastrow_een_new[nw] - jastrow_een_old[nw]) - delta_een[nw]) < 1.e-12);
    for (int i = 0; i < elec_num; i++){
      for (int k = 0; k < 4; k++){
        //printf("ee_gl_new[%d][%d][%d] = %f\n", nw, k, i, ee_gl_new[nw][k][i]);
        //printf("ee_gl_old[%d][%d][%d] = %f\n", nw, k, i, ee_gl_old[nw][k][i]);
        //printf("delta_ee_gl[%d][%d][%d] = %f\n", nw, i, k, delta_ee_gl[nw][i][k]);
        assert(fabs((en_gl_new[nw][k][i] - en_gl_old[nw][k][i]) - delta_en_gl[nw][i][k]) < 1.e-12);
        assert(fabs((ee_gl_new[nw][k][i] - ee_gl_old[nw][k][i]) - delta_ee_gl[nw][i][k]) < 1.e-12);
        assert(fabs((een_gl_new[nw][k][i] - een_gl_old[nw][k][i]) - delta_een_gl[nw][k][i]) < 1.e-12);
      }
    }
  }

  for (int nw = 0; nw < walk_num; nw++){
    for (int i = 0; i < elec_num; i++) {
      if (i == 1) continue;
      //printf("ee_distance[%d][elec][%d] = %f\n", nw, i, ee_distance[nw][elec][i]);
      //printf("single_ee_distance[%d][%d] = %f\n", nw, i, single_ee_distance[nw][i]);
      assert(fabs((ee_distance[nw][elec][i]-single_ee_distance[nw][i])) < 1.e-12);
    }
  }

  for (int nw = 0; nw < walk_num; nw++) {
    for (int a = 0; a < nucl_num; a++){
          assert(fabs((en_distance[nw][elec][a]-single_en_distance[nw][a])) < 1.e-12);
    }
  }
  for (int nw = 0; nw < walk_num; nw++) {
    for (int i = 0; i < elec_num; i++){
      assert(fabs(ee_rescaled[nw][elec][i]-single_ee_rescaled[nw][i]) < 1.e-12);
    }
  }

  for (int nw = 0; nw < walk_num; nw++) {
    for (int i = 0; i < elec_num; i++) {
      for (int m = 0; m < 4; m++) {
        if (i == elec) continue;
        //printf("%f\n", ee_rescaled_gl[nw][elec][i][m]);
        //printf("%f\n", single_ee_rescaled_gl[nw][i][m]);
        assert(fabs(ee_rescaled_gl[nw][elec][i][m] - single_ee_rescaled_gl[nw][i][m]) < 1.e-12);
      
      }
    }
  }

}

printf("OK\n");

rc = qmckl_context_destroy(context);
    assert (rc == QMCKL_SUCCESS);

    return 0;
}
