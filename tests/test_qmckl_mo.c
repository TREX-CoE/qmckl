#include "qmckl.h"
#include "assert.h"
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <math.h>
#include "chbrclf.h"
#include "qmckl_electron_private_func.h"
#include "qmckl_ao_private_func.h"
#include "qmckl_mo_private_func.h"

int main() {
    qmckl_context context;
    context = qmckl_context_create();

    qmckl_exit_code rc;

{
#define walk_num chbrclf_walk_num
#define elec_num chbrclf_elec_num
#define shell_num chbrclf_shell_num
#define ao_num chbrclf_ao_num

  int64_t elec_up_num   = chbrclf_elec_up_num;
  int64_t elec_dn_num   = chbrclf_elec_dn_num;
  double* elec_coord    = &(chbrclf_elec_coord[0][0][0]);
  int64_t   nucl_num      = chbrclf_nucl_num;
  const double*   nucl_charge   = chbrclf_charge;
  const double*   nucl_coord    = &(chbrclf_nucl_coord[0][0]);

  int64_t point_num = walk_num*elec_num;

  rc = qmckl_set_electron_num (context, elec_up_num, elec_dn_num);
  assert (rc == QMCKL_SUCCESS);

  assert(qmckl_electron_provided(context));

  rc = qmckl_set_point(context, 'N', point_num, elec_coord, point_num*3);
  assert(rc == QMCKL_SUCCESS);

  rc = qmckl_set_nucleus_num (context, nucl_num);
  assert(rc == QMCKL_SUCCESS);

  rc = qmckl_set_nucleus_coord (context, 'T', &(nucl_coord[0]), nucl_num*3);
  assert(rc == QMCKL_SUCCESS);

  rc = qmckl_set_nucleus_charge(context, nucl_charge, nucl_num);
  assert(rc == QMCKL_SUCCESS);

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
  assert(rc == QMCKL_SUCCESS);
  assert(!qmckl_ao_basis_provided(context));

  rc = qmckl_set_ao_basis_shell_num (context, chbrclf_shell_num);
  assert(rc == QMCKL_SUCCESS);
  assert(!qmckl_ao_basis_provided(context));

  rc = qmckl_set_ao_basis_prim_num (context, chbrclf_prim_num);
  assert(rc == QMCKL_SUCCESS);
  assert(!qmckl_ao_basis_provided(context));

  rc = qmckl_set_ao_basis_nucleus_index (context, nucleus_index, nucl_num);
  assert(rc == QMCKL_SUCCESS);
  assert(!qmckl_ao_basis_provided(context));

  rc = qmckl_set_ao_basis_nucleus_shell_num (context, nucleus_shell_num, nucl_num);
  assert(rc == QMCKL_SUCCESS);
  assert(!qmckl_ao_basis_provided(context));

  rc = qmckl_set_ao_basis_shell_ang_mom (context, shell_ang_mom, chbrclf_shell_num);
  assert(rc == QMCKL_SUCCESS);
  assert(!qmckl_ao_basis_provided(context));

  rc = qmckl_set_ao_basis_shell_factor  (context, shell_factor, chbrclf_shell_num);
  assert(rc == QMCKL_SUCCESS);
  assert(!qmckl_ao_basis_provided(context));

  rc = qmckl_set_ao_basis_shell_prim_num (context, shell_prim_num, chbrclf_shell_num);
  assert(rc == QMCKL_SUCCESS);
  assert(!qmckl_ao_basis_provided(context));

  rc = qmckl_set_ao_basis_shell_prim_index (context, shell_prim_index, chbrclf_shell_num);
  assert(rc == QMCKL_SUCCESS);
  assert(!qmckl_ao_basis_provided(context));

  rc = qmckl_set_ao_basis_exponent      (context, exponent, chbrclf_prim_num);
  assert(rc == QMCKL_SUCCESS);
  assert(!qmckl_ao_basis_provided(context));

  rc = qmckl_set_ao_basis_coefficient   (context, coefficient, chbrclf_prim_num);
  assert(rc == QMCKL_SUCCESS);
  assert(!qmckl_ao_basis_provided(context));

  rc = qmckl_set_ao_basis_prim_factor (context, prim_factor, chbrclf_prim_num);
  assert(rc == QMCKL_SUCCESS);

  rc = qmckl_set_ao_basis_ao_num(context, chbrclf_ao_num);
  assert(rc == QMCKL_SUCCESS);

  rc = qmckl_set_ao_basis_ao_factor (context, ao_factor, chbrclf_ao_num);
  assert(rc == QMCKL_SUCCESS);

  assert(qmckl_ao_basis_provided(context));


  double ao_vgl[point_num][5][chbrclf_ao_num];

  rc = qmckl_get_ao_basis_ao_vgl(context, &(ao_vgl[0][0][0]),
                                 (int64_t) 5*point_num*chbrclf_ao_num);
  assert (rc == QMCKL_SUCCESS);

  /* Set up MO data */
  int64_t mo_num = chbrclf_mo_num;
  rc = qmckl_set_mo_basis_mo_num(context, mo_num);
  assert (rc == QMCKL_SUCCESS);

  const double  * mo_coefficient          =  &(chbrclf_mo_coef[0]);

  rc = qmckl_set_mo_basis_coefficient(context, mo_coefficient, chbrclf_mo_num*chbrclf_ao_num);
  assert (rc == QMCKL_SUCCESS);

  assert(qmckl_mo_basis_provided(context));

  rc = qmckl_context_touch(context);
  assert (rc == QMCKL_SUCCESS);

  double mo_value[point_num][chbrclf_mo_num];
  rc = qmckl_get_mo_basis_mo_value(context, &(mo_value[0][0]), point_num*chbrclf_mo_num);
  assert (rc == QMCKL_SUCCESS);

  double mo_vgl[point_num][5][chbrclf_mo_num];
  rc = qmckl_get_mo_basis_mo_vgl(context, &(mo_vgl[0][0][0]), point_num*5*chbrclf_mo_num);
  assert (rc == QMCKL_SUCCESS);

  for (int i=0 ; i< point_num; ++i) {
    for (int k=0 ; k< chbrclf_mo_num ; ++k) {
      assert(fabs(mo_vgl[i][0][k] - mo_value[i][k]) < 1.e-12) ;
    }
  }

  rc = qmckl_context_touch(context);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_mo_basis_mo_value(context, &(mo_value[0][0]), point_num*chbrclf_mo_num);
  assert (rc == QMCKL_SUCCESS);

  for (int i=0 ; i< point_num; ++i) {
    for (int k=0 ; k< chbrclf_mo_num ; ++k) {
      assert(fabs(mo_vgl[i][0][k] - mo_value[i][k]) < 1.e-12) ;
    }
  }

  rc = qmckl_mo_basis_rescale(context, 0.);
  assert (rc != QMCKL_SUCCESS);

  rc = qmckl_mo_basis_rescale(context, 2.);
  assert (rc == QMCKL_SUCCESS);


  rc = qmckl_get_mo_basis_mo_value(context, &(mo_value[0][0]), point_num*chbrclf_mo_num);
  assert (rc == QMCKL_SUCCESS);

  for (int i=0 ; i< point_num; ++i) {
    for (int k=0 ; k< chbrclf_mo_num ; ++k) {
      assert(fabs(2.*mo_vgl[i][0][k] - mo_value[i][k]) < 1.e-12) ;
    }
  }

  rc = qmckl_mo_basis_rescale(context, 0.5);
  assert (rc == QMCKL_SUCCESS);


  printf("\n");
  printf(" mo_vgl mo_vgl[0][26][219] %25.15e\n", mo_vgl[2][0][3]);
  printf(" mo_vgl mo_vgl[1][26][219] %25.15e\n", mo_vgl[2][1][3]);
  printf(" mo_vgl mo_vgl[0][26][220] %25.15e\n", mo_vgl[2][0][3]);
  printf(" mo_vgl mo_vgl[1][26][220] %25.15e\n", mo_vgl[2][1][3]);
  printf(" mo_vgl mo_vgl[0][26][221] %25.15e\n", mo_vgl[2][0][3]);
  printf(" mo_vgl mo_vgl[1][26][221] %25.15e\n", mo_vgl[2][1][3]);
  printf(" mo_vgl mo_vgl[0][26][222] %25.15e\n", mo_vgl[2][0][3]);
  printf(" mo_vgl mo_vgl[1][26][222] %25.15e\n", mo_vgl[2][1][3]);
  printf(" mo_vgl mo_vgl[0][26][223] %25.15e\n", mo_vgl[2][0][3]);
  printf(" mo_vgl mo_vgl[1][26][223] %25.15e\n", mo_vgl[2][1][3]);
  printf(" mo_vgl mo_vgl[0][26][224] %25.15e\n", mo_vgl[2][0][3]);
  printf(" mo_vgl mo_vgl[1][26][224] %25.15e\n", mo_vgl[2][1][3]);
  printf("\n");

  /* Check single precision */
  rc = qmckl_set_numprec_precision(context,53);
  rc = qmckl_context_touch(context);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_mo_basis_mo_value(context, &(mo_value[0][0]), point_num*chbrclf_mo_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_mo_basis_mo_vgl(context, &(mo_vgl[0][0][0]), point_num*5*chbrclf_mo_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_set_numprec_precision(context,23);
  rc = qmckl_context_touch(context);
  assert (rc == QMCKL_SUCCESS);

  double mo_value_sp[point_num][chbrclf_mo_num];
  rc = qmckl_get_mo_basis_mo_value(context, &(mo_value_sp[0][0]), point_num*chbrclf_mo_num);
  assert (rc == QMCKL_SUCCESS);

  double mo_vgl_sp[point_num][5][chbrclf_mo_num];
  rc = qmckl_get_mo_basis_mo_vgl(context, &(mo_vgl_sp[0][0][0]), point_num*5*chbrclf_mo_num);
  assert (rc == QMCKL_SUCCESS);


  uint64_t average_prec = 0;
  int32_t nbits;
  for (int i=0 ; i< point_num; ++i) {
    for (int k=0 ; k< chbrclf_mo_num ; ++k) {
      nbits = qmckl_test_precision_64(mo_value_sp[i][k], mo_value[i][k]);
//      printf("%d %d %25.15e %25.15e %d\n", i, k, mo_value_sp[i][k], mo_value[i][k], nbits);
      average_prec += nbits;
    }
  }
  printf("Average precision for %d: %d\n", qmckl_get_numprec_precision(context),
            (int) (average_prec/(point_num*chbrclf_mo_num)));
  assert(nbits > 12);
  fflush(stdout);

  average_prec = 0;
  for (int i=0 ; i< point_num; ++i) {
    for (int k=0 ; k< chbrclf_mo_num ; ++k) {
//      printf("%d %d\n", i, k);
      nbits = qmckl_test_precision_64(mo_vgl_sp[i][0][k], mo_vgl[i][0][k]);
      average_prec += nbits;
//      printf("%25.15e %25.15e %d\n", mo_vgl_sp[i][0][k], mo_vgl[i][0][k], nbits);
      nbits = qmckl_test_precision_64(mo_vgl_sp[i][1][k], mo_vgl[i][1][k]);
      average_prec += nbits;
//      printf("%25.15e %25.15e %d\n", mo_vgl_sp[i][1][k], mo_vgl[i][1][k], nbits);
      nbits = qmckl_test_precision_64(mo_vgl_sp[i][2][k], mo_vgl[i][2][k]);
      average_prec += nbits;
//      printf("%25.15e %25.15e %d\n", mo_vgl_sp[i][2][k], mo_vgl[i][2][k], nbits);
      nbits = qmckl_test_precision_64(mo_vgl_sp[i][3][k], mo_vgl[i][3][k]);
      average_prec += nbits;
//      printf("%25.15e %25.15e %d\n", mo_vgl_sp[i][3][k], mo_vgl[i][3][k], nbits);
      nbits = qmckl_test_precision_64(mo_vgl_sp[i][4][k], mo_vgl[i][4][k]);
      average_prec += nbits;
//      printf("%25.15e %25.15e %d\n", mo_vgl_sp[i][4][k], mo_vgl[i][4][k], nbits);
    }
  }
  nbits = (int) (average_prec/(point_num*chbrclf_mo_num*5));
  printf("Average precision for %d: %d\n", qmckl_get_numprec_precision(context), nbits);
  assert(nbits > 11);
  fflush(stdout);
  rc = qmckl_set_numprec_precision(context,53);


  /* Check selection of MOs */

  int32_t keep[mo_num];
  for (int i=0 ; i<mo_num ; ++i) {
    keep[i] = 0;
  }
  keep[2] = 1;
  keep[5] = 1;

  rc = qmckl_mo_basis_select_mo(context, &(keep[0]), mo_num);
  assert (rc == QMCKL_SUCCESS);

  rc = qmckl_get_mo_basis_mo_num(context, &mo_num);
  printf(" mo_num: %ld\n", (long) mo_num);
  assert(mo_num == 2);

  double mo_coefficient_new[mo_num][ao_num];
  rc = qmckl_get_mo_basis_coefficient (context, &(mo_coefficient_new[0][0]), mo_num*ao_num);
  for (int i=0 ; i<ao_num ; ++i) {
    assert(mo_coefficient_new[0][i] == mo_coefficient[i + ao_num*2]);
    assert(mo_coefficient_new[1][i] == mo_coefficient[i + ao_num*5]);
  }


}

rc = qmckl_context_destroy(context);
    assert (rc == QMCKL_SUCCESS);

    return 0;
}
