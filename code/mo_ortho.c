#include <qmckl.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>

int main(int argc, char** argv)
{
  const char* trexio_filename = "..//share/qmckl/test_data/h2o_5z.h5";
  qmckl_exit_code rc = QMCKL_SUCCESS;

qmckl_context context = qmckl_context_create();

rc = qmckl_trexio_read(context, trexio_filename, strlen(trexio_filename));

if (rc != QMCKL_SUCCESS) {
  fprintf(stderr, "Error reading TREXIO file:\n%s\n", qmckl_string_of_error(rc));
  exit(1);
}

int64_t nucl_num;

rc = qmckl_get_nucleus_num(context, &nucl_num);

if (rc != QMCKL_SUCCESS) {
  fprintf(stderr, "Error getting nucl_num:\n%s\n", qmckl_string_of_error(rc));
  exit(1);
}


double nucl_charge[nucl_num];

rc = qmckl_get_nucleus_charge(context, &(nucl_charge[0]), nucl_num);

if (rc != QMCKL_SUCCESS) {
  fprintf(stderr, "Error getting nucl_charge:\n%s\n", qmckl_string_of_error(rc));
  exit(1);
}


double r_cusp[nucl_num];

for (size_t i=0 ; i<nucl_num ; ++i) {

  switch ((int) nucl_charge[i]) {

  case 1:
    r_cusp[i] = 0.5;
    break;

  case 8:
    r_cusp[i] = 0.1;
    break;
  }

}


rc = qmckl_set_mo_basis_r_cusp(context, &(r_cusp[0]), nucl_num);

if (rc != QMCKL_SUCCESS) {
  fprintf(stderr, "Error setting r_cusp:\n%s\n", qmckl_string_of_error(rc));
  exit(1);
}

double nucl_coord[nucl_num][3];

rc = qmckl_get_nucleus_coord(context, 'N', &(nucl_coord[0][0]), nucl_num*3);

if (rc != QMCKL_SUCCESS) {
  fprintf(stderr, "Error getting nucl_coord:\n%s\n", qmckl_string_of_error(rc));
  exit(1);
}

for (size_t i=0 ; i<nucl_num ; ++i) {
  printf("%d  %+f %+f %+f\n",
         (int32_t) nucl_charge[i],
         nucl_coord[i][0],
         nucl_coord[i][1],
         nucl_coord[i][2]);
}

size_t nx[3] = { 120, 120, 120 };
double shift[3] = {5.,5.,5.};
int64_t point_num = nx[0] * nx[1] * nx[2];

double rmin[3] = { nucl_coord[0][0], nucl_coord[0][1], nucl_coord[0][2] } ;
double rmax[3] = { nucl_coord[0][0], nucl_coord[0][1], nucl_coord[0][2] } ;

for (size_t i=0 ; i<nucl_num ; ++i) {
  for (int j=0 ; j<3 ; ++j) {
    rmin[j] = nucl_coord[i][j] < rmin[j] ? nucl_coord[i][j] : rmin[j];
    rmax[j] = nucl_coord[i][j] > rmax[j] ? nucl_coord[i][j] : rmax[j];
  }
}

rmin[0] -= shift[0]; rmin[1] -= shift[1]; rmin[2] -= shift[2];
rmax[0] += shift[0]; rmax[1] += shift[1]; rmax[2] += shift[2];

double step[3];

double* linspace[3];
for (int i=0 ; i<3 ; ++i) {

  linspace[i] = (double*) calloc( sizeof(double), nx[i] );

  if (linspace[i] == NULL) {
    fprintf(stderr, "Allocation failed (linspace)\n");
    exit(1);
  }

  step[i] = (rmax[i] - rmin[i]) / ((double) (nx[i]-1));

  for (size_t j=0 ; j<nx[i] ; ++j) {
    linspace[i][j] = rmin[i] + j*step[i];
  }

}

double dr = step[0] * step[1] * step[2];

double* point = (double*) calloc(sizeof(double), 3*point_num);

if (point == NULL) {
  fprintf(stderr, "Allocation failed (point)\n");
  exit(1);
}

size_t m = 0;
for (size_t i=0 ; i<nx[0] ; ++i) {
  for (size_t j=0 ; j<nx[1] ; ++j) {
    for (size_t k=0 ; k<nx[2] ; ++k) {

      point[m] = linspace[0][i];
      m++;

      point[m] = linspace[1][j];
      m++;

      point[m] = linspace[2][k];
      m++;

    }
  }
}

rc = qmckl_set_point(context, 'N', point_num, point, (point_num*3));

if (rc != QMCKL_SUCCESS) {
  fprintf(stderr, "Error setting points:\n%s\n", qmckl_string_of_error(rc));
  exit(1);
}

int64_t mo_num;
rc = qmckl_get_mo_basis_mo_num(context, &mo_num);

long before, after;
struct timeval timecheck;

double* mo_value = (double*) calloc(sizeof(double), point_num*mo_num);
if (mo_value == NULL) {
  fprintf(stderr, "Allocation failed (mo_value)\n");
  exit(1);
}

gettimeofday(&timecheck, NULL);
before = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;

rc = qmckl_get_mo_basis_mo_value(context, mo_value, point_num*mo_num);
if (rc != QMCKL_SUCCESS) {
  fprintf(stderr, "Error getting mo_value)\n");
  exit(1);
}

gettimeofday(&timecheck, NULL);
after = (long)timecheck.tv_sec * 1000 + (long)timecheck.tv_usec / 1000;

printf("Number of MOs: %ld\n", (long) mo_num);
printf("Number of grid points: %ld\n", (long) point_num);
printf("Execution time : %f seconds\n", (after - before)*1.e-3);

double* overlap = (double*) malloc (mo_num*mo_num*sizeof(double));

  rc = qmckl_dgemm(context, 'N', 'T', mo_num, mo_num, point_num, dr,
                   mo_value, mo_num, mo_value, mo_num, 0.0,
                   overlap, mo_num);

  for (size_t i=0 ; i<mo_num ; ++i) {
    printf("%4ld", i);
    for (size_t j=0 ; j<mo_num ; ++j) {
      printf(" %f", overlap[i*mo_num+j]);
    }
    printf("\n");
  }

  // Clean-up and exit
  free(mo_value);
  free(overlap);

  rc = qmckl_context_destroy(context);
  if (rc != QMCKL_SUCCESS) {
    fprintf(stderr, "Error destroying context)\n");
    exit(1);
  }

  return 0;
}
