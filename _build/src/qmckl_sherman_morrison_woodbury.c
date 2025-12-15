/* C sources */
/* Common includes and macros used by all the Sherman-Morrison-Woodbury kernels. */

#include <stdbool.h>
#include <math.h>
#include "qmckl.h"
#include "config.h"
#include "assert.h"
#include "stdio.h"



/* ~qmckl_sm_naive_hpc~ is a high performance variation of */
/* ~qmckl_sm_naive~ written in C. It is used in cases when ~Dim~ is */
/* smaller than the leading dimension ~LDS~, irrespective of whetether ~LDS~ */
/* includes zero padding to benefit from SIMD instructions or not. Cases like this */
/* include situations where one wants to apply updates to a square submatrix of the */
/* full matrix. */
/* It takes advantage of memory aligned data and assumes no data dependencies */
/* inside the loops. The loops are fully vectorised whenever ~Dim~ is an integer */
/* multiple of ~SIMD_LENGTH~. */

qmckl_exit_code qmckl_sm_naive_hpc(
    const qmckl_context context,
    const uint64_t LDS,
    const uint64_t Dim,
    const uint64_t N_updates,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
      return qmckl_failwith( context,
          QMCKL_NULL_CONTEXT,
          "qmckl_sm_naive_hpc",
          NULL);
  }

  double __attribute__((aligned(8))) C[Dim];
  double __attribute__((aligned(8))) D[LDS];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = S^{-1} x u_l
    for (uint64_t i = 0; i < Dim; i++) {
      C[i] = 0.0f;
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < Dim; j++) {
        C[i] += Slater_inv[i * LDS + j] * Updates[l * LDS + j];
      }
    }

    // Denominator: v_l^T * C
    const int cui = Updates_index[l] - 1;
    double den = 1.0f + C[cui];

    if (fabs(den) < breakdown)
      return QMCKL_FAILURE;

    double iden = 1.0f / den;

    // Update det(A)
    if (determinant)
      *determinant *= den;

    // selecting column: v_l^T * S_inv
    IVDEP
    ALIGNED
    for (uint64_t j = 0; j < Dim; j++) {
      D[j] = Slater_inv[cui * LDS + j];
    }

    // A^{-1} = A^{-1} - C x D / den
    for (uint64_t i = 0; i < Dim; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < Dim; j++) {
        const double update = C[i] * D[j] * iden;
        Slater_inv[i * LDS + j] -= update;
      }
    }
    l += 1;
  }
  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_sm_naive_2(
    const qmckl_context context,
    const uint64_t N_updates,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_naive_2",
        NULL);
  }

  #define D2_P ((1+(2-1)/SIMD_LENGTH)*SIMD_LENGTH)

  double __attribute__((aligned(8))) C[2];
  double __attribute__((aligned(8))) D[D2_P];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = A^{-1} x U_l
    for (uint64_t i = 0; i < 2; i++) {
      C[i] = 0;
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D2_P; j++) {
        C[i] += Slater_inv[i * D2_P + j] * Updates[l * D2_P + j];
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0f + C[cui];

    if (fabs(den) < breakdown) {
      return QMCKL_FAILURE;
    }
    double iden = 1.0f / den;

    // Update det(A)
    if (determinant)
      *determinant *= den;

    // selecting column: D = v_l^T * S_inv
    IVDEP
    ALIGNED
    for (uint64_t j = 0; j < D2_P; j++) {
      D[j] = Slater_inv[cui * D2_P + j];
    }

    // A^{-1} = A^{-1} - C x D / den
    for (uint64_t i = 0; i < 2; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D2_P; j++) {
        double update = C[i] * D[j] * iden;
        Slater_inv[i * D2_P + j] -= update;
      }
    }

    l += 1;
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_sm_naive_3(
    const qmckl_context context,
    const uint64_t N_updates,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_naive_3",
        NULL);
  }

  #define D3_P ((1+(3-1)/SIMD_LENGTH)*SIMD_LENGTH)

  double __attribute__((aligned(8))) C[3];
  double __attribute__((aligned(8))) D[D3_P];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = A^{-1} x U_l
    for (uint64_t i = 0; i < 3; i++) {
      C[i] = 0;
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D3_P; j++) {
        C[i] += Slater_inv[i * D3_P + j] * Updates[l * D3_P + j];
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0f + C[cui];

    if (fabs(den) < breakdown) {
      return QMCKL_FAILURE;
    }
    double iden = 1.0f / den;

    // Update det(A)
    if (determinant)
      *determinant *= den;

    // selecting column: D = v_l^T * S_inv
    IVDEP
    ALIGNED
    for (uint64_t j = 0; j < D3_P; j++) {
      D[j] = Slater_inv[cui * D3_P + j];
    }

    // A^{-1} = A^{-1} - C x D / den
    for (uint64_t i = 0; i < 3; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D3_P; j++) {
        double update = C[i] * D[j] * iden;
        Slater_inv[i * D3_P + j] -= update;
      }
    }

    l += 1;
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_sm_naive_4(
    const qmckl_context context,
    const uint64_t N_updates,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_naive_4",
        NULL);
  }

  #define D4_P ((1+(4-1)/SIMD_LENGTH)*SIMD_LENGTH)

  double __attribute__((aligned(8))) C[4];
  double __attribute__((aligned(8))) D[D4_P];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = A^{-1} x U_l
    for (uint64_t i = 0; i < 4; i++) {
      C[i] = 0;
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D4_P; j++) {
        C[i] += Slater_inv[i * D4_P + j] * Updates[l * D4_P + j];
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0f + C[cui];

    if (fabs(den) < breakdown) {
      return QMCKL_FAILURE;
    }
    double iden = 1.0f / den;

    // Update det(A)
    if (determinant)
      *determinant *= den;

    // selecting column: D = v_l^T * S_inv
    IVDEP
    ALIGNED
    for (uint64_t j = 0; j < D4_P; j++) {
      D[j] = Slater_inv[cui * D4_P + j];
    }

    // A^{-1} = A^{-1} - C x D / den
    for (uint64_t i = 0; i < 4; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D4_P; j++) {
        double update = C[i] * D[j] * iden;
        Slater_inv[i * D4_P + j] -= update;
      }
    }

    l += 1;
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_sm_naive_5(
    const qmckl_context context,
    const uint64_t N_updates,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_naive_5",
        NULL);
  }

  #define D5_P ((1+(5-1)/SIMD_LENGTH)*SIMD_LENGTH)

  double __attribute__((aligned(8))) C[5];
  double __attribute__((aligned(8))) D[D5_P];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = A^{-1} x U_l
    for (uint64_t i = 0; i < 5; i++) {
      C[i] = 0;
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D5_P; j++) {
        C[i] += Slater_inv[i * D5_P + j] * Updates[l * D5_P + j];
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0f + C[cui];

    if (fabs(den) < breakdown) {
      return QMCKL_FAILURE;
    }
    double iden = 1.0f / den;

    // Update det(A)
    if (determinant)
      *determinant *= den;

    // selecting column: D = v_l^T * S_inv
    IVDEP
    ALIGNED
    for (uint64_t j = 0; j < D5_P; j++) {
      D[j] = Slater_inv[cui * D5_P + j];
    }

    // A^{-1} = A^{-1} - C x D / den
    for (uint64_t i = 0; i < 5; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D5_P; j++) {
        double update = C[i] * D[j] * iden;
        Slater_inv[i * D5_P + j] -= update;
      }
    }

    l += 1;
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_sm_naive_6(
    const qmckl_context context,
    const uint64_t N_updates,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_naive_6",
        NULL);
  }

  #define D6_P ((1+(6-1)/SIMD_LENGTH)*SIMD_LENGTH)

  double __attribute__((aligned(8))) C[6];
  double __attribute__((aligned(8))) D[D6_P];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = A^{-1} x U_l
    for (uint64_t i = 0; i < 6; i++) {
      C[i] = 0;
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D6_P; j++) {
        C[i] += Slater_inv[i * D6_P + j] * Updates[l * D6_P + j];
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0f + C[cui];

    if (fabs(den) < breakdown) {
      return QMCKL_FAILURE;
    }
    double iden = 1.0f / den;

    // Update det(A)
    if (determinant)
      *determinant *= den;

    // selecting column: D = v_l^T * S_inv
    IVDEP
    ALIGNED
    for (uint64_t j = 0; j < D6_P; j++) {
      D[j] = Slater_inv[cui * D6_P + j];
    }

    // A^{-1} = A^{-1} - C x D / den
    for (uint64_t i = 0; i < 6; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D6_P; j++) {
        double update = C[i] * D[j] * iden;
        Slater_inv[i * D6_P + j] -= update;
      }
    }

    l += 1;
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_sm_naive_7(
    const qmckl_context context,
    const uint64_t N_updates,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_naive_7",
        NULL);
  }

  #define D7_P ((1+(7-1)/SIMD_LENGTH)*SIMD_LENGTH)

  double __attribute__((aligned(8))) C[7];
  double __attribute__((aligned(8))) D[D7_P];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = A^{-1} x U_l
    for (uint64_t i = 0; i < 7; i++) {
      C[i] = 0;
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D7_P; j++) {
        C[i] += Slater_inv[i * D7_P + j] * Updates[l * D7_P + j];
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0f + C[cui];

    if (fabs(den) < breakdown) {
      return QMCKL_FAILURE;
    }
    double iden = 1.0f / den;

    // Update det(A)
    if (determinant)
      *determinant *= den;

    // selecting column: D = v_l^T * S_inv
    IVDEP
    ALIGNED
    for (uint64_t j = 0; j < D7_P; j++) {
      D[j] = Slater_inv[cui * D7_P + j];
    }

    // A^{-1} = A^{-1} - C x D / den
    for (uint64_t i = 0; i < 7; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D7_P; j++) {
        double update = C[i] * D[j] * iden;
        Slater_inv[i * D7_P + j] -= update;
      }
    }

    l += 1;
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_sm_naive_8(
    const qmckl_context context,
    const uint64_t N_updates,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_naive_8",
        NULL);
  }

  #define D8_P ((1+(8-1)/SIMD_LENGTH)*SIMD_LENGTH)

  double __attribute__((aligned(8))) C[8];
  double __attribute__((aligned(8))) D[D8_P];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = A^{-1} x U_l
    for (uint64_t i = 0; i < 8; i++) {
      C[i] = 0;
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D8_P; j++) {
        C[i] += Slater_inv[i * D8_P + j] * Updates[l * D8_P + j];
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0f + C[cui];

    if (fabs(den) < breakdown) {
      return QMCKL_FAILURE;
    }
    double iden = 1.0f / den;

    // Update det(A)
    if (determinant)
      *determinant *= den;

    // selecting column: D = v_l^T * S_inv
    IVDEP
    ALIGNED
    for (uint64_t j = 0; j < D8_P; j++) {
      D[j] = Slater_inv[cui * D8_P + j];
    }

    // A^{-1} = A^{-1} - C x D / den
    for (uint64_t i = 0; i < 8; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D8_P; j++) {
        double update = C[i] * D[j] * iden;
        Slater_inv[i * D8_P + j] -= update;
      }
    }

    l += 1;
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_sm_naive_9(
    const qmckl_context context,
    const uint64_t N_updates,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_naive_9",
        NULL);
  }

  #define D9_P ((1+(9-1)/SIMD_LENGTH)*SIMD_LENGTH)

  double __attribute__((aligned(8))) C[9];
  double __attribute__((aligned(8))) D[D9_P];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = A^{-1} x U_l
    for (uint64_t i = 0; i < 9; i++) {
      C[i] = 0;
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D9_P; j++) {
        C[i] += Slater_inv[i * D9_P + j] * Updates[l * D9_P + j];
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0f + C[cui];

    if (fabs(den) < breakdown) {
      return QMCKL_FAILURE;
    }
    double iden = 1.0f / den;

    // Update det(A)
    if (determinant)
      *determinant *= den;

    // selecting column: D = v_l^T * S_inv
    IVDEP
    ALIGNED
    for (uint64_t j = 0; j < D9_P; j++) {
      D[j] = Slater_inv[cui * D9_P + j];
    }

    // A^{-1} = A^{-1} - C x D / den
    for (uint64_t i = 0; i < 9; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D9_P; j++) {
        double update = C[i] * D[j] * iden;
        Slater_inv[i * D9_P + j] -= update;
      }
    }

    l += 1;
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_sm_naive_10(
    const qmckl_context context,
    const uint64_t N_updates,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_naive_10",
        NULL);
  }

  #define D10_P ((1+(10-1)/SIMD_LENGTH)*SIMD_LENGTH)

  double __attribute__((aligned(8))) C[10];
  double __attribute__((aligned(8))) D[D10_P];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = A^{-1} x U_l
    for (uint64_t i = 0; i < 10; i++) {
      C[i] = 0;
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D10_P; j++) {
        C[i] += Slater_inv[i * D10_P + j] * Updates[l * D10_P + j];
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0f + C[cui];

    if (fabs(den) < breakdown) {
      return QMCKL_FAILURE;
    }
    double iden = 1.0f / den;

    // Update det(A)
    if (determinant)
      *determinant *= den;

    // selecting column: D = v_l^T * S_inv
    IVDEP
    ALIGNED
    for (uint64_t j = 0; j < D10_P; j++) {
      D[j] = Slater_inv[cui * D10_P + j];
    }

    // A^{-1} = A^{-1} - C x D / den
    for (uint64_t i = 0; i < 10; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D10_P; j++) {
        double update = C[i] * D[j] * iden;
        Slater_inv[i * D10_P + j] -= update;
      }
    }

    l += 1;
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_sm_naive_11(
    const qmckl_context context,
    const uint64_t N_updates,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_naive_11",
        NULL);
  }

  #define D11_P ((1+(11-1)/SIMD_LENGTH)*SIMD_LENGTH)

  double __attribute__((aligned(8))) C[11];
  double __attribute__((aligned(8))) D[D11_P];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = A^{-1} x U_l
    for (uint64_t i = 0; i < 11; i++) {
      C[i] = 0;
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D11_P; j++) {
        C[i] += Slater_inv[i * D11_P + j] * Updates[l * D11_P + j];
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0f + C[cui];

    if (fabs(den) < breakdown) {
      return QMCKL_FAILURE;
    }
    double iden = 1.0f / den;

    // Update det(A)
    if (determinant)
      *determinant *= den;

    // selecting column: D = v_l^T * S_inv
    IVDEP
    ALIGNED
    for (uint64_t j = 0; j < D11_P; j++) {
      D[j] = Slater_inv[cui * D11_P + j];
    }

    // A^{-1} = A^{-1} - C x D / den
    for (uint64_t i = 0; i < 11; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D11_P; j++) {
        double update = C[i] * D[j] * iden;
        Slater_inv[i * D11_P + j] -= update;
      }
    }

    l += 1;
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_sm_naive_12(
    const qmckl_context context,
    const uint64_t N_updates,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_naive_12",
        NULL);
  }

  #define D12_P ((1+(12-1)/SIMD_LENGTH)*SIMD_LENGTH)

  double __attribute__((aligned(8))) C[12];
  double __attribute__((aligned(8))) D[D12_P];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = A^{-1} x U_l
    for (uint64_t i = 0; i < 12; i++) {
      C[i] = 0;
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D12_P; j++) {
        C[i] += Slater_inv[i * D12_P + j] * Updates[l * D12_P + j];
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0f + C[cui];

    if (fabs(den) < breakdown) {
      return QMCKL_FAILURE;
    }
    double iden = 1.0f / den;

    // Update det(A)
    if (determinant)
      *determinant *= den;

    // selecting column: D = v_l^T * S_inv
    IVDEP
    ALIGNED
    for (uint64_t j = 0; j < D12_P; j++) {
      D[j] = Slater_inv[cui * D12_P + j];
    }

    // A^{-1} = A^{-1} - C x D / den
    for (uint64_t i = 0; i < 12; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D12_P; j++) {
        double update = C[i] * D[j] * iden;
        Slater_inv[i * D12_P + j] -= update;
      }
    }

    l += 1;
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_sm_naive_13(
    const qmckl_context context,
    const uint64_t N_updates,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_naive_13",
        NULL);
  }

  #define D13_P ((1+(13-1)/SIMD_LENGTH)*SIMD_LENGTH)

  double __attribute__((aligned(8))) C[13];
  double __attribute__((aligned(8))) D[D13_P];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = A^{-1} x U_l
    for (uint64_t i = 0; i < 13; i++) {
      C[i] = 0;
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D13_P; j++) {
        C[i] += Slater_inv[i * D13_P + j] * Updates[l * D13_P + j];
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0f + C[cui];

    if (fabs(den) < breakdown) {
      return QMCKL_FAILURE;
    }
    double iden = 1.0f / den;

    // Update det(A)
    if (determinant)
      *determinant *= den;

    // selecting column: D = v_l^T * S_inv
    IVDEP
    ALIGNED
    for (uint64_t j = 0; j < D13_P; j++) {
      D[j] = Slater_inv[cui * D13_P + j];
    }

    // A^{-1} = A^{-1} - C x D / den
    for (uint64_t i = 0; i < 13; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D13_P; j++) {
        double update = C[i] * D[j] * iden;
        Slater_inv[i * D13_P + j] -= update;
      }
    }

    l += 1;
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_sm_naive_14(
    const qmckl_context context,
    const uint64_t N_updates,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_naive_14",
        NULL);
  }

  #define D14_P ((1+(14-1)/SIMD_LENGTH)*SIMD_LENGTH)

  double __attribute__((aligned(8))) C[14];
  double __attribute__((aligned(8))) D[D14_P];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = A^{-1} x U_l
    for (uint64_t i = 0; i < 14; i++) {
      C[i] = 0;
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D14_P; j++) {
        C[i] += Slater_inv[i * D14_P + j] * Updates[l * D14_P + j];
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0f + C[cui];

    if (fabs(den) < breakdown) {
      return QMCKL_FAILURE;
    }
    double iden = 1.0f / den;

    // Update det(A)
    if (determinant)
      *determinant *= den;

    // selecting column: D = v_l^T * S_inv
    IVDEP
    ALIGNED
    for (uint64_t j = 0; j < D14_P; j++) {
      D[j] = Slater_inv[cui * D14_P + j];
    }

    // A^{-1} = A^{-1} - C x D / den
    for (uint64_t i = 0; i < 14; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D14_P; j++) {
        double update = C[i] * D[j] * iden;
        Slater_inv[i * D14_P + j] -= update;
      }
    }

    l += 1;
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_sm_naive_15(
    const qmckl_context context,
    const uint64_t N_updates,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_naive_15",
        NULL);
  }

  #define D15_P ((1+(15-1)/SIMD_LENGTH)*SIMD_LENGTH)

  double __attribute__((aligned(8))) C[15];
  double __attribute__((aligned(8))) D[D15_P];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = A^{-1} x U_l
    for (uint64_t i = 0; i < 15; i++) {
      C[i] = 0;
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D15_P; j++) {
        C[i] += Slater_inv[i * D15_P + j] * Updates[l * D15_P + j];
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0f + C[cui];

    if (fabs(den) < breakdown) {
      return QMCKL_FAILURE;
    }
    double iden = 1.0f / den;

    // Update det(A)
    if (determinant)
      *determinant *= den;

    // selecting column: D = v_l^T * S_inv
    IVDEP
    ALIGNED
    for (uint64_t j = 0; j < D15_P; j++) {
      D[j] = Slater_inv[cui * D15_P + j];
    }

    // A^{-1} = A^{-1} - C x D / den
    for (uint64_t i = 0; i < 15; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D15_P; j++) {
        double update = C[i] * D[j] * iden;
        Slater_inv[i * D15_P + j] -= update;
      }
    }

    l += 1;
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_sm_naive_16(
    const qmckl_context context,
    const uint64_t N_updates,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_naive_16",
        NULL);
  }

  #define D16_P ((1+(16-1)/SIMD_LENGTH)*SIMD_LENGTH)

  double __attribute__((aligned(8))) C[16];
  double __attribute__((aligned(8))) D[D16_P];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = A^{-1} x U_l
    for (uint64_t i = 0; i < 16; i++) {
      C[i] = 0;
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D16_P; j++) {
        C[i] += Slater_inv[i * D16_P + j] * Updates[l * D16_P + j];
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0f + C[cui];

    if (fabs(den) < breakdown) {
      return QMCKL_FAILURE;
    }
    double iden = 1.0f / den;

    // Update det(A)
    if (determinant)
      *determinant *= den;

    // selecting column: D = v_l^T * S_inv
    IVDEP
    ALIGNED
    for (uint64_t j = 0; j < D16_P; j++) {
      D[j] = Slater_inv[cui * D16_P + j];
    }

    // A^{-1} = A^{-1} - C x D / den
    for (uint64_t i = 0; i < 16; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D16_P; j++) {
        double update = C[i] * D[j] * iden;
        Slater_inv[i * D16_P + j] -= update;
      }
    }

    l += 1;
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_sm_naive_17(
    const qmckl_context context,
    const uint64_t N_updates,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_naive_17",
        NULL);
  }

  #define D17_P ((1+(17-1)/SIMD_LENGTH)*SIMD_LENGTH)

  double __attribute__((aligned(8))) C[17];
  double __attribute__((aligned(8))) D[D17_P];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = A^{-1} x U_l
    for (uint64_t i = 0; i < 17; i++) {
      C[i] = 0;
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D17_P; j++) {
        C[i] += Slater_inv[i * D17_P + j] * Updates[l * D17_P + j];
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0f + C[cui];

    if (fabs(den) < breakdown) {
      return QMCKL_FAILURE;
    }
    double iden = 1.0f / den;

    // Update det(A)
    if (determinant)
      *determinant *= den;

    // selecting column: D = v_l^T * S_inv
    IVDEP
    ALIGNED
    for (uint64_t j = 0; j < D17_P; j++) {
      D[j] = Slater_inv[cui * D17_P + j];
    }

    // A^{-1} = A^{-1} - C x D / den
    for (uint64_t i = 0; i < 17; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D17_P; j++) {
        double update = C[i] * D[j] * iden;
        Slater_inv[i * D17_P + j] -= update;
      }
    }

    l += 1;
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_sm_naive_18(
    const qmckl_context context,
    const uint64_t N_updates,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_naive_18",
        NULL);
  }

  #define D18_P ((1+(18-1)/SIMD_LENGTH)*SIMD_LENGTH)

  double __attribute__((aligned(8))) C[18];
  double __attribute__((aligned(8))) D[D18_P];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = A^{-1} x U_l
    for (uint64_t i = 0; i < 18; i++) {
      C[i] = 0;
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D18_P; j++) {
        C[i] += Slater_inv[i * D18_P + j] * Updates[l * D18_P + j];
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0f + C[cui];

    if (fabs(den) < breakdown) {
      return QMCKL_FAILURE;
    }
    double iden = 1.0f / den;

    // Update det(A)
    if (determinant)
      *determinant *= den;

    // selecting column: D = v_l^T * S_inv
    IVDEP
    ALIGNED
    for (uint64_t j = 0; j < D18_P; j++) {
      D[j] = Slater_inv[cui * D18_P + j];
    }

    // A^{-1} = A^{-1} - C x D / den
    for (uint64_t i = 0; i < 18; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D18_P; j++) {
        double update = C[i] * D[j] * iden;
        Slater_inv[i * D18_P + j] -= update;
      }
    }

    l += 1;
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_sm_naive_19(
    const qmckl_context context,
    const uint64_t N_updates,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_naive_19",
        NULL);
  }

  #define D19_P ((1+(19-1)/SIMD_LENGTH)*SIMD_LENGTH)

  double __attribute__((aligned(8))) C[19];
  double __attribute__((aligned(8))) D[D19_P];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = A^{-1} x U_l
    for (uint64_t i = 0; i < 19; i++) {
      C[i] = 0;
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D19_P; j++) {
        C[i] += Slater_inv[i * D19_P + j] * Updates[l * D19_P + j];
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0f + C[cui];

    if (fabs(den) < breakdown) {
      return QMCKL_FAILURE;
    }
    double iden = 1.0f / den;

    // Update det(A)
    if (determinant)
      *determinant *= den;

    // selecting column: D = v_l^T * S_inv
    IVDEP
    ALIGNED
    for (uint64_t j = 0; j < D19_P; j++) {
      D[j] = Slater_inv[cui * D19_P + j];
    }

    // A^{-1} = A^{-1} - C x D / den
    for (uint64_t i = 0; i < 19; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D19_P; j++) {
        double update = C[i] * D[j] * iden;
        Slater_inv[i * D19_P + j] -= update;
      }
    }

    l += 1;
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_sm_naive_20(
    const qmckl_context context,
    const uint64_t N_updates,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_naive_20",
        NULL);
  }

  #define D20_P ((1+(20-1)/SIMD_LENGTH)*SIMD_LENGTH)

  double __attribute__((aligned(8))) C[20];
  double __attribute__((aligned(8))) D[D20_P];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = A^{-1} x U_l
    for (uint64_t i = 0; i < 20; i++) {
      C[i] = 0;
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D20_P; j++) {
        C[i] += Slater_inv[i * D20_P + j] * Updates[l * D20_P + j];
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0f + C[cui];

    if (fabs(den) < breakdown) {
      return QMCKL_FAILURE;
    }
    double iden = 1.0f / den;

    // Update det(A)
    if (determinant)
      *determinant *= den;

    // selecting column: D = v_l^T * S_inv
    IVDEP
    ALIGNED
    for (uint64_t j = 0; j < D20_P; j++) {
      D[j] = Slater_inv[cui * D20_P + j];
    }

    // A^{-1} = A^{-1} - C x D / den
    for (uint64_t i = 0; i < 20; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D20_P; j++) {
        double update = C[i] * D[j] * iden;
        Slater_inv[i * D20_P + j] -= update;
      }
    }

    l += 1;
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_sm_naive_21(
    const qmckl_context context,
    const uint64_t N_updates,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_naive_21",
        NULL);
  }

  #define D21_P ((1+(21-1)/SIMD_LENGTH)*SIMD_LENGTH)

  double __attribute__((aligned(8))) C[21];
  double __attribute__((aligned(8))) D[D21_P];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = A^{-1} x U_l
    for (uint64_t i = 0; i < 21; i++) {
      C[i] = 0;
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D21_P; j++) {
        C[i] += Slater_inv[i * D21_P + j] * Updates[l * D21_P + j];
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0f + C[cui];

    if (fabs(den) < breakdown) {
      return QMCKL_FAILURE;
    }
    double iden = 1.0f / den;

    // Update det(A)
    if (determinant)
      *determinant *= den;

    // selecting column: D = v_l^T * S_inv
    IVDEP
    ALIGNED
    for (uint64_t j = 0; j < D21_P; j++) {
      D[j] = Slater_inv[cui * D21_P + j];
    }

    // A^{-1} = A^{-1} - C x D / den
    for (uint64_t i = 0; i < 21; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D21_P; j++) {
        double update = C[i] * D[j] * iden;
        Slater_inv[i * D21_P + j] -= update;
      }
    }

    l += 1;
  }

  return QMCKL_SUCCESS;
}



/* ~qmckl_sm_naive~ is a generic function that contains decision making logic that calls the proper kernel based on the used library configuration (~--enable-doc~ and ~--enable-hpc~) and the passed array dimensions ~LDS~ and ~Dim~. */

qmckl_exit_code qmckl_sm_naive(const qmckl_context context,
    const uint64_t LDS,
    const uint64_t Dim,
    const uint64_t N_updates,
    const double* Updates,
    const uint64_t* Updates_index,
    const double breakdown,
    double* Slater_inv,
    double* determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(
        context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_naive",
        NULL);
  }

  #ifdef HAVE_HPC__BROKEN_WITH_CRAY
  if (LDS == (1+(Dim-1)/SIMD_LENGTH)*SIMD_LENGTH) { // Most cases
    switch (Dim) {
      case 2:  
        return qmckl_sm_naive_2(context,
          N_updates,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 3:  
        return qmckl_sm_naive_3(context,
          N_updates,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 4:  
        return qmckl_sm_naive_4(context,
          N_updates,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 5:  
        return qmckl_sm_naive_5(context,
          N_updates,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 6:  
        return qmckl_sm_naive_6(context,
          N_updates,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 7:  
        return qmckl_sm_naive_7(context,
          N_updates,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 8:  
        return qmckl_sm_naive_8(context,
          N_updates,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 9:  
        return qmckl_sm_naive_9(context,
          N_updates,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 10:  
        return qmckl_sm_naive_10(context,
          N_updates,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 11:  
        return qmckl_sm_naive_11(context,
          N_updates,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 12:  
        return qmckl_sm_naive_12(context,
          N_updates,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 13:  
        return qmckl_sm_naive_13(context,
          N_updates,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 14:  
        return qmckl_sm_naive_14(context,
          N_updates,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 15:  
        return qmckl_sm_naive_15(context,
          N_updates,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 16:  
        return qmckl_sm_naive_16(context,
          N_updates,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 17:  
        return qmckl_sm_naive_17(context,
          N_updates,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 18:  
        return qmckl_sm_naive_18(context,
          N_updates,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 19:  
        return qmckl_sm_naive_19(context,
          N_updates,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 20:  
        return qmckl_sm_naive_20(context,
          N_updates,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 21:  
        return qmckl_sm_naive_21(context,
          N_updates,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
    }
  }
  else
  { // Updating smaller sub-matrix
    return qmckl_sm_naive_hpc(
        context,
        LDS,
        Dim,
        N_updates,
        Updates,
        Updates_index,
        breakdown,
        Slater_inv,
        determinant);
  }
  #else
  return qmckl_sm_naive_doc(
      context,
      LDS,
      Dim,
      N_updates,
      Updates,
      Updates_index,
      breakdown,
      Slater_inv,
      determinant);
  #endif
  
  return QMCKL_FAILURE;
}

/* C sources */

qmckl_exit_code qmckl_sm_splitting_core_hpc(
    const qmckl_context context,
    uint64_t LDS,
    uint64_t Dim,
    uint64_t N_updates,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict later_updates,
    uint64_t* restrict later_index,
    uint64_t* restrict later,
    double* restrict determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(
        context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_splitting_core_hpc",
        NULL);
  }

  double __attribute__((aligned(8))) C[LDS];
  double __attribute__((aligned(8))) D[LDS];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = S^{-1} x U_l
    for (uint64_t i = 0; i < Dim; i++) {
      C[i] = 0.0f;
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < LDS; j++) {
        C[i] += Slater_inv[i * LDS + j] * Updates[l * LDS + j];
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0f + C[cui];
    if (fabs(den) < breakdown) {
      // U_l = U_l / 2: split the update in 2 equal halves and save the
      // second halve in later_updates
      IVDEP
      ALIGNED
      for (uint64_t i = 0; i < LDS; i++) {
        later_updates[*later * LDS + i] = Updates[l * LDS + i] * 0.5f;
        C[i] *= 0.5f;
      }
      later_index[*later] = Updates_index[l];
      (*later)++;

      den = 1.0f + C[cui];
    } // From here onwards we continue with applying the first halve of the
      // update to Slater_inv
    double iden = 1.0f / den;

    if (determinant)
      *determinant *= den;

    // D = v^T x S^{-1} : 1 x LDS
    IVDEP
    ALIGNED
    for (uint64_t j = 0; j < LDS; j++) {
      D[j] = Slater_inv[cui * LDS + j];
    }

    // S^{-1} = S^{-1} - C x D / den
    for (uint64_t i = 0; i < Dim; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < LDS; j++) {
        const double update = C[i] * D[j] * iden;
        Slater_inv[i * LDS + j] -= update;
      }
    }
    l += 1;
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_sm_splitting_core_2(
    const qmckl_context context,
    uint64_t N_updates,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict later_updates,
    uint64_t* restrict later_index,
    uint64_t* restrict later,
    double* restrict determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(
        context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_splitting_core_2",
        NULL);
  }

  double __attribute__((aligned(8))) C[D2_P];
  double __attribute__((aligned(8))) D[D2_P];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = S^{-1} x U_l
    for (uint64_t i = 0; i < 2; i++) {
      C[i] = 0.0f;
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D2_P; j++) {
        C[i] += Slater_inv[i * D2_P + j] * Updates[l * D2_P + j];
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0f + C[cui];
    if (fabs(den) < breakdown) {
      // U_l = U_l / 2: split the update in 2 equal halves and save the
      // second halve in later_updates
      IVDEP
      ALIGNED
      for (uint64_t i = 0; i < D2_P; i++) {
        later_updates[*later * D2_P + i] = Updates[l * D2_P + i] * 0.5f;
        C[i] *= 0.5f;
      }
      later_index[*later] = Updates_index[l];
      (*later)++;

      den = 1.0f + C[cui];
    } // From here onwards we continue with applying the first halve of the
      // update to Slater_inv
    double iden = 1.0f / den;

    if (determinant)
      *determinant *= den;

    // D = v^T x S^{-1} : 1 x D2_P
    IVDEP
    ALIGNED
    for (uint64_t j = 0; j < D2_P; j++) {
      D[j] = Slater_inv[cui * D2_P + j];
    }

    // S^{-1} = S^{-1} - C x D / den
    for (uint64_t i = 0; i < 2; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D2_P; j++) {
        const double update = C[i] * D[j] * iden;
        Slater_inv[i * D2_P + j] -= update;
      }
    }
    l += 1;
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_sm_splitting_core_3(
    const qmckl_context context,
    uint64_t N_updates,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict later_updates,
    uint64_t* restrict later_index,
    uint64_t* restrict later,
    double* restrict determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(
        context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_splitting_core_3",
        NULL);
  }

  double __attribute__((aligned(8))) C[D3_P];
  double __attribute__((aligned(8))) D[D3_P];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = S^{-1} x U_l
    for (uint64_t i = 0; i < 3; i++) {
      C[i] = 0.0f;
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D3_P; j++) {
        C[i] += Slater_inv[i * D3_P + j] * Updates[l * D3_P + j];
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0f + C[cui];
    if (fabs(den) < breakdown) {
      // U_l = U_l / 2: split the update in 2 equal halves and save the
      // second halve in later_updates
      IVDEP
      ALIGNED
      for (uint64_t i = 0; i < D3_P; i++) {
        later_updates[*later * D3_P + i] = Updates[l * D3_P + i] * 0.5f;
        C[i] *= 0.5f;
      }
      later_index[*later] = Updates_index[l];
      (*later)++;

      den = 1.0f + C[cui];
    } // From here onwards we continue with applying the first halve of the
      // update to Slater_inv
    double iden = 1.0f / den;

    if (determinant)
      *determinant *= den;

    // D = v^T x S^{-1} : 1 x D3_P
    IVDEP
    ALIGNED
    for (uint64_t j = 0; j < D3_P; j++) {
      D[j] = Slater_inv[cui * D3_P + j];
    }

    // S^{-1} = S^{-1} - C x D / den
    for (uint64_t i = 0; i < 3; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D3_P; j++) {
        const double update = C[i] * D[j] * iden;
        Slater_inv[i * D3_P + j] -= update;
      }
    }
    l += 1;
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_sm_splitting_core_4(
    const qmckl_context context,
    uint64_t N_updates,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict later_updates,
    uint64_t* restrict later_index,
    uint64_t* restrict later,
    double* restrict determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(
        context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_splitting_core_4",
        NULL);
  }

  double __attribute__((aligned(8))) C[D4_P];
  double __attribute__((aligned(8))) D[D4_P];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = S^{-1} x U_l
    for (uint64_t i = 0; i < 4; i++) {
      C[i] = 0.0f;
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D4_P; j++) {
        C[i] += Slater_inv[i * D4_P + j] * Updates[l * D4_P + j];
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0f + C[cui];
    if (fabs(den) < breakdown) {
      // U_l = U_l / 2: split the update in 2 equal halves and save the
      // second halve in later_updates
      IVDEP
      ALIGNED
      for (uint64_t i = 0; i < D4_P; i++) {
        later_updates[*later * D4_P + i] = Updates[l * D4_P + i] * 0.5f;
        C[i] *= 0.5f;
      }
      later_index[*later] = Updates_index[l];
      (*later)++;

      den = 1.0f + C[cui];
    } // From here onwards we continue with applying the first halve of the
      // update to Slater_inv
    double iden = 1.0f / den;

    if (determinant)
      *determinant *= den;

    // D = v^T x S^{-1} : 1 x D4_P
    IVDEP
    ALIGNED
    for (uint64_t j = 0; j < D4_P; j++) {
      D[j] = Slater_inv[cui * D4_P + j];
    }

    // S^{-1} = S^{-1} - C x D / den
    for (uint64_t i = 0; i < 4; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D4_P; j++) {
        const double update = C[i] * D[j] * iden;
        Slater_inv[i * D4_P + j] -= update;
      }
    }
    l += 1;
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_sm_splitting_core_5(
    const qmckl_context context,
    uint64_t N_updates,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict later_updates,
    uint64_t* restrict later_index,
    uint64_t* restrict later,
    double* restrict determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(
        context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_splitting_core_5",
        NULL);
  }

  double __attribute__((aligned(8))) C[D5_P];
  double __attribute__((aligned(8))) D[D5_P];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = S^{-1} x U_l
    for (uint64_t i = 0; i < 5; i++) {
      C[i] = 0.0f;
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D5_P; j++) {
        C[i] += Slater_inv[i * D5_P + j] * Updates[l * D5_P + j];
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0f + C[cui];
    if (fabs(den) < breakdown) {
      // U_l = U_l / 2: split the update in 2 equal halves and save the
      // second halve in later_updates
      IVDEP
      ALIGNED
      for (uint64_t i = 0; i < D5_P; i++) {
        later_updates[*later * D5_P + i] = Updates[l * D5_P + i] * 0.5f;
        C[i] *= 0.5f;
      }
      later_index[*later] = Updates_index[l];
      (*later)++;

      den = 1.0f + C[cui];
    } // From here onwards we continue with applying the first halve of the
      // update to Slater_inv
    double iden = 1.0f / den;

    if (determinant)
      *determinant *= den;

    // D = v^T x S^{-1} : 1 x D5_P
    IVDEP
    ALIGNED
    for (uint64_t j = 0; j < D5_P; j++) {
      D[j] = Slater_inv[cui * D5_P + j];
    }

    // S^{-1} = S^{-1} - C x D / den
    for (uint64_t i = 0; i < 5; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D5_P; j++) {
        const double update = C[i] * D[j] * iden;
        Slater_inv[i * D5_P + j] -= update;
      }
    }
    l += 1;
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_sm_splitting_core_6(
    const qmckl_context context,
    uint64_t N_updates,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict later_updates,
    uint64_t* restrict later_index,
    uint64_t* restrict later,
    double* restrict determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(
        context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_splitting_core_6",
        NULL);
  }

  double __attribute__((aligned(8))) C[D6_P];
  double __attribute__((aligned(8))) D[D6_P];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = S^{-1} x U_l
    for (uint64_t i = 0; i < 6; i++) {
      C[i] = 0.0f;
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D6_P; j++) {
        C[i] += Slater_inv[i * D6_P + j] * Updates[l * D6_P + j];
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0f + C[cui];
    if (fabs(den) < breakdown) {
      // U_l = U_l / 2: split the update in 2 equal halves and save the
      // second halve in later_updates
      IVDEP
      ALIGNED
      for (uint64_t i = 0; i < D6_P; i++) {
        later_updates[*later * D6_P + i] = Updates[l * D6_P + i] * 0.5f;
        C[i] *= 0.5f;
      }
      later_index[*later] = Updates_index[l];
      (*later)++;

      den = 1.0f + C[cui];
    } // From here onwards we continue with applying the first halve of the
      // update to Slater_inv
    double iden = 1.0f / den;

    if (determinant)
      *determinant *= den;

    // D = v^T x S^{-1} : 1 x D6_P
    IVDEP
    ALIGNED
    for (uint64_t j = 0; j < D6_P; j++) {
      D[j] = Slater_inv[cui * D6_P + j];
    }

    // S^{-1} = S^{-1} - C x D / den
    for (uint64_t i = 0; i < 6; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D6_P; j++) {
        const double update = C[i] * D[j] * iden;
        Slater_inv[i * D6_P + j] -= update;
      }
    }
    l += 1;
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_sm_splitting_core_7(
    const qmckl_context context,
    uint64_t N_updates,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict later_updates,
    uint64_t* restrict later_index,
    uint64_t* restrict later,
    double* restrict determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(
        context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_splitting_core_7",
        NULL);
  }

  double __attribute__((aligned(8))) C[D7_P];
  double __attribute__((aligned(8))) D[D7_P];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = S^{-1} x U_l
    for (uint64_t i = 0; i < 7; i++) {
      C[i] = 0.0f;
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D7_P; j++) {
        C[i] += Slater_inv[i * D7_P + j] * Updates[l * D7_P + j];
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0f + C[cui];
    if (fabs(den) < breakdown) {
      // U_l = U_l / 2: split the update in 2 equal halves and save the
      // second halve in later_updates
      IVDEP
      ALIGNED
      for (uint64_t i = 0; i < D7_P; i++) {
        later_updates[*later * D7_P + i] = Updates[l * D7_P + i] * 0.5f;
        C[i] *= 0.5f;
      }
      later_index[*later] = Updates_index[l];
      (*later)++;

      den = 1.0f + C[cui];
    } // From here onwards we continue with applying the first halve of the
      // update to Slater_inv
    double iden = 1.0f / den;

    if (determinant)
      *determinant *= den;

    // D = v^T x S^{-1} : 1 x D7_P
    IVDEP
    ALIGNED
    for (uint64_t j = 0; j < D7_P; j++) {
      D[j] = Slater_inv[cui * D7_P + j];
    }

    // S^{-1} = S^{-1} - C x D / den
    for (uint64_t i = 0; i < 7; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D7_P; j++) {
        const double update = C[i] * D[j] * iden;
        Slater_inv[i * D7_P + j] -= update;
      }
    }
    l += 1;
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_sm_splitting_core_8(
    const qmckl_context context,
    uint64_t N_updates,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict later_updates,
    uint64_t* restrict later_index,
    uint64_t* restrict later,
    double* restrict determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(
        context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_splitting_core_8",
        NULL);
  }

  double __attribute__((aligned(8))) C[D8_P];
  double __attribute__((aligned(8))) D[D8_P];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = S^{-1} x U_l
    for (uint64_t i = 0; i < 8; i++) {
      C[i] = 0.0f;
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D8_P; j++) {
        C[i] += Slater_inv[i * D8_P + j] * Updates[l * D8_P + j];
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0f + C[cui];
    if (fabs(den) < breakdown) {
      // U_l = U_l / 2: split the update in 2 equal halves and save the
      // second halve in later_updates
      IVDEP
      ALIGNED
      for (uint64_t i = 0; i < D8_P; i++) {
        later_updates[*later * D8_P + i] = Updates[l * D8_P + i] * 0.5f;
        C[i] *= 0.5f;
      }
      later_index[*later] = Updates_index[l];
      (*later)++;

      den = 1.0f + C[cui];
    } // From here onwards we continue with applying the first halve of the
      // update to Slater_inv
    double iden = 1.0f / den;

    if (determinant)
      *determinant *= den;

    // D = v^T x S^{-1} : 1 x D8_P
    IVDEP
    ALIGNED
    for (uint64_t j = 0; j < D8_P; j++) {
      D[j] = Slater_inv[cui * D8_P + j];
    }

    // S^{-1} = S^{-1} - C x D / den
    for (uint64_t i = 0; i < 8; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D8_P; j++) {
        const double update = C[i] * D[j] * iden;
        Slater_inv[i * D8_P + j] -= update;
      }
    }
    l += 1;
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_sm_splitting_core_9(
    const qmckl_context context,
    uint64_t N_updates,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict later_updates,
    uint64_t* restrict later_index,
    uint64_t* restrict later,
    double* restrict determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(
        context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_splitting_core_9",
        NULL);
  }

  double __attribute__((aligned(8))) C[D9_P];
  double __attribute__((aligned(8))) D[D9_P];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = S^{-1} x U_l
    for (uint64_t i = 0; i < 9; i++) {
      C[i] = 0.0f;
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D9_P; j++) {
        C[i] += Slater_inv[i * D9_P + j] * Updates[l * D9_P + j];
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0f + C[cui];
    if (fabs(den) < breakdown) {
      // U_l = U_l / 2: split the update in 2 equal halves and save the
      // second halve in later_updates
      IVDEP
      ALIGNED
      for (uint64_t i = 0; i < D9_P; i++) {
        later_updates[*later * D9_P + i] = Updates[l * D9_P + i] * 0.5f;
        C[i] *= 0.5f;
      }
      later_index[*later] = Updates_index[l];
      (*later)++;

      den = 1.0f + C[cui];
    } // From here onwards we continue with applying the first halve of the
      // update to Slater_inv
    double iden = 1.0f / den;

    if (determinant)
      *determinant *= den;

    // D = v^T x S^{-1} : 1 x D9_P
    IVDEP
    ALIGNED
    for (uint64_t j = 0; j < D9_P; j++) {
      D[j] = Slater_inv[cui * D9_P + j];
    }

    // S^{-1} = S^{-1} - C x D / den
    for (uint64_t i = 0; i < 9; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D9_P; j++) {
        const double update = C[i] * D[j] * iden;
        Slater_inv[i * D9_P + j] -= update;
      }
    }
    l += 1;
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_sm_splitting_core_10(
    const qmckl_context context,
    uint64_t N_updates,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict later_updates,
    uint64_t* restrict later_index,
    uint64_t* restrict later,
    double* restrict determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(
        context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_splitting_core_10",
        NULL);
  }

  double __attribute__((aligned(8))) C[D10_P];
  double __attribute__((aligned(8))) D[D10_P];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = S^{-1} x U_l
    for (uint64_t i = 0; i < 10; i++) {
      C[i] = 0.0f;
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D10_P; j++) {
        C[i] += Slater_inv[i * D10_P + j] * Updates[l * D10_P + j];
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0f + C[cui];
    if (fabs(den) < breakdown) {
      // U_l = U_l / 2: split the update in 2 equal halves and save the
      // second halve in later_updates
      IVDEP
      ALIGNED
      for (uint64_t i = 0; i < D10_P; i++) {
        later_updates[*later * D10_P + i] = Updates[l * D10_P + i] * 0.5f;
        C[i] *= 0.5f;
      }
      later_index[*later] = Updates_index[l];
      (*later)++;

      den = 1.0f + C[cui];
    } // From here onwards we continue with applying the first halve of the
      // update to Slater_inv
    double iden = 1.0f / den;

    if (determinant)
      *determinant *= den;

    // D = v^T x S^{-1} : 1 x D10_P
    IVDEP
    ALIGNED
    for (uint64_t j = 0; j < D10_P; j++) {
      D[j] = Slater_inv[cui * D10_P + j];
    }

    // S^{-1} = S^{-1} - C x D / den
    for (uint64_t i = 0; i < 10; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D10_P; j++) {
        const double update = C[i] * D[j] * iden;
        Slater_inv[i * D10_P + j] -= update;
      }
    }
    l += 1;
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_sm_splitting_core_11(
    const qmckl_context context,
    uint64_t N_updates,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict later_updates,
    uint64_t* restrict later_index,
    uint64_t* restrict later,
    double* restrict determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(
        context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_splitting_core_11",
        NULL);
  }

  double __attribute__((aligned(8))) C[D11_P];
  double __attribute__((aligned(8))) D[D11_P];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = S^{-1} x U_l
    for (uint64_t i = 0; i < 11; i++) {
      C[i] = 0.0f;
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D11_P; j++) {
        C[i] += Slater_inv[i * D11_P + j] * Updates[l * D11_P + j];
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0f + C[cui];
    if (fabs(den) < breakdown) {
      // U_l = U_l / 2: split the update in 2 equal halves and save the
      // second halve in later_updates
      IVDEP
      ALIGNED
      for (uint64_t i = 0; i < D11_P; i++) {
        later_updates[*later * D11_P + i] = Updates[l * D11_P + i] * 0.5f;
        C[i] *= 0.5f;
      }
      later_index[*later] = Updates_index[l];
      (*later)++;

      den = 1.0f + C[cui];
    } // From here onwards we continue with applying the first halve of the
      // update to Slater_inv
    double iden = 1.0f / den;

    if (determinant)
      *determinant *= den;

    // D = v^T x S^{-1} : 1 x D11_P
    IVDEP
    ALIGNED
    for (uint64_t j = 0; j < D11_P; j++) {
      D[j] = Slater_inv[cui * D11_P + j];
    }

    // S^{-1} = S^{-1} - C x D / den
    for (uint64_t i = 0; i < 11; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D11_P; j++) {
        const double update = C[i] * D[j] * iden;
        Slater_inv[i * D11_P + j] -= update;
      }
    }
    l += 1;
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_sm_splitting_core_12(
    const qmckl_context context,
    uint64_t N_updates,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict later_updates,
    uint64_t* restrict later_index,
    uint64_t* restrict later,
    double* restrict determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(
        context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_splitting_core_12",
        NULL);
  }

  double __attribute__((aligned(8))) C[D12_P];
  double __attribute__((aligned(8))) D[D12_P];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = S^{-1} x U_l
    for (uint64_t i = 0; i < 12; i++) {
      C[i] = 0.0f;
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D12_P; j++) {
        C[i] += Slater_inv[i * D12_P + j] * Updates[l * D12_P + j];
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0f + C[cui];
    if (fabs(den) < breakdown) {
      // U_l = U_l / 2: split the update in 2 equal halves and save the
      // second halve in later_updates
      IVDEP
      ALIGNED
      for (uint64_t i = 0; i < D12_P; i++) {
        later_updates[*later * D12_P + i] = Updates[l * D12_P + i] * 0.5f;
        C[i] *= 0.5f;
      }
      later_index[*later] = Updates_index[l];
      (*later)++;

      den = 1.0f + C[cui];
    } // From here onwards we continue with applying the first halve of the
      // update to Slater_inv
    double iden = 1.0f / den;

    if (determinant)
      *determinant *= den;

    // D = v^T x S^{-1} : 1 x D12_P
    IVDEP
    ALIGNED
    for (uint64_t j = 0; j < D12_P; j++) {
      D[j] = Slater_inv[cui * D12_P + j];
    }

    // S^{-1} = S^{-1} - C x D / den
    for (uint64_t i = 0; i < 12; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D12_P; j++) {
        const double update = C[i] * D[j] * iden;
        Slater_inv[i * D12_P + j] -= update;
      }
    }
    l += 1;
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_sm_splitting_core_13(
    const qmckl_context context,
    uint64_t N_updates,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict later_updates,
    uint64_t* restrict later_index,
    uint64_t* restrict later,
    double* restrict determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(
        context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_splitting_core_13",
        NULL);
  }

  double __attribute__((aligned(8))) C[D13_P];
  double __attribute__((aligned(8))) D[D13_P];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = S^{-1} x U_l
    for (uint64_t i = 0; i < 13; i++) {
      C[i] = 0.0f;
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D13_P; j++) {
        C[i] += Slater_inv[i * D13_P + j] * Updates[l * D13_P + j];
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0f + C[cui];
    if (fabs(den) < breakdown) {
      // U_l = U_l / 2: split the update in 2 equal halves and save the
      // second halve in later_updates
      IVDEP
      ALIGNED
      for (uint64_t i = 0; i < D13_P; i++) {
        later_updates[*later * D13_P + i] = Updates[l * D13_P + i] * 0.5f;
        C[i] *= 0.5f;
      }
      later_index[*later] = Updates_index[l];
      (*later)++;

      den = 1.0f + C[cui];
    } // From here onwards we continue with applying the first halve of the
      // update to Slater_inv
    double iden = 1.0f / den;

    if (determinant)
      *determinant *= den;

    // D = v^T x S^{-1} : 1 x D13_P
    IVDEP
    ALIGNED
    for (uint64_t j = 0; j < D13_P; j++) {
      D[j] = Slater_inv[cui * D13_P + j];
    }

    // S^{-1} = S^{-1} - C x D / den
    for (uint64_t i = 0; i < 13; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D13_P; j++) {
        const double update = C[i] * D[j] * iden;
        Slater_inv[i * D13_P + j] -= update;
      }
    }
    l += 1;
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_sm_splitting_core_14(
    const qmckl_context context,
    uint64_t N_updates,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict later_updates,
    uint64_t* restrict later_index,
    uint64_t* restrict later,
    double* restrict determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(
        context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_splitting_core_14",
        NULL);
  }

  double __attribute__((aligned(8))) C[D14_P];
  double __attribute__((aligned(8))) D[D14_P];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = S^{-1} x U_l
    for (uint64_t i = 0; i < 14; i++) {
      C[i] = 0.0f;
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D14_P; j++) {
        C[i] += Slater_inv[i * D14_P + j] * Updates[l * D14_P + j];
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0f + C[cui];
    if (fabs(den) < breakdown) {
      // U_l = U_l / 2: split the update in 2 equal halves and save the
      // second halve in later_updates
      IVDEP
      ALIGNED
      for (uint64_t i = 0; i < D14_P; i++) {
        later_updates[*later * D14_P + i] = Updates[l * D14_P + i] * 0.5f;
        C[i] *= 0.5f;
      }
      later_index[*later] = Updates_index[l];
      (*later)++;

      den = 1.0f + C[cui];
    } // From here onwards we continue with applying the first halve of the
      // update to Slater_inv
    double iden = 1.0f / den;

    if (determinant)
      *determinant *= den;

    // D = v^T x S^{-1} : 1 x D14_P
    IVDEP
    ALIGNED
    for (uint64_t j = 0; j < D14_P; j++) {
      D[j] = Slater_inv[cui * D14_P + j];
    }

    // S^{-1} = S^{-1} - C x D / den
    for (uint64_t i = 0; i < 14; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D14_P; j++) {
        const double update = C[i] * D[j] * iden;
        Slater_inv[i * D14_P + j] -= update;
      }
    }
    l += 1;
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_sm_splitting_core_15(
    const qmckl_context context,
    uint64_t N_updates,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict later_updates,
    uint64_t* restrict later_index,
    uint64_t* restrict later,
    double* restrict determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(
        context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_splitting_core_15",
        NULL);
  }

  double __attribute__((aligned(8))) C[D15_P];
  double __attribute__((aligned(8))) D[D15_P];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = S^{-1} x U_l
    for (uint64_t i = 0; i < 15; i++) {
      C[i] = 0.0f;
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D15_P; j++) {
        C[i] += Slater_inv[i * D15_P + j] * Updates[l * D15_P + j];
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0f + C[cui];
    if (fabs(den) < breakdown) {
      // U_l = U_l / 2: split the update in 2 equal halves and save the
      // second halve in later_updates
      IVDEP
      ALIGNED
      for (uint64_t i = 0; i < D15_P; i++) {
        later_updates[*later * D15_P + i] = Updates[l * D15_P + i] * 0.5f;
        C[i] *= 0.5f;
      }
      later_index[*later] = Updates_index[l];
      (*later)++;

      den = 1.0f + C[cui];
    } // From here onwards we continue with applying the first halve of the
      // update to Slater_inv
    double iden = 1.0f / den;

    if (determinant)
      *determinant *= den;

    // D = v^T x S^{-1} : 1 x D15_P
    IVDEP
    ALIGNED
    for (uint64_t j = 0; j < D15_P; j++) {
      D[j] = Slater_inv[cui * D15_P + j];
    }

    // S^{-1} = S^{-1} - C x D / den
    for (uint64_t i = 0; i < 15; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D15_P; j++) {
        const double update = C[i] * D[j] * iden;
        Slater_inv[i * D15_P + j] -= update;
      }
    }
    l += 1;
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_sm_splitting_core_16(
    const qmckl_context context,
    uint64_t N_updates,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict later_updates,
    uint64_t* restrict later_index,
    uint64_t* restrict later,
    double* restrict determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(
        context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_splitting_core_16",
        NULL);
  }

  double __attribute__((aligned(8))) C[D16_P];
  double __attribute__((aligned(8))) D[D16_P];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = S^{-1} x U_l
    for (uint64_t i = 0; i < 16; i++) {
      C[i] = 0.0f;
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D16_P; j++) {
        C[i] += Slater_inv[i * D16_P + j] * Updates[l * D16_P + j];
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0f + C[cui];
    if (fabs(den) < breakdown) {
      // U_l = U_l / 2: split the update in 2 equal halves and save the
      // second halve in later_updates
      IVDEP
      ALIGNED
      for (uint64_t i = 0; i < D16_P; i++) {
        later_updates[*later * D16_P + i] = Updates[l * D16_P + i] * 0.5f;
        C[i] *= 0.5f;
      }
      later_index[*later] = Updates_index[l];
      (*later)++;

      den = 1.0f + C[cui];
    } // From here onwards we continue with applying the first halve of the
      // update to Slater_inv
    double iden = 1.0f / den;

    if (determinant)
      *determinant *= den;

    // D = v^T x S^{-1} : 1 x D16_P
    IVDEP
    ALIGNED
    for (uint64_t j = 0; j < D16_P; j++) {
      D[j] = Slater_inv[cui * D16_P + j];
    }

    // S^{-1} = S^{-1} - C x D / den
    for (uint64_t i = 0; i < 16; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D16_P; j++) {
        const double update = C[i] * D[j] * iden;
        Slater_inv[i * D16_P + j] -= update;
      }
    }
    l += 1;
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_sm_splitting_core_17(
    const qmckl_context context,
    uint64_t N_updates,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict later_updates,
    uint64_t* restrict later_index,
    uint64_t* restrict later,
    double* restrict determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(
        context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_splitting_core_17",
        NULL);
  }

  double __attribute__((aligned(8))) C[D17_P];
  double __attribute__((aligned(8))) D[D17_P];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = S^{-1} x U_l
    for (uint64_t i = 0; i < 17; i++) {
      C[i] = 0.0f;
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D17_P; j++) {
        C[i] += Slater_inv[i * D17_P + j] * Updates[l * D17_P + j];
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0f + C[cui];
    if (fabs(den) < breakdown) {
      // U_l = U_l / 2: split the update in 2 equal halves and save the
      // second halve in later_updates
      IVDEP
      ALIGNED
      for (uint64_t i = 0; i < D17_P; i++) {
        later_updates[*later * D17_P + i] = Updates[l * D17_P + i] * 0.5f;
        C[i] *= 0.5f;
      }
      later_index[*later] = Updates_index[l];
      (*later)++;

      den = 1.0f + C[cui];
    } // From here onwards we continue with applying the first halve of the
      // update to Slater_inv
    double iden = 1.0f / den;

    if (determinant)
      *determinant *= den;

    // D = v^T x S^{-1} : 1 x D17_P
    IVDEP
    ALIGNED
    for (uint64_t j = 0; j < D17_P; j++) {
      D[j] = Slater_inv[cui * D17_P + j];
    }

    // S^{-1} = S^{-1} - C x D / den
    for (uint64_t i = 0; i < 17; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D17_P; j++) {
        const double update = C[i] * D[j] * iden;
        Slater_inv[i * D17_P + j] -= update;
      }
    }
    l += 1;
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_sm_splitting_core_18(
    const qmckl_context context,
    uint64_t N_updates,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict later_updates,
    uint64_t* restrict later_index,
    uint64_t* restrict later,
    double* restrict determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(
        context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_splitting_core_18",
        NULL);
  }

  double __attribute__((aligned(8))) C[D18_P];
  double __attribute__((aligned(8))) D[D18_P];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = S^{-1} x U_l
    for (uint64_t i = 0; i < 18; i++) {
      C[i] = 0.0f;
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D18_P; j++) {
        C[i] += Slater_inv[i * D18_P + j] * Updates[l * D18_P + j];
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0f + C[cui];
    if (fabs(den) < breakdown) {
      // U_l = U_l / 2: split the update in 2 equal halves and save the
      // second halve in later_updates
      IVDEP
      ALIGNED
      for (uint64_t i = 0; i < D18_P; i++) {
        later_updates[*later * D18_P + i] = Updates[l * D18_P + i] * 0.5f;
        C[i] *= 0.5f;
      }
      later_index[*later] = Updates_index[l];
      (*later)++;

      den = 1.0f + C[cui];
    } // From here onwards we continue with applying the first halve of the
      // update to Slater_inv
    double iden = 1.0f / den;

    if (determinant)
      *determinant *= den;

    // D = v^T x S^{-1} : 1 x D18_P
    IVDEP
    ALIGNED
    for (uint64_t j = 0; j < D18_P; j++) {
      D[j] = Slater_inv[cui * D18_P + j];
    }

    // S^{-1} = S^{-1} - C x D / den
    for (uint64_t i = 0; i < 18; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D18_P; j++) {
        const double update = C[i] * D[j] * iden;
        Slater_inv[i * D18_P + j] -= update;
      }
    }
    l += 1;
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_sm_splitting_core_19(
    const qmckl_context context,
    uint64_t N_updates,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict later_updates,
    uint64_t* restrict later_index,
    uint64_t* restrict later,
    double* restrict determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(
        context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_splitting_core_19",
        NULL);
  }

  double __attribute__((aligned(8))) C[D19_P];
  double __attribute__((aligned(8))) D[D19_P];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = S^{-1} x U_l
    for (uint64_t i = 0; i < 19; i++) {
      C[i] = 0.0f;
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D19_P; j++) {
        C[i] += Slater_inv[i * D19_P + j] * Updates[l * D19_P + j];
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0f + C[cui];
    if (fabs(den) < breakdown) {
      // U_l = U_l / 2: split the update in 2 equal halves and save the
      // second halve in later_updates
      IVDEP
      ALIGNED
      for (uint64_t i = 0; i < D19_P; i++) {
        later_updates[*later * D19_P + i] = Updates[l * D19_P + i] * 0.5f;
        C[i] *= 0.5f;
      }
      later_index[*later] = Updates_index[l];
      (*later)++;

      den = 1.0f + C[cui];
    } // From here onwards we continue with applying the first halve of the
      // update to Slater_inv
    double iden = 1.0f / den;

    if (determinant)
      *determinant *= den;

    // D = v^T x S^{-1} : 1 x D19_P
    IVDEP
    ALIGNED
    for (uint64_t j = 0; j < D19_P; j++) {
      D[j] = Slater_inv[cui * D19_P + j];
    }

    // S^{-1} = S^{-1} - C x D / den
    for (uint64_t i = 0; i < 19; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D19_P; j++) {
        const double update = C[i] * D[j] * iden;
        Slater_inv[i * D19_P + j] -= update;
      }
    }
    l += 1;
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_sm_splitting_core_20(
    const qmckl_context context,
    uint64_t N_updates,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict later_updates,
    uint64_t* restrict later_index,
    uint64_t* restrict later,
    double* restrict determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(
        context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_splitting_core_20",
        NULL);
  }

  double __attribute__((aligned(8))) C[D20_P];
  double __attribute__((aligned(8))) D[D20_P];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = S^{-1} x U_l
    for (uint64_t i = 0; i < 20; i++) {
      C[i] = 0.0f;
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D20_P; j++) {
        C[i] += Slater_inv[i * D20_P + j] * Updates[l * D20_P + j];
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0f + C[cui];
    if (fabs(den) < breakdown) {
      // U_l = U_l / 2: split the update in 2 equal halves and save the
      // second halve in later_updates
      IVDEP
      ALIGNED
      for (uint64_t i = 0; i < D20_P; i++) {
        later_updates[*later * D20_P + i] = Updates[l * D20_P + i] * 0.5f;
        C[i] *= 0.5f;
      }
      later_index[*later] = Updates_index[l];
      (*later)++;

      den = 1.0f + C[cui];
    } // From here onwards we continue with applying the first halve of the
      // update to Slater_inv
    double iden = 1.0f / den;

    if (determinant)
      *determinant *= den;

    // D = v^T x S^{-1} : 1 x D20_P
    IVDEP
    ALIGNED
    for (uint64_t j = 0; j < D20_P; j++) {
      D[j] = Slater_inv[cui * D20_P + j];
    }

    // S^{-1} = S^{-1} - C x D / den
    for (uint64_t i = 0; i < 20; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D20_P; j++) {
        const double update = C[i] * D[j] * iden;
        Slater_inv[i * D20_P + j] -= update;
      }
    }
    l += 1;
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_sm_splitting_core_21(
    const qmckl_context context,
    uint64_t N_updates,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict later_updates,
    uint64_t* restrict later_index,
    uint64_t* restrict later,
    double* restrict determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(
        context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_splitting_core_21",
        NULL);
  }

  double __attribute__((aligned(8))) C[D21_P];
  double __attribute__((aligned(8))) D[D21_P];

  uint64_t l = 0;
  // For each update
  while (l < N_updates) {
    // C = S^{-1} x U_l
    for (uint64_t i = 0; i < 21; i++) {
      C[i] = 0.0f;
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D21_P; j++) {
        C[i] += Slater_inv[i * D21_P + j] * Updates[l * D21_P + j];
      }
    }

    // Denominator
    const int cui = Updates_index[l] - 1;
    double den = 1.0f + C[cui];
    if (fabs(den) < breakdown) {
      // U_l = U_l / 2: split the update in 2 equal halves and save the
      // second halve in later_updates
      IVDEP
      ALIGNED
      for (uint64_t i = 0; i < D21_P; i++) {
        later_updates[*later * D21_P + i] = Updates[l * D21_P + i] * 0.5f;
        C[i] *= 0.5f;
      }
      later_index[*later] = Updates_index[l];
      (*later)++;

      den = 1.0f + C[cui];
    } // From here onwards we continue with applying the first halve of the
      // update to Slater_inv
    double iden = 1.0f / den;

    if (determinant)
      *determinant *= den;

    // D = v^T x S^{-1} : 1 x D21_P
    IVDEP
    ALIGNED
    for (uint64_t j = 0; j < D21_P; j++) {
      D[j] = Slater_inv[cui * D21_P + j];
    }

    // S^{-1} = S^{-1} - C x D / den
    for (uint64_t i = 0; i < 21; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D21_P; j++) {
        const double update = C[i] * D[j] * iden;
        Slater_inv[i * D21_P + j] -= update;
      }
    }
    l += 1;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_sm_splitting_core(
    const qmckl_context context,
    const uint64_t LDS,
    const uint64_t Dim,
    const uint64_t N_updates,
    const double* Updates,
    const uint64_t* Updates_index,
    const double breakdown,
    double* Slater_inv,
    double* later_updates,
    uint64_t* later_index,
    uint64_t* later,
    double* determinant) {

  #ifdef HAVE_HPC__BROKEN_WITH_CRAY
  if (LDS == (1+(Dim-1)/SIMD_LENGTH)*SIMD_LENGTH) { // Most cases
    switch (Dim) {
      case 2: {
        return qmckl_sm_splitting_core_2(
            context,
            N_updates,
            Updates,
            Updates_index,
            breakdown,
            Slater_inv,
            later_updates,
            later_index,
            later,
            determinant);
      }
      case 3: {
        return qmckl_sm_splitting_core_3(
            context,
            N_updates,
            Updates,
            Updates_index,
            breakdown,
            Slater_inv,
            later_updates,
            later_index,
            later,
            determinant);
      }
      case 4: {
        return qmckl_sm_splitting_core_4(
            context,
            N_updates,
            Updates,
            Updates_index,
            breakdown,
            Slater_inv,
            later_updates,
            later_index,
            later,
            determinant);
      }
      case 5: {
        return qmckl_sm_splitting_core_5(
            context,
            N_updates,
            Updates,
            Updates_index,
            breakdown,
            Slater_inv,
            later_updates,
            later_index,
            later,
            determinant);
      }
      case 6: {
        return qmckl_sm_splitting_core_6(
            context,
            N_updates,
            Updates,
            Updates_index,
            breakdown,
            Slater_inv,
            later_updates,
            later_index,
            later,
            determinant);
      }
      case 7: {
        return qmckl_sm_splitting_core_7(
            context,
            N_updates,
            Updates,
            Updates_index,
            breakdown,
            Slater_inv,
            later_updates,
            later_index,
            later,
            determinant);
      }
      case 8: {
        return qmckl_sm_splitting_core_8(
            context,
            N_updates,
            Updates,
            Updates_index,
            breakdown,
            Slater_inv,
            later_updates,
            later_index,
            later,
            determinant);
      }
      case 9: {
        return qmckl_sm_splitting_core_9(
            context,
            N_updates,
            Updates,
            Updates_index,
            breakdown,
            Slater_inv,
            later_updates,
            later_index,
            later,
            determinant);
      }
      case 10: {
        return qmckl_sm_splitting_core_10(
            context,
            N_updates,
            Updates,
            Updates_index,
            breakdown,
            Slater_inv,
            later_updates,
            later_index,
            later,
            determinant);
      }
      case 11: {
        return qmckl_sm_splitting_core_11(
            context,
            N_updates,
            Updates,
            Updates_index,
            breakdown,
            Slater_inv,
            later_updates,
            later_index,
            later,
            determinant);
      }
      case 12: {
        return qmckl_sm_splitting_core_12(
            context,
            N_updates,
            Updates,
            Updates_index,
            breakdown,
            Slater_inv,
            later_updates,
            later_index,
            later,
            determinant);
      }
      case 13: {
        return qmckl_sm_splitting_core_13(
            context,
            N_updates,
            Updates,
            Updates_index,
            breakdown,
            Slater_inv,
            later_updates,
            later_index,
            later,
            determinant);
      }
      case 14: {
        return qmckl_sm_splitting_core_14(
            context,
            N_updates,
            Updates,
            Updates_index,
            breakdown,
            Slater_inv,
            later_updates,
            later_index,
            later,
            determinant);
      }
      case 15: {
        return qmckl_sm_splitting_core_15(
            context,
            N_updates,
            Updates,
            Updates_index,
            breakdown,
            Slater_inv,
            later_updates,
            later_index,
            later,
            determinant);
      }
      case 16: {
        return qmckl_sm_splitting_core_16(
            context,
            N_updates,
            Updates,
            Updates_index,
            breakdown,
            Slater_inv,
            later_updates,
            later_index,
            later,
            determinant);
      }
      case 17: {
        return qmckl_sm_splitting_core_17(
            context,
            N_updates,
            Updates,
            Updates_index,
            breakdown,
            Slater_inv,
            later_updates,
            later_index,
            later,
            determinant);
      }
      case 18: {
        return qmckl_sm_splitting_core_18(
            context,
            N_updates,
            Updates,
            Updates_index,
            breakdown,
            Slater_inv,
            later_updates,
            later_index,
            later,
            determinant);
      }
      case 19: {
        return qmckl_sm_splitting_core_19(
            context,
            N_updates,
            Updates,
            Updates_index,
            breakdown,
            Slater_inv,
            later_updates,
            later_index,
            later,
            determinant);
      }
      case 20: {
        return qmckl_sm_splitting_core_20(
            context,
            N_updates,
            Updates,
            Updates_index,
            breakdown,
            Slater_inv,
            later_updates,
            later_index,
            later,
            determinant);
      }
      case 21: {
        return qmckl_sm_splitting_core_21(
            context,
            N_updates,
            Updates,
            Updates_index,
            breakdown,
            Slater_inv,
            later_updates,
            later_index,
            later,
            determinant);
      }
      default: {
        assert(0 == 1 && "TEMPLATE NOT IMPLEMENTED!");
        break;
      }
    }
  }
  else { // Updating smaller sub-matrix
    return qmckl_sm_splitting_core_hpc(
        context,
        LDS,
        Dim,
        N_updates,
        Updates,
        Updates_index,
        breakdown,
        Slater_inv,
        later_updates,
        later_index,
        later,
        determinant);
  }
  #else
  return qmckl_sm_splitting_core_doc(
      context,
      LDS,
      Dim,
      N_updates,
      Updates,
      Updates_index,
      breakdown,
      Slater_inv,
      later_updates,
      later_index,
      later,
      determinant);
  #endif

  return QMCKL_FAILURE;
}

/* C sources */

qmckl_exit_code qmckl_woodbury_2x2_hpc(const qmckl_context context,
    const uint64_t LDS,
    const uint64_t Dim,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {
/*
    C := S^{-1} * U,    dim x 2
    B := 1 + V * C,     2 x 2
    D := V * S^{-1},    2 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
                          QMCKL_NULL_CONTEXT,
                          "qmckl_woodbury_2x2_hpc",
                          NULL);
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);

  // Compute C = (S^T)^{-1}U : Dim x 2
  double __attribute__((aligned(8))) C[2 * Dim];
  for (uint64_t i = 0; i < Dim; i++) {
      C[i * 2] = 0;
      C[i * 2 + 1] = 0;
      for (uint64_t k = 0; k < LDS; k++) {
          C[i * 2] += Slater_inv[i * LDS + k] * Updates[k];
          C[i * 2 + 1] += Slater_inv[i * LDS + k] * Updates[LDS + k];
      }
  }

  // Compute B = 1 + VC : 2 x 2
  const double B0 = C[row1 * 2] + 1;
  const double B1 = C[row1 * 2 + 1];
  const double B2 = C[row2 * 2];
  const double B3 = C[row2 * 2 + 1] + 1;

  // Check if determinant of inverted matrix is not zero
  double det = B0 * B3 - B1 * B2;
  if (fabs(det) < breakdown) {
      return QMCKL_FAILURE;
  }

  // Update det(S) when passed
  if (determinant)
      *determinant *= det;

  // Compute B^{-1} with explicit formula for 2 x 2 inversion
  double __attribute__((aligned(8))) Binv[4], idet = 1.0 / det;
  Binv[0] = idet * B3;
  Binv[1] = -1.0 * idet * B1;
  Binv[2] = -1.0 * idet * B2;
  Binv[3] = idet * B0;

  // tmp = B^{-1}D : 2 x LDS
  double __attribute__((aligned(8))) tmp[2 * LDS];
  double* r1dim = &(Slater_inv[row1 * LDS]);
  double* r2dim = &(Slater_inv[row2 * LDS]);
  for (uint64_t j = 0; j < LDS; j++) {
      tmp[j] = Binv[0] * r1dim[j] + Binv[1] * r2dim[j];
      tmp[LDS + j] = Binv[2] * r1dim[j] + Binv[3] * r2dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : Dim x LDS
  for (uint64_t i = 0; i < Dim; i++) {
      for (uint64_t j = 0; j < LDS; j++) {
          Slater_inv[i * LDS + j] -= C[i * 2] * tmp[j];
          Slater_inv[i * LDS + j] -= C[i * 2 + 1] * tmp[LDS + j];
      }
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_woodbury_2x2_2(
    const qmckl_context context,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {
/*
    C := S^{-1} * U,    dim x 2
    B := 1 + V * C,     2 x 2
    D := V * S^{-1},    2 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
                          QMCKL_NULL_CONTEXT,
                          "qmckl_woodbury_2x2_2",
                          NULL);
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);

  // Compute C = (S^T)^{-1}U : 2 x 2
  double __attribute__((aligned(8))) C[2 * 2];
  for (uint64_t i = 0; i < 2; i++) {
      C[i * 2] = 0;
      C[i * 2 + 1] = 0;
      IVDEP
      ALIGNED
      for (uint64_t k = 0; k < D2_P; k++) {
          C[i * 2] += Slater_inv[i * D2_P + k] * Updates[k];
          C[i * 2 + 1] += Slater_inv[i * D2_P + k] * Updates[D2_P + k];
      }
  }

  // Compute B = 1 + VC : 2 x 2
  const double B0 = C[row1 * 2] + 1;
  const double B1 = C[row1 * 2 + 1];
  const double B2 = C[row2 * 2];
  const double B3 = C[row2 * 2 + 1] + 1;

  // Check if determinant of inverted matrix is not zero
  double det = B0 * B3 - B1 * B2;
  if (fabs(det) < breakdown) {
      return QMCKL_FAILURE;
  }

  // Update det(S) when passed
  if (determinant)
      *determinant *= det;

  // Compute B^{-1} with explicit formula for 2 x 2 inversion
  double __attribute__((aligned(8))) Binv[4], idet = 1.0 / det;
  Binv[0] = idet * B3;
  Binv[1] = -1.0 * idet * B1;
  Binv[2] = -1.0 * idet * B2;
  Binv[3] = idet * B0;

  // tmp = B^{-1}D : 2 x D2_P
  double __attribute__((aligned(8))) tmp[2 * D2_P];
  double* r1dim = &(Slater_inv[row1 * D2_P]);
  double* r2dim = &(Slater_inv[row2 * D2_P]);
  IVDEP
  ALIGNED
  for (uint64_t j = 0; j < D2_P; j++) {
      tmp[j] = Binv[0] * r1dim[j] + Binv[1] * r2dim[j];
      tmp[D2_P + j] = Binv[2] * r1dim[j] + Binv[3] * r2dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : 2 x D2_P
  for (uint64_t i = 0; i < 2; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D2_P; j++) {
          Slater_inv[i * D2_P + j] -= C[i * 2] * tmp[j];
          Slater_inv[i * D2_P + j] -= C[i * 2 + 1] * tmp[D2_P + j];
      }
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_woodbury_2x2_3(
    const qmckl_context context,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {
/*
    C := S^{-1} * U,    dim x 2
    B := 1 + V * C,     2 x 2
    D := V * S^{-1},    2 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
                          QMCKL_NULL_CONTEXT,
                          "qmckl_woodbury_2x2_3",
                          NULL);
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);

  // Compute C = (S^T)^{-1}U : 3 x 2
  double __attribute__((aligned(8))) C[2 * 3];
  for (uint64_t i = 0; i < 3; i++) {
      C[i * 2] = 0;
      C[i * 2 + 1] = 0;
      IVDEP
      ALIGNED
      for (uint64_t k = 0; k < D3_P; k++) {
          C[i * 2] += Slater_inv[i * D3_P + k] * Updates[k];
          C[i * 2 + 1] += Slater_inv[i * D3_P + k] * Updates[D3_P + k];
      }
  }

  // Compute B = 1 + VC : 2 x 2
  const double B0 = C[row1 * 2] + 1;
  const double B1 = C[row1 * 2 + 1];
  const double B2 = C[row2 * 2];
  const double B3 = C[row2 * 2 + 1] + 1;

  // Check if determinant of inverted matrix is not zero
  double det = B0 * B3 - B1 * B2;
  if (fabs(det) < breakdown) {
      return QMCKL_FAILURE;
  }

  // Update det(S) when passed
  if (determinant)
      *determinant *= det;

  // Compute B^{-1} with explicit formula for 2 x 2 inversion
  double __attribute__((aligned(8))) Binv[4], idet = 1.0 / det;
  Binv[0] = idet * B3;
  Binv[1] = -1.0 * idet * B1;
  Binv[2] = -1.0 * idet * B2;
  Binv[3] = idet * B0;

  // tmp = B^{-1}D : 2 x D3_P
  double __attribute__((aligned(8))) tmp[2 * D3_P];
  double* r1dim = &(Slater_inv[row1 * D3_P]);
  double* r2dim = &(Slater_inv[row2 * D3_P]);
  IVDEP
  ALIGNED
  for (uint64_t j = 0; j < D3_P; j++) {
      tmp[j] = Binv[0] * r1dim[j] + Binv[1] * r2dim[j];
      tmp[D3_P + j] = Binv[2] * r1dim[j] + Binv[3] * r2dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : 3 x D3_P
  for (uint64_t i = 0; i < 3; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D3_P; j++) {
          Slater_inv[i * D3_P + j] -= C[i * 2] * tmp[j];
          Slater_inv[i * D3_P + j] -= C[i * 2 + 1] * tmp[D3_P + j];
      }
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_woodbury_2x2_4(
    const qmckl_context context,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {
/*
    C := S^{-1} * U,    dim x 2
    B := 1 + V * C,     2 x 2
    D := V * S^{-1},    2 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
                          QMCKL_NULL_CONTEXT,
                          "qmckl_woodbury_2x2_4",
                          NULL);
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);

  // Compute C = (S^T)^{-1}U : 4 x 2
  double __attribute__((aligned(8))) C[2 * 4];
  for (uint64_t i = 0; i < 4; i++) {
      C[i * 2] = 0;
      C[i * 2 + 1] = 0;
      IVDEP
      ALIGNED
      for (uint64_t k = 0; k < D4_P; k++) {
          C[i * 2] += Slater_inv[i * D4_P + k] * Updates[k];
          C[i * 2 + 1] += Slater_inv[i * D4_P + k] * Updates[D4_P + k];
      }
  }

  // Compute B = 1 + VC : 2 x 2
  const double B0 = C[row1 * 2] + 1;
  const double B1 = C[row1 * 2 + 1];
  const double B2 = C[row2 * 2];
  const double B3 = C[row2 * 2 + 1] + 1;

  // Check if determinant of inverted matrix is not zero
  double det = B0 * B3 - B1 * B2;
  if (fabs(det) < breakdown) {
      return QMCKL_FAILURE;
  }

  // Update det(S) when passed
  if (determinant)
      *determinant *= det;

  // Compute B^{-1} with explicit formula for 2 x 2 inversion
  double __attribute__((aligned(8))) Binv[4], idet = 1.0 / det;
  Binv[0] = idet * B3;
  Binv[1] = -1.0 * idet * B1;
  Binv[2] = -1.0 * idet * B2;
  Binv[3] = idet * B0;

  // tmp = B^{-1}D : 2 x D4_P
  double __attribute__((aligned(8))) tmp[2 * D4_P];
  double* r1dim = &(Slater_inv[row1 * D4_P]);
  double* r2dim = &(Slater_inv[row2 * D4_P]);
  IVDEP
  ALIGNED
  for (uint64_t j = 0; j < D4_P; j++) {
      tmp[j] = Binv[0] * r1dim[j] + Binv[1] * r2dim[j];
      tmp[D4_P + j] = Binv[2] * r1dim[j] + Binv[3] * r2dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : 4 x D4_P
  for (uint64_t i = 0; i < 4; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D4_P; j++) {
          Slater_inv[i * D4_P + j] -= C[i * 2] * tmp[j];
          Slater_inv[i * D4_P + j] -= C[i * 2 + 1] * tmp[D4_P + j];
      }
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_woodbury_2x2_5(
    const qmckl_context context,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {
/*
    C := S^{-1} * U,    dim x 2
    B := 1 + V * C,     2 x 2
    D := V * S^{-1},    2 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
                          QMCKL_NULL_CONTEXT,
                          "qmckl_woodbury_2x2_5",
                          NULL);
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);

  // Compute C = (S^T)^{-1}U : 5 x 2
  double __attribute__((aligned(8))) C[2 * 5];
  for (uint64_t i = 0; i < 5; i++) {
      C[i * 2] = 0;
      C[i * 2 + 1] = 0;
      IVDEP
      ALIGNED
      for (uint64_t k = 0; k < D5_P; k++) {
          C[i * 2] += Slater_inv[i * D5_P + k] * Updates[k];
          C[i * 2 + 1] += Slater_inv[i * D5_P + k] * Updates[D5_P + k];
      }
  }

  // Compute B = 1 + VC : 2 x 2
  const double B0 = C[row1 * 2] + 1;
  const double B1 = C[row1 * 2 + 1];
  const double B2 = C[row2 * 2];
  const double B3 = C[row2 * 2 + 1] + 1;

  // Check if determinant of inverted matrix is not zero
  double det = B0 * B3 - B1 * B2;
  if (fabs(det) < breakdown) {
      return QMCKL_FAILURE;
  }

  // Update det(S) when passed
  if (determinant)
      *determinant *= det;

  // Compute B^{-1} with explicit formula for 2 x 2 inversion
  double __attribute__((aligned(8))) Binv[4], idet = 1.0 / det;
  Binv[0] = idet * B3;
  Binv[1] = -1.0 * idet * B1;
  Binv[2] = -1.0 * idet * B2;
  Binv[3] = idet * B0;

  // tmp = B^{-1}D : 2 x D5_P
  double __attribute__((aligned(8))) tmp[2 * D5_P];
  double* r1dim = &(Slater_inv[row1 * D5_P]);
  double* r2dim = &(Slater_inv[row2 * D5_P]);
  IVDEP
  ALIGNED
  for (uint64_t j = 0; j < D5_P; j++) {
      tmp[j] = Binv[0] * r1dim[j] + Binv[1] * r2dim[j];
      tmp[D5_P + j] = Binv[2] * r1dim[j] + Binv[3] * r2dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : 5 x D5_P
  for (uint64_t i = 0; i < 5; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D5_P; j++) {
          Slater_inv[i * D5_P + j] -= C[i * 2] * tmp[j];
          Slater_inv[i * D5_P + j] -= C[i * 2 + 1] * tmp[D5_P + j];
      }
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_woodbury_2x2_6(
    const qmckl_context context,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {
/*
    C := S^{-1} * U,    dim x 2
    B := 1 + V * C,     2 x 2
    D := V * S^{-1},    2 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
                          QMCKL_NULL_CONTEXT,
                          "qmckl_woodbury_2x2_6",
                          NULL);
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);

  // Compute C = (S^T)^{-1}U : 6 x 2
  double __attribute__((aligned(8))) C[2 * 6];
  for (uint64_t i = 0; i < 6; i++) {
      C[i * 2] = 0;
      C[i * 2 + 1] = 0;
      IVDEP
      ALIGNED
      for (uint64_t k = 0; k < D6_P; k++) {
          C[i * 2] += Slater_inv[i * D6_P + k] * Updates[k];
          C[i * 2 + 1] += Slater_inv[i * D6_P + k] * Updates[D6_P + k];
      }
  }

  // Compute B = 1 + VC : 2 x 2
  const double B0 = C[row1 * 2] + 1;
  const double B1 = C[row1 * 2 + 1];
  const double B2 = C[row2 * 2];
  const double B3 = C[row2 * 2 + 1] + 1;

  // Check if determinant of inverted matrix is not zero
  double det = B0 * B3 - B1 * B2;
  if (fabs(det) < breakdown) {
      return QMCKL_FAILURE;
  }

  // Update det(S) when passed
  if (determinant)
      *determinant *= det;

  // Compute B^{-1} with explicit formula for 2 x 2 inversion
  double __attribute__((aligned(8))) Binv[4], idet = 1.0 / det;
  Binv[0] = idet * B3;
  Binv[1] = -1.0 * idet * B1;
  Binv[2] = -1.0 * idet * B2;
  Binv[3] = idet * B0;

  // tmp = B^{-1}D : 2 x D6_P
  double __attribute__((aligned(8))) tmp[2 * D6_P];
  double* r1dim = &(Slater_inv[row1 * D6_P]);
  double* r2dim = &(Slater_inv[row2 * D6_P]);
  IVDEP
  ALIGNED
  for (uint64_t j = 0; j < D6_P; j++) {
      tmp[j] = Binv[0] * r1dim[j] + Binv[1] * r2dim[j];
      tmp[D6_P + j] = Binv[2] * r1dim[j] + Binv[3] * r2dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : 6 x D6_P
  for (uint64_t i = 0; i < 6; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D6_P; j++) {
          Slater_inv[i * D6_P + j] -= C[i * 2] * tmp[j];
          Slater_inv[i * D6_P + j] -= C[i * 2 + 1] * tmp[D6_P + j];
      }
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_woodbury_2x2_7(
    const qmckl_context context,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {
/*
    C := S^{-1} * U,    dim x 2
    B := 1 + V * C,     2 x 2
    D := V * S^{-1},    2 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
                          QMCKL_NULL_CONTEXT,
                          "qmckl_woodbury_2x2_7",
                          NULL);
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);

  // Compute C = (S^T)^{-1}U : 7 x 2
  double __attribute__((aligned(8))) C[2 * 7];
  for (uint64_t i = 0; i < 7; i++) {
      C[i * 2] = 0;
      C[i * 2 + 1] = 0;
      IVDEP
      ALIGNED
      for (uint64_t k = 0; k < D7_P; k++) {
          C[i * 2] += Slater_inv[i * D7_P + k] * Updates[k];
          C[i * 2 + 1] += Slater_inv[i * D7_P + k] * Updates[D7_P + k];
      }
  }

  // Compute B = 1 + VC : 2 x 2
  const double B0 = C[row1 * 2] + 1;
  const double B1 = C[row1 * 2 + 1];
  const double B2 = C[row2 * 2];
  const double B3 = C[row2 * 2 + 1] + 1;

  // Check if determinant of inverted matrix is not zero
  double det = B0 * B3 - B1 * B2;
  if (fabs(det) < breakdown) {
      return QMCKL_FAILURE;
  }

  // Update det(S) when passed
  if (determinant)
      *determinant *= det;

  // Compute B^{-1} with explicit formula for 2 x 2 inversion
  double __attribute__((aligned(8))) Binv[4], idet = 1.0 / det;
  Binv[0] = idet * B3;
  Binv[1] = -1.0 * idet * B1;
  Binv[2] = -1.0 * idet * B2;
  Binv[3] = idet * B0;

  // tmp = B^{-1}D : 2 x D7_P
  double __attribute__((aligned(8))) tmp[2 * D7_P];
  double* r1dim = &(Slater_inv[row1 * D7_P]);
  double* r2dim = &(Slater_inv[row2 * D7_P]);
  IVDEP
  ALIGNED
  for (uint64_t j = 0; j < D7_P; j++) {
      tmp[j] = Binv[0] * r1dim[j] + Binv[1] * r2dim[j];
      tmp[D7_P + j] = Binv[2] * r1dim[j] + Binv[3] * r2dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : 7 x D7_P
  for (uint64_t i = 0; i < 7; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D7_P; j++) {
          Slater_inv[i * D7_P + j] -= C[i * 2] * tmp[j];
          Slater_inv[i * D7_P + j] -= C[i * 2 + 1] * tmp[D7_P + j];
      }
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_woodbury_2x2_8(
    const qmckl_context context,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {
/*
    C := S^{-1} * U,    dim x 2
    B := 1 + V * C,     2 x 2
    D := V * S^{-1},    2 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
                          QMCKL_NULL_CONTEXT,
                          "qmckl_woodbury_2x2_8",
                          NULL);
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);

  // Compute C = (S^T)^{-1}U : 8 x 2
  double __attribute__((aligned(8))) C[2 * 8];
  for (uint64_t i = 0; i < 8; i++) {
      C[i * 2] = 0;
      C[i * 2 + 1] = 0;
      IVDEP
      ALIGNED
      for (uint64_t k = 0; k < D8_P; k++) {
          C[i * 2] += Slater_inv[i * D8_P + k] * Updates[k];
          C[i * 2 + 1] += Slater_inv[i * D8_P + k] * Updates[D8_P + k];
      }
  }

  // Compute B = 1 + VC : 2 x 2
  const double B0 = C[row1 * 2] + 1;
  const double B1 = C[row1 * 2 + 1];
  const double B2 = C[row2 * 2];
  const double B3 = C[row2 * 2 + 1] + 1;

  // Check if determinant of inverted matrix is not zero
  double det = B0 * B3 - B1 * B2;
  if (fabs(det) < breakdown) {
      return QMCKL_FAILURE;
  }

  // Update det(S) when passed
  if (determinant)
      *determinant *= det;

  // Compute B^{-1} with explicit formula for 2 x 2 inversion
  double __attribute__((aligned(8))) Binv[4], idet = 1.0 / det;
  Binv[0] = idet * B3;
  Binv[1] = -1.0 * idet * B1;
  Binv[2] = -1.0 * idet * B2;
  Binv[3] = idet * B0;

  // tmp = B^{-1}D : 2 x D8_P
  double __attribute__((aligned(8))) tmp[2 * D8_P];
  double* r1dim = &(Slater_inv[row1 * D8_P]);
  double* r2dim = &(Slater_inv[row2 * D8_P]);
  IVDEP
  ALIGNED
  for (uint64_t j = 0; j < D8_P; j++) {
      tmp[j] = Binv[0] * r1dim[j] + Binv[1] * r2dim[j];
      tmp[D8_P + j] = Binv[2] * r1dim[j] + Binv[3] * r2dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : 8 x D8_P
  for (uint64_t i = 0; i < 8; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D8_P; j++) {
          Slater_inv[i * D8_P + j] -= C[i * 2] * tmp[j];
          Slater_inv[i * D8_P + j] -= C[i * 2 + 1] * tmp[D8_P + j];
      }
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_woodbury_2x2_9(
    const qmckl_context context,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {
/*
    C := S^{-1} * U,    dim x 2
    B := 1 + V * C,     2 x 2
    D := V * S^{-1},    2 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
                          QMCKL_NULL_CONTEXT,
                          "qmckl_woodbury_2x2_9",
                          NULL);
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);

  // Compute C = (S^T)^{-1}U : 9 x 2
  double __attribute__((aligned(8))) C[2 * 9];
  for (uint64_t i = 0; i < 9; i++) {
      C[i * 2] = 0;
      C[i * 2 + 1] = 0;
      IVDEP
      ALIGNED
      for (uint64_t k = 0; k < D9_P; k++) {
          C[i * 2] += Slater_inv[i * D9_P + k] * Updates[k];
          C[i * 2 + 1] += Slater_inv[i * D9_P + k] * Updates[D9_P + k];
      }
  }

  // Compute B = 1 + VC : 2 x 2
  const double B0 = C[row1 * 2] + 1;
  const double B1 = C[row1 * 2 + 1];
  const double B2 = C[row2 * 2];
  const double B3 = C[row2 * 2 + 1] + 1;

  // Check if determinant of inverted matrix is not zero
  double det = B0 * B3 - B1 * B2;
  if (fabs(det) < breakdown) {
      return QMCKL_FAILURE;
  }

  // Update det(S) when passed
  if (determinant)
      *determinant *= det;

  // Compute B^{-1} with explicit formula for 2 x 2 inversion
  double __attribute__((aligned(8))) Binv[4], idet = 1.0 / det;
  Binv[0] = idet * B3;
  Binv[1] = -1.0 * idet * B1;
  Binv[2] = -1.0 * idet * B2;
  Binv[3] = idet * B0;

  // tmp = B^{-1}D : 2 x D9_P
  double __attribute__((aligned(8))) tmp[2 * D9_P];
  double* r1dim = &(Slater_inv[row1 * D9_P]);
  double* r2dim = &(Slater_inv[row2 * D9_P]);
  IVDEP
  ALIGNED
  for (uint64_t j = 0; j < D9_P; j++) {
      tmp[j] = Binv[0] * r1dim[j] + Binv[1] * r2dim[j];
      tmp[D9_P + j] = Binv[2] * r1dim[j] + Binv[3] * r2dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : 9 x D9_P
  for (uint64_t i = 0; i < 9; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D9_P; j++) {
          Slater_inv[i * D9_P + j] -= C[i * 2] * tmp[j];
          Slater_inv[i * D9_P + j] -= C[i * 2 + 1] * tmp[D9_P + j];
      }
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_woodbury_2x2_10(
    const qmckl_context context,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {
/*
    C := S^{-1} * U,    dim x 2
    B := 1 + V * C,     2 x 2
    D := V * S^{-1},    2 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
                          QMCKL_NULL_CONTEXT,
                          "qmckl_woodbury_2x2_10",
                          NULL);
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);

  // Compute C = (S^T)^{-1}U : 10 x 2
  double __attribute__((aligned(8))) C[2 * 10];
  for (uint64_t i = 0; i < 10; i++) {
      C[i * 2] = 0;
      C[i * 2 + 1] = 0;
      IVDEP
      ALIGNED
      for (uint64_t k = 0; k < D10_P; k++) {
          C[i * 2] += Slater_inv[i * D10_P + k] * Updates[k];
          C[i * 2 + 1] += Slater_inv[i * D10_P + k] * Updates[D10_P + k];
      }
  }

  // Compute B = 1 + VC : 2 x 2
  const double B0 = C[row1 * 2] + 1;
  const double B1 = C[row1 * 2 + 1];
  const double B2 = C[row2 * 2];
  const double B3 = C[row2 * 2 + 1] + 1;

  // Check if determinant of inverted matrix is not zero
  double det = B0 * B3 - B1 * B2;
  if (fabs(det) < breakdown) {
      return QMCKL_FAILURE;
  }

  // Update det(S) when passed
  if (determinant)
      *determinant *= det;

  // Compute B^{-1} with explicit formula for 2 x 2 inversion
  double __attribute__((aligned(8))) Binv[4], idet = 1.0 / det;
  Binv[0] = idet * B3;
  Binv[1] = -1.0 * idet * B1;
  Binv[2] = -1.0 * idet * B2;
  Binv[3] = idet * B0;

  // tmp = B^{-1}D : 2 x D10_P
  double __attribute__((aligned(8))) tmp[2 * D10_P];
  double* r1dim = &(Slater_inv[row1 * D10_P]);
  double* r2dim = &(Slater_inv[row2 * D10_P]);
  IVDEP
  ALIGNED
  for (uint64_t j = 0; j < D10_P; j++) {
      tmp[j] = Binv[0] * r1dim[j] + Binv[1] * r2dim[j];
      tmp[D10_P + j] = Binv[2] * r1dim[j] + Binv[3] * r2dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : 10 x D10_P
  for (uint64_t i = 0; i < 10; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D10_P; j++) {
          Slater_inv[i * D10_P + j] -= C[i * 2] * tmp[j];
          Slater_inv[i * D10_P + j] -= C[i * 2 + 1] * tmp[D10_P + j];
      }
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_woodbury_2x2_11(
    const qmckl_context context,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {
/*
    C := S^{-1} * U,    dim x 2
    B := 1 + V * C,     2 x 2
    D := V * S^{-1},    2 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
                          QMCKL_NULL_CONTEXT,
                          "qmckl_woodbury_2x2_11",
                          NULL);
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);

  // Compute C = (S^T)^{-1}U : 11 x 2
  double __attribute__((aligned(8))) C[2 * 11];
  for (uint64_t i = 0; i < 11; i++) {
      C[i * 2] = 0;
      C[i * 2 + 1] = 0;
      IVDEP
      ALIGNED
      for (uint64_t k = 0; k < D11_P; k++) {
          C[i * 2] += Slater_inv[i * D11_P + k] * Updates[k];
          C[i * 2 + 1] += Slater_inv[i * D11_P + k] * Updates[D11_P + k];
      }
  }

  // Compute B = 1 + VC : 2 x 2
  const double B0 = C[row1 * 2] + 1;
  const double B1 = C[row1 * 2 + 1];
  const double B2 = C[row2 * 2];
  const double B3 = C[row2 * 2 + 1] + 1;

  // Check if determinant of inverted matrix is not zero
  double det = B0 * B3 - B1 * B2;
  if (fabs(det) < breakdown) {
      return QMCKL_FAILURE;
  }

  // Update det(S) when passed
  if (determinant)
      *determinant *= det;

  // Compute B^{-1} with explicit formula for 2 x 2 inversion
  double __attribute__((aligned(8))) Binv[4], idet = 1.0 / det;
  Binv[0] = idet * B3;
  Binv[1] = -1.0 * idet * B1;
  Binv[2] = -1.0 * idet * B2;
  Binv[3] = idet * B0;

  // tmp = B^{-1}D : 2 x D11_P
  double __attribute__((aligned(8))) tmp[2 * D11_P];
  double* r1dim = &(Slater_inv[row1 * D11_P]);
  double* r2dim = &(Slater_inv[row2 * D11_P]);
  IVDEP
  ALIGNED
  for (uint64_t j = 0; j < D11_P; j++) {
      tmp[j] = Binv[0] * r1dim[j] + Binv[1] * r2dim[j];
      tmp[D11_P + j] = Binv[2] * r1dim[j] + Binv[3] * r2dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : 11 x D11_P
  for (uint64_t i = 0; i < 11; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D11_P; j++) {
          Slater_inv[i * D11_P + j] -= C[i * 2] * tmp[j];
          Slater_inv[i * D11_P + j] -= C[i * 2 + 1] * tmp[D11_P + j];
      }
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_woodbury_2x2_12(
    const qmckl_context context,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {
/*
    C := S^{-1} * U,    dim x 2
    B := 1 + V * C,     2 x 2
    D := V * S^{-1},    2 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
                          QMCKL_NULL_CONTEXT,
                          "qmckl_woodbury_2x2_12",
                          NULL);
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);

  // Compute C = (S^T)^{-1}U : 12 x 2
  double __attribute__((aligned(8))) C[2 * 12];
  for (uint64_t i = 0; i < 12; i++) {
      C[i * 2] = 0;
      C[i * 2 + 1] = 0;
      IVDEP
      ALIGNED
      for (uint64_t k = 0; k < D12_P; k++) {
          C[i * 2] += Slater_inv[i * D12_P + k] * Updates[k];
          C[i * 2 + 1] += Slater_inv[i * D12_P + k] * Updates[D12_P + k];
      }
  }

  // Compute B = 1 + VC : 2 x 2
  const double B0 = C[row1 * 2] + 1;
  const double B1 = C[row1 * 2 + 1];
  const double B2 = C[row2 * 2];
  const double B3 = C[row2 * 2 + 1] + 1;

  // Check if determinant of inverted matrix is not zero
  double det = B0 * B3 - B1 * B2;
  if (fabs(det) < breakdown) {
      return QMCKL_FAILURE;
  }

  // Update det(S) when passed
  if (determinant)
      *determinant *= det;

  // Compute B^{-1} with explicit formula for 2 x 2 inversion
  double __attribute__((aligned(8))) Binv[4], idet = 1.0 / det;
  Binv[0] = idet * B3;
  Binv[1] = -1.0 * idet * B1;
  Binv[2] = -1.0 * idet * B2;
  Binv[3] = idet * B0;

  // tmp = B^{-1}D : 2 x D12_P
  double __attribute__((aligned(8))) tmp[2 * D12_P];
  double* r1dim = &(Slater_inv[row1 * D12_P]);
  double* r2dim = &(Slater_inv[row2 * D12_P]);
  IVDEP
  ALIGNED
  for (uint64_t j = 0; j < D12_P; j++) {
      tmp[j] = Binv[0] * r1dim[j] + Binv[1] * r2dim[j];
      tmp[D12_P + j] = Binv[2] * r1dim[j] + Binv[3] * r2dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : 12 x D12_P
  for (uint64_t i = 0; i < 12; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D12_P; j++) {
          Slater_inv[i * D12_P + j] -= C[i * 2] * tmp[j];
          Slater_inv[i * D12_P + j] -= C[i * 2 + 1] * tmp[D12_P + j];
      }
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_woodbury_2x2_13(
    const qmckl_context context,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {
/*
    C := S^{-1} * U,    dim x 2
    B := 1 + V * C,     2 x 2
    D := V * S^{-1},    2 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
                          QMCKL_NULL_CONTEXT,
                          "qmckl_woodbury_2x2_13",
                          NULL);
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);

  // Compute C = (S^T)^{-1}U : 13 x 2
  double __attribute__((aligned(8))) C[2 * 13];
  for (uint64_t i = 0; i < 13; i++) {
      C[i * 2] = 0;
      C[i * 2 + 1] = 0;
      IVDEP
      ALIGNED
      for (uint64_t k = 0; k < D13_P; k++) {
          C[i * 2] += Slater_inv[i * D13_P + k] * Updates[k];
          C[i * 2 + 1] += Slater_inv[i * D13_P + k] * Updates[D13_P + k];
      }
  }

  // Compute B = 1 + VC : 2 x 2
  const double B0 = C[row1 * 2] + 1;
  const double B1 = C[row1 * 2 + 1];
  const double B2 = C[row2 * 2];
  const double B3 = C[row2 * 2 + 1] + 1;

  // Check if determinant of inverted matrix is not zero
  double det = B0 * B3 - B1 * B2;
  if (fabs(det) < breakdown) {
      return QMCKL_FAILURE;
  }

  // Update det(S) when passed
  if (determinant)
      *determinant *= det;

  // Compute B^{-1} with explicit formula for 2 x 2 inversion
  double __attribute__((aligned(8))) Binv[4], idet = 1.0 / det;
  Binv[0] = idet * B3;
  Binv[1] = -1.0 * idet * B1;
  Binv[2] = -1.0 * idet * B2;
  Binv[3] = idet * B0;

  // tmp = B^{-1}D : 2 x D13_P
  double __attribute__((aligned(8))) tmp[2 * D13_P];
  double* r1dim = &(Slater_inv[row1 * D13_P]);
  double* r2dim = &(Slater_inv[row2 * D13_P]);
  IVDEP
  ALIGNED
  for (uint64_t j = 0; j < D13_P; j++) {
      tmp[j] = Binv[0] * r1dim[j] + Binv[1] * r2dim[j];
      tmp[D13_P + j] = Binv[2] * r1dim[j] + Binv[3] * r2dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : 13 x D13_P
  for (uint64_t i = 0; i < 13; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D13_P; j++) {
          Slater_inv[i * D13_P + j] -= C[i * 2] * tmp[j];
          Slater_inv[i * D13_P + j] -= C[i * 2 + 1] * tmp[D13_P + j];
      }
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_woodbury_2x2_14(
    const qmckl_context context,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {
/*
    C := S^{-1} * U,    dim x 2
    B := 1 + V * C,     2 x 2
    D := V * S^{-1},    2 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
                          QMCKL_NULL_CONTEXT,
                          "qmckl_woodbury_2x2_14",
                          NULL);
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);

  // Compute C = (S^T)^{-1}U : 14 x 2
  double __attribute__((aligned(8))) C[2 * 14];
  for (uint64_t i = 0; i < 14; i++) {
      C[i * 2] = 0;
      C[i * 2 + 1] = 0;
      IVDEP
      ALIGNED
      for (uint64_t k = 0; k < D14_P; k++) {
          C[i * 2] += Slater_inv[i * D14_P + k] * Updates[k];
          C[i * 2 + 1] += Slater_inv[i * D14_P + k] * Updates[D14_P + k];
      }
  }

  // Compute B = 1 + VC : 2 x 2
  const double B0 = C[row1 * 2] + 1;
  const double B1 = C[row1 * 2 + 1];
  const double B2 = C[row2 * 2];
  const double B3 = C[row2 * 2 + 1] + 1;

  // Check if determinant of inverted matrix is not zero
  double det = B0 * B3 - B1 * B2;
  if (fabs(det) < breakdown) {
      return QMCKL_FAILURE;
  }

  // Update det(S) when passed
  if (determinant)
      *determinant *= det;

  // Compute B^{-1} with explicit formula for 2 x 2 inversion
  double __attribute__((aligned(8))) Binv[4], idet = 1.0 / det;
  Binv[0] = idet * B3;
  Binv[1] = -1.0 * idet * B1;
  Binv[2] = -1.0 * idet * B2;
  Binv[3] = idet * B0;

  // tmp = B^{-1}D : 2 x D14_P
  double __attribute__((aligned(8))) tmp[2 * D14_P];
  double* r1dim = &(Slater_inv[row1 * D14_P]);
  double* r2dim = &(Slater_inv[row2 * D14_P]);
  IVDEP
  ALIGNED
  for (uint64_t j = 0; j < D14_P; j++) {
      tmp[j] = Binv[0] * r1dim[j] + Binv[1] * r2dim[j];
      tmp[D14_P + j] = Binv[2] * r1dim[j] + Binv[3] * r2dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : 14 x D14_P
  for (uint64_t i = 0; i < 14; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D14_P; j++) {
          Slater_inv[i * D14_P + j] -= C[i * 2] * tmp[j];
          Slater_inv[i * D14_P + j] -= C[i * 2 + 1] * tmp[D14_P + j];
      }
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_woodbury_2x2_15(
    const qmckl_context context,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {
/*
    C := S^{-1} * U,    dim x 2
    B := 1 + V * C,     2 x 2
    D := V * S^{-1},    2 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
                          QMCKL_NULL_CONTEXT,
                          "qmckl_woodbury_2x2_15",
                          NULL);
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);

  // Compute C = (S^T)^{-1}U : 15 x 2
  double __attribute__((aligned(8))) C[2 * 15];
  for (uint64_t i = 0; i < 15; i++) {
      C[i * 2] = 0;
      C[i * 2 + 1] = 0;
      IVDEP
      ALIGNED
      for (uint64_t k = 0; k < D15_P; k++) {
          C[i * 2] += Slater_inv[i * D15_P + k] * Updates[k];
          C[i * 2 + 1] += Slater_inv[i * D15_P + k] * Updates[D15_P + k];
      }
  }

  // Compute B = 1 + VC : 2 x 2
  const double B0 = C[row1 * 2] + 1;
  const double B1 = C[row1 * 2 + 1];
  const double B2 = C[row2 * 2];
  const double B3 = C[row2 * 2 + 1] + 1;

  // Check if determinant of inverted matrix is not zero
  double det = B0 * B3 - B1 * B2;
  if (fabs(det) < breakdown) {
      return QMCKL_FAILURE;
  }

  // Update det(S) when passed
  if (determinant)
      *determinant *= det;

  // Compute B^{-1} with explicit formula for 2 x 2 inversion
  double __attribute__((aligned(8))) Binv[4], idet = 1.0 / det;
  Binv[0] = idet * B3;
  Binv[1] = -1.0 * idet * B1;
  Binv[2] = -1.0 * idet * B2;
  Binv[3] = idet * B0;

  // tmp = B^{-1}D : 2 x D15_P
  double __attribute__((aligned(8))) tmp[2 * D15_P];
  double* r1dim = &(Slater_inv[row1 * D15_P]);
  double* r2dim = &(Slater_inv[row2 * D15_P]);
  IVDEP
  ALIGNED
  for (uint64_t j = 0; j < D15_P; j++) {
      tmp[j] = Binv[0] * r1dim[j] + Binv[1] * r2dim[j];
      tmp[D15_P + j] = Binv[2] * r1dim[j] + Binv[3] * r2dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : 15 x D15_P
  for (uint64_t i = 0; i < 15; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D15_P; j++) {
          Slater_inv[i * D15_P + j] -= C[i * 2] * tmp[j];
          Slater_inv[i * D15_P + j] -= C[i * 2 + 1] * tmp[D15_P + j];
      }
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_woodbury_2x2_16(
    const qmckl_context context,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {
/*
    C := S^{-1} * U,    dim x 2
    B := 1 + V * C,     2 x 2
    D := V * S^{-1},    2 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
                          QMCKL_NULL_CONTEXT,
                          "qmckl_woodbury_2x2_16",
                          NULL);
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);

  // Compute C = (S^T)^{-1}U : 16 x 2
  double __attribute__((aligned(8))) C[2 * 16];
  for (uint64_t i = 0; i < 16; i++) {
      C[i * 2] = 0;
      C[i * 2 + 1] = 0;
      IVDEP
      ALIGNED
      for (uint64_t k = 0; k < D16_P; k++) {
          C[i * 2] += Slater_inv[i * D16_P + k] * Updates[k];
          C[i * 2 + 1] += Slater_inv[i * D16_P + k] * Updates[D16_P + k];
      }
  }

  // Compute B = 1 + VC : 2 x 2
  const double B0 = C[row1 * 2] + 1;
  const double B1 = C[row1 * 2 + 1];
  const double B2 = C[row2 * 2];
  const double B3 = C[row2 * 2 + 1] + 1;

  // Check if determinant of inverted matrix is not zero
  double det = B0 * B3 - B1 * B2;
  if (fabs(det) < breakdown) {
      return QMCKL_FAILURE;
  }

  // Update det(S) when passed
  if (determinant)
      *determinant *= det;

  // Compute B^{-1} with explicit formula for 2 x 2 inversion
  double __attribute__((aligned(8))) Binv[4], idet = 1.0 / det;
  Binv[0] = idet * B3;
  Binv[1] = -1.0 * idet * B1;
  Binv[2] = -1.0 * idet * B2;
  Binv[3] = idet * B0;

  // tmp = B^{-1}D : 2 x D16_P
  double __attribute__((aligned(8))) tmp[2 * D16_P];
  double* r1dim = &(Slater_inv[row1 * D16_P]);
  double* r2dim = &(Slater_inv[row2 * D16_P]);
  IVDEP
  ALIGNED
  for (uint64_t j = 0; j < D16_P; j++) {
      tmp[j] = Binv[0] * r1dim[j] + Binv[1] * r2dim[j];
      tmp[D16_P + j] = Binv[2] * r1dim[j] + Binv[3] * r2dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : 16 x D16_P
  for (uint64_t i = 0; i < 16; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D16_P; j++) {
          Slater_inv[i * D16_P + j] -= C[i * 2] * tmp[j];
          Slater_inv[i * D16_P + j] -= C[i * 2 + 1] * tmp[D16_P + j];
      }
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_woodbury_2x2_17(
    const qmckl_context context,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {
/*
    C := S^{-1} * U,    dim x 2
    B := 1 + V * C,     2 x 2
    D := V * S^{-1},    2 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
                          QMCKL_NULL_CONTEXT,
                          "qmckl_woodbury_2x2_17",
                          NULL);
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);

  // Compute C = (S^T)^{-1}U : 17 x 2
  double __attribute__((aligned(8))) C[2 * 17];
  for (uint64_t i = 0; i < 17; i++) {
      C[i * 2] = 0;
      C[i * 2 + 1] = 0;
      IVDEP
      ALIGNED
      for (uint64_t k = 0; k < D17_P; k++) {
          C[i * 2] += Slater_inv[i * D17_P + k] * Updates[k];
          C[i * 2 + 1] += Slater_inv[i * D17_P + k] * Updates[D17_P + k];
      }
  }

  // Compute B = 1 + VC : 2 x 2
  const double B0 = C[row1 * 2] + 1;
  const double B1 = C[row1 * 2 + 1];
  const double B2 = C[row2 * 2];
  const double B3 = C[row2 * 2 + 1] + 1;

  // Check if determinant of inverted matrix is not zero
  double det = B0 * B3 - B1 * B2;
  if (fabs(det) < breakdown) {
      return QMCKL_FAILURE;
  }

  // Update det(S) when passed
  if (determinant)
      *determinant *= det;

  // Compute B^{-1} with explicit formula for 2 x 2 inversion
  double __attribute__((aligned(8))) Binv[4], idet = 1.0 / det;
  Binv[0] = idet * B3;
  Binv[1] = -1.0 * idet * B1;
  Binv[2] = -1.0 * idet * B2;
  Binv[3] = idet * B0;

  // tmp = B^{-1}D : 2 x D17_P
  double __attribute__((aligned(8))) tmp[2 * D17_P];
  double* r1dim = &(Slater_inv[row1 * D17_P]);
  double* r2dim = &(Slater_inv[row2 * D17_P]);
  IVDEP
  ALIGNED
  for (uint64_t j = 0; j < D17_P; j++) {
      tmp[j] = Binv[0] * r1dim[j] + Binv[1] * r2dim[j];
      tmp[D17_P + j] = Binv[2] * r1dim[j] + Binv[3] * r2dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : 17 x D17_P
  for (uint64_t i = 0; i < 17; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D17_P; j++) {
          Slater_inv[i * D17_P + j] -= C[i * 2] * tmp[j];
          Slater_inv[i * D17_P + j] -= C[i * 2 + 1] * tmp[D17_P + j];
      }
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_woodbury_2x2_18(
    const qmckl_context context,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {
/*
    C := S^{-1} * U,    dim x 2
    B := 1 + V * C,     2 x 2
    D := V * S^{-1},    2 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
                          QMCKL_NULL_CONTEXT,
                          "qmckl_woodbury_2x2_18",
                          NULL);
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);

  // Compute C = (S^T)^{-1}U : 18 x 2
  double __attribute__((aligned(8))) C[2 * 18];
  for (uint64_t i = 0; i < 18; i++) {
      C[i * 2] = 0;
      C[i * 2 + 1] = 0;
      IVDEP
      ALIGNED
      for (uint64_t k = 0; k < D18_P; k++) {
          C[i * 2] += Slater_inv[i * D18_P + k] * Updates[k];
          C[i * 2 + 1] += Slater_inv[i * D18_P + k] * Updates[D18_P + k];
      }
  }

  // Compute B = 1 + VC : 2 x 2
  const double B0 = C[row1 * 2] + 1;
  const double B1 = C[row1 * 2 + 1];
  const double B2 = C[row2 * 2];
  const double B3 = C[row2 * 2 + 1] + 1;

  // Check if determinant of inverted matrix is not zero
  double det = B0 * B3 - B1 * B2;
  if (fabs(det) < breakdown) {
      return QMCKL_FAILURE;
  }

  // Update det(S) when passed
  if (determinant)
      *determinant *= det;

  // Compute B^{-1} with explicit formula for 2 x 2 inversion
  double __attribute__((aligned(8))) Binv[4], idet = 1.0 / det;
  Binv[0] = idet * B3;
  Binv[1] = -1.0 * idet * B1;
  Binv[2] = -1.0 * idet * B2;
  Binv[3] = idet * B0;

  // tmp = B^{-1}D : 2 x D18_P
  double __attribute__((aligned(8))) tmp[2 * D18_P];
  double* r1dim = &(Slater_inv[row1 * D18_P]);
  double* r2dim = &(Slater_inv[row2 * D18_P]);
  IVDEP
  ALIGNED
  for (uint64_t j = 0; j < D18_P; j++) {
      tmp[j] = Binv[0] * r1dim[j] + Binv[1] * r2dim[j];
      tmp[D18_P + j] = Binv[2] * r1dim[j] + Binv[3] * r2dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : 18 x D18_P
  for (uint64_t i = 0; i < 18; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D18_P; j++) {
          Slater_inv[i * D18_P + j] -= C[i * 2] * tmp[j];
          Slater_inv[i * D18_P + j] -= C[i * 2 + 1] * tmp[D18_P + j];
      }
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_woodbury_2x2_19(
    const qmckl_context context,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {
/*
    C := S^{-1} * U,    dim x 2
    B := 1 + V * C,     2 x 2
    D := V * S^{-1},    2 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
                          QMCKL_NULL_CONTEXT,
                          "qmckl_woodbury_2x2_19",
                          NULL);
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);

  // Compute C = (S^T)^{-1}U : 19 x 2
  double __attribute__((aligned(8))) C[2 * 19];
  for (uint64_t i = 0; i < 19; i++) {
      C[i * 2] = 0;
      C[i * 2 + 1] = 0;
      IVDEP
      ALIGNED
      for (uint64_t k = 0; k < D19_P; k++) {
          C[i * 2] += Slater_inv[i * D19_P + k] * Updates[k];
          C[i * 2 + 1] += Slater_inv[i * D19_P + k] * Updates[D19_P + k];
      }
  }

  // Compute B = 1 + VC : 2 x 2
  const double B0 = C[row1 * 2] + 1;
  const double B1 = C[row1 * 2 + 1];
  const double B2 = C[row2 * 2];
  const double B3 = C[row2 * 2 + 1] + 1;

  // Check if determinant of inverted matrix is not zero
  double det = B0 * B3 - B1 * B2;
  if (fabs(det) < breakdown) {
      return QMCKL_FAILURE;
  }

  // Update det(S) when passed
  if (determinant)
      *determinant *= det;

  // Compute B^{-1} with explicit formula for 2 x 2 inversion
  double __attribute__((aligned(8))) Binv[4], idet = 1.0 / det;
  Binv[0] = idet * B3;
  Binv[1] = -1.0 * idet * B1;
  Binv[2] = -1.0 * idet * B2;
  Binv[3] = idet * B0;

  // tmp = B^{-1}D : 2 x D19_P
  double __attribute__((aligned(8))) tmp[2 * D19_P];
  double* r1dim = &(Slater_inv[row1 * D19_P]);
  double* r2dim = &(Slater_inv[row2 * D19_P]);
  IVDEP
  ALIGNED
  for (uint64_t j = 0; j < D19_P; j++) {
      tmp[j] = Binv[0] * r1dim[j] + Binv[1] * r2dim[j];
      tmp[D19_P + j] = Binv[2] * r1dim[j] + Binv[3] * r2dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : 19 x D19_P
  for (uint64_t i = 0; i < 19; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D19_P; j++) {
          Slater_inv[i * D19_P + j] -= C[i * 2] * tmp[j];
          Slater_inv[i * D19_P + j] -= C[i * 2 + 1] * tmp[D19_P + j];
      }
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_woodbury_2x2_20(
    const qmckl_context context,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {
/*
    C := S^{-1} * U,    dim x 2
    B := 1 + V * C,     2 x 2
    D := V * S^{-1},    2 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
                          QMCKL_NULL_CONTEXT,
                          "qmckl_woodbury_2x2_20",
                          NULL);
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);

  // Compute C = (S^T)^{-1}U : 20 x 2
  double __attribute__((aligned(8))) C[2 * 20];
  for (uint64_t i = 0; i < 20; i++) {
      C[i * 2] = 0;
      C[i * 2 + 1] = 0;
      IVDEP
      ALIGNED
      for (uint64_t k = 0; k < D20_P; k++) {
          C[i * 2] += Slater_inv[i * D20_P + k] * Updates[k];
          C[i * 2 + 1] += Slater_inv[i * D20_P + k] * Updates[D20_P + k];
      }
  }

  // Compute B = 1 + VC : 2 x 2
  const double B0 = C[row1 * 2] + 1;
  const double B1 = C[row1 * 2 + 1];
  const double B2 = C[row2 * 2];
  const double B3 = C[row2 * 2 + 1] + 1;

  // Check if determinant of inverted matrix is not zero
  double det = B0 * B3 - B1 * B2;
  if (fabs(det) < breakdown) {
      return QMCKL_FAILURE;
  }

  // Update det(S) when passed
  if (determinant)
      *determinant *= det;

  // Compute B^{-1} with explicit formula for 2 x 2 inversion
  double __attribute__((aligned(8))) Binv[4], idet = 1.0 / det;
  Binv[0] = idet * B3;
  Binv[1] = -1.0 * idet * B1;
  Binv[2] = -1.0 * idet * B2;
  Binv[3] = idet * B0;

  // tmp = B^{-1}D : 2 x D20_P
  double __attribute__((aligned(8))) tmp[2 * D20_P];
  double* r1dim = &(Slater_inv[row1 * D20_P]);
  double* r2dim = &(Slater_inv[row2 * D20_P]);
  IVDEP
  ALIGNED
  for (uint64_t j = 0; j < D20_P; j++) {
      tmp[j] = Binv[0] * r1dim[j] + Binv[1] * r2dim[j];
      tmp[D20_P + j] = Binv[2] * r1dim[j] + Binv[3] * r2dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : 20 x D20_P
  for (uint64_t i = 0; i < 20; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D20_P; j++) {
          Slater_inv[i * D20_P + j] -= C[i * 2] * tmp[j];
          Slater_inv[i * D20_P + j] -= C[i * 2 + 1] * tmp[D20_P + j];
      }
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_woodbury_2x2_21(
    const qmckl_context context,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {
/*
    C := S^{-1} * U,    dim x 2
    B := 1 + V * C,     2 x 2
    D := V * S^{-1},    2 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
                          QMCKL_NULL_CONTEXT,
                          "qmckl_woodbury_2x2_21",
                          NULL);
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);

  // Compute C = (S^T)^{-1}U : 21 x 2
  double __attribute__((aligned(8))) C[2 * 21];
  for (uint64_t i = 0; i < 21; i++) {
      C[i * 2] = 0;
      C[i * 2 + 1] = 0;
      IVDEP
      ALIGNED
      for (uint64_t k = 0; k < D21_P; k++) {
          C[i * 2] += Slater_inv[i * D21_P + k] * Updates[k];
          C[i * 2 + 1] += Slater_inv[i * D21_P + k] * Updates[D21_P + k];
      }
  }

  // Compute B = 1 + VC : 2 x 2
  const double B0 = C[row1 * 2] + 1;
  const double B1 = C[row1 * 2 + 1];
  const double B2 = C[row2 * 2];
  const double B3 = C[row2 * 2 + 1] + 1;

  // Check if determinant of inverted matrix is not zero
  double det = B0 * B3 - B1 * B2;
  if (fabs(det) < breakdown) {
      return QMCKL_FAILURE;
  }

  // Update det(S) when passed
  if (determinant)
      *determinant *= det;

  // Compute B^{-1} with explicit formula for 2 x 2 inversion
  double __attribute__((aligned(8))) Binv[4], idet = 1.0 / det;
  Binv[0] = idet * B3;
  Binv[1] = -1.0 * idet * B1;
  Binv[2] = -1.0 * idet * B2;
  Binv[3] = idet * B0;

  // tmp = B^{-1}D : 2 x D21_P
  double __attribute__((aligned(8))) tmp[2 * D21_P];
  double* r1dim = &(Slater_inv[row1 * D21_P]);
  double* r2dim = &(Slater_inv[row2 * D21_P]);
  IVDEP
  ALIGNED
  for (uint64_t j = 0; j < D21_P; j++) {
      tmp[j] = Binv[0] * r1dim[j] + Binv[1] * r2dim[j];
      tmp[D21_P + j] = Binv[2] * r1dim[j] + Binv[3] * r2dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : 21 x D21_P
  for (uint64_t i = 0; i < 21; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D21_P; j++) {
          Slater_inv[i * D21_P + j] -= C[i * 2] * tmp[j];
          Slater_inv[i * D21_P + j] -= C[i * 2 + 1] * tmp[D21_P + j];
      }
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_woodbury_2x2(const qmckl_context context,
    const uint64_t LDS,
    const uint64_t Dim,
    const double* Updates,
    const uint64_t* Updates_index,
    const double breakdown,
    double* Slater_inv,
    double* determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(
        context,
        QMCKL_NULL_CONTEXT,
        "qmckl_woodbury_2x2",
        NULL);
  }

  #ifdef HAVE_HPC__BROKEN_WITH_CRAY
  if (LDS == (1+(Dim-1)/SIMD_LENGTH)*SIMD_LENGTH) { // Most cases
    switch (Dim) {
      case 2:  
        return qmckl_woodbury_2x2_2(context,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 3:  
        return qmckl_woodbury_2x2_3(context,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 4:  
        return qmckl_woodbury_2x2_4(context,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 5:  
        return qmckl_woodbury_2x2_5(context,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 6:  
        return qmckl_woodbury_2x2_6(context,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 7:  
        return qmckl_woodbury_2x2_7(context,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 8:  
        return qmckl_woodbury_2x2_8(context,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 9:  
        return qmckl_woodbury_2x2_9(context,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 10:  
        return qmckl_woodbury_2x2_10(context,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 11:  
        return qmckl_woodbury_2x2_11(context,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 12:  
        return qmckl_woodbury_2x2_12(context,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 13:  
        return qmckl_woodbury_2x2_13(context,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 14:  
        return qmckl_woodbury_2x2_14(context,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 15:  
        return qmckl_woodbury_2x2_15(context,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 16:  
        return qmckl_woodbury_2x2_16(context,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 17:  
        return qmckl_woodbury_2x2_17(context,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 18:  
        return qmckl_woodbury_2x2_18(context,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 19:  
        return qmckl_woodbury_2x2_19(context,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 20:  
        return qmckl_woodbury_2x2_20(context,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 21:  
        return qmckl_woodbury_2x2_21(context,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
    }
  }
  else { // Updating smaller sub-matrix
    return qmckl_woodbury_2x2_hpc(
        context,
        LDS,
        Dim,
        Updates,
        Updates_index,
        breakdown,
        Slater_inv,
        determinant);
  }
  #else
  return qmckl_woodbury_2x2_doc(
      context,
      LDS,
      Dim,
      Updates,
      Updates_index,
      breakdown,
      Slater_inv,
      determinant);
  // return qmckl_woodbury_2x2_hpc(
  //     context,
  //     LDS,
  //     Dim,
  //     Updates,
  //     Updates_index,
  //     breakdown,
  //     Slater_inv,
  //     determinant);
  #endif

  return QMCKL_FAILURE;
}

/* C sources */

qmckl_exit_code qmckl_woodbury_3x3_hpc(const qmckl_context context,
    const uint64_t LDS,
    const uint64_t Dim,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {
/*
    C := S^{-1} * U,    dim x 3
    B := 1 + V * C,     3 x 3
    D := V * S^{-1},    3 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
                          QMCKL_NULL_CONTEXT,
                          "qmckl_woodbury_3x3_hpc",
                          NULL);
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);
  const uint64_t row3 = (Updates_index[2] - 1);

  // Compute C = (S^T)^{-1}U : Dim x 3
  double __attribute__((aligned(8))) C[3 * Dim];
  for (uint64_t i = 0; i < Dim; i++) {
      C[i * 3] = 0;
      C[i * 3 + 1] = 0;
      C[i * 3 + 2] = 0;
      IVDEP
      ALIGNED
      for (uint64_t k = 0; k < LDS; k++) {
          C[i * 3]     += Slater_inv[i * LDS + k] * Updates[k];
          C[i * 3 + 1] += Slater_inv[i * LDS + k] * Updates[LDS + k];
          C[i * 3 + 2] += Slater_inv[i * LDS + k] * Updates[2 * LDS + k];
      }
  }

  // Compute B = 1 + VC : 3 x 3
  const double B0 = C[row1 * 3] + 1;
  const double B1 = C[row1 * 3 + 1];
  const double B2 = C[row1 * 3 + 2];
  const double B3 = C[row2 * 3];
  const double B4 = C[row2 * 3 + 1] + 1;
  const double B5 = C[row2 * 3 + 2];
  const double B6 = C[row3 * 3];
  const double B7 = C[row3 * 3 + 1];
  const double B8 = C[row3 * 3 + 2] + 1;

  // Check if determinant of inverted matrix is not zero
  double det;
  det = B0 * (B4 * B8 - B5 * B7) - B1 * (B3 * B8 - B5 * B6) +
        B2 * (B3 * B7 - B4 * B6);
  if (fabs(det) < breakdown) {
      return QMCKL_FAILURE;
  }

  // Update det(S) when passed
  if (determinant)
      *determinant *= det;

  // Compute B^{-1} with explicit formula for 2 x 2 inversion
  double __attribute__((aligned(8))) Binv[9], idet = 1.0 / det;
  Binv[0] =  (B4 * B8 - B7 * B5) * idet;
  Binv[1] = -(B1 * B8 - B7 * B2) * idet;
  Binv[2] =  (B1 * B5 - B4 * B2) * idet;
  Binv[3] = -(B3 * B8 - B6 * B5) * idet;
  Binv[4] =  (B0 * B8 - B6 * B2) * idet;
  Binv[5] = -(B0 * B5 - B3 * B2) * idet;
  Binv[6] =  (B3 * B7 - B6 * B4) * idet;
  Binv[7] = -(B0 * B7 - B6 * B1) * idet;
  Binv[8] =  (B0 * B4 - B3 * B1) * idet;

  // tmp = B^{-1}D : 2 x LDS
  double __attribute__((aligned(8))) tmp[3 * LDS];
  double* r1dim = &(Slater_inv[row1 * LDS]);
  double* r2dim = &(Slater_inv[row2 * LDS]);
  double* r3dim = &(Slater_inv[row3 * LDS]);
  IVDEP
  ALIGNED
  for (uint64_t j = 0; j < LDS; j++) {
      tmp[j] = Binv[0] * r1dim[j] + Binv[1] * r2dim[j] + Binv[2] * r3dim[j];
      tmp[LDS + j] =
          Binv[3] * r1dim[j] + Binv[4] * r2dim[j] + Binv[5] * r3dim[j];
      tmp[2 * LDS + j] =
          Binv[6] * r1dim[j] + Binv[7] * r2dim[j] + Binv[8] * r3dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : Dim x LDS
  for (uint64_t i = 0; i < Dim; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < LDS; j++) {
          Slater_inv[i * LDS + j] -= C[i * 3] * tmp[j];
          Slater_inv[i * LDS + j] -= C[i * 3 + 1] * tmp[LDS + j];
          Slater_inv[i * LDS + j] -= C[i * 3 + 2] * tmp[2 * LDS + j];
      }
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_woodbury_3x3_2(
    const qmckl_context context,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {
/*
    C := S^{-1} * U,    dim x 3
    B := 1 + V * C,     3 x 3
    D := V * S^{-1},    3 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
                          QMCKL_NULL_CONTEXT,
                          "qmckl_woodbury_3x3_2",
                          NULL);
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);
  const uint64_t row3 = (Updates_index[2] - 1);

  // Compute C = (S^T)^{-1}U : 2 x 3
  double __attribute__((aligned(8))) C[3 * 2];
  for (uint64_t i = 0; i < 2; i++) {
      C[i * 3] = 0;
      C[i * 3 + 1] = 0;
      C[i * 3 + 2] = 0;
      IVDEP
      ALIGNED
      for (uint64_t k = 0; k < D2_P; k++) {
          C[i * 3] += Slater_inv[i * D2_P + k] * Updates[k];
          C[i * 3 + 1] += Slater_inv[i * D2_P + k] * Updates[D2_P + k];
          C[i * 3 + 2] += Slater_inv[i * D2_P + k] * Updates[2 * D2_P + k];
      }
  }

  // Compute B = 1 + VC : 3 x 3
  const double B0 = C[row1 * 3] + 1;
  const double B1 = C[row1 * 3 + 1];
  const double B2 = C[row1 * 3 + 2];
  const double B3 = C[row2 * 3];
  const double B4 = C[row2 * 3 + 1] + 1;
  const double B5 = C[row2 * 3 + 2];
  const double B6 = C[row3 * 3];
  const double B7 = C[row3 * 3 + 1];
  const double B8 = C[row3 * 3 + 2] + 1;

  // Check if determinant of B is not too close to zero
  double det;
  det = B0 * (B4 * B8 - B5 * B7) - B1 * (B3 * B8 - B5 * B6) +
        B2 * (B3 * B7 - B4 * B6);
  if (fabs(det) < breakdown) {
      return QMCKL_FAILURE;
  }

  // Update det(Slater) if passed
  if (determinant)
      *determinant *= det;

  // Compute B^{-1} with explicit formula for 3 x 3 inversion
  double __attribute__((aligned(8))) Binv[9], idet = 1.0 / det;
  Binv[0] = (B4 * B8 - B7 * B5) * idet;
  Binv[1] = -(B1 * B8 - B7 * B2) * idet;
  Binv[2] = (B1 * B5 - B4 * B2) * idet;
  Binv[3] = -(B3 * B8 - B6 * B5) * idet;
  Binv[4] = (B0 * B8 - B6 * B2) * idet;
  Binv[5] = -(B0 * B5 - B3 * B2) * idet;
  Binv[6] = (B3 * B7 - B6 * B4) * idet;
  Binv[7] = -(B0 * B7 - B6 * B1) * idet;
  Binv[8] = (B0 * B4 - B3 * B1) * idet;

  // tmp = B^{-1}D : 3 x D2_P
  double __attribute__((aligned(8))) tmp[3 * D2_P];
  double* r1dim = &(Slater_inv[row1 * D2_P]);
  double* r2dim = &(Slater_inv[row2 * D2_P]);
  double* r3dim = &(Slater_inv[row3 * D2_P]);
      IVDEP
      ALIGNED
  for (uint64_t j = 0; j < D2_P; j++) {
      tmp[j] = Binv[0] * r1dim[j] + Binv[1] * r2dim[j] + Binv[2] * r3dim[j];
      tmp[D2_P + j] =
          Binv[3] * r1dim[j] + Binv[4] * r2dim[j] + Binv[5] * r3dim[j];
      tmp[2 * D2_P + j] =
          Binv[6] * r1dim[j] + Binv[7] * r2dim[j] + Binv[8] * r3dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : 2 x D2_P
  for (uint64_t i = 0; i < 2; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D2_P; j++) {
          Slater_inv[i * D2_P + j] -= C[i * 3] * tmp[j];
          Slater_inv[i * D2_P + j] -= C[i * 3 + 1] * tmp[D2_P + j];
          Slater_inv[i * D2_P + j] -= C[i * 3 + 2] * tmp[2 * D2_P + j];
      }
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_woodbury_3x3_3(
    const qmckl_context context,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {
/*
    C := S^{-1} * U,    dim x 3
    B := 1 + V * C,     3 x 3
    D := V * S^{-1},    3 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
                          QMCKL_NULL_CONTEXT,
                          "qmckl_woodbury_3x3_3",
                          NULL);
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);
  const uint64_t row3 = (Updates_index[2] - 1);

  // Compute C = (S^T)^{-1}U : 3 x 3
  double __attribute__((aligned(8))) C[3 * 3];
  for (uint64_t i = 0; i < 3; i++) {
      C[i * 3] = 0;
      C[i * 3 + 1] = 0;
      C[i * 3 + 2] = 0;
      IVDEP
      ALIGNED
      for (uint64_t k = 0; k < D3_P; k++) {
          C[i * 3] += Slater_inv[i * D3_P + k] * Updates[k];
          C[i * 3 + 1] += Slater_inv[i * D3_P + k] * Updates[D3_P + k];
          C[i * 3 + 2] += Slater_inv[i * D3_P + k] * Updates[2 * D3_P + k];
      }
  }

  // Compute B = 1 + VC : 3 x 3
  const double B0 = C[row1 * 3] + 1;
  const double B1 = C[row1 * 3 + 1];
  const double B2 = C[row1 * 3 + 2];
  const double B3 = C[row2 * 3];
  const double B4 = C[row2 * 3 + 1] + 1;
  const double B5 = C[row2 * 3 + 2];
  const double B6 = C[row3 * 3];
  const double B7 = C[row3 * 3 + 1];
  const double B8 = C[row3 * 3 + 2] + 1;

  // Check if determinant of B is not too close to zero
  double det;
  det = B0 * (B4 * B8 - B5 * B7) - B1 * (B3 * B8 - B5 * B6) +
        B2 * (B3 * B7 - B4 * B6);
  if (fabs(det) < breakdown) {
      return QMCKL_FAILURE;
  }

  // Update det(Slater) if passed
  if (determinant)
      *determinant *= det;

  // Compute B^{-1} with explicit formula for 3 x 3 inversion
  double __attribute__((aligned(8))) Binv[9], idet = 1.0 / det;
  Binv[0] = (B4 * B8 - B7 * B5) * idet;
  Binv[1] = -(B1 * B8 - B7 * B2) * idet;
  Binv[2] = (B1 * B5 - B4 * B2) * idet;
  Binv[3] = -(B3 * B8 - B6 * B5) * idet;
  Binv[4] = (B0 * B8 - B6 * B2) * idet;
  Binv[5] = -(B0 * B5 - B3 * B2) * idet;
  Binv[6] = (B3 * B7 - B6 * B4) * idet;
  Binv[7] = -(B0 * B7 - B6 * B1) * idet;
  Binv[8] = (B0 * B4 - B3 * B1) * idet;

  // tmp = B^{-1}D : 3 x D3_P
  double __attribute__((aligned(8))) tmp[3 * D3_P];
  double* r1dim = &(Slater_inv[row1 * D3_P]);
  double* r2dim = &(Slater_inv[row2 * D3_P]);
  double* r3dim = &(Slater_inv[row3 * D3_P]);
      IVDEP
      ALIGNED
  for (uint64_t j = 0; j < D3_P; j++) {
      tmp[j] = Binv[0] * r1dim[j] + Binv[1] * r2dim[j] + Binv[2] * r3dim[j];
      tmp[D3_P + j] =
          Binv[3] * r1dim[j] + Binv[4] * r2dim[j] + Binv[5] * r3dim[j];
      tmp[2 * D3_P + j] =
          Binv[6] * r1dim[j] + Binv[7] * r2dim[j] + Binv[8] * r3dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : 3 x D3_P
  for (uint64_t i = 0; i < 3; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D3_P; j++) {
          Slater_inv[i * D3_P + j] -= C[i * 3] * tmp[j];
          Slater_inv[i * D3_P + j] -= C[i * 3 + 1] * tmp[D3_P + j];
          Slater_inv[i * D3_P + j] -= C[i * 3 + 2] * tmp[2 * D3_P + j];
      }
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_woodbury_3x3_4(
    const qmckl_context context,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {
/*
    C := S^{-1} * U,    dim x 3
    B := 1 + V * C,     3 x 3
    D := V * S^{-1},    3 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
                          QMCKL_NULL_CONTEXT,
                          "qmckl_woodbury_3x3_4",
                          NULL);
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);
  const uint64_t row3 = (Updates_index[2] - 1);

  // Compute C = (S^T)^{-1}U : 4 x 3
  double __attribute__((aligned(8))) C[3 * 4];
  for (uint64_t i = 0; i < 4; i++) {
      C[i * 3] = 0;
      C[i * 3 + 1] = 0;
      C[i * 3 + 2] = 0;
      IVDEP
      ALIGNED
      for (uint64_t k = 0; k < D4_P; k++) {
          C[i * 3] += Slater_inv[i * D4_P + k] * Updates[k];
          C[i * 3 + 1] += Slater_inv[i * D4_P + k] * Updates[D4_P + k];
          C[i * 3 + 2] += Slater_inv[i * D4_P + k] * Updates[2 * D4_P + k];
      }
  }

  // Compute B = 1 + VC : 3 x 3
  const double B0 = C[row1 * 3] + 1;
  const double B1 = C[row1 * 3 + 1];
  const double B2 = C[row1 * 3 + 2];
  const double B3 = C[row2 * 3];
  const double B4 = C[row2 * 3 + 1] + 1;
  const double B5 = C[row2 * 3 + 2];
  const double B6 = C[row3 * 3];
  const double B7 = C[row3 * 3 + 1];
  const double B8 = C[row3 * 3 + 2] + 1;

  // Check if determinant of B is not too close to zero
  double det;
  det = B0 * (B4 * B8 - B5 * B7) - B1 * (B3 * B8 - B5 * B6) +
        B2 * (B3 * B7 - B4 * B6);
  if (fabs(det) < breakdown) {
      return QMCKL_FAILURE;
  }

  // Update det(Slater) if passed
  if (determinant)
      *determinant *= det;

  // Compute B^{-1} with explicit formula for 3 x 3 inversion
  double __attribute__((aligned(8))) Binv[9], idet = 1.0 / det;
  Binv[0] = (B4 * B8 - B7 * B5) * idet;
  Binv[1] = -(B1 * B8 - B7 * B2) * idet;
  Binv[2] = (B1 * B5 - B4 * B2) * idet;
  Binv[3] = -(B3 * B8 - B6 * B5) * idet;
  Binv[4] = (B0 * B8 - B6 * B2) * idet;
  Binv[5] = -(B0 * B5 - B3 * B2) * idet;
  Binv[6] = (B3 * B7 - B6 * B4) * idet;
  Binv[7] = -(B0 * B7 - B6 * B1) * idet;
  Binv[8] = (B0 * B4 - B3 * B1) * idet;

  // tmp = B^{-1}D : 3 x D4_P
  double __attribute__((aligned(8))) tmp[3 * D4_P];
  double* r1dim = &(Slater_inv[row1 * D4_P]);
  double* r2dim = &(Slater_inv[row2 * D4_P]);
  double* r3dim = &(Slater_inv[row3 * D4_P]);
      IVDEP
      ALIGNED
  for (uint64_t j = 0; j < D4_P; j++) {
      tmp[j] = Binv[0] * r1dim[j] + Binv[1] * r2dim[j] + Binv[2] * r3dim[j];
      tmp[D4_P + j] =
          Binv[3] * r1dim[j] + Binv[4] * r2dim[j] + Binv[5] * r3dim[j];
      tmp[2 * D4_P + j] =
          Binv[6] * r1dim[j] + Binv[7] * r2dim[j] + Binv[8] * r3dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : 4 x D4_P
  for (uint64_t i = 0; i < 4; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D4_P; j++) {
          Slater_inv[i * D4_P + j] -= C[i * 3] * tmp[j];
          Slater_inv[i * D4_P + j] -= C[i * 3 + 1] * tmp[D4_P + j];
          Slater_inv[i * D4_P + j] -= C[i * 3 + 2] * tmp[2 * D4_P + j];
      }
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_woodbury_3x3_5(
    const qmckl_context context,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {
/*
    C := S^{-1} * U,    dim x 3
    B := 1 + V * C,     3 x 3
    D := V * S^{-1},    3 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
                          QMCKL_NULL_CONTEXT,
                          "qmckl_woodbury_3x3_5",
                          NULL);
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);
  const uint64_t row3 = (Updates_index[2] - 1);

  // Compute C = (S^T)^{-1}U : 5 x 3
  double __attribute__((aligned(8))) C[3 * 5];
  for (uint64_t i = 0; i < 5; i++) {
      C[i * 3] = 0;
      C[i * 3 + 1] = 0;
      C[i * 3 + 2] = 0;
      IVDEP
      ALIGNED
      for (uint64_t k = 0; k < D5_P; k++) {
          C[i * 3] += Slater_inv[i * D5_P + k] * Updates[k];
          C[i * 3 + 1] += Slater_inv[i * D5_P + k] * Updates[D5_P + k];
          C[i * 3 + 2] += Slater_inv[i * D5_P + k] * Updates[2 * D5_P + k];
      }
  }

  // Compute B = 1 + VC : 3 x 3
  const double B0 = C[row1 * 3] + 1;
  const double B1 = C[row1 * 3 + 1];
  const double B2 = C[row1 * 3 + 2];
  const double B3 = C[row2 * 3];
  const double B4 = C[row2 * 3 + 1] + 1;
  const double B5 = C[row2 * 3 + 2];
  const double B6 = C[row3 * 3];
  const double B7 = C[row3 * 3 + 1];
  const double B8 = C[row3 * 3 + 2] + 1;

  // Check if determinant of B is not too close to zero
  double det;
  det = B0 * (B4 * B8 - B5 * B7) - B1 * (B3 * B8 - B5 * B6) +
        B2 * (B3 * B7 - B4 * B6);
  if (fabs(det) < breakdown) {
      return QMCKL_FAILURE;
  }

  // Update det(Slater) if passed
  if (determinant)
      *determinant *= det;

  // Compute B^{-1} with explicit formula for 3 x 3 inversion
  double __attribute__((aligned(8))) Binv[9], idet = 1.0 / det;
  Binv[0] = (B4 * B8 - B7 * B5) * idet;
  Binv[1] = -(B1 * B8 - B7 * B2) * idet;
  Binv[2] = (B1 * B5 - B4 * B2) * idet;
  Binv[3] = -(B3 * B8 - B6 * B5) * idet;
  Binv[4] = (B0 * B8 - B6 * B2) * idet;
  Binv[5] = -(B0 * B5 - B3 * B2) * idet;
  Binv[6] = (B3 * B7 - B6 * B4) * idet;
  Binv[7] = -(B0 * B7 - B6 * B1) * idet;
  Binv[8] = (B0 * B4 - B3 * B1) * idet;

  // tmp = B^{-1}D : 3 x D5_P
  double __attribute__((aligned(8))) tmp[3 * D5_P];
  double* r1dim = &(Slater_inv[row1 * D5_P]);
  double* r2dim = &(Slater_inv[row2 * D5_P]);
  double* r3dim = &(Slater_inv[row3 * D5_P]);
      IVDEP
      ALIGNED
  for (uint64_t j = 0; j < D5_P; j++) {
      tmp[j] = Binv[0] * r1dim[j] + Binv[1] * r2dim[j] + Binv[2] * r3dim[j];
      tmp[D5_P + j] =
          Binv[3] * r1dim[j] + Binv[4] * r2dim[j] + Binv[5] * r3dim[j];
      tmp[2 * D5_P + j] =
          Binv[6] * r1dim[j] + Binv[7] * r2dim[j] + Binv[8] * r3dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : 5 x D5_P
  for (uint64_t i = 0; i < 5; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D5_P; j++) {
          Slater_inv[i * D5_P + j] -= C[i * 3] * tmp[j];
          Slater_inv[i * D5_P + j] -= C[i * 3 + 1] * tmp[D5_P + j];
          Slater_inv[i * D5_P + j] -= C[i * 3 + 2] * tmp[2 * D5_P + j];
      }
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_woodbury_3x3_6(
    const qmckl_context context,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {
/*
    C := S^{-1} * U,    dim x 3
    B := 1 + V * C,     3 x 3
    D := V * S^{-1},    3 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
                          QMCKL_NULL_CONTEXT,
                          "qmckl_woodbury_3x3_6",
                          NULL);
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);
  const uint64_t row3 = (Updates_index[2] - 1);

  // Compute C = (S^T)^{-1}U : 6 x 3
  double __attribute__((aligned(8))) C[3 * 6];
  for (uint64_t i = 0; i < 6; i++) {
      C[i * 3] = 0;
      C[i * 3 + 1] = 0;
      C[i * 3 + 2] = 0;
      IVDEP
      ALIGNED
      for (uint64_t k = 0; k < D6_P; k++) {
          C[i * 3] += Slater_inv[i * D6_P + k] * Updates[k];
          C[i * 3 + 1] += Slater_inv[i * D6_P + k] * Updates[D6_P + k];
          C[i * 3 + 2] += Slater_inv[i * D6_P + k] * Updates[2 * D6_P + k];
      }
  }

  // Compute B = 1 + VC : 3 x 3
  const double B0 = C[row1 * 3] + 1;
  const double B1 = C[row1 * 3 + 1];
  const double B2 = C[row1 * 3 + 2];
  const double B3 = C[row2 * 3];
  const double B4 = C[row2 * 3 + 1] + 1;
  const double B5 = C[row2 * 3 + 2];
  const double B6 = C[row3 * 3];
  const double B7 = C[row3 * 3 + 1];
  const double B8 = C[row3 * 3 + 2] + 1;

  // Check if determinant of B is not too close to zero
  double det;
  det = B0 * (B4 * B8 - B5 * B7) - B1 * (B3 * B8 - B5 * B6) +
        B2 * (B3 * B7 - B4 * B6);
  if (fabs(det) < breakdown) {
      return QMCKL_FAILURE;
  }

  // Update det(Slater) if passed
  if (determinant)
      *determinant *= det;

  // Compute B^{-1} with explicit formula for 3 x 3 inversion
  double __attribute__((aligned(8))) Binv[9], idet = 1.0 / det;
  Binv[0] = (B4 * B8 - B7 * B5) * idet;
  Binv[1] = -(B1 * B8 - B7 * B2) * idet;
  Binv[2] = (B1 * B5 - B4 * B2) * idet;
  Binv[3] = -(B3 * B8 - B6 * B5) * idet;
  Binv[4] = (B0 * B8 - B6 * B2) * idet;
  Binv[5] = -(B0 * B5 - B3 * B2) * idet;
  Binv[6] = (B3 * B7 - B6 * B4) * idet;
  Binv[7] = -(B0 * B7 - B6 * B1) * idet;
  Binv[8] = (B0 * B4 - B3 * B1) * idet;

  // tmp = B^{-1}D : 3 x D6_P
  double __attribute__((aligned(8))) tmp[3 * D6_P];
  double* r1dim = &(Slater_inv[row1 * D6_P]);
  double* r2dim = &(Slater_inv[row2 * D6_P]);
  double* r3dim = &(Slater_inv[row3 * D6_P]);
      IVDEP
      ALIGNED
  for (uint64_t j = 0; j < D6_P; j++) {
      tmp[j] = Binv[0] * r1dim[j] + Binv[1] * r2dim[j] + Binv[2] * r3dim[j];
      tmp[D6_P + j] =
          Binv[3] * r1dim[j] + Binv[4] * r2dim[j] + Binv[5] * r3dim[j];
      tmp[2 * D6_P + j] =
          Binv[6] * r1dim[j] + Binv[7] * r2dim[j] + Binv[8] * r3dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : 6 x D6_P
  for (uint64_t i = 0; i < 6; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D6_P; j++) {
          Slater_inv[i * D6_P + j] -= C[i * 3] * tmp[j];
          Slater_inv[i * D6_P + j] -= C[i * 3 + 1] * tmp[D6_P + j];
          Slater_inv[i * D6_P + j] -= C[i * 3 + 2] * tmp[2 * D6_P + j];
      }
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_woodbury_3x3_7(
    const qmckl_context context,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {
/*
    C := S^{-1} * U,    dim x 3
    B := 1 + V * C,     3 x 3
    D := V * S^{-1},    3 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
                          QMCKL_NULL_CONTEXT,
                          "qmckl_woodbury_3x3_7",
                          NULL);
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);
  const uint64_t row3 = (Updates_index[2] - 1);

  // Compute C = (S^T)^{-1}U : 7 x 3
  double __attribute__((aligned(8))) C[3 * 7];
  for (uint64_t i = 0; i < 7; i++) {
      C[i * 3] = 0;
      C[i * 3 + 1] = 0;
      C[i * 3 + 2] = 0;
      IVDEP
      ALIGNED
      for (uint64_t k = 0; k < D7_P; k++) {
          C[i * 3] += Slater_inv[i * D7_P + k] * Updates[k];
          C[i * 3 + 1] += Slater_inv[i * D7_P + k] * Updates[D7_P + k];
          C[i * 3 + 2] += Slater_inv[i * D7_P + k] * Updates[2 * D7_P + k];
      }
  }

  // Compute B = 1 + VC : 3 x 3
  const double B0 = C[row1 * 3] + 1;
  const double B1 = C[row1 * 3 + 1];
  const double B2 = C[row1 * 3 + 2];
  const double B3 = C[row2 * 3];
  const double B4 = C[row2 * 3 + 1] + 1;
  const double B5 = C[row2 * 3 + 2];
  const double B6 = C[row3 * 3];
  const double B7 = C[row3 * 3 + 1];
  const double B8 = C[row3 * 3 + 2] + 1;

  // Check if determinant of B is not too close to zero
  double det;
  det = B0 * (B4 * B8 - B5 * B7) - B1 * (B3 * B8 - B5 * B6) +
        B2 * (B3 * B7 - B4 * B6);
  if (fabs(det) < breakdown) {
      return QMCKL_FAILURE;
  }

  // Update det(Slater) if passed
  if (determinant)
      *determinant *= det;

  // Compute B^{-1} with explicit formula for 3 x 3 inversion
  double __attribute__((aligned(8))) Binv[9], idet = 1.0 / det;
  Binv[0] = (B4 * B8 - B7 * B5) * idet;
  Binv[1] = -(B1 * B8 - B7 * B2) * idet;
  Binv[2] = (B1 * B5 - B4 * B2) * idet;
  Binv[3] = -(B3 * B8 - B6 * B5) * idet;
  Binv[4] = (B0 * B8 - B6 * B2) * idet;
  Binv[5] = -(B0 * B5 - B3 * B2) * idet;
  Binv[6] = (B3 * B7 - B6 * B4) * idet;
  Binv[7] = -(B0 * B7 - B6 * B1) * idet;
  Binv[8] = (B0 * B4 - B3 * B1) * idet;

  // tmp = B^{-1}D : 3 x D7_P
  double __attribute__((aligned(8))) tmp[3 * D7_P];
  double* r1dim = &(Slater_inv[row1 * D7_P]);
  double* r2dim = &(Slater_inv[row2 * D7_P]);
  double* r3dim = &(Slater_inv[row3 * D7_P]);
      IVDEP
      ALIGNED
  for (uint64_t j = 0; j < D7_P; j++) {
      tmp[j] = Binv[0] * r1dim[j] + Binv[1] * r2dim[j] + Binv[2] * r3dim[j];
      tmp[D7_P + j] =
          Binv[3] * r1dim[j] + Binv[4] * r2dim[j] + Binv[5] * r3dim[j];
      tmp[2 * D7_P + j] =
          Binv[6] * r1dim[j] + Binv[7] * r2dim[j] + Binv[8] * r3dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : 7 x D7_P
  for (uint64_t i = 0; i < 7; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D7_P; j++) {
          Slater_inv[i * D7_P + j] -= C[i * 3] * tmp[j];
          Slater_inv[i * D7_P + j] -= C[i * 3 + 1] * tmp[D7_P + j];
          Slater_inv[i * D7_P + j] -= C[i * 3 + 2] * tmp[2 * D7_P + j];
      }
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_woodbury_3x3_8(
    const qmckl_context context,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {
/*
    C := S^{-1} * U,    dim x 3
    B := 1 + V * C,     3 x 3
    D := V * S^{-1},    3 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
                          QMCKL_NULL_CONTEXT,
                          "qmckl_woodbury_3x3_8",
                          NULL);
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);
  const uint64_t row3 = (Updates_index[2] - 1);

  // Compute C = (S^T)^{-1}U : 8 x 3
  double __attribute__((aligned(8))) C[3 * 8];
  for (uint64_t i = 0; i < 8; i++) {
      C[i * 3] = 0;
      C[i * 3 + 1] = 0;
      C[i * 3 + 2] = 0;
      IVDEP
      ALIGNED
      for (uint64_t k = 0; k < D8_P; k++) {
          C[i * 3] += Slater_inv[i * D8_P + k] * Updates[k];
          C[i * 3 + 1] += Slater_inv[i * D8_P + k] * Updates[D8_P + k];
          C[i * 3 + 2] += Slater_inv[i * D8_P + k] * Updates[2 * D8_P + k];
      }
  }

  // Compute B = 1 + VC : 3 x 3
  const double B0 = C[row1 * 3] + 1;
  const double B1 = C[row1 * 3 + 1];
  const double B2 = C[row1 * 3 + 2];
  const double B3 = C[row2 * 3];
  const double B4 = C[row2 * 3 + 1] + 1;
  const double B5 = C[row2 * 3 + 2];
  const double B6 = C[row3 * 3];
  const double B7 = C[row3 * 3 + 1];
  const double B8 = C[row3 * 3 + 2] + 1;

  // Check if determinant of B is not too close to zero
  double det;
  det = B0 * (B4 * B8 - B5 * B7) - B1 * (B3 * B8 - B5 * B6) +
        B2 * (B3 * B7 - B4 * B6);
  if (fabs(det) < breakdown) {
      return QMCKL_FAILURE;
  }

  // Update det(Slater) if passed
  if (determinant)
      *determinant *= det;

  // Compute B^{-1} with explicit formula for 3 x 3 inversion
  double __attribute__((aligned(8))) Binv[9], idet = 1.0 / det;
  Binv[0] = (B4 * B8 - B7 * B5) * idet;
  Binv[1] = -(B1 * B8 - B7 * B2) * idet;
  Binv[2] = (B1 * B5 - B4 * B2) * idet;
  Binv[3] = -(B3 * B8 - B6 * B5) * idet;
  Binv[4] = (B0 * B8 - B6 * B2) * idet;
  Binv[5] = -(B0 * B5 - B3 * B2) * idet;
  Binv[6] = (B3 * B7 - B6 * B4) * idet;
  Binv[7] = -(B0 * B7 - B6 * B1) * idet;
  Binv[8] = (B0 * B4 - B3 * B1) * idet;

  // tmp = B^{-1}D : 3 x D8_P
  double __attribute__((aligned(8))) tmp[3 * D8_P];
  double* r1dim = &(Slater_inv[row1 * D8_P]);
  double* r2dim = &(Slater_inv[row2 * D8_P]);
  double* r3dim = &(Slater_inv[row3 * D8_P]);
      IVDEP
      ALIGNED
  for (uint64_t j = 0; j < D8_P; j++) {
      tmp[j] = Binv[0] * r1dim[j] + Binv[1] * r2dim[j] + Binv[2] * r3dim[j];
      tmp[D8_P + j] =
          Binv[3] * r1dim[j] + Binv[4] * r2dim[j] + Binv[5] * r3dim[j];
      tmp[2 * D8_P + j] =
          Binv[6] * r1dim[j] + Binv[7] * r2dim[j] + Binv[8] * r3dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : 8 x D8_P
  for (uint64_t i = 0; i < 8; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D8_P; j++) {
          Slater_inv[i * D8_P + j] -= C[i * 3] * tmp[j];
          Slater_inv[i * D8_P + j] -= C[i * 3 + 1] * tmp[D8_P + j];
          Slater_inv[i * D8_P + j] -= C[i * 3 + 2] * tmp[2 * D8_P + j];
      }
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_woodbury_3x3_9(
    const qmckl_context context,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {
/*
    C := S^{-1} * U,    dim x 3
    B := 1 + V * C,     3 x 3
    D := V * S^{-1},    3 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
                          QMCKL_NULL_CONTEXT,
                          "qmckl_woodbury_3x3_9",
                          NULL);
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);
  const uint64_t row3 = (Updates_index[2] - 1);

  // Compute C = (S^T)^{-1}U : 9 x 3
  double __attribute__((aligned(8))) C[3 * 9];
  for (uint64_t i = 0; i < 9; i++) {
      C[i * 3] = 0;
      C[i * 3 + 1] = 0;
      C[i * 3 + 2] = 0;
      IVDEP
      ALIGNED
      for (uint64_t k = 0; k < D9_P; k++) {
          C[i * 3] += Slater_inv[i * D9_P + k] * Updates[k];
          C[i * 3 + 1] += Slater_inv[i * D9_P + k] * Updates[D9_P + k];
          C[i * 3 + 2] += Slater_inv[i * D9_P + k] * Updates[2 * D9_P + k];
      }
  }

  // Compute B = 1 + VC : 3 x 3
  const double B0 = C[row1 * 3] + 1;
  const double B1 = C[row1 * 3 + 1];
  const double B2 = C[row1 * 3 + 2];
  const double B3 = C[row2 * 3];
  const double B4 = C[row2 * 3 + 1] + 1;
  const double B5 = C[row2 * 3 + 2];
  const double B6 = C[row3 * 3];
  const double B7 = C[row3 * 3 + 1];
  const double B8 = C[row3 * 3 + 2] + 1;

  // Check if determinant of B is not too close to zero
  double det;
  det = B0 * (B4 * B8 - B5 * B7) - B1 * (B3 * B8 - B5 * B6) +
        B2 * (B3 * B7 - B4 * B6);
  if (fabs(det) < breakdown) {
      return QMCKL_FAILURE;
  }

  // Update det(Slater) if passed
  if (determinant)
      *determinant *= det;

  // Compute B^{-1} with explicit formula for 3 x 3 inversion
  double __attribute__((aligned(8))) Binv[9], idet = 1.0 / det;
  Binv[0] = (B4 * B8 - B7 * B5) * idet;
  Binv[1] = -(B1 * B8 - B7 * B2) * idet;
  Binv[2] = (B1 * B5 - B4 * B2) * idet;
  Binv[3] = -(B3 * B8 - B6 * B5) * idet;
  Binv[4] = (B0 * B8 - B6 * B2) * idet;
  Binv[5] = -(B0 * B5 - B3 * B2) * idet;
  Binv[6] = (B3 * B7 - B6 * B4) * idet;
  Binv[7] = -(B0 * B7 - B6 * B1) * idet;
  Binv[8] = (B0 * B4 - B3 * B1) * idet;

  // tmp = B^{-1}D : 3 x D9_P
  double __attribute__((aligned(8))) tmp[3 * D9_P];
  double* r1dim = &(Slater_inv[row1 * D9_P]);
  double* r2dim = &(Slater_inv[row2 * D9_P]);
  double* r3dim = &(Slater_inv[row3 * D9_P]);
      IVDEP
      ALIGNED
  for (uint64_t j = 0; j < D9_P; j++) {
      tmp[j] = Binv[0] * r1dim[j] + Binv[1] * r2dim[j] + Binv[2] * r3dim[j];
      tmp[D9_P + j] =
          Binv[3] * r1dim[j] + Binv[4] * r2dim[j] + Binv[5] * r3dim[j];
      tmp[2 * D9_P + j] =
          Binv[6] * r1dim[j] + Binv[7] * r2dim[j] + Binv[8] * r3dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : 9 x D9_P
  for (uint64_t i = 0; i < 9; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D9_P; j++) {
          Slater_inv[i * D9_P + j] -= C[i * 3] * tmp[j];
          Slater_inv[i * D9_P + j] -= C[i * 3 + 1] * tmp[D9_P + j];
          Slater_inv[i * D9_P + j] -= C[i * 3 + 2] * tmp[2 * D9_P + j];
      }
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_woodbury_3x3_10(
    const qmckl_context context,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {
/*
    C := S^{-1} * U,    dim x 3
    B := 1 + V * C,     3 x 3
    D := V * S^{-1},    3 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
                          QMCKL_NULL_CONTEXT,
                          "qmckl_woodbury_3x3_10",
                          NULL);
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);
  const uint64_t row3 = (Updates_index[2] - 1);

  // Compute C = (S^T)^{-1}U : 10 x 3
  double __attribute__((aligned(8))) C[3 * 10];
  for (uint64_t i = 0; i < 10; i++) {
      C[i * 3] = 0;
      C[i * 3 + 1] = 0;
      C[i * 3 + 2] = 0;
      IVDEP
      ALIGNED
      for (uint64_t k = 0; k < D10_P; k++) {
          C[i * 3] += Slater_inv[i * D10_P + k] * Updates[k];
          C[i * 3 + 1] += Slater_inv[i * D10_P + k] * Updates[D10_P + k];
          C[i * 3 + 2] += Slater_inv[i * D10_P + k] * Updates[2 * D10_P + k];
      }
  }

  // Compute B = 1 + VC : 3 x 3
  const double B0 = C[row1 * 3] + 1;
  const double B1 = C[row1 * 3 + 1];
  const double B2 = C[row1 * 3 + 2];
  const double B3 = C[row2 * 3];
  const double B4 = C[row2 * 3 + 1] + 1;
  const double B5 = C[row2 * 3 + 2];
  const double B6 = C[row3 * 3];
  const double B7 = C[row3 * 3 + 1];
  const double B8 = C[row3 * 3 + 2] + 1;

  // Check if determinant of B is not too close to zero
  double det;
  det = B0 * (B4 * B8 - B5 * B7) - B1 * (B3 * B8 - B5 * B6) +
        B2 * (B3 * B7 - B4 * B6);
  if (fabs(det) < breakdown) {
      return QMCKL_FAILURE;
  }

  // Update det(Slater) if passed
  if (determinant)
      *determinant *= det;

  // Compute B^{-1} with explicit formula for 3 x 3 inversion
  double __attribute__((aligned(8))) Binv[9], idet = 1.0 / det;
  Binv[0] = (B4 * B8 - B7 * B5) * idet;
  Binv[1] = -(B1 * B8 - B7 * B2) * idet;
  Binv[2] = (B1 * B5 - B4 * B2) * idet;
  Binv[3] = -(B3 * B8 - B6 * B5) * idet;
  Binv[4] = (B0 * B8 - B6 * B2) * idet;
  Binv[5] = -(B0 * B5 - B3 * B2) * idet;
  Binv[6] = (B3 * B7 - B6 * B4) * idet;
  Binv[7] = -(B0 * B7 - B6 * B1) * idet;
  Binv[8] = (B0 * B4 - B3 * B1) * idet;

  // tmp = B^{-1}D : 3 x D10_P
  double __attribute__((aligned(8))) tmp[3 * D10_P];
  double* r1dim = &(Slater_inv[row1 * D10_P]);
  double* r2dim = &(Slater_inv[row2 * D10_P]);
  double* r3dim = &(Slater_inv[row3 * D10_P]);
      IVDEP
      ALIGNED
  for (uint64_t j = 0; j < D10_P; j++) {
      tmp[j] = Binv[0] * r1dim[j] + Binv[1] * r2dim[j] + Binv[2] * r3dim[j];
      tmp[D10_P + j] =
          Binv[3] * r1dim[j] + Binv[4] * r2dim[j] + Binv[5] * r3dim[j];
      tmp[2 * D10_P + j] =
          Binv[6] * r1dim[j] + Binv[7] * r2dim[j] + Binv[8] * r3dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : 10 x D10_P
  for (uint64_t i = 0; i < 10; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D10_P; j++) {
          Slater_inv[i * D10_P + j] -= C[i * 3] * tmp[j];
          Slater_inv[i * D10_P + j] -= C[i * 3 + 1] * tmp[D10_P + j];
          Slater_inv[i * D10_P + j] -= C[i * 3 + 2] * tmp[2 * D10_P + j];
      }
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_woodbury_3x3_11(
    const qmckl_context context,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {
/*
    C := S^{-1} * U,    dim x 3
    B := 1 + V * C,     3 x 3
    D := V * S^{-1},    3 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
                          QMCKL_NULL_CONTEXT,
                          "qmckl_woodbury_3x3_11",
                          NULL);
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);
  const uint64_t row3 = (Updates_index[2] - 1);

  // Compute C = (S^T)^{-1}U : 11 x 3
  double __attribute__((aligned(8))) C[3 * 11];
  for (uint64_t i = 0; i < 11; i++) {
      C[i * 3] = 0;
      C[i * 3 + 1] = 0;
      C[i * 3 + 2] = 0;
      IVDEP
      ALIGNED
      for (uint64_t k = 0; k < D11_P; k++) {
          C[i * 3] += Slater_inv[i * D11_P + k] * Updates[k];
          C[i * 3 + 1] += Slater_inv[i * D11_P + k] * Updates[D11_P + k];
          C[i * 3 + 2] += Slater_inv[i * D11_P + k] * Updates[2 * D11_P + k];
      }
  }

  // Compute B = 1 + VC : 3 x 3
  const double B0 = C[row1 * 3] + 1;
  const double B1 = C[row1 * 3 + 1];
  const double B2 = C[row1 * 3 + 2];
  const double B3 = C[row2 * 3];
  const double B4 = C[row2 * 3 + 1] + 1;
  const double B5 = C[row2 * 3 + 2];
  const double B6 = C[row3 * 3];
  const double B7 = C[row3 * 3 + 1];
  const double B8 = C[row3 * 3 + 2] + 1;

  // Check if determinant of B is not too close to zero
  double det;
  det = B0 * (B4 * B8 - B5 * B7) - B1 * (B3 * B8 - B5 * B6) +
        B2 * (B3 * B7 - B4 * B6);
  if (fabs(det) < breakdown) {
      return QMCKL_FAILURE;
  }

  // Update det(Slater) if passed
  if (determinant)
      *determinant *= det;

  // Compute B^{-1} with explicit formula for 3 x 3 inversion
  double __attribute__((aligned(8))) Binv[9], idet = 1.0 / det;
  Binv[0] = (B4 * B8 - B7 * B5) * idet;
  Binv[1] = -(B1 * B8 - B7 * B2) * idet;
  Binv[2] = (B1 * B5 - B4 * B2) * idet;
  Binv[3] = -(B3 * B8 - B6 * B5) * idet;
  Binv[4] = (B0 * B8 - B6 * B2) * idet;
  Binv[5] = -(B0 * B5 - B3 * B2) * idet;
  Binv[6] = (B3 * B7 - B6 * B4) * idet;
  Binv[7] = -(B0 * B7 - B6 * B1) * idet;
  Binv[8] = (B0 * B4 - B3 * B1) * idet;

  // tmp = B^{-1}D : 3 x D11_P
  double __attribute__((aligned(8))) tmp[3 * D11_P];
  double* r1dim = &(Slater_inv[row1 * D11_P]);
  double* r2dim = &(Slater_inv[row2 * D11_P]);
  double* r3dim = &(Slater_inv[row3 * D11_P]);
      IVDEP
      ALIGNED
  for (uint64_t j = 0; j < D11_P; j++) {
      tmp[j] = Binv[0] * r1dim[j] + Binv[1] * r2dim[j] + Binv[2] * r3dim[j];
      tmp[D11_P + j] =
          Binv[3] * r1dim[j] + Binv[4] * r2dim[j] + Binv[5] * r3dim[j];
      tmp[2 * D11_P + j] =
          Binv[6] * r1dim[j] + Binv[7] * r2dim[j] + Binv[8] * r3dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : 11 x D11_P
  for (uint64_t i = 0; i < 11; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D11_P; j++) {
          Slater_inv[i * D11_P + j] -= C[i * 3] * tmp[j];
          Slater_inv[i * D11_P + j] -= C[i * 3 + 1] * tmp[D11_P + j];
          Slater_inv[i * D11_P + j] -= C[i * 3 + 2] * tmp[2 * D11_P + j];
      }
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_woodbury_3x3_12(
    const qmckl_context context,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {
/*
    C := S^{-1} * U,    dim x 3
    B := 1 + V * C,     3 x 3
    D := V * S^{-1},    3 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
                          QMCKL_NULL_CONTEXT,
                          "qmckl_woodbury_3x3_12",
                          NULL);
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);
  const uint64_t row3 = (Updates_index[2] - 1);

  // Compute C = (S^T)^{-1}U : 12 x 3
  double __attribute__((aligned(8))) C[3 * 12];
  for (uint64_t i = 0; i < 12; i++) {
      C[i * 3] = 0;
      C[i * 3 + 1] = 0;
      C[i * 3 + 2] = 0;
      IVDEP
      ALIGNED
      for (uint64_t k = 0; k < D12_P; k++) {
          C[i * 3] += Slater_inv[i * D12_P + k] * Updates[k];
          C[i * 3 + 1] += Slater_inv[i * D12_P + k] * Updates[D12_P + k];
          C[i * 3 + 2] += Slater_inv[i * D12_P + k] * Updates[2 * D12_P + k];
      }
  }

  // Compute B = 1 + VC : 3 x 3
  const double B0 = C[row1 * 3] + 1;
  const double B1 = C[row1 * 3 + 1];
  const double B2 = C[row1 * 3 + 2];
  const double B3 = C[row2 * 3];
  const double B4 = C[row2 * 3 + 1] + 1;
  const double B5 = C[row2 * 3 + 2];
  const double B6 = C[row3 * 3];
  const double B7 = C[row3 * 3 + 1];
  const double B8 = C[row3 * 3 + 2] + 1;

  // Check if determinant of B is not too close to zero
  double det;
  det = B0 * (B4 * B8 - B5 * B7) - B1 * (B3 * B8 - B5 * B6) +
        B2 * (B3 * B7 - B4 * B6);
  if (fabs(det) < breakdown) {
      return QMCKL_FAILURE;
  }

  // Update det(Slater) if passed
  if (determinant)
      *determinant *= det;

  // Compute B^{-1} with explicit formula for 3 x 3 inversion
  double __attribute__((aligned(8))) Binv[9], idet = 1.0 / det;
  Binv[0] = (B4 * B8 - B7 * B5) * idet;
  Binv[1] = -(B1 * B8 - B7 * B2) * idet;
  Binv[2] = (B1 * B5 - B4 * B2) * idet;
  Binv[3] = -(B3 * B8 - B6 * B5) * idet;
  Binv[4] = (B0 * B8 - B6 * B2) * idet;
  Binv[5] = -(B0 * B5 - B3 * B2) * idet;
  Binv[6] = (B3 * B7 - B6 * B4) * idet;
  Binv[7] = -(B0 * B7 - B6 * B1) * idet;
  Binv[8] = (B0 * B4 - B3 * B1) * idet;

  // tmp = B^{-1}D : 3 x D12_P
  double __attribute__((aligned(8))) tmp[3 * D12_P];
  double* r1dim = &(Slater_inv[row1 * D12_P]);
  double* r2dim = &(Slater_inv[row2 * D12_P]);
  double* r3dim = &(Slater_inv[row3 * D12_P]);
      IVDEP
      ALIGNED
  for (uint64_t j = 0; j < D12_P; j++) {
      tmp[j] = Binv[0] * r1dim[j] + Binv[1] * r2dim[j] + Binv[2] * r3dim[j];
      tmp[D12_P + j] =
          Binv[3] * r1dim[j] + Binv[4] * r2dim[j] + Binv[5] * r3dim[j];
      tmp[2 * D12_P + j] =
          Binv[6] * r1dim[j] + Binv[7] * r2dim[j] + Binv[8] * r3dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : 12 x D12_P
  for (uint64_t i = 0; i < 12; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D12_P; j++) {
          Slater_inv[i * D12_P + j] -= C[i * 3] * tmp[j];
          Slater_inv[i * D12_P + j] -= C[i * 3 + 1] * tmp[D12_P + j];
          Slater_inv[i * D12_P + j] -= C[i * 3 + 2] * tmp[2 * D12_P + j];
      }
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_woodbury_3x3_13(
    const qmckl_context context,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {
/*
    C := S^{-1} * U,    dim x 3
    B := 1 + V * C,     3 x 3
    D := V * S^{-1},    3 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
                          QMCKL_NULL_CONTEXT,
                          "qmckl_woodbury_3x3_13",
                          NULL);
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);
  const uint64_t row3 = (Updates_index[2] - 1);

  // Compute C = (S^T)^{-1}U : 13 x 3
  double __attribute__((aligned(8))) C[3 * 13];
  for (uint64_t i = 0; i < 13; i++) {
      C[i * 3] = 0;
      C[i * 3 + 1] = 0;
      C[i * 3 + 2] = 0;
      IVDEP
      ALIGNED
      for (uint64_t k = 0; k < D13_P; k++) {
          C[i * 3] += Slater_inv[i * D13_P + k] * Updates[k];
          C[i * 3 + 1] += Slater_inv[i * D13_P + k] * Updates[D13_P + k];
          C[i * 3 + 2] += Slater_inv[i * D13_P + k] * Updates[2 * D13_P + k];
      }
  }

  // Compute B = 1 + VC : 3 x 3
  const double B0 = C[row1 * 3] + 1;
  const double B1 = C[row1 * 3 + 1];
  const double B2 = C[row1 * 3 + 2];
  const double B3 = C[row2 * 3];
  const double B4 = C[row2 * 3 + 1] + 1;
  const double B5 = C[row2 * 3 + 2];
  const double B6 = C[row3 * 3];
  const double B7 = C[row3 * 3 + 1];
  const double B8 = C[row3 * 3 + 2] + 1;

  // Check if determinant of B is not too close to zero
  double det;
  det = B0 * (B4 * B8 - B5 * B7) - B1 * (B3 * B8 - B5 * B6) +
        B2 * (B3 * B7 - B4 * B6);
  if (fabs(det) < breakdown) {
      return QMCKL_FAILURE;
  }

  // Update det(Slater) if passed
  if (determinant)
      *determinant *= det;

  // Compute B^{-1} with explicit formula for 3 x 3 inversion
  double __attribute__((aligned(8))) Binv[9], idet = 1.0 / det;
  Binv[0] = (B4 * B8 - B7 * B5) * idet;
  Binv[1] = -(B1 * B8 - B7 * B2) * idet;
  Binv[2] = (B1 * B5 - B4 * B2) * idet;
  Binv[3] = -(B3 * B8 - B6 * B5) * idet;
  Binv[4] = (B0 * B8 - B6 * B2) * idet;
  Binv[5] = -(B0 * B5 - B3 * B2) * idet;
  Binv[6] = (B3 * B7 - B6 * B4) * idet;
  Binv[7] = -(B0 * B7 - B6 * B1) * idet;
  Binv[8] = (B0 * B4 - B3 * B1) * idet;

  // tmp = B^{-1}D : 3 x D13_P
  double __attribute__((aligned(8))) tmp[3 * D13_P];
  double* r1dim = &(Slater_inv[row1 * D13_P]);
  double* r2dim = &(Slater_inv[row2 * D13_P]);
  double* r3dim = &(Slater_inv[row3 * D13_P]);
      IVDEP
      ALIGNED
  for (uint64_t j = 0; j < D13_P; j++) {
      tmp[j] = Binv[0] * r1dim[j] + Binv[1] * r2dim[j] + Binv[2] * r3dim[j];
      tmp[D13_P + j] =
          Binv[3] * r1dim[j] + Binv[4] * r2dim[j] + Binv[5] * r3dim[j];
      tmp[2 * D13_P + j] =
          Binv[6] * r1dim[j] + Binv[7] * r2dim[j] + Binv[8] * r3dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : 13 x D13_P
  for (uint64_t i = 0; i < 13; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D13_P; j++) {
          Slater_inv[i * D13_P + j] -= C[i * 3] * tmp[j];
          Slater_inv[i * D13_P + j] -= C[i * 3 + 1] * tmp[D13_P + j];
          Slater_inv[i * D13_P + j] -= C[i * 3 + 2] * tmp[2 * D13_P + j];
      }
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_woodbury_3x3_14(
    const qmckl_context context,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {
/*
    C := S^{-1} * U,    dim x 3
    B := 1 + V * C,     3 x 3
    D := V * S^{-1},    3 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
                          QMCKL_NULL_CONTEXT,
                          "qmckl_woodbury_3x3_14",
                          NULL);
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);
  const uint64_t row3 = (Updates_index[2] - 1);

  // Compute C = (S^T)^{-1}U : 14 x 3
  double __attribute__((aligned(8))) C[3 * 14];
  for (uint64_t i = 0; i < 14; i++) {
      C[i * 3] = 0;
      C[i * 3 + 1] = 0;
      C[i * 3 + 2] = 0;
      IVDEP
      ALIGNED
      for (uint64_t k = 0; k < D14_P; k++) {
          C[i * 3] += Slater_inv[i * D14_P + k] * Updates[k];
          C[i * 3 + 1] += Slater_inv[i * D14_P + k] * Updates[D14_P + k];
          C[i * 3 + 2] += Slater_inv[i * D14_P + k] * Updates[2 * D14_P + k];
      }
  }

  // Compute B = 1 + VC : 3 x 3
  const double B0 = C[row1 * 3] + 1;
  const double B1 = C[row1 * 3 + 1];
  const double B2 = C[row1 * 3 + 2];
  const double B3 = C[row2 * 3];
  const double B4 = C[row2 * 3 + 1] + 1;
  const double B5 = C[row2 * 3 + 2];
  const double B6 = C[row3 * 3];
  const double B7 = C[row3 * 3 + 1];
  const double B8 = C[row3 * 3 + 2] + 1;

  // Check if determinant of B is not too close to zero
  double det;
  det = B0 * (B4 * B8 - B5 * B7) - B1 * (B3 * B8 - B5 * B6) +
        B2 * (B3 * B7 - B4 * B6);
  if (fabs(det) < breakdown) {
      return QMCKL_FAILURE;
  }

  // Update det(Slater) if passed
  if (determinant)
      *determinant *= det;

  // Compute B^{-1} with explicit formula for 3 x 3 inversion
  double __attribute__((aligned(8))) Binv[9], idet = 1.0 / det;
  Binv[0] = (B4 * B8 - B7 * B5) * idet;
  Binv[1] = -(B1 * B8 - B7 * B2) * idet;
  Binv[2] = (B1 * B5 - B4 * B2) * idet;
  Binv[3] = -(B3 * B8 - B6 * B5) * idet;
  Binv[4] = (B0 * B8 - B6 * B2) * idet;
  Binv[5] = -(B0 * B5 - B3 * B2) * idet;
  Binv[6] = (B3 * B7 - B6 * B4) * idet;
  Binv[7] = -(B0 * B7 - B6 * B1) * idet;
  Binv[8] = (B0 * B4 - B3 * B1) * idet;

  // tmp = B^{-1}D : 3 x D14_P
  double __attribute__((aligned(8))) tmp[3 * D14_P];
  double* r1dim = &(Slater_inv[row1 * D14_P]);
  double* r2dim = &(Slater_inv[row2 * D14_P]);
  double* r3dim = &(Slater_inv[row3 * D14_P]);
      IVDEP
      ALIGNED
  for (uint64_t j = 0; j < D14_P; j++) {
      tmp[j] = Binv[0] * r1dim[j] + Binv[1] * r2dim[j] + Binv[2] * r3dim[j];
      tmp[D14_P + j] =
          Binv[3] * r1dim[j] + Binv[4] * r2dim[j] + Binv[5] * r3dim[j];
      tmp[2 * D14_P + j] =
          Binv[6] * r1dim[j] + Binv[7] * r2dim[j] + Binv[8] * r3dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : 14 x D14_P
  for (uint64_t i = 0; i < 14; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D14_P; j++) {
          Slater_inv[i * D14_P + j] -= C[i * 3] * tmp[j];
          Slater_inv[i * D14_P + j] -= C[i * 3 + 1] * tmp[D14_P + j];
          Slater_inv[i * D14_P + j] -= C[i * 3 + 2] * tmp[2 * D14_P + j];
      }
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_woodbury_3x3_15(
    const qmckl_context context,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {
/*
    C := S^{-1} * U,    dim x 3
    B := 1 + V * C,     3 x 3
    D := V * S^{-1},    3 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
                          QMCKL_NULL_CONTEXT,
                          "qmckl_woodbury_3x3_15",
                          NULL);
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);
  const uint64_t row3 = (Updates_index[2] - 1);

  // Compute C = (S^T)^{-1}U : 15 x 3
  double __attribute__((aligned(8))) C[3 * 15];
  for (uint64_t i = 0; i < 15; i++) {
      C[i * 3] = 0;
      C[i * 3 + 1] = 0;
      C[i * 3 + 2] = 0;
      IVDEP
      ALIGNED
      for (uint64_t k = 0; k < D15_P; k++) {
          C[i * 3] += Slater_inv[i * D15_P + k] * Updates[k];
          C[i * 3 + 1] += Slater_inv[i * D15_P + k] * Updates[D15_P + k];
          C[i * 3 + 2] += Slater_inv[i * D15_P + k] * Updates[2 * D15_P + k];
      }
  }

  // Compute B = 1 + VC : 3 x 3
  const double B0 = C[row1 * 3] + 1;
  const double B1 = C[row1 * 3 + 1];
  const double B2 = C[row1 * 3 + 2];
  const double B3 = C[row2 * 3];
  const double B4 = C[row2 * 3 + 1] + 1;
  const double B5 = C[row2 * 3 + 2];
  const double B6 = C[row3 * 3];
  const double B7 = C[row3 * 3 + 1];
  const double B8 = C[row3 * 3 + 2] + 1;

  // Check if determinant of B is not too close to zero
  double det;
  det = B0 * (B4 * B8 - B5 * B7) - B1 * (B3 * B8 - B5 * B6) +
        B2 * (B3 * B7 - B4 * B6);
  if (fabs(det) < breakdown) {
      return QMCKL_FAILURE;
  }

  // Update det(Slater) if passed
  if (determinant)
      *determinant *= det;

  // Compute B^{-1} with explicit formula for 3 x 3 inversion
  double __attribute__((aligned(8))) Binv[9], idet = 1.0 / det;
  Binv[0] = (B4 * B8 - B7 * B5) * idet;
  Binv[1] = -(B1 * B8 - B7 * B2) * idet;
  Binv[2] = (B1 * B5 - B4 * B2) * idet;
  Binv[3] = -(B3 * B8 - B6 * B5) * idet;
  Binv[4] = (B0 * B8 - B6 * B2) * idet;
  Binv[5] = -(B0 * B5 - B3 * B2) * idet;
  Binv[6] = (B3 * B7 - B6 * B4) * idet;
  Binv[7] = -(B0 * B7 - B6 * B1) * idet;
  Binv[8] = (B0 * B4 - B3 * B1) * idet;

  // tmp = B^{-1}D : 3 x D15_P
  double __attribute__((aligned(8))) tmp[3 * D15_P];
  double* r1dim = &(Slater_inv[row1 * D15_P]);
  double* r2dim = &(Slater_inv[row2 * D15_P]);
  double* r3dim = &(Slater_inv[row3 * D15_P]);
      IVDEP
      ALIGNED
  for (uint64_t j = 0; j < D15_P; j++) {
      tmp[j] = Binv[0] * r1dim[j] + Binv[1] * r2dim[j] + Binv[2] * r3dim[j];
      tmp[D15_P + j] =
          Binv[3] * r1dim[j] + Binv[4] * r2dim[j] + Binv[5] * r3dim[j];
      tmp[2 * D15_P + j] =
          Binv[6] * r1dim[j] + Binv[7] * r2dim[j] + Binv[8] * r3dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : 15 x D15_P
  for (uint64_t i = 0; i < 15; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D15_P; j++) {
          Slater_inv[i * D15_P + j] -= C[i * 3] * tmp[j];
          Slater_inv[i * D15_P + j] -= C[i * 3 + 1] * tmp[D15_P + j];
          Slater_inv[i * D15_P + j] -= C[i * 3 + 2] * tmp[2 * D15_P + j];
      }
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_woodbury_3x3_16(
    const qmckl_context context,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {
/*
    C := S^{-1} * U,    dim x 3
    B := 1 + V * C,     3 x 3
    D := V * S^{-1},    3 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
                          QMCKL_NULL_CONTEXT,
                          "qmckl_woodbury_3x3_16",
                          NULL);
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);
  const uint64_t row3 = (Updates_index[2] - 1);

  // Compute C = (S^T)^{-1}U : 16 x 3
  double __attribute__((aligned(8))) C[3 * 16];
  for (uint64_t i = 0; i < 16; i++) {
      C[i * 3] = 0;
      C[i * 3 + 1] = 0;
      C[i * 3 + 2] = 0;
      IVDEP
      ALIGNED
      for (uint64_t k = 0; k < D16_P; k++) {
          C[i * 3] += Slater_inv[i * D16_P + k] * Updates[k];
          C[i * 3 + 1] += Slater_inv[i * D16_P + k] * Updates[D16_P + k];
          C[i * 3 + 2] += Slater_inv[i * D16_P + k] * Updates[2 * D16_P + k];
      }
  }

  // Compute B = 1 + VC : 3 x 3
  const double B0 = C[row1 * 3] + 1;
  const double B1 = C[row1 * 3 + 1];
  const double B2 = C[row1 * 3 + 2];
  const double B3 = C[row2 * 3];
  const double B4 = C[row2 * 3 + 1] + 1;
  const double B5 = C[row2 * 3 + 2];
  const double B6 = C[row3 * 3];
  const double B7 = C[row3 * 3 + 1];
  const double B8 = C[row3 * 3 + 2] + 1;

  // Check if determinant of B is not too close to zero
  double det;
  det = B0 * (B4 * B8 - B5 * B7) - B1 * (B3 * B8 - B5 * B6) +
        B2 * (B3 * B7 - B4 * B6);
  if (fabs(det) < breakdown) {
      return QMCKL_FAILURE;
  }

  // Update det(Slater) if passed
  if (determinant)
      *determinant *= det;

  // Compute B^{-1} with explicit formula for 3 x 3 inversion
  double __attribute__((aligned(8))) Binv[9], idet = 1.0 / det;
  Binv[0] = (B4 * B8 - B7 * B5) * idet;
  Binv[1] = -(B1 * B8 - B7 * B2) * idet;
  Binv[2] = (B1 * B5 - B4 * B2) * idet;
  Binv[3] = -(B3 * B8 - B6 * B5) * idet;
  Binv[4] = (B0 * B8 - B6 * B2) * idet;
  Binv[5] = -(B0 * B5 - B3 * B2) * idet;
  Binv[6] = (B3 * B7 - B6 * B4) * idet;
  Binv[7] = -(B0 * B7 - B6 * B1) * idet;
  Binv[8] = (B0 * B4 - B3 * B1) * idet;

  // tmp = B^{-1}D : 3 x D16_P
  double __attribute__((aligned(8))) tmp[3 * D16_P];
  double* r1dim = &(Slater_inv[row1 * D16_P]);
  double* r2dim = &(Slater_inv[row2 * D16_P]);
  double* r3dim = &(Slater_inv[row3 * D16_P]);
      IVDEP
      ALIGNED
  for (uint64_t j = 0; j < D16_P; j++) {
      tmp[j] = Binv[0] * r1dim[j] + Binv[1] * r2dim[j] + Binv[2] * r3dim[j];
      tmp[D16_P + j] =
          Binv[3] * r1dim[j] + Binv[4] * r2dim[j] + Binv[5] * r3dim[j];
      tmp[2 * D16_P + j] =
          Binv[6] * r1dim[j] + Binv[7] * r2dim[j] + Binv[8] * r3dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : 16 x D16_P
  for (uint64_t i = 0; i < 16; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D16_P; j++) {
          Slater_inv[i * D16_P + j] -= C[i * 3] * tmp[j];
          Slater_inv[i * D16_P + j] -= C[i * 3 + 1] * tmp[D16_P + j];
          Slater_inv[i * D16_P + j] -= C[i * 3 + 2] * tmp[2 * D16_P + j];
      }
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_woodbury_3x3_17(
    const qmckl_context context,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {
/*
    C := S^{-1} * U,    dim x 3
    B := 1 + V * C,     3 x 3
    D := V * S^{-1},    3 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
                          QMCKL_NULL_CONTEXT,
                          "qmckl_woodbury_3x3_17",
                          NULL);
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);
  const uint64_t row3 = (Updates_index[2] - 1);

  // Compute C = (S^T)^{-1}U : 17 x 3
  double __attribute__((aligned(8))) C[3 * 17];
  for (uint64_t i = 0; i < 17; i++) {
      C[i * 3] = 0;
      C[i * 3 + 1] = 0;
      C[i * 3 + 2] = 0;
      IVDEP
      ALIGNED
      for (uint64_t k = 0; k < D17_P; k++) {
          C[i * 3] += Slater_inv[i * D17_P + k] * Updates[k];
          C[i * 3 + 1] += Slater_inv[i * D17_P + k] * Updates[D17_P + k];
          C[i * 3 + 2] += Slater_inv[i * D17_P + k] * Updates[2 * D17_P + k];
      }
  }

  // Compute B = 1 + VC : 3 x 3
  const double B0 = C[row1 * 3] + 1;
  const double B1 = C[row1 * 3 + 1];
  const double B2 = C[row1 * 3 + 2];
  const double B3 = C[row2 * 3];
  const double B4 = C[row2 * 3 + 1] + 1;
  const double B5 = C[row2 * 3 + 2];
  const double B6 = C[row3 * 3];
  const double B7 = C[row3 * 3 + 1];
  const double B8 = C[row3 * 3 + 2] + 1;

  // Check if determinant of B is not too close to zero
  double det;
  det = B0 * (B4 * B8 - B5 * B7) - B1 * (B3 * B8 - B5 * B6) +
        B2 * (B3 * B7 - B4 * B6);
  if (fabs(det) < breakdown) {
      return QMCKL_FAILURE;
  }

  // Update det(Slater) if passed
  if (determinant)
      *determinant *= det;

  // Compute B^{-1} with explicit formula for 3 x 3 inversion
  double __attribute__((aligned(8))) Binv[9], idet = 1.0 / det;
  Binv[0] = (B4 * B8 - B7 * B5) * idet;
  Binv[1] = -(B1 * B8 - B7 * B2) * idet;
  Binv[2] = (B1 * B5 - B4 * B2) * idet;
  Binv[3] = -(B3 * B8 - B6 * B5) * idet;
  Binv[4] = (B0 * B8 - B6 * B2) * idet;
  Binv[5] = -(B0 * B5 - B3 * B2) * idet;
  Binv[6] = (B3 * B7 - B6 * B4) * idet;
  Binv[7] = -(B0 * B7 - B6 * B1) * idet;
  Binv[8] = (B0 * B4 - B3 * B1) * idet;

  // tmp = B^{-1}D : 3 x D17_P
  double __attribute__((aligned(8))) tmp[3 * D17_P];
  double* r1dim = &(Slater_inv[row1 * D17_P]);
  double* r2dim = &(Slater_inv[row2 * D17_P]);
  double* r3dim = &(Slater_inv[row3 * D17_P]);
      IVDEP
      ALIGNED
  for (uint64_t j = 0; j < D17_P; j++) {
      tmp[j] = Binv[0] * r1dim[j] + Binv[1] * r2dim[j] + Binv[2] * r3dim[j];
      tmp[D17_P + j] =
          Binv[3] * r1dim[j] + Binv[4] * r2dim[j] + Binv[5] * r3dim[j];
      tmp[2 * D17_P + j] =
          Binv[6] * r1dim[j] + Binv[7] * r2dim[j] + Binv[8] * r3dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : 17 x D17_P
  for (uint64_t i = 0; i < 17; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D17_P; j++) {
          Slater_inv[i * D17_P + j] -= C[i * 3] * tmp[j];
          Slater_inv[i * D17_P + j] -= C[i * 3 + 1] * tmp[D17_P + j];
          Slater_inv[i * D17_P + j] -= C[i * 3 + 2] * tmp[2 * D17_P + j];
      }
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_woodbury_3x3_18(
    const qmckl_context context,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {
/*
    C := S^{-1} * U,    dim x 3
    B := 1 + V * C,     3 x 3
    D := V * S^{-1},    3 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
                          QMCKL_NULL_CONTEXT,
                          "qmckl_woodbury_3x3_18",
                          NULL);
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);
  const uint64_t row3 = (Updates_index[2] - 1);

  // Compute C = (S^T)^{-1}U : 18 x 3
  double __attribute__((aligned(8))) C[3 * 18];
  for (uint64_t i = 0; i < 18; i++) {
      C[i * 3] = 0;
      C[i * 3 + 1] = 0;
      C[i * 3 + 2] = 0;
      IVDEP
      ALIGNED
      for (uint64_t k = 0; k < D18_P; k++) {
          C[i * 3] += Slater_inv[i * D18_P + k] * Updates[k];
          C[i * 3 + 1] += Slater_inv[i * D18_P + k] * Updates[D18_P + k];
          C[i * 3 + 2] += Slater_inv[i * D18_P + k] * Updates[2 * D18_P + k];
      }
  }

  // Compute B = 1 + VC : 3 x 3
  const double B0 = C[row1 * 3] + 1;
  const double B1 = C[row1 * 3 + 1];
  const double B2 = C[row1 * 3 + 2];
  const double B3 = C[row2 * 3];
  const double B4 = C[row2 * 3 + 1] + 1;
  const double B5 = C[row2 * 3 + 2];
  const double B6 = C[row3 * 3];
  const double B7 = C[row3 * 3 + 1];
  const double B8 = C[row3 * 3 + 2] + 1;

  // Check if determinant of B is not too close to zero
  double det;
  det = B0 * (B4 * B8 - B5 * B7) - B1 * (B3 * B8 - B5 * B6) +
        B2 * (B3 * B7 - B4 * B6);
  if (fabs(det) < breakdown) {
      return QMCKL_FAILURE;
  }

  // Update det(Slater) if passed
  if (determinant)
      *determinant *= det;

  // Compute B^{-1} with explicit formula for 3 x 3 inversion
  double __attribute__((aligned(8))) Binv[9], idet = 1.0 / det;
  Binv[0] = (B4 * B8 - B7 * B5) * idet;
  Binv[1] = -(B1 * B8 - B7 * B2) * idet;
  Binv[2] = (B1 * B5 - B4 * B2) * idet;
  Binv[3] = -(B3 * B8 - B6 * B5) * idet;
  Binv[4] = (B0 * B8 - B6 * B2) * idet;
  Binv[5] = -(B0 * B5 - B3 * B2) * idet;
  Binv[6] = (B3 * B7 - B6 * B4) * idet;
  Binv[7] = -(B0 * B7 - B6 * B1) * idet;
  Binv[8] = (B0 * B4 - B3 * B1) * idet;

  // tmp = B^{-1}D : 3 x D18_P
  double __attribute__((aligned(8))) tmp[3 * D18_P];
  double* r1dim = &(Slater_inv[row1 * D18_P]);
  double* r2dim = &(Slater_inv[row2 * D18_P]);
  double* r3dim = &(Slater_inv[row3 * D18_P]);
      IVDEP
      ALIGNED
  for (uint64_t j = 0; j < D18_P; j++) {
      tmp[j] = Binv[0] * r1dim[j] + Binv[1] * r2dim[j] + Binv[2] * r3dim[j];
      tmp[D18_P + j] =
          Binv[3] * r1dim[j] + Binv[4] * r2dim[j] + Binv[5] * r3dim[j];
      tmp[2 * D18_P + j] =
          Binv[6] * r1dim[j] + Binv[7] * r2dim[j] + Binv[8] * r3dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : 18 x D18_P
  for (uint64_t i = 0; i < 18; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D18_P; j++) {
          Slater_inv[i * D18_P + j] -= C[i * 3] * tmp[j];
          Slater_inv[i * D18_P + j] -= C[i * 3 + 1] * tmp[D18_P + j];
          Slater_inv[i * D18_P + j] -= C[i * 3 + 2] * tmp[2 * D18_P + j];
      }
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_woodbury_3x3_19(
    const qmckl_context context,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {
/*
    C := S^{-1} * U,    dim x 3
    B := 1 + V * C,     3 x 3
    D := V * S^{-1},    3 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
                          QMCKL_NULL_CONTEXT,
                          "qmckl_woodbury_3x3_19",
                          NULL);
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);
  const uint64_t row3 = (Updates_index[2] - 1);

  // Compute C = (S^T)^{-1}U : 19 x 3
  double __attribute__((aligned(8))) C[3 * 19];
  for (uint64_t i = 0; i < 19; i++) {
      C[i * 3] = 0;
      C[i * 3 + 1] = 0;
      C[i * 3 + 2] = 0;
      IVDEP
      ALIGNED
      for (uint64_t k = 0; k < D19_P; k++) {
          C[i * 3] += Slater_inv[i * D19_P + k] * Updates[k];
          C[i * 3 + 1] += Slater_inv[i * D19_P + k] * Updates[D19_P + k];
          C[i * 3 + 2] += Slater_inv[i * D19_P + k] * Updates[2 * D19_P + k];
      }
  }

  // Compute B = 1 + VC : 3 x 3
  const double B0 = C[row1 * 3] + 1;
  const double B1 = C[row1 * 3 + 1];
  const double B2 = C[row1 * 3 + 2];
  const double B3 = C[row2 * 3];
  const double B4 = C[row2 * 3 + 1] + 1;
  const double B5 = C[row2 * 3 + 2];
  const double B6 = C[row3 * 3];
  const double B7 = C[row3 * 3 + 1];
  const double B8 = C[row3 * 3 + 2] + 1;

  // Check if determinant of B is not too close to zero
  double det;
  det = B0 * (B4 * B8 - B5 * B7) - B1 * (B3 * B8 - B5 * B6) +
        B2 * (B3 * B7 - B4 * B6);
  if (fabs(det) < breakdown) {
      return QMCKL_FAILURE;
  }

  // Update det(Slater) if passed
  if (determinant)
      *determinant *= det;

  // Compute B^{-1} with explicit formula for 3 x 3 inversion
  double __attribute__((aligned(8))) Binv[9], idet = 1.0 / det;
  Binv[0] = (B4 * B8 - B7 * B5) * idet;
  Binv[1] = -(B1 * B8 - B7 * B2) * idet;
  Binv[2] = (B1 * B5 - B4 * B2) * idet;
  Binv[3] = -(B3 * B8 - B6 * B5) * idet;
  Binv[4] = (B0 * B8 - B6 * B2) * idet;
  Binv[5] = -(B0 * B5 - B3 * B2) * idet;
  Binv[6] = (B3 * B7 - B6 * B4) * idet;
  Binv[7] = -(B0 * B7 - B6 * B1) * idet;
  Binv[8] = (B0 * B4 - B3 * B1) * idet;

  // tmp = B^{-1}D : 3 x D19_P
  double __attribute__((aligned(8))) tmp[3 * D19_P];
  double* r1dim = &(Slater_inv[row1 * D19_P]);
  double* r2dim = &(Slater_inv[row2 * D19_P]);
  double* r3dim = &(Slater_inv[row3 * D19_P]);
      IVDEP
      ALIGNED
  for (uint64_t j = 0; j < D19_P; j++) {
      tmp[j] = Binv[0] * r1dim[j] + Binv[1] * r2dim[j] + Binv[2] * r3dim[j];
      tmp[D19_P + j] =
          Binv[3] * r1dim[j] + Binv[4] * r2dim[j] + Binv[5] * r3dim[j];
      tmp[2 * D19_P + j] =
          Binv[6] * r1dim[j] + Binv[7] * r2dim[j] + Binv[8] * r3dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : 19 x D19_P
  for (uint64_t i = 0; i < 19; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D19_P; j++) {
          Slater_inv[i * D19_P + j] -= C[i * 3] * tmp[j];
          Slater_inv[i * D19_P + j] -= C[i * 3 + 1] * tmp[D19_P + j];
          Slater_inv[i * D19_P + j] -= C[i * 3 + 2] * tmp[2 * D19_P + j];
      }
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_woodbury_3x3_20(
    const qmckl_context context,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {
/*
    C := S^{-1} * U,    dim x 3
    B := 1 + V * C,     3 x 3
    D := V * S^{-1},    3 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
                          QMCKL_NULL_CONTEXT,
                          "qmckl_woodbury_3x3_20",
                          NULL);
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);
  const uint64_t row3 = (Updates_index[2] - 1);

  // Compute C = (S^T)^{-1}U : 20 x 3
  double __attribute__((aligned(8))) C[3 * 20];
  for (uint64_t i = 0; i < 20; i++) {
      C[i * 3] = 0;
      C[i * 3 + 1] = 0;
      C[i * 3 + 2] = 0;
      IVDEP
      ALIGNED
      for (uint64_t k = 0; k < D20_P; k++) {
          C[i * 3] += Slater_inv[i * D20_P + k] * Updates[k];
          C[i * 3 + 1] += Slater_inv[i * D20_P + k] * Updates[D20_P + k];
          C[i * 3 + 2] += Slater_inv[i * D20_P + k] * Updates[2 * D20_P + k];
      }
  }

  // Compute B = 1 + VC : 3 x 3
  const double B0 = C[row1 * 3] + 1;
  const double B1 = C[row1 * 3 + 1];
  const double B2 = C[row1 * 3 + 2];
  const double B3 = C[row2 * 3];
  const double B4 = C[row2 * 3 + 1] + 1;
  const double B5 = C[row2 * 3 + 2];
  const double B6 = C[row3 * 3];
  const double B7 = C[row3 * 3 + 1];
  const double B8 = C[row3 * 3 + 2] + 1;

  // Check if determinant of B is not too close to zero
  double det;
  det = B0 * (B4 * B8 - B5 * B7) - B1 * (B3 * B8 - B5 * B6) +
        B2 * (B3 * B7 - B4 * B6);
  if (fabs(det) < breakdown) {
      return QMCKL_FAILURE;
  }

  // Update det(Slater) if passed
  if (determinant)
      *determinant *= det;

  // Compute B^{-1} with explicit formula for 3 x 3 inversion
  double __attribute__((aligned(8))) Binv[9], idet = 1.0 / det;
  Binv[0] = (B4 * B8 - B7 * B5) * idet;
  Binv[1] = -(B1 * B8 - B7 * B2) * idet;
  Binv[2] = (B1 * B5 - B4 * B2) * idet;
  Binv[3] = -(B3 * B8 - B6 * B5) * idet;
  Binv[4] = (B0 * B8 - B6 * B2) * idet;
  Binv[5] = -(B0 * B5 - B3 * B2) * idet;
  Binv[6] = (B3 * B7 - B6 * B4) * idet;
  Binv[7] = -(B0 * B7 - B6 * B1) * idet;
  Binv[8] = (B0 * B4 - B3 * B1) * idet;

  // tmp = B^{-1}D : 3 x D20_P
  double __attribute__((aligned(8))) tmp[3 * D20_P];
  double* r1dim = &(Slater_inv[row1 * D20_P]);
  double* r2dim = &(Slater_inv[row2 * D20_P]);
  double* r3dim = &(Slater_inv[row3 * D20_P]);
      IVDEP
      ALIGNED
  for (uint64_t j = 0; j < D20_P; j++) {
      tmp[j] = Binv[0] * r1dim[j] + Binv[1] * r2dim[j] + Binv[2] * r3dim[j];
      tmp[D20_P + j] =
          Binv[3] * r1dim[j] + Binv[4] * r2dim[j] + Binv[5] * r3dim[j];
      tmp[2 * D20_P + j] =
          Binv[6] * r1dim[j] + Binv[7] * r2dim[j] + Binv[8] * r3dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : 20 x D20_P
  for (uint64_t i = 0; i < 20; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D20_P; j++) {
          Slater_inv[i * D20_P + j] -= C[i * 3] * tmp[j];
          Slater_inv[i * D20_P + j] -= C[i * 3 + 1] * tmp[D20_P + j];
          Slater_inv[i * D20_P + j] -= C[i * 3 + 2] * tmp[2 * D20_P + j];
      }
  }

  return QMCKL_SUCCESS;
}

static inline qmckl_exit_code qmckl_woodbury_3x3_21(
    const qmckl_context context,
    const double* restrict Updates,
    const uint64_t* restrict Updates_index,
    const double breakdown,
    double* restrict Slater_inv,
    double* restrict determinant) {
/*
    C := S^{-1} * U,    dim x 3
    B := 1 + V * C,     3 x 3
    D := V * S^{-1},    3 x dim
*/

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(context,
                          QMCKL_NULL_CONTEXT,
                          "qmckl_woodbury_3x3_21",
                          NULL);
  }

  const uint64_t row1 = (Updates_index[0] - 1);
  const uint64_t row2 = (Updates_index[1] - 1);
  const uint64_t row3 = (Updates_index[2] - 1);

  // Compute C = (S^T)^{-1}U : 21 x 3
  double __attribute__((aligned(8))) C[3 * 21];
  for (uint64_t i = 0; i < 21; i++) {
      C[i * 3] = 0;
      C[i * 3 + 1] = 0;
      C[i * 3 + 2] = 0;
      IVDEP
      ALIGNED
      for (uint64_t k = 0; k < D21_P; k++) {
          C[i * 3] += Slater_inv[i * D21_P + k] * Updates[k];
          C[i * 3 + 1] += Slater_inv[i * D21_P + k] * Updates[D21_P + k];
          C[i * 3 + 2] += Slater_inv[i * D21_P + k] * Updates[2 * D21_P + k];
      }
  }

  // Compute B = 1 + VC : 3 x 3
  const double B0 = C[row1 * 3] + 1;
  const double B1 = C[row1 * 3 + 1];
  const double B2 = C[row1 * 3 + 2];
  const double B3 = C[row2 * 3];
  const double B4 = C[row2 * 3 + 1] + 1;
  const double B5 = C[row2 * 3 + 2];
  const double B6 = C[row3 * 3];
  const double B7 = C[row3 * 3 + 1];
  const double B8 = C[row3 * 3 + 2] + 1;

  // Check if determinant of B is not too close to zero
  double det;
  det = B0 * (B4 * B8 - B5 * B7) - B1 * (B3 * B8 - B5 * B6) +
        B2 * (B3 * B7 - B4 * B6);
  if (fabs(det) < breakdown) {
      return QMCKL_FAILURE;
  }

  // Update det(Slater) if passed
  if (determinant)
      *determinant *= det;

  // Compute B^{-1} with explicit formula for 3 x 3 inversion
  double __attribute__((aligned(8))) Binv[9], idet = 1.0 / det;
  Binv[0] = (B4 * B8 - B7 * B5) * idet;
  Binv[1] = -(B1 * B8 - B7 * B2) * idet;
  Binv[2] = (B1 * B5 - B4 * B2) * idet;
  Binv[3] = -(B3 * B8 - B6 * B5) * idet;
  Binv[4] = (B0 * B8 - B6 * B2) * idet;
  Binv[5] = -(B0 * B5 - B3 * B2) * idet;
  Binv[6] = (B3 * B7 - B6 * B4) * idet;
  Binv[7] = -(B0 * B7 - B6 * B1) * idet;
  Binv[8] = (B0 * B4 - B3 * B1) * idet;

  // tmp = B^{-1}D : 3 x D21_P
  double __attribute__((aligned(8))) tmp[3 * D21_P];
  double* r1dim = &(Slater_inv[row1 * D21_P]);
  double* r2dim = &(Slater_inv[row2 * D21_P]);
  double* r3dim = &(Slater_inv[row3 * D21_P]);
      IVDEP
      ALIGNED
  for (uint64_t j = 0; j < D21_P; j++) {
      tmp[j] = Binv[0] * r1dim[j] + Binv[1] * r2dim[j] + Binv[2] * r3dim[j];
      tmp[D21_P + j] =
          Binv[3] * r1dim[j] + Binv[4] * r2dim[j] + Binv[5] * r3dim[j];
      tmp[2 * D21_P + j] =
          Binv[6] * r1dim[j] + Binv[7] * r2dim[j] + Binv[8] * r3dim[j];
  }

  // Compute (S^T)^{-1} - C * tmp : 21 x D21_P
  for (uint64_t i = 0; i < 21; i++) {
      IVDEP
      ALIGNED
      for (uint64_t j = 0; j < D21_P; j++) {
          Slater_inv[i * D21_P + j] -= C[i * 3] * tmp[j];
          Slater_inv[i * D21_P + j] -= C[i * 3 + 1] * tmp[D21_P + j];
          Slater_inv[i * D21_P + j] -= C[i * 3 + 2] * tmp[2 * D21_P + j];
      }
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_woodbury_3x3(const qmckl_context context,
    const uint64_t LDS,
    const uint64_t Dim,
    const double* Updates,
    const uint64_t* Updates_index,
    const double breakdown,
    double* Slater_inv,
    double* determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(
        context,
        QMCKL_NULL_CONTEXT,
        "qmckl_woodbury_3x3",
        NULL);
  }

  #ifdef HAVE_HPC__BROKEN_WITH_CRAY
  if (LDS == (1+(Dim-1)/SIMD_LENGTH)*SIMD_LENGTH) { // Most cases
    switch (Dim) {
      case 2:  
        return qmckl_woodbury_3x3_2(context,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 3:  
        return qmckl_woodbury_3x3_3(context,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 4:  
        return qmckl_woodbury_3x3_4(context,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 5:  
        return qmckl_woodbury_3x3_5(context,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 6:  
        return qmckl_woodbury_3x3_6(context,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 7:  
        return qmckl_woodbury_3x3_7(context,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 8:  
        return qmckl_woodbury_3x3_8(context,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 9:  
        return qmckl_woodbury_3x3_9(context,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 10:  
        return qmckl_woodbury_3x3_10(context,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 11:  
        return qmckl_woodbury_3x3_11(context,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 12:  
        return qmckl_woodbury_3x3_12(context,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 13:  
        return qmckl_woodbury_3x3_13(context,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 14:  
        return qmckl_woodbury_3x3_14(context,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 15:  
        return qmckl_woodbury_3x3_15(context,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 16:  
        return qmckl_woodbury_3x3_16(context,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 17:  
        return qmckl_woodbury_3x3_17(context,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 18:  
        return qmckl_woodbury_3x3_18(context,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 19:  
        return qmckl_woodbury_3x3_19(context,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 20:  
        return qmckl_woodbury_3x3_20(context,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
      case 21:  
        return qmckl_woodbury_3x3_21(context,
          Updates,
          Updates_index,
          breakdown,
          Slater_inv,
          determinant);
    }
  }
  else { // Updating smaller sub-matrix
    return qmckl_woodbury_3x3_hpc(
        context,
        LDS,
        Dim,
        Updates,
        Updates_index,
        breakdown,
        Slater_inv,
        determinant);
  }
  #else
  return qmckl_woodbury_3x3_doc(
      context,
      LDS,
      Dim,
      Updates,
      Updates_index,
      breakdown,
      Slater_inv,
      determinant);
  // return qmckl_woodbury_3x3_hpc(
  //     context,
  //     LDS,
  //     Dim,
  //     Updates,
  //     Updates_index,
  //     breakdown,
  //     Slater_inv,
  //     determinant);
  #endif

  return QMCKL_FAILURE;
}

qmckl_exit_code qmckl_sm_splitting_hpc(
    const qmckl_context context,
    const uint64_t LDS,
    const uint64_t Dim,
    const uint64_t N_updates,
    const double* Updates,
    const uint64_t* Updates_index,
    const double breakdown,
    double* Slater_inv,
    double* determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(
        context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_splitting_hpc",
        NULL);
  }

  double __attribute__((aligned(8))) later_updates[LDS * N_updates];
  uint64_t later_index[N_updates];
  uint64_t later = 0;

  qmckl_exit_code rc;
  if (LDS == (1+(Dim-1)/SIMD_LENGTH)*SIMD_LENGTH) {
    switch (Dim) {
      case 2: {
        rc = qmckl_sm_splitting_core_2(
            context,
            N_updates,
            Updates,
            Updates_index,
            breakdown,
            Slater_inv,
            later_updates,
            later_index, &later, determinant);
        break;
      }
      
      
      case 3: {
        rc = qmckl_sm_splitting_core_3(
            context,
            N_updates,
            Updates,
            Updates_index,
            breakdown,
            Slater_inv,
            later_updates,
            later_index, &later, determinant);
        break;
      }
      
      
      case 4: {
        rc = qmckl_sm_splitting_core_4(
            context,
            N_updates,
            Updates,
            Updates_index,
            breakdown,
            Slater_inv,
            later_updates,
            later_index, &later, determinant);
        break;
      }
      
      
      case 5: {
        rc = qmckl_sm_splitting_core_5(
            context,
            N_updates,
            Updates,
            Updates_index,
            breakdown,
            Slater_inv,
            later_updates,
            later_index, &later, determinant);
        break;
      }
      
      
      case 6: {
        rc = qmckl_sm_splitting_core_6(
            context,
            N_updates,
            Updates,
            Updates_index,
            breakdown,
            Slater_inv,
            later_updates,
            later_index, &later, determinant);
        break;
      }
      
      
      case 7: {
        rc = qmckl_sm_splitting_core_7(
            context,
            N_updates,
            Updates,
            Updates_index,
            breakdown,
            Slater_inv,
            later_updates,
            later_index, &later, determinant);
        break;
      }
      
      
      case 8: {
        rc = qmckl_sm_splitting_core_8(
            context,
            N_updates,
            Updates,
            Updates_index,
            breakdown,
            Slater_inv,
            later_updates,
            later_index, &later, determinant);
        break;
      }
      
      
      case 9: {
        rc = qmckl_sm_splitting_core_9(
            context,
            N_updates,
            Updates,
            Updates_index,
            breakdown,
            Slater_inv,
            later_updates,
            later_index, &later, determinant);
        break;
      }
      
      
      case 10: {
        rc = qmckl_sm_splitting_core_10(
            context,
            N_updates,
            Updates,
            Updates_index,
            breakdown,
            Slater_inv,
            later_updates,
            later_index, &later, determinant);
        break;
      }
      
      
      case 11: {
        rc = qmckl_sm_splitting_core_11(
            context,
            N_updates,
            Updates,
            Updates_index,
            breakdown,
            Slater_inv,
            later_updates,
            later_index, &later, determinant);
        break;
      }
      
      
      case 12: {
        rc = qmckl_sm_splitting_core_12(
            context,
            N_updates,
            Updates,
            Updates_index,
            breakdown,
            Slater_inv,
            later_updates,
            later_index, &later, determinant);
        break;
      }
      
      
      case 13: {
        rc = qmckl_sm_splitting_core_13(
            context,
            N_updates,
            Updates,
            Updates_index,
            breakdown,
            Slater_inv,
            later_updates,
            later_index, &later, determinant);
        break;
      }
      
      
      case 14: {
        rc = qmckl_sm_splitting_core_14(
            context,
            N_updates,
            Updates,
            Updates_index,
            breakdown,
            Slater_inv,
            later_updates,
            later_index, &later, determinant);
        break;
      }
      
      
      case 15: {
        rc = qmckl_sm_splitting_core_15(
            context,
            N_updates,
            Updates,
            Updates_index,
            breakdown,
            Slater_inv,
            later_updates,
            later_index, &later, determinant);
        break;
      }
      
      
      case 16: {
        rc = qmckl_sm_splitting_core_16(
            context,
            N_updates,
            Updates,
            Updates_index,
            breakdown,
            Slater_inv,
            later_updates,
            later_index, &later, determinant);
        break;
      }
      
      
      case 17: {
        rc = qmckl_sm_splitting_core_17(
            context,
            N_updates,
            Updates,
            Updates_index,
            breakdown,
            Slater_inv,
            later_updates,
            later_index, &later, determinant);
        break;
      }
      
      
      case 18: {
        rc = qmckl_sm_splitting_core_18(
            context,
            N_updates,
            Updates,
            Updates_index,
            breakdown,
            Slater_inv,
            later_updates,
            later_index, &later, determinant);
        break;
      }
      
      
      case 19: {
        rc = qmckl_sm_splitting_core_19(
            context,
            N_updates,
            Updates,
            Updates_index,
            breakdown,
            Slater_inv,
            later_updates,
            later_index, &later, determinant);
        break;
      }
      
      
      case 20: {
        rc = qmckl_sm_splitting_core_20(
            context,
            N_updates,
            Updates,
            Updates_index,
            breakdown,
            Slater_inv,
            later_updates,
            later_index, &later, determinant);
        break;
      }
      
      
      case 21: {
        rc = qmckl_sm_splitting_core_21(
            context,
            N_updates,
            Updates,
            Updates_index,
            breakdown,
            Slater_inv,
            later_updates,
            later_index, &later, determinant);
        break;
      }
      default: {
        assert(0 == 1 && "TEMPLATE NOT IMPLEMENTED!");
        break;
      }
    }
  } else {
    rc = qmckl_sm_splitting_core_hpc(
            context, LDS, Dim, N_updates, Updates, Updates_index,
            breakdown, Slater_inv, later_updates,
            later_index, &later, determinant);
  }
  if (rc != QMCKL_SUCCESS) return QMCKL_FAILURE;

  if (later > 0) {
    qmckl_exit_code rc = qmckl_sm_splitting_hpc(
                            context, LDS, Dim, later,
                            later_updates, later_index,
                            breakdown, Slater_inv, determinant);
    if (rc != QMCKL_SUCCESS) return QMCKL_FAILURE;
  }

  return QMCKL_SUCCESS;
}

qmckl_exit_code qmckl_sm_splitting(
    const qmckl_context context,
    const uint64_t LDS,
    const uint64_t Dim,
    const uint64_t N_updates,
    const double* Updates,
    const uint64_t* Updates_index,
    const double breakdown,
    double* Slater_inv,
    double* determinant) {

  if (qmckl_context_check(context) == QMCKL_NULL_CONTEXT) {
    return qmckl_failwith(
        context,
        QMCKL_NULL_CONTEXT,
        "qmckl_sm_splitting",
        NULL);
  }
  #ifdef HAVE_HPC__BROKEN_WITH_CRAY
    return qmckl_sm_splitting_hpc(
        context,
        LDS,
        Dim,
        N_updates,
        Updates,
        Updates_index,
        breakdown,
        Slater_inv,
        determinant);
  #else
    return qmckl_sm_splitting_doc(
        context,
        LDS,
        Dim,
        N_updates,
        Updates,
        Updates_index,
        breakdown,
        Slater_inv,
        determinant);
  #endif 

}
