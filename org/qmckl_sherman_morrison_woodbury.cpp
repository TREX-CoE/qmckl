#include <cmath>

// Sherman-Morrison-Woodbury break-down threshold
#ifndef THRESHOLD
#define THRESHOLD 1e-3
#endif
static double threshold();

// Naïve Sherman Morrison
bool qmckl_sherman_morrison(double *Slater_inv, unsigned int Dim, unsigned int N_updates,
         double *Updates, unsigned int *Updates_index);

// Woodbury 2x2 kernel
bool qmckl_woodbury_2(double *Slater_inv, const unsigned int Dim, double *Updates,
         const unsigned int *Updates_index);

// Woodbury 3x3 kernel
bool qmckl_woodbury_3(double *Slater_inv, const unsigned int Dim, double *Updates,
         const unsigned int *Updates_index);

// Sherman Morrison, with J. Slagel splitting (caller function)
void qmckl_sherman_morrison_splitting(double *Slater_inv, unsigned int Dim, unsigned int N_updates,
         double *Updates, unsigned int *Updates_index);

// Sherman Morrison, with J. Slagel splitting
// http://hdl.handle.net/10919/52966
static void slagel_splitting(double *Slater_inv, unsigned int Dim, unsigned int N_updates,
             double *Updates, unsigned int *Updates_index,
             double *later_updates, unsigned int *later_index,
             unsigned int *later);

// Mixed Sherman-Morrison-Woodbury kernel using
// Woodbury 2x2 and Sherman-Morrison with update-splitting
void qmckl_sherman_morrison_woodbury_2(double *Slater_inv, const unsigned int Dim,
           const unsigned int N_updates, double *Updates,
           unsigned int *Updates_index);

// Mixed Sherman-Morrison-Woodbury kernel using
// Woodbury 3x3, Woodbury 2x2 and Sherman-Morrison with update-splitting
void qmckl_sherman_morrison_woodbury_3(double *Slater_inv, const unsigned int Dim,
           const unsigned int N_updates, double *Updates,
           unsigned int *Updates_index);



// Sherman-Morrison-Woodbury break-down threshold
static double threshold() {
  const double threshold = THRESHOLD;
#ifdef DEBUG
  std::cerr << "Break-down threshold set to: " << threshold << std::endl;
#endif
  return threshold;
}

// Naïve Sherman Morrison
bool qmckl_sherman_morrison(double *Slater_inv, unsigned int Dim, unsigned int N_updates,
         double *Updates, unsigned int *Updates_index) {
#ifdef DEBUG
  std::cerr << "Called qmckl_sherman_morrison with " << N_updates << " updates" << std::endl;
#endif

  double C[Dim];
  double D[Dim];

  unsigned int l = 0;
  // For each update
  while (l < N_updates) {
    // C = A^{-1} x U_l
    for (unsigned int i = 0; i < Dim; i++) {
      C[i] = 0;
      for (unsigned int j = 0; j < Dim; j++) {
        C[i] += Slater_inv[i * Dim + j] * Updates[l * Dim + j];
      }
    }

    // Denominator
    double den = 1 + C[Updates_index[l] - 1];
    if (std::fabs(den) < threshold()) {
      return false;
    }
    double iden = 1 / den;

    // D = v^T x A^{-1}
    for (unsigned int j = 0; j < Dim; j++) {
      D[j] = Slater_inv[(Updates_index[l] - 1) * Dim + j];
    }

    // A^{-1} = A^{-1} - C x D / den
    for (unsigned int i = 0; i < Dim; i++) {
      for (unsigned int j = 0; j < Dim; j++) {
        double update = C[i] * D[j] * iden;
        Slater_inv[i * Dim + j] -= update;
      }
    }

    l += 1;
  }
  return true;
}

// Woodbury 2x2 kernel
bool qmckl_woodbury_2(double *Slater_inv, const unsigned int Dim, double *Updates,
         const unsigned int *Updates_index) {
/*
    C := S^{-1} * U,    dim x 2
    B := 1 + V * C,     2 x 2
    D := V * S^{-1},    2 x dim
*/
#ifdef DEBUG
  std::cerr << "Called Woodbury 2x2 kernel" << std::endl;
#endif

  const unsigned int row1 = (Updates_index[0] - 1);
  const unsigned int row2 = (Updates_index[1] - 1);

  // Compute C = S_inv * U  !! NON-STANDARD MATRIX MULTIPLICATION BECAUSE
  // OF LAYOUT OF 'Updates' !!
  double C[2 * Dim];
  for (unsigned int i = 0; i < Dim; i++) {
    for (unsigned int j = 0; j < 2; j++) {
      C[i * 2 + j] = 0;
      for (unsigned int k = 0; k < Dim; k++) {
        C[i * 2 + j] += Slater_inv[i * Dim + k] * Updates[Dim * j + k];
      }
    }
  }

  // Compute B = 1 + V * C
  const double B0 = C[row1 * 2] + 1;
  const double B1 = C[row1 * 2 + 1];
  const double B2 = C[row2 * 2];
  const double B3 = C[row2 * 2 + 1] + 1;

  // Check if determinant of inverted matrix is not zero
  double det = B0 * B3 - B1 * B2;
  if (std::fabs(det) < threshold()) {
    return false;
  }

  // Compute B^{-1} with explicit formula for 2x2 inversion
  double Binv[4], idet = 1.0 / det;
  Binv[0] = idet * B3;
  Binv[1] = -1.0 * idet * B1;
  Binv[2] = -1.0 * idet * B2;
  Binv[3] = idet * B0;

  // Compute tmp = B^{-1} x (V.S^{-1})
  double tmp[2 * Dim];
  for (unsigned int i = 0; i < 2; i++) {
    for (unsigned int j = 0; j < Dim; j++) {
      tmp[i * Dim + j] = Binv[i * 2] * Slater_inv[row1 * Dim + j];
      tmp[i * Dim + j] += Binv[i * 2 + 1] * Slater_inv[row2 * Dim + j];
    }
  }

  // Compute (S + U V)^{-1} = S^{-1} - C x tmp
  for (unsigned int i = 0; i < Dim; i++) {
    for (unsigned int j = 0; j < Dim; j++) {
      Slater_inv[i * Dim + j] -= C[i * 2] * tmp[j];
      Slater_inv[i * Dim + j] -= C[i * 2 + 1] * tmp[Dim + j];
    }
  }

  return true;
}

// Woodbury 3x3 kernel
bool qmckl_woodbury_3(double *Slater_inv, const unsigned int Dim, double *Updates,
         const unsigned int *Updates_index) {
/*
    C := S^{-1} * U,    dim x 3
    B := 1 + V * C,     3 x 3
    D := V * S^{-1},    3 x dim
*/
#ifdef DEBUG
  std::cerr << "Called Woodbury 3x3 kernel" << std::endl;
#endif

  const unsigned int row1 = (Updates_index[0] - 1);
  const unsigned int row2 = (Updates_index[1] - 1);
  const unsigned int row3 = (Updates_index[2] - 1);

  // Compute C = S_inv * U  !! NON-STANDARD MATRIX MULTIPLICATION BECAUSE
  // OF LAYOUT OF 'Updates' !!
  double C[3 * Dim];
  for (unsigned int i = 0; i < Dim; i++) {
    for (unsigned int j = 0; j < 3; j++) {
      C[i * 3 + j] = 0;
      for (unsigned int k = 0; k < Dim; k++) {
        C[i * 3 + j] += Slater_inv[i * Dim + k] * Updates[Dim * j + k];
      }
    }
  }

  // Compute B = 1 + V.C
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

  if (std::fabs(det) < threshold()) {
    return false;
  }

  // Compute B^{-1} with explicit formula for 3x3 inversion
  double Binv[9], idet = 1.0 / det;
  Binv[0] = (B4 * B8 - B7 * B5) * idet;
  Binv[1] = -(B1 * B8 - B7 * B2) * idet;
  Binv[2] = (B1 * B5 - B4 * B2) * idet;
  Binv[3] = -(B3 * B8 - B6 * B5) * idet;
  Binv[4] = (B0 * B8 - B6 * B2) * idet;
  Binv[5] = -(B0 * B5 - B3 * B2) * idet;
  Binv[6] = (B3 * B7 - B6 * B4) * idet;
  Binv[7] = -(B0 * B7 - B6 * B1) * idet;
  Binv[8] = (B0 * B4 - B3 * B1) * idet;

  // Compute tmp = B^{-1} x (V.S^{-1})
  double tmp[3 * Dim];
  for (unsigned int i = 0; i < 3; i++) {
    for (unsigned int j = 0; j < Dim; j++) {
      tmp[i * Dim + j] = Binv[i * 3] * Slater_inv[row1 * Dim + j];
      tmp[i * Dim + j] += Binv[i * 3 + 1] * Slater_inv[row2 * Dim + j];
      tmp[i * Dim + j] += Binv[i * 3 + 2] * Slater_inv[row3 * Dim + j];
    }
  }

  // Compute (S + U V)^{-1} = S^{-1} - C x tmp
  for (unsigned int i = 0; i < Dim; i++) {
    for (unsigned int j = 0; j < Dim; j++) {
      Slater_inv[i * Dim + j] -= C[i * 3] * tmp[j];
      Slater_inv[i * Dim + j] -= C[i * 3 + 1] * tmp[Dim + j];
      Slater_inv[i * Dim + j] -= C[i * 3 + 2] * tmp[2 * Dim + j];
    }
  }

  return true;
}

// Sherman Morrison, with J. Slagel splitting (caller function)
// http://hdl.handle.net/10919/52966
void qmckl_sherman_morrison_splitting(double *Slater_inv, unsigned int Dim, unsigned int N_updates,
         double *Updates, unsigned int *Updates_index) {
#ifdef DEBUG
  std::cerr << "Called qmckl_sherman_morrison_splitting with " << N_updates << " updates" << std::endl;
#endif

  double later_updates[Dim * N_updates];
  unsigned int later_index[N_updates];
  unsigned int later = 0;

  slagel_splitting(Slater_inv, Dim, N_updates, Updates, Updates_index, later_updates,
          later_index, &later);

  if (later > 0) {
    qmckl_sherman_morrison_splitting(Slater_inv, Dim, later, later_updates, later_index);
  }
}

// Sherman Morrison, with J. Slagel splitting
// http://hdl.handle.net/10919/52966
static void slagel_splitting(double *Slater_inv, unsigned int Dim, unsigned int N_updates,
             double *Updates, unsigned int *Updates_index,
             double *later_updates, unsigned int *later_index,
             unsigned int *later) {
#ifdef DEBUG
  std::cerr << "Called slagel_splitting with " << N_updates << " updates" << std::endl;
#endif

  double C[Dim];
  double D[Dim];

  unsigned int l = 0;
  // For each update
  while (l < N_updates) {
    // C = S^{-1} x U_l
    for (unsigned int i = 0; i < Dim; i++) {
      C[i] = 0;
      for (unsigned int j = 0; j < Dim; j++) {
        C[i] += Slater_inv[i * Dim + j] * Updates[l * Dim + j];
      }
    }

    // Denominator
    double den = 1 + C[Updates_index[l] - 1];
    if (std::fabs(den) < threshold()) {

      // U_l = U_l / 2 (do the split)
      for (unsigned int i = 0; i < Dim; i++) {
        later_updates[*later * Dim + i] = Updates[l * Dim + i] / 2.0;
        C[i] /= 2.0;
      }
      later_index[*later] = Updates_index[l];
      (*later)++;

      den = 1 + C[Updates_index[l] - 1];
    }
    double iden = 1 / den;

    // D = v^T x S^{-1}
    for (unsigned int j = 0; j < Dim; j++) {
      D[j] = Slater_inv[(Updates_index[l] - 1) * Dim + j];
    }

    // S^{-1} = S^{-1} - C x D / den
    for (unsigned int i = 0; i < Dim; i++) {
      for (unsigned int j = 0; j < Dim; j++) {
        double update = C[i] * D[j] * iden;
        Slater_inv[i * Dim + j] -= update;
      }
    }
    l += 1;
  }
}

// Sherman-Morrison-Woodbury 2x2 kernel
// qmckl_woodbury_2, slagel_splitting mixing scheme
void qmckl_sherman_morrison_woodbury_2(double *Slater_inv, const unsigned int Dim,
           const unsigned int N_updates, double *Updates,
           unsigned int *Updates_index) {
#ifdef DEBUG
  std::cerr << "Called qmckl_sherman_morrison_woodbury_2 with " << N_updates
            << " updates" << std::endl;
#endif

  unsigned int n_of_2blocks = N_updates / 2;
  unsigned int remainder = N_updates % 2;
  unsigned int length_2block = 2 * Dim;
  unsigned int length_1block = 1 * Dim;

  // Apply first 2*n_of_2blocks updates in n_of_2blocks blocks of 2 updates with
  // Woodbury 2x2 kernel
  double later_updates[Dim * N_updates];
  unsigned int later_index[N_updates];
  unsigned int later = 0;
  if (n_of_2blocks > 0) {
    for (unsigned int i = 0; i < n_of_2blocks; i++) {
      double *Updates_2block = &Updates[i * length_2block];
      unsigned int *Updates_index_2block = &Updates_index[i * 2];
      bool ok;
      ok = qmckl_woodbury_2(Slater_inv, Dim, Updates_2block, Updates_index_2block);
      if (!ok) { // Send the entire block to slagel_splitting
        unsigned int l = 0;
        slagel_splitting(Slater_inv, Dim, 2, Updates_2block, Updates_index_2block,
                later_updates + (Dim * later), later_index + later, &l);
        later = later + l;
      }
    }
  }

  if (remainder == 1) { // Apply last remaining update with slagel_splitting
    double *Updates_1block = &Updates[n_of_2blocks * length_2block];
    unsigned int *Updates_index_1block = &Updates_index[2 * n_of_2blocks];
    unsigned int l = 0;
    slagel_splitting(Slater_inv, Dim, 1, Updates_1block, Updates_index_1block,
            later_updates + (Dim * later), later_index + later, &l);
    later = later + l;
  }

  if (later > 0) {
    qmckl_sherman_morrison_splitting(Slater_inv, Dim, later, later_updates, later_index);
  }
}

// Sherman-Morrison-Woodbury 3x3 kernel
// qmckl_woodbury_2, qmckl_woodbury_3, slagel_splitting mixing scheme
void qmckl_sherman_morrison_woodbury_3(double *Slater_inv, const unsigned int Dim,
           const unsigned int N_updates, double *Updates,
           unsigned int *Updates_index) {
#ifdef DEBUG
  std::cerr << "Called qmckl_sherman_morrison_woodbury_3 with " << N_updates
            << " updates" << std::endl;
#endif

  unsigned int n_of_3blocks = N_updates / 3;
  unsigned int remainder = N_updates % 3;
  unsigned int length_3block = 3 * Dim;
  unsigned int length_2block = 2 * Dim;
  unsigned int length_1block = 1 * Dim;

  // Apply first 3*n_of_3blocks updates in n_of_3blocks blocks of 3 updates with
  // Woodbury 3x3 kernel
  double later_updates[Dim * N_updates];
  unsigned int later_index[N_updates];
  unsigned int later = 0;
  if (n_of_3blocks > 0) {
    for (unsigned int i = 0; i < n_of_3blocks; i++) {
      double *Updates_3block = &Updates[i * length_3block];
      unsigned int *Updates_index_3block = &Updates_index[i * 3];
      bool ok;
      ok = qmckl_woodbury_3(Slater_inv, Dim, Updates_3block, Updates_index_3block);
      if (!ok) { // Send the entire block to slagel_splitting
        unsigned int l = 0;
        slagel_splitting(Slater_inv, Dim, 3, Updates_3block, Updates_index_3block,
                later_updates + (Dim * later), later_index + later, &l);
        later = later + l;
      }
    }
  }

  if (remainder == 2) { // Apply last remaining block of 2 updates with Woodbury 2x2 kernel
    double *Updates_2block = &Updates[n_of_3blocks * length_3block];
    unsigned int *Updates_index_2block = &Updates_index[3 * n_of_3blocks];
    bool ok;
    ok = qmckl_woodbury_2(Slater_inv, Dim, Updates_2block, Updates_index_2block);
    if (!ok) { // Send the entire block to slagel_splitting
      unsigned int l = 0;
      slagel_splitting(Slater_inv, Dim, 2, Updates_2block, Updates_index_2block,
              later_updates + (Dim * later), later_index + later, &l);
      later = later + l;
    }
  }
  else if (remainder == 1) { // Apply last remaining update with slagel_splitting
    double *Updates_1block = &Updates[n_of_3blocks * length_3block];
    unsigned int *Updates_index_1block = &Updates_index[3 * n_of_3blocks];
    unsigned int l = 0;
    slagel_splitting(Slater_inv, Dim, 1, Updates_1block, Updates_index_1block,
            later_updates + (Dim * later), later_index + later, &l);
    later = later + l;
  }

  if (later > 0) {
    qmckl_sherman_morrison_splitting(Slater_inv, Dim, later, later_updates, later_index);
  }
}
