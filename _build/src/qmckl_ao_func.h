

/* To set the basis set, all the following functions need to be */
/* called. */


qmckl_exit_code
qmckl_set_ao_basis_type (qmckl_context context,
                         const char basis_type);

qmckl_exit_code
qmckl_set_ao_basis_shell_num (qmckl_context context,
                              const int64_t shell_num);

qmckl_exit_code
qmckl_set_ao_basis_prim_num (qmckl_context context,
                             const int64_t prim_num);

qmckl_exit_code
qmckl_set_ao_basis_nucleus_shell_num (qmckl_context context,
                                      const int64_t* nucleus_shell_num,
                                      const int64_t size_max);

qmckl_exit_code
qmckl_set_ao_basis_nucleus_index (qmckl_context context,
                                  const int64_t* nucleus_index,
                                  const int64_t size_max);

qmckl_exit_code
qmckl_set_ao_basis_shell_ang_mom (qmckl_context context,
                                  const int32_t* shell_ang_mom,
                                  const int64_t size_max);

qmckl_exit_code
qmckl_set_ao_basis_shell_prim_num (qmckl_context context,
                                   const int64_t* shell_prim_num,
                                   const int64_t size_max);

qmckl_exit_code
qmckl_set_ao_basis_shell_prim_index (qmckl_context context,
                                     const int64_t* shell_prim_index,
                                     const int64_t size_max);

qmckl_exit_code
qmckl_set_ao_basis_shell_factor (qmckl_context context,
                                 const double* shell_factor,
                                 const int64_t size_max);

qmckl_exit_code
qmckl_set_ao_basis_exponent (qmckl_context context,
                             const double* exponent,
                             const int64_t size_max);

qmckl_exit_code
qmckl_set_ao_basis_coefficient (qmckl_context context,
                                const double* coefficient,
                                const int64_t size_max);

qmckl_exit_code
qmckl_set_ao_basis_prim_factor (qmckl_context context,
                                const double* prim_factor,
                                const int64_t size_max);

qmckl_exit_code
qmckl_set_ao_basis_ao_num (qmckl_context context,
                           const int64_t ao_num);

qmckl_exit_code
qmckl_set_ao_basis_ao_factor (qmckl_context context,
                              const double* ao_factor,
                              const int64_t size_max);

qmckl_exit_code
qmckl_set_ao_basis_r_power (qmckl_context context,
                            const int32_t* r_power,
                            const int64_t size_max);

qmckl_exit_code
qmckl_set_ao_basis_cartesian (qmckl_context context,
                              const bool cartesian);

/* C interface */


qmckl_exit_code
qmckl_get_ao_basis_type (const qmckl_context context,
                         char* const basis_type);

qmckl_exit_code
qmckl_get_ao_basis_shell_num (const qmckl_context context,
                              int64_t* const shell_num);

qmckl_exit_code
qmckl_get_ao_basis_prim_num (const qmckl_context context,
                             int64_t* const prim_num);

qmckl_exit_code
qmckl_get_ao_basis_nucleus_shell_num (const qmckl_context context,
                                      int64_t* const nucleus_shell_num,
                                      const int64_t size_max);

qmckl_exit_code
qmckl_get_ao_basis_nucleus_index (const qmckl_context context,
                                  int64_t* const nucleus_index,
                                  const int64_t size_max);

qmckl_exit_code
qmckl_get_ao_basis_r_power (const qmckl_context context,
                            int32_t* const r_power,
                            const int64_t size_max);

qmckl_exit_code
qmckl_get_ao_basis_shell_ang_mom (const qmckl_context context,
                                  int32_t* const shell_ang_mom,
                                  const int64_t size_max);

qmckl_exit_code
qmckl_get_ao_basis_shell_prim_num (const qmckl_context context,
                                   int64_t* const shell_prim_num,
                                   const int64_t size_max);

qmckl_exit_code
qmckl_get_ao_basis_shell_prim_index (const qmckl_context context,
                                     int64_t* const shell_prim_index,
                                     const int64_t size_max);

qmckl_exit_code
qmckl_get_ao_basis_shell_factor (const qmckl_context context,
                                 double*  const shell_factor,
                                 const int64_t size_max);

qmckl_exit_code
qmckl_get_ao_basis_exponent (const qmckl_context context,
                             double*  const exponent,
                             const int64_t size_max);

qmckl_exit_code
qmckl_get_ao_basis_coefficient (const qmckl_context context,
                                double*  const coefficient,
                                const int64_t size_max);

qmckl_exit_code
qmckl_get_ao_basis_prim_factor (const qmckl_context context,
                                double*  const prim_factor,
                                const int64_t size_max);

qmckl_exit_code
qmckl_get_ao_basis_ao_num (const qmckl_context context,
                           int64_t* const ao_num);

qmckl_exit_code
qmckl_get_ao_basis_ao_factor (const qmckl_context context,
                              double* const ao_factor,
                              const int64_t size_max);




/* When all the data for the AOs have been provided, the following */
/* function returns ~true~. */


bool qmckl_ao_basis_provided (const qmckl_context context);

/* Access functions */


qmckl_exit_code
qmckl_get_ao_basis_primitive_vgl (qmckl_context context,
                                  double* const primitive_vgl,
                                  const int64_t size_max);

qmckl_exit_code
qmckl_get_ao_basis_shell_vgl (qmckl_context context,
                              double* const shell_vgl,
                              const int64_t size_max);

qmckl_exit_code
qmckl_get_ao_basis_ao_vgl (qmckl_context context,
                           double* const ao_vgl,
                           const int64_t size_max);



/* Uses the given array to compute the VGL. */


qmckl_exit_code
qmckl_get_ao_basis_ao_vgl_inplace (qmckl_context context,
                                   double* const ao_vgl,
                                   const int64_t size_max);

qmckl_exit_code
qmckl_get_ao_basis_ao_value (qmckl_context context,
                             double* const ao_value,
                             const int64_t size_max);



/* Uses the given array to compute the value. */


qmckl_exit_code
qmckl_get_ao_basis_ao_value_inplace (qmckl_context context,
                                     double* const ao_value,
                                     const int64_t size_max);

qmckl_exit_code
qmckl_ao_gaussian_vgl(const qmckl_context context,
                      const double *X,
                      const double *R,
                      const int64_t *n,
                      const int64_t *A,
                      const double *VGL,
                      const int64_t ldv);

qmckl_exit_code
qmckl_ao_slater_vgl(const qmckl_context context,
                    const double *X,
                    const double *R,
                    const int64_t *num_slater,
                    const int64_t *N,
                    const double *A,
                    const double *VGL,
                    const int64_t ldv);

qmckl_exit_code
qmckl_ao_slater_vgl_hpc(const qmckl_context context,
                         const double *X,
                         const double *R,
                         const int64_t *num_slater,
                         const int64_t *N,
                         const double *A,
                         double *VGL,
                         const int64_t ldv);

/* Computation of primitives */
/*    :PROPERTIES: */
/*    :Name:     qmckl_compute_ao_basis_primitive_gaussian_vgl */
/*    :CRetType: qmckl_exit_code */
/*    :FRetType: qmckl_exit_code */
/*    :END: */

/*    #+NAME: qmckl_ao_basis_primitive_gaussian_vgl_args */
/*    | Variable             | Type                             | In/Out | Description                                      | */
/*    |----------------------+----------------------------------+--------+--------------------------------------------------| */
/*    | ~context~            | ~qmckl_context~                  | in     | Global state                                     | */
/*    | ~prim_num~           | ~int64_t~                        | in     | Number of primitives                             | */
/*    | ~point_num~          | ~int64_t~                        | in     | Number of points considered                      | */
/*    | ~nucl_num~           | ~int64_t~                        | in     | Number of nuclei                                 | */
/*    | ~nucleus_prim_index~ | ~int64_t[nucl_num+1]~            | in     | Index of the 1st primitive of each nucleus       | */
/*    | ~coord~              | ~double[3][point_num]~           | in     | Coordinates                                      | */
/*    | ~nucl_coord~         | ~double[3][nucl_num]~            | in     | Nuclear coordinates                              | */
/*    | ~expo~               | ~double[prim_num]~               | in     | Exponents of the primitives                      | */
/*    | ~primitive_vgl~      | ~double[point_num][5][prim_num]~ | out    | Value, gradients and Laplacian of the primitives | */

/*    #+CALL: generate_private_c_header(table=qmckl_ao_basis_primitive_gaussian_vgl_args,rettyp=get_value("CRetType"),fname="qmckl_compute_ao_basis_primitive_gaussian_vgl")) */

/*    #+RESULTS: */

qmckl_exit_code qmckl_compute_ao_basis_primitive_gaussian_vgl (
      const qmckl_context context,
      const int64_t prim_num,
      const int64_t point_num,
      const int64_t nucl_num,
      const int64_t* nucleus_prim_index,
      const double* coord,
      const double* nucl_coord,
      const double* expo,
      double* const primitive_vgl );

/* Gaussians */
/*    :PROPERTIES: */
/*    :Name:     qmckl_compute_ao_basis_shell_gaussian_vgl */
/*    :CRetType: qmckl_exit_code */
/*    :FRetType: qmckl_exit_code */
/*    :END: */

/*    #+NAME: qmckl_ao_basis_shell_gaussian_vgl_args */
/*    | Variable            | Type                              | In/Out | Description                                  | */
/*    |---------------------+-----------------------------------+--------+----------------------------------------------| */
/*    | ~context~           | ~qmckl_context~                   | in     | Global state                                 | */
/*    | ~prim_num~          | ~int64_t~                         | in     | Number of primitives                         | */
/*    | ~shell_num~         | ~int64_t~                         | in     | Number of shells                             | */
/*    | ~point_num~         | ~int64_t~                         | in     | Number of points                             | */
/*    | ~nucl_num~          | ~int64_t~                         | in     | Number of nuclei                             | */
/*    | ~nucleus_shell_num~ | ~int64_t[nucl_num]~               | in     | Number of shells for each nucleus            | */
/*    | ~nucleus_index~     | ~int64_t[nucl_num]~               | in     | Index of the 1st shell of each nucleus       | */
/*    | ~nucleus_range~     | ~double[nucl_num]~                | in     | Range of the nucleus                         | */
/*    | ~shell_prim_index~  | ~int64_t[shell_num]~              | in     | Index of the 1st primitive of each shell     | */
/*    | ~shell_prim_num~    | ~int64_t[shell_num]~              | in     | Number of primitives per shell               | */
/*    | ~coord~             | ~double[3][point_num]~            | in     | Coordinates                                  | */
/*    | ~nucl_coord~        | ~double[3][nucl_num]~             | in     | Nuclear coordinates                          | */
/*    | ~expo~              | ~double[prim_num]~                | in     | Exponents of the primitives                  | */
/*    | ~coef_normalized~   | ~double[prim_num]~                | in     | Coefficients of the primitives               | */
/*    | ~shell_vgl~         | ~double[point_num][5][shell_num]~ | out    | Value, gradients and Laplacian of the shells | */

/*    #+CALL: generate_private_c_header(table=qmckl_ao_basis_shell_gaussian_vgl_args,rettyp=get_value("CRetType"),fname="qmckl_compute_ao_basis_shell_gaussian_vgl")) */

/*    #+RESULTS: */

qmckl_exit_code qmckl_compute_ao_basis_shell_gaussian_vgl (
      const qmckl_context context,
      const int64_t prim_num,
      const int64_t shell_num,
      const int64_t point_num,
      const int64_t nucl_num,
      const int64_t* nucleus_shell_num,
      const int64_t* nucleus_index,
      const double* nucleus_range,
      const int64_t* shell_prim_index,
      const int64_t* shell_prim_num,
      const double* coord,
      const double* nucl_coord,
      const double* expo,
      const double* coef_normalized,
      double* const shell_vgl );

/* Slater */
/*    :PROPERTIES: */
/*    :Name:     qmckl_compute_ao_basis_shell_slater_vgl */
/*    :CRetType: qmckl_exit_code */
/*    :FRetType: qmckl_exit_code */
/*    :END: */

/*    #+NAME: qmckl_ao_basis_shell_slater_vgl_args */
/*    | Variable            | Type                              | In/Out | Description                                  | */
/*    |---------------------+-----------------------------------+--------+----------------------------------------------| */
/*    | ~context~           | ~qmckl_context~                   | in     | Global state                                 | */
/*    | ~prim_num~          | ~int64_t~                         | in     | Number of primitives                         | */
/*    | ~shell_num~         | ~int64_t~                         | in     | Number of shells                             | */
/*    | ~point_num~         | ~int64_t~                         | in     | Number of points                             | */
/*    | ~nucl_num~          | ~int64_t~                         | in     | Number of nuclei                             | */
/*    | ~nucleus_shell_num~ | ~int64_t[nucl_num]~               | in     | Number of shells per nucleus                 | */
/*    | ~nucleus_index~     | ~int64_t[nucl_num]~               | in     | Index of the first shell of each nucleus     | */
/*    | ~nucleus_range~     | ~double[nucl_num]~                | in     | Range beyond which all shells are zero       | */
/*    | ~shell_prim_index~  | ~int64_t[shell_num]~              | in     | Index of the first primitive of each shell   | */
/*    | ~shell_prim_num~    | ~int64_t[shell_num]~              | in     | Number of primitives per shell               | */
/*    | ~r_power~           | ~int32_t[shell_num]~              | in     | Power of r prefix                            | */
/*    | ~coord~             | ~double[3][point_num]~            | in     | Coordinates of the points                    | */
/*    | ~nucl_coord~        | ~double[3][nucl_num]~             | in     | Nuclear coordinates                          | */
/*    | ~expo~              | ~double[prim_num]~                | in     | Exponents of the primitives                  | */
/*    | ~coef_normalized~   | ~double[prim_num]~                | in     | Coefficients of the primitives               | */
/*    | ~shell_vgl~         | ~double[point_num][5][shell_num]~ | out    | Value, gradients and Laplacian of the shells | */

/*    #+CALL: generate_private_c_header(table=qmckl_ao_basis_shell_slater_vgl_args,rettyp=get_value("CRetType"),fname="qmckl_compute_ao_basis_shell_slater_vgl")) */

/*    #+RESULTS: */

qmckl_exit_code qmckl_compute_ao_basis_shell_slater_vgl (
      const qmckl_context context,
      const int64_t prim_num,
      const int64_t shell_num,
      const int64_t point_num,
      const int64_t nucl_num,
      const int64_t* nucleus_shell_num,
      const int64_t* nucleus_index,
      const double* nucleus_range,
      const int64_t* shell_prim_index,
      const int64_t* shell_prim_num,
      const int32_t* r_power,
      const double* coord,
      const double* nucl_coord,
      const double* expo,
      const double* coef_normalized,
      double* const shell_vgl );

/* Get */


qmckl_exit_code
qmckl_get_ao_basis_shell_hessian (qmckl_context context,
                              double* const shell_hessian,
                              const int64_t size_max);

/* Computation of shell Hessians with Gaussian functions */
/*    :PROPERTIES: */
/*    :Name:     qmckl_compute_ao_basis_shell_gaussian_hessian */
/*    :CRetType: qmckl_exit_code */
/*    :FRetType: qmckl_exit_code */
/*    :END: */

/*    #+NAME: qmckl_ao_basis_shell_gaussian_hessian_args */
/*    | Variable            | Type                                 | In/Out | Description                              | */
/*    |---------------------+--------------------------------------+--------+------------------------------------------| */
/*    | ~context~           | ~qmckl_context~                      | in     | Global state                             | */
/*    | ~prim_num~          | ~int64_t~                            | in     | Number of primitives                     | */
/*    | ~shell_num~         | ~int64_t~                            | in     | Number of shells                         | */
/*    | ~point_num~         | ~int64_t~                            | in     | Number of points                         | */
/*    | ~nucl_num~          | ~int64_t~                            | in     | Number of nuclei                         | */
/*    | ~nucleus_shell_num~ | ~int64_t[nucl_num]~                  | in     | Number of shells for each nucleus        | */
/*    | ~nucleus_index~     | ~int64_t[nucl_num]~                  | in     | Index of the 1st shell of each nucleus   | */
/*    | ~nucleus_range~     | ~double[nucl_num]~                   | in     | Range of the nucleus                     | */
/*    | ~shell_prim_index~  | ~int64_t[shell_num]~                 | in     | Index of the 1st primitive of each shell | */
/*    | ~shell_prim_num~    | ~int64_t[shell_num]~                 | in     | Number of primitives per shell           | */
/*    | ~coord~             | ~double[3][point_num]~               | in     | Coordinates                              | */
/*    | ~nucl_coord~        | ~double[3][nucl_num]~                | in     | Nuclear coordinates                      | */
/*    | ~expo~              | ~double[prim_num]~                   | in     | Exponents of the primitives              | */
/*    | ~coef_normalized~   | ~double[prim_num]~                   | in     | Coefficients of the primitives           | */
/*    | ~shell_hessian~     | ~double[point_num][3][4][shell_num]~ | out    | Hessian of the shells                    | */
/*    |---------------------+--------------------------------------+--------+------------------------------------------| */


qmckl_exit_code qmckl_compute_ao_basis_shell_gaussian_hessian (
      const qmckl_context context,
      const int64_t prim_num,
      const int64_t shell_num,
      const int64_t point_num,
      const int64_t nucl_num,
      const int64_t* nucleus_shell_num,
      const int64_t* nucleus_index,
      const double* nucleus_range,
      const int64_t* shell_prim_index,
      const int64_t* shell_prim_num,
      const double* coord,
      const double* nucl_coord,
      const double* expo,
      const double* coef_normalized,
      double* const shell_hessian );

/* TODO Computation of shell Hessians with Slater functions */
/*    :PROPERTIES: */
/*    :Name:     qmckl_compute_ao_basis_shell_slater_hessian */
/*    :CRetType: qmckl_exit_code */
/*    :FRetType: qmckl_exit_code */
/*    :END: */

/*    #+NAME: qmckl_ao_basis_shell_slater_hessian_args */
/*    | Variable            | Type                                 | In/Out | Description                                | */
/*    |---------------------+--------------------------------------+--------+--------------------------------------------| */
/*    | ~context~           | ~qmckl_context~                      | in     | Global state                               | */
/*    | ~prim_num~          | ~int64_t~                            | in     | Number of primitives                       | */
/*    | ~shell_num~         | ~int64_t~                            | in     | Number of shells                           | */
/*    | ~point_num~         | ~int64_t~                            | in     | Number of points                           | */
/*    | ~nucl_num~          | ~int64_t~                            | in     | Number of nuclei                           | */
/*    | ~nucleus_shell_num~ | ~int64_t[nucl_num]~                  | in     | Number of shells per nucleus               | */
/*    | ~nucleus_index~     | ~int64_t[nucl_num]~                  | in     | Index of the first shell of each nucleus   | */
/*    | ~nucleus_range~     | ~double[nucl_num]~                   | in     | Range beyond which all shells are zero     | */
/*    | ~shell_prim_index~  | ~int64_t[shell_num]~                 | in     | Index of the first primitive of each shell | */
/*    | ~shell_prim_num~    | ~int64_t[shell_num]~                 | in     | Number of primitives per shell             | */
/*    | ~r_power~           | ~int32_t[shell_num]~                 | in     | Power of the r prefix                      | */
/*    | ~coord~             | ~double[3][point_num]~               | in     | Coordinates of the points                  | */
/*    | ~nucl_coord~        | ~double[3][nucl_num]~                | in     | Nuclear coordinates                        | */
/*    | ~expo~              | ~double[prim_num]~                   | in     | Exponents of the primitives                | */
/*    | ~coef_normalized~   | ~double[prim_num]~                   | in     | Coefficients of the primitives             | */
/*    | ~shell_hessian~     | ~double[point_num][3][4][shell_num]~ | out    | Hessian of the shells                      | */
/*    |---------------------+--------------------------------------+--------+--------------------------------------------| */


qmckl_exit_code qmckl_compute_ao_basis_shell_slater_hessian (
      const qmckl_context context,
      const int64_t prim_num,
      const int64_t shell_num,
      const int64_t point_num,
      const int64_t nucl_num,
      const int64_t* nucleus_shell_num,
      const int64_t* nucleus_index,
      const double* nucleus_range,
      const int64_t* shell_prim_index,
      const int64_t* shell_prim_num,
      const int32_t* r_power,
      const double* coord,
      const double* nucl_coord,
      const double* expo,
      const double* coef_normalized,
      double* const shell_hessian );

/* General functions for Powers of $x-X_i$ */
/*    :PROPERTIES: */
/*    :Name:     qmckl_ao_power */
/*    :CRetType: qmckl_exit_code */
/*    :FRetType: qmckl_exit_code */
/*    :END: */

/*    The ~qmckl_ao_power~ function computes all the powers of the ~n~ */
/*    input data up to the given maximum value given in input for each of */
/*    the $n$ points: */

/*    \[ P_{ik} = X_i^k \] */

/*    #+NAME: qmckl_ao_power_args */
/*    | Variable  | Type            | In/Out | Description                                       | */
/*    |-----------+-----------------+--------+---------------------------------------------------| */
/*    | ~context~ | ~qmckl_context~ | in     | Global state                                      | */
/*    | ~n~       | int64_t         | in     | Number of values                                  | */
/*    | ~X~       | double[n]       | in     | Array containing the input values                 | */
/*    | ~LMAX~    | int32_t[n]      | in     | Array containing the maximum power for each value | */
/*    | ~P~       | double[n][ldp]  | out    | Array containing all the powers of ~X~            | */
/*    | ~ldp~     | int64_t         | in     | Leading dimension of array ~P~                    | */

/*    Requirements: */

/*    - ~context~ is not ~QMCKL_NULL_CONTEXT~ */
/*    - ~n~ > 0 */
/*    - ~X~ is allocated with at least $n \times 8$ bytes */
/*    - ~LMAX~ is allocated with at least $n \times 4$ bytes */
/*    - ~P~ is allocated with at least $n \times \max_i \text{LMAX}_i \times 8$ bytes */
/*    - ~LDP~ >= $\max_i$ ~LMAX[i]~ */

/*    #+CALL: generate_c_header(table=qmckl_ao_power_args,rettyp=get_value("CRetType"),fname="qmckl_ao_power") */

/*    #+RESULTS: */

qmckl_exit_code qmckl_ao_power (
      const qmckl_context context,
      const int64_t n,
      const double* X,
      const int32_t* LMAX,
      double* const P,
      const int64_t ldp );

/* General functions for Value, Gradient and Laplacian of a polynomial */
/*    :PROPERTIES: */
/*    :Name:     qmckl_ao_polynomial_vgl */
/*    :CRetType: qmckl_exit_code */
/*    :FRetType: qmckl_exit_code */
/*    :END: */

/*    A polynomial is centered on a nucleus $\mathbf{R}_i$ */

/*    \[ */
/*    P_l(\mathbf{r},\mathbf{R}_i)  =   (x-X_i)^a (y-Y_i)^b (z-Z_i)^c */
/*    \] */

/*    The gradients with respect to electron coordinates are */

/*    \begin{eqnarray*} */
/*    \frac{\partial }{\partial x} P_l\left(\mathbf{r},\mathbf{R}_i \right) & */
/*                   = & a (x-X_i)^{a-1} (y-Y_i)^b (z-Z_i)^c \\ */
/*    \frac{\partial }{\partial y} P_l\left(\mathbf{r},\mathbf{R}_i \right) & */
/*                   = & b (x-X_i)^a (y-Y_i)^{b-1} (z-Z_i)^c \\ */
/*    \frac{\partial }{\partial z} P_l\left(\mathbf{r},\mathbf{R}_i \right) & */
/*                   = & c (x-X_i)^a (y-Y_i)^b (z-Z_i)^{c-1} \\ */
/*    \end{eqnarray*} */

/*    and the Laplacian is */

/*    \begin{eqnarray*} */
/*    \left( \frac{\partial }{\partial x^2} + */
/*               \frac{\partial }{\partial y^2} + */
/*               \frac{\partial }{\partial z^2} \right) P_l */
/*               \left(\mathbf{r},\mathbf{R}_i \right) &  = & */
/*             a(a-1) (x-X_i)^{a-2} (y-Y_i)^b (z-Z_i)^c + \\ */
/*          && b(b-1) (x-X_i)^a (y-Y_i)^{b-2} (z-Z_i)^c + \\ */
/*          && c(c-1) (x-X_i)^a (y-Y_i)^b (z-Z_i)^{c-2}. */
/*    \end{eqnarray*} */

/*    ~qmckl_ao_polynomial_vgl~ computes the values, gradients and */
/*    Laplacians at a given point in space, of all polynomials with an */
/*    angular momentum up to ~lmax~. */

/*    #+NAME: qmckl_ao_polynomial_vgl_args */
/*    | Variable  | Type              | In/Out | Description                                          | */
/*    |-----------+-------------------+--------+------------------------------------------------------| */
/*    | ~context~ | ~qmckl_context~   | in     | Global state                                         | */
/*    | ~X~       | ~double[3]~       | in     | Array containing the coordinates of the points       | */
/*    | ~R~       | ~double[3]~       | in     | Array containing the x,y,z coordinates of the center | */
/*    | ~lmax~    | ~int32_t~         | in     | Maximum angular momentum                             | */
/*    | ~n~       | ~int64_t~         | inout  | Number of computed polynomials                       | */
/*    | ~L~       | ~int32_t[n][ldl]~ | out    | Contains a,b,c for all ~n~ results                   | */
/*    | ~ldl~     | ~int64_t~         | in     | Leading dimension of ~L~                             | */
/*    | ~VGL~     | ~double[n][ldv]~  | out    | Value, gradients and Laplacian of the polynomials    | */
/*    | ~ldv~     | ~int64_t~         | in     | Leading dimension of array ~VGL~                     | */
/*    |-----------+-------------------+--------+------------------------------------------------------| */

/*    Requirements: */

/*     - ~context~ \ne ~QMCKL_NULL_CONTEXT~ */
/*     - ~n~ > 0 */
/*     - ~lmax~ >= 0 */
/*     - ~ldl~ >= 3 */
/*     - ~ldv~ >= 5 */
/*     - ~X~ is allocated with at least $3 \times 8$ bytes */
/*     - ~R~ is allocated with at least $3 \times 8$ bytes */
/*     - ~n~ >= ~(lmax+1)(lmax+2)(lmax+3)/6~ */
/*     - ~L~ is allocated with at least $3 \times n \times 4$ bytes */
/*     - ~VGL~ is allocated with at least $5 \times n \times 8$ bytes */
/*     - On output, ~n~ should be equal to ~(lmax+1)(lmax+2)(lmax+3)/6~ */
/*     - On output, the powers are given in the following order (l=a+b+c): */
/*       - Increasing values of ~l~ */
/*       - Within a given value of ~l~, alphabetical order of the */
/*         string made by a*"x" + b*"y" + c*"z" (in Python notation). */
/*         For example, with a=0, b=2 and c=1 the string is "yyz" */

/*     #+CALL: generate_c_header(table=qmckl_ao_polynomial_vgl_args,rettyp=get_value("CRetType"),fname=get_value("Name")) */

/*     #+RESULTS: */

qmckl_exit_code qmckl_ao_polynomial_vgl (
      const qmckl_context context,
      const double* X,
      const double* R,
      const int32_t lmax,
      int64_t* n,
      int32_t* const L,
      const int64_t ldl,
      double* const VGL,
      const int64_t ldv );



/* #+CALL: generate_c_header(table=qmckl_ao_polynomial_vgl_args,rettyp=get_value("CRetType"),fname="qmckl_ao_polynomial_vgl_doc") */

/* #+RESULTS: */

qmckl_exit_code qmckl_ao_polynomial_vgl_doc (
      const qmckl_context context,
      const double* X,
      const double* R,
      const int32_t lmax,
      int64_t* n,
      int32_t* const L,
      const int64_t ldl,
      double* const VGL,
      const int64_t ldv );



/* #+CALL: generate_c_header(table=qmckl_ao_polynomial_vgl_args,rettyp=get_value("CRetType"),fname="qmckl_ao_polynomial_transp_vgl") */

/* #+RESULTS: */

qmckl_exit_code qmckl_ao_polynomial_transp_vgl (
      const qmckl_context context,
      const double* X,
      const double* R,
      const int32_t lmax,
      int64_t* n,
      int32_t* const L,
      const int64_t ldl,
      double* const VGL,
      const int64_t ldv );



/* #+CALL: generate_c_header(table=qmckl_ao_polynomial_vgl_args,rettyp=get_value("CRetType"),fname="qmckl_ao_polynomial_transp_vgl_doc") */

/* #+RESULTS: */

qmckl_exit_code qmckl_ao_polynomial_transp_vgl_doc (
      const qmckl_context context,
      const double* X,
      const double* R,
      const int32_t lmax,
      int64_t* n,
      int32_t* const L,
      const int64_t ldl,
      double* const VGL,
      const int64_t ldv );



/* #   #+CALL: generate_c_header(table=qmckl_ao_polynomial_vgl_args,rettyp=get_value("CRetType"),fname="qmckl_ao_polynomial_transp_vgl_hpc") */

/*     #+RESULTS: */

qmckl_exit_code qmckl_ao_polynomial_transp_vgl_hpc (
      const qmckl_context context,
      const double* X,
      const double* R,
      const int32_t lmax,
      int64_t* n,
      int32_t* const L,
      const int64_t ldl,
      double* const VGL,
      const int64_t ldv );

/* Hessian */

/* Compute the Hessian of the polynomial part of the atomic orbitals. Similarly to the shells, the ~hessian[:][3][:]~ component is reserved for the derivative of the Laplacian. */



/*   :PROPERTIES: */
/*    :Name:     qmckl_compute_ao_polynomial_hessian */
/*    :CRetType: qmckl_exit_code */
/*    :FRetType: qmckl_exit_code */
/*    :END: */

/*    #+NAME: qmckl_compute_ao_polynomial_hessian_args */
/*    | Variable  | Type                | In/Out | Description                                          | */
/*    |-----------+---------------------+--------+------------------------------------------------------| */
/*    | ~context~ | ~qmckl_context~     | in     | Global state                                         | */
/*    | ~X~       | ~double[3]~         | in     | Array containing the coordinates of the points       | */
/*    | ~R~       | ~double[3]~         | in     | Array containing the x,y,z coordinates of the center | */
/*    | ~lmax~    | ~int32_t~           | in     | Maximum angular momentum                             | */
/*    | ~n~       | ~int64_t~           | inout  | Number of computed polynomials                       | */
/*    | ~hessian~ | ~double[ldv][4][3]~ | out    | Hessian of the polynomials                           | */
/*    |-----------+---------------------+--------+------------------------------------------------------| */


qmckl_exit_code qmckl_compute_ao_polynomial_hessian_doc (
      const qmckl_context context,
      const double* X,
      const double* R,
      const int32_t lmax,
      int64_t* n,
      double* const hessian);

qmckl_exit_code qmckl_ao_polynomial_hessian (
      const qmckl_context context,
      const double* X,
      const double* R,
      const int32_t lmax,
      int64_t* n,
      double* const hessian );

qmckl_exit_code qmckl_compute_ao_polynomial_hessian_doc (
      const qmckl_context context,
      const double* X,
      const double* R,
      const int32_t lmax,
      int64_t* n,
      double* const hessian);

/* Get */


qmckl_exit_code
qmckl_get_ao_basis_ao_hessian(qmckl_context context,
                          double* const ao_hessian,
                          const int64_t size_max);
