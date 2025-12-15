bool qmckl_electron_provided (const qmckl_context context);

/* Initialization functions */

/*    To set the data relative to the electrons in the context, the */
/*    following functions need to be called. When the data structure is */
/*    initialized, the internal ~coord_new~ and ~coord_old~ arrays are */
/*    both not allocated. */


qmckl_exit_code qmckl_set_electron_num      (qmckl_context context, const int64_t up_num, const int64_t down_num);
qmckl_exit_code qmckl_set_electron_coord    (qmckl_context context, const char transp, const int64_t walk_num, const double* coord, const int64_t size_max);

/* Number of electrons */


qmckl_exit_code qmckl_get_electron_num        (const qmckl_context context, int64_t* const num);
qmckl_exit_code qmckl_get_electron_up_num     (const qmckl_context context, int64_t* const up_num);
qmckl_exit_code qmckl_get_electron_down_num   (const qmckl_context context, int64_t* const down_num);

/* Number of walkers */

/*     A walker is a set of electron coordinates that are arguments of */
/*     the wave function. ~walk_num~ is the number of walkers. */


qmckl_exit_code qmckl_get_electron_walk_num   (const qmckl_context context, int64_t* const walk_num);

/* Electron coordinates */

/*     Returns the current electron coordinates. The pointer is assumed */
/*     to point on a memory block of size ~size_max~ \ge ~3 * elec_num * walker.num~. */
/*     The order of the indices is: */

/*     |         | Normal                     | Transposed                 | */
/*     |---------+----------------------------+----------------------------| */
/*     | C       | ~[walker.num*elec_num][3]~ | ~[3][walker.num*elec_num]~ | */
/*     | Fortran | ~(3,walker.num*elec_num)~  | ~(walker.num*elec_num, 3)~ | */



qmckl_exit_code
qmckl_get_electron_coord (const qmckl_context context,
                          const char transp,
                          double* const coord,
                          const int64_t size_max);

/* Get */


qmckl_exit_code
qmckl_get_electron_ee_distance(qmckl_context context,
                               double* const distance,
                               const int64_t size_max);

/* Get */


qmckl_exit_code qmckl_get_electron_ee_potential(qmckl_context context, double* const ee_potential);

/* Get */


qmckl_exit_code
qmckl_get_electron_en_distance(qmckl_context context,
                               double* const distance,
                               const int64_t size_max);

/* Get */


qmckl_exit_code qmckl_get_electron_en_potential(qmckl_context context, double* const en_potential);

qmckl_exit_code qmckl_compute_en_potential (
      const qmckl_context context,
      const int64_t elec_num,
      const int64_t nucl_num,
      const int64_t walk_num,
      const double* charge,
      const double* en_distance,
      double* const en_potential );
