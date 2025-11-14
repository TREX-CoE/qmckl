! Compute kinetic enregy
!    :PROPERTIES:
!    :Name:     qmckl_compute_kinetic_energy
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!    #+NAME: qmckl_compute_kinetic_energy_args
!    | ~qmckl_context~ | ~context~                                                             | in  | Global state                              |
!    | ~int64_t~       | ~walk_num~                                                            | in  | Number of walkers                         |
!    | ~int64_t~       | ~det_num_alpha~                                                       | in  | Number of determinants                    |
!    | ~int64_t~       | ~det_num_beta~                                                        | in  | Number of determinants                    |
!    | ~int64_t~       | ~alpha_num~                                                           | in  | Number of electrons                       |
!    | ~int64_t~       | ~beta_num~                                                            | in  | Number of electrons                       |
!    | ~int64_t~       | ~elec_num~                                                            | in  | Number of electrons                       |
!    | ~int64_t~       | ~mo_index_alpha[det_num_alpha][walk_num][alpha_num]~                  | in  | MO indices for electrons                  |
!    | ~int64_t~       | ~mo_index_beta[det_num_beta][walk_num][beta_num]~                     | in  | MO indices for electrons                  |
!    | ~int64_t~       | ~mo_num~                                                              | in  | Number of MOs                             |
!    | ~double~        | ~mo_vgl[5][elec_num][mo_num]~                                         | in  | Value, gradients and Laplacian of the MOs |
!    | ~double~        | ~det_value_alpha[det_num_alpha][walk_num]~                            | in  | Det of wavefunction                       |
!    | ~double~        | ~det_value_beta[det_num_beta][walk_num]~                              | in  | Det of wavefunction                       |
!    | ~double~        | ~det_inv_matrix_alpha[det_num_alpha][walk_num][alpha_num][alpha_num]~ | in  | Value, gradients and Laplacian of the Det |
!    | ~double~        | ~det_inv_matrix_beta[det_num_beta][walk_num][beta_num][beta_num]~     | in  | Value, gradients and Laplacian of the Det |
!    | ~double~        | ~e_kin[walk_num]~                                                     | out | Kinetic energy                            |


integer function qmckl_compute_kinetic_energy_f(context, walk_num, &
     det_num_alpha, det_num_beta, alpha_num, beta_num, elec_num, mo_index_alpha, mo_index_beta, &
     mo_num, mo_vgl, det_value_alpha, det_value_beta, det_inv_matrix_alpha, det_inv_matrix_beta, e_kin) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context)  , intent(in)  :: context
  integer*8, intent(in)             :: walk_num
  integer*8, intent(in)             :: det_num_alpha
  integer*8, intent(in)             :: det_num_beta
  integer*8, intent(in)             :: alpha_num
  integer*8, intent(in)             :: beta_num
  integer*8, intent(in)             :: elec_num
  integer*8, intent(in)             :: mo_num
  integer*8, intent(in)             :: mo_index_alpha(alpha_num, walk_num, det_num_alpha)
  integer*8, intent(in)             :: mo_index_beta(beta_num, walk_num, det_num_beta)
  double precision, intent(in)      :: mo_vgl(mo_num, elec_num, 5)
  double precision, intent(in)      :: det_value_alpha(walk_num, det_num_alpha)
  double precision, intent(in)      :: det_value_beta(walk_num, det_num_beta)
  double precision, intent(in)      :: det_inv_matrix_alpha(alpha_num, alpha_num, walk_num, det_num_alpha)
  double precision, intent(in)      :: det_inv_matrix_beta(beta_num, beta_num, walk_num, det_num_beta)
  double precision, intent(inout)   :: e_kin(walk_num)
  double precision                  :: tmp_e
  integer*8 :: idet, iwalk, ielec, mo_id, imo

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (alpha_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  if (beta_num < 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  if (elec_num <= 0) then
     info = QMCKL_INVALID_ARG_5
     return
  endif

  e_kin = 0.0d0
  do idet = 1, det_num_alpha
  do iwalk = 1, walk_num
    ! Alpha part
    do imo = 1, alpha_num
    do ielec = 1, alpha_num
      mo_id = mo_index_alpha(imo, iwalk, idet)
      e_kin(iwalk) = e_kin(iwalk) - 0.5d0 * det_inv_matrix_alpha(imo, ielec, iwalk, idet) * &
                                    mo_vgl(mo_id, ielec, 5)
    end do
    end do
    ! Beta part
    do imo = 1, beta_num
    do ielec = 1, beta_num
      mo_id = mo_index_beta(imo, iwalk, idet)
      e_kin(iwalk) = e_kin(iwalk) - 0.5d0 * det_inv_matrix_beta(imo, ielec, iwalk, idet) * &
                                     mo_vgl(mo_id, alpha_num + ielec, 5)
    end do
    end do
  end do
  end do

end function qmckl_compute_kinetic_energy_f



! #+CALL: generate_c_interface(table=qmckl_compute_kinetic_energy_args,rettyp=get_value("CRetType"),fname="qmckl_compute_kinetic_energy"))

! #+RESULTS:

integer(c_int32_t) function qmckl_compute_kinetic_energy &
    (context, &
	     walk_num, &
	     det_num_alpha, &
	     det_num_beta, &
	     alpha_num, &
	     beta_num, &
	     elec_num, &
	     mo_index_alpha, &
	     mo_index_beta, &
	     mo_num, &
	     mo_vgl, &
	     det_value_alpha, &
	     det_value_beta, &
	     det_inv_matrix_alpha, &
	     det_inv_matrix_beta, &
	     e_kin) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: walk_num
  integer (c_int64_t) , intent(in)  , value :: det_num_alpha
  integer (c_int64_t) , intent(in)  , value :: det_num_beta
  integer (c_int64_t) , intent(in)  , value :: alpha_num
  integer (c_int64_t) , intent(in)  , value :: beta_num
  integer (c_int64_t) , intent(in)  , value :: elec_num
  integer (c_int64_t) , intent(in)          :: mo_index_alpha(alpha_num,walk_num,det_num_alpha)
  integer (c_int64_t) , intent(in)          :: mo_index_beta(beta_num,walk_num,det_num_beta)
  integer (c_int64_t) , intent(in)  , value :: mo_num
  real    (c_double ) , intent(in)          :: mo_vgl(mo_num,elec_num,5)
  real    (c_double ) , intent(in)          :: det_value_alpha(walk_num,det_num_alpha)
  real    (c_double ) , intent(in)          :: det_value_beta(walk_num,det_num_beta)
  real    (c_double ) , intent(in)          :: det_inv_matrix_alpha(alpha_num,alpha_num,walk_num,det_num_alpha)
  real    (c_double ) , intent(in)          :: det_inv_matrix_beta(beta_num,beta_num,walk_num,det_num_beta)
  real    (c_double ) , intent(out)         :: e_kin(walk_num)

  integer(c_int32_t), external :: qmckl_compute_kinetic_energy_f
  info = qmckl_compute_kinetic_energy_f &
	 (context, &
	     walk_num, &
	     det_num_alpha, &
	     det_num_beta, &
	     alpha_num, &
	     beta_num, &
	     elec_num, &
	     mo_index_alpha, &
	     mo_index_beta, &
	     mo_num, &
	     mo_vgl, &
	     det_value_alpha, &
	     det_value_beta, &
	     det_inv_matrix_alpha, &
	     det_inv_matrix_beta, &
	     e_kin)

end function qmckl_compute_kinetic_energy

! Compute potential enregy
!    :PROPERTIES:
!    :Name:     qmckl_compute_potential_energy
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!     #+NAME: qmckl_compute_potential_energy_args
!    | ~qmckl_context~ | ~context~                | in  | Global state        |
!    | ~int64_t~       | ~walk_num~               | in  | Number of walkers   |
!    | ~int64_t~       | ~elec_num~               | in  | Number of electrons |
!    | ~int64_t~       | ~nucl_num~               | in  | Number of MOs       |
!    | ~double~        | ~ee_potential[walk_num]~ | in  | ee potential        |
!    | ~double~        | ~en_potential[walk_num]~ | in  | en potential        |
!    | ~double~        | ~repulsion~              | in  | en potential        |
!    | ~double~        | ~e_pot[walk_num]~        | out | Potential energy    |


integer function qmckl_compute_potential_energy_f(context, walk_num, &
     elec_num, nucl_num, ee_potential, en_potential, repulsion, e_pot) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context)  , intent(in)  :: context
  integer*8, intent(in)             :: walk_num
  integer*8, intent(in)             :: elec_num
  integer*8, intent(in)             :: nucl_num
  double precision, intent(in)      :: ee_potential(walk_num)
  double precision, intent(in)      :: en_potential(walk_num)
  double precision, intent(in)      :: repulsion
  double precision, intent(inout)   :: e_pot(walk_num)
  integer*8 :: idet, iwalk, ielec, mo_id, imo

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (elec_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  do iwalk = 1, walk_num
    e_pot(iwalk) = ee_potential(iwalk) + en_potential(iwalk) + repulsion
  end do

end function qmckl_compute_potential_energy_f



! #+CALL: generate_c_interface(table=qmckl_compute_potential_energy_args,rettyp=get_value("CRetType"),fname="qmckl_compute_potential_energy"))

! #+RESULTS:

integer(c_int32_t) function qmckl_compute_potential_energy &
    (context, walk_num, elec_num, nucl_num, ee_potential, en_potential, repulsion, e_pot) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: walk_num
  integer (c_int64_t) , intent(in)  , value :: elec_num
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  real    (c_double ) , intent(in)          :: ee_potential(walk_num)
  real    (c_double ) , intent(in)          :: en_potential(walk_num)
  real    (c_double ) , intent(in)  , value :: repulsion
  real    (c_double ) , intent(out)         :: e_pot(walk_num)

  integer(c_int32_t), external :: qmckl_compute_potential_energy_f
  info = qmckl_compute_potential_energy_f &
	 (context, walk_num, elec_num, nucl_num, ee_potential, en_potential, repulsion, e_pot)

end function qmckl_compute_potential_energy

! Compute local enregy
!    :PROPERTIES:
!    :Name:     qmckl_compute_local_energy
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!     #+NAME: qmckl_compute_local_energy_args
!    | ~qmckl_context~ | ~context~           | in  | Global state      |
!    | ~int64_t~       | ~walk_num~          | in  | Number of walkers |
!    | ~double~        | ~e_kin[walk_num]~   | in  | e kinetic         |
!    | ~double~        | ~e_pot[walk_num]~   | in  | e potential       |
!    | ~double~        | ~e_local[walk_num]~ | out | local energy      |


integer function qmckl_compute_local_energy_f(context, walk_num, &
     e_kin, e_pot, e_local) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context)  , intent(in)  :: context
  integer*8, intent(in)             :: walk_num
  double precision, intent(in)      :: e_kin(walk_num)
  double precision, intent(in)      :: e_pot(walk_num)
  double precision, intent(inout)   :: e_local(walk_num)
  integer*8 :: idet, iwalk, ielec, mo_id, imo

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  e_local = 0.0d0
  do iwalk = 1, walk_num
    e_local(iwalk) = e_local(iwalk) + e_kin(iwalk) + e_pot(iwalk)
  end do

end function qmckl_compute_local_energy_f



! #+CALL: generate_c_interface(table=qmckl_compute_local_energy_args,rettyp=get_value("CRetType"),fname="qmckl_compute_local_energy"))

!  #+RESULTS:

integer(c_int32_t) function qmckl_compute_local_energy &
    (context, walk_num, e_kin, e_pot, e_local) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: walk_num
  real    (c_double ) , intent(in)          :: e_kin(walk_num)
  real    (c_double ) , intent(in)          :: e_pot(walk_num)
  real    (c_double ) , intent(out)         :: e_local(walk_num)

  integer(c_int32_t), external :: qmckl_compute_local_energy_f
  info = qmckl_compute_local_energy_f &
	 (context, walk_num, e_kin, e_pot, e_local)

end function qmckl_compute_local_energy

! Compute drift vector
!    :PROPERTIES:
!    :Name:     qmckl_compute_drift_vector
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!     #+NAME: qmckl_compute_drift_vector_args
!    | ~qmckl_context~ | ~context~                                                             | in  | Global state                              |
!    | ~int64_t~       | ~walk_num~                                                            | in  | Number of walkers                         |
!    | ~int64_t~       | ~det_num_alpha~                                                       | in  | Number of determinants                    |
!    | ~int64_t~       | ~det_num_beta~                                                        | in  | Number of determinants                    |
!    | ~int64_t~       | ~alpha_num~                                                           | in  | Number of electrons                       |
!    | ~int64_t~       | ~beta_num~                                                            | in  | Number of electrons                       |
!    | ~int64_t~       | ~elec_num~                                                            | in  | Number of electrons                       |
!    | ~int64_t~       | ~mo_index_alpha[det_num_alpha][walk_num][alpha_num]~                  | in  | MO indices for electrons                  |
!    | ~int64_t~       | ~mo_index_beta[det_num_beta][walk_num][beta_num]~                     | in  | MO indices for electrons                  |
!    | ~int64_t~       | ~mo_num~                                                              | in  | Number of MOs                             |
!    | ~double~        | ~mo_vgl[5][elec_num][mo_num]~                                         | in  | Value, gradients and Laplacian of the MOs |
!    | ~double~        | ~det_inv_matrix_alpha[det_num_alpha][walk_num][alpha_num][alpha_num]~ | in  | Value, gradients and Laplacian of the Det |
!    | ~double~        | ~det_inv_matrix_beta[det_num_beta][walk_num][beta_num][beta_num]~     | in  | Value, gradients and Laplacian of the Det |
!    | ~double~        | ~r_drift[walk_num][elec_num][3]~                                      | out | Kinetic energy                            |


integer function qmckl_compute_drift_vector_f(context, walk_num, &
     det_num_alpha, det_num_beta, alpha_num, beta_num, elec_num, mo_index_alpha, mo_index_beta, &
     mo_num, mo_vgl, det_inv_matrix_alpha, det_inv_matrix_beta, r_drift) &
     result(info)
  use qmckl
  implicit none
  integer(qmckl_context)  , intent(in)  :: context
  integer*8, intent(in)             :: walk_num
  integer*8, intent(in)             :: det_num_alpha
  integer*8, intent(in)             :: det_num_beta
  integer*8, intent(in)             :: alpha_num
  integer*8, intent(in)             :: beta_num
  integer*8, intent(in)             :: elec_num
  integer*8, intent(in)             :: mo_num
  integer*8, intent(in)             :: mo_index_alpha(alpha_num, walk_num, det_num_alpha)
  integer*8, intent(in)             :: mo_index_beta(beta_num, walk_num, det_num_beta)
  double precision, intent(in)      :: mo_vgl(mo_num, elec_num, 5)
  double precision, intent(in)      :: det_inv_matrix_alpha(alpha_num, alpha_num, walk_num, det_num_alpha)
  double precision, intent(in)      :: det_inv_matrix_beta(beta_num, beta_num, walk_num, det_num_beta)
  double precision, intent(inout)   :: r_drift(3,elec_num,walk_num)
  integer*8 :: idet, iwalk, ielec, mo_id, imo

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (alpha_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  if (beta_num < 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  if (elec_num <= 0) then
     info = QMCKL_INVALID_ARG_5
     return
  endif

  r_drift = 0.0d0
  do idet = 1, det_num_alpha
  do iwalk = 1, walk_num
    ! Alpha part
    do imo = 1, alpha_num
    do ielec = 1, alpha_num
      mo_id = mo_index_alpha(imo, iwalk, idet)
      r_drift(1,ielec,iwalk) = r_drift(1,ielec,iwalk) + 2.0d0 * det_inv_matrix_alpha(imo, ielec, iwalk, idet) * &
                                    mo_vgl(mo_id, ielec, 2)
      r_drift(2,ielec,iwalk) = r_drift(2,ielec,iwalk) + 2.0d0 * det_inv_matrix_alpha(imo, ielec, iwalk, idet) * &
                                    mo_vgl(mo_id, ielec, 3)
      r_drift(3,ielec,iwalk) = r_drift(3,ielec,iwalk) + 2.0d0 * det_inv_matrix_alpha(imo, ielec, iwalk, idet) * &
                                    mo_vgl(mo_id, ielec, 4)
    end do
    end do
    ! Beta part
    do imo = 1, beta_num
    do ielec = 1, beta_num
      mo_id = mo_index_beta(imo, iwalk, idet)
      r_drift(1,alpha_num + ielec,iwalk) = r_drift(1,alpha_num + ielec,iwalk) + &
                                    2.0d0 * det_inv_matrix_beta(imo, ielec, iwalk, idet) * &
                                    mo_vgl(mo_id, alpha_num + ielec, 2)
      r_drift(2,alpha_num + ielec,iwalk) = r_drift(2,alpha_num + ielec,iwalk) + &
                                    2.0d0 * det_inv_matrix_beta(imo, ielec, iwalk, idet) * &
                                    mo_vgl(mo_id, alpha_num + ielec, 3)
      r_drift(3,alpha_num + ielec,iwalk) = r_drift(3,alpha_num + ielec,iwalk) + &
                                    2.0d0 * det_inv_matrix_beta(imo, ielec, iwalk, idet) * &
                                    mo_vgl(mo_id, alpha_num + ielec, 4)
    end do
    end do
  end do
  end do

end function qmckl_compute_drift_vector_f



! #+CALL: generate_c_interface(table=qmckl_compute_drift_vector_args,rettyp=get_value("CRetType"),fname="qmckl_compute_drift_vector"))

!  #+RESULTS:

integer(c_int32_t) function qmckl_compute_drift_vector &
    (context, &
	     walk_num, &
	     det_num_alpha, &
	     det_num_beta, &
	     alpha_num, &
	     beta_num, &
	     elec_num, &
	     mo_index_alpha, &
	     mo_index_beta, &
	     mo_num, &
	     mo_vgl, &
	     det_inv_matrix_alpha, &
	     det_inv_matrix_beta, &
	     r_drift) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: walk_num
  integer (c_int64_t) , intent(in)  , value :: det_num_alpha
  integer (c_int64_t) , intent(in)  , value :: det_num_beta
  integer (c_int64_t) , intent(in)  , value :: alpha_num
  integer (c_int64_t) , intent(in)  , value :: beta_num
  integer (c_int64_t) , intent(in)  , value :: elec_num
  integer (c_int64_t) , intent(in)          :: mo_index_alpha(alpha_num,walk_num,det_num_alpha)
  integer (c_int64_t) , intent(in)          :: mo_index_beta(beta_num,walk_num,det_num_beta)
  integer (c_int64_t) , intent(in)  , value :: mo_num
  real    (c_double ) , intent(in)          :: mo_vgl(mo_num,elec_num,5)
  real    (c_double ) , intent(in)          :: det_inv_matrix_alpha(alpha_num,alpha_num,walk_num,det_num_alpha)
  real    (c_double ) , intent(in)          :: det_inv_matrix_beta(beta_num,beta_num,walk_num,det_num_beta)
  real    (c_double ) , intent(out)         :: r_drift(3,elec_num,walk_num)

  integer(c_int32_t), external :: qmckl_compute_drift_vector_f
  info = qmckl_compute_drift_vector_f &
	 (context, &
	     walk_num, &
	     det_num_alpha, &
	     det_num_beta, &
	     alpha_num, &
	     beta_num, &
	     elec_num, &
	     mo_index_alpha, &
	     mo_index_beta, &
	     mo_num, &
	     mo_vgl, &
	     det_inv_matrix_alpha, &
	     det_inv_matrix_beta, &
	     r_drift)

end function qmckl_compute_drift_vector
