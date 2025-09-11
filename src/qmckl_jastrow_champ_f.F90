! Compute
!      :PROPERTIES:
!      :Name:     qmckl_compute_jastrow_champ_asymp_jasb
!      :CRetType: qmckl_exit_code
!      :FRetType: qmckl_exit_code
!      :END:

!      #+NAME: qmckl_asymp_jasb_args
!      |---------------------+----------------------+--------+---------------------------------------------------------------|
!      | Variable            | Type                 | In/Out | Description                                                   |
!      |---------------------+----------------------+--------+---------------------------------------------------------------|
!      | ~context~           | ~qmckl_context~      | in     | Global state                                                  |
!      | ~bord_num~          | ~int64_t~            | in     | Order of the polynomial                                       |
!      | ~b_vector~          | ~double[bord_num+1]~ | in     | Values of b                                                   |
!      | ~rescale_factor_ee~ | ~double~             | in     | Electron coordinates                                          |
!      | ~spin_independent~  | ~int32_t~            | in     | If 1, same parameters for parallel and anti-parallel pairs    |
!      | ~asymp_jasb~        | ~double[2]~          | out    | Asymptotic value                                              |
!      |---------------------+----------------------+--------+---------------------------------------------------------------|


function qmckl_compute_jastrow_champ_asymp_jasb_doc(context, &
     bord_num, b_vector, rescale_factor_ee, spin_independent, asymp_jasb) &
     bind(C) result(info)
  use qmckl
  implicit none

  integer (qmckl_context) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: bord_num
  real    (c_double ) , intent(in)          :: b_vector(bord_num+1)
  real    (c_double ) , intent(in)  , value :: rescale_factor_ee
  integer (c_int32_t) , intent(in)  , value :: spin_independent
  real    (c_double ) , intent(out)         :: asymp_jasb(2)
  integer(qmckl_exit_code)                  :: info

  integer*8 :: i, p
  double precision   :: kappa_inv, x, asym_one
  kappa_inv = 1.0d0 / rescale_factor_ee

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (bord_num < 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  asym_one = b_vector(1) * kappa_inv / (1.0d0 + b_vector(2) * kappa_inv)
  if (spin_independent == 1) then
     asymp_jasb(:) = (/asym_one, asym_one/)
  else
     asymp_jasb(:) = (/0.5d0*asym_one, asym_one/)
  end if

  x = kappa_inv
  do p = 2, bord_num
     x = x * kappa_inv
     do i = 1, 2
        asymp_jasb(i) = asymp_jasb(i) + b_vector(p + 1) * x
     end do
  end do

end function qmckl_compute_jastrow_champ_asymp_jasb_doc

! Compute
!      :PROPERTIES:
!      :Name:     qmckl_compute_ee_distance_rescaled
!      :CRetType: qmckl_exit_code
!      :FRetType: qmckl_exit_code
!      :END:

!      #+NAME: qmckl_ee_distance_rescaled_args
!      |---------------------+----------------------------------------+--------+--------------------------------------|
!      | Variable            | Type                                   | In/Out | Description                          |
!      |---------------------+----------------------------------------+--------+--------------------------------------|
!      | ~context~           | ~qmckl_context~                        | in     | Global state                         |
!      | ~elec_num~          | ~int64_t~                              | in     | Number of electrons                  |
!      | ~rescale_factor_ee~ | ~double~                               | in     | Factor to rescale ee distances       |
!      | ~walk_num~          | ~int64_t~                              | in     | Number of walkers                    |
!      | ~coord~             | ~double[3][walk_num][elec_num]~        | in     | Electron coordinates                 |
!      | ~ee_distance~       | ~double[walk_num][elec_num][elec_num]~ | out    | Electron-electron rescaled distances |
!      |---------------------+----------------------------------------+--------+--------------------------------------|


function qmckl_compute_ee_distance_rescaled_doc(context, &
     elec_num, rescale_factor_ee, walk_num, &
     coord, ee_distance_rescaled) &
     bind(C) result(info)
  use qmckl
  implicit none

  integer(qmckl_context), intent(in), value :: context
  integer (c_int64_t) , intent(in)  , value :: elec_num
  real    (c_double ) , intent(in)  , value :: rescale_factor_ee
  integer (c_int64_t) , intent(in)  , value :: walk_num
  real    (c_double ) , intent(in)          :: coord(elec_num,walk_num,3)
  real    (c_double ) , intent(out)         :: ee_distance_rescaled(elec_num,elec_num,walk_num)
  integer(qmckl_exit_code)                  :: info

  integer*8 :: k

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (elec_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  do k=1,walk_num
     info = qmckl_distance_rescaled(context, 'T', 'T', elec_num, elec_num, &
          coord(1,k,1), elec_num * walk_num, &
          coord(1,k,1), elec_num * walk_num, &
          ee_distance_rescaled(1,1,k), elec_num, rescale_factor_ee)
     if (info /= QMCKL_SUCCESS) then
        exit
     endif
  end do

end function qmckl_compute_ee_distance_rescaled_doc

! Compute
!      :PROPERTIES:
!      :Name:     qmckl_compute_ee_distance_rescaled_gl
!      :CRetType: qmckl_exit_code
!      :FRetType: qmckl_exit_code
!      :END:

!      #+NAME: qmckl_ee_distance_rescaled_gl_args
!      |---------------------------+-------------------------------------------+--------+-------------------------------------------------|
!      | Variable                  | Type                                      | In/Out | Description                                     |
!      |---------------------------+-------------------------------------------+--------+-------------------------------------------------|
!      | ~context~                 | ~qmckl_context~                           | in     | Global state                                    |
!      | ~elec_num~                | ~int64_t~                                 | in     | Number of electrons                             |
!      | ~rescale_factor_ee~       | ~double~                                  | in     | Factor to rescale ee distances                  |
!      | ~walk_num~                | ~int64_t~                                 | in     | Number of walkers                               |
!      | ~coord~                   | ~double[3][walk_num][elec_num]~           | in     | Electron coordinates                            |
!      | ~ee_distance_rescaled_gl~ | ~double[walk_num][elec_num][elec_num][4]~ | out    | Electron-electron rescaled distance derivatives |
!      |---------------------------+-------------------------------------------+--------+-------------------------------------------------|


function qmckl_compute_ee_distance_rescaled_gl_doc(context,  &
     elec_num, rescale_factor_ee, walk_num, coord, ee_distance_rescaled_gl) &
     bind(C) result(info)
  use qmckl
  implicit none

  integer(qmckl_context), intent(in), value :: context
  integer (c_int64_t) , intent(in)  , value :: elec_num
  real    (c_double ) , intent(in)  , value :: rescale_factor_ee
  integer (c_int64_t) , intent(in)  , value :: walk_num
  real    (c_double ) , intent(in)          :: coord(elec_num,walk_num,3)
  real    (c_double ) , intent(out)         :: ee_distance_rescaled_gl(4,elec_num,elec_num,walk_num)
  integer(qmckl_exit_code)                  :: info

  integer*8 :: k

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (elec_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  do k=1,walk_num
     info = qmckl_distance_rescaled_gl(context, 'T', 'T', elec_num, elec_num, &
          coord(1,k,1), elec_num*walk_num, &
          coord(1,k,1), elec_num*walk_num, &
          ee_distance_rescaled_gl(1,1,1,k), elec_num, rescale_factor_ee)
     if (info /= QMCKL_SUCCESS) then
        exit
     endif
  end do

end function qmckl_compute_ee_distance_rescaled_gl_doc

! Compute
!      :PROPERTIES:
!      :Name:     qmckl_compute_jastrow_champ_factor_ee
!      :CRetType: qmckl_exit_code
!      :FRetType: qmckl_exit_code
!      :END:

!      #+NAME: qmckl_factor_ee_args
!      |------------------------+----------------------------------------+--------+---------------------------------|
!      | Variable               | Type                                   | In/Out | Description                     |
!      |------------------------+----------------------------------------+--------+---------------------------------|
!      | ~context~              | ~qmckl_context~                        | in     | Global state                    |
!      | ~walk_num~             | ~int64_t~                              | in     | Number of walkers               |
!      | ~elec_num~             | ~int64_t~                              | in     | Number of electrons             |
!      | ~up_num~               | ~int64_t~                              | in     | Number of alpha electrons       |
!      | ~bord_num~             | ~int64_t~                              | in     | Number of coefficients          |
!      | ~b_vector~             | ~double[bord_num+1]~                   | in     | List of coefficients            |
!      | ~ee_distance_rescaled~ | ~double[walk_num][elec_num][elec_num]~ | in     | Electron-electron distances     |
!      | ~asymp_jasb~           | ~double[2]~                            | in     | Asymptotic value of the Jastrow |
!      | ~factor_ee~            | ~double[walk_num]~                     | out    | $f_{ee}$                        |
!      |------------------------+----------------------------------------+--------+---------------------------------|


function qmckl_compute_jastrow_champ_factor_ee_doc(context, &
     walk_num, elec_num, up_num, bord_num, b_vector, &
     ee_distance_rescaled, asymp_jasb, spin_independent, factor_ee) &
     bind(C) result(info)
  use qmckl
  implicit none

  integer (qmckl_context), intent(in), value :: context
  integer (c_int64_t)    , intent(in), value :: walk_num
  integer (c_int64_t)    , intent(in), value :: elec_num
  integer (c_int64_t)    , intent(in), value :: up_num
  integer (c_int64_t)    , intent(in), value :: bord_num
  real    (c_double )    , intent(in)        :: b_vector(bord_num+1)
  real    (c_double )    , intent(in)        :: ee_distance_rescaled(elec_num,elec_num,walk_num)
  real    (c_double )    , intent(in)        :: asymp_jasb(2)
  integer (c_int32_t)    , intent(in), value :: spin_independent
  real    (c_double )    , intent(out)       :: factor_ee(walk_num)
  integer(qmckl_exit_code)                   :: info

  integer*8 :: i, j, k, nw
  double precision   :: x, xk

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

  if (bord_num < 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif



  do nw =1, walk_num

     factor_ee(nw) = 0.0d0
     do j=1,elec_num
        do i=1,j-1
           x = ee_distance_rescaled(i,j,nw)
           if (spin_independent == 1) then
              factor_ee(nw) = factor_ee(nw) + b_vector(1) * x / (1.d0 + b_vector(2) * x) - asymp_jasb(2)
           else
              if ( (j <= up_num).or.(i > up_num) ) then
                 factor_ee(nw) = factor_ee(nw) + 0.5d0 * b_vector(1) * x / (1.d0 + b_vector(2) * x) - asymp_jasb(1)
              else
                 factor_ee(nw) = factor_ee(nw) + b_vector(1) * x / (1.d0 + b_vector(2) * x) - asymp_jasb(2)
              endif
           endif

           xk = x
           do k=2,bord_num
              xk = xk * x
              factor_ee(nw) = factor_ee(nw) + b_vector(k+1)* xk
           end do
        end do
     end do

  end do

end function qmckl_compute_jastrow_champ_factor_ee_doc

! Compute
!      :PROPERTIES:
!      :Name:     qmckl_compute_jastrow_champ_factor_ee_gl
!      :CRetType: qmckl_exit_code
!      :FRetType: qmckl_exit_code
!      :END:

!      #+NAME: qmckl_factor_ee_gl_args
!      |---------------------------+-------------------------------------------+--------+-----------------------------------------------------------|
!      | Variable                  | Type                                      | In/Out | Description                                               |
!      |---------------------------+-------------------------------------------+--------+-----------------------------------------------------------|
!      | ~context~                 | ~qmckl_context~                           | in     | Global state                                              |
!      | ~walk_num~                | ~int64_t~                                 | in     | Number of walkers                                         |
!      | ~elec_num~                | ~int64_t~                                 | in     | Number of electrons                                       |
!      | ~up_num~                  | ~int64_t~                                 | in     | Number of alpha electrons                                 |
!      | ~bord_num~                | ~int64_t~                                 | in     | Number of coefficients                                    |
!      | ~b_vector~                | ~double[bord_num+1]~                      | in     | List of coefficients                                      |
!      | ~ee_distance_rescaled~    | ~double[walk_num][elec_num][elec_num]~    | in     | Electron-electron distances                               |
!      | ~ee_distance_rescaled_gl~ | ~double[walk_num][elec_num][elec_num][4]~ | in     | Electron-electron distances                               |
!      | ~spin_independent~        | ~int32_t~                                 | in     | If 1, same parameters for parallel and antiparallel spins |
!      | ~factor_ee_gl~            | ~double[walk_num][4][elec_num]~           | out    | Electron-electron distances                               |
!      |---------------------------+-------------------------------------------+--------+-----------------------------------------------------------|


function qmckl_compute_jastrow_champ_factor_ee_gl_doc( &
     context, walk_num, elec_num, up_num, bord_num, &
     b_vector, ee_distance_rescaled, ee_distance_rescaled_gl,  &
     spin_independent, factor_ee_gl) &
     bind(C) result(info)
  use qmckl
  implicit none

  integer (qmckl_context), intent(in), value :: context
  integer (c_int64_t) , intent(in)  , value :: walk_num
  integer (c_int64_t) , intent(in)  , value :: elec_num
  integer (c_int64_t) , intent(in)  , value :: up_num
  integer (c_int64_t) , intent(in)  , value :: bord_num
  real    (c_double ) , intent(in)          :: b_vector(bord_num+1)
  real    (c_double ) , intent(in)          :: ee_distance_rescaled(elec_num,elec_num,walk_num)
  real    (c_double ) , intent(in)          :: ee_distance_rescaled_gl(4,elec_num,elec_num,walk_num)
  integer (c_int32_t) , intent(in)  , value :: spin_independent
  real    (c_double ) , intent(out)         :: factor_ee_gl(elec_num,4,walk_num)
  integer(qmckl_exit_code)                  :: info

  integer*8 :: i, j, k, nw, ii
  double precision   :: x, x1, kf
  double precision   :: denom, invdenom, invdenom2, f
  double precision   :: grad_c2
  double precision   :: dx(4)

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

  if (bord_num < 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  if ((spin_independent < 0).or.(spin_independent > 1)) then
     info = QMCKL_INVALID_ARG_8
     return
  endif

  do nw =1, walk_num
     factor_ee_gl(:,:,nw) = 0.0d0

     do j = 1, elec_num
        do i = 1, elec_num
           if (i == j) cycle

           x = ee_distance_rescaled(i,j,nw)

           denom         = 1.0d0 + b_vector(2) * x
           invdenom      = 1.0d0 / denom
           invdenom2     = invdenom * invdenom

           dx(1) = ee_distance_rescaled_gl(1, i, j, nw)
           dx(2) = ee_distance_rescaled_gl(2, i, j, nw)
           dx(3) = ee_distance_rescaled_gl(3, i, j, nw)
           dx(4) = ee_distance_rescaled_gl(4, i, j, nw)

           grad_c2 = dx(1)*dx(1) + dx(2)*dx(2) + dx(3)*dx(3)

           if (spin_independent == 1) then
              f = b_vector(1) * invdenom2
           else
              if((i <= up_num .and. j <= up_num ) .or. (i >  up_num .and. j >  up_num)) then
                 f = 0.5d0 * b_vector(1) * invdenom2
              else
                 f = b_vector(1) * invdenom2
              end if
           end if

           factor_ee_gl(i,1,nw) = factor_ee_gl(i,1,nw) + f * dx(1)
           factor_ee_gl(i,2,nw) = factor_ee_gl(i,2,nw) + f * dx(2)
           factor_ee_gl(i,3,nw) = factor_ee_gl(i,3,nw) + f * dx(3)
           factor_ee_gl(i,4,nw) = factor_ee_gl(i,4,nw) &
                + f * (dx(4) - 2.d0 * b_vector(2) * grad_c2 * invdenom)


           kf = 2.d0
           x1 = x
           x = 1.d0
           do k=2, bord_num
              f = b_vector(k+1) * kf * x
              factor_ee_gl(i,1,nw) = factor_ee_gl(i,1,nw) + f * x1 * dx(1)
              factor_ee_gl(i,2,nw) = factor_ee_gl(i,2,nw) + f * x1 * dx(2)
              factor_ee_gl(i,3,nw) = factor_ee_gl(i,3,nw) + f * x1 * dx(3)
              factor_ee_gl(i,4,nw) = factor_ee_gl(i,4,nw) &
                   + f * (x1 * dx(4) + (kf-1.d0) * grad_c2)
              x = x*x1
              kf = kf + 1.d0
           end do

        end do
     end do

  end do

end function qmckl_compute_jastrow_champ_factor_ee_gl_doc

! Compute
!      :PROPERTIES:
!      :Name:     qmckl_compute_jastrow_champ_asymp_jasa
!      :CRetType: qmckl_exit_code
!      :FRetType: qmckl_exit_code
!      :END:

!      #+NAME: qmckl_asymp_jasa_args
!      |---------------------+-------------------------------------+--------+----------------------------|
!      | Variable            | Type                                | In/Out | Description                |
!      |---------------------+-------------------------------------+--------+----------------------------|
!      | ~context~           | ~qmckl_context~                     | in     | Global state               |
!      | ~aord_num~          | ~int64_t~                           | in     | Order of the polynomial    |
!      | ~type_nucl_num~     | ~int64_t~                           | in     | Number of nucleus types    |
!      | ~a_vector~          | ~double[type_nucl_num][aord_num+1]~ | in     | Values of a                |
!      | ~rescale_factor_en~ | ~double[type_nucl_num]~             | in     | Electron nucleus distances |
!      | ~asymp_jasa~        | ~double[type_nucl_num]~             | out    | Asymptotic value           |
!      |---------------------+-------------------------------------+--------+----------------------------|


function qmckl_compute_jastrow_champ_asymp_jasa(context, aord_num, type_nucl_num, a_vector, &
     rescale_factor_en, asymp_jasa) &
     result(info) bind(C)
  use, intrinsic :: iso_c_binding
  use qmckl
  implicit none
  integer(qmckl_context), intent(in), value :: context
  integer (c_int64_t) , intent(in)  , value :: aord_num
  integer (c_int64_t) , intent(in)  , value :: type_nucl_num
  real    (c_double ) , intent(in)          :: a_vector(aord_num+1,type_nucl_num)
  real    (c_double ) , intent(in)          :: rescale_factor_en(type_nucl_num)
  real    (c_double ) , intent(out)         :: asymp_jasa(type_nucl_num)
  integer(qmckl_exit_code)                  :: info

  integer*8 :: i, j, p
  double precision   :: kappa_inv, x, asym_one

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (aord_num < 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  do i=1,type_nucl_num

     kappa_inv = 1.0d0 / rescale_factor_en(i)

     asymp_jasa(i) = a_vector(1,i) * kappa_inv / (1.0d0 + a_vector(2,i) * kappa_inv)

     x = kappa_inv
     do p = 2, aord_num
        x = x * kappa_inv
        asymp_jasa(i) = asymp_jasa(i) + a_vector(p+1, i) * x
     end do

  end do

end function qmckl_compute_jastrow_champ_asymp_jasa

! Compute
!      :PROPERTIES:
!      :Name:     qmckl_compute_en_distance_rescaled
!      :CRetType: qmckl_exit_code
!      :FRetType: qmckl_exit_code
!      :END:

!     #+NAME: qmckl_en_distance_rescaled_args
!      |------------------------+----------------------------------------+--------+-----------------------------------|
!      | Variable               | Type                                   | In/Out | Description                       |
!      |------------------------+----------------------------------------+--------+-----------------------------------|
!      | ~context~              | ~qmckl_context~                        | in     | Global state                      |
!      | ~elec_num~             | ~int64_t~                              | in     | Number of electrons               |
!      | ~nucl_num~             | ~int64_t~                              | in     | Number of nuclei                  |
!      | ~type_nucl_num~        | ~int64_t~                              | in     | Number of types of nuclei         |
!      | ~type_nucl_vector~     | ~int64_t[nucl_num]~                    | in     | Number of types of nuclei         |
!      | ~rescale_factor_en~    | ~double[type_nucl_num]~                | in     | The factor for rescaled distances |
!      | ~walk_num~             | ~int64_t~                              | in     | Number of walkers                 |
!      | ~elec_coord~           | ~double[3][walk_num][elec_num]~        | in     | Electron coordinates              |
!      | ~nucl_coord~           | ~double[3][nucl_num]~                  | in     | Nuclear coordinates               |
!      | ~en_distance_rescaled~ | ~double[walk_num][nucl_num][elec_num]~ | out    | Electron-nucleus distances        |
!      |------------------------+----------------------------------------+--------+-----------------------------------|


function qmckl_compute_en_distance_rescaled_doc(context, &
     elec_num, nucl_num, type_nucl_num, &
     type_nucl_vector, rescale_factor_en, walk_num, elec_coord, &
     nucl_coord, en_distance_rescaled) &
     bind(C) result(info)
  use qmckl
  implicit none
  integer (qmckl_context), intent(in), value :: context
  integer (c_int64_t) , intent(in)  , value :: elec_num
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  integer (c_int64_t) , intent(in)  , value :: type_nucl_num
  integer (c_int64_t) , intent(in)          :: type_nucl_vector(nucl_num)
  real    (c_double ) , intent(in)          :: rescale_factor_en(type_nucl_num)
  integer (c_int64_t) , intent(in)  , value :: walk_num
  real    (c_double ) , intent(in)          :: elec_coord(elec_num,walk_num,3)
  real    (c_double ) , intent(in)          :: nucl_coord(nucl_num,3)
  real    (c_double ) , intent(out)         :: en_distance_rescaled(elec_num,nucl_num,walk_num)
  integer(qmckl_exit_code)                  :: info

  integer*8 :: i, k
  double precision      :: coord(3)

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (elec_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (nucl_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_5
     return
  endif

  do i=1, nucl_num
     coord(1:3) = nucl_coord(i,1:3)
     do k=1,walk_num
        info = qmckl_distance_rescaled(context, 'T', 'N', elec_num, 1_8, &
             elec_coord(1,k,1), elec_num*walk_num, coord, 3_8, &
             en_distance_rescaled(1,i,k), elec_num, rescale_factor_en(type_nucl_vector(i)+1))
        if (info /= QMCKL_SUCCESS) then
           return
        endif
     end do
  end do

end function qmckl_compute_en_distance_rescaled_doc

! Compute
!      :PROPERTIES:
!      :Name:     qmckl_compute_en_distance_rescaled_gl
!      :CRetType: qmckl_exit_code
!      :FRetType: qmckl_exit_code
!      :END:

!      #+NAME: qmckl_en_distance_rescaled_gl_args
!      |---------------------------+-------------------------------------------+--------+---------------------------------------|
!      | Variable                  | Type                                      | In/Out | Description                           |
!      |---------------------------+-------------------------------------------+--------+---------------------------------------|
!      | ~context~                 | ~qmckl_context~                           | in     | Global state                          |
!      | ~elec_num~                | ~int64_t~                                 | in     | Number of electrons                   |
!      | ~nucl_num~                | ~int64_t~                                 | in     | Number of nuclei                      |
!      | ~type_nucl_num~           | ~int64_t~                                 | in     | Number of nucleus types               |
!      | ~type_nucl_vector~        | ~int64_t[nucl_num]~                       | in     | Array of nucleus types                |
!      | ~rescale_factor_en~       | ~double[nucl_num]~                        | in     | The factors for rescaled distances    |
!      | ~walk_num~                | ~int64_t~                                 | in     | Number of walkers                     |
!      | ~elec_coord~              | ~double[3][walk_num][elec_num]~           | in     | Electron coordinates                  |
!      | ~nucl_coord~              | ~double[3][nucl_num]~                     | in     | Nuclear coordinates                   |
!      | ~en_distance_rescaled_gl~ | ~double[walk_num][nucl_num][elec_num][4]~ | out    | Electron-nucleus distance derivatives |
!      |---------------------------+-------------------------------------------+--------+---------------------------------------|


function qmckl_compute_en_distance_rescaled_gl_doc(context, elec_num, nucl_num, &
     type_nucl_num, type_nucl_vector, rescale_factor_en, walk_num, elec_coord, &
     nucl_coord, en_distance_rescaled_gl) &
     result(info) bind(C)
  use, intrinsic :: iso_c_binding
  use qmckl
  implicit none
  integer(qmckl_context), intent(in)  :: context
  integer (c_int64_t) , intent(in)  , value :: elec_num
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  integer (c_int64_t) , intent(in)  , value :: type_nucl_num
  integer (c_int64_t) , intent(in)          :: type_nucl_vector(nucl_num)
  real    (c_double ) , intent(in)          :: rescale_factor_en(nucl_num)
  integer (c_int64_t) , intent(in)  , value :: walk_num
  real    (c_double ) , intent(in)          :: elec_coord(elec_num,walk_num,3)
  real    (c_double ) , intent(in)          :: nucl_coord(nucl_num,3)
  real    (c_double ) , intent(out)         :: en_distance_rescaled_gl(4,elec_num,nucl_num,walk_num)

  integer(qmckl_exit_code)                  :: info
  integer*8 :: i, k
  double precision :: coord(3)

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (elec_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (nucl_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_5
     return
  endif

  do i=1, nucl_num
     coord(1:3) = nucl_coord(i,1:3)
     do k=1,walk_num
        info = qmckl_distance_rescaled_gl(context, 'T', 'T', elec_num, 1_8, &
             elec_coord(1,k,1), elec_num*walk_num, coord, 1_8, &
             en_distance_rescaled_gl(1,1,i,k), elec_num, rescale_factor_en(type_nucl_vector(i)+1))
        if (info /= QMCKL_SUCCESS) then
           return
        endif
     end do
  end do

end function qmckl_compute_en_distance_rescaled_gl_doc

! Compute
!      :PROPERTIES:
!      :Name:     qmckl_compute_jastrow_champ_factor_en_doc
!      :CRetType: qmckl_exit_code
!      :FRetType: qmckl_exit_code
!      :END:

!      #+NAME: qmckl_factor_en_args
!      |------------------------+----------------------------------------+--------+----------------------------|
!      | Variable               | Type                                   | In/Out | Description                |
!      |------------------------+----------------------------------------+--------+----------------------------|
!      | ~context~              | ~qmckl_context~                        | in     | Global state               |
!      | ~walk_num~             | ~int64_t~                              | in     | Number of walkers          |
!      | ~elec_num~             | ~int64_t~                              | in     | Number of electrons        |
!      | ~nucl_num~             | ~int64_t~                              | in     | Number of nuclei           |
!      | ~type_nucl_num~        | ~int64_t~                              | in     | Number of unique nuclei    |
!      | ~type_nucl_vector~     | ~int64_t[nucl_num]~                    | in     | IDs of unique nuclei       |
!      | ~aord_num~             | ~int64_t~                              | in     | Number of coefficients     |
!      | ~a_vector~             | ~double[type_nucl_num][aord_num+1]~    | in     | List of coefficients       |
!      | ~en_distance_rescaled~ | ~double[walk_num][nucl_num][elec_num]~ | in     | Electron-nucleus distances |
!      | ~asymp_jasa~           | ~double[type_nucl_num]~                | in     | Type of nuclei             |
!      | ~factor_en~            | ~double[walk_num]~                     | out    | Electron-nucleus jastrow   |
!      |------------------------+----------------------------------------+--------+----------------------------|


function qmckl_compute_jastrow_champ_factor_en_doc( &
     context, walk_num, elec_num, nucl_num, type_nucl_num, &
     type_nucl_vector, aord_num, a_vector, &
     en_distance_rescaled, asymp_jasa, factor_en) &
     bind(C) result(info)
  use qmckl
  implicit none

  integer (qmckl_context), intent(in), value :: context
  integer (c_int64_t) , intent(in)  , value :: walk_num
  integer (c_int64_t) , intent(in)  , value :: elec_num
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  integer (c_int64_t) , intent(in)  , value :: type_nucl_num
  integer (c_int64_t) , intent(in)          :: type_nucl_vector(nucl_num)
  integer (c_int64_t) , intent(in)  , value :: aord_num
  real    (c_double ) , intent(in)          :: a_vector(aord_num+1,type_nucl_num)
  real    (c_double ) , intent(in)          :: en_distance_rescaled(elec_num,nucl_num,walk_num)
  real    (c_double ) , intent(in)          :: asymp_jasa(type_nucl_num)
  real    (c_double ) , intent(out)         :: factor_en(walk_num)
  integer(qmckl_exit_code)                  :: info

  integer*8 :: i, a, p, nw
  double precision   :: x, power_ser

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

  if (nucl_num <= 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  if (type_nucl_num <= 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  if (aord_num < 0) then
     info = QMCKL_INVALID_ARG_7
     return
  endif


  do nw =1, walk_num
     factor_en(nw) = 0.0d0
     do a = 1, nucl_num
        do i = 1, elec_num
           x = en_distance_rescaled(i, a, nw)

           factor_en(nw) = factor_en(nw) + a_vector(1, type_nucl_vector(a)+1) * x / &
                (1.0d0 + a_vector(2, type_nucl_vector(a)+1) * x) - asymp_jasa(type_nucl_vector(a)+1)

           do p = 2, aord_num
              x = x * en_distance_rescaled(i, a, nw)
              factor_en(nw) = factor_en(nw) + a_vector(p + 1, type_nucl_vector(a)+1) * x
           end do

        end do
     end do
  end do

end function qmckl_compute_jastrow_champ_factor_en_doc

! Compute
!      :PROPERTIES:
!      :Name:     qmckl_compute_jastrow_champ_factor_en_gl
!      :CRetType: qmckl_exit_code
!      :FRetType: qmckl_exit_code
!      :END:

!      #+NAME: qmckl_factor_en_gl_args
!      |---------------------------+-------------------------------------------+--------+---------------------------------------|
!      | Variable                  | Type                                      | In/Out | Description                           |
!      |---------------------------+-------------------------------------------+--------+---------------------------------------|
!      | ~context~                 | ~qmckl_context~                           | in     | Global state                          |
!      | ~walk_num~                | ~int64_t~                                 | in     | Number of walkers                     |
!      | ~elec_num~                | ~int64_t~                                 | in     | Number of electrons                   |
!      | ~nucl_num~                | ~int64_t~                                 | in     | Number of nuclei                      |
!      | ~type_nucl_num~           | ~int64_t~                                 | in     | Number of unique nuclei               |
!      | ~type_nucl_vector~        | ~int64_t[nucl_num]~                       | in     | IDs of unique nuclei                  |
!      | ~aord_num~                | ~int64_t~                                 | in     | Number of coefficients                |
!      | ~a_vector~                | ~double[type_nucl_num][aord_num+1]~       | in     | List of coefficients                  |
!      | ~en_distance_rescaled~    | ~double[walk_num][nucl_num][elec_num]~    | in     | Electron-nucleus distances            |
!      | ~en_distance_rescaled_gl~ | ~double[walk_num][nucl_num][elec_num][4]~ | in     | Electron-nucleus distance derivatives |
!      | ~factor_en_gl~            | ~double[walk_num][4][elec_num]~           | out    | Electron-nucleus jastrow              |
!      |---------------------------+-------------------------------------------+--------+---------------------------------------|


function qmckl_compute_jastrow_champ_factor_en_gl_doc( &
     context, walk_num, elec_num, nucl_num, type_nucl_num, &
     type_nucl_vector, aord_num, a_vector, &
     en_distance_rescaled, en_distance_rescaled_gl, factor_en_gl) &
     bind(C) result(info)
  use qmckl
  implicit none

  integer (qmckl_context), intent(in), value :: context
  integer (c_int64_t) , intent(in)  , value :: walk_num
  integer (c_int64_t) , intent(in)  , value :: elec_num
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  integer (c_int64_t) , intent(in)  , value :: type_nucl_num
  integer (c_int64_t) , intent(in)          :: type_nucl_vector(nucl_num)
  integer (c_int64_t) , intent(in)  , value :: aord_num
  real    (c_double ) , intent(in)          :: a_vector(aord_num+1,type_nucl_num)
  real    (c_double ) , intent(in)          :: en_distance_rescaled(elec_num,nucl_num,walk_num)
  real    (c_double ) , intent(in)          :: en_distance_rescaled_gl(4, elec_num,nucl_num,walk_num)
  real    (c_double ) , intent(out)         :: factor_en_gl(elec_num,4,walk_num)
  integer(qmckl_exit_code)                  :: info

  integer*8 :: i, a, k, nw, ii
  double precision   :: x, x1, kf
  double precision   :: denom, invdenom, invdenom2, f
  double precision   :: grad_c2
  double precision   :: dx(4)

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

  if (nucl_num <= 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  if (aord_num < 0) then
     info = QMCKL_INVALID_ARG_7
     return
  endif

  do nw =1, walk_num
     factor_en_gl(:,:,nw) = 0.0d0
     do a = 1, nucl_num
        do i = 1, elec_num

           x = en_distance_rescaled(i,a,nw)
           if(abs(x) < 1.d-12) continue

           denom = 1.0d0 + a_vector(2, type_nucl_vector(a)+1) * x
           invdenom = 1.0d0 / denom
           invdenom2 = invdenom*invdenom

           dx(1) = en_distance_rescaled_gl(1,i,a,nw)
           dx(2) = en_distance_rescaled_gl(2,i,a,nw)
           dx(3) = en_distance_rescaled_gl(3,i,a,nw)
           dx(4) = en_distance_rescaled_gl(4,i,a,nw)

           f = a_vector(1, type_nucl_vector(a)+1) * invdenom2
           grad_c2 = dx(1)*dx(1) + dx(2)*dx(2) + dx(3)*dx(3)

           factor_en_gl(i,1,nw) = factor_en_gl(i,1,nw) + f * dx(1)
           factor_en_gl(i,2,nw) = factor_en_gl(i,2,nw) + f * dx(2)
           factor_en_gl(i,3,nw) = factor_en_gl(i,3,nw) + f * dx(3)
           factor_en_gl(i,4,nw) = factor_en_gl(i,4,nw) &
                + f * (dx(4) - 2.d0 * a_vector(2, type_nucl_vector(a)+1) * grad_c2 * invdenom)


           kf = 2.d0
           x1 = x
           x = 1.d0
           do k=2, aord_num
              f = a_vector(k+1,type_nucl_vector(a)+1) * kf * x
              factor_en_gl(i,1,nw) = factor_en_gl(i,1,nw) + f * x1 * dx(1)
              factor_en_gl(i,2,nw) = factor_en_gl(i,2,nw) + f * x1 * dx(2)
              factor_en_gl(i,3,nw) = factor_en_gl(i,3,nw) + f * x1 * dx(3)
              factor_en_gl(i,4,nw) = factor_en_gl(i,4,nw) &
                   + f * (x1 * dx(4) + (kf-1.d0) * grad_c2)
              x = x*x1
              kf = kf + 1.d0
           end do

        end do
     end do
  end do

end function qmckl_compute_jastrow_champ_factor_en_gl_doc

! Compute
!      :PROPERTIES:
!      :Name:     qmckl_compute_een_rescaled_e
!      :CRetType: qmckl_exit_code
!      :FRetType: qmckl_exit_code
!      :END:

!      #+NAME: qmckl_factor_een_rescaled_e_args
!      |---------------------+----------------------------------------------------+--------+------------------------------------------------------|
!      | Variable            | Type                                               | In/Out | Description                                          |
!      |---------------------+----------------------------------------------------+--------+------------------------------------------------------|
!      | ~context~           | ~qmckl_context~                                    | in     | Global state                                         |
!      | ~walk_num~          | ~int64_t~                                          | in     | Number of walkers                                    |
!      | ~elec_num~          | ~int64_t~                                          | in     | Number of electrons                                  |
!      | ~cord_num~          | ~int64_t~                                          | in     | Order of polynomials                                 |
!      | ~rescale_factor_ee~ | ~double~                                           | in     | Factor to rescale ee distances                       |
!      | ~ee_distance~       | ~double[walk_num][elec_num][elec_num]~             | in     | Electron-electron distances for each walker          |
!      | ~een_rescaled_e~    | ~double[walk_num][0:cord_num][elec_num][elec_num]~ | out    | Electron-electron rescaled distances for each walker |
!      |---------------------+----------------------------------------------------+--------+------------------------------------------------------|


function qmckl_compute_een_rescaled_e_doc( &
     context, walk_num, elec_num, cord_num, rescale_factor_ee,  &
     ee_distance, een_rescaled_e) &
     result(info) bind(C)
  use, intrinsic :: iso_c_binding
  use qmckl
  implicit none
  integer(qmckl_context), intent(in), value :: context
  integer (c_int64_t) , intent(in)  , value :: walk_num
  integer (c_int64_t) , intent(in)  , value :: elec_num
  integer (c_int64_t) , intent(in)  , value :: cord_num
  real    (c_double ) , intent(in)  , value :: rescale_factor_ee
  real    (c_double ) , intent(in)          :: ee_distance(elec_num,elec_num,walk_num)
  real    (c_double ) , intent(out)         :: een_rescaled_e(elec_num,elec_num,0:cord_num,walk_num)
  integer(qmckl_exit_code)                  :: info

  double precision                    :: x
  integer*8                           :: i, j, k, l, nw

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

  if (cord_num < 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  do nw = 1, walk_num
     do l = 0, cord_num
        do j = 1, elec_num
           do i = 1, elec_num
              een_rescaled_e(i, j, l, nw) = dexp(-rescale_factor_ee * ee_distance(i, j, nw))**l
           end do
           een_rescaled_e(j, j, l, nw) = 0.d0
        end do
     end do
  end do

end function qmckl_compute_een_rescaled_e_doc

! Compute
!      :PROPERTIES:
!      :Name:     qmckl_compute_jastrow_champ_factor_een_rescaled_e_gl
!      :CRetType: qmckl_exit_code
!      :FRetType: qmckl_exit_code
!      :END:

!      #+NAME: qmckl_factor_een_rescaled_e_gl_args
!      |---------------------+-------------------------------------------------------+--------+--------------------------------------|
!      | Variable            | Type                                                  | In/Out | Description                          |
!      |---------------------+-------------------------------------------------------+--------+--------------------------------------|
!      | ~context~           | ~qmckl_context~                                       | in     | Global state                         |
!      | ~walk_num~          | ~int64_t~                                             | in     | Number of walkers                    |
!      | ~elec_num~          | ~int64_t~                                             | in     | Number of electrons                  |
!      | ~cord_num~          | ~int64_t~                                             | in     | Order of polynomials                 |
!      | ~rescale_factor_ee~ | ~double~                                              | in     | Factor to rescale ee distances       |
!      | ~coord_ee~          | ~double[3][walk_num][elec_num]~                       | in     | Electron coordinates                 |
!      | ~ee_distance~       | ~double[walk_num][elec_num][elec_num]~                | in     | Electron-electron distances          |
!      | ~een_rescaled_e~    | ~double[walk_num][0:cord_num][elec_num][elec_num]~    | in     | Electron-electron distances          |
!      | ~een_rescaled_e_gl~ | ~double[walk_num][0:cord_num][elec_num][4][elec_num]~ | out    | Electron-electron rescaled distances |
!      |---------------------+-------------------------------------------------------+--------+--------------------------------------|


function qmckl_compute_jastrow_champ_factor_een_rescaled_e_gl_doc( &
     context, walk_num, elec_num, cord_num, rescale_factor_ee,  &
     coord_ee, ee_distance, een_rescaled_e, een_rescaled_e_gl) &
     result(info) bind(C)
  use, intrinsic :: iso_c_binding
  use qmckl
  implicit none
  integer(qmckl_context), intent(in), value  :: context
  integer(c_int64_t)    , intent(in), value  :: walk_num
  integer(c_int64_t)    , intent(in), value  :: elec_num
  integer(c_int64_t)    , intent(in), value  :: cord_num
  real(c_double)        , intent(in), value  :: rescale_factor_ee
  real(c_double)        , intent(in)  :: coord_ee(elec_num,walk_num,3)
  real(c_double)        , intent(in)  :: ee_distance(elec_num,elec_num,walk_num)
  real(c_double)        , intent(in)  :: een_rescaled_e(elec_num,elec_num,0:cord_num,walk_num)
  real(c_double)        , intent(out) :: een_rescaled_e_gl(elec_num,4,elec_num,0:cord_num,walk_num)
  integer(qmckl_exit_code)            :: info

  double precision, allocatable       :: elec_dist_gl(:,:,:)
  double precision                    :: x, kappa_l
  integer*8                           :: i, j, k, l, nw, ii

  double precision  :: rij_inv(elec_num)


  allocate(elec_dist_gl(elec_num, 4, elec_num))
  elec_dist_gl = 0.d0
  een_rescaled_e_gl = 0.d0

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

  if (cord_num < 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  do nw = 1, walk_num

     ! Prepare table of exponentiated distances raised to appropriate power
     do j = 1, elec_num
        do i = 1, j-1
           rij_inv(i) = 1.0d0 / ee_distance(i, j, nw)
        enddo
        rij_inv(j) = 0.0d0
        do i = j+1, elec_num
           rij_inv(i) = 1.0d0 / ee_distance(i, j, nw)
        enddo
        do i = 1, elec_num
           do ii = 1, 3
              elec_dist_gl(i, ii, j) = (coord_ee(i, nw, ii) - coord_ee(j, nw, ii)) * rij_inv(i)
           end do
           elec_dist_gl(i, 4, j) = 2.0d0 * rij_inv(i)
        end do
     end do

     !   Not necessary: should be set to zero by qmckl_malloc
     een_rescaled_e_gl(:,:,:,0,nw) = 0.d0

     do l = 1, cord_num
        kappa_l = -dble(l) * rescale_factor_ee
        do j = 1, elec_num
           do i = 1, elec_num
              if (i /= j) then
                 een_rescaled_e_gl(i, 1, j, l, nw) = kappa_l *  elec_dist_gl(i, 1, j) * een_rescaled_e(i,j,l,nw)
                 een_rescaled_e_gl(i, 2, j, l, nw) = kappa_l *  elec_dist_gl(i, 2, j) * een_rescaled_e(i,j,l,nw)
                 een_rescaled_e_gl(i, 3, j, l, nw) = kappa_l *  elec_dist_gl(i, 3, j) * een_rescaled_e(i,j,l,nw)
                 een_rescaled_e_gl(i, 4, j, l, nw) = kappa_l * (elec_dist_gl(i, 4, j) + kappa_l) * een_rescaled_e(i,j,l,nw)
              else
                 een_rescaled_e_gl(i, 1, j, l, nw) = 0.d0
                 een_rescaled_e_gl(i, 2, j, l, nw) = 0.d0
                 een_rescaled_e_gl(i, 3, j, l, nw) = 0.d0
                 een_rescaled_e_gl(i, 4, j, l, nw) = 0.d0
              end if
           end do
        end do
     end do
  end do

end function qmckl_compute_jastrow_champ_factor_een_rescaled_e_gl_doc

! Compute
!      :PROPERTIES:
!      :Name:     qmckl_compute_een_rescaled_n
!      :CRetType: qmckl_exit_code
!      :FRetType: qmckl_exit_code
!      :END:

!      #+NAME: qmckl_factor_een_rescaled_n_args
!      |---------------------+----------------------------------------------------+--------+-------------------------------------|
!      | Variable            | Type                                               | In/Out | Description                         |
!      |---------------------+----------------------------------------------------+--------+-------------------------------------|
!      | ~context~           | ~qmckl_context~                                    | in     | Global state                        |
!      | ~walk_num~          | ~int64_t~                                          | in     | Number of walkers                   |
!      | ~elec_num~          | ~int64_t~                                          | in     | Number of electrons                 |
!      | ~nucl_num~          | ~int64_t~                                          | in     | Number of atoms                     |
!      | ~type_nucl_num~     | ~int64_t~                                          | in     | Number of atom types                |
!      | ~type_nucl_vector~  | ~int64_t[nucl_num]~                                | in     | Types of atoms                      |
!      | ~cord_num~          | ~int64_t~                                          | in     | Order of polynomials                |
!      | ~rescale_factor_en~ | ~double[nucl_num]~                                 | in     | Factor to rescale ee distances      |
!      | ~en_distance~       | ~double[walk_num][elec_num][nucl_num]~             | in     | Electron-nucleus distances          |
!      | ~een_rescaled_n~    | ~double[walk_num][0:cord_num][nucl_num][elec_num]~ | out    | Electron-nucleus rescaled distances |
!      |---------------------+----------------------------------------------------+--------+-------------------------------------|


function qmckl_compute_een_rescaled_n( &
     context, walk_num, elec_num, nucl_num, &
     type_nucl_num, type_nucl_vector, cord_num, rescale_factor_en,  &
     en_distance, een_rescaled_n) &
     result(info) bind(C)
  use, intrinsic :: iso_c_binding
  use qmckl
  implicit none
  integer(qmckl_context), intent(in), value :: context
  integer (c_int64_t) , intent(in)  , value :: walk_num
  integer (c_int64_t) , intent(in)  , value :: elec_num
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  integer (c_int64_t) , intent(in)  , value :: type_nucl_num
  integer (c_int64_t) , intent(in)          :: type_nucl_vector(nucl_num)
  integer (c_int64_t) , intent(in)  , value :: cord_num
  real    (c_double ) , intent(in)          :: rescale_factor_en(nucl_num)
  real    (c_double ) , intent(in)          :: en_distance(nucl_num,elec_num,walk_num)
  real    (c_double ) , intent(out)         :: een_rescaled_n(elec_num,nucl_num,0:cord_num,walk_num)
  integer(qmckl_exit_code)                  :: info

  double precision                    :: x
  integer*8                           :: i, a, k, l, nw

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

  if (nucl_num <= 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  if (cord_num < 0) then
     info = QMCKL_INVALID_ARG_5
     return
  endif

  do nw = 1, walk_num

     ! prepare the actual een table
     een_rescaled_n(:, :, 0, nw) = 1.0d0

     do a = 1, nucl_num
        do i = 1, elec_num
           een_rescaled_n(i, a, 1, nw) = dexp(-rescale_factor_en(type_nucl_vector(a)+1) * en_distance(a, i, nw))
        end do
     end do

     do l = 2, cord_num
        do a = 1, nucl_num
           do i = 1, elec_num
              een_rescaled_n(i, a, l, nw) = een_rescaled_n(i, a, l - 1, nw) * een_rescaled_n(i, a, 1, nw)
           end do
        end do
     end do

  end do

end function qmckl_compute_een_rescaled_n

! Compute
!      :PROPERTIES:
!      :Name:     qmckl_compute_jastrow_champ_factor_een_rescaled_n_gl
!      :CRetType: qmckl_exit_code
!      :FRetType: qmckl_exit_code
!      :END:

!      #+NAME: qmckl_compute_jastrow_champ_factor_een_rescaled_n_gl_args
!      |---------------------+-------------------------------------------------------+--------+-------------------------------------|
!      | Variable            | Type                                                  | In/Out | Description                         |
!      |---------------------+-------------------------------------------------------+--------+-------------------------------------|
!      | ~context~           | ~qmckl_context~                                       | in     | Global state                        |
!      | ~walk_num~          | ~int64_t~                                             | in     | Number of walkers                   |
!      | ~elec_num~          | ~int64_t~                                             | in     | Number of electrons                 |
!      | ~nucl_num~          | ~int64_t~                                             | in     | Number of atoms                     |
!      | ~type_nucl_num~     | ~int64_t~                                             | in     | Number of atom types                |
!      | ~type_nucl_vector~  | ~int64_t[nucl_num]~                                   | in     | Types of atoms                      |
!      | ~cord_num~          | ~int64_t~                                             | in     | Order of polynomials                |
!      | ~rescale_factor_en~ | ~double[nucl_num]~                                    | in     | Factor to rescale ee distances      |
!      | ~coord_ee~          | ~double[3][walk_num][elec_num]~                       | in     | Electron coordinates                |
!      | ~coord_n~           | ~double[3][nucl_num]~                                 | in     | Nuclear coordinates                 |
!      | ~en_distance~       | ~double[walk_num][elec_num][nucl_num]~                | in     | Electron-nucleus distances          |
!      | ~een_rescaled_n~    | ~double[walk_num][0:cord_num][nucl_num][elec_num]~    | in     | Electron-nucleus distances          |
!      | ~een_rescaled_n_gl~ | ~double[walk_num][0:cord_num][nucl_num][4][elec_num]~ | out    | Electron-nucleus rescaled distances |
!      |---------------------+-------------------------------------------------------+--------+-------------------------------------|


function qmckl_compute_jastrow_champ_factor_een_rescaled_n_gl( &
     context, walk_num, elec_num, nucl_num, type_nucl_num, type_nucl_vector, &
     cord_num, rescale_factor_en, &
     coord_ee, coord_n, en_distance, een_rescaled_n, een_rescaled_n_gl) &
     result(info) bind(C)
  use, intrinsic :: iso_c_binding
  use qmckl
  implicit none
  integer(qmckl_context), intent(in), value :: context
  integer (c_int64_t) , intent(in)  , value :: walk_num
  integer (c_int64_t) , intent(in)  , value :: elec_num
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  integer (c_int64_t) , intent(in)  , value :: type_nucl_num
  integer (c_int64_t) , intent(in)          :: type_nucl_vector(nucl_num)
  integer (c_int64_t) , intent(in)  , value :: cord_num
  real    (c_double ) , intent(in)          :: rescale_factor_en(nucl_num)
  real    (c_double ) , intent(in)          :: coord_ee(elec_num,walk_num,3)
  real    (c_double ) , intent(in)          :: coord_n(nucl_num,3)
  real    (c_double ) , intent(in)          :: en_distance(nucl_num,elec_num,walk_num)
  real    (c_double ) , intent(in)          :: een_rescaled_n(elec_num,nucl_num,0:cord_num,walk_num)
  real    (c_double ) , intent(out)         :: een_rescaled_n_gl(elec_num,4,nucl_num,0:cord_num,walk_num)
  integer(qmckl_exit_code)                  :: info

  double precision,allocatable        :: elnuc_dist_gl(:,:,:)
  double precision                    :: x, ria_inv, kappa_l
  integer*8                           :: i, a, k, l, nw, ii

  allocate(elnuc_dist_gl(elec_num, 4, nucl_num))

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

  if (nucl_num <= 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  if (cord_num < 0) then
     info = QMCKL_INVALID_ARG_5
     return
  endif

  ! Prepare table of exponentiated distances raised to appropriate power
  een_rescaled_n_gl             = 0.0d0
  do nw = 1, walk_num

     ! prepare the actual een table
     do a = 1, nucl_num
        do i = 1, elec_num
           ria_inv = 1.0d0 / en_distance(a, i, nw)
           do ii = 1, 3
              elnuc_dist_gl(i, ii, a) = (coord_ee(i, nw, ii) - coord_n(a, ii)) * ria_inv
           end do
           elnuc_dist_gl(i, 4, a) = 2.0d0 * ria_inv
        end do
     end do

     do l = 0, cord_num
        do a = 1, nucl_num
           kappa_l = - dble(l) * rescale_factor_en(type_nucl_vector(a)+1)
           do i = 1, elec_num
              een_rescaled_n_gl(i, 1, a, l, nw) = kappa_l * elnuc_dist_gl(i, 1, a)
              een_rescaled_n_gl(i, 2, a, l, nw) = kappa_l * elnuc_dist_gl(i, 2, a)
              een_rescaled_n_gl(i, 3, a, l, nw) = kappa_l * elnuc_dist_gl(i, 3, a)
              een_rescaled_n_gl(i, 4, a, l, nw) = kappa_l * elnuc_dist_gl(i, 4, a)

              een_rescaled_n_gl(i, 4, a, l, nw) = een_rescaled_n_gl(i, 4, a, l, nw)           &
                   + een_rescaled_n_gl(i, 1, a, l, nw) * een_rescaled_n_gl(i, 1, a, l, nw) &
                   + een_rescaled_n_gl(i, 2, a, l, nw) * een_rescaled_n_gl(i, 2, a, l, nw) &
                   + een_rescaled_n_gl(i, 3, a, l, nw) * een_rescaled_n_gl(i, 3, a, l, nw)

              een_rescaled_n_gl(i, 1, a, l, nw) = een_rescaled_n_gl(i, 1, a, l, nw) * &
                   een_rescaled_n(i, a, l, nw)
              een_rescaled_n_gl(i, 2, a, l, nw) = een_rescaled_n_gl(i, 2, a, l, nw) * &
                   een_rescaled_n(i, a, l, nw)
              een_rescaled_n_gl(i, 3, a, l, nw) = een_rescaled_n_gl(i, 3, a, l, nw) * &
                   een_rescaled_n(i, a, l, nw)
              een_rescaled_n_gl(i, 4, a, l, nw) = een_rescaled_n_gl(i, 4, a, l, nw) * &
                   een_rescaled_n(i, a, l, nw)
           end do
        end do
     end do
  end do

end function qmckl_compute_jastrow_champ_factor_een_rescaled_n_gl



! #+RESULTS:
! | $N_{ord}$ | Number of parameters |
! |           |                      |
! |         1 |                    0 |
! |         2 |                    2 |
! |         3 |                    6 |
! |         4 |                   13 |
! |         5 |                   23 |
! |         6 |                   37 |
! |         7 |                   55 |
! |         8 |                   78 |
! |         9 |                  106 |
! |        10 |                  140 |

! #+NAME: qmckl_factor_dim_c_vector_args
! |----------------+-----------------+--------+------------------------------------|
! | Variable       | Type            | In/Out | Description                        |
! |----------------+-----------------+--------+------------------------------------|
! | ~context~      | ~qmckl_context~ | in     | Global state                       |
! | ~cord_num~     | ~int64_t~       | in     | Order of polynomials               |
! | ~dim_c_vector~ | ~int64_t~       | out    | Number of parameters per atom type |
! |----------------+-----------------+--------+------------------------------------|


function qmckl_compute_dim_c_vector_doc( &
     context, cord_num, dim_c_vector) &
     result(info) bind(C)
  use, intrinsic :: iso_c_binding
  use qmckl
  implicit none
  integer(qmckl_context), intent(in), value :: context
  integer (c_int64_t) , intent(in)  , value :: cord_num
  integer (c_int64_t) , intent(out)         :: dim_c_vector
  integer(qmckl_exit_code)                  :: info

  double precision                    :: x
  integer*8                           :: i, a, k, l, p, lmax

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (cord_num < 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  dim_c_vector = 0

  do p = 2, cord_num
    do k = p - 1, 0, -1
      if (k .ne. 0) then
        lmax = p - k
      else
        lmax = p - k - 2
      endif
      do l = lmax, 0, -1
        if (iand(p - k - l, 1_8) == 1) cycle
        dim_c_vector = dim_c_vector + 1
      end do
    end do
  end do

end function qmckl_compute_dim_c_vector_doc

! Compute c_vector_full
!      :PROPERTIES:
!      :Name:     qmckl_compute_c_vector_full
!      :CRetType: qmckl_exit_code
!      :FRetType: qmckl_exit_code
!      :END:

!      #+NAME: qmckl_factor_c_vector_full_args
!      |--------------------+---------------------------------------+--------+------------------------------|
!      | Variable           | Type                                  | In/Out | Description                  |
!      |--------------------+---------------------------------------+--------+------------------------------|
!      | ~context~          | ~qmckl_context~                       | in     | Global state                 |
!      | ~nucl_num~         | ~int64_t~                             | in     | Number of atoms              |
!      | ~dim_c_vector~     | ~int64_t~                             | in     | dimension of cord full table |
!      | ~type_nucl_num~    | ~int64_t~                             | in     | dimension of cord full table |
!      | ~type_nucl_vector~ | ~int64_t[nucl_num]~                   | in     | dimension of cord full table |
!      | ~c_vector~         | ~double[dim_c_vector][type_nucl_num]~ | in     | dimension of cord full table |
!      | ~c_vector_full~    | ~double[dim_c_vector][nucl_num]~      | out    | Full list of coefficients    |
!      |--------------------+---------------------------------------+--------+------------------------------|


function qmckl_compute_c_vector_full_doc( &
     context, nucl_num, dim_c_vector, type_nucl_num,  &
     type_nucl_vector, c_vector, c_vector_full) &
     result(info) bind(C)
  use, intrinsic :: iso_c_binding
  use qmckl
  implicit none
  integer(qmckl_context), intent(in), value :: context
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  integer (c_int64_t) , intent(in)  , value :: dim_c_vector
  integer (c_int64_t) , intent(in)  , value :: type_nucl_num
  integer (c_int64_t) , intent(in)          :: type_nucl_vector(nucl_num)
  real    (c_double ) , intent(in)          :: c_vector(dim_c_vector, type_nucl_num)
  real    (c_double ) , intent(out)         :: c_vector_full(nucl_num,dim_c_vector)
  integer(qmckl_exit_code)                  :: info

  double precision                    :: x
  integer*8                           :: i, a, k, l, nw

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) info = QMCKL_INVALID_CONTEXT
  if (nucl_num <= 0)                 info = QMCKL_INVALID_ARG_2
  if (dim_c_vector < 0)              info = QMCKL_INVALID_ARG_3
  if (type_nucl_num <= 0)            info = QMCKL_INVALID_ARG_4
  if (info /= QMCKL_SUCCESS)         return

  do a = 1, nucl_num
    c_vector_full(a,1:dim_c_vector) = c_vector(1:dim_c_vector, type_nucl_vector(a)+1)
  end do

end function qmckl_compute_c_vector_full_doc

! Compute lkpm_combined_index
!      :PROPERTIES:
!      :Name:     qmckl_compute_lkpm_combined_index
!      :CRetType: qmckl_exit_code
!      :FRetType: qmckl_exit_code
!      :END:

!      #+NAME: lkpm_combined_index_args
!      |-----------------------+----------------------------+--------+-------------------------------|
!      | Variable              | Type                       | In/Out | Description                   |
!      |-----------------------+----------------------------+--------+-------------------------------|
!      | ~context~             | ~qmckl_context~            | in     | Global state                  |
!      | ~cord_num~            | ~int64_t~                  | in     | Order of polynomials          |
!      | ~dim_c_vector~        | ~int64_t~                  | in     | dimension of cord full table  |
!      | ~lkpm_combined_index~ | ~int64_t[4][dim_c_vector]~ | out    | Full list of combined indices |
!      |-----------------------+----------------------------+--------+-------------------------------|


function qmckl_compute_lkpm_combined_index_doc( &
     context, cord_num, dim_c_vector,  lkpm_combined_index) &
     result(info) bind(C)
  use, intrinsic :: iso_c_binding
  use qmckl
  implicit none
  integer(qmckl_context), intent(in), value :: context
  integer (c_int64_t) , intent(in)  , value :: cord_num
  integer (c_int64_t) , intent(in)  , value :: dim_c_vector
  integer (c_int64_t) , intent(out)         :: lkpm_combined_index(dim_c_vector,4)
  integer(qmckl_exit_code)                  :: info

  double precision                    :: x
  integer*8                           :: i, a, k, l, kk, p, lmax, m

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) info = QMCKL_INVALID_CONTEXT
  if (cord_num < 0)                  info = QMCKL_INVALID_ARG_2
  if (dim_c_vector < 0)              info = QMCKL_INVALID_ARG_3
  if (info /= QMCKL_SUCCESS)         return

  kk = 0
  do p = 2, cord_num
    do k = p - 1, 0, -1
      if (k /= 0) then
        lmax = p - k
      else
        lmax = p - k - 2
      end if
      do l = lmax, 0, -1
        if (iand(p - k - l, 1_8) .eq. 1_8) cycle
        m = (p - k - l)/2
        kk = kk + 1
        lkpm_combined_index(kk, 1) = l
        lkpm_combined_index(kk, 2) = k
        lkpm_combined_index(kk, 3) = p
        lkpm_combined_index(kk, 4) = m
      end do
    end do
  end do

end function qmckl_compute_lkpm_combined_index_doc

function qmckl_compute_tmp_c_doc( &
     context, cord_num, elec_num, nucl_num, &
     walk_num, een_rescaled_e, een_rescaled_n, tmp_c) &
     result(info) bind(C)
  use, intrinsic :: iso_c_binding
  use qmckl
  implicit none
  integer(qmckl_context), intent(in), value  :: context
  integer (c_int64_t) , intent(in)  , value :: cord_num
  integer (c_int64_t) , intent(in)  , value :: elec_num
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  integer (c_int64_t) , intent(in)  , value :: walk_num
  real    (c_double ) , intent(in)          :: een_rescaled_e(elec_num,elec_num,0:cord_num,walk_num)
  real    (c_double ) , intent(in)          :: een_rescaled_n(elec_num,nucl_num,0:cord_num,walk_num)
  real    (c_double ) , intent(out)         :: tmp_c(elec_num,nucl_num,0:cord_num,0:cord_num-1,walk_num)
  integer(qmckl_exit_code)                  :: info

  double precision                    :: x
  integer*8                           :: i, j, a, l, kk, p, lmax, nw
  character                           :: TransA, TransB
  double precision                    :: alpha, beta
  integer*8                           :: M, N, K, LDA, LDB, LDC

  TransA = 'N'
  TransB = 'N'
  alpha = 1.0d0
  beta  = 0.0d0

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) info = QMCKL_INVALID_CONTEXT
  if (cord_num <  0)                 info = QMCKL_INVALID_ARG_2
  if (elec_num <= 0)                 info = QMCKL_INVALID_ARG_3
  if (nucl_num <= 0)                 info = QMCKL_INVALID_ARG_4
  if (walk_num <= 0)                 info = QMCKL_INVALID_ARG_5
  if (info /= QMCKL_SUCCESS)         return


  M = elec_num
  N = nucl_num*(cord_num + 1)
  K = elec_num
  LDA = size(een_rescaled_e,1)
  LDB = size(een_rescaled_n,1)
  LDC = size(tmp_c,1)

  do nw=1, walk_num
    do i=0, cord_num-1
      info = qmckl_dgemm(context, TransA, TransB, M, N, K, alpha,     &
                         een_rescaled_e(1,1,i,nw),LDA*1_8,                     &
                         een_rescaled_n(1,1,0,nw),LDB*1_8,                     &
                         beta,                                       &
                         tmp_c(1,1,0,i,nw),LDC)
    end do
  end do

end function qmckl_compute_tmp_c_doc

function qmckl_compute_dtmp_c_doc( &
     context, cord_num, elec_num, nucl_num, &
     walk_num, een_rescaled_e_gl, een_rescaled_n, dtmp_c) &
     result(info) bind(C)
  use, intrinsic :: iso_c_binding
  use qmckl
  implicit none
  integer(qmckl_context), intent(in), value :: context
  integer (c_int64_t) , intent(in)  , value :: cord_num
  integer (c_int64_t) , intent(in)  , value :: elec_num
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  integer (c_int64_t) , intent(in)  , value :: walk_num
  real    (c_double ) , intent(in)          :: een_rescaled_e_gl(elec_num,4,elec_num,0:cord_num,walk_num)
  real    (c_double ) , intent(in)          :: een_rescaled_n(elec_num,nucl_num,0:cord_num,walk_num)
  real    (c_double ) , intent(out)         :: dtmp_c(elec_num,4,nucl_num,0:cord_num,0:cord_num-1,walk_num)
  integer(qmckl_exit_code)                  :: info

  double precision                    :: x
  integer*8                           :: i, j, a, l, kk, p, lmax, nw, ii
  character                           :: TransA, TransB
  double precision                    :: alpha, beta
  integer*8                           :: M, N, K, LDA, LDB, LDC

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) info = QMCKL_INVALID_CONTEXT
  if (cord_num <  0)                 info = QMCKL_INVALID_ARG_2
  if (elec_num <= 0)                 info = QMCKL_INVALID_ARG_3
  if (nucl_num <= 0)                 info = QMCKL_INVALID_ARG_4
  if (walk_num <= 0)                 info = QMCKL_INVALID_ARG_5
  if (info /= QMCKL_SUCCESS)         return

  TransA = 'N'
  TransB = 'N'
  alpha = 1.0d0
  beta  = 0.0d0

  M = 4*elec_num
  N = nucl_num*(cord_num + 1)
  K = elec_num
  LDA = 4*size(een_rescaled_e_gl,1)
  LDB = size(een_rescaled_n,1)
  LDC = 4*size(dtmp_c,1)

  do nw=1, walk_num
     do i=0, cord_num-1
        info = qmckl_dgemm(context,TransA, TransB, M, N, K, alpha,  &
             een_rescaled_e_gl(1,1,1,i,nw),LDA*1_8,            &
             een_rescaled_n(1,1,0,nw),LDB*1_8,                      &
             beta,                                                  &
             dtmp_c(1,1,1,0,i,nw),LDC)
     end do
  end do

end function qmckl_compute_dtmp_c_doc

! Compute naive
!      :PROPERTIES:
!      :Name:     qmckl_compute_jastrow_champ_factor_een_naive
!      :CRetType: qmckl_exit_code
!      :FRetType: qmckl_exit_code
!      :END:

!      #+NAME: qmckl_factor_een_naive_args
!      |-----------------------+----------------------------------------------------+--------+--------------------------------------|
!      | Variable              | Type                                               | In/Out | Description                          |
!      |-----------------------+----------------------------------------------------+--------+--------------------------------------|
!      | ~context~             | ~qmckl_context~                                    | in     | Global state                         |
!      | ~walk_num~            | ~int64_t~                                          | in     | Number of walkers                    |
!      | ~elec_num~            | ~int64_t~                                          | in     | Number of electrons                  |
!      | ~nucl_num~            | ~int64_t~                                          | in     | Number of nuclei                     |
!      | ~cord_num~            | ~int64_t~                                          | in     | order of polynomials                 |
!      | ~dim_c_vector~        | ~int64_t~                                          | in     | dimension of full coefficient vector |
!      | ~c_vector_full~       | ~double[dim_c_vector][nucl_num]~                   | in     | full coefficient vector              |
!      | ~lkpm_combined_index~ | ~int64_t[4][dim_c_vector]~                         | in     | combined indices                     |
!      | ~een_rescaled_e~      | ~double[walk_num][0:cord_num][elec_num][elec_num]~ | in     | Electron-nucleus rescaled            |
!      | ~een_rescaled_n~      | ~double[walk_num][0:cord_num][nucl_num][elec_num]~ | in     | Electron-nucleus rescaled factor     |
!      | ~factor_een~          | ~double[walk_num]~                                 | out    | Electron-nucleus jastrow             |
!      |-----------------------+----------------------------------------------------+--------+--------------------------------------|


function qmckl_compute_jastrow_champ_factor_een_naive( &
     context, walk_num, elec_num, nucl_num, cord_num,&
     dim_c_vector, c_vector_full, lkpm_combined_index, &
     een_rescaled_e, een_rescaled_n, factor_een) &
     result(info) bind(C)
  use, intrinsic :: iso_c_binding
  use qmckl
  implicit none
  integer(qmckl_context), intent(in), value :: context
  integer (c_int64_t) , intent(in)  , value :: walk_num
  integer (c_int64_t) , intent(in)  , value :: elec_num
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  integer (c_int64_t) , intent(in)  , value :: cord_num
  integer (c_int64_t) , intent(in)  , value :: dim_c_vector
  real    (c_double ) , intent(in)          :: c_vector_full(nucl_num,dim_c_vector)
  integer (c_int64_t) , intent(in)          :: lkpm_combined_index(dim_c_vector,4)
  real    (c_double ) , intent(in)          :: een_rescaled_e(elec_num,elec_num,0:cord_num,walk_num)
  real    (c_double ) , intent(in)          :: een_rescaled_n(elec_num,nucl_num,0:cord_num,walk_num)
  real    (c_double ) , intent(out)         :: factor_een(walk_num)
  integer(qmckl_exit_code)                  :: info

  integer*8 :: i, a, j, l, k, m, n, p, nw
  double precision :: cn

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) info = QMCKL_INVALID_CONTEXT
  if (walk_num <= 0)                 info = QMCKL_INVALID_ARG_2
  if (elec_num <= 0)                 info = QMCKL_INVALID_ARG_3
  if (nucl_num <= 0)                 info = QMCKL_INVALID_ARG_4
  if (cord_num <  0)                 info = QMCKL_INVALID_ARG_5
  if (info /= QMCKL_SUCCESS)         return

  do nw =1, walk_num
     factor_een(nw) = 0.d0
     do n = 1, dim_c_vector
        l = lkpm_combined_index(n, 1)
        k = lkpm_combined_index(n, 2)
        p = lkpm_combined_index(n, 3)
        m = lkpm_combined_index(n, 4)

        do a = 1, nucl_num
           cn = c_vector_full(a, n)
           if (cn == 0.d0) cycle
           do j = 1, elec_num
              do i = 1, j-1
                 factor_een(nw) = factor_een(nw) + cn*( &
                      een_rescaled_e(i,j,k,nw) *       &
                      (een_rescaled_n(i,a,l,nw) + een_rescaled_n(j,a,l,nw)) * &
                      (een_rescaled_n(i,a,m,nw) * een_rescaled_n(j,a,m,nw)) )
              end do
           end do
        end do
     end do
  end do

end function qmckl_compute_jastrow_champ_factor_een_naive

! Compute
!      :PROPERTIES:
!      :Name:     qmckl_compute_jastrow_champ_factor_een_doc
!      :CRetType: qmckl_exit_code
!      :FRetType: qmckl_exit_code
!      :END:

!      #+NAME: qmckl_factor_een_args
!      |-----------------------+------------------------------------------------------------------+--------+--------------------------------------|
!      | Variable              | Type                                                             | In/Out | Description                          |
!      |-----------------------+------------------------------------------------------------------+--------+--------------------------------------|
!      | ~context~             | ~qmckl_context~                                                  | in     | Global state                         |
!      | ~walk_num~            | ~int64_t~                                                        | in     | Number of walkers                    |
!      | ~elec_num~            | ~int64_t~                                                        | in     | Number of electrons                  |
!      | ~nucl_num~            | ~int64_t~                                                        | in     | Number of nuclei                     |
!      | ~cord_num~            | ~int64_t~                                                        | in     | order of polynomials                 |
!      | ~dim_c_vector~        | ~int64_t~                                                        | in     | dimension of full coefficient vector |
!      | ~c_vector_full~       | ~double[dim_c_vector][nucl_num]~                                 | in     | full coefficient vector              |
!      | ~lkpm_combined_index~ | ~int64_t[4][dim_c_vector]~                                       | in     | combined indices                     |
!      | ~tmp_c~               | ~double[walk_num][0:cord_num-1][0:cord_num][nucl_num][elec_num]~ | in     | vector of non-zero coefficients      |
!      | ~een_rescaled_n~      | ~double[walk_num][0:cord_num][nucl_num][elec_num]~               | in     | Electron-nucleus rescaled distances  |
!      | ~factor_een~          | ~double[walk_num]~                                               | out    | Electron-nucleus jastrow             |
!      |-----------------------+------------------------------------------------------------------+--------+--------------------------------------|


function qmckl_compute_jastrow_champ_factor_een_doc( &
     context, walk_num, elec_num, nucl_num, cord_num,   &
     dim_c_vector, c_vector_full, lkpm_combined_index, &
     tmp_c, een_rescaled_n, factor_een) &
     result(info) bind(C)
  use, intrinsic :: iso_c_binding
  use qmckl
  implicit none
  integer(qmckl_context), intent(in), value :: context
  integer (c_int64_t) , intent(in)  , value :: walk_num
  integer (c_int64_t) , intent(in)  , value :: elec_num
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  integer (c_int64_t) , intent(in)  , value :: cord_num
  integer (c_int64_t) , intent(in)  , value :: dim_c_vector
  real    (c_double ) , intent(in)          :: c_vector_full(nucl_num,dim_c_vector)
  integer (c_int64_t) , intent(in)          :: lkpm_combined_index(dim_c_vector,4)
  real    (c_double ) , intent(in)          :: tmp_c(elec_num,nucl_num,0:cord_num,0:cord_num-1,walk_num)
  real    (c_double ) , intent(in)          :: een_rescaled_n(elec_num,nucl_num,0:cord_num,walk_num)
  real    (c_double ) , intent(out)         :: factor_een(walk_num)
  integer(qmckl_exit_code)                  :: info

  integer*8 :: i, a, j, l, k, p, m, n, nw
  double precision :: accu, accu2, cn

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) info = QMCKL_INVALID_CONTEXT
  if (walk_num <= 0)                 info = QMCKL_INVALID_ARG_2
  if (elec_num <= 0)                 info = QMCKL_INVALID_ARG_3
  if (nucl_num <= 0)                 info = QMCKL_INVALID_ARG_4
  if (cord_num <  0)                 info = QMCKL_INVALID_ARG_5
  if (info /= QMCKL_SUCCESS)         return

  factor_een = 0.0d0

  if (cord_num == 0) return

  do nw =1, walk_num
     do n = 1, dim_c_vector
        l = lkpm_combined_index(n, 1)
        k = lkpm_combined_index(n, 2)
        m = lkpm_combined_index(n, 4)

        do a = 1, nucl_num
           cn = c_vector_full(a, n)
           if(cn == 0.d0) cycle

           accu = 0.0d0
           do j = 1, elec_num
              accu = accu + een_rescaled_n(j,a,m,nw) * tmp_c(j,a,m+l,k,nw)
           end do
           factor_een(nw) = factor_een(nw) + accu * cn
        end do
     end do
  end do

end function qmckl_compute_jastrow_champ_factor_een_doc

! Compute Naive
!      :PROPERTIES:
!      :Name:     qmckl_compute_jastrow_champ_factor_een_gl_naive
!      :CRetType: qmckl_exit_code
!      :FRetType: qmckl_exit_code
!      :END:

!      #+NAME: qmckl_factor_een_gl_naive_args
!      |-----------------------+-------------------------------------------------------+--------+--------------------------------------|
!      | Variable              | Type                                                  | In/Out | Description                          |
!      |-----------------------+-------------------------------------------------------+--------+--------------------------------------|
!      | ~context~             | ~qmckl_context~                                       | in     | Global state                         |
!      | ~walk_num~            | ~int64_t~                                             | in     | Number of walkers                    |
!      | ~elec_num~            | ~int64_t~                                             | in     | Number of electrons                  |
!      | ~nucl_num~            | ~int64_t~                                             | in     | Number of nuclei                     |
!      | ~cord_num~            | ~int64_t~                                             | in     | order of polynomials                 |
!      | ~dim_c_vector~        | ~int64_t~                                             | in     | dimension of full coefficient vector |
!      | ~c_vector_full~       | ~double[dim_c_vector][nucl_num]~                      | in     | full coefficient vector              |
!      | ~lkpm_combined_index~ | ~int64_t[4][dim_c_vector]~                            | in     | combined indices                     |
!      | ~een_rescaled_e~      | ~double[walk_num][0:cord_num][elec_num][elec_num]~    | in     | Electron-nucleus rescaled            |
!      | ~een_rescaled_n~      | ~double[walk_num][0:cord_num][nucl_num][elec_num]~    | in     | Electron-nucleus rescaled factor     |
!      | ~een_rescaled_e_gl~   | ~double[walk_num][0:cord_num][elec_num][4][elec_num]~ | in     | Electron-nucleus rescaled            |
!      | ~een_rescaled_n_gl~   | ~double[walk_num][0:cord_num][nucl_num][4][elec_num]~ | in     | Electron-nucleus rescaled factor     |
!      | ~factor_een_gl~       | ~double[walk_num][4][elec_num]~                       | out    | Electron-nucleus jastrow             |
!      |-----------------------+-------------------------------------------------------+--------+--------------------------------------|


function qmckl_compute_jastrow_champ_factor_een_gl_naive( &
     context, walk_num, elec_num, nucl_num, cord_num, dim_c_vector, &
     c_vector_full, lkpm_combined_index, een_rescaled_e, een_rescaled_n, &
     een_rescaled_e_gl, een_rescaled_n_gl, factor_een_gl)&
     result(info) bind(C)
  use, intrinsic :: iso_c_binding
  use qmckl
  implicit none
  integer(qmckl_context), intent(in), value :: context
  integer (c_int64_t) , intent(in)  , value :: walk_num
  integer (c_int64_t) , intent(in)  , value :: elec_num
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  integer (c_int64_t) , intent(in)  , value :: cord_num
  integer (c_int64_t) , intent(in)  , value :: dim_c_vector
  real    (c_double ) , intent(in)          :: c_vector_full(nucl_num,dim_c_vector)
  integer (c_int64_t) , intent(in)          :: lkpm_combined_index(dim_c_vector,4)
  real    (c_double ) , intent(in)          :: een_rescaled_e(elec_num,elec_num,0:cord_num,walk_num)
  real    (c_double ) , intent(in)          :: een_rescaled_n(elec_num,nucl_num,0:cord_num,walk_num)
  real    (c_double ) , intent(in)          :: een_rescaled_e_gl(elec_num,4,elec_num,0:cord_num,walk_num)
  real    (c_double ) , intent(in)          :: een_rescaled_n_gl(elec_num,4,nucl_num,0:cord_num,walk_num)
  real    (c_double ) , intent(out)         :: factor_een_gl(elec_num,4,walk_num)
  integer(qmckl_exit_code)                  :: info

  integer*8 :: i, a, j, l, k, m, n, nw
  double precision :: accu, accu2, cn
  double precision :: daccu(1:4), daccu2(1:4)

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) info = QMCKL_INVALID_CONTEXT
  if (walk_num <= 0)                 info = QMCKL_INVALID_ARG_2
  if (elec_num <= 0)                 info = QMCKL_INVALID_ARG_3
  if (nucl_num <= 0)                 info = QMCKL_INVALID_ARG_4
  if (cord_num <  0)                 info = QMCKL_INVALID_ARG_5
  if (info /= QMCKL_SUCCESS)         return

  factor_een_gl = 0.0d0

  do nw =1, walk_num
     do n = 1, dim_c_vector
        l = lkpm_combined_index(n, 1)
        k = lkpm_combined_index(n, 2)
        m = lkpm_combined_index(n, 4)

        do a = 1, nucl_num
           cn = c_vector_full(a, n)
           do j = 1, elec_num
              accu = 0.0d0
              accu2 = 0.0d0
              daccu = 0.0d0
              daccu2 = 0.0d0
              do i = 1, elec_num
                 accu = accu + een_rescaled_e(i, j, k, nw) * een_rescaled_n(i, a, m, nw)
                 accu2 = accu2 + een_rescaled_e(i, j, k, nw) * een_rescaled_n(i, a, m + l, nw)
                 daccu(1:4) = daccu(1:4) + een_rescaled_e_gl(j, 1:4, i, k, nw) *   &
                      een_rescaled_n(i, a, m, nw)
                 daccu2(1:4) = daccu2(1:4) + een_rescaled_e_gl(j, 1:4, i, k, nw) * &
                      een_rescaled_n(i, a, m + l, nw)
              end do
              factor_een_gl(j, 1:4, nw) = factor_een_gl(j, 1:4, nw) +   &
                   (accu * een_rescaled_n_gl(j, 1:4, a, m + l, nw)      &
                   + daccu(1:4) * een_rescaled_n(j, a, m + l, nw)       &
                   + daccu2(1:4) * een_rescaled_n(j, a, m, nw)          &
                   + accu2 * een_rescaled_n_gl(j, 1:4, a, m, nw)) * cn

              factor_een_gl(j, 4, nw) = factor_een_gl(j, 4, nw) + 2.0d0 * ( &
                   daccu (1) * een_rescaled_n_gl(j, 1, a, m + l, nw) +      &
                   daccu (2) * een_rescaled_n_gl(j, 2, a, m + l, nw) +      &
                   daccu (3) * een_rescaled_n_gl(j, 3, a, m + l, nw) +      &
                   daccu2(1) * een_rescaled_n_gl(j, 1, a, m, nw    ) +      &
                   daccu2(2) * een_rescaled_n_gl(j, 2, a, m, nw    ) +      &
                   daccu2(3) * een_rescaled_n_gl(j, 3, a, m, nw    ) ) * cn

           end do
        end do
     end do
  end do

end function qmckl_compute_jastrow_champ_factor_een_gl_naive

! Compute GL
!      :PROPERTIES:
!      :Name:     qmckl_compute_jastrow_champ_factor_een_gl
!      :CRetType: qmckl_exit_code
!      :FRetType: qmckl_exit_code
!      :END:

!      #+NAME: qmckl_factor_een_gl_args
!      |-----------------------+---------------------------------------------------------------------+--------+------------------------------------------------|
!      | Variable              | Type                                                                | In/Out | Description                                    |
!      |-----------------------+---------------------------------------------------------------------+--------+------------------------------------------------|
!      | ~context~             | ~qmckl_context~                                                     | in     | Global state                                   |
!      | ~walk_num~            | ~int64_t~                                                           | in     | Number of walkers                              |
!      | ~elec_num~            | ~int64_t~                                                           | in     | Number of electrons                            |
!      | ~nucl_num~            | ~int64_t~                                                           | in     | Number of nuclei                               |
!      | ~cord_num~            | ~int64_t~                                                           | in     | order of polynomials                           |
!      | ~dim_c_vector~        | ~int64_t~                                                           | in     | dimension of full coefficient vector           |
!      | ~c_vector_full~       | ~double[dim_c_vector][nucl_num]~                                    | in     | full coefficient vector                        |
!      | ~lkpm_combined_index~ | ~int64_t[4][dim_c_vector]~                                          | in     | combined indices                               |
!      | ~tmp_c~               | ~double[walk_num][0:cord_num-1][0:cord_num][nucl_num][elec_num]~    | in     | Temporary intermediate tensor                  |
!      | ~dtmp_c~              | ~double[walk_num][0:cord_num-1][0:cord_num][nucl_num][4][elec_num]~ | in     | vector of non-zero coefficients                |
!      | ~een_rescaled_n~      | ~double[walk_num][0:cord_num][nucl_num][elec_num]~                  | in     | Electron-nucleus rescaled factor               |
!      | ~een_rescaled_n_gl~   | ~double[walk_num][0:cord_num][nucl_num][4][elec_num]~               | in     | Derivative of Electron-nucleus rescaled factor |
!      | ~factor_een_gl~       | ~double[walk_num][4][elec_num]~                                     | out    | Derivative of Electron-nucleus jastrow         |
!      |-----------------------+---------------------------------------------------------------------+--------+------------------------------------------------|



function qmckl_compute_jastrow_champ_factor_een_gl_doc( &
     context, walk_num, elec_num, nucl_num, &
     cord_num, dim_c_vector, c_vector_full, lkpm_combined_index, &
     tmp_c, dtmp_c, een_rescaled_n, een_rescaled_n_gl, factor_een_gl)&
     result(info) bind(C)
  use, intrinsic :: iso_c_binding
  use qmckl
  implicit none
  integer(qmckl_context), intent(in), value :: context
  integer (c_int64_t) , intent(in)  , value :: walk_num
  integer (c_int64_t) , intent(in)  , value :: elec_num
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  integer (c_int64_t) , intent(in)  , value :: cord_num
  integer (c_int64_t) , intent(in)  , value :: dim_c_vector
  real    (c_double ) , intent(in)          :: c_vector_full(nucl_num,dim_c_vector)
  integer (c_int64_t) , intent(in)          :: lkpm_combined_index(dim_c_vector,4)
  real    (c_double ) , intent(in)          :: tmp_c(elec_num,nucl_num,0:cord_num,0:cord_num-1,walk_num)
  real    (c_double ) , intent(in)          :: dtmp_c(elec_num,4,nucl_num,0:cord_num,0:cord_num-1,walk_num)
  real    (c_double ) , intent(in)          :: een_rescaled_n(elec_num,nucl_num,0:cord_num,walk_num)
  real    (c_double ) , intent(in)          :: een_rescaled_n_gl(elec_num,4,nucl_num,0:cord_num,walk_num)
  real    (c_double ) , intent(out)         :: factor_een_gl(elec_num,4,walk_num)
  integer(qmckl_exit_code)                  :: info

  integer*8 :: i, a, j, l, k, m, n, nw, ii
  double precision :: accu, accu2, cn

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) info = QMCKL_INVALID_CONTEXT
  if (walk_num <= 0)                 info = QMCKL_INVALID_ARG_2
  if (elec_num <= 0)                 info = QMCKL_INVALID_ARG_3
  if (nucl_num <= 0)                 info = QMCKL_INVALID_ARG_4
  if (cord_num <  0)                 info = QMCKL_INVALID_ARG_5
  if (info /= QMCKL_SUCCESS)         return

  if (cord_num == 0) then
     factor_een_gl = 0.0d0
     return
  end if

  do nw =1, walk_num
     factor_een_gl(:,:,nw) = 0.0d0
     do n = 1, dim_c_vector
        l = lkpm_combined_index(n, 1)
        k = lkpm_combined_index(n, 2)
        m = lkpm_combined_index(n, 4)

        do a = 1, nucl_num
           cn = c_vector_full(a, n)
           if(cn == 0.d0) cycle

           do ii = 1, 4
              do j = 1, elec_num
                 factor_een_gl(j,ii,nw) = factor_een_gl(j,ii,nw) + ( &
                      tmp_c (j,   a,m  ,k,nw) * een_rescaled_n_gl(j,ii,a,m+l,nw) + &
                      tmp_c (j,   a,m+l,k,nw) * een_rescaled_n_gl(j,ii,a,m  ,nw) + &
                      dtmp_c(j,ii,a,m  ,k,nw) * een_rescaled_n   (j,   a,m+l,nw) + &
                      dtmp_c(j,ii,a,m+l,k,nw) * een_rescaled_n   (j,   a,m  ,nw)   &
                      ) * cn
              end do
           end do

           cn = cn + cn
           do j = 1, elec_num
              factor_een_gl(j,4,nw) = factor_een_gl(j,4,nw) +  ( &
                   dtmp_c(j,1,a,m  ,k,nw) * een_rescaled_n_gl(j,1,a,m+l,nw)  + &
                   dtmp_c(j,2,a,m  ,k,nw) * een_rescaled_n_gl(j,2,a,m+l,nw)  + &
                   dtmp_c(j,3,a,m  ,k,nw) * een_rescaled_n_gl(j,3,a,m+l,nw)  + &
                   dtmp_c(j,1,a,m+l,k,nw) * een_rescaled_n_gl(j,1,a,m  ,nw)  + &
                   dtmp_c(j,2,a,m+l,k,nw) * een_rescaled_n_gl(j,2,a,m  ,nw)  + &
                   dtmp_c(j,3,a,m+l,k,nw) * een_rescaled_n_gl(j,3,a,m  ,nw)    &
                   ) * cn
           end do
        end do
     end do
  end do

end function qmckl_compute_jastrow_champ_factor_een_gl_doc

! Compute Gradient only
!      :PROPERTIES:
!      :Name:     qmckl_compute_jastrow_champ_factor_een_grad
!      :CRetType: qmckl_exit_code
!      :FRetType: qmckl_exit_code
!      :END:

!      #+NAME: qmckl_factor_een_grad_args
!      |-----------------------+---------------------------------------------------------------------+--------+------------------------------------------------|
!      | Variable              | Type                                                                | In/Out | Description                                    |
!      |-----------------------+---------------------------------------------------------------------+--------+------------------------------------------------|
!      | ~context~             | ~qmckl_context~                                                     | in     | Global state                                   |
!      | ~walk_num~            | ~int64_t~                                                           | in     | Number of walkers                              |
!      | ~elec_num~            | ~int64_t~                                                           | in     | Number of electrons                            |
!      | ~nucl_num~            | ~int64_t~                                                           | in     | Number of nuclei                               |
!      | ~cord_num~            | ~int64_t~                                                           | in     | order of polynomials                           |
!      | ~dim_c_vector~        | ~int64_t~                                                           | in     | dimension of full coefficient vector           |
!      | ~c_vector_full~       | ~double[dim_c_vector][nucl_num]~                                    | in     | full coefficient vector                        |
!      | ~lkpm_combined_index~ | ~int64_t[4][dim_c_vector]~                                          | in     | combined indices                               |
!      | ~tmp_c~               | ~double[walk_num][0:cord_num-1][0:cord_num][nucl_num][elec_num]~    | in     | Temporary intermediate tensor                  |
!      | ~dtmp_c~              | ~double[walk_num][0:cord_num-1][0:cord_num][nucl_num][4][elec_num]~ | in     | vector of non-zero coefficients                |
!      | ~een_rescaled_n~      | ~double[walk_num][0:cord_num][nucl_num][elec_num]~                  | in     | Electron-nucleus rescaled factor               |
!      | ~een_rescaled_n_gl~   | ~double[walk_num][0:cord_num][nucl_num][4][elec_num]~               | in     | Derivative of Electron-nucleus rescaled factor |
!      | ~factor_een_grad~     | ~double[walk_num][3][elec_num]~                                     | out    | Derivative of Electron-nucleus jastrow         |
!      |-----------------------+---------------------------------------------------------------------+--------+------------------------------------------------|



function qmckl_compute_jastrow_champ_factor_een_grad_doc( &
     context, walk_num, elec_num, nucl_num, &
     cord_num, dim_c_vector, c_vector_full, lkpm_combined_index, &
     tmp_c, dtmp_c, een_rescaled_n, een_rescaled_n_gl, factor_een_grad) &
     bind(C) result(info)

  use, intrinsic :: iso_c_binding
  use qmckl
  implicit none
  integer(qmckl_context), intent(in), value :: context
  integer (c_int64_t) , intent(in)  , value :: walk_num
  integer (c_int64_t) , intent(in)  , value :: elec_num
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  integer (c_int64_t) , intent(in)  , value :: cord_num
  integer (c_int64_t) , intent(in)  , value :: dim_c_vector
  real    (c_double ) , intent(in)          :: c_vector_full(nucl_num,dim_c_vector)
  integer (c_int64_t) , intent(in)          :: lkpm_combined_index(dim_c_vector,4)
  real    (c_double ) , intent(in)          :: tmp_c(elec_num,nucl_num,0:cord_num,0:cord_num-1,walk_num)
  real    (c_double ) , intent(in)          :: dtmp_c(elec_num,4,nucl_num,0:cord_num,0:cord_num-1,walk_num)
  real    (c_double ) , intent(in)          :: een_rescaled_n(elec_num,nucl_num,0:cord_num,walk_num)
  real    (c_double ) , intent(in)          :: een_rescaled_n_gl(elec_num,4,nucl_num,0:cord_num,walk_num)
  real    (c_double ) , intent(out)         :: factor_een_grad(elec_num,3,walk_num)
  integer(qmckl_exit_code)                  :: info

  integer*8 :: i, a, j, l, k, m, n, nw, ii
  double precision :: accu, accu2, cn

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) info = QMCKL_INVALID_CONTEXT
  if (walk_num <= 0)                 info = QMCKL_INVALID_ARG_2
  if (elec_num <= 0)                 info = QMCKL_INVALID_ARG_3
  if (nucl_num <= 0)                 info = QMCKL_INVALID_ARG_4
  if (cord_num <  0)                 info = QMCKL_INVALID_ARG_5
  if (info /= QMCKL_SUCCESS)         return


  if (cord_num == 0) then
     factor_een_grad = 0.0d0
     return
  end if

  do nw =1, walk_num
     factor_een_grad(:,:,nw) = 0.0d0
     do n = 1, dim_c_vector
        l = lkpm_combined_index(n, 1)
        k = lkpm_combined_index(n, 2)
        m = lkpm_combined_index(n, 4)

        do a = 1, nucl_num
           cn = c_vector_full(a, n)
           if(cn == 0.d0) cycle

           do ii = 1, 3
              do j = 1, elec_num
                 factor_een_grad(j,ii,nw) = factor_een_grad(j,ii,nw) + (           &
                      dtmp_c(j,ii,a,m  ,k,nw) * een_rescaled_n   (j,   a,m+l,nw) + &
                      dtmp_c(j,ii,a,m+l,k,nw) * een_rescaled_n   (j,   a,m  ,nw) + &
                      tmp_c(j,a,m  ,k,nw)     * een_rescaled_n_gl(j,ii,a,m+l,nw) + &
                      tmp_c(j,a,m+l,k,nw)     * een_rescaled_n_gl(j,ii,a,m  ,nw)   &
                      ) * cn
              end do
           end do

        end do
     end do
  end do

end function qmckl_compute_jastrow_champ_factor_een_grad_doc

! Compute
!      :PROPERTIES:
!      :Name:     qmckl_compute_jastrow_champ_value_doc
!      :CRetType: qmckl_exit_code
!      :FRetType: qmckl_exit_code
!      :END:

!      #+NAME: qmckl_jastrow_champ_value_args
!      |------------+--------------------+--------+----------------------|
!      | Variable   | Type               | In/Out | Description          |
!      |------------+--------------------+--------+----------------------|
!      | ~context~  | ~qmckl_context~    | in     | Global state         |
!      | ~walk_num~ | ~int64_t~          | in     | Number of walkers    |
!      | ~f_ee~     | ~double[walk_num]~ | in     | ee component         |
!      | ~f_en~     | ~double[walk_num]~ | in     | eN component         |
!      | ~f_een~    | ~double[walk_num]~ | in     | eeN component        |
!      | ~value~    | ~double[walk_num]~ | out    | Total Jastrow factor |
!      |------------+--------------------+--------+----------------------|


function qmckl_compute_jastrow_champ_value_doc(context, &
     walk_num, f_ee, f_en, f_een, value) &
     result(info) bind(C)
  use, intrinsic :: iso_c_binding
  use qmckl
  implicit none
  integer(qmckl_context), intent(in), value :: context
  integer (c_int64_t) , intent(in)  , value :: walk_num
  real    (c_double ) , intent(in)          :: f_ee(walk_num)
  real    (c_double ) , intent(in)          :: f_en(walk_num)
  real    (c_double ) , intent(in)          :: f_een(walk_num)
  real    (c_double ) , intent(out)         :: value(walk_num)
  integer(qmckl_exit_code)                  :: info

  integer*8 :: i

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  do i = 1, walk_num
     value(i) = f_ee(i) + f_en(i) + f_een(i)
  end do

  do i = 1, walk_num
     ! Flush to zero to avoid floating-point exception
     if (value(i) < -100.d0) then
       value(i) = 0.d0
     else
       value(i) = dexp(value(i))
     endif
  end do

end function qmckl_compute_jastrow_champ_value_doc

! Compute GL
!      :PROPERTIES:
!      :Name:     qmckl_compute_jastrow_champ_gl_doc
!      :CRetType: qmckl_exit_code
!      :FRetType: qmckl_exit_code
!      :END:

!      #+NAME: qmckl_jastrow_champ_gl_args
!      |------------+---------------------------------+--------+----------------------|
!      | Variable   | Type                            | In/Out | Description          |
!      |------------+---------------------------------+--------+----------------------|
!      | ~context~  | ~qmckl_context~                 | in     | Global state         |
!      | ~walk_num~ | ~int64_t~                       | in     | Number of walkers    |
!      | ~elec_num~ | ~int64_t~                       | in     | Number of electrons  |
!      | ~value~    | ~double[walk_num]~              | in     | Total Jastrow        |
!      | ~gl_ee~    | ~double[walk_num][4][elec_num]~ | in     | ee component         |
!      | ~gl_en~    | ~double[walk_num][4][elec_num]~ | in     | eN component         |
!      | ~gl_een~   | ~double[walk_num][4][elec_num]~ | in     | eeN component        |
!      | ~gl~       | ~double[walk_num][4][elec_num]~ | out    | Total Jastrow factor |
!      |------------+---------------------------------+--------+----------------------|


function qmckl_compute_jastrow_champ_gl_doc(context, &
     walk_num, elec_num, value, gl_ee, gl_en, gl_een, gl) &
     bind(C) result(info)
  use qmckl
  use, intrinsic :: iso_c_binding
  implicit none
  integer(qmckl_context), intent(in), value :: context
  integer (c_int64_t) , intent(in)  , value :: walk_num
  integer (c_int64_t) , intent(in)  , value :: elec_num
  real    (c_double ) , intent(in)          :: value(walk_num)
  real    (c_double ) , intent(in)          :: gl_ee(elec_num,4,walk_num)
  real    (c_double ) , intent(in)          :: gl_en(elec_num,4,walk_num)
  real    (c_double ) , intent(in)          :: gl_een(elec_num,4,walk_num)
  real    (c_double ) , intent(out)         :: gl(elec_num,4,walk_num)

  integer(qmckl_exit_code)                 :: info
  integer*8 :: i, j, k

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  do k = 1, walk_num
     do j=1,4
        do i = 1, elec_num
           gl(i,j,k) = gl_ee(i,j,k) + gl_en(i,j,k) + gl_een(i,j,k)
        end do
     end do
     do i = 1, elec_num
        gl(i,4,k) = gl(i,4,k) + &
             gl(i,1,k) * gl(i,1,k) + &
             gl(i,2,k) * gl(i,2,k) + &
             gl(i,3,k) * gl(i,3,k)
     end do
     gl(:,:,k) = gl(:,:,k) * value(k)
  end do


end function qmckl_compute_jastrow_champ_gl_doc

! Compute Gradient only
!      :PROPERTIES:
!      :Name:     qmckl_compute_jastrow_champ_grad_doc
!      :CRetType: qmckl_exit_code
!      :FRetType: qmckl_exit_code
!      :END:

!      #+NAME: qmckl_jastrow_champ_grad_args
!      |------------+---------------------------------+--------+----------------------|
!      | Variable   | Type                            | In/Out | Description          |
!      |------------+---------------------------------+--------+----------------------|
!      | ~context~  | ~qmckl_context~                 | in     | Global state         |
!      | ~walk_num~ | ~int64_t~                       | in     | Number of walkers    |
!      | ~elec_num~ | ~int64_t~                       | in     | Number of electrons  |
!      | ~value~    | ~double[walk_num]~              | in     | Total Jastrow        |
!      | ~gl_ee~    | ~double[walk_num][4][elec_num]~ | in     | ee component         |
!      | ~gl_en~    | ~double[walk_num][4][elec_num]~ | in     | eN component         |
!      | ~grad_een~ | ~double[walk_num][3][elec_num]~ | in     | eeN component        |
!      | ~grad~     | ~double[walk_num][3][elec_num]~ | out    | Total Jastrow factor |
!      |------------+---------------------------------+--------+----------------------|


function qmckl_compute_jastrow_champ_grad_doc(context, &
     walk_num, elec_num, value, gl_ee, gl_en, grad_een, grad) &
     result(info) bind(C)
  use, intrinsic :: iso_c_binding
  use qmckl
  implicit none
  integer(qmckl_context), intent(in), value :: context
  integer (c_int64_t) , intent(in)  , value :: walk_num
  integer (c_int64_t) , intent(in)  , value :: elec_num
  real    (c_double ) , intent(in)          :: value(walk_num)
  real    (c_double ) , intent(in)          :: gl_ee(elec_num,4,walk_num)
  real    (c_double ) , intent(in)          :: gl_en(elec_num,4,walk_num)
  real    (c_double ) , intent(in)          :: grad_een(elec_num,3,walk_num)
  real    (c_double ) , intent(out)         :: grad(elec_num,3,walk_num)

  integer(qmckl_exit_code) :: info
  integer*8 :: i, j, k

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  do k = 1, walk_num
     do j=1,3
        do i = 1, elec_num
           grad(i,j,k) = gl_ee(i,j,k) + gl_en(i,j,k) + grad_een(i,j,k)
        end do
     end do
     grad(:,:,k) = grad(:,:,k) * value(k)
  end do


end function qmckl_compute_jastrow_champ_grad_doc
