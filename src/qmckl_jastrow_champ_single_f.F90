! Compute
!    :PROPERTIES:
!    :Name:     qmckl_compute_single_ee_distance
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!     #+NAME: qmckl_single_ee_distance_args
!     |----------------------+---------------------------------+--------+-------------------------------------------------|
!     | Variable             | Type                            | In/Out | Description                                     |
!     |----------------------+---------------------------------+--------+-------------------------------------------------|
!     | ~context~            | ~qmckl_context~                 | in     | Global state                                    |
!     | ~num~                | ~int64_t~                       | in     | Index of single electron                        |
!     | ~elec_num~           | ~int64_t~                       | in     | Number of electrons                             |
!     | ~walk_num~           | ~int64_t~                       | in     | Number of walkers                               |
!     | ~coord~              | ~double[3][walk_num][elec_num]~ | in     | Electron coordinates                            |
!     | ~single_coord~       | ~double[walk_num][3]~           | in     | Single electron coordinates                     |
!     | ~single_ee_distance~ | ~double[walk_num][elec_num]~    | out    | Electron-electron distances for single electron |
!     |----------------------+---------------------------------+--------+-------------------------------------------------|


integer(qmckl_exit_code) function qmckl_compute_single_ee_distance(context, &
     num_in, elec_num, walk_num, coord, single_coord, single_ee_distance) &
     result(info) bind(C)
  use, intrinsic :: iso_c_binding
  use qmckl
  implicit none
  integer(qmckl_context), intent(in), value  :: context
  integer (c_int64_t) , intent(in)  , value :: elec_num, num_in
  integer (c_int64_t) , intent(in)  , value :: walk_num
  real    (c_double ) , intent(in)          :: coord(elec_num,walk_num,3)
  real    (c_double ) , intent(in)          :: single_coord(3,walk_num)
  real    (c_double ) , intent(out)         :: single_ee_distance(elec_num,walk_num)

  integer*8 :: k, i, j, num
  double precision :: x, y, z

  info = QMCKL_SUCCESS

  num = num_in + 1

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
     info = qmckl_distance(context, 'T', 'N', elec_num, 1_8, &
          coord(1,k,1), elec_num*walk_num, &
          single_coord(1,k), 3_8, &
          single_ee_distance(1,k), elec_num)
     if (info /= QMCKL_SUCCESS) then
        exit
     endif
     single_ee_distance(num,k) = 0.0d0
  end do



end function qmckl_compute_single_ee_distance

! Compute
!    :PROPERTIES:
!    :Name:     qmckl_compute_single_en_distance
!    :CRetType: qmckl_exit_code
!    :FRetType: qmckl_exit_code
!    :END:

!     #+NAME: qmckl_single_en_distance_args
!     |----------------------+------------------------------+--------+------------------------------------------------|
!     | Variable             | Type                         | In/Out | Description                                    |
!     |----------------------+------------------------------+--------+------------------------------------------------|
!     | ~context~            | ~qmckl_context~              | in     | Global state                                   |
!     | ~nucl_num~           | ~int64_t~                    | in     | Number of nuclei                               |
!     | ~walk_num~           | ~int64_t~                    | in     | Number of walkers                              |
!     | ~elec_coord~         | ~double[3][walk_num]~        | in     | Electron coordinates                           |
!     | ~nucl_coord~         | ~double[3][nucl_num]~        | in     | Nuclear coordinates                            |
!     | ~single_en_distance~ | ~double[walk_num][nucl_num]~ | out    | Electron-nucleus distances for single-electron |
!     |----------------------+------------------------------+--------+------------------------------------------------|


integer function qmckl_compute_single_en_distance(context, nucl_num, walk_num,  &
     elec_coord, nucl_coord, single_en_distance) result(info) bind(C)
  use, intrinsic :: iso_c_binding
  use qmckl
  implicit none
  integer(qmckl_context), intent(in), value :: context
  integer (c_int64_t) , intent(in)  , value :: nucl_num, walk_num
  real    (c_double ) , intent(in)          :: elec_coord(3,walk_num)
  real    (c_double ) , intent(in)          :: nucl_coord(nucl_num,3)
  real    (c_double ) , intent(out)         :: single_en_distance(nucl_num, walk_num)

  integer*8 :: k

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (nucl_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  info = qmckl_distance(context, 'T', 'N', nucl_num, walk_num, &
          nucl_coord, nucl_num, &
          elec_coord, 3_8, &
          single_en_distance, nucl_num)

end function qmckl_compute_single_en_distance

! Compute
!      :PROPERTIES:
!      :Name:     qmckl_compute_een_rescaled_single_e
!      :CRetType: qmckl_exit_code
!      :FRetType: qmckl_exit_code
!      :END:

!      #+NAME: qmckl_factor_een_rescaled_single_e_args
!      |-------------------------+----------------------------------------------------+--------+-------------------------------------------------------------|
!      | Variable                | Type                                               | In/Out | Description                                                 |
!      |-------------------------+----------------------------------------------------+--------+-------------------------------------------------------------|
!      | ~context~               | ~qmckl_context~                                    | in     | Global state                                                |
!      | ~num~                   | ~int64_t~                                          | in     | Number of single electron                                   |
!      | ~walk_num~              | ~int64_t~                                          | in     | Number of walkers                                           |
!      | ~elec_num~              | ~int64_t~                                          | in     | Number of electrons                                         |
!      | ~cord_num~              | ~int64_t~                                          | in     | Order of polynomials                                        |
!      | ~rescale_factor_ee~     | ~double~                                           | in     | Factor to rescale ee distances                              |
!      | ~single_ee_distance~    | ~double[walk_num][elec_num]~                       | in     | Single electron-electron distances for each walker          |
!      | ~een_rescaled_e~        | ~double[walk_num][0:cord_num][elec_num][elec_num]~ | in     | Rescaled electron-electron distances for each walker        |
!      | ~een_rescaled_single_e~ | ~double[walk_num][0:cord_num][elec_num]~           | out    | Single electron-electron rescaled distances for each walker |
!      |-------------------------+----------------------------------------------------+--------+-------------------------------------------------------------|


integer function qmckl_compute_een_rescaled_single_e_doc( &
     context, num_in, walk_num, elec_num, cord_num, rescale_factor_ee,  &
     single_ee_distance, een_rescaled_e, een_rescaled_single_e) &
     result(info) bind(C)
  use, intrinsic :: iso_c_binding
  use qmckl
  implicit none
  integer(qmckl_context), intent(in), value  :: context
  integer(c_int64_t)    , intent(in), value  :: num_in
  integer(c_int64_t)    , intent(in), value  :: walk_num
  integer(c_int64_t)    , intent(in), value  :: elec_num
  integer(c_int64_t)    , intent(in), value  :: cord_num
  real(c_double)        , intent(in), value  :: rescale_factor_ee
  real(c_double)        , intent(in)         :: single_ee_distance(elec_num,walk_num)
  real(c_double)        , intent(in)         :: een_rescaled_e(elec_num,elec_num,0:cord_num,walk_num)
  real(c_double)        , intent(out)        :: een_rescaled_single_e(elec_num,0:cord_num,walk_num)

  double precision,allocatable        :: een_rescaled_single_e_ij(:,:)
  double precision                    :: x
  integer*8                           :: i, j, k, l, nw, num

  num = num_in + 1
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

  allocate(een_rescaled_single_e_ij(elec_num, cord_num + 1))

  ! Prepare table of exponentiated distances raised to appropriate power
  do nw = 1, walk_num
     een_rescaled_single_e_ij(:, 1) = 1.0d0


     do j = 1, elec_num
        een_rescaled_single_e_ij(j, 2) = dexp(-rescale_factor_ee * single_ee_distance(j, nw))
     end do


     do l = 2, cord_num
        do k = 1, elec_num
           een_rescaled_single_e_ij(k, l + 1) = een_rescaled_single_e_ij(k, l) * een_rescaled_single_e_ij(k, 2)
        end do
     end do

     ! prepare the actual een table
     een_rescaled_single_e(:,0,nw) = 1.0d0

     do l = 1, cord_num
        do j = 1, elec_num
            x = een_rescaled_single_e_ij(j, l + 1)
            een_rescaled_single_e(j, l, nw) = x
        end do
     end do

     een_rescaled_single_e(num, :, :) = 0.0d0

  end do

end function qmckl_compute_een_rescaled_single_e_doc

! Compute
!      :PROPERTIES:
!      :Name:     qmckl_compute_een_rescaled_single_n
!      :CRetType: qmckl_exit_code
!      :FRetType: qmckl_exit_code
!      :END:

!      #+NAME: qmckl_factor_een_rescaled_single_n_args
!      | Variable                | Type                                               | In/Out | Description                                |
!      |-------------------------+----------------------------------------------------+--------+--------------------------------------------|
!      | ~context~               | ~qmckl_context~                                    | in     | Global state                               |
!      | ~num~                   | ~int64_t~                                          | in     | Number of single electron                  |
!      | ~walk_num~              | ~int64_t~                                          | in     | Number of walkers                          |
!      | ~elec_num~              | ~int64_t~                                          | in     | Number of atoms                            |
!      | ~nucl_num~              | ~int64_t~                                          | in     | Number of atoms                            |
!      | ~type_nucl_num~         | ~int64_t~                                          | in     | Number of atom types                       |
!      | ~type_nucl_vector~      | ~int64_t[nucl_num]~                                | in     | Types of atoms                             |
!      | ~cord_num~              | ~int64_t~                                          | in     | Order of polynomials                       |
!      | ~rescale_factor_en~     | ~double[nucl_num]~                                 | in     | Factor to rescale ee distances             |
!      | ~single_en_distance~    | ~double[walk_num][nucl_num]~                       | in     | Electron-nucleus distances                 |
!      | ~een_rescaled_n~        | ~double[walk_num][0:cord_num][nucl_num][elec_num]~ | in     | Electron-nucleus rescaled distances        |
!      | ~een_rescaled_single_n~ | ~double[walk_num][0:cord_num][nucl_num]~           | out    | Single electron-nucleus rescaled distances |
!      |-------------------------+----------------------------------------------------+--------+--------------------------------------------|


integer function qmckl_compute_een_rescaled_single_n( &
     context, num_in, walk_num, elec_num, nucl_num, &
     type_nucl_num, type_nucl_vector, cord_num, rescale_factor_en,  &
     single_en_distance, een_rescaled_n, een_rescaled_single_n) &
     result(info) bind(C)
  use, intrinsic :: iso_c_binding
  use qmckl
  implicit none
  integer(qmckl_context), intent(in), value  :: context
  integer(c_int64_t)    , intent(in), value  :: num_in
  integer(c_int64_t)    , intent(in), value  :: walk_num
  integer(c_int64_t)    , intent(in), value  :: elec_num
  integer(c_int64_t)    , intent(in), value  :: nucl_num
  integer(c_int64_t)    , intent(in), value  :: type_nucl_num
  integer(c_int64_t)    , intent(in)         :: type_nucl_vector(nucl_num)
  integer(c_int64_t)    , intent(in), value  :: cord_num
  real(c_double)        , intent(in)         :: rescale_factor_en(type_nucl_num)
  real(c_double)        , intent(in)         :: single_en_distance(nucl_num,walk_num)
  real(c_double)        , intent(in)         :: een_rescaled_n(elec_num,nucl_num,0:cord_num,walk_num)
  real(c_double)        , intent(out)        :: een_rescaled_single_n(nucl_num,0:cord_num,walk_num)

  double precision                    :: x
  integer*8                           :: i, a, k, l, nw, num

  num = num_in + 1

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (nucl_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  if (cord_num < 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  do nw = 1, walk_num

     ! prepare the actual een table
     een_rescaled_single_n(:, 0, nw) = 1.0d0

     do a = 1, nucl_num
        een_rescaled_single_n(a, 1, nw) = dexp(-rescale_factor_en(type_nucl_vector(a)+1) * single_en_distance(a, nw))
     end do

     do l = 2, cord_num
        do a = 1, nucl_num
           een_rescaled_single_n(a, l, nw) = een_rescaled_single_n(a, l - 1, nw) * een_rescaled_single_n(a, 1, nw)
        end do
     end do

  end do

end function qmckl_compute_een_rescaled_single_n

! Compute
!      :PROPERTIES:
!      :Name:     qmckl_compute_jastrow_champ_delta_p_doc
!      :CRetType: qmckl_exit_code
!      :FRetType: qmckl_exit_code
!      :END:

!      #+NAME: qmckl_factor_delta_p_args
!      | Variable                | Type                                                             | In/Out | Description                                 |
!      |-------------------------+------------------------------------------------------------------+--------+---------------------------------------------|
!      | ~context~               | ~qmckl_context~                                                  | in     | Global state                                |
!      | ~num~                   | ~int64_t~                                                        | in     | Single point index                          |
!      | ~walk_num~              | ~int64_t~                                                        | in     | Number of walkers                           |
!      | ~elec_num~              | ~int64_t~                                                        | in     | Number of electrons                         |
!      | ~nucl_num~              | ~int64_t~                                                        | in     | Number of nuclei                            |
!      | ~cord_num~              | ~int64_t~                                                        | in     | order of polynomials                        |
!      | ~een_rescaled_n~        | ~double[walk_num][0:cord_num][nucl_num][elec_num]~               | in     | Electron-nucleus rescaled distances         |
!      | ~een_rescaled_e~        | ~double[walk_num][0:cord_num][elec_num][elec_num]~               | in     | Electron-electron rescaled distances        |
!      | ~een_rescaled_single_n~ | ~double[walk_num][0:cord_num][nucl_num]~                         | in     | Electron-nucleus single rescaled distances  |
!      | ~een_rescaled_single_e~ | ~double[walk_num][0:cord_num][elec_num]~                         | in     | Electron-electron single rescaled distances |
!      | ~delta_p~               | ~double[walk_num][0:cord_num-1][0:cord_num][nucl_num][elec_num]~ | out    | Delta P matrix                              |
!      |-------------------------+------------------------------------------------------------------+--------+---------------------------------------------|


integer function qmckl_compute_jastrow_champ_delta_p_doc( &
     context, num_in, walk_num, elec_num, nucl_num, cord_num,   &
     een_rescaled_n, een_rescaled_e, een_rescaled_single_n, een_rescaled_single_e, delta_p) &
     result(info) bind(C)
  use, intrinsic :: iso_c_binding
  use qmckl
  implicit none
  integer(qmckl_context), intent(in), value  :: context
  integer(c_int64_t), intent(in), value  :: num_in, walk_num, elec_num, cord_num, nucl_num
  real(c_double) , intent(in)            :: een_rescaled_n(elec_num, nucl_num, 0:cord_num, walk_num)
  real(c_double) , intent(in)            :: een_rescaled_e(elec_num, elec_num, 0:cord_num, walk_num)
  real(c_double) , intent(in)            :: een_rescaled_single_n(nucl_num, 0:cord_num, walk_num)
  real(c_double) , intent(in)            :: een_rescaled_single_e(elec_num, 0:cord_num, walk_num)
  real(c_double) , intent(out)           :: delta_p(elec_num, nucl_num,0:cord_num, 0:cord_num-1,  walk_num)

  double precision        :: een_rescaled_delta_e(elec_num)

  integer*8 :: i, a, c, j, l, k, p, m, n, nw, num
  double precision :: dn, dn2
  integer*8                           :: LDA, LDB, LDC

  num = num_in + 1

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) info = QMCKL_INVALID_CONTEXT
  if (walk_num <= 0)                 info = QMCKL_INVALID_ARG_3
  if (elec_num <= 0)                 info = QMCKL_INVALID_ARG_4
  if (nucl_num <= 0)                 info = QMCKL_INVALID_ARG_5
  if (cord_num <  0)                 info = QMCKL_INVALID_ARG_6
  if (info /= QMCKL_SUCCESS)         return

  if (cord_num == 0) return

  do nw=1, walk_num

     do i=0, cord_num-1

        een_rescaled_delta_e(:) = een_rescaled_single_e(:,i,nw) - een_rescaled_e(:,num,i,nw)

        do c=0,cord_num
           do a=1,nucl_num
              dn = een_rescaled_single_n(a,c,nw) - een_rescaled_n(num,a,c,nw)
              dn2 = een_rescaled_single_n(a,c,nw)
              do j=1,elec_num
                 delta_p(j,a,c,i,nw) = een_rescaled_e(j,num,i,nw)*dn + een_rescaled_delta_e(j) * dn2
              enddo
           end do
        end do

        info = qmckl_dgemm(context, 'T', 'N', 1_8, nucl_num * (cord_num+1_8), elec_num, 1.0d0,     &
             een_rescaled_delta_e,elec_num, &
             een_rescaled_n(1,1,0,nw),elec_num,     &
             1.0d0,                                 &
             delta_p(num,1,0,i,nw),elec_num)

     enddo

  end do

end function qmckl_compute_jastrow_champ_delta_p_doc

! Compute
!      :PROPERTIES:
!      :Name:     qmckl_compute_jastrow_champ_factor_single_een_doc
!      :CRetType: qmckl_exit_code
!      :FRetType: qmckl_exit_code
!      :END:

!      #+NAME: qmckl_factor_single_een_args
!      | Variable                | Type                                                             | In/Out | Description                                       |
!      |-------------------------+------------------------------------------------------------------+--------+---------------------------------------------------|
!      | ~context~               | ~qmckl_context~                                                  | in     | Global state                                      |
!      | ~num~                   | ~int64_t~                                                        | in     | Single point number                               |
!      | ~walk_num~              | ~int64_t~                                                        | in     | Number of walkers                                 |
!      | ~elec_num~              | ~int64_t~                                                        | in     | Number of electrons                               |
!      | ~nucl_num~              | ~int64_t~                                                        | in     | Number of nuclei                                  |
!      | ~cord_num~              | ~int64_t~                                                        | in     | order of polynomials                              |
!      | ~dim_c_vector~          | ~int64_t~                                                        | in     | dimension of full coefficient vector              |
!      | ~c_vector_full~         | ~double[dim_c_vector][nucl_num]~                                 | in     | full coefficient vector                           |
!      | ~lkpm_combined_index~   | ~int64_t[4][dim_c_vector]~                                       | in     | combined indices                                  |
!      | ~tmp_c~                 | ~double[walk_num][0:cord_num-1][0:cord_num][nucl_num][elec_num]~ | in     | P matrix                                          |
!      | ~delta_p~               | ~double[walk_num][0:cord_num-1][0:cord_num][nucl_num][elec_num]~ | in     | Delta P matrix                                    |
!      | ~een_rescaled_n~        | ~double[walk_num][0:cord_num][nucl_num][elec_num]~               | in     | Electron-nucleus rescaled distances               |
!      | ~een_rescaled_e~        | ~double[walk_num][0:cord_num][elec_num][elec_num]~               | in     | Electron-electron rescaled distances              |
!      | ~een_rescaled_single_n~ | ~double[walk_num][0:cord_num][nucl_num]~                         | in     | Electron-nucleus single rescaled distances        |
!      | ~een_rescaled_single_e~ | ~double[walk_num][0:cord_num][elec_num]~                         | in     | Electron-electron single rescaled distances       |
!      | ~delta_een~             | ~double[walk_num]~                                               | out    | Single electron electron-electron-nucleus Jastrow |
!      |-------------------------+------------------------------------------------------------------+--------+---------------------------------------------------|


integer function qmckl_compute_jastrow_champ_factor_single_een_doc( &
     context, num_in, walk_num, elec_num, nucl_num, cord_num,   &
     dim_c_vector, c_vector_full, lkpm_combined_index, &
     tmp_c, delta_p, een_rescaled_n, een_rescaled_e, een_rescaled_single_n, &
     een_rescaled_single_e, delta_een) &
     result(info) bind(C)
  use, intrinsic :: iso_c_binding
  use qmckl
  implicit none
  integer(qmckl_context), intent(in), value  :: context
  integer(c_int64_t)    , intent(in), value  :: num_in, walk_num, elec_num, cord_num, nucl_num, dim_c_vector
  integer(c_int64_t)    , intent(in)  :: lkpm_combined_index(dim_c_vector,4)
  real(c_double)        , intent(in)  :: c_vector_full(nucl_num, dim_c_vector)
  real(c_double)        , intent(in)  :: tmp_c(elec_num, nucl_num,0:cord_num, 0:cord_num-1,  walk_num)
  real(c_double)        , intent(in)  :: delta_p(elec_num, nucl_num,0:cord_num, 0:cord_num-1,  walk_num)
  real(c_double)        , intent(in)  :: een_rescaled_n(elec_num, nucl_num, 0:cord_num, walk_num)
  real(c_double)        , intent(in)  :: een_rescaled_e(elec_num, elec_num, 0:cord_num, walk_num)
  real(c_double)        , intent(in)  :: een_rescaled_single_n(nucl_num, 0:cord_num, walk_num)
  real(c_double)        , intent(in)  :: een_rescaled_single_e(elec_num, 0:cord_num, walk_num)
  real(c_double)        , intent(out) :: delta_een(walk_num)

  double precision        :: delta_c(nucl_num,0:cord_num, 0:cord_num-1,  walk_num)
  double precision        :: delta_c2(elec_num, nucl_num,0:cord_num, 0:cord_num-1,  walk_num)

  double precision        :: een_rescaled_delta_n(nucl_num, 0:cord_num)

  integer*8 :: i, a, j, l, k, p, m, n, nw, num
  double precision :: accu, accu2, cn
  integer*8                           :: LDA, LDB, LDC

  num = num_in + 1

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) info = QMCKL_INVALID_CONTEXT
  if (walk_num <= 0)                 info = QMCKL_INVALID_ARG_3
  if (elec_num <= 0)                 info = QMCKL_INVALID_ARG_4
  if (nucl_num <= 0)                 info = QMCKL_INVALID_ARG_5
  if (cord_num <  0)                 info = QMCKL_INVALID_ARG_6
  if (info /= QMCKL_SUCCESS)         return

  delta_een = 0.0d0

  if (cord_num == 0) return

  do nw =1, walk_num
     een_rescaled_delta_n(:,:) = een_rescaled_single_n(:,:,nw) - een_rescaled_n(num,:,:,nw)
     do n = 1, dim_c_vector
        l = lkpm_combined_index(n, 1)
        k = lkpm_combined_index(n, 2)
        p = lkpm_combined_index(n, 3)
        m = lkpm_combined_index(n, 4)

        do a = 1, nucl_num
           cn = c_vector_full(a, n)
           if(cn == 0.d0) cycle

           accu = 0.0d0
           do j = 1, elec_num
              accu = accu + een_rescaled_n(j,a,m,nw) * delta_p(j,a,m+l,k,nw)
           end do
           accu = accu + een_rescaled_delta_n(a,m) * (tmp_c(num,a,m+l,k,nw) + delta_p(num,a,m+l,k,nw))
           delta_een(nw) = delta_een(nw) + accu * cn
        end do
     end do
  end do

end function qmckl_compute_jastrow_champ_factor_single_een_doc

! Compute
!      :PROPERTIES:
!      :Name:     qmckl_compute_een_rescaled_single_n_gl
!      :CRetType: qmckl_exit_code
!      :FRetType: qmckl_exit_code
!      :END:

!      #+NAME: qmckl_compute_een_rescaled_single_n_gl_args
!      |----------------------------+---------------------------------------------+--------+-------------------------------------------------------|
!      | Variable                   | Type                                        | In/Out | Description                                           |
!      |----------------------------+---------------------------------------------+--------+-------------------------------------------------------|
!      | ~context~                  | ~qmckl_context~                             | in     | Global state                                          |
!      | ~walk_num~                 | ~int64_t~                                   | in     | Number of walkers                                     |
!      | ~nucl_num~                 | ~int64_t~                                   | in     | Number of atoms                                       |
!      | ~type_nucl_num~            | ~int64_t~                                   | in     | Number of atom types                                  |
!      | ~type_nucl_vector~         | ~int64_t[nucl_num]~                         | in     | Types of atoms                                        |
!      | ~cord_num~                 | ~int64_t~                                   | in     | Order of polynomials                                  |
!      | ~rescale_factor_en~        | ~double[nucl_num]~                          | in     | Factor to rescale ee distances                        |
!      | ~coord_ee~                 | ~double[walk_num][3]~                       | in     | Electron coordinates                                  |
!      | ~coord_n~                  | ~double[3][nucl_num]~                       | in     | Nuclear coordinates                                   |
!      | ~single_en_distance~       | ~double[walk_num][nucl_num]~                | in     | Electron-nucleus single distances                     |
!      | ~een_rescaled_single_n~    | ~double[walk_num][0:cord_num][nucl_num]~    | in     | Electron-nucleus rescaled single distances            |
!      | ~een_rescaled_single_n_gl~ | ~double[walk_num][0:cord_num][nucl_num][4]~ | out    | Electron-nucleus rescaled single distances derivative |
!      |----------------------------+---------------------------------------------+--------+-------------------------------------------------------|


integer function qmckl_compute_een_rescaled_single_n_gl( &
     context, walk_num, nucl_num, type_nucl_num, type_nucl_vector, &
     cord_num, rescale_factor_en, coord_ee, coord_n, single_en_distance, &
     een_rescaled_single_n, een_rescaled_single_n_gl) &
     result(info) bind(C)
  use, intrinsic :: iso_c_binding
  use qmckl
  implicit none
  integer(qmckl_context), intent(in), value  :: context
  integer(c_int64_t)    , intent(in), value  :: walk_num
  integer(c_int64_t)    , intent(in), value  :: nucl_num
  integer(c_int64_t)    , intent(in), value  :: type_nucl_num
  integer(c_int64_t)    , intent(in)         :: type_nucl_vector(nucl_num)
  integer(c_int64_t)    , intent(in), value  :: cord_num
  real(c_double)        , intent(in)         :: rescale_factor_en(type_nucl_num)
  real(c_double)        , intent(in)         :: coord_ee(3,walk_num)
  real(c_double)        , intent(in)         :: coord_n(nucl_num,3)
  real(c_double)        , intent(in)         :: single_en_distance(nucl_num,walk_num)
  real(c_double)        , intent(in)         :: een_rescaled_single_n(nucl_num,0:cord_num,walk_num)
  real(c_double)        , intent(out)        :: een_rescaled_single_n_gl(4,nucl_num,0:cord_num,walk_num)

  double precision,allocatable   :: elnuc_dist_gl(:,:)
  double precision               :: x, ria_inv, kappa_l
  integer*8                      :: i, a, k, l, nw, ii

  allocate(elnuc_dist_gl(4, nucl_num))

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif


  if (nucl_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  if (cord_num < 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  ! Prepare table of exponentiated distances raised to appropriate power
  een_rescaled_single_n_gl             = 0.0d0
  do nw = 1, walk_num
     ! prepare the actual een table
     do a = 1, nucl_num
        ria_inv = 1.0d0 / single_en_distance(a, nw)
        do ii = 1, 3
           elnuc_dist_gl(ii, a) = (coord_ee(ii,nw) - coord_n(a, ii)) * ria_inv
        end do
        elnuc_dist_gl(4, a) = 2.0d0 * ria_inv
     end do

     do l = 0, cord_num
        do a = 1, nucl_num
           kappa_l = - dble(l) * rescale_factor_en(type_nucl_vector(a)+1)
           een_rescaled_single_n_gl(1, a, l, nw) = kappa_l * elnuc_dist_gl(1, a)
           een_rescaled_single_n_gl(2, a, l, nw) = kappa_l * elnuc_dist_gl(2, a)
           een_rescaled_single_n_gl(3, a, l, nw) = kappa_l * elnuc_dist_gl(3, a)
           een_rescaled_single_n_gl(4, a, l, nw) = kappa_l * (elnuc_dist_gl(4, a) + kappa_l)

           een_rescaled_single_n_gl(1, a, l, nw) = een_rescaled_single_n_gl(1, a, l, nw) * &
                een_rescaled_single_n(a, l, nw)
           een_rescaled_single_n_gl(2, a, l, nw) = een_rescaled_single_n_gl(2, a, l, nw) * &
                een_rescaled_single_n(a, l, nw)
           een_rescaled_single_n_gl(3, a, l, nw) = een_rescaled_single_n_gl(3, a, l, nw) * &
                een_rescaled_single_n(a, l, nw)
           een_rescaled_single_n_gl(4, a, l, nw) = een_rescaled_single_n_gl(4, a, l, nw) * &
                een_rescaled_single_n(a, l, nw)
        end do
     end do
  end do

end function qmckl_compute_een_rescaled_single_n_gl

! Compute
!      :PROPERTIES:
!      :Name:     qmckl_compute_een_rescaled_single_e_gl
!      :CRetType: qmckl_exit_code
!      :FRetType: qmckl_exit_code
!      :END:

!      #+NAME: qmckl_een_rescaled_single_e_gl_args
!      |----------------------------+---------------------------------------------+--------+--------------------------------------------------------|
!      | Variable                   | Type                                        | In/Out | Description                                            |
!      |----------------------------+---------------------------------------------+--------+--------------------------------------------------------|
!      | ~context~                  | ~qmckl_context~                             | in     | Global state                                           |
!      | ~num~                      | ~int64_t~                                   | in     | Index of single electron                               |
!      | ~walk_num~                 | ~int64_t~                                   | in     | Number of walkers                                      |
!      | ~elec_num~                 | ~int64_t~                                   | in     | Number of electrons                                    |
!      | ~cord_num~                 | ~int64_t~                                   | in     | Order of polynomials                                   |
!      | ~rescale_factor_ee~        | ~double~                                    | in     | Factor to rescale ee distances                         |
!      | ~coord~                    | ~double[walk_num][3]~                       | in     | Single electron coordinates                            |
!      | ~coord_ee~                 | ~double[3][walk_num][elec_num]~             | in     | Electron coordinates                                   |
!      | ~single_ee_distance~       | ~double[walk_num][elec_num]~                | in     | Electron-electron single distances                     |
!      | ~een_rescaled_single_e~    | ~double[walk_num][0:cord_num][elec_num]~    | in     | Electron-electron rescaled single distances            |
!      | ~een_rescaled_single_e_gl~ | ~double[walk_num][0:cord_num][elec_num][4]~ | out    | Electron-electron rescaled single distances derivative |
!      |----------------------------+---------------------------------------------+--------+--------------------------------------------------------|


integer function qmckl_compute_een_rescaled_single_e_gl_doc( &
     context, num_in, walk_num, elec_num, cord_num, rescale_factor_ee,  &
     coord, coord_ee, single_ee_distance, een_rescaled_single_e, een_rescaled_single_e_gl) &
     result(info) bind(C)
  use, intrinsic :: iso_c_binding
  use qmckl
  implicit none
  integer(qmckl_context), intent(in), value  :: context
  integer(c_int64_t)    , intent(in), value :: num_in
  integer(c_int64_t)    , intent(in), value :: walk_num
  integer(c_int64_t)    , intent(in), value :: elec_num
  integer(c_int64_t)    , intent(in), value :: cord_num
  real(c_double)        , intent(in), value :: rescale_factor_ee
  real(c_double)        , intent(in)        :: coord(3,walk_num)
  real(c_double)        , intent(in)        :: coord_ee(elec_num,walk_num,3)
  real(c_double)        , intent(in)        :: single_ee_distance(elec_num,walk_num)
  real(c_double)        , intent(in)        :: een_rescaled_single_e(elec_num,0:cord_num,walk_num)
  real(c_double)        , intent(out)       :: een_rescaled_single_e_gl(4,elec_num,0:cord_num,walk_num)

  double precision,allocatable   :: elec_dist_gl(:,:)
  double precision               :: x, rij_inv, kappa_l
  integer*8                      :: i, j, k, l, nw, ii, num

  num = num_in + 1

  allocate(elec_dist_gl(4, elec_num))

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  if (elec_num <= 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  if (cord_num < 0) then
     info = QMCKL_INVALID_ARG_5
     return
  endif

  !   Not necessary: should be set to zero by qmckl_malloc
  !   een_rescaled_single_e_gl = 0.0d0

  ! Prepare table of exponentiated distances raised to appropriate power

  do nw = 1, walk_num
     do i = 1, elec_num
        if (i == num) cycle
        rij_inv = 1.0d0 / single_ee_distance(i, nw)
        do ii = 1, 3
           elec_dist_gl(ii, i) = (coord(ii, nw) - coord_ee(i, nw, ii)) * rij_inv
        end do
        elec_dist_gl(4, i) = 2.0d0 * rij_inv
     end do

     elec_dist_gl(:, num) = 0.0d0

     do l = 1, cord_num
        kappa_l = - dble(l) * rescale_factor_ee
        do i = 1, elec_num
           een_rescaled_single_e_gl(1, i, l, nw) = kappa_l * elec_dist_gl(1, i)
           een_rescaled_single_e_gl(2, i, l, nw) = kappa_l * elec_dist_gl(2, i)
           een_rescaled_single_e_gl(3, i, l, nw) = kappa_l * elec_dist_gl(3, i)
           een_rescaled_single_e_gl(4, i, l, nw) = kappa_l * (elec_dist_gl(4, i) + kappa_l)

           een_rescaled_single_e_gl(1,i,l,nw) = een_rescaled_single_e_gl(1,i,l,nw) * een_rescaled_single_e(i,l,nw)
           een_rescaled_single_e_gl(2,i,l,nw) = een_rescaled_single_e_gl(2,i,l,nw) * een_rescaled_single_e(i,l,nw)
           een_rescaled_single_e_gl(3,i,l,nw) = een_rescaled_single_e_gl(3,i,l,nw) * een_rescaled_single_e(i,l,nw)
           een_rescaled_single_e_gl(4,i,l,nw) = een_rescaled_single_e_gl(4,i,l,nw) * een_rescaled_single_e(i,l,nw)

        end do
     end do
  end do



end function qmckl_compute_een_rescaled_single_e_gl_doc

! Compute
!      :PROPERTIES:
!      :Name:     qmckl_compute_jastrow_champ_delta_p_gl_doc
!      :CRetType: qmckl_exit_code
!      :FRetType: qmckl_exit_code
!      :END:

!      #+NAME: qmckl_factor_delta_p_gl_args
!      | Variable                   | Type                                                                | In/Out | Description                                             |
!      |----------------------------+---------------------------------------------------------------------+--------+---------------------------------------------------------|
!      | ~context~                  | ~qmckl_context~                                                     | in     | Global state                                            |
!      | ~num~                      | ~int64_t~                                                           | in     | Index of single electron                                |
!      | ~walk_num~                 | ~int64_t~                                                           | in     | Number of walkers                                       |
!      | ~elec_num~                 | ~int64_t~                                                           | in     | Number of electrons                                     |
!      | ~nucl_num~                 | ~int64_t~                                                           | in     | Number of nuclei                                        |
!      | ~cord_num~                 | ~int64_t~                                                           | in     | order of polynomials                                    |
!      | ~een_rescaled_n~           | ~double[walk_num][0:cord_num][nucl_num][elec_num]~                  | in     | Electron-nucleus rescaled distances                     |
!      | ~een_rescaled_e~           | ~double[walk_num][0:cord_num][elec_num][elec_num]~                  | in     | Electron-electron rescaled distances                    |
!      | ~een_rescaled_single_n~    | ~double[walk_num][0:cord_num][nucl_num]~                            | in     | Electron-nucleus single rescaled distances              |
!      | ~een_rescaled_single_e~    | ~double[walk_num][0:cord_num][elec_num]~                            | in     | Electron-electron single rescaled distances             |
!      | ~een_rescaled_n_gl~        | ~double[walk_num][0:cord_num][nucl_num][4][elec_num]~               | in     | Electron-nucleus rescaled distances derivatives         |
!      | ~een_rescaled_e_gl~        | ~double[walk_num][0:cord_num][elec_num][4][elec_num]~               | in     | Electron-electron rescaled distances derivatives        |
!      | ~een_rescaled_single_n_gl~ | ~double[walk_num][0:cord_num][nucl_num][4]~                         | in     | Electron-nucleus single rescaled distances derivatives  |
!      | ~een_rescaled_single_e_gl~ | ~double[walk_num][0:cord_num][elec_num][4]~                         | in     | Electron-electron single rescaled distances derivatives |
!      | ~delta_p_gl~               | ~double[walk_num][0:cord_num-1][0:cord_num][4][nucl_num][elec_num]~ | out    | Delta P matrix gradient and Laplacian                   |
!      |----------------------------+---------------------------------------------------------------------+--------+---------------------------------------------------------|


integer(qmckl_exit_code) function qmckl_compute_jastrow_champ_delta_p_gl_doc( &
     context, num_in, walk_num, elec_num, nucl_num, cord_num,   &
     een_rescaled_n, een_rescaled_e, een_rescaled_single_n, een_rescaled_single_e, &
     een_rescaled_n_gl, een_rescaled_e_gl, een_rescaled_single_n_gl, een_rescaled_single_e_gl, delta_p_gl) &
     result(info) bind(C)
  use, intrinsic :: iso_c_binding
  use qmckl
  implicit none
  integer(qmckl_context), intent(in), value  :: context
  integer(c_int64_t)    , intent(in), value  :: num_in, walk_num, elec_num, cord_num, nucl_num
  real(c_double)        , intent(in)  :: een_rescaled_n(elec_num, nucl_num, 0:cord_num, walk_num)
  real(c_double)        , intent(in)  :: een_rescaled_e(elec_num, elec_num, 0:cord_num, walk_num)
  real(c_double)        , intent(in)  :: een_rescaled_single_n(nucl_num, 0:cord_num, walk_num)
  real(c_double)        , intent(in)  :: een_rescaled_single_e(elec_num, 0:cord_num, walk_num)
  real(c_double)        , intent(in)  :: een_rescaled_n_gl(elec_num, 4, nucl_num, 0:cord_num, walk_num)
  real(c_double)        , intent(in)  :: een_rescaled_e_gl(elec_num, 4, elec_num, 0:cord_num, walk_num)
  real(c_double)        , intent(in)  :: een_rescaled_single_n_gl(4, nucl_num, 0:cord_num, walk_num)
  real(c_double)        , intent(in)  :: een_rescaled_single_e_gl(4,elec_num, 0:cord_num, walk_num)
  real(c_double)        , intent(out)  :: delta_p_gl(elec_num,nucl_num,4,0:cord_num, 0:cord_num-1,  walk_num)

  double precision        :: delta_e_gl(elec_num,4)

  double precision        :: een_rescaled_delta_n, een_re_n, een_re_single_n

  integer*8 :: i, a, j, l, k, p, m, n, nw, num
  double precision :: tmp, cummu
  integer*8        :: LDA, LDB, LDC

  num = num_in + 1

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) info = QMCKL_INVALID_CONTEXT
  if (walk_num <= 0)                 info = QMCKL_INVALID_ARG_3
  if (elec_num <= 0)                 info = QMCKL_INVALID_ARG_4
  if (nucl_num <= 0)                 info = QMCKL_INVALID_ARG_5
  if (cord_num <  0)                 info = QMCKL_INVALID_ARG_6
  if (info /= QMCKL_SUCCESS)         return


  if (cord_num == 0) then
     delta_p_gl = 0.d0
     return
  endif

  do nw=1, walk_num
     do m=1, cord_num-1
        do k = 1, 4
          do j = 1, elec_num
            delta_e_gl(j,k) = een_rescaled_single_e_gl(k,j,m,nw) - een_rescaled_e_gl(num, k, j, m, nw)
          end do
        end do
        do k = 1, 4
          delta_e_gl(num, k) = 0.0d0
        end do

        do l=0, cord_num
          do k = 1, 3
            do a = 1, nucl_num

              een_re_n = een_rescaled_n(num, a, l, nw)
              een_re_single_n = een_rescaled_single_n(a,l,nw)

              cummu = 0.0d0
              do i = 1, elec_num

                delta_p_gl(i,a,k,l,m,nw) = -een_rescaled_e_gl(i,k,num,m,nw) * een_re_n &
                  - een_rescaled_single_e_gl(k,i,m,nw) * een_re_single_n

                cummu = cummu + delta_e_gl(i,k) * een_rescaled_n(i,a,l,nw)
              end do
              delta_p_gl(num,a,k,l,m,nw) = delta_p_gl(num,a,k,l,m,nw) + cummu

            end do
          end do
          do a = 1, nucl_num
            een_rescaled_delta_n = een_rescaled_single_n(a,l,nw) - een_rescaled_n(num, a, l, nw)
            cummu = 0.0d0
            een_re_single_n = een_rescaled_single_n(a,l,nw)
            
            do i = 1, elec_num
              delta_p_gl(i,a,4,l,m,nw) = een_rescaled_e_gl(i,4,num,m,nw) * een_rescaled_delta_n &
                +delta_e_gl(i,4) * een_re_single_n

              cummu = cummu + delta_e_gl(i,4) * een_rescaled_n(i,a,l,nw)
            end do
            delta_p_gl(num,a,4,l,m,nw) = delta_p_gl(num,a,4,l,m,nw) + cummu

          end do
        end do
     end do
  end do

end function qmckl_compute_jastrow_champ_delta_p_gl_doc

! Compute
!      :PROPERTIES:
!      :Name:     qmckl_compute_jastrow_champ_factor_single_een_gl_doc
!      :CRetType: qmckl_exit_code
!      :FRetType: qmckl_exit_code
!      :END:

!      #+NAME: qmckl_factor_single_een_gl_args
!      | Variable                   | Type                                                                | In/Out | Description                                                    |
!      |----------------------------+---------------------------------------------------------------------+--------+----------------------------------------------------------------|
!      | ~context~                  | ~qmckl_context~                                                     | in     | Global state                                                   |
!      | ~num~                      | ~int64_t~                                                           | in     | Index of single electron                                       |
!      | ~walk_num~                 | ~int64_t~                                                           | in     | Number of walkers                                              |
!      | ~elec_num~                 | ~int64_t~                                                           | in     | Number of electrons                                            |
!      | ~nucl_num~                 | ~int64_t~                                                           | in     | Number of nuclei                                               |
!      | ~cord_num~                 | ~int64_t~                                                           | in     | order of polynomials                                           |
!      | ~dim_c_vector~             | ~int64_t~                                                           | in     | dimension of full coefficient vector                           |
!      | ~c_vector_full~            | ~double[dim_c_vector][nucl_num]~                                    | in     | full coefficient vector                                        |
!      | ~lkpm_combined_index~      | ~int64_t[4][dim_c_vector]~                                          | in     | combined indices                                               |
!      | ~tmp_c~                    | ~double[walk_num][0:cord_num-1][0:cord_num][nucl_num][elec_num]~    | in     | P matrix                                                       |
!      | ~dtmp_c~                   | ~double[walk_num][0:cord_num-1][0:cord_num][nucl_num][4][elec_num]~ | in     | P matrix derivative                                            |
!      | ~delta_p~                  | ~double[walk_num][0:cord_num-1][0:cord_num][nucl_num][elec_num]~    | in     | Delta P matrix                                                 |
!      | ~delta_p_gl~               | ~double[walk_num][0:cord_num-1][0:cord_num][4][nucl_num][elec_num]~ | in     | Delta P matrix derivative                                      |
!      | ~een_rescaled_n~           | ~double[walk_num][0:cord_num][nucl_num][elec_num]~                  | in     | Electron-nucleus rescaled distances                            |
!      | ~een_rescaled_single_n~    | ~double[walk_num][0:cord_num][nucl_num]~                            | in     | Electron-nucleus single rescaled distances                     |
!      | ~een_rescaled_n_gl~        | ~double[walk_num][0:cord_num][nucl_num][4][elec_num]~               | in     | Electron-nucleus rescaled distances derivatives                |
!      | ~een_rescaled_single_n_gl~ | ~double[walk_num][0:cord_num][nucl_num][4]~                         | in     | Electron-nucleus single rescaled distances derivatives         |
!      | ~delta_een_gl~             | ~double[walk_num][4][elec_num]~                                     | out    | Delta electron-electron-nucleus jastrow gradient and Laplacian |
!      |----------------------------+---------------------------------------------------------------------+--------+----------------------------------------------------------------|


integer(qmckl_exit_code) function qmckl_compute_jastrow_champ_factor_single_een_gl_doc( &
     context, num_in, walk_num, elec_num, nucl_num, cord_num,   &
     dim_c_vector, c_vector_full, lkpm_combined_index, &
     tmp_c, dtmp_c, delta_p, delta_p_gl, een_rescaled_n, een_rescaled_single_n, &
     een_rescaled_n_gl, een_rescaled_single_n_gl, delta_een_gl) &
     result(info) bind(C)
  use, intrinsic :: iso_c_binding
  use qmckl
  implicit none
  integer(qmckl_context), intent(in), value  :: context
  integer(c_int64_t)    , intent(in), value  :: num_in, walk_num, elec_num, cord_num, nucl_num, dim_c_vector
  integer(c_int64_t)    , intent(in)  :: lkpm_combined_index(dim_c_vector,4)
  real(c_double)        , intent(in)  :: c_vector_full(nucl_num, dim_c_vector)
  real(c_double)        , intent(in)  :: tmp_c(elec_num, nucl_num,0:cord_num, 0:cord_num-1,  walk_num)
  real(c_double)        , intent(in)  :: dtmp_c(elec_num, 4, nucl_num,0:cord_num, 0:cord_num-1,  walk_num)
  real(c_double)        , intent(in)  :: delta_p(elec_num, nucl_num,0:cord_num, 0:cord_num-1,  walk_num)
  real(c_double)        , intent(in)  :: delta_p_gl(elec_num, nucl_num, 4, 0:cord_num, 0:cord_num-1,  walk_num)
  real(c_double)        , intent(in)  :: een_rescaled_n(elec_num, nucl_num, 0:cord_num, walk_num)
  real(c_double)        , intent(in)  :: een_rescaled_single_n(nucl_num, 0:cord_num, walk_num)
  real(c_double)        , intent(in)  :: een_rescaled_n_gl(elec_num, 4, nucl_num, 0:cord_num, walk_num)
  real(c_double)        , intent(in)  :: een_rescaled_single_n_gl(4, nucl_num, 0:cord_num, walk_num)
  real(c_double)        , intent(out) :: delta_een_gl(elec_num, 4, walk_num)

  integer*8 :: i, a, j, l, k, p, m, n, nw, kk, num
  double precision :: accu, accu2, cn
  integer*8                           :: LDA, LDB, LDC

  double precision  :: een_rescaled_delta_n(nucl_num, 0:cord_num)
  double precision  :: een_rescaled_delta_n_gl(4, nucl_num, 0:cord_num)
  double precision :: dpg1_m, dpg1_ml, dp_m, dp_ml, een_r_m, een_r_ml, een_r_gl_m, een_r_gl_ml
  num = num_in + 1

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) info = QMCKL_INVALID_CONTEXT
  if (walk_num <= 0)                 info = QMCKL_INVALID_ARG_3
  if (elec_num <= 0)                 info = QMCKL_INVALID_ARG_4
  if (nucl_num <= 0)                 info = QMCKL_INVALID_ARG_5
  if (cord_num <  0)                 info = QMCKL_INVALID_ARG_6
  if (info /= QMCKL_SUCCESS)         return

  delta_een_gl = 0.0d0

  if (cord_num == 0) return


  do nw =1, walk_num

     een_rescaled_delta_n(:,:) = een_rescaled_single_n(:,:,nw) - een_rescaled_n(num, :, :, nw)
     een_rescaled_delta_n_gl(:,:,:) = een_rescaled_single_n_gl(:,:,:,nw) - een_rescaled_n_gl(num, :,:,:,nw)

     do n = 1, dim_c_vector
        l = lkpm_combined_index(n, 1)
        k = lkpm_combined_index(n, 2)
        p = lkpm_combined_index(n, 3)
        m = lkpm_combined_index(n, 4)

        do kk = 1, 4
           do a = 1, nucl_num
              cn = c_vector_full(a, n)
              if(cn == 0.d0) cycle
              do i = 1, elec_num
                 delta_een_gl(i,kk,nw) = delta_een_gl(i,kk,nw) + ( &
                      delta_p_gl(i,a,kk,m  ,k,nw) * een_rescaled_n(i,a,m+l,nw) + &
                      delta_p_gl(i,a,kk,m+l,k,nw) * een_rescaled_n(i,a,m  ,nw) + &
                      delta_p(i,a,m  ,k,nw) * een_rescaled_n_gl(i,kk,a,m+l,nw) + &
                      delta_p(i,a,m+l,k,nw) * een_rescaled_n_gl(i,kk,a,m  ,nw) ) * cn
              end do

              delta_een_gl(num,kk,nw) = delta_een_gl(num,kk,nw) + ( &
                   (dtmp_c(num,kk,a,m  ,k,nw) + delta_p_gl(num,a,kk,m  ,k,nw)) * een_rescaled_delta_n(a,m+l) + &
                   (dtmp_c(num,kk,a,m+l,k,nw) + delta_p_gl(num,a,kk,m+l,k,nw)) * een_rescaled_delta_n(a,m  ) + &
                   (tmp_c(num,a,m  ,k,nw) + delta_p(num,a,m  ,k,nw)) * een_rescaled_delta_n_gl(kk,a,m+l)  + &
                   (tmp_c(num,a,m+l,k,nw) + delta_p(num,a,m+l,k,nw)) * een_rescaled_delta_n_gl(kk,a,m  ) )* cn
           end do
        end do
        do a = 1, nucl_num
           cn = c_vector_full(a, n)
           if(cn == 0.d0) cycle
           cn = cn + cn
           do i = 1, elec_num
              delta_een_gl(i,4,nw) = delta_een_gl(i,4,nw) + ( &
                   delta_p_gl(i,a,1,m  ,k,nw) * een_rescaled_n_gl(i,1,a,m+l,nw) + &
                   delta_p_gl(i,a,1,m+l,k,nw) * een_rescaled_n_gl(i,1,a,m  ,nw) + &
                   delta_p_gl(i,a,2,m  ,k,nw) * een_rescaled_n_gl(i,2,a,m+l,nw) + &
                   delta_p_gl(i,a,2,m+l,k,nw) * een_rescaled_n_gl(i,2,a,m  ,nw) + &
                   delta_p_gl(i,a,3,m  ,k,nw) * een_rescaled_n_gl(i,3,a,m+l,nw) + &
                   delta_p_gl(i,a,3,m+l,k,nw) * een_rescaled_n_gl(i,3,a,m  ,nw) ) * cn
           end do
          delta_een_gl(num,4,nw) = delta_een_gl(num,4,nw) + ( &
               (delta_p_gl(num,a,1,m  ,k,nw) + dtmp_c(num,1,a,m  ,k,nw)) * een_rescaled_delta_n_gl(1,a,m+l) + &
               (delta_p_gl(num,a,1,m+l,k,nw) + dtmp_c(num,1,a,m+l,k,nw)) * een_rescaled_delta_n_gl(1,a,m  ) + &
               (delta_p_gl(num,a,2,m  ,k,nw) + dtmp_c(num,2,a,m  ,k,nw)) * een_rescaled_delta_n_gl(2,a,m+l) + &
               (delta_p_gl(num,a,2,m+l,k,nw) + dtmp_c(num,2,a,m+l,k,nw)) * een_rescaled_delta_n_gl(2,a,m  ) + &
               (delta_p_gl(num,a,3,m  ,k,nw) + dtmp_c(num,3,a,m  ,k,nw)) * een_rescaled_delta_n_gl(3,a,m+l) + &
               (delta_p_gl(num,a,3,m+l,k,nw) + dtmp_c(num,3,a,m+l,k,nw)) * een_rescaled_delta_n_gl(3,a,m  ) ) * cn
        end do
     end do
  end do

end function qmckl_compute_jastrow_champ_factor_single_een_gl_doc

! Compute
!      :PROPERTIES:
!      :Name:     qmckl_compute_jastrow_champ_delta_p_g_doc
!      :CRetType: qmckl_exit_code
!      :FRetType: qmckl_exit_code
!      :END:

!      #+NAME: qmckl_factor_delta_p_g_args
!      | Variable                   | Type                                                                | In/Out | Description                                             |
!      |----------------------------+---------------------------------------------------------------------+--------+---------------------------------------------------------|
!      | ~context~                  | ~qmckl_context~                                                     | in     | Global state                                            |
!      | ~num~                      | ~int64_t~                                                           | in     | Index of single electron                                |
!      | ~walk_num~                 | ~int64_t~                                                           | in     | Number of walkers                                       |
!      | ~elec_num~                 | ~int64_t~                                                           | in     | Number of electrons                                     |
!      | ~nucl_num~                 | ~int64_t~                                                           | in     | Number of nuclei                                        |
!      | ~cord_num~                 | ~int64_t~                                                           | in     | order of polynomials                                    |
!      | ~een_rescaled_n~           | ~double[walk_num][0:cord_num][nucl_num][elec_num]~                  | in     | Electron-nucleus rescaled distances                     |
!      | ~een_rescaled_e~           | ~double[walk_num][0:cord_num][elec_num][elec_num]~                  | in     | Electron-electron rescaled distances                    |
!      | ~een_rescaled_single_n~    | ~double[walk_num][0:cord_num][nucl_num]~                            | in     | Electron-nucleus single rescaled distances              |
!      | ~een_rescaled_single_e~    | ~double[walk_num][0:cord_num][elec_num]~                            | in     | Electron-electron single rescaled distances             |
!      | ~een_rescaled_n_gl~        | ~double[walk_num][0:cord_num][nucl_num][4][elec_num]~               | in     | Electron-nucleus rescaled distances derivatives         |
!      | ~een_rescaled_e_gl~        | ~double[walk_num][0:cord_num][elec_num][4][elec_num]~               | in     | Electron-electron rescaled distances derivatives        |
!      | ~een_rescaled_single_n_gl~ | ~double[walk_num][0:cord_num][nucl_num][4]~                         | in     | Electron-nucleus single rescaled distances derivatives  |
!      | ~een_rescaled_single_e_gl~ | ~double[walk_num][0:cord_num][elec_num][4]~                         | in     | Electron-electron single rescaled distances derivatives |
!      | ~delta_p_g~                | ~double[walk_num][0:cord_num-1][0:cord_num][4][nucl_num][elec_num]~ | out    | Delta P matrix gradient                                 |
!      |----------------------------+---------------------------------------------------------------------+--------+---------------------------------------------------------|


integer(qmckl_exit_code) function qmckl_compute_jastrow_champ_delta_p_g_doc( &
     context, num_in, walk_num, elec_num, nucl_num, cord_num,   &
     een_rescaled_n, een_rescaled_e, een_rescaled_single_n, een_rescaled_single_e, &
     een_rescaled_n_gl, een_rescaled_e_gl, een_rescaled_single_n_gl, een_rescaled_single_e_gl, delta_p_g) &
     result(info) bind(C)
  use, intrinsic :: iso_c_binding
  use qmckl
  implicit none
  integer(qmckl_context), intent(in), value  :: context
  integer(c_int64_t)    , intent(in), value  :: num_in, walk_num, elec_num, cord_num, nucl_num
  real(c_double)        , intent(in)  :: een_rescaled_n(elec_num, nucl_num, 0:cord_num, walk_num)
  real(c_double)        , intent(in)  :: een_rescaled_e(elec_num, elec_num, 0:cord_num, walk_num)
  real(c_double)        , intent(in)  :: een_rescaled_single_n(nucl_num, 0:cord_num, walk_num)
  real(c_double)        , intent(in)  :: een_rescaled_single_e(elec_num, 0:cord_num, walk_num)
  real(c_double)        , intent(in)  :: een_rescaled_n_gl(elec_num, 4, nucl_num, 0:cord_num, walk_num)
  real(c_double)        , intent(in)  :: een_rescaled_e_gl(elec_num, 4, elec_num, 0:cord_num, walk_num)
  real(c_double)        , intent(in)  :: een_rescaled_single_n_gl(4, nucl_num, 0:cord_num, walk_num)
  real(c_double)        , intent(in)  :: een_rescaled_single_e_gl(4,elec_num, 0:cord_num, walk_num)
  real(c_double)        , intent(out)  :: delta_p_g(elec_num,nucl_num,4,0:cord_num, 0:cord_num-1,  walk_num)

  double precision        :: delta_e_gl(elec_num,4)

  double precision        :: een_rescaled_delta_n, een_re_n, een_re_single_n

  integer*8 :: i, a, j, l, k, p, m, n, nw, num
  double precision :: tmp, cummu
  integer*8        :: LDA, LDB, LDC

  num = num_in + 1

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) info = QMCKL_INVALID_CONTEXT
  if (num <= 0)                      info = QMCKL_INVALID_ARG_2
  if (walk_num <= 0)                 info = QMCKL_INVALID_ARG_3
  if (elec_num <= 0)                 info = QMCKL_INVALID_ARG_4
  if (nucl_num <= 0)                 info = QMCKL_INVALID_ARG_5
  if (cord_num <  0)                 info = QMCKL_INVALID_ARG_6
  if (info /= QMCKL_SUCCESS)         return


  if (cord_num == 0) then
     delta_p_g = 0.d0
     return
  endif

  do nw=1, walk_num
     do m=1, cord_num-1
        do j = 1, elec_num
          do k = 1, 4
            delta_e_gl(j,k) = een_rescaled_single_e_gl(k,j,m,nw) - een_rescaled_e_gl(num, k, j, m, nw)
          end do
        end do
        do k = 1, 4
          delta_e_gl(num, k) = 0.0d0
        end do

        do l=0, cord_num
          do k = 1, 3
            do a = 1, nucl_num

              cummu = 0.0d0
              do i = 1, elec_num
                cummu = cummu + delta_e_gl(i,k) * een_rescaled_n(i,a,l,nw)
              end do
              delta_p_g(num,a,k,l,m,nw) = cummu


            end do
          end do
        end do
     end do
  end do

end function qmckl_compute_jastrow_champ_delta_p_g_doc

integer(qmckl_exit_code) function qmckl_compute_jastrow_champ_delta_p_g_hpc( &
     context, num_in, walk_num, elec_num, nucl_num, cord_num,   &
     een_rescaled_n, een_rescaled_e, een_rescaled_single_n, een_rescaled_single_e, &
     een_rescaled_n_gl, een_rescaled_e_gl, een_rescaled_single_n_gl, een_rescaled_single_e_gl, delta_p_g) &
     result(info) bind(C)
  use, intrinsic :: iso_c_binding
  use qmckl
  implicit none
  integer(qmckl_context), intent(in), value  :: context
  integer(c_int64_t)    , intent(in), value  :: num_in, walk_num, elec_num, cord_num, nucl_num
  real(c_double)        , intent(in)  :: een_rescaled_n(elec_num, nucl_num, 0:cord_num, walk_num)
  real(c_double)        , intent(in)  :: een_rescaled_e(elec_num, elec_num, 0:cord_num, walk_num)
  real(c_double)        , intent(in)  :: een_rescaled_single_n(nucl_num, 0:cord_num, walk_num)
  real(c_double)        , intent(in)  :: een_rescaled_single_e(elec_num, 0:cord_num, walk_num)
  real(c_double)        , intent(in)  :: een_rescaled_n_gl(elec_num, 4, nucl_num, 0:cord_num, walk_num)
  real(c_double)        , intent(in)  :: een_rescaled_e_gl(elec_num, 4, elec_num, 0:cord_num, walk_num)
  real(c_double)        , intent(in)  :: een_rescaled_single_n_gl(4, nucl_num, 0:cord_num, walk_num)
  real(c_double)        , intent(in)  :: een_rescaled_single_e_gl(4,elec_num, 0:cord_num, walk_num)
  real(c_double)        , intent(out)  :: delta_p_g(elec_num,nucl_num,4,0:cord_num, 0:cord_num-1,  walk_num)

  double precision        :: delta_e_gl(3,elec_num)

  double precision        :: een_rescaled_delta_n, een_re_n, een_re_single_n

  integer*8 :: i, a, j, l, k, p, m, n, nw, num
  double precision, allocatable :: tmp(:,:,:)
  integer*8        :: LDA, LDB, LDC

  num = num_in + 1

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) info = QMCKL_INVALID_CONTEXT
  if (num <= 0)                      info = QMCKL_INVALID_ARG_2
  if (walk_num <= 0)                 info = QMCKL_INVALID_ARG_3
  if (elec_num <= 0)                 info = QMCKL_INVALID_ARG_4
  if (nucl_num <= 0)                 info = QMCKL_INVALID_ARG_5
  if (cord_num <  0)                 info = QMCKL_INVALID_ARG_6
  if (info /= QMCKL_SUCCESS)         return


  if (cord_num == 0) then
     delta_p_g = 0.d0
     return
  endif

  allocate( tmp(3,nucl_num,0:cord_num) )
  do nw=1, walk_num
     do m=1, cord_num-1
        delta_e_gl(1:3,1:elec_num) = een_rescaled_single_e_gl(1:3,1:elec_num,m,nw) - &
             een_rescaled_e_gl(num, 1:3, 1:elec_num, m, nw)
        delta_e_gl(1:3,num) = 0.0d0

        call dgemm('N','N', 3, int(nucl_num*(cord_num+1),4), int(elec_num,4), 1.d0, &
             delta_e_gl(1,1), 3, een_rescaled_n(1,1,0,nw), int(elec_num,4), 0.d0, &
             tmp, 3)

        delta_p_g(num,1:nucl_num,1,0:cord_num,m,nw) = tmp(1,1:nucl_num,0:cord_num)
        delta_p_g(num,1:nucl_num,2,0:cord_num,m,nw) = tmp(2,1:nucl_num,0:cord_num)
        delta_p_g(num,1:nucl_num,3,0:cord_num,m,nw) = tmp(3,1:nucl_num,0:cord_num)
     end do
  end do
  deallocate(tmp)

end function qmckl_compute_jastrow_champ_delta_p_g_hpc

! Compute
!      :PROPERTIES:
!      :Name:     qmckl_compute_jastrow_champ_factor_single_een_g_doc
!      :CRetType: qmckl_exit_code
!      :FRetType: qmckl_exit_code
!      :END:

!      #+NAME: qmckl_factor_single_een_g_args
!      | Variable                   | Type                                                                | In/Out | Description                                                    |
!      |----------------------------+---------------------------------------------------------------------+--------+----------------------------------------------------------------|
!      | ~context~                  | ~qmckl_context~                                                     | in     | Global state                                                   |
!      | ~num~                      | ~int64_t~                                                           | in     | Index of single electron                                       |
!      | ~walk_num~                 | ~int64_t~                                                           | in     | Number of walkers                                              |
!      | ~elec_num~                 | ~int64_t~                                                           | in     | Number of electrons                                            |
!      | ~nucl_num~                 | ~int64_t~                                                           | in     | Number of nuclei                                               |
!      | ~cord_num~                 | ~int64_t~                                                           | in     | order of polynomials                                           |
!      | ~dim_c_vector~             | ~int64_t~                                                           | in     | dimension of full coefficient vector                           |
!      | ~c_vector_full~            | ~double[dim_c_vector][nucl_num]~                                    | in     | full coefficient vector                                        |
!      | ~lkpm_combined_index~      | ~int64_t[4][dim_c_vector]~                                          | in     | combined indices                                               |
!      | ~tmp_c~                    | ~double[walk_num][0:cord_num-1][0:cord_num][nucl_num][elec_num]~    | in     | P matrix                                                       |
!      | ~dtmp_c~                   | ~double[walk_num][0:cord_num-1][0:cord_num][nucl_num][4][elec_num]~ | in     | P matrix derivative                                            |
!      | ~delta_p~                  | ~double[walk_num][0:cord_num-1][0:cord_num][nucl_num][elec_num]~    | in     | Delta P matrix                                                 |
!      | ~delta_p_gl~               | ~double[walk_num][0:cord_num-1][0:cord_num][4][nucl_num][elec_num]~ | in     | Delta P matrix derivative                                      |
!      | ~een_rescaled_n~           | ~double[walk_num][0:cord_num][nucl_num][elec_num]~                  | in     | Electron-nucleus rescaled distances                            |
!      | ~een_rescaled_single_n~    | ~double[walk_num][0:cord_num][nucl_num]~                            | in     | Electron-nucleus single rescaled distances                     |
!      | ~een_rescaled_n_gl~        | ~double[walk_num][0:cord_num][nucl_num][4][elec_num]~               | in     | Electron-nucleus rescaled distances derivatives                |
!      | ~een_rescaled_single_n_gl~ | ~double[walk_num][0:cord_num][nucl_num][4]~                         | in     | Electron-nucleus single rescaled distances derivatives         |
!      | ~delta_een_g~              | ~double[walk_num][4][elec_num]~                                     | out    | Delta electron-electron-nucleus jastrow gradient               |
!      |----------------------------+---------------------------------------------------------------------+--------+----------------------------------------------------------------|


integer(qmckl_exit_code) function qmckl_compute_jastrow_champ_factor_single_een_g_doc( &
     context, num_in, walk_num, elec_num, nucl_num, cord_num,   &
     dim_c_vector, c_vector_full, lkpm_combined_index, &
     tmp_c, dtmp_c, delta_p, delta_p_gl, een_rescaled_n, een_rescaled_single_n, &
     een_rescaled_n_gl, een_rescaled_single_n_gl, delta_een_g) &
     result(info) bind(C)
  use, intrinsic :: iso_c_binding
  use qmckl
  implicit none
  integer(qmckl_context), intent(in), value  :: context
  integer(c_int64_t)    , intent(in), value  :: num_in, walk_num, elec_num, cord_num, nucl_num, dim_c_vector
  integer(c_int64_t)    , intent(in)  :: lkpm_combined_index(dim_c_vector,4)
  real(c_double)        , intent(in)  :: c_vector_full(nucl_num, dim_c_vector)
  real(c_double)        , intent(in)  :: tmp_c(elec_num, nucl_num,0:cord_num, 0:cord_num-1,  walk_num)
  real(c_double)        , intent(in)  :: dtmp_c(elec_num, 4, nucl_num,0:cord_num, 0:cord_num-1,  walk_num)
  real(c_double)        , intent(in)  :: delta_p(elec_num, nucl_num,0:cord_num, 0:cord_num-1,  walk_num)
  real(c_double)        , intent(in)  :: delta_p_gl(elec_num, nucl_num, 4, 0:cord_num, 0:cord_num-1,  walk_num)
  real(c_double)        , intent(in)  :: een_rescaled_n(elec_num, nucl_num, 0:cord_num, walk_num)
  real(c_double)        , intent(in)  :: een_rescaled_single_n(nucl_num, 0:cord_num, walk_num)
  real(c_double)        , intent(in)  :: een_rescaled_n_gl(elec_num, 4, nucl_num, 0:cord_num, walk_num)
  real(c_double)        , intent(in)  :: een_rescaled_single_n_gl(4, nucl_num, 0:cord_num, walk_num)
  real(c_double)        , intent(out) :: delta_een_g(elec_num, 4, walk_num)

  integer*8 :: i, a, j, l, k, p, m, n, nw, kk, num
  double precision :: accu, accu2, cn
  integer*8                           :: LDA, LDB, LDC

  double precision  :: een_rescaled_delta_n_gl(4, nucl_num, 0:cord_num, walk_num)
  double precision  :: een_rescaled_delta_n(nucl_num, 0:cord_num, walk_num)

  num = num_in + 1

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) info = QMCKL_INVALID_CONTEXT
  if (walk_num <= 0)                 info = QMCKL_INVALID_ARG_3
  if (elec_num <= 0)                 info = QMCKL_INVALID_ARG_4
  if (nucl_num <= 0)                 info = QMCKL_INVALID_ARG_5
  if (cord_num <  0)                 info = QMCKL_INVALID_ARG_6
  if (info /= QMCKL_SUCCESS)         return

  delta_een_g = 0.0d0

  if (cord_num == 0) return

  een_rescaled_delta_n(:,:,:) = een_rescaled_single_n(:,:,:) - een_rescaled_n(num, :, :, :)
  een_rescaled_delta_n_gl(:,:,:,:) = een_rescaled_single_n_gl(:,:,:,:) - een_rescaled_n_gl(num, :,:,:,:)


  do nw =1, walk_num
     do n = 1, dim_c_vector
        l = lkpm_combined_index(n, 1)
        k = lkpm_combined_index(n, 2)
        p = lkpm_combined_index(n, 3)
        m = lkpm_combined_index(n, 4)

        do kk = 1, 3
           do a = 1, nucl_num
              cn = c_vector_full(a, n)
              if(cn == 0.d0) cycle

              !   delta_een_g(num,kk,nw) = delta_een_g(num,kk,nw) + ( &
              !        delta_p_gl(num,a,kk,m  ,k,nw) * een_rescaled_n(num,a,m+l,nw) + &
              !        delta_p_gl(num,a,kk,m+l,k,nw) * een_rescaled_n(num,a,m  ,nw) + &
              !        delta_p(num,a,m  ,k,nw) * een_rescaled_n_gl(num,kk,a,m+l,nw) + &
              !        delta_p(num,a,m+l,k,nw) * een_rescaled_n_gl(num,kk,a,m  ,nw) ) * cn

              !delta_een_g(num,kk,nw) = delta_een_g(num,kk,nw) + ( &
              !     (dtmp_c(num,kk,a,m  ,k,nw) + delta_p_gl(num,a,kk,m  ,k,nw)) * een_rescaled_delta_n(a,m+l,nw) + &
              !     (dtmp_c(num,kk,a,m+l,k,nw) + delta_p_gl(num,a,kk,m+l,k,nw)) * een_rescaled_delta_n(a,m  ,nw) + &
              !     (tmp_c(num,a,m  ,k,nw) + delta_p(num,a,m  ,k,nw)) * een_rescaled_delta_n_gl(kk,a,m+l,nw)  + &
              !     (tmp_c(num,a,m+l,k,nw) + delta_p(num,a,m+l,k,nw)) * een_rescaled_delta_n_gl(kk,a,m  ,nw) )* cn

    

              delta_een_g(num,kk,nw) = delta_een_g(num,kk,nw) + ( &
                   dtmp_c(num,kk,a,m  ,k,nw) * een_rescaled_delta_n(a,m+l,nw) + &
                   dtmp_c(num,kk,a,m+l,k,nw) * een_rescaled_delta_n(a,m  ,nw) + &
                   tmp_c(num,a,m  ,k,nw) * een_rescaled_delta_n_gl(kk,a,m+l,nw)  + &
                   tmp_c(num,a,m+l,k,nw) * een_rescaled_delta_n_gl(kk,a,m  ,nw) + &
                   delta_p_gl(num,a,kk,m  ,k,nw) * een_rescaled_single_n(a,m+l,nw) + &
                   delta_p_gl(num,a,kk,m+l,k,nw) * een_rescaled_single_n(a,m  ,nw) + &
                   delta_p(num,a,m  ,k,nw) * een_rescaled_single_n_gl(kk,a,m+l,nw) + & 
                   delta_p(num,a,m+l,k,nw) * een_rescaled_single_n_gl(kk,a,m  ,nw) )* cn
           end do
        end do
     end do
  end do

end function qmckl_compute_jastrow_champ_factor_single_een_g_doc

! Compute
!      :PROPERTIES:
!      :Name:     qmckl_compute_ee_rescaled_single
!      :CRetType: qmckl_exit_code
!      :FRetType: qmckl_exit_code
!      :END:

!      #+NAME: qmckl_ee_rescaled_single_args
!      | Variable             | Type                         | In/Out | Description                          |
!      |----------------------+------------------------------+--------+--------------------------------------|
!      | ~context~            | ~qmckl_context~              | in     | Global state                         |
!      | ~elec_num~           | ~int64_t~                    | in     | Number of electrons                  |
!      | ~rescale_factor_ee~  | ~double~                     | in     | Factor to rescale ee distances       |
!      | ~walk_num~           | ~int64_t~                    | in     | Number of walkers                    |
!      | ~single_ee_distance~ | ~double[walk_num][elec_num]~ | in     | Single electron-electron distances   |
!      | ~ee_rescaled_single~ | ~double[walk_num][elec_num]~ | out    | Electron-electron rescaled distances |
!      |----------------------+------------------------------+--------+--------------------------------------|


function qmckl_compute_ee_rescaled_single_doc(context, &
     elec_num, rescale_factor_ee, walk_num, &
     single_ee_distance, ee_rescaled_single) &
     bind(C) result(info)
  use qmckl
  implicit none

  integer(qmckl_context), intent(in), value :: context
  integer (c_int64_t) , intent(in)  , value :: elec_num
  real    (c_double ) , intent(in)  , value :: rescale_factor_ee
  integer (c_int64_t) , intent(in)  , value :: walk_num
  real    (c_double ) , intent(in)          :: single_ee_distance(elec_num,walk_num)
  real    (c_double ) , intent(out)         :: ee_rescaled_single(elec_num,walk_num)
  integer(qmckl_exit_code)                  :: info

  integer*8 :: k, i
  real (c_double) :: inverse_rescale_factor_ee

  inverse_rescale_factor_ee = 1.0d0 / rescale_factor_ee

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
      do i=1,elec_num
     ee_rescaled_single(i,k) = (1.0d0 - dexp(-rescale_factor_ee * single_ee_distance(i,k))) * inverse_rescale_factor_ee
     enddo
  end do

end function qmckl_compute_ee_rescaled_single_doc

! Compute
!      :PROPERTIES:
!      :Name:     qmckl_compute_jastrow_champ_single_ee
!      :CRetType: qmckl_exit_code
!      :FRetType: qmckl_exit_code
!      :END:

!      #+NAME: qmckl_single_ee_args
!      | Variable               | Type                                   | In/Out | Description                                 |
!      |------------------------+----------------------------------------+--------+---------------------------------------------|
!      | ~context~              | ~qmckl_context~                        | in     | Global state                                |
!      | ~num~                  | ~int64_t~                              | in     | Index of single point                       |
!      | ~walk_num~             | ~int64_t~                              | in     | Number of walkers                           |
!      | ~elec_num~             | ~int64_t~                              | in     | Number of electrons                         |
!      | ~up_num~               | ~int64_t~                              | in     | Number of alpha electrons                   |
!      | ~bord_num~             | ~int64_t~                              | in     | Number of coefficients                      |
!      | ~b_vector~             | ~double[bord_num+1]~                   | in     | List of coefficients                        |
!      | ~ee_distance_rescaled~ | ~double[walk_num][elec_num][elec_num]~ | in     | Electron-electron rescaled distances        |
!      | ~ee_rescaled_single~   | ~double[walk_num][elec_num]~           | in     | Electron-electron rescaled single distances |
!      | ~delta_ee~             | ~double[walk_num]~                     | out    | Single electron-electron Jastrow            |
!      |------------------------+----------------------------------------+--------+---------------------------------------------|


function qmckl_compute_jastrow_champ_single_ee_doc(context, &
     num_in, walk_num, elec_num, up_num, bord_num, b_vector, &
     ee_distance_rescaled, ee_rescaled_single, spin_independent, delta_ee) &
     bind(C) result(info)
  use qmckl
  implicit none

  integer (qmckl_context), intent(in), value :: context
  integer (c_int64_t)    , intent(in), value :: num_in
  integer (c_int64_t)    , intent(in), value :: walk_num
  integer (c_int64_t)    , intent(in), value :: elec_num
  integer (c_int64_t)    , intent(in), value :: up_num
  integer (c_int64_t)    , intent(in), value :: bord_num
  real    (c_double )    , intent(in)        :: b_vector(bord_num+1)
  real    (c_double )    , intent(in)        :: ee_distance_rescaled(elec_num,elec_num,walk_num)
  real    (c_double )    , intent(in)        :: ee_rescaled_single(elec_num,walk_num)
  integer (c_int32_t)    , intent(in), value :: spin_independent
  real    (c_double )    , intent(out)       :: delta_ee(walk_num)
  integer(qmckl_exit_code)                   :: info

  integer*8 :: i, j, k, nw, num
  double precision   :: x, xk, y, yk

  num = num_in + 1

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  if (elec_num <= 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  if (bord_num < 0) then
     info = QMCKL_INVALID_ARG_5
     return
  endif

  do nw =1, walk_num

     delta_ee(nw) = 0.0d0
     do i=1,elec_num
         !print *,i, ee_rescaled_single(i,nw)
         !print *, i, ee_distance_rescaled(i,num,nw)
         !print *, '   '
        if (i.ne.num) then
           x = ee_distance_rescaled(i,num,nw)
           y = ee_rescaled_single(i,nw)
           if (spin_independent == 1) then
              delta_ee(nw) = delta_ee(nw) - (b_vector(1) * x / (1.d0 + b_vector(2) * x)) &
                   + (b_vector(1) * y / (1.d0 + b_vector(2) * y))
           else
              if ((i <= up_num .and. num <= up_num ) .or. (i >  up_num .and. num >  up_num)) then
                 delta_ee(nw) = delta_ee(nw) - (0.5d0 * b_vector(1) * x / (1.d0 + b_vector(2) * x)) &
                      + (0.5d0 * b_vector(1) * y / (1.d0 + b_vector(2) * y))
              else
                 delta_ee(nw) = delta_ee(nw) - (b_vector(1) * x / (1.d0 + b_vector(2) * x)) &
                      + (b_vector(1) * y / (1.d0 + b_vector(2) * y))
              endif
           endif

           xk = x
           yk = y
           do k=2,bord_num
              xk = xk * x
              yk = yk * y
              delta_ee(nw) = delta_ee(nw) - (b_vector(k+1) * xk) + (b_vector(k+1) * yk)
           end do
        endif
     end do

  end do

end function qmckl_compute_jastrow_champ_single_ee_doc

! Compute
!      :PROPERTIES:
!      :Name:     qmckl_compute_ee_rescaled_single_gl
!      :CRetType: qmckl_exit_code
!      :FRetType: qmckl_exit_code
!      :END:

!      #+NAME: qmckl_ee_rescaled_single_gl_args
!      | Variable                | Type                            | In/Out | Description                                            |
!      |-------------------------+---------------------------------+--------+--------------------------------------------------------|
!      | ~context~               | ~qmckl_context~                 | in     | Global state                                           |
!      | ~num~                   | ~int64_t~                       | in     | Index of single electron                               |
!      | ~elec_num~              | ~int64_t~                       | in     | Number of electrons                                    |
!      | ~rescale_factor_ee~     | ~double~                        | in     | Factor to rescale ee distances                         |
!      | ~walk_num~              | ~int64_t~                       | in     | Number of walkers                                      |
!      | ~single_ee_distance~    | ~double[elec_num][walk_num]~    | in     | Single electron-electron distances                     |
!      | ~elec_coord~            | ~double[3][walk_num][elec_num]~ | in     | Electron coordinates                                   |
!      | ~coord~                 | ~double[walk_num][3]~           | in     | Single electron coordinates                            |
!      | ~ee_rescaled_single_gl~ | ~double[walk_num][elec_num][4]~ | out    | Electron-electron rescaled single distance derivatives |
!      |-------------------------+---------------------------------+--------+--------------------------------------------------------|


function qmckl_compute_ee_rescaled_single_gl_doc(context, num_in,  &
     elec_num, rescale_factor_ee, walk_num, single_ee_distance, elec_coord, coord, ee_rescaled_single_gl) &
     bind(C) result(info)
  use qmckl
  implicit none

  integer(qmckl_context), intent(in), value :: context
   integer (c_int64_t) , intent(in)  , value :: num_in
  integer (c_int64_t) , intent(in)  , value :: elec_num
  real    (c_double ) , intent(in)  , value :: rescale_factor_ee
  integer (c_int64_t) , intent(in)  , value :: walk_num
  real    (c_double ) , intent(in)          :: single_ee_distance(elec_num,walk_num)
  real    (c_double ) , intent(in)          :: elec_coord(elec_num,walk_num,3)
  real    (c_double ) , intent(in)          :: coord(3,walk_num)
  real    (c_double ) , intent(out)         :: ee_rescaled_single_gl(4,elec_num,walk_num)
  integer(qmckl_exit_code)                  :: info

  integer*8 :: nw, i, ii, num
  double precision :: rij_inv, elel_dist_gl(4, elec_num), kappa_l

  num = num_in + 1

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


  ee_rescaled_single_gl = 0.0d0
  do nw = 1, walk_num

     ! prepare the actual een table
     do i = 1, elec_num
        rij_inv = 1.0d0 / single_ee_distance(i, nw)
        do ii = 1, 3
          elel_dist_gl(ii, i) = (elec_coord(i,nw, ii) - coord(ii,nw)) * rij_inv
        end do
        elel_dist_gl(4, i) = 2.0d0 * rij_inv
     end do

      do i = 1, elec_num
        kappa_l = -1 * rescale_factor_ee
        ee_rescaled_single_gl(1, i, nw) = elel_dist_gl(1, i)
        ee_rescaled_single_gl(2, i, nw) = elel_dist_gl(2, i)
        ee_rescaled_single_gl(3, i, nw) = elel_dist_gl(3, i)
        ee_rescaled_single_gl(4, i, nw) = elel_dist_gl(4, i)

        ee_rescaled_single_gl(4, i, nw) = ee_rescaled_single_gl(4, i, nw) + kappa_l

        ee_rescaled_single_gl(1, i, nw) = ee_rescaled_single_gl(1, i, nw)  * dexp(kappa_l * single_ee_distance(i,nw))
        ee_rescaled_single_gl(2, i, nw) = ee_rescaled_single_gl(2, i, nw)  * dexp(kappa_l * single_ee_distance(i,nw))
        ee_rescaled_single_gl(3, i, nw) = ee_rescaled_single_gl(3, i, nw)  * dexp(kappa_l * single_ee_distance(i,nw))
        ee_rescaled_single_gl(4, i, nw) = ee_rescaled_single_gl(4, i, nw)  * dexp(kappa_l * single_ee_distance(i,nw))
     end do

     ee_rescaled_single_gl(1, num, nw) = 0.0d0
     ee_rescaled_single_gl(2, num, nw) = 0.0d0
     ee_rescaled_single_gl(3, num, nw) = 0.0d0
     ee_rescaled_single_gl(4, num, nw) = 0.0d0
  end do


end function qmckl_compute_ee_rescaled_single_gl_doc

! Compute
!      :PROPERTIES:
!      :Name:     qmckl_compute_jastrow_champ_single_ee_gl
!      :CRetType: qmckl_exit_code
!      :FRetType: qmckl_exit_code
!      :END:

!      #+NAME: qmckl_single_ee_gl_args
!      | Variable                  | Type                                      | In/Out | Description                                               |
!      |---------------------------+-------------------------------------------+--------+-----------------------------------------------------------|
!      | ~context~                 | ~qmckl_context~                           | in     | Global state                                              |
!      | ~num~                     | ~int64_t~                                 | in     | Index of single electron                                  |
!      | ~walk_num~                | ~int64_t~                                 | in     | Number of walkers                                         |
!      | ~elec_num~                | ~int64_t~                                 | in     | Number of electrons                                       |
!      | ~up_num~                  | ~int64_t~                                 | in     | Number of alpha electrons                                 |
!      | ~bord_num~                | ~int64_t~                                 | in     | Number of coefficients                                    |
!      | ~b_vector~                | ~double[bord_num+1]~                      | in     | List of coefficients                                      |
!      | ~ee_distance_rescaled~    | ~double[walk_num][elec_num][elec_num]~    | in     | Electron-electron rescaled distances                      |
!      | ~ee_distance_rescaled_gl~ | ~double[walk_num][4][elec_num][elec_num]~ | in     | Electron-electron rescaled distances derivatives          |
!      | ~ee_rescaled_single~      | ~double[walk_num][elec_num]~              | in     | Electron-electron rescaled single distances               |
!      | ~ee_rescaled_single_gl~   | ~double[walk_num][4][elec_num]~           | in     | Electron-electron rescaled single distances derivatives   |
!      | ~spin_independent~        | ~int32_t~                                 | in     | If 1, same parameters for parallel and antiparallel spins |
!      | ~delta_ee_gl~             | ~double[walk_num][elec_num][4]~           | out    | Single electron-electron jastrow gradients and Laplacian  |
!      |---------------------------+-------------------------------------------+--------+-----------------------------------------------------------|


function qmckl_compute_jastrow_champ_single_ee_gl_doc( &
     context, num_in, walk_num, elec_num, up_num, bord_num, &
     b_vector, ee_distance_rescaled, ee_distance_rescaled_gl,  &
     ee_rescaled_single, ee_rescaled_single_gl,  &
     spin_independent, delta_ee_gl) &
     bind(C) result(info)
  use qmckl
  implicit none

  integer (qmckl_context), intent(in), value :: context
  integer (c_int64_t) , intent(in)  , value :: num_in
  integer (c_int64_t) , intent(in)  , value :: walk_num
  integer (c_int64_t) , intent(in)  , value :: elec_num
  integer (c_int64_t) , intent(in)  , value :: up_num
  integer (c_int64_t) , intent(in)  , value :: bord_num
  real    (c_double ) , intent(in)          :: b_vector(bord_num+1)
  real    (c_double ) , intent(in)          :: ee_distance_rescaled(elec_num,elec_num,walk_num)
  real    (c_double ) , intent(in)          :: ee_distance_rescaled_gl(4,elec_num,elec_num,walk_num)
  real    (c_double ) , intent(in)          :: ee_rescaled_single(elec_num,walk_num)
  real    (c_double ) , intent(in)          :: ee_rescaled_single_gl(4,elec_num,walk_num)
  integer (c_int32_t) , intent(in)  , value :: spin_independent
  real    (c_double ) , intent(out)         :: delta_ee_gl(4,elec_num,walk_num)
  integer(qmckl_exit_code)                  :: info

  integer*8 :: i, j, k, nw, ii, num
  double precision   :: x, x1, kf, x_old, x1_old
  double precision   :: denom, invdenom, invdenom2, f
  double precision   :: denom_old, invdenom_old, invdenom2_old, f_old
  double precision   :: grad_c2, grad_c2_old
  double precision   :: dx(4), dx_old(4)

  num = num_in + 1

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  if (elec_num <= 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  if (bord_num < 0) then
     info = QMCKL_INVALID_ARG_5
     return
  endif

  if ((spin_independent < 0).or.(spin_independent > 1)) then
     info = QMCKL_INVALID_ARG_8
     return
  endif

  do nw =1, walk_num
     delta_ee_gl(:,:,nw) = 0.0d0
        do i = 1, elec_num
           if (i == num) cycle

           x = ee_rescaled_single(i,nw)
           x_old = ee_distance_rescaled(i,num,nw)

           denom         = 1.0d0 + b_vector(2) * x
           invdenom      = 1.0d0 / denom
           invdenom2     = invdenom * invdenom

           denom_old         = 1.0d0 + b_vector(2) * x_old
           invdenom_old      = 1.0d0 / denom_old
           invdenom2_old     = invdenom_old * invdenom_old

           dx(1) = ee_rescaled_single_gl(1, i, nw)
           dx(2) = ee_rescaled_single_gl(2, i, nw)
           dx(3) = ee_rescaled_single_gl(3, i, nw)
           dx(4) = ee_rescaled_single_gl(4, i, nw)

           dx_old(1) = ee_distance_rescaled_gl(1, i, num, nw)
           dx_old(2) = ee_distance_rescaled_gl(2, i, num, nw)
           dx_old(3) = ee_distance_rescaled_gl(3, i, num, nw)
           dx_old(4) = ee_distance_rescaled_gl(4, i, num, nw)

           grad_c2 = dx(1)*dx(1) + dx(2)*dx(2) + dx(3)*dx(3)
           grad_c2_old = dx_old(1)*dx_old(1) + dx_old(2)*dx_old(2) + dx_old(3)*dx_old(3)

           if (spin_independent == 1) then
              f = b_vector(1) * invdenom2
              f_old = b_vector(1) * invdenom2_old
           else
              if((i <= up_num .and. num <= up_num ) .or. (i >  up_num .and. num >  up_num)) then
                 f = 0.5d0 * b_vector(1) * invdenom2
                 f_old = 0.5d0 * b_vector(1) * invdenom2_old
              else
                 f = b_vector(1) * invdenom2
                 f_old = b_vector(1) * invdenom2_old
              end if
           end if

           delta_ee_gl(1,i,nw) = delta_ee_gl(1,i,nw) + f * dx(1) - f_old * dx_old(1)
           delta_ee_gl(2,i,nw) = delta_ee_gl(2,i,nw) + f * dx(2) - f_old * dx_old(2)
           delta_ee_gl(3,i,nw) = delta_ee_gl(3,i,nw) + f * dx(3) - f_old * dx_old(3)
           delta_ee_gl(4,i,nw) = delta_ee_gl(4,i,nw) &
                + f * (dx(4) - 2.d0 * b_vector(2) * grad_c2 * invdenom) &
                - f_old * (dx_old(4) - 2.d0 * b_vector(2) * grad_c2_old * invdenom_old)

           delta_ee_gl(1,num,nw) = delta_ee_gl(1,num,nw) - f * dx(1) + f_old * dx_old(1)
           delta_ee_gl(2,num,nw) = delta_ee_gl(2,num,nw) - f * dx(2) + f_old * dx_old(2)
           delta_ee_gl(3,num,nw) = delta_ee_gl(3,num,nw) - f * dx(3) + f_old * dx_old(3)
           delta_ee_gl(4,num,nw) = delta_ee_gl(4,num,nw) &
                + f * (dx(4) - 2.d0 * b_vector(2) * grad_c2 * invdenom) &
                - f_old * (dx_old(4) - 2.d0 * b_vector(2) * grad_c2_old * invdenom_old)


           kf = 2.d0
           x1 = x
           x1_old = x_old
           x = 1.d0
           x_old = 1.d0
           do k=2, bord_num
              f = b_vector(k+1) * kf * x
              f_old = b_vector(k+1) * kf * x_old
              delta_ee_gl(1,i,nw) = delta_ee_gl(1,i,nw) + f * x1 * dx(1) - f_old * x1_old * dx_old(1)
              delta_ee_gl(2,i,nw) = delta_ee_gl(2,i,nw) + f * x1 * dx(2) - f_old * x1_old * dx_old(2)
              delta_ee_gl(3,i,nw) = delta_ee_gl(3,i,nw) + f * x1 * dx(3) - f_old * x1_old * dx_old(3)
              delta_ee_gl(4,i,nw) = delta_ee_gl(4,i,nw) &
                   + f * (x1 * dx(4) + (kf-1.d0) * grad_c2) &
                   - f_old * (x1_old * dx_old(4) + (kf-1.d0) * grad_c2_old)

              delta_ee_gl(1,num,nw) = delta_ee_gl(1,num,nw) - f * x1 * dx(1) + f_old * x1_old * dx_old(1)
              delta_ee_gl(2,num,nw) = delta_ee_gl(2,num,nw) - f * x1 * dx(2) + f_old * x1_old * dx_old(2)
              delta_ee_gl(3,num,nw) = delta_ee_gl(3,num,nw) - f * x1 * dx(3) + f_old * x1_old * dx_old(3)
              delta_ee_gl(4,num,nw) = delta_ee_gl(4,num,nw) &
                   + f * (x1 * dx(4) + (kf-1.d0) * grad_c2) &
                   - f_old * (x1_old * dx_old(4) + (kf-1.d0) * grad_c2_old)
              x = x*x1
              x_old = x_old*x1_old
              kf = kf + 1.d0
           end do


        end do

  end do

end function qmckl_compute_jastrow_champ_single_ee_gl_doc

! Compute
!      :PROPERTIES:
!      :Name:     qmckl_compute_en_rescaled_single
!      :CRetType: qmckl_exit_code
!      :FRetType: qmckl_exit_code
!      :END:

!      #+NAME: qmckl_en_rescaled_single_args
!      | Variable             | Type                         | In/Out | Description                         |
!      |----------------------+------------------------------+--------+-------------------------------------|
!      | ~context~            | ~qmckl_context~              | in     | Global state                        |
!      | ~elec_num~           | ~int64_t~                    | in     | Number of electrons                 |
!      | ~nucl_num~           | ~int64_t~                    | in     | Number of nuclei                    |
!      | ~type_nucl_num~      | ~int64_t~                    | in     | Number of types of nuclei           |
!      | ~type_nucl_vector~   | ~int64_t[nucl_num]~          | in     | Number of types of nuclei           |
!      | ~rescale_factor_en~  | ~double[type_nucl_num]~      | in     | The factor for rescaled distances   |
!      | ~walk_num~           | ~int64_t~                    | in     | Number of walkers                   |
!      | ~single_en_distance~ | ~double[walk_num][nucl_num]~ | in     | Single electron-nucleus distances   |
!      | ~en_rescaled_single~ | ~double[walk_num][nucl_num]~ | out    | Electron-nucleus rescaled distances |
!      |----------------------+------------------------------+--------+-------------------------------------|


function qmckl_compute_en_rescaled_single_doc(context, &
     nucl_num, type_nucl_num, type_nucl_vector, rescale_factor_en, &
     walk_num, single_en_distance, en_rescaled_single) &
     bind(C) result(info)
  use qmckl
  implicit none
  integer (qmckl_context), intent(in), value :: context
  integer (c_int64_t) , intent(in)  , value  :: nucl_num
  integer (c_int64_t) , intent(in)  , value  :: type_nucl_num
  integer (c_int64_t) , intent(in)           :: type_nucl_vector(nucl_num)
  real    (c_double ) , intent(in)           :: rescale_factor_en(type_nucl_num)
  integer (c_int64_t) , intent(in)  , value  :: walk_num
  real    (c_double ) , intent(in)           :: single_en_distance(nucl_num,walk_num)
  real    (c_double ) , intent(out)          :: en_rescaled_single(nucl_num,walk_num)
  integer(qmckl_exit_code)                   :: info

  integer*8 :: i, k
  double precision      :: coord(3)

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif


  if (nucl_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  do i=1, nucl_num
     do k=1,walk_num
        en_rescaled_single(i,k) = (1.0d0 - dexp(-rescale_factor_en(type_nucl_vector(i)+1) * &
             single_en_distance(i,k))) / rescale_factor_en(type_nucl_vector(i)+1)
     end do
  end do

end function qmckl_compute_en_rescaled_single_doc

! Compute
!      :PROPERTIES:
!      :Name:     qmckl_compute_jastrow_champ_single_en_doc
!      :CRetType: qmckl_exit_code
!      :FRetType: qmckl_exit_code
!      :END:

!      #+NAME: qmckl_single_en_args
!      |------------------------+----------------------------------------+--------+--------------------------------------------|
!      | Variable               | Type                                   | In/Out | Description                                |
!      |------------------------+----------------------------------------+--------+--------------------------------------------|
!      | ~context~              | ~qmckl_context~                        | in     | Global state                               |
!      | ~num~                  | ~int64_t~                              | in     | Index of single point                      |
!      | ~walk_num~             | ~int64_t~                              | in     | Number of walkers                          |
!      | ~elec_num~             | ~int64_t~                              | in     | Number of electrons                        |
!      | ~nucl_num~             | ~int64_t~                              | in     | Number of nuclei                           |
!      | ~type_nucl_num~        | ~int64_t~                              | in     | Number of unique nuclei                    |
!      | ~type_nucl_vector~     | ~int64_t[nucl_num]~                    | in     | IDs of unique nuclei                       |
!      | ~aord_num~             | ~int64_t~                              | in     | Number of coefficients                     |
!      | ~a_vector~             | ~double[type_nucl_num][aord_num+1]~    | in     | List of coefficients                       |
!      | ~en_distance_rescaled~ | ~double[walk_num][nucl_num][elec_num]~ | in     | Electron-nucleus rescaled distances        |
!      | ~en_rescaled_single ~  | ~double[walk_num][nucl_num]~           | in     | Electron-nucleus rescaled single distances |
!      | ~delta_en~             | ~double[walk_num]~                     | out    | Single electron-nucleus jastrow            |
!      |------------------------+----------------------------------------+--------+--------------------------------------------|


function qmckl_compute_jastrow_champ_single_en_doc( &
     context, num_in, walk_num, elec_num, nucl_num, type_nucl_num, &
     type_nucl_vector, aord_num, a_vector, &
     en_distance_rescaled, en_rescaled_single, delta_en) &
     bind(C) result(info)
  use qmckl
  implicit none

  integer (qmckl_context), intent(in), value :: context
  integer (c_int64_t) , intent(in)  , value :: num_in
  integer (c_int64_t) , intent(in)  , value :: walk_num
  integer (c_int64_t) , intent(in)  , value :: elec_num
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  integer (c_int64_t) , intent(in)  , value :: type_nucl_num
  integer (c_int64_t) , intent(in)          :: type_nucl_vector(nucl_num)
  integer (c_int64_t) , intent(in)  , value :: aord_num
  real    (c_double ) , intent(in)          :: a_vector(aord_num+1,type_nucl_num)
  real    (c_double ) , intent(in)          :: en_distance_rescaled(elec_num,nucl_num,walk_num)
  real    (c_double ) , intent(in)          :: en_rescaled_single(nucl_num,walk_num)
  real    (c_double ) , intent(out)         :: delta_en(walk_num)
  integer(qmckl_exit_code)                  :: info

  integer*8 :: i, a, p, nw, num
  double precision   :: x, power_ser, y

  num = num_in + 1

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
     delta_en(nw) = 0.0d0
     do a = 1, nucl_num
          x = en_distance_rescaled(num, a, nw)
          y = en_rescaled_single(a, nw)

          delta_en(nw) = delta_en(nw) - a_vector(1, type_nucl_vector(a)+1) * x / (1.0d0 + a_vector(2, type_nucl_vector(a)+1) * x)
          delta_en(nw) = delta_en(nw) + a_vector(1, type_nucl_vector(a)+1) * y / (1.0d0 + a_vector(2, type_nucl_vector(a)+1) * y)

          do p = 2, aord_num
            x = x * en_distance_rescaled(num, a, nw)
            y = y * en_rescaled_single(a, nw)
            delta_en(nw) = delta_en(nw) - a_vector(p + 1, type_nucl_vector(a)+1) * x + a_vector(p + 1, type_nucl_vector(a)+1) * y
          end do

     end do
  end do

end function qmckl_compute_jastrow_champ_single_en_doc

! Compute
!      :PROPERTIES:
!      :Name:     qmckl_compute_en_rescaled_single_gl
!      :CRetType: qmckl_exit_code
!      :FRetType: qmckl_exit_code
!      :END:

!      #+NAME: qmckl_en_rescaled_single_gl_args
!      | Variable                | Type                            | In/Out | Description                                           |
!      |-------------------------+---------------------------------+--------+-------------------------------------------------------|
!      | ~context~               | ~qmckl_context~                 | in     | Global state                                          |
!      | ~nucl_num~              | ~int64_t~                       | in     | Number of nuclei                                      |
!      | ~type_nucl_num~         | ~int64_t~                       | in     | Number of nucleus types                               |
!      | ~type_nucl_vector~      | ~int64_t[nucl_num]~             | in     | Array of nucleus types                                |
!      | ~rescale_factor_en~     | ~double[nucl_num]~              | in     | The factors for rescaled distances                    |
!      | ~walk_num~              | ~int64_t~                       | in     | Number of walkers                                     |
!      | ~single_en_distance~    | ~double[walk_num][nucl_num]~    | in     | Single electorn distances                             |
!      | ~coord~                 | ~double[walk_num][3]~           | in     | Single electron coordinates                           |
!      | ~nucl_coord~            | ~double[3][nucl_num]~           | in     | Nucleus coordinates                                   |
!      | ~en_rescaled_single_gl~ | ~double[walk_num][nucl_num][4]~ | out    | Electron-nucleus rescaled single distance derivatives |
!      |-------------------------+---------------------------------+--------+-------------------------------------------------------|


integer function qmckl_compute_en_rescaled_single_gl_doc(context, nucl_num, &
     type_nucl_num, type_nucl_vector, rescale_factor_en, walk_num, &
     single_en_distance, coord, nucl_coord, en_rescaled_single_gl) &
     result(info) bind(C)
  use qmckl
  implicit none
  integer(qmckl_context), intent(in), value :: context
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  integer (c_int64_t) , intent(in)  , value :: type_nucl_num
  integer (c_int64_t) , intent(in)          :: type_nucl_vector(nucl_num)
  real    (c_double ) , intent(in)          :: rescale_factor_en(nucl_num)
  integer (c_int64_t) , intent(in)  , value :: walk_num
  real    (c_double ) , intent(in)          :: single_en_distance(nucl_num,walk_num)
  real    (c_double ) , intent(in)          :: coord(3,walk_num)
  real    (c_double ) , intent(in)          :: nucl_coord(nucl_num,3)
  real    (c_double ) , intent(out)         :: en_rescaled_single_gl(4,nucl_num,walk_num)

  integer*8 :: nw, a, ii
  double precision :: ria_inv, elnuc_dist_gl(4, nucl_num), kappa_l

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif


  if (nucl_num <= 0) then
     info = QMCKL_INVALID_ARG_2
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif


  en_rescaled_single_gl = 0.0d0
  do nw = 1, walk_num

     ! prepare the actual een table
     do a = 1, nucl_num
        ria_inv = 1.0d0 / single_en_distance(a, nw)
        do ii = 1, 3
          elnuc_dist_gl(ii, a) = (coord(ii,nw) - nucl_coord(a, ii)) * ria_inv
        end do
        elnuc_dist_gl(4, a) = 2.0d0 * ria_inv
     end do

      do a = 1, nucl_num
        kappa_l = -1 * rescale_factor_en(type_nucl_vector(a)+1)
        en_rescaled_single_gl(1, a, nw) = elnuc_dist_gl(1, a)
        en_rescaled_single_gl(2, a, nw) = elnuc_dist_gl(2, a)
        en_rescaled_single_gl(3, a, nw) = elnuc_dist_gl(3, a)
        en_rescaled_single_gl(4, a, nw) = elnuc_dist_gl(4, a)

        en_rescaled_single_gl(4, a, nw) = en_rescaled_single_gl(4, a, nw) + kappa_l

        en_rescaled_single_gl(1, a, nw) = en_rescaled_single_gl(1, a, nw)  * dexp(kappa_l * single_en_distance(a,nw))
        en_rescaled_single_gl(2, a, nw) = en_rescaled_single_gl(2, a, nw)  * dexp(kappa_l * single_en_distance(a,nw))
        en_rescaled_single_gl(3, a, nw) = en_rescaled_single_gl(3, a, nw)  * dexp(kappa_l * single_en_distance(a,nw))
        en_rescaled_single_gl(4, a, nw) = en_rescaled_single_gl(4, a, nw)  * dexp(kappa_l * single_en_distance(a,nw))
     end do

  end do


end function qmckl_compute_en_rescaled_single_gl_doc

! Compute
!      :PROPERTIES:
!      :Name:     qmckl_compute_jastrow_champ_single_en_gl
!      :CRetType: qmckl_exit_code
!      :FRetType: qmckl_exit_code
!      :END:

!      #+NAME: qmckl_single_en_gl_args
!      | Variable                  | Type                                      | In/Out | Description                                             |
!      |---------------------------+-------------------------------------------+--------+---------------------------------------------------------|
!      | ~context~                 | ~qmckl_context~                           | in     | Global state                                            |
!      | ~num~                     | ~int64_t~                                 | in     | Index of single electron                                |
!      | ~walk_num~                | ~int64_t~                                 | in     | Number of walkers                                       |
!      | ~elec_num~                | ~int64_t~                                 | in     | Number of electrons                                     |
!      | ~nucl_num~                | ~int64_t~                                 | in     | Number of nuclei                                        |
!      | ~type_nucl_num~           | ~int64_t~                                 | in     | Number of unique nuclei                                 |
!      | ~type_nucl_vector~        | ~int64_t[nucl_num]~                       | in     | IDs of unique nuclei                                    |
!      | ~aord_num~                | ~int64_t~                                 | in     | Number of coefficients                                  |
!      | ~a_vector~                | ~double[type_nucl_num][aord_num+1]~       | in     | List of coefficients                                    |
!      | ~en_distance_rescaled~    | ~double[walk_num][nucl_num][elec_num]~    | in     | Electron-nucleus rescaled distances                     |
!      | ~en_distance_rescaled_gl~ | ~double[walk_num][nucl_num][elec_num][4]~ | in     | Electron-nucleus rescaled distance derivatives          |
!      | ~en_rescaled_single~      | ~double[walk_num][nucl_num]~              | in     | Electron-nucleus rescaled single distances              |
!      | ~en_rescaled_single_gl~   | ~double[walk_num][nucl_num][4]~           | in     | Electron-nucleus rescaled single distance derivatives   |
!      | ~delta_en_gl~             | ~double[walk_num][elec_num][4]~           | out    | Single electron-nucleus Jastrow gradients and Laplacian |
!      |---------------------------+-------------------------------------------+--------+---------------------------------------------------------|


function qmckl_compute_jastrow_champ_single_en_gl_doc( &
     context, num_in, walk_num, elec_num, nucl_num, type_nucl_num, &
     type_nucl_vector, aord_num, a_vector, &
     en_distance_rescaled, en_distance_rescaled_gl, en_rescaled_single, en_rescaled_single_gl, delta_en_gl) &
     bind(C) result(info)
  use qmckl
  implicit none

  integer (qmckl_context), intent(in), value :: context
  integer (c_int64_t) , intent(in)  , value :: num_in
  integer (c_int64_t) , intent(in)  , value :: walk_num
  integer (c_int64_t) , intent(in)  , value :: elec_num
  integer (c_int64_t) , intent(in)  , value :: nucl_num
  integer (c_int64_t) , intent(in)  , value :: type_nucl_num
  integer (c_int64_t) , intent(in)          :: type_nucl_vector(nucl_num)
  integer (c_int64_t) , intent(in)  , value :: aord_num
  real    (c_double ) , intent(in)          :: a_vector(aord_num+1,type_nucl_num)
  real    (c_double ) , intent(in)          :: en_distance_rescaled(elec_num,nucl_num,walk_num)
  real    (c_double ) , intent(in)          :: en_distance_rescaled_gl(4, elec_num,nucl_num,walk_num)
  real    (c_double ) , intent(in)          :: en_rescaled_single(nucl_num,walk_num)
  real    (c_double ) , intent(in)          :: en_rescaled_single_gl(4, nucl_num,walk_num)
  real    (c_double ) , intent(out)         :: delta_en_gl(4,elec_num,walk_num)
  integer(qmckl_exit_code)                   :: info

  integer*8 :: i, a, k, nw, ii, num
  double precision   :: x, x1, kf, x_old, x1_old
  double precision   :: denom, invdenom, invdenom2, f
  double precision   :: denom_old, invdenom_old, invdenom2_old, f_old
  double precision   :: grad_c2, grad_c2_old
  double precision   :: dx(4), dx_old(4)

  num = num_in + 1

  info = QMCKL_SUCCESS

  if (context == QMCKL_NULL_CONTEXT) then
     info = QMCKL_INVALID_CONTEXT
     return
  endif

  if (walk_num <= 0) then
     info = QMCKL_INVALID_ARG_3
     return
  endif

  if (elec_num <= 0) then
     info = QMCKL_INVALID_ARG_4
     return
  endif

  if (nucl_num <= 0) then
     info = QMCKL_INVALID_ARG_5
     return
  endif

  if (aord_num < 0) then
     info = QMCKL_INVALID_ARG_8
     return
  endif


  do nw =1, walk_num
    delta_en_gl(:,:,nw) = 0.0d0

    do a = 1, nucl_num

          x_old = en_distance_rescaled(num,a,nw)
          x = en_rescaled_single(a,nw)


          denom = 1.0d0 + a_vector(2, type_nucl_vector(a)+1) * x
          invdenom = 1.0d0 / denom
          invdenom2 = invdenom*invdenom

          denom_old = 1.0d0 + a_vector(2, type_nucl_vector(a)+1) * x_old
          invdenom_old = 1.0d0 / denom_old
          invdenom2_old = invdenom_old*invdenom_old

          dx(1) = en_rescaled_single_gl(1,a,nw)
          dx(2) = en_rescaled_single_gl(2,a,nw)
          dx(3) = en_rescaled_single_gl(3,a,nw)
          dx(4) = en_rescaled_single_gl(4,a,nw)

          dx_old(1) = en_distance_rescaled_gl(1,num,a,nw)
          dx_old(2) = en_distance_rescaled_gl(2,num,a,nw)
          dx_old(3) = en_distance_rescaled_gl(3,num,a,nw)
          dx_old(4) = en_distance_rescaled_gl(4,num,a,nw)

          f = a_vector(1, type_nucl_vector(a)+1) * invdenom2
          grad_c2 = dx(1)*dx(1) + dx(2)*dx(2) + dx(3)*dx(3)

          f_old = a_vector(1, type_nucl_vector(a)+1) * invdenom2_old
          grad_c2_old = dx_old(1)*dx_old(1) + dx_old(2)*dx_old(2) + dx_old(3)*dx_old(3)

          delta_en_gl(1,num,nw) = delta_en_gl(1,num,nw) + f * dx(1) - f_old * dx_old(1)
          delta_en_gl(2,num,nw) = delta_en_gl(2,num,nw) + f * dx(2) - f_old * dx_old(2)
          delta_en_gl(3,num,nw) = delta_en_gl(3,num,nw) + f * dx(3) - f_old * dx_old(3)
          delta_en_gl(4,num,nw) = delta_en_gl(4,num,nw) &
               + f * (dx(4) - 2.d0 * a_vector(2, type_nucl_vector(a)+1) * grad_c2 * invdenom) &
               - f_old * (dx_old(4) - 2.d0 * a_vector(2, type_nucl_vector(a)+1) * grad_c2_old * invdenom_old)


           kf = 2.d0
           x1 = x
           x = 1.d0
           x1_old = x_old
           x_old = 1.d0
           do k=2, aord_num
              f = a_vector(k+1,type_nucl_vector(a)+1) * kf * x
              f_old = a_vector(k+1,type_nucl_vector(a)+1) * kf * x_old
              delta_en_gl(1,num,nw) = delta_en_gl(1,num,nw) + f * x1 * dx(1) - f_old * x1_old * dx_old(1)
              delta_en_gl(2,num,nw) = delta_en_gl(2,num,nw) + f * x1 * dx(2) - f_old * x1_old * dx_old(2)
              delta_en_gl(3,num,nw) = delta_en_gl(3,num,nw) + f * x1 * dx(3) - f_old * x1_old * dx_old(3)
              delta_en_gl(4,num,nw) = delta_en_gl(4,num,nw) &
                   + f * (x1 * dx(4) + (kf-1.d0) * grad_c2) &
                   - f_old * (x1_old * dx_old(4) + (kf-1.d0) * grad_c2_old)
              x = x*x1
              x_old = x_old*x1_old
              kf = kf + 1.d0
           end do
     end do
  end do

end function qmckl_compute_jastrow_champ_single_en_gl_doc
