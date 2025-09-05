subroutine convert(upds, s_inv, Updates, Inverse, nupdates, lds, dim)
  implicit none
  integer*8 , intent(in)                             :: lds, dim, nupdates
  real*8    , intent(in)                             :: upds(nupdates * lds)
  real*8    , intent(in)                             :: s_inv(dim * lds)  
  real*8    , intent(out) , dimension(dim, nupdates) :: Updates
  real*8    , intent(out) , dimension(dim, dim)      :: Inverse

  integer*8                                          :: i, j

  ! Construct Updates: lds x nupdates
  do i = 1, nupdates
    do j = 1, dim
      Updates(j, i) = upds((i - 1) * lds + j)
    end do
  end do

  ! Construct Inverse: dim x lds
  do i = 1, dim
    do j = 1, dim
      Inverse(i, j) = s_inv((i - 1) * lds + j)
    end do
  end do
end subroutine convert

subroutine copy_back_inv(Inverse, s_inv, lds, dim)
  implicit none
  integer*8 , intent(in)                       :: lds, dim
  real*8    , intent(in) , dimension(dim, dim) :: Inverse
  real*8    , intent(out)                      :: s_inv(dim * lds)  

  integer*8                                    :: i, j

  ! Copy updated inverse back to s_inv
  do i = 1, dim
    do j = 1, dim
      s_inv((i - 1) * lds + j) = Inverse(i, j)
    end do
  end do
end subroutine copy_back_inv

subroutine copy_back_lu(Later_updates, later_upds, lds, dim, nupdates)
  implicit none
  integer*8 , intent(in)                            :: lds, dim, nupdates
  real*8    , intent(in) , dimension(dim, nupdates) :: Later_updates
  real*8    , intent(out)                           :: later_upds(nupdates * lds)  

  integer*8                                    :: i, j

  ! Copy updated inverse back to s_inv
  do i = 1, nupdates
    do j = 1, dim
      later_upds((i - 1) * lds + j) = Later_updates(j, i)
    end do
  end do
end subroutine copy_back_lu

integer function qmckl_sm_naive_doc_f(context, &
    lds, dim, &
    nupdates, &
    upds, &
    updates_index, &
    breakdown, &
    s_inv, &
    determinant) result(info)

  use qmckl
  implicit none
  integer*8 , intent(in)    :: context
  integer*8 , intent(in)    :: lds, dim
  integer*8 , intent(in)    :: nupdates
  integer*8 , intent(in)    :: updates_index(nupdates)
  real*8    , intent(in)    :: upds(nupdates * lds)
  real*8    , intent(in)    :: breakdown
  real*8    , intent(inout) :: s_inv(dim * lds)
  real*8    , intent(inout) :: determinant

  real*8 , dimension(dim, nupdates) :: Updates
  real*8 , dimension(dim, dim)      :: Inverse
  real*8 , dimension(dim)           :: C
  real*8 , dimension(dim)           :: D
  real*8                            :: denominator, idenominator, update
  integer*8                         :: i, j, l, row

  info = QMCKL_FAILURE

  if (context == QMCKL_NULL_CONTEXT) then
    info = QMCKL_INVALID_CONTEXT
    return
  endif

  ! Convert 'upds' and 's_inv' into the more easily readable Fortran
  ! matrices 'Updates' and 'Inverse'.
  call convert(upds, s_inv, Updates, Inverse, nupdates, lds, dim)

  l = 1;
  ! For each update do...
  do while (l < nupdates + 1)

    ! Compute C = S^{-1}U(l)
    do i = 1, dim
      C(i) = 0
      do j = 1, dim
        C(i) = C(i) + Inverse(i, j) * Updates(j, l)
      end do
    end do

    ! Compute denominator = 1 + V(l)^TC
    row = updates_index(l)
    denominator = 1 + C(row)

    ! Return early if denominator is too small
    if (abs(denominator) < breakdown) return
    idenominator = 1 / denominator

    ! Update det(S)
    determinant = determinant * denominator

    ! selecting column: v_l^T * S_inv
    D = Inverse(row, :)

    ! A^{-1} = A^{-1} - C x D / denominator
    do i = 1, dim
      do j = 1, dim
        update = C(i) * D(j) * idenominator
        Inverse(i, j) = Inverse(i, j) - update
      end do
    end do

    l = l + 1
  end do

  ! Copy updated inverse back to s_inv
  call copy_back_inv(Inverse, s_inv, lds, dim)

  info = QMCKL_SUCCESS

end function qmckl_sm_naive_doc_f

! C interface (not directly exposed)
! The following Fortran function ~qmckl_sm_naive_doc~ makes sure
! that the pedagogical kernel ~qmckl_sm_naive_doc_f~, written in
! Fortran, can be called from C using the ~ISO_C_BINDING~. The Fortran function ~qmckl_sm_naive_doc~ will be exposed in the header file 'qmckl.h'
! for C users and in the module file 'qmckl_f.F90' for Fortran users.

! #+CALL: generate_c_interface(table=qmckl_sm_naive_args,rettyp=get_value("CRetType"),fname="qmckl_sm_naive_doc")

! #+RESULTS:

integer(c_int32_t) function qmckl_sm_naive_doc &
    (context, LDS, Dim, N_updates, Updates, Updates_index, breakdown, Slater_inv, determinant) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: LDS
  integer (c_int64_t) , intent(in)  , value :: Dim
  integer (c_int64_t) , intent(in)  , value :: N_updates
  real    (c_double ) , intent(in)          :: Updates(N_updates*LDS)
  integer (c_int64_t) , intent(in)          :: Updates_index(N_updates)
  real    (c_double ) , intent(in)  , value :: breakdown
  real    (c_double ) , intent(inout)        :: Slater_inv(Dim*LDS)
  real    (c_double ) , intent(inout)        :: determinant

  integer(c_int32_t), external :: qmckl_sm_naive_doc_f
  info = qmckl_sm_naive_doc_f &
         (context, LDS, Dim, N_updates, Updates, Updates_index, breakdown, Slater_inv, determinant)

end function qmckl_sm_naive_doc

integer function qmckl_sm_splitting_core_doc_f( &
    context, &
    lds, dim, &
    nupdates, &
    upds, &
    updates_index, &
    breakdown, &
    s_inv, &
    later_upds, &
    Later_index, &
    Later, &
    determinant) result(info)

  use qmckl
  implicit none
  integer*8 , intent(in)            :: context
  integer*8 , intent(in)            :: lds, dim
  integer*8 , intent(in)            :: nupdates
  integer*8 , intent(in)            :: updates_index(nupdates)
  real*8    , intent(in)            :: upds(lds * nupdates)
  real*8    , intent(in)            :: breakdown
  real*8    , intent(inout)         :: s_inv(dim * lds)
  real*8    , intent(inout)         :: determinant
  integer*8 , intent(inout)         :: Later
  integer*8 , intent(inout)         :: Later_index(nupdates)
  real*8    , intent(inout)         :: later_upds(lds * nupdates)

  real*8 , dimension(dim, nupdates) :: Updates
  real*8 , dimension(dim, nupdates) :: Later_updates
  real*8 , dimension(dim, dim)      :: Inverse
  real*8 , dimension(dim)           :: C
  real*8 , dimension(dim)           :: D
  real*8                            :: denominator, idenominator, update
  integer*8                         :: i, j, l, row

  info = QMCKL_FAILURE

  if (context == QMCKL_NULL_CONTEXT) then
    info = QMCKL_INVALID_CONTEXT
    return
  endif

  ! Convert 'upds' and 's_inv' into the more easily readable Fortran
  ! matrices 'Updates' and 'Inverse'.
  call convert(upds, s_inv, Updates, Inverse, nupdates, lds, dim)

  l = 1;
  ! For each update do...
  do while (l < nupdates + 1)

    ! Compute C = S^{-1}U(l)
    do i = 1, dim
      C(i) = 0
      do j = 1, dim
        C(i) = C(i) + Inverse(i, j) * Updates(j, l)
      end do
    end do

    ! Compute denominator = 1 + V(l)^TC
    row = updates_index(l)
    denominator = 1 + C(row)

    ! If denominator is too close to zero:
    ! - Split update in 2 before storing in Later_updates
    ! - Split previously computed vector C in 2
    ! - Recompute the denominator
    if (abs(denominator) < breakdown) then
      do i = 1, dim
        Later_updates(i, l) = Updates(i, l) / 2
        C(i) = C(i) / 2
      end do
      Later_index(Later + 1) = updates_index(l)
      Later = Later + 1
      denominator = 1 + C(row)
    end if

    idenominator = 1 / denominator

    ! Update det(S)
    determinant = determinant * denominator

    ! selecting column: v_l^T * S_inv
    D = Inverse(row, :)

    ! A^{-1} = A^{-1} - C x D / denominator
    do i = 1, dim
      do j = 1, dim
        update = C(i) * D(j) * idenominator
        Inverse(i, j) = Inverse(i, j) - update
      end do
    end do

    l = l + 1
  end do

  ! Copy updated inverse and later updates
  ! back to s_inv and later_upds
  call copy_back_inv(Inverse, s_inv, lds, dim)
  call copy_back_lu(Later_Updates, later_upds, lds, dim, nupdates)

  info = QMCKL_SUCCESS

end function qmckl_sm_splitting_core_doc_f

! C interface to the pedagogical kernel (not directly exposed)
! The function ~qmckl_sm_splitting_core_doc~ makes sure that
! ~qmckl_sm_splitting_core_doc_f~ can be called from C using the
! ~ISO_C_BINDING~. Function ~qmckl_sm_splitting_core_doc~ will be
! exposed in ~qmckl.h~ and ~qmckl_f.F90~, but
! ~qmckl_sm_splitting_core_doc_f~ will not.

! #+CALL: generate_c_interface(table=qmckl_sm_splitting_core_args,rettyp=get_value("CRetType"),fname="qmckl_sm_splitting_core_doc")

! #+RESULTS:

integer(c_int32_t) function qmckl_sm_splitting_core_doc &
    (context, &
     LDS, &
     Dim, &
     N_updates, &
     Updates, &
     Updates_index, &
     breakdown, &
     Slater_inv, &
     later_updates, &
     later_index, &
     later, &
     determinant) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: LDS
  integer (c_int64_t) , intent(in)  , value :: Dim
  integer (c_int64_t) , intent(in)  , value :: N_updates
  real    (c_double ) , intent(in)          :: Updates(LDS*N_updates)
  integer (c_int64_t) , intent(in)          :: Updates_index(N_updates)
  real    (c_double ) , intent(in)  , value :: breakdown
  real    (c_double ) , intent(inout)        :: Slater_inv(Dim*LDS)
  real    (c_double ) , intent(inout)        :: later_updates(LDS*N_updates)
  integer (c_int64_t) , intent(inout)        :: later_index(N_updates)
  integer (c_int64_t) , intent(inout)        :: later
  real    (c_double ) , intent(inout)        :: determinant

  integer(c_int32_t), external :: qmckl_sm_splitting_core_doc_f
  info = qmckl_sm_splitting_core_doc_f &
         (context, &
     LDS, &
     Dim, &
     N_updates, &
     Updates, &
     Updates_index, &
     breakdown, &
     Slater_inv, &
     later_updates, &
     later_index, &
     later, &
     determinant)

end function qmckl_sm_splitting_core_doc

integer function qmckl_woodbury_2x2_doc_f(&
    context, &
    lds, dim, &
    upds, &
    updates_index, &
    breakdown, &
    s_inv, &
    determinant) result(info)

  use qmckl
  implicit none
  integer*8 , intent(in)            :: context
  integer*8 , intent(in)            :: lds, dim
  integer*8 , intent(in)            :: updates_index(2)
  real*8    , intent(in)            :: upds(2 * lds)
  real*8    , intent(in)            :: breakdown
  real*8    , intent(inout)         :: s_inv(dim * lds)
  real*8    , intent(inout)         :: determinant

  integer*8 , dimension(2, dim)     :: V
  integer*8 , dimension(2, 2)       :: Id
  real*8    , dimension(dim, dim)   :: Inverse
  real*8    , dimension(dim, 2)     :: Updates, C
  real*8    , dimension(2, 2)       :: D, invD
  real*8    , dimension(2, dim)     :: E, F
  
  real*8                            :: detD, idenominator, update
  integer*8                         :: i, j, k, l

  info = QMCKL_FAILURE

  if (context == QMCKL_NULL_CONTEXT) then
    info = QMCKL_INVALID_CONTEXT
    return
  endif

  ! Construct V(2, dim) matrix
  V = 0
  V(1, updates_index(1)) = 1
  V(2, updates_index(2)) = 1

  ! Construct Id(2, 2) matrix
  Id = 0
  Id(1, 1) = 1
  Id(2, 2) = 1

  ! Convert 'upds' and 's_inv' into the more easily readable Fortran
  ! matrices 'Updates' and 'Inverse'.
  call convert(upds, s_inv, Updates, Inverse, int(2,8), lds, dim)

  ! Compute C(dim, 2) = Inverse(dim, dim) x Updates(dim, 2)
  C = 0
  do i = 1, dim
    do j = 1, 2
      do k = 1, dim
        C(i, j) = C(i, j) + Inverse(i, k) * Updates(k, j)
      end do
    end do
  end do

  ! Construct matrix D(2, 2) := I(2, 2) + V(2, dim) x C(dim, 2)
  D = 0
  do i = 1, 2 
     do j = 1, 2
        do k = 2, dim
           D(i, j) = D(i, j) + V(i, k) * C(k, j)
        end do
      end do
  end do
  D = Id + D
  
  ! Compute determinant := det(D) explicitly
  detD = D(1,1) * D(2,2) - D(1,2) * D(2,1)

    ! Return early if det(D) is too small
  if (abs(detD) < breakdown) return

  ! Update det(S)
  determinant = determinant * detD

  ! Compute inv(D) explicitly
  invD(1,1) =   D(2,2)
  invD(1,2) = - D(1,2)
  invD(2,1) = - D(2,1)
  invD(2,2) =   D(1,1)
  invD = invD / detD

  ! Compute E(2, dim) := V(2, dim) x Inverse(dim, dim)
  E = 0
  do i = 1, 2 
    do j = 1, dim
      do k = 1, dim
        E(i, j) = E(i, j) + V(i, k) * Inverse(k, j)
      end do
    end do
  end do

  ! Compute F(2, dim) := invD(2, 2) x E(2, dim)
  F = 0
  do i = 1, 2
    do j = 1, dim
      do k = 1, 2
        F(i, j) = F(i, j) + invD(i, k) * E(k, j)
      end do
    end do
  end do

  ! Compute Inverse(dim, dim) := Inverse(dim, dim) - C(dim, 2) x F(2, dim)
  do i = 1, dim
    do j = 1, dim
      do k = 1, 2
        Inverse(i, j) = Inverse(i, j) - C(i, k) * F(k, j)
      end do
    end do
  end do

  ! Copy updated inverse and later updates
  ! back to s_inv and later_upds
  call copy_back_inv(Inverse, s_inv, lds, dim)

  info = QMCKL_SUCCESS

end function qmckl_woodbury_2x2_doc_f

! C interface (not directly exposed)
! The function ~qmckl_sm_splitting_core_doc~ makes sure that
! ~qmckl_sm_splitting_core_doc_f~ can be called from C using the
! ~ISO_C_BINDING~. Function ~qmckl_sm_splitting_core_doc~ will be
! exposed in ~qmckl.h~ and ~qmckl_f.F90~, but
! ~qmckl_sm_splitting_core_doc_f~ will not.

! #+CALL: generate_c_interface(table=qmckl_woodbury_2x2_args,rettyp=get_value("CRetType"),fname="qmckl_woodbury_2x2_doc")

! #+RESULTS:

integer(c_int32_t) function qmckl_woodbury_2x2_doc &
    (context, LDS, Dim, Updates, Updates_index, breakdown, Slater_inv, determinant) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: LDS
  integer (c_int64_t) , intent(in)  , value :: Dim
  real    (c_double ) , intent(in)          :: Updates(2*Dim)
  integer (c_int64_t) , intent(in)          :: Updates_index(2)
  real    (c_double ) , intent(in)  , value :: breakdown
  real    (c_double ) , intent(inout)        :: Slater_inv(LDS*Dim)
  real    (c_double ) , intent(inout)        :: determinant

  integer(c_int32_t), external :: qmckl_woodbury_2x2_doc_f
  info = qmckl_woodbury_2x2_doc_f &
         (context, LDS, Dim, Updates, Updates_index, breakdown, Slater_inv, determinant)

end function qmckl_woodbury_2x2_doc

integer function qmckl_woodbury_3x3_doc_f(&
    context, &
    lds, dim, &
    upds, &
    updates_index, &
    breakdown, &
    s_inv, &
    determinant) result(info)

  use qmckl
  implicit none
  integer*8 , intent(in)            :: context
  integer*8 , intent(in)            :: lds, dim
  integer*8 , intent(in)            :: updates_index(3)
  real*8    , intent(in)            :: upds(3 * lds)
  real*8    , intent(in)            :: breakdown
  real*8    , intent(inout)         :: s_inv(dim * lds)
  real*8    , intent(inout)         :: determinant

  integer*8 , dimension(3, dim)     :: V
  integer*8 , dimension(3, 3)       :: Id
  real*8    , dimension(dim, dim)   :: Inverse
  real*8    , dimension(dim, 3)     :: Updates, C
  real*8    , dimension(3, 3)       :: D, invD
  real*8    , dimension(3, dim)     :: E, F
  
  real*8                            :: detD, idetD, idenominator, update
  integer*8                         :: i, j, k, l

  info = QMCKL_FAILURE

  if (context == QMCKL_NULL_CONTEXT) then
    info = QMCKL_INVALID_CONTEXT
    return
  endif

  ! Construct V(3, dim) matrix
  V = 0
  V(1, updates_index(1)) = 1
  V(2, updates_index(2)) = 1
  V(3, updates_index(3)) = 1

  ! Construct Id(3, 3) matrix
  Id = 0
  Id(1, 1) = 1
  Id(2, 2) = 1
  Id(3, 3) = 1

  ! Convert 'upds' and 's_inv' into the more easily readable Fortran
  ! matrices 'Updates' and 'Inverse'.
  call convert(upds, s_inv, Updates, Inverse, int(3,8), lds, dim)

  ! Compute C(dim, 3) = Inverse(dim, dim) x Updates(dim, 3)
  C = 0
  do i = 1, dim
    do j = 1, 3
      do k = 1, dim
        C(i, j) = C(i, j) + Inverse(i, k) * Updates(k, j)
      end do
    end do
  end do

  ! Construct matrix D(3, 3) := I(3, 3) + V(3, dim) x C(dim, 3)
  D = 0
  do i = 1, 3 
     do j = 1, 3
        do k = 3, dim
           D(i, j) = D(i, j) + V(i, k) * C(k, j)
        end do
      end do
  end do
  D = Id + D
  
  ! Compute determinant := det(D) explicitly
  detD = D(1,1) * (D(2,2) * D(3,3) - D(2,3) * D(3,2)) - &
         D(1,2) * (D(2,1) * D(3,3) - D(2,3) * D(3,1)) + &
         D(1,3) * (D(2,1) * D(3,2) - D(2,2) * D(3,1))
         
  ! Return early if det(D) is too small
  if (abs(detD) < breakdown) return

  ! Update det(S)
  determinant = determinant * detD

  idetD = 1.0d0 / detD
  ! Compute inv(D) explicitly
  invD(1,1) =  (D(2,2) * D(3,3) - D(3,2) * D(2,3)) * idetD
  invD(1,2) = -(D(1,2) * D(3,3) - D(3,2) * D(1,3)) * idetD
  invD(1,3) =  (D(1,2) * D(2,3) - D(2,2) * D(1,3)) * idetD
  invD(2,1) = -(D(2,1) * D(3,3) - D(3,1) * D(2,3)) * idetD
  invD(2,2) =  (D(1,1) * D(3,3) - D(3,1) * D(1,3)) * idetD
  invD(2,3) = -(D(1,1) * D(2,3) - D(2,1) * D(1,3)) * idetD
  invD(3,1) =  (D(2,1) * D(3,2) - D(3,1) * D(2,2)) * idetD
  invD(3,2) = -(D(1,1) * D(3,2) - D(3,1) * D(1,2)) * idetD
  invD(3,3) =  (D(1,1) * D(2,2) - D(2,1) * D(1,2)) * idetD
  
  ! Compute E(3, dim) := V(3, dim) x Inverse(dim, dim)
  E = 0
  do i = 1, 3 
    do j = 1, dim
      do k = 1, dim
        E(i, j) = E(i, j) + V(i, k) * Inverse(k, j)
      end do
    end do
  end do

  ! Compute F(3, dim) := invD(3, 3) x E(3, dim)
  F = 0
  do i = 1, 3
    do j = 1, dim
      do k = 1, 3
        F(i, j) = F(i, j) + invD(i, k) * E(k, j)
      end do
    end do
  end do

  ! Compute Inverse(dim, dim) := Inverse(dim, dim) - C(dim, 3) x F(3, dim)
  do i = 1, dim
    do j = 1, dim
      do k = 1, 3
        Inverse(i, j) = Inverse(i, j) - C(i, k) * F(k, j)
      end do
    end do
  end do

  ! Copy updated inverse and later updates
  ! back to s_inv and later_upds
  call copy_back_inv(Inverse, s_inv, lds, dim)

  info = QMCKL_SUCCESS

end function qmckl_woodbury_3x3_doc_f

! C interface (not directly exposed)
! The function ~qmckl_sm_splitting_core_doc~ makes sure that
! ~qmckl_sm_splitting_core_doc_f~ can be called from C using the
! ~ISO_C_BINDING~. Function ~qmckl_sm_splitting_core_doc~ will be
! exposed in ~qmckl.h~ and ~qmckl_f.F90~, but
! ~qmckl_sm_splitting_core_doc_f~ will not.

! #+CALL: generate_c_interface(table=qmckl_woodbury_3x3_args,rettyp=get_value("CRetType"),fname="qmckl_woodbury_3x3_doc")

! #+RESULTS:

integer(c_int32_t) function qmckl_woodbury_3x3_doc &
    (context, LDS, Dim, Updates, Updates_index, breakdown, Slater_inv, determinant) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: LDS
  integer (c_int64_t) , intent(in)  , value :: Dim
  real    (c_double ) , intent(in)          :: Updates(3*Dim)
  integer (c_int64_t) , intent(in)          :: Updates_index(3)
  real    (c_double ) , intent(in)  , value :: breakdown
  real    (c_double ) , intent(inout)        :: Slater_inv(LDS*Dim)
  real    (c_double ) , intent(inout)        :: determinant

  integer(c_int32_t), external :: qmckl_woodbury_3x3_doc_f
  info = qmckl_woodbury_3x3_doc_f &
         (context, LDS, Dim, Updates, Updates_index, breakdown, Slater_inv, determinant)

end function qmckl_woodbury_3x3_doc

integer recursive function qmckl_sm_splitting_doc_f( &
    context, &
    lds, dim, &
    nupdates, &
    upds, &
    updates_index, &
    breakdown, &
    s_inv, &
    determinant) result(info)

  use qmckl
  implicit none
  integer*8 , intent(in)                :: context
  integer*8 , intent(in)                :: lds, dim
  integer*8 , intent(in)                :: nupdates
  integer*8 , intent(in)                :: updates_index(nupdates)
  real*8    , intent(in)                :: upds(lds * nupdates)
  real*8    , intent(in)                :: breakdown
  real*8    , intent(inout)             :: s_inv(dim * lds)
  real*8    , intent(inout)             :: determinant

  integer   , external                  :: qmckl_sm_splitting_core_doc_f
  
  integer*8                             :: Later
  integer*8 , dimension(nupdates)       :: Later_index
  real*8    , dimension(lds * nupdates) :: Later_updates

  info = QMCKL_FAILURE

  if (context == QMCKL_NULL_CONTEXT) then
    info = QMCKL_INVALID_CONTEXT
    return
  endif

  Later = 0
  Later_index = 0
  Later_updates = 0

  info = qmckl_sm_splitting_core_doc_f( &
    context, &
    lds, dim, &
    nupdates, &
    upds, &
    updates_index, &
    breakdown, &
    s_inv, &
    Later_updates, &
    Later_index, &
    Later, &
    determinant)
  
  if (Later > 0) then
    info = qmckl_sm_splitting_doc_f( &
      context, &
      lds, dim, &
      Later, &
      Later_updates, &
      Later_index, &
      breakdown, &
      s_inv, &
      determinant)
  end if

  info = QMCKL_SUCCESS

end function qmckl_sm_splitting_doc_f

! C interface to the pedagogical kernel (not directly exposed)
! The following Fortran function ~qmckl_sm_splitting_core_doc~ makes sure
! that the pedagogical kernel ~qmckl_sm_splitting_core_doc_f~, written in
! Fortran, can be called from C using the ~ISO_C_BINDING~. The Fortran function
! ~qmckl_sm_splitting_core_doc~ will be exposed in the header file 'qmckl.h'
! for C users and in the module file 'qmckl_f.F90' for Fortran users.

! #+CALL: generate_c_interface(table=qmckl_sm_splitting_args,rettyp=get_value("CRetType"),fname="qmckl_sm_splitting_doc")

! #+RESULTS:

integer(c_int32_t) function qmckl_sm_splitting_doc &
    (context, LDS, Dim, N_updates, Updates, Updates_index, breakdown, Slater_inv, determinant) &
    bind(C) result(info)

  use, intrinsic :: iso_c_binding
  implicit none

  integer (c_int64_t) , intent(in)  , value :: context
  integer (c_int64_t) , intent(in)  , value :: LDS
  integer (c_int64_t) , intent(in)  , value :: Dim
  integer (c_int64_t) , intent(in)  , value :: N_updates
  real    (c_double ) , intent(in)          :: Updates(N_updates*LDS)
  integer (c_int64_t) , intent(in)          :: Updates_index(N_updates)
  real    (c_double ) , intent(in)  , value :: breakdown
  real    (c_double ) , intent(inout)        :: Slater_inv(Dim*LDS)
  real    (c_double ) , intent(inout)        :: determinant

  integer(c_int32_t), external :: qmckl_sm_splitting_doc_f
  info = qmckl_sm_splitting_doc_f &
         (context, LDS, Dim, N_updates, Updates, Updates_index, breakdown, Slater_inv, determinant)

end function qmckl_sm_splitting_doc
