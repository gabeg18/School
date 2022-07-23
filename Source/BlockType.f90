MODULE block_type
  USE ISO_FORTRAN_ENV
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: grid_block
  INTEGER, PARAMETER :: rDef = REAL64

  ! This module contains a derived data type for the blocks
  ! to be used in a multi-block grid generator

  TYPE grid_block                                 

    INTEGER :: block_number,  &
               xMin,          &
               xMax,          &
               iMin,          &
               iMax

    REAL(KIND=rDef) :: sigma,    &
                       mu,       &
                       dX,       &
                       dT,       & 
                       uMin,     &
                       uMax,     &
                       uInitial, &
                       l2Err,    &
                       dU

    REAL(KIND=rDef), DIMENSION(:), ALLOCATABLE :: u,      &
                                                  x,      &
                                                  uNew,   &
                                                  uExact, &
                                                  fVec0,  & ! For the RK Scheme
                                                  fVecS,  & ! For the RK Scheme
                                                  rhs,    & ! For the I1 Scheme
                                                  dm1,    & ! For the I1 Scheme 
                                                  d0,     & ! For the I1 Scheme
                                                  dp1       ! For the I1 Scheme
                                                  

  END TYPE grid_block

END MODULE block_type


