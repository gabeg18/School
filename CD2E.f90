MODULE CD2E

  USE, INTRINSIC :: ISO_FORTRAN_ENV     ! ISO Fortran Environment
  IMPLICIT NONE
  INTEGER, PARAMETER :: rDef = REAL64   ! Defining precision

  ! A module to calculate the time rate of change of the heat
  ! equation using 2nd Order Central Spatial Differencing and
  ! 1st Order Explicit Time Marching schemes

  INTERFACE GetCD2E1Update
    MODULE PROCEDURE CalculateCD2E1Update
  END INTERFACE

  CONTAINS

  SUBROUTINE CalculateCD2E1Update(fVec,     &
                                  mu,       &
                                  dX,       &
                                  dT,       &
                                  fNewMin,  &
                                  fNewMax,  & 
                                  fNewVec)

    ! Dummy Variable Declaration
    REAL(KIND = rDef), DIMENSION(:), INTENT(IN) :: fVec     ! Provided function

    REAL(KIND = rDef), INTENT(IN) :: dX,      &             ! Delta x
                                     mu,      &             ! (Diffusion?) Coefficient
                                     dT,      &             ! Delta t
                                     fNewMin, &             ! Lower Boundary Condition
                                     fNewMax                ! Upper Boundary Condition

    REAL(KIND = rDef), DIMENSION(:), INTENT(OUT) :: fNewVec ! Output Function 

    ! Local Variable Declaration
    REAL(KIND = rDef) :: sigma
    INTEGER :: i, iMax

    iMax  = UBOUND(fNewVec,1)
    sigma = mu*dT/(dX*dX)

    ! Boundary conditions are imposed at the upper and lower boundaries
    ! Therefore, we will compute the solution from the values:
    ! i = 2 , iMax - 1

    ! Calculate Resulting Function
    DO i = 2, iMax-1
      fNewVec(i) = fVec(i) + sigma*(fVec(i+1)  &
                           - 2.0_rDef*(fVec(i)) &
                           + fVec(i-1))
    END DO

    RETURN

  END SUBROUTINE CalculateCD2E1Update

END MODULE CD2E
