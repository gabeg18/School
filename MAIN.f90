PROGRAM MAIN
  USE, INTRINSIC :: ISO_FORTRAN_ENV  ! Using Intrinsic Fortran Environment

  IMPLICIT NONE

  ! Variable Declaration
  INTEGER, PARAMETER :: rDef       = REAL64,  &
                        numSchemes = 3,       &
                        maxItn     = 1000000  ! In case it blows up

  REAL(KIND=rDef), PARAMETER :: xMin  =     0.0_rDef, &
                                xMax  =   100.0_rDef, &
                                l2Tol = 1.0e-10_rDef, &
                                l2Max =     1.0_rDef, & ! To stop the code
                                mu    =     0.1_rDef

  REAL(KIND=rDef), DIMENSION(3), PARAMETER :: sigmaMax = (/ &
                  0.5_rDef, &     ! Max Stability for E1
               1.0e15_rDef, &     ! Max Stability for I1
              0.6963234_rDef/)    ! Max Stability for RK4

  INTEGER :: j, i, &
             outFileUnit,  &
             outDataFileUnit,  &
             nG,nStep

  REAL(KIND=rDef) :: dX,        &   ! Delta x
                     dT,        &   ! Delta time
                     uMin,      &   ! Minimum value of function
                     uMax,      &   ! Maximum value of function
                     uInitial,  &   ! Initial value of function
                     l2Err,     &   ! L2 Error
                     dU             ! Change in function

  REAL(KIND=rDef), DIMENSION(:) :: u,       &   ! Function
                                   uNew,    &   ! Updated Function value
                                   uExact,  &   ! Exact Function value
                                   x            ! X value

  INTEGER :: iMax

  CHARACTER(LEN=30) :: outFileName, &
                       outDataFileName

  CHARACTER(LEN=5), DIMENSION(3) :: dType

  iMax = 101        ! For 0 <= x <= 100

  dType(1) = 'CD2E1'
  dType(2) = 'CD2I1'
  dType(3) = 'CD2RK'

  uMin  = 1.0_rDef  ! Lower Boundary Condition
  uMax  = 3.0_rDef  ! Upper Boundary Condition
  uInit = 2.0_rDef  ! Initial Condition

  !DO nG = 1,numGrids

  ! Determine value of delta x
  dX = (xMax - xMin)/REAL(iMax-1,rDef)  


  ! Calculate Exact Solution
  DO i = 1,iMax
  x(i) = xMin + REAL(i-1,rDef)*dX                         ! Initialize Space Domain
  uExact(i) = uMin + (uMax-uMin)*(x(i)-xMin)/(xMax-xMin)  ! Calculate solution
  END DO


  ! Calculation of Each Scheme
  DO j = 1, numSchemes    ! ADD BACK IN LATER FOR DIFFERENT SCHEMES
    uNew(j) = uMin
    uNew(iMax) = uMax

    ! Initialize Data
    DO i = 1, iMax
      u(i) = uInitial
    END DO

    ! Create Output File Name
    WRITE(outFileName,'(a,a,a,i3.3,a)') 'RunData',  &
                                        dType(j),   &
                                        'iMax',iMax,'.dat'

    ! Open File
    OPEN(NEWUNIT = outFileUnit,       & ! Specifying file unit value
         FILE    = TRIM(outFileName), & ! Calling aforementioned name
         FORM    = 'FORMATTED')         ! Formatted

    ! Set dT
    dT = sigmaMax(j)*(dX*dX)/mu         ! Determined by Stability Limit (sigmaMax)

    nStep = 0
    99 CONTINUE

    nStep = nStep + 1

    IF(j == 1) THEN ! E1 Scheme
      CALL GetCD2E1Update(fVec    = u,    &
                          dX      = dX,   &
                          mu      = mu,   &
                          fNewMin = uMin, &
                          fNewMax = uMax, &
                          fnewVec = uNew)

    ELSE
      CONTINUE
    END IF

    ! L2 Calculations

    l2Err = 0.0_rDef
    DO i = 1, iMax
      u(i)  = uNew(i)
      dU    = u(i) - uExact(i)
      l2Err = l2Err + (dU*dU)
    END DO
    l2Err = SQRT(l2Err/REAL(iMax,rDef))

    ! Write Output

    WRITE(outFileUnit,*) nStep,l2Err
    WRITE(6,*) nStep,l2Err

    ! Continue?
    IF ((l2Err < l2Max) .AND. &
        (l2Err > l2tol) .AND. &
        (nStep < maxItn)) THEN ! Keep going
      GO TO 99
    ELSE  ! Done
      CONTINUE
    END IF

    CLOSE(outFileUnit)

    WRITE(outFileDataName,'(a,a,a,i3.3,a)') 'SolutionData', &
                                             dType(j),      &
                                             'iMax',iMax,'.dat'

    OPEN(NEWUNIT = outFileDataName,        &
         FILE    = TRIM(outDataFileName),  &
         FORM    = 'FORMATTED')

    DO i = 1, iMax
      WRITE(outDataFileUnit,*) x(i),u(i),uExact(i)
    END DO

    CLOSE(OutDataFileUnit)

  END DO  ! Scheme Loop

  STOP

END PROGRAM MAIN

