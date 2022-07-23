MODULE Input_Processing
  USE, INTRINSIC :: ISO_FORTRAN_ENV
  USE block_type
  IMPLICIT NONE
  PRIVATE
  PUBLIC Sort_Data
  INTEGER, PARAMETER :: rDef = REAL64


  ! This module is for processing the input file
  ! supplied for the multi-block grid generator

  ! NOTE: The input file format is as follows:

!|  Block Number  | xMin  | xMax  | iMin  | iMax  |
!|----------------|-------|-------|-------|-------|


  CONTAINS

    ! Subroutine to sort the Input File Data 
    SUBROUTINE Sort_Data(grid_blocks,num_blocks)

    ! Dummy Variable Declaration(s)
    INTEGER, INTENT(IN) :: num_blocks
    TYPE(grid_block), DIMENSION(:), INTENT(INOUT) :: grid_blocks

    ! Local Varaiable Declaration(s)
    INTEGER :: i,j                    ! DO-loop integers

    TYPE(grid_block), DIMENSION(SIZE(grid_blocks)) :: temp  ! For Sorting   

    !! Sorting Data
    !------------------------------------------------------------------

    ! Sorting Algorithm:
    
    ! 1.) Utilize intrinsic FINDLOC to determine location
    !     where block number of temporary matches the 
    !     incremental i integer
    
    ! 2.) Take the values of the "temporary" array at this 
    !     location
    
    ! 3.) Set value of block array equal to these values


    ! Transfer Data to "Temp" Array
    DO i = 1,num_blocks
      
      temp(i) = grid_blocks(i) 

    END DO                    
    
    WRITE (0,'(/"Sorted Data:")')
    DO i = 1,num_blocks
      grid_blocks(i) = temp(FINDLOC(temp % block_number,i,1)) 
      WRITE (0,'(4(I3,2X),I3)')    grid_blocks(i) % block_number, &
                                   grid_blocks(i) % xMin,         &
                                   grid_blocks(i) % xMax,         &
                                   grid_blocks(i) % iMin,         &
                                   grid_blocks(i) % iMax
    END DO

    END SUBROUTINE Sort_Data
END MODULE Input_Processing


MODULE SetParams
  USE, INTRINSIC :: ISO_FORTRAN_ENV
  USE block_type
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: GetParams
  
  INTEGER, PARAMETER :: rDef = REAL64

INTERFACE GetParams
  MODULE PROCEDURE Calculate_Parameters
END INTERFACE

  CONTAINS

    ! This subroutine calculates the dX, mu, and dT for each block when provided
    ! the block values. These parameters are then passed back for use in the 
    ! actual solver

    SUBROUTINE Calculate_Parameters(grid_blocks,num_blocks)

      ! Dummy Argument Declarations
      TYPE(grid_block), DIMENSION(:), INTENT(INOUT) :: grid_blocks
      INTEGER, INTENT(IN) :: num_blocks

      ! Local Variables
      INTEGER :: i

      !num_blocks = 1

      WRITE(0,'(/"Calculated Parameters: "/)')

      DO i = 1,num_blocks

        grid_blocks(i) % dX =     (grid_blocks(i) % xMax - grid_blocks(i) % xMin)/ &
                                  REAL(grid_blocks(i) % iMax - 1,rDef)

        grid_blocks(i) % mu = 1.0_rDef/10.0_rDef

        grid_blocks(i) % dT = 0.5_rDef*(grid_blocks(i) % dX * grid_blocks(i) % dX)/ &
                                            grid_blocks(i) % mu



        WRITE(0,'("Block: ",I3/,                 &
          "dX: ",F5.3,2X,                        &
          "mu: ",F5.3,2X/)') grid_blocks(i) % block_number,  &
                             grid_blocks(i) % dX,            &
                             grid_blocks(i) % mu            
      END DO

    END SUBROUTINE Calculate_Parameters

END MODULE SetParams

MODULE Exact_Solution
  USE, INTRINSIC :: ISO_FORTRAN_ENV
  USE block_type
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: GetExactSoln

  INTEGER, PARAMETER :: rDef = REAL64

INTERFACE GetExactSoln
  MODULE PROCEDURE GlobalExactSoln
END INTERFACE

  CONTAINS 
    SUBROUTINE GlobalExactSoln(grid_blocks,num_blocks)

      ! Dummy Variable Declaration
      TYPE(grid_block), DIMENSION(:), INTENT(INOUT) :: grid_blocks
      INTEGER, INTENT(IN) ::  num_blocks

      ! Local Variable Declaration
      INTEGER :: i,j,  &
                 nStep  ! For "running total"

      ! Initialize Variables
      grid_blocks(num_blocks) % uMin     = 1.0_rDef
      grid_blocks(num_blocks) % uMax     = 3.0_rDef
      grid_blocks(num_blocks) % uInitial = 2.0_rDef

      ! Calculate Global Exact Solution
      DO i = 1,grid_blocks(num_blocks) % iMax
        grid_blocks(num_blocks) % uExact(i) =  grid_blocks(num_blocks) % uMin   &
                                            + (grid_blocks(num_blocks) % uMax   &
                                            -  grid_blocks(num_blocks) % uMin)  &
                                            * (grid_blocks(num_blocks) % x(i)   & 
                                            -  grid_blocks(num_blocks) % xMin)  &
                                            / (grid_blocks(num_blocks) % xMax   &
                                            -  grid_blocks(num_blocks) % xMin)
      END DO

      ! Map Exact Solution to Each Block
      nStep = 0
      DO i = 1,num_blocks
        DO j = 1,grid_blocks(i) % iMax
          nStep = nStep + 1

          ! Only Map to "Actual" Blocks (Block 5 is "Master" Block)
          IF (i /= num_blocks) THEN

            ! Take uExact(nStep) since j resets each loop
            grid_blocks(i) % uExact(j) = grid_blocks(num_blocks) % uExact(nStep)

          END IF

        END DO
      END DO

    END SUBROUTINE GlobalExactSoln

END MODULE Exact_Solution
  



MODULE CD2E
  USE, INTRINSIC :: ISO_FORTRAN_ENV
  USE block_type
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: GetCD2E1Update, &
            GetCD2RKUpdate

  INTEGER, PARAMETER :: rDef = REAL64

 INTERFACE GetCD2E1Update 
   MODULE PROCEDURE CalculateCD2E1Update
 END INTERFACE

 INTERFACE GetCD2RKUpdate 
   MODULE PROCEDURE CalculateCD2RKUpdate
 END INTERFACE

  CONTAINS 
    SUBROUTINE CalculateCD2E1Update(grid_blocks,num_blocks)

    ! Dummy Variables
    TYPE(grid_block), DIMENSION(:), INTENT(INOUT) :: grid_blocks
    INTEGER, INTENT(IN) :: num_blocks

    ! Local Variables
    INTEGER :: i,j
    
    ! Exact Solution

    ! Set Min, Max, and Initial Conditions
    grid_blocks(1) % uMin          = 1.0_rDef   ! Global Minimum for FIRST block
    grid_blocks(num_blocks) % uMax = 3.0_rDef   ! Global Maximum for LAST block

    ! Initialize Data
    DO i = 1,num_blocks

      grid_blocks(i) % sigma = grid_blocks(i) % mu * grid_blocks(i) % dT/ &
                              (grid_blocks(i) % dX * grid_blocks(i) % dX)

    ! Flow Solution (Starting with 1 block)

      ! Apply B.C.
      grid_blocks(i) % uNew(grid_blocks(i) % iMin) = grid_blocks(i) % &
                                                     uExact(grid_blocks(i) % iMin)

      grid_blocks(i) % uNew(grid_blocks(i) % iMax) = grid_blocks(i) % &
                                                     uExact(grid_blocks(i) % iMax)

      ! Solving for Interior Domain
      DO j = 2,(grid_blocks(i) % iMax) - 1
        grid_blocks(i) % uNew(j) = grid_blocks(i) % u(j) + grid_blocks(i) % sigma *      &
                                  (grid_blocks(i) % u(j+1)                               &
                       -  2.0_rDef*grid_blocks(i) % u(j)                                 &
                       +           grid_blocks(i) % u(j-1))
      END DO
    END DO

    RETURN

    END SUBROUTINE CalculateCD2E1Update

    SUBROUTINE CalculateCD2RKUpdate(grid_blocks,num_blocks)

      ! Dummy Variable Declarations
      TYPE(grid_block), DIMENSION(:), INTENT(INOUT) :: grid_blocks
      INTEGER, INTENT(IN) :: num_blocks

      ! Local Variable Declarations
      REAL(KIND=rDef), DIMENSION(4), PARAMETER :: dTFrac = (/ 1.0_rDef/4.0_rDef,  & ! 1st Stage
                                                              1.0_rDef/3.0_rDef,  & ! 2nd Stage
                                                              1.0_rDef/2.0_rDef,  & ! 3rd Stage
                                                              1.0_rDef  /)          ! 4th Stage

      INTEGER :: i,j, &
                 nS   ! RK Stage 

      !
      ! Solving from i=2 to i=iMax-1 -- B.C. on i=1 and i=iMax
      !

      ! Effect of BC (and removing i-1 and i=iMax from the update, since we are specifying
      ! the values)

      DO i = 1,num_blocks
        
        grid_blocks(i) % sigma = grid_blocks(i) % mu * grid_blocks(i) % dT/ &
                                (grid_blocks(i) % dX * grid_blocks(i) % dX)

        ! Set Min, Max, and Initial Conditions
        grid_blocks(1) % uMin          = 1.0_rDef   ! Global Minimum for FIRST block
        grid_blocks(num_blocks) % uMax = 3.0_rDef   ! Global Maximum for LAST block

        ! Apply B.C.
        grid_blocks(i) % uNew(grid_blocks(i) % iMin) = grid_blocks(i) % &
                                                       uExact(grid_blocks(i) % iMin)

        grid_blocks(i) % uNew(grid_blocks(i) % iMax) = grid_blocks(i) % &
                                                       uExact(grid_blocks(i) % iMax)

        ! Initializing RK Arrays
        DO j = 1,grid_blocks(i) % iMax
          grid_blocks(i) % fVec0(j) = grid_blocks(i) % u(j)
          grid_blocks(i) % fVecS(j) = grid_blocks(i) % u(j)
        END DO

      END DO

      ! RK Scheme (Beginning with 1 Block)
      DO i = 1,num_blocks
        
        ! Loop Through Each Stage
        DO nS = 1,4

          ! Loop Through Interior Domain
          DO j = 2,grid_blocks(i) % iMax-1
            grid_blocks(i) % uNew(j) = grid_blocks(i) % fVec0(j) + dTFrac(nS) &
                                     * grid_blocks(i) % sigma                 &
                                    * (grid_blocks(i) % fVecS(j+1)            &
                          - 2.0_rDef * grid_blocks(i) % fVecS(j)              &
                                     + grid_blocks(i) % fVecS(j-1))
          END DO
          IF (nS < 4) THEN
            DO j = 1,grid_blocks(i) % iMax
              grid_blocks(i) % fVecS(j) = grid_blocks(i) % uNew(j)
            END DO
          ELSE
            CONTINUE
          END IF
        END DO
      END DO

    RETURN
    END SUBROUTINE CalculateCD2RKUpdate

END MODULE CD2E

MODULE CD2I
  USE, INTRINSIC :: ISO_FORTRAN_ENV
  USE block_type
  IMPLICIT NONE
  PRIVATE
  PUBLIC GetCD2I1Update

  INTEGER, PARAMETER :: rDef = REAL64

  INTERFACE GetCD2I1Update
    MODULE PROCEDURE CalculateCD2I1Update
  END INTERFACE

  CONTAINS
    SUBROUTINE CalculateCD2I1Update(grid_blocks,num_blocks)

      ! Dummy Variable Declaration
      TYPE(grid_block), DIMENSION(:), INTENT(INOUT) :: grid_blocks
      INTEGER, INTENT(IN) :: num_blocks

      ! Local Variable Declaration
      INTEGER :: i,j
      REAL(KIND=rDef) :: fac

      !
      ! Solving from i=2 to i=iMax-1 -- B.C. on i=1 and i=iMax
      !

      ! Loop Through Each Block
      DO i = 1,num_blocks
        grid_blocks(i) % sigma = grid_blocks(i) % mu * grid_blocks(i) % dT/ &
                                (grid_blocks(i) % dX * grid_blocks(i) % dX)
       
        ! Loop Through Interior Domain (B.C. on i=1 and i=iMax)
        DO j = 2,grid_blocks(i) % iMax-1

          grid_blocks(i) % rhs(j) =  grid_blocks(i) % u(j)  
          grid_blocks(i) % dm1(j) = -grid_blocks(i) % sigma                   ! Lower (i-1)
          grid_blocks(i) % d0(j)  = 1.0_rDef+2.0_rDef*grid_blocks(i) % sigma  ! Main  (i)
          grid_blocks(i) % dp1(j) = -grid_blocks(i) % sigma                   ! Upper (i+1)

        END DO

        ! Effect of B.C. (and removing i=1 and i=iMax from the matrix, since we're specifying
        !                 the values)

        grid_blocks(i) % uNew(grid_blocks(i) % iMin) = grid_blocks(i) % &
                                                       uExact(grid_blocks(i) % iMin)

        grid_blocks(i) % uNew(grid_blocks(i) % iMax) = grid_blocks(i) % &
                                                       uExact(grid_blocks(i) % iMax)
        j = 2
        grid_blocks(i) % rhs(j) = grid_blocks(i) % rhs(j) &
                                - grid_blocks(i) % dm1(j) * grid_blocks(i) % uNew(j-1)

        j = grid_blocks(i) % iMax-1
        grid_blocks(i) % rhs(j) = grid_blocks(i) % rhs(j) &
                                - grid_blocks(i) % dp1(j) * grid_blocks(i) % uNew(j+1)

        ! Thomas Algorithm Solve
        DO j = 3,grid_blocks(i) % iMax-1

          fac = grid_blocks(i) % dm1(j)/grid_blocks(i) % d0(j-1)
          grid_blocks(i) % d0(j) = grid_blocks(i) % d0(j) - (grid_blocks(i) % dp1(j-1)*fac)
          grid_blocks(i) % rhs(j) = grid_blocks(i) % rhs(j) - (grid_blocks(i) % rhs(j-1)*fac)

        END DO

        j = grid_blocks(i) % iMax-1
        grid_blocks(i) % uNew(j) = grid_blocks(i) % rhs(j)/grid_blocks(i) % d0(j)

        ! Backsubstitution
        DO j = grid_blocks(i) % iMax-2,2,-1

          grid_blocks(i) % uNew(j) = (grid_blocks(i) % rhs(j) - grid_blocks(i) % dp1(j) &
                                   *  grid_blocks(i) % uNew(j+1))/grid_blocks(i) % d0(j)

        END DO

      END DO


      RETURN
    END SUBROUTINE CalculateCD2I1Update

END MODULE CD2I


PROGRAM MAIN
  USE block_type
  USE Input_Processing
  USE SetParams
  USE Exact_Solution
  USE CD2E
  USE CD2I
  USE ISO_FORTRAN_ENV
  IMPLICIT NONE

    ! Constant Declarations
    INTEGER, PARAMETER :: rDef           = REAL64, & 
                          num_grids      = 2,      &  
                          num_schemes    = 3,      &  
                          maxItn = 1000000    ! Fail-Safe

    REAL(KIND=rDef), PARAMETER :: x_min   = 0.0_rDef,    &
                                  x_max   = 100.0_rDef,  &
                                  l2Tol  = 1.0e-10_rDef, &
                                  l2Max  = 1.0_rDef,     & ! Stops the code
                                  mu      = 0.1_rDef

    REAL(KIND=rDef), DIMENSION(3), PARAMETER :: sigmaMax = (/ &
                        0.5_rDef,                &  ! E1
                        0.6963233908512847_rDef, &  ! RK
                        1.0e15_rDef              /) ! I1

    CHARACTER(LEN=20) :: InputFileName

    ! Variable Declaration(s)
    INTEGER :: num_blocks,      &
               ios,             &
               i,j,             &
               inFileUnit,      &
               outFileUnit,     &
               outDataFileUnit, &
               nScheme,nG,nStep

    TYPE(grid_block), DIMENSION(:), ALLOCATABLE :: grid_blocks


    CHARACTER(LEN=30) :: outFileName, &
                         outDataFileName

    CHARACTER(LEN=5), DIMENSION(3) :: dType 

    ! Stencil "Names"
    dType(1) = 'CD2E1'
    dType(2) = 'CD2RK'
    dType(3) = 'CD2I1'


    !! MAIN Code Starts Here
    !--------------------------------------------------------------------
      
    !! READ DATA
    !--------------------------------------------------------------------
    DO

     ! Open data file on unit number "in" and read data
     OPEN (15, FILE='input.txt', STATUS="OLD", IOSTAT=ios)

     ! Determine Number of Blocks from Input File
     READ (15,'(I3)') num_blocks
     WRITE(0,'(/"Number of blocks (specified in file):",I3/)') num_blocks

     ! Allocate Appropriate Number of Blocks
     ALLOCATE (grid_blocks(num_blocks))

     ! Read Data from Input and Write Results to Screen
     WRITE(0,'("Read Data:")') 

     DO i = 1,num_blocks
       READ (15,'(4(I3,2X),I3)') grid_blocks(i) % block_number, &
                                         grid_blocks(i) % xMin,         &
                                         grid_blocks(i) % xMax,         &
                                         grid_blocks(i) % iMin,         &
                                         grid_blocks(i) % iMax

      ! Allocate Flow Value Arrays
      ALLOCATE(grid_blocks(i) % u(grid_blocks(i) % iMax),        &
               grid_blocks(i) % x(grid_blocks(i) % iMax),        &
               grid_blocks(i) % uNew(grid_blocks(i) % iMax),     &
               grid_blocks(i) % uExact(grid_blocks(i) % iMax),   &
               grid_blocks(i) % fVec0(grid_blocks(i) % iMax),    &
               grid_blocks(i) % fVecS(grid_blocks(i) % iMax),    &
               grid_blocks(i) % rhs(grid_blocks(i) % iMax),      &
               grid_blocks(i) % dm1(grid_blocks(i) % iMax),      &
               grid_blocks(i) % d0(grid_blocks(i) % iMax),       &
               grid_blocks(i) % dp1(grid_blocks(i) % iMax))     





      ! Write Read Data to Screen
       WRITE (0,'(4(I3,2X),I3)') grid_blocks(i) % block_number, &
                                 grid_blocks(i) % xMin,         &
                                 grid_blocks(i) % xMax,         &
                                 grid_blocks(i) % iMin,         &
                                 grid_blocks(i) % iMax
     END DO

     ! Repeat request if file not opened satisfactorily
     IF (ios == 0) EXIT
       PRINT *,"Unable to open file -- please try again"
    END DO  ! Open File
    CLOSE(15) ! Input
    
    ! Sort Data
    CALL Sort_Data(grid_blocks,num_blocks)

    ! Calculate/Set Parameters
    CALL GetParams(grid_blocks,num_blocks)


    ! EXACT SOLUTION
    !--------------------------------------------------------------------

    ! Initialize Flow Values
    DO i = 1,num_blocks
      DO j = 1,grid_blocks(i)%iMax
        grid_blocks(i) % u(j)      = 0.0_rDef
        grid_blocks(i) % x(j)      = 0.0_rDef 
        grid_blocks(i) % uNew(j)   = 0.0_rDef
        grid_blocks(i) % uExact(j) = 0.0_rDef
        grid_blocks(i) % fVec0(j)  = 0.0_rDef
        grid_blocks(i) % fVecS(j)  = 0.0_rDef
        grid_blocks(i) % rhs(j)    = 0.0_rDef 
        grid_blocks(i) % dm1(j)    = 0.0_rDef 
        grid_blocks(i) % d0(j)     = 0.0_rDef 
        grid_blocks(i) % dp1(j)    = 0.0_rDef 
      END DO

    ! Apply GLOBAL IC
    grid_blocks(i) % uInitial     = 2.0_rDef   
    END DO

    ! Create Space Domain 
    DO i = 1,num_blocks
      DO j = 1,grid_blocks(i) % iMax
        grid_blocks(i) % x(j) = grid_blocks(i) %  xMin + REAL(j-1,rDef) * grid_blocks(i) % dX
      END DO
    END DO

    ! Calculate Exact Solution
    CALL GetExactSoln(grid_blocks,num_blocks)


    ! NUMERICAL APPROXIMATION
    !---------------------------------------------------------------------

    ! Initialize Data For Each Grid Block

    ! Loop Through Each Block
    DO i = 1,num_blocks
      
      ! Loop Through Each Block's Domain
      DO j = 1, grid_blocks(i) % iMax
        grid_blocks(i) % u(j) = grid_blocks(i) % uInitial
      END DO
    END DO

    ! Loop Through Each Problem
    DO nScheme = 1,num_schemes

      WRITE(0,*) nScheme

      WRITE(outFileName, '(A,A,A,I3.3,A)') 'RunData',                       &
                                            dType(nScheme),                 &
                                           'iMax',                          &
                                            grid_blocks(num_blocks) % iMax, &
                                            '.dat'

      OPEN(NEWUNIT = outFileUnit,       &
           FILE    = TRIM(outFileName), &
           FORM    = 'FORMATTED')

      ! Set dT
      grid_blocks(i) % dT = sigmaMax(nScheme)                           &
                          * (grid_blocks(i) % dX * grid_blocks(i) % dX) &
                          /  grid_blocks(i) % mu
                                           
      ! Loop Through Each Block
      DO i = 1,num_blocks

        ! Initialize Step
        nStep = 0
        99 CONTINUE

        ! Increment Step
        nStep = nStep + 1

        ! Calculate Solution(s)
        IF (nScheme == 1) THEN
          CALL GetCD2E1Update(grid_blocks,num_blocks) ! Explicit Scheme

        ELSE IF (nScheme == 2) THEN
          CALL GetCD2RKUpdate(grid_blocks,num_blocks) ! RK Scheme

        ELSE IF (nScheme == 3) THEN
          CALL GetCD2I1Update(grid_blocks,num_blocks) ! Implicit Scheme
        ELSE 
          CONTINUE
        END IF

        ! Calculate L2 Error
        
        ! Initialize L2 Error
        grid_blocks(i) % l2Err = 0.0_rDef

        ! Sum Value of L2 Within Each Block
        DO j = 1,grid_blocks(i) % iMax

          ! Update Flux
          grid_blocks(i) % u(j)   = grid_blocks(i) % uNew(j)

          ! Calculate dU
          grid_blocks(i) % dU     = grid_blocks(i) % u(j) - grid_blocks(i) % uExact(j)

          ! Calculate Numerator of L2 Function
          grid_blocks(i) % l2Err  = grid_blocks(i) % l2Err +                            &
                                   (grid_blocks(i) % dU * grid_blocks(i) % dU)
        END DO ! L2 Loop

        ! Calculate (Actual) L2
        grid_blocks(i) % l2Err = SQRT(grid_blocks(i) % l2Err/REAL(grid_blocks(i)%iMax))

        ! Write Ouput
        WRITE(outFileUnit,*) nStep,grid_blocks(i) % l2Err ! To file
        !WRITE(0,*) nStep,grid_blocks(i) % l2Err           ! To screen

        ! Determine Whether to Continue or Not
        IF ((grid_blocks(i) % l2Err < l2Max) .AND.  &
            (grid_blocks(i) % l2Err > l2Tol) .AND.  &
            (nStep < maxItn)) THEN  ! Keep going
        GO TO 99  
        ELSE  ! Converged
          WRITE(0,'(/"Block: ",I3," converged!"/)') grid_blocks(i) % block_number
        END IF

      END DO  ! num_blocks Loop

      CLOSE(outFileUnit)

      ! Write Results to File
      WRITE(outDataFileName,'(A,A,A,I3.3,A)') 'SolutionData',                  &
                                               dType(nScheme),                 &
                                              'iMax',                          &
                                               grid_blocks(num_blocks) % iMax, &
                                              '.dat'
      OPEN(NEWUNIT = outDataFileUnit,      &
           FILE    = TRIM(outDataFileName),&
           FORM    = 'FORMATTED')

      ! (Actually) Write Results
      DO i = 1,num_blocks
        DO j = grid_blocks(i) % iMin,grid_blocks(i) % iMax
          WRITE(outDataFileUnit,*) grid_blocks(i) % x(j),          &
                                   grid_blocks(i) % uExact(j),     &
                                   grid_blocks(i) % uNew(j)
        END DO
      END DO


      CLOSE(outDataFileUnit) ! Output
    END DO ! nScheme Loop


    ! De-allocate Block Arrays
    DEALLOCATE(grid_blocks)

    STOP
END PROGRAM MAIN
