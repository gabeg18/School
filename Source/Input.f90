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

        grid_blocks(i) % sigma = grid_blocks(i) % mu * grid_blocks(i) % dT/ &
                                    (grid_blocks(i) % dX * grid_blocks(i) % dX)

        WRITE(0,'("Block: ",I3/,                 &
          "sigma: ",F5.3,2X,                     &
          "dX: ",F5.3,2X,                        &
          "mu: ",F5.3,2X,                        &
          "dT: ",F5.3/)')  grid_blocks(i) % block_number,  &
                           grid_blocks(i) % sigma,         &
                           grid_blocks(i) % dX,            &
                           grid_blocks(i) % mu,            &
                           grid_blocks(i) % dT
      END DO

    END SUBROUTINE Calculate_Parameters

END MODULE SetParams



MODULE CD2E
  USE, INTRINSIC :: ISO_FORTRAN_ENV
  USE block_type
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: GetCD2E1Update

  INTEGER, PARAMETER :: rDef = REAL64

 INTERFACE GetCD2E1Update 
   MODULE PROCEDURE CalculateCD2E1Update
 END INTERFACE

  CONTAINS 
    SUBROUTINE CalculateCD2E1Update(grid_blocks)

    ! Dummy Variables
    TYPE(grid_block), DIMENSION(:), INTENT(INOUT) :: grid_blocks

    ! Local Variables
    INTEGER :: i,j,   &
               num_blocks

    REAL(KIND=rDef) :: xMin,  &
                       xMax
    
    num_blocks = 1

    ! Beginning with a single block...

    ! Exact Solution

    DO i = 1,num_blocks
!
!     ! Allocate Memory
!     ALLOCATE(solved_blocks(i) % u(solved_blocks(i)%iMax),        &
!              solved_blocks(i) % uNew(solved_blocks(i)%iMax),     &
!              solved_blocks(i) % uExact(solved_blocks(i)%iMax),   &
!              solved_blocks(i) % x(solved_blocks(i)%iMax))


    ! Flow Solution (Starting with 1 block)

      ! Apply B.C.
      grid_blocks(i) % uNew(grid_blocks(i) % iMin) = grid_blocks(i) % uMin
      grid_blocks(i) % uNew(grid_blocks(i) % iMax) = grid_blocks(i) % uMax

      ! Solving for Interior Domain
      DO j = 2,(grid_blocks(i) % iMax) - 1
        grid_blocks(i) % uNew(j) = grid_blocks(i) % u(j) + grid_blocks(i) % sigma *  &
                                  (grid_blocks(i) % u(j+1)                               &
                       -  2.0_rDef*grid_blocks(i) % u(j)                                 &
                       +           grid_blocks(i) % u(j-1))
      END DO
    END DO

    END SUBROUTINE CalculateCD2E1Update

END MODULE CD2E


PROGRAM MAIN
  USE block_type
  USE Input_Processing
  USE SetParams
  USE CD2E
  USE ISO_FORTRAN_ENV
  IMPLICIT NONE

    ! Constant Declarations
    INTEGER, PARAMETER :: rDef           = REAL64, &  ! Precision
                          num_grids      = 2,      &  ! Grids
                          num_schemes    = 3,      &  ! Stencils
                          maxItn = 1000000    ! Fail-Safe

    REAL(KIND=rDef), PARAMETER :: x_min   = 0.0_rDef,     &
                                  x_max   = 100.0_rDef,   &
                                  l2Tol  = 1.0e-10_rDef, &
                                  l2Max  = 1.0_rDef,     & ! Stops the code
                                  nu      = 0.1_rDef

    REAL(KIND=rDef), DIMENSION(3), PARAMETER :: sigmaMax = (/ &
                           0.5_rDef, &              ! E1
                        1.0e15_rDef, &              ! I1
                        0.6963233908512847_rDef /)  ! RK

    CHARACTER(LEN=20) :: InputFileName

    ! Variable Declaration(s)
    INTEGER :: num_blocks,      &
               ios,             &
               i,j,             &
               outFileUnit,     &
               outDataFileUnit, &
               nG,nStep

    TYPE(grid_block), DIMENSION(:), ALLOCATABLE :: grid_blocks


    CHARACTER(LEN=30) :: OutputFileName, &
                        OutputDataFileName

    CHARACTER(LEN=5), DIMENSION(3) :: dType 

    integer :: test

    ! Stencil "Names"
    dType(1) = 'CD2E1'
    dType(2) = 'CD2I1'
    dType(3) = 'CD2RK'


    !! MAIN Code Starts Here
    !--------------------------------------------------------------------
      
    DO
     !PRINT *,"Please give name of data file"
     !READ *,InputFileName

     ! Open data file on unit number "in" and read data
     OPEN (15, FILE='input.txt', STATUS="OLD", IOSTAT=ios)

     !! Reading Data
     !------------------------------------------------------------------
     ! Determine Number of Blocks from Input File
     READ (15,'(I3)') num_blocks
     WRITE(0,'(/"Number of blocks (specified in file):",I3/)') num_blocks

     ! Allocate Appropriate Number of Blocks
     ALLOCATE (grid_blocks(num_blocks))      ! Actual solved blocks      

     WRITE(0,'("Read Data:")') 


     ! Read Data from Input and Write Results to Screen
     DO i = 1,num_blocks
       READ (15,'(4(I3,2X),I3)') grid_blocks(i) % block_number, &
                                 grid_blocks(i) % xMin,         &
                                 grid_blocks(i) % xMax,         &
                                 grid_blocks(i) % iMin,         &
                                 grid_blocks(i) % iMax

      ! Allocate Flow Value Arrays
      ALLOCATE(grid_blocks(i) % u(grid_blocks(i)%iMax),        &
               grid_blocks(i) % uNew(grid_blocks(i)%iMax),     &
               grid_blocks(i) % uExact(grid_blocks(i)%iMax),   &
               grid_blocks(i) % x(grid_blocks(i)%iMax))

      ! Initialize Flow Values
      DO j = 1,grid_blocks(i)%iMax
        grid_blocks(i) % u(j)      = 0.0_rDef
        grid_blocks(i) % uNew(j)   = 0.0_rDef
        grid_blocks(i) % uExact(j) = 0.0_rDef
        grid_blocks(i) % x(j)      = 0.0_rDef 
      END DO

      ! Set Min, Max, and Initial Conditions
      grid_blocks(i) % uMin     = 1.0_rDef
      grid_blocks(i) % uMax     = 3.0_rDef
      grid_blocks(i) % uInitial = 2.0_rDef

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
    END DO

    
    CALL Sort_Data(grid_blocks,num_blocks)
    !CLOSE (15)

    ! Set Parameters
    CALL GetParams(grid_blocks,num_blocks)

    ! Exact Solution
    ! Create Space Domain and Map Exact Solution
    DO i = 1,num_blocks
      DO j = 1,grid_blocks(i) % iMax

        grid_blocks(i) % x(j) = grid_blocks(i) %  xMin + REAL(j-1,rDef) * grid_blocks(i) % dX
        grid_blocks(i) % uExact(j) = grid_blocks(i) % uMin + (grid_blocks(i) % uMax   &
                                                               -  grid_blocks(i) % uMin)  &
                                     * (grid_blocks(i) % x(j) - grid_blocks(i) % xMin)/ &
                                       (grid_blocks(i) % xMax - grid_blocks(i) % xMin)
      END DO
    END DO

    ! Initialize Data
    DO i = 1, grid_blocks(1) % iMax
      grid_blocks(1) % u(i) = grid_blocks(1) % uInitial
    END DO

    nStep = 0
    99 CONTINUE

    ! Calculate Solution(s)
    nStep = nStep + 1

    CALL GetCD2E1Update(grid_blocks)

    ! Calculate l2Err
    grid_blocks(1) % l2Err = 0.0_rDef
    DO i = 1,grid_blocks(1) % iMax
      grid_blocks(1) % u(i)   = grid_blocks(1) % uNew(i)
      grid_blocks(1) % dU     = grid_blocks(1) % u(i) - grid_blocks(1) % uExact(i)
      grid_blocks(1) % l2Err  = grid_blocks(1) % l2Err + &
                               (grid_blocks(1) % dU * grid_blocks(1) % dU)
    END DO

    grid_blocks(1) % l2Err = SQRT(grid_blocks(1) % l2Err/REAL(grid_blocks(1)%iMax))

    WRITE(0,*) nStep,grid_blocks(1) % l2Err
!
!   IF ((solved_blocks(1) % l2Err < l2Max) .AND.  &
!      (solved_blocks(1) % l2Err > l2Tol) .AND.  &
!      (nStep < maxItn)) THEN  ! Keep going
!    GO TO 99
!  ELSE
!    CONTINUE  ! Done
!  END IF
!
    OPEN(55, FILE='test.dat', STATUS='NEW')
!
    ! Write Results
    DO i = grid_blocks(1) % iMin,grid_blocks(1) % iMax
      WRITE(55,*) grid_blocks(1) % x(i),      &
                  grid_blocks(1) % uExact(i), &
                  grid_blocks(1) % uNew(i)
    END DO
!
    CLOSE(55)

END PROGRAM MAIN



  
