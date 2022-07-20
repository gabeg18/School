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

    SUBROUTINE Calculate_Parameters(unsolved_blocks,num_blocks)

      ! Dummy Argument Declarations
      TYPE(grid_block), DIMENSION(:), INTENT(INOUT) :: unsolved_blocks
      INTEGER, INTENT(IN) :: num_blocks

      ! Local Variables
      INTEGER :: i

      !num_blocks = 1

      WRITE(0,'(/"Calculated Parameters: "/)')

      DO i = 1,num_blocks

        unsolved_blocks(i) % dX =     (unsolved_blocks(i) % xMax - unsolved_blocks(i) % xMin)/ &
                                  REAL(unsolved_blocks(i) % iMax - 1,rDef)

        unsolved_blocks(i) % mu = 1.0_rDef/10.0_rDef

        unsolved_blocks(i) % dT = 0.5_rDef*(unsolved_blocks(i) % dX * unsolved_blocks(i) % dX)/ &
                                            unsolved_blocks(i) % mu

        unsolved_blocks(i) % sigma = unsolved_blocks(i) % mu * unsolved_blocks(i) % dT/ &
                                    (unsolved_blocks(i) % dX * unsolved_blocks(i) % dX)

        WRITE(0,'("Block: ",I3/,                 &
          "sigma: ",F5.3,2X,                     &
          "dX: ",F5.3,2X,                        &
          "mu: ",F5.3,2X,                        &
          "dT: ",F5.3/)')  unsolved_blocks(i) % block_number,  &
                           unsolved_blocks(i) % sigma,         &
                           unsolved_blocks(i) % dX,            &
                           unsolved_blocks(i) % mu,            &
                           unsolved_blocks(i) % dT
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
    SUBROUTINE CalculateCD2E1Update(unsolved_blocks,solved_blocks)

    ! Dummy Variables
    TYPE(grid_block), DIMENSION(:), INTENT(IN) :: unsolved_blocks
    TYPE(grid_block), DIMENSION(:), INTENT(OUT) ::   solved_blocks

    ! Local Variables
    INTEGER :: i,j,   &
               iMin,  &
               iMax,  &
               num_blocks

    REAL(KIND=rDef) :: xMin,  &
                       xMax
    
    num_blocks = 1

    ! Beginning with a single block...

    DO i = 1,num_blocks

      ! Allocate Memory
      ALLOCATE(solved_blocks(i) % u(unsolved_blocks(i)%iMax),        &
               solved_blocks(i) % uNew(unsolved_blocks(i)%iMax),     &
               solved_blocks(i) % uExact(unsolved_blocks(i)%iMax),   &
               solved_blocks(i) % x(unsolved_blocks(i)%iMax))
             
      ! "Transfer" Data to "Solved" Blocks
      solved_blocks(i) % block_number = unsolved_blocks(i) % block_number
      solved_blocks(i) % xMin         = unsolved_blocks(i) % xMin
      solved_blocks(i) % xMax         = unsolved_blocks(i) % xMax
      solved_blocks(i) % iMin         = unsolved_blocks(i) % iMin
      solved_blocks(i) % iMax         = unsolved_blocks(i) % iMax
      solved_blocks(i) % block_number = unsolved_blocks(i) % block_number
      solved_blocks(i) % sigma        = unsolved_blocks(i) % sigma
      solved_blocks(i) % dX           = unsolved_blocks(i) % dX
      solved_blocks(i) % mu           = unsolved_blocks(i) % mu
      solved_blocks(i) % dT           = unsolved_blocks(i) % dT

      ! Initialize Flow Variables
      DO j = 1,solved_blocks(i)%iMax
        solved_blocks(i) % u(j)      = 0.0_rDef
        solved_blocks(i) % uNew(j)   = 0.0_rDef
        solved_blocks(i) % uExact(j) = 0.0_rDef
        solved_blocks(i) % x(j)      = 0.0_rDef 
      END DO

      solved_blocks(i) % uMin     = 1.0_rDef
      solved_blocks(i) % uMax     = 3.0_rDef
      solved_blocks(i) % uInitial = 2.0_rDef

      ! Create Space Domain and Map Exact Solution
      DO j = 1,solved_blocks(i) % iMax

        solved_blocks(i) % x(j) = solved_blocks(i) %  xMin + REAL(j-1,rDef) * solved_blocks(i) % dX
        solved_blocks(i) % uExact(j) = solved_blocks(i) % uMin + (solved_blocks(i) % uMax   &
                                                               -  solved_blocks(i) % uMin)  &
                                     * (solved_blocks(i) % x(j) - solved_blocks(i) % xMin)/ &
                                       (solved_blocks(i) % xMax - solved_blocks(i) % xMin)
      END DO
      ! Set Local Values
      iMin = solved_blocks(i) % iMin
      iMax = solved_blocks(i) % iMax

      xMin = solved_blocks(i) % xMin
      xMax = solved_blocks(i) % xMax

    END DO


    ! Flow Solution (Starting with 1 block)

    ! Apply B.C.
    solved_blocks(1) % uNew(1) = solved_blocks(1) % uMin
    solved_blocks(1) % uNew(iMax) = solved_blocks(1) % uMax

    ! Solving for Interior Domain
    DO i = 2,iMax-1
      solved_blocks(1) % uNew(i) = solved_blocks(1) % u(i)   + solved_blocks(1) % sigma *  &
                                  (solved_blocks(1) % u(i+1)                               &
                       -  2.0_rDef*solved_blocks(1) % u(i)                                 &
                       +           solved_blocks(1) % u(i-1))
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

    TYPE(grid_block), DIMENSION(:), ALLOCATABLE :: unsolved_blocks, &
                                                     solved_blocks


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
     ALLOCATE (unsolved_blocks(num_blocks), &    ! For performing calculations
                  solved_blocks(num_blocks))      ! Actual solved blocks      

     WRITE(0,'("Read Data:")') 


     ! Read Data from Input and Write Results to Screen
     DO i = 1,num_blocks
       READ (15,'(4(I3,2X),I3)') unsolved_blocks(i) % block_number, &
                                 unsolved_blocks(i) % xMin,         &
                                 unsolved_blocks(i) % xMax,         &
                                 unsolved_blocks(i) % iMin,         &
                                 unsolved_blocks(i) % iMax

      ALLOCATE(solved_blocks(i) % u(unsolved_blocks(i)%iMax),        &
               solved_blocks(i) % uNew(unsolved_blocks(i)%iMax),     &
               solved_blocks(i) % uExact(unsolved_blocks(i)%iMax),   &
               solved_blocks(i) % x(unsolved_blocks(i)%iMax))

       WRITE (0,'(4(I3,2X),I3)') unsolved_blocks(i) % block_number, &
                                 unsolved_blocks(i) % xMin,         &
                                 unsolved_blocks(i) % xMax,         &
                                 unsolved_blocks(i) % iMin,         &
                                 unsolved_blocks(i) % iMax
     END DO

     ! Repeat request if file not opened satisfactorily
     IF (ios == 0) EXIT
       PRINT *,"Unable to open file -- please try again"
    END DO

    
    CALL Sort_Data(unsolved_blocks,num_blocks)
    !CLOSE (15)

    ! Set Parameters
    CALL GetParams(unsolved_blocks,num_blocks)

    ! Initialize Data
    DO i = 1, solved_blocks(1) % iMax
      solved_blocks(1) % u(i) = solved_blocks(1) % uInitial
    END DO

    nStep = 0
    99 CONTINUE

    ! Calculate Solution(s)
    nStep = nStep + 1

    CALL GetCD2E1Update(unsolved_blocks,solved_blocks)

    ! Calculate l2Err
    solved_blocks(1) % l2Err = 0.0_rDef
    DO i = 1,solved_blocks(1) % iMax
      solved_blocks(1) % u(i)   = solved_blocks(1) % uNew(i)
      solved_blocks(1) % dU     = solved_blocks(1) % u(i) - solved_blocks(1) % uExact(i)
      solved_blocks(1) % l2Err  = solved_blocks(1) % l2Err + &
                                 (solved_blocks(1) % dU * solved_blocks(1) % dU)
    END DO

    solved_blocks(1) % l2Err = SQRT(solved_blocks(1) % l2Err/REAL(solved_blocks(1)%iMax))
!
!   WRITE(0,*) nStep,solved_blocks(1) % l2Err
!
!   IF ((solved_blocks(1) % l2Err < l2Max) .AND.  &
!      (solved_blocks(1) % l2Err > l2Tol) .AND.  &
!      (nStep < maxItn)) THEN  ! Keep going
!    GO TO 99
!  ELSE
!    CONTINUE  ! Done
!  END IF

    OPEN(55, FILE='test.dat', STATUS='NEW')

    ! Write Results
    DO i = solved_blocks(1) % iMin,solved_blocks(1) % iMax
      WRITE(55,*) solved_blocks(1) % x(i),      &
                  solved_blocks(1) % uExact(i), &
                  solved_blocks(1) % uNew(i)
    END DO

     CLOSE(55)

END PROGRAM MAIN



  
