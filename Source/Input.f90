MODULE Input_Processing
  USE, INTRINSIC :: ISO_FORTRAN_ENV
  USE block_type
  IMPLICIT NONE
  PRIVATE
  PUBLIC Process_Data
  INTEGER, PARAMETER :: rDef = REAL64

 INTERFACE Process_Data 
   MODULE PROCEDURE Read_Data
 END INTERFACE Process_Data

  ! This module is for processing the input file
  ! supplied for the multi-block grid generator

  ! NOTE: The input file format is as follows:

!|  Block Number  | xMin  | xMax  | iMin  | iMax  |
!|----------------|-------|-------|-------|-------|


  CONTAINS

    ! Subroutine to sort the Input File Data 
    SUBROUTINE Read_Data(InputFileName,  &
                         blocks,         &
                         num_blocks,     &
                         global_iMax)

    ! Dummy Variable Declaration(s)
    CHARACTER(LEN=20), INTENT(IN) :: InputFileName

    TYPE(grid_block), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: blocks
    INTEGER, INTENT(OUT) :: global_iMax,  &
                            num_blocks

    ! Local Varaiable Declaration(s)
    INTEGER :: i,j                    ! DO-loop integers

    TYPE(grid_block), DIMENSION(:), ALLOCATABLE :: temp   ! For sorting


    !! Reading Data
    !------------------------------------------------------------------
    ! Determine Number of Blocks from Input File
    READ (15,'(I3)') num_blocks
    WRITE(0,'("Number of blocks (specified in file):",I3/)') num_blocks
    WRITE(0,'("Read Data:")') 

    ! Allocate Appropriate Number of Blocks
    ALLOCATE (blocks(num_blocks), &         ! Allocate for "temp" array
              temp(num_blocks))             ! too. Will use this for
                                            ! sorting data
                                            

    ! Read Data from Input and Write Results to Screen
    DO i = 1,num_blocks
      READ (15,'(4(I3,2X),I3)') blocks(i) % block_number, &
                                blocks(i) % xMin,         &
                                blocks(i) % xMax,         &
                                blocks(i) % iMin,         &
                                blocks(i) % iMax

      WRITE (0,'(4(I3,2X),I3)') blocks(i) % block_number, &
                                blocks(i) % xMin,         &
                                blocks(i) % xMax,         &
                                blocks(i) % iMin,         &
                                blocks(i) % iMax
    END DO

    !! Sorting Data
    !------------------------------------------------------------------

    ! Transfer Data to "Temp" Array
    DO i = 1,num_blocks
      
      temp(i) = blocks(i) 

    END DO                    

    PRINT *,""

    !Sanity Check
    ! Write Original Array to Screen
!   WRITE(0,'("Original Array (blocks):")')
!
!   ! Initialize Values Array Values to Zero (For Security)
!   blocks % block_number = 0
!   blocks % xMin         = 0
!   blocks % xMax         = 0
!   blocks % iMin         = 0
!   blocks % iMax         = 0
!
!   DO i = 1,num_blocks
!
!     WRITE (0,*) blocks(i) % block_number, &
!                 blocks(i) % xMin,         &
!                 blocks(i) % xMax,         &
!                 blocks(i) % iMin,         &
!                 blocks(i) % iMax
!
!   END DO                    
!
!   Print *,""
!
!   WRITE(6,'("New Array (temp):")')
!   DO i = 1,num_blocks
!     WRITE (0,*) temp(i) % block_number, &
!                 temp(i) % xMin,         &
!                 temp(i) % xMax,         &
!                 temp(i) % iMin,         &
!                 temp(i) % iMax
!   END DO
!
!   PRINT *,""

    ! (Actually) Sort Data 
    
    ! Sorting Algorithm:
    
    ! 1.) Utilize intrinsic FINDLOC to determine location
    !     where block number of temporary matches the 
    !     incremental i integer
    
    ! 2.) Take the values of the "temporary" array at this 
    !     location
    
    ! 3.) Set value of block array equal to these values
    
!   WRITE (0,'("Sorted Data:")')
    DO i = 1,num_blocks
      blocks(i) = temp(FINDLOC(temp % block_number,i,1)) 

    END DO

    ! Calculate Global Parameters

    ! iMax
    global_iMax = blocks(1) % iMax
    DO i=2,num_blocks-1
      global_iMax = global_iMax + blocks(i) % iMax
    END DO
    WRITE (0,'("Global iMax: ",I3,/)') global_iMax

    ! Allocate and Initialize Flow Variables
!   WRITE(0,'("Block ",3X,"i",6X,"u",8X,"uNew",6X,"uExact",8X,"x")')

    DO i = 1,num_blocks
      ALLOCATE(blocks(i) % u(blocks(i)%iMax),      &
               blocks(i) % uNew(blocks(i)%iMax),   &
               blocks(i) % uExact(blocks(i)%iMax), &
               blocks(i) % x(blocks(i)%iMax))


      DO j = 1,blocks(i)%iMax
        blocks(i) % u(j)      = 0.0_rDef
        blocks(i) % uNew(j)   = 0.0_rDef
        blocks(i) % uExact(j) = 0.0_rDef
        blocks(i) % x(j)      = 0.0_rDef 

      ! SANITY CHECK: Write to Screen
!     WRITE(0,'(2I5,3(F9.3,2X),F9.3)') i,j,                            &  
!                                      blocks(i) % u(j),      &     
!                                      blocks(i) % uNew(j),   &
!                                      blocks(i) % uExact(j), &
!                                      blocks(i) % x(j)

      END DO

    END DO


    END SUBROUTINE Read_Data
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

    SUBROUTINE Calculate_Parameters(unsolved_blocks,  &
                                    sigma,            &
                                    dX,               &
                                    mu,               &
                                    dT)

      ! Dummy Argument Declarations
      TYPE(grid_block), DIMENSION(:), INTENT(IN) :: unsolved_blocks


      REAL(KIND=rDef), DIMENSION(:), ALLOCATABLE,  INTENT(OUT) :: sigma,  &
                                                                  dX,     &
                                                                  mu,     & 
                                                                  dT

      ! Local Variables
      INTEGER :: num_blocks,  &
                 i

      !num_blocks = SIZE(unsolved_blocks,1)

      num_blocks = 1

      ! Allocate Memory
      ALLOCATE(sigma(num_blocks), &
                  dx(num_blocks), &
                  mu(num_blocks), &
                  dT(num_blocks))

      DO i = 1,num_blocks

        dX(i) =     (unsolved_blocks(i) % xMax - unsolved_blocks(1) % xMin)/ &
                REAL(unsolved_blocks(i) % iMax - 1,rDef)

        mu(i) = 1.0_rDef/10.0_rDef

        dT(i) = 0.5_rDef*(dX(i)*dX(i))/mu(i)

        sigma(i) = mu(i)*dT(i)/(dX(i)*dX(i))

        WRITE(6,'("Calculated Parameters: ",/,  &
          "Block: ",I3,/,                       &
          "sigma: ",F5.3,2X,"dX: ",F5.3,2X,"mu: ",F5.3,2X,"dT: ",F5.3)') unsolved_blocks(i) % block_number,  &
                                                                         sigma(i),       &
                                                                         dX(i),          &
                                                                         mu(i),          &
                                                                         dT(i)
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
    SUBROUTINE CalculateCD2E1Update(sigma,dX,mu,dT,   &
                                    unsolved_blocks,  &
                                    solved_blocks)

    ! Dummy Variables
    REAL(KIND=rDef), DIMENSION(:), INTENT(IN) :: sigma, dX, mu, dT
    TYPE(grid_block), DIMENSION(:), INTENT(IN) :: unsolved_blocks

    TYPE(grid_block), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: solved_blocks


    ! Local Variables
    REAL(KIND=rDef), DIMENSION(:), ALLOCATABLE :: x,  &
                                                  uExact, &
                                                  u,  &
                                                  uNew


    INTEGER :: i,j,   &
               iMin,  &
               iMax

    REAL(KIND=rDef) :: xMin,  &
                       xMax
    
    ! Beginning with a single block...

    ! Allocate Memory
    ALLOCATE(solved_blocks(1))

    ! Set Local Values
    iMin = unsolved_blocks(1)%iMin
    iMax = unsolved_blocks(1)%iMax

    ALLOCATE(x(iMax),       &
             uExact(iMax),  &
             u(iMax),       &
             uNew(iMax))

    xMin = unsolved_blocks(1)%xMin
    xMax = unsolved_blocks(1)%xMax


    solved_blocks(1)%uMin     = 1.0_rDef
    solved_blocks(1)%uMax     = 3.0_rDef
    solved_blocks(1)%uInitial = 2.0_rDef


    ! Flow Solution (Starting with 1 block)

!    OPEN(55, FILE='test.dat', STATUS='NEW')
!

    ! Initialize and Map Exact Solution
!   DO i = iMin,iMax
!     u(i) = solved_blocks(1) % uInitial
!     x(i) = xMin + REAL(i-1,rDef)*dX(1)
!     uExact(i) = solved_blocks(1) % uMin + (solved_blocks(1) % uMax  &
!                                         -  solved_blocks(1) % uMin) &
!                                         * (x(i)-xMin)/(xMax - xMin)
!   END DO
!
!   ! Apply B.C.
!   uNew(1) = solved_blocks(1) % uMin
!   uNew(iMax) = solved_blocks(1) % uMax
!
!   ! Solving for Interior Domain
!   DO i = 2,iMax-1
!    uNew(i) = u(i) + sigma(1)*(u(i+1)   &
!                   -  2.0_rDef*u(i)     &
!                   +           u(i-1))
!   END DO
!
!
!
!
!   ! Write Results
!   DO i = iMin,iMax
!      WRITE(55,*) x(i),uExact(i),uNew(i)
!   END DO
!
!   CLOSE(55)

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
                          max_iterations = 1000000    ! Fail-Safe

    REAL(KIND=rDef), PARAMETER :: x_min   = 0.0_rDef,     &
                                  x_max   = 100.0_rDef,   &
                                  l2_tol  = 1.0e-10_rDef, &
                                  l2_max  = 1.0_rDef,     & ! Stops the code
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
               nG,nStep,        &
               iMax

   REAL(KIND=rDef), DIMENSION(:), ALLOCATABLE :: sigma      
!                                                 dX,        &
!                                                 mu,        &
!                                                 dT,        &
!                                                 uMin,      & 
!                                                 uMax,      &
!                                                 uInitial,  & 
!                                                 l2Err,     &
!                                                 dU
!
!   REAL(KIND=rDef), DIMENSION(:), ALLOCATABLE :: u,      &
!                                                 uNew,   &
!                                                 uExact, &
!                                                 x
!
    TYPE(grid_block), DIMENSION(:), ALLOCATABLE :: unsolved_blocks, &
                                                     solved_blocks, &
                                                     block_parameters

!
!   CHARACTER(LEN=30) :: OutputFileName, &
!                        OutputDataFileName
!
!   CHARACTER(LEN=5), DIMENSION(3) :: dType 
!
!   integer :: test
!
!   ! Stencil "Names"
!   dType(1) = 'CD2E1'
!   dType(2) = 'CD2I1'
!   dType(3) = 'CD2RK'
!
!   ! Boundary and Initial Condtions
!   uMin     = 1.0_rDef
!   uMax     = 3.0_rDef
!   uInitial = 2.0_rDef


    !! MAIN Code Starts Here
    !--------------------------------------------------------------------
      
    DO
     PRINT *,"Please give name of data file"
     READ *,InputFileName

     ! Open data file on unit number "in" and call subroutine "Read_Data"
     OPEN (15, FILE=InputFileName, STATUS="OLD", IOSTAT=ios)

     ! Repeat request if file not opened satisfactorily
     IF (ios == 0) EXIT
       PRINT *,"Unable to open file -- please try again"
    END DO

    ! Read data
    CALL Process_Data(InputFileName,unsolved_blocks,num_blocks,iMax)
    CLOSE (15)

    ! Set Parameters
!   CALL GetParams(unsolved_blocks,  &
!                  sigma,            &
!                  dX,               &
!                  mu,               &
!                  dT)
!
!
!
!   ! Calculate Solution(s)
!   CALL GetCD2E1Update(sigma,dX,mu,dT,   &
!                       unsolved_blocks,  &
!                       solved_blocks)
!
!   ! L2 Calculation
!   l2Err = 0.0_rDef
!   DO i = 1,iMax
!     u(i)  = uNew(i)
!     dU    = u(i) - uNew(i)
!     l2Err = l2Err + (dU*dU)
!   END DO
!   l2Err = SQRT(l2Err/REAL(iMax,rDef))
  

END PROGRAM MAIN



  
