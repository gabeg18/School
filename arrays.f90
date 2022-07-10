MODULE grid_block_type
  IMPLICIT NONE
  SAVE

  ! This module defines a derived data type for each block of a multi-block
  ! grid generator

  TYPE grid_block

    INTEGER :: lower_boundary, &    ! Grid Block Lower Boundary
               upper_boundary, &    ! Grid Block Upper Boundary
               block_number         ! Respective Block Number

    INTEGER, DIMENSION (:), ALLOCATABLE :: num_grid_points ! Number of points
                                                           ! per each grid 
                                                           ! block

    LOGICAL :: boundary_condition   ! T = Boundary Block
                                    ! F = Interior Block

  END TYPE grid_block

END MODULE grid_block_type


MODULE multiblock_grid
  USE grid_block_type
  IMPLICIT NONE
  PRIVATE
  PUBLIC gen_multiblock_grid


  CONTAINS
  SUBROUTINE gen_multiblock_grid(num_blocks, &
                                  num_points, &
                                  blocks)

  ! Dummy Variable Declaration(s)
  INTEGER, INTENT (IN) :: num_blocks,  & 
                          num_points

  ! Utilize block defined data type
  TYPE (grid_block), DIMENSION (:), ALLOCATABLE, INTENT (OUT) :: blocks
  
  ! Local Variable Declaration(s)
  INTEGER :: i, j  ! Count integers


  ! Allocate num_blocks number of grid blocks
  ALLOCATE ( blocks (num_blocks) )

  ! Allocate num_blocks number of grid blocks, each containing num_points
  ! number of grid points
  DO i=1, num_blocks
    ALLOCATE ( blocks (i) % num_grid_points (num_points) )
  END DO

  ! Block Array Creation
  DO i=1, num_blocks                              

    ! Block Array Population
    DO j=1,num_points                             

      ! Fill First Block
      IF (i==1) THEN
      blocks (i) % num_grid_points (j) = (i-1)+j      ! Fill with num_points number of
                                                      ! grid points

      ! Determine Lower and Upper Boundaries of First Block
      blocks (i) % lower_boundary = blocks(i) % num_grid_points(1)            ! Lower
      blocks (i) % upper_boundary = blocks(i) % num_grid_points(num_points)   ! Upper
      
      ! Number Block  
      blocks (i) % block_number = i

      ! Fill Successive Blocks
      ELSE
      blocks(i) % num_grid_points(j) = (blocks(i-1) % num_grid_points(num_points))  &
                                        + (j-1)       ! Reference previous block and 
                                                      ! fill with successive points
                                                      ! (share common point at boundaries)

      ! Determine Lower and Upper Boundaries of Each Block
      blocks (i) % lower_boundary = blocks(i) % num_grid_points(1)            ! Lower
      blocks(i) % upper_boundary = blocks(i) % num_grid_points(num_points)    ! Upper

      ! Number Block
      blocks (i) % block_number = i

      END IF
    ! End Block Population
    END DO               

    ! Write Results to Screen
    WRITE (*, *)            
    WRITE (*,'("Block:",I3)') blocks(i) % block_number
    WRITE (*, '(10I5)') blocks (i) % num_grid_points
    WRITE (*, '("Lower Boundary:",I3,T25,"Upper Boundary:",I3)') blocks(i) % lower_boundary, &
                                                                 blocks(i) % upper_boundary

  ! End Block Array Creation
  END DO

  END SUBROUTINE gen_multiblock_grid
END MODULE multiblock_grid

PROGRAM arrays
  USE grid_block_type
  USE multiblock_grid

  IMPLICIT NONE

  LOGICAL, PARAMETER :: Yes = .TRUE., &
                        No  = .FALSE. 
  INTEGER :: num_blocks,  &
             num_points

  ! Utilize block defined data type
  TYPE (grid_block), DIMENSION (:), ALLOCATABLE :: blocks

  LOGICAL :: boundary_blocks

  ! Data Input
  PRINT *,"Specify number of blocks"
  READ *,num_blocks

  PRINT *,"Specify number of grid points per block"
  READ *,num_points

  PRINT *,"Specify the boundary blocks"
  READ *,boundary_blocks

  CALL gen_multiblock_grid(num_blocks, &
                            num_points, &
                            blocks)

END PROGRAM arrays
