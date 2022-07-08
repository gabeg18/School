PROGRAM MAIN
  
  ! Program to test the use of the mblock_gen and block_data modules

  USE, INTRINSIC :: ISO_FORTRAN_ENV
  !USE block_data
  USE mblock_gen

  IMPLICIT NONE

  INTEGER, PARAMETER :: iMin       = 1,      & 
                        iMax       = 101

  INTEGER :: num_blocks

  REAL(KIND=rDef), PARAMETER :: x_min   =   0.0, &
                                x_max   = 100.0, &
                                delta_x =   1.0

  REAL(KIND=rDef), DIMENSION(:), ALLOCATABLE :: block_min, & 
                                                block_max


  ! Data Input
  PRINT *,"Please input the number of blocks you wish to use"
  READ *,num_blocks

  ALLOCATE( &
    block_min(num_blocks), &
    block_max(num_blocks))

  ! Call Multi-Block Generator from mblock_gen Module
  CALL blocks(x_min,      &
              x_max,      &
              delta_x,    &
              num_blocks, &
              block_min,  &
              block_max)

  DEALLOCATE( &
    block_min,&
    block_max)

END PROGRAM MAIN
