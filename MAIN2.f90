PROGRAM MAIN
  
  ! Program to test the use of the block_gen module

  USE block_gen

  IMPLICIT NONE

  INTEGER :: iMin = 1,  &
             iMax = 101

  REAL, PARAMETER :: x_min   =   0.0, &
                     x_max   = 100.0, &
                     delta_x =   1.0

  REAL, DIMENSION(6) :: block_min, block_max

  CALL blocks(x_min,x_max,delta_x,block_min,block_max)

END PROGRAM MAIN
