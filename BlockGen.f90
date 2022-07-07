MODULE block_gen
  USE, INTRINSIC :: ISO_FORTRAN_ENV
  IMPLICIT NONE

  INTEGER, PARAMETER :: rDef = REAL64

  CONTAINS

  SUBROUTINE blocks(x_min,      & ! Global Minimum
                    x_max,      & ! Global Maximum
                    delta_x,    & 
                    block_min,  & ! Individual Block Minimum
                    block_max)    ! Individual Block Maximum

    ! Dummy Variable Declaration
    REAL(KIND=rDef), INTENT(IN) :: x_min, x_max, delta_x
    REAL(KIND=rDef), DIMENSION(6), INTENT(OUT) :: block_min, block_max

    ! Local Variable Declaration
    INTEGER :: num_blocks, i

    num_blocks = 6                              ! Modify later for diff number of blocks?

    ! Calculate Block Minimums and Maximums
    block_min(1) = x_min                        ! Set minimum for first block
    block_max(1) = block_min(1) + x_max/6.0

    WRITE(UNIT=6,FMT=201) 1,block_min(1),block_max(1)   ! Write results to screen


    ! Loop to Generate Grid Blocks and Boundaries
    DO i = 2,num_blocks

      block_min(i) = x_min + block_max(i-1)
      block_max(i) = block_min(i) + x_max/6.0           
      WRITE(UNIT=6,FMT=201) i,block_min(i),block_max(i) ! Write results to screen

    END DO

    
    ! Format Statement(s)
    201 FORMAT("Block",I2,T12,"Lower Bound:",F8.3,T39,"Upper Bound:",F15.3)
    

  END SUBROUTINE blocks

END MODULE block_gen

