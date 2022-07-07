MODULE block_select
  USE :: ISO_FORTRAN_ENV
  IMPLICIT NONE
  ! This module defines the blocks for use in a 
  ! multi-block solver of the Heat Equation

  ! Beginning with putting the blocks in numerical
  ! order

  ! Variable Declaration
  INTEGER, PARAMETER :: rDef       = REAL64, &
                        num_blocks = 6,      &
                        iMin       = 1,      &  ! Global Minimum
                        iMax       = 101        ! Global Maximum

  REAL(KIND=rDef), PARAMETER :: x_min   =   0.0_rDef,  &
                                x_max   = 100.0_rDef,  &
                                delta_x =   1.0_rDef

  CONTAINS
    SUBROUTINE block_generator(x_min,     &
                               x_max,     &
                               delta_x,   &
                               block_min, &
                               block_max)

    ! This subroutine generates the grid for each block

    ! Dummy Variable Declaration

    REAL(KIND=rDef),INTENT(IN) :: x_min,  &
                                  x_max,  &
                                  delta_x

    REAL(KIND=rDef),INTENT(OUT) :: block_min,  & 
                                   block_max,  &

    ! Local Variable Declaration

    INTEGER :: i

                        
    ! Determining iMax from provided values
    !iMax = (x_max - x_min)/delta_x + 1
    !PRINT *,iMax

    ! Initialize Entire Space Domain
    DO i = 1, iMax

    ! Defining Block Grid
    DO i = 1, num_blocks
      

    END DO

    END SUBROUTINE block_generator

END MODULE block_select

PROGRAM test
  USE block_select
  IMPLICIT NONE

END PROGRAM test
