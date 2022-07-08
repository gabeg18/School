MODULE block_data
  USE, INTRINSIC :: ISO_FORTRAN_ENV
  IMPLICIT NONE

  INTEGER, PARAMETER :: rDef = REAL64

  !
  ! This module defines a derived data type for each block
  ! 

  TYPE grid_block

    REAL(KIND=rDef) :: lower_boundary, &  ! Local lower boundary
                       upper_boundary     ! Local upper boundary

    INTEGER :: block_number, &            ! Respective block number (i.e. 1,2,3, etc.)
               num_grid_points            ! Number of grid points per block

    LOGICAL :: boundary_condition         ! Boundary condition of block

  END TYPE grid_block

END MODULE block_data
