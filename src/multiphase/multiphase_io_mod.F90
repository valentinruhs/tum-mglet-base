!====================================================================
!  Module: multiphase_io_mod
!
!  Responsibilities:
!     - Reads initial Color-Function field c
!
!  Author:      Valentin Ruhs
!  Created:     2026-02
!  Last update: 2026-02
!
!====================================================================

MODULE multiphase_io_mod

    USE precision_mod, ONLY: intk
    USE field_mod, ONLY: field_t

    IMPLICIT NONE
    PRIVATE 

    PUBLIC :: read_color_function

CONTAINS

    SUBROUTINE init_multiphase_io()

        continue

    END SUBROUTINE init_multiphase_io

    !================================================================

    SUBROUTINE finish_multiphase_io()

        continue

    END SUBROUTINE finish_multiphase_io

    !================================================================

    SUBROUTINE read_color_function(c)
    !----------------------------------------------------------------
    !   What it does:
    !   This subroutine reads the initial values of the Color-
    !   Function field from a csv-file generated in column-major 
    !   order. The grid size needs to be specified via Nx, Ny and Nz.
    !----------------------------------------------------------------

        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: c
        
        ! Local variables
        INTEGER(intk) :: Nx, Ny, Nz
        INTEGER(intk) :: i, unit

        Nx = 20
        Ny = 20
        Nz = 20

        OPEN(newunit=unit,file="c.csv",status="old",action="read")
        DO i = 1, Nx*Ny*Nz
            READ(unit,*) c%arr(i)
        END DO
        CLOSE(unit)

    END SUBROUTINE read_color_function

END MODULE multiphase_io_mod