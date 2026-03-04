!====================================================================
!  Module: multiphase_io_mod
!
!  Responsibilities:
!     - Reads initial color-function field c
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

        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: c
        
        ! Local variables
        INTEGER(intk) :: Nx, Ny, Nz
        INTEGER(intk) :: i, unit

        Nx = 54
        Ny = 54
        Nz = 54

        OPEN(newunit=unit,file="c.csv",status="old",action="read")
        DO i = 1, Nx*Ny*Nz
            READ(unit,*) c%arr(i)
        END DO
        CLOSE(unit)

    END SUBROUTINE read_color_function

END MODULE multiphase_io_mod