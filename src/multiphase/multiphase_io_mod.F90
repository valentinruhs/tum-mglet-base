!====================================================================
!  Module: multiphase_io_mod
!
!  Responsibilities:
!     - Reads initial volume fraction field vff
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

    PUBLIC :: init_multiphase_io, finish_multiphase_io, read_vff

CONTAINS

    SUBROUTINE init_multiphase_io()

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        continue
    END SUBROUTINE init_multiphase_io

    !================================================================

    SUBROUTINE finish_multiphase_io()

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        continue
    END SUBROUTINE finish_multiphase_io

    !================================================================

    SUBROUTINE read_vff(vff)
    !----------------------------------------------------------------
    !   What it does:
    !   This subroutine reads the initial values of the volume
    !   fraction field from a csv-file generated in column-major 
    !   order. The grid size needs to be specified via Nx, Ny and Nz.
    !----------------------------------------------------------------

        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: vff
        
        ! Local variables
        INTEGER(intk) :: Nx, Ny, Nz
        INTEGER(intk) :: i, unit

        Nx = 20
        Ny = 20
        Nz = 20

        OPEN(newunit=unit,file="vff.csv",status="old",action="read")
        DO i = 1, Nx*Ny*Nz
            READ(unit,*) vff%arr(i)
        END DO
        CLOSE(unit)

    END SUBROUTINE read_vff

END MODULE multiphase_io_mod