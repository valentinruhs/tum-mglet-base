!====================================================================
!  Module: multiphasecore_mod
!
!  Description:
!     Provides data structures and basic operations for the
!     multiphase model.
!
!  Responsibilities:
!     - Stores volume fraction field C
!     - Computes mixture density and viscosity
!     - Provides interface geometry utilities
!
!  Author:      Valentin Ruhs
!  Created:     2026-02
!  Last update: 2026-02
!
!====================================================================

MODULE multiphasecore_mod

    USE fort7_mod, ONLY: fort7, dread, dwrite
    USE config_mod, ONLY: config_t
    USE fields_mod, ONLY: set_field
    USE precision_mod, ONLY: intk
    USE comms_mod, ONLY: myid
    
    IMPLICIT NONE(type, external)
    PRIVATE

    ! Control parameters
    LOGICAL, PROTECTED :: has_multiphase, solve_multiphase

    ! Physical parameters

    PUBLIC :: init_multiphasecore, finish_multiphasecore, has_multiphase, solve_multiphase

CONTAINS

    SUBROUTINE init_multiphasecore()

        ! Subroutine arguments
        ! None

        ! Local variables
        TYPE(config_t) :: multiphaseconf
        INTEGER(intk), PARAMETER :: units_c(7) = [0, 0, 0, 0, 0, 0, 0]
        CHARACTER(len=*), PARAMETER :: description_c = "Color-function field c is 0 for fluid 1 and 1 for fluid 2"

        ! Read configuration values for multiphase flow
        has_multiphase = .FALSE.
        IF (.NOT. fort7%exists("/multiphase")) THEN
            IF (myid == 0) THEN
                WRITE(*, '("NO MULTIPHASE FLOW")')
                WRITE(*, '()')
            END IF
            RETURN
        END IF
        has_multiphase = .TRUE.

        ! Initialize multiphaseconf
        CALL fort7%get(multiphaseconf, "/multiphase")
        CALL multiphaseconf%get_value("/solve", solve_multiphase, .TRUE.)

        CALL set_field("C", description=description_c , units=units_c, &
            dread=dread, required=dread, dwrite=dwrite, buffers=.TRUE.)

    END SUBROUTINE init_multiphasecore

    !================================================================

    SUBROUTINE finish_multiphasecore
        continue
    END SUBROUTINE  finish_multiphasecore

END MODULE multiphasecore_mod