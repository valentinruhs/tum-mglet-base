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

    USE core_mod
    
    IMPLICIT NONE
    PRIVATE

CONTAINS

    SUBROUTINE init_multiphasecore()

        ! Subroutine arguments
        ! none

        ! Local variables
        TYPE(config_t) :: multiphaseconf

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

        ! CALL flowconf%get_value("/", gmol)




        

    END SUBROUTINE init_multiphasecore

    SUBROUTINE finish_multiphasecore
        continue
    END SUBROUTINE  finish_multiphasecore

END MODULE multiphasecore_mod