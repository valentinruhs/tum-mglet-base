!====================================================================
!  Module: mph_mod
!
!   Responsibilities:
!   - Initializes all multi-phase modules
!   - Finishes all multi-phase modules
!   - Serves as public multi-phase gateway
!
!   Author:      Valentin Ruhs
!   e-Mail:      valentin.ruhs@gmx.de
!   Created:     2026-02
!   Last update: 2026-09
!
!====================================================================

MODULE mph_mod

    USE mphcore_mod, ONLY: init_mphcore, finish_mphcore, hasMph
    USE mph_vof_mod, ONLY: init_mph_vof, finish_mph_vof
    USE mph_plic_mod, ONLY: init_mph_plic, finish_mph_plic
    USE mph_props_mod, ONLY: init_mph_props, finish_mph_props
    USE mph_utils_mod, ONLY: init_mph_utils, finish_mph_utils
    USE mph_test_mod, ONLY: init_mph_test, finish_mph_test

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: init_mph, finish_mph

CONTAINS

    SUBROUTINE init_mph()

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        CALL init_mphcore()
        IF ( .NOT. hasMph ) RETURN

        CALL init_mph_test()
        CALL init_mph_props()
        CALL init_mph_plic()
        CALL init_mph_vof()
        CALL init_mph_utils()

    END SUBROUTINE init_mph

    !================================================================

    SUBROUTINE finish_mph()

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        IF ( .NOT. has_mph ) RETURN

        CALL finish_mph_utils()
        CALL finish_mph_vof()
        CALL finish_mph_plic()
        CALL finish_mph_props()
        CALL finish_mph_test()

    END SUBROUTINE finish_mph

END MODULE mph_mod