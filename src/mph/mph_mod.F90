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

    USE precision_mod, ONLY: intk, realk
    USE field_mod, ONLY: field_t

    USE mphcore_mod, ONLY: init_mphcore, finish_mphcore, hasMph
    USE mph_vof_mod, ONLY: init_mph_vof, finish_mph_vof, cpy_flds, &
        adve_operator, diff_operator, pres_operator, exte_operator
    USE mph_plic_mod, ONLY: init_mph_plic, finish_mph_plic
    USE mph_props_mod, ONLY: init_mph_props, finish_mph_props, &
        comp_props
    USE mph_utils_mod, ONLY: init_mph_utils, finish_mph_utils
    USE mph_test_mod, ONLY: init_mph_test, finish_mph_test
    USE mph_pois_mod, ONLY: init_mph_pois, finish_mph_pois
    USE mph_chk_mod, ONLY: init_mph_chk, finish_mph_chk, final_vol_chk

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: init_mph, finish_mph, mph_step

CONTAINS

    SUBROUTINE init_mph()

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        CALL init_mphcore()

        IF ( hasMph ) THEN
            CALL init_mph_props()
            CALL init_mph_plic()
            CALL init_mph_vof()
            CALL init_mph_utils()
            CALL init_mph_pois()
            CALL init_mph_test()
            CALL init_mph_chk()
        END IF

    END SUBROUTINE init_mph

    !================================================================

    SUBROUTINE finish_mph()

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        CALL final_vol_chk()

        IF ( hasMph ) THEN
            CALL finish_mph_chk()
            CALL finish_mph_test()
            CALL finish_mph_pois()
            CALL finish_mph_utils()
            CALL finish_mph_vof()
            CALL finish_mph_plic()
            CALL finish_mph_props()
        END IF

        CALL finish_mphcore()

    END SUBROUTINE finish_mph

    !================================================================

    SUBROUTINE mph_step(uo_f, vo_f, wo_f, dt, itstep)
    !----------------------------------------------------------------
    !   What it does:
    !    
    !----------------------------------------------------------------

        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: uo_f, vo_f, wo_f
        REAL(realk), INTENT(in) :: dt
        INTEGER(intk), INTENT(in) :: itstep

        ! Local variables
        ! None

        CALL cpy_flds()
        CALL adve_operator(dt, itstep)
        CALL comp_props()
        CALL diff_operator(uo_f, vo_f, wo_f)
        CALL pres_operator(uo_f, vo_f, wo_f)
        CALL exte_operator(uo_f, vo_f, wo_f)

    END SUBROUTINE mph_step

END MODULE mph_mod