!====================================================================
!  Module: mph_mod
!
!   Responsibilities:
!   - Serves as public gateway of the multi-phase solver
!   - Initializes and finishes all multi-phase modules
!   - Advances one RK stage of the multi-phase momentum equation
!     (mph_step): geometric VOF advection and property update in
!     the first stage, viscous, pressure and gravity terms in every
!     stage
!
!   Author:      Valentin Ruhs
!   e-Mail:      valentin.ruhs@gmx.de
!   Created:     2026-02
!
!====================================================================

MODULE mph_mod

    USE precision_mod, ONLY: intk, realk
    USE field_mod, ONLY: field_t

    USE mphcore_mod, ONLY: init_mphcore, finish_mphcore, hasMph
    USE mph_vof_mod, ONLY: init_mph_vof, finish_mph_vof, cpy_flds, &
        adve_operator, diff_operator, pres_operator, exte_operator
    USE mph_plic_mod, ONLY: init_mph_plic, finish_mph_plic
    USE mph_props_mod, ONLY: init_mph_props, finish_mph_props, comp_props
    USE mph_utils_mod, ONLY: init_mph_utils, finish_mph_utils
    USE mph_test_mod, ONLY: init_mph_test, finish_mph_test, frc_vel_fld, frcVelFld
    USE mph_pois_mod, ONLY: init_mph_pois, finish_mph_pois
    USE mph_chk_mod, ONLY: init_mph_chk, finish_mph_chk, final_chk, chk_ifc_par
    USE mph_xfer_mod, ONLY: init_mph_xfer, finish_mph_xfer
    USE mph_slc_mod, ONLY: init_mph_slc, finish_mph_slc

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: init_mph, finish_mph, mph_step

CONTAINS

    SUBROUTINE init_mph()
    !----------------------------------------------------------------
    !   Short: Initialize all multi-phase modules.
    !
    !   init_mphcore is always called and sets hasMph. The other
    !   modules are only initialized if the parameter file contains a
    !   /multiphase section. The order matters: init_mph_test fills c
    !   and reconstructs the interface, init_mph_chk afterwards
    !   stores the initial volume and error norms.
    !----------------------------------------------------------------

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        CALL init_mphcore()

        IF ( hasMph ) THEN
            CALL init_mph_slc()
            CALL init_mph_xfer()
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
    !----------------------------------------------------------------
    !   Short: Run final check and finalize all multi-phase modules.
    !----------------------------------------------------------------

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        IF ( hasMph ) THEN
            CALL final_chk()

            CALL finish_mph_chk()
            CALL finish_mph_test()
            CALL finish_mph_pois()
            CALL finish_mph_utils()
            CALL finish_mph_vof()
            CALL finish_mph_plic()
            CALL finish_mph_props()
            CALL finish_mph_xfer()
            CALL finish_mph_slc()
        END IF

        CALL finish_mphcore()

    END SUBROUTINE finish_mph

    !================================================================

    SUBROUTINE mph_step(uo_f, vo_f, wo_f, dt, irk, itstep, timeph, gradp)
    !----------------------------------------------------------------
    !   Short: Advance one RK stage: advection and properties (first
    !          stage only), viscous, pressure and gravity terms.
    !
    !   Called by timeintegrate_flow instead of the single-phase
    !   right-hand side (tstle4).
    !
    !   The geometric VOF advection is a one-step method. Hence, it
    !   is performed only in the first RK stage over the full time
    !   step dt with the velocities at the beginning of the step (UP,
    !   VP, WP). If frcVelFld is set, the prescribed velocity field
    !   of the test case is imposed at t + dt/2 beforehand.
    !   Afterwards, the PAR face check (chk_ifc_par) is done and the
    !   material properties are updated (comp_props).
    !
    !   The viscous, pressure and gravity terms are added to the
    !   right-hand side uo, vo, wo in every RK stage. gradp is the
    !   prescribed driving pressure gradient.
    !----------------------------------------------------------------

        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: uo_f, vo_f, wo_f
        REAL(realk), INTENT(in) :: dt, timeph
        INTEGER(intk), INTENT(in) :: irk, itstep
        REAL(realk), INTENT(in) :: gradp(3)

        IF ( irk == 1 ) THEN
            IF ( frcVelFld ) CALL frc_vel_fld(timeph + 0.5_realk*dt)
            CALL cpy_flds()
            CALL adve_operator(dt, itstep)
            CALL chk_ifc_par()
            CALL comp_props()
        END IF

        CALL diff_operator(uo_f, vo_f, wo_f)
        CALL pres_operator(uo_f, vo_f, wo_f, gradp)
        CALL exte_operator(uo_f, vo_f, wo_f)

    END SUBROUTINE mph_step

END MODULE mph_mod