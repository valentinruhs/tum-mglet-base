!====================================================================
!  Module: mphcore_mod
!
!   Responsibilities:
!   - Reads the /multiphase section of the parameter file and sets
!     hasMph
!   - Provides the control parameters (test case, skip switches,
!     tolerances, advection scheme, ...) and the physical parameters
!     (densities, viscosities, gravity)
!   - Defines the phase convention: index 1 refers to the phase with
!     c = 1, index 2 to the phase with c = 0
!   - Provides the error codes of the multi-phase modules
!   - Registers the basic multi-phase fields: C, CP, D, UP, VP, WP
!
!   Author:      Valentin Ruhs
!   e-Mail:      valentin.ruhs@gmx.de
!   Created:     2026-02
!
!====================================================================

MODULE mphcore_mod

    USE fort7_mod, ONLY: fort7
    USE config_mod, ONLY: config_t
    USE fields_mod, ONLY: set_field
    USE precision_mod, ONLY: intk, realk
    USE comms_mod, ONLY: myid
    USE err_mod, ONLY: errr
    USE fort7_mod, ONLY: dcont

    IMPLICIT NONE(type, external)
    PRIVATE

    ! Control parameters
    LOGICAL, PROTECTED :: hasMph
    CHARACTER(len=30), PROTECTED :: mphTst
    LOGICAL, PROTECTED :: skpAdv, skpDif, skpPre, skpExt
    REAL(realk), PROTECTED :: vofTol, volTol
    CHARACTER(len=5), PROTECTED :: advScm
    LOGICAL, PROTECTED :: donCen, denStg
    INTEGER(intk), PROTECTED :: parIfcBand
    LOGICAL, PROTECTED :: finChk

    ! Physical parameters
    REAL(realk), PROTECTED :: rho1, rho2
    REAL(realk), PROTECTED :: gmol1, gmol2
    REAL(realk), PROTECTED :: grav(3)

    ! Error codes
    INTEGER(intk), PARAMETER :: mphInitErr = 124, propsErr = 125, vofErr = 126, plicErr = 127

    PUBLIC :: init_mphcore, finish_mphcore, hasMph, &
        mphTst, skpAdv, skpDif, skpPre, skpExt, vofTol, &
        volTol, advScm, donCen, denStg, parIfcBand, finChk, rho1, rho2, &
        gmol1, gmol2, grav, propsErr, vofErr, plicErr, mphInitErr

CONTAINS

    SUBROUTINE init_mphcore()
    !----------------------------------------------------------------
    !   Short: Read control and physical parameters; register c,
    !          density and previous-step fields.
    !
    !   If the parameter file has no /multiphase section, hasMph is
    !   set to .FALSE. and the routine returns. Otherwise, the
    !   control parameters (test case, skip switches, tolerances,
    !   advection scheme, ...) and the physical parameters are read:
    !   densities rho1/rho2, dynamic viscosities gmol1/gmol2 and
    !   gravity. Index 1 refers to the phase with c = 1.
    !
    !   Registered fields: UP, VP, WP (velocities of the previous
    !   time step), C (volume fraction, read on restart and written
    !   to the output), CP (volume fraction of the previous time
    !   step) and D (cell density).
    !----------------------------------------------------------------

        ! Subroutine arguments
        ! None

        ! Local variables
        TYPE(config_t) :: mphConf
        CHARACTER(len=*), PARAMETER :: descVel = "prev. vel. fld."
        CHARACTER(len=*), PARAMETER :: descD = "density fld."
        CHARACTER(len=*), PARAMETER :: descC = "vol. frac. fld."
        CHARACTER(len=*), PARAMETER :: descCp = "prev. vol. frac. fld."

        hasMph = .FALSE.
        IF (.NOT. fort7%exists("/multiphase")) THEN
            IF (myid == 0) THEN
                WRITE(*, '("NO MULTIPHASE FLOW")')
                WRITE(*, '()')
            END IF
            RETURN
        END IF
        hasMph = .TRUE.

        ! Initialize mphConf
        CALL fort7%get(mphConf, "/multiphase")

        ! Read steering input
        CALL mphConf%get_value("/test", mphTst, "None")
        CALL mphConf%get_value("/skipAdvection", skpAdv, .FALSE.)
        CALL mphConf%get_value("/skipDiffusion", skpDif, .FALSE.)
        CALL mphConf%get_value("/skipPressure", skpPre, .FALSE.)
        CALL mphConf%get_value("/skipExternal", skpExt, .FALSE.)
        CALL mphConf%get_value("/vofTolerance", vofTol, 1.0E-12_realk)
        CALL mphConf%get_value("/volTolerance", volTol, 1.0E-8_realk)
        CALL mphConf%get_value("/advectionScheme", advScm, "QUICK")
        CALL mphConf%get_value("/donatingCentered", donCen, .FALSE.)
        CALL mphConf%get_value("/densityFaceStag", denStg, .FALSE.)
        CALL mphConf%get_value("/parentInterfaceBand", parIfcBand, 4)
        CALL mphConf%get_value("/finalCheck", finChk, .FALSE.)

        ! Read densities
        CALL mphConf%get_value("/rho1", rho1, 1.0_realk)
        CALL mphConf%get_value("/rho2", rho2, 1.0_realk)
        IF (rho1 <= 0.0_realk .OR. rho2 <= 0.0_realk) THEN
            WRITE(*, *) "Densities must be positive. rho1 = ", rho1, ", rho2 = ", rho2
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Read viscosities
        CALL mphConf%get_value("/gmol1", gmol1, 1.0E-3_realk)
        CALL mphConf%get_value("/gmol2", gmol2, 1.0E-3_realk)
        IF (gmol1 <= 0.0_realk .OR. gmol2 <= 0.0_realk) THEN
            WRITE(*, *) "Viscosities must be positive. gmol1 = ", gmol1, ", gmol2 = ", gmol2
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Read gravity
        CALL mphConf%get_array("/gravity", grav)

        ! Initialize multi-phase fields
        CALL set_field("UP", description=descVel, istag=1, &
            dread=.FALSE., required=.TRUE., dwrite=.FALSE., buffers=.TRUE.)
        CALL set_field("VP", description=descVel, jstag=1, &
            dread=.FALSE., required=.TRUE., dwrite=.FALSE., buffers=.TRUE.)
        CALL set_field("WP", description=descVel, kstag=1, &
            dread=.FALSE., required=.TRUE., dwrite=.FALSE., buffers=.TRUE.)

        CALL set_field("C", description=descC, &
            dread=dcont, required=.TRUE., dwrite=.TRUE., buffers=.TRUE.)
        CALL set_field("CP", description=descCp, &
            dread=.FALSE., required=.TRUE., dwrite=.FALSE., buffers=.TRUE.)
        CALL set_field("D", description=descD, &
            dread=.FALSE., required=.TRUE., dwrite=.FALSE., buffers=.TRUE.)

    END SUBROUTINE init_mphcore

    !================================================================

    SUBROUTINE finish_mphcore()
    !----------------------------------------------------------------
    !   Short: Finalize module (no operation).
    !----------------------------------------------------------------

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        CONTINUE

    END SUBROUTINE  finish_mphcore

END MODULE mphcore_mod