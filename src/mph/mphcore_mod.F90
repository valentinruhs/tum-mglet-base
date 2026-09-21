!====================================================================
!   Module: mphcore_mod
!
!   Responsibilities:
!   - reads and provides control and physical parameters
!   - initializes multi-phase fields
!
!   Author:      Valentin Ruhs
!   e-Mail:      valentin.ruhs@gmx.de
!   Created:     2026-02
!   Last update: 2026-09
!
!====================================================================

MODULE mphcore_mod

    USE fort7_mod, ONLY: fort7
    USE config_mod, ONLY: config_t
    USE fields_mod, ONLY: set_field
    USE precision_mod, ONLY: intk, realk
    USE comms_mod, ONLY: myid
    USE err_mod, ONLY: errr
    
    IMPLICIT NONE(type, external)
    PRIVATE

    ! Control parameters
    LOGICAL, PROTECTED :: hasMph
    CHARACTER(len=30), PROTECTED :: mphTst
    INTEGER(intk), PROTECTED :: splPer
    LOGICAL, PROTECTED :: skpAdv, skpDif, skpPre, skpExt
    REAL(realk), PROTECTED :: vofTol, divTol, volTol
    CHARACTER(len=5), PROTECTED :: advScm
    LOGICAL, PROTECTED :: donCen
    LOGICAL, PROTECTED :: vofChk, volChk, divChk

    ! Physical parameters
    REAL(realk), PROTECTED :: rho1, rho2
    REAL(realk), PROTECTED :: gmol1, gmol2
    REAL(realk), PROTECTED :: grav(3)

    ! Error codes
    INTEGER(intk), PARAMETER :: mphInitErr = 124, propsErr = 125, vofErr = 126, plicErr = 127

    PUBLIC :: init_mphcore, finish_mphcore, hasMph, &
        mphTst, skpAdv, skpDif, skpPre, skpExt, splPer, vofTol, divTol, &
        voltol, advScm, donCen, vofChk, volChk, divChk, rho1, rho2, &
        gmol1, gmol2, grav, propsErr, vofErr, plicErr, mphInitErr

CONTAINS

    SUBROUTINE init_mphcore()

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
        CALL mphConf%get_value("/splitPermutation", splPer, 3_intk)
        CALL mphConf%get_value("/vofTolerance", vofTol, 1.0E-12_realk)
        CALL mphConf%get_value("/divTolerance", divTol, 1.0E-8_realk)
        CALL mphConf%get_value("/volTolerance", voltol, 1.0E-8_realk)
        CALL mphConf%get_value("/advectionScheme", advScm, "QUICK")
        CALL mphConf%get_value("/donatingCentered", donCen, .FALSE.)
        CALL mphConf%get_value("/vofCheck", vofChk, .FALSE.)
        CALL mphConf%get_value("/divCheck", divChk, .FALSE.)
        CALL mphConf%get_value("/volCheck", volChk, .FALSE.)

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

        ! Initialize mph-flds
        CALL set_field("UP", description=descVel, istag=1, &
            dread=.FALSE., required=.TRUE., dwrite=.FALSE., buffers=.TRUE.)
        CALL set_field("VP", description=descVel, jstag=1, &
            dread=.FALSE., required=.TRUE., dwrite=.FALSE., buffers=.TRUE.)
        CALL set_field("WP", description=descVel, kstag=1, &
            dread=.FALSE., required=.TRUE., dwrite=.FALSE., buffers=.TRUE.)

        CALL set_field("C", description=descC, &
            dread=.FALSE., required=.TRUE., dwrite=.TRUE., buffers=.TRUE.)
        CALL set_field("CP", description=descCp, &
            dread=.FALSE., required=.TRUE., dwrite=.FALSE., buffers=.TRUE.)
        CALL set_field("D", description=descD, &
            dread=.FALSE., required=.TRUE., dwrite=.FALSE., buffers=.TRUE.)

    END SUBROUTINE init_mphcore

    !================================================================

    SUBROUTINE finish_mphcore()

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        CONTINUE

    END SUBROUTINE  finish_mphcore

END MODULE mphcore_mod