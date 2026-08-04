!====================================================================
!  Module: multiphasecore_mod
!
!  Description:
!     Provides data structures and basic operations for the
!     multiphase model.
!
!  Responsibilities:
!     - Initializes steering values
!     - Initializes multiphase fields
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
    USE precision_mod, ONLY: intk, realk
    USE comms_mod, ONLY: myid
    USE err_mod, ONLY: errr
    
    IMPLICIT NONE(type, external)
    PRIVATE

    ! Control parameters
    LOGICAL, PROTECTED :: has_multiphase, solve_multiphase
    CHARACTER(len=6), PROTECTED :: test_multiphase
    INTEGER(intk), PROTECTED :: permutation_multiphase
    LOGICAL, PROTECTED :: omitAdve, omitDiff, omitExte
    REAL(realk), PROTECTED :: tol
    CHARACTER(len=5), PROTECTED :: fluxLimiter
    LOGICAL, PROTECTED :: checkContinuity, checkSolenoidality, checkBalance

    ! Physical parameters
    REAL(realk), PROTECTED :: rho1, rho2
    REAL(realk), PROTECTED :: gmol1, gmol2
    REAL(realk), PROTECTED :: grav(3)

    PUBLIC :: init_multiphasecore, finish_multiphasecore, &
        has_multiphase, solve_multiphase, test_multiphase, &
        permutation_multiphase, omitAdve, omitDiff, omitExte, &
        tol, fluxLimiter, checkContinuity, checkSolenoidality, &
        checkBalance, rho1, rho2, gmol1, gmol2, grav

CONTAINS

    SUBROUTINE init_multiphasecore()

        ! Subroutine arguments
        ! None

        ! Local variables
        TYPE(config_t) :: multiphaseconf
        INTEGER(intk), PARAMETER :: unitsvp(7) = [0, 1, -1, 0, 0, 0, 0]
        INTEGER(intk), PARAMETER :: unitsvff(7) = [0, 0, 0, 0, 0, 0, 0]
        INTEGER(intk), PARAMETER :: unitsnorm(7) = [0, 0, 0, 0, 0, 0, 0]
        INTEGER(intk), PARAMETER :: unitsalpha(7) = [0, 1, 0, 0, 0, 0, 0]
        CHARACTER(len=*), PARAMETER :: descriptionvp = "Prev. Velocity"
        CHARACTER(len=*), PARAMETER :: descriptionvff = "Volume Fraction Field"
        CHARACTER(len=*), PARAMETER :: descriptionvffp = "Prev. Volume Fraction Field"
        CHARACTER(len=*), PARAMETER :: descriptionnorm = "Norm"
        CHARACTER(len=*), PARAMETER :: descriptionalpha = "Alpha"

        ! Decide wether multiphase is used or not
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

        ! Read steering input
        CALL multiphaseconf%get_value("/solve", solve_multiphase, .TRUE.)
        CALL multiphaseconf%get_value("/test", test_multiphase, 'none')
        CALL multiphaseconf%get_value("/omitAdve", omitAdve, .FALSE.)
        CALL multiphaseconf%get_value("/omitDiff", omitDiff, .FALSE.)
        CALL multiphaseconf%get_value("/omitExte", omitExte, .FALSE.)
        CALL multiphaseconf%get_value("/permutation", permutation_multiphase, 3_intk)
        CALL multiphaseconf%get_value("/tolerance", tol, 1.0E-12_realk)
        CALL multiphaseconf%get_value("/fluxLimiter", fluxLimiter, 'QUICK')
        CALL multiphaseconf%get_value("/checkContinuity", checkContinuity, .FALSE.)
        CALL multiphaseconf%get_value("/checkSolenoidality", checkSolenoidality, .FALSE.)
        CALL multiphaseconf%get_value("/checkBalance", checkBalance, .FALSE.)

        ! Read densities
        CALL multiphaseconf%get_value("/rho1", rho1, 1.0_realk)
        CALL multiphaseconf%get_value("/rho2", rho2, 1.0_realk)
        IF (rho1 <= 0.0 .OR. rho2 <= 0.0) THEN
            WRITE(*, *) "Densities must be positive. rho1 = ", rho1, ", rho2 = ", rho2
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Read viscosities
        CALL multiphaseconf%get_value("/gmol1", gmol1)
        CALL multiphaseconf%get_value("/gmol2", gmol2)
        IF (gmol1 <= 0.0_realk .OR. gmol2 <= 0.0_realk) THEN
            WRITE(*, *) "Viscosities must be positive. gmol1 = ", gmol1, ", gmol2 = ", gmol2
            CALL errr(__FILE__, __LINE__)
        END IF

        ! Read gravity
        CALL multiphaseconf%get_array("/gravity", grav)

        ! Initialize multiphase fields
        CALL set_field("UP", description=descriptionvp, istag=1, units=unitsvp, &
            dread=.FALSE., required=dread, dwrite=.FALSE., buffers=.TRUE.)
        CALL set_field("VP", description=descriptionvp, jstag=1, units=unitsvp, &
            dread=.FALSE., required=dread, dwrite=.FALSE., buffers=.TRUE.)
        CALL set_field("WP", description=descriptionvp, kstag=1, units=unitsvp, &
            dread=.FALSE., required=dread, dwrite=.FALSE., buffers=.TRUE.)
        CALL set_field("VFF", description=descriptionvff , units=unitsvff, &
            dread=dread, required=dread, dwrite=dwrite, buffers=.TRUE.)
        CALL set_field("VFFP", description=descriptionvffp , units=unitsvff, &
            dread=dread, required=dread, dwrite=.FALSE., buffers=.TRUE.)
        CALL set_field("NORMX", description=descriptionnorm , units=unitsnorm, &
            dread=.FALSE., required=dread, dwrite=dwrite, buffers=.TRUE.)
        CALL set_field("NORMY", description=descriptionnorm , units=unitsnorm, &
            dread=.FALSE., required=dread, dwrite=dwrite, buffers=.TRUE.)
        CALL set_field("NORMZ", description=descriptionnorm , units=unitsnorm, &
            dread=.FALSE., required=dread, dwrite=dwrite, buffers=.TRUE.)
        CALL set_field("ALPHA", description=descriptionalpha , units=unitsalpha, &
            dread=.FALSE., required=dread, dwrite=dwrite, buffers=.TRUE.)

    END SUBROUTINE init_multiphasecore

    !================================================================

    SUBROUTINE finish_multiphasecore

        
        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        continue
    END SUBROUTINE  finish_multiphasecore

END MODULE multiphasecore_mod