!====================================================================
!  Module: mph_utils_mod
!
!   Responsibilities:
!   - Provides helper routines shared by the multi-phase modules:
!     - index offsets, cell spacings and velocity components of a
!       split direction l or a staggered component q
!     - loop bounds of the unknowns of the main and staggered grids
!       depending on the boundary types (CON faces)
!     - clipping of the volume fraction to [0, 1]
!     - integer to character conversion for field names
!
!   Author:      Valentin Ruhs
!   e-Mail:      valentin.ruhs@gmx.de
!   Created:     2026-02
!
!====================================================================

MODULE mph_utils_mod

    USE precision_mod, ONLY: intk, realk
    USE mphcore_mod, ONLY: vofErr, vofTol
    USE err_mod, ONLY: err_abort
    USE grids_mod, ONLY: get_mgdims, get_mgbasb

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: init_mph_utils, finish_mph_utils, sel_ind, sel_ext, &
        sel_vel, clp_tol, clp, int2char, get_lp_mdf, get_lp_bnds

CONTAINS

    SUBROUTINE init_mph_utils()
    !----------------------------------------------------------------
    !   Short: Initialize module (no operation).
    !----------------------------------------------------------------

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        CONTINUE

    END SUBROUTINE init_mph_utils

    !================================================================

    SUBROUTINE finish_mph_utils()
    !----------------------------------------------------------------
    !   Short: Finalize module (no operation).
    !----------------------------------------------------------------

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        CONTINUE

    END SUBROUTINE finish_mph_utils

    !================================================================

    SUBROUTINE sel_ind(lOrq, io, jo, ko)
    !----------------------------------------------------------------
    !   Short: Define spatial index offsets for direction l or
    !          component q.
    !
    !   (io, jo, ko) is the unit vector of direction lOrq = 1, 2, 3.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: lOrq
        INTEGER(intk), INTENT(out) :: io, jo, ko

        ! Local variables
        ! None

        io = 0
        jo = 0
        ko = 0

        SELECT CASE ( lOrq )
        CASE ( 1 )
            io = 1
        CASE ( 2 )
            jo = 1
        CASE ( 3 )
            ko = 1
        CASE DEFAULT
            CALL err_abort(vofErr, "invalid direction lOrq.", __FILE__, __LINE__)
        END SELECT

    END SUBROUTINE sel_ind

    !================================================================

    SUBROUTINE sel_ext(kk, jj, ii, q, l, dx, dy, dz, ddx, ddy, ddz, dsx, dsy, dsz)
    !----------------------------------------------------------------
    !   Short: Select face-to-face or centre-to-centre spacing for q
    !          and l.
    !
    !   ds(.) = dd(.) (face-to-face) in all directions, except for
    !   direction l if q == l: there, the staggered control volume
    !   spans from centre to centre (ds(.) = d(.)). q = 0 denotes the
    !   main grid.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, q, l
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(out) :: dsx(ii), dsy(jj), dsz(kk)

        ! Local variables
        ! None

        IF ( l < 1 .OR. l > 3 .OR. q < 0 .OR. q > 3 ) THEN
            CALL err_abort(vofErr, "invalid direction l or q.", __FILE__, __LINE__)
        END IF

        dsx = ddx ; dsy = ddy ; dsz = ddz

        IF ( q == l ) THEN
            SELECT CASE ( l )
            CASE ( 1 )
                dsx = dx
            CASE ( 2 )
                dsy = dy
            CASE ( 3 )
                dsz = dz
            END SELECT
        END IF

    END SUBROUTINE sel_ext

    !================================================================

    SUBROUTINE sel_vel(lOrq, u, v, w, vel)
    !----------------------------------------------------------------
    !   Short: Point to the velocity component of direction l or q.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: lOrq
        REAL(realk), POINTER, CONTIGUOUS, INTENT(in) :: u(:,:,:), v(:,:,:), w(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS, INTENT(out) :: vel(:,:,:)

        ! Local variables
        ! None

        vel => NULL()

        SELECT CASE ( lOrq )
        CASE ( 1 )
            vel => u
        CASE ( 2 )
            vel => v
        CASE ( 3 )
            vel => w
        CASE DEFAULT
            CALL err_abort(vofErr, "invalid direction lOrq.", __FILE__, __LINE__)
        END SELECT

    END SUBROUTINE sel_vel

    !================================================================

    SUBROUTINE get_lp_mdf(igrid, nfro, nbac, nrgt, nlft, nbot, ntop, &
        nfu, nbu, nrv, nlv, nbw, ntw)
    !----------------------------------------------------------------
    !   Short: Get boundary codes and loop modifiers at CON faces;
    !          abort for OP1.
    !
    !   nbu, nlv, ntw = 1 at CON faces (code 7) on the back, left and
    !   top side: there, the velocity at the boundary face is an
    !   unknown of this grid and included in the loops. nfu, nrv, nbw
    !   are always 0, since the face at the lower boundary belongs to
    !   the neighbouring grid.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid
        INTEGER(intk), INTENT(out) :: nfro, nbac, nrgt, nlft, nbot, ntop
        INTEGER(intk), INTENT(out) :: nfu, nbu, nrv, nlv, nbw, ntw

        ! Local variables
        ! None

        CALL get_mgbasb(nfro, nbac, nrgt, nlft, nbot, ntop, igrid)

        IF ( ANY([nfro, nbac, nrgt, nlft, nbot, ntop] == 3) ) THEN
            CALL err_abort(vofErr, "OP1 not supported in multiphase solver.", __FILE__, __LINE__)
        END IF

        nfu = 0
        nbu = 0
        nrv = 0
        nlv = 0
        nbw = 0
        ntw = 0

        ! CON = 7
        IF (nbac == 7) nbu = 1
        IF (nlft == 7) nlv = 1
        IF (ntop == 7) ntw = 1

    END SUBROUTINE get_lp_mdf

    !================================================================

    SUBROUTINE get_lp_bnds(igrid, q, ista, iend, jsta, jend, ksta, kend)
    !----------------------------------------------------------------
    !   Short: Get loop bounds of the unknowns of component q.
    !
    !   q = 0: pressure cell, q = 1, 2, 3: u-, v-, w-cell. Only
    !   the bounds in q-direction depend on the boundary type
    !   (see get_lp_mdf).
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, q
        INTEGER(intk), INTENT(out) :: ista, iend, jsta, jend, ksta, kend

        ! Local variables
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: nfro, nbac, nrgt, nlft, nbot, ntop
        INTEGER(intk) :: nfu, nbu, nrv, nlv, nbw, ntw

        CALL get_mgdims(kk, jj, ii, igrid)
        CALL get_lp_mdf(igrid, nfro, nbac, nrgt, nlft, nbot, ntop, &
            nfu, nbu, nrv, nlv, nbw, ntw)

        ista = 3 ; iend = ii-2
        jsta = 3 ; jend = jj-2
        ksta = 3 ; kend = kk-2

        SELECT CASE ( q )
        CASE ( 0 )
            CONTINUE
        CASE ( 1 )
            ista = 3-nfu ; iend = ii-3+nbu
        CASE ( 2 )
            jsta = 3-nrv ; jend = jj-3+nlv
        CASE ( 3 )
            ksta = 3-nbw ; kend = kk-3+ntw
        CASE DEFAULT
            CALL err_abort(vofErr, "invalid component q.", __FILE__, __LINE__)
        END SELECT

    END SUBROUTINE get_lp_bnds

    !================================================================

    ELEMENTAL FUNCTION clp_tol(c) RESULT(cc)
    !----------------------------------------------------------------
    !   Short: Volume fraction clipped to [0, 1] and snapped to 0 or
    !          1 within vofTol.
    !
    !   Values that comp_isIfc_grd does not treat as interface
    !   (c <= vofTol, c >= 1-vofTol) are set to 0 or 1.
    !----------------------------------------------------------------

        REAL(realk), INTENT(in) :: c
        REAL(realk) :: cc
        cc = MAX(MIN(c, 1.0_realk), 0.0_realk)
        IF ( cc <= vofTol ) cc = 0.0_realk
        IF ( cc >= 1.0_realk - vofTol ) cc = 1.0_realk

    END FUNCTION clp_tol

    !================================================================

    ELEMENTAL FUNCTION clp(c) RESULT(cc)
    !----------------------------------------------------------------
    !   Short: Volume fraction clipped to [0, 1].
    !----------------------------------------------------------------

        REAL(realk), INTENT(in) :: c
        REAL(realk) :: cc
        cc = MAX(MIN(c, 1.0_realk), 0.0_realk)

    END FUNCTION clp

    !================================================================

    FUNCTION int2char(i) RESULT(c)
    !----------------------------------------------------------------
    !   Short: Integer converted to character.
    !
    !   Single digit only (len = 1).
    !----------------------------------------------------------------
        INTEGER(intk), INTENT(in) :: i
        CHARACTER(len=1) :: c
        WRITE(c, '(I0)') i
    END FUNCTION int2char

END MODULE mph_utils_mod