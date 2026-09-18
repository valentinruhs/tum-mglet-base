
!====================================================================
!  Module: mph_utils_mod
!
!   Responsibilities:
!   - 
!
!   Author:      Valentin Ruhs
!   e-Mail:      valentin.ruhs@gmx.de
!   Created:     2026-02
!   Last update: 2026-09
!
!====================================================================

MODULE mph_utils_mod

    USE grids_mod, ONLY: nmygrids, mygrids, get_mgdims, iparent, &
        idprocofgrd, iposition, jposition, kposition, ngrid
    USE comms_mod, ONLY: myid
    USE precision_mod, ONLY: intk, realk
    USE mphcore_mod, ONLY: divTol, volTol, vofErr
    USE fields_mod, ONLY: get_fieldptr, set_field
    USE err_mod, ONLY: err_abort

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: init_mph_utils, finish_mph_utils, sel_ind, sel_extent, &
        sel_vel, clp, int2char

CONTAINS

    SUBROUTINE init_mph_utils()

        ! Subroutine arguments
        ! None

        ! Local variables
        CHARACTER(len=*), PARAMETER :: descGrdmask = "uncov. cells"
        INTEGER(intk) :: n, igrid, igridf, ipar
        INTEGER(intk) :: kk, jj, ii, kc0, jc0, ic0
        REAL(realk), POINTER, CONTIGUOUS :: grdMask(:,:,:)

        CALL set_field("GRDMASK", description=descGrdmask, &
            dread=.FALSE., required=.TRUE., dwrite=.FALSE., buffers=.TRUE.)

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_fieldptr(grdMask, "GRDMASK", igrid)
            grdMask = 1.0_realk
        END DO

        DO igridf = 1, ngrid
            ipar = iparent(igridf)

            IF (ipar == 0) CYCLE
            IF (idprocofgrd(ipar) /= myid) CYCLE

            CALL get_fieldptr(grdMask, "GRDMASK", ipar)
            CALL get_mgdims(kk, jj, ii, igridf)

            ic0 = iposition(igridf)
            jc0 = jposition(igridf)
            kc0 = kposition(igridf)

            grdMask(kc0:kc0+(kk-4)/2-1, &
                    jc0:jc0+(jj-4)/2-1, &
                    ic0:ic0+(ii-4)/2-1) = 0.0_realk
        END DO

    END SUBROUTINE init_mph_utils

    !================================================================

    SUBROUTINE finish_mph_utils()

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        CONTINUE

    END SUBROUTINE finish_mph_utils

    !================================================================

    SUBROUTINE sel_ind(lOrq, io, jo, ko)
    !----------------------------------------------------------------
    !   What it does:
    !   Depending on the direction (l) of component (q) the 
    !   indices io, jo and ko are set to 0 or 1. 
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

    SUBROUTINE sel_extent(kk, jj, ii, q, l, dx, dy, dz, ddx, ddy, ddz, dsx, dsy, dsz)
    !----------------------------------------------------------------
    !   What it does:
    !   Depending on the direction (l) of component (q) the 
    !   extents ds(.) are set to d(.) or dd(.). The term ds(.) stands
    !   for spacing in (.)-direction and is a neutral specification
    !   for face-to-face or center-to-center distance.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, q, l
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(out) :: dsx(ii), dsy(jj), dsz(kk)

        ! Local variables
        ! None

        IF ( l < 1 .OR. l > 3 .OR. q < 1 .OR. q > 3 ) THEN
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

    END SUBROUTINE sel_extent

    !================================================================

    SUBROUTINE sel_vel(lOrq, u, v, w, vel)
    !----------------------------------------------------------------
    !   What it does:
    !   Select the for lOrq relevant velocity.
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

    SUBROUTINE comp_vol()
    !----------------------------------------------------------------
    !   What it does:
    !
    !----------------------------------------------------------------

        ! Subroutine arguments

        ! Local variables


    END SUBROUTINE comp_vol

    !================================================================

    SUBROUTINE comp_vol_grd(kk, jj, ii, c, grdMask, ddx, ddy, ddz)
    !----------------------------------------------------------------
    !   What it does:
    !
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: c(kk, jj, ii), grdMask(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)

        ! Local variables
        INTEGER(intk) :: k, j, i

        volFld1 = 0.0_realk
        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    vol = ddx(i)*ddy(j)*ddz(k)
                    volFld1 = volFld1 + c(k,j,i)*vol*grdMask(k,j,i)
                END DO
            END DO
        END DO

    END SUBROUTINE comp_vol_grd

    !================================================================

    ELEMENTAL FUNCTION clp(c) RESULT(cc)
    !----------------------------------------------------------------
    !   What it does:
    !   Clips c to its boundaries [0, 1].
    !
    !   Source:
    !   T. Arrufat et al., “A mass-momentum consistent, 
    !   Volume-of-Fluid method for incompressible flow on staggered 
    !   grids,” Computers & Fluids, vol. 215, p. 104785, Jan. 2021, 
    !   doi: 10.1016/j.compfluid.2020.104785.
    !----------------------------------------------------------------

        REAL(realk), INTENT(in) :: c
        REAL(realk) :: cc
        cc = MAX(MIN(c, 1.0_realk), 0.0_realk)

    END FUNCTION clp

    !================================================================

    FUNCTION int2char(i) RESULT(c)
        INTEGER(intk), INTENT(in) :: i
        CHARACTER(len=1) :: c
        WRITE(c, '(I0)') i
    END FUNCTION int2char

END MODULE mph_utils_mod