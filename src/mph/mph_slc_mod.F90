!====================================================================
!  Module: mph_slc_mod
!
!   Responsibilities:
!   - Registers and fills the reduced snapshot fields of the open
!     channel flow for video
!   - Values are cell-centred and set in the interior only
!
!   Fixed grid: finest level 128 x 96 x 128 uniform cells, origin
!   x = y = z = 0 (Lx = 2 pi, Ly = 1.2, Lz = pi).
!
!====================================================================

MODULE mph_slc_mod

    USE precision_mod, ONLY: intk, realk
    USE fields_mod, ONLY: set_field, get_fieldptr
    USE grids_mod, ONLY: get_mgdims, get_bbox, get_level, get_jj, &
        minlevel, maxlevel, mygrids, nmygrids

    IMPLICIT NONE(type, external)
    PRIVATE

    ! Global cell indices of the planes on the finest level
    INTEGER(intk), PARAMETER :: jWal = 9     ! y = 0.106 (y+ = 16)
    INTEGER(intk), PARAMETER :: jSur = 77    ! y = 0.956
    INTEGER(intk), PARAMETER :: iMid = 65    ! x = Lx/2 + dx/2
    INTEGER(intk), PARAMETER :: kMid = 65    ! z = Lz/2 + dz/2

    PUBLIC :: init_mph_slc, finish_mph_slc, mph_fields_to_write

CONTAINS

    SUBROUTINE init_mph_slc()
    !----------------------------------------------------------------
    !   Short: Register the slice fields on the finest level.
    !----------------------------------------------------------------

        ! Local variables
        INTEGER(intk) :: l
        LOGICAL :: fin(minlevel:maxlevel)

        fin = [(l == maxlevel, l = minlevel, maxlevel)]

        CALL set_field("UXZWAL", ndim=2, normal_direction=2, active_level=fin)
        CALL set_field("VXZSUR", ndim=2, normal_direction=2, active_level=fin)
        CALL set_field("ETAXZ", ndim=2, normal_direction=2, active_level=fin)
        CALL set_field("UXYMID", ndim=2, normal_direction=3, active_level=fin)
        CALL set_field("CXYMID", ndim=2, normal_direction=3, active_level=fin)
        CALL set_field("UYZMID", ndim=2, normal_direction=1, active_level=fin)
        CALL set_field("VYZMID", ndim=2, normal_direction=1, active_level=fin)
        CALL set_field("WYZMID", ndim=2, normal_direction=1, active_level=fin)
        CALL set_field("CYZMID", ndim=2, normal_direction=1, active_level=fin)
        CALL set_field("UPROF", ndim=1, get_len=get_jj, active_level=fin)

    END SUBROUTINE init_mph_slc

    !================================================================

    SUBROUTINE finish_mph_slc()
    !----------------------------------------------------------------
    !   Short: Finalize module (no operation).
    !----------------------------------------------------------------

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        CONTINUE

    END SUBROUTINE finish_mph_slc

    !================================================================

    SUBROUTINE mph_fields_to_write()
    !----------------------------------------------------------------
    !   Short: Fill the slice fields from U, V, W and C.
    !
    !   Call directly before the snapshot is written (ghost layers of
    !   U, V, W valid after mgpoisl).
    !----------------------------------------------------------------

        ! Local variables
        INTEGER(intk) :: n, igrid, ilevel, kk, jj, ii, i, jW, jS, k
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz

        DO n = 1, nmygrids
            igrid = mygrids(n)
            CALL get_level(ilevel, igrid)
            IF ( ilevel /= maxlevel ) CYCLE

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

            ! Local indices of the planes (0: not in this grid)
            i = loc_ind(iMid, minx, maxx, ii)
            jW = loc_ind(jWal, miny, maxy, jj)
            jS = loc_ind(jSur, miny, maxy, jj)
            k = loc_ind(kMid, minz, maxz, kk)

            CALL fill_slc(igrid, 2, jW, "U", 1, "UXZWAL")
            CALL fill_slc(igrid, 2, jS, "V", 2, "VXZSUR")
            CALL fill_slc(igrid, 3, k, "U", 1, "UXYMID")
            CALL fill_slc(igrid, 3, k, "C", 0, "CXYMID")
            CALL fill_slc(igrid, 1, i, "U", 1, "UYZMID")
            CALL fill_slc(igrid, 1, i, "V", 2, "VYZMID")
            CALL fill_slc(igrid, 1, i, "W", 3, "WYZMID")
            CALL fill_slc(igrid, 1, i, "C", 0, "CYZMID")
            CALL fill_sum(igrid)
        END DO

    END SUBROUTINE mph_fields_to_write

    !================================================================

    PURE FUNCTION loc_ind(g, lo, hi, n) RESULT(p)
    !----------------------------------------------------------------
    !   Short: Local index of global cell g in a grid from lo to hi
    !          with n cells incl. ghosts; 0 if not in the interior.
    !----------------------------------------------------------------

        ! Function arguments
        INTEGER(intk), INTENT(in) :: g, n
        REAL(realk), INTENT(in) :: lo, hi
        INTEGER(intk) :: p

        p = g - NINT(lo*(n-4)/(hi-lo)) + 2
        IF ( p < 3 .OR. p > n-2 ) p = 0

    END FUNCTION loc_ind

    !================================================================

    SUBROUTINE fill_slc(igrid, l, p, fldNam, q, slcNam)
    !----------------------------------------------------------------
    !   Short: Fill slice slcNam with fldNam in plane p normal to l.
    !
    !   q = 0: fldNam cell-centred, q = 1, 2, 3: staggered in q.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, l, p, q
        CHARACTER(len=*), INTENT(in) :: fldNam, slcNam

        ! Local variables
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: f(:,:,:), slc(:,:)

        CALL get_mgdims(kk, jj, ii, igrid)
        CALL get_fieldptr(f, fldNam, igrid)
        CALL get_fieldptr(slc, slcNam, igrid)

        CALL fill_slc_grd(kk, jj, ii, l, p, q, f, slc)

    END SUBROUTINE fill_slc

    !================================================================

    SUBROUTINE fill_slc_grd(kk, jj, ii, l, p, q, f, slc)
    !----------------------------------------------------------------
    !   Short: Grid kernel of fill_slc().
    !
    !   Cell-centred value (f(k,j,i) + f(k-kq,j-jq,i-iq))/2 in plane
    !   p; slc(k,j), slc(k,i), slc(j,i) for l = 1, 2, 3.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, l, p, q
        REAL(realk), INTENT(in) :: f(kk, jj, ii)
        REAL(realk), INTENT(out) :: slc(:,:)

        ! Local variables
        INTEGER(intk) :: k, j, i, iq, jq, kq

        slc = 0.0_realk
        IF ( p == 0 ) RETURN

        iq = MERGE(1, 0, q == 1)
        jq = MERGE(1, 0, q == 2)
        kq = MERGE(1, 0, q == 3)

        SELECT CASE ( l )
        CASE ( 1 )
            DO j = 3, jj-2
                DO k = 3, kk-2
                    slc(k, j) = 0.5_realk*(f(k, j, p) + f(k-kq, j-jq, p-iq))
                END DO
            END DO
        CASE ( 2 )
            DO i = 3, ii-2
                DO k = 3, kk-2
                    slc(k, i) = 0.5_realk*(f(k, p, i) + f(k-kq, p-jq, i-iq))
                END DO
            END DO
        CASE ( 3 )
            DO i = 3, ii-2
                DO j = 3, jj-2
                    slc(j, i) = 0.5_realk*(f(p, j, i) + f(p-kq, j-jq, i-iq))
                END DO
            END DO
        END SELECT

    END SUBROUTINE fill_slc_grd

    !================================================================

    SUBROUTINE fill_sum(igrid)
    !----------------------------------------------------------------
    !   Short: Fill ETAXZ (column sums of c) and UPROF (plane sums
    !          of u) of grid igrid.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid

        ! Local variables
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: c(:,:,:), u(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:), eta(:,:), prf(:)

        CALL get_mgdims(kk, jj, ii, igrid)
        CALL get_fieldptr(c, "C", igrid)
        CALL get_fieldptr(u, "U", igrid)
        CALL get_fieldptr(ddx, "DDX", igrid)
        CALL get_fieldptr(ddy, "DDY", igrid)
        CALL get_fieldptr(ddz, "DDZ", igrid)
        CALL get_fieldptr(eta, "ETAXZ", igrid)
        CALL get_fieldptr(prf, "UPROF", igrid)

        CALL fill_sum_grd(kk, jj, ii, c, u, ddx, ddy, ddz, eta, prf)

    END SUBROUTINE fill_sum

    !================================================================

    SUBROUTINE fill_sum_grd(kk, jj, ii, c, u, ddx, ddy, ddz, eta, prf)
    !----------------------------------------------------------------
    !   Short: Grid kernel of fill_sum().
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: c(kk, jj, ii), u(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(out) :: eta(kk, ii), prf(jj)

        ! Local variables
        INTEGER(intk) :: k, j, i

        eta = 0.0_realk
        prf = 0.0_realk

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    eta(k, i) = eta(k, i) + c(k, j, i)*ddy(j)
                    prf(j) = prf(j) + 0.5_realk*(u(k, j, i) + u(k, j, i-1))*ddx(i)*ddz(k)
                END DO
            END DO
        END DO

    END SUBROUTINE fill_sum_grd

END MODULE mph_slc_mod