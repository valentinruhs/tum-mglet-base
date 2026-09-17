!====================================================================
!  Module: mph_test_mod
!
!   Responsibilities:
!   - 
!
!   Author:      Valentin Ruhs
!   e-Mail:      valentin.ruhs@gmx.de
!   Created:     2026-09
!   Last update: 2026-09
!
!====================================================================

MODULE mph_test_mod

    USE mphcore_mod, ONLY: mphTst, mphInitErr
    USE err_mod, ONLY: err_abort
    USE precision_mod, ONLY: realk, intk, pi

    IMPLICIT NONE(type, external)
    PRIVATE

    INTEGER(intk), PARAMETER :: shpCircle = 1, shpZalesak = 2, &
        hpEllipse = 3, shpPlane = 4

    PUBLIC :: init_mph_test, finish_mph_test

    ABSTRACT INTERFACE
        PURE FUNCTION dist_func_interface(x, y, z) RESULT(phi)
            REAL(realk), INTENT(in) :: x, y, z
            REAL(realk) :: phi
        END FUNCTION dist_func_interface
    END INTERFACE

CONTAINS

    SUBROUTINE init_mph_test()

        ! Subroutine arguments
        ! None

        ! Local variables
        INTEGER(intk) :: shp
        REAL(realk) :: xc, yc, zc, ra, slotDx, slotDy, rb, theta, lvl

        SELECT CASE ( mphTst )
        CASE ( "Cylinder Translation" )
            shp = shpCircle
            xc = 0.5_realk
            yc = 0.5_realk
            ra = 0.15_realk
        CASE ( "Zalesak Disk" )
            shp = shpZalesak
            xc = 0.5_realk
            yc = 0.5_realk
            ra = 0.15_realk
            slotDx = 0.05_realk
            slotDy = 0.25_realk
        CASE ( "Rider-Kothe Vortex" )
            shp = shpCircle
            xc = 0.5_realk
            yc = 0.75_realk
            ra = 0.15_realk
        CASE ( "Uniform Cylinder Advection" )
            shp = shpCircle
            xc = 0.2_realk
            yc = 0.2_realk
            ra = 0.1_realk
        CASE ( "Abrupt Cylinder Advection" )
            shp = shpCircle
            xc = 0.2_realk
            yc = 0.2_realk
            ra = 0.1_realk
        CASE ( "Ellipse Reconstruction" )
            shp = shpEllipse
            xc = 0.5_realk
            yc = 0.5_realk
            ra = 0.3464_realk
            rb = 0.1414_realk
            theta = pi
        CASE ( "Open Channel Flow" )
            shp = shpPlane
            lvl = 0.0_realk
        CASE DEFAULT
            CALL err_abort(mphInitErr, "unknown test case.", __FILE__, __LINE__)
        END SELECT

        CALL fill_c_dom(dist_func)
        CALL fill_c_bou(shp)

    CONTAINS

        PURE FUNCTION dist_func(x, y, z) RESULT(phi)
            REAL(realk), INTENT(in) :: x, y, z
            REAL(realk) :: phi

            SELECT CASE ( shp )
            CASE ( shpCircle )
                phi = dist_func_circle(x, y, xc, yc, ra)
            CASE ( shpZalesak )
                phi = MAX(dist_func_circle(x, y, xc, yc, ra), &
                         -dist_func_box(x, y, xc, yc - ra + 0.5_realk*slotH, 0.5_realk*slotW, 0.5_realk*slotH))
            CASE ( shpEllipse )
                phi = dist_func_ellipse(x, y, xc, yc, ra, rb, theta)
            CASE ( shpPlane )
                phi = y - lvl
            CASE DEFAULT
                phi = 1.0_realk
            END SELECT
        END FUNCTION dist_func

    END SUBROUTINE init_mph_test

    !================================================================

    SUBROUTINE finish_mph_test()

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        CONTINUE

    END SUBROUTINE finish_mph_test

    !================================================================

    SUBROUTINE fill_c_dom(dist_func)
    !----------------------------------------------------------------
    !   What it does:
    !   Fill the volume fraction field c depending on the selected
    !   distance function. Only the inside of the domain is filled.
    !----------------------------------------------------------------

        ! Subroutine arguments
        PROCEDURE(dist_func_interface) :: dist_func

        ! Local variables
        INTEGER(intk) :: n, igrid, kk, jj, ii
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz
        REAL(realk), POINTER, CONTIGUOUS :: c(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        INTEGER(intk) :: i, j, k
        REAL(realk) :: xm, x, ym, y, zm, z, halfDiag, dist
        INTEGER(intk) :: ins, is, js, ks
        INTEGER, PARAMETER :: nSub = 16
        REAL(realk) :: xs, ys, zs

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

            CALL get_fieldptr(c, "C", igrid)
            CALL get_fieldptr(dx, "DX", igrid)
            CALL get_fieldptr(dy, "DY", igrid)
            CALL get_fieldptr(dz, "DZ", igrid)
            CALL get_fieldptr(ddx, "DDX", igrid)
            CALL get_fieldptr(ddy, "DDY", igrid)
            CALL get_fieldptr(ddz, "DDZ", igrid)

            c = 0.0_realk
            xm = minx
            DO i = 3, ii-2
                x = xm + 0.5_realk*ddx(i)
                ym = miny
                DO j = 3, jj-2
                    y = ym + 0.5_realk*ddy(j)
                    zm = minz
                    DO k = 3, kk-2
                        z = zm + 0.5_realk*ddz(k)

                        halfDiag = 0.5_realk*SQRT(ddx(i)**2 + ddy(j)**2 + ddz(k)**2)
                        dist = dist_func(x, y, z)

                        IF ( dist <= -halfDiag ) THEN
                            c(k,j,i) = 1.0_realk
                        ELSE IF ( dist >= halfDiag ) THEN
                            c(k,j,i) = 0.0_realk
                        ELSE 
                            ins = 0
                            DO is = 1, nSub
                                xs = xm + (REAL(is, realk) - 0.5_realk)*ddx(i)/REAL(nSub, realk)
                                DO js = 1, nSub
                                    ys = ym + (REAL(js, realk) - 0.5_realk)*ddy(j)/REAL(nSub, realk)
                                    DO ks = 1, nSub
                                        zs = zm + (REAL(ks, realk) - 0.5_realk)*ddz(k)/REAL(nSub, realk)
                                        IF ( dist_func(xs, ys, zs) < 0.0_realk ) THEN
                                            ins = ins + 1
                                        END IF
                                    END DO
                                END DO
                            END DO
                            c(k,j,i) = REAL(ins, realk)/REAL(nSub, realk)**3
                        END IF
                        zm = zm + ddz(k)
                    END DO
                    ym = ym + ddy(j)
                END DO
                xm = xm + ddx(i)
            END DO
        END DO

    END SUBROUTINE fill_c_dom

    !================================================================

    SUBROUTINE fill_c_bou(shp)
    !----------------------------------------------------------------
    !   What it does:
    !   Fill the volume fraction field c depending on the selected
    !   distance function. Only the boundaries of the domain are
    !   filled.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: shp

        ! Local variabels
        ! None

        CONTINUE

    END SUBROUTINE fill_c_bou

    !================================================================

    PURE FUNCTION dist_func_circle(x, y, xc, yc, r) RESULT(phi)
        REAL(realk), INTENT(in) :: x, y, xc, yc, r
        REAL(realk) :: phi
        phi = SQRT((x - xc)**2 + (y - yc)**2) - r
    END FUNCTION dist_func_circle

    !================================================================

    PURE FUNCTION dist_func_box(x, y, xc, yc, hx, hy) RESULT(phi)
        REAL(realk), INTENT(in) :: x, y, xc, yc, hx, hy
        REAL(realk) :: phi
        phi = MAX(ABS(x - xc) - hx, ABS(y - yc) - hy)
    END FUNCTION dist_func_box

    !================================================================

    PURE FUNCTION dist_func_ellipse(x, y, xc, yc, ra, rb, theta) RESULT(phi)
        REAL(realk), INTENT(in) :: x, y, xc, yc, ra, rb, theta
        REAL(realk) :: phi
        REAL(realk) :: ct, st, xr, yr
        ct = COS(theta)
        st = SIN(theta)
        xr =  ct*(x - xc) + st*(y - yc)
        yr = -st*(x - xc) + ct*(y - yc)
        phi = MIN(ra, rb)*(SQRT((xr/ra)**2 + (yr/rb)**2) - 1.0_realk)
    END FUNCTION dist_func_ellipse

END MODULE mph_test_mod