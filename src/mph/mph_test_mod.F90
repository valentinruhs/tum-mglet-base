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

    USE err_mod, ONLY: err_abort
    USE precision_mod, ONLY: realk, intk, pi
    USE grids_mod, ONLY: nmygrids, mygrids, get_mgdims, get_bbox
    USE fields_mod, ONLY: get_fieldptr
    USE mphcore_mod, ONLY: mphTst, mphInitErr, vofTol
    USE mph_plic_mod, ONLY: comp_c_stg

    IMPLICIT NONE(type, external)
    PRIVATE

    INTEGER(intk), PARAMETER :: tstCylTra=1, tstZalDis=2, tstRKoVor=3, tstUCylAd=4, tstACylAd=5, tstEllRec=6, tstOpCFl=7
    INTEGER(intk), PARAMETER :: shpCircle=1, shpZalesak=2, shpEllipse=3, shpPlane=4
    LOGICAL, PROTECTED :: frcVelFld

    PUBLIC :: init_mph_test, finish_mph_test

    TYPE :: shape_t
        INTEGER(intk) :: shp = shpCircle
        REAL(realk) :: xc = 0.0_realk, yc = 0.0_realk, zc = 0.0_realk
        REAL(realk) :: ra = 0.0_realk, rb = 0.0_realk, theta = 0.0_realk
        REAL(realk) :: lvl = 0.0_realk
        REAL(realk) :: slotW = 0.0_realk, slotH = 0.0_realk
    END TYPE shape_t

    ABSTRACT INTERFACE
        PURE FUNCTION dist_func_interface(x, y, z, s) RESULT(phi)
            IMPORT :: realk, shape_t
            REAL(realk), INTENT(in) :: x, y, z
            TYPE(shape_t), INTENT(in) :: s
            REAL(realk) :: phi
        END FUNCTION dist_func_interface
    END INTERFACE

CONTAINS

    SUBROUTINE init_mph_test()

        ! Subroutine arguments
        ! None

        ! Local variables
        INTEGER(intk) :: tst
        TYPE(shape_t) :: shape

        SELECT CASE ( mphTst )
        CASE ( "Cylinder Translation" )
            tst = tstCylTra
            shape = shape_t(shp=shpCircle, &
                xc=0.5_realk, yc=0.5_realk, ra=0.15_realk)
            frcVelFld = .TRUE.

        CASE ( "Zalesak Disk" )
            tst = tstZalDis
            shape = shape_t(shp=shpZalesak, &
                xc=0.5_realk, yc=0.5_realk, ra=0.15_realk, &
                slotW=0.05_realk, slotH=0.25_realk)
            frcVelFld = .TRUE.

        CASE ( "Rider-Kothe Vortex" )
            tst = tstRKoVor
            shape = shape_t(shp=shpCircle, &
                xc=0.5_realk, yc=0.75_realk, ra=0.15_realk)
            frcVelFld = .TRUE.

        CASE ( "Uniform Cylinder Advection" )
            tst = tstUCylAd
            shape = shape_t(shp=shpCircle, &
                xc=0.2_realk, yc=0.2_realk, ra=0.1_realk)
            frcVelFld = .FALSE.

        CASE ( "Abrupt Cylinder Advection" )
            tst = tstACylAd
            shape = shape_t(shp=shpCircle, &
                xc=0.2_realk, yc=0.2_realk, ra=0.1_realk)
            frcVelFld = .FALSE.

        CASE ( "Ellipse Reconstruction" )
            tst = tstEllRec
            shape = shape_t(shp=shpEllipse, &
                xc=0.5_realk, yc=0.5_realk, &
                ra=0.3464_realk, rb=0.1414_realk, theta=pi)
            frcVelFld = .TRUE.

        CASE ( "Open Channel Flow" )
            tst = tstOpCFl
            shape = shape_t(shp=shpPlane, lvl=0.0_realk)
            frcVelFld = .FALSE.

        CASE DEFAULT
            CALL err_abort(mphInitErr, "unknown test case.", __FILE__, __LINE__)
        END SELECT

        CALL fill_c_dom(dist_func, shape)
        CALL fill_c_bou(shape%shp)
        CALL init_vel(tst)

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

    SUBROUTINE fill_c_dom(dist_func, s)
    !----------------------------------------------------------------
    !   What it does:
    !   Fill the volume fraction field c depending on the selected
    !   distance function. Only the inside of the domain is filled.
    !----------------------------------------------------------------

        ! Subroutine arguments
        PROCEDURE(dist_func_interface) :: dist_func
        TYPE(shape_t), INTENT(in) :: s

        ! Local variables
        INTEGER(intk) :: n, igrid, kk, jj, ii
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz
        REAL(realk), POINTER, CONTIGUOUS :: c(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        INTEGER(intk) :: i, j, k
        REAL(realk) :: xm, x, ym, y, zm, z, halfDiag, dist
        INTEGER(intk) :: ins, is, js, ks
        INTEGER, PARAMETER :: nSub = 64
        REAL(realk) :: xs, ys, zs

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

            CALL get_fieldptr(c, "C", igrid)
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
                        dist = dist_func(x, y, z, s)

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
                                        IF ( dist_func(xs, ys, zs, s) < 0.0_realk ) THEN
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

    SUBROUTINE init_vel(tst)
    !----------------------------------------------------------------
    !   What it does:
    !   Imposes the prescribed velocity field of the selected test
    !   case. Does nothing for tests whose velocity field is produced
    !   by the flow solver.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: tst

        ! Local variables
        ! None

        SELECT CASE ( tst )
        CASE ( tstUCylAd )
            CALL set_vel_uni(0.016_realk, 0.016_realk, 0.0_realk)
        CASE ( tstACylAd )
            CALL set_vel_c(0.016_realk, 0.016_realk, 0.0_realk, h=2)
        CASE ( tstEllRec )
            CALL set_vel_uni(0.0_realk, 0.0_realk, 0.0_realk)
        CASE DEFAULT
            CALL err_abort(mphInitErr, "no velocity field for this test.", __FILE__, __LINE__)
        END SELECT

    END SUBROUTINE init_vel

    !================================================================

    SUBROUTINE set_vel_uni(uc, vc, wc)
    !----------------------------------------------------------------
    !   What it does:
    !   Sets u, v and w to a uniform value on every grid. The whole
    !   array is written, ghost layers included, hence no connect or
    !   parent is needed afterwards.
    !----------------------------------------------------------------

        ! Subroutine arguments
        REAL(realk), INTENT(in) :: uc, vc, wc

        ! Local variables
        INTEGER(intk) :: n, igrid
        REAL(realk), POINTER, CONTIGUOUS :: u(:,:,:), v(:,:,:), w(:,:,:)

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_fieldptr(u, "U", igrid)
            CALL get_fieldptr(v, "V", igrid)
            CALL get_fieldptr(w, "W", igrid)

            u = uc
            v = vc
            w = wc
        END DO

    END SUBROUTINE set_vel_uni

    !================================================================

    SUBROUTINE set_vel_c(uc, vc, wc, h)
    !----------------------------------------------------------------
    !   What it does:
    !   Sets u, v and w to a uniform value on every grid. The whole
    !   array is written, ghost layers included, hence no connect or
    !   parent is needed afterwards.
    !----------------------------------------------------------------

        ! Subroutine arguments
        REAL(realk), INTENT(in) :: uc, vc, wc
        INTEGER(intk), INTENT(in) :: h

        ! Local variables
        INTEGER(intk) :: n, igrid
        INTEGER(intk) :: kk, jj, ii, q
        REAL(realk), POINTER, CONTIGUOUS :: cS1(:,:,:), cS2(:,:,:), cS3(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: u(:,:,:), v(:,:,:), w(:,:,:)

        DO n = 1, nmygrids
            igrid = mygrids(n)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_fieldptr(cS1, "CS1", igrid)
            CALL get_fieldptr(cS2, "CS2", igrid)
            CALL get_fieldptr(cS3, "CS3", igrid)
            CALL get_fieldptr(u, "U", igrid)
            CALL get_fieldptr(v, "V", igrid)
            CALL get_fieldptr(w, "W", igrid)

            DO q = 1, 3
                CALL comp_c_stg(q)
            END DO
            CALL set_vel_c_grd(kk, jj, ii, cS1, cS2, cS3, u, v, w, uc, vc, wc, h)
        END DO

    END SUBROUTINE set_vel_c

    !================================================================

    SUBROUTINE set_vel_c_grd(kk, jj, ii, cS1, cS2, cS3, u, v, w, uc, vc, wc, h)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: cS1(kk, jj, ii), cS2(kk, jj, ii), cS3(kk, jj, ii)
        REAl(realk), INTENT(inout) :: u(kk, jj, ii), w(kk, jj, ii), v(kk, jj, ii)
        REAL(realk), INTENT(in) :: uc, vc, wc
        INTEGER(intk), INTENT(in) :: h

        ! Local variables
        INTEGER(intk) :: k, j, I

        DO i = 2, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    IF ( cS1(k,j,i) > vofTol ) THEN
                        u(k-h:k+h,j-h:j+h,i-h:i+h) = cS1(k-h:k+h,j-h:j+h,i-h:i+h)*uc
                    END IF
                END DO
            END DO
        END DO

        DO i = 3, ii-2
            DO j = 2, jj-2
                DO k = 3, kk-2
                    IF ( cS2(k,j,i) > vofTol ) THEN
                        v(k-h:k+h,j-h:j+h,i-h:i+h) = cS2(k-h:k+h,j-h:j+h,i-h:i+h)*vc
                    END IF
                END DO
            END DO
        END DO

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 2, kk-2
                    IF ( cS3(k,j,i) > vofTol ) THEN
                        w(k-h:k+h,j-h:j+h,i-h:i+h) = cS3(k-h:k+h,j-h:j+h,i-h:i+h)*wc
                    END IF
                END DO
            END DO
        END DO

    END SUBROUTINE set_vel_c_grd

    !================================================================

    ! SUBROUTINE frc_vel_fld()
    ! !----------------------------------------------------------------
    ! !   What it does:
    ! !   Re-imposes the prescribed velocity field. Called once per
    ! !   time step so that neither the flow solver nor upd_vel_stg
    ! !   can alter it.
    ! !----------------------------------------------------------------

    !     ! Subroutine arguments
    !     ! None

    !     ! Local variables
    !     ! None

    !     CALL init_vel(mphTstId)

    ! END SUBROUTINE frc_vel_fld

    !================================================================

    PURE FUNCTION dist_func(x, y, z, s) RESULT(phi)
        REAL(realk), INTENT(in) :: x, y, z
        TYPE(shape_t), INTENT(in) :: s
        REAL(realk) :: phi

        SELECT CASE ( s%shp )
        CASE ( shpCircle )
            phi = dist_func_circle(x, y, s%xc, s%yc, s%ra)
        CASE ( shpZalesak )
            phi = MAX(dist_func_circle(x, y, s%xc, s%yc, s%ra), &
                     -dist_func_box(x, y, s%xc, s%yc - s%ra + 0.5_realk*s%slotH, &
                                    0.5_realk*s%slotW, 0.5_realk*s%slotH))
        CASE ( shpEllipse )
            phi = dist_func_ellipse(x, y, s%xc, s%yc, s%ra, s%rb, s%theta)
        CASE ( shpPlane )
            phi = y - s%lvl
        CASE DEFAULT
            phi = 1.0_realk
        END SELECT
    END FUNCTION dist_func

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