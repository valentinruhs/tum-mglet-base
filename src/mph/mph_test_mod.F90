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
    USE mph_plic_mod, ONLY: comp_c_loc

    IMPLICIT NONE(type, external)
    PRIVATE

    INTEGER(intk), PARAMETER :: tstCylTra=1, tstZalDis=2, tstRKoVor=3, tstUCylAd=4, tstACylAd=5, tstEllRec=6, tstOpCFl=7
    INTEGER(intk), PARAMETER :: shpCircle=1, shpZalesak=2, shpEllipse=3, shpPlane=4
    LOGICAL, PROTECTED :: frcVelFld, isRevTst
    INTEGER(intk), PROTECTED :: tstId
    REAL(realk), PROTECTED :: circumf

    PUBLIC :: init_mph_test, finish_mph_test, frc_vel_fld, isRevTst, comp_eGeo, comp_eIfc, circumf, shape, frcVelFld

    TYPE :: shape_t
        INTEGER(intk) :: shp = shpCircle
        REAL(realk) :: xc = 0.0_realk, yc = 0.0_realk, zc = 0.0_realk
        REAL(realk) :: ra = 0.0_realk, rb = 0.0_realk, theta = 0.0_realk
        REAL(realk) :: lvl = 0.0_realk
        REAL(realk) :: slotW = 0.0_realk, slotH = 0.0_realk
    END TYPE shape_t

    TYPE(shape_t), PROTECTED :: shape

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
        ! None

        SELECT CASE ( mphTst )
        CASE ( "Cylinder Translation" )
            tstId = tstCylTra
            shape = shape_t(shp=shpCircle, &
                xc=0.5_realk, yc=0.5_realk, ra=0.15_realk)
            frcVelFld = .TRUE.
            isRevTst = .TRUE.

        CASE ( "Zalesak Disk" )
            tstId = tstZalDis
            shape = shape_t(shp=shpZalesak, &
                xc=0.5_realk, yc=0.75_realk, ra=0.15_realk, &
                slotW=0.05_realk, slotH=0.25_realk)
            frcVelFld = .TRUE.
            isRevTst = .TRUE.

        CASE ( "Rider-Kothe Vortex" )
            tstId = tstRKoVor
            shape = shape_t(shp=shpCircle, &
                xc=0.5_realk, yc=0.75_realk, ra=0.15_realk)
            frcVelFld = .TRUE.
            isRevTst = .TRUE.

        CASE ( "Uniform Cylinder Advection" )
            tstId = tstUCylAd
            shape = shape_t(shp=shpCircle, &
                xc=0.2_realk, yc=0.2_realk, ra=0.1_realk)
            frcVelFld = .FALSE.

        CASE ( "Abrupt Cylinder Advection" )
            tstId = tstACylAd
            shape = shape_t(shp=shpCircle, &
                xc=0.2_realk, yc=0.2_realk, ra=0.1_realk)
            frcVelFld = .FALSE.

        CASE ( "Ellipse Reconstruction" )
            tstId = tstEllRec
            shape = shape_t(shp=shpEllipse, &
                xc=0.5_realk, yc=0.5_realk, &
                ra=0.3464_realk, rb=0.1414_realk, theta=pi)
            frcVelFld = .TRUE.

        CASE ( "Open Channel Flow" )
            tstId = tstOpCFl
            shape = shape_t(shp=shpPlane, lvl=0.0_realk)
            frcVelFld = .FALSE.
            isRevTst = .FALSE.

        CASE DEFAULT
            CALL err_abort(mphInitErr, "unknown test case.", __FILE__, __LINE__)

        END SELECT

        circumf = circ_func(shape)
        CALL fill_c_dom(dist_func, shape)
        CALL fill_c_bou(shape%shp)
        CALL init_vel()

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
        REAL(realk) :: xMi, x, yMi, y, zMi, z, halfDiag, dist
        INTEGER(intk) :: ins, is, js
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
            xMi = minx
            DO i = 3, ii-2
                x = xMi + 0.5_realk*ddx(i)
                yMi = miny
                DO j = 3, jj-2
                    y = yMi + 0.5_realk*ddy(j)
                    zMi = minz
                    DO k = 3, kk-2
                        z = zMi + 0.5_realk*ddz(k)

                        halfDiag = 0.5_realk*SQRT(ddx(i)**2 + ddy(j)**2 + ddz(k)**2)
                        dist = dist_func(x, y, z, s)

                        IF ( dist <= -halfDiag ) THEN
                            c(k,j,i) = 1.0_realk
                        ELSE IF ( dist >= halfDiag ) THEN
                            c(k,j,i) = 0.0_realk
                        ELSE 
                            ins = 0
                            DO is = 1, nSub
                                xs = xMi + (REAL(is, realk) - 0.5_realk)*ddx(i)/REAL(nSub, realk)
                                DO js = 1, nSub
                                    ys = yMi + (REAL(js, realk) - 0.5_realk)*ddy(j)/REAL(nSub, realk)
                                    zs = zMi + 0.5_realk*ddz(k)
                                    IF ( dist_func(xs, ys, zs, s) < 0.0_realk ) THEN
                                        ins = ins + 1
                                    END IF
                                END DO
                            END DO
                            c(k,j,i) = REAL(ins, realk)/REAL(nSub, realk)**2
                        END IF
                        zMi = zMi + ddz(k)
                    END DO
                    yMi = yMi + ddy(j)
                END DO
                xMi = xMi + ddx(i)
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

    SUBROUTINE init_vel(dt, itstep)
    !----------------------------------------------------------------
    !   What it does:
    !   Imposes the prescribed velocity field of the selected test
    !   case. Does nothing for tests whose velocity field is produced
    !   by the flow solver.
    !----------------------------------------------------------------

        ! Subroutine arguments
        REAL(realk), INTENT(in), OPTIONAL :: dt
        INTEGER(intk), INTENT(in), OPTIONAL :: itstep

        ! Local variables
        REAL(realk) :: t

        IF ( PRESENT(itstep) ) THEN
            t = dt*itstep
        ELSE
            t = 0.0_realk
        END IF

        SELECT CASE ( tstId )
        CASE ( tstZalDis )
            CALL set_vel_cav()
        CASE ( tstRKoVor )
            CALL set_vel_vtx(t)
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
    !   
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

    SUBROUTINE set_vel_vtx(t)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        REAL(realk), INTENT(in) :: t

        ! Local variables
        INTEGER(intk) :: n, igrid
        INTEGER(intk) :: kk, jj, ii
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz
        REAL(realk), POINTER, CONTIGUOUS :: u(:,:,:), v(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:)

        DO n = 1, nmygrids
            igrid = mygrids(n)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)
            CALL get_fieldptr(u, "U", igrid)
            CALL get_fieldptr(v, "V", igrid)
            CALL get_fieldptr(ddx, "DDX", igrid)
            CALL get_fieldptr(ddy, "DDY", igrid)

            CALL set_vel_vtx_grid(kk, jj, ii, u, v, ddx, ddy, minx, miny, t)
        END DO

    END SUBROUTINE set_vel_vtx

    !================================================================

    SUBROUTINE set_vel_vtx_grid(kk, jj, ii, u, v, ddx, ddy, minx, miny, t)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAl(realk), INTENT(inout) :: u(kk, jj, ii), v(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj)
        REAL(realk), INTENT(in) :: minx, miny, t

        ! Local variables
        INTEGER(intk) :: k, j, i
        REAL(realk) :: psi(jj, ii)
        REAL(realk) :: x, y, xMi, yMi

        xMi = minx - 0.5*ddx(2)
        DO i = 2, ii-2
            x = xMi + 0.5_realk*ddx(i)
            yMi = miny - 0.5*ddy(2)
            DO j = 2, jj-2
                y = yMi + 0.5_realk*ddy(j)

                psi(j,i) = 1.0_realk/pi*COS(pi*t/2.0_realk)* &
                    SIN(pi*x)**2 * SIN(pi*y)**2

                yMi = yMi + ddy(j)
            END DO
            xMi = xMi + ddx(i)
        END DO

        DO i = 2, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    u(k,j,i) = -(psi(j,i) - psi(j-1,i))/ddy(j)
                END DO
            END DO
        END DO

        DO i = 3, ii-2
            DO j = 2, jj-2
                DO k = 3, kk-2
                    v(k,j,i) = (psi(j,i) - psi(j,i-1))/ddx(i)
                END DO
            END DO
        END DO

    END SUBROUTINE set_vel_vtx_grid

    !================================================================

    SUBROUTINE set_vel_cav()
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        ! None

        ! Local variables
        INTEGER(intk) :: n, igrid
        INTEGER(intk) :: kk, jj, ii
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz
        REAL(realk), POINTER, CONTIGUOUS :: u(:,:,:), v(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:)

        DO n = 1, nmygrids
            igrid = mygrids(n)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)
            CALL get_fieldptr(u, "U", igrid)
            CALL get_fieldptr(v, "V", igrid)
            CALL get_fieldptr(ddx, "DDX", igrid)
            CALL get_fieldptr(ddy, "DDY", igrid)

            CALL set_vel_cav_grid(kk, jj, ii, u, v, ddx, ddy, minx, miny)
        END DO

    END SUBROUTINE set_vel_cav

    !================================================================

    SUBROUTINE set_vel_cav_grid(kk, jj, ii, u, v, ddx, ddy, minx, miny)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAl(realk), INTENT(inout) :: u(kk, jj, ii), v(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj)
        REAL(realk), INTENT(in) :: minx, miny

        ! Local variables
        INTEGER(intk) :: k, j, i
        REAL(realk) :: omega
        REAL(realk) :: x, y, xMi, yMi

        omega = 2.0_realk*pi/6.28_realk

        DO i = 2, ii-2
            yMi = miny
            DO j = 3, jj-2
                y = yMi + 0.5_realk*ddy(j)
                DO k = 3, kk-2
                    u(k,j,i) = -omega*(y - 0.5_realk)
                END DO
                yMi = yMi + ddy(j)
            END DO
        END DO

        xMi = minx
        DO i = 3, ii-2
            x = xMi + 0.5_realk*ddx(i)
            DO j = 2, jj-2
                DO k = 3, kk-2
                    v(k,j,i) = omega*(x - 0.5_realk)
                END DO
            END DO
            xMi = xMi + ddx(i)
        END DO

    END SUBROUTINE set_vel_cav_grid

    !================================================================

    SUBROUTINE frc_vel_fld(dt, itstep)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        REAL(realk), INTENT(in) :: dt
        INTEGER(intk), INTENT(in) :: itstep

        ! Local variables
        ! None

        CALL init_vel(dt, itstep)

    END SUBROUTINE frc_vel_fld

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

    PURE FUNCTION circ_func(s) RESULT(phi)
        TYPE(shape_t), INTENT(in) :: s
        REAL(realk) :: phi

        SELECT CASE ( s%shp )
        CASE ( shpCircle )
            phi = circ_func_circle(s%ra)
        CASE ( shpZalesak )
            phi = circ_func_circle(s%ra) - 2.0_realk*s%ra*ASIN(0.5_realk*s%slotW/s%ra) &
                + 2.0_realk*(s%slotH - s%ra + SQRT(s%ra**2 - 0.25_realk*s%slotW**2)) + s%slotW
        CASE ( shpEllipse )
            phi = circ_func_ellipse(s%ra, s%rb)
        CASE DEFAULT
            phi = 1.0_realk
        END SELECT
    END FUNCTION circ_func

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

    !================================================================

    PURE FUNCTION circ_func_circle(ra) RESULT(phi)
        REAL(realk), INTENT(in) :: ra
        REAL(realk) :: phi
        phi = 2.0_realk*pi*ra
    END FUNCTION circ_func_circle

    !================================================================

    PURE FUNCTION circ_func_ellipse(ra, rb) RESULT(phi)
        REAL(realk), INTENT(in) :: ra, rb
        REAL(realk) :: phi
        REAL(realk) :: lbd
        lbd = (ra - rb)/(ra + rb)
        phi = (ra + rb)*pi*(1.0_realk + (3.0_realk*lbd**2)/ &
            (10.0_realk + SQRT(4.0_realk - 3.0_realk*lbd**2)))
    END FUNCTION circ_func_ellipse

    !================================================================

    FUNCTION comp_eGeo(ddx, ddy, ddz, cAct, cRef) RESULT(eGeo)
        REAL(realk) :: ddx, ddy, ddz, cAct, cRef
        REAL(realk) :: eGeo

        eGeo = ddx*ddy*ddz*ABS(cAct - cRef)

    END FUNCTION comp_eGeo

    !================================================================

    RECURSIVE FUNCTION comp_eIfc(xMi, yMi, zMi, ddx, ddy, ddz, &
        normx, normy, normz, alpha, c, isIfc, s) RESULT(eIfc)
        REAL(realk), INTENT(in) :: xMi, yMi, zMi, ddx, ddy, ddz
        REAL(realk), INTENT(in) :: normx, normy, normz, alpha, c
        REAL(realk), INTENT(in) :: isIfc
        TYPE(shape_t), INTENT(in) :: s
        REAL(realk) :: eIfc

        INTEGER(intk), PARAMETER :: nSub = 256
        REAL(realk) :: vol, cAppr, dist, halfDiag, xs, ys, zs
        INTEGER(intk) :: is, js, ins
        LOGICAL :: inExac, inAppr

        vol = ddx*ddy*ddz
        cAppr = c

        halfDiag = 0.5_realk*SQRT(ddx**2 + ddy**2 + ddz**2)
        dist = dist_func(xMi + 0.5_realk*ddx, yMi + 0.5_realk*ddy, zMi + 0.5_realk*ddz, s)

        IF ( dist <= -halfDiag ) THEN
            eIfc = (1.0_realk - cAppr)*vol
        ELSE IF ( dist >= halfDiag ) THEN
            eIfc = cAppr*vol
        ELSE
            ins = 0
            DO is = 1, nSub
                xs = (REAL(is, realk) - 0.5_realk)*ddx/REAL(nSub, realk)
                DO js = 1, nSub
                    ys = (REAL(js, realk) - 0.5_realk)*ddy/REAL(nSub, realk)
                    zs = 0.5_realk*ddz
                    inExac = dist_func(xMi+xs, yMi+ys, zMi+zs, s) < 0.0_realk
                    IF ( isIfc > 0.0_realk ) THEN
                        inAppr = normx*xs + normy*ys + normz*zs < alpha
                    ELSE
                        inAppr = c > 0.5_realk
                    END IF
                    IF (inExac .NEQV. inAppr) ins = ins + 1
                END DO
            END DO
            eIfc = vol*REAL(ins,realk)/REAL(nSub,realk)**2
        END IF

    END FUNCTION comp_eIfc

END MODULE mph_test_mod