!====================================================================
!  Module: mph_test_mod
!
!   Responsibilities:
!   - Defines the test cases of the multi-phase solver, selected by
!     /multiphase/test:
!     - Cylinder Translation, Zalesak Disk, Rider-Kothe Vortex
!       (prescribed velocity field)
!     - Uniform and Abrupt Cylinder Advection (initial velocity
!       field)
!     - Open Channel Flow (initial velocity field and hydrostatic
!       pressure)
!   - Initializes c from the signed distance function of the test
!     shape as well as velocity and pressure (not on restart)
!   - Re-imposes prescribed velocity fields during the simulation
!     (frcVelFld)
!   - Provides the exact shape (distance function, interface
!     length, area) for the error norms in mph_chk_mod
!
!   Author:      Valentin Ruhs
!   e-Mail:      valentin.ruhs@gmx.de
!   Created:     2026-09
!
!   Sources:
!   H. Reichardt, “Vollständige Darstellung der turbulenten
!   Geschwindigkeitsverteilung in glatten Leitungen,” Zeitschrift für
!   Angewandte Mathematik und Mechanik, vol. 31, no. 7, pp. 208–219,
!   Jan. 1951, doi: 10.1002/zamm.19510310704.
!
!   W. Schoppa and F. Hussain, “Coherent structure dynamics in
!   near-wall turbulence,” Fluid Dynamics Research, vol. 26, no. 2,
!   pp. 119–139, Feb. 2000, doi: 10.1016/S0169-5983(99)00018-0.
!
!   E. de Villiers, “The Potential of Large Eddy Simulation for the
!   Modelling of Wall Bounded Flows,” Degree of Doctor of Philosophy
!   of the University of London and Diploma of Imperial College,
!   Imperial College London, 2006.
!
!====================================================================

MODULE mph_test_mod

    USE err_mod, ONLY: err_abort
    USE precision_mod, ONLY: realk, intk, pi
    USE grids_mod, ONLY: nmygrids, mygrids, get_mgdims, get_bbox
    USE fields_mod, ONLY: get_fieldptr, set_field
    USE fort7_mod, ONLY: dcont

    USE mphcore_mod, ONLY: mphTst, mphInitErr, vofTol
    USE mph_plic_mod, ONLY: comp_c_stg, comp_ifc
    USE mph_plic_mod, ONLY: comp_c_loc
    USE mphcore_mod, ONLY: gmol1, rho1, rho2, grav
    USE mph_xfer_mod, ONLY: rstr, prlg


    IMPLICIT NONE(type, external)
    PRIVATE

    INTEGER(intk), PARAMETER :: tstCylTra=1, tstZalDis=2, tstRKoVor=3, tstUCylAd=4, tstACylAd=5, tstOpCFl=6
    INTEGER(intk), PARAMETER :: shpCircle=1, shpZalesak=2, shpPlane=3, shpSphere=4
    INTEGER, PARAMETER :: nSub = 64
    LOGICAL, PROTECTED :: frcVelFld, isRevTst
    INTEGER(intk), PROTECTED :: tstId
    REAL(realk), PROTECTED :: circumf, area

    PUBLIC :: init_mph_test, finish_mph_test, frc_vel_fld, isRevTst, comp_abs_res, circumf, &
        area, shape, frcVelFld

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
    !----------------------------------------------------------------
    !   Short: Select test case; initialize c, interface, velocity
    !          and pressure.
    !
    !   The test case is selected by /multiphase/test: Cylinder
    !   Translation, Zalesak Disk, Rider-Kothe Vortex, Uniform
    !   Cylinder Advection, Abrupt Cylinder Advection or Open Channel
    !   Flow. On a restart (dcont), c, velocity and pressure are
    !   taken from the restart file; only restriction, prolongation
    !   and interface reconstruction are performed.
    !----------------------------------------------------------------

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        SELECT CASE ( mphTst )
        CASE ( "Cylinder Translation" )
            tstId = tstCylTra
            shape = shape_t(shp=shpCircle, &
                xc=0.5_realk, yc=0.5_realk, ra=5.5_realk/80.0_realk)
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
            isRevTst = .TRUE.

        CASE ( "Abrupt Cylinder Advection" )
            tstId = tstACylAd
            shape = shape_t(shp=shpSphere, &
                xc=0.2_realk, yc=0.2_realk, zc=0.5_realk, ra=0.15_realk)
            frcVelFld = .FALSE.
            isRevTst = .FALSE.

        CASE ( "Open Channel Flow" )
            tstId = tstOpCFl
            shape = shape_t(shp=shpPlane, lvl=1.0_realk)
            frcVelFld = .FALSE.
            isRevTst = .FALSE.

        CASE DEFAULT
            CALL err_abort(mphInitErr, "unknown test case.", __FILE__, __LINE__)
        END SELECT

        circumf = circf_func(shape)
        area = area_func(shape)
        IF ( .NOT. dcont ) CALL fill_c_dom(dist_func, shape)
        CALL rstr(fldName="C", flag="D")
        CALL prlg(fldName="C")
        CALL comp_ifc()
        IF ( .NOT. dcont ) CALL init_vel()
        IF ( .NOT. dcont ) CALL init_pre()
        CALL rstr(fldName="P", flag="P")
        CALL prlg(fldName="P")

    END SUBROUTINE init_mph_test

    !================================================================

    SUBROUTINE finish_mph_test()
    !----------------------------------------------------------------
    !   Short: Finalize module (no operation).
    !----------------------------------------------------------------

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        CONTINUE

    END SUBROUTINE finish_mph_test

    !================================================================

    SUBROUTINE fill_c_dom(dist_func, s)
    !----------------------------------------------------------------
    !   Short: Initialize c from a signed distance function
    !          by sub-sampling.
    !
    !   Only the interior cells are filled. Cells further than half a
    !   cell diagonal from the interface are set to 0 or 1; cells cut
    !   by the interface are sub-sampled with nSub x nSub points in
    !   the x-y plane (2D shapes).
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
        INTEGER(intk) :: ins, is, js, ks
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
                                    DO ks = 1, nSub
                                        zs = zMi + (REAL(ks, realk) - 0.5_realk)*ddz(k)/REAL(nSub, realk)
                                        IF ( dist_func(xs, ys, zs, s) < 0.0_realk ) THEN
                                            ins = ins + 1
                                        END IF
                                    END DO
                                END DO
                            END DO
                            c(k,j,i) = REAL(ins, realk)/REAL(nSub, realk)**3
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

    SUBROUTINE init_vel(time)
    !----------------------------------------------------------------
    !   Short: Impose the prescribed velocity field of the selected
    !          test case.
    !
    !   time: evaluation time of time-dependent fields (default 0).
    !   For the open channel flow only the initial field is set.
    !----------------------------------------------------------------

        ! Subroutine arguments
        REAL(realk), INTENT(in), OPTIONAL :: time

        ! Local variables
        REAL(realk) :: t

        IF ( PRESENT(time) ) THEN
            t = time
        ELSE
            t = 0.0_realk
        END IF

        SELECT CASE ( tstId )
        CASE ( tstCylTra )
            CALL set_vel_tra(0.0125_realk, t)
        CASE ( tstZalDis )
            CALL set_vel_cav()
        CASE ( tstRKoVor )
            CALL set_vel_vtx(t)
        CASE ( tstUCylAd )
            CALL set_vel_uni(0.016_realk, 0.016_realk, 0.0_realk)
        CASE ( tstACylAd )
            CALL set_vel_c(10.0_realk, 10.0_realk, 0.0_realk, h=1)
        CASE ( tstOpCFl )
            CALL set_vel_cha(150.0_realk, 1.0_realk, 2.0_realk*pi, pi, "TURB")
        CASE DEFAULT
            CALL err_abort(mphInitErr, "no velocity field for this test.", __FILE__, __LINE__)
        END SELECT

    END SUBROUTINE init_vel

    !================================================================

    SUBROUTINE set_vel_tra(mag, t)
    !----------------------------------------------------------------
    !   Short: Set piecewise constant translation velocity
    !          (cylinder translation).
    !----------------------------------------------------------------

        ! Subroutine arguments
        REAL(realk), INTENT(in) :: mag, t

        ! Local variables
        INTEGER(intk) :: n, igrid
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: u(:,:,:), v(:,:,:)

        DO n = 1, nmygrids
            igrid = mygrids(n)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_fieldptr(u, "U", igrid)
            CALL get_fieldptr(v, "V", igrid)

            CALL set_vel_tra_grid(kk, jj, ii, u, v, mag, t)
        END DO

    END SUBROUTINE set_vel_tra

    !================================================================

    SUBROUTINE set_vel_tra_grid(kk, jj, ii, u, v, mag, t)
    !----------------------------------------------------------------
    !   Short: Grid kernel of set_vel_tra.
    !
    !   Seven phases of duration tPha: six with magnitude mag and
    !   fixed directions, the seventh with the negative sum of the
    !   first six, so that the cylinder returns to its initial
    !   position. Zero velocity afterwards.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: u(kk, jj, ii), v(kk, jj, ii)
        REAL(realk), INTENT(in) :: mag, t

        ! Local variables
        REAL(realk), PARAMETER :: tPha = 200.0_realk
        REAL(realk) :: fac(6), uPha(7), vPha(7)
        INTEGER(intk) :: iPha

        fac = [0.7686000_realk, 0.5968000_realk, 0.1035700_realk, &
               0.5514000_realk, 0.2242010_realk, 0.2512981_realk]

        uPha(1:6) = COS(2.0_realk*pi*fac)*mag
        vPha(1:6) = SIN(2.0_realk*pi*fac)*mag
        uPha(7) = -SUM(uPha(1:6))
        vPha(7) = -SUM(vPha(1:6))
        iPha = FLOOR(t/tPha) + 1

        IF ( iPha <= 7 ) THEN
            u = uPha(iPha)
            v = vPha(iPha)
        ELSE
            u = 0.0_realk
            v = 0.0_realk
        END IF

    END SUBROUTINE set_vel_tra_grid

    !================================================================

    SUBROUTINE set_vel_uni(uc, vc, wc)
    !----------------------------------------------------------------
    !   Short: Set uniform velocity on all grids incl. ghost layers.
    !
    !   The whole array is written, ghost layers included, hence no
    !   connect or parent is needed afterwards.
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
    !   Short: Set velocity only in and around cells containing
    !          phase 1.
    !
    !   The velocity (uc, vc, wc) is set in a cube of half width h
    !   around every staggered cell with cS > vofTol; elsewhere the
    !   velocity is unchanged.
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
    !   Short: Grid kernel of set_vel_c.
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
                        u(k-h:k+h,j-h:j+h,i-h:i+h) = uc
                    END IF
                END DO
            END DO
        END DO

        DO i = 3, ii-2
            DO j = 2, jj-2
                DO k = 3, kk-2
                    IF ( cS2(k,j,i) > vofTol ) THEN
                        v(k-h:k+h,j-h:j+h,i-h:i+h) = vc
                    END IF
                END DO
            END DO
        END DO

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 2, kk-2
                    IF ( cS3(k,j,i) > vofTol ) THEN
                        w(k-h:k+h,j-h:j+h,i-h:i+h) = wc
                    END IF
                END DO
            END DO
        END DO

    END SUBROUTINE set_vel_c_grd

    !================================================================

    SUBROUTINE set_vel_vtx(t)
    !----------------------------------------------------------------
    !   Short: Set time-dependent Rider-Kothe vortex velocity.
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
    !   Short: Grid kernel: velocities from the stream function.
    !
    !   psi = 1/pi*cos(pi*t/8)*sin(pi*x)^2*sin(pi*y)^2 is evaluated
    !   at the cell corners of the x-y plane (exact positions on
    !   equidistant grids). u and v follow from differences of psi
    !   and are therefore discretely divergence-free.
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

                psi(j,i) = 1.0_realk/pi*COS(pi*t/8.0_realk)* &
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
    !   Short: Set solid body rotation velocity (Zalesak disk).
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
    !   Short: Grid kernel of set_vel_cav.
    !
    !   Rotation about (0.5, 0.5) with omega = 2*pi/6.28, i.e. one
    !   revolution in t = 6.28.
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

    SUBROUTINE set_vel_cha(ReTau, h, lx, lz, flowStateFlag)
    !----------------------------------------------------------------
    !   Short: Set initial open channel velocity on all grids with
    !          grid-wise random seed.
    !
    !   lx, lz are the (periodic) domain lengths in x and z. The seed
    !   depends on the grid number, so that the perturbations are
    !   reproducible and independent of the distribution of the grids
    !   on the processes.
    !----------------------------------------------------------------

        ! Subroutine arguments
        REAL(realk), INTENT(in) :: ReTau, h, lx, lz
        CHARACTER(len=4), INTENT(in) :: flowStateFlag

        ! Local variables
        INTEGER(intk) :: n, igrid
        INTEGER(intk) :: kk, jj, ii
        INTEGER :: nSeed, s
        INTEGER, ALLOCATABLE :: seed(:)
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz
        REAL(realk), POINTER, CONTIGUOUS :: u(:,:,:), w(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)

        CALL RANDOM_SEED(size=nSeed)
        ALLOCATE(seed(nSeed))

        DO n = 1, nmygrids
            igrid = mygrids(n)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)
            CALL get_fieldptr(u, "U", igrid)
            CALL get_fieldptr(w, "W", igrid)
            CALL get_fieldptr(ddx, "DDX", igrid)
            CALL get_fieldptr(ddy, "DDY", igrid)
            CALL get_fieldptr(ddz, "DDZ", igrid)

            seed = [(123457 + 7919*INT(igrid) + 104729*s, s = 1, nSeed)]
            CALL RANDOM_SEED(put=seed)

            CALL set_vel_cha_grid(kk, jj, ii, u, w, ddx, ddy, ddz, &
                minx, miny, minz, ReTau, h, lx, lz, flowStateFlag)
        END DO

    END SUBROUTINE set_vel_cha

    !================================================================

    SUBROUTINE set_vel_cha_grid(kk, jj, ii, u, w, ddx, ddy, ddz, &
            minx, miny, minz, ReTau, h, lx, lz, flowStateFlag)
    !----------------------------------------------------------------
    !   Short: Grid kernel: turbulent (Reichardt profile and
    !          streaks), laminar or zero velocity.
    !
    !   "TURB": Reichardt mean profile plus near-wall streaks (u) and
    !           spanwise perturbations (w), both with a random
    !           amplitude variation of +-devAmp.
    !   "LAMI": laminar profile u+ = y+ - y+^2/(2*Re_tau).
    !   "NONE": u = 0.
    !
    !   y is limited to the interface height h, i.e. the gas phase
    !   gets the velocity at the surface.
    !
    !   Sources:
    !   W. Schoppa and F. Hussain, “Coherent structure dynamics in
    !   near-wall turbulence,” Fluid Dynamics Research, vol. 26,
    !   no. 2, pp. 119–139, Feb. 2000,
    !   doi: 10.1016/S0169-5983(99)00018-0.
    !
    !   E. de Villiers, “The Potential of Large Eddy Simulation for
    !   the Modelling of Wall Bounded Flows,” Degree of Doctor of
    !   Philosophy of the University of London and Diploma of
    !   Imperial College, Imperial College London, 2006.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: u(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: minx, miny, minz
        REAL(realk), INTENT(in) :: ReTau, h, lx, lz
        CHARACTER(len=4), INTENT(in) :: flowStateFlag

        ! Local parameters
        REAL(realk), PARAMETER :: duPl = 11.2_realk
        REAL(realk), PARAMETER :: epsPl = 0.072_realk
        REAL(realk), PARAMETER :: yMaxPl = 30.0_realk
        REAL(realk), PARAMETER :: devAmp = 0.2_realk
        INTEGER(intk), PARAMETER :: nStr = 4
        INTEGER(intk), PARAMETER :: nWav = 2

        ! Local variables
        INTEGER(intk) :: k, j, i
        REAL(realk) :: x, y, z, xMi, yMi, zMi, yPl
        REAL(realk) :: uTau, nu, alphaPl, betaPl, r

        nu = gmol1/rho1
        uTau = ReTau*nu/h

        alphaPl = 2.0_realk*pi*REAL(nWav, realk)/lx * nu/uTau
        betaPl = 2.0_realk*pi*REAL(nStr, realk)/lz * nu/uTau

        IF ( flowStateFlag == "TURB" ) THEN
            DO i = 2, ii-2
                yMi = miny
                DO j = 3, jj-2
                    y = MIN(yMi + 0.5_realk*ddy(j), h)
                    yPl = y*uTau/nu
                    zMi = minz
                    DO k = 3, kk-2
                        z = zMi + 0.5_realk*ddz(k)
                        CALL RANDOM_NUMBER(r)
                        u(k,j,i) = uTau*(reichardt(yPl) &
                            + (1.0_realk + devAmp*(2.0_realk*r - 1.0_realk)) &
                            * uPert(yPl, z*uTau/nu))
                        zMi = zMi + ddz(k)
                    END DO
                    yMi = yMi + ddy(j)
                END DO
            END DO

            xMi = minx
            DO i = 3, ii-2
                x = xMi + 0.5_realk*ddx(i)
                yMi = miny
                DO j = 3, jj-2
                    y = MIN(yMi + 0.5_realk*ddy(j), h)
                    DO k = 2, kk-2
                        CALL RANDOM_NUMBER(r)
                        w(k,j,i) = uTau*(1.0_realk + devAmp*(2.0_realk*r - 1.0_realk)) &
                            * wPert(x*uTau/nu, y*uTau/nu)
                    END DO
                    yMi = yMi + ddy(j)
                END DO
                xMi = xMi + ddx(i)
            END DO

        ELSE IF ( flowStateFlag == "LAMI" ) THEN
            yMi = miny
            DO j = 3, jj-2
                y = MIN(yMi + 0.5_realk*ddy(j), h)
                u(3:kk-2,j,2:ii-2) = uTau*parabola(y*uTau/nu, ReTau)
                yMi = yMi + ddy(j)
            END DO
        ELSE IF ( flowStateFlag == "NONE" ) THEN
            u = 0.0_realk
        ELSE
            CALL err_abort(mphInitErr, "unknown flow state.", __FILE__, __LINE__)
        END IF

    CONTAINS

        PURE REAL(realk) FUNCTION reichardt(yPl) RESULT(uPl)
        !------------------------------------------------------------
        !   Short: Mean velocity profile u+(y+) after Reichardt.
        !
        !   Source:
        !   H. Reichardt, “Vollständige Darstellung der turbulenten
        !   Geschwindigkeitsverteilung in glatten Leitungen,”
        !   Zeitschrift für Angewandte Mathematik und Mechanik,
        !   vol. 31, no. 7, pp. 208–219, Jan. 1951,
        !   doi: 10.1002/zamm.19510310704.
        !------------------------------------------------------------
            REAL(realk), INTENT(in) :: yPl
            REAL(realk), PARAMETER :: kappa = 0.41_realk, C = 7.8_realk
            uPl = LOG(1.0_realk + kappa*yPl)/kappa + &
                C*(1.0_realk - EXP(-yPl/11.0_realk) - yPl/11.0_realk*EXP(-yPl/3.0_realk))
        END FUNCTION reichardt

        PURE REAL(realk) FUNCTION parabola(yPl, ReTau) RESULT(uPl)
        !------------------------------------------------------------
        !   Short: Laminar open channel profile u+(y+).
        !------------------------------------------------------------
            REAL(realk), INTENT(in) :: yPl, ReTau
            uPl = yPl - 0.5_realk*yPl**2/ReTau
        END FUNCTION parabola

        PURE REAL(realk) FUNCTION uPert(yPl, zPl) RESULT(uPl)
        !------------------------------------------------------------
        !   Short: Near-wall streak perturbation u+(y+, z+).
        !
        !   Source:
        !   E. de Villiers, “The Potential of Large Eddy Simulation
        !   for the Modelling of Wall Bounded Flows,” Degree of
        !   Doctor of Philosophy of the University of London and
        !   Diploma of Imperial College, Imperial College London,
        !   2006.
        !------------------------------------------------------------
            REAL(realk), INTENT(in) :: yPl, zPl
            uPl = 0.5_realk*duPl*COS(betaPl*zPl)*(yPl/yMaxPl) &
                * EXP(0.5_realk*(1.0_realk - (yPl/yMaxPl)**2))
        END FUNCTION uPert

        PURE REAL(realk) FUNCTION wPert(xPl, yPl) RESULT(wPl)
        !------------------------------------------------------------
        !   Short: Spanwise perturbation w+(x+, y+).
        !
        !   Source:
        !   E. de Villiers, “The Potential of Large Eddy Simulation
        !   for the Modelling of Wall Bounded Flows,” Degree of
        !   Doctor of Philosophy of the University of London and
        !   Diploma of Imperial College, Imperial College London,
        !   2006.
        !------------------------------------------------------------
            REAL(realk), INTENT(in) :: xPl, yPl
            wPl = epsPl*SIN(alphaPl*xPl)*yPl*EXP(-0.5_realk*(yPl/yMaxPl)**2)
        END FUNCTION wPert

    END SUBROUTINE set_vel_cha_grid

    !================================================================

    SUBROUTINE frc_vel_fld(time)
    !----------------------------------------------------------------
    !   Short: Re-impose the prescribed velocity at a given time.
    !----------------------------------------------------------------

        ! Subroutine arguments
        REAL(realk), INTENT(in) :: time

        ! Local variables
        ! None

        CALL init_vel(time)

    END SUBROUTINE frc_vel_fld

    !================================================================

    SUBROUTINE init_pre()
    !----------------------------------------------------------------
    !   Short: Initialize the pressure field of the selected
    !          test case.
    !
    !   Only the open channel flow has an initial (hydrostatic)
    !   pressure; otherwise the pressure is unchanged.
    !----------------------------------------------------------------

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        SELECT CASE ( tstId )
        CASE ( tstOpCFl )
            CALL set_pre_cha(shape)
        CASE DEFAULT
            CONTINUE
        END SELECT

    END SUBROUTINE init_pre

    !================================================================

    SUBROUTINE set_pre_cha(s)
    !----------------------------------------------------------------
    !   Short: Set hydrostatic pressure for the open channel flow.
    !----------------------------------------------------------------

        ! Subroutine arguments
        TYPE(shape_t), INTENT(in) :: s

        ! Local variables
        INTEGER(intk) :: n, igrid
        INTEGER(intk) :: kk, jj, ii
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz
        REAL(realk), POINTER, CONTIGUOUS :: p(:,:,:), c(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: ddy(:)

        DO n = 1, nmygrids
            igrid = mygrids(n)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)
            CALL get_fieldptr(p, "P", igrid)
            CALL get_fieldptr(c, "C", igrid)
            CALL get_fieldptr(ddy, "DDY", igrid)

            CALL set_pre_cha_grid(kk, jj, ii, p, c, ddy, miny, s)
        END DO

    END SUBROUTINE set_pre_cha

    !================================================================

    SUBROUTINE set_pre_cha_grid(kk, jj, ii, p, c, ddy, miny, s)
    !----------------------------------------------------------------
    !   Short: Grid kernel of set_pre_cha.
    !
    !   p = rho(c)*g_y*(y - lvl) with the local mixture density, i.e.
    !   p = 0 at the interface and linear in both phases.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAl(realk), INTENT(inout) :: p(kk, jj, ii)
        REAL(realk), INTENT(in) :: c(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddy(jj)
        REAL(realk), INTENT(in) :: miny
        TYPE(shape_t), INTENT(in) :: s

        ! Local variables
        INTEGER(intk) :: j
        REAL(realk) :: y, yMi

        yMi = miny
        DO j = 3, jj-2
            y = yMi + 0.5_realk*ddy(j)
            p(3:kk-2, j, 3:ii-2) = (c(3:kk-2, j, 3:ii-2)*(rho1 - rho2) + rho2)* &
                grav(2)*(y - s%lvl)
            yMi = yMi + ddy(j)
        END DO

    END SUBROUTINE set_pre_cha_grid

    !================================================================

    PURE FUNCTION dist_func(x, y, z, s) RESULT(phi)
    !----------------------------------------------------------------
    !   Short: Signed distance to the selected shape.
    !
    !   Negative inside phase 1.
    !----------------------------------------------------------------
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
        CASE ( shpPlane )
            phi = y - s%lvl
        CASE ( shpSphere )
            phi = dist_func_sphere(x, y, z, s%xc, s%yc, s%zc, s%ra)
        CASE DEFAULT
            phi = 1.0_realk
        END SELECT
    END FUNCTION dist_func

    !================================================================

    PURE FUNCTION circf_func(s) RESULT(phi)
    !----------------------------------------------------------------
    !   Short: Exact interface length of the selected shape.
    !----------------------------------------------------------------
        TYPE(shape_t), INTENT(in) :: s
        REAL(realk) :: phi

        SELECT CASE ( s%shp )
        CASE ( shpCircle )
            phi = circf_func_circle(s%ra)
        CASE ( shpZalesak )
            phi = circf_func_circle(s%ra) - 2.0_realk*s%ra*ASIN(0.5_realk*s%slotW/s%ra) &
                + 2.0_realk*(s%slotH - s%ra + SQRT(s%ra**2 - 0.25_realk*s%slotW**2)) + s%slotW
        CASE DEFAULT
            phi = 1.0_realk
        END SELECT
    END FUNCTION circf_func

    !================================================================

    PURE FUNCTION area_func(s) RESULT(phi)
    !----------------------------------------------------------------
    !   Short: Exact area of the selected shape.
    !----------------------------------------------------------------
        TYPE(shape_t), INTENT(in) :: s
        REAL(realk) :: phi

        SELECT CASE ( s%shp )
        CASE ( shpSphere )
            phi = area_func_sphere(s%ra)
        CASE DEFAULT
            phi = 1.0_realk
        END SELECT
    END FUNCTION area_func

    !================================================================

    PURE FUNCTION dist_func_circle(x, y, xc, yc, r) RESULT(phi)
    !----------------------------------------------------------------
    !   Short: Signed distance to a circle.
    !----------------------------------------------------------------
        REAL(realk), INTENT(in) :: x, y, xc, yc, r
        REAL(realk) :: phi
        phi = SQRT((x - xc)**2 + (y - yc)**2) - r
    END FUNCTION dist_func_circle

    !================================================================

    PURE FUNCTION dist_func_box(x, y, xc, yc, hx, hy) RESULT(phi)
    !----------------------------------------------------------------
    !   Short: Signed distance to a rectangle (maximum norm).
    !----------------------------------------------------------------
        REAL(realk), INTENT(in) :: x, y, xc, yc, hx, hy
        REAL(realk) :: phi
        phi = MAX(ABS(x - xc) - hx, ABS(y - yc) - hy)
    END FUNCTION dist_func_box

    !================================================================

    PURE FUNCTION dist_func_sphere(x, y, z, xc, yc, zc, r) RESULT(phi)
    !----------------------------------------------------------------
    !   Short: Signed distance to a circle.
    !----------------------------------------------------------------
        REAL(realk), INTENT(in) :: x, y, z, xc, yc, zc, r
        REAL(realk) :: phi
        phi = SQRT((x - xc)**2 + (y - yc)**2 + (z - zc)**2) - r
    END FUNCTION dist_func_sphere

    !================================================================

    PURE FUNCTION circf_func_circle(ra) RESULT(phi)
    !----------------------------------------------------------------
    !   Short: Circumference of a circle.
    !----------------------------------------------------------------
        REAL(realk), INTENT(in) :: ra
        REAL(realk) :: phi
        phi = 2.0_realk*pi*ra
    END FUNCTION circf_func_circle

    !================================================================

    PURE FUNCTION area_func_sphere(ra) RESULT(phi)
    !----------------------------------------------------------------
    !   Short: Area of a circle.
    !----------------------------------------------------------------
        REAL(realk), INTENT(in) :: ra
        REAL(realk) :: phi
        phi = 4.0_realk*pi*ra**2
    END FUNCTION area_func_sphere

    !================================================================

    FUNCTION comp_abs_res(xMi, yMi, zMi, ddx, ddy, ddz, &
        normx, normy, normz, alpha, c, isIfc, s) RESULT(L1)
    !----------------------------------------------------------------
    !   Short: Cell-wise volume between exact and
    !          reconstructed interface.
    !
    !   Cells further than half a cell diagonal from the exact
    !   interface give (1 - c) or c times the cell volume. All other
    !   cells are sub-sampled (nSub x nSub points, 2D); the
    !   approximate phase follows from the PLIC plane in interface
    !   cells and from c > 0.5 otherwise.
    !----------------------------------------------------------------
        REAL(realk), INTENT(in) :: xMi, yMi, zMi, ddx, ddy, ddz
        REAL(realk), INTENT(in) :: normx, normy, normz, alpha, c
        REAL(realk), INTENT(in) :: isIfc
        TYPE(shape_t), INTENT(in) :: s
        REAL(realk) :: L1

        REAL(realk) :: vol, cAppr, dist, halfDiag, xs, ys, zs
        INTEGER(intk) :: is, js, inRibbon
        LOGICAL :: inExac, inAppr

        vol = ddx*ddy*ddz
        cAppr = c

        halfDiag = 0.5_realk*SQRT(ddx**2 + ddy**2 + ddz**2)
        dist = dist_func(xMi + 0.5_realk*ddx, yMi + 0.5_realk*ddy, zMi + 0.5_realk*ddz, s)

        IF ( dist <= -halfDiag ) THEN
            L1 = (1.0_realk - cAppr)*vol
        ELSE IF ( dist >= halfDiag ) THEN
            L1 = cAppr*vol
        ELSE
            inRibbon = 0
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
                    IF (inExac .NEQV. inAppr) inRibbon = inRibbon + 1
                END DO
            END DO
            L1 = vol*REAL(inRibbon,realk)/REAL(nSub,realk)**2
        END IF

    END FUNCTION comp_abs_res

END MODULE mph_test_mod