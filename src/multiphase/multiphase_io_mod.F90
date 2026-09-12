!====================================================================
!  Module: multiphase_io_mod
!
!  Responsibilities:
!     - Reads initial volume fraction field vff
!
!  Author:      Valentin Ruhs
!  Created:     2026-02
!  Last update: 2026-02
!
!====================================================================

MODULE multiphase_io_mod

    USE MPI_f08
    USE precision_mod, ONLY: intk, mglet_mpi_real
    USE field_mod, ONLY: field_t
    USE grids_mod, ONLY: nmygrids, mygrids
    USE field_mod, ONLY: field_t
    USE grids_mod, ONLY: get_mgdims, get_mgbasb, get_bbox
    USE comms_mod, ONLY: myid
    USE precision_mod, ONLY: intk, realk
    USE fields_mod, ONLY: get_field
    USE multiphasecore_mod, ONLY: test_multiphase
    USE connect2_mod, ONLY: connect
    USE grids_mod, ONLY: minlevel, maxlevel
    USE multiphasecore_mod, ONLY: gmol1, gmol2, rho1, rho2

    IMPLICIT NONE
    PRIVATE 

    REAL(realk), PARAMETER :: pi = 4.0_realk * ATAN(1.0_realk)
    REAL(realk), PROTECTED :: trueVol = 0.0_realk
    REAL(realk), PROTECTED :: initErr = 0.0_realk
    REAL(realk), PROTECTED :: initVol = 0.0_realk

    PUBLIC :: init_multiphase_io, finish_multiphase_io, update_velocity, trueVol, initErr, initVol

CONTAINS

    SUBROUTINE init_multiphase_io()

        ! Subroutine arguments
        ! None

        ! Local variables
        TYPE(field_t), POINTER :: u_f, v_f, w_f, vff_f
        TYPE(field_t), POINTER :: dx_f, dy_f, dz_f, ddx_f, ddy_f, ddz_f
        TYPE(field_t), POINTER :: grdMask_f
        INTEGER(intk) :: n, igrid
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: i, j, k, di, dj, dk, r, ilevel, halo
        INTEGER(intk) :: iSub, jSub, kSub
        REAL(realk), POINTER, CONTIGUOUS :: u(:,:,:), v(:,:,:), w(:,:,:), vff(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        REAL(realk), POINTER, CONTIGUOUS :: grdMask(:,:,:)
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz
        REAL(realk), ALLOCATABLE :: xR(:), yR(:), zR(:), psi(:,:,:)
        REAL(realk) :: centx, centy, centz, rad, dh, x, y, z, inside, a, b, e, maxTrans, theta
        REAL(realk) :: randx(10), randy(10), randt(10)
        REAL(realk) :: vol

        CALL get_field(u_f, "U"); CALL get_field(v_f, "V"); CALL get_field(w_f, "W"); CALL get_field(vff_f, "VFF")
        CALL get_field(dx_f, "DX"); CALL get_field(dy_f, "DY"); CALL get_field(dz_f, "DZ")
        CALL get_field(ddx_f, "DDX"); CALL get_field(ddy_f, "DDY"); CALL get_field(ddz_f, "DDZ")
        CALL get_field(grdMask_f, "GRDMASK")

        WRITE(*,'(A,I0,A)') "Initializing volume fraction field in ", nmygrids, " grids ..."
        WRITE(*,*) ""

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

            CALL u_f%get_ptr(u, igrid); CALL v_f%get_ptr(v, igrid); cALL w_f%get_ptr(w, igrid); CALL vff_f%get_ptr(vff, igrid)
            CALL dx_f%get_ptr(dx, igrid); CALL dy_f%get_ptr(dy, igrid); CALL dz_f%get_ptr(dz, igrid)
            CALL ddx_f%get_ptr(ddx, igrid); CALL ddy_f%get_ptr(ddy, igrid); CALL ddz_f%get_ptr(ddz, igrid)
            CALL grdMask_f%get_ptr(grdMask, igrid)

            ALLOCATE(xR(ii))
            ALLOCATE(yR(jj))
            ALLOCATE(zR(kk))

            CALL get_top_right_corner(xR, ddx, minx, ii)
            CALL get_top_right_corner(yR, ddy, miny, jj)
            CALL get_top_right_corner(zR, ddz, minz, kk)

            SELECT CASE( test_multiphase )
            CASE ( 'Sphere Translation' )
                !----------------------------------------------------
                centx = 0.5_realk
                centy = 0.5_realk
                centz = 0.5_realk
                rad = 0.06875_realk
                trueVol = 4.0_realk / 3.0_realk * pi * rad**3.0_realk
                iSub = 128
                jSub = 128
                DO i = 3, ii-2
                    DO j = 3, jj-2
                        DO k = 3, kk-2
                            inside = 0.0_realk
                            DO di = 0, iSub-1
                                x = xR(i) + (di + 0.5_realk)/iSub*ddx(i)
                                DO dj = 0, jSub-1
                                    y = yR(j) + (dj + 0.5_realk)/jSub*ddy(j)
                                    DO dk = 0, kSub-1
                                        z = zR(k) + (dk + 0.5_realk)/kSub*ddz(k)
                                        IF ((x - centx)**2 + (y - centy)**2 + (z - centz)**2 <= rad**2) THEN
                                            inside = inside + 1.0_realk
                                        ENDIF
                                    ENDDO
                                ENDDO
                            ENDDO
                            vff(3:kk-2,j,i) = inside / (iSub*jSub)
                        ENDDO
                    ENDDO
                ENDDO
                !----------------------------------------------------
            CASE ( 'Vortex in a Box' )
                !----------------------------------------------------
                centx = 0.5_realk
                centy = 0.75_realk
                centz = 0.0_realk
                rad = 0.15_realk
                trueVol = pi * rad**2.0_realk * (maxz - minz)
                iSub = 128
                jSub = 128
                DO i = 3, ii-2
                    DO j = 3, jj-2
                        inside = 0.0_realk
                        DO di = 0, iSub-1
                            x = xR(i) + (di + 0.5_realk)/iSub*ddx(i)
                            DO dj = 0, jSub-1
                                y = yR(j) + (dj + 0.5_realk)/jSub*ddy(j)
                                IF ((x - centx)**2 + (y - centy)**2 <= rad**2) THEN
                                    inside = inside + 1.0_realk
                                ENDIF
                            ENDDO
                        ENDDO
                        vff(3:kk-2,j,i) = inside / (iSub*jSub)
                    ENDDO
                ENDDO

                ALLOCATE(psi(kk,jj,ii))
                DO i = 1, ii
                    DO j = 1, jj
                        DO k = 1, kk
                            psi(k,j,i) = 1/pi * &
                                SIN(pi*(xR(i)+ddx(i)))**2 * &
                                SIN(pi*(yR(j)+ddy(j)))**2
                        END DO
                    END DO 
                END DO
                DO i = 2, ii-1
                    DO j = 2, jj-1
                        DO k = 2, kk-1
                            u(k,j,i) = (psi(k,j,i) - psi(k,j-1,i))/ddy(j)
                            v(k,j,i) = - (psi(k,j,i) - psi(k,j,i-1))/ddx(i)
                            w(k,j,i) = 0.0_realk
                        END DO
                    END DO
                END DO
                DEALLOCATE(psi)
                !----------------------------------------------------
            CASE ( 'Cylinder Advection' ) 
                !----------------------------------------------------
                centx = 0.2_realk
                centy = 0.2_realk
                centz = 0.0_realk
                rad = 0.1_realk
                trueVol = pi * rad**2.0_realk * (maxz - minz)
                iSub = 128
                jSub = 128
                DO i = 3, ii-2
                    DO j = 3, jj-2
                        inside = 0.0_realk
                        DO di = 0, iSub-1
                            x = xR(i) + (di + 0.5_realk)/iSub*ddx(i)
                            DO dj = 0, jSub-1
                                y = yR(j) + (dj + 0.5_realk)/jSub*ddy(j)
                                IF ((x - centx)**2 + (y - centy)**2 <= rad**2) THEN
                                    inside = inside + 1.0_realk
                                ENDIF
                            ENDDO
                        ENDDO
                        vff(3:kk-2,j,i) = inside / (iSub*jSub)
                    ENDDO
                ENDDO

                u = 0.016_realk
                v = 0.016_realk
                w = 0.0_realk
                !----------------------------------------------------
            CASE ( 'Sudden Cylinder Acceleration' )
                !----------------------------------------------------
                centx = 0.2_realk
                centy = 0.2_realk
                centz = 0.0_realk
                rad = 0.1_realk
                trueVol = pi * rad**2.0_realk * (maxz - minz)
                iSub = 128
                jSub = 128
                DO i = 3, ii-2
                    DO j = 3, jj-2
                        inside = 0.0_realk
                        DO di = 0, iSub-1
                            x = xR(i) + (di + 0.5_realk)/iSub*ddx(i)
                            DO dj = 0, jSub-1
                                y = yR(j) + (dj + 0.5_realk)/jSub*ddy(j)
                                IF ((x - centx)**2 + (y - centy)**2 <= rad**2) THEN
                                    inside = inside + 1.0_realk
                                ENDIF
                            ENDDO
                        ENDDO
                        vff(3:kk-2,j,i) = inside / (iSub*jSub)
                    ENDDO
                ENDDO

                halo = 1
                DO i = 3, ii-2
                    DO j = 3, jj-2
                        DO k = 3, kk-2
                            IF ( vff(k,j,i) > 0.0_realk .OR. &
                                 vff(k,j,i+halo) > 0.0_realk .OR. vff(k,j,i-halo) > 0.0_realk .OR. &
                                 vff(k,j+halo,i) > 0.0_realk .OR. vff(k,j-halo,i) > 0.0_realk .OR. &
                                 vff(k+halo,j,i) > 0.0_realk .OR. vff(k-halo,j,i) > 0.0_realk ) THEN
                                u(k,j,i) = 0.016_realk
                                v(k,j,i) = 0.016_realk
                                w(k,j,i) = 0.0_realk
                            ENDIF
                        ENDDO
                    ENDDO
                ENDDO
                !----------------------------------------------------
            CASE ( 'Open Channel Flow' )
                !----------------------------------------------------
                dh = 1.0_realk/64.0_realk
                trueVol = 2.0_realk * pi**2 * ( 1.0_realk - dh )
                kSub = 2
                DO k = 3, kk-2
                    inside = 0.0_realk
                    DO dk = 0, kSub-1
                        z = zR(k) + (dk + 0.5_realk)/kSub*ddz(k)
                        IF ( z <= - dh ) THEN
                            inside = inside + 1.0_realk
                        ENDIF
                    ENDDO
                    vff(k,3:jj-2,3:ii-2) = inside / (kSub)
                ENDDO

                u = 0.0_realk
                v = 0.0_realk
                w = 0.0_realk
                !----------------------------------------------------
            CASE ( 'PLIC Ellipse' )
                !----------------------------------------------------
                DO r = 1, 1
                    randx = [0.2047, 0.8682, 0.8612, 0.9178,  -0.42, 0.7934, 0.9832, 0.0922, -0.348, 0.2323]
                    randy = [0.4138, 0.7814, 0.3582, -0.841, 0.9691, 0.8432, 0.2889, -0.527, 0.2153, -0.781]
                    randt = [0.7019, 0.5916, 0.6264, 0.9859, 0.4035, 0.1262, 0.3234, 0.1483, 0.8225, 0.4155]
                    maxTrans = 0.15_realk
                    theta = randt(r) * 2.0_realk * pi
                    centx = 0.5_realk + randx(r) * maxTrans
                    centy = 0.5_realk + randy(r) * maxTrans
                    a = 0.3464_realk ; b = 0.1414_realk ; e = SQRT(a**2.0_realk - b**2.0_realk)
                    iSub = 128 ; jSub = 128 ; kSub = 1
                    ! Outer loop over cells
                    DO i = 3, ii-2 ; DO j = 3, jj-2 ; DO k = 3, kk-2
                        inside = 0.0_realk
                        ! Inner loop over (.)Sub for refinement
                        DO di = 0, iSub-1 ; DO dj = 0, jSub-1
                            x = minx + ( i-3 + (di + 0.5_realk)/iSub ) * dx(1) ! assume equidistance
                            y = miny + ( j-3 + (dj + 0.5_realk)/jSub ) * dy(1) ! assume equidistance
                            IF ( SQRT( ((x - centx) - e*COS(theta))**2.0_realk + ((y - centy) - e*SIN(theta))**2.0_realk ) + &
                                 SQRT( ((x - centx) + e*COS(theta))**2.0_realk + ((y - centy) + e*SIN(theta))**2.0_realk ) <= 2*a ) THEN
                                inside = inside + 1.0_realk
                            ENDIF
                        ENDDO ; ENDDO
                        vff(k,j,i) = inside / (iSub * jSub * kSub)
                    ENDDO ; ENDDO ; ENDDO
                    trueVol = pi * a * b * ( maxz - minz )
                ENDDO
            END SELECT

            DO k = 3, kk-2
                DO j = 3, jj-2
                    DO i = 3, ii-2
                        vol = grdMask(k,j,i) * ddx(i) * ddy(j) * ddz(k)
                        initVol = initVol + vff(k,j,i) * vol
                    ENDDO
                ENDDO
            ENDDO

            DEALLOCATE(xR)
            DEALLOCATE(yR)
            DEALLOCATE(zR)
        ENDDO

        CALL MPI_Allreduce(MPI_IN_PLACE, initVol, 1, mglet_mpi_real, MPI_SUM, MPI_COMM_WORLD)

        initErr = trueVol - initVol

        IF ( myid == 0 ) THEN
            WRITE(*,'(A,ES24.16,A)') "True volume is ", trueVol, " (trueVol)"
            WRITE(*,'(A,ES24.16,A)') "Initial volume is ", initVol, " (initVol)"
            WRITE(*,'(A,ES24.16,A)') "Initial volume error is ", initErr, " (trueVol - initVol)"
            WRITE(*,*) ""
        ENDIF

        DO ilevel = minlevel, maxlevel
            CALL connect(ilevel, layers=2, s1=vff_f, corners=.TRUE.)
        ENDDO

    CONTAINS

        PURE REAL(realk) FUNCTION reichardt(yPl) RESULT(uPl)
            REAL(realk), INTENT(in) :: yPl
            REAL(realk), PARAMETER :: C=7.8_realk, kappa=0.41
            uPl = (LOG(1.0_realk + kappa*yPl)/kappa +  & 
                C*(1.0_realk - EXP(-yPl/11.0_realk) - yPl/11.0_realk*EXP(-yPl/3.0_realk)))
        END FUNCTION reichardt

    END SUBROUTINE init_multiphase_io

    !================================================================

    SUBROUTINE finish_multiphase_io()

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        continue
    END SUBROUTINE finish_multiphase_io

    !================================================================

    SUBROUTINE update_velocity(u_f, v_f, w_f, vff_f, itstep, dt)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        TYPE(field_t), INTENT(in) :: u_f
        TYPE(field_t), INTENT(in) :: v_f
        TYPE(field_t), INTENT(in) :: w_f
        TYPE(field_t), INTENT(in) :: vff_f
        INTEGER(intk), INTENT(in) :: itstep
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        TYPE(field_t), POINTER :: dx_f, dy_f, dz_f, ddx_f, ddy_f, ddz_f
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: u, v, w, vff
        INTEGER(intk) :: n, igrid
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        INTEGER(intk) :: kk, jj, ii, k, j, i
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz
        REAL(realk) :: magnitude, fac(6)
        REAL(realk), ALLOCATABLE :: psi(:,:,:), xR(:), yR(:), zR(:)

        CALL get_field(dx_f, "DX")
        CALL get_field(dy_f, "DY")
        CALL get_field(dz_f, "DZ")

        CALL get_field(ddx_f, "DDX")
        CALL get_field(ddy_f, "DDY")
        CALL get_field(ddz_f, "DDZ")

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

            CALL u_f%get_ptr(u, igrid)
            CALL v_f%get_ptr(v, igrid)
            CALL w_f%get_ptr(w, igrid)
            CALL vff_f%get_ptr(vff, igrid)

            CALL dx_f%get_ptr(dx, igrid)
            CALL dy_f%get_ptr(dy, igrid)
            CALL dz_f%get_ptr(dz, igrid)

            CALL ddx_f%get_ptr(ddx, igrid)
            CALL ddy_f%get_ptr(ddy, igrid)
            CALL ddz_f%get_ptr(ddz, igrid)

            ALLOCATE(xR(ii))
            ALLOCATE(yR(jj))
            ALLOCATE(zR(kk))

            CALL get_top_right_corner(xR, ddx, minx, ii)
            CALL get_top_right_corner(yR, ddy, miny, jj)
            CALL get_top_right_corner(zR, ddz, minz, kk)

            SELECT CASE( test_multiphase )
            CASE ( 'Sphere Translation' )
                fac = [0.7686000, 0.5968000, 0.1035700, &
                       0.5514000, 0.2242010, 0.2512981]
                magnitude = 0.0125_realk

                ! In the first 1200 steps the sphere is translated in
                ! 6 random directions (see fac). In steps 1201-1400
                ! the sphere is translated back to [0.5, 0.5].
                IF ( itstep <= 1200 ) THEN 
                    u = COS(fac(floor((itstep-1)/200.0_realk) + 1)*2.0_realk*pi) * magnitude
                    v = SIN(fac(floor((itstep-1)/200.0_realk) + 1)*2.0_realk*pi) * magnitude
                ELSE IF ( itstep <= 1400 ) THEN
                    u =   0.008793826_realk
                    v = - 0.008883616_realk
                ELSE 
                    u = 0.0_realk
                    v = 0.0_realk
                END IF
                w = 0.0_realk
            CASE ( 'Vortex in a Box' )
                ALLOCATE(psi(kk,jj,ii))
                DO i = 1, ii
                    DO j = 1, jj
                        DO k = 1, kk
                            psi(k,j,i) = 1/pi * COS(pi*dt*itstep/2.0_realk) * &
                                SIN(pi*(xR(i)+ddx(i)))**2 * &
                                SIN(pi*(yR(j)+ddy(j)))**2
                        END DO
                    END DO 
                END DO
                DO i = 2, ii-1
                    DO j = 2, jj-1
                        DO k = 2, kk-1
                            u(k,j,i) = (psi(k,j,i) - psi(k,j-1,i))/ddy(j)
                            v(k,j,i) = - (psi(k,j,i) - psi(k,j,i-1))/ddx(i)
                            w(k,j,i) = 0.0_realk
                        END DO
                    END DO
                END DO
                DEALLOCATE(psi)
            END SELECT

            DEALLOCATE(xR)
            DEALLOCATE(yR)
            DEALLOCATE(zR)
        END DO

    END SUBROUTINE

    !================================================================

    SUBROUTINE get_top_right_corner(c, ddn, cmin, nmax)
    !----------------------------------------------------------------
    !   What it does:
    !    
    !----------------------------------------------------------------

        ! Subroutine arguments
        REAL(realk), INTENT(out) :: c(:)
        REAL(realk), INTENT(in)  :: ddn(:), cmin
        INTEGER, INTENT(in)      :: nmax

        ! Local variables
        INTEGER :: n

        c(3) = cmin
        DO n = 4, nmax
            c(n) = c(n-1) + ddn(n-1)
        END DO

        DO n = 2, 1, -1
            c(n) = c(n+1) - ddn(n)
        END DO
    END SUBROUTINE

END MODULE multiphase_io_mod