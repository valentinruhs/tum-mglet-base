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

    USE precision_mod, ONLY: intk
    USE field_mod, ONLY: field_t
    USE grids_mod, ONLY: nmygrids, mygrids
    USE field_mod, ONLY: field_t
    USE grids_mod, ONLY: get_mgdims, get_mgbasb, get_bbox
    USE precision_mod, ONLY: intk, realk
    USE fields_mod, ONLY: get_field
    USE multiphasecore_mod, ONLY: test_multiphase
    USE connect2_mod, ONLY: connect
    USE grids_mod, ONLY: minlevel, maxlevel

    IMPLICIT NONE
    PRIVATE 

    REAL(realk), PARAMETER :: pi = 4.0_realk * atan(1.0_realk)
    REAL(realk), PROTECTED :: initErr
    REAL(realk), PROTECTED :: apprVol

    PUBLIC :: init_multiphase_io, finish_multiphase_io, update_velocity, apprVol

CONTAINS

    SUBROUTINE init_multiphase_io()

        ! Subroutine arguments
        ! None

        ! Local variables
        TYPE(field_t), POINTER :: vff_f
        TYPE(field_t), POINTER :: dx_f, dy_f, dz_f, ddx_f, ddy_f, ddz_f
        INTEGER(intk) :: n, igrid
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: i, j, k, di, dj, dk, r, ilevel
        INTEGER(intk) :: iSub, jSub, kSub
        REAL(realk), POINTER, CONTIGUOUS :: vff(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz
        REAL(realk) :: centx, centy, centz, rad, x, y, z, inside, a, b, e, maxTrans, theta
        REAL(realk) :: randx(10), randy(10), randt(10)
        REAL(realk) :: trueVol

        CALL get_field(vff_f, "VFF")
        CALL get_field(dx_f, "DX"); CALL get_field(dy_f, "DY"); CALL get_field(dz_f, "DZ")
        CALL get_field(ddx_f, "DDX"); CALL get_field(ddy_f, "DDY"); CALL get_field(ddz_f, "DDZ")

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

            CALL vff_f%get_ptr(vff, igrid)
            CALL dx_f%get_ptr(dx, igrid); CALL dy_f%get_ptr(dy, igrid); CALL dz_f%get_ptr(dz, igrid)
            CALL ddx_f%get_ptr(ddx, igrid); CALL ddy_f%get_ptr(ddy, igrid); CALL ddz_f%get_ptr(ddz, igrid)

            WRITE(*,'(A19,1X,I12)')   "INITIALIZE GRID    " , igrid

            SELECT CASE( test_multiphase )
            CASE ( 'SphTrF' ) ! Sphere Translation Fine
                !----------------------------------------------------
                centx = 0.5_realk ; centy = 0.5_realk ; centz = 0.5_realk ; rad = 0.06875_realk
                iSub = 32 ; jSub = 32 ; kSub = 32
                ! Outer loop over cells
                DO i = 3, ii-2 ; DO j = 3, jj-2 ; DO k = 3, kk-2
                    inside = 0.0_realk
                    ! Inner loop over (.)Sub for refinement
                    DO di = 0, iSub-1 ; DO dj = 0, jSub-1 ; DO dk = 0, kSub-1
                        x = minx + ( i-3 + (di + 0.5_realk)/iSub ) * dx(1) ! assume equidistance
                        y = miny + ( j-3 + (dj + 0.5_realk)/jSub ) * dy(1) ! assume equidistance
                        z = minz + ( k-3 + (dk + 0.5_realk)/kSub ) * dz(1) ! assume equidistance
                        IF ( (x - centx)**2.0_realk + &
                             (y - centy)**2.0_realk + &
                             (z - centz)**2.0_realk <= rad**2.0_realk ) THEN
                            inside = inside + 1.0_realk
                        ENDIF
                    ENDDO ; ENDDO ; ENDDO
                    vff(k,j,i) = inside / (iSub * jSub * kSub)
                ENDDO ; ENDDO ; ENDDO
                trueVol = 4.0_realk / 3.0_realk * pi * rad**3.0_realk
                apprVol = apprVol + SUM(vff) * ( dx(1) * dy(1) * dz(1) )
                CALL print_statistics(iSub, jSub, kSub, trueVol, apprVol)
                !----------------------------------------------------
            CASE ( 'SphTrC' ) ! Sphere Translation Coarse
                !----------------------------------------------------
                centx = 0.5_realk ; centy = 0.5_realk ; centz = 0.5_realk ; rad = 0.06875_realk
                iSub = 32 ; jSub = 32 ; kSub = 32
                ! Outer loop over cells
                DO i = 3, ii-2 ; DO j = 3, jj-2 ; DO k = 3, kk-2
                    inside = 0.0_realk
                    ! Inner loop over (.)Sub for refinement
                    DO di = 0, iSub-1 ; DO dj = 0, jSub-1 ; DO dk = 0, kSub-1
                        x = minx + ( i-3 + (di + 0.5_realk)/iSub ) * dx(1) ! assume equidistance
                        y = miny + ( j-3 + (dj + 0.5_realk)/jSub ) * dy(1) ! assume equidistance
                        z = minz + ( k-3 + (dk + 0.5_realk)/kSub ) * dz(1) ! assume equidistance
                        IF ( (x - centx)**2.0_realk + &
                             (y - centy)**2.0_realk + &
                             (z - centz)**2.0_realk <= rad**2.0_realk ) THEN
                            inside = inside + 1.0_realk
                        ENDIF
                    ENDDO ; ENDDO ; ENDDO
                    vff(k,j,i) = inside / (iSub * jSub * kSub)
                ENDDO ; ENDDO ; ENDDO
                trueVol = 4.0_realk / 3.0_realk * pi * rad**3.0_realk
                apprVol = apprVol + SUM(vff) * ( dx(1) * dy(1) * dz(1) )
                CALL print_statistics(iSub, jSub, kSub, trueVol, apprVol)
                !----------------------------------------------------
            CASE ( 'VorBoF' ) ! Vortex-in-a-Box Fine
                !----------------------------------------------------
                centx = 0.5_realk ; centy = 0.75_realk ; centz = 0.0_realk ; rad = 0.15_realk
                iSub = 512 ; jSub = 512 ; kSub = 1
                ! Outer loop over cells
                DO i = 3, ii-2 ; DO j = 3, jj-2 ; DO k = 3, kk-2
                    inside = 0.0_realk
                    ! Inner loop over (.)Sub for refinement
                    DO di = 0, iSub-1 ; DO dj = 0, jSub-1
                        x = minx + ( i-3 + (di + 0.5_realk)/iSub ) * dx(1) ! assume equidistance
                        y = miny + ( j-3 + (dj + 0.5_realk)/jSub ) * dy(1) ! assume equidistance
                        IF ( (x - centx)**2.0_realk + &
                             (y - centy)**2.0_realk <= rad**2.0_realk ) THEN
                            inside = inside + 1.0_realk
                        ENDIF
                    ENDDO ; ENDDO
                    vff(k,j,i) = inside / (iSub * jSub * kSub)
                ENDDO ; ENDDO ; ENDDO
                trueVol = pi * rad**2.0_realk * ( maxz - minz )
                apprVol = apprVol + SUM(vff) * ( dx(1) * dy(1) * dz(1) )
                CALL print_statistics(iSub, jSub, kSub, trueVol, apprVol)
                !----------------------------------------------------
            CASE ( 'VorBoC' ) ! Vortex-in-a-Box Coarse
                !----------------------------------------------------
                centx = 0.5_realk ; centy = 0.75_realk ; centz = 0.0_realk ; rad = 0.15_realk
                iSub = 512 ; jSub = 512 ; kSub = 1
                ! Outer loop over cells
                DO i = 3, ii-2 ; DO j = 3, jj-2 ; DO k = 3, kk-2
                    inside = 0.0_realk
                    ! Inner loop over (.)Sub for refinement
                    DO di = 0, iSub-1 ; DO dj = 0, jSub-1
                        x = minx + ( i-3 + (di + 0.5_realk)/iSub ) * dx(1) ! assume equidistance
                        y = miny + ( j-3 + (dj + 0.5_realk)/jSub ) * dy(1) ! assume equidistance
                        IF ( (x - centx)**2.0_realk + &
                             (y - centy)**2.0_realk <= rad**2.0_realk ) THEN
                            inside = inside + 1.0_realk
                        ENDIF
                    ENDDO ; ENDDO
                    vff(k,j,i) = inside / (iSub * jSub * kSub)
                ENDDO ; ENDDO ; ENDDO
                trueVol = pi * rad**2.0_realk * ( maxz - minz )
                apprVol = apprVol + SUM(vff) * ( dx(1) * dy(1) * dz(1) )
                CALL print_statistics(iSub, jSub, kSub, trueVol, apprVol)
                !----------------------------------------------------
            CASE ( 'CylAdF' ) ! Cylinder Advection Fine
                !----------------------------------------------------
                centx = 0.2_realk ; centy = 0.2_realk ; centz = 0.0_realk ; rad = 0.1_realk
                iSub = 512 ; jSub = 512 ; kSub = 1
                ! Outer loop over cells
                DO i = 3, ii-2 ; DO j = 3, jj-2 ; DO k = 3, kk-2
                    inside = 0.0_realk
                    ! Inner loop over (.)Sub for refinement
                    DO di = 0, iSub-1 ; DO dj = 0, jSub-1
                        x = minx + ( i-3 + (di + 0.5_realk)/iSub ) * dx(1) ! assume equidistance
                        y = miny + ( j-3 + (dj + 0.5_realk)/jSub ) * dy(1) ! assume equidistance
                        IF ( (x - centx)**2.0_realk + &
                             (y - centy)**2.0_realk <= rad**2.0_realk ) THEN
                            inside = inside + 1.0_realk
                        ENDIF
                    ENDDO ; ENDDO
                    vff(k,j,i) = inside / (iSub * jSub * kSub)
                ENDDO ; ENDDO ; ENDDO
                trueVol = pi * rad**2.0_realk * ( maxy - miny )
                apprVol = apprVol + SUM(vff) * ( dx(1) * dy(1) * dz(1) )
                CALL print_statistics(iSub, jSub, kSub, trueVol, apprVol)
                !----------------------------------------------------
            CASE ( 'CylAdC' ) ! Cylinder Advection Coarse
                !----------------------------------------------------
                centx = 0.2_realk ; centy = 0.2_realk ; centz = 0.0_realk ; rad = 0.1_realk
                iSub = 512 ; jSub = 512 ; kSub = 1
                ! Outer loop over cells
                DO i = 3, ii-2 ; DO j = 3, jj-2 ; DO k = 3, kk-2
                    inside = 0.0_realk
                    ! Inner loop over (.)Sub for refinement
                    DO di = 0, iSub-1 ; DO dj = 0, jSub-1
                        x = minx + ( i-3 + (di + 0.5_realk)/iSub ) * dx(1) ! assume equidistance
                        y = miny + ( j-3 + (dj + 0.5_realk)/jSub ) * dy(1) ! assume equidistance
                        IF ( (x - centx)**2.0_realk + &
                             (y - centy)**2.0_realk <= rad**2.0_realk ) THEN
                            inside = inside + 1.0_realk
                        ENDIF
                    ENDDO ; ENDDO
                    vff(k,j,i) = inside / (iSub * jSub * kSub)
                ENDDO ; ENDDO ; ENDDO
                trueVol = pi * rad**2.0_realk * ( maxz - minz )
                apprVol = apprVol + SUM(vff) * ( dx(1) * dy(1) * dz(1) )
                CALL print_statistics(iSub, jSub, kSub, trueVol, apprVol)
                !----------------------------------------------------
            CASE ( 'PlicEl' ) ! PLIC Ellipse
                !----------------------------------------------------
                DO r = 1, 1
                    randx = [0.2047, 0.8682, 0.8612, 0.9178,  -0.42, 0.7934, 0.9832, 0.0922, -0.348, 0.2323]
                    randy = [0.4138, 0.7814, 0.3582, -0.841, 0.9691, 0.8432, 0.2889, -0.527, 0.2153, -0.781]
                    randt = [0.7019, 0.5916, 0.6264, 0.9859, 0.4035, 0.1262, 0.3234, 0.1483, 0.8225, 0.4155]
                    maxTrans = 0.15_realk
                    theta = randt(r) * 2.0_realk * pi
                    centx = 0.5_realk !+ randx(r) * maxTrans
                    centy = 0.5_realk !+ randy(r) * maxTrans
                    a = 0.3464_realk ; b = 0.1414_realk ; e = SQRT(a**2.0_realk - b**2.0_realk)
                    iSub = 256 ; jSub = 256 ; kSub = 1
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
                    apprVol = apprVol + SUM(vff(3:kk-2,3:jj-2,3:ii-2)) * ( dx(1) * dy(1) * dz(1) )
                    WRITE(*,*) r
                    CALL print_statistics(iSub, jSub, kSub, trueVol, apprVol)
                ENDDO
                !----------------------------------------------------
            END SELECT
        ENDDO

        DO ilevel = minlevel, maxlevel
            CALL connect(ilevel, layers=2, s1=vff_f, corners=.TRUE.)
        ENDDO

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

    SUBROUTINE print_statistics(iSub, jSub, kSub, trueVol, apprVol)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: iSub, jSub, kSub
        REAL(realk), INTENT(in) :: trueVol, apprVol

        ! Local variables
        ! None
        
        initErr = ABS( trueVol - apprVol )

        WRITE(*,'(A32)')          "----------------------------------"
        WRITE(*,'(A19,1X,I12)')   "iSub:              " , iSub
        WRITE(*,'(A19,1X,I12)')   "jSub:              " , jSub
        WRITE(*,'(A19,1X,I12)')   "kSub:              " , kSub
        WRITE(*,'(A10,1X,E21.15)') "trueVol:  " , trueVol
        WRITE(*,'(A10,1X,E21.15)') "apprVol:  " , apprVol
        WRITE(*,'(A10,1X,E21.15)') "initErr:  " , initErr
        WRITE(*,'(A23)') ""

    END SUBROUTINE print_statistics

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
        INTEGER(intk) :: n, igrid!, i, j, k
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        INTEGER(intk) :: kk, jj, ii, k, j, i
        REAL(realk) :: magnitude, fac(6)
        REAL(realk), ALLOCATABLE :: psi(:,:,:)

        CALL get_field(dx_f, "DX")
        CALL get_field(dy_f, "DY")
        CALL get_field(dz_f, "DZ")

        CALL get_field(ddx_f, "DDX")
        CALL get_field(ddy_f, "DDY")
        CALL get_field(ddz_f, "DDZ")

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_mgdims(kk, jj, ii, igrid)

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

            IF ( test_multiphase == 'SphTrF' .OR. test_multiphase == 'SphTrC' ) THEN
                fac = [0.7686000, 0.5968000, 0.1035700, &
                       0.5514000, 0.2242010, 0.2512981]
                magnitude = 0.0125_realk

                ! In the first 1200 steps the sphere is translated in
                ! 6 random directions (see fac). In steps 1201-1400
                ! the sphere is translated back to [0.5, 0.5].
                IF ( itstep <= 1200 ) THEN 
                    u = cos(fac(floor((itstep-1)/200.0_realk) + 1)*2.0_realk*pi) * magnitude
                    v = sin(fac(floor((itstep-1)/200.0_realk) + 1)*2.0_realk*pi) * magnitude
                ELSE IF ( itstep <= 1400 ) THEN
                    u =   0.008793826_realk
                    v = - 0.008883616_realk
                ELSE 
                    u = 0.0_realk
                    v = 0.0_realk
                END IF
                w = 0.0_realk
            ELSE IF ( test_multiphase == 'VorBoF' .OR. test_multiphase == 'VorBoC' ) THEN
                IF (.NOT. ALLOCATED(psi)) ALLOCATE(psi(kk,jj,ii))
                DO i = 1, ii
                    DO j = 1, jj
                        DO k = 1, kk
                            psi(k,j,i) = - 1/pi * cos(pi * itstep * dt / 2.0_realk) * &
                                sin(pi*(- 2.0_realk * ddx(1) + i * ddx(1)))**2.0_realk * &
                                sin(pi*(- 2.0_realk * ddy(1) + j * ddy(1)))**2.0_realk
                        END DO
                    END DO 
                END DO
                DO i = 2, ii-1
                    DO j = 2, jj-1
                        DO k = 2, kk-1
                            u(k,j,i) = ( psi(k,j,i) - psi(k,j-1,i) ) / ddy(j)
                            v(k,j,i) = - ( psi(k,j,i) - psi(k,j,i-1) ) / ddx(i)
                            w(k,j,i) = 0.0_realk
                        END DO
                    END DO
                END DO   
            ELSE IF ( test_multiphase == 'CylAdF' .OR. test_multiphase == 'CylAdC' ) THEN
                IF ( itstep == 1 ) THEN
                    DO i = 1, ii
                        DO j = 1, jj
                            DO k = 1, kk
                                IF ( vff(k,j,i) > 0.0_realk ) THEN
                                    u(k,j,i) = 0.016_realk
                                    v(k,j,i) = 0.016_realk
                                    w(k,j,i) = 0.0_realk
                                END IF
                            END DO
                        END DO
                    END DO  
                END IF
            END IF

        END DO

    END SUBROUTINE

END MODULE multiphase_io_mod