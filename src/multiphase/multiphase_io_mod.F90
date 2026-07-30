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
    ! USE flowcore_mod, ONLY: uinf
    USE multiphasecore_mod, ONLY: gmol1, gmol2, rho1, rho2

    IMPLICIT NONE
    PRIVATE 

    REAL(realk), PARAMETER :: pi = 4.0_realk * ATAN(1.0_realk)
    REAL(realk), PROTECTED :: trueVol = 0.0_realk
    REAL(realk), PROTECTED :: initErr = 0.0_realk
    REAL(realk), PROTECTED :: initVol = 0.0_realk

    PUBLIC :: init_multiphase_io, finish_multiphase_io, update_velocity, trueVol, initErr, initVol, validate_velocity

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
        REAL(realk) :: vol

        CALL get_field(vff_f, "VFF")
        CALL get_field(dx_f, "DX"); CALL get_field(dy_f, "DY"); CALL get_field(dz_f, "DZ")
        CALL get_field(ddx_f, "DDX"); CALL get_field(ddy_f, "DDY"); CALL get_field(ddz_f, "DDZ")

        WRITE(*,'(A,I0,A)') "Initializing volume fraction field in ", nmygrids, " grids ..."
        WRITE(*,*) ""

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)

            CALL vff_f%get_ptr(vff, igrid)
            CALL dx_f%get_ptr(dx, igrid); CALL dy_f%get_ptr(dy, igrid); CALL dz_f%get_ptr(dz, igrid)
            CALL ddx_f%get_ptr(ddx, igrid); CALL ddy_f%get_ptr(ddy, igrid); CALL ddz_f%get_ptr(ddz, igrid)

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
                !----------------------------------------------------
            CASE ( 'SCylAc' ) ! Sudden Cylinder Accerleration
                !----------------------------------------------------
                centx = 0.2_realk ; centy = 0.2_realk ; centz = 0.0_realk ; rad = 0.1_realk
                iSub = 1024 ; jSub = 1024 ; kSub = 1
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
                !----------------------------------------------------
            CASE ( 'PlicEl' ) ! PLIC Ellipse
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
                    iSub = 1024 ; jSub = 1024 ; kSub = 1
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
                !----------------------------------------------------
                CASE ( 'StFstP' ) ! Stokes First Problem
                !----------------------------------------------------
                iSub = 1 ; jSub = 1024 ; kSub = 1
                ! Outer loop over cells
                DO i = 3, ii-2 ; DO j = 3, jj-2 ; DO k = 3, kk-2
                    inside = 0.0_realk
                    ! Inner loop over (.)Sub for refinement
                    DO dj = 0, jSub-1
                        y = miny + ( j-3 + (dj + 0.5_realk)/jSub ) * dy(1)
                        IF ( y <= 0.5_realk ) THEN
                            inside = inside + 1.0_realk
                        ENDIF
                    ENDDO
                    vff(k,j,i) = inside / (iSub * jSub * kSub)
                ENDDO ; ENDDO ; ENDDO
                trueVol = 0.5_realk * 1.0_realk * ( maxz - minz )
                !----------------------------------------------------
            END SELECT

            DO k = 3, kk-2
                DO j = 3, jj-2
                    DO i = 3, ii-2
                        vol = ddx(i) * ddy(j) * ddz(k)
                        initVol = initVol + vff(k,j,i) * vol
                    ENDDO
                ENDDO
            ENDDO
        ENDDO

        CALL MPI_Allreduce(MPI_IN_PLACE, initVol, 1, mglet_mpi_real, MPI_SUM, MPI_COMM_WORLD)

        initErr = trueVol - initVol

        IF ( myid == 0 ) THEN
            WRITE(*,'(A,ES14.6)') "Initial volume error is ", initErr, " (trueVol - initVol)"
            WRITE(*,*) ""
        ENDIF

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
        REAL(realk) :: magnitude, fac(6)
        REAL(realk), ALLOCATABLE :: psi(:,:,:)
        INTEGER(intk) :: halo

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
            ELSE IF ( test_multiphase == 'VorBoF' .OR. test_multiphase == 'VorBoC' ) THEN
                IF (.NOT. ALLOCATED(psi)) ALLOCATE(psi(kk,jj,ii))
                DO i = 1, ii
                    DO j = 1, jj
                        DO k = 1, kk
                            psi(k,j,i) = 1/pi * COS(pi * itstep * dt / 2.0_realk) * &
                                SIN(pi*(- 1.5_realk * ddx(1) + (i-1) * ddx(1)))**2.0_realk * &
                                SIN(pi*(- 1.5_realk * ddy(1) + (j-1) * ddy(1)))**2.0_realk
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
                    u = 0.016_realk
                    v = 0.016_realk
                    w = 0.0_realk
                END IF
            ELSE IF ( test_multiphase == 'SCylAc' ) THEN
                halo = 3
                IF ( itstep == 1 ) THEN
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
                END IF
            ELSE IF ( test_multiphase == 'StFstP' ) THEN
                IF ( itstep == 1 ) THEN
                    u = 0.0_realk
                    v = 0.0_realk
                    w = 0.0_realk
                ENDIF
            END IF

        END DO

    END SUBROUTINE

    !================================================================

    SUBROUTINE validate_velocity(u_f, v_f, w_f, itstep, dt)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        TYPE(field_t), INTENT(in) :: u_f
        TYPE(field_t), INTENT(in) :: v_f
        TYPE(field_t), INTENT(in) :: w_f
        INTEGER(intk), INTENT(in) :: itstep
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        INTEGER(intk) :: n, igrid, m
        TYPE(field_t), POINTER :: dx_f, dy_f, dz_f, ddx_f, ddy_f, ddz_f
        INTEGER(intk) :: kk, jj, ii, k, j, i
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: u, v, w
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        REAL(realk), ALLOCATABLE :: trueVel(:,:,:), trueVel1(:,:,:), trueVel2(:,:,:)
        REAL(realk) :: beta(100), betaHat(100), alpha, kappa, a, aHat, h, sigma, yHat, tHat, sum1, sum2

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

            CALL dx_f%get_ptr(dx, igrid)
            CALL dy_f%get_ptr(dy, igrid)
            CALL dz_f%get_ptr(dz, igrid)

            CALL ddx_f%get_ptr(ddx, igrid)
            CALL ddy_f%get_ptr(ddy, igrid)
            CALL ddz_f%get_ptr(ddz, igrid)

            IF ( test_multiphase == 'StFstP' ) THEN
                IF (.NOT. ALLOCATED(trueVel)) ALLOCATE(trueVel(kk,jj,ii))
                IF (.NOT. ALLOCATED(trueVel1)) ALLOCATE(trueVel1(kk,jj,ii))
                IF (.NOT. ALLOCATED(trueVel2)) ALLOCATE(trueVel2(kk,jj,ii))
                ! the betas are the positive roots of: cot(beta) + sigma * cot(kappa * beta * a)
                ! see: C.-O. Ng, “Starting flow in channels with boundary slip,” Meccanica, 
                !      vol. 52, no. 1–2, pp. 45–67, Jan. 2017, doi: 10.1007/s11012-016-0384-4.
                beta = [ 1.57079633,  4.71238898,  7.85398163, 10.99557429, 14.13716694, &
                        17.27875959, 20.42035225, 23.56194490, 26.70353756, 29.84513021, &
                        32.98672286, 36.12831552, 39.26990817, 42.41150082, 45.55309348, &
                        48.69468613, 51.83627878, 54.97787144, 58.11946409, 61.26105675, &
                        64.40264940, 67.54424205, 70.68583471, 73.82742736, 76.96902001, &
                        80.11061267, 83.25220532, 86.39379797, 89.53539063, 92.67698328, &
                        95.81857593, 98.96016859, 102.10176124, 105.24335390, 108.38494655, &
                        111.52653920, 114.66813186, 117.80972451, 120.95131716, 124.09290982, &
                        127.23450247, 130.37609512, 133.51768778, 136.65928043, 139.80087308, &
                        142.94246574, 146.08405839, 149.22565105, 152.36724370, 155.50883635, &
                        158.65042901, 161.79202166, 164.93361431, 168.07520697, 171.21679962, &
                        174.35839227, 177.49998493, 180.64157758, 183.78317024, 186.92476289, &
                        190.06635554, 193.20794820, 196.34954085, 199.49113350, 202.63272616, &
                        205.77431881, 208.91591146, 212.05750412, 215.19909677, 218.34068942, &
                        221.48228208, 224.62387473, 227.76546739, 230.90706004, 234.04865269, &
                        237.19024535, 240.33183800, 243.47343065, 246.61502331, 249.75661596, &
                        252.89820861, 256.03980127, 259.18139392, 262.32298657, 265.46457923, &
                        268.60617188, 271.74776454, 274.88935719, 278.03094984, 281.17254250, &
                        284.31413515, 287.45572780, 290.59732046, 293.73891311, 296.88050576, &
                        300.02209842, 303.16369107, 306.30528373, 309.44687638, 312.58846903 ]

                alpha = gmol2 / gmol1
                kappa = SQRT(gmol1 / rho1 * rho2 / gmol2) ; sigma = kappa * ( gmol2 / gmol1 )
                a = 0.5_realk ; h = 0.5_realk
                aHat = a / h ; betaHat = beta * h
                
                trueVel1 = 0.0_realk
                trueVel2 = 0.0_realk
                DO i = 3, ii-2
                    DO j = 3, jj-2
                        DO k = 3, kk-2
                            trueVel(k,j,i) = 0.0 !uinf(1) - uinf(1) * ERF(( ABS(ddy(1)) / 2 + (j-3) * ABS(ddy(1)) )/( SQRT(4.0_realk * gmol1 / rho1 * itstep * dt) ))
                            yHat = ( ABS(ddy(1)) / 2 + (j-3) * ABS(ddy(1)) ) / h
                            tHat = itstep * dt * ( ( gmol1 / rho1 ) / h**2.0_realk )

                            sum1 = 0.0_realk
                            sum2 = 0.0_realk
                            DO m = 1, 100
                                sum1 = sum1 + ( sin(kappa * betaHat(m) * aHat)**2.0_realk * sin(betaHat(m) * ( 1.0_realk + yHat )) ) / &
                                    ( betaHat(m) * ( sin(kappa * betaHat(m) * aHat)**2.0_realk + sigma * kappa * aHat * sin(betaHat(m))**2.0_realk ) ) * &
                                    EXP(-betaHat(m)**2.0_realk * tHat)
                                sum2 = sum2 + ( sin(betaHat(m)) * sin(kappa * betaHat(m) * aHat) * sin(kappa * betaHat(m) * (aHat - yHat)) ) / &
                                    ( betaHat(m) * ( sin(kappa * betaHat(m) * aHat)**2.0_realk + sigma * kappa * aHat * sin(betaHat(m))**2.0_realk ) ) * &
                                    EXP(-betaHat(m)**2.0_realk * tHat)
                            ENDDO

                            IF ( yHat <= 0.0 ) THEN
                                trueVel1(k,j,i) = ( alpha * aHat - yhat ) / ( alpha * aHat + 1.0_realk ) - 2.0_realk * sum1
                            ELSE
                                trueVel2(k,j,i) = ( alpha * aHat - yhat ) / ( alpha * aHat + 1.0_realk ) - 2.0_realk * sum2
                            ENDIF
                        ENDDO
                    ENDDO
                ENDDO
                DO j = 3, jj-2
                    ! WRITE(*,*) "ERROR at j = ", j, ": ", ABS(u(7,j,7) - trueVel(7,j,7))
                    ! WRITE(*,*) u(7,j,7) / uinf(1), trueVel(7,j,7) / uinf(1), trueVel1(7,j,7), trueVel2(7,j,7)
                ENDDO
            ELSE 
                return
            ENDIF

        ENDDO

    END SUBROUTINE validate_velocity

    !================================================================

    ! SUBROUTINE validate_ellipse()
    ! !----------------------------------------------------------------
    ! !   What it does:
    ! !   
    ! !----------------------------------------------------------------

    !     ! Subroutine arguments
    !     TYPE(field_t), INTENT(in) :: vff_f

    !     ! Local variables
    !     INTEGER(intk) :: l
    !     REAL(realk) :: centx, centy, x, y, a, b, e, maxTrans, theta, dlinSpac
    !     REAL(realk) :: randx(10), randy(10), randt(10)
    !     REAL(realk) :: linSpac(100000), ellx(100000), elly(100000)

    !     randx = [0.2047, 0.8682, 0.8612, 0.9178,  -0.42, 0.7934, 0.9832, 0.0922, -0.348, 0.2323]
    !     randy = [0.4138, 0.7814, 0.3582, -0.841, 0.9691, 0.8432, 0.2889, -0.527, 0.2153, -0.781]
    !     randt = [0.7019, 0.5916, 0.6264, 0.9859, 0.4035, 0.1262, 0.3234, 0.1483, 0.8225, 0.4155]
    !     maxTrans = 0.15_realk
    !     theta = randt(r) * 2.0_realk * pi
    !     centx = 0.5_realk + randx(r) * maxTrans
    !     centy = 0.5_realk + randy(r) * maxTrans
    !     a = 0.3464_realk ; b = 0.1414_realk ; e = SQRT(a**2.0_realk - b**2.0_realk)
    !     dlinSpac = 2.0_realk * pi / SIZE(linSpac)

    !     DO l = 1, 100000
    !         linSpac(l) = l * dlinSpac
    !     ENDDO

    !     ellx = centx + a * COS(linSpac) * COS(theta) - b * SIN(linSpac) * SIN(theta)
    !     elly = centy + a * COS(linSpac) * SIN(theta) + b * SIN(linSpac) * COS(theta)

    !     CALL get_field(normx_f, "NORMX")
    !     CALL get_field(normy_f, "NORMY")
    !     CALL get_field(normz_f, "NORMZ")
    !     CALL get_field(alpha_f, "ALPHA")

    !     DO n = 1, nmygrids
    !         igrid = mygrids(n)

    !         CALL get_mgdims(kk, jj, ii, igrid)

    !         CALL normx_f%get_ptr(normx, igrid)
    !         CALL normy_f%get_ptr(normy, igrid)
    !         CALL normz_f%get_ptr(normz, igrid)
    !         CALL alpha_f%get_ptr(alpha, igrid)

    !         DO i = 3, ii-2
    !             DO j = 3, jj-2
    !                 DO k = 3, kk-2
    !                     IF ( normx(k,j,i) >= 1.0E-12_realk .OR. &
    !                          normy(k,j,i) >= 1.0E-12_realk .OR. &
    !                          normz(k,j,i) >= 1.0E-12_realk ) THEN

    !                     ENDIF
    !                 ENDDO
    !             ENDDO
    !         ENDDO

    !     ENDDO

    ! END SUBROUTINE validate_ellipse

END MODULE multiphase_io_mod