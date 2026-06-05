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

    IMPLICIT NONE
    PRIVATE 

    PUBLIC :: init_multiphase_io, finish_multiphase_io, update_velocity

CONTAINS

    SUBROUTINE init_multiphase_io()

        ! Subroutine arguments
        ! None

        ! Local variables
        TYPE(field_t), POINTER :: vff_f
        TYPE(field_t), POINTER :: dx_f, dy_f, dz_f, ddx_f, ddy_f, ddz_f
        INTEGER(intk) :: n, igrid
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: i, j, k, di, dj, dk
        INTEGER(intk) :: iSub, jSub, kSub
        INTEGER(intk) :: nxDom, nyDom, nzDom
        REAL(realk), POINTER, CONTIGUOUS :: vff(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz
        REAL(realk) :: centx, centy, centz, rad, x, y, z, inside

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

            SELECT CASE( test_multiphase )
            CASE ( 'SphTra' )

                centx = 0.2_realk ; centy = 0.2_realk ; rad = 0.1_realk
                iSub = 128 ; jSub = 128 ; kSub = 128
                nxDom = 16 ; nyDom = 16 ; nzDom = 16

                ! Outer loop over cells
                DO i = 1, ii ; DO j = 1, jj ; DO k = 1, kk

                    inside = 0

                    ! Inner loop over (.)Sub for refinement
                    DO di = 1, iSub ; DO dj = 1, jSub ; DO dk = 1, kSub

                        x = ( i + (di + 0.5_realk)/iSub ) / nxDom
                        y = ( j + (dj + 0.5_realk)/jSub ) / nyDom
                        z = ( k + (dk + 0.5_realk)/kSub ) / nzDom

                        IF ( (x - centx)**2_realk + (y - centy)**2_realk + (z - centz)**2_realk <= rad**2_realk ) THEN
                            inside = inside + 1
                        ENDIF
                        
                    ENDDO ; ENDDO ; ENDDO

                    vff(k,j,i) = inside / (iSub * jSub * kSub)

                ENDDO ; ENDDO ; ENDDO

            CASE ( 'VorBoF' )
                continue
            CASE ( 'VorBoC' )
                                
                centx = 0.5_realk ; centy = 0.75_realk ; rad = 0.15_realk
                iSub = 256 ; jSub = 256
                nxDom = 32 ; nyDom = 32
                
                ! Outer loop over cells
                DO i = 3, ii-2 ; DO j = 3, jj-2 ; DO k = 3, kk-2

                    inside = 0.0_realk

                    ! Inner loop over (.)Sub for refinement
                    DO di = 0, iSub-1 ; DO dj = 0, jSub-1

                        x = minx + ( i-3 + (di + 0.5_realk)/iSub ) / nxDom
                        y = miny + ( j-3 + (dj + 0.5_realk)/jSub ) / nyDom

                        IF ( (x - centx)**2.0_realk + (y - centy)**2.0_realk <= rad**2.0_realk ) THEN
                            inside = inside + 1.0_realk
                        ENDIF
                        
                    ENDDO ; ENDDO
                    
                    vff(k,j,i) = inside / (iSub * jSub)

                ENDDO ; ENDDO ; ENDDO

            CASE ( 'CylAdF' )
                continue
            CASE ( 'CylAdC' )
                continue
            END SELECT
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
        INTEGER(intk) :: n, igrid!, i, j, k
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        INTEGER(intk) :: kk, jj, ii, k, j, i
        REAL(realk) :: magnitude, fac(6)
        REAL(realk), PARAMETER :: pi = 4.0_realk * atan(1.0_realk)
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

            IF ( test_multiphase == 'SphTra' ) THEN
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
                    DO i = 2, ii-1
                        DO j = 2, jj-1
                            DO k = 2, kk-1
                                IF ( vff(k,j,i) >= 0.0_realk ) THEN
                                    u(k,j,i) = 0.016_realk
                                    v(k,j,i) = 0.0_realk
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