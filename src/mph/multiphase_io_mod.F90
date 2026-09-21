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
        REAL(realk), ALLOCATABLE :: psi(:,:,:), xPl(:), yPl(:), zPl(:)

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

            ALLOCATE(xPl(ii))
            ALLOCATE(yPl(jj))
            ALLOCATE(zPl(kk))

            CALL get_top_right_corner(xPl, ddx, minx, ii)
            CALL get_top_right_corner(yPl, ddy, miny, jj)
            CALL get_top_right_corner(zPl, ddz, minz, kk)

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
                                SIN(pi*(xPl(i)+ddx(i)))**2 * &
                                SIN(pi*(yPl(j)+ddy(j)))**2
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

            DEALLOCATE(xPl)
            DEALLOCATE(yPl)
            DEALLOCATE(zPl)
        END DO

    END SUBROUTINE

END MODULE multiphase_io_mod