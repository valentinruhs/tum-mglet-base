!====================================================================
!  Module: multiphase_vof_transport_mod
!
!  Responsibilities:
!     - Solves the volume fraction transport equation
!
!  Author:      Valentin Ruhs
!  Created:     2026-02
!  Last update: 2026-02
!
!====================================================================

MODULE multiphase_vof_transport_mod

    USE precision_mod, ONLY: intk, realk
    USE grids_mod, ONLY: nmygrids, mygrids
    USE field_mod, ONLY: field_t
    USE fields_mod, ONLY: get_field
    USE grids_mod, ONLY: get_mgdims, get_mgbasb
    
    IMPLICIT NONE
    PRIVATE 

    PUBLIC :: multiphase_vof_transport

CONTAINS

    SUBROUTINE init_multiphase_vof_transport()

        continue

    END SUBROUTINE init_multiphase_vof_transport

    !================================================================

    SUBROUTINE finish_multiphase_vof_transport()

        continue

    END SUBROUTINE finish_multiphase_vof_transport

    !================================================================

    SUBROUTINE multiphase_vof_transport(c_f, u_f, v_f, w_f, dtfu)

        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: c_f
        TYPE(field_t), INTENT(in) :: u_f
        TYPE(field_t), INTENT(in) :: v_f
        TYPE(field_t), INTENT(in) :: w_f
        REAL(realk), INTENT(in) :: dtfu

        ! Local variables
        TYPE(field_t), POINTER :: dx_f, dy_f, dz_f, ddx_f, ddy_f, ddz_f
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: c, u, v, w
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        INTEGER(intk) :: i, igrid
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: nfro, nbac, nrgt, nlft, nbot, ntop

        CALL get_field(ddx_f, "DDX")
        CALL get_field(ddy_f, "DDY")
        CALL get_field(ddz_f, "DDZ")

        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_mgbasb(nfro, nbac, nrgt, nlft, nbot, ntop, igrid)

            CALL c_f%get_ptr(c, igrid)
            CALL u_f%get_ptr(u, igrid)
            CALL v_f%get_ptr(v, igrid)
            CALL w_f%get_ptr(w, igrid)

            CALL ddx_f%get_ptr(ddx, igrid)
            CALL ddy_f%get_ptr(ddy, igrid)
            CALL ddz_f%get_ptr(ddz, igrid)

            CALL multiphase_vof_transport_advection(kk, jj, ii, c, u, v, w, & 
                ddx, ddy, ddz, dtfu, nfro, nbac, nrgt, nlft, nbot, ntop)
        END DO

    END SUBROUTINE multiphase_vof_transport

    !================================================================
    
    SUBROUTINE multiphase_vof_transport_advection(kk, jj, ii, c, u, v, w, & 
        ddx, ddy, ddz, dtfu, nfro, nbac, nrgt, nlft, nbot, ntop)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: c(kk, jj, ii)
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: dtfu
        INTEGER, INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop

        ! Local variables
        REAL(realk) :: fluxx(kk, jj, ii), fluxy(kk, jj, ii), fluxz(kk, jj, ii)
        REAL(realk), PARAMETER :: tol = 1.0E-15
        LOGICAL :: adv_x = .TRUE., adv_y = .FALSE., adv_z = .FALSE.

        ! Move c in x-direction--------------------------------------
        CALL compute_fluxx(fluxx, kk, jj, ii, c, u, ddy, ddz, tol, & 
            nfro, nbac, nrgt, nlft, nbot, ntop)
        CALL update_color_function(kk, jj, ii, c, fluxx, fluxy, fluxz, & 
            adv_x, adv_y, adv_z, ddx, ddy, ddz, dtfu, nfro, nbac, nrgt, nlft, nbot, ntop)

        ! Move c in y-direction--------------------------------------
        CALL compute_fluxy(fluxy, kk, jj, ii, c, v, ddx, ddz, tol, & 
            nfro, nbac, nrgt, nlft, nbot, ntop)
        CALL update_color_function(kk, jj, ii, c, fluxx, fluxy, fluxz, & 
            adv_x, adv_y, adv_z, ddx, ddy, ddz, dtfu, nfro, nbac, nrgt, nlft, nbot, ntop)

        ! Move c in z-direction--------------------------------------
        CALL compute_fluxz(fluxz, kk, jj, ii, c, w, ddx, ddy, tol, & 
            nfro, nbac, nrgt, nlft, nbot, ntop)
        CALL update_color_function(kk, jj, ii, c, fluxx, fluxy, fluxz, & 
            adv_x, adv_y, adv_z, ddx, ddy, ddz, dtfu, nfro, nbac, nrgt, nlft, nbot, ntop)

    END SUBROUTINE multiphase_vof_transport_advection

    !================================================================

    SUBROUTINE compute_fluxx(fluxx, kk, jj, ii, c, u, ddy, ddz, tol, & 
        nfro, nbac, nrgt, nlft, nbot, ntop)

        ! Subroutine arguments
        REAL(realk), INTENT(out) :: fluxx(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: c(kk, jj, ii)
        REAL(realk), INTENT(in) :: u(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: tol
        INTEGER, INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: nbu, nfu, nrv, nbw, ntw, nlv
        REAL(realk) :: c_flux

        nfu = 0
        nbu = 0
        nrv = 0
        nlv = 0
        nbw = 0
        ntw = 0

        ! CON = 7
        IF (nbac == 7) nbu = 1
        IF (nlft == 7) nlv = 1
        IF (ntop == 7) ntw = 1

        ! OP1 = 3
        IF (nfro == 3) nfu = 1
        IF (nbac == 3) nbu = 1
        IF (nrgt == 3) nrv = 1
        IF (nlft == 3) nlv = 1
        IF (nbot == 3) nbw = 1
        IF (ntop == 3) ntw = 1

        ! Calculate flux in x-direction------------------------------
        DO i = 3-nfu, ii-3+nbu
            DO j = 3, jj-2
                DO k = 3, kk-2
                    IF ( u(k,j,i) > tol ) THEN
                        c_flux = c(k,j,i) * abs( u(k,j,i) ) * ddy(j) * ddz(k)
                    ELSE IF ( u(k,j,i) < -tol ) THEN
                        c_flux = c(k,j,i+1) * abs( u(k,j,i) ) * ddy(j+1) * ddz(k+1)
                    ELSE
                        c_flux = 0.0
                    END IF
                    fluxx(k,j,i) = sign( 1.0, u(k,j,i) ) * c_flux
                END DO
            END DO
        END DO

    END SUBROUTINE compute_fluxx

    !================================================================

    SUBROUTINE compute_fluxy(fluxy, kk, jj, ii, c, v, ddx, ddz, tol, & 
        nfro, nbac, nrgt, nlft, nbot, ntop)

        ! Subroutine arguments
        REAL(realk), INTENT(out) :: fluxy(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: c(kk, jj, ii)
        REAL(realk), INTENT(in) :: v(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddz(kk)
        REAL(realk), INTENT(in) :: tol
        INTEGER, INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: nbu, nfu, nrv, nbw, ntw, nlv
        REAL(realk) :: c_flux

        nfu = 0
        nbu = 0
        nrv = 0
        nlv = 0
        nbw = 0
        ntw = 0

        ! CON = 7
        IF (nbac == 7) nbu = 1
        IF (nlft == 7) nlv = 1
        IF (ntop == 7) ntw = 1

        ! OP1 = 3
        IF (nfro == 3) nfu = 1
        IF (nbac == 3) nbu = 1
        IF (nrgt == 3) nrv = 1
        IF (nlft == 3) nlv = 1
        IF (nbot == 3) nbw = 1
        IF (ntop == 3) ntw = 1

        ! Calculate flux in y-direction------------------------------
        DO i = 3, ii-2
            DO j = 3-nrv, jj-3+nlv
                DO k = 3, kk-2
                    IF ( v(k,j,i) > tol ) THEN
                        c_flux = c(k,j,i) * abs( v(k,j,i) ) * ddx(i) * ddz(k)
                    ELSE IF ( v(k,j,i) < -tol ) THEN
                        c_flux = c(k,j,i+1) * abs( v(k,j,i) ) * ddx(i+1) * ddz(k+1)
                    ELSE
                        c_flux = 0.0
                    END IF
                    fluxy(k,j,i) = sign( 1.0, v(k,j,i) ) * c_flux
                END DO
            END DO
        END DO

    END SUBROUTINE compute_fluxy

    !================================================================

    SUBROUTINE compute_fluxz(fluxz, kk, jj, ii, c, w, ddx, ddy, tol, & 
        nfro, nbac, nrgt, nlft, nbot, ntop)

        ! Subroutine arguments
        REAL(realk), INTENT(out) :: fluxz(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: c(kk, jj, ii)
        REAL(realk), INTENT(in) :: w(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj)
        REAL(realk), INTENT(in) :: tol
        INTEGER, INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: nbu, nfu, nrv, nbw, ntw, nlv
        REAL(realk) :: c_flux

        nfu = 0
        nbu = 0
        nrv = 0
        nlv = 0
        nbw = 0
        ntw = 0

        ! CON = 7
        IF (nbac == 7) nbu = 1
        IF (nlft == 7) nlv = 1
        IF (ntop == 7) ntw = 1

        ! OP1 = 3
        IF (nfro == 3) nfu = 1
        IF (nbac == 3) nbu = 1
        IF (nrgt == 3) nrv = 1
        IF (nlft == 3) nlv = 1
        IF (nbot == 3) nbw = 1
        IF (ntop == 3) ntw = 1

        ! Calculate flux in z-direction------------------------------
        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3-nbw, kk-3+ntw
                    IF ( w(k,j,i) > tol ) THEN
                        c_flux = c(k,j,i) * abs( w(k,j,i) ) * ddx(i) * ddy(j)
                    ELSE IF ( w(k,j,i) < -tol ) THEN
                        c_flux = c(k,j,i+1) * abs( w(k,j,i) ) * ddx(i+1) * ddy(j+1)
                    ELSE
                        c_flux = 0.0
                    END IF
                    fluxz(k,j,i) = sign( 1.0, w(k,j,i) ) * c_flux
                END DO
            END DO
        END DO

    END SUBROUTINE compute_fluxz

    !================================================================

    SUBROUTINE update_color_function(kk, jj, ii, c, fluxx, fluxy, fluxz, & 
            adv_x, adv_y, adv_z, ddx, ddy, ddz, dtfu, nfro, nbac, nrgt, nlft, nbot, ntop)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: c(kk, jj, ii)
        REAL(realk), INTENT(in) :: fluxx(kk, jj, ii), fluxy(kk, jj, ii), fluxz(kk, jj, ii)
        LOGICAL, INTENT(inout) :: adv_x, adv_y, adv_z
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: dtfu
        INTEGER, INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: nbu, nfu, nrv, nbw, ntw, nlv

        nfu = 0
        nbu = 0
        nrv = 0
        nlv = 0
        nbw = 0
        ntw = 0

        ! CON = 7
        IF (nbac == 7) nbu = 1
        IF (nlft == 7) nlv = 1
        IF (ntop == 7) ntw = 1

        ! OP1 = 3
        IF (nfro == 3) nfu = 1
        IF (nbac == 3) nbu = 1
        IF (nrgt == 3) nrv = 1
        IF (nlft == 3) nlv = 1
        IF (nbot == 3) nbw = 1
        IF (ntop == 3) ntw = 1

        IF ( adv_x ) THEN 

            DO i = 3-nfu, ii-3+nbu
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        c(k,j,i) = c(k,j,i) + dtfu * ( fluxx(k,j,i-1) - fluxx(k,j,i) ) / ( ddx(i) * ddy(j) * ddz(k) )
                    END DO 
                END DO 
            END DO

            adv_x = .false.
            adv_y = .true.
            adv_z = .false.

        ELSEIF ( adv_y ) THEN

            DO i = 3, ii-2
                DO j = 3-nrv, jj-3+nlv
                    DO k = 3, kk-2
                        c(k,j,i) = c(k,j,i) + dtfu * ( fluxy(k,j-1,i) - fluxy(k,j,i) ) / ( ddx(i) * ddy(j) * ddz(k) )
                    END DO 
                END DO 
            END DO

            adv_x = .false.
            adv_y = .false.
            adv_z = .true.

        ELSEIF ( adv_z ) THEN

            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3-nbw, kk-3+ntw
                        c(k,j,i) = c(k,j,i) + dtfu * ( fluxz(k-1,j,i) - fluxz(k,j,i) ) / ( ddx(i) * ddy(j) * ddz(k) )
                    END DO 
                END DO 
            END DO

            adv_x = .true.
            adv_y = .false.
            adv_z = .false.

        END IF

    END SUBROUTINE update_color_function

END MODULE multiphase_vof_transport_mod
