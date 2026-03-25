!====================================================================
!  Module: multiphase_vof_transport_mod
!
!  Responsibilities:
!     - Solves the incompressible volume fraction transport equation
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
    USE err_mod, ONLY: errr
    USE multiphase_plic_mod, ONLY: track_interface, compute_normal_vector, compute_alpha, compute_cFluxVol
    
    IMPLICIT NONE
    PRIVATE 

    PUBLIC :: multiphase_vof_transport

CONTAINS

    ! SUBROUTINE init_multiphase_vof_transport()

    !     continue

    ! END SUBROUTINE init_multiphase_vof_transport

    ! !================================================================

    ! SUBROUTINE finish_multiphase_vof_transport()

    !     continue

    ! END SUBROUTINE finish_multiphase_vof_transport

    ! !================================================================

    SUBROUTINE multiphase_vof_transport(vff_f, u_f, v_f, w_f, dtfu, itstep)
    !----------------------------------------------------------------
    !   What it does:
    !   The subroutine gets all nescessary input variables to perform
    !   the Volume-of-Fluid advection. It also manages the
    !   calculation on all grids.
    !----------------------------------------------------------------

        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: vff_f
        TYPE(field_t), INTENT(in) :: u_f
        TYPE(field_t), INTENT(in) :: v_f
        TYPE(field_t), INTENT(in) :: w_f
        REAL(realk), INTENT(in) :: dtfu
        INTEGER(intk), INTENT(in) :: itstep

        ! Local variables
        TYPE(field_t), POINTER :: ddx_f, ddy_f, ddz_f
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: vff, u, v, w
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

            CALL vff_f%get_ptr(vff, igrid)
            CALL u_f%get_ptr(u, igrid)
            CALL v_f%get_ptr(v, igrid)
            CALL w_f%get_ptr(w, igrid)

            CALL ddx_f%get_ptr(ddx, igrid)
            CALL ddy_f%get_ptr(ddy, igrid)
            CALL ddz_f%get_ptr(ddz, igrid)

            CALL multiphase_vof_transport_advection(kk, jj, ii, vff, u, v, w, & 
                ddx, ddy, ddz, dtfu, nfro, nbac, nrgt, nlft, nbot, ntop, itstep)
        END DO

    END SUBROUTINE multiphase_vof_transport

    !================================================================
    
    SUBROUTINE multiphase_vof_transport_advection(kk, jj, ii, vff, u, v, w, & 
        ddx, ddy, ddz, dtfu, nfro, nbac, nrgt, nlft, nbot, ntop, itstep)
    !----------------------------------------------------------------
    !   What it does:
    !   Coordinates the Volume-of-Fluid advection.
    !   1. Track interface cells
    !   2. Compute gradient unit normal vectors
    !   3. Compute interface distance
    !      => Interface is located
    !   4. Calculate fluxes
    !   5. Time integration of volume fraction field
    !   6. Clip volume fraction field to its bounds 
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: vff(kk, jj, ii)
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: dtfu
        INTEGER, INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop
        INTEGER(intk), INTENT(in) :: itstep

        ! Local variables
        LOGICAL :: isInterface(kk, jj, ii)
        REAL(realk) :: normx(kk, jj, ii), normy(kk, jj, ii), normz(kk, jj, ii)
        REAL(realk) :: alpha(kk, jj, ii)
        REAL(realk) :: fluxx(kk, jj, ii), fluxy(kk, jj, ii), fluxz(kk, jj, ii)
        REAL(realk), PARAMETER :: tol = 1.0E-15
        LOGICAL :: adv_x, adv_y, adv_z
        INTEGER(intk) :: permutation_index

        ! permutation_index only changes in a new time-step
        permutation_index = mod(itstep-1, 3)

        ! Select permutation of split advection
        SELECT CASE (permutation_index)
            CASE (0)
                adv_x = .TRUE.
                adv_y = .FALSE.
                adv_z = .FALSE.
            CASE (1)
                adv_x = .FALSE.
                adv_y = .TRUE.
                adv_z = .FALSE.
            CASE (2)
                adv_x = .FALSE.
                adv_y = .FALSE.
                adv_z = .TRUE.
        END SELECT

        ! Move vff in x-direction
        CALL track_interface(isInterface, kk, jj, ii, vff, tol)
        CALL compute_normal_vector(normx, normy, normz, kk, jj, ii, vff, ddx, ddy, ddz, tol)
        CALL compute_alpha(alpha, kk, jj, ii, vff, isInterface, ddx, ddy, ddz, normx, normy, normz, tol)
        CALL compute_fluxx(fluxx, kk, jj, ii, vff, isInterface, u, alpha, normx, normy, normz, ddy, ddz, tol, & 
            nfro, nbac, nrgt, nlft, nbot, ntop)
        CALL update_volume_fraction_field(kk, jj, ii, vff, fluxx, fluxy, fluxz, & 
            adv_x, adv_y, adv_z, tol, ddx, ddy, ddz, dtfu, nfro, nbac, nrgt, nlft, nbot, ntop)
        CALL clip_volume_fraction_field(kk, ii, jj, vff, tol)

        ! Move vff in y-direction
        CALL track_interface(isInterface, kk, jj, ii, vff, tol)
        CALL compute_normal_vector(normx, normy, normz, kk, jj, ii, vff, ddx, ddy, ddz, tol)
        CALL compute_alpha(alpha, kk, jj, ii, vff, isInterface, ddx, ddy, ddz, normx, normy, normz, tol)
        CALL compute_fluxy(fluxy, kk, jj, ii, vff, isInterface, v, alpha, normx, normy, normz, ddx, ddz, tol, & 
            nfro, nbac, nrgt, nlft, nbot, ntop)
        CALL update_volume_fraction_field(kk, jj, ii, vff, fluxx, fluxy, fluxz, & 
            adv_x, adv_y, adv_z, tol, ddx, ddy, ddz, dtfu, nfro, nbac, nrgt, nlft, nbot, ntop)
        CALL clip_volume_fraction_field(kk, ii, jj, vff, tol)

        ! Move vff in z-direction
        CALL track_interface(isInterface, kk, jj, ii, vff, tol)
        CALL compute_normal_vector(normx, normy, normz, kk, jj, ii, vff, ddx, ddy, ddz, tol)
        CALL compute_alpha(alpha, kk, jj, ii, vff, isInterface, ddx, ddy, ddz, normx, normy, normz, tol)
        CALL compute_fluxz(fluxz, kk, jj, ii, vff, isInterface, w, alpha, normx, normy, normz, ddx, ddy, tol, & 
            nfro, nbac, nrgt, nlft, nbot, ntop)
        CALL update_volume_fraction_field(kk, jj, ii, vff, fluxx, fluxy, fluxz, & 
            adv_x, adv_y, adv_z, tol, ddx, ddy, ddz, dtfu, nfro, nbac, nrgt, nlft, nbot, ntop)
        CALL clip_volume_fraction_field(kk, ii, jj, vff, tol)

    END SUBROUTINE multiphase_vof_transport_advection

    !================================================================

    SUBROUTINE compute_fluxx(fluxx, kk, jj, ii, vff, isInterface, u, alpha, normx, normy, normz, ddy, ddz, tol, & 
        nfro, nbac, nrgt, nlft, nbot, ntop)
    !----------------------------------------------------------------
    !   What it does:
    !   This subroutine calculates the flux of the volume fraction
    !   field vff in the x direction. The interface area is taken 
    !   into account. Therefore the flux has a unit of L^3/T. 
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(out) :: fluxx(kk, jj, ii)
        REAL(realk), INTENT(in) :: vff(kk, jj, ii)
        LOGICAL, INTENT(in) :: isInterface(kk, jj, ii)
        REAL(realk), INTENT(in) :: u(kk, jj, ii)
        REAL(realk), INTENT(in) :: alpha(kk, jj, ii)
        REAL(realk), INTENT(in) :: normx(kk, jj, ii), normy(kk, jj, ii), normz(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: tol
        INTEGER, INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: nbu, nfu, nrv, nbw, ntw, nlv
        REAL(realk) :: cFlux

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

        ! Calculate flux in x-direction
        DO i = 3-nfu, ii-3+nbu
            DO j = 3, jj-2
                DO k = 3, kk-2
                    IF ( u(k,j,i) > tol ) THEN
                        IF ( isInterface(k,j,i) ) THEN
                            ! Calculate flux for multiphase cell
                            cFlux = 0.0
                        ELSE
                            ! Calculate flux for singlephase cell
                            cFlux = vff(k,j,i) * abs( u(k,j,i) ) * ddy(j) * ddz(k)
                        END IF
                    ELSE IF ( u(k,j,i) < -tol ) THEN
                        IF ( isInterface(k,j,i) ) THEN
                            ! Calculate flux for multiphase cell
                            cFlux = 0.0
                        ELSE
                            ! Calculate flux for singlephase cell
                            cFlux = vff(k,j,i+1) * abs( u(k,j,i) ) * ddy(j) * ddz(k)
                        END IF
                    ELSE
                        cFlux = 0.0
                    END IF
                    fluxx(k,j,i) = sign( 1.0, u(k,j,i) ) * cFlux
                END DO
            END DO
        END DO

    END SUBROUTINE compute_fluxx

    !================================================================

    SUBROUTINE compute_fluxy(fluxy, kk, jj, ii, vff, isInterface, v, alpha, normx, normy, normz, ddx, ddz, tol, & 
        nfro, nbac, nrgt, nlft, nbot, ntop)
    !----------------------------------------------------------------
    !   What it does:
    !   This subroutine calculates the flux of the volume fraction
    !   field vff in the y direction. The interface area is taken 
    !   into account. Therefore the flux has a unit of L^3/T. 
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(out) :: fluxy(kk, jj, ii)
        REAL(realk), INTENT(in) :: vff(kk, jj, ii)
        LOGICAL, INTENT(in) :: isInterface(kk, jj, ii)
        REAL(realk), INTENT(in) :: v(kk, jj, ii)
        REAL(realk), INTENT(in) :: alpha(kk, jj, ii)
        REAL(realk), INTENT(in) :: normx(kk, jj, ii), normy(kk, jj, ii), normz(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddz(kk)
        REAL(realk), INTENT(in) :: tol
        INTEGER, INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: nbu, nfu, nrv, nbw, ntw, nlv
        REAL(realk) :: cFlux

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

        ! Calculate flux in y-direction
        DO i = 3, ii-2
            DO j = 3-nrv, jj-3+nlv
                DO k = 3, kk-2
                    IF ( v(k,j,i) > tol ) THEN
                        IF ( isInterface(k,j,i) ) THEN
                            ! Calculate flux for multiphase cell
                            cFlux = 0.0
                        ELSE
                            ! Calculate flux for singlephase cell
                            cFlux = vff(k,j,i) * abs( v(k,j,i) ) * ddx(i) * ddz(k)
                        END IF
                    ELSE IF ( v(k,j,i) < -tol ) THEN
                        IF ( isInterface(k,j,i) ) THEN
                            ! Calculate flux for multiphase cell
                            cFlux = 0.0
                        ELSE
                            ! Calculate flux for singlephase cell
                            cFlux = vff(k,j,i+1) * abs( v(k,j,i) ) * ddx(i) * ddz(k)
                        END IF
                    ELSE
                        cFlux = 0.0
                    END IF
                    fluxy(k,j,i) = sign( 1.0, v(k,j,i) ) * cFlux
                END DO
            END DO
        END DO

    END SUBROUTINE compute_fluxy

    !================================================================

    SUBROUTINE compute_fluxz(fluxz, kk, jj, ii, vff, isInterface, w, alpha, normx, normy, normz, ddx, ddy, tol, & 
        nfro, nbac, nrgt, nlft, nbot, ntop)
    !----------------------------------------------------------------
    !   What it does:
    !   This subroutine calculates the flux of the volume fraction
    !   field vff in the z direction. The interface area is taken 
    !   into account. Therefore the flux has a unit of L^3/T. 
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(out) :: fluxz(kk, jj, ii)
        REAL(realk), INTENT(in) :: vff(kk, jj, ii)
        LOGICAL, INTENT(in) :: isInterface(kk, jj, ii)
        REAL(realk), INTENT(in) :: w(kk, jj, ii)
        REAL(realk), INTENT(in) :: alpha(kk, jj, ii)
        REAL(realk), INTENT(in) :: normx(kk, jj, ii), normy(kk, jj, ii), normz(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj)
        REAL(realk), INTENT(in) :: tol
        INTEGER, INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: nbu, nfu, nrv, nbw, ntw, nlv
        REAL(realk) :: cFlux

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

        ! Calculate flux in z-direction
        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3-nbw, kk-3+ntw
                    IF ( w(k,j,i) > tol ) THEN
                        IF ( isInterface(k,j,i) ) THEN
                            ! Calculate flux for multiphase cell
                            cFlux = 0.0
                        ELSE
                            ! Calculate flux for singlephase cell
                            cFlux = vff(k,j,i) * abs( w(k,j,i) ) * ddx(i) * ddy(j)
                        END IF
                    ELSE IF ( w(k,j,i) < -tol ) THEN
                        IF ( isInterface(k,j,i) ) THEN
                            ! Calculate flux for multiphase cell
                            cFlux = 0.0
                        ELSE
                            ! Calculate flux for singlephase cell
                            cFlux = vff(k,j,i+1) * abs( w(k,j,i) ) * ddx(i) * ddy(j)
                        END IF
                    ELSE
                        cFlux = 0.0
                    END IF
                    fluxz(k,j,i) = sign( 1.0, w(k,j,i) ) * cFlux
                END DO
            END DO
        END DO

    END SUBROUTINE compute_fluxz

    !================================================================

    SUBROUTINE update_volume_fraction_field(kk, jj, ii, vff, fluxx, fluxy, fluxz, & 
            adv_x, adv_y, adv_z, tol, ddx, ddy, ddz, dtfu, nfro, nbac, nrgt, nlft, nbot, ntop)
    !----------------------------------------------------------------
    !   What it does:
    !   The subroutine performs the time integration of the volume 
    !   fraction field vff. In each Runge-Kutta step the subroutine 
    !   is called three times. Each time vff is updated taking into 
    !   account one spatial dimension. Each time-step the order of 
    !   the dimensional splitting is permuted by the calling 
    !   subroutine. The variables adv_(.) track, which dimension will 
    !   be integrated. At the moment this time integration performes 
    !   one Euler step with the length of the Runge-Kutta sub-step. 
    !   The previous Runge-Kutta stages are not considered.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: vff(kk, jj, ii)
        REAL(realk), INTENT(in) :: fluxx(kk, jj, ii), fluxy(kk, jj, ii), fluxz(kk, jj, ii)
        LOGICAL, INTENT(inout) :: adv_x, adv_y, adv_z
        REAL(realk), INTENT(in) :: tol
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
        
        ! Update volume fraction field with x fluxes
        IF ( adv_x ) THEN 

            DO i = 3-nfu, ii-3+nbu
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        vff(k,j,i) = vff(k,j,i) + dtfu * ( fluxx(k,j,i-1) - fluxx(k,j,i) ) / ( ddx(i) * ddy(j) * ddz(k) )
                    END DO 
                END DO 
            END DO

            adv_x = .false.
            adv_y = .true.
            adv_z = .false.

        ! Update volume fraction field with y fluxes
        ELSEIF ( adv_y ) THEN

            DO i = 3, ii-2
                DO j = 3-nrv, jj-3+nlv
                    DO k = 3, kk-2
                        vff(k,j,i) = vff(k,j,i) + dtfu * ( fluxy(k,j-1,i) - fluxy(k,j,i) ) / ( ddx(i) * ddy(j) * ddz(k) )
                    END DO 
                END DO 
            END DO

            adv_x = .false.
            adv_y = .false.
            adv_z = .true.

        ! Update volume fraction field with z fluxes
        ELSEIF ( adv_z ) THEN

            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3-nbw, kk-3+ntw
                        vff(k,j,i) = vff(k,j,i) + dtfu * ( fluxz(k-1,j,i) - fluxz(k,j,i) ) / ( ddx(i) * ddy(j) * ddz(k) )
                    END DO 
                END DO 
            END DO

            adv_x = .true.
            adv_y = .false.
            adv_z = .false.

        END IF

        ! ! DEBUG
        ! IF ( minval(vff) < -tol .OR. maxval(vff) > 1+tol ) THEN
        !     WRITE(*, *) "volume fraction field vff must be in bounds [0,1]. c_min = ", minval(vff), " c_max = ", maxval(vff)
        !     CALL errr(__FILE__, __LINE__)
        ! END IF

    END SUBROUTINE update_volume_fraction_field

    !================================================================

    SUBROUTINE clip_volume_fraction_field(kk, ii, jj, vff, tol)
    !----------------------------------------------------------------
    !   What it does:
    !   This subroutine ensures that the volume fraction field vff stays
    !   within its bounds of [0,1].
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: vff(kk, jj, ii)
        REAL(realk), INTENT(in) :: tol

        ! Local variables


        continue

    END SUBROUTINE clip_volume_fraction_field

END MODULE multiphase_vof_transport_mod
