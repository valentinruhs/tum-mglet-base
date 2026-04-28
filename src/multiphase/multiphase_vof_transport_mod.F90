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
    USE multiphase_plic_mod, ONLY: track_interface, compute_normal_vector, compute_alpha, compute_cell_proportion, compute_iStag_vff, compute_jStag_vff, compute_kStag_vff, interface_reconstruction_wrapper
    USE rungekutta_mod, ONLY: rk_2n_t
    
    IMPLICIT NONE
    PRIVATE 

    PUBLIC :: init_multiphase_vof_transport, finish_multiphase_vof_transport, field_flux_wrapper, compute_density_flux, get_advection_sequence, compression_term_wrapper, update_field, clip_volume_fraction_field

    INTERFACE field_flux_wrapper
        MODULE PROCEDURE field_flux_wrapper_pres
        MODULE PROCEDURE field_flux_wrapper_stag
    END INTERFACE

    INTERFACE compute_normal_strain_rates
        MODULE PROCEDURE compute_normal_strain_rates_pres
        MODULE PROCEDURE compute_normal_strain_rates_stag
    END INTERFACE

    INTERFACE compression_term_wrapper
        MODULE PROCEDURE compression_term_wrapper_pres
        MODULE PROCEDURE compression_term_wrapper_stag
    END INTERFACE

    INTERFACE compute_non_directional_compression_coeffiecient
        MODULE PROCEDURE compute_non_directional_compression_coeffiecient_pres
        MODULE PROCEDURE compute_non_directional_compression_coeffiecient_stag
    END INTERFACE

CONTAINS

    SUBROUTINE init_multiphase_vof_transport()

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        continue
    END SUBROUTINE init_multiphase_vof_transport

    !================================================================

    SUBROUTINE finish_multiphase_vof_transport()

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None
        
        continue
    END SUBROUTINE finish_multiphase_vof_transport

    !================================================================

    SUBROUTINE compute_normal_strain_rates_pres(kk, jj, ii, strainRateX, strainRateY, strainRateZ, u, v, w, ddx, ddy, ddz)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(out) :: strainRateX(kk, jj, ii), strainRateY(kk, jj, ii), strainRateZ(kk, jj, ii)
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)

        ! Local variables
        INTEGER(intk) :: k, j, i

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    strainRateX(k,j,i) = ( u(k,j,i) - u(k,j,i-1) ) / ddx(i)
                    strainRateY(k,j,i) = ( v(k,j,i) - v(k,j-1,i) ) / ddy(j)
                    strainRateZ(k,j,i) = ( w(k,j,i) - w(k-1,j,i) ) / ddz(k)
                END DO
            END DO
        END DO

    END SUBROUTINE compute_normal_strain_rates_pres

    !================================================================

    SUBROUTINE compute_normal_strain_rates_stag(kk, jj, ii, component, strainRateX, strainRateY, strainRateZ, u, v, w, dx, dy, dz, ddx, ddy, ddz)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: component
        REAL(realk), INTENT(out) :: strainRateX(kk, jj, ii, 3), strainRateY(kk, jj, ii, 3), strainRateZ(kk, jj, ii, 3)
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)

        ! Local variables
        INTEGER(intk) :: k, j, i
        REAL(realk) :: uE, uW, vN, vS, wT, wB
        REAL(realk) :: iStag, jStag, kStag

        IF ( component == 1 ) THEN
            iStag = 1.0
            jStag = 0.0
            kStag = 0.0
        ELSE IF ( component == 2 ) THEN
            iStag = 0.0
            jStag = 1.0
            kStag = 0.0
        ELSE IF ( component == 3 ) THEN
            iStag = 0.0
            jStag = 0.0
            kStag = 1.0
        END IF

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2

                    uE = 0.5 * ( iStag * ( u(k,j,i) + u(k,j,i+1) ) + jStag * ( u(k,j,i) + u(k,j+1,i) ) + kStag * ( u(k,j,i) + u(k+1,j,i) ) ) 
                    uW = 0.5 * ( iStag * ( u(k,j,i-1) + u(k,j,i) ) + jStag * ( u(k,j,i-1) + u(k,j+1,i-1) ) + kStag * ( u(k,j,i-1) + u(k+1,j,i-1) ) )
                    vN = 0.5 * ( iStag * ( v(k,j,i) + v(k,j,i+1) ) + jStag * ( v(k,j,i) + v(k,j+1,i) ) + kStag * ( v(k,j,i) + v(k+1,j,i) ) )
                    vS = 0.5 * ( iStag * ( v(k,j-1,i) + v(k,j-1,i+1) ) + jStag * ( v(k,j-1,i) + v(k,j,i) ) + kStag * ( v(k+1,j-1,i) + v(k,j-1,i) ) )
                    wT = 0.5 * ( iStag * ( w(k,j,i) + w(k,j,i+1) ) + jStag * ( w(k,j,i) + w(k,j+1,i) ) + kStag * ( w(k,j,i) + w(k+1,j,i) ) ) 
                    wB = 0.5 * ( iStag * ( w(k-1,j,i) + w(k-1,j,i+1) ) + jStag * ( w(k-1,j,i) + w(k-1,j+1,i) ) + kStag * ( w(k-1,j,i) + w(k,j,i) ) ) 

                    strainRateX(k,j,i,component) = ( uE - uW ) / ( iStag * dx(i) + jStag * ddx(i) + kStag * ddx(i) )
                    strainRateY(k,j,i,component) = ( vN - vS ) / ( iStag * ddy(j) + jStag * dy(j) + kStag * ddy(j) )
                    strainRateZ(k,j,i,component) = ( wT - wB ) / ( iStag * ddz(k) + jStag * ddz(k) + kStag * dz(k) )

                END DO
            END DO
        END DO

    END SUBROUTINE compute_normal_strain_rates_stag
    
    !================================================================

    SUBROUTINE compute_non_directional_compression_coeffiecient_pres(kk, jj, ii, nonDirectionalCompressionCoefficient, vff)
    !----------------------------------------------------------------
    !   What it does:
    !   This subroutine computes the nondirectional compression 
    !   coefficient c for Weymouth and Yue's advection scheme.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(out) :: nonDirectionalCompressionCoefficient(kk, jj, ii)
        REAL(realk), INTENT(in) :: vff(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    IF ( vff(k,j,i) >= 0.5 ) THEN
                        nonDirectionalCompressionCoefficient(k,j,i) = 1.0
                    ELSE
                        nonDirectionalCompressionCoefficient(k,j,i) = 0.0
                    END IF
                END DO
            END DO
        END DO

    END SUBROUTINE compute_non_directional_compression_coeffiecient_pres

    !================================================================

    SUBROUTINE compute_non_directional_compression_coeffiecient_stag(kk, jj, ii, component, nonDirectionalCompressionCoefficient, vff)
    !----------------------------------------------------------------
    !   What it does:
    !   This subroutine computes the nondirectional compression 
    !   coefficient c for Weymouth and Yue's advection scheme.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: component
        REAL(realk), INTENT(out) :: nonDirectionalCompressionCoefficient(kk, jj, ii, 3)
        REAL(realk), INTENT(in) :: vff(kk, jj, ii, 3)

        ! Local variables
        INTEGER(intk) :: k, j, i

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    IF ( vff(k,j,i,component) >= 0.5 ) THEN
                        nonDirectionalCompressionCoefficient(k,j,i,component) = 1.0
                    ELSE
                        nonDirectionalCompressionCoefficient(k,j,i,component) = 0.0
                    END IF
                END DO
            END DO
        END DO

    END SUBROUTINE compute_non_directional_compression_coeffiecient_stag

    !================================================================

    SUBROUTINE compute_fluxx(fluxx, kk, jj, ii, field, isInterface, u, alpha, dt, normx, normy, normz, ddx, ddy, ddz, tol, & 
        nfro, nbac, nrgt, nlft, nbot, ntop)
    !----------------------------------------------------------------
    !   What it does:
    !   This subroutine calculates the flux of a field in the x
    !   direction.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(out) :: fluxx(kk, jj, ii)
        REAL(realk), INTENT(in) :: field(kk, jj, ii)
        LOGICAL, INTENT(in) :: isInterface(kk, jj, ii)
        REAL(realk), INTENT(in) :: u(kk, jj, ii)
        REAL(realk), INTENT(in) :: alpha(kk, jj, ii)
        REAL(realk), INTENT(in) :: dt
        REAL(realk), INTENT(in) :: normx(kk, jj, ii), normy(kk, jj, ii), normz(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: tol
        INTEGER, INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: nbu, nfu, nrv, nbw, ntw, nlv
        REAL(realk) :: fieldFlux, fluxedProportion, eulerianFluxWidth, eulerianFluxAlpha

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
                            ! Calculate the width of the fluxed volume and the alpha value for this subcell of the investigated cell
                            eulerianFluxWidth = abs( u(k,j,i) ) * dt
                            eulerianFluxAlpha = alpha(k,j,i) - normx(k,j,i) * ( ddx(i) - eulerianFluxWidth )

                            ! Caluculate the volume fraction in the fluxed volume subcell of the investigated cell 
                            CALL compute_cell_proportion(fluxedProportion, eulerianFluxAlpha, field(k,j,i), eulerianFluxWidth, ddy(j), ddz(k), normx(k,j,i), normy(k,j,i), normz(k,j,i), tol)
                            
                            ! Calculate flux for multiphase cell
                            fieldFlux = fluxedProportion  * ( abs( u(k,j,i) ) * dt / ddx(i) )
                        ELSE
                            ! Calculate flux for singlephase cell
                            fieldFlux = field(k,j,i) * ( abs( u(k,j,i) ) * dt / ddx(i) )
                        END IF
                    ELSE IF ( u(k,j,i) < -tol ) THEN
                        IF ( isInterface(k,j,i) ) THEN
                            ! Calculate the width of the fluxed volume and the alpha value for this subcell of the eastern cell
                            eulerianFluxWidth = abs( u(k,j,i) ) * dt
                            eulerianFluxAlpha = alpha(k,j,i+1)

                            ! Caluculate the volume fraction in the fluxed volume subcell of the eastern cell
                            CALL compute_cell_proportion(fluxedProportion, eulerianFluxAlpha, field(k,j,i+1), eulerianFluxWidth, ddy(j), ddz(k), normx(k,j,i+1), normy(k,j,i+1), normz(k,j,i+1), tol)

                            ! Calculate flux for multiphase cell
                            fieldFlux = fluxedProportion * ( abs( u(k,j,i) ) * dt / ddx(i+1) )
                        ELSE
                            ! Calculate flux for singlephase cell
                            fieldFlux = field(k,j,i+1) * ( abs( u(k,j,i) ) * dt / ddx(i+1) )
                        END IF
                    ELSE
                        fieldFlux = 0.0
                    END IF
                    fluxx(k,j,i) = sign( 1.0, u(k,j,i) ) * fieldFlux / dt
                END DO
            END DO
        END DO

    END SUBROUTINE compute_fluxx

    !================================================================

    SUBROUTINE compute_fluxy(fluxy, kk, jj, ii, field, isInterface, v, alpha, dt, normx, normy, normz, ddx, ddy, ddz, tol, & 
        nfro, nbac, nrgt, nlft, nbot, ntop)
    !----------------------------------------------------------------
    !   What it does:
    !   This subroutine calculates the flux of a field in the y 
    !   direction.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(out) :: fluxy(kk, jj, ii)
        REAL(realk), INTENT(in) :: field(kk, jj, ii)
        LOGICAL, INTENT(in) :: isInterface(kk, jj, ii)
        REAL(realk), INTENT(in) :: v(kk, jj, ii)
        REAL(realk), INTENT(in) :: alpha(kk, jj, ii)
        REAL(realk), INTENT(in) :: dt
        REAL(realk), INTENT(in) :: normx(kk, jj, ii), normy(kk, jj, ii), normz(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: tol
        INTEGER, INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: nbu, nfu, nrv, nbw, ntw, nlv
        REAL(realk) :: fieldFlux, fluxedProportion, eulerianFluxWidth, eulerianFluxAlpha

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
                            ! Calculate the width of the fluxed volume and the alpha value for this subcell of the investigated cell
                            eulerianFluxWidth = abs( v(k,j,i) ) * dt
                            eulerianFluxAlpha = alpha(k,j,i) - normy(k,j,i) * ( ddy(j) - eulerianFluxWidth )

                            ! Caluculate the volume fraction in the fluxed volume subcell of the investigated cell 
                            CALL compute_cell_proportion(fluxedProportion, eulerianFluxAlpha, field(k,j,i), ddx(i), eulerianFluxWidth, ddz(k), normx(k,j,i), normy(k,j,i), normz(k,j,i), tol)

                            ! Calculate flux for multiphase cell
                            fieldFlux = fluxedProportion * ( abs( v(k,j,i) ) * dt / ddy(j) )
                        ELSE
                            ! Calculate flux for singlephase cell
                            fieldFlux = field(k,j,i) * ( abs( v(k,j,i) ) * dt / ddy(j) )
                        END IF
                    ELSE IF ( v(k,j,i) < -tol ) THEN
                        IF ( isInterface(k,j,i) ) THEN
                            ! Calculate the width of the fluxed volume and the alpha value for this subcell of the northern cell
                            eulerianFluxWidth = abs( v(k,j,i) ) * dt
                            eulerianFluxAlpha = alpha(k,j+1,i)

                            ! Caluculate the volume fraction in the fluxed volume subcell of the northern cell 
                            CALL compute_cell_proportion(fluxedProportion, eulerianFluxAlpha, field(k,j+1,i), ddx(i), eulerianFluxWidth, ddz(k), normx(k,j+1,i), normy(k,j+1,i), normz(k,j+1,i), tol)

                            ! Calculate flux for multiphase cell
                            fieldFlux = fluxedProportion * ( abs( v(k,j,i) ) * dt / ddy(j+1) )
                        ELSE
                            ! Calculate flux for singlephase cell
                            fieldFlux = field(k,j+1,i) * ( abs( v(k,j,i) ) * dt / ddy(j+1) )
                        END IF
                    ELSE
                        fieldFlux = 0.0
                    END IF
                    fluxy(k,j,i) = sign( 1.0, v(k,j,i) ) * fieldFlux / dt
                END DO
            END DO
        END DO

    END SUBROUTINE compute_fluxy

    !================================================================

    SUBROUTINE compute_fluxz(fluxz, kk, jj, ii, field, isInterface, w, alpha, dt, normx, normy, normz, ddx, ddy, ddz, tol, & 
        nfro, nbac, nrgt, nlft, nbot, ntop)
    !----------------------------------------------------------------
    !   What it does:
    !   This subroutine calculates the flux of a field in the z 
    !   direction.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(out) :: fluxz(kk, jj, ii)
        REAL(realk), INTENT(in) :: field(kk, jj, ii)
        LOGICAL, INTENT(in) :: isInterface(kk, jj, ii)
        REAL(realk), INTENT(in) :: w(kk, jj, ii)
        REAL(realk), INTENT(in) :: alpha(kk, jj, ii)
        REAL(realk), INTENT(in) :: dt
        REAL(realk), INTENT(in) :: normx(kk, jj, ii), normy(kk, jj, ii), normz(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: tol
        INTEGER, INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: nbu, nfu, nrv, nbw, ntw, nlv
        REAL(realk) :: fieldFlux, fluxedProportion, eulerianFluxWidth, eulerianFluxAlpha

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
                            ! Calculate the width of the fluxed volume and the alpha value for this subcell of the investigated cell
                            eulerianFluxWidth = abs( w(k,j,i) ) * dt
                            eulerianFluxAlpha = alpha(k,j,i) - normz(k,j,i) * ( ddz(k) - eulerianFluxWidth )

                            ! Caluculate the volume fraction in the fluxed volume subcell of the investigated cell 
                            CALL compute_cell_proportion(fluxedProportion, eulerianFluxAlpha, field(k,j,i), ddx(i), ddy(j), eulerianFluxWidth, normx(k,j,i), normy(k,j,i), normz(k,j,i), tol)

                            ! Calculate flux for multiphase cell
                            fieldFlux = fluxedProportion * ( abs( w(k,j,i) ) * dt / ddz(k) )
                        ELSE
                            ! Calculate flux for singlephase cell
                            fieldFlux = field(k,j,i) * ( abs( w(k,j,i) ) * dt / ddz(k) )
                        END IF
                    ELSE IF ( w(k,j,i) < -tol ) THEN
                        IF ( isInterface(k,j,i) ) THEN
                            ! Calculate the width of the fluxed volume and the alpha value for this subcell of the investigated cell
                            eulerianFluxWidth = abs( w(k,j,i) ) * dt
                            eulerianFluxAlpha = alpha(k+1,j,i)

                            ! Caluculate the volume fraction in the fluxed volume subcell of the investigated cell 
                            CALL compute_cell_proportion(fluxedProportion, eulerianFluxAlpha, field(k+1,j,i), ddx(i), ddy(j), eulerianFluxWidth, normx(k+1,j,i), normy(k+1,j,i), normz(k+1,j,i), tol)

                            ! Calculate flux for multiphase cell
                            fieldFlux = fluxedProportion * ( abs( w(k,j,i) ) * dt / ddz(k+1) )
                        ELSE
                            ! Calculate flux for singlephase cell
                            fieldFlux = field(k+1,j,i) * ( abs( w(k,j,i) ) * dt / ddz(k+1) )
                        END IF
                    ELSE
                        fieldFlux = 0.0
                    END IF
                    fluxz(k,j,i) = sign( 1.0, w(k,j,i) ) * fieldFlux
                END DO
            END DO
        END DO

    END SUBROUTINE compute_fluxz

    !================================================================

    SUBROUTINE update_field(kk, jj, ii, splitDir, field, fluxx, fluxy, fluxz, xCompressionTerm, yCompressionTerm, zCompressionTerm, & 
            dt, nfro, nbac, nrgt, nlft, nbot, ntop)
    !----------------------------------------------------------------
    !   What it does:
    !   The subroutine performs the summation of fluxes, updating a 
    !   field. In each Runge-Kutta step the subroutine is called 
    !   three times. Each time field is updated taking into account 
    !   one spatial dimension. Each time-step the order of the 
    !   dimensional splitting is permuted by the calling subroutine. 
    !   The variables adv_(.) track, which dimension will be 
    !   integrated. 
    !   
    !   The time integration performed here is of 
    !   "geometrical" nature. It consideres real volumes, which are
    !   moved by the underlying velocity field and therefore is 
    !   exact. There is no need of a Runge-Kutta like canceling of an
    !   error term, hence previous Runge-Kutta stages are not 
    !   considered.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: splitDir
        REAL(realk), INTENT(inout) :: field(kk, jj, ii)
        REAL(realk), INTENT(in) :: fluxx(kk, jj, ii), fluxy(kk, jj, ii), fluxz(kk, jj, ii)
        REAL(realk), INTENT(in) :: xCompressionTerm(kk, jj, ii), yCompressionTerm(kk, jj, ii), zCompressionTerm(kk, jj, ii)
        REAL(realk), INTENT(in) :: dt
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
        
        ! Update field with x fluxes
        IF ( splitDir == 1 ) THEN 

            DO i = 3-nfu, ii-3+nbu
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        field(k,j,i) = ( field(k,j,i) + dt * ( fluxx(k,j,i-1) - fluxx(k,j,i) + xCompressionTerm(k,j,i) ) ) 
                    END DO 
                END DO 
            END DO

        ! Update field with y fluxes
        ELSEIF ( splitDir == 2 ) THEN

            DO i = 3, ii-2
                DO j = 3-nrv, jj-3+nlv
                    DO k = 3, kk-2
                        field(k,j,i) = ( field(k,j,i) + dt * ( fluxy(k,j-1,i) - fluxy(k,j,i) + yCompressionTerm(k,j,i) ) )
                    END DO 
                END DO 
            END DO

        ! Update field with z fluxes
        ELSEIF ( splitDir == 3 ) THEN

            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3-nbw, kk-3+ntw
                        field(k,j,i) = ( field(k,j,i) + dt * ( fluxz(k-1,j,i) - fluxz(k,j,i) + zCompressionTerm(k,j,i) ) )
                    END DO 
                END DO 
            END DO

        END IF

    END SUBROUTINE update_field

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
        INTEGER(intk) :: k, j, i
        
        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    IF ( vff(k,j,i) < tol ) THEN
                        vff(k,j,i) = 0.0
                    ELSE IF ( vff(k,j,i) > 1 - tol ) THEN
                        vff(k,j,i) = 1.0
                    END IF
                END DO
            END DO
        END DO

    END SUBROUTINE clip_volume_fraction_field

        !================================================================

    SUBROUTINE field_flux_wrapper_pres(kk, jj, ii, splitDir, field, isInterface, u, v, w, alpha, dt, normx, normy, normz, ddx, ddy, ddz, tol, nfro, nbac, nrgt, nlft, nbot, ntop, fieldFlux, complementFieldFlux)
    !----------------------------------------------------------------
    !   What it does:
    !   The subroutine is just a wrapper for the subroutines, which
    !   are used to compute the flux of a volume fraction field and
    !   its complement field in advection direction.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: splitDir
        REAL(realk), INTENT(in) :: field(kk, jj, ii)
        LOGICAL, INTENT(in) :: isInterface(kk, jj, ii)
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: alpha(kk, jj, ii)
        REAL(realk), INTENT(in) :: dt
        REAL(realk), INTENT(in) :: normx(kk, jj, ii), normy(kk, jj, ii), normz(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: tol
        INTEGER, INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop
        REAL(realk), INTENT(out) :: fieldFlux(kk, jj, ii)
        REAL(realk), INTENT(out) :: complementFieldFlux(kk, jj, ii)
 
        ! Local variables
        REAL(realk) :: complementField(kk, jj, ii)

        complementField = 1.0 - field
        
        IF ( splitDir == 1 ) THEN 

            CALL compute_fluxx(fieldFlux, kk, jj, ii, field, isInterface, u, alpha, dt, normx, normy, normz, ddx, ddy, ddz, tol, & 
                nfro, nbac, nrgt, nlft, nbot, ntop)
            CALL compute_fluxx(complementFieldFlux, kk, jj, ii, complementField, isInterface, u, alpha, dt, normx, normy, normz, ddx, ddy, ddz, tol, & 
                nfro, nbac, nrgt, nlft, nbot, ntop)

        ELSE IF ( splitDir == 2 ) THEN

            CALL compute_fluxy(fieldFlux, kk, jj, ii, field, isInterface, v, alpha, dt, normx, normy, normz, ddx, ddy, ddz, tol, & 
                nfro, nbac, nrgt, nlft, nbot, ntop)
            CALL compute_fluxy(complementFieldFlux, kk, jj, ii, complementField, isInterface, v, alpha, dt, normx, normy, normz, ddx, ddy, ddz, tol, & 
                nfro, nbac, nrgt, nlft, nbot, ntop)

        ELSE IF ( splitDir == 3 ) THEN

            CALL compute_fluxz(fieldFlux, kk, jj, ii, field, isInterface, w, alpha, dt, normx, normy, normz, ddx, ddy, ddz, tol, & 
                nfro, nbac, nrgt, nlft, nbot, ntop)
            CALL compute_fluxz(complementFieldFlux, kk, jj, ii, complementField, isInterface, w, alpha, dt, normx, normy, normz, ddx, ddy, ddz, tol, & 
                nfro, nbac, nrgt, nlft, nbot, ntop)

        END IF


    END SUBROUTINE field_flux_wrapper_pres

    !================================================================

    SUBROUTINE field_flux_wrapper_stag(kk, jj, ii, component, splitDir, field, isInterface, u, v, w, alpha, dt, normx, normy, normz, dx, dy, dz, ddx, ddy, ddz, tol, nfro, nbac, nrgt, nlft, nbot, ntop, fieldFlux, complementFieldFlux)
    !----------------------------------------------------------------
    !   What it does:
    !   The subroutine is just a wrapper for the subroutines, which
    !   are used to compute the flux of a volume fraction field and
    !   its complement field in advection direction.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: component, splitDir
        REAL(realk), INTENT(in) :: field(kk, jj, ii, 3)
        LOGICAL, INTENT(in) :: isInterface(kk, jj, ii, 3)
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: alpha(kk, jj, ii, 3)
        REAL(realk), INTENT(in) :: dt
        REAL(realk), INTENT(in) :: normx(kk, jj, ii, 3), normy(kk, jj, ii, 3), normz(kk, jj, ii, 3)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: tol
        INTEGER, INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop
        REAL(realk), INTENT(out) :: fieldFlux(kk, jj, ii, 3)
        REAL(realk), INTENT(out) :: complementFieldFlux(kk, jj, ii, 3)
 
        ! Local variables
        REAL(realk) :: complementField(kk, jj, ii, 3)
        REAL(realk) :: deltaX(ii), deltaY(jj), deltaZ(kk)

        IF ( component == 1 ) THEN
            deltaX = dx
            deltaY = ddy
            deltaZ = ddz
        ELSE IF ( component == 2 ) THEN
            deltaX = ddx
            deltaY = dy
            deltaZ = ddz
        ELSE IF ( component == 3 ) THEN
            deltaX = ddx
            deltaY = ddy
            deltaZ = dz
        END IF

        complementField = 1.0 - field
        
        IF ( splitDir == 1 ) THEN 

            CALL compute_fluxx(fieldFlux(:,:,:,component), kk, jj, ii, field(:,:,:,component), isInterface(:,:,:,component), u, alpha(:,:,:,component), dt, normx(:,:,:,component), normy(:,:,:,component), normz(:,:,:,component), deltaX, deltaY, deltaZ, tol, & 
                nfro, nbac, nrgt, nlft, nbot, ntop)
            CALL compute_fluxx(complementFieldFlux(:,:,:,component), kk, jj, ii, complementField(:,:,:,component), isInterface(:,:,:,component), u, alpha(:,:,:,component), dt, normx(:,:,:,component), normy(:,:,:,component), normz(:,:,:,component), deltaX, deltaY, deltaZ, tol, & 
                nfro, nbac, nrgt, nlft, nbot, ntop)

        ELSE IF ( splitDir == 2 ) THEN

            CALL compute_fluxy(fieldFlux(:,:,:,component), kk, jj, ii, field(:,:,:,component), isInterface(:,:,:,component), v, alpha(:,:,:,component), dt, normx(:,:,:,component), normy(:,:,:,component), normz(:,:,:,component), deltaX, deltaY, deltaZ, tol, & 
                nfro, nbac, nrgt, nlft, nbot, ntop)
            CALL compute_fluxy(complementFieldFlux(:,:,:,component), kk, jj, ii, complementField(:,:,:,component), isInterface(:,:,:,component), v, alpha(:,:,:,component), dt, normx(:,:,:,component), normy(:,:,:,component), normz(:,:,:,component), deltaX, deltaY, deltaZ, tol, & 
                nfro, nbac, nrgt, nlft, nbot, ntop)

        ELSE IF ( splitDir == 3 ) THEN

            CALL compute_fluxz(fieldFlux(:,:,:,component), kk, jj, ii, field(:,:,:,component), isInterface(:,:,:,component), w, alpha(:,:,:,component), dt, normx(:,:,:,component), normy(:,:,:,component), normz(:,:,:,component), deltaX, deltaY, deltaZ, tol, & 
                nfro, nbac, nrgt, nlft, nbot, ntop)
            CALL compute_fluxz(complementFieldFlux(:,:,:,component), kk, jj, ii, complementField(:,:,:,component), isInterface(:,:,:,component), w, alpha(:,:,:,component), dt, normx(:,:,:,component), normy(:,:,:,component), normz(:,:,:,component), deltaX, deltaY, deltaZ, tol, & 
                nfro, nbac, nrgt, nlft, nbot, ntop)

        END IF


    END SUBROUTINE field_flux_wrapper_stag

    !================================================================

    SUBROUTINE compute_density_flux(kk, jj, ii, component, flux, fluxComp, rho1, rho2, densityFlux)
    !----------------------------------------------------------------
    !   What it does:
    !   The subroutine computes the density fluxes using the volume
    !   fraction field fluxes and their complements.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: component
        REAL(realk), INTENT(in) :: flux(kk, jj, ii, 3), fluxComp(kk, jj, ii, 3)
        REAL(realk), INTENT(in) :: rho1, rho2
        REAL(realk), INTENT(out) :: densityFlux(kk, jj, ii, 3)

        ! Local variables
        INTEGER(intk) :: k, j, i
        
        DO i = 1, ii
            DO j = 1, jj
                DO k = 1, kk
                    densityFlux(k,j,i,component) = rho1 * flux(k,j,i,component) + rho2 * fluxComp(k,j,i,component)
                END DO
            END DO
        END DO

    END SUBROUTINE compute_density_flux

    !================================================================
    
    SUBROUTINE get_advection_sequence(iteration, advSeq)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: iteration
        INTEGER(intk), INTENT(inout) :: advSeq(3)

        ! Local variables
        INTEGER(intk) :: permutationIndex

        ! permutationIndex only changes in a new time-step
        permutationIndex = mod(iteration-1, 3)

        ! Select permutation of split advection
        SELECT CASE (permutationIndex)
            CASE (0)
                advSeq = [1, 2, 3]
            CASE (1)
                advSeq = [3, 1, 2]
            CASE (2)
                advSeq = [2, 3, 1]
        END SELECT

    END SUBROUTINE get_advection_sequence

    !================================================================

    SUBROUTINE compression_term_wrapper_pres(kk, jj, ii, u, v, w, vff, ddx, ddy, ddz, compressionTermX, compressionTermY, compressionTermZ, propertyFluid1, propertyFluid2)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii), vff(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(out) :: compressionTermX(kk, jj, ii), compressionTermY(kk, jj, ii), compressionTermZ(kk, jj, ii)
        REAL(realk), INTENT(in), OPTIONAL :: propertyFluid1, propertyFluid2

        ! Local variables
        REAL(realk) :: strainRateX(kk, jj, ii), strainRateY(kk, jj, ii), strainRateZ(kk, jj, ii)
        REAL(realk) :: nonDirectionalCompressionCoefficient(kk, jj, ii)

        CALL compute_normal_strain_rates(kk, jj, ii, strainRateX, strainRateY, strainRateZ, u, v, w, ddx, ddy, ddz)
        CALL compute_non_directional_compression_coeffiecient(kk, jj, ii, nonDirectionalCompressionCoefficient, vff)

        IF ( .NOT. PRESENT(propertyFluid1) .OR. .NOT. PRESENT(propertyFluid2) ) THEN
            compressionTermX = nonDirectionalCompressionCoefficient * strainRateX
            compressionTermY = nonDirectionalCompressionCoefficient * strainRateY
            compressionTermZ = nonDirectionalCompressionCoefficient * strainRateZ
        ELSE IF ( PRESENT(propertyFluid1) .AND. PRESENT(propertyFluid2) ) THEN
            compressionTermX = ( nonDirectionalCompressionCoefficient * propertyFluid1 + ( 1.0 - nonDirectionalCompressionCoefficient ) * propertyFluid2 ) * strainRateX
            compressionTermY = ( nonDirectionalCompressionCoefficient * propertyFluid1 + ( 1.0 - nonDirectionalCompressionCoefficient ) * propertyFluid2 ) * strainRateY
            compressionTermZ = ( nonDirectionalCompressionCoefficient * propertyFluid1 + ( 1.0 - nonDirectionalCompressionCoefficient ) * propertyFluid2 ) * strainRateZ
        END IF

    END SUBROUTINE compression_term_wrapper_pres

    !================================================================

    SUBROUTINE compression_term_wrapper_stag(kk, jj, ii, component, u, v, w, vff, dx, dy, dz, ddx, ddy, ddz, compressionTermX, compressionTermY, compressionTermZ, propertyFluid1, propertyFluid2)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: component
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii), vff(kk, jj, ii, 3)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(out) :: compressionTermX(kk, jj, ii, 3), compressionTermY(kk, jj, ii, 3), compressionTermZ(kk, jj, ii, 3)
        REAL(realk), INTENT(in), OPTIONAL :: propertyFluid1, propertyFluid2

        ! Local variables
        REAL(realk) :: strainRateX(kk, jj, ii, 3), strainRateY(kk, jj, ii, 3), strainRateZ(kk, jj, ii, 3)
        REAL(realk) :: nonDirectionalCompressionCoefficient(kk, jj, ii, 3)

        CALL compute_normal_strain_rates(kk, jj, ii, component, strainRateX, strainRateY, strainRateZ, u, v, w, dx, dy, dz, ddx, ddy, ddz)
        CALL compute_non_directional_compression_coeffiecient(kk, jj, ii, component, nonDirectionalCompressionCoefficient, vff)

        IF ( .NOT. PRESENT(propertyFluid1) .OR. .NOT. PRESENT(propertyFluid2) ) THEN
            compressionTermX(:,:,:,component) = nonDirectionalCompressionCoefficient(:,:,:,component) * strainRateX(:,:,:,component)
            compressionTermY(:,:,:,component) = nonDirectionalCompressionCoefficient(:,:,:,component) * strainRateY(:,:,:,component)
            compressionTermZ(:,:,:,component) = nonDirectionalCompressionCoefficient(:,:,:,component) * strainRateZ(:,:,:,component)
        ELSE IF ( PRESENT(propertyFluid1) .AND. PRESENT(propertyFluid2) ) THEN
            compressionTermX(:,:,:,component) = ( nonDirectionalCompressionCoefficient(:,:,:,component) * propertyFluid1 + ( 1.0 - nonDirectionalCompressionCoefficient(:,:,:,component) ) * propertyFluid2 ) * strainRateX(:,:,:,component)
            compressionTermY(:,:,:,component) = ( nonDirectionalCompressionCoefficient(:,:,:,component) * propertyFluid1 + ( 1.0 - nonDirectionalCompressionCoefficient(:,:,:,component) ) * propertyFluid2 ) * strainRateY(:,:,:,component)
            compressionTermZ(:,:,:,component) = ( nonDirectionalCompressionCoefficient(:,:,:,component) * propertyFluid1 + ( 1.0 - nonDirectionalCompressionCoefficient(:,:,:,component) ) * propertyFluid2 ) * strainRateZ(:,:,:,component)
        END IF

    END SUBROUTINE compression_term_wrapper_stag

END MODULE multiphase_vof_transport_mod
