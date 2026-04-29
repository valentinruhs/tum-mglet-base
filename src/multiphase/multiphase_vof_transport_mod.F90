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

    PUBLIC :: init_multiphase_vof_transport, finish_multiphase_vof_transport, compute_flux, compute_density_flux, get_advection_sequence, compression_term_wrapper, update_field, clip_volume_fraction_field

    INTERFACE compute_normal_strain_rate
        MODULE PROCEDURE compute_normal_strain_rate_pres
        MODULE PROCEDURE compute_normal_strain_rate_stag
    END INTERFACE

    INTERFACE compression_term_wrapper
        MODULE PROCEDURE compression_term_wrapper_pres
        MODULE PROCEDURE compression_term_wrapper_stag
    END INTERFACE

    INTERFACE compute_non_directional_compression_coeffiecient
        MODULE PROCEDURE compute_non_directional_compression_coeffiecient_pres
        MODULE PROCEDURE compute_non_directional_compression_coeffiecient_stag
    END INTERFACE

    INTERFACE update_field
        MODULE PROCEDURE update_field_pres
        MODULE PROCEDURE update_field_stag
    END INTERFACE

    INTERFACE compute_flux
        MODULE PROCEDURE compute_flux_pres
        MODULE PROCEDURE compute_flux_stag
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

    SUBROUTINE compute_normal_strain_rate_pres(kk, jj, ii, splitDir, strainRate, u, v, w, ddx, ddy, ddz)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: splitDir
        REAL(realk), INTENT(out) :: strainRate(kk, jj, ii)
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)

        ! Local variables
        INTEGER(intk) :: k, j, i

        IF ( splitDir == 1 ) THEN
            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        strainRate(k,j,i) = ( u(k,j,i) - u(k,j,i-1) ) / ddx(i)
                    END DO
                END DO
            END DO
        ELSE IF ( splitDir == 2 ) THEN
            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        strainRate(k,j,i) = ( v(k,j,i) - v(k,j-1,i) ) / ddy(j)
                    END DO
                END DO
            END DO
        ELSE IF ( splitDir == 3 ) THEN 
            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        strainRate(k,j,i) = ( w(k,j,i) - w(k-1,j,i) ) / ddz(k)
                    END DO
                END DO
            END DO
        END IF

    END SUBROUTINE compute_normal_strain_rate_pres

    !================================================================

    SUBROUTINE compute_normal_strain_rate_stag(kk, jj, ii, component, splitDir, strainRate, u, v, w, dx, dy, dz, ddx, ddy, ddz)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: component
        INTEGER(intk), INTENT(in) :: splitDir
        REAL(realk), INTENT(out) :: strainRate(kk, jj, ii, 3)
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

        IF ( splitDir == 1 ) THEN
            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        uE = 0.5 * ( iStag * ( u(k,j,i) + u(k,j,i+1) ) + jStag * ( u(k,j,i) + u(k,j+1,i) ) + kStag * ( u(k,j,i) + u(k+1,j,i) ) ) 
                        uW = 0.5 * ( iStag * ( u(k,j,i-1) + u(k,j,i) ) + jStag * ( u(k,j,i-1) + u(k,j+1,i-1) ) + kStag * ( u(k,j,i-1) + u(k+1,j,i-1) ) )
                        strainRate(k,j,i,component) = ( uE - uW ) / ( iStag * dx(i) + jStag * ddx(i) + kStag * ddx(i) )
                    END DO
                END DO
            END DO
        ELSE IF ( splitDir == 2 ) THEN
            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        vN = 0.5 * ( iStag * ( v(k,j,i) + v(k,j,i+1) ) + jStag * ( v(k,j,i) + v(k,j+1,i) ) + kStag * ( v(k,j,i) + v(k+1,j,i) ) )
                        vS = 0.5 * ( iStag * ( v(k,j-1,i) + v(k,j-1,i+1) ) + jStag * ( v(k,j-1,i) + v(k,j,i) ) + kStag * ( v(k+1,j-1,i) + v(k,j-1,i) ) )
                        strainRate(k,j,i,component) = ( vN - vS ) / ( iStag * ddy(j) + jStag * dy(j) + kStag * ddy(j) )
                    END DO
                END DO
            END DO
        ELSE IF ( splitDir == 3 ) THEN 
            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        wT = 0.5 * ( iStag * ( w(k,j,i) + w(k,j,i+1) ) + jStag * ( w(k,j,i) + w(k,j+1,i) ) + kStag * ( w(k,j,i) + w(k+1,j,i) ) ) 
                        wB = 0.5 * ( iStag * ( w(k-1,j,i) + w(k-1,j,i+1) ) + jStag * ( w(k-1,j,i) + w(k-1,j+1,i) ) + kStag * ( w(k-1,j,i) + w(k,j,i) ) ) 
                        strainRate(k,j,i,component) = ( wT - wB ) / ( iStag * ddz(k) + jStag * ddz(k) + kStag * dz(k) )                    
                    END DO
                END DO
            END DO
        END IF

    END SUBROUTINE compute_normal_strain_rate_stag
    
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

    SUBROUTINE compute_flux_pres(kk, jj, ii, splitDir, fieldFlux, field, isInterface, u, v, w, alpha, dt, normx, normy, normz, ddx, ddy, ddz, tol, & 
        nfro, nbac, nrgt, nlft, nbot, ntop)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: splitDir
        REAL(realk), INTENT(out) :: fieldFlux(kk, jj, ii)
        REAL(realk), INTENT(in) :: field(kk, jj, ii)
        LOGICAL, INTENT(in) :: isInterface(kk, jj, ii)
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: alpha(kk, jj, ii)
        REAL(realk), INTENT(in) :: dt
        REAL(realk), INTENT(in) :: normx(kk, jj, ii), normy(kk, jj, ii), normz(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: tol
        INTEGER, INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: nbu, nfu, nrv, nbw, ntw, nlv
        REAL(realk) :: flux, fluxedProportion, eulerianFluxWidth, eulerianFluxAlpha

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

        IF ( splitDir == 1 ) THEN 
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
                                flux = fluxedProportion  * ( abs( u(k,j,i) ) * dt / ddx(i) )
                            ELSE
                                ! Calculate flux for singlephase cell
                                flux = field(k,j,i) * ( abs( u(k,j,i) ) * dt / ddx(i) )
                            END IF
                        ELSE IF ( u(k,j,i) < -tol ) THEN
                            IF ( isInterface(k,j,i) ) THEN
                                ! Calculate the width of the fluxed volume and the alpha value for this subcell of the eastern cell
                                eulerianFluxWidth = abs( u(k,j,i) ) * dt
                                eulerianFluxAlpha = alpha(k,j,i+1)

                                ! Caluculate the volume fraction in the fluxed volume subcell of the eastern cell
                                CALL compute_cell_proportion(fluxedProportion, eulerianFluxAlpha, field(k,j,i+1), eulerianFluxWidth, ddy(j), ddz(k), normx(k,j,i+1), normy(k,j,i+1), normz(k,j,i+1), tol)

                                ! Calculate flux for multiphase cell
                                flux = fluxedProportion * ( abs( u(k,j,i) ) * dt / ddx(i+1) )
                            ELSE
                                ! Calculate flux for singlephase cell
                                flux = field(k,j,i+1) * ( abs( u(k,j,i) ) * dt / ddx(i+1) )
                            END IF
                        ELSE
                            flux = 0.0
                        END IF
                        fieldFlux(k,j,i) = sign( 1.0, u(k,j,i) ) * flux / dt
                    END DO
                END DO
            END DO
        ELSE IF ( splitDir == 2 ) THEN 
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
                                flux = fluxedProportion * ( abs( v(k,j,i) ) * dt / ddy(j) )
                            ELSE
                                ! Calculate flux for singlephase cell
                                flux = field(k,j,i) * ( abs( v(k,j,i) ) * dt / ddy(j) )
                            END IF
                        ELSE IF ( v(k,j,i) < -tol ) THEN
                            IF ( isInterface(k,j,i) ) THEN
                                ! Calculate the width of the fluxed volume and the alpha value for this subcell of the northern cell
                                eulerianFluxWidth = abs( v(k,j,i) ) * dt
                                eulerianFluxAlpha = alpha(k,j+1,i)

                                ! Caluculate the volume fraction in the fluxed volume subcell of the northern cell 
                                CALL compute_cell_proportion(fluxedProportion, eulerianFluxAlpha, field(k,j+1,i), ddx(i), eulerianFluxWidth, ddz(k), normx(k,j+1,i), normy(k,j+1,i), normz(k,j+1,i), tol)

                                ! Calculate flux for multiphase cell
                                flux = fluxedProportion * ( abs( v(k,j,i) ) * dt / ddy(j+1) )
                            ELSE
                                ! Calculate flux for singlephase cell
                                flux = field(k,j+1,i) * ( abs( v(k,j,i) ) * dt / ddy(j+1) )
                            END IF
                        ELSE
                            flux = 0.0
                        END IF
                        fieldFlux(k,j,i) = sign( 1.0, v(k,j,i) ) * flux / dt
                    END DO
                END DO
            END DO
        ELSE IF ( splitDir == 3 ) THEN
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
                                flux = fluxedProportion * ( abs( w(k,j,i) ) * dt / ddz(k) )
                            ELSE
                                ! Calculate flux for singlephase cell
                                flux = field(k,j,i) * ( abs( w(k,j,i) ) * dt / ddz(k) )
                            END IF
                        ELSE IF ( w(k,j,i) < -tol ) THEN
                            IF ( isInterface(k,j,i) ) THEN
                                ! Calculate the width of the fluxed volume and the alpha value for this subcell of the investigated cell
                                eulerianFluxWidth = abs( w(k,j,i) ) * dt
                                eulerianFluxAlpha = alpha(k+1,j,i)

                                ! Caluculate the volume fraction in the fluxed volume subcell of the investigated cell 
                                CALL compute_cell_proportion(fluxedProportion, eulerianFluxAlpha, field(k+1,j,i), ddx(i), ddy(j), eulerianFluxWidth, normx(k+1,j,i), normy(k+1,j,i), normz(k+1,j,i), tol)

                                ! Calculate flux for multiphase cell
                                flux = fluxedProportion * ( abs( w(k,j,i) ) * dt / ddz(k+1) )
                            ELSE
                                ! Calculate flux for singlephase cell
                                flux = field(k+1,j,i) * ( abs( w(k,j,i) ) * dt / ddz(k+1) )
                            END IF
                        ELSE
                            flux = 0.0
                        END IF
                        fieldFlux(k,j,i) = sign( 1.0, w(k,j,i) ) * flux / dt
                    END DO
                END DO
            END DO
        END IF

    END SUBROUTINE compute_flux_pres

    !================================================================

    SUBROUTINE compute_flux_stag(kk, jj, ii, component, splitDir, fieldFlux, complementFieldFlux, field, isInterface, u, v, w, alpha, dt, normx, normy, normz, dx, dy, dz, ddx, ddy, ddz, tol, & 
        nfro, nbac, nrgt, nlft, nbot, ntop)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: component
        INTEGER(intk), INTENT(in) :: splitDir
        REAL(realk), INTENT(out) :: fieldFlux(kk, jj, ii, 3), complementFieldFlux(kk, jj, ii, 3)
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

        ! Local variables
        REAL(realk) :: iStag, jStag, kStag
        REAL(realk) :: deltaX(ii), deltaY(jj), deltaZ(kk)
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: nbu, nfu, nrv, nbw, ntw, nlv
        REAL(realk) :: flux, complementFlux, fluxedProportion, eulerianFluxWidth, eulerianFluxAlpha, uE, vN, wT

        IF ( component == 1 ) THEN
            iStag = 1.0
            jStag = 0.0
            kStag = 0.0
            deltaX = dx
            deltaY = ddy
            deltaZ = ddz
        ELSE IF ( component == 2 ) THEN
            iStag = 0.0
            jStag = 1.0
            kStag = 0.0
            deltaX = ddx
            deltaY = dy
            deltaZ = ddz
        ELSE IF ( component == 3 ) THEN
            iStag = 0.0
            jStag = 0.0
            kStag = 1.0
            deltaX = ddx
            deltaY = ddy
            deltaZ = dz
        END IF

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

        IF ( splitDir == 1 ) THEN
            DO i = 3-nfu, ii-3+nbu
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        uE = 0.5 * ( iStag * ( u(k,j,i) + u(k,j,i+1) ) + jStag * ( u(k,j,i) + u(k,j+1,i) ) + kStag * ( u(k,j,i) + u(k+1,j,i) ) ) 
                        IF ( uE > tol ) THEN
                            IF ( isInterface(k,j,i,component) ) THEN
                                ! Calculate the width of the fluxed volume and the alpha value for this subcell of the investigated cell
                                eulerianFluxWidth = abs( uE ) * dt
                                eulerianFluxAlpha = alpha(k,j,i,component) - normx(k,j,i,component) * ( deltaX(i) - eulerianFluxWidth )

                                ! Caluculate the volume fraction in the fluxed volume subcell of the investigated cell 
                                CALL compute_cell_proportion(fluxedProportion, eulerianFluxAlpha, field(k,j,i,component), eulerianFluxWidth, ddy(j), ddz(k), normx(k,j,i,component), normy(k,j,i,component), normz(k,j,i,component), tol)
                                
                                ! Calculate flux for multiphase cell
                                flux = fluxedProportion  * ( abs( uE ) * dt / deltaX(i) )
                                complementFlux = ( 1 - fluxedProportion ) * ( abs( uE ) * dt / deltaX(i) )
                            ELSE
                                ! Calculate flux for singlephase cell
                                flux = field(k,j,i,component) * ( abs( uE ) * dt / deltaX(i) )
                                complementFlux = ( 1 - field(k,j,i,component) ) * ( abs( uE ) * dt / deltaX(i) )
                            END IF
                        ELSE IF ( uE < -tol ) THEN
                            IF ( isInterface(k,j,i,component) ) THEN
                                ! Calculate the width of the fluxed volume and the alpha value for this subcell of the eastern cell
                                eulerianFluxWidth = abs( uE ) * dt
                                eulerianFluxAlpha = alpha(k,j,i+1,component)

                                ! Caluculate the volume fraction in the fluxed volume subcell of the eastern cell
                                CALL compute_cell_proportion(fluxedProportion, eulerianFluxAlpha, field(k,j,i+1,component), eulerianFluxWidth, ddy(j), ddz(k), normx(k,j,i+1,component), normy(k,j,i+1,component), normz(k,j,i+1,component), tol)

                                ! Calculate flux for multiphase cell
                                flux = fluxedProportion * ( abs( uE ) * dt / deltaX(i+1) )
                                complementFlux = ( 1 - fluxedProportion ) * ( abs( uE ) * dt / deltaX(i+1) )
                            ELSE
                                ! Calculate flux for singlephase cell
                                flux = field(k,j,i+1,component) * ( abs( uE ) * dt / deltaX(i+1) )
                                complementFlux = ( 1 - field(k,j,i+1,component) ) * ( abs( uE ) * dt / deltaX(i+1) )
                            END IF
                        ELSE
                            flux = 0.0
                            complementFlux = 0.0
                        END IF
                        fieldFlux(k,j,i,component) = sign( 1.0, uE ) * flux / dt
                        complementFieldFlux(k,j,i,component) = sign( 1.0, uE ) * complementFlux / dt
                    END DO
                END DO
            END DO
        ELSE IF ( splitDir == 2 ) THEN
            DO i = 3, ii-2
                DO j = 3-nrv, jj-3+nlv
                    DO k = 3, kk-2
                        vN = 0.5 * ( iStag * ( v(k,j,i) + v(k,j,i+1) ) + jStag * ( v(k,j,i) + v(k,j+1,i) ) + kStag * ( v(k,j,i) + v(k+1,j,i) ) )
                        IF ( vN > tol ) THEN
                            IF ( isInterface(k,j,i,component) ) THEN
                                ! Calculate the width of the fluxed volume and the alpha value for this subcell of the investigated cell
                                eulerianFluxWidth = abs( vN ) * dt
                                eulerianFluxAlpha = alpha(k,j,i,component) - normy(k,j,i,component) * ( deltaY(j) - eulerianFluxWidth )

                                ! Caluculate the volume fraction in the fluxed volume subcell of the investigated cell 
                                CALL compute_cell_proportion(fluxedProportion, eulerianFluxAlpha, field(k,j,i,component), ddx(i), eulerianFluxWidth, ddz(k), normx(k,j,i,component), normy(k,j,i,component), normz(k,j,i,component), tol)

                                ! Calculate flux for multiphase cell
                                flux = fluxedProportion * ( abs( vN ) * dt / deltaY(j) )
                                complementFlux = ( 1 - fluxedProportion ) * ( abs( vN ) * dt / deltaY(j) )
                            ELSE
                                ! Calculate flux for singlephase cell
                                flux = field(k,j,i,component) * ( abs( vN ) * dt / deltaY(j) )
                                complementFlux = ( 1 - field(k,j,i,component) ) * ( abs( vN ) * dt / deltaY(j) )
                            END IF
                        ELSE IF ( vN < -tol ) THEN
                            IF ( isInterface(k,j,i,component) ) THEN
                                ! Calculate the width of the fluxed volume and the alpha value for this subcell of the northern cell
                                eulerianFluxWidth = abs( vN ) * dt
                                eulerianFluxAlpha = alpha(k,j+1,i,component)

                                ! Caluculate the volume fraction in the fluxed volume subcell of the northern cell 
                                CALL compute_cell_proportion(fluxedProportion, eulerianFluxAlpha, field(k,j+1,i,component), ddx(i), eulerianFluxWidth, ddz(k), normx(k,j+1,i,component), normy(k,j+1,i,component), normz(k,j+1,i,component), tol)

                                ! Calculate flux for multiphase cell
                                flux = fluxedProportion * ( abs( vN ) * dt / deltaY(j+1) )
                                complementFlux = ( 1 - fluxedProportion ) * ( abs( vN ) * dt / deltaY(j+1) )
                            ELSE
                                ! Calculate flux for singlephase cell
                                flux = field(k,j+1,i,component) * ( abs( vN ) * dt / deltaY(j+1) )
                                complementFlux = ( 1 - field(k,j+1,i,component) ) * ( abs( vN ) * dt / deltaY(j+1) )
                            END IF
                        ELSE
                            flux = 0.0
                            complementFlux = 0.0
                        END IF
                        fieldFlux(k,j,i,component) = sign( 1.0, vN ) * flux / dt
                        complementFieldFlux(k,j,i,component) = sign( 1.0, vN ) * complementFlux / dt
                    END DO
                END DO
            END DO
        ELSE IF ( splitDir == 3 ) THEN
            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3-nbw, kk-3+ntw
                        wT = 0.5 * ( iStag * ( w(k,j,i) + w(k,j,i+1) ) + jStag * ( w(k,j,i) + w(k,j+1,i) ) + kStag * ( w(k,j,i) + w(k+1,j,i) ) )
                        IF ( wT > tol ) THEN
                            IF ( isInterface(k,j,i,component) ) THEN
                                ! Calculate the width of the fluxed volume and the alpha value for this subcell of the investigated cell
                                eulerianFluxWidth = abs( wT ) * dt
                                eulerianFluxAlpha = alpha(k,j,i,component) - normz(k,j,i,component) * ( deltaZ(k) - eulerianFluxWidth )

                                ! Caluculate the volume fraction in the fluxed volume subcell of the investigated cell 
                                CALL compute_cell_proportion(fluxedProportion, eulerianFluxAlpha, field(k,j,i,component), ddx(i), ddy(j), eulerianFluxWidth, normx(k,j,i,component), normy(k,j,i,component), normz(k,j,i,component), tol)

                                ! Calculate flux for multiphase cell
                                flux = fluxedProportion * ( abs( wT ) * dt / deltaZ(k) )
                                complementFlux = ( 1 - fluxedProportion ) * ( abs( wT ) * dt / deltaZ(k) )
                            ELSE
                                ! Calculate flux for singlephase cell
                                flux = field(k,j,i,component) * ( abs( wT ) * dt / deltaZ(k) )
                                complementFlux = ( 1 - field(k,j,i,component) ) * ( abs( wT ) * dt / deltaZ(k) )
                            END IF
                        ELSE IF ( wT < -tol ) THEN
                            IF ( isInterface(k,j,i,component) ) THEN
                                ! Calculate the width of the fluxed volume and the alpha value for this subcell of the investigated cell
                                eulerianFluxWidth = abs( wT ) * dt
                                eulerianFluxAlpha = alpha(k+1,j,i,component)

                                ! Caluculate the volume fraction in the fluxed volume subcell of the investigated cell 
                                CALL compute_cell_proportion(fluxedProportion, eulerianFluxAlpha, field(k+1,j,i,component), ddx(i), ddy(j), eulerianFluxWidth, normx(k+1,j,i,component), normy(k+1,j,i,component), normz(k+1,j,i,component), tol)

                                ! Calculate flux for multiphase cell
                                flux = fluxedProportion * ( abs( wT ) * dt / deltaZ(k+1) )
                                complementFlux = ( 1 - fluxedProportion ) * ( abs( wT ) * dt / deltaZ(k+1) )
                            ELSE
                                ! Calculate flux for singlephase cell
                                flux = field(k+1,j,i,component) * ( abs( wT ) * dt / deltaZ(k+1) )
                                complementFlux = ( 1 - field(k+1,j,i,component) ) * ( abs( wT ) * dt / deltaZ(k+1) )
                            END IF
                        ELSE
                            flux = 0.0
                            complementFLux = 0.0
                        END IF
                        fieldFlux(k,j,i,component) = sign( 1.0, wT ) * flux / dt
                        complementFieldFlux(k,j,i,component) = sign( 1.0, wT ) * complementFlux / dt
                    END DO
                END DO
            END DO
        END IF

    END SUBROUTINE compute_flux_stag

    !================================================================

    SUBROUTINE update_field_pres(kk, jj, ii, splitDir, field, flux, compressionTerm, & 
            dt, nfro, nbac, nrgt, nlft, nbot, ntop)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: splitDir
        REAL(realk), INTENT(inout) :: field(kk, jj, ii)
        REAL(realk), INTENT(in) :: flux(kk, jj, ii)
        REAL(realk), INTENT(in) :: compressionTerm(kk, jj, ii)
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
        
        IF ( splitDir == 1 ) THEN 
            DO i = 3-nfu, ii-3+nbu
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        field(k,j,i) = ( field(k,j,i) + dt * ( flux(k,j,i-1) - flux(k,j,i) + compressionTerm(k,j,i) ) ) 
                    END DO 
                END DO 
            END DO
        ELSEIF ( splitDir == 2 ) THEN
            DO i = 3, ii-2
                DO j = 3-nrv, jj-3+nlv
                    DO k = 3, kk-2
                        field(k,j,i) = ( field(k,j,i) + dt * ( flux(k,j-1,i) - flux(k,j,i) + compressionTerm(k,j,i) ) )
                    END DO 
                END DO 
            END DO
        ELSEIF ( splitDir == 3 ) THEN
            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3-nbw, kk-3+ntw
                        field(k,j,i) = ( field(k,j,i) + dt * ( flux(k-1,j,i) - flux(k,j,i) + compressionTerm(k,j,i) ) )
                    END DO 
                END DO 
            END DO
        END IF

    END SUBROUTINE update_field_pres

    !================================================================

    SUBROUTINE update_field_stag(kk, jj, ii, component, splitDir, field, flux, compressionTerm, & 
            dt, nfro, nbac, nrgt, nlft, nbot, ntop)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: component
        INTEGER(intk), INTENT(in) :: splitDir
        REAL(realk), INTENT(inout) :: field(kk, jj, ii, 3)
        REAL(realk), INTENT(in) :: flux(kk, jj, ii, 3)
        REAL(realk), INTENT(in) :: compressionTerm(kk, jj, ii, 3)
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
        
        IF ( splitDir == 1 ) THEN 
            DO i = 3-nfu, ii-3+nbu
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        field(k,j,i,component) = ( field(k,j,i,component) + dt * ( flux(k,j,i-1,component) - flux(k,j,i,component) + compressionTerm(k,j,i,component) ) ) 
                    END DO 
                END DO 
            END DO
        ELSEIF ( splitDir == 2 ) THEN
            DO i = 3, ii-2
                DO j = 3-nrv, jj-3+nlv
                    DO k = 3, kk-2
                        field(k,j,i,component) = ( field(k,j,i,component) + dt * ( flux(k,j-1,i,component) - flux(k,j,i,component) + compressionTerm(k,j,i,component) ) )
                    END DO 
                END DO 
            END DO
        ELSEIF ( splitDir == 3 ) THEN
            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3-nbw, kk-3+ntw
                        field(k,j,i,component) = ( field(k,j,i,component) + dt * ( flux(k-1,j,i,component) - flux(k,j,i,component) + compressionTerm(k,j,i,component) ) )
                    END DO 
                END DO 
            END DO
        END IF

    END SUBROUTINE update_field_stag

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

    SUBROUTINE compression_term_wrapper_pres(kk, jj, ii, splitDir, u, v, w, vff, ddx, ddy, ddz, compressionTerm, propertyFluid1, propertyFluid2)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: splitDir
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii), vff(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(out) :: compressionTerm(kk, jj, ii)
        REAL(realk), INTENT(in), OPTIONAL :: propertyFluid1, propertyFluid2

        ! Local variables
        REAL(realk) :: strainRate(kk, jj, ii)
        REAL(realk) :: nonDirectionalCompressionCoefficient(kk, jj, ii)

        CALL compute_normal_strain_rate(kk, jj, ii, splitDir, strainRate, u, v, w, ddx, ddy, ddz)
        CALL compute_non_directional_compression_coeffiecient(kk, jj, ii, nonDirectionalCompressionCoefficient, vff)

        IF ( .NOT. PRESENT(propertyFluid1) .OR. .NOT. PRESENT(propertyFluid2) ) THEN
            compressionTerm = nonDirectionalCompressionCoefficient * strainRate
        ELSE IF ( PRESENT(propertyFluid1) .AND. PRESENT(propertyFluid2) ) THEN
            compressionTerm = ( nonDirectionalCompressionCoefficient * propertyFluid1 + ( 1.0 - nonDirectionalCompressionCoefficient ) * propertyFluid2 ) * strainRate
        END IF

    END SUBROUTINE compression_term_wrapper_pres

    !================================================================

    SUBROUTINE compression_term_wrapper_stag(kk, jj, ii, component, splitDir, u, v, w, vff, dx, dy, dz, ddx, ddy, ddz, compressionTerm, propertyFluid1, propertyFluid2)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: component
        INTEGER(intk), INTENT(in) :: splitDir
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii), vff(kk, jj, ii, 3)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(out) :: compressionTerm(kk, jj, ii, 3)
        REAL(realk), INTENT(in), OPTIONAL :: propertyFluid1, propertyFluid2

        ! Local variables
        REAL(realk) :: strainRate(kk, jj, ii, 3)
        REAL(realk) :: nonDirectionalCompressionCoefficient(kk, jj, ii, 3)

        CALL compute_normal_strain_rate(kk, jj, ii, component, splitDir, strainRate, u, v, w, dx, dy, dz, ddx, ddy, ddz)
        CALL compute_non_directional_compression_coeffiecient(kk, jj, ii, component, nonDirectionalCompressionCoefficient, vff)

        IF ( .NOT. PRESENT(propertyFluid1) .OR. .NOT. PRESENT(propertyFluid2) ) THEN
            compressionTerm(:,:,:,component) = nonDirectionalCompressionCoefficient(:,:,:,component) * strainRate(:,:,:,component)
        ELSE IF ( PRESENT(propertyFluid1) .AND. PRESENT(propertyFluid2) ) THEN
            compressionTerm(:,:,:,component) = ( nonDirectionalCompressionCoefficient(:,:,:,component) * propertyFluid1 + ( 1.0 - nonDirectionalCompressionCoefficient(:,:,:,component) ) * propertyFluid2 ) * strainRate(:,:,:,component)
        END IF

    END SUBROUTINE compression_term_wrapper_stag

END MODULE multiphase_vof_transport_mod
