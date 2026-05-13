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
    USE multiphase_plic_mod, ONLY: track_interface, compute_normal_vector, compute_alpha, compute_cell_proportion, compute_iStag_vff, compute_jStag_vff, compute_kStag_vff, interface_reconstruction_wrapper, staggered_fractions_wrapper
    USE rungekutta_mod, ONLY: rk_2n_t
    USE multiphasecore_mod, ONLY: gmol1, gmol2, rho1, rho2
    USE multiphase_material_mod, ONLY: compute_material_property_field
    
    IMPLICIT NONE
    PRIVATE 

    PUBLIC :: init_multiphase_vof_transport, finish_multiphase_vof_transport, compute_flux, compute_density_flux, get_advection_sequence, compression_term_wrapper, update_field, clip_volume_fraction_field, multiphase_split_advection

    INTERFACE compute_normal_strain_rate
        MODULE PROCEDURE comp_normal_strain_rate_pres
        MODULE PROCEDURE comp_normal_strain_rate_stag
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

    SUBROUTINE comp_normal_strain_rate_pres(kk, jj, ii, splitDir, &
    u, v, w, ddx, ddy, ddz, strainRate)
    !----------------------------------------------------------------
    !   What it does:
    !   Computation of the normal strain rates dependend on the 
    !   directional split. 
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: splitDir
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(out) :: strainRate(kk, jj, ii)

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

    END SUBROUTINE comp_normal_strain_rate_pres

    !================================================================

    SUBROUTINE comp_normal_strain_rate_stag(kk, jj, ii, q, splitDir, &
        u, v, w, dx, dy, dz, ddx, ddy, ddz, strainRate)
    !----------------------------------------------------------------
    !   What it does:
    !   Computation of the normal strain rates dependend on the 
    !   directional split for the three staggered grids.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, q, splitDir
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(out) :: strainRate(kk, jj, ii, 3)

        ! Local variables
        INTEGER(intk) :: k, j, i
        REAL(realk) :: advrE(kk, jj, ii), advrW(kk, jj, ii)
        REAL(realk) :: advrN(kk, jj, ii), advrS(kk, jj, ii)
        REAL(realk) :: advrT(kk, jj, ii), advrB(kk, jj, ii)
        REAL(realk) :: iStag, jStag, kStag
        REAL(realk) :: deltaX(ii), deltaY(jj), deltaZ(kk)

        CALL get_component_specifics(kk, jj, ii, q, dx=dx, dy=dy, dz=dz, ddx=ddx, ddy=ddy, ddz=ddz, &
            iStag=iStag, jStag=jStag, kStag=kStag, deltaX=deltaX, deltaY=deltaY, deltaZ=deltaZ)
        CALL comp_advr_centr(kk, jj, ii, u, v, w, iStag, jStag, kStag, advrE, advrW, advrN, advrS, advrT, advrB)

        IF ( splitDir == 1 ) THEN
            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        strainRate(k,j,i,q) = ( advrE(k,j,i) - advrW(k,j,i) ) / deltaX(i)
                    END DO
                END DO
            END DO
        ELSE IF ( splitDir == 2 ) THEN
            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        strainRate(k,j,i,q) = ( advrN(k,j,i) - advrS(k,j,i) ) / deltaY(j)
                    END DO
                END DO
            END DO
        ELSE IF ( splitDir == 3 ) THEN 
            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        strainRate(k,j,i,q) = ( advrT(k,j,i) - advrB(k,j,i) ) / deltaZ(k)                 
                    END DO
                END DO
            END DO
        END IF

    END SUBROUTINE comp_normal_strain_rate_stag
    
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
                    IF ( vff(k,j,i) >= 0.5_realk ) THEN
                        nonDirectionalCompressionCoefficient(k,j,i) = 1.0_realk
                    ELSE
                        nonDirectionalCompressionCoefficient(k,j,i) = 0.0_realk
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
                    IF ( vff(k,j,i,component) >= 0.5_realk ) THEN
                        nonDirectionalCompressionCoefficient(k,j,i,component) = 1.0_realk
                    ELSE
                        nonDirectionalCompressionCoefficient(k,j,i,component) = 0.0_realk
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
                            IF ( isInterface(k,j,i+1) ) THEN
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
                            flux = 0.0_realk
                        END IF
                        fieldFlux(k,j,i) = sign( 1.0_realk, u(k,j,i) ) * flux / dt
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
                            IF ( isInterface(k,j+1,i) ) THEN
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
                            flux = 0.0_realk
                        END IF
                        fieldFlux(k,j,i) = sign( 1.0_realk, v(k,j,i) ) * flux / dt
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
                            IF ( isInterface(k+1,j,i) ) THEN
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
                            flux = 0.0_realk
                        END IF
                        fieldFlux(k,j,i) = sign( 1.0_realk, w(k,j,i) ) * flux / dt
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
        REAL(realk), INTENT(inout) :: fieldFlux(kk, jj, ii, 3), complementFieldFlux(kk, jj, ii, 3)
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
        REAL(realk) :: flux, complementFlux, fluxedProportion, eulerianFluxWidth, eulerianFluxAlpha
        REAL(realk) :: advrE(kk, jj, ii), advrW(kk, jj, ii), advrN(kk, jj, ii), advrS(kk, jj, ii), advrT(kk, jj, ii), advrB(kk, jj, ii)

        IF ( component == 1 ) THEN
            iStag = 1.0_realk
            jStag = 0.0_realk
            kStag = 0.0_realk
            deltaX = dx
            deltaY = ddy
            deltaZ = ddz
        ELSE IF ( component == 2 ) THEN
            iStag = 0.0_realk
            jStag = 1.0_realk
            kStag = 0.0_realk
            deltaX = ddx
            deltaY = dy
            deltaZ = ddz
        ELSE IF ( component == 3 ) THEN
            iStag = 0.0_realk
            jStag = 0.0_realk
            kStag = 1.0_realk
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

        CALL comp_advr_centr(kk, jj, ii, u, v, w, iStag, jStag, kStag, advrE, advrW, advrN, advrS, advrT, advrB)

        IF ( splitDir == 1 ) THEN
            DO i = 3-nfu, ii-3+nbu
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        IF ( advrE(k,j,i) > tol ) THEN
                            IF ( isInterface(k,j,i,component) ) THEN
                                ! Calculate the width of the fluxed volume and the alpha value for this subcell of the investigated cell
                                eulerianFluxWidth = abs( advrE(k,j,i) ) * dt
                                eulerianFluxAlpha = alpha(k,j,i,component) - normx(k,j,i,component) * ( deltaX(i) - eulerianFluxWidth )

                                ! Caluculate the volume fraction in the fluxed volume subcell of the investigated cell 
                                CALL compute_cell_proportion(fluxedProportion, eulerianFluxAlpha, field(k,j,i,component), eulerianFluxWidth, ddy(j), ddz(k), normx(k,j,i,component), normy(k,j,i,component), normz(k,j,i,component), tol)
                                
                                ! Calculate flux for multiphase cell
                                flux = fluxedProportion  * ( abs( advrE(k,j,i) ) * dt / deltaX(i) )
                                complementFlux = ( 1.0_realk - fluxedProportion ) * ( abs( advrE(k,j,i) ) * dt / deltaX(i) )
                            ELSE
                                ! Calculate flux for singlephase cell
                                flux = field(k,j,i,component) * ( abs( advrE(k,j,i) ) * dt / deltaX(i) )
                                complementFlux = ( 1.0_realk - field(k,j,i,component) ) * ( abs( advrE(k,j,i) ) * dt / deltaX(i) )
                            END IF
                        ELSE IF ( advrE(k,j,i) < -tol ) THEN
                            IF ( isInterface(k,j,i+1,component) ) THEN
                                ! Calculate the width of the fluxed volume and the alpha value for this subcell of the eastern cell
                                eulerianFluxWidth = abs( advrE(k,j,i) ) * dt
                                eulerianFluxAlpha = alpha(k,j,i+1,component)

                                ! Caluculate the volume fraction in the fluxed volume subcell of the eastern cell
                                CALL compute_cell_proportion(fluxedProportion, eulerianFluxAlpha, field(k,j,i+1,component), eulerianFluxWidth, ddy(j), ddz(k), normx(k,j,i+1,component), normy(k,j,i+1,component), normz(k,j,i+1,component), tol)

                                ! Calculate flux for multiphase cell
                                flux = fluxedProportion * ( abs( advrE(k,j,i) ) * dt / deltaX(i+1) )
                                complementFlux = ( 1.0_realk - fluxedProportion ) * ( abs( advrE(k,j,i) ) * dt / deltaX(i+1) )
                            ELSE
                                ! Calculate flux for singlephase cell
                                flux = field(k,j,i+1,component) * ( abs( advrE(k,j,i) ) * dt / deltaX(i+1) )
                                complementFlux = ( 1.0_realk - field(k,j,i+1,component) ) * ( abs( advrE(k,j,i) ) * dt / deltaX(i+1) )
                            END IF
                        ELSE
                            flux = 0.0_realk
                            complementFlux = 0.0_realk
                        END IF
                        fieldFlux(k,j,i,component) = sign( 1.0_realk, advrE(k,j,i) ) * flux / dt
                        complementFieldFlux(k,j,i,component) = sign( 1.0_realk, advrE(k,j,i) ) * complementFlux / dt
                    END DO
                END DO
            END DO
        ELSE IF ( splitDir == 2 ) THEN
            DO i = 3, ii-2
                DO j = 3-nrv, jj-3+nlv
                    DO k = 3, kk-2
                        IF ( advrN(k,j,i) > tol ) THEN
                            IF ( isInterface(k,j,i,component) ) THEN
                                ! Calculate the width of the fluxed volume and the alpha value for this subcell of the investigated cell
                                eulerianFluxWidth = abs( advrN(k,j,i) ) * dt
                                eulerianFluxAlpha = alpha(k,j,i,component) - normy(k,j,i,component) * ( deltaY(j) - eulerianFluxWidth )

                                ! Caluculate the volume fraction in the fluxed volume subcell of the investigated cell 
                                CALL compute_cell_proportion(fluxedProportion, eulerianFluxAlpha, field(k,j,i,component), ddx(i), eulerianFluxWidth, ddz(k), normx(k,j,i,component), normy(k,j,i,component), normz(k,j,i,component), tol)

                                ! Calculate flux for multiphase cell
                                flux = fluxedProportion * ( abs( advrN(k,j,i) ) * dt / deltaY(j) )
                                complementFlux = ( 1.0_realk - fluxedProportion ) * ( abs( advrN(k,j,i) ) * dt / deltaY(j) )
                            ELSE
                                ! Calculate flux for singlephase cell
                                flux = field(k,j,i,component) * ( abs( advrN(k,j,i) ) * dt / deltaY(j) )
                                complementFlux = ( 1.0_realk - field(k,j,i,component) ) * ( abs( advrN(k,j,i) ) * dt / deltaY(j) )
                            END IF
                        ELSE IF ( advrN(k,j,i) < -tol ) THEN
                            IF ( isInterface(k,j+1,i,component) ) THEN
                                ! Calculate the width of the fluxed volume and the alpha value for this subcell of the northern cell
                                eulerianFluxWidth = abs( advrN(k,j,i) ) * dt
                                eulerianFluxAlpha = alpha(k,j+1,i,component)

                                ! Caluculate the volume fraction in the fluxed volume subcell of the northern cell 
                                CALL compute_cell_proportion(fluxedProportion, eulerianFluxAlpha, field(k,j+1,i,component), ddx(i), eulerianFluxWidth, ddz(k), normx(k,j+1,i,component), normy(k,j+1,i,component), normz(k,j+1,i,component), tol)

                                ! Calculate flux for multiphase cell
                                flux = fluxedProportion * ( abs( advrN(k,j,i) ) * dt / deltaY(j+1) )
                                complementFlux = ( 1.0_realk - fluxedProportion ) * ( abs( advrN(k,j,i) ) * dt / deltaY(j+1) )
                            ELSE
                                ! Calculate flux for singlephase cell
                                flux = field(k,j+1,i,component) * ( abs( advrN(k,j,i) ) * dt / deltaY(j+1) )
                                complementFlux = ( 1.0_realk - field(k,j+1,i,component) ) * ( abs( advrN(k,j,i) ) * dt / deltaY(j+1) )
                            END IF
                        ELSE
                            flux = 0.0_realk
                            complementFlux = 0.0_realk
                        END IF
                        fieldFlux(k,j,i,component) = sign( 1.0_realk, advrN(k,j,i) ) * flux / dt
                        complementFieldFlux(k,j,i,component) = sign( 1.0_realk, advrN(k,j,i) ) * complementFlux / dt
                    END DO
                END DO
            END DO
        ELSE IF ( splitDir == 3 ) THEN
            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3-nbw, kk-3+ntw
                        IF ( advrT(k,j,i) > tol ) THEN
                            IF ( isInterface(k,j,i,component) ) THEN
                                ! Calculate the width of the fluxed volume and the alpha value for this subcell of the investigated cell
                                eulerianFluxWidth = abs( advrT(k,j,i) ) * dt
                                eulerianFluxAlpha = alpha(k,j,i,component) - normz(k,j,i,component) * ( deltaZ(k) - eulerianFluxWidth )

                                ! Caluculate the volume fraction in the fluxed volume subcell of the investigated cell 
                                CALL compute_cell_proportion(fluxedProportion, eulerianFluxAlpha, field(k,j,i,component), ddx(i), ddy(j), eulerianFluxWidth, normx(k,j,i,component), normy(k,j,i,component), normz(k,j,i,component), tol)

                                ! Calculate flux for multiphase cell
                                flux = fluxedProportion * ( abs( advrT(k,j,i) ) * dt / deltaZ(k) )
                                complementFlux = ( 1.0_realk - fluxedProportion ) * ( abs( advrT(k,j,i) ) * dt / deltaZ(k) )
                            ELSE
                                ! Calculate flux for singlephase cell
                                flux = field(k,j,i,component) * ( abs( advrT(k,j,i) ) * dt / deltaZ(k) )
                                complementFlux = ( 1.0_realk - field(k,j,i,component) ) * ( abs( advrT(k,j,i) ) * dt / deltaZ(k) )
                            END IF
                        ELSE IF ( advrT(k,j,i) < -tol ) THEN
                            IF ( isInterface(k+1,j,i,component) ) THEN
                                ! Calculate the width of the fluxed volume and the alpha value for this subcell of the investigated cell
                                eulerianFluxWidth = abs( advrT(k,j,i) ) * dt
                                eulerianFluxAlpha = alpha(k+1,j,i,component)

                                ! Caluculate the volume fraction in the fluxed volume subcell of the investigated cell 
                                CALL compute_cell_proportion(fluxedProportion, eulerianFluxAlpha, field(k+1,j,i,component), ddx(i), ddy(j), eulerianFluxWidth, normx(k+1,j,i,component), normy(k+1,j,i,component), normz(k+1,j,i,component), tol)

                                ! Calculate flux for multiphase cell
                                flux = fluxedProportion * ( abs( advrT(k,j,i) ) * dt / deltaZ(k+1) )
                                complementFlux = ( 1.0_realk - fluxedProportion ) * ( abs( advrT(k,j,i) ) * dt / deltaZ(k+1) )
                            ELSE
                                ! Calculate flux for singlephase cell
                                flux = field(k+1,j,i,component) * ( abs( advrT(k,j,i) ) * dt / deltaZ(k+1) )
                                complementFlux = ( 1.0_realk - field(k+1,j,i,component) ) * ( abs( advrT(k,j,i) ) * dt / deltaZ(k+1) )
                            END IF
                        ELSE
                            flux = 0.0_realk
                            complementFLux = 0.0_realk
                        END IF
                        fieldFlux(k,j,i,component) = sign( 1.0_realk, advrT(k,j,i) ) * flux / dt
                        complementFieldFlux(k,j,i,component) = sign( 1.0_realk, advrT(k,j,i) ) * complementFlux / dt
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
                        vff(k,j,i) = 0.0_realk
                    ELSE IF ( vff(k,j,i) > 1.0_realk - tol ) THEN
                        vff(k,j,i) = 1.0_realk
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
        REAL(realk), INTENT(inout) :: densityFlux(kk, jj, ii, 3)

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

        CALL compute_normal_strain_rate(kk, jj, ii, splitDir, u, v, w, ddx, ddy, ddz, strainRate)
        CALL compute_non_directional_compression_coeffiecient(kk, jj, ii, nonDirectionalCompressionCoefficient, vff)

        IF ( .NOT. PRESENT(propertyFluid1) .AND. .NOT. PRESENT(propertyFluid2) ) THEN
            compressionTerm = nonDirectionalCompressionCoefficient * strainRate
        ELSE IF ( PRESENT(propertyFluid1) .AND. PRESENT(propertyFluid2) ) THEN
            compressionTerm = ( nonDirectionalCompressionCoefficient * propertyFluid1 + ( 1.0_realk - nonDirectionalCompressionCoefficient ) * propertyFluid2 ) * strainRate
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

        CALL compute_normal_strain_rate(kk, jj, ii, component, splitDir, u, v, w, dx, dy, dz, ddx, ddy, ddz, strainRate)
        CALL compute_non_directional_compression_coeffiecient(kk, jj, ii, component, nonDirectionalCompressionCoefficient, vff)

        IF ( .NOT. PRESENT(propertyFluid1) .OR. .NOT. PRESENT(propertyFluid2) ) THEN
            compressionTerm(:,:,:,component) = nonDirectionalCompressionCoefficient(:,:,:,component) * strainRate(:,:,:,component)
        ELSE IF ( PRESENT(propertyFluid1) .AND. PRESENT(propertyFluid2) ) THEN
            compressionTerm(:,:,:,component) = ( nonDirectionalCompressionCoefficient(:,:,:,component) * propertyFluid1 + ( 1.0_realk - nonDirectionalCompressionCoefficient(:,:,:,component) ) * propertyFluid2 ) * strainRate(:,:,:,component)
        END IF

    END SUBROUTINE compression_term_wrapper_stag

    !================================================================

    SUBROUTINE multiphase_split_advection(uo_f, vo_f, wo_f, u_f, v_f, w_f, ut_f, vt_f, wt_f, &
        vff_f, p_f, g_f, d_f, dtrki, itstep)
    !----------------------------------------------------------------
    !   What it does:
    !    
    !----------------------------------------------------------------

        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: uo_f
        TYPE(field_t), INTENT(inout) :: vo_f
        TYPE(field_t), INTENT(inout) :: wo_f
        TYPE(field_t), INTENT(in) :: u_f
        TYPE(field_t), INTENT(in) :: v_f
        TYPE(field_t), INTENT(in) :: w_f
        TYPE(field_t), INTENT(in) :: ut_f
        TYPE(field_t), INTENT(in) :: vt_f
        TYPE(field_t), INTENT(in) :: wt_f
        TYPE(field_t), INTENT(in) :: vff_f
        TYPE(field_t), INTENT(in) :: p_f
        TYPE(field_t), INTENT(in) :: g_f
        TYPE(field_t), INTENT(in) :: d_f
        REAL(realk), INTENT(in) :: dtrki
        INTEGER(intk), INTENT(in) :: itstep

        ! Local variables
        TYPE(field_t), POINTER :: dx_f, dy_f, dz_f, ddx_f, ddy_f, ddz_f
        TYPE(field_t), POINTER :: rdx_f, rdy_f, rdz_f, rddx_f, rddy_f, rddz_f
        TYPE(field_t), POINTER :: normx_f, normy_f, normz_f, alpha_f
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: uo, vo, wo
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: u, v, w
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: ut, vt, wt
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: vff, p, g, d
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        REAL(realk), POINTER, CONTIGUOUS :: rdx(:), rdy(:), rdz(:)
        REAL(realk), POINTER, CONTIGUOUS :: rddx(:), rddy(:), rddz(:)
        REAL(realk), POINTER, CONTIGUOUS :: normxField(:,:,:), normyField(:,:,:), normzField(:,:,:), alphaField(:,:,:)
        INTEGER(intk) :: advSeq(3)
        INTEGER(intk) :: i, igrid, l, q, splitDir
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: nfro, nbac, nrgt, nlft, nbot, ntop
        REAL(realk), PARAMETER :: tol = 1.0E-15_realk
        REAL(realk), ALLOCATABLE :: normx(:,:,:), normy(:,:,:), normz(:,:,:)
        REAL(realk), ALLOCATABLE :: normxStag(:,:,:,:), normyStag(:,:,:,:), normzStag(:,:,:,:)
        REAL(realk), ALLOCATABLE :: alpha(:,:,:), alphaStag(:,:,:,:)
        LOGICAL, ALLOCATABLE :: isInterface(:,:,:), isInterfaceStag(:,:,:,:)
        LOGICAL, ALLOCATABLE :: isNearInterface(:,:,:), isNearInterfaceStag(:,:,:,:)
        REAL(realk), ALLOCATABLE :: vffStag(:,:,:,:), vffOld(:,:,:), VffStagOld(:,:,:,:)
        REAL(realk), ALLOCATABLE :: densityFieldStag(:,:,:,:), densityFieldStagOld(:,:,:,:)
        REAL(realk), ALLOCATABLE :: vffFlux(:,:,:), vffFluxStag(:,:,:,:)
        REAL(realk), ALLOCATABLE :: complementvffFlux(:,:,:), complementvffFluxStag(:,:,:,:)
        REAL(realk), ALLOCATABLE :: densityFieldFluxStag(:,:,:,:)
        REAL(realk), ALLOCATABLE :: vffCompressionTerm(:,:,:)
        REAL(realk), ALLOCATABLE :: densityCompressionTermStag(:,:,:,:)
        REAL(realk), ALLOCATABLE :: uNew(:,:,:), vNew(:,:,:), wNew(:,:,:)

        ! Set all the output to zero everywhere before we start!
        uo_f = 0.0_realk
        vo_f = 0.0_realk
        wo_f = 0.0_realk

        CALL get_field(dx_f, "DX")
        CALL get_field(dy_f, "DY")
        CALL get_field(dz_f, "DZ")

        CALL get_field(ddx_f, "DDX")
        CALL get_field(ddy_f, "DDY")
        CALL get_field(ddz_f, "DDZ")

        CALL get_field(rdx_f, "RDX")
        CALL get_field(rdy_f, "RDY")
        CALL get_field(rdz_f, "RDZ")

        CALL get_field(rddx_f, "RDDX")
        CALL get_field(rddy_f, "RDDY")
        CALL get_field(rddz_f, "RDDZ")

        CALL get_field(normx_f, "NORMX")
        CALL get_field(normy_f, "NORMY")
        CALL get_field(normz_f, "NORMZ")

        CALL get_field(alpha_f, "ALPHA")

        DO i = 1, nmygrids
            igrid = mygrids(i)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_mgbasb(nfro, nbac, nrgt, nlft, nbot, ntop, igrid)

            CALL uo_f%get_ptr(uo, igrid)
            CALL vo_f%get_ptr(vo, igrid)
            CALL wo_f%get_ptr(wo, igrid)

            CALL u_f%get_ptr(u, igrid)
            CALL v_f%get_ptr(v, igrid)
            CALL w_f%get_ptr(w, igrid)

            CALL ut_f%get_ptr(ut, igrid)
            CALL vt_f%get_ptr(vt, igrid)
            CALL wt_f%get_ptr(wt, igrid)

            CALL vff_f%get_ptr(vff, igrid)
            CALL p_f%get_ptr(p, igrid)
            CALL g_f%get_ptr(g, igrid)
            CALL d_f%get_ptr(d, igrid)

            CALL dx_f%get_ptr(dx, igrid)
            CALL dy_f%get_ptr(dy, igrid)
            CALL dz_f%get_ptr(dz, igrid)

            CALL ddx_f%get_ptr(ddx, igrid)
            CALL ddy_f%get_ptr(ddy, igrid)
            CALL ddz_f%get_ptr(ddz, igrid)

            CALL rdx_f%get_ptr(rdx, igrid)
            CALL rdy_f%get_ptr(rdy, igrid)
            CALL rdz_f%get_ptr(rdz, igrid)

            CALL normx_f%get_ptr(normxField, igrid)
            CALL normy_f%get_ptr(normyField, igrid)
            CALL normz_f%get_ptr(normzField, igrid)
            CALL alpha_f%get_ptr(alphaField, igrid)

            CALL rddx_f%get_ptr(rddx, igrid)
            CALL rddy_f%get_ptr(rddy, igrid)
            CALL rddz_f%get_ptr(rddz, igrid)

            IF (.NOT. ALLOCATED(normx))     ALLOCATE(normx(kk,jj,ii))
            IF (.NOT. ALLOCATED(normy))     ALLOCATE(normy(kk,jj,ii))
            IF (.NOT. ALLOCATED(normz))     ALLOCATE(normz(kk,jj,ii))

            IF (.NOT. ALLOCATED(normxStag)) ALLOCATE(normxStag(kk,jj,ii,3))
            IF (.NOT. ALLOCATED(normyStag)) ALLOCATE(normyStag(kk,jj,ii,3))
            IF (.NOT. ALLOCATED(normzStag)) ALLOCATE(normzStag(kk,jj,ii,3))

            IF (.NOT. ALLOCATED(alpha))     ALLOCATE(alpha(kk,jj,ii))
            IF (.NOT. ALLOCATED(alphaStag)) ALLOCATE(alphaStag(kk,jj,ii,3))

            IF (.NOT. ALLOCATED(isInterface))     ALLOCATE(isInterface(kk,jj,ii))
            IF (.NOT. ALLOCATED(isInterfaceStag)) ALLOCATE(isInterfaceStag(kk,jj,ii,3))

            IF (.NOT. ALLOCATED(isNearInterface))     ALLOCATE(isNearInterface(kk,jj,ii))
            IF (.NOT. ALLOCATED(isNearInterfaceStag)) ALLOCATE(isNearInterfaceStag(kk,jj,ii,3))

            IF (.NOT. ALLOCATED(vffStag)) ALLOCATE(vffStag(kk,jj,ii,3))
            IF (.NOT. ALLOCATED(vffStagOld)) ALLOCATE(vffStagOld(kk,jj,ii,3))
            IF (.NOT. ALLOCATED(vffOld)) ALLOCATE(vffOld(kk,jj,ii))

            IF (.NOT. ALLOCATED(densityFieldStag)) ALLOCATE(densityFieldStag(kk,jj,ii,3))
            IF (.NOT. ALLOCATED(densityFieldStagOld)) ALLOCATE(densityFieldStagOld(kk,jj,ii,3))

            IF (.NOT. ALLOCATED(vffFlux)) ALLOCATE(vffFlux(kk,jj,ii))
            IF (.NOT. ALLOCATED(vffFluxStag)) ALLOCATE(vffFluxStag(kk,jj,ii,3))
            
            IF (.NOT. ALLOCATED(complementvffFlux)) ALLOCATE(complementvffFlux(kk,jj,ii))
            IF (.NOT. ALLOCATED(complementvffFluxStag)) ALLOCATE(complementvffFluxStag(kk,jj,ii,3))

            IF (.NOT. ALLOCATED(densityFieldFluxStag)) ALLOCATE(densityFieldFluxStag(kk,jj,ii,3))

            IF (.NOT. ALLOCATED(vffCompressionTerm)) ALLOCATE(vffCompressionTerm(kk,jj,ii))

            IF (.NOT. ALLOCATED(densityCompressionTermStag)) ALLOCATE(densityCompressionTermStag(kk,jj,ii,3))

            IF (.NOT. ALLOCATED(uNew)) ALLOCATE(uNew(kk,jj,ii))
            IF (.NOT. ALLOCATED(vNew)) ALLOCATE(vNew(kk,jj,ii))
            IF (.NOT. ALLOCATED(wNew)) ALLOCATE(wNew(kk,jj,ii))

            CALL get_advection_sequence(itstep, advSeq)
            CALL interface_reconstruction_wrapper(kk, jj, ii, vff, ddx, ddy, ddz, tol, normx, normy, normz, alpha, isInterface)

            DO q = 1, 3             ! Loop over staggered components u, v and w
                
                CALL staggered_fractions_wrapper(kk, jj, ii, q, alpha, vff, isInterface, ddx, ddy, ddz, normx, normy, normz, tol, vffStag)
                CALL interface_reconstruction_wrapper(kk, jj, ii, q, vffStag, dx, dy, dz, ddx, ddy, ddz, tol, normxStag, normyStag, normzStag, alphaStag, isInterfaceStag, isNearInterfaceStag)
                CALL compute_material_property_field(kk, jj, ii, q, densityFieldStag, vffStag, rho1, rho2)

            END DO

            densityFieldStagOld = densityFieldStag
            vffOld = vff
            vffStagOld = vffStag

            DO l = 1, 3             ! Loop over dimensions x, y and z for split-advection

                splitDir = advSeq(l)

                CALL interface_reconstruction_wrapper(kk, jj, ii, vff, ddx, ddy, ddz, tol, normx, normy, normz, alpha, isInterface)
                CALL compute_flux(kk, jj, ii, splitDir, vffFlux, vff, isInterface, u, v, w, alpha, dtrki, normx, normy, normz, ddx, ddy, ddz, tol, nfro, nbac, nrgt, nlft, nbot, ntop)                
                CALL compression_term_wrapper(kk, jj, ii, splitDir, u, v, w, vffOld, ddx, ddy, ddz, vffCompressionTerm)
                CALL update_field(kk, jj, ii, splitDir, vff, vffFlux, vffCompressionTerm, dtrki, nfro, nbac, nrgt, nlft, nbot, ntop)
                CALL clip_volume_fraction_field(kk, ii, jj, vff, tol) 

                DO q = 1, 3         ! Loop over staggered components u, v and w
                    
                    CALL compute_flux(kk, jj, ii, q, splitDir, vffFluxStag, complementvffFluxStag, vffStag, isInterfaceStag, u, v, w, alphaStag, dtrki, normxStag, normyStag, normzStag, dx, dy, dz, ddx, ddy, ddz, tol, nfro, nbac, nrgt, nlft, nbot, ntop)
                    CALL compute_density_flux(kk, jj, ii, q, vffFluxStag, complementvffFluxStag, rho1, rho2, densityFieldFluxStag)
                    CALL compression_term_wrapper(kk, jj, ii, q, splitDir, u, v, w, vffStagOld, dx, dy, dz, ddx, ddy, ddz, densityCompressionTermStag, rho1, rho2)
                    CALL update_field(kk, jj, ii, q, splitDir, densityFieldStag, densityFieldFluxStag, densityCompressionTermStag, dtrki, nfro, nbac, nrgt, nlft, nbot, ntop)
                    CALL multiphase_advect_momentum(kk, jj, ii, q, splitDir, u, v, w, uNew, vNew, wNew, densityFieldFluxStag, &
                        densityCompressionTermStag, densityFieldStagOld, densityFieldStag, &
                        isNearInterfaceStag, dx, dy, dz, ddx, ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz, &
                        nfro, nbac, nrgt, nlft, nbot, ntop)

                END DO
            END DO
            u(3:-3,3:-3,3:-3) = uNew(3:-3,3:-3,3:-3)
            v(3:-3,3:-3,3:-3) = vNew(3:-3,3:-3,3:-3)
            w(3:-3,3:-3,3:-3) = wNew(3:-3,3:-3,3:-3)
        END DO

        normxField = normx
        normyField = normy
        normzField = normz
        alphaField = alpha

    END SUBROUTINE multiphase_split_advection

    !================================================================

    SUBROUTINE multiphase_advect_momentum(kk, jj, ii, component, splitDir, u, v, w, uNew, vNew, wNew, densityFieldFluxStag, &
        densityCompressionTermStag, densityFieldStagOld, densityFieldStag, &
        isNearInterfaceStag, dx, dy, dz, ddx, ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz, &
        nfro, nbac, nrgt, nlft, nbot, ntop)
    !----------------------------------------------------------------
    !   What it does:
    !    
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: splitDir
        INTEGER(intk), INTENT(in) :: component
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(out) :: uNew(kk, jj, ii), vNew(kk, jj, ii), wNew(kk, jj, ii)
        REAL(realk), INTENT(in) :: densityFieldStagOld(kk, jj, ii, 3), densityFieldFluxStag(kk, jj, ii, 3) 
        REAL(realk), INTENT(in) :: densityCompressionTermStag(kk, jj, ii, 3)
        REAL(realk), INTENT(in) :: densityFieldStag(kk, jj, ii, 3)
        LOGICAL, INTENT(in) :: isNearInterfaceStag(kk, jj, ii, 3)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: rdx(ii), rdy(jj), rdz(kk)
        REAL(realk), INTENT(in) :: rddx(ii), rddy(jj), rddz(kk)
        INTEGER, INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: nbu, nfu, nrv, nlv, nbw, ntw
        REAL(realk) :: advrE(kk, jj, ii), advrW(kk, jj, ii), advrN(kk, jj, ii), advrS(kk, jj, ii), advrT(kk, jj, ii), advrB(kk, jj, ii)
        REAL(realk) :: adveE(kk, jj, ii), adveW(kk, jj, ii), adveN(kk, jj, ii), adveS(kk, jj, ii), adveT(kk, jj, ii), adveB(kk, jj, ii)
        REAL(realk) :: iStag, jStag, kStag
        REAL(realk) :: velocity(kk, jj, ii)
        REAL(realk) :: dVelocity, dMomentum

        ! return
        ! WRITE(*,*) "Mom. Adv. Running"

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

        IF ( component == 1 ) THEN
            iStag = 1.0_realk
            jStag = 0.0_realk
            kStag = 0.0_realk
            velocity = u
        ELSE IF ( component == 2 ) THEN
            iStag = 0.0_realk
            jStag = 1.0_realk
            kStag = 0.0_realk
            velocity = v
        ELSE IF ( component == 3 ) THEN
            iStag = 0.0_realk
            jStag = 0.0_realk
            kStag = 1.0_realk
            velocity = w
        END IF

        CALL comp_advr_centr(kk, jj, ii, u, v, w, iStag, jStag, kStag, advrE, advrW, advrN, advrS, advrT, advrB)
        CALL comp_adve_quick(kk, jj, ii, velocity, advrE, advrW, advrN, advrS, advrT, advrB, adveE, adveW, adveN, adveS, adveT, adveB)

        IF ( splitDir == 1 ) THEN
            DO i = 3-nfu, ii-3+nbu
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        IF ( isNearInterfaceStag(k,j,i,component) ) THEN
                            dMomentum = - ( adveE(k,j,i) * densityFieldFluxStag(k,j,i,component) - adveW(k,j,i) * densityFieldFluxStag(k,j,i-1,component) ) + velocity(k,j,i) * densityCompressionTermStag(k,j,i,component)
                            velocity(k,j,i) = 1.0_realk / densityFieldStag(k,j,i,component) * ( densityFieldStagOld(k,j,i,component) * velocity(k,j,i) + dMomentum )
                        ELSE
                            dVelocity = - ( ( adveE(k,j,i) * advrE(k,j,i) - adveW(k,j,i) * advrW(k,j,i) ) * rdx(i) )
                            velocity(k,j,i) = velocity(k,j,i) + dVelocity
                        END IF
                    END DO
                END DO
            END DO
        ELSE IF ( splitDir == 2 ) THEN
            DO i = 3, ii-2
                DO j = 3-nrv, jj-3+nlv
                    DO k = 3, kk-2
                        IF ( isNearInterfaceStag(k,j,i,component) ) THEN
                            dMomentum = - ( adveN(k,j,i) * densityFieldFluxStag(k,j,i,component) - adveS(k,j,i) * densityFieldFluxStag(k,j-1,i,component) ) + velocity(k,j,i) * densityCompressionTermStag(k,j,i,component)
                            velocity(k,j,i) = 1.0_realk / densityFieldStag(k,j,i,component) * ( densityFieldStagOld(k,j,i,component) * velocity(k,j,i) + dMomentum )
                        ELSE
                            dVelocity = - ( ( adveN(k,j,i) * advrN(k,j,i) - adveS(k,j,i) * advrS(k,j,i) ) * rdy(j) )
                            velocity(k,j,i) = velocity(k,j,i) + dVelocity
                        END IF
                    END DO
                END DO
            END DO
        ELSE IF ( splitDir == 3 ) THEN
            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3-nbw, kk-3+ntw
                        IF ( isNearInterfaceStag(k,j,i,component) ) THEN
                            dMomentum = - ( adveT(k,j,i) * densityFieldFluxStag(k,j,i,component) - adveB(k,j,i) * densityFieldFluxStag(k-1,j,i,component) ) + velocity(k,j,i) * densityCompressionTermStag(k,j,i,component)
                            velocity(k,j,i) = 1.0_realk / densityFieldStag(k,j,i,component) * ( densityFieldStagOld(k,j,i,component) * velocity(k,j,i) + dMomentum )
                        ELSE
                            dVelocity = - ( ( adveT(k,j,i) * advrT(k,j,i) - adveB(k,j,i) * advrB(k,j,i) ) * rdz(k) )
                            velocity(k,j,i) = velocity(k,j,i) + dVelocity
                        END IF
                    END DO
                END DO
            END DO
        END IF

        IF ( component == 1 ) THEN
            uNew = velocity
        ELSE IF ( component == 2 ) THEN
            vNew = velocity
        ELSE IF ( component == 3 ) THEN
            wNew = velocity
        END IF

    END SUBROUTINE multiphase_advect_momentum

    !================================================================

    PURE SUBROUTINE comp_adve_quick(kk, jj, ii, adveField, &
        advrE, advrW, advrN, advrS, advrT, advrB, &
        adveE, adveW, adveN, adveS, adveT, adveB)
    !----------------------------------------------------------------
    !   What it does:
    !   QUICK interpolation to compute the advected veloctiy
    !   (advectee) on staggered grid cells.
    !   adve = advected component (advectee)
    !   advr = advecting component (advector)
    !   An indicator function is used to avoid if-statements within
    !   loops.
    !   
    !   Source: 
    !   T. Arrufat et al., “A mass-momentum consistent, 
    !   Volume-of-Fluid method for incompressible flow on staggered 
    !   grids,” Computers & Fluids, vol. 215, p. 104785, Jan. 2021, 
    !   doi: 10.1016/j.compfluid.2020.104785.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: adveField(kk, jj, ii)
        REAL(realk), INTENT(in) :: advrE(kk, jj, ii), advrW(kk, jj, ii)
        REAL(realk), INTENT(in) :: advrN(kk, jj, ii), advrS(kk, jj, ii)
        REAL(realk), INTENT(in) :: advrT(kk, jj, ii), advrB(kk, jj, ii)
        REAL(realk), INTENT(out) :: adveE(kk, jj, ii), adveW(kk, jj, ii)
        REAL(realk), INTENT(out) :: adveN(kk, jj, ii), adveS(kk, jj, ii)
        REAL(realk), INTENT(out) :: adveT(kk, jj, ii), adveB(kk, jj, ii)

        ! Loval variables
        INTEGER(intk) :: k, j, i
        REAL(realk) :: indE(2), indW(2), indN(2), indS(2), indT(2), indB(2)

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    indE(1) = merge(1.0_realk, 0.0_realk, advrE(k,j,i) >= 0.0_realk)
                    indW(1) = merge(1.0_realk, 0.0_realk, advrW(k,j,i) >= 0.0_realk)
                    indN(1) = merge(1.0_realk, 0.0_realk, advrN(k,j,i) >= 0.0_realk)
                    indS(1) = merge(1.0_realk, 0.0_realk, advrS(k,j,i) >= 0.0_realk)
                    indT(1) = merge(1.0_realk, 0.0_realk, advrT(k,j,i) >= 0.0_realk)
                    indB(1) = merge(1.0_realk, 0.0_realk, advrB(k,j,i) >= 0.0_realk)

                    indE(2) = 1.0_realk - indE(1)
                    indW(2) = 1.0_realk - indW(1)
                    indN(2) = 1.0_realk - indN(1)
                    indS(2) = 1.0_realk - indS(1)
                    indT(2) = 1.0_realk - indT(1)
                    indB(2) = 1.0_realk - indB(1)

                    adveE(k,j,i) = indE(1) * ( 0.75_realk * adveField(k,j,i) + 0.375_realk * adveField(k,j,i+1) - 0.125_realk * adveField(k,j,i-1) ) + &
                                   indE(2) * ( 0.75_realk * adveField(k,j,i+1) + 0.375_realk * adveField(k,j,i) - 0.125_realk * adveField(k,j,i+2) )
                    adveW(k,j,i) = indW(1) * ( 0.75_realk * adveField(k,j,i-1) + 0.375_realk * adveField(k,j,i) - 0.125_realk * adveField(k,j,i-2) ) + &
                                   indW(2) * ( 0.75_realk * adveField(k,j,i) + 0.375_realk * adveField(k,j,i-1) - 0.125_realk * adveField(k,j,i+1) )
                    adveN(k,j,i) = indN(1) * ( 0.75_realk * adveField(k,j,i) + 0.375_realk * adveField(k,j+1,i) - 0.125_realk * adveField(k,j-1,i) ) + &
                                   indN(2) * ( 0.75_realk * adveField(k,j+1,i) + 0.375_realk * adveField(k,j,i) - 0.125_realk * adveField(k,j+2,i) )
                    adveS(k,j,i) = indS(1) * ( 0.75_realk * adveField(k,j-1,i) + 0.375_realk * adveField(k,j,i) - 0.125_realk * adveField(k,j-2,i) ) + &
                                   indS(2) * ( 0.75_realk * adveField(k,j,i) + 0.375_realk * adveField(k,j-1,i) - 0.125_realk * adveField(k,j+1,i) )
                    adveT(k,j,i) = indT(1) * ( 0.75_realk * adveField(k,j,i) + 0.375_realk * adveField(k+1,j,i) - 0.125_realk * adveField(k-1,j,i) ) + &
                                   indT(2) * ( 0.75_realk * adveField(k+1,j,i) + 0.375_realk * adveField(k,j,i) - 0.125_realk * adveField(k+2,j,i) )
                    adveB(k,j,i) = indB(1) * ( 0.75_realk * adveField(k-1,j,i) + 0.375_realk * adveField(k,j,i) - 0.125_realk * adveField(k-2,j,i) ) + &
                                   indB(2) * ( 0.75_realk * adveField(k,j,i) + 0.375_realk * adveField(k-1,j,i) - 0.125_realk * adveField(k+1,j,i) )
                END DO
            END DO
        END DO

    END SUBROUTINE comp_adve_quick

    !================================================================

    PURE SUBROUTINE comp_advr_centr(kk, jj, ii, u, v, w, &
        iStag, jStag, kStag, advrE, advrW, advrN, advrS, advrT, advrB)
    !----------------------------------------------------------------
    !   What it does:
    !   Central difference to compute the advecting velocity
    !   (advector) on staggered grid cells.
    !   advr = advecting component (advector)
    !   
    !   Source: 
    !   T. Arrufat et al., “A mass-momentum consistent, 
    !   Volume-of-Fluid method for incompressible flow on staggered 
    !   grids,” Computers & Fluids, vol. 215, p. 104785, Jan. 2021, 
    !   doi: 10.1016/j.compfluid.2020.104785.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: iStag, jStag, kStag
        REAL(realk), INTENT(out) :: advrE(kk, jj, ii), advrW(kk, jj, ii)
        REAL(realk), INTENT(out) :: advrN(kk, jj, ii), advrS(kk, jj, ii)
        REAL(realk), INTENT(out) :: advrT(kk, jj, ii), advrB(kk, jj, ii)

        ! Loval variables
        INTEGER(intk) :: k, j, i

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    advrE(k,j,i) = 0.5_realk * ( iStag * ( u(k,j,i) + u(k,j,i+1) ) + &
                                                 jStag * ( u(k,j,i) + u(k,j+1,i) ) + &
                                                 kStag * ( u(k,j,i) + u(k+1,j,i) ) )
                    advrW(k,j,i) = 0.5_realk * ( iStag * ( u(k,j,i-1) + u(k,j,i) ) + &
                                                 jStag * ( u(k,j,i-1) + u(k,j+1,i-1) ) + &
                                                 kStag * ( u(k,j,i-1) + u(k+1,j,i-1) ) )
                    advrN(k,j,i) = 0.5_realk * ( iStag * ( v(k,j,i) + v(k,j,i+1) ) + &
                                                 jStag * ( v(k,j,i) + v(k,j+1,i) ) + &
                                                 kStag * ( v(k,j,i) + v(k+1,j,i) ) )
                    advrS(k,j,i) = 0.5_realk * ( iStag * ( v(k,j-1,i) + v(k,j-1,i+1) ) + &
                                                 jStag * ( v(k,j-1,i) + v(k,j,i) ) + &
                                                 kStag * ( v(k,j-1,i) + v(k+1,j-1,i) ) )
                    advrT(k,j,i) = 0.5_realk * ( iStag * ( w(k,j,i) + w(k,j,i+1) ) + &
                                                 jStag * ( w(k,j,i) + w(k,j+1,i) ) + &
                                                 kStag * ( w(k,j,i) + w(k+1,j,i) ) )
                    advrB(k,j,i) = 0.5_realk * ( iStag * ( w(k-1,j,i) + w(k-1,j,i+1) ) + &
                                                 jStag * ( w(k-1,j,i) + w(k-1,j+1,i) ) + &
                                                 kStag * ( w(k-1,j,i) + w(k,j,i) ) )
                END DO
            END DO
        END DO

    END SUBROUTINE comp_advr_centr

    !================================================================

    SUBROUTINE get_component_specifics(kk, jj, ii, q, u, v, w, &
        dx, dy, dz, ddx, ddy, ddz, iStag, jStag, kStag, &
        deltaX, deltaY, deltaZ, veloFld)
    !----------------------------------------------------------------
    !   What it does:
    !   Select the specifics for the qth staggered grid.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, q
        REAL(realk), INTENT(in), OPTIONAL :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in), OPTIONAL :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in), OPTIONAL :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(out), OPTIONAL :: iStag, jStag, kStag
        REAL(realk), INTENT(out), OPTIONAL :: deltaX(ii), deltaY(jj), deltaZ(kk)
        REAL(realk), INTENT(out), OPTIONAL :: veloFld(kk, jj, ii)

        ! Local variables
        ! None

        IF ( q == 1 ) THEN
            IF ( PRESENT(iStag) ) iStag = 1.0_realk
            IF ( PRESENT(jStag) ) jStag = 0.0_realk
            IF ( PRESENT(kStag) ) kStag = 0.0_realk
            IF ( PRESENT(deltaX) ) deltaX = dx
            IF ( PRESENT(deltaY) ) deltaY = ddy
            IF ( PRESENT(deltaZ) ) deltaZ = ddz
            IF ( PRESENT(veloFld) ) veloFld = u
        ELSE IF ( q == 2 ) THEN
            IF ( PRESENT(iStag) ) iStag = 0.0_realk
            IF ( PRESENT(jStag) ) jStag = 1.0_realk
            IF ( PRESENT(kStag) ) kStag = 0.0_realk
            IF ( PRESENT(deltaX) ) deltaX = ddx
            IF ( PRESENT(deltaY) ) deltaY = dy
            IF ( PRESENT(deltaZ) ) deltaZ = ddz
            IF ( PRESENT(veloFld) ) veloFld = v
        ELSE IF ( q == 3 ) THEN
            IF ( PRESENT(iStag) ) iStag = 0.0_realk
            IF ( PRESENT(jStag) ) jStag = 0.0_realk
            IF ( PRESENT(kStag) ) kStag = 1.0_realk
            IF ( PRESENT(deltaX) ) deltaX = ddx
            IF ( PRESENT(deltaY) ) deltaY = ddy
            IF ( PRESENT(deltaZ) ) deltaZ = dz
            IF ( PRESENT(veloFld) ) veloFld = w
        END IF
    
    END SUBROUTINE get_component_specifics

END MODULE multiphase_vof_transport_mod
