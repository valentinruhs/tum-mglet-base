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
    USE multiphase_material_mod, ONLY: comp_material_property_field
    USE multiphasecore_mod, ONLY: mom_multiphase
    
    IMPLICIT NONE
    PRIVATE 

    PUBLIC :: init_multiphase_vof_transport, finish_multiphase_vof_transport, multiphase_solve

    ! INTERFACE compute_normal_strain_rate
    !     MODULE PROCEDURE comp_normal_strain_rate_pres
    !     MODULE PROCEDURE comp_normal_strain_rate_stag
    ! END INTERFACE

    ! INTERFACE update_field
    !     MODULE PROCEDURE update_field_pres
    !     MODULE PROCEDURE update_field_stag
    ! END INTERFACE

    INTERFACE comp_flux
        MODULE PROCEDURE comp_flux_pres
        MODULE PROCEDURE comp_flux_stag
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

    ! SUBROUTINE comp_normal_strain_rate_pres(kk, jj, ii, splitDir, &
    ! u, v, w, ddx, ddy, ddz, strainRate)
    ! !----------------------------------------------------------------
    ! !   What it does:
    ! !   Computation of the normal strain rates dependend on the 
    ! !   directional split. 
    ! !----------------------------------------------------------------

    !     ! Subroutine arguments
    !     INTEGER(intk), INTENT(in) :: kk, jj, ii
    !     INTEGER(intk), INTENT(in) :: splitDir
    !     REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
    !     REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
    !     REAL(realk), INTENT(out) :: strainRate(kk, jj, ii)

    !     ! Local variables
    !     INTEGER(intk) :: k, j, i

    !     IF ( splitDir == 1 ) THEN
    !         DO i = 3, ii-2
    !             DO j = 3, jj-2
    !                 DO k = 3, kk-2
    !                     strainRate(k,j,i) = ( u(k,j,i) - u(k,j,i-1) ) / ddx(i)
    !                 END DO
    !             END DO
    !         END DO
    !     ELSE IF ( splitDir == 2 ) THEN
    !         DO i = 3, ii-2
    !             DO j = 3, jj-2
    !                 DO k = 3, kk-2
    !                     strainRate(k,j,i) = ( v(k,j,i) - v(k,j-1,i) ) / ddy(j)
    !                 END DO
    !             END DO
    !         END DO
    !     ELSE IF ( splitDir == 3 ) THEN
    !         DO i = 3, ii-2
    !             DO j = 3, jj-2
    !                 DO k = 3, kk-2
    !                     strainRate(k,j,i) = ( w(k,j,i) - w(k-1,j,i) ) / ddz(k)
    !                 END DO
    !             END DO
    !         END DO
    !     END IF

    ! END SUBROUTINE comp_normal_strain_rate_pres

    ! !================================================================

    ! SUBROUTINE comp_normal_strain_rate_stag(kk, jj, ii, q, splitDir, &
    !     u, v, w, dx, dy, dz, ddx, ddy, ddz, strainRate)
    ! !----------------------------------------------------------------
    ! !   What it does:
    ! !   Computation of the normal strain rates dependend on the 
    ! !   directional split for the three staggered grids.
    ! !----------------------------------------------------------------

    !     ! Subroutine arguments
    !     INTEGER(intk), INTENT(in) :: kk, jj, ii, q, splitDir
    !     REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
    !     REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
    !     REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
    !     REAL(realk), INTENT(out) :: strainRate(kk, jj, ii, 3)

    !     ! Local variables
    !     INTEGER(intk) :: k, j, i
    !     REAL(realk) :: advrE(kk, jj, ii), advrW(kk, jj, ii)
    !     REAL(realk) :: advrN(kk, jj, ii), advrS(kk, jj, ii)
    !     REAL(realk) :: advrT(kk, jj, ii), advrB(kk, jj, ii)
    !     REAL(realk) :: iStag, jStag, kStag
    !     REAL(realk) :: deltaX(ii), deltaY(jj), deltaZ(kk)

    !     CALL get_component_specifics(kk, jj, ii, q, dx=dx, dy=dy, dz=dz, &
    !         ddx=ddx, ddy=ddy, ddz=ddz, iStag=iStag, jStag=jStag, kStag=kStag, &
    !         deltaX=deltaX, deltaY=deltaY, deltaZ=deltaZ)
    !     CALL comp_advr_centr(kk, jj, ii, u, v, w, iStag, jStag, kStag, &
    !         advrE, advrW, advrN, advrS, advrT, advrB)

    !     IF ( splitDir == 1 ) THEN
    !         DO i = 3, ii-2
    !             DO j = 3, jj-2
    !                 DO k = 3, kk-2
    !                     strainRate(k,j,i,q) = ( advrE(k,j,i) - advrW(k,j,i) ) / deltaX(i)
    !                 END DO
    !             END DO
    !         END DO
    !     ELSE IF ( splitDir == 2 ) THEN
    !         DO i = 3, ii-2
    !             DO j = 3, jj-2
    !                 DO k = 3, kk-2
    !                     strainRate(k,j,i,q) = ( advrN(k,j,i) - advrS(k,j,i) ) / deltaY(j)
    !                 END DO
    !             END DO
    !         END DO
    !     ELSE IF ( splitDir == 3 ) THEN 
    !         DO i = 3, ii-2
    !             DO j = 3, jj-2
    !                 DO k = 3, kk-2
    !                     strainRate(k,j,i,q) = ( advrT(k,j,i) - advrB(k,j,i) ) / deltaZ(k)                 
    !                 END DO
    !             END DO
    !         END DO
    !     END IF

    ! END SUBROUTINE comp_normal_strain_rate_stag
    
    !================================================================

    SUBROUTINE comp_cWY(kk, jj, ii, vff, cWY)
    !----------------------------------------------------------------
    !   What it does:
    !   Computes the nondirectional compression coefficient c for 
    !   Weymouth and Yue's advection scheme.
    !
    !   Source:
    !   T. Arrufat et al., “A mass-momentum consistent, 
    !   Volume-of-Fluid method for incompressible flow on staggered 
    !   grids,” Computers & Fluids, vol. 215, p. 104785, Jan. 2021, 
    !   doi: 10.1016/j.compfluid.2020.104785.
    !   
    !   G. D. Weymouth and D. K.-P. Yue, “Conservative 
    !   Volume-of-Fluid method for free-surface simulations on 
    !   Cartesian-grids,” Journal of Computational Physics, vol. 229,
    !   no. 8, pp. 2853–2865, Apr. 2010, 
    !   doi: 10.1016/j.jcp.2009.12.018.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: vff(kk, jj, ii)
        REAL(realk), INTENT(out) :: cWY(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i

        DO i = 1, ii
            DO j = 1, jj
                DO k = 1, kk
                    IF ( vff(k,j,i) > 0.5_realk ) THEN
                        cWY(k,j,i) = 1.0_realk
                    ELSE
                        cWY(k,j,i) = 0.0_realk
                    END IF
                END DO
            END DO
        END DO

    END SUBROUTINE comp_cWY

    !================================================================

    SUBROUTINE comp_flux_pres(kk, jj, ii, splitDir, field, isInterface, u, v, w, alpha, dt, normx, normy, normz, ddx, ddy, ddz, tol, fieldFlux)
    !----------------------------------------------------------------
    !   What it does:
    !   Computes the volume fraction fluxes depending on the current
    !   split direction. It is distinguished between several cases 
    !   depending on the velocity direction and volume fraction 
    !   field.
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
        REAL(realk), INTENT(out) :: fieldFlux(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i
        REAL(realk) :: flux, fluxedProp, fluxWidth, fluxAlpha

        IF ( splitDir == 1 ) THEN 
            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        IF ( u(k,j,i) > tol ) THEN
                            IF ( isInterface(k,j,i) ) THEN
                                ! CASE 1
                                ! vff:      <1     E       
                                !       +----------+----------+
                                !       |     |    |          |
                                !       |     |   --> u       |
                                !       |     |    |          |
                                !       +----------+----------+
                                !             <---->
                                !             u * dt
                                ! 
                                ! The proportion of the cell, which 
                                ! is relevant for volume transport is
                                ! calculated. The alpha value is
                                ! transformed to match the new node
                                ! of the cell and the fluxed
                                ! proportion is computed.

                                fluxWidth = abs( u(k,j,i) ) * dt
                                fluxAlpha = alpha(k,j,i) - normx(k,j,i) * ( ddx(i) - fluxWidth )

                                CALL compute_cell_proportion(fluxedProp, fluxAlpha, field(k,j,i), &
                                    fluxWidth, ddy(j), ddz(k), normx(k,j,i), normy(k,j,i), normz(k,j,i), tol)
                                
                                flux = fluxedProp  * ( abs( u(k,j,i) ) * dt / ddx(i) )
                            ELSE
                                ! CASE 2
                                ! vff:      =1     E       
                                !       +----------+----------+
                                !       |          |          |
                                !       |         --> u       |
                                !       |          |          |
                                !       +----------+----------+
                                ! 
                                ! Since there is no interface in the 
                                ! upwind cell there is no need to 
                                ! compute the proportion.

                                flux = field(k,j,i) * ( abs( u(k,j,i) ) * dt / ddx(i) )
                            END IF
                        ELSE IF ( u(k,j,i) < -tol ) THEN
                            IF ( isInterface(k,j,i+1) ) THEN
                                ! CASE 3
                                ! vff:             E       
                                !       +----------+----------+
                                !       |          |    |     |
                                !       |       u <--   |     |
                                !       |          |    |     |
                                !       +----------+----------+
                                !                  <---->
                                !                  u * dt
                                ! 
                                ! The proportion of the cell, which 
                                ! is relevant for volume transport is
                                ! calculated. The alpha value is
                                ! transformed to match the new node
                                ! of the cell and the fluxed
                                ! proportion is computed.

                                fluxWidth = abs( u(k,j,i) ) * dt
                                fluxAlpha = alpha(k,j,i+1)

                                CALL compute_cell_proportion(fluxedProp, fluxAlpha, field(k,j,i+1), fluxWidth, ddy(j), ddz(k), normx(k,j,i+1), normy(k,j,i+1), normz(k,j,i+1), tol)

                                flux = fluxedProp * ( abs( u(k,j,i) ) * dt / ddx(i+1) )
                            ELSE
                                ! CASE 4
                                ! vff:             E    =1 
                                !       +----------+----------+
                                !       |          |          |
                                !       |       u <--         |
                                !       |          |          |
                                !       +----------+----------+
                                ! 
                                ! Since there is no interface in the  
                                ! upwind cell there is no need to 
                                ! compute the proportion.

                                flux = field(k,j,i+1) * ( abs( u(k,j,i) ) * dt / ddx(i+1) )
                            END IF
                        ELSE
                            ! CASE 5
                            ! Velocity is smaller than tolerance. 
                            ! Hence, there is no flux.

                            flux = 0.0_realk
                        END IF
                        fieldFlux(k,j,i) = sign( 1.0_realk, u(k,j,i) ) * flux / dt
                    END DO
                END DO
            END DO
        ELSE IF ( splitDir == 2 ) THEN 
            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        IF ( v(k,j,i) > tol ) THEN
                            IF ( isInterface(k,j,i) ) THEN
                                ! See explanation above.
                                fluxWidth = abs( v(k,j,i) ) * dt
                                fluxAlpha = alpha(k,j,i) - normy(k,j,i) * ( ddy(j) - fluxWidth )

                                CALL compute_cell_proportion(fluxedProp, fluxAlpha, field(k,j,i), ddx(i), fluxWidth, ddz(k), normx(k,j,i), normy(k,j,i), normz(k,j,i), tol)

                                flux = fluxedProp * ( abs( v(k,j,i) ) * dt / ddy(j) )
                            ELSE
                                ! See explanation above.
                                flux = field(k,j,i) * ( abs( v(k,j,i) ) * dt / ddy(j) )
                            END IF
                        ELSE IF ( v(k,j,i) < -tol ) THEN
                            IF ( isInterface(k,j+1,i) ) THEN
                                ! See explanation above.
                                fluxWidth = abs( v(k,j,i) ) * dt
                                fluxAlpha = alpha(k,j+1,i)

                                CALL compute_cell_proportion(fluxedProp, fluxAlpha, field(k,j+1,i), ddx(i), fluxWidth, ddz(k), normx(k,j+1,i), normy(k,j+1,i), normz(k,j+1,i), tol)

                                flux = fluxedProp * ( abs( v(k,j,i) ) * dt / ddy(j+1) )
                            ELSE
                                ! See explanation above.
                                flux = field(k,j+1,i) * ( abs( v(k,j,i) ) * dt / ddy(j+1) )
                            END IF
                        ELSE
                            ! See explanation above.
                            flux = 0.0_realk
                        END IF
                        fieldFlux(k,j,i) = sign( 1.0_realk, v(k,j,i) ) * flux / dt
                    END DO
                END DO
            END DO
        ELSE IF ( splitDir == 3 ) THEN
            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        IF ( w(k,j,i) > tol ) THEN
                            IF ( isInterface(k,j,i) ) THEN
                                ! See explanation above.
                                fluxWidth = abs( w(k,j,i) ) * dt
                                fluxAlpha = alpha(k,j,i) - normz(k,j,i) * ( ddz(k) - fluxWidth )

                                CALL compute_cell_proportion(fluxedProp, fluxAlpha, field(k,j,i), ddx(i), ddy(j), fluxWidth, normx(k,j,i), normy(k,j,i), normz(k,j,i), tol)

                                flux = fluxedProp * ( abs( w(k,j,i) ) * dt / ddz(k) )
                            ELSE
                                ! See explanation above.
                                flux = field(k,j,i) * ( abs( w(k,j,i) ) * dt / ddz(k) )
                            END IF
                        ELSE IF ( w(k,j,i) < -tol ) THEN
                            IF ( isInterface(k+1,j,i) ) THEN
                                ! See explanation above.
                                fluxWidth = abs( w(k,j,i) ) * dt
                                fluxAlpha = alpha(k+1,j,i)

                                CALL compute_cell_proportion(fluxedProp, fluxAlpha, field(k+1,j,i), ddx(i), ddy(j), fluxWidth, normx(k+1,j,i), normy(k+1,j,i), normz(k+1,j,i), tol)

                                flux = fluxedProp * ( abs( w(k,j,i) ) * dt / ddz(k+1) )
                            ELSE
                                ! See explanation above.
                                flux = field(k+1,j,i) * ( abs( w(k,j,i) ) * dt / ddz(k+1) )
                            END IF
                        ELSE
                            ! See explanation above.
                            flux = 0.0_realk
                        END IF
                        fieldFlux(k,j,i) = sign( 1.0_realk, w(k,j,i) ) * flux / dt
                    END DO
                END DO
            END DO
        END IF

    END SUBROUTINE comp_flux_pres

    !================================================================

    SUBROUTINE comp_flux_stag(kk, jj, ii, q, splitDir, field, isInterface, u, v, w, alpha, dt, normx, normy, normz, dx, dy, dz, ddx, ddy, ddz, tol, fieldFlux)
    !----------------------------------------------------------------
    !   What it does:
    !   Computes the volume fraction fluxes depending on the current
    !   split direction for the three staggered girds. It is 
    !   distinguished between several cases depending on the velocity 
    !   direction and volume fraction field.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, q
        INTEGER(intk), INTENT(in) :: splitDir
        REAL(realk), INTENT(in) :: field(kk, jj, ii)
        LOGICAL, INTENT(in) :: isInterface(kk, jj, ii)
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: alpha(kk, jj, ii)
        REAL(realk), INTENT(in) :: dt
        REAL(realk), INTENT(in) :: normx(kk, jj, ii), normy(kk, jj, ii), normz(kk, jj, ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: tol
        REAL(realk), INTENT(inout) :: fieldFlux(kk, jj, ii)

        ! Local variables
        REAL(realk) :: iStag, jStag, kStag
        REAL(realk) :: deltaX(ii), deltaY(jj), deltaZ(kk)
        INTEGER(intk) :: k, j, i
        REAL(realk) :: flux, complementFlux, fluxedProp, fluxWidth, fluxAlpha
        REAL(realk) :: advrE(kk, jj, ii), advrW(kk, jj, ii), advrN(kk, jj, ii), advrS(kk, jj, ii), advrT(kk, jj, ii), advrB(kk, jj, ii)
        REAL(realk) :: complementFieldFlux(kk, jj, ii)

        CALL get_component_specifics(kk, jj, ii, q, dx=dx, dy=dy, dz=dz, &
            ddx=ddx, ddy=ddy, ddz=ddz, iStag=iStag, jStag=jStag, kStag=kStag, &
            deltaX=deltaX, deltaY=deltaY, deltaZ=deltaZ)
        CALL comp_advr_centr(kk, jj, ii, u, v, w, iStag, jStag, kStag, &
            advrE, advrW, advrN, advrS, advrT, advrB)

        IF ( splitDir == 1 ) THEN
            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        IF ( advrE(k,j,i) > tol ) THEN
                            IF ( isInterface(k,j,i) ) THEN
                                ! See explanation above.
                                fluxWidth = abs( advrE(k,j,i) ) * dt
                                fluxAlpha = alpha(k,j,i) - normx(k,j,i) * ( deltaX(i)/2.0_realk - fluxWidth )

                                CALL compute_cell_proportion(fluxedProp, fluxAlpha, field(k,j,i), fluxWidth, ddy(j), ddz(k), normx(k,j,i), normy(k,j,i), normz(k,j,i), tol)
                                
                                flux = fluxedProp  * ( abs( advrE(k,j,i) ) * dt / deltaX(i) )
                                complementFlux = ( 1.0_realk - fluxedProp ) * ( abs( advrE(k,j,i) ) * dt / deltaX(i) )
                            ELSE
                                ! See explanation above.
                                flux = field(k,j,i) * ( abs( advrE(k,j,i) ) * dt / deltaX(i) )
                                complementFlux = max( 1.0_realk - field(k,j,i), 0.0_realk ) * ( abs( advrE(k,j,i) ) * dt / deltaX(i) )
                            END IF
                        ELSE IF ( advrE(k,j,i) < -tol ) THEN
                            IF ( isInterface(k,j,i+1) ) THEN
                                ! See explanation above.
                                fluxWidth = abs( advrE(k,j,i) ) * dt
                                fluxAlpha = alpha(k,j,i+1) - normx(k,j,i+1) * deltaX(i+1)/2.0_realk

                                CALL compute_cell_proportion(fluxedProp, fluxAlpha, field(k,j,i+1), fluxWidth, ddy(j), ddz(k), normx(k,j,i+1), normy(k,j,i+1), normz(k,j,i+1), tol)

                                flux = fluxedProp * ( abs( advrE(k,j,i) ) * dt / deltaX(i+1) )
                                complementFlux = ( 1.0_realk - fluxedProp ) * ( abs( advrE(k,j,i) ) * dt / deltaX(i+1) )
                            ELSE
                                ! See explanation above.
                                flux = field(k,j,i+1) * ( abs( advrE(k,j,i) ) * dt / deltaX(i+1) )
                                complementFlux = max( 1.0_realk - field(k,j,i+1), 0.0_realk ) * ( abs( advrE(k,j,i) ) * dt / deltaX(i+1) )
                            END IF
                        ELSE
                            ! See explanation above.
                            flux = 0.0_realk
                            complementFlux = 0.0_realk
                        END IF
                        fieldFlux(k,j,i) = sign( 1.0_realk, advrE(k,j,i) ) * flux / dt
                        complementFieldFlux(k,j,i) = sign( 1.0_realk, advrE(k,j,i) ) * complementFlux / dt
                    END DO
                END DO
            END DO
        ELSE IF ( splitDir == 2 ) THEN
            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        IF ( advrN(k,j,i) > tol ) THEN
                            IF ( isInterface(k,j,i) ) THEN
                                ! See explanation above.
                                fluxWidth = abs( advrN(k,j,i) ) * dt
                                fluxAlpha = alpha(k,j,i) - normy(k,j,i) * ( deltaY(j)/2.0_realk - fluxWidth )

                                CALL compute_cell_proportion(fluxedProp, fluxAlpha, field(k,j,i), ddx(i), fluxWidth, ddz(k), normx(k,j,i), normy(k,j,i), normz(k,j,i), tol)

                                flux = fluxedProp * ( abs( advrN(k,j,i) ) * dt / deltaY(j) )
                                complementFlux = ( 1.0_realk - fluxedProp ) * ( abs( advrN(k,j,i) ) * dt / deltaY(j) )
                            ELSE
                                ! See explanation above.
                                flux = field(k,j,i) * ( abs( advrN(k,j,i) ) * dt / deltaY(j) )
                                complementFlux = max( 1.0_realk - field(k,j,i), 0.0_realk ) * ( abs( advrN(k,j,i) ) * dt / deltaY(j) )
                            END IF
                        ELSE IF ( advrN(k,j,i) < -tol ) THEN
                            IF ( isInterface(k,j+1,i) ) THEN
                                ! See explanation above.
                                fluxWidth = abs( advrN(k,j,i) ) * dt
                                fluxAlpha = alpha(k,j+1,i) - normy(k,j+1,i) * deltaY(j+1)/2.0_realk

                                CALL compute_cell_proportion(fluxedProp, fluxAlpha, field(k,j+1,i), ddx(i), fluxWidth, ddz(k), normx(k,j+1,i), normy(k,j+1,i), normz(k,j+1,i), tol)

                                flux = fluxedProp * ( abs( advrN(k,j,i) ) * dt / deltaY(j+1) )
                                complementFlux = ( 1.0_realk - fluxedProp ) * ( abs( advrN(k,j,i) ) * dt / deltaY(j+1) )
                            ELSE
                                ! See explanation above.
                                flux = field(k,j+1,i) * ( abs( advrN(k,j,i) ) * dt / deltaY(j+1) )
                                complementFlux = max( 1.0_realk - field(k,j+1,i), 0.0_realk ) * ( abs( advrN(k,j,i) ) * dt / deltaY(j+1) )
                            END IF
                        ELSE
                            ! See explanation above.
                            flux = 0.0_realk
                            complementFlux = 0.0_realk
                        END IF
                        fieldFlux(k,j,i) = sign( 1.0_realk, advrN(k,j,i) ) * flux / dt
                        complementFieldFlux(k,j,i) = sign( 1.0_realk, advrN(k,j,i) ) * complementFlux / dt
                    END DO
                END DO
            END DO
        ELSE IF ( splitDir == 3 ) THEN
            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        IF ( advrT(k,j,i) > tol ) THEN
                            IF ( isInterface(k,j,i) ) THEN
                                ! See explanation above.
                                fluxWidth = abs( advrT(k,j,i) ) * dt
                                fluxAlpha = alpha(k,j,i) - normz(k,j,i) * ( deltaZ(k)/2.0_realk - fluxWidth )

                                CALL compute_cell_proportion(fluxedProp, fluxAlpha, field(k,j,i), ddx(i), ddy(j), fluxWidth, normx(k,j,i), normy(k,j,i), normz(k,j,i), tol)

                                flux = fluxedProp * ( abs( advrT(k,j,i) ) * dt / deltaZ(k) )
                                complementFlux = ( 1.0_realk - fluxedProp ) * ( abs( advrT(k,j,i) ) * dt / deltaZ(k) )
                            ELSE
                                ! See explanation above.
                                flux = field(k,j,i) * ( abs( advrT(k,j,i) ) * dt / deltaZ(k) )
                                complementFlux = max( 1.0_realk - field(k,j,i), 0.0_realk ) * ( abs( advrT(k,j,i) ) * dt / deltaZ(k) )
                            END IF
                        ELSE IF ( advrT(k,j,i) < -tol ) THEN
                            IF ( isInterface(k+1,j,i) ) THEN
                                ! See explanation above.
                                fluxWidth = abs( advrT(k,j,i) ) * dt
                                fluxAlpha = alpha(k+1,j,i) - normz(k+1,j,i) * deltaZ(k+1)/2.0_realk

                                CALL compute_cell_proportion(fluxedProp, fluxAlpha, field(k+1,j,i), ddx(i), ddy(j), fluxWidth, normx(k+1,j,i), normy(k+1,j,i), normz(k+1,j,i), tol)

                                flux = fluxedProp * ( abs( advrT(k,j,i) ) * dt / deltaZ(k+1) )
                                complementFlux = ( 1.0_realk - fluxedProp ) * ( abs( advrT(k,j,i) ) * dt / deltaZ(k+1) )
                            ELSE
                                ! See explanation above.
                                flux = field(k+1,j,i) * ( abs( advrT(k,j,i) ) * dt / deltaZ(k+1) )
                                complementFlux = max( 1.0_realk - field(k+1,j,i), 0.0_realk ) * ( abs( advrT(k,j,i) ) * dt / deltaZ(k+1) )
                            END IF
                        ELSE
                            ! See explanation above.
                            flux = 0.0_realk
                            complementFLux = 0.0_realk
                        END IF
                        fieldFlux(k,j,i) = sign( 1.0_realk, advrT(k,j,i) ) * flux / dt
                        complementFieldFlux(k,j,i) = sign( 1.0_realk, advrT(k,j,i) ) * complementFlux / dt
                    END DO
                END DO
            END DO
        END IF

    END SUBROUTINE comp_flux_stag

    !================================================================

    ! SUBROUTINE update_field_pres(kk, jj, ii, splitDir, flux, comprTerm, & 
    !         dt, nfro, nbac, nrgt, nlft, nbot, ntop, fld)
    ! !----------------------------------------------------------------
    ! !   What it does:
    ! !   Update the input field with corresponding fluxes and 
    ! !   compression term.
    ! !----------------------------------------------------------------

    !     ! Subroutine arguments
    !     INTEGER(intk), INTENT(in) :: kk, jj, ii
    !     INTEGER(intk), INTENT(in) :: splitDir
    !     REAL(realk), INTENT(in) :: flux(kk, jj, ii)
    !     REAL(realk), INTENT(in) :: comprTerm(kk, jj, ii)
    !     REAL(realk), INTENT(in) :: dt
    !     INTEGER, INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop
    !     REAL(realk), INTENT(inout) :: fld(kk, jj, ii)

    !     ! Local variables
    !     INTEGER(intk) :: k, j, i
    !     INTEGER(intk) :: nbu, nfu, nrv, nbw, ntw, nlv

    !     nfu = 0
    !     nbu = 0
    !     nrv = 0
    !     nlv = 0
    !     nbw = 0
    !     ntw = 0

    !     ! CON = 7
    !     IF (nbac == 7) nbu = 1
    !     IF (nlft == 7) nlv = 1
    !     IF (ntop == 7) ntw = 1

    !     ! OP1 = 3
    !     IF (nfro == 3) nfu = 1
    !     IF (nbac == 3) nbu = 1
    !     IF (nrgt == 3) nrv = 1
    !     IF (nlft == 3) nlv = 1
    !     IF (nbot == 3) nbw = 1
    !     IF (ntop == 3) ntw = 1
        
    !     IF ( splitDir == 1 ) THEN 
    !         DO i = 3-nfu, ii-3+nbu
    !             DO j = 3, jj-2
    !                 DO k = 3, kk-2
    !                     fld(k,j,i) = fld(k,j,i) + dt * ( flux(k,j,i-1) - flux(k,j,i) + comprTerm(k,j,i) )
    !                 END DO 
    !             END DO 
    !         END DO
    !     ELSEIF ( splitDir == 2 ) THEN
    !         DO i = 3, ii-2
    !             DO j = 3-nrv, jj-3+nlv
    !                 DO k = 3, kk-2
    !                     fld(k,j,i) = fld(k,j,i) + dt * ( flux(k,j-1,i) - flux(k,j,i) + comprTerm(k,j,i) )
    !                 END DO 
    !             END DO 
    !         END DO
    !     ELSEIF ( splitDir == 3 ) THEN
    !         DO i = 3, ii-2
    !             DO j = 3, jj-2
    !                 DO k = 3-nbw, kk-3+ntw
    !                     fld(k,j,i) = fld(k,j,i) + dt * ( flux(k-1,j,i) - flux(k,j,i) + comprTerm(k,j,i) )
    !                 END DO 
    !             END DO 
    !         END DO
    !     END IF

    ! END SUBROUTINE update_field_pres

    ! !================================================================

    ! SUBROUTINE update_field_stag(kk, jj, ii, q, splitDir, flux, comprTerm, & 
    !         dt, nfro, nbac, nrgt, nlft, nbot, ntop, fld)
    ! !----------------------------------------------------------------
    ! !   What it does:
    ! !   Update the input field with corresponding fluxes and 
    ! !   compression term for the three staggered girds.
    ! !----------------------------------------------------------------

    !     ! Subroutine arguments
    !     INTEGER(intk), INTENT(in) :: kk, jj, ii
    !     INTEGER(intk), INTENT(in) :: q
    !     INTEGER(intk), INTENT(in) :: splitDir
    !     REAL(realk), INTENT(in) :: flux(kk, jj, ii, 3)
    !     REAL(realk), INTENT(in) :: comprTerm(kk, jj, ii, 3)
    !     REAL(realk), INTENT(in) :: dt
    !     INTEGER, INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop
    !     REAL(realk), INTENT(inout) :: fld(kk, jj, ii, 3)

    !     ! Local variables
    !     INTEGER(intk) :: k, j, i
    !     INTEGER(intk) :: nbu, nfu, nrv, nbw, ntw, nlv

    !     nfu = 0
    !     nbu = 0
    !     nrv = 0
    !     nlv = 0
    !     nbw = 0
    !     ntw = 0

    !     ! CON = 7
    !     IF (nbac == 7) nbu = 1
    !     IF (nlft == 7) nlv = 1
    !     IF (ntop == 7) ntw = 1

    !     ! OP1 = 3
    !     IF (nfro == 3) nfu = 1
    !     IF (nbac == 3) nbu = 1
    !     IF (nrgt == 3) nrv = 1
    !     IF (nlft == 3) nlv = 1
    !     IF (nbot == 3) nbw = 1
    !     IF (ntop == 3) ntw = 1
        
    !     IF ( splitDir == 1 ) THEN 
    !         DO i = 3-nfu, ii-3+nbu
    !             DO j = 3, jj-2
    !                 DO k = 3, kk-2
    !                     fld(k,j,i,q) = fld(k,j,i,q) + dt * ( flux(k,j,i-1,q) - flux(k,j,i,q) + comprTerm(k,j,i,q) )
    !                 END DO 
    !             END DO 
    !         END DO
    !     ELSEIF ( splitDir == 2 ) THEN
    !         DO i = 3, ii-2
    !             DO j = 3-nrv, jj-3+nlv
    !                 DO k = 3, kk-2
    !                     fld(k,j,i,q) = fld(k,j,i,q) + dt * ( flux(k,j-1,i,q) - flux(k,j,i,q) + comprTerm(k,j,i,q) )
    !                 END DO 
    !             END DO 
    !         END DO
    !     ELSEIF ( splitDir == 3 ) THEN
    !         DO i = 3, ii-2
    !             DO j = 3, jj-2
    !                 DO k = 3-nbw, kk-3+ntw
    !                     fld(k,j,i,q) = fld(k,j,i,q) + dt * ( flux(k-1,j,i,q) - flux(k,j,i,q) + comprTerm(k,j,i,q) )
    !                 END DO 
    !             END DO 
    !         END DO
    !     END IF

    ! END SUBROUTINE update_field_stag

    ! !================================================================

    ! SUBROUTINE clip_volume_fraction_field(kk, ii, jj, vff, tol)
    ! !----------------------------------------------------------------
    ! !   What it does:
    ! !   This subroutine ensures that the volume fraction field vff 
    ! !   stays within its bounds of [0,1].
    ! !----------------------------------------------------------------

    !     ! Subroutine arguments
    !     INTEGER(intk), INTENT(in) :: kk, jj, ii
    !     REAL(realk), INTENT(inout) :: vff(kk, jj, ii)
    !     REAL(realk), INTENT(in) :: tol

    !     ! Local variables
    !     INTEGER(intk) :: k, j, i
        
    !     DO i = 3, ii-2
    !         DO j = 3, jj-2
    !             DO k = 3, kk-2
    !                 IF ( vff(k,j,i) < tol ) THEN
    !                     vff(k,j,i) = 0.0_realk
    !                 ELSE IF ( vff(k,j,i) > 1.0_realk - tol ) THEN
    !                     vff(k,j,i) = 1.0_realk
    !                 END IF
    !             END DO
    !         END DO
    !     END DO

    ! END SUBROUTINE clip_volume_fraction_field

    ! !================================================================

    ! SUBROUTINE compute_density_flux(kk, jj, ii, q, flux, fluxComp, rho1, rho2, densityFlux)
    ! !----------------------------------------------------------------
    ! !   What it does:
    ! !   The subroutine computes the density fluxes using the volume
    ! !   fraction field fluxes and their complements.
    ! !----------------------------------------------------------------

    !     ! Subroutine arguments
    !     INTEGER(intk), INTENT(in) :: kk, jj, ii
    !     INTEGER(intk), INTENT(in) :: q
    !     REAL(realk), INTENT(in) :: flux(kk, jj, ii, 3), fluxComp(kk, jj, ii, 3)
    !     REAL(realk), INTENT(in) :: rho1, rho2
    !     REAL(realk), INTENT(inout) :: densityFlux(kk, jj, ii, 3)

    !     ! Local variables
    !     INTEGER(intk) :: k, j, i
        
    !     DO i = 1, ii
    !         DO j = 1, jj
    !             DO k = 1, kk
    !                 densityFlux(k,j,i,q) = rho1 * flux(k,j,i,q) + rho2 * fluxComp(k,j,i,q)
    !             END DO
    !         END DO
    !     END DO

    ! END SUBROUTINE compute_density_flux

    ! !================================================================
    
    ! SUBROUTINE get_advection_sequence(iteration, advSeq)
    ! !----------------------------------------------------------------
    ! !   What it does:
    ! !   
    ! !----------------------------------------------------------------

    !     ! Subroutine arguments
    !     INTEGER(intk), INTENT(in) :: iteration
    !     INTEGER(intk), INTENT(inout) :: advSeq(3)

    !     ! Local variables
    !     INTEGER(intk) :: permutationIndex

    !     ! permutationIndex only changes in a new time-step
    !     permutationIndex = mod(iteration-1, 3)

    !     ! Select permutation of split advection
    !     SELECT CASE (permutationIndex)
    !         CASE (0)
    !             advSeq = [1, 2, 3]
    !         CASE (1)
    !             advSeq = [3, 1, 2]
    !         CASE (2)
    !             advSeq = [2, 3, 1]
    !     END SELECT

    ! END SUBROUTINE get_advection_sequence

    !================================================================

    SUBROUTINE multiphase_solve(u_f, v_f, w_f, f_f, p_f, g_f, d_f, dtrki, itstep)
    !----------------------------------------------------------------
    !   What it does:
    !    
    !----------------------------------------------------------------

        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: u_f
        TYPE(field_t), INTENT(inout) :: v_f
        TYPE(field_t), INTENT(inout) :: w_f
        TYPE(field_t), INTENT(in) :: f_f
        TYPE(field_t), INTENT(in) :: p_f
        TYPE(field_t), INTENT(in) :: g_f
        TYPE(field_t), INTENT(in) :: d_f
        REAL(realk), INTENT(in) :: dtrki
        INTEGER(intk), INTENT(in) :: itstep

        ! Local variables
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: u, v, w
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: vff, p, g, d

        TYPE(field_t), POINTER :: dx_f, dy_f, dz_f, ddx_f, ddy_f, ddz_f
        TYPE(field_t), POINTER :: rdx_f, rdy_f, rdz_f, rddx_f, rddy_f, rddz_f
        TYPE(field_t), POINTER :: vffiStag_f, vffjStag_f, vffkStag_f
        TYPE(field_t), POINTER :: diStag_f, djStag_f, dkStag_f
        TYPE(field_t), POINTER :: normx_f, normy_f, normz_f
        TYPE(field_t), POINTER :: normxiStag_f, normyiStag_f, normziStag_f
        TYPE(field_t), POINTER :: normxjStag_f, normyjStag_f, normzjStag_f
        TYPE(field_t), POINTER :: normxkStag_f, normykStag_f, normzkStag_f
        TYPE(field_t), POINTER :: alpha_f
        TYPE(field_t), POINTER :: alphaiStag_f, alphajStag_f, alphakStag_f

        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:), ddx(:), ddy(:), ddz(:)
        REAL(realk), POINTER, CONTIGUOUS :: rdx(:), rdy(:), rdz(:), rddx(:), rddy(:), rddz(:)
        REAL(realk), POINTER, CONTIGUOUS :: vffiStag(:,:,:), vffjStag(:,:,:), vffkStag(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: diStag(:,:,:), djStag(:,:,:), dkStag(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: normx(:,:,:), normy(:,:,:), normz(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: normxiStag(:,:,:), normyiStag(:,:,:), normziStag(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: normxjStag(:,:,:), normyjStag(:,:,:), normzjStag(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: normxkStag(:,:,:), normykStag(:,:,:), normzkStag(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: alpha(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: alphaiStag(:,:,:), alphajStag(:,:,:), alphakStag(:,:,:)

        INTEGER(intk) :: i, igrid
        INTEGER(intk) :: kk, jj, ii, q
        INTEGER(intk) :: nfro, nbac, nrgt, nlft, nbot, ntop
        REAL(realk), ALLOCATABLE :: vffStag(:,:,:)
        REAL(realk), ALLOCATABLE :: dStag(:,:,:)
        REAL(realk), ALLOCATABLE :: normxStag(:,:,:), normyStag(:,:,:), normzStag(:,:,:)
        REAL(realk), ALLOCATABLE :: alphaStag(:,:,:)
        LOGICAL, ALLOCATABLE :: isInterface(:,:,:)
        LOGICAL, ALLOCATABLE :: isInterfaceStag(:,:,:)
        LOGICAL, ALLOCATABLE :: isNearInterface(:,:,:)
        LOGICAL, ALLOCATABLE :: isNearInterfaceStag(:,:,:)
        REAL(realk), ALLOCATABLE :: vffFlux(:,:,:), mom(:,:,:), cWY(:,:,:), veloWY(:,:,:)
        
        REAL(realk), PARAMETER :: tol = 1.0E-15_realk

        CALL get_field(dx_f, "DX"); CALL get_field(dy_f, "DY"); CALL get_field(dz_f, "DZ")
        CALL get_field(ddx_f, "DDX"); CALL get_field(ddy_f, "DDY"); CALL get_field(ddz_f, "DDZ")

        CALL get_field(rdx_f, "RDX"); CALL get_field(rdy_f, "RDY"); CALL get_field(rdz_f, "RDZ")
        CALL get_field(rddx_f, "RDDX"); CALL get_field(rddy_f, "RDDY"); CALL get_field(rddz_f, "RDDZ")

        CALL get_field(vffiStag_f, "VFFiStag"); CALL get_field(vffjStag_f, "VFFjStag"); CALL get_field(vffkStag_f, "VFFkStag")
        CALL get_field(diStag_f, "DiStag"); CALL get_field(djStag_f, "DjStag"); CALL get_field(dkStag_f, "DkStag")
        CALL get_field(normx_f, "NORMX"); CALL get_field(normy_f, "NORMY"); CALL get_field(normz_f, "NORMZ")
        CALL get_field(normxiStag_f, "NORMXiStag"); CALL get_field(normyiStag_f, "NORMYiStag"); CALL get_field(normziStag_f, "NORMZiStag")
        CALL get_field(normxjStag_f, "NORMXjStag"); CALL get_field(normyjStag_f, "NORMYjStag"); CALL get_field(normzjStag_f, "NORMZjStag")
        CALL get_field(normxkStag_f, "NORMXkStag"); CALL get_field(normykStag_f, "NORMYkStag"); CALL get_field(normzkStag_f, "NORMZkStag")
        CALL get_field(alpha_f, "ALPHA")
        CALL get_field(alphaiStag_f, "ALPHAiStag"); CALL get_field(alphajStag_f, "ALPHAjStag"); CALL get_field(alphakStag_f, "ALPHAkStag")

        DO i = 1, nmygrids
            igrid = mygrids(i)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_mgbasb(nfro, nbac, nrgt, nlft, nbot, ntop, igrid)

            CALL u_f%get_ptr(u, igrid)
            CALL v_f%get_ptr(v, igrid)
            CALL w_f%get_ptr(w, igrid)
            CALL f_f%get_ptr(vff, igrid)
            CALL p_f%get_ptr(p, igrid)
            CALL g_f%get_ptr(g, igrid)
            CALL d_f%get_ptr(d, igrid)

            CALL dx_f%get_ptr(dx, igrid); CALL dy_f%get_ptr(dy, igrid); CALL dz_f%get_ptr(dz, igrid)
            CALL ddx_f%get_ptr(ddx, igrid); CALL ddy_f%get_ptr(ddy, igrid); CALL ddz_f%get_ptr(ddz, igrid)

            CALL rdx_f%get_ptr(rdx, igrid); CALL rdy_f%get_ptr(rdy, igrid); CALL rdz_f%get_ptr(rdz, igrid)
            CALL rddx_f%get_ptr(rddx, igrid); CALL rddy_f%get_ptr(rddy, igrid); CALL rddz_f%get_ptr(rddz, igrid)

            CALL vffiStag_f%get_ptr(vffiStag, igrid); CALL vffjStag_f%get_ptr(vffjStag, igrid); CALL vffkStag_f%get_ptr(vffkStag, igrid)
            CALL diStag_f%get_ptr(diStag, igrid); CALL djStag_f%get_ptr(djStag, igrid); CALL dkStag_f%get_ptr(dkStag, igrid)
            CALL normx_f%get_ptr(normx, igrid); CALL normy_f%get_ptr(normy, igrid); CALL normz_f%get_ptr(normz, igrid)
            CALL normxiStag_f%get_ptr(normxiStag, igrid); CALL normyiStag_f%get_ptr(normyiStag, igrid); CALL normziStag_f%get_ptr(normziStag, igrid)
            CALL normxjStag_f%get_ptr(normxjStag, igrid); CALL normyjStag_f%get_ptr(normyjStag, igrid); CALL normzjStag_f%get_ptr(normzjStag, igrid)
            CALL normxkStag_f%get_ptr(normxkStag, igrid); CALL normykStag_f%get_ptr(normykStag, igrid); CALL normzkStag_f%get_ptr(normzkStag, igrid)
            CALL alpha_f%get_ptr(alpha, igrid)
            CALL alphaiStag_f%get_ptr(alphaiStag, igrid); CALL alphajStag_f%get_ptr(alphajStag, igrid); CALL alphakStag_f%get_ptr(alphakStag, igrid)

            IF (.NOT. ALLOCATED(vffStag))             ALLOCATE(vffStag(kk,jj,ii))
            IF (.NOT. ALLOCATED(dStag))               ALLOCATE(dStag(kk,jj,ii))
            IF (.NOT. ALLOCATED(normxStag))           ALLOCATE(normxStag(kk,jj,ii))
            IF (.NOT. ALLOCATED(normyStag))           ALLOCATE(normyStag(kk,jj,ii))
            IF (.NOT. ALLOCATED(normzStag))           ALLOCATE(normzStag(kk,jj,ii))
            IF (.NOT. ALLOCATED(alphaStag))           ALLOCATE(alphaStag(kk,jj,ii))
            IF (.NOT. ALLOCATED(isInterface))         ALLOCATE(isInterface(kk,jj,ii))
            IF (.NOT. ALLOCATED(isInterfaceStag))     ALLOCATE(isInterfaceStag(kk,jj,ii))
            IF (.NOT. ALLOCATED(isNearInterface))     ALLOCATE(isNearInterface(kk,jj,ii))
            IF (.NOT. ALLOCATED(isNearInterfaceStag)) ALLOCATE(isNearInterfaceStag(kk,jj,ii))
            IF (.NOT. ALLOCATED(vffFlux))             ALLOCATE(vffFlux(kk,jj,ii))
            IF (.NOT. ALLOCATED(mom))                 ALLOCATE(mom(kk,jj,ii))
            IF (.NOT. ALLOCATED(cWY))                 ALLOCATE(cWY(kk,jj,ii))
            IF (.NOT. ALLOCATED(veloWY))              ALLOCATE(veloWY(kk,jj,ii))
            
            IF ( mom_multiphase == "RUDMAN" ) THEN

                CALL interface_reconstruction_wrapper(kk, jj, ii, vff, ddx, ddy, ddz, tol, normx, normy, normz, alpha, isInterface, isNearInterface)
                
                DO q = 1, 3

                    vffStag = 0.0_realk; normxStag = 0.0_realk; normyStag = 0.0_realk; normzStag = 0.0_realk
                    alphaStag = 0.0_realk; isInterfaceStag = .False.; isNearInterfaceStag = .False.
                    dStag = 0.0_realk; mom = 0.0_realk; cWY = 0.0_realk; vffFlux = 0.0_realk
                    
                    CALL staggered_fractions_wrapper(kk, jj, ii, q, alpha, vff, isInterface, ddx, ddy, ddz, normx, normy, normz, tol, vffStag)
                    CALL interface_reconstruction_wrapper(kk, jj, ii, vffStag, ddx, ddy, ddz, tol, normxStag, normyStag, normzStag, alphaStag, isInterfaceStag, isNearInterfaceStag)
                    CALL comp_material_property_field(kk, jj, ii, vffStag, rho1, rho2, dStag)
                    CALL comp_momentum(kk, jj, ii, q, dStag, u, v, w, mom)
                    CALL comp_cWY(kk, jj, ii, vffStag, cWY)

                    IF ( MOD(itstep, 3) == 0 ) THEN ! z, x, y
                        
                        CALL comp_flux(kk, jj, ii, q, 3, vff, isInterfaceStag, u, v, w, alphaStag, dtrki, normxStag, normyStag, normzStag, dx, dy, dz, ddx, ddy, ddz, tol, vffFlux)
                        CALL adv_vof(kk, jj, ii, q, 3, vffFlux, cWY, veloWY, dtrki, vffStag)
                        
                        CALL comp_flux(kk, jj, ii, q, 1, vff, isInterfaceStag, u, v, w, alphaStag, dtrki, normxStag, normyStag, normzStag, dx, dy, dz, ddx, ddy, ddz, tol, vffFlux)
                        CALL adv_vof(kk, jj, ii, q, 1, vffFlux, cWY, veloWY, dtrki, vffStag)

                        CALL comp_flux(kk, jj, ii, q, 2, vff, isInterfaceStag, u, v, w, alphaStag, dtrki, normxStag, normyStag, normzStag, dx, dy, dz, ddx, ddy, ddz, tol, vffFlux)
                        CALL adv_vof(kk, jj, ii, q, 2, vffFlux, cWY, veloWY, dtrki, vffStag)

                    ELSE IF ( MOD(itstep, 3) == 1 ) THEN ! y, z, x
                        
                        CALL comp_flux(kk, jj, ii, q, 2, vff, isInterfaceStag, u, v, w, alphaStag, dtrki, normxStag, normyStag, normzStag, dx, dy, dz, ddx, ddy, ddz, tol, vffFlux)
                        CALL adv_vof(kk, jj, ii, q, 2, vffFlux, cWY, veloWY, dtrki, vffStag)

                        CALL comp_flux(kk, jj, ii, q, 3, vff, isInterfaceStag, u, v, w, alphaStag, dtrki, normxStag, normyStag, normzStag, dx, dy, dz, ddx, ddy, ddz, tol, vffFlux)
                        CALL adv_vof(kk, jj, ii, q, 3, vffFlux, cWY, veloWY, dtrki, vffStag)
                        
                        CALL comp_flux(kk, jj, ii, q, 1, vff, isInterfaceStag, u, v, w, alphaStag, dtrki, normxStag, normyStag, normzStag, dx, dy, dz, ddx, ddy, ddz, tol, vffFlux)
                        CALL adv_vof(kk, jj, ii, q, 1, vffFlux, cWY, veloWY, dtrki, vffStag)

                    ELSE IF ( MOD(itstep, 3) == 2 ) THEN ! x, y, z
                                                
                        CALL comp_flux(kk, jj, ii, q, 1, vff, isInterfaceStag, u, v, w, alphaStag, dtrki, normxStag, normyStag, normzStag, dx, dy, dz, ddx, ddy, ddz, tol, vffFlux)
                        CALL adv_vof(kk, jj, ii, q, 1, vffFlux, cWY, veloWY, dtrki, vffStag)

                        CALL comp_flux(kk, jj, ii, q, 2, vff, isInterfaceStag, u, v, w, alphaStag, dtrki, normxStag, normyStag, normzStag, dx, dy, dz, ddx, ddy, ddz, tol, vffFlux)
                        CALL adv_vof(kk, jj, ii, q, 2, vffFlux, cWY, veloWY, dtrki, vffStag)

                        CALL comp_flux(kk, jj, ii, q, 3, vff, isInterfaceStag, u, v, w, alphaStag, dtrki, normxStag, normyStag, normzStag, dx, dy, dz, ddx, ddy, ddz, tol, vffFlux)
                        CALL adv_vof(kk, jj, ii, q, 3, vffFlux, cWY, veloWY, dtrki, vffStag)

                    END IF

                    IF ( q == 1 ) THEN
                        vffiStag = vffStag
                        diStag = dStag
                        normxiStag = normxStag
                        normyiStag = normyStag
                        normziStag = normzStag
                        alphaiStag = alphaStag
                    ELSE IF ( q == 2 ) THEN
                        vffjStag = vffStag
                        djStag = dStag
                        normxjStag = normxStag
                        normyjStag = normyStag
                        normzjStag = normzStag
                        alphajStag = alphaStag
                    ELSE IF ( q == 3 ) THEN
                        vffkStag = vffStag
                        dkStag = dStag
                        normxkStag = normxStag
                        normykStag = normyStag
                        normzkStag = normzStag
                        alphakStag = alphaStag
                    END IF

                END DO

                CALL comp_cWY(kk, jj, ii, vff, cWY)

                IF ( MOD(itstep, 3) == 0 ) THEN ! z, x, y
                    
                    CALL comp_flux(kk, jj, ii, 3, vff, isInterface, u, v, w, alpha, dtrki, normx, normy, normz, ddx, ddy, ddz, tol, vffFlux)
                    CALL adv_vof(kk, jj, ii, q, 3, vffFlux, cWY, veloWY, dtrki, vff)
                    
                    CALL comp_flux(kk, jj, ii, 1, vff, isInterface, u, v, w, alpha, dtrki, normx, normy, normz, ddx, ddy, ddz, tol, vffFlux)
                    CALL adv_vof(kk, jj, ii, q, 1, vffFlux, cWY, veloWY, dtrki, vff)

                    CALL comp_flux(kk, jj, ii, 2, vff, isInterface, u, v, w, alpha, dtrki, normx, normy, normz, ddx, ddy, ddz, tol, vffFlux)
                    CALL adv_vof(kk, jj, ii, q, 2, vffFlux, cWY, veloWY, dtrki, vff)

                ELSE IF ( MOD(itstep, 3) == 1 ) THEN ! y, z, x
                    
                    CALL comp_flux(kk, jj, ii, 2, vff, isInterface, u, v, w, alpha, dtrki, normx, normy, normz, ddx, ddy, ddz, tol, vffFlux)
                    CALL adv_vof(kk, jj, ii, q, 2, vffFlux, cWY, veloWY, dtrki, vff)

                    CALL comp_flux(kk, jj, ii, 3, vff, isInterface, u, v, w, alpha, dtrki, normx, normy, normz, ddx, ddy, ddz, tol, vffFlux)
                    CALL adv_vof(kk, jj, ii, q, 3, vffFlux, cWY, veloWY, dtrki, vff)
                    
                    CALL comp_flux(kk, jj, ii, 1, vff, isInterface, u, v, w, alpha, dtrki, normx, normy, normz, ddx, ddy, ddz, tol, vffFlux)
                    CALL adv_vof(kk, jj, ii, q, 1, vffFlux, cWY, veloWY, dtrki, vff)

                ELSE IF ( MOD(itstep, 3) == 2 ) THEN ! x, y, z
                                            
                    CALL comp_flux(kk, jj, ii, 1, vff, isInterface, u, v, w, alpha, dtrki, normx, normy, normz, ddx, ddy, ddz, tol, vffFlux)
                    CALL adv_vof(kk, jj, ii, q, 1, vffFlux, cWY, veloWY, dtrki, vff)

                    CALL comp_flux(kk, jj, ii, 2, vff, isInterface, u, v, w, alpha, dtrki, normx, normy, normz, ddx, ddy, ddz, tol, vffFlux)
                    CALL adv_vof(kk, jj, ii, q, 2, vffFlux, cWY, veloWY, dtrki, vff)

                    CALL comp_flux(kk, jj, ii, 3, vff, isInterface, u, v, w, alpha, dtrki, normx, normy, normz, ddx, ddy, ddz, tol, vffFlux)
                    CALL adv_vof(kk, jj, ii, q, 3, vffFlux, cWY, veloWY, dtrki, vff)

                END IF

            ELSE IF ( mom_multiphase == "ARRUFAT" ) THEN
                WRITE(*,*) "Not yet implemented!"
                continue
            END IF

        END DO

    END SUBROUTINE multiphase_solve

    !================================================================

    SUBROUTINE adv_vof(kk, jj, ii, q, splitDir, vffFlux, cWY, veloWY, dtrki, vff)
    !----------------------------------------------------------------
    !   What it does:
    !    
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: q, splitDir
        REAL(realk), INTENT(in) :: vffFlux(:,:,:)
        REAL(realk), INTENT(in) :: cWY(:,:,:), veloWY(:,:,:)
        REAL(realk), INTENT(in) :: dtrki
        REAL(realk), INTENT(inout) :: vff(:,:,:)

        ! Local variables
        INTEGER(intk) :: i, j, k
        INTEGER(intk) :: i0, j0, k0

        CALL get_component_specifics(kk, jj, ii, q=q, i0=i0, j0=j0, k0=k0)
        WRITE(*,*) i0, j0, k0

        DO i = 2, ii-1
            DO j = 2, jj-1
                DO k = 2, kk-1
                    vff(k,j,i) = vff(k,j,i) - dtrki * ( vffFlux(k,j,i) - vffFlux(k-i0,j-j0,i-k0) ) + dtrki * cWY(k,j,i)
                END DO
            END DO
        END DO

    END SUBROUTINE adv_vof

    SUBROUTINE comp_momentum(kk, jj, ii, q, dStag, u, v, w, mom)
    !----------------------------------------------------------------
    !   What it does:
    !    
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: q
        REAL(realk), INTENT(in) :: dStag(kk, jj, ii)
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(out) :: mom(kk,jj,ii)

        ! Local variables
        INTEGER(intk) :: i, j, k
        REAL(realk) :: veloFld(kk, jj, ii)

        CALL get_component_specifics(kk, jj, ii, q=q, u=u, v=v, w=w, veloFld=veloFld)

        DO i = 1, ii
            DO j = 1, jj
                DO k = 1, kk
                    mom(k,j,i) = veloFld(k,j,i) * dStag(k,j,i)
                END DO
            END DO
        END DO

    END SUBROUTINE comp_momentum

    !================================================================

    ! SUBROUTINE comp_velocity(kk, jj, ii, q, mom, densityFieldStag, tol, velocity)
    ! !----------------------------------------------------------------
    ! !   What it does:
    ! !    
    ! !----------------------------------------------------------------

    !     ! Subroutine arguments
    !     INTEGER(intk), INTENT(in) :: kk, jj, ii
    !     INTEGER(intk), INTENT(in) :: q
    !     REAL(realk), INTENT(in) :: mom(kk, jj, ii, 3), densityFieldStag(kk, jj, ii, 3)
    !     REAL(realk), INTENT(in) :: tol
    !     REAL(realk), INTENT(out) :: velocity(kk,jj,ii)

    !     ! Local variables
    !     INTEGER(intk) :: i, j, k
    
    !     DO i = 1, ii
    !         DO j = 1, jj
    !             DO k = 1, kk
    !                 IF (densityFieldStag(k,j,i,q) > tol) THEN
    !                     velocity(k,j,i) = mom(k,j,i,q) / densityFieldStag(k,j,i,q)
    !                 ELSE
    !                     velocity(k,j,i) = 0.0_realk
    !                 END IF
    !             END DO
    !         END DO
    !     END DO

    ! END SUBROUTINE

    ! !================================================================

    ! SUBROUTINE multiphase_momentum_advection(kk, jj, ii, q, splitDir, mom, densityFieldStag, densityFieldFluxStag, &
    !                     densityCompressionTermStag, dt, dx, dy, dz, ddx, ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz, &
    !                     nfro, nbac, nrgt, nlft, nbot, ntop)
    ! !----------------------------------------------------------------
    ! !   What it does:
    ! !    
    ! !----------------------------------------------------------------

    !     ! Subroutine arguments
    !     INTEGER(intk), INTENT(in) :: kk, jj, ii
    !     INTEGER(intk), INTENT(in) :: splitDir
    !     INTEGER(intk), INTENT(in) :: q
    !     REAL(realk), INTENT(inout) :: mom(kk, jj, ii, 3)
    !     REAL(realk), INTENT(in) :: densityFieldStag(kk, jj, ii, 3)
    !     REAL(realk), INTENT(in) :: densityFieldFluxStag(kk, jj, ii, 3) 
    !     REAL(realk), INTENT(in) :: densityCompressionTermStag(kk, jj, ii, 3)
    !     REAL(realk), INTENT(in) :: dt
    !     REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
    !     REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
    !     REAL(realk), INTENT(in) :: rdx(ii), rdy(jj), rdz(kk)
    !     REAL(realk), INTENT(in) :: rddx(ii), rddy(jj), rddz(kk)
    !     INTEGER, INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop

    !     ! Local variables
    !     INTEGER(intk) :: k, j, i
    !     INTEGER(intk) :: nbu, nfu, nrv, nlv, nbw, ntw
    !     REAL(realk) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
    !     REAL(realk) :: advrE(kk, jj, ii), advrW(kk, jj, ii), advrN(kk, jj, ii), advrS(kk, jj, ii), advrT(kk, jj, ii), advrB(kk, jj, ii)
    !     REAL(realk) :: adveE(kk, jj, ii), adveW(kk, jj, ii), adveN(kk, jj, ii), adveS(kk, jj, ii), adveT(kk, jj, ii), adveB(kk, jj, ii)
    !     REAL(realk) :: iStag, jStag, kStag
    !     REAL(realk) :: velocity(kk, jj, ii)

    !     nfu = 0
    !     nbu = 0
    !     nrv = 0
    !     nlv = 0
    !     nbw = 0
    !     ntw = 0

    !     ! CON = 7
    !     IF (nbac == 7) nbu = 1
    !     IF (nlft == 7) nlv = 1
    !     IF (ntop == 7) ntw = 1

    !     ! OP1 = 3
    !     IF (nfro == 3) nfu = 1
    !     IF (nbac == 3) nbu = 1
    !     IF (nrgt == 3) nrv = 1
    !     IF (nlft == 3) nlv = 1
    !     IF (nbot == 3) nbw = 1
    !     IF (ntop == 3) ntw = 1

    !     IF ( q == 1 ) THEN
    !         iStag = 1.0_realk
    !         jStag = 0.0_realk
    !         kStag = 0.0_realk
    !         velocity = mom(:,:,:,q) / densityFieldStag(:,:,:,q)
    !     ELSE IF ( q == 2 ) THEN
    !         iStag = 0.0_realk
    !         jStag = 1.0_realk
    !         kStag = 0.0_realk
    !         velocity = mom(:,:,:,q) / densityFieldStag(:,:,:,q)
    !     ELSE IF ( q == 3 ) THEN
    !         iStag = 0.0_realk
    !         jStag = 0.0_realk
    !         kStag = 1.0_realk
    !         velocity = mom(:,:,:,q) / densityFieldStag(:,:,:,q)
    !     END IF

    !     u = mom(:,:,:,1) / densityFieldStag(:,:,:,1)
    !     v = mom(:,:,:,2) / densityFieldStag(:,:,:,2)
    !     w = mom(:,:,:,3) / densityFieldStag(:,:,:,3)
        
    !     CALL comp_advr_centr(kk, jj, ii, u, v, w, iStag, jStag, kStag, advrE, advrW, advrN, advrS, advrT, advrB)
    !     CALL comp_adve_quick(kk, jj, ii, velocity, advrE, advrW, advrN, advrS, advrT, advrB, adveE, adveW, adveN, adveS, adveT, adveB)

    !     IF ( splitDir == 1 ) THEN
    !         DO i = 4-nfu, ii-4+nbu
    !             DO j = 4, jj-3
    !                 DO k = 4, kk-4
    !                     mom(k,j,i,q) = mom(k,j,i,q) - dt * ( adveE(k,j,i) * densityFieldFluxStag(k,j,i,q) - adveW(k,j,i) * densityFieldFluxStag(k,j,i-1,q) - velocity(k,j,i) * densityCompressionTermStag(k,j,i,q) )
    !                 END DO
    !             END DO
    !         END DO
    !     ELSE IF ( splitDir == 2 ) THEN
    !         DO i = 4, ii-3
    !             DO j = 4-nrv, jj-4+nlv
    !                 DO k = 4, kk-3
    !                     mom(k,j,i,q) = mom(k,j,i,q) - dt * ( adveN(k,j,i) * densityFieldFluxStag(k,j,i,q) - adveS(k,j,i) * densityFieldFluxStag(k,j-1,i,q) - velocity(k,j,i) * densityCompressionTermStag(k,j,i,q) )
    !                 END DO
    !             END DO
    !         END DO
    !     ELSE IF ( splitDir == 3 ) THEN
    !         DO i = 4, ii-3
    !             DO j = 4, jj-3
    !                 DO k = 4-nbw, kk-4+ntw
    !                     mom(k,j,i,q) = mom(k,j,i,q) - dt * ( adveT(k,j,i) * densityFieldFluxStag(k,j,i,q) - adveB(k,j,i) * densityFieldFluxStag(k-1,j,i,q) - velocity(k,j,i) * densityCompressionTermStag(k,j,i,q) )
    !                 END DO
    !             END DO
    !         END DO
    !     END IF

    ! END SUBROUTINE multiphase_momentum_advection

    !================================================================

    ! SUBROUTINE multiphase_momentum_diffusion(kk, jj, ii, uo, vo, wo, u, v, w, g, &
    !         dx, dy, dz, ddx, ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz, &
    !         nfro, nbac, nrgt, nlft, nbot, ntop, densityFieldiStag, densityFieldjStag, densityFieldkStag)
    ! !----------------------------------------------------------------
    ! !   What it does:
    ! !    
    ! !----------------------------------------------------------------
    
    !     ! Subroutine arguments
    !     INTEGER(intk), INTENT(in) :: kk, jj, ii
    !     REAL(realk), INTENT(inout) :: uo(kk, jj, ii), vo(kk, jj, ii), &
    !         wo(kk, jj, ii)
    !     REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
    !     REAL(realk), INTENT(in) :: g(kk, jj, ii)
    !     REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
    !     REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
    !     REAL(realk), INTENT(in) :: rdx(ii), rdy(jj), rdz(kk)
    !     REAL(realk), INTENT(in) :: rddx(ii), rddy(jj), rddz(kk)
    !     INTEGER, INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop
    !     REAL(realk), INTENT(in) :: densityFieldiStag(kk, jj, ii), densityFieldjStag(kk, jj, ii), densityFieldkStag(kk, jj, ii)

    !     ! Local variables
    !     INTEGER(intk) :: k, j, i
    !     INTEGER(intk) :: nbu, nfu, nrv, nbw, ntw, nlv
    !     REAL(realk) :: ge, gw, gn, gs, gt, gb
    !     REAL(realk) :: tauxxe, tauxxw, tauyxn, tauyxs, tauzxt, tauzxb
    !     REAL(realk) :: tauxye, tauxyw, tauyyn, tauyys, tauzyt, tauzyb
    !     REAL(realk) :: tauxze, tauxzw, tauyzn, tauyzs, tauzzt, tauzzb
    !     REAL(realk) :: duo, dvo, dwo

    !     nfu = 0
    !     nbu = 0
    !     nrv = 0
    !     nlv = 0
    !     nbw = 0
    !     ntw = 0

    !     ! CON = 7
    !     IF (nbac == 7) nbu = 1
    !     IF (nlft == 7) nlv = 1
    !     IF (ntop == 7) ntw = 1

    !     ! OP1 = 3
    !     IF (nfro == 3) nfu = 1
    !     IF (nbac == 3) nbu = 1
    !     IF (nrgt == 3) nrv = 1
    !     IF (nlft == 3) nlv = 1
    !     IF (nbot == 3) nbw = 1
    !     IF (ntop == 3) ntw = 1

    !     ! CALL swcle3d(kk, jj, ii, uo, vo, wo, u, v, w, &
    !     !     ddx, ddy, ddz, nfro, nbac, nrgt, nlft, nbot, ntop)

    !     DO i = 3-nfu, ii-3+nbu
    !         DO j = 3, jj-2
    !             DO k = 3, kk-2
    !                 ! Face values of dynamic viscosity on u-momentum cell
    !                 ! Harmonic mean for a more physical treatment at interfaces
    !                 ge = g(k, j, i+1)
    !                 gw = g(k, j, i)
    !                 gn = g(k, j, i)*g(k, j+1, i) &
    !                     / MAX(g(k, j, i) + g(k, j+1, i), MIN(gmol1,gmol2)) &
    !                     + g(k, j, i+1)*g(k, j+1, i+1) &
    !                     / MAX(g(k, j, i+1) + g(k, j+1, i+1), MIN(gmol1,gmol2))
    !                 gs = g(k, j-1, i)*g(k, j, i) &
    !                     / MAX(g(k, j-1, i) + g(k, j, i), MIN(gmol1,gmol2)) &
    !                     + g(k, j-1, i+1)*g(k, j, i+1) &
    !                     / MAX(g(k, j-1, i+1) + g(k, j, i+1), MIN(gmol1,gmol2))
    !                 gt = g(k, j, i)*g(k+1, j, i) &
    !                     / MAX(g(k, j, i) + g(k+1, j, i), MIN(gmol1,gmol2)) &
    !                     + g(k, j, i+1)*g(k+1, j, i+1) &
    !                     /MAX(g(k, j, i+1) + g(k+1, j, i+1), MIN(gmol1,gmol2))
    !                 gb = g(k-1, j, i)*g(k, j, i) &
    !                     / MAX(g(k-1, j, i) + g(k, j, i), MIN(gmol1,gmol2)) &
    !                     + g(k-1, j, i+1)*g(k, j, i+1) &
    !                     / MAX(g(k-1, j, i+1) + g(k, j, i+1), MIN(gmol1,gmol2))

    !                 ! Normal stresses
    !                 !             ---------------inner derivatives---------------
    !                 tauxxe = ge * 2.0_realk * (u(k,j,i+1) - u(k,j,i)) * rddx(i+1)
    !                 tauxxw = gw * 2.0_realk * (u(k,j,i) - u(k,j,i-1)) * rddx(i)

    !                 ! Shear stresses
    !                 !             ------------------------------inner derivatives------------------------------
    !                 tauyxn = gn * ( (u(k,j+1,i) - u(k,j,i)) * rdy(j)   + (v(k,j,i+1) - v(k,j,i))     * rdx(i) )
    !                 tauyxs = gs * ( (u(k,j,i) - u(k,j-1,i)) * rdy(j-1) + (v(k,j-1,i+1) - v(k,j-1,i)) * rdx(i) )
    !                 tauzxt = gt * ( (u(k+1,j,i) - u(k,j,i)) * rdz(k)   + (w(k,j,i+1) - w(k,j,i))     * rdx(i) )
    !                 tauzxb = gb * ( (u(k,j,i) - u(k-1,j,i)) * rdz(k-1) + (w(k-1,j,i+1) - w(k-1,j,i)) * rdx(i) )

    !                 ! Change due to diffusion
    !                 !                                  ---------------------------------------outer derivatives----------------------------------------
    !                 duo = 1/densityFieldiStag(k,j,i) * ( ( tauxxe - tauxxw ) * rdx(i) + ( tauyxn - tauyxs ) * rddy(j) + ( tauzxt - tauzxb ) * rddz(k) )

    !                 ! Addition
    !                 uo(k, j, i) = uo(k, j, i) + duo
    !             END DO
    !         END DO
    !     END DO

    !     DO i = 3, ii-2
    !         DO j = 3-nrv, jj-3+nlv
    !             DO k = 3, kk-2
    !                 ! Face values of dynamic viscosity on v-momentum cell
    !                 ! Harmonic mean for a more physical treatment at interfaces
    !                 ge = g(k, j, i)*g(k, j, i+1) &
    !                     / MAX(g(k, j, i) + g(k, j, i+1), MIN(gmol1,gmol2)) &
    !                     + g(k, j+1, i)*g(k, j+1, i+1) &
    !                     / MAX(g(k, j+1, i) + g(k, j+1, i+1), MIN(gmol1,gmol2))
    !                 gw = g(k, j, i-1)*g(k, j, i) &
    !                     / MAX(g(k, j, i-1) + g(k, j, i), MIN(gmol1,gmol2)) &
    !                     + g(k, j+1, i-1)*g(k, j+1, i) &
    !                     / MAX(g(k, j+1, i-1) + g(k, j+1, i), MIN(gmol1,gmol2))
    !                 gn = g(k, j+1, i)
    !                 gs = g(k, j, i)
    !                 gt = g(k, j, i)*g(k+1, j, i) &
    !                     / MAX(g(k, j, i) + g(k+1, j, i), MIN(gmol1,gmol2)) &
    !                     + g(k, j+1, i)*g(k+1, j+1, i) &
    !                     / MAX(g(k, j+1, i) + g(k+1, j+1, i), MIN(gmol1,gmol2))
    !                 gb = g(k-1, j, i)*g(k, j, i) &
    !                     / MAX(g(k-1, j, i) + g(k, j, i), MIN(gmol1,gmol2)) &
    !                     + g(k-1, j+1, i)*g(k, j+1, i) &
    !                     / MAX(g(k-1, j+1, i) + g(k, j+1, i), MIN(gmol1,gmol2))

    !                 ! Shear stresses
    !                 !             ------------------------------inner derivatives------------------------------
    !                 tauxye = ge * ( (u(k,j+1,i) - u(k,j,i))     * rdy(j) + (v(k,j,i+1) - v(k,j,i)) * rdx(i)   )
    !                 tauxyw = gw * ( (u(k,j+1,i-1) - u(k,j,i-1)) * rdy(j) + (v(k,j,i) - v(k,j,i-1)) * rdx(i-1) )

    !                 ! Normal stresses
    !                 !             ---------------inner derivatives---------------
    !                 tauyyn = gn * 2.0_realk * (v(k,j+1,i) - v(k,j,i)) * rddy(j+1)
    !                 tauyys = gs * 2.0_realk * (v(k,j,i) - v(k,j-1,i)) * rddy(j)
                    
    !                 ! Shear stresses
    !                 !             ------------------------------inner derivatives------------------------------
    !                 tauzyt = gt * ( (v(k+1,j,i) - v(k,j,i)) * rdz(k)   + (w(k,j+1,i) - w(k,j,i))     * rdy(j) )
    !                 tauzyb = gb * ( (v(k,j,i) - v(k-1,j,i)) * rdz(k-1) + (w(k-1,j+1,i) - w(k-1,j,i)) * rdy(j) )

    !                 ! Change due to diffusion
    !                 !                                  ---------------------------------------outer derivatives----------------------------------------
    !                 dvo = 1/densityFieldjStag(k,j,i) * ( ( tauxye - tauxyw ) * rddx(i) + ( tauyyn - tauyys ) * rdy(j) + ( tauzyt - tauzyb ) * rddz(k) )

    !                 ! Addition
    !                 vo(k, j, i) = vo(k, j, i) + dvo
    !             END DO
    !         END DO
    !     END DO

    !     DO i = 3, ii-2
    !         DO j = 3, jj-2
    !             DO k = 3-nbw, kk-3+ntw
    !                 ! Face values of dynamic viscosity on w-momentum cell
    !                 ! Harmonic mean for a more physical treatment at interfaces
    !                 ge = g(k, j, i)*g(k, j, i+1) &
    !                     / MAX(g(k, j, i) + g(k, j, i+1), MIN(gmol1,gmol2)) &
    !                     + g(k+1, j, i)*g(k+1, j, i+1) &
    !                     / MAX(g(k+1, j, i) + g(k+1, j, i+1), MIN(gmol1,gmol2))
    !                 gw = g(k, j, i-1)*g(k, j, i) &
    !                     / MAX(g(k, j, i-1) + g(k, j, i), MIN(gmol1,gmol2)) &
    !                     + g(k+1, j, i-1)*g(k+1, j, i) &
    !                     / MAX(g(k+1, j, i-1) + g(k+1, j, i), MIN(gmol1,gmol2))
    !                 gn = g(k, j, i)*g(k, j+1, i) &
    !                     / MAX(g(k, j, i) + g(k, j+1, i), MIN(gmol1,gmol2)) &
    !                     + g(k+1, j, i)*g(k+1, j+1, i) &
    !                     / MAX(g(k+1, j, i) + g(k+1, j+1, i), MIN(gmol1,gmol2))
    !                 gs = g(k, j-1, i)*g(k, j, i) &
    !                     / MAX(g(k, j-1, i) + g(k, j, i), MIN(gmol1,gmol2)) &
    !                     + g(k+1, j-1, i)*g(k+1, j, i) &
    !                     / MAX(g(k+1, j-1, i) + g(k+1, j, i), MIN(gmol1,gmol2))
    !                 gt = g(k+1, j, i)
    !                 gb = g(k, j, i)

    !                 ! Shear stresses
    !                 !             ------------------------------inner derivatives------------------------------
    !                 tauxze = ge * ( (u(k+1,j,i) - u(k,j,i))     * rdz(k) + (w(k,j,i+1) - w(k,j,i)) * rdx(i)   )
    !                 tauxzw = gw * ( (u(k+1,j,i-1) - u(k,j,i-1)) * rdz(k) + (w(k,j,i) - w(k,j,i-1)) * rdx(i-1) )
    !                 tauyzn = gn * ( (v(k+1,j,i) - v(k,j,i))     * rdz(k) + (w(k,j+1,i) - w(k,j,i)) * rdy(j)   )
    !                 tauyzs = gs * ( (v(k+1,j-1,i) - v(k,j-1,i)) * rdz(k) + (w(k,j,i) - w(k,j-1,i)) * rdy(j-1) )
                    
    !                 ! Normal stresses
    !                 !             ---------------inner derivatives---------------
    !                 tauzzt = gt * 2.0_realk * (w(k+1,j,i) - w(k,j,i)) * rddz(k+1)
    !                 tauzzb = gb * 2.0_realk * (w(k,j,i) - w(k-1,j,i)) * rddz(k)

    !                 ! Change due to diffusion
    !                 !                                  ---------------------------------------outer derivatives----------------------------------------
    !                 dwo = 1/densityFieldkStag(k,j,i) * ( ( tauxze - tauxzw ) * rddx(i) + ( tauyzn - tauyzs ) * rddy(j) + ( tauzzt - tauzzb ) * rdz(k) )

    !                 ! Addition
    !                 wo(k, j, i) = wo(k, j, i) + dwo
    !             END DO
    !         END DO
    !     END DO 

    ! END SUBROUTINE multiphase_momentum_diffusion

    !================================================================

    ! PURE SUBROUTINE comp_adve_quick(kk, jj, ii, adveField, &
    !     advrE, advrW, advrN, advrS, advrT, advrB, &
    !     adveE, adveW, adveN, adveS, adveT, adveB)
    ! !----------------------------------------------------------------
    ! !   What it does:
    ! !   QUICK interpolation to compute the advected veloctiy
    ! !   (advectee) on staggered grid cells.
    ! !   adve = advected q (advectee)
    ! !   advr = advecting q (advector)
    ! !   An indicator function is used to avoid if-statements within
    ! !   loops.
    ! !   
    ! !   Source: 
    ! !   T. Arrufat et al., “A mass-momentum consistent, 
    ! !   Volume-of-Fluid method for incompressible flow on staggered 
    ! !   grids,” Computers & Fluids, vol. 215, p. 104785, Jan. 2021, 
    ! !   doi: 10.1016/j.compfluid.2020.104785.
    ! !----------------------------------------------------------------

    !     ! Subroutine arguments
    !     INTEGER(intk), INTENT(in) :: kk, jj, ii
    !     REAL(realk), INTENT(in) :: adveField(kk, jj, ii)
    !     REAL(realk), INTENT(in) :: advrE(kk, jj, ii), advrW(kk, jj, ii)
    !     REAL(realk), INTENT(in) :: advrN(kk, jj, ii), advrS(kk, jj, ii)
    !     REAL(realk), INTENT(in) :: advrT(kk, jj, ii), advrB(kk, jj, ii)
    !     REAL(realk), INTENT(out) :: adveE(kk, jj, ii), adveW(kk, jj, ii)
    !     REAL(realk), INTENT(out) :: adveN(kk, jj, ii), adveS(kk, jj, ii)
    !     REAL(realk), INTENT(out) :: adveT(kk, jj, ii), adveB(kk, jj, ii)

    !     ! Loval variables
    !     INTEGER(intk) :: k, j, i
    !     REAL(realk) :: indE(2), indW(2), indN(2), indS(2), indT(2), indB(2)

    !     DO i = 3, ii-2
    !         DO j = 3, jj-2
    !             DO k = 3, kk-2
    !                 indE(1) = merge(1.0_realk, 0.0_realk, advrE(k,j,i) >= 0.0_realk)
    !                 indW(1) = merge(1.0_realk, 0.0_realk, advrW(k,j,i) >= 0.0_realk)
    !                 indN(1) = merge(1.0_realk, 0.0_realk, advrN(k,j,i) >= 0.0_realk)
    !                 indS(1) = merge(1.0_realk, 0.0_realk, advrS(k,j,i) >= 0.0_realk)
    !                 indT(1) = merge(1.0_realk, 0.0_realk, advrT(k,j,i) >= 0.0_realk)
    !                 indB(1) = merge(1.0_realk, 0.0_realk, advrB(k,j,i) >= 0.0_realk)

    !                 indE(2) = 1.0_realk - indE(1)
    !                 indW(2) = 1.0_realk - indW(1)
    !                 indN(2) = 1.0_realk - indN(1)
    !                 indS(2) = 1.0_realk - indS(1)
    !                 indT(2) = 1.0_realk - indT(1)
    !                 indB(2) = 1.0_realk - indB(1)

    !                 adveE(k,j,i) = indE(1) * ( 0.75_realk * adveField(k,j,i) + 0.375_realk * adveField(k,j,i+1) - 0.125_realk * adveField(k,j,i-1) ) + &
    !                                indE(2) * ( 0.75_realk * adveField(k,j,i+1) + 0.375_realk * adveField(k,j,i) - 0.125_realk * adveField(k,j,i+2) )
    !                 adveW(k,j,i) = indW(1) * ( 0.75_realk * adveField(k,j,i-1) + 0.375_realk * adveField(k,j,i) - 0.125_realk * adveField(k,j,i-2) ) + &
    !                                indW(2) * ( 0.75_realk * adveField(k,j,i) + 0.375_realk * adveField(k,j,i-1) - 0.125_realk * adveField(k,j,i+1) )
    !                 adveN(k,j,i) = indN(1) * ( 0.75_realk * adveField(k,j,i) + 0.375_realk * adveField(k,j+1,i) - 0.125_realk * adveField(k,j-1,i) ) + &
    !                                indN(2) * ( 0.75_realk * adveField(k,j+1,i) + 0.375_realk * adveField(k,j,i) - 0.125_realk * adveField(k,j+2,i) )
    !                 adveS(k,j,i) = indS(1) * ( 0.75_realk * adveField(k,j-1,i) + 0.375_realk * adveField(k,j,i) - 0.125_realk * adveField(k,j-2,i) ) + &
    !                                indS(2) * ( 0.75_realk * adveField(k,j,i) + 0.375_realk * adveField(k,j-1,i) - 0.125_realk * adveField(k,j+1,i) )
    !                 adveT(k,j,i) = indT(1) * ( 0.75_realk * adveField(k,j,i) + 0.375_realk * adveField(k+1,j,i) - 0.125_realk * adveField(k-1,j,i) ) + &
    !                                indT(2) * ( 0.75_realk * adveField(k+1,j,i) + 0.375_realk * adveField(k,j,i) - 0.125_realk * adveField(k+2,j,i) )
    !                 adveB(k,j,i) = indB(1) * ( 0.75_realk * adveField(k-1,j,i) + 0.375_realk * adveField(k,j,i) - 0.125_realk * adveField(k-2,j,i) ) + &
    !                                indB(2) * ( 0.75_realk * adveField(k,j,i) + 0.375_realk * adveField(k-1,j,i) - 0.125_realk * adveField(k+1,j,i) )
    !             END DO
    !         END DO
    !     END DO

    ! END SUBROUTINE comp_adve_quick

    !================================================================

    PURE SUBROUTINE comp_advr_centr(kk, jj, ii, u, v, w, &
        iStag, jStag, kStag, advrE, advrW, advrN, advrS, advrT, advrB)
    !----------------------------------------------------------------
    !   What it does:
    !   Central difference to compute the advecting velocity
    !   (advector) on staggered grid cells.
    !   advr = advecting q (advector)
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
        dx, dy, dz, ddx, ddy, ddz, i0, j0, k0, iStag, jStag, kStag, &
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
        INTEGER(intk), INTENT(out), OPTIONAL :: i0, j0, k0
        REAL(realk), INTENT(out), OPTIONAL :: iStag, jStag, kStag
        REAL(realk), INTENT(out), OPTIONAL :: deltaX(ii), deltaY(jj), deltaZ(kk)
        REAL(realk), INTENT(out), OPTIONAL :: veloFld(kk, jj, ii)

        ! Local variables
        ! None

        IF ( q == 1 ) THEN
            IF ( PRESENT(iStag) ) iStag = 1.0_realk
            IF ( PRESENT(jStag) ) jStag = 0.0_realk
            IF ( PRESENT(kStag) ) kStag = 0.0_realk
            IF ( PRESENT(i0) ) i0 = 1
            IF ( PRESENT(j0) ) j0 = 0
            IF ( PRESENT(k0) ) k0 = 0
            IF ( PRESENT(deltaX) ) deltaX = dx
            IF ( PRESENT(deltaY) ) deltaY = ddy
            IF ( PRESENT(deltaZ) ) deltaZ = ddz
            IF ( PRESENT(veloFld) ) veloFld = u
        ELSE IF ( q == 2 ) THEN
            IF ( PRESENT(iStag) ) iStag = 0.0_realk
            IF ( PRESENT(jStag) ) jStag = 1.0_realk
            IF ( PRESENT(kStag) ) kStag = 0.0_realk
            IF ( PRESENT(i0) ) i0 = 0
            IF ( PRESENT(j0) ) j0 = 1
            IF ( PRESENT(k0) ) k0 = 0
            IF ( PRESENT(deltaX) ) deltaX = ddx
            IF ( PRESENT(deltaY) ) deltaY = dy
            IF ( PRESENT(deltaZ) ) deltaZ = ddz
            IF ( PRESENT(veloFld) ) veloFld = v
        ELSE IF ( q == 3 ) THEN
            IF ( PRESENT(iStag) ) iStag = 0.0_realk
            IF ( PRESENT(jStag) ) jStag = 0.0_realk
            IF ( PRESENT(kStag) ) kStag = 1.0_realk
            IF ( PRESENT(i0) ) i0 = 0
            IF ( PRESENT(j0) ) j0 = 0
            IF ( PRESENT(k0) ) k0 = 1
            IF ( PRESENT(deltaX) ) deltaX = ddx
            IF ( PRESENT(deltaY) ) deltaY = ddy
            IF ( PRESENT(deltaZ) ) deltaZ = dz
            IF ( PRESENT(veloFld) ) veloFld = w
        END IF
    
    END SUBROUTINE get_component_specifics

END MODULE multiphase_vof_transport_mod
