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
    USE multiphase_plic_mod, ONLY: comp_prop, iface_recon_wrap, comp_stag_frac_wrap
    USE rungekutta_mod, ONLY: rk_2n_t
    USE multiphasecore_mod, ONLY: gmol1, gmol2, rho1, rho2, splitting_multiphase, permutation_multiphase
    USE multiphase_material_mod, ONLY: comp_material_property_field
    
    IMPLICIT NONE
    PRIVATE 

    PUBLIC :: init_multiphase_vof_transport, finish_multiphase_vof_transport, multiphase_solve

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

    SUBROUTINE comp_flux_cent(kk, jj, ii, splitDir, field, isInterface, u, v, w, alpha, dt, normx, normy, normz, ddx, ddy, ddz, tol, fieldFlux)
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

                                CALL comp_prop(fluxedProp, fluxAlpha, field(k,j,i), &
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

                                CALL comp_prop(fluxedProp, fluxAlpha, field(k,j,i+1), fluxWidth, ddy(j), ddz(k), normx(k,j,i+1), normy(k,j,i+1), normz(k,j,i+1), tol)

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

                                CALL comp_prop(fluxedProp, fluxAlpha, field(k,j,i), ddx(i), fluxWidth, ddz(k), normx(k,j,i), normy(k,j,i), normz(k,j,i), tol)

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

                                CALL comp_prop(fluxedProp, fluxAlpha, field(k,j+1,i), ddx(i), fluxWidth, ddz(k), normx(k,j+1,i), normy(k,j+1,i), normz(k,j+1,i), tol)

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

                                CALL comp_prop(fluxedProp, fluxAlpha, field(k,j,i), ddx(i), ddy(j), fluxWidth, normx(k,j,i), normy(k,j,i), normz(k,j,i), tol)

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

                                CALL comp_prop(fluxedProp, fluxAlpha, field(k+1,j,i), ddx(i), ddy(j), fluxWidth, normx(k+1,j,i), normy(k+1,j,i), normz(k+1,j,i), tol)

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

    END SUBROUTINE comp_flux_cent

    !================================================================

    SUBROUTINE comp_flux_stag(kk, jj, ii, q, splitDir, field, isInterface, u, v, w, alpha, dt, normx, normy, normz, dx, dy, dz, ddx, ddy, ddz, tol, fieldFlux, complementFieldFlux)
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
        REAL(realk), INTENT(inout) :: fieldFlux(kk, jj, ii), complementFieldFlux(kk, jj, ii)

        ! Local variables
        REAL(realk) :: deltaX(ii), deltaY(jj), deltaZ(kk)
        INTEGER(intk) :: k, j, i
        REAL(realk) :: flux, complementFlux, fluxedProp, fluxWidth, fluxAlpha
        REAL(realk) :: advrE(kk, jj, ii), advrN(kk, jj, ii), advrT(kk, jj, ii)

        CALL get_component_specifics(kk, jj, ii, q, dx=dx, dy=dy, dz=dz, ddx=ddx, ddy=ddy, ddz=ddz, deltaX=deltaX, deltaY=deltaY, deltaZ=deltaZ)
        CALL comp_advr_centr(kk, jj, ii, q, u, v, w, advrE, advrN, advrT)

        IF ( splitDir == 1 ) THEN
            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        IF ( advrE(k,j,i) > tol ) THEN
                            IF ( isInterface(k,j,i) ) THEN
                                ! See explanation above.
                                fluxWidth = abs( advrE(k,j,i) ) * dt
                                fluxAlpha = alpha(k,j,i) - normx(k,j,i) * ( deltaX(i)/2.0_realk - fluxWidth )

                                CALL comp_prop(fluxedProp, fluxAlpha, field(k,j,i), fluxWidth, ddy(j), ddz(k), normx(k,j,i), normy(k,j,i), normz(k,j,i), tol)
                                
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

                                CALL comp_prop(fluxedProp, fluxAlpha, field(k,j,i+1), fluxWidth, ddy(j), ddz(k), normx(k,j,i+1), normy(k,j,i+1), normz(k,j,i+1), tol)

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

                                CALL comp_prop(fluxedProp, fluxAlpha, field(k,j,i), ddx(i), fluxWidth, ddz(k), normx(k,j,i), normy(k,j,i), normz(k,j,i), tol)

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

                                CALL comp_prop(fluxedProp, fluxAlpha, field(k,j+1,i), ddx(i), fluxWidth, ddz(k), normx(k,j+1,i), normy(k,j+1,i), normz(k,j+1,i), tol)

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

                                CALL comp_prop(fluxedProp, fluxAlpha, field(k,j,i), ddx(i), ddy(j), fluxWidth, normx(k,j,i), normy(k,j,i), normz(k,j,i), tol)

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

                                CALL comp_prop(fluxedProp, fluxAlpha, field(k+1,j,i), ddx(i), ddy(j), fluxWidth, normx(k+1,j,i), normy(k+1,j,i), normz(k+1,j,i), tol)

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
    
    PURE SUBROUTINE def_advection_sequence(iteration, advSeq)
    !----------------------------------------------------------------
    !   What it does:
    !   Defines the sequence, in which the volume fraction field and
    !   the momentum are advected. 
    !
    !   The PARIS solver uses a cyclic periodicity of three. Hence,
    !   the sequence can be
    !   x -> y -> z,
    !   y -> z -> x or
    !   z -> x -> y.
    !   
    !   Permutation can also be see as the sum of all possible 
    !   sequences. For three dimensions we get six sequences.
    !   x -> y -> z,
    !   y -> z -> x,
    !   z -> x -> y,
    !   x -> z -> y,
    !   y -> x -> z or
    !   z -> y -> x.
    !
    !   The permutation_multiphase variable controls which version 
    !   is used.
    !
    !   Source:
    !   T. Arrufat et al., “A mass-momentum consistent, 
    !   Volume-of-Fluid method for incompressible flow on staggered 
    !   grids,” Computers & Fluids, vol. 215, p. 104785, Jan. 2021, 
    !   doi: 10.1016/j.compfluid.2020.104785.
    !
    !   W. Aniszewski et al., “PArallel, Robust, Interface Simulator
    !   (PARIS),” Computer Physics Communications, vol. 263, 
    !   p. 107849, Jun. 2021, doi: 10.1016/j.cpc.2021.107849.
    !   
    !   PARIS source code (accessed: Mai 2026)
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: iteration
        INTEGER(intk), INTENT(inout) :: advSeq(3)

        ! Local variables
        INTEGER(intk) :: permutationIndex

        ! permutationIndex only changes in a new time-step
        permutationIndex = mod(iteration-1, permutation_multiphase)

        ! Select permutation of split advection
        SELECT CASE (permutationIndex)
            CASE (0)
                advSeq = [1, 2, 3]
            CASE (1)
                advSeq = [3, 1, 2]
            CASE (2)
                advSeq = [2, 3, 1]
            CASE (3)
                advSeq = [1, 3, 2]
            CASE (4)
                advSeq = [3, 2, 1]
            CASE (5)
                advSeq = [2, 1, 3]
        END SELECT

    END SUBROUTINE def_advection_sequence

    !================================================================

    SUBROUTINE multiphase_solve(u_f, v_f, w_f, f_f, p_f, g_f, d_f, dtrki, itstep, uo_f, vo_f, wo_f)
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
        TYPE(field_t), INTENT(inout) :: uo_f, vo_f, wo_f

        ! Local variables
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: u, v, w
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: vff, p, g, d
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: uo, vo, wo

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

        REAL(realk), ALLOCATABLE :: vffPrev(:,:,:)

        INTEGER(intk) :: i, igrid
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: nfro, nbac, nrgt, nlft, nbot, ntop
        REAL(realk), PARAMETER :: tol = 1.0E-12_realk

        uo_f = 0.0_realk
        vo_f = 0.0_realk
        wo_f = 0.0_realk

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
            CALL uo_f%get_ptr(uo, igrid)
            CALL vo_f%get_ptr(vo, igrid)
            CALL wo_f%get_ptr(wo, igrid)

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

            IF ( .NOT. ALLOCATED(vffPrev) ) ALLOCATE(vffPrev(kk, jj, ii))
            vffPrev = vff

            CALL adve_operator(kk, jj, ii, u, v, w, vff, dx, dy, dz, ddx, ddy, ddz, dtrki, itstep, tol, &
                normx, normy, normz, alpha, &
                vffiStag, vffjStag, vffkStag, diStag, djStag, dkStag, &
                normxiStag, normyiStag, normziStag, & 
                normxjStag, normyjStag, normzjStag, &
                normxkStag, normykStag, normzkStag, &
                alphaiStag, alphajStag, alphakStag, uo, vo, wo)

            CALL diff_operator(kk, jj, ii, u, v, w, vffPrev, g, rdx, rdy, rdz, rddx, rddy, rddz, uo, vo, wo)

            ! CALL exte_operator()

            u = u + uo * dtrki
            v = v + vo * dtrki
            w = w + wo * dtrki

        END DO

    END SUBROUTINE multiphase_solve

    !================================================================

    SUBROUTINE adve_operator(kk, jj, ii, u, v, w, vff, dx, dy, dz, ddx, ddy, ddz, dtrki, itstep, tol, &
        normx, normy, normz, alpha, &
        vffiStag, vffjStag, vffkStag, diStag, djStag, dkStag, &
        normxiStag, normyiStag, normziStag, & 
        normxjStag, normyjStag, normzjStag, &
        normxkStag, normykStag, normzkStag, &
        alphaiStag, alphajStag, alphakStag, uo, vo, wo)
    !----------------------------------------------------------------
    !   What it does:
    !    
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(inout) :: vff(kk, jj, ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: tol, dtrki
        INTEGER(intk), INTENT(in) :: itstep
        REAL(realk), INTENT(out) :: normx(kk, jj, ii), normy(kk, jj, ii), normz(kk, jj, ii), alpha(kk, jj, ii)
        REAL(realk), INTENT(out) :: vffiStag(kk, jj, ii), vffjStag(kk, jj, ii), vffkStag(kk, jj, ii)
        REAL(realk), INTENT(out) :: diStag(kk, jj, ii), djStag(kk, jj, ii), dkStag(kk, jj, ii)
        REAL(realk), INTENT(out) :: normxiStag(kk, jj, ii), normyiStag(kk, jj, ii), normziStag(kk, jj, ii)
        REAL(realk), INTENT(out) :: normxjStag(kk, jj, ii), normyjStag(kk, jj, ii), normzjStag(kk, jj, ii)
        REAL(realk), INTENT(out) :: normxkStag(kk, jj, ii), normykStag(kk, jj, ii), normzkStag(kk, jj, ii)
        REAL(realk), INTENT(out) :: alphaiStag(kk, jj, ii), alphajStag(kk, jj, ii), alphakStag(kk, jj, ii)
        REAL(realk), INTENT(inout) :: uo(kk, jj, ii), vo(kk, jj, ii), wo(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: q, advSeq(3), l, splitDir
        REAL(realk), ALLOCATABLE :: vffStag(:,:,:)
        REAL(realk), ALLOCATABLE :: dStag(:,:,:)
        REAL(realk), ALLOCATABLE :: normxStag(:,:,:), normyStag(:,:,:), normzStag(:,:,:)
        REAL(realk), ALLOCATABLE :: alphaStag(:,:,:)
        LOGICAL, ALLOCATABLE :: isInterface(:,:,:)
        LOGICAL, ALLOCATABLE :: isInterfaceStag(:,:,:)
        LOGICAL, ALLOCATABLE :: isNearInterface(:,:,:)
        LOGICAL, ALLOCATABLE :: isNearInterfaceStag(:,:,:)
        REAL(realk), ALLOCATABLE :: vffFlux(:,:,:), complVffFlux(:,:,:), mom(:,:,:), cWY(:,:,:), cWYStag(:,:,:)
        REAL(realk), ALLOCATABLE :: advrE(:,:,:), advrN(:,:,:), advrT(:,:,:), advrSplitDir(:,:,:,:)
        REAL(realk), ALLOCATABLE :: uNew(:,:,:), vNew(:,:,:), wNew(:,:,:)
        REAL(realk), ALLOCATABLE :: mom4D(:,:,:,:), vffStag4D(:,:,:,:), cWYStag4D(:,:,:,:)
        
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
        IF (.NOT. ALLOCATED(complVffFlux))        ALLOCATE(complVffFlux(kk,jj,ii))
        IF (.NOT. ALLOCATED(mom))                 ALLOCATE(mom(kk,jj,ii))
        IF (.NOT. ALLOCATED(cWY))                 ALLOCATE(cWY(kk,jj,ii))
        IF (.NOT. ALLOCATED(cWYStag))             ALLOCATE(cWYStag(kk,jj,ii))
        IF (.NOT. ALLOCATED(advrE))               ALLOCATE(advrE(kk,jj,ii))
        IF (.NOT. ALLOCATED(advrN))               ALLOCATE(advrN(kk,jj,ii))
        IF (.NOT. ALLOCATED(advrT))               ALLOCATE(advrT(kk,jj,ii))
        IF (.NOT. ALLOCATED(advrSplitDir))        ALLOCATE(advrSplitDir(kk,jj,ii,3))
        IF (.NOt. ALLOCATED(uNew))                ALLOCATE(uNew(kk,jj,ii))
        IF (.NOt. ALLOCATED(vNew))                ALLOCATE(vNew(kk,jj,ii))
        IF (.NOt. ALLOCATED(wNew))                ALLOCATE(wNew(kk,jj,ii))
        IF (.NOt. ALLOCATED(mom4D))               ALLOCATE(mom4D(kk,jj,ii,3))
        IF (.NOt. ALLOCATED(vffStag4D))           ALLOCATE(vffStag4D(kk,jj,ii,3))
        IF (.NOT. ALLOCATED(cWYStag4D))           ALLOCATE(cWYStag4D(kk,jj,ii,3))
        
        IF ( splitting_multiphase == "component-wise" ) THEN

            CALL check_solenoidality(kk, jj, ii, u, v, w, dx, dy, dz)
            CALL def_advection_sequence(itstep, advSeq)
            CALL iface_recon_wrap(kk, jj, ii, vff, ddx, ddy, ddz, tol, normx, normy, normz, alpha, isInterface, isNearInterface)

            DO q = 1, 3
                
                CALL comp_stag_frac_wrap(kk, jj, ii, q, alpha, vff, isInterface, ddx, ddy, ddz, normx, normy, normz, tol, vffStag)
                CALL iface_recon_wrap(kk, jj, ii, q, vffStag, dx, dy, dz, ddx, ddy, ddz, tol, normxStag, normyStag, normzStag, alphaStag, isInterfaceStag, isNearInterfaceStag)
                CALL comp_material_property_field(kk, jj, ii, vffStag, rho1, rho2, dStag)
                CALL comp_momentum(kk, jj, ii, q, dStag, u, v, w, mom)
                CALL comp_cWY(kk, jj, ii, vffStag, cWYStag)
                CALL comp_advr_centr(kk, jj, ii, q, u, v, w, advrE, advrN, advrT)
                advrSplitDir(:,:,:,1) = advrE
                advrSplitDir(:,:,:,2) = advrN
                advrSplitDir(:,:,:,3) = advrT

                DO l = 1, 3

                    splitDir = advSeq(l)
                                            
                    CALL comp_flux_stag(kk, jj, ii, q, splitDir, vff, isInterface, u, v, w, alpha, dtrki, normx, normy, normz, dx, dy, dz, ddx, ddy, ddz, tol, vffFlux, complVffFlux)
                    CALL adv_mom(kk, jj, ii, splitDir, vffStag, vffFlux, complVffFlux, cWYStag, advrSplitDir(:,:,:,splitDir), dx, dy, dz, ddx, ddy, ddz, dtrki, mom)
                    CALL adv_vof(kk, jj, ii, splitDir, vffFlux, cWYStag, advrE, advrN, advrT, dx, dy, dz, ddx, ddy, ddz, dtrki, tol, vffStag)

                END DO
                
                IF ( q == 1 ) THEN
                    CALL comp_velocity_change(kk, jj, ii, q, u, vffStag, mom, dtrki, uo)
                ELSE IF ( q == 2 ) THEN
                    CALL comp_velocity_change(kk, jj, ii, q, v, vffStag, mom, dtrki, vo)
                ELSE IF ( q == 3 ) THEN
                    CALL comp_velocity_change(kk, jj, ii, q, w, vffStag, mom, dtrki, wo)
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

            DO l = 1, 3

                splitDir = advSeq(l)
                
                CALL iface_recon_wrap(kk, jj, ii, vff, ddx, ddy, ddz, tol, normx, normy, normz, alpha, isInterface, isNearInterface)
                CALL comp_flux_cent(kk, jj, ii, splitDir, vff, isInterface, u, v, w, alpha, dtrki, normx, normy, normz, ddx, ddy, ddz, tol, vffFlux)
                CALL adv_vof(kk, jj, ii, splitDir, vffFlux, cWY, u, v, w, dx, dy, dz, ddx, ddy, ddz, dtrki, tol, vff)
                CALL clip_vff(kk, jj, ii, tol, vff)

            END DO

        ELSE IF ( splitting_multiphase == "direction-wise" ) THEN
            
            CALL check_solenoidality(kk, jj, ii, u, v, w, dx, dy, dz)
            CALL def_advection_sequence(itstep, advSeq)
            CALL clip_vff(kk, jj, ii, tol, vff)
            CALL iface_recon_wrap(kk, jj, ii, vff, ddx, ddy, ddz, tol, normx, normy, normz, alpha, isInterface, isNearInterface)

            DO q = 1, 3

                CALL comp_stag_frac_wrap(kk, jj, ii, q, alpha, vff, isInterface, ddx, ddy, ddz, normx, normy, normz, tol, vffStag)
                CALL iface_recon_wrap(kk, jj, ii, vffStag, ddx, ddy, ddz, tol, normxStag, normyStag, normzStag, alphaStag, isInterfaceStag, isNearInterfaceStag)
                CALL comp_material_property_field(kk, jj, ii, vffStag, rho1, rho2, dStag)
                CALL comp_momentum(kk, jj, ii, q, dStag, u, v, w, mom)
                CALL comp_cWY(kk, jj, ii, vffStag, cWYStag)

                vffStag4D(:,:,:,q) = vffStag
                mom4D(:,:,:,q) = mom
                cWYStag4D(:,:,:,q) = cWYStag
                
            END DO

            CALL comp_cWY(kk, jj, ii, vff, cWY)

            DO l = 1, 3

                splitDir = advSeq(l)

                DO q = 1, 3
                    
                    vffStag = vffStag4D(:,:,:,q)
                    mom = mom4D(:,:,:,q)
                    cWYStag = cWYStag4D(:,:,:,q)

                    CALL comp_advr_centr(kk, jj, ii, q, u, v, w, advrE, advrN, advrT)
                    advrSplitDir(:,:,:,1) = advrE
                    advrSplitDir(:,:,:,2) = advrN
                    advrSplitDir(:,:,:,3) = advrT

                    CALL comp_flux_stag(kk, jj, ii, q, splitDir, vff, isInterface, u, v, w, alpha, dtrki, normx, normy, normz, dx, dy, dz, ddx, ddy, ddz, tol, vffFlux, complVffFlux)
                    CALL adv_mom(kk, jj, ii, splitDir, vffStag, vffFlux, complVffFlux, cWYStag, advrSplitDir(:,:,:,splitDir), dx, dy, dz, ddx, ddy, ddz, dtrki, mom)
                    CALL adv_vof(kk, jj, ii, splitDir, vffFlux, cWYStag, advrE, advrN, advrT, dx, dy, dz, ddx, ddy, ddz, dtrki, tol, vffStag)

                    IF ( q == 1 ) THEN
                        CALL comp_velocity_change(kk, jj, ii, q, u, vffStag, mom, dtrki, uo)
                    ELSE IF ( q == 2 ) THEN
                        CALL comp_velocity_change(kk, jj, ii, q, v, vffStag, mom, dtrki, vo)
                    ELSE IF ( q == 3 ) THEN
                        CALL comp_velocity_change(kk, jj, ii, q, w, vffStag, mom, dtrki, wo)
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

                CALL iface_recon_wrap(kk, jj, ii, vff, ddx, ddy, ddz, tol, normx, normy, normz, alpha, isInterface, isNearInterface)
                CALL comp_flux_cent(kk, jj, ii, splitDir, vff, isInterface, u, v, w, alpha, dtrki, normx, normy, normz, ddx, ddy, ddz, tol, vffFlux)
                CALL adv_vof(kk, jj, ii, splitDir, vffFlux, cWY, u, v, w, dx, dy, dz, ddx, ddy, ddz, dtrki, tol, vff)

            END DO

        END IF

    END SUBROUTINE

    !================================================================

    SUBROUTINE adv_vof(kk, jj, ii, splitDir, vffFlux, cWY, velWY1, velWY2, velWY3, dx, dy, dz, ddx, ddy, ddz, dtrki, tol, vff)
    !----------------------------------------------------------------
    !   What it does:
    !    
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: splitDir
        REAL(realk), INTENT(in) :: vffFlux(kk, jj, ii)
        REAL(realk), INTENT(in) :: cWY(kk, jj, ii), velWY1(kk, jj, ii), velWY2(kk, jj, ii), velWY3(kk, jj, ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: dtrki, tol
        REAL(realk), INTENT(inout) :: vff(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: i, j, k
        INTEGER(intk) :: i0, j0, k0
        REAL(realk) :: deltaX(ii), deltaY(jj), deltaZ(kk)
        REAL(realk) :: velWY(kk, jj, ii)
        REAL(realk) :: div(kk, jj, ii)
        
        i0 = 0_intk; j0 = 0_intk; k0 = 0_intk

        CALL get_component_specifics(kk, jj, ii, splitDir, vel1=velWY1, vel2=velWY2, vel3=velWY3, dx=dx, dy=dy, dz=dz, ddx=ddx, ddy=ddy, ddz=ddz, &
            i0=i0, j0=j0, k0=k0, deltaX=deltaX, deltaY=deltaY, deltaZ=deltaZ, velFld=velWY)

        DO i = 2, ii-1
            DO j = 2, jj-1
                DO k = 2, kk-1
                    div(k,j,i) = ( velWY(k,j,i) - velWY(k-k0,j-j0,i-i0) ) / deltaX(ii)
                    vff(k,j,i) = vff(k,j,i) - dtrki * ( vffFlux(k,j,i) - vffFlux(k-k0,j-j0,i-i0) ) + dtrki * cWY(k,j,i) * div(k,j,i)
                END DO
            END DO
        END DO

    END SUBROUTINE adv_vof

    !================================================================

    SUBROUTINE adv_mom(kk, jj, ii, splitDir, vff, vffFlux, complVffFlux, cWY, velWY, dx, dy, dz, ddx, ddy, ddz, dtrki, mom)
    !----------------------------------------------------------------
    !   What it does:
    !    
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: splitDir
        REAL(realk), INTENT(in) :: vff(kk, jj, ii), vffFlux(kk, jj, ii), complVffFlux(kk, jj, ii)
        REAL(realk), INTENT(in) :: cWY(kk, jj, ii), velWY(kk, jj, ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: dtrki
        REAL(realk), INTENT(inout) :: mom(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: i, j, k
        INTEGER(intk) :: i0, j0, k0
        REAL(realk) :: deltaX(ii), deltaY(jj), deltaZ(kk)
        REAL(realk) :: dStag(kk, jj, ii)
        REAL(realk) :: vel1, vel2, vel3, a1, a2, advrU, advrD
        REAL(realk) :: momFlux(kk, jj, ii)
        REAL(realk) :: div(kk, jj, ii)

        i0 = 0_intk; j0 = 0_intk; k0 = 0_intk

        CALL get_component_specifics(kk, jj, ii, splitDir, dx=dx, dy=dy, dz=dz, ddx=ddx, ddy=ddy, ddz=ddz, &
            i0=i0, j0=j0, k0=k0, deltaX=deltaX, deltaY=deltaY, deltaZ=deltaZ)

        CALL comp_material_property_field(kk, jj, ii, vff, rho1, rho2, dStag)

        DO i = 2, ii-1
            DO j = 2, jj-1
                DO k = 2, kk-1
                    vel1 = mom(k-k0,j-j0,i-i0) / dStag(k-k0,j-j0,i-i0)
                    vel2 = mom(k,j,i) / dStag(k,j,i)
                    vel3 = mom(k+k0,j+j0,i+i0) / dStag(k+k0,j+j0,i+i0)

                    a1 = velWY(k-k0,j-j0,i-i0)*dtrki/deltaX(i-i0)
                    a2 = velWY(k,j,i)*dtrki/deltaX(i)

                    CALL comp_advr_inter(vel1, vel2, vel3, -0.5_realk*(1.0_realk + a1), "ENO", advrU)
                    CALL comp_advr_inter(vel1, vel2, vel3,  0.5_realk*(1.0_realk - a2), "ENO", advrD)
                    
                    momFlux(k,j,i) = ( rho1 * vffFlux(k,j,i) + rho2 * complVffFlux(k,j,i) ) * advrD
                END DO
            END DO 
        END DO

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    div(k,j,i) = ( velWY(k,j,i) - velWY(k-k0,j-j0,i-i0) ) / deltaX(ii)
                    mom(k,j,i) = mom(k,j,i) - dtrki * ( momFlux(k,j,i) - momFlux(k-k0,j-j0,i-i0) ) + dtrki * (rho1-rho2) * velWY(k,j,i) * cWY(k,j,i) * div(k,j,i)
                END DO
            END DO
        END DO

    END SUBROUTINE adv_mom

    !================================================================

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
        REAL(realk) :: velFld(kk, jj, ii)

        CALL get_component_specifics(kk, jj, ii, q, vel1=u, vel2=v, vel3=w, velFld=velFld)

        DO i = 1, ii
            DO j = 1, jj
                DO k = 1, kk
                    mom(k,j,i) = velFld(k,j,i) * dStag(k,j,i)
                END DO
            END DO
        END DO

    END SUBROUTINE comp_momentum

    !================================================================

    SUBROUTINE comp_velocity_change(kk, jj, ii, q, vel, vff, mom, dt, velo)
    !----------------------------------------------------------------
    !   What it does:
    !    
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: q
        REAL(realk), INTENT(in) :: vel(kk, jj, ii), vff(kk, jj, ii), mom(kk, jj, ii)
        REAL(realk), INTENT(in) :: dt
        REAL(realk), INTENT(out) :: velo(kk,jj,ii)

        ! Local variables
        INTEGER(intk) :: i, j, k
        REAL(realk) :: dStag(kk, jj, ii)
        
        CALL comp_material_property_field(kk, jj, ii, vff, rho1, rho2, dStag)
    
        DO i = 4, ii-3
            DO j = 4, jj-3
                DO k = 4, kk-3
                    IF ( vff(k,j,i) >= 0.0_realk .AND. vff(k,j,i) <= 1.0_realk ) THEN
                        velo(k,j,i) = ( mom(k,j,i) / dStag(k,j,i) - vel(k,j,i) ) / dt
                    END IF
                END DO
            END DO
        END DO

    END SUBROUTINE comp_velocity_change

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

    SUBROUTINE clip_vff(kk, jj, ii, tol, vff)
    !----------------------------------------------------------------
    !   What it does:
    !   Clips the volume fraction field to its boundaries [0, 1].
    !    
    !   Source:
    !   T. Arrufat et al., “A mass-momentum consistent, 
    !   Volume-of-Fluid method for incompressible flow on staggered 
    !   grids,” Computers & Fluids, vol. 215, p. 104785, Jan. 2021, 
    !   doi: 10.1016/j.compfluid.2020.104785.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: tol
        REAL(realk), INTENT(inout) :: vff(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i
        return
        DO i = 1, ii
            DO j = 1, jj
                DO k = 1, kk
                    IF ( vff(k,j,i) <= tol ) THEN
                        vff(k,j,i) = 0.0_realk
                    ELSE IF ( vff(k,j,i) >= ( 1.0_realk - tol ) ) THEN
                        vff(k,j,i) = 1.0_realk
                    END IF
                END DO 
            END DO
        END DO

    END SUBROUTINE clip_vff

    !================================================================

    SUBROUTINE check_solenoidality(kk, jj, ii, u, v, w, dx, dy, dz)
    !----------------------------------------------------------------
    !   What it does:
    !    
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)

        ! Local variables
        INTEGER(intk) :: k, j, i
        REAL(realk) :: div(kk, jj, ii)
        LOGICAL :: isSolenoidal
        REAL(realk), PARAMETER :: eps = 1.0E-10_realk
        REAL(realk) :: uChar, lChar, L1Eps, L2Eps, LinfEps, L1Norm, L2Norm, LinfNorm

        isSolenoidal = .TRUE.

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2

                    div(k,j,i) = ( u(k,j,i+1) - u(k,j,i) )/dx(i) + &
                                 ( v(k,j+1,i) - v(k,j,i) )/dy(j) + &
                                 ( w(k+1,j,i) - w(k,j,i) )/dz(k)

                END DO
            END DO
        END DO

        uChar = MAX(MAXVAL(ABS(u(3:kk-2,3:jj-2,3:ii-2))), &
                    MAXVAL(ABS(v(3:kk-2,3:jj-2,3:ii-2))), &
                    MAXVAL(ABS(w(3:kk-2,3:jj-2,3:ii-2))))
        lChar = MAX(MAXVAL(ABS(dx(3:ii-2))), &
                    MAXVAL(ABS(dy(3:jj-2))), &
                    MAXVAL(ABS(dz(3:kk-2))))

        L1Eps = eps * uChar / lChar
        L2Eps = eps * uChar / lChar * ( SIZE(div(3:kk-2,3:jj-2,3:ii-2)) )**(1.0_realk/2.0_realk)
        LinfEps = eps * uChar / lChar * SIZE(div(3:kk-2,3:jj-2,3:ii-2))
        
        L1Norm = SUM( ABS(div(3:kk-2,3:jj-2,3:ii-2)) )
        L2Norm = ( SUM( div(3:kk-2,3:jj-2,3:ii-2)**2.0_realk ) )**(1.0_realk/2.0_realk)
        LinfNorm = MAXVAL( ABS(div(3:kk-2,3:jj-2,3:ii-2)) )

        IF ( L1Norm > L1Eps .OR. L2Norm > L2Eps .OR. LinfNorm > LinfEps ) THEN
            isSolenoidal = .FALSE.
            WRITE(*,*) "Velocity not solenoidal: L1 =", L1Norm, " L2 =", L2Norm, " Linf =", LinfNorm
        END IF

    END SUBROUTINE check_solenoidality

    !================================================================

    SUBROUTINE diff_operator(kk, jj, ii, u, v, w, vff, g, rdx, rdy, rdz, rddx, rddy, rddz, uo, vo, wo)
    !----------------------------------------------------------------
    !   What it does:
    !    
    !----------------------------------------------------------------
    
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii), vff(kk, jj, ii), g(kk, jj, ii)
        REAL(realk), INTENT(in) :: rdx(ii), rdy(jj), rdz(kk)
        REAL(realk), INTENT(in) :: rddx(ii), rddy(jj), rddz(kk)
        REAL(realk), INTENT(inout) :: uo(kk, jj, ii), vo(kk, jj, ii), wo(kk, jj, ii)

        ! Local variables
        REAL(realk) :: d(kk, jj, ii)
        INTEGER(intk) :: k, j, i
        REAL(realk) :: ge, gw, gn, gs, gt, gb
        REAL(realk) :: tauxxe, tauxxw, tauyxn, tauyxs, tauzxt, tauzxb
        REAL(realk) :: tauxye, tauxyw, tauyyn, tauyys, tauzyt, tauzyb
        REAL(realk) :: tauxze, tauxzw, tauyzn, tauyzs, tauzzt, tauzzb
        return
        CALL comp_material_property_field(kk, jj, ii, vff, rho1, rho2, d)

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    ! Face values of dynamic viscosity on u-momentum cell
                    ! Harmonic mean for a more physical treatment at interfaces
                    ge = g(k, j, i+1)
                    gw = g(k, j, i)
                    gn = g(k, j, i)*g(k, j+1, i) &
                        / MAX(g(k, j, i) + g(k, j+1, i), MIN(gmol1,gmol2)) &
                        + g(k, j, i+1)*g(k, j+1, i+1) &
                        / MAX(g(k, j, i+1) + g(k, j+1, i+1), MIN(gmol1,gmol2))
                    gs = g(k, j-1, i)*g(k, j, i) &
                        / MAX(g(k, j-1, i) + g(k, j, i), MIN(gmol1,gmol2)) &
                        + g(k, j-1, i+1)*g(k, j, i+1) &
                        / MAX(g(k, j-1, i+1) + g(k, j, i+1), MIN(gmol1,gmol2))
                    gt = g(k, j, i)*g(k+1, j, i) &
                        / MAX(g(k, j, i) + g(k+1, j, i), MIN(gmol1,gmol2)) &
                        + g(k, j, i+1)*g(k+1, j, i+1) &
                        /MAX(g(k, j, i+1) + g(k+1, j, i+1), MIN(gmol1,gmol2))
                    gb = g(k-1, j, i)*g(k, j, i) &
                        / MAX(g(k-1, j, i) + g(k, j, i), MIN(gmol1,gmol2)) &
                        + g(k-1, j, i+1)*g(k, j, i+1) &
                        / MAX(g(k-1, j, i+1) + g(k, j, i+1), MIN(gmol1,gmol2))

                    ! Normal stresses
                    tauxxe = ge * 2.0_realk * (u(k,j,i+1) - u(k,j,i)) * rddx(i+1)
                    tauxxw = gw * 2.0_realk * (u(k,j,i) - u(k,j,i-1)) * rddx(i)

                    ! Shear stresses
                    tauyxn = gn * ( (u(k,j+1,i) - u(k,j,i)) * rdy(j)   + (v(k,j,i+1) - v(k,j,i))     * rdx(i) )
                    tauyxs = gs * ( (u(k,j,i) - u(k,j-1,i)) * rdy(j-1) + (v(k,j-1,i+1) - v(k,j-1,i)) * rdx(i) )
                    tauzxt = gt * ( (u(k+1,j,i) - u(k,j,i)) * rdz(k)   + (w(k,j,i+1) - w(k,j,i))     * rdx(i) )
                    tauzxb = gb * ( (u(k,j,i) - u(k-1,j,i)) * rdz(k-1) + (w(k-1,j,i+1) - w(k-1,j,i)) * rdx(i) )

                    ! Change due to diffusion
                    uo(k,j,i) = uo(k,j,i) + 1.0_realk/( 0.5_realk * ( d(k,j,i) + d(k,j,i+1) ) ) * &
                        ( ( tauxxe - tauxxw ) * rdx(i) + &
                          ( tauyxn - tauyxs ) * rddy(j) + &
                          ( tauzxt - tauzxb ) * rddz(k) )
                END DO
            END DO
        END DO

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    ! Face values of dynamic viscosity on v-momentum cell
                    ! Harmonic mean for a more physical treatment at interfaces
                    ge = g(k, j, i)*g(k, j, i+1) &
                        / MAX(g(k, j, i) + g(k, j, i+1), MIN(gmol1,gmol2)) &
                        + g(k, j+1, i)*g(k, j+1, i+1) &
                        / MAX(g(k, j+1, i) + g(k, j+1, i+1), MIN(gmol1,gmol2))
                    gw = g(k, j, i-1)*g(k, j, i) &
                        / MAX(g(k, j, i-1) + g(k, j, i), MIN(gmol1,gmol2)) &
                        + g(k, j+1, i-1)*g(k, j+1, i) &
                        / MAX(g(k, j+1, i-1) + g(k, j+1, i), MIN(gmol1,gmol2))
                    gn = g(k, j+1, i)
                    gs = g(k, j, i)
                    gt = g(k, j, i)*g(k+1, j, i) &
                        / MAX(g(k, j, i) + g(k+1, j, i), MIN(gmol1,gmol2)) &
                        + g(k, j+1, i)*g(k+1, j+1, i) &
                        / MAX(g(k, j+1, i) + g(k+1, j+1, i), MIN(gmol1,gmol2))
                    gb = g(k-1, j, i)*g(k, j, i) &
                        / MAX(g(k-1, j, i) + g(k, j, i), MIN(gmol1,gmol2)) &
                        + g(k-1, j+1, i)*g(k, j+1, i) &
                        / MAX(g(k-1, j+1, i) + g(k, j+1, i), MIN(gmol1,gmol2))

                    ! Shear stresses
                    tauxye = ge * ( (u(k,j+1,i) - u(k,j,i))     * rdy(j) + (v(k,j,i+1) - v(k,j,i)) * rdx(i)   )
                    tauxyw = gw * ( (u(k,j+1,i-1) - u(k,j,i-1)) * rdy(j) + (v(k,j,i) - v(k,j,i-1)) * rdx(i-1) )

                    ! Normal stresses
                    tauyyn = gn * 2.0_realk * (v(k,j+1,i) - v(k,j,i)) * rddy(j+1)
                    tauyys = gs * 2.0_realk * (v(k,j,i) - v(k,j-1,i)) * rddy(j)
                    
                    ! Shear stresses
                    tauzyt = gt * ( (v(k+1,j,i) - v(k,j,i)) * rdz(k)   + (w(k,j+1,i) - w(k,j,i))     * rdy(j) )
                    tauzyb = gb * ( (v(k,j,i) - v(k-1,j,i)) * rdz(k-1) + (w(k-1,j+1,i) - w(k-1,j,i)) * rdy(j) )

                    ! Change due to diffusion
                    vo(k,j,i) = vo(k,j,i) + 1.0_realk/( 0.5_realk * ( d(k,j,i) + d(k,j+1,i) ) ) * &
                        ( ( tauxye - tauxyw ) * rddx(i) + &
                          ( tauyyn - tauyys ) * rdy(j) + &
                          ( tauzyt - tauzyb ) * rddz(k) )
                END DO
            END DO
        END DO

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    ! Face values of dynamic viscosity on w-momentum cell
                    ! Harmonic mean for a more physical treatment at interfaces
                    ge = g(k, j, i)*g(k, j, i+1) &
                        / MAX(g(k, j, i) + g(k, j, i+1), MIN(gmol1,gmol2)) &
                        + g(k+1, j, i)*g(k+1, j, i+1) &
                        / MAX(g(k+1, j, i) + g(k+1, j, i+1), MIN(gmol1,gmol2))
                    gw = g(k, j, i-1)*g(k, j, i) &
                        / MAX(g(k, j, i-1) + g(k, j, i), MIN(gmol1,gmol2)) &
                        + g(k+1, j, i-1)*g(k+1, j, i) &
                        / MAX(g(k+1, j, i-1) + g(k+1, j, i), MIN(gmol1,gmol2))
                    gn = g(k, j, i)*g(k, j+1, i) &
                        / MAX(g(k, j, i) + g(k, j+1, i), MIN(gmol1,gmol2)) &
                        + g(k+1, j, i)*g(k+1, j+1, i) &
                        / MAX(g(k+1, j, i) + g(k+1, j+1, i), MIN(gmol1,gmol2))
                    gs = g(k, j-1, i)*g(k, j, i) &
                        / MAX(g(k, j-1, i) + g(k, j, i), MIN(gmol1,gmol2)) &
                        + g(k+1, j-1, i)*g(k+1, j, i) &
                        / MAX(g(k+1, j-1, i) + g(k+1, j, i), MIN(gmol1,gmol2))
                    gt = g(k+1, j, i)
                    gb = g(k, j, i)

                    ! Shear stresses
                    tauxze = ge * ( (u(k+1,j,i) - u(k,j,i))     * rdz(k) + (w(k,j,i+1) - w(k,j,i)) * rdx(i)   )
                    tauxzw = gw * ( (u(k+1,j,i-1) - u(k,j,i-1)) * rdz(k) + (w(k,j,i) - w(k,j,i-1)) * rdx(i-1) )
                    tauyzn = gn * ( (v(k+1,j,i) - v(k,j,i))     * rdz(k) + (w(k,j+1,i) - w(k,j,i)) * rdy(j)   )
                    tauyzs = gs * ( (v(k+1,j-1,i) - v(k,j-1,i)) * rdz(k) + (w(k,j,i) - w(k,j-1,i)) * rdy(j-1) )
                    
                    ! Normal stresses
                    tauzzt = gt * 2.0_realk * (w(k+1,j,i) - w(k,j,i)) * rddz(k+1)
                    tauzzb = gb * 2.0_realk * (w(k,j,i) - w(k-1,j,i)) * rddz(k)

                    ! Change due to diffusion
                    wo(k,j,i) = wo(k,j,i) + 1.0_realk/( 0.5_realk * ( d(k,j,i) + d(k,j,i+1) ) ) * &
                        ( ( tauxze - tauxzw ) * rddx(i) + &
                          ( tauyzn - tauyzs ) * rddy(j) + &
                          ( tauzzt - tauzzb ) * rdz(k) )
                END DO
            END DO
        END DO 

    END SUBROUTINE diff_operator

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

    PURE SUBROUTINE comp_advr_centr(kk, jj, ii, q, u, v, w, advrE, advrN, advrT)
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
        INTEGER(intk), INTENT(in) :: kk, jj, ii, q
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(out) :: advrE(kk, jj, ii), advrN(kk, jj, ii), advrT(kk, jj, ii)

        ! Loval variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: k0, j0, i0

        i0 = 0_intk; j0 = 0_intk; k0 = 0_intk

        CALL get_component_specifics(kk, jj, ii, q, i0=i0, j0=j0, k0=k0)        

        DO i = 2, ii-1
            DO j = 2, jj-1
                DO k = 2, kk-1
                    advrE(k,j,i) = 0.5_realk * ( u(k,j,i) + u(k+k0,j+j0,i+i0) )
                    advrN(k,j,i) = 0.5_realk * ( v(k,j,i) + v(k+k0,j+j0,i+i0) )
                    advrT(k,j,i) = 0.5_realk * ( w(k,j,i) + w(k+k0,j+j0,i+i0) )
                END DO
            END DO
        END DO

    END SUBROUTINE comp_advr_centr

    !================================================================

    PURE SUBROUTINE comp_advr_inter(vel1, vel2, vel3, length, scheme, advrD)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !   
    !   Source: 
    !   W. Aniszewski et al., “PArallel, Robust, Interface Simulator
    !   (PARIS),” Computer Physics Communications, vol. 263, 
    !   p. 107849, Jun. 2021, doi: 10.1016/j.cpc.2021.107849.
    !   
    !   PARIS source code function slope_lim (accessed: Mai 2026)
    !----------------------------------------------------------------

        ! Subroutine arguments
        REAL(realk), INTENT(in) :: vel1, vel2, vel3
        REAL(realk), INTENT(in) :: length
        CHARACTER(len=*), INTENT(in) :: scheme
        REAL(realk), INTENT(out) :: advrD 

        ! Loval variables
        REAL(realk) :: limiter
        REAL(realk) :: a, a1, a2

        IF ( scheme == 'ENO' ) THEN
            IF ( abs(vel3-vel2) < abs(vel2-vel1) ) THEN
                limiter = vel3-vel2
            ELSE 
                limiter = vel2-vel1
            END IF
        ELSE IF ( scheme == 'WENO' ) THEN
            a1 = 1.0_realk / ( ( vel2 - vel1 )**2.0_realk + 1.0E-16_realk )
            a2 = 1.0_realk / ( ( vel3 - vel2 )**2.0_realk + 1.0E-16_realk )
            a = a1 + a2
            limiter = ( a1 * ( vel2 - vel1 ) + a2 * (vel3 - vel2) ) / a
        END IF
        
        advrD = vel2 + limiter * length

    END SUBROUTINE comp_advr_inter

    !================================================================

    PURE SUBROUTINE get_component_specifics(kk, jj, ii, q, vel1, vel2, vel3, &
        dx, dy, dz, ddx, ddy, ddz, i0, j0, k0, iStag, jStag, kStag, &
        deltaX, deltaY, deltaZ, velFld)
    !----------------------------------------------------------------
    !   What it does:
    !   Select the specifics for the qth staggered grid.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, q
        REAL(realk), INTENT(in), OPTIONAL :: vel1(kk, jj, ii), vel2(kk, jj, ii), vel3(kk, jj, ii)
        REAL(realk), INTENT(in), OPTIONAL :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in), OPTIONAL :: ddx(ii), ddy(jj), ddz(kk)
        INTEGER(intk), INTENT(out), OPTIONAL :: i0, j0, k0
        REAL(realk), INTENT(out), OPTIONAL :: iStag, jStag, kStag
        REAL(realk), INTENT(out), OPTIONAL :: deltaX(ii), deltaY(jj), deltaZ(kk)
        REAL(realk), INTENT(out), OPTIONAL :: velFld(kk, jj, ii)

        ! Local variables
        ! None

        IF ( q == 1 ) THEN
            IF ( PRESENT(iStag) ) iStag = 1.0_realk
            IF ( PRESENT(jStag) ) jStag = 0.0_realk
            IF ( PRESENT(kStag) ) kStag = 0.0_realk
            IF ( PRESENT(i0) ) i0 = 1_intk
            IF ( PRESENT(j0) ) j0 = 0_intk
            IF ( PRESENT(k0) ) k0 = 0_intk
            IF ( PRESENT(deltaX) ) deltaX = dx
            IF ( PRESENT(deltaY) ) deltaY = ddy
            IF ( PRESENT(deltaZ) ) deltaZ = ddz
            IF ( PRESENT(velFld) ) velFld = vel1
        ELSE IF ( q == 2 ) THEN
            IF ( PRESENT(iStag) ) iStag = 0.0_realk
            IF ( PRESENT(jStag) ) jStag = 1.0_realk
            IF ( PRESENT(kStag) ) kStag = 0.0_realk
            IF ( PRESENT(i0) ) i0 = 0_intk
            IF ( PRESENT(j0) ) j0 = 1_intk
            IF ( PRESENT(k0) ) k0 = 0_intk
            IF ( PRESENT(deltaX) ) deltaX = ddx
            IF ( PRESENT(deltaY) ) deltaY = dy
            IF ( PRESENT(deltaZ) ) deltaZ = ddz
            IF ( PRESENT(velFld) ) velFld = vel2
        ELSE IF ( q == 3 ) THEN
            IF ( PRESENT(iStag) ) iStag = 0.0_realk
            IF ( PRESENT(jStag) ) jStag = 0.0_realk
            IF ( PRESENT(kStag) ) kStag = 1.0_realk
            IF ( PRESENT(i0) ) i0 = 0_intk
            IF ( PRESENT(j0) ) j0 = 0_intk
            IF ( PRESENT(k0) ) k0 = 1_intk
            IF ( PRESENT(deltaX) ) deltaX = ddx
            IF ( PRESENT(deltaY) ) deltaY = ddy
            IF ( PRESENT(deltaZ) ) deltaZ = dz
            IF ( PRESENT(velFld) ) velFld = vel3
        END IF
    
    END SUBROUTINE get_component_specifics

END MODULE multiphase_vof_transport_mod
