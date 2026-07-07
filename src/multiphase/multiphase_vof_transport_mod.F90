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
    USE grids_mod, ONLY: get_mgdims, get_mgbasb, get_gradpxflag
    USE err_mod, ONLY: errr
    USE multiphase_plic_mod, ONLY: comp_frac, iface_recon_wrap, comp_stag_frac_wrap, track_iface, track_iface_vic
    USE rungekutta_mod, ONLY: rk_2n_t
    USE multiphasecore_mod, ONLY: gmol1, gmol2, rho1, rho2, grav, splitting_multiphase, permutation_multiphase
    USE multiphase_material_mod, ONLY: comp_material_property_field
    USE flowcore_mod, ONLY: gradp
    USE connect2_mod, ONLY: connect
    USE parent_mod, ONLY: parent
    USE grids_mod, ONLY: minlevel, maxlevel
    USE multiphase_io_mod, ONLY: apprVol
    USE err_mod, ONLY: errr
    
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
            DO i = 1, ii-1
                DO j = 1, jj-1
                    DO k = 1, kk-1
                        IF ( u(k,j,i) > tol ) THEN
                            IF ( isInterface(k,j,i) ) THEN
                                ! CASE 1
                                !       W   <1     E       
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

                                CALL comp_frac(fluxedProp, fluxAlpha, field(k,j,i), &
                                    fluxWidth, ddy(j), ddz(k), normx(k,j,i), normy(k,j,i), normz(k,j,i), tol)
                                
                                flux = fluxedProp  * ( abs( u(k,j,i) ) * dt / ddx(i) )
                            ELSE
                                ! CASE 2
                                !       W   =1     E
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
                                !       W          E     <1
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

                                CALL comp_frac(fluxedProp, fluxAlpha, field(k,j,i+1), fluxWidth, ddy(j), ddz(k), normx(k,j,i+1), normy(k,j,i+1), normz(k,j,i+1), tol)

                                flux = fluxedProp * ( abs( u(k,j,i) ) * dt / ddx(i+1) )
                            ELSE
                                ! CASE 4
                                !       W          E    =1 
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
            DO i = 1, ii-1
                DO j = 1, jj-1
                    DO k = 1, kk-1
                        IF ( v(k,j,i) > tol ) THEN
                            IF ( isInterface(k,j,i) ) THEN
                                ! See explanation above.
                                fluxWidth = abs( v(k,j,i) ) * dt
                                fluxAlpha = alpha(k,j,i) - normy(k,j,i) * ( ddy(j) - fluxWidth )

                                CALL comp_frac(fluxedProp, fluxAlpha, field(k,j,i), ddx(i), fluxWidth, ddz(k), normx(k,j,i), normy(k,j,i), normz(k,j,i), tol)

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

                                CALL comp_frac(fluxedProp, fluxAlpha, field(k,j+1,i), ddx(i), fluxWidth, ddz(k), normx(k,j+1,i), normy(k,j+1,i), normz(k,j+1,i), tol)

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
            DO i = 1, ii-1
                DO j = 1, jj-1
                    DO k = 1, kk-1
                        IF ( w(k,j,i) > tol ) THEN
                            IF ( isInterface(k,j,i) ) THEN
                                ! See explanation above.
                                fluxWidth = abs( w(k,j,i) ) * dt
                                fluxAlpha = alpha(k,j,i) - normz(k,j,i) * ( ddz(k) - fluxWidth )

                                CALL comp_frac(fluxedProp, fluxAlpha, field(k,j,i), ddx(i), ddy(j), fluxWidth, normx(k,j,i), normy(k,j,i), normz(k,j,i), tol)

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

                                CALL comp_frac(fluxedProp, fluxAlpha, field(k+1,j,i), ddx(i), ddy(j), fluxWidth, normx(k+1,j,i), normy(k+1,j,i), normz(k+1,j,i), tol)

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
        REAL(realk) :: dsx(ii), dsy(jj), dsz(kk)
        INTEGER(intk) :: k, j, i
        REAL(realk) :: flux, complementFlux, fluxedProp, fluxWidth, fluxAlpha
        REAL(realk) :: advr(kk, jj, ii)

        CALL get_spatial_extents(kk, jj, ii, q, dx, dy, dz, ddx, ddy, ddz, dsx, dsy, dsz)
        CALL comp_advr_linear_interpolation(kk, jj, ii, splitDir, u, v, w, advr)

        IF ( splitDir == 1 ) THEN
            DO i = 1, ii-1
                DO j = 1, jj-1
                    DO k = 1, kk-1
                        IF ( advr(k,j,i) > tol ) THEN
                            IF ( isInterface(k,j,i) ) THEN
                                ! See explanation above.
                                fluxWidth = abs( advr(k,j,i) ) * dt
                                fluxAlpha = alpha(k,j,i) - normx(k,j,i) * ( dsx(i)/2.0_realk - fluxWidth )

                                CALL comp_frac(fluxedProp, fluxAlpha, field(k,j,i), fluxWidth, ddy(j), ddz(k), normx(k,j,i), normy(k,j,i), normz(k,j,i), tol)
                                
                                flux = fluxedProp  * ( abs( advr(k,j,i) ) * dt / dsx(i) )
                                complementFlux = ( 1.0_realk - fluxedProp ) * ( abs( advr(k,j,i) ) * dt / dsx(i) )
                            ELSE
                                ! See explanation above.
                                flux = field(k,j,i) * ( abs( advr(k,j,i) ) * dt / dsx(i) )
                                complementFlux = max( 1.0_realk - field(k,j,i), 0.0_realk ) * ( abs( advr(k,j,i) ) * dt / dsx(i) )
                            END IF
                        ELSE IF ( advr(k,j,i) < -tol ) THEN
                            IF ( isInterface(k,j,i+1) ) THEN
                                ! See explanation above.
                                fluxWidth = abs( advr(k,j,i) ) * dt
                                fluxAlpha = alpha(k,j,i+1) - normx(k,j,i+1) * dsx(i+1)/2.0_realk

                                CALL comp_frac(fluxedProp, fluxAlpha, field(k,j,i+1), fluxWidth, ddy(j), ddz(k), normx(k,j,i+1), normy(k,j,i+1), normz(k,j,i+1), tol)

                                flux = fluxedProp * ( abs( advr(k,j,i) ) * dt / dsx(i+1) )
                                complementFlux = ( 1.0_realk - fluxedProp ) * ( abs( advr(k,j,i) ) * dt / dsx(i+1) )
                            ELSE
                                ! See explanation above.
                                flux = field(k,j,i+1) * ( abs( advr(k,j,i) ) * dt / dsx(i+1) )
                                complementFlux = max( 1.0_realk - field(k,j,i+1), 0.0_realk ) * ( abs( advr(k,j,i) ) * dt / dsx(i+1) )
                            END IF
                        ELSE
                            ! See explanation above.
                            flux = 0.0_realk
                            complementFlux = 0.0_realk
                        END IF
                        fieldFlux(k,j,i) = sign( 1.0_realk, advr(k,j,i) ) * flux / dt
                        complementFieldFlux(k,j,i) = sign( 1.0_realk, advr(k,j,i) ) * complementFlux / dt
                    END DO
                END DO
            END DO
        ELSE IF ( splitDir == 2 ) THEN
            DO i = 1, ii-1
                DO j = 1, jj-1
                    DO k = 1, kk-1
                        IF ( advr(k,j,i) > tol ) THEN
                            IF ( isInterface(k,j,i) ) THEN
                                ! See explanation above.
                                fluxWidth = abs( advr(k,j,i) ) * dt
                                fluxAlpha = alpha(k,j,i) - normy(k,j,i) * ( dsy(j)/2.0_realk - fluxWidth )

                                CALL comp_frac(fluxedProp, fluxAlpha, field(k,j,i), ddx(i), fluxWidth, ddz(k), normx(k,j,i), normy(k,j,i), normz(k,j,i), tol)

                                flux = fluxedProp * ( abs( advr(k,j,i) ) * dt / dsy(j) )
                                complementFlux = ( 1.0_realk - fluxedProp ) * ( abs( advr(k,j,i) ) * dt / dsy(j) )
                            ELSE
                                ! See explanation above.
                                flux = field(k,j,i) * ( abs( advr(k,j,i) ) * dt / dsy(j) )
                                complementFlux = max( 1.0_realk - field(k,j,i), 0.0_realk ) * ( abs( advr(k,j,i) ) * dt / dsy(j) )
                            END IF
                        ELSE IF ( advr(k,j,i) < -tol ) THEN
                            IF ( isInterface(k,j+1,i) ) THEN
                                ! See explanation above.
                                fluxWidth = abs( advr(k,j,i) ) * dt
                                fluxAlpha = alpha(k,j+1,i) - normy(k,j+1,i) * dsy(j+1)/2.0_realk

                                CALL comp_frac(fluxedProp, fluxAlpha, field(k,j+1,i), ddx(i), fluxWidth, ddz(k), normx(k,j+1,i), normy(k,j+1,i), normz(k,j+1,i), tol)

                                flux = fluxedProp * ( abs( advr(k,j,i) ) * dt / dsy(j+1) )
                                complementFlux = ( 1.0_realk - fluxedProp ) * ( abs( advr(k,j,i) ) * dt / dsy(j+1) )
                            ELSE
                                ! See explanation above.
                                flux = field(k,j+1,i) * ( abs( advr(k,j,i) ) * dt / dsy(j+1) )
                                complementFlux = max( 1.0_realk - field(k,j+1,i), 0.0_realk ) * ( abs( advr(k,j,i) ) * dt / dsy(j+1) )
                            END IF
                        ELSE
                            ! See explanation above.
                            flux = 0.0_realk
                            complementFlux = 0.0_realk
                        END IF
                        fieldFlux(k,j,i) = sign( 1.0_realk, advr(k,j,i) ) * flux / dt
                        complementFieldFlux(k,j,i) = sign( 1.0_realk, advr(k,j,i) ) * complementFlux / dt
                    END DO
                END DO
            END DO
        ELSE IF ( splitDir == 3 ) THEN
            DO i = 1, ii-1
                DO j = 1, jj-1
                    DO k = 1, kk-1
                        IF ( advr(k,j,i) > tol ) THEN
                            IF ( isInterface(k,j,i) ) THEN
                                ! See explanation above.
                                fluxWidth = abs( advr(k,j,i) ) * dt
                                fluxAlpha = alpha(k,j,i) - normz(k,j,i) * ( dsz(k)/2.0_realk - fluxWidth )

                                CALL comp_frac(fluxedProp, fluxAlpha, field(k,j,i), ddx(i), ddy(j), fluxWidth, normx(k,j,i), normy(k,j,i), normz(k,j,i), tol)

                                flux = fluxedProp * ( abs( advr(k,j,i) ) * dt / dsz(k) )
                                complementFlux = ( 1.0_realk - fluxedProp ) * ( abs( advr(k,j,i) ) * dt / dsz(k) )
                            ELSE
                                ! See explanation above.
                                flux = field(k,j,i) * ( abs( advr(k,j,i) ) * dt / dsz(k) )
                                complementFlux = max( 1.0_realk - field(k,j,i), 0.0_realk ) * ( abs( advr(k,j,i) ) * dt / dsz(k) )
                            END IF
                        ELSE IF ( advr(k,j,i) < -tol ) THEN
                            IF ( isInterface(k+1,j,i) ) THEN
                                ! See explanation above.
                                fluxWidth = abs( advr(k,j,i) ) * dt
                                fluxAlpha = alpha(k+1,j,i) - normz(k+1,j,i) * dsz(k+1)/2.0_realk

                                CALL comp_frac(fluxedProp, fluxAlpha, field(k+1,j,i), ddx(i), ddy(j), fluxWidth, normx(k+1,j,i), normy(k+1,j,i), normz(k+1,j,i), tol)

                                flux = fluxedProp * ( abs( advr(k,j,i) ) * dt / dsz(k+1) )
                                complementFlux = ( 1.0_realk - fluxedProp ) * ( abs( advr(k,j,i) ) * dt / dsz(k+1) )
                            ELSE
                                ! See explanation above.
                                flux = field(k+1,j,i) * ( abs( advr(k,j,i) ) * dt / dsz(k+1) )
                                complementFlux = max( 1.0_realk - field(k+1,j,i), 0.0_realk ) * ( abs( advr(k,j,i) ) * dt / dsz(k+1) )
                            END IF
                        ELSE
                            ! See explanation above.
                            flux = 0.0_realk
                            complementFLux = 0.0_realk
                        END IF
                        fieldFlux(k,j,i) = sign( 1.0_realk, advr(k,j,i) ) * flux / dt
                        complementFieldFlux(k,j,i) = sign( 1.0_realk, advr(k,j,i) ) * complementFlux / dt
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

    SUBROUTINE multiphase_solve(u_f, v_f, w_f, vff_f, p_f, g_f, d_f, dt, itstep, uo_f, vo_f, wo_f)
    !----------------------------------------------------------------
    !   What it does:
    !    
    !----------------------------------------------------------------

        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: u_f
        TYPE(field_t), INTENT(inout) :: v_f
        TYPE(field_t), INTENT(inout) :: w_f
        TYPE(field_t), INTENT(in) :: vff_f
        TYPE(field_t), INTENT(in) :: p_f
        TYPE(field_t), INTENT(in) :: g_f
        TYPE(field_t), INTENT(in) :: d_f
        REAL(realk), INTENT(in) :: dt
        INTEGER(intk), INTENT(in) :: itstep
        TYPE(field_t), INTENT(inout) :: uo_f, vo_f, wo_f

        ! Local variables
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: u, v, w
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: vff, p, g, d
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: uo, vo, wo

        TYPE(field_t), POINTER :: dx_f, dy_f, dz_f, ddx_f, ddy_f, ddz_f
        TYPE(field_t), POINTER :: rdx_f, rdy_f, rdz_f, rddx_f, rddy_f, rddz_f
        TYPE(field_t), POINTER :: normx_f, normy_f, normz_f
        TYPE(field_t), POINTER :: alpha_f

        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:), ddx(:), ddy(:), ddz(:)
        REAL(realk), POINTER, CONTIGUOUS :: rdx(:), rdy(:), rdz(:), rddx(:), rddy(:), rddz(:)
        REAL(realk), POINTER, CONTIGUOUS :: normx(:,:,:), normy(:,:,:), normz(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: alpha(:,:,:)

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

        CALL get_field(normx_f, "NORMX"); CALL get_field(normy_f, "NORMY"); CALL get_field(normz_f, "NORMZ")
         CALL get_field(alpha_f, "ALPHA")

        DO i = 1, nmygrids
            igrid = mygrids(i)

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_mgbasb(nfro, nbac, nrgt, nlft, nbot, ntop, igrid)

            CALL u_f%get_ptr(u, igrid)
            CALL v_f%get_ptr(v, igrid)
            CALL w_f%get_ptr(w, igrid)
            CALL vff_f%get_ptr(vff, igrid)
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

            CALL normx_f%get_ptr(normx, igrid); CALL normy_f%get_ptr(normy, igrid); CALL normz_f%get_ptr(normz, igrid)
            CALL alpha_f%get_ptr(alpha, igrid)

            IF ( .NOT. ALLOCATED(vffPrev) ) ALLOCATE(vffPrev(kk, jj, ii))
            vffPrev = vff

            CALL adve_operator(kk, jj, ii, u, v, w, vff, dx, dy, dz, ddx, ddy, ddz, dt, itstep, tol, &
                normx, normy, normz, alpha, uo, vo, wo)

            CALL diff_operator(kk, jj, ii, u, v, w, vffPrev, g, rdx, rdy, rdz, rddx, rddy, rddz, uo, vo, wo)

            CALL pres_operator(kk, jj, ii, vff, p, rdx, rdy, rdz, igrid, uo, vo, wo)

            CALL exte_operator(kk, jj, ii, vff, dx, dy, dz, ddx, ddy, ddz, uo, vo, wo)

            u = u + uo * dt
            v = v + vo * dt
            w = w + wo * dt

        END DO

        CALL check_continuity(tol)
        CALL check_solenoidality(tol)

    END SUBROUTINE multiphase_solve

    !================================================================

    SUBROUTINE adve_operator(kk, jj, ii, u, v, w, vff, dx, dy, dz, ddx, ddy, ddz, dt, itstep, tol, &
        normx, normy, normz, alpha, uo, vo, wo)
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
        REAL(realk), INTENT(in) :: tol, dt
        INTEGER(intk), INTENT(in) :: itstep
        REAL(realk), INTENT(out) :: normx(kk, jj, ii), normy(kk, jj, ii), normz(kk, jj, ii), alpha(kk, jj, ii)
        REAL(realk), INTENT(inout) :: uo(kk, jj, ii), vo(kk, jj, ii), wo(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: q, advSeq(3), l, splitDir
        REAL(realk), ALLOCATABLE :: vffStag(:,:,:)
        REAL(realk), ALLOCATABLE :: dStag(:,:,:)
        LOGICAL, ALLOCATABLE :: isInterface(:,:,:)
        LOGICAL, ALLOCATABLE :: isNearInterface(:,:,:)
        REAL(realk), ALLOCATABLE :: vffFlux(:,:,:), complVffFlux(:,:,:), mom(:,:,:), cWY(:,:,:), cWYStag(:,:,:)
        REAL(realk), ALLOCATABLE :: advr(:,:,:)
        REAL(realk), ALLOCATABLE :: adve(:,:,:)
        REAL(realk), ALLOCATABLE :: vel(:,:,:), velo(:,:,:,:)
        REAL(realk), ALLOCATABLE :: mom4D(:,:,:,:), vffStag4D(:,:,:,:), cWYStag4D(:,:,:,:)

        IF (.NOT. ALLOCATED(vffStag))             ALLOCATE(vffStag(kk,jj,ii))
        IF (.NOT. ALLOCATED(dStag))               ALLOCATE(dStag(kk,jj,ii))
        IF (.NOT. ALLOCATED(isInterface))         ALLOCATE(isInterface(kk,jj,ii))
        IF (.NOT. ALLOCATED(isNearInterface))     ALLOCATE(isNearInterface(kk,jj,ii))
        IF (.NOT. ALLOCATED(vffFlux))             ALLOCATE(vffFlux(kk,jj,ii))
        IF (.NOT. ALLOCATED(complVffFlux))        ALLOCATE(complVffFlux(kk,jj,ii))
        IF (.NOT. ALLOCATED(mom))                 ALLOCATE(mom(kk,jj,ii))
        IF (.NOT. ALLOCATED(cWY))                 ALLOCATE(cWY(kk,jj,ii))
        IF (.NOT. ALLOCATED(cWYStag))             ALLOCATE(cWYStag(kk,jj,ii))
        IF (.NOT. ALLOCATED(advr))                ALLOCATE(advr(kk,jj,ii))
        IF (.NOT. ALLOCATED(adve))                ALLOCATE(adve(kk,jj,ii))
        IF (.NOT. ALLOCATED(vel))                 ALLOCATE(vel(kk,jj,ii))
        IF (.NOt. ALLOCATED(mom4D))               ALLOCATE(mom4D(kk,jj,ii,3))
        IF (.NOt. ALLOCATED(vffStag4D))           ALLOCATE(vffStag4D(kk,jj,ii,3))
        IF (.NOT. ALLOCATED(cWYStag4D))           ALLOCATE(cWYStag4D(kk,jj,ii,3))
        IF (.NOT. ALLOCATED(velo))                ALLOCATE(velo(kk,jj,ii,3))
        
        IF ( splitting_multiphase == "component-wise" ) THEN

            CALL def_advection_sequence(itstep, advSeq)
            CALL iface_recon_wrap(kk, jj, ii, vff, ddx, ddy, ddz, tol, normx, normy, normz, alpha, isInterface, isNearInterface)

            DO q = 1, 3
                
                CALL comp_stag_frac_wrap(kk, jj, ii, q, alpha, vff, isInterface, ddx, ddy, ddz, normx, normy, normz, tol, vffStag)
                CALL comp_material_property_field(kk, jj, ii, vffStag, rho1, rho2, dStag)
                CALL comp_momentum(kk, jj, ii, q, dStag, u, v, w, mom)
                CALL comp_cWY(kk, jj, ii, vffStag, cWYStag)

                DO splitDir = 1, 3

                    l = advSeq(splitDir)

                    CALL comp_advr_linear_interpolation(kk, jj, ii, l, u, v, w, advr)
                    CALL comp_adve_quick(kk, jj, ii, q, l, u, v, w, advr, vffStag, tol, adve)
                    CALL comp_flux_stag(kk, jj, ii, q, l, vff, isInterface, u, v, w, alpha, dt, normx, normy, normz, dx, dy, dz, ddx, ddy, ddz, tol, vffFlux, complVffFlux)
                    CALL adv_mom(kk, jj, ii, q, l, u, v, w, advr, adve, vffStag, vffFlux, complVffFlux, cWYStag, dx, dy, dz, ddx, ddy, ddz, dt, mom)
                    CALL adv_vof(kk, jj, ii, l, vffFlux, cWYStag, advr, dx, dy, dz, ddx, ddy, ddz, dt, tol, vffStag)

                END DO

                CALL comp_velocity_change(kk, jj, ii, q, u, v, w, vffStag, mom, dt, velo)

            END DO

            uo = velo(:,:,:,1)
            vo = velo(:,:,:,2)
            wo = velo(:,:,:,3)

            CALL comp_cWY(kk, jj, ii, vff, cWY)

            DO splitDir = 1, 3

                l = advSeq(splitDir)

                CALL iface_recon_wrap(kk, jj, ii, vff, ddx, ddy, ddz, tol, normx, normy, normz, alpha, isInterface, isNearInterface)
                CALL comp_flux_cent(kk, jj, ii, l, vff, isInterface, u, v, w, alpha, dt, normx, normy, normz, ddx, ddy, ddz, tol, vffFlux)
                CALL get_condit_velocity(kk, jj, ii, l, u, v, w, vel)
                CALL adv_vof(kk, jj, ii, l, vffFlux, cWY, vel, dx, dy, dz, ddx, ddy, ddz, dt, tol, vff)
                CALL clip_vff(kk, jj, ii, tol, vff)
                CALL app_bcon()

            END DO

        ELSE IF ( splitting_multiphase == "direction-wise" ) THEN

            CALL def_advection_sequence(itstep, advSeq)
            CALL iface_recon_wrap(kk, jj, ii, vff, ddx, ddy, ddz, tol, normx, normy, normz, alpha, isInterface, isNearInterface)

            DO q = 1, 3

                CALL comp_stag_frac_wrap(kk, jj, ii, q, alpha, vff, isInterface, ddx, ddy, ddz, normx, normy, normz, tol, vffStag)
                CALL comp_material_property_field(kk, jj, ii, vffStag, rho1, rho2, dStag)
                CALL comp_momentum(kk, jj, ii, q, dStag, u, v, w, mom)
                CALL comp_cWY(kk, jj, ii, vffStag, cWYStag)

                vffStag4D(:,:,:,q) = vffStag
                mom4D(:,:,:,q) = mom
                cWYStag4D(:,:,:,q) = cWYStag

            END DO

            CALL comp_cWY(kk, jj, ii, vff, cWY)

            DO splitDir = 1, 3

                l = advSeq(splitDir)

                DO q = 1, 3

                    vffStag = vffStag4D(:,:,:,q)
                    mom = mom4D(:,:,:,q)
                    cWYStag = cWYStag4D(:,:,:,q)

                    CALL comp_advr_linear_interpolation(kk, jj, ii, l, u, v, w, advr)
                    CALL comp_adve_quick(kk, jj, ii, q, l, u, v, w, advr, vffStag, tol, adve)
                    CALL comp_flux_stag(kk, jj, ii, q, l, vff, isInterface, u, v, w, alpha, dt, normx, normy, normz, dx, dy, dz, ddx, ddy, ddz, tol, vffFlux, complVffFlux)
                    CALL adv_mom(kk, jj, ii, q, l, u, v, w, advr, adve, vffStag, vffFlux, complVffFlux, cWYStag, dx, dy, dz, ddx, ddy, ddz, dt, mom)
                    CALL adv_vof(kk, jj, ii, l, vffFlux, cWYStag, advr, dx, dy, dz, ddx, ddy, ddz, dt, tol, vffStag)
                    CALL comp_velocity_change(kk, jj, ii, q, u, v, w, vffStag, mom, dt, velo)

                END DO

                uo = velo(:,:,:,1)
                vo = velo(:,:,:,2)
                wo = velo(:,:,:,3)

                CALL iface_recon_wrap(kk, jj, ii, vff, ddx, ddy, ddz, tol, normx, normy, normz, alpha, isInterface, isNearInterface)
                CALL comp_flux_cent(kk, jj, ii, l, vff, isInterface, u, v, w, alpha, dt, normx, normy, normz, ddx, ddy, ddz, tol, vffFlux)
                CALL get_condit_velocity(kk, jj, ii, l, u, v, w, vel)
                CALL adv_vof(kk, jj, ii, l, vffFlux, cWY, vel, dx, dy, dz, ddx, ddy, ddz, dt, tol, vff)
                CALL clip_vff(kk, jj, ii, tol, vff)
                CALL app_bcon()

            END DO

        END IF

    END SUBROUTINE adve_operator

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

        CALL comp_material_property_field(kk, jj, ii, vff, rho1, rho2, d)

        DO i = 2, ii-1
            DO j = 2, jj-1
                DO k = 2, kk-1
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
                    uo(k,j,i) = uo(k,j,i) + 2.0_realk/( d(k,j,i) + d(k,j,i+1) ) * &
                        ( ( tauxxe - tauxxw ) * rdx(i) + &
                          ( tauyxn - tauyxs ) * rddy(j) + &
                          ( tauzxt - tauzxb ) * rddz(k) )
                END DO
            END DO
        END DO

        DO i = 2, ii-1
            DO j = 2, jj-1
                DO k = 2, kk-1
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
                    vo(k,j,i) = vo(k,j,i) + 2.0_realk/( d(k,j,i) + d(k,j+1,i) ) * &
                        ( ( tauxye - tauxyw ) * rddx(i) + &
                          ( tauyyn - tauyys ) * rdy(j) + &
                          ( tauzyt - tauzyb ) * rddz(k) )
                END DO
            END DO
        END DO

        DO i = 2, ii-1
            DO j = 2, jj-1
                DO k = 2, kk-1
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
                    wo(k,j,i) = wo(k,j,i) + 2.0_realk/( d(k,j,i) + d(k,j,i+1) ) * &
                        ( ( tauxze - tauxzw ) * rddx(i) + &
                          ( tauyzn - tauyzs ) * rddy(j) + &
                          ( tauzzt - tauzzb ) * rdz(k) )
                END DO
            END DO
        END DO 

    END SUBROUTINE diff_operator

    !================================================================

    SUBROUTINE pres_operator(kk, jj, ii, vff, p, rdx, rdy, rdz, igrid, uo, vo, wo)
    !----------------------------------------------------------------
    !   What it does:
    !    
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: vff(kk, jj, ii), p(kk, jj, ii)
        REAL(realk), INTENT(in) :: rdx(ii), rdy(jj), rdz(kk)
        INTEGER(intk), INTENT(in) :: igrid
        REAL(realk), INTENT(inout) :: uo(kk, jj, ii), vo(kk, jj, ii), wo(kk, jj, ii)

        ! Local variables
        REAL(realk) :: d(kk, jj, ii)
        INTEGER(intk) :: gradpflag
        REAL(realk) :: gpx, gpy, gpz
        INTEGER(intk) :: i, j, k

        CALL comp_material_property_field(kk, jj, ii, vff, rho1, rho2, d)

        CALL get_gradpxflag(gradpflag, igrid)
        gpx = gradp(1)*gradpflag
        gpy = gradp(2)*gradpflag
        gpz = gradp(3)*gradpflag

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    uo(k,j,i) = uo(k,j,i) - 2.0_realk / ( d(k,j,i) + d(k,j,i+1) ) * ( p(k,j,i+1) - p(k,j,i) + gpx / rdx(i) ) * rdx(i)
                END DO
            END DO
        END DO

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    vo(k,j,i) = vo(k,j,i) - 2.0_realk / ( d(k,j,i) + d(k,j+1,i) ) * ( p(k,j+1,i) - p(k,j,i) + gpy / rdy(j) ) * rdy(j)
                END DO
            END DO
        END DO

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    wo(k,j,i) = wo(k,j,i) - 2.0_realk / ( d(k,j,i) + d(k+1,j,i) ) * ( p(k+1,j,i) - p(k,j,i) + gpz / rdz(k) ) * rdz(k)
                END DO
            END DO
        END DO

    END SUBROUTINE pres_operator

    !================================================================

    SUBROUTINE exte_operator(kk, jj, ii, vff, dx, dy, dz, ddx, ddy, ddz, uo, vo, wo)
    !----------------------------------------------------------------
    !   What it does:
    !    
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: vff(kk, jj, ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(inout) :: uo(kk, jj, ii), vo(kk, jj, ii), wo(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: i, j, k
        REAL(realk) :: d(kk, jj, ii)
        return
        CALL comp_material_property_field(kk, jj, ii, vff, rho1, rho2, d)

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    uo(i,j,i) = uo(k,j,i) - grav(1) ! * 0.5 * ( d(k,j,i) + d(k,j,i+1) )
                    vo(i,j,i) = vo(k,j,i) - grav(2) ! * 0.5 * ( d(k,j,i) + d(k,j+1,i) )
                    wo(i,j,i) = wo(k,j,i) - grav(3) ! * 0.5 * ( d(k,j,i) + d(k+1,j,i) )
                ENDDO
            ENDDO
        ENDDO

    END SUBROUTINE exte_operator

    !================================================================

    SUBROUTINE adv_vof(kk, jj, ii, l, vffFlux, cWY, vel, dx, dy, dz, ddx, ddy, ddz, dt, tol, vff)
    !----------------------------------------------------------------
    !   What it does:
    !    
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: l
        REAL(realk), INTENT(in) :: vffFlux(kk, jj, ii)
        REAL(realk), INTENT(in) :: cWY(kk, jj, ii), vel(kk,jj,ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: dt, tol
        REAL(realk), INTENT(inout) :: vff(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: i, j, k
        INTEGER(intk) :: il, jl, kl
        REAL(realk) :: dsx(ii), dsy(jj), dsz(kk)
        REAL(realk) :: div(kk, jj, ii)

        CALL get_spatial_indices(kk, jj, ii, l, il, jl, kl)
        CALL get_spatial_extents(kk, jj, ii, l, dx, dy, dz, ddx, ddy, ddz, dsx, dsy, dsz)

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    div(k,j,i) = ( vel(k,j,i) - vel(k-kl,j-jl,i-il) ) / ( il * dsx(i) + jl * dsy(j) + kl * dsz(k) )
                    vff(k,j,i) = vff(k,j,i) - dt * ( vffFlux(k,j,i) - vffFlux(k-kl,j-jl,i-il) ) + dt * cWY(k,j,i) * div(k,j,i)
                END DO
            END DO
        END DO

    END SUBROUTINE adv_vof

    !================================================================

    SUBROUTINE adv_mom(kk, jj, ii, q, l, u, v, w, advr, adve, vff, vffFlux, complVffFlux, cWY, dx, dy, dz, ddx, ddy, ddz, dt, mom)
    !----------------------------------------------------------------
    !   What it does:
    !    
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, q, l
        REAL(realk), INTENT(in) :: u(kk,jj,ii), v(kk,jj,ii), w(kk,jj,ii)
        REAL(realk), INTENT(in) :: advr(kk,jj,ii), adve(kk,jj,ii)
        REAL(realk), INTENT(in) :: vff(kk, jj, ii), vffFlux(kk, jj, ii), complVffFlux(kk, jj, ii)
        REAL(realk), INTENT(in) :: cWY(kk, jj, ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: dt
        REAL(realk), INTENT(inout) :: mom(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: i, j, k
        INTEGER(intk) :: il, jl, kl
        REAL(realk) :: dsx(ii), dsy(jj), dsz(kk)
        REAL(realk) :: vel(kk,jj,ii)
        REAL(realk) :: dStag(kk, jj, ii)
        REAL(realk) :: momFlux(kk, jj, ii)
        REAL(realk) :: div, com

        CALL get_spatial_indices(kk, jj, ii, l, il, jl, kl)
        CALL get_spatial_extents(kk, jj, ii, l, dx, dy, dz, ddx, ddy, ddz, dsx, dsy, dsz)
        CALL get_condit_velocity(kk, jj, ii, q, u, v, w, vel)

        CALL comp_material_property_field(kk, jj, ii, vff, rho1, rho2, dStag)

        DO i = 2, ii-1
            DO j = 2, jj-1
                DO k = 2, kk-1
                    momFlux(k,j,i) = adve(k,j,i) * ( rho1 * vffFlux(k,j,i) + rho2 * complVffFlux(k,j,i) )
                END DO
            END DO 
        END DO

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    div = ( advr(k,j,i) - advr(k-kl,j-jl,i-il) ) / ( il * dsx(i) + jl * dsy(j) + kl * dsz(k) )
                    com = ( rho1 * cWY(k,j,i) + rho2 * (1.0_realk - cWY(k,j,i)) ) * div
                    mom(k,j,i) = mom(k,j,i) - dt * ( momFlux(k,j,i) - momFlux(k-kl,j-jl,i-il) ) + dt * vel(k,j,i) * com
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
        REAL(realk) :: vel(kk, jj, ii)

        CALL get_condit_velocity(kk, jj, ii, q, u, v, w, vel)

        DO i = 1, ii
            DO j = 1, jj
                DO k = 1, kk
                    mom(k,j,i) = vel(k,j,i) * dStag(k,j,i)
                END DO
            END DO
        END DO

    END SUBROUTINE comp_momentum

    !================================================================

    SUBROUTINE comp_velocity_change(kk, jj, ii, q, u, v, w, vff, mom, dt, velo)
    !----------------------------------------------------------------
    !   What it does:
    !    
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: q
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii), vff(kk, jj, ii), mom(kk, jj, ii)
        REAL(realk), INTENT(in) :: dt
        REAL(realk), INTENT(out) :: velo(kk,jj,ii,q)

        ! Local variables
        INTEGER(intk) :: i, j, k
        REAL(realk) :: vel(kk,jj,ii)
        REAL(realk) :: dStag(kk, jj, ii)

        CALL get_condit_velocity(kk, jj, ii, q, u, v, w, vel)
        CALL comp_material_property_field(kk, jj, ii, vff, rho1, rho2, dStag)
    
        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    IF ( vff(k,j,i) >= 0.0_realk .AND. vff(k,j,i) <= 1.0_realk ) THEN
                        velo(k,j,i,q) = ( mom(k,j,i) / dStag(k,j,i) - vel(k,j,i) ) / dt
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

    SUBROUTINE check_solenoidality(tol)
    !----------------------------------------------------------------
    !   What it does:
    !   Checks, if the field is solenoidal. Only used during coding.
    !   Solenoidality is assured by pressure correction.
    !----------------------------------------------------------------

        ! Subroutine arguments
        REAL(realk) :: tol

        ! Local variables
        TYPE(field_t), POINTER :: u_f, v_f, w_f
        TYPE(field_t), POINTER :: ddx_f, ddy_f, ddz_f
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        REAL(realk), POINTER, CONTIGUOUS :: u(:,:,:), v(:,:,:), w(:,:,:)
        INTEGER(intk) :: kk, jj, ii, k, j, i, n, igrid
        REAL(realk), ALLOCATABLE :: div(:,:,:)
        REAL(realk) :: uChar, lChar, L1tol, L2tol, Linftol, L1Norm, L2Norm, LinfNorm

        return

        CALL get_field(u_f, "U")
        CALL get_field(v_f, "V")
        CALL get_field(w_f, "W")
        CALL get_field(ddx_f, "DDX")
        CALL get_field(ddy_f, "DDY")
        CALL get_field(ddz_f, "DDZ")

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL u_f%get_ptr(u, igrid)
            CALL v_f%get_ptr(v, igrid)
            CALL w_f%get_ptr(w, igrid)
            CALL ddx_f%get_ptr(ddx, igrid)
            CALL ddy_f%get_ptr(ddy, igrid)
            CALL ddz_f%get_ptr(ddz, igrid)

            IF ( .NOT. ALLOCATED(div) ) ALLOCATE(div(kk,jj,ii))

            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        div(k,j,i) = ( u(k,j,i+1) - u(k,j,i) )/ddx(i) + &
                                     ( v(k,j+1,i) - v(k,j,i) )/ddy(j) + &
                                     ( w(k+1,j,i) - w(k,j,i) )/ddz(k)
                    END DO
                END DO
            END DO

            uChar = MAX(MAXVAL(ABS(u(3:kk-2,3:jj-2,3:ii-2))), &
                        MAXVAL(ABS(v(3:kk-2,3:jj-2,3:ii-2))), &
                        MAXVAL(ABS(w(3:kk-2,3:jj-2,3:ii-2))))
            lChar = MAX(MAXVAL(ABS(ddx(3:ii-2))), &
                        MAXVAL(ABS(ddy(3:jj-2))), &
                        MAXVAL(ABS(ddz(3:kk-2))))

            L1tol = tol * uChar / lChar
            L2tol = tol * uChar / lChar * ( SIZE(div(3:kk-2,3:jj-2,3:ii-2)) )**(1.0_realk/2.0_realk)
            Linftol = tol * uChar / lChar * SIZE(div(3:kk-2,3:jj-2,3:ii-2))
            
            L1Norm = SUM( ABS(div(3:kk-2,3:jj-2,3:ii-2)) )
            L2Norm = ( SUM( div(3:kk-2,3:jj-2,3:ii-2)**2.0_realk ) )**(1.0_realk/2.0_realk)
            LinfNorm = MAXVAL( ABS(div(3:kk-2,3:jj-2,3:ii-2)) )

            IF ( L1Norm > L1tol .OR. L2Norm > L2tol .OR. LinfNorm > Linftol ) THEN
                WRITE(*,*) "Velocity not solenoidal: L1 =", L1Norm, " L2 =", L2Norm, " Linf =", LinfNorm
            END IF

        ENDDO

    END SUBROUTINE check_solenoidality

    !================================================================

    SUBROUTINE check_continuity(tol)
    !----------------------------------------------------------------
    !   What it does:
    !    
    !----------------------------------------------------------------

        ! Subroutine arguments
        REAL(realk) :: tol

        ! Local variables
        TYPE(field_t), POINTER :: vff_f
        TYPE(field_t), POINTER :: ddx_f, ddy_f, ddz_f
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        REAL(realk), POINTER, CONTIGUOUS :: vff(:,:,:)
        INTEGER(intk) :: kk, jj, ii, k, j, i, n, igrid
        REAL(realk) :: volFl1

        CALL get_field(vff_f, "VFF")
        CALL get_field(ddx_f, "DDX")
        CALL get_field(ddy_f, "DDY")
        CALL get_field(ddz_f, "DDZ")

        volFl1 = 0.0_realk

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL vff_f%get_ptr(vff, igrid)
            CALL ddx_f%get_ptr(ddx, igrid)
            CALL ddy_f%get_ptr(ddy, igrid)
            CALL ddz_f%get_ptr(ddz, igrid)

            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        volFl1 = volFl1 + vff(k,j,i) * ddx(i) * ddy(j) * ddz(k)
                    ENDDO
                ENDDO
            ENDDO
            
        ENDDO

        IF ( ABS(apprVol - volFl1) >= tol ) THEN
            WRITE(*,*) "Continuity equation is violated: volumeError = ", ABS(apprVol - volFl1)
        ENDIF

    END SUBROUTINE check_continuity

    !================================================================

    SUBROUTINE app_bcon()
    !----------------------------------------------------------------
    !   What it does:
    !    
    !----------------------------------------------------------------

        ! Subroutine arguments
        ! None

        ! Local variables
        TYPE(field_t), POINTER :: vff_f
        INTEGER(intk) :: ilevel

        CALL get_field(vff_f, "VFF")

        DO ilevel = minlevel, maxlevel
            CALL connect(ilevel, layers=2, s1=vff_f, corners=.TRUE.)
        ENDDO

    END SUBROUTINE app_bcon

    !================================================================

    SUBROUTINE comp_adve_quick(kk, jj, ii, q, l, u, v, w, &
        advr, vff, tol, adve)
    !----------------------------------------------------------------
    !   What it does:
    !   QUICK interpolation to compute the advected veloctiy
    !   (advectee) on staggered grid cells.
    !   adve = advected q (advectee)
    !   advr = advecting q (advector)
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
        INTEGER(intk), INTENT(in) :: kk, jj, ii, q, l
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: advr(kk, jj, ii)
        REAL(realk), INTENT(in) :: vff(kk,jj,ii)
        REAL(realk), INTENT(in) :: tol
        REAL(realk), INTENT(out) :: adve(kk, jj, ii)

        ! Loval variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: kl, jl, il
        REAL(realk) :: vel(kk,jj,ii)
        LOGICAL :: isIface(kk, jj, ii), isIfaceVic(kk, jj, ii)
        REAL(realk) :: signInd(2), iFacInd(2)

        CALL get_spatial_indices(kk, jj, ii, l, il, jl, kl)
        CALL get_condit_velocity(kk, jj, ii, q, u, v, w, vel)
        CALL track_iface(isIface, kk, jj, ii, vff, tol)
        CALL track_iface_vic(isIfaceVic, kk, jj, ii, isIface)

        DO i = 2, ii-2
            DO j = 2, jj-2
                DO k = 2, kk-2
                    signInd(1) = merge(1.0_realk, 0.0_realk, advr(k,j,i) >= 0.0_realk)
                    signInd(2) = 1.0_realk - signInd(1)
                    iFacInd(1) = merge(1.0_realk, 0.0_realk, isIfaceVic(k,j,i))
                    iFacInd(2) = 1.0_realk - iFacInd(1)

                    adve(k,j,i) = iFacInd(1) * ( signInd(1) * vel(k,j,i) + &
                                                 signInd(2) * vel(k+kl,j+jl,i+il) ) + &
                                  iFacInd(2) * ( signInd(1) * ( 0.750_realk * vel(k,j,i) + &
                                                                0.375_realk * vel(k+kl,j+jl,i+il) - &
                                                                0.125_realk * vel(k-kl,j-jl,i-il) ) + &
                                                 signInd(2) * ( 0.750_realk * vel(k+kl,j+jl,i+il) + &
                                                                0.375_realk * vel(k,j,i) - &
                                                                0.125_realk * vel(k+2*kl,j+2*jl,i+2*il) ) )
                END DO
            END DO
        END DO

    END SUBROUTINE comp_adve_quick

    !================================================================

    SUBROUTINE comp_advr_linear_interpolation(kk, jj, ii, l, u, v, w, advr)
    !----------------------------------------------------------------
    !   What it does:
    !   Linear interpolation to compute the advecting velocity
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
        INTEGER(intk), INTENT(in) :: kk, jj, ii, l
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(out) :: advr(kk, jj, ii)

        ! Loval variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: kl, jl, il
        REAL(realk) :: vel(kk,jj,ii)

        CALL get_spatial_indices(kk, jj, ii, l, il, jl, kl)
        CALL get_condit_velocity(kk, jj, ii, l, u, v, w, vel)

        DO i = 1, ii-1
            DO j = 1, jj-1
                DO k = 1, kk-1
                    advr(k,j,i) = 0.5_realk * ( vel(k,j,i) + vel(k+kl,j+jl,i+il) )
                END DO
            END DO
        END DO

    END SUBROUTINE comp_advr_linear_interpolation

    !================================================================

    SUBROUTINE get_spatial_indices(kk, jj, ii, lOrq, io, jo, ko)
    !----------------------------------------------------------------
    !   What it does:
    !   Depending on the direction (l) of component (q) the 
    !   indices io, jo and ko are set to 0 or 1. 
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, lOrq
        INTEGER(intk), INTENT(out) :: io, jo, ko

        ! Local variables
        ! None

        io = 0_intk ; jo = 0_intk ; ko = 0_intk

        IF ( lOrq == 1 ) THEN
            io = 1_intk
        ELSE IF ( lOrq == 2 ) THEN
            jo = 1_intk
        ELSE IF ( lOrq == 3 ) THEN
            ko = 1_intk
        ELSE
            CALL errr(__FILE__, __LINE__)
        END IF

    END SUBROUTINE get_spatial_indices

    !================================================================

    SUBROUTINE get_spatial_extents(kk, jj, ii, lOrq, dx, dy, dz, ddx, ddy, ddz, dsx, dsy, dsz)
    !----------------------------------------------------------------
    !   What it does:
    !   Depending on the direction (l) of component (q) the 
    !   extents ds(.) are set to d(.) or dd(.). The term ds(.) stands
    !   for spacing in (.)-direction and is a neutral specification
    !   for face-to-face or center-to-center distance.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, lOrq
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(out) :: dsx(ii), dsy(jj), dsz(kk)

        ! Local variables
        ! None

        dsx = ddx ; dsy = ddy ; dsz = ddz

        IF ( lOrq == 1 ) THEN
            dsx = dx
        ELSE IF ( lOrq == 2 ) THEN
            dsy = dy
        ELSE IF ( lOrq == 3 ) THEN
            dsz = dz
        ELSE
            CALL errr(__FILE__, __LINE__)
        END IF

    END SUBROUTINE get_spatial_extents

    !================================================================

    SUBROUTINE get_condit_velocity(kk, jj, ii, lOrq, u, v, w, vel)
    !----------------------------------------------------------------
    !   What it does:
    !   Depending on the direction (l) of component (q) the 
    !   specific velocity (vel) is set to u, v or w.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, lOrq
        REAL(realk), INTENT(in) :: u(kk,jj,ii), v(kk,jj,ii), w(kk,jj,ii)
        REAL(realk), INTENT(out) :: vel(kk,jj,ii)

        ! Local variables
        ! None

        vel = 0.0_realk

        IF ( lOrq == 1 ) THEN
            vel = u
        ELSE IF ( lOrq == 2 ) THEN
            vel = v
        ELSE IF ( lOrq == 3 ) THEN
            vel = w
        ELSE
            CALL errr(__FILE__, __LINE__)
        END IF

    END SUBROUTINE get_condit_velocity

END MODULE multiphase_vof_transport_mod
