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

    USE MPI_f08
    USE precision_mod, ONLY: intk, realk, mglet_mpi_real
    USE grids_mod, ONLY: nmygrids, mygrids
    USE field_mod, ONLY: field_t
    USE fields_mod, ONLY: get_field
    USE grids_mod, ONLY: get_mgdims, get_mgbasb, get_gradpxflag
    USE pointers_mod, ONLY: get_ip3
    USE multiphase_plic_mod, ONLY: comp_frac, iface_reconstruction, comp_stag_frac, track_iface, track_iface_vic
    USE multiphasecore_mod, ONLY: gmol1, gmol2, rho1, rho2, grav, permutation_multiphase, omitAdve, omitDiff, omitExte, tol, checkContinuity, checkSolenoidality, checkBalance, fluxLimiter, fluxCentered
    USE multiphase_material_mod, ONLY: comp_material_property_field, comp_property_face_value_cent, comp_property_face_value_stag
    USE flowcore_mod, ONLY: gradp
    USE connect2_mod, ONLY: connect
    USE parent_mod, ONLY: parent
    USE ftoc_mod, ONLY: ftoc
    USE grids_mod, ONLY: minlevel, maxlevel
    USE err_mod, ONLY: errr, err_abort
    USE multiphase_utils_mod, ONLY: get_spatial_indices, get_spatial_extents, get_condit_velocity, clip_vff, check_continuity, check_solenoidality, comp_vol_phase1, sanity_check

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

    SUBROUTINE comp_flux_cent(kk, jj, ii, l, vff, vel, vffFlux1, &
        ddx, ddy, ddz, normx, normy, normz, alpha, isIface, dt)
    !----------------------------------------------------------------
    !   What it does:
    !   Computes the volume fraction fluxes depending on the current
    !   split direction. It is distinguished between several cases 
    !   depending on the velocity direction and volume fraction 
    !   field.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: l
        REAL(realk), INTENT(in) :: vff(kk, jj, ii)
        REAL(realk), INTENT(in) :: vel(kk, jj, ii)
        REAL(realk), INTENT(out) :: vffFlux1(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: normx(kk, jj, ii), normy(kk, jj, ii), normz(kk, jj, ii)
        REAL(realk), INTENT(in) :: alpha(kk, jj, ii)
        LOGICAL, INTENT(in) :: isIface(kk, jj, ii)
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: kl, jl, il
        REAL(realk) :: dds, norms
        REAL(realk) :: dimx, dimy, dimz
        REAL(realk) :: fluxedProp, fluxWidth, fluxAlpha

        CALL get_spatial_indices(l, il, jl, kl)

        DO i = 2, ii-2
            DO j = 2, jj-2
                DO k = 2, kk-2
                    IF ( vel(k,j,i) > tol ) THEN

                        ! Compute face flux width and characteristic length
                        fluxWidth = abs( vel(k,j,i) ) * dt
                        dds = il * ddx(i) + jl * ddy(j) + kl * ddz(k)
                        IF ( fluxWidth > 0.5_realk*dds  ) THEN
                            CALL err_abort(155, "fluxWidth > 0.5*cellWidth! Hint: reduce dt", __FILE__, __LINE__)
                        ENDIF

                        IF ( isIface(k,j,i) ) THEN
                            ! Compute norm
                            norms = il * normx(k,j,i) + jl * normy(k,j,i) + kl * normz(k,j,i)

                            ! Compute proper alpha
                            fluxAlpha = alpha(k,j,i) - norms * ( dds - fluxWidth )

                            ! Compute dimensions of fluxed cuboid
                            dimx = il*fluxWidth + (1-il)*ddx(i)
                            dimy = jl*fluxWidth + (1-jl)*ddy(j)
                            dimz = kl*fluxWidth + (1-kl)*ddz(k)

                            ! Compute vff in fluxed cuboid
                            CALL comp_frac(fluxedProp, fluxAlpha, dimx, dimy, dimz, &
                                normx(k,j,i), normy(k,j,i), normz(k,j,i))
                        ELSE
                            fluxedProp = vff(k,j,i)
                        END IF
                    ELSE IF ( vel(k,j,i) < -tol ) THEN

                        ! Compute face flux width and characteristic length
                        fluxWidth = abs( vel(k,j,i) ) * dt
                        dds = il * ddx(i+il) + jl * ddy(j+jl) + kl * ddz(k+kl)
                        IF ( fluxWidth > 0.5_realk*dds  ) THEN
                            CALL err_abort(155, "fluxWidth > 0.5*cellWidth! Hint: reduce dt", __FILE__, __LINE__)
                        ENDIF

                        IF ( isIface(k+kl,j+jl,i+il) ) THEN
                            ! Compute proper alpha
                            fluxAlpha = alpha(k+kl,j+jl,i+il)

                            ! Compute dimensions of fluxed cuboid
                            dimx = il*fluxWidth + (1-il)*ddx(i+il)
                            dimy = jl*fluxWidth + (1-jl)*ddy(j+jl)
                            dimz = kl*fluxWidth + (1-kl)*ddz(k+kl)

                            ! Compute vff in fluxed cuboid
                            CALL comp_frac(fluxedProp, fluxAlpha, dimx, dimy, dimz, &
                                normx(k+kl,j+jl,i+il), normy(k+kl,j+jl,i+il), normz(k+kl,j+jl,i+il))
                        ELSE
                            fluxedProp = vff(k+kl,j+jl,i+il)
                        END IF
                    ELSE
                        fluxedProp = 0.0_realk
                    END IF
                    vffFlux1(k,j,i) = vel(k,j,i) * fluxedProp
                END DO
            END DO
        END DO

    END SUBROUTINE comp_flux_cent

    !================================================================

    SUBROUTINE comp_flux_stag(kk, jj, ii, q, l, vff, advr, &
        vffFlux1, vffFlux2, ddx, ddy, ddz, &
        normx, normy, normz, alpha, isIface, dt)
    !----------------------------------------------------------------
    !   What it does:
    !   Computes the volume fraction fluxes depending on the current
    !   split direction. It is distinguished between several cases 
    !   depending on the velocity direction and volume fraction 
    !   field.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, q, l
        REAL(realk), INTENT(in) :: vff(kk, jj, ii)
        REAL(realk), INTENT(in) :: advr(kk, jj, ii)
        REAL(realk), INTENT(out) :: vffFlux1(kk, jj, ii), vffFlux2(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: normx(kk, jj, ii), normy(kk, jj, ii), normz(kk, jj, ii)
        REAL(realk), INTENT(in) :: alpha(kk, jj, ii)
        LOGICAL, INTENT(in) :: isIface(kk, jj, ii)
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        INTEGER(intk) :: k, j, i, kl, jl, il, kq, jq, iq
        INTEGER(intk) :: kDonMi, jDonMi, iDonMi, kDonPl, jDonPl, iDonPl
        REAL(realk) :: ddsDon, normDon, fluxWidth, fluxAlpha, fluxDimx, fluxDimy, fluxDimz, fluxedProp
        REAL(realk) :: farEnd, ddslMi, ddslPl, ddsqMi, ddsqPl
        REAL(realk) :: normlMi, normlPl, normqMi, normqPl, fluxAlphaMi, fluxAlphaPl
        REAL(realk) :: fluxDimxMi, fluxDimyMi, fluxDimzMi, fluxDimxPl, fluxDimyPl, fluxDimzPl, fracMi, fracPl


        CALL get_spatial_indices(l, il, jl, kl)
        CALL get_spatial_indices(q, iq, jq, kq)

        IF ( q == l ) THEN
            DO i = 2, ii-2
                DO j = 2, jj-2
                    DO k = 2, kk-2
                        IF ( ABS(advr(k,j,i)) < tol ) THEN
                            vffFlux1(k,j,i) = 0.0_realk
                            vffFlux2(k,j,i) = 0.0_realk
                            CYCLE
                        ENDIF

                        ddsDon  = il*ddx(i+il) + jl*ddy(j+jl) + kl*ddz(k+kl)
                        normDon = il*normx(k+kl,j+jl,i+il) &
                                + jl*normy(k+kl,j+jl,i+il) &
                                + kl*normz(k+kl,j+jl,i+il)

                        fluxWidth = ABS( advr(k,j,i) ) * dt

                        IF ( fluxWidth > 0.5_realk*ddsDon ) THEN
                            CALL err_abort(155, "fluxWidth > 0.5*cellWidth! Hint: reduce dt", __FILE__, __LINE__)
                        ENDIF

                        IF ( isIface(k+kl,j+jl,i+il) ) THEN
                            fluxAlpha = alpha(k+kl,j+jl,i+il) - normDon * &
                                ( 0.5_realk*ddsDon &
                                - MERGE(fluxWidth, 0.0_realk, advr(k,j,i) > 0.0_realk) )
                            fluxDimx = il*fluxWidth + (1-il)*ddx(i)
                            fluxDimy = jl*fluxWidth + (1-jl)*ddy(j)
                            fluxDimz = kl*fluxWidth + (1-kl)*ddz(k)
                            CALL comp_frac(fluxedProp, fluxAlpha, fluxDimx, fluxDimy, fluxDimz, &
                                normx(k+kl,j+jl,i+il), normy(k+kl,j+jl,i+il), &
                                normz(k+kl,j+jl,i+il))
                        ELSE
                            fluxedProp = vff(k+kl,j+jl,i+il)
                        ENDIF

                        vffFlux1(k,j,i) = advr(k,j,i) * fluxedProp
                        vffFlux2(k,j,i) = advr(k,j,i) * ( 1.0_realk - fluxedProp )
                    END DO
                END DO
            END DO
        ELSE
            DO i = 2, ii-2
                DO j = 2, jj-2
                    DO k = 2, kk-2

                    IF ( ABS(advr(k,j,i)) < tol ) THEN
                        vffFlux1(k,j,i) = 0.0_realk
                        vffFlux2(k,j,i) = 0.0_realk
                        CYCLE
                    ENDIF

                    ! Indizes of the half-pressure-cells overlapping the velocity-cells
                    IF ( advr(k,j,i) > 0.0_realk ) THEN
                        iDonMi = i ; jDonMi = j ; kDonMi = k ; farEnd = 1.0_realk
                    ELSE
                        iDonMi = i+il ; jDonMi = j+jl ; kDonMi = k+kl ; farEnd = 0.0_realk
                    ENDIF

                    iDonPl = iDonMi+iq ; jDonPl = jDonMi+jq ; kDonPl = kDonMi+kq

                    fluxWidth = ABS( advr(k,j,i) ) * dt

                    ddslMi = il*ddx(iDonMi) + jl*ddy(jDonMi) + kl*ddz(kDonMi)
                    ddslPl = il*ddx(iDonPl) + jl*ddy(jDonPl) + kl*ddz(kDonPl)

                    IF ( fluxWidth > 0.5_realk*ddslMi .OR. fluxWidth > 0.5_realk*ddslPl ) THEN
                        CALL err_abort(155, "fluxWidth > 0.5*cellWidth! Hint: reduce dt", __FILE__, __LINE__)
                    ENDIF

                    ddsqMi = iq*ddx(iDonMi) + jq*ddy(jDonMi) + kq*ddz(kDonMi)
                    ddsqPl = iq*ddx(iDonPl) + jq*ddy(jDonPl) + kq*ddz(kDonPl)

                    normlMi = il*normx(kDonMi,jDonMi,iDonMi) + jl*normy(kDonMi,jDonMi,iDonMi) + kl*normz(kDonMi,jDonMi,iDonMi)
                    normlPl = il*normx(kDonPl,jDonPl,iDonPl) + jl*normy(kDonPl,jDonPl,iDonPl) + kl*normz(kDonPl,jDonPl,iDonPl)
                    normqMi = iq*normx(kDonMi,jDonMi,iDonMi) + jq*normy(kDonMi,jDonMi,iDonMi) + kq*normz(kDonMi,jDonMi,iDonMi)
                    normqPl = iq*normx(kDonPl,jDonPl,iDonPl) + jq*normy(kDonPl,jDonPl,iDonPl) + kq*normz(kDonPl,jDonPl,iDonPl)

                    IF ( isIface(kDonMi,jDonMi,iDonMi) ) THEN
                        fluxAlphaMi = alpha(kDonMi,jDonMi,iDonMi) &
                                    - normlMi * farEnd * ( ddslMi - fluxWidth ) &
                                    - normqMi * 0.5_realk * ddsqMi

                        fluxDimxMi = il*fluxWidth + iq*ddx(iDonMi)/2.0_realk + (1-il-iq)*ddx(iDonMi)
                        fluxDimyMi = jl*fluxWidth + jq*ddy(jDonMi)/2.0_realk + (1-jl-jq)*ddy(jDonMi)
                        fluxDimzMi = kl*fluxWidth + kq*ddz(kDonMi)/2.0_realk + (1-kl-kq)*ddz(kDonMi)

                        CALL comp_frac(fracMi, fluxAlphaMi, fluxDimxMi, fluxDimyMi, fluxDimzMi, &
                            normx(kDonMi,jDonMi,iDonMi), normy(kDonMi,jDonMi,iDonMi), normz(kDonMi,jDonMi,iDonMi))
                    ELSE
                        fracMi = vff(kDonMi,jDonMi,iDonMi)
                    ENDIF

                    IF ( isIface(kDonPl,jDonPl,iDonPl) ) THEN
                        fluxAlphaPl = alpha(kDonPl,jDonPl,iDonPl) &
                                    - normlPl * farEnd * ( ddslPl - fluxWidth )

                        fluxDimxPl = il*fluxWidth + iq*ddx(iDonPl)/2.0_realk + (1-il-iq)*ddx(iDonPl)
                        fluxDimyPl = jl*fluxWidth + jq*ddy(jDonPl)/2.0_realk + (1-jl-jq)*ddy(jDonPl)
                        fluxDimzPl = kl*fluxWidth + kq*ddz(kDonPl)/2.0_realk + (1-kl-kq)*ddz(kDonPl)

                        CALL comp_frac(fracPl, fluxAlphaPl, fluxDimxPl, fluxDimyPl, fluxDimzPl, &
                            normx(kDonPl,jDonPl,iDonPl), normy(kDonPl,jDonPl,iDonPl), &
                            normz(kDonPl,jDonPl,iDonPl))
                    ELSE
                        fracPl = vff(kDonPl,jDonPl,iDonPl)
                    ENDIF

                    fluxedProp = ( fracMi*ddsqMi + fracPl*ddsqPl ) / ( ddsqMi + ddsqPl )

                    vffFlux1(k,j,i) = advr(k,j,i) * fluxedProp
                    vffFlux2(k,j,i) = advr(k,j,i) * ( 1.0_realk - fluxedProp )

                    END DO
                END DO
            END DO
        ENDIF

    END SUBROUTINE comp_flux_stag

    !================================================================

    PURE SUBROUTINE def_advection_sequence(iteration, advSeq)
    !----------------------------------------------------------------
    !   What it does:
    !   Defines the sequence, in which the volume fraction vff and
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

    SUBROUTINE multiphase_solve(u_f, v_f, w_f, vff_f, p_f, dt, &
        itstep, uo_f, vo_f, wo_f)
    !----------------------------------------------------------------
    !   What it does:
    !    
    !----------------------------------------------------------------

        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: u_f
        TYPE(field_t), INTENT(inout) :: v_f
        TYPE(field_t), INTENT(inout) :: w_f
        TYPE(field_t), INTENT(inout) :: vff_f
        TYPE(field_t), INTENT(in) :: p_f
        REAL(realk), INTENT(in) :: dt
        INTEGER(intk), INTENT(in) :: itstep
        TYPE(field_t), INTENT(inout) :: uo_f, vo_f, wo_f

        ! Local variables
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: vff, p
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: uo, vo, wo

        TYPE(field_t), POINTER :: dx_f, dy_f, dz_f, ddx_f, ddy_f, ddz_f
        TYPE(field_t), POINTER :: rdx_f, rdy_f, rdz_f, rddx_f, rddy_f, rddz_f
        TYPE(field_t), POINTER :: up_f, vp_f, wp_f, vffp_f

        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:), ddx(:), ddy(:), ddz(:)
        REAL(realk), POINTER, CONTIGUOUS :: rdx(:), rdy(:), rdz(:), rddx(:), rddy(:), rddz(:)
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: up, vp, wp, vffp

        INTEGER(intk) :: i, igrid
        INTEGER(intk) :: kk, jj, ii

        IF ( checkContinuity ) CALL check_continuity(itstep)
        IF ( checkSolenoidality ) CALL check_solenoidality(itstep, dt)

        uo_f = 0.0_realk
        vo_f = 0.0_realk
        wo_f = 0.0_realk

        CALL get_field(dx_f, "DX"); CALL get_field(dy_f, "DY"); CALL get_field(dz_f, "DZ")
        CALL get_field(ddx_f, "DDX"); CALL get_field(ddy_f, "DDY"); CALL get_field(ddz_f, "DDZ")

        CALL get_field(rdx_f, "RDX"); CALL get_field(rdy_f, "RDY"); CALL get_field(rdz_f, "RDZ")
        CALL get_field(rddx_f, "RDDX"); CALL get_field(rddy_f, "RDDY"); CALL get_field(rddz_f, "RDDZ")

        CALL get_field(up_f, "UP"); CALL get_field(vp_f, "VP"); CALL get_field(wp_f, "WP")
        CALL get_field(vffp_f, "VFFP")

        up_f%arr = u_f%arr
        vp_f%arr = v_f%arr
        wp_f%arr = w_f%arr
        vffp_f%arr = vff_f%arr

        CALL adve_operator(u_f, v_f, w_f, vff_f, dt, itstep)

        DO i = 1, nmygrids
            igrid = mygrids(i)
            CALL get_mgdims(kk, jj, ii, igrid)

            CALL vff_f%get_ptr(vff, igrid)
            CALL p_f%get_ptr(p, igrid)
            CALL uo_f%get_ptr(uo, igrid)
            CALL vo_f%get_ptr(vo, igrid)
            CALL wo_f%get_ptr(wo, igrid)

            CALL up_f%get_ptr(up, igrid)
            CALL vp_f%get_ptr(vp, igrid)
            CALL wp_f%get_ptr(wp, igrid)
            CALL vffp_f%get_ptr(vffp, igrid)

            CALL dx_f%get_ptr(dx, igrid)
            CALL dy_f%get_ptr(dy, igrid)
            CALL dz_f%get_ptr(dz, igrid)
            CALL ddx_f%get_ptr(ddx, igrid)
            CALL ddy_f%get_ptr(ddy, igrid)
            CALL ddz_f%get_ptr(ddz, igrid)
            CALL rdx_f%get_ptr(rdx, igrid)
            CALL rdy_f%get_ptr(rdy, igrid)
            CALL rdz_f%get_ptr(rdz, igrid)
            CALL rddx_f%get_ptr(rddx, igrid)
            CALL rddy_f%get_ptr(rddy, igrid)
            CALL rddz_f%get_ptr(rddz, igrid)

            CALL diff_operator(kk, jj, ii, up, vp, wp, vffp, vff, ddx, ddy, ddz, &
                rdx, rdy, rdz, rddx, rddy, rddz, uo, vo, wo)
            CALL pres_operator(kk, jj, ii, vff, p, ddx, ddy, ddz, &
                rdx, rdy, rdz, igrid, uo, vo, wo)
            CALL exte_operator(kk, jj, ii, uo, vo, wo)
        END DO

    END SUBROUTINE multiphase_solve

    !================================================================

    SUBROUTINE adve_operator(u_f, v_f, w_f, vff_f, dt, itstep)
    !----------------------------------------------------------------
    !   What it does:
    !   Performs the spatial and temporal integration of the
    !   advection operator. The integration is combined, since
    !   VOF/PLIC is "exact" up to the order of accuracy of the 
    !   interface reconstruction in PLIC.
    !   In this routine four types of fields exist.
    !
    !   Type 1, truly time-persistant fields:
    !       - exist over multiple time-steps
    !       - previous time-steps status is important
    !       - initialized in core/multiphasecore with set_field
    !       Examples: u, v, w, p, vff
    !   Type 2, falsely time-persistant fields:
    !       - exist over multiple time-steps
    !       - previous time-steps status is overwritten
    !       - final result is needed for post-processing
    !       - initialized in core/multiphasecore with set_field
    !       Examples: normx, normy, normz and alpha
    !   Type 3, sweep-persistant fields:
    !       - exist over multiple sweeps
    !       - previous time-steps status is overwritten
    !       - initialized in this routine
    !       Examples: vffStag(q), mom(q), cWyStag(q), cWy, ...
    !   Type 4, sweep-temporary fields:
    !       - exist only during one sweep
    !       - previous sweeps status is overwritten
    !       - initialized for each grid
    !       Examples: dStag, vffFlux1, vffFlux1Stag, ...
    !
    !   The reason for this "clumsy" approach is, the connection of
    !   multiple grids. MGELT only uses two hollow cells at the grids
    !   boundaries. For the multi-phase algorithm this is not enough.
    !   Hence, after each sweep, the hollow cells need to be
    !   connected. The introduction of bigger buffer layers could
    !   improve the performance of advection algorithm.
    !----------------------------------------------------------------

        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: u_f
        TYPE(field_t), INTENT(inout) :: v_f
        TYPE(field_t), INTENT(inout) :: w_f
        TYPE(field_t), INTENT(inout) :: vff_f
        REAL(realk), INTENT(in) :: dt
        INTEGER(intk), INTENT(in) :: itstep

        ! Local variables
        INTEGER(intk) :: q, l, advSeq(3), splitDir
        TYPE(field_t), POINTER :: dx_f, dy_f, dz_f, ddx_f, ddy_f, ddz_f
        TYPE(field_t), POINTER :: normx_f, normy_f, normz_f, alpha_f
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:), ddx(:), ddy(:), ddz(:)
        REAL(realk), POINTER, CONTIGUOUS :: normx(:,:,:), normy(:,:,:), normz(:,:,:), alpha(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: u(:,:,:), v(:,:,:), w(:,:,:), vff(:,:,:)
        TYPE(field_t) :: vffStag(3), mom(3), cWyStag(3), cWy
        CHARACTER(len=1), PARAMETER :: component(3) = ['X','Y','Z']
        INTEGER(intk) :: n, igrid, kk, jj, ii, ip3, ilevel
        REAL(realk), ALLOCATABLE :: dStag(:,:,:), advr(:,:,:), adve(:,:,:)
        REAL(realk), ALLOCATABLE :: vffFlux1Stag(:,:,:), vffFlux2Stag(:,:,:), vel(:,:,:), vffFlux1(:,:,:)
        LOGICAL, ALLOCATABLE :: isIface(:,:,:)
        REAL(realk) :: volPhase1r, volPhase1r1, volFluxPhase1, volCompPhase1, volClipPhase1

        IF ( omitAdve ) RETURN

        volPhase1r = 0.0_realk
        volPhase1r1 = 0.0_realk
        volFluxPhase1 = 0.0_realk
        volCompPhase1 = 0.0_realk
        volClipPhase1 = 0.0_realk

        ! Get missing truly time-persistant fields
        CALL get_field(dx_f, "DX")
        CALL get_field(dy_f, "DY")
        CALL get_field(dz_f, "DZ")
        CALL get_field(ddx_f, "DDX")
        CALL get_field(ddy_f, "DDY")
        CALL get_field(ddz_f, "DDZ")

        ! Get falsely time-persistant fields
        CALL get_field(normx_f, "NORMX")
        CALL get_field(normy_f, "NORMY")
        CALL get_field(normz_f, "NORMZ")
        CALL get_field(alpha_f, "ALPHA")

        ! Initialize sweep persistant fields
        DO q = 1, 3
            CALL vffStag(q)%init("VFFSTAG"//component(q))
            CALL vffStag(q)%init_buffers()
            CALL mom(q)%init("MOM"//component(q))
            CALL mom(q)%init_buffers()
            CALL cWyStag(q)%init("CWYSTAG"//component(q))
        END DO
        CALL cWy%init("CWY")

        DO n = 1, nmygrids
            igrid = mygrids(n)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)

            ! Get pointers to truly time-persistant fields
            CALL vff_f%get_ptr(vff, igrid)
            CALL dx_f%get_ptr(dx, igrid)
            CALL dy_f%get_ptr(dy, igrid)
            CALL dz_f%get_ptr(dz, igrid)
            CALL ddx_f%get_ptr(ddx, igrid)
            CALL ddy_f%get_ptr(ddy, igrid)
            CALL ddz_f%get_ptr(ddz, igrid)

            ! Get pointers to falsely time-persistant fields
            CALL normx_f%get_ptr(normx, igrid)
            CALL normy_f%get_ptr(normy, igrid)
            CALL normz_f%get_ptr(normz, igrid)
            CALL alpha_f%get_ptr(alpha, igrid)

            ! Interface reconstruction and Weymouth-Yue-Coefficient on all grids
            CALL iface_reconstruction(kk, jj, ii, vff, dx, dy, dz, &
                ddx, ddy, ddz, normx, normy, normz, alpha)
            CALL comp_cWy(kk, jj, ii, vff, cWy%arr(ip3))

            ! Compute Volume of Phase 1 at rk-step r
            CALL comp_vol_phase1(kk, jj, ii, vff, ddx, ddy, ddz, volPhase1r)
        ENDDO

        CALL MPI_Allreduce(MPI_IN_PLACE, volPhase1r, 1, mglet_mpi_real, MPI_SUM, MPI_COMM_WORLD)

        DO n = 1, nmygrids
            igrid = mygrids(n)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)

            ! Get pointers to truly time-persistant fields
            CALL u_f%get_ptr(u, igrid)
            CALL v_f%get_ptr(v, igrid)
            CALL w_f%get_ptr(w, igrid)
            CALL vff_f%get_ptr(vff, igrid)
            CALL ddx_f%get_ptr(ddx, igrid)
            CALL ddy_f%get_ptr(ddy, igrid)
            CALL ddz_f%get_ptr(ddz, igrid)

            ! Get pointers to falsely time-persistant fields
            CALL normx_f%get_ptr(normx, igrid)
            CALL normy_f%get_ptr(normy, igrid)
            CALL normz_f%get_ptr(normz, igrid)
            CALL alpha_f%get_ptr(alpha, igrid)

            ! Allocate sweep-temporary fields
            ALLOCATE(dStag(kk, jj, ii))
            ALLOCATE(isIface(kk, jj, ii))

            CALL track_iface(isIface, kk, jj, ii, vff)

            DO q = 1, 3
                CALL comp_stag_frac(kk, jj, ii, q, vff, vffStag(q)%arr(ip3), &
                    ddx, ddy, ddz, normx, normy, normz, alpha, isIface)
                CALL comp_material_property_field(kk, jj, ii, vffStag(q)%arr(ip3), &
                    dStag, rho1, rho2)
                CALL comp_momentum(kk, jj, ii, q, dStag, u, v, w, mom(q)%arr(ip3))
                CALL comp_cWy(kk, jj, ii, vffStag(q)%arr(ip3), cWyStag(q)%arr(ip3))
            END DO

            ! Deallocate sweep-temporary fields
            DEALLOCATE(isIface)
            DEALLOCATE(dStag)
        ENDDO

        ! Fine to coarse within the domain (ftoc)
        DO ilevel = maxlevel, minlevel, -1
            CALL ftoc(ilevel, vffStag(1)%arr, vffStag(1)%arr, 'A')
            CALL ftoc(ilevel, vffStag(2)%arr, vffStag(2)%arr, 'B')
            CALL ftoc(ilevel, vffStag(3)%arr, vffStag(3)%arr, 'C')
            CALL ftoc(ilevel, vff_f%arr, vff_f%arr, 'D')
            CALL ftoc(ilevel, mom(1)%arr, mom(1)%arr, 'A')
            CALL ftoc(ilevel, mom(2)%arr, mom(2)%arr, 'B')
            CALL ftoc(ilevel, mom(3)%arr, mom(3)%arr, 'C')
        ENDDO

        ! Coarse to fine in the buffer (parent)
        ! Fill buffers of vff and vffStag with initial values (connect)
        DO ilevel = minlevel, maxlevel
            CALL parent(ilevel, s1=vff_f)
            CALL parent(ilevel, v1=vffStag(1), v2=vffStag(2), v3=vffStag(3))
            CALL parent(ilevel, v1=mom(1), v2=mom(2), v3=mom(3))
            CALL connect(ilevel, 2, s1=vff_f, corners=.TRUE.)
            CALL connect(ilevel, 2, v1=vffStag(1), v2=vffStag(2), v3=vffStag(3), corners=.TRUE.)
            CALL connect(ilevel, 2, v1=mom(1), v2=mom(2), v3=mom(3), corners=.TRUE.)
        END DO

        CALL def_advection_sequence(itstep, advSeq)

        DO splitDir = 1, 3
            l = advSeq(splitDir)
            DO n = 1, nmygrids
                igrid = mygrids(n)
                CALL get_mgdims(kk, jj, ii, igrid)
                CALL get_ip3(ip3, igrid)

                ! Get pointers to truly time-persistant fields
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

                ! Get pointers to falsely time-persistant fields
                CALL normx_f%get_ptr(normx, igrid)
                CALL normy_f%get_ptr(normy, igrid)
                CALL normz_f%get_ptr(normz, igrid)
                CALL alpha_f%get_ptr(alpha, igrid)

                ! Allocate sweep-temporary fields
                ALLOCATE(advr(kk, jj, ii))
                ALLOCATE(adve(kk, jj, ii))
                ALLOCATE(vffFlux1Stag(kk, jj, ii))
                ALLOCATE(vffFlux2Stag(kk, jj, ii))
                ALLOCATE(isIface(kk, jj, ii))
                ALLOCATE(vel(kk, jj, ii))
                ALLOCATE(vffFlux1(kk, jj, ii))

                CALL track_iface(isIface, kk, jj, ii, vff)

                DO q = 1, 3
                    CALL comp_advr_linear_interpolation(kk, jj, ii, q, l, u, v, w, &
                        advr, ddx, ddy, ddz)
                    CALL comp_adve(kk, jj, ii, q, l, vffStag(q)%arr(ip3), u, v, w, &
                        advr, adve, dx, dy, dz, ddx, ddy, ddz, dt)
                    CALL comp_flux_stag(kk, jj, ii, q, l, vff, advr, vffFlux1Stag, vffFlux2Stag, &
                        ddx, ddy, ddz, normx, normy, normz, alpha, isIface, dt)
                    CALL adv_mom(kk, jj, ii, q, l, vffStag(q)%arr(ip3), cWyStag(q)%arr(ip3), u, v, w, &
                        advr, adve, mom(q)%arr(ip3), vffFlux1Stag, vffFlux2Stag, dx, dy, dz, ddx, ddy, ddz, dt)
                    CALL adv_vof(kk, jj, ii, q, l, vffStag(q)%arr(ip3), cWyStag(q)%arr(ip3), &
                        advr, vffFlux1Stag, dx, dy, dz, ddx, ddy, ddz, dt)
                END DO

                CALL get_condit_velocity(kk, jj, ii, l, u, v, w, vel)
                CALL comp_flux_cent(kk, jj, ii, l, vff, vel, vffFlux1, &
                    ddx, ddy, ddz, normx, normy, normz, alpha, isIface, dt)
                CALL adv_vof(kk, jj, ii, 0, l, vff, cWy%arr(ip3), vel, vffFlux1, &
                    dx, dy, dz, ddx, ddy, ddz, dt, volFluxPhase1, volCompPhase1)

                ! Deallocate sweep-temporary fields
                DEALLOCATE(vffFlux1)
                DEALLOCATE(vel)
                DEALLOCATE(isIface)
                DEALLOCATE(vffFlux2Stag)
                DEALLOCATE(vffFlux1Stag)
                DEALLOCATE(adve)
                DEALLOCATE(advr)
            ENDDO

            ! Fine to coarse within the domain (ftoc)
            DO ilevel = maxlevel, minlevel, -1
                CALL ftoc(ilevel, vffStag(1)%arr, vffStag(1)%arr, 'A')
                CALL ftoc(ilevel, vffStag(2)%arr, vffStag(2)%arr, 'B')
                CALL ftoc(ilevel, vffStag(3)%arr, vffStag(3)%arr, 'C')
                CALL ftoc(ilevel, vff_f%arr, vff_f%arr, 'D')
                CALL ftoc(ilevel, mom(1)%arr, mom(1)%arr, 'A')
                CALL ftoc(ilevel, mom(2)%arr, mom(2)%arr, 'B')
                CALL ftoc(ilevel, mom(3)%arr, mom(3)%arr, 'C')
            ENDDO

            ! Coarse to fine in the buffer (parent)
            ! Fill buffers of vff and vffStag after l sweep (connect)
            DO ilevel = minlevel, maxlevel
                CALL parent(ilevel, s1=vff_f)
                CALL parent(ilevel, v1=vffStag(1), v2=vffStag(2), v3=vffStag(3))
                CALL parent(ilevel, v1=mom(1), v2=mom(2), v3=mom(3))
                CALL connect(ilevel, 2, s1=vff_f, corners=.TRUE.)
                CALL connect(ilevel, 2, v1=vffStag(1), v2=vffStag(2), v3=vffStag(3), corners=.TRUE.)
                CALL connect(ilevel, 2, v1=mom(1), v2=mom(2), v3=mom(3), corners=.TRUE.)
            END DO

            DO n = 1, nmygrids
                igrid = mygrids(n)
                CALL get_mgdims(kk, jj, ii, igrid)

                ! Get pointers to truly time-persistant fields
                CALL vff_f%get_ptr(vff, igrid)
                CALL dx_f%get_ptr(dx, igrid)
                CALL dy_f%get_ptr(dy, igrid)
                CALL dz_f%get_ptr(dz, igrid)
                CALL ddx_f%get_ptr(ddx, igrid)
                CALL ddy_f%get_ptr(ddy, igrid)
                CALL ddz_f%get_ptr(ddz, igrid)

                ! Get pointers to falsely time-persistant fields
                CALL normx_f%get_ptr(normx, igrid)
                CALL normy_f%get_ptr(normy, igrid)
                CALL normz_f%get_ptr(normz, igrid)
                CALL alpha_f%get_ptr(alpha, igrid)

                ! Interface reconstruction
                CALL iface_reconstruction(kk, jj, ii, vff, dx, dy, dz, ddx, ddy, ddz, normx, normy, normz, alpha)
            ENDDO
        END DO

        CALL MPI_Allreduce(MPI_IN_PLACE, volFluxPhase1, 1, mglet_mpi_real, MPI_SUM, MPI_COMM_WORLD)
        CALL MPI_Allreduce(MPI_IN_PLACE, volCompPhase1, 1, mglet_mpi_real, MPI_SUM, MPI_COMM_WORLD)

        DO n = 1, nmygrids
            igrid = mygrids(n)
            CALL get_mgdims(kk, jj, ii, igrid)

            ! Get pointers to truly time-persistant fields
            CALL vff_f%get_ptr(vff, igrid)
            CALL ddx_f%get_ptr(ddx, igrid)
            CALL ddy_f%get_ptr(ddy, igrid)
            CALL ddz_f%get_ptr(ddz, igrid)

            ! Clip vff at rk-step r+1
            CALL clip_vff(kk, jj, ii, vff, ddx, ddy, ddz, volClipPhase1)

            ! Compute Volume of Phase 1 at rk-step r+1
            CALL comp_vol_phase1(kk, jj, ii, vff, ddx, ddy, ddz, volPhase1r1)
        ENDDO

        CALL MPI_Allreduce(MPI_IN_PLACE, volClipPhase1, 1, mglet_mpi_real, MPI_SUM, MPI_COMM_WORLD)
        CALL MPI_Allreduce(MPI_IN_PLACE, volPhase1r1, 1, mglet_mpi_real, MPI_SUM, MPI_COMM_WORLD)

        IF ( checkBalance ) CALL sanity_check(volPhase1r, volPhase1r1, volFluxPhase1, volCompPhase1, volClipPhase1)

        DO n = 1, nmygrids
            igrid = mygrids(n)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_ip3(ip3, igrid)

            CALL u_f%get_ptr(u, igrid)
            CALL v_f%get_ptr(v, igrid)
            CALL w_f%get_ptr(w, igrid)

            DO q = 1, 3
                CALL update_velocity(kk, jj, ii, q, vffStag(q)%arr(ip3), u, v, w, mom(q)%arr(ip3), dt)
            END DO
        ENDDO

        ! Finish sweep persistant fields
        CALL cWy%finish()
        DO q = 1, 3
            CALL cWyStag(q)%finish()
            CALL mom(q)%finish()
            CALL vffStag(q)%finish()
        END DO

    END SUBROUTINE adve_operator

    !================================================================

    SUBROUTINE diff_operator(kk, jj, ii, u, v, w, vffp, vff, &
        ddx, ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz, uo, vo, wo)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------
    
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: vffp(kk, jj, ii), vff(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: rdx(ii), rdy(jj), rdz(kk)
        REAL(realk), INTENT(in) :: rddx(ii), rddy(jj), rddz(kk)
        REAL(realk), INTENT(inout) :: uo(kk, jj, ii), vo(kk, jj, ii), wo(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i
        REAL(realk) :: g(kk, jj, ii)
        REAL(realk) :: gxy(kk, jj, ii), gxz(kk, jj, ii), gyz(kk, jj, ii)
        REAL(realk) :: rhoe(kk, jj, ii), rhon(kk, jj, ii), rhot(kk, jj, ii)
        REAL(realk) :: tauxxe, tauxxw, tauxyn, tauxys, tauxzt, tauxzb
        REAL(realk) :: tauyxe, tauyxw, tauyyn, tauyys, tauyzt, tauyzb
        REAL(realk) :: tauzxe, tauzxw, tauzyn, tauzys, tauzzt, tauzzb

        IF ( omitDiff ) RETURN

        CALL comp_material_property_field(kk, jj, ii, vffp, g, gmol1, gmol2)
        CALL comp_property_face_value_stag(kk, jj, ii, vffp, gmol1, gmol2, &
            gxy, gxz, gyz, ddx, ddy, ddz)
        CALL comp_property_face_value_cent(kk, jj, ii, vff, rho1, rho2, 'ARI', &
            rhoe, rhon, rhot, ddx, ddy, ddz)

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    ! Stresses
                    tauxxe = g(k,j,i+1) * 2.0_realk * (u(k,j,i+1) - u(k,j,i))*rddx(i+1)
                    tauxxw = g(k,j,i) * 2.0_realk * (u(k,j,i) - u(k,j,i-1))*rddx(i)
                    tauxyn = gxy(k,j,i) * ((u(k,j+1,i) - u(k,j,i))*rdy(j) + (v(k,j,i+1) - v(k,j,i))*rdx(i))
                    tauxys = gxy(k,j-1,i) * ((u(k,j,i) - u(k,j-1,i))*rdy(j-1) + (v(k,j-1,i+1) - v(k,j-1,i))*rdx(i))
                    tauxzt = gxz(k,j,i) * ((u(k+1,j,i) - u(k,j,i))*rdz(k) + (w(k,j,i+1) - w(k,j,i))*rdx(i))
                    tauxzb = gxz(k-1,j,i) * ((u(k,j,i) - u(k-1,j,i))*rdz(k-1) + (w(k-1,j,i+1) - w(k-1,j,i))*rdx(i))

                    ! Change due to diffusion
                    uo(k,j,i) = uo(k,j,i) + 1.0_realk/rhoe(k,j,i)* &
                        ((tauxxe - tauxxw)*rdx(i) + &
                        (tauxyn - tauxys)*rddy(j) + &
                        (tauxzt - tauxzb)*rddz(k))
                END DO
            END DO
        END DO

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    ! Stresses
                    tauyxe = gxy(k,j,i) * ((u(k,j+1,i) - u(k,j,i))*rdy(j) + (v(k,j,i+1) - v(k,j,i))*rdx(i))
                    tauyxw = gxy(k,j,i-1) * ((u(k,j+1,i-1) - u(k,j,i-1))*rdy(j) + (v(k,j,i) - v(k,j,i-1))*rdx(i-1))
                    tauyyn = g(k,j+1,i) * 2.0_realk * (v(k,j+1,i) - v(k,j,i))*rddy(j+1)
                    tauyys = g(k,j,i) * 2.0_realk * (v(k,j,i) - v(k,j-1,i))*rddy(j)
                    tauyzt = gyz(k,j,i) * ((v(k+1,j,i) - v(k,j,i))*rdz(k) + (w(k,j+1,i) - w(k,j,i))*rdy(j))
                    tauyzb = gyz(k-1,j,i) * ((v(k,j,i) - v(k-1,j,i))*rdz(k-1) + (w(k-1,j+1,i) - w(k-1,j,i))*rdy(j))

                    ! Change due to diffusion
                    vo(k,j,i) = vo(k,j,i) + 1.0_realk/rhon(k,j,i)* &
                        ((tauyxe - tauyxw)*rddx(i) + &
                        (tauyyn - tauyys)*rdy(j) + &
                        (tauyzt - tauyzb)*rddz(k))
                END DO
            END DO
        END DO

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    ! Stresses
                    tauzxe = gxz(k,j,i) * ((u(k+1,j,i) - u(k,j,i))*rdz(k) + (w(k,j,i+1) - w(k,j,i))*rdx(i))
                    tauzxw = gxz(k,j,i-1) * ((u(k+1,j,i-1) - u(k,j,i-1))*rdz(k) + (w(k,j,i) - w(k,j,i-1))*rdx(i-1))
                    tauzyn = gyz(k,j,i) * ((v(k+1,j,i) - v(k,j,i))*rdz(k) + (w(k,j+1,i) - w(k,j,i))*rdy(j))
                    tauzys = gyz(k,j-1,i) * ((v(k+1,j-1,i) - v(k,j-1,i))*rdz(k) + (w(k,j,i) - w(k,j-1,i))*rdy(j-1))
                    tauzzt = g(k+1,j,i) * 2.0_realk * (w(k+1,j,i) - w(k,j,i))*rddz(k+1)
                    tauzzb = g(k,j,i) * 2.0_realk * (w(k,j,i) - w(k-1,j,i))*rddz(k)

                    ! Change due to diffusion
                    wo(k,j,i) = wo(k,j,i) + 1.0_realk/rhot(k,j,i)*&
                        ((tauzxe - tauzxw)*rddx(i) + &
                        (tauzyn - tauzys)*rddy(j) + &
                        (tauzzt - tauzzb)*rdz(k))
                END DO
            END DO
        END DO

    END SUBROUTINE diff_operator

    !================================================================

    SUBROUTINE pres_operator(kk, jj, ii, vff, p, ddx, ddy, ddz, rdx, rdy, rdz, &
        igrid, uo, vo, wo)
    !----------------------------------------------------------------
    !   What it does:
    !    
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: vff(kk, jj, ii), p(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: rdx(ii), rdy(jj), rdz(kk)
        INTEGER(intk), INTENT(in) :: igrid
        REAL(realk), INTENT(inout) :: uo(kk, jj, ii), vo(kk, jj, ii), wo(kk, jj, ii)

        ! Local variables
        REAL(realk) :: rhoe(kk, jj, ii), rhon(kk, jj, ii), rhot(kk, jj, ii)
        INTEGER(intk) :: gradpflag
        REAL(realk) :: gpx(kk, jj, ii), gpy(kk, jj, ii), gpz(kk, jj, ii)
        INTEGER(intk) :: i, j, k

        CALL comp_property_face_value_cent(kk, jj, ii, vff, rho1, rho2, 'ARI', &
            rhoe, rhon, rhot, ddx, ddy, ddz)

        CALL get_gradpxflag(gradpflag, igrid)

        gpx = 0.0_realk
        gpy = 0.0_realk
        gpz = 0.0_realk
        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    gpx(k,j,i) = gradp(1)*gradpflag*MERGE(1.0_realk, 0.0_realk, vff(k,j,i) > tol)
                    gpy(k,j,i) = gradp(2)*gradpflag*MERGE(1.0_realk, 0.0_realk, vff(k,j,i) > tol)
                    gpz(k,j,i) = gradp(3)*gradpflag*MERGE(1.0_realk, 0.0_realk, vff(k,j,i) > tol)
                ENDDO
            ENDDO
        ENDDO

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    uo(k,j,i) = uo(k,j,i) - 1.0_realk / rhoe(k,j,i) * ( ( p(k,j,i+1) - p(k,j,i) ) * rdx(i) + gpx(k,j,i) )
                END DO
            END DO
        END DO

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    vo(k,j,i) = vo(k,j,i) - 1.0_realk / rhon(k,j,i) * ( ( p(k,j+1,i) - p(k,j,i) ) * rdy(j) + gpy(k,j,i) )
                END DO
            END DO
        END DO

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    wo(k,j,i) = wo(k,j,i) - 1.0_realk / rhot(k,j,i) * ( ( p(k+1,j,i) - p(k,j,i) ) * rdz(k) + gpz(k,j,i) )
                END DO
            END DO
        END DO

    END SUBROUTINE pres_operator

    !================================================================

    SUBROUTINE exte_operator(kk, jj, ii, uo, vo, wo)
    !----------------------------------------------------------------
    !   What it does:
    !    
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: uo(kk, jj, ii), vo(kk, jj, ii), wo(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: i, j, k

        IF ( omitExte ) RETURN

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    uo(k,j,i) = uo(k,j,i) + grav(1)
                    vo(k,j,i) = vo(k,j,i) + grav(2)
                    wo(k,j,i) = wo(k,j,i) + grav(3)
                ENDDO
            ENDDO
        ENDDO

    END SUBROUTINE exte_operator

    !================================================================

    SUBROUTINE adv_vof(kk, jj, ii, q, l, vff, cWy, vel, vffFlux1, &
        dx, dy, dz, ddx, ddy, ddz, dt, volFluxPhase1, volCompPhase1)
    !----------------------------------------------------------------
    !   What it does:
    !    
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, q, l
        REAL(realk), INTENT(inout) :: vff(kk, jj, ii)
        REAL(realk), INTENT(in) :: cWy(kk, jj, ii)
        REAL(realk), INTENT(in) :: vel(kk,jj,ii)
        REAL(realk), INTENT(in) :: vffFlux1(kk, jj, ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: dt
        REAL(realk), INTENT(inout), OPTIONAL :: volFluxPhase1, volCompPhase1

        ! Local variables
        INTEGER(intk) :: i, j, k
        INTEGER(intk) :: il, jl, kl
        REAL(realk) :: dsx(ii), dsy(jj), dsz(kk), dsCV, dV
        REAL(realk) :: div(kk, jj, ii)

        CALL get_spatial_indices(l, il, jl, kl)
        CALL get_spatial_extents(kk, jj, ii, q, l, dx, dy, dz, ddx, ddy, ddz, dsx, dsy, dsz)

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    dsCV = il * dsx(i) + jl * dsy(j) + kl * dsz(k)
                    div(k,j,i) = ( vel(k,j,i) - vel(k-kl,j-jl,i-il) ) / dsCV

                    IF ( PRESENT(volFluxPhase1) .AND. PRESENT(volCompPhase1) ) THEN
                        dV   = ddx(i)*ddy(j)*ddz(k)
                        volFluxPhase1 = volFluxPhase1 - dt/dsCV * ( vffFlux1(k,j,i) - vffFlux1(k-kl,j-jl,i-il) ) * dV
                        volCompPhase1 = volCompPhase1 + dt * cWy(k,j,i) * div(k,j,i) * dV
                    ENDIF

                    vff(k,j,i) = vff(k,j,i) &
                        - dt/dsCV * ( vffFlux1(k,j,i) - vffFlux1(k-kl,j-jl,i-il) ) &
                        + dt * cWy(k,j,i) * div(k,j,i)
                END DO
            END DO
        END DO

    END SUBROUTINE adv_vof

    !================================================================

    SUBROUTINE adv_mom(kk, jj, ii, q, l, vff, cWy, u, v, w, &
        advr, adve, mom, vffFlux1, vffFlux2, dx, dy, dz, &
        ddx, ddy, ddz, dt)
    !----------------------------------------------------------------
    !   What it does:
    !    
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, q, l
        REAL(realk), INTENT(in) :: vff(kk, jj, ii), cWy(kk, jj, ii)
        REAL(realk), INTENT(in) :: u(kk,jj,ii), v(kk,jj,ii), w(kk,jj,ii)
        REAL(realk), INTENT(in) :: advr(kk,jj,ii), adve(kk,jj,ii)
        REAL(realk), INTENT(inout) :: mom(kk, jj, ii)
        REAL(realk), INTENT(in) :: vffFlux1(kk, jj, ii), vffFlux2(kk, jj, ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        INTEGER(intk) :: i, j, k
        INTEGER(intk) :: il, jl, kl
        REAL(realk) :: dsx(ii), dsy(jj), dsz(kk), dsCV
        REAL(realk) :: vel(kk,jj,ii)
        REAL(realk) :: momFlux(kk, jj, ii)
        REAL(realk) :: div, com

        CALL get_spatial_indices(l, il, jl, kl)
        CALL get_spatial_extents(kk, jj, ii, q, l, dx, dy, dz, ddx, ddy, ddz, dsx, dsy, dsz)
        CALL get_condit_velocity(kk, jj, ii, q, u, v, w, vel)

        DO i = 2, ii-2
            DO j = 2, jj-2
                DO k = 2, kk-2
                    momFlux(k,j,i) = adve(k,j,i) * ( rho1 * vffFlux1(k,j,i) + rho2 * vffFlux2(k,j,i) )
                END DO
            END DO 
        END DO

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    dsCV = il * dsx(i) + jl * dsy(j) + kl * dsz(k)
                    div = ( advr(k,j,i) - advr(k-kl,j-jl,i-il) ) / dsCV
                    com = ( rho1 * cWy(k,j,i) + rho2 * ( 1.0_realk - cWy(k,j,i) ) ) * div

                    mom(k,j,i) = mom(k,j,i) &
                        - dt/dsCV * ( momFlux(k,j,i) - momFlux(k-kl,j-jl,i-il) ) &
                        + dt * vel(k,j,i) * com
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

        DO i = 2, ii-2
            DO j = 2, jj-2
                DO k = 2, kk-2
                    mom(k,j,i) = vel(k,j,i) * dStag(k,j,i)
                END DO
            END DO
        END DO

    END SUBROUTINE comp_momentum

    !================================================================

    SUBROUTINE update_velocity(kk, jj, ii, q, vffStag, u, v, w, mom, &
        dt)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, q
        REAL(realk), INTENT(in) :: vffStag(kk, jj, ii)
        REAL(realk), INTENT(inout) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: mom(kk, jj, ii)
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        REAL(realk) :: vel(kk, jj, ii), dStag(kk, jj, ii)
        INTEGER(intk) :: k, j, i

        CALL get_condit_velocity(kk, jj, ii, q, u, v, w, vel)
        CALL comp_material_property_field(kk, jj, ii, vffStag, dStag, rho1, rho2)
    
        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    vel(k,j,i) = mom(k,j,i) / dStag(k,j,i)
                END DO
            END DO
        END DO

        IF (q == 1) THEN
            u = vel
        ELSEIF (q == 2) THEN
            v = vel
        ELSEIF (q == 3) THEN
            w = vel
        ENDIF

    END SUBROUTINE update_velocity

    !================================================================

    SUBROUTINE comp_cWy(kk, jj, ii, vff, cWy)
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
        REAL(realk), INTENT(out) :: cWy(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    IF ( vff(k,j,i) > 0.5_realk ) THEN
                        cWy(k,j,i) = 1.0_realk
                    ELSE
                        cWy(k,j,i) = 0.0_realk
                    END IF
                END DO
            END DO
        END DO

    END SUBROUTINE comp_cWy

    !================================================================

    SUBROUTINE comp_adve(kk, jj, ii, q, l, vff, u, v, w, &
        advr, adve, dx, dy, dz, ddx, ddy, ddz, dt)
    !----------------------------------------------------------------
    !   What it does:
    !   Selects interpolation scheme.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, q, l
        REAL(realk), INTENT(in) :: vff(kk, jj, ii)
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: advr(kk, jj, ii)
        REAL(realk), INTENT(out) :: adve(kk, jj, ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        ! None

        IF ( fluxLimiter == 'QUICK' ) THEN
            CALL comp_adve_quick(kk, jj, ii, q, l, vff, u, v, w, advr, adve, dx, dy, dz, ddx, ddy, ddz, dt)
        ELSEIF ( fluxLimiter == 'ENO' ) THEN
            CALL comp_adve_eno(kk, jj, ii, q, l, u, v, w, advr, adve, dx, dy, dz, ddx, ddy, ddz, dt)
        ELSE
            CALL err_abort(155, "Unknown flux limiter!", __FILE__, __LINE__)
        ENDIF

    END SUBROUTINE comp_adve

    !================================================================

    SUBROUTINE comp_adve_quick(kk, jj, ii, q, l, vff, u, v, w, &
        advr, adve, dx, dy, dz, ddx, ddy, ddz, dt)
    !----------------------------------------------------------------
    !   What it does:
    !   QUICK interpolation to compute the advected veloctiy
    !   (advectee) on staggered grid cells.
    !   adve = advected q (advectee)
    !   advr = advecting q (advector)
    !   An indicator function is used to avoid if-statements within
    !   loops.
    !
    !   Sources: 
    !   T. Arrufat et al., “A mass-momentum consistent, 
    !   Volume-of-Fluid method for incompressible flow on staggered 
    !   grids,” Computers & Fluids, vol. 215, p. 104785, Jan. 2021, 
    !   doi: 10.1016/j.compfluid.2020.104785.
    !
    !   B. P. Leonard, “A stable and accurate convective modelling 
    !   procedure based on quadratic upstream interpolation,”
    !   Computer Methods in Applied Mechanics and Engineering,
    !   vol. 19, no. 1, pp. 59–98, Jun. 1979,
    !   doi: 10.1016/0045-7825(79)90034-3.
    !
    !   W. Aniszewski et al., “PArallel, Robust, Interface Simulator
    !   (PARIS),” Computer Physics Communications, vol. 263,
    !   p. 107849, Jun. 2021, doi: 10.1016/j.cpc.2021.107849.
    !   PARIS source code, grep "interpole_quad" (accessed: Mai 2026)
    !
    !   D. Herrmann, Numerische Mathematik — 40 BASIC-Programme. 
    !   Wiesbaden: Vieweg+Teubner Verlag, 1983. 
    !   doi: 10.1007/978-3-322-96321-5.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, q, l
        REAL(realk), INTENT(in) :: vff(kk,jj,ii)
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: advr(kk, jj, ii)
        REAL(realk), INTENT(out) :: adve(kk, jj, ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: kl, jl, il
        REAL(realk) :: vel(kk,jj,ii)
        LOGICAL :: isIface(kk, jj, ii), isIfaceVic(kk, jj, ii)
        REAL(realk) :: dnx(ii), dny(jj), dnz(kk)
        REAL(realk) :: dcx(ii), dcy(jj), dcz(kk)
        REAL(realk) :: signInd(2), iFacInd(2)
        REAL(realk) :: dnslL, dnslC, dnslR, dcslMi, dcslPl
        REAL(realk) :: adveMi, advePl
        REAL(realk) :: velMi, velCe, velPl, velFP

        CALL get_spatial_indices(l, il, jl, kl)
        CALL get_condit_velocity(kk, jj, ii, q, u, v, w, vel)
        CALL track_iface(isIface, kk, jj, ii, vff)
        CALL track_iface_vic(isIfaceVic, kk, jj, ii, isIface)

        IF ( q == l ) THEN
            dnx(1:ii-1) = ddx(2:ii) ; dnx(ii) = 0.0_realk
            dny(1:jj-1) = ddy(2:jj) ; dny(jj) = 0.0_realk
            dnz(1:kk-1) = ddz(2:kk) ; dnz(kk) = 0.0_realk
            dcx = dnx ; dcy = dny ; dcz = dnz
        ELSE
            dnx = dx  ; dny = dy  ; dnz = dz
            dcx = ddx ; dcy = ddy ; dcz = ddz
        ENDIF

        DO i = 2, ii-2
            DO j = 2, jj-2
                DO k = 2, kk-2
                    signInd(1) = MERGE(1.0_realk, 0.0_realk, advr(k,j,i) >= 0.0_realk)
                    signInd(2) = 1.0_realk - signInd(1)
                    iFacInd(1) = MERGE(1.0_realk, 0.0_realk, isIfaceVic(k,j,i))
                    iFacInd(2) = 1.0_realk - iFacInd(1)

                    dnslL = il*dnx(i-il) + jl*dny(j-jl) + kl*dnz(k-kl)
                    dnslC = il*dnx(i) + jl*dny(j) + kl*dnz(k)
                    dnslR = il*dnx(i+il) + jl*dny(j+jl) + kl*dnz(k+kl)

                    dcslMi = 0.5_realk * (il*dcx(i) + jl*dcy(j) + kl*dcz(k))
                    dcslPl  = dnslC - dcslMi

                    IF ( fluxCentered ) THEN
                        dcslMi = dcslMi - ABS(advr(k,j,i))*dt/2.0_realk
                        dcslPl = dcslPl - ABS(advr(k,j,i))*dt/2.0_realk
                    ENDIF

                    velMi = vel(k-kl,j-jl,i-il)
                    velCe = vel(k,j,i)
                    velPl = vel(k+kl,j+jl,i+il)
                    velFP = vel(k+2*kl,j+2*jl,i+2*il)

                    adveMi = newton_interpolation(velMi, velCe, velPl, dnslL, dnslC, dcslMi)
                    advePl = newton_interpolation(velFP, velPl, velCe, dnslR, dnslC, dcslPl)

                    adve(k,j,i) = iFacInd(1) * (signInd(1)*velCe + signInd(2)*velPl) + &
                                  iFacInd(2) * (signInd(1)*adveMi + signInd(2)*advePl)
                END DO
            END DO
        END DO

    CONTAINS

        !------------------------------------------------------------
        ! Both, Lagrange and Newton form yield the same results.
        ! The Lagrange form was more familiar for me. Hence, I 
        ! implemented it first. The Newton form is closer to 
        ! Leonard's formulation and easier to adjust in the future
        ! (see QUICKEST scheme).
        !------------------------------------------------------------
        ! PURE REAL(realk) FUNCTION lagrange_interpolation(phiUU, phiU, phiD, dsUU, dsD, dsr) RESULT(r)
        !     REAL(realk), INTENT(in) :: phiUU, phiU, phiD, dsUU, dsD, dsr
        !     r = phiUU*(dsr*(dsr-dsD))/(dsUU*(dsUU+dsD)) &
        !         - phiU*((dsr+dsUU)*(dsr-dsD))/(dsUU*dsD) &
        !         + phiD*(dsr*(dsr+dsUU))/(dsD*(dsUU+dsD))
        ! END FUNCTION lagrange_interpolation

        PURE REAL(realk) FUNCTION newton_interpolation(phiUU, phiU, phiD, dsUU, dsD, dsr) RESULT(r)
            REAL(realk), INTENT(in) :: phiUU, phiU, phiD, dsUU, dsD, dsr
            REAL(realk) :: dd1U, dd1D, curv
            dd1D = (phiD - phiU)/dsD
            dd1U = (phiU - phiUU)/dsUU
            curv = (dd1D - dd1U)/(dsUU + dsD)
            r = phiU + dsr*dd1D + dsr*(dsr - dsD)*curv
        END FUNCTION newton_interpolation

    END SUBROUTINE comp_adve_quick

    !================================================================

    SUBROUTINE comp_adve_eno(kk, jj, ii, q, l, u, v, w, &
        advr, adve, dx, dy, dz, ddx, ddy, ddz, dt)
    !----------------------------------------------------------------
    !   What it does:
    !   ENO interpolation to compute the advected veloctiy
    !   (advectee) on staggered grid cells.
    !   adve = advected q (advectee)
    !   advr = advecting q (advector)
    !   An indicator function is used to avoid if-statements within
    !   loops.
    !   
    !   Sources: 
    !   G. Tryggvason, R. Scardovelli, and S. Zaleski, Direct
    !   Numerical Simulations of Gas–Liquid Multiphase Flows,
    !   1st ed. Cambridge University Press, 2011.
    !   doi: 10.1017/CBO9780511975264.
    !
    !   P. K. Sweby, “High Resolution Schemes Using Flux Limiters
    !   for Hyperbolic Conservation Laws,” SIAM J. Numer. Anal.,
    !   vol. 21, no. 5, pp. 995–1011, Oct. 1984,
    !   doi: 10.1137/0721062.
    !
    !   W. Aniszewski et al., “PArallel, Robust, Interface Simulator
    !   (PARIS),” Computer Physics Communications, vol. 263,
    !   p. 107849, Jun. 2021, doi: 10.1016/j.cpc.2021.107849.
    !   PARIS source code, grep "interpole3" (accessed: Mai 2026)
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, q, l
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: advr(kk, jj, ii)
        REAL(realk), INTENT(out) :: adve(kk, jj, ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: dt

        ! Loval variables
        INTEGER(intk) :: k, j, i, kl, jl, il
        REAL(realk) :: vel(kk,jj,ii)
        REAL(realk) :: dnx(ii), dny(jj), dnz(kk)
        REAL(realk) :: dcx(ii), dcy(jj), dcz(kk)
        REAL(realk) :: signInd(2)
        REAL(realk) :: dnslMi, dnslCe, dnslPl, dcslMi, dcslPl
        REAL(realk) :: slopeMi, slopeCe, slopePl, s

        CALL get_spatial_indices(l, il, jl, kl)
        CALL get_condit_velocity(kk, jj, ii, q, u, v, w, vel)

        IF ( q == l ) THEN
            dnx(1:ii-1) = ddx(2:ii) ; dnx(ii) = 0.0_realk
            dny(1:jj-1) = ddy(2:jj) ; dny(jj) = 0.0_realk
            dnz(1:kk-1) = ddz(2:kk) ; dnz(kk) = 0.0_realk
            dcx = dnx ; dcy = dny ; dcz = dnz
        ELSE
            dnx = dx  ; dny = dy  ; dnz = dz
            dcx = ddx ; dcy = ddy ; dcz = ddz
        ENDIF

        DO i = 2, ii-2
            DO j = 2, jj-2
                DO k = 2, kk-2
                    signInd(1) = MERGE(1.0_realk, 0.0_realk, advr(k,j,i) >= 0.0_realk)
                    signInd(2) = 1.0_realk - signInd(1)

                    dnslMi = il*dnx(i-1) + jl*dny(j-1) + kl*dnz(k-1)
                    dnslCe = il*dnx(i) + jl*dny(j) + kl*dnz(k)
                    dnslPl = il*dnx(i+1) + jl*dny(j+1) + kl*dnz(k+1)

                    slopeMi = (vel(k,j,i) - vel(k-kl,j-jl,i-il))/dnslMi
                    slopeCe = (vel(k+kl,j+jl,i+il) - vel(k,j,i))/dnslCe
                    slopePl = (vel(k+2*kl,j+2*jl,i+2*il) - vel(k+kl,j+jl,i+il))/dnslPl

                    s = signInd(1) * minmod(slopeMi, slopeCe) + &
                        signInd(2) * minmod(slopeCe, slopePl)

                    dcslMi = 0.5_realk * (il*dcx(i) + jl*dcy(j) + kl*dcz(k))
                    dcslPl = dnslCe - dcslMi

                    IF ( fluxCentered ) THEN
                        dcslMi = dcslMi - ABS(advr(k,j,i))*dt/2.0_realk
                        dcslPl = dcslPl - ABS(advr(k,j,i))*dt/2.0_realk
                    ENDIF

                    adve(k,j,i) = signInd(1) * (vel(k,j,i) + s*dcslMi) + &
                                  signInd(2) * (vel(k+kl,j+jl,i+il) - s*dcslPl)
                END DO
            END DO
        END DO

    CONTAINS

        PURE REAL(realk) FUNCTION minmod(a, b) RESULT(res)
            REAL(realk), INTENT(in) :: a, b
            res = 0.5_realk * ( SIGN(1.0_realk, a) + SIGN(1.0_realk, b) ) * MIN(ABS(a), ABS(b))
        END FUNCTION minmod

    END SUBROUTINE comp_adve_eno

    !================================================================

    SUBROUTINE comp_advr_linear_interpolation(kk, jj, ii, q, l, &
        u, v, w, advr, ddx, ddy, ddz)
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
        INTEGER(intk), INTENT(in) :: kk, jj, ii, q, l
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(out) :: advr(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)

        ! Loval variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: kq, jq, iq
        REAL(realk) :: vel(kk,jj,ii)
        REAL(realk) :: ddnqMi, ddnqPl

        CALL get_spatial_indices(q, iq, jq, kq)
        CALL get_condit_velocity(kk, jj, ii, l, u, v, w, vel)

        IF ( q == l ) THEN
            DO i = 2, ii-2
                DO j = 2, jj-2
                    DO k = 2, kk-2
                        advr(k,j,i) = 0.5_realk*(vel(k,j,i) + vel(k+kq,j+jq,i+iq))
                    ENDDO
                ENDDO
            ENDDO
        ELSE
            DO i = 2, ii-2
                DO j = 2, jj-2
                    DO k = 2, kk-2
                        ddnqMi = iq*ddx(i) + jq*ddy(j) + kq*ddz(k)
                        ddnqPl = iq*ddx(i+iq) + jq*ddy(j+jq) + kq*ddz(k+kq)
                        advr(k,j,i) = (vel(k,j,i)*ddnqMi + vel(k+kq,j+jq,i+iq)*ddnqPl)/(ddnqMi + ddnqPl)
                    ENDDO
                ENDDO
            ENDDO
        ENDIF

    END SUBROUTINE comp_advr_linear_interpolation

END MODULE multiphase_vof_transport_mod
