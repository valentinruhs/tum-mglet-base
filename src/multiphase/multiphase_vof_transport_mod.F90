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
    USE multiphase_plic_mod, ONLY: comp_frac, iface_reconstruction, comp_stag_frac, track_iface, track_iface_vic
    USE rungekutta_mod, ONLY: rk_2n_t
    USE multiphasecore_mod, ONLY: gmol1, gmol2, rho1, rho2, grav, permutation_multiphase
    USE multiphase_material_mod, ONLY: comp_material_property_field, comp_property_face_value
    USE flowcore_mod, ONLY: gradp
    USE connect2_mod, ONLY: connect
    USE parent_mod, ONLY: parent
    USE grids_mod, ONLY: minlevel, maxlevel
    USE multiphase_io_mod, ONLY: apprVol
    USE err_mod, ONLY: errr
    USE multiphase_utils_mod, ONLY: get_spatial_indices, get_spatial_extents, get_condit_velocity
    
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

    SUBROUTINE comp_flux_cent(kk, jj, ii, l, vff, isIface, u, v, w, alpha, dt, normx, normy, normz, ddx, ddy, ddz, tol, vffFlux)
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
        LOGICAL, INTENT(in) :: isIface(kk, jj, ii)
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: alpha(kk, jj, ii)
        REAL(realk), INTENT(in) :: dt
        REAL(realk), INTENT(in) :: normx(kk, jj, ii), normy(kk, jj, ii), normz(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: tol
        REAL(realk), INTENT(out) :: vffFlux(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: kl, jl, il
        REAL(realk) :: vel(kk, jj, ii), dds, norms
        REAL(realk) :: dimx, dimy, dimz
        REAL(realk) :: fluxedProp, fluxWidth, fluxAlpha

        CALL get_spatial_indices(kk, jj, ii, l, il, jl, kl)
        CALL get_condit_velocity(kk, jj, ii, l, u, v, w, vel)

        DO i = 2, ii-2
            DO j = 2, jj-2
                DO k = 2, kk-2
                    IF ( vel(k,j,i) > tol ) THEN
                        IF ( isIface(k,j,i) ) THEN
                            ! Compute characteristic length and norm
                            dds = il * ddx(i) + jl * ddy(j) + kl * ddz(k)
                            norms = il * normx(k,j,i) + jl * normy(k,j,i) + kl * normz(k,j,i)

                            ! Compute face fluxwidth and proper alpha
                            fluxWidth = abs( vel(k,j,i) ) * dt
                            fluxAlpha = alpha(k,j,i) - norms * ( dds - fluxWidth )

                            ! Compute dimensions of fluxed cuboid
                            dimx = il * fluxWidth + jl * ddy(j) + kl * ddz(k)
                            dimy = il * ddx(i) + jl * fluxWidth + kl * ddz(k)
                            dimz = il * ddx(i) + jl * ddy(j) + kl * fluxWidth

                            ! Compute vff in fluxed cuboid
                            CALL comp_frac(fluxedProp, fluxAlpha, vff(k,j,i), dimx, dimy, dimz, normx(k,j,i), normy(k,j,i), normz(k,j,i), tol)
                        ELSE
                            ! Compute characteristic length
                            dds = il * ddx(i) + jl * ddy(j) + kl * ddz(k)

                            ! Compute face fluxwidth
                            fluxWidth = abs( vel(k,j,i) ) * dt
                            fluxedProp = vff(k,j,i)
                        END IF
                    ELSE IF ( vel(k,j,i) < -tol ) THEN
                        IF ( isIface(k+kl,j+jl,i+il) ) THEN
                            ! Compute characteristic length
                            dds = il * ddx(i+1) + jl * ddy(j+1) + kl * ddz(k+1)

                            ! Compute face fluxwidth and proper alpha
                            fluxWidth = abs( vel(k,j,i) ) * dt
                            fluxAlpha = alpha(k+kl,j+jl,i+il)

                            ! Compute dimensions of fluxed cuboid
                            dimx = il * fluxWidth + jl * ddy(j+1) + kl * ddz(k+1)
                            dimy = il * ddx(i+1) + jl * fluxWidth + kl * ddz(k+1)
                            dimz = il * ddx(i+1) + jl * ddy(j+1) + kl * fluxWidth

                            ! Compute vff in fluxed cuboid
                            CALL comp_frac(fluxedProp, fluxAlpha, vff(k+kl,j+jl,i+il), dimx, dimy, dimz, normx(k+kl,j+jl,i+il), normy(k+kl,j+jl,i+il), normz(k+kl,j+jl,i+il), tol)
                        ELSE
                            ! Compute characteristic length
                            dds = il * ddx(i+1) + jl * ddy(j+1) + kl * ddz(k+1)

                            ! Compute face fluxwidth
                            fluxWidth = abs( vel(k,j,i) ) * dt
                            fluxedProp = vff(k+kl,j+jl,i+il)
                        END IF
                    ELSE
                        fluxedProp = 0.0_realk
                    END IF
                    vffFlux(k,j,i) = vel(k,j,i) * fluxedProp
                END DO
            END DO
        END DO

    END SUBROUTINE comp_flux_cent

    !================================================================

    SUBROUTINE comp_flux_stag(kk, jj, ii, q, l, vff, isIface, advr, alpha, dt, normx, normy, normz, dx, dy, dz, ddx, ddy, ddz, tol, vffFlux, complvffFlux)
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
        LOGICAL, INTENT(in) :: isIface(kk, jj, ii)
        REAL(realk), INTENT(in) :: advr(kk, jj, ii)
        REAL(realk), INTENT(in) :: alpha(kk, jj, ii)
        REAL(realk), INTENT(in) :: dt
        REAL(realk), INTENT(in) :: normx(kk, jj, ii), normy(kk, jj, ii), normz(kk, jj, ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: tol
        REAL(realk), INTENT(out) :: vffFlux(kk, jj, ii), complvffFlux(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i, kl, jl, il, kq, jq, iq
        INTEGER(intk) :: kDonMi, jDonMi, iDonMi, kDonPl, jDonPl, iDonPl
        REAL(realk) :: ddsDon, normDon, fluxWidth, fluxAlpha, fluxDimx, fluxDimy, fluxDimz, fluxedProp
        REAL(realk) :: farEnd, ddslMi, ddslPl, ddsqMi, ddsqPl, normlMi, normlPl, normqMi, normqPl, fluxAlphaMi, fluxAlphaPl, fluxDimxMi, fluxDimyMi, fluxDimzMi, fluxDimxPl, fluxDimyPl, fluxDimzPl, fracMi, fracPl


        CALL get_spatial_indices(kk, jj, ii, l, il, jl, kl)
        CALL get_spatial_indices(kk, jj, ii, q, iq, jq, kq)

        IF ( q == l ) THEN
            DO i = 2, ii-2
                DO j = 2, jj-2
                    DO k = 2, kk-2
                        IF ( ABS(advr(k,j,i)) < tol ) THEN
                            vffFlux(k,j,i)      = 0.0_realk
                            complVffFlux(k,j,i) = 0.0_realk
                            CYCLE
                        ENDIF

                        ddsDon  = il*ddx(i+il) + jl*ddy(j+jl) + kl*ddz(k+kl)
                        normDon = il*normx(k+kl,j+jl,i+il) &
                                + jl*normy(k+kl,j+jl,i+il) &
                                + kl*normz(k+kl,j+jl,i+il)

                        fluxWidth = ABS( advr(k,j,i) ) * dt
                        IF ( fluxWidth > 0.5_realk*ddsDon ) CALL errr(__FILE__, __LINE__)

                        IF ( isIface(k+kl,j+jl,i+il) ) THEN
                            fluxAlpha = alpha(k+kl,j+jl,i+il) - normDon * &
                                ( 0.5_realk*ddsDon &
                                - MERGE(fluxWidth, 0.0_realk, advr(k,j,i) > 0.0_realk) )
                            fluxDimx = il*fluxWidth + (1-il)*ddx(i)
                            fluxDimy = jl*fluxWidth + (1-jl)*ddy(j)
                            fluxDimz = kl*fluxWidth + (1-kl)*ddz(k)
                            CALL comp_frac(fluxedProp, fluxAlpha, vff(k+kl,j+jl,i+il), &
                                fluxDimx, fluxDimy, fluxDimz, &
                                normx(k+kl,j+jl,i+il), normy(k+kl,j+jl,i+il), &
                                normz(k+kl,j+jl,i+il), tol)
                        ELSE
                            fluxedProp = vff(k+kl,j+jl,i+il)
                        ENDIF

                        vffFlux(k,j,i) = advr(k,j,i) * fluxedProp
                        complVffFlux(k,j,i) = advr(k,j,i) * ( 1.0_realk - fluxedProp )
                    END DO
                END DO
            END DO
        ELSE
            DO i = 2, ii-2
                DO j = 2, jj-2
                    DO k = 2, kk-2

                    IF ( ABS(advr(k,j,i)) < tol ) THEN
                        vffFlux(k,j,i)      = 0.0_realk
                        complVffFlux(k,j,i) = 0.0_realk
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
                    IF ( fluxWidth > ddslMi .OR. fluxWidth > ddslPl ) CALL errr(__FILE__, __LINE__)

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

                        CALL comp_frac(fracMi, fluxAlphaMi, vff(kDonMi,jDonMi,iDonMi), &
                            fluxDimxMi, fluxDimyMi, fluxDimzMi, &
                            normx(kDonMi,jDonMi,iDonMi), normy(kDonMi,jDonMi,iDonMi), normz(kDonMi,jDonMi,iDonMi), tol)
                    ELSE
                        fracMi = vff(kDonMi,jDonMi,iDonMi)
                    ENDIF

                    IF ( isIface(kDonPl,jDonPl,iDonPl) ) THEN
                        fluxAlphaPl = alpha(kDonPl,jDonPl,iDonPl) &
                                    - normlPl * farEnd * ( ddslPl - fluxWidth )

                        fluxDimxPl = il*fluxWidth + iq*ddx(iDonPl)/2.0_realk + (1-il-iq)*ddx(iDonPl)
                        fluxDimyPl = jl*fluxWidth + jq*ddy(jDonPl)/2.0_realk + (1-jl-jq)*ddy(jDonPl)
                        fluxDimzPl = kl*fluxWidth + kq*ddz(kDonPl)/2.0_realk + (1-kl-kq)*ddz(kDonPl)

                        CALL comp_frac(fracPl, fluxAlphaPl, vff(kDonPl,jDonPl,iDonPl), &
                            fluxDimxPl, fluxDimyPl, fluxDimzPl, &
                            normx(kDonPl,jDonPl,iDonPl), normy(kDonPl,jDonPl,iDonPl), &
                            normz(kDonPl,jDonPl,iDonPl), tol)
                    ELSE
                        fracPl = vff(kDonPl,jDonPl,iDonPl)
                    ENDIF

                    fluxedProp = ( fracMi*ddsqMi + fracPl*ddsqPl ) / ( ddsqMi + ddsqPl )

                    vffFlux(k,j,i) = advr(k,j,i) * fluxedProp
                    complVffFlux(k,j,i) = advr(k,j,i) * ( 1.0_realk - fluxedProp )

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

        CALL check_continuity(tol, itstep)
        CALL check_solenoidality(tol, itstep)

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

            CALL adve_operator(kk, jj, ii, u, v, w, vff, d, dx, dy, dz, ddx, ddy, ddz, dt, itstep, tol, &
                normx, normy, normz, alpha)

            CALL diff_operator(kk, jj, ii, u, v, w, vffPrev, ddx, ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz, tol, uo, vo, wo)

            CALL pres_operator(kk, jj, ii, vff, p, rdx, rdy, rdz, igrid, uo, vo, wo)

            CALL exte_operator(kk, jj, ii, vff, dx, dy, dz, ddx, ddy, ddz, uo, vo, wo)

        END DO

    END SUBROUTINE multiphase_solve

    !================================================================

    SUBROUTINE adve_operator(kk, jj, ii, u, v, w, vff, d, dx, dy, dz, ddx, ddy, ddz, dt, itstep, tol, &
        normx, normy, normz, alpha)
    !----------------------------------------------------------------
    !   What it does:
    !    
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(inout) :: vff(kk, jj, ii), d(kk, jj, ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: tol, dt
        INTEGER(intk), INTENT(in) :: itstep
        REAL(realk), INTENT(out) :: normx(kk, jj, ii), normy(kk, jj, ii), normz(kk, jj, ii), alpha(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: q, advSeq(3), l, splitDir
        REAL(realk) :: vffStag(kk, jj, ii)
        REAL(realk) :: dStag(kk, jj, ii)
        LOGICAL :: isIface(kk, jj, ii)
        LOGICAL :: isIfaceVic(kk, jj, ii)
        REAL(realk) :: vffFlux(kk, jj, ii), complVffFlux(kk, jj, ii), mom(kk, jj, ii), cWY(kk, jj, ii), cWYStag(kk, jj, ii)
        REAL(realk) :: advr(kk, jj, ii)
        REAL(realk) :: adve(kk, jj, ii)
        REAL(realk) :: vel(kk, jj, ii), velo(kk, jj, ii)
        REAL(realk) :: mom4D(kk, jj, ii, 3), vffStag4D(kk, jj, ii, 3), cWYStag4D(kk, jj, ii, 3)
        REAL(realk) :: uPrev(kk, jj, ii), vPrev(kk, jj, ii), wPrev(kk, jj, ii)

        uPrev = u
        vPrev = v
        wPrev = w

        CALL def_advection_sequence(itstep, advSeq)
        CALL iface_reconstruction(kk, jj, ii, vff, ddx, ddy, ddz, tol, normx, normy, normz, alpha, isIface, isIfaceVic)

        DO q = 1, 3

            CALL comp_stag_frac(kk, jj, ii, q, vff, alpha, isIface, ddx, ddy, ddz, normx, normy, normz, tol, vffStag)
            CALL comp_material_property_field(kk, jj, ii, vffStag, rho1, rho2, 'ARI', dStag)
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

                CALL comp_advr_linear_interpolation(kk, jj, ii, q, l, u, v, w, advr)
                CALL comp_adve_quick(kk, jj, ii, q, l, u, v, w, advr, vffStag, tol, adve)
                ! CALL comp_adve_eno(kk, jj, ii, q, l, u, v, w, advr, dx, dy, dz, ddx, ddy, ddz, tol, adve)
                CALL comp_flux_stag(kk, jj, ii, q, l, vff, isIface, advr, alpha, dt, normx, normy, normz, dx, dy, dz, ddx, ddy, ddz, tol, vffFlux, complVffFlux)
                CALL adv_mom(kk, jj, ii, q, l, u, v, w, advr, adve, vffStag, vffFlux, complVffFlux, cWYStag, dx, dy, dz, ddx, ddy, ddz, dt, mom)
                CALL adv_vof(kk, jj, ii, q, l, vffFlux, cWYStag, advr, dx, dy, dz, ddx, ddy, ddz, dt, tol, vffStag)
                CALL comp_velocity_change(kk, jj, ii, q, u, v, w, vffStag, mom, dt, velo)
                CALL update_velocity(kk, jj, ii, q, u, v, w, velo, dt)

                mom4D(:,:,:,q) = mom
                vffStag4D(:,:,:,q) = vffStag

            END DO

            CALL iface_reconstruction(kk, jj, ii, vff, ddx, ddy, ddz, tol, normx, normy, normz, alpha, isIface, isIfaceVic)
            CALL comp_flux_cent(kk, jj, ii, l, vff, isIface, uPrev, vPrev, wPrev, alpha, dt, normx, normy, normz, ddx, ddy, ddz, tol, vffFlux)
            CALL get_condit_velocity(kk, jj, ii, l, uPrev, vPrev, wPrev, vel)
            CALL adv_vof(kk, jj, ii, 0_intk, l, vffFlux, cWY, vel, dx, dy, dz, ddx, ddy, ddz, dt, tol, vff)
            CALL clip_vff(kk, jj, ii, tol, vff)
            CALL app_bcon()
            CALL comp_material_property_field(kk, jj, ii, vff, rho1, rho2, 'ARI', d)

        END DO

    END SUBROUTINE adve_operator

    !================================================================

    SUBROUTINE diff_operator(kk, jj, ii, u, v, w, vff, ddx, ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz, tol, uo, vo, wo)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------
    
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii), vff(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: rdx(ii), rdy(jj), rdz(kk)
        REAL(realk), INTENT(in) :: rddx(ii), rddy(jj), rddz(kk)
        REAL(realk), INTENT(in) :: tol
        REAL(realk), INTENT(inout) :: uo(kk, jj, ii), vo(kk, jj, ii), wo(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i, q
        REAL(realk) :: d(kk, jj, ii), g(kk, jj, ii)
        REAL(realk) :: normx(kk, jj, ii), normy(kk, jj, ii), normz(kk, jj, ii), alpha(kk, jj, ii)
        LOGICAL :: isIface(kk, jj, ii)
        REAL(realk) :: vffStag(kk, jj, ii)
        REAL(realk) :: ge(kk, jj, ii, 3), gn(kk, jj, ii, 3), gt(kk, jj, ii, 3)
        REAL(realk) :: tauxxe, tauxxw, tauxyn, tauxys, tauxzt, tauxzb
        REAL(realk) :: tauyxe, tauyxw, tauyyn, tauyys, tauyzt, tauyzb
        REAL(realk) :: tauzxe, tauzxw, tauzyn, tauzys, tauzzt, tauzzb

        CALL comp_material_property_field(kk, jj, ii, vff, rho1, rho2, 'ARI', d)
        CALL comp_material_property_field(kk, jj, ii, vff, gmol1, gmol2, 'HAR', g)
        CALL iface_reconstruction(kk, jj, ii, vff, ddx, ddy, ddz, tol, normx, normy, normz, alpha, isIface)

        DO q = 1, 3

            CALL comp_stag_frac(kk, jj, ii, q, vff, alpha, isIface, ddx, ddy, ddz, normx, normy, normz, tol, vffStag)
            CALL comp_mean_harm(kk, jj, ii, q, vffStag, g, ge(:,:,:,q), gn(:,:,:,q), gt(:,:,:,q))
            ! CALL comp_mean_arit(kk, jj, ii, q, vffStag, g, ge, gn, gt)

        ENDDO

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    ! Normal stresses
                    tauxxe = ge(k,j,i,1) * 2.0_realk * (u(k,j,i+1) - u(k,j,i)) * rddx(i+1)
                    tauxxw = ge(k,j,i-1,1) * 2.0_realk * (u(k,j,i) - u(k,j,i-1)) * rddx(i)

                    ! Shear stresses
                    tauxyn = gn(k,j,i,1) * ( (u(k,j+1,i) - u(k,j,i)) * rdy(j)   + (v(k,j,i+1) - v(k,j,i))     * rdx(i) )
                    tauxys = gn(k,j-1,i,1) * ( (u(k,j,i) - u(k,j-1,i)) * rdy(j-1) + (v(k,j-1,i+1) - v(k,j-1,i)) * rdx(i) )
                    tauxzt = gt(k,j,i,1) * ( (u(k+1,j,i) - u(k,j,i)) * rdz(k)   + (w(k,j,i+1) - w(k,j,i))     * rdx(i) )
                    tauxzb = gt(k-1,j,i,1) * ( (u(k,j,i) - u(k-1,j,i)) * rdz(k-1) + (w(k-1,j,i+1) - w(k-1,j,i)) * rdx(i) )

                    ! Change due to diffusion
                    uo(k,j,i) = uo(k,j,i) + 2.0_realk/( d(k,j,i) + d(k,j,i+1) ) * &
                        ( ( tauxxe - tauxxw ) * rdx(i) + &
                          ( tauxyn - tauxys ) * rddy(j) + &
                          ( tauxzt - tauxzb ) * rddz(k) )
                END DO
            END DO
        END DO

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    ! Shear stresses
                    tauyxe = ge(k,j,i,2) * ( (u(k,j+1,i) - u(k,j,i))     * rdy(j) + (v(k,j,i+1) - v(k,j,i)) * rdx(i)   )
                    tauyxw = ge(k,j,i-1,2) * ( (u(k,j+1,i-1) - u(k,j,i-1)) * rdy(j) + (v(k,j,i) - v(k,j,i-1)) * rdx(i-1) )

                    ! Normal stresses
                    tauyyn = gn(k,j,i,2) * 2.0_realk * (v(k,j+1,i) - v(k,j,i)) * rddy(j+1)
                    tauyys = gn(k,j-1,i,2) * 2.0_realk * (v(k,j,i) - v(k,j-1,i)) * rddy(j)
                    
                    ! Shear stresses
                    tauyzt = gt(k,j,i,2) * ( (v(k+1,j,i) - v(k,j,i)) * rdz(k)   + (w(k,j+1,i) - w(k,j,i))     * rdy(j) )
                    tauyzb = gt(k-1,j,i,2) * ( (v(k,j,i) - v(k-1,j,i)) * rdz(k-1) + (w(k-1,j+1,i) - w(k-1,j,i)) * rdy(j) )

                    ! Change due to diffusion
                    vo(k,j,i) = vo(k,j,i) + 2.0_realk/( d(k,j,i) + d(k,j+1,i) ) * &
                        ( ( tauyxe - tauyxw ) * rddx(i) + &
                          ( tauyyn - tauyys ) * rdy(j) + &
                          ( tauyzt - tauyzb ) * rddz(k) )
                END DO
            END DO
        END DO

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    ! Shear stresses
                    tauzxe = ge(k,j,i,3) * ( (u(k+1,j,i) - u(k,j,i))     * rdz(k) + (w(k,j,i+1) - w(k,j,i)) * rdx(i)   )
                    tauzxw = ge(k,j,i-1,3) * ( (u(k+1,j,i-1) - u(k,j,i-1)) * rdz(k) + (w(k,j,i) - w(k,j,i-1)) * rdx(i-1) )
                    tauzyn = gn(k,j,i,3) * ( (v(k+1,j,i) - v(k,j,i))     * rdz(k) + (w(k,j+1,i) - w(k,j,i)) * rdy(j)   )
                    tauzys = gn(k,j-1,i,3) * ( (v(k+1,j-1,i) - v(k,j-1,i)) * rdz(k) + (w(k,j,i) - w(k,j-1,i)) * rdy(j-1) )
                    
                    ! Normal stresses
                    tauzzt = gt(k,j,i,3) * 2.0_realk * (w(k+1,j,i) - w(k,j,i)) * rddz(k+1)
                    tauzzb = gt(k-1,j,i,3) * 2.0_realk * (w(k,j,i) - w(k-1,j,i)) * rddz(k)

                    ! Change due to diffusion
                    wo(k,j,i) = wo(k,j,i) + 2.0_realk/( d(k,j,i) + d(k+1,j,i) ) * &
                        ( ( tauzxe - tauzxw ) * rddx(i) + &
                          ( tauzyn - tauzys ) * rddy(j) + &
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
        REAL(realk) :: rhoe(kk, jj, ii), rhon(kk, jj, ii), rhot(kk, jj, ii)
        INTEGER(intk) :: gradpflag
        REAL(realk) :: gpx, gpy, gpz
        INTEGER(intk) :: i, j, k

        CALL comp_property_face_value(kk, jj, ii, vff, rho1, rho2, 'ARI', rhoe, rhon, rhot)

        CALL get_gradpxflag(gradpflag, igrid)
        gpx = gradp(1)*gradpflag
        gpy = gradp(2)*gradpflag
        gpz = gradp(3)*gradpflag

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    uo(k,j,i) = uo(k,j,i) - 1.0_realk / rhoe(k,j,i) * ( ( p(k,j,i+1) - p(k,j,i) ) * rdx(i) + gpx )
                END DO
            END DO
        END DO

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    vo(k,j,i) = vo(k,j,i) - 1.0_realk / rhon(k,j,i) * ( ( p(k,j+1,i) - p(k,j,i) ) * rdy(j) + gpy )
                END DO
            END DO
        END DO

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    wo(k,j,i) = wo(k,j,i) - 1.0_realk / rhot(k,j,i) * ( ( p(k+1,j,i) - p(k,j,i) ) * rdz(k) + gpz )
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
        CALL comp_material_property_field(kk, jj, ii, vff, rho1, rho2, 'ARI', d)

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    uo(k,j,i) = uo(k,j,i) - grav(1) ! * 0.5 * ( d(k,j,i) + d(k,j,i+1) )
                    vo(k,j,i) = vo(k,j,i) - grav(2) ! * 0.5 * ( d(k,j,i) + d(k,j+1,i) )
                    wo(k,j,i) = wo(k,j,i) - grav(3) ! * 0.5 * ( d(k,j,i) + d(k+1,j,i) )
                ENDDO
            ENDDO
        ENDDO

    END SUBROUTINE exte_operator

    !================================================================

    SUBROUTINE adv_vof(kk, jj, ii, q, l, vffFlux, cWY, vel, dx, dy, dz, ddx, ddy, ddz, dt, tol, vff)
    !----------------------------------------------------------------
    !   What it does:
    !    
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, q, l
        REAL(realk), INTENT(in) :: vffFlux(kk, jj, ii)
        REAL(realk), INTENT(in) :: cWY(kk, jj, ii), vel(kk,jj,ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: dt, tol
        REAL(realk), INTENT(inout) :: vff(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: i, j, k
        INTEGER(intk) :: il, jl, kl
        REAL(realk) :: dsx(ii), dsy(jj), dsz(kk), dsCV
        REAL(realk) :: div(kk, jj, ii)

        CALL get_spatial_indices(kk, jj, ii, l, il, jl, kl)
        CALL get_spatial_extents(kk, jj, ii, q, l, dx, dy, dz, ddx, ddy, ddz, dsx, dsy, dsz)

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    dsCV = il * dsx(i) + jl * dsy(j) + kl * dsz(k)
                    div(k,j,i) = ( vel(k,j,i) - vel(k-kl,j-jl,i-il) ) / dsCV
                    vff(k,j,i) = vff(k,j,i) - dt/dsCV * ( vffFlux(k,j,i) - vffFlux(k-kl,j-jl,i-il) ) + dt * cWY(k,j,i) * div(k,j,i)
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
        REAL(realk) :: dsx(ii), dsy(jj), dsz(kk), dsCV
        REAL(realk) :: vel(kk,jj,ii)
        REAL(realk) :: dStag(kk, jj, ii)
        REAL(realk) :: momFlux(kk, jj, ii)
        REAL(realk) :: div, com

        CALL get_spatial_indices(kk, jj, ii, l, il, jl, kl)
        CALL get_spatial_extents(kk, jj, ii, q, l, dx, dy, dz, ddx, ddy, ddz, dsx, dsy, dsz)
        CALL get_condit_velocity(kk, jj, ii, q, u, v, w, vel)

        CALL comp_material_property_field(kk, jj, ii, vff, rho1, rho2, 'ARI', dStag)

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
                    dsCV = il * dsx(i) + jl * dsy(j) + kl * dsz(k)
                    div = ( advr(k,j,i) - advr(k-kl,j-jl,i-il) ) / dsCV
                    com = ( rho1 * cWY(k,j,i) + rho2 * (1.0_realk - cWY(k,j,i)) ) * div
                    mom(k,j,i) = mom(k,j,i) - dt/dsCV * ( momFlux(k,j,i) - momFlux(k-kl,j-jl,i-il) ) + dt * vel(k,j,i) * com
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
        REAL(realk), INTENT(out) :: velo(kk,jj,ii)

        ! Local variables
        INTEGER(intk) :: i, j, k
        REAL(realk) :: vel(kk,jj,ii)
        REAL(realk) :: dStag(kk, jj, ii)

        CALL get_condit_velocity(kk, jj, ii, q, u, v, w, vel)
        CALL comp_material_property_field(kk, jj, ii, vff, rho1, rho2, 'ARI', dStag)
    
        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    IF ( vff(k,j,i) >= 0.0_realk .AND. vff(k,j,i) <= 1.0_realk ) THEN
                        velo(k,j,i) = ( mom(k,j,i) / dStag(k,j,i) - vel(k,j,i) ) / dt
                    END IF
                END DO
            END DO
        END DO

    END SUBROUTINE comp_velocity_change

    !================================================================

    SUBROUTINE update_velocity(kk, jj, ii, q, u, v, w, velo, dt)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, q
        REAL(realk), INTENT(inout) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: velo(kk, jj, ii)
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        ! None

        IF (q == 1) THEN
            u = u + velo * dt
        ELSEIF (q == 2) THEN
            v = v + velo * dt
        ELSEIF (q == 3) THEN
            w = w + velo * dt
        ENDIF

    END SUBROUTINE

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

    SUBROUTINE check_solenoidality(tol, itstep)
    !----------------------------------------------------------------
    !   What it does:
    !   Checks, if the field is solenoidal. Only used during coding.
    !   Solenoidality is assured by pressure correction.
    !----------------------------------------------------------------

        ! Subroutine arguments
        REAL(realk) :: tol
        INTEGER(intk) :: itstep

        ! Local variables
        TYPE(field_t), POINTER :: u_f, v_f, w_f
        TYPE(field_t), POINTER :: ddx_f, ddy_f, ddz_f
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        REAL(realk), POINTER, CONTIGUOUS :: u(:,:,:), v(:,:,:), w(:,:,:)
        INTEGER(intk) :: kk, jj, ii, k, j, i, n, igrid
        REAL(realk), ALLOCATABLE :: div(:,:,:)

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
            div = 0.0_realk

            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        div(k,j,i) = ( u(k,j,i) - u(k,j,i-1) )/ddx(i) + &
                                     ( v(k,j,i) - v(k,j-1,i) )/ddy(j) + &
                                     ( w(k,j,i) - w(k-1,j,i) )/ddz(k)
                    END DO
                END DO
            END DO

            IF ( SUM(ABS(div)) >= tol ) THEN
                WRITE(*,*) "Solenoidality violated at itstep ", itstep, ": L1-Norm Error = ", SUM(ABS(div))
            ENDIF

        ENDDO

    END SUBROUTINE check_solenoidality

    !================================================================

    SUBROUTINE check_continuity(tol, itstep)
    !----------------------------------------------------------------
    !   What it does:
    !    
    !----------------------------------------------------------------

        ! Subroutine arguments
        REAL(realk) :: tol
        INTEGER(intk) :: itstep

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
            WRITE(*,*) "Conti. equ. violated at itstep ", itstep, ": vol. err. [%] = ", ABS(apprVol - volFl1) / apprVol * 100.0_realk
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

    SUBROUTINE comp_mean_harm(kk, jj, ii, q, vff, g, ge, gn, gt)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !   
    !   Source: 
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, q
        REAL(realk), INTENT(in) :: vff(kk, jj, ii), g(kk, jj, ii)
        REAL(realk), INTENT(out) :: ge(kk, jj, ii), gn(kk, jj, ii), gt(kk, jj, ii)

        ! Loval variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: kq, jq, iq

        CALL get_spatial_indices(kk, jj, ii, q, iq, jq, kq)

        DO i = 2, ii-2
            DO j = 2, jj-2
                DO k = 2, kk-2

                    ge(k,j,i) = iq * g(k,j,i+1) + &
                                jq * 1.0_realk / ( vff(k,j,i+1) / gmol1 + ( 1.0_realk - vff(k,j,i+1) ) / gmol2 ) + &
                                kq * 1.0_realk / ( vff(k,j,i+1) / gmol1 + ( 1.0_realk - vff(k,j,i+1) ) / gmol2 )
                    gn(k,j,i) = iq * 1.0_realk / ( vff(k,j+1,i) / gmol1 + ( 1.0_realk - vff(k,j+1,i) ) / gmol2 ) + &
                                jq * g(k,j+1,i) + &
                                kq * 1.0_realk / ( vff(k,j+1,i) / gmol1 + ( 1.0_realk - vff(k,j+1,i) ) / gmol2 )
                    gt(k,j,i) = iq * 1.0_realk / ( vff(k+1,j,i) / gmol1 + ( 1.0_realk - vff(k+1,j,i) ) / gmol2 ) + &
                                jq * 1.0_realk / ( vff(k+1,j,i) / gmol1 + ( 1.0_realk - vff(k+1,j,i) ) / gmol2 ) + &
                                kq * g(k+1,j,i)

                ENDDO
            ENDDO
        ENDDO

    END SUBROUTINE comp_mean_harm

    !================================================================

    ! SUBROUTINE comp_mean_arit(kk, jj, ii, q, vff, g, ge, gn, gt)
    ! !----------------------------------------------------------------
    ! !   What it does:
    ! !   
    ! !   
    ! !   Source: 
    ! !   
    ! !----------------------------------------------------------------

    !     ! Subroutine arguments
    !     INTEGER(intk), INTENT(in) :: kk, jj, ii, q
    !     REAL(realk), INTENT(in) :: vff(kk, jj, ii), g(kk, jj, ii)
    !     REAL(realk), INTENT(out) :: ge(kk, jj, ii), gn(kk, jj, ii), gt(kk, jj, ii)

    !     ! Loval variables
    !     INTEGER(intk) :: k, j, i
    !     INTEGER(intk) :: kq, jq, iq

    !     CALL get_spatial_indices(kk, jj, ii, q, iq, jq, kq)

    !     DO i = 2, ii-2
    !         DO j = 2, jj-2
    !             DO k = 2, kk-2

    !                 ge(k,j,i) = fx * g(k,j,i+1) + &
    !                             jl * ( vff(k,j,i+1) * gmol1 + ( 1.0_realk - vff(k,j,i+1) ) * gmol2 ) + &
    !                             kl * ( vff(k,j,i+1) * gmol1 + ( 1.0_realk - vff(k,j,i+1) ) * gmol2 )
    !                 gn(k,j,i) = fx * ( vff(k,j+1,i) * gmol1 + ( 1.0_realk - vff(k,j+1,i) ) * gmol2 ) + &
    !                             jl * g(k,j+1,i) + &
    !                             kl * ( vff(k,j+1,i) * gmol1 + ( 1.0_realk - vff(k,j+1,i) ) * gmol2 )
    !                 gt(k,j,i) = fx * ( vff(k+1,j,i) * gmol1 + ( 1.0_realk - vff(k+1,j,i) ) * gmol2 ) + &
    !                             jl * ( vff(k+1,j,i) * gmol1 + ( 1.0_realk - vff(k+1,j,i) ) * gmol2 ) + &
    !                             kl * g(k+1,j,i)

    !             ENDDO
    !         ENDDO
    !     ENDDO

    ! END SUBROUTINE comp_mean_arit

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
                    signInd(1) = MERGE(1.0_realk, 0.0_realk, advr(k,j,i) >= 0.0_realk)
                    signInd(2) = 1.0_realk - signInd(1)
                    iFacInd(1) = MERGE(1.0_realk, 0.0_realk, isIfaceVic(k,j,i))
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

    ! SUBROUTINE comp_adve_eno(kk, jj, ii, q, l, u, v, w, &
    !     advr, dx, dy, dz, ddx, ddy, ddz, tol, adve)
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
    !     INTEGER(intk), INTENT(in) :: kk, jj, ii, q, l
    !     REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
    !     REAL(realk), INTENT(in) :: advr(kk, jj, ii)
    !     REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
    !     REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
    !     REAL(realk), INTENT(in) :: tol
    !     REAL(realk), INTENT(out) :: adve(kk, jj, ii)

    !     ! Loval variables
    !     INTEGER(intk) :: k, j, i
    !     INTEGER(intk) :: kl, jl, il
    !     REAL(realk) :: vel(kk,jj,ii)
    !     REAL(realk) :: signInd(2)
    !     REAL(realk) :: slopeLimLe, slopeLimRi

    !     CALL get_spatial_indices(kk, jj, ii, l, il, jl, kl)
    !     CALL get_condit_velocity(kk, jj, ii, q, u, v, w, vel)

    !     DO i = 2, ii-2
    !         DO j = 2, jj-2
    !             DO k = 2, kk-2
    !                 signInd(1) = MERGE(1.0_realk, 0.0_realk, advr(k,j,i) >= 0.0_realk)
    !                 signInd(2) = 1.0_realk - signInd(1)

    !                 slopeLimRi = 0.5_realk * ( SIGN(1.0_realk, ABS(vel(k+kl,j+jl,i+il))-ABS(vel(k,j,i)) * ( vel(k,j,i) - vel(k+kl,j+jl,i+il) ) + vel(k,j,i) + vel(k+kl,j+jl,i+il) ) )
    !                 slopeLimLe = 0.5_realk * ( SIGN(1.0_realk, ABS(vel(k,j,i))-ABS(vel(k-kl,j-jl,i-il)) * ( vel(k-kl,j-jl,i-il) - vel(k,j,i) ) + vel(k-kl,j-jl,i-il) + vel(k,j,i) ) )

    !                 adve(k,j,i) = signInd(1) * vel(k,j,i) + slopeLimRi * 0.5_realk + &
    !                               signInd(2) * vel(k+kl,j+jl,i+il) - slopeLimLe * 0.5_realk
    !             END DO
    !         END DO
    !     END DO

    ! END SUBROUTINE comp_adve_eno

    !================================================================

    SUBROUTINE comp_advr_linear_interpolation(kk, jj, ii, q, l, u, v, w, advr)
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

        ! Loval variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: kq, jq, iq
        REAL(realk) :: vel(kk,jj,ii)

        CALL get_spatial_indices(kk, jj, ii, q, iq, jq, kq)
        CALL get_condit_velocity(kk, jj, ii, l, u, v, w, vel)

        DO i = 2, ii-2
            DO j = 2, jj-2
                DO k = 2, kk-2
                    advr(k,j,i) = 0.5_realk * ( vel(k,j,i) + vel(k+kq,j+jq,i+iq) )
                END DO
            END DO
        END DO

    END SUBROUTINE comp_advr_linear_interpolation

END MODULE multiphase_vof_transport_mod
