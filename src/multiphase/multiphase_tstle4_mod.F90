!====================================================================
!  Module: multiphase_tstle4_mod
!
!  Responsibilities:
!     - 
!
!  Author:      Valentin Ruhs
!  Created:     2026-02
!  Last update: 2026-02
!
!====================================================================

MODULE multiphase_tstle4_mod

    USE precision_mod, ONLY: intk, realk
    USE grids_mod, ONLY: nmygrids, mygrids
    USE field_mod, ONLY: field_t
    USE fields_mod, ONLY: get_field
    USE grids_mod, ONLY: get_mgdims, get_mgbasb
    USE err_mod, ONLY: errr
    USE lesmodel_mod, ONLY: ilesmodel
    USE multiphasecore_mod, ONLY: gmol1, gmol2, rho1, rho2
    USE multiphase_plic_mod, ONLY: interface_reconstruction_wrapper, staggered_fractions_wrapper
    USE multiphase_vof_transport_mod, ONLY: field_flux_wrapper, get_density_flux, get_advection_direction, compression_term_wrapper

    IMPLICIT NONE
    PRIVATE

    PUBLIC :: multiphase_tstle4

CONTAINS
    
    SUBROUTINE multiphase_tstle4(uo_f, vo_f, wo_f, u_f, v_f, w_f, ut_f, vt_f, wt_f, &
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
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: uo, vo, wo
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: u, v, w
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: ut, vt, wt
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: vff, p, g, d
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        REAL(realk), POINTER, CONTIGUOUS :: rdx(:), rdy(:), rdz(:)
        REAL(realk), POINTER, CONTIGUOUS :: rddx(:), rddy(:), rddz(:)
        INTEGER(intk) :: i, igrid, advectionDirection
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: nfro, nbac, nrgt, nlft, nbot, ntop
        REAL(realk), PARAMETER :: tol = 1.0E-15
        LOGICAL :: advX, advY, advZ
        REAL(realk), ALLOCATABLE :: normx(:,:,:), normy(:,:,:), normz(:,:,:)
        REAL(realk), ALLOCATABLE :: normxiStag(:,:,:), normyiStag(:,:,:), normziStag(:,:,:)
        REAL(realk), ALLOCATABLE :: normxjStag(:,:,:), normyjStag(:,:,:), normzjStag(:,:,:)
        REAL(realk), ALLOCATABLE :: normxkStag(:,:,:), normykStag(:,:,:), normzkStag(:,:,:)
        REAL(realk), ALLOCATABLE :: alpha(:,:,:), alphaiStag(:,:,:), alphajStag(:,:,:), alphakStag(:,:,:)
        LOGICAL, ALLOCATABLE :: isInterface(:,:,:), isInterfaceiStag(:,:,:), isInterfacejStag(:,:,:), isInterfacekStag(:,:,:)
        LOGICAL, ALLOCATABLE :: isNearInterface(:,:,:), isNearInterfaceiStag(:,:,:), isNearInterfacejStag(:,:,:), isNearInterfacekStag(:,:,:)
        REAL(realk), ALLOCATABLE :: vffiStag(:,:,:), vffjStag(:,:,:), vffkStag(:,:,:)
        REAL(realk), ALLOCATABLE :: densityFieldiStag(:,:,:), densityFieldjStag(:,:,:), densityFieldkStag(:,:,:)
        REAL(realk), ALLOCATABLE :: vffiStagFlux(:,:,:), vffjStagFlux(:,:,:), vffkStagFlux(:,:,:)
        REAL(realk), ALLOCATABLE :: complementvffiStagFlux(:,:,:), complementvffjStagFlux(:,:,:), complementvffkStagFlux(:,:,:)
        REAL(realk), ALLOCATABLE :: densityFluxiStag(:,:,:), densityFluxjStag(:,:,:), densityFluxkStag(:,:,:)
        REAL(realk), ALLOCATABLE :: densityCompressionTermXiStag(:,:,:), densityCompressionTermXjStag(:,:,:), densityCompressionTermXkStag(:,:,:)
        REAL(realk), ALLOCATABLE :: densityCompressionTermYiStag(:,:,:), densityCompressionTermYjStag(:,:,:), densityCompressionTermYkStag(:,:,:)
        REAL(realk), ALLOCATABLE :: densityCompressionTermZiStag(:,:,:), densityCompressionTermZjStag(:,:,:), densityCompressionTermZkStag(:,:,:)

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

            CALL rddx_f%get_ptr(rddx, igrid)
            CALL rddy_f%get_ptr(rddy, igrid)
            CALL rddz_f%get_ptr(rddz, igrid)

            IF (.NOT. ALLOCATED(normx)) ALLOCATE(normx(kk,jj,ii))
            IF (.NOT. ALLOCATED(normy)) ALLOCATE(normy(kk,jj,ii))
            IF (.NOT. ALLOCATED(normz)) ALLOCATE(normz(kk,jj,ii))

            IF (.NOT. ALLOCATED(normxiStag)) ALLOCATE(normxiStag(kk,jj,ii))
            IF (.NOT. ALLOCATED(normyiStag)) ALLOCATE(normyiStag(kk,jj,ii))
            IF (.NOT. ALLOCATED(normziStag)) ALLOCATE(normziStag(kk,jj,ii))

            IF (.NOT. ALLOCATED(normxjStag)) ALLOCATE(normxjStag(kk,jj,ii))
            IF (.NOT. ALLOCATED(normyjStag)) ALLOCATE(normyjStag(kk,jj,ii))
            IF (.NOT. ALLOCATED(normzjStag)) ALLOCATE(normzjStag(kk,jj,ii))

            IF (.NOT. ALLOCATED(normxkStag)) ALLOCATE(normxkStag(kk,jj,ii))
            IF (.NOT. ALLOCATED(normykStag)) ALLOCATE(normykStag(kk,jj,ii))
            IF (.NOT. ALLOCATED(normzkStag)) ALLOCATE(normzkStag(kk,jj,ii))

            IF (.NOT. ALLOCATED(alpha))      ALLOCATE(alpha(kk,jj,ii))
            IF (.NOT. ALLOCATED(alphaiStag)) ALLOCATE(alphaiStag(kk,jj,ii))
            IF (.NOT. ALLOCATED(alphajStag)) ALLOCATE(alphajStag(kk,jj,ii))
            IF (.NOT. ALLOCATED(alphakStag)) ALLOCATE(alphakStag(kk,jj,ii))

            IF (.NOT. ALLOCATED(isInterface))      ALLOCATE(isInterface(kk,jj,ii))
            IF (.NOT. ALLOCATED(isInterfaceiStag)) ALLOCATE(isInterfaceiStag(kk,jj,ii))
            IF (.NOT. ALLOCATED(isInterfacejStag)) ALLOCATE(isInterfacejStag(kk,jj,ii))
            IF (.NOT. ALLOCATED(isInterfacekStag)) ALLOCATE(isInterfacekStag(kk,jj,ii))

            IF (.NOT. ALLOCATED(isNearInterface))      ALLOCATE(isNearInterface(kk,jj,ii))
            IF (.NOT. ALLOCATED(isNearInterfaceiStag)) ALLOCATE(isNearInterfaceiStag(kk,jj,ii))
            IF (.NOT. ALLOCATED(isNearInterfacejStag)) ALLOCATE(isNearInterfacejStag(kk,jj,ii))
            IF (.NOT. ALLOCATED(isNearInterfacekStag)) ALLOCATE(isNearInterfacekStag(kk,jj,ii))

            IF (.NOT. ALLOCATED(vffiStag)) ALLOCATE(vffiStag(kk,jj,ii))
            IF (.NOT. ALLOCATED(vffjStag)) ALLOCATE(vffjStag(kk,jj,ii))
            IF (.NOT. ALLOCATED(vffkStag)) ALLOCATE(vffkStag(kk,jj,ii))

            IF (.NOT. ALLOCATED(densityFieldiStag)) ALLOCATE(densityFieldiStag(kk,jj,ii))
            IF (.NOT. ALLOCATED(densityFieldjStag)) ALLOCATE(densityFieldjStag(kk,jj,ii))
            IF (.NOT. ALLOCATED(densityFieldkStag)) ALLOCATE(densityFieldkStag(kk,jj,ii))

            IF (.NOT. ALLOCATED(vffiStagFlux)) ALLOCATE(vffiStagFlux(kk,jj,ii))
            IF (.NOT. ALLOCATED(vffjStagFlux)) ALLOCATE(vffjStagFlux(kk,jj,ii))
            IF (.NOT. ALLOCATED(vffkStagFlux)) ALLOCATE(vffkStagFlux(kk,jj,ii))

            IF (.NOT. ALLOCATED(complementvffiStagFlux)) ALLOCATE(complementvffiStagFlux(kk,jj,ii))
            IF (.NOT. ALLOCATED(complementvffjStagFlux)) ALLOCATE(complementvffjStagFlux(kk,jj,ii))
            IF (.NOT. ALLOCATED(complementvffkStagFlux)) ALLOCATE(complementvffkStagFlux(kk,jj,ii))

            IF (.NOT. ALLOCATED(densityFluxiStag)) ALLOCATE(densityFluxiStag(kk,jj,ii))
            IF (.NOT. ALLOCATED(densityFluxjStag)) ALLOCATE(densityFluxjStag(kk,jj,ii))
            IF (.NOT. ALLOCATED(densityFluxkStag)) ALLOCATE(densityFluxkStag(kk,jj,ii))

            IF (.NOT. ALLOCATED(densityCompressionTermXiStag)) ALLOCATE(densityCompressionTermXiStag(kk,jj,ii))
            IF (.NOT. ALLOCATED(densityCompressionTermXjStag)) ALLOCATE(densityCompressionTermXjStag(kk,jj,ii))
            IF (.NOT. ALLOCATED(densityCompressionTermXkStag)) ALLOCATE(densityCompressionTermXkStag(kk,jj,ii))

            IF (.NOT. ALLOCATED(densityCompressionTermYiStag)) ALLOCATE(densityCompressionTermYiStag(kk,jj,ii))
            IF (.NOT. ALLOCATED(densityCompressionTermYjStag)) ALLOCATE(densityCompressionTermYjStag(kk,jj,ii))
            IF (.NOT. ALLOCATED(densityCompressionTermYkStag)) ALLOCATE(densityCompressionTermYkStag(kk,jj,ii))

            IF (.NOT. ALLOCATED(densityCompressionTermZiStag)) ALLOCATE(densityCompressionTermZiStag(kk,jj,ii))
            IF (.NOT. ALLOCATED(densityCompressionTermZjStag)) ALLOCATE(densityCompressionTermZjStag(kk,jj,ii))
            IF (.NOT. ALLOCATED(densityCompressionTermZkStag)) ALLOCATE(densityCompressionTermZkStag(kk,jj,ii))

            CALL get_advection_direction(itstep, advX, advY, advZ)

            DO advectionDirection = 2, 4

                CALL interface_reconstruction_wrapper(kk, jj, ii, vff, ddx, ddy, ddz, tol, normx, normy, normz, alpha, isInterface, isNearInterface)

                CALL staggered_fractions_wrapper(kk, jj, ii, alpha, vff, isInterface, ddx, ddy, ddz, normx, normy, normz, tol, rho1, rho2, &
                    vffiStag, vffjStag, vffkStag, densityFieldiStag, densityFieldjStag, densityFieldkStag)

                CALL interface_reconstruction_wrapper(kk, jj, ii, vffiStag, dx, ddy, ddz, tol, normxiStag, normyiStag, normziStag, alphaiStag, isInterfaceiStag, isNearInterfaceiStag)
                CALL interface_reconstruction_wrapper(kk, jj, ii, vffjStag, ddx, dy, ddz, tol, normxjStag, normyjStag, normzjStag, alphajStag, isInterfacejStag, isNearInterfacejStag)
                CALL interface_reconstruction_wrapper(kk, jj, ii, vffkStag, ddx, ddy, dz, tol, normxkStag, normykStag, normzkStag, alphakStag, isInterfacekStag, isNearInterfacekStag)

                CALL field_flux_wrapper(kk, jj, ii, advX, advY, advZ, vffiStag, isInterfaceiStag, u, v, w, alphaiStag, dtrki, normxiStag, normyiStag, normziStag, dx, ddy, ddz, tol, nfro, nbac, nrgt, nlft, nbot, ntop, vffiStagFlux, complementvffiStagFlux)
                CALL field_flux_wrapper(kk, jj, ii, advX, advY, advZ, vffjStag, isInterfacejStag, u, v, w, alphajStag, dtrki, normxjStag, normyjStag, normzjStag, dx, ddy, ddz, tol, nfro, nbac, nrgt, nlft, nbot, ntop, vffjStagFlux, complementvffjStagFlux)
                CALL field_flux_wrapper(kk, jj, ii, advX, advY, advZ, vffkStag, isInterfacekStag, u, v, w, alphakStag, dtrki, normxkStag, normykStag, normzkStag, dx, ddy, ddz, tol, nfro, nbac, nrgt, nlft, nbot, ntop, vffkStagFlux, complementvffkStagFlux)

                CALL get_density_flux(kk, jj, ii, vffiStagFlux, complementvffiStagFlux, rho1, rho2, densityFluxiStag)
                CALL get_density_flux(kk, jj, ii, vffjStagFlux, complementvffjStagFlux, rho1, rho2, densityFluxjStag)
                CALL get_density_flux(kk, jj, ii, vffkStagFlux, complementvffkStagFlux, rho1, rho2, densityFluxkStag)

                CALL compression_term_wrapper(kk, jj, ii, u, v, w, vffiStag, dx, dy, dz, ddx, ddy, ddz, 1.0, 0.0, 0.0, rho1, rho2, densityCompressionTermXiStag, densityCompressionTermYiStag, densityCompressionTermZiStag)
                CALL compression_term_wrapper(kk, jj, ii, u, v, w, vffjStag, dx, dy, dz, ddx, ddy, ddz, 0.0, 1.0, 0.0, rho1, rho2, densityCompressionTermXjStag, densityCompressionTermYjStag, densityCompressionTermZjStag)
                CALL compression_term_wrapper(kk, jj, ii, u, v, w, vffjStag, dx, dy, dz, ddx, ddy, ddz, 0.0, 0.0, 1.0, rho1, rho2, densityCompressionTermXkStag, densityCompressionTermYkStag, densityCompressionTermZkStag)

                CALL multiphase_tstle4_kon()

                CALL get_advection_direction(advectionDirection, advX, advY, advZ)

            END DO

            CALL multiphase_tstle4_diff(kk, jj, ii, uo, vo, wo, u, v, w, g, &
                dx, dy, dz, ddx, ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz, &
                nfro, nbac, nrgt, nlft, nbot, ntop, densityFieldiStag, densityFieldjStag, densityFieldkStag)

            CALL multiphase_tstle4_gradp()

        END DO

    END SUBROUTINE multiphase_tstle4

    !================================================================

    SUBROUTINE multiphase_tstle4_kon()
    !----------------------------------------------------------------
    !   What it does:
    !    
    !----------------------------------------------------------------

        ! Subroutine arguments


        ! Local variables


        

    END SUBROUTINE multiphase_tstle4_kon

    !================================================================

    SUBROUTINE multiphase_tstle4_diff(kk, jj, ii, uo, vo, wo, u, v, w, g, &
            dx, dy, dz, ddx, ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz, &
            nfro, nbac, nrgt, nlft, nbot, ntop, densityFieldiStag, densityFieldjStag, densityFieldkStag)
    !----------------------------------------------------------------
    !   What it does:
    !    
    !----------------------------------------------------------------
    
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: uo(kk, jj, ii), vo(kk, jj, ii), &
            wo(kk, jj, ii)
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: g(kk, jj, ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: rdx(ii), rdy(jj), rdz(kk)
        REAL(realk), INTENT(in) :: rddx(ii), rddy(jj), rddz(kk)
        INTEGER, INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop
        REAL(realk), INTENT(in) :: densityFieldiStag(kk, jj, ii), densityFieldjStag(kk, jj, ii), densityFieldkStag(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: nbu, nfu, nrv, nbw, ntw, nlv
        INTEGER(intk) :: iles
        REAL(realk) :: ge, gw, gn, gs, gt, gb
        REAL(realk) :: tauxxe, tauxxw, tauyxn, tauyxs, tauzxt, tauzxb
        REAL(realk) :: tauxye, tauxyw, tauyyn, tauyys, tauzyt, tauzyb
        REAL(realk) :: tauxze, tauxzw, tauyzn, tauyzs, tauzzt, tauzzb
        REAL(realk) :: duo, dvo, dwo

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

        iles = 1
        IF (ilesmodel == 0) iles = 0

        ! CALL swcle3d(kk, jj, ii, uo, vo, wo, u, v, w, &
        !     ddx, ddy, ddz, nfro, nbac, nrgt, nlft, nbot, ntop)

        DO i = 3-nfu, ii-3+nbu
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
                    !             -----------inner derivatives-----------
                    tauxxe = ge * 2 * (u(k,j,i+1) - u(k,j,i)) * rddx(i+1)
                    tauxxw = gw * 2 * (u(k,j,i) - u(k,j,i-1)) * rddx(i)

                    ! Shear stresses
                    !             ------------------------------inner derivatives------------------------------
                    tauyxn = gn * ( (u(k,j+1,i) - u(k,j,i)) * rdy(j)   + (v(k,j,i+1) - v(k,j,i))     * rdx(i) )
                    tauyxs = gs * ( (u(k,j,i) - u(k,j-1,i)) * rdy(j-1) + (v(k,j-1,i+1) - v(k,j-1,i)) * rdx(i) )
                    tauzxt = gt * ( (u(k+1,j,i) - u(k,j,i)) * rdz(k)   + (w(k,j,i+1) - w(k,j,i))     * rdx(i) )
                    tauzxb = gb * ( (u(k,j,i) - u(k-1,j,i)) * rdz(k-1) + (w(k-1,j,i+1) - w(k-1,j,i)) * rdx(i) )

                    ! Change due to diffusion
                    !                                  ---------------------------------------outer derivatives----------------------------------------
                    duo = 1/densityFieldiStag(k,j,i) * ( ( tauxxe - tauxxw ) * rdx(i) + ( tauyxn - tauyxs ) * rddy(j) + ( tauzxt - tauzxb ) * rddz(k) )

                    IF ( duo > 0.0 ) THEN
                        WRITE(*,*) duo
                    END IF

                    ! Addition
                    uo(k, j, i) = uo(k, j, i) + duo
                END DO
            END DO
        END DO

        DO i = 3, ii-2
            DO j = 3-nrv, jj-3+nlv
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
                    !             ------------------------------inner derivatives------------------------------
                    tauxye = ge * ( (u(k,j+1,i) - u(k,j,i))     * rdy(j) + (v(k,j,i+1) - v(k,j,i)) * rdx(i)   )
                    tauxyw = gw * ( (u(k,j+1,i-1) - u(k,j,i-1)) * rdy(j) + (v(k,j,i) - v(k,j,i-1)) * rdx(i-1) )

                    ! Normal stresses
                    !             -----------inner derivatives-----------
                    tauyyn = gn * 2 * (v(k,j+1,i) - v(k,j,i)) * rddy(j+1)
                    tauyys = gs * 2 * (v(k,j,i) - v(k,j-1,i)) * rddy(j)
                    
                    ! Shear stresses
                    !             ------------------------------inner derivatives------------------------------
                    tauzyt = gt * ( (v(k+1,j,i) - v(k,j,i)) * rdz(k)   + (w(k,j+1,i) - w(k,j,i))     * rdy(j) )
                    tauzyb = gb * ( (v(k,j,i) - v(k-1,j,i)) * rdz(k-1) + (w(k-1,j+1,i) - w(k-1,j,i)) * rdy(j) )

                    ! Change due to diffusion
                    !                                  ---------------------------------------outer derivatives----------------------------------------
                    dvo = 1/densityFieldjStag(k,j,i) * ( ( tauxye - tauxyw ) * rddx(i) + ( tauyyn - tauyys ) * rdy(j) + ( tauzyt - tauzyb ) * rddz(k) )

                    ! Addition
                    vo(k, j, i) = vo(k, j, i) + dvo
                END DO
            END DO
        END DO

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3-nbw, kk-3+ntw
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
                    !             ------------------------------inner derivatives------------------------------
                    tauxze = ge * ( (u(k+1,j,i) - u(k,j,i))     * rdz(k) + (w(k,j,i+1) - w(k,j,i)) * rdx(i)   )
                    tauxzw = gw * ( (u(k+1,j,i-1) - u(k,j,i-1)) * rdz(k) + (w(k,j,i) - w(k,j,i-1)) * rdx(i-1) )
                    tauyzn = gn * ( (v(k+1,j,i) - v(k,j,i))     * rdz(k) + (w(k,j+1,i) - w(k,j,i)) * rdy(j)   )
                    tauyzs = gs * ( (v(k+1,j-1,i) - v(k,j-1,i)) * rdz(k) + (w(k,j,i) - w(k,j-1,i)) * rdy(j-1) )
                    
                    ! Normal stresses
                    !             -----------inner derivatives-----------
                    tauzzt = gt * 2 * (w(k+1,j,i) - w(k,j,i)) * rddz(k+1)
                    tauzzb = gb * 2 * (w(k,j,i) - w(k-1,j,i)) * rddz(k)

                    ! Change due to diffusion
                    !                                  ---------------------------------------outer derivatives----------------------------------------
                    dwo = 1/densityFieldkStag(k,j,i) * ( ( tauxze - tauxzw ) * rddx(i) + ( tauyzn - tauyzs ) * rddy(j) + ( tauzzt - tauzzb ) * rdz(k) )

                    ! Addition
                    wo(k, j, i) = wo(k, j, i) + dwo
                END DO
            END DO
        END DO 

    END SUBROUTINE multiphase_tstle4_diff

    !================================================================

    SUBROUTINE multiphase_tstle4_gradp()
    !----------------------------------------------------------------
    !   What it does:
    !    
    !----------------------------------------------------------------

        ! Subroutine arguments


        ! Local variables


        

    END SUBROUTINE multiphase_tstle4_gradp

    !================================================================

    ! PURE SUBROUTINE quick_interpolation_scheme(kk, jj, ii, k, j, i, adveField, &
    !     adveE, adveW, adveN, adveS, adveT, adveB, &
    !     advrE, advrW, advrN, advrS, advrT, advrB)
    ! !----------------------------------------------------------------
    ! !   What it does:
    ! !   The subroutine performes a QUICK interpolation for the 
    ! !   advected components of the momentum calculation.
    ! !   adve = advected component (advectee)
    ! !   advr = advecting component (advector)
    ! !----------------------------------------------------------------

    !     ! Subroutine arguments
    !     INTEGER(intk), INTENT(in) :: kk, jj, ii
    !     INTEGER(intk), INTENT(in) :: k, j, i
    !     REAL(realk), INTENT(in) :: adveField(kk, jj, ii)
    !     REAL(realk), INTENT(out) :: adveE, adveW, adveN, adveS, adveT, adveB
    !     REAL(realk), INTENT(in) :: advrE, advrW, advrN, advrS, advrT, advrB

    !     ! Loval variables
    !     ! None

    !     !       -----indicator-function----   --------------------------QUICK 3^rd order interpolation-------------------------
    !     adveE = 0.5 * ( 1.0 + SIGN(1.0,advrE) ) * 0.75 * adveField(k,j,i) + 0.375 * adveField(k,j,i+1) - 0.125 * adveField(k,j,i-1) + &
    !             0.5 * ( 1.0 - SIGN(1.0,advrE) ) * 0.75 * adveField(k,j,i+1) + 0.375 * adveField(k,j,i) - 0.125 * adveField(k,j,i+2)
    !     adveW = 0.5 * ( 1.0 + SIGN(1.0,advrW) ) * 0.75 * adveField(k,j,i-1) + 0.375 * adveField(k,j,i) - 0.125 * adveField(k,j,i-2) + &
    !             0.5 * ( 1.0 - SIGN(1.0,advrW) ) * 0.75 * adveField(k,j,i) + 0.375 * adveField(k,j,i-1) - 0.125 * adveField(k,j,i+1)
    !     adveN = 0.5 * ( 1.0 + SIGN(1.0,advrN) ) * 0.75 * adveField(k,j,i) + 0.375 * adveField(k,j+1,i) - 0.125 * adveField(k,j-1,i) + &
    !             0.5 * ( 1.0 - SIGN(1.0,advrN) ) * 0.75 * adveField(k,j+1,i) + 0.375 * adveField(k,j,i) - 0.125 * adveField(k,j+2,i)
    !     adveS = 0.5 * ( 1.0 + SIGN(1.0,advrS) ) * 0.75 * adveField(k,j-1,i) + 0.375 * adveField(k,j,i) - 0.125 * adveField(k,j-2,i) + &
    !             0.5 * ( 1.0 - SIGN(1.0,advrS) ) * 0.75 * adveField(k,j,i) + 0.375 * adveField(k,j-1,i) - 0.125 * adveField(k,j+1,i)
    !     adveT = 0.5 * ( 1.0 + SIGN(1.0,advrT) ) * 0.75 * adveField(k,j,i) + 0.375 * adveField(k+1,j,i) - 0.125 * adveField(k-1,j,i) + &
    !             0.5 * ( 1.0 - SIGN(1.0,advrT) ) * 0.75 * adveField(k+1,j,i) + 0.375 * adveField(k,j,i) - 0.125 * adveField(k+2,j,i)
    !     adveB = 0.5 * ( 1.0 + SIGN(1.0,advrB) ) * 0.75 * adveField(k-1,j,i) + 0.375 * adveField(k,j,i) - 0.125 * adveField(k-2,j,i) + &
    !             0.5 * ( 1.0 - SIGN(1.0,advrB) ) * 0.75 * adveField(k,j,i) + 0.375 * adveField(k-1,j,i) - 0.125 * adveField(k+1,j,i)

    ! END SUBROUTINE quick_interpolation_scheme

    !================================================================

    ! PURE SUBROUTINE advecting_interpolation_scheme(kk, jj, ii, k, j, i, adveField, &
    !     advrE, advrW, advrN, advrS, advrT, advrB)
    ! !----------------------------------------------------------------
    ! !   What it does:
    ! !   
    ! !----------------------------------------------------------------

    !     ! Subroutine arguments
    !     INTEGER(intk), INTENT(in) :: kk, jj, ii
    !     INTEGER(intk), INTENT(in) :: k, j, i
    !     REAL(realk), INTENT(in) :: adveField(kk, jj, ii)
    !     REAL(realk), INTENT(out) :: advrE, advrW, advrN, advrS, advrT, advrB

    !     ! Loval variables
    !     ! None

    !     advrE = 0.5 * ( u(k,j,i) + u(k,j,i+1) )
    !     advrW = 0.5 * ( u(k,j,i-1) + u(k,j,i) )
    !     advrN = 0.5 * ( v(k,j,i) + v(k,j,i+1) )
    !     advrS = 0.5 * ( v(k,j-1,i) + v(k,j-1,i+1) )
    !     advrT = 0.5 * ( w(k,j,i) + w(k,j,i+1) )
    !     advrB = 0.5 * ( w(k-1,j,i) + w(k-1,j,i+1) )

    ! END SUBROUTINE advecting_interpolation_scheme


END MODULE multiphase_tstle4_mod