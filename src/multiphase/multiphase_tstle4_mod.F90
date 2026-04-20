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

    IMPLICIT NONE
    PRIVATE

    ! PUBLIC :: multiphase_tstle4

CONTAINS
    
    ! SUBROUTINE multiphase_tstle4(uo_f, vo_f, wo_f, u_f, v_f, w_f, ut_f, vt_f, wt_f, &
    !             p_f, g_f, vff_f, itstep)
    ! !----------------------------------------------------------------
    ! !   What it does:
    ! !    
    ! !----------------------------------------------------------------

    !     ! Subroutine arguments
    !     TYPE(field_t), INTENT(inout) :: uo_f
    !     TYPE(field_t), INTENT(inout) :: vo_f
    !     TYPE(field_t), INTENT(inout) :: wo_f
    !     TYPE(field_t), INTENT(in) :: u_f
    !     TYPE(field_t), INTENT(in) :: v_f
    !     TYPE(field_t), INTENT(in) :: w_f
    !     TYPE(field_t), INTENT(in) :: ut_f
    !     TYPE(field_t), INTENT(in) :: vt_f
    !     TYPE(field_t), INTENT(in) :: wt_f
    !     TYPE(field_t), INTENT(in) :: p_f
    !     TYPE(field_t), INTENT(in) :: g_f
    !     TYPE(field_t), INTENT(in) :: vff_f
    !     INTEGER(intk), INTENT(in) :: itstep

    !     ! Local variables
    !     TYPE(field_t), POINTER :: dx_f, dy_f, dz_f, ddx_f, ddy_f, ddz_f
    !     TYPE(field_t), POINTER :: rdx_f, rdy_f, rdz_f, rddx_f, rddy_f, rddz_f
    !     REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: uo, vo, wo
    !     REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: vff, u, v, w
    !     REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: ut, vt, wt
    !     REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: p, g
    !     REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)
    !     REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
    !     REAL(realk), POINTER, CONTIGUOUS :: rdx(:), rdy(:), rdz(:)
    !     REAL(realk), POINTER, CONTIGUOUS :: rddx(:), rddy(:), rddz(:)
    !     INTEGER(intk) :: i, igrid, dim
    !     INTEGER(intk) :: kk, jj, ii
    !     INTEGER(intk) :: nfro, nbac, nrgt, nlft, nbot, ntop
    !     REAL(realk), ALLOCATABLE ::  densityFieldiStag(:,:,:), densityFieldjStag(:,:,:), densityFieldkStag(:,:,:)
    !     LOGICAL :: adv_x, adv_y, adv_z
    !     INTEGER(intk) :: permutationIndex

    !     ! Set all the output to zero everywhere before we start!
    !     uo_f = 0.0_realk
    !     vo_f = 0.0_realk
    !     wo_f = 0.0_realk

    !     CALL get_field(dx_f, "DX")
    !     CALL get_field(dy_f, "DY")
    !     CALL get_field(dz_f, "DZ")

    !     CALL get_field(ddx_f, "DDX")
    !     CALL get_field(ddy_f, "DDY")
    !     CALL get_field(ddz_f, "DDZ")

    !     CALL get_field(rdx_f, "RDX")
    !     CALL get_field(rdy_f, "RDY")
    !     CALL get_field(rdz_f, "RDZ")

    !     CALL get_field(rddx_f, "RDDX")
    !     CALL get_field(rddy_f, "RDDY")
    !     CALL get_field(rddz_f, "RDDZ")

    !     DO i = 1, nmygrids
    !         igrid = mygrids(i)

    !         CALL get_mgdims(kk, jj, ii, igrid)
    !         CALL get_mgbasb(nfro, nbac, nrgt, nlft, nbot, ntop, igrid)

    !         CALL uo_f%get_ptr(uo, igrid)
    !         CALL vo_f%get_ptr(vo, igrid)
    !         CALL wo_f%get_ptr(wo, igrid)

    !         CALL u_f%get_ptr(u, igrid)
    !         CALL v_f%get_ptr(v, igrid)
    !         CALL w_f%get_ptr(w, igrid)

    !         CALL ut_f%get_ptr(ut, igrid)
    !         CALL vt_f%get_ptr(vt, igrid)
    !         CALL wt_f%get_ptr(wt, igrid)

    !         CALL p_f%get_ptr(p, igrid)
    !         CALL g_f%get_ptr(g, igrid)
    !         CALL vff_f%get_ptr(vff, igrid)

    !         CALL dx_f%get_ptr(dx, igrid)
    !         CALL dy_f%get_ptr(dy, igrid)
    !         CALL dz_f%get_ptr(dz, igrid)

    !         CALL ddx_f%get_ptr(ddx, igrid)
    !         CALL ddy_f%get_ptr(ddy, igrid)
    !         CALL ddz_f%get_ptr(ddz, igrid)

    !         CALL rdx_f%get_ptr(rdx, igrid)
    !         CALL rdy_f%get_ptr(rdy, igrid)
    !         CALL rdz_f%get_ptr(rdz, igrid)

    !         CALL rddx_f%get_ptr(rddx, igrid)
    !         CALL rddy_f%get_ptr(rddy, igrid)
    !         CALL rddz_f%get_ptr(rddz, igrid)

    !         ! permutationIndex only changes in a new time-step
    !         permutationIndex = mod(itstep-1, 3)

    !         ! Select permutation of split advection
    !         SELECT CASE (permutationIndex)
    !             CASE (0)
    !                 adv_x = .TRUE.
    !                 adv_y = .FALSE.
    !                 adv_z = .FALSE.
    !             CASE (1)
    !                 adv_x = .FALSE.
    !                 adv_y = .TRUE.
    !                 adv_z = .FALSE.
    !             CASE (2)
    !                 adv_x = .FALSE.
    !                 adv_y = .FALSE.
    !                 adv_z = .TRUE.
    !         END SELECT

    !         IF ( .NOT. ALLOCATED(densityFieldiStag) .OR. .NOT. ALLOCATED(densityFieldjStag) .OR. .NOT. ALLOCATED(densityFieldkStag) ) THEN
    !             ALLOCATE(densityFieldiStag(kk, jj, ii), densityFieldjStag(kk, jj, ii), densityFieldkStag(kk, jj, ii))
    !         END IF

    !         CALL compute_normal_strain_rates(kk, jj, ii, strainRatex, strainRatey, strainRatez, vff, u, v, w, ddx, ddy, ddz)
    !         CALL compute_non_directional_compression_coeffiecient(kk, jj, ii, nonDirectionalCompressionCoefficient, vff)
    !         uCompressionTerm = nonDirectionalCompressionCoefficient * strainRatex
    !         vCompressionTerm = nonDirectionalCompressionCoefficient * strainRatey
    !         zCompressionTerm = nonDirectionalCompressionCoefficient * strainRatez

    !         xDensityCompressionTerm = ( nonDirectionalCompressionCoefficient * rho1 + ( 1 - nonDirectionalCompressionCoefficient ) * rho2 ) * strainRatex
    !         yDensityCompressionTerm = ( nonDirectionalCompressionCoefficient * rho1 + ( 1 - nonDirectionalCompressionCoefficient ) * rho2 ) * strainRatey
    !         zDensityCompressionTerm = ( nonDirectionalCompressionCoefficient * rho1 + ( 1 - nonDirectionalCompressionCoefficient ) * rho2 ) * strainRatez

    !         DO dim = 1, 3
    !             CALL interface_reconstruction_wrapper(kk, jj, ii, vff, ddx, ddy, ddz, tol, normx, normy, normz, alpha, isInterface)

    !             CALL compute_iStag_vff(kk, jj, ii, vffiStag, alpha, vff, isInterface, ddx, ddy, ddz, normx, normy, normz, tol)
    !             CALL compute_jStag_vff(kk, jj, ii, vffjStag, alpha, vff, isInterface, ddx, ddy, ddz, normx, normy, normz, tol)
    !             CALL compute_kStag_vff(kk, jj, ii, vffkStag, alpha, vff, isInterface, ddx, ddy, ddz, normx, normy, normz, tol)

    !             ! Get the old density fields
    !             CALL get_material_property_field(kk, jj, ii, densityFieldiStagOld, vffiStag, rho1, rho2)
    !             CALL get_material_property_field(kk, jj, ii, densityFieldjStagOld, vffjStag, rho1, rho2)
    !             CALL get_material_property_field(kk, jj, ii, densityFieldkStagOld, vffkStag, rho1, rho2)

    !             CALL interface_reconstruction_wrapper(kk, jj, ii, vffiStag, dx, ddy, ddz, tol, normxiStag, normyiStag, normziStag, alphaiStag, isInterfaceiStag, isNearInterfaceiStag)
    !             CALL interface_reconstruction_wrapper(kk, jj, ii, vffjStag, dx, ddy, ddz, tol, normxjStag, normyjStag, normzjStag, alphajStag, isInterfacejStag, isNearInterfaceiStag)
    !             CALL interface_reconstruction_wrapper(kk, jj, ii, vffkStag, dx, ddy, ddz, tol, normxkStag, normykStag, normzkStag, alphakStag, isInterfacekStag, isNearInterfaceiStag)

    !             IF ( adv_x ) THEN
    !                 ! Compute fluxes of volume fraction field
    !                 CALL compute_fluxx(vffFluxXiStag, kk, jj, ii, vffiStag, isInterfaceiStag, u, alphaiStag, dt, normxiStag, normyiStag, normziStag, dx, ddy, ddz, tol, & 
    !                                 nfro, nbac, nrgt, nlft, nbot, ntop)
    !                 CALL compute_fluxx(vffFluxXjStag, kk, jj, ii, vffjStag, isInterfacejStag, u, alphajStag, dt, normxjStag, normyjStag, normzjStag, dx, ddy, ddz, tol, & 
    !                                 nfro, nbac, nrgt, nlft, nbot, ntop)
    !                 CALL compute_fluxx(vffFluxXkStag, kk, jj, ii, vffkStag, isInterfacekStag, u, alphakStag, dt, normxkStag, normykStag, normzkStag, dx, ddy, ddz, tol, & 
    !                                 nfro, nbac, nrgt, nlft, nbot, ntop)
                    
    !                 ! Compute fluxes of complementary volume fraction field
    !                 vffCompiStag = 1 - vffiStag
    !                 vffCompjStag = 1 - vffjStag
    !                 vffCompkStag = 1 - vffkStag
    !                 CALL compute_fluxx(vffCompFluxXiStag, kk, jj, ii, vffCompiStag, isInterfaceiStag, u, alphaiStag, dt, normxiStag, normyiStag, normziStag, dx, ddy, ddz, tol, & 
    !                                 nfro, nbac, nrgt, nlft, nbot, ntop)
    !                 CALL compute_fluxx(vffCompFluxXjStag, kk, jj, ii, vffCompjStag, isInterfacejStag, u, alphajStag, dt, normxjStag, normyjStag, normzjStag, dx, ddy, ddz, tol, & 
    !                                 nfro, nbac, nrgt, nlft, nbot, ntop)
    !                 CALL compute_fluxx(vffCompFluxXkStag, kk, jj, ii, vffCompkStag, isInterfacekStag, u, alphakStag, dt, normxkStag, normykStag, normzkStag, dx, ddy, ddz, tol, & 
    !                                 nfro, nbac, nrgt, nlft, nbot, ntop)

    !                 ! Compute density fluxes from vff and vffComp fluxes
    !                 CALL get_density_flux(kk, jj, ii, vffFluxXiStag, vffCompFluxXiStag, rho1, rho2, densityFluxiStag)
    !                 CALL get_density_flux(kk, jj, ii, vffFluxXjStag, vffCompFluxXjStag, rho1, rho2, densityFluxjStag)
    !                 CALL get_density_flux(kk, jj, ii, vffFluxXkStag, vffCompFluxXkStag, rho1, rho2, densityFluxkStag)

    !                 ! Update vffs with fluxes
    !                 CALL update_field(kk, jj, ii, vffiStag, vffFluxXiStag, fluxy, fluxz, uCompressionTerm, vCompressionTerm, wCompressionTerm, & 
    !                                 adv_x, adv_y, adv_z, dt, nfro, nbac, nrgt, nlft, nbot, ntop)
    !                 CALL update_field(kk, jj, ii, vffjStag, vffFluxXjStag, fluxy, fluxz, uCompressionTerm, vCompressionTerm, wCompressionTerm, & 
    !                                 adv_x, adv_y, adv_z, dt, nfro, nbac, nrgt, nlft, nbot, ntop)
    !                 CALL update_field(kk, jj, ii, vffkStag, vffFluxXkStag, fluxy, fluxz, uCompressionTerm, vCompressionTerm, wCompressionTerm, & 
    !                                 adv_x, adv_y, adv_z, dt, nfro, nbac, nrgt, nlft, nbot, ntop)

    !                 ! Get the updated density fields
    !                 CALL get_material_property_field(kk, jj, ii, densityFieldiStag, vffiStag, rho1, rho2)
    !                 CALL get_material_property_field(kk, jj, ii, densityFieldjStag, vffjStag, rho1, rho2)
    !                 CALL get_material_property_field(kk, jj, ii, densityFieldkStag, vffkStag, rho1, rho2)
                    
    !             ELSE IF ( adv_y ) THEN
                    
    !             ELSE IF ( adv_z ) THEN
                    
    !             END IF

    !             ! Momentum advection
    !             CALL tstle4_kon(kk, jj, ii, uo, vo, wo, u, v, w, ut, vt, wt, &
    !                 dx, dy, dz, ddx, ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz, &
    !                 nfro, nbac, nrgt, nlft, nbot, ntop, itstep)
    !             CALL update_field(kk, jj, ii, field, fluxx, fluxy, fluxz, uDivergence, vDivergence, wDivergence, & 
    !                 adv_x, adv_y, adv_z, tol, ddx, ddy, ddz, dt, nfro, nbac, nrgt, nlft, nbot, ntop)
    !         END DO

    !         CALL tstle4_diff(kk, jj, ii, uo, vo, wo, u, v, w, g, &
    !             dx, dy, dz, ddx, ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz, &
    !             nfro, nbac, nrgt, nlft, nbot, ntop, densityFieldiStag, densityFieldjStag, densityFieldkStag)

    !         ! CALL tstle4_gradp(kk, jj, ii, uo, vo, wo, p, dx, dy, dz, &
    !         !     nfro, nbac, nrgt, nlft, nbot, ntop, igrid)

    !     END DO
        
    ! END SUBROUTINE multiphase_tstle4

    ! !================================================================

    ! SUBROUTINE multiphase_tstle4_kon()
    ! !----------------------------------------------------------------
    ! !   What it does:
    ! !    
    ! !----------------------------------------------------------------

    ! END SUBROUTINE multiphase_tstle4_kon

    ! !================================================================

    ! SUBROUTINE multiphase_tstle4_diff()
    ! !----------------------------------------------------------------
    ! !   What it does:
    ! !    
    ! !----------------------------------------------------------------

    ! END SUBROUTINE multiphase_tstle4_diff

    ! !================================================================

    ! SUBROUTINE multiphase_tstle4_gradp()
    ! !----------------------------------------------------------------
    ! !   What it does:
    ! !    
    ! !----------------------------------------------------------------

    ! END SUBROUTINE multiphase_tstle4_gradp

END MODULE multiphase_tstle4_mod