!====================================================================
!  Module: multiphase_advection_mod
!
!  Responsibilities:
!     - 
!
!  Author:      Valentin Ruhs
!  Created:     2026-02
!  Last update: 2026-02
!
!====================================================================

MODULE multiphase_advection_mod

    USE precision_mod, ONLY: intk, realk
    USE grids_mod, ONLY: nmygrids, mygrids
    USE field_mod, ONLY: field_t
    USE fields_mod, ONLY: get_field
    USE grids_mod, ONLY: get_mgdims, get_mgbasb
    USE err_mod, ONLY: errr
    USE multiphasecore_mod, ONLY: gmol1, gmol2, rho1, rho2
    USE multiphase_material_mod, ONLY: compute_material_property_field
    USE multiphase_plic_mod, ONLY: interface_reconstruction_wrapper, staggered_fractions_wrapper
    USE multiphase_vof_transport_mod, ONLY: compute_flux, compute_density_flux, get_advection_sequence, compression_term_wrapper, update_field, clip_volume_fraction_field

    IMPLICIT NONE
    PRIVATE

    PUBLIC :: multiphase_split_advection

CONTAINS
    
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
        REAL(realk) :: advrE, advrW, advrN, advrS, advrT, advrB
        REAL(realk) :: adveE, adveW, adveN, adveS, adveT, adveB
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

        IF ( splitDir == 1 ) THEN
            DO i = 3-nfu, ii-3+nbu
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        CALL advecting_interpolation_scheme(kk, jj, ii, k, j, i, u, v, w, &
                            advrE, advrW, advrN, advrS, advrT, advrB, iStag, jStag, kStag)

                        CALL quick_advected_interpolation_scheme(kk, jj, ii, k, j, i, velocity, &
                            adveE, adveW, adveN, adveS, adveT, adveB, &
                            advrE, advrW, advrN, advrS, advrT, advrB)

                        IF ( isNearInterfaceStag(k,j,i,component) ) THEN
                            dMomentum = - ( adveE * densityFieldFluxStag(k,j,i,component) - adveW * densityFieldFluxStag(k,j,i-1,component) ) + velocity(k,j,i) * densityCompressionTermStag(k,j,i,component)
                            velocity(k,j,i) = 1.0_realk / densityFieldStag(k,j,i,component) * ( densityFieldStagOld(k,j,i,component) * velocity(k,j,i) + dMomentum )
                        ELSE
                            dVelocity = - ( ( adveE * advrE - adveW * advrW ) * rdx(i) )
                            velocity(k,j,i) = velocity(k,j,i) + dVelocity
                        END IF
                    END DO
                END DO
            END DO
        ELSE IF ( splitDir == 2 ) THEN
            DO i = 3, ii-2
                DO j = 3-nrv, jj-3+nlv
                    DO k = 3, kk-2
                        CALL advecting_interpolation_scheme(kk, jj, ii, k, j, i, u, v, w, &
                            advrE, advrW, advrN, advrS, advrT, advrB, iStag, jStag, kStag)

                        CALL quick_advected_interpolation_scheme(kk, jj, ii, k, j, i, velocity, &
                            adveE, adveW, adveN, adveS, adveT, adveB, &
                            advrE, advrW, advrN, advrS, advrT, advrB)

                        IF ( isNearInterfaceStag(k,j,i,component) ) THEN
                            dMomentum = - ( adveN * densityFieldFluxStag(k,j,i,component) - adveS * densityFieldFluxStag(k,j-1,i,component) ) + velocity(k,j,i) * densityCompressionTermStag(k,j,i,component)
                            velocity(k,j,i) = 1.0_realk / densityFieldStag(k,j,i,component) * ( densityFieldStagOld(k,j,i,component) * velocity(k,j,i) + dMomentum )
                        ELSE
                            dVelocity = - ( ( adveN * advrN - adveS * advrS ) * rdy(j) )
                            velocity(k,j,i) = velocity(k,j,i) + dVelocity
                        END IF
                    END DO
                END DO
            END DO
        ELSE IF ( splitDir == 3 ) THEN
            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3-nbw, kk-3+ntw
                        CALL advecting_interpolation_scheme(kk, jj, ii, k, j, i, u, v, w, &
                            advrE, advrW, advrN, advrS, advrT, advrB, iStag, jStag, kStag)

                        CALL quick_advected_interpolation_scheme(kk, jj, ii, k, j, i, velocity, &
                            adveE, adveW, adveN, adveS, adveT, adveB, &
                            advrE, advrW, advrN, advrS, advrT, advrB)

                        IF ( isNearInterfaceStag(k,j,i,component) ) THEN
                            dMomentum = - ( adveT * densityFieldFluxStag(k,j,i,component) - adveB * densityFieldFluxStag(k-1,j,i,component) ) + velocity(k,j,i) * densityCompressionTermStag(k,j,i,component)
                            velocity(k,j,i) = 1.0_realk / densityFieldStag(k,j,i,component) * ( densityFieldStagOld(k,j,i,component) * velocity(k,j,i) + dMomentum )
                        ELSE
                            dVelocity = - ( ( adveT * advrT - adveB * advrB ) * rdz(k) )
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

    PURE SUBROUTINE quick_advected_interpolation_scheme(kk, jj, ii, k, j, i, adveVelocity, &
        adveE, adveW, adveN, adveS, adveT, adveB, &
        advrE, advrW, advrN, advrS, advrT, advrB)
    !----------------------------------------------------------------
    !   What it does:
    !   The subroutine performes a QUICK interpolation for the 
    !   advected components of the momentum calculation.
    !   adve = advected component (advectee)
    !   advr = advecting component (advector)
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: k, j, i
        REAL(realk), INTENT(in) :: adveVelocity(kk, jj, ii)
        REAL(realk), INTENT(out) :: adveE, adveW, adveN, adveS, adveT, adveB
        REAL(realk), INTENT(in) :: advrE, advrW, advrN, advrS, advrT, advrB

        ! Loval variables
        ! None

        !                   ----------indicator-function---------   -----------------------------------------QUICK 3^rd order interpolation-----------------------------------------
        adveE = 0.5_realk * ( 1.0_realk + SIGN(1.0_realk,advrE) ) * ( 0.75_realk * adveVelocity(k,j,i) + 0.375_realk * adveVelocity(k,j,i+1) - 0.125_realk * adveVelocity(k,j,i-1) ) + &
                0.5_realk * ( 1.0_realk - SIGN(1.0_realk,advrE) ) * ( 0.75_realk * adveVelocity(k,j,i+1) + 0.375_realk * adveVelocity(k,j,i) - 0.125_realk * adveVelocity(k,j,i+2) )
        adveW = 0.5_realk * ( 1.0_realk + SIGN(1.0_realk,advrW) ) * ( 0.75_realk * adveVelocity(k,j,i-1) + 0.375_realk * adveVelocity(k,j,i) - 0.125_realk * adveVelocity(k,j,i-2) ) + &
                0.5_realk * ( 1.0_realk - SIGN(1.0_realk,advrW) ) * ( 0.75_realk * adveVelocity(k,j,i) + 0.375_realk * adveVelocity(k,j,i-1) - 0.125_realk * adveVelocity(k,j,i+1) )
        adveN = 0.5_realk * ( 1.0_realk + SIGN(1.0_realk,advrN) ) * ( 0.75_realk * adveVelocity(k,j,i) + 0.375_realk * adveVelocity(k,j+1,i) - 0.125_realk * adveVelocity(k,j-1,i) ) + &
                0.5_realk * ( 1.0_realk - SIGN(1.0_realk,advrN) ) * ( 0.75_realk * adveVelocity(k,j+1,i) + 0.375_realk * adveVelocity(k,j,i) - 0.125_realk * adveVelocity(k,j+2,i) )
        adveS = 0.5_realk * ( 1.0_realk + SIGN(1.0_realk,advrS) ) * ( 0.75_realk * adveVelocity(k,j-1,i) + 0.375_realk * adveVelocity(k,j,i) - 0.125_realk * adveVelocity(k,j-2,i) ) + &
                0.5_realk * ( 1.0_realk - SIGN(1.0_realk,advrS) ) * ( 0.75_realk * adveVelocity(k,j,i) + 0.375_realk * adveVelocity(k,j-1,i) - 0.125_realk * adveVelocity(k,j+1,i) )
        adveT = 0.5_realk * ( 1.0_realk + SIGN(1.0_realk,advrT) ) * ( 0.75_realk * adveVelocity(k,j,i) + 0.375_realk * adveVelocity(k+1,j,i) - 0.125_realk * adveVelocity(k-1,j,i) ) + &
                0.5_realk * ( 1.0_realk - SIGN(1.0_realk,advrT) ) * ( 0.75_realk * adveVelocity(k+1,j,i) + 0.375_realk * adveVelocity(k,j,i) - 0.125_realk * adveVelocity(k+2,j,i) )
        adveB = 0.5_realk * ( 1.0_realk + SIGN(1.0_realk,advrB) ) * ( 0.75_realk * adveVelocity(k-1,j,i) + 0.375_realk * adveVelocity(k,j,i) - 0.125_realk * adveVelocity(k-2,j,i) ) + &
                0.5_realk * ( 1.0_realk - SIGN(1.0_realk,advrB) ) * ( 0.75_realk * adveVelocity(k,j,i) + 0.375_realk * adveVelocity(k-1,j,i) - 0.125_realk * adveVelocity(k+1,j,i) )

    END SUBROUTINE quick_advected_interpolation_scheme

    !================================================================

    ! PURE SUBROUTINE quick_advected_interpolation_scheme_test(kk, jj, ii, k, j, i, adveVelocity, &
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
    !     REAL(realk), INTENT(in) :: adveVelocity(kk, jj, ii)
    !     REAL(realk), INTENT(out) :: adveE, adveW, adveN, adveS, adveT, adveB
    !     REAL(realk), INTENT(in) :: advrE, advrW, advrN, advrS, advrT, advrB

    !     ! Loval variables
    !     ! None

    !     IF ( advrE >= 0.0_realk ) THEN
    !         adveE = ( 0.75_realk * adveVelocity(k,j,i) + 0.375_realk * adveVelocity(k,j,i+1) - 0.125_realk * adveVelocity(k,j,i-1) )
    !     ELSE IF ( advrE > 0.0_realk ) THEN
    !         adveE = ( 0.75_realk * adveVelocity(k,j,i+1) + 0.375_realk * adveVelocity(k,j,i) - 0.125_realk * adveVelocity(k,j,i+2) )
    !     END IF

    !     IF ( advrW >= 0.0_realk ) THEN
    !         adveW = ( 0.75_realk * adveVelocity(k,j,i-1) + 0.375_realk * adveVelocity(k,j,i) - 0.125_realk * adveVelocity(k,j,i-2) )
    !     ELSE IF ( advrW < 0.0_realk ) THEN
    !         adveW = ( 0.75_realk * adveVelocity(k,j,i) + 0.375_realk * adveVelocity(k,j,i-1) - 0.125_realk * adveVelocity(k,j,i+1) )
    !     END IF

    !     IF ( advrN >= 0.0_realk ) THEN
    !         adveN = ( 0.75_realk * adveVelocity(k,j,i) + 0.375_realk * adveVelocity(k,j+1,i) - 0.125_realk * adveVelocity(k,j-1,i) )
    !     ELSE IF ( advrN < 0.0_realk ) THEN
    !         adveN = ( 0.75_realk * adveVelocity(k,j+1,i) + 0.375_realk * adveVelocity(k,j,i) - 0.125_realk * adveVelocity(k,j+2,i) )
    !     END IF

    !     IF ( advrS >= 0.0_realk ) THEN
    !         adveS = ( 0.75_realk * adveVelocity(k,j-1,i) + 0.375_realk * adveVelocity(k,j,i) - 0.125_realk * adveVelocity(k,j-2,i) )
    !     ELSE IF ( advrS < 0.0_realk ) THEN
    !         adveS = ( 0.75_realk * adveVelocity(k,j,i) + 0.375_realk * adveVelocity(k,j-1,i) - 0.125_realk * adveVelocity(k,j+1,i) )
    !     END IF

    !     IF ( advrT >= 0.0_realk ) THEN
    !         adveT = ( 0.75_realk * adveVelocity(k,j,i) + 0.375_realk * adveVelocity(k+1,j,i) - 0.125_realk * adveVelocity(k-1,j,i) )
    !     ELSE IF ( advrT < 0.0_realk ) THEN
    !         adveT = ( 0.75_realk * adveVelocity(k+1,j,i) + 0.375_realk * adveVelocity(k,j,i) - 0.125_realk * adveVelocity(k+2,j,i) )
    !     END IF

    !     IF ( advrB >= 0.0_realk ) THEN
    !         adveB = ( 0.75_realk * adveVelocity(k-1,j,i) + 0.375_realk * adveVelocity(k,j,i) - 0.125_realk * adveVelocity(k-2,j,i) )
    !     ELSE IF ( advrB < 0.0_realk ) THEN
    !         adveB = ( 0.75_realk * adveVelocity(k,j,i) + 0.375_realk * adveVelocity(k-1,j,i) - 0.125_realk * adveVelocity(k+1,j,i) )
    !     END IF

    ! END SUBROUTINE quick_advected_interpolation_scheme_test

    !================================================================

    PURE SUBROUTINE advecting_interpolation_scheme(kk, jj, ii, k, j, i, u, v, w, &
        advrE, advrW, advrN, advrS, advrT, advrB, iStag, jStag, kStag)
    !----------------------------------------------------------------
    !   What it does:
    !   The subroutine performes the average interpolation for the 
    !   advecting components of the momentum calculation.
    !   advr = advecting component (advector)
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: k, j, i
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(out) :: advrE, advrW, advrN, advrS, advrT, advrB
        REAL(realk), INTENT(in) :: iStag, jStag, kStag

        ! Loval variables
        ! None

        advrE = 0.5_realk * ( iStag * ( u(k,j,i) + u(k,j,i+1) ) + jStag * ( u(k,j,i) + u(k,j+1,i) ) + kStag * ( u(k,j,i) + u(k+1,j,i) ) )
        advrW = 0.5_realk * ( iStag * ( u(k,j,i-1) + u(k,j,i) ) + jStag * ( u(k,j,i-1) + u(k,j+1,i-1) ) + kStag * ( u(k,j,i-1) + u(k+1,j,i-1) ) )
        advrN = 0.5_realk * ( iStag * ( v(k,j,i) + v(k,j,i+1) ) + jStag * ( v(k,j,i) + v(k,j+1,i) ) + kStag * ( v(k,j,i) + v(k+1,j,i) ) )
        advrS = 0.5_realk * ( iStag * ( v(k,j-1,i) + v(k,j-1,i+1) ) + jStag * ( v(k,j-1,i) + v(k,j,i) ) + kStag * ( v(k,j-1,i) + v(k+1,j-1,i) ) )
        advrT = 0.5_realk * ( iStag * ( w(k,j,i) + w(k,j,i+1) ) + jStag * ( w(k,j,i) + w(k,j+1,i) ) + kStag * ( w(k,j,i) + w(k+1,j,i) ) )
        advrB = 0.5_realk * ( iStag * ( w(k-1,j,i) + w(k-1,j,i+1) ) + jStag * ( w(k-1,j,i) + w(k-1,j+1,i) ) + kStag * ( w(k-1,j,i) + w(k,j,i) ) )

    END SUBROUTINE advecting_interpolation_scheme

END MODULE multiphase_advection_mod