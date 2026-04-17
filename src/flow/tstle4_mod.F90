MODULE tstle4_mod
    USE core_mod
    USE flowcore_mod
    USE lesmodel_mod, ONLY: ilesmodel
    USE wernerwengle_mod, ONLY: tauwin
    USE multiphasecore_mod, ONLY: gmol1, gmol2, rho1, rho2
    USE multiphase_material_mod, ONLY: get_material_property_field

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: tstle4

CONTAINS
    SUBROUTINE tstle4(uo_f, vo_f, wo_f, u_f, v_f, w_f, ut_f, vt_f, wt_f, &
            p_f, g_f, vff_f, itstep)
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
        TYPE(field_t), INTENT(in) :: p_f
        TYPE(field_t), INTENT(in) :: g_f
        TYPE(field_t), INTENT(in) :: vff_f
        INTEGER(intk), INTENT(in) :: itstep

        ! Local variables
        TYPE(field_t), POINTER :: dx_f, dy_f, dz_f, ddx_f, ddy_f, ddz_f
        TYPE(field_t), POINTER :: rdx_f, rdy_f, rdz_f, rddx_f, rddy_f, rddz_f
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: uo, vo, wo
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: vff, u, v, w
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: ut, vt, wt
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: p, g
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        REAL(realk), POINTER, CONTIGUOUS :: rdx(:), rdy(:), rdz(:)
        REAL(realk), POINTER, CONTIGUOUS :: rddx(:), rddy(:), rddz(:)
        INTEGER(intk) :: i, igrid, dim
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: nfro, nbac, nrgt, nlft, nbot, ntop
        REAL(realk), ALLOCATABLE ::  densityFieldiStag(:,:,:), densityFieldjStag(:,:,:), densityFieldkStag(:,:,:)
        LOGICAL :: adv_x, adv_y, adv_z
        INTEGER(intk) :: permutationIndex

        CALL start_timer(310)

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

            CALL p_f%get_ptr(p, igrid)
            CALL g_f%get_ptr(g, igrid)
            CALL vff_f%get_ptr(vff, igrid)

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

            ! permutationIndex only changes in a new time-step
            permutationIndex = mod(itstep-1, 3)

            ! Select permutation of split advection
            SELECT CASE (permutationIndex)
                CASE (0)
                    adv_x = .TRUE.
                    adv_y = .FALSE.
                    adv_z = .FALSE.
                CASE (1)
                    adv_x = .FALSE.
                    adv_y = .TRUE.
                    adv_z = .FALSE.
                CASE (2)
                    adv_x = .FALSE.
                    adv_y = .FALSE.
                    adv_z = .TRUE.
            END SELECT

            CALL get_material_property_field(kk, jj, ii, g, vff, gmol1, gmol2)

            IF ( .NOT. ALLOCATED(densityFieldiStag) .OR. .NOT. ALLOCATED(densityFieldjStag) .OR. .NOT. ALLOCATED(densityFieldkStag) ) THEN
                ALLOCATE(densityFieldiStag(kk, jj, ii), densityFieldjStag(kk, jj, ii), densityFieldkStag(kk, jj, ii))
            END IF

            CALL compute_normal_strain_rates(kk, jj, ii, strainRatex, strainRatey, strainRatez, vff, u, v, w, ddx, ddy, ddz)
            CALL compute_non_directional_compression_coeffiecient(kk, jj, ii, nonDirectionalCompressionCoefficient, vff)
            uCompressionTerm = nonDirectionalCompressionCoefficient * strainRatex
            vCompressionTerm = nonDirectionalCompressionCoefficient * strainRatey
            zCompressionTerm = nonDirectionalCompressionCoefficient * strainRatez

            xDensityCompressionTerm = ( nonDirectionalCompressionCoefficient * rho1 + ( 1 - nonDirectionalCompressionCoefficient ) * rho2 ) * strainRatex
            yDensityCompressionTerm = ( nonDirectionalCompressionCoefficient * rho1 + ( 1 - nonDirectionalCompressionCoefficient ) * rho2 ) * strainRatey
            zDensityCompressionTerm = ( nonDirectionalCompressionCoefficient * rho1 + ( 1 - nonDirectionalCompressionCoefficient ) * rho2 ) * strainRatez

            DO dim = 1, 3
                CALL interface_reconstruction_wrapper(kk, jj, ii, vff, ddx, ddy, ddz, tol, normx, normy, normz, alpha, isInterface)

                CALL compute_iStag_vff(kk, jj, ii, vffiStag, alpha, vff, isInterface, ddx, ddy, ddz, normx, normy, normz, tol)
                CALL compute_jStag_vff(kk, jj, ii, vffjStag, alpha, vff, isInterface, ddx, ddy, ddz, normx, normy, normz, tol)
                CALL compute_kStag_vff(kk, jj, ii, vffkStag, alpha, vff, isInterface, ddx, ddy, ddz, normx, normy, normz, tol)

                ! Get the old density fields
                CALL get_material_property_field(kk, jj, ii, densityFieldiStagOld, vffiStag, rho1, rho2)
                CALL get_material_property_field(kk, jj, ii, densityFieldjStagOld, vffjStag, rho1, rho2)
                CALL get_material_property_field(kk, jj, ii, densityFieldkStagOld, vffkStag, rho1, rho2)

                CALL interface_reconstruction_wrapper(kk, jj, ii, vffiStag, dx, ddy, ddz, tol, normxiStag, normyiStag, normziStag, alphaiStag, isInterfaceiStag, isNearInterfaceiStag)
                CALL interface_reconstruction_wrapper(kk, jj, ii, vffjStag, dx, ddy, ddz, tol, normxjStag, normyjStag, normzjStag, alphajStag, isInterfacejStag, isNearInterfaceiStag)
                CALL interface_reconstruction_wrapper(kk, jj, ii, vffkStag, dx, ddy, ddz, tol, normxkStag, normykStag, normzkStag, alphakStag, isInterfacekStag, isNearInterfaceiStag)

                IF ( adv_x ) THEN
                    ! Compute fluxes of volume fraction field
                    CALL compute_fluxx(vffFluxXiStag, kk, jj, ii, vffiStag, isInterfaceiStag, u, alphaiStag, dt, normxiStag, normyiStag, normziStag, dx, ddy, ddz, tol, & 
                                       nfro, nbac, nrgt, nlft, nbot, ntop)
                    CALL compute_fluxx(vffFluxXjStag, kk, jj, ii, vffjStag, isInterfacejStag, u, alphajStag, dt, normxjStag, normyjStag, normzjStag, dx, ddy, ddz, tol, & 
                                       nfro, nbac, nrgt, nlft, nbot, ntop)
                    CALL compute_fluxx(vffFluxXkStag, kk, jj, ii, vffkStag, isInterfacekStag, u, alphakStag, dt, normxkStag, normykStag, normzkStag, dx, ddy, ddz, tol, & 
                                       nfro, nbac, nrgt, nlft, nbot, ntop)
                    
                    ! Compute fluxes of complementary volume fraction field
                    vffCompiStag = 1 - vffiStag
                    vffCompjStag = 1 - vffjStag
                    vffCompkStag = 1 - vffkStag
                    CALL compute_fluxx(vffCompFluxXiStag, kk, jj, ii, vffCompiStag, isInterfaceiStag, u, alphaiStag, dt, normxiStag, normyiStag, normziStag, dx, ddy, ddz, tol, & 
                                       nfro, nbac, nrgt, nlft, nbot, ntop)
                    CALL compute_fluxx(vffCompFluxXjStag, kk, jj, ii, vffCompjStag, isInterfacejStag, u, alphajStag, dt, normxjStag, normyjStag, normzjStag, dx, ddy, ddz, tol, & 
                                       nfro, nbac, nrgt, nlft, nbot, ntop)
                    CALL compute_fluxx(vffCompFluxXkStag, kk, jj, ii, vffCompkStag, isInterfacekStag, u, alphakStag, dt, normxkStag, normykStag, normzkStag, dx, ddy, ddz, tol, & 
                                       nfro, nbac, nrgt, nlft, nbot, ntop)

                    ! Compute density fluxes from vff and vffComp fluxes
                    CALL get_density_flux(kk, jj, ii, vffFluxXiStag, vffCompFluxXiStag, rho1, rho2, densityFluxiStag)
                    CALL get_density_flux(kk, jj, ii, vffFluxXjStag, vffCompFluxXjStag, rho1, rho2, densityFluxjStag)
                    CALL get_density_flux(kk, jj, ii, vffFluxXkStag, vffCompFluxXkStag, rho1, rho2, densityFluxkStag)

                    ! Update vffs with fluxes
                    CALL update_field(kk, jj, ii, vffiStag, vffFluxXiStag, fluxy, fluxz, uCompressionTerm, vCompressionTerm, wCompressionTerm, & 
                                      adv_x, adv_y, adv_z, dt, nfro, nbac, nrgt, nlft, nbot, ntop)
                    CALL update_field(kk, jj, ii, vffjStag, vffFluxXjStag, fluxy, fluxz, uCompressionTerm, vCompressionTerm, wCompressionTerm, & 
                                      adv_x, adv_y, adv_z, dt, nfro, nbac, nrgt, nlft, nbot, ntop)
                    CALL update_field(kk, jj, ii, vffkStag, vffFluxXkStag, fluxy, fluxz, uCompressionTerm, vCompressionTerm, wCompressionTerm, & 
                                      adv_x, adv_y, adv_z, dt, nfro, nbac, nrgt, nlft, nbot, ntop)

                    ! Get the updated density fields
                    CALL get_material_property_field(kk, jj, ii, densityFieldiStag, vffiStag, rho1, rho2)
                    CALL get_material_property_field(kk, jj, ii, densityFieldjStag, vffjStag, rho1, rho2)
                    CALL get_material_property_field(kk, jj, ii, densityFieldkStag, vffkStag, rho1, rho2)
                    
                ELSE IF ( adv_y ) THEN
                    
                ELSE IF ( adv_z ) THEN
                    
                END IF

                ! Momentum advection
                CALL tstle4_kon(kk, jj, ii, uo, vo, wo, u, v, w, ut, vt, wt, &
                    dx, dy, dz, ddx, ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz, &
                    nfro, nbac, nrgt, nlft, nbot, ntop, itstep)
                CALL update_field(kk, jj, ii, field, fluxx, fluxy, fluxz, uDivergence, vDivergence, wDivergence, & 
                    adv_x, adv_y, adv_z, tol, ddx, ddy, ddz, dt, nfro, nbac, nrgt, nlft, nbot, ntop)
            END DO

            CALL tstle4_diff(kk, jj, ii, uo, vo, wo, u, v, w, g, &
                dx, dy, dz, ddx, ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz, &
                nfro, nbac, nrgt, nlft, nbot, ntop, densityFieldiStag, densityFieldjStag, densityFieldkStag)

            ! CALL tstle4_gradp(kk, jj, ii, uo, vo, wo, p, dx, dy, dz, &
            !     nfro, nbac, nrgt, nlft, nbot, ntop, igrid)

            ! CALL tstle4_par(kk, jj, ii, uo, vo, wo, u, v, w, ut, vt, wt, &
            !     dx, dy, dz, ddx, ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz, &
            !     nfro, nbac, nrgt, nlft, nbot, ntop)
        END DO

        CALL stop_timer(310)
    END SUBROUTINE tstle4


    ! The convective terms are computed in two steps:
    ! First, the mass fluxes (transporting velocities) are interpolated to
    ! the faces of the momentum cell. This interpolation is performed in a
    ! way which ensures mass conservation at the momentum cell if the
    ! adveField field is divergence-free on the adjacent pressure cells.
    ! Second, the transported velocities are interpolated in a symmetry-
    ! preserving manner (the convective term has to be skew-symmetric in
    ! order to conserve energy).
    !
    ! Details can be found in:
    ! [1] Heinz Werner, Grobstruktursimulation der turbulenten Strömung
    !     über eine querliegende Rippe in einem Plattenkanal bei hoher
    !     Reynolds-Zahl, PhD Thesis, Technical University of Munich, 1991
    ! [2] Verstappen et al., SYMMETRY-PRESERVING DISCRETIZATIONS OF THE
    !     INCOMPRESSIBLE NAVIER-STOKES EQUATIONS, European Conference on
    !     Computational Fluid Dynamics, ECCOMAS CFD 2006
    SUBROUTINE tstle4_kon(kk, jj, ii, uo, vo, wo, u, v, w, ut, vt, wt, &
            dx, dy, dz, ddx, ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz, &
            nfro, nbac, nrgt, nlft, nbot, ntop, bu, bv, bw, itstep)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: uo(kk, jj, ii), vo(kk, jj, ii), &
            wo(kk, jj, ii)
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: ut(kk, jj, ii), vt(kk, jj, ii), &
            wt(kk, jj, ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: rdx(ii), rdy(jj), rdz(kk)
        REAL(realk), INTENT(in) :: rddx(ii), rddy(jj), rddz(kk)
        INTEGER, INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop
        REAL(realk), INTENT(in), OPTIONAL :: bu(kk, jj, ii), bv(kk, jj, ii), &
            bw(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: itstep

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: nbu, nfu, nrv, nbw, ntw, nlv
        REAL(realk) :: ax, ay, az
        REAL(realk) :: fw, fe, ft, fb, fn, fs
        REAL(realk) :: qw, qe, qt, qb, qn, qs

        ! Sanity check
        IF (PRESENT(bu) .NEQV. PRESENT(bv) .OR. &
                PRESENT(bu) .NEQV. PRESENT(bw)) THEN
            CALL errr(__FILE__, __LINE__)
        END IF


        DO i = 3-nfu, ii-3+nbu
            DO j = 3, jj-2
                DO k = 3, kk-2
                    ax = ddy(j)*ddz(k)
                    ay = dx(i)*ddz(k)
                    az = dx(i)*ddy(j)

                    uAdvectingE = 0.5 * ( u(k,j,i) + u(k,j,i+1) )
                    uAdvectingW = 0.5 * ( u(k,j,i-1) + u(k,j,i) )
                    vAdvectingN = 0.5 * ( v(k,j,i) + v(k,j,i+1) )
                    vAdvectingS = 0.5 * ( v(k,j-1,i) + v(k,j-1,i+1) )
                    wAdvectingT = 0.5 * ( w(k,j,i) + w(k,j,i+1) )
                    wAdvectingB = 0.5 * ( w(k-1,j,i) + w(k-1,j,i+1) )

                    CALL quick_interpolation_scheme(kk, jj, ii, u, &
                        uAdvectedE, uAdvectedW, uAdvectedN, uAdvectedS, uAdvectedT, uAdvectedB, &
                        uAdvectingE, uAdvectingW, vAdvectingN, vAdvectingS, wAdvectingT, wAdvectingB)
                    
                    IF ( isNearInterfaceiStag ) THEN
                        duo = - ( uAdvectedE * densityFluxiStag(k,j,i) - uAdvectedW * densityFluxiStag(k,j,i-1) ) * rdx(i) + &
                                ( uAdvectedN * densityFluxjStag(k,j,i) - uAdvectedS * densityFluxjStag(k,j-1,i) ) * rddy(j) + &
                                ( uAdvectedT * densityFluxkStag(k,j,i) - uAdvectedB * densityFluxkStag(k-1,j,i) ) * rddz(k) + &
                                u(k,j,i) * xDensityCompressionTerm(k,j,i)

                        uo(k,j,i) = 1 / densityFieldiStag(k,j,i) * ( uo(k,j,i) * densityFieldiStagOld(k,j,i) + dou )
                    ELSE
                        duo = - ( ( uAdvectedE * uAdvectingE - uAdvectedW * uAdvectingW ) * rdx(i) + &
                                  ( uAdvectedN * vAdvectingN - uAdvectedS * vAdvectingS ) * rddy(j) + &
                                  ( uAdvectedT * wAdvectingT - uAdvectedB * wAdvectingB ) * rddz(k) )

                        uo(k,j,i) = uo(k,j,i) + dou
                    END IF

                END DO
            END DO
        END DO

        ! nfu = 0
        ! nbu = 0
        ! nrv = 0
        ! nlv = 0
        ! nbw = 0
        ! ntw = 0

        ! ! CON = 7
        ! IF (nbac == 7) nbu = 1
        ! IF (nlft == 7) nlv = 1
        ! IF (ntop == 7) ntw = 1

        ! ! OP1 = 3
        ! IF (nfro == 3) nfu = 1
        ! IF (nbac == 3) nbu = 1
        ! IF (nrgt == 3) nrv = 1
        ! IF (nlft == 3) nlv = 1
        ! IF (nbot == 3) nbw = 1
        ! IF (ntop == 3) ntw = 1

        ! DO i = 3-nfu, ii-3+nbu
        !     DO j = 3, jj-2
        !         DO k = 3, kk-2
        !             ax = ddy(j)*ddz(k)
        !             ay = dx(i)*ddz(k)
        !             az = dx(i)*ddy(j)

        !             fe = ax*(ut(k, j, i) + (ut(k, j, i+1) - ut(k, j, i)) &
        !                 * 0.5*dx(i)/ddx(i+1))
        !             fw = ax*(ut(k, j, i-1) + (ut(k, j, i) - ut(k, j, i-1)) &
        !                 * 0.5*dx(i-1)/ddx(i))
        !             fn = ay*(vt(k, j, i) + vt(k, j, i+1))*0.5
        !             fs = ay*(vt(k, j-1, i) + vt(k, j-1, i+1))*0.5
        !             ft = az*(wt(k, j, i) + wt(k, j, i+1))*0.5
        !             fb = az*(wt(k-1, j, i) + wt(k-1, j, i+1))*0.5

        !             qe = 0.5*fe*(u(k, j, i) + u(k, j, i+1))
        !             qw = 0.5*fw*(u(k, j, i-1) + u(k, j, i))
        !             qn = 0.5*fn*(u(k, j, i) + u(k, j+1, i))
        !             qs = 0.5*fs*(u(k, j-1, i) + u(k, j, i))
        !             qt = 0.5*ft*(u(k, j, i) + u(k+1, j, i))
        !             qb = 0.5*fb*(u(k-1, j, i) + u(k, j, i))

        !             uo(k, j, i) = -(qe-qw+qn-qs+qt-qb)
        !         END DO

        !         IF (PRESENT(bu)) THEN
        !             DO k = 3, kk-2
        !                 uo(k, j, i) = bu(k, j, i)*uo(k, j, i)
        !             END DO
        !         ELSE
        !             DO k = 3, kk-2
        !                 uo(k, j, i) = rdx(i)*rddy(j)*rddz(k)*uo(k, j, i)
        !             END DO
        !         END IF
        !     END DO
        ! END DO

        ! DO i = 3, ii-2
        !     DO j = 3-nrv, jj-3+nlv
        !         DO k = 3, kk-2
        !             ax = dy(j)*ddz(k)
        !             ay = ddx(i)*ddz(k)
        !             az = ddx(i)*dy(j)

        !             fe = ax*(ut(k, j, i) + ut(k, j+1, i))*0.5
        !             fw = ax*(ut(k, j, i-1) + ut(k, j+1, i-1))*0.5
        !             fn = ay*(vt(k, j, i) + (vt(k, j+1, i) - vt(k, j, i)) &
        !                 * 0.5*dy(j)/ddy(j+1))
        !             fs = ay*(vt(k, j-1, i) + (vt(k, j, i) -vt(k, j-1, i)) &
        !                 * 0.5*dy(j-1)/ddy(j))
        !             ft = az*(wt(k, j, i) + wt(k, j+1, i))*0.5
        !             fb = az*(wt(k-1, j, i) + wt(k-1, j+1, i))*0.5

        !             qe = 0.5*fe*(v(k, j, i) + v(k, j, i+1))
        !             qw = 0.5*fw*(v(k, j, i-1) + v(k, j, i))
        !             qn = 0.5*fn*(v(k, j, i) + v(k, j+1, i))
        !             qs = 0.5*fs*(v(k, j-1, i) + v(k, j, i))
        !             qt = 0.5*ft*(v(k, j, i) + v(k+1, j, i))
        !             qb = 0.5*fb*(v(k-1, j, i) + v(k, j, i))

        !             vo(k, j, i) = -(qe-qw+qn-qs+qt-qb)
        !         END DO

        !         IF (PRESENT(bv)) THEN
        !             DO k = 3, kk-2
        !                 vo(k, j, i) = bv(k, j, i)*vo(k, j, i)
        !             END DO
        !         ELSE
        !             DO k = 3, kk-2
        !                 vo(k, j, i) = rddx(i)*rdy(j)*rddz(k)*vo(k, j, i)
        !             END DO
        !         END IF
        !     END DO
        ! END DO

        ! DO i = 3, ii-2
        !     DO j = 3, jj-2
        !         DO k = 3-nbw, kk-3+ntw
        !             ax = ddy(j)*dz(k)
        !             ay = ddx(i)*dz(k)
        !             az = ddx(i)*ddy(j)

        !             fe = ax*(ut(k, j, i) + ut(k+1, j, i))*0.5
        !             fw = ax*(ut(k, j, i-1)+ ut(k+1, j, i-1))*0.5
        !             fn = ay*(vt(k, j, i) + vt(k+1, j, i))*0.5
        !             fs = ay*(vt(k, j-1, i)+ vt(k+1, j-1, i))*0.5
        !             ft = az*(wt(k, j, i) + (wt(k+1, j, i) - wt(k, j, i)) &
        !                 * 0.5*dz(k)/ddz(k+1))
        !             fb = az*(wt(k-1, j, i) + (wt(k, j, i) - wt(k-1, j, i)) &
        !                 * 0.5*dz(k-1)/ddz(k))

        !             qe = 0.5*fe*(w(k, j, i) + w(k, j, i+1))
        !             qw = 0.5*fw*(w(k, j, i-1) + w(k, j, i))
        !             qn = 0.5*fn*(w(k, j, i) + w(k, j+1, i))
        !             qs = 0.5*fs*(w(k, j-1, i) + w(k, j, i))
        !             qt = 0.5*ft*(w(k, j, i) + w(k+1, j, i))
        !             qb = 0.5*fb*(w(k-1, j, i) + w(k, j, i))

        !             wo(k, j, i) = -(qe-qw+qn-qs+qt-qb)
        !         END DO

        !         IF (PRESENT(bw)) THEN
        !             DO k = 3-nbw, kk-3+ntw
        !                 wo(k, j, i) = bw(k, j, i)*wo(k, j, i)
        !             END DO
        !         ELSE
        !             DO k = 3-nbw, kk-3+ntw
        !                 wo(k, j, i) = rddx(i)*rddy(j)*rdz(k)*wo(k, j, i)
        !             END DO
        !         END IF
        !     END DO
        ! END DO
    END SUBROUTINE tstle4_kon


    SUBROUTINE tstle4_diff(kk, jj, ii, uo, vo, wo, u, v, w, g, &
            dx, dy, dz, ddx, ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz, &
            nfro, nbac, nrgt, nlft, nbot, ntop, densityFieldiStag, densityFieldjStag, densityFieldkStag)
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

        CALL swcle3d(kk, jj, ii, uo, vo, wo, u, v, w, &
            ddx, ddy, ddz, nfro, nbac, nrgt, nlft, nbot, ntop)

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
                    !                                    ---------------------------------------outer derivatives----------------------------------------
                    duo = - 1/densityFieldiStag(k,j,i) * ( ( tauxxe - tauxxw ) * rdx(i) + ( tauyxn - tauyxs ) * rddy(j) + ( tauzxt - tauzxb ) * rddz(k) )

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
                    !                                    ---------------------------------------outer derivatives----------------------------------------
                    dvo = - 1/densityFieldjStag(k,j,i) * ( ( tauxye - tauxyw ) * rddx(i) + ( tauyyn - tauyys ) * rdy(j) + ( tauzyt - tauzyb ) * rddz(k) )

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
                    !                                    ---------------------------------------outer derivatives----------------------------------------
                    dwo = - 1/densityFieldkStag(k,j,i) * ( ( tauxze - tauxzw ) * rddx(i) + ( tauyzn - tauyzs ) * rddy(j) + ( tauzzt - tauzzb ) * rdz(k) )

                    ! Addition
                    wo(k, j, i) = wo(k, j, i) + dwo
                END DO
            END DO
        END DO
    END SUBROUTINE tstle4_diff


    ! SUBROUTINE tstle4_gradp(kk, jj, ii, uo, vo, wo, p, dx, dy, dz, &
    !         nfro, nbac, nrgt, nlft, nbot, ntop, igrid)
    !     ! Subroutine arguments
    !     INTEGER(intk), INTENT(in) :: kk, jj, ii
    !     REAL(realk), INTENT(inout) :: uo(kk, jj, ii), vo(kk, jj, ii), &
    !         wo(kk, jj, ii)
    !     REAL(realk), INTENT(in) :: p(kk, jj, ii)
    !     REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
    !     INTEGER, INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop
    !     INTEGER, INTENT(in) :: igrid

    !     ! Local variables
    !     INTEGER(intk) :: k, j, i
    !     INTEGER(intk) :: nbu, nfu, nrv, nbw, ntw, nlv
    !     INTEGER(intk) :: gradpflag
    !     REAL(realk) :: gpx, gpy, gpz

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

    !     CALL get_gradpxflag(gradpflag, igrid)
    !     gpx = gradp(1)*gradpflag
    !     gpy = gradp(2)*gradpflag
    !     gpz = gradp(3)*gradpflag

    !     DO i = 3-nfu, ii-3+nbu
    !         DO j = 3, jj-2
    !             DO k = 3, kk-2
    !                 uo(k, j, i) = uo(k, j, i) - 1.0/(rho*dx(i)) &
    !                     *(p(k, j, i+1) - p(k, j, i) + gpx*dx(i))
    !             END DO
    !         END DO
    !     END DO

    !     DO i = 3, ii-2
    !         DO j = 3-nrv, jj-3+nlv
    !             DO k = 3, kk-2
    !                 vo(k, j, i) = vo(k, j, i) - 1.0/(rho*dy(j)) &
    !                     *(p(k, j+1, i) - p(k, j, i) + gpy*dy(j))
    !             END DO
    !         END DO
    !     END DO

    !     DO i = 3, ii-2
    !         DO j = 3, jj-2
    !             DO k = 3-nbw, kk-3+ntw
    !                 wo(k, j, i) = wo(k, j, i) - 1.0/(rho*dz(k)) &
    !                     *(p(k+1, j, i) - p(k, j, i) + gpz*dz(k))
    !             END DO
    !         END DO
    !     END DO
    ! END SUBROUTINE tstle4_gradp


    ! SUBROUTINE tstle4_par(kk, jj, ii, uo, vo, wo, u, v, w, ut, vt, wt, &
    !         dx, dy, dz, ddx, ddy, ddz, rdx, rdy, rdz, rddx, rddy, rddz, &
    !         nfro, nbac, nrgt, nlft, nbot, ntop)
    !     ! Subroutine arguments
    !     INTEGER(intk), INTENT(in) :: kk, jj, ii
    !     REAL(realk), INTENT(inout) :: uo(kk, jj, ii), vo(kk, jj, ii), &
    !         wo(kk, jj, ii)
    !     REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
    !     REAL(realk), INTENT(in) :: ut(kk, jj, ii), vt(kk, jj, ii), &
    !         wt(kk, jj, ii)
    !     REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
    !     REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
    !     REAL(realk), INTENT(in) :: rdx(ii), rdy(jj), rdz(kk)
    !     REAL(realk), INTENT(in) :: rddx(ii), rddy(jj), rddz(kk)
    !     INTEGER, INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop

    !     ! Local variables
    !     INTEGER(intk) :: k, j, i
    !     REAL(realk) :: fkdtu, fkdtv, fkdtw
    !     REAL(realk) :: qkubadd, qkusadd, qkvbadd
    !     REAL(realk) :: qkvwadd, qkwsadd, qkwwadd

    !     REAL(realk) :: qkut, qkub, qkun, qkus, qkvw, qkve, qkvt, qkvb, &
    !         qkww, qkwe, qkwn, qkws, &
    !         fut, fub, fun, fus, auy, auz, &
    !         fvw, fve, fvt, fvb, avx, avz, &
    !         fww, fwe, fwn, fws, awx, awy
    !     REAL(realk) :: dxi, ddxi, dyj, ddyj, dzk, ddzk, rdzk, rddzk
    !     REAL(realk), PARAMETER :: wkon = 1.0

    !     ! Temporary storage
    !     ! TODO: Can this be replaced with a single tmp(:, :, :) array??
    !     ! I see no places where three arrays are needed at the same time...
    !     REAL(realk), ALLOCATABLE :: wcu(:, :, :), wcv(:, :, :), wcw(:, :, :)

    !     ALLOCATE(wcu(kk, jj, ii))
    !     ALLOCATE(wcv(kk, jj, ii))
    !     ALLOCATE(wcw(kk, jj, ii))

    !     ! Upwind in vorletzter schicht bei PAR-randbedingung
    !     IF (nfro == 8) THEN
    !         i = 4
    !         DO j = 3, jj-2
    !             dyj = dy(j)
    !             ddyj = ddy(j)
    !             fkdtu = -1.0*rddy(j)*rdx(i)*wkon
    !             fkdtv = -1.0*rddx(i)*rdy(j)*wkon
    !             fkdtw = -1.0*rddx(i)*rddy(j)*wkon

    !             DO k = 3, kk-2
    !                 dzk = dz(k)
    !                 rdzk = rdz(k)
    !                 ddzk = ddz(k)
    !                 rddzk = rddz(k)
    !                 avx = dyj*ddzk
    !                 awx = ddyj*dzk
    !                 fvw = avx*(ut(k, j, i-1) + ut(k, j+1, i-1))*0.5
    !                 qkvwadd = 0.5*(fvw-ABS(fvw)) &
    !                     *(0.5*v(k, j, i)-1.5*v(k, j, i-1)+v(k, j, i-2))*0.5*0.5

    !                 vo(k, j, i) = vo(k, j, i) + fkdtv * rddzk * (-qkvwadd)
    !                 vo(k, j, i-1) = vo(k, j, i-1) + fkdtv * rddzk * (+qkvwadd)

    !                 fww = awx*(ut(k, j, i-1)+ ut(k+1, j, i-1))*0.5
    !                 qkwwadd = 0.5 *(fww-ABS(fww)) &
    !                     *(0.5*w(k, j, i)-1.5*w(k, j, i-1)+w(k, j, i-2))*0.5*0.5

    !                 wo(k, j, i) = wo(k, j, i) + fkdtw * rdzk * (-qkwwadd)
    !                 wo(k, j, i-1) = wo(k, j, i-1) + fkdtw * rdzk * (+qkwwadd)
    !             END DO
    !         END DO
    !     END IF

    !     IF (nbac == 8) THEN
    !         i = ii-2
    !         DO j = 3, jj-2
    !             dyj = dy(j)
    !             ddyj = ddy(j)
    !             fkdtu = -1.0*rddy(j)*rdx(i)*wkon
    !             fkdtv = -1.0*rddx(i)*rdy(j)*wkon
    !             fkdtw = -1.0*rddx(i)*rddy(j)*wkon

    !             DO k = 3, kk-2
    !                 dzk = dz(k)
    !                 rdzk = rdz(k)
    !                 ddzk = ddz(k)
    !                 rddzk = rddz(k)
    !                 avx = dyj*ddzk
    !                 awx = ddyj*dzk
    !                 fvw = avx*(ut(k, j, i-1) + ut(k, j+1, i-1))*0.5

    !                 qkvwadd = 0.5*(fvw+ABS(fvw)) &
    !                     *(0.5*v(k, j, i-1)-1.5*v(k, j, i)+v(k, j, i+1))*0.5*0.5

    !                 vo(k, j, i) = vo(k, j, i) + fkdtv * rddzk * (-qkvwadd)
    !                 vo(k, j, i-1) = vo(k, j, i-1) + fkdtv * rddzk * (+qkvwadd)

    !                 fww = awx*(ut(k, j, i-1)+ ut(k+1, j, i-1))*0.5
    !                 qkwwadd = 0.5*(fww+ABS(fww)) &
    !                     *(0.5*w(k, j, i-1)-1.5*w(k, j, i)+w(k, j, i+1))*0.5*0.5

    !                 wo(k, j, i) = wo(k, j, i) + fkdtw * rdzk * (-qkwwadd)
    !                 wo(k, j, i-1) = wo(k, j, i-1) + fkdtw * rdzk * (+qkwwadd)
    !             END DO
    !         END DO
    !     END IF

    !     IF (nrgt == 8) THEN
    !         DO i = 3, ii-2
    !             j = 4
    !             dxi = dx(i)
    !             ddxi = ddx(i)
    !             fkdtu = -1.0*rddy(j)*rdx(i)*wkon
    !             fkdtv = -1.0*rddx(i)*rdy(j)*wkon
    !             fkdtw = -1.0*rddx(i)*rddy(j)*wkon

    !             DO k = 3, kk-2
    !                 dzk = dz(k)
    !                 rdzk = rdz(k)
    !                 ddzk = ddz(k)
    !                 rddzk = rddz(k)
    !                 auy = dxi*ddzk
    !                 awy = ddxi*dzk

    !                 fus = auy*(vt(k, j-1, i) + vt(k, j-1, i+1))*0.5
    !                 qkusadd = 0.5 *(fus-ABS(fus)) &
    !                     *(0.5*u(k, j, i)-1.5*u(k, j-1, i)+u(k, j-2, i))*0.5*0.5

    !                 uo(k, j, i) = uo(k, j, i) + fkdtu * rddzk * (-qkusadd)
    !                 uo(k, j-1, i) = uo(k, j-1, i) + fkdtu * rddzk * (+qkusadd)

    !                 fws = awy*(vt(k, j-1, i)+ vt(k+1, j-1, i))*0.5
    !                 qkwsadd = 0.5 *(fws-ABS(fws)) &
    !                     *(0.5*w(k, j, i)-1.5*w(k, j-1, i)+w(k, j-2, i))*0.5*0.5

    !                 wo(k, j, i) = wo(k, j, i) + fkdtw * rdzk * (-qkwsadd)
    !                 wo(k, j-1, i) = wo(k, j-1, i) + fkdtw * rdzk * (+qkwsadd)
    !             END DO
    !         END DO
    !     END IF

    !     IF (nlft == 8) THEN
    !         DO i = 3, ii-2
    !             j = jj-2
    !             dxi = dx(i)
    !             ddxi = ddx(i)
    !             fkdtu = -1.0*rddy(j)*rdx(i)*wkon
    !             fkdtv = -1.0*rddx(i)*rdy(j)*wkon
    !             fkdtw = -1.0*rddx(i)*rddy(j)*wkon

    !             DO k = 3, kk-2
    !                 dzk = dz(k)
    !                 rdzk = rdz(k)
    !                 ddzk = ddz(k)
    !                 rddzk = rddz(k)
    !                 auy = dxi*ddzk
    !                 awy = ddxi*dzk

    !                 fus = auy*(vt(k, j-1, i) + vt(k, j-1, i+1))*0.5
    !                 qkusadd = 0.5*(fus+ABS(fus)) &
    !                     *(0.5*u(k, j-1, i)-1.5*u(k, j, i)+u(k, j+1, i))*0.5*0.5

    !                 uo(k, j, i) = uo(k, j, i) + fkdtu * rddzk * (-qkusadd)
    !                 uo(k, j-1, i) = uo(k, j-1, i) + fkdtu * rddzk * (+qkusadd)
    !                 fws = awy*(vt(k, j-1, i) + vt(k+1, j-1, i))*0.5
    !                 qkwsadd = 0.5*(fws+ABS(fws)) &
    !                     *(0.5*w(k, j-1, i)-1.5*w(k, j, i)+w(k, j+1, i))*0.5*0.5

    !                 wo(k, j, i) = wo(k, j, i) + fkdtw * rdzk * (-qkwsadd)
    !                 wo(k, j-1, i) = wo(k, j-1, i) + fkdtw * rdzk * (+qkwsadd)
    !             END DO
    !         END DO
    !     END IF

    !     IF (nbot == 8) THEN
    !         DO i = 3, ii-2
    !             dxi = dx(i)
    !             ddxi = ddx(i)

    !             DO j = 3, jj-2
    !                 k = 4

    !                 fkdtu = -1.0*rddy(j)*rdx(i)*wkon
    !                 fkdtv = -1.0*rddx(i)*rdy(j)*wkon
    !                 fkdtw = -1.0*rddx(i)*rddy(j)*wkon

    !                 dyj = dy(j)
    !                 ddyj = ddy(j)
    !                 auz = dxi*ddyj
    !                 avz = ddxi*dyj
    !                 rdzk = rdz(k)
    !                 rddzk = rddz(k)

    !                 fub = auz*(wt(k-1, j, i)+ wt(k-1, j, i+1))*0.5
    !                 qkubadd = 0.5*(fub-ABS(fub)) &
    !                     *(0.5*u(k, j, i)-1.5*u(k-1, j, i)+u(k-2, j, i))*0.5*0.5
    !                 uo(k, j, i) = uo(k, j, i) + fkdtu * rddzk * (-qkubadd)
    !                 uo(k-1, j, i) = uo(k-1, j, i) + fkdtu * rddzk * (+qkubadd)

    !                 fvb = avz*(wt(k-1, j, i) + wt(k-1, j+1, i))*0.5
    !                 qkvbadd = 0.5*(fvb-ABS(fvb)) &
    !                     *(0.5*v(k, j, i)-1.5*v(k-1, j, i)+v(k-2, j, i))*0.5*0.5

    !                 vo(k, j, i) = vo(k, j, i) + fkdtv * rddzk * (-qkvbadd)
    !                 vo(k-1, j, i) = vo(k-1, j, i) + fkdtv * rddzk * (+qkvbadd)
    !             END DO
    !         END DO
    !     END IF

    !     IF (ntop == 8) THEN
    !         DO i = 3, ii-2
    !             dxi = dx(i)
    !             ddxi = ddx(i)

    !             DO j = 3, jj-2
    !                 k = kk-2

    !                 fkdtu = -1.0*rddy(j)*rdx(i)*wkon
    !                 fkdtv = -1.0*rddx(i)*rdy(j)*wkon
    !                 fkdtw = -1.0*rddx(i)*rddy(j)*wkon

    !                 dyj = dy(j)
    !                 ddyj = ddy(j)
    !                 auz = dxi*ddyj
    !                 avz = ddxi* dyj

    !                 rdzk = rdz(k)
    !                 rddzk = rddz(k)

    !                 fub = auz*(wt(k-1, j, i) + wt(k-1, j, i+1))*0.5
    !                 qkubadd = 0.5 *(fub+ABS(fub)) &
    !                     *(0.5*u(k-1, j, i)-1.5*u(k, j, i)+u(k+1, j, i))*0.5*0.5
    !                 uo(k, j, i) = uo(k, j, i) + fkdtu * rddzk * (-qkubadd)
    !                 uo(k-1, j, i) = uo(k-1, j, i) + fkdtu * rddzk * (+qkubadd)

    !                 fvb = avz*(wt(k-1, j, i)+ wt(k-1, j+1, i))*0.5
    !                 qkvbadd = 0.5 *(fvb+ABS(fvb)) &
    !                     *(0.5*v(k-1, j, i)-1.5*v(k, j, i)+v(k+1, j, i))*0.5*0.5

    !                 vo(k, j, i) = vo(k, j, i) + fkdtv * rddzk * (-qkvbadd)
    !                 vo(k-1, j, i) = vo(k-1, j, i) + fkdtv * rddzk * (+qkvbadd)
    !             END DO
    !         END DO

    !     END IF

    !     ! PAR-RB Impulserhaltend BACK
    !     IF (nbac == 8) THEN
    !         ! W-Impulszelle
    !         ! STANDART-BERECHNUNG DES KONVEKTIVEN FLUSSES
    !         i = ii-2
    !         DO j = 3, jj-2
    !             ddyj = ddy(j)
    !             DO k = 2, kk-2
    !                 dzk = dz(k)
    !                 awx = ddyj*dzk
    !                 fwe = awx*(ut(k, j, i) + ut(k+1, j, i))*0.5
    !                 qkwe = 0.5*fwe*(w(k, j, i) + w(k, j, i+1))
    !                 wcw(k, j, i+1) = qkwe
    !             END DO
    !         END DO

    !         ! NEUE BERECHNUNG DES KONVEKTIVEN FLUSSES
    !         ! FUER JEDE GROBGITTERGESCHWINDIGKEITT
    !         i = ii-2
    !         DO j = 3, jj-2, 2
    !             ddyj = ddy(j)
    !             k = 2
    !             dzk = dz(k)
    !             awx = ddyj*dzk
    !             fwe = awx*(2*ut(k, j, i) + 2*ut(k, j+1, i) &
    !                 + ut(k+1, j, i) + ut(k+1, j+1, i) &
    !                 + ut(k+2, j, i) + ut(k+2, j+1, i))*0.125
    !             qkwe = 0.5*fwe*(w(k, j, i) + w(k, j, i+1))
    !             wcw(k, j, i) = qkwe
    !             qkwe = 0.5*fwe*(w(k, j+1, i) + w(k, j+1, i+1))
    !             wcw(k, j+1, i) = qkwe
    !             DO k = 4, kk-4, 2
    !                 dzk = dz(k)
    !                 awx = ddyj*dzk
    !                 fwe = awx*(ut(k-1, j, i) + ut(k-1, j+1, i) &
    !                     + ut(k, j, i) + ut(k, j+1, i) &
    !                     + ut(k+1, j, i) + ut(k+1, j+1, i) &
    !                     + ut(k+2, j, i)+ ut(k+2, j+1, i))*0.125
    !                 qkwe = 0.5*fwe*(w(k, j, i) + w(k, j, i+1))
    !                 wcw(k, j, i) = qkwe
    !                 qkwe = 0.5*fwe*(w(k, j+1, i) + w(k, j+1, i+1))
    !                 wcw(k, j+1, i) = qkwe
    !             END DO
    !             k = kk-2
    !             dzk = dz(k)
    !             awx = ddyj*dzk
    !             fwe = awx*(ut(k-1, j, i) + ut(k-1, j+1, i) &
    !                 + ut(k, j, i) + ut(k, j+1, i)&
    !                 + 2*ut(k+1, j, i) + 2*ut(k+1, j+1, i)) *0.125
    !             qkwe = 0.5*fwe*(w(k, j, i) + w(k, j, i+1))
    !             wcw(k, j, i) = qkwe
    !             qkwe = 0.5*fwe*(w(k, j+1, i) + w(k, j+1, i+1))
    !             wcw(k, j+1, i) = qkwe
    !         END DO

    !         ! VERTEILUNG
    !         i = ii-2
    !         DO j = 3, jj-2
    !             DO k = 2, kk-4, 2
    !                 wcw(k+1, j, i) = 0.5*(wcw(k, j, i) + wcw(k+2, j, i))
    !             END DO
    !         END DO

    !         ! AUF WO SCHREIBEN
    !         DO j = 3, jj-2
    !             fkdtw = -1.0*rddx(i)*rddy(j)*wkon
    !             DO k = 2, kk-2
    !                 rdzk = rdz(k)
    !                 wo(k, j, i) = wo(k, j, i) &
    !                     + fkdtw*rdzk*(-wcw(k, j, i+1) + wcw(k, j, i))
    !             END DO
    !         END DO

    !         ! V-Impulszelle
    !         ! STANDARD-BERECHNUNG DES KONVEKTIVEN FLUSSES
    !         i = ii-2
    !         DO j = 2, jj-2
    !             dyj = dy(j)
    !             DO k = 3, kk-2
    !                 ddzk = ddz(k)
    !                 avx = ddzk*dyj
    !                 fve = avx *(ut(k, j, i) + ut(k, j+1, i))*0.5
    !                 qkve = 0.5*fve*(v(k, j, i) + v(k, j, i+1))
    !                 wcv(k, j, i+1) = qkve
    !             END DO
    !         END DO

    !         ! NEUE BERECHNUNG DES KONVEKTIVEN FLUSSES
    !         ! FUER JEDE GROBGITTERGESCHWINDIGKEIT
    !         i = ii-2

    !         ! YM-RAND
    !         j = 2
    !         dyj = dy(j)
    !         DO k = 3, kk-2, 2
    !             ddzk = ddz(k)
    !             avx = ddzk*dyj
    !             fve = avx*(2*ut(k, j, i) + 2*ut(k+1, j, i) &
    !                 + ut(k, j+1, i) + ut(k+1, j+1, i) &
    !                 + ut(k, j+2, i) + ut(k+1, j+2, i))*0.125
    !             qkve = 0.5*fve*(v(k, j, i) + v(k, j, i+1))
    !             wcv(k, j, i) = qkve
    !             qkve = 0.5*fve*(v(k+1, j, i) + v(k+1, j, i+1))
    !             wcv(k+1, j, i) = qkve
    !         END DO

    !         ! IM GEBIET
    !         DO j = 4, jj-4, 2
    !             dyj = dy(j)
    !             DO k = 3, kk-2, 2
    !                 ddzk = ddz(k)
    !                 avx = ddzk*dyj
    !                 fve = avx*(ut(k, j-1, i) + ut(k+1, j-1, i) &
    !                     + ut(k, j, i) + ut(k+1, j, i) &
    !                     + ut(k, j+1, i) + ut(k+1, j+1, i) &
    !                     + ut(k, j+2, i) + ut(k+1, j+2, i))*0.125
    !                 qkve = 0.5*fve*(v(k, j, i) + v(k, j, i+1))
    !                 wcv(k, j, i) = qkve
    !                 qkve = 0.5*fve*(v(k+1, j, i) + v(k+1, j, i+1))
    !                 wcv(k+1, j, i) = qkve
    !             END DO
    !         END DO

    !         ! YP-RAND
    !         j = jj-2
    !         dyj = dy(j)
    !         DO k = 3, kk-2, 2
    !             ddzk = ddz(k)
    !             avx = ddzk*dyj
    !             fve = avx*(ut(k, j-1, i) + ut(k+1, j-1, i) &
    !                 + ut(k, j, i) + ut(k+1, j, i) &
    !                 + 2*ut(k, j+1, i) + 2*ut(k+1, j+1, i))*0.125
    !             qkve = 0.5*fve*(v(k, j, i) + v(k, j, i+1))
    !             wcv(k, j, i) = qkve
    !             qkve = 0.5*fve*(v(k+1, j, i) + v(k+1, j, i+1))
    !             wcv(k+1, j, i) = qkve
    !         END DO

    !         ! VERTEILUNG
    !         i = ii-2
    !         DO j = 2, jj-4, 2
    !             DO k = 3, kk-2
    !                 wcv(k, j+1, i) = 0.5*(wcv(k, j, i) + wcv(k, j+2, i))
    !             END DO
    !         END DO

    !         ! AUF VO SCHREIBEN
    !         DO j = 2, jj-2
    !             fkdtv = -1.0*rddx(i)* rdy(j)*wkon
    !             DO k = 3, kk-2
    !                 rddzk = rddz(k)
    !                 vo(k, j, i) = vo(k, j, i) &
    !                     + fkdtv*rddzk*(wcv(k, j, i) - wcv(k, j, i+1))
    !             END DO
    !         END DO
    !     END IF

    !     ! PAR-RB Impulserhaltend FRONT
    !     IF (nfro == 8) THEN
    !         ! W-Impulszelle
    !         ! STANDART-BERECHNUNG DES KONVEKTIVEN FLUSSES
    !         i = 3
    !         DO j = 3, jj-2
    !             ddyj = ddy(j)
    !             DO k = 2, kk-2
    !                 dzk = dz(k)
    !                 awx = ddyj*dzk
    !                 fww = awx*(ut(k, j, i-1) + ut(k+1, j, i-1))*0.5
    !                 qkww = 0.5*fww*(w(k, j, i) + w(k, j, i-1))
    !                 wcw(k, j, i-1) = qkww
    !             END DO
    !         END DO

    !         ! NEUE BERECHNUNG DES KONVEKTIVEN FLUSSES
    !         ! FUER JEDE GROBGITTERGESCHWINDIGKEIT
    !         i = 3
    !         DO j = 3, jj-2, 2
    !             ddyj = ddy(j)
    !             k = 2
    !             dzk = dz(k)
    !             awx = ddyj*dzk
    !             fww = awx*(2*ut(k, j, i-1) + 2*ut(k, j+1, i-1) &
    !                 + ut(k+1, j, i-1) + ut(k+1, j+1, i-1) &
    !                 + ut(k+2, j, i-1) + ut(k+2, j+1, i-1))*0.125
    !             qkww = 0.5*fww*(w(k, j, i) + w(k, j, i-1))
    !             wcw(k, j, i) = qkww
    !             qkww = 0.5*fww*(w(k, j+1, i) + w(k, j+1, i-1))
    !             wcw(k, j+1, i) = qkww
    !             DO k = 4, kk-4, 2
    !                 dzk = dz(k)
    !                 awx = ddyj*dzk
    !                 fww = awx*(ut(k-1, j, i-1) + ut(k-1, j+1, i-1) &
    !                     + ut(k, j, i-1) + ut(k, j+1, i-1) &
    !                     + ut(k+1, j, i-1) + ut(k+1, j+1, i-1) &
    !                     + ut(k+2, j, i-1) + ut(k+2, j+1, i-1))*0.125
    !                 qkww = 0.5*fww*(w(k, j, i) + w(k, j, i-1))
    !                 wcw(k, j, i) = qkww
    !                 qkww = 0.5*fww*(w(k, j+1, i) + w(k, j+1, i-1))
    !                 wcw(k, j+1, i) = qkww
    !             END DO
    !             k = kk-2
    !             dzk = dz(k)
    !             awx = ddyj*dzk
    !             fww = awx*(ut(k-1, j, i-1) + ut(k-1, j+1, i-1) &
    !                 + ut(k, j, i-1) + ut(k, j+1, i-1) &
    !                 + 2*ut(k+1, j, i-1) + 2*ut(k+1, j+1, i-1))*0.125
    !             qkww = 0.5*fww*(w(k, j, i) + w(k, j, i-1))
    !             wcw(k, j, i) = qkww
    !             qkww = 0.5*fww*(w(k, j+1, i) + w(k, j+1, i-1))
    !             wcw(k, j+1, i) = qkww
    !         END DO

    !         ! VERTEILUNG
    !         i = 3
    !         DO j = 3, jj-2
    !             DO k = 2, kk-4, 2
    !                 wcw(k+1, j, i) = 0.5*(wcw(k, j, i) + wcw(k+2, j, i))
    !             END DO
    !         END DO

    !         ! AUF WO SCHREIBEN
    !         DO j = 3, jj-2
    !             fkdtw = -1.0*rddx(i)*rddy(j)*wkon
    !             DO k = 2, kk-2
    !                 rdzk = rdz(k)
    !                 wo(k, j, i) = wo(k, j, i) &
    !                     + fkdtw*rdzk*(wcw(k, j, i-1) - wcw(k, j, i))
    !             END DO
    !         END DO

    !         ! V-Impulszelle
    !         ! STANDARD-BERECHNUNG DES KONVEKTIVEN FLUSSES
    !         i = 3
    !         DO j = 2, jj-2
    !             dyj = dy(j)
    !             DO k = 3, kk-2
    !                 ddzk = ddz(k)
    !                 avx = ddzk*dyj
    !                 fvw = avx *(ut(k, j, i-1) + ut(k, j+1, i-1))*0.5
    !                 qkvw = 0.5*fvw*(v(k, j, i) + v(k, j, i-1))
    !                 wcv(k, j, i-1) = qkvw
    !             END DO
    !         END DO

    !         ! NEUE BERECHNUNG DES KONVEKTIVEN FLUSSES
    !         ! FUER JEDE GROBGITTERGESCHWINDIGKEIT
    !         i = 3

    !         ! YM-RAND
    !         j = 2
    !         dyj = dy(j)
    !         DO k = 3, kk-2, 2
    !             ddzk = ddz(k)
    !             avx = ddzk*dyj
    !             fvw = avx*(2*ut(k, j, i-1) + 2*ut(k+1, j, i-1) &
    !                 + ut(k, j+1, i-1) + ut(k+1, j+1, i-1) &
    !                 + ut(k, j+2, i-1) + ut(k+1, j+2, i-1))*0.125
    !             qkvw = 0.5*fvw*(v(k, j, i) + v(k, j, i-1))
    !             wcv(k, j, i) = qkvw
    !             qkvw = 0.5*fvw*(v(k+1, j, i) + v(k+1, j, i-1))
    !             wcv(k+1, j, i) = qkvw
    !         END DO

    !         ! IM GEBIET
    !         DO j = 4, jj-4, 2
    !             dyj = dy(j)
    !             DO k = 3, kk-2, 2
    !                 ddzk = ddz(k)
    !                 avx = ddzk*dyj
    !                 fvw = avx*(ut(k, j-1, i-1) + ut(k+1, j-1, i-1) &
    !                     + ut(k, j, i-1) + ut(k+1, j, i-1) &
    !                     + ut(k, j+1, i-1) + ut(k+1, j+1, i-1) &
    !                     + ut(k, j+2, i-1) + ut(k+1, j+2, i-1))*0.125
    !                 qkvw = 0.5*fvw*(v(k, j, i) + v(k, j, i-1))
    !                 wcv(k, j, i) = qkvw
    !                 qkvw = 0.5*fvw*(v(k+1, j, i) + v(k+1, j, i-1))
    !                 wcv(k+1, j, i) = qkvw
    !             END DO
    !         END DO

    !         ! YP-RAND
    !         j = jj-2
    !         dyj = dy(j)
    !         DO k = 3, kk-2, 2
    !             ddzk = ddz(k)
    !             avx = ddzk*dyj
    !             fvw = avx*(ut(k, j-1, i-1) + ut(k+1, j-1, i-1) &
    !                 + ut(k, j, i-1) + ut(k+1, j, i-1) &
    !                 + 2*ut(k, j+1, i-1) + 2*ut(k+1, j+1, i-1))*0.125
    !             qkvw = 0.5*fvw*(v(k, j, i) + v(k, j, i-1))
    !             wcv(k, j, i) = qkvw
    !             qkvw = 0.5*fvw*(v(k+1, j, i) + v(k+1, j, i-1))
    !             wcv(k+1, j, i) = qkvw
    !         END DO

    !         ! VERTEILUNG
    !         i = 3
    !         DO j = 2, jj-4, 2
    !             DO k = 3, kk-2
    !                 wcv(k, j+1, i) = 0.5*(wcv(k, j, i) + wcv(k, j+2, i))
    !             END DO
    !         END DO

    !         ! AUF VO SCHREIBEN
    !         DO j = 2, jj-2
    !             fkdtv = -1.0*rddx(i)*rdy(j)*wkon
    !             DO k = 3, kk-2
    !                 rddzk = rddz(k)
    !                 vo(k, j, i) = vo(k, j, i) &
    !                     + fkdtv*rddzk*(-wcv(k, j, i) + wcv(k, j, i-1))
    !             END DO
    !         END DO
    !     END IF

    !     ! PAR-RB Impulserhaltend TOP
    !     IF (ntop == 8) THEN

    !         ! U-Impulszelle
    !         ! STANDARD-BERECHNUNG DES KONVEKTIVEN FLUSSES
    !         k = kk-2
    !         DO i = 2, ii-2
    !             dxi = dx(i)
    !             DO j = 3, jj-2
    !                 ddyj = ddy(j)
    !                 auz = dxi*ddyj
    !                 fut = auz*(wt(k, j, i) + wt(k, j, i+1))*0.5
    !                 qkut = 0.5*fut*(u(k, j, i) + u(k+1, j, i))
    !                 wcu(k+1, j, i) = qkut
    !             END DO
    !         END DO

    !         ! NEUE BERECHNUNG DES KONVEKTIVEN FLUSSES
    !         ! FUER JEDE GROBGITTERGESCHWINDIGKEIT
    !         k = kk-2

    !         ! XM-RAND
    !         i = 2
    !         dxi = dx(i)
    !         DO j = 3, jj-2, 2
    !             ddyj = ddy(j)
    !             auz = dxi*ddyj
    !             fut = auz*(2*wt(k, j, i) + 2*wt(k, j+1, i) &
    !                 + wt(k, j, i+1) + wt(k, j+1, i+1) &
    !                 + wt(k, j, i+2) + wt(k, j+1, i+2))*0.125
    !             qkut = 0.5*fut*(u(k, j, i) + u(k+1, j, i))
    !             wcu(k, j, i) = qkut
    !             qkut = 0.5*fut*(u(k, j+1, i) + u(k+1, j+1, i))
    !             wcu(k, j+1, i) = qkut
    !         END DO

    !         ! IM GEBIET
    !         DO i = 4, ii-4, 2
    !             dxi = dx(i)
    !             DO j = 3, jj-2, 2
    !                 ddyj = ddy(j)
    !                 auz = dxi*ddyj
    !                 fut = auz*(wt(k, j, i-1) + wt(k, j+1, i-1) &
    !                     + wt(k, j, i) + wt(k, j+1, i) &
    !                     + wt(k, j, i+1) + wt(k, j+1, i+1) &
    !                     + wt(k, j, i+2) + wt(k, j+1, i+2))*0.125
    !                 qkut = 0.5*fut*(u(k, j, i) + u(k+1, j, i))
    !                 wcu(k, j, i) = qkut
    !                 qkut = 0.5*fut*(u(k, j+1, i) + u(k+1, j+1, i))
    !                 wcu(k, j+1, i) = qkut
    !             END DO
    !         END DO

    !         ! XP-RAND
    !         i = ii-2
    !         dxi = dx(i)
    !         DO j = 3, jj-2, 2
    !             ddyj = ddy(j)
    !             auz = dxi*ddyj
    !             fut = auz*(wt(k, j, i-1) + wt(k, j+1, i-1) &
    !                 + wt(k, j, i) + wt(k, j+1, i) &
    !                 + 2*wt(k, j, i+1) + 2*wt(k, j+1, i+1))*0.125
    !             qkut = 0.5*fut*(u(k, j, i) + u(k+1, j, i))
    !             wcu(k, j, i) = qkut
    !             qkut = 0.5*fut*(u(k, j+1, i) + u(k+1, j+1, i))
    !             wcu(k, j+1, i) = qkut
    !         END DO

    !         ! VERTEILUNG
    !         k = kk-2
    !         DO i = 2, ii-4, 2
    !             DO j = 3, jj-2
    !                 wcu(k, j, i+1) = 0.5*(wcu(k, j, i) + wcu(k, j, i+2))
    !             END DO
    !         END DO

    !         ! AUF U0 SCHREIBEN
    !         rddzk = rddz(k)
    !         DO i = 2, ii-2
    !             DO j = 3, jj-2
    !                 fkdtu = -1.0*rddy(j)*rdx(i)*wkon
    !                 uo(k, j, i) = uo(k, j, i) &
    !                     + fkdtu*rddzk*(wcu(k, j, i) - wcu(k+1, j, i))
    !             END DO
    !         END DO

    !         ! V-Impulszelle
    !         ! STANDARD-BERECHNUNG DES KONVEKTIVEN FLUSSES
    !         k = kk-2
    !         DO i = 3, ii-2
    !             ddxi = ddx(i)
    !             DO j = 2, jj-2
    !                 dyj = dy(j)
    !                 avz = ddxi*dyj
    !                 fvt = avz*(wt(k, j, i) + wt(k, j+1, i))*0.5
    !                 qkvt = 0.5*fvt*(v(k, j, i) + v(k+1, j, i))
    !                 wcv(k+1, j, i) = qkvt
    !             END DO
    !         END DO

    !         ! NEUE BERECHNUNG DES KONVEKTIVEN FLUSSES
    !         ! FUER JEDE GROBGITTERGESCHWINDIGKEIT
    !         k = kk-2
    !         DO i = 3, ii-2, 2
    !             ddxi = ddx(i)

    !             ! ym-rand
    !             j = 2
    !             dyj = dy(j)
    !             avz = ddxi*dyj
    !             fvt = avz*(2*wt(k, j, i) + 2*wt(k, j, i+1) &
    !                 + wt(k, j+1, i) + wt(k, j+1, i+1) &
    !                 + wt(k, j+2, i) + wt(k, j+2, i+1))*0.125
    !             qkvt = 0.5*fvt*(v(k, j, i) + v(k+1, j, i))
    !             wcv(k, j, i) = qkvt
    !             qkvt = 0.5*fvt*(v(k, j, i+1) + v(k+1, j, i+1))
    !             wcv(k, j, i+1) = qkvt

    !             ! IM GEBIET
    !             DO j = 4, jj-4, 2
    !                 dyj = dy(j)
    !                 avz = ddxi*dyj
    !                 fvt = avz*(wt(k, j-1, i) + wt(k, j-1, i+1) &
    !                     + wt(k, j, i) + wt(k, j, i+1) &
    !                     + wt(k, j+1, i) + wt(k, j+1, i+1) &
    !                     + wt(k, j+2, i) + wt(k, j+2, i+1))*0.125
    !                 qkvt = 0.5*fvt*(v(k, j, i) + v(k+1, j, i))
    !                 wcv(k, j, i) = qkvt
    !                 qkvt = 0.5*fvt*(v(k, j, i+1) + v(k+1, j, i+1))
    !                 wcv(k, j, i+1) = qkvt
    !             END DO

    !             ! yp-rand
    !             j = jj-2
    !             dyj = dy(j)
    !             avz = ddxi*dyj
    !             fvt = avz*(wt(k, j-1, i) + wt(k, j-1, i+1) &
    !                 + wt(k, j, i) + wt(k, j, i+1) &
    !                 + 2*wt(k, j+1, i) + 2*wt(k, j+1, i+1))*0.125
    !             qkvt = 0.5*fvt*(v(k, j, i) + v(k+1, j, i))
    !             wcv(k, j, i) = qkvt
    !             qkvt = 0.5*fvt*(v(k, j, i+1) + v(k+1, j, i+1))
    !             wcv(k, j, i+1) = qkvt
    !         END DO

    !         ! VERTEILUNG
    !         k = kk-2
    !         DO i = 3, ii-2
    !             DO j = 2, jj-4, 2
    !                 wcv(k, j+1, i) = 0.5*(wcv(k, j, i) + wcv(k, j+2, i))
    !             END DO
    !         END DO

    !         ! AUF VO SCHREIBEN
    !         rddzk = rddz(k)
    !         DO i = 3, ii-2
    !             DO j = 2, jj-2
    !                 fkdtv = -1.0*rddx(i)*rdy(j)*wkon
    !                 vo(k, j, i) = vo(k, j, i) &
    !                     + fkdtv*rddzk*(wcv(k, j, i) - wcv(k+1, j, i))
    !             END DO
    !         END DO
    !     END IF

    !     ! PAR-RB Impulserhaltend BOTTOM
    !     IF (nbot == 8) THEN
    !         ! U-Impulszelle
    !         ! STANDARD-BERECHNUNG DES KONVEKTIVEN FLUSSES
    !         k = 3
    !         DO i = 2, ii-2
    !             dxi = dx(i)
    !             DO j = 3, jj-2
    !                 ddyj = ddy(j)
    !                 auz = dxi*ddyj
    !                 fub = auz*(wt(k-1, j, i) + wt(k-1, j, i+1))*0.5
    !                 qkub = 0.5*fub*(u(k, j, i) + u(k-1, j, i))
    !                 wcu(k-1, j, i) = qkub
    !             END DO
    !         END DO

    !         ! NEUE BERECHNUNG DES KONVEKTIVEN FLUSSES
    !         ! FUER JEDE GROBGITTERGESCHWINDIGKEIT
    !         k = 3

    !         ! XM-RAND
    !         i = 2
    !         dxi = dx(i)
    !         DO j = 3, jj-2, 2
    !             ddyj = ddy(j)
    !             auz = dxi*ddyj
    !             fub = auz*(2*wt(k-1, j, i) + 2*wt(k-1, j+1, i) &
    !                 + wt(k-1, j, i+1) + wt(k-1, j+1, i+1) &
    !                 + wt(k-1, j, i+2) + wt(k-1, j+1, i+2))*0.125
    !             qkub = 0.5*fub*(u(k, j, i) + u(k-1, j, i))
    !             wcu(k, j, i) = qkub
    !             qkub = 0.5*fub*(u(k, j+1, i) + u(k-1, j+1, i))
    !             wcu(k, j+1, i) = qkub
    !         END DO

    !         ! IM GEBIET
    !         DO i = 4, ii-4, 2
    !             dxi = dx(i)
    !             DO j = 3, jj-2, 2
    !                 ddyj = ddy(j)
    !                 auz = dxi*ddyj
    !                 fub = auz*(wt(k-1, j, i-1) + wt(k-1, j+1, i-1) &
    !                     + wt(k-1, j, i) + wt(k-1, j+1, i) &
    !                     + wt(k-1, j, i+1) + wt(k-1, j+1, i+1) &
    !                     + wt(k-1, j, i+2) + wt(k-1, j+1, i+2))*0.125
    !                 qkub = 0.5*fub*(u(k, j, i) + u(k-1, j, i))
    !                 wcu(k, j, i) = qkub
    !                 qkub = 0.5*fub*(u(k, j+1, i) + u(k-1, j+1, i))
    !                 wcu(k, j+1, i) = qkub
    !             END DO
    !         END DO

    !         ! XP-RAND
    !         i = ii-2
    !         dxi = dx(i)
    !         DO j = 3, jj-2, 2
    !             ddyj = ddy(j)
    !             auz = dxi*ddyj
    !             fub = auz*(wt(k-1, j, i-1) + wt(k-1, j+1, i-1) &
    !                 + wt(k-1, j, i) + wt(k-1, j+1, i) &
    !                 + 2*wt(k-1, j, i+1) + 2*wt(k-1, j+1, i+1))*0.125
    !             qkub = 0.5*fub*(u(k, j, i) + u(k-1, j, i))
    !             wcu(k, j, i) = qkub
    !             qkub = 0.5*fub*(u(k, j+1, i) + u(k-1, j+1, i))
    !             wcu(k, j+1, i) = qkub
    !         END DO

    !         ! VERTEILUNG
    !         k = 3
    !         DO i = 2, ii-4, 2
    !             DO j = 3, jj-2
    !                 wcu(k, j, i+1) = 0.5*(wcu(k, j, i) + wcu(k, j, i+2))
    !             END DO
    !         END DO

    !         ! AUF U0 SCHREIBEN
    !         rddzk = rddz(k)
    !         DO i = 2, ii-2
    !             DO j = 3, jj-2
    !                 fkdtu = -1.0*rddy(j)*rdx(i)*wkon
    !                 uo(k, j, i) = uo(k, j, i) &
    !                     + fkdtu*rddzk*(-wcu(k, j, i) + wcu(k-1, j, i))
    !             END DO
    !         END DO

    !         ! V-Impulszelle
    !         ! STANDARD-BERECHNUNG DES KONVEKTIVEN FLUSSES
    !         k = 3
    !         DO i = 3, ii-2
    !             ddxi = ddx(i)
    !             DO j = 2, jj-2
    !                 dyj = dy(j)
    !                 avz = ddxi*dyj
    !                 fvb = avz*(wt(k-1, j, i) + wt(k-1, j+1, i))*0.5
    !                 qkvb = 0.5*fvb*(v(k, j, i) + v(k-1, j, i))
    !                 wcv(k-1, j, i) = qkvb
    !             END DO
    !         END DO

    !         ! NEUE BERECHNUNG DES KONVEKTIVEN FLUSSES
    !         ! FUER JEDE GROBGITTERGESCHWINDIGKEIT
    !         k = 3
    !         DO i = 3, ii-2, 2
    !             ddxi = ddx(i)

    !             ! ym-rand
    !             j = 2
    !             dyj = dy(j)
    !             avz = ddxi* dyj
    !             fvb = avz*(2*wt(k-1, j, i) + 2*wt(k-1, j, i+1) &
    !                 + wt(k-1, j+1, i) + wt(k-1, j+1, i+1) &
    !                 + wt(k-1, j+2, i) + wt(k-1, j+2, i+1))*0.125
    !             qkvb = 0.5*fvb*(v(k, j, i) + v(k-1, j, i))
    !             wcv(k, j, i) = qkvb
    !             qkvb = 0.5*fvb*(v(k, j, i+1) + v(k-1, j, i+1))
    !             wcv(k, j, i+1) = qkvb

    !             ! IM GEBIET
    !             DO j = 4, jj-4, 2
    !                 dyj = dy(j)
    !                 avz = ddxi* dyj
    !                 fvb = avz*(wt(k-1, j-1, i) + wt(k-1, j-1, i+1) &
    !                     + wt(k-1, j, i) + wt(k-1, j, i+1) &
    !                     + wt(k-1, j+1, i) + wt(k-1, j+1, i+1) &
    !                     + wt(k-1, j+2, i) + wt(k-1, j+2, i+1))*0.125
    !                 qkvb = 0.5*fvb*(v(k, j, i) + v(k-1, j, i))
    !                 wcv(k, j, i) = qkvb
    !                 qkvb = 0.5*fvb*(v(k, j, i+1) + v(k-1, j, i+1))
    !                 wcv(k, j, i+1) = qkvb

    !             END DO

    !             ! yp-rand
    !             j = jj-2
    !             dyj = dy(j)
    !             avz = ddxi* dyj
    !             fvb = avz*(wt(k-1, j-1, i) + wt(k-1, j-1, i+1) &
    !                 + wt(k-1, j, i) + wt(k-1, j, i+1) &
    !                 + 2*wt(k-1, j+1, i) + 2*wt(k-1, j+1, i+1))*0.125
    !             qkvb = 0.5*fvb*(v(k, j, i) + v(k-1, j, i))
    !             wcv(k, j, i) = qkvb
    !             qkvb = 0.5*fvb*(v(k, j, i+1) + v(k-1, j, i+1))
    !             wcv(k, j, i+1) = qkvb
    !         END DO

    !         ! VERTEILUNG
    !         k = 3
    !         DO i = 3, ii-2
    !             DO j = 2, jj-4, 2
    !                 wcv(k, j+1, i) = 0.5*(wcv(k, j, i) + wcv(k, j+2, i))
    !             END DO
    !         END DO

    !         ! AUF VO SCHREIBEN
    !         rddzk = rddz(k)
    !         DO i = 3, ii-2
    !             DO j = 2, jj-2
    !                 fkdtv = -1.0*rddx(i)*rdy(j)*wkon
    !                 vo(k, j, i) = vo(k, j, i) &
    !                     + fkdtv*rddzk*(-wcv(k, j, i) + wcv(k-1, j, i))
    !             END DO
    !         END DO
    !     END IF

    !     ! PAR-RB Impulserhaltend LEFT
    !     IF (nlft == 8) THEN
    !         ! U-Impulszelle
    !         ! STANDARD-BERECHNUNG DES KONVEKTIVEN FLUSSES
    !         j = jj-2
    !         DO i = 2, ii-2
    !             dxi = dx(i)
    !             DO k = 3, kk-2
    !                 ddzk = ddz(k)
    !                 auy = dxi*ddzk
    !                 fun = auy*(vt(k, j, i) + vt(k, j, i+1))*0.5
    !                 qkun = 0.5*fun*(u(k, j, i) + u(k, j+1, i))
    !                 wcu(k, j+1, i) = qkun
    !             END DO
    !         END DO

    !         ! NEUE BERECHNUNG DES KONVEKTIVEN FLUSSES
    !         ! FUER JEDE GROBGITTERGESCHWINDIGKEIT
    !         j = jj-2

    !         ! XM-RAND
    !         i = 2
    !         dxi = dx(i)
    !         DO k = 3, kk-2, 2
    !             ddzk = ddz(k)
    !             auy = dxi*ddzk
    !             fun = auy*(2*vt(k, j, i) + 2*vt(k+1, j, i) &
    !                 + vt(k, j, i+1) + vt(k+1, j, i+1) &
    !                 + vt(k, j, i+2) + vt(k+1, j, i+2))*0.125
    !             qkun = 0.5*fun*(u(k, j, i) + u(k, j+1, i))
    !             wcu(k, j, i) = qkun
    !             qkun = 0.5*fun*(u(k+1, j, i) + u(k+1, j+1, i))
    !             wcu(k+1, j, i) = qkun
    !         END DO

    !         ! IM GEBIET
    !         DO i = 4, ii-4, 2
    !             dxi = dx(i)
    !             DO k = 3, kk-2, 2
    !                 ddzk = ddz(k)
    !                 auy = dxi*ddzk
    !                 fun = auy*(vt(k, j, i-1) + vt(k+1, j, i-1) &
    !                     + vt(k, j, i) + vt(k+1, j, i) &
    !                     + vt(k, j, i+1) + vt(k+1, j, i+1) &
    !                     + vt(k, j, i+2) + vt(k+1, j, i+2))*0.125
    !                 qkun = 0.5*fun*(u(k, j, i) + u(k, j+1, i))
    !                 wcu(k, j, i) = qkun
    !                 qkun = 0.5*fun*(u(k+1, j, i) + u(k+1, j+1, i))
    !                 wcu(k+1, j, i) = qkun
    !             END DO
    !         END DO

    !         ! XP-RAND
    !         i = ii-2
    !         dxi = dx(i)
    !         DO k = 3, kk-2, 2
    !             ddzk = ddz(k)
    !             auy = dxi*ddzk
    !             fun = auy*(vt(k, j, i-1) + vt(k+1, j, i-1) &
    !                 + vt(k, j, i) + vt(k+1, j, i) &
    !                 + 2*vt(k, j, i+1) + 2*vt(k+1, j, i+1))*0.125
    !             qkun = 0.5*fun*(u(k, j, i) + u(k, j+1, i))
    !             wcu(k, j, i) = qkun
    !             qkun = 0.5*fun*(u(k+1, j, i) + u(k+1, j+1, i))
    !             wcu(k+1, j, i) = qkun
    !         END DO

    !         ! VERTEILUNG
    !         j = jj-2
    !         DO i = 2, ii-4, 2
    !             DO k = 3, kk-2
    !                 wcu(k, j, i+1) = 0.5*(wcu(k, j, i) + wcu(k, j, i+2))
    !             END DO
    !         END DO

    !         ! AUF U0 SCHREIBEN
    !         j = jj-2
    !         DO i = 2, ii-2
    !             fkdtu = -1.0*rddy(j)*rdx(i)*wkon
    !             DO k = 3, kk-2
    !                 rddzk = rddz(k)
    !                 uo(k, j, i) = uo(k, j, i) &
    !                     + fkdtu*rddzk*(wcu(k, j, i) - wcu(k, j+1, i))
    !             END DO
    !         END DO

    !         ! W-IMPULSZELLE
    !         ! STANDARD-BERECHNUNG DES KONVEKTIVEN FLUSSES
    !         j = jj-2
    !         DO i = 3, ii-2
    !             ddxi = ddx(i)
    !             DO k = 2, kk-2
    !                 dzk = dz(k)
    !                 awy = ddxi*dzk
    !                 fwn = awy*(vt(k, j, i) + vt(k+1, j, i))*0.5
    !                 qkwn = 0.5*fwn*(w(k, j, i) + w(k, j+1, i))
    !                 wcw(k, j+1, i) = qkwn
    !             END DO
    !         END DO

    !         ! NEUE BERECHNUNG DES KONVEKTIVEN FLUSSES
    !         ! FUER JEDE GROBGITTERGESCHWINDIGKEIT
    !         j = jj-2
    !         DO i = 3, ii-2, 2
    !             ddxi = ddx(i)

    !             ! ZM-RAND
    !             k = 2
    !             dzk = dz(k)
    !             awy = ddxi*dzk
    !             fwn = awy*(2*vt(k, j, i) + 2*vt(k, j, i+1) &
    !                 + vt(k+1, j, i) + vt(k+1, j, i+1) &
    !                 + vt(k+2, j, i) + vt(k+2, j, i+1))*0.125
    !             qkwn = 0.5*fwn*(w(k, j, i) + w(k, j+1, i))
    !             wcw(k, j, i) = qkwn
    !             qkwn = 0.5*fwn*(w(k, j, i+1) + w(k, j+1, i+1))
    !             wcw(k, j, i+1) = qkwn

    !             ! IM-GEBIET
    !             DO k = 4, kk-4
    !                 dzk = dz(k)
    !                 awy = ddxi*dzk
    !                 fwn = awy*(vt(k-1, j, i) + vt(k-1, j, i+1) &
    !                     + vt(k, j, i) + vt(k, j, i+1) &
    !                     + vt(k+1, j, i) + vt(k+1, j, i+1) &
    !                     + vt(k+2, j, i) + vt(k+2, j, i+1))*0.125
    !                 qkwn = 0.5*fwn*(w(k, j, i) + w(k, j+1, i))
    !                 wcw(k, j, i) = qkwn
    !                 qkwn = 0.5*fwn*(w(k, j, i+1) + w(k, j+1, i+1))
    !                 wcw(k, j, i+1) = qkwn
    !             END DO

    !             ! ZP-RAND
    !             k = kk-2
    !             dzk = dz(k)
    !             awy = ddxi*dzk
    !             fwn = awy*(vt(k-1, j, i) + vt(k-1, j, i+1) &
    !                 + vt(k, j, i) + vt(k, j, i+1) &
    !                 + 2*vt(k+1, j, i) + 2*vt(k+1, j, i+1))*0.125
    !             qkwn = 0.5*fwn*(w(k, j, i) + w(k, j+1, i))
    !             wcw(k, j, i) = qkwn
    !             qkwn = 0.5*fwn*(w(k, j, i+1) + w(k, j+1, i+1))
    !             wcw(k, j, i+1) = qkwn
    !         END DO

    !         ! VERTEILUNG
    !         j = jj-2
    !         DO i = 3, ii-2
    !             DO k = 2, kk-4, 2
    !                 wcw(k+1, j, i) = 0.5*(wcw(k, j, i) + wcw(k+2, j, i))
    !             END DO
    !         END DO

    !         ! AUF W0 SCHREIBEN
    !         j = jj-2
    !         DO i = 3, ii-2
    !             fkdtw = -1.0*rddx(i)*rddy(j)*wkon
    !             DO k = 2, kk-2
    !                 rdzk = rdz(k)
    !                 wo(k, j, i) = wo(k, j, i) &
    !                     + fkdtw*rdzk*(-wcw(k, j+1, i) + wcw(k, j, i))
    !             END DO
    !         END DO
    !     END IF

    !     ! PAR-RB Impulserhaltend RIGHT
    !     IF (nrgt == 8) THEN
    !         ! U-Impulszelle
    !         ! STANDARD-BERECHNUNG DES KONVEKTIVEN FLUSSES
    !         j = 3
    !         DO i = 2, ii-2
    !             dxi = dx(i)
    !             DO k = 3, kk-2
    !                 ddzk = ddz(k)
    !                 auy = dxi*ddzk
    !                 fus = auy*(vt(k, j-1, i) + vt(k, j-1, i+1))*0.5
    !                 qkus = 0.5*fus*(u(k, j, i) + u(k, j-1, i))
    !                 wcu(k, j-1, i) = qkus
    !             END DO
    !         END DO

    !         ! NEUE BERECHNUNG DES KONVEKTIVEN FLUSSES
    !         ! FUER JEDE GROBGITTERGESCHWINDIGKEIT
    !         j = 3

    !         ! XM-RAND
    !         i = 2
    !         dxi = dx(i)
    !         DO k = 3, kk-2, 2
    !             ddzk = ddz(k)
    !             auy = dxi*ddzk
    !             fus = auy*(2*vt(k, j-1, i) + 2*vt(k+1, j-1, i) &
    !                 + vt(k, j-1, i+1) + vt(k+1, j-1, i+1) &
    !                 + vt(k, j-1, i+2) + vt(k+1, j-1, i+2))*0.125
    !             qkus = 0.5*fus*(u(k, j, i) + u(k, j-1, i))
    !             wcu(k, j, i) = qkus
    !             qkus = 0.5*fus*(u(k+1, j, i) + u(k+1, j-1, i))
    !             wcu(k+1, j, i) = qkus
    !         END DO

    !         ! IM GEBIET
    !         DO i = 4, ii-4, 2
    !             dxi = dx(i)
    !             DO k = 3, kk-2, 2
    !                 ddzk = ddz(k)
    !                 auy = dxi*ddzk
    !                 fus = auy*(vt(k, j-1, i-1) + vt(k+1, j-1, i-1) &
    !                     + vt(k, j-1, i) + vt(k+1, j-1, i) &
    !                     + vt(k, j-1, i+1) + vt(k+1, j-1, i+1) &
    !                     + vt(k, j-1, i+2) + vt(k+1, j-1, i+2))*0.125
    !                 qkus = 0.5*fus*(u(k, j, i) + u(k, j-1, i))
    !                 wcu(k, j, i) = qkus
    !                 qkus = 0.5*fus*(u(k+1, j, i) + u(k+1, j-1, i))
    !                 wcu(k+1, j, i) = qkus
    !             END DO
    !         END DO

    !         ! XP-RAND
    !         i = ii-2
    !         dxi = dx(i)
    !         DO k = 3, kk-2, 2
    !             ddzk = ddz(k)
    !             auy = dxi*ddzk
    !             fus = auy*(vt(k, j-1, i-1) + vt(k+1, j-1, i-1) &
    !                 + vt(k, j-1, i) + vt(k+1, j-1, i) &
    !                 +2 *vt(k, j-1, i+1) + 2*vt(k+1, j-1, i+1))*0.125
    !             qkus = 0.5*fus*(u(k, j, i) + u(k, j-1, i))
    !             wcu(k, j, i) = qkus
    !             qkus = 0.5*fus*(u(k+1, j, i) + u(k+1, j-1, i))
    !             wcu(k+1, j, i) = qkus
    !         END DO

    !         ! VERTEILUNG
    !         j = 3
    !         DO i = 2, ii-4, 2
    !             DO k = 3, kk-2
    !                 wcu(k, j, i+1) = 0.5*(wcu(k, j, i) + wcu(k, j, i+2))
    !             END DO
    !         END DO

    !         ! AUF U0 SCHREIBEN
    !         j = 3
    !         DO i = 2, ii-2
    !             fkdtu = -1.0*rddy(j)*rdx(i)*wkon
    !             DO k = 3, kk-2
    !                 rddzk = rddz(k)
    !                 uo(k, j, i) = uo(k, j, i) &
    !                     + fkdtu*rddzk*(-wcu(k, j, i) + wcu(k, j-1, i))
    !             END DO
    !         END DO

    !         ! W-IMPULSZELLE
    !         ! STANDARD-BERECHNUNG DES KONVEKTIVEN FLUSSES
    !         j = 3
    !         DO i = 3, ii-2
    !             ddxi = ddx(i)
    !             DO k = 2, kk-2
    !                 dzk = dz(k)
    !                 awy = ddxi*dzk
    !                 fws = awy*(vt(k, j-1, i) + vt(k+1, j-1, i))*0.5
    !                 qkws = 0.5*fws*(w(k, j, i) + w(k, j-1, i))
    !                 wcw(k, j-1, i) = qkws
    !             END DO
    !         END DO

    !         ! NEUE BERECHNUNG DES KONVEKTIVEN FLUSSES
    !         ! FUER JEDE GROBGITTERGESCHWINDIGKEIT
    !         j = 3
    !         DO i = 3, ii-2, 2
    !             ddxi = ddx(i)

    !             ! zm-rand
    !             k = 2
    !             dzk = dz(k)
    !             awy = ddxi*dzk
    !             fws = awy*(2*vt(k, j-1, i) + 2*vt(k, j-1, i+1) &
    !                 + vt(k+1, j-1, i) + vt(k+1, j-1, i+1) &
    !                 + vt(k+2, j-1, i) + vt(k+2, j-1, i+1))*0.125
    !             qkws = 0.5*fws*(w(k, j, i) + w(k, j-1, i))
    !             wcw(k, j, i) = qkws
    !             qkws = 0.5*fws*(w(k, j, i+1) + w(k, j-1, i+1))
    !             wcw(k, j, i+1) = qkws

    !             ! IM-GEBIET
    !             DO k = 4, kk-4
    !                 dzk = dz(k)
    !                 awy = ddxi*dzk
    !                 fws = awy*(vt(k-1, j-1, i) + vt(k-1, j-1, i+1) &
    !                     + vt(k, j-1, i) + vt(k, j-1, i+1) &
    !                     + vt(k+1, j-1, i) + vt(k+1, j-1, i+1) &
    !                     + vt(k+2, j-1, i) + vt(k+2, j-1, i+1))*0.125
    !                 qkws = 0.5*fws*(w(k, j, i) + w(k, j-1, i))
    !                 wcw(k, j, i) = qkws
    !                 qkws = 0.5*fws*(w(k, j, i+1) + w(k, j-1, i+1))
    !                 wcw(k, j, i+1) = qkws
    !             END DO

    !             ! ZP-RAND
    !             k = kk-2
    !             dzk = dz(k)
    !             awy = ddxi*dzk
    !             fws = awy*(vt(k-1, j-1, i) + vt(k-1, j-1, i+1) &
    !                 + vt(k, j-1, i) + vt(k, j-1, i+1) &
    !                 + 2*vt(k+1, j-1, i) + 2*vt(k+1, j-1, i+1))*0.125
    !             qkws = 0.5*fws*(w(k, j, i) + w(k, j-1, i))
    !             wcw(k, j, i) = qkws
    !             qkws = 0.5*fws*(w(k, j, i+1) + w(k, j-1, i+1))
    !             wcw(k, j, i+1) = qkws
    !         END DO

    !         ! VERTEILUNG
    !         j = 3
    !         DO i = 3, ii-2
    !             DO k = 2, kk-4, 2
    !                 wcw(k+1, j, i) = 0.5*(wcw(k, j, i) + wcw(k+2, j, i))
    !             END DO
    !         END DO

    !         ! AUF W0 SCHREIBEN
    !         j = 3
    !         DO i = 3, ii-2
    !             fkdtw = -1.0*rddx(i)*rddy(j)*wkon
    !             DO k = 2, kk-2
    !                 rdzk = rdz(k)
    !                 wo(k, j, i) = wo(k, j, i) &
    !                     + fkdtw*rdzk*(wcw(k, j-1, i) - wcw(k, j, i))
    !             END DO
    !         END DO
    !     END IF

    !     DEALLOCATE(wcu)
    !     DEALLOCATE(wcv)
    !     DEALLOCATE(wcw)
    ! END SUBROUTINE tstle4_par


    SUBROUTINE swcle3d(kk, jj, ii, uo, vo, wo, u, v, w, ddx, ddy, ddz, &
            nfro, nbac, nrgt, nlft, nbot, ntop)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: uo(kk, jj, ii), vo(kk, jj, ii), &
            wo(kk, jj, ii)
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        INTEGER(intk), INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop

        ! Local variables
        INTEGER(intk) :: k, j, i

        IF (nfro == 5) THEN
            i = 3
            DO j = 2, jj-1
                DO k = 2, kk-1
                    vo(k, j, i) = vo(k, j, i) - swcle3d_one(ddx(i), v(k, j, i))
                    wo(k, j, i) = wo(k, j, i) - swcle3d_one(ddx(i), w(k, j, i))
                END DO
            END DO
        END IF

        IF (nbac == 5) THEN
            i = ii-2
            DO j = 2, jj-1
                DO k = 2, kk-1
                    vo(k, j, i) = vo(k, j, i) - swcle3d_one(ddx(i), v(k, j, i))
                    wo(k, j, i) = wo(k, j, i) - swcle3d_one(ddx(i), w(k, j, i))
                END DO
            END DO
        END IF

        IF (nrgt == 5) THEN
            j = 3
            DO i = 2, ii-1
                DO k = 2, kk-1
                    uo(k, j, i) = uo(k, j, i) - swcle3d_one(ddy(j), u(k, j, i))
                    wo(k, j, i) = wo(k, j, i) - swcle3d_one(ddy(j), w(k, j, i))
                END DO
            END DO
        END IF

        IF (nlft == 5) THEN
            j = jj-2
            DO i = 2, ii-1
                DO k = 2, kk-1
                    uo(k, j, i) = uo(k, j, i) - swcle3d_one(ddy(j), u(k, j, i))
                    wo(k, j, i) = wo(k, j, i) - swcle3d_one(ddy(j), w(k, j, i))
                END DO
            END DO
        END IF

        IF (nbot == 5) THEN
            k = 3
            DO i = 2, ii-1
                DO j = 2, jj-1
                    uo(k, j, i) = uo(k, j, i) - swcle3d_one(ddz(k), u(k, j, i))
                    vo(k, j, i) = vo(k, j, i) - swcle3d_one(ddz(k), v(k, j, i))
                END DO
            END DO
        END IF

        IF (ntop == 5) THEN
            k = kk-2
            DO i = 2, ii-1
                DO j = 2, jj-1
                    uo(k, j, i) = uo(k, j, i) - swcle3d_one(ddz(k), u(k, j, i))
                    vo(k, j, i) = vo(k, j, i) - swcle3d_one(ddz(k), v(k, j, i))
                END DO
            END DO
        END IF
    END SUBROUTINE swcle3d

    PURE SUBROUTINE quick_interpolation_scheme(kk, jj, ii, adveField, &
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
        REAL(realk), INTENT(in) :: adveField(kk, jj, ii)
        REAL(realk), INTENT(out) :: adveE, adveW, adveN, adveS, adveT, adveB
        REAL(realk), INTENT(in) :: advrE, advrW, advrN, advrS, advrT, advrB

        ! Loval variables
        ! None

        !       -----indicator-function----   --------------------------QUICK 3^rd order interpolation-------------------------
        adveE = 0.5 * ( 1 + SIGN(1,advrE) ) * 0.75 * adveField(k,j,i) + 0.375 * adveField(k,j,i+1) - 0.125 * adveField(k,j,i-1) + &
                0.5 * ( 1 - SIGN(1,advrE) ) * 0.75 * adveField(k,j,i+1) + 0.375 * adveField(k,j,i) - 0.125 * adveField(k,j,i+2)
        adveW = 0.5 * ( 1 + SIGN(1,advrW) ) * 0.75 * adveField(k,j,i-1) + 0.375 * adveField(k,j,i) - 0.125 * adveField(k,j,i-2) + &
                0.5 * ( 1 - SIGN(1,advrW) ) * 0.75 * adveField(k,j,i) + 0.375 * adveField(k,j,i-1) - 0.125 * adveField(k,j,i+1)
        adveN = 0.5 * ( 1 + SIGN(1,advrN) ) * 0.75 * adveField(k,j,i) + 0.375 * adveField(k,j+1,i) - 0.125 * adveField(k,j-1,i) + &
                0.5 * ( 1 - SIGN(1,advrN) ) * 0.75 * adveField(k,j+1,i) + 0.375 * adveField(k,j,i) - 0.125 * adveField(k,j+2,i)
        adveS = 0.5 * ( 1 + SIGN(1,advrS) ) * 0.75 * adveField(k,j-1,i) + 0.375 * adveField(k,j,i) - 0.125 * adveField(k,j-2,i) + &
                0.5 * ( 1 - SIGN(1,advrS) ) * 0.75 * adveField(k,j,i) + 0.375 * adveField(k,j-1,i) - 0.125 * adveField(k,j+1,i)
        adveT = 0.5 * ( 1 + SIGN(1,advrT) ) * 0.75 * adveField(k,j,i) + 0.375 * adveField(k+1,j,i) - 0.125 * adveField(k-1,j,i) + &
                0.5 * ( 1 - SIGN(1,advrT) ) * 0.75 * adveField(k+1,j,i) + 0.375 * adveField(k,j,i) - 0.125 * adveField(k+2,j,i)
        adveB = 0.5 * ( 1 + SIGN(1,advrB) ) * 0.75 * adveField(k-1,j,i) + 0.375 * adveField(k,j,i) - 0.125 * adveField(k-2,j,i) + &
                0.5 * ( 1 - SIGN(1,advrB) ) * 0.75 * adveField(k,j,i) + 0.375 * adveField(k-1,j,i) - 0.125 * adveField(k+1,j,i)

    END SUBROUTINE quick_interpolation_scheme

    PURE ELEMENTAL REAL(realk) FUNCTION swcle3d_one(ddz, u) RESULT(uo)
        !$omp declare simd(swcle3d_one)

        ! Function arguments
        REAL(realk), INTENT(in) :: ddz  ! wall normal
        REAL(realk), INTENT(in) :: u    ! adveField

        ! Local variables
        ! none...

        uo = tauwin(u, ddz)/rho/ddz
    END FUNCTION swcle3d_one
END MODULE tstle4_mod
