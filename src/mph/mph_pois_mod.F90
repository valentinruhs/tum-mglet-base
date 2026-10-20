!====================================================================
!  Module: mph_pois_mod
!
!   Responsibilities:
!   - 
!
!   Author:      Valentin Ruhs
!   e-Mail:      valentin.ruhs@gmx.de
!   Created:     2026-09
!   Last update: 2026-10
!
!====================================================================

MODULE mph_pois_mod

    USE precision_mod, ONLY: intk, realk
    USE grids_mod, ONLY: nmygrids, mygrids, get_mgdims, get_mgbasb, &
        minlevel, maxlevel
    USE field_mod, ONLY: field_t
    USE fields_mod, ONLY: get_fieldptr, get_field, set_field
    USE ftoc_mod, ONLY: ftoc

    USE mph_xfer_mod, ONLY: rstr

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: init_mph_pois, finish_mph_pois, comp_mat_coeff_mph

CONTAINS

    SUBROUTINE init_mph_pois()

        ! Subroutine arguments
        ! None

        ! Local variables
        CHARACTER(len=*), PARAMETER :: descBeta = "face cond. 1/rho"

        CALL set_field("BETABA", description=descBeta)
        CALL set_field("BETALE", description=descBeta)
        CALL set_field("BETATO", description=descBeta)

    END SUBROUTINE init_mph_pois

    !================================================================

    SUBROUTINE finish_mph_pois()

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        CONTINUE

    END SUBROUTINE finish_mph_pois

    !================================================================

    SUBROUTINE comp_mat_coeff_mph()
    !----------------------------------------------------------------
    !   What it does:
    !   Computes the coefficients of the pressure Poisson equation
    !   from the cell-wise face conductivities (see comp_beta). On the
    !   finest level they equal 1/rho_face, i.e. the coefficients are
    !   identical to the former formulation with DBA/DLE/DTO.
    !----------------------------------------------------------------

        ! Subroutine arguments
        ! None

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: n, igrid
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ae(:,:,:), aw(:,:,:), &
                                            an(:,:,:), as(:,:,:), &
                                            at(:,:,:), ab(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: ap(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: betaBa(:,:,:), betaLe(:,:,:), betaTo(:,:,:)
        INTEGER(intk) :: nfro, nbac, nrgt, nlft, nbot, ntop

        CALL comp_beta()

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_fieldptr(dx, "DX", igrid)
            CALL get_fieldptr(dy, "DY", igrid)
            CALL get_fieldptr(dz, "DZ", igrid)
            CALL get_fieldptr(ddx, "DDX", igrid)
            CALL get_fieldptr(ddy, "DDY", igrid)
            CALL get_fieldptr(ddz, "DDZ", igrid)

            CALL get_fieldptr(aw, "GSAW", igrid)
            CALL get_fieldptr(ae, "GSAE", igrid)
            CALL get_fieldptr(as, "GSAS", igrid)
            CALL get_fieldptr(an, "GSAN", igrid)
            CALL get_fieldptr(ab, "GSAB", igrid)
            CALL get_fieldptr(at, "GSAT", igrid)
            CALL get_fieldptr(ap, "GSAP", igrid)

            CALL get_fieldptr(betaBa, "BETABA", igrid)
            CALL get_fieldptr(betaLe, "BETALE", igrid)
            CALL get_fieldptr(betaTo, "BETATO", igrid)

            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        aw(k,j,i) = (2.0_realk*betaBa(k,j,i-1))/((dx(i-1)+dx(i))*dx(i-1))
                        ae(k,j,i) = (2.0_realk*betaBa(k,j,i))/((dx(i-1)+dx(i))*dx(i))
                        as(k,j,i) = (2.0_realk*betaLe(k,j-1,i))/((dy(j-1)+dy(j))*dy(j-1))
                        an(k,j,i) = (2.0_realk*betaLe(k,j,i))/((dy(j-1)+dy(j))*dy(j))
                        ab(k,j,i) = (2.0_realk*betaTo(k-1,j,i))/((dz(k-1)+dz(k))*dz(k-1))
                        at(k,j,i) = (2.0_realk*betaTo(k,j,i))/((dz(k-1)+dz(k))*dz(k))
                    ENDDO
                ENDDO
            ENDDO

            CALL get_mgbasb(nfro, nbac, nrgt, nlft, nbot, ntop, igrid)
            IF (ANY(nfro == [2, 5, 6, 19])) aw(:, :, 3) = 0.0_realk
            IF (ANY(nbac == [2, 5, 6, 19])) ae(:, :, ii-2) = 0.0_realk
            IF (ANY(nrgt == [2, 5, 6, 19])) as(:, 3, :) = 0.0_realk
            IF (ANY(nlft == [2, 5, 6, 19])) an(:, jj-2, :) = 0.0_realk
            IF (ANY(nbot == [2, 5, 6, 19])) ab(3, :, :) = 0.0_realk
            IF (ANY(ntop == [2, 5, 6, 19])) at(kk-2, :, :) = 0.0_realk

            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        ap(k, j, i) = -(ae(k,j,i) + aw(k,j,i) + an(k,j,i) + &
                                        as(k,j,i) + at(k,j,i) + ab(k,j,i))
                    END DO
                END DO
            END DO
        ENDDO

    END SUBROUTINE comp_mat_coeff_mph

    !================================================================

    SUBROUTINE comp_beta()
    !----------------------------------------------------------------
    !   What it does:
    !   see PARIS 
    !       fill_coefficients
    !       get_MGmatrix_coef
    !       NewSolverMG
    !----------------------------------------------------------------

        ! Subroutine arguments
        ! None

        ! Local variables
        INTEGER(intk) :: n, igrid
        REAL(realk), POINTER, CONTIGUOUS :: betaBa(:,:,:), betaLe(:,:,:), betaTo(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: dBa(:,:,:), dLe(:,:,:), dTo(:,:,:)

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_fieldptr(betaBa, "BETABA", igrid)
            CALL get_fieldptr(betaLe, "BETALE", igrid)
            CALL get_fieldptr(betaTo, "BETATO", igrid)
            CALL get_fieldptr(dBa, "DBA", igrid)
            CALL get_fieldptr(dLe, "DLE", igrid)
            CALL get_fieldptr(dTo, "DTO", igrid)

            betaBa = inv(dBa)
            betaLe = inv(dLe)
            betaTo = inv(dTo)
        END DO

        CALL rstr("BETABA", "U")
        CALL rstr("BETALE", "V")
        CALL rstr("BETATO", "W")

    CONTAINS

        PURE ELEMENTAL REAL(realk) FUNCTION inv(d) RESULT(r)
            REAL(realk), INTENT(in) :: d
            IF (d > 0.0_realk) THEN
                r = 1.0_realk/d
            ELSE
                r = 0.0_realk
            END IF
        END FUNCTION inv

    END SUBROUTINE comp_beta

END MODULE mph_pois_mod