!====================================================================
!  Module: multiphase_mod
!
!  Responsibilities:
!     - 
!
!  Coordinates:
!     - 
!
!  Author:      Valentin Ruhs
!  Created:     2026-02
!  Last update: 2026-02
!
!====================================================================

MODULE multiphase_mod

    USE multiphasecore_mod, ONLY: init_multiphasecore, finish_multiphasecore, has_multiphase, solve_multiphase
    USE multiphase_vof_transport_mod, ONLY: init_multiphase_vof_transport, finish_multiphase_vof_transport
    USE multiphase_plic_mod, ONLY: init_multiphase_plic, finish_multiphase_plic
    USE multiphase_material_mod, ONLY: init_multiphase_material, finish_multiphase_material, comp_property_face_value
    USE multiphase_io_mod, ONLY: init_multiphase_io, finish_multiphase_io
    USE multiphasecore_mod, ONLY: gmol1, gmol2, rho1, rho2
    USE multiphase_utils_mod, ONLY: init_multiphase_utils, finish_multiphase_utils
    USE precision_mod, ONLY: intk, realk
    USE grids_mod, ONLY: nmygrids, mygrids, get_mgdims, get_mgbasb
    USE fields_mod, ONLY: get_field, set_field, get_fieldptr
    USE field_mod, ONLY: field_t
    
    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: init_multiphase, finish_multiphase, comp_matrix_coeff_multiphase, comp_factor_coeff_multiphase

CONTAINS

    SUBROUTINE init_multiphase()
        
        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        CALL init_multiphasecore()
        CALL init_multiphase_utils()
        CALL init_multiphase_material()
        CALL init_multiphase_vof_transport()
        CALL init_multiphase_plic()
        CALL init_multiphase_io()
        IF(.NOT. has_multiphase) RETURN
        IF(.NOT. solve_multiphase) RETURN


    END SUBROUTINE init_multiphase

    !================================================================

    SUBROUTINE finish_multiphase()

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        IF(.NOT. has_multiphase) RETURN

        CALL finish_multiphase_plic()
        CALL finish_multiphase_vof_transport()
        CALL finish_multiphase_material()
        CALL finish_multiphase_io()
        CALL finish_multiphase_utils()
        CALL finish_multiphasecore()

    END SUBROUTINE finish_multiphase

    !================================================================

    SUBROUTINE comp_matrix_coeff_multiphase()
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        ! None

        ! Local variables
        TYPE(field_t), POINTER :: vff_f
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: igr, igrid

        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ae(:,:,:), aw(:,:,:), &
                                            an(:,:,:), as(:,:,:), &
                                            at(:,:,:), ab(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: ap(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: vff(:, :, :)
        REAL(realk), ALLOCATABLE :: rhoe(:, :, :), rhon(:, :, :), rhot(:, :, :)

        CALL get_field(vff_f, "VFF")

        DO igr = 1, nmygrids
            igrid = mygrids(igr)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_fieldptr(dx, "DX", igrid)
            CALL get_fieldptr(dy, "DY", igrid)
            CALL get_fieldptr(dz, "DZ", igrid)

            CALL get_fieldptr(aw, "GSAW", igrid)
            CALL get_fieldptr(ae, "GSAE", igrid)
            CALL get_fieldptr(as, "GSAS", igrid)
            CALL get_fieldptr(an, "GSAN", igrid)
            CALL get_fieldptr(ab, "GSAB", igrid)
            CALL get_fieldptr(at, "GSAT", igrid)

            CALL get_fieldptr(ap, "GSAP", igrid)

            CALL vff_f%get_ptr(vff, igrid)

            IF ( .NOT. ALLOCATED(rhoe)) ALLOCATE(rhoe(kk, jj, ii))
            IF ( .NOT. ALLOCATED(rhon)) ALLOCATE(rhon(kk, jj, ii))
            IF ( .NOT. ALLOCATED(rhot)) ALLOCATE(rhot(kk, jj, ii))

            CALL comp_property_face_value(kk, jj, ii, vff, rho1, rho2, 'ARI', rhoe, rhon, rhot)

            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        ae(k,j,i) = 2.0/((dx(i-1)+dx(i))*dx(i)*rhoe(k,j,i))
                        aw(k,j,i) = 2.0/((dx(i-1)+dx(i))*dx(i-1)*rhoe(k,j,i-1))
                        an(k,j,i) = 2.0/((dy(j-1)+dy(j))*dy(j)*rhon(k,j,i))
                        as(k,j,i) = 2.0/((dy(j-1)+dy(j))*dy(j-1)*rhon(k,j-1,i))
                        at(k,j,i) = 2.0/((dz(k-1)+dz(k))*dz(k)*rhot(k,j,i))
                        ab(k,j,i) = 2.0/((dz(k-1)+dz(k))*dz(k-1)*rhot(k-1,j,i))
                    ENDDO
                ENDDO
            ENDDO

            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        ap(k, j, i) = -( ae(k,j,i) + aw(k,j,i) + an(k,j,i) &
                                       + as(k,j,i) + at(k,j,i) + ab(k,j,i) )
                    END DO
                END DO
            END DO

            IF ( ALLOCATED(rhot)) DEALLOCATE(rhot)
            IF ( ALLOCATED(rhon)) DEALLOCATE(rhon)
            IF ( ALLOCATED(rhoe)) DEALLOCATE(rhoe)
            
        ENDDO

    END SUBROUTINE comp_matrix_coeff_multiphase

    !================================================================

    SUBROUTINE comp_factor_coeff_multiphase()
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        ! None

        ! Local variables
        INTEGER(intk) :: igr, igrid
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), PARAMETER :: alfa = 0.92
        REAL(realk) :: p1, p2, p3
        REAL(realk), POINTER, CONTIGUOUS :: bp(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: lb(:, :, :), lw(:, :, :), &
            ls(:, :, :), ue(:, :, :), un(:, :, :), ut(:, :, :), lpr(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: ae(:,:,:), aw(:,:,:), &
                                            an(:,:,:), as(:,:,:), &
                                            at(:,:,:), ab(:,:,:), ap(:, :, :)

        DO igr = 1, nmygrids
            igrid = mygrids(igr)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_fieldptr(lw, "SIPLW", igrid)
            CALL get_fieldptr(ls, "SIPLS", igrid)
            CALL get_fieldptr(lb, "SIPLB", igrid)
            CALL get_fieldptr(ue, "SIPUE", igrid)
            CALL get_fieldptr(un, "SIPUN", igrid)
            CALL get_fieldptr(ut, "SIPUT", igrid)
            CALL get_fieldptr(lpr, "SIPLPR", igrid)

            CALL get_fieldptr(aw, "GSAW", igrid)
            CALL get_fieldptr(ae, "GSAE", igrid)
            CALL get_fieldptr(as, "GSAS", igrid)
            CALL get_fieldptr(an, "GSAN", igrid)
            CALL get_fieldptr(ab, "GSAB", igrid)
            CALL get_fieldptr(at, "GSAT", igrid)

            CALL get_fieldptr(ap, "GSAP", igrid)

            CALL get_fieldptr(bp, "BP", igrid)

            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        lw(k, j, i) = aw(k,j,i)*bu(k, j, i-1) &
                            /(1.0 + alfa*(un(k, j, i-1) + ut(k, j, i-1)))

                        ls(k, j, i) = as(k,j,i)*bv(k, j-1, i) &
                            /(1.0 + alfa*(ue(k, j-1, i) + ut(k, j-1, i)))

                        lb(k, j, i) = ab(k,j,i)*bw(k-1, j, i) &
                            /(1.0 + alfa*(un(k-1, j, i) + ue(k-1, j, i)))

                        p1 = alfa*(lb(k, j, i)*ue(k-1, j, i) &
                            + ls(k, j, i)*ue(k, j-1, i))
                        p2 = alfa*(lb(k, j, i)*un(k-1, j, i) &
                            + lw(k, j, i)*un(k, j, i-1))
                        p3 = alfa*(lw(k, j, i)*ut(k, j, i-1) &
                            + ls(k, j, i)*ut(k, j-1, i))

                        lpr(k, j, i) = 1.0/(ap(k, j, i) + p1 + p2 + p3 &
                            - lb(k, j, i)*ut(k-1, j, i) &
                            - lw(k, j, i)*ue(k, j, i-1) &
                            - ls(k, j, i)*un(k, j-1, i) &
                            + 1.0e-20)

                        ue(k, j, i) = (ae(k,j,i)*bu(k, j, i) - p1)*lpr(k, j, i)
                        un(k, j, i) = (an(k,j,i)*bv(k, j, i) - p2)*lpr(k, j, i)
                        ut(k, j, i) = (at(k,j,i)*bw(k, j, i) - p3)*lpr(k, j, i)
                    END DO
                END DO
            END DO

            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        lw(k, j, i) = lw(k, j, i)*lpr(k, j, i)
                    END DO
                    DO k = 3, kk-2
                        ls(k, j, i) = ls(k, j, i)*lpr(k, j, i)
                    END DO
                    DO k = 3, kk-2
                        lb(k, j, i) = lb(k, j, i)*lpr(k, j, i)
                    END DO
                END DO
            END DO
        END DO

    CONTAINS
        ! These routines are not protected against out-of-bounds and not valid
        ! for the edges ii, jj, kk - but in the scope above that is OK
        PURE REAL(realk) FUNCTION bu(k, j, i)
            INTEGER(intk), INTENT(in) :: k, j, i
            bu = bp(k, j, i)*bp(k, j, i+1)
        END FUNCTION bu

        PURE REAL(realk) FUNCTION bv(k, j, i)
            INTEGER(intk), INTENT(in) :: k, j, i
            bv = bp(k, j, i)*bp(k, j+1, i)
        END FUNCTION bv

        PURE REAL(realk) FUNCTION bw(k, j, i)
            INTEGER(intk), INTENT(in) :: k, j, i
            bw = bp(k, j, i)*bp(k+1, j, i)
        END FUNCTION bw

    END SUBROUTINE comp_factor_coeff_multiphase

END MODULE multiphase_mod