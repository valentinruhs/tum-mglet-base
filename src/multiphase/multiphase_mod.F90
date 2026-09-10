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
    USE multiphase_material_mod, ONLY: init_multiphase_material, finish_multiphase_material, comp_property_face_value_cent
    USE multiphase_io_mod, ONLY: init_multiphase_io, finish_multiphase_io
    USE multiphasecore_mod, ONLY: gmol1, gmol2, rho1, rho2
    USE multiphase_utils_mod, ONLY: init_multiphase_utils, finish_multiphase_utils
    USE precision_mod, ONLY: intk, realk
    USE grids_mod, ONLY: nmygrids, mygrids, get_mgdims, get_mgbasb, idprocofgrd, iposition, jposition, kposition, iparent, ngrid
    USE fields_mod, ONLY: get_field, set_field, get_fieldptr
    USE field_mod, ONLY: field_t
    USE comms_mod, ONLY: myid
    
    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: init_multiphase, finish_multiphase, comp_matrix_coeff_multiphase

CONTAINS

    SUBROUTINE init_multiphase()
        
        ! Subroutine arguments
        ! None

        ! Local variables
        INTEGER(intk) :: n, igrid, igridf, ipar
        INTEGER(intk) :: kk, jj, ii, kc0, jc0, ic0
        REAL(realk), POINTER, CONTIGUOUS :: grdMask(:,:,:)

        IF(.NOT. has_multiphase) RETURN

        CALL init_multiphasecore()
        CALL init_multiphase_utils()
        CALL init_multiphase_material()
        CALL init_multiphase_vof_transport()
        CALL init_multiphase_plic()

        DO n = 1, nmygrids
            igrid = mygrids(n)
            CALL get_fieldptr(grdMask, "GRDMASK", igrid)
            grdMask = 1.0_realk
        END DO

        DO igridf = 1, ngrid
            ipar = iparent(igridf)
            IF (ipar == 0) CYCLE
            IF (idprocofgrd(ipar) /= myid) CYCLE

            CALL get_fieldptr(grdMask, "GRDMASK", ipar)
            CALL get_mgdims(kk, jj, ii, igridf)

            ic0 = iposition(igridf)
            jc0 = jposition(igridf)
            kc0 = kposition(igridf)

            grdMask(kc0:kc0+(kk-4)/2-1, &
                    jc0:jc0+(jj-4)/2-1, &
                    ic0:ic0+(ii-4)/2-1) = 0.0_realk
        END DO

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
        INTEGER(intk) :: n, igrid

        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ae(:,:,:), aw(:,:,:), &
                                            an(:,:,:), as(:,:,:), &
                                            at(:,:,:), ab(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: ap(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: vff(:, :, :)
        REAL(realk), ALLOCATABLE :: rhoe(:, :, :), rhon(:, :, :), rhot(:, :, :)

        CALL get_field(vff_f, "VFF")

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

            CALL vff_f%get_ptr(vff, igrid)

            IF ( .NOT. ALLOCATED(rhoe)) ALLOCATE(rhoe(kk, jj, ii))
            IF ( .NOT. ALLOCATED(rhon)) ALLOCATE(rhon(kk, jj, ii))
            IF ( .NOT. ALLOCATED(rhot)) ALLOCATE(rhot(kk, jj, ii))

            CALL comp_property_face_value_cent(kk, jj, ii, vff, rho1, rho2, 'ARI', rhoe, rhon, rhot, ddx, ddy, ddz)

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

END MODULE multiphase_mod