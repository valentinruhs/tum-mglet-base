!====================================================================
!  Module: multiphase_material_mod
!
!  Responsibilities:
!     - Defines mixture laws for material properties
!     - Computes density as function of volume fraction
!     - Computes viscosity as function of volume fraction
!
!  Author:      Valentin Ruhs
!  Created:     2026-02
!  Last update: 2026-02
!
!====================================================================

MODULE multiphase_material_mod

    USE precision_mod, ONLY: intk, realk
    USE grids_mod, ONLY: nmygrids, mygrids
    USE field_mod, ONLY: field_t
    USE fields_mod, ONLY: get_field
    USE grids_mod, ONLY: get_mgdims, get_mgbasb
    USE multiphasecore_mod, ONLY: gmol1, gmol2, rho1, rho2
    USE err_mod, ONLY: err_abort
    USE multiphase_utils_mod, ONLY: clip_vff

    IMPLICIT NONE
    PRIVATE

    PUBLIC :: init_multiphase_material, finish_multiphase_material, &
        comp_material_property_field, comp_property_face_value_cent, &
        comp_property_face_value_stag

CONTAINS

    SUBROUTINE init_multiphase_material()

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        continue
    END SUBROUTINE init_multiphase_material

    !================================================================

    SUBROUTINE finish_multiphase_material()

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        continue
    END SUBROUTINE finish_multiphase_material

    !================================================================

    SUBROUTINE comp_material_property_field(kk, jj, ii, vff, propFld, prop1, prop2)
    !----------------------------------------------------------------
    !   What it does:
    !   This subroutine computes the weighted material property for
    !   all cells.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: vff(kk, jj, ii)
        REAL(realk), INTENT(out) :: propFld(kk, jj, ii)
        REAL(realk), INTENT(in) :: prop1, prop2

        ! Local variables
        INTEGER(intk) :: k, j, i
        REAL(realk) :: vffClip(kk, jj, ii)

        vffClip = vff
        CALL clip_vff(kk, jj, ii, vffClip)

        propFld = vffClip * ( prop1 - prop2 ) + prop2

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    IF ( propFld(k,j,i) > MAX(prop1, prop2) .OR. propFld(k,j,i) < MIN(prop1, prop2) ) THEN
                        CALL err_abort(155, "material property out of bounds!", __FILE__, __LINE__)
                    ENDIF
                ENDDO
            ENDDO
        ENDDO

    END SUBROUTINE comp_material_property_field

    !================================================================

    SUBROUTINE comp_property_face_value_cent(kk, jj, ii, vff, prop1, prop2, average, pe, pn, pt)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: vff(kk, jj, ii)
        REAL(realk), INTENT(in) :: prop1, prop2
        CHARACTER(len=3), INTENT(in) :: average
        REAL(realk), INTENT(out) :: pe(kk, jj, ii), pn(kk, jj, ii), pt(kk, jj, ii)

        ! Local variables
        REAL(realk) :: propFld(kk, jj, ii)
        INTEGER(intk) :: k, j, i

        CALL comp_material_property_field(kk, jj, ii, vff, propFld, prop1, prop2)

        IF ( average == 'ARI' ) THEN
            DO i = 2, ii-2
                DO j = 2, jj-2
                    DO k = 2, kk-2
                        pe(k,j,i) = ( propFld(k,j,i) + propFld(k,j,i+1) ) / 2.0_realk
                        pn(k,j,i) = ( propFld(k,j,i) + propFld(k,j+1,i) ) / 2.0_realk
                        pt(k,j,i) = ( propFld(k,j,i) + propFld(k+1,j,i) ) / 2.0_realk
                    ENDDO
                ENDDO
            ENDDO
        ELSEIF ( average == 'HAR' ) THEN
            DO i = 2, ii-2
                DO j = 2, jj-2
                    DO k = 2, kk-2
                        pe(k,j,i) = 2.0_realk / ( 1.0_realk / propFld(k,j,i) + 1.0_realk / propFld(k,j,i+1) )
                        pn(k,j,i) = 2.0_realk / ( 1.0_realk / propFld(k,j,i) + 1.0_realk / propFld(k,j+1,i) )
                        pt(k,j,i) = 2.0_realk / ( 1.0_realk / propFld(k,j,i) + 1.0_realk / propFld(k+1,j,i) )
                    ENDDO
                ENDDO
            ENDDO
        ENDIF

    ENDSUBROUTINE comp_property_face_value_cent

    !================================================================

    SUBROUTINE comp_property_face_value_stag(kk, jj, ii, vff, prop1, prop2, pxy, pxz, pyz)
    !----------------------------------------------------------------
    !   What it does:
    !   Compute the "matrix" of a property on the staggered cells
    !   faces. Each cell has 6 faces, hence there are 18 values to 
    !   compute. Since (k,j,i)+ = (k+kq,j+jq,i+iq)- this reduces to
    !   9 values. 
    !
    !   Of these 9, 6 are the same:
    !   xStagN = yStagE, xStagW = yStagS,
    !   xStagT = zStagE, xStagW = zStagB,
    !   yStagT = zStagN, yStagS = zStagB.
    !
    !   Additionally, the diagonal values are trivial, since they
    !   are centered on the pressure grid.
    !
    !   -> 3 different values to compute!
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: vff(kk, jj, ii)
        REAL(realk), INTENT(in) :: prop1, prop2
        REAL(realk), INTENT(out) :: pxy(kk, jj, ii), pxz(kk, jj, ii), pyz(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i
        REAL(realk) :: phi

        DO i = 2, ii-2
            DO j = 2, jj-2
                DO k = 2, kk-2
                    phi = MIN(MAX(0.5_realk*(vff(k,j,i) + vff(k,j,i+1) + vff(k,j+1,i) + vff(k,j+1,i+1)) - 0.5_realk, 0.0_realk), 1.0_realk)
                    pxy(k,j,i) = 1.0_realk / ( phi/prop1 + (1.0_realk - phi)/prop2 )

                    phi = MIN(MAX(0.5_realk*(vff(k,j,i) + vff(k,j,i+1) + vff(k+1,j,i) + vff(k+1,j,i+1)) - 0.5_realk, 0.0_realk), 1.0_realk)
                    pxz(k,j,i) = 1.0_realk / ( phi/prop1 + (1.0_realk - phi)/prop2 )

                    phi = MIN(MAX(0.5_realk*(vff(k,j,i) + vff(k,j+1,i) + vff(k+1,j,i) + vff(k+1,j+1,i)) - 0.5_realk, 0.0_realk), 1.0_realk)
                    pyz(k,j,i) = 1.0_realk / ( phi/prop1 + (1.0_realk - phi)/prop2 )
                ENDDO
            ENDDO
        ENDDO

    ENDSUBROUTINE comp_property_face_value_stag

END MODULE multiphase_material_mod