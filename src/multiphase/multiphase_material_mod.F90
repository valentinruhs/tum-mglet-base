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
    USE err_mod, ONLY: errr

    IMPLICIT NONE
    PRIVATE

    PUBLIC :: init_multiphase_material, finish_multiphase_material, comp_material_property_field, comp_property_face_value

CONTAINS

    SUBROUTINE init_multiphase_material()

        ! Subroutine arguments
        ! None

        ! Local variables
        TYPE(field_t), POINTER :: d_f
        TYPE(field_t), POINTER :: vff_f
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: i, igrid
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: d, vff

        CALL get_field(d_f, "D")
        CALL get_field(vff_f, "VFF")

        DO i = 1, nmygrids
            igrid = mygrids(i)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL d_f%get_ptr(d, igrid)
            CALL vff_f%get_ptr(vff, igrid)

            CALL comp_material_property_field(kk, jj, ii, vff, d, rho1, rho2, 'ARI')
        END DO

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

    SUBROUTINE comp_material_property_field(kk, jj, ii, vff, propertyField, propertyFluid1, propertyFluid2, average)
    !----------------------------------------------------------------
    !   What it does:
    !   This subroutine computes the weighted material property for
    !   all cells.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: vff(kk, jj, ii)
        REAL(realk), INTENT(out) :: propertyField(kk, jj, ii)
        REAL(realk), INTENT(in) :: propertyFluid1, propertyFluid2
        CHARACTER(len=3), INTENT(in) :: average

        ! Local variables
        ! None

        IF ( average == 'ARI' ) THEN
            propertyField = vff * ( propertyFluid1 - propertyFluid2 ) + propertyFluid2
        ELSEIF ( average == 'HAR' ) THEN
            propertyField = 1.0_realk / ( vff * ( 1.0_realk / propertyFluid1 - 1.0_realk / propertyFluid2 ) + 1.0_realk / propertyFluid2 )
        ELSE 
            CALL errr(__FILE__, __LINE__)
        ENDIF

    END SUBROUTINE comp_material_property_field

    !================================================================

    SUBROUTINE comp_property_face_value(kk, jj, ii, vff, propertyFluid1, propertyFluid2, average, pe, pn, pt)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: vff(kk, jj, ii)
        REAL(realk), INTENT(in) :: propertyFluid1, propertyFluid2
        CHARACTER(len=3), INTENT(in) :: average
        REAL(realk), INTENT(out) :: pe(kk, jj, ii), pn(kk, jj, ii), pt(kk, jj, ii)

        ! Local variables
        REAL(realk) :: propertyField(kk, jj, ii)
        INTEGER(intk) :: k, j, i

        CALL comp_material_property_field(kk, jj, ii, vff, propertyField, propertyFluid1, propertyFluid2, 'ARI')

        IF ( average == 'ARI' ) THEN
            DO i = 2, ii-2
                DO j = 2, jj-2
                    DO k = 2, kk-2
                        pe(k,j,i) = ( propertyField(k,j,i) + propertyField(k,j,i+1) ) / 2.0_realk
                        pn(k,j,i) = ( propertyField(k,j,i) + propertyField(k,j+1,i) ) / 2.0_realk
                        pt(k,j,i) = ( propertyField(k,j,i) + propertyField(k+1,j,i) ) / 2.0_realk
                    ENDDO
                ENDDO
            ENDDO
        ELSEIF ( average == 'HAR' ) THEN
            DO i = 2, ii-2
                DO j = 2, jj-2
                    DO k = 2, kk-2
                        pe(k,j,i) = 2.0_realk / ( 1.0_realk / propertyField(k,j,i) + 1.0_realk / propertyField(k,j,i+1) )
                        pn(k,j,i) = 2.0_realk / ( 1.0_realk / propertyField(k,j,i) + 1.0_realk / propertyField(k,j+1,i) )
                        pt(k,j,i) = 2.0_realk / ( 1.0_realk / propertyField(k,j,i) + 1.0_realk / propertyField(k+1,j,i) )
                    ENDDO
                ENDDO
            ENDDO
        ELSE 
            CALL errr(__FILE__, __LINE__)
        ENDIF

    ENDSUBROUTINE comp_property_face_value

END MODULE multiphase_material_mod