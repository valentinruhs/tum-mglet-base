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

    PUBLIC :: init_multiphase_material, finish_multiphase_material, comp_material_property_field

CONTAINS

    SUBROUTINE init_multiphase_material()

        ! Subroutine arguments
        ! None

        ! Local variables
        TYPE(field_t), POINTER :: d_f
        TYPE(field_t), POINTER :: g_f
        TYPE(field_t), POINTER :: vff_f
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: i, igrid
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: d, g, vff

        CALL get_field(d_f, "D")
        CALL get_field(g_f, "G")
        CALL get_field(vff_f, "VFF")

        DO i = 1, nmygrids
            igrid = mygrids(i)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL d_f%get_ptr(d, igrid)
            CALL g_f%get_ptr(g, igrid)
            CALL vff_f%get_ptr(vff, igrid)

            CALL comp_material_property_field(kk, jj, ii, vff, rho1, rho2, d)
            CALL comp_material_property_field(kk, jj, ii, vff, gmol1, gmol2, g)
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

    SUBROUTINE comp_material_property_field(kk, jj, ii, vff, propertyFluid1, propertyFluid2, propertyField, harmonic)
    !----------------------------------------------------------------
    !   What it does:
    !   This subroutine computes the weighted material property for
    !   all cells.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: vff(kk, jj, ii)
        REAL(realk), INTENT(in) :: propertyFluid1, propertyFluid2
        REAL(realk), INTENT(inout) :: propertyField(kk, jj, ii)
        LOGICAL, OPTIONAL :: harmonic

        ! Local variables
        ! None

        IF ( .NOT. PRESENT(harmonic) ) THEN
            propertyField = vff * ( propertyFluid1 - propertyFluid2 ) + propertyFluid2
        ELSEIF ( harmonic ) THEN
            propertyField = 1.0_realk / ( vff * ( 1.0_realk / propertyFluid1 - 1.0_realk / propertyFluid2 ) + 1.0_realk / propertyFluid2 )
        ELSE 
            CALL errr(__FILE__, __LINE__)
        ENDIF

    END SUBROUTINE comp_material_property_field

END MODULE multiphase_material_mod