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

    IMPLICIT NONE
    PRIVATE

    PUBLIC :: init_multiphase_material, finish_multiphase_material

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

    SUBROUTINE get_material_property_field(kk, jj, ii, propertyField, vff, propertyFluid1, propertyFluid2, dx, dy, dz)
    !----------------------------------------------------------------
    !   What it does:
    !   This subroutine computes the weighted material property for
    !   all cells.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(out) :: propertyField(kk, jj, ii)
        REAL(realk), INTENT(in) :: vff(kk, jj, ii)
        REAL(realk), INTENT(in) :: propertyFluid1, propertyFluid2
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)

        ! Local variables
        INTEGER(intk) :: k, j, i

        DO i = 1, ii
            DO j = 1, jj
                DO k = 1, kk
                    propertyField(k,j,i) = propertyFluid1 * vff(k,j,i) + propertyFluid2 * ( 1 - vff(k,j,i) )
                END DO
            END DO
        END DO

    END SUBROUTINE get_material_property_field

END MODULE multiphase_material_mod