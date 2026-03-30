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

    USE multiphasecore_mod, ONLY: rho1, rho2, gmol1, gmol2
    
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

END MODULE multiphase_material_mod