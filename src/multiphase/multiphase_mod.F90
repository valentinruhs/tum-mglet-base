!====================================================================
!  Module: multiphase_mod
!
!  Responsibilities:
!     - Acts as the high-level driver for the multiphase model
!     - Couples multiphase physics to the flow solver
!     - Advances the volume fraction field in time
!     - Updates material properties (rho, mu) based on C
!
!  Coordinates:
!     - multiphase_transport_mod
!     - multiphase_material_mod
!     - multiphase_interface_mod
!
!  Author:      Valentin Ruhs
!  Created:     2026-02
!  Last update: 2026-02
!
!====================================================================

MODULE multiphase_mod

    USE multiphasecore_mod, ONLY: init_multiphasecore, finish_multiphasecore, has_multiphase, solve_multiphase
    USE fields_mod, ONLY: get_field
    
    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: init_multiphase, finish_multiphase

CONTAINS

    SUBROUTINE init_multiphase()
        
        ! Subroutine arguments
        ! None

        ! Local variables
        

        CALL init_multiphasecore()
        IF(.NOT. has_multiphase) RETURN

        CALL init_c()
        

        IF(.NOT. solve_multiphase) RETURN


    END SUBROUTINE init_multiphase

    !================================================================

    SUBROUTINE finish_multiphase()

        IF(.NOT. has_multiphase) RETURN

        CALL finish_multiphasecore()

    END SUBROUTINE finish_multiphase

    !================================================================

    SUBROUTINE init_c()
        
        CALL get_field(c, "C")

    END SUBROUTINE init_c

END MODULE multiphase_mod