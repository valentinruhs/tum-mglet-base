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
    USE multiphase_material_mod, ONLY: init_multiphase_material, finish_multiphase_material
    USE multiphase_io_mod, ONLY: init_multiphase_io, finish_multiphase_io
    USE multiphasecore_mod, ONLY: gmol1, gmol2, rho1, rho2
    USE multiphase_utils_mod, ONLY: init_multiphase_utils, finish_multiphase_utils
    USE precision_mod, ONLY: intk, realk
    USE grids_mod, ONLY: nmygrids, mygrids, get_mgdims, get_mgbasb
    USE fields_mod, ONLY: get_field
    USE field_mod, ONLY: field_t
    
    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: init_multiphase, finish_multiphase

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

END MODULE multiphase_mod