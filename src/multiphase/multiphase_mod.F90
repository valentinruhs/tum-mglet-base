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
!     - multiphase_vof_transport_mod
!     - multiphase_levelset_transport_mod
!     - multiphase_material_mod
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
    USE multiphase_io_mod, ONLY: init_multiphase_io, finish_multiphase_io, read_vff
    USE fields_mod, ONLY: get_field
    USE field_mod, ONLY: field_t
    
    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: init_multiphase, finish_multiphase, init_vff

CONTAINS

    SUBROUTINE init_multiphase()
        
        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        CALL init_multiphasecore()
        CALL init_multiphase_io()
        CALL init_multiphase_material()
        CALL init_multiphase_vof_transport()
        CALL init_multiphase_plic()
        IF(.NOT. has_multiphase) RETURN

        CALL init_vff()
        

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
        CALL finish_multiphasecore()

    END SUBROUTINE finish_multiphase

    !================================================================

    SUBROUTINE init_vff()
    !----------------------------------------------------------------
    !   What it does:
    !   The subroutine manages the allocation of the volume fraction
    !   field. 
    !----------------------------------------------------------------

        TYPE(field_t), POINTER :: vff
        
        CALL get_field(vff, "VFF")

        CALL read_vff(vff)

    END SUBROUTINE init_vff

END MODULE multiphase_mod