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
    USE multiphase_plic_mod, ONLY: init_multiphase_plic, finish_multiphase_plic, compute_iStag_vff, compute_jStag_vff, compute_kStag_vff, track_interface, compute_normal_vector, compute_alpha
    USE multiphase_material_mod, ONLY: init_multiphase_material, finish_multiphase_material, get_material_property_field
    USE multiphase_io_mod, ONLY: init_multiphase_io, finish_multiphase_io, read_vff
    USE multiphasecore_mod, ONLY: gmol1, gmol2, rho1, rho2
    USE precision_mod, ONLY: intk, realk
    USE grids_mod, ONLY: nmygrids, mygrids, get_mgdims, get_mgbasb
    USE fields_mod, ONLY: get_field
    USE field_mod, ONLY: field_t
    
    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: init_multiphase, finish_multiphase, init_vff, compute_shifted_volume_properties

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

    SUBROUTINE compute_shifted_volume_properties(kk, jj, ii, vff, ddx, ddy, ddz, denistyFieldiStag, densityFieldjStag, densityFieldkStag)

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: vff(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(out) :: denistyFieldiStag(kk, jj, ii), densityFieldjStag(kk, jj, ii), densityFieldkStag(kk, jj, ii)

        ! Local variables
        LOGICAL :: isInterface(kk, jj, ii)
        REAL(realk) :: normx(kk, jj, ii), normy(kk, jj, ii), normz(kk, jj, ii)
        REAL(realk) :: alpha(kk, jj, ii)
        REAL(realk) :: vffiStag(kk, jj, ii), vffjStag(kk, jj, ii), vffkStag(kk, jj, ii)
        REAL(realk), PARAMETER :: tol = 1.0E-15

        CALL track_interface(isInterface, kk, jj, ii, vff, tol)
        CALL compute_normal_vector(normx, normy, normz, kk, jj, ii, vff, ddx, ddy, ddz, tol)
        CALL compute_alpha(alpha, kk, jj, ii, vff, isInterface, ddx, ddy, ddz, normx, normy, normz, tol)

        CALL compute_iStag_vff(kk, jj, ii, vffiStag, alpha, vff, isInterface, ddx, ddy, ddz, normx, normy, normz, tol)
        CALL compute_jStag_vff(kk, jj, ii, vffjStag, alpha, vff, isInterface, ddx, ddy, ddz, normx, normy, normz, tol)
        CALL compute_kStag_vff(kk, jj, ii, vffkStag, alpha, vff, isInterface, ddx, ddy, ddz, normx, normy, normz, tol)

        CALL get_material_property_field(kk, jj, ii, denistyFieldiStag, vffiStag, rho1, rho2)
        CALL get_material_property_field(kk, jj, ii, densityFieldjStag, vffjStag, rho1, rho2)
        CALL get_material_property_field(kk, jj, ii, densityFieldkStag, vffkStag, rho1, rho2)

    END SUBROUTINE compute_shifted_volume_properties

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