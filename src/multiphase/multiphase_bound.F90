!====================================================================
!  Module: multiphase_bound_mod
!
!  Responsibilities:
!     - 
!
!  Author:      Valentin Ruhs
!  Created:     2026-02
!  Last update: 2026-02
!
!====================================================================

MODULE multiphase_bound_mod

    USE bound_mod, ONLY: bound_t

    IMPLICIT NONE
    PRIVATE 

    ! Bound operation 'T' operate on vff
    TYPE, EXTENDS(bound_t) :: bound_multiphase_t
    CONTAINS
        PROCEDURE, NOPASS :: front => bfront
        PROCEDURE, NOPASS :: back => bfront
        PROCEDURE, NOPASS :: right => bright
        PROCEDURE, NOPASS :: left => bright
        PROCEDURE, NOPASS :: bottom => bbottom
        PROCEDURE, NOPASS :: top => bbottom
    END TYPE bound_multiphase_t
    TYPE(bound_multiphase_t) :: bound_multiphase

    PUBLIC :: bound_multiphase

CONTAINS

    SUBROUTINE init_multiphase_bound()

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        continue
    END SUBROUTINE init_multiphase_bound

    !================================================================

    SUBROUTINE finish_multiphase_bound()

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        continue
    END SUBROUTINE finish_multiphase_bound

    !================================================================

    

END MODULE multiphase_bound_mod