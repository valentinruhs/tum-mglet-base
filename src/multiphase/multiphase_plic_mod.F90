!====================================================================
!  Module: multiphase_plic_mod
!
!  Responsibilities:
!     - 
!
!  Author:      Valentin Ruhs
!  Created:     2026-02
!  Last update: 2026-02
!
!====================================================================

MODULE multiphase_plic_mod

    USE precision_mod, ONLY: intk, realk
    
    IMPLICIT NONE

CONTAINS

    SUBROUTINE init_multiphase_plic()

        continue

    END SUBROUTINE init_multiphase_plic

    !================================================================

    SUBROUTINE finish_multiphase_plic()

        continue

    END SUBROUTINE finish_multiphase_plic

    !================================================================

    SUBROUTINE track_interface(is_interface, kk, jj, ii, c)

        ! Subroutine arguments
        LOGICAL, INTENT(out) :: is_interface(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: c(kk, jj, ii)

        ! Local variables


        continue

    END SUBROUTINE track_interface

    !================================================================

    SUBROUTINE compute_normal_vector(normx, normy, normz, kk, jj, ii, c, ddx, ddy, ddz, tol)

        ! Subroutine arguments
        REAL(realk), INTENT(out) :: normx(kk, jj, ii), normy(kk, jj, ii), normz(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: c(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: tol

        ! Local variables
        

        continue

    END SUBROUTINE compute_normal_vector

    !================================================================

    SUBROUTINE compute_alpha(alpha, kk, jj, ii, c, is_interface, ddx, ddy, ddz, normx, normy, normz)

        ! Subroutine arguments
        REAL(realk), INTENT(out) :: alpha(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: c(kk, jj, ii)
        LOGICAL, INTENT(out) :: is_interface(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(out) :: normx(kk, jj, ii), normy(kk, jj, ii), normz(kk, jj, ii)

        ! Local variables
        

        continue

    END SUBROUTINE compute_alpha

END MODULE multiphase_plic_mod