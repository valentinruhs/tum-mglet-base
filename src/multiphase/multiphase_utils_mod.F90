
!====================================================================
!  Module: multiphase_utils_mod
!
!  Responsibilities:
!     - 
!
!  Author:      Valentin Ruhs
!  Created:     2026-02
!  Last update: 2026-02
!
!====================================================================

MODULE multiphase_utils_mod

    USE precision_mod, ONLY: intk, realk
    USE err_mod, ONLY: errr
    
    IMPLICIT NONE
    PRIVATE 

    PUBLIC :: init_multiphase_utils, finish_multiphase_utils, get_spatial_indices, get_spatial_extents, get_condit_velocity

CONTAINS

    SUBROUTINE init_multiphase_utils()

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        continue
    END SUBROUTINE init_multiphase_utils

    !================================================================

    SUBROUTINE finish_multiphase_utils()

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None
        
        continue
    END SUBROUTINE finish_multiphase_utils

    !================================================================

    SUBROUTINE get_spatial_indices(kk, jj, ii, lOrq, io, jo, ko)
    !----------------------------------------------------------------
    !   What it does:
    !   Depending on the direction (l) of component (q) the 
    !   indices io, jo and ko are set to 0 or 1. 
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, lOrq
        INTEGER(intk), INTENT(out) :: io, jo, ko

        ! Local variables
        ! None

        io = 0_intk ; jo = 0_intk ; ko = 0_intk

        IF ( lOrq == 1 ) THEN
            io = 1_intk
        ELSE IF ( lOrq == 2 ) THEN
            jo = 1_intk
        ELSE IF ( lOrq == 3 ) THEN
            ko = 1_intk
        ELSE
            CALL errr(__FILE__, __LINE__)
        END IF

    END SUBROUTINE get_spatial_indices

    !================================================================

    SUBROUTINE get_spatial_extents(kk, jj, ii, q, l, dx, dy, dz, ddx, ddy, ddz, dsx, dsy, dsz)
    !----------------------------------------------------------------
    !   What it does:
    !   Depending on the direction (l) of component (q) the 
    !   extents ds(.) are set to d(.) or dd(.). The term ds(.) stands
    !   for spacing in (.)-direction and is a neutral specification
    !   for face-to-face or center-to-center distance.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, q, l
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(out) :: dsx(ii), dsy(jj), dsz(kk)

        ! Local variables
        ! None

        dsx = ddx ; dsy = ddy ; dsz = ddz

        IF ( q == l ) THEN
            IF ( l == 1 ) THEN
                dsx = dx
            ELSEIF ( l == 2 ) THEN
                dsy = dy
            ELSEIF ( l == 3 ) THEN
                dsz = dz
            ELSE
                CALL errr(__FILE__, __LINE__)
            ENDIF
        ENDIF

    END SUBROUTINE get_spatial_extents

    !================================================================

    SUBROUTINE get_condit_velocity(kk, jj, ii, lOrq, u, v, w, vel)
    !----------------------------------------------------------------
    !   What it does:
    !   Depending on the direction (l) of component (q) the 
    !   specific velocity (vel) is set to u, v or w.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, lOrq
        REAL(realk), INTENT(in) :: u(kk,jj,ii), v(kk,jj,ii), w(kk,jj,ii)
        REAL(realk), INTENT(out) :: vel(kk,jj,ii)

        ! Local variables
        ! None

        vel = 0.0_realk

        IF ( lOrq == 1 ) THEN
            vel = u
        ELSE IF ( lOrq == 2 ) THEN
            vel = v
        ELSE IF ( lOrq == 3 ) THEN
            vel = w
        ELSE
            CALL errr(__FILE__, __LINE__)
        END IF

    END SUBROUTINE get_condit_velocity

END MODULE multiphase_utils_mod