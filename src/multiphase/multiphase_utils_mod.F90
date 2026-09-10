
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

    USE MPI_f08
    USE grids_mod, ONLY: nmygrids, mygrids
    USE field_mod, ONLY: field_t
    USE fields_mod, ONLY: get_field
    USE grids_mod, ONLY: get_mgdims
    USE comms_mod, ONLY: myid
    USE precision_mod, ONLY: intk, realk, mglet_mpi_real
    USE multiphasecore_mod, ONLY: tol
    USE multiphase_io_mod, ONLY: initVol, initErr, trueVol

    IMPLICIT NONE
    PRIVATE

    REAL(realk), PROTECTED :: currErr = 0.0_realk
    REAL(realk), PROTECTED :: relaErr = 0.0_realk

    REAL(realk), PROTECTED :: cumFlux = 0.0_realk, cumComp = 0.0_realk
    REAL(realk), PROTECTED :: cumClip = 0.0_realk, cumResi = 0.0_realk
    REAL(realk), PROTECTED :: volRef  = 0.0_realk
    INTEGER(intk), PROTECTED :: nBal  = 0

    PUBLIC :: init_multiphase_utils, finish_multiphase_utils, &
        get_spatial_indices, get_spatial_extents, get_condit_velocity, clip_vff, &
        check_continuity, check_solenoidality, comp_vol_phase1, sanity_check

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

    SUBROUTINE get_spatial_indices(lOrq, io, jo, ko)
    !----------------------------------------------------------------
    !   What it does:
    !   Depending on the direction (l) of component (q) the 
    !   indices io, jo and ko are set to 0 or 1. 
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: lOrq
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
        ENDIF

    END SUBROUTINE get_condit_velocity

    !================================================================

    SUBROUTINE clip_vff(kk, jj, ii, vff, ddx, ddy, ddz, volClipPhase1)
    !----------------------------------------------------------------
    !   What it does:
    !   Clips the volume fraction field to its boundaries [0, 1].
    !    
    !   Source:
    !   T. Arrufat et al., “A mass-momentum consistent, 
    !   Volume-of-Fluid method for incompressible flow on staggered 
    !   grids,” Computers & Fluids, vol. 215, p. 104785, Jan. 2021, 
    !   doi: 10.1016/j.compfluid.2020.104785.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: vff(kk, jj, ii)
        REAL(realk), INTENT(in), OPTIONAL :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(inout), OPTIONAL :: volClipPhase1

        ! Local variables
        INTEGER(intk) :: k, j, i
        REAL(realk) :: volBeforClip, volAfterClip, vol

        volBeforClip = 0.0_realk
        volAfterClip = 0.0_realk

        IF ( PRESENT(volClipPhase1) ) THEN
            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        vol = ddx(i) * ddy(j) * ddz(k)
                        volBeforClip = volBeforClip + vff(k,j,i) * vol
                    ENDDO
                ENDDO
            ENDDO
        ENDIF

        DO i = 1, ii
            DO j = 1, jj
                DO k = 1, kk
                    IF ( vff(k,j,i) <= tol ) THEN
                        vff(k,j,i) = 0.0_realk
                    ELSE IF ( vff(k,j,i) >= ( 1.0_realk - tol ) ) THEN
                        vff(k,j,i) = 1.0_realk
                    ENDIF
                ENDDO 
            ENDDO
        ENDDO

        IF ( PRESENT(volClipPhase1) ) THEN
            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        vol = ddx(i) * ddy(j) * ddz(k)
                        volAfterClip = volAfterClip + vff(k,j,i) * vol
                    ENDDO
                ENDDO
            ENDDO
            volClipPhase1 = volClipPhase1 + ( volAfterClip - volBeforClip )
        ENDIF

    END SUBROUTINE clip_vff

    !================================================================

    SUBROUTINE check_solenoidality(itstep, dt)
    !----------------------------------------------------------------
    !   What it does:
    !    
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: itstep
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        TYPE(field_t), POINTER :: u_f, v_f, w_f
        TYPE(field_t), POINTER :: ddx_f, ddy_f, ddz_f
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        REAL(realk), POINTER, CONTIGUOUS :: u(:,:,:), v(:,:,:), w(:,:,:)
        INTEGER(intk) :: kk, jj, ii, k, j, i, n, igrid
        REAL(realk) :: div, divMax, divMaxGlob

        div = 0.0_realk
        divMax = 0.0_realk
        divMaxGlob = 0.0_realk

        CALL get_field(u_f, "U")
        CALL get_field(v_f, "V")
        CALL get_field(w_f, "W")
        CALL get_field(ddx_f, "DDX")
        CALL get_field(ddy_f, "DDY")
        CALL get_field(ddz_f, "DDZ")

        DO n = 1, nmygrids
            igrid = mygrids(n)
            CALL get_mgdims(kk, jj, ii, igrid)

            CALL u_f%get_ptr(u, igrid)
            CALL v_f%get_ptr(v, igrid)
            CALL w_f%get_ptr(w, igrid)
            CALL ddx_f%get_ptr(ddx, igrid)
            CALL ddy_f%get_ptr(ddy, igrid)
            CALL ddz_f%get_ptr(ddz, igrid)

            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        div = ( u(k,j,i) - u(k,j,i-1) ) / ddx(i) + &
                              ( v(k,j,i) - v(k,j-1,i) ) / ddy(j) + &
                              ( w(k,j,i) - w(k-1,j,i) ) / ddz(k)
                        IF ( ABS(div) > divMax ) THEN
                            divMax = ABS(div)
                        ENDIF
                    ENDDO
                ENDDO
            ENDDO
        ENDDO

        CALL MPI_Allreduce(divMax, divMaxGlob, 1, mglet_mpi_real, MPI_MAX, MPI_COMM_WORLD)

        IF ( divMaxGlob * dt >= tol ) THEN
            IF ( myid == 0 ) THEN
                WRITE(*,'(A,ES14.6,A,ES14.6)') "Solenoidality violated! max|div| = ", divMaxGlob, &
                    "  max|div|*dt = ", divMaxGlob*dt
            ENDIF
        ENDIF

    END SUBROUTINE check_solenoidality

    !================================================================

    SUBROUTINE check_continuity(itstep)
    !----------------------------------------------------------------
    !   What it does:
    !    
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: itstep

        ! Local variables
        TYPE(field_t), POINTER :: vff_f
        TYPE(field_t), POINTER :: ddx_f, ddy_f, ddz_f
        TYPE(field_t), POINTER :: grdMask_f
        REAL(realk), POINTER, CONTIGUOUS :: vff(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        REAL(realk), POINTER, CONTIGUOUS :: grdMask(:,:,:)
        INTEGER(intk) :: kk, jj, ii, k, j, i, n, igrid
        REAL(realk) :: currVol, domaVol

        currVol = 0.0_realk
        domaVol = 0.0_realk

        CALL get_field(vff_f, "VFF")
        CALL get_field(ddx_f, "DDX")
        CALL get_field(ddy_f, "DDY")
        CALL get_field(ddz_f, "DDZ")
        CALL get_field(grdMask_f, "GRDMASK")

        DO n = 1, nmygrids
            igrid = mygrids(n)
            CALL get_mgdims(kk, jj, ii, igrid)

            CALL vff_f%get_ptr(vff, igrid)
            CALL ddx_f%get_ptr(ddx, igrid)
            CALL ddy_f%get_ptr(ddy, igrid)
            CALL ddz_f%get_ptr(ddz, igrid)
            CALL grdMask_f%get_ptr(grdMask, igrid)

            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        domaVol = domaVol + grdMask(k,j,i) * ddx(i) * ddy(j) * ddz(k)
                        currVol = currVol + grdMask(k,j,i) * vff(k,j,i) * ddx(i) * ddy(j) * ddz(k)
                    ENDDO
                ENDDO
            ENDDO
        ENDDO

        CALL MPI_Allreduce(MPI_IN_PLACE, currVol, 1, mglet_mpi_real, MPI_SUM, MPI_COMM_WORLD)
        CALL MPI_Allreduce(MPI_IN_PLACE, domaVol, 1, mglet_mpi_real, MPI_SUM, MPI_COMM_WORLD)

        currErr = trueVol - currVol
        relaErr = ( currVol - initVol ) / initVol

        IF ( myid == 0 .AND. ABS(relaErr) >= tol ) THEN
            WRITE(*,'(A,ES14.6)') "Continuity violated! Relative volume error of ", relaErr
        ENDIF

    END SUBROUTINE check_continuity

    !================================================================

    SUBROUTINE comp_vol_phase1(kk, jj, ii, vff, ddx, ddy, ddz, volPhase1)
    !----------------------------------------------------------------
    !   What it does:
    !    
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: vff(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(inout) :: volPhase1

        ! Local variables
        INTEGER(intk) :: k, j, i

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    volPhase1 = volPhase1 + vff(k,j,i) * ddx(i) * ddy(j) * ddz(k)
                ENDDO
            ENDDO
        ENDDO

    END SUBROUTINE comp_vol_phase1

    !================================================================

    SUBROUTINE sanity_check(volPhase1r, volPhase1r1, volFluxPhase1, volCompPhase1, volClipPhase1)
    !----------------------------------------------------------------
    !   What it does:
    !    
    !----------------------------------------------------------------

        ! Subroutine arguments
        REAL(realk), INTENT(in) :: volPhase1r, volPhase1r1
        REAL(realk), INTENT(in) :: volFluxPhase1, volCompPhase1, volClipPhase1

        ! Local variables
        REAL(realk) :: volResi, volTres

        IF ( nBal == 0 ) volRef = volPhase1r

        volResi = ( volPhase1r1 - volPhase1r ) - ( volFluxPhase1 + volCompPhase1 + volClipPhase1 )

        cumFlux = cumFlux + volFluxPhase1
        cumComp = cumComp + volCompPhase1
        cumClip = cumClip + volClipPhase1
        cumResi = cumResi + volResi
        nBal    = nBal + 1

        volTres = ( volPhase1r1 - volRef ) - ( cumFlux + cumComp + cumClip + cumResi )

        IF ( myid == 0 ) THEN
            WRITE(*,'(A,6(A,ES14.6))') "VOF acc. budget: ", &
            "  dV/V = ", ( volPhase1r1 - volRef ) / volRef, &
            "  flux/V = ", cumFlux / volRef, &
            "  comp/V = ", cumComp / volRef, &
            "  clip/V = ", cumClip / volRef, &
            "  resi/V = ", cumResi / volRef, &
            "  tres/V = ", volTres / volRef
        ENDIF

    END SUBROUTINE sanity_check

END MODULE multiphase_utils_mod