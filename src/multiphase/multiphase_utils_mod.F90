
!====================================================================
!  Module: mph_utils_mod
!
!   Responsibilities:
!   - 
!
!   Author:      Valentin Ruhs
!   e-Mail:      valentin.ruhs@gmx.de
!   Created:     2026-02
!   Last update: 2026-09
!
!====================================================================

MODULE mph_utils_mod

    USE MPI_f08
    USE grids_mod, ONLY: nmygrids, mygrids
    USE field_mod, ONLY: field_t
    USE fields_mod, ONLY: get_field
    USE grids_mod, ONLY: get_mgdims
    USE comms_mod, ONLY: myid
    USE precision_mod, ONLY: intk, realk, mglet_mpi_real
    USE mphcore_mod, ONLY: divTol, volTol
    USE mph_io_mod, ONLY: initVol, initErr, trueVol

    IMPLICIT NONE
    PRIVATE

    PUBLIC :: init_mph_utils, finish_mph_utils

CONTAINS

    SUBROUTINE init_mph_utils()

        ! Subroutine arguments
        ! None

        ! Local variables
        INTEGER(intk) :: n, igrid, igridf, ipar
        INTEGER(intk) :: kk, jj, ii, kc0, jc0, ic0
        REAL(realk), POINTER, CONTIGUOUS :: grdMask(:,:,:)

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_fieldptr(grdMask, "GRDMASK", igrid)
            grdMask = 1.0_realk
        END DO

        DO igridf = 1, ngrid
            ipar = iparent(igridf)

            IF (ipar == 0) CYCLE
            IF (idprocofgrd(ipar) /= myid) CYCLE

            CALL get_fieldptr(grdMask, "GRDMASK", ipar)
            CALL get_mgdims(kk, jj, ii, igridf)

            ic0 = iposition(igridf)
            jc0 = jposition(igridf)
            kc0 = kposition(igridf)

            grdMask(kc0:kc0+(kk-4)/2-1, &
                    jc0:jc0+(jj-4)/2-1, &
                    ic0:ic0+(ii-4)/2-1) = 0.0_realk
        END DO

    END SUBROUTINE init_mph_utils

    !================================================================

    SUBROUTINE finish_mph_utils()

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        CONTINUE

    END SUBROUTINE finish_mph_utils

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

    SUBROUTINE clip(kk, jj, ii, c, cClip)
    !----------------------------------------------------------------
    !   What it does:
    !   Clips c to its boundaries [0, 1].
    !
    !   Source:
    !   T. Arrufat et al., “A mass-momentum consistent, 
    !   Volume-of-Fluid method for incompressible flow on staggered 
    !   grids,” Computers & Fluids, vol. 215, p. 104785, Jan. 2021, 
    !   doi: 10.1016/j.compfluid.2020.104785.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: c(kk, jj, ii)
        REAL(realk), INTENT(out) :: cClip(kk, jj, ii)

        ! Local variables
        ! None

        cClip = MAX(MIN(c, 1.0_realk), 0.0_realk)

    END SUBROUTINE clip

    !================================================================

    SUBROUTINE check_solenoidality(dt)
    !----------------------------------------------------------------
    !   What it does:
    !    
    !----------------------------------------------------------------

        ! Subroutine arguments
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        REAL(realk), POINTER, CONTIGUOUS :: rddx(:), rddy(:), rddz(:)
        REAL(realk), POINTER, CONTIGUOUS :: u(:,:,:), v(:,:,:), w(:,:,:)
        INTEGER(intk) :: kk, jj, ii, k, j, i, n, igrid
        REAL(realk) :: div, divMax, divMaxGlob

        divMax = 0.0_realk
        divMaxGlob = 0.0_realk

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_fieldptr(u, "U", igrid)
            CALL get_fieldptr(v, "V", igrid)
            CALL get_fieldptr(w, "W", igrid)
            CALL get_fieldptr(rddx, "RDDX", igrid)
            CALL get_fieldptr(rddy, "RDDY", igrid)
            CALL get_fieldptr(rddz, "RDDZ", igrid)

            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        div = ( u(k,j,i) - u(k,j,i-1) )*rddx(i) + &
                              ( v(k,j,i) - v(k,j-1,i) )*rddy(j) + &
                              ( w(k,j,i) - w(k-1,j,i) )*rddz(k)
                        IF ( ABS(div) > divMax ) THEN
                            divMax = ABS(div)
                        ENDIF
                    ENDDO
                ENDDO
            ENDDO
        ENDDO

        CALL MPI_Allreduce(divMax, divMaxGlob, 1, mglet_mpi_real, MPI_MAX, MPI_COMM_WORLD)

        IF ( divMaxGlob*dt >= divTol ) THEN
            IF ( myid == 0 ) THEN
                WRITE(*,'(A,ES14.6,A,ES14.6)') "Solenoidality violated! max|div| = ", divMaxGlob, &
                    "  max|div|*dt = ", divMaxGlob*dt
            ENDIF
        ENDIF

    END SUBROUTINE check_solenoidality

    !================================================================

    SUBROUTINE check_continuity()
    !----------------------------------------------------------------
    !   What it does:
    !    
    !----------------------------------------------------------------

        ! Subroutine arguments

        ! Local variables
        REAL(realk), POINTER, CONTIGUOUS :: c(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        REAL(realk), POINTER, CONTIGUOUS :: grdMask(:,:,:)
        INTEGER(intk) :: kk, jj, ii, k, j, i, n, igrid
        REAL(realk) :: currVol, domaVol

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_fieldptr(c, "C", igrid)
            CALL get_fieldptr(ddx, "DDX", igrid)
            CALL get_fieldptr(ddy, "DDY", igrid)
            CALL get_fieldptr(ddz, "DDZ", igrid)
            CALL get_fieldptr(ddz, "DDZ", igrid)

            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        domaVol = domaVol + grdMask(k,j,i) * ddx(i) * ddy(j) * ddz(k)
                        currVol = currVol + grdMask(k,j,i) * c(k,j,i) * ddx(i) * ddy(j) * ddz(k)
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

    SUBROUTINE comp_vol_phase1(kk, jj, ii, c, ddx, ddy, ddz, volPhase1)
    !----------------------------------------------------------------
    !   What it does:
    !    
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: c(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(inout) :: volPhase1

        ! Local variables
        INTEGER(intk) :: k, j, i

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    volPhase1 = volPhase1 + c(k,j,i) * ddx(i) * ddy(j) * ddz(k)
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

    !================================================================

    SUBROUTINE comp_matrix_coeff_mph()
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        ! None

        ! Local variables
        TYPE(field_t), POINTER :: c_f
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: n, igrid

        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ae(:,:,:), aw(:,:,:), &
                                            an(:,:,:), as(:,:,:), &
                                            at(:,:,:), ab(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: ap(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: c(:, :, :)
        REAL(realk), ALLOCATABLE :: rhoe(:, :, :), rhon(:, :, :), rhot(:, :, :)

        CALL get_field(c_f, "c")

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_fieldptr(dx, "DX", igrid)
            CALL get_fieldptr(dy, "DY", igrid)
            CALL get_fieldptr(dz, "DZ", igrid)

            CALL get_fieldptr(rdx, "RDX", igrid)
            CALL get_fieldptr(rdy, "RDY", igrid)
            CALL get_fieldptr(rdz, "RDZ", igrid)

            CALL get_fieldptr(aw, "GSAW", igrid)
            CALL get_fieldptr(ae, "GSAE", igrid)
            CALL get_fieldptr(as, "GSAS", igrid)
            CALL get_fieldptr(an, "GSAN", igrid)
            CALL get_fieldptr(ab, "GSAB", igrid)
            CALL get_fieldptr(at, "GSAT", igrid)

            CALL get_fieldptr(ap, "GSAP", igrid)

            CALL c_f%get_ptr(c, igrid)

            IF ( .NOT. ALLOCATED(rhoe)) ALLOCATE(rhoe(kk, jj, ii))
            IF ( .NOT. ALLOCATED(rhon)) ALLOCATE(rhon(kk, jj, ii))
            IF ( .NOT. ALLOCATED(rhot)) ALLOCATE(rhot(kk, jj, ii))

            CALL comp_prop_face(kk, jj, ii, c, rhoe, rhon, rhot, rho1, rho2, rdx, rdy, rdz)

            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        ae(k,j,i) = 2.0/((dx(i-1)+dx(i))*dx(i)*rhoe(k,j,i))
                        aw(k,j,i) = 2.0/((dx(i-1)+dx(i))*dx(i-1)*rhoe(k,j,i-1))
                        an(k,j,i) = 2.0/((dy(j-1)+dy(j))*dy(j)*rhon(k,j,i))
                        as(k,j,i) = 2.0/((dy(j-1)+dy(j))*dy(j-1)*rhon(k,j-1,i))
                        at(k,j,i) = 2.0/((dz(k-1)+dz(k))*dz(k)*rhot(k,j,i))
                        ab(k,j,i) = 2.0/((dz(k-1)+dz(k))*dz(k-1)*rhot(k-1,j,i))
                    ENDDO
                ENDDO
            ENDDO

            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        ap(k, j, i) = -( ae(k,j,i) + aw(k,j,i) + an(k,j,i) &
                                       + as(k,j,i) + at(k,j,i) + ab(k,j,i) )
                    END DO
                END DO
            END DO

            IF ( ALLOCATED(rhot)) DEALLOCATE(rhot)
            IF ( ALLOCATED(rhon)) DEALLOCATE(rhon)
            IF ( ALLOCATED(rhoe)) DEALLOCATE(rhoe)
            
        ENDDO

    END SUBROUTINE comp_matrix_coeff_mph

END MODULE mph_utils_mod