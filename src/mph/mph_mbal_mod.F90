!====================================================================
!  Module: mph_mbal_mod
!
!   Responsibilities:
!   - 
!
!   Author:      Valentin Ruhs
!   e-Mail:      valentin.ruhs@gmx.de
!   Created:     2026-09
!   Last update: 2026-09
!
!====================================================================

MODULE mph_mbal_mod

    USE MPI_f08
    USE precision_mod, ONLY: intk, realk, mglet_mpi_real
    USE comms_mod, ONLY: myid
    USE grids_mod, ONLY: nmygrids, mygrids, get_mgdims, get_mgbasb, &
        level, maxlevel
    USE fields_mod, ONLY: get_fieldptr
    USE fort7_mod, ONLY: fort7
    USE config_mod, ONLY: config_t

    USE mphcore_mod, ONLY: rho1, grav
    USE mph_utils_mod, ONLY: get_lp_bnds
    USE mph_plic_mod, ONLY: comp_c_stg
    USE mph_props_mod, ONLY: comp_d_stg

    IMPLICIT NONE(type, external)
    PRIVATE

    LOGICAL, PROTECTED :: mbalOn = .FALSE.

    REAL(realk) :: m0 = 0.0_realk, m1 = 0.0_realk, vLiq = 1.0_realk
    LOGICAL :: havePrev = .FALSE.
    INTEGER(intk) :: stepPrev = 0
    REAL(realk) :: dtPrev = 0.0_realk, m2Prev = 0.0_realk
    REAL(realk) :: dmResPrev = 0.0_realk, dmAdvPrev = 0.0_realk
    REAL(realk) :: fDrvPrev = 0.0_realk, fGrvPrev = 0.0_realk, &
        fWalPrev = 0.0_realk
    REAL(realk) :: cumRes = 0.0_realk, cumAdv = 0.0_realk, &
        cumRk = 0.0_realk, cumT = 0.0_realk
    LOGICAL :: headerDone = .FALSE.
    LOGICAL :: firstDone = .FALSE.

    PUBLIC :: init_mph_mbal, finish_mph_mbal, mbal_pre, mbal_post, mbalOn

CONTAINS

    SUBROUTINE init_mph_mbal()

        ! Subroutine arguments
        ! None

        ! Local variables
        TYPE(config_t) :: mphConf

        CALL fort7%get(mphConf, "/multiphase")
        CALL mphConf%get_value("/momentumBalance", mbalOn, .FALSE.)

        IF (mbalOn .AND. myid == 0) THEN
            WRITE(*, '(A)') "MBAL: x-momentum budget active (finest level only)"
            WRITE(*, '()')
        END IF

    END SUBROUTINE init_mph_mbal

    !================================================================

    SUBROUTINE finish_mph_mbal()

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        CONTINUE

    END SUBROUTINE finish_mph_mbal

    !================================================================

    SUBROUTINE mbal_pre()
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        ! None

        ! Local variables
        REAL(realk) :: fWal, dmRk, aNorm, cNorm

        IF (.NOT. mbalOn) RETURN

        m0 = sum_u("DBA", "U")
        vLiq = sum_c()
        fWal = wall_force()

        IF (havePrev) THEN
            dmRk = (m0 - m2Prev) - dtPrev*(fDrvPrev + fGrvPrev &
                - 0.5_realk*(fWalPrev + fWal))

            cumRes = cumRes + dmResPrev
            cumAdv = cumAdv + dmAdvPrev
            cumRk = cumRk + dmRk
            cumT = cumT + dtPrev

            IF (myid == 0) THEN
                aNorm = 1.0_realk/(dtPrev*rho1*vLiq)
                cNorm = 1.0_realk/(cumT*rho1*vLiq)
                IF (.NOT. headerDone) THEN
                    WRITE(*, '(A)') "MBAL     step         aRes         aAdv" // &
                        "          aRk       aDrive        aWall" // &
                        "       <aRes>       <aAdv>        <aRk>"
                    headerDone = .TRUE.
                END IF
                WRITE(*, '(A,I8,8ES13.4)') "MBAL ", stepPrev, &
                    dmResPrev*aNorm, dmAdvPrev*aNorm, dmRk*aNorm, &
                    fDrvPrev/(rho1*vLiq), &
                    0.5_realk*(fWalPrev + fWal)/(rho1*vLiq), &
                    cumRes*cNorm, cumAdv*cNorm, cumRk*cNorm
            END IF
        END IF

        CALL comp_c_stg(1)
        CALL comp_d_stg(1)
        m1 = sum_u("DS1", "U")

    END SUBROUTINE mbal_pre

    !================================================================

    SUBROUTINE mbal_post(itstep, dt, gradp)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: itstep
        REAL(realk), INTENT(in) :: dt
        REAL(realk), INTENT(in) :: gradp(3)

        ! Local variables
        REAL(realk) :: m2

        IF (.NOT. mbalOn) RETURN

        m2 = sum_u("DBA", "U")

        IF (firstDone) THEN
            dmResPrev = m1 - m0
        ELSE
            dmResPrev = 0.0_realk
            firstDone = .TRUE.
        END IF
        dmAdvPrev = m2 - m1

        fDrvPrev = -gradp(1)*sum_u("CS1")
        fGrvPrev = grav(1)*sum_u("DBA")
        fWalPrev = wall_force()

        m2Prev = m2
        dtPrev = dt
        stepPrev = itstep
        havePrev = .TRUE.

    END SUBROUTINE mbal_post

    !================================================================

    FUNCTION sum_u(name1, name2) RESULT(s)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Function arguments
        CHARACTER(len=*), INTENT(in) :: name1
        CHARACTER(len=*), INTENT(in), OPTIONAL :: name2
        REAL(realk) :: s

        ! Local variables
        INTEGER(intk) :: n, igrid, i, j, k
        INTEGER(intk) :: ista, iend, jsta, jend, ksta, kend
        REAL(realk) :: sLoc(1), sGlo(1)
        REAL(realk), POINTER, CONTIGUOUS :: f1(:,:,:), f2(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), ddy(:), ddz(:)

        sLoc = 0.0_realk
        DO n = 1, nmygrids
            igrid = mygrids(n)
            IF (level(igrid) /= maxlevel) CYCLE

            CALL get_fieldptr(f1, name1, igrid)
            IF (PRESENT(name2)) CALL get_fieldptr(f2, name2, igrid)
            CALL get_fieldptr(dx, "DX", igrid)
            CALL get_fieldptr(ddy, "DDY", igrid)
            CALL get_fieldptr(ddz, "DDZ", igrid)

            CALL get_lp_bnds(igrid, 1, ista, iend, jsta, jend, ksta, kend)

            IF (PRESENT(name2)) THEN
                DO i = ista, iend
                    DO j = jsta, jend
                        DO k = ksta, kend
                            sLoc(1) = sLoc(1) + f1(k,j,i)*f2(k,j,i) &
                                *dx(i)*ddy(j)*ddz(k)
                        END DO
                    END DO
                END DO
            ELSE
                DO i = ista, iend
                    DO j = jsta, jend
                        DO k = ksta, kend
                            sLoc(1) = sLoc(1) + f1(k,j,i)*dx(i)*ddy(j)*ddz(k)
                        END DO
                    END DO
                END DO
            END IF
        END DO

        CALL MPI_Allreduce(sLoc, sGlo, 1, mglet_mpi_real, MPI_SUM, MPI_COMM_WORLD)
        s = sGlo(1)

    END FUNCTION sum_u

    !================================================================

    FUNCTION sum_c() RESULT(s)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Function arguments
        REAL(realk) :: s

        ! Local variables
        INTEGER(intk) :: n, igrid, i, j, k, kk, jj, ii
        REAL(realk) :: sLoc(1), sGlo(1)
        REAL(realk), POINTER, CONTIGUOUS :: c(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)

        sLoc = 0.0_realk
        DO n = 1, nmygrids
            igrid = mygrids(n)
            IF (level(igrid) /= maxlevel) CYCLE

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_fieldptr(c, "C", igrid)
            CALL get_fieldptr(ddx, "DDX", igrid)
            CALL get_fieldptr(ddy, "DDY", igrid)
            CALL get_fieldptr(ddz, "DDZ", igrid)

            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        sLoc(1) = sLoc(1) + c(k,j,i)*ddx(i)*ddy(j)*ddz(k)
                    END DO
                END DO
            END DO
        END DO

        CALL MPI_Allreduce(sLoc, sGlo, 1, mglet_mpi_real, MPI_SUM, MPI_COMM_WORLD)
        s = sGlo(1)

    END FUNCTION sum_c

    !================================================================

    FUNCTION wall_force() RESULT(f)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Function arguments
        REAL(realk) :: f

        ! Local variables
        INTEGER(intk) :: n, igrid, i, j, k, kk, jj, ii
        INTEGER(intk) :: nfro, nbac, nrgt, nlft, nbot, ntop
        INTEGER(intk) :: ista, iend, jsta, jend, ksta, kend
        REAL(realk) :: fLoc(1), fGlo(1)
        REAL(realk), POINTER, CONTIGUOUS :: u(:,:,:), gUv(:,:,:), gUw(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), ddy(:), ddz(:)

        fLoc = 0.0_realk
        DO n = 1, nmygrids
            igrid = mygrids(n)
            IF (level(igrid) /= maxlevel) CYCLE

            CALL get_mgbasb(nfro, nbac, nrgt, nlft, nbot, ntop, igrid)
            IF (ALL([nrgt, nlft, nbot, ntop] /= 5)) CYCLE

            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_fieldptr(u, "U", igrid)
            CALL get_fieldptr(gUv, "GUV", igrid)
            CALL get_fieldptr(gUw, "GUW", igrid)
            CALL get_fieldptr(dx, "DX", igrid)
            CALL get_fieldptr(ddy, "DDY", igrid)
            CALL get_fieldptr(ddz, "DDZ", igrid)

            CALL get_lp_bnds(igrid, 1, ista, iend, jsta, jend, ksta, kend)

            ! Right
            IF (nrgt == 5) THEN
                j = 3
                DO i = ista, iend
                    DO k = ksta, kend
                        fLoc(1) = fLoc(1) + gUv(k,j-1,i)*2.0_realk*u(k,j,i) &
                            /ddy(j)*dx(i)*ddz(k)
                    END DO
                END DO
            END IF

            ! Left
            IF (nlft == 5) THEN
                j = jj-2
                DO i = ista, iend
                    DO k = ksta, kend
                        fLoc(1) = fLoc(1) + gUv(k,j,i)*2.0_realk*u(k,j,i) &
                            /ddy(j)*dx(i)*ddz(k)
                    END DO
                END DO
            END IF

            ! Bottom
            IF (nbot == 5) THEN
                k = 3
                DO i = ista, iend
                    DO j = jsta, jend
                        fLoc(1) = fLoc(1) + gUw(k-1,j,i)*2.0_realk*u(k,j,i) &
                            /ddz(k)*dx(i)*ddy(j)
                    END DO
                END DO
            END IF

            ! Top
            IF (ntop == 5) THEN
                k = kk-2
                DO i = ista, iend
                    DO j = jsta, jend
                        fLoc(1) = fLoc(1) + gUw(k,j,i)*2.0_realk*u(k,j,i) &
                            /ddz(k)*dx(i)*ddy(j)
                    END DO
                END DO
            END IF
        END DO

        CALL MPI_Allreduce(fLoc, fGlo, 1, mglet_mpi_real, MPI_SUM, MPI_COMM_WORLD)
        f = fGlo(1)

    END FUNCTION wall_force

END MODULE mph_mbal_mod