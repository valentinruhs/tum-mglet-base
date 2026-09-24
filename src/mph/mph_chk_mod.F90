!====================================================================
!  Module: mph_chk_mod
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

MODULE mph_chk_mod

    USE MPI_f08
    USE grids_mod, ONLY: nmygrids, mygrids, get_mgdims, iparent, &
        idprocofgrd, iposition, jposition, kposition, ngrid, get_bbox
    USE comms_mod, ONLY: myid
    USE precision_mod, ONLY: intk, realk, mglet_mpi_real
    USE fields_mod, ONLY: get_fieldptr, set_field
    USE err_mod, ONLY: err_abort
    USE comms_mod, ONLY: myid
    USE grids_mod, ONLY: minlevel, maxlevel

    USE mphcore_mod, ONLY: volChk, vofChk, volTol, vofTol, hasMph
    USE mph_test_mod, ONLY: comp_abs_res, circumf, area, shape

    IMPLICIT NONE(type, external)
    PRIVATE

    REAL(realk) :: volInit

    PUBLIC :: init_mph_chk, finish_mph_chk, comp_vol, final_chk, itinfo_mph

CONTAINS

    SUBROUTINE init_mph_chk()

        ! Subroutine arguments
        ! None

        ! Local variables
        CHARACTER(len=*), PARAMETER :: descGrdmask = "uncov. cells"
        CHARACTER(len=*), PARAMETER :: descCInit = "init. C"
        INTEGER(intk) :: n, igrid, igridf, ipar
        INTEGER(intk) :: kk, jj, ii, kc0, jc0, ic0
        REAL(realk), POINTER, CONTIGUOUS :: grdMask(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: c(:,:,:), cInit(:,:,:)

        CALL set_field("GRDMASK", description=descGrdmask, buffers=.TRUE.)
        CALL set_field("CINIT", description=descCInit, buffers=.TRUE.)

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_fieldptr(grdMask, "GRDMASK", igrid)
            CALL get_fieldptr(c, "C", igrid)
            CALL get_fieldptr(cInit, "CINIT", igrid)

            grdMask = 1.0_realk
            cInit = c
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

        CALL comp_vol(volInit)
        WRITE(*,'(A,E11.5)') "Initial volume is: ", volInit
        WRITE(*, '()')

    END SUBROUTINE init_mph_chk

    !================================================================

    SUBROUTINE finish_mph_chk()

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        CONTINUE

    END SUBROUTINE finish_mph_chk

    !================================================================

    SUBROUTINE itinfo_mph()
    !----------------------------------------------------------------
    !   What it does:
    !
    !----------------------------------------------------------------

        ! Subroutine arguments
        ! None

        ! Local variables
        REAL(realk) :: volCurr
        INTEGER(intk) :: i

        IF (.NOT. hasMph) RETURN

        CALL comp_vol(volCurr)

        IF (myid == 0) THEN
            DO i = minlevel, maxlevel
                WRITE(*, '(A,A,E20.10,E20.10)') &
                    "ABSVOLERR, ", "RELVOLERR: ", &
                    ABS(volInit-volCurr), ABS(volInit-volCurr)/volInit
            END DO
        END IF

    END SUBROUTINE itinfo_mph

    !================================================================

    SUBROUTINE final_chk()
    !----------------------------------------------------------------
    !   What it does:
    !
    !----------------------------------------------------------------

        ! Subroutine arguments
        ! None

        ! Local variables
        REAL(realk) :: volFini, L1, L2

        CALL comp_vol(volFini)
        CALL comp_L_norm_err(L1, L2)

        IF ( myid == 0 ) THEN
            WRITE(*, '()')
            WRITE(*,'(A,E11.5)') "Final volume is: ", volFini
            WRITE(*,'(A,E11.5)') "Absolute volume error: ", ABS(volInit-volFini)
            WRITE(*,'(A,E11.5)') "Relative volume error: ", ABS(volInit-volFini)/volInit
            WRITE(*,'(A,E11.5)') "E_{l=1,B=1}: ", L1
            WRITE(*,'(A,E11.5)') "E_{l=1,B=L}: ", L1/circumf
            WRITE(*, '()')
        END IF

    END SUBROUTINE final_chk

    !================================================================

    SUBROUTINE comp_L_norm_err(L1, L2)
    !----------------------------------------------------------------
    !   What it does:
    !
    !----------------------------------------------------------------

        ! Subroutine arguments
        REAL(realk), INTENT(out) :: L1, L2

        ! Local variables
        INTEGER(intk) :: n, igrid
        INTEGER(intk) :: kk, jj, ii
        REAL(realk) :: minx, maxx, miny, maxy, minz, maxz
        REAL(realk), POINTER, CONTIGUOUS :: c(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        REAL(realk), POINTER, CONTIGUOUS :: normx(:,:,:), normy(:,:,:), normz(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: alpha(:,:,:)
        REAl(realk), POINTER, CONTIGUOUS :: isIfc(:,:,:)
        REAL(realk) :: L1Loc, L2Loc, L1Grd, L2Grd

        L1Loc = 0.0_realk
        L2Loc = 0.0_realk
        DO n = 1, nmygrids
            igrid = mygrids(n)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_bbox(minx, maxx, miny, maxy, minz, maxz, igrid)
            CALL get_fieldptr(c, "C", igrid)
            CALL get_fieldptr(ddx, "DDX", igrid)
            CALL get_fieldptr(ddy, "DDY", igrid)
            CALL get_fieldptr(ddz, "DDZ", igrid)
            CALL get_fieldptr(normx, "NORMX", igrid)
            CALL get_fieldptr(normy, "NORMY", igrid)
            CALL get_fieldptr(normz, "NORMZ", igrid)
            CALL get_fieldptr(alpha, "ALPHA", igrid)
            CALL get_fieldptr(isIfc, "ISIFC", igrid)

            CALL comp_L_norm_err_grd(kk, jj, ii, minx, miny, minz, c, &
                normx, normy, normz, alpha, isIfc, ddx, ddy, ddz, L1Grd, L2Grd)
            L1Loc = L1Loc + L1Grd
            L2Loc = L2Loc + L2Grd
        END DO

        CALL MPI_Allreduce(L1Loc, L1, 1, mglet_mpi_real, MPI_SUM, MPI_COMM_WORLD)
        CALL MPI_Allreduce(L2Loc, L2, 1, mglet_mpi_real, MPI_SUM, MPI_COMM_WORLD)

        L1 = L1
        L2 = SQRT(L2)

    END SUBROUTINE comp_L_norm_err

    !================================================================

    SUBROUTINE comp_L_norm_err_grd(kk, jj, ii, minx, miny, minz, c, &
        normx, normy, normz, alpha, isIfc, ddx, ddy, ddz, L1, L2)
    !----------------------------------------------------------------
    !   What it does:
    !
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: minx, miny, minz
        REAL(realk), INTENT(in) :: c(kk, jj, ii)
        REAL(realk), INTENT(in) :: normx(kk, jj, ii), normy(kk, jj, ii), normz(kk, jj, ii)
        REAL(realk), INTENT(in) :: alpha(kk, jj, ii), isIfc(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(out) :: L1, L2

        ! Local variables
        INTEGER(intk) :: k, j, i
        REAL(realk) :: xMi, yMi, zMi
        REAL(realk) :: vol, absRes, eps

        L1 = 0.0_realk
        L2 = 0.0_realk
        xMi = minx
        DO i = 3, ii-2
            yMi = miny
            DO j = 3, jj-2
                zMi = minz
                DO k = 3, kk-2
                    vol = ddx(i)*ddy(j)*ddz(k)
                    absRes = comp_abs_res(xMi, yMi, zMi, ddx(i), ddy(j), ddz(k), &
                        normx(k,j,i), normy(k,j,i), normz(k,j,i), alpha(k,j,i), &
                        c(k,j,i), isIfc(k,j,i), shape)
                    eps = absRes/vol
                    L1 = L1 + eps*vol
                    L2 = L2 + eps**2*vol
                zMi = zMi + ddz(k)
                END DO
                yMi = yMi + ddy(j)
            END DO
            xMi = xMi + ddx(i)
        END DO

    END SUBROUTINE comp_L_norm_err_grd

    !================================================================

    SUBROUTINE comp_vol(vol)
    !----------------------------------------------------------------
    !   What it does:
    !
    !----------------------------------------------------------------

        ! Subroutine arguments
        REAL(realk), INTENT(out) :: vol

        ! Local variables
        INTEGER(intk) :: n, igrid
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: c(:,:,:), grdMask(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        REAL(realk) :: volLoc, volFld1

        volLoc = 0.0_realk
        DO n = 1, nmygrids
            igrid = mygrids(n)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_fieldptr(c, "C", igrid)
            CALL get_fieldptr(grdMask, "GRDMASK", igrid)
            CALL get_fieldptr(ddx, "DDX", igrid)
            CALL get_fieldptr(ddy, "DDY", igrid)
            CALL get_fieldptr(ddz, "DDZ", igrid)

            CALL comp_vol_grd(kk, jj, ii, c, grdMask, ddx, ddy, ddz, volFld1)
            volLoc =volLoc + volFld1
        END DO

        CALL MPI_Allreduce(volLoc, vol, 1, mglet_mpi_real, MPI_SUM, MPI_COMM_WORLD)

    END SUBROUTINE comp_vol

    !================================================================

    SUBROUTINE comp_vol_grd(kk, jj, ii, c, grdMask, ddx, ddy, ddz, volFld1)
    !----------------------------------------------------------------
    !   What it does:
    !
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: c(kk, jj, ii), grdMask(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(out) :: volFld1

        ! Local variables
        INTEGER(intk) :: k, j, i
        REAL(realk) :: volCell

        volFld1 = 0.0_realk
        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    volCell = ddx(i)*ddy(j)*ddz(k)
                    volFld1 = volFld1 + c(k,j,i)*volCell*grdMask(k,j,i)
                END DO
            END DO
        END DO

    END SUBROUTINE comp_vol_grd

END MODULE mph_chk_mod