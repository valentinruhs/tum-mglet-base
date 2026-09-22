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
        idprocofgrd, iposition, jposition, kposition, ngrid
    USE comms_mod, ONLY: myid
    USE precision_mod, ONLY: intk, realk, mglet_mpi_real
    USE fields_mod, ONLY: get_fieldptr, set_field
    USE err_mod, ONLY: err_abort
    USE comms_mod, ONLY: myid
    USE grids_mod, ONLY: minlevel, maxlevel

    USE mphcore_mod, ONLY: volChk, divChk, vofChk, volTol, divTol, vofTol, hasMph

    IMPLICIT NONE(type, external)
    PRIVATE

    REAL(realk) :: volInit

    PUBLIC :: init_mph_chk, finish_mph_chk, comp_vol, final_vol_chk, itinfo_mph

CONTAINS

    SUBROUTINE init_mph_chk()

        ! Subroutine arguments
        ! None

        ! Local variables
        CHARACTER(len=*), PARAMETER :: descGrdmask = "uncov. cells"
        INTEGER(intk) :: n, igrid, igridf, ipar
        INTEGER(intk) :: kk, jj, ii, kc0, jc0, ic0
        REAL(realk), POINTER, CONTIGUOUS :: grdMask(:,:,:)

        CALL set_field("GRDMASK", description=descGrdmask, &
            dread=.FALSE., required=.TRUE., dwrite=.FALSE., buffers=.TRUE.)

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

    SUBROUTINE final_vol_chk()
    !----------------------------------------------------------------
    !   What it does:
    !
    !----------------------------------------------------------------

        ! Subroutine arguments
        ! None

        ! Local variables
        REAL(realk) :: volFini

        IF ( myid /= 0 ) RETURN

        CALL comp_vol(volFini)
        WRITE(*, '()')
        WRITE(*,'(A,E11.5)') "Final volume is: ", volFini
        WRITE(*,'(A,E11.5)') "Absolute volume error: ", ABS(volInit-volFini)
        WRITE(*,'(A,E11.5)') "Relative volume error: ", ABS(volInit-volFini)/volInit
        WRITE(*, '()')

    END SUBROUTINE final_vol_chk

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

    !================================================================

END MODULE mph_chk_mod