
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

    IMPLICIT NONE(type, external)
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

    SUBROUTINE sel_index(lOrq, io, jo, ko)
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

        io = 0
        jo = 0
        ko = 0

        SELECT CASE ( lOrq )
        CASE ( 1 )
            io = 1
        CASE ( 2 )
            jo = 1
        CASE ( 3 )
            ko = 1
        CASE DEFAULT
            CALL err_abort(vofErr, "invalid direction lOrq.", __FILE__, __LINE__)
        END SELECT

    END SUBROUTINE sel_index

    !================================================================

    SUBROUTINE sel_extent(kk, jj, ii, q, l, dx, dy, dz, ddx, ddy, ddz, dsx, dsy, dsz)
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

        IF ( l < 1 .OR. l > 3 .OR. q < 1 .OR. q > 3 ) THEN
            CALL err_abort(vofErr, "invalid direction l or q.", __FILE__, __LINE__)
        END IF

        dsx = ddx ; dsy = ddy ; dsz = ddz

        IF ( q == l ) THEN
            SELECT CASE ( l )
            CASE ( 1 )
                dsx = dx
            CASE ( 2 )
                dsy = dy
            CASE ( 3 )
                dsz = dz
            END SELECT
        END IF

    END SUBROUTINE sel_extent

    !================================================================

    SUBROUTINE sel_velocity(kk, jj, ii, lOrq, u, v, w, vel)
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

        SELECT CASE ( lOrq )
        CASE ( 1 )
            vel = u
        CASE ( 2 )
            vel = v
        CASE ( 3 )
            vel = w
        CASE DEFAULT
            CALL err_abort(vofErr, "invalid direction lOrq.", __FILE__, __LINE__)
        END SELECT

    END SUBROUTINE sel_velocity

    !================================================================

    SUBROUTINE clip(kk, jj, ii, c)
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
        REAL(realk), INTENT(inout) :: c(kk, jj, ii)

        ! Local variables
        ! None

        c = MAX(MIN(c, 1.0_realk), 0.0_realk)

    END SUBROUTINE clip

    !================================================================

    SUBROUTINE comp_mat_coeff_mph()
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        ! None

        ! Local variables
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
        REAL(realk), POINTER, CONTIGUOUS :: dBa(:, :, :), dLe(:, :, :), dTo(:, :, :)


        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_fieldptr(c, "C", igrid)
            CALL get_fieldptr(dx, "DX", igrid)
            CALL get_fieldptr(dy, "DY", igrid)
            CALL get_fieldptr(dz, "DZ", igrid)
            CALL get_fieldptr(ddx, "DDX", igrid)
            CALL get_fieldptr(ddy, "DDY", igrid)
            CALL get_fieldptr(ddz, "DDZ", igrid)

            CALL get_fieldptr(aw, "GSAW", igrid)
            CALL get_fieldptr(ae, "GSAE", igrid)
            CALL get_fieldptr(as, "GSAS", igrid)
            CALL get_fieldptr(an, "GSAN", igrid)
            CALL get_fieldptr(ab, "GSAB", igrid)
            CALL get_fieldptr(at, "GSAT", igrid)
            CALL get_fieldptr(ap, "GSAP", igrid)

            CALL get_fieldptr(dBa, "dBa", igrid)
            CALL get_fieldptr(dLe, "dLe", igrid)
            CALL get_fieldptr(dTo, "dTo", igrid)

            CALL comp_prop_face(kk, jj, ii, c, dBa, dLe, dTo, rho1, rho2, ddx, ddy, ddz)

            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        ae(k,j,i) = 2.0/((dx(i-1)+dx(i))*dx(i)*dBa(k,j,i))
                        aw(k,j,i) = 2.0/((dx(i-1)+dx(i))*dx(i-1)*dBa(k,j,i-1))
                        an(k,j,i) = 2.0/((dy(j-1)+dy(j))*dy(j)*dLe(k,j,i))
                        as(k,j,i) = 2.0/((dy(j-1)+dy(j))*dy(j-1)*dLe(k,j-1,i))
                        at(k,j,i) = 2.0/((dz(k-1)+dz(k))*dz(k)*dTo(k,j,i))
                        ab(k,j,i) = 2.0/((dz(k-1)+dz(k))*dz(k-1)*dTo(k-1,j,i))
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
        ENDDO

    END SUBROUTINE comp_mat_coeff_mph

END MODULE mph_utils_mod