!====================================================================
!  Module: mph_props_mod
!
!   Responsibilities:
!   - Fill property fields initially
!   - Provide routines related to material properties
!
!   Author:      Valentin Ruhs
!   e-Mail:      valentin.ruhs@gmx.de
!   Created:     2026-02
!   Last update: 2026-09
!
!====================================================================

MODULE mph_props_mod

    USE mphcore_mod, ONLY: gmol1, gmol2, rho1, rho2, propsErr
    USE mph_utils_mod, ONLY: clp, int2char
    USE precision_mod, ONLY: intk, realk
    USE grids_mod, ONLY: nmygrids, mygrids, get_mgdims
    USE fields_mod, ONLY: get_fieldptr, set_field
    USE err_mod, ONLY: err_abort

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: init_mph_props, finish_mph_props, comp_prop, &
        comp_prop_face, comp_prop_face_stag, comp_d_stag, comp_props

CONTAINS

    SUBROUTINE init_mph_props()

        ! Subroutine arguments
        ! None

        ! Local variables
        CHARACTER(len=*), PARAMETER :: descDBa = "density back"
        CHARACTER(len=*), PARAMETER :: descDLe = "density left"
        CHARACTER(len=*), PARAMETER :: descDTo = "density top"
        CHARACTER(len=*), PARAMETER :: descGUv = "visc. coupling u/v"
        CHARACTER(len=*), PARAMETER :: descGUw = "visc. coupling u/w"
        CHARACTER(len=*), PARAMETER :: descGVw = "visc. coupling w/v"
        INTEGER(intk) :: n, igrid
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: c(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: d(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: g(:,:,:)

        CALL set_field("DBA", description=descDBa)
        CALL set_field("DLE", description=descDLe)
        CALL set_field("DTO", description=descDTo)
        CALL set_field("GUV", description=descGUv)
        CALL set_field("GUW", description=descGUw)
        CALL set_field("GVW", description=descGVw)

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_fieldptr(c, "C", igrid)
            CALL get_fieldptr(d, "D", igrid)
            CALL get_fieldptr(g, "G", igrid)

            CALL comp_prop(kk, jj, ii, c, d, rho1, rho2, meanFlag="arit")
            CALL comp_prop(kk, jj, ii, c, g, gmol1, gmol2, meanFlag="harm")

        END DO

    END SUBROUTINE init_mph_props

    !================================================================

    SUBROUTINE finish_mph_props()

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        CONTINUE

    END SUBROUTINE finish_mph_props

    !================================================================

    SUBROUTINE comp_prop(kk, jj, ii, c, propFld, prop1, prop2, meanFlag)
    !----------------------------------------------------------------
    !   What it does:
    !   Compute property field values depending on c using one of two
    !   means.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: c(kk, jj, ii)
        REAL(realk), INTENT(out) :: propFld(kk, jj, ii)
        REAL(realk), INTENT(in) :: prop1, prop2
        CHARACTER(len=4), INTENT(in), OPTIONAL :: meanFlag

        ! Local variables
        ! None

        IF ( PRESENT(meanFlag) ) THEN
            SELECT CASE ( meanFlag )
            CASE ( "arit" )
                CALL comp_prop_arit(kk, jj, ii, c, propFld, prop1, prop2)
            CASE ( "harm" )
                CALL comp_prop_harm(kk, jj, ii, c, propFld, prop1, prop2)
            CASE DEFAULT
                CALL err_abort(propsErr, "unknown meanFlag.", __FILE__, __LINE__)
            END SELECT
        ELSE
            CALL comp_prop_arit(kk, jj, ii, c, propFld, prop1, prop2)
        END IF

    END SUBROUTINE comp_prop

    !================================================================

    SUBROUTINE comp_prop_arit(kk, jj, ii, c, propFld, prop1, prop2)
    !----------------------------------------------------------------
    !   What it does:
    !   Compute property field values depending on c using the 
    !   arithmetic mean.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: c(kk, jj, ii)
        REAL(realk), INTENT(out) :: propFld(kk, jj, ii)
        REAL(realk), INTENT(in) :: prop1, prop2

        ! Local variables
        ! None

        propFld = clp(c)*( prop1 - prop2 ) + prop2

        IF ( MAXVAL(propFld) > MAX(prop1, prop2) .OR. &
             MINVAL(propFld) < MIN(prop1, prop2) ) THEN
            CALL err_abort(propsErr, "property out of bounds.", __FILE__, __LINE__)
        END IF

    END SUBROUTINE comp_prop_arit

    !================================================================

    SUBROUTINE comp_prop_harm(kk, jj, ii, c, propFld, prop1, prop2)
    !----------------------------------------------------------------
    !   What it does:
    !   Compute property field values depending on c using the 
    !   harmonic mean.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: c(kk, jj, ii)
        REAL(realk), INTENT(out) :: propFld(kk, jj, ii)
        REAL(realk), INTENT(in) :: prop1, prop2

        ! Local variables
        REAL(realk) :: invProp1, invProp2

        invProp1 = 1.0_realk/prop1
        invProp2 = 1.0_realk/prop2

        propFld = 1.0_realk/( clp(c)*( invProp1 - invProp2 ) + invProp2 )

        IF ( MAXVAL(propFld) > MAX(prop1, prop2) .OR. &
             MINVAL(propFld) < MIN(prop1, prop2) ) THEN
            CALL err_abort(propsErr, "property out of bounds.", __FILE__, __LINE__)
        END IF

    END SUBROUTINE comp_prop_harm

    !================================================================

    SUBROUTINE comp_prop_face(kk, jj, ii, c, pBa, pLe, pTo, prop1, prop2, ddx, ddy, ddz, meanFlag)
    !----------------------------------------------------------------
    !   What it does:
    !   Compute property face values depending on c using one of two 
    !   means.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: c(kk, jj, ii)
        REAL(realk), INTENT(out) :: pBa(kk, jj, ii), pLe(kk, jj, ii), pTo(kk, jj, ii)
        REAL(realk), INTENT(in) :: prop1, prop2
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        CHARACTER(len=4), INTENT(in), OPTIONAL :: meanFlag

        ! Local variables
        ! None

        IF ( PRESENT(meanFlag) ) THEN
            SELECT CASE ( meanFlag )
            CASE ( "arit" )
                CALL comp_prop_arit_face(kk, jj, ii, c, pBa, pLe, pTo, prop1, prop2, ddx, ddy, ddz)
            CASE ( "harm" )
                CALL comp_prop_harm_face(kk, jj, ii, c, pBa, pLe, pTo, prop1, prop2, ddx, ddy, ddz)
            CASE DEFAULT
                CALL err_abort(propsErr, "unknown meanFlag.", __FILE__, __LINE__)
            END SELECT
        ELSE
            CALL comp_prop_arit_face(kk, jj, ii, c, pBa, pLe, pTo, prop1, prop2, ddx, ddy, ddz)
        END IF

    END SUBROUTINE comp_prop_face

    !================================================================

    SUBROUTINE comp_prop_arit_face(kk, jj, ii, c, pBa, pLe, pTo, prop1, prop2, ddx, ddy, ddz)
    !----------------------------------------------------------------
    !   What it does:
    !   Compute property face values depending on c using the 
    !   arithmetic mean.
    !
    !   Source:
    !   G. Tryggvason, R. Scardovelli, und S. Zaleski, &
    !   Direct Numerical Simulations of Gas–Liquid Multiphase &
    !   Flows, 1. Aufl. Cambridge University Press, 2011. &
    !   doi: 10.1017/CBO9780511975264.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: c(kk, jj, ii)
        REAL(realk), INTENT(out) :: pBa(kk, jj, ii), pLe(kk, jj, ii), pTo(kk, jj, ii)
        REAL(realk), INTENT(in) :: prop1, prop2
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)

        ! Local variables
        INTEGER(intk) :: k, j, i
        REAL(realk) :: w(2), phi

        DO i = 2, ii-2
            DO j = 2, jj-2
                DO k = 2, kk-2
                    w(1) = ddx(i)
                    w(2) = ddx(i+1)
                    phi = MIN(MAX((c(k,j,i)*w(1) + c(k,j,i+1)*w(2))/SUM(w), 0.0_realk), 1.0_realk)
                    pBa(k,j,i) = phi*( prop1 - prop2 ) + prop2

                    w(1) = ddy(j)
                    w(2) = ddy(j+1)
                    phi = MIN(MAX((c(k,j,i)*w(1) + c(k,j+1,i)*w(2))/SUM(w), 0.0_realk), 1.0_realk)
                    pLe(k,j,i) = phi*( prop1 - prop2 ) + prop2

                    w(1) = ddz(k)
                    w(2) = ddz(k+1)
                    phi = MIN(MAX((c(k,j,i)*w(1) + c(k+1,j,i)*w(2))/SUM(w), 0.0_realk), 1.0_realk)
                    pTo(k,j,i) = phi*( prop1 - prop2 ) + prop2
                ENDDO
            ENDDO
        ENDDO

    END SUBROUTINE comp_prop_arit_face

    !================================================================

    SUBROUTINE comp_prop_harm_face(kk, jj, ii, c, pBa, pLe, pTo, prop1, prop2, ddx, ddy, ddz)
    !----------------------------------------------------------------
    !   What it does:
    !   Compute property face values depending on c using the 
    !   harmonic mean. Care, only valid for 
    !   equidistant grids.
    !
    !   Source:
    !   G. Tryggvason, R. Scardovelli, und S. Zaleski, &
    !   Direct Numerical Simulations of Gas–Liquid Multiphase &
    !   Flows, 1. Aufl. Cambridge University Press, 2011. &
    !   doi: 10.1017/CBO9780511975264.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: c(kk, jj, ii)
        REAL(realk), INTENT(out) :: pBa(kk, jj, ii), pLe(kk, jj, ii), pTo(kk, jj, ii)
        REAL(realk), INTENT(in) :: prop1, prop2
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)

        ! Local variables
        INTEGER(intk) :: k, j, i
        REAL(realk) :: w(2), phi
        REAL(realk) :: invProp1, invProp2

        invProp1 = 1.0_realk/prop1
        invProp2 = 1.0_realk/prop2

        DO i = 2, ii-2
            DO j = 2, jj-2
                DO k = 2, kk-2
                    w(1) = ddx(i)
                    w(2) = ddx(i+1)
                    phi = MIN(MAX(2.0_realk*(c(k,j,i)*w(1) + c(k,j,i+1)*w(2))/SUM(w) - &
                                  0.5_realk, 0.0_realk), 1.0_realk)
                    pBa(k,j,i) = 1.0_realk/( phi*( invProp1 - invProp2 ) + invProp2 )

                    w(1) = ddy(j)
                    w(2) = ddy(j+1)
                    phi = MIN(MAX(2.0_realk*(c(k,j,i)*w(1) + c(k,j+1,i)*w(2))/SUM(w) - &
                                  0.5_realk, 0.0_realk), 1.0_realk)
                    pLe(k,j,i) = 1.0_realk/( phi*( invProp1 - invProp2 ) + invProp2 )

                    w(1) = ddz(k)
                    w(2) = ddz(k+1)
                    phi = MIN(MAX(2.0_realk*(c(k,j,i)*w(1) + c(k+1,j,i)*w(2))/SUM(w) - &
                                  0.5_realk, 0.0_realk), 1.0_realk)
                    pTo(k,j,i) = 1.0_realk/( phi*( invProp1 - invProp2 ) + invProp2 )
                ENDDO
            ENDDO
        ENDDO

    END SUBROUTINE comp_prop_harm_face

    !================================================================

    SUBROUTINE comp_prop_face_stag(kk, jj, ii, c, pxy, pxz, pyz, prop1, prop2, ddx, ddy, ddz, meanFlag)
    !----------------------------------------------------------------
    !   What it does:
    !   Compute property face values on the staggered grid depending 
    !   on c using one of two means.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: c(kk, jj, ii)
        REAL(realk), INTENT(out) :: pxy(kk, jj, ii), pxz(kk, jj, ii), pyz(kk, jj, ii)
        REAL(realk), INTENT(in) :: prop1, prop2
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        CHARACTER(len=4), INTENT(in), OPTIONAL :: meanFlag

        ! Local variables
        ! None

        IF ( PRESENT(meanFlag) ) THEN
            SELECT CASE ( meanFlag )
            CASE ( "arit" )
                CALL comp_prop_arit_face_stag(kk, jj, ii, c, pxy, pxz, pyz, prop1, prop2, ddx, ddy, ddz)
            CASE ( "harm" )
                CALL comp_prop_harm_face_stag(kk, jj, ii, c, pxy, pxz, pyz, prop1, prop2, ddx, ddy, ddz)
            CASE DEFAULT
                CALL err_abort(propsErr, "unknown meanFlag.", __FILE__, __LINE__)
            END SELECT
        ELSE
            CALL comp_prop_arit_face_stag(kk, jj, ii, c, pxy, pxz, pyz, prop1, prop2, ddx, ddy, ddz)
        END IF

    END SUBROUTINE comp_prop_face_stag

    !================================================================

    SUBROUTINE comp_prop_arit_face_stag(kk, jj, ii, c, pxy, pxz, pyz, prop1, prop2, ddx, ddy, ddz)
    !----------------------------------------------------------------
    !   What it does:
    !   Compute property face values on the staggered grid depending 
    !   on c using the arithmetic mean.
    !
    !   Source:
    !   G. Tryggvason, R. Scardovelli, und S. Zaleski, &
    !   Direct Numerical Simulations of Gas–Liquid Multiphase &
    !   Flows, 1. Aufl. Cambridge University Press, 2011. &
    !   doi: 10.1017/CBO9780511975264.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: c(kk, jj, ii)
        REAL(realk), INTENT(out) :: pxy(kk, jj, ii), pxz(kk, jj, ii), pyz(kk, jj, ii)
        REAL(realk), INTENT(in) :: prop1, prop2
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)

        ! Local variables
        INTEGER(intk) :: k, j, i
        REAL(realk) :: w(4), phi

        DO i = 2, ii-2
            DO j = 2, jj-2
                DO k = 2, kk-2
                    w(1) = ddx(i)*ddy(j)
                    w(2) = ddx(i+1)*ddy(j)
                    w(3) = ddx(i)*ddy(j+1)
                    w(4) = ddx(i+1)*ddy(j+1)
                    phi = MIN(MAX((c(k,j,i)*w(1) + c(k,j,i+1)*w(2) + &
                                   c(k,j+1,i)*w(3) + c(k,j+1,i+1)*w(4))/SUM(w), 0.0_realk), 1.0_realk)
                    pxy(k,j,i) = phi*( prop1 - prop2 ) + prop2

                    w(1) = ddx(i)*ddz(k)
                    w(2) = ddx(i+1)*ddz(k)
                    w(3) = ddx(i)*ddz(k+1)
                    w(4) = ddx(i+1)*ddz(k+1)
                    phi = MIN(MAX((c(k,j,i)*w(1) + c(k,j,i+1)*w(2) + &
                                   c(k+1,j,i)*w(3) + c(k+1,j,i+1)*w(4))/SUM(w), 0.0_realk), 1.0_realk)
                    pxz(k,j,i) = phi*( prop1 - prop2 ) + prop2

                    w(1) = ddy(j)*ddz(k)
                    w(2) = ddy(j+1)*ddz(k)
                    w(3) = ddy(j)*ddz(k+1)
                    w(4) = ddy(j+1)*ddz(k+1)
                    phi = MIN(MAX((c(k,j,i)*w(1) + c(k,j+1,i)*w(2) + &
                                   c(k+1,j,i)*w(3) + c(k+1,j+1,i)*w(4))/SUM(w), 0.0_realk), 1.0_realk)
                    pyz(k,j,i) = phi*( prop1 - prop2 ) + prop2
                ENDDO
            ENDDO
        ENDDO

    END SUBROUTINE comp_prop_arit_face_stag

    !================================================================

    SUBROUTINE comp_prop_harm_face_stag(kk, jj, ii, c, pxy, pxz, pyz, prop1, prop2, ddx, ddy, ddz)
    !----------------------------------------------------------------
    !   What it does:
    !   Compute property face values on the staggered grid depending 
    !   on c using the harmonic mean. Care, only valid for 
    !   equidistant grids.
    !
    !   Source:
    !   G. Tryggvason, R. Scardovelli, und S. Zaleski, &
    !   Direct Numerical Simulations of Gas–Liquid Multiphase &
    !   Flows, 1. Aufl. Cambridge University Press, 2011. &
    !   doi: 10.1017/CBO9780511975264.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: c(kk, jj, ii)
        REAL(realk), INTENT(out) :: pxy(kk, jj, ii), pxz(kk, jj, ii), pyz(kk, jj, ii)
        REAL(realk), INTENT(in) :: prop1, prop2
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)

        ! Local variables
        INTEGER(intk) :: k, j, i
        REAL(realk) :: w(4), phi
        REAL(realk) :: invProp1, invProp2

        invProp1 = 1.0_realk/prop1
        invProp2 = 1.0_realk/prop2

        DO i = 2, ii-2
            DO j = 2, jj-2
                DO k = 2, kk-2
                    w(1) = ddx(i)*ddy(j)
                    w(2) = ddx(i+1)*ddy(j)
                    w(3) = ddx(i)*ddy(j+1)
                    w(4) = ddx(i+1)*ddy(j+1)
                    phi = MIN(MAX(2.0_realk*(c(k,j,i)*w(1) + c(k,j,i+1)*w(2) + &
                                             c(k,j+1,i)*w(3) + c(k,j+1,i+1)*w(4))/SUM(w) - &
                                             0.5_realk, 0.0_realk), 1.0_realk)
                    pxy(k,j,i) = 1.0_realk/( phi*( invProp1 - invProp2 ) + invProp2 )

                    w(1) = ddx(i)*ddz(k)
                    w(2) = ddx(i+1)*ddz(k)
                    w(3) = ddx(i)*ddz(k+1)
                    w(4) = ddx(i+1)*ddz(k+1)
                    phi = MIN(MAX(2.0_realk*(c(k,j,i)*w(1) + c(k,j,i+1)*w(2) + &
                                             c(k+1,j,i)*w(3) + c(k+1,j,i+1)*w(4))/SUM(w) - &
                                             0.5_realk, 0.0_realk), 1.0_realk)
                    pxz(k,j,i) = 1.0_realk/( phi*( invProp1 - invProp2 ) + invProp2 )

                    w(1) = ddy(j)*ddz(k)
                    w(2) = ddy(j+1)*ddz(k)
                    w(3) = ddy(j)*ddz(k+1)
                    w(4) = ddy(j+1)*ddz(k+1)
                    phi = MIN(MAX(2.0_realk*(c(k,j,i)*w(1) + c(k,j+1,i)*w(2) + &
                                             c(k+1,j,i)*w(3) + c(k+1,j+1,i)*w(4))/SUM(w) - &
                                             0.5_realk, 0.0_realk), 1.0_realk)
                    pyz(k,j,i) = 1.0_realk/( phi*( invProp1 - invProp2 ) + invProp2 )
                ENDDO
            ENDDO
        ENDDO

    END SUBROUTINE comp_prop_harm_face_stag

    !================================================================

    SUBROUTINE comp_d_stag(q)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: q

        ! Local variables
        CHARACTER(len=3) :: cFldName, dFldName
        INTEGER(intk) :: n, igrid
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: cSq(:,:,:), dSq(:,:,:)

        cFldName = "CS"//int2char(q)
        dFldName = "DS"//int2char(q)

        DO n = 1, nmygrids
            igrid = mygrids(n)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_fieldptr(cSq, cFldName, igrid)
            CALL get_fieldptr(dSq, dFldName, igrid)

            CALL comp_prop(kk, jj, ii, cSq, dSq, rho1, rho2, meanFlag="arit")
        END DO

    END SUBROUTINE comp_d_stag

    !================================================================

    SUBROUTINE comp_props()
    !----------------------------------------------------------------
    !   What it does:
    !   Compute all nescessary properies 
    !----------------------------------------------------------------

        ! Subroutine arguments
        ! None

        ! Local variables
        INTEGER(intk) :: n, igrid
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: cp(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: g(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: gUv(:,:,:), gUw(:,:,:), gVw(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: dBa(:,:,:), dLe(:,:,:), dTo(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)

        DO n = 1, nmygrids
            igrid = mygrids(n)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_fieldptr(cp, "CP", igrid)
            CALL get_fieldptr(g, "G", igrid)
            CALL get_fieldptr(gUv, "GUV", igrid)
            CALL get_fieldptr(gUw, "GUW", igrid)
            CALL get_fieldptr(gVw, "GVW", igrid)
            CALL get_fieldptr(dBa, "DBA", igrid)
            CALL get_fieldptr(dLe, "DLE", igrid)
            CALL get_fieldptr(dTo, "DTO", igrid)
            CALL get_fieldptr(ddx, "DDX", igrid)
            CALL get_fieldptr(ddy, "DDY", igrid)
            CALL get_fieldptr(ddz, "DDZ", igrid)

            CALL comp_prop(kk, jj, ii, cp, g, gmol1, gmol2, meanFlag="harm")
            CALL comp_prop_face_stag(kk, jj, ii, cp, gUv, gUw, gVw, &
                gmol1, gmol2, ddx, ddy, ddz, meanFlag="harm")
            CALL comp_prop_face(kk, jj, ii, cp, dBa, dLe, dTo, &
                rho1, rho2, ddx, ddy, ddz, meanFlag="arit")
        END DO

    END SUBROUTINE comp_props

END MODULE mph_props_mod