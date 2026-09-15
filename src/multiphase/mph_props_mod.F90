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

    USE mphcore_mod, ONLY: gmol1, gmol2, rho1, rho2
    USE mph_utils_mod, ONLY: clip
    USE precision_mod, ONLY: intk, realk
    USE grids_mod, ONLY: nmygrids, mygrids
    USE fields_mod, ONLY: get_fieldptr
    USE grids_mod, ONLY: get_mgdims
    USE err_mod, ONLY: err_abort

    IMPLICIT NONE
    PRIVATE

    PUBLIC :: init_mph_props, finish_mph_props, comp_prop, &
        comp_prop_face, comp_prop_face_stag

CONTAINS

    SUBROUTINE init_mph_props()

        ! Subroutine arguments
        ! None

        ! Local variables
        INTEGER(intk) :: n, igrid
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: c(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: d(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: g(:,:,:)

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_fieldptr(c, "C", igrid)
            CALL get_fieldptr(d, "D", igrid)
            CALL get_fieldptr(g, "G", igrid)

            CALL comp_prop(kk, jj, ii, c, d, rho1, rho2, meanFlag='arit')
            CALL comp_prop(kk, jj, ii, c, g, gmol1, gmol2, meanFlag='harm')

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

        IF ( .NOT. PRESENT(meanFlag) ) THEN
            comp_prop_arit(kk, jj, ii, c, propFld, prop1, prop2)
            RETURN
        END IF

        IF ( meanFlag == 'arit' ) THEN
            comp_prop_arit(kk, jj, ii, c, propFld, prop1, prop2)
        ELSE IF ( meanFlag == 'harm' ) THEN
            comp_prop_harm(kk, jj, ii, c, propFld, prop1, prop2)
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
        REAL(realk) :: cClip(kk, jj, ii)

        CALL clip(kk, jj, ii, c, cClip)

        propFld = cClip*( prop1 - prop2 ) + prop2

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
        REAL(realk) :: cClip(kk, jj, ii)
        REAL(realk) :: invProp1, invProp2

        CALL clip(kk, jj, ii, c, cClip)

        invProp1 = 1.0_realk/( prop1 + 1.0E-16_realk )
        invProp2 = 1.0_realk/( prop2 + 1.0E-16_realk )

        propFld = 1.0_realk/( cClip*( invProp1 - invProp2 ) + invProp2 )

        IF ( MAXVAL(propFld) > MAX(prop1, prop2) .OR. &
             MINVAL(propFld) < MIN(prop1, prop2) ) THEN
            CALL err_abort(propsErr, "property out of bounds.", __FILE__, __LINE__)
        END IF

    END SUBROUTINE comp_prop_harm

    !================================================================

    SUBROUTINE comp_prop_face(kk, jj, ii, propFld, pe, pn, pt, prop1, prop2, rdx, rdy, rdz, meanFlag)
    !----------------------------------------------------------------
    !   What it does:
    !   Compute property face values depending on a property field
    !   using one of two means.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: propFld(kk, jj, ii)
        REAL(realk), INTENT(out) :: pe(kk, jj, ii), pn(kk, jj, ii), pt(kk, jj, ii)
        REAL(realk), INTENT(in) :: prop1, prop2
        REAL(realk), INTENT(in) :: rdx(ii), rdy(jj), rdz(kk)
        CHARACTER(len=4), INTENT(in), OPTIONAL :: meanFlag

        ! Local variables
        ! None

        IF ( .NOT. PRESENT(meanFlag) ) THEN
            CALL comp_prop_arit_face(kk, jj, ii, propFld, pe, pn, pt, prop1, prop2, rdx, rdy, rdz)
            RETURN
        END IF

        IF ( meanFlag == 'arit' ) THEN
            CALL comp_prop_arit_face(kk, jj, ii, propFld, pe, pn, pt, prop1, prop2, rdx, rdy, rdz)
        ELSE IF ( meanFlag == 'harm' ) THEN
            CALL comp_prop_harm_face(kk, jj, ii, propFld, pe, pn, pt, prop1, prop2, rdx, rdy, rdz)
        END IF

    END SUBROUTINE comp_prop_face

    !================================================================

    SUBROUTINE comp_prop_arit_face(kk, jj, ii, propFld, pe, pn, pt, prop1, prop2, rdx, rdy, rdz)
    !----------------------------------------------------------------
    !   What it does:
    !   Compute property face values depending on a property field
    !   using the arithmetic mean.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: propFld(kk, jj, ii)
        REAL(realk), INTENT(out) :: pe(kk, jj, ii), pn(kk, jj, ii), pt(kk, jj, ii)
        REAL(realk), INTENT(in) :: prop1, prop2
        REAL(realk), INTENT(in) :: rdx(ii), rdy(jj), rdz(kk)

        ! Local variables
        INTEGER(intk) :: k, j, i

        DO i = 2, ii-2
            DO j = 2, jj-2
                DO k = 2, kk-2
                    pe(k,j,i) = (propFld(k,j,i)*rdx(i) + propFld(k,j,i+1)*rdx(i+1))/(rdx(i) + rdx(i+1))
                    pn(k,j,i) = (propFld(k,j,i)*rdy(j) + propFld(k,j+1,i)*rdy(j+1))/(rdy(j) + rdy(j+1))
                    pt(k,j,i) = (propFld(k,j,i)*rdz(k) + propFld(k+1,j,i)*rdz(k+1))/(rdz(k) + rdz(k+1))
                ENDDO
            ENDDO
        ENDDO

    ENDSUBROUTINE comp_prop_arit_face

    !================================================================

    SUBROUTINE comp_prop_harm_face(kk, jj, ii, propFld, pe, pn, pt, prop1, prop2, rdx, rdy, rdz)
    !----------------------------------------------------------------
    !   What it does:
    !   Compute property face values depending on a property field
    !   using the harmonic mean.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: propFld(kk, jj, ii)
        REAL(realk), INTENT(out) :: pe(kk, jj, ii), pn(kk, jj, ii), pt(kk, jj, ii)
        REAL(realk), INTENT(in) :: prop1, prop2
        REAL(realk), INTENT(in) :: rdx(ii), rdy(jj), rdz(kk)

        ! Local variables
        INTEGER(intk) :: k, j, i
        REAL(realk) :: invPropFld(kk, jj, ii)

        invPropFld = 1.0_realk/( propFld + 1.0E-16_realk )

        DO i = 2, ii-2
            DO j = 2, jj-2
                DO k = 2, kk-2
                    pe(k,j,i) = (rdx(i) + rdx(i+1))/(invPropFld(k,j,i)*rdx(i) + invPropFld(k,j,i+1)*rdx(i+1))
                    pn(k,j,i) = (rdy(j) + rdy(j+1))/(invPropFld(k,j,i)*rdy(j) + invPropFld(k,j+1,i)*rdy(j+1))
                    pt(k,j,i) = (rdz(k) + rdz(k+1))/(invPropFld(k,j,i)*rdz(k) + invPropFld(k+1,j,i)*rdz(k+1))
                ENDDO
            ENDDO
        ENDDO

    ENDSUBROUTINE comp_prop_harm_face

    !================================================================

    SUBROUTINE comp_prop_face_stag(kk, jj, ii, c, pxy, pxz, pyz, prop1, prop2, ddx, ddy, ddz, meanFlag)
    !----------------------------------------------------------------
    !   What it does:
    !   Compute the "matrix" of a property on the staggered cells
    !   faces. Each cell has 6 faces, hence there are 18 values to 
    !   compute. Since (k,j,i)+ = (k+kq,j+jq,i+iq)- this reduces to
    !   9 values. 
    !
    !   Of these 9, 6 are the same:
    !   xStagN = yStagE, xStagW = yStagS,
    !   xStagT = zStagE, xStagW = zStagB,
    !   yStagT = zStagN, yStagS = zStagB.
    !
    !   Additionally, the diagonal values are trivial, since they
    !   are centered on the pressure grid.
    !
    !   -> 3 different values to compute!
    !
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

        IF ( .NOT. PRESENT(meanFlag) ) THEN
            comp_prop_arit_face_stag(kk, jj, ii, c, pxy, pxz, pyz, prop1, prop2, ddx, ddy, ddz)
            RETURN
        END IF

        IF ( meanFlag == 'arit' ) THEN
            comp_prop_arit_face_stag(kk, jj, ii, c, pxy, pxz, pyz, prop1, prop2, ddx, ddy, ddz)
        ELSE IF ( meanFlag == 'harm' ) THEN
            comp_prop_harm_face_stag(kk, jj, ii, c, pxy, pxz, pyz, prop1, prop2, ddx, ddy, ddz)
        END IF

    ENDSUBROUTINE comp_prop_face_stag

    !================================================================

    SUBROUTINE comp_prop_arit_face_stag(kk, jj, ii, c, pxy, pxz, pyz, prop1, prop2, ddx, ddy, ddz)
    !----------------------------------------------------------------
    !   What it does:
    !   Compute property face values on the staggered grid depending 
    !   on c using the arithmetic mean.
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
                    phi = MIN(MAX(2.0_realk*(c(k,j,i)*w(1) + c(k,j,i+1)*w(2) + &
                                                    c(k,j+1,i)*w(3) + c(k,j+1,i+1)*w(4))/SUM(w) - &
                                                    0.5_realk, 0.0_realk), 1.0_realk)
                    pxy(k,j,i) = phi*( prop1 - prop2 ) + prop2

                    w(1) = ddx(i)*ddz(k)
                    w(2) = ddx(i+1)*ddz(k)
                    w(3) = ddx(i)*ddz(k+1)
                    w(4) = ddx(i+1)*ddz(k+1)
                    phi = MIN(MAX(2.0_realk*(c(k,j,i)*w(1) + c(k,j,i+1)*w(2) + &
                                                    c(k+1,j,i)*w(3) + c(k+1,j,i+1)*w(4))/SUM(w) - &
                                                    0.5_realk, 0.0_realk), 1.0_realk)
                    pxz(k,j,i) = phi*( prop1 - prop2 ) + prop2

                    w(1) = ddy(j)*ddz(k)
                    w(2) = ddy(j+1)*ddz(k)
                    w(3) = ddy(j)*ddz(k+1)
                    w(4) = ddy(j+1)*ddz(k+1)
                    phi = MIN(MAX(2.0_realk*(c(k,j,i)*w(1) + c(k,j+1,i)*w(2) + &
                                                    c(k+1,j,i)*w(3) + c(k+1,j+1,i)*w(4))/SUM(w) - &
                                                    0.5_realk, 0.0_realk), 1.0_realk)
                    pyz(k,j,i) = phi*( prop1 - prop2 ) + prop2
                ENDDO
            ENDDO
        ENDDO

    ENDSUBROUTINE comp_prop_arit_face_stag

    !================================================================

    SUBROUTINE comp_prop_harm_face_stag(kk, jj, ii, c, pxy, pxz, pyz, prop1, prop2, ddx, ddy, ddz)
    !----------------------------------------------------------------
    !   What it does:
    !   Compute property face values on the staggered grid depending 
    !   on c using the harmonic mean.
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

        invProp1 = 1.0_realk/( prop1 + 1.0E-16_realk )
        invProp2 = 1.0_realk/( prop2 + 1.0E-16_realk )

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

    ENDSUBROUTINE comp_prop_harm_face_stag

END MODULE mph_props_mod