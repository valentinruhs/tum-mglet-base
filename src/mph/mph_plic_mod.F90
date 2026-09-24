!====================================================================
!  Module: mph_plic_mod
!
!   Responsibilities:
!   - Provides routines related to piecewise linear interface 
!     calculation
!
!   Author:      Valentin Ruhs
!   e-Mail:      valentin.ruhs@gmx.de
!   Created:     2026-02
!   Last update: 2026-09
!
!   Source:
!   R. Scardovelli und S. Zaleski, „Analytical Relations 
!   Connecting Linear Interfaces and Volume Fractions in 
!   Rectangular Grids“, Journal of Computational Physics, 
!   Bd. 164, Nr. 1, S. 228–237, Okt. 2000, 
!   doi: 10.1006/jcph.2000.6567.
!
!====================================================================

MODULE mph_plic_mod

    USE precision_mod, ONLY: intk, realk
    USE mph_utils_mod, ONLY: sel_ind, int2char
    USE grids_mod, ONLY: nmygrids, mygrids, get_mgdims
    USE mphcore_mod, ONLY: vofTol
    USE fields_mod, ONLY: set_field, get_fieldptr

    IMPLICIT NONE(type, external)
    PRIVATE 

    PUBLIC :: init_mph_plic, finish_mph_plic, &
        comp_ifc, comp_c_stg, comp_isIfc_stg, comp_c_loc

CONTAINS

    SUBROUTINE init_mph_plic()

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        CHARACTER(len=*), PARAMETER :: descNorm = "ifc. norm. vec."
        CHARACTER(len=*), PARAMETER :: descAlpha = "ifc. plane const."
        CHARACTER(len=*), PARAMETER :: descIsIfc = "cell with ifc."
        CHARACTER(len=*), PARAMETER :: descIsIfcVic = "cell in vic. of ifc."

        CALL set_field("NORMX", description=descNorm, &
            dread=.FALSE., required=.TRUE., dwrite=.TRUE., buffers=.TRUE.)
        CALL set_field("NORMY", description=descNorm, &
            dread=.FALSE., required=.TRUE., dwrite=.TRUE., buffers=.TRUE.)
        CALL set_field("NORMZ", description=descNorm, &
            dread=.FALSE., required=.TRUE., dwrite=.TRUE., buffers=.TRUE.)
        CALL set_field("ALPHA", description=descAlpha, &
            dread=.FALSE., required=.TRUE., dwrite=.TRUE., buffers=.TRUE.)
        CALL set_field("ISIFC", description=descIsIfc, &
            dread=.FALSE., required=.TRUE., dwrite=.TRUE., buffers=.TRUE.)
        CALL set_field("ISIFCS1", description=descIsIfc, &
            dread=.FALSE., required=.TRUE., dwrite=.TRUE., buffers=.TRUE.)
        CALL set_field("ISIFCS2", description=descIsIfc, &
            dread=.FALSE., required=.TRUE., dwrite=.TRUE., buffers=.TRUE.)
        CALL set_field("ISIFCS3", description=descIsIfc, &
            dread=.FALSE., required=.TRUE., dwrite=.TRUE., buffers=.TRUE.)
        CALL set_field("ISIFCVICS1", description=descIsIfcVic, &
            dread=.FALSE., required=.TRUE., dwrite=.TRUE., buffers=.TRUE.)
        CALL set_field("ISIFCVICS2", description=descIsIfcVic, &
            dread=.FALSE., required=.TRUE., dwrite=.TRUE., buffers=.TRUE.)
        CALL set_field("ISIFCVICS3", description=descIsIfcVic, &
            dread=.FALSE., required=.TRUE., dwrite=.TRUE., buffers=.TRUE.)

    END SUBROUTINE init_mph_plic

    !================================================================

    SUBROUTINE finish_mph_plic()

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        CONTINUE

    END SUBROUTINE finish_mph_plic

    !================================================================

    SUBROUTINE comp_isIfc_stg(q)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: q

        ! Local variables
        CHARACTER(len=10) :: cFldName, isIfcFldName, isIfcVicFldName
        INTEGER(intk) :: n, igrid
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: cSq(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: isIfcSq(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: isIfcVicSq(:,:,:)

        cFldName = "CS"//int2char(q)
        isIfcFldName = "ISIFCS"//int2char(q)
        isIfcVicFldName = "ISIFCVICS"//int2char(q)

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_fieldptr(cSq, cFldName, igrid)
            CALL get_fieldptr(isIfcSq, isIfcFldName, igrid)
            CALL get_fieldptr(isIfcVicSq, isIfcVicFldName, igrid)

            CALL comp_isIfc_grd(kk, jj, ii, cSq, isIfcSq, isIfcVicSq)
        END DO

    END SUBROUTINE comp_isIfc_stg

    !================================================================

    SUBROUTINE comp_isIfc_grd(kk, jj, ii, c, isIfc, isIfcVic)
    !----------------------------------------------------------------
    !   What it does:
    !   Track cells containing an interface or in vicinity of an
    !   interface.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: c(kk, jj, ii)
        REAL(realk), INTENT(out) :: isIfc(kk, jj, ii)
        REAL(realk), INTENT(out), OPTIONAL :: isIfcVic(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i, vic

        isIfc = -1.0_realk
        DO i = 1, ii
            DO j = 1, jj
                DO k = 1, kk
                    IF ( c(k,j,i) > vofTol .AND. c(k,j,i) < 1.0_realk - vofTol ) THEN
                        isIfc(k,j,i) = 1.0_realk
                    ENDIF
                ENDDO
            ENDDO
        ENDDO

        IF ( .NOT. PRESENT(isIfcVic) ) RETURN

        vic = 2
        isIfcVic = -1.0_realk
        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    IF ( isIfc(k,j,i) > 0.0_realk ) THEN
                        isIfcVic(k-vic:k+vic,j-vic:j+vic,i-vic:i+vic) = 1.0_realk
                    ENDIF
                ENDDO
            ENDDO
        ENDDO

    END SUBROUTINE comp_isIfc_grd

    !================================================================

    SUBROUTINE comp_norm_vec(normx, normy, normz, kk, jj, ii, c, dx, dy, dz)
    !----------------------------------------------------------------
    !   What it does:
    !   Computes the normal vector components normx, normy and normz
    !   of the gradient of the volume fraction function c. The 
    !   components are normalized by the length to get the unit
    !   normal components. The gradient in each cell is calculated
    !   by taking into account its eight surrounding cells weighted
    !   with the three-dimensional stcl operator:
    !            / 1  2  1 \
    !   stcl =   | 2  4  2 |
    !            \ 1  2  1 /
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(out) :: normx(kk, jj, ii), normy(kk, jj, ii), normz(kk, jj, ii)
        REAL(realk), INTENT(in) :: c(kk, jj, ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)

        ! Local variables
        INTEGER(intk) :: k, j, i, d1, d2
        INTEGER(intk), PARAMETER :: stcl(-1:1, -1:1) = &
            RESHAPE([1, 2, 1, 2, 4, 2, 1, 2, 1], [3,3])
        REAL(realk) :: length
        REAL(realk) :: sumx, sumy, sumz

        DO i = 2, ii-1
            DO j = 2, jj-1
                DO k = 2, kk-1

                    sumx = 0.0_realk
                    sumy = 0.0_realk
                    sumz = 0.0_realk

                    DO d2 = -1, 1
                        DO d1 = -1, 1
                            sumx = sumx + stcl(d1, d2) * &
                                ( c(k+d1, j+d2, i+1) - c(k+d1, j+d2, i-1) )
                            sumy = sumy + stcl(d1, d2) * &
                                ( c(k+d1, j+1, i+d2) - c(k+d1, j-1, i+d2) )
                            sumz = sumz + stcl(d1, d2) * &
                                ( c(k+1, j+d1, i+d2) - c(k-1, j+d1, i+d2) )
                        ENDDO
                    ENDDO

                    normx(k,j,i) = sumx / ( SUM(stcl) * (dx(i-1)+dx(i)) )
                    normy(k,j,i) = sumy / ( SUM(stcl) * (dy(j-1)+dy(j)) )
                    normz(k,j,i) = sumz / ( SUM(stcl) * (dz(k-1)+dz(k)) )

                    ! Calculate normal vector length
                    length = SQRT( normx(k,j,i)**2 + &
                                   normy(k,j,i)**2 + &
                                   normz(k,j,i)**2 )

                    ! Normalize with direction from high c to low c
                    IF ( length > vofTol ) THEN
                        normx(k,j,i) = - normx(k,j,i) / length
                        normy(k,j,i) = - normy(k,j,i) / length
                        normz(k,j,i) = - normz(k,j,i) / length
                    ELSE
                        normx(k,j,i) = 0.0_realk
                        normy(k,j,i) = 0.0_realk
                        normz(k,j,i) = 0.0_realk
                    ENDIF

                ENDDO
            ENDDO
        ENDDO

    END SUBROUTINE comp_norm_vec

    !================================================================

    SUBROUTINE comp_ifc()
    !----------------------------------------------------------------
    !   What it does:
    !   Reconstruct interface on multi-grid level.
    !----------------------------------------------------------------

        ! Subroutine arguments
        ! None

        ! Local variables
        INTEGER(intk) :: n, igrid
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: c(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        REAL(realk), POINTER, CONTIGUOUS :: normx(:,:,:), normy(:,:,:), normz(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: alpha(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: isIfc(:,:,:)

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
            CALL get_fieldptr(normx, "NORMX", igrid)
            CALL get_fieldptr(normy, "NORMY", igrid)
            CALL get_fieldptr(normz, "NORMZ", igrid)
            CALL get_fieldptr(alpha, "ALPHA", igrid)
            CALL get_fieldptr(isIfc, "ISIFC", igrid)

            CALL comp_ifc_grd(kk, jj, ii, c, dx, dy, dz, &
                ddx, ddy, ddz, normx, normy, normz, alpha, isIfc)
        END DO

    END SUBROUTINE comp_ifc

    !================================================================

    SUBROUTINE comp_ifc_grd(kk, jj, ii, c, dx, dy, dz, ddx, ddy, ddz, normx, normy, normz, alpha, isIfc)
    !----------------------------------------------------------------
    !   What it does:
    !   Reconstruct interface on grid level.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: c(kk, jj, ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(out) :: normx(kk, jj, ii), normy(kk, jj, ii), normz(kk, jj, ii)
        REAL(realk), INTENT(out) :: alpha(kk, jj, ii)
        REAL(realk), INTENT(out) :: isIfc(kk, jj, ii)

        ! Local variables
        ! None

        CALL comp_isIfc_grd(kk, jj, ii, c, isIfc)
        CALL comp_norm_vec(normx, normy, normz, kk, jj, ii, c, dx, dy, dz)
        ! Fix for cells with no normal vector e.g. flotsam
        WHERE ( normx**2 + normy**2 + normz**2 < 0.5_realk ) isIfc = -1.0_realk
        CALL comp_alpha(alpha, kk, jj, ii, c, isIfc, ddx, ddy, ddz, normx, normy, normz)

    END SUBROUTINE comp_ifc_grd

    !================================================================

    SUBROUTINE comp_c_loc(cLoc, alpha, ddx, ddy, ddz, normx, normy, normz)
    !----------------------------------------------------------------
    !   What it does:
    !   Compute volume fraction given the interface parameters of a
    !   cell.
    !----------------------------------------------------------------

        ! Subroutine arguments
        REAL(realk), INTENT(out) :: cLoc
        REAL(realk), INTENT(in) :: alpha
        REAL(realk), INTENT(in) :: ddx, ddy, ddz
        REAL(realk), INTENT(in) :: normx, normy, normz

        ! Local variables
        REAL(realk) :: m1, m2, m3, c1, c2, c3
        REAL(realk) :: alphaStd

        CALL get_order(m1, m2, m3, c1, c2, c3, normx, normy, normz, ddx, ddy, ddz)
        alphaStd = alpha - mirror_shift(normx, normy, normz, ddx, ddy, ddz)
        CALL comp_c_std(m1, m2, m3, c1, c2, c3, alphaStd, cLoc)

    END SUBROUTINE comp_c_loc

    !================================================================

    SUBROUTINE comp_c_stg(q)
    !----------------------------------------------------------------
    !   What it does:
    !   Compute the volume fraction field for the staggered cells
    !   depending on q. The staggered cells are either moved by
    !   1/2 ddx, 1/2 ddy or 1/2 ddz.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: q

        ! Local variables
        CHARACTER(len=3) :: cFldName
        INTEGER(intk) :: n, igrid
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: c(:,:,:), cSq(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        REAL(realk), POINTER, CONTIGUOUS :: normx(:,:,:), normy(:,:,:), normz(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: alpha(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: isIfc(:,:,:)

        cFldName = "CS"//int2char(q)

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_fieldptr(c, "C", igrid)
            CALL get_fieldptr(cSq, cFldName, igrid)
            CALL get_fieldptr(ddx, "DDX", igrid)
            CALL get_fieldptr(ddy, "DDY", igrid)
            CALL get_fieldptr(ddz, "DDZ", igrid)
            CALL get_fieldptr(normx, "NORMX", igrid)
            CALL get_fieldptr(normy, "NORMY", igrid)
            CALL get_fieldptr(normz, "NORMZ", igrid)
            CALL get_fieldptr(alpha, "ALPHA", igrid)
            CALL get_fieldptr(isIfc, "ISIFC", igrid)

            CALL comp_c_stg_grd(kk, jj, ii, q, c, cSq, &
                ddx, ddy, ddz, normx, normy, normz, alpha, isIfc)
        END DO

    END SUBROUTINE comp_c_stg

    !================================================================

    SUBROUTINE comp_c_stg_grd(kk, jj, ii, q, c, cSq, &
        ddx, ddy, ddz, normx, normy, normz, alpha, isIfc)
    !----------------------------------------------------------------
    !   What it does:
    !   Compute the volume fraction field for the staggered cells
    !   depending on q. The staggered cells are either moved by
    !   1/2ddx, 1/2ddy or 1/2ddz.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, q
        REAL(realk), INTENT(in) :: c(kk, jj, ii)
        REAL(realk), INTENT(inout) :: cSq(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: normx(kk, jj, ii), normy(kk, jj, ii), normz(kk, jj, ii)
        REAL(realk), INTENT(in) :: alpha(kk, jj, ii)
        REAL(realk), INTENT(in) :: isIfc(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: kq, jq, iq
        REAL(realk) :: alphaOffset
        REAL(realk) :: alphaMi, ddxMi, ddyMi, ddzMi, normxMi, normyMi, normzMi, ddsMi, halfFractionMi
        REAL(realk) :: alphaPl, ddxPl, ddyPl, ddzPl, normxPl, normyPl, normzPl, ddsPl, halfFractionPl

        cSq = 0.0_realk

        CALL sel_ind(q, iq, jq, kq)

        DO i = 2, ii-2
            DO j = 2, jj-2
                DO k = 2, kk-2
                    alphaOffset = ( iq * normx(k,j,i) * ddx(i) + jq * normy(k,j,i) * ddy(j) + kq * normz(k,j,i) * ddz(k) ) / 2.0_realk

                    alphaMi = alpha(k,j,i) - alphaOffset
                    ddxMi = ( 1.0_realk - iq ) * ddx(i) + iq * ddx(i) / 2.0_realk
                    ddyMi = ( 1.0_realk - jq ) * ddy(j) + jq * ddy(j) / 2.0_realk
                    ddzMi = ( 1.0_realk - kq ) * ddz(k) + kq * ddz(k) / 2.0_realk
                    normxMi = normx(k,j,i)
                    normyMi = normy(k,j,i)
                    normzMi = normz(k,j,i)
                    ddsMi = iq * ddxMi + jq * ddyMi + kq * ddzMi

                    alphaPl = alpha(k+kq,j+jq,i+iq)
                    ddxPl = ( 1.0_realk - iq ) * ddx(i+iq) + iq * ddx(i+iq) / 2.0_realk
                    ddyPl = ( 1.0_realk - jq ) * ddy(j+jq) + jq * ddy(j+jq) / 2.0_realk
                    ddzPl = ( 1.0_realk - kq ) * ddz(k+kq) + kq * ddz(k+kq) / 2.0_realk
                    normxPl = normx(k+kq,j+jq,i+iq)
                    normyPl = normy(k+kq,j+jq,i+iq)
                    normzPl = normz(k+kq,j+jq,i+iq)
                    ddsPl = iq * ddxPl + jq * ddyPl + kq * ddzPl

                    IF ( isIfc(k,j,i) > 0.0_realk .AND. isIfc(k+kq,j+jq,i+iq) > 0.0_realk ) THEN
                        CALL comp_c_loc(halfFractionMi, alphaMi, ddxMi, ddyMi, ddzMi, normxMi, normyMi, normzMi)
                        CALL comp_c_loc(halfFractionPl, alphaPl, ddxPl, ddyPl, ddzPl, normxPl, normyPl, normzPl)
                        cSq(k,j,i) = ( halfFractionMi * ddsMi + halfFractionPl * ddsPl ) / ( ddsMi + ddsPl )
                    ELSE IF ( isIfc(k,j,i) > 0.0_realk .AND. isIfc(k+kq,j+jq,i+iq) < 0.0_realk ) THEN
                        CALL comp_c_loc(halfFractionMi, alphaMi, ddxMi, ddyMi, ddzMi, normxMi, normyMi, normzMi)
                        cSq(k,j,i) = ( halfFractionMi * ddsMi + c(k+kq,j+jq,i+iq) * ddsPl ) / ( ddsMi + ddsPl )
                    ELSE IF ( isIfc(k,j,i) < 0.0_realk .AND. isIfc(k+kq,j+jq,i+iq) > 0.0_realk ) THEN
                        CALL comp_c_loc(halfFractionPl, alphaPl, ddxPl, ddyPl, ddzPl, normxPl, normyPl, normzPl)
                        cSq(k,j,i) = ( c(k,j,i) * ddsMi + halfFractionPl * ddsPl ) / ( ddsMi + ddsPl )
                    ELSE
                        cSq(k,j,i) = ( c(k,j,i) * ddsMi + c(k+kq,j+jq,i+iq) * ddsPl ) / ( ddsMi + ddsPl )
                    ENDIF
                ENDDO
            ENDDO
        ENDDO

    END SUBROUTINE comp_c_stg_grd

    !================================================================

    PURE SUBROUTINE get_order(m1, m2, m3, c1, c2, c3, normx, normy, normz, ddx, ddy, ddz)
    !----------------------------------------------------------------
    !   What it does:
    !   Determines the order of the scaled normal values
    !   and writes it to m1-m3 and c1-c3.
    !----------------------------------------------------------------

        ! Subroutine arguments
        REAL(realk), INTENT(out) :: m1, m2, m3, c1, c2, c3
        REAL(realk), INTENT(in) :: normx, normy, normz
        REAL(realk), INTENT(in) :: ddx, ddy, ddz

        ! Local variables
        REAL(realk) :: scaledNorm1, scaledNorm2, scaledNorm3, tmp
        INTEGER(intk) :: i1, i2, i3, tmpi

        ! Set initial order
        scaledNorm1 = ABS(normx)*ddx
        scaledNorm2 = ABS(normy)*ddy
        scaledNorm3 = ABS(normz)*ddz

        i1 = 1; i2 = 2; i3 = 3

        ! normx*ddx >= normy*ddy => normx is m2 or m3 and normy is m1 or m2
        ! normx*ddx <  normy*ddy => normx is m1 or m2 and normy is m2 or m3
        IF (scaledNorm1 >= scaledNorm2) THEN
            tmp = scaledNorm1; scaledNorm1 = scaledNorm2; scaledNorm2 = tmp
            tmpi = i1; i1 = i2; i2 = tmpi
        ENDIF

        ! normx*ddx >= normz*ddz => normx is m3 and normz is m1 or m2
        ! normx*ddx <  normz*ddz => normz is m3 and normx is m1 (or m2)
        ! --OR--
        ! normy*ddy >= normz*ddz => normy is m3 and normz is m1 or m2
        ! normy*ddy <  normz*ddz => normz is m3 and normy is m1 (or m2)
        IF (scaledNorm2 >= scaledNorm3) THEN
            tmp = scaledNorm2; scaledNorm2 = scaledNorm3; scaledNorm3 = tmp
            tmpi = i2; i2 = i3; i3 = tmpi
        ENDIF

        ! normx*ddx >= normz*ddz => normx is m2 and normz is m1
        ! normx*ddx <  normz*ddz => normz is m2 and normx is m1
        ! --OR--
        ! normy*ddy >= normz*ddz => normy is m2 and normz is m1
        ! normy*ddy <  normz*ddz => normz is m2 and normy is m1
        ! --OR--
        ! normy*ddy >= normx*ddx => normy is m2 and normz is m1
        ! normy*ddy <  normx*ddx => normx is m2 and normy is m1
        IF (scaledNorm1 >= scaledNorm2) THEN
            tmp = scaledNorm1; scaledNorm1 = scaledNorm2; scaledNorm2 = tmp
            tmpi = i1; i1 = i2; i2 = tmpi
        ENDIF

        ! Assign new order to m1-m3 and c1-c3
        SELECT CASE (i1)
        CASE (1); m1 = ABS(normx); c1 = ddx
        CASE (2); m1 = ABS(normy); c1 = ddy
        CASE (3); m1 = ABS(normz); c1 = ddz
        END SELECT

        SELECT CASE (i2)
        CASE (1); m2 = ABS(normx); c2 = ddx
        CASE (2); m2 = ABS(normy); c2 = ddy
        CASE (3); m2 = ABS(normz); c2 = ddz
        END SELECT

        SELECT CASE (i3)
        CASE (1); m3 = ABS(normx); c3 = ddx
        CASE (2); m3 = ABS(normy); c3 = ddy
        CASE (3); m3 = ABS(normz); c3 = ddz
        END SELECT

    END SUBROUTINE get_order

    !================================================================

    SUBROUTINE comp_alpha(alpha, kk, jj, ii, c, isIfc, ddx, ddy, ddz, normx, normy, normz)
    !----------------------------------------------------------------
    !   What it does:
    !   Compute alpha given the normal vector and volume fraction of
    !   a cell.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(out) :: alpha(kk, jj, ii)
        REAL(realk), INTENT(in) :: c(kk, jj, ii)
        REAL(realk), INTENT(in) :: isIfc(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: normx(kk, jj, ii), normy(kk, jj, ii), normz(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i
        REAL(realk) :: m1, m2, m3, c1, c2, c3
        REAL(realk) :: alphaStd

        alpha = 0.0_realk
        DO i = 2, ii-1
            DO j = 2, jj-1
                DO k = 2, kk-1
                    IF ( isIfc(k,j,i) < 0.0_realk ) CYCLE
                    CALL get_order(m1, m2, m3, c1, c2, c3, normx(k,j,i), normy(k,j,i), normz(k,j,i), ddx(i), ddy(j), ddz(k))
                    CALL comp_alpha_std(m1, m2, m3, c1, c2, c3, alphaStd, c(k,j,i))
                    alpha(k,j,i) = alphaStd + mirror_shift(normx(k,j,i), normy(k,j,i), normz(k,j,i), ddx(i), ddy(j), ddz(k))
                ENDDO
            ENDDO
        ENDDO

    END SUBROUTINE comp_alpha

    !================================================================

    PURE SUBROUTINE comp_alpha_std(m1, m2, m3, c1, c2, c3, &
        alphaStd, cLoc)
    !----------------------------------------------------------------
    !   What it does:
    !   Solves the standart case for the plane constant alpha.
    !   Mirrors comp_c_std.
    !----------------------------------------------------------------

        ! Subroutine arguments
        REAL(realk), INTENT(in) :: m1, m2, m3, c1, c2, c3
        REAL(realk), INTENT(out) :: alphaStd
        REAL(realk), INTENT(in) :: cLoc

        ! Local variables
        LOGICAL :: isUpper
        REAL(realk) :: mc1, mc2, mc3
        REAL(realk) :: alphaMax, alphaLow, volLow
        REAL(realk) :: area, areaCrit
        REAL(realk) :: V1, V2, V3
        REAL(realk) :: a0, a1, a2
        REAL(realk) :: qo, po
        REAL(realk) :: theta

        ! Scaled normal components
        mc1 = m1*c1
        mc2 = m2*c2
        mc3 = m3*c3
        alphaMax = comp_alpha_max(mc1, mc2, mc3)

        ! Fold input onto the lower half
        isUpper = ( cLoc > 0.5_realk )
        volLow = MIN(cLoc, 1.0_realk - cLoc)*c1*c2*c3

        IF ( mc1 < vofTol ) THEN
            IF ( mc2 < vofTol ) THEN
                ! 1D
                alphaLow = volLow/(c1*c2)
            ELSE
                ! 2D
                areaCrit = 0.5_realk*c2**2*m2/m3
                area = volLow/c1
                IF ( area < areaCrit ) THEN
                    alphaLow = SQRT(2.0_realk*area*m2*m3)
                ELSE
                    alphaLow = m3/c2*area + mc2*0.5_realk
                ENDIF
            ENDIF
        ELSE
            ! 3D
            V1 = mc1**2*c1/(MAX(6.0_realk*m2*m3, vofTol))
            V2 = V1 + c1*c2*(mc2 - mc1)/(2.0_realk*m3)
            IF ( mc3 < mc1 + mc2 ) THEN
                V3 = (mc3**2*(3.0_realk*(mc1 + mc2) - mc3) + &
                    mc1**2*(mc1 - 3.0_realk*mc3) + &
                    mc2**2*(mc2 - 3.0_realk*mc3))/ &
                    (6.0_realk*m1*m2*m3)
            ELSE
                V3 = c1*c2*(mc1 + mc2)/(2.0_realk*m3)
            ENDIF

            IF ( volLow < V1 ) THEN
                alphaLow = (6.0_realk*m1*m2*m3*volLow)** &
                    (1.0_realk/3.0_realk)
            ELSE IF ( volLow < V2 ) THEN
                alphaLow = 0.5_realk*(mc1 + SQRT(mc1**2 + &
                    8.0_realk*m2*m3*(volLow - V1)/c1))
            ELSE IF ( volLow < V3 ) THEN
                a2 = -3.0_realk*(mc1 + mc2)
                a1 = 3.0_realk*(mc1**2 + mc2**2)
                a0 = -(mc1**3 + mc2**3) + 6.0_realk*m1*m2*m3*volLow
                po = a1/3.0_realk - a2**2/9.0_realk
                qo = (a1*a2 - 3.0_realk*a0)/6.0_realk - &
                    a2**3/27.0_realk
                theta = ACOS(qo/(-po*SQRT(-po)))/3.0_realk
                alphaLow = SQRT(-po)*(SQRT(3.0_realk)*SIN(theta) - &
                    COS(theta)) - a2/3.0_realk
            ELSE IF ( volLow >= V3 .AND. mc3 <= mc1 + mc2 ) THEN
                a2 = -3.0_realk/2.0_realk*(mc1 + mc2 + mc3)
                a1 = 3.0_realk/2.0_realk*(mc1**2 + mc2**2 + mc3**2)
                a0 = -0.5_realk*(mc1**3 + mc2**3 + mc3**3) + &
                    3.0_realk*m1*m2*m3*volLow
                po = a1/3.0_realk - a2**2/9.0_realk
                qo = (a1*a2 - 3.0_realk*a0)/6.0_realk - &
                    a2**3/27.0_realk
                theta = ACOS(qo/(-po*SQRT(-po)))/3.0_realk
                alphaLow = SQRT(-po)*(SQRT(3.0_realk)*SIN(theta) - &
                    COS(theta)) - a2/3.0_realk
            ELSE IF ( volLow >= V3 .AND. mc3 > mc1 + mc2 ) THEN
                alphaLow = m3*volLow/(c1*c2) + (mc1 + mc2)*0.5_realk
            ENDIF
        ENDIF

        ! Unfold result onto the full range
        alphaStd = MERGE(alphaMax - alphaLow, alphaLow, isUpper)

    END SUBROUTINE comp_alpha_std

    !================================================================

    PURE SUBROUTINE comp_c_std(m1, m2, m3, c1, c2, c3, &
        alphaStd, cLoc)
    !----------------------------------------------------------------
    !   What it does:
    !   Solves the standart case for the volume fraction.
    !   Mirrors comp_alpha_std.
    !----------------------------------------------------------------

        ! Subroutine arguments
        REAL(realk), INTENT(in) :: m1, m2, m3, c1, c2, c3
        REAL(realk), INTENT(in) :: alphaStd
        REAL(realk), INTENT(out) :: cLoc

        ! Local variables
        LOGICAL :: isUpper
        REAL(realk) :: mc1, mc2, mc3
        REAL(realk) :: alphaMax, alphaLow, volLow
        REAL(realk) :: area, areaCrit
        REAL(realk) :: V1
        REAL(realk) :: cLow

        ! Scaled normal components
        mc1 = m1*c1
        mc2 = m2*c2
        mc3 = m3*c3
        alphaMax = comp_alpha_max(mc1, mc2, mc3)

        ! Fold input onto the lower half
        isUpper = ( alphaStd > 0.5_realk*alphaMax )
        alphaLow = MAX(MIN(alphaStd, alphaMax - alphaStd), 0.0_realk)

        IF ( mc1 < vofTol ) THEN
            IF ( mc2 < vofTol ) THEN
                ! 1D
                volLow = alphaLow*(c1*c2)
            ELSE
                ! 2D
                areaCrit = 0.5_realk*c2**2*m2/m3
                IF ( alphaLow < mc2 ) THEN
                    area = 0.5_realk*alphaLow**2/(m2*m3)
                ELSE
                    area = c2*alphaLow/m3 - areaCrit
                ENDIF
                volLow = area*c1
            ENDIF
        ELSE
            ! 3D
            V1 = mc1**2*c1/(MAX(6.0_realk*m2*m3, vofTol))
            IF ( alphaLow < mc1 ) THEN
                volLow = alphaLow**3/(6.0_realk*m1*m2*m3)
            ELSE IF ( alphaLow < mc2 ) THEN
                volLow = (alphaLow*c1*(alphaLow - mc1))/ &
                    (2.0_realk*m2*m3) + V1
            ELSE IF ( alphaLow < MIN(mc1 + mc2, mc3) ) THEN
                volLow = (alphaLow**2*(3.0_realk* &
                    (mc1 + mc2) - alphaLow) + &
                    mc1**2*(mc1 - 3.0_realk*alphaLow) + &
                    mc2**2*(mc2 - 3.0_realk*alphaLow))/ &
                    (6.0_realk*m1*m2*m3)
            ELSE IF ( alphaLow >= MIN(mc1 + mc2, mc3) .AND. &
                      mc3 <= mc1 + mc2 ) THEN
                volLow = (alphaLow**2*(3.0_realk* &
                    (mc1 + mc2 + mc3) - 2.0_realk*alphaLow) + &
                    mc1**2*(mc1 - 3.0_realk*alphaLow) + &
                    mc2**2*(mc2 - 3.0_realk*alphaLow) + &
                    mc3**2*(mc3 - 3.0_realk*alphaLow))/ &
                    (6.0_realk*m1*m2*m3)
            ELSE IF ( alphaLow >= MIN(mc1 + mc2, mc3) .AND. &
                      mc3 > mc1 + mc2 ) THEN
                volLow = (c1*c2*(2.0_realk*alphaLow - &
                    (mc1 + mc2)))/(2.0_realk*m3)
            ENDIF
        ENDIF

        ! Unfold result onto the full range
        cLow = volLow/(c1*c2*c3)
        cLoc = MERGE(1.0_realk - cLow, cLow, isUpper)

    END SUBROUTINE comp_c_std

    !================================================================

    PURE FUNCTION mirror_shift(normx, normy, normz, ddx, ddy, ddz) RESULT(shift)
    !----------------------------------------------------------------
    !   What it does:
    !   Each direction with a negative normal component
    !   mirrors the cell and shifts the plane constant by 
    !   norm(.)*dd(.).
    !----------------------------------------------------------------

        REAL(realk), INTENT(in) :: normx, normy, normz
        REAL(realk), INTENT(in) :: ddx, ddy, ddz
        REAL(realk) :: shift

        shift = MIN(normx*ddx, 0.0_realk) + &
                MIN(normy*ddy, 0.0_realk) + &
                MIN(normz*ddz, 0.0_realk)

    END FUNCTION mirror_shift

    !================================================================

    PURE FUNCTION comp_alpha_max(mc1, mc2, mc3) RESULT(alphaMax)
    !----------------------------------------------------------------
    !   What it does:
    !   Compute theoretical maximum alpha value.
    !----------------------------------------------------------------

        REAL(realk), INTENT(in) :: mc1, mc2, mc3
        REAL(realk) :: alphaMax

        alphaMax = MERGE(mc1, 0.0_realk, mc1 >= vofTol) + &
                   MERGE(mc2, 0.0_realk, mc2 >= vofTol) + &
                   MERGE(mc3, 0.0_realk, mc3 >= vofTol)

    END FUNCTION comp_alpha_max

END MODULE mph_plic_mod