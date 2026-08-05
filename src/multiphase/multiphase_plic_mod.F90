    !====================================================================
    !  Module: multiphase_plic_mod
    !
    !  Responsibilities:
    !     - Tracks the cells containing an interfac
    !     - Computes the interface normal vector within a cell
    !     - Computes the distance of the interface alpha within a cell
    !
    !  Author:      Valentin Ruhs
    !  Created:     2026-02
    !  Last update: 2026-02
    !
    !====================================================================

    MODULE multiphase_plic_mod

    USE precision_mod, ONLY: intk, realk
    USE multiphase_utils_mod, ONLY: get_spatial_indices
    USE multiphasecore_mod, ONLY: tol
        
    IMPLICIT NONE
    PRIVATE 

    PUBLIC :: init_multiphase_plic, finish_multiphase_plic, iface_reconstruction, &
            comp_frac, comp_stag_frac, track_iface, track_iface_vic

    CONTAINS

    SUBROUTINE init_multiphase_plic()

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        continue
    END SUBROUTINE init_multiphase_plic

    !================================================================

    SUBROUTINE finish_multiphase_plic()

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None
        
        continue
    END SUBROUTINE finish_multiphase_plic

    !================================================================

    SUBROUTINE track_iface(isIface, kk, jj, ii, vff)
    !----------------------------------------------------------------
    !   What it does:
    !   Identifies which of the cells in the domain contains a volume
    !   fraction of two fluids. These cells have to be taken into
    !   account when reconstructing interfaces. Therefore the
    !   variable containign this information is called isIface.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        LOGICAL, INTENT(out) :: isIface(kk, jj, ii)
        REAL(realk), INTENT(in) :: vff(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i

        isIface = .FALSE.
        DO i = 1, ii
            DO j = 1, jj
                DO k = 1, kk
                    IF ( vff(k,j,i) > tol .AND. vff(k,j,i) < 1.0_realk - tol ) THEN
                        isIface(k,j,i) = .TRUE.
                    ENDIF
                ENDDO
            ENDDO
        ENDDO

    END SUBROUTINE track_iface

    !================================================================

    SUBROUTINE track_iface_vic(isIfaceVic, kk, jj, ii, isIface)
    !----------------------------------------------------------------
    !   What it does:
    !   Creates the isIfaceVic array, which stores TRUE
    !   when a cell is considered "near" an interface cell. Since
    !   indices with (.)-2 and (.)+2 are used in the momentum
    !   advection a 5x5x5 volume is considered to be "near" to an 
    !   interface cell. 
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        LOGICAL, INTENT(out) :: isIfaceVic(kk, jj, ii)
        LOGICAL, INTENT(in) :: isIface(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i, vic

        vic = 2
        isIfaceVic = .FALSE.

        ! Find cells with interface
        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    ! If there is an interface, set the 5x5x5 neighbouring cells to near interface
                    IF ( isIface(k,j,i) ) THEN
                        isIfaceVic(k-vic:k+vic,j-vic:j+vic,i-vic:i+vic) = .TRUE.
                    ENDIF
                ENDDO
            ENDDO
        ENDDO

    END SUBROUTINE track_iface_vic

    !================================================================

    SUBROUTINE comp_norm_vec(normx, normy, normz, kk, jj, ii, vff, dx, dy, dz)
    !----------------------------------------------------------------
    !   What it does:
    !   Computes the normal vector components normx, normy and normz
    !   of the gradient of the volume fraction function vff. The 
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
        REAL(realk), INTENT(in) :: vff(kk, jj, ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)

        ! Local variables
        INTEGER(intk) :: k, j, i, d1, d2
        INTEGER(intk), PARAMETER :: stcl(-1:1, -1:1) = &
            RESHAPE([1, 2, 1, 2, 4, 2, 1, 2, 1], [3,3])
        REAL(realk) :: length
        REAL(realk) :: sumx, sumy, sumz

        ! Loop over cells
        DO i = 2, ii-1
            DO j = 2, jj-1
                DO k = 2, kk-1

                    sumx = 0.0_realk
                    sumy = 0.0_realk
                    sumz = 0.0_realk

                    ! Loop over stcl Stcl
                    DO d2 = -1, 1
                        DO d1 = -1, 1
                            sumx = sumx + stcl(d1, d2) * &
                                ( vff(k+d1, j+d2, i+1) - vff(k+d1, j+d2, i-1) )
                            sumy = sumy + stcl(d1, d2) * &
                                ( vff(k+d1, j+1, i+d2) - vff(k+d1, j-1, i+d2) )
                            sumz = sumz + stcl(d1, d2) * &
                                ( vff(k+1, j+d1, i+d2) - vff(k-1, j+d1, i+d2) )
                        ENDDO
                    ENDDO

                    normx(k,j,i) = sumx / ( SUM(stcl) * (dx(i-1)+dx(i)) )
                    normy(k,j,i) = sumy / ( SUM(stcl) * (dy(j-1)+dy(j)) )
                    normz(k,j,i) = sumz / ( SUM(stcl) * (dz(k-1)+dz(k)) )

                    ! Calculate normal vector length
                    length = SQRT( normx(k,j,i)**2.0_realk + &
                                   normy(k,j,i)**2.0_realk + &
                                   normz(k,j,i)**2.0_realk )

                    ! Normalize with direction from high vff to low vff
                    IF ( length > tol ) THEN
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

    SUBROUTINE iface_reconstruction(kk, jj, ii, vff, dx, dy, dz, ddx, ddy, ddz, normx, normy, normz, alpha)
    !----------------------------------------------------------------
    !   What it does:
    !   The subroutine is just a wrapper for the subroutines, which
    !   are used to reconstrunct the interface in cells.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: vff(kk, jj, ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(out) :: normx(kk, jj, ii), normy(kk, jj, ii), normz(kk, jj, ii)
        REAL(realk), INTENT(out) :: alpha(kk, jj, ii)

        ! Local variables
        LOGICAL :: isIface(kk, jj, ii)

        CALL track_iface(isIface, kk, jj, ii, vff)
        CALL comp_norm_vec(normx, normy, normz, kk, jj, ii, vff, dx, dy, dz)
        CALL comp_alph(alpha, kk, jj, ii, vff, isIface, ddx, ddy, ddz, normx, normy, normz)

    END SUBROUTINE iface_reconstruction

    !================================================================

    SUBROUTINE comp_alph(alpha, kk, jj, ii, vff, isIface, ddx, ddy, ddz, normx, normy, normz)
    !----------------------------------------------------------------
    !   What it does:
    !   This subroutine calculates the alpha value for PLIC. The 
    !   alpha value describes the distance of the interface in a cell
    !   from a defined reference (left bottom front corner).
    !   1. Assign norm(.) to m1, m2 and m3 and c1-c3 respectively
    !   2. Transform vff to a actual volume in bounds [0,0.5] * dV
    !   3. Solve the standart case for alpha
    !   4. If necessary, transform alpha to its conjugate 
    !      alphaMax - alpha
    !   5. If necessary, transform alpha regarding to its negative
    !      normal vector components
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(out) :: alpha(kk, jj, ii)
        REAL(realk), INTENT(in) :: vff(kk, jj, ii)
        LOGICAL, INTENT(in) :: isIface(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: normx(kk, jj, ii), normy(kk, jj, ii), normz(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i
        REAL(realk) :: m1, m2, m3, c1, c2, c3
        REAL(realk) :: alphaMax(kk, jj, ii)


        ! Loop over cells
        DO i = 2, ii-1
            DO j = 2, jj-1
                DO k = 2, kk-1

                    ! Only calculate interface for intersected cells
                    IF ( .NOT. isIface(k,j,i) ) THEN 
                        CYCLE
                    ENDIF

                    ! 1. Assign norm(.) to m1, m2 and m3 and c1-c3 respectively
                    ! To enhance performance consider inlining
                    CALL get_order(m1, m2, m3, c1, c2, c3, &
                        normx(k,j,i), normy(k,j,i), normz(k,j,i), ddx(i), ddy(j), ddz(k))

                    ! 2. Transform vff to a actual volume in bounds [0,0.5] * dV
                    ! 3. Solve the standart cases for alpha
                    ! Source: R. Scardovelli und S. Zaleski, „Analytical Relations Connecting Linear Interfaces and Volume Fractions in Rectangular Grids“,
                    !         Journal of Computational Physics, Bd. 164, Nr. 1, S. 228–237, Okt. 2000, doi: 10.1006/jcph.2000.6567.
                    ! To enhance performance consider inlining
                    CALL comp_alph_std(m1, m2, m3, c1, c2, c3, alpha(k,j,i), &
                        alphaMax(k,j,i), vff(k,j,i), ddx(i), ddy(j), ddz(k))

                    ! 4. If necessary, transform alpha back to volume bounds [0,1] * dV
                    ! If the volume fraction function has a value above 0.5 the "inverse problem" is solved. Therefore, the result is no longer 
                    ! alpha, but alphaMax - alpha. It can be seen as a rotation of the voxel. This is the inverse rotation (see comp_alph_std)
                    IF ( vff(k,j,i) > 0.5_realk ) THEN
                        alpha(k,j,i) = alphaMax(k,j,i) - alpha(k,j,i)
                    ENDIF

                    ! 5. If necessary, transform alpha regarding to its negative normal vector components
                    ! If one of the normal vector components is negative, a mirrored case is solved. Therefore, the solution
                    ! has to be transformed back (see get_order)
                    IF ( normx(k,j,i) < 0.0_realk ) THEN
                        alpha(k,j,i) = alpha(k,j,i) + ddx(i)*normx(k,j,i)
                    ENDIF

                    IF ( normy(k,j,i) < 0.0_realk ) THEN
                        alpha(k,j,i) = alpha(k,j,i) + ddy(j)*normy(k,j,i)
                    ENDIF

                    IF ( normz(k,j,i) < 0.0_realk ) THEN
                        alpha(k,j,i) = alpha(k,j,i) + ddz(k)*normz(k,j,i)
                    ENDIF

                ENDDO
            ENDDO
        ENDDO

    END SUBROUTINE comp_alph

    !================================================================

    SUBROUTINE comp_frac(cellProportion, alpha, vff, ddx, ddy, ddz, normx, normy, normz)
    !----------------------------------------------------------------
    !   What it does:
    !   This subroutine calculates the value of the volume fraction 
    !   function vff in the voxel ddx*ddy*ddz, given alpha. 
    !   1. Assign norm(.) to m1, m2 and m3 and c1-c3 respectively
    !   2. If necessary, transform alpha regarding to its negative
    !      normal vector components
    !   3. If necessary, transform alpha to its conjugate 
    !      alphaMax - alpha 
    !   4. Solve the standart case for vol
    !----------------------------------------------------------------

        ! Subroutine arguments
        REAL(realk), INTENT(out) :: cellProportion
        REAL(realk), INTENT(in) :: alpha
        REAL(realk), INTENT(in) :: vff
        REAL(realk), INTENT(in) :: ddx, ddy, ddz
        REAL(realk), INTENT(in) :: normx, normy, normz

        ! Input variables
        REAL(realk) :: m1, m2, m3, c1, c2, c3
        REAL(realk) :: alphaStd

        ! 1. Assign norm(.) to m1, m2 and m3 and c1-c3 respectively
        ! To enhance performance consider inlining
        CALL get_order(m1, m2, m3, c1, c2, c3, normx, normy, normz, ddx, ddy, ddz)

        ! 2. If necessary, transform alpha regarding to its negative normal vector components
        alphaStd = alpha
        IF ( normx < 0.0_realk ) THEN
            alphaStd = alphaStd - ddx*normx
        ENDIF

        IF ( normy < 0.0_realk ) THEN
            alphaStd = alphaStd - ddy*normy
        ENDIF

        IF ( normz < 0.0_realk ) THEN
            alphaStd = alphaStd - ddz*normz
        ENDIF

        ! 3. If necessary, transform alpha to its conjugate alphaMax - alpha
        ! 4. Solve the standart case for vol
        CALL comp_frac_std(m1, m2, m3, c1, c2, c3, alphaStd, cellProportion)

    END SUBROUTINE comp_frac

    !================================================================

    SUBROUTINE comp_stag_frac(kk, jj, ii, q, vff, vffStag, ddx, ddy, ddz, normx, normy, normz, alpha, isIface)
    !----------------------------------------------------------------
    !   What it does:
    !   Compute the volume fraction field for the staggered cells
    !   depending on q. The staggered cells are either moved by
    !   1/2 ddx, 1/2 ddy or 1/2 ddz.
    !----------------------------------------------------------------

    ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, q
        REAL(realk), INTENT(in) :: vff(kk, jj, ii)
        REAL(realk), INTENT(out) :: vffStag(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: normx(kk, jj, ii), normy(kk, jj, ii), normz(kk, jj, ii)
        REAL(realk), INTENT(in) :: alpha(kk, jj, ii)
        LOGICAL, INTENT(in) :: isIface(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: kq, jq, iq
        REAL(realk) :: alphaOffset
        REAL(realk) :: alphaMi, vffMi, ddxMi, ddyMi, ddzMi, normxMi, normyMi, normzMi, ddsMi, halfFractionMi
        REAL(realk) :: alphaPl, vffPl, ddxPl, ddyPl, ddzPl, normxPl, normyPl, normzPl, ddsPl, halfFractionPl

        vffStag = 0.0_realk

        CALL get_spatial_indices(kk, jj, ii, q, iq, jq, kq)

        DO i = 2, ii-2
            DO j = 2, jj-2
                DO k = 2, kk-2
                    alphaOffset = ( iq * normx(k,j,i) * ddx(i) + jq * normy(k,j,i) * ddy(j) + kq * normz(k,j,i) * ddz(k) ) / 2.0_realk

                    alphaMi = alpha(k,j,i) - alphaOffset
                    vffMi = vff(k,j,i)
                    ddxMi = ( 1.0_realk - iq ) * ddx(i) + iq * ddx(i) / 2.0_realk
                    ddyMi = ( 1.0_realk - jq ) * ddy(j) + jq * ddy(j) / 2.0_realk
                    ddzMi = ( 1.0_realk - kq ) * ddz(k) + kq * ddz(k) / 2.0_realk
                    normxMi = normx(k,j,i)
                    normyMi = normy(k,j,i)
                    normzMi = normz(k,j,i)
                    ddsMi = iq * ddxMi + jq * ddyMi + kq * ddzMi

                    alphaPl = alpha(k+kq,j+jq,i+iq)
                    vffPl = vff(k+kq,j+jq,i+iq)
                    ddxPl = ( 1.0_realk - iq ) * ddx(i+iq) + iq * ddx(i+iq) / 2.0_realk
                    ddyPl = ( 1.0_realk - jq ) * ddy(j+jq) + jq * ddy(j+jq) / 2.0_realk
                    ddzPl = ( 1.0_realk - kq ) * ddz(k+kq) + kq * ddz(k+kq) / 2.0_realk
                    normxPl = normx(k+kq,j+jq,i+iq)
                    normyPl = normy(k+kq,j+jq,i+iq)
                    normzPl = normz(k+kq,j+jq,i+iq)
                    ddsPl = iq * ddxPl + jq * ddyPl + kq * ddzPl

                    IF ( isIface(k,j,i) .AND. isIface(k+kq,j+jq,i+iq) ) THEN
                        CALL comp_frac(halfFractionMi, alphaMi, vffMi, ddxMi, ddyMi, ddzMi, normxMi, normyMi, normzMi)
                        CALL comp_frac(halfFractionPl, alphaPl, vffPl, ddxPl, ddyPl, ddzPl, normxPl, normyPl, normzPl)
                        vffStag(k,j,i) = ( halfFractionMi * ddsMi + halfFractionPl * ddsPl ) / ( ddsMi + ddsPl )
                    ELSE IF ( isIface(k,j,i) .AND. .NOT. isIface(k+kq,j+jq,i+iq) ) THEN
                        CALL comp_frac(halfFractionMi, alphaMi, vffMi, ddxMi, ddyMi, ddzMi, normxMi, normyMi, normzMi)
                        vffStag(k,j,i) = ( halfFractionMi * ddsMi + vff(k+kq,j+jq,i+iq) * ddsPl ) / ( ddsMi + ddsPl )
                    ELSE IF ( .NOT. isIface(k,j,i) .AND. isIface(k+kq,j+jq,i+iq) ) THEN
                        CALL comp_frac(halfFractionPl, alphaPl, vffPl, ddxPl, ddyPl, ddzPl, normxPl, normyPl, normzPl)
                        vffStag(k,j,i) = ( vff(k,j,i) * ddsMi + halfFractionPl * ddsPl ) / ( ddsMi + ddsPl )
                    ELSE
                        vffStag(k,j,i) = ( vff(k,j,i) * ddsMi + vff(k+kq,j+jq,i+iq) * ddsPl ) / ( ddsMi + ddsPl )
                    ENDIF
                ENDDO
            ENDDO
        ENDDO

    END SUBROUTINE comp_stag_frac

    !================================================================

    PURE SUBROUTINE get_order(m1, m2, m3, c1, c2, c3, normx, normy, normz, ddx, ddy, ddz)
    !----------------------------------------------------------------
    !   What it does:
    !   This is a pure subroutine to enhance the performance by 
    !   inlining. It determines the order of the scaled normal values
    !   and writes it to m1-m3 and c1-c3 respectively.
    !----------------------------------------------------------------

        ! Subroutine arguments
        REAL(realk), INTENT(out) :: m1, m2, m3, c1, c2, c3
        REAL(realk), INTENT(in) :: normx, normy, normz
        REAL(realk), INTENT(in) :: ddx, ddy, ddz

        ! Local variables
        REAL(realk) :: scaledNorm1, scaledNorm2, scaledNorm3, tmp
        INTEGER(intk) :: i1, i2, i3, tmpi

        ! Set initial order
        scaledNorm1 = abs(normx)*ddx
        scaledNorm2 = abs(normy)*ddy
        scaledNorm3 = abs(normz)*ddz

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

        ! Assign new order to m1-m3 and c1-c3 respectively
        ! The absolute value of norm(.) is a mirror transform (see comp_alph 5.).
        SELECT CASE (i1)
        CASE (1); m1 = abs(normx); c1 = ddx
        CASE (2); m1 = abs(normy); c1 = ddy
        CASE (3); m1 = abs(normz); c1 = ddz
        END SELECT

        SELECT CASE (i2)
        CASE (1); m2 = abs(normx); c2 = ddx
        CASE (2); m2 = abs(normy); c2 = ddy
        CASE (3); m2 = abs(normz); c2 = ddz
        END SELECT

        SELECT CASE (i3)
        CASE (1); m3 = abs(normx); c3 = ddx
        CASE (2); m3 = abs(normy); c3 = ddy
        CASE (3); m3 = abs(normz); c3 = ddz
        END SELECT

    END SUBROUTINE get_order

    !================================================================

    SUBROUTINE comp_alph_std(m1, m2, m3, c1, c2, c3, alphaStd, alphaMax, vff, ddx, ddy, ddz)
    !----------------------------------------------------------------
    !   What it does:
    !   This is a pure subroutine to enhance the performance by 
    !   inlining. It solves the cubic equation
    !
    !   vol = 1 / (6 * m1 * m2 * m3) * [alphaStd^3
    !   - sum_{j=1..3} H(alphaStd - m_j * c_j) * (alphaStd - m_j * c_j)^3
    !   + sum_{j=1..3} H(alphaStd - alphaMax + m_j * c_j) * 
    !   (alphaStd - alphaMax + m_j * c_j)^3]
    !
    !   where alphaMax = m1*c1 + m2*c2 + m3*c3
    !
    !   for alphaStd. This is done for the standart cases:
    !       - 0 = mc1 = mc2 < mc3 (one-dimensional)
    !       - 0 = mc1 < mc3 < mc3 (two-dimensional)
    !       - mc1 < mc2 < mc3 (three-dimensional)
    !----------------------------------------------------------------

        ! Subroutine arguments
        REAL(realk), INTENT(in) :: m1, m2, m3, c1, c2, c3
        REAL(realk), INTENT(out) :: alphaStd, alphaMax
        REAL(realk), INTENT(in) :: vff
        REAL(realk), INTENT(in) :: ddx, ddy, ddz
        
        ! Local variables
        REAL(realk) :: mc1, mc2, mc3
        REAL(realk) :: vol
        REAL(realk) :: baseArea, criticalBaseArea
        REAL(realk) :: V1, V2, V3
        REAL(realk) :: a0, a1, a2
        REAL(realk) :: qo, po
        REAL(realk) :: theta

        ! Transform vff to a actual volume in bounds [0,0.5] * dV
        ! Rotate voxel into standart configuration (see comp_alph 4.)
        vol = MIN(vff, 1.0_realk - vff) * ddx * ddy * ddz
        
        ! Solve the standart cases for alphaStd
        ! Source: R. Scardovelli und S. Zaleski, „Analytical Relations Connecting Linear Interfaces and Volume Fractions in Rectangular Grids“,
        !         Journal of Computational Physics, Bd. 164, Nr. 1, S. 228–237, Okt. 2000, doi: 10.1006/jcph.2000.6567.
        mc1 = m1*c1
        mc2 = m2*c2
        mc3 = m3*c3
        
        IF ( mc1 < tol ) THEN
            IF ( mc2 < tol ) THEN
                ! One-dimensional case
                alphaMax = mc3
                alphaStd = vol / (c1*c2)
            ELSE
                ! Two-dimensional cases
                alphaMax = mc2 + mc3
                
                ! actual base area
                baseArea = vol / c1
                
                ! When the critical base area is exceeded the volume shape transforms to a chamfered rectangle prism instead of triangular prism
                criticalBaseArea = 1.0_realk/2.0_realk * c2**2.0_realk * m2/m3
                
                IF ( baseArea < criticalBaseArea ) THEN
                    ! Here both interception lines of the interface with the coordinate axis are within the cell => triangular prism
                    alphaStd = SQRT(2.0_realk * baseArea * m2 * m3)
                ELSE
                    ! Here one interception line (with the c2 axis) is outside the cell => chamfered rectangle prism
                    alphaStd = (m3) / (c2) * baseArea + (mc2) / 2.0_realk
                ENDIF
            ENDIF
        ELSE
            ! Three-dimensional cases
            alphaMax = mc1 + mc2 + mc3

            ! Define interval boundaries V1, V2, V3
            V1 = mc1**2.0_realk * c1 / ( MAX(6.0_realk * m2 * m3, tol) )
            V2 = V1 + c1 * c2 * ( mc2 - mc1 ) / ( 2.0_realk * m3 )
            IF ( mc3 < mc1 + mc2 ) THEN
                V3 = ( mc3**2.0_realk * ( 3.0_realk * ( mc1 + mc2 ) - mc3 ) + &
                        mc1**2.0_realk * ( mc1 - 3.0_realk * mc3 ) + &
                        mc2**2.0_realk * ( mc2 - 3.0_realk * mc3 ) ) / &
                        ( 6.0_realk * m1 * m2 * m3 )
            ELSE
                V3 = c1 * c2 * ( mc1 + mc2 ) / ( 2.0_realk * m3 )
            ENDIF
            
            ! Calculate alphaStd dependent on V1, V2 and V3
            IF ( vol < V1 ) THEN
                alphaStd = ( 6.0_realk * m1 * m2 * m3 * vol )**( 1.0_realk/3.0_realk )
            ELSE IF ( vol < V2 ) THEN
                alphaStd = 1.0_realk/2.0_realk * ( mc1 + SQRT(mc1**2.0_realk + 8.0_realk * m2 * m3 * (vol - V1) / c1) )
            ELSE IF ( vol < V3 ) THEN
                a2 = - 3.0_realk * ( mc1 + mc2 )
                a1 = 3.0_realk * ( mc1**2.0_realk + mc2**2.0_realk )
                a0 = - (mc1**3.0_realk + mc2**3.0_realk) + 6.0_realk * m1 * m2 * m3 * vol
                po = a1 / 3.0_realk - a2**2.0_realk / 9.0_realk
                qo = ( a1 * a2 - 3.0_realk * a0 ) / 6.0_realk - a2**3.0_realk / 27.0_realk
                
                theta = ACOS(qo / (-po)**1.5_realk) / 3.0_realk
                alphaStd = SQRT(-po) * ( SQRT(3.0_realk) * SIN(theta) - COS(theta) ) - a2 / 3.0_realk
            ELSE IF ( vol >= V3 .AND. mc3 <= mc1 + mc2 ) THEN
                a2 = - 3.0_realk/2.0_realk * ( mc1 + mc2 + mc3 )
                a1 = 3.0_realk/2.0_realk * ( mc1**2.0_realk + mc2**2.0_realk + mc3**2.0_realk )
                a0 = - 1.0_realk/2.0_realk * ( mc1**3.0_realk + mc2**3.0_realk + mc3**3.0_realk ) + 3.0_realk * m1 * m2 * m3 * vol
                po = a1 / 3.0_realk - a2**2.0_realk / 9.0_realk
                qo = ( a1 * a2 - 3.0_realk * a0 ) / 6.0_realk - a2**3.0_realk / 27.0_realk
                
                theta = ACOS(qo / (-po)**1.5_realk) / 3.0_realk
                alphaStd = SQRT(-po) * ( SQRT(3.0_realk) * SIN(theta) - COS(theta) ) - a2 / 3.0_realk
            ELSE IF ( vol >= V3 .AND. mc3 > mc1 + mc2 ) THEN
                alphaStd = m3 * vol / ( c1 * c2 ) + ( mc1 + mc2 ) / 2.0_realk
            ENDIF
        ENDIF

    END SUBROUTINE comp_alph_std

    !================================================================

    PURE SUBROUTINE comp_frac_std(m1, m2, m3, c1, c2, c3, alphaLoc, cellProportion)
    !----------------------------------------------------------------
    !   What it does:
    !   This is a pure subroutine to enhance the performance by 
    !   inlining. It solves the cubic equation
    !
    !   vol = 1 / (6 * m1 * m2 * m3) * [alphaStd^3
    !   - sum_{j=1..3} H(alphaStd - m_j * c_j) * (alphaStd - m_j * c_j)^3
    !   + sum_{j=1..3} H(alphaStd - alphaMax + m_j * c_j) * 
    !   (alphaStd - alphaMax + m_j * c_j)^3]
    !
    !   where alphaMax = m1*c1 + m2*c2 + m3*c3
    !
    !   for vol. This is done for the standart cases:
    !       - 0 = mc1 = mc2 < mc3 (one-dimensional)
    !       - 0 = mc1 < mc3 < mc3 (two-dimensional)
    !       - mc1 < mc2 < mc3 (three-dimensional)
    !----------------------------------------------------------------

        ! Subroutine arguments
        REAL(realk), INTENT(in) :: m1, m2, m3, c1, c2, c3
        REAL(realk), INTENT(in) :: alphaLoc
        REAL(realk), INTENT(out) :: cellProportion

        ! Local variables
        REAL(realk) :: mc1, mc2, mc3
        REAL(realk) :: vol
        REAL(realk) :: alphaMax
        REAL(realk) :: alphaStd
        REAL(realk) :: V1
        REAL(realk) :: baseArea, chamferedRectangleArea, triangularArea

        ! Solve the standart cases for vol
        ! Source: R. Scardovelli und S. Zaleski, „Analytical Relations Connecting Linear Interfaces and Volume Fractions in Rectangular Grids“,
        !         Journal of Computational Physics, Bd. 164, Nr. 1, S. 228–237, Okt. 2000, doi: 10.1006/jcph.2000.6567.
        mc1 = m1*c1
        mc2 = m2*c2
        mc3 = m3*c3

        IF ( mc1 < tol ) THEN
            IF ( mc2 < tol ) THEN
                ! One-dimensional case
                alphaMax = mc3
                alphaStd = MIN(alphaLoc, alphaMax - alphaLoc)

                IF ( alphaStd <= 0.0_realk ) THEN
                    IF ( alphaLoc >= alphaMax ) THEN
                        vol = c1 * c2 * c3
                    ELSE
                        vol = 0.0_realk
                    ENDIF 
                ELSE
                    vol = alphaStd * (c1*c2)
                ENDIF
            ELSE
                ! Two-dimensional cases
                alphaMax = mc2 + mc3
                alphaStd = MIN(alphaLoc, alphaMax - alphaLoc)

                IF ( alphaStd <= 0.0_realk ) THEN
                    IF ( alphaLoc >= alphaMax ) THEN
                        vol = c1 * c2 * c3
                    ELSE
                        vol = 0.0_realk
                    ENDIF
                ! Calculate vol dependent on mc2 and mc3
                ELSEIF ( alphaStd < mc2 ) THEN
                    baseArea = 1.0_realk/2.0_realk * alphaStd**2.0_realk / ( m2 * m3 )
                    vol = baseArea * c1
                ELSE
                    triangularArea = 1.0_realk/2.0_realk * c2**2.0_realk * m2 / m3
                    chamferedRectangleArea = c2 * alphaStd / m3 - triangularArea
                    vol = chamferedRectangleArea * c1
                ENDIF
            ENDIF
        ELSE
            ! Three-dimensional cases
            alphaMax = mc1 + mc2 + mc3
            alphaStd = MIN(alphaLoc, alphaMax - alphaLoc)

            V1 = mc1**2 * c1 / ( MAX(6.0_realk * m2 * m3, tol) )

            IF ( alphaStd <= 0.0_realk ) THEN
                IF ( alphaLoc >= alphaMax ) THEN
                    vol = c1 * c2 * c3
                ELSE
                    vol = 0.0_realk
                ENDIF
            ! Calculate vol dependent on mc1, mc2 and mc3
            ELSEIF ( alphaStd < mc1 ) THEN
                vol = alphaStd**3.0_realk / ( 6.0_realk * m1 * m2 * m3 )
            ELSE IF ( alphaStd < mc2 ) THEN
                vol = ( alphaStd * c1 * ( alphaStd - mc1 ) ) / ( 2.0_realk * m2 * m3 ) + V1
            ELSE IF ( alphaStd < MIN(mc1 + mc2, mc3) ) THEN
                vol = ( alphaStd**2.0_realk * ( 3.0_realk * ( mc1 + mc2 ) - alphaStd ) + mc1**2.0_realk * ( mc1 - 3.0_realk * alphaStd ) + mc2**2.0_realk * ( mc2 - 3.0_realk * alphaStd ) ) / ( 6.0_realk * m1 * m2 * m3 )
            ELSE IF ( alphaStd >= MIN(mc1 + mc2, mc3) .AND. mc3 <= mc1 + mc2 ) THEN
                vol = ( alphaStd**2.0_realk * ( 3.0_realk * ( mc1 + mc2 + mc3 ) - 2.0_realk * alphaStd ) &
                    + mc1**2.0_realk * ( mc1 - 3.0_realk * alphaStd ) &
                    + mc2**2.0_realk * ( mc2 - 3.0_realk * alphaStd ) &
                    + mc3**2.0_realk * ( mc3 - 3.0_realk * alphaStd ) ) / ( 6.0_realk * m1 * m2 * m3 )
            ELSE IF ( alphaStd >= MIN(mc1 + mc2, mc3) .AND. mc3 > mc1 + mc2 ) THEN
                vol = ( c1 * c2 * ( 2.0_realk * alphaStd - ( mc1 + mc2 ) ) ) / ( 2.0_realk * m3 )
            ENDIF
        ENDIF

        cellProportion = vol / ( c1 * c2 * c3 )

        IF ( alphaLoc > 0.5_realk * alphaMax .AND. alphaLoc < alphaMax ) THEN
            cellProportion = 1.0_realk - cellProportion
        ENDIF

    END SUBROUTINE comp_frac_std

    END MODULE multiphase_plic_mod