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
    USE err_mod, ONLY: errr
    
    IMPLICIT NONE
    PRIVATE 

    PUBLIC :: init_multiphase_plic, finish_multiphase_plic, track_interface, compute_normal_vector, compute_alpha, compute_cell_proportion

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

    SUBROUTINE track_interface(isInterface, kk, jj, ii, vff, tol)
    !----------------------------------------------------------------
    !   What it does:
    !   Identifies which of the cells in the domain contains a volume
    !   fraction of two fluids. These cells have to be taken into
    !   account when reconstructing interfaces. Therefore the
    !   variable containign this information is called isInterface.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        LOGICAL, INTENT(out) :: isInterface(kk, jj, ii)
        REAL(realk), INTENT(in) :: vff(kk, jj, ii)
        REAL(realk), INTENT(in) :: tol

        ! Local variables
        INTEGER(intk) :: k, j, i

        isInterface = .FALSE.

        ! Find cells with interface
        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    IF ( vff(k,j,i) > tol .AND. vff(k,j,i) < 1.0 - tol ) THEN
                        isInterface(k,j,i) = .TRUE.
                    END IF
                END DO
            END DO
        END DO

    END SUBROUTINE track_interface

    !================================================================

    SUBROUTINE compute_normal_vector(normx, normy, normz, kk, jj, ii, vff, ddx, ddy, ddz, tol)
    !----------------------------------------------------------------
    !   What it does:
    !   Computes the normal vector components normx, normy and normz
    !   of the gradient of the volume fraction function vff. The 
    !   components are normalized by the length to get the unit
    !   normal components. The gradient in each cell is calculated
    !   by taking into account its eight surrounding cells weighted
    !   with the three-dimensional sobel operator:
    !            / 1  2  1 \
    !   sobel =  | 2  4  2 |
    !            \ 1  2  1 /
    !   The gradient is approximated by a central-difference scheme:
    !   norm(.) = (upwind sum - downwind sum) / 2 * dd(.) * sobel sum
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(out) :: normx(kk, jj, ii), normy(kk, jj, ii), normz(kk, jj, ii)
        REAL(realk), INTENT(in) :: vff(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: tol

        ! Local variables
        INTEGER(intk) :: k, j, i, dDim1, dDim2
        INTEGER(intk), PARAMETER :: sobel(-1:1, -1:1) = reshape([1, 2, 1, 2, 4, 2, 1, 2, 1], [3,3])
        REAL(realk) :: normLength
        REAL(realk) :: sumStencilx, sumStencily, sumStencilz

        ! Loop over cells
        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2

                    sumStencilx = 0.0_realk
                    sumStencily = 0.0_realk
                    sumStencilz = 0.0_realk

                    ! Loop over sobel stencil
                    DO dDim2 = -1, 1
                        DO dDim1 = -1, 1
                            sumStencilx = sumStencilx + sobel(dDim1, dDim2) * ( vff(k+dDim1, j+dDim2, i+1) - vff(k+dDim1, j+dDim2, i-1) )
                            sumStencily = sumStencily + sobel(dDim1, dDim2) * ( vff(k+dDim1, j+1, i+dDim2) - vff(k+dDim1, j-1, i+dDim2) )
                            sumStencilz = sumStencilz + sobel(dDim1, dDim2) * ( vff(k+1, j+dDim1, i+dDim2) - vff(k-1, j+dDim1, i+dDim2) )
                        END DO
                    END DO

                    normx(k,j,i) = sumStencilx / ( 2 * sum(sobel) * ddx(i) )
                    normy(k,j,i) = sumStencily / ( 2 * sum(sobel) * ddy(j) )
                    normz(k,j,i) = sumStencilz / ( 2 * sum(sobel) * ddz(k) )

                    ! Calculate normal vector length
                    normLength = sqrt( normx(k,j,i)**2 + normy(k,j,i)**2 + normz(k,j,i)**2 )

                    ! Normalize with direction from high vff to low vff
                    IF ( normLength > tol ) THEN
                        normx(k,j,i) = - normx(k,j,i) / normLength
                        normy(k,j,i) = - normy(k,j,i) / normLength
                        normz(k,j,i) = - normz(k,j,i) / normLength
                    ELSE
                        normx(k,j,i) = 0.0
                        normy(k,j,i) = 0.0
                        normz(k,j,i) = 0.0
                    END IF

                END DO
            END DO
        END DO

    END SUBROUTINE compute_normal_vector

    !================================================================

    SUBROUTINE compute_alpha(alpha, kk, jj, ii, vff, isInterface, ddx, ddy, ddz, normx, normy, normz, tol)
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
        LOGICAL, INTENT(out) :: isInterface(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: normx(kk, jj, ii), normy(kk, jj, ii), normz(kk, jj, ii)
        REAL(realk), INTENT(in) :: tol

        ! Local variables
        INTEGER(intk) :: k, j, i
        REAL(realk) :: m1, m2, m3, c1, c2, c3
        REAL(realk) :: alphaStd(kk, jj, ii)
        REAL(realk) :: alphaMax(kk, jj, ii)


        ! Loop over cells
        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2

                    ! Only calculate interface for intersected cells
                    IF ( .NOT. isInterface(k,j,i) ) THEN 
                        CYCLE
                    END IF

                    ! 1. Assign norm(.) to m1, m2 and m3 and c1-c3 respectively
                    ! To enhance performance consider inlining
                    CALL get_corner_crossing_order(m1, m2, m3, c1, c2, c3, normx(k,j,i), normy(k,j,i), normz(k,j,i), ddx(i), ddy(j), ddz(k))

                    ! 2. Transform vff to a actual volume in bounds [0,0.5] * dV
                    ! 3. Solve the standart cases for alpha
                    ! Source: R. Scardovelli und S. Zaleski, „Analytical Relations Connecting Linear Interfaces and Volume Fractions in Rectangular Grids“,
                    !         Journal of Computational Physics, Bd. 164, Nr. 1, S. 228–237, Okt. 2000, doi: 10.1006/jcph.2000.6567.
                    ! To enhance performance consider inlining
                    CALL solve_alpha_standart_cases(m1, m2, m3, c1, c2, c3, alphaStd(k,j,i), alphaMax(k,j,i), vff(k,j,i), ddx(i), ddy(j), ddz(k), tol)

                    ! 4. If necessary, transform alpha back to volume bounds [0,1] * dV
                    ! If the volume fraction function has a value above 0.5 the "inverse problem" is solved. Therefore, the result is no longer 
                    ! alpha, but alphaMax - alpha. It can be seen as a rotation of the voxel. This is the inverse rotation (see solve_alpha_standart_cases)
                    alpha(k,j,i) = alphaStd(k,j,i)
                    IF ( vff(k,j,i) > 1.0/2.0 ) THEN
                        alpha(k,j,i) = alphaMax(k,j,i) - alpha(k,j,i)
                    END IF

                    ! 5. If necessary, transform alpha regarding to its negative normal vector components
                    ! If one of the normal vector components is negative, a mirrored case is solved. Therefore, the solution
                    ! has to be transformed back (see get_corner_crossing_order)
                    IF ( normx(k,j,i) < 0.0 ) THEN
                        alpha(k,j,i) = alpha(k,j,i) + ddx(i)*normx(k,j,i)
                    END IF

                    IF ( normy(k,j,i) < 0.0 ) THEN
                        alpha(k,j,i) = alpha(k,j,i) + ddy(j)*normy(k,j,i)
                    END IF

                    IF ( normz(k,j,i) < 0.0 ) THEN
                        alpha(k,j,i) = alpha(k,j,i) + ddz(k)*normz(k,j,i)
                    END IF

                END DO
            END DO
        END DO

    END SUBROUTINE compute_alpha

    !================================================================

    SUBROUTINE compute_cell_proportion(cellProportion, alpha, vff, ddx, ddy, ddz, normx, normy, normz, tol)
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
        REAL(realk), INTENT(in) :: tol

        ! Input variables
        REAL(realk) :: m1, m2, m3, c1, c2, c3
        REAL(realk) :: alphaLoc

        ! 1. Assign norm(.) to m1, m2 and m3 and c1-c3 respectively
        ! To enhance performance consider inlining
        CALL get_corner_crossing_order(m1, m2, m3, c1, c2, c3, normx, normy, normz, ddx, ddy, ddz)

        ! 2. If necessary, transform alpha regarding to its negative normal vector components
        alphaLoc = alpha
        IF ( normx < 0.0 ) THEN
            alphaLoc = alphaLoc - ddx*normx
        END IF

        IF ( normy < 0.0 ) THEN
            alphaLoc = alphaLoc - ddy*normy
        END IF

        IF ( normz < 0.0 ) THEN
            alphaLoc = alphaLoc - ddz*normz
        END IF

        ! 3. If necessary, transform alpha to its conjugate alphaMax - alpha
        ! 4. Solve the standart case for vol
        CALL solve_vol_standart_cases(m1, m2, m3, c1, c2, c3, alphaLoc, cellProportion, tol)

    END SUBROUTINE compute_cell_proportion

    !================================================================

    PURE SUBROUTINE get_corner_crossing_order(m1, m2, m3, c1, c2, c3, normx, normy, normz, ddx, ddy, ddz)
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
        END IF

        ! normx*ddx >= normz*ddz => normx is m3 and normz is m1 or m2
        ! normx*ddx <  normz*ddz => normz is m3 and normx is m1 (or m2)
        ! --OR--
        ! normy*ddy >= normz*ddz => normy is m3 and normz is m1 or m2
        ! normy*ddy <  normz*ddz => normz is m3 and normy is m1 (or m2)
        IF (scaledNorm2 >= scaledNorm3) THEN
            tmp = scaledNorm2; scaledNorm2 = scaledNorm3; scaledNorm3 = tmp
            tmpi = i2; i2 = i3; i3 = tmpi
        END IF

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
        END IF

        ! Assign new order to m1-m3 and c1-c3 respectively
        ! The absolute value of norm(.) is a mirror transform (see compute_alpha 5.).
        SELECT CASE (i1)
        CASE (1); m1 = abs(normx); c1 = ddx
        CASE (2); m1 = abs(normy); c1 = ddy
        CASE (3); m1 = abs(normz); c1 = ddz
        end SELECT

        SELECT CASE (i2)
        CASE (1); m2 = abs(normx); c2 = ddx
        CASE (2); m2 = abs(normy); c2 = ddy
        CASE (3); m2 = abs(normz); c2 = ddz
        end SELECT

        SELECT CASE (i3)
        CASE (1); m3 = abs(normx); c3 = ddx
        CASE (2); m3 = abs(normy); c3 = ddy
        CASE (3); m3 = abs(normz); c3 = ddz
        end SELECT

    END SUBROUTINE get_corner_crossing_order

    !================================================================

    PURE SUBROUTINE solve_alpha_standart_cases(m1, m2, m3, c1, c2, c3, alphaStd, alphaMax, vff, ddx, ddy, ddz, tol)
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
        REAL(realk), INTENT(in) :: tol
        
        ! Local variables
        REAL(realk) :: mc1, mc2, mc3
        REAL(realk) :: vol
        REAL(realk) :: baseArea, criticalBaseArea
        REAL(realk) :: V1, V2, V3
        REAL(realk) :: a0, a1, a2
        REAL(realk) :: qo, po
        REAL(realk) :: theta

        ! Transform vff to a actual volume in bounds [0,0.5] * dV
        ! Rotate voxel into standart configuration (see compute_alpha 4.)
        vol = min(vff, 1.0 - vff) * ddx * ddy * ddz
        
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
                criticalBaseArea = 1.0/2.0 * c2**2 * m2/m3
                
                IF ( baseArea <= criticalBaseArea ) THEN
                    ! Here both interception lines of the interface with the coordinate axis are within the cell => triangular prism
                    alphaStd = sqrt(2.0 * baseArea * m2 * m3)
                ELSE
                    ! Here one interception line (with the c2 axis) is outside the cell => chamfered rectangle prism
                    alphaStd = (m3) / (c2) * baseArea + (mc2) / 2.0
                END IF
            END IF
        ELSE
            ! Three-dimensional cases
            alphaMax = mc1 + mc2 + mc3

            ! Define interval boundaries V1, V2, V3
            V1 = mc1**2 * c1 / ( max(6.0 * m2 * m3, tol) )
            V2 = V1 + c1 * c2 * ( mc2 - mc1 ) / ( 2.0 * m3 )
            IF ( mc3 < mc1 + mc2 ) THEN
                V3 = ( mc3**2 * ( 3.0 * ( mc1 + mc2 ) - mc3 ) + mc1**2 * ( mc1 - 3.0 * mc3 ) + mc2**2 * ( mc2 - 3.0 * mc3 ) ) / ( 6.0 * m1 * m2 * m3 )
            ELSE
                V3 = c1 * c2 * ( mc1 + mc2 ) / ( 2.0 * m3 )
            END IF
            
            ! Calculate alphaStd dependent on V1, V2 and V3
            IF ( vol < V1 ) THEN
                alphaStd = ( 6.0 * m1 * m2 * m3 * vol )**( 1.0/3.0 )
            ELSE IF ( vol < V2 ) THEN
                alphaStd = 1.0/2.0 * ( mc1 + sqrt(mc1**2 + 8.0 * m2 * m3 * (vol - V1) / c1) )
            ELSE IF ( vol < V3 ) THEN
                a2 = - 3.0 * ( mc1 + mc2 )
                a1 = 3.0 * ( mc1**2 + mc2**2 )
                a0 = - (mc1**3 + mc2**3) + 6.0 * m1 * m2 * m3 * vol
                po = a1 / 3.0 - a2**2 / 9.0
                qo = ( a1 * a2 - 3.0 * a0 ) / 6.0 - a2**3 / 27.0
                
                ! ! Debug
                ! IF ( po**3 + qo**2 > 0 ) THEN
                !     WRITE(*, *) "No real roots for alphaStd. po^3 + qo^2 = ", po**3 + qo**2, " po = ", po
                !     CALL errr(__FILE__, __LINE__)
                ! END IF
                
                theta = acos(qo / sqrt((-po)**3)) / 3.0
                alphaStd = sqrt(-po) * ( sqrt(3.0) * sin(theta) - cos(theta) ) - a2 / 3.0
            ELSE IF ( vol >= V3 .AND. mc3 <= mc1 + mc2 ) THEN
                a2 = - 3.0/2.0 * ( mc1 + mc2 + mc3 )
                a1 = 3.0/2.0 * ( mc1**2 + mc2**2 + mc3**2 )
                a0 = - 1.0/2.0 * ( mc1**3 + mc2**3 + mc3**3 ) + 3.0 * m1 * m2 * m3 * vol
                po = a1 / 3.0 - a2**2 / 9.0
                qo = ( a1 * a2 - 3.0 * a0 ) / 6.0 - a2**3 / 27.0
                
                ! ! Debug
                ! IF ( po**3 + qo**2 > 0 ) THEN
                !     WRITE(*, *) "No real roots for alphaStd. po^3 + qo^2 = ", po**3 + qo**2, " po = ", po
                !     CALL errr(__FILE__, __LINE__)
                ! END IF
                
                theta = acos(qo / sqrt((-po)**3)) / 3.0
                alphaStd = sqrt(-po) * ( sqrt(3.0) * sin(theta) - cos(theta) ) - a2 / 3.0
            ELSE IF ( vol >= V3 .AND. mc3 > mc1 + mc2 ) THEN
                alphaStd = m3 * vol / ( c1 * c2 ) + ( mc1 + mc2 ) / 2.0
            END IF
        END IF

    END SUBROUTINE solve_alpha_standart_cases

    !================================================================

    PURE SUBROUTINE solve_vol_standart_cases(m1, m2, m3, c1, c2, c3, alphaLoc, vff, tol)
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
        REAL(realk), INTENT(in) :: tol
        
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
                alphaStd = min(alphaLoc, alphaMax - alphaLoc)

                IF ( alphaStd <= 0.0 ) THEN
                    IF ( alphaLoc >= alphaMax ) THEN
                        vol = c1 * c2 * c3
                    ELSE
                        vol = 0.0
                    END IF 
                ELSE
                    vol = alphaStd * (c1*c2)
                END IF
            ELSE
                ! Two-dimensional cases
                alphaMax = mc2 + mc3
                alphaStd = min(alphaLoc, alphaMax - alphaLoc)

                IF ( alphaStd <= 0.0 ) THEN
                    IF ( alphaLoc >= alphaMax ) THEN
                        vol = c1 * c2 * c3
                    ELSE
                        vol = 0.0
                    END IF
                ! Calculate vol dependent on mc2 and mc3
                ELSEIF ( alphaStd < mc2 ) THEN
                    baseArea = 1.0/2.0 * alphaStd**2 / ( m2 * m3 )
                    vol = baseArea * c1
                ELSE
                    triangularArea = 1.0/2.0 * c2**2 * m2 / m3
                    chamferedRectangleArea = c2 * alphaStd / m3 - triangularArea
                    vol = chamferedRectangleArea * c1
                END IF
            END IF
        ELSE
            ! Three-dimensional cases
            alphaMax = mc1 + mc2 + mc3
            alphaStd = min(alphaLoc, alphaMax - alphaLoc)

            V1 = mc1**2 * c1 / ( max(6.0 * m2 * m3, tol) )

            IF ( alphaStd <= 0.0 ) THEN
                IF ( alphaLoc >= alphaMax ) THEN
                    vol = c1 * c2 * c3
                ELSE
                    vol = 0.0
                END IF
            ! Calculate vol dependent on mc1, mc2 and mc3
            ELSEIF ( alphaStd < mc1 ) THEN
                vol = alphaStd**3 / ( 6.0 * m1 * m2 * m3 )
            ELSE IF ( alphaStd < mc2 ) THEN
                vol = ( alphaStd * c1 * ( alphaStd - mc1 ) ) / ( 2.0 * m2 * m3 ) + V1
            ELSE IF ( alphaStd < min(mc1 + mc2, mc3) ) THEN
                vol = ( alphaStd**2 * ( 3.0 * ( mc1 + mc2 ) - alphaStd ) + mc1**2 * ( mc1 - 3.0 * alphaStd ) + mc2**2 * ( mc2 - 3.0 * alphaStd ) ) / ( 6.0 * m1 * m2 * m3 )
            ELSE IF ( alphaStd >= min(mc1 + mc2, mc3) .AND. mc3 <= mc1 + mc2 ) THEN
                vol = ( alphaStd**2 * ( 3.0 * ( mc1 + mc2 + mc3 ) - 2.0 * alphaStd ) &
                      + mc1**2 * ( mc1 - 3.0 * alphaStd ) &
                      + mc2**2 * ( mc2 - 3.0 * alphaStd ) &
                      + mc3**2 * ( mc3 - 3.0 * alphaStd ) ) / ( 6.0 * m1 * m2 * m3 )
            ELSE IF ( alphaStd >= min(mc1 + mc2, mc3) .AND. mc3 > mc1 + mc2 ) THEN
                vol = ( c1 * c2 * ( 2.0 * alphaStd - ( mc1 + mc2 ) ) ) / ( 2.0 * m3 )
            END IF
        END IF

        cellProportion = vol / ( c1 * c2 * c3 )

        IF ( alphaLoc > 0.5 * alphaMax .AND. alphaLoc < alphaMax ) THEN
            cellProportion = 1 - cellProportion
        END IF

    END SUBROUTINE solve_vol_standart_cases

END MODULE multiphase_plic_mod