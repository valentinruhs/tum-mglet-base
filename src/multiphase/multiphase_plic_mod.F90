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

    PUBLIC :: track_interface, compute_normal_vector, compute_alpha

CONTAINS

    SUBROUTINE init_multiphase_plic()

        continue

    END SUBROUTINE init_multiphase_plic

    !================================================================

    SUBROUTINE finish_multiphase_plic()

        continue

    END SUBROUTINE finish_multiphase_plic

    !================================================================

    SUBROUTINE track_interface(is_interface, kk, jj, ii, c, tol)
    !----------------------------------------------------------------
    !   What it does:
    !   Identifies which of the cells in the domain contains a volume
    !   fraction of two fluids. These cells have to be taken into
    !   account when reconstructing interfaces. Therefore the
    !   variable containign this information is called is_interface.
    !----------------------------------------------------------------

        ! Subroutine arguments
        LOGICAL, INTENT(out) :: is_interface(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: c(kk, jj, ii)
        REAL(realk), INTENT(in) :: tol

        ! Local variables
        INTEGER(intk) :: k, j, i

        is_interface = .FALSE.

        ! Find cells with interface
        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    IF ( c(k,j,i) > tol .AND. c(k,j,i) < 1.0 - tol ) THEN
                        is_interface(k,j,i) = .TRUE.
                    END IF
                END DO
            END DO
        END DO

    END SUBROUTINE track_interface

    !================================================================

    SUBROUTINE compute_normal_vector(normx, normy, normz, kk, jj, ii, c, ddx, ddy, ddz, tol)
    !----------------------------------------------------------------
    !   What it does:
    !   Computes the normal vector components normx, normy and normz
    !   of the gradient of the Color-Function c. The components are 
    !   normalized by the length to get the unit normal components. 
    !   The gradient in each cell is calculated by taking into 
    !   account its eight surrounding cells weighted with the 
    !   three-dimensional sobel operator:
    !           1  2  1
    !   sobel = 2  4  2
    !           1  2  1
    !   The gradient is approximated by a central-difference scheme:
    !   norm(.) = (upwind sum - downwind sum) / 2 * dd(.) * sobel sum
    !----------------------------------------------------------------

        ! Subroutine arguments
        REAL(realk), INTENT(out) :: normx(kk, jj, ii), normy(kk, jj, ii), normz(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: c(kk, jj, ii)
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
                            sumStencilx = sumStencilx + sobel(dDim1, dDim2) * ( c(k+dDim1, j+dDim2, i+1) - c(k+dDim1, j+dDim2, i-1) )
                            sumStencily = sumStencily + sobel(dDim1, dDim2) * ( c(k+dDim1, j+1, i+dDim2) - c(k+dDim1, j-1, i+dDim2) )
                            sumStencilz = sumStencilz + sobel(dDim1, dDim2) * ( c(k+1, j+dDim1, i+dDim2) - c(k-1, j+dDim1, i+dDim2) )
                        END DO
                    END DO

                    normx(k,j,i) = sumStencilx / ( 2 * sum(sobel) * ddx(i) )
                    normy(k,j,i) = sumStencily / ( 2 * sum(sobel) * ddy(j) )
                    normz(k,j,i) = sumStencilz / ( 2 * sum(sobel) * ddz(k) )

                    ! Calculate normal vector length
                    normLength = sqrt( normx(k,j,i)**2 + normy(k,j,i)**2 + normz(k,j,i)**2 )

                    ! Normalize with direction from high c to low c
                    IF ( normLength > tol ) THEN
                        normx(k,j,i) = - normx(k,j,i) / normLength
                        normy(k,j,i) = - normy(k,j,i) / normLength
                        normz(k,j,i) = - normz(k,j,i) / normLength     
                    END IF

                END DO
            END DO
        END DO

    END SUBROUTINE compute_normal_vector

    !================================================================

    SUBROUTINE compute_alpha(alpha, kk, jj, ii, c, is_interface, ddx, ddy, ddz, normx, normy, normz, tol)
    !----------------------------------------------------------------
    !   What it does:
    !   This subroutine calculates the alpha value for PLIC. The 
    !   alpha value describes the distance of the interface in a cell
    !   from a defined reference (left bottom front corner).
    !   1. Assign norm(.) to m1, m2 and m3 and c1-c3 respectively
    !   2. Transform c to a actual volume in bounds [0,0.5] * dV
    !   3. Solve the standart case for alpha
    !   4. If necessary, transform alpha back to volume bounds 
    !      [0,1] * dV
    !   5. If necessary, transform alpha regarding to its negative
    !      normal vector components
    !----------------------------------------------------------------

        ! Subroutine arguments
        REAL(realk), INTENT(out) :: alpha(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: c(kk, jj, ii)
        LOGICAL, INTENT(out) :: is_interface(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: normx(kk, jj, ii), normy(kk, jj, ii), normz(kk, jj, ii)
        REAL(realk), INTENT(in) :: tol

        ! Local variables
        INTEGER(intk) :: k, j, i
        REAL(realk) :: m1, m2, m3, c1, c2, c3, mc1, mc2, mc3
        REAL(realk) :: alpha_max(kk, jj, ii)
        REAL(realk) :: vol
        REAL(realk) :: base_area, critical_base_area
        REAL(realk) :: V1, V2, V3
        REAL(realk) :: a0, a1, a2
        REAL(realk) :: qo, po
        REAL(realk) :: theta


        ! Loop over cells
        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2

                    ! Only calculate interface for intersected cells
                    IF ( .NOT. is_interface(k,j,i) ) THEN 
                        CYCLE
                    END IF

                    ! 1. Assign norm(.) to m1, m2 and m3 and c1-c3 respectively
                    ! To enhance performance consider inlining
                    CALL get_corner_crossing_order(m1, m2, m3, c1, c2, c3, normx(k,j,i), normy(k,j,i), normz(k,j,i), ddx(i), ddy(j), ddz(k))

                    ! 2. Transform c to a actual volume in bounds [0,0.5] * dV
                    vol = min(c(k,j,i), 1 - c(k,j,i)) * ddx(i) * ddy(j) * ddz(k)

                    ! 3. Solve the standart cases for alpha
                    ! Source: R. Scardovelli und S. Zaleski, „Analytical Relations Connecting Linear Interfaces and Volume Fractions in Rectangular Grids“,
                    !         Journal of Computational Physics, Bd. 164, Nr. 1, S. 228–237, Okt. 2000, doi: 10.1006/jcph.2000.6567.
                    mc1 = m1*c1
                    mc2 = m2*c2
                    mc3 = m3*c3
                    
                    IF ( mc1 < tol ) THEN
                        IF ( mc2 < tol ) THEN
                            ! One-dimensional case
                            alpha_max(k,j,i) = mc3
                            alpha(k,j,i) = vol / (c1*c2)
                        ELSE
                            ! Two-dimensional cases
                            alpha_max(k,j,i) = mc2 + mc3

                            ! actual base area
                            base_area = vol / c1

                            ! When the critical base area is exceeded the volume shape transforms to a chamfered rectangle prism instead of triangular prism
                            critical_base_area = 1.0/2.0 * c2**2 * m2/m3

                            IF ( base_area <= critical_base_area ) THEN
                                ! Here both interception points of the interface are within the cell => triangular prism
                                alpha(k,j,i) = sqrt(2 * base_area * m2 * m3)
                            ELSE
                                ! Here one interception point (on the c2 axis) is outside the cell => chamfered rectangle prism
                                alpha(k,j,i) = (m3) / (c2) * base_area + (mc2) / 2 
                            END IF
                        END IF
                    ELSE
                        ! Three-dimensional cases
                        alpha_max(k,j,i) = mc1 + mc2 + mc3
                        ! Define interval boundaries V1, V2, V3
                        V1 = mc1**2 * c1 / ( max(6*m2*m3, tol) )
                        V2 = V1 + c1 * c2 * ( mc2 - mc1 ) / ( 2*m3 )
                        IF ( mc3 < mc1 + mc2 ) THEN
                            V3 = ( mc3**2 * ( 3 * ( mc1 + mc2 ) - mc3 ) + mc1**2 * ( mc1 - 3 * mc3 ) + mc2**2 * ( mc2 - 3 * mc3 ) ) / ( 6 * m1 * m2 * m3 )
                        ELSE
                            V3 = c1 * c2 * ( mc1 + mc2 ) / ( 2 * m3 )
                        END IF

                        ! Calculate alpha dependent on V1, V2 and V3
                        IF ( vol < V1 ) THEN
                            alpha(k,j,i) = ( 6 * m1 * m2 * m3 * vol )**( 1.0/3.0 )
                        ELSE IF ( vol < V2 ) THEN
                            alpha(k,j,i) = 1.0/2.0 * ( mc1 + sqrt(mc1**2 + 8 * m2 * m3 * (vol - V1) / c1) )
                        ELSE IF ( vol < V3 ) THEN
                            a2 = - 3 * ( mc1 + mc2 )
                            a1 = 3 * ( mc1**2 + mc2**2 )
                            a0 = - (mc1**3 + mc2**3) + 6 * m1 * m2 * m3 * vol
                            po = a1 / 3 - a2**2 / 9
                            qo = ( a1 * a2 - 3 * a0 ) / 6 - a2**3 / 27

                            ! ! Debug
                            ! IF ( po**3 + qo**2 > 0 .OR. po > 0 ) THEN
                            !     WRITE(*, *) "No real roots for alpha. po^3 + qo^2 = ", po**3 + qo**2, " po = ", po
                            !     CALL errr(__FILE__, __LINE__)
                            ! END IF

                            theta = acos(qo / sqrt((-po)**3)) / 3
                            alpha(k,j,i) = sqrt(-po) * ( sqrt(3.0) * sin(theta) - cos(theta) ) - a2 / 3
                        ELSE IF ( vol >= V3 .AND. mc3 <= mc1 + mc2 ) THEN
                            a2 = - 3.0/2.0
                            a1 = 3.0/2.0 * ( mc1**2 + mc2**2 + mc3**2 )
                            a0 = - 1.0/2.0 * (mc1**3 + mc2**3) + 3 * m1 * m2 * m3 * vol
                            po = a1 / 3 - a2**2 / 9
                            qo = ( a1 * a2 - 3 * a0 ) / 6 - a2**3 / 27

                            ! ! Debug
                            ! IF ( po**3 + qo**2 > 0 .OR. po > 0 ) THEN
                            !     WRITE(*, *) "No real roots for alpha. po^3 + qo^2 = ", po**3 + qo**2, " po = ", po
                            !     CALL errr(__FILE__, __LINE__)
                            ! END IF

                            theta = acos(qo / sqrt((-po)**3)) / 3
                            alpha(k,j,i) = sqrt(-po) * ( sqrt(3.0) * sin(theta) - cos(theta) ) - a2 / 3
                        ELSE IF ( vol >= V3 .AND. mc3 > mc1 + mc2 ) THEN
                            alpha(k,j,i) = m3 * vol / ( c1 * c2 ) + ( mc1 + mc2 ) / 2
                        END IF
                    END IF

                    ! 4. If necessary, transform alpha back to volume bounds [0,1] * dV
                    ! If the Color-Function has a value above 0.5 the "inverse problem" is solved. Therefore, the result is no longer 
                    ! alpha, but alpha_max - alpha
                    IF ( c(k,j,i) > 1.0/2.0 ) THEN
                        alpha(k,j,i) = alpha_max(k,j,i) - alpha(k,j,i)
                    END IF

                    ! 5. If necessary, transform alpha regarding to its negative normal vector components
                    ! If one of the normal vector components is negative, a mirrored case is solved. Therefore, the solution
                    ! has to be transformed back
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

END MODULE multiphase_plic_mod