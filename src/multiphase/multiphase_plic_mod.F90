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
    
    IMPLICIT NONE

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

    SUBROUTINE compute_alpha(alpha, kk, jj, ii, c, is_interface, ddx, ddy, ddz, normx, normy, normz)

        ! Subroutine arguments
        REAL(realk), INTENT(out) :: alpha(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: c(kk, jj, ii)
        LOGICAL, INTENT(out) :: is_interface(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(out) :: normx(kk, jj, ii), normy(kk, jj, ii), normz(kk, jj, ii)

        ! Local variables
        

        continue

    END SUBROUTINE compute_alpha

END MODULE multiphase_plic_mod