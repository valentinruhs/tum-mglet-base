!====================================================================
!  Module: mph_pois_mod
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

MODULE mph_pois_mod

    USE precision_mod, ONLY: intk, realk
    USE grids_mod, ONLY: nmygrids, mygrids, get_mgdims
    USE fields_mod, ONLY: get_fieldptr
        
    IMPLICIT NONE(type, external)
    PRIVATE 

    PUBLIC :: init_mph_pois, finish_mph_pois, comp_mat_coeff_mph

CONTAINS

    SUBROUTINE init_mph_pois()

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        CONTINUE

    END SUBROUTINE init_mph_pois

    !================================================================

    SUBROUTINE finish_mph_pois()

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        CONTINUE

    END SUBROUTINE finish_mph_pois

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
        REAL(realk), POINTER, CONTIGUOUS :: dBa(:, :, :), dLe(:, :, :), dTo(:, :, :)


        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_mgdims(kk, jj, ii, igrid)

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

            CALL get_fieldptr(dBa, "DBA", igrid)
            CALL get_fieldptr(dLe, "DLE", igrid)
            CALL get_fieldptr(dTo, "DTO", igrid)

            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        ae(k,j,i) = 1.0_realk/(ddx(i)*dx(i)*dBa(k,j,i))
                        aw(k,j,i) = 1.0_realk/(ddx(i)*dx(i-1)*dBa(k,j,i-1))
                        an(k,j,i) = 1.0_realk/(ddy(j)*dy(j)*dLe(k,j,i))
                        as(k,j,i) = 1.0_realk/(ddy(j)*dy(j-1)*dLe(k,j-1,i))
                        at(k,j,i) = 1.0_realk/(ddz(k)*dz(k)*dTo(k,j,i))
                        ab(k,j,i) = 1.0_realk/(ddz(k)*dz(k-1)*dTo(k-1,j,i))
                    ENDDO
                ENDDO
            ENDDO

            DO i = 3, ii-2
                DO j = 3, jj-2
                    DO k = 3, kk-2
                        ap(k, j, i) = -(ae(k,j,i) + aw(k,j,i) + an(k,j,i) + &
                                        as(k,j,i) + at(k,j,i) + ab(k,j,i))
                    END DO
                END DO
            END DO
        ENDDO

    END SUBROUTINE comp_mat_coeff_mph

END MODULE mph_pois_mod