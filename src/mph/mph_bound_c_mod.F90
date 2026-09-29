MODULE mph_bound_c_mod

    USE precision_mod, ONLY: realk, intk
    USE grids_mod, ONLY: get_mgdims
    USE field_mod, ONLY: field_t
    USE bound_mod, ONLY: bound_t
    USE err_mod, ONLY: errr

    IMPLICIT NONE(type, external)
    PRIVATE

    ! Bound operation 'T' operate on C
    TYPE, EXTENDS(bound_t) :: bound_c_t
    CONTAINS
        PROCEDURE, NOPASS :: front => bfront
        PROCEDURE, NOPASS :: back => bfront
        PROCEDURE, NOPASS :: right => bright
        PROCEDURE, NOPASS :: left => bright
        PROCEDURE, NOPASS :: bottom => bbottom
        PROCEDURE, NOPASS :: top => bbottom
    END TYPE bound_c_t
    TYPE(bound_c_t) :: bound_c

    PUBLIC :: bound_c

CONTAINS
    SUBROUTINE bfront(igrid, iface, ibocd, ctyp, f1, f2, f3, f4, timeph)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, iface, ibocd
        CHARACTER(len=*), INTENT(in) :: ctyp
        TYPE(field_t), INTENT(inout) :: f1
        TYPE(field_t), INTENT(inout), OPTIONAL :: f2, f3, f4
        REAL(realk), INTENT(in), OPTIONAL :: timeph

        ! Local variables
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: k, j, i1, i2, i3, i4
        REAL(realk), POINTER, CONTIGUOUS :: c(:, :, :)

        ! Return early when no action is to be taken
        SELECT CASE (ctyp)
        CASE ("FIX", "NOS", "SLI")
            CONTINUE
        CASE DEFAULT
            RETURN
        END SELECT

        ! Fetch pointers
        CALL f1%get_ptr(c, igrid)
        CALL get_mgdims(kk, jj, ii, igrid)

        SELECT CASE (iface)
        CASE (1)
            ! Front
            i1 = 1
            i2 = 2
            i3 = 3
            i4 = 4
        CASE (2)
            ! Back
            i1 = ii
            i2 = ii - 1
            i3 = ii - 2
            i4 = ii - 3
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        SELECT CASE ( ctyp )
        CASE ("FIX", "NOS", "SLI")
            DO j = 1, jj
                DO k = 1, kk
                    c(k, j, i2) = c(k, j, i3)
                    c(k, j, i1) = c(k, j, i4)
                END DO
            END DO
        ! extendable
        END SELECT
    
    END SUBROUTINE bfront


    SUBROUTINE bright(igrid, iface, ibocd, ctyp, f1, f2, f3, f4, timeph)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, iface, ibocd
        CHARACTER(len=*), INTENT(in) :: ctyp
        TYPE(field_t), INTENT(inout) :: f1
        TYPE(field_t), INTENT(inout), OPTIONAL :: f2, f3, f4
        REAL(realk), INTENT(in), OPTIONAL :: timeph

        ! Local variables
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: k, i, j1, j2, j3, j4
        REAL(realk), POINTER, CONTIGUOUS :: c(:, :, :)

        ! Return early when no action is to be taken
        SELECT CASE (ctyp)
        CASE ("FIX", "NOS", "SLI")
            CONTINUE
        CASE DEFAULT
            RETURN
        END SELECT

        ! Fetch pointers
        CALL f1%get_ptr(c, igrid)
        CALL get_mgdims(kk, jj, ii, igrid)

        SELECT CASE (iface)
        CASE (3)
            ! Right
            j1 = 1
            j2 = 2
            j3 = 3
            j4 = 4
        CASE (4)
            ! Left
            j1 = jj
            j2 = jj - 1
            j3 = jj - 2
            j4 = jj - 3
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        SELECT CASE ( ctyp )
        CASE ("FIX", "NOS", "SLI")
            DO i = 1, ii
                DO k = 1, kk
                    c(k, j2, i) = c(k, j3, i)
                    c(k, j1, i) = c(k, j4, i)
                END DO
            END DO
        ! extendable
        END SELECT

    END SUBROUTINE bright


    SUBROUTINE bbottom(igrid, iface, ibocd, ctyp, f1, f2, f3, f4, timeph)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, iface, ibocd
        CHARACTER(len=*), INTENT(in) :: ctyp
        TYPE(field_t), INTENT(inout) :: f1
        TYPE(field_t), INTENT(inout), OPTIONAL :: f2, f3, f4
        REAL(realk), INTENT(in), OPTIONAL :: timeph

        ! Local variables
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: j, i, k1, k2, k3, k4
        REAL(realk), POINTER, CONTIGUOUS :: c(:, :, :)

        ! Return early when no action is to be taken
        SELECT CASE (ctyp)
        CASE ("FIX", "NOS", "SLI")
            CONTINUE
        CASE DEFAULT
            RETURN
        END SELECT

        ! Fetch pointers
        CALL f1%get_ptr(c, igrid)
        CALL get_mgdims(kk, jj, ii, igrid)

        SELECT CASE (iface)
        CASE (5)
            ! Bottom
            k1 = 1
            k2 = 2
            k3 = 3
            k4 = 4
        CASE (6)
            ! Top
            k1 = kk
            k2 = kk - 1
            k3 = kk - 2
            k4 = kk - 3
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        SELECT CASE ( ctyp )
        CASE ("FIX", "NOS", "SLI")
            DO i = 1, ii
                DO j = 1, jj
                    c(k2, j, i) = c(k3, j, i)
                    c(k1, j, i) = c(k4, j, i)
                END DO
            END DO
        ! extendable
        END SELECT

    END SUBROUTINE bbottom

END MODULE mph_bound_c_mod
