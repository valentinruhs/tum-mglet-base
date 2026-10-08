MODULE mph_bound_c_mod

    USE precision_mod, ONLY: realk, intk
    USE grids_mod, ONLY: get_mgdims
    USE field_mod, ONLY: field_t
    USE bound_mod, ONLY: bound_t
    USE err_mod, ONLY: errr

    IMPLICIT NONE(type, external)
    PRIVATE

    ! Bound operation on cell-centred VOF scalars (C)
    !   FIX, NOS, SLI: zero-gradient copy
    !   PAR:           injection of the coarse value from the
    !                  parent buffer into both ghost layers
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

    ! Bound operation on staggered triples (CS1..3, DS1..3, MS1..3,
    ! U, V, W). f1/f2/f3 are staggered in x/y/z. Only PAR is handled:
    !   normal component:      PAR face and outer face = coarse face
    !                          value (as bound_flow does for u)
    !   tangential components: both ghost layers = coarse value of the
    !                          first coarse cell outside (injection)
    TYPE, EXTENDS(bound_t) :: bound_stg_t
    CONTAINS
        PROCEDURE, NOPASS :: front => sfront
        PROCEDURE, NOPASS :: back => sfront
        PROCEDURE, NOPASS :: right => sright
        PROCEDURE, NOPASS :: left => sright
        PROCEDURE, NOPASS :: bottom => sbottom
        PROCEDURE, NOPASS :: top => sbottom
    END TYPE bound_stg_t
    TYPE(bound_stg_t) :: bound_stg

    PUBLIC :: bound_c, bound_stg

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
        REAL(realk), POINTER, CONTIGUOUS :: c(:, :, :), buf(:, :, :)

        ! Return early when no action is to be taken
        SELECT CASE (ctyp)
        CASE ("FIX", "NOS", "SLI", "PAR")
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
        CASE ("PAR")
            ! Both ghost layers lie inside the first coarse cell
            ! outside the grid: inject its value (bounded and
            ! consistent with the volume-weighted restriction "D")
            CALL f1%buffers%get_buffer(buf, igrid, iface)
            DO j = 1, jj
                DO k = 1, kk
                    c(k, j, i2) = buf(k, j, 1)
                    c(k, j, i1) = buf(k, j, 1)
                END DO
            END DO
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
        REAL(realk), POINTER, CONTIGUOUS :: c(:, :, :), buf(:, :, :)

        ! Return early when no action is to be taken
        SELECT CASE (ctyp)
        CASE ("FIX", "NOS", "SLI", "PAR")
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
        CASE ("PAR")
            CALL f1%buffers%get_buffer(buf, igrid, iface)
            DO i = 1, ii
                DO k = 1, kk
                    c(k, j2, i) = buf(k, i, 1)
                    c(k, j1, i) = buf(k, i, 1)
                END DO
            END DO
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
        REAL(realk), POINTER, CONTIGUOUS :: c(:, :, :), buf(:, :, :)

        ! Return early when no action is to be taken
        SELECT CASE (ctyp)
        CASE ("FIX", "NOS", "SLI", "PAR")
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
        CASE ("PAR")
            CALL f1%buffers%get_buffer(buf, igrid, iface)
            DO i = 1, ii
                DO j = 1, jj
                    c(k2, j, i) = buf(j, i, 1)
                    c(k1, j, i) = buf(j, i, 1)
                END DO
            END DO
        END SELECT

    END SUBROUTINE bbottom


    SUBROUTINE sfront(igrid, iface, ibocd, ctyp, f1, f2, f3, f4, timeph)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, iface, ibocd
        CHARACTER(len=*), INTENT(in) :: ctyp
        TYPE(field_t), INTENT(inout) :: f1
        TYPE(field_t), INTENT(inout), OPTIONAL :: f2, f3, f4
        REAL(realk), INTENT(in), OPTIONAL :: timeph

        ! Local variables
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: k, j, ig1, ig2, inf, ino
        REAL(realk), POINTER, CONTIGUOUS :: fn(:, :, :), ft1(:, :, :), ft2(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: bn(:, :, :), bt1(:, :, :), bt2(:, :, :)

        IF (ctyp /= "PAR") RETURN
        IF (.NOT. (PRESENT(f2) .AND. PRESENT(f3))) CALL errr(__FILE__, __LINE__)

        CALL get_mgdims(kk, jj, ii, igrid)

        ! ig1/ig2: ghost cell layers, inf: PAR face, ino: outer face
        SELECT CASE (iface)
        CASE (1)
            ! Front
            ig1 = 1
            ig2 = 2
            inf = 2
            ino = 1
        CASE (2)
            ! Back
            ig1 = ii
            ig2 = ii - 1
            inf = ii - 2
            ino = ii - 1
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        CALL f1%get_ptr(fn, igrid)
        CALL f2%get_ptr(ft1, igrid)
        CALL f3%get_ptr(ft2, igrid)
        CALL f1%buffers%get_buffer(bn, igrid, iface)
        CALL f2%buffers%get_buffer(bt1, igrid, iface)
        CALL f3%buffers%get_buffer(bt2, igrid, iface)

        DO j = 1, jj
            DO k = 1, kk
                fn(k, j, inf) = bn(k, j, 1)
                fn(k, j, ino) = bn(k, j, 1)
                ft1(k, j, ig2) = bt1(k, j, 1)
                ft1(k, j, ig1) = bt1(k, j, 1)
                ft2(k, j, ig2) = bt2(k, j, 1)
                ft2(k, j, ig1) = bt2(k, j, 1)
            END DO
        END DO

    END SUBROUTINE sfront


    SUBROUTINE sright(igrid, iface, ibocd, ctyp, f1, f2, f3, f4, timeph)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, iface, ibocd
        CHARACTER(len=*), INTENT(in) :: ctyp
        TYPE(field_t), INTENT(inout) :: f1
        TYPE(field_t), INTENT(inout), OPTIONAL :: f2, f3, f4
        REAL(realk), INTENT(in), OPTIONAL :: timeph

        ! Local variables
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: k, i, jg1, jg2, jnf, jno
        REAL(realk), POINTER, CONTIGUOUS :: fn(:, :, :), ft1(:, :, :), ft2(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: bn(:, :, :), bt1(:, :, :), bt2(:, :, :)

        IF (ctyp /= "PAR") RETURN
        IF (.NOT. (PRESENT(f2) .AND. PRESENT(f3))) CALL errr(__FILE__, __LINE__)

        CALL get_mgdims(kk, jj, ii, igrid)

        SELECT CASE (iface)
        CASE (3)
            ! Right
            jg1 = 1
            jg2 = 2
            jnf = 2
            jno = 1
        CASE (4)
            ! Left
            jg1 = jj
            jg2 = jj - 1
            jnf = jj - 2
            jno = jj - 1
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        ! Normal component is f2 (y-staggered)
        CALL f2%get_ptr(fn, igrid)
        CALL f1%get_ptr(ft1, igrid)
        CALL f3%get_ptr(ft2, igrid)
        CALL f2%buffers%get_buffer(bn, igrid, iface)
        CALL f1%buffers%get_buffer(bt1, igrid, iface)
        CALL f3%buffers%get_buffer(bt2, igrid, iface)

        DO i = 1, ii
            DO k = 1, kk
                fn(k, jnf, i) = bn(k, i, 1)
                fn(k, jno, i) = bn(k, i, 1)
                ft1(k, jg2, i) = bt1(k, i, 1)
                ft1(k, jg1, i) = bt1(k, i, 1)
                ft2(k, jg2, i) = bt2(k, i, 1)
                ft2(k, jg1, i) = bt2(k, i, 1)
            END DO
        END DO

    END SUBROUTINE sright


    SUBROUTINE sbottom(igrid, iface, ibocd, ctyp, f1, f2, f3, f4, timeph)
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: igrid, iface, ibocd
        CHARACTER(len=*), INTENT(in) :: ctyp
        TYPE(field_t), INTENT(inout) :: f1
        TYPE(field_t), INTENT(inout), OPTIONAL :: f2, f3, f4
        REAL(realk), INTENT(in), OPTIONAL :: timeph

        ! Local variables
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: j, i, kg1, kg2, knf, kno
        REAL(realk), POINTER, CONTIGUOUS :: fn(:, :, :), ft1(:, :, :), ft2(:, :, :)
        REAL(realk), POINTER, CONTIGUOUS :: bn(:, :, :), bt1(:, :, :), bt2(:, :, :)

        IF (ctyp /= "PAR") RETURN
        IF (.NOT. (PRESENT(f2) .AND. PRESENT(f3))) CALL errr(__FILE__, __LINE__)

        CALL get_mgdims(kk, jj, ii, igrid)

        SELECT CASE (iface)
        CASE (5)
            ! Bottom
            kg1 = 1
            kg2 = 2
            knf = 2
            kno = 1
        CASE (6)
            ! Top
            kg1 = kk
            kg2 = kk - 1
            knf = kk - 2
            kno = kk - 1
        CASE DEFAULT
            CALL errr(__FILE__, __LINE__)
        END SELECT

        ! Normal component is f3 (z-staggered)
        CALL f3%get_ptr(fn, igrid)
        CALL f1%get_ptr(ft1, igrid)
        CALL f2%get_ptr(ft2, igrid)
        CALL f3%buffers%get_buffer(bn, igrid, iface)
        CALL f1%buffers%get_buffer(bt1, igrid, iface)
        CALL f2%buffers%get_buffer(bt2, igrid, iface)

        DO i = 1, ii
            DO j = 1, jj
                fn(knf, j, i) = bn(j, i, 1)
                fn(kno, j, i) = bn(j, i, 1)
                ft1(kg2, j, i) = bt1(j, i, 1)
                ft1(kg1, j, i) = bt1(j, i, 1)
                ft2(kg2, j, i) = bt2(j, i, 1)
                ft2(kg1, j, i) = bt2(j, i, 1)
            END DO
        END DO

    END SUBROUTINE sbottom

END MODULE mph_bound_c_mod