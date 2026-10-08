!====================================================================
!  Module: mph_xfer_mod
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

MODULE mph_xfer_mod

    USE precision_mod, ONLY: intk, realk
    USE grids_mod, ONLY: minlevel, maxlevel
    USE connect2_mod, ONLY: connect
    USE parent_mod, ONLY: parent
    USE ftoc_mod, ONLY: ftoc
    USE field_mod, ONLY: field_t
    USE fields_mod, ONLY: get_field
    USE mph_bound_c_mod, ONLY: bound_c, bound_stg
    USE mph_utils_mod, ONLY: int2char

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: init_mph_xfer, finish_mph_xfer, rstr, rstr_stg, &
        prlg, prlg_stg, cnct

CONTAINS

    SUBROUTINE init_mph_xfer()

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        CONTINUE

    END SUBROUTINE init_mph_xfer

    !================================================================

    SUBROUTINE finish_mph_xfer()

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        CONTINUE

    END SUBROUTINE finish_mph_xfer
    
    !================================================================

    SUBROUTINE rstr(fldName, flag)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        CHARACTER(len=*), INTENT(in) :: fldName
        CHARACTER(len=*), INTENT(in) :: flag

        ! Local variables
        TYPE(field_t), POINTER :: fld_p
        INTEGER(intk) :: ilevel

        CALL get_field(fld_p, fldName)
        DO ilevel = maxlevel, minlevel, -1
            CALL ftoc(ilevel, fld_p%arr, fld_p%arr, flag)
        END DO

    END SUBROUTINE rstr

    !================================================================

    SUBROUTINE rstr_stg(q, fldName)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: q
        CHARACTER(len=*), INTENT(in) :: fldName

        ! Local variables
        CHARACTER(len=3) :: name
        CHARACTER(len=1) :: flag
        TYPE(field_t), POINTER :: fld_p
        INTEGER(intk) :: ilevel

        name = fldName//int2Char(q)

        SELECT CASE ( q )
        CASE ( 1 ); flag = "A"
        CASE ( 2 ); flag = "B"
        CASE ( 3 ); flag = "C"
        END SELECT

        CALL get_field(fld_p, name)
        DO ilevel = maxlevel, minlevel, -1
            CALL connect(ilevel, 2, s1=fld_p)
            CALL ftoc(ilevel, fld_p%arr, fld_p%arr, flag)
        END DO

    END SUBROUTINE rstr_stg

    !================================================================

    SUBROUTINE prlg(fldName)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        CHARACTER(len=*), INTENT(in) :: fldName

        ! Local variables
        TYPE(field_t), POINTER :: fld_p
        INTEGER(intk) :: ilevel

        CALL get_field(fld_p, fldName)

        DO ilevel = minlevel, maxlevel
            CALL parent(ilevel, s1=fld_p)
            CALL bound_c%bound(ilevel, f1=fld_p)
            CALL connect(ilevel, 2, s1=fld_p, corners=.TRUE.)
        END DO

    END SUBROUTINE prlg

    !================================================================

    SUBROUTINE prlg_stg(fldName1, fldName2, fldName3)
    !----------------------------------------------------------------
    !   What it does:
    !   Prolongation of a staggered triple to the boundaries of the
    !   finer grids: parent fills the face buffers, bound_stg writes
    !   them into PAR faces/ghost layers, connect fills CON layers.
    !----------------------------------------------------------------

        ! Subroutine arguments
        CHARACTER(len=*), INTENT(in) :: fldName1, fldName2, fldName3

        ! Local variables
        TYPE(field_t), POINTER :: fld1_p, fld2_p, fld3_p
        INTEGER(intk) :: ilevel

        CALL get_field(fld1_p, fldName1)
        CALL get_field(fld2_p, fldName2)
        CALL get_field(fld3_p, fldName3)

        DO ilevel = minlevel, maxlevel
            CALL parent(ilevel, v1=fld1_p, v2=fld2_p, v3=fld3_p)
            CALL bound_stg%bound(ilevel, f1=fld1_p, f2=fld2_p, f3=fld3_p)
            CALL connect(ilevel, 2, v1=fld1_p, v2=fld2_p, v3=fld3_p, corners=.TRUE.)
        END DO

    END SUBROUTINE prlg_stg

    !================================================================

    SUBROUTINE cnct(fldName)
    !----------------------------------------------------------------
    !   What it does:
    !   Same-level exchange (connect) of one field on all levels.
    !----------------------------------------------------------------

        ! Subroutine arguments
        CHARACTER(len=*), INTENT(in) :: fldName

        ! Local variables
        TYPE(field_t), POINTER :: fld_p
        INTEGER(intk) :: ilevel

        CALL get_field(fld_p, fldName)
        DO ilevel = minlevel, maxlevel
            CALL connect(ilevel, 2, s1=fld_p, corners=.TRUE.)
        END DO

    END SUBROUTINE cnct

END MODULE mph_xfer_mod
