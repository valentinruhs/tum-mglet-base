!====================================================================
!  Module: mph_xfer_mod
!
!   Responsibilities:
!   - Provides the grid transfers of the multi-phase fields between
!     the levels of the locally refined grid:
!     - restriction (fine to coarse) of cell-centred, face and
!       staggered fields, the latter with PLIC half-cell values
!     - prolongation of cell-centred fields and staggered field
!       triples to the boundaries of the finer grids (parent,
!       bound, connect)
!     - exchange of ghost layers on the same level (connect)
!   - Registers the half-cell fields HLFFRCMI and HLFFRCPL
!
!   Author:      Valentin Ruhs
!   e-Mail:      valentin.ruhs@gmx.de
!   Created:     2026-09
!
!====================================================================

MODULE mph_xfer_mod

    USE precision_mod, ONLY: intk, realk
    USE grids_mod, ONLY: minlevel, maxlevel
    USE connect2_mod, ONLY: connect
    USE parent_mod, ONLY: parent
    USE ftoc_mod, ONLY: ftoc
    USE field_mod, ONLY: field_t
    USE fields_mod, ONLY: get_field, set_field
    USE bound_mod, ONLY: bound_t

    USE mph_utils_mod, ONLY: int2char
    USE mph_plic_mod, ONLY: comp_hlf_frc
    USE mph_bound_mod, ONLY: bound_c

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: init_mph_xfer, finish_mph_xfer, rstr, rstr_stg, &
        prlg, prlg_stg, cnct

CONTAINS

    SUBROUTINE init_mph_xfer()
    !----------------------------------------------------------------
    !   Short: Register half-cell fields for the
    !          staggered restriction.
    !----------------------------------------------------------------

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        CALL set_field("HLFFRCMI")
        CALL set_field("HLFFRCPL")

    END SUBROUTINE init_mph_xfer

    !================================================================

    SUBROUTINE finish_mph_xfer()
    !----------------------------------------------------------------
    !   Short: Finalize module (no operation).
    !----------------------------------------------------------------

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        CONTINUE

    END SUBROUTINE finish_mph_xfer

    !================================================================

    SUBROUTINE rstr(fldName, flag)
    !----------------------------------------------------------------
    !   Short: Restrict a cell-centred or face field from fine to
    !          coarse levels.
    !
    !   ftoc from the finest to the coarsest level with the given
    !   flag. For the cell flags (D, P) the coarse value is the
    !   volume-weighted mean of the eight fine cells, for the face
    !   flags (U, V, W) the area-weighted mean of the four fine faces
    !   that make up the coarse face (noib_restrict_mod).
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
    !   Short: Restrict a field on staggered grid q using PLIC
    !          half-cell values.
    !
    !   The ghost layers are exchanged first (connect), since the
    !   outermost coarse cells are restricted using CON values from
    !   the ghost cells. Then the half-cell values are computed
    !   (comp_hlf_frc) and ftoc is called with flag A, B or C for
    !   q = 1, 2, 3.
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
        ! Perform connect first, since the outermost cells are
        ! restricted using CON values from the ghost cells
        DO ilevel = maxlevel, minlevel, -1
            CALL connect(ilevel, 2, s1=fld_p)
        END DO
        CALL comp_hlf_frc(q, name)
        DO ilevel = maxlevel, minlevel, -1
            CALL ftoc(ilevel, fld_p%arr, fld_p%arr, flag)
        END DO

    END SUBROUTINE rstr_stg

    !================================================================

    SUBROUTINE prlg(fldName)
    !----------------------------------------------------------------
    !   Short: Prolongate a cell-centred field to fine-grid
    !          boundaries (parent, bound, connect).
    !
    !   From the coarsest to the finest level: parent fills the
    !   boundary buffers of the finer grids, bound_c writes them into
    !   the ghost layers at PAR faces (zero gradient at NOS/SLI),
    !   connect fills the CON ghost layers incl. corners.
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

    SUBROUTINE prlg_stg(fldName1, fldName2, fldName3, bnd)
    !----------------------------------------------------------------
    !   Short: Prolongate a staggered field triple to
    !          fine-grid boundaries.
    !
    !   parent fills the face buffers, bnd (bound_stg or bound_flow)
    !   writes them into PAR faces/ghost layers, connect fills the
    !   CON layers incl. corners.
    !----------------------------------------------------------------

        ! Subroutine arguments
        CHARACTER(len=*), INTENT(in) :: fldName1, fldName2, fldName3
        CLASS(bound_t), INTENT(inout) :: bnd

        ! Local variables
        TYPE(field_t), POINTER :: fld1_p, fld2_p, fld3_p
        INTEGER(intk) :: ilevel

        CALL get_field(fld1_p, fldName1)
        CALL get_field(fld2_p, fldName2)
        CALL get_field(fld3_p, fldName3)

        DO ilevel = minlevel, maxlevel
            CALL parent(ilevel, v1=fld1_p, v2=fld2_p, v3=fld3_p)
            CALL bnd%bound(ilevel, f1=fld1_p, f2=fld2_p, f3=fld3_p)
            CALL connect(ilevel, 2, v1=fld1_p, v2=fld2_p, v3=fld3_p, corners=.TRUE.)
        END DO

    END SUBROUTINE prlg_stg

    !================================================================

    SUBROUTINE cnct(fldName)
    !----------------------------------------------------------------
    !   Short: Exchange ghost layers of a field on all
    !          levels (connect).
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
