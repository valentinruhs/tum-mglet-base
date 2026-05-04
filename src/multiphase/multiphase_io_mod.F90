!====================================================================
!  Module: multiphase_io_mod
!
!  Responsibilities:
!     - Reads initial volume fraction field vff
!
!  Author:      Valentin Ruhs
!  Created:     2026-02
!  Last update: 2026-02
!
!====================================================================

MODULE multiphase_io_mod

    USE precision_mod, ONLY: intk
    USE field_mod, ONLY: field_t
    USE grids_mod, ONLY: nmygrids, mygrids
    USE field_mod, ONLY: field_t
    USE grids_mod, ONLY: get_mgdims, get_mgbasb
    USE precision_mod, ONLY: intk, realk
    USE fields_mod, ONLY: get_field
    USE multiphasecore_mod, ONLY: test_multiphase

    IMPLICIT NONE
    PRIVATE 

    PUBLIC :: init_multiphase_io, finish_multiphase_io, read_vff, update_velocity

CONTAINS

    SUBROUTINE init_multiphase_io()

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        continue
    END SUBROUTINE init_multiphase_io

    !================================================================

    SUBROUTINE finish_multiphase_io()

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        continue
    END SUBROUTINE finish_multiphase_io

    !================================================================

    SUBROUTINE read_vff(vff)
    !----------------------------------------------------------------
    !   What it does:
    !   This subroutine reads the initial values of the volume
    !   fraction field from a csv-file generated in column-major 
    !   order. The grid size needs to be specified via Nx, Ny and Nz.
    !----------------------------------------------------------------

        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: vff
        
        ! Local variables
        INTEGER(intk) :: Nx, Ny, Nz
        INTEGER(intk) :: i, unit

        IF ( test_multiphase == 'CylTra' ) THEN
            Nx = 84
            Ny = 84
            Nz = 8
        END IF

        OPEN(newunit=unit,file="vff.csv",status="old",action="read")
        DO i = 1, Nx*Ny*Nz
            READ(unit,*) vff%arr(i)
        END DO
        CLOSE(unit)

    END SUBROUTINE read_vff

    !================================================================

    SUBROUTINE update_velocity(u_f, v_f, w_f, vff_f, itstep)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        TYPE(field_t), INTENT(in) :: u_f
        TYPE(field_t), INTENT(in) :: v_f
        TYPE(field_t), INTENT(in) :: w_f
        TYPE(field_t), INTENT(in) :: vff_f
        INTEGER(intk), INTENT(in) :: itstep

        ! Local variables
        TYPE(field_t), POINTER :: dx_f, dy_f, dz_f, ddx_f, ddy_f, ddz_f
        REAL(realk), POINTER, CONTIGUOUS, DIMENSION(:, :, :) :: u, v, w, vff
        INTEGER(intk) :: n, igrid!, i, j, k
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk), ALLOCATABLE :: seed(:)
        REAL(realk) :: magnitude, fac = 0.25

        CALL get_field(dx_f, "DX")
        CALL get_field(dy_f, "DY")
        CALL get_field(dz_f, "DZ")

        CALL get_field(ddx_f, "DDX")
        CALL get_field(ddy_f, "DDY")
        CALL get_field(ddz_f, "DDZ")

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL u_f%get_ptr(u, igrid)
            CALL v_f%get_ptr(v, igrid)
            CALL w_f%get_ptr(w, igrid)
            CALL vff_f%get_ptr(vff, igrid)

            CALL dx_f%get_ptr(dx, igrid)
            CALL dy_f%get_ptr(dy, igrid)
            CALL dz_f%get_ptr(dz, igrid)

            CALL ddx_f%get_ptr(ddx, igrid)
            CALL ddy_f%get_ptr(ddy, igrid)
            CALL ddz_f%get_ptr(ddz, igrid)

            IF ( itstep == 1 ) THEN
                CALL random_seed(size = n)
                IF (.NOT. ALLOCATED(seed)) ALLOCATE(seed(n))
                seed = 1234
                CALL random_seed(put = seed)
            END IF

            IF ( test_multiphase == 'CylTra' ) THEN

                IF ( MOD(itstep, 10) == 1 .AND. itstep /= 1 ) THEN
                    CALL random_number(fac)
                END IF

                magnitude = 0.08                       
                u = cos(fac*2.0*3.1415) * magnitude
                v = sin(fac*2.0*3.1415) * magnitude
                w = 0.0

            END IF

        END DO

    END SUBROUTINE

END MODULE multiphase_io_mod