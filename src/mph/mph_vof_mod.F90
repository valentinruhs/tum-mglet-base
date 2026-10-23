!====================================================================
!  Module: mph_vof_mod
!
!   Responsibilities:
!   - Provides the operators of the multi-phase momentum equation
!     in conservative form on the staggered grid:
!     - Advection: directionally split geometric VOF/PLIC advection
!       of c on the main grid and of mass and momentum on the three
!       staggered grids with the same phase fluxes (mass-momentum
!       consistent), Weymouth-Yue compression term; performed once
!       per time step
!     - Diffusion: full viscous stress tensor with variable
!       viscosity, correction at no-slip walls
!     - Pressure: pressure gradient and c-weighted driving pressure
!       gradient, divided by the face density
!     - External forces: gravity
!   - Updates the velocities from the advected staggered momentum
!     and density
!   - Stores c and the velocities of the previous time step
!
!   Author:      Valentin Ruhs
!   e-Mail:      valentin.ruhs@gmx.de
!   Created:     2026-02
!
!   Sources:
!   W. Aniszewski et al., “PArallel, Robust, Interface Simulator
!   (PARIS),” Computer Physics Communications, vol. 263, p. 107849,
!   Jun. 2021, doi: 10.1016/j.cpc.2021.107849.
!   See PARIS (accessed: May 2026) source code subroutine(s):
!       - vofsweeps (vofmodules.f90)
!       - interpole_quad (rudman_modules.f90)
!       - interpole3 (modules.f90)
!       - slope_lim (modules.f90)
!
!   T. Arrufat et al., “A mass-momentum consistent, Volume-of-Fluid
!   method for incompressible flow on staggered grids,” Computers &
!   Fluids, vol. 215, p. 104785, Jan. 2021,
!   doi: 10.1016/j.compfluid.2020.104785.
!
!   D. Herrmann, Numerische Mathematik — 40 BASIC-Programme.
!   Wiesbaden: Vieweg+Teubner Verlag, 1983.
!   doi: 10.1007/978-3-322-96321-5.
!
!   B. P. Leonard, “A stable and accurate convective modelling
!   procedure based on quadratic upstream interpolation,” Computer
!   Methods in Applied Mechanics and Engineering, vol. 19, no. 1,
!   pp. 59–98, Jun. 1979, doi: 10.1016/0045-7825(79)90034-3.
!
!   P. K. Sweby, “High Resolution Schemes Using Flux Limiters for
!   Hyperbolic Conservation Laws,” SIAM Journal on Numerical
!   Analysis, vol. 21, no. 5, pp. 995–1011, Oct. 1984,
!   doi: 10.1137/0721062.
!
!   G. Tryggvason, R. Scardovelli, and S. Zaleski, Direct Numerical
!   Simulations of Gas–Liquid Multiphase Flows, 1st ed. Cambridge
!   University Press, 2011. doi: 10.1017/CBO9780511975264.
!
!   H. Werner and H. Wengle, “Large-eddy simulation of turbulent flow
!   over and around a cube in a plate channel,” in Turbulent Shear
!   Flows 8. Berlin, Heidelberg: Springer, 1993, pp. 155–168,
!   doi: 10.1007/978-3-642-77674-8_12.
!
!   G. D. Weymouth and D. K.-P. Yue, “Conservative Volume-of-Fluid
!   method for free-surface simulations on Cartesian-grids,” Journal
!   of Computational Physics, vol. 229, no. 8, pp. 2853–2865,
!   Apr. 2010, doi: 10.1016/j.jcp.2009.12.018.
!
!   See MGLET source code subroutine(s):
!       - swcle3d (tstle4.F90)
!
!====================================================================

MODULE mph_vof_mod

    USE precision_mod, ONLY: intk, realk
    USE field_mod, ONLY: field_t
    USE fields_mod, ONLY: set_field, get_fieldptr
    USE grids_mod, ONLY: nmygrids, mygrids, get_mgdims
    USE err_mod, ONLY: err_abort
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: ERROR_UNIT
    USE bound_flow_mod, ONLY: bound_flow

    USE mphcore_mod, ONLY: rho1, rho2, gmol1, gmol2, grav, &
        skpAdv, skpDif, skpPre, skpExt, vofTol, advScm, donCen, vofErr
    USE mph_utils_mod, ONLY: sel_ind, sel_ext, sel_vel, clp_tol, int2char, get_lp_mdf, get_lp_bnds
    USE mph_plic_mod, ONLY: comp_ifc, comp_c_stg, comp_isIfc_stg, comp_c_loc
    USE mph_props_mod, ONLY: comp_d_stg
    USE mph_xfer_mod, ONLY: rstr, rstr_stg, prlg, prlg_stg, cnct
    USE mph_bound_mod, ONLY: bound_stg

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: init_mph_vof, finish_mph_vof, cpy_flds, &
        adve_operator, diff_operator, pres_operator, exte_operator

CONTAINS

    SUBROUTINE init_mph_vof()
    !----------------------------------------------------------------
    !   Short: Register staggered c, rho, rho*u, compression
    !          coefficient and flux fields.
    !
    !   CS1-3, DS1-3, MS1-3: volume fraction, density and momentum of
    !   the staggered u-, v- and w-cells. CWY, CWYS1-3: Weymouth-Yue
    !   compression coefficients of the main and staggered grids.
    !   ADVR, ADVE (advector, advectee) and CFLX1, CFLX2 (phase
    !   fluxes) are work arrays of the current split.
    !----------------------------------------------------------------

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        ! Initialize vof fields
        ! Staggered volume fraction fields
        CALL set_field("CS1", istag=1, buffers=.TRUE.)
        CALL set_field("CS2", jstag=1, buffers=.TRUE.)
        CALL set_field("CS3", kstag=1, buffers=.TRUE.)

        ! Staggered density fields
        CALL set_field("DS1", istag=1, buffers=.TRUE.)
        CALL set_field("DS2", jstag=1, buffers=.TRUE.)
        CALL set_field("DS3", kstag=1, buffers=.TRUE.)

        ! Staggered momentum fields
        CALL set_field("MS1", istag=1, buffers=.TRUE.)
        CALL set_field("MS2", jstag=1, buffers=.TRUE.)
        CALL set_field("MS3", kstag=1, buffers=.TRUE.)

        ! Main compression coefficient field
        CALL set_field("CWY", buffers=.TRUE.)

        ! Staggered compression coefficient fields
        CALL set_field("CWYS1", istag=1, buffers=.TRUE.)
        CALL set_field("CWYS2", jstag=1, buffers=.TRUE.)
        CALL set_field("CWYS3", kstag=1, buffers=.TRUE.)

        ! Advector/advectee and flux fields
        CALL set_field("ADVR", buffers=.TRUE.)
        CALL set_field("ADVE", buffers=.TRUE.)
        CALL set_field("CFLX1", buffers=.TRUE.)
        CALL set_field("CFLX2", buffers=.TRUE.)

    END SUBROUTINE init_mph_vof

    !================================================================

    SUBROUTINE finish_mph_vof()
    !----------------------------------------------------------------
    !   Short: Finalize module (no operation).
    !----------------------------------------------------------------

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        CONTINUE

    END SUBROUTINE finish_mph_vof

    !================================================================

    SUBROUTINE cpy_flds()
    !----------------------------------------------------------------
    !   Short: Store c and velocities of the previous time step.
    !
    !   UP, VP and WP are the advecting velocities of all three
    !   directional splits, since U, V and W are updated after each
    !   split (upd_vel_stg).
    !----------------------------------------------------------------

        ! Subroutine arguments
        ! None

        ! Local variables
        INTEGER(intk) :: n, igrid
        REAL(realk), POINTER, CONTIGUOUS :: c(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: u(:,:,:), v(:,:,:), w(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: cp(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: up(:,:,:), vp(:,:,:), wp(:,:,:)

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_fieldptr(c, "C", igrid)
            CALL get_fieldptr(u, "U", igrid)
            CALL get_fieldptr(v, "V", igrid)
            CALL get_fieldptr(w, "W", igrid)
            CALL get_fieldptr(cp, "CP", igrid)
            CALL get_fieldptr(up, "UP", igrid)
            CALL get_fieldptr(vp, "VP", igrid)
            CALL get_fieldptr(wp, "WP", igrid)

            cp = c
            up = u
            vp = v
            wp = w

        END DO

    END SUBROUTINE cpy_flds

    !================================================================

    SUBROUTINE comp_flx(l, dt)
    !----------------------------------------------------------------
    !   Short: Compute VOF fluxes on the main grid in split
    !          direction l.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: l
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        INTEGER(intk) :: n, igrid
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: c(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: up(:,:,:), vp(:,:,:), wp(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: vel(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: cFlx1(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        REAL(realk), POINTER, CONTIGUOUS :: normx(:,:,:), normy(:,:,:), normz(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: alpha(:,:,:)
        REAl(realk), POINTER, CONTIGUOUS :: isIfc(:,:,:)
        INTEGER(intk) :: ista, iend, jsta, jend, ksta, kend

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_fieldptr(c, "C", igrid)
            CALL get_fieldptr(up, "UP", igrid)
            CALL get_fieldptr(vp, "VP", igrid)
            CALL get_fieldptr(wp, "WP", igrid)
            CALL get_fieldptr(cFlx1, "CFLX1", igrid)
            CALL get_fieldptr(ddx, "DDX", igrid)
            CALL get_fieldptr(ddy, "DDY", igrid)
            CALL get_fieldptr(ddz, "DDZ", igrid)
            CALL get_fieldptr(normx, "NORMX", igrid)
            CALL get_fieldptr(normy, "NORMY", igrid)
            CALL get_fieldptr(normz, "NORMZ", igrid)
            CALL get_fieldptr(alpha, "ALPHA", igrid)
            CALL get_fieldptr(isIfc, "ISIFC", igrid)

            CALL sel_vel(l, up, vp, wp, vel)
            CALL get_lp_bnds(igrid, 0, ista, iend, jsta, jend, ksta, kend)
            CALL comp_flx_grd(kk, jj, ii, l, c, vel, cFlx1, &
                ddx, ddy, ddz, normx, normy, normz, alpha, isIfc, dt, &
                ista, iend, jsta, jend, ksta, kend)
        END DO

    END SUBROUTINE comp_flx

    !================================================================

    SUBROUTINE comp_flx_grd(kk, jj, ii, l, c, vel, cFlx1, &
        ddx, ddy, ddz, normx, normy, normz, alpha, isIfc, dt, &
        ista, iend, jsta, jend, ksta, kend)
    !----------------------------------------------------------------
    !   Short: Grid kernel: geometric PLIC flux from the upwind cell,
    !          CFL check.
    !
    !   For every face in direction l, the volume of phase 1 crossing
    !   the face within dt is computed geometrically. The donating
    !   region is the cuboid of width |u|*dt inside the upwind cell
    !   adjacent to the face. In interface cells its volume fraction
    !   follows from the PLIC plane (comp_c_loc); for u > 0 the plane
    !   constant is shifted to the origin of the cuboid. In all other
    !   cells the cell value of c is used. The flux is
    !   cFlx1 = u*c_don, i.e. a volume flux per unit face area.
    !
    !   Faces with |u| <= vofTol carry no flux. Aborts if |u|*dt
    !   exceeds the width of the donor cell (CFL > 1).
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: l
        REAL(realk), INTENT(in) :: c(kk, jj, ii), vel(kk, jj, ii)
        REAL(realk), INTENT(inout) :: cFlx1(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: normx(kk, jj, ii), normy(kk, jj, ii), normz(kk, jj, ii)
        REAL(realk), INTENT(in) :: alpha(kk, jj, ii)
        REAL(realk), INTENT(in) :: isIfc(kk, jj, ii)
        REAL(realk), INTENT(in) :: dt
        INTEGER(intk), INTENT(in) :: ista, iend, jsta, jend, ksta, kend

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: kl, jl, il
        REAL(realk) :: dds, norms
        REAL(realk) :: dimx, dimy, dimz
        REAL(realk) :: flxedProp, flxWidth, flxAlpha

        CALL sel_ind(l, il, jl, kl)

        DO i = ista-il, iend
            DO j = jsta-jl, jend
                DO k = ksta-kl, kend
                    IF ( vel(k,j,i) > vofTol ) THEN

                        ! Compute face flux width and characteristic length
                        flxWidth = abs( vel(k,j,i) )*dt
                        dds = il*ddx(i) + jl*ddy(j) + kl*ddz(k)
                        IF ( flxWidth > (1.0_realk + 1.0E-12_realk)*dds ) THEN
                            WRITE(ERROR_UNIT,'(A,1I2,3I4,4ES12.4)') "VOFCFL l,k,j,i,advr,dds,c,cfl: ", l, &
                                k, j, i, vel(k,j,i), dds, c(k+kl,j+jl,i+il), flxWidth/dds
                            FLUSH(ERROR_UNIT)
                            CALL err_abort(vofErr, "flxWidth > cellWidth.", __FILE__, __LINE__)
                        ENDIF

                        IF ( isIfc(k,j,i) > 0.0_realk ) THEN
                            ! Compute norm
                            norms = il*normx(k,j,i) + jl*normy(k,j,i) + kl*normz(k,j,i)

                            ! Compute proper alpha
                            flxAlpha = alpha(k,j,i) - norms*( dds - flxWidth )

                            ! Compute dimensions of fluxed cuboid
                            dimx = il*flxWidth + (1-il)*ddx(i)
                            dimy = jl*flxWidth + (1-jl)*ddy(j)
                            dimz = kl*flxWidth + (1-kl)*ddz(k)

                            ! Compute c in fluxed cuboid
                            CALL comp_c_loc(flxedProp, flxAlpha, dimx, dimy, dimz, &
                                normx(k,j,i), normy(k,j,i), normz(k,j,i))
                        ELSE
                            flxedProp = c(k,j,i)
                        END IF
                    ELSE IF ( vel(k,j,i) < -vofTol ) THEN

                        ! Compute face flux width and characteristic length
                        flxWidth = abs( vel(k,j,i) )*dt
                        dds = il*ddx(i+il) + jl*ddy(j+jl) + kl*ddz(k+kl)
                        IF ( flxWidth > (1.0_realk + 1.0E-12_realk)*dds ) THEN
                            WRITE(ERROR_UNIT,'(A,1I2,3I4,4ES12.4)') "VOFCFL l,k,j,i,advr,dds,c,cfl: ", l, &
                                k, j, i, vel(k,j,i), dds, c(k+kl,j+jl,i+il), flxWidth/dds
                            FLUSH(ERROR_UNIT)
                            CALL err_abort(vofErr, "flxWidth > cellWidth.", __FILE__, __LINE__)
                        ENDIF

                        IF ( isIfc(k+kl,j+jl,i+il) > 0.0_realk ) THEN
                            ! Compute proper alpha
                            flxAlpha = alpha(k+kl,j+jl,i+il)

                            ! Compute dimensions of fluxed cuboid
                            dimx = il*flxWidth + (1-il)*ddx(i+il)
                            dimy = jl*flxWidth + (1-jl)*ddy(j+jl)
                            dimz = kl*flxWidth + (1-kl)*ddz(k+kl)

                            ! Compute c in fluxed cuboid
                            CALL comp_c_loc(flxedProp, flxAlpha, dimx, dimy, dimz, &
                                normx(k+kl,j+jl,i+il), normy(k+kl,j+jl,i+il), normz(k+kl,j+jl,i+il))
                        ELSE
                            flxedProp = c(k+kl,j+jl,i+il)
                        END IF
                    ELSE
                        flxedProp = 0.0_realk
                    END IF
                    cFlx1(k,j,i) = vel(k,j,i)*flxedProp
                END DO
            END DO
        END DO

    END SUBROUTINE comp_flx_grd

    !================================================================

    SUBROUTINE comp_flx_stg(q, l, dt)
    !----------------------------------------------------------------
    !   Short: Compute phase fluxes on staggered grid q in split
    !          direction l.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: q, l
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        INTEGER(intk) :: kk, jj, ii, n, igrid
        REAL(realk), POINTER, CONTIGUOUS :: c(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: advr(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: cFlx1(:,:,:), cFlx2(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        REAL(realk), POINTER, CONTIGUOUS :: normx(:,:,:), normy(:,:,:), normz(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: alpha(:,:,:)
        REAl(realk), POINTER, CONTIGUOUS :: isIfc(:,:,:)
        INTEGER(intk) :: ista, iend, jsta, jend, ksta, kend

        DO n = 1, nmygrids
            igrid = mygrids(n)
            CALL get_mgdims(kk, jj, ii, igrid)
            CALL get_fieldptr(c, "C", igrid)
            CALL get_fieldptr(advr, "ADVR", igrid)
            CALL get_fieldptr(cFlx1, "CFLX1", igrid)
            CALL get_fieldptr(cFlx2, "CFLX2", igrid)
            CALL get_fieldptr(ddx, "DDX", igrid)
            CALL get_fieldptr(ddy, "DDY", igrid)
            CALL get_fieldptr(ddz, "DDZ", igrid)
            CALL get_fieldptr(normx, "NORMX", igrid)
            CALL get_fieldptr(normy, "NORMY", igrid)
            CALL get_fieldptr(normz, "NORMZ", igrid)
            CALL get_fieldptr(alpha, "ALPHA", igrid)
            CALL get_fieldptr(isIfc, "ISIFC", igrid)

            CALL get_lp_bnds(igrid, q, ista, iend, jsta, jend, ksta, kend)
            CALL comp_flx_stg_grd(kk, jj, ii, q, l, c, advr, cFlx1, cFlx2, &
                ddx, ddy, ddz, normx, normy, normz, alpha, isIfc, dt, &
                ista, iend, jsta, jend, ksta, kend)
        END DO

    END SUBROUTINE comp_flx_stg

    !================================================================

    SUBROUTINE comp_flx_stg_grd(kk, jj, ii, q, l, c, advr, cFlx1, cFlx2, &
        ddx, ddy, ddz, normx, normy, normz, alpha, isIfc, dt, &
        ista, iend, jsta, jend, ksta, kend)
    !----------------------------------------------------------------
    !   Short: Grid kernel: PLIC phase fluxes from donor (half)
    !          cells, CFL check.
    !
    !   Phase fluxes cFlx1 (phase 1) and cFlx2 (phase 2) through the
    !   faces of the staggered q-cells in direction l, based on the
    !   PLIC interface of the main grid and the advector advr.
    !
    !   q == l: the face lies in the centre of the main cell
    !           (k+kl,j+jl,i+il). The donating region of width
    !           |advr|*dt lies in the upstream half of this cell,
    !           adjacent to its centre.
    !   q /= l: the face coincides with a main grid face in direction
    !           l, but the staggered cell covers the upper q-half of
    !           the main cell Mi and the lower q-half of the main
    !           cell Pl. The donating region is split into these two
    !           half cells; their volume fractions are weighted with
    !           the cell widths in q-direction.
    !
    !   Faces with |advr| < vofTol carry no flux. Aborts if |advr|*dt
    !   exceeds half of the donor cell width.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, q, l
        REAL(realk), INTENT(in) :: c(kk, jj, ii)
        REAL(realk), INTENT(in) :: advr(kk, jj, ii)
        REAL(realk), INTENT(inout) :: cFlx1(kk, jj, ii), cFlx2(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: normx(kk, jj, ii), normy(kk, jj, ii), normz(kk, jj, ii)
        REAL(realk), INTENT(in) :: alpha(kk, jj, ii)
        REAL(realk), INTENT(in) :: isIfc(kk, jj, ii)
        REAL(realk), INTENT(in) :: dt
        INTEGER(intk), INTENT(in) :: ista, iend, jsta, jend, ksta, kend

        ! Local variables
        INTEGER(intk) :: k, j, i, kl, jl, il, kq, jq, iq
        INTEGER(intk) :: kDonMi, jDonMi, iDonMi, kDonPl, jDonPl, iDonPl
        REAL(realk) :: ddsDon, normDon, flxWidth, flxAlpha, flxDimx, flxDimy, flxDimz, flxedProp
        REAL(realk) :: farEnd, ddslMi, ddslPl, ddsqMi, ddsqPl
        REAL(realk) :: normlMi, normlPl, normqMi, normqPl, flxAlphaMi, flxAlphaPl
        REAL(realk) :: flxDimxMi, flxDimyMi, flxDimzMi, flxDimxPl, flxDimyPl, flxDimzPl, fracMi, fracPl

        CALL sel_ind(l, il, jl, kl)
        CALL sel_ind(q, iq, jq, kq)

        IF ( q == l ) THEN
            DO i = ista-il, iend
                DO j = jsta-jl, jend
                    DO k = ksta-kl, kend
                        IF ( ABS(advr(k,j,i)) < vofTol ) THEN
                            cFlx1(k,j,i) = 0.0_realk
                            cFlx2(k,j,i) = 0.0_realk
                            CYCLE
                        ENDIF

                        ddsDon  = il*ddx(i+il) + jl*ddy(j+jl) + kl*ddz(k+kl)
                        normDon = il*normx(k+kl,j+jl,i+il) &
                                + jl*normy(k+kl,j+jl,i+il) &
                                + kl*normz(k+kl,j+jl,i+il)

                        flxWidth = ABS( advr(k,j,i) )*dt

                        IF ( flxWidth > 0.5_realk*ddsDon ) THEN
                            WRITE(ERROR_UNIT,'(A,2I2,3I4,4ES12.4)') "VOFCFL q,l,k,j,i,advr,dds,c,cfl: ", q, l, &
                                k, j, i, advr(k,j,i), ddsDon, c(k+kl,j+jl,i+il), flxWidth/ddsDon
                            FLUSH(ERROR_UNIT)
                            CALL err_abort(vofErr, "flxWidth > 0.5*cellWidth.", __FILE__, __LINE__)
                        ENDIF

                        IF ( isIfc(k+kl,j+jl,i+il) > 0.0_realk ) THEN
                            flxAlpha = alpha(k+kl,j+jl,i+il) - normDon*&
                                ( 0.5_realk*ddsDon &
                                - MERGE(flxWidth, 0.0_realk, advr(k,j,i) > 0.0_realk) )
                            flxDimx = il*flxWidth + (1-il)*ddx(i)
                            flxDimy = jl*flxWidth + (1-jl)*ddy(j)
                            flxDimz = kl*flxWidth + (1-kl)*ddz(k)
                            CALL comp_c_loc(flxedProp, flxAlpha, flxDimx, flxDimy, flxDimz, &
                                normx(k+kl,j+jl,i+il), normy(k+kl,j+jl,i+il), &
                                normz(k+kl,j+jl,i+il))
                        ELSE
                            flxedProp = c(k+kl,j+jl,i+il)
                        ENDIF

                        cFlx1(k,j,i) = advr(k,j,i)*flxedProp
                        cFlx2(k,j,i) = advr(k,j,i)*( 1.0_realk - flxedProp )
                    END DO
                END DO
            END DO
        ELSE
            DO i = ista-il, iend
                DO j = jsta-jl, jend
                    DO k = ksta-kl, kend

                    IF ( ABS(advr(k,j,i)) < vofTol ) THEN
                        cFlx1(k,j,i) = 0.0_realk
                        cFlx2(k,j,i) = 0.0_realk
                        CYCLE
                    ENDIF

                    ! Indices of the half pressure cells overlapping the velocity cells
                    IF ( advr(k,j,i) > 0.0_realk ) THEN
                        iDonMi = i ; jDonMi = j ; kDonMi = k ; farEnd = 1.0_realk
                    ELSE
                        iDonMi = i+il ; jDonMi = j+jl ; kDonMi = k+kl ; farEnd = 0.0_realk
                    ENDIF

                    iDonPl = iDonMi+iq ; jDonPl = jDonMi+jq ; kDonPl = kDonMi+kq

                    flxWidth = ABS( advr(k,j,i) )*dt

                    ddslMi = il*ddx(iDonMi) + jl*ddy(jDonMi) + kl*ddz(kDonMi)
                    ddslPl = il*ddx(iDonPl) + jl*ddy(jDonPl) + kl*ddz(kDonPl)

                    IF ( flxWidth > 0.5_realk*ddslMi .OR. flxWidth > 0.5_realk*ddslPl ) THEN
                            WRITE(ERROR_UNIT,'(A,2I2,3I4,4ES12.4)') "VOFCFL q,l,k,j,i,advr,dds,c,cfl: ", q, l, &
                                k, j, i, advr(k,j,i), ddslMi, c(k+kl,j+jl,i+il), flxWidth/ddslMi
                            FLUSH(ERROR_UNIT)
                            CALL err_abort(vofErr, "flxWidth > 0.5*cellWidth.", __FILE__, __LINE__)
                        ENDIF

                    ddsqMi = iq*ddx(iDonMi) + jq*ddy(jDonMi) + kq*ddz(kDonMi)
                    ddsqPl = iq*ddx(iDonPl) + jq*ddy(jDonPl) + kq*ddz(kDonPl)

                    normlMi = il*normx(kDonMi,jDonMi,iDonMi) + jl*normy(kDonMi,jDonMi,iDonMi) + kl*normz(kDonMi,jDonMi,iDonMi)
                    normlPl = il*normx(kDonPl,jDonPl,iDonPl) + jl*normy(kDonPl,jDonPl,iDonPl) + kl*normz(kDonPl,jDonPl,iDonPl)
                    normqMi = iq*normx(kDonMi,jDonMi,iDonMi) + jq*normy(kDonMi,jDonMi,iDonMi) + kq*normz(kDonMi,jDonMi,iDonMi)
                    normqPl = iq*normx(kDonPl,jDonPl,iDonPl) + jq*normy(kDonPl,jDonPl,iDonPl) + kq*normz(kDonPl,jDonPl,iDonPl)

                    IF ( isIfc(kDonMi,jDonMi,iDonMi) > 0.0_realk ) THEN
                        flxAlphaMi = alpha(kDonMi,jDonMi,iDonMi) &
                                    - normlMi*farEnd*( ddslMi - flxWidth ) &
                                    - normqMi*0.5_realk*ddsqMi

                        flxDimxMi = il*flxWidth + iq*ddx(iDonMi)/2.0_realk + (1-il-iq)*ddx(iDonMi)
                        flxDimyMi = jl*flxWidth + jq*ddy(jDonMi)/2.0_realk + (1-jl-jq)*ddy(jDonMi)
                        flxDimzMi = kl*flxWidth + kq*ddz(kDonMi)/2.0_realk + (1-kl-kq)*ddz(kDonMi)

                        CALL comp_c_loc(fracMi, flxAlphaMi, flxDimxMi, flxDimyMi, flxDimzMi, &
                            normx(kDonMi,jDonMi,iDonMi), normy(kDonMi,jDonMi,iDonMi), normz(kDonMi,jDonMi,iDonMi))
                    ELSE
                        fracMi = c(kDonMi,jDonMi,iDonMi)
                    ENDIF

                    IF ( isIfc(kDonPl,jDonPl,iDonPl) > 0.0_realk ) THEN
                        flxAlphaPl = alpha(kDonPl,jDonPl,iDonPl) &
                                    - normlPl*farEnd*( ddslPl - flxWidth )

                        flxDimxPl = il*flxWidth + iq*ddx(iDonPl)/2.0_realk + (1-il-iq)*ddx(iDonPl)
                        flxDimyPl = jl*flxWidth + jq*ddy(jDonPl)/2.0_realk + (1-jl-jq)*ddy(jDonPl)
                        flxDimzPl = kl*flxWidth + kq*ddz(kDonPl)/2.0_realk + (1-kl-kq)*ddz(kDonPl)

                        CALL comp_c_loc(fracPl, flxAlphaPl, flxDimxPl, flxDimyPl, flxDimzPl, &
                            normx(kDonPl,jDonPl,iDonPl), normy(kDonPl,jDonPl,iDonPl), &
                            normz(kDonPl,jDonPl,iDonPl))
                    ELSE
                        fracPl = c(kDonPl,jDonPl,iDonPl)
                    ENDIF

                    flxedProp = ( fracMi*ddsqMi + fracPl*ddsqPl ) / ( ddsqMi + ddsqPl )

                    cFlx1(k,j,i) = advr(k,j,i)*flxedProp
                    cFlx2(k,j,i) = advr(k,j,i)*( 1.0_realk - flxedProp )

                    END DO
                END DO
            END DO
        ENDIF

    END SUBROUTINE comp_flx_stg_grd

    !================================================================

    PURE SUBROUTINE def_adv_seq(iteration, advSeq)
    !----------------------------------------------------------------
    !   Short: Define cyclic permutation of the split
    !          advection directions.
    !
    !   The sequence in which c and the momentum are advected is
    !   permuted cyclically with the time step (period of three, as
    !   in PARIS):
    !
    !   itstep = 1, 4, ...: x -> y -> z
    !   itstep = 2, 5, ...: z -> x -> y
    !   itstep = 3, 6, ...: y -> z -> x
    !
    !   Sources:
    !   W. Aniszewski et al., “PArallel, Robust, Interface Simulator
    !   (PARIS),” Computer Physics Communications, vol. 263,
    !   p. 107849, Jun. 2021, doi: 10.1016/j.cpc.2021.107849.
    !   See PARIS (accessed: May 2026) source code subroutine(s):
    !       - vofsweeps (vofmodules.f90)
    !
    !   T. Arrufat et al., “A mass-momentum consistent,
    !   Volume-of-Fluid method for incompressible flow on staggered
    !   grids,” Computers & Fluids, vol. 215, p. 104785, Jan. 2021,
    !   doi: 10.1016/j.compfluid.2020.104785.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: iteration
        INTEGER(intk), INTENT(inout) :: advSeq(3)

        ! Local variables
        INTEGER(intk) :: perInd

        ! perInd only changes in a new time-step
        perInd = mod(iteration-1, 3)

        ! Select permutation of split advection
        SELECT CASE ( perInd )
            CASE (0); advSeq = [1, 2, 3]
            CASE (1); advSeq = [3, 1, 2]
            CASE (2); advSeq = [2, 3, 1]
        END SELECT

    END SUBROUTINE def_adv_seq

    !================================================================

    SUBROUTINE adve_operator(dt, itstep)
    !----------------------------------------------------------------
    !   Short: Advect c and the staggered mass and momentum by
    !          directional splitting, incl. grid transfers.
    !
    !   Spatial and temporal integration of the advection terms in
    !   one step. Geometric VOF/PLIC is "exact" up to the accuracy of
    !   the interface reconstruction, hence no RK sub-steps are used.
    !
    !   1) Compression coefficient of the main grid; staggered c,
    !      density, momentum and compression coefficients from the
    !      current PLIC interface.
    !   2) Restriction and prolongation of the staggered fields and
    !      of the velocities.
    !   3) For each split direction l (sequence from def_adv_seq):
    !      advect momentum, c and density of the three staggered
    !      grids with the same phase fluxes, then c of the main grid;
    !      restrict and prolongate c and reconstruct the interface;
    !      restrict and prolongate the staggered fields and update
    !      the velocities u_q = (rho*u)_q/rho_q.
    !
    !   Each directional split has to be completed on the whole
    !   domain (including restriction and prolongation) before
    !   the next split starts. Otherwise, errors occur at the
    !   grid boundaries.
    !
    !   With skpAdv only c of the main grid is advected; momentum and
    !   staggered fields are not advected (velocity unchanged).
    !
    !   Sources:
    !   T. Arrufat et al., “A mass-momentum consistent,
    !   Volume-of-Fluid method for incompressible flow on staggered
    !   grids,” Computers & Fluids, vol. 215, p. 104785, Jan. 2021,
    !   doi: 10.1016/j.compfluid.2020.104785.
    !
    !   G. D. Weymouth and D. K.-P. Yue, “Conservative
    !   Volume-of-Fluid method for free-surface simulations on
    !   Cartesian-grids,” Journal of Computational Physics, vol. 229,
    !   no. 8, pp. 2853–2865, Apr. 2010,
    !   doi: 10.1016/j.jcp.2009.12.018.
    !----------------------------------------------------------------

        ! Subroutine arguments
        REAL(realk), INTENT(in) :: dt
        INTEGER(intk), INTENT(in) :: itstep

        ! Local variables
        INTEGER(intk) :: q, advSeq(3), dirLoop, l

        ! Initialize cWy
        CALL comp_cWy()

        ! Initialize stg. grid cS(.), dS(.), mS(.), cWyS(.)
        DO q = 1, 3
            CALL comp_c_stg(q)
            CALL comp_d_stg(q)
            IF ( skpAdv ) CYCLE
            CALL comp_m_stg(q)
            CALL comp_cWy_stg(q)
        END DO

        ! Restriction
        DO q = 1, 3
            IF ( skpAdv ) CYCLE
            CALL rstr_stg(q=q, fldName="CS")
            CALL rstr_stg(q=q, fldName="DS")
            CALL rstr_stg(q=q, fldName="MS")
        END DO

        ! Prolongation
        IF ( .NOT. skpAdv ) THEN
            CALL prlg_stg(fldName1="CS1", fldName2="CS2", fldName3="CS3", bnd=bound_stg)
            CALL prlg_stg(fldName1="DS1", fldName2="DS2", fldName3="DS3", bnd=bound_stg)
            CALL prlg_stg(fldName1="U", fldName2="V", fldName3="W", bnd=bound_flow)
        END IF

        CALL def_adv_seq(itstep, advSeq)
        DO dirLoop = 1, 3
            l = advSeq(dirLoop)
            ! Advect stg. grid mS(.), cS(.), dS(.)
            DO q = 1, 3
                IF ( skpAdv ) CYCLE
                CALL comp_isIfc_stg(q)
                CALL cnct("ISIFCVICS"//int2char(q))
                CALL comp_advr_stg(q, l)
                CALL comp_adve_stg(q, l, dt)
                CALL comp_flx_stg(q, l, dt)
                CALL adv_m_stg(q, l, dt)
                CALL adv_c_stg(q, l, dt)
                CALL comp_d_stg(q)
            END DO
            ! Advect main grid c
            CALL comp_flx(l, dt)
            CALL adv_c(l, dt)

            CALL rstr(fldName="C", flag="D")
            CALL prlg(fldName="C")
            CALL comp_ifc()

            DO q = 1, 3
                IF ( skpAdv ) CYCLE
                CALL rstr_stg(q=q, fldName="CS")
                CALL rstr_stg(q=q, fldName="DS")
                CALL rstr_stg(q=q, fldName="MS")
            END DO

            IF ( .NOT. skpAdv ) THEN
                CALL prlg_stg(fldName1="CS1", fldName2="CS2", fldName3="CS3", bnd=bound_stg)
                CALL prlg_stg(fldName1="DS1", fldName2="DS2", fldName3="DS3", bnd=bound_stg)
                DO q = 1, 3
                    CALL upd_vel_stg(q)
                END DO
                CALL prlg_stg(fldName1="U", fldName2="V", fldName3="W", bnd=bound_flow)
            END IF
        END DO

    END SUBROUTINE adve_operator

    !================================================================

    SUBROUTINE diff_operator(uo_f, vo_f, wo_f)
    !----------------------------------------------------------------
    !   Short: Add viscous term with variable viscosity to the
    !          momentum RHS.
    !
    !   For each grid, the wall correction at no-slip walls
    !   (swcle3d_grid) and the divergence of the viscous stress
    !   tensor (diff_operator_grd) are added to uo, vo, wo. Uses
    !   the cell viscosity G, the edge viscosities GUV, GUW, GVW
    !   and the face densities DBA, DLE, DTO from comp_props.
    !   Skipped if skpDif.
    !----------------------------------------------------------------

        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: uo_f, vo_f, wo_f

        ! Local variables
        INTEGER(intk) :: n, igrid
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: u(:,:,:), v(:,:,:), w(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: uo(:,:,:), vo(:,:,:), wo(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: g(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: gUv(:,:,:), gUw(:,:,:), gVw(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: dBa(:,:,:), dLe(:,:,:), dTo(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        REAL(realk), POINTER, CONTIGUOUS :: rdx(:), rdy(:), rdz(:)
        REAL(realk), POINTER, CONTIGUOUS :: rddx(:), rddy(:), rddz(:)
        INTEGER(intk) :: nfro, nbac, nrgt, nlft, nbot, ntop
        INTEGER(intk) :: nfu, nbu, nrv, nlv, nbw, ntw

        IF ( skpDif ) RETURN

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_fieldptr(u, "U", igrid)
            CALL get_fieldptr(v, "V", igrid)
            CALL get_fieldptr(w, "W", igrid)
            CALL uo_f%get_ptr(uo, igrid)
            CALL vo_f%get_ptr(vo, igrid)
            CALL wo_f%get_ptr(wo, igrid)
            CALL get_fieldptr(g, "G", igrid)
            CALL get_fieldptr(gUv, "GUV", igrid)
            CALL get_fieldptr(gUw, "GUW", igrid)
            CALL get_fieldptr(gVw, "GVW", igrid)
            CALL get_fieldptr(dBa, "DBA", igrid)
            CALL get_fieldptr(dLe, "DLE", igrid)
            CALL get_fieldptr(dTo, "DTO", igrid)
            CALL get_fieldptr(ddx, "DDX", igrid)
            CALL get_fieldptr(ddy, "DDY", igrid)
            CALL get_fieldptr(ddz, "DDZ", igrid)
            CALL get_fieldptr(rdx, "RDX", igrid)
            CALL get_fieldptr(rdy, "RDY", igrid)
            CALL get_fieldptr(rdz, "RDZ", igrid)
            CALL get_fieldptr(rddx, "RDDX", igrid)
            CALL get_fieldptr(rddy, "RDDY", igrid)
            CALL get_fieldptr(rddz, "RDDZ", igrid)

            CALL get_lp_mdf(igrid, nfro, nbac, nrgt, nlft, nbot, ntop, &
                nfu, nbu, nrv, nlv, nbw, ntw)
            CALL swcle3d_grid(kk, jj, ii, nfro, nbac, nrgt, nlft, nbot, ntop, &
                u, v, w, gUv, gUw, gVw, dBa, dLe, dTo, &
                ddx, ddy, ddz, rdx, rdy, rdz, uo, vo, wo)
            CALL diff_operator_grd(kk, jj, ii, u, v, w, &
                g, gUv, gUw, gVw, dBa, dLe, dTo, &
                rdx, rdy, rdz, rddx, rddy, rddz, uo, vo, wo, &
                nfu, nbu, nrv, nlv, nbw, ntw)
        END DO

    END SUBROUTINE diff_operator

    !================================================================

    SUBROUTINE diff_operator_grd(kk, jj, ii, u, v, w, g, &
        gUv, gUw, gVw, dBa, dLe, dTo, rdx, rdy, rdz, &
        rddx, rddy, rddz, uo, vo, wo, nfu, nbu, nrv, nlv, nbw, ntw)
    !----------------------------------------------------------------
    !   Short: Grid kernel: divergence of the full viscous
    !          stress tensor.
    !
    !   uo_i += 1/rho_face*d(tau_ij)/dx_j with
    !   tau_ij = mu*(du_i/dx_j + du_j/dx_i). The full tensor (not the
    !   Laplacian form) is required for variable viscosity. Normal
    !   stresses use the cell viscosity g, shear stresses the edge
    !   viscosities gUv, gUw, gVw. At CON faces the loops in normal
    !   direction include the boundary face (loop modifiers nfu, ...,
    !   ntw, see get_lp_mdf).
    !
    !   Source:
    !   G. Tryggvason, R. Scardovelli, and S. Zaleski, Direct
    !   Numerical Simulations of Gas–Liquid Multiphase Flows, 1st ed.
    !   Cambridge University Press, 2011.
    !   doi: 10.1017/CBO9780511975264.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: g(kk, jj, ii)
        REAL(realk), INTENT(in) :: gUv(kk, jj, ii), gUw(kk, jj, ii), gVw(kk, jj, ii)
        REAL(realk), INTENT(in) :: dBa(kk, jj, ii), dLe(kk, jj, ii), dTo(kk, jj, ii)
        REAL(realk), INTENT(in) :: rdx(ii), rdy(jj), rdz(kk)
        REAL(realk), INTENT(in) :: rddx(ii), rddy(jj), rddz(kk)
        REAL(realk), INTENT(inout) :: uo(kk, jj, ii), vo(kk, jj, ii), wo(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: nfu, nbu, nrv, nlv, nbw, ntw

        ! Local variables
        INTEGER(intk) :: k, j, i
        REAL(realk) :: tauXxPl, tauXxMi, tauXyPl, tauXyMi, tauXzPl, tauXzMi
        REAL(realk) :: tauYxPl, tauYxMi, tauYyPl, tauYyMi, tauYzPl, tauYzMi
        REAL(realk) :: tauZxPl, tauZxMi, tauZyPl, tauZyMi, tauZzPl, tauZzMi

        DO i = 3-nfu, ii-3+nbu
            DO j = 3, jj-2
                DO k = 3, kk-2
                    tauXxPl = g(k,j,i+1)*2.0_realk*(u(k,j,i+1) - u(k,j,i))*rddx(i+1)
                    tauXxMi = g(k,j,i)*2.0_realk*(u(k,j,i) - u(k,j,i-1))*rddx(i)
                    tauXyPl = gUv(k,j,i)*((u(k,j+1,i) - u(k,j,i))*rdy(j) + (v(k,j,i+1) - v(k,j,i))*rdx(i))
                    tauXyMi = gUv(k,j-1,i)*((u(k,j,i) - u(k,j-1,i))*rdy(j-1) + (v(k,j-1,i+1) - v(k,j-1,i))*rdx(i))
                    tauXzPl = gUw(k,j,i)*((u(k+1,j,i) - u(k,j,i))*rdz(k) + (w(k,j,i+1) - w(k,j,i))*rdx(i))
                    tauXzMi = gUw(k-1,j,i)*((u(k,j,i) - u(k-1,j,i))*rdz(k-1) + (w(k-1,j,i+1) - w(k-1,j,i))*rdx(i))

                    uo(k,j,i) = uo(k,j,i) + 1.0_realk/dBa(k,j,i)* &
                        ((tauXxPl - tauXxMi)*rdx(i) + &
                         (tauXyPl - tauXyMi)*rddy(j) + &
                         (tauXzPl - tauXzMi)*rddz(k))
                END DO
            END DO
        END DO

        DO i = 3, ii-2
            DO j = 3-nrv, jj-3+nlv
                DO k = 3, kk-2
                    tauYxPl = gUv(k,j,i)*((u(k,j+1,i) - u(k,j,i))*rdy(j) + (v(k,j,i+1) - v(k,j,i))*rdx(i))
                    tauYxMi = gUv(k,j,i-1)*((u(k,j+1,i-1) - u(k,j,i-1))*rdy(j) + (v(k,j,i) - v(k,j,i-1))*rdx(i-1))
                    tauYyPl = g(k,j+1,i)*2.0_realk*(v(k,j+1,i) - v(k,j,i))*rddy(j+1)
                    tauYyMi = g(k,j,i)*2.0_realk*(v(k,j,i) - v(k,j-1,i))*rddy(j)
                    tauYzPl = gVw(k,j,i)*((v(k+1,j,i) - v(k,j,i))*rdz(k) + (w(k,j+1,i) - w(k,j,i))*rdy(j))
                    tauYzMi = gVw(k-1,j,i)*((v(k,j,i) - v(k-1,j,i))*rdz(k-1) + (w(k-1,j+1,i) - w(k-1,j,i))*rdy(j))

                    vo(k,j,i) = vo(k,j,i) + 1.0_realk/dLe(k,j,i)* &
                        ((tauYxPl - tauYxMi)*rddx(i) + &
                         (tauYyPl - tauYyMi)*rdy(j) + &
                         (tauYzPl - tauYzMi)*rddz(k))
                END DO
            END DO
        END DO

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3-nbw, kk-3+ntw
                    tauZxPl = gUw(k,j,i)*((u(k+1,j,i) - u(k,j,i))*rdz(k) + (w(k,j,i+1) - w(k,j,i))*rdx(i))
                    tauZxMi = gUw(k,j,i-1)*((u(k+1,j,i-1) - u(k,j,i-1))*rdz(k) + (w(k,j,i) - w(k,j,i-1))*rdx(i-1))
                    tauZyPl = gVw(k,j,i)*((v(k+1,j,i) - v(k,j,i))*rdz(k) + (w(k,j+1,i) - w(k,j,i))*rdy(j))
                    tauZyMi = gVw(k,j-1,i)*((v(k+1,j-1,i) - v(k,j-1,i))*rdz(k) + (w(k,j,i) - w(k,j-1,i))*rdy(j-1))
                    tauZzPl = g(k+1,j,i)*2.0_realk*(w(k+1,j,i) - w(k,j,i))*rddz(k+1)
                    tauZzMi = g(k,j,i)*2.0_realk*(w(k,j,i) - w(k-1,j,i))*rddz(k)

                    wo(k,j,i) = wo(k,j,i) + 1.0_realk/dTo(k,j,i)*&
                        ((tauZxPl - tauZxMi)*rddx(i) + &
                         (tauZyPl - tauZyMi)*rddy(j) + &
                         (tauZzPl - tauZzMi)*rdz(k))
                END DO
            END DO
        END DO

    END SUBROUTINE diff_operator_grd

    !================================================================

    SUBROUTINE swcle3d_grid(kk, jj, ii, nfro, nbac, nrgt, nlft, nbot, ntop, &
        u, v, w, gUv, gUw, gVw, dBa, dLe, dTo, &
        ddx, ddy, ddz, rdx, rdy, rdz, uo, vo, wo)
    !----------------------------------------------------------------
    !   Short: Grid kernel: correct the viscous wall flux at
    !          no-slip walls.
    !
    !   Face by face wall correction of the tangential velocity
    !   components at no-slip walls (boundary code 5 = NOS). wall_acc
    !   removes the shear stress computed with the ghost cell value
    !   in diff_operator_grd and adds the wall shear stress
    !   2*mu*u/dds with the wall half a cell away (equivalent to the
    !   Werner-Wengle laminar sublayer case).
    !
    !   The wall edge of the viscosity fields and of rd* lies between
    !   ghost and first interior cell: index 2 at the lower faces,
    !   index *-2 at the upper faces.
    !
    !   Sources:
    !   H. Werner and H. Wengle, “Large-eddy simulation of turbulent
    !   flow over and around a cube in a plate channel,” in Turbulent
    !   Shear Flows 8. Berlin, Heidelberg: Springer, 1993,
    !   pp. 155–168, doi: 10.1007/978-3-642-77674-8_12.
    !
    !   See MGLET source code subroutine(s):
    !       - swcle3d (tstle4.F90)
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        INTEGER(intk), INTENT(in) :: nfro, nbac, nrgt, nlft, nbot, ntop
        REAL(realk), INTENT(in) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: gUv(kk, jj, ii), gUw(kk, jj, ii), gVw(kk, jj, ii)
        REAL(realk), INTENT(in) :: dBa(kk, jj, ii), dLe(kk, jj, ii), dTo(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: rdx(ii), rdy(jj), rdz(kk)
        REAL(realk), INTENT(inout) :: uo(kk, jj, ii), vo(kk, jj, ii), wo(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i

        ! Front
        IF ( nfro == 5 ) THEN
            i = 3
            DO j = 2, jj-1
                DO k = 2, kk-1
                    vo(k,j,i) = vo(k,j,i) - wall_acc(gUv(k,j,i-1), dLe(k,j,i), &
                        ddx(i), rdx(i-1), v(k,j,i), v(k,j,i-1))
                    wo(k,j,i) = wo(k,j,i) - wall_acc(gUw(k,j,i-1), dTo(k,j,i), &
                        ddx(i), rdx(i-1), w(k,j,i), w(k,j,i-1))
                END DO
            END DO
        END IF

        ! Back
        IF ( nbac == 5 ) THEN
            i = ii-2
            DO j = 2, jj-1
                DO k = 2, kk-1
                    vo(k,j,i) = vo(k,j,i) - wall_acc(gUv(k,j,i), dLe(k,j,i), &
                        ddx(i), rdx(i), v(k,j,i), v(k,j,i+1))
                    wo(k,j,i) = wo(k,j,i) - wall_acc(gUw(k,j,i), dTo(k,j,i), &
                        ddx(i), rdx(i), w(k,j,i), w(k,j,i+1))
                END DO
            END DO
        END IF

        ! Right
        IF ( nrgt == 5 ) THEN
            j = 3
            DO i = 2, ii-1
                DO k = 2, kk-1
                    uo(k,j,i) = uo(k,j,i) - wall_acc(gUv(k,j-1,i), dBa(k,j,i), &
                        ddy(j), rdy(j-1), u(k,j,i), u(k,j-1,i))
                    wo(k,j,i) = wo(k,j,i) - wall_acc(gVw(k,j-1,i), dTo(k,j,i), &
                        ddy(j), rdy(j-1), w(k,j,i), w(k,j-1,i))
                END DO
            END DO
        END IF

        ! Left
        IF ( nlft == 5 ) THEN
            j = jj-2
            DO i = 2, ii-1
                DO k = 2, kk-1
                    uo(k,j,i) = uo(k,j,i) - wall_acc(gUv(k,j,i), dBa(k,j,i), &
                        ddy(j), rdy(j), u(k,j,i), u(k,j+1,i))
                    wo(k,j,i) = wo(k,j,i) - wall_acc(gVw(k,j,i), dTo(k,j,i), &
                        ddy(j), rdy(j), w(k,j,i), w(k,j+1,i))
                END DO
            END DO
        END IF

        ! Bottom
        IF ( nbot == 5 ) THEN
            k = 3
            DO i = 2, ii-1
                DO j = 2, jj-1
                    uo(k,j,i) = uo(k,j,i) - wall_acc(gUw(k-1,j,i), dBa(k,j,i), &
                        ddz(k), rdz(k-1), u(k,j,i), u(k-1,j,i))
                    vo(k,j,i) = vo(k,j,i) - wall_acc(gVw(k-1,j,i), dLe(k,j,i), &
                        ddz(k), rdz(k-1), v(k,j,i), v(k-1,j,i))
                END DO
            END DO
        END IF

        ! Top
        IF ( ntop == 5 ) THEN
            k = kk-2
            DO i = 2, ii-1
                DO j = 2, jj-1
                    uo(k,j,i) = uo(k,j,i) - wall_acc(gUw(k,j,i), dBa(k,j,i), &
                        ddz(k), rdz(k), u(k,j,i), u(k+1,j,i))
                    vo(k,j,i) = vo(k,j,i) - wall_acc(gVw(k,j,i), dLe(k,j,i), &
                        ddz(k), rdz(k), v(k,j,i), v(k+1,j,i))
                END DO
            END DO
        END IF

    CONTAINS

        PURE ELEMENTAL REAL(realk) FUNCTION wall_acc(mu, d, dds, rds, vel, velGhost) RESULT(acc)
        !------------------------------------------------------------
        !   Short: Acceleration that replaces the ghost-cell shear
        !          stress by the wall shear stress.
        !
        !   Equivalent to the Werner-Wengle laminar sublayer case,
        !   but with multi-phase fluid properties.
        !------------------------------------------------------------
            REAL(realk), INTENT(in) :: mu, d, dds, rds, vel, velGhost
            acc = (2.0_realk*mu*vel/dds - mu*(vel - velGhost)*rds)/(d*dds)
        END FUNCTION wall_acc

    END SUBROUTINE swcle3d_grid

    !================================================================

    SUBROUTINE pres_operator(uo_f, vo_f, wo_f, gradp)
    !----------------------------------------------------------------
    !   Short: Add pressure gradient and c-weighted driving gradient,
    !          divided by the face density, to the momentum RHS.
    !
    !   uo_q -= (dp/dx_q + cS_q*gradp_q)/rho_face. The prescribed
    !   driving pressure gradient gradp (e.g. open channel flow) acts
    !   on phase 1 only, weighted with the staggered volume fraction
    !   cS_q. Skipped if skpPre.
    !----------------------------------------------------------------

        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: uo_f, vo_f, wo_f
        REAL(realk), INTENT(in) :: gradp(3)

        ! Local variables
        INTEGER(intk) :: n, igrid
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: cS1(:,:,:), cS2(:,:,:), cS3(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: p(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: uo(:,:,:), vo(:,:,:), wo(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: dBa(:,:,:), dLe(:,:,:), dTo(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: rdx(:), rdy(:), rdz(:)
        INTEGER(intk) :: nfro, nbac, nrgt, nlft, nbot, ntop
        INTEGER(intk) :: nfu, nbu, nrv, nlv, nbw, ntw

        IF ( skpPre ) RETURN

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_fieldptr(cS1, "CS1", igrid)
            CALL get_fieldptr(cS2, "CS2", igrid)
            CALL get_fieldptr(cS3, "CS3", igrid)
            CALL get_fieldptr(p, "P", igrid)
            CALL uo_f%get_ptr(uo, igrid)
            CALL vo_f%get_ptr(vo, igrid)
            CALL wo_f%get_ptr(wo, igrid)
            CALL get_fieldptr(dBa, "DBA", igrid)
            CALL get_fieldptr(dLe, "DLE", igrid)
            CALL get_fieldptr(dTo, "DTO", igrid)
            CALL get_fieldptr(rdx, "RDX", igrid)
            CALL get_fieldptr(rdy, "RDY", igrid)
            CALL get_fieldptr(rdz, "RDZ", igrid)

            CALL get_lp_mdf(igrid, nfro, nbac, nrgt, nlft, nbot, ntop, &
                nfu, nbu, nrv, nlv, nbw, ntw)
            CALL pres_operator_grd(kk, jj, ii, cS1, cS2, cS3, p, dBa, dLe, dTo, &
                rdx, rdy, rdz, uo, vo, wo, gradp, nfu, nbu, nrv, nlv, nbw, ntw)
        END DO

    END SUBROUTINE pres_operator

    !================================================================

    SUBROUTINE pres_operator_grd(kk, jj, ii, cS1, cS2, cS3, p, dBa, dLe, dTo, &
        rdx, rdy, rdz, uo, vo, wo, gradp, nfu, nbu, nrv, nlv, nbw, ntw)
    !----------------------------------------------------------------
    !   Short: Grid kernel of pres_operator.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: cS1(kk, jj, ii), cS2(kk, jj, ii), cS3(kk, jj, ii)
        REAL(realk), INTENT(in) :: p(kk, jj, ii)
        REAL(realk), INTENT(in) :: dBa(kk, jj, ii), dLe(kk, jj, ii), dTo(kk, jj, ii)
        REAL(realk), INTENT(in) :: rdx(ii), rdy(jj), rdz(kk)
        REAL(realk), INTENT(inout) :: uo(kk, jj, ii), vo(kk, jj, ii), wo(kk, jj, ii)
        REAL(realk), INTENT(in) :: gradp(3)
        INTEGER(intk), INTENT(in) :: nfu, nbu, nrv, nlv, nbw, ntw

        ! Local variables
        INTEGER(intk) :: i, j, k
        REAL(realk) :: gpx, gpy, gpz

        DO i = 3-nfu, ii-3+nbu
            DO j = 3, jj-2
                DO k = 3, kk-2
                    gpx = gradp(1)*cS1(k, j, i)
                    uo(k,j,i) = uo(k,j,i) - 1.0_realk/dBa(k,j,i)*((p(k,j,i+1) - p(k,j,i))*rdx(i) + gpx)
                END DO
            END DO
        END DO

        DO i = 3, ii-2
            DO j = 3-nrv, jj-3+nlv
                DO k = 3, kk-2
                    gpy = gradp(2)*cS2(k, j, i)
                    vo(k,j,i) = vo(k,j,i) - 1.0_realk/dLe(k,j,i)*((p(k,j+1,i) - p(k,j,i))*rdy(j) + gpy)
                END DO
            END DO
        END DO

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3-nbw, kk-3+ntw
                    gpz = gradp(3)*cS3(k, j, i)
                    wo(k,j,i) = wo(k,j,i) - 1.0_realk/dTo(k,j,i)*((p(k+1,j,i) - p(k,j,i))*rdz(k) + gpz)
                END DO
            END DO
        END DO

    END SUBROUTINE pres_operator_grd

    !================================================================

    SUBROUTINE exte_operator(uo_f, vo_f, wo_f)
    !----------------------------------------------------------------
    !   Short: Add gravitational acceleration to the momentum RHS.
    !
    !   Skipped if skpExt.
    !----------------------------------------------------------------

        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: uo_f, vo_f, wo_f

        ! Local variables
        INTEGER(intk) :: n, igrid
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: uo(:,:,:), vo(:,:,:), wo(:,:,:)
        INTEGER(intk) :: nfro, nbac, nrgt, nlft, nbot, ntop
        INTEGER(intk) :: nfu, nbu, nrv, nlv, nbw, ntw

        IF ( skpExt ) RETURN

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL uo_f%get_ptr(uo, igrid)
            CALL vo_f%get_ptr(vo, igrid)
            CALL wo_f%get_ptr(wo, igrid)

            CALL get_lp_mdf(igrid, nfro, nbac, nrgt, nlft, nbot, ntop, &
                nfu, nbu, nrv, nlv, nbw, ntw)
            CALL exte_operator_grd(kk, jj, ii, uo, vo, wo, &
                nfu, nbu, nrv, nlv, nbw, ntw)
        END DO

    END SUBROUTINE exte_operator

    !================================================================

    SUBROUTINE exte_operator_grd(kk, jj, ii, uo, vo, wo, &
        nfu, nbu, nrv, nlv, nbw, ntw)
    !----------------------------------------------------------------
    !   Short: Grid kernel of exte_operator.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: uo(kk, jj, ii), vo(kk, jj, ii), wo(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: nfu, nbu, nrv, nlv, nbw, ntw

        ! Local variables
        INTEGER(intk) :: i, j, k

        DO i = 3-nfu, ii-3+nbu
            DO j = 3, jj-2
                DO k = 3, kk-2
                    uo(k,j,i) = uo(k,j,i) + grav(1)
                ENDDO
            ENDDO
        ENDDO

        DO i = 3, ii-2
            DO j = 3-nrv, jj-3+nlv
                DO k = 3, kk-2
                    vo(k,j,i) = vo(k,j,i) + grav(2)
                ENDDO
            ENDDO
        ENDDO

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3-nbw, kk-3+ntw
                    wo(k,j,i) = wo(k,j,i) + grav(3)
                ENDDO
            ENDDO
        ENDDO

    END SUBROUTINE exte_operator_grd

    !================================================================

    SUBROUTINE adv_c(l, dt)
    !----------------------------------------------------------------
    !   Short: Advect c on the main grid in split direction l.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: l
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        INTEGER(intk) :: n, igrid
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: c(:,:,:), cWy(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: up(:,:,:), vp(:,:,:), wp(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: vel(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: cFlx1(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        INTEGER(intk) :: ista, iend, jsta, jend, ksta, kend

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_fieldptr(c, "C", igrid)
            CALL get_fieldptr(cWy, "CWY", igrid)
            CALL get_fieldptr(up, "UP", igrid)
            CALL get_fieldptr(vp, "VP", igrid)
            CALL get_fieldptr(wp, "WP", igrid)
            CALL get_fieldptr(cFlx1, "CFLX1", igrid)
            CALL get_fieldptr(dx, "DX", igrid)
            CALL get_fieldptr(dy, "DY", igrid)
            CALL get_fieldptr(dz, "DZ", igrid)
            CALL get_fieldptr(ddx, "DDX", igrid)
            CALL get_fieldptr(ddy, "DDY", igrid)
            CALL get_fieldptr(ddz, "DDZ", igrid)

            CALL sel_vel(l, up, vp, wp, vel)
            CALL get_lp_bnds(igrid, 0, ista, iend, jsta, jend, ksta, kend)
            CALL adv_c_grd(kk, jj, ii, 0, l, c, cWy, vel, cFlx1, &
                dx, dy, dz, ddx, ddy, ddz, dt, &
                ista, iend, jsta, jend, ksta, kend)
        END DO

    END SUBROUTINE adv_c

    !================================================================

    SUBROUTINE adv_c_stg(q, l, dt)
    !----------------------------------------------------------------
    !   Short: Advect c on staggered grid q in split direction l;
    !          clip it and correct the momentum by the clipped mass.
    !
    !   After the flux divergence (adv_c_grd), cS_q is clipped
    !   with clp_tol. The mass added or removed by the clipping,
    !   (rho1 - rho2)*dc, is added to the staggered momentum
    !   with the current velocity u_q, so that mass and momentum
    !   stay consistent.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: q, l
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        CHARACTER(len=5) :: cFldName, cWyFldName, mFldName
        INTEGER(intk) :: n, igrid
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: cSq(:,:,:), cWySq(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: advr(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: cFlx1(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: mSq(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: u(:,:,:), v(:,:,:), w(:,:,:), vel(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        INTEGER(intk) :: ista, iend, jsta, jend, ksta, kend

        cFldName = "CS"//int2Char(q)
        cWyFldName = "CWYS"//int2Char(q)
        mFldName = "MS"//int2Char(q)

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_fieldptr(cSq, cFldName, igrid)
            CALL get_fieldptr(cWySq, cWyFldName, igrid)
            CALL get_fieldptr(advr, "ADVR", igrid)
            CALL get_fieldptr(cFlx1, "CFLX1", igrid)
            CALL get_fieldptr(mSq, mFldName, igrid)
            CALL get_fieldptr(u, "U", igrid)
            CALL get_fieldptr(v, "V", igrid)
            CALL get_fieldptr(w, "W", igrid)
            CALL get_fieldptr(dx, "DX", igrid)
            CALL get_fieldptr(dy, "DY", igrid)
            CALL get_fieldptr(dz, "DZ", igrid)
            CALL get_fieldptr(ddx, "DDX", igrid)
            CALL get_fieldptr(ddy, "DDY", igrid)
            CALL get_fieldptr(ddz, "DDZ", igrid)

            CALL get_lp_bnds(igrid, q, ista, iend, jsta, jend, ksta, kend)
            CALL adv_c_grd(kk, jj, ii, q, l, cSq, cWySq, advr, cFlx1, &
                dx, dy, dz, ddx, ddy, ddz, dt, ista, iend, jsta, jend, ksta, kend)

            ! Clip cS and correct mS by clipped amount
            CALL sel_vel(q, u, v, w, vel)
            mSq = mSq + (rho1 - rho2)*( clp_tol(cSq) - cSq )*vel
            cSq = clp_tol(cSq)
        END DO

    END SUBROUTINE adv_c_stg

    !================================================================

    SUBROUTINE adv_c_grd(kk, jj, ii, q, l, c, cWy, vel, cFlx1, &
        dx, dy, dz, ddx, ddy, ddz, dt, ista, iend, jsta, jend, ksta, kend)
    !----------------------------------------------------------------
    !   Short: Grid kernel: flux divergence plus Weymouth-Yue
    !          compression term.
    !
    !   c = c - dt/ds*(F(+) - F(-)) + dt*cWy*du_l/dx_l with the
    !   extent ds of the control volume in direction l (sel_ext) and
    !   the velocity vel at its faces. The compression term with the
    !   coefficient cWy, fixed at the beginning of the time step,
    !   makes the directionally split scheme conservative (see
    !   comp_cWy_grd). The main grid (q = 0) is clipped with clp_tol
    !   here, the staggered grids in adv_c_stg.
    !
    !   Source:
    !   G. D. Weymouth and D. K.-P. Yue, “Conservative
    !   Volume-of-Fluid method for free-surface simulations on
    !   Cartesian-grids,” Journal of Computational Physics, vol. 229,
    !   no. 8, pp. 2853–2865, Apr. 2010,
    !   doi: 10.1016/j.jcp.2009.12.018.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, q, l
        REAL(realk), INTENT(inout) :: c(kk, jj, ii)
        REAL(realk), INTENT(in) :: cWy(kk, jj, ii)
        REAL(realk), INTENT(in) :: vel(kk,jj,ii)
        REAL(realk), INTENT(in) :: cFlx1(kk, jj, ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: dt
        INTEGER(intk), INTENT(in) :: ista, iend, jsta, jend, ksta, kend

        ! Local variables
        INTEGER(intk) :: i, j, k
        INTEGER(intk) :: il, jl, kl
        REAL(realk) :: dsx(ii), dsy(jj), dsz(kk), dsCV
        REAL(realk) :: div

        CALL sel_ind(l, il, jl, kl)
        CALL sel_ext(kk, jj, ii, q, l, dx, dy, dz, ddx, ddy, ddz, dsx, dsy, dsz)

        DO i = ista, iend
            DO j = jsta, jend
                DO k = ksta, kend
                    dsCV = il*dsx(i) + jl*dsy(j) + kl*dsz(k)
                    div = ( vel(k,j,i) - vel(k-kl,j-jl,i-il) ) / dsCV
                    c(k,j,i) = c(k,j,i) &
                        - dt/dsCV*( cFlx1(k,j,i) - cFlx1(k-kl,j-jl,i-il) ) &
                        + dt*cWy(k,j,i)*div
                END DO
            END DO
        END DO

        IF ( q == 0 ) c = clp_tol(c)

    END SUBROUTINE adv_c_grd

    !================================================================

    SUBROUTINE adv_m_stg(q, l, dt)
    !----------------------------------------------------------------
    !   Short: Advect momentum on staggered grid q in split
    !          direction l.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: q, l
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        CHARACTER(len=5) :: cWyFldName, mFldName
        INTEGER(intk) :: n, igrid
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: cWySq(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: u(:,:,:), v(:,:,:), w(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: up(:,:,:), vp(:,:,:), wp(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: vel(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: advr(:,:,:), adve(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: mSq(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: cFlx1(:,:,:), cFlx2(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        INTEGER(intk) :: ista, iend, jsta, jend, ksta, kend

        cWyFldName = "CWYS"//int2Char(q)
        mFldName = "MS"//int2Char(q)

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_fieldptr(cWySq, cWyFldName, igrid)
            CALL get_fieldptr(u, "U", igrid)
            CALL get_fieldptr(v, "V", igrid)
            CALL get_fieldptr(w, "W", igrid)
            CALL get_fieldptr(up, "UP", igrid)
            CALL get_fieldptr(vp, "VP", igrid)
            CALL get_fieldptr(wp, "WP", igrid)
            CALL get_fieldptr(advr, "ADVR", igrid)
            CALL get_fieldptr(adve, "ADVE", igrid)
            CALL get_fieldptr(mSq, mFldName, igrid)
            CALL get_fieldptr(cFlx1, "CFLX1", igrid)
            CALL get_fieldptr(cFlx2, "CFLX2", igrid)
            CALL get_fieldptr(dx, "DX", igrid)
            CALL get_fieldptr(dy, "DY", igrid)
            CALL get_fieldptr(dz, "DZ", igrid)
            CALL get_fieldptr(ddx, "DDX", igrid)
            CALL get_fieldptr(ddy, "DDY", igrid)
            CALL get_fieldptr(ddz, "DDZ", igrid)

            CALL sel_vel(q, up, vp, wp, vel)
            CALL get_lp_bnds(igrid, q, ista, iend, jsta, jend, ksta, kend)
            CALL adv_m_grd(kk, jj, ii, q, l, cWySq, vel, &
                advr, adve, mSq, cFlx1, cFlx2, dx, dy, dz, &
                ddx, ddy, ddz, dt, ista, iend, jsta, jend, ksta, kend)
        END DO

    END SUBROUTINE adv_m_stg

    !================================================================

    SUBROUTINE adv_m_grd(kk, jj, ii, q, l, cWy, vel, &
        advr, adve, mom, cFlx1, cFlx2, dx, dy, dz, &
        ddx, ddy, ddz, dt, ista, iend, jsta, jend, ksta, kend)
    !----------------------------------------------------------------
    !   Short: Grid kernel: mass-consistent momentum flux divergence
    !          plus compression term.
    !
    !   The momentum flux is adve*(rho1*cFlx1 + rho2*cFlx2), i.e.
    !   it uses the same phase fluxes as the advection of cS_q,
    !   which makes mass and momentum transport consistent (Arrufat
    !   et al.). The compression term is the density of the
    !   compression coefficient, rho1*cWy + rho2*(1 - cWy), times
    !   the velocity at the beginning of the time step and the
    !   divergence of the advector.
    !
    !   Source:
    !   T. Arrufat et al., “A mass-momentum consistent,
    !   Volume-of-Fluid method for incompressible flow on staggered
    !   grids,” Computers & Fluids, vol. 215, p. 104785, Jan. 2021,
    !   doi: 10.1016/j.compfluid.2020.104785.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, q, l
        REAL(realk), INTENT(in) :: cWy(kk, jj, ii)
        REAL(realk), INTENT(in) :: vel(kk,jj,ii)
        REAL(realk), INTENT(in) :: advr(kk,jj,ii), adve(kk,jj,ii)
        REAL(realk), INTENT(inout) :: mom(kk, jj, ii)
        REAL(realk), INTENT(in) :: cFlx1(kk, jj, ii), cFlx2(kk, jj, ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: dt
        INTEGER(intk), INTENT(in) :: ista, iend, jsta, jend, ksta, kend

        ! Local variables
        INTEGER(intk) :: i, j, k
        INTEGER(intk) :: il, jl, kl
        REAL(realk) :: dsx(ii), dsy(jj), dsz(kk), dsCV
        REAL(realk) :: momFlx(kk, jj, ii)
        REAL(realk) :: div, com

        CALL sel_ind(l, il, jl, kl)
        CALL sel_ext(kk, jj, ii, q, l, dx, dy, dz, ddx, ddy, ddz, dsx, dsy, dsz)

        DO i = ista-il, iend
            DO j = jsta-jl, jend
                DO k = ksta-kl, kend
                    momFlx(k,j,i) = adve(k,j,i)*( rho1*cFlx1(k,j,i) + rho2*cFlx2(k,j,i) )
                END DO
            END DO
        END DO

        DO i = ista, iend
            DO j = jsta, jend
                DO k = ksta, kend
                    dsCV = il*dsx(i) + jl*dsy(j) + kl*dsz(k)
                    div = ( advr(k,j,i) - advr(k-kl,j-jl,i-il) ) / dsCV
                    com = ( rho1*cWy(k,j,i) + rho2*( 1.0_realk - cWy(k,j,i) ) )*div

                    mom(k,j,i) = mom(k,j,i) &
                        - dt/dsCV*( momFlx(k,j,i) - momFlx(k-kl,j-jl,i-il) ) &
                        + dt*vel(k,j,i)*com
                END DO
            END DO
        END DO

    END SUBROUTINE adv_m_grd

    !================================================================

    SUBROUTINE comp_m_stg(q)
    !----------------------------------------------------------------
    !   Short: Compute momentum rho*u on staggered grid q from
    !          density and velocity.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: q

        ! Local variables
        CHARACTER(len=3) :: dFldName, mFldName
        INTEGER(intk) :: n, igrid
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: u(:,:,:), v(:,:,:), w(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: vel(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: dSq(:,:,:), mSq(:,:,:)

        dFldName = "DS"//int2Char(q)
        mFldName = "MS"//int2Char(q)

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_fieldptr(dSq, dFldName, igrid)
            CALL get_fieldptr(u, "U", igrid)
            CALL get_fieldptr(v, "V", igrid)
            CALL get_fieldptr(w, "W", igrid)
            CALL get_fieldptr(mSq, mFldName, igrid)

            CALL sel_vel(q, u, v, w, vel)
            CALL comp_m_stg_grd(kk, jj, ii, q, dSq, vel, mSq)
        END DO

    END SUBROUTINE comp_m_stg

    !================================================================

    SUBROUTINE comp_m_stg_grd(kk, jj, ii, q, dSq, vel, mSq)
    !----------------------------------------------------------------
    !   Short: Grid kernel of comp_m_stg.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, q
        REAL(realk), INTENT(in) :: vel(kk, jj, ii)
        REAL(realk), INTENT(in) :: dSq(kk, jj, ii)
        REAL(realk), INTENT(out) :: mSq(kk,jj,ii)

        ! Local variables
        ! None

        mSq = vel*dSq

    END SUBROUTINE comp_m_stg_grd

    !================================================================

    SUBROUTINE upd_vel_stg(q)
    !----------------------------------------------------------------
    !   Short: Update velocity component q from advected momentum
    !          and density.
    !
    !   u_q = mS_q/dS_q. Only the unknowns owned by the grid are
    !   updated, i.e. with the same loop bounds as in adv_m_stg (see
    !   get_lp_bnds). Faces on walls and refinement boundaries as
    !   well as ghost layers are not touched here; the refinement
    !   boundaries and ghost layers are filled afterwards by
    !   prlg_stg("U", "V", "W") in adve_operator.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: q

        ! Local variables
        CHARACTER(len=3) :: dFldName, mFldName
        INTEGER(intk) :: n, igrid
        INTEGER(intk) :: kk, jj, ii
        INTEGER(intk) :: ista, iend, jsta, jend, ksta, kend
        REAL(realk), POINTER, CONTIGUOUS :: u(:,:,:), v(:,:,:), w(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: vel(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: dSq(:,:,:), mSq(:,:,:)

        dFldName = "DS"//int2char(q)
        mFldName = "MS"//int2char(q)

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_fieldptr(dSq, dFldName, igrid)
            CALL get_fieldptr(mSq, mFldName, igrid)
            CALL get_fieldptr(u, "U", igrid)
            CALL get_fieldptr(v, "V", igrid)
            CALL get_fieldptr(w, "W", igrid)

            CALL sel_vel(q, u, v, w, vel)
            CALL get_lp_bnds(igrid, q, ista, iend, jsta, jend, ksta, kend)
            CALL upd_vel_stg_grd(kk, jj, ii, dSq, mSq, vel, &
                ista, iend, jsta, jend, ksta, kend)
        END DO

    END SUBROUTINE upd_vel_stg

    !================================================================

    SUBROUTINE upd_vel_stg_grd(kk, jj, ii, dSq, mSq, vel, &
        ista, iend, jsta, jend, ksta, kend)
    !----------------------------------------------------------------
    !   Short: Grid kernel: u_q = (rho*u)_q/rho_q on owned unknowns.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: dSq(kk, jj, ii), mSq(kk, jj, ii)
        REAL(realk), INTENT(inout) :: vel(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: ista, iend, jsta, jend, ksta, kend

        ! Local variables
        INTEGER(intk) :: i, j, k

        DO i = ista, iend
            DO j = jsta, jend
                DO k = ksta, kend
                    vel(k, j, i) = mSq(k, j, i)/dSq(k, j, i)
                END DO
            END DO
        END DO

    END SUBROUTINE upd_vel_stg_grd

    !================================================================

    SUBROUTINE comp_cWy()
    !----------------------------------------------------------------
    !   Short: Compute Weymouth-Yue compression coefficient on the
    !          main grid.
    !----------------------------------------------------------------

        ! Subroutine arguments
        ! None

        ! Local variables
        INTEGER(intk) :: n, igrid
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: c(:,:,:), cWy(:,:,:)

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_fieldptr(c, "C", igrid)
            CALL get_fieldptr(cWy, "CWY", igrid)

            CALL comp_cWy_grd(kk, jj, ii, c, cWy)
        END DO

    END SUBROUTINE comp_cWy

    !================================================================

    SUBROUTINE comp_cWy_stg(q)
    !----------------------------------------------------------------
    !   Short: Compute Weymouth-Yue compression coefficient on
    !          staggered grid q.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: q

        ! Local variables
        CHARACTER(len=5) :: cFldName, cWyFldName
        INTEGER(intk) :: n, igrid
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: cSq(:,:,:), cWySq(:,:,:)

        cFldName = "CS"//int2Char(q)
        cWyFldName = "CWYS"//int2Char(q)

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_fieldptr(cSq, cFldName, igrid)
            CALL get_fieldptr(cWySq, cWyFldName, igrid)

            CALL comp_cWy_grd(kk, jj, ii, cSq, cWySq)
        END DO

    END SUBROUTINE comp_cWy_stg

    !================================================================

    SUBROUTINE comp_cWy_grd(kk, jj, ii, c, cWy)
    !----------------------------------------------------------------
    !   Short: Grid kernel: coefficient is 1 for c > 0.5, else 0.
    !
    !   The coefficient is non-directional. It is evaluated once per
    !   time step before the first split and kept constant for all
    !   splits, which makes the directionally split advection
    !   exactly conservative.
    !
    !   Sources:
    !   T. Arrufat et al., “A mass-momentum consistent,
    !   Volume-of-Fluid method for incompressible flow on staggered
    !   grids,” Computers & Fluids, vol. 215, p. 104785, Jan. 2021,
    !   doi: 10.1016/j.compfluid.2020.104785.
    !
    !   G. D. Weymouth and D. K.-P. Yue, “Conservative
    !   Volume-of-Fluid method for free-surface simulations on
    !   Cartesian-grids,” Journal of Computational Physics, vol. 229,
    !   no. 8, pp. 2853–2865, Apr. 2010,
    !   doi: 10.1016/j.jcp.2009.12.018.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: c(kk, jj, ii)
        REAL(realk), INTENT(inout) :: cWy(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    IF ( c(k,j,i) > 0.5_realk ) THEN
                        cWy(k,j,i) = 1.0_realk
                    ELSE
                        cWy(k,j,i) = 0.0_realk
                    END IF
                END DO
            END DO
        END DO

    END SUBROUTINE comp_cWy_grd

    !================================================================

    SUBROUTINE comp_adve_stg(q, l, dt)
    !----------------------------------------------------------------
    !   Short: Compute advected velocity (advectee) on staggered
    !          grid q.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: q, l
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        CHARACTER(len=10) :: isIfcVicFldName
        INTEGER(intk) :: n, igrid
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: u(:,:,:), v(:,:,:), w(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: vel(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: advr(:,:,:), adve(:,:,:)
        REAl(realk), POINTER, CONTIGUOUS :: isIfcVicSq(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        INTEGER(intk) :: ista, iend, jsta, jend, ksta, kend

        isIfcVicFldName = "ISIFCVICS"//int2char(q)

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_fieldptr(u, "U", igrid)
            CALL get_fieldptr(v, "V", igrid)
            CALL get_fieldptr(w, "W", igrid)
            CALL get_fieldptr(advr, "ADVR", igrid)
            CALL get_fieldptr(adve, "ADVE", igrid)
            CALL get_fieldptr(isIfcVicSq, isIfcVicFldName, igrid)
            CALL get_fieldptr(dx, "DX", igrid)
            CALL get_fieldptr(dy, "DY", igrid)
            CALL get_fieldptr(dz, "DZ", igrid)
            CALL get_fieldptr(ddx, "DDX", igrid)
            CALL get_fieldptr(ddy, "DDY", igrid)
            CALL get_fieldptr(ddz, "DDZ", igrid)

            CALL sel_vel(q, u, v, w, vel)
            CALL get_lp_bnds(igrid, q, ista, iend, jsta, jend, ksta, kend)
            CALL comp_adve_stg_grd(kk, jj, ii, q, l, vel, advr, adve, &
                isIfcVicSq, dx, dy, dz, ddx, ddy, ddz, dt, &
                ista, iend, jsta, jend, ksta, kend)
        END DO

    END SUBROUTINE comp_adve_stg

    !================================================================

    SUBROUTINE comp_adve_stg_grd(kk, jj, ii, q, l, vel, advr, adve, &
        isIfcVic, dx, dy, dz, ddx, ddy, ddz, dt, &
        ista, iend, jsta, jend, ksta, kend)
    !----------------------------------------------------------------
    !   Short: Grid kernel: select interpolation scheme (QUICK or
    !          ENO) for the advectee.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, q, l
        REAL(realk), INTENT(in) :: vel(kk, jj, ii)
        REAL(realk), INTENT(in) :: advr(kk, jj, ii)
        REAL(realk), INTENT(inout) :: adve(kk, jj, ii)
        REAL(realk), INTENT(in) :: isIfcVic(kk, jj, ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: dt
        INTEGER(intk), INTENT(in) :: ista, iend, jsta, jend, ksta, kend

        ! Local variables
        ! None

        IF ( advScm == "QUICK" ) THEN
            CALL comp_adve_stg_quick(kk, jj, ii, q, l, vel, advr, &
                adve, isIfcVic, dx, dy, dz, ddx, ddy, ddz, dt, &
                ista, iend, jsta, jend, ksta, kend)
        ELSEIF ( advScm == "ENO" ) THEN
            CALL comp_adve_stg_eno(kk, jj, ii, q, l, vel, advr, &
                adve, dx, dy, dz, ddx, ddy, ddz, dt, &
                ista, iend, jsta, jend, ksta, kend)
        ELSE
            CALL err_abort(vofErr, "Unknown interpolation.", __FILE__, __LINE__)
        ENDIF

    END SUBROUTINE comp_adve_stg_grd

    !================================================================

    SUBROUTINE comp_adve_stg_quick(kk, jj, ii, q, l, vel, &
        advr, adve, isIfcVic, dx, dy, dz, ddx, ddy, ddz, dt, &
        ista, iend, jsta, jend, ksta, kend)
    !----------------------------------------------------------------
    !   Short: Interpolate the advectee with QUICK; upwind value near
    !          the interface.
    !
    !   Interpolates the advected velocity component q (advectee,
    !   adve) to the faces of the staggered q-cells in direction l;
    !   the upwind side follows from the advecting velocity
    !   (advector, advr). The quadratic upstream interpolation is
    !   evaluated in Newton form on non-uniform grids. In the
    !   vicinity of the interface (isIfcVic > 0) the first-order
    !   upwind value is used.
    !
    !   With donCen, the interpolation point is moved by |advr|*dt/2
    !   upstream, i.e. to the centre of the donating region.
    !   Indicator functions are used to avoid if-statements within
    !   the loops.
    !
    !   Sources:
    !   W. Aniszewski et al., “PArallel, Robust, Interface Simulator
    !   (PARIS),” Computer Physics Communications, vol. 263,
    !   p. 107849, Jun. 2021, doi: 10.1016/j.cpc.2021.107849.
    !   See PARIS (accessed: May 2026) source code subroutine(s):
    !       - interpole_quad (rudman_modules.f90)
    !
    !   T. Arrufat et al., “A mass-momentum consistent,
    !   Volume-of-Fluid method for incompressible flow on staggered
    !   grids,” Computers & Fluids, vol. 215, p. 104785, Jan. 2021,
    !   doi: 10.1016/j.compfluid.2020.104785.
    !
    !   D. Herrmann, Numerische Mathematik — 40 BASIC-Programme.
    !   Wiesbaden: Vieweg+Teubner Verlag, 1983.
    !   doi: 10.1007/978-3-322-96321-5.
    !
    !   B. P. Leonard, “A stable and accurate convective modelling
    !   procedure based on quadratic upstream interpolation,”
    !   Computer Methods in Applied Mechanics and Engineering,
    !   vol. 19, no. 1, pp. 59–98, Jun. 1979,
    !   doi: 10.1016/0045-7825(79)90034-3.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, q, l
        REAL(realk), INTENT(in) :: vel(kk, jj, ii)
        REAL(realk), INTENT(in) :: advr(kk, jj, ii)
        REAL(realk), INTENT(inout) :: adve(kk, jj, ii)
        REAL(realk), INTENT(in) :: isIfcVic(kk, jj, ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: dt
        INTEGER(intk), INTENT(in) :: ista, iend, jsta, jend, ksta, kend

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: kl, jl, il
        REAL(realk) :: dnx(ii), dny(jj), dnz(kk)
        REAL(realk) :: dcx(ii), dcy(jj), dcz(kk)
        REAL(realk) :: signInd(2), iFacInd(2)
        REAL(realk) :: dnslL, dnslC, dnslR, dcslMi, dcslPl
        REAL(realk) :: adveMi, advePl
        REAL(realk) :: velMi, velCe, velPl, velFP

        CALL sel_ind(l, il, jl, kl)

        IF ( q == l ) THEN
            dnx(1:ii-1) = ddx(2:ii) ; dnx(ii) = 0.0_realk
            dny(1:jj-1) = ddy(2:jj) ; dny(jj) = 0.0_realk
            dnz(1:kk-1) = ddz(2:kk) ; dnz(kk) = 0.0_realk
            dcx = dnx ; dcy = dny ; dcz = dnz
        ELSE
            dnx = dx  ; dny = dy  ; dnz = dz
            dcx = ddx ; dcy = ddy ; dcz = ddz
        ENDIF

        DO i = ista-il, iend
            DO j = jsta-jl, jend
                DO k = ksta-kl, kend
                    signInd(1) = MERGE(1.0_realk, 0.0_realk, advr(k,j,i) >= 0.0_realk)
                    signInd(2) = 1.0_realk - signInd(1)
                    iFacInd(1) = MERGE(1.0_realk, 0.0_realk, isIfcVic(k,j,i) > 0.0_realk)
                    iFacInd(2) = 1.0_realk - iFacInd(1)

                    dnslL = il*dnx(i-il) + jl*dny(j-jl) + kl*dnz(k-kl)
                    dnslC = il*dnx(i) + jl*dny(j) + kl*dnz(k)
                    dnslR = il*dnx(i+il) + jl*dny(j+jl) + kl*dnz(k+kl)

                    dcslMi = 0.5_realk*(il*dcx(i) + jl*dcy(j) + kl*dcz(k))
                    dcslPl  = dnslC - dcslMi

                    IF ( donCen ) THEN
                        dcslMi = dcslMi - ABS(advr(k,j,i))*dt/2.0_realk
                        dcslPl = dcslPl - ABS(advr(k,j,i))*dt/2.0_realk
                    ENDIF

                    velMi = vel(k-kl,j-jl,i-il)
                    velCe = vel(k,j,i)
                    velPl = vel(k+kl,j+jl,i+il)
                    velFP = vel(k+2*kl,j+2*jl,i+2*il)

                    adveMi = newton_interpolation(velMi, velCe, velPl, dnslL, dnslC, dcslMi)
                    advePl = newton_interpolation(velFP, velPl, velCe, dnslR, dnslC, dcslPl)

                    adve(k,j,i) = iFacInd(1)*(signInd(1)*velCe + signInd(2)*velPl) + &
                                  iFacInd(2)*(signInd(1)*adveMi + signInd(2)*advePl)
                END DO
            END DO
        END DO

    CONTAINS

        !------------------------------------------------------------
        ! Both the Lagrange and the Newton form yield the same
        ! results. The Lagrange form was more familiar to me, hence
        ! I implemented it first. The Newton form is closer to the
        ! formulation of Leonard and easier to adjust in the future
        ! (see QUICKEST scheme).
        !------------------------------------------------------------
        ! PURE REAL(realk) FUNCTION lagrange_interpolation(phiFU, phiU, phiD, dsFU, dsD, dsr) RESULT(r)
        !     REAL(realk), INTENT(in) :: phiFU, phiU, phiD, dsFU, dsD, dsr
        !     r = phiFU*(dsr*(dsr-dsD))/(dsFU*(dsFU+dsD)) &
        !         - phiU*((dsr+dsFU)*(dsr-dsD))/(dsFU*dsD) &
        !         + phiD*(dsr*(dsr+dsFU))/(dsD*(dsFU+dsD))
        ! END FUNCTION lagrange_interpolation

        PURE REAL(realk) FUNCTION newton_interpolation(phiFU, phiU, phiD, dsFU, dsD, dsr) RESULT(r)
        !------------------------------------------------------------
        !   Short: Quadratic upstream interpolation (Newton form) at
        !          distance dsr from the upwind node.
        !------------------------------------------------------------
            REAL(realk), INTENT(in) :: phiFU, phiU, phiD, dsFU, dsD, dsr
            REAL(realk) :: dd1U, dd1D, curv
            dd1D = (phiD - phiU)/dsD
            dd1U = (phiU - phiFU)/dsFU
            curv = (dd1D - dd1U)/(dsFU + dsD)
            r = phiU + dsr*dd1D + dsr*(dsr - dsD)*curv
        END FUNCTION newton_interpolation

    END SUBROUTINE comp_adve_stg_quick

    !================================================================

    SUBROUTINE comp_adve_stg_eno(kk, jj, ii, q, l, vel, &
        advr, adve, dx, dy, dz, ddx, ddy, ddz, dt, &
        ista, iend, jsta, jend, ksta, kend)
    !----------------------------------------------------------------
    !   Short: Interpolate the advectee with ENO
    !          (minmod-limited slope).
    !
    !   Second-order upwind extrapolation of the advectee from the
    !   upwind node with the minmod-limited slope at this node.
    !   donCen as in comp_adve_stg_quick. Indicator functions are
    !   used to avoid if-statements within the loops.
    !
    !   Sources:
    !   W. Aniszewski et al., “PArallel, Robust, Interface Simulator
    !   (PARIS),” Computer Physics Communications, vol. 263,
    !   p. 107849, Jun. 2021, doi: 10.1016/j.cpc.2021.107849.
    !   See PARIS (accessed: May 2026) source code subroutine(s):
    !       - interpole3 (modules.f90)
    !       - slope_lim (modules.f90)
    !
    !   P. K. Sweby, “High Resolution Schemes Using Flux Limiters for
    !   Hyperbolic Conservation Laws,” SIAM Journal on Numerical
    !   Analysis, vol. 21, no. 5, pp. 995–1011, Oct. 1984,
    !   doi: 10.1137/0721062.
    !
    !   G. Tryggvason, R. Scardovelli, and S. Zaleski, Direct
    !   Numerical Simulations of Gas–Liquid Multiphase Flows, 1st ed.
    !   Cambridge University Press, 2011.
    !   doi: 10.1017/CBO9780511975264.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, q, l
        REAL(realk), INTENT(in) :: vel(kk, jj, ii)
        REAL(realk), INTENT(in) :: advr(kk, jj, ii)
        REAL(realk), INTENT(inout) :: adve(kk, jj, ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: dt
        INTEGER(intk), INTENT(in) :: ista, iend, jsta, jend, ksta, kend

        ! Local variables
        INTEGER(intk) :: k, j, i, kl, jl, il
        REAL(realk) :: dnx(ii), dny(jj), dnz(kk)
        REAL(realk) :: dcx(ii), dcy(jj), dcz(kk)
        REAL(realk) :: signInd(2)
        REAL(realk) :: dnslMi, dnslCe, dnslPl, dcslMi, dcslPl
        REAL(realk) :: slopeMi, slopeCe, slopePl, s

        CALL sel_ind(l, il, jl, kl)

        IF ( q == l ) THEN
            dnx(1:ii-1) = ddx(2:ii) ; dnx(ii) = 0.0_realk
            dny(1:jj-1) = ddy(2:jj) ; dny(jj) = 0.0_realk
            dnz(1:kk-1) = ddz(2:kk) ; dnz(kk) = 0.0_realk
            dcx = dnx ; dcy = dny ; dcz = dnz
        ELSE
            dnx = dx  ; dny = dy  ; dnz = dz
            dcx = ddx ; dcy = ddy ; dcz = ddz
        ENDIF

        DO i = ista-il, iend
            DO j = jsta-jl, jend
                DO k = ksta-kl, kend
                    signInd(1) = MERGE(1.0_realk, 0.0_realk, advr(k,j,i) >= 0.0_realk)
                    signInd(2) = 1.0_realk - signInd(1)

                    dnslMi = il*dnx(i-1) + jl*dny(j-1) + kl*dnz(k-1)
                    dnslCe = il*dnx(i) + jl*dny(j) + kl*dnz(k)
                    dnslPl = il*dnx(i+1) + jl*dny(j+1) + kl*dnz(k+1)

                    slopeMi = (vel(k,j,i) - vel(k-kl,j-jl,i-il))/dnslMi
                    slopeCe = (vel(k+kl,j+jl,i+il) - vel(k,j,i))/dnslCe
                    slopePl = (vel(k+2*kl,j+2*jl,i+2*il) - vel(k+kl,j+jl,i+il))/dnslPl

                    s = signInd(1)*minmod(slopeMi, slopeCe) + &
                        signInd(2)*minmod(slopeCe, slopePl)

                    dcslMi = 0.5_realk*(il*dcx(i) + jl*dcy(j) + kl*dcz(k))
                    dcslPl = dnslCe - dcslMi

                    IF ( donCen ) THEN
                        dcslMi = dcslMi - ABS(advr(k,j,i))*dt/2.0_realk
                        dcslPl = dcslPl - ABS(advr(k,j,i))*dt/2.0_realk
                    ENDIF

                    adve(k,j,i) = signInd(1)*(vel(k,j,i) + s*dcslMi) + &
                                  signInd(2)*(vel(k+kl,j+jl,i+il) - s*dcslPl)
                END DO
            END DO
        END DO

    CONTAINS

        PURE REAL(realk) FUNCTION minmod(a, b) RESULT(res)
        !------------------------------------------------------------
        !   Short: Minmod limiter.
        !------------------------------------------------------------
            REAL(realk), INTENT(in) :: a, b
            res = 0.5_realk*( SIGN(1.0_realk, a) + SIGN(1.0_realk, b) )*MIN(ABS(a), ABS(b))
        END FUNCTION minmod

    END SUBROUTINE comp_adve_stg_eno

    !================================================================

    SUBROUTINE comp_advr_stg(q, l)
    !----------------------------------------------------------------
    !   Short: Compute advecting velocity (advector) on staggered
    !          grid q.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: q, l

        ! Local variables
        INTEGER(intk) :: n, igrid
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: up(:,:,:), vp(:,:,:), wp(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: vel(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: advr(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        INTEGER(intk) :: ista, iend, jsta, jend, ksta, kend

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_fieldptr(up, "UP", igrid)
            CALL get_fieldptr(vp, "VP", igrid)
            CALL get_fieldptr(wp, "WP", igrid)
            CALL get_fieldptr(advr, "ADVR", igrid)
            CALL get_fieldptr(ddx, "DDX", igrid)
            CALL get_fieldptr(ddy, "DDY", igrid)
            CALL get_fieldptr(ddz, "DDZ", igrid)

            CALL sel_vel(l, up, vp, wp, vel)
            CALL get_lp_bnds(igrid, q, ista, iend, jsta, jend, ksta, kend)
            CALL comp_advr_stg_grd(kk, jj, ii, q, l, vel, advr, ddx, ddy, ddz, &
                ista, iend, jsta, jend, ksta, kend)
        END DO

    END SUBROUTINE comp_advr_stg

    !================================================================

    SUBROUTINE comp_advr_stg_grd(kk, jj, ii, q, l, vel, advr, ddx, ddy, ddz, &
        ista, iend, jsta, jend, ksta, kend)
    !----------------------------------------------------------------
    !   Short: Grid kernel: advector as (area-weighted) mean of the
    !          main-grid face velocities.
    !
    !   The advector is based on the velocity at the beginning of the
    !   time step (UP, VP, WP).
    !
    !   q == l: the staggered face lies at the centre of a main cell;
    !           advr is the mean of the velocities at its two faces.
    !   q /= l: the staggered face covers one half of two adjacent
    !           main faces; advr is the mean of their velocities
    !           weighted with the face areas, consistent with the
    !           half-cell fluxes in comp_flx_stg_grd.
    !
    !   Source:
    !   T. Arrufat et al., “A mass-momentum consistent,
    !   Volume-of-Fluid method for incompressible flow on staggered
    !   grids,” Computers & Fluids, vol. 215, p. 104785, Jan. 2021,
    !   doi: 10.1016/j.compfluid.2020.104785.
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, q, l
        REAL(realk), INTENT(in) :: vel(kk, jj, ii)
        REAL(realk), INTENT(inout) :: advr(kk, jj, ii)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        INTEGER(intk), INTENT(in) :: ista, iend, jsta, jend, ksta, kend

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: kl, jl, il, kq, jq, iq
        REAL(realk) :: ddnqMi, ddnqPl

        CALL sel_ind(l, il, jl, kl)
        CALL sel_ind(q, iq, jq, kq)

        IF ( q == l ) THEN
            DO i = ista-il, iend
                DO j = jsta-jl, jend
                    DO k = ksta-kl, kend
                        advr(k,j,i) = 0.5_realk*(vel(k,j,i) + vel(k+kq,j+jq,i+iq))
                    ENDDO
                ENDDO
            ENDDO
        ELSE
            DO i = ista-il, iend
                DO j = jsta-jl, jend
                    DO k = ksta-kl, kend
                        ddnqMi = iq*ddx(i) + jq*ddy(j) + kq*ddz(k)
                        ddnqPl = iq*ddx(i+iq) + jq*ddy(j+jq) + kq*ddz(k+kq)
                        advr(k,j,i) = (vel(k,j,i)*ddnqMi + vel(k+kq,j+jq,i+iq)*ddnqPl)/(ddnqMi + ddnqPl)
                    ENDDO
                ENDDO
            ENDDO
        ENDIF

    END SUBROUTINE comp_advr_stg_grd

END MODULE mph_vof_mod
