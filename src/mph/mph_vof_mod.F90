!====================================================================
!  Module: mph_vof_mod
!
!   Responsibilities:
!   - 
!
!   Author:      Valentin Ruhs
!   e-Mail:      valentin.ruhs@gmx.de
!   Created:     2026-02
!   Last update: 2026-09
!
!====================================================================

MODULE mph_vof_mod

    USE precision_mod, ONLY: intk, realk
    USE field_mod, ONLY: field_t
    USE fields_mod, ONLY: get_field, set_field, get_fieldptr
    USE grids_mod, ONLY: nmygrids, mygrids, get_mgdims, get_gradpxflag, &
        minlevel, maxlevel
    USE flowcore_mod, ONLY: gradp
    USE connect2_mod, ONLY: connect
    USE parent_mod, ONLY: parent
    USE ftoc_mod, ONLY: ftoc
    USE err_mod, ONLY: err_abort

    USE mphcore_mod, ONLY: rho1, rho2, gmol1, gmol2, grav, splPer, &
        skpAdv, skpDif, skpPre, skpExt, vofTol, advScm, donCen, volChk, &
        divChk, vofErr
    USE mph_utils_mod, ONLY: sel_ind, sel_ext, sel_vel, clp, int2char
    USE mph_plic_mod, ONLY: comp_ifc, comp_c_stg, comp_isIfc_stg, comp_c_loc
    USE mph_props_mod, ONLY: comp_prop, comp_prop_face, comp_prop_face_stg, &
        comp_d_stg

    IMPLICIT NONE(type, external)
    PRIVATE

    PUBLIC :: init_mph_vof, finish_mph_vof, cpy_flds, &
        adve_operator, diff_operator, pres_operator, exte_operator

CONTAINS

    SUBROUTINE init_mph_vof()

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

        ! Adver/adve and flux fields
        CALL set_field("ADVR", buffers=.TRUE.)
        CALL set_field("ADVE", buffers=.TRUE.)
        CALL set_field("CFLX1", buffers=.TRUE.)
        CALL set_field("CFLX2", buffers=.TRUE.)

    END SUBROUTINE init_mph_vof

    !================================================================

    SUBROUTINE finish_mph_vof()

        ! Subroutine arguments
        ! None

        ! Local variables
        ! None

        CONTINUE

    END SUBROUTINE finish_mph_vof

    !================================================================

    SUBROUTINE cpy_flds()
    !----------------------------------------------------------------
    !   What it does:
    !   Copy the current fields into the previous fields.
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
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: l
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        INTEGER(intk) :: n, igrid
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: c(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: u(:,:,:), v(:,:,:), w(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: vel(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: cFlx1(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        REAL(realk), POINTER, CONTIGUOUS :: normx(:,:,:), normy(:,:,:), normz(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: alpha(:,:,:)
        REAl(realk), POINTER, CONTIGUOUS :: isIfc(:,:,:)

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_fieldptr(c, "C", igrid)
            CALL get_fieldptr(u, "U", igrid)
            CALL get_fieldptr(v, "V", igrid)
            CALL get_fieldptr(w, "W", igrid)
            CALL get_fieldptr(cFlx1, "CFLX1", igrid)
            CALL get_fieldptr(ddx, "DDX", igrid)
            CALL get_fieldptr(ddy, "DDY", igrid)
            CALL get_fieldptr(ddz, "DDZ", igrid)
            CALL get_fieldptr(normx, "NORMX", igrid)
            CALL get_fieldptr(normy, "NORMY", igrid)
            CALL get_fieldptr(normz, "NORMZ", igrid)
            CALL get_fieldptr(alpha, "ALPHA", igrid)
            CALL get_fieldptr(isIfc, "ISIFC", igrid)

            CALL sel_vel(l, u, v, w, vel)
            CALL comp_flx_grd(kk, jj, ii, l, c, vel, cFlx1, &
                ddx, ddy, ddz, normx, normy, normz, alpha, isIfc, dt)
        END DO

    END SUBROUTINE comp_flx

    !================================================================

    SUBROUTINE comp_flx_grd(kk, jj, ii, l, c, vel, cFlx1, &
        ddx, ddy, ddz, normx, normy, normz, alpha, isIfc, dt)
    !----------------------------------------------------------------
    !   What it does:
    !   Computes the volume fraction flxes depending on the current
    !   split direction. It is distinguished between several cases 
    !   depending on the velocity direction and volume fraction 
    !   field.
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

        ! Local variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: kl, jl, il
        REAL(realk) :: dds, norms
        REAL(realk) :: dimx, dimy, dimz
        REAL(realk) :: flxedProp, flxWidth, flxAlpha

        CALL sel_ind(l, il, jl, kl)

        DO i = 2, ii-2
            DO j = 2, jj-2
                DO k = 2, kk-2
                    IF ( vel(k,j,i) > vofTol ) THEN

                        ! Compute face flx width and characteristic length
                        flxWidth = abs( vel(k,j,i) )*dt
                        dds = il*ddx(i) + jl*ddy(j) + kl*ddz(k)
                        IF ( flxWidth > 0.5_realk*dds  ) THEN
                            CALL err_abort(vofErr, "flxWidth > 0.5*cellWidth.", __FILE__, __LINE__)
                        ENDIF

                        IF ( isIfc(k,j,i) > 0.0_realk ) THEN
                            ! Compute norm
                            norms = il*normx(k,j,i) + jl*normy(k,j,i) + kl*normz(k,j,i)

                            ! Compute proper alpha
                            flxAlpha = alpha(k,j,i) - norms*( dds - flxWidth )

                            ! Compute dimensions of flxed cuboid
                            dimx = il*flxWidth + (1-il)*ddx(i)
                            dimy = jl*flxWidth + (1-jl)*ddy(j)
                            dimz = kl*flxWidth + (1-kl)*ddz(k)

                            ! Compute c in flxed cuboid
                            CALL comp_c_loc(flxedProp, flxAlpha, dimx, dimy, dimz, &
                                normx(k,j,i), normy(k,j,i), normz(k,j,i))
                        ELSE
                            flxedProp = c(k,j,i)
                        END IF
                    ELSE IF ( vel(k,j,i) < -vofTol ) THEN

                        ! Compute face flx width and characteristic length
                        flxWidth = abs( vel(k,j,i) )*dt
                        dds = il*ddx(i+il) + jl*ddy(j+jl) + kl*ddz(k+kl)
                        IF ( flxWidth > 0.5_realk*dds  ) THEN
                            CALL err_abort(vofErr, "flxWidth > 0.5*cellWidth.", __FILE__, __LINE__)
                        ENDIF

                        IF ( isIfc(k+kl,j+jl,i+il) > 0.0_realk ) THEN
                            ! Compute proper alpha
                            flxAlpha = alpha(k+kl,j+jl,i+il)

                            ! Compute dimensions of flxed cuboid
                            dimx = il*flxWidth + (1-il)*ddx(i+il)
                            dimy = jl*flxWidth + (1-jl)*ddy(j+jl)
                            dimz = kl*flxWidth + (1-kl)*ddz(k+kl)

                            ! Compute c in flxed cuboid
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
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: q, l
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        INTEGER(intk) :: n, igrid
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: c(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: advr(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: cFlx1(:,:,:), cFlx2(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)
        REAL(realk), POINTER, CONTIGUOUS :: normx(:,:,:), normy(:,:,:), normz(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: alpha(:,:,:)
        REAl(realk), POINTER, CONTIGUOUS :: isIfc(:,:,:)

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

            CALL comp_flx_stg_grd(kk, jj, ii, q, l, c, advr, &
                cFlx1, cFlx2, ddx, ddy, ddz, &
                normx, normy, normz, alpha, isIfc, dt)
        END DO

    END SUBROUTINE comp_flx_stg

    !================================================================

    SUBROUTINE comp_flx_stg_grd(kk, jj, ii, q, l, c, advr, &
        cFlx1, cFlx2, ddx, ddy, ddz, &
        normx, normy, normz, alpha, isIfc, dt)
    !----------------------------------------------------------------
    !   What it does:
    !   Computes the volume fraction flxes depending on the current
    !   split direction. It is distinguished between several cases 
    !   depending on the velocity direction and volume fraction 
    !   field.
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
            DO i = 2, ii-2
                DO j = 2, jj-2
                    DO k = 2, kk-2
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
                            CALL err_abort(vofErr, "flxWidth > 0.5*cellWidt.", __FILE__, __LINE__)
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
            DO i = 2, ii-2
                DO j = 2, jj-2
                    DO k = 2, kk-2

                    IF ( ABS(advr(k,j,i)) < vofTol ) THEN
                        cFlx1(k,j,i) = 0.0_realk
                        cFlx2(k,j,i) = 0.0_realk
                        CYCLE
                    ENDIF

                    ! Indizes of the half-pressure-cells overlapping the velocity-cells
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
    !   What it does:
    !   Defines the sequence, in which the volume fraction c and
    !   the momentum are advected. 
    !
    !   The PARIS solver uses a cyclic periodicity of three. Hence,
    !   the sequence can be
    !   x -> y -> z,
    !   y -> z -> x or
    !   z -> x -> y.
    !   
    !   Permutation can also be see as the sum of all possible 
    !   sequences. For three dimensions we get six sequences.
    !   x -> y -> z,
    !   y -> z -> x,
    !   z -> x -> y,
    !   x -> z -> y,
    !   y -> x -> z or
    !   z -> y -> x.
    !
    !   The splPer variable controls which version 
    !   is used.
    !
    !   Source:
    !   T. Arrufat et al., “A mass-momentum consistent, 
    !   Volume-of-Fluid method for incompressible flow on staggered 
    !   grids,” Computers & Fluids, vol. 215, p. 104785, Jan. 2021, 
    !   doi: 10.1016/j.compfluid.2020.104785.
    !
    !   W. Aniszewski et al., “PArallel, Robust, Interface Simulator
    !   (PARIS),” Computer Physics Communications, vol. 263, 
    !   p. 107849, Jun. 2021, doi: 10.1016/j.cpc.2021.107849.
    !   PARIS source code (accessed: Mai 2026)
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: iteration
        INTEGER(intk), INTENT(inout) :: advSeq(3)

        ! Local variables
        INTEGER(intk) :: perInd

        ! perInd only changes in a new time-step
        perInd = mod(iteration-1, splPer)

        ! Select permutation of split advection
        SELECT CASE ( perInd )
            CASE (0); advSeq = [1, 2, 3]
            CASE (1); advSeq = [3, 1, 2]
            CASE (2); advSeq = [2, 3, 1]
            CASE (3); advSeq = [1, 3, 2]
            CASE (4); advSeq = [3, 2, 1]
            CASE (5); advSeq = [2, 1, 3]
        END SELECT

    END SUBROUTINE def_adv_seq

    !================================================================

    SUBROUTINE adve_operator(dt, itstep)
    !----------------------------------------------------------------
    !   What it does:
    !   Performs the spatial and temporal integration of the
    !   advection operator. The integration is combined, since
    !   VOF/PLIC is "exact" up to the order of accuracy of the 
    !   interface reconstruction in PLIC.
    !   The structure is needed due to the split direction advection.
    !   Each directional split has to be performed on the hole
    !   domain (including restriction and prologation) before the
    !   next directional split can be performed. Otherwise, there are
    !   errors at the grids boundaries.
    !----------------------------------------------------------------

        ! Subroutine arguments
        REAL(realk), INTENT(in) :: dt
        INTEGER(intk), INTENT(in) :: itstep

        ! Local variables
        INTEGER(intk) :: q, advSeq(3), dirLoop, l

        IF ( skpAdv ) RETURN

        ! Initialize alpha, norm(.), cWy
        CALL comp_ifc()
        CALL comp_cWy()

        ! Initialize stg. grid cS(.), dS(.), mS(.), cWyS(.)
        DO q = 1, 3
            CALL comp_c_stg(q)
            CALL comp_d_stg(q)
            CALL comp_m_stg(q)
            CALL comp_cWy_stg(q)
        END DO

        ! Restriction
        CALL rstr(fldName="C", flag="D")
        DO q = 1, 3
            CALL rstr_stg(q=q, fldName="CS")
            CALL rstr_stg(q=q, fldName="DS")
            CALL rstr_stg(q=q, fldName="MS")
        END DO

        ! Prologation
        CALL prlg(fldName="C")
        CALL prlg_stg(fldName1="CS1", fldName2="CS2", fldName3="CS3")
        CALL prlg_stg(fldName1="DS1", fldName2="DS2", fldName3="DS3")
        CALL prlg_stg(fldName1="MS1", fldName2="MS2", fldName3="MS3")

        CALL def_adv_seq(itstep, advSeq)
        DO dirLoop = 1, 3
            l = advSeq(dirLoop)
                ! Advect stg. grid mS(.), cS(.), dS(.)
                DO q = 1, 3
                    CALL comp_isIfc_stg(q)
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

            ! Restriction
            CALL rstr(fldName="C", flag="D")
            DO q = 1, 3
                CALL rstr_stg(q=q, fldName="CS")
                CALL rstr_stg(q=q, fldName="DS")
                CALL rstr_stg(q=q, fldName="MS")
            END DO

            ! Prologation
            CALL prlg(fldName="C")
            CALL prlg_stg(fldName1="CS1", fldName2="CS2", fldName3="CS3")
            CALL prlg_stg(fldName1="DS1", fldName2="DS2", fldName3="DS3")
            CALL prlg_stg(fldName1="MS1", fldName2="MS2", fldName3="MS3")

            CALL comp_ifc()
        END DO

        ! Update u, v, w with mS(.) and dS(.)
        DO q = 1, 3
            CALL upd_vel_stg(q)
        END DO

    END SUBROUTINE adve_operator

    !================================================================

    SUBROUTINE diff_operator(uo_f, vo_f, wo_f)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: uo_f, vo_f, wo_f

        ! Local variables
        INTEGER(intk) :: n, igrid
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: up(:,:,:), vp(:,:,:), wp(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: uo(:,:,:), vo(:,:,:), wo(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: g(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: gUv(:,:,:), gUw(:,:,:), gVw(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: dBa(:,:,:), dLe(:,:,:), dTo(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: rdx(:), rdy(:), rdz(:)
        REAL(realk), POINTER, CONTIGUOUS :: rddx(:), rddy(:), rddz(:)

        IF ( skpDif ) RETURN

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_fieldptr(up, "UP", igrid)
            CALL get_fieldptr(vp, "VP", igrid)
            CALL get_fieldptr(wp, "WP", igrid)
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
            CALL get_fieldptr(rdx, "RDX", igrid)
            CALL get_fieldptr(rdy, "RDY", igrid)
            CALL get_fieldptr(rdz, "RDZ", igrid)
            CALL get_fieldptr(rddx, "RDDX", igrid)
            CALL get_fieldptr(rddy, "RDDY", igrid)
            CALL get_fieldptr(rddz, "RDDZ", igrid)

            CALL diff_operator_grd(kk, jj, ii, up, vp, wp, &
                g, gUv, gUw, gVw, dBa, dLe, dTo, &
                rdx, rdy, rdz, rddx, rddy, rddz, uo, vo, wo)
        END DO

    END SUBROUTINE diff_operator

    !================================================================

    SUBROUTINE diff_operator_grd(kk, jj, ii, up, vp, wp, g, &
        gUv, gUw, gVw, dBa, dLe, dTo, rdx, rdy, rdz, &
        rddx, rddy, rddz, uo, vo, wo)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------
    
        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: up(kk, jj, ii), vp(kk, jj, ii), wp(kk, jj, ii)
        REAL(realk), INTENT(in) :: g(kk, jj, ii)
        REAL(realk), INTENT(in) :: gUv(kk, jj, ii), gUw(kk, jj, ii), gVw(kk, jj, ii)
        REAL(realk), INTENT(in) :: dBa(kk, jj, ii), dLe(kk, jj, ii), dTo(kk, jj, ii)
        REAL(realk), INTENT(in) :: rdx(ii), rdy(jj), rdz(kk)
        REAL(realk), INTENT(in) :: rddx(ii), rddy(jj), rddz(kk)
        REAL(realk), INTENT(inout) :: uo(kk, jj, ii), vo(kk, jj, ii), wo(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: k, j, i
        REAL(realk) :: tauXxPl, tauXxMi, tauXyPl, tauXyMi, tauXzPl, tauXzMi
        REAL(realk) :: tauYxPl, tauYxMi, tauYyPl, tauYyMi, tauYzPl, tauYzMi
        REAL(realk) :: tauZxPl, tauZxMi, tauZyPl, tauZyMi, tauZzPl, tauZzMi

        DO i = 2, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    tauXxPl = g(k,j,i+1)*2.0_realk*(up(k,j,i+1) - up(k,j,i))*rddx(i+1)
                    tauXxMi = g(k,j,i)*2.0_realk*(up(k,j,i) - up(k,j,i-1))*rddx(i)
                    tauXyPl = gUv(k,j,i)*((up(k,j+1,i) - up(k,j,i))*rdy(j) + (vp(k,j,i+1) - vp(k,j,i))*rdx(i))
                    tauXyMi = gUv(k,j-1,i)*((up(k,j,i) - up(k,j-1,i))*rdy(j-1) + (vp(k,j-1,i+1) - vp(k,j-1,i))*rdx(i))
                    tauXzPl = gUw(k,j,i)*((up(k+1,j,i) - up(k,j,i))*rdz(k) + (wp(k,j,i+1) - wp(k,j,i))*rdx(i))
                    tauXzMi = gUw(k-1,j,i)*((up(k,j,i) - up(k-1,j,i))*rdz(k-1) + (wp(k-1,j,i+1) - wp(k-1,j,i))*rdx(i))

                    uo(k,j,i) = uo(k,j,i) + 1.0_realk/dBa(k,j,i)* &
                        ((tauXxPl - tauXxMi)*rdx(i) + &
                         (tauXyPl - tauXyMi)*rddy(j) + &
                         (tauXzPl - tauXzMi)*rddz(k))
                END DO
            END DO
        END DO

        DO i = 3, ii-2
            DO j = 2, jj-2
                DO k = 3, kk-2
                    tauYxPl = gUv(k,j,i)*((up(k,j+1,i) - up(k,j,i))*rdy(j) + (vp(k,j,i+1) - vp(k,j,i))*rdx(i))
                    tauYxMi = gUv(k,j,i-1)*((up(k,j+1,i-1) - up(k,j,i-1))*rdy(j) + (vp(k,j,i) - vp(k,j,i-1))*rdx(i-1))
                    tauYyPl = g(k,j+1,i)*2.0_realk*(vp(k,j+1,i) - vp(k,j,i))*rddy(j+1)
                    tauYyMi = g(k,j,i)*2.0_realk*(vp(k,j,i) - vp(k,j-1,i))*rddy(j)
                    tauYzPl = gVw(k,j,i)*((vp(k+1,j,i) - vp(k,j,i))*rdz(k) + (wp(k,j+1,i) - wp(k,j,i))*rdy(j))
                    tauYzMi = gVw(k-1,j,i)*((vp(k,j,i) - vp(k-1,j,i))*rdz(k-1) + (wp(k-1,j+1,i) - wp(k-1,j,i))*rdy(j))

                    vo(k,j,i) = vo(k,j,i) + 1.0_realk/dLe(k,j,i)* &
                        ((tauYxPl - tauYxMi)*rddx(i) + &
                         (tauYyPl - tauYyMi)*rdy(j) + &
                         (tauYzPl - tauYzMi)*rddz(k))
                END DO
            END DO
        END DO

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 2, kk-2
                    tauZxPl = gUw(k,j,i)*((up(k+1,j,i) - up(k,j,i))*rdz(k) + (wp(k,j,i+1) - wp(k,j,i))*rdx(i))
                    tauZxMi = gUw(k,j,i-1)*((up(k+1,j,i-1) - up(k,j,i-1))*rdz(k) + (wp(k,j,i) - wp(k,j,i-1))*rdx(i-1))
                    tauZyPl = gVw(k,j,i)*((vp(k+1,j,i) - vp(k,j,i))*rdz(k) + (wp(k,j+1,i) - wp(k,j,i))*rdy(j))
                    tauZyMi = gVw(k,j-1,i)*((vp(k+1,j-1,i) - vp(k,j-1,i))*rdz(k) + (wp(k,j,i) - wp(k,j-1,i))*rdy(j-1))
                    tauZzPl = g(k+1,j,i)*2.0_realk*(wp(k+1,j,i) - wp(k,j,i))*rddz(k+1)
                    tauZzMi = g(k,j,i)*2.0_realk*(wp(k,j,i) - wp(k-1,j,i))*rddz(k)

                    wo(k,j,i) = wo(k,j,i) + 1.0_realk/dTo(k,j,i)*&
                        ((tauZxPl - tauZxMi)*rddx(i) + &
                         (tauZyPl - tauZyMi)*rddy(j) + &
                         (tauZzPl - tauZzMi)*rdz(k))
                END DO
            END DO
        END DO

    END SUBROUTINE diff_operator_grd

    !================================================================

    SUBROUTINE pres_operator(uo_f, vo_f, wo_f)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: uo_f, vo_f, wo_f

        ! Local variables
        INTEGER(intk) :: n, igrid
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: c(:,:,:), p(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: uo(:,:,:), vo(:,:,:), wo(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: dBa(:,:,:), dLe(:,:,:), dTo(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: rdx(:), rdy(:), rdz(:)
        INTEGER(intk) :: gradpflag

        IF ( skpPre ) RETURN

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_fieldptr(c, "C", igrid)
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

            CALL get_gradpxflag(gradpflag, igrid)
            CALL pres_operator_grd(kk, jj, ii, c, p, dBa, dLe, dTo, &
                rdx, rdy, rdz, uo, vo, wo, gradpflag)
        END DO

    END SUBROUTINE pres_operator

    !================================================================

    SUBROUTINE pres_operator_grd(kk, jj, ii, c, p, dBa, dLe, dTo, &
        rdx, rdy, rdz, uo, vo, wo, gradpflag)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(in) :: c(kk, jj, ii), p(kk, jj, ii)
        REAL(realk), INTENT(in) :: dBa(kk, jj, ii), dLe(kk, jj, ii), dTo(kk, jj, ii)
        REAL(realk), INTENT(in) :: rdx(ii), rdy(jj), rdz(kk)
        REAL(realk), INTENT(inout) :: uo(kk, jj, ii), vo(kk, jj, ii), wo(kk, jj, ii)
        INTEGER(intk), INTENT(in) :: gradpflag

        ! Local variables
        INTEGER(intk) :: i, j, k
        REAL(realk) :: gpx(kk, jj, ii), gpy(kk, jj, ii), gpz(kk, jj, ii)

        gpx = 0.0_realk
        gpy = 0.0_realk
        gpz = 0.0_realk
        DO i = 2, ii-2
            DO j = 2, jj-2
                DO k = 2, kk-2
                    gpx(k,j,i) = gradp(1)*gradpflag*MERGE(1.0_realk, 0.0_realk, c(k,j,i) > vofTol)
                    gpy(k,j,i) = gradp(2)*gradpflag*MERGE(1.0_realk, 0.0_realk, c(k,j,i) > vofTol)
                    gpz(k,j,i) = gradp(3)*gradpflag*MERGE(1.0_realk, 0.0_realk, c(k,j,i) > vofTol)
                ENDDO
            ENDDO
        ENDDO

        DO i = 2, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    uo(k,j,i) = uo(k,j,i) - 1.0_realk/dBa(k,j,i)*((p(k,j,i+1) - p(k,j,i))*rdx(i) + gpx(k,j,i))
                END DO
            END DO
        END DO

        DO i = 3, ii-2
            DO j = 2, jj-2
                DO k = 3, kk-2
                    vo(k,j,i) = vo(k,j,i) - 1.0_realk/dLe(k,j,i)*((p(k,j+1,i) - p(k,j,i))*rdy(j) + gpy(k,j,i))
                END DO
            END DO
        END DO

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 2, kk-2
                    wo(k,j,i) = wo(k,j,i) - 1.0_realk/dTo(k,j,i)*((p(k+1,j,i) - p(k,j,i))*rdz(k) + gpz(k,j,i))
                END DO
            END DO
        END DO

    END SUBROUTINE pres_operator_grd

    !================================================================

    SUBROUTINE exte_operator(uo_f, vo_f, wo_f)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        TYPE(field_t), INTENT(inout) :: uo_f, vo_f, wo_f

        ! Local variables
        INTEGER(intk) :: n, igrid
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: uo(:,:,:), vo(:,:,:), wo(:,:,:)
        

        IF ( skpExt ) RETURN

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL uo_f%get_ptr(uo, igrid)
            CALL vo_f%get_ptr(vo, igrid)
            CALL wo_f%get_ptr(wo, igrid)

            CALL exte_operator_grd(kk, jj, ii, uo, vo, wo)
        END DO

    END SUBROUTINE exte_operator

    !================================================================

    SUBROUTINE exte_operator_grd(kk, jj, ii, uo, vo, wo)
    !----------------------------------------------------------------
    !   What it does:
    !    
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii
        REAL(realk), INTENT(inout) :: uo(kk, jj, ii), vo(kk, jj, ii), wo(kk, jj, ii)

        ! Local variables
        INTEGER(intk) :: i, j, k

        DO i = 2, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    uo(k,j,i) = uo(k,j,i) + grav(1)
                ENDDO
            ENDDO
        ENDDO

        DO i = 3, ii-2
            DO j = 2, jj-2
                DO k = 3, kk-2
                    vo(k,j,i) = vo(k,j,i) + grav(2)
                ENDDO
            ENDDO
        ENDDO

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 2, kk-2
                    wo(k,j,i) = wo(k,j,i) + grav(3)
                ENDDO
            ENDDO
        ENDDO

    END SUBROUTINE exte_operator_grd

    !================================================================

    SUBROUTINE adv_c(l, dt)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: l
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        INTEGER(intk) :: n, igrid
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: c(:,:,:), cWy(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: u(:,:,:), v(:,:,:), w(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: vel(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: cFlx1(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_fieldptr(c, "C", igrid)
            CALL get_fieldptr(cWy, "CWY", igrid)
            CALL get_fieldptr(u, "U", igrid)
            CALL get_fieldptr(v, "V", igrid)
            CALL get_fieldptr(w, "W", igrid)
            CALL get_fieldptr(cFlx1, "CFLX1", igrid)
            CALL get_fieldptr(dx, "DX", igrid)
            CALL get_fieldptr(dy, "DY", igrid)
            CALL get_fieldptr(dz, "DZ", igrid)
            CALL get_fieldptr(ddx, "DDX", igrid)
            CALL get_fieldptr(ddy, "DDY", igrid)
            CALL get_fieldptr(ddz, "DDZ", igrid)

            CALL sel_vel(l, u, v, w, vel)
            CALL adv_c_grd(kk, jj, ii, 0, l, c, cWy, vel, cFlx1, &
                dx, dy, dz, ddx, ddy, ddz, dt)
        END DO

    END SUBROUTINE adv_c

    !================================================================

    SUBROUTINE adv_c_stg(q, l, dt)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: q, l
        REAL(realk), INTENT(in) :: dt

        ! Local variables
        CHARACTER(len=5) :: cFldName, cWyFldName
        INTEGER(intk) :: n, igrid
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: cSq(:,:,:), cWySq(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: advr(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: cFlx1(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)

        cFldName = "CS"//int2Char(q)
        cWyFldName = "CWYS"//int2Char(q)

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_fieldptr(cSq, cFldName, igrid)
            CALL get_fieldptr(cWySq, cWyFldName, igrid)
            CALL get_fieldptr(advr, "ADVR", igrid)
            CALL get_fieldptr(cFlx1, "CFLX1", igrid)
            CALL get_fieldptr(dx, "DX", igrid)
            CALL get_fieldptr(dy, "DY", igrid)
            CALL get_fieldptr(dz, "DZ", igrid)
            CALL get_fieldptr(ddx, "DDX", igrid)
            CALL get_fieldptr(ddy, "DDY", igrid)
            CALL get_fieldptr(ddz, "DDZ", igrid)

            CALL adv_c_grd(kk, jj, ii, q, l, cSq, cWySq, advr, cFlx1, &
                dx, dy, dz, ddx, ddy, ddz, dt)
        END DO

    END SUBROUTINE adv_c_stg

    !================================================================

    SUBROUTINE adv_c_grd(kk, jj, ii, q, l, c, cWy, vel, cFlx1, &
        dx, dy, dz, ddx, ddy, ddz, dt)
    !----------------------------------------------------------------
    !   What it does:
    !    
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

        ! Local variables
        INTEGER(intk) :: i, j, k
        INTEGER(intk) :: il, jl, kl
        REAL(realk) :: dsx(ii), dsy(jj), dsz(kk), dsCV
        REAL(realk) :: div

        CALL sel_ind(l, il, jl, kl)
        CALL sel_ext(kk, jj, ii, q, l, dx, dy, dz, ddx, ddy, ddz, dsx, dsy, dsz)

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
                    dsCV = il*dsx(i) + jl*dsy(j) + kl*dsz(k)
                    div = ( vel(k,j,i) - vel(k-kl,j-jl,i-il) ) / dsCV
                    c(k,j,i) = c(k,j,i) &
                        - dt/dsCV*( cFlx1(k,j,i) - cFlx1(k-kl,j-jl,i-il) ) &
                        + dt*cWy(k,j,i)*div
                END DO
            END DO
        END DO

    END SUBROUTINE adv_c_grd

    !================================================================

    SUBROUTINE adv_m_stg(q, l, dt)
    !----------------------------------------------------------------
    !   What it does:
    !   
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
        REAL(realk), POINTER, CONTIGUOUS :: vel(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: advr(:,:,:), adve(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: mSq(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: cFlx1(:,:,:), cFlx2(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: dx(:), dy(:), dz(:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)

        cWyFldName = "CWYS"//int2Char(q)
        mFldName = "MS"//int2Char(q)

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_fieldptr(cWySq, cWyFldName, igrid)
            CALL get_fieldptr(u, "U", igrid)
            CALL get_fieldptr(v, "V", igrid)
            CALL get_fieldptr(w, "W", igrid)
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

            CALL sel_vel(q, u, v, w, vel)
            CALL adv_m_grd(kk, jj, ii, q, l, cWySq, vel, &
                advr, adve, mSq, cFlx1, cFlx2, dx, dy, dz, &
                ddx, ddy, ddz, dt)
        END DO

    END SUBROUTINE adv_m_stg

    !================================================================

    SUBROUTINE adv_m_grd(kk, jj, ii, q, l, cWy, vel, &
        advr, adve, mom, cFlx1, cFlx2, dx, dy, dz, &
        ddx, ddy, ddz, dt)
    !----------------------------------------------------------------
    !   What it does:
    !    
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

        ! Local variables
        INTEGER(intk) :: i, j, k
        INTEGER(intk) :: il, jl, kl
        REAL(realk) :: dsx(ii), dsy(jj), dsz(kk), dsCV
        REAL(realk) :: momFlx(kk, jj, ii)
        REAL(realk) :: div, com

        CALL sel_ind(l, il, jl, kl)
        CALL sel_ext(kk, jj, ii, q, l, dx, dy, dz, ddx, ddy, ddz, dsx, dsy, dsz)

        DO i = 2, ii-2
            DO j = 2, jj-2
                DO k = 2, kk-2
                    momFlx(k,j,i) = adve(k,j,i)*( rho1*cFlx1(k,j,i) + rho2*cFlx2(k,j,i) )
                END DO
            END DO 
        END DO

        DO i = 3, ii-2
            DO j = 3, jj-2
                DO k = 3, kk-2
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
    !   What it does:
    !   Compute the volume fraction field for the staggered cells
    !   depending on q. The staggered cells are either moved by
    !   1/2 ddx, 1/2 ddy or 1/2 ddz.
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
    !   What it does:
    !    
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
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: q

        ! Local variables
        CHARACTER(len=3) :: dFldName, mFldName
        INTEGER(intk) :: n, igrid
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: u(:,:,:), v(:,:,:), w(:,:,:)
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

            CALL upd_vel_stg_grd(kk, jj, ii, q, dSq, u, v, w, mSq)
        END DO

    END SUBROUTINE upd_vel_stg

    !================================================================

    SUBROUTINE upd_vel_stg_grd(kk, jj, ii, q, dSq, u, v, w, mSq)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, q
        REAL(realk), INTENT(in) :: dSq(kk, jj, ii)
        REAL(realk), INTENT(inout) :: u(kk, jj, ii), v(kk, jj, ii), w(kk, jj, ii)
        REAL(realk), INTENT(in) :: mSq(kk, jj, ii)

        ! Local variables
        ! None

        IF (q == 1) THEN
            u = mSq/dSq
        ELSEIF (q == 2) THEN
            v = mSq/dSq
        ELSEIF (q == 3) THEN
            w = mSq/dSq
        ENDIF

    END SUBROUTINE upd_vel_stg_grd

    !================================================================

    SUBROUTINE comp_cWy()
    !----------------------------------------------------------------
    !   What it does:
    !   Comput nondirectional compression coefficient cWy on 
    !   multi-grid level.
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
    !   What it does:
    !   Comput nondirectional compression coefficient cWy on 
    !   multi-grid level.
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
    !   What it does:
    !   Computes the nondirectional compression coefficient c for 
    !   Weymouth and Yue"s advection scheme.
    !
    !   Source:
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

        CALL get_field(fld_p, fldName)
        CALL rstr_grd(fld_p, flag)

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

        name = fldName//int2Char(q)

        SELECT CASE ( q )
        CASE ( 1 ); flag = "A"
        CASE ( 2 ); flag = "B"
        CASE ( 3 ); flag = "C"
        END SELECT

        CALL get_field(fld_p, name)
        CALL rstr_grd(fld_p, flag)

    END SUBROUTINE rstr_stg

    !================================================================

    SUBROUTINE rstr_grd(fld_p, flag)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        TYPE(field_t), POINTER :: fld_p
        CHARACTER(len=*), INTENT(in) :: flag

        ! Local variables
        INTEGER(intk) :: ilevel

        DO ilevel = maxlevel, minlevel, -1
            CALL ftoc(ilevel, fld_p%arr, fld_p%arr, flag)
        END DO

    END SUBROUTINE rstr_grd

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
            CALL connect(ilevel, 2, s1=fld_p, corners=.TRUE.)
        END DO

    END SUBROUTINE prlg

    !================================================================

    SUBROUTINE prlg_stg(fldName1, fldName2, fldName3)
    !----------------------------------------------------------------
    !   What it does:
    !   
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
            CALL connect(ilevel, 2, v1=fld1_p, v2=fld2_p, v3=fld3_p, corners=.TRUE.)
        END DO

    END SUBROUTINE prlg_stg

    !================================================================

    SUBROUTINE comp_adve_stg(q, l, dt)
    !----------------------------------------------------------------
    !   What it does:
    !   
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
            CALL comp_adve_stg_grd(kk, jj, ii, q, l, vel, advr, adve, &
                isIfcVicSq, dx, dy, dz, ddx, ddy, ddz, dt)
        END DO

    END SUBROUTINE comp_adve_stg

    !================================================================

    SUBROUTINE comp_adve_stg_grd(kk, jj, ii, q, l, vel, advr, adve, &
        isIfcVic, dx, dy, dz, ddx, ddy, ddz, dt)
    !----------------------------------------------------------------
    !   What it does:
    !   Selects interpolation scheme.
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

        ! Local variables
        ! None

        IF ( advScm == "QUICK" ) THEN
            CALL comp_adve_stg_quick(kk, jj, ii, q, l, vel, advr, &
                adve, isIfcVic, dx, dy, dz, ddx, ddy, ddz, dt)
        ELSEIF ( advScm == "ENO" ) THEN
            CALL comp_adve_stg_eno(kk, jj, ii, q, l, vel, advr, &
                adve, dx, dy, dz, ddx, ddy, ddz, dt)
        ELSE
            CALL err_abort(vofErr, "Unknown interpolation.", __FILE__, __LINE__)
        ENDIF

    END SUBROUTINE comp_adve_stg_grd

    !================================================================

    SUBROUTINE comp_adve_stg_quick(kk, jj, ii, q, l, vel, &
        advr, adve, isIfcVic, dx, dy, dz, ddx, ddy, ddz, dt)
    !----------------------------------------------------------------
    !   What it does:
    !   QUICK interpolation to compute the advected veloctiy
    !   (advectee) on staggered grid cells.
    !   adve = advected q (advectee)
    !   advr = advecting q (advector)
    !   An indicator function is used to avoid if-statements within
    !   loops.
    !
    !   Sources: 
    !   T. Arrufat et al., “A mass-momentum consistent, 
    !   Volume-of-Fluid method for incompressible flow on staggered 
    !   grids,” Computers & Fluids, vol. 215, p. 104785, Jan. 2021, 
    !   doi: 10.1016/j.compfluid.2020.104785.
    !
    !   B. P. Leonard, “A stable and accurate convective modelling 
    !   procedure based on quadratic upstream interpolation,”
    !   Computer Methods in Applied Mechanics and Engineering,
    !   vol. 19, no. 1, pp. 59–98, Jun. 1979,
    !   doi: 10.1016/0045-7825(79)90034-3.
    !
    !   W. Aniszewski et al., “PArallel, Robust, Interface Simulator
    !   (PARIS),” Computer Physics Communications, vol. 263,
    !   p. 107849, Jun. 2021, doi: 10.1016/j.cpc.2021.107849.
    !   PARIS source code, grep "interpole_quad" (accessed: Mai 2026)
    !
    !   D. Herrmann, Numerische Mathematik — 40 BASIC-Programme. 
    !   Wiesbaden: Vieweg+Teubner Verlag, 1983. 
    !   doi: 10.1007/978-3-322-96321-5.
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

        DO i = 2, ii-2
            DO j = 2, jj-2
                DO k = 2, kk-2
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
        ! Both, Lagrange and Newton form yield the same results.
        ! The Lagrange form was more familiar for me. Hence, I 
        ! implemented it first. The Newton form is closer to 
        ! Leonard"s formulation and easier to adjust in the future
        ! (see QUICKEST scheme).
        !------------------------------------------------------------
        ! PURE REAL(realk) FUNCTION lagrange_interpolation(phiUU, phiU, phiD, dsUU, dsD, dsr) RESULT(r)
        !     REAL(realk), INTENT(in) :: phiUU, phiU, phiD, dsUU, dsD, dsr
        !     r = phiUU*(dsr*(dsr-dsD))/(dsUU*(dsUU+dsD)) &
        !         - phiU*((dsr+dsUU)*(dsr-dsD))/(dsUU*dsD) &
        !         + phiD*(dsr*(dsr+dsUU))/(dsD*(dsUU+dsD))
        ! END FUNCTION lagrange_interpolation

        PURE REAL(realk) FUNCTION newton_interpolation(phiUU, phiU, phiD, dsUU, dsD, dsr) RESULT(r)
            REAL(realk), INTENT(in) :: phiUU, phiU, phiD, dsUU, dsD, dsr
            REAL(realk) :: dd1U, dd1D, curv
            dd1D = (phiD - phiU)/dsD
            dd1U = (phiU - phiUU)/dsUU
            curv = (dd1D - dd1U)/(dsUU + dsD)
            r = phiU + dsr*dd1D + dsr*(dsr - dsD)*curv
        END FUNCTION newton_interpolation

    END SUBROUTINE comp_adve_stg_quick

    !================================================================

    SUBROUTINE comp_adve_stg_eno(kk, jj, ii, q, l, vel, &
        advr, adve, dx, dy, dz, ddx, ddy, ddz, dt)
    !----------------------------------------------------------------
    !   What it does:
    !   ENO interpolation to compute the advected veloctiy
    !   (advectee) on staggered grid cells.
    !   adve = advected q (advectee)
    !   advr = advecting q (advector)
    !   An indicator function is used to avoid if-statements within
    !   loops.
    !   
    !   Sources: 
    !   G. Tryggvason, R. Scardovelli, and S. Zaleski, Direct
    !   Numerical Simulations of Gas–Liquid mph Flows,
    !   1st ed. Cambridge University Press, 2011.
    !   doi: 10.1017/CBO9780511975264.
    !
    !   P. K. Sweby, “High Resolution Schemes Using Flx Limiters
    !   for Hyperbolic Conservation Laws,” SIAM J. Numer. Anal.,
    !   vol. 21, no. 5, pp. 995–1011, Oct. 1984,
    !   doi: 10.1137/0721062.
    !
    !   W. Aniszewski et al., “PArallel, Robust, Interface Simulator
    !   (PARIS),” Computer Physics Communications, vol. 263,
    !   p. 107849, Jun. 2021, doi: 10.1016/j.cpc.2021.107849.
    !   PARIS source code, grep "interpole3" (accessed: Mai 2026)
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: kk, jj, ii, q, l
        REAL(realk), INTENT(in) :: vel(kk, jj, ii)
        REAL(realk), INTENT(in) :: advr(kk, jj, ii)
        REAL(realk), INTENT(inout) :: adve(kk, jj, ii)
        REAL(realk), INTENT(in) :: dx(ii), dy(jj), dz(kk)
        REAL(realk), INTENT(in) :: ddx(ii), ddy(jj), ddz(kk)
        REAL(realk), INTENT(in) :: dt

        ! Loval variables
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

        DO i = 2, ii-2
            DO j = 2, jj-2
                DO k = 2, kk-2
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
            REAL(realk), INTENT(in) :: a, b
            res = 0.5_realk*( SIGN(1.0_realk, a) + SIGN(1.0_realk, b) )*MIN(ABS(a), ABS(b))
        END FUNCTION minmod

    END SUBROUTINE comp_adve_stg_eno

    !================================================================

    SUBROUTINE comp_advr_stg(q, l)
    !----------------------------------------------------------------
    !   What it does:
    !   
    !----------------------------------------------------------------

        ! Subroutine arguments
        INTEGER(intk), INTENT(in) :: q, l

        ! Local variables
        INTEGER(intk) :: n, igrid
        INTEGER(intk) :: kk, jj, ii
        REAL(realk), POINTER, CONTIGUOUS :: u(:,:,:), v(:,:,:), w(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: vel(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: advr(:,:,:)
        REAL(realk), POINTER, CONTIGUOUS :: ddx(:), ddy(:), ddz(:)

        DO n = 1, nmygrids
            igrid = mygrids(n)

            CALL get_mgdims(kk, jj, ii, igrid)

            CALL get_fieldptr(u, "U", igrid)
            CALL get_fieldptr(v, "V", igrid)
            CALL get_fieldptr(w, "W", igrid)
            CALL get_fieldptr(advr, "ADVR", igrid)
            CALL get_fieldptr(ddx, "DDX", igrid)
            CALL get_fieldptr(ddy, "DDY", igrid)
            CALL get_fieldptr(ddz, "DDZ", igrid)

            CALL sel_vel(l, u, v, w, vel)
            CALL comp_advr_stg_grd(kk, jj, ii, q, l, vel, advr, &
                ddx, ddy, ddz)
        END DO

    END SUBROUTINE comp_advr_stg

    !================================================================

    SUBROUTINE comp_advr_stg_grd(kk, jj, ii, q, l, &
        vel, advr, ddx, ddy, ddz)
    !----------------------------------------------------------------
    !   What it does:
    !   Linear interpolation to compute the advecting velocity
    !   (advector) on staggered grid cells.
    !   advr = advecting q (advector)
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

        ! Loval variables
        INTEGER(intk) :: k, j, i
        INTEGER(intk) :: kq, jq, iq
        REAL(realk) :: ddnqMi, ddnqPl

        CALL sel_ind(q, iq, jq, kq)

        IF ( q == l ) THEN
            DO i = 2, ii-2
                DO j = 2, jj-2
                    DO k = 2, kk-2
                        advr(k,j,i) = 0.5_realk*(vel(k,j,i) + vel(k+kq,j+jq,i+iq))
                    ENDDO
                ENDDO
            ENDDO
        ELSE
            DO i = 2, ii-2
                DO j = 2, jj-2
                    DO k = 2, kk-2
                        ddnqMi = iq*ddx(i) + jq*ddy(j) + kq*ddz(k)
                        ddnqPl = iq*ddx(i+iq) + jq*ddy(j+jq) + kq*ddz(k+kq)
                        advr(k,j,i) = (vel(k,j,i)*ddnqPl + vel(k+kq,j+jq,i+iq)*ddnqMi)/(ddnqMi + ddnqPl)
                    ENDDO
                ENDDO
            ENDDO
        ENDIF

    END SUBROUTINE comp_advr_stg_grd

END MODULE mph_vof_mod
