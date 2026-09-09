! --------------------------------------------------------------------
! Copyright (C) 1991 - 2026 - EDF - www.code-aster.org
! This file is part of code_aster.
!
! code_aster is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! code_aster is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
! --------------------------------------------------------------------
!
subroutine dkqrig(plateCara, plateOrie, &
                  xyzl, option, pgl, &
                  matrRigi_, ener_)
!
    use plate_type
    implicit none
!
#include "asterc/r8gaem.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/bsthpl.h"
#include "asterfort/dkqbf.h"
#include "asterfort/dkqshp.h"
#include "asterfort/dxmate.h"
#include "asterfort/dxqbm.h"
#include "asterfort/dxqgm.h"
#include "asterfort/dxqloc.h"
#include "asterfort/dxqlocdri1.h"
#include "asterfort/dxqlocdri2.h"
#include "asterfort/dxqlocdri3.h"
#include "asterfort/dxqlocdri4.h"
#include "asterfort/dxqloe_NV.h"
#include "asterfort/dxqloe.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/gquad4.h"
#include "asterfort/jevech.h"
#include "asterfort/jquad4.h"
#include "asterfort/r8inir.h"
#include "asterfort/utbtab.h"
#include "asterfort/utctab.h"
#include "asterfort/utpvgl.h"
#include "blas/dcopy.h"
#include "blas/dscal.h"
#include "jeveux.h"
!
    type(plateOrie_Para), intent(in) :: plateOrie
    type(plateCara_Para), intent(in) :: plateCara
    real(kind=8), intent(in) :: xyzl(3, *), pgl(*)
    character(len=16), intent(in) :: option
    real(kind=8), optional, intent(out) :: matrRigi_(300), ener_(3)
!
! --------------------------------------------------------------------------------------------------
!
!     MATRICE DE RIGIDITE DE L'ELEMENT DE PLAQUE DKQ
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: npg, ndim, ipoids, icoopg
    integer(kind=8) :: multic, i, jvDisp
    real(kind=8), parameter :: un = 1.d0
    real(kind=8) :: wgt
    real(kind=8) :: df(9), dm(9), dmf(9), dc(4), dci(4)
    real(kind=8) :: df2(9), dm2(9), dmf2(9)
    real(kind=8) :: dmc(3, 2), dfc(3, 2)
    real(kind=8) :: bf(3, 12), bm(3, 8)
    real(kind=8) :: xab1(3, 12), depl(24), caraq4(25), jacob(5)
    real(kind=8) :: qsi, eta
    real(kind=8) :: flex(144), memb(64), mefl(96)
    real(kind=8) :: bsigth(24), excent, ctor
    aster_logical :: coupmf, exce
    integer(kind=8) :: j, ii, jj, irot
    integer(kind=8), parameter :: npgmx = 9
    real(kind=8) :: shp(3, 4, npgmx), shpr1(3, 4, npgmx), shpr2(3, 4, npgmx), bb(12, npgmx)
    real(kind=8) :: gshp1(3, 4), gshp2(3, 4)
    real(kind=8) :: dArea, gam, epais, fact, gm(3, 4)
    real(kind=8) :: gmemb(4, 4), btgmemb(8, 4), gmefl(4, 12)
    real(kind=8) :: bxb(12, 12)
    aster_logical :: dri
    blas_int :: b_incx, b_incy, b_n
    real(kind=8) :: matrRigi(300), ener(3)
!
! --------------------------------------------------------------------------------------------------
!
    df = 0.d0
    dm = 0.d0
    dmf = 0.d0
    dc = 0.d0
    dci = 0.0
    df2 = 0.d0
    dm2 = 0.d0
    dmf2 = 0.0
    dmc = 0.d0
    dfc = 0.0
    bf = 0.d0
    bm = 0.0
    xab1 = 0.d0
    depl = 0.d0
    caraq4 = 0.d0
    jacob = 0.0
    qsi = 0.d0
    eta = 0.0
    bsigth = 0.d0
    excent = 0.d0
    ctor = 0.0
    dArea = 0.d0
    gam = 0.d0
    epais = 0.d0
    fact = 0.0
    coupmf = ASTER_FALSE
    ener = 0.d0
    matrRigi = 0.d0
!
    call elrefe_info(fami='RIGI', npg=npg, ndim=ndim, &
                     jpoids=ipoids, jcoopg=icoopg)

! - Get parameters
    ctor = plateCara%coefRigiDRZ
    excent = plateCara%offset
    epais = plateCara%thick

! - Flags
    dri = ctor .lt. 0.0d0
    exce = (abs(excent) .gt. un/r8gaem())

! - Geometric properties
    call gquad4(xyzl, caraq4)

! - Get elementary matrix of rigidity
    call dxmate(plateCara, plateOrie, &
                'RIGI', df, dm, dmf, dc, &
                dci, dmc, dfc, &
                multic, coupmf)

    memb = 0.d0
    flex = 0.d0
    mefl = 0.d0
!
    if (dri) then
        call r8inir(12*npgmx, 0.d0, shp, 1)
        call r8inir(12*npgmx, 0.d0, shpr1, 1)
        call r8inir(12*npgmx, 0.d0, shpr2, 1)
        call r8inir(12*npgmx, 0.d0, bb, 1)
!
        dArea = 0.0d0
!
        call r8inir(12, 0.d0, gshp1, 1)
        call r8inir(12, 0.d0, gshp2, 1)
        call r8inir(12, 0.d0, gm, 1)
        call r8inir(16, 0.d0, gmemb, 1)
        call r8inir(32, 0.d0, btgmemb, 1)
        call r8inir(48, 0.d0, gmefl, 1)
        call r8inir(144, 0.d0, bxb, 1)
!
        gam = abs(ctor)*dm(1)
        do ii = 1, npg
!
!        ----- COORDINATES :
            qsi = zr(icoopg-1+ndim*(ii-1)+1)
            eta = zr(icoopg-1+ndim*(ii-1)+2)
!
!        ----- JACOBIAN AND WEIGHT :
            call jquad4(xyzl, qsi, eta, jacob)
            wgt = zr(ipoids+ii-1)*jacob(1)

!        -- ELEMENT AREA :
            dArea = dArea+wgt
!
!        -- COMPUTE LINEAR AND ROTATIONAL SHAPE FUNCTIONS AND DERIVATIVES :
            call dkqshp(qsi, eta, caraq4, jacob, shp(1, 1, ii), &
                        shpr1(1, 1, ii), shpr2(1, 1, ii))
!
            do j = 1, 4
                do i = 1, 3
                    gshp1(i, j) = gshp1(i, j)+shpr1(i, j, ii)*wgt
                    gshp2(i, j) = gshp2(i, j)+shpr2(i, j, ii)*wgt
                end do
            end do
        end do
!
        do ii = 1, npg
            do j = 1, 4
                do i = 1, 3
                    shpr1(i, j, ii) = shpr1(i, j, ii)-gshp1(i, j)/dArea
                    shpr2(i, j, ii) = shpr2(i, j, ii)-gshp2(i, j)/dArea
                end do
            end do
            do i = 1, 4
                j = 3*(i-1)
                bb(1+j, ii) = bb(1+j, ii)-shp(2, i, ii)
                bb(2+j, ii) = bb(2+j, ii)+shp(1, i, ii)
                bb(3+j, ii) = bb(3+j, ii)-2.d0*shp(3, i, ii)
                bb(3+j, ii) = bb(3+j, ii)-shpr1(2, i, ii)+shpr2(1, i, ii)
            end do
        end do
    end if

    do i = 1, npg
        qsi = zr(icoopg-1+ndim*(i-1)+1)
        eta = zr(icoopg-1+ndim*(i-1)+2)
!        ----- CALCUL DU JACOBIEN SUR LE QUADRANGLE --------------------
        call jquad4(xyzl, qsi, eta, jacob)
        wgt = zr(ipoids+i-1)*jacob(1)
!
!        -- FLEXION :
        call dkqbf(qsi, eta, jacob(2), caraq4, bf)
!        ----- CALCUL DU PRODUIT BFT.DF.BF -----------------------------
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, df, b_incx, df2, b_incy)
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        call dscal(b_n, wgt, df2, b_incx)
        call utbtab('CUMU', 3, 12, df2, bf, &
                    xab1, flex)
!
!        -- MEMBRANE :
        call dxqbm(qsi, eta, jacob(2), bm)
!        ----- CALCUL DU PRODUIT BMT.DM.BM -----------------------------
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, dm, b_incx, dm2, b_incy)
        b_n = to_blas_int(9)
        b_incx = to_blas_int(1)
        call dscal(b_n, wgt, dm2, b_incx)
        call utbtab('CUMU', 3, 8, dm2, bm, &
                    xab1, memb)
!
!   compute rotational part of membrane B matrix
!
        if (dri) then
!=====================================================================
! ---  CALCUL DE LA PARTIE MEMBRANE DE LA MATRICE DE MASSE =
! ---  LES TERMES SONT EN NK*NP                                      =
!=====================================================================
!
!        -- MEMBRANE (DRILLING PART) Gm:
            call dxqgm(shpr1(1, 1, i), shpr2(1, 1, i), gm)
!
!        ----- CALCUL DU PRODUIT GMT.DM.GM
            b_n = to_blas_int(9)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, dm, b_incx, dm2, b_incy)
            b_n = to_blas_int(9)
            b_incx = to_blas_int(1)
            call dscal(b_n, wgt, dm2, b_incx)
            call utbtab('CUMU', 3, 4, dm2, gm, &
                        xab1, gmemb)
!        ----- CALCUL DU PRODUIT BMT.DM.GM
            b_n = to_blas_int(9)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, dm, b_incx, dm2, b_incy)
            b_n = to_blas_int(9)
            b_incx = to_blas_int(1)
            call dscal(b_n, wgt, dm2, b_incx)
            call utctab('CUMU', 3, 4, 8, dm2, &
                        gm, bm, xab1, btgmemb)
!
!        ----- CALCUL DU PRODUIT gam/Omega*b(x)b
!
            do irot = 1, 12
                fact = wgt*gam/dArea*bb(irot, i)
                do jj = 1, 12
                    bxb(irot, jj) = bxb(irot, jj)+fact*bb(jj, i)
                end do
            end do
        end if
!
!
!
!        -- COUPLAGE :
        if (coupmf .or. exce) then
            if (dri) then
!           ----- CALCUL DU PRODUIT BMT.DMF.BF -------------------------
                b_n = to_blas_int(9)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call dcopy(b_n, dmf, b_incx, dmf2, b_incy)
                b_n = to_blas_int(9)
                b_incx = to_blas_int(1)
                call dscal(b_n, wgt, dmf2, b_incx)
                call utctab('CUMU', 3, 12, 8, dmf2, &
                            bf, bm, xab1, mefl)
!
!   compute product Gmt.Dmf.Bf
!           ----- CALCUL DU PRODUIT GMT.DMF.BF -------------------------
                b_n = to_blas_int(9)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call dcopy(b_n, dmf, b_incx, dmf2, b_incy)
                b_n = to_blas_int(9)
                b_incx = to_blas_int(1)
                call dscal(b_n, wgt, dmf2, b_incx)
                call utctab('CUMU', 3, 12, 4, dmf2, &
                            bf, gm, xab1, gmefl)
            else if (.not. dri) then
!           ----- CALCUL DU PRODUIT BMT.DMF.BF -------------------------
                b_n = to_blas_int(9)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call dcopy(b_n, dmf, b_incx, dmf2, b_incy)
                b_n = to_blas_int(9)
                b_incx = to_blas_int(1)
                call dscal(b_n, wgt, dmf2, b_incx)
                call utctab('CUMU', 3, 12, 8, dmf2, &
                            bf, bm, xab1, mefl)
            else
                ASSERT(ASTER_FALSE)
            end if
        end if
    end do
!
    if (option .eq. 'RIGI_MECA') then
        if (.not. dri) then
            call dxqloc(flex, memb, mefl, ctor, matrRigi)
        else if (dri) then
!     Add rotational to stiffness matrix
            ctor = 0.d0
            call dxqloc(flex, memb, mefl, ctor, matrRigi)
            call dxqlocdri1(gmemb, matrRigi)
            call dxqlocdri2(btgmemb, matrRigi)
            call dxqlocdri3(gmefl, matrRigi)
            call dxqlocdri4(bxb, matrRigi)
        else
            ASSERT(ASTER_FALSE)
        end if
    else if (option .eq. 'EPOT_ELEM') then
        call jevech('PDEPLAR', 'L', jvDisp)
        call utpvgl(4, 6, pgl, zr(jvDisp), depl)
        if (.not. dri) then
            call dxqloe(flex, memb, mefl, ctor, coupmf, &
                        depl, ener)
        else if (dri) then
!     Add rotational to stiffness matrix
            ctor = 0.d0
            call dxqloc(flex, memb, mefl, ctor, matrRigi)
            call dxqlocdri1(gmemb, matrRigi)
            call dxqlocdri2(btgmemb, matrRigi)
            call dxqlocdri3(gmefl, matrRigi)
            call dxqlocdri4(bxb, matrRigi)
            call dxqloe_NV(coupmf, matrRigi, depl, ener)
        else
            ASSERT(ASTER_FALSE)
        end if
        ! call bsthpl(plateCara, plateOrie, &
        !             jvGeom, nomte, xyzl, &
        !             bsigth)
        ! if (indith) then
        !     enerTher = 0.d0
        !     do i = 1, 24
        !         enerTher = enerTher+depl(i)*bsigth(i)
        !     end do
        !     ener(1) = ener(1)-enerTher
        ! end if
    end if
!
    if (present(matrRigi_)) then
        matrRigi_ = matrRigi
    end if
    if (present(ener_)) then
        ener_ = ener
    end if
!
end subroutine
