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
! aslint: disable=W0413
!
subroutine q4grig(plateCara, plateOrie, &
                  xyzl, option, pgl, &
                  matrRigi_, ener_)
!
    use plate_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/bsthpl.h"
#include "asterfort/dsqbfb.h"
#include "asterfort/dxmate.h"
#include "asterfort/dxqbm.h"
#include "asterfort/dxqloc.h"
#include "asterfort/dxqloe.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/gquad4.h"
#include "asterfort/jevech.h"
#include "asterfort/jquad4.h"
#include "asterfort/q4gbc.h"
#include "asterfort/utbtab.h"
#include "asterfort/utctab.h"
#include "asterfort/utdtab.h"
#include "asterfort/utpvgl.h"
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
!     MATRICE DE RIGIDITE DE L'ELEMENT Q4GAMMA (AVEC CISAILLEMENT)
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), parameter :: zero = 0.d0
    integer(kind=8) :: npg, ndim, ipoids, icoopg
    integer(kind=8) :: multic
    real(kind=8) :: wgt, depl(24)
    real(kind=8) :: df(3, 3), dm(3, 3), dmf(3, 3), dc(2, 2), dci(2, 2)
    real(kind=8) :: dmc(3, 2), dfc(3, 2)
    real(kind=8) :: bf(3, 12), bc(2, 12), bm(3, 8)
    real(kind=8) :: xab1(3, 12), xab2(2, 12), xab3(3, 8)
    real(kind=8) :: xab4(3, 12)
    real(kind=8) :: kf(144), kc(144)
    real(kind=8) :: flexi(144), flex(144)
    real(kind=8) :: membi(64), memb(64)
    real(kind=8) :: mefli(96), mefl(96), kmc(96), kfc(144)
    real(kind=8) :: enerTher, caraq4(25)
    real(kind=8) :: jacob(5), qsi, eta
    aster_logical :: coupmf
    integer(kind=8) :: i, jvDisp, k
    real(kind=8) :: ctor, excent
    real(kind=8) :: matrRigi(300), ener(3)
!
! --------------------------------------------------------------------------------------------------
!
    enerTher = zero

!
    call elrefe_info(fami='RIGI', npg=npg, ndim=ndim, &
                     jpoids=ipoids, jcoopg=icoopg)

! - Get parameters
    ctor = plateCara%coefRigiDRZ
    excent = plateCara%offset
    if (excent .ne. zero) then
        coupmf = .true.
    end if

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
    kmc = 0.d0
    kfc = 0.d0
!
    do i = 1, npg
        qsi = zr(icoopg-1+ndim*(i-1)+1)
        eta = zr(icoopg-1+ndim*(i-1)+2)
!        ---------------------------------------------------------------
!        CALCUL DE LA MATRICE DE RIGIDITE DE L'ELEMENT EN FLEXION
!        ---------------------------------------------------------------
!        ----- CALCUL DU JACOBIEN SUR LE QUADRANGLE --------------------
        call jquad4(xyzl, qsi, eta, jacob)
!        ---- CALCUL DE LA MATRICE BF ----------------------------------
        call dsqbfb(qsi, eta, jacob(2), bf)
!        ---- CALCUL DU PRODUIT BFT.DF.BF ------------------------------
        call utbtab('ZERO', 3, 12, df, bf, &
                    xab1, kf)
!        ---- CALCUL DE LA MATRICE BC ----------------------------------
        call q4gbc(qsi, eta, jacob(2), caraq4, bc)
!        ---- CALCUL DU PRODUIT BCT.DC.BC -----------------------------
        call utbtab('ZERO', 2, 12, dc, bc, &
                    xab2, kc)
!        ----- CALCUL DU PRODUIT BFT.DFC.BC ----------------------
        call utdtab('ZERO', 3, 2, 12, 12, &
                    dfc, bc, bf, xab4, kfc)
!        ----- CALCUL DE LA SOMME KF + KC = FLEXI ----------------------
        do k = 1, 144
            flexi(k) = kf(k)+kc(k)+kfc(k)
        end do
        wgt = zr(ipoids+i-1)*jacob(1)
        do k = 1, 144
            flex(k) = flex(k)+flexi(k)*wgt
        end do
!        ---------------------------------------------------------------
!        CALCUL DE LA MATRICE DE RIGIDITE DE L'ELEMENT EN MEMBRANE
!        ---------------------------------------------------------------
!        ----- CALCUL DE LA MATRICE BM ---------------------------------
        call dxqbm(qsi, eta, jacob(2), bm)
!        ----- CALCUL DU PRODUIT BMT.DM.BM -----------------------------
        call utbtab('ZERO', 3, 8, dm, bm, &
                    xab3, membi)
!        ----- CALCUL DE LA MATRICE DE RIGIDITE EN MEMBRANE ------------
        do k = 1, 64
            memb(k) = memb(k)+membi(k)*wgt
        end do
!        ----- CALCUL DU PRODUIT BMT.DMC.BC ----------------------
        call utdtab('ZERO', 3, 2, 12, 8, &
                    dmc, bc, bm, xab4, kmc)
!
        if (coupmf) then
!           ------------------------------------------------------------
!           CALCUL DES MATRICES DE COUPLAGE MEMBRANE/FLEXION
!           ------------------------------------------------------------
!           ----- CALCUL DU PRODUIT BMT.DMF.BF -------------------------
            call utctab('ZERO', 3, 12, 8, dmf, &
                        bf, bm, xab1, mefli)
            do k = 1, 96
                mefl(k) = mefl(k)+(mefli(k)+kmc(k))*wgt
            end do
        end if
    end do
!
    if (option .eq. 'RIGI_MECA') then
        call dxqloc(flex, memb, mefl, ctor, matrRigi)
    else if (option .eq. 'EPOT_ELEM') then
        call jevech('PDEPLAR', 'L', jvDisp)
        call utpvgl(4, 6, pgl, zr(jvDisp), depl)
        call dxqloe(flex, memb, mefl, ctor, coupmf, &
                    depl, ener)
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
