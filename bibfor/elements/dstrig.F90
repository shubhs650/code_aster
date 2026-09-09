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
subroutine dstrig(plateCara, plateOrie, &
                  xyzl, option, pgl, &
                  rig, ener)
!
    use plate_type
    implicit none
!
#include "asterc/r8gaem.h"
#include "asterf_types.h"
#include "asterfort/bsthpl.h"
#include "asterfort/dstbfa.h"
#include "asterfort/dstbfb.h"
#include "asterfort/dstci2.h"
#include "asterfort/dstcis.h"
#include "asterfort/dsxhft.h"
#include "asterfort/dxmate.h"
#include "asterfort/dxtbm.h"
#include "asterfort/dxtloc.h"
#include "asterfort/dxtloe.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/gtria3.h"
#include "asterfort/jevech.h"
#include "asterfort/promat.h"
#include "asterfort/r8inir.h"
#include "asterfort/utbtab.h"
#include "asterfort/utctab.h"
#include "asterfort/utdtab.h"
#include "asterfort/utmess.h"
#include "asterfort/utpvgl.h"
#include "jeveux.h"
!
    type(plateOrie_Para), intent(in) :: plateOrie
    type(plateCara_Para), intent(in) :: plateCara
    real(kind=8) :: xyzl(3, *), pgl(*), rig(*), ener(*)
    character(len=16) :: option
!
! --------------------------------------------------------------------------------------------------
!
!     MATRICE DE RIGIDITE DE L'ELEMENT DE PLAQUE DST (AVEC CISAILLEMENT)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ndim, nnos, npg, ipoids, icoopg, ivf, idfdx, idfd2, jgano
    integer(kind=8) :: int, multic, perm(9), perm2(36)
    integer(kind=8) :: i, j, jdepg, k, k1, k2
    real(kind=8) :: wgt, aire
    real(kind=8) :: df(3, 3), dm(3, 3), dmf(3, 3), dc(2, 2), dci(2, 2)
    real(kind=8) :: dmc(3, 2), dfc(3, 2)
    real(kind=8) :: bfb(3, 9), bfa(3, 3), hft2(2, 6)
    real(kind=8) :: bca(2, 3), pb(3, 9), bm(3, 6), pm(3, 6)
    real(kind=8) :: xab1(3, 6), xab2(3, 9), xab3(2, 3), xab4(3, 3), xab5(3, 2)
    real(kind=8) :: xab6(9, 2)
    real(kind=8) :: xab7(3, 2), xab8(6, 2), xab9(3, 6)
!                   ----(9,9)  ---(9,9)  ---(9,9)
    real(kind=8) :: kf11(81), kfc(81), kfb(81), kf12(9, 3)
    real(kind=8) :: kmf11(6, 9), kmf(6, 9)
!                   ----(3,3) --(3,3) ---(3,3)
    real(kind=8) :: kf22(9), ka(9), kaa(9)
!                   ----(9,9)  -----(9,9)  ----(6,6)
    real(kind=8) :: flex(81), flexi(81), memb(36), memexc(36)
!                   ----(6,9)  -----(6,9)
!
    real(kind=8) :: mefl(6, 9), mefli(6, 9), depl(18)
    real(kind=8) :: kfc11(9, 3), kfc21(9), kmc(6, 3)
    real(kind=8) :: kmf12(6, 3), kmf12a(36)
    real(kind=8) :: enerth, excent, un, zero
    real(kind=8) :: qsi, eta, carat3(21)
    aster_logical :: coupmf, exce, ismultic
!     ------------------------------------------------------------------
    real(kind=8) :: ctor
    data perm/1, 4, 7, 2, 5, 8, 3, 6, 9/
    data perm2/1, 7, 13, 19, 25, 31,&
     &             2, 8, 14, 20, 26, 32,&
     &             3, 9, 15, 21, 27, 33,&
     &             4, 10, 16, 22, 28, 34,&
     &             5, 11, 17, 23, 29, 35,&
     &             6, 12, 18, 24, 30, 36/
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', ndim=ndim, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jcoopg=icoopg, jvf=ivf, jdfde=idfdx, jdfd2=idfd2, &
                     jgano=jgano)
!
    zero = 0.0d0
    un = 1.0d0
    call r8inir(27, zero, kfc11, 1)
    call r8inir(9, zero, kfc21, 1)
    call r8inir(18, zero, kmc, 1)
    call r8inir(18, zero, pm, 1)
    call r8inir(18, zero, kmf12, 1)
    call r8inir(36, zero, kmf12a, 1)
    enerth = zero
!
!
    ctor = plateCara%coefRigiDRZ
    excent = plateCara%offset
    exce = (abs(excent) .gt. un/r8gaem())
    ismultic = .false.
!
!     ----- CALCUL DES GRANDEURS GEOMETRIQUES SUR LE TRIANGLE --------
    call gtria3(xyzl, carat3)
!
!     ----- CALCUL DES MATRICES DE RIGIDITE DU MATERIAU EN FLEXION,
!           MEMBRANE ET CISAILLEMENT INVERSEE
    call dxmate(plateCara, plateOrie, &
                'RIGI', df, dm, dmf, dc, &
                dci, dmc, dfc, &
                multic, coupmf)
!   VERIFICATION CAS EXCENTREMENT MULTICOUCHES
    if (multic .gt. 0) ismultic = .true.
    if (exce .and. ismultic) then
        call utmess('F', 'PLATE1_15')
    end if
!     ------------------------------------------------------------------
!     CALCUL DE LA MATRICE DE RIGIDITE DE L'ELEMENT EN MEMBRANE
!     ------------------------------------------------------------------
!
!     ------ CALCUL DE LA MATRICE BM -----------------------------------
    call dxtbm(carat3(9), bm)
!     ------ CALCUL DU PRODUIT BMT.DM.BM -------------------------------
    call utbtab('ZERO', 3, 6, dm, bm, &
                xab1, memb)
    aire = carat3(8)
    do k = 1, 36
        memb(k) = memb(k)*aire
    end do
!
!     ------------------------------------------------------------------
!     CALCUL DES MATRICES DE RIGIDITE DE L'ELEMENT EN FLEXION ET
!     COUPLAGE MEMBRANE/FLEXION
!     ------------------------------------------------------------------
!
!     ------- CALCUL DE LA MATRICE BFB -------------------------------
    call dstbfb(carat3(9), bfb)
!
!     ------- CALCUL DU PRODUIT BFBT.DF.BFB --------------------------
    call utbtab('ZERO', 3, 9, df, bfb, &
                xab2, kf11)
!
    if (coupmf) then
!        ----- CALCUL DU PRODUIT BMT.DMF.BFB ------------------------
        call utctab('ZERO', 3, 9, 6, dmf, &
                    bfb, bm, xab2, kmf11)
    end if
!
!     ------- CALCUL DU PRODUIT HF.T2 ----------------------------------
    call dsxhft(df, carat3(9), hft2)
!
    if (exce) then
!
! ---   CALCUL DES MATRICES BCA ,PB ET PM DANS LE CAS DE L'EXCENTREMENT:
!       ---------------------------------------------------------------
        call dstci2(dci, carat3, hft2, dfc, dmc, &
                    bca, pb, pm)
    else
!
! ---   CALCUL DES MATRICES BCA ET PB DANS LE CAS NON EXCENTRE :
!       ------------------------------------------------------
        call dstcis(dci, carat3, hft2, bca, pb)
!
    end if
!
!     ------- CALCUL DU PRODUIT BCAT.DCI.BCA ---------------------------
    call utbtab('ZERO', 2, 3, dci, bca, &
                xab3, kaa)
!
    do k = 1, 81
        flex(k) = 0.d0
    end do
    do i = 1, 6
        do j = 1, 9
            mefl(i, j) = 0.d0
        end do
    end do
!
    do int = 1, npg
!
! ---   COORDONNEES DU POINT D'INTEGRATION COURANT :
!       ------------------------------------------
        qsi = zr(icoopg-1+ndim*(int-1)+1)
        eta = zr(icoopg-1+ndim*(int-1)+2)
!
! ---   CALCUL DE LA MATRICE BFA AU POINT QSI ETA :
!       -----------------------------------------
        call dstbfa(qsi, eta, carat3, bfa)
!
! ---   CALCUL DU PRODUIT BFBT.DF.BFA :
!       -----------------------------
        call utctab('ZERO', 3, 3, 9, df, &
                    bfa, bfb, xab4, kf12)
!
! ---   CALCUL DU PRODUIT BFAT.DF.BFA :
!       -----------------------------
        call utbtab('ZERO', 3, 3, df, bfa, &
                    xab4, kf22)
!
!=========================================
! ---   CAS DU COMPORTEMENT ELAS_COQUE   =
!=========================================
!
! ---   CALCUL DU PRODUIT BFBT.DFC.DCI.BCA :
!       ----------------------------------
        call utdtab('ZERO', 3, 2, 2, 9, &
                    dfc, dci, bfb, xab5, xab6)
        call promat(xab6, 9, 9, 2, bca, &
                    2, 2, 3, kfc11)
!
! ---   CALCUL DU PRODUIT BFAT.DFC.DCI.BCA :
!       ----------------------------------
        call utdtab('ZERO', 3, 2, 2, 3, &
                    dfc, dci, bfa, xab5, xab7)
        call promat(xab7, 3, 3, 2, bca, &
                    2, 2, 3, kfc21)
!
! ---   CALCUL DU PRODUIT BMT.DMC.DCI.BCA :
!       ---------------------------------
        call utdtab('ZERO', 3, 2, 2, 6, &
                    dmc, dci, bm, xab5, xab8)
        call promat(xab8, 6, 6, 2, bca, &
                    2, 2, 3, kmc)
!
!==============================================================
! ---   CALCUL DE LA PARTIE FLEXION DE LA MATRICE DE RIGIDITE =
! ---   FLEXI =  KF11 + KF12 * PB + PB T * KF12 T             =
! ---          + PB T * ( KF22 + KAA ) * PB                   =
!==============================================================
!
! ---   CALCUL DE LA SOMME KF22 + KAA :
!       -----------------------------
        do k = 1, 9
            ka(k) = kf22(k)+kaa(k)+kfc21(k)+kfc21(perm(k))
        end do
!
! ---   CALCUL DU PRODUIT PBT.KA.PB :
!       ---------------------------
        call utbtab('ZERO', 3, 9, ka, pb, &
                    xab2, flexi)
        do k1 = 1, 81
            kfb(k1) = 0.d0
        end do
        do i = 1, 9
            do j = 1, 9
                k1 = 9*(j-1)+i
                k2 = 9*(i-1)+j
                do k = 1, 3
                    kfb(k1) = kfb(k1)+(kf12(i, k)+kfc11(i, k))*pb(k, j)
                end do
                kfc(k2) = kfb(k1)
            end do
        end do
        do k = 1, 81
            flexi(k) = flexi(k)+kf11(k)+kfb(k)+kfc(k)
        end do
!
        wgt = zr(ipoids+int-1)*carat3(7)
        do k = 1, 81
            flex(k) = flex(k)+flexi(k)*wgt
        end do
!
        if (coupmf .or. exce) then
!
! ---     DETERMINATION DE LA MATRICE [KMF12] = [BM]T*[DMF]*[BFA]
! ---     CETTE MATRICE INTERVIENT DANS LA MATRICE DE RIGIDITE
! ---     DE MEMBRANE ET DANS LA PARTIE COUPLAGE MEMBRANE-FLEXION
! ---     DE LA RIGIDITE :
!         --------------
            call utctab('ZERO', 3, 3, 6, dmf, &
                        bfa, bm, xab4, kmf12)
!
        end if
!
        if (exce) then
!=====================================================================
! ---     L'EXCENTREMENT FAIT QUE L'ON DOIT AJOUTER A LA RIGIDITE    =
! ---     DE MEMBRANE CLASSIQUE, LES TERMES                          =
! ---     [PM]T*([KF22]+[KAA])*[PM] + [KMF12]*[PM] + [PM]T*[KMF12]T  =
!=====================================================================
!
! ---     CALCUL DU TERME [PM]T*([KF22] + [KAA])*[PM] :
!         ------------------------------------------
            call utbtab('ZERO', 3, 6, ka, pm, &
                        xab9, memexc)
!
! ---     CALCUL DU TERME [KMF12]*[PM] :
!         ----------------------------
            call promat(kmf12, 6, 6, 3, pm, &
                        3, 3, 6, kmf12a)
!
! ---     AJOUT DE CES NOUVELLES RIGIDITES DE MEMBRANE A LA
! ---     RIGIDITE DE MEMBRANE CLASSIQUE :
!         -----------------------------
            do i = 1, 36
                memb(i) = memb(i)+(memexc(i)+kmf12a(i)+kmf12a(perm2(i)))*wgt
            end do
!
        end if
!
!=====================================================================
! ---     TRAITEMENT DU COUPLAGE MEMBRANE-FLEXION                    =
! ---     L'EXCENTREMENT FAIT QUE L'ON DOIT AJOUTER A LA RIGIDITE    =
! ---     DE COUPLAGE MEMBRANE-FLEXION LES TERMES                    =
! ---     [PM]T*([KF22]+[KAA])*[PB] + [PM]T*[KF12]T + [KMF12]*[PB]   =
!=====================================================================
!
        if (coupmf .or. exce) then
!
            kmf = 0.d0
!
            if (exce) then
!
! ---       AJOUT DANS LA RIGIDITE DE COUPLAGE MEMBRANE-FLEXION
! ---       DU TERME  [PM]T*([KF22]+[KAA])*[PB] :
!           ----------------------------------
                call utctab('ZERO', 3, 9, 6, ka, &
                            pb, pm, xab2, kmf)
!
            end if
!
! ---     AJOUT DANS LA RIGIDITE DE COUPLAGE MEMBRANE-FLEXION
! ---     DES TERMES [PM]T*[KF12]T + [KMF12]*[PB]   :
!         --------------------------------------
            do i = 1, 9
                do j = 1, 6
                    do k = 1, 3
                        kmf(j, i) = kmf(j, i)+(kmf12(j, k)+kmc(j, k))*pb(k, i)+pm(k, j)*kf12(i, k)
                    end do
                    k = 6*(i-1)+j
                    mefli(j, i) = kmf11(j, i)+kmf(j, i)
                end do
            end do
            do i = 1, 6
                do j = 1, 9
                    mefl(i, j) = mefl(i, j)+mefli(i, j)*wgt
                end do
            end do
        end if
!
    end do
!
    if (option .eq. 'RIGI_MECA') then
        call dxtloc(flex, memb, mefl, ctor, rig)
!
    else if (option .eq. 'EPOT_ELEM') then
        call jevech('PDEPLAR', 'L', jdepg)
        call utpvgl(3, 6, pgl, zr(jdepg), depl)
        call dxtloe(flex, memb, mefl, ctor, coupmf, &
                    depl, ener)
        ! call bsthpl(plateCara, plateOrie, &
        !             jvGeom, nomte, xyzl, &
        !             bsigth)
        ! if (indith) then
        !     do i = 1, 18
        !         enerth = enerth+depl(i)*bsigth(i)
        !     end do
        !     ener(1) = ener(1)-enerth
        ! end if
    end if
!
end subroutine
