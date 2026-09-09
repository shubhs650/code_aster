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
subroutine te0412(option, nomte)
!
    use MaterialPara_module
    use MaterialPara_type
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystPara, compCoorSystPlate, &
                                isPlateTria, isPlateQuad, isPlateQ4GG, isPlateDKTG
    implicit none
!
#include "asterc/r8dgrd.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/cosiro.h"
#include "asterfort/dkqbf.h"
#include "asterfort/dkqedg.h"
#include "asterfort/dktbf.h"
#include "asterfort/dktedg.h"
#include "asterfort/dsqedg.h"
#include "asterfort/dstedg.h"
#include "asterfort/dxeffi.h"
#include "asterfort/dxefro.h"
#include "asterfort/dxmate.h"
#include "asterfort/dxqbm.h"
#include "asterfort/dxtbm.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/gquad4.h"
#include "asterfort/gtria3.h"
#include "asterfort/jevech.h"
#include "asterfort/jquad4.h"
#include "asterfort/pmrvec.h"
#include "asterfort/q4gedg.h"
#include "asterfort/r8inir.h"
#include "asterfort/t3gedg.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/utpvgl.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: DKT/DKTG/DST/Q4G/Q4GG
!
! Options: ENEL_ELGA/ENEL_ELEM
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nnomx = 4, nbsm = 3, npgmx = 4
    integer(kind=8) :: nbsig, ndim
    real(kind=8) :: pgl(3, 3)
    real(kind=8) :: eps(3), khi(3), gam(2)
    real(kind=8) :: bf(3, 3*nnomx), bm(3, 2*nnomx), um(2, nnomx), uf(3, nnomx)
    real(kind=8) :: ul(6, nnomx), qsi, eta, xyzl(3, 4), jacob(5), poids
    real(kind=8) :: cara(25)
    real(kind=8) :: nmm(nbsm), mff(nbsm)
    real(kind=8) :: enelm(npgmx), enelf(npgmx)
    real(kind=8) :: enelt(npgmx), enelc(npgmx), enemf(npgmx)
    real(kind=8) :: ent, enm, enf, enc, enmf
    real(kind=8) :: effint(32), effort(32), degpg(32)
    real(kind=8) :: dmeps(3), dfkhi(3), dcgam(3)
    real(kind=8) :: df(9), dm(9), dmf(9), dc(4), dci(4)
    real(kind=8) :: dmc(3, 2), dfc(3, 2)
    integer(kind=8) :: nno, nnoel, npg, ipoids, icoopg
    integer(kind=8) :: jvGeom, kpg, ino, jvDisp, isig, jsig, jvEner, iret
    integer(kind=8) :: icompo, icontp, jvari, nbvar, ivpg
    integer(kind=8) :: multic
    character(len=16), pointer :: compor(:) => null()
    character(len=16) :: optio2, relaComp, relaFlua
    aster_logical ::  lKitDDI, coupmf
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Get plate parameters
    call getCara(plateCara, plateOrie)

    nbsig = 6
    if (isPlateDKTG(plateCara) .or. isPlateQ4GG(plateCara)) then
        nbsig = 8
    end if
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnoel, npg=npg, &
                     jpoids=ipoids, jcoopg=icoopg)

! - Calculate the transformation: global coordinate system/intrinsic coordinate system
    call compCoorSystPara(plateCara, zr(jvGeom), pgl)

! - Compute coordinate system for plate
    call compCoorSystPlate(pgl, plateCara, plateOrie)

! - Change coordinates of geometry
    call utpvgl(nno, 3, pgl, zr(jvGeom), xyzl)

! - Compute geometric parametres of plate
    if (isPlateQuad(plateCara)) then
        call gquad4(xyzl, cara)
    elseif (isPlateTria(plateCara)) then
        call gtria3(xyzl, cara)
    else
        ASSERT(ASTER_FALSE)
    end if

! --- INITIALISATION
    call r8inir(npgmx, 0.d0, enelt, 1)
    call r8inir(npgmx, 0.d0, enelm, 1)
    call r8inir(npgmx, 0.d0, enelf, 1)
    call r8inir(npgmx, 0.d0, enelc, 1)
    call r8inir(npgmx, 0.d0, enemf, 1)
    ent = 0.0d0
    enm = 0.0d0
    enf = 0.0d0
    enc = 0.0d0
    enmf = 0.0d0

! - ON REGARDE SI ON EST EN LINEAIRE OU ENN NON-LINEAIRE
    call tecach('NNO', 'PCOMPOR', 'L', iret, iad=icompo)
    if (iret .eq. 0) then
        call jevech('PCOMPOR', 'L', vk16=compor)
        relaComp = compor(RELA_NAME)
        relaFlua = compor(CREEP_NAME)
        lKitDDI = relaComp(1:7) .eq. 'KIT_DDI'
        if (relaComp(1:4) .eq. 'ELAS' .or. relaComp(1:4) .eq. 'ENDO' .or. &
            relaComp(1:6) .eq. 'MAZARS' .or. relaComp(1:7) .eq. 'GLRC_DM' .or. &
            relaComp(1:11) .eq. 'GLRC_DAMAGE' .or. &
            (lKitDDI .and. relaFlua .eq. 'GLRC_DM')) then
            if (option .eq. 'ENEL_ELGA') then
                call jevech('PDEPLAR', 'L', jvDisp)
                if (.not. (isPlateDKTG(plateCara) .or. isPlateQ4GG(plateCara))) then
                    call cosiro(plateCara, plateOrie, &
                                'PCONTRR', 'L', 'UI', 'G', &
                                icontp)
                else
                    call jevech('PCONTRR', 'L', icontp)
                end if
            else if (option .eq. 'ENEL_ELEM') then
                call jevech('PDEPLR', 'L', jvDisp)
                if (.not. (isPlateDKTG(plateCara) .or. isPlateQ4GG(plateCara))) then
                    call cosiro(plateCara, plateOrie, &
                                'PCONTPR', 'L', 'UI', 'G', &
                                icontp)
                else
                    call jevech('PCONTPR', 'L', icontp)
                end if
            else
                ASSERT(ASTER_FALSE)
            end if
            if (((isPlateDKTG(plateCara) .or. isPlateQ4GG(plateCara)) .and. lKitDDI) .or. &
                relaComp .eq. 'GLRC_DAMAGE') then
                if (option .eq. 'ENEL_ELGA') then
                    call jevech('PVARIGR', 'L', jvari)
                else if (option .eq. 'ENEL_ELEM') then
                    call jevech('PVARIPR', 'L', jvari)
                end if
            end if

            if ((.not. lKitDDI) .or. &
                (.not. (isPlateDKTG(plateCara) .or. isPlateQ4GG(plateCara)))) then
! ------------- Change coordinates of displacements
                call utpvgl(nno, 6, pgl, zr(jvDisp), ul)
!
!       -- PARTITION DU DEPLACEMENT EN MEMBRANE/FLEXION :
!       -------------------------------------------------
                do ino = 1, nnoel
                    um(1, ino) = ul(1, ino)
                    um(2, ino) = ul(2, ino)
                    uf(1, ino) = ul(3, ino)
                    uf(2, ino) = ul(5, ino)
                    uf(3, ino) = -ul(4, ino)
                end do
            end if
!
!     -- CALCUL DES CONTRAINTES GENERALISEES :
!     -------------------------------------------------
            if (isPlateDKTG(plateCara) .or. isPlateQ4GG(plateCara)) then
                do kpg = 1, npg
                    do isig = 1, nbsig
                        effort((kpg-1)*nbsig+isig) = zr(icontp-1+(kpg-1)*8+isig)
                    end do
                end do
                call dxefro(npg, plateOrie%t2ui, effort, effint)
            else
                call dxeffi(plateCara, plateOrie, &
                            option, nomte, zr(icontp), nbsig, &
                            effint)
            end if

            do kpg = 1, npg
!
                qsi = zr(icoopg-1+ndim*(kpg-1)+1)
                eta = zr(icoopg-1+ndim*(kpg-1)+2)
                if (isPlateQuad(plateCara)) then
                    call jquad4(xyzl, qsi, eta, jacob)
                    poids = zr(ipoids+kpg-1)*jacob(1)
                    call dxqbm(qsi, eta, jacob(2), bm)
                    call dkqbf(qsi, eta, jacob(2), cara, bf)
                else
                    poids = zr(ipoids+kpg-1)*cara(7)
                    call dxtbm(cara(9), bm)
                    call dktbf(qsi, eta, cara, bf)
                end if
!
                if ((isPlateDKTG(plateCara) .or. isPlateQ4GG(plateCara)) .and. lKitDDI) then
                    read (zk16(icompo-1+NVAR), '(I16)') nbvar
                    ivpg = jvari+(kpg-1)*nbvar+24
                    do isig = 1, nbsm
                        eps(isig) = zr(ivpg+isig)
                        khi(isig) = zr(ivpg+isig+3)
                    end do
                else
!
!         -- CALCUL DE EPS, KHI :
!         -----------------------------------
                    call pmrvec('ZERO', 3, 2*nnoel, bm, um, eps)
                    call pmrvec('ZERO', 3, 3*nnoel, bf, uf, khi)
                    if (relaComp .eq. 'GLRC_DAMAGE') then
                        read (zk16(icompo-1+NVAR), '(I16)') nbvar
                        ivpg = jvari+(kpg-1)*nbvar-1
                        do isig = 1, nbsm
                            eps(isig) = eps(isig)-zr(ivpg+isig)
                            khi(isig) = khi(isig)-zr(ivpg+isig+3)
                        end do
                    end if
                end if
!
!  --    CALCUL DE LA DENSITE D'ENERGIE POTENTIELLE ELASTIQUE :
!        ==========================================================
                if ((option .eq. 'ENEL_ELGA') .or. (option .eq. 'ENEL_ELEM')) then
!
!  --      DENSITE D'ENERGIE POTENTIELLE ELASTIQUE AU POINT
!  --      D'INTEGRATION COURANT
!          ---------------------
                    call r8inir(nbsm, 0.d0, nmm, 1)
                    call r8inir(nbsm, 0.d0, mff, 1)
!
                    do isig = 1, nbsm
                        nmm(isig) = effint((kpg-1)*nbsig+isig)
                        mff(isig) = effint((kpg-1)*nbsig+isig+3)
                    end do
!
                    do jsig = 1, nbsm
                        enelm(kpg) = enelm(kpg)+0.5d0*nmm(jsig)*eps(jsig)
                        enelf(kpg) = enelf(kpg)+0.5d0*mff(jsig)*khi(jsig)
                    end do
                    enelt(kpg) = enelm(kpg)+enelf(kpg)
!
                    enm = enm+enelm(kpg)*poids
                    enf = enf+enelf(kpg)*poids
                    ent = ent+enelt(kpg)*poids
                end if
            end do
        end if

    else
!
        if (option .eq. 'ENEL_ELGA') then
            call jevech('PDEPLAR', 'L', jvDisp)
        else if (option .eq. 'ENEL_ELEM') then
            call jevech('PDEPLR', 'L', jvDisp)
        end if
! ----- Change coordinates of displacements
        call utpvgl(nno, 6, pgl, zr(jvDisp), ul)
!
        call dxmate(plateCara, plateOrie, &
                    'RIGI', df, dm, dmf, dc, &
                    dci, dmc, dfc, &
                    multic, coupmf)
!
!     -- CALCUL DES DEFORMATIONS GENERALISEES AUX POINTS DE GAUSS
!     -----------------------------------------------------------
        optio2 = 'DEGE_ELGA'
        if (nomte .eq. 'MEDKTR3' .or. nomte .eq. 'MEDKTG3') then
            call dktedg(plateCara, plateOrie, &
                        xyzl, optio2, ul, &
                        degpg, multic)
        else if (nomte .eq. 'MEDSTR3') then
            call dstedg(plateCara, plateOrie, &
                        xyzl, optio2, ul, &
                        degpg)
        else if (nomte .eq. 'MEDKQU4' .or. nomte .eq. 'MEDKQG4') then
            call dkqedg(plateCara, plateOrie, &
                        xyzl, optio2, ul, &
                        degpg)
        else if (nomte .eq. 'MEDSQU4') then
            call dsqedg(plateCara, plateOrie, &
                        xyzl, optio2, ul, &
                        degpg)
        else if (nomte .eq. 'MEQ4QU4' .or. nomte .eq. 'MEQ4GG4') then
            call q4gedg(plateCara, plateOrie, &
                        xyzl, optio2, ul, &
                        degpg)
        else if (nomte .eq. 'MET3TR3' .or. nomte .eq. 'MET3GG3') then
            call t3gedg(plateCara, plateOrie, &
                        xyzl, optio2, ul, &
                        degpg)
        end if

        do kpg = 1, npg
            qsi = zr(icoopg-1+ndim*(kpg-1)+1)
            eta = zr(icoopg-1+ndim*(kpg-1)+2)
            if (isPlateQuad(plateCara)) then
                call jquad4(xyzl, qsi, eta, jacob)
                poids = zr(ipoids+kpg-1)*jacob(1)
            else
                poids = zr(ipoids+kpg-1)*cara(7)
            end if

!  --    CALCUL DE LA DENSITE D'ENERGIE POTENTIELLE ELASTIQUE
            if ((option .eq. 'ENEL_ELGA') .or. (option .eq. 'ENEL_ELEM')) then
                do isig = 1, nbsm
                    eps(isig) = degpg((kpg-1)*8+isig)
                    khi(isig) = degpg((kpg-1)*8+isig+3)
                end do
                do isig = 1, 2
                    gam(isig) = degpg((kpg-1)*8+isig+6)
                end do
!
! --- CALCUL DES PRODUITS :
!           MEMBRANE     : [DM]{EPSI}
!           FLEXION      : [DF]{KHI}
!           CISAILLEMENT : [DC]{GAM}
!
                eps(3) = eps(3)*2.d0
                khi(3) = khi(3)*2.d0
!
                call pmrvec('ZERO', 3, 3, dm, eps, dmeps)
                call pmrvec('ZERO', 3, 3, df, khi, dfkhi)
                call pmrvec('ZERO', 2, 2, dc, gam, dcgam)
!
                do isig = 1, nbsm
                    enelm(kpg) = enelm(kpg)+0.5d0*eps(isig)*dmeps(isig)
                    enelf(kpg) = enelf(kpg)+0.5d0*khi(isig)*dfkhi(isig)
                end do
                do isig = 1, 2
                    enelc(kpg) = enelc(kpg)+0.5d0*gam(isig)*dcgam(isig)
                end do
!
! --- COUPLAGE MEMBRANE - FLEXION (ELAS_COQUE)
!
                if (coupmf) then
                    call pmrvec('ZERO', 3, 3, dmf, eps, dmeps)
                    call pmrvec('ZERO', 3, 3, dmf, khi, dfkhi)
!
                    do isig = 1, nbsm
                        enemf(kpg) = enemf(kpg)+ &
                                     0.5d0*(eps(isig)*dfkhi(isig)+khi(isig)*dmeps(isig))
                    end do
                end if
!
                enelt(kpg) = enelm(kpg)+enelf(kpg)+enelc(kpg)+enemf(kpg)
                enm = enm+enelm(kpg)*poids
                enf = enf+enelf(kpg)*poids
                enc = enc+enelc(kpg)*poids
                enmf = enmf+enemf(kpg)*poids
                ent = ent+enelt(kpg)*poids
            end if
        end do
    end if
!
! ---- RECUPERATION DU CHAMP DES DENSITES D'ENERGIE DE DEFORMATION
! ---- ELASTIQUE EN SORTIE
!      -------------------
    if (option .eq. 'ENEL_ELGA') then
        call jevech('PENERDR', 'E', jvEner)
    else if (option .eq. 'ENEL_ELEM') then
        call jevech('PENERD1', 'E', jvEner)
    end if

    if (option .eq. 'ENEL_ELGA') then
        do kpg = 1, npg
            zr(jvEner-1+(kpg-1)*5+1) = enelt(kpg)
            zr(jvEner-1+(kpg-1)*5+2) = enelm(kpg)
            zr(jvEner-1+(kpg-1)*5+3) = enelf(kpg)
            zr(jvEner-1+(kpg-1)*5+4) = enelc(kpg)
            zr(jvEner-1+(kpg-1)*5+5) = enemf(kpg)
        end do
    else if (option .eq. 'ENEL_ELEM') then
        zr(jvEner) = ent
        zr(jvEner+1) = enm
        zr(jvEner+2) = enf
        zr(jvEner+3) = enc
        zr(jvEner+4) = enmf
    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
