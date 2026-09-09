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
subroutine te0338(option, nomte)
    implicit none
#include "jeveux.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/assert.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/epdcp.h"
#include "asterfort/fgequi.h"
#include "asterfort/jevech.h"
#include "asterfort/psvari.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/tecach.h"
!
    character(len=*) :: option, nomte
!     FONCTION REALISEE :
!
!         CALCUL DU CHAM_ELEM DE WEIBULL
!         COMPORTEMENT NON-LINEAIRE.
!         ELEMENTS ISOPARAMETRIQUES 3D.
!
!         OPTION : 'WEIBULL'
!
! ENTREE  --->  OPTION : NOM DE L'OPTION DE CALCUL
!         --->  NOMTE  : NOM DU TYPE D'ELEMENT
!
!
!-DEL CHARACTER*32 JEXNUM,JEXNOM,JEXATR,JEXR8
!
    integer(kind=8) :: icodre(4)
    integer(kind=8) :: codres
    character(len=8) :: fami
    character(len=32) :: phenom
    character(len=16) :: optcal(12), nomres(4)
!
    real(kind=8) :: sigm(6), sig1, sigwk, valres(4), epsg(6), eps1
    real(kind=8) :: m, vref, sref, seuil, dvpg, poids, vkp
    real(kind=8) :: equi(6), pp, ppt, vkpact
    real(kind=8) :: sigold, signew, tg, tmoy
!
    integer(kind=8) :: i, kp, ndim, nbvari, ipopp, ipoppt
    integer(kind=8) :: jgano, ipoids, ivf, idfde, npg, nno, nnos
    integer(kind=8) :: imate, igeom, icong, ivarig, issopt, iweib, idefg, nbvp
    integer(kind=8) :: isigie, isigis, jtab(7), iret
    character(len=16), pointer :: compor(:) => null()
    character(len=16) :: rela_comp
!
!======================== CORPS DU PROGRAMME ===========================
!
!
    fami = 'RIGI'
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
    nbvp = 3
!
!     1.2 CHAMPS IN
!     -------------
    call jevech('PMATERC', 'L', imate)
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PCONTRG', 'L', icong)
    call jevech('PVARIPG', 'L', ivarig)
    call jevech('PSOUSOP', 'L', issopt)
    call jevech('PDOMMAG', 'L', isigie)
    call tecach('OOO', 'PVARIPG', 'L', iret, nval=7, &
                itab=jtab)
    nbvari = max(jtab(6), 1)*jtab(7)
    call jevech('PCOMPOR', 'L', vk16=compor)
    rela_comp = compor(RELA_NAME)
!
    call psvari(rela_comp, nbvari, ipopp, ipoppt)
!
!     1.3 CHAMPS OUT
!     --------------
    call jevech('PWEIBUL', 'E', iweib)
    call jevech('PSIGISG', 'E', isigis)
!
!     1.4 OPTIONS DE CALCUL
!     ---------------------
    optcal(1) = zk24(issopt) (1:16)
    optcal(2) = zk24(issopt) (17:19)
!
!     1.5 DONNES DE LA RC DU MATERIAU
!     -------------------------------
    nomres(1) = 'M'
    nomres(2) = 'VOLU_REFE'
    nomres(3) = 'SEUIL_EPSP_CUMU'
    nomres(4) = 'SIGM_REFE'
!
    call rccoma(zi(imate), 'WEIBULL', 1, phenom, codres)
!
!     --- S'IL N Y A PAS DE CHAMP DE TEMPERATURE
!     ARRET
!
    if (optcal(1) .eq. 'SIGM_ELMOY') then
        tmoy = 0.d0
    end if
!
    call rcvalb(fami, 1, 1, '+', zi(imate), &
                ' ', phenom, 0, ' ', [0.d0], &
                3, nomres, valres, icodre, 1)
    call rcvalb(fami, 1, 1, '+', zi(imate), &
                ' ', phenom, 0, ' ', [0.d0], &
                1, nomres(3), valres(3), icodre(3), 1)
    if (icodre(3) .ne. 0) valres(3) = 1.0d-6
    m = valres(1)
    vref = valres(2)
    seuil = valres(3)
!
!     1.6 INITIALISATION
!     ------------------
    poids = 0.d0
    dvpg = 0.d0
    vkp = 0.d0
    vkpact = 0.d0
    sigwk = 0.d0
    sigm = 0.d0
    epsg = 0.d0
! -CRITERE PLASTIQUE
    ppt = 0.d0
    pp = 0.d0
!
!
!     2. BOUCLE SUR POINTS DE GAUSS SUIVANT OPTIONS DE CALCUL
!     -------------------------------------------------------
!     2.1 SIGM_W A PARTIR DE SIGM MOYENNE SANS CORRECTION PLASTIQUE
!     -------------------------------------------------------------
    if ((optcal(1) .eq. 'SIGM_ELMOY') .and. (optcal(2) .eq. 'NON')) then
!        2.1.1 INTEGRATION DE SIGM SUR LA PARTIE PLASTIFIEE
!        --------------------------------------------------
        do kp = 1, npg, 1
! VOLUME PLASTIFIE
            pp = zr(ivarig+nbvari*(kp-1)+ipopp-1)
            if (pp .ge. seuil) then
                call dfdm3d(nno, kp, ipoids, idfde, zr(igeom), &
                            poids)
                dvpg = poids
                vkp = vkp+dvpg
                do i = 1, 6, 1
                    sigm(i) = sigm(i)+dvpg*zr(icong+6*kp+i-7)
                end do
!           --- TEMPERATURE MOYENNE
                call rcvarc(' ', 'TEMP', '+', 'RIGI', kp, &
                            1, tg, iret)
                if (iret .ne. 0) tg = 0.d0
                tmoy = tmoy+tg*dvpg
            end if
! VOLUME PLASTIQUE ACTIF
            if (rela_comp .eq. 'LEMAITRE' .and. (pp .ge. seuil)) then
                ppt = 1.d0
            else
                ppt = zr(ivarig+nbvari*(kp-1)+ipoppt-1)
            end if
            if (ppt .ge. (1.d0)) then
                dvpg = poids
                vkpact = vkpact+dvpg
            end if
!
        end do
!
        sig1 = 0.d0
        if ((vkp .ne. 0.0d0) .and. (vkpact .ne. 0.d0)) then
!           2.1.2 CALCUL DE LA VALEUR MOYENNE DE SIGM SUR LA  MAILLE
!           --------------------------------------------------------
            do i = 1, 6, 1
                sigm(i) = sigm(i)/vkp
            end do
!
            tmoy = tmoy/vkp
            call rcvalb('RIGI', 1, 1, '+', zi(imate), &
                        ' ', phenom, 1, 'TEMP', [tmoy], &
                        1, nomres(4), valres(4), icodre(4), 1)
            sref = valres(4)
!
!           2.1.3 CALCUL DE SIGM_W
!           ----------------------
            call fgequi(sigm, 'SIGM', nbvp, equi)
            sig1 = max(equi(3), equi(4), equi(5))
            sig1 = sig1/sref
        end if
!
        sigold = zr(isigie)
        if (sig1 .gt. sigold) then
            zr(isigis) = sig1
        else
            zr(isigis) = zr(isigie)
        end if
        sig1 = zr(isigis)
!
        sigwk = (vkp/vref)*(sig1**m)
!
!     2.2 SIGM_W A PARTIR DE SIGM MOYENNE AVEC CORRECTION PLASTIQUE
!     -------------------------------------------------------------
    else if ((optcal(1) .eq. 'SIGM_ELMOY') .and. (optcal(2) .eq. 'OUI')) &
        then
!        2.2.1 INTEGRATION DE SIGM SUR LA PARTIE PLASTIFIEE
!        --------------------------------------------------
        call jevech('PDEFORR', 'L', idefg)
        do kp = 1, npg, 1
! VOLUME PLASTIFIE
            pp = zr(ivarig+nbvari*(kp-1)+ipopp-1)
            if (pp .ge. seuil) then
                call dfdm3d(nno, kp, ipoids, idfde, zr(igeom), &
                            poids)
                dvpg = poids
                vkp = vkp+dvpg
                do i = 1, 6, 1
                    sigm(i) = sigm(i)+dvpg*zr(icong+6*kp+i-7)
                    epsg(i) = epsg(i)+dvpg*zr(idefg+6*kp+i-7)
                end do
!           --- TEMPERATURE AU PG
                call rcvarc(' ', 'TEMP', '+', 'RIGI', kp, &
                            1, tg, iret)
                if (iret .ne. 0) tg = 0.d0
                tmoy = tmoy+tg*dvpg
            end if
! VOLUME PLASTIQUE ACTIF
            if (rela_comp .eq. 'LEMAITRE' .and. (pp .ge. seuil)) then
                ppt = 1.d0
            else
                ppt = zr(ivarig+nbvari*(kp-1)+ipoppt-1)
            end if
            if (ppt .ge. (1.d0)) then
                dvpg = poids
                vkpact = vkpact+dvpg
            end if
        end do
!
        signew = 0.d0
        if ((vkp .ne. 0.0d0) .and. (vkpact .ne. 0.d0)) then
!           2.2.2 CALCUL DE LA VALEUR MOYENNE DE SIGM SUR LA  MAILLE
!           --------------------------------------------------------
            do i = 1, 6, 1
                sigm(i) = sigm(i)/vkp
                epsg(i) = epsg(i)/vkp
            end do
!
            tmoy = tmoy/vkp
            call rcvalb('RIGI', 1, 1, '+', zi(imate), &
                        ' ', phenom, 1, 'TEMP', [tmoy], &
                        1, nomres(4), valres(4), icodre(4), 1)
            sref = valres(4)
!
!           2.2.3 CALCUL DE SIGM_W
!           ----------------------
            call epdcp(sigm, epsg, sig1, eps1)
            signew = (sig1/sref)*exp(-eps1*0.5d0)
        end if
!
        sigold = zr(isigie)
        if (signew .gt. sigold) then
            zr(isigis) = signew
        else
            zr(isigis) = zr(isigie)
        end if
        signew = zr(isigis)
!
        sigwk = (vkp/vref)*(signew**m)
!
!     2.3 SIGM_W A PARTIR DE SIGM ORIGINAL AVEC CORRECTION PLASTIQUE
!     -------------------------------------------------------------
    else if ((optcal(1) .eq. 'SIGM_ELGA') .and. (optcal(2) .eq. 'OUI')) &
        then
        call jevech('PDEFORR', 'L', idefg)
        do kp = 1, npg, 1
            pp = zr(ivarig+nbvari*(kp-1)+ipopp-1)
            signew = 0.d0
            if (pp .ge. seuil) then
                call dfdm3d(nno, kp, ipoids, idfde, zr(igeom), &
                            poids)
                dvpg = poids
                if ((rela_comp .eq. 'LEMAITRE') .and. (pp .ge. seuil)) then
                    ppt = 1.d0
                else
                    ppt = zr(ivarig+nbvari*(kp-1)+ipoppt-1)
                end if
                if (ppt .ge. (1.d0)) then
                    do i = 1, 6, 1
                        sigm(i) = zr(icong+6*kp+i-7)
                        epsg(i) = zr(idefg+6*kp+i-7)
                    end do
                    call epdcp(sigm, epsg, sig1, eps1)
                    call rcvalb(fami, kp, 1, '+', zi(imate), &
                                ' ', phenom, 0, ' ', [0.d0], &
                                1, nomres(4), valres(4), icodre(4), 1)
                    sref = valres(4)
                    signew = (sig1/sref)*exp(-eps1*0.5d0)
                end if
            end if
            sigold = zr(isigie+kp-1)
            if (signew .gt. sigold) then
                zr(isigis+kp-1) = signew
            else
                zr(isigis+kp-1) = zr(isigie+kp-1)
            end if
            signew = zr(isigis+kp-1)
            sigwk = sigwk+(dvpg/vref)*(signew**m)
        end do
!
!     2.4 SIGM_W A PARTIR DE SIGM ORIGINAL SANS CORRECTION
!     ----------------------------------------------------
    else if ((optcal(1) .eq. 'SIGM_ELGA') .and. (optcal(2) .eq. 'NON')) &
        then
        do kp = 1, npg, 1
            pp = zr(ivarig+nbvari*(kp-1)+ipopp-1)
            if (pp .ge. seuil) then
                call dfdm3d(nno, kp, ipoids, idfde, zr(igeom), &
                            poids)
                dvpg = poids
                if ((rela_comp .eq. 'LEMAITRE') .and. (pp .ge. seuil)) then
                    ppt = 1.d0
                else
                    ppt = zr(ivarig+nbvari*(kp-1)+ipoppt-1)
                end if
                sig1 = 0.d0
                if (ppt .ge. (1.d0)) then
                    do i = 1, 6, 1
                        sigm(i) = zr(icong+6*kp+i-7)
                    end do
                    call fgequi(sigm, 'SIGM', nbvp, equi)
                    sig1 = max(equi(3), equi(4), equi(5))
                    call rcvalb(fami, kp, 1, '+', zi(imate), &
                                ' ', phenom, 0, ' ', [0.d0], &
                                1, nomres(4), valres(4), icodre(4), 1)
                    sref = valres(4)
                    sig1 = sig1/sref
                end if
            else
                sig1 = 0.d0
            end if
            sigold = zr(isigie+kp-1)
            if (sig1 .gt. sigold) then
                zr(isigis+kp-1) = sig1
            else
                zr(isigis+kp-1) = zr(isigie+kp-1)
            end if
            sig1 = zr(isigis+kp-1)
!
            sigwk = sigwk+(dvpg/vref)*(sig1**m)
        end do
!
!     2.5 TRAITEMENT DES OPTIONS INVALIDES
!     ------------------------------------
    else
!       OPTION DE CALCUL NON VALIDE
        ASSERT(.false.)
    end if
!
!
!     3. ECRITURE DU CHAM_ELEM DE WEIBULL
!     -----------------------------------
    zr(iweib) = sigwk
!
end subroutine
