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
subroutine te0331(option, nomte)
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/epdcp.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/psvari.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/tecach.h"
#include "asterfort/vpri2d.h"
#include "asterfort/Behaviour_type.h"
!
    character(len=*) :: option, nomte
!     FONCTION REALISEE :
!
!         CALCUL DE LA CONTRAINTE DE WEIBULL D'UNE STRUCTURE
!         EN COMPORTEMENT NON-LINEAIRE.
!         ELEMENTS ISOPARAMETRIQUES 2D.
!
!         OPTION : 'WEIBULL'
!
! ENTREE  --->  OPTION : NOM DE L'OPTION DE CALCUL
!         --->  NOMTE  : NOM DU TYPE D'ELEMENT
!
!     ------------------------------------------------------------------
!
    integer(kind=8) :: icodre(4)
    integer(kind=8) :: codres
    character(len=8) :: fami
    character(len=16) :: nomres(4), optcal(12)
    character(len=32) :: phenom
    real(kind=8) :: sig(6), sigi, dsigwb, valres(4), epsgi
    real(kind=8) :: poids, r, volume, volact, dvol, seuil, m, v0
    real(kind=8) :: cong(4), epsq(4), dfdx(9), dfdy(9), pp, ppt
    real(kind=8) :: tc(6), tdp(6), sigold, signew, sref, tg, tmoy
!
    integer(kind=8) :: nno, kp, npg, k, ii, iweib, jtab(7), nnos, jgano, ndim
    integer(kind=8) :: idefg, issopt, ipopp, ipoppt
    integer(kind=8) :: ipoids, ivf, idfde, imate
    integer(kind=8) :: igeom, icong, ivarig
    integer(kind=8) :: isigie, isigis, nbvari
    aster_logical :: laxi
    character(len=16), pointer :: compor(:) => null()
    character(len=16) :: rela_comp
!     ------------------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: iret
!-----------------------------------------------------------------------
    fami = 'RIGI'
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
    poids = 0.d0
    dsigwb = 0.d0
    volume = 0.d0
    volact = 0.d0
    dvol = 0.d0
    laxi = .false.
    if (lteatt('AXIS', 'OUI')) laxi = .true.
!
    nomres(1) = 'M'
    nomres(2) = 'VOLU_REFE'
    nomres(3) = 'SEUIL_EPSP_CUMU'
    nomres(4) = 'SIGM_REFE'
!
    call jevech('PMATERC', 'L', imate)
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PCONTRG', 'L', icong)
    call jevech('PVARIPG', 'L', ivarig)
    call jevech('PSOUSOP', 'L', issopt)
    call jevech('PDOMMAG', 'L', isigie)
    call jevech('PWEIBUL', 'E', iweib)
    call jevech('PSIGISG', 'E', isigis)
!
    call tecach('OOO', 'PVARIPG', 'L', iret, nval=7, itab=jtab)
    nbvari = max(jtab(6), 1)*jtab(7)
    call jevech('PCOMPOR', 'L', vk16=compor)
    rela_comp = compor(RELA_NAME)
!
    call psvari(rela_comp, nbvari, ipopp, ipoppt)
!
    optcal(1) = zk24(issopt) (1:16)
    optcal(2) = zk24(issopt) (17:19)
!
    cong = 0.d0
    epsq = 0.d0

! -FONCTION SEUIL
    ppt = 0.d0
    pp = 0.d0
!
!
!     --- CAS SIGU DEPEND DE LA TEMPERATURE (WEIBULL_FO)?
!     SI OUI ET QU IL N Y A PAS DE CHAMP DE TEMPERATURE
!     ARRET
!
    call rccoma(zi(imate), 'WEIBULL', 1, phenom, codres)
!
    if (optcal(1) .eq. 'SIGM_ELMOY') then
        tmoy = 0.d0
    end if
!
!
!     --- RECUPERATION DES DONNEES MATERIAU ---
!
    call rcvalb(fami, 1, 1, '+', zi(imate), &
                ' ', phenom, 0, ' ', [0.d0], &
                3, nomres, valres, icodre, 1)
    call rcvalb(fami, 1, 1, '+', zi(imate), &
                ' ', phenom, 0, ' ', [0.d0], &
                1, nomres(3), valres(3), icodre(3), 1)
    if (icodre(3) .ne. 0) valres(3) = 1.d-6
    m = valres(1)
    v0 = valres(2)
    seuil = valres(3)
!
!     --- BOUCLE SUR POINTS DE GAUSS SUIVANT OPTIONS DE CALCUL ---
!
!=================================================================
!=================================================================
    if ((optcal(1) .eq. 'SIGM_ELMOY') .and. (optcal(2) .eq. 'NON')) then
        do kp = 1, npg
            k = (kp-1)*nno
            r = 0.d0
            call dfdm2d(nno, kp, ipoids, idfde, zr(igeom), &
                        poids, dfdx, dfdy)
            if (laxi) then
                do ii = 1, nno
                    r = r+zr(igeom+2*ii-2)*zr(ivf+k+ii-1)
                end do
                poids = poids*r
            end if
! VOLUME PLASTIFIE
            pp = zr(ivarig+nbvari*(kp-1)+ipopp-1)
            if (pp .ge. seuil) then
                dvol = poids
                volume = volume+dvol
                do ii = 1, 4
                    cong(ii) = cong(ii)+dvol*zr(icong+4*kp+ii-5)
                end do
!           --- TEMPERATURE MOYENNE
                call rcvarc(' ', 'TEMP', '+', 'RIGI', kp, 1, tg, iret)
                if (iret .ne. 0) tg = 0.d0
                tmoy = tmoy+tg*dvol
            end if
! VOLUME PLASTIQUE ACTIF
            if (rela_comp .eq. 'LEMAITRE' .and. (pp .ge. seuil)) then
                ppt = 1.d0
            else
                ppt = zr(ivarig+nbvari*(kp-1)+ipoppt-1)
            end if
            if (ppt .ge. (1.0d0)) then
                dvol = poids
                volact = volact+dvol
            end if
        end do
!
        sigi = 0.d0
        if ((volact .ne. 0.d0) .and. (volume .ne. 0.d0)) then
            sig(1) = cong(1)/volume
            sig(2) = cong(2)/volume
            sig(3) = cong(3)/volume
            sig(4) = cong(4)/volume
            call vpri2d(sig, sigi)
!
            tmoy = tmoy/volume
            call rcvalb(fami, 1, 1, '+', zi(imate), &
                        ' ', phenom, 1, 'TEMP', [tmoy], &
                        1, nomres(4), valres(4), icodre(4), 1)
            sref = valres(4)
            sigi = sigi/sref
        end if
        sigold = zr(isigie)
        if (sigi .gt. sigold) then
            zr(isigis) = sigi
        else
            zr(isigis) = zr(isigie)
        end if
        sigi = zr(isigis)
!
        dsigwb = volume/v0*(sigi**m)
!=================================================================
!=================================================================
    elseif ((optcal(1) .eq. 'SIGM_ELGA') .and. (optcal(2) .eq. 'OUI')) then
        do kp = 1, npg
            r = 0.d0
            k = (kp-1)*nno
            call dfdm2d(nno, kp, ipoids, idfde, zr(igeom), &
                        poids, dfdx, dfdy)
            if (laxi) then
                do ii = 1, nno
                    r = r+zr(igeom+2*ii-2)*zr(ivf+k+ii-1)
                end do
                poids = poids*r
            end if
            volume = poids
            call jevech('PDEFORR', 'L', idefg)
            do ii = 1, 4
                cong(ii) = zr(icong+4*kp+ii-5)
                epsq(ii) = zr(idefg+4*kp+ii-5)
            end do
            pp = zr(ivarig+nbvari*(kp-1)+ipopp-1)
            if (rela_comp .eq. 'LEMAITRE' .and. (pp .ge. seuil)) then
                ppt = 1.d0
            else
                ppt = zr(ivarig+nbvari*(kp-1)+ipoppt-1)
            end if
!
            signew = 0.d0
            if (ppt .ge. (1.d0)) then
!
!              ------CALCUL DE SIGI ET EPSI---------
!
                tc(1) = cong(1)
                tc(2) = cong(2)
                tc(3) = cong(3)
                tc(4) = cong(4)
                tc(5) = 0.d0
                tc(6) = 0.d0
!
                tdp(1) = epsq(1)
                tdp(2) = epsq(2)
                tdp(3) = epsq(3)
                tdp(4) = epsq(4)
                tdp(5) = 0.d0
                tdp(6) = 0.d0
                call epdcp(tc, tdp, sigi, epsgi)
                call rcvalb(fami, kp, 1, '+', zi(imate), &
                            ' ', phenom, 0, ' ', [0.d0], &
                            1, nomres(4), valres(4), icodre(4), 1)
                sref = valres(4)
!
                signew = exp((-epsgi/2.d0))*sigi/sref
            end if
            sigold = zr(isigie+kp-1)
            if (signew .gt. sigold) then
                zr(isigis+kp-1) = signew
            else
                zr(isigis+kp-1) = zr(isigie+kp-1)
            end if
            signew = zr(isigis+kp-1)
            dsigwb = dsigwb+volume*(signew**m)/v0
!
        end do
!=================================================================
!=================================================================
    elseif ((optcal(1) .eq. 'SIGM_ELMOY') .and. (optcal(2) .eq. 'OUI')) then
!
        do kp = 1, npg
            r = 0.d0
            k = (kp-1)*nno
            call dfdm2d(nno, kp, ipoids, idfde, zr(igeom), &
                        poids, dfdx, dfdy)
            if (laxi) then
                do ii = 1, nno
                    r = r+zr(igeom+2*ii-2)*zr(ivf+k+ii-1)
                end do
                poids = poids*r
            end if
! VOL PLASTIFIE
            pp = zr(ivarig+nbvari*(kp-1)+ipopp-1)
            call jevech('PDEFORR', 'L', idefg)
            if (pp .ge. seuil) then
                dvol = poids
                volume = volume+dvol
                do ii = 1, 4
                    cong(ii) = cong(ii)+dvol*zr(icong+4*kp+ii-5)
                    epsq(ii) = epsq(ii)+dvol*zr(idefg+4*kp+ii-5)
                end do
!           --- TEMPERATURE MOYENNE
                call rcvarc(' ', 'TEMP', '+', 'RIGI', kp, &
                            1, tg, iret)
                if (iret .ne. 0) tg = 0.d0
                tmoy = tmoy+tg*dvol
            end if
! VOL PLASTIQUE ACTIF
            if (rela_comp .eq. 'LEMAITRE' .and. (pp .ge. seuil)) then
                ppt = 1.d0
            else
                ppt = zr(ivarig+nbvari*(kp-1)+ipoppt-1)
            end if
            if (ppt .ge. (1.0d0)) then
                dvol = poids
                volact = volact+dvol
            end if
        end do
!
        signew = 0.d0
        if ((volact .ne. (0.d0)) .and. (volume .ne. 0.d0)) then
            tc(1) = cong(1)/volume
            tc(2) = cong(2)/volume
            tc(3) = cong(3)/volume
            tc(4) = cong(4)/volume
            tc(5) = 0.d0
            tc(6) = 0.d0
!
            tdp(1) = epsq(1)/volume
            tdp(2) = epsq(2)/volume
            tdp(3) = epsq(3)/volume
            tdp(4) = epsq(4)/volume
            tdp(5) = 0.d0
            tdp(6) = 0.d0
            call epdcp(tc, tdp, sigi, epsgi)
            tmoy = tmoy/volume
            call rcvalb(fami, 1, 1, '+', zi(imate), &
                        ' ', phenom, 1, 'TEMP', [tmoy], &
                        1, nomres(4), valres(4), icodre(4), 1)
            sref = valres(4)
            signew = exp((-epsgi/2.d0))*sigi/sref
        end if
        sigold = zr(isigie)
        if (signew .gt. sigold) then
            zr(isigis) = signew
        else
            zr(isigis) = zr(isigie)
        end if
        signew = zr(isigis)
        dsigwb = volume*(signew**m)/v0
!=================================================================
!=================================================================
    elseif ((optcal(1) .eq. 'SIGM_ELGA') .and. (optcal(2) .eq. 'NON')) then
        do kp = 1, npg
            k = (kp-1)*nno
            r = 0.d0
            do ii = 1, 4
                cong(ii) = zr(icong+(4*kp)-5+ii)
            end do
            call dfdm2d(nno, kp, ipoids, idfde, zr(igeom), &
                        poids, dfdx, dfdy)
            if (laxi) then
                do ii = 1, nno
                    r = r+zr(igeom+2*ii-2)*zr(ivf+k+ii-1)
                end do
                poids = poids*r
            end if
            volume = poids
!
            sigi = 0.d0
            pp = zr(ivarig+nbvari*(kp-1)+ipopp-1)
            if (rela_comp .eq. 'LEMAITRE' .and. (pp .ge. seuil)) then
                ppt = 1.d0
            else
                ppt = zr(ivarig+nbvari*(kp-1)+ipoppt-1)
            end if
            if (ppt .ge. (1.d0)) then
! CALCUL DE SIGI
                call vpri2d(cong, sigi)
                call rcvalb(fami, kp, 1, '+', zi(imate), &
                            ' ', phenom, 0, ' ', [0.d0], &
                            1, nomres(4), valres(4), icodre(4), 1)
                sref = valres(4)
                sigi = sigi/sref
            end if
            sigold = zr(isigie+kp-1)
            if (sigi .gt. sigold) then
                zr(isigis+kp-1) = sigi
            else
                zr(isigis+kp-1) = zr(isigie+kp-1)
            end if
            sigi = zr(isigis+kp-1)
            dsigwb = dsigwb+volume*(sigi**m)/v0
        end do
    else
!        OPTION DE CALCUL NON VALIDE
        ASSERT(.false.)
    end if
!=================================================================
!=================================================================
    zr(iweib) = dsigwb
!
!     DESTRUCTION DES OBJETS CREES DANS LA BASE
!
end subroutine
