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
subroutine te0576(option, nomte)
!
    use Behaviour_module
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterc/r8vide.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/enelpg.h"
#include "asterfort/eps1mc.h"
#include "asterfort/jevech.h"
#include "asterfort/nbsigm.h"
#include "asterfort/nmgeom.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
#include "MeshTypes_type.h"
#include "FE_module.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 3D
!
! Options: ENEL_ELGA, ETOT_ELEM, ETOT_ELGA
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = "RIGI"
    integer(kind=8), parameter :: ksp = 1
    integer(kind=8), parameter :: mxcmel = 162, nbSigx = 6, nbEpsi = 6
    real(kind=8) :: epsiCurr(mxcmel), epsiPrev(mxcmel)
    real(kind=8), parameter :: zero = 0.d0, undemi = 0.5d0, deux = 2.d0
    real(kind=8), parameter :: nharm = 0.d0
    aster_logical, parameter :: axi = ASTER_FALSE
    integer(kind=8) :: jvDispCurr, jvDispPrev, jvVari
    integer(kind=8) :: jvSigmCurr, jvSigmPrev, jvGeom, jvMaterc, jvTime
    integer(kind=8) :: jvGaussWeight, jvBaseFunc, jvDBaseFunc
    integer(kind=8) :: nbsig, ndim, nno, npg, nbvari
    integer(kind=8) :: idenem, idener
    real(kind=8) :: enerElem, poids
    integer(kind=8) :: jtab(7), iret
    integer(kind=8) :: kpg, isigm, iEpsi
    real(kind=8) :: epsiKpgCurr(nbEpsi), epsiKpgPrev(nbEpsi), epsiKpgDelta(nbEpsi)
    real(kind=8) :: time
    real(kind=8) :: enerKpg(MT_NNOMAX), r
    real(kind=8) :: f(3, 3)
    real(kind=8) :: sigmKpgPrev(nbSigx), sigmKpgCurr(nbSigx)
    real(kind=8) :: integ1, integ2, integ
    real(kind=8) :: epsiDumm(6), dfdbid(MAX_BV_CG)
    character(len=16) :: relaName, defoComp
    aster_logical :: largeStrain
    type(Material_Para) :: materPara
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, npg=npg, &
                     jpoids=jvGaussWeight, jvf=jvBaseFunc, jdfde=jvDBaseFunc)
!
    nbsig = nbsigm()
    ASSERT(nno .le. MT_NNOMAX)
    ASSERT(nbsig .le. 6)
    ASSERT(npg .le. 27)
    enerElem = zero
    enerKpg = zero

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

    if (option(1:4) .eq. 'ENEL') then
! ----- Material parameters
        call jevech('PMATERC', 'L', jvMaterc)

! ----- Initializations of material parameters on current cell
        call initParaCell(fami, zi(jvMaterc), materPara)

! ----- Set local coordinate system from user
        call getUserLCS(ndim, nno, jvGeom, materPara%lcsPara)

! ----- Current displacements
        call jevech('PDEPLAR', 'L', jvDispCurr)

! ----- Current stresses
        call jevech('PCONTRR', 'L', jvSigmCurr)

! ----- Get time
        time = r8vide()
        call tecach('ONO', 'PINSTR', 'L', iret, iad=jvTime)
        if (jvTime .ne. 0) then
            time = zr(jvTime)
        end if

! ----- Internal state variables (only for non-linear)
        jvVari = 1
        nbvari = 0
        call tecach('ONO', 'PVARIGR', 'L', iret, nval=7, itab=jtab)
        if (iret .eq. 0) then
            jvVari = jtab(1)
            nbvari = max(jtab(6), 1)*jtab(7)
        end if

    end if

! - Get parameters from behaviour (linear and non-linear cases)
    call behaviourGetParameters(relaName, defoComp)
    if ((defoComp .eq. 'SIMO_MIEHE') .or. (defoComp .eq. 'GDEF_LOG')) then
        largeStrain = ASTER_TRUE
    else
        largeStrain = ASTER_FALSE
    end if

! - CAS DU CALCUL DE LA DENSITE D'ENERGIE TOTALE
    if (option(1:4) .eq. 'ETOT') then
        if (largeStrain) then
            call utmess('F', 'ENERGY1_7', sk=defoComp)
        end if

! ----- RECUPERATION DU CHAMP DE DEPLACEMENT A L'INSTANT COURANT
        call jevech('PDEPLR', 'L', jvDispCurr)

! ----- RECUPERATION EVENTUELLE DU CHAMP DE DEPLACEMENT A L'INSTANT PRECEDENT
        call tecach('NNO', 'PDEPLM', 'L', iret, iad=jvDispPrev)

! ----- RECUPERATION DU CHAMP DE CONTRAINTES A L'INSTANT COURANT
        call jevech('PCONTPR', 'L', jvSigmCurr)

! ----- RECUPERATION EVENTUELLE DU CHAMP DE CONTRAINTES A L'INSTANT PRECEDENT
        call tecach('NNO', 'PCONTMR', 'L', iret, iad=jvSigmPrev)

! ----- CALCUL DU CHAMP DE DEFORMATIONS AU PREMIER ORDRE (current)
        call eps1mc(nno, ndim, nbsig, npg, &
                    jvGaussWeight, jvBaseFunc, jvDBaseFunc, &
                    zr(jvGeom), zr(jvDispCurr), nharm, &
                    epsiCurr)

! ----- CALCUL DU CHAMP DE DEFORMATIONS AU PREMIER ORDRE (previous)
        if (jvDispPrev .ne. 0) then
            call eps1mc(nno, ndim, nbsig, npg, &
                        jvGaussWeight, jvBaseFunc, jvDBaseFunc, &
                        zr(jvGeom), zr(jvDispPrev), nharm, &
                        epsiPrev)
        end if
    end if

    do kpg = 1, npg
! ----- CALCUL DU JACOBIEN AU POINT D'INTEGRATION COURANT
        call dfdm3d(nno, kpg, jvGaussWeight, jvDBaseFunc, zr(jvGeom), poids)
        epsiKpgCurr = 0.d0

        if (option(1:4) .eq. 'ENEL') then
! --------- Initializations of material parameters on current integration point
            call initParaPoin(kpg, ksp, materPara)

! --------- TENSEUR DES CONTRAINTES AU POINT D'INTEGRATION COURANT :
            sigmKpgCurr = 0.d0
            do iSigm = 1, nbsig
                sigmKpgCurr(iSigm) = zr(jvSigmCurr+(kpg-1)*nbsig+iSigm-1)
            end do

! --------- CALCUL DU JACOBIEN AU POINT D'INTEGRATION COURANT
            call nmgeom(3, nno, axi, largeStrain, zr(jvGeom), kpg, &
                        jvGaussWeight, jvBaseFunc, jvDBaseFunc, &
                        zr(jvDispCurr), &
                        ASTER_TRUE, poids, dfdbid, f, epsiDumm, &
                        r)

! --------- CALCUL DE L'ENERGIE DE DEFORMATION ELASTIQUE
            call enelpg(materPara, time, &
                        relaName, defoComp, &
                        f, sigmKpgCurr, &
                        nbvari, zr(jvVari+(kpg-1)*nbvari), &
                        enerKpg(kpg))

        else if (option(1:4) .eq. 'ETOT') then
! --------- TENSEURS DES DEFORMATIONS  ET DES CONTRAINTES AU PAS DE
! --------- TEMPS COURANT ET AU PAS DE TEMPS PRECEDENT S'IL Y A LIEU,
! --------- AU POINT D'INTEGRATION COURANT
            do iEpsi = 1, nbEpsi
                epsiKpgCurr(iEpsi) = epsiCurr(iEpsi+(kpg-1)*nbEpsi)
                if (jvDispPrev .ne. 0) then
                    epsiKpgPrev(iEpsi) = epsiPrev(iEpsi+(kpg-1)*nbEpsi)
                end if
            end do
            do iSigm = 1, nbSig
                sigmKpgCurr(iSigm) = zr(jvSigmCurr+(kpg-1)*nbsig+iSigm-1)
                if (jvSigmPrev .ne. 0) then
                    sigmKpgPrev(iSigm) = zr(jvSigmPrev+(kpg-1)*nbsig+iSigm-1)
                end if
            end do
            epsiKpgDelta = 0.d0
            if (jvDispPrev .ne. 0) then
                do iEpsi = 1, nbEpsi
                    epsiKpgDelta(iEpsi) = epsiKpgCurr(iEpsi)-epsiKpgPrev(iEpsi)
                end do
            end if

! --------- CALCUL DES TERMES A SOMMER POUR OBTENIR LA DENSITE D'ENERGIE TOTALE
            integ1 = sigmKpgCurr(1)*epsiKpgDelta(1)+ &
                     sigmKpgCurr(2)*epsiKpgDelta(2)+ &
                     sigmKpgCurr(3)*epsiKpgDelta(3)+ &
                     deux*sigmKpgCurr(4)*epsiKpgDelta(4)+ &
                     deux*sigmKpgCurr(5)*epsiKpgDelta(5)+ &
                     deux*sigmKpgCurr(6)*epsiKpgDelta(6)

            if (jvDispPrev .ne. 0 .and. jvSigmPrev .ne. 0) then
                integ2 = sigmKpgPrev(1)*epsiKpgDelta(1)+ &
                         sigmKpgPrev(2)*epsiKpgDelta(2)+ &
                         sigmKpgPrev(3)*epsiKpgDelta(3)+ &
                         deux*sigmKpgPrev(4)*epsiKpgDelta(4)+ &
                         deux*sigmKpgPrev(5)*epsiKpgDelta(5)+ &
                         deux*sigmKpgPrev(6)*epsiKpgDelta(6)
                enerKpg(kpg) = undemi*(integ1+integ2)
            else
                integ = sigmKpgCurr(1)*epsiKpgCurr(1)+ &
                        sigmKpgCurr(2)*epsiKpgCurr(2)+ &
                        sigmKpgCurr(3)*epsiKpgCurr(3)+ &
                        deux*sigmKpgCurr(4)*epsiKpgCurr(4)+ &
                        deux*sigmKpgCurr(5)*epsiKpgCurr(5)+ &
                        deux*sigmKpgCurr(6)*epsiKpgCurr(6)
                enerKpg(kpg) = undemi*integ
            end if
            enerElem = enerElem+enerKpg(kpg)*poids
        end if
    end do

! - RECUPERATION DU CHAMP DES DENSITES D'ENERGIE DE DEFORMATION ELASTIQUE EN SORTIE
    call jevech('PENERDR', 'E', idener)
    if (option(1:4) .eq. 'ETOT') then
        call jevech('PENERDM', 'L', idenem)
        if (option .eq. 'ETOT_ELGA') then
            do kpg = 1, npg
                zr(idener+kpg-1) = zr(idenem+kpg-1)+enerKpg(kpg)
            end do
        else if (option .eq. 'ETOT_ELEM') then
            zr(idener) = zr(idenem)+enerElem
        end if
    else
        if (option .eq. 'ENEL_ELGA') then
            do kpg = 1, npg
                zr(idener+kpg-1) = enerKpg(kpg)
            end do
        end if
    end if
!
end subroutine
