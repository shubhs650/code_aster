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
subroutine ccvrrl(mesh, model, caraElem, &
                  lRestCell, nbRestCell, restCellJv, &
                  fieldElnoS, codret)
!
    implicit none
!
#include "asterc/indik8.h"
#include "asterc/r8rddg.h"
#include "asterf_types.h"
#include "asterfort/angvec.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/carces.h"
#include "asterfort/cccmcr.h"
#include "asterfort/cncinv.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/utpvlg.h"
#include "jeveux.h"
!
    character(len=8), intent(in) :: mesh, model, caraElem
    aster_logical, intent(in) :: lRestCell
    integer(kind=8), intent(in) :: nbRestCell
    character(len=24), intent(in) :: restCellJv
    character(len=19), intent(in) :: fieldElnoS
    integer(kind=8), intent(out) :: codret
!
! --------------------------------------------------------------------------------------------------
!
!  CALC_CHAMP - VERIFICATION DES REPERES LOCAUX
!
! --------------------------------------------------------------------------------------------------
!
!  ROUTINE SERVANT A VERIFIER L'ORIENTATION DES REPERES LOCAUX
!    LORS DU PASSAGE ELNO -> NOEU
!
! IN  :
!   NOMMAI  K8   NOM DU MAILLAGE A VERIFIER
!   MODELE  K8   NOM DU MODELE
!   CARAEL  K8   NOM DU CARAELEM
!   MESMAI  K24  NOM DU VECTEUR CONTENANT LES MAILLES SUR LESQUELLES
!                LE CALCUL EST DEMANDE
!   CHAMES  K19  NOM DU CHAM_ELEM_S POUR LEQUEL ON VERIFIE LES REPERES
!   CMPERR  K1   COMPORTEMENT EN CAS DE PROBLEME
!     'F' : EMISSION D'UNE ERREUR <F>
!     'A' : EMISSION D'UNE ALARME POUR PREVENIR L'UTILISATEUR
!     ' ' : SILENCE => CODE RETOUR
!
! OUT :
!   CODRET  I    CODE RETOUR
!     0 SI OK
!     1 EN CAS DE PROBLEME
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), parameter :: maxtol = 8.7266463d-2
    integer(kind=8) :: jvRestCell, nbCell, ier, nbNodeMesh, jconi1, jconi2, iCell, ncmax
    integer(kind=8) :: posit, posma, cellNume, iCell2, cellNume2, iNodeMesh, jcesd, jcesl
    integer(kind=8) :: jcesv, iexori, jrepe, iepais, nbCellMesh, idir
    integer(kind=8) :: jconx1, jconx2, jcoord, jcesdc, ialpha, ibeta
    integer(kind=8) :: jalpha, jbeta, jgamma, jcesc, jcescc, jcesdd, jceslc, jcesvc
    integer(kind=8) :: adcar1(3), adcar2(3)
    real(kind=8) :: maxdif, angle1, angle2
    real(kind=8) :: pgl(3, 3), vl(3), vg1(3), vg2(3), vg3(3), vg4(3)
    character(len=16) :: modeli
    character(len=19) :: cnxinv, carsd, modelLigrel, carcc
    aster_logical :: lprobm
    integer(kind=8), pointer :: dime(:) => null()
    real(kind=8), pointer :: workVect(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    codret = 0
    carsd = '&&CCVRRL.CARORIEN'
    carcc = '&&CCVRRL.CARCOQUE'
    cnxinv = '&&CCVRRL.CNCINV'
    call dismoi("NOM_LIGREL", model, "MODELE", repk=modelLigrel)
    call jeveuo(modelLigrel//'.REPE', 'L', jrepe)

! - Access to mesh
    call jeveuo(mesh//'.DIME', 'L', vi=dime)
    nbCellMesh = dime(3)

! - Access to restricted list of cells
    if (lRestCell) then
        call jeveuo(restCellJv, 'L', jvRestCell)
    else
        jvRestCell = 1
        ASSERT(nbRestCell .eq. 0)
    end if
    call jeveuo(fieldElnoS//'.CESD', 'L', jcesdd)

!   CONVERSION DE LA CARTE D'ORIENTATION EN UN CHAM_ELEM_S
    call exisd('CARTE', caraElem//'.CARORIEN', iexori)
    if (iexori .eq. 1) then
        call carces(caraElem//'.CARORIEN', 'ELEM', ' ', 'V', carsd, ' ', ier)
        call jeveuo(carsd//'.CESD', 'L', jcesd)
        call jeveuo(carsd//'.CESL', 'L', jcesl)
        call jeveuo(carsd//'.CESV', 'L', jcesv)
        call jeveuo(carsd//'.CESC', 'L', jcesc)
        call jelira(carsd//'.CESC', 'LONMAX', ncmax)
        jalpha = indik8(zk8(jcesc), 'ALPHA   ', 1, ncmax)
        jbeta = indik8(zk8(jcesc), 'BETA    ', 1, ncmax)
        jgamma = indik8(zk8(jcesc), 'GAMMA   ', 1, ncmax)
    else
        jcesd = 0
        jcesl = 0
        jcesv = 0
        jcesc = 0
    end if

!   CONVERSION DE LA CARTE CARACTERISTIQUE DES COQUES EN UN CHAM_ELEM_S
    call exisd('CARTE', caraElem//'.CARCOQUE', iexori)
    if (iexori .eq. 1) then
        call carces(caraElem//'.CARCOQUE', 'ELEM', ' ', 'V', carcc, ' ', ier)
        call jeveuo(carcc//'.CESD', 'L', jcesdc)
        call jeveuo(carcc//'.CESL', 'L', jceslc)
        call jeveuo(carcc//'.CESV', 'L', jcesvc)
        call jeveuo(carcc//'.CESC', 'L', jcescc)
        call jelira(carcc//'.CESC', 'LONMAX', ncmax)
        iepais = indik8(zk8(jcescc), 'EP      ', 1, ncmax)
        ialpha = indik8(zk8(jcescc), 'ALPHA   ', 1, ncmax)
        ibeta = indik8(zk8(jcescc), 'BETA    ', 1, ncmax)
    else
        jcesdc = 0
        jceslc = 0
        jcesvc = 0
        jcescc = 0
    end if

! - Working vector
    AS_ALLOCATE(vr=workVect, size=6*nbCellMesh)

! - Inverse connectivity
    if (lRestCell) then
        call cncinv(mesh, zi(jvRestCell), nbRestCell, 'V', cnxinv)
    else
        call cncinv(mesh, [0], 0, 'V', cnxinv)
    end if
    call jelira(cnxinv, 'NUTIOC', nbNodeMesh)
    call jeveuo(jexnum(cnxinv, 1), 'L', jconi1)
    call jeveuo(jexatr(cnxinv, 'LONCUM'), 'L', jconi2)
    call jeveuo(jexnum(mesh//'.CONNEX', 1), 'L', jconx1)
    call jeveuo(jexatr(mesh//'.CONNEX', 'LONCUM'), 'L', jconx2)
    call jeveuo(mesh//'.COORDO    .VALE', 'L', jcoord)
!
    adcar1(1) = jcesd
    adcar1(2) = jcesl
    adcar1(3) = jcesv
    adcar2(1) = jcesdc
    adcar2(2) = jceslc
    adcar2(3) = jcesvc
    lprobm = .false.

!   BOUCLE SUR TOUS LES NOEUDS DU MAILLAGE
    maxdif = 0.d0
    do iNodeMesh = 1, nbNodeMesh
        nbCell = zi(jconi2+iNodeMesh)-zi(jconi2+iNodeMesh-1)
        posit = zi(jconi2+iNodeMesh-1)

        do iCell = 1, nbCell
            posma = zi(jconi1+posit+iCell-2)
            if (posma .eq. 0) cycle

! --------- Get index of cell
            if (lRestCell) then
                cellNume = zi(jvRestCell+posma-1)
            else
                cellNume = posma
            end if
            if (cellNume .eq. 0) cycle
!
            if ((workVect(6*(cellNume-1)+1) .eq. 0.d0) .and. &
                (workVect(6*(cellNume-1)+2) .eq. 0.d0) .and. &
                (workVect(6*(cellNume-1)+3) .eq. 0.d0) .and. &
                (workVect(6*(cellNume-1)+4) .eq. 0.d0) .and. &
                (workVect(6*(cellNume-1)+5) .eq. 0.d0) .and. &
                (workVect(6*(cellNume-1)+6) .eq. 0.d0)) then
                call cccmcr(jcesdd, cellNume, jrepe, jconx2, jconx1, &
                            jcoord, adcar1, adcar2, ialpha, ibeta, &
                            iepais, jalpha, jbeta, jgamma, modelLigrel, &
                            iNodeMesh, pgl, modeli, ier)
                if (ier .eq. 3) cycle
                if (ier .eq. 1) then
                    call utmess('A', 'MODELISA10_5', si=cellNume)
                end if
!
                vl(1) = 1.d0
                vl(2) = 0.d0
                vl(3) = 0.d0
                call utpvlg(1, 3, pgl, vl, vg1)
                vl(1) = 0.d0
                vl(2) = 1.d0
                vl(3) = 0.d0
                call utpvlg(1, 3, pgl, vl, vg2)
!               sauvegarde de la valeur trouvee sauf pour les coques 3d
                if (modeli .ne. 'CQ3') then
                    do idir = 1, 3
                        workVect(6*(cellNume-1)+idir) = vg1(idir)
                        workVect(3+6*(cellNume-1)+idir) = vg2(idir)
                    end do
                end if
            else
                do idir = 1, 3
                    vg1(idir) = workVect(6*(cellNume-1)+idir)
                    vg2(idir) = workVect(3+6*(cellNume-1)+idir)
                end do
            end if
!
!           Compare les repères des autres mailles liées au noeud iNodeMesh
            do iCell2 = iCell+1, nbCell
                posma = zi(jconi1+posit+iCell2-2)
                if (posma .eq. 0) cycle
! ------------- Get index of cell
                if (lRestCell) then
                    cellNume2 = zi(jvRestCell+posma-1)
                else
                    cellNume2 = posma
                end if
                if (cellNume2 .eq. 0) cycle
!
                if ((workVect(6*(cellNume2-1)+1) .eq. 0.d0) .and. &
                    (workVect(6*(cellNume2-1)+2) .eq. 0.d0) .and. &
                    (workVect(6*(cellNume2-1)+3) .eq. 0.d0) .and. &
                    (workVect(6*(cellNume2-1)+4) .eq. 0.d0) .and. &
                    (workVect(6*(cellNume2-1)+5) .eq. 0.d0) .and. &
                    (workVect(6*(cellNume2-1)+6) .eq. 0.d0)) then
!
                    call cccmcr(jcesdd, cellNume2, jrepe, jconx2, jconx1, &
                                jcoord, adcar1, adcar2, ialpha, ibeta, &
                                iepais, jalpha, jbeta, jgamma, modelLigrel, &
                                iNodeMesh, pgl, modeli, ier)
                    if (ier .eq. 3) cycle
                    if (ier .eq. 1) then
                        call utmess('A', 'MODELISA10_5', si=cellNume2)
                    end if
!
                    vl(1) = 1.d0
                    vl(2) = 0.d0
                    vl(3) = 0.d0
                    call utpvlg(1, 3, pgl, vl, vg3)
                    vl(1) = 0.d0
                    vl(2) = 1.d0
                    vl(3) = 0.d0
                    call utpvlg(1, 3, pgl, vl, vg4)
!                   sauvegarde de la valeur trouvee sauf pour les coques3d
                    if (modeli .ne. 'COQUE_3D') then
                        do idir = 1, 3
                            workVect(6*(cellNume2-1)+idir) = vg3(idir)
                            workVect(3+6*(cellNume2-1)+idir) = vg4(idir)
                        end do
                    end if
!
                else
                    do idir = 1, 3
                        vg3(idir) = workVect(6*(cellNume2-1)+idir)
                        vg4(idir) = workVect(3+6*(cellNume2-1)+idir)
                    end do
                end if
!
                call angvec(vg1, vg3, angle1)
                call angvec(vg2, vg4, angle2)
                if (angle1 .gt. maxtol .or. angle2 .gt. maxtol) then
                    maxdif = max(angle1, maxdif)
                    maxdif = max(angle2, maxdif)
                    lprobm = .true.
                end if
            end do
        end do
    end do
!
    if (lprobm) then
        call utmess('A', 'UTILITAI_4', sr=maxdif*r8rddg())
        codret = 1
    end if
!
    AS_DEALLOCATE(vr=workVect)
    call jedetr(cnxinv)
    call detrsd('CHAM_ELEM_S', carsd)
    call detrsd('CHAM_ELEM_S', carcc)
!
    call jedema()
!
end subroutine
