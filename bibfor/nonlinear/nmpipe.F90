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
subroutine nmpipe(model, ligrpi, cartyp, careta, ds_material, &
                  ds_constitutive, valinc, depdel, ddepl0, &
                  ddepl1, tau, nbeffe, eta, pilcvg, &
                  typpil, caraElem)
!
    use NonLin_Datastructure_type
    use HHO_precalc_module, only: hhoAddInputField
    use coorSyst_module, only: setOrieFields
    implicit none
!
#include "asterc/r8vide.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/celces.h"
#include "asterfort/cesexi.h"
#include "asterfort/dbgcal.h"
#include "asterfort/detrsd.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mecact.h"
#include "asterfort/megeom.h"
#include "asterfort/nmchex.h"
#include "asterfort/pipere.h"
#include "asterfort/sdmpic.h"
#include "asterfort/setStructFields.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
!
    integer(kind=8) :: pilcvg, nbeffe
    real(kind=8) :: tau, eta(2)
    character(len=24) :: typpil
    character(len=19) :: ddepl0, ddepl1
    character(len=19) :: ligrpi, cartyp, careta
    character(len=24) :: model, caraElem
    type(NL_DS_Material), intent(in) :: ds_material
    type(NL_DS_Constitutive), intent(in) :: ds_constitutive
    character(len=19) :: depdel, valinc(*)
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME - PILOTAGE)
!
! RESOLUTION DE L'EQUATION DE PILOTAGE PAR PREDICTION ELASTIQUE OU
! DEFORMATION
!
! --------------------------------------------------------------------------------------------------
!
! IN  MODELE : MODELE
! IN  LIGRPI : LIGREL DES MAILLES CONTROLEES PAR LE PILOTAGE
! IN  CARTYP : CARTE CONTENANT LE TYPE DE PILOTAGE
! In  ds_material      : datastructure for material parameters
! IN  CARELE : CARACTERISTIQUES DES ELEMENTS DE STRUCTURE
! In  ds_constitutive  : datastructure for constitutive laws management
! IN  DEPDEL : INCREMENT DE DEPLACEMENT
! IN  DDEPL0 : VARIATION DE DEPLACEMENT K-1.F0
! IN  DDEPL1 : VARIATION DE DEPLACEMENT K-1.F1
! IN  VALINC : VARIABLE CHAPEAU POUR INCREMENTS VARIABLES
! IN  TAU    : SECOND MEMBRE DE L'EQUATION DE PILOTAGE
! IN  TYPPIL : TYPE PILOTAGE : PRED_ELAS OU DEFORMATION
! OUT NBEFFE : NOMBRE DE SOLUTIONS EFFECTIVES
! OUT ETA    : ETA_PILOTAGE
! OUT PILCVG : CODE DE CONVERGENCE POUR LE PILOTAGE
!                -1 : PAS DE CALCUL DU PILOTAGE
!                 0 : CAS DU FONCTIONNEMENT NORMAL
!                 1 : PAS DE SOLUTION
!                 2 : BORNE ATTEINTE -> FIN DU CALCUL
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter ::  nbFieldOut = 1, nbFieldInMax = 100
    character(len=8) :: lpaout(nbFieldOut), lpain(nbFieldInMax)
    character(len=19) :: lchout(nbFieldOut), lchin(nbFieldInMax)
!
    integer(kind=8) :: nbCell, npg, icmp, iCell, kpg, nbgmax
    integer(kind=8) :: jcesd, jcesl, ja0a1, ja0, ja1, ja2, ja3, jtrav
    integer(kind=8) :: iret, ja4, nbFieldIn
    character(len=8), parameter :: cmpName = 'A0'
    character(len=24) :: a0a1, trav
    character(len=19) :: chgeom
    character(len=19) :: dispPrev, sigmPrev, variPrev, varcPrev
    character(len=19) :: dispCurr
    character(len=16) :: option
    real(kind=8), pointer :: cesv(:) => null()
    character(len=19), parameter :: ctau = '&&NMPIPE.CTAU'
    character(len=19), parameter :: copilo = '&&NMPIPE.COPILO'
    character(len=19), parameter :: copils = '&&NMPIPE.COPILS'
    data a0a1, trav/'&&NMPIPE.A0A1', '&&NMPIPE.TRAV'/
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    lpain = " "
    lpaout = " "
    lchin = " "
    lchout = " "
    pilcvg = 0
    if (typpil .eq. 'PRED_ELAS') then
        option = 'PILO_PRED_ELAS'
    else if (typpil .eq. 'DEFORMATION') then
        option = 'PILO_PRED_DEFO'
    else
        ASSERT(ASTER_FALSE)
    end if

! - DECOMPACTION VARIABLES CHAPEAUX
    call nmchex(valinc, 'VALINC', 'DEPMOI', dispPrev)
    call nmchex(valinc, 'VALINC', 'SIGMOI', sigmPrev)
    call nmchex(valinc, 'VALINC', 'VARMOI', variPrev)
    call nmchex(valinc, 'VALINC', 'COMMOI', varcPrev)
    call nmchex(valinc, 'VALINC', 'DEPPLU', dispCurr)
!
    call sdmpic('CHAM_ELEM', sigmPrev)
    call sdmpic('CHAM_ELEM', variPrev)

! - Get geometry field
    call megeom(model, chgeom)

! - Create field for PILOTAGE
    call detrsd('CARTE', ctau)
    call mecact('V', ctau, 'LIGREL', ligrpi, 'PILO_R', &
                ncmp=1, nomcmp=cmpName, sr=tau)

! - Add input fields
    lpain(1) = 'PGEOMER'
    lchin(1) = chgeom
    lpain(2) = 'PMATERC'
    lchin(2) = ds_material%mateco(1:19)
    lpain(3) = 'PCOMPOR'
    lchin(3) = ds_constitutive%compor(1:19)
    lpain(4) = 'PDEPLMR'
    lchin(4) = dispPrev
    lpain(5) = 'PCONTMR'
    lchin(5) = sigmPrev
    lpain(6) = 'PVARIMR'
    lchin(6) = variPrev
    lpain(7) = 'PDDEPLR'
    lchin(7) = depdel
    lpain(8) = 'PDEPL0R'
    lchin(8) = ddepl0
    lpain(9) = 'PDEPL1R'
    lchin(9) = ddepl1
    lpain(10) = 'PTYPEPI'
    lchin(10) = cartyp
    lpain(11) = 'PBORNPI'
    lchin(11) = careta
    lpain(12) = 'PCDTAU'
    lchin(12) = ctau
    lpain(14) = 'PCARCRI'
    lchin(14) = ds_constitutive%carcri(1:19)
    lpain(15) = 'PCAGNBA'
    lchin(15) = caraElem(1:8)//'.CARGENBA'
    nbFieldIn = 15

! - Add fields for orientation
    call setOrieFields(nbFieldInMax, lpain, lchin, &
                       nbFieldIn, caraElem)

! - Add HHO field
    call hhoAddInputField(model, nbFieldInMax, lchin, lpain, nbFieldIn)

! - Add output field
    lpaout(1) = 'PCOPILO'
    lchout(1) = copilo

! - Compute
    call calcul('S', option, ligrpi, &
                nbFieldIn, lchin, lpain, &
                nbFieldOut, lchout, lpaout, &
                'V', 'OUI')
    call sdmpic('CHAM_ELEM', copilo)

! - TRANSFORMATION EN CHAM_ELEM_S
    call celces(copilo, 'V', copils)
    call jeveuo(copils//'.CESD', 'L', jcesd)
    call jeveuo(copils//'.CESL', 'L', jcesl)
    call jeveuo(copils//'.CESV', 'L', vr=cesv)
    nbCell = zi(jcesd-1+1)
    npg = zi(jcesd-1+3)
    nbgmax = nbCell*npg

! - ESPACE MEMOIRE POUR LE TABLEAU A0,A1
    call jeexin(a0a1, iret)
    if (iret .eq. 0) then
        call wkvect(a0a1, 'V V R', 4*nbgmax, ja0a1)
        call wkvect(trav, 'V V I', 4*(nbgmax+1), jtrav)
    else
        call jeveuo(a0a1, 'E', ja0a1)
        call jeveuo(trav, 'E', jtrav)
    end if

! - LECTURE DES COMPOSANTES DU CHAM_ELEM_S
    icmp = 0
    do iCell = 1, nbCell
        do kpg = 1, npg
            call cesexi('C', jcesd, jcesl, iCell, kpg, 1, 1, ja0)
            call cesexi('C', jcesd, jcesl, iCell, kpg, 1, 2, ja1)
            call cesexi('C', jcesd, jcesl, iCell, kpg, 1, 3, ja2)
            call cesexi('C', jcesd, jcesl, iCell, kpg, 1, 4, ja3)
            call cesexi('C', jcesd, jcesl, iCell, kpg, 1, 5, ja4)

! ---     LECTURE DU CODE RETOUR
!
            if (ja4 .ne. 0) then
                if (cesv(ja4) .ne. r8vide()) then
! ---         A T ON REMPLI CODE-RETOUR ? OUI -> PAS DE SOLUTION
                    pilcvg = 1
                    goto 999
                end if
            end if
!
! ---     COEFFICIENTS DE LA OU DES DROITES
!
            if (ja0 .ne. 0) then
                if (cesv(ja0) .ne. r8vide()) then
                    zr(ja0a1+icmp) = cesv(ja0)
                    zr(ja0a1+icmp+1) = cesv(ja1)
                    icmp = icmp+2
                    if (cesv(ja2) .ne. r8vide()) then
                        zr(ja0a1+icmp) = cesv(ja2)
                        zr(ja0a1+icmp+1) = cesv(ja3)
                        icmp = icmp+2
                    end if
                end if
            end if
        end do
    end do
!
    npg = icmp/2

! - RESOLUTION DE L'EQUATION DE PILOTAGE P(U(ETA)) = TAU
    call pipere(npg, zr(ja0a1), tau, nbeffe, eta)
!
    if (nbeffe .eq. 0) then
        pilcvg = 1
    end if
!
999 continue
!
    call jedema()
end subroutine
