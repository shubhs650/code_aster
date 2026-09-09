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
subroutine nmleeb(sderro, loopName, loopState)
!
    implicit none
!
#include "asterfort/NonLinear_type.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
!
    character(len=24), intent(in) :: sderro
    character(len=4), intent(in) :: loopName
    character(len=4), intent(out) :: loopState
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME)
!
! ETAT DE LA CONVERGENCE
!
! --------------------------------------------------------------------------------------------------
!
! IN  SDERRO : SD GESTION DES ERREURS
! IN  NOMBCL : NOM DE LA BOUCLE
!               'RESI' - BOUCLE SUR LES RESIDUS D'EQUILIBRE
!               'NEWT' - BOUCLE DE NEWTON
!               'FIXE' - BOUCLE DE POINT FIXE
!               'INST' - BOUCLE SUR LES PAS DE TEMPS
!               'CALC' - CALCUL
! OUT ETABCL : ETAT DE LA BOUCLE
!               'CONT' - ON CONTINUE LA BOUCLE
!               'CTCD' - ON CONTINUE LA BOUCLE APRES LA PREDICTION
!               'CONV' - ON STOPPE LA BOUCLE : CONVERGEE
!               'EVEN' - EVENEMENT PENDANT LA BOUCLE
!               'ERRE' - ON STOPPE LA BOUCLE : ERREUR TRAITEE
!               'STOP' - ON STOPPE LA BOUCLE : ERREUR NON TRAITEE
!               'ARCH' - ON ARCHIVE LE DERNIER PAS DE TEMPS ET ON ARRETE LE CALCUL
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24) :: eventCONVJv
    integer(kind=8), pointer :: eventCONV(:) => null()
    integer(kind=8) :: convState
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - INITIALISATIONS
    convState = LOOP_STATE_CONTINUE
    loopState = ' '

! - Access to datastructure
    eventCONVJv = sderro(1:19)//'.CONV'
    call jeveuo(eventCONVJv, 'L', vi=eventCONV)

! - Get state of convergence for loop
    if (loopName .eq. 'RESI') then
        convState = eventCONV(LOOP_RESI)
    else if (loopName .eq. 'NEWT') then
        convState = eventCONV(LOOP_NEWT)
    else if (loopName .eq. 'FIXE') then
        convState = eventCONV(LOOP_FIXE)
    else if (loopName .eq. 'INST') then
        convState = eventCONV(LOOP_INST)
    else if (loopName .eq. 'CALC') then
        convState = eventCONV(LOOP_CALC)
    else
        ASSERT(ASTER_FALSE)
    end if

! --- SELON ETAT
    if (convState .eq. LOOP_STATE_CONTINUE) then
        loopState = 'CONT'
    else if (convState .eq. LOOP_STATE_CONVERGE) then
        loopState = 'CONV'
    else if (convState .eq. LOOP_STATE_EVENT) then
        loopState = 'EVEN'
    else if (convState .eq. LOOP_STATE_ERROR) then
        loopState = 'ERRE'
    else if (convState .eq. LOOP_STATE_STOP) then
        loopState = 'STOP'
    else if (convState .eq. LOOP_STATE_CTCD) then
        loopState = 'CTCD'
    else if (convState .eq. LOOP_STATE_ARCH) then
        loopState = 'ARCH'
    else
        ASSERT(ASTER_FALSE)
    end if
!
    call jedema()
end subroutine
