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
subroutine nmeceb(sderro, loopName, loopState)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/NonLinear_type.h"
#include "jeveux.h"
!
    character(len=24), intent(in) :: sderro
    character(len=4), intent(in) :: loopName, loopState
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME)
!
! ECRITURE DE L'ETAT DE LA BOUCLE
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
! IN  ETABCL : ETAT DE LA BOUCLE
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

! - Access to datastructure
    eventCONVJv = sderro(1:19)//'.CONV'
    call jeveuo(eventCONVJv, 'E', vi=eventCONV)
!
! --- SELON ETAT
!
    if (loopState .eq. 'CONT') then
        convState = LOOP_STATE_CONTINUE
    else if (loopState .eq. 'CONV') then
        convState = LOOP_STATE_CONVERGE
    else if (loopState .eq. 'EVEN') then
        convState = LOOP_STATE_EVENT
    else if (loopState .eq. 'ERRE') then
        convState = LOOP_STATE_ERROR
    else if (loopState .eq. 'STOP') then
        convState = LOOP_STATE_STOP
    else if (loopState .eq. 'CTCD') then
        convState = LOOP_STATE_CTCD
    else if (loopState .eq. 'ARCH') then
        convState = LOOP_STATE_ARCH
    else
        ASSERT(ASTER_FALSE)
    end if

! - Set state of convergence for loop
    if (loopName .eq. 'RESI') then
        eventCONV(LOOP_RESI) = convState
    else if (loopName .eq. 'NEWT') then
        eventCONV(LOOP_NEWT) = convState
    else if (loopName .eq. 'FIXE') then
        eventCONV(LOOP_FIXE) = convState
    else if (loopName .eq. 'INST') then
        eventCONV(LOOP_INST) = convState
    else if (loopName .eq. 'CALC') then
        eventCONV(LOOP_CALC) = convState
    else
        ASSERT(ASTER_FALSE)
    end if
!
    call jedema()
end subroutine
