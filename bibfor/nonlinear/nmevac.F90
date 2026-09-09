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
subroutine nmevac(sddisc, sderro, i_fail_acti, nume_inst, iterat, &
                  retact, ds_print_, ds_contact_)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/getFailAction.h"
#include "asterfort/getFailEvent.h"
#include "asterfort/nmdeco.h"
#include "asterfort/nmecev.h"
#include "asterfort/nmeraz.h"
#include "asterfort/nmerge.h"
#include "asterfort/utmess.h"
#include "event_def.h"
!
    character(len=19), intent(in) :: sddisc
    character(len=24), intent(in) :: sderro
    integer(kind=8), intent(in) :: i_fail_acti
    integer(kind=8), intent(in) :: nume_inst
    integer(kind=8), intent(in) :: iterat
    integer(kind=8), intent(out) :: retact
    type(NL_DS_Print), optional, intent(in) :: ds_print_
    type(NL_DS_Contact), optional, intent(in) :: ds_contact_
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME)
!
! ACTIONS SUITE A UN EVENEMENT
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_print         : datastructure for printing parameters
! In  sddisc           : datastructure for time discretization
! IN  SDERRO : SD ERREUR
! In  ds_contact       : datastructure for contact management
! IN  IEVDAC : INDICE DE L'EVENEMENT ACTIF
! IN  NUMINS : NUMERO D'INSTANT
! IN  ITERAT : NUMERO D'ITERATION DE NEWTON
! OUT RETACT : CODE RETOUR
!     0 - ON NE FAIT RIEN
!     1 - ON REFAIT LE PAS DE TEMPS
!     2 - ON CONTINUE LA BOUCLE DE NEWTON (ITERATIONS EN PLUS)
!     3 - L'ACTION A ECHOUE
!     5 - ON ARRETE LE CALCUL ET ON ARCHIVE LE DERNIER PAS DE TEMPS
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: retdec, failType, actionType
    aster_logical :: trydec, litmax
!
! --------------------------------------------------------------------------------------------------
!
    retact = 3
    actionType = FAIL_ACT_STOP
    trydec = ASTER_FALSE
    ASSERT(i_fail_acti .ne. 0)

! --- RECUPERATION ERREURS PARTICULIERES
    litmax = ASTER_FALSE
    if (sderro .ne. ' ') then
        call nmerge(sderro, 'ITER_MAXI', litmax)
    end if

! - Get event and action
    call getFailEvent(sddisc, i_fail_acti, failType)
    call getFailAction(sddisc, i_fail_acti, actionType)

! - Action
    if (actionType .eq. FAIL_ACT_STOP) then
        call utmess('I', 'MECANONLINE10_30')
        retact = 3
        trydec = ASTER_FALSE
    else if (actionType .eq. FAIL_ACT_CUT) then
        trydec = ASTER_TRUE
    else if (actionType .eq. FAIL_ACT_CONTINUE) then
        retact = 0
    else if (actionType .eq. FAIL_ACT_ARCHI) then
        retact = 5
    else
        ASSERT(ASTER_FALSE)
    end if

! - For step cutting
    if (trydec) then
        call utmess('I', 'MECANONLINE10_33')
        call nmdeco(sddisc, nume_inst, iterat, i_fail_acti, retdec)
        if (retdec .eq. 0) then
            retact = 3
        else if (retdec .eq. 1) then
            retact = 1
        else if (retdec .eq. 2) then
            retact = 0
        else
            ASSERT(ASTER_FALSE)
        end if
    end if

! - ECHEC DE L'ACTION -> EVENEMENT ERREUR FATALE
    if (retact .eq. 3) then
        call nmecev(sderro, 'E', failType, actionType)
    end if
!
! --- ON DESACTIVE LES EVENEMENTS
!
    call nmeraz(sderro, 'EVEN')
!
end subroutine
