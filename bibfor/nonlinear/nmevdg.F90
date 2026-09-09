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
subroutine nmevdg(sddisc, hvalIncr, iEvent, i_echec_acti)
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "event_def.h"
#include "asterfort/assert.h"
#include "asterfort/extdch.h"
#include "asterfort/infdbg.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utdidt.h"
!
    integer(kind=8) :: iEvent, i_echec_acti
    character(len=19) :: sddisc, hvalIncr(*)
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME - EVENEMENTS)
!
! GESTION DE L'EVENEMENT DELTA_GRANDEUR
!
! --------------------------------------------------------------------------------------------------
!
! In  sddisc           : datastructure for time discretization TEMPORELLE
! IN  VALE   : INCREMENTS DES VARIABLES
!               OP0070: VARIABLE CHAPEAU
!               OP0033: TABLE
! IN  iEvent : OCCURRENCE DE L'ECHEC
! OUT i_echec_acti : VAUT iEvent SI EVENEMENT DECLENCHE
!                   0 SINON
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: typeExtr = 'MAX_ABS'
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: deb, fin, etat_loca
    integer(kind=8), pointer:: loca(:) => null()
    real(kind=8) :: valeRefe, vale
    character(len=8) :: crit
    character(len=16) :: fieldType, cmpName
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        write (ifm, *) '<MECANONLINE> ... DELTA_GRANDEUR'
    end if

! - INITIALISATIONS
    i_echec_acti = 0

! - PARAMETRES
    call utdidt('L', sddisc, 'ECHE', 'NOM_CHAM', index_=iEvent, valk_=fieldType)
    call utdidt('L', sddisc, 'ECHE', 'NOM_CMP', index_=iEvent, valk_=cmpName)
    call utdidt('L', sddisc, 'ECHE', 'VALE_REF', index_=iEvent, valr_=valeRefe)
    call utdidt('L', sddisc, 'ECHE', 'CRIT_COMP', index_=iEvent, valk_=crit)
    ASSERT(crit .eq. 'GT')

! - Extraction du filtre sur la liste des mailles
    call jeveuo(sddisc//'.ELOC', 'L', vi=loca)
    etat_loca = loca(SIZE_LELOCA*(iEvent-1)+1)

    if (etat_loca .eq. LOCA_VIDE) then
        vale = 0
    else if (etat_loca .eq. LOCA_PARTIEL) then
        deb = loca(SIZE_LELOCA*(iEvent-1)+2)
        fin = loca(SIZE_LELOCA*(iEvent-1)+3)
        call extdch(typeExtr, hvalIncr, fieldType, cmpName, vale, lst_loca=loca(deb:fin))
    else if (etat_loca .eq. LOCA_TOUT) then
        call extdch(typeExtr, hvalIncr, fieldType, cmpName, vale)
    else
        ASSERT(ASTER_FALSE)
    end if

!
    if (vale .gt. valeRefe) then
        i_echec_acti = iEvent
    end if
!
    call jedema()
end subroutine
