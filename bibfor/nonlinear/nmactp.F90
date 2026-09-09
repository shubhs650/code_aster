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
subroutine nmactp(ds_print, sddisc, sderro, ds_contact, &
                  ds_conv, nbiter, numins)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "event_def.h"
#include "asterfort/assert.h"
#include "asterfort/nmacto.h"
#include "asterfort/nmeceb.h"
#include "asterfort/nmevac.h"
#include "asterfort/nmleeb.h"
#include "asterfort/utmess.h"
!
    type(NL_DS_Print), intent(in) :: ds_print
    character(len=24), intent(in) :: sderro
    type(NL_DS_Contact), intent(in) :: ds_contact
    character(len=19), intent(in) :: sddisc
    type(NL_DS_Conv), intent(in) :: ds_conv
    integer(kind=8), intent(in) :: nbiter
    integer(kind=8), intent(in) :: numins
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME)
!
! GESTION DES ACTIONS A LA FIN D'UN PAS DE TEMPS
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_print         : datastructure for printing parameters
! In  sddisc           : datastructure for time discretization
! IN  SDERRO : SD GESTION DES ERREURS
! In  ds_contact       : datastructure for contact management
! In  ds_conv          : datastructure for convergence management
! IN  NBITER : NOMBRE D'ITERATIONS DE NEWTON
! IN  NUMINS : NUMERO D'INSTANT
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: retact, i_echec_acti, actpas, iterat
    character(len=4) :: etinst
!
! --------------------------------------------------------------------------------------------------
!
    retact = 4
    actpas = 3
    iterat = nbiter-1

! - ETAT DE LA BOUCLE EN TEMPS ?
    call nmleeb(sderro, 'INST', etinst)

! - ACTIONS SUITE A UN EVENEMENT
    if (etinst .eq. 'CONV') then
        retact = 0
    else if (etinst .eq. 'EVEN') then
        call nmacto(sddisc, i_echec_acti)
        call nmevac(sddisc, sderro, i_echec_acti, numins, iterat, &
                    retact, ds_print, ds_contact)
    else if (etinst .eq. 'ERRE') then
        retact = 1
    else if (etinst .eq. 'STOP') then
        retact = 4
    else
        ASSERT(.false.)
    end if

! --- TRAITEMENT DE L'ACTION
    if (retact .eq. 0) then
! ----- TOUT EST OK -> ON PASSE A LA SUITE
        actpas = 0

    else if (retact .eq. 1) then
! ----- ON REFAIT LE PAS DE TEMPS
        actpas = 1

    else if (retact .eq. 2) then
! ----- PAS D'ITERATION EN PLUS ICI
        ASSERT(ASTER_FALSE)

    else if (retact .eq. 3) then
! ----- ECHEC DE L'ACTION
        if (.not. ds_conv%l_stop) then
! ------- CONVERGENCE FORCEE -> ON PASSE A LA SUITE
            call utmess('A', 'MECANONLINE2_37')
            actpas = 0
        else
! ------- ARRET DU CALCUL
            actpas = 3
        end if

    else if (retact .eq. 4) then
! ----- ARRET DU CALCUL
        actpas = 3

    else if (retact .eq. 5) then
!
! ----- ARRET DU CALCUL AVEC ARCHIVAGE DERNIER PAS DE TEMPS
!
        actpas = 4
    else
        ASSERT(ASTER_FALSE)
    end if

! - CHANGEMENT DE STATUT DE LA BOUCLE
    if (actpas .eq. 0) then
        call nmeceb(sderro, 'INST', 'CONV')
    else if (actpas .eq. 1) then
        call nmeceb(sderro, 'INST', 'ERRE')
    else if (actpas .eq. 3) then
        call nmeceb(sderro, 'INST', 'STOP')
    else if (actpas .eq. 4) then
        call nmeceb(sderro, 'INST', 'ARCH')
    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
