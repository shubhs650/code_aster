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
subroutine vechmp(model, materField, materCode, caraElem, &
                  varplu, lXFEM, partps, &
                  nbFieldInMax, lpain, lchin, nbFieldIn)
!
    use coorSyst_module, only: setOrieFields
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/mecact.h"
#include "asterfort/mecoor.h"
#include "asterfort/setStructFields.h"
#include "asterfort/xajcin.h"
#include "jeveux.h"
!
    character(len=8), intent(in) :: model
    character(len=24), intent(in) :: materField, caraElem, materCode
    character(len=19), intent(in) :: varplu
    aster_logical, intent(in) :: lXFEM
    real(kind=8), intent(in) :: partps(3)
    integer(kind=8), intent(in) :: nbFieldInMax
    character(len=8), intent(inout) :: lpain(nbFieldInMax)
    character(len=19), intent(inout) :: lchin(nbFieldInMax)
    integer(kind=8), intent(inout) :: nbFieldIn
!
! --------------------------------------------------------------------------------------------------
!
! CALCUL DES VECTEURS ELEMENTAIRES DES CHARGEMENTS MECANIQUES
! DE NEUMANN
!
! PREPARATION DES CHAMPS D'ENTREE STANDARDS
!
! --------------------------------------------------------------------------------------------------
!
! IN  NOMO   : NOM DU MODELE
! IN  PARTPS : TABLEAU DONNANT T+, DELTAT ET THETA (POUR LE THM)
! IN  CARELE : CARACTERISTIQUES DES POUTRES ET COQUES
! IN  MATE   : MATERIAU CODE
! IN  VARPLU : VARIABLES DE COMMANDE A L'INSTANT T+
! IN  LXFEM  : .TRUE. SI XFEM
! IN  nbFieldInMax   : NOMBRE MAXI DE CHAMPS D'ENTREE
! OUT LPAIN  : LISTE DES PARAMETRES IN
! OUT LCHIN  : LISTE DES CHAMPS IN
! OUT LASTIN : NOMBRE EFFECTIF DE CHAMPS IN
!
! --------------------------------------------------------------------------------------------------
!
    character(len=19), parameter :: chtime = '&&VECHMP.CH_INST_R'
    character(len=8), parameter :: cmpName(3) = (/'INST  ', 'DELTAT', 'THETA '/)
    character(len=19) :: modelLigrel
    character(len=19) :: chgeom
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    lpain = " "
    lchin = " "
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)

! - CHAMP DE GEOMETRIE
    call mecoor(modelLigrel, chgeom)

! - CREATION DE LA CARTE DES INSTANTS
    call mecact('V', chtime, 'LIGREL', modelLigrel, 'INST_R', &
                ncmp=3, lnomcmp=cmpName, vr=partps)

! - Set input fields
    lpain(1) = 'PGEOMER'
    lchin(1) = chgeom
    lpain(2) = 'PINSTR'
    lchin(2) = chtime
    lpain(3) = 'PMATERC'
    lchin(3) = materCode(1:19)
    lpain(4) = 'PVARCPR'
    lchin(4) = varplu
    lpain(5) = 'PCOMPOR'
    lchin(5) = materField(1:8)//'.COMPOR'
    nbFieldIn = 5

! - Add fields for structural elements
    call setStructFields(caraElem, nbFieldInMax, lchin, lpain, nbFieldIn)

! - Add fields for XFEM
    if (lXFEM) then
        call xajcin(model, 'CHAR_MECA_NEUM', nbFieldInMax, lchin, lpain, nbFieldIn)
    end if

! - Add fields for orientation
    call setOrieFields(nbFieldInMax, lpain, lchin, &
                       nbFieldIn, caraElem)
!
    call jedema()
end subroutine
