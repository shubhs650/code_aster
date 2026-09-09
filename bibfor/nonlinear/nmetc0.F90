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

subroutine nmetc0(model, caraElem, compor, ds_inout)
!
    use coorSyst_module, only: setOrieFields
    use NonLin_Datastructure_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/calcul.h"
#include "asterfort/dismoi.h"
#include "asterfort/alchml.h"
!
!
    character(len=24), intent(in) :: model
    character(len=24), intent(in) :: caraElem
    character(len=19), intent(in) :: compor
    type(NL_DS_InOut), intent(in) :: ds_inout
!
! --------------------------------------------------------------------------------------------------
!
! *_NON_LINE - Input/output datastructure
!
! Compute initial fields if necessary
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : name of model
! In  caraElem        : name of datastructure for elementary parameters (CARTE)
! In  compor           : name of <CARTE> COMPOR
! In  ds_inout         : datastructure for input/output management
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbFieldInMax = 100, nbFieldOut = 1
    character(len=8) :: lpain(nbFieldInMax), lpaout(nbFieldOut)
    character(len=19) :: lchin(nbFieldInMax), lchout(nbFieldOut)
!
    integer(kind=8) :: nbFieldIn
    character(len=24) :: field_type, init_name
    character(len=24) :: sief_init, vari_init, strx_init
    integer(kind=8) :: i_field, nb_field, iret
    character(len=19) :: modelLigrel
    aster_logical :: l_sief, l_vari, l_strx, l_acti
!
! --------------------------------------------------------------------------------------------------
!
    nb_field = ds_inout%nb_field
!
! - Create initial fields or not ?
!
    l_sief = ASTER_FALSE
    l_vari = ASTER_FALSE
    l_strx = ASTER_FALSE
    do i_field = 1, nb_field
        field_type = ds_inout%field(i_field)%type
        init_name = ds_inout%field(i_field)%init_name
        l_acti = ds_inout%l_field_acti(i_field)
        if (field_type .eq. 'SIEF_ELGA' .and. l_acti) then
            l_sief = .true.
            sief_init = init_name
        end if
        if (field_type .eq. 'VARI_ELGA' .and. l_acti) then
            l_vari = .true.
            vari_init = init_name
        end if
        if (field_type .eq. 'STRX_ELGA' .and. l_acti) then
            l_strx = .true.
            strx_init = init_name
        end if
    end do
!
! - Initial fields: compute stress and internal variables
!
    if (l_vari .or. l_sief) then
        call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)
        call alchml(modelLigrel, 'TOU_INI_ELGA', 'PSIEF_R', 'V', sief_init, iret, compor)
        call alchml(modelLigrel, 'TOU_INI_ELGA', 'PVARI_R', 'V', vari_init, iret, compor)
    end if

! - Initial fields: special multifibers field
    if (l_strx) then
        nbFieldIn = 0
! ----- Add fields for orientation
        call setOrieFields(nbFieldInMax, lpain, lchin, &
                           nbFieldIn, caraElem)

        lpaout(1) = 'PSTRX_R'
        lchout(1) = strx_init(1:19)
        call calcul('S', 'INI_STRX', modelLigrel, 1, lchin, &
                    lpain, 1, lchout, lpaout, 'V', &
                    'OUI')
    end if
!
end subroutine
