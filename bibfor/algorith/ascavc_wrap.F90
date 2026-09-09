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

subroutine ascavc_wrap(model, list_load, numedd, inst, vci, base)
!
    use HHO_type
    use HHO_Dirichlet_module
    use listLoad_module, only: hasDiriElimBC
!
    implicit none
#include "asterfort/as_deallocate.h"
#include "asterfort/ascavc.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"

    character(len=8) :: model
    character(len=19) :: list_load, vci
    character(len=14) :: numedd
    real(kind=8) :: inst
    character(len=1) :: base
!
    type(HHO_Field) :: hhoField
    character(len=8) :: answer
    character(len=24) :: lchar, infcha, fomult

! Wrapper to call AVCAVC from C++

    call dismoi('EXI_HHO', model, 'MODELE', repk=answer)
!
    lchar = list_load//'.LCHA'
    infcha = list_load//'.INFC'
    fomult = list_load//'.FCHA'
!
    if (answer .eq. 'OUI') then
        if (hasDiriElimBC(list_load)) then
!
! --- Prepare fields for Dirichlet loads
!
            hhoField%fieldCineFunc = '&&HHO.CINEFUNC'
            hhoField%fieldCineVale = '&&HHO.CINEVALE'

            call hhoDiriFuncPrepare(model, list_load, hhoField)

            if (hhoField%l_cine_f) then
                call hhoDiriFuncCompute(model, hhoField, inst)
            end if

            call ascavc(lchar, infcha, fomult, numedd, inst, vci, &
                        l_hho_=ASTER_TRUE, hhoField_=hhoField, basez=base)
!
! --- Cleaning
            call detrsd("CARTE", hhoField%fieldCineFunc)
            call detrsd('CHAM_ELEM', hhoField%fieldCineVale)
            if (hhoField%l_cine_f) then
                AS_DEALLOCATE(vi=hhoField%v_info_cine)
            end if
        end if
    else
        call ascavc(lchar, infcha, fomult, numedd, inst, vci, basez=base)
    end if

end subroutine
