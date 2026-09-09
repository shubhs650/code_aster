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
!
subroutine addGrpMa(mesh, group_ma, listCells, nbCells, l_added_grpma)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/existGrpMa.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jecreo.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
!
    character(len=8), intent(in)  :: mesh
    character(len=24), intent(in) :: group_ma
    integer(kind=8), intent(in)           :: listCells(*)
    integer(kind=8), intent(in)           :: nbCells
    aster_logical, intent(out), optional :: l_added_grpma
!
!---------------------------------------------------------------------------------------------------
!   But :
!     Add a group of Cells in a group_ma of the mesh
!
!   IN:
!     mesh      : name of the mesh
!     group_ma  : name of the group of Cells to test
!     listCells : list of Cells of group_ma
!     nbCells   : number of Cells
!
!  OUT:
!     l_added_grpma : the group_ma has been added to the mesh ?
!
!---------------------------------------------------------------------------------------------------
    character(len=24) :: grmama, grmamap, nomgrp
    integer(kind=8) :: nbGrp, iaux, iret
    aster_logical :: l_parallel_mesh, l_exi_in_grp, l_exi_in_grp_p, l_added
    integer(kind=8), pointer :: v_cells(:) => null()
    character(len=24), pointer :: v_grpp(:) => null()
!-----------------------------------------------------------------------
!
    grmama = mesh//'.GROUPEMA'
    grmamap = mesh//'.PAR_GRPMAI'
!
    call jemarq()
!
    l_parallel_mesh = isParallelMesh(mesh)
!
    call existGrpMa(mesh, group_ma, l_exi_in_grp, l_exi_in_grp_p)
!
! --- The group already exists - do nothing
!
    if (l_exi_in_grp) then
        l_added = ASTER_FALSE
        call utmess('A', 'SOUSTRUC_88', sk=group_ma)
    elseif (nbCells <= 0) then
        l_added = ASTER_FALSE
        call utmess('A', 'SOUSTRUC_36', sk=group_ma)
    else
        l_added = ASTER_TRUE
        ASSERT(.not. l_exi_in_grp)
!
! --- Add group_ma
!
        call jecroc(jexnom(grmama, group_ma))
        call jeecra(jexnom(grmama, group_ma), 'LONMAX', max(nbCells, 1))
        call jeecra(jexnom(grmama, group_ma), 'LONUTI', nbCells)
        call jeveuo(jexnom(grmama, group_ma), 'E', vi=v_cells)
!
! --- Add Cells
!
        v_cells(1:nbCells) = listCells(1:nbCells)
!
! --- For a ParallelMesh, add new group_ma to global group_ma
!
        if (l_parallel_mesh) then
            if (.not. l_exi_in_grp_p) then
                call jeexin(grmamap, iret)
                if (iret .ne. 0) then
                    call jelira(grmamap, 'NOMMAX', nbGrp)
                    nbGrp = abs(nbGrp)
                else
                    nbGrp = 0
                end if
                call wkvect('&&TMP', 'V V K24', nbGrp+1, vk24=v_grpp)
                if (iret .ne. 0) then
                    do iaux = 1, nbGrp
                        call jenuno(jexnum(grmamap, iaux), nomgrp)
                        v_grpp(iaux) = nomgrp
                    end do
                end if
                v_grpp(nbGrp+1) = group_ma
                call jedetr(grmamap)
                call jecreo(grmamap, 'G N K24')
                call jeecra(grmamap, 'NOMMAX', nbGrp+1)
                do iaux = 1, nbGrp+1
                    call jecroc(jexnom(grmamap, v_grpp(iaux)))
                end do
                call jedetr('&&TMP')
!
! --- Becarefull, there are not yet shared by all meshes
!
            end if
        end if
    end if
!
    if (present(l_added_grpma)) then
        l_added_grpma = l_added
    end if
!
    call jedema()
end subroutine
