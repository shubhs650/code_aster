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
subroutine mtdorc(factorKeyword, model, comporMeta, lEmptyMap)
!
    use Metallurgy_type
    implicit none
!
#include "asterfort/Behaviour_type.h"
#include "asterfort/comp_init.h"
#include "asterfort/comp_meta_clean.h"
#include "asterfort/comp_meta_info.h"
#include "asterfort/comp_meta_prnt.h"
#include "asterfort/comp_meta_pvar.h"
#include "asterfort/comp_meta_read.h"
#include "asterfort/comp_meta_save.h"
#include "asterfort/dismoi.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nocart.h"
!
    character(len=16), intent(in) :: factorKeyword
    character(len=8), intent(in) :: model
    character(len=24), intent(in) :: comporMeta
    aster_logical, intent(out) :: lEmptyMap
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (metallurgy)
!
! Construct map for behaviour in metallurgy
!
! --------------------------------------------------------------------------------------------------
!
! In  factorKeyword    : factor keyword to read
! In  model            : name of model
! In  comporMeta       : name of map for behaviour in metallurgy
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nbCmp
    character(len=8) :: mesh
    character(len=19), parameter :: comporMetaInfo = '&&MTDORC.INFO'
    type(META_PrepBehaviour) :: metaPrepBehaviour
    character(len=16), pointer :: comporValv(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    lEmptyMap = ASTER_TRUE
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)

! - Create datastructure to prepare comportement
    call comp_meta_info(factorKeyword, metaPrepBehaviour)
    lEmptyMap = metaPrepBehaviour%nbFactorKeyword .eq. 0

! - Create COMPOR <CARTE>
    call comp_init(mesh, comporMeta, 'V', nbCmp)

! - Access to map
    call jeveuo(comporMeta(1:19)//'.VALV', 'E', vk16=comporValv)

! - Init <CARTE>
    comporValv(1) = "VIDE"
    comporValv(2) = '0'
    comporValv(3) = "VIDE"
    comporValv(4) = '0'

! - Create <CARTE>
    call nocart(comporMeta, 1, nbCmp)

! - Read informations from command file
    call comp_meta_read(metaPrepBehaviour)

! - Save informations in COMPOR <CARTE>
    call comp_meta_save(mesh, comporMeta, nbCmp, metaPrepBehaviour)

! - Prepare informations about internal variables
    call comp_meta_pvar(model, comporMeta, comporMetaInfo)

! - Print informations about internal variables:
    call comp_meta_prnt(metaPrepBehaviour%hasTemper, comporMetaInfo)

! - Clean
    deallocate (metaPrepBehaviour%paraBehaviour)
    call comp_meta_clean(comporMetaInfo)
!
end subroutine
