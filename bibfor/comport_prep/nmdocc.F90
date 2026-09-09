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
subroutine nmdocc(model, chmate, lInitialState, compor, base, l_verbose)
!
    use BehaviourPrepare_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/comp_info.h"
#include "asterfort/comp_init.h"
#include "asterfort/comp_meca_chck.h"
#include "asterfort/comp_meca_cvar.h"
#include "asterfort/comp_meca_elas.h"
#include "asterfort/comp_meca_full.h"
#include "asterfort/comp_meca_info.h"
#include "asterfort/comp_meca_read.h"
#include "asterfort/comp_meca_save.h"
#include "asterfort/dismoi.h"
#include "asterfort/infniv.h"
#include "asterfort/utmess.h"
!
    character(len=8), intent(in) :: model, chmate
    aster_logical, intent(in) :: lInitialState
    character(len=19), intent(in) :: compor
    character(len=1), intent(in) :: base
    aster_logical, intent(in), optional :: l_verbose
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of behaviours (mechanics)
!
! Get parameters from COMPORTEMENT keyword and prepare COMPOR <CARTE>
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : model
! In  chmate           : material field
! In  lInitialState      : .true. if initial state is defined
! In  compor           : map for parameters of constitutive laws
! In  base             : permanent or temporary allocation
! In  l_verbose        : .true. to enable verbose mode
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    aster_logical :: verbose
    character(len=8) :: mesh
    character(len=19), parameter :: fullElemField = '&&NMDOCC.FULL_ELEM'
    type(BehaviourPrep_MapCompor) :: prepMapCompor
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE12_4')
    end if

! - Initialisations
    verbose = ASTER_FALSE
    if (present(l_verbose)) then
        verbose = l_verbose
    end if
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)

! - Create datastructure to prepare comportement
    call comp_meca_info(prepMapCompor)

! - Create COMPOR <CARTE>
    call comp_init(mesh, compor, base)

! - Set ELASTIQUE COMPOR
    call comp_meca_elas(compor, lInitialState)

! - Read informations from command file
    call comp_meca_read(lInitialState, prepMapCompor, model)

! - Create <CARTE> of FULL_MECA option for checking
    call comp_meca_full(model, compor, fullElemField)

! - Check informations in COMPOR <CARTE>
    call comp_meca_chck(model, mesh, chmate, fullElemField, lInitialState, prepMapCompor)

! - Count internal variables
    call comp_meca_cvar(prepMapCompor)

! - Save informations in COMPOR <CARTE>
    call comp_meca_save(model, mesh, chmate, compor, prepMapCompor)

! - Verbose mode
    if (verbose) then
        call comp_info(model, compor)
    end if

! - Clean
    deallocate (prepMapCompor%prepPara)
    deallocate (prepMapCompor%prepExte)
!
end subroutine
