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
subroutine diinit(meshZ, modelZ, ds_inout, materField, materCode, caraElem, &
                  listFuncActi, sddyna, ds_conv, ds_algopara, solver, &
                  ds_contact, sddisc)
!
    use NonLin_Datastructure_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/crsvsi.h"
#include "asterfort/getvid.h"
#include "asterfort/isfonc.h"
#include "asterfort/ndxcfl.h"
#include "asterfort/ndynlo.h"
#include "asterfort/nmcrar.h"
#include "asterfort/nmcrli.h"
#include "asterfort/nmcrsu.h"
#include "asterfort/utdidt.h"
#include "asterfort/utmess.h"
!
    character(len=*), intent(in) :: meshZ, modelZ
    character(len=19), intent(in) :: sddisc, sddyna
    character(len=24), intent(in) :: caraElem
    character(len=24), intent(in) :: materField, materCode
    type(NL_DS_Conv), intent(in) :: ds_conv
    type(NL_DS_AlgoPara), intent(in) :: ds_algopara
    type(NL_DS_InOut), intent(in) :: ds_inout
    character(len=19), intent(in) :: solver
    type(NL_DS_Contact), intent(in) :: ds_contact
    integer(kind=8), intent(in) :: listFuncActi(*)
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Datastructures
!
! Time discretization and storing datastructures
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  model            : name of model
! In  materField       : name of material characteristics (field)
! In  caraElem         : name of elementary characteristics (field)
! In  sddyna           : name of dynamic parameters
! In  ds_conv          : datastructure for convergence management
! In  ds_algo          : datastructure for algorithm management
! In  listFuncActi     : active functionnalities vector (see nmfonc)
! In  ds_contact       : datastructure for contact management
! In  solver           : name of solver parameters
! In  sddisc           : datastructure for time discretization
! In  ds_inout         : datastructure for input/output management
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: factorKeyword = "INCREMENT"
    aster_logical :: l_expl, l_implex, lCutStep
    character(len=19) :: listInst, sdarch
    character(len=8) :: model, mesh, result
    character(len=16) :: answer
!
! --------------------------------------------------------------------------------------------------
!
    model = modelZ
    mesh = meshZ

! - Get object for time list from user
    call getvid(factorKeyword, 'LIST_INST', iocc=1, scal=listInst)

! - Get parameters
    result = ds_inout%result

! - Active functionnalities
    l_expl = ndynlo(sddyna, 'EXPLICITE')
    l_implex = isfonc(listFuncActi, 'IMPLEX')

! - Create time discretization datastructure
    call nmcrli(listInst, sddisc)

! - Courant condition
    if (l_expl) then
        call ndxcfl(materField, materCode, caraElem, sddyna, sddisc)
    end if

! - Create storing objects
    sdarch = sddisc(1:14)//'.ARCH'
    call nmcrar(result, sddisc, sdarch, listFuncActi)

! - Automatic management of time stepping
    call nmcrsu(sddisc, listInst, ds_conv, &
                l_implex, ds_contact)

! - Solver management for time step cut
    call utdidt('L', sddisc, 'LIST', 'EXIS_DECOUPE', valk_=answer)
    lCutStep = answer .eq. 'OUI'
    if (lCutStep) then
        call crsvsi(solver)
    end if

! - Restriction of cut step
    if (lCutStep) then
        if (ds_algopara%matrix_pred .eq. 'DEPL_CALCULE') then
            call utmess('F', 'SUBDIVISE_99')
        end if
    end if
!
end subroutine
