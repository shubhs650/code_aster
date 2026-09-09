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
subroutine nonlinLoadDirichletCompute(list_load, model, nume_dof, &
                                      ds_measure, matr_asse, disp, &
                                      hval_veelem, hval_veasse)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assvec.h"
#include "asterfort/conlag.h"
#include "asterfort/dismoi.h"
#include "asterfort/infdbg.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nmchex.h"
#include "asterfort/nmdebg.h"
#include "asterfort/nmtime.h"
#include "asterfort/utmess.h"
#include "asterfort/vebume.h"
#include "asterfort/vector_update_ghost_values.h"
#include "asterfort/vector_zero_late_values.h"
!
    character(len=19), intent(in) :: list_load
    character(len=24), intent(in) :: model, nume_dof
    type(NL_DS_Measure), intent(inout) :: ds_measure
    character(len=19), intent(in) :: matr_asse
    character(len=19), intent(in) :: disp
    character(len=19), intent(in) :: hval_veelem(*), hval_veasse(*)
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Algorithm
!
! Compute values of Dirichlet conditions
!
! --------------------------------------------------------------------------------------------------
!
! In  list_load        : name of datastructure for list of loads
! In  model            : name of model
! In  nume_dof         : name of numbering object (NUME_DDL)
! IO  ds_measure       : datastructure for measure and statistics management
! In  matr_asse        : matrix
! In  disp             : displacements
! In  hval_veelem      : hat-variable for elementary vectors
! In  hval_veasse      : hat-variable for vectors (node fields)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    character(len=19) :: vebudi, cnbudi, nume_equa
    character(len=8) :: mesh
    real(kind=8) :: alpha
    real(kind=8), pointer :: vale(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE11_12')
    end if
!
! - Get hat variables
!
    call nmchex(hval_veelem, 'VEELEM', 'CNBUDI', vebudi)
    call nmchex(hval_veasse, 'VEASSE', 'CNBUDI', cnbudi)
!
! - Launch timer
!
    call nmtime(ds_measure, 'Init', '2nd_Member')
    call nmtime(ds_measure, 'Launch', '2nd_Member')
!
! - Compute
!
    call conlag(matr_asse, alpha)
    call vebume(model, disp, list_load, vebudi, alpha, 'V')
    call assvec('V', cnbudi, 1, vebudi, [1.d0], nume_dof)
!
! - Pour les lagrange distribués, il faut assembler le vecteur pour avoir les contrib
!   des autres procs (ceci arrive quand un Lagrange est relié à plusieurs sous-domaines)
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
    if (isParallelMesh(mesh)) then
        call jeveuo(cnbudi//'.VALE', 'E', vr=vale)
        call dismoi('NUME_EQUA', cnbudi, 'CHAM_NO', repk=nume_equa)
        call vector_update_ghost_values(vale, nume_equa, "BIDIR")
        call vector_zero_late_values(vale, nume_equa)
    end if
!
! - Stop timer
!
    call nmtime(ds_measure, 'Stop', '2nd_Member')
!
! - Debug
!
    if (niv .ge. 2) then
        call nmdebg('VECT', cnbudi, 6)
    end if
!
end subroutine
