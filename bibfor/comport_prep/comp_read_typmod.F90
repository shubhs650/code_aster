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
subroutine comp_read_typmod(mesh, modelCell, &
                            cellAffeJv, lAllCellAffe, nbCellAffe, &
                            relaComp, relaCompPY, &
                            factorKeyword, iFactorKeyword, &
                            modelMGIS, cplaMGIS)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/BehaviourMGIS_type.h"
#include "asterfort/comp_mfront_modelem.h"
#include "asterfort/dismoi.h"
#include "asterfort/getMFrontPlaneStress.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/utmess.h"
!
    character(len=8), intent(in) :: mesh
    integer(kind=8), pointer :: modelCell(:)
    character(len=24), intent(in) :: cellAffeJv
    aster_logical, intent(in) :: lAllCellAffe
    integer(kind=8), intent(in):: nbCellAffe
    character(len=16), intent(in) :: factorKeyword
    integer(kind=8), intent(in) :: iFactorKeyword
    character(len=16), intent(in) :: relaComp, relaCompPY
    integer(kind=8), intent(out) :: modelMGIS
    character(len=16), intent(out) :: cplaMGIS
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (mechanics)
!
! Find dimension and type of modelisation for MFront
!
! --------------------------------------------------------------------------------------------------
!
! In  mesh             : name of mesh
! In  modelCell        : pointer to list of elements in model
! In  cellAffeJv       : name of object for affected cells
! In  lAllCellAffe     : flag if all cells on mesh have been affected
! In  nbCellAffe       : number of affected cells
! In  relaComp         : behaviour (RELATION keyword)
! In  relaCompPY       : behaviour (RELATION keyword) - For Python
! In  factorKeyword    : factor keyword to read (COMPORTEMENT)
! In  iFactorKeyword   : index of factor keyword
! Out modelMGIS        : type of modelisation MFront
! Out cplaMGIS         : stress plane hypothesis (for Deborst)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nbCell, iCell, cellNume, modelMGISSave
    integer(kind=8) :: elemTypeNume, codret
    aster_logical :: l_mfront_cp
    character(len=16) :: elemTypeName
    integer(kind=8), pointer :: cellAffe(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    modelMGIS = MGIS_MODEL_UNSET
    cplaMGIS = 'VIDE'

! - Access to list of cells
    if (lAllCellAffe) then
        call dismoi('NB_MA_MAILLA', mesh, 'MAILLAGE', repi=nbCell)
    else
        call jeveuo(cellAffeJv, 'L', vi=cellAffe)
        nbCell = nbCellAffe
    end if

! - Get plane stress hypothesis
    call getMFrontPlaneStress(relaComp, relaCompPY, &
                              factorKeyword, iFactorKeyword, &
                              l_mfront_cp)

! - Loop on cells
    modelMGISSave = MGIS_MODEL_UNSET
    do iCell = 1, nbCell
! ----- Get current cell
        if (lAllCellAffe) then
            cellNume = iCell
        else
            cellNume = cellAffe(iCell)
        end if
        elemTypeNume = modelCell(cellNume)

! ----- Select type of modelisation for MFront
        if (elemTypeNume .ne. 0) then
            call jenuno(jexnum('&CATA.TE.NOMTE', elemTypeNume), elemTypeName)
            call comp_mfront_modelem(elemTypeName, l_mfront_cp, &
                                     modelMGIS, cplaMGIS, codret)

            if (modelMGIS .ne. MGIS_MODEL_UNSET) then
                if (modelMGISSave .eq. MGIS_MODEL_UNSET) then
                    modelMGISSave = modelMGIS
                else
                    if ((modelMGISSave .ne. modelMGIS)) then
                        codret = 1
                    end if
                end if
            end if
            if (codret .eq. 1) then
                call utmess('F', 'COMPOR4_13', ni=2, &
                            vali=[modelMGISSave, modelMGIS], &
                            sk="MGISBehaviourFort.h")
            end if
            if (codret .eq. 2) then
                call utmess('F', 'COMPOR4_14', si=modelMGIS, &
                            sk="MGISBehaviourFort.h")
            end if
        end if
    end do

! - Final model for MFront
    modelMGIS = modelMGISSave
!
end subroutine
