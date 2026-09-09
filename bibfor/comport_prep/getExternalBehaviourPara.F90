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
subroutine getExternalBehaviourPara(mesh, modelCell, &
                                    cellAffeJv, lAllCellAffe, nbCellAffe, &
                                    relaComp, relaCompPY, relaMeca, defoComp, &
                                    factorKeyword, iFactorKeyword, &
                                    prepExte)
!
    use BehaviourPrepare_type
    implicit none
!
#include "asterc/mgis_debug.h"
#include "asterc/mgis_load_library.h"
#include "asterc/umat_get_function.h"
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/BehaviourMGIS_type.h"
#include "asterfort/comp_meca_l.h"
#include "asterfort/comp_read_exte.h"
#include "asterfort/comp_read_mfront.h"
#include "asterfort/comp_read_typmod.h"
#include "asterfort/getExternalStrainModel.h"
!
    character(len=8), intent(in) :: mesh
    integer(kind=8), pointer :: modelCell(:)
    character(len=24), intent(in) :: cellAffeJv
    aster_logical, intent(in) :: lAllCellAffe
    integer(kind=8), intent(in):: nbCellAffe
    character(len=16), intent(in) :: relaComp, relaCompPY, relaMeca, defoComp
    character(len=16), intent(in) :: factorKeyword
    integer(kind=8), intent(in) :: iFactorKeyword
    type(BehaviourPrep_Exte), intent(inout) :: prepExte
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of behaviour (mechanics)
!
! Get parameters for external programs (MFRONT/UMAT)
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
! In  relaMeca         : mechanical part of behaviour
! In  defoComp         : model of strain (DEFORMATION keyword)
! In  factorKeyword    : factor keyword to read (COMPORTEMENT)
! In  iFactorKeyword   : index of factor keyword
! IO  prepExte         : external behaviours parameters
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_mfront_proto, l_mfront_offi, l_umat
    character(len=255) :: librNameUMAT, subrNameUMAT
    integer(kind=8) :: modelMGIS, nbVariUMAT
    character(len=16) :: adrsMGIS, cplaMGIS
    integer(kind=8) :: adrsUMAT, solvBehavType, strainMGIS
!
! --------------------------------------------------------------------------------------------------
!
    l_umat = ASTER_FALSE
    l_mfront_proto = ASTER_FALSE
    l_mfront_offi = ASTER_FALSE
    cplaMGIS = "VIDE"

! - Detect external integrator
    call comp_meca_l(relaMeca, 'UMAT', l_umat)
    call comp_meca_l(relaMeca, 'MFRONT_OFFI', l_mfront_offi)
    call comp_meca_l(relaMeca, 'MFRONT_PROTO', l_mfront_proto)
    solvBehavType = SOLV_BEHAV_ASTER
    if (l_mfront_offi) then
        solvBehavType = SOLV_BEHAV_MGIS_OFFI
    elseif (l_mfront_proto) then
        solvBehavType = SOLV_BEHAV_MGIS_PROTO
    elseif (l_umat) then
        solvBehavType = SOLV_BEHAV_UMAT
    end if

! - Get parameters for UMAT
    adrsUMAT = 0
    librNameUMAT = ' '
    subrNameUMAT = ' '
    nbVariUMAT = 0
    if (solvBehavType == SOLV_BEHAV_UMAT) then
        if (iFactorKeyword .ne. 0) then
            call comp_read_exte(factorKeyword, iFactorKeyword, &
                                librNameUMAT, subrNameUMAT, nbVariUMAT)
            call umat_get_function(librNameUMAT, subrNameUMAT, adrsUMAT)
        end if
    end if

! - Get parameters for MFRONT
    adrsMGIS = " "
    modelMGIS = MGIS_MODEL_UNSET
    strainMGIS = MGIS_STRAIN_UNSET
    if (solvBehavType == SOLV_BEHAV_MGIS_OFFI .or. &
        solvBehavType == SOLV_BEHAV_MGIS_PROTO) then
        call comp_read_mfront(factorKeyword, iFactorKeyword, adrsMGIS)

        if (adrsMGIS .ne. " ") then
! --------- Finite element support
            call comp_read_typmod(mesh, modelCell, &
                                  cellAffeJv, lAllCellAffe, nbCellAffe, &
                                  relaComp, relaCompPY, &
                                  factorKeyword, iFactorKeyword, &
                                  modelMGIS, cplaMGIS)

! --------- Get model of strains and load library
            call getExternalStrainModel(defoComp, strainMGIS)

! --------- Load library
            call mgis_load_library(adrsMGIS, modelMGIS, strainMGIS)
            ! call mgis_debug(adrsMGIS, "Loaded behaviour:")

        end if
    end if

! - Save
    prepExte%solvBehavType = solvBehavType
    prepExte%librNameUMAT = librNameUMAT
    prepExte%subrNameUMAT = subrNameUMAT
    prepExte%adrsUMAT = adrsUMAT
    prepExte%nbVariUMAT = nbVariUMAT
    prepExte%adrsMGIS = adrsMGIS
    prepExte%modelMGIS = modelMGIS
    prepExte%strainMGIS = strainMGIS
    prepExte%cplaMGIS = cplaMGIS
!
end subroutine
