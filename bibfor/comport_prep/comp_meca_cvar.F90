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
subroutine comp_meca_cvar(prepMapCompor)
!
    use BehaviourPrepare_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/comp_nbvari.h"
!
    type(BehaviourPrep_MapCompor), intent(inout) :: prepMapCompor
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of constitutive laws (mechanics)
!
! Count all internal state variables
!
! --------------------------------------------------------------------------------------------------
!
! IO  prepMapCompor    : datastructure to construct COMPOR map
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iFactorKeyword, nbFactorKeyword
    character(len=16) :: relaComp, defoComp, multComp, kitComp(4), typeCpla
    character(len=16) :: reguVisc, postIter, postIncr
    integer(kind=8) :: numeLawKit(4), nbVari, nbVariKit(4), numeLaw, nbVariUMAT
    integer(kind=8) :: solvBehavType
    character(len=16) :: adrsMGIS
!
! --------------------------------------------------------------------------------------------------
!
    nbFactorKeyword = prepMapCompor%nb_comp
    do iFactorKeyword = 1, nbFactorKeyword
        nbVari = 0
        numeLaw = 0
        nbVariKit = 0
        numeLawKit = 0

! ----- Get parameters
        relaComp = prepMapCompor%prepPara(iFactorKeyword)%rela_comp
        defoComp = prepMapCompor%prepPara(iFactorKeyword)%defo_comp
        typeCpla = prepMapCompor%prepPara(iFactorKeyword)%type_cpla
        kitComp = prepMapCompor%prepPara(iFactorKeyword)%kit_comp
        multComp = prepMapCompor%prepPara(iFactorKeyword)%mult_comp
        postIter = prepMapCompor%prepPara(iFactorKeyword)%post_iter
        reguVisc = prepMapCompor%prepPara(iFactorKeyword)%regu_visc
        postIncr = prepMapCompor%prepPara(iFactorKeyword)%post_incr
        nbVariUMAT = prepMapCompor%prepExte(iFactorKeyword)%nbVariUMAT
        adrsMGIS = prepMapCompor%prepExte(iFactorKeyword)%adrsMGIS
        solvBehavType = prepMapCompor%prepExte(iFactorKeyword)%solvBehavType

! ----- Count the number of internal state variables and index of behaviours
        call comp_nbvari(relaComp, defoComp, typeCpla, kitComp, &
                         postIter, multComp, reguVisc, postIncr, &
                         solvBehavType, adrsMGIS, &
                         nbVariUMAT, &
                         nbVari, numeLaw, nbVariKit, numeLawKit)

! ----- Save informations
        prepMapCompor%prepPara(iFactorKeyword)%nbVari = nbVari
        prepMapCompor%prepPara(iFactorKeyword)%nbVariKit = nbVariKit
        prepMapCompor%prepPara(iFactorKeyword)%numeLaw = numeLaw
        prepMapCompor%prepPara(iFactorKeyword)%numeLawKit = numeLawKit

        if (prepMapCompor%lDebug) then
            WRITE (6, *) "- Occurrence : ", iFactorKeyword
            WRITE (6, *) "--- nbVari : ", prepMapCompor%prepPara(iFactorKeyword)%nbVari
            WRITE (6, *) "--- nbVariKit : ", prepMapCompor%prepPara(iFactorKeyword)%nbVariKit
            WRITE (6, *) "--- numeLaw : ", prepMapCompor%prepPara(iFactorKeyword)%numeLaw
            WRITE (6, *) "--- numeLawKit : ", prepMapCompor%prepPara(iFactorKeyword)%numeLawKit
        end if

    end do
!
end subroutine
