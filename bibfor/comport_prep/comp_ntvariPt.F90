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
subroutine comp_ntvariPt(comporList, comporInfo, &
                         ntVari, nbVariMaxi, prepExte)
!
    use BehaviourPrepare_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/compGetMecaPart.h"
#include "asterfort/getExternalBehaviourParaPtAdr.h"
#include "asterfort/jeveuo.h"
!
    character(len=16), intent(in) :: comporList(COMPOR_SIZE)
    character(len=19), intent(in) :: comporInfo
    integer(kind=8), intent(out) :: ntVari, nbVariMaxi
    type(BehaviourPrep_Exte), intent(out) :: prepExte
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of constitutive laws (mechanics)
!
! Count total of internal state variables
!
! --------------------------------------------------------------------------------------------------
!
! In  comporList       : list for parameters of constitutive laws
! In  comporInfo       : object for information about internal state variables and behaviour
! Out ntVari           : total number of internal variables (on all <CARTE> COMPOR)
! Out nbVariMaxi       : maximum number of internal variables on all comportments"
! Out prepExte         : pointer to external behaviours parameters
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), pointer :: comporInfoZone(:) => null()
    integer(kind=8) :: nbVari
    character(len=16) :: relaComp, defoComp, kitComp(4), relaMeca
    character(len=16) :: adrsMGIS
!
! --------------------------------------------------------------------------------------------------
!
    ntVari = 0
    nbVariMaxi = 0

! - Create list of zones: for each zone (in CARTE), how many elements
    call jeveuo(comporInfo(1:19)//'.ZONE', 'L', vi=comporInfoZone)

! - Get parameters
    adrsMGIS = comporList(MGIS_ADDR)
    relaComp = comporList(RELA_NAME)
    defoComp = comporList(DEFO)
    kitComp(1) = comporList(KIT1_NAME)
    kitComp(2) = comporList(KIT2_NAME)
    kitComp(3) = comporList(KIT3_NAME)
    kitComp(4) = comporList(KIT4_NAME)

! - Get mechanical part of behaviour
    call compGetMecaPart(relaComp, kitComp, relaMeca)

! - Get parameters for external programs (MFRONT/UMAT)
    call getExternalBehaviourParaPtAdr(relaMeca, defoComp, &
                                       adrsMGIS, prepExte)

! - Get number of internal variables
    read (comporList(NVAR), '(I16)') nbVari
    ntVari = ntVari+nbVari
    nbVariMaxi = max(nbVariMaxi, nbVari)
!
end subroutine
