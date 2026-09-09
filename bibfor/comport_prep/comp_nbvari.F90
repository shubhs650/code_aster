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
subroutine comp_nbvari(relaComp, defoComp, typeCpla, kitComp, &
                       postIter, multComp, reguVisc, postIncr, &
                       solvBehavType, adrsMGIS, &
                       nbVariUMAT, &
                       nbVari, numeLaw, nbVariKit, numeLawKit)
!
    implicit none
!
#include "asterc/mgis_get_sizeof_isvs.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/comp_meca_l.h"
#include "asterfort/comp_nbvari_kit.h"
#include "asterfort/comp_nbvari_std.h"
#include "asterfort/jeveuo.h"
!
    character(len=16), intent(in) :: relaComp, defoComp, typeCpla, kitComp(4)
    character(len=16), intent(in) :: postIter, multComp, reguVisc, postIncr
    integer(kind=8), intent(in) :: solvBehavType
    character(len=16), intent(in) :: adrsMGIS
    integer(kind=8), intent(in) :: nbVariUMAT
    integer(kind=8), intent(out) :: nbVari, numeLaw, nbVariKit(4), numeLawKit(4)
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of constitutive laws (mechanics)
!
! Count the number of internal state variables and index of behaviours
!
! --------------------------------------------------------------------------------------------------
!
! In  relaComp         : behaviour (RELATION keyword)
! In  defoComp         : model of strain (DEFORMATION keyword)
! In  typeCpla         : plane stress method (analytical or De Borst algorithm)
! In  kitComp          : KIT behaviour
! In  postIter         : type of post_treatment at each Newton iteration (POST_ITER keyword)
! In  multComp         : name of map for multi-behaviours (DEFI_COMPOR)
! In  reguVisc         : keyword for viscuous regularization (REGU_VISC keyword)
! In  postIncr         : type of post-treatment at end of time step (POST_INCR keyword)
! In  solvBehavType    : type of type of integration (internal, official, proto, umat)
! In  adrsMGIS         : address (hexadecimal) for the MGIS Behaviour
! In  nbVariUMAT       : number of internal state variables for UMAT
! Out nbVari           : number of internal state variables
! Out numeLaw          : index of subroutine for behaviour
! Out nbVariKit        : number of internal state variables for components in kit
! Out numeLawKit       : index of subroutine for components in kit
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nbVariExte, nbVariFromKit, nbVariCrystal
    aster_logical :: l_cristal, l_kit_meta, l_kit_thm, l_kit_ddi, l_kit_cg, l_kit
    integer(kind=8), pointer :: cpri(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    nbVari = 0
    numeLaw = 0
    nbVariKit = 0
    numeLawKit = 0

! - Detection of specific cases
    call comp_meca_l(relaComp, 'KIT', l_kit)
    call comp_meca_l(relaComp, 'CRISTAL', l_cristal)
    call comp_meca_l(relaComp, 'KIT_META', l_kit_meta)
    call comp_meca_l(relaComp, 'KIT_THM', l_kit_thm)
    call comp_meca_l(relaComp, 'KIT_DDI', l_kit_ddi)
    call comp_meca_l(relaComp, 'KIT_CG', l_kit_cg)

! - Get number of internal state variables for KIT
    nbVariFromKit = 0
    if (l_kit) then
        call comp_nbvari_kit(kitComp, &
                             l_kit_meta, l_kit_thm, l_kit_ddi, l_kit_cg, &
                             nbVariFromKit, nbVariKit, numeLawKit)
    end if

! - Special for CRISTAL
    nbVariCrystal = 0
    if (l_cristal) then
        call jeveuo(multComp(1:8)//'.CPRI', 'L', vi=cpri)
        nbVariCrystal = cpri(3)
        if (defoComp .eq. 'SIMO_MIEHE') then
            nbVariCrystal = nbVariCrystal+3+9
        end if
    end if

! - Get number of internal state variables
    call comp_nbvari_std(relaComp, defoComp, typeCpla, &
                         kitComp, postIter, reguVisc, postIncr, &
                         nbVari, numeLaw)

! - Get number of internal state variables for external behaviours
    if (solvBehavType .eq. SOLV_BEHAV_ASTER) then
        nbVariExte = 0
    elseif (solvBehavType .eq. SOLV_BEHAV_MGIS_OFFI .or. &
            solvBehavType .eq. SOLV_BEHAV_MGIS_PROTO) then
        call mgis_get_sizeof_isvs(adrsMGIS, nbVariExte)
        if (nbVariExte .eq. 0) then
            nbVariExte = 1
        end if
        nbVariKit(4) = nbVariExte
    elseif (solvBehavType .eq. SOLV_BEHAV_UMAT) then
        nbVariExte = nbVariUMAT
        nbVariKit(4) = nbVariExte
    else
        ASSERT(ASTER_FALSE)
    end if

! - Total number of internal state variables
    nbVari = nbVariFromKit+nbVari
    nbVari = nbVariCrystal+nbVari
    nbVari = nbVariExte+nbVari
!
end subroutine
