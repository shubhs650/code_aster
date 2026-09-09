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
subroutine epstmc(materPara, &
                  poum, time, ndim, &
                  indxVarcStrain, allVarcStrain, &
                  epsiVarc_)
!
    use BehaviourStrain_module
    use BehaviourStrain_type
    use MaterialPara_type
    use MaterialPara_module
    implicit none
!
#include "asterc/r8vide.h"
#include "asterfort/ElasticityMaterial_type.h"
#include "asterfort/lteatt.h"
#include "asterfort/matrot.h"
#include "asterfort/utmess.h"
#include "asterfort/utpslg.h"
#include "jeveux.h"
!
    type(Material_Para), intent(inout) :: materPara
    character(len=*), intent(in) :: poum
    real(kind=8), intent(in) :: time
    integer(kind=8), intent(in) :: ndim
    integer(kind=8), intent(in) :: indxVarcStrain
    type(All_Varc_Strain), intent(inout) :: allVarcStrain
    real(kind=8), optional, intent(out) :: epsiVarc_(6)
!
! --------------------------------------------------------------------------------------------------
!
! Compute inelastic strains from external state variables
!
! For isoparametric elements
!
! --------------------------------------------------------------------------------------------------
!
! IO  materPara        : parameters of material
! In  poum             : '-'  '+' or 'T' (previous, current and both)
! In  time             : given time
! In  ndim             : dimension of space
! In  indxVarcStrain   : index of external state variable
! IO  allVarcStrain    : all external state variables for anelastic strains
! Out epsiVarc         : anelastic strains from all external state variables
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: lTHM
    real(kind=8) :: epsiVarcLoca(6), epsiVarcLocaIn(6), epsiVarcLocaOut(6)
    real(kind=8) :: pgl(3, 3)
!
! --------------------------------------------------------------------------------------------------
!
    lTHM = lteatt('TYPMOD2', 'THM')

! - Current time
    allVarcStrain%time = time
    if (time .eq. r8vide()) then
        allVarcStrain%hasTime = ASTER_FALSE
    else
        allVarcStrain%hasTime = ASTER_TRUE
    end if

! - Detect external state variable
    call strainDetectVarc(poum, lTHM, materPara, &
                          allVarcStrain, indxVarcStrain)

! - Compute non-mechanical strains
    call compVarcStrain(poum, materPara, allVarcStrain)

! - Return values if required
    if (present(epsiVarc_)) then
        epsiVarc_ = 0.d0
        call getVarcStrain(poum, indxVarcStrain, allVarcStrain, epsiVarcLoca)

! ----- Non-isotropic elasticity: rotate strains
        if (materPara%elasID .eq. ELAS_ISOT) then
            epsiVarc_ = epsiVarcLoca
        else
            if (.not. chckLCSDefine(materPara%lcsPara)) then
                call utmess("F", "ALGORITH8_20")
            end if
            call matrot(materPara%lcsPara%lcsAngle, pgl)
            epsiVarcLocaIn(1) = epsiVarcLoca(1)
            epsiVarcLocaIn(2) = epsiVarcLoca(4)
            epsiVarcLocaIn(3) = epsiVarcLoca(2)
            epsiVarcLocaIn(4) = epsiVarcLoca(5)
            epsiVarcLocaIn(5) = epsiVarcLoca(6)
            epsiVarcLocaIn(6) = epsiVarcLoca(3)
            call utpslg(1, 3, pgl, epsiVarcLocaIn, epsiVarcLocaOut)
            epsiVarc_(1) = epsiVarcLocaOut(1)
            epsiVarc_(2) = epsiVarcLocaOut(3)
            epsiVarc_(3) = epsiVarcLocaOut(6)
            epsiVarc_(4) = epsiVarcLocaOut(2)
            epsiVarc_(5) = epsiVarcLocaOut(4)
            epsiVarc_(6) = epsiVarcLocaOut(5)
            if (ndim .eq. 2) then
                epsiVarc_(3) = epsiVarcLoca(3)
            end if
        end if
    end if
!
end subroutine
