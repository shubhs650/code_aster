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
subroutine epthmc(materPara, time, &
                  nbEpsi, npg, ndim, &
                  indxVarcStrain, epsiVarc)
!
    use BehaviourStrain_module
    use BehaviourStrain_type
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/epstmc.h"
!
    type(Material_Para), intent(inout) :: materPara
    real(kind=8), intent(in) :: time
    integer(kind=8), intent(in) :: nbEpsi, npg, ndim
    integer(kind=8), intent(in) :: indxVarcStrain
    real(kind=8), intent(out) :: epsiVarc(nbEpsi*npg)
!
! --------------------------------------------------------------------------------------------------
!
! Compute anelastic strains from external state variables
!
! --------------------------------------------------------------------------------------------------
!
! IO  materPara        : parameters of material
! In  time             : given time
! In  ndim             : dimension of space
! In  nbEpsi           : number of strain tensor components
! In  npg              : number of Gauss points
! In  indxVarcStrain   : index of external state variable
! Out epsiVarc         : anelastic strains from all external state variables
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: ksp = 1
    integer(kind=8) :: kpg, iEpsi
    real(kind=8) :: epsiVarcKpg(6)
    type(All_Varc_Strain) :: allVarcStrain
!
! --------------------------------------------------------------------------------------------------
!
    epsiVarc = 0.d0
    ASSERT(nbEpsi .le. 6)

! - Current time
    allVarcStrain%time = time
    if (time .eq. r8vide()) then
        allVarcStrain%hasTime = ASTER_FALSE
    else
        allVarcStrain%hasTime = ASTER_TRUE
    end if

! - Loop on Gauss points
    do kpg = 1, npg
! ----- Initializations of material parameters on current integration point
        call initParaPoin(kpg, ksp, materPara)

! ----- Compute inelastic strains from external state variables on current integration point
        epsiVarcKpg = 0.d0
        call epstmc(materPara, "+", time, ndim, &
                    indxVarcStrain, allVarcStrain, &
                    epsiVarcKpg)
        do iEpsi = 1, nbEpsi
            epsiVarc(nbEpsi*(kpg-1)+iEpsi) = epsiVarcKpg(iEpsi)
        end do
    end do
!
end subroutine
