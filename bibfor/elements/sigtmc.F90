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
subroutine sigtmc(materPara, time, &
                  nbsig, npg, ndim, &
                  indxVarcStrain, sigmVarc)
!
    use BehaviourStrain_module
    use BehaviourStrain_type
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterc/r8miem.h"
#include "asterfort/assert.h"
#include "asterfort/dmatmc.h"
#include "asterfort/epstmc.h"
!
    type(Material_Para), intent(inout) :: materPara
    real(kind=8), intent(in) :: time
    integer(kind=8), intent(in) :: nbsig, npg, ndim
    integer(kind=8), intent(in) :: indxVarcStrain
    real(kind=8), intent(out) :: sigmVarc(nbsig*npg)
!
! --------------------------------------------------------------------------------------------------
!
! Compute stresses from external state variables
!
! --------------------------------------------------------------------------------------------------
!
! IO  materPara        : parameters of material
! In  time             : current time
! In  nbsig            : number of components for stress tensors (4 or 6)
! In  npg              : number of Gauss points
! In  indxVarcStrain   : index of external state variable
! Out sigmVarc         : Stresses from external state variables
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: ksp = 1
    real(kind=8) :: d(36), epsiVarc(6)
    integer(kind=8) :: k
    aster_logical :: hasEpsiVarc
    integer(kind=8) :: iSigm, kpg, jSigm
    type(All_Varc_Strain) :: allVarcStrain
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nbsig .le. 6)
    sigmVarc = 0.d0

! - Loop on Gauss points
    do kpg = 1, npg
! ----- Initializations of material parameters on current integration point
        call initParaPoin(kpg, ksp, materPara)

! ----- Compute inelastic strains from external state variables
        call epstmc(materPara, "+", time, ndim, &
                    indxVarcStrain, allVarcStrain, &
                    epsiVarc)

! ----- Compute
        hasEpsiVarc = ASTER_FALSE
        do k = 1, 6
            if (abs(epsiVarc(k)) .gt. r8miem()) then
                hasEpsiVarc = ASTER_TRUE
            end if
        end do
        if (hasEpsiVarc) then
            epsiVarc(4:6) = 2.d0*epsiVarc(4:6)

! --------- Compute Hooke matrix [D]
            call dmatmc(materPara, "+", time, nbsig, d)

! --------- Compute stresses from external state variables
            do iSigm = 1, nbsig
                do jSigm = 1, nbsig
                    sigmVarc(iSigm+nbsig*(kpg-1)) = sigmVarc(iSigm+nbsig*(kpg-1))+ &
                                                    d(jSigm+(iSigm-1)*nbsig)*epsiVarc(jSigm)
                end do
            end do
        end if
    end do
!
end subroutine
