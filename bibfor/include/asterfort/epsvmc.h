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
interface
    subroutine epsvmc(nno, ndim, nbEpsi, npg, &
                      jvGaussWeight, jvBaseFunc, jvDBaseFunc, &
                      nodeCoor, nodeDisp, &
                      time, nharm, &
                      strainType, lStrainMeca, materPara, &
                      epsi)
        use MaterialPara_type
        integer(kind=8), intent(in) :: nno, ndim, nbEpsi, npg
        integer(kind=8), intent(in) :: jvGaussWeight, jvBaseFunc, jvDBaseFunc
        real(kind=8), intent(in) :: nodeCoor(ndim*nno), nodeDisp(ndim*nno)
        real(kind=8), intent(in) :: time, nharm
        integer(kind=8), intent(in) :: strainType
        aster_logical, intent(in) :: lStrainMeca
        type(Material_Para), intent(inout) :: materPara
        real(kind=8), intent(out) :: epsi(nbEpsi*npg)
    end subroutine epsvmc
end interface
