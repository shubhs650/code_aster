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
#include "asterfort/Behaviour_type.h"
!
interface
    subroutine nmiclb(materPara, &
                      option, relaComp, carcri, &
                      xlong0, aire, tmoins, tplus, &
                      dlong0, effnom, vim, effnop, vip, &
                      klv, fono, epsm, codret)
        use MaterialPara_type
        type(Material_Para), intent(in) :: materPara
        character(len=16), intent(in) :: option, relaComp
        real(kind=8) :: carcri(CARCRI_SIZE)
        real(kind=8) :: xlong0
        real(kind=8) :: aire
        real(kind=8) :: tmoins
        real(kind=8) :: tplus
        real(kind=8) :: dlong0
        real(kind=8) :: effnom
        real(kind=8) :: vim(*)
        real(kind=8) :: effnop
        real(kind=8) :: vip(*)
        real(kind=8) :: klv(21)
        real(kind=8) :: fono(6)
        real(kind=8) :: epsm
        integer(kind=8) :: codret
    end subroutine nmiclb
end interface
