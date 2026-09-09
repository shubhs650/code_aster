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
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
!
interface
    subroutine nmfi2d(BEHInteg, &
                      npg, lgpg, option, geom, &
                      deplm, ddepl, sigmo, sigma, fint, &
                      kTang, vim, vip, &
                      tm, tp, &
                      carcri, compor, typmod, &
                      lMatr, lVect, lSigm, &
                      codret)
        use Behaviour_type
        type(Behaviour_Integ), intent(inout) :: BEHinteg
        integer(kind=8) :: lgpg
        integer(kind=8) :: npg
        character(len=8), intent(in) :: typmod(2)
        character(len=16), intent(in) :: option, compor(COMPOR_SIZE)
        real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
        real(kind=8) :: geom(2, 4)
        real(kind=8) :: deplm(8)
        real(kind=8) :: ddepl(8)
        real(kind=8) :: sigmo(6, npg)
        real(kind=8) :: sigma(6, npg)
        real(kind=8) :: fint(8)
        real(kind=8) :: kTang(8, 8)
        real(kind=8) :: vim(lgpg, npg)
        real(kind=8) :: vip(lgpg, npg)
        real(kind=8) :: tm
        real(kind=8) :: tp
        aster_logical, intent(in) :: lMatr, lVect, lSigm
        integer(kind=8) :: codret
    end subroutine nmfi2d
end interface
