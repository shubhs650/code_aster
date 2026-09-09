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
    subroutine varcCalcPrep(modelZ, caraElemZ, materCodeZ, &
                            poum, &
                            l_temp, l_meta, &
                            varcRefeZ, varcPrevZ, varcCurrZ, &
                            comporZ, multCompZ, chsithz, &
                            sigmz, variz, &
                            nbFieldInMax, nbFieldOutMax, &
                            nbFieldIn, nbFieldOut, &
                            lpain, lchin, &
                            lpaout, lchout)
        character(len=*), intent(in) :: modelZ, caraElemZ, materCodeZ
        aster_logical, intent(in) :: l_temp, l_meta
        character(len=1), intent(in) :: poum
        character(len=*), intent(in) :: varcRefeZ, varcPrevZ, varcCurrZ
        character(len=*), intent(in) :: comporZ, multCompZ, chsithz
        character(len=*), intent(in) :: sigmz, variz
        integer(kind=8), intent(in) :: nbFieldInMax, nbFieldOutMax
        integer(kind=8), intent(out) :: nbFieldIn, nbFieldOut
        character(len=8), intent(out) :: lpaout(nbFieldOutMax), lpain(nbFieldInMax)
        character(len=19), intent(out) :: lchout(nbFieldOutMax), lchin(nbFieldInMax)
    end subroutine varcCalcPrep
end interface
