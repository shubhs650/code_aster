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
    subroutine nmgpfi(BEHInteg, &
                      typmod, option, &
                      nno, npg, ndim, geomInit, &
                      compor, carcri, multComp, &
                      instm, instp, &
                      dispPrev, dispIncr, &
                      lgpg, sigmPrev, vim, sigmCurr, &
                      vip, fint, &
                      matr, codret)
        use Behaviour_type
        type(Behaviour_Integ), intent(inout) :: BEHInteg
        integer(kind=8) :: ndim, nno, npg, lgpg
        character(len=8), intent(in) :: typmod(2)
        character(len=16), intent(in) :: option
        character(len=16), intent(in) :: compor(COMPOR_SIZE)
        real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
        character(len=16), intent(in) :: multComp
        real(kind=8) :: geomInit(*), instm, instp
        real(kind=8) :: dispPrev(*), dispIncr(*), sigmPrev(2*ndim, npg)
        real(kind=8) :: vim(lgpg, npg), sigmCurr(2*ndim, npg), vip(lgpg, npg)
        real(kind=8) :: matr(*), fint(*)
        integer(kind=8), intent(out) :: codret
    end subroutine nmgpfi
end interface
