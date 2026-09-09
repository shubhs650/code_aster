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
    subroutine nmas2d(BEHinteg, &
                      nno, npg, &
                      ipoids, ivf, idfde, &
                      geom, typmod, option, &
                      compor, mult_comp, lgpg, carcri, instam, instap, &
                      deplm, deplp, sigm, vim, &
                      dfdi, def, sigp, vip, matuu, &
                      vectu, codret)
        use Behaviour_type
        type(Behaviour_Integ), intent(inout) :: BEHinteg
        integer(kind=8) :: lgpg
        integer(kind=8) :: npg
        integer(kind=8) :: nno
        integer(kind=8) :: ipoids
        integer(kind=8) :: ivf
        integer(kind=8) :: idfde
        real(kind=8) :: geom(2, nno)
        character(len=8) :: typmod(2)
        character(len=16) :: option
        character(len=16), intent(in) :: compor(COMPOR_SIZE)
        character(len=16), intent(in) :: mult_comp
        real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
        real(kind=8) :: instam
        real(kind=8) :: instap
        real(kind=8) :: deplm(2, nno), deplp(2, nno)
        real(kind=8) :: sigm(10, npg)
        real(kind=8) :: vim(lgpg, npg)
        real(kind=8) :: dfdi(nno, 2)
        real(kind=8) :: def(4, nno, 2)
        real(kind=8) :: sigp(10, npg)
        real(kind=8) :: vip(lgpg, npg)
        real(kind=8) :: matuu(*)
        real(kind=8) :: vectu(2, nno)
        integer(kind=8) :: codret
    end subroutine nmas2d
end interface
