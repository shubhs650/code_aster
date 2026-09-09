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
    subroutine assesu(ds_thm, &
                      lMatr, lVect, lSigm, &
                      lVari, lMatrPred, &
                      option, typmod, &
                      ndim, nbvari, &
                      nno, nnos, nface, &
                      dimdef, dimcon, dimuel, &
                      mecani, press1, press2, tempe, &
                      compor, carcri, &
                      elem_coor, &
                      dispm, dispp, &
                      defgem, defgep, &
                      congem, congep, &
                      vintm, vintp, &
                      time_prev, time_curr, &
                      matuu, vectu)
        use THM_type
        use MaterialPara_type
        type(THM_DS), intent(inout) :: ds_thm
        aster_logical, intent(in) :: lVect, lMatr, lVari, lSigm, lMatrPred
        integer(kind=8), parameter :: maxfa = 6
        character(len=16), intent(in) :: option
        character(len=8), intent(in) :: typmod(2)
        integer(kind=8), intent(in) :: ndim, nbvari
        integer(kind=8), intent(in) :: nno, nnos, nface
        integer(kind=8), intent(in) :: dimdef, dimcon, dimuel
        integer(kind=8), intent(in) :: mecani(5), press1(7), press2(7), tempe(5)
        character(len=16), intent(in) :: compor(COMPOR_SIZE)
        real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
        real(kind=8), intent(in) :: elem_coor(ndim, nno)
        real(kind=8), intent(in) :: dispm(dimuel), dispp(dimuel)
        real(kind=8), intent(inout) :: defgem(dimdef), defgep(dimdef)
        real(kind=8), intent(in) :: congem(dimcon, maxfa+1)
        real(kind=8), intent(inout) :: congep(dimcon, maxfa+1)
        real(kind=8), intent(in) :: vintm(nbvari, maxfa+1)
        real(kind=8), intent(inout) :: vintp(nbvari, maxfa+1)
        real(kind=8), intent(in) :: time_curr, time_prev
        real(kind=8), intent(inout) :: matuu(dimuel*dimuel)
        real(kind=8), intent(inout) :: vectu(dimuel)
    end subroutine assesu
end interface
