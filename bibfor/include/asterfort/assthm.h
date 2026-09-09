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
!
interface
    subroutine assthm(ds_thm, &
                      lMatr, lSigm, lVect, &
                      lVari, lMatrPred, l_axi, &
                      option, typmod, inte_type, &
                      ndim, nbvari, &
                      nno, nnos, npg, npi, &
                      nddls, nddlm, nddl_meca, nddl_p1, nddl_p2, nddl_2nd, &
                      dimdef, dimcon, dimuel, &
                      mecani, press1, press2, tempe, second, &
                      compor, carcri, &
                      jv_poids, jv_poids2, &
                      jv_func, jv_func2, &
                      jv_dfunc, jv_dfunc2, &
                      elem_coor, &
                      dispm, dispp, &
                      congem, congep, &
                      vintm, vintp, &
                      time_prev, time_curr, &
                      matuu, vectu, codret)
        use THM_type
        type(THM_DS), intent(inout) :: ds_thm

        aster_logical, intent(in) :: lMatr, lSigm, lVari, lMatrPred, lVect
        aster_logical, intent(in)  :: l_axi
        character(len=16), intent(in) :: option
        character(len=8), intent(in) :: typmod(2)
        character(len=3), intent(in) :: inte_type
        integer(kind=8), intent(in) :: nbvari, ndim
        integer(kind=8), intent(in) :: nno, nnos
        integer(kind=8), intent(in) :: npg, npi
        integer(kind=8), intent(in) :: nddls, nddlm, nddl_meca, nddl_p1, nddl_p2, nddl_2nd
        integer(kind=8), intent(in) :: dimuel, dimdef, dimcon
        integer(kind=8), intent(in) :: mecani(5), press1(7), press2(7), tempe(5), second(5)
        character(len=16), intent(in)  :: compor(*)
        real(kind=8), intent(in) :: carcri(*)
        integer(kind=8), intent(in) :: jv_poids, jv_poids2
        integer(kind=8), intent(in) :: jv_func, jv_func2
        integer(kind=8), intent(in) :: jv_dfunc, jv_dfunc2
        real(kind=8), intent(in) :: elem_coor(ndim, nno)
        real(kind=8), intent(in) :: dispm(dimuel), dispp(dimuel)
        real(kind=8), intent(inout) :: congem(dimcon*npi)
        real(kind=8), intent(inout) :: congep(dimcon*npi)
        real(kind=8), intent(in) :: vintm(nbvari*npi)
        real(kind=8), intent(inout) :: vintp(nbvari*npi)
        real(kind=8), intent(in) :: time_prev, time_curr
        real(kind=8), intent(inout) :: matuu(dimuel*dimuel)
        real(kind=8), intent(inout) :: vectu(dimuel)
        integer(kind=8), intent(out) :: codret
    end subroutine assthm
end interface
