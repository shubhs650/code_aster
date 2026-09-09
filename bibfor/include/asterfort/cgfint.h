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
    subroutine cgfint(BEHInteg, &
                      typmod, option, &
                      comporKit, carcrikit, &
                      ndim, nno1, nno2, npg, wref, &
                      vff1, vff2, dffr1, &
                      geom, tang, &
                      instam, instap, &
                      ddlm, ddld, &
                      iu, iuc, im, a, &
                      lgpg, sigm, vim, &
                      sigp, vip, &
                      matr, vect, &
                      codret)
        use Behaviour_type
        type(Behaviour_Integ), intent(inout) :: BEHInteg
        character(len=8), intent(in) :: typmod(2)
        character(len=16), intent(in) :: option, comporKit(COMPOR_SIZE)
        real(kind=8), intent(in) :: carcrikit(CARCRI_SIZE)
        integer(kind=8) :: lgpg
        integer(kind=8) :: npg
        integer(kind=8) :: nno2
        integer(kind=8) :: nno1
        integer(kind=8) :: ndim
        real(kind=8) :: wref(npg)
        real(kind=8) :: vff1(nno1, npg)
        real(kind=8) :: vff2(nno2, npg)
        real(kind=8) :: dffr1(nno1, npg)
        real(kind=8) :: geom(ndim, nno1)
        real(kind=8) :: tang(*)
        real(kind=8) :: instam
        real(kind=8) :: instap
        real(kind=8) :: ddlm(nno1*(ndim+1)+nno2)
        real(kind=8) :: ddld(nno1*(ndim+1)+nno2)
        integer(kind=8) :: iu(3, 3)
        integer(kind=8) :: iuc(3)
        integer(kind=8) :: im(3)
        real(kind=8) :: a
        real(kind=8) :: sigm(3, npg)
        real(kind=8) :: vim(lgpg, npg)
        real(kind=8) :: sigp(3, npg)
        real(kind=8) :: vip(lgpg, npg)
        real(kind=8) :: matr(*)
        real(kind=8) :: vect(nno1*(ndim+1)+nno2)
        integer(kind=8) :: codret
    end subroutine cgfint
end interface
