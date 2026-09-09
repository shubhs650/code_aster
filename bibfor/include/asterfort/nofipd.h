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
    subroutine nofipd(BEHInteg, &
                      ndim, nnod, nnop, nnog, npg, &
                      iw, vffd, vffp, vffg, idffd, &
                      vu, vp, vpi, &
                      geomi, typmod, option, nomte, compor, &
                      lgpg, carcri, instm, instp, &
                      ddlm, ddld, &
                      sigm, vim, sigp, vip, &
                      vect, matr, codret, &
                      lSigm, lVect, lMatr)
        use Behaviour_type
        type(Behaviour_Integ), intent(inout) :: BEHinteg
        character(len=8), intent(in) :: typmod(2)
        character(len=16), intent(in) :: compor(COMPOR_SIZE)
        real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
        integer(kind=8) :: ndim, nnod, nnop, nnog, npg, iw, idffd, lgpg
        integer(kind=8) :: vu(3, 27), vp(27), vpi(3, 27)
        integer(kind=8) :: codret
        real(kind=8) :: vffd(nnod, npg), vffp(nnop, npg), vffg(nnog, npg)
        real(kind=8) :: instm, instp
        real(kind=8) :: geomi(ndim, nnod), ddlm(*), ddld(*)
        real(kind=8) :: sigm(2*ndim+1, npg), sigp(2*ndim+1, npg)
        real(kind=8) :: vim(lgpg, npg), vip(lgpg, npg)
        real(kind=8) :: vect(*), matr(*)
        character(len=16), intent(in) :: option
        character(len=16) :: nomte
        aster_logical, intent(in) :: lSigm, lVect, lMatr
    end subroutine nofipd
end interface
