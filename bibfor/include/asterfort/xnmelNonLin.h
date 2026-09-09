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
    subroutine xnmelNonLin(BEHInteg, &
                           option, typmod, &
                           compor, carcri, &
                           nnop, nfh, nfe, &
                           ddlc, ddlm, jvGeom, &
                           lgpg, jpintt, cnset, heavt, lonch, basloc, &
                           instam, instap, idepl, &
                           lsn, lst, sig, vi, matuu, ivectu, &
                           codret, jpmilt, nfiss, jheavn, jstno, &
                           lMatr, lVect, lSigm)
        use Behaviour_type
        type(Behaviour_Integ), intent(inout) :: BEHInteg
        character(len=8), intent(in) :: typmod(2)
        character(len=16), intent(in) :: option, compor(COMPOR_SIZE)
        real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
        integer(kind=8) :: nnop, imate, lgpg, codret, jvGeom, nfiss, jheavn
        integer(kind=8) :: cnset(4*32), heavt(*), lonch(10), ndim
        integer(kind=8) :: nfh, nfe, ddlc, ddlm
        integer(kind=8) :: ivectu, idepl, jpintt, jpmilt
        integer(kind=8) :: jstno
        real(kind=8) :: instam, instap
        real(kind=8) :: vi(*), crit2(1), vi2(1), sig2(1)
        real(kind=8) :: lsn(nnop)
        real(kind=8) :: lst(nnop), matuu(*), sig(*), basloc(*)
        aster_logical, intent(in) :: lMatr, lVect, lSigm
    end subroutine xnmelNonLin
end interface
