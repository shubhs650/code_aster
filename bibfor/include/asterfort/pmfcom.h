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
    subroutine pmfcom(materPara, &
                      option, carcri, &
                      kpg, debsp, pmfCompor, &
                      nf, instam, instap, nbvalc, &
                      defam, defap, varim, varimp, contm, &
                      defm, ddefp, epsm, modf, sigf, &
                      varip, codret)
        use MaterialPara_type
        type(Material_Para), intent(inout) :: materPara
        character(len=16), intent(in) :: option
        real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
        integer(kind=8) :: nbvalc
        integer(kind=8) :: nf
        integer(kind=8) :: kpg
        integer(kind=8) :: debsp
        character(len=24) :: pmfCompor(*)
        real(kind=8) :: instam
        real(kind=8) :: instap
        real(kind=8) :: defam(*)
        real(kind=8) :: defap(*)
        real(kind=8) :: varim(nbvalc*nf)
        real(kind=8) :: varimp(nbvalc*nf)
        real(kind=8) :: contm(nf)
        real(kind=8) :: defm(nf)
        real(kind=8) :: ddefp(nf)
        real(kind=8) :: epsm
        real(kind=8) :: modf(nf)
        real(kind=8) :: sigf(nf)
        real(kind=8) :: varip(nbvalc*nf)
        integer(kind=8) :: codret
    end subroutine pmfcom
end interface
