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
    subroutine pmf_vmis(for_pmf, nf, nbvalc, &
                        pmfCompor, materPara, &
                        varim, contm, &
                        ddefp, modf, &
                        sigf, varip, codret)
        use pmfcom_type
        use MaterialPara_type
        type(pmfcom_user), intent(in) :: for_pmf
        integer(kind=8) :: nf, nbvalc
        character(len=24) :: pmfCompor(*)
        type(Material_Para), intent(inout) :: materPara
        real(kind=8) :: varim(nbvalc*nf), contm(nf), ddefp(nf), modf(nf)
        real(kind=8) :: sigf(nf), varip(nbvalc*nf)
        integer(kind=8) :: codret
    end subroutine pmf_vmis
end interface
