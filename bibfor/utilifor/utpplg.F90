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

subroutine utpplg(nn, nc, p, sl, sg)
    implicit none
#include "asterfort/mapvec.h"
#include "asterfort/mavec.h"
#include "asterfort/upletr.h"
#include "asterfort/uplstr.h"
#include "asterfort/utpalg.h"
#include "asterfort/utpslg.h"
#include "asterfort/vecmap.h"
    integer(kind=8), intent(in) :: nn, nc
    real(kind=8), intent(in) :: p(3, 3), sl(*)
    real(kind=8), intent(out) :: sg(*)
!
    integer(kind=8) :: n, n1, nddl
    real(kind=8), dimension(nn*nc, nn*nc) :: matsy1, matsy2, matas2
    real(kind=8), dimension(nn*nc, nn*nc) :: matsym, matasy
    real(kind=8), dimension(nn*nc, nn*nc) :: parsmg, parayg, matril, matrig
    real(kind=8), dimension(nn*nc*(nn*nc+1)/2) :: parsym, parasy, vecsym, vecasy
!     ------------------------------------------------------------------
!     PASSAGE D'UNE MATRICE PLEINE NONSYMETRIQUE DE NN*NC LIGNES
!     DU REPERE LOCAL AU REPERE GLOBAL (3D)
!     ------------------------------------------------------------------
!IN   I   NN   NOMBRE DE NOEUDS
!IN   I   N    NOMBRE DE NOEUDS
!IN   I   NC   NOMBRE DE COMPOSANTES
!IN   R   P    MATRICE DE PASSAGE 3D DE GLOBAL A LOCAL
!IN   R   SL   NN*NC COMPOSANTES DE LA MATRICE PLEINE SL DANS LOCAL
!OUT  R   SG   NN*NC COMPOSANTES DE LA MATRICE PLEINE SG DANS GLOBAL
!     ------------------------------------------------------------------
!
    nddl = nn*nc
    n = nddl*nddl
    n1 = (nddl+1)*nddl/2
!
    call vecmap(sl, n, matril, nddl)
    matsy1 = transpose(matril)
    matsy2 = matril+matsy1
    matas2 = matril-matsy1
    matsym = 0.5d0*matsy2
    call mavec(matsym, nddl, vecsym, n1)
    matasy = 0.5d0*matas2
    call mavec(matasy, nddl, vecasy, n1)
    call utpslg(nn, nc, p, vecsym, parsym)
    call uplstr(nddl, parsmg, parsym)
    call utpalg(nn, nc, p, vecasy, parasy)
    call upletr(nddl, parayg, parasy)
    matrig = parsmg+parayg
    call mapvec(matrig, nddl, sg, n)
!
!
end subroutine
