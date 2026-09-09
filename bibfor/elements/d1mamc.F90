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
!
subroutine d1mamc(materPara, poum, time, nbSig, d1)
!
    use MaterialPara_type
    implicit none
!
#include "asterfort/d1ma3d.h"
#include "asterfort/d1macp.h"
#include "asterfort/d1madp.h"
#include "asterfort/lteatt.h"
#include "asterfort/utmess.h"
!
    type(Material_Para), intent(in) :: materPara
    character(len=*), intent(in) :: poum
    real(kind=8), intent(in) :: time
    integer(kind=8), intent(in) :: nbsig
    real(kind=8), intent(out) :: d1(nbsig, 1)
!
! --------------------------------------------------------------------------------------------------
!
!      D1MAMC :   CALCUL DE L'INVERSE DE LA MATRICE DE HOOKE
!                 POUR LES ELEMENTS ISOPARAMETRIQUES POUR DES
!                 MATERIAUX ISOTROPE, ORTHOTROPE ET ISOTROPE TRANSVERSE
!
! --------------------------------------------------------------------------------------------------
!
    if (lteatt('DIM_TOPO_MAILLE', '3') .or. lteatt('FOURIER', 'OUI')) then
        call d1ma3d(materPara, poum, time, d1)

    elseif (lteatt('D_PLAN', 'OUI') .or. lteatt('AXIS', 'OUI')) then
        call d1madp(materPara, poum, time, d1)

    else if (lteatt('C_PLAN', 'OUI')) then
        call d1macp(materPara, poum, time, d1)

    else
        call utmess('F', 'ELEMENTS_11')
    end if
!
end subroutine
