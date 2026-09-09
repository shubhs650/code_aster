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
subroutine te0119(option, nomte)
!
    use pipeElem_module
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystNone
    implicit none
!
#include "asterc/r8prem.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/teattr.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
!  VERI_CARA_ELEM
!
!  Vérification du contenu des cartes sur les éléments
!
!  Remarque : Si la vérification est faite içi, il est inutile de la faire dans "veri_affe_carte"
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ibid, jvCodret
    real(kind=8) :: excent
    character(len=3) :: cmod
    character(len=8) :: alias8
    aster_logical :: lPipe
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!
    call teattr('S', 'ALIAS8', alias8, ibid)
    cmod = alias8(3:5)
    lPipe = lteatt('TUYAU', 'OUI')

! - Vérification que l'excentrement est nul pour COQUE_3D
    if (cmod .eq. 'CQ3') then
        call getCara(plateCara, plateOrie)
        call compCoorSystNone(plateOrie)
        call jevech('PCODRET', "E", jvCodret)
        excent = plateCara%offset
        if (abs(excent) .ge. r8prem()) then
            zi(jvCodret-1+1) = 1
        end if
    end if

! - Checks for TUYAU
    if (lPipe) then
        call pipeCheckMetric(nomte)
    end if
!
end subroutine
