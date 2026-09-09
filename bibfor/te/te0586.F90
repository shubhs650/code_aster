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
subroutine te0586(option, nomte)
!
    use pipeElem_module
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterfort/jevech.h"
#include "asterfort/pipeElem_type.h"
#include "asterfort/tufull.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: TUYAU
!
! Options: FULL_MECA_*, RIGI_MECA_*, RAPH_MECA
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = 'RIGI'
    integer(kind=8) :: nbFourier, nbDof, nbNode
    type(Material_Para) :: materPara
    integer(kind=8) :: jvMaterc
!
! --------------------------------------------------------------------------------------------------
!
    call pipeGetDime(nomte, fami, &
                     nbNode, nbFourier, nbDof)

! - Material parameters
    call jevech('PMATERC', 'L', jvMaterc)

! - Initializations of material parameters on current cell
    call initParaCell(fami, zi(jvMaterc), materPara)

!   Angle du mot clef MASSIF de AFFE_CARA_ELEM, initialisé à 0, nécessaire pour les LdC
! - LEMAITRE_IRRA et VISC_IRRA_LOG (voir ssnl121c)
    call initLCSZero(materPara)

! - Compute option
    call tufull(materPara, option, nbFourier, nbDof)
!
end subroutine
