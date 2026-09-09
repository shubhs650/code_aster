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
! aslint: disable=W0413
!
subroutine dxroep(plateCara, rho, epais)
!
    use plate_type
    implicit none
!
#include "asterc/r8maem.h"
#include "asterfort/jevech.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvala.h"
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    type(plateCara_Para), intent(in) :: plateCara
    real(kind=8), intent(out) :: rho, epais
!
! --------------------------------------------------------------------------------------------------
!
!     APPEL DES MASSE VOLUMIQUE DU MATERIAU ET EPAISSEUR DE LA PLAQUE
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbPropMaxi = 2
    real(kind=8) :: propVale(nbPropMaxi)
    integer(kind=8) :: propCode(nbPropMaxi)
    character(len=16) :: propName(nbPropMaxi)
    integer(kind=8) :: jvMaterc, nbProp, iadzi, iazk24
    real(kind=8) :: r8bid
    character(len=32) :: elasKeyword
!
! --------------------------------------------------------------------------------------------------
!
    r8bid = 0.d0
    call jevech('PMATERC', 'L', jvMaterc)
    call rccoma(zi(jvMaterc), 'ELAS', 1, elasKeyword)
!
    if (elasKeyword .eq. 'ELAS_COQMU') then
        propName(1) = 'HOM_19'
        propName(2) = 'HOM_20'
        nbProp = 2
        call rcvala(zi(jvMaterc), ' ', elasKeyword, &
                    0, ' ', [r8bid], &
                    nbProp, propName, &
                    propVale, propCode, 1)
        epais = propVale(1)
        rho = propVale(2)
        if (rho .eq. r8maem()) then
            call tecael(iadzi, iazk24)
            call utmess('F', 'ELEMENTS4_81', sk='RHO', si=zi(iadzi-1+1))
        end if

    elseif (elasKeyword .eq. 'ELAS' .or. elasKeyword .eq. 'ELAS_COQUE' .or. &
            elasKeyword .eq. 'ELAS_ISTR' .or. elasKeyword .eq. 'ELAS_ORTH' .or. &
            elasKeyword .eq. 'ELAS_GLRC' .or. elasKeyword .eq. 'ELAS_DHRC') then
        propName(1) = 'RHO'
        nbProp = 1
        call rcvala(zi(jvMaterc), ' ', elasKeyword, &
                    0, ' ', [r8bid], &
                    nbProp, propName, &
                    propVale, propCode, 1)
        rho = propVale(1)
        epais = plateCara%thick

    else
        call utmess('F', 'ELEMENTS_50')
    end if
!
end subroutine
