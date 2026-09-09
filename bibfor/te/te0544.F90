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
subroutine te0544(option, nomte)

    implicit none
!
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/te0544_implement.h"
#include "asterfort/teattr.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
!  CALCUL DES COEFFICIENTS A0 ET A1 POUR LE PILOTAGE PAR INCREMENT DE DEFORMATION
!  POUR LES ELEMENTS A VARIABLES LOCALES 2D et 3D
!
! --------------------------------------------------------------------------------------------------
    character(len=8), parameter :: fami = 'RIGI'
! --------------------------------------------------------------------------------------------------
    character(len=8):: typmod(2)
    character(len=16), pointer :: compor(:) => null()
    integer(kind=8) :: ndim, nno, npg
    integer(kind=8) :: jv_poids, jv_vff, jv_dfde, jv_geom
    integer(kind=8) :: jv_deplm, jv_ddepl, jv_depl0, jv_depl1
    integer(kind=8) :: jv_copil, jv_dtau
! --------------------------------------------------------------------------------------------------
!
! - Type of modelling
    call teattr('S', 'TYPMOD', typmod(1))
    call teattr('C', 'TYPMOD2', typmod(2), vattr_missing=' ')

! - Get parameters of element
    call elrefe_info(fami=fami, nno=nno, npg=npg, ndim=ndim, jpoids=jv_poids, &
                     jvf=jv_vff, jdfde=jv_dfde)

! - Common parameters
    call jevech('PGEOMER', 'L', jv_geom)
    call jevech('PDEPLMR', 'L', jv_deplm)
    call jevech('PDDEPLR', 'L', jv_ddepl)
    call jevech('PDEPL0R', 'L', jv_depl0)
    call jevech('PDEPL1R', 'L', jv_depl1)
    call jevech('PCOMPOR', 'L', vk16=compor)
    call jevech('PCDTAU', 'L', jv_dtau)
    call jevech('PCOPILO', 'E', jv_copil)

! - Main subroutine to compute coefficients
    call te0544_implement(typmod, compor, ndim, nno, npg, &
                          jv_poids, jv_vff, jv_dfde, zr(jv_geom), &
                          zr(jv_deplm), zr(jv_ddepl), zr(jv_depl0), zr(jv_depl1), &
                          zr(jv_dtau), zr(jv_copil))

end subroutine
