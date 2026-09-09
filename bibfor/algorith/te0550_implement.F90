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
! aslint: disable=W0104

subroutine te0550_implement(option, fami, nno, npg, ndim_sp, &
                            wref, vff, dxi_ff, aire, geom, pesa, &
                            mate, fext)
!
! --------------------------------------------------------------------------------------------------
!     CALCUL DES FORCES ELEMENTAIRES LINEIQUES POUR LES ELEMENTS BARRE
! --------------------------------------------------------------------------------------------------
!
! option : nom de l'option à calculer
!       CHAR_MECA_PESA_R    : charges de pesanteur
! --------------------------------------------------------------------------------------------------
!
    use tenseur_dime_module, only: proten
    implicit none

#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/dfdm1d_xd.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
!
    character(len=*) :: option
    character(len=8) :: fami
    integer(kind=8), intent(in):: nno, npg, ndim_sp, mate
    real(kind=8), intent(in):: geom(ndim_sp, nno), wref(npg), vff(nno, npg), dxi_ff(nno, npg)
    real(kind=8), intent(in):: aire, pesa(0:ndim_sp)
    real(kind=8), intent(out):: fext(ndim_sp, nno)
! --------------------------------------------------------------------------------------------------
    integer(kind=8):: cret(1), g
    real(kind=8):: w(npg), ds_ff(nno), rho(npg), grav(ndim_sp), fvol(ndim_sp, npg)
! --------------------------------------------------------------------------------------------------

! - Gauss points weight
    do g = 1, npg
        call dfdm1d_xd(dxi_ff(:, g), wref(g), geom, ds_ff, w(g))
    end do
    w = w*aire

! - Accélération de la pesanteur
    grav(:) = pesa(0)*pesa(1:ndim_sp)

! - density
    do g = 1, npg
        call rcvalb(fami, g, 1, '+', mate, ' ', 'ELAS', 0, ' ', [0.d0], 1, 'RHO', rho(g), cret, 1)
    end do

! - Volumetric force
    fvol(:, :) = proten(grav, rho)

! - External forces
    fext = matmul(fvol*spread(w, 1, ndim_sp), transpose(vff))

end subroutine
