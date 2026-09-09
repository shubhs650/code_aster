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
subroutine ffpoutimo(x, xl, mate, materi, ff)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/matela.h"
#include "asterfort/poutre_modloc.h"
!
    real(kind=8), intent(in)       :: x(3), xl
    integer(kind=8), intent(in)    :: mate
    character(len=8), intent(in)   :: materi
    real(kind=8), intent(out)      :: ff(18)
!
! --------------------------------------------------------------------------------------------------
!
! Finite elements management
!
! Value of shape functions at given point for Timoshenko Beam
!
! --------------------------------------------------------------------------------------------------
!
! In  x       : coordinates in parametric space to evaluate shape function
! In  xl      : beam element length
! In  mate    : address of coded material
! In  materi  : name of coded material
! Out ff      : value of shape functions at point x
!
! --------------------------------------------------------------------------------------------------
!
! --------------------------------------------------------------------------------------------------
!
! parameters depuis porigi
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: e, g, nu
    real(kind=8) :: a, alfay, alfaz
    real(kind=8) :: xiy, xiz
    real(kind=8) :: xx, xi
    real(kind=8) :: xl2, eiy, eiz, phiy, phiz
! --------------------------------------------------------------------------------------------------

    integer(kind=8), parameter :: nb_cara = 11
    real(kind=8) :: vale_cara(nb_cara)
    character(len=8) :: noms_cara(nb_cara)
    data noms_cara/'A1', 'IY1', 'IZ1', 'AY1', 'AZ1', 'EY1', 'EZ1', 'JX1', 'JG1', &
        'IYR21', 'IZR21'/

!
! --------------------------------------------------------------------------------------------------
!

!   Récupération des aramètres matériau sans prise en compte de la température
    call matela(mate, materi, 0, 0.d0, e, nu)
    g = e/(2.0d0*(1.0d0+nu))

!   Récuperation des caracteristiques generales des sections
    call poutre_modloc('CAGNP1', noms_cara, nb_cara, lvaleur=vale_cara)
    a = vale_cara(1)
    xiy = vale_cara(2)
    xiz = vale_cara(3)
    alfay = vale_cara(4)
    alfaz = vale_cara(5)

!   Calcul de constantes
    xl2 = xl*xl
    eiy = e*xiy
    eiz = e*xiz
    phiy = (12.d0*eiz*alfay)/(g*a*xl2)
    phiz = (12.d0*eiy*alfaz)/(g*a*xl2)

!   Evaluation des fonctions de forme aux points de Gauss
!       elles sont donc écrites en xi telles que : N_xx((1+xi)*Lp/2)
    xi = x(1)
    xx = (1+xi)*xl/2

    ff(1) = 1-xx/xl
    ff(2) = 1/(1+phiz)*(2*(xx/xl)**3-3*(xx/xl)**2-phiz*(xx/xl)+1+phiz)
    ff(3) = xl/(1+phiz)*((xx/xl)**3-(2+phiz/2)*(xx/xl)**2+(1+phiz/2)*(xx/xl))
    ff(4) = xx/xl
    ff(5) = -1/(1+phiz)*(2*(xx/xl)**3-3*(xx/xl)**2-phiz*(xx/xl))
    ff(6) = xl/(1+phiz)*((xx/xl)**3-(1-phiz/2)*(xx/xl)**2-phiz/2*(xx/xl))
    ff(7) = 6/((1+phiz)*xl)*((xx/xl)**2-xx/xl)
    ff(8) = 1/(1+phiz)*(3*(xx/xl)**2-(4+phiz)*xx/xl+1+phiz)
    ff(9) = -6/((1+phiz)*xl)*((xx/xl)**2-xx/xl)
    ff(10) = 1/(1+phiz)*(3*(xx/xl)**2-(2-phiz)*xx/xl)

    ff(11) = 1/(1+phiy)*(2*(xx/xl)**3-3*(xx/xl)**2-phiy*(xx/xl)+1+phiy)
    ff(12) = -xl/(1+phiy)*((xx/xl)**3-(2+phiy/2)*(xx/xl)**2+(1+phiy/2)*(xx/xl))
    ff(13) = -1/(1+phiy)*(2*(xx/xl)**3-3*(xx/xl)**2-phiy*(xx/xl))
    ff(14) = -xl/(1+phiy)*((xx/xl)**3-(1-phiy/2)*(xx/xl)**2-phiy/2*(xx/xl))
    ff(15) = 6/((1+phiy)*xl)*((xx/xl)**2-xx/xl)
    ff(16) = 1/(1+phiy)*(3*(xx/xl)**2-(4+phiy)*xx/xl+1+phiy)
    ff(17) = -6/((1+phiy)*xl)*((xx/xl)**2-xx/xl)
    ff(18) = 1/(1+phiy)*(3*(xx/xl)**2-(2-phiy)*xx/xl)
!
end subroutine
