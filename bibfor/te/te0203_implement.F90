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
subroutine te0203_implement(BEHInteg, typmod, compor, ndim, nno, npg, geom, &
                            deplm, ddepl, ddepl0, ddepl1, &
                            lgpg, sigm, vim, etamin, etamax, dtau, &
                            copilo)

    use Behaviour_type
    implicit none

#include "asterf_types.h"
#include "asterc/r8vide.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/nmfisa.h"
#include "asterfort/pi0000.h"

    type(Behaviour_Integ) :: BEHInteg
    character(len=8) :: typmod(2)
    character(len=16) ::compor(COMPOR_SIZE)
    integer(kind=8), intent(in):: ndim
    integer(kind=8), intent(in):: nno
    integer(kind=8), intent(in):: npg
    real(kind=8), intent(in) :: geom(ndim, nno)
    real(kind=8), intent(in) :: deplm(ndim*nno)
    real(kind=8), intent(in) :: ddepl(ndim*nno)
    real(kind=8), intent(in) :: ddepl0(ndim*nno)
    real(kind=8), intent(in) :: ddepl1(ndim*nno)
    integer(kind=8), intent(in) :: lgpg
    real(kind=8), intent(in) :: sigm(ndim, npg)
    real(kind=8), intent(in) :: vim(lgpg, npg)
    real(kind=8), intent(in) :: etamin
    real(kind=8), intent(in) :: etamax
    real(kind=8), intent(in) :: dtau
    real(kind=8), intent(out) :: copilo(5, npg)
! --------------------------------------------------------------------------------------------------
!  PILOTAGE PRED_ELAS POUR LES ELEMENTS DE JOINT 2D
! --------------------------------------------------------------------------------------------------
    aster_logical :: axi
    integer(kind=8) :: kpg
    real(kind=8) :: epsm(ndim), epsd_cste(ndim), epsd_pilo(ndim), b(ndim, ndim*nno), poids
! --------------------------------------------------------------------------------------------------

    axi = typmod(1) .eq. 'AXIS'
    copilo = r8vide()

    do kpg = 1, npg

        ! Kinematics
        call nmfisa(axi, geom, kpg, poids, b)
        epsm = matmul(b, deplm)
        epsd_cste = matmul(b, ddepl+ddepl0)
        epsd_pilo = matmul(b, ddepl1)

        ! Path-following coefficients
        call pi0000(BEHInteg, compor, typmod, ndim, &
                    epsm, epsd_cste, epsd_pilo, &
                    sigm(:, kpg), vim(:, kpg), dtau, etamin, etamax, &
                    copilo(:, kpg))

    end do

end subroutine
