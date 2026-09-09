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
subroutine ethdst(materPara, &
                  nno, ndim, nbsig, npg, &
                  jvGaussWeight, jvBaseFunc, jvDBaseFunc, &
                  nodeCoor, time, &
                  enerTherTher)
!
    use BehaviourStrain_type
    use MaterialPara_type
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/epthmc.h"
#include "asterfort/lteatt.h"
#include "asterfort/sigtmc.h"
#include "MeshTypes_type.h"
#include "jeveux.h"
!
    type(Material_Para), intent(inout) :: materPara
    integer(kind=8), intent(in) :: nno, ndim, nbsig, npg
    integer(kind=8), intent(in) :: jvGaussWeight, jvBaseFunc, jvDBaseFunc
    real(kind=8), intent(in) :: nodeCoor(ndim*nno)
    real(kind=8), intent(in) :: time
    real(kind=8), intent(out) :: enerTherTher

! --------------------------------------------------------------------------------------------------
!
!      ETHDST   -- CALCUL DU TERME EPSTHT*D*EPSTH RENTRANT
!                  DANS LE CALCUL DE L'ENERGIE POTENTIELLE
!                  (I.E.  1/2*UT*K*U - UT*FTH + 1/2*EPSTHT*D*EPSTH)
!                  POUR LES ELEMENTS ISOPARAMETRIQUES
!
! --------------------------------------------------------------------------------------------------
!
! IO  materPara        : parameters of material
! In  nno              : number of nodes
! In  ndim             : dimension of space
! In  nbsig            : number of stress tensor components
! In  npg              : number of Gauss points
! In  jvGaussWeight    : adresse to Gauss point weight
! In  jvBaseFunc       : adresse to shape functions
! In  jvDBaseFunc      : adresse to derivative of shape functions
! In  nodeCoor         : coordinates of nodes
! In  time             : given time
! Out enerTherTher     : SOMME(EPSTH_T*D*EPSTH)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ino, isig, kpg, nbEpsi
    real(kind=8) :: epsiTher(162), sigmTher(162)
    real(kind=8) :: rayon
    real(kind=8) :: enerTherKpg, dfdx(MT_NNOMAX3D), dfdy(MT_NNOMAX3D), dfdz(MT_NNOMAX3D)
    real(kind=8) :: jacobKpg
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nbsig .le. 6)
    ASSERT(npg .le. 27)
    nbEpsi = nbsig
    enerTherTher = 0.d0

! - Calcul des déformations thermiques
    call epthmc(materPara, time, &
                nbEpsi, npg, ndim, &
                VARC_STRAIN_TEMP, epsiTher)

! - Calcul des contraintes thermiques
    call sigtmc(materPara, time, &
                nbsig, npg, ndim, &
                VARC_STRAIN_TEMP, sigmTher)

! - Loop on Gauss points
    do kpg = 1, npg
        enerTherKpg = 0.d0

! ----- Get jacobian at current Gauss point
        if (lteatt('DIM_TOPO_MAILLE', '3')) then
            call dfdm3d(nno, kpg, &
                        jvGaussWeight, jvDBaseFunc, nodeCoor, &
                        jacobKpg, dfdx, dfdy, dfdz)
        else
            call dfdm2d(nno, kpg, &
                        jvGaussWeight, jvDBaseFunc, nodeCoor, &
                        jacobKpg, dfdx, dfdy)
            if (lteatt('AXIS', 'OUI')) then
                rayon = 0.d0
                do ino = 1, nno
                    rayon = rayon+zr(jvBaseFunc+(kpg-1)*nno-1+ino)*nodeCoor(2*(ino-1)+1)
                end do
                jacobKpg = jacobKpg*rayon
            end if
        end if
        do isig = 1, nbsig
            enerTherKpg = enerTherKpg+ &
                          epsiTher(isig+nbsig*(kpg-1))*sigmTher(isig+nbsig*(kpg-1))
        end do
        enerTherTher = enerTherTher+(enerTherKpg*jacobKpg)
    end do
!
end subroutine
