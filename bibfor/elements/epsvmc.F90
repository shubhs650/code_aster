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
subroutine epsvmc(nno, ndim, nbEpsi, npg, &
                  jvGaussWeight, jvBaseFunc, jvDBaseFunc, &
                  nodeCoor, nodeDisp, &
                  time, nharm, &
                  strainType, lStrainMeca, materPara, &
                  epsi)
!
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/dmatmc.h"
#include "asterfort/eps1mc.h"
#include "asterfort/eps2mc.h"
#include "asterfort/epslmc.h"
#include "asterfort/epthmc.h"
#include "asterfort/lteatt.h"
#include "jeveux.h"
!
    integer(kind=8), intent(in) :: nno, ndim, nbEpsi, npg
    integer(kind=8), intent(in) :: jvGaussWeight, jvBaseFunc, jvDBaseFunc
    real(kind=8), intent(in) :: nodeCoor(ndim*nno), nodeDisp(ndim*nno)
    real(kind=8), intent(in) :: time, nharm
    integer(kind=8), intent(in) :: strainType
    aster_logical, intent(in) :: lStrainMeca
    type(Material_Para), intent(inout) :: materPara
    real(kind=8), intent(out) :: epsi(nbEpsi*npg)
!
! --------------------------------------------------------------------------------------------------
!
! Compute mechanical strains or total strains
!
! Mechanical strains = total strains - command variables strains
!
! --------------------------------------------------------------------------------------------------
!
! In  nno              : number of nodes
! In  ndim             : dimension of space
! In  nbEpsi           : number of strain tensor components
! In  npg              : number of Gauss points
! In  jvGaussWeight    : adresse to Gauss point weight
! In  jvBaseFunc       : adresse to shape functions
! In  jvDBaseFunc      : adresse to derivative of shape functions
! In  nodeCoor         : coordinates of nodes
! In  nodeDisp         : displacements of nodes
! In  time             : given time
! In  nharm            : Fourier mode
! In  strainType       : type of strain (small, Green, log, etc.)
! In  lStrainMeca      : flag to compute mechanical strains
! IO  materPara        : parameters of material
! Out epsi             : mechanical strains or total strains
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: ksp = 1
    aster_logical :: l_modi_cp
    real(kind=8) :: epsiVarc(162), epsiMeca(162), epsiTota(162)
    real(kind=8) :: epsiLine(162), epsiNlin(162)
    real(kind=8) :: d(4, 4)
    integer(kind=8) :: kpg
    aster_logical :: lStrainVarc
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nbEpsi .le. 6)
    ASSERT(npg .le. 27)
    epsiVarc = 0.d0
    epsiMeca = 0.d0
    epsiTota = 0.d0
    epsiLine = 0.d0
    epsiNlin = 0.d0

! - When compute strains from external state variables ?
    lStrainVarc = lStrainMeca .or. lteatt('C_PLAN', 'OUI')

! - Compute total strain
    if (strainType .eq. STRAIN_TYPE_LOG) then
        call epslmc(nno, ndim, nbEpsi, &
                    npg, jvGaussWeight, jvBaseFunc, &
                    jvDBaseFunc, nodeCoor, nodeDisp, &
                    epsiTota)

    elseif (strainType .eq. STRAIN_TYPE_SMALL) then
        call eps1mc(nno, ndim, nbEpsi, npg, &
                    jvGaussWeight, jvBaseFunc, jvDBaseFunc, &
                    nodeCoor, nodeDisp, nharm, &
                    epsiTota)

    elseif (strainType .eq. STRAIN_TYPE_GREEN) then
        call eps1mc(nno, ndim, nbEpsi, npg, &
                    jvGaussWeight, jvBaseFunc, jvDBaseFunc, &
                    nodeCoor, nodeDisp, nharm, &
                    epsiLine)
        call eps2mc(nno, ndim, nbEpsi, npg, &
                    jvGaussWeight, jvBaseFunc, jvDBaseFunc, &
                    nodeCoor, nodeDisp, &
                    epsiNLin)
        epsiTota = epsiLine+epsiNLin

    else
        ASSERT(ASTER_FALSE)
    end if

! - Compute anelastic strains from external state variables
    if (lStrainVarc) then
        call epthmc(materPara, time, &
                    nbEpsi, npg, ndim, &
                    VARC_STRAIN_ALL, epsiVarc)
    end if

! - Compute mechanical strains
    if (lStrainMeca) then
        epsiMeca = epsiTota-epsiVarc
    end if

! - Modification for plaste stress hypothesis
    if (lteatt('C_PLAN', 'OUI')) then
        do kpg = 1, npg
            l_modi_cp = ASTER_TRUE

! --------- Initializations of material parameters on current integration point
            call initParaPoin(kpg, ksp, materPara)

! --------- Hooke matrix for iso-parametric elements
            call dmatmc(materPara, '+', time, &
                        nbEpsi, d, l_modi_cp)

! --------- Modification of strains
            if (lStrainMeca) then
                epsiMeca(nbEpsi*(kpg-1)+3) = -1.d0/d(3, 3)* &
                                             (d(3, 1)*epsiMeca(nbEpsi*(kpg-1)+1)+ &
                                              d(3, 2)*epsiMeca(nbEpsi*(kpg-1)+2)+ &
                                              d(3, 4)*epsiMeca(nbEpsi*(kpg-1)+4)*2.d0)
            else
                epsiTota(nbEpsi*(kpg-1)+3) = -1.d0/d(3, 3)* &
                                             (d(3, 1)*(epsiTota(nbEpsi*(kpg-1)+1)- &
                                                       epsiVarc(nbEpsi*(kpg-1)+1))+ &
                                              d(3, 2)*(epsiTota(nbEpsi*(kpg-1)+2)- &
                                                       epsiVarc(nbEpsi*(kpg-1)+2))+ &
                                              d(3, 4)*(epsiTota(nbEpsi*(kpg-1)+4)- &
                                                       epsiVarc(nbEpsi*(kpg-1)+4))*2.d0)+ &
                                             epsiVarc(nbEpsi*(kpg-1)+3)
            end if
        end do
    end if

! - Select output strains
    epsi = 0.d0
    if (lStrainMeca) then
        epsi(1:nbEpsi*npg) = epsiMeca(1:nbEpsi*npg)
    else
        epsi(1:nbEpsi*npg) = epsiTota(1:nbEpsi*npg)
    end if
    if (lteatt('D_PLAN', 'OUI')) then
        do kpg = 1, npg
            epsi(nbEpsi*(kpg-1)+3) = 0.d0
        end do
    end if
!
end subroutine
