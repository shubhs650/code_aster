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
subroutine mat_proto(BEHInteg, &
                     fami, kpg, ksp, poum, jvMaterCode, relaComp, &
                     nprops, props)
!
    use Behaviour_type
    implicit none
!
#include "asterc/r8nnem.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/metaGetPhase.h"
#include "asterfort/metaGetType.h"
#include "asterfort/rcadlv.h"
#include "jeveux.h"
!
    type(Behaviour_Integ), intent(in) :: BEHInteg
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: kpg, ksp
    character(len=1), intent(in) :: poum
    integer(kind=8), intent(in) :: jvMaterCode
    character(len=*), intent(in) :: relaComp
    integer(kind=8), intent(inout) :: nprops
    real(kind=8), intent(out) :: props(*)
!
! --------------------------------------------------------------------------------------------------
!
! Behaviour (MFront/UMAT)
!
! Calculer les coef. materiau pour l'interface umat et mfront
!
! --------------------------------------------------------------------------------------------------
!
! In  BEHInteg         : parameters for integration of behaviour
!       in   fami      : famille de point de gauss (rigi,mass,...)
!       in   kpg,ksp   : numero du (sous)point de gauss
!       in   poum      : '+' /'-'
!       in   nprops    : en entree : dimension du tableau props
!       in   relaComp    : nom de l'interface de prototypage
!       out  nprops    : en sortie : nombre de coefficients recuperes dans props
!            props(*)  : coeficients du materiau
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i, jadr, icodre, ncoef, nbPara
    real(kind=8) :: phase(5), zalpha
    integer(kind=8) :: metaType, nbPhases
    integer(kind=8), parameter :: nbParaMaxi = 4
    real(kind=8) :: paraVale(nbParaMaxi)
    character(len=16) :: paraName(nbParaMaxi)
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(relaComp .eq. 'UMAT' .or. relaComp .eq. 'MFRONT')
    ASSERT(nprops .gt. 0 .and. nprops .le. 197)
    props(1:nprops) = r8nnem()

! - Get material properties
    zalpha = 0.d0
    nbPara = 0
    paraName = " "
    paraVale = 0.d0
    if (BEHInteg%materPara%lElasIsMeta) then
        call metaGetType(metaType, nbPhases)
        if (nbPhases .ne. 0) then
            call metaGetPhase(fami, poum, kpg, ksp, metaType, &
                              nbPhases, phase, zcold_=zalpha)
            nbPara = nbPara+1
            paraName(nbPara) = "META"
            paraVale(nbPara) = zalpha
        end if
    end if

    if (.not. BEHInteg%behavESVA%lGeomInESVA) then
        ASSERT(fami .ne. 'XFEM')
        nbPara = nbPara+1
        paraName(nbPara) = "X"
        paraVale(nbPara) = BEHInteg%behavESVA%behavESVAGeom%coorElga(kpg, 1)
        nbPara = nbPara+1
        paraName(nbPara) = "Y"
        paraVale(nbPara) = BEHInteg%behavESVA%behavESVAGeom%coorElga(kpg, 2)
        nbPara = nbPara+1
        paraName(nbPara) = "Z"
        paraVale(nbPara) = BEHInteg%behavESVA%behavESVAGeom%coorElga(kpg, 3)
    end if
    ASSERT(nbPara .le. nbParaMaxi)
    call rcadlv(fami, kpg, ksp, poum, jvMaterCode, ' ', relaComp, &
                'LISTE_COEF', nbPara, paraName, paraVale, jadr, ncoef, icodre, 1)
    ASSERT(icodre .eq. 0)

! - Copy properties
    ASSERT(ncoef .le. nprops)
    do i = 1, ncoef
        props(i) = zr(jadr-1+i)
    end do
    nprops = ncoef
!
end subroutine
