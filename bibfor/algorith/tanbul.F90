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
subroutine tanbul(materPara, relaComp, &
                  ndim, mini, &
                  alpha, dsbdep, &
                  lVect_, traceEpsiTher_)
!
    use MaterialPara_type
    use BehaviourStrain_type
    implicit none
!
#include "asterc/r8miem.h"
#include "asterc/r8vide.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/ElasticityMaterial_type.h"
#include "asterfort/epstmc.h"
#include "asterfort/get_elas_para.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    type(Material_Para), intent(inout) :: materPara
    character(len=16), intent(in) :: relaComp
    integer(kind=8), intent(in) :: ndim
    aster_logical, intent(in) :: mini
    real(kind=8), intent(out) :: alpha, dsbdep(2*ndim, 2*ndim)
    aster_logical, optional, intent(in) :: lVect_
    real(kind=8), optional, intent(out) :: traceEpsiTher_
!
! --------------------------------------------------------------------------------------------------
!
!          CALCUL DE LA MATRICE TANGENTE BULLE
!
! --------------------------------------------------------------------------------------------------
!
! IN  RESI    : CALCUL DES FORCES INTERNES
! IN  MINI    : STABILISATION BULLE - MINI ELEMENT
! IN  OPTION  : OPTION DE CALCUL
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  G       : NUMERO DU POINT DE GAUSS
! IN  MATE    : NUMERO DU MATERIAU
! OUT ALPHA   : INVERSE DE KAPPA
! OUT DSBDEP  : MATRICE TANGENTE BULLE
! OUT TREPST  : TRACE DU TENSEUR DES DEFORMATIONS THERMIQUES
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: k, jvTime, iret
    real(kind=8) :: e, nu, time
    real(kind=8) :: epsiTher(6)
    real(kind=8) :: coef, coef1, coef2, coef3
    aster_logical :: hasEpsiTher, lVect
    type(All_Varc_Strain) :: allVarcStrain
!
! --------------------------------------------------------------------------------------------------
!
    lVect = ASTER_FALSE
    if (present(lVect_)) then
        lVect = lVect_
    end if
    dsbdep = 0.d0
    alpha = 0.d0
    if (.not. (relaComp(1:4) .eq. 'ELAS' .or. relaComp(1:9) .eq. 'VMIS_ISOT')) then
        call utmess('F', 'ELEMENTSINCO_1')
    end if

! - Get current time
    call tecach('NNO', 'PINSTR', 'L', iret, iad=jvTime)
    if (jvTime .ne. 0) then
        time = zr(jvTime)
    else
        time = r8vide()
    end if
    allVarcStrain%time = time

! - Get elastic parameters
    if (materPara%elasID .ne. ELAS_ISOT) then
        call utmess("F", "ELEMENTSINCO_2")
    end if
    if (allVarcStrain%hasTime) then
        call get_elas_para(materPara%schemePara%fami, materPara%jvMaterCode, '+', &
                           materPara%schemePara%kpg, materPara%schemePara%ksp, &
                           materPara%elasID, materPara%elasKeyword, &
                           time=time, &
                           e_=e, nu_=nu)
    else
        call get_elas_para(materPara%schemePara%fami, materPara%jvMaterCode, '+', &
                           materPara%schemePara%kpg, materPara%schemePara%ksp, &
                           materPara%elasID, materPara%elasKeyword, &
                           e_=e, nu_=nu)
    end if
    alpha = (3.d0*(1.d0-2.d0*nu))/e

! - Compute matrix
    if (mini) then
        coef = 1.d0/((1.d0+nu)*(1.d0-2.d0*nu))
        coef1 = e*(1.d0-nu)*coef
        coef2 = e*nu*coef
        coef3 = 2.d0*e/(1.d0+nu)
        dsbdep(1, 1) = coef1
        dsbdep(1, 2) = coef2
        dsbdep(1, 3) = coef2
        dsbdep(2, 1) = coef2
        dsbdep(2, 2) = coef1
        dsbdep(2, 3) = coef2
        dsbdep(3, 1) = coef2
        dsbdep(3, 2) = coef2
        dsbdep(3, 3) = coef1
        dsbdep(4, 4) = coef3
        if (ndim .eq. 3) then
            dsbdep(5, 5) = coef3
            dsbdep(6, 6) = coef3
        end if
    end if

! - Compute residual
    if (lVect) then
        time = r8vide()
        epsiTher = 0.d0
        ASSERT(present(traceEpsiTher_))
        traceEpsiTher_ = 0.d0

! ----- Compute thermal strains
        call epstmc(materPara, '+', time, ndim, &
                    VARC_STRAIN_TEMP, allVarcStrain, &
                    epsiTher)

! ----- Compute trace
        hasEpsiTher = ASTER_FALSE
        do k = 1, 6
            if (abs(epsiTher(k)) .gt. r8miem()) then
                hasEpsiTher = ASTER_TRUE
            end if
        end do
        if (hasEpsiTher) then
            do k = 1, 3
                traceEpsiTher_ = traceEpsiTher_+epsiTher(k)
            end do
        end if
    end if
!
end subroutine
