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
!
subroutine lc0054(fami, kpg, ksp, ndim, imate, &
                  compor, carcri, instam, instap, epsm, &
                  deps, sigm, nvi, vim, option, &
                  sigp, vip, typmod, &
                  dsidep, codret)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/nmtevp.h"
!
    integer(kind=8) :: imate, ndim, kpg, ksp, codret
    character(len=16), intent(in) :: compor(COMPOR_SIZE), option
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    real(kind=8) :: instam, instap
    integer(kind=8) :: nvi
    real(kind=8) :: epsm(6), deps(6)
    real(kind=8) :: sigm(6), sigp(6)
    real(kind=8) :: vim(nvi), vip(nvi)
    real(kind=8) :: dsidep(6, 6)
    character(len=8) :: typmod(*)
    character(len=*) :: fami
!
!
!.......................................................................
!
!     BUT: LOI DE COMPORTEMENT A ECROUISSAGE ISOTROPE DEPENDANT DE LA
!          VITESSE DE DEFORMATION
!
!          RELATIONS : 'VMIS_JOHN_COOK'
!
!       IN      FAMI    FAMILLE DE POINT DE GAUSS (RIGI,MASS,...)
!       IN      KPG,KSP NUMERO DU (SOUS)POINT DE GAUSS
!       IN      NDIM    DIMENSION DE L ESPACE (3D=3,2D=2,1D=1)
!               TYPMOD  TYPE DE MODELISATION
!               IMATE    ADRESSE DU MATERIAU CODE
!               COMPOR    COMPORTEMENT DE L ELEMENT
!               INSTAM   INSTANT T
!               INSTAP   INSTANT T+DT
!               EPSM   DEFORMATION TOTALE A T
!               DEPS   INCREMENT DE DEFORMATION TOTALE
!               SIGM    CONTRAINTE A T
!               VIM    VARIABLES INTERNES A T    + INDICATEUR ETAT T
!    ATTENTION  VIM    VARIABLES INTERNES A T MODIFIEES SI REDECOUPAGE
!               OPTION     OPTION DE CALCUL A FAIRE
!               ANGMAS
!       OUT     SIGP    CONTRAINTE A T+DT
!               VIP    VARIABLES INTERNES A T+DT + INDICATEUR ETAT T+DT
!               DSIDEP    MATRICE DE COMPORTEMENT TANGENT A T+DT OU T
!.......................................................................
!               CODRET
    character(len=16) :: relaComp
    ASSERT(nvi .eq. 5)
    relaComp = compor(RELA_NAME)

    call nmtevp(fami, kpg, ksp, ndim, typmod, &
                imate, relaComp, carcri, instam, instap, &
                deps, sigm, vim, option, sigp, &
                vip, dsidep, codret)

!
end subroutine
