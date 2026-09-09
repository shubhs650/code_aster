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
subroutine lc0004(fami, kpg, ksp, ndim, imate, &
                  compor, carcri, instam, instap, &
                  deps, sigm, vim, option, &
                  sigp, vip, typmod, &
                  nvi, dsidep, codret)
!
    implicit none
!
#include "asterfort/nmchab.h"
#include "asterfort/Behaviour_type.h"
!
    integer(kind=8) :: kpg, ksp, ndim, imate, codret, nvi
    character(len=*) :: fami
    character(len=8) :: typmod(*)
    character(len=16), intent(in) :: compor(COMPOR_SIZE)
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    character(len=16) :: option
    real(kind=8) :: instam, instap, deps(6)
    real(kind=8) :: sigm(6), vim(nvi), sigp(6), vip(nvi), dsidep(6, 6)
!
!

!
! ======================================================================
!
!     INTEGRATION LOCALE DES LOIS DE COMPORTEMENT DE CHABOCHE
!          RELATIONS : 'VMIS_CIN1_CHAB' 'VMIS_CIN2_CHAB'
!          RELATIONS : 'VISC_CIN1_CHAB' 'VISC_CIN2_CHAB'
!          RELATIONS : 'VMIS_CIN2_MEMO' 'VISC_CIN2_MEMO'
!
!     ARGUMENTS :
!       IN      FAMI    FAMILLE DE POINT DE GAUSS (RIGI,MASS,...)
!       IN      KPG,KSP NUMERO DU (SOUS)POINT DE GAUSS
!       IN      NDIM    DIMENSION DE L ESPACE (3D=3,2D=2,1D=1)
!       IN      TYPMOD  TYPE DE MODELISATION
!               IMATE   ADRESSE DU MATERIAU CODE
!               COMPOR  COMPORTEMENT DE L ELEMENT
!                       COMPOR(1) = RELATION DE COMPORTEMENT
!                       COMPOR(2) = NB DE VARIABLES INTERNES
!                       COMPOR(3) = TYPE DE DEFORMATION
!               CRIT    CRITERES  LOCAUX, EN PARTICULIER
!                       CRIT(1) = NOMBRE D ITERATIONS MAXI A CONVERGENCE
!                                 (ITER_INTE_MAXI == ITECREL)
!                       CRIT(3) = VALEUR DE LA TOLERANCE DE CONVERGENCE
!                                 (RESI_INTE == RESCREL)
!               INSTAM  INSTANT T
!               INSTAP  INSTANT T+DT
!               DEPS    INCREMENT DE DEFORMATION TOTALE
!               SIGM    CONTRAINTE A T
!               VIM     VARIABLES INTERNES A T    + INDICATEUR ETAT T
!    ATTENTION  VIM     VARIABLES INTERNES A T MODIFIEES SI REDECOUPAGE
!               OPTION     OPTION DE CALCUL A FAIRE
!                             'RIGI_MECA_TANG'> DSIDEP(T)
!                             'FULL_MECA'     > DSIDEP(T+DT) , SIG(T+DT)
!                             'RAPH_MECA'     > SIG(T+DT)
!       OUT     SIGP    CONTRAINTE A T+DT
!               VIP     VARIABLES INTERNES A T+DT + INDICATEUR ETAT T+DT
!               DSIDEP  MATRICE DE COMPORTEMENT TANGENT A T+DT OU T
!               IRET    CODE RETOUR DE  L'INTEGRATION DE LA LDC
!                              IRET=0 => PAS DE PROBLEME
!                              IRET=1 => ABSENCE DE CONVERGENCE
!
!               ATTENTION LES TENSEURS ET MATRICES SONT RANGES DANS
!               L'ORDRE :  XX,YY,ZZ,SQRT(2)*XY,SQRT(2)*XZ,SQRT(2)*YZ
!               -----------------------------------------------------
!
    character(len=16) :: relaComp
!

    relaComp = compor(RELA_NAME)
    call nmchab(fami, kpg, ksp, ndim, &
                imate, relaComp, carcri, instam, instap, &
                deps, sigm, vim, option, sigp, &
                vip, dsidep, codret)
end subroutine
