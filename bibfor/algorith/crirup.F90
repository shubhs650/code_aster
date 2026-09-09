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
subroutine crirup(materPara, &
                  ndim, npg, lgpg, &
                  option, sigp, vip, vim, &
                  instam, instap)
!
    use MaterialPara_type
    implicit none
!
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/fgequi.h"
#include "asterfort/rcvalb.h"
!
    type(Material_Para), intent(in) :: materPara
    character(len=16), intent(in) :: option
!
! --------------------------------------------------------------------------------------------------
!
! SOUS ROUTINE PERMETTANT L ELVALUATION DU CRITERE
! METHODE MISE EN OEUVRE : MOYENNE SUR LES POINTS DE GAUSS PUIS
! EVALUATION DE LA CONTRAINTE PRINCIPALE MAXIMALE
!
! --------------------------------------------------------------------------------------------------
!
!    NDIM   :  DIMENSION DE L'ESPACE
!    NPG    :  NOMBRE DE POINTS DE GAUSS
!    LGPG   :  NOMBRE TOTAL DE VI APRES AJOUT DES VI RUPTURE
!   OPTION  : OPTION DEMANDEE : RIGI_MECA_TANG ,FULL_MECA ,RAPH_MECA
!   COMPOR  : CARTE DECRIVANT LES PARAMETRES K16 DU COMPORTEMENT
!     VIM   : VARIABLES INTERNES A L'INSTANT DU CALCUL PRECEDENT
!    INSTAM : INSTANT PRECEDENT
! OUT  SIGP : CONTRAINTES DE CAUCHY (RAPH_MECA ET FULL_MECA)
!     VIP   : VARIABLES INTERNES    (RAPH_MECA ET FULL_MECA)
!    INSTAP : INSTANT DE CALCUL
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: npg, kpg, i, ndim, lgpg, ivp, cerr(1)
    real(kind=8) :: sigmoy(6), sigp(2*ndim, npg), equi(20), sc(1), prin1
    real(kind=8) :: vip(lgpg, npg), vim(lgpg, npg), pm, pp, dp, instam, instap
    real(kind=8) :: dt
!
! --------------------------------------------------------------------------------------------------
!
    if (option(1:9) .eq. 'FULL_MECA' .or. option(1:9) .eq. 'RAPH_MECA') then

! ----- Get parameters
        call rcvalb(materPara%schemePara%fami, &
                    materPara%schemePara%kpg, &
                    materPara%schemePara%ksp, &
                    '+', materPara%jvMaterCode, &
                    ' ', 'CRIT_RUPT', &
                    0, ' ', [0.d0], &
                    1, 'SIGM_C', sc, cerr, 1)

!       CALCUL DU TENSEUR MOYEN
        sigmoy = 0.d0
        do i = 1, 2*ndim
            do kpg = 1, npg
                sigmoy(i) = sigmoy(i)+sigp(i, kpg)/npg
            end do
        end do

!       EVALUATION DE LA CONTRAINTE PRINCIPALE MAXIMALE
        call fgequi(sigmoy, 'SIGM', ndim, equi)
        prin1 = max(equi(3), equi(4))
        prin1 = max(equi(5), prin1)

!       CALCUL DE P MOYEN
        pp = 0
        pm = 0
        if (ndim .eq. 2) then
            ivp = 9
        elseif (ndim .eq. 3) then
            ivp = 13
        else
            ASSERT(ASTER_FALSE)
        end if
        do kpg = 1, npg
            pm = pm+vim(ivp, kpg)/npg
            pp = pp+vip(ivp, kpg)/npg
        end do
        dp = pp-pm
        dt = instap-instam
!
!     EVALUATION DE LA VITESSE DE DEFORMATION PLASTIQUE : DP/DT
!     ET DE DIFFERENTES ENERGIES :
!     V(LGPG-4) : ENERGIE DISSIPEE, DP*SIGMOY_EG
!     V(LGPG-3) : ENREGIE DISSIPEE CUMULEE A CHAQUE PAS,
!     V(LGPG-2) : PUISSANCE DISSIPEE, DP/DT*SIGMOY_EG
!     V(LGPG-1) : PUISSANCE DISSIPEE CUMULEE A CHAQUE PAS,
        do kpg = 1, npg
            vip(lgpg-5, kpg) = dp/dt
            vip(lgpg-4, kpg) = dp*equi(1)
            vip(lgpg-3, kpg) = dp*equi(1)+vim(lgpg-3, kpg)
            vip(lgpg-2, kpg) = dp*equi(1)/dt
            vip(lgpg-1, kpg) = dp*equi(1)/dt+vim(lgpg-1, kpg)
        end do
!
!     CRITERE DE RUPTURE
        if (prin1 .gt. sc(1)) then
! LA CONTRAINTE PRINCIPALE MAXI DEPASSE LE SEUIL
            do kpg = 1, npg
                vip(lgpg, kpg) = 1.d0
            end do
        else if (abs(vim(lgpg, 1)-1.d0) .lt. r8prem()) then
! LA MAILLE ETAIT DEJA CASSEE. ELLE LE RESTE
            do kpg = 1, npg
                vip(lgpg, kpg) = 1.d0
            end do
        else
! MAILLE SAINE
            do kpg = 1, npg
                vip(lgpg, kpg) = 0.d0
            end do
        end if
!
    end if
end subroutine
