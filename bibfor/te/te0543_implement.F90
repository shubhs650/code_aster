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
subroutine te0543_implement(BEHInteg, &
                            typmod, compor, &
                            ndim, nno, npg, &
                            jv_poids, jv_vff, jv_dfde, geom, &
                            lgpg, deplm, sigm, &
                            vim, ddepl, depl0, depl1, &
                            bornes, dtau, copilo)
!
    use Behaviour_type
    use Behaviour_module
    use tenseur_dime_module, only: voigt
    implicit none
!
#include "asterc/matfpe.h"
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/pi0000.h"
#include "asterfort/pipdef.h"
!
    type(Behaviour_Integ), intent(inout) :: BEHInteg
    character(len=8), intent(in) :: typmod(:)
    character(len=16), intent(in) :: compor(COMPOR_SIZE)
    integer(kind=8) :: ndim, nno, npg
    integer(kind=8) :: jv_poids, jv_vff, jv_dfde
    integer(kind=8) :: lgpg
    real(kind=8) :: geom(ndim, nno), deplm(ndim, nno), ddepl(ndim, nno)
    real(kind=8) :: sigm(BEHInteg%behavPara%ndimsi, npg), vim(lgpg, npg)
    real(kind=8) :: depl0(ndim, nno), depl1(ndim, nno)
    real(kind=8) :: copilo(5, npg)
    real(kind=8) :: bornes(2), dtau
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (PILOTAGE)
!
! CALCUL DES COEFFICIENTS DE PILOTAGE POUR PRED_ELAS/DEFORMATION
!
! --------------------------------------------------------------------------------------------------
! in  ndim   : dimension de l'espace
! in  nno    : nombre de noeuds de l'element
! in  npg    : nombre de points de gauss
! in  jv_poids : poids des points de gauss
! in  jv_vff    : valeur des fonctions de forme
! in  jv_dfde  : derivee des fonctions de forme element de reference
! in  geom   : coordonees des noeuds
! in  lgpg   : "longueur" des variables internes pour 1 point de gauss
!             cette longueur est un majorant du nbre reel de var. int.
! in  deplm  : deplacement en t-
! in  ddepl  : increment de deplacement a l'iteration newton courante
! in  sigm   : contraintes de cauchy en t-
! in  vim    : variables internes en t-
! in  depl0  : correction de deplacement pour forces fixes
! in  depl1  : correction de deplacement pour forces pilotees
! in  bornes : bornes pilotage [etamax, etamin]
! in  dtau   : parametre pilotage
! out copilo : coefficients a0 et a1 pour chaque point de gauss
! --------------------------------------------------------------------------------------------------
    integer(kind=8), parameter :: ksp = 1
! --------------------------------------------------------------------------------------------------
    integer(kind=8) :: kpg, ndimsi
    real(kind=8) :: epsm(BEHInteg%behavPara%ndimsi)
    real(kind=8) :: epsd_cste(BEHInteg%behavPara%ndimsi)
    real(kind=8) :: epsd_pilo(BEHInteg%behavPara%ndimsi)
    real(kind=8) :: sigm_ldc(BEHInteg%behavPara%ndimsi)
! --------------------------------------------------------------------------------------------------
!
    call matfpe(-1)

! - initialisations
    ndimsi = BEHInteg%behavPara%ndimsi
    ASSERT(ndimsi .eq. 4 .or. ndimsi .eq. 6)

    sigm_ldc = 0
    copilo = r8vide()

! - traitement de chaque point de gauss
    do kpg = 1, npg

        ! Set main parameters for behaviour (on point)
        call behaviourSetParaPoin(kpg, ksp, BEHInteg)

        ! Calcul des déformations pour le pilotage
        call pipdef(typmod, &
                    ndim, nno, kpg, jv_poids, jv_vff, &
                    jv_dfde, geom, deplm, &
                    ddepl, depl0, depl1, &
                    epsm, epsd_cste, epsd_pilo)

        ! préparation des contraintes en t-
        sigm_ldc = sigm(1:ndimsi, kpg)*voigt(ndimsi)

        ! Calcul des coefficients de pilotage
        call pi0000(BEHInteg, compor, typmod, ndim, &
                    epsm, epsd_cste, epsd_pilo, &
                    sigm_ldc, vim(:, kpg), dtau, bornes(2), bornes(1), &
                    copilo(:, kpg))

    end do
    call matfpe(1)

end subroutine
