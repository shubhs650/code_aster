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
subroutine te0544_implement(typmod, compor, ndim, nno, npg, &
                            jv_poids, jv_vff, jv_dfde, geom_i, &
                            deplm, ddepl, depl0, depl1, dtau, copilo)

    use tenseur_dime_module, only: NDIM_TO_NDIMSI, matrix_to_voigt, identity
    implicit none

#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterc/matfpe.h"
#include "asterc/r8vide.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/nmgeom.h"
#include "asterfort/pidefo.h"
#include "asterfort/pipdef.h"
!
    character(len=8), intent(in):: typmod(:)
    character(len=16), intent(in) :: compor(COMPOR_SIZE)
    integer(kind=8) :: ndim, nno, npg
    integer(kind=8) :: jv_poids, jv_vff, jv_dfde
    real(kind=8) :: geom_i(ndim, nno), deplm(ndim, nno), ddepl(ndim, nno)
    real(kind=8) :: depl0(ndim, nno), depl1(ndim, nno)
    real(kind=8) :: dtau
    real(kind=8) :: copilo(5, npg)
!
! --------------------------------------------------------------------------------------------------
!
! routine meca_non_line (pilotage)
!
! calcul des coefficients de pilotage pour pred_elas/deformation
!
! --------------------------------------------------------------------------------------------------
! in  typmod : type de modélisation (pour distinguer 1D, AXIS et 2D/3D)
! in  compor : carte comportement (pour grandes déformations ou pas)
! in  ndim   : dimension de l'espace
! in  nno    : nombre de noeuds de l'element
! in  npg    : nombre de points de gauss
! in  jv_poids : poids des points de gauss
! in  jv_vff    : valeur des fonctions de forme
! in  jv_dfde  : derivee des fonctions de forme element de reference
! in  geom   : coordonees des noeuds
! in  deplm  : deplacement en t-
! in  ddepl  : increment de deplacement a l'iteration newton courante
! in  depl0  : correction de deplacement pour forces fixes
! in  depl1  : correction de deplacement pour forces pilotees
! out copilo : coefficients a0 et a1 pour chaque point de gauss
! --------------------------------------------------------------------------------------------------
    aster_logical:: axi, grand
    integer(kind=8) :: kpg, ndimsi
    real(kind=8):: epsm(NDIM_TO_NDIMSI(ndim))
    real(kind=8):: epsd_cste(NDIM_TO_NDIMSI(ndim))
    real(kind=8):: epsd_pilo(NDIM_TO_NDIMSI(ndim))
    real(kind=8):: geom_m(ndim, nno), dfdi(ndim, nno), poids, r
    real(kind=8):: fm(3, 3), finger(3, 3)
! --------------------------------------------------------------------------------------------------

    call matfpe(-1)

! - initialisations
    ndimsi = NDIM_TO_NDIMSI(ndim)
    copilo = r8vide()
    axi = typmod(1) .eq. 'AXIS'
    grand = compor(DEFO) .ne. 'PETIT'

    ! Géométrie actualisée en grandes déformation
    geom_m = merge(geom_i+deplm, geom_i, grand)

! - traitement de chaque point de gauss
    do kpg = 1, npg

        ! Deformations (linéarisées autour de la configuration u- si grandes déf)
        call pipdef(typmod, &
                    ndim, nno, kpg, jv_poids, jv_vff, &
                    jv_dfde, geom_m, deplm, &
                    ddepl, depl0, depl1, &
                    epsm, epsd_cste, epsd_pilo)

        ! En grandes déformations, on prend pour epsm la déformation de Finger en U-
        if (grand) then

            call nmgeom(ndim, nno, axi, ASTER_TRUE, geom_i, &
                        kpg, jv_poids, jv_vff, jv_dfde, deplm, &
                        ASTER_TRUE, poids, dfdi, fm, epsm, r)

            finger = 0.5d0*(matmul(fm, transpose(fm))-identity(3))
            epsm = matrix_to_voigt(finger, ndimsi)

        end if

        ! Pilotage en déformation
        call pidefo(epsm, epsd_cste, epsd_pilo, dtau, copilo(:, kpg))

    end do
    call matfpe(1)

end subroutine
