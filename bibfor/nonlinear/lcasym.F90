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

subroutine lcasym(materPara, option, sigm, vim, deps, sigp, vip, dsde)
!
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterfort/nm1das.h"
#include "asterfort/rcvalb.h"
!
    type(Material_Para), intent(in) :: materPara
    character(len=*), intent(in) :: option
    real(kind=8), intent(in):: sigm, vim(:), deps
    real(kind=8), intent(out) :: sigp, vip(:), dsde
!
! --------------------------------------------------------------------------------------------------
!
!    TRAITEMENT DE LA RELATION DE COMPORTEMENT -ELASTOPLASTICITE-
!    ECROUISSAGE ISOTROPE ASYMETRIQUE LINEAIRE - VON MISES-
!    POUR UN MODELE BARRE ELEMENT MECA_BARRE
!
! --------------------------------------------------------------------------------------------------
    real(kind=8) :: syc, syt, etc, ett, e
    real(kind=8) :: dsdem, dsdep
    character(len=16), parameter :: propElas = "E"
    integer(kind=8), parameter :: nbProp = 4
    integer(kind=8) :: propCode(nbProp)
    real(kind=8) :: propVale(nbProp)
    character(len=16), parameter :: propName(nbProp) = &
                                    (/'SY_C        ', 'DC_SIGM_EPSI', &
                                      'SY_T        ', 'DT_SIGM_EPSI'/)
! --------------------------------------------------------------------------------------------------

! - CARACTERISTIQUES ELASTIQUES
    call rcvalb(materPara%schemePara%fami, &
                materPara%schemePara%kpg, &
                materPara%schemePara%ksp, &
                '+', &
                materPara%jvMaterCode, &
                ' ', 'ELAS', &
                0, ' ', [0.d0], &
                1, propElas, propVale, &
                propCode, 1)
    e = propVale(1)

! - CARACTERISTIQUES ECROUISSAGE LINEAIRE ASYMETRIQUE
    call rcvalb(materPara%schemePara%fami, &
                materPara%schemePara%kpg, &
                materPara%schemePara%ksp, &
                '+', &
                materPara%jvMaterCode, &
                ' ', 'ECRO_ASYM_LINE', &
                0, ' ', [0.d0], &
                nbProp, propName, propVale, &
                propCode, 1)
    syc = propVale(1)
    etc = propVale(2)
    syt = propVale(3)
    ett = propVale(4)
    call nm1das(materPara, &
                e, syc, &
                syt, etc, ett, &
                sigm, deps, vim, &
                sigp, vip, dsdem, dsdep)

! - CALCUL DU COEFFICIENT NON NUL DE LA MATRICE TANGENTE
    if (option(1:10) .eq. 'RIGI_MECA_' .or. option(1:9) .eq. 'FULL_MECA') then
        if (option(11:14) .eq. 'ELAS') then
            dsde = e
        else
            if (option(1:14) .eq. 'RIGI_MECA_TANG') then
                dsde = dsdem
            else
                dsde = dsdep
            end if
        end if

    end if

end subroutine
