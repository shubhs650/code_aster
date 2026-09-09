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

subroutine nmasym(materPara, option, &
                  xlong0, a, dlong0, &
                  effnom, vim, effnop, vip, klv, &
                  fono)
!
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterfort/nm1das.h"
#include "asterfort/r8inir.h"
#include "asterfort/rcvalb.h"
!
    integer(kind=8), parameter :: neq = 6, nbt = 21, nvar = 4
    type(Material_Para), intent(in) :: materPara
    character(len=*) :: option
    real(kind=8) :: xlong0, a, syc, syt, etc, ett
    real(kind=8) :: e, dlong0
    real(kind=8) :: effnom, vim(nvar)
    real(kind=8) :: effnop, vip(nvar), fono(neq), klv(nbt)
!
! --------------------------------------------------------------------------------------------------
!
!    TRAITEMENT DE LA RELATION DE COMPORTEMENT -ELASTOPLASTICITE-
!    ECROUISSAGE ISOTROPE ASYMETRIQUE LINEAIRE - VON MISES-
!    POUR UN MODELE BARRE ELEMENT MECA_BARRE
!
! --------------------------------------------------------------------------------------------------
!
!       XLONG0 : LONGUEUR DE L'ELEMENT DE BARRE AU REPOS
!       A      : SECTION DE LA BARRE
!       XLONGM : LONGEUR DE L'ELEMENT AU TEMPS MOINS
!       DLONG0 : INCREMENT D'ALLONGEMENT DE L'ELEMENT
!       EFFNOM : EFFORT NORMAL PRECEDENT
!       OPTION : OPTION DEMANDEE (R_M_T,FULL OU RAPH_MECA)
!
! OUT : EFFNOP : CONTRAINTE A L'INSTANT ACTUEL
!       VIP    : VARIABLE INTERNE A L'INSTANT ACTUEL
!       FONO   : FORCES NODALES COURANTES
!       KLV    : MATRICE TANGENTE
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: kpgFPG1 = 1, kspFPG1 = 1
    character(len=8), parameter :: famiFPG1 = "FPG1"
    type(Material_Para) :: materParaFPG1
    real(kind=8) :: sigm, deps, dsdem, dsdep, sigp, xrig
    character(len=16), parameter :: propElas = "E"
    integer(kind=8), parameter :: nbProp = 4
    integer(kind=8) :: propCode(nbProp)
    real(kind=8) :: propVale(nbProp)
    character(len=16), parameter :: propName(nbProp) = &
                                    (/'SY_C        ', 'DC_SIGM_EPSI', &
                                      'SY_T        ', 'DT_SIGM_EPSI'/)
!
! --------------------------------------------------------------------------------------------------
!
    call r8inir(nbt, 0.d0, klv, 1)
    call r8inir(neq, 0.d0, fono, 1)
!
!----------RECUPERATION DES CARACTERISTIQUES
!
    deps = dlong0/xlong0
    sigm = effnom/a

! - Copy material parameters with other scheme parameters
    call copyMaterPara(materPara, famiFPG1, kpgFPG1, kspFPG1, &
                       materParaFPG1)

! - CARACTERISTIQUES ELASTIQUES
    call rcvalb(materParaFPG1%schemePara%fami, &
                materParaFPG1%schemePara%kpg, &
                materParaFPG1%schemePara%ksp, &
                '+', &
                materParaFPG1%jvMaterCode, &
                ' ', 'ELAS', &
                0, ' ', [0.d0], &
                1, propElas, propVale, &
                propCode, 1)
    e = propVale(1)

! - CARACTERISTIQUES ECROUISSAGE LINEAIRE ASYMETRIQUE
    call rcvalb(materParaFPG1%schemePara%fami, &
                materParaFPG1%schemePara%kpg, &
                materParaFPG1%schemePara%ksp, &
                '+', &
                materParaFPG1%jvMaterCode, &
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
    effnop = sigp*a

! - CALCUL DU COEFFICIENT NON NUL DE LA MATRICE TANGENTE
    if (option(1:10) .eq. 'RIGI_MECA_' .or. option(1:9) .eq. 'FULL_MECA') then
!
        if (option(11:14) .eq. 'ELAS') then
            xrig = e*a/xlong0
        else
            if (option(1:14) .eq. 'RIGI_MECA_TANG') then
                xrig = dsdem*a/xlong0
            else
                xrig = dsdep*a/xlong0
            end if
        end if
        klv(1) = xrig
        klv(7) = -xrig
        klv(10) = xrig
    end if

! --- CALCUL DES FORCES NODALES
    fono(1) = -effnop
    fono(4) = effnop
!
end subroutine
