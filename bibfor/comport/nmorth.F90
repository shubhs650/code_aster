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
subroutine nmorth(materPara, ndim, poum, &
                  dEpsiIn, sigmPrev, option, &
                  sigmCurr, dsidep)
!
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterc/r8vide.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/d1ma3d.h"
#include "asterfort/d1mamc.h"
#include "asterfort/dmat3d.h"
#include "asterfort/dmatmc.h"
#include "asterfort/ElasticityMaterial_type.h"
#include "asterfort/lteatt.h"
#include "asterfort/matrot.h"
#include "asterfort/utbtab.h"
#include "asterfort/utmess.h"
#include "asterfort/verifepsa.h"
#include "asterfort/verifh.h"
#include "asterfort/verifs.h"
#include "asterfort/verift.h"
!
    type(Material_Para), intent(in) :: materPara
    integer(kind=8), intent(in) :: ndim
    character(len=*), intent(in) :: poum
    real(kind=8), intent(in) :: dEpsiIn(2*ndim)
    real(kind=8), intent(in) :: sigmPrev(2*ndim)
    character(len=16), intent(in):: option
    real(kind=8), intent(out) :: sigmCurr(2*ndim)
    real(kind=8), intent(out) :: dsidep(2*ndim, 2*ndim)
!
! --------------------------------------------------------------------------------------------------
!
!  IN    FAMI   : FAMILLE DE POINT DE GAUSS
!  IN    KPG    : NUMERO DU POINT DE GAUSS
!  IN    KSP    : NUMERO DU SOUS POINT DE GAUSS
!  IN    NDIM   : DIMENSION DU PROBLEME
!  IN    PHENOM : PHENOMENE (ELAS_ORTH OU ELAS_ISTR)
!  IN    TYPMOD : TYPE DE MODELISATION
!  IN    IMATE  : ADRESSE DU MATERIAU
!  IN    POUM   : '+' INSTANT SUIVANT OU '-' INSTANT COURANT
!                 OU 'T' = '+' - '-' INCREMENT
!  IN    EPSM   : DEFORMATION A L INSTANT T-
!  IN    DESPS  : INCREMENT DE DEFORMATION
!  IN    SIGM   : CONTRAINTE A L INSTANT T-
!  IN    OPTION : OPTION A CALCULER
!  IN    ANGL_NAUT : ANGLE DU REPERE LOCAL D ORTHOTROPIE
!  OUT   SIGP   : CONTRAINTE A L INSTANT T+
!                 CAR IL EN EXISTE FORCEMENT UNE)
!  OUT   DSIDEP : MATRICE DE RIGIDITE TANGENTE
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), parameter :: rac2 = sqrt(2.d0)
    real(kind=8) :: p(3, 3)
    real(kind=8) :: timeNaN, hooke(36), ekooh(36)
    real(kind=8) :: dEpsi(6)
    real(kind=8) :: epsiTherAnis(3), dEpsiTherLoca(6), dEpsiTherGlob(6)
    real(kind=8) :: dEpsiHydr, dEpsiSech, dEpsiAnel(6)
    real(kind=8) :: dEpsiMeca(6), epsm2(6)
    real(kind=8) :: work(3, 3), deplth_mat(3, 3), depgth_mat(3, 3)
    real(kind=8) :: sigm2(2*ndim)
    integer(kind=8) :: ndimsi, i, j
    aster_logical :: lLegitModel
!
! --------------------------------------------------------------------------------------------------
!
    timeNaN = r8vide()
!
    if (materPara%elasID .eq. ELAS_ISTR .and. ndim .eq. 2) then
        call utmess('F', 'ELEMENTS3_2')
    end if

    ndimsi = ndim*2
    dsidep = 0
    dEpsiTherGlob = 0
!
    if (option .eq. 'FULL_MECA' .or. option .eq. 'RAPH_MECA') then
        do i = 1, ndimsi
            if (i .le. 3) then
                dEpsi(i) = dEpsiIn(i)
            else
                dEpsi(i) = dEpsiIn(i)*rac2
            end if
        end do
    end if
!
    call chckLCSValid(materPara)

! - Check legit model
    lLegitModel = ASTER_FALSE
    if (materPara%schemePara%fami .eq. 'PMAT') then
        lLegitModel = ASTER_TRUE
    else
        if (lteatt('DIM_TOPO_MAILLE', '3')) then
            lLegitModel = ASTER_TRUE
        else if (lteatt('C_PLAN', 'OUI')) then
            lLegitModel = ASTER_TRUE
        else if (lteatt('D_PLAN', 'OUI')) then
            lLegitModel = ASTER_TRUE
        else if (lteatt('AXIS', 'OUI')) then
            lLegitModel = ASTER_TRUE
        end if
    end if
    if (.not. lLegitModel) then
        call utmess('F', 'ALGORITH8_22')
    end if

! - MATRICES TANGENTES
    if (materPara%schemePara%fami .eq. 'PMAT') then
        if (option .eq. 'RIGI_MECA_TANG') then
            call dmat3d(materPara, '-', timeNaN, hooke)
        else
            call d1ma3d(materPara, '-', timeNaN, ekooh)
            call dmat3d(materPara, '+', timeNaN, hooke)
        end if
    else
        if (option .eq. 'RIGI_MECA_TANG') then
            call dmatmc(materPara, "-", timeNaN, ndimsi, hooke)
        else
            call d1mamc(materPara, '-', timeNaN, ndimsi, ekooh)
            call dmatmc(materPara, '+', timeNaN, ndimsi, hooke)
        end if
    end if
!
    if (option(1:10) .eq. 'RIGI_MECA_' .or. option(1:9) .eq. 'FULL_MECA') then
        do i = 1, ndimsi
            do j = 1, ndimsi
                dsidep(i, j) = hooke(ndimsi*(j-1)+i)
            end do
        end do
    end if

! - INTEGRATION
    if (option .eq. 'FULL_MECA' .or. option .eq. 'RAPH_MECA') then
! ----- Thermal strains (local)
        dEpsiTherLoca = 0.d0
        if (materPara%elasID .eq. ELAS_ORTH) then
            call verift(materPara%schemePara%fami, &
                        materPara%schemePara%kpg, materPara%schemePara%ksp, &
                        poum, materPara%jvMaterCode, &
                        epsth_anis_=epsiTherAnis)
            dEpsiTherLoca(1) = epsiTherAnis(1)
            dEpsiTherLoca(2) = epsiTherAnis(2)
            dEpsiTherLoca(3) = epsiTherAnis(3)

        else if (materPara%elasID .eq. ELAS_ISTR) then
            call verift(materPara%schemePara%fami, &
                        materPara%schemePara%kpg, materPara%schemePara%ksp, &
                        poum, materPara%jvMaterCode, &
                        epsth_anis_=epsiTherAnis)
            dEpsiTherLoca(1) = epsiTherAnis(1)
            dEpsiTherLoca(2) = epsiTherAnis(1)
            dEpsiTherLoca(3) = epsiTherAnis(2)

        else
            ASSERT(ASTER_FALSE)
        end if

!       RECUPERATION DE LA MATRICE DE PASSAGE
        if (chckLCSDefine(materPara%lcsPara)) then
            call matrot(materPara%lcsPara%lcsAngle, p)
        else
            call utmess("F", "ALGORITH8_20")
        end if

! ----- Thermal strains (global)
        do i = 1, 3
            deplth_mat(i, i) = dEpsiTherLoca(i)
        end do
        deplth_mat(1, 2) = dEpsiTherLoca(4)
        deplth_mat(1, 3) = dEpsiTherLoca(5)
        deplth_mat(2, 3) = dEpsiTherLoca(6)
        deplth_mat(2, 1) = deplth_mat(1, 2)
        deplth_mat(3, 1) = deplth_mat(1, 3)
        deplth_mat(3, 2) = deplth_mat(2, 3)
        call utbtab('ZERO', 3, 3, deplth_mat, p, work, depgth_mat)
        dEpsiTherGlob(1) = depgth_mat(1, 1)
        dEpsiTherGlob(2) = depgth_mat(2, 2)
        dEpsiTherGlob(3) = depgth_mat(3, 3)
        dEpsiTherGlob(4) = depgth_mat(1, 2)
        dEpsiTherGlob(5) = depgth_mat(1, 3)
        dEpsiTherGlob(6) = depgth_mat(2, 3)

! ----- Get increment of external state variables
        call verifh(materPara%schemePara%fami, &
                    materPara%schemePara%kpg, materPara%schemePara%ksp, &
                    poum, materPara%jvMaterCode, dEpsiHydr)
        call verifs(materPara%schemePara%fami, &
                    materPara%schemePara%kpg, materPara%schemePara%ksp, &
                    poum, materPara%jvMaterCode, dEpsiSech)
        call verifepsa(materPara%schemePara%fami, &
                       materPara%schemePara%kpg, materPara%schemePara%ksp, poum, dEpsiAnel)

        do i = 1, ndimsi
            if (i .le. 3) then
                dEpsiMeca(i) = dEpsi(i)-dEpsiTherGlob(i)-dEpsiHydr-dEpsiSech-dEpsiAnel(i)
            else
                dEpsiMeca(i) = dEpsi(i)-2.0*dEpsiTherGlob(i)-2.0*dEpsiAnel(i)
            end if
        end do

        sigm2 = sigmPrev
        do i = 4, ndimsi
            sigm2(i) = sigm2(i)/rac2
        end do

! MODIFICATION DE SIGM POUR PRENDRE EN COMPTE LA VARIATION DE
! COEF ELASTIQUES AVEC LA TEMPERATURE
!
        do i = 1, ndimsi
            epsm2(i) = 0.d0
            do j = 1, ndimsi
                epsm2(i) = epsm2(i)+ekooh(ndimsi*(j-1)+i)*sigm2(j)
            end do
        end do
!
        do i = 1, ndimsi
            sigmCurr(i) = 0.d0
            do j = 1, ndimsi
                sigmCurr(i) = sigmCurr(i)+ &
                              hooke(ndimsi*(j-1)+i)*(dEpsiMeca(j)+epsm2(j))
            end do
        end do
        do i = 4, ndimsi
            sigmCurr(i) = sigmCurr(i)*rac2
        end do
    end if
!
    if (option(1:10) .eq. 'RIGI_MECA_' .or. option(1:9) .eq. 'FULL_MECA') then
        do i = 1, ndimsi
            do j = 4, ndimsi
                dsidep(i, j) = dsidep(i, j)*sqrt(2.d0)
            end do
        end do
        do i = 4, ndimsi
            do j = 1, ndimsi
                dsidep(i, j) = dsidep(i, j)*sqrt(2.d0)
            end do
        end do
    end if
!
end subroutine
