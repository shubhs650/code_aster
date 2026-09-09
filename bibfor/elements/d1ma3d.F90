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
subroutine d1ma3d(materPara, poum, time, d1)
!
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/d1pa3d.h"
#include "asterfort/ElasticityMaterial_type.h"
#include "asterfort/MaterialPara_type.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utbtab.h"
#include "asterfort/utmess.h"
!
    type(Material_Para), intent(in) :: materPara
    character(len=*), intent(in) :: poum
    real(kind=8), intent(in) :: time
    real(kind=8), intent(out) :: d1(6, 6)
!
! --------------------------------------------------------------------------------------------------
!
!     D1MA3D  --   CALCUL DE L'INVERSE DE LA MATRICE DE HOOKE
!                  POUR LES ELEMENTS MASSIFS EN 3D OU EN SERIE DE
!                  FOURIER POUR DES MATERIAUX ISOTROPE, ORTHOTROPE
!                  ET ISOTROPE TRANSVERSE
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), parameter :: zero = 0.d0, un = 1.d0, deux = 2.d0
    integer(kind=8), parameter :: nbPropMaxi = 9
    integer(kind=8) :: nbProp
    integer(kind=8) :: propCode(nbPropMaxi)
    character(len=16) :: propName(nbPropMaxi)
    real(kind=8) :: propVale(nbPropMaxi)
    integer(kind=8), parameter :: nbParaMaxi = 1
    integer(kind=8) :: nbPara
    character(len=8) :: paraName(nbParaMaxi)
    real(kind=8) :: paraVale(nbParaMaxi)
    integer(kind=8) :: i, irep, j
    real(kind=8) :: coef1, coef2, coef3, e, e1, e2
    real(kind=8) :: e3
    real(kind=8) :: passag(6, 6), d1orth(6, 6), work(6, 6)
    real(kind=8) :: nu, nu12, nu21, nu13, nu23, nu31, nu32
!
! --------------------------------------------------------------------------------------------------
!
    d1 = zero
!
    if (time .eq. r8vide()) then
        nbPara = 0
    else
        nbPara = 1
        paraName(1) = 'INST'
        paraVale(1) = time
    end if
    ASSERT(nbPara .le. nbParaMaxi)
!
    d1orth = zero
    work = zero

    if (materPara%elasID == ELAS_ISOT) then
        propName(1) = 'E'
        propName(2) = 'NU'
        nbProp = 2
        call rcvalb(materPara%schemePara%fami, &
                    materPara%schemePara%kpg, &
                    materPara%schemePara%ksp, &
                    poum, &
                    materPara%jvMaterCode, ' ', materPara%elasKeyword, &
                    nbPara, paraName, [paraVale], &
                    nbProp, propName, propVale, propCode, 1)
        e = propVale(1)
        nu = propVale(2)
        coef1 = un/e
        coef2 = -nu/e
        coef3 = deux*(un+nu)/e
        d1(1, 1) = coef1
        d1(1, 2) = coef2
        d1(1, 3) = coef2
        d1(2, 1) = coef2
        d1(2, 2) = coef1
        d1(2, 3) = coef2
        d1(3, 1) = coef2
        d1(3, 2) = coef2
        d1(3, 3) = coef1
        d1(4, 4) = coef3
        d1(5, 5) = coef3
        d1(6, 6) = coef3

    else if (materPara%elasID == ELAS_ORTH) then
        if (.not. chckLCSDefine(materPara%lcsPara)) then
            call utmess("F", "ALGORITH8_20")
        end if
        propName(1) = 'E_L'
        propName(2) = 'E_T'
        propName(3) = 'E_N'
        propName(4) = 'NU_LT'
        propName(5) = 'NU_LN'
        propName(6) = 'NU_TN'
        propName(7) = 'G_LT'
        propName(8) = 'G_LN'
        propName(9) = 'G_TN'
        nbProp = 9
        call rcvalb(materPara%schemePara%fami, &
                    materPara%schemePara%kpg, &
                    materPara%schemePara%ksp, &
                    poum, &
                    materPara%jvMaterCode, ' ', materPara%elasKeyword, &
                    nbPara, paraName, [paraVale], &
                    nbProp, propName, propVale, propCode, 1)
        e1 = propVale(1)
        e2 = propVale(2)
        e3 = propVale(3)
        nu12 = propVale(4)
        nu13 = propVale(5)
        nu23 = propVale(6)
        nu21 = e2*nu12/e1
        nu31 = e3*nu13/e1
        nu32 = e3*nu23/e2
!
        d1orth(1, 1) = un/e1
        d1orth(1, 2) = -nu21/e2
        d1orth(1, 3) = -nu31/e3
        d1orth(2, 2) = un/e2
        d1orth(2, 3) = -nu32/e3
        d1orth(3, 3) = un/e3
        d1orth(2, 1) = d1orth(1, 2)
        d1orth(3, 1) = d1orth(1, 3)
        d1orth(3, 2) = d1orth(2, 3)
!
        d1orth(4, 4) = un/propVale(7)
        d1orth(5, 5) = un/propVale(8)
        d1orth(6, 6) = un/propVale(9)
!
! ----   CALCUL DE LA MATRICE DE PASSAGE DU REPERE D'ORTHOTROPIE AU
! ----   REPERE GLOBAL POUR L'INVERSE DE LA MATRICE DE HOOKE
!        ---------------------------------------------------
        call d1pa3d(materPara%lcsPara%lcsAngle, irep, passag)
!
! ----   'INVERSE' DU TENSEUR D'ELASTICITE DANS LE REPERE GLOBAL :
! ----    D1_GLOB = PASSAG_T * D1_ORTH * PASSAG
! ----    (ON NE FAIT REELLEMENT LE PRODUIT QUE SI LA MATRICE
! ----     DE PASSAGE N'EST PAS L'IDENTITE)
!        ----------------------------------
        ASSERT((irep .eq. 1) .or. (irep .eq. 0))
        if (irep .eq. 1) then
            call utbtab('ZERO', 6, 6, d1orth, passag, work, d1)
        else if (irep .eq. 0) then
            do i = 1, 6
                do j = 1, 6
                    d1(i, j) = d1orth(i, j)
                end do
            end do
        end if

    else if (materPara%elasID == ELAS_ISTR) then
        if (.not. chckLCSDefine(materPara%lcsPara)) then
            call utmess("F", "ALGORITH8_20")
        end if
        propName(1) = 'E_L'
        propName(2) = 'E_N'
        propName(3) = 'NU_LT'
        propName(4) = 'NU_LN'
        propName(5) = 'G_LN'
        nbProp = 5
!
! ----   INTERPOLATION DES COEFFICIENTS EN FONCTION DE LA TEMPERATURE
! ----   ET DU TEMPS
!        -----------
        call rcvalb(materPara%schemePara%fami, &
                    materPara%schemePara%kpg, &
                    materPara%schemePara%ksp, &
                    poum, &
                    materPara%jvMaterCode, ' ', materPara%elasKeyword, &
                    nbPara, paraName, [paraVale], &
                    nbProp, propName, propVale, propCode, 1)
!
        e1 = propVale(1)
        e3 = propVale(2)
        nu12 = propVale(3)
        nu13 = propVale(4)
        nu31 = e3*nu13/e1
!
        d1orth(1, 1) = un/e1
        d1orth(1, 2) = -nu12/e1
        d1orth(1, 3) = -nu31/e3
        d1orth(2, 1) = d1orth(1, 2)
        d1orth(2, 2) = un/e1
        d1orth(2, 3) = -nu31/e3
        d1orth(3, 1) = d1orth(1, 3)
        d1orth(3, 2) = d1orth(2, 3)
        d1orth(3, 3) = un/e3
        d1orth(4, 4) = deux*(un+nu12)/e1
        d1orth(5, 5) = un/propVale(5)
        d1orth(6, 6) = d1orth(5, 5)
!
! ----   CALCUL DE LA MATRICE DE PASSAGE DU REPERE D'ORTHOTROPIE AU
! ----   REPERE GLOBAL POUR L'INVERSE DE LA MATRICE DE HOOKE
!        ---------------------------------------------------
        call d1pa3d(materPara%lcsPara%lcsAngle, irep, passag)
!
! ----   'INVERSE' DU TENSEUR D'ELASTICITE DANS LE REPERE GLOBAL :
! ----    D_GLOB = PASSAG_T * D_ORTH * PASSAG
! ----    (ON NE FAIT REELLEMENT LE PRODUIT QUE SI LA MATRICE
! ----     DE PASSAGE N'EST PAS L'IDENTITE)
!        ----------------------------------
        ASSERT((irep .eq. 1) .or. (irep .eq. 0))
        if (irep .eq. 1) then
            call utbtab('ZERO', 6, 6, d1orth, passag, &
                        work, d1)
        else if (irep .eq. 0) then
            do i = 1, 6
                do j = 1, 6
                    d1(i, j) = d1orth(i, j)
                end do
            end do
        end if

    else
        call utmess('F', 'ELEMENTS_15', sk=materPara%elasKeyword)
    end if
!
end subroutine
