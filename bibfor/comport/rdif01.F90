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
! aslint: disable=W1306,W1504
!
subroutine rdif01(materPara, &
                  relaComp, typmod1, &
                  matcst, nbcomm, cpmono, nfs, &
                  nsg, toutms, nvi, nmat, vini, &
                  cothe, coeff, dcothe, dcoeff, pgl, &
                  nbphas, coel, x, dtime, neps, &
                  epsd, detot, dvin, nhsr, numhsr, &
                  hsr, itmax, toler, iret)
!
    use MaterialPara_type
    implicit none
!
#include "asterfort/calsig.h"
#include "asterfort/coefft.h"
#include "asterfort/lcdvin.h"
#include "asterfort/lcmmon.h"
#include "asterfort/lcmmop.h"
!
    type(Material_Para), intent(in) :: materPara
    character(len=16), intent(in) :: relaComp
    character(len=8), intent(in) :: typmod1
!
! --------------------------------------------------------------------------------------------------
!
!     INTEGRATION DE LOIS DE COMPORTEMENT ELASTO-VISCOPLASTIQUE
!     PAR UNE METHODE DE RUNGE KUTTA
!
! --------------------------------------------------------------------------------------------------
!
!         MATCST  :  NATURE DES PARAMETRES INELASTIQUES
!         NVI     :  NOMBRE DE VARIABLES INTERNES
!         NMAT    :  NOMBRE DE PARAMETRES MATERIAU INELASTIQUE
!         VINI    :  VARIABLES INTERNES A T
!         COTHE   :  COEFFICIENTS MATERIAU ELASTIQUE A T
!         COEFF   :  COEFFICIENTS MATERIAU INELASTIQUE A T
!         DCOTHE  :  DELTA COEFFICIENTS MATERIAU ELASTIQUE A T+DT
!         DCOEFF  :  DELTA COEFFICIENTS MATERIAU INELASTIQUE A T+DT
!         COEL    :  COEFFICIENTS D'ELASTICITE
!         X       :  INTERVALE DE TEMPS ADAPTATIF
!         DTIME   :  INTERVALE DE TEMPS
!         EPSD    :  DEFORMATION TOTALE A T
!         DETOT   :  INCREMENT DE DEFORMATION TOTALE
!         DVIN    :  DERIVEES DES VARIABLES INTERNES A T
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nmat, nvi, nbcomm(nmat, 3), itens
    integer(kind=8) :: nbphas, nfs, iret, itmax, nsg, nhsr, numhsr(*), neps
    character(len=24) :: cpmono(5*nmat+1)
    character(len=3) :: matcst
    real(kind=8) :: pgl(3, 3), toler, x, dtime, coel(nmat)
    real(kind=8) :: cothe(nmat), dcothe(nmat), coeff(nmat), dcoeff(nmat)
    real(kind=8) :: epsd(6), detot(6), coeft(nmat), xm, sigi(6)
    real(kind=8) :: vini(nvi), dvin(nvi), hsr(nsg, nsg, nhsr), evi(6)
!     POUR GAGNER EN TEMPS CPU
    real(kind=8) :: toutms(*)
    character(len=8) :: fami
    integer(kind=8) :: jvMaterCode, kpg, ksp
!
! --------------------------------------------------------------------------------------------------
!
    jvMaterCode = materPara%jvMaterCode
    fami = materPara%schemePara%fami
    kpg = materPara%schemePara%kpg
    ksp = materPara%schemePara%ksp

    if (relaComp .eq. 'MONOCRISTAL') then
!       PAS DE VARIATION DES COEF AVEC LA TEMPERATURE
        xm = 0.d0
        call coefft(cothe, coeff, dcothe, dcoeff, xm, &
                    dtime, coeft, nmat, coel)
        call lcmmon(fami, kpg, ksp, relaComp, nbcomm, &
                    cpmono, nmat, nvi, vini, x, &
                    dtime, pgl, typmod1, coeft, neps, &
                    epsd, detot, coel, dvin, nfs, &
                    nsg, toutms, hsr(1, 1, 1), itmax, toler, &
                    iret)
!
    else if (relaComp .eq. 'POLYCRISTAL') then
!       PAS DE VARIATION DES COEF AVEC LA TEMPERATURE
        xm = 0.d0
        call coefft(cothe, coeff, dcothe, dcoeff, xm, &
                    dtime, coeft, nmat, coel)
        call lcmmop(fami, kpg, ksp, relaComp, nbcomm, &
                    cpmono, nmat, nvi, vini, x, &
                    dtime, typmod1, coeft, epsd, detot, &
                    coel, nbphas, nfs, nsg, toutms, &
                    dvin, nhsr, numhsr, hsr, itmax, &
                    toler, iret)

    else
        do itens = 1, 6
            evi(itens) = vini(itens)
        end do
        call coefft(cothe, coeff, dcothe, dcoeff, x, &
                    dtime, coeft, nmat, coel)
        call calsig(fami, kpg, ksp, evi, typmod1, &
                    relaComp, vini, x, dtime, epsd, &
                    detot, nmat, coel, sigi)
        call lcdvin(fami, kpg, ksp, relaComp, typmod1, &
                    jvMaterCode, matcst, nvi, nmat, vini, &
                    coeft, x, dtime, sigi, dvin, &
                    iret)
!
    end if
end subroutine
