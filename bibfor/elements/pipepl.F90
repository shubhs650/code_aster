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
! aslint: disable=W0413
!
subroutine pipepl(materPara, ndim, relaComp, typmod, &
                  tau, &
                  sigm, vim, epsp, epsd, a0, &
                  a1, a2, a3, etas)
!
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterc/r8vide.h"
#include "asterf_types.h"
#include "asterfort/rcfonc.h"
#include "asterfort/rctrac.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
#include "asterfort/zerop2.h"
#include "blas/ddot.h"
!
    type(Material_Para), intent(in) :: materPara
    character(len=8), intent(in) :: typmod(2)
    character(len=16), intent(in) :: relaComp
    integer(kind=8), intent(in) :: ndim
    real(kind=8) :: epsp(6), epsd(6), tau
    real(kind=8) :: vim(2), sigm(6)
    real(kind=8) :: a0, a1, a2, a3, etas
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (PILOTAGE - PRED_ELAS)
!
! LOI DE COMPORTEMENT PLASTIQUE VMIS_ISOT_*
!
! --------------------------------------------------------------------------------------------------
!
! IN  TAU    : 2ND MEMBRE DE L'EQUATION F(ETA)=TAU
! IN  SIGM   : CONTRAINTE EN T-
! IN  VIM    : VARIABLES INTERNES EN T-
! IN  EPSP   : CORRECTION DE DEFORMATIONS DUES AUX CHARGES FIXES
! IN  EPSD   : CORRECTION DE DEFORMATIONS DUES AUX CHARGES PILOTEES
! OUT A0     : LINEARISATION DU CRITERE : FEL = A0 + A1*ETA
! OUT A1     : IDEM A0
! OUT A2     : IDEM A0 POUR LA 2E SOLUTION EVENTUELLE. R8VIDE SINON
! OUT A3     : IDEM A1 POUR LA 2E SOLUTION EVENTUELLE. R8VIDE SINON
! OUT ETAS   : SI PAS DE SOLUTION : LE MINIMUM. R8VIDE SINON
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: kpgFPG1 = 1, kspFPG1 = 1
    character(len=8), parameter :: famiFPG1 = "FPG1"
    type(Material_Para) :: materParaFPG1
    character(len=8), parameter :: poum = "+"
    integer(kind=8), parameter :: nbProp = 4
    integer(kind=8) :: propCode(nbProp)
    character(len=16) :: propName(nbProp)
    real(kind=8) :: propVale(nbProp)
    aster_logical :: cplan
    integer(kind=8) :: ndimsi, k, nrac, jprol, jvale, nbvale
    real(kind=8) :: sigmh, epsph, epsdh, s0h, s1h, s0(6), s1(6)
    real(kind=8) :: p0, p1, p2, eta, rac(2)
    real(kind=8) :: young, nu, deuxmu, rp, h, et, sy
    blas_int :: b_incx, b_incy, b_n
    real(kind=8), parameter :: kron(6) = (/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0/)
!
! --------------------------------------------------------------------------------------------------
!

    ndimsi = 2*ndim
    cplan = (typmod(1) .eq. 'C_PLAN  ')
    if (cplan) then
        call utmess('F', 'PILOTAGE_1')
    end if

! - Copy material parameters with other scheme parameters
    call copyMaterPara(materPara, famiFPG1, kpgFPG1, kspFPG1, &
                       materParaFPG1)

! - LECTURE DES CARACTERISTIQUES
    if (relaComp .eq. 'VMIS_ISOT_LINE') then
        propName(1) = 'E'
        propName(2) = 'NU'
        propName(3) = 'SY'
        propName(4) = 'D_SIGM_EPSI'
        call rcvalb(materParaFPG1%schemePara%fami, &
                    materParaFPG1%schemePara%kpg, &
                    materParaFPG1%schemePara%ksp, &
                    poum, &
                    materParaFPG1%jvMaterCode, ' ', 'ELAS', &
                    0, ' ', [0.d0], &
                    2, propName, propVale, &
                    propCode, 2)
        call rcvalb(materParaFPG1%schemePara%fami, &
                    materParaFPG1%schemePara%kpg, &
                    materParaFPG1%schemePara%ksp, &
                    poum, &
                    materParaFPG1%jvMaterCode, ' ', 'ECRO_LINE', &
                    0, ' ', [0.d0], &
                    2, propName(3), propVale(3), &
                    propCode(3), 2)
        young = propVale(1)
        nu = propVale(2)
        sy = propVale(3)
        et = propVale(4)
        h = young*et/(young-et)
        rp = sy+h*vim(1)
!
    else
        call rcvalb(materParaFPG1%schemePara%fami, &
                    materParaFPG1%schemePara%kpg, &
                    materParaFPG1%schemePara%ksp, &
                    poum, &
                    materParaFPG1%jvMaterCode, ' ', 'ELAS', &
                    0, ' ', [0.d0], &
                    1, 'NU', propVale, &
                    propCode, 2)
        nu = propVale(1)
        call rctrac(materParaFPG1%jvMaterCode, 1, 'SIGM', 0.d0, jprol, &
                    jvale, nbvale, young)
        call rcfonc('V', 1, jprol, jvale, nbvale, &
                    p=vim(1), rp=rp)
    end if
!
    deuxmu = young/(1.d0+nu)
!
! ======================================================================
!                CALCUL DES DEFORMATIONS POUR LINEARISATION
! ======================================================================
! - PARTITION TRACE / DEVIATEUR
    sigmh = (sigm(1)+sigm(2)+sigm(3))/3
    epsph = (epsp(1)+epsp(2)+epsp(3))/3
    epsdh = (epsd(1)+epsd(2)+epsd(3))/3
!
    s0h = deuxmu*epsph+sigmh
    s1h = deuxmu*epsdh
    do k = 1, ndimsi
        s0(k) = sigm(k)+deuxmu*epsp(k)-s0h*kron(k)
        s1(k) = deuxmu*epsd(k)-s1h*kron(k)
    end do

! - COEFFICIENTS DE LA FORME QUADRATIQUE DU CRITERE
! - FEL = SQRT(P0 + 2P1 ETA + P2 ETA**2) - 1
    b_n = to_blas_int(ndimsi)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    p0 = ddot(b_n, s0, b_incx, s0, b_incy)*(1.5d0/rp**2)
    b_n = to_blas_int(ndimsi)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    p1 = ddot(b_n, s0, b_incx, s1, b_incy)*(1.5d0/rp**2)
    b_n = to_blas_int(ndimsi)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    p2 = ddot(b_n, s1, b_incx, s1, b_incy)*(1.5d0/rp**2)

    if (p2 .eq. 0) then
! ----- POINT A DEVIATEUR NUL : PAS DE PILOTAGE POSSIBLE
        a0 = 0.d0
        a1 = 0.d0
        a2 = 0.d0
        a3 = 0.d0
    else
! ----- RECHERCHE DES INTERSECTIONS ELLIPSE / DROITE
        call zerop2(2*p1/p2, (p0-(1+tau)**2)/p2, rac, nrac)

! ----- PAS DE SOLUTION : POINT LE PLUS PROCHE
        if (nrac .eq. 0) then
            etas = -p1/p2

! ----- UNE OU DEUX SOLUTIONS : ON LINEARISE AUTOUR DES DEUX
        else if (nrac .eq. 1) then
            eta = rac(1)
            a1 = (p2*eta+p1)/(1+tau)
            a0 = tau-a1*eta
            a2 = r8vide()
            a3 = r8vide()
        else
            eta = rac(1)
            a1 = (p2*eta+p1)/(1+tau)
            a0 = tau-a1*eta
            eta = rac(2)
            a3 = (p2*eta+p1)/(1+tau)
            a2 = tau-a3*eta
        end if
    end if
end subroutine
