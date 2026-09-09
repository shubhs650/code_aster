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
subroutine betjpl(BEHInteg, &
                  mod, nmat, mater, sig, vin, &
                  dsde)
!
    use Behaviour_type
    implicit none
!
!       BETON_DOUBLE_DP: LOI ELASTO PLASTIQUE AVEC DOUBLE CRITERE DE
!       PLASTICITE AVEC UN SEUIL EN COMPRESSION ET UN SEUIL EN TRACTION
!       MATRICE SYMETRIQUE DE COMPORTEMENT TANGENT ELASTO_PLASTIQUE
!       EN VITESSE A T OU T+DT
!       ----------------------------------------------------------------
!       IN  MOD    :  TYPE DE MODELISATION
!           NMAT   :  DIMENSION MATER
!           TEMP   :  TEMPERATURE
!           MATER  :  COEFFICIENTS MATERIAU
!           SIG    :  CONTRAINTES
!           VIN    :  VARIABLES INTERNES
!       OUT DSDE   :  MATRICE DE COMPORTEMENT TANGENT = DSIG/DEPS
!       ----------------------------------------------------------------
#include "jeveux.h"
#include "asterfort/betfpp.h"
#include "asterfort/lcdevi.h"
#include "asterfort/lcopli.h"
#include "asterfort/lcprte.h"
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"
    type(Behaviour_Integ), intent(in) :: BEHInteg
    integer(kind=8) :: nmat, nseuil
    real(kind=8) :: un, zero, rac2, deux, trois
    parameter(deux=2.d0)
    parameter(trois=3.d0)
    parameter(un=1.d0)
    parameter(zero=0.d0)
!
    real(kind=8) :: vin(*), sig(6)
    real(kind=8) :: hook(6, 6), dsde(6, 6), vtmp(6)
!
    real(kind=8) :: mater(nmat, 2)
!
    character(len=8) :: mod
    real(kind=8) :: trav1(6), trav2(6), pi0(6), dev(6)
    real(kind=8) :: sigeq, matr1(6, 6)
    real(kind=8) :: fc, ft, beta, kuc, kut, ke
    real(kind=8) :: a, b, c, d
    real(kind=8) :: pc, pt, dfcdlc, dftdlt, dfcds(6), dftds(6)
    real(kind=8) :: coef1, coef2, hdfcds(6), hdftds(6)
    real(kind=8) :: cc, ccc, tt, ttt, ct, tc, discr
    integer(kind=8) :: iadzi, iazk24
    character(len=8) :: nomail
!       ----------------------------------------------------------------
    integer(kind=8) :: ndt, ndi
!     ------------------------------------------------------------------
    common/tdim/ndt, ndi
!       ----------------------------------------------------------------
    data pi0/un, un, un, zero, zero, zero/
!       ----------------------------------------------------------------
!
!
! --- INITIALISATION
!
    kuc = 0
    kut = 0
    rac2 = sqrt(deux)
    beta = mater(3, 2)
!
    a = rac2*(beta-un)/(deux*beta-un)
    b = rac2/trois*beta/(deux*beta-un)
    c = rac2
    d = deux*rac2/trois
!
    call lcopli('ISOTROPE', mod, mater(1, 1), hook)
!
    pc = vin(1)
    pt = vin(2)
    nseuil = int(vin(4)+0.5d0)
!
! --- CONTRAINTE EQUIVALENTE
!
    call lcdevi(sig, dev)
    sigeq = sqrt(1.5d0*dot_product(dev(1:ndt), dev(1:ndt)))
    if (sigeq .eq. zero) then
        call tecael(iadzi, iazk24)
        nomail = zk24(iazk24-1+3) (1:8)
        call utmess('A', 'ALGORITH_48', sk=nomail)
        sigeq = 1.d0
    end if
!
! --  CALCUL DES ECROUISSAGES ET DERIVES DES COURBES D'ADOUCISSEMENT
!
    call betfpp(BEHInteg, &
                mater, nmat, pc, pt, &
                nseuil, fc, ft, dfcdlc, dftdlt, &
                kuc, kut, ke)
!
! --- DERIVEES DU CRITERE EN COMPRESSION
!
    if (nseuil .eq. 1 .or. nseuil .eq. 3) then
!
        coef1 = un/(rac2*b*sigeq)
        coef2 = a/(trois*b)
        trav1(1:ndt) = coef1*dev(1:ndt)
        trav2(1:ndt) = coef2*pi0(1:ndt)
        dfcds(1:ndt) = trav1(1:ndt)+trav2(1:ndt)
    end if
!
! --- DERIVEES DU CRITERE EN TRACTION
!
    if (nseuil .eq. 2 .or. nseuil .eq. 3) then
!
        coef1 = un/(rac2*d*sigeq)
        coef2 = c/(trois*d)
        trav1(1:ndt) = coef1*dev(1:ndt)
        trav2(1:ndt) = coef2*pi0(1:ndt)
        dftds(1:ndt) = trav1(1:ndt)+trav2(1:ndt)
    end if
!
! --- DERIVEES DU CRITERE EN TRACTION AVEC PROJECTION AU SOMMET DES
!     CONES DE TRACTION ET DE COMPRESSION
!
    if (nseuil .eq. 11 .or. nseuil .eq. 33) then
        coef2 = a/(trois*b)
        dfcds(1:ndt) = coef2*pi0(1:ndt)
    end if
    if (nseuil .eq. 22 .or. nseuil .eq. 33) then
        coef2 = c/(trois*d)
        dftds(1:ndt) = coef2*pi0(1:ndt)
    end if
!
! --- MATRICE DE COMPORTEMENT TANGENT
!
    if (nseuil .eq. 3 .or. nseuil .eq. 33) then
        hdfcds(1:ndt) = matmul(hook(1:ndt, 1:ndt), dfcds(1:ndt))
        hdftds(1:ndt) = matmul(hook(1:ndt, 1:ndt), dftds(1:ndt))
        cc = dot_product(dfcds(1:ndt), hdfcds(1:ndt))
        tc = dot_product(dftds(1:ndt), hdfcds(1:ndt))
        ct = dot_product(dfcds(1:ndt), hdftds(1:ndt))
        tt = dot_product(dftds(1:ndt), hdftds(1:ndt))
        ccc = cc+dfcdlc
        ttt = tt+dftdlt
        discr = -un/(ccc*ttt-ct*tc)
        vtmp(1:ndt) = ((discr*ttt))*hdfcds(1:ndt)
        call lcprte(hdfcds, vtmp, dsde)
        vtmp(1:ndt) = ((discr*ccc))*hdftds(1:ndt)
        call lcprte(hdftds, vtmp, matr1)
        dsde(1:ndt, 1:ndt) = matr1(1:ndt, 1:ndt)+dsde(1:ndt, 1:ndt)
        vtmp(1:ndt) = (discr*ct)*hdftds(1:ndt)
        call lcprte(hdfcds, vtmp, matr1)
        dsde(1:ndt, 1:ndt) = dsde(1:ndt, 1:ndt)-matr1(1:ndt, 1:ndt)
        vtmp(1:ndt) = (discr*tc)*hdfcds(1:ndt)
        call lcprte(hdftds, vtmp, matr1)
        dsde(1:ndt, 1:ndt) = dsde(1:ndt, 1:ndt)-matr1(1:ndt, 1:ndt)
        dsde(1:ndt, 1:ndt) = hook(1:ndt, 1:ndt)+dsde(1:ndt, 1:ndt)
    end if
!
    if (nseuil .eq. 2 .or. nseuil .eq. 22) then
        hdftds(1:ndt) = matmul(hook(1:ndt, 1:ndt), dftds(1:ndt))
        tt = dot_product(dftds(1:ndt), hdftds(1:ndt))
        ttt = tt+dftdlt
        discr = -un/ttt
        vtmp(1:ndt) = discr*hdftds(1:ndt)
        call lcprte(hdftds, vtmp, dsde)
        dsde(1:ndt, 1:ndt) = hook(1:ndt, 1:ndt)+dsde(1:ndt, 1:ndt)
    end if
!
    if (nseuil .eq. 1 .or. nseuil .eq. 11) then
        hdfcds(1:ndt) = matmul(hook(1:ndt, 1:ndt), dfcds(1:ndt))
        cc = dot_product(dfcds(1:ndt), hdfcds(1:ndt))
        ccc = cc+dfcdlc
        discr = -un/ccc
        vtmp(1:ndt) = discr*hdfcds(1:ndt)
        call lcprte(hdfcds, vtmp, dsde)
        dsde(1:ndt, 1:ndt) = hook(1:ndt, 1:ndt)+dsde(1:ndt, 1:ndt)
    end if
!
end subroutine
