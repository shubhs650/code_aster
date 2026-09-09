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
! aslint: disable=W1504
!
subroutine lcplnf(BEHInteg, &
                  rela_comp, vind, nbcomm, nmat, cpmono, &
                  materf, iter, nvi, itmax, &
                  toler, pgl, nfs, nsg, toutms, &
                  hsr, dt, dy, yd, yf, &
                  vinf, sigd, sigf, &
                  deps, nr, mod, timef, &
                  codret)
!
    use Behaviour_type
!
    implicit none
!
    type(Behaviour_Integ), intent(in) :: BEHInteg
!   POST-TRAITEMENTS SPECIFIQUES AUX LOIS
!
!   CORRESPONDANCE ENTRE LES VARIABLES INTERNES ET LES EQUATIONS
!          DU SYSTEME DIFFERENTIEL APRES INTEGRATION
!
!   CAS GENERAL :
!      COPIE DES YF DANS VINF
!      LA DERNIERE C'EST TOUJOURS L'INDICATEUR PLASTIQUE
!
!   CAS PARTICULIER DU  MONOCRISTAL  :
!       ON GARDE 1 VARIABLE INTERNE PAR SYSTEME DE GLISSEMENT SUR 3
!       DEFORMATION PLASTIQUE EQUIVALENTE CUMULEE MACROSCOPIQUE
! ----------------------------------------------------------------
!  IN
!     LOI    :  NOM DE LA LOI
!     VIND   :  VARIABLE INTERNES A T
!     MATERF :  COEF MATERIAU A T+DT
!     NBCOMM :  INCIDES DES COEF MATERIAU
!     NMAT   :  DIMENSION MATER ET DE NBCOMM
!     NVI    :  NOMBRE DE VARIABLES INTERNES
!     DT     : INCREMENT DE TEMPS
!     NR     : DIMENSION VECTEUR INCONNUES (YF/DY)
!     YF     : EQUATIONS DU COMPORTEMENT INTEGRES A T+DT
!     DY     : INCREMENT DES VARIABLES INTERNES
!     TIMED  : INSTANT T
!     TIMEF  : INSTANT T+DT
!  OUT
!     VINF   :  VARIABLES INTERNES A T+DT
! ----------------------------------------------------------------
#include "asterfort/irrlnf.h"
#include "asterfort/lcdpec.h"
#include "asterfort/lcopli.h"
#include "asterfort/lkilnf.h"
#include "asterfort/srilnf.h"
    integer(kind=8) :: ndt, nvi, nmat, ndi, nbcomm(nmat, 3), iter, itmax, nr, codret
    integer(kind=8) :: nfs, nsg, i
    real(kind=8) :: materf(nmat, 2), timef
    real(kind=8) :: pkc, m13, dtot, hookf(6, 6)
    real(kind=8) :: yd(*), vind(*), toler, pgl(3, 3), dt
    real(kind=8) :: toutms(nfs, nsg, 6), hsr(nsg, nsg), dy(*), yf(*), vinf(*)
    character(len=16) :: rela_comp
    character(len=24) :: cpmono(5*nmat+1)
    character(len=8) :: mod
    real(kind=8) :: sigf(6), deps(*), sigd(6)
!
    common/tdim/ndt, ndi
! --- -------------------------------------------------------------
!
!     MISE A JOUR DE SIGF , VINF
    sigf(1:ndt) = yf(1:ndt)
!
    if (rela_comp .eq. 'MONOCRISTAL') then
! ---    DEFORMATION PLASTIQUE EQUIVALENTE CUMULEE MACROSCOPIQUE
        call lcdpec(BEHInteg, &
                    vind, nbcomm, nmat, ndt, cpmono, &
                    materf, iter, nvi, itmax, toler, &
                    pgl, nfs, nsg, toutms, hsr, &
                    dt, dy, yd, vinf, &
                    sigf, deps, nr, mod, &
                    codret)
!
    else if (rela_comp .eq. 'IRRAD3M') then
        call irrlnf(nmat, materf, yf(ndt+1), 1.0d0, vinf)
    else if (rela_comp .eq. 'LETK') then
        call lkilnf(nvi, vind, nmat, materf, dt, &
                    sigd, nr, yd, yf, deps, &
                    vinf)
    else if (rela_comp .eq. 'LKR') then
        call srilnf(nvi, vind, nmat, materf, dt, &
                    nr, yf, deps, vinf)
    else if (rela_comp .eq. 'HAYHURST') then
!        DEFORMATION PLASTIQUE CUMULEE
        vinf(7) = yf(ndt+1)
!        H1
        vinf(8) = yf(ndt+2)
!        H2
        vinf(9) = yf(ndt+3)
!        PHI
        pkc = materf(11, 2)
        m13 = -1.d0/3.d0
        vinf(10) = 1.d0-(1.d0+pkc*timef)**m13
!        DEFORMATION PLASTIQUE
!        D
        vinf(11) = yf(ndt+4)
        dtot = (1.d0-vinf(11))
        call lcopli('ISOTROPE', mod, materf(1, 1), hookf)
        sigf(1:ndt) = matmul(hookf(1:ndt, 1:ndt), yf(1:ndt))
        sigf(1:ndt) = dtot*sigf(1:ndt)
        do i = 1, ndt
            vinf(i) = yf(i)
        end do
        vinf(nvi) = iter
    else
!        CAS GENERAL :
        vinf(1:nvi-1) = yf(ndt+1:ndt+nvi-1)
        vinf(nvi) = iter
    end if
!
end subroutine
