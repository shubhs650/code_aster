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
subroutine lcconv(rela_comp, yd, dy, ddy, &
                  nr, itmax, toler, iter, intg, &
                  nmat, mater, r, rini, epstr, &
                  typess, essai, cutLevel, nvi, &
                  vinf, &
                  iret)

    implicit none
!     ROUTINE D AIGUILLAGE
!     ----------------------------------------------------------------
!     CONTROLE DE LA CONVERGENCE DE LA METHODE DE NEWTON (LCPLNL):
!
!                     - CONTROLE DU NOMBRE D ITERATIONS
!                     - CONTROLE DE LA PRECISION DE CONVERGENCE
!                     - CONTROLE DE LA VALIDITE SOLUTION A CONVERGENCE
!                     - CONTROLE DES RE-INTEGRATIONS EVENTUELLES
!                     - CONTROLE DU REDECOUPAGE DU PAS DE TEMPS
!
!     ----------------------------------------------------------------
!     IN  LOI    :  MODELE DE COMPORTEMENT
!         TYPESS :  TYPE DE SOLUTION D ESSAI POUR DY(DEPEND DU MODELE)
!                    > VOIR XXXCVG ET XXXINI
!         ESSAI  :  VALEUR SOLUTION D ESSAI
!         ITMAX  :  NB MAXI D ITERATIONS LOCALES
!         TOLER  :  TOLERANCE A CONVERGENCE
!         ITER   :  NUMERO ITERATION COURANTE
!         INTG   :  NUMERO INTEGRATION COURANTE
!         NR     :  DIMENSION DY DDY
!         DY     :  VECTEUR SOLUTION = ( DSIG DVIN (DEPS3) )
!         DDY    :  VECTEUR CORRECTION SUR LA SOLUTION
!         NVI    :  NOMBRE DE VARIABLES INTERNES
!         VINF   :  VARIABLES INTERNES A L'INSTANT T+DT
!
!     OUT IRET = 0:  CONVERGENCE
!         IRET = 1:  ITERATION SUIVANTE
!         IRET = 2:  RE-INTEGRATION
!         IRET = 3:  REDECOUPAGE DU PAS DE TEMPS
!         (VINF) UNIQUEMENT POUR LETK  - ETAT PLASTIQUE DESACTIVE?
!     ----------------------------------------------------------------
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/cvmcvg.h"
#include "asterfort/irrcvg.h"
#include "asterfort/lccong.h"
#include "asterfort/lcmmcv.h"
#include "asterfort/lkicvg.h"
#include "asterfort/sricvg.h"
    integer(kind=8) :: typess, itmax, iter, intg, nr, cutLevel
    integer(kind=8) :: iret, nmat, nvi
    real(kind=8) :: toler, essai, ddy(*), dy(*), r(*), rini(*), yd(*)
    real(kind=8) :: mater(nmat, 2), epstr(6), vinf(nvi)
    character(len=16), intent(in) :: rela_comp
!     ----------------------------------------------------------------
!
    if (rela_comp .eq. 'VISCOCHAB') then
!
        call cvmcvg(dy, ddy, nr, itmax, toler, &
                    iter, intg, typess, essai, cutLevel, &
                    iret)
!
    else if (rela_comp .eq. 'MONOCRISTAL') then
!
        call lcmmcv(yd, dy, ddy, nr, itmax, &
                    toler, iter, r, rini, epstr, &
                    iret)
!
    else if (rela_comp .eq. 'IRRAD3M') then
!
        call irrcvg(dy, ddy, nr, nmat, mater, &
                    itmax, toler, iter, r, rini, &
                    iret)
!
    else if (rela_comp .eq. 'LETK') then
!
        call lkicvg(nr, itmax, toler, iter, r, &
                    nvi, vinf, dy, iret)
!
    else if (rela_comp .eq. 'LKR') then
!
        call sricvg(nr, itmax, toler, iter, r, &
                    nvi, vinf, dy, iret)
!
    else
!
        call lccong(nr, itmax, toler, iter, r, &
                    rini, yd, dy, iret)
!
    end if
!
    ASSERT(iret .ge. 0)
    ASSERT(iret .le. 3)
!
end subroutine
