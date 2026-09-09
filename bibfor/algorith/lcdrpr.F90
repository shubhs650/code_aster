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
subroutine lcdrpr(fami, typmod, option, imate, sigm, &
                  depsm, vim, &
                  vip, sig, dsidep, iret)
!
    implicit none
!
#include "asterfort/dpmate.h"
#include "asterfort/dpvpdi.h"
#include "asterfort/lcdpli.h"
#include "asterfort/lcdppa.h"
#include "asterfort/get_varc.h"
!
!
    character(len=*), intent(in) :: fami
    integer(kind=8) :: imate, iret
    real(kind=8) :: depsm(6), vim(*), vip(*), sig(6), dsidep(6, 6)
    real(kind=8) :: sigm(6)
    character(len=8) :: typmod(*)
    character(len=16) :: option
! ======================================================================
! --- LOI DE COMPORTEMENT DE TYPE DRUCKER PRAGER -----------------------
! --- ELASTICITE ISOTROPE ----------------------------------------------
! --- PLASTICITE DE VON MISES + TERME DE TRACE -------------------------
! --- ECROUISSAGE ISOTROPE LINEAIRE ------------------------------------
! ======================================================================
! IN  OPTION  OPTION DE CALCUL (RAPH_MECA, RIGI_MECA_TANG OU FULL_MECA)
! IN  IMATE   NATURE DU MATERIAU
! IN  EPSM    CHAMP DE DEFORMATION EN T-
! IN  DEPS    INCREMENT DU CHAMP DE DEFORMATION
! IN  VIM     VARIABLES INTERNES EN T-
!               1   : ENDOMMAGEMENT (D)
!               2   : INDICATEUR DISSIPATIF (1) OU ELASTIQUE (0)
! VAR VIP     VARIABLES INTERNES EN T+
!              IN  ESTIMATION (ITERATION PRECEDENTE)
!              OUT CALCULEES
! OUT SIGP    CONTRAINTES EN T+
! OUT DSIDEP  MATRICE TANGENTE
! OUT IRET    CODE RETOUR (0 = OK)
! ======================================================================
    integer(kind=8) :: nbmat, typedp, ndt, ndi, nvi
    real(kind=8) :: td, tf, tr
    parameter(nbmat=5)
    real(kind=8) :: materf(nbmat, 2), deps(6)
    character(len=8) :: mod
! ======================================================================
    common/tdim/ndt, ndi
!
! - Get temperatures
!
    call get_varc(fami, 1, 1, 'T', &
                  td, tf, tr)
! ======================================================================
! --- RECUPERATION DU TYPE DE LOI DE COMPORTEMENT DP -------------------
! ======================================================================
    mod = typmod(1)
    call dpmate(mod, imate, materf, ndt, ndi, &
                nvi, typedp)
! ======================================================================
! --- RETRAIT DE LA DEFORMATION DUE A LA DILATATION THERMIQUE ----------
! ======================================================================
    call dpvpdi(nbmat, materf, td, tf, tr, &
                depsm, deps)
! ======================================================================
    if (typedp .eq. 1) then
! ======================================================================
! --- CAS LINEAIRE -----------------------------------------------------
! ======================================================================
        call lcdpli(mod, nvi, option, materf, sigm, &
                    deps, vim, vip, sig, dsidep, &
                    iret)
    else if (typedp .eq. 2) then
! ======================================================================
! --- CAS PARABOLIQUE --------------------------------------------------
! ======================================================================
        call lcdppa(mod, nvi, option, materf, &
                    sigm, deps, vim, vip, sig, &
                    dsidep, iret)
    end if
! ======================================================================
end subroutine
