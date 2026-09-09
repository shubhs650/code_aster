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
! aslint: disable=W0104
!
subroutine lc0025(fami, kpg, ksp, ndim, imate, &
                  instam, instap, &
                  epsm, deps, sigm, nvi, vim, option, &
                  sigp, vip, typmod, &
                  numlc, dsidep, codret)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/lcrank.h"
#include "asterfort/mctgel.h"
#include "asterfort/rcvala.h"
#include "asterfort/rcvarc.h"
#include "asterfort/Behaviour_type.h"
!
    character(len=*) :: fami
    integer(kind=8) :: kpg
    integer(kind=8) :: ksp
    integer(kind=8) :: ndim
    integer(kind=8) :: imate
    real(kind=8) :: instam
    real(kind=8) :: instap
    real(kind=8) :: epsm(6)
    real(kind=8) :: deps(6)
    real(kind=8) :: sigm(6)
    integer(kind=8), intent(in) :: nvi
    real(kind=8) :: vim(nvi)
    character(len=16) :: option
    real(kind=8) :: sigp(6)
    real(kind=8) :: vip(nvi)
    character(len=8) :: typmod(*)
    integer(kind=8) :: numlc
    real(kind=8) :: dsidep(6, 6)
    integer(kind=8) :: codret
! ======================================================================
!
    integer(kind=8) :: iret, icode(3)
    real(kind=8) :: tp, tm, tref, rprops(3), r0
    character(len=16) :: nomres(3)
    aster_logical     :: epflag
    data r0/0.0d0/
!
    common/debug/epflag
!
! Remarque: Utilise RESI_INTE = instant_post > 0.
! --------  pour activation affichage detaille
!
!
! ======================================================================
!
!     BUT: LOI DE COMPORTEMENT DE RANKINE
!
!
!       IN      FAMI    FAMILLE DE POINT DE GAUSS (RIGI,MASS,...)
!       IN      KPG,KSP NUMERO DU (SOUS)POINT DE GAUSS
!       IN      NDIM    DIMENSION DE L ESPACE (3D=3,2D=2,1D=1)
!               TYPMOD  TYPE DE MODELISATION
!               IMATE    ADRESSE DU MATERIAU CODE
!               COMPOR    COMPORTEMENT DE L ELEMENT
!                     COMPOR(1) = RELATION DE COMPORTEMENT (CHABOCHE...)
!                     COMPOR(2) = NB DE VARIABLES INTERNES
!                     COMPOR(3) = TYPE DE DEFORMATION (PETIT,JAUMANN...)
!               CRIT    CRITERES  LOCAUX
!                       CRIT(1) = NOMBRE D ITERATIONS MAXI A CONVERGENCE
!                                 (ITER_INTE_MAXI == ITECREL)
!                       CRIT(2) = TYPE DE JACOBIEN A T+DT
!                                 (TYPE_MATR_COMP == MACOMP)
!                                 0 = EN VITESSE     > SYMETRIQUE
!                                 1 = EN INCREMENTAL > NON-SYMETRIQUE
!                       CRIT(3) = VALEUR DE LA TOLERANCE DE CONVERGENCE
!                                 (RESI_INTE == RESCREL)
!                       CRIT(5) = NOMBRE D'INCREMENTS POUR LE
!                                 REDECOUPAGE LOCAL DU PAS DE TEMPS
!                                 (ITER_INTE_PAS == ITEDEC)
!                                 0 = PAS DE REDECOUPAGE
!                                 N = NOMBRE DE PALIERS
!               INSTAM   INSTANT T
!               INSTAP   INSTANT T+DT
!               EPSM   DEFORMATION TOTALE A T
!               DEPS   INCREMENT DE DEFORMATION TOTALE
!               SIGM    CONTRAINTE A T
!               VIM    VARIABLES INTERNES A T    + INDICATEUR ETAT T
!    ATTENTION  VIM    VARIABLES INTERNES A T MODIFIEES SI REDECOUPAGE
!               OPTION     OPTION DE CALCUL A FAIRE
!                             'RIGI_MECA_TANG'> DSIDEP(T)
!                             'FULL_MECA'     > DSIDEP(T+DT) , SIG(T+DT)
!                             'RAPH_MECA'     > SIG(T+DT)
!
!       OUT     SIGP    CONTRAINTE A T+DT
!               VIP    VARIABLES INTERNES A T+DT + INDICATEUR ETAT T+DT
!               DSIDEP    MATRICE DE COMPORTEMENT TANGENT A T+DT OU T
!               CODRET
!
! ======================================================================
!
!       APPEL DE RCVARC POUR LA RECUPERATION DE LA TEMPERATURE
!       RAISON: CETTE ROUTINE EST APPELEE EN THM AUSSI... (CALCME)
    call rcvarc(' ', 'TEMP', '-', fami, kpg, &
                ksp, tm, iret)
    call rcvarc(' ', 'TEMP', '+', fami, kpg, &
                ksp, tp, iret)
    call rcvarc(' ', 'TEMP', 'REF', fami, kpg, &
                ksp, tref, iret)
!
    dsidep(:, :) = r0
!
    if (option(1:14) .eq. 'RIGI_MECA_TANG') then

!       write(6,'(A)')
!       write(6,'(A)')'> LC0101 :: entering MCTGEL'
!       write(6,'(A,6(1X,E15.8))')'! * DEPS =',(deps(i),i=1,6)
!
        nomres(1) = 'E       '
        nomres(2) = 'NU      '
        call rcvala(imate, ' ', 'ELAS', 0, '   ', &
                    [tp], 2, nomres, rprops(2), icode, 2)
!
        call mctgel(dsidep, rprops)
!
    else
!
!       write(6,'(A)')
!       write(6,'(A)')'> LC0101 :: entering LCRANK'
!       write(6,'(A,6(1X,E15.8))')'! * DEPS =',(deps(i),i=1,6)
        call lcrank(ndim, typmod, imate, option, tm, tp, &
                    deps, sigm, sigp, vim, vip, dsidep, codret)
!
    end if
!
end subroutine
