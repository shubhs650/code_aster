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
subroutine lc0059(BEHInteg, &
                  fami, kpg, ksp, jvMaterCode, &
                  compor, carcri, instam, instap, neps, epsm, &
                  deps, nsig, sigm, nvi, vim, option, &
                  sigp, vip, &
                  typmod, dsidep, codret)
!
    use Behaviour_type
    implicit none
!
#include "asterfort/Behaviour_type.h"
#include "asterfort/plasti.h"
#include "asterfort/srcomp.h"
#include "asterfort/utlcal.h"
!
    type(Behaviour_Integ), intent(in) :: BEHInteg
    character(len=*), intent(in) :: fami
    integer(kind=8), intent(in) :: kpg
    integer(kind=8), intent(in) :: ksp
    integer(kind=8), intent(in) :: jvMaterCode
    character(len=16), intent(in) :: compor(COMPOR_SIZE), option
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    real(kind=8), intent(in) :: instam
    real(kind=8), intent(in) :: instap
    integer(kind=8), intent(in) :: neps
    integer(kind=8), intent(in) :: nsig
    real(kind=8), intent(in) :: epsm(neps)
    real(kind=8), intent(in) :: deps(neps)
    real(kind=8), intent(in) :: sigm(nsig)
    integer(kind=8), intent(in) :: nvi
    real(kind=8), intent(in) :: vim(nvi)
    real(kind=8), intent(out) :: sigp(nsig)
    real(kind=8), intent(out) :: vip(nvi)
    character(len=8), intent(in) :: typmod(2)
    real(kind=8), intent(out) :: dsidep(6, 6)
    integer(kind=8), intent(out) :: codret
!
! --------------------------------------------------------------------------------------------------
!
! Behaviour
!
! LKR
!
! --------------------------------------------------------------------------------------------------
!
! In  BEHInteg       : parameters for integration of behaviour
!
! VARIABLES INTERNES DU MODELE :
!         1.  RXIP      : VARIABLE D ECROUISSAGE MECA. PLASTIQUE
!         2.  RGAMMAP   : DISTORSION PLASTIQUE
!         3.  RXIVP     : VARIABLE D ECROUISSAGE DU MECANISME
!                         VISCOPLASTIQUE
!         4.  RGAMMAVP  : DISTORSION VISCOPLASTIQUE
!         5.  RINDICDIL : INDICATEUR DE DILATANCE : 1 SI DIL    0 SINON
!         6.  INDIVISC  : INDICATEUR DE VISCO.    : 1 SI VISCO. 0 SINON
!         7.  INDIPLAS  : INDICATEUR DE PLAST.    : 1 SI PLAST. 0 SINON
!         8.  RDVME     : DEF. VOL. ELAS. MECA.
!         9.  RDVTE     : DEF. VOL. ELAS. THER.
!         10. RDVPL     : DEF. VOL. PLAS.
!         11. RDVVP     : DEF. VOL. VISCOPLAS.
!         12. DOMAINES  : DOMAINE EN FONCTION DES VALEURS DE XIP
!                          ELAS                      ---> DOMAINE = 0
!                          PLAS PRE-PIC              ---> DOMAINE = 1
!                          PLAS POST-PIC AVT CLIVAGE ---> DOMAINE = 2
!                          PLAS POST-PIC AP CLIVAGE  ---> DOMAINE = 3
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: algoInte
!
! --------------------------------------------------------------------------------------------------
!
    call utlcal('VALE_NOM', algoInte, carcri(6))
    if ((algoInte(1:10) .eq. 'SPECIFIQUE') .or. (option(1:14) .eq. 'RIGI_MECA_TANG')) then
        call srcomp(typmod, jvMaterCode, instam, instap, deps, sigm, vim, &
                    option, sigp, vip, dsidep, codret, nvi)
    else
        call plasti(BEHInteg, &
                    option, typmod, &
                    fami, kpg, ksp, jvMaterCode, &
                    compor, carcri, instam, instap, &
                    epsm, deps, &
                    sigm, &
                    nvi, vim, &
                    sigp, vip, &
                    dsidep, codret)
    end if
end subroutine
