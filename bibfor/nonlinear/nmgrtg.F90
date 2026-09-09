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
!
subroutine nmgrtg(FEBasis, coorpg, weight, BGSEval, &
                  lVect, lMatr, lMatrPred, &
                  fPrev, fCurr, dsidep, sigmPrev, &
                  sigmCurr, matsym, matuu, vectu)
!
    use FE_basis_module
    use FE_stiffness_module
    use FE_mechanics_module
!
    implicit none
!
#include "asterf_types.h"
#include "FE_module.h"
!
    type(FE_basis), intent(in) :: FEBasis
    real(kind=8), intent(in) :: dsidep(6, 6), weight, coorpg(3), BGSEval(3, MAX_BS_CG)
    real(kind=8), intent(in) :: sigmCurr(6), sigmPrev(6), fPrev(3, 3), fCurr(3, 3)
    real(kind=8), intent(inout) :: matuu(*), vectu(*)
    aster_logical, intent(in) :: matsym, lVect, lMatr, lMatrPred
!
! --------------------------------------------------------------------------------------------------
!
!     BUT:  CALCUL DE LA MATRICE TANGENTE EN CONFIGURATION LAGRANGIENNE
!           OPTIONS RIGI_MECA_TANG ET FULL_MECA
!
! --------------------------------------------------------------------------------------------------
!
! IN  NNO     : NOMBRE DE NOEUDS DE L'ELEMENT
! IN  NDIM    : DIMENSION DU PB
! IN  POIDS   : POIDS DES POINTS DE GAUSS
! IN  KPG     : NUMERO DU POINT DE GAUSS
! IN  VFF     : VALEUR  DES FONCTIONS DE FORME
! IN  DEF     : PRODUIT DE F PAR LA DERIVEE DES FONCTIONS DE FORME
! IN  PFF     : PRODUIT DES FONCTIONS DE FORME
! IN  OPTION  : OPTION DE CALCUL
! IN  AXI     : .TRUE. SI AXIS
! IN  R       : RAYON DU POINT DE GAUSS COURANT (EN AXI)
! IN  DSIDEP  : OPERATEUR TANGENT ISSU DU COMPORTEMENT
! IN  SIGN    : CONTRAINTES PK2 A L'INSTANT PRECEDENT (AVEC RAC2)
! IN  SIGMA   : CONTRAINTES PK2 A L'INSTANT ACTUEL    (AVEC RAC2)
! IN  MATSYM  : VRAI SI LA MATRICE DE RIGIDITE EST SYMETRIQUE
! OUT MATUU   : MATRICE DE RIGIDITE PROFIL (RIGI_MECA_TANG ET FULL_MECA)
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) ::  pff(6, MAX_BS_CG, MAX_BS_CG), def(6, MAX_BS_CG, 3)
!
! --------------------------------------------------------------------------------------------------
!
! ----- Compute internal forces vector and rigidity matrix
    if (lVect) then
        call FEMatFB(FEBasis, coorpg, BGSEval, fCurr, def)
    else
        call FEMatFB(FEBasis, coorpg, BGSEval, fPrev, def)
    end if
! ----- Rigidity matrix
    if (lMatr) then
        call FEStiffJacoVectSymAdd(FEBasis, def, weight, dsidep, matsym, matuu)
        call FEMatBB(FEBasis, BGSEval, pff)
        if (lMatrPred) then
            call FEStiffGeomVectSymAdd(FEBasis, pff, weight, sigmPrev, matsym, &
                                       matuu)
        else
            call FEStiffGeomVectSymAdd(FEBasis, pff, weight, sigmCurr, matsym, &
                                       matuu)
        end if
    end if
! ----- Internal forces
    if (lVect) then
        call FEStiffResiVectSymAdd(FEBasis, def, weight, sigmCurr, vectu)
    end if
!
end subroutine
