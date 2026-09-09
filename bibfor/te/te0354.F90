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
subroutine te0354(option, nomte)
!
    use FE_topo_module
    use FE_quadrature_module
    use FE_basis_module
    use FE_rhs_module
    use FE_eval_module
    use FE_mass_module
!
    implicit none
#include "asterfort/jevech.h"
#include "asterfort/fointe.h"
#include "asterfort/foderi.h"
#include "asterfort/writeVector.h"
#include "asterfort/writeMatrix.h"
#include "FE_module.h"
#include "jeveux.h"
!
    character(len=16) :: option, nomte
!
! ----------------------------------------------------------------------
!                   SOURCE THERMIQUE NON LINEAIRE
!
!    - FONCTION REALISEE:  CALCUL DES MATRICES ELEMENTAIRES
!                          OPTION : 'CHAR_THER_SOURNL'
!                                   'RESI_THER_SOURNL'
!                                   'MTAN_THER_SOURNL'
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
    type(FE_Cell) :: FECell
    type(FE_Quadrature) :: FEQuadCell
    type(FE_basis) :: FEBasis
!
    integer(kind=8) :: kp, iret, itemps, isour
    real(kind=8) :: valQP(MAX_QP), tg, sour, theta, coefop, dsdt
    real(kind=8) :: resi(MAX_BS_CG), mass(MAX_BS_CG, MAX_BS_CG)
    real(kind=8), pointer :: tempi(:) => null()
    aster_logical :: l_resi
!
    call FECell%init()
    call FEQuadCell%initCell(FECell, "RIGI")
    call FEBasis%initCell(FECell)
!
    call jevech('PINSTR', 'L', itemps)
    call jevech('PSOURNL', 'L', isour)
    if (zk8(isour) (1:7) .eq. '&FOZERO') goto 999
    theta = zr(itemps+2)
!
!    LECTURE DES PARAMETRES SPECIFIQUES A CHAQUE OPTION
    if (option(1:4) .eq. 'CHAR') then
        l_resi = .true.
        coefop = 1-theta
        call jevech('PTEMPER', 'L', vr=tempi)
    else if (option(1:4) .eq. 'RESI') then
        l_resi = .true.
        coefop = -1.d0
        call jevech('PTEMPEI', 'L', vr=tempi)
    else
        l_resi = .false.
        coefop = -1.d0
        call jevech('PTEMPEI', 'L', vr=tempi)
    end if
!
    valQP = 0.0
    do kp = 1, FEQuadCell%nbQuadPoints
        tg = FEEvalFuncRScal(FEBasis, tempi, FEQuadCell%points_param(1:3, kp))
!
        if (l_resi) then
            call fointe('FM', zk8(isour), 1, ['TEMP'], [tg], &
                        sour, iret)

            if (theta < -0.5) then
                valQP(kp) = sour
            else
                valQP(kp) = sour*coefop
            end if
        else
            call foderi(zk8(isour), tg, sour, dsdt)
            valQP(kp) = dsdt*coefop
        end if
    end do
!
    if (option(1:4) .eq. 'CHAR') then
        call FeMakeRhsScal(FEQuadCell, FEBasis, valQP, resi)
        call writeVector("PVECTTR", FEBasis%size, resi)
    else if (option(1:4) .eq. 'RESI') then
        call FeMakeRhsScal(FEQuadCell, FEBasis, valQP, resi)
        call writeVector("PRESIDU", FEBasis%size, resi)
    else
        call FEMassMatScal(FEQuadCell, FEBasis, mass, valQP)
        call writeMatrix("PMATTTR", FEBasis%size, FEBasis%size, ASTER_TRUE, mass)
    end if
!
999 continue
end subroutine
