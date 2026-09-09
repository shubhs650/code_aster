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

subroutine te0080(option, nomte)
    use FE_topo_module
    use FE_quadrature_module
    use FE_basis_module
    use FE_rhs_module
!
    implicit none
#include "FE_module.h"
#include "asterfort/fointe_varc.h"
#include "asterfort/jevech.h"
#include "asterfort/assert.h"
#include "asterfort/writeVector.h"
#include "jeveux.h"
!
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DES CHARGEMENTS ELEMENTAIRES
!                          OPTION : 'CHAR_THER_SOUR'
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
    integer(kind=8) ::  kp, itemps, isour, icode
    integer(kind=8) :: nbDof, nbres
    parameter(nbres=4)
    character(len=8) :: nompar(nbres)
    real(kind=8) :: valpar(nbres), sour, valQP(MAX_QP)
    real(kind=8) :: theta, soun, sounp1, rhs(MAX_BS_CG)
!
    call FECell%init()
    call FEQuadCell%initCell(FECell, "RIGI")
    call FEBasis%initCell(FECell)
    nbDof = FEBasis%size
!
    valQP = 0.d0
!
    if (option .eq. 'CHAR_THER_SOUR_R') then
        call jevech('PSOURCR', 'L', isour)
        do kp = 1, FEQuadCell%nbQuadPoints
            valQP(kp) = zr(isour+kp-1)
        end do
    elseif (option .eq. 'CHAR_THER_SOUR_F') then
        call jevech('PINSTR', 'L', itemps)
        call jevech('PSOURCF', 'L', isour)
!
        theta = zr(itemps+2)
        nompar(1:3) = ['X', 'Y', 'Z']
        nompar(4) = 'INST'
!
        do kp = 1, FEQuadCell%nbQuadPoints
            valpar(1:3) = FEQuadCell%points(1:3, kp)
            valpar(4) = zr(itemps)
            call fointe_varc('FM', 'RIGI', 1, 1, '+', &
                             zk8(isour), 4, nompar, valpar, &
                             sounp1, icode)
            if (theta .ne. 1.0d0) then
                valpar(4) = zr(itemps)-zr(itemps+1)
                call fointe_varc('FM', 'RIGI', 1, 1, '+', &
                                 zk8(isour), 4, nompar, valpar, &
                                 soun, icode)
            else
                soun = 0.d0
            end if
!
            if (theta < -0.5) then
                sour = sounp1
            else
                sour = theta*sounp1+(1.0d0-theta)*soun
            end if

            valQP(kp) = sour
        end do
    else
        ASSERT(ASTER_FALSE)
    end if
!
    call FeMakeRhsScal(FEQuadCell, FEBasis, valQP, rhs)
!
    call writeVector("PVECTTR", nbDof, rhs)
!
end subroutine
