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
subroutine te0217(option, nomte)
!
    use FE_topo_module
    use FE_quadrature_module
    use FE_basis_module
    use FE_stiffness_module
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/fointe.h"
#include "asterfort/jevech.h"
#include "asterfort/rcvalb.h"
#include "asterfort/writeVector.h"
#include "FE_module.h"
!
    character(len=16) :: option, nomte
!.......................................................................
!
!     BUT: CALCUL DU SECOND MEMBRE ELEMENTAIRE EN THERMIQUE CORRESPON-
!          DANT A UN GRADIENT IMPOSE DE TEMPERATURE
!          ELEMENTS ISOPARAMETRIQUES
!
!          OPTION : 'CHAR_THER_GRAI_R/F'
!
!     ENTREES  ---> OPTION : OPTION DE CALCUL
!              ---> NOMTE  : NOM DU TYPE ELEMENT
!.......................................................................
!
    type(FE_Cell) :: FECell
    type(FE_Quadrature) :: FEQuadCell
    type(FE_basis) :: FEBasis
!
    integer(kind=8) :: icodre(1), kp, igrai, imate, ier, itemps
    character(len=8) :: nompar(4), grxf, gryf, grzf
    character(len=8), parameter :: fami = "RIGI", poum = "+"
!
    real(kind=8) :: valres(1), valpar(4), lambda
    real(kind=8) :: grx, gry, grz, time
    real(kind=8) :: load(MAX_BS_CG), valQP(3, MAX_QP)
!
    aster_logical :: fonc
!-----------------------------------------------------------------------
!
    call FECell%init()
    call FEBasis%initCell(FECell)
    call FEQuadCell%initCell(FECell, fami)
!
    call jevech('PMATERC', 'L', imate)
!
    if (option .eq. 'CHAR_THER_GRAI_R') then
        fonc = .false.
        call jevech('PGRAINR', 'L', igrai)
        grx = zr(igrai)
        gry = zr(igrai+1)
        grz = zr(igrai+2)
        time = 0.d0
    else if (option .eq. 'CHAR_THER_GRAI_F') then
        fonc = .true.
        call jevech('PINSTR', 'L', itemps)
        call jevech('PGRAINF', 'L', igrai)
        grxf = zk8(igrai)
        gryf = zk8(igrai+1)
        grzf = zk8(igrai+2)
        nompar(1) = 'X'
        nompar(2) = 'Y'
        nompar(3) = 'Z'
        nompar(4) = 'INST'
        time = zr(itemps)
        valpar(4) = time
    end if
!
    do kp = 1, FEQuadCell%nbQuadPoints
!
        call rcvalb(fami, kp, 1, poum, zi(imate), &
                    ' ', 'THER', 1, 'INST', [time], &
                    1, 'LAMBDA', valres, icodre, 1)
        lambda = valres(1)
!
        if (fonc) then
            valpar(1:3) = FEQuadCell%points_param(1:3, kp)
            call fointe('FM', grxf, 4, nompar, valpar, grx, ier)
            call fointe('FM', gryf, 4, nompar, valpar, gry, ier)
            if (FECell%ndim == 3) then
                call fointe('FM', grzf, 4, nompar, valpar, grz, ier)
            else
                grz = 0.d0
            end if
        end if
        valQP(1:3, kp) = lambda*[grx, gry, grz]
!
    end do
!
    call FEStiffResiScal(FEQuadCell, FEBasis, valQP, load)
!
    call writeVector("PVECTTR", FEBasis%size, load)
!
end subroutine
