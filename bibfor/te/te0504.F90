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
subroutine te0504(option, nomte)
!
    use HHO_basis_module
    use HHO_eval_module
    use HHO_init_module
    use HHO_matrix_module
    use HHO_Meca_module
    use HHO_quadrature_module
    use HHO_size_module
    use HHO_SmallStrainMeca_module
    use HHO_type
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterc/r8vide.h"
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/dmatmc.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/jevech.h"
#include "asterfort/nbsigm.h"
#include "asterfort/readVector.h"
#include "asterfort/sigtmc.h"
#include "asterfort/writeVector.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: HHO
!
! Option: SIEF_ELGA
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = 'RIGI'
    integer(kind=8), parameter :: ksp = 1
    type(HHO_basis_cell) :: hhoBasisCell
    type(HHO_Cell) :: hhoCell
    type(HHO_Data) :: hhoData
    type(HHO_matrix) :: gradsym
    type(HHO_Quadrature) :: hhoQuadCellRigi
    integer(kind=8) :: cbs, fbs, total_dofs, npg, kpg, gbs, gbs_sym
    integer(kind=8) :: nbsig, jvMaterc, i, j, gbs_axis
    real(kind=8) :: time, sigma(6), weight, coorpg(3)
    real(kind=8) :: E_coeff(MSIZE_CELL_MAT), Eps(6)
    real(kind=8) :: dmat(6, 6), BSCEval(MSIZE_CELL_SCAL), dmat3d(6, 6)
    real(kind=8) :: sigmVarc(6*MSIZE_QP_CELL), sief(6*MSIZE_QP_CELL)
    real(kind=8) :: depl(MSIZE_TDOFS_VEC), sigm_mat(6)
    type(Material_Para) :: materPara
!
! --------------------------------------------------------------------------------------------------
!

! - Get HHO informations
    call elrefe_info(fami=fami, npg=npg)
    call hhoInfoInitCell(hhoCell, hhoData, npg, hhoQuadCellRigi)

! - Number of dofs
    call hhoMecaNLDofs(hhoCell, hhoData, cbs, fbs, total_dofs, gbs, gbs_sym, gbs_axis)
    nbsig = nbsigm()

! - Type of finite element
    call hhoBasisCell%initialize(hhoCell)

! - Compute Operators
    if (hhoData%precompute()) then
        call hhoReloadPreCalcMeca(hhoCell, hhoData, ASTER_FALSE, gradsym)
    else
        call hhoCalcOpMeca(hhoCell, hhoData, ASTER_FALSE, gradsym)
    end if

! - Get displacements
    call readVector("PDEPLAR", total_dofs, depl)

! - Get current time
    time = r8vide()

! - Get material parameters
    call jevech('PMATERC', 'L', jvMaterc)

! - Initializations of material parameters on current cell
    call initParaCell(fami, zi(jvMaterc), materPara)

! - Set local coordinate system from user
    call getUserLCSWithBaryCenter(hhoCell%ndim, hhoCell%barycenter, materPara%lcsPara)

! - Compute stresses from external state variables
    call sigtmc(materPara, time, &
                nbsig, npg, hhoCell%ndim, &
                VARC_STRAIN_ALL, sigmVarc)

! - Compute E = gradrec_sym * depl
    call gradsym%dot(depl, E_coeff)

! - Compute total stress
    sief = 0.d0
    do kpg = 1, hhoQuadCellRigi%nbQuadPoints
        coorpg(1:3) = hhoQuadCellRigi%points(1:3, kpg)
        weight = hhoQuadCellRigi%weights(kpg)

! ----- Initializations of material parameters on current integration point
        call initParaPoin(kpg, ksp, materPara)

! ----- Eval basis function at the quadrature point
        call hhoBasisCell%BSEval(coorpg(1:3), 0, hhoData%grad_degree(), BSCEval)

! ----- Eval deformations
        Eps = hhoEvalSymMatCell(hhoCell%ndim, gbs_sym, BSCEval, E_coeff)

! ----- Compute elasticity matrix
        call dmatmc(materPara, '+', time, &
                    nbsig, dmat)
        call tranfoTensToSym(nbsig, dmat, dmat3d)
!
! ----- Compute SIGM_ELGA
        sigma = 0.d0
        do i = 1, nbsig
            do j = 1, nbsig
                sigma(i) = sigma(i)+Eps(j)*dmat3d(i, j)
            end do
        end do
        call tranfoSymToMat(hhoCell%ndim, sigma, sigm_mat)
!
! ----- Compute SIEF_ELGA
        sief((kpg-1)*nbsig+1:kpg*nbsig) = sigm_mat(1:nbsig)-sigmVarc((kpg-1)*nbsig+1:kpg*nbsig)
    end do

! - Set output vector
    call writeVector("PCONTRR", npg*nbsig, sief)
!
    call gradsym%free()
!
end subroutine
