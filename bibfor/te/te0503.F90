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
subroutine te0503(option, nomte)
!
    use HHO_algebra_module
    use HHO_basis_module
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
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/nbsigm.h"
#include "asterfort/sigtmc.h"
#include "asterfort/tecach.h"
#include "asterfort/writeVector.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: HHO - External state variable
!
! Option: CHAR_MECA_TEMP_R / CHAR_MECA_HYDR_R / CHAR_MECA_PTOT_R / CHAR_MECA_SECH_R
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
    integer(kind=8) :: cbs, fbs, total_dofs, npg, kpg, gbs, gbs_sym, cbs_cmp, gbs_axis
    integer(kind=8) :: nbsig, indxVarcStrain, jvMaterc, jvInstr, faces_dofs, gbs_cmp, iret
    real(kind=8) :: time, Cauchy_curr(6), weight, coorpg(3)
    real(kind=8) :: rhs(MSIZE_TDOFS_VEC), BSCEval(MSIZE_CELL_SCAL)
    real(kind=8), dimension(MSIZE_CELL_MAT) :: bT
    real(kind=8) :: sigmVarc(6*MSIZE_QP_CELL)
    aster_logical :: l_axis
    type(Material_Para) :: materPara
!
! --------------------------------------------------------------------------------------------------
!

! - Get HHO informations
    call elrefe_info(fami=fami, npg=npg)
    call hhoInfoInitCell(hhoCell, hhoData, npg, hhoQuadCellRigi)

! - Number of dofs
    call hhoMecaNLDofs(hhoCell, hhoData, cbs, fbs, total_dofs, gbs, gbs_sym, gbs_axis)
    faces_dofs = total_dofs-cbs
    gbs_cmp = (gbs-gbs_axis)/(hhoCell%ndim*hhoCell%ndim)
    cbs_cmp = cbs/hhoCell%ndim
    nbsig = nbsigm()

! - Type of finite element
    call hhoBasisCell%initialize(hhoCell)
    l_axis = lteatt('AXIS', 'OUI')

! - Compute Operators
    if (hhoData%precompute()) then
        call hhoReloadPreCalcMeca(hhoCell, hhoData, ASTER_FALSE, gradsym)
    else
        call hhoCalcOpMeca(hhoCell, hhoData, ASTER_FALSE, gradsym)
    end if

! - Choice of external state variable
    if (option .eq. 'CHAR_MECA_TEMP_R') then
        indxVarcStrain = VARC_STRAIN_TEMP
    elseif (option .eq. 'CHAR_MECA_HYDR_R') then
        indxVarcStrain = VARC_STRAIN_HYDR
    elseif (option .eq. 'CHAR_MECA_SECH_R') then
        indxVarcStrain = VARC_STRAIN_SECH
    elseif (option .eq. 'CHAR_MECA_PTOT_R') then
        indxVarcStrain = VARC_STRAIN_PTOT
        ASSERT(hhoCell%ndim == 2)
    else
        ASSERT(ASTER_FALSE)
    end if

! - Get current time
    time = r8vide()
    call tecach('ONO', 'PINSTR', 'L', iret, iad=jvInstr)
    if (jvInstr .ne. 0) then
        time = zr(jvInstr)
    end if

! - Get material parameters
    call jevech('PMATERC', 'L', jvMaterc)

! - Initializations of material parameters on current cell
    call initParaCell(fami, zi(jvMaterc), materPara)

! - Set local coordinate system from user
    call getUserLCSWithBaryCenter(hhoCell%ndim, hhoCell%barycenter, materPara%lcsPara)

! - Compute stresses from external state variables
    call sigtmc(materPara, time, &
                nbsig, npg, hhoCell%ndim, &
                indxVarcStrain, sigmVarc)

    rhs = 0.d0
    bT = 0.d0

! - Loop on quadrature point
    do kpg = 1, hhoQuadCellRigi%nbQuadPoints
        coorpg(1:3) = hhoQuadCellRigi%points(1:3, kpg)
        weight = hhoQuadCellRigi%weights(kpg)

! ----- Initializations of material parameters on current integration point
        call initParaPoin(kpg, ksp, materPara)

! ----- tranform sigm in symmetric form
        call tranfoMatToSym(hhoCell%ndim, sigmVarc((kpg-1)*nbsig+1:kpg*nbsig), Cauchy_curr)

! ----- Eval basis function at the quadrature point
        call hhoBasisCell%BSEval(coorpg(1:3), 0, hhoData%grad_degree(), BSCEval)
!
        call hhoComputeRhsSmall(hhoCell, Cauchy_curr, weight, BSCEval, gbs_cmp, bT)
    end do
    call hho_dgemv_T(1.d0, gradsym, bT, 1.d0, rhs)

! - Set output vector
    call writeVector("PVECTUR", total_dofs, rhs)
!
    call gradsym%free()
!
end subroutine
