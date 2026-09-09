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
subroutine te0450(option, nomte)
!
    use HHO_basis_module
    use HHO_compor_module
    use HHO_eval_module
    use HHO_init_module, only: hhoInfoInitCell
    use HHO_Meca_module
    use HHO_quadrature_module
    use HHO_size_module
    use HHO_type
    use HHO_utils_module
    use resi_refe_module, only: RESI_REFE
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/writeVector.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: nomte, option
!
! --------------------------------------------------------------------------------------------------
!
! HHO
! Mechanics - FORC_NODA and REFE_FORC_NODA
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = "RIGI", typmod2 = "HHO"
    type(HHO_Data) :: hhoData
    type(HHO_Cell) :: hhoCell
    type(HHO_Meca_State) :: hhoMecaState
    type(HHO_Quadrature) :: hhoQuadCellRigi
    type(HHO_Compor_State) :: hhoCS
    integer(kind=8) :: cbs, fbs, total_dofs, npg, i, isig, cbs_cmp, j, fbs_cmp, ipg
    integer(kind=8) :: faces_dofs
    aster_logical :: l_largestrains
    real(kind=8) :: rhs(MSIZE_TDOFS_VEC), refe_rhs(MSIZE_TDOFS_VEC)
    real(kind=8) :: stress(6*MSIZE_QP_CELL), sigm_refe, val_refe(3)
    type(RESI_REFE):: refe
!
! --------------------------------------------------------------------------------------------------
!
    if (option /= "FORC_NODA" .and. option /= "REFE_FORC_NODA") then
        ASSERT(ASTER_FALSE)
    end if

! - Get element parameters
    call elrefe_info(fami=fami, npg=npg)

! - Get HHO data on the modelisation
    call hhoInfoInitCell(hhoCell, hhoData, npg, hhoQuadCellRigi)

! - Number of dofs
    call hhoMecaDofs(hhoCell, hhoData, cbs, fbs, total_dofs)
    faces_dofs = total_dofs-cbs

! - Type of finite element
    call hhoCS%initialize(fami, option, hhoCell%ndim, hhoCell%barycenter, typmod2)
    call hhoMecaState%initialize(hhoCell, hhoData, hhoCS)

! - Large strains ?
    l_largestrains = hhoCS%l_largestrain

! - Compute Operators
    if (hhoData%precompute()) then
        call hhoReloadPreCalcMeca(hhoCell, hhoData, l_largestrains, &
                                  hhoMecaState%grad, hhoMecaState%stab)
    else
        call hhoCalcOpMeca(hhoCell, hhoData, l_largestrains, hhoMecaState%grad, hhoMecaState%stab)
    end if
!
    if (option == "FORC_NODA") then
        call hhoLocalForcNoda(hhoCell, hhoData, hhoQuadCellRigi, hhoMecaState, &
                              hhoCS, hhoCS%sig_prev, rhs)

    elseif (option == "REFE_FORC_NODA") then
        call refe%Init(nomte)
        sigm_refe = refe%GetRef('SIGM')
        call refe%Check()

        stress = 0.d0
        refe_rhs = 0.d0
        do isig = 1, hhoCS%nbsigm
            do ipg = 1, hhoQuadCellRigi%nbQuadPoints
                stress((ipg-1)*hhoCS%nbsigm+isig) = sigm_refe
            end do
            call hhoLocalForcNoda(hhoCell, hhoData, hhoQuadCellRigi, hhoMecaState, &
                                  hhoCS, stress, rhs)
            do i = 1, total_dofs
                refe_rhs(i) = refe_rhs(i)+abs(rhs(i))
            end do
            do ipg = 1, hhoQuadCellRigi%nbQuadPoints
                stress((ipg-1)*hhoCS%nbsigm+isig) = 0.d0
            end do
        end do
        cbs_cmp = cbs/hhoCell%ndim
        fbs_cmp = fbs/hhoCell%ndim
        ! Set mean value
        do i = 1, hhoCell%ndim
            val_refe(i) = refe_rhs(faces_dofs+(i-1)*cbs_cmp+1)
            do j = 1, hhoCell%nbfaces
                val_refe(i) = max(val_refe(i), refe_rhs((j-1)*fbs+(i-1)*fbs_cmp+1))
            end do
        end do
        rhs = 0.d0
        do i = 1, hhoCell%ndim
            rhs(faces_dofs+(i-1)*cbs_cmp+1:faces_dofs+i*cbs_cmp) = val_refe(i)
            do j = 1, hhoCell%nbfaces
                rhs((j-1)*fbs+(i-1)*fbs_cmp+1:(j-1)*fbs+i*fbs_cmp) = val_refe(i)
            end do
        end do
    else
        ASSERT(ASTER_FALSE)
    end if

! - Save rhs
    call writeVector('PVECTUR', total_dofs, rhs)
!
    call hhoMecaState%free()
!
end subroutine
