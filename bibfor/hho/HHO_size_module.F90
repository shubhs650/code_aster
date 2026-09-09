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
module HHO_size_module
!
    use HHO_type
!
    implicit none
!
    private
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/binomial.h"
!
! --------------------------------------------------------------------------------------------------
!
! HHO
!
! Define size function that are shared by the different modules
! Include astefort/HHO_size_module.h to define the sizes of different table
!
! --------------------------------------------------------------------------------------------------
!
!
    public :: hhoTherDofs, hhoTherNLDofs, hhoMecaDofs, hhoMecaNLDofs, hhoMecaGradDofs
    public :: hhoTherFaceDofs, hhoMecaFaceDofs, hhoTherCellDofs, hhoMecaCellDofs
    private :: hhoMecaFaceDofs2, hhoTherFaceDofs2
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoTherDofs(hhoCell, hhoData, cbs, fbs, total_dofs)
!
        implicit none
!
        type(HHO_Cell), intent(in)          :: hhoCell
        type(HHO_Data), intent(in)          :: hhoData
        integer(kind=8), intent(out)                :: cbs
        integer(kind=8), intent(out)                :: fbs
        integer(kind=8), intent(out)                :: total_dofs
!
! --------------------------------------------------------------------------------------------------
!   HHO - thermics
!
!   Compute the number of dofs for thermics
!   In hhoCell      : the current HHO Cell
!   In hhoDta       : information on HHO methods
!   Out cbs         : number of cell dofs
!   Out fbs         : number of face dofs
!   Out total_dofs  : number of total dofs
!
! --------------------------------------------------------------------------------------------------
!
! ---- number of dofs
        call hhoTherCellDofs(hhoCell, hhoData, cbs)
        call hhoTherFaceDofs2(hhoCell, hhoData, fbs)
        total_dofs = cbs+hhoCell%nbfaces*fbs
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoTherNLDofs(hhoCell, hhoData, cbs, fbs, total_dofs, gbs)
!
        implicit none
!
        type(HHO_Cell), intent(in)          :: hhoCell
        type(HHO_Data), intent(in)          :: hhoData
        integer(kind=8), intent(out)                :: cbs
        integer(kind=8), intent(out)                :: fbs
        integer(kind=8), intent(out)                :: total_dofs
        integer(kind=8), intent(out)                :: gbs
!
! --------------------------------------------------------------------------------------------------
!   HHO - thermics
!
!   Compute the number of dofs for non-linear thermics
!   In hhoCell      : the current HHO Cell
!   In hhoDta       : information on HHO methods
!   Out cbs         : number of cell dofs
!   Out fbs         : number of face dofs
!   Out total_dofs  : number of total dofs
!   Out gbs         : number of gradient dofs
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: ndim
! --------------------------------------------------------------------------------------------------
!
        ndim = hhoCell%ndim
!
! ---- number of dofs
!
        call hhoTherDofs(hhoCell, hhoData, cbs, fbs, total_dofs)
        gbs = ndim*binomial(hhoData%grad_degree()+ndim, hhoData%grad_degree())
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoMecaDofs(hhoCell, hhoData, cbs, fbs, total_dofs)
!
        implicit none
!
        type(HHO_Cell), intent(in)          :: hhoCell
        type(HHO_Data), intent(in)          :: hhoData
        integer(kind=8), intent(out)                :: cbs
        integer(kind=8), intent(out)                :: fbs
        integer(kind=8), intent(out)                :: total_dofs
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the number of dofs for mechanics
!   In hhoCell      : the current HHO Cell
!   In hhoDta       : information on HHO methods
!   Out cbs         : number of cell dofs
!   Out fbs         : number of face dofs
!   Out total_dofs  : number of total dofs
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: ndim
! --------------------------------------------------------------------------------------------------
!
        ndim = hhoCell%ndim
!
! ---- number of dofs
!
        call hhoMecaCellDofs(hhoCell, hhoData, cbs)
        call hhoMecaFaceDofs2(hhoCell, hhoData, fbs)
        total_dofs = cbs+hhoCell%nbfaces*fbs
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoMecaNLDofs(hhoCell, hhoData, cbs, fbs, total_dofs, gbs, gbs_sym, gbs_axis)
!
        implicit none
!
        type(HHO_Cell), intent(in)          :: hhoCell
        type(HHO_Data), intent(in)          :: hhoData
        integer(kind=8), intent(out)        :: cbs
        integer(kind=8), intent(out)        :: fbs
        integer(kind=8), intent(out)        :: total_dofs
        integer(kind=8), intent(out)        :: gbs
        integer(kind=8), intent(out)        :: gbs_sym
        integer(kind=8), intent(out)        :: gbs_axis
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the number of dofs for non-linear mechanics
!   In hhoCell      : the current HHO Cell
!   In hhoDta       : information on HHO methods
!   Out cbs         : number of cell dofs
!   Out fbs         : number of face dofs
!   Out total_dofs  : number of total dofs
!   Out gbs         : number of gradient dofs
!   Out gbs_sym     : number of symmetric gradient dofs
!   Out gbs_axis     : number of axis gradient dofs
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: ndim, gbs_comp
! --------------------------------------------------------------------------------------------------
!
        ndim = hhoCell%ndim
        gbs_comp = binomial(hhoData%grad_degree()+ndim, hhoData%grad_degree())
! ---- number of dofs
        call hhoMecaDofs(hhoCell, hhoData, cbs, fbs, total_dofs)
        call hhoMecaGradDofs(hhoCell, hhoData, gbs, gbs_sym, gbs_axis)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoTherFaceDofs(hhoFace, hhoData, fbs)
!
        implicit none
!
        type(HHO_Face), intent(in)  :: hhoFace
        type(HHO_Data), intent(in)  :: hhoData
        integer(kind=8), intent(out)        :: fbs
!
! --------------------------------------------------------------------------------------------------
!   HHO - thermic
!
!   Compute the number of dofs for thermic
!   In hhoFace      : the current HHO Face
!   In hhoData      : information on HHO methods
!   Out fbs         : number of face dofs
!
! --------------------------------------------------------------------------------------------------
!
        fbs = binomial(hhoData%face_degree()+hhoFace%ndim, hhoData%face_degree())
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoTherFaceDofs2(hhoCell, hhoData, fbs)
!
        implicit none
!
        type(HHO_Cell), intent(in)  :: hhoCell
        type(HHO_Data), intent(in)  :: hhoData
        integer(kind=8), intent(out)        :: fbs
!
! --------------------------------------------------------------------------------------------------
!   HHO - thermic
!
!   Compute the number of dofs for thermic
!   In hhoCell      : the current HHO cell
!   In hhoData      : information on HHO methods
!   Out fbs         : number of face dofs
!
! --------------------------------------------------------------------------------------------------
!
        fbs = binomial(hhoData%face_degree()+hhoCell%ndim-1, hhoData%face_degree())
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoMecaFaceDofs(hhoFace, hhoData, fbs)
!
        implicit none
!
        type(HHO_Face), intent(in)  :: hhoFace
        type(HHO_Data), intent(in)  :: hhoData
        integer(kind=8), intent(out)        :: fbs
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the number of dofs for mechanics
!   In hhoFace      : the current HHO Face
!   In hhoData      : information on HHO methods
!   Out fbs         : number of face dofs
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: fbs_ther
! --------------------------------------------------------------------------------------------------
!
        call hhoTherFaceDofs(hhoFace, hhoData, fbs_ther)
        fbs = (hhoFace%ndim+1)*fbs_ther
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoMecaFaceDofs2(hhoCell, hhoData, fbs)
!
        implicit none
!
        type(HHO_Cell), intent(in)  :: hhoCell
        type(HHO_Data), intent(in)  :: hhoData
        integer(kind=8), intent(out)        :: fbs
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the number of dofs for mechanics
!   In hhoCell      : the current HHO Face
!   In hhoData      : information on HHO methods
!   Out fbs         : number of face dofs
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: fbs_ther
! --------------------------------------------------------------------------------------------------
!
        call hhoTherFaceDofs2(hhoCell, hhoData, fbs_ther)
        fbs = hhoCell%ndim*fbs_ther
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoTherCellDofs(hhoCell, hhoData, cbs)
!
        implicit none
!
        type(HHO_Cell), intent(in)  :: hhoCell
        type(HHO_Data), intent(in)  :: hhoData
        integer(kind=8), intent(out)        :: cbs
!
! --------------------------------------------------------------------------------------------------
!   HHO - thermic
!
!   Compute the number of dofs for thermic
!   In hhoCell      : the current HHO Cell
!   In hhoData      : information on HHO methods
!   Out cbs         : number of cell dofs
!
! --------------------------------------------------------------------------------------------------
!
        cbs = binomial(hhoData%cell_degree()+hhoCell%ndim, hhoData%cell_degree())
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoMecaCellDofs(hhoCell, hhoData, cbs)
!
        implicit none
!
        type(HHO_Cell), intent(in)  :: hhoCell
        type(HHO_Data), intent(in)  :: hhoData
        integer(kind=8), intent(out)        :: cbs
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the number of dofs for mechanics
!   In hhoCell      : the current HHO Cell
!   In hhoData      : information on HHO methods
!   Out cbs         : number of cell dofs
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: cbs_ther
! --------------------------------------------------------------------------------------------------
!
        call hhoTherCellDofs(hhoCell, hhoData, cbs_ther)
        cbs = hhoCell%ndim*cbs_ther
!
    end subroutine
!
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoMecaGradDofs(hhoCell, hhoData, gbs, gbs_sym, gbs_axis)
!
        implicit none
!
        type(HHO_Cell), intent(in)  :: hhoCell
        type(HHO_Data), intent(in)  :: hhoData
        integer(kind=8), intent(out) :: gbs
        integer(kind=8), intent(out) :: gbs_sym
        integer(kind=8), intent(out) :: gbs_axis
!
! --------------------------------------------------------------------------------------------------
!   HHO - mechanics
!
!   Compute the number of dofs for mechanics
!   In hhoCell      : the current HHO Cell
!   In hhoData      : information on HHO methods
!   Out gbs         : number of grad dofs
!   Out gbs_sym     : number of symmetric grad dofs
!   Out gbs_axis     : number of axis gradient dofs
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: gbs_comp, ndim
! --------------------------------------------------------------------------------------------------
!
        ndim = hhoCell%ndim
        gbs_comp = binomial(hhoData%grad_degree()+ndim, hhoData%grad_degree())
!
        gbs = ndim*ndim*gbs_comp
!
        if (ndim == 3) then
            gbs_sym = 6*gbs_comp
        else if (ndim == 2) then
            gbs_sym = 3*gbs_comp
        else
            ASSERT(ASTER_FALSE)
        end if
!
        if (hhoCell%l_axis) then
            gbs_axis = gbs_comp
        else
            gbs_axis = 0
        end if
        gbs = gbs+gbs_axis
        gbs_sym = gbs_sym+gbs_axis
!
    end subroutine
!
end module
