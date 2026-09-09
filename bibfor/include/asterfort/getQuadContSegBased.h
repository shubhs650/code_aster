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
#include "asterf_types.h"
#include "MeshTypes_type.h"
#include "asterfort/mesh_pairing_type.h"
!
interface
    subroutine getQuadContSegBased(elem_dime, &
                           elem_slav_code, elem_mast_code, &
                           nbPoinInte, poinInteSlav, &
                           nb_qp, coor_qp, &
                           l_axis_, nb_node_slav_, elem_slav_coor_, &
                           weight_qp_)
        integer(kind=8), intent(in) :: elem_dime
        character(len=8), intent(in) :: elem_slav_code, elem_mast_code
        integer(kind=8), intent(in) :: nbPoinInte
        real(kind=8), intent(in) :: poinInteSlav(2, MAX_NB_INTE)
        real(kind=8), intent(out) :: coor_qp(2, MAX_NB_QUAD)
        integer(kind=8), intent(out) :: nb_qp
        integer(kind=8), optional, intent(in) :: nb_node_slav_
        real(kind=8), optional, intent(in) :: elem_slav_coor_(3, MT_NNOMAX2D)
        aster_logical, optional, intent(in) :: l_axis_
        real(kind=8), optional, intent(out) :: weight_qp_(MAX_NB_QUAD)
    end subroutine getQuadContSegBased
end interface
