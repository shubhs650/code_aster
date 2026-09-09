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
subroutine getQuadCont(parameters, elem_dime, &
                       elem_slav_code, elem_mast_code, &
                       nbPoinInte, poinInteSlav, &
                       nb_qp, coor_qp, &
                       l_axis_, nb_node_slav_, elem_slav_coor_, &
                       weight_qp_)
    use contact_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/mesh_pairing_type.h"
#include "asterfort/getQuadContSegBased.h"
#include "asterfort/getQuadContEleBased.h"
#include "MeshTypes_type.h"
#include "jeveux.h"
!
    type(ContactParameters), intent(in) :: parameters
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
!
! --------------------------------------------------------------------------------------------------
!
! Contact - Quadrature
!
! Compute quadrature in slave slide
!
! --------------------------------------------------------------------------------------------------
!
! In  modelDime        : dimension of model
! In  elem_slav_code   : code element for slave side from contact element
! In  elem_slav_coor   : coordinates from slave side of contact element
! In  nbPoinInte       : number of intersection points
! In  poinInteSlav     : coordinates of intersection points (in slave parametric space)
! Out nb_qp            : number of quadrature points
! Out coor_qp          : coordinates of quadrature points (parametric slave space)
! Out weight_qp        : weight of quadrature points
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: use_segbased
!
! --------------------------------------------------------------------------------------------------
!
    use_segbased = (parameters%inte_type .gt. 0)
    if (use_segbased) then
        call getQuadContSegBased(elem_dime, &
                                 elem_slav_code, elem_mast_code, &
                                 nbPoinInte, poinInteSlav, &
                                 nb_qp, coor_qp, &
                                 l_axis_, nb_node_slav_, elem_slav_coor_, &
                                 weight_qp_)
    else
        call getQuadContEleBased(elem_dime, &
                                 elem_slav_code, &
                                 nbPoinInte, poinInteSlav, &
                                 nb_qp, coor_qp, &
                                 l_axis_, nb_node_slav_, elem_slav_coor_, &
                                 weight_qp_)
    end if
end subroutine
