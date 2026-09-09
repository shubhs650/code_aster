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
subroutine getQuadContEleBased(elem_dime, &
                               elem_slav_code, &
                               nbPoinInte, poinInteSlav, &
                               nb_qp, coor_qp, &
                               l_axis_, nb_node_slav_, elem_slav_coor_, &
                               weight_qp_)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/mesh_pairing_type.h"
#include "asterfort/elraga.h"
#include "asterfort/mmdonf.h"
#include "asterfort/mmmjac.h"
#include "asterfort/mmnonf.h"
#include "MeshTypes_type.h"
#include "asterfort/assert.h"
#include "jeveux.h"
!
    integer(kind=8), intent(in) :: elem_dime
    character(len=8), intent(in) :: elem_slav_code
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
! Compute quadrature in slave slide by clipping (element-based approach)
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
    character(len=8) :: elga_fami
    aster_logical :: isInside
    integer(kind=8) :: model_ndim
    real(kind=8) :: gausWeightSlav(12), gausCoorSlav(2, 12)
    real(kind=8) :: tol, xi, yi, xj, yj, ai
    integer(kind=8) :: iqp, ipt, jpt, s_init, s_i, nbGaussSlav
    real(kind=8) :: shape_func(MT_NNOMAX2D), shape_dfunc(2, MT_NNOMAX2D), jacobian_sl
!
! - Hyperparameters of the numerical method
    if (elem_slav_code .eq. 'SE2') then
        elga_fami = 'FPG3'
    else if (elem_slav_code .eq. 'TR3') then
        elga_fami = 'FPG6'
    else if (elem_slav_code .eq. 'QU4') then
        elga_fami = 'FPG9'
    else if (elem_slav_code .eq. 'SE3') then
        elga_fami = 'FPG4'
    else if (elem_slav_code .eq. 'TR6') then
        elga_fami = 'FPG12'
    else if (elem_slav_code .eq. 'TR7') then
        elga_fami = 'FPG12'
    else if (elem_slav_code .eq. 'QU8') then
        elga_fami = 'FPG16'
    else if (elem_slav_code .eq. 'QU9') then
        elga_fami = 'FPG16'
    else
        ASSERT(ASTER_FALSE)
    end if
    tol = 1.0e-15
! - Initialisation of parameters
    nb_qp = 0
    coor_qp = 0
    if (present(weight_qp_)) then
        weight_qp_ = 0.d0
    end if
    model_ndim = elem_dime-1

! - Get quadrature points for the slave cell
    call elraga(elem_slav_code, elga_fami, model_ndim, nbGaussSlav, &
                gausCoorSlav, gausWeightSlav)
! - Loop over the quadrature points
    do iqp = 1, nbGaussSlav
        isInside = .True.

        s_init = 0

        do ipt = 1, nbPoinInte
            xi = poinInteSlav(1, ipt)-gausCoorSlav(1, iqp)
            yi = poinInteSlav(2, ipt)-gausCoorSlav(2, iqp)

            jpt = ipt+1
            ! i+1 modulo n
            if (jpt > nbPoinInte) jpt = 1

            xj = poinInteSlav(1, jpt)-gausCoorSlav(1, iqp)
            yj = poinInteSlav(2, jpt)-gausCoorSlav(2, iqp)

            ai = xj*yi-xi*yj
            if (abs(ai) <= tol) cycle
            if (ai > tol) then
                s_i = +1
            else
                s_i = -1
            end if

            if (s_init == 0) then
                s_init = s_i
            else if (s_i /= s_init) then
                isInside = .False.
            end if
        end do

        ! Add quadrature points if it is inside the computed intersection
        if (isInside) then
            nb_qp = nb_qp+1
            coor_qp(1:2, nb_qp) = gausCoorSlav(1:2, iqp)
            if (present(weight_qp_)) then
                ! ------------- Get shape functions and first derivative only (for perf)
                call mmnonf(elem_slav_code, &
                            gausCoorSlav(1, iqp), gausCoorSlav(2, iqp), &
                            shape_func)
                call mmdonf(elem_slav_code, &
                            gausCoorSlav(1, iqp), gausCoorSlav(2, iqp), &
                            shape_dfunc)
! ------------- Compute jacobian
                call mmmjac(l_axis_, nb_node_slav_, elem_dime, &
                            elem_slav_code, elem_slav_coor_, &
                            shape_func, shape_dfunc, jacobian_sl)
                weight_qp_(nb_qp) = jacobian_sl*gausWeightSlav(iqp)
            end if
        end if
    end do
end subroutine
