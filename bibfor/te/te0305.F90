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
subroutine te0305(option, nomte)
!
    use FE_topo_module
    use FE_quadrature_module
    use FE_basis_module
    use FE_rhs_module
    use FE_eval_module
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/jevecd.h"
#include "asterfort/jevech.h"
#include "asterfort/tefrep.h"
#include "asterfort/writeVector.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "FE_module.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 3D/D_PLAN (skin elements)
!
! Options: CHAR_MECA_FRSU** (for skin elements only)
!          RIGI_MECA_FRSU_* (for skin elements only)
!
! --------------------------------------------------------------------------------------------------
!
    type(FE_Skin) :: FESkin
    type(FE_Quadrature) :: FEQuad
    type(FE_Basis) :: FEBasis
!
    character(len=8) :: nompar(4)
    aster_logical :: l_func, l_suiv
    integer(kind=8) :: jv_time, jv_forc, ideplm, ideplp
    integer(kind=8) :: kp, nbpar
    real(kind=8) :: rhs(MAX_BV_CG), disp_curr(MAX_BV_CG)
    real(kind=8) :: forc_pg(3, MAX_QP), valpar(4)
    blas_int :: b_incx, b_incy, b_n

!
! --------------------------------------------------------------------------------------------------
!
    l_func = (option .eq. 'CHAR_MECA_FFSU23') .or. (option .eq. 'CHAR_MECA_FFSU12') .or. &
             (option .eq. 'CHAR_MECA_FF2D3D') .or. (option .eq. 'CHAR_MECA_FF1D2D') .or. &
             (option .eq. 'RIGI_MECA_FFSU23') .or. (option .eq. 'RIGI_MECA_FFSU12')
    l_suiv = option(13:14) .eq. "SU"
!
    call FESkin%init()
    if (l_suiv) then
        call jevech('PDEPLMR', 'L', ideplm)
        call jevech('PDEPLPR', 'L', ideplp)
        b_n = to_blas_int((FESkin%ndim+1)*FESkin%nbnodes)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, zr(ideplm), b_incx, disp_curr, b_incy)
        call daxpy(b_n, 1.d0, zr(ideplp), b_incx, disp_curr, b_incy)
        call FESkin%updateCoordinates(disp_curr)
    end if
    call FEQuad%initFace(FESkin, "RIGI")
    call FEBasis%initFace(FESkin)
!
! - Input fields: for force, no node affected -> 0
!
    if (l_func) then
        ASSERT(.not. l_suiv)
        nompar(1) = 'X'
        nompar(2) = 'Y'
        if (FESkin%ndim == 2) then
            call jevecd('PFF2D3D', jv_forc, 0.d0)
            nbpar = 4
            nompar(3) = 'Z'
        else
            nbpar = 3
            call jevecd('PFF1D2D', jv_forc, 0.d0)
        end if
        nompar(nbpar) = 'INST'
        call jevech('PINSTR', 'L', jv_time)
        valpar(nbpar) = zr(jv_time)
    else
        if (FESkin%ndim == 2) then
            call tefrep(option, 'PFR2D3D', jv_forc)
        else
            call tefrep(option, 'PFR1D2D', jv_forc)
        end if
    end if
!
! - Evaluation of force at Gauss points (from nodes)
!
    forc_pg = 0.d0
    do kp = 1, FEQuad%nbQuadPoints
        if (l_func) then
            valpar(1:FESkin%ndim+1) = FEQuad%points(1:FESkin%ndim+1, kp)
            forc_pg(1:3, kp) = FEEvalFuncFVec(zk8(jv_forc), nbpar, nompar, valpar, FEBasis%ndim)
        else
            forc_pg(1:3, kp) = FEEvalFuncRVec(FEBasis, zr(jv_forc), &
                                              FEQuad%points_param(1:3, kp))
        end if
    end do
!
! - Output
!
    if (option(1:9) .eq. 'CHAR_MECA') then
        call FeMakeRhsVec(FEQuad, FEBasis, forc_pg, rhs)
        call writeVector("PVECTUR", FEBasis%ndim*FEBasis%size, rhs)
    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
