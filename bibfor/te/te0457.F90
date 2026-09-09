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
subroutine te0457(option, nomte)
!
    use HHO_type
    use HHO_basis_module
    use HHO_size_module
    use HHO_quadrature_module
    use HHO_Neumann_module
    use HHO_init_module, only: hhoInfoInitFace
    use HHO_eval_module
    use HHO_utils_module
!
    implicit none
!
#include "asterc/r8t0.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/foderi.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/jevech.h"
#include "asterfort/readVector.h"
#include "asterfort/writeMatrix.h"
#include "blas/dsyr.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
!---------------------------------------------------------------------------------------------------
!
!  HHO METHODS
!     BUT: CALCUL DES MATRICES ELEMENTAIRES EN THERMIQUE
!          CORRESPONDANT A ECHANGE PAROI POUR HHO
!          (LE CHARGEMENT PEUT ETRE DONNE SOUS FORME D'UNE FONCTION)
!
!          OPTIONS : 'RIGI_THER_ECHA_R'
!                    'RIGI_THER_ECHA_F'
!                    'MTAN_THER_FLUNL'
!                    'MTAN_THER_RAYO_R'
!                    'MTAN_THER_RAYO_F'
!
!  ENTREES  ---> OPTION : OPTION DE CALCUL
!           ---> NOMTE  : NOM DU TYPE ELEMENT
!
!---------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: maxpara = 4
    real(kind=8) :: valpar(maxpara)
    character(len=8) :: nompar(maxpara)
    type(HHO_Data) :: hhoData
    type(HHO_Face) :: hhoFace
    type(HHO_Quadrature) :: hhoQuadFace
    type(HHO_basis_face) :: hhoBasisFace
    real(kind=8), dimension(MSIZE_FACE_SCAL) :: basisScalEval, temp_F_curr
    real(kind=8), dimension(MSIZE_FACE_SCAL, MSIZE_FACE_SCAL) :: lhs
    real(kind=8) :: CoeffQP_curr(MSIZE_QP_FACE), coeff, d_alpha, temp_eval_curr
    real(kind=8) :: sigmEner(MSIZE_QP_FACE), epsil(MSIZE_QP_FACE), rbid, tz0, time_curr
    integer(kind=8) :: fbs, celldim, ipg, nbpara, npg
    integer(kind=8) :: j_time, j_coefh, j_para
    blas_int :: b_incx, b_lda, b_n
!
!
! -- Get number of Gauss points
!
    call elrefe_info(fami='RIGI', npg=npg)
!
! -- Retrieve HHO informations
!
    call hhoInfoInitFace(hhoFace, hhoData, npg, hhoQuadFace)
    call hhoBasisFace%initialize(hhoFace)
!
! ---- number of dofs
!
    call hhoTherFaceDofs(hhoFace, hhoData, fbs)
!
    ASSERT(hhoQuadFace%nbQuadPoints <= MSIZE_QP_FACE)
!
    celldim = hhoFace%ndim+1
    CoeffQP_curr = 0.d0
    nompar(:) = 'XXXXXXXX'
    valpar(:) = 0.d0
!
    call jevech('PINSTR', 'L', j_time)
    time_curr = zr(j_time)
!
! ---- Which option ?
!
    if (option .eq. 'RIGI_THER_ECHA_R') then
!
! ----- Get real value COEF_H
!
        call jevech('PCOEFHR', 'L', j_coefh)
        CoeffQP_curr(1:hhoQuadFace%nbQuadPoints) = zr(j_coefh)
!
    else if (option .eq. 'RIGI_THER_ECHA_F') then
!
! ---- Get Function Parameters
!
        if (celldim == 3) then
            nbpara = 4
            nompar(1:3) = (/'X', 'Y', 'Z'/)
        else if (celldim == 2) then
            nbpara = 3
            nompar(1:2) = (/'X', 'Y'/)
        else
            ASSERT(ASTER_FALSE)
        end if
!
! ---- Time
!
        nompar(nbpara) = 'INST'
        valpar(nbpara) = time_curr
!
! ----- Evaluate the analytical function COEF_H
!
        call jevech('PCOEFHF', 'L', j_coefh)
        call hhoFuncFScalEvalQp(hhoQuadFace, zk8(j_coefh), nbpara, nompar, valpar, &
                                celldim, CoeffQP_curr)
!
    else if (option .eq. 'MTAN_THER_RAYO_R') then
!
! ----- Get real value (sigmEner, epsil, temp_inf)
!
        call jevech('PRAYONR', 'L', j_para)
        sigmEner = zr(j_para)
        epsil = zr(j_para+1)
        tz0 = r8t0()
!
        call readVector('PTEMPEI', fbs, temp_F_curr)
!
        do ipg = 1, hhoQuadFace%nbQuadPoints
            temp_eval_curr = hhoEvalScalFace( &
                             hhoBasisFace, hhoData%face_degree(), hhoQuadFace%points(1:3, ipg), &
                             temp_F_curr)
!
            CoeffQP_curr(ipg) = 4.d0*sigmEner(ipg)*epsil(ipg)*(temp_eval_curr+tz0)**3
        end do
    else if (option .eq. 'MTAN_THER_RAYO_F') then
        call jevech('PRAYONF', 'L', j_para)
!
! ---- Get Function Parameters (sigmEner, epsil, temp_inf)
!
        if (celldim == 3) then
            nbpara = 4
            nompar(1:3) = (/'X', 'Y', 'Z'/)
        else if (celldim == 2) then
            nbpara = 3
            nompar(1:2) = (/'X', 'Y'/)
        else
            ASSERT(ASTER_FALSE)
        end if
!
! ---- Time +
!
        nompar(nbpara) = 'INST'
        valpar(nbpara) = time_curr
!
! ----- Evaluate the analytical function at T+
!
        call hhoFuncFScalEvalQp(hhoQuadFace, zk8(j_para), nbpara, nompar, valpar, &
                                celldim, sigmEner)
        call hhoFuncFScalEvalQp(hhoQuadFace, zk8(j_para+1), nbpara, nompar, valpar, &
                                celldim, epsil)
!
        tz0 = r8t0()
        call readVector('PTEMPEI', fbs, temp_F_curr)
!
        do ipg = 1, hhoQuadFace%nbQuadPoints
            temp_eval_curr = hhoEvalScalFace( &
                             hhoBasisFace, hhoData%face_degree(), hhoQuadFace%points(1:3, ipg), &
                             temp_F_curr)
!
            CoeffQP_curr(ipg) = 4.d0*sigmEner(ipg)*epsil(ipg)*(temp_eval_curr+tz0)**3
        end do
!
    else if (option .eq. 'MTAN_THER_FLUXNL') then
        call jevech('PFLUXNL', 'L', j_para)
!
        call readVector('PTEMPEI', fbs, temp_F_curr)
!
        do ipg = 1, hhoQuadFace%nbQuadPoints
            temp_eval_curr = hhoEvalScalFace( &
                             hhoBasisFace, hhoData%face_degree(), hhoQuadFace%points(1:3, ipg), &
                             temp_F_curr)
            call foderi(zk8(j_para), temp_eval_curr, rbid, d_alpha)
            CoeffQP_curr(ipg) = -d_alpha
        end do
    else
        ASSERT(ASTER_FALSE)
    end if
!
! ---- Compute mass matrix
!
    lhs = 0.d0
!
    call hhoBasisFace%initialize(hhoFace)
!
! ----- Loop on quadrature point
    do ipg = 1, hhoQuadFace%nbQuadPoints
! --------- Eval basis function at the quadrature point
        call hhoBasisFace%BSEval(hhoQuadFace%points(1:3, ipg), 0, hhoData%face_degree(), &
                                 basisScalEval)
! --------  Eval massMat
        coeff = CoeffQP_curr(ipg)*hhoQuadFace%weights(ipg)
        b_n = to_blas_int(fbs)
        b_incx = to_blas_int(1)
        b_lda = to_blas_int(MSIZE_FACE_SCAL)
        call dsyr('U', b_n, coeff, basisScalEval, b_incx, &
                  lhs, b_lda)
    end do
!
! ----- Copy the lower part
!
    call hhoCopySymPartMat('U', lhs(1:fbs, 1:fbs))
!
! ---- save result
!
    call writeMatrix('PMATTTR', fbs, fbs, ASTER_TRUE, lhs)
!
end subroutine
