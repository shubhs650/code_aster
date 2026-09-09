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
subroutine pidefo_hho(compor, &
                  ndim, kpg, fm, &
                  epsm, epsp, epsd, dtau, copilo)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/r8inir.h"
#include "blas/dcopy.h"
#include "blas/ddot.h"
#include "blas/dnrm2.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: compor(COMPOR_SIZE)
    integer(kind=8) :: ndim, kpg
    real(kind=8) :: epsm(6), epsp(6), epsd(6)
    real(kind=8) :: fm(3, 3)
    real(kind=8) :: dtau
    real(kind=8) :: copilo(:,:)
!
! --------------------------------------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (PILOTAGE - PRED_ELAS/DEFORMATION)
!
! PILOTAGE PAR DEFORMATION
!
! --------------------------------------------------------------------------------------------------
!
! IN  KPG    : NUMERO DU POINT DE GAUSS
! IN  FM     : GRADIENT DE LA TRANSFORMATION AU TEMPS MOINS
! IN  EPSM   : DEFORMATIONS AU TEMPS MOINS
! IN  EPSP   : CORRECTION DE DEFORMATIONS DUES AUX CHARGES FIXES
! IN  EPSD   : CORRECTION DE DEFORMATIONS DUES AUX CHARGES PILOTEES
! OUT COPILO : COEFFICIENTS A0 ET A1 POUR CHAQUE POINT DE GAUSS
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), parameter :: rac2 = sqrt(2.d0)
    aster_logical :: grand
    integer(kind=8) :: ndimsi
    real(kind=8) :: ff, em(6), epsmno
    integer(kind=8) :: ij, kl, i, j, k, l
    blas_int :: b_incx, b_incy, b_n
    integer(kind=8), parameter :: indi(6) = (/1, 2, 3, 2, 3, 3/)
    integer(kind=8), parameter :: indj(6) = (/1, 2, 3, 1, 1, 2/)
    integer(kind=8), parameter :: prac(6) = (/0, 0, 0, 1, 1, 1/)
!
! --------------------------------------------------------------------------------------------------
!
    grand = compor(DEFO) .ne. 'PETIT'
    ndimsi = 2*ndim

! - TRANSPORT DU TENSEUR DES DEFORMATIONS E := F E FT
    if (grand) then
        b_n = to_blas_int(ndimsi)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, epsm, b_incx, em, b_incy)
        call r8inir(ndimsi, 0.d0, epsm, 1)
!
        do ij = 1, ndimsi
            do kl = 1, ndimsi
                i = indi(ij)
                j = indj(ij)
                k = indi(kl)
                l = indj(kl)
                ff = (fm(i, k)*fm(j, l)+fm(i, l)*fm(j, k))/2
                ff = ff*rac2**prac(ij)*rac2**prac(kl)
                epsm(ij) = epsm(ij)+ff*em(kl)
            end do
        end do
    end if

! - INCREMENT DE DEFORMATION PROJETE
    b_n = to_blas_int(ndimsi)
    b_incx = to_blas_int(1)
    epsmno = dnrm2(b_n, epsm, b_incx)
    b_n = to_blas_int(ndimsi)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    copilo(1, kpg) = ddot(b_n, epsm, b_incx, epsp, b_incy)/epsmno
    b_n = to_blas_int(ndimsi)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    copilo(2, kpg) = ddot(b_n, epsm, b_incx, epsd, b_incy)/epsmno
!
end subroutine
