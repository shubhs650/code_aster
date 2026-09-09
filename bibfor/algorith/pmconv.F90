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
subroutine pmconv(resi, resiInit, resiEval, &
                  ds_conv, &
                  timeCurr, iterNewt, &
                  coefAdim, sigmCurr, &
                  conver, lIterNewtMaxi)
!
    use NonLin_Datastructure_type
    implicit none
!
#include "asterc/r8prem.h"
#include "asterf_types.h"
#include "asterfort/GetResi.h"
#include "asterfort/pmimpr.h"
#include "asterfort/utmess.h"
#include "blas/dcopy.h"
#include "blas/dscal.h"
!
    real(kind=8), intent(in) :: resi(12), resiInit(12)
    real(kind=8), intent(inout) :: resiEval(12)
    type(NL_DS_Conv), intent(in) :: ds_conv
    real(kind=8), intent(in) :: timeCurr
    integer(kind=8), intent(in) :: iterNewt
    real(kind=8), intent(in) :: coefAdim, sigmCurr(6)
    aster_logical, intent(out) :: conver, lIterNewtMaxi
!
! --------------------------------------------------------------------------------------------------
!
! SIMU_POINT_MAT
!
! Management of convergence
!
! --------------------------------------------------------------------------------------------------
!
! IN   R      : RESIDU ACTUEL
! IN   RINI   : RESIDU INITIAL
! IN/OUT R1   : RESIDU PREMIERE ITERATION
! IN   INST   : INSTANT ACTUEL
! IN   SIGP   : CONTRAINTES ACTUELLES (POUR CONSTRUIRE LE DENOMINATEUR)
! IN   COEF   : COEF POUR ADIMENSIONNALISER LE PB
! IN   ITER   : NUMERO D'ITERATION
! In  ds_conv          : datastructure for convergence management
! OUT  ITEMAX : .TRUE. SI ITERATION MAXIMUM ATTEINTE
! OUT  CONVER : .TRUE. SI CONVERGENCE REALISEE
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: l_rela
    integer(kind=8) :: prtLevel, i, iterGlobMaxi
    real(kind=8) :: resiGlobRela, resiGlobMaxi
    real(kind=8) :: ee, e1, e2, toler, e1ini, e2ini, er1, eini
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------------------------------------
!
    conver = ASTER_FALSE
    lIterNewtMaxi = ASTER_FALSE

! - Get parameters
    call GetResi(ds_conv, type='RESI_GLOB_RELA', user_para_=resiGlobRela, l_resi_test_=l_rela)
    call GetResi(ds_conv, type='RESI_GLOB_MAXI', user_para_=resiGlobMaxi)
    iterGlobMaxi = ds_conv%iter_glob_maxi
!
    e1 = 0.d0
    e2 = 0.d0
    e1ini = 0.d0
    e2ini = 0.d0
    er1 = 0.d0

    if (iterNewt .eq. 1) then
        b_n = to_blas_int(6)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, sigmCurr, b_incx, resiEval(1), b_incy)
        b_n = to_blas_int(6)
        b_incx = to_blas_int(1)
        call dscal(b_n, 1.d0/coefAdim, resiEval(1), b_incx)
        b_n = to_blas_int(6)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, resi(7), b_incx, resiEval(7), b_incy)
        do i = 1, 12
            er1 = max(er1, abs(resiEval(i)))
        end do
        if (er1 .le. r8prem()) then
            ee = er1
            prtLevel = 4
            conver = ASTER_TRUE
            goto 999
        end if
    end if
!
    do i = 1, 6
        e1 = max(e1, abs(resi(i)))
        e1ini = max(e1ini, abs(resiInit(i)))
        e1ini = max(e1ini, abs(resiEval(i)))
    end do
    do i = 7, 12
        e2 = max(e2, abs(resi(i)))
        e2ini = max(e2ini, abs(resiInit(i)))
        e2ini = max(e2ini, abs(resiEval(i)))
    end do
    eini = max(e1ini, e2ini)

!   TEST RELATIF OU ABSOLU
    if (l_rela) then
        toler = resiGlobRela
        if (eini .gt. r8prem()) then
            e1 = e1/eini
            e2 = e2/eini
            ee = max(e1, e2)
            prtLevel = 3
        end if
    else
        toler = resiGlobMaxi
        ee = max(e1, e2)
        prtLevel = 4
    end if

    if (iterNewt .lt. iterGlobMaxi) then
        if (ee .gt. toler) then
            conver = ASTER_FALSE
        else
            conver = ASTER_TRUE
        end if
    else
        conver = ASTER_FALSE
        lIterNewtMaxi = ASTER_TRUE
        call utmess('I', 'COMPOR2_5')
    end if
!
999 continue
!
    call pmimpr(prtLevel, &
                timeCurr, iterNewt, &
                ee_=ee, eini_=eini)
end subroutine
