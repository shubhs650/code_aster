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
! aslint: disable=W1504
!
subroutine pmstab(sigmPrev, sigmCurr, epsiPrev, epsiIncr, &
                  nbVari, vim, vip, &
                  timePrev, timeCurr, iterNewt, &
                  tablName, tablType, tablNbPara, &
                  tablParaName, tablVale, &
                  lLoadGrad, valeImpo, lPrintMatr, dsidep, variName, &
                  nbVariTabl)
!
    implicit none
!
#include "asterfort/fgequi.h"
#include "asterfort/tbajli.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "blas/dscal.h"
!
    real(kind=8), intent(in) :: sigmCurr(6), epsiIncr(9)
    real(kind=8), intent(out) :: sigmPrev(6), epsiPrev(9)
    integer(kind=8), intent(in) :: nbVari
    real(kind=8), intent(out) :: vim(nbVari)
    real(kind=8), intent(in) :: vip(nbVari)
    real(kind=8), intent(out) :: timePrev
    real(kind=8), intent(in) :: timeCurr
    integer(kind=8), intent(in) :: iterNewt
    character(len=8), intent(in) :: tablName
    integer(kind=8), intent(in) :: tablType
    integer(kind=8), intent(in) :: tablNbPara
    character(len=16), dimension(:), intent(in) :: tablParaName
    real(kind=8), dimension(:), intent(inout) :: tablVale
    aster_logical, intent(in) :: lLoadGrad
    real(kind=8), intent(in) :: valeImpo(9)
    aster_logical, intent(in) :: lPrintMatr
    real(kind=8), intent(in) :: dsidep(*)
    character(len=8), intent(in) :: variName(nbVari)
!
! --------------------------------------------------------------------------------------------------
!
! SIMU_POINT_MAT
!
! Save values in table
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), parameter  :: rac2 = sqrt(2.d0)
    complex(kind=8), parameter :: c16Dummy = (0.d0, 0.d0)
    character(len=4), parameter :: epsiName(6) = (/'EPXX', 'EPYY', 'EPZZ', &
                                                   'EPXY', 'EPXZ', 'EPYZ'/)
    character(len=4), parameter :: sigmName(6) = (/'SIXX', 'SIYY', 'SIZZ', &
                                                   'SIXY', 'SIXZ', 'SIYZ'/)
    character(len=4), parameter :: gradName(9) = (/'F11', 'F12', 'F13', &
                                                   'F21', 'F22', 'F23', &
                                                   'F31', 'F32', 'F33'/)
    integer(kind=8) :: i, nbCmpEpsi, nbVariTabl
    character(len=8) :: k8b, vk8(2)
    real(kind=8) ::  epsiCurr(9), epsiTabl(9)
    real(kind=8) :: sigmEqui(17), sigmTabl(6)
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------------------------------------
!
    if (lLoadGrad) then
        nbCmpEpsi = 9
    else
        nbCmpEpsi = 6
    end if

! - Update time
    timePrev = timeCurr

! - Update strains
    if (nbCmpEpsi .eq. 6) then
        b_n = to_blas_int(nbCmpEpsi)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, epsiPrev, b_incx, epsiCurr, b_incy)
        b_n = to_blas_int(nbCmpEpsi)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0, epsiIncr, b_incx, epsiCurr, b_incy)
        b_n = to_blas_int(nbCmpEpsi)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, epsiCurr, b_incx, epsiPrev, b_incy)
    else
        b_n = to_blas_int(nbCmpEpsi)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, valeImpo, b_incx, epsiPrev, b_incy)
    end if

! - Create strains for table
    if (nbCmpEpsi .eq. 6) then
        b_n = to_blas_int(nbCmpEpsi)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, epsiCurr, b_incx, epsiTabl, b_incy)
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        call dscal(b_n, 1.d0/rac2, epsiTabl(4), b_incx)
    else
        b_n = to_blas_int(nbCmpEpsi)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, valeImpo, b_incx, epsiTabl, b_incy)
    end if

! - Update internal state variables
    b_n = to_blas_int(nbVari)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, vip, b_incx, vim, b_incy)

! - Update stress
    b_n = to_blas_int(6)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, sigmCurr, b_incx, sigmPrev, b_incy)

! - Create stresses for table
    b_n = to_blas_int(6)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, sigmCurr, b_incx, sigmTabl, b_incy)

! - Compute invariants of increment of stress
    b_n = to_blas_int(3)
    b_incx = to_blas_int(1)
    call dscal(b_n, 1.d0/rac2, sigmTabl(4), b_incx)
    call fgequi(sigmTabl, 'SIGM_DIR', 3, sigmEqui)

    if (tablType .eq. 0) then
! ----- Save time
        tablVale(1) = timeCurr

! ----- Save strains in table
        tablVale(2:2+nbCmpEpsi-1) = epsiTabl(1:nbCmpEpsi)

! ----- Save stresses in table
        tablVale(2+nbCmpEpsi:2+nbCmpEpsi+6-1) = sigmTabl(1:6)

! ----- Save invariants of stress in table
        tablVale(nbCmpEpsi+8) = sigmEqui(16)
        tablVale(nbCmpEpsi+9) = sigmEqui(1)

! ----- Save internal state variables
        tablVale(2+nbCmpEpsi+8:2+nbCmpEpsi+8+nbVariTabl-1) = vip(1:nbVariTabl)

! ----- Save iteration of Newton
        tablVale(tablNbPara) = iterNewt

! ----- Save tangent matrix
        if (lPrintMatr) then
            tablVale(2+nbCmpEpsi+8+nbVariTabl:2+nbCmpEpsi+8+nbVariTabl+36-1) = dsidep(1:36)
        end if

        call tbajli(tablName, tablNbPara, tablParaName, [0], tablVale, &
                    [c16Dummy], k8b, 0)
    else
! ----- Save time
        tablVale(1) = timeCurr

! ----- Save strains
        vk8(1) = 'EPSI'
        do i = 1, nbCmpEpsi
            tablVale(2) = epsiTabl(i)
            if (lLoadGrad) then
                vk8(2) = gradName(i)
            else
                vk8(2) = epsiName(i)
            end if
            call tbajli(tablName, tablNbPara, tablParaName, [0], tablVale, &
                        [c16Dummy], vk8, 0)
        end do

! ----- Save stresses
        vk8(1) = 'SIGM'
        do i = 1, 6
            tablVale(2) = sigmTabl(i)
            vk8(2) = sigmName(i)
            call tbajli(tablName, tablNbPara, tablParaName, [0], tablVale, &
                        [c16Dummy], vk8, 0)
        end do

! ----- Save invariants of stress in table
        vk8(1) = 'SIEQ'
        tablVale(2) = sigmEqui(1)
        vk8(2) = 'VMIS'
        call tbajli(tablName, tablNbPara, tablParaName, [0], tablVale, &
                    [c16Dummy], vk8, 0)
        vk8(2) = 'VMIS'
        tablVale(2) = sigmEqui(16)
        vk8(2) = 'TRACE'
        call tbajli(tablName, tablNbPara, tablParaName, [0], tablVale, &
                    [c16Dummy], vk8, 0)

! ----- Save internal state variables
        vk8(1) = 'VARI'
        do i = 1, nbVariTabl
            tablVale(2) = vip(i)
            vk8(2) = variName(i)
            call tbajli(tablName, tablNbPara, tablParaName, [0], tablVale, &
                        [c16Dummy], vk8, 0)
        end do
    end if
!
end subroutine
