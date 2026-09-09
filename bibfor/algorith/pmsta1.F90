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
subroutine pmsta1(sigmPrev, sigmCurr, epsiIncr, &
                  nbVari, nbVariTabl, &
                  vim, vip, &
                  tablNbPara, tablType, &
                  tablParaName, tablParaType, tablVale, &
                  lLoadGrad, variName, sddisc, &
                  liccvg, lIterNewtMaxi, conver, newtLoopAction)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/detrsd.h"
#include "asterfort/fgequi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/pmevdr.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbcrsd.h"
#include "asterfort/wkvect.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
!
    real(kind=8), intent(in) :: sigmPrev(6), sigmCurr(6), epsiIncr(9)
    integer(kind=8), intent(in)  :: nbVari, nbVariTabl
    real(kind=8), intent(in) :: vim(nbVari), vip(nbVari)
    integer(kind=8), intent(in) :: tablNbPara, tablType
    character(len=16), dimension(:), intent(in) :: tablParaName
    character(len=8), dimension(:), intent(in) :: tablParaType
    real(kind=8), dimension(:), intent(inout) :: tablVale
    aster_logical, intent(in) :: lLoadGrad
    character(len=8), intent(in) :: variName(nbVari)
    character(len=19), intent(in) :: sddisc
    integer(kind=8), intent(in) :: liccvg(5)
    aster_logical, intent(in) :: lIterNewtMaxi, conver
    integer(kind=8), intent(out) :: newtLoopAction
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
    character(len=24), parameter :: tablIncr = '&&OP0033.TABINC'
    character(len=8), parameter :: k8Dummy = " "
    integer(kind=8) :: i, nbCmpEpsi
    character(len=24), parameter :: variIncrJv = "&&OP0033.VARI"
    integer(kind=8) :: jvVariIncr
    character(len=8) :: vk8(2)
    real(kind=8) :: sigmIncr(6), epsiIncrTabl(9), sigmEqui(17)
    blas_int :: b_incx, b_incy, b_n
    character(len=4), parameter :: epsiName(6) = (/'EPXX', 'EPYY', 'EPZZ', &
                                                   'EPXY', 'EPXZ', 'EPYZ'/)
    character(len=4), parameter :: sigmName(6) = (/'SIXX', 'SIYY', 'SIZZ', &
                                                   'SIXY', 'SIXZ', 'SIYZ'/)
    character(len=4), parameter :: gradName(9) = (/'F11', 'F12', 'F13', &
                                                   'F21', 'F22', 'F23', &
                                                   'F31', 'F32', 'F33'/)
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
    if (lLoadGrad) then
        nbCmpEpsi = 9
    else
        nbCmpEpsi = 6
    end if

! - Transform increment of strains for event
    epsiIncrTabl(1:nbCmpEpsi) = epsiIncr(1:nbCmpEpsi)
    if (.not. lLoadGrad) then
        epsiIncrTabl(4:6) = (1.d0/rac2)*epsiIncrTabl(4:6)
    end if

! - Compute increment of stress
    sigmIncr = sigmCurr-sigmPrev

! - Compute invariants of increment of stress
    sigmIncr(1:3) = (1.d0/rac2)*sigmIncr(1:3)
    call fgequi(sigmIncr, 'SIGM_DIR', 3, sigmEqui)

! - Create table for incremental values
    call detrsd('TABLE', tablIncr)
    call tbcrsd(tablIncr, 'V')
    call tbajpa(tablIncr, tablNbPara, tablParaName, tablParaType)
!
    if (tablType .eq. 0) then
! ----- Save increment of strains in table
        tablVale(2:2+nbCmpEpsi-1) = epsiIncrTabl(1:nbCmpEpsi)

! ----- Save increment of stresses in table
        tablVale(2+nbCmpEpsi:2+nbCmpEpsi+6-1) = sigmIncr(1:6)

! ----- Save invariants of increment of stress in table
        tablVale(nbCmpEpsi+8) = sigmEqui(16)
        tablVale(nbCmpEpsi+9) = sigmEqui(1)

! ----- Compute increment of internal state variables and save them in table
        tablVale(2+nbCmpEpsi+8:2+nbCmpEpsi+8+nbVariTabl-1) = vip(1:nbVariTabl)-vim(1:nbVariTabl)

! ----- Add line in table
        call tbajli(tablIncr, tablNbPara, tablParaName, [0], tablVale, &
                    [c16Dummy], k8Dummy, 0)
!
    else
        call wkvect(variIncrJv, 'V V R8', nbVariTabl, jvVariIncr)

! ----- Compute increment of internal state variables
        b_n = to_blas_int(nbVariTabl)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, vip, b_incx, zr(jvVariIncr), b_incy)
        b_n = to_blas_int(nbVariTabl)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, -1.d0, vim, b_incx, zr(jvVariIncr), b_incy)
!
        tablVale(1) = 0.d0
        vk8(1) = 'EPSI'
        do i = 1, nbCmpEpsi
            tablVale(2) = epsiIncrTabl(i)
            if (lLoadGrad) then
                vk8(2) = gradName(i)
            else
                vk8(2) = epsiName(i)
            end if
            call tbajli(tablIncr, tablNbPara, tablParaName, [0], tablVale, &
                        [c16Dummy], vk8, 0)
!
        end do
        vk8(1) = 'SIGM'
        do i = 1, 6
            tablVale(2) = sigmIncr(i)
            vk8(2) = sigmName(i)
            call tbajli(tablIncr, tablNbPara, tablParaName, [0], tablVale, &
                        [c16Dummy], vk8, 0)
!
        end do
        vk8(1) = 'SIEQ'
        tablVale(2) = sigmEqui(1)
        vk8(2) = 'VMIS'
        call tbajli(tablIncr, tablNbPara, tablParaName, [0], tablVale, &
                    [c16Dummy], vk8, 0)
!
        tablVale(2) = sigmEqui(16)
        vk8(2) = 'TRACE'
        call tbajli(tablIncr, tablNbPara, tablParaName, [0], tablVale, &
                    [c16Dummy], vk8, 0)
!
        vk8(1) = 'VARI'
        do i = 1, nbVariTabl
            tablVale(2) = zr(jvVariIncr-1+i)
            vk8(2) = variName(i)
            call tbajli(tablIncr, tablNbPara, tablParaName, [0], tablVale, &
                        [c16Dummy], vk8, 0)
        end do
!
    end if

! - Check event-driven
    call pmevdr(sddisc, tablType, tablIncr, &
                liccvg, lIterNewtMaxi, conver, &
                newtLoopAction)
!
    call jedetr(variIncrJv)
!
    call jedema()
end subroutine
