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
! aslint: disable=W0413
!
subroutine te0154(option, nomte)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/ElasticityMaterial_type.h"
#include "asterfort/get_elas_para.h"
#include "asterfort/getDensity.h"
#include "asterfort/jevech.h"
#include "asterfort/lonele.h"
#include "asterfort/matrot.h"
#include "asterfort/pmavec.h"
#include "asterfort/ptenci.h"
#include "asterfort/ptenpo.h"
#include "asterfort/ptenth.h"
#include "asterfort/tecach.h"
#include "asterfort/utpvgl.h"
#include "asterfort/verift.h"
#include "jeveux.h"
!
    character(len=*), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: MECA_BARRE, MECA_2D_BARRE
!
! Options: ECIN_ELEM, EPOT_ELEM, EPSI_ELGA, SIEF_ELGA
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = 'RIGI'
    character(len=3), parameter :: stopz = 'ONO'
    integer(kind=8), parameter :: nc = 3, nno = 2, kpg = 1, ksp = 1
    integer(kind=8), parameter :: elasID = ELAS_ISOT
    character(len=16), parameter :: elasKeyword = "ELAS"
    real(kind=8) :: pgl(3, 3), klc(6, 6)
    real(kind=8) :: ugr(6), ulr(12), flr(6)
    real(kind=8) :: enerTher, e, rho
    aster_logical :: hasTemp, lRod3d
    real(kind=8) :: aire, epsiTher, xfl1, xfl4, xl, xmas, xrig
    integer(kind=8) :: ii, iif, itype, kanl, iret
    integer(kind=8) :: jvSigm, jvEner, jvPuls, jvEpsi
    integer(kind=8) :: jvMater, jvCaorien, jvSect
    integer(kind=8) :: jvDisp, jvVite
!
! --------------------------------------------------------------------------------------------------
!
    lRod3d = nomte .eq. 'MECA_BARRE'

! - Get material parameters
    epsiTher = 0.d0
    hasTemp = ASTER_FALSE
    e = 0.d0
    if (option .ne. "EPSI_ELGA") then
        call jevech('PMATERC', 'L', jvMater)

! ----- Compute thermal strains
        call verift(fami, 1, 1, '+', zi(jvMater), epsth_=epsiTher)
        if (epsiTher .ne. 0.d0) then
            hasTemp = ASTER_TRUE
        end if

! ----- Get young modulus
        call get_elas_para(fami, zi(jvMater), '+', kpg, ksp, &
                           elasID, elasKeyword, e_=e)

    end if

! - Length of rod
    if (nomte .eq. 'MECA_BARRE') then
        xl = lonele()
    else if (nomte .eq. 'MECA_2D_BARRE') then
        xl = lonele(dime=2)
    end if

! - Get area
    aire = 0.d0
    if (option .ne. 'EPSI_ELGA') then
        call jevech('PCAGNBA', 'L', jvSect)
        aire = zr(jvSect)
    end if

! - Get orientation
    call jevech('PCAORIE', 'L', jvCaorien)
    call matrot(zr(jvCaorien), pgl)

! - Get displacements or speeds
    ugr = 0.d0
    if (option .ne. 'ECIN_ELEM') then
        call jevech('PDEPLAR', 'L', jvDisp)
        if (lRod3d) then
            do ii = 1, 6
                ugr(ii) = zr(jvDisp+ii-1)
            end do
        else
            ugr(1) = zr(jvDisp+1-1)
            ugr(2) = zr(jvDisp+2-1)
            ugr(4) = zr(jvDisp+3-1)
            ugr(5) = zr(jvDisp+4-1)
        end if
    else
        call tecach(stopz, 'PVITESR', 'L', iret, iad=jvVite)
        if (iret .eq. 0) then
            if (lRod3d) then
                do ii = 1, 6
                    ugr(ii) = zr(jvVite+ii-1)
                end do
            else
                ugr(1) = zr(jvVite+1-1)
                ugr(2) = zr(jvVite+2-1)
                ugr(4) = zr(jvVite+3-1)
                ugr(5) = zr(jvVite+4-1)
            end if
        else
            call tecach(stopz, 'PDEPLAR', 'L', iret, iad=jvDisp)
            if (iret .eq. 0) then
                if (lRod3d) then
                    do ii = 1, 6
                        ugr(ii) = zr(jvDisp+ii-1)
                    end do
                else
                    ugr(1) = zr(jvDisp+1-1)
                    ugr(2) = zr(jvDisp+2-1)
                    ugr(4) = zr(jvDisp+3-1)
                    ugr(5) = zr(jvDisp+4-1)
                end if
            else
                ASSERT(ASTER_FALSE)
            end if
        end if
    end if

! - Change coordinate system from global to local
    call utpvgl(nno, nc, pgl, ugr, ulr)
!
    klc = 0.d0
    if (option .eq. 'EPOT_ELEM') then
        call jevech('PENERDR', 'E', jvEner)
        xrig = e*aire/xl
        klc(1, 1) = xrig
        klc(1, 4) = -xrig
        klc(4, 1) = -xrig
        klc(4, 4) = xrig
        iif = 0
        call ptenpo(6, ulr, klc, zr(jvEner), iif, iif)
        if (hasTemp) then
            call ptenth(ulr, xl, epsiTher, 6, klc, enerTher)
            zr(jvEner) = zr(jvEner)-enerTher
        end if

    else if (option .eq. 'ECIN_ELEM') then
        call getDensity(zi(jvMater), rho, elasKeywordZ_=elasKeyword)
        call jevech('PENERCR', 'E', jvEner)
        call jevech('POMEGA2', 'L', jvPuls)
        xmas = rho*aire*xl/6.d0
        klc(1, 1) = xmas*2.d0
        klc(2, 2) = xmas*2.d0
        klc(3, 3) = xmas*2.d0
        klc(4, 4) = xmas*2.d0
        klc(5, 5) = xmas*2.d0
        klc(6, 6) = xmas*2.d0
        klc(1, 4) = xmas
        klc(4, 1) = xmas
        klc(2, 5) = xmas
        klc(5, 2) = xmas
        klc(3, 6) = xmas
        klc(6, 3) = xmas
        iif = 0
        itype = 50
        kanl = 1
        call ptenci(6, ulr, klc, zr(jvPuls), zr(jvEner), itype, kanl, iif)
!
    else if (option .eq. 'EPSI_ELGA') then
        call jevech('PDEFOPG', 'E', jvEpsi)
        zr(jvEpsi-1+1) = (ulr(4)-ulr(1))/xl

    else
        xrig = e*aire/xl
        klc(1, 1) = xrig
        klc(1, 4) = -xrig
        klc(4, 1) = -xrig
        klc(4, 4) = xrig
        call pmavec('ZERO', 6, klc, ulr, flr)
        if (hasTemp) then
            xfl1 = -epsiTher*e*aire
            xfl4 = -xfl1
            flr(1) = flr(1)-xfl1
            flr(4) = flr(4)-xfl4
        end if
        if (option .eq. 'SIEF_ELGA') then
            call jevech('PCONTRR', 'E', jvSigm)
            zr(jvSigm) = -flr(1)
        end if
    end if
!
end subroutine
