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
subroutine lcjohm(materParaFPG1, &
                  lSigm, lMatr, lVari, &
                  kpg, npg, &
                  addeme, advico, ndim, dimdef, &
                  dimcon, nbvari, defgem, defgep, varim, &
                  varip, sigm, sigp, drde, ouvh, &
                  retcom)
!
    use MaterialPara_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
!
    type(Material_Para), intent(in) :: materParaFPG1
    integer(kind=8) :: kpg, npg, addeme, advico, ndim, dimdef, dimcon, nbvari
    real(kind=8) :: defgem(dimdef), varim(nbvari), sigm(dimcon)
    aster_logical, intent(in) :: lSigm, lMatr, lVari
    integer(kind=8) :: retcom
    real(kind=8) :: defgep(dimdef), varip(nbvari), sigp(dimcon)
    real(kind=8) :: drde(dimdef, dimdef), ouvh
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: poum = "+"
    integer(kind=8) :: i
    real(kind=8) :: kni, umc, gamma, kt, clo, valr(2), tmecn, tmecs
    integer(kind=8), parameter :: nbProp = 4
    integer(kind=8) :: propCode(nbProp)
    real(kind=8) :: propVale(nbProp)
    character(len=8), parameter :: propName(nbProp) = (/'K    ', 'DMAX ', 'GAMMA', 'KT   '/)
!
! --------------------------------------------------------------------------------------------------
!
    call rcvalb(materParaFPG1%schemePara%fami, &
                materParaFPG1%schemePara%kpg, &
                materParaFPG1%schemePara%ksp, &
                poum, &
                materParaFPG1%jvMaterCode, &
                ' ', 'JOINT_BANDIS', &
                0, ' ', [0.d0], &
                nbProp, propName, propVale, &
                propCode, 1)
    kni = propVale(1)
    umc = propVale(2)
    gamma = propVale(3)
    kt = propVale(4)

! - MISE A JOUR FERMETURE
    clo = 0.d0
    ouvh = varim(advico)
    clo = umc-ouvh
    clo = clo-defgep(addeme)+defgem(addeme)

! - Internal state variable
    if (lVari) then
        ouvh = umc-clo
        varip(advico) = ouvh
    end if
!
! - Stresses
!
    if (lSigm) then
        if ((clo .gt. umc) .or. (clo .lt. -1.d-3)) then
            valr(1) = clo
            valr(2) = umc
            call utmess('A', 'ALGORITH17_11', nr=2, valr=valr)
            retcom = 1
            goto 999
        end if
        do i = 1, dimcon
            sigp(i) = 0.d0
        end do
        sigp(1) = sigm(1)-kni/(1-clo/umc)**gamma*(varim(advico)-varip(advico))
        do i = 2, ndim
            sigp(i) = sigm(i)+kt*(defgep(addeme+1)-defgem(addeme+1))
        end do
    end if

! - CALCUP OPERATEUR TANGENT
    if (lMatr .and. (kpg .le. npg)) then
        tmecn = kni/(1-clo/umc)**gamma
        tmecs = kt
        drde(addeme, addeme) = tmecn
        do i = 2, ndim
            drde(addeme+i-1, addeme+i-1) = tmecs
        end do
    end if
!
999 continue
!
end subroutine
