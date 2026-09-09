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
interface
    subroutine rk21co(materPara, &
                      relaComp, typmod1, &
                      matcst, nbcomm, cpmono, nfs, &
                      nsg, toutms, nvi, nmat, y, &
                      kp, ee, a, h, pgl, &
                      nbphas, cothe, coeff, dcothe, dcoeff, &
                      coel, x, pas, neps, epsd, &
                      detot, nhsr, numhsr, hsr, itmax, &
                      toler, iret)
        use MaterialPara_type
        type(Material_Para), intent(in) :: materPara
        character(len=16), intent(in) :: relaComp
        character(len=8), intent(in) :: typmod1
        integer(kind=8) :: nhsr
        integer(kind=8) :: nmat
        integer(kind=8) :: nvi
        integer(kind=8) :: nsg
        character(len=3) :: matcst
        integer(kind=8) :: nbcomm(nmat, 3)
        character(len=24) :: cpmono(5*nmat+1)
        integer(kind=8) :: nfs
        real(kind=8) :: toutms(*)
        real(kind=8) :: y(nvi)
        integer(kind=8) :: kp
        real(kind=8) :: ee(nvi)
        real(kind=8) :: a(nvi)
        real(kind=8) :: h
        real(kind=8) :: pgl(3, 3)
        integer(kind=8) :: nbphas
        real(kind=8) :: cothe(nmat)
        real(kind=8) :: coeff(nmat)
        real(kind=8) :: dcothe(nmat)
        real(kind=8) :: dcoeff(nmat)
        real(kind=8) :: coel(nmat)
        real(kind=8) :: x
        real(kind=8) :: pas
        integer(kind=8) :: neps
        real(kind=8) :: epsd(6)
        real(kind=8) :: detot(6)
        integer(kind=8) :: numhsr(*)
        real(kind=8) :: hsr(nsg, nsg, nhsr)
        integer(kind=8) :: itmax
        real(kind=8) :: toler
        integer(kind=8) :: iret
    end subroutine rk21co
end interface
