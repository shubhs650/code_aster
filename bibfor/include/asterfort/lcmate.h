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
#include "asterfort/Behaviour_type.h"
!
interface
    subroutine lcmate(materPara, &
                      carcri, relaComp, typmod1, &
                      nmat, tempd, tempf, tref, rungeKutta, &
                      typma, hsr, materd, materf, matcst, &
                      nbcomm, cpmono, pgl, itmax, &
                      toler, ndt, ndi, nr, &
                      nvi, vind, nfs, nsg, toutms, &
                      nhsr, numhsr, sigd, multComp_)
        use MaterialPara_type
        type(Material_Para), intent(in) :: materPara
        real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
        character(len=16), intent(in) :: relaComp
        character(len=8), intent(in) :: typmod1
        integer(kind=8), intent(in) :: nvi
        character(len=16), optional, intent(in) :: multComp_
        integer(kind=8) :: nmat
        real(kind=8) :: tempd
        real(kind=8) :: tempf
        real(kind=8) :: tref
        integer(kind=8) :: rungeKutta
        character(len=8) :: typma
        real(kind=8) :: hsr(*)
        real(kind=8) :: materd(nmat, 2)
        real(kind=8) :: materf(nmat, 2)
        character(len=3) :: matcst
        integer(kind=8) :: nbcomm(*)
        character(len=24) :: cpmono(*)
        real(kind=8) :: pgl(3, 3)
        integer(kind=8) :: itmax
        real(kind=8) :: toler
        integer(kind=8) :: ndt
        integer(kind=8) :: ndi
        integer(kind=8) :: nr
        real(kind=8) :: vind(*)
        integer(kind=8) :: nfs
        integer(kind=8) :: nsg
        real(kind=8) :: toutms(*)
        integer(kind=8) :: nhsr
        integer(kind=8) :: numhsr(*)
        real(kind=8) :: sigd(6)
    end subroutine lcmate
end interface
