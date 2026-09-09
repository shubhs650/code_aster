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
#include "asterf_types.h"
!
interface
    subroutine xfnohm(ds_thm, &
                      nno, npg, &
                      ipoids, ivf, idfde, &
                      geom, congem, b, dfdi, dfdi2, &
                      r, vectu, mecani, press1, &
                      dimcon, nddls, nddlm, dimuel, nmec, &
                      np1, ndim, axi, dimenr, nnop, &
                      nnops, nnopm, igeom, jpintt, jpmilt, &
                      jheavn, lonch, cnset, heavt, enrmec, enrhyd, &
                      nfiss, nfh, jfisno)
        use THM_type
        type(THM_DS), intent(inout) :: ds_thm
        integer(kind=8) :: nnops
        integer(kind=8) :: nnop
        integer(kind=8) :: dimenr
        integer(kind=8) :: ndim
        integer(kind=8) :: dimuel
        aster_logical :: fnoevo
        real(kind=8) :: deltat
        integer(kind=8) :: nno
        integer(kind=8) :: npg
        integer(kind=8) :: ipoids
        integer(kind=8) :: ivf
        integer(kind=8) :: idfde
        real(kind=8) :: geom(ndim, nnop)
        real(kind=8) :: congem(*)
        real(kind=8) :: b(dimenr, dimuel)
        real(kind=8) :: dfdi(nnop, ndim)
        real(kind=8) :: dfdi2(nnops, ndim)
        real(kind=8) :: r(1:dimenr)
        real(kind=8) :: vectu(dimuel)
        integer(kind=8) :: mecani(5)
        integer(kind=8) :: press1(7)
        integer(kind=8) :: dimcon
        integer(kind=8) :: nddls
        integer(kind=8) :: nddlm
        integer(kind=8) :: nmec
        integer(kind=8) :: np1
        integer(kind=8) :: jheavn
        aster_logical :: axi
        integer(kind=8) :: nnopm
        integer(kind=8) :: igeom
        integer(kind=8) :: jpintt
        integer(kind=8) :: jpmilt
        integer(kind=8) :: lonch(10)
        integer(kind=8) :: cnset(*)
        integer(kind=8) :: heavt(*)
        integer(kind=8) :: enrmec(3)
        integer(kind=8) :: enrhyd(3)
        integer(kind=8) :: nfiss
        integer(kind=8) :: nfh
        integer(kind=8) :: jfisno
    end subroutine xfnohm
end interface
