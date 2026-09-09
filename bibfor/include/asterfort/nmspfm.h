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
    subroutine nmspfm(BEHinteg, typmod, ndim, nno, nddl, nddlsym, &
                      nno_p, nno_s, nddl_p, nddl_s, &
                      npg, lgpg, wref, &
                      vff_s, vf_p, pgl, geom, nsigm, mate, matint, matpou, &
                      option, deplm, ddepl, sigm, sigp, fint, &
                      ktan, vim, vip, carcri, compor, &
                      tm, tp, coopg, matsym, lMatr, lVect, lSigm, lElas, &
                      codret)
        use Behaviour_type
        type(Behaviour_Integ) :: BEHinteg
        character(len=8) :: typmod(2)
        integer(kind=8) :: ndim
        integer(kind=8) :: nno
        integer(kind=8) :: nddl
        integer(kind=8) :: nddlsym
        integer(kind=8) :: nno_p
        integer(kind=8) :: nno_s
        integer(kind=8) :: nddl_p
        integer(kind=8) :: nddl_s
        integer(kind=8) :: npg
        integer(kind=8) :: lgpg
        real(kind=8) :: wref(npg)
        real(kind=8) :: vf_p(ndim, npg)
        real(kind=8) :: vff_s(nno, npg)
        real(kind=8) :: pgl(3, 3)
        real(kind=8) :: geom(3*nno)
        integer(kind=8) :: nsigm
        integer(kind=8) :: mate
        character(len=8)  :: matint
        character(len=8)  :: matpou
        character(len=16) :: option
        real(kind=8) :: deplm(nddl)
        real(kind=8) :: ddepl(nddl)
        real(kind=8) :: sigm(nsigm, npg)
        real(kind=8) :: sigp(nsigm, npg)
        real(kind=8) :: fint(nddl)
        real(kind=8) :: ktan(nddlsym)
        real(kind=8) :: vim(lgpg, npg)
        real(kind=8) :: vip(lgpg, npg)
        real(kind=8) :: carcri(*)
        character(len=16) :: compor(COMPOR_SIZE)
        real(kind=8) :: tm
        real(kind=8) :: tp
        real(kind=8) :: coopg(4, npg)
        aster_logical, intent(in) :: matsym, lMatr, lVect, lSigm, lElas
        integer(kind=8) :: codret
    end subroutine nmspfm
end interface
