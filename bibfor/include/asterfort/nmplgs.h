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
    subroutine nmplgs(BEHInteg, &
                      ndim, nno1, nno2, npg, &
                      vff1, idfde1, &
                      vff2, idfde2, &
                      iw, geom, &
                      typmod, option, compor, carcri, &
                      instam, instap, &
                      ddlm, ddld, &
                      lgpg, sigm, vim, &
                      sigp, vip, &
                      matr, vect, codret, &
                      livois, &
                      nbvois, numa, lisoco, nbsoco, &
                      lVari, lSigm, lMatr, lVect)
        use Behaviour_type
        type(Behaviour_Integ), intent(inout) :: BEHInteg
        character(len=8), intent(in) :: typmod(2)
        character(len=16), intent(in) :: option, compor(COMPOR_SIZE)
        real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
        aster_logical, intent(in) :: lVari, lSigm, lMatr, lVect
        integer(kind=8), intent(in) :: ndim, nno1, nno2, npg
        real(kind=8), intent(in) :: vff1(nno1, npg)
        integer(kind=8), intent(in) :: idfde1
        real(kind=8), intent(in) :: vff2(nno2, npg)
        integer(kind=8), intent(in) :: idfde2
        integer(kind=8), intent(in) :: iw
        real(kind=8) :: geom(ndim, nno1)
        real(kind=8) :: instam
        real(kind=8) :: instap
        real(kind=8) :: ddlm(*)
        real(kind=8) :: ddld(*)
        integer(kind=8), intent(in) :: lgpg
        real(kind=8) :: sigm(2*ndim, npg)
        real(kind=8) :: vim(lgpg, npg)
        real(kind=8) :: sigp(2*ndim, npg)
        real(kind=8) :: vip(lgpg, npg)
        real(kind=8) :: matr(*)
        real(kind=8) :: vect(*)
        integer(kind=8) :: codret
        integer(kind=8), parameter :: nvoima = 12, nscoma = 4
        integer(kind=8) :: livois(1:nvoima)
        integer(kind=8) :: nbvois
        integer(kind=8) :: numa
        integer(kind=8) :: lisoco(1:nvoima, 1:nscoma, 1:2)
        integer(kind=8) :: nbsoco(1:nvoima)
    end subroutine nmplgs
end interface
