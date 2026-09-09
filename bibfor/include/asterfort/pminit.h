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
    subroutine pminit(tablName, tablNbPara, tablType, &
                      tablParaName, tablParaType, tablVale, &
                      pgl, lRota, &
                      epsiPrev, sigmPrev, &
                      nbVari, vim, vip, &
                      loadEpsiType, loadType, loadFunc, coefImpo, &
                      coefAdim, typeMatrPred, lMatrElas, matrElas, lPrintMatr, option, &
                      variName, nbVariTabl, &
                      sddisc, ds_conv, ds_algopara, sderro, materPara)
        use NonLin_Datastructure_type
        use MaterialPara_type
        character(len=8), intent(out) :: tablName
        integer(kind=8), intent(out) :: tablNbPara, tablType
        character(len=16), allocatable, intent(out) :: tablParaName(:)
        character(len=8), allocatable, intent(out) :: tablParaType(:)
        real(kind=8), allocatable, intent(out) :: tablVale(:)
        real(kind=8), intent(out) :: pgl(3, 3)
        aster_logical, intent(out) :: lRota
        real(kind=8), intent(out) :: epsiPrev(9), sigmPrev(6)
        integer(kind=8), intent(in) :: nbVari
        real(kind=8), intent(out) :: vim(nbVari), vip(nbVari)
        integer(kind=8), intent(out) :: loadEpsiType, loadType(9)
        character(len=8), intent(out) :: loadFunc(9)
        real(kind=8), intent(out) :: coefImpo(6, 12), coefAdim
        integer(kind=8), intent(out) :: typeMatrPred
        aster_logical, intent(out) :: lMatrElas
        real(kind=8), intent(out) :: matrElas(6, 6)
        aster_logical, intent(out) :: lPrintMatr
        character(len=16), intent(out) :: option
        character(len=8), intent(out) :: variName(nbVari)
        integer(kind=8), intent(out) :: nbVariTabl
        character(len=19), intent(out) :: sddisc
        type(NL_DS_Conv), intent(inout) :: ds_conv
        type(NL_DS_AlgoPara), intent(inout) :: ds_algopara
        character(len=24), intent(out) :: sderro
        type(Material_Para), intent(inout) :: materPara
    end subroutine pminit
end interface
