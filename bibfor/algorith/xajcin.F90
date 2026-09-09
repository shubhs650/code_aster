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

subroutine xajcin(model, option, nbFieldInMax, lchin, lpain, &
                  nbFieldIn)
!
    implicit none
!
#include "asterfort/assert.h"
!
    integer(kind=8), intent(in) :: nbFieldInMax
    character(len=*), intent(in) :: model
    character(len=*), intent(in) :: option
    character(len=*), intent(inout) :: lpain(nbFieldInMax)
    character(len=*), intent(inout) :: lchin(nbFieldInMax)
    integer(kind=8), intent(inout) :: nbFieldIn
!
! --------------------------------------------------------------------------------------------------
!
! Add XFEM fields for input fields
!
!  -> OPTIONS : - CHAR_MECA_TEMP_R
!               - CHAR_THER_PARO_F
!               - CHAR_THER_PARO_R
!               - FULL_MECA
!               - RIGI_MECA_*
!               - RAPH_MECA
!               - CHAR_MECA_NEUM
!               - MASS_THER
!
! --------------------------------------------------------------------------------------------------
!
! In  model  : name of model
! In  option : option to select input fields
! In  nbFieldInMax : maximum number of input fields
! IO  lpain  : list of parameters
! IO  lchin  : list of fields
! IO  nbin   : number of input fields
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nbadd
!
! --------------------------------------------------------------------------------------------------
!
    if ((option .eq. 'CHAR_MECA_TEMP_R') .or. &
        (option(1:9) .eq. 'FULL_MECA') .or. &
        (option(1:9) .eq. 'RAPH_MECA') .or. &
        (option(1:9) .eq. 'RIGI_MECA') .or. &
        (option(1:9) .eq. 'RIGI_GEOM')) then
!
        nbadd = 12
        ASSERT(nbFieldIn+nbadd .le. nbFieldInMax)
        lpain(nbFieldIn+1) = 'PPINTTO'
        lchin(nbFieldIn+1) = model(1:8)//'.TOPOSE.PIN'
        lpain(nbFieldIn+2) = 'PCNSETO'
        lchin(nbFieldIn+2) = model(1:8)//'.TOPOSE.CNS'
        lpain(nbFieldIn+3) = 'PHEAVTO'
        lchin(nbFieldIn+3) = model(1:8)//'.TOPOSE.HEA'
        lpain(nbFieldIn+4) = 'PLONCHA'
        lchin(nbFieldIn+4) = model(1:8)//'.TOPOSE.LON'
        lpain(nbFieldIn+5) = 'PBASLOR'
        lchin(nbFieldIn+5) = model(1:8)//'.BASLOC'
        lpain(nbFieldIn+6) = 'PLSN'
        lchin(nbFieldIn+6) = model(1:8)//'.LNNO'
        lpain(nbFieldIn+7) = 'PLST'
        lchin(nbFieldIn+7) = model(1:8)//'.LTNO'
        lpain(nbFieldIn+8) = 'PSTANO'
        lchin(nbFieldIn+8) = model(1:8)//'.STNO'
        lpain(nbFieldIn+9) = 'PPMILTO'
        lchin(nbFieldIn+9) = model(1:8)//'.TOPOSE.PMI'
        lpain(nbFieldIn+10) = 'PFISNO'
        lchin(nbFieldIn+10) = model(1:8)//'.FISSNO'
        lpain(nbFieldIn+11) = 'PHEA_NO'
        lchin(nbFieldIn+11) = model(1:8)//'.TOPONO.HNO'
        lpain(nbFieldIn+12) = 'PHEA_FA'
        lchin(nbFieldIn+12) = model(1:8)//'.TOPONO.HFA'
        nbFieldIn = nbFieldIn+nbadd
!
    elseif (option .eq. 'CHAR_MECA_NEUM') then
!
        nbadd = 19
        ASSERT(nbFieldIn+nbadd .le. nbFieldInMax)
        lpain(nbFieldIn+1) = 'PPINTTO'
        lchin(nbFieldIn+1) = model(1:8)//'.TOPOSE.PIN'
        lpain(nbFieldIn+2) = 'PCNSETO'
        lchin(nbFieldIn+2) = model(1:8)//'.TOPOSE.CNS'
        lpain(nbFieldIn+3) = 'PHEAVTO'
        lchin(nbFieldIn+3) = model(1:8)//'.TOPOSE.HEA'
        lpain(nbFieldIn+4) = 'PLONCHA'
        lchin(nbFieldIn+4) = model(1:8)//'.TOPOSE.LON'
        lpain(nbFieldIn+5) = 'PLSN'
        lchin(nbFieldIn+5) = model(1:8)//'.LNNO'
        lpain(nbFieldIn+6) = 'PLST'
        lchin(nbFieldIn+6) = model(1:8)//'.LTNO'
        lpain(nbFieldIn+7) = 'PSTANO'
        lchin(nbFieldIn+7) = model(1:8)//'.STNO'
        lpain(nbFieldIn+8) = 'PPMILTO'
        lchin(nbFieldIn+8) = model(1:8)//'.TOPOSE.PMI'
        lpain(nbFieldIn+9) = 'PFISNO'
        lchin(nbFieldIn+9) = model(1:8)//'.FISSNO'
        lpain(nbFieldIn+10) = 'PPINTER'
        lchin(nbFieldIn+10) = model(1:8)//'.TOPOFAC.OE'
        lpain(nbFieldIn+11) = 'PAINTER'
        lchin(nbFieldIn+11) = model(1:8)//'.TOPOFAC.AI'
        lpain(nbFieldIn+12) = 'PCFACE'
        lchin(nbFieldIn+12) = model(1:8)//'.TOPOFAC.CF'
        lpain(nbFieldIn+13) = 'PLONGCO'
        lchin(nbFieldIn+13) = model(1:8)//'.TOPOFAC.LO'
        lpain(nbFieldIn+14) = 'PBASECO'
        lchin(nbFieldIn+14) = model(1:8)//'.TOPOFAC.BA'
        lpain(nbFieldIn+15) = 'PHEA_NO'
        lchin(nbFieldIn+15) = model(1:8)//'.TOPONO.HNO'
        lpain(nbFieldIn+16) = 'PHEA_SE'
        lchin(nbFieldIn+16) = model(1:8)//'.TOPONO.HSE'
        lpain(nbFieldIn+17) = 'PHEA_FA'
        lchin(nbFieldIn+17) = model(1:8)//'.TOPONO.HFA'
        lpain(nbFieldIn+18) = 'PBASLOR'
        lchin(nbFieldIn+18) = model(1:8)//'.BASLOC'
        lpain(nbFieldIn+19) = 'PHEAVNO'
        lchin(nbFieldIn+19) = model(1:8)//'.HEAVNO'
        nbFieldIn = nbFieldIn+nbadd
!
    elseif (option .eq. 'REFE_FORC_NODA') then
!
        nbadd = 14
        ASSERT(nbFieldIn+nbadd .le. nbFieldInMax)
        lpain(nbFieldIn+1) = 'PPINTTO'
        lchin(nbFieldIn+1) = model(1:8)//'.TOPOSE.PIN'
        lpain(nbFieldIn+2) = 'PCNSETO'
        lchin(nbFieldIn+2) = model(1:8)//'.TOPOSE.CNS'
        lpain(nbFieldIn+3) = 'PHEAVTO'
        lchin(nbFieldIn+3) = model(1:8)//'.TOPOSE.HEA'
        lpain(nbFieldIn+4) = 'PLONCHA'
        lchin(nbFieldIn+4) = model(1:8)//'.TOPOSE.LON'
        lpain(nbFieldIn+5) = 'PLSN'
        lchin(nbFieldIn+5) = model(1:8)//'.LNNO'
        lpain(nbFieldIn+6) = 'PLST'
        lchin(nbFieldIn+6) = model(1:8)//'.LTNO'
        lpain(nbFieldIn+7) = 'PPMILTO'
        lchin(nbFieldIn+7) = model(1:8)//'.TOPOSE.PMI'
        lpain(nbFieldIn+8) = 'PPINTER'
        lchin(nbFieldIn+8) = model(1:8)//'.TOPOFAC.OE'
        lpain(nbFieldIn+9) = 'PAINTER'
        lchin(nbFieldIn+9) = model(1:8)//'.TOPOFAC.AI'
        lpain(nbFieldIn+10) = 'PCFACE'
        lchin(nbFieldIn+10) = model(1:8)//'.TOPOFAC.CF'
        lpain(nbFieldIn+11) = 'PBASECO'
        lchin(nbFieldIn+11) = model(1:8)//'.TOPOFAC.BA'
        lpain(nbFieldIn+12) = 'PHEA_NO'
        lchin(nbFieldIn+12) = model(1:8)//'.TOPONO.HNO'
        lpain(nbFieldIn+13) = 'PLONFA'
        lchin(nbFieldIn+13) = model(1:8)//'.TOPOFAC.LO'
        lpain(nbFieldIn+14) = 'PBASLOR'
        lchin(nbFieldIn+14) = model(1:8)//'.BASLOC'
        nbFieldIn = nbFieldIn+nbadd
!
    elseif (option .eq. 'FORC_NODA') then
!
        nbadd = 11
        ASSERT(nbFieldIn+nbadd .le. nbFieldInMax)
        lpain(nbFieldIn+1) = 'PPINTTO'
        lchin(nbFieldIn+1) = model(1:8)//'.TOPOSE.PIN'
        lpain(nbFieldIn+2) = 'PCNSETO'
        lchin(nbFieldIn+2) = model(1:8)//'.TOPOSE.CNS'
        lpain(nbFieldIn+3) = 'PHEAVTO'
        lchin(nbFieldIn+3) = model(1:8)//'.TOPOSE.HEA'
        lpain(nbFieldIn+4) = 'PLONCHA'
        lchin(nbFieldIn+4) = model(1:8)//'.TOPOSE.LON'
        lpain(nbFieldIn+5) = 'PPMILTO'
        lchin(nbFieldIn+5) = model(1:8)//'.TOPOSE.PMI'
        lpain(nbFieldIn+6) = 'PBASLOR'
        lchin(nbFieldIn+6) = model(1:8)//'.BASLOC'
        lpain(nbFieldIn+7) = 'PLSN'
        lchin(nbFieldIn+7) = model(1:8)//'.LNNO'
        lpain(nbFieldIn+8) = 'PLST'
        lchin(nbFieldIn+8) = model(1:8)//'.LTNO'
        lpain(nbFieldIn+9) = 'PSTANO'
        lchin(nbFieldIn+9) = model(1:8)//'.STNO'
        lpain(nbFieldIn+10) = 'PFISNO'
        lchin(nbFieldIn+10) = model(1:8)//'.FISSNO'
        lpain(nbFieldIn+11) = 'PHEA_NO'
        lchin(nbFieldIn+11) = model(1:8)//'.TOPONO.HNO'
        nbFieldIn = nbFieldIn+nbadd
!
    elseif (option(6:14) .eq. 'THER_PARO') then
!
        nbadd = 9
        ASSERT(nbFieldIn+nbadd .le. nbFieldInMax)
        lpain(nbFieldIn+1) = 'PPINTER'
        lchin(nbFieldIn+1) = model(1:8)//'.TOPOFAC.OE'
        lpain(nbFieldIn+2) = 'PAINTER'
        lchin(nbFieldIn+2) = model(1:8)//'.TOPOFAC.AI'
        lpain(nbFieldIn+3) = 'PCFACE'
        lchin(nbFieldIn+3) = model(1:8)//'.TOPOFAC.CF'
        lpain(nbFieldIn+4) = 'PLONGCO'
        lchin(nbFieldIn+4) = model(1:8)//'.TOPOFAC.LO'
        lpain(nbFieldIn+5) = 'PLST'
        lchin(nbFieldIn+5) = model(1:8)//'.LTNO'
        lpain(nbFieldIn+6) = 'PSTANO'
        lchin(nbFieldIn+6) = model(1:8)//'.STNO'
        lpain(nbFieldIn+7) = 'PBASECO'
        lchin(nbFieldIn+7) = model(1:8)//'.TOPOFAC.BA'
        lpain(nbFieldIn+8) = 'PLSN'
        lchin(nbFieldIn+8) = model(1:8)//'.LNNO'
        lpain(nbFieldIn+9) = 'PHEA_NO'
        lchin(nbFieldIn+9) = model(1:8)//'.TOPONO.HNO'
        nbFieldIn = nbFieldIn+nbadd
!
    elseif (option .eq. 'RIGI_THER') then
!
        nbadd = 9
        ASSERT(nbFieldIn+nbadd .le. nbFieldInMax)
        lpain(nbFieldIn+1) = 'PSTANO'
        lchin(nbFieldIn+1) = model(1:8)//'.STNO'
        lpain(nbFieldIn+2) = 'PPINTTO'
        lchin(nbFieldIn+2) = model(1:8)//'.TOPOSE.PIN'
        lpain(nbFieldIn+3) = 'PCNSETO'
        lchin(nbFieldIn+3) = model(1:8)//'.TOPOSE.CNS'
        lpain(nbFieldIn+4) = 'PHEAVTO'
        lchin(nbFieldIn+4) = model(1:8)//'.TOPOSE.HEA'
        lpain(nbFieldIn+5) = 'PLONCHA'
        lchin(nbFieldIn+5) = model(1:8)//'.TOPOSE.LON'
        lpain(nbFieldIn+6) = 'PBASLOR'
        lchin(nbFieldIn+6) = model(1:8)//'.BASLOC'
        lpain(nbFieldIn+7) = 'PLSN'
        lchin(nbFieldIn+7) = model(1:8)//'.LNNO'
        lpain(nbFieldIn+8) = 'PLST'
        lchin(nbFieldIn+8) = model(1:8)//'.LTNO'
        lpain(nbFieldIn+9) = 'PHEA_NO'
        lchin(nbFieldIn+9) = model(1:8)//'.TOPONO.HNO'
        nbFieldIn = nbFieldIn+nbadd
!
    elseif (option .eq. 'CHAR_THER_EVOL') then
!
        nbadd = 9
        ASSERT(nbFieldIn+nbadd .le. nbFieldInMax)
        lpain(nbFieldIn+1) = 'PSTANO'
        lchin(nbFieldIn+1) = model(1:8)//'.STNO'
        lpain(nbFieldIn+2) = 'PPINTTO'
        lchin(nbFieldIn+2) = model(1:8)//'.TOPOSE.PIN'
        lpain(nbFieldIn+3) = 'PCNSETO'
        lchin(nbFieldIn+3) = model(1:8)//'.TOPOSE.CNS'
        lpain(nbFieldIn+4) = 'PHEAVTO'
        lchin(nbFieldIn+4) = model(1:8)//'.TOPOSE.HEA'
        lpain(nbFieldIn+5) = 'PLONCHA'
        lchin(nbFieldIn+5) = model(1:8)//'.TOPOSE.LON'
        lpain(nbFieldIn+6) = 'PBASLOR'
        lchin(nbFieldIn+6) = model(1:8)//'.BASLOC'
        lpain(nbFieldIn+7) = 'PLSN'
        lchin(nbFieldIn+7) = model(1:8)//'.LNNO'
        lpain(nbFieldIn+8) = 'PLST'
        lchin(nbFieldIn+8) = model(1:8)//'.LTNO'
        lpain(nbFieldIn+9) = 'PHEA_NO'
        lchin(nbFieldIn+9) = model(1:8)//'.TOPONO.HNO'
        nbFieldIn = nbFieldIn+nbadd
    elseif (option .eq. 'MASS_THER') then
!
        nbadd = 9
        ASSERT(nbFieldIn+nbadd .le. nbFieldInMax)
        lpain(nbFieldIn+1) = 'PSTANO'
        lchin(nbFieldIn+1) = model(1:8)//'.STNO'
        lpain(nbFieldIn+2) = 'PPINTTO'
        lchin(nbFieldIn+2) = model(1:8)//'.TOPOSE.PIN'
        lpain(nbFieldIn+3) = 'PCNSETO'
        lchin(nbFieldIn+3) = model(1:8)//'.TOPOSE.CNS'
        lpain(nbFieldIn+4) = 'PHEAVTO'
        lchin(nbFieldIn+4) = model(1:8)//'.TOPOSE.HEA'
        lpain(nbFieldIn+5) = 'PLONCHA'
        lchin(nbFieldIn+5) = model(1:8)//'.TOPOSE.LON'
        lpain(nbFieldIn+6) = 'PBASLOR'
        lchin(nbFieldIn+6) = model(1:8)//'.BASLOC'
        lpain(nbFieldIn+7) = 'PLSN'
        lchin(nbFieldIn+7) = model(1:8)//'.LNNO'
        lpain(nbFieldIn+8) = 'PLST'
        lchin(nbFieldIn+8) = model(1:8)//'.LTNO'
        lpain(nbFieldIn+9) = 'PHEA_NO'
        lchin(nbFieldIn+9) = model(1:8)//'.TOPONO.HNO'
        nbFieldIn = nbFieldIn+nbadd
    elseif (option(1:5) .eq. 'MASS_') then
!
        nbadd = 10
        ASSERT(nbFieldIn+nbadd .le. nbFieldInMax)
        lpain(nbFieldIn+1) = 'PPINTTO'
        lchin(nbFieldIn+1) = model(1:8)//'.TOPOSE.PIN'
        lpain(nbFieldIn+2) = 'PHEAVTO'
        lchin(nbFieldIn+2) = model(1:8)//'.TOPOSE.HEA'
        lpain(nbFieldIn+3) = 'PLONCHA'
        lchin(nbFieldIn+3) = model(1:8)//'.TOPOSE.LON'
        lpain(nbFieldIn+4) = 'PCNSETO'
        lchin(nbFieldIn+4) = model(1:8)//'.TOPOSE.CNS'
        lpain(nbFieldIn+5) = 'PBASLOR'
        lchin(nbFieldIn+5) = model(1:8)//'.BASLOC'
        lpain(nbFieldIn+6) = 'PLSN'
        lchin(nbFieldIn+6) = model(1:8)//'.LNNO'
        lpain(nbFieldIn+7) = 'PLST'
        lchin(nbFieldIn+7) = model(1:8)//'.LTNO'
        lpain(nbFieldIn+8) = 'PSTANO'
        lchin(nbFieldIn+8) = model(1:8)//'.STNO'
        lpain(nbFieldIn+9) = 'PHEA_NO'
        lchin(nbFieldIn+9) = model(1:8)//'.TOPONO.HNO'
        lpain(nbFieldIn+10) = 'PPMILTO'
        lchin(nbFieldIn+10) = model(1:8)//'.TOPOSE.PMI'
        nbFieldIn = nbFieldIn+nbadd

    elseif (option .eq. 'ENEL_ELEM') then
        nbadd = 4
        ASSERT(nbFieldIn+nbadd .le. nbFieldInMax)

        lpain(nbFieldIn+1) = 'PPINTTO'
        lchin(nbFieldIn+1) = model(1:8)//'.TOPOSE.PIN'
        lpain(nbFieldIn+2) = 'PPMILTO'
        lchin(nbFieldIn+2) = model(1:8)//'.TOPOSE.PMI'
        lpain(nbFieldIn+3) = 'PCNSETO'
        lchin(nbFieldIn+3) = model(1:8)//'.TOPOSE.CNS'
        lpain(nbFieldIn+4) = 'PLONCHA'
        lchin(nbFieldIn+4) = model(1:8)//'.TOPOSE.LON'
        nbFieldIn = nbFieldIn+nbadd

    elseif (option .eq. 'COOR_ELGA') then
        nbadd = 4
        ASSERT(nbFieldIn+nbadd .le. nbFieldInMax)

        lpain(nbFieldIn+1) = 'PPINTTO'
        lchin(nbFieldIn+1) = model(1:8)//'.TOPOSE.PIN'
        lpain(nbFieldIn+2) = 'PPMILTO'
        lchin(nbFieldIn+2) = model(1:8)//'.TOPOSE.PMI'
        lpain(nbFieldIn+3) = 'PCNSETO'
        lchin(nbFieldIn+3) = model(1:8)//'.TOPOSE.CNS'
        lpain(nbFieldIn+4) = 'PLONCHA'
        lchin(nbFieldIn+4) = model(1:8)//'.TOPOSE.LON'
        nbFieldIn = nbFieldIn+nbadd

!
    else
        ASSERT(.false.)
    end if
!
end subroutine
