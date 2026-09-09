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
subroutine ddi_kit_read(factorKeyword, iFactorKeyword, &
                        rela_flua, rela_plas, rela_cpla, rela_coup)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/getvtx.h"
#include "asterfort/assert.h"
#include "asterfort/utmess.h"
!
    character(len=16), intent(in) :: factorKeyword
    integer(kind=8), intent(in) :: iFactorKeyword
    character(len=16), intent(out) :: rela_flua, rela_plas, rela_cpla, rela_coup
!
! --------------------------------------------------------------------------------------------------
!
! KIT_DDI
!
! Read informations for KIT
!
! --------------------------------------------------------------------------------------------------
!
! In  factorKeyword    : factor keyword to read (COMPORTEMENT)
! In  iFactorKeyword   : index of factor keyword
! Out rela_flua        : comportment relation for fluage
! Out rela_plas        : comportment relation for plasticity
! Out rela_cpla        : comportment relation for plane stress (GLRC)
! Out rela_coup        : comportment relation for coupling (GLRC)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter ::  dmflua = 5, dmplas = 9
    integer(kind=8) :: ikit, ii, nbFactorKeyword
    character(len=16) :: rela_kit(2)
    character(len=16), parameter :: poflua(dmflua) = (/'BETON_GRANGER  ', 'BETON_GRANGER_V', &
                                                       'BETON_UMLV     ', 'GLRC_DM        ', &
                                                       'GLRC_DAMAGE    '/)
    character(len=16), parameter :: poplas(dmplas) = (/'ELAS           ', 'VMIS_ISOT_TRAC ', &
                                                       'VMIS_ISOT_PUIS ', 'VMIS_ISOT_LINE ', &
                                                       'VMIS_CINE_LINE ', 'ROUSS_PR       ', &
                                                       'BETON_DOUBLE_DP', 'ENDO_ISOT_BETON', &
                                                       'MAZARS         '/)
!
! --------------------------------------------------------------------------------------------------
!
    rela_flua = 'VIDE'
    rela_plas = 'VIDE'
    rela_cpla = 'VIDE'
    rela_coup = 'VIDE'

! - Read command file
    call getvtx(factorKeyword, 'RELATION_KIT', iocc=iFactorKeyword, nbval=0, &
                nbret=nbFactorKeyword)
    nbFactorKeyword = -nbFactorKeyword
    ASSERT(nbFactorKeyword .le. 2)
    call getvtx(factorKeyword, 'RELATION_KIT', iocc=iFactorKeyword, nbval=nbFactorKeyword, &
                vect=rela_kit)

! - Get relations for kit
    do ikit = 1, nbFactorKeyword

! ----- Fluage
        do ii = 1, dmflua
            if (rela_kit(ikit) .eq. poflua(ii)) then
                rela_flua = rela_kit(ikit)
                cycle
            end if
        end do

! ----- Elasto-plastic
        do ii = 1, dmplas
            if (rela_kit(ikit) .eq. poplas(ii)) then
                rela_plas = rela_kit(ikit)
                cycle
            end if
        end do
    end do

! - Compatibility
    if (rela_flua(1:13) .eq. 'BETON_GRANGER') then
        if (rela_plas .ne. 'ELAS' .and. rela_plas .ne. 'VMIS_ISOT_TRAC' .and. rela_plas &
            .ne. 'VMIS_ISOT_PUIS' .and. rela_plas .ne. 'VMIS_ISOT_LINE' .and. rela_plas &
            .ne. 'ROUSS_PR' .and. rela_plas .ne. 'BETON_DOUBLE_DP') then
            call utmess('F', 'COMPOR3_2', sk=rela_plas)
        end if
    else if (rela_flua .eq. 'BETON_UMLV') then
        if (rela_plas .ne. 'ENDO_ISOT_BETON' .and. rela_plas .ne. 'MAZARS') then
            call utmess('F', 'COMPOR3_3', sk=rela_plas)
        end if
    else if (rela_flua(1:4) .eq. 'GLRC') then
        if (rela_plas .ne. 'VMIS_ISOT_TRAC' .and. rela_plas .ne. 'VMIS_ISOT_LINE' .and. &
            rela_plas .ne. 'VMIS_CINE_LINE') then
            call utmess('F', 'COMPOR3_4', sk=rela_plas)
        end if
    else
        call utmess('F', 'COMPOR3_6', sk=rela_flua)
    end if

! - For GLRC: internal Deborst Algorithm and special internal variables
    if (rela_flua(1:4) .eq. 'GLRC') then
        rela_coup = 'DDI_PLAS_ENDO'
        rela_cpla = 'DEBORST'
    end if
!
end subroutine
