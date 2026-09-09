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
subroutine nmdocn(ds_conv, mesh)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infdbg.h"
#include "asterfort/utmess.h"
#include "asterfort/SetResi.h"
#include "asterfort/SetResiRefe.h"
!
    type(NL_DS_Conv), intent(inout) :: ds_conv
    character(len=8), optional:: mesh
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Convergence management
!
! Read parameters for convergence management
!
! --------------------------------------------------------------------------------------------------
!
! IO  ds_conv          : datastructure for convergence management
! IN  mesh             : mesh name
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    character(len=16), parameter :: factorKeyword = 'CONVERGENCE'
    integer(kind=8) :: iret, iret_rela, iret_maxi, iret_refe, iret_comp, para_inte, iret_iter
    real(kind=8) :: para_real
    character(len=24) :: answer
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE12_8')
    end if

! - Initializations
    iret_refe = 0
    iret_comp = 0

! - Get convergence parameters (maximum iterations)
    call getvis(factorKeyword, 'ITER_GLOB_MAXI', iocc=1, scal=para_inte)
    ds_conv%iter_glob_maxi = para_inte
    call getvis(factorKeyword, 'ITER_GLOB_ELAS', iocc=1, scal=para_inte, nbret=iret_iter)
    if (iret_iter .ne. 0) then
        call getvis(factorKeyword, 'ITER_GLOB_ELAS', iocc=1, scal=para_inte)
        ds_conv%iter_glob_elas = para_inte
    end if

! - Get convergence parameters (residuals)
    call getvr8(factorKeyword, 'RESI_GLOB_RELA', iocc=1, scal=para_real, nbret=iret_rela)
    if (iret_rela .eq. 1) then
        call SetResi(ds_conv, type_='RESI_GLOB_RELA', &
                     user_para_=para_real, l_resi_test_=.true._1)
    end if
    call getvr8(factorKeyword, 'RESI_GLOB_MAXI', iocc=1, scal=para_real, nbret=iret_maxi)
    if (iret_maxi .eq. 1) then
        call SetResi(ds_conv, type_='RESI_GLOB_MAXI', &
                     user_para_=para_real, l_resi_test_=.true._1)
    end if
    call getvr8(factorKeyword, 'RESI_COMP_RELA', iocc=1, scal=para_real, nbret=iret_comp)
    if (iret_comp .eq. 1) then
        call SetResi(ds_conv, type_='RESI_COMP_RELA', &
                     user_para_=para_real, l_resi_test_=.true._1)
    end if
    call getvr8(factorKeyword, 'RESI_REFE_RELA', iocc=1, scal=para_real, nbret=iret_refe)
    if (iret_refe .eq. 1) then
        call SetResi(ds_conv, type_='RESI_REFE_RELA', &
                     user_para_=para_real, l_resi_test_=.true._1)
    end if
    call getvr8(factorKeyword, 'RESI_GLOB_RELA', iocc=1, scal=para_real, nbret=iret_rela)
    if (iret_rela .eq. 1) then
        call SetResi(ds_conv, type_='RESI_GLOB_RELA', &
                     user_para_=para_real, l_resi_test_=.true._1)
    end if
    call getvtx(factorKeyword, 'VERIF', iocc=1, scal=answer, nbret=iret_rela)
    ds_conv%lCritereOr = ASTER_FALSE
    if (iret_rela .eq. 1) then
        if (answer == "TOUT") then
            ds_conv%lCritereOr = ASTER_FALSE
        elseif (answer == "AU_MOINS_UN") then
            ds_conv%lCritereOr = ASTER_TRUE
        else
            write (6, *) "answer:", answer
            ASSERT(ASTER_FALSE)
        end if
    end if

! - Reference residuals
    if (iret_refe .eq. 1) then
        ASSERT(present(mesh))
        call SetResiRefe(mesh, ds_conv)
    end if

! - Forced convergence
    call getvtx(factorKeyword, 'ARRET', iocc=1, scal=answer, nbret=iret)
    if (iret .gt. 0) then
        ds_conv%l_stop = answer .eq. 'OUI'
    end if
!
end subroutine
