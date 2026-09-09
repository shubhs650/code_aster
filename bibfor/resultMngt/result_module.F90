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
! ==================================================================================================
!
! Module for management of result datastructures
!
! ==================================================================================================
!
module result_module
! ==================================================================================================
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: rsCopyPara
! ==================================================================================================
    private
#include "asterf_types.h"
#include "asterfort/copisd.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnopa.h"
#include "jeveux.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! rsCopyPara
!
! Copy parameters from one result to another one
!
! In  resultIn          : name of datastructure for input results
! In  resultOut         : name of datastructure for output results
! In  nbStore           : number of storing indexes
! Ptr listStore         : pointer to list of storing indexes
!
! --------------------------------------------------------------------------------------------------
    subroutine rsCopyPara(resultInZ, resultOutZ, nbStore, listStore, &
                          copyField_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: resultInZ, resultOutZ
        integer(kind=8), intent(in) :: nbStore
        integer(kind=8), pointer :: listStore(:)
        aster_logical, optional, intent(in) :: copyField_
! ----- Local
        character(len=24), parameter :: paraJvName = '&&CCBCOP.NOMS_PARA'
        integer(kind=8) :: nbParaAccess, nbPara, nbParaTotal, iParaTotal
        integer(kind=8) :: numeStore, iStore, iret
        character(len=19) :: resultIn, resultOut
        character(len=8) :: paraType
        integer(kind=8) :: jvParaIn, jvParaOut
        aster_logical :: copyField
        character(len=19) :: comporToCopy, comporToSave
        character(len=24) :: paraIn, paraOut
        character(len=16) :: paraName
        character(len=16), pointer :: listParaName(:) => null()
!   ------------------------------------------------------------------------------------------------
!
        call jemarq()

! ----- Initializations
        resultIn = resultInZ
        resultOut = resultOutZ
        copyField = ASTER_FALSE
        if (present(copyField_)) then
            copyField = copyField_
        end if

! ----- Acces to parameters
        call rsnopa(resultIn, 2, paraJvName, nbParaAccess, nbPara)
        nbParaTotal = nbParaAccess+nbPara
        call jeveuo(paraJvName, 'L', vk16=listParaName)

! ----- Copy parameters
        do iStore = 1, nbStore
            numeStore = listStore(iStore)
            do iParaTotal = 1, nbParaTotal
                paraName = listParaName(iParaTotal)
                call rsadpa(resultIn, 'L', 1, paraName, numeStore, &
                            1, sjv=jvParaIn, styp=paraType, istop=0)
                call rsadpa(resultOut, 'E', 1, paraName, numeStore, &
                            1, sjv=jvParaOut, styp=paraType)
                if (copyField) then
                    call rsexch(' ', resultIn, 'COMPORTEMENT', numeStore, comporToCopy, iret)
                    if (iret .eq. 0) then
                        call rsexch(' ', resultOut, 'COMPORTEMENT', numeStore, comporToSave, iret)
                        call copisd('CHAMP_GD', 'G', comporToCopy, comporToSave)
                    end if
                end if
                if (paraType(1:1) .eq. 'I') then
                    zi(jvParaOut) = zi(jvParaIn)
                else if (paraType(1:1) .eq. 'R') then
                    zr(jvParaOut) = zr(jvParaIn)
                else if (paraType(1:1) .eq. 'C') then
                    zc(jvParaOut) = zc(jvParaIn)
                else if (paraType(1:3) .eq. 'K80') then
                    zk80(jvParaOut) = zk80(jvParaIn)
                else if (paraType(1:3) .eq. 'K32') then
                    zk32(jvParaOut) = zk32(jvParaIn)
                else if (paraType(1:3) .eq. 'K24') then
                    zk24(jvParaOut) = zk24(jvParaIn)
                    if (copyField) then
                        paraIn = zk24(jvParaIn)
                        if (paraName .eq. 'EXCIT' .and. paraIn(1:2) .ne. '  ') then
                            paraOut = resultOut(1:8)//paraIn(9:)
                            call copisd('LISTE_CHARGES', 'G', paraIn(1:19), paraOut(1:19))
                            zk24(jvParaOut) = paraOut
                        end if
                    end if
                else if (paraType(1:3) .eq. 'K16') then
                    zk16(jvParaOut) = zk16(jvParaIn)
                else if (paraType(1:2) .eq. 'K8') then
                    zk8(jvParaOut) = zk8(jvParaIn)
                end if
            end do
        end do
        call jedetr(paraJvName)
!
        call jedema()
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module result_module
