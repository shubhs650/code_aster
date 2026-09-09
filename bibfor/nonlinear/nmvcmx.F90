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
subroutine nmvcmx(materField, varcRefe, varc)
!
    implicit none
!
#include "jeveux.h"
#include "asterc/r8maem.h"
#include "asterfort/celces.h"
#include "asterfort/cesexi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nmvcex.h"
#include "asterfort/utmess.h"
!
    character(len=24), intent(in) :: materField, varcRefe
    character(len=19), intent(in) :: varc
!
! --------------------------------------------------------------------------------------------------
!
! Material - External state variables (VARC)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nbcmp, nbcmp2
    character(len=8) :: valk(2)
    character(len=19), parameter :: chscom = '&&NMVCMX.COMVAL_SIM'
    character(len=19), parameter :: chsref = '&&NMVCMX.COMREF_SIM'
    character(len=24) :: varcAllCurr, varcAllRefe
    integer(kind=8) :: jcesd, jcesl, nbCell, nbpt, nbsp, icmp
    integer(kind=8) :: jcrsd, jcrsl, iCell, ipt, isp, iad, iad2
    integer(kind=8) :: cellNumeMaxi, cellNumeMini, iref
    real(kind=8) :: valeMini, valeMaxi, valr(2)
    real(kind=8) :: valeur, valref
    character(len=8), pointer :: cvrcvarc(:) => null()
    real(kind=8), pointer :: cesv(:) => null()
    real(kind=8), pointer :: crsv(:) => null()
    character(len=8), pointer :: cvrcnom(:) => null()
!
! ----------------------------------------------------------------------
!
    call jemarq()

! - Get fields for external state variables

    call nmvcex('TOUT', varcRefe, varcAllRefe)
    call nmvcex('TOUT', varc, varcAllCurr)

! - TRANSFO. EN CHAM_NO_S
    call celces(varcAllCurr, 'V', chscom)
    call celces(varcAllRefe, 'V', chsref)
!
    call utmess('A+', 'MECANONLINE2_97')
!
!     CALCUL DU MIN / MAX
!
!     DESCRIPTEUR
    call jeveuo(chscom//'.CESD', 'L', jcesd)
    call jeveuo(chsref//'.CESD', 'L', jcrsd)
!     PRESENCE DES CMP (R)
    call jeveuo(chscom//'.CESL', 'L', jcesl)
    call jeveuo(chsref//'.CESL', 'L', jcrsl)
!     VALEUR DES CMP (R)
    call jeveuo(chscom//'.CESV', 'L', vr=cesv)
    call jeveuo(chsref//'.CESV', 'L', vr=crsv)
!
!     RECUPERATION DES NOMS DES VARC
    call jelira(materField(1:8)//'.CVRCNOM', 'LONMAX', ival=nbcmp2)
    call jeveuo(materField(1:8)//'.CVRCNOM', 'L', vk8=cvrcnom)
    call jeveuo(materField(1:8)//'.CVRCVARC', 'L', vk8=cvrcvarc)
!
    nbCell = zi(jcesd-1+1)
!
    do icmp = 1, nbcmp2
        valeMaxi = -r8maem()
        valeMini = r8maem()
        cellNumeMini = 0
        cellNumeMaxi = 0
        iref = 0
        if (cvrcvarc(icmp) .eq. 'TEMP' .or. cvrcvarc(icmp) .eq. 'SECH') then
            iref = 1
        end if
!
        do iCell = 1, nbCell
            nbcmp = zi(jcesd-1+5+4*(iCell-1)+3)
            if (nbcmp .eq. 0) cycle
            call cesexi('C', jcrsd, jcrsl, iCell, 1, &
                        1, icmp, iad2)
            if (iad2 .le. 0) cycle
!
!           VALEURS DE REFERENCE
            if (iref .eq. 1) then
                call cesexi('C', jcrsd, jcrsl, iCell, 1, &
                            1, icmp, iad2)
                valref = crsv(iad2)
            end if
            nbpt = zi(jcesd-1+5+4*(iCell-1)+1)
            nbsp = zi(jcesd-1+5+4*(iCell-1)+2)
            do ipt = 1, nbpt
                do isp = 1, nbsp
                    call cesexi('C', jcesd, jcesl, iCell, ipt, &
                                isp, icmp, iad)
                    if (iad .gt. 0) then
                        valeur = cesv(iad)
                        if (isnan(valeur)) cycle
!
                        if (iref .eq. 1) then
                            valeur = abs(valeur-valref)
                        end if
                        if (valeur .gt. valeMaxi) then
                            cellNumeMaxi = iCell
                            valeMaxi = valeur
                        end if
                        if (valeur .lt. valeMini) then
                            cellNumeMini = iCell
                            valeMini = valeur
                        end if
                    end if
                end do
            end do
        end do
        if (cellNumeMaxi .gt. 0) then
            valk(1) = cvrcvarc(icmp)
            valk(2) = cvrcnom(icmp)
            valr(1) = valeMaxi
            valr(2) = valeMini
            if (iref .eq. 1) then
                call utmess('A+', 'MECANONLINE2_95', &
                            nk=2, valk=valk, &
                            nr=2, valr=valr, &
                            ni=2, vali=[cellNumeMaxi, cellNumeMini])
            else
                call utmess('A+', 'MECANONLINE2_94', &
                            nk=2, valk=valk, &
                            nr=2, valr=valr, &
                            ni=2, vali=[cellNumeMaxi, cellNumeMini])
            end if
        end if
    end do
    call utmess('A', 'MECANONLINE2_93')
!
    call jedetr(chscom)
    call jedetr(chsref)
    call jedema()
end subroutine
