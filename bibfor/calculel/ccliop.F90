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
subroutine ccliop(option, jvBaseName, listOptEffJv, nbOptEff)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
!
    character(len=8), intent(in) :: jvBaseName
    character(len=16), intent(in) :: option
    character(len=24), intent(out) :: listOptEffJv
    integer(kind=8), intent(out) :: nbOptEff
!
! --------------------------------------------------------------------------------------------------
!
!  CALC_CHAMP
!
!  Prepare list of options to compute (dependencies)
!
! --------------------------------------------------------------------------------------------------
!
! In  option            : option to compute
! In  jvBaseName        : string base for names of JEVEUX objects
! Out nbOptEff          : number of options to compute (after dependencies)
! Out listOptEffJv      : list of options to compute (after dependencies)
!
!       NOLIOP(N)   = OPTION (EN ARGUEMENT DE LA ROUTINE)
!       NOLIOP(N-1) = OPTIO2 (DONT DEPEND OPTION)
!       ...
!       NOLIOP(M)   = OPTIOM (DONT DEPEND OPTION)
!       ...
!       NOLIOP(P)   = OPTIOP (DONT DEPEND OPTIO2)
!       ...
!
!  D'OU NOLORI DE TAILLE 2*N
!       NOLORI(N)   = N-1
!       NOLORI(N-1) = M
!       NOLORI(3) = P
!       ...
!  LA LISTE NOLDEP RAPPELLE LA DEPENDANCE TEMPORELLE
!  EX : SI OPTION DEPEND DE L'INSTANT N+1 D'OPTIO2 ALORS :
!       NOLDEP(N-1) = 'NP1'
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbOptMax = 100
    integer(kind=8) :: lopor1(nbOptMax), lopor2(nbOptMax)
    character(len=1) :: isodep(nbOptMax)
    character(len=4) :: lopdep(nbOptMax)
    character(len=16) :: loptio(nbOptMax)
    integer(kind=8) :: cataNbIn, iCataIn
    integer(kind=8) :: optionNume, iopdeb, iOpt, cataInOptionNume
    integer(kind=8) :: nopous
    integer(kind=8) :: jliori, jlidep, jlnoin, jlisde
    aster_logical :: optionHasBeenAdded
    character(len=16) :: optionCurr, cataInStep, cataInOption
    character(len=24) :: cataInName
    character(len=24) :: nolori, noldep, noliin, nolisd
    integer(kind=8), pointer :: cataDescopt(:) => null()
    character(len=24), pointer :: cataLocalis(:) => null()
    character(len=16), pointer :: listOptEff(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    nbOptEff = 0

! - Access to catalog of options
    if (option(6:9) .ne. 'NOEU') then
        call jenonu(jexnom('&CATA.OP.NOMOPT', option), optionNume)
        ASSERT(optionNume .ne. 0)
        call jeveuo(jexnum('&CATA.OP.DESCOPT', optionNume), 'L', vi=cataDescopt)
        cataNbIn = cataDescopt(2)
        if (cataNbIn .eq. 0) goto 999
        call jeveuo(jexnum('&CATA.OP.LOCALIS', optionNume), 'L', vk24=cataLocalis)
    end if

! - INITIALISATION DES ENTIERS
    iopdeb = 1
    nbOptEff = 1
    nopous = nbOptEff

! - Initialization of list of options (dependencies: the first one)
    loptio(1) = option
    lopdep(1) = 'NSP'
!
40  continue

! - Create temporary list of dependencies
    do iOpt = iopdeb, nopous
        isodep(iOpt) = ' '
        optionCurr = loptio(iOpt)
        optionHasBeenAdded = ASTER_FALSE

        if (optionCurr(6:9) .eq. 'NOEU') then
! --------- Dependency from ELNO field
            cataInOption = option(1:5)//'ELNO'
            call jenonu(jexnom('&CATA.OP.NOMOPT', cataInOption), optionNume)
            if (optionNume .ne. 0) then
                nbOptEff = nbOptEff+1
                loptio(nbOptEff) = cataInOption
                lopor1(iOpt) = nbOptEff
                lopor2(iOpt) = 0
                lopor1(nbOptEff) = 0
                lopor2(nbOptEff) = 0
                lopdep(nbOptEff) = 'N'
                optionHasBeenAdded = ASTER_TRUE
            end if
            cataNbIn = 0
        else if ((optionCurr .eq. 'DEPL') .or. (optionCurr .eq. 'VITE') .or. &
                 (optionCurr .eq. 'ACCE') .or. (optionCurr .eq. 'TEMP') .or. &
                 (optionCurr .eq. 'VARI_ELGA')) then
! --------- No dependencies
            cataNbIn = 0
        else
            call jenonu(jexnom('&CATA.OP.NOMOPT', optionCurr), optionNume)
            call jeveuo(jexnum('&CATA.OP.DESCOPT', optionNume), 'L', vi=cataDescopt)
            cataNbIn = cataDescopt(2)
            if (cataNbIn .eq. 0) cycle
            call jeveuo(jexnum('&CATA.OP.LOCALIS', optionNume), 'L', vk24=cataLocalis)
        end if

! ----- Loop on input parameters for this option
        do iCataIn = 1, cataNbIn
! --------- Get properties of input parameter
            cataInName = cataLocalis(3*iCataIn-1)
            cataInStep = cataLocalis(3*iCataIn) (1:16)

! --------- Is cataInName is an option ?
            cataInOption = cataInName(1:16)
            call jenonu(jexnom('&CATA.OP.NOMOPT', cataInOption), cataInOptionNume)

! --------- Depend on previous or next option, not itself
            if (optionCurr .eq. cataInName(1:16)) then
                if (cataInStep .eq. 'NP1') then
                    isodep(iOpt) = '+'
                else if (cataInStep .eq. 'NM1') then
                    isodep(iOpt) = '-'
                end if
                cycle
            end if
!
            if (cataInOptionNume .ne. 0) then
! ------------- Dependency is another option
                nbOptEff = nbOptEff+1
                loptio(nbOptEff) = cataInOption
                if (.not. optionHasBeenAdded) then
                    lopor1(iOpt) = nbOptEff
                    lopor2(iOpt) = 0
                end if
                lopor1(nbOptEff) = 0
                lopor2(nbOptEff) = 0
                lopdep(nbOptEff) = cataInStep(1:4)
                optionHasBeenAdded = ASTER_TRUE
                isodep(nbOptEff) = ' '
            elseif (((cataInOption .eq. 'DEPL') .or. &
                     (cataInOption .eq. 'SIEF_ELGA') .or. &
                     (cataInOption .eq. 'VARI_ELGA')) .and. &
                    ((cataInStep .eq. 'NP1') .or. (cataInStep(1:3) .eq. 'NM1'))) then
! ------------- Dependency is field at time step NP1 or NM1
                nbOptEff = nbOptEff+1
                loptio(nbOptEff) = cataInOption
                if (.not. optionHasBeenAdded) then
                    lopor1(iOpt) = nbOptEff
                    lopor2(iOpt) = 0
                end if
                lopor1(nbOptEff) = 0
                lopor2(nbOptEff) = 0
                lopdep(nbOptEff) = cataInStep(1:4)
                optionHasBeenAdded = ASTER_TRUE
                isodep(nbOptEff) = ' '
            end if
        end do
        if (.not. optionHasBeenAdded) then
            lopor1(iOpt) = 0
            lopor2(iOpt) = 0
        else
            lopor2(iOpt) = nbOptEff
        end if
    end do

!     SI ON A AJOUTE UNE OPTION LORS DE LA DERNIERE PASSE, ON
!     DOIT CHERCHER SES DEPENDANCES
    if (optionHasBeenAdded) then
        iopdeb = nopous+1
        nopous = nbOptEff
        goto 40
    end if
    ASSERT(nbOptEff .le. nbOptMax)

! - Create objects for final list of dependencies
    listOptEffJv = jvBaseName//'.LISOPT'
    nolori = jvBaseName//'.LISORI'
    noldep = jvBaseName//'.LISDEP'
    noliin = jvBaseName//'.LNOINS'
    nolisd = jvBaseName//'.ISODEP'
    call wkvect(listOptEffJv, 'V V K16', nbOptEff, vk16=listOptEff)
    call wkvect(nolori, 'V V I', 2*nbOptEff, jliori)
    call wkvect(noldep, 'V V K8', nbOptEff, jlidep)
    call wkvect(noliin, 'V V K24', nbOptEff, jlnoin)
    call wkvect(nolisd, 'V V K8', nbOptEff, jlisde)

! - Create final list of dependencies
    do iOpt = 1, nbOptEff
!       ON PARCOURT LA LISTE A L'ENVERS PUISQUE PAR CONSTRUCTION
!       LES OPTIONS 'D'EN HAUT' DEPENDENT DES OPTIONS 'D'EN BAS'
        cataInOption = loptio(nbOptEff-iOpt+1)

! ----- Add this option to compute
        listOptEff(iOpt) = cataInOption

!       ON REGARDE DE QUELLES OPTIONS DEPEND L'OPTION COURANTE
        if (lopor1(nbOptEff-iOpt+1) .ne. 0) then
            zi(jliori+2*iOpt-2) = nbOptEff-lopor2(nbOptEff-iOpt+1)+1
            zi(jliori+2*iOpt-1) = nbOptEff-lopor1(nbOptEff-iOpt+1)+1
        else
            zi(jliori+2*iOpt-2) = 0
            zi(jliori+2*iOpt-1) = 0
        end if
        zk8(jlidep+iOpt-1) = lopdep(nbOptEff-iOpt+1)
        zk8(jlisde+iOpt-1) = isodep(nbOptEff-iOpt+1)
    end do
!
999 continue
!
    call jedema()
!
end subroutine
