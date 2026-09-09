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
subroutine crvarc()
!
    use coorSyst_module, only: setOrieFields
    implicit none
!
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/cesvar.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exlima.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mecact.h"
#include "asterfort/megeom.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rscrsd.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnoch.h"
#include "asterfort/rsorac.h"
#include "asterfort/setStructFields.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
! ------------------------------------------------------------------------------
!
!                           COMMANDE    :   CREA_RESU
!
! ------------------------------------------------------------------------------
!
!       Création d'une SD de type "EVOL_THER" contenant la température sur les
!       couches des coques multicouche à partir :
!           d'un CHAM_GD de fonctions [ INST, [EPAIS|X,Y,Z] ]
!           d'un EVOL_THER contenant TEMP/TEMP_INF/TEMP_SUP
!
! ------------------------------------------------------------------------------
!
    character(len=16), parameter :: factorKeyw = "PREP_VARC"
    integer(kind=8), parameter :: nbFieldOut = 1, nbFieldInMax = 100
    character(len=8) :: lpaout(nbFieldOut), lpain(nbFieldInMax)
    character(len=19) :: lchout(nbFieldOut), lchin(nbFieldInMax)
    integer(kind=8) :: nbFieldIn
    integer(kind=8) :: nbfac, jvPara, jvTime, iret, iexi
    integer(kind=8) :: n1, n2, iTime, nbTime, numeStore, nbStore
    real(kind=8) :: timeCurr, prec
    character(len=5) :: LesInst, LeCas
    character(len=8) :: resultOut, caraElem, thermalField
    character(len=8) :: model2, model
    character(len=16) :: resultType, cmdName, crit
    character(len=19) :: modelLigrel, ligrelCalc, listInst, thermalResult
    character(len=19), parameter ::  chinst = '&&CRVRC1.CHINST'
    character(len=24) :: chtemp, chgeom, tempField
    integer(kind=8) :: lesordres(1)
    integer(kind=8) :: ibib
    character(len=8) :: k8bid
    complex(kind=8) :: c16bid
    integer(kind=8), pointer :: storeList(:) => null()
    real(kind=8), pointer :: timeList(:) => null()
    character(len=24), pointer :: celk(:) => null()
!
! ------------------------------------------------------------------------------
    call jemarq()
!
    call getfac(factorKeyw, nbfac)
    if (nbfac .eq. 0) goto 20
    ASSERT(nbfac .eq. 1)

! - Initializations
    call getres(resultOut, resultType, cmdName)
    lpain = " "
    lpaout = " "
    lchin = " "
    lchout = " "

!   C'est CHAM_GD ou EVOL_THER
    LeCas = ' '
    call getvid(factorKeyw, 'CHAM_GD', iocc=1, nbval=0, nbret=n1)
    call getvid(factorKeyw, 'EVOL_THER', iocc=1, nbval=0, nbret=n2)
    if (n1 .lt. 0) then
        call getvid(factorKeyw, 'CHAM_GD', iocc=1, scal=thermalField, nbret=n1)
        LeCas = 'CHAMP'
    else if (n2 .lt. 0) then
        call getvid(factorKeyw, 'EVOL_THER', iocc=1, scal=thermalResult, nbret=n2)
        LeCas = 'THERM'
    else
        ASSERT(ASTER_FALSE)
    end if
!
    call getvid(factorKeyw, 'MODELE', iocc=1, scal=model, nbret=n1)
    call getvid(factorKeyw, 'CARA_ELEM', iocc=1, scal=caraElem, nbret=n1)

!   On vérifie que le cara_elem s'appuie bien sur le modèle
    call jeexin(caraElem//'.CANBSP    .CELK', iexi)
    if (iexi .eq. 0) then
        call utmess('F', 'CALCULEL4_14')
    end if
    call jeveuo(caraElem//'.CANBSP    .CELK', 'L', vk24=celk)
    model2 = celk(1) (1:8)
    if (model2 .ne. model) then
        call utmess('F', 'CALCULEL4_15')
    end if
!
!   Les instants
    LesInst = ' '
    call getvr8(factorKeyw, 'INST', iocc=1, nbval=0, nbret=n1)
    call getvr8(factorKeyw, 'LIST_INST', iocc=1, nbval=0, nbret=n2)
    if (n1 .lt. 0) then
        nbTime = -n1
        AS_ALLOCATE(vr=timeList, size=nbTime)
        call getvr8(factorKeyw, 'INST', iocc=1, nbval=nbTime, vect=timeList, nbret=n1)
        ASSERT(nbTime .eq. n1)
        LesInst = 'LINST'
    else if (n2 .lt. 0) then
        call getvid(factorKeyw, 'LIST_INST', iocc=1, scal=listInst, nbret=n2)
        ASSERT(n2 .eq. 1)
        call jeveuo(listInst(1:19)//'.VALE', 'L', jvTime)
        call jelira(listInst(1:19)//'.VALE', 'LONMAX', nbTime)
        AS_ALLOCATE(vr=timeList, size=nbTime)
        do iTime = 1, nbTime
            timeList(iTime) = zr(jvTime+iTime-1)
        end do
        LesInst = 'LINST'
    end if
!
!   Si EVOL_THER
    prec = 0.0; crit = ' '
    if (LeCas .eq. 'THERM') then
        ! Si INST ou LIST_INST : pour recherche de INST dans la SD
        if (LesInst .eq. 'LINST') then
            call getvr8(factorKeyw, 'PRECISION', iocc=1, scal=prec, nbret=n1)
            call getvtx(factorKeyw, 'CRITERE', iocc=1, scal=crit, nbret=n2)
        else
            LesInst = 'TINST'
        end if
    end if
!
!   Si EVOL_THER, on vérifie qu'il y a des NUME_ORDRE
    if (LeCas .eq. 'THERM') then
        call jelira(thermalResult//'.ORDR', 'LONUTI', nbStore)
        call jeveuo(thermalResult//'.ORDR', 'L', vi=storeList)
        ASSERT(nbStore .gt. 0)
        if (LesInst .eq. 'TINST') then
            nbTime = nbStore
        end if
    end if

! - Create output result (no resuse !)
    call jeexin(resultOut//'           .DESC', iret)
    if (iret .ne. 0) then
        call utmess('F', 'CALCULEL7_6')
    end if
    call rscrsd('G', resultOut, 'EVOL_THER', nbTime)

! - Add fields for structural elements
    nbFieldIn = 0
    call setStructFields(caraElem, nbFieldInMax, lchin, lpain, nbFieldIn)

! - Add fields for orientation
    call setOrieFields(nbFieldInMax, lpain, lchin, &
                       nbFieldIn, caraElem)

! - Set output fields
    lpaout(1) = 'PTEMPCR'

! - Select FED
    if (LeCas .eq. 'CHAMP') then
        call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)

    else if (LeCas .eq. 'THERM') then
        call exlima(factorKeyw, 1, 'G', model, ligrelCalc)

    end if
!
    if (LeCas .eq. 'CHAMP') then
! ----- Add input fields
        call megeom(model, chgeom)
        nbFieldIn = nbFieldin+1
        lpain(nbFieldIn) = 'PTEMPEF'
        lchin(nbFieldIn) = thermalField
        nbFieldIn = nbFieldin+1
        lpain(nbFieldIn) = 'PINST_R'
        lchin(nbFieldIn) = chinst
        nbFieldIn = nbFieldin+1
        lpain(nbFieldIn) = 'PGEOMER'
        lchin(nbFieldIn) = chgeom(1:19)

        do iTime = 1, nbTime
! --------- Get current time and create map
            timeCurr = timeList(iTime)
            call mecact('V', chinst, 'MODELE', modelLigrel, 'INST_R', &
                        ncmp=1, nomcmp='INST', sr=timeCurr)

! --------- Get temperature field in output resultOut
            call rsexch(' ', resultOut, 'TEMP', iTime, tempField, iret)

! --------- Prepare field for dynamic components (VARi_R)
            call cesvar(caraElem, ' ', modelLigrel, tempField)

! --------- Compute
            lchout(1) = tempField(1:19)
            call calcul('S', 'PREP_VRC', modelLigrel, &
                        nbFieldIn, lchin, lpain, &
                        nbFieldOut, lchout, lpaout, &
                        'G', 'OUI')
            call detrsd('CHAM_ELEM_S', tempField)

! --------- Indexation of field in resultOut and save parameters
            call rsnoch(resultOut, 'TEMP', iTime)
            call rsadpa(resultOut, 'E', 1, 'INST', iTime, 0, sjv=jvPara)
            zr(jvPara) = timeCurr
            call rsadpa(resultOut, 'E', 1, 'CARAELEM', iTime, 0, sjv=jvPara)
            zk8(jvPara) = caraElem
            call detrsd('CHAMP', chinst)

        end do

    else if (LeCas .eq. 'THERM') then
! ----- Add input fields
        nbFieldIn = nbFieldin+1
        lpain(nbFieldIn) = 'PTEMPER'

        do iTime = 1, nbTime
            if (LesInst .eq. 'TINST') then
! ------------- Get current time
                numeStore = storeList(iTime)
                call rsadpa(thermalResult, 'L', 1, 'INST', numeStore, 0, sjv=jvPara)
                timeCurr = zr(jvPara)

            else
! ------------- Get current time
                timeCurr = timeList(iTime)
                call rsorac(thermalResult, 'INST', ibib, timeCurr, k8bid, c16bid, &
                            prec, crit, lesordres, 1, n1)
                if (n1 .lt. 0) then
                    call utmess('F', 'ALGORITH12_83', si=-n1, sr=timeCurr)
                else if (n1 .eq. 0) then
                    call utmess('F', 'ALGORITH12_84', sr=timeCurr)
                end if
                numeStore = lesordres(1)
            end if

! --------- Get temperature field in thermal resultOut
            call rsexch('F', thermalResult, 'TEMP', numeStore, chtemp, iret)
            lchin(nbFieldIn) = chtemp(1:19)

! --------- Get temperature field in output resultOut
            call rsexch(' ', resultOut, 'TEMP', numeStore, tempField, iret)

! --------- Prepare field for dynamic components (VARi_R)
            call cesvar(caraElem, ' ', ligrelCalc, tempField)

! --------- Compute
            lchout(1) = tempField(1:19)
            call calcul('S', 'PREP_VRC', ligrelCalc, &
                        nbFieldIn, lchin, lpain, &
                        nbFieldOut, lchout, lpaout, &
                        'G', 'OUI')
            call detrsd('CHAM_ELEM_S', tempField)

! --------- Indexation of field in resultOut and save parameters
            call rsnoch(resultOut, 'TEMP', numeStore)
            call rsadpa(resultOut, 'E', 1, 'INST', numeStore, 0, sjv=jvPara)
            zr(jvPara) = timeCurr
            call rsadpa(resultOut, 'E', 1, 'CARAELEM', numeStore, 0, sjv=jvPara)
            zk8(jvPara) = caraElem
        end do
    end if
!
    AS_DEALLOCATE(vr=timeList)
!
20  continue
    call jedema()
end subroutine
