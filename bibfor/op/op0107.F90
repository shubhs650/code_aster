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

subroutine op0107()
    implicit none
!     OPERATEUR   POST_ELEM
!     ------------------------------------------------------------------
!
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterfort/assert.h"
#include "asterfort/chpve2.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/medomp.h"
#include "asterfort/peaire.h"
#include "asterfort/pecage.h"
#include "asterfort/pecapo.h"
#include "asterfort/pechli.h"
#include "asterfort/peecin.h"
#include "asterfort/peeint.h"
#include "asterfort/peepot.h"
#include "asterfort/peingl.h"
#include "asterfort/pemain.h"
#include "asterfort/pemima.h"
#include "asterfort/penorm.h"
#include "asterfort/peritr.h"
#include "asterfort/pevolu.h"
#include "asterfort/peweib.h"
#include "asterfort/pewext.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsutnu.h"
#include "asterfort/titre.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    integer(kind=8) :: numeHarm, iret, jordr, n1, n2, nbFactorKeyword, nbordr, nc, np, nr, ier
    real(kind=8) :: prec
    character(len=8) :: model, caraElem, deform, result, crit, mesh, k8b
    character(len=16) :: concep, nomcmd
    character(len=19) :: tablOut, knum, tabtyp(3)
    character(len=24) :: materField, materCode, chdef
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infmaj()
    k8b = " "

! - Get output table
    call getres(tablOut, concep, nomcmd)

! - Get input result
    call getvid(' ', 'RESULTAT', scal=result, nbret=nr)
    if (nr .eq. 0) then
        result = ' '
    end if

!
    call getfac('TRAV_EXT', nbFactorKeyword)
    if (nbFactorKeyword .ne. 0) then
        call pewext(tablOut)
    end if
!
    call getfac('CHAR_LIMITE', nbFactorKeyword)
    if (nbFactorKeyword .ne. 0) then
        call medomp(result, model, mateco=materCode)
        call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
        ASSERT(.not. isParallelMesh(mesh))
        call pechli(tablOut, model, materCode)
    end if
!
    call getfac('AIRE_INTERNE', nbFactorKeyword)
    if (nbFactorKeyword .ne. 0) then
        call medomp(result, model)
        call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
        ASSERT(.not. isParallelMesh(mesh))
        call peaire(tablOut, mesh, nbFactorKeyword)
    end if
!
    call getfac('MASS_INER', nbFactorKeyword)
    if (nbFactorKeyword .ne. 0) then
        call medomp(result, model, materField, materCode, caraElem, numeHarm)
        call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
        ASSERT(.not. isParallelMesh(mesh))
        chdef = ' '
        call getvtx(' ', 'GEOMETRIE', scal=deform, nbret=n1)
        if (deform .eq. 'DEFORMEE') then
            call getvid(' ', 'CHAM_GD', scal=chdef, nbret=n2)
            if (n2 .eq. 0) then
                tabtyp(1) = 'NOEU#DEPL_R'
                tabtyp(2) = 'NOEU#TEMP_R'
                tabtyp(3) = 'ELEM#ENER_R'
                knum = '&&OP0107.NUME_ORDRE'
                call getvid(' ', 'RESULTAT', scal=result, nbret=nr)
                call getvr8(' ', 'PRECISION', scal=prec, nbret=np)
                call getvtx(' ', 'CRITERE', scal=crit, nbret=nc)
                call rsutnu(result, ' ', 0, knum, nbordr, &
                            prec, crit, iret)
                if (nbordr .ne. 1) then
                    call utmess('F', 'POSTELEM_10')
                end if
                if (iret .ne. 0) goto 999
                call jeveuo(knum, 'L', jordr)
                call rsexch('F', result, 'DEPL', zi(jordr), chdef, &
                            iret)
                call chpve2(chdef, 3, tabtyp, ier)
            end if
        end if
        call pemain(tablOut, &
                    model, materField, materCode, caraElem, numeHarm, &
                    nbFactorKeyword, chdef)
    end if
!
    call getfac('ENER_POT', nbFactorKeyword)
    if (nbFactorKeyword .ne. 0) then
        call medomp(result, model, materField, materCode, caraElem, numeHarm)
        call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
        ASSERT(.not. isParallelMesh(mesh))
        call peepot(tablOut, &
                    model, materField, materCode, caraElem, numeHarm, &
                    nbFactorKeyword)
    end if
!
    call getfac('ENER_CIN', nbFactorKeyword)
    if (nbFactorKeyword .ne. 0) then
        call medomp(result, model, materField, materCode, caraElem, numeHarm)
        call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
        ASSERT(.not. isParallelMesh(mesh))
        call peecin(tablOut, &
                    model, materField, materCode, caraElem, numeHarm, &
                    nbFactorKeyword)
!
    end if
!
    call getfac('INTEGRALE', nbFactorKeyword)
    if (nbFactorKeyword .ne. 0) then
        call medomp(result, model)
        call peeint(tablOut, model, nbFactorKeyword)
    end if
!
    call getfac('NORME', nbFactorKeyword)
    if (nbFactorKeyword .ne. 0) then
!         --- ON RECUPERE LE MODELE
        call getvid('NORME', 'CHAM_GD', iocc=1, scal=chdef, nbret=n1)
        if (n1 .ne. 0) then
            call getvid('NORME', 'MODELE', iocc=1, scal=model, nbret=n2)
        else
            call getvid('NORME', 'RESULTAT', iocc=1, scal=result, nbret=nr)
            call medomp(result, model)
        end if
        call penorm(tablOut, model)
    end if
!
    call getfac('VOLUMOGRAMME', nbFactorKeyword)
    if (nbFactorKeyword .ne. 0) then
        call medomp(result, model, carele=caraElem)
        call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
        ASSERT(.not. isParallelMesh(mesh))
        call pevolu(tablOut, model, caraElem, nbFactorKeyword)
    end if
!
    call getfac('MINMAX', nbFactorKeyword)
    if (nbFactorKeyword .ne. 0) then
        call getvid('MINMAX', 'CHAM_GD', iocc=1, scal=chdef, nbret=n1)
        if (n1 .ne. 0) then
            call getvid('MINMAX', 'MODELE', iocc=1, scal=model, nbret=n2)
        else
            call getvid('MINMAX', 'RESULTAT', iocc=1, scal=result, nbret=nr)
            call medomp(result, model)
        end if
        call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
        ! ASSERT(.not. isParallelMesh(mesh))
        call pemima(n1, chdef, tablOut, model, nbFactorKeyword)
    end if
!
    call getfac('WEIBULL', nbFactorKeyword)
    if (nbFactorKeyword .ne. 0) then
        call medomp(result, model, materField, materCode, caraElem, numeHarm)
        call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
        ASSERT(.not. isParallelMesh(mesh))
        call peweib(tablOut, model, materField, materCode, caraElem, k8b, &
                    numeHarm, nbFactorKeyword, 0, nomcmd)
    end if
!
    call getfac('RICE_TRACEY', nbFactorKeyword)
    if (nbFactorKeyword .ne. 0) then
        call medomp(result, model, nh=numeHarm)
        call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
        ASSERT(.not. isParallelMesh(mesh))
        call peritr(tablOut, model, numeHarm, nbFactorKeyword)
    end if
!
    call getfac('CARA_GEOM', nbFactorKeyword)
    if (nbFactorKeyword .ne. 0) then
        call medomp(result, model)
        call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
        ASSERT(.not. isParallelMesh(mesh))
        call pecage(tablOut, model, nbFactorKeyword)
    end if
!
    call getfac('CARA_POUTRE', nbFactorKeyword)
    if (nbFactorKeyword .ne. 0) then
        call medomp(result, model, carele=caraElem, nh=numeHarm)
        call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
        ASSERT(.not. isParallelMesh(mesh))
        call pecapo(tablOut, model, numeHarm)
    end if
!
    call getfac('INDIC_ENER', nbFactorKeyword)
    if (nbFactorKeyword .ne. 0) then
        call medomp(result, model, materField, materCode, caraElem, numeHarm)
        call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
        ASSERT(.not. isParallelMesh(mesh))
        call peingl(tablOut, model, materField, materCode, caraElem, numeHarm, &
                    nbFactorKeyword, 'INDIC_ENER')
    end if
!
    call getfac('INDIC_SEUIL', nbFactorKeyword)
    if (nbFactorKeyword .ne. 0) then
        call medomp(result, model, materField, materCode, caraElem, numeHarm)
        call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
        ASSERT(.not. isParallelMesh(mesh))
        call peingl(tablOut, model, materField, materCode, caraElem, numeHarm, &
                    nbFactorKeyword, 'INDIC_SEUIL')
    end if
!
    call getfac('ENER_ELAS', nbFactorKeyword)
    if (nbFactorKeyword .ne. 0) then
        call medomp(result, model, materField, materCode, caraElem, numeHarm)
        call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
        ASSERT(.not. isParallelMesh(mesh))
        call peingl(tablOut, model, materField, materCode, caraElem, numeHarm, &
                    nbFactorKeyword, 'ENER_ELAS')
    end if
!
    call getfac('ENER_ELTR', nbFactorKeyword)
    if (nbFactorKeyword .ne. 0) then
        call medomp(result, model, materField, materCode, caraElem, numeHarm)
        call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
        ASSERT(.not. isParallelMesh(mesh))
        call peingl(tablOut, model, materField, materCode, caraElem, numeHarm, &
                    nbFactorKeyword, 'ENER_ELTR')
    end if

!
    call getfac('ENER_TOTALE', nbFactorKeyword)
    if (nbFactorKeyword .ne. 0) then
        call medomp(result, model, materField, materCode, caraElem, numeHarm)
        call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
        ASSERT(.not. isParallelMesh(mesh))
        call peingl(tablOut, model, materField, materCode, caraElem, numeHarm, &
                    nbFactorKeyword, 'ENER_TOTALE')
    end if
!
    call getfac('ENER_DISS', nbFactorKeyword)
    if (nbFactorKeyword .ne. 0) then
        call medomp(result, model, materField, materCode, caraElem, numeHarm)
        call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
        ASSERT(.not. isParallelMesh(mesh))
        call peingl(tablOut, model, materField, materCode, caraElem, numeHarm, &
                    nbFactorKeyword, 'ENER_DISS')
    end if
!
999 continue
    call titre()
!
    call jedema()
!
end subroutine
