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
subroutine op0038()
!
    use HHO_precalc_module, only: hhoAddInputField
    use coorSyst_module, only: setOrieFields
    implicit none
!
#include "asterc/getres.h"
#include "asterc/r8vide.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/cesvar.h"
#include "asterfort/chpver.h"
#include "asterfort/copisd.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exlima.h"
#include "asterfort/getvid.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/mecham.h"
#include "asterfort/mechti.h"
#include "asterfort/rcmfmc.h"
#include "asterfort/sdmpic.h"
#include "asterfort/setStructFields.h"
#include "asterfort/vrcins.h"
#include "jeveux.h"
!
! --------------------------------------------------------------------------------------------------
!
! Command: CALC_CHAM_ELEM
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbFieldInMax = 100, nbFieldOut = 1
    character(len=8) :: lpaout(nbFieldOut), lpain(nbFieldInMax)
    character(len=19) :: lchout(nbFieldOut), lchin(nbFieldInMax)
    character(len=1), parameter :: jvBase = 'G'
    integer(kind=8) :: ierd, iret, numeHarm, nbRet, nbFieldIn
    real(kind=8) :: time, rundf
    character(len=2) :: chdret
    character(len=8) :: model, caraElem, temp, mesh, kmpic, chmate
    character(len=16) :: type, oper, option, phenom
    character(len=19) :: chelem, press, ligrel
    character(len=24) :: chgeom, chharm, materCode
    character(len=24) :: chtemp, chtime, chpres
    character(len=24), parameter :: chflug = '&&OP0038.FLUXGAUSS'
    character(len=24), parameter :: chvarc = '&&OP0038.CHVARC'
    aster_logical :: exitim, l_ther
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infmaj()

! - Initializations
    rundf = r8vide()

! - Output result
    call getres(chelem, type, oper)

! - Get main parameters
    model = ' '
    materCode = ' '
    caraElem = ' '
    chmate = ' '
    chtemp = ' '
    chpres = ' '
    call getvid(' ', 'MODELE', scal=model, nbret=nbRet)
    ASSERT(nbRet .eq. 1)
    call dismoi('PHENOMENE', model, 'MODELE', repk=phenom)
    l_ther = phenom .eq. 'THERMIQUE'
    call getvid(' ', 'CARA_ELEM', scal=caraElem, nbret=nbRet)
    call getvid(' ', 'CHAM_MATER', scal=chmate, nbret=nbRet)
    materCode = ' '
    if (nbRet .ne. 0) then
        call rcmfmc(chmate, materCode, l_ther_=l_ther)
    end if
    call getvtx(' ', 'OPTION', scal=option, nbret=nbRet)
    temp = ' '
    call getvid(' ', 'TEMP', scal=temp, nbret=nbRet)
    if (nbRet .ne. 0) then
        chtemp = temp
        call chpver('F', chtemp, 'NOEU', 'TEMP_R', ierd)
    end if
    press = ' '
    call getvid(' ', 'PRES', scal=press, nbret=nbRet)
    if (nbRet .ne. 0) then
        chpres = press
        call chpver('F', chpres, 'NOEU', 'PRES_C', ierd)
    end if
    call getvr8(' ', 'INST', scal=time, nbret=nbRet)
    exitim = nbRet .ne. 0
    call getvis(' ', 'MODE_FOURIER', scal=numeHarm, nbret=nbRet)
    if (nbRet .eq. 0) then
        numeHarm = 0
    end if
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)

! - List of cells for computation: all model
    call exlima(' ', 0, jvBase, model, ligrel)

! - Prepare input field
    call mecham(option, model, numeHarm, &
                chgeom, chharm, iret)

! - Create field for time
    chtime = ' '
    if (exitim) then
        call mechti(mesh, time, rundf, rundf, chtime)
    end if

    if (iret .ne. 0) goto 10

! - Compute
    if (option(1:7) .eq. 'FLUX_EL') then
        lchin(1) = chgeom(1:19)
        lpain(1) = 'PGEOMER'
        lchin(2) = materCode(1:19)
        lpain(2) = 'PMATERC'
        lchin(5) = chtemp(1:19)
        lpain(5) = 'PTEMPER'
        lchin(6) = chtime(1:19)
        lpain(6) = 'PINSTR'
        lchin(7) = chharm(1:19)
        lpain(7) = 'PHARMON'
        nbFieldIn = 7

! ----- Add fields for orientation
        call setOrieFields(nbFieldInMax, lpain, lchin, &
                           nbFieldIn, caraElem)

! ----- Add HHO field
        call hhoAddInputField(model, nbFieldInMax, lchin, lpain, nbFieldIn)

        lchout(1) = chflug(1:19)
        lpaout(1) = 'PFLUXPG'
        call calcul('S', 'FLUX_ELGA', ligrel, &
                    nbFieldIn, lchin, lpain, &
                    nbFieldOut, lchout, lpaout, &
                    jvBase, 'OUI')

        if (option .eq. 'FLUX_ELNO') then
            lchin(1) = chflug(1:19)
            lpain(1) = 'PFLUXPG'
            lchout(1) = chelem(1:19)
            lpaout(1) = 'PFLUXNO'
            nbFieldIn = 1
            call calcul('S', option, ligrel, &
                        nbFieldIn, lchin, lpain, &
                        nbFieldOut, lchout, lpaout, &
                        jvBase, 'OUI')

        else if (option .eq. 'FLUX_ELGA') then
            call copisd('CHAMP', jvBase, chflug, chelem)

        else
            ASSERT(ASTER_FALSE)

        end if

    else if (option .eq. 'COOR_ELGA') then
        lchin(1) = chgeom(1:19)
        lpain(1) = 'PGEOMER'
        nbFieldIn = 1

! ----- Add fields for orientation
        call setOrieFields(nbFieldInMax, lpain, lchin, &
                           nbFieldIn, caraElem)

! ----- Add fields for structural elements
        call setStructFields(caraElem, nbFieldInMax, lchin, lpain, nbFieldIn)
        lchout(1) = chelem(1:19)
        lpaout(1) = 'PCOORPG'
        call cesvar(caraElem, ' ', ligrel, lchout(1))
        call calcul('S', option, ligrel, &
                    nbFieldIn, lchin, lpain, &
                    nbFieldOut, lchout, lpaout, &
                    jvBase, 'OUI')

    else if (option .eq. 'ROCH_ELNO') then
        call vrcins(model, chmate, caraElem, time, chvarc(1:19), &
                    chdret)

        lchin(1) = materCode(1:19)
        lpain(1) = 'PMATERC'
        lchin(2) = chvarc(1:19)
        lpain(2) = 'PVARCPR'
        nbFieldIn = 2
! ----- Add fields for structural elements
        call setStructFields(caraElem, nbFieldInMax, lchin, lpain, nbFieldIn)

        lchout(1) = chelem
        lpaout(1) = 'PROCHRR'
        call calcul('S', option, ligrel, &
                    nbFieldIn, lchin, lpain, &
                    nbFieldOut, lchout, lpaout, &
                    jvBase, 'OUI')
        call detrsd('CHAMP_GD', chvarc)

    else if (option .eq. 'PRAC_ELNO') then
        lpain(1) = 'PPRESSC'
        lchin(1) = chpres(1:19)
        lchout(1) = chelem(1:19)
        lpaout(1) = 'PPRAC_R'
        nbFieldIn = 1
        call calcul('S', option, ligrel, &
                    nbFieldOut, lchin, lpain, &
                    nbFieldOut, lchout, lpaout, &
                    jvBase, 'OUI')

    else
        ASSERT(ASTER_FALSE)

    end if
!
10  continue

! - SI CHELEM N'EST PAS MPI_COMPLET, ON LE COMPLETE :
    call dismoi('MPI_COMPLET', chelem, 'CHAM_ELEM', repk=kmpic)
    if (kmpic .eq. 'NON') call sdmpic('CHAM_ELEM', chelem)
!
    call jedema()
end subroutine
