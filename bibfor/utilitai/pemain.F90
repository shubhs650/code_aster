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
subroutine pemain(tablOutZ, &
                  modelZ, materFieldZ, materCodeZ, caraElemZ, numeHarm, &
                  nbFactorKeyword, deformZ)
!
    use coorSyst_module, only: setOrieFields
    implicit none
!
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/char8_to_int.h"
#include "asterfort/detrsd.h"
#include "asterfort/exlim3.h"
#include "asterfort/getvem.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/mecham.h"
#include "asterfort/pemica.h"
#include "asterfort/setStructFields.h"
#include "asterfort/tbajli.h"
#include "asterfort/tbajpa.h"
#include "asterfort/tbcrsd.h"
#include "asterfort/umalma.h"
#include "asterfort/utmess.h"
#include "asterfort/vtgpld.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
!
    character(len=*), intent(in) :: tablOutZ
    character(len=*), intent(in) :: modelZ, materFieldZ, materCodeZ, caraElemZ
    integer(kind=8), intent(in) :: numeHarm, nbFactorKeyword
    character(len=*), intent(in) :: deformZ
!
! --------------------------------------------------------------------------------------------------
!
! POST_ELEM
!
! TRAITEMENT DU MOT CLE-FACTEUR "MASS_INER"
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: option = "MASS_INER"
    integer(kind=8), parameter :: nbFieldInMax = 100, nbFieldOut = 1
    character(len=8) :: lpain(nbFieldInMax), lpaout(nbFieldOut)
    character(len=19) :: lchin(nbFieldInMax), lchout(nbFieldOut)
!
    integer(kind=8) :: nbFieldIn
    integer(kind=8) :: ibid, iret, iFactorKeyword, nt, ng, nr, nm, nbgrma, jgr, ig, nbma, jad
    integer(kind=8) :: nbmail, jma, im, nume, ifm, niv, iorig
    integer(kind=8) :: icage, nbtot, nbMaiT, nre
    integer(kind=8), parameter :: nbValeR1 = 16, nbValeR2 = 25
    integer(kind=8) :: nbValeR
    real(kind=8) :: zero, orig(3), r8b
    character(len=8) :: k8b, mesh, valk(2)
    character(len=19), parameter :: chelem = '&&PEMAIN.MASS_INER'
    character(len=24), parameter :: chgeo2 = '&&PEMAIN.CH_GEOMER'
    character(len=19) :: chdef
    character(len=24) :: valk2(2)
    character(len=24) :: chgeom, chharm, ligrel
    complex(kind=8) :: c16b
    real(kind=8), pointer :: trav1(:) => null()
    integer(kind=8), pointer :: v_allma(:) => null()
!
    integer(kind=8), parameter :: nbParaResu1 = 18, nbParaResu2 = 27
    integer(kind=8) :: nbParaResu
    integer(kind=8), parameter :: nbParaResuMax = 27
    character(len=16), parameter :: tablParaName(nbParaResuMax) = &
                                    (/'LIEU     ', 'ENTITE   ', 'MASSE    ', 'CDG_X    ', &
                                      'CDG_Y    ', 'CDG_Z    ', 'IX_G     ', 'IY_G     ', &
                                      'IZ_G     ', 'IXY_G    ', 'IXZ_G    ', 'IYZ_G    ', &
                                      'IX_PRIN_G', 'IY_PRIN_G', 'IZ_PRIN_G', 'ALPHA    ', &
                                      'BETA     ', 'GAMMA    ', 'X_P      ', 'Y_P      ', &
                                      'Z_P      ', 'IX_P     ', 'IY_P     ', 'IZ_P     ', &
                                      'IXY_P    ', 'IXZ_P    ', 'IYZ_P    '/)
    character(len=8), parameter :: tablParaType(nbParaResuMax) = &
                                   (/'K24', 'K8 ', 'R  ', 'R  ', &
                                     'R  ', 'R  ', 'R  ', 'R  ', &
                                     'R  ', 'R  ', 'R  ', 'R  ', &
                                     'R  ', 'R  ', 'R  ', 'R  ', &
                                     'R  ', 'R  ', 'R  ', 'R  ', &
                                     'R  ', 'R  ', 'R  ', 'R  ', &
                                     'R  ', 'R  ', 'R  '/)
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infniv(ifm, niv)

! - Initializations
    ibid = 0
    c16b = (0.d0, 0.d0)
    icage = 0
    zero = 0.0d0
    r8b = 0.0d0
    chdef = deformZ
    lpain = ' '
    lchin = ' '
    lpaout = ' '
    lchout = ' '

! - Check and prepare input fields
    call mecham(option, modelZ, numeHarm, &
                chgeom, chharm, iret)

    if (iret .ne. 0) goto 60
    mesh = chgeom(1:8)
!
    call exlim3(option, 'V', modelZ, ligrel)

! - Add input fields
    lpain(1) = 'PGEOMER'
    if (chdef .ne. ' ') then
        call vtgpld('CUMU', 1.d0, chgeom, chdef, 'V', chgeo2)
        lchin(1) = chgeo2(1:19)
    else
        lchin(1) = chgeom(1:19)
    end if
    lpain(2) = 'PMATERC'
    lchin(2) = materCodeZ
    lpain(3) = 'PCOMPOR'
    lchin(3) = materFieldZ(1:8)//'.COMPOR'
    nbFieldIn = 3

! - Add fields for structural elements
    call setStructFields(caraElemZ, nbFieldInMax, lchin, lpain, nbFieldIn)

! - Add fields for orientation
    call setOrieFields(nbFieldInMax, lpain, lchin, &
                       nbFieldIn, caraElemZ)

! - Set output field
    lpaout(1) = 'PMASSINE'
    lchout(1) = chelem
!
    call calcul('S', option, ligrel, &
                nbFieldIn, lchin, lpain, &
                nbFieldOut, lchout, lpaout, &
                'V', 'OUI')

! - Prepare list of parameters in table
    nbValeR = nbValeR1
    nbParaResu = nbParaResu1
    do iFactorKeyword = 1, nbFactorKeyword
        call getvr8(option, 'ORIG_INER', iocc=iFactorKeyword, nbval=0, nbret=nr)
        if (nr .ne. 0) then
            nbValeR = nbValeR2
            nbParaResu = nbParaResu2
            exit
        end if
    end do

! - Create output table
    call tbcrsd(tablOutZ, 'G')
    call tbajpa(tablOutZ, nbParaResu, tablParaName, tablParaType)
!
    AS_ALLOCATE(vr=trav1, size=nbValeR)
    do iFactorKeyword = 1, nbFactorKeyword
        iorig = 0
        orig = zero
        call getvtx(option, 'TOUT', iocc=iFactorKeyword, nbval=0, nbret=nt)
        call getvem(mesh, 'GROUP_MA', option, 'GROUP_MA', iFactorKeyword, &
                    0, k8b, ng)
        call getvem(mesh, 'MAILLE', option, 'MAILLE', iFactorKeyword, &
                    0, k8b, nm)

        call getvr8(option, 'ORIG_INER', iocc=iFactorKeyword, nbval=0, nbret=nr)
        if (nr .ne. 0) then
            iorig = 1
            nre = -nr
            call getvr8(option, 'ORIG_INER', iocc=iFactorKeyword, nbval=nre, vect=orig, &
                        nbret=nr)
        end if
        if (nt .ne. 0) then
            call pemica(chelem, nbValeR, trav1, 0, [ibid], &
                        orig, iorig, icage)
            valk(1) = mesh
            valk(2) = 'TOUT'
            call tbajli(tablOutZ, nbParaResu, tablParaName, [ibid], trav1, &
                        [c16b], valk, 0)
        end if
        if (ng .ne. 0) then
            nbgrma = -ng
            call wkvect('&&PEMAIN_GROUPM', 'V V K24', nbgrma, jgr)
            call getvem(mesh, 'GROUP_MA', option, 'GROUP_MA', iFactorKeyword, &
                        nbgrma, zk24(jgr), ng)
            valk2(2) = 'GROUP_MA'
            do ig = 1, nbgrma
                call jeexin(jexnom(mesh//'.GROUPEMA', zk24(jgr+ig-1)), iret)
                if (iret .eq. 0) then
                    call utmess('A', 'UTILITAI3_46', sk=zk24(jgr+ig-1))
                    cycle
                end if
                call jelira(jexnom(mesh//'.GROUPEMA', zk24(jgr+ig-1)), 'LONUTI', nbma)
                if (nbma .eq. 0) then
                    call utmess('A', 'UTILITAI3_47', sk=zk24(jgr+ig-1))
                    cycle
                end if
                call jeveuo(jexnom(mesh//'.GROUPEMA', zk24(jgr+ig-1)), 'L', jad)
                call pemica(chelem, nbValeR, trav1, nbma, zi(jad), &
                            orig, iorig, icage)
                valk2(1) = zk24(jgr+ig-1)
                call tbajli(tablOutZ, nbParaResu, tablParaName, [ibid], trav1, &
                            [c16b], valk2, 0)
            end do
!
!
! --- UNION
            if (nbgrma > 1) then
                call umalma(mesh, zk24(jgr), nbgrma, v_allma, nbtot)
                ASSERT(nbtot > 0)
                !
                call pemica(chelem, nbValeR, trav1, nbtot, v_allma, orig, iorig, icage)
                valk2(1) = "UNION_GROUP_MA"
                call tbajli(tablOutZ, nbParaResu, tablParaName, [ibid], trav1, [c16b], valk2, 0)
                !
                AS_DEALLOCATE(vi=v_allma)
            end if
            call jedetr('&&PEMAIN_GROUPM')
        end if
        if (nm .ne. 0) then
            nbmail = -nm
            call wkvect('&&PEMAIN_MAILLE', 'V V K8', nbmail, jma)
            call getvem(mesh, 'MAILLE', option, 'MAILLE', iFactorKeyword, &
                        nbmail, zk8(jma), nm)
            valk(2) = 'MAILLE'
            call jelira(mesh//'.TYPMAIL', 'LONMAX', nbMaiT)
            do im = 1, nbmail
                nume = char8_to_int(zk8(jma+im-1))
                if ((nume .gt. nbMaiT) .or. (nume .le. 0)) then
                    call utmess('A', 'UTILITAI3_49', sk=zk8(jma+im-1))
                    cycle
                end if
                call pemica(chelem, nbValeR, trav1, 1, [nume], &
                            orig, iorig, icage)
                valk(1) = zk8(jma+im-1)
                call tbajli(tablOutZ, nbParaResu, tablParaName, [ibid], trav1, &
                            [c16b], valk, 0)
            end do
            call jedetr('&&PEMAIN_MAILLE')
        end if
    end do

! - MENAGE
    call detrsd('CHAM_ELEM', chelem)
    call detrsd('CHAMP_GD', chgeo2)
    AS_DEALLOCATE(vr=trav1)
!
60  continue
!
    call jedema()
end subroutine
