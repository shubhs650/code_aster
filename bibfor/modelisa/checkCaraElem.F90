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
subroutine checkCaraElem(modelZ, caraElemZ)
!
    use coorSyst_module, only: setOrieFields
    implicit none
!
#include "asterf_types.h"
#include "asterfort/calcul.h"
#include "asterfort/dismoi.h"
#include "asterfort/mecact.h"
#include "asterfort/nmiret.h"
#include "asterfort/pipeElem_type.h"
#include "asterfort/setStructFields.h"
#include "asterfort/utmess.h"
!
    character(len=*), intent(in) :: modelZ, caraElemZ
!
! --------------------------------------------------------------------------------------------------
!
! AFFE_CARA_ELEM
!
! Some checks for structural elements
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : name of the model
! In  caraElem         : name of elementary characteristics (field)
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: option = "VERI_CARA_ELEM"
    integer(kind=8), parameter :: nbFieldInMax = 100, nbFieldOut = 2
    character(len=8) :: lpain(nbFieldInMax), lpaout(nbFieldOut)
    character(len=19) :: lchin(nbFieldInMax), lchout(nbFieldOut)
!
    integer(kind=8) :: nbFieldIn
    character(len=8) :: model, caraElem, mesh, answer
    character(len=24) :: modelLigrel
    character(len=24), parameter :: paraCheck = "&&OP0019.PARACHECK"
    integer(kind=8), parameter :: nbCmp = 2
    character(len=8), parameter :: cmpName(nbCmp) = (/"X1", "X2"/)
    real(kind=8) :: cmpVale(nbCmp)
    character(len=24), parameter :: codret = '&&OP0019.CODRET', indicr = '&&OP0019.INDICR'
    aster_logical :: tabret(0:10), lparallel_mesh
!
! --------------------------------------------------------------------------------------------------
!
    model = modelZ
    caraElem = caraElemZ
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
    call dismoi('PARALLEL_MESH', mesh, 'MAILLAGE', repk=answer)
    lparallel_mesh = answer .eq. 'OUI'

    if (.not. lparallel_mesh) then

! ----- Create map for parameters
        cmpVale(1) = PIPE_METRIC_LIMIT
        cmpVale(2) = 1.d0
        call mecact('V', paraCheck, 'LIGREL', modelLigrel, 'NEUT_R', &
                    ncmp=nbCmp, lnomcmp=cmpName, vr=cmpVale)

! ----- Input fields
        lpain(1) = 'PCHCKPR'
        lchin(1) = paraCheck(1:19)
        nbFieldIn = 1

! ----- Add fields for structural elements
        call setStructFields(caraElem, nbFieldInMax, lchin, lpain, nbFieldIn)

! ----- Add fields for orientation
        call setOrieFields(nbFieldInMax, lpain, lchin, &
                           nbFieldIn, caraElem)

! ----- Output fields
        lchout(1) = codret(1:19)
        lpaout(1) = 'PCODRET'
        lchout(2) = indicr(1:19)
        lpaout(2) = 'PINDICR'

! ----- Computation
        call calcul('C', option, modelLigrel, &
                    nbFieldIn, lchin, lpain, &
                    nbFieldOut, lchout, lpaout, &
                    'V', 'OUI')

! ----- Return code
        call nmiret(codret, tabret)
        if (tabret(1)) then
            call utmess('F', 'CALCULEL2_31')
        end if
        if (tabret(2)) then
            call utmess('A', 'PIPE1_54', sr=PIPE_METRIC_LIMIT)
        end if
    end if
!
! --------------------------------------------------------------------------------------------------
!
end subroutine
