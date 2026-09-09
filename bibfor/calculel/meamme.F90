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
subroutine meamme(modelZ, &
                  materFieldZ, materCodeZ, caraElemZ, &
                  time, jvBaseZ, &
                  matrRigiZ, matrMassZ, &
                  matrElemz, &
                  variZ, comporZ, sddyna)
!
    use coorSyst_module, only: setOrieFields
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mecact.h"
#include "asterfort/mecham.h"
#include "asterfort/memare.h"
#include "asterfort/ndynkk.h"
#include "asterfort/reajre.h"
#include "asterfort/redetr.h"
#include "asterfort/setStructFields.h"
#include "asterfort/vrcins.h"
!
    character(len=*), intent(in) :: modelZ
    character(len=*), intent(in) :: materFieldZ, materCodeZ, caraElemZ
    real(kind=8), intent(in) :: time
    character(len=*), intent(in) :: jvBaseZ
    character(len=*), intent(in) :: matrRigiZ, matrMassZ, matrElemz
    character(len=*), intent(in) :: variZ, comporZ
    character(len=19), intent(in) :: sddyna
!
! --------------------------------------------------------------------------------------------------
!
! Elementary matrix for AMOR_MECA / RIGI_MECA_HYST
!
! NB: careful, compute Dirichlet [B] matrix too when RIGI_MECA_HYST
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : option to compute
! In  model            : name of model
! In  nbLoad           : number of loads
! In  listLoadK24      : pointer to the name of loads
! In  materField       : name of material characteristics (field)
! In  materCode        : name of coded material
! In  caraElem         : name of elementary characteristics (field)
! In  time             : current time
! In  jvBase           : JEVEUX base to create matrElem
! In  matrRigi         : elementary rigidity matrix
! In  matrMass         : elementary rigidity mass
! In  listElemCalc     : list of element (LIGREL) where matrElem is computed
! In  matrElem         : elementary matrix
! In  numeHarm         : index of Fourier mode
! In  vari             : internal state variables
! In  compor           : field of behaviour (non-linear cases)
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: option = 'AMOR_MECA'
    integer(kind=8), parameter :: nbFieldInMax = 100, nbFieldOutMax = 2
    character(len=8) :: lpain(nbFieldInMax), lpaout(nbFieldOutMax)
    character(len=19) :: lchin(nbFieldInMax), lchout(nbFieldOutMax)
!
    integer(kind=8) :: nbFieldIn, nbFieldOut
    character(len=2) :: codret
    integer(kind=8) :: iret
    integer(kind=8), parameter :: numeHarm = 0
    character(len=24), parameter :: chvarc = '&&MEAMME.CHVARC'
    character(len=24) :: compor, vari
    character(len=8) :: physQuantityName
    character(len=24) :: matrRigi, matrMass
    character(len=24) :: resuElemRigi, resuElemMass
    character(len=24) :: chgeom, chharm
    character(len=1) :: jvBase
    character(len=8) :: model, caraElem, mesh
    character(len=24) :: materField, materCode, amor_flui
    character(len=19) :: matrElem
    integer(kind=8) :: nbResuElem, iResuElem, idxResuElemRigi
    integer(kind=8) :: nbSubstruct
    character(len=24), pointer :: listResuElem(:) => null()
    character(len=19) :: modelLigrel, resuLigrel
    character(len=24), parameter :: nonLinearMap = "&&MEAMMA.NONLIN"
    integer(kind=8), parameter :: nbCmp = 1
    character(len=8), parameter :: cmpName = ('X1')
    integer(kind=8), parameter :: cmpVale = 1
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    model = modelZ
    caraElem = caraElemZ
    materField = materFieldZ
    materCode = materCodeZ
    matrElem = matrElemz
    jvBase = jvBaseZ
    matrRigi = matrRigiZ
    matrMass = matrMassZ
    compor = comporZ
    vari = variZ
    lpain = ' '
    lchin = ' '
    lpaout = ' '
    lchout = ' '

! - Get parameters
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)
    call dismoi('NB_SS_ACTI', model, 'MODELE', repi=nbSubstruct)
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)

! - Preparation of input fields
    call mecham(option, model, numeHarm, &
                chgeom, chharm, iret)

! - Special map for non-linear cases
    call jedetr(nonLinearMap)
    call mecact('V', nonLinearMap, 'MAILLA', mesh, 'NEUT_I', &
                ncmp=nbCmp, nomcmp=cmpName, si=cmpVale)

! - Special map for fluid damping
    call ndynkk(sddyna, 'AMOR_FLUI', amor_flui)

! - Field for external state variables
    call vrcins(model, materField, caraElem, time, chvarc, codret)

! - Get RESU_ELEM from rigidity matrix
    resuElemRigi = ' '
    idxResuElemRigi = 0
    if (matrRigi(1:1) .ne. ' ') then
        call jeexin(matrRigi(1:19)//'.RELR', iret)
        if (iret .gt. 0) then
            call jeveuo(matrRigi(1:19)//'.RELR', 'L', vk24=listResuElem)
            call jelira(matrRigi(1:19)//'.RELR', 'LONUTI', nbResuElem)
            do iResuElem = 1, nbResuElem
                resuElemRigi = listResuElem(iResuElem)
                idxResuElemRigi = iResuElem
                call dismoi('NOM_LIGREL', resuElemRigi, 'RESUELEM', repk=resuLigrel)
                if (resuLigrel .eq. modelLigrel) then
                    goto 20
                end if
            end do
            ASSERT(ASTER_FALSE)
20          continue
        end if
    end if

! - Get RESU_ELEM from mass matrix
    resuElemMass = ' '
    if (matrMass(1:1) .ne. ' ') then
        call jeexin(matrMass(1:19)//'.RELR', iret)
        if (iret .gt. 0) then
            call jeveuo(matrMass(1:19)//'.RELR', 'L', vk24=listResuElem)
            call jelira(matrMass(1:19)//'.RELR', 'LONUTI', nbResuElem)
            do iResuElem = 1, nbResuElem
                resuElemMass = listResuElem(iResuElem)
                call dismoi('NOM_LIGREL', resuElemMass, 'RESUELEM', repk=resuLigrel)
                if (resuLigrel .eq. modelLigrel) then
                    goto 40
                end if
            end do
            ASSERT(ASTER_FALSE)
40          continue
        end if
    end if

! - Prepare RESU_ELEM objects
    call memare(jvBase, matrElem, model, 'AMOR_MECA', to_aster_logical(nbSubstruct > 0))
    call jedetr(matrElem//'.RELR')

! - Add input fields
    lpain(1) = 'PGEOMER'
    lchin(1) = chgeom(1:19)
    lpain(2) = 'PMATERC'
    lchin(2) = materCodeZ(1:19)
    lpain(3) = 'PMASSEL'
    lchin(3) = resuElemMass(1:19)
    lpain(4) = 'PCOMPOR'
    lchin(4) = compor(1:19)
    lpain(5) = 'PNONLIN'
    lchin(5) = nonLinearMap(1:19)
    lpain(6) = 'PVARIPG'
    lchin(6) = vari(1:19)
    lpain(7) = 'PAMORFL'
    lchin(7) = amor_flui(1:19)
    lpain(8) = 'PVARCPR'
    lchin(8) = chvarc(1:19)
    nbFieldIn = 8

! - Add fields for structural elements
    call setStructFields(caraElem, nbFieldInMax, lchin, lpain, nbFieldIn)

! - Add fields for orientation
    call setOrieFields(nbFieldInMax, lpain, lchin, &
                       nbFieldIn, caraElem)

! - Get symmetric or unsymmetric rigidity matrix
    if (resuElemRigi .ne. ' ') then
        nbFieldIn = nbFieldIn+1
        lchin(nbFieldIn) = resuElemRigi(1:19)
        call dismoi('NOM_GD', resuElemRigi, 'RESUELEM', repk=physQuantityName)
        if (physQuantityName .eq. 'MDNS_R') then
            lpain(nbFieldIn) = 'PRIGINS'
        else
            lpain(nbFieldIn) = 'PRIGIEL'
            call jeveuo(matrRigi(1:19)//'.RELR', 'L', vk24=listResuElem)
            call jelira(matrRigi(1:19)//'.RELR', 'LONUTI', nbResuElem)
            if (idxResuElemRigi .lt. nbResuElem) then
                resuElemRigi = listResuElem(idxResuElemRigi+1)
                call dismoi('NOM_GD', resuElemRigi, 'RESUELEM', repk=physQuantityName)
                if (physQuantityName .eq. 'MDNS_R') then
                    nbFieldIn = nbFieldIn+1
                    lpain(nbFieldIn) = 'PRIGINS'
                    lchin(nbFieldIn) = resuElemRigi(1:19)
                end if
            end if
        end if
    end if

! - Output fields
    lpaout(1) = 'PMATUUR'
    lpaout(2) = 'PMATUNS'
    lchout(1) = matrElem(1:8)//'.ME001'
    lchout(2) = matrElem(1:8)//'.ME002'
    nbFieldOut = 2

! - Compute
    ASSERT(nbFieldIn .le. nbFieldInMax)
    ASSERT(nbFieldOut .le. nbFieldOutMax)
    call calcul('S', &
                option, modelLigrel, &
                nbFieldIn, lchin, lpain, &
                nbFieldOut, lchout, lpaout, &
                jvBase, 'OUI')

! - Save RESU_ELEM
    call reajre(matrElem, lchout(1), jvBase)
    call reajre(matrElem, lchout(2), jvBase)

! - Clean
    call redetr(matrElem)
    call detrsd('CHAMP_GD', chvarc)
!
    call jedema()
end subroutine
