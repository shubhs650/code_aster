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
subroutine merimp(lXFEM, lDyna, &
                  model, caraElem, sddyna, iterNewt, &
                  ds_constitutive, ds_material, &
                  hval_incr, hval_algo, caco3d, &
                  nbFieldInMax, lpain, lchin, nbFieldIn)
!
    use NonLin_Datastructure_type
    use HHO_precalc_module, only: hhoAddInputField
    use coorSyst_module, only: setOrieFields
    implicit none
!
#include "asterf_types.h"
#include "asterfort/cesvar.h"
#include "asterfort/copisd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/mecact.h"
#include "asterfort/mecara.h"
#include "asterfort/megeom.h"
#include "asterfort/ndynkk.h"
#include "asterfort/nmchex.h"
#include "asterfort/nmvcex.h"
#include "asterfort/setStructFields.h"
#include "asterfort/vtzero.h"
#include "asterfort/xajcin.h"
!
    aster_logical, intent(in) :: lXFEM, lDyna
    character(len=24), intent(in) :: model, caraElem
    character(len=19), intent(in) :: sddyna
    integer(kind=8), intent(in) :: iterNewt
    type(NL_DS_Constitutive), intent(in) :: ds_constitutive
    type(NL_DS_Material), intent(in) :: ds_material
    character(len=19), intent(in) :: hval_incr(*), hval_algo(*)
    character(len=24), intent(in) :: caco3d
    integer(kind=8), intent(in) :: nbFieldInMax
    character(len=8), intent(inout) :: lpain(nbFieldInMax)
    character(len=19), intent(inout) :: lchin(nbFieldInMax)
    integer(kind=8), intent(out) :: nbFieldIn
!
! --------------------------------------------------------------------------------------------------
!
! Nonlinear mechanics (algorithm)
!
! Computation of rigidity matrix and internal forces - Input fields
!
! --------------------------------------------------------------------------------------------------
!
! In  lXFEM            : flag for XFEM elements
! In  lDyna            : flag for dynamic
! In  l_hho            : flag for HHO elements
! In  model            : name of model
! In  caraElem         : name of elementary characteristics (field)
! In  sddyna           : datastructure for dynamic
! In  iterNewt         : index of current Newton iteration
! In  ds_constitutive  : datastructure for constitutive laws management
! In  ds_material      : datastructure for material parameters
! In  hval_incr        : hat-variable for incremental values fields
! In  hval_algo        : hat-variable for algorithms fields
! In  hhoField         : datastructure for HHO method
! In  caco3d           : name of field for COQUE_3D (field of normals)
! In  nbFieldInMax     : maximum number of input fields
! IO  lpain            : list of input parameters
! IO  lchin            : list of input fields
! Out nbFieldIn        : number of input fields
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24), parameter :: chiter = '&&MERIMO.CH_ITERAT'
    character(len=24), parameter :: variIter = '&&MERIMO.VARMOJ', strxIter = '&&MERIMO.STRMOJ'
    character(len=16), parameter :: option = 'FULL_MECA'
    integer(kind=8) :: iret
    character(len=24) :: chgeom
    character(len=19) :: stadyn, depent, vitent
    character(len=19) :: dispPrev, sigmPrev, variPrev, varcPrev, strxPrev
    character(len=19) :: dispCurr, sigmCurr, variCurr, varcCurr, strxCurr
    character(len=19) :: timePrev, varcAllPrev
    character(len=19) :: timeCurr, varcAllCurr
    character(len=24) :: varcRefe
    character(len=19) :: depkm1, vitkm1, acckm1
    character(len=19) :: viteCurr, acceCurr, vitePrev, accePrev
    character(len=19) :: romkm1, romk
    character(len=24) :: modelLigrel
    character(len=19) :: dispIter, dispCumuInst
!
! --------------------------------------------------------------------------------------------------
!
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)
    nbFieldIn = 0

! - Get fields from hat-variables - Begin of time step
    call nmchex(hval_incr, 'VALINC', 'DEPMOI', dispPrev)
    call nmchex(hval_incr, 'VALINC', 'VITMOI', vitePrev)
    call nmchex(hval_incr, 'VALINC', 'ACCMOI', accePrev)
    call nmchex(hval_incr, 'VALINC', 'SIGMOI', sigmPrev)
    call nmchex(hval_incr, 'VALINC', 'VARMOI', variPrev)
    call nmchex(hval_incr, 'VALINC', 'COMMOI', varcPrev)
    call nmchex(hval_incr, 'VALINC', 'STRMOI', strxPrev)

! - Get fields from hat-variables - End of time step
    call nmchex(hval_incr, 'VALINC', 'DEPPLU', dispCurr)
    call nmchex(hval_incr, 'VALINC', 'VITPLU', viteCurr)
    call nmchex(hval_incr, 'VALINC', 'ACCPLU', acceCurr)
    call nmchex(hval_incr, 'VALINC', 'SIGPLU', sigmCurr)
    call nmchex(hval_incr, 'VALINC', 'VARPLU', variCurr)
    call nmchex(hval_incr, 'VALINC', 'COMPLU', varcCurr)
    call nmchex(hval_incr, 'VALINC', 'STRPLU', strxCurr)
!
    call nmchex(hval_incr, 'VALINC', 'DEPKM1', depkm1)
    call nmchex(hval_incr, 'VALINC', 'VITKM1', vitkm1)
    call nmchex(hval_incr, 'VALINC', 'ACCKM1', acckm1)
    call nmchex(hval_incr, 'VALINC', 'ROMKM1', romkm1)
    call nmchex(hval_incr, 'VALINC', 'ROMK  ', romk)
!
    call nmchex(hval_algo, 'SOLALG', 'DEPDEL', dispCumuInst)
    call nmchex(hval_algo, 'SOLALG', 'DDEPLA', dispIter)

! - Dynamic fields
    if (lDyna) then
        call ndynkk(sddyna, 'DEPENT', depent)
        call ndynkk(sddyna, 'VITENT', vitent)
        call ndynkk(sddyna, 'STADYN', stadyn)
    end if

! - Get external state variables
    call nmvcex('TOUT', varcPrev, varcAllPrev)
    call nmvcex('INST', varcPrev, timePrev)
    call nmvcex('TOUT', varcCurr, varcAllCurr)
    call nmvcex('INST', varcCurr, timeCurr)
    call nmvcex('TOUT', ds_material%varc_refe, varcRefe)

! - Get internal state variables from previous iteration
    call exisd('CHAMP_GD', variCurr, iret)
    if (iret .ne. 0) then
        call copisd('CHAMP_GD', 'V', variCurr, variIter)
    else
        call copisd('CHAMP_GD', 'V', variPrev, variIter)
    end if

! - Get structural variables from previous iteration
    call exisd('CHAMP_GD', strxPrev, iret)
    if (iret .ne. 0 .and. iterNewt .lt. 2) then
        call copisd('CHAMP_GD', 'V', strxPrev, strxIter)
        call vtzero(strxIter, 'CHAM_ELEM')
    elseif (iterNewt .ge. 2) then
        call copisd('CHAMP_GD', 'V', strxCurr, strxIter)
    end if

! - Extend elementary field for internal variables
    call exisd('CHAM_ELEM_S', ds_constitutive%compor, iret)
    if (iret .eq. 0) then
        call cesvar(caraElem, ds_constitutive%compor, modelLigrel, ds_constitutive%compor)
    end if
    call copisd('CHAM_ELEM_S', 'V', ds_constitutive%compor, variCurr)
    call copisd('CHAM_ELEM_S', 'V', ds_constitutive%compor, sigmCurr)

! - Get geometry field
    call megeom(model, chgeom)

! - Create field for iteration number
    call mecact('V', chiter, 'MODELE', modelLigrel, 'NEUT_I', &
                ncmp=1, nomcmp='X1', si=iterNewt)

! - Set input fields
    lpain(1) = 'PGEOMER'
    lchin(1) = chgeom(1:19)
    lpain(2) = 'PMATERC'
    lchin(2) = ds_material%mateco(1:19)
    lpain(3) = 'PCONTMR'
    lchin(3) = sigmPrev(1:19)
    lpain(4) = 'PVARIMR'
    lchin(4) = variPrev(1:19)
    lpain(5) = 'PCOMPOR'
    lchin(5) = ds_constitutive%compor(1:19)
    lpain(6) = 'PDEPLMR'
    lchin(6) = dispPrev(1:19)
    lpain(7) = 'PDEPLPR'
    lchin(7) = dispCumuInst(1:19)
    lpain(8) = 'PINSTMR'
    lchin(8) = timePrev(1:19)
    lpain(9) = 'PINSTPR'
    lchin(9) = timeCurr(1:19)
    lpain(10) = 'PCARCRI'
    lchin(10) = ds_constitutive%carcri(1:19)
    lpain(11) = 'PITERAT'
    lchin(11) = chiter(1:19)
    lpain(12) = 'PDDEPLA'
    lchin(12) = dispIter(1:19)
    lpain(13) = 'PDEPKM1'
    lchin(13) = depkm1(1:19)
    lpain(14) = 'PVITKM1'
    lchin(14) = vitkm1(1:19)
    lpain(15) = 'PACCKM1'
    lchin(15) = acckm1(1:19)
    lpain(16) = 'PROMKM1'
    lchin(16) = romkm1(1:19)
    lpain(17) = 'PROMK'
    lchin(17) = romk(1:19)
    lpain(18) = 'PVARIMP'
    lchin(18) = variIter(1:19)
    lpain(19) = 'PVARCMR'
    lchin(19) = varcAllPrev(1:19)
    lpain(20) = 'PVARCPR'
    lchin(20) = varcAllCurr(1:19)
    lpain(21) = 'PVARCRR'
    lchin(21) = varcRefe(1:19)
    lpain(22) = 'PCACO3D'
    lchin(22) = caco3d(1:19)
    lpain(23) = 'PSTRXMR'
    lchin(23) = strxPrev(1:19)
    lpain(24) = 'PSTRXMP'
    lchin(24) = strxIter(1:19)
    lpain(25) = 'PMULCOM'
    lchin(25) = ds_constitutive%mult_comp(1:19)
    nbFieldIn = 25

! - Add fields for structural elements
    call setStructFields(caraElem, nbFieldInMax, lchin, lpain, nbFieldIn)

! - Add fields for XFEM
    if (lXFEM) then
        call xajcin(model, option, nbFieldInMax, lchin, lpain, nbFieldIn)
    end if

! - Add fields for dynamic
    if (lDyna) then
        nbFieldIn = nbFieldIn+1
        lpain(nbFieldIn) = 'PDEPENT'
        lchin(nbFieldIn) = depent(1:19)
        nbFieldIn = nbFieldIn+1
        lpain(nbFieldIn) = 'PVITENT'
        lchin(nbFieldIn) = vitent(1:19)
        nbFieldIn = nbFieldIn+1
        lpain(nbFieldIn) = 'PSTADYN'
        lchin(nbFieldIn) = stadyn(1:19)
        nbFieldIn = nbFieldIn+1
        lpain(nbFieldIn) = 'PVITPLU'
        lchin(nbFieldIn) = viteCurr(1:19)
        nbFieldIn = nbFieldIn+1
        lpain(nbFieldIn) = 'PACCPLU'
        lchin(nbFieldIn) = acceCurr(1:19)
        nbFieldIn = nbFieldIn+1
        lpain(nbFieldIn) = 'PVITMOI'
        lchin(nbFieldIn) = vitePrev(1:19)
        nbFieldIn = nbFieldIn+1
        lpain(nbFieldIn) = 'PACCMOI'
        lchin(nbFieldIn) = accePrev(1:19)
    end if

! - Add fields for HHO
    call hhoAddInputField(model, nbFieldInMax, lchin, lpain, nbFieldIn)

! - Add fields for orientation
    call setOrieFields(nbFieldInMax, lpain, lchin, &
                       nbFieldIn, caraElem)
!
end subroutine
