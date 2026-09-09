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
subroutine merimo(jvBase, &
                  lXFEM, lMacrElem, &
                  model, caraElem, iterNewt, &
                  ds_constitutive, ds_material, ds_system, &
                  hval_incr, hval_algo, &
                  optionZ, ldccvg, sddynaZ)
!
    use NonLin_Datastructure_type
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
#include "asterfort/jemarq.h"
#include "asterfort/memare.h"
#include "asterfort/merimp.h"
#include "asterfort/ndynlo.h"
#include "asterfort/nmchex.h"
#include "asterfort/nmiret.h"
#include "asterfort/reajre.h"
#include "asterfort/redetr.h"
#include "asterfort/vemare.h"
#include "jeveux.h"
!
    character(len=1), intent(in) :: jvBase
    aster_logical, intent(in) :: lXFEM, lMacrElem
    character(len=24), intent(in) :: model, caraElem
    integer(kind=8), intent(in) :: iterNewt
    type(NL_DS_Constitutive), intent(in) :: ds_constitutive
    type(NL_DS_Material), intent(in) :: ds_material
    type(NL_DS_System), intent(in) :: ds_system
    character(len=19), intent(in) :: hval_incr(*), hval_algo(*)
    character(len=*), intent(in) :: optionZ
    integer(kind=8), intent(out) :: ldccvg
    character(len=*), optional, intent(in) :: sddynaZ
!
! --------------------------------------------------------------------------------------------------
!
! Nonlinear mechanics (algorithm)
!
! Computation of rigidity matrix and internal forces
!
! --------------------------------------------------------------------------------------------------
!
! In  jvBase           : JEVEUX jvBase to create objects
! In  lXFEM            : flag for XFEM elements
! In  lMacrElem        : flag for macro-elements
! In  model            : name of model
! In  caraElem         : name of elementary characteristics (field)
! In  iterNewt         : index of current Newton iteration
! In  ds_constitutive  : datastructure for constitutive laws management
! In  ds_material      : datastructure for material parameters
! In  ds_system        : datastructure for non-linear system management
! In  hval_incr        : hat-variable for incremental values fields
! In  hval_algo        : hat-variable for algorithms fields
! In  option           : name of option to compute
! Out ldccvg           : return code from integration of behaviour
!                       -1 - no integration of behaviour
!                        0 - OK
!                        1 - Failure
!                        2 - Failure but not fatal during Newton's iterations
!                        3 - De Borst not converged
!                        4 - Using law outside bounds (VERI_BORNE)
! In  sddyna           : datastructure for dynamic
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbFieldOutMax = 11, nbFieldInMax = 100
    character(len=8) :: lpaout(nbFieldOutMax), lpain(nbFieldInMax)
    character(len=19) :: lchout(nbFieldOutMax), lchin(nbFieldInMax)
    aster_logical :: l_merigi, l_veinte, l_sigmex
    aster_logical :: l_codret, l_codpre, l_dyna
    integer(kind=8) :: ires, iret, nbFieldIn, nbFieldOut
    character(len=24) :: modelLigrel
    character(len=19) :: sigmExtr, sigmCurr, variCurr, strxCurr, sddyna
    character(len=16) :: option
    integer(kind=8) :: ich_matrixs, ich_matrixn, ich_veinte, ich_codret, ich_copred
    aster_logical :: tabret(0:10)
    character(len=24), parameter :: caco3d = '&&MERIMO.CARA_ROTAF'
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    option = optionZ
    sddyna = ' '
    if (present(sddynaZ)) then
        sddyna = sddynaZ
    end if
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelLigrel)
    tabret = ASTER_FALSE
    ldccvg = 0
    l_dyna = ndynlo(sddyna, 'DYNAMIQUE')
    lpain = " "
    lpaout = " "
    lchin = " "
    lchout = " "

! - Get fields from hat-variables
    call nmchex(hval_incr, 'VALINC', 'SIGEXT', sigmExtr)
    call nmchex(hval_incr, 'VALINC', 'SIGPLU', sigmCurr)
    call nmchex(hval_incr, 'VALINC', 'VARPLU', variCurr)
    call nmchex(hval_incr, 'VALINC', 'STRPLU', strxCurr)

! - Set input fields
    call merimp(lXFEM, l_dyna, &
                model, caraElem, sddyna, iterNewt, &
                ds_constitutive, ds_material, &
                hval_incr, hval_algo, caco3d, &
                nbFieldInMax, lpain, lchin, nbFieldIn)

! - Prepare flags
    if (option(1:9) .eq. 'FULL_MECA') then
        l_merigi = ASTER_TRUE
        l_veinte = ASTER_TRUE
        l_codret = ASTER_TRUE
        l_sigmex = ASTER_FALSE
        l_codpre = ASTER_FALSE
    else if (option(1:10) .eq. 'RIGI_MECA ') then
        l_merigi = ASTER_TRUE
        l_veinte = ASTER_FALSE
        l_codret = ASTER_FALSE
        l_sigmex = ASTER_FALSE
        l_codpre = ASTER_FALSE
    else if (option(1:16) .eq. 'RIGI_MECA_IMPLEX') then
        l_merigi = ASTER_TRUE
        l_veinte = ASTER_FALSE
        l_codret = ASTER_FALSE
        l_sigmex = ASTER_TRUE
        l_codpre = ASTER_FALSE
    else if (option(1:10) .eq. 'RIGI_MECA_') then
        l_merigi = ASTER_TRUE
        l_veinte = ASTER_TRUE
        l_codret = ASTER_FALSE
        l_sigmex = ASTER_FALSE
        l_codpre = ASTER_TRUE
        if (option .eq. 'RIGI_MECA_TANG') then
            call detrsd('CHAM_ELEM', ds_constitutive%comp_error)
            l_codret = ASTER_TRUE
        end if
    else if (option(1:9) .eq. 'RAPH_MECA') then
        l_merigi = ASTER_FALSE
        l_veinte = ASTER_TRUE
        l_codret = ASTER_TRUE
        l_sigmex = ASTER_FALSE
        l_codpre = ASTER_FALSE
    else
        ASSERT(ASTER_FALSE)
    end if

! - Prepare vector and matrix
    if (l_merigi) then
        call detrsd('MATR_ELEM', ds_system%merigi)
        call jeexin(ds_system%merigi//'.RERR', ires)
        if (ires .eq. 0) then
            call memare(jvBase, ds_system%merigi, model, 'RIGI_MECA', lMacrElem)
        end if
        call reajre(ds_system%merigi, ' ', jvBase)
    end if
!
    if (l_veinte) then
        call jeexin(ds_system%veinte//'.RELR', iret)
        if (iret .eq. 0) then
            call vemare(jvBase, ds_system%veinte, model)
        end if
        call jedetr(ds_system%veinte//'.RELR')
        call reajre(ds_system%veinte, ' ', jvBase)
    end if

! - Set output fields
    lpaout(1) = 'PVARIPR'
    lchout(1) = variCurr(1:19)
    lpaout(2) = 'PCACO3D'
    lchout(2) = caco3d(1:19)
    lpaout(3) = 'PSTRXPR'
    lchout(3) = strxCurr(1:19)
    nbFieldOut = 3
    if (l_merigi) then
        nbFieldOut = nbFieldOut+1
        lpaout(nbFieldOut) = 'PMATUUR'
        lchout(nbFieldOut) = ds_system%merigi(1:15)//'.M01'
        ich_matrixs = nbFieldOut
        nbFieldOut = nbFieldOut+1
        lpaout(nbFieldOut) = 'PMATUNS'
        lchout(nbFieldOut) = ds_system%merigi(1:15)//'.M02'
        ich_matrixn = nbFieldOut
    end if
    if (l_veinte) then
        nbFieldOut = nbFieldOut+1
        lpaout(nbFieldOut) = 'PVECTUR'
        lchout(nbFieldOut) = ds_system%veinte(1:15)//'.R01'
        ich_veinte = nbFieldOut
    end if
    if (l_sigmex) then
        nbFieldOut = nbFieldOut+1
        lpaout(nbFieldOut) = 'PCONTXR'
        lchout(nbFieldOut) = sigmExtr(1:19)
    end if
    if (l_codret) then
        nbFieldOut = nbFieldOut+1
        lpaout(nbFieldOut) = 'PCODRET'
        lchout(nbFieldOut) = ds_constitutive%comp_error(1:19)
        ich_codret = nbFieldOut
    end if
    if (l_codpre) then
        nbFieldOut = nbFieldOut+1
        lpaout(nbFieldOut) = 'PCOPRED'
        lchout(nbFieldOut) = ds_constitutive%code_pred(1:19)
        ich_copred = nbFieldOut
        nbFieldOut = nbFieldOut+1
        lpaout(nbFieldOut) = 'PCONTPR'
        lchout(nbFieldOut) = ds_constitutive%sigm_pred(1:19)
    else
        nbFieldOut = nbFieldOut+1
        lpaout(nbFieldOut) = 'PCONTPR'
        lchout(nbFieldOut) = sigmCurr(1:19)
    end if
!
    ASSERT(nbFieldOut .le. nbFieldOutMax)
    ASSERT(nbFieldIn .le. nbFieldInMax)
!
! - Compute
!
    call calcul('S', option, modelLigrel, nbFieldIn, lchin, &
                lpain, nbFieldOut, lchout, lpaout, jvBase, &
                'NON')
!
! - Save
!
    if (l_merigi) then
        call reajre(ds_system%merigi, lchout(ich_matrixs), jvBase)
        call reajre(ds_system%merigi, lchout(ich_matrixn), jvBase)
        call redetr(ds_system%merigi)
    end if
    if (l_veinte) then
        call reajre(ds_system%veinte, lchout(ich_veinte), jvBase)
    end if
!
! - Errors
!
    if (l_codret) then
        call nmiret(lchout(ich_codret), tabret)
        if (tabret(0)) then
            if (tabret(1)) then
                ldccvg = 1
            else if (tabret(3)) then
                ldccvg = 3
            else if (tabret(2)) then
                ldccvg = 2
            else if (tabret(4)) then
                ldccvg = 4
            else
                ldccvg = 1
            end if
        else
            ldccvg = 0
        end if
    end if
!
    call jedema()
end subroutine
