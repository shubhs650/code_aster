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
subroutine resi_ther(l_stat, &
                     modelZ, caraElemZ, materCodeZ, &
                     timePara, timeMapZ, varcCurrZ, &
                     comporTherZ, tempIterZ, &
                     tempPrevZ, hydrPrevZ, hydrCurrZ, &
                     resuElemZ, vectElemZ, jvBase)
!
    use coorSyst_module, only: setOrieFields
    implicit none
!
#include "asterf_types.h"
#include "asterfort/calcul.h"
#include "asterfort/corich.h"
#include "asterfort/dismoi.h"
#include "asterfort/gcnco2.h"
#include "asterfort/megeom.h"
#include "asterfort/multResuElem.h"
#include "asterfort/reajre.h"
!
    aster_logical, intent(in) :: l_stat
    character(len=*), intent(in) :: modelZ, caraElemZ, materCodeZ
    real(kind=8), intent(in) :: timePara(2)
    character(len=*), intent(in) :: tempIterZ, comporTherZ, varcCurrZ
    character(len=*), intent(in) :: tempPrevZ, hydrPrevZ, hydrCurrZ, timeMapZ
    character(len=*), intent(inout) :: resuElemZ
    character(len=*), intent(in) :: vectElemZ
    character(len=1), intent(in) :: jvBase
!
! --------------------------------------------------------------------------------------------------
!
! Thermic
!
! Residuals from non-linear laws
!
! --------------------------------------------------------------------------------------------------
!
! In  l_stat           : flag for stationnary computation (no mass term)
! In  model            : name of the model
! In  caraElem         : name of elementary characteristics (field)
! In  materCode        : name of coding material characteristics (field)
! In  timePara         : timePara(1) = theta
!                        timePara(2) = deltat
! In  varcCurr         : command variable for current time
! In  comporTher       : name of comportment definition (field)
! In  tempIter         : temperature field at current Newton iteration
! In  tempPrev         : previous temperature
! In  hydrPrev         : previous hydration
! In  hydrCurr         : current hydration
! IO  resuElem         : name of resu_elem
! In  vectElem         : name of vect_elem result
! In  jvBase           : JEVEUX base for object
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: optionRigi = 'RAPH_THER', optionMass = 'MASS_THER_RESI'
    character(len=16), parameter :: optionHydr = "HYDR_ELGA"
    integer(kind=8), parameter :: nbFieldInMax = 100, nbFieldOutMax = 2
    character(len=8) :: lpain(nbFieldInMax), lpaout(nbFieldOutMax)
    character(len=24) :: lchin(nbFieldInMax), lchout(nbFieldOutMax)
    character(len=24) :: modelLigrel, chgeom
    character(len=19) :: resuElem
    real(kind=8) :: theta, deltat
    integer(kind=8) :: nbFieldIn, nbFieldOut
    character(len=8) :: newnom, answer
    aster_logical :: l_dry
!
! --------------------------------------------------------------------------------------------------
!
    call dismoi('NOM_LIGREL', modelZ, 'MODELE', repk=modelLigrel)
    call dismoi('EXI_SECH', modelZ, 'MODELE', repk=answer)
    l_dry = answer .eq. 'OUI'
    theta = timePara(1)
    deltat = timePara(2)
    lpain = " "
    lchin = " "
    lpaout = " "
    lchout = " "
    resuElem = resuElemZ(1:19)

! - Get geometry field
    call megeom(modelZ, chgeom)

! - Add input fields
    lpain(1) = 'PGEOMER'
    lchin(1) = chgeom
    lpain(2) = 'PMATERC'
    lchin(2) = materCodeZ
    lpain(3) = 'PTEMPEI'
    lchin(3) = tempIterZ
    lpain(4) = 'PCOMPOR'
    lchin(4) = comporTherZ
    lpain(5) = 'PVARCPR'
    lchin(5) = varcCurrZ
    nbFieldIn = 5

! - Add fields for orientation
    call setOrieFields(nbFieldInMax, lpain, lchin, &
                       nbFieldIn, caraElemZ)

! - Add output fields
    lpaout(1) = 'PRESIDU'
    lchout(1) = resuElemZ
    lpaout(2) = 'PFLUXPR'
    lchout(2) = "&&RESI_THER.FLUXPR"
    call corich('E', lchout(1), ichin_=-1)
    nbFieldOut = 2

! - Compute rigidity term
    call calcul("S", optionRigi, modelLigrel, &
                nbFieldIn, lchin, lpain, &
                nbFieldOut, lchout, lpaout, &
                jvBase, 'OUI')

! - Multiply values by theta
    call multResuElem(resuElem, theta)

! - Add RESU_ELEM in VECT_ELEM
    call reajre(vectElemZ, resuElem, jvBase)

! - Compute hydratation
    if (.not. l_stat .and. .not. l_dry) then
! ----- Input fields
        lpain = " "
        lchin = " "
        lpain(1) = 'PMATERC'
        lchin(1) = materCodeZ
        lpain(2) = 'PCOMPOR'
        lchin(2) = comporTherZ
        lpain(3) = 'PINSTR'
        lchin(3) = timeMapZ
        lpain(4) = 'PTEMPMR'
        lchin(4) = tempPrevZ
        lpain(5) = 'PTEMPPR'
        lchin(5) = tempIterZ
        lpain(6) = 'PHYDRMR'
        lchin(6) = hydrPrevZ
        lpain(7) = 'PGEOMER'
        lchin(7) = chgeom
        nbFieldIn = 7

! - --- Output fields
        lpaout(1) = 'PHYDRPR'
        lchout(1) = hydrCurrZ
        nbFieldOut = 1

! - --- Compute
        call calcul("S", optionHydr, modelLigrel, &
                    nbFieldIn, lchin, lpain, &
                    nbFieldOut, lchout, lpaout, &
                    jvBase, 'OUI')
    end if

! - Compute mass term
    if (.not. l_stat) then
! ----- Input fields
        lpain = " "
        lchin = " "
        lpain(1) = 'PGEOMER'
        lchin(1) = chgeom
        lpain(2) = 'PMATERC'
        lchin(2) = materCodeZ
        lpain(3) = 'PTEMPEI'
        lchin(3) = tempIterZ
        lpain(4) = 'PCOMPOR'
        lchin(4) = comporTherZ
        lpain(5) = 'PVARCPR'
        lchin(5) = varcCurrZ
        lpain(6) = 'PHYDRPR'
        lchin(6) = hydrCurrZ
        nbFieldIn = 6

! - --- Output fields
        newnom = resuElem(9:16)
        call gcnco2(newnom)
        resuElem(10:16) = newnom(2:8)
        lpaout(1) = 'PRESIDU'
        lchout(1) = resuElem
        call corich('E', lchout(1), ichin_=-1)
        nbFieldOut = 1

! - --- Compute
        call calcul("S", optionMass, modelLigrel, &
                    nbFieldIn, lchin, lpain, &
                    nbFieldOut, lchout, lpaout, &
                    jvBase, 'OUI')

! - --- Multiply values by 1/dt
        call multResuElem(resuElem, 1.d0/deltat)

! - --- Add RESU_ELEM in VECT_ELEM
        call reajre(vectElemZ, resuElem, jvBase)
    end if
!
    resuElemZ = resuElem
!
end subroutine
