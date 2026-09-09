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
subroutine ther_mtan(l_stat, &
                     modelZ, caraElemZ, materCodeZ, &
                     timePara, varcCurrZ, &
                     comporTherZ, tempIterZ, &
                     resuElemZ, matrElemZ, jvBase)
!
    use coorSyst_module, only: setOrieFields
    implicit none
!
#include "asterf_types.h"
#include "asterfort/calcul.h"
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
    character(len=*), intent(inout) :: resuElemZ
    character(len=*), intent(in) :: matrElemZ
    character(len=1), intent(in) :: jvBase
!
! --------------------------------------------------------------------------------------------------
!
! Thermic
!
! Tangent matrix (volumic terms)
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
! IO  resuElem         : name of resu_elem
! In  matrElem         : name of matr_elem result
! In  jvBase           : JEVEUX base for object
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: optionRigi = 'RIGI_THER_TANG', optionMass = 'MASS_THER_TANG'
    integer(kind=8), parameter :: nbFieldInMax = 100, nbFieldOut = 1
    character(len=8) :: lpain(nbFieldInMax), lpaout(nbFieldOut)
    character(len=24) :: lchin(nbFieldInMax), lchout(nbFieldOut)
    character(len=24) :: modelLigrel, chgeom
    character(len=19) :: resuElem
    real(kind=8) :: theta, deltat
    integer(kind=8) :: nbFieldIn
    character(len=8) :: newnom
!
! --------------------------------------------------------------------------------------------------
!
    call dismoi('NOM_LIGREL', modelZ, 'MODELE', repk=modelLigrel)
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

! - Add output field
    lpaout(1) = 'PMATTSR'
    lchout(1) = resuElemZ

! - Compute rigidity term
    call calcul("S", optionRigi, modelLigrel, &
                nbFieldIn, lchin, lpain, &
                nbFieldOut, lchout, lpaout, &
                jvBase, 'OUI')

! - Multiply values by theta
    call multResuElem(resuElem, theta)

! - Add RESU_ELEM in MATR_ELEM
    call reajre(matrElemZ, resuElem, jvBase)

! - Compute mass term
    if (.not. l_stat) then
! - --- Output fields
        newnom = resuElem(9:16)
        call gcnco2(newnom)
        resuElem(10:16) = newnom(2:8)
        lpaout(1) = 'PMATTTR'
        lchout(1) = resuElem

! - --- Compute
        call calcul("S", optionMass, modelLigrel, &
                    nbFieldIn, lchin, lpain, &
                    nbFieldOut, lchout, lpaout, &
                    jvBase, 'OUI')

! - --- Multiply values by 1/dt
        call multResuElem(resuElem, 1.d0/deltat)

! - --- Add RESU_ELEM in MATR_ELEM
        call reajre(matrElemZ, resuElem, jvBase)

    end if
!
    resuElemZ = resuElem
!
end subroutine
