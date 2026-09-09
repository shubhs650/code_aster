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
subroutine varcCalcComp(modelZ, chsithz, &
                        l_temp, l_hydr, l_ptot, &
                        l_sech, l_epsa, &
                        nbFieldIn, nbFieldOut, &
                        lpain, lchin, &
                        lpaout, lchout, &
                        vectElemZ)
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/calcul.h"
#include "asterfort/corich.h"
#include "asterfort/dismoi.h"
#include "asterfort/gcnco2.h"
#include "asterfort/reajre.h"
!
    character(len=*), intent(in) :: modelZ, chsithz
    aster_logical, intent(in)  :: l_temp, l_hydr, l_ptot, l_sech, l_epsa
    integer(kind=8), intent(in) :: nbFieldIn, nbFieldOut
    character(len=8), intent(in) :: lpain(*), lpaout(*)
    character(len=19), intent(in) :: lchin(*)
    character(len=19), intent(inout) :: lchout(*)
    character(len=*), intent(in) :: vectElemZ
!
! --------------------------------------------------------------------------------------------------
!
! Material - External state variables (VARC)
!
! Calls to CALCUL
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : name of model
! In  chsith           : commande variable for temperature in XFEM
! In  l_temp           : .true. if temperature exists
! In  l_hydr           : .true. if hydration exists
! In  l_ptot           : .true. if total pressure (THM) exists
! In  l_sech           : .true. if drying exists
! In  l_epsa           : .true. if non-elastic strain exists
! In  nbFieldIn        : effective number of input fields
! In  nbFieldOut       : effective number of output fields
! In  lpain            : list of input parameters
! In  lchin            : list of input fields
! In  lpaout           : list of output parameters
! IO  lchout           : list of output fields
! In  vect_elem        : name of elementary vectors
!
! --------------------------------------------------------------------------------------------------
!
    character(len=1), parameter :: jvBase = "V"
    character(len=8) :: newnom
    character(len=16) :: option
    character(len=19) :: resuElem, modelLigrel
!
! --------------------------------------------------------------------------------------------------
!
    call dismoi('NOM_LIGREL', modelZ, 'MODELE', repk=modelLigrel)
    newnom = '.0000000'
    resuElem = vectElemZ(1:8)//'.0000000'

! - Temperature
    if (l_temp) then
        call gcnco2(newnom)
        resuElem(10:16) = newnom(2:8)
        call corich('E', resuElem, ichin_=-1)
        lchout(1) = resuElem
        option = 'CHAR_MECA_TEMP_R'
        if (nbFieldOut .eq. 2) then
            lchout(2) = chsithz(1:19)
        end if
        call calcul('C', option, modelLigrel, &
                    nbFieldIn, lchin, lpain, &
                    nbFieldOut, lchout, lpaout, &
                    jvBase, 'OUI')
        call reajre(vectElemZ, resuElem, jvBase)
    end if

! - Hydratation
    if (l_hydr) then
        call gcnco2(newnom)
        resuElem(10:16) = newnom(2:8)
        call corich('E', resuElem, ichin_=-1)
        lchout(1) = resuElem
        option = 'CHAR_MECA_HYDR_R'
        call calcul('C', option, modelLigrel, &
                    nbFieldIn, lchin, lpain, &
                    nbFieldOut, lchout, lpaout, &
                    jvBase, 'OUI')
        call reajre(vectElemZ, resuElem, jvBase)
    end if

! - Total pressure (THM)
    if (l_ptot) then
        call gcnco2(newnom)
        resuElem(10:16) = newnom(2:8)
        call corich('E', resuElem, ichin_=-1)
        lchout(1) = resuElem
        option = 'CHAR_MECA_PTOT_R'
        call calcul('C', option, modelLigrel, &
                    nbFieldIn, lchin, lpain, &
                    nbFieldOut, lchout, lpaout, &
                    jvBase, 'OUI')
        call reajre(vectElemZ, resuElem, jvBase)
    end if

! - Drying
    if (l_sech) then
        call gcnco2(newnom)
        resuElem(10:16) = newnom(2:8)
        call corich('E', resuElem, ichin_=-1)
        lchout(1) = resuElem
        option = 'CHAR_MECA_SECH_R'
        call calcul('C', option, modelLigrel, &
                    nbFieldIn, lchin, lpain, &
                    nbFieldOut, lchout, lpaout, &
                    jvBase, 'OUI')
        call reajre(vectElemZ, resuElem, jvBase)
    end if

! - Non-elastic strain
    if (l_epsa) then
        call gcnco2(newnom)
        resuElem(10:16) = newnom(2:8)
        call corich('E', resuElem, ichin_=-1)
        lchout(1) = resuElem
        option = 'CHAR_MECA_EPSA_R'
        call calcul('C', option, modelLigrel, &
                    nbFieldIn, lchin, lpain, &
                    nbFieldOut, lchout, lpaout, &
                    jvBase, 'OUI')
        call reajre(vectElemZ, resuElem, jvBase)
    end if
!
end subroutine
