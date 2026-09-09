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
subroutine dilcar(option, compor, icontm, ivarim, ideplm, ideplp, &
                  igeom, imatuu, ivectu, icontp, &
                  ivarip, ichg, ichn, jcret, jvCarcri, iinstm, iinstp)
!
    use Behaviour_module, only: behaviourOption
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/ismaem.h"
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
!
    integer(kind=8) :: icontm, ivarim, ideplm, ideplp, igeom, jcret
    integer(kind=8) :: imatuu, ivectu, icontp, ichg, ichn, ivarip, jvCarcri, iinstm, iinstp
    character(len=16) :: option
    character(len=16), pointer :: compor(:)
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! RECUPERATION DES ADRESSES DES CHAMPS DE L'ELEMENT POUR LES MODELES SECOND GRADIENT
!
! --------------------------------------------------------------------------------------------------
!
    aster_logical :: lVect, lMatr, lVari, lSigm
!
! --------------------------------------------------------------------------------------------------
!
    lVect = ASTER_FALSE
    lMatr = ASTER_FALSE
    lSigm = ASTER_FALSE
    lVari = ASTER_FALSE
    icontm = ismaem()
    ivarim = ismaem()
    ideplm = ismaem()
    ideplp = ismaem()
    imatuu = ismaem()
    ivectu = ismaem()
    icontp = ismaem()
    ivarip = ismaem()
    ichg = ismaem()
    ichn = ismaem()
    jcret = ismaem()
    jvCarcri = ismaem()
    iinstm = ismaem()
    iinstp = ismaem()
!
! - Input fields
!
    if (option(1:9) .eq. 'RIGI_MECA') then
        call jevech('PCONTMR', 'L', icontm)
        call jevech('PVARIMR', 'L', ivarim)
        call jevech('PDEPLMR', 'L', ideplm)
        call jevech('PDEPLPR', 'L', ideplp)
        call jevech('PGEOMER', 'L', igeom)
        call jevech('PCARCRI', 'L', jvCarcri)
        call jevech('PINSTMR', 'L', iinstm)
        call jevech('PINSTPR', 'L', iinstp)
    else if (option .eq. 'RAPH_MECA') then
        call jevech('PCONTMR', 'L', icontm)
        call jevech('PVARIMR', 'L', ivarim)
        call jevech('PDEPLMR', 'L', ideplm)
        call jevech('PDEPLPR', 'L', ideplp)
        call jevech('PGEOMER', 'L', igeom)
        call jevech('PCARCRI', 'L', jvCarcri)
        call jevech('PINSTMR', 'L', iinstm)
        call jevech('PINSTPR', 'L', iinstp)
    else if (option(1:9) .eq. 'FULL_MECA') then
        call jevech('PCONTMR', 'L', icontm)
        call jevech('PVARIMR', 'L', ivarim)
        call jevech('PDEPLMR', 'L', ideplm)
        call jevech('PDEPLPR', 'L', ideplp)
        call jevech('PGEOMER', 'L', igeom)
        call jevech('PCARCRI', 'L', jvCarcri)
        call jevech('PINSTMR', 'L', iinstm)
        call jevech('PINSTPR', 'L', iinstp)
    else if (option .eq. 'FORC_NODA') then
        call jevech('PSIEFR', 'L', icontm)
        call jevech('PGEOMER', 'L', igeom)
    else if (option .eq. 'REFE_FORC_NODA') then
        call jevech('PGEOMER', 'L', igeom)
    else
        ASSERT(ASTER_FALSE)
    end if

! - Select objects to construct from option name
    call jevech('PCOMPOR', 'L', vk16=compor)
    call behaviourOption(option, compor, &
                         lMatr, lVect, &
                         lVari, lSigm)

! - Output fields
    if ((option .eq. 'FORC_NODA') .or. (option .eq. 'REFE_FORC_NODA')) then
        call jevech('PVECTUR', 'E', ivectu)
    else
        if (lMatr) then
            call jevech('PMATUNS', 'E', imatuu)
        end if
        if (lVect) then
            call jevech('PVECTUR', 'E', ivectu)
        end if
        if (lSigm) then
            call jevech('PCONTPR', 'E', icontp)
            call jevech('PCODRET', 'E', jcret)
        end if
        if (lVari) then
            call jevech('PVARIPR', 'E', ivarip)
        end if
    end if
!
end subroutine
