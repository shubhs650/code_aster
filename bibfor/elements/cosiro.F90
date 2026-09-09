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
subroutine cosiro(plateCara, plateOrie, &
                  paraNameZ, loue, sens, goun, &
                  jtens_)
!
    use plateGeom_module, only: isShell3D, isPlate
    use plate_type
    implicit none
!
#include "asterc/r8dgrd.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/dxsiro.h"
#include "asterfort/tecach.h"
#include "asterfort/vdsiro.h"
#include "asterfort/plate_type.h"
#include "jeveux.h"
!
    type(plateCara_Para), intent(in) :: plateCara
    type(plateOrie_Para), intent(in) :: plateOrie
    character(len=*), intent(in) :: paraNameZ
    character(len=1), intent(in) :: loue
    character(len=2), intent(in) :: sens
    character(len=1), intent(in) :: goun
    integer(kind=8), optional, intent(out) :: jtens_
!
! --------------------------------------------------------------------------------------------------
!
! CHANGER LE REPERE : INTRINSEQUE <-> UTILISATEUR  POUR UN CHAMP LOCAL
!       (MODELISATIONS : DKT/DST/Q4G/COQUE_3D)
!
! --------------------------------------------------------------------------------------------------
!
!  PARAM  IN : NOM DU CHAMP LOCAL A MODIFIER
!  LOUE   IN :  / 'L' : PARAMETRE EN LECTURE
!               / 'E' : PARAMETRE EN ECRITURE
!  SENS   IN :  / 'IU' : INTRINSEQUE -> UTILISATEUR
!               / 'UI' : UTILISATEUR -> INTRINSEQUE
!  GOUN   IN :  / 'G' : CHAMP ELGA
!               / 'N' : CHAMP ELNO
!  JTENS  OUT : ADRESSE DU CHAMP LOCAL (QUE L'ON A MODIFIE)
!
!  REMARQUE : UTILISER SOUR='R' PEUT FAIRE GAGNER UN PEU DE TEMPS MAIS
!             CELA PERMET SURTOUT DE SE PROTEGER DES TE00IJ QUI
!             MODIFIENT LA GEOMETRIE INITIALE (EX : TE0031)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: itab(7), iret, nbpt, nbsp
    integer(kind=8) :: jtens
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(loue .eq. 'L' .or. loue .eq. 'E')
    ASSERT(sens .eq. 'UI' .or. sens .eq. 'IU')
    ASSERT(goun .eq. 'G' .or. goun .eq. 'N')

! - Access to field (tensor)
    call tecach('NNO', paraNameZ, loue, iret, nval=7, itab=itab)
    ASSERT(iret .eq. 0 .or. iret .eq. 1)

! - Frame is OK !
    ASSERT(plateOrie%lUpdate)

    if (iret .ne. 1) then
! ----- Parameters of field
        jtens = itab(1)
        nbpt = itab(3)
        nbsp = itab(7)

! ----- Apply local coordinate system on tensor
        if (isShell3D(plateCara)) then
            if (goun .eq. 'G') then
                call vdsiro(nbpt, nbsp, plateOrie%matevg, sens, goun, &
                            zr(jtens), zr(jtens))
            else
                call vdsiro(nbpt, nbsp, plateOrie%matevn, sens, goun, &
                            zr(jtens), zr(jtens))
            end if

        elseif (isPlate(plateCara)) then
            if (sens .eq. 'UI') then
                call dxsiro(nbpt*nbsp, plateOrie%t2ui, zr(jtens), zr(jtens))
            else
                call dxsiro(nbpt*nbsp, plateOrie%t2iu, zr(jtens), zr(jtens))
            end if
        else
            ASSERT(ASTER_FALSE)
        end if
    end if

    if (present(jtens_)) then
        jtens_ = jtens
    end if
end subroutine
