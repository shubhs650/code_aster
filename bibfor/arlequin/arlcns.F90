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
subroutine arlcns(nummai, connex, loncum, nbno, cxno)
!
!
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "MeshTypes_type.h"
!
    integer(kind=8) :: nummai, connex(*), loncum(*)
    integer(kind=8) :: nbno
    integer(kind=8) :: cxno(nbno)
!
! ----------------------------------------------------------------------
!
! CONSTRUCTION DE BOITES ENGLOBANTES POUR UN GROUPE DE MAILLES
!
! RETOURNE LES NUMEROS ABSOLUS DES NOEUDS DE LA CONNECTIVITE DE LA
! MAILLE POUR LES ELEMENTS SOLIDES
!
! ----------------------------------------------------------------------
!
!
! IN  NUMMAI : NUMERO ABSOLU DE LA MAILLE DANS LE MAILLAGE
! IN  CONNEX : CONNEXITE DES MAILLES
! IN  LONCUM : LONGUEUR CUMULEE DE CONNEX
! IN  NBNO   : NOMBRE DE NOEUDS DE LA MAILLE
! OUT CXNO   : CONNECTIVITE DE LA MAILLE
!
! ----------------------------------------------------------------------
!
    integer(kind=8) :: ino, jdec
!
! ----------------------------------------------------------------------
!
    jdec = loncum(nummai)
!
    if ((nbno < 1) .or. (nbno > MT_NNOMAX)) then
        ASSERT(.false.)
    end if
!
    do ino = 1, nbno
        cxno(ino) = connex(jdec-1+ino)
    end do
!
end subroutine
