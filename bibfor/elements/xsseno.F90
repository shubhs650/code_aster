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
subroutine xsseno(nno, nbsig, nse, npg, jgano, &
                  jsigpg, siseno)
    implicit none
!
!    - FONCTIONS REALISEES :  ROUTINE X-FEM
!
!     CALCUL DES CONTRAINTES PAR SOUS-ELEMENTS AUX NOEUDS (SENO)
!
! ......................................................................
!
!
!
!
#include "MeshTypes_type.h"
#include "jeveux.h"
#include "asterfort/ppgan2.h"
    integer(kind=8) :: mxval
    parameter(mxval=32*10*6)
!     EN 2D :
!     MXVAL =  6 (NBSE MAX) * 6 (NBNOSE MAX) * 4 (NBCMP MAX)
!     EN 3D :
!     MXVAL = 32 (NBSE MAX) * 10 (NBNOSE MAX) * 6 (NBCMP MAX)
!
    integer(kind=8) :: nno, npg, jgano
    integer(kind=8) :: nbsig
    integer(kind=8) :: jsigpg
    integer(kind=8) :: idecpg
    integer(kind=8) :: nse, ise, in, kpg, ic
!
    real(kind=8) :: vpg(15), vno(MT_NNOMAX3D)
!
    real(kind=8) :: siseno(mxval)
!
!
! ----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!     CALCUL DES CONTRAINTES PAR SOUS-ELEMENTS AUX NOEUDS (SENO)
!-----------------------------------------------------------------------
!
!   BOUCLE SUR LES NSE SOUS-ELEMENTS
    do ise = 1, nse
!
!       DEBUT DE LA ZONE MEMOIRE DE SIG  CORRESPONDANTE
        idecpg = npg*(ise-1)
!
!       BOUCLE NCMP DES CONTRAINTES
        do ic = 1, nbsig
!
            do kpg = 1, npg
                vpg(kpg) = zr(jsigpg+(kpg-1+idecpg)*nbsig+ic-1)
            end do
!
            call ppgan2(jgano, 1, 1, vpg, vno)
!
            do in = 1, nno
                siseno(nbsig*nno*(ise-1)+nbsig*(in-1)+ic) = vno(in)
            end do
!
        end do
!
    end do
!
!
end subroutine
