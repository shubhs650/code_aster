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
subroutine majour(neq, lgrot, sdnume, chaini, &
                  chadel, coef, chamaj, ordre)
!
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/nmgrot.h"
    character(len=19) :: sdnume
    aster_logical :: lgrot
    integer(kind=8) :: neq, ordre
    real(kind=8) :: chaini(*), chadel(*), chamaj(*), coef
!
! ----------------------------------------------------------------------
!
! ROUTINE MECA_NON_LINE (ALGORITHME - UTILITAIRE - STATIQUE)
!
! MET A JOUR LES CHAM_NO DES DEPLACEMENTS
!
! ----------------------------------------------------------------------
!
!
! CHAMAJ = CHAINI + COEF*CHADEL.
!   POUR LES TRANSLATIONS ET LES PETITES ROTATIONS, ON APPLIQUE
!   LA FORMULE PRECEDENTE A LA LETTRE.
!   POUR LES GRANDES ROTATIONS, LE VECTEUR-ROTATION DE CHAMAJ
!   EST CELUI DU PRODUIT DE LA ROTATION DEFINIE DANS CHAINI PAR
!   COEF FOIS L'INCREMENT DE ROTATION DEFINI DANS CHADEL.
!
! IN  NEQ    : LONGUEUR DES CHAM_NO
! IN  SDNUME : SD NUMEROTATION
! IN  LGROT  : TRUE  S'IL Y A DES DDL DE GRDE ROTATION
!                       FALSE SINON
! IN  CHAINI : CHAM_NO DONNE
! IN  CHADEL : CHAM_NO DONNE
! IN  COEF   : REEL DONNE
! IN  ORDRE  : 0 -> MAJ INCREMENTS
!              1 -> MAJ DEPL
! OUT CHAMAJ : CHAM_NO MIS A JOUR
!
!
!
!
    integer(kind=8) :: iran(3), i, icomp
    real(kind=8) :: theta(3), deldet(3)
    integer(kind=8) :: ptdo, indic1, indic2
    integer(kind=8), pointer :: ndro(:) => null()
!
! ----------------------------------------------------------------------
!
    ptdo = 0
    indic1 = 0
    indic2 = 0
!
    call jemarq()
!
    if (lgrot) then
        call jeveuo(sdnume//'.NDRO', 'L', vi=ndro)
    end if
!
    if (.not. lgrot) then
        do i = 1, neq
            chamaj(i) = chaini(i)+coef*chadel(i)
        end do
    else
        icomp = 0
        do i = 1, neq
            if (ndro(i) .eq. 0) then
                chamaj(i) = chaini(i)+coef*chadel(i)
            else if (ndro(i) .eq. 1) then
                icomp = icomp+1
                iran(icomp) = i
                theta(icomp) = chaini(i)
                deldet(icomp) = coef*chadel(i)
                if (icomp .eq. 3) then
                    icomp = 0
                    call nmgrot(iran, deldet, theta, chamaj)
                end if
            else
                ASSERT(.false.)
            end if
        end do
    end if
!
    call jedema()
end subroutine
