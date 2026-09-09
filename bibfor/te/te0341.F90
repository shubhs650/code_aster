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

subroutine te0341(option, nomte)
    use resi_refe_module, only: RESI_REFE
!
!
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/cgfono.h"
#include "asterfort/cgfore.h"
#include "asterfort/cginit.h"
#include "asterfort/cgtang.h"
#include "asterfort/elref2.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/Behaviour_type.h"
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:  FORC_NODA ET REFE_FORC_NODA
!                          POUR ELEMENTS GAINE/CABLE
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
!      INSPIRE DE TE0361
! ......................................................................
!
    character(len=8) :: lielrf(10)
    integer(kind=8) :: nno1, nno2, npg, ivf2, idf2, nnos, jgn
    integer(kind=8) :: iw, ivf1, idf1, igeom, jvSief, ivectu, ndim, ntrou
    integer(kind=8) :: npgn, iwn, ivf1n, idf1n, jgnn, ino, i, nddl1
    integer(kind=8) :: iu(3, 3), iuc(3), im(3), isect, jvDisp
    character(len=16), pointer :: compor(:) => null()
    real(kind=8) :: tang(3, 3), forref, sigref, depref, a
    real(kind=8) :: geom(3, 3)
    aster_logical :: reactu
    type(RESI_REFE):: refe
!
!
    call elref2(nomte, 2, lielrf, ntrou)
    call elrefe_info(elrefe=lielrf(1), fami='RIGI', ndim=ndim, nno=nno1, nnos=nnos, &
                     npg=npg, jpoids=iw, jvf=ivf1, jdfde=idf1, jgano=jgn)
    call elrefe_info(elrefe=lielrf(1), fami='NOEU', ndim=ndim, nno=nno1, nnos=nnos, &
                     npg=npgn, jpoids=iwn, jvf=ivf1n, jdfde=idf1n, jgano=jgnn)
    call elrefe_info(elrefe=lielrf(2), fami='RIGI', ndim=ndim, nno=nno2, nnos=nnos, &
                     npg=npg, jpoids=iw, jvf=ivf2, jdfde=idf2, jgano=jgn)
    ndim = 3
    nddl1 = 5
!
! - DECALAGE D'INDICE POUR LES ELEMENTS D'INTERFACE
    call cginit(nomte, iu, iuc, im)
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PVECTUR', 'E', ivectu)
!
!     MISE A JOUR EVENTUELLE DE LA GEOMETRIE
!
    reactu = .false.
    if (option .eq. 'FORC_NODA') then
        call jevech('PCOMPOR', 'L', vk16=compor)
        if (compor(DEFO) .eq. 'PETIT_REAC') reactu = .true.
    end if
    if (.not. reactu) then
        do ino = 1, nno1
            do i = 1, ndim
                geom(i, ino) = zr(igeom-1+(ino-1)*ndim+i)
            end do
        end do
    else
        call jevech('PDEPLAR', 'L', jvDisp)
        do ino = 1, nno1
            do i = 1, ndim
                geom(i, ino) = zr(igeom-1+(ino-1)*ndim+i)+zr(jvDisp-1+(ino-1)*nddl1+i)
            end do
        end do
    end if
!     DEFINITION DES TANGENTES
!
    call cgtang(3, nno1, npgn, geom, zr(idf1n), &
                tang)
!
!      OPTIONS FORC_NODA ET REFE_FORC_NODA
!
    if (option .eq. 'FORC_NODA') then
!
        call jevech('PSIEFR', 'L', jvSief)
        call cgfono(ndim, nno1, nno2, npg, zr(iw), &
                    zr(ivf1), zr(ivf2), zr(idf1), geom, tang, &
                    iu, iuc, im, zr(jvSief), zr(ivectu))
!
    else
!
        call jevech('PCAGNBA', 'L', isect)
        a = zr(isect)
!
        call refe%Init(nomte)
        sigref = refe%GetRef('SIGM')
        depref = refe%GetRef('DEPL')
        forref = refe%GetRef('EFFORT')
        call refe%Check()
!
        call cgfore(ndim, nno1, nno2, npg, zr(iw), &
                    zr(ivf1), zr(ivf2), zr(idf1), a, geom, &
                    tang, iu, iuc, im, forref, &
                    sigref, depref, zr(ivectu))
!
    end if
!
end subroutine
