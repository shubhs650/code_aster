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
subroutine te0395(option, nomte)
    use resi_refe_module, only: RESI_REFE
    implicit none
#include "jeveux.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/nmasf3.h"
#include "asterfort/r8inir.h"
#include "asterfort/tecach.h"
#include "blas/daxpy.h"
!
    character(len=16) :: option, nomte
! ----------------------------------------------------------------------
! FONCTION REALISEE:  CALCUL DE L'OPTION FORC_NODA ELEMENT HEXAS8
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
    real(kind=8) :: bsigm(3, 8), geo(24), sigtmp(6), ftemp(24), sigref
    integer(kind=8) :: jgano, nno, k, npg1, i, j, ivectu, ndim, nnos
    integer(kind=8) :: ipoids, ivf, idfde, igeom, jvSief, imate, jvDisp
    integer(kind=8) :: icomp, ii, iretc, iretd
    blas_int :: b_incx, b_incy, b_n
    type(RESI_REFE):: refe
! DEB ------------------------------------------------------------------
!
! ---- CARACTERISTIQUES DU TYPE D'ELEMENT :
! ---- GEOMETRIE ET INTEGRATION
!      ------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg1, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
! ---- PARAMETRES EN ENTREE
! ----     COORDONNEES DES CONNECTIVITES
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imate)
    do i = 1, ndim*nno
        geo(i) = zr(igeom-1+i)
    end do
!
! ---- PARAMETRES EN SORTIE
!      --------------------
! ----     VECTEUR DES FORCES INTERNES (BT*SIGMA)
    call jevech('PVECTUR', 'E', ivectu)
!
! ---- CALCUL DE FORC_NODA
    call tecach('ONO', 'PCOMPOR', 'L', iretc, iad=icomp)
!
    if (option .eq. 'FORC_NODA') then
!      --------------------
!         CHAMPS POUR LA REACTUALISATION DE LA GEOMETRIE
        call jevech('PDEPLAR', 'L', jvDisp)
!
! ----     CONTRAINTES AUX POINTS D'INTEGRATION
        call jevech('PSIEFR', 'L', jvSief)
!
! ---- CALCUL DU VECTEUR DES FORCES INTERNES (BT*SIGMA) :
!      --------------------------------------------------
        call nmasf3(nno, npg1, ipoids, ivf, idfde, &
                    zi(imate), geo, zr(jvDisp), zr(jvSief), zr(ivectu), &
                    zk16(icomp))
!
!
    else if (option .eq. 'REFE_FORC_NODA') then
        call refe%Init(nomte)
        sigref = refe%GetRef('SIGM')
        call refe%Check()
!
        call tecach('ONO', 'PDEPLMR', 'L', iretd, iad=jvDisp)
!
        call r8inir(6*npg1, 0.d0, sigtmp, 1)
        call r8inir(3*nno, 0.d0, ftemp, 1)
        do i = 1, 6*npg1
!
            sigtmp(i) = sigref
            call nmasf3(nno, npg1, ipoids, ivf, idfde, &
                        zi(imate), geo, zr(jvDisp), sigtmp, bsigm, &
                        zk16(icomp))
!
            do j = 1, nno
                ii = 3*(j-1)
                do k = 1, 3
                    ftemp(ii+k) = ftemp(ii+k)+abs(bsigm(k, j))
                end do
            end do
!
        end do
!
        b_n = to_blas_int(ndim*nno)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0/npg1, ftemp, b_incx, zr(ivectu), &
                   b_incy)
!
    end if
!
! FIN ------------------------------------------------------------------
end subroutine
