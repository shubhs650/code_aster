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
subroutine te0169(option, nomte)
! SUPPRESSION D'INSTRUCTIONS INUTILES
    use resi_refe_module, only: RESI_REFE
    implicit none
#include "jeveux.h"
#include "asterfort/jevech.h"
#include "blas/ddot.h"
!
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DES FORCES NODALES DE MEPOULI
!                          REFE_FORC_NODA
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
    real(kind=8) :: w(9), l1(3), l2(3), forref
    real(kind=8) :: norml1, norml2, coef1, coef2
    integer(kind=8) :: jefint, jvSief, igeom, jvDisp, ivectu, nno, nc
    integer(kind=8) :: ino, i, kc
    blas_int :: b_incx, b_incy, b_n
    type(RESI_REFE):: refe
! ----------------------------------------------------------------------
!
    if (option .eq. 'REFE_FORC_NODA') then
        nno = 3
        nc = 3
        call refe%Init(nomte)
        forref = refe%GetRef('EFFORT')
        call refe%Check()
        call jevech('PVECTUR', 'E', ivectu)
        do ino = 1, nno
            do i = 1, nc
                zr(ivectu+(ino-1)*nc+i-1) = forref
            end do
        end do
!
    else if (option .eq. 'FORC_NODA') then
        call jevech('PGEOMER', 'L', igeom)
        call jevech('PDEPLAR', 'L', jvDisp)
        call jevech('PSIEFR', 'L', jvSief)
        call jevech('PVECTUR', 'E', jefint)
!
        do i = 1, 9
            w(i) = zr(jvDisp-1+i)
        end do
!
        do kc = 1, 3
            l1(kc) = w(kc)+zr(igeom-1+kc)-w(6+kc)-zr(igeom+5+kc)
        end do
        do kc = 1, 3
            l2(kc) = w(3+kc)+zr(igeom+2+kc)-w(6+kc)-zr(igeom+5+kc)
        end do
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        norml1 = ddot(b_n, l1, b_incx, l1, b_incy)
        b_n = to_blas_int(3)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        norml2 = ddot(b_n, l2, b_incx, l2, b_incy)
        norml1 = sqrt(norml1)
        norml2 = sqrt(norml2)
!
        coef1 = zr(jvSief)/norml1
        coef2 = zr(jvSief)/norml2
!
        do i = 1, 3
            zr(jefint+i-1) = coef1*l1(i)
            zr(jefint+i+2) = coef2*l2(i)
            zr(jefint+i+5) = -zr(jefint+i-1)-zr(jefint+i+2)
        end do
    end if
end subroutine
