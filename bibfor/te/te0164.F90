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
subroutine te0164(option, nomte)
    use resi_refe_module, only: RESI_REFE
    implicit none
#include "jeveux.h"
#include "asterfort/biline.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/matvec.h"
!
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DES FORCES NODALES DE MECABL2
!                          OPTION : 'FORC_NODA', 'REFE_FORC_NODA'
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
    real(kind=8) :: coef, jacobi, nx, ytywpq(9), w(9), forref
    integer(kind=8) :: nno, kp, i, ipoids, ivf, igeom, nc, nordre, k
    integer(kind=8) :: ivectu, ino, ndim, nnos, npg
    integer(kind=8) :: idfdk, jgano, iyty, jvDisp, jvSief, jefint
    type(RESI_REFE):: refe
! ----------------------------------------------------------------------
!
    if (option .eq. 'REFE_FORC_NODA') then
        nno = 2
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
        call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                         jpoids=ipoids, jvf=ivf, jdfde=idfdk, jgano=jgano)
        call jevete('&INEL.CABPOU.YTY', 'L', iyty)
        nordre = 3*nno
        call jevech('PGEOMER', 'L', igeom)
        call jevech('PDEPLAR', 'L', jvDisp)
        call jevech('PSIEFR', 'L', jvSief)
        call jevech('PVECTUR', 'E', jefint)
!

        do i = 1, 3*nno
            w(i) = zr(jvDisp-1+i)
        end do

        do kp = 1, npg
            k = (kp-1)*nordre*nordre
            jacobi = sqrt(biline(nordre, zr(igeom), zr(iyty+k), zr(igeom)))
            nx = zr(jvSief-1+kp)
            call matvec(nordre, zr(iyty+k), 2, zr(igeom), w, &
                        ytywpq)
            coef = nx*zr(ipoids-1+kp)/jacobi
            do i = 1, nordre
                zr(jefint-1+i) = zr(jefint-1+i)+coef*ytywpq(i)
            end do
        end do
    end if
end subroutine
