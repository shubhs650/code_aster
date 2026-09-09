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
subroutine te0008(option, nomte)
    use resi_refe_module, only: RESI_REFE
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/bsigmc.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/nbsigm.h"
#include "asterfort/tecach.h"
#include "MeshTypes_type.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "FE_module.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! ELEMENTS ISOPARAMETRIQUES 2D ET 3D
!
! CALCUL DES OPTIONS FORC_NODA ET REFE_FORC_NODA
!
! --------------------------------------------------------------------------------------------------
!
! IN  OPTION : OPTION DE CALCUL
! IN  NOMTE  : NOM DU TYPE ELEMENT
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: sigref, sigtmp(6*MT_NNOMAX)
    real(kind=8) :: nharm, bsigm(MAX_BV_CG), geo(3*MT_NNOMAX), ftemp(MAX_BV_CG)
    integer(kind=8) :: nbsig, ndim, nno, npg
    integer(kind=8) :: ipoids, ivf, idfde
    integer(kind=8) :: igeom, ivectu
    integer(kind=8) :: jvDisp, jvSief, jvCompor
    integer(kind=8) :: i, j, nbinco, iretCompor, iretDisp
    blas_int :: b_incx, b_incy, b_n
    type(RESI_REFE):: refe
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, npg=npg, jpoids=ipoids, &
                     jvf=ivf, jdfde=idfde)
!
! - Initializations
    nharm = 0.d0
    nbinco = nno*ndim
    ASSERT(nbinco .le. MAX_BV_CG)
    nbsig = nbsigm()
    ASSERT(nbsig .le. 6)
!
! - Output vector
    call jevech('PVECTUR', 'E', ivectu)
!
! - Get input field
    call jevech('PGEOMER', 'L', igeom)
    b_n = to_blas_int(nbinco)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, zr(igeom), b_incx, geo, b_incy)
!
! - Compute
    if (option .eq. 'FORC_NODA') then
        call jevech('PSIEFR', 'L', jvSief)
        call tecach('ONO', 'PCOMPOR', 'L', iretCompor, iad=jvCompor)
        if (iretCompor .eq. 0) then
            if (zk16(jvCompor+2) (1:6) .ne. 'PETIT ') then
                call tecach('ONO', 'PDEPLAR', 'L', iretDisp, iad=jvDisp)
                if (iretDisp .eq. 0) then
                    call jevech('PDEPLAR', 'L', jvDisp)
                    b_n = to_blas_int(nbinco)
                    b_incx = to_blas_int(1)
                    b_incy = to_blas_int(1)
                    call daxpy(b_n, 1.d0, zr(jvDisp), b_incx, geo, &
                               b_incy)
                end if
            end if
        end if
        call bsigmc(nno, ndim, nbsig, npg, ipoids, &
                    ivf, idfde, geo, nharm, zr(jvSief), &
                    bsigm)
        b_n = to_blas_int(nbinco)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, bsigm, b_incx, zr(ivectu), b_incy)
!
    else if (option .eq. 'REFE_FORC_NODA') then
        call refe%Init(nomte)
        sigref = refe%GetRef('SIGM')
        call refe%Check()
        sigtmp = 0.d0
        ftemp = 0.d0
        do i = 1, nbsig*npg
            sigtmp(i) = sigref
            call bsigmc(nno, ndim, nbsig, npg, ipoids, &
                        ivf, idfde, geo, nharm, sigtmp, &
                        bsigm)
            do j = 1, nbinco
                ftemp(j) = ftemp(j)+abs(bsigm(j))
            end do
            sigtmp(i) = 0.d0
        end do
        b_n = to_blas_int(nbinco)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0/npg, ftemp, b_incx, zr(ivectu), &
                   b_incy)
    end if
!
end subroutine
