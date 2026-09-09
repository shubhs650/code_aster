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
subroutine te0401(option, nomte)
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystCO3D
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/bsthco.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/matpgl.h"
#include "asterfort/r8inir.h"
#include "asterfort/tranlg.h"
#include "asterfort/utvtsv.h"
#include "asterfort/vdxrig.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
!     CALCUL DES OPTIONS DES ELEMENTS DE COQUE : COQUE_3D
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: indm, indf
    integer(kind=8) :: nb1, nb2, nddlet
    integer(kind=8) :: lzr
    integer(kind=8) :: jvGeom, jener
    integer(kind=8) :: i, j, kompt
    integer(kind=8) :: jvDisp, imatuu
    real(kind=8) :: matrRigiLoca(51, 51), plg(9, 3, 3)
    real(kind=8) :: vrs(1326)
    real(kind=8) :: bsigth(51), enerth
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!
    enerth = 0.0d0

! - Get plate parameters
    call getCara(plateCara, plateOrie)

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Compute global<=>local transformation
    call compCoorSystCO3D(nomte, jvGeom, &
                          plateCara, plateOrie)

! - Access to static objects of COQUE_3D
    call jevete('&INEL.'//nomte(1:8)//'.DESR', ' ', lzr)

! - Compute elastic matrix
    if (option .eq. 'RIGI_MECA' .or. option .eq. 'EPOT_ELEM') then
        indm = 0
        indf = 0
        call vdxrig(plateCara, plateOrie, &
                    nomte, zr(jvGeom), matrRigiLoca, nb1, &
                    indm, indf)

! ----- Get matrix
        nb2 = nb1+1
        call matpgl(nb2, zr(lzr), plg)

! ----- Frame modification of matrix
        nddlet = 6*nb1+3
        vrs = 0.d0
        call tranlg(nb1, 51, nddlet, plg, matrRigiLoca, vrs)
!
    else
        ASSERT(ASTER_FALSE)
    end if
!
!
    if (option .eq. 'RIGI_MECA') then
        call jevech('PMATUUR', 'E', imatuu)
        kompt = 0
        do j = 1, 6*nb1+3
            do i = 1, j
                kompt = kompt+1
                zr(imatuu-1+kompt) = vrs(kompt)
            end do
        end do
!
    end if

!---- ENERGIES DE DEFORMATION ELASTIQUE
    if (option .eq. 'EPOT_ELEM') then
        call jevech('PENERDR', 'E', jener)
        call jevech('PDEPLAR', 'L', jvDisp)

!------ ENERGIE DE DEFORMATION TOTALE
        call utvtsv('ZERO', 6*nb1+3, vrs, zr(jvDisp), zr(jener))
!
        zr(jener) = 0.5d0*zr(jener)

!
        call bsthco(plateCara, plateOrie, &
                    nomte, bsigth)
!
        do i = 1, 6*nb1+3
            enerth = enerth+bsigth(i)*zr(jvDisp+i-1)
        end do
        zr(jener) = zr(jener)-enerth
!
        if (abs(zr(jener)) .gt. 1.d-6) then
!--------- ENERGIE DE DEFORMATION DE MEMBRANE
            indm = 1
            indf = 0
            call vdxrig(plateCara, plateOrie, &
                        nomte, zr(jvGeom), matrRigiLoca, nb1, &
                        indm, indf)

            vrs = 0.d0
            call tranlg(nb1, 51, nddlet, plg, matrRigiLoca, vrs)
            call utvtsv('ZERO', 6*nb1+3, vrs, zr(jvDisp), zr(jener+1))
            zr(jener+1) = 0.5d0*zr(jener+1)

!--------- ENERGIE DE DEFORMATION DE FLEXION
            indm = 0
            indf = 1
            call vdxrig(plateCara, plateOrie, &
                        nomte, zr(jvGeom), matrRigiLoca, nb1, &
                        indm, indf)
            call r8inir(1326, 0.d0, vrs, 1)
            vrs = 0.d0
            call tranlg(nb1, 51, nddlet, plg, matrRigiLoca, vrs)
            call utvtsv('ZERO', 6*nb1+3, vrs, zr(jvDisp), zr(jener+2))
!
            zr(jener+2) = 0.5d0*zr(jener+2)
            zr(jener+1) = zr(jener+1)/zr(jener)
            zr(jener+2) = zr(jener+2)/zr(jener)
!
        else
            call r8inir(2, 0.d0, zr(jener+1), 1)
        end if
    end if
!
end subroutine
