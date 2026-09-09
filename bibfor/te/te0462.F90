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
subroutine te0462(option, nomte)
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystPara, compCoorSystNone, &
                                isPlateQuad, isPlateTria
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fmater.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/plate_type.h"
#include "asterfort/tecach.h"
#include "asterfort/utpvlg.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: DKT, GRILLE_EXCENTRE
!
! Option: COOR_ELGA_MATER
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nfpgmx = 10, nbniv = 3, ndim = 3
    real(kind=8), parameter :: gm1(3) = (/0.d0, 0.d0, 1.d0/)
    integer(kind=8) :: nno, npg, ivf
    integer(kind=8) :: jvGeom, jtab(7), icopg, iret, decpo, iad
    integer(kind=8) :: nbsp, nbLayer, nfpg, decfpg
    integer(kind=8) :: ifpg, kpg, iLayer, iniv, ino
    real(kind=8) :: pgl(3, 3), xx, yy, zz
    real(kind=8) :: epais, excen, gm2(3), epc, bas, hh
    aster_logical :: grille
    character(len=8) :: fami(nfpgmx)
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!

! - Get plate parameters
    call getCara(plateCara, plateOrie)
    ASSERT(plateCara%type .eq. PLATE_DKT .or. plateCara%type .eq. PLATE_GRID)
    grille = lteatt('MODELI', 'GRC')
    if (isPlateQuad(plateCara)) then
        nno = 4
    else if (isPlateTria(plateCara)) then
        nno = 3
    else
        ASSERT(ASTER_FALSE)
    end if

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

!    ZR(ICOPG) : COORDONNEES DE SOUS-POINTS DE GAUSS
    call tecach('OOO', 'PCOOPGM', 'E', iret, nval=7, itab=jtab)
    icopg = jtab(1)
    nbsp = jtab(7)
    ASSERT(nbsp .gt. 0)

! - Get shell parameters
    nbLayer = plateCara%nbLayer
    epais = plateCara%thick
    excen = plateCara%offset
    if (grille) then
        ASSERT(nbLayer .eq. 1)
    else
        bas = -epais/2.d0+excen
        epc = epais/nbLayer
    end if

! - Calculate the transformation: global coordinate system/intrinsic coordinate system
    call compCoorSystPara(plateCara, zr(jvGeom), pgl)

! - No coordinate system from user
    call compCoorSystNone(plateOrie)

    call utpvlg(1, 3, pgl, gm1, gm2)

!
    call fmater(nfpgmx, nfpg, fami)
    decfpg = 0
    do ifpg = 1, nfpg
        call elrefe_info(fami=fami(ifpg), npg=npg, jvf=ivf)
        do kpg = 1, npg
! --------- Coordinates of Gauss point
            xx = 0.d0
            yy = 0.d0
            zz = 0.d0
            do ino = 1, nno
                xx = xx+zr(jvGeom+3*(ino-1)+0)*zr(ivf+(kpg-1)*nno+ino-1)
                yy = yy+zr(jvGeom+3*(ino-1)+1)*zr(ivf+(kpg-1)*nno+ino-1)
                zz = zz+zr(jvGeom+3*(ino-1)+2)*zr(ivf+(kpg-1)*nno+ino-1)
            end do

            if (grille) then
                decpo = ndim*(decfpg+kpg-1)
                iad = icopg+decpo
                zr(iad+0) = xx+excen*gm2(1)
                zr(iad+1) = yy+excen*gm2(2)
                zr(iad+2) = zz+excen*gm2(3)
            else
                decpo = nbLayer*nbniv*ndim*(decfpg+kpg-1)
                do iLayer = 1, nbLayer
                    do iniv = 1, nbniv
                        hh = bas+dble(iLayer-1)*epc+dble(iniv-1)*epc/2.d0
                        iad = icopg+decpo+(iLayer-1)*nbniv*ndim+(iniv-1)*ndim
                        zr(iad+0) = xx+hh*gm2(1)
                        zr(iad+1) = yy+hh*gm2(2)
                        zr(iad+2) = zz+hh*gm2(3)
                    end do
                end do
            end if
        end do
        decfpg = decfpg+npg
    end do
!
end subroutine
