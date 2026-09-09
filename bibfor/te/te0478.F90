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
subroutine te0478(option, nomte)
!
    use pipeElem_module
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystNone
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/dfdm1d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/matrot.h"
#include "asterfort/pmfinfo.h"
#include "asterfort/ppga1d.h"
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
! Element: lineic
!
! Option: COOR_ELGA
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: spaceDime, nbNode, npg, jvCoorSupp, idfde, ipoids, ivf, jvGeom
    integer(kind=8) :: tab(2), iret
    integer(kind=8) :: jacf, iorien, nbsp, nbLayer, nbptcou
    integer(kind=8) :: iLayer, isp, kpg, ifi, iNode, jvCoorPg
    real(kind=8) :: copg(4, 4), copg2(3, 4), pgl(3, 3), gm1(3), gm2(3), airesp
    real(kind=8) :: layerThickness, thickness, hh, radius
    real(kind=8) :: dfdx(3), cour, jacp, cosa, sina, spoid
    aster_logical :: gauss_support
    integer(kind=8) :: nbfibr, nbgrfi, tygrfi, nbcarm, nug(10)
    integer(kind=8) :: nbFourier
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', nno=nbNode, &
                     npg=npg, jpoids=ipoids, jvf=ivf, jdfde=idfde)
    ASSERT(npg .le. 4)

! - Access to geometry
    call tecach('OOO', 'PGEOMER', 'L', iret, nval=2, itab=tab)
    spaceDime = tab(2)/nbNode
    jvGeom = tab(1)

! - Access to integration scheme
    call jevech('PCOORPG', 'E', jvCoorPg)

! - Calcul des coordonnées des points de Gauss du support si besoin
    call tecach('NNN', 'PCOORSU', 'E', iret, iad=jvCoorSupp)
    gauss_support = (iret .eq. 0)
!
! ------------------------------------------------------------------------------
!   POUTRES MULTIFIBRES
    if ((nomte .eq. 'MECA_POU_D_EM') .or. (nomte .eq. 'MECA_POU_D_TGM') .or. &
        (nomte .eq. 'MECA_POU_D_SQUE')) then
! ----- Properties of fibers
        call pmfinfo(nbfibr, nbgrfi, tygrfi, nbcarm, nug, jacf=jacf)
        ASSERT(nbgrfi .le. 10)
        call jevech('PCAORIE', 'L', iorien)
        call matrot(zr(iorien), pgl)

! ----- Points of integration scheme for SEG support
        call ppga1d(spaceDime, nbNode, npg, zr(ipoids), zr(ivf), zr(idfde), zr(jvGeom), copg)
        if (gauss_support) then
            do kpg = 1, npg
! ------------- Coordinates
                zr(jvCoorSupp+(kpg-1)*4+0) = copg(1, kpg)
                zr(jvCoorSupp+(kpg-1)*4+1) = copg(2, kpg)
                zr(jvCoorSupp+(kpg-1)*4+2) = copg(3, kpg)
! ------------- Weight
                zr(jvCoorSupp+(kpg-1)*4+3) = copg(4, kpg)
            end do
        end if

! ----- Sub-Points of integration scheme
        gm1 = 0.d0
!       boucle sur les fibres/sous-points
!           données   : nbcarm valeurs par fibre <yf,zf,Aire> + <yp,zp,Numgr>
!           résultats : 4 valeurs par fibre  <x,y,z,w>
        do ifi = 1, nbfibr
            gm1(2) = zr(jacf+(ifi-1)*nbcarm)
            gm1(3) = zr(jacf+(ifi-1)*nbcarm+1)
            call utpvlg(1, 3, pgl, gm1, gm2)
            airesp = zr(jacf+(ifi-1)*nbcarm+2)
            do kpg = 1, npg
! ------------- Coordinates
                zr(jvCoorPg+(nbfibr*(kpg-1)+(ifi-1))*4+0) = copg(1, kpg)+gm2(1)
                zr(jvCoorPg+(nbfibr*(kpg-1)+(ifi-1))*4+1) = copg(2, kpg)+gm2(2)
                zr(jvCoorPg+(nbfibr*(kpg-1)+(ifi-1))*4+2) = copg(3, kpg)+gm2(3)
! ------------- Weight
                zr(jvCoorPg+(nbfibr*(kpg-1)+(ifi-1))*4+3) = copg(4, kpg)*airesp
            end do
        end do
!
! ------------------------------------------------------------------------------
!   TUYAUX
    else if (lteatt('TUYAU', 'OUI')) then
        call pipeGetDime(nomte, 'RIGI', &
                         nbNode, nbFourier)
        call pipeCoorElga(gauss_support, &
                          nbNode, npg, nbFourier, &
                          ipoids, ivf, idfde, jvGeom, &
                          jvCoorPg, jvCoorSupp)
!
! --------------------------------------------------------------------------------------------------
!   COQUE_AXIS
    else if (nomte .eq. 'MECXSE3') then
        ASSERT(spaceDime .eq. 2)
        call getCara(plateCara, plateOrie)
        call compCoorSystNone(plateOrie)

! ----- Get layers
        nbLayer = plateCara%nbLayer
        thickness = plateCara%thick
        layerThickness = thickness/nbLayer

! ----- Points of integration scheme for SEG support
        call ppga1d(spaceDime, nbNode, npg, zr(ipoids), zr(ivf), zr(idfde), zr(jvGeom), copg2)
        if (gauss_support) then
            do kpg = 1, npg
                zr(jvCoorSupp+(kpg-1)*3+0) = copg2(1, kpg)
                zr(jvCoorSupp+(kpg-1)*3+1) = copg2(2, kpg)
                zr(jvCoorSupp+(kpg-1)*3+2) = copg2(3, kpg)
            end do
        end if

! ----- Sub-Points of integration scheme
        nbptcou = 3
        nbsp = nbptcou*nbLayer
        do kpg = 1, npg
!           Calcul du vecteur normal unitaire au point de gauss
            call dfdm1d(nbNode, zr(ipoids+kpg-1), zr(idfde+(kpg-1)*nbNode), zr(jvGeom), dfdx, &
                        cour, jacp, cosa, sina)
            radius = 0.d0
            do iNode = 1, nbNode
                radius = radius+zr(jvGeom+2*(iNode-1))*zr(ivf+(kpg-1)*nbNode+iNode-1)
            end do
            jacp = jacp*radius
            gm2(1) = cosa
            gm2(2) = sina
!
            do iLayer = 1, nbLayer
                do isp = 1, nbptcou
                    hh = -thickness/2.0+(iLayer-1+0.5d0*(isp-1))*layerThickness
                    zr(jvCoorPg+((kpg-1)*nbsp+(iLayer-1)*nbptcou+(isp-1))*3+0) = &
                        copg2(1, kpg)+hh*gm2(1)
                    zr(jvCoorPg+((kpg-1)*nbsp+(iLayer-1)*nbptcou+(isp-1))*3+1) = &
                        copg2(2, kpg)+hh*gm2(2)
                    if (isp .eq. 2) then
                        spoid = 2.0d0/3.0d0
                    else
                        spoid = 1.0d0/6.0d0
                    end if
!                   pour le poids, on multiplie par l'epaisseur par couche
                    zr(jvCoorPg+((kpg-1)*nbsp+(iLayer-1)*nbptcou+(isp-1))*3+2) = &
                        jacp*spoid*layerThickness
                end do
            end do
        end do
!
! --------------------------------------------------------------------------------------------------
!   autres éléments
    else
        call ppga1d(spaceDime, nbNode, npg, zr(ipoids), zr(ivf), zr(idfde), zr(jvGeom), &
                    zr(jvCoorPg))
    end if
!
end subroutine
