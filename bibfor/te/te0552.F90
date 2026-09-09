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
subroutine te0552(option, nomte)
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystPara, compCoorSystGrid, compCoorSystPlate
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/plate_type.h"
#include "asterfort/utmess.h"
#include "asterfort/utpvgl.h"
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
! Options: DEPL_ELGA
!
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter  :: mxnoeu = 4, mxnpg = 4, nbcompo = 6, nbdepl = 3
    integer(kind=8), parameter  :: sp_couche_dkt = 3, sp_couche_gri = 1
    integer(kind=8) :: icmp, indga, indno, indx, ino, ipg, iptc, iLayer
    integer(kind=8) :: nno, npg, nbLayer
    integer(kind=8) :: jvDispElga, jvf
    integer(kind=8) :: jvDispNode, jvGeom
    aster_logical :: elem_dkt, elem_gri
    real(kind=8) :: pgl(3, 3)
    real(kind=8) :: deplno(nbcompo*mxnoeu), deplga(nbcompo*mxnpg), deplsp(nbdepl)
    real(kind=8) :: epcoqu, hcouche, zic, zmin, excent
    character(len=4), parameter :: fami = 'RIGI'
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(option .eq. 'DEPL_ELGA')
    call elrefe_info(fami=fami, nno=nno, npg=npg, jvf=jvf)
    ASSERT(nno .le. mxnoeu)
    ASSERT(npg .le. mxnpg)

! - Get plate parameters
    call getCara(plateCara, plateOrie)
    elem_dkt = plateCara%type .eq. PLATE_DKT
    elem_gri = plateCara%type .eq. PLATE_GRID

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Calculate the transformation: global coordinate system/intrinsic coordinate system
    call compCoorSystPara(plateCara, zr(jvGeom), pgl)

    if (elem_dkt) then
        call compCoorSystPlate(pgl, plateCara, plateOrie)
    elseif (elem_gri) then
        call compCoorSystGrid(pgl, plateCara, plateOrie)
    else
        ASSERT(ASTER_FALSE)
    end if

! - Displacements (at nodes)
    call jevech('PDEPLAR', 'L', jvDispNode)

! - Passage des déplacements dans le repère local
    call utpvgl(nno, nbcompo, pgl, zr(jvDispNode), deplno)

! - Output field: displacements ELGA
    call jevech('PDEPLGA', 'E', jvDispElga)

! - Calcul des déplacements aux points de gauss de l'élément support
    deplga = 0.0
    do ipg = 1, npg
        do ino = 1, nno
            do icmp = 1, nbcompo
                indga = (ipg-1)*nbcompo+icmp
                indno = (ino-1)*nbcompo+icmp
                deplga(indga) = deplga(indga)+deplno(indno)*zr(jvf+(ipg-1)*nno+ino-1)
            end do
        end do
    end do

! - Calcul des déplacements aux sous-points
    if (elem_dkt) then
        epcoqu = plateCara%thick
        excent = plateCara%offset
        nbLayer = plateCara%nbLayer
        if (nbLayer .le. 0) then
            call utmess('F', 'PLATE1_10')
        end if
        hcouche = epcoqu/nbLayer
        zmin = excent-epcoqu*0.50

        do iLayer = 1, nbLayer
            do iptc = 1, sp_couche_dkt
                if (iptc .eq. 1) then
                    zic = zmin+(iLayer-1)*hcouche
                else if (iptc .eq. 2) then
                    zic = zmin+(iLayer-1)*hcouche+hcouche*0.50
                else
                    zic = zmin+(iLayer-1)*hcouche+hcouche
                end if
                do ipg = 1, npg
                    indga = (ipg-1)*nbcompo
                    deplsp(1) = deplga(indga+1)+deplga(indga+5)*zic
                    deplsp(2) = deplga(indga+2)-deplga(indga+4)*zic
                    deplsp(3) = deplga(indga+3)
                    indx = jvDispElga+nbdepl*(sp_couche_dkt*(nbLayer*(ipg-1)+iLayer-1)+iptc-1)
                    call utpvlg(1, nbdepl, pgl, deplsp, zr(indx))
                end do
            end do
        end do
    else if (elem_gri) then
        epcoqu = plateCara%thick
        nbLayer = plateCara%nbLayer
        excent = plateCara%offset
        ASSERT(nbLayer .eq. 1)
        iLayer = 1
        iptc = 1
        zic = excent

        do ipg = 1, npg
            indga = (ipg-1)*nbcompo
            deplsp(1) = deplga(indga+1)+deplga(indga+5)*zic
            deplsp(2) = deplga(indga+2)-deplga(indga+4)*zic
            deplsp(3) = deplga(indga+3)
            indx = jvDispElga+nbdepl*(sp_couche_gri*(nbLayer*(ipg-1)+iLayer-1)+iptc-1)
            ! Passage des déplacements dans le repère global
            call utpvlg(1, nbdepl, pgl, deplsp, zr(indx))
        end do
    end if
!
end subroutine
