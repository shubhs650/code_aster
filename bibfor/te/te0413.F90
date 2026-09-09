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
subroutine te0413(option, nomte)
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystPara, compCoorSystPlate, &
                                isPlateTria, isPlateQuad
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/glrc_recup_mate.h"
#include "asterfort/gquad4.h"
#include "asterfort/gtria3.h"
#include "asterfort/jevech.h"
#include "asterfort/jquad4.h"
#include "asterfort/utmess.h"
#include "asterfort/utpvgl.h"
#include "jeveux.h"
!
    character(len=16) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: DKTG
!
! Options: DISS_ELEM, DISS_ELGA
!
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: npgmx = 4
    real(kind=8) :: pgl(3, 3)
    real(kind=8) :: qsi, eta, xyzl(3, 4), jacob(5), poids, cara(25)
    real(kind=8) :: disse(npgmx), dse
    real(kind=8) :: ep, seuil
    integer(kind=8) :: ndim, nno, npg, ipoids, icoopg
    integer(kind=8) :: jvGeom, kpg, jvDiss, jvMaterc
    integer(kind=8) :: jvVari, nbvar
    character(len=16), pointer :: compor(:) => null()
    character(len=16) :: valk(2)
    aster_logical :: lkit
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, npg=npg, &
                     jpoids=ipoids, jcoopg=icoopg)
    ASSERT(npg .le. npgmx)

! - Get plate parameters
    call getCara(plateCara, plateOrie)

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Calculate the transformation: global coordinate system/intrinsic coordinate system
    call compCoorSystPara(plateCara, zr(jvGeom), pgl)

! - Compute coordinate system for plate
    call compCoorSystPlate(pgl, plateCara, plateOrie)

! - Change coordinates of geometry
    call utpvgl(nno, 3, pgl, zr(jvGeom), xyzl)

! - Compute geometric parametrs of plate
    if (isPlateQuad(plateCara)) then
        call gquad4(xyzl, cara)
    elseif (isPlateTria(plateCara)) then
        call gtria3(xyzl, cara)
    else
        ASSERT(ASTER_FALSE)
    end if

! - Get plate parameters
    ep = plateCara%thick

! - Non-linear behaviour
    call jevech('PCOMPOR', 'L', vk16=compor)
    lkit = compor(RELA_NAME) (1:7) .eq. 'KIT_DDI'
    read (compor(NVAR), '(I16)') nbvar

! - Get internal state variables
    if (option .eq. 'DISS_ELGA') then
        call jevech('PVARIGR', 'L', jvVari)
    else if (option .eq. 'DISS_ELEM') then
        call jevech('PVARIPR', 'L', jvVari)
    end if

! - Get output field
    if (option .eq. 'DISS_ELGA') then
        call jevech('PDISSPG', 'E', jvDiss)
    else if (option .eq. 'DISS_ELEM') then
        call jevech('PDISSD1', 'E', jvDiss)
    end if
!
    if ((compor(RELA_NAME) (1:7) .eq. 'GLRC_DM') .or. &
        (lkit .and. (compor(CREEP_NAME) (1:7) .eq. 'GLRC_DM'))) then

        disse = 0.d0
        dse = 0.0d0
        do kpg = 1, npg
            qsi = zr(icoopg-1+ndim*(kpg-1)+1)
            eta = zr(icoopg-1+ndim*(kpg-1)+2)
            if (isPlateQuad(plateCara)) then
                call jquad4(xyzl, qsi, eta, jacob)
                poids = zr(ipoids+kpg-1)*jacob(1)
            else
                poids = zr(ipoids+kpg-1)*cara(7)
            end if
            call jevech('PMATERC', 'L', jvMaterc)
            call glrc_recup_mate(zi(jvMaterc), compor(RELA_NAME), .false._1, ep, seuil=seuil)
            if ((option .eq. 'DISS_ELGA') .or. (option .eq. 'DISS_ELEM')) then
                disse(kpg) = (zr(jvVari-1+(kpg-1)*nbvar+1)+ &
                              zr(jvVari-1+(kpg-1)*nbvar+2))*seuil
                dse = dse+disse(kpg)*poids
            end if
        end do

        if (option .eq. 'DISS_ELGA') then
            do kpg = 1, npg
                zr(jvDiss-1+(kpg-1)*1+1) = disse(kpg)
            end do
        else if (option .eq. 'DISS_ELEM') then
            zr(jvDiss-1+1) = dse
        end if

    elseif (compor(RELA_NAME) (1:4) .eq. 'DHRC') then
        disse = 0.D0
        dse = 0.0d0

        do kpg = 1, npg
            qsi = zr(icoopg-1+ndim*(kpg-1)+1)
            eta = zr(icoopg-1+ndim*(kpg-1)+2)
            if (isPlateQuad(plateCara)) then
                call jquad4(xyzl, qsi, eta, jacob)
                poids = zr(ipoids+kpg-1)*jacob(1)
            else
                poids = zr(ipoids+kpg-1)*cara(7)
            end if

            if ((option .eq. 'DISS_ELGA') .or. (option .eq. 'DISS_ELEM')) then
                disse(kpg) = zr(jvVari-1+(kpg-1)*nbvar+9)
                dse = dse+disse(kpg)*poids
            end if
        end do

        if (option .eq. 'DISS_ELGA') then
            do kpg = 1, npg
                zr(jvDiss-1+(kpg-1)*1+1) = disse(kpg)
            end do
        else if (option .eq. 'DISS_ELEM') then
            zr(jvDiss-1+1) = dse
        end if

    else
        valk(1) = option
        valk(2) = compor(RELA_NAME) (1:7)
        call utmess('A', 'ELEMENTS4_63', nk=2, valk=valk)
    end if
!
end subroutine
