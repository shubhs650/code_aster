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
! ==================================================================================================
!
! Module for geometry of plates
!
! ==================================================================================================
!
module plateGeom_module
! ==================================================================================================
    use plate_type
    use calcul_module, only: ca_nomte_, ca_nomtm_
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: checkPlaneity, getCara, compCoorSystPara
    public :: compCoorSystPlate, compCoorSystCO3D, updateCoorSystCO3D, compCoorSystNone
    public :: compCoorSystGrid, compCoorSystMemb
    public :: getManifoldBase, creaCaraMini
    public :: isPlateTria, isPlateQuad, isShell3D, isPlate, isPlateDKT, isPlateDKTG, isPlateQ4GG
    private :: getType, getNbLayer
! ==================================================================================================
    private
#include "asterc/r8dgrd.h"
#include "asterc/r8miem.h"
#include "asterc/r8vide.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/coqrep.h"
#include "asterfort/dxqpgl.h"
#include "asterfort/dxtpgl.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/plate_type.h"
#include "asterfort/plateGeom_module.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/vdrepe.h"
#include "asterfort/vdxrep.h"
#include "asterfort/vectgt.h"
#include "jeveux.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! checkPlaneity
!
! Check planeity of quadrangular cell
!
! In  xyzg             : coordinates of vertices
! In  errorTole        : tolerance for check
! Out errorCode        : error return code
! Out distAbso         : distance to plane (absolute)
! Out distRela         : distance to plane (relative)
!
! --------------------------------------------------------------------------------------------------
    subroutine checkPlaneity(xyzg, errorTole, errorCode, distAbso, distRela)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        real(kind=8), intent(in) :: xyzg(3, 4)
        integer(kind=8), intent(out) :: errorCode
        real(kind=8), intent(in) :: errorTole
        real(kind=8), intent(out) :: distAbso, distRela
! ----- Local
        real(kind=8) :: x12, y12, z12, x13, y13, z13, x14, y14, z14
        real(kind=8) :: ux, uy, uz, pscal, normu, norm4, dist
!   ------------------------------------------------------------------------------------------------
!
        errorCode = BASE_NO_ERROR
        distAbso = 0.d0
        distRela = 0.d0

! ----- First vector
        x12 = xyzg(1, 2)-xyzg(1, 1)
        y12 = xyzg(2, 2)-xyzg(2, 1)
        z12 = xyzg(3, 2)-xyzg(3, 1)

! ----- Second vector
        x13 = xyzg(1, 3)-xyzg(1, 1)
        y13 = xyzg(2, 3)-xyzg(2, 1)
        z13 = xyzg(3, 3)-xyzg(3, 1)

! ----- Third vector
        x14 = xyzg(1, 4)-xyzg(1, 1)
        y14 = xyzg(2, 4)-xyzg(2, 1)
        z14 = xyzg(3, 4)-xyzg(3, 1)

! ----- Define plane on node 1-2-3
        ux = (y12*z13)-(y13*z12)
        uy = (z12*x13)-(z13*x12)
        uz = (x12*y13)-(x13*y12)

        pscal = (ux*x14)+(uy*y14)+(uz*z14)

! ----- Compute normal to plane
        normu = sqrt((ux*ux)+(uy*uy)+(uz*uz))
        if (normu .lt. r8miem()) then
! --------- Something wrong: degenerated normal
            errorCode = BASE_CELL_DEGE
        else
            norm4 = sqrt((x14*x14)+(y14*y14)+(z14*z14))
            dist = pscal/normu
            pscal = dist/norm4
            if (abs(pscal) .gt. errorTole) then
                errorCode = BASE_QUAD_NOPLANE
                distAbso = abs(dist)
                distRela = pscal*100
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getCara
!
! Get caracteristics of plates (thickness, shear coefficient, etc.)
!
! --------------------------------------------------------------------------------------------------
    subroutine getCara(plateCara, plateOrie)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(plateCara_Para), intent(inout) :: plateCara
        type(plateOrie_Para), intent(inout) :: plateOrie
! ----- Local
        integer(kind=8) :: jvCacoqu, iret, jTab(8)
        integer(kind=8), parameter :: indxC3D(10) = (/1, 4, 0, 5, 6, 7, 2, 3, 0, 0/)
        integer(kind=8), parameter :: indxDkt(10) = (/1, 0, 0, 4, 5, 6, 2, 3, 0, 0/)
        integer(kind=8), parameter :: indxMITC(10) = (/1, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        integer(kind=8), parameter :: indxAxi(10) = (/1, 2, 3, 0, 0, 0, 0, 0, 0, 0/)
        integer(kind=8), parameter :: indxGri(10) = (/0, 0, 0, 5, 4, 0, 2, 3, 1, 0/)
        integer(kind=8), parameter :: indxMem(10) = (/1, 0, 0, 0, 0, 0, 2, 3, 0, 4/)
        integer(kind=8), parameter :: indxTh1(10) = (/1, 0, 0, 0, 0, 0, 0, 0, 0, 0/)
        integer(kind=8), parameter :: indxTh2(10) = (/1, 0, 0, 0, 0, 0, 2, 3, 0, 0/)
        integer(kind=8) :: indxSelect(10)
!   ------------------------------------------------------------------------------------------------
!
        call tecach('NNO', 'PCACOQU', 'L', iret, nval=8, itab=jTab)

! ----- Detect type
        call getType(plateCara)

! ----- Number of layers
        call getNbLayer(plateCara)

! ----- Get parameters
        if (iret .eq. 0) then
            jvCacoqu = jTab(1)
            indxSelect = 0
            if (plateCara%type .eq. PLATE_CO3D) then
                indxSelect = indxC3D
            elseif (plateCara%type .eq. PLATE_DKT) then
                indxSelect = indxDkt
            elseif (plateCara%type .eq. PLATE_DKTG) then
                indxSelect = indxDkt
            elseif (plateCara%type .eq. PLATE_MITC) then
                indxSelect = indxMITC
            elseif (plateCara%type .eq. PLATE_COAX) then
                indxSelect = indxAxi
            elseif (plateCara%type .eq. PLATE_GRID) then
                indxSelect = indxGri
            elseif (plateCara%type .eq. PLATE_MEMB) then
                indxSelect = indxMem
            elseif (plateCara%type .eq. PLATE_COTH1) then
                indxSelect = indxTh1
            elseif (plateCara%type .eq. PLATE_COTH2) then
                indxSelect = indxTh2
            elseif (plateCara%type .eq. PLATE_DST) then
                indxSelect = indxDkt
            elseif (plateCara%type .eq. PLATE_Q4G) then
                indxSelect = indxDkt
            elseif (plateCara%type .eq. PLATE_Q4GG) then
                indxSelect = indxDkt
            elseif (plateCara%type .eq. PLATE_SOLID_LINK) then
                indxSelect = indxDkt
            elseif (plateCara%type .eq. PLATE_BOUND) then
                indxSelect = indxDkt
            else
                ASSERT(ASTER_FALSE)
            end if

            if (indxSelect(1) .ne. 0) plateCara%thick = zr(jvCacoqu-1+indxSelect(1))
            if (indxSelect(2) .ne. 0) plateCara%shearCoef = zr(jvCacoqu-1+indxSelect(2))
            if (indxSelect(3) .ne. 0) plateCara%metric = zr(jvCacoqu-1+indxSelect(3))
            if (indxSelect(4) .ne. 0) plateCara%coefRigiDRZ = zr(jvCacoqu-1+indxSelect(4))
            if (indxSelect(5) .ne. 0) plateCara%offset = zr(jvCacoqu-1+indxSelect(5))
            if (indxSelect(6) .ne. 0) plateCara%inerRota = zr(jvCacoqu-1+indxSelect(6))
            if (indxSelect(9) .ne. 0) plateCara%section = zr(jvCacoqu-1+indxSelect(9))
            if (indxSelect(10) .ne. 0) plateCara%tension = zr(jvCacoqu-1+indxSelect(10))
            if (indxSelect(7) .ne. 0) then
                plateOrie%alpha = zr(jvCacoqu-1+indxSelect(7))*r8dgrd()
                ASSERT(indxSelect(8) .ne. 0)
                plateOrie%beta = zr(jvCacoqu-1+indxSelect(8))*r8dgrd()
                plateOrie%lRead = ASTER_TRUE
            end if
            plateCara%lRead = ASTER_TRUE
        else
            plateCara%lRead = ASTER_FALSE
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compCoorSystPlate
!
! Compute coordinate system for plate
!
! --------------------------------------------------------------------------------------------------
    subroutine compCoorSystPlate(pgl, plateCara, plateOrie)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        real(kind=8), intent(in) :: pgl(3, 3)
        type(plateCara_Para), intent(inout) :: plateCara
        type(plateOrie_Para), intent(inout) :: plateOrie
! ----- Local
        real(kind=8) :: alpha, beta, t1ve(9)
        real(kind=8) :: t2iu(4), t2ui(4), c, s
        character(len=16) :: elasKeyword
        integer(kind=8) :: jvMaterc, iret
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(plateCara%lRead)
        ASSERT(plateOrie%lRead)
        ASSERT(isPlate(plateCara))
        alpha = plateOrie%alpha
        beta = plateOrie%beta

! ----- Compute operators for coordinate transformation
        call coqrep(pgl, alpha, beta, &
                    t2iu, t2ui, c, s)
        plateOrie%t2iu = t2iu
        plateOrie%t2ui = t2ui
        plateOrie%c = c
        plateOrie%s = s

        call tecach('NNO', 'PMATERC', 'L', iret, iad=jvMaterc)
        elasKeyword = "ELAS"
        if (iret .eq. 0) then
            call rccoma(zi(jvMaterc), 'ELAS', 1, elasKeyword)
        end if
        t1ve = 0.d0
        if (elasKeyword .ne. "ELAS") then
            t1ve(1) = plateOrie%c*plateOrie%c
            t1ve(4) = plateOrie%s*plateOrie%s
            t1ve(7) = plateOrie%c*plateOrie%s
            t1ve(2) = t1ve(4)
            t1ve(5) = t1ve(1)
            t1ve(8) = -t1ve(7)
            t1ve(3) = -t1ve(7)-t1ve(7)
            t1ve(6) = t1ve(7)+t1ve(7)
            t1ve(9) = t1ve(1)-t1ve(4)
        end if
        plateOrie%t1ve = t1ve
        plateOrie%lUpdate = ASTER_TRUE
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getType
!
! Get coordinate system of plate
!
! --------------------------------------------------------------------------------------------------
    subroutine getType(plateCara)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(plateCara_Para), intent(inout) :: plateCara
! ----- Local
        integer(kind=8) :: plateType
!   ------------------------------------------------------------------------------------------------
!
        plateType = PLATE_UNKW
        if (lteatt('MODELI', 'CQ3')) then
            plateType = PLATE_CO3D
        elseif (lteatt('MODELI', 'DKT')) then
            plateType = PLATE_DKT
        elseif (lteatt('MODELI', 'MIT')) then
            plateType = PLATE_MITC
        elseif (lteatt('MODELI', 'DTG')) then
            plateType = PLATE_DKTG
        elseif (lteatt('MODELI', 'DST')) then
            plateType = PLATE_DST
        elseif (lteatt('MODELI', 'Q4G')) then
            plateType = PLATE_Q4G
        elseif (lteatt('MODELI', 'Q4S')) then
            plateType = PLATE_Q4GG
        elseif (lteatt('MODELI', 'RC3')) then
            plateType = PLATE_SOLID_LINK
        elseif (lteatt('MODELI', 'CQA')) then
            plateType = PLATE_COAX
        elseif (lteatt('MODELI', 'GRC') .or. lteatt('MODELI', 'GRM')) then
            plateType = PLATE_GRID
        elseif (lteatt('MODELI', 'MMB')) then
            plateType = PLATE_MEMB
        elseif (lteatt('MODELI', 'CTA') .or. lteatt('MODELI', 'CQP')) then
            plateType = PLATE_COTH1
        elseif (lteatt('MODELI', 'CQ_')) then
            plateType = PLATE_COTH2
        elseif (ca_nomte_ .eq. "MEBODKT" .or. &
                ca_nomte_ .eq. "MEBOQ4G" .or. &
                ca_nomte_ .eq. "MEBODST") then
            plateType = PLATE_BOUND
        else
            ASSERT(ASTER_FALSE)
        end if
        plateCara%type = plateType
        plateCara%geom = PLATE_GEOM_UNKW
        if (ca_nomtm_ .eq. 'TRIA3' .or. ca_nomtm_ .eq. 'TRIA6' .or. &
            ca_nomtm_ .eq. 'TRIA7') then
            plateCara%geom = PLATE_GEOM_TRIA
        else if (ca_nomtm_ .eq. 'QUAD4' .or. ca_nomtm_ .eq. 'QUAD8' .or. &
                 ca_nomtm_ .eq. 'QUAD9') then
            plateCara%geom = PLATE_GEOM_QUAD
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getNbLayer
!
! Get number of layers
!
! --------------------------------------------------------------------------------------------------
    subroutine getNbLayer(plateCara)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(plateCara_Para), intent(inout) :: plateCara
! ----- Local
        integer(kind=8) :: jvNbsp, nbLayer, iret
!   ------------------------------------------------------------------------------------------------
!
        nbLayer = 0
        ASSERT(plateCara%type .ne. PLATE_UNKW)
        call tecach('NNO', 'PNBSP_I', 'L', iret, iad=jvNbsp)
        if (iret .eq. 0) then
            nbLayer = zi(jvNbsp-1+1)
        else
            nbLayer = 1
        end if
        if (plateCara%type .eq. PLATE_DKTG) then
            ASSERT(nbLayer .eq. 1)
        end if
        if (nbLayer .le. 0) then
            call utmess('F', 'PLATE1_10')
        end if
        plateCara%nbLayer = nbLayer
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compCoorSystCO3D
!
! Compute coordinate system for COQUE_3D
!
! --------------------------------------------------------------------------------------------------
    subroutine compCoorSystCO3D(nomte, jvGeom, &
                                plateCara, plateOrie)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=16), intent(in) :: nomte
        integer(kind=8), intent(in) :: jvGeom
        type(plateCara_Para), intent(in) :: plateCara
        type(plateOrie_Para), intent(inout) :: plateOrie
! ----- Local
        real(kind=8) :: thick
        integer(kind=8), parameter :: npgt = 10
        real(kind=8) :: matevn(2, 2, npgt), matevg(2, 2, npgt)
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(plateCara%lRead)
        ASSERT(plateCara%type .eq. PLATE_CO3D)
        ASSERT(plateOrie%lRead)
        thick = plateCara%thick

! ----- Compute local basis on cell (at nodes)
        call vdxrep(plateOrie, &
                    nomte, thick, zr(jvGeom))

! ----- Compute global<=>local matrices (at nodes and integration points)
        call vdrepe(plateOrie, nomte, matevn, matevg)
        plateOrie%matevn = matevn
        plateOrie%matevg = matevg
        plateOrie%lUpdate = ASTER_TRUE
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! updateCoorSystCO3D
!
! Update coordinate system for COQUE_3D
!
! --------------------------------------------------------------------------------------------------
    subroutine updateCoorSystCO3D(plateCara, plateOrie, &
                                  nomte, nodeCoor, &
                                  npgsr, nb1, lzr, &
                                  matevn, matevg)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(plateCara_Para), intent(in) :: plateCara
        type(plateOrie_Para), intent(in) :: plateOrie
        character(len=16), intent(in) :: nomte
        real(kind=8), intent(in) :: nodeCoor(3, 9)
        integer(kind=8), intent(in) :: npgsr, nb1, lzr
        real(kind=8), intent(out) :: matevn(2, 2, npgsr), matevg(2, 2, npgsr)
! ----- Local
        integer(kind=8), parameter :: ptTypeRedu = 0
        real(kind=8), parameter :: ptMiddleCoor = 0.d0
        integer(kind=8) :: kpgsr, i, j, k
        real(kind=8) :: thick
        real(kind=8) :: vectBaseKpg(3, 3)
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(plateCara%lRead)
        ASSERT(plateOrie%lRead)
        ASSERT(plateOrie%lUpdate)
        thick = plateCara%thick
        k = 0
        do kpgsr = 1, npgsr
            call vectgt(plateOrie, ptTypeRedu, nb1, &
                        nodeCoor, ptMiddleCoor, kpgsr, &
                        thick, zr(lzr), &
                        vectBaseKpg)
            do j = 1, 3
                do i = 1, 3
                    k = k+1
                    zr(lzr+2000+k-1) = vectBaseKpg(i, j)
                end do
            end do
        end do

! ----- Compute global<=>local matrices (at nodes and integration points)
        call vdrepe(plateOrie, nomte, matevn, matevg)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compCoorSystNone
!
! No coordinate system !
!
! --------------------------------------------------------------------------------------------------
    subroutine compCoorSystNone(plateOrie)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(plateOrie_Para), intent(inout) :: plateOrie
!   ------------------------------------------------------------------------------------------------
!
        plateOrie%lUpdate = ASTER_FALSE
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getManifoldBase
!
! --------------------------------------------------------------------------------------------------
    subroutine getManifoldBase(nb1, nb2, desr, &
                               nodeCoor, &
                               vecta)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(out) :: nb1, nb2
        real(kind=8), intent(inout) :: desr(*)
        real(kind=8), intent(in) :: nodeCoor(3, 9)
        real(kind=8), intent(out) :: vecta(9, 2, 3)
! ----- Locals
        integer(kind=8), parameter :: jvDFuncKsi = 828, jvDFuncEta = 900
        integer(kind=8) :: i, j, k
!   ------------------------------------------------------------------------------------------------
!
        vecta = 0.d0
        do i = 1, nb2
            do k = 1, 3
                vecta(i, 1, k) = 0.d0
                vecta(i, 2, k) = 0.d0
                do j = 1, nb1
                    vecta(i, 1, k) = vecta(i, 1, k)+ &
                                     desr(jvDFuncKsi+8*(i-1)+j)*nodeCoor(k, j)
                    vecta(i, 2, k) = vecta(i, 2, k)+ &
                                     desr(jvDFuncEta+8*(i-1)+j)*nodeCoor(k, j)
                end do
            end do
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compCoorSystGrid
!
! Compute coordinate system for grids
!
! --------------------------------------------------------------------------------------------------
    subroutine compCoorSystGrid(pgl, plateCara, plateOrie)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        real(kind=8), intent(in) :: pgl(3, 3)
        type(plateCara_Para), intent(in) :: plateCara
        type(plateOrie_Para), intent(inout) :: plateOrie
! ----- Local
        real(kind=8) :: alpha, beta
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(plateCara%lRead)
        ASSERT(plateOrie%lRead)
        ASSERT(plateCara%type .eq. PLATE_GRID)
        alpha = plateOrie%alpha
        beta = plateOrie%beta
        plateOrie%gridDir11(1) = cos(beta)*cos(alpha)
        plateOrie%gridDir11(2) = cos(beta)*sin(alpha)
        plateOrie%gridDir11(3) = -sin(beta)
        if (lteatt('MODELI', 'GRC')) then
            plateOrie%gridNorm(1:3) = plateCara%offset*pgl(3, 1:3)
        else
            plateOrie%gridNorm = r8vide()
        end if
        plateOrie%lUpdate = ASTER_TRUE
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compCoorSystMemb
!
! Compute coordinate system for membrane
!
! --------------------------------------------------------------------------------------------------
    subroutine compCoorSystMemb(plateOrie)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(plateOrie_Para), intent(inout) :: plateOrie
!   ------------------------------------------------------------------------------------------------
!
        plateOrie%lUpdate = ASTER_TRUE
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! isShell3D
!
! --------------------------------------------------------------------------------------------------
    function isShell3D(plateCara)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        aster_logical :: isShell3D
        type(plateCara_Para), intent(in) :: plateCara
! --------------------------------------------------------------------------------------------------
!
        isShell3D = (plateCara%type .eq. PLATE_CO3D)
!
!   ------------------------------------------------------------------------------------------------
    end function
! --------------------------------------------------------------------------------------------------
!
! isPlateDKT
!
! --------------------------------------------------------------------------------------------------
    function isPlateDKT(plateCara)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        aster_logical :: isPlateDKT
        type(plateCara_Para), intent(in) :: plateCara
! --------------------------------------------------------------------------------------------------
!
        isPlateDKT = (plateCara%type .eq. PLATE_DKT)
!
!   ------------------------------------------------------------------------------------------------
    end function
! --------------------------------------------------------------------------------------------------
!
! isPlateDKTG
!
! --------------------------------------------------------------------------------------------------
    function isPlateDKTG(plateCara)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        aster_logical :: isPlateDKTG
        type(plateCara_Para), intent(in) :: plateCara
! --------------------------------------------------------------------------------------------------
!
        isPlateDKTG = (plateCara%type .eq. PLATE_DKTG)
!
!   ------------------------------------------------------------------------------------------------
    end function
! --------------------------------------------------------------------------------------------------
!
! isPlateQ4GG
!
! --------------------------------------------------------------------------------------------------
    function isPlateQ4GG(plateCara)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        aster_logical :: isPlateQ4GG
        type(plateCara_Para), intent(in) :: plateCara
! --------------------------------------------------------------------------------------------------
!
        isPlateQ4GG = (plateCara%type .eq. PLATE_Q4GG)
!
!   ------------------------------------------------------------------------------------------------
    end function
! --------------------------------------------------------------------------------------------------
!
! isPlate
!
! --------------------------------------------------------------------------------------------------
    function isPlate(plateCara)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        aster_logical :: isPlate
        type(plateCara_Para), intent(in) :: plateCara
! --------------------------------------------------------------------------------------------------
!
        isPlate = (plateCara%type .eq. PLATE_DKT .or. plateCara%type .eq. PLATE_DKTG .or. &
                   plateCara%type .eq. PLATE_Q4G .or. plateCara%type .eq. PLATE_Q4GG .or. &
                   plateCara%type .eq. PLATE_DST)
!
!   ------------------------------------------------------------------------------------------------
    end function
! --------------------------------------------------------------------------------------------------
!
! isPlateTria
!
! --------------------------------------------------------------------------------------------------
    function isPlateTria(plateCara)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        aster_logical :: isPlateTria
        type(plateCara_Para), intent(in) :: plateCara
! --------------------------------------------------------------------------------------------------
!
        isPlateTria = plateCara%geom .eq. PLATE_GEOM_TRIA
!
!   ------------------------------------------------------------------------------------------------
    end function
! --------------------------------------------------------------------------------------------------
!
! isPlateQuad
!
! --------------------------------------------------------------------------------------------------
    function isPlateQuad(plateCara)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        aster_logical :: isPlateQuad
        type(plateCara_Para), intent(in) :: plateCara
! --------------------------------------------------------------------------------------------------
!
        isPlateQuad = plateCara%geom .eq. PLATE_GEOM_QUAD
!
!   ------------------------------------------------------------------------------------------------
    end function
! --------------------------------------------------------------------------------------------------
!
! compCoorSystPara
!
! Calculate the coordinate transformation between the global coordinate system and the
! plate's intrinsic coordinate system.
!
! --------------------------------------------------------------------------------------------------
    subroutine compCoorSystPara(plateCara, xyzg, pgl)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(plateCara_Para), intent(in) :: plateCara
        real(kind=8), intent(in) :: xyzg(3, *)
        real(kind=8), intent(out) :: pgl(3, 3)
!   ------------------------------------------------------------------------------------------------
!
        pgl = 0.d0
        if (isPlateTria(plateCara)) then
            call dxtpgl(xyzg, pgl)
        elseif (isPlateQuad(plateCara)) then
            call dxqpgl(xyzg, pgl)
        else
            ASSERT(ASTER_FALSE)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! creaCaraMini
!
! --------------------------------------------------------------------------------------------------
    subroutine creaCaraMini(plateCara)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(plateCara_Para), intent(inout) :: plateCara
!   ------------------------------------------------------------------------------------------------
!
        call getType(plateCara)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module plateGeom_module
