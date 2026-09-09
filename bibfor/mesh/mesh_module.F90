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
! aslint: disable=W0413
!
module mesh_module
! ==================================================================================================
    implicit none
! ==================================================================================================
    private :: getSignNormalSkinToSupport
    public :: getPropertiesOfListOfCells, getPropertiesOfCell, getSkinCellSupport
    public :: getSkinCell
    public :: checkNormalOnSkinCell, checkInclude, getCellOptionForName, createNameOfCell
    public :: getNodeOptionForName, createNameOfNode, getGroupsFromCell, getMeshDimension
    public :: getFirstNodeFromNodeGroup, getListOfCellGroup, checkCellsAreSkin, checkCellsAreVolume
    public :: compMinMaxEdges
    public :: getCellsFromGroup
    public :: hasCellsDefinedFromCmd, allCellsDefinedFromCmd
! ==================================================================================================
    private
#include "asterc/r8gaem.h"
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/asmpi_comm_vect.h"
#include "asterfort/assert.h"
#include "asterfort/char8_to_int.h"
#include "asterfort/codent.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvem.h"
#include "asterfort/getvtx.h"
#include "asterfort/int_to_char8.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/lxlgut.h"
#include "asterfort/normev.h"
#include "asterfort/provec.h"
#include "asterfort/rminsp.h"
#include "asterfort/utmasu.h"
#include "asterfort/utmess.h"
#include "asterfort/utnono.h"
#include "blas/ddot.h"
#include "jeveux.h"
#include "MeshTypes_type.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! getPropertiesOfListOfCells
!
! Get properties of list of cells
!
! In  mesh             : name of mesh
! In  nbCell           : number of skin cells
! In  listCell         : list of skin cells (number)
! IO  cellNbNode       : number of nodes for each cell
! IO  cellNodeIndx     : index of first node for each cell
! Out lCellSurf        : 2D cells (surfacic) exist
! Out lCellLine        : 1D cells (lineic) exist
!
! --------------------------------------------------------------------------------------------------
    subroutine getPropertiesOfListOfCells(meshz, nbCell, listCell, cellNbNode, cellNodeIndx, &
                                          lCellSurf, lCellLine)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=*), intent(in) :: meshz
        integer(kind=8), intent(in) :: nbCell, listCell(*)
        integer(kind=8), intent(inout) :: cellNbNode(nbCell), cellNodeIndx(nbCell)
        aster_logical, intent(out) :: lCellSurf, lCellLine
! - Local
        integer(kind=8) :: iCell, cellNume, cellTypeNume
        character(len=8) :: mesh, cellTypeName
        character(len=4) :: cellTopo
        integer(kind=8), pointer :: skinCellType(:) => null()
        integer(kind=8), pointer :: meshConnexLen(:) => null()
!   ------------------------------------------------------------------------------------------------
        mesh = meshz
        lCellSurf = ASTER_FALSE
        lCellLine = ASTER_FALSE
        call jeveuo(mesh//'.TYPMAIL', 'L', vi=skinCellType)
        call jeveuo(jexatr(mesh//'.CONNEX', 'LONCUM'), 'L', vi=meshConnexLen)
        do iCell = 1, nbCell
            cellNume = listCell(iCell)
            cellNbNode(iCell) = meshConnexLen(cellNume+1)-meshConnexLen(cellNume)
            cellNodeIndx(iCell) = meshConnexLen(cellNume)
            cellTypeNume = skinCellType(cellNume)
            call jenuno(jexnum('&CATA.TM.NOMTM', cellTypeNume), cellTypeName)
            call getPropertiesOfCell(cellTypeName, cellTopo)
            if (cellTopo .eq. 'SURF') then
                lCellSurf = ASTER_TRUE
            else if (cellTopo .eq. 'LINE') then
                lCellLine = ASTER_TRUE
            else
! --------- Ignore
            end if
        end do
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getSkinCellSupport
!
! Get "volumic" cells support of skin cells
!
! In  mesh             : name of mesh
! In  nbCellSkin       : number of skin cells
! Ptr cellSkinNume     : list of skin cells (number)
! In  lCell2d          : flag if 2d cells exist
! In  lCell1d          : flag if 1D cells exist
! Ptr cellSupp         : volumic cells support of skin cells (pre-allocated before)
! In  nbCell           : number of cells where to look for support cells
! Ptr listCell_        : list of cells where to look for support cells
!
! --------------------------------------------------------------------------------------------------
    subroutine getSkinCellSupport(meshz, nbCellSkin, cellSkinNume, lCell2d, lCell1d, &
                                  cellSupp, nbCell_, listCell_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: meshz
        integer(kind=8), intent(in) :: nbCellSkin
        integer(kind=8), pointer :: cellSkinNume(:)
        aster_logical, intent(in) :: lCell2d, lCell1d
        integer(kind=8), pointer :: cellSupp(:)
        integer(kind=8), optional, intent(in) :: nbCell_
        integer(kind=8), optional, pointer :: listCell_(:)
! ----- Local
        character(len=8) :: mesh
        character(len=2) :: kdim
        real(kind=8), pointer :: meshNodeCoor(:) => null()
        integer(kind=8), parameter :: zero = 0
        aster_logical, parameter :: skinInsideVolume = ASTER_TRUE
        integer(kind=8) :: ibid(1), iSkinCell
        character(len=24), parameter :: cellSuppJv = '&&UTMASU'
        integer(kind=8), pointer :: vCellVolume(:) => null()
!   ------------------------------------------------------------------------------------------------
        mesh = meshz
        kdim = '  '
        if (lCell1d) then
            kdim = '2D'
            ASSERT(.not. lCell2d)
        end if
        if (lCell2d) then
            kdim = '3D'
            ASSERT(.not. lCell1d)
        end if
!
        call jeveuo(mesh//'.COORDO    .VALE', 'L', vr=meshNodeCoor)
        if (present(nbCell_)) then
            call utmasu(mesh, kdim, nbCellSkin, cellSkinNume, cellSuppJv, &
                        meshNodeCoor, nbCell_, listCell_, skinInsideVolume)
        else
            call utmasu(mesh, kdim, nbCellSkin, cellSkinNume, cellSuppJv, &
                        meshNodeCoor, zero, ibid, skinInsideVolume)
        end if
        call jeveuo(cellSuppJv, 'L', vi=vCellVolume)
        do iSkinCell = 1, nbCellSkin
            cellSupp(iSkinCell) = vCellVolume(iSkinCell)
        end do
        call jedetr(cellSuppJv)
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! checkNormalOnSkinCell
!
! Get "volumic" cells support of skin cells
!
! In  mesh             : name of mesh
! In  modelDime        : global dimension of model
! In  nbCellSkin       : number of skin cells
! In  cellSkinNume     : index for skin cells
! In  cellSkinNbNode   : number of nodes for skin cells
! In  cellSkinNodeIndx : index of first node for skin cells
! In  cellSupp         : index for "volumic" cells support of skin cells
!
! --------------------------------------------------------------------------------------------------
    subroutine checkNormalOnSkinCell(meshz, modelDime, nbCellSkin, cellSkinNume, cellSkinNbNode, &
                                     cellSkinNodeIndx, cellSupp, lMisoriented)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: meshz
        integer(kind=8), intent(in) :: modelDime
        integer(kind=8), intent(in) :: nbCellSkin, cellSkinNume(nbCellSkin)
        integer(kind=8), intent(in) :: cellSupp(nbCellSkin)
        integer(kind=8), intent(in) :: cellSkinNbNode(nbCellSkin), cellSkinNodeIndx(nbCellSkin)
        aster_logical, intent(out) :: lMisoriented
! ----- Local
        character(len=8) :: mesh, skinCellType, suppCellType, cellName
        integer(kind=8) :: iSkinCell, suppNume, nodeFirst, iNode
        integer(kind=8) :: skinNume
        integer(kind=8) :: skinNbNode, skinNodeNume(27), suppNbNode, suppNodeNume(27)
        real(kind=8) :: signNorm
        real(kind=8), pointer :: meshNodeCoor(:) => null()
        integer(kind=8), pointer :: meshConnexLen(:) => null()
        integer(kind=8), pointer :: meshConnex(:) => null()
        integer(kind=8), pointer :: meshCellType(:) => null()
!   ------------------------------------------------------------------------------------------------
        mesh = meshz
        call jeveuo(mesh//'.TYPMAIL        ', 'L', vi=meshCellType)
        call jeveuo(mesh//'.COORDO    .VALE', 'L', vr=meshNodeCoor)
        call jeveuo(jexatr(mesh//'.CONNEX', 'LONCUM'), 'L', vi=meshConnexLen)
        call jeveuo(mesh//'.CONNEX', 'E', vi=meshConnex)
        lMisoriented = ASTER_FALSE
!
        do iSkinCell = 1, nbCellSkin
! ----- Support cell of skin
            suppNume = cellSupp(iSkinCell)
            if (suppNume .eq. 0) then
                call utmess('A', 'FLUID1_3')
            else
! --------- Parameters of skin cell
                skinNume = cellSkinNume(iSkinCell)
                skinNbNode = cellSkinNbNode(iSkinCell)
                ASSERT(skinNbNode .le. 27)
                nodeFirst = cellSkinNodeIndx(iSKinCell)
                skinNodeNume = 0
                do iNode = 1, skinNbNode
                    skinNodeNume(iNode) = meshConnex(nodeFirst+iNode-1)
                end do
                call jenuno(jexnum('&CATA.TM.NOMTM', meshCellType(skinNume)), skinCellType)
! --------- Parameters of support cell
                nodeFirst = meshConnexLen(suppNume)
                suppNbNode = meshConnexLen(suppNume+1)-nodeFirst
                ASSERT(suppNbNode .le. 27)
                suppNodeNume = 0
                do iNode = 1, suppNbNode
                    suppNodeNume(iNode) = meshConnex(nodeFirst+iNode-1)
                end do
                call jenuno(jexnum('&CATA.TM.NOMTM', meshCellType(suppNume)), suppCellType)
                cellName = int_to_char8(skinNume)
!
! --------- Get sign of normal from skin cell to its volumic support
                call getSignNormalSkinToSupport(modelDime, skinNodeNume, suppNodeNume, &
                                                skinCellType, suppCellType, meshNodeCoor, &
                                                signNorm)
                if (signNorm .gt. 0.d0) then
                    lMisoriented = ASTER_TRUE
                end if
            end if
!
        end do
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getSignNormalSkinToSupport
!
! get sign of normal from skin cell to its volumic support
!
! In  modelDime        : global dimension of model
! In  skinNodeNume     : index of nodes in mesh for skin cell
! In  suppNodeNume     : index of nodes in mesh for support of skin cell
! In  skinCellType     : type of skin cell
! In  suppCellType     : type of support of skin cell
! In  meshNodeCoor     : coordinates of nodes
! Out signNorm         : sign of normal from skin to support
!                          >0 : normal if from skin to inside support cell
!
! --------------------------------------------------------------------------------------------------
    subroutine getSignNormalSkinToSupport(modelDime, skinNodeNume, suppNodeNume, skinCellType, &
                                          suppCellType, meshNodeCoor, signNorm)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        integer(kind=8), intent(in) :: modelDime
        integer(kind=8), intent(in) :: skinNodeNume(27), suppNodeNume(27)
        character(len=8), intent(in) :: skinCellType, suppCellType
        real(kind=8), intent(in) :: meshNodeCoor(*)
        real(kind=8), intent(out) :: signNorm
! - Local
        integer(kind=8) :: node1, node2, node3, iDime, iNode
        integer(kind=8) :: skinNbNode, suppNbNode
        real(kind=8) :: node1Coor(3), node2Coor(3), node3Coor(3)
        real(kind=8) :: n1n3(3), n1n2(3), normal(3), n1g(3), norm
        real(kind=8) :: xgm(3), xg3d(3)
        blas_int :: b_incx, b_incy, b_n
!
!   ------------------------------------------------------------------------------------------------
        signNorm = 0.d0
        normal = 0.d0
!
! - Compute normal to skin element
!
        node1 = skinNodeNume(1)
        node2 = skinNodeNume(2)
        if (modelDime .eq. 3) then
            node3 = skinNodeNume(3)
        end if
        node1Coor = 0.d0
        node2Coor = 0.d0
        do iDime = 1, modelDime
            node1Coor(iDime) = meshNodeCoor(3*(node1-1)+iDime)
            node2Coor(iDime) = meshNodeCoor(3*(node2-1)+iDime)
        end do
        n1n2 = node2Coor-node1Coor
        if (modelDime .eq. 2) then
            normal(1) = n1n2(2)
            normal(2) = -n1n2(1)
        else if (modelDime .eq. 3) then
            do iDime = 1, 3
                node3Coor(iDime) = meshNodeCoor(3*(node3-1)+iDime)
            end do
            n1n3 = node3Coor-node1Coor
            call provec(n1n2, n1n3, normal)
        else
            ASSERT(ASTER_FALSE)
        end if
        call normev(normal, norm)
!
! - Compute barycentric for skin element
!
        if (skinCellType(1:4) .eq. 'QUAD') then
            skinNbNode = 4
        else if (skinCellType(1:4) .eq. 'TRIA') then
            skinNbNode = 3
        else if (skinCellType(1:3) .eq. 'SEG') then
            skinNbNode = 2
            ASSERT(modelDime .eq. 2)
        else
            ASSERT(ASTER_FALSE)
        end if
        xgm = 0.d0
        do iNode = 1, skinNbNode
            xgm(1) = xgm(1)+meshNodeCoor(3*(skinNodeNume(iNode)-1)+1)
            xgm(2) = xgm(2)+meshNodeCoor(3*(skinNodeNume(iNode)-1)+2)
            if (modelDime .eq. 3) then
                xgm(3) = xgm(3)+meshNodeCoor(3*(skinNodeNume(iNode)-1)+3)
            end if
        end do
        xgm(1) = xgm(1)/skinNbNode
        xgm(2) = xgm(2)/skinNbNode
        xgm(3) = xgm(3)/skinNbNode
!
! - Compute barycentric for support element
!
        if (suppCellType(1:4) .eq. 'HEXA') then
            suppNbNode = 8
        else if (suppCellType(1:4) .eq. 'PENT') then
            suppNbNode = 6
        else if (suppCellType(1:4) .eq. 'PYRA') then
            suppNbNode = 5
        else if (suppCellType(1:4) .eq. 'TETR') then
            suppNbNode = 4
        else if (suppCellType(1:4) .eq. 'QUAD') then
            suppNbNode = 4
        else if (suppCellType(1:4) .eq. 'TRIA') then
            suppNbNode = 3
        else
            ASSERT(ASTER_FALSE)
        end if
        xg3d = 0.d0
        do iNode = 1, suppNbNode
            xg3d(1) = xg3d(1)+meshNodeCoor(3*(suppNodeNume(iNode)-1)+1)
            xg3d(2) = xg3d(2)+meshNodeCoor(3*(suppNodeNume(iNode)-1)+2)
            if (modelDime .eq. 3) then
                xg3d(3) = xg3d(3)+meshNodeCoor(3*(suppNodeNume(iNode)-1)+3)
            end if
        end do
        xg3d(1) = xg3d(1)/suppNbNode
        xg3d(2) = xg3d(2)/suppNbNode
        xg3d(3) = xg3d(3)/suppNbNode
!
! - Compute sign of normal from skin to support
!
        n1g = xg3d-xgm
        call normev(n1g, norm)
        b_n = to_blas_int(modelDime)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        signNorm = ddot(b_n, n1g, b_incx, normal, b_incy)
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! checkInclude
!
! check include for mesh
!
! --------------------------------------------------------------------------------------------------
    subroutine checkInclude()
!   ------------------------------------------------------------------------------------------------
! - Local
        integer(kind=8) :: cellTypeNume, nbCell
!   ------------------------------------------------------------------------------------------------
        nbCell = 0
        call jenonu(jexnom('&CATA.TM.NOMTM', 'POI1'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_POI1)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'SEG2'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_SEG2)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'SEG3'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_SEG3)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'SEG4'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_SEG4)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'TRIA3'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_TRIA3)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'TRIA6'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_TRIA6)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'TRIA7'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_TRIA7)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'TRIA10'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_TRIA10)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'QUAD4'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_QUAD4)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'QUAD8'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_QUAD8)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'QUAD9'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_QUAD9)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'QUAD12'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_QUAD12)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'TETRA4'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_TETRA4)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'TETRA10'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_TETRA10)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'TETRA15'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_TETRA15)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'TETRA20'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_TETRA20)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'PENTA6'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_PENTA6)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'PENTA15'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_PENTA15)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'PENTA18'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_PENTA18)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'PENTA21'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_PENTA21)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'PYRAM5'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_PYRAM5)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'PYRAM13'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_PYRAM13)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'PYRAM19'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_PYRAM19)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'HEXA8'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_HEXA8)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'HEXA20'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_HEXA20)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'HEXA27'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_HEXA27)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'HEXA9'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_HEXA9)
        nbCell = nbCell+1
        call jenonu(jexnom('&CATA.TM.NOMTM', 'PENTA7'), cellTypeNume)
        ASSERT(cellTypeNume .eq. MT_PENTA7)
        nbCell = nbCell+1
!
! - Total number of "physical" cells - See MeshTypes_type.h
        ASSERT(nbCell .eq. MT_NPHMAX)
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getPropertiesOfCell
!
! Get properties of cell
!
! --------------------------------------------------------------------------------------------------
    subroutine getPropertiesOfCell(cellTypeName, cellTopo_, cellOrder_)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=8), intent(in) :: cellTypeName
        character(len=4), optional, intent(out) :: cellTopo_
        integer(kind=8), optional, intent(out) :: cellOrder_
! - Local
        integer(kind=8) :: cellOrder
        character(len=4) :: cellTopo
!   ------------------------------------------------------------------------------------------------
        if (cellTypeName .eq. 'POI1') then
            cellTopo = 'PONC'
            cellOrder = 0
        else if (cellTypeName .eq. 'SEG2') then
            cellTopo = 'LINE'
            cellOrder = 1
        else if (cellTypeName .eq. 'SEG3') then
            cellTopo = 'LINE'
            cellOrder = 2
        else if (cellTypeName .eq. 'SEG4') then
            cellTopo = 'LINE'
            cellOrder = 3
        else if (cellTypeName .eq. 'TRIA3') then
            cellTopo = 'SURF'
            cellOrder = 1
        else if (cellTypeName .eq. 'TRIA6') then
            cellTopo = 'SURF'
            cellOrder = 2
        else if (cellTypeName .eq. 'TRIA7') then
            cellTopo = 'SURF'
            cellOrder = 2
        else if (cellTypeName .eq. 'TRIA10') then
            cellTopo = 'SURF'
            cellOrder = 3
        else if (cellTypeName .eq. 'QUAD4') then
            cellTopo = 'SURF'
            cellOrder = 1
        else if (cellTypeName .eq. 'QUAD8') then
            cellTopo = 'SURF'
            cellOrder = 2
        else if (cellTypeName .eq. 'QUAD9') then
            cellTopo = 'SURF'
            cellOrder = 2
        else if (cellTypeName .eq. 'QUAD12') then
            cellTopo = 'SURF'
            cellOrder = 3
        else if (cellTypeName .eq. 'TETRA4') then
            cellTopo = 'VOLU'
            cellOrder = 1
        else if (cellTypeName .eq. 'TETRA10') then
            cellTopo = 'VOLU'
            cellOrder = 2
        else if (cellTypeName .eq. 'TETRA15') then
            cellTopo = 'VOLU'
            cellOrder = 2
        else if (cellTypeName .eq. 'TETRA20') then
            cellTopo = 'VOLU'
            cellOrder = 3
        else if (cellTypeName .eq. 'PENTA6') then
            cellTopo = 'VOLU'
            cellOrder = 1
        else if (cellTypeName .eq. 'PENTA15') then
            cellTopo = 'VOLU'
            cellOrder = 2
        else if (cellTypeName .eq. 'PENTA18') then
            cellTopo = 'VOLU'
            cellOrder = 2
        else if (cellTypeName .eq. 'PENTA21') then
            cellTopo = 'VOLU'
            cellOrder = 2
        else if (cellTypeName .eq. 'PYRAM5') then
            cellTopo = 'VOLU'
            cellOrder = 1
        else if (cellTypeName .eq. 'PYRAM13') then
            cellTopo = 'VOLU'
            cellOrder = 2
        else if (cellTypeName .eq. 'PYRAM19') then
            cellTopo = 'VOLU'
            cellOrder = 2
        else if (cellTypeName .eq. 'HEXA8') then
            cellTopo = 'VOLU'
            cellOrder = 1
        else if (cellTypeName .eq. 'HEXA20') then
            cellTopo = 'VOLU'
            cellOrder = 2
        else if (cellTypeName .eq. 'HEXA27') then
            cellTopo = 'VOLU'
            cellOrder = 2
        else if (cellTypeName .eq. 'HEXA9') then
            cellTopo = 'VOLU'
            cellOrder = 1
        else if (cellTypeName .eq. 'PENTA7') then
            cellTopo = 'VOLU'
            cellOrder = 1
        else
            ASSERT(ASTER_FALSE)
        end if
        if (present(cellTopo_)) then
            cellTopo_ = cellTopo
        end if
        if (present(cellOrder_)) then
            cellOrder_ = cellOrder
        end if
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getCellOptionForName
!
! Get options from user for name of cell
!
! Out lPrefCellName    : flag if a string was given for name of cell
! Out lPrefCellNume    : flag if a integer was given for name of cell
! Out prefCellName     : string for name of cell
! Out prefCellNume     : integer for name of cell
!
! --------------------------------------------------------------------------------------------------
    subroutine getCellOptionForName(lPrefCellName, lPrefCellNume, prefCellName, &
                                    prefCellNume)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        aster_logical, intent(out) :: lPrefCellName, lPrefCellNume
        character(len=8), intent(out) :: prefCellName
        integer(kind=8), intent(out) :: prefCellNume
!   ------------------------------------------------------------------------------------------------
        lPrefCellName = ASTER_FALSE
        lPrefCellNume = ASTER_FALSE
        prefCellName = ' '
        prefCellNume = -1
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! createNameOfCell
!
! Create name for cell
!
! IO  cellName         : name of cell
! In  prefCellName     : string for name of cell
! IO  prefCellNume     : integer for name of cell
!
! --------------------------------------------------------------------------------------------------
    subroutine createNameOfCell(cellName, prefCellName, &
                                prefCellNume)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=8), intent(inout) :: cellName
        character(len=8), intent(in) :: prefCellName
        integer(kind=8), intent(inout) :: prefCellNume
! - Local
        integer(kind=8) :: lenPrefName, lenPrefPrev
        character(len=8) :: knume
!   ------------------------------------------------------------------------------------------------
        lenPrefName = lxlgut(prefCellName)
        call codent(prefCellNume, 'G', knume)
        prefCellNume = prefCellNume+1
        lenPrefPrev = lxlgut(knume)
        if (lenPrefName+lenPrefPrev .gt. 8) then
            call utmess('F', 'MESH2_1')
        end if
        cellName = knume
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getNodeOptionForName
!
! Get options from user for name of node
!
! Out lPrefNodeNume    : flag if a integer was given for name of node
! Out prefNodeName     : string for name of node
! Out prefNodeNume     : integer for name of node
!
! --------------------------------------------------------------------------------------------------
    subroutine getNodeOptionForName(lPrefNodeNume, prefNodeName, &
                                    prefNodeNume)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        aster_logical, intent(out) :: lPrefNodeNume
        character(len=8), intent(out) :: prefNodeName
        integer(kind=8), intent(out) :: prefNodeNume
!   ------------------------------------------------------------------------------------------------
        lPrefNodeNume = ASTER_FALSE
        prefNodeName = ' '
        prefNodeNume = -1
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! createNameOfNode
!
! Create name for node
!
! IO  nodeName         : name of node
! In  lPrefNodeNume    : flag if a integer was given for name of node
! In  prefNodeName     : string for name of node
! IO  prefNodeNume     : integer for name of node
!
! --------------------------------------------------------------------------------------------------
    subroutine createNameOfNode(nodeName, lPrefNodeNume, prefNodeName, &
                                prefNodeNume)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=8), intent(inout) :: nodeName
        aster_logical, intent(in) :: lPrefNodeNume
        character(len=8), intent(in) :: prefNodeName
        integer(kind=8), intent(inout) :: prefNodeNume
! - Local
        integer(kind=8) :: lenPrefName, lenPrefPrev
        character(len=8) :: knume
!   ------------------------------------------------------------------------------------------------
        lenPrefName = lxlgut(prefNodeName)
        if (lPrefNodeNume) then
            call codent(prefNodeNume, 'G', knume)
            prefNodeNume = prefNodeNume+1
            lenPrefPrev = lxlgut(knume)
            if (lenPrefName+lenPrefPrev .gt. 8) then
                call utmess('F', 'MESH2_1')
            end if
            nodeName = knume
        else
            ASSERT(ASTER_FALSE)
        end if
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getGroupsFromCell
!
! Get some groups where cells are
!
! In  mesh             : mesh
! In  cellNume         : index of current cell
! Out groupCell        : name of groups of cells
! Out nbGroupCell      : number of group
!
! --------------------------------------------------------------------------------------------------
    subroutine getGroupsFromCell(meshZ, cellNume, groupCell, nbGroupCell)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=*), intent(in) :: meshZ
        integer(kind=8), intent(in) :: cellNume
        character(len=24), intent(out) :: groupCell(4)
        integer(kind=8), intent(out) :: nbGroupCell
! - Local
        character(len=8) :: mesh
        integer(kind=8) :: iexi, nbGroup, iGroup, iCell, nbCell
        character(len=24) :: groupName
        integer(kind=8), pointer :: cellList(:) => null()
!   ------------------------------------------------------------------------------------------------
        mesh = meshZ
        nbGroupCell = 0
        groupCell = " "
        call jeexin(mesh//'.GROUPEMA', iexi)
        if (iexi .gt. 0) then
            call jelira(mesh//'.GROUPEMA', 'NMAXOC', nbGroup)
            do iGroup = 1, nbGroup
                call jeexin(jexnum(mesh//'.GROUPEMA', iGroup), iexi)
                if (iexi .eq. 0) cycle
                call jenuno(jexnum(mesh//'.GROUPEMA', iGroup), groupName)
                call jelira(jexnum(mesh//'.GROUPEMA', iGroup), 'LONUTI', nbCell)
                call jeveuo(jexnum(mesh//'.GROUPEMA', iGroup), 'L', vi=cellList)
                do iCell = 1, nbCell
                    if (cellList(iCell) .eq. cellNume) then
                        nbGroupCell = nbGroupCell+1
                        groupCell(nbGroupCell) = groupName
                        if (nbGroupCell .eq. 4) goto 100
                    end if
                end do
            end do
        end if
100     continue
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getMeshDimension
!
! Get dimension of mesh
!
! In  mesh             : mesh
! Out meshDime         : dimension of mesh
!
! --------------------------------------------------------------------------------------------------
    subroutine getMeshDimension(meshZ, meshDime)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: meshZ
        integer(kind=8), intent(out) :: meshDime
! ----- Local
        character(len=8) :: mesh
        character(len=24) :: answer
!   ------------------------------------------------------------------------------------------------
!
        mesh = meshZ
        meshDime = 0
        call dismoi('Z_CST', mesh, 'MAILLAGE', repk=answer)
        if (answer .eq. 'OUI') then
            meshDime = 2
        else
            meshDime = 3
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getFirstNodeFromNodeGroup
!
! Get first node from group of nodes
!
! In  mesh             : mesh
! In  nodeGroup        ; name of group of node
! Out nodeName         : name of first node from group of nodes
! Out nodeNume         : index of first node from group of nodes
!
! --------------------------------------------------------------------------------------------------
    subroutine getFirstNodeFromNodeGroup(meshZ, nodeGroupZ, nodeName, nodeNume)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: meshZ, nodeGroupZ
        character(len=8), intent(out) :: nodeName
        integer(kind=8), intent(out) :: nodeNume
! ----- Local
        character(len=8) :: mesh
        character(len=24) :: nodeGroup
        integer(kind=8) :: iret
!   ------------------------------------------------------------------------------------------------
!
        mesh = meshZ
        nodeGroup = nodeGroupZ
        nodeName = " "
        nodeNume = 0
        call utnono(" ", mesh, 'NOEUD', nodeGroup, nodeName, &
                    iret)
        if (iret .eq. 0) then
            nodeNume = char8_to_int(nodeName)
        else if (iret .eq. 10) then
            call utmess('F', 'MESH3_1', sk=nodeGroup)
        else if (iret .eq. 1) then
            call utmess('A', 'MESH3_2', sk=nodeGroup)
        else
            ASSERT(ASTER_FALSE)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getListOfCellGroup
!
! Get list of groups of cells
!
! In  mesh             : mesh
!
! --------------------------------------------------------------------------------------------------
    subroutine getListOfCellGroup(meshZ, factorKeywordZ, iFactorKeyword, nbGroup, listOfGroups)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: meshZ, factorKeywordZ
        integer(kind=8), intent(in) :: iFactorKeyword
        character(len=24), pointer :: listOfGroups(:)
        integer(kind=8), intent(out) :: nbGroup
! ----- Local
        character(len=8) :: mesh
        character(len=24) :: k24Dummy
        integer(kind=8) :: nbRet
!   ------------------------------------------------------------------------------------------------
!
        mesh = meshZ
        nbGroup = 0
        call getvem(mesh, 'GROUP_MA', factorKeywordZ, 'GROUP_MA', iFactorKeyword, &
                    0, k24Dummy, nbGroup)
        nbGroup = -nbGroup
        if (nbGroup .ne. 0) then
            allocate (listOfGroups(nbGroup))
            call getvem(mesh, 'GROUP_MA', factorKeywordZ, 'GROUP_MA', iFactorKeyword, &
                        nbGroup, listOfGroups, nbRet)
            ASSERT(nbRet .eq. nbGroup)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! checkCellsAreSkin
!
! Check if cells are skins' ones
!
! In  mesh             : mesh
!
! --------------------------------------------------------------------------------------------------
    subroutine checkCellsAreSkin(meshZ, nbCell, listCellNume, onlySkin1D, listCellType, &
                                 hasSkin1D, hasSkin2D)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: meshZ
        integer(kind=8), intent(in) :: nbCell
        integer(kind=8), pointer :: listCellNume(:)
        aster_logical, intent(in) :: onlySkin1D
        character(len=8), pointer :: listCellType(:)
        aster_logical, intent(out) :: hasSkin1D, hasSkin2D
! ----- Local
        character(len=8) :: mesh
        integer(kind=8) :: cellNume, iCell, cellTypeNume
        character(len=8) :: cellTypeName
        integer(kind=8), pointer :: meshTypmail(:) => null()
!
!   ------------------------------------------------------------------------------------------------
!
        mesh = meshZ
        hasSkin1D = ASTER_FALSE
        hasSkin2D = ASTER_FALSE
!
! ----- Access to mesh datastructures
        call jeveuo(mesh//'.TYPMAIL', 'L', vi=meshTypmail)
!
! ----- Check cells
        do iCell = 1, nbCell
            cellNume = listCellNume(iCell)
!
! --------- Get type of cell
            cellTypeNume = meshTypmail(cellNume)
            call jenuno(jexnum('&CATA.TM.NOMTM', cellTypeNume), cellTypeName)
            listCellType(iCell) = cellTypeName
!
! --------- Detect type
            if (cellTypeName(1:4) .eq. 'QUAD') then
                hasSkin2D = ASTER_TRUE
            else if (cellTypeName(1:4) .eq. 'TRIA') then
                hasSkin2D = ASTER_TRUE
            else if (cellTypeName(1:3) .eq. 'SEG') then
                hasSkin1D = ASTER_TRUE
            else
                call utmess('F', 'MESH3_94', sk=cellTypeName)
            end if
            if (hasSkin1D .and. hasSkin2D) then
                call utmess('F', 'MESH3_98')
            end if
        end do
!
        if (onlySkin1D .and. hasSkin2D) then
            call utmess('F', 'MESH3_92')
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! checkCellsAreVolume
!
! Check if cells are volumes' ones
!
! In  mesh             : mesh
!
! --------------------------------------------------------------------------------------------------
    subroutine checkCellsAreVolume(meshZ, nbCell, listCellNume, onlyHexaBiQuad, listCellType, &
                                   hasVolume, hasVoluNotHexaBiQ)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: meshZ
        integer(kind=8), intent(in) :: nbCell
        integer(kind=8), pointer :: listCellNume(:)
        aster_logical, intent(in) :: onlyHexaBiQuad
        character(len=8), pointer :: listCellType(:)
        aster_logical, intent(out) :: hasVolume, hasVoluNotHexaBiQ
! ----- Local
        character(len=8) :: mesh
        integer(kind=8) :: cellNume, iCell, cellTypeNume
        character(len=8) :: cellTypeName
        integer(kind=8), pointer :: meshTypmail(:) => null()
!
!   ------------------------------------------------------------------------------------------------
!
        mesh = meshZ
        hasVolume = ASTER_FALSE
        hasVoluNotHexaBiQ = ASTER_FALSE
!
! ----- Access to mesh datastructures
        call jeveuo(mesh//'.TYPMAIL', 'L', vi=meshTypmail)
!
! ----- Check cells
        do iCell = 1, nbCell
            cellNume = listCellNume(iCell)
!
! --------- Get type of cell
            cellTypeNume = meshTypmail(cellNume)
            call jenuno(jexnum('&CATA.TM.NOMTM', cellTypeNume), cellTypeName)
            listCellType(iCell) = cellTypeName
!
! --------- Detect type
            if (cellTypeName(1:4) .eq. 'HEXA') then
                hasVolume = ASTER_TRUE
                if (cellTypeName(5:6) .ne. '27') then
                    hasVoluNotHexaBiQ = ASTER_TRUE
                end if
            else if (cellTypeName(1:4) .eq. 'PENTA') then
                hasVolume = ASTER_TRUE
                hasVoluNotHexaBiQ = ASTER_TRUE
            else
                call utmess('F', 'MESH3_84', sk=cellTypeName)
            end if
        end do
!
        if (onlyHexaBiQuad .and. hasVoluNotHexaBiQ) then
            call utmess('F', 'MESH3_82')
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compMinMaxEdges
!
! Compute min and max edges on all cells in mesh (default) or only in a list of cells in mesh
!
! In  mesh             : mesh
! Out edgeMin          : length of smallest edge
! Out edgeMax          : length of greatest edge
! In  nbCell           : number of cells where to compute ()
! Ptr listCellNume     : pointer to index of cells where to compute ()
!
! --------------------------------------------------------------------------------------------------
    subroutine compMinMaxEdges(meshZ, edgeMin, edgeMax, &
                               nbCell_, listCellNume_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: meshZ
        real(kind=8), intent(out) :: edgeMin, edgeMax
        integer(kind=8), optional, intent(in) :: nbCell_
        integer(kind=8), optional, pointer :: listCellNume_(:)
! ----- Local
        character(len=8) :: mesh, cellTypeName
        integer(kind=8) :: nbCell, nbCellMesh
        integer(kind=8) :: cellNume, cellTypeNume, nodeNume
        integer(kind=8) :: iNode, iCell
        real(kind=8) :: x(MT_MAX_NBNODE_LINE), y(MT_MAX_NBNODE_LINE), z(MT_MAX_NBNODE_LINE)
        real(kind=8) :: d1, d2, d3, d4
        integer(kind=8), pointer :: meshDime(:) => null(), meshTypeCell(:) => null()
        integer(kind=8), pointer :: meshConnex(:) => null()
        real(kind=8), pointer :: meshNodeCoor(:) => null()
        aster_logical :: l_pmesh, onAllMesh
!
!   ------------------------------------------------------------------------------------------------
!
        mesh = meshZ
        edgeMin = r8gaem()
        edgeMax = 0.d0
        if (present(nbCell_)) then
            ASSERT(present(listCellNume_))
            onAllMesh = ASTER_FALSE
        else
            onAllMesh = ASTER_TRUE
        end if

! ----- Access to mesh
        call jeveuo(mesh//'.DIME', 'L', vi=meshDime)
        call jeveuo(mesh//'.TYPMAIL', 'L', vi=meshTypeCell)
        call jeveuo(mesh//'.COORDO    .VALE', 'L', vr=meshNodeCoor)
        nbCellMesh = meshDime(3)
        l_pmesh = isParallelMesh(mesh)

! ----- Set list
        if (onAllMesh) then
            nbCell = nbCellMesh
        else
            nbCell = nbCell_
            ASSERT(nbCell .le. nbCellMesh)
        end if

        do iCell = 1, nbCell
! --------- Current cell
            if (onAllMesh) then
                cellNume = iCell
            else
                cellNume = listCellNume_(iCell)
            end if

            ASSERT(cellNume .le. nbCellMesh)
            cellTypeNume = meshTypeCell(cellNume)
            call jenuno(jexnum('&CATA.TM.NOMTM', cellTypeNume), cellTypeName)
            call jeveuo(jexnum(mesh//'.CONNEX', cellNume), 'L', vi=meshConnex)

! --------- Compute max and min edges on this cell
            if (cellTypeName(1:3) .eq. 'POI') then
            else if (cellTypeName(1:3) .eq. 'SEG') then
                do iNode = 1, 2
                    nodeNume = meshConnex(iNode)
                    x(iNode) = meshNodeCoor(3*(nodeNume-1)+1)
                    y(iNode) = meshNodeCoor(3*(nodeNume-1)+2)
                    z(iNode) = meshNodeCoor(3*(nodeNume-1)+3)
                end do
                d1 = (x(2)-x(1))**2+(y(2)-y(1))**2+(z(2)-z(1))**2
                edgeMin = rminsp(edgeMin, d1, 0.d0, 0.d0, 0.d0)
                edgeMax = max(edgeMax, d1)
            else if (cellTypeName(1:4) .eq. 'TRIA') then
                do iNode = 1, 3
                    nodeNume = meshConnex(iNode)
                    x(iNode) = meshNodeCoor(3*(nodeNume-1)+1)
                    y(iNode) = meshNodeCoor(3*(nodeNume-1)+2)
                    z(iNode) = meshNodeCoor(3*(nodeNume-1)+3)
                end do
                d1 = (x(2)-x(1))**2+(y(2)-y(1))**2+(z(2)-z(1))**2
                d2 = (x(3)-x(2))**2+(y(3)-y(2))**2+(z(3)-z(2))**2
                d3 = (x(1)-x(3))**2+(y(1)-y(3))**2+(z(1)-z(3))**2
                edgeMin = rminsp(edgeMin, d1, d2, d3, 0.d0)
                edgeMax = max(edgeMax, d1, d2, d3)
            else if (cellTypeName(1:4) .eq. 'QUAD') then
                do iNode = 1, 4
                    nodeNume = meshConnex(iNode)
                    x(iNode) = meshNodeCoor(3*(nodeNume-1)+1)
                    y(iNode) = meshNodeCoor(3*(nodeNume-1)+2)
                    z(iNode) = meshNodeCoor(3*(nodeNume-1)+3)
                end do
                d1 = (x(2)-x(1))**2+(y(2)-y(1))**2+(z(2)-z(1))**2
                d2 = (x(3)-x(2))**2+(y(3)-y(2))**2+(z(3)-z(2))**2
                d3 = (x(4)-x(3))**2+(y(4)-y(3))**2+(z(4)-z(3))**2
                d4 = (x(1)-x(4))**2+(y(1)-y(4))**2+(z(1)-z(4))**2
                edgeMin = rminsp(edgeMin, d1, d2, d3, d4)
                edgeMax = max(edgeMax, d1, d2, d3, d4)
            else if (cellTypeName(1:4) .eq. 'HEXA') then
                do iNode = 1, 8
                    nodeNume = meshConnex(iNode)
                    x(iNode) = meshNodeCoor(3*(nodeNume-1)+1)
                    y(iNode) = meshNodeCoor(3*(nodeNume-1)+2)
                    z(iNode) = meshNodeCoor(3*(nodeNume-1)+3)
                end do
                d1 = (x(2)-x(1))**2+(y(2)-y(1))**2+(z(2)-z(1))**2
                d2 = (x(3)-x(2))**2+(y(3)-y(2))**2+(z(3)-z(2))**2
                d3 = (x(4)-x(3))**2+(y(4)-y(3))**2+(z(4)-z(3))**2
                d4 = (x(1)-x(4))**2+(y(1)-y(4))**2+(z(1)-z(4))**2
                edgeMin = rminsp(edgeMin, d1, d2, d3, d4)
                edgeMax = max(edgeMax, d1, d2, d3, d4)
                d1 = (x(6)-x(5))**2+(y(6)-y(5))**2+(z(6)-z(5))**2
                d2 = (x(7)-x(6))**2+(y(7)-y(6))**2+(z(7)-z(6))**2
                d3 = (x(8)-x(7))**2+(y(8)-y(7))**2+(z(8)-z(7))**2
                d4 = (x(5)-x(8))**2+(y(5)-y(8))**2+(z(5)-z(8))**2
                edgeMin = rminsp(edgeMin, d1, d2, d3, d4)
                edgeMax = max(edgeMax, d1, d2, d3, d4)
                d1 = (x(1)-x(5))**2+(y(1)-y(5))**2+(z(1)-z(5))**2
                d2 = (x(2)-x(6))**2+(y(2)-y(6))**2+(z(2)-z(6))**2
                d3 = (x(3)-x(7))**2+(y(3)-y(7))**2+(z(3)-z(7))**2
                d4 = (x(4)-x(8))**2+(y(4)-y(8))**2+(z(4)-z(8))**2
                edgeMin = rminsp(edgeMin, d1, d2, d3, d4)
                edgeMax = max(edgeMax, d1, d2, d3, d4)
            else if (cellTypeName(1:5) .eq. 'PENTA') then
                do iNode = 1, 6
                    nodeNume = meshConnex(iNode)
                    x(iNode) = meshNodeCoor(3*(nodeNume-1)+1)
                    y(iNode) = meshNodeCoor(3*(nodeNume-1)+2)
                    z(iNode) = meshNodeCoor(3*(nodeNume-1)+3)
                end do
                d1 = (x(2)-x(1))**2+(y(2)-y(1))**2+(z(2)-z(1))**2
                d2 = (x(3)-x(2))**2+(y(3)-y(2))**2+(z(3)-z(2))**2
                d3 = (x(1)-x(3))**2+(y(1)-y(3))**2+(z(1)-z(3))**2
                edgeMin = rminsp(edgeMin, d1, d2, d3, 0.d0)
                edgeMax = max(edgeMax, d1, d2, d3)
                d1 = (x(5)-x(4))**2+(y(5)-y(4))**2+(z(5)-z(4))**2
                d2 = (x(6)-x(5))**2+(y(6)-y(5))**2+(z(6)-z(5))**2
                d3 = (x(4)-x(6))**2+(y(4)-y(6))**2+(z(4)-z(6))**2
                edgeMin = rminsp(edgeMin, d1, d2, d3, 0.d0)
                edgeMax = max(edgeMax, d1, d2, d3)
                d1 = (x(1)-x(4))**2+(y(1)-y(4))**2+(z(1)-z(4))**2
                d2 = (x(2)-x(5))**2+(y(2)-y(5))**2+(z(2)-z(5))**2
                d3 = (x(3)-x(6))**2+(y(3)-y(6))**2+(z(3)-z(6))**2
                edgeMin = rminsp(edgeMin, d1, d2, d3, 0.d0)
                edgeMax = max(edgeMax, d1, d2, d3)
            else if (cellTypeName(1:5) .eq. 'TETRA') then
                do iNode = 1, 4
                    nodeNume = meshConnex(iNode)
                    x(iNode) = meshNodeCoor(3*(nodeNume-1)+1)
                    y(iNode) = meshNodeCoor(3*(nodeNume-1)+2)
                    z(iNode) = meshNodeCoor(3*(nodeNume-1)+3)
                end do
                d1 = (x(2)-x(1))**2+(y(2)-y(1))**2+(z(2)-z(1))**2
                d2 = (x(3)-x(2))**2+(y(3)-y(2))**2+(z(3)-z(2))**2
                d3 = (x(1)-x(3))**2+(y(1)-y(3))**2+(z(1)-z(3))**2
                edgeMin = rminsp(edgeMin, d1, d2, d3, 0.d0)
                edgeMax = max(edgeMax, d1, d2, d3)
                d1 = (x(4)-x(1))**2+(y(4)-y(1))**2+(z(4)-z(1))**2
                d2 = (x(4)-x(2))**2+(y(4)-y(2))**2+(z(4)-z(2))**2
                d3 = (x(4)-x(3))**2+(y(4)-y(3))**2+(z(4)-z(3))**2
                edgeMin = rminsp(edgeMin, d1, d2, d3, 0.d0)
                edgeMax = max(edgeMax, d1, d2, d3)
            else if (cellTypeName(1:4) .eq. 'PYRA') then
                do iNode = 1, 5
                    nodeNume = meshConnex(iNode)
                    x(iNode) = meshNodeCoor(3*(nodeNume-1)+1)
                    y(iNode) = meshNodeCoor(3*(nodeNume-1)+2)
                    z(iNode) = meshNodeCoor(3*(nodeNume-1)+3)
                end do
                d1 = (x(2)-x(1))**2+(y(2)-y(1))**2+(z(2)-z(1))**2
                d2 = (x(3)-x(2))**2+(y(3)-y(2))**2+(z(3)-z(2))**2
                d3 = (x(4)-x(3))**2+(y(4)-y(3))**2+(z(4)-z(3))**2
                d4 = (x(1)-x(4))**2+(y(1)-y(4))**2+(z(1)-z(4))**2
                edgeMin = rminsp(edgeMin, d1, d2, d3, d4)
                edgeMax = max(edgeMax, d1, d2, d3, d4)
                d1 = (x(5)-x(1))**2+(y(5)-y(1))**2+(z(5)-z(1))**2
                d2 = (x(5)-x(2))**2+(y(5)-y(2))**2+(z(5)-z(2))**2
                d3 = (x(5)-x(3))**2+(y(5)-y(3))**2+(z(5)-z(3))**2
                d4 = (x(5)-x(4))**2+(y(5)-y(4))**2+(z(5)-z(4))**2
                edgeMin = rminsp(edgeMin, d1, d2, d3, d4)
                edgeMax = max(edgeMax, d1, d2, d3, d4)
            else
                call utmess('F', 'MESH3_10', sk=cellTypeName)
            end if

        end do
        if (l_pmesh) then
            call asmpi_comm_vect("MPI_MIN", 'R', scr=edgeMin)
            call asmpi_comm_vect("MPI_MAX", 'R', scr=edgeMax)
        end if
        if (edgeMin .eq. r8gaem()) then
            edgeMin = 0.d0
        end if
        edgeMin = sqrt(edgeMin)
        edgeMax = sqrt(edgeMax)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getCellsFromGroup
!
! Get list of cells from group
!
! In  mesh             : mesh
! In  groupName        : name of group of cell
! Out edgeMax          : length of greatest edge
! Out nbCell           : number of cells
! Ptr listCellNume     : pointer to index of cells
!
! --------------------------------------------------------------------------------------------------
    subroutine getCellsFromGroup(meshZ, groupNameZ, nbCell, listCellNume)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: meshZ, groupNameZ
        integer(kind=8), intent(out) :: nbCell
        integer(kind=8), pointer :: listCellNume(:)
! ----- Local
        character(len=8) :: mesh
        character(len=24) :: groupName
        integer(kind=8) :: iexi, iret
!
!   ------------------------------------------------------------------------------------------------
!
        mesh = meshZ
        groupName = groupNameZ

        call jeexin(mesh//'.GROUPEMA', iexi)
        if (iexi .gt. 0) then
            call jeexin(jexnom(mesh//'.GROUPEMA', groupName), iret)
            if (iret .eq. 0) then
                call utmess('F', 'MESH3_11', sk=groupName)
            else
                call jelira(jexnom(mesh//'.GROUPEMA', groupName), 'LONMAX', nbCell)
                if (nbCell .ne. 0) then
                    call jeveuo(jexnom(mesh//'.GROUPEMA', groupName), 'L', vi=listCellNume)
                end if
            end if
        else
            call utmess('F', 'MESH3_12')
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! hasCellsDefinedFromCmd
!
! Has cells defined from command file ?
!
! In  factorKeyword    : factor keyword to read
! In  iOcc             : index of factor keyword
!
! --------------------------------------------------------------------------------------------------
    function hasCellsDefinedFromCmd(factorKeywordZ, iOcc)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        aster_logical :: hasCellsDefinedFromCmd
        character(len=*), intent(in) :: factorKeywordZ
        integer(kind=8), intent(in) :: iOcc
! ----- Local
        character(len=24) :: factorKeyword
        integer(kind=8) :: iret
!   ------------------------------------------------------------------------------------------------
!
        hasCellsDefinedFromCmd = ASTER_FALSE
        factorKeyword = factorKeywordZ
        call getvtx(factorKeyword, 'GROUP_MA', iocc=iOcc, nbval=0, nbret=iret)
        if (iret .ne. 0) then
            hasCellsDefinedFromCmd = ASTER_TRUE
        end if
!
!   ------------------------------------------------------------------------------------------------
    end function
! --------------------------------------------------------------------------------------------------
!
! allCellsDefinedFromCmd
!
! All cells are selected in command file ?
!
! In  factorKeyword    : factor keyword to read
! In  iOcc             : index of factor keyword
!
! --------------------------------------------------------------------------------------------------
    function allCellsDefinedFromCmd(factorKeywordZ, iOcc)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        aster_logical :: allCellsDefinedFromCmd
        character(len=*), intent(in) :: factorKeywordZ
        integer(kind=8), intent(in) :: iOcc
! ----- Local
        character(len=24) :: factorKeyword
        character(len=8) :: answer
        integer(kind=8) :: iret
!   ------------------------------------------------------------------------------------------------
!
        allCellsDefinedFromCmd = ASTER_FALSE
        factorKeyword = factorKeywordZ

        if (factorKeyword .eq. " ") then
            call getvtx(' ', 'TOUT', nbval=0, nbret=iret)
            if (iret .ne. 0) then
                allCellsDefinedFromCmd = ASTER_TRUE
            end if
        else
            call getvtx(factorKeyword, 'TOUT', iocc=iOcc, nbval=0, nbret=iret)
            if (iret .ne. 0) then
                iret = abs(iret)
                ASSERT(iret .eq. 1)
                call getvtx(factorKeyword, 'TOUT', iocc=iOcc, scal=answer)
                ASSERT(answer .eq. 'OUI')
                allCellsDefinedFromCmd = ASTER_TRUE
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end function
! --------------------------------------------------------------------------------------------------
!
! getSkinCell
!
! Get all skin cells
!
! In  mesh             : name of mesh
! In  cellSkinDime     : dimension of skin cell
! In  nbCell           : number of cells
! Ptr listCell         : list of cells
! Out nbCellSkin       : number of skin cells
! Ptr listCellSkin     : list of skin cells
!
! --------------------------------------------------------------------------------------------------
    subroutine getSkinCell(meshZ, cellSkinDime, &
                           nbCell, listCell, &
                           nbCellSkin, listCellSkin)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: meshZ
        integer(kind=8), intent(in) :: nbCell, cellSkinDime
        integer(kind=8), pointer :: listCell(:)
        integer(kind=8), intent(out) :: nbCellSkin
        integer(kind=8), pointer :: listCellSkin(:)
! ----- Local
        character(len=8) :: mesh
        integer(kind=8) :: iCellType, nbCellType, iCell
        integer(kind=8) :: nbWorkCell, iWorkCell
        character(len=8), pointer :: cataTypeCell(:) => null()
        integer(kind=8), pointer :: cataDimeCell(:) => null()
        integer(kind=8), pointer :: meshTypmail(:) => null()
        integer(kind=8), pointer :: workCell(:) => null()
        integer(kind=8), pointer :: listCellType(:) => null()
!   ------------------------------------------------------------------------------------------------
!
        mesh = meshZ
        ASSERT(nbCell .gt. 0)
        ASSERT(cellSkinDime .ge. 1)

! ----- Access to mesh
        call jeveuo(mesh(1:8)//'.TYPMAIL        ', 'L', vi=meshTypmail)

! ----- Access to catalog of list of type cell
        call jelira('&CATA.TM.NOMTM', 'NOMMAX', nbCellType)

! ----- Working vectors
        AS_ALLOCATE(vk8=cataTypeCell, size=nbCellType)
        AS_ALLOCATE(vi=cataDimeCell, size=nbCellType)
        AS_ALLOCATE(vi=listCellType, size=nbCell)

! ----- Get all type cells in catalog
        do iCellType = 1, nbCellType
            call jenuno(jexnum('&CATA.TM.NOMTM', iCellType), cataTypeCell(iCellType))
            call dismoi('DIM_TOPO', cataTypeCell(iCellType), 'TYPE_MAILLE', &
                        repi=cataDimeCell(iCellType))
        end do

! ----- Get all skin cells
        AS_ALLOCATE(vi=workCell, size=nbCell)
        nbWorkCell = 0
        iWorkCell = 1
        do iCell = 1, nbCell
            listCellType(iCell) = meshTypmail(listCell(iCell))
            if (cataDimeCell(listCellType(iCell)) .eq. cellSkinDime) then
                workCell(iWorkCell) = listCell(iCell)
                nbWorkCell = nbWorkCell+1
                iWorkCell = iWorkCell+1
            end if
        end do

! ----- Final copy
        if (nbWorkCell .ne. 0) then
            nbCellSkin = nbWorkCell
            AS_ALLOCATE(vi=listCellSkin, size=nbCellSkin)
            do iWorkCell = 1, nbWorkCell
                listCellSkin(iWorkCell) = workCell(iWorkCell)
            end do
        end if
!
        AS_DEALLOCATE(vk8=cataTypeCell)
        AS_DEALLOCATE(vi=cataDimeCell)
        AS_DEALLOCATE(vi=workCell)
        AS_DEALLOCATE(vi=listCellType)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module mesh_module
