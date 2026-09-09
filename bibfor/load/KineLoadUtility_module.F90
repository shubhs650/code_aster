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
! Module for kinematic loads
!
! ==================================================================================================
!
module KineLoadUtility_module
! ==================================================================================================
    use KineListRela_type
    use KineListRela_module
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: kineLoadGetSlaveNodes, kineLoadGetSlaveCells, kineLoadGetSlaveNodesFromSkin
    public :: kineLoadGetMasterNodes, kineLoadGetMasterCells
    public :: kineLoadNormFromSkinCells, kineLoadNormAtNodes
    public :: kineLoadSetLCS, kineLoadElimMult, kineLoadReadTransf, kineLoadApplyTransf
    public :: kineLoadGetPhysQuanInfo, kineLoadDeleteNodePair
    public :: kineLoadCheckCmpOnNode, kineLoadApplyEccentricity
! ==================================================================================================
    private
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/indik8.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/assert.h"
#include "asterfort/calirg.h"
#include "asterfort/canorm.h"
#include "asterfort/canort.h"
#include "asterfort/char_read_tran.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisdg.h"
#include "asterfort/fointe.h"
#include "asterfort/getelem.h"
#include "asterfort/getnode.h"
#include "asterfort/irmail.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jelira.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/nbnlma.h"
#include "asterfort/normev.h"
#include "asterfort/pacoor.h"
#include "asterfort/reliem.h"
#include "asterfort/ulaffe.h"
#include "asterfort/ulnume.h"
#include "asterfort/ulopen.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
#include "MeshTypes_type.h"
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! kineLoadGetSlaveNodesFromSkin
!
! Get list of slave nodes from skin
!
! In  model            : model
! In  mesh             : mesh
! In  meshNbNode       : total number of nodes of mesh
! In  geomDime         : geometric dimension (2 or 3)
! In  factorKeyword    : factor keyword
! In  iOcc             : index of factor keyword
! Out nodeSlavJv       : name of JEVEUX object for list of slave nodes
! Out nbNodeSlav       : number of slave nodes
! Ptr nodeSlav         : pointer to list of slave nodes
! Out nbNodeSlav       : number of slave cells
! Ptr cellSlav         : pointer to list of slave cells
! Ptr nodeSkinToBody   : pointer between body and skin nodes
!
! --------------------------------------------------------------------------------------------------
    subroutine kineLoadGetSlaveNodesFromSkin(modelZ, meshZ, meshNbNode, geomDime, &
                                             factorKeywordZ, iOcc, &
                                             nodeSlavJv, nbNodeSlav, nodeSlav, &
                                             nbCellSlav, cellSlav, &
                                             nodeSkinToBody)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=*), intent(in) :: modelZ, meshZ, factorKeywordZ
        integer(kind=8), intent(in) :: iOcc, meshNbNode, geomDime
        character(len=24), intent(out) :: nodeSlavJv
        integer(kind=8), intent(out) :: nbNodeSlav, nbCellSlav
        integer(kind=8), pointer :: nodeSlav(:), cellSlav(:)
        integer(kind=8), pointer :: nodeSkinToBody(:)
! - Local
        character(len=24) :: cellSlavJv
        aster_logical :: lError
        character(len=8) :: cellNameError, mesh, model
        integer(kind=8) :: iNodeSlav
        integer(kind=8), parameter :: nbCellSkin2D = 3, nbCellSkin3D = 8
        character(len=8), parameter :: cellSkin2D(nbCellSkin2D) = (/'SEG2', 'SEG3', 'SEG4'/)
        character(len=8), parameter :: cellSkin3D(nbCellSkin3D) = (/'TRIA3', 'TRIA6', &
                                                                    'QUAD4', 'QUAD8', &
                                                                    'QUAD9', 'SEG2 ', &
                                                                    'SEG3 ', 'SEG4 '/)
!   ------------------------------------------------------------------------------------------------
!
        mesh = meshZ
        model = modelZ
        nbNodeSlav = 0

! - Get slave cells
        call kineLoadGetSlaveCells(modelZ, meshZ, &
                                   factorKeywordZ, iOcc, &
                                   cellSlavJv, nbCellSlav, cellSlav)

! - Get nodes from skin cells
        nodeSlavJv = '&&NBNLMA.LN'
        call jedetr(nodeSlavJv)
        if (geomDime .eq. 2) then
            call nbnlma(mesh, &
                        nbCellSlav, cellSlav, &
                        nbCellSkin2D, cellSkin2D, &
                        nbNodeSlav, lError, cellNameError)
        elseif (geomDime .eq. 3) then
            call nbnlma(mesh, &
                        nbCellSlav, cellSlav, &
                        nbCellSkin3D, cellSkin3D, &
                        nbNodeSlav, lError, cellNameError)
        else
            ASSERT(ASTER_FALSE)
        end if
        call jedetr('&&NBNLMA.NBN')
        if (lError) then
            call utmess('F', 'CHARGES7_4')
        end if
        call jeveuo(nodeSlavJv, 'L', vi=nodeSlav)
        call jelira(nodeSlavJv, 'LONUTI', nbNodeSlav)

! - Indirection between volumic and skin nodes
        AS_ALLOCATE(vi=nodeSkinToBody, size=meshNbNode)
        do iNodeSlav = 1, nbNodeSlav
            nodeSkinToBody(nodeSlav(iNodeSlav)) = iNodeSlav
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! kineLoadNormFromSkinCells
!
! Compute normals at nodes from skin cells
!
! In  mesh             : mesh
! In  geomDime         : geometric dimension (2 or 3)
! In  nbCellSkin       : number of skin cells
! Ptr cellSkin         : pointer to skin cells
! In  nbNode           : number of nodes
! Ptr node             : pointer to nodes
! Ptr norm             : pointer to list of normals
!
! --------------------------------------------------------------------------------------------------
    subroutine kineLoadNormFromSkinCells(meshZ, geomDime, &
                                         nbCellSkin, cellSkin, &
                                         nbNode, node, norm)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=*), intent(in) :: meshZ
        integer(kind=8), intent(in) :: geomDime
        integer(kind=8), intent(in) :: nbCellSkin
        integer(kind=8), pointer :: cellSkin(:)
        integer(kind=8), intent(in) :: nbNode
        integer(kind=8), pointer :: node(:)
        real(kind=8), pointer :: norm(:)
! - Local
        character(len=24), parameter :: normJv = '&&CANORT.NORMALE'
        character(len=8) :: mesh
!   ------------------------------------------------------------------------------------------------
!
        mesh = meshZ
        call jedetr(normJv)
        call canort(mesh, nbCellSkin, cellSkin, geomDime, nbNode, node, 1)
        call jeveuo(normJv, 'L', vr=norm)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! kineLoadGetSlaveNodes
!
! Get list of slave nodes
!
! In  mesh             : mesh
! In  factorKeyword    : factor keyword
! In  iOcc             : index of factor keyword
! Out nodeSlavJv       : name of JEVEUX object for list of slave nodes
! Out nbNodeSlav       : number of slave nodes
! Ptr nodeSlav         : pointer to list of slave nodes
! In  keywordSuffix    : suffix to keyword to read slave nodes
!
! --------------------------------------------------------------------------------------------------
    subroutine kineLoadGetSlaveNodes(meshZ, &
                                     factorKeywordZ, iOcc, &
                                     nodeSlavJv, nbNodeSlav, nodeSlav, &
                                     keywordSuffixZ_)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=*), intent(in) :: meshZ, factorKeywordZ
        integer(kind=8), intent(in) :: iOcc
        character(len=24), intent(out) :: nodeSlavJv
        integer(kind=8), intent(out) :: nbNodeSlav
        integer(kind=8), pointer :: nodeSlav(:)
        character(len=*), optional, intent(in) :: keywordSuffixZ_
! - Local
        character(len=8) :: mesh, keywordSuffix
!   ------------------------------------------------------------------------------------------------
!
        mesh = meshZ
        nbNodeSlav = 0
        nodeSlavJv = '&&CALIRC.LINONU2'
        keywordSuffix = '_ESCL'
        if (present(keywordSuffixZ_)) then
            keywordSuffix = keywordSuffixZ_
        end if
        call jedetr(nodeSlavJv)
        call getnode(mesh, factorKeywordZ, iOcc, 'F', nodeSlavJv, nbNodeSlav, &
                     suffixZ_=keywordSuffix)
        call jeveuo(nodeSlavJv, 'L', vi=nodeSlav)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! kineLoadGetMasterCells
!
! Get list of master cells
!
! In  model            : model
! In  mesh             : mesh
! In  factorKeyword    : factor keyword
! In  iOcc             : index of factor keyword
! Out cellMastJv       : name of JEVEUX object for list of master cells
! Out nbCellMast       : number of master cells
! Ptr cellMast         : pointer to list of master cells
! In  keywordSuffix    : suffix to keyword to read master cells
!
! --------------------------------------------------------------------------------------------------
    subroutine kineLoadGetMasterCells(modelZ, meshZ, &
                                      factorKeywordZ, iOcc, &
                                      cellMastJv, nbCellMast, cellMast, &
                                      keywordSuffixZ_)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=*), intent(in) :: modelZ, meshZ, factorKeywordZ
        integer(kind=8), intent(in) :: iOcc
        character(len=24), intent(out) :: cellMastJv
        integer(kind=8), intent(out) :: nbCellMast
        integer(kind=8), pointer :: cellMast(:)
        character(len=*), optional, intent(in) :: keywordSuffixZ_
! - Local
        character(len=8) :: mesh, model, keywordSuffix
!   ------------------------------------------------------------------------------------------------
!
        mesh = meshZ
        model = modelZ
        nbCellMast = 0
        cellMastJv = '&&CALIRC.LIMANU1'
        keywordSuffix = '_MAIT'
        if (present(keywordSuffixZ_)) then
            keywordSuffix = keywordSuffixZ_
        end if
        call jedetr(cellMastJv)
        call getelem(mesh, factorKeywordZ, iOcc, 'F', cellMastJv, nbCellMast, keywordSuffix, model)
        call jeveuo(cellMastJv, 'L', vi=cellMast)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! kineLoadGetSlaveCells
!
! Get list of slave cells
!
! In  model            : model
! In  mesh             : mesh
! In  factorKeyword    : factor keyword
! In  iOcc             : index of factor keyword
! Out nbCellSlav       : number of slave cells
! Out cellSlavJv       : name of JEVEUX object for list of slave cells
! Ptr cellSlav         : pointer to list of slave cells
! In  keywordSuffix    : suffix to keyword to read slave cells
!
! --------------------------------------------------------------------------------------------------
    subroutine kineLoadGetSlaveCells(modelZ, meshZ, &
                                     factorKeywordZ, iOcc, &
                                     cellSlavJv, nbCellSlav, cellSlav, &
                                     keywordSuffixZ_)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=*), intent(in) :: modelZ, meshZ, factorKeywordZ
        integer(kind=8), intent(in) :: iOcc
        character(len=24), intent(out) :: cellSlavJv
        integer(kind=8), intent(out) :: nbCellSlav
        integer(kind=8), pointer :: cellSlav(:)
        character(len=*), optional, intent(in) :: keywordSuffixZ_
! - Local
        character(len=8) :: mesh, model, keywordSuffix
!   ------------------------------------------------------------------------------------------------
!
        mesh = meshZ
        model = modelZ
        cellSlavJv = '&&CALIRC.LIMANU2'
        keywordSuffix = '_ESCL'
        if (present(keywordSuffixZ_)) then
            keywordSuffix = keywordSuffixZ_
        end if
        call jedetr(cellSlavJv)
        call getelem(mesh, factorKeywordZ, iOcc, ' ', cellSlavJv, nbCellSlav, keywordSuffix, model)
        if (nbCellSlav .eq. 0) then
            call utmess('F', 'CHARGES7_49')
        end if
        call jeveuo(cellSlavJv, 'L', vi=cellSlav)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! kineLoadGetMasterNodes
!
! Get list of master nodes
!
! In  mesh             : mesh
! In  factorKeyword    : factor keyword
! In  iOcc             : index of factor keyword
! Out nodeMastJv       : name of JEVEUX object for list of master nodes
! Out nbNodeMast       : number of master nodes
! Ptr nodeMast         : pointer to list of master nodes
! In  keywordSuffix    : suffix to keyword to read master nodes
!
! --------------------------------------------------------------------------------------------------
    subroutine kineLoadGetMasterNodes(meshZ, &
                                      factorKeywordZ, iOcc, &
                                      nodeMastJv, nbNodeMast, nodeMast, &
                                      keywordSuffixZ_)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=*), intent(in) :: meshZ, factorKeywordZ
        integer(kind=8), intent(in) :: iOcc
        character(len=24), intent(out) :: nodeMastJv
        integer(kind=8), intent(out) :: nbNodeMast
        integer(kind=8), pointer :: nodeMast(:)
        character(len=*), optional, intent(in) :: keywordSuffixZ_
! - Local
        character(len=8) :: mesh, keywordSuffix
!   ------------------------------------------------------------------------------------------------
!
        mesh = meshZ
        nbNodeMast = 0
        nodeMastJv = '&&CALIRC.LINONU1'
        keywordSuffix = '_MAIT'
        if (present(keywordSuffixZ_)) then
            keywordSuffix = keywordSuffixZ_
        end if
        call jedetr(nodeMastJv)
        call getnode(mesh, factorKeywordZ, iOcc, ' ', nodeMastJv, nbNodeMast, &
                     suffixZ_=keywordSuffix)
        call jeveuo(nodeMastJv, 'L', vi=nodeMast)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! kineLoadElimMult
!
! Suppress duplicates nodes from list of nodes
!
! In  iOcc             : index of factor keyword
! In  nodeJv           : name of JEVEUX object for list of nodes
! IO  nbNode           : number of nodes
! Ptr node             : pointer to list of nodes
! Ptr nodeElim         : list of duplicate nodes
!
! --------------------------------------------------------------------------------------------------
    subroutine kineLoadElimMult(iOcc, nodeJv, nbNode, node, nodeElim)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        integer(kind=8), intent(in) :: iOcc
        character(len=24), intent(in) :: nodeJv
        integer(kind=8), intent(inout) :: nbNode
        integer(kind=8), pointer :: node(:), nodeElim(:)
! - Local
        integer(kind=8) :: iNode, elimNodeIndx, nodeNume
        integer(kind=8), pointer :: nodeCopy(:) => null()
!   ------------------------------------------------------------------------------------------------
!
        elimNodeIndx = 0
        AS_ALLOCATE(vi=nodeCopy, size=nbNode)
        do iNode = 1, nbNode
            nodeNume = node(iNode)
            if (nodeElim(nodeNume) .eq. 0) then
                nodeElim(nodeNume) = 1
                elimNodeIndx = elimNodeIndx+1
                nodeCopy(elimNodeIndx) = nodeNume
            end if
        end do

! - New number of nodes
        nbNode = elimNodeIndx
        if (nbNode .eq. 0) then
            call utmess('F', 'CHARGES7_48', si=iOcc)
        end if

! - New list of nodes
        call jedetr(nodeJv)
        call wkvect(nodeJv, 'V V I', nbNode, vi=node)
        do iNode = 1, nbNode
            node(iNode) = nodeCopy(iNode)
        end do
        AS_DEALLOCATE(vi=nodeCopy)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! kineLoadReadTransf
!
! Read transformation for translation/rotation
!
! In  mesh             : mesh
! In  geomDime         : geometric dimension (2 or 3)
! In  factorKeyword    : factor keyword
! In  iOcc             : index of factor keyword
! In  nodeJv           : name of JEVEUX object for list of nodes
! In  nbNode           : number of nodes
! Out geomJv           : name of JEVEUX object for coordinates of mesh after transformation
!                          WARNING: defined on ALL mesh nodes (>= nbNode)
! Out lApplyRota       : flag when rotation is applied
! Out rotaMatr         : rotation matrix
!
! --------------------------------------------------------------------------------------------------
    subroutine kineLoadReadTransf(meshZ, geomDime, &
                                  factorKeywordZ, iOcc, &
                                  nodeJv, nbNode, geomJv, &
                                  lApplyRota_, rotaMatr_)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=*), intent(in) :: meshZ, factorKeywordZ
        integer(kind=8), intent(in) :: iOcc, geomDime, nbNode
        character(len=24), intent(in) :: nodeJv
        character(len=24), intent(out) :: geomJv
        aster_logical, optional, intent(out) :: lApplyRota_
        real(kind=8), optional, intent(out) :: rotaMatr_(3, 3)
! - Local
        character(len=8) :: mesh
        aster_logical :: lTran, lCent, lAnglNaut
        real(kind=8) :: tranPara(3), centPara(3), anglNautPara(3)
        aster_logical :: lApplyRota
        real(kind=8) :: rotaMatr(3, 3)
!   ------------------------------------------------------------------------------------------------
!
        lApplyRota = ASTER_FALSE
        mesh = meshZ
        rotaMatr = 0.d0
        geomJv = '&&CALIRC.GEOM_TRANS2'

! - Get transformation from user
        call char_read_tran(factorKeywordZ, iOcc, geomDime, &
                            lTran, tranPara, &
                            lCent, centPara, &
                            lAnglNaut, anglNautPara)

! - Apply transformation on nodes
        call calirg(mesh, nbNode, nodeJv, &
                    tranPara, centPara, &
                    lAnglNaut, anglNautPara, &
                    geomJv, lApplyRota, rotaMatr)
        if (present(lApplyRota_)) then
            lApplyRota_ = lApplyRota
            rotaMatr_ = rotaMatr
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! kineLoadApplyTransf
!
! Apply transformation for translation/rotation for shells
!
! In  outputFile       : name of output file for advanced debug print
! In  mesh             : mesh
! In  meshNbNode       : total number of nodes of mesh
! In  factorKeyword    : factor keyword
! In  iOcc             : index of factor keyword
! In  lMastTransf      : apply geometric transformation to master side for COQUE
! In  mastTransf       : geometric transformation to master side for COQUE
! In  lSlavTransf      : apply geometric transformation to slave side for COQUE
! In  slavTransf       : geometric transformation to slave side for COQUE
! In  nbNodeSlav       : number of slave nodes
! Ptr nodeSlav         : pointer to list of slave nodes
! In  lVerbose         : flag to advanced debug print
! In  meshDebugJv      : name of JEVEUX object for advanced debug print
! In  geomSlavJv       : name of JEVEUX object for coordinates of slave nodes after transformation
! Out geomMastJv       : name of JEVEUX object for coordinates of master nodes after transformation
!
! --------------------------------------------------------------------------------------------------
    subroutine kineLoadApplyTransf(outputFileZ, meshZ, meshNbNode, &
                                   factorKeywordZ, iOcc, &
                                   lMastTransf, mastTransf, &
                                   lSlavTransf, slavTransf, &
                                   nbNodeSlav, nodeSlav, &
                                   lVerbose, meshDebugJv, &
                                   geomSlavJv, geomMastJv)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=*), intent(in) :: outputFileZ, meshZ
        integer(kind=8), intent(in) :: meshNbNode
        character(len=*), intent(in) :: factorKeywordZ
        integer(kind=8), intent(in) :: iOcc
        aster_logical, intent(in) :: lMastTransf, lSlavTransf
        character(len=8), intent(in) :: mastTransf(3), slavTransf(3)
        character(len=24), intent(out) :: geomMastJv
        integer(kind=8), intent(in) :: nbNodeSlav
        integer(kind=8), pointer :: nodeSlav(:)
        character(len=24), intent(in) :: geomSlavJv
        aster_logical, intent(in) :: lVerbose
        character(len=8), intent(in) :: meshDebugJv
! - Local
        integer(kind=8) :: iNode, iFunc
        integer(kind=8) :: ier, jvGeomInit, nodeNume
        character(len=24) :: nodeMastJv
        character(len=8) :: mesh, outputFile
        character(len=8), parameter :: funcParaName(3) = (/'X', 'Y', 'Z'/)
        real(kind=8) :: funcEval(3)
        integer(kind=8) :: nbNodeMast
        integer(kind=8), pointer :: nodeMast(:) => null()
        real(kind=8), pointer :: geomMast(:) => null(), geomSlav(:) => null()
        real(kind=8), pointer :: geomVerbose(:) => null()
!   ------------------------------------------------------------------------------------------------
!
        mesh = meshZ
        outputFile = outputFileZ

! - Access to meshes
        call jeveuo(mesh//'.COORDO    .VALE', 'L', jvGeomInit)
        if (lVerbose) then
            call jeveuo(meshDebugJv(1:8)//'.COORDO    .VALE', 'E', vr=geomVerbose)
        end if

! - Transformation on master side
        if (lMastTransf) then
            geomMastJv = '&&CALIRC.GEOM_TRANS1'
            call kineLoadGetMasterNodes(meshZ, &
                                        factorKeywordZ, iOcc, &
                                        nodeMastJv, nbNodeMast, nodeMast)
            call wkvect(geomMastJv, 'V V R', 3*meshNbNode, vr=geomMast)
            do iNode = 1, nbNodeMast
                nodeNume = nodeMast(iNode)
                do iFunc = 1, 3
                    call fointe('F', mastTransf(iFunc), &
                                3, funcParaName, &
                                zr(jvGeomInit+3*(nodeNume-1)), funcEval(iFunc), ier)
                    ASSERT(ier .eq. 0)
                end do
                geomMast(3*(nodeNume-1)+1) = funcEval(1)
                geomMast(3*(nodeNume-1)+2) = funcEval(2)
                geomMast(3*(nodeNume-1)+3) = funcEval(3)
                if (lVerbose) then
                    geomVerbose(3*(nodeNume-1)+1) = funcEval(1)
                    geomVerbose(3*(nodeNume-1)+2) = funcEval(2)
                    geomVerbose(3*(nodeNume-1)+3) = funcEval(3)
                end if
            end do
        end if

! - Transformation on slave side
        if (lSlavTransf) then
            call jeveuo(geomSlavJv, 'E', vr=geomSlav)
            do iNode = 1, nbNodeSlav
                nodeNume = nodeSlav(iNode)
                do iFunc = 1, 3
                    call fointe('F', slavTransf(iFunc), &
                                3, funcParaName, &
                                zr(jvGeomInit+3*(nodeNume-1)), funcEval(iFunc), ier)
                    ASSERT(ier .eq. 0)
                end do
                geomSlav(3*(nodeNume-1)+1) = funcEval(1)
                geomSlav(3*(nodeNume-1)+2) = funcEval(2)
                geomSlav(3*(nodeNume-1)+3) = funcEval(3)
                if (lVerbose) then
                    geomVerbose(3*(nodeNume-1)+1) = funcEval(1)
                    geomVerbose(3*(nodeNume-1)+2) = funcEval(2)
                    geomVerbose(3*(nodeNume-1)+3) = funcEval(3)
                end if
            end do
        end if

! - Print mesh for debug
        if (lVerbose) then
            call kineLoadPrintMesh(meshDebugJv, outputFileZ, iOcc)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! kineLoadSetLCS
!
! Compute locale coordinate system on two sides: slave and master
!
! In  geomDime         : geometric dimension (2 or 3)
! In  nbNodeMast       : number of master nodes
! In  iLink            : current link
! Ptr normSlav         : pointer to list of slave normals
! Ptr nodeSkinToBody   : pointer between body and skin nodes
! In  lApplyRota       : flag when rotation is applied
! In  rotaMatr         : rotation matrix
! IO  kineListRela     : object for list of linear relations
!
! --------------------------------------------------------------------------------------------------
    subroutine kineLoadSetLCS(geomDime, nbNodeMast, &
                              iLink, &
                              normSlav, nodeSkinToBody, &
                              lApplyRota, rotaMatr, &
                              kineListRela)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        integer(kind=8), intent(in) :: geomDime, nbNodeMast
        integer(kind=8), intent(in) :: iLink
        real(kind=8), pointer :: normSlav(:)
        integer(kind=8), pointer :: nodeSkinToBody(:)
        aster_logical, intent(in) :: lApplyRota
        real(kind=8), intent(in) :: rotaMatr(3, 3)
        type(KINE_LIST_RELA), intent(inout) :: kineListRela
! - Local
        integer(kind=8) :: iNodeMast, iGeomDime, jGeomDime
        real(kind=8) :: normal(3)
!   ------------------------------------------------------------------------------------------------
!
        kineListRela%LCSType(1:nbNodeMast+1) = geomDime
        if (lApplyRota) then

! ----- Compute normal with rotation between slave and master sides
            normal = 0.d0
            do iGeomDime = 1, geomDime
                do jGeomDime = 1, geomDime
                    normal(iGeomDime) = normal(iGeomDime)+ &
                                        rotaMatr(jGeomDime, iGeomDime)* &
                                        normSlav(geomDime*(nodeSkinToBody(iLink)-1)+jGeomDime)
                end do
            end do

! ----- Save normal for slave side
            do iGeomDime = 1, geomDime
                kineListRela%LCSVale(iGeomDime) = &
                    normSlav(geomDime*(nodeSkinToBody(iLink)-1)+iGeomDime)
            end do

! ----- Save normal for master side
            do iNodeMast = 2, nbNodeMast+1
                do iGeomDime = 1, geomDime
                    kineListRela%LCSVale(3*(iNodeMast-1)+iGeomDime) = -normal(iGeomDime)
                end do
            end do
        else

! ----- Save normal for slave and master sideS
            do iNodeMast = 1, nbNodeMast+1
                do iGeomDime = 1, geomDime
                    kineListRela%LCSVale(3*(iNodeMast-1)+iGeomDime) = &
                        normSlav(geomDime*(nodeSkinToBody(iLink)-1)+iGeomDime)
                end do
            end do
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! kineLoadPrintMesh
!
! Print mesh in MED file
!
! In  meshDebugJv      : name of JEVEUX object for advanced debug print
! In  fileBaseName     : base name of output file
! In  fileIndx         : index of output file
!
! --------------------------------------------------------------------------------------------------
    subroutine kineLoadPrintMesh(meshDebugJv, fileBaseNameZ, fileIndx)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=8), intent(in) :: meshDebugJv
        character(len=*), intent(in) :: fileBaseNameZ
        integer(kind=8), intent(in) :: fileIndx
! - Local
        character(len=8) :: fileBaseName
        integer(kind=8) :: fileUnit
        character(len=4) :: fileIndxStr
        character(len=80) :: fileName
        character(len=8), parameter :: k8dummy = ' '
        character(len=16), parameter :: k16dummy = ' '
        character(len=19), parameter :: k19dummy = ' '
!   ------------------------------------------------------------------------------------------------
!
        fileBaseName = fileBaseNameZ
        write (fileIndxStr, '(A1,I3.3)') '_', fileIndx
        fileName = 'REPE_OUT/'//fileBaseName(1:8)//fileIndxStr//'_transf_geom.med'
        fileUnit = ulnume()
        if (fileUnit .le. 0) call utmess('F', 'UTILITAI5_10')
        call ulaffe(fileUnit, fileName, ' ', 'N', 'O')
        call irmail('MED', fileUnit, 0, meshDebugJv, ASTER_FALSE, k19dummy, 1, k16dummy)
        call ulopen(-fileUnit, k8dummy, k8dummy, k8dummy, k8dummy)
        call jedetr(meshDebugJv)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! kineLoadGetPhysQuanInfo
!
! Get information about physical quantity
!
! In  physQuanName     : name of physical quantity
! Out nbEc             : number of "coded integers"
! Out dofDZIndx        : index of dof "DZ" in physical quantity
!
! --------------------------------------------------------------------------------------------------
    subroutine kineLoadGetPhysQuanInfo(physQuanNameZ, nbEc, dofDZIndx_)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=*), intent(in) :: physQuanNameZ
        integer(kind=8), intent(out) :: nbEc
        integer(kind=8), optional, intent(out) :: dofDZIndx_
! - Local
        integer(kind=8), parameter :: nbDofMaxi = 320
        character(len=8) :: physQuanDofName(nbDofMaxi)
        character(len=8) :: physQuanName
        integer(kind=8) :: dofDZIndx, jvPhysQuanDofName
        integer(kind=8) :: physQuanSize, physQuanNbDof, iQuanDof
!   ------------------------------------------------------------------------------------------------
!
        physQuanName = physQuanNameZ
        nbEc = 0
        dofDZIndx = 0

! - Access to list of DOF
        physQuanDofName = ' '
        call dismoi('NB_EC', physQuanName, 'GRANDEUR', repi=nbec)
        ASSERT(nbec .le. 11)
        call jelira(jexnom('&CATA.GD.NOMCMP', physQuanName), 'LONMAX', physQuanSize)
        physQuanNbDof = physQuanSize-1
        ASSERT(physQuanSize .le. nbDofMaxi)
        call jeveuo(jexnom('&CATA.GD.NOMCMP', physQuanName), 'L', jvPhysQuanDofName)
        do iQuanDof = 1, physQuanNbDof
            physQuanDofName(iQuanDof) = zk8(jvPhysQuanDofName-1+iQuanDof)
        end do

! - Index of dof DZ in physical quantity
        dofDZIndx = indik8(physQuanDofName, 'DZ', 1, physQuanNbDof)

        if (present(dofDZIndx_)) then
            dofDZIndx_ = dofDZIndx
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! kineLoadNormAtNodes
!
! Compute normals at nodes
!
! In  mesh             : mesh
! In  model            : model
! In  geomDime         : geometric dimension (2 or 3)
! In  factorKeyword    : factor keyword
! In  iOcc             : index of factor keyword
! Ptr coni             : pointer to the pairing of nodes
! In  conrJv           : name of JEVEUX object for normals
!
! --------------------------------------------------------------------------------------------------
    subroutine kineLoadNormAtNodes(meshZ, modelZ, geomDime, &
                                   factorKeywordZ, iOcc, &
                                   coni, conrJvZ)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=*), intent(in) :: meshZ, modelZ
        integer(kind=8), intent(in) :: geomDime
        character(len=*), intent(in) :: factorKeywordZ
        integer(kind=8), intent(in) :: iOcc
        integer(kind=8), pointer :: coni(:)
        character(len=*), intent(in) :: conrJvZ
! - Local
        character(len=8) :: mesh, model
        character(len=24) :: conrJv
        integer(kind=8) :: jvGeom
        integer(kind=8) :: iNode, iDime, nbNode, inoma
        integer(kind=8) :: nodeMastNume, nodeSlavNume
        real(kind=8) :: normMast(3), normSlav(3)
        real(kind=8) :: normDumm, jeu
        real(kind=8), pointer :: normNode(:) => null()
        character(len=24) :: cellMastJv, cellSlavJv
        integer(kind=8) :: nbCellMast, nbCellSlav
        integer(kind=8), pointer :: cellMast(:) => null(), cellSlav(:) => null()

!   ------------------------------------------------------------------------------------------------
!
        mesh = meshZ
        model = modelZ
        conrJv = conrJvZ

! - Access to mesh
        call jeveuo(mesh//'.COORDO    .VALE', 'L', jvGeom)

! - Access to pairing
        nbNode = coni(1)

! - Create object for normal
        call jecroc(jexnum(conrJv, iOcc))
        if (geomDime .eq. 2) then
            call jeecra(jexnum(conrJv, iOcc), 'LONMAX', 12*nbNode)
        elseif (geomDime .eq. 3) then
            call jeecra(jexnum(conrJv, iOcc), 'LONMAX', 22*nbNode)
        else
            ASSERT(ASTER_FALSE)
        end if
        call jeveuo(jexnum(conrJv, iocc), 'E', vr=normNode)

        do iNode = 1, nbNode
! ----- Current pair
            nodeMastNume = coni(1+2*(iNode-1)+1)
            nodeSlavNume = coni(1+2*(iNode-1)+2)

! ----- Get list of master cells
            call kineLoadGetMasterCells(model, mesh, &
                                        factorKeywordZ, iOcc, &
                                        cellMastJv, nbCellMast, cellMast, &
                                        '_1')

! ----- Get list of slave cells
            call kineLoadGetSlaveCells(model, mesh, &
                                       factorKeywordZ, iOcc, &
                                       cellSlavJv, nbCellSlav, cellSlav, &
                                       '_2')

! ----- Compute average normals on cells
            normMast = 0.d0
            normSlav = 0.d0
            call kineLoadAverageNormal(mesh, geomDime, &
                                       nodeMastNume, nbCellMast, cellMast, &
                                       nodeSlavNume, nbCellSlav, cellSlav, &
                                       inoma, normMast, normSlav)

! ----- For cells = 'POI1'
            if (inoma .ne. -1) then
                call normev(normMast, normDumm)
            end if
            if (inoma .ne. -2) then
                call normev(normSlav, normDumm)
            end if
            jeu = 0.d0

!       ANGLE ENTRE NORMALES ET MOYENNE DES NORMALES UNITILE SI POI1
            if ((inoma .ne. -1) .and. (inoma .ne. -2)) then
                normMast = (normMast-normSlav)/2.d0

            else if (inoma .eq. -1) then
                normMast = -normSlav

            else if (inoma .eq. -2) then

            end if

! ----- Compute gap
            jeu = 0.d0
            do iDime = 1, geomDime
                jeu = jeu-zr(jvGeom-1+3*(nodeMastNume-1)+iDime)*normMast(iDime)+ &
                      zr(jvGeom-1+3*(nodeSlavNume-1)+iDime)*normMast(iDime)
            end do

! ----- Save
            do iDime = 1, geomDime
                normNode((2*geomDime+1)*(iNode-1)+iDime) = normMast(iDime)
                normNode((2*geomDime+1)*(iNode-1)+iDime+geomDime) = normSlav(iDime)
            end do
            normNode((2*geomDime+1)*iNode) = jeu
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! kineLoadDeleteNodePair
!
! Delete node pair from SANS_* keyword
!
! In  mesh             : mesh
! In  factorKeyword    : factor keyword
! In  iOcc             : index of factor keyword
! Ptr coni             : pointer to the pairing of nodes
!
! --------------------------------------------------------------------------------------------------
    subroutine kineLoadDeleteNodePair(meshZ, factorKeywordZ, iOcc, coni)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=*), intent(in) :: meshZ, factorKeywordZ
        integer(kind=8), intent(in) :: iOcc
        integer(kind=8), pointer :: coni(:)
! - Local
        character(len=8) :: mesh
        integer(kind=8), parameter :: nbKeywordExcl = 2
        character(len=24), parameter :: keywordExcl(nbKeywordExcl) = (/'SANS_GROUP_NO   ', &
                                                                       'SANS_NOEUD      '/)
        character(len=16), parameter :: keywordExclType(nbKeywordExcl) = (/'GROUP_NO        ', &
                                                                           'NOEUD           '/)
        integer(kind=8) :: iNodeInit, iNodeExcl
        character(len=24), parameter :: nodeExclJv = '&&CAEXNO.LISTENOEUD'
        integer(kind=8), pointer :: nodeExcl(:) => null()
        integer(kind=8) :: nbNodeInit, nbNodeExcl, nbNode
!   ------------------------------------------------------------------------------------------------
!
        mesh = meshZ

! - Get list of exclude nodes
        call reliem(' ', mesh, 'NU_NOEUD', factorKeywordZ, iocc, &
                    nbKeywordExcl, keywordExcl, keywordExclType, nodeExclJv, nbNodeExcl)

! - Change pairing
        if (nbNodeExcl .ne. 0) then
            call jeveuo(nodeExclJv, 'L', vi=nodeExcl)
            nbNodeInit = coni(1)
            nbNode = 0
            do iNodeInit = 1, nbNodeInit
                nbNode = nbNode+1
                do iNodeExcl = 1, nbNodeExcl
                    if ((coni(1+2*(iNodeInit-1)+1) .eq. nodeExcl(iNodeExcl)) .or. &
                        (coni(1+2*(iNodeInit-1)+2) .eq. nodeExcl(iNodeExcl))) then
                        nbNode = nbNode-1
                        goto 2
                    end if
                end do
                coni(1+2*(nbNode-1)+1) = coni(1+2*(iNodeInit-1)+1)
                coni(1+2*(nbNode-1)+2) = coni(1+2*(iNodeInit-1)+2)
2               continue
            end do
            coni(1) = nbNode
        end if
!
        call jedetr(nodeExclJv)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! kineLoadAverageNormals
!
! Compute average normal on cell
!
! In  mesh             : mesh
! In  geomDime         : geometric dimension (2 or 3)
! In  nodeMastNume     : index of master node where to compute normal
! In  nbCellMast       : number of master cells
! Ptr cellMast         : pointer to master cells
! In  nodeSlavNume     : index of slave node where to compute normal
! In  nbCellSlav       : number of slave cells
! Ptr cellSlav         : pointer to slave cells
! Out inoma
! Out normMast         : normal for master side
! Out normSlav         : normal for slave side
!
! --------------------------------------------------------------------------------------------------
    subroutine kineLoadAverageNormal(meshZ, geomDime, &
                                     nodeMastNume, nbCellMast, cellMast, &
                                     nodeSlavNume, nbCellSlav, cellSlav, &
                                     inoma, normMast, normSlav)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=*), intent(in) :: meshZ
        integer(kind=8), intent(in) :: geomDime
        integer(kind=8), intent(in) :: nodeMastNume, nbCellMast
        integer(kind=8), pointer :: cellMast(:)
        integer(kind=8), intent(in) :: nodeSlavNume, nbCellSlav
        integer(kind=8), pointer :: cellSlav(:)
        integer(kind=8), intent(out) :: inoma
        real(kind=8), intent(out) :: normMast(3), normSlav(3)
! - Local
        character(len=8) :: mesh
        integer(kind=8), pointer :: connex(:) => null()
        integer(kind=8), pointer ::cellNbNode(:) => null(), typmail(:) => null()
        integer(kind=8), parameter :: normNorm = 1
        integer(kind=8) :: cellTypeNume, iCellMast, iCellSlav, iNode, nbNode
        integer(kind=8) :: cellMastNume, cellSlavNume
        aster_logical :: mastHasPOI1, slavHasPOI1
        real(kind=8) :: normCell(3), cellCoor(3*MT_NNOMAX2D)
!   ------------------------------------------------------------------------------------------------
!
        mesh = meshZ
        inoma = 0
        normMast = 0.d0
        normSlav = 0.d0

! - Access to mesh
        call jeveuo(mesh//'.TYPMAIL', 'L', vi=typmail)

! - Compute normal for master side
        mastHasPOI1 = ASTER_FALSE
        do iCellMast = 1, nbCellMast
! ----- Current cell
            cellMastNume = cellMast(iCellMast)

! ----- Get properties of current cell
            cellTypeNume = typmail(cellMastNume)
            call jeveuo(jexnum('&CATA.TM.NBNO', cellTypeNume), 'L', vi=cellNbNode)
            nbNode = cellNbNode(1)
            call jeveuo(jexnum(mesh//'.CONNEX', cellMastNume), 'L', vi=connex)
            do iNode = 1, nbNode
                if (connex(iNode) .eq. nodeMastNume) then
                    if (nbNode .eq. 1) then
                        mastHasPOI1 = ASTER_TRUE
                        exit
                    end if
                    inoma = 1
                    ASSERT(nbNode .le. MT_NNOMAX2D)
                    call pacoor(mesh, cellMastNume, nbNode, cellCoor)
                    call canorm(cellCoor, normCell, geomDime, cellTypeNume, normNorm)
                    normMast = normMast+normCell
                    exit
                end if
            end do
        end do

! - Compute normal for slave side
        slavHasPOI1 = ASTER_FALSE
        do iCellSlav = 1, nbCellSlav
! ----- Current cell
            cellSlavNume = cellSlav(iCellSlav)

! ----- Get properties of current cell
            cellTypeNume = typmail(cellSlavNume)
            call jeveuo(jexnum('&CATA.TM.NBNO', cellTypeNume), 'L', vi=cellNbNode)
            nbNode = cellNbNode(1)
            call jeveuo(jexnum(mesh//'.CONNEX', cellSlavNume), 'L', vi=connex)
            do iNode = 1, nbNode
                if (connex(iNode) .eq. nodeSlavNume) then
                    if (nbNode .eq. 1) then
                        slavHasPOI1 = ASTER_TRUE
                        exit
                    end if
                    inoma = 1
                    ASSERT(nbNode .le. MT_NNOMAX2D)
                    call pacoor(mesh, cellSlavNume, nbNode, cellCoor)
                    call canorm(cellCoor, normCell, geomDime, cellTypeNume, normNorm)
                    normSlav = normSlav+normCell
                    exit
                end if
            end do
        end do

! - Treatment of POI1 cells
        if (mastHasPOI1 .and. slavHasPOI1) then
            call utmess('F', 'CHARGES7_9')
        end if
        if (mastHasPOI1) then
            inoma = -1
        else if (slavHasPOI1) then
            inoma = -2
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! kineLoadCheckCmpOnNode
!
! Check components on node
!
! --------------------------------------------------------------------------------------------------
    subroutine kineLoadCheckCmpOnNode(jvPrnm, nodeNume, &
                                      physQuanNbCmp, jvPhysQuanCmpName, &
                                      nbDof, nbec, dofName, &
                                      dofExist, oneDofDoesntExist)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        integer(kind=8), intent(in) :: jvPrnm, nodeNume
        integer(kind=8), intent(in) :: physQuanNbCmp, nbec, nbDof
        character(len=8), pointer :: dofName(:)
        integer(kind=8), intent(in) :: jvPhysQuanCmpName
        aster_logical, pointer :: dofExist(:)
        aster_logical, intent(out) :: oneDofDoesntExist
! - Local
        integer(kind=8) :: iDof, idxCmp
!   ------------------------------------------------------------------------------------------------
!
        dofExist = ASTER_TRUE
        oneDofDoesntExist = ASTER_FALSE
        do iDof = 1, nbDof
            idxCmp = indik8(zk8(jvPhysQuanCmpName), dofName(iDof), 1, physQuanNbCmp)
            ASSERT(idxCmp .gt. 0)
            if (.not. exisdg(zi(jvPrnm-1+(nodeNume-1)*nbec+1), idxCmp)) then
                dofExist(iDof) = ASTER_FALSE
                oneDofDoesntExist = ASTER_TRUE
            end if
        end do
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! kineLoadApplyEccentricity
!
! Apply eccentricity on coefficients
!
! --------------------------------------------------------------------------------------------------
    subroutine kineLoadApplyEccentricity(iDime, nodeNameZ, rotaCmpNameZ, &
                                         coef, coefZero, xyzom, &
                                         iTerm, kineListRela)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        integer(kind=8), intent(in) :: iDime
        character(len=*), intent(in) :: nodeNameZ, rotaCmpNameZ
        real(kind=8), intent(in) :: coef, coefZero
        real(kind=8), intent(in) :: xyzom(3)
        integer(kind=8), intent(inout) :: iTerm
        type(KINE_LIST_RELA), intent(inout) :: kineListRela
! - Local
        character(len=8) :: rotaCmpName, nodeName
!   ------------------------------------------------------------------------------------------------
!
        nodeName = nodeNameZ
        rotaCmpName = rotaCmpNameZ
        if (abs(coef*xyzom(iDime)) .gt. coefZero) then
            iTerm = iTerm+1
            kineListRela%nodeName(iTerm) = nodeName
            kineListRela%coefMultReal(iTerm) = coef*xyzom(iDime)
            kineListRela%dofName(iTerm) = rotaCmpName
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module KineLoadUtility_module
