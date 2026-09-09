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
! Module for elementary computation for solid-shells elements
!
! ==================================================================================================
!
module SolidShell_Elementary_module
! ==================================================================================================
    use Behaviour_module
    use MaterialPara_module
    use MaterialPara_type
    use resi_refe_module, only: RESI_REFE
    use SolidShell_Debug_module
    use SolidShell_Elementary_Hexa_module
    use SolidShell_Geometry_Hexa_module
    use SolidShell_Geometry_module
    use SolidShell_NonLinear_Hexa_module
    use SolidShell_type
    use SolidShell_Utilities_module
! ==================================================================================================
    implicit none
! ==================================================================================================
    public  :: compRigiMatr, compSiefElga, compForcNoda, compNonLinear, &
               compEpsiElga, compEpslElga, &
               compLoad, compMassMatr, compRigiGeomMatr, &
               compRefeForcNoda, compLoadExteStatVari, compEpvcElga
    private :: compElemElasMatrix, &
               initGeomCell, initMatePara, initElemProp, initBehaPara
! ==================================================================================================
    private
#include "asterc/r8vide.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/dmat3d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/SolidShell_type.h"
#include "asterfort/tecach.h"
#include "jeveux.h"
#include "MeshTypes_type.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! initElemProp
!
! Initialization of general properties of finite element
!
! In  inteFami         : name of integration scheme family
! Out elemProp         : general properties of element
!
! --------------------------------------------------------------------------------------------------
    subroutine initElemProp(inteFami, elemProp)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=4), intent(in)     :: inteFami
        type(SSH_ELEM_PROP), intent(out) :: elemProp
! ----- Local
        integer(kind=8) :: npg, nno
        integer(kind=8) :: jvWeight, jvCoor, jvShape, jvDShape
!   ------------------------------------------------------------------------------------------------
!
        if (SSH_DBG_ELEM) SSH_DBG_STRG('> initElemProp')

! ----- Get element parameters
        call elrefe_info(fami=inteFami, npg=npg, nno=nno, &
                         jpoids=jvWeight, jcoopg=jvCoor, &
                         jvf=jvShape, jdfde=jvDShape)

! ----- Set parameters for integration scheme
        elemProp%elemInte%inteFami = inteFami
        elemProp%elemInte%nbIntePoint = npg
        elemProp%elemInte%jvCoor = jvCoor
        elemProp%elemInte%jvWeight = jvWeight
        elemProp%elemInte%jvShape = jvShape
        elemProp%elemInte%jvDShape = jvDShape

! ----- Set main parameters of finite element
        if (nno .eq. 9) then
            elemProp%cellType = SSH_CELL_HEXA
        else
            ASSERT(ASTER_FALSE)
        end if
        elemProp%nbNode = nno
        elemProp%nbNodeGeom = nno-1
        elemProp%nbDofGeom = 3*elemProp%nbNodeGeom
        elemProp%nbDof = elemProp%nbDofGeom+1
!
        if (SSH_DBG_ELEM) SSH_DBG_STRG('< initElemProp')
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! initGeomCell
!
! Initialization of geometric properties of cell
!
! In  elemProp         : general properties of element
! Out cellGeom         : general geometric properties of cell
!
! --------------------------------------------------------------------------------------------------
    subroutine initGeomCell(elemProp, cellGeom)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(SSH_ELEM_PROP), intent(in)  :: elemProp
        type(SSH_CELL_GEOM), intent(out) :: cellGeom
! ----- Local
        integer(kind=8)      :: iDofGeom, iNodeGeom
        real(kind=8) :: detJac0
        real(kind=8) :: area0
!   ------------------------------------------------------------------------------------------------
!
        if (SSH_DBG_ELEM) SSH_DBG_STRG('> initGeomCell')

! ----- Access to field of coordinates
        call jevech('PGEOMER', 'L', cellGeom%jvGeom)

! ----- Set initial geometry
        do iDofGeom = 1, elemProp%nbDofGeom
            cellGeom%geomInit(iDofGeom) = zr(cellGeom%jvGeom+iDofGeom-1)
        end do
        do iNodeGeom = 1, elemProp%nbNodeGeom
            cellGeom%geomInitX(iNodeGeom) = cellGeom%geomInit(3*(iNodeGeom-1)+1)
            cellGeom%geomInitY(iNodeGeom) = cellGeom%geomInit(3*(iNodeGeom-1)+2)
            cellGeom%geomInitZ(iNodeGeom) = cellGeom%geomInit(3*(iNodeGeom-1)+3)
        end do

! ----- Set center
        if (elemProp%cellType == SSH_CELL_HEXA) then
            cellGeom%cellCenterCova = hexaCovaCenter
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- Compute Jacobian matrix at center of element on initial configuration
        call compJacoMatr(elemProp, &
                          cellGeom%geomInit, cellGeom%cellCenterCova, &
                          cellGeom%Jac0, cellGeom%JacInv0, &
                          detJac0)
        cellGeom%detJac0 = abs(detJac0)
! ----- Compute average Thickness
        call compAhmadFrame(cellGeom, area0)
        cellGeom%h0 = abs(detJac0)*8/area0
!
        if (SSH_DBG_ELEM) SSH_DBG_STRG('< initGeomCell')
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! initMatePara
!
! Initialization of properties of material
!
! In  elemProp         : general properties of element
! In  cellGeom         : general geometric properties of cell
! In  timeCurr         : current time
! Out sshMatePara         : parameters of material
!
! --------------------------------------------------------------------------------------------------
    subroutine initMatePara(elemProp, cellGeom, timeCurr, sshMatePara)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(SSH_ELEM_PROP), intent(in)  :: elemProp
        type(SSH_CELL_GEOM), intent(in)  :: cellGeom
        real(kind=8), intent(in)         :: timeCurr
        type(SSH_MATE_PARA), intent(out) :: sshMatePara
! ----- Local
        character(len=8), parameter :: fami = 'RIGI'
        integer(kind=8), parameter :: ndim = 3, kpg = 1, ksp = 1
        integer(kind=8) :: jvMaterc
!   ------------------------------------------------------------------------------------------------
!
        if (SSH_DBG_ELEM) SSH_DBG_STRG('> initMatePara')

! ----- Get material parameters
        call jevech('PMATERC', 'L', jvMaterc)

! ----- Initializations of material parameters on current cell
        call initParaCell(fami, zi(jvMaterc), sshMatePara%materPara)

! ----- Set local coordinate system from user
        call getUserLCS(ndim, elemProp%nbNodeGeom, cellGeom%jvGeom, sshMatePara%materPara%lcsPara)

! ----- At center of element
        call initParaPoin(kpg, ksp, sshMatePara%materPara)

! ----- Compute elasticity matrix at middle of cell
        call compElemElasMatrix(timeCurr, sshMatePara)

!
        if (SSH_DBG_ELEM) SSH_DBG_STRG('< initMatePara')
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! initBehaPara
!
! Initialization of properties of behaviour
!
! In  option           : name of option to compute
! In  elemProp         : general properties of element
! In  cellGeom         : general geometric properties of cell
! IO  sshMatePara      : parameters of material
! Out behaPara         : parameters of behaviour
!
! --------------------------------------------------------------------------------------------------
    subroutine initBehaPara(option, elemProp, cellGeom, sshMatePara, behaPara)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=16), intent(in) :: option
        type(SSH_ELEM_PROP), intent(in) :: elemProp
        type(SSH_CELL_GEOM), intent(in) :: cellGeom
        type(SSH_MATE_PARA), intent(inout) :: sshMatePara
        type(SSH_BEHA_PARA), intent(out) :: behaPara
! ----- Local
        integer(kind=8) :: nno, npg
        integer(kind=8) :: jvWeight, jvShape, jvDShape
        integer(kind=8) :: jvTimeM, jvTimeP
        aster_logical :: lMatrSyme
!   ------------------------------------------------------------------------------------------------
!
        if (SSH_DBG_ELEM) SSH_DBG_STRG('> initBehaPara')

! ----- Properties of finite element
        nno = elemProp%nbNodeGeom
        npg = elemProp%elemInte%nbIntePoint
        jvWeight = elemProp%elemInte%jvWeight
        jvShape = elemProp%elemInte%jvShape
        jvDShape = elemProp%elemInte%jvDShape

! ----- Access to fields time
        call jevech('PINSTMR', 'L', jvTimeM)
        call jevech('PINSTPR', 'L', jvTimeP)

! ----- Access to fields of behaviours parameters
        call jevech('PCOMPOR', 'L', vk16=behaPara%compor)
        call jevech('PCARCRI', 'L', vr=behaPara%carcri)

! ----- Main parameters
        behaPara%relaComp = behaPara%compor(RELA_NAME)
        behaPara%typeComp = behaPara%compor(INCRELAS)
        behaPara%defoComp = behaPara%compor(DEFO)
        if (behaPara%defoComp .eq. 'PETIT') then
            behaPara%lLarge = ASTER_FALSE
        elseif (behaPara%defoComp .eq. 'GDEF_LOG') then
            behaPara%lLarge = ASTER_TRUE
        else
            ASSERT(ASTER_FALSE)
        end if
        lMatrSyme = ASTER_TRUE
        if (nint(behaPara%carcri(CARCRI_MATRSYME)) .gt. 0) then
            lMatrSyme = ASTER_FALSE
        end if
        behaPara%lMatrSyme = lMatrSyme

! ----- Select objects to construct from option name
        call behaviourOption(option, behaPara%compor, &
                             behaPara%lMatr, behaPara%lVect, &
                             behaPara%lVari, behaPara%lSigm)

! ----- Initialisation of behaviour datastructure
        call behaviourInit(behaPara%BEHInteg)

! ----- Set main parameters for behaviour (on cell)
        call behaviourSetParaCell(typmod, option, &
                                  behaPara%compor, behaPara%carcri, &
                                  zr(jvTimeM), zr(jvTimeP), &
                                  sshMatePara%materPara, behaPara%BEHInteg)

! ----- Prepare external state variables (geometry)
        call behaviourPrepESVAGeom(nno, npg, SSH_NDIM, &
                                   jvWeight, jvShape, jvDShape, &
                                   cellGeom%geomInit, &
                                   behaPara%BEHInteg)
!
        if (SSH_DBG_ELEM) SSH_DBG_STRG('< initBehaPara')
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compElemElasMatrix
!
! Compute elasticity matrix at middle of cell
!
! In  timeCurr         : current time
! IO  sshMatePara      : parameters of material
!
! --------------------------------------------------------------------------------------------------
    subroutine compElemElasMatrix(timeCurr, sshMatePara)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        real(kind=8), intent(in) :: timeCurr
        type(SSH_MATE_PARA), intent(inout) :: sshMatePara
!   ------------------------------------------------------------------------------------------------
!
        call dmat3d(sshMatePara%materPara, '+', timeCurr, sshMatePara%elemHookeMatrix)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compRigiMatr
!
! Compute rigidity matrix - RIGI_MECA
!
! --------------------------------------------------------------------------------------------------
    subroutine compRigiMatr()
!   ------------------------------------------------------------------------------------------------
! ----- Local
        character(len=4), parameter :: inteFami = 'RIGI'
        type(SSH_CELL_GEOM) :: cellGeom
        type(SSH_ELEM_PROP) :: elemProp
        type(SSH_MATE_PARA) :: sshMatePara
        real(kind=8) :: matrRigi(SSH_NBDOF_MAX, SSH_NBDOF_MAX), timeCurr
        integer(kind=8) :: jvMatr, i, j, k
!   ------------------------------------------------------------------------------------------------
!
        matrRigi = 0.d0

! ----- Non-sense ! To suppress (see issue30887)
        timeCurr = r8vide()

! ----- Initialization of general properties of finite element
        call initElemProp(inteFami, elemProp)
        if (SSH_DBG_ELEM) call dbgObjElemProp(elemProp)

! ----- Initialization of geometric properties of cell
        call initGeomCell(elemProp, cellGeom)
        if (SSH_DBG_GEOM) call dbgObjCellGeom(cellGeom)

! ----- Initialization of properties of material
        call initMatePara(elemProp, cellGeom, timeCurr, sshMatePara)
        if (SSH_DBG_MATE) call dbgObjMatePara(sshMatePara)

! ----- Compute rigidity matrix
        if (elemProp%cellType .eq. SSH_CELL_HEXA) then
            call compRigiMatrHexa(elemProp, cellGeom, sshMatePara, matrRigi)
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- Save matrix
        call jevech('PMATUUR', 'E', jvMatr)
        k = 0
        do i = 1, elemProp%nbDof
            do j = 1, i
                k = k+1
                zr(jvMatr-1+k) = matrRigi(i, j)
            end do
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compSiefElga
!
! Compute stresses - SIEF_ELGA
!
! --------------------------------------------------------------------------------------------------
    subroutine compSiefElga()
!   ------------------------------------------------------------------------------------------------
! ----- Local
        character(len=4), parameter :: inteFami = 'RIGI'
        type(SSH_CELL_GEOM) :: cellGeom
        type(SSH_ELEM_PROP) :: elemProp
        type(SSH_MATE_PARA) :: sshMatePara
        real(kind=8) :: siefElga(SSH_SIZE_TENS*SSH_NBPG_MAX), timeCurr
        integer(kind=8) :: jvSigm, jvDisp, i
!   ------------------------------------------------------------------------------------------------
!
        siefElga = 0.d0

! ----- Non-sense ! To suppress (see issue30887)
        timeCurr = r8vide()

! ----- Initialization of general properties of finite element
        call initElemProp(inteFami, elemProp)
        if (SSH_DBG_ELEM) call dbgObjElemProp(elemProp)

! ----- Initialization of geometric properties of cell
        call initGeomCell(elemProp, cellGeom)
        if (SSH_DBG_GEOM) call dbgObjCellGeom(cellGeom)

! ----- Initialization of properties of material
        call initMatePara(elemProp, cellGeom, timeCurr, sshMatePara)
        if (SSH_DBG_MATE) call dbgObjMatePara(sshMatePara)

! ----- Get displacements
        call jevech('PDEPLAR', 'L', jvDisp)

! ----- Compute stresses
        if (elemProp%cellType .eq. SSH_CELL_HEXA) then
            call compSiefElgaHexa(elemProp, cellGeom, sshMatePara, zr(jvDisp), &
                                  siefElga)
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- Save stress
        call jevech('PCONTRR', 'E', jvSigm)
        do i = 1, SSH_SIZE_TENS*elemProp%elemInte%nbIntePoint
            zr(jvSigm-1+i) = siefElga(i)
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compForcNoda
!
! Compute nodal forces - FORC_NODA
!
! --------------------------------------------------------------------------------------------------
    subroutine compForcNoda()
!   ------------------------------------------------------------------------------------------------
! ----- Local
        character(len=4), parameter :: inteFami = 'RIGI'
        type(SSH_CELL_GEOM) :: cellGeom
        type(SSH_ELEM_PROP) :: elemProp
        real(kind=8) :: forcNoda(SSH_NBDOF_MAX)
        integer(kind=8) :: jvSigm, jvVect, i
!   ------------------------------------------------------------------------------------------------
!
        forcNoda = 0.d0

! ----- Initialization of general properties of finite element
        call initElemProp(inteFami, elemProp)
        if (SSH_DBG_ELEM) call dbgObjElemProp(elemProp)

! ----- Initialization of geometric properties of cell
        call initGeomCell(elemProp, cellGeom)
        if (SSH_DBG_GEOM) call dbgObjCellGeom(cellGeom)

! ----- Get stresses
        call jevech('PSIEFR', 'L', jvSigm)

! ----- Compute nodal forces
        if (elemProp%cellType .eq. SSH_CELL_HEXA) then
            call compForcNodaHexa(elemProp, cellGeom, &
                                  zr(jvSigm), forcNoda)
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- Save vector
        call jevech('PVECTUR', 'E', jvVect)
        do i = 1, elemProp%nbDof
            zr(jvVect-1+i) = forcNoda(i)
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compNonLinear
!
! Compute non-linear options
!
! In  option           : name of option to compute
!
! --------------------------------------------------------------------------------------------------
    subroutine compNonLinear(option)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=16), intent(in)    :: option
! ----- Local
        character(len=4), parameter :: inteFami = 'RIGI'
        type(SSH_CELL_GEOM) :: cellGeom
        type(SSH_ELEM_PROP) :: elemProp
        type(SSH_MATE_PARA) :: sshMatePara
        type(SSH_BEHA_PARA) :: behaPara
        real(kind=8) :: timeCurr
!   ------------------------------------------------------------------------------------------------
!

! ----- Non-sense ! To suppress (see issue30887)
        timeCurr = r8vide()

! ----- Initialization of general properties of finite element
        call initElemProp(inteFami, elemProp)
        if (SSH_DBG_ELEM) call dbgObjElemProp(elemProp)

! ----- Initialization of geometric properties of cell
        call initGeomCell(elemProp, cellGeom)
        if (SSH_DBG_GEOM) call dbgObjCellGeom(cellGeom)

! ----- Initialization of properties of material
        call initMatePara(elemProp, cellGeom, timeCurr, sshMatePara)
        if (SSH_DBG_MATE) call dbgObjMatePara(sshMatePara)

! ----- Initialization of properties of behaviour
        call initBehaPara(option, elemProp, cellGeom, sshMatePara, behaPara)
        if (SSH_DBG_BEHA) call dbgObjBehaPara(behaPara)

! ----- Compute non-linear options
        if (elemProp%cellType .eq. SSH_CELL_HEXA) then
            call compNonLinearHexa(option, elemProp, cellGeom, behaPara)
        else
            ASSERT(ASTER_FALSE)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compEpsiElga
!
! Compute small strains - EPSI_ELGA
!
! --------------------------------------------------------------------------------------------------
    subroutine compEpsiElga()
!   ------------------------------------------------------------------------------------------------
! ----- Local
        character(len=4), parameter :: inteFami = 'RIGI'
        type(SSH_CELL_GEOM) :: cellGeom
        type(SSH_ELEM_PROP) :: elemProp
        real(kind=8) :: epsiElga(SSH_SIZE_TENS*SSH_NBPG_MAX)
        integer(kind=8) :: jvEpsi, jvDisp, i
!   ------------------------------------------------------------------------------------------------
!
        epsiElga = 0.d0

! ----- Initialization of general properties of finite element
        call initElemProp(inteFami, elemProp)
        if (SSH_DBG_ELEM) call dbgObjElemProp(elemProp)

! ----- Initialization of geometric properties of cell
        call initGeomCell(elemProp, cellGeom)
        if (SSH_DBG_GEOM) call dbgObjCellGeom(cellGeom)

! ----- Get displacements
        call jevech('PDEPLAR', 'L', jvDisp)

! ----- Compute strains
        if (elemProp%cellType .eq. SSH_CELL_HEXA) then
            call compEpsiElgaHexa(elemProp, cellGeom, zr(jvDisp), epsiElga)
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- Save strains
        call jevech('PDEFOPG', 'E', jvEpsi)
        do i = 1, SSH_SIZE_TENS*elemProp%elemInte%nbIntePoint
            zr(jvEpsi-1+i) = epsiElga(i)
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compEpslElga
!
! Compute logarithmic strains - EPSL_ELGA
!
! --------------------------------------------------------------------------------------------------
    subroutine compEpslElga()
!   ------------------------------------------------------------------------------------------------
! ----- Local
        character(len=4), parameter :: inteFami = 'RIGI'
        type(SSH_CELL_GEOM) :: cellGeom
        type(SSH_ELEM_PROP) :: elemProp
        real(kind=8) :: epslElga(SSH_SIZE_TENS*SSH_NBPG_MAX)
        integer(kind=8) :: jvEpsi, jvDisp, i
!   ------------------------------------------------------------------------------------------------
!
        epslElga = 0.d0

! ----- Initialization of general properties of finite element
        call initElemProp(inteFami, elemProp)
        if (SSH_DBG_ELEM) call dbgObjElemProp(elemProp)

! ----- Initialization of geometric properties of cell
        call initGeomCell(elemProp, cellGeom)
        if (SSH_DBG_GEOM) call dbgObjCellGeom(cellGeom)

! ----- Get displacements
        call jevech('PDEPLAR', 'L', jvDisp)

! ----- Compute stresses
        if (elemProp%cellType .eq. SSH_CELL_HEXA) then
            call compEpslElgaHexa(elemProp, cellGeom, zr(jvDisp), epslElga)
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- Save stress
        call jevech('PDEFOPG', 'E', jvEpsi)
        do i = 1, SSH_SIZE_TENS*elemProp%elemInte%nbIntePoint
            zr(jvEpsi-1+i) = epslElga(i)
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compLoad
!
! Compute load - CHAR_MECA_PRES_R / CHAR_MECA_PESA_R / CHAR_MECA_FF3D3D / CHAR_MECA_FR3D3D
!
! In  option           : name of option to compute
!
! --------------------------------------------------------------------------------------------------
    subroutine compLoad(option)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=16), intent(in)   :: option
! ----- Local
        character(len=4), parameter :: inteFami = 'RIGI'
        type(SSH_CELL_GEOM) :: cellGeom
        type(SSH_ELEM_PROP) :: elemProp
        type(SSH_MATE_PARA) :: sshMatePara
        real(kind=8) :: loadNoda(SSH_NBDOF_MAX), timeCurr
        integer(kind=8) :: jvVect, i
!   ------------------------------------------------------------------------------------------------
!
        loadNoda = 0.d0

! ----- Non-sense ! To suppress (see issue30887)
        timeCurr = r8vide()

! ----- Initialization of general properties of finite element
        call initElemProp(inteFami, elemProp)
        if (SSH_DBG_ELEM) call dbgObjElemProp(elemProp)

! ----- Initialization of geometric properties of cell
        call initGeomCell(elemProp, cellGeom)
        if (SSH_DBG_GEOM) call dbgObjCellGeom(cellGeom)

! ----- Initialization of properties of material
        if (option .eq. 'CHAR_MECA_PESA_R') then
            call initMatePara(elemProp, cellGeom, timeCurr, sshMatePara)
            if (SSH_DBG_MATE) call dbgObjMatePara(sshMatePara)
        end if

! ----- Compute load
        if (elemProp%cellType .eq. SSH_CELL_HEXA) then
            call compLoadHexa(elemProp, cellGeom, sshMatePara, option, loadNoda)
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- Save vector
        call jevech('PVECTUR', 'E', jvVect)
        do i = 1, elemProp%nbDof
            zr(jvVect-1+i) = loadNoda(i)
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compMassMatr
!
! Compute mass matrix - MASS_MECA
!
! --------------------------------------------------------------------------------------------------
    subroutine compMassMatr()
!   ------------------------------------------------------------------------------------------------
! ----- Local
        character(len=4), parameter :: inteFami = 'MASS'
        type(SSH_CELL_GEOM) :: cellGeom
        type(SSH_ELEM_PROP) :: elemProp
        type(SSH_MATE_PARA) :: sshMatePara
        real(kind=8) :: matrMass(SSH_NBDOF_MAX, SSH_NBDOF_MAX), timeCurr
        integer(kind=8) :: jvMatr, i, j, k
!   ------------------------------------------------------------------------------------------------
!
        matrMass = 0.d0

! ----- Non-sense ! To suppress (see issue30887)
        timeCurr = r8vide()

! ----- Initialization of general properties of finite element
        call initElemProp(inteFami, elemProp)
        if (SSH_DBG_ELEM) call dbgObjElemProp(elemProp)

! ----- Initialization of geometric properties of cell
        call initGeomCell(elemProp, cellGeom)
        if (SSH_DBG_GEOM) call dbgObjCellGeom(cellGeom)

! ----- Initialization of properties of material
        call initMatePara(elemProp, cellGeom, timeCurr, sshMatePara)
        if (SSH_DBG_MATE) call dbgObjMatePara(sshMatePara)

! ----- Compute mass matrix
        if (elemProp%cellType .eq. SSH_CELL_HEXA) then
            call compMassMatrHexa(elemProp, cellGeom, sshMatePara, matrMass)
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- Save matrix
        call jevech('PMATUUR', 'E', jvMatr)
        k = 0
        do i = 1, elemProp%nbDof
            do j = 1, i
                k = k+1
                zr(jvMatr-1+k) = matrMass(i, j)
            end do
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compRigiGeomMatr
!
! Compute rigidity geometric matrix - RIGI_GEOM
!
! --------------------------------------------------------------------------------------------------
    subroutine compRigiGeomMatr()
!   ------------------------------------------------------------------------------------------------
! ----- Local
        character(len=4), parameter :: inteFami = 'RIGI'
        type(SSH_CELL_GEOM) :: cellGeom
        type(SSH_ELEM_PROP) :: elemProp
        real(kind=8) :: matrRigiGeom(SSH_NBDOF_MAX, SSH_NBDOF_MAX)
        integer(kind=8) :: nbIntePoint
        integer(kind=8) :: jvMatr, jvSigm, i, j, k
!   ------------------------------------------------------------------------------------------------
!
        matrRigiGeom = 0.d0

! ----- Initialization of general properties of finite element
        call initElemProp(inteFami, elemProp)
        if (SSH_DBG_ELEM) call dbgObjElemProp(elemProp)
        nbIntePoint = elemProp%elemInte%nbIntePoint

! ----- Initialization of geometric properties of cell
        call initGeomCell(elemProp, cellGeom)
        if (SSH_DBG_GEOM) call dbgObjCellGeom(cellGeom)

! ----- Get stress tensor
        call jevech('PCONTRR', 'L', jvSigm)

! ----- Compute geometric rigidity matrix
        if (elemProp%cellType .eq. SSH_CELL_HEXA) then
            call compRigiGeomMatrHexa(elemProp, cellGeom, nbIntePoint, zr(jvSigm), matrRigiGeom)
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- Save matrix
        call jevech('PMATUUR', 'E', jvMatr)
        k = 0
        do i = 1, elemProp%nbDof
            do j = 1, i
                k = k+1
                zr(jvMatr-1+k) = matrRigiGeom(i, j)
            end do
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compRefeForcNoda
!
! Compute reference force - REFE_FORC_NODA
!
! --------------------------------------------------------------------------------------------------
    subroutine compRefeForcNoda()
!   ------------------------------------------------------------------------------------------------
! ----- Local
        character(len=4), parameter :: inteFami = 'RIGI'
        type(SSH_CELL_GEOM) :: cellGeom
        type(SSH_ELEM_PROP) :: elemProp
        real(kind=8) :: refeForcNoda(SSH_NBDOF_MAX), sigmRefe
        integer(kind=8) :: jvVect, iDof
        type(RESI_REFE):: refe
!   ------------------------------------------------------------------------------------------------
!
        refeForcNoda = 0.d0

! ----- Initialization of general properties of finite element
        call initElemProp(inteFami, elemProp)
        if (SSH_DBG_ELEM) call dbgObjElemProp(elemProp)

! ----- Initialization of geometric properties of cell
        call initGeomCell(elemProp, cellGeom)
        if (SSH_DBG_GEOM) call dbgObjCellGeom(cellGeom)

! ----- Get reference stress
        call refe%Init('SolidShell')
        sigmRefe = refe%GetRef('SIGM')
        call refe%Check()

! ----- Compute reference force
        if (elemProp%cellType .eq. SSH_CELL_HEXA) then
            call compRefeForcNodaHexa(elemProp, cellGeom, sigmRefe, refeForcNoda)
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- Save vector
        call jevech('PVECTUR', 'E', jvVect)
        do iDof = 1, elemProp%nbDof
            zr(jvVect-1+iDof) = refeForcNoda(iDof)
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compLoadExteStatVariHexa
!
! Compute external state variable load - CHAR_MECA_TEMP_R
!
! --------------------------------------------------------------------------------------------------
    subroutine compLoadExteStatVari(indxVarcStrain)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: indxVarcStrain
! ----- Local
        character(len=4), parameter :: inteFami = 'RIGI'
        type(SSH_CELL_GEOM) :: cellGeom
        type(SSH_ELEM_PROP) :: elemProp
        type(SSH_MATE_PARA) :: sshMatePara
        real(kind=8) :: loadNoda(SSH_NBDOF_MAX), timeCurr
        integer(kind=8) :: jvVect, jvTime, i, iret
!   ------------------------------------------------------------------------------------------------
!
        loadNoda = 0.d0

! ----- Non-sense ! To suppress (see issue30887)
        timeCurr = r8vide()

! ----- Initialization of general properties of finite element
        call initElemProp(inteFami, elemProp)
        if (SSH_DBG_ELEM) call dbgObjElemProp(elemProp)

! ----- Initialization of geometric properties of cell
        call initGeomCell(elemProp, cellGeom)
        if (SSH_DBG_GEOM) call dbgObjCellGeom(cellGeom)

! ----- Get current time
        call tecach('ONO', 'PINSTR', 'L', iret, iad=jvTime)
        if (jvTime .ne. 0) then
            timeCurr = zr(jvTime)
        end if

! ----- Initialization of properties of material
        call initMatePara(elemProp, cellGeom, timeCurr, sshMatePara)
        if (SSH_DBG_MATE) call dbgObjMatePara(sshMatePara)

! ----- Compute external state variable load
        if (elemProp%cellType .eq. SSH_CELL_HEXA) then
            call compLoadExteStatVariHexa(elemProp, cellGeom, sshMatePara, indxVarcStrain, &
                                          loadNoda)
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- Save vector
        call jevech('PVECTUR', 'E', jvVect)
        do i = 1, elemProp%nbDof
            zr(jvVect-1+i) = loadNoda(i)
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compEpvcElga
!
! Compute strains from external state variables - EPVC_ELGA
!
! --------------------------------------------------------------------------------------------------
    subroutine compEpvcElga()
!   ------------------------------------------------------------------------------------------------
! ----- Local
        character(len=4), parameter :: inteFami = 'RIGI'
        integer(kind=8), parameter :: nbCmp = 6
        character(len=16) :: option
        type(SSH_CELL_GEOM) :: cellGeom
        type(SSH_ELEM_PROP) :: elemProp
        type(SSH_MATE_PARA) :: sshMatePara
        real(kind=8) :: epvcElga(SSH_NBPG_MAX, SSH_SIZE_TENS), timeCurr
        real(kind=8) :: epvcElgaAllCmp(SSH_NBPG_MAX, nbCmp)
        integer(kind=8) :: jvEpsi, iIntePoint, iCmp
!   ------------------------------------------------------------------------------------------------
!
        epvcElga = 0.d0

! ----- Non-sense ! To suppress (see issue30887)
        timeCurr = r8vide()

! ----- Initialization of general properties of finite element
        call initElemProp(inteFami, elemProp)
        if (SSH_DBG_ELEM) call dbgObjElemProp(elemProp)

! ----- Initialization of geometric properties of cell
        call initGeomCell(elemProp, cellGeom)
        if (SSH_DBG_GEOM) call dbgObjCellGeom(cellGeom)

! ----- Initialization of properties of material (to suppress, see 30888)
        call initMatePara(elemProp, cellGeom, timeCurr, sshMatePara)
        if (SSH_DBG_MATE) call dbgObjMatePara(sshMatePara)

! ----- Compute strains
        if (elemProp%cellType .eq. SSH_CELL_HEXA) then
            option = 'EPVC_ELGA_TEMP'
            call compEpvcElgaHexa(elemProp, sshMatePara, epvcElga)
            epvcElgaAllCmp(:, 1) = epvcElga(:, 1)
            epvcElgaAllCmp(:, 2) = epvcElga(:, 2)
            epvcElgaAllCmp(:, 3) = epvcElga(:, 3)
            option = 'EPVC_ELGA_SECH'
            call compEpvcElgaHexa(elemProp, sshMatePara, epvcElga)
            epvcElgaAllCmp(:, 4) = epvcElga(:, 1)
            option = 'EPVC_ELGA_HYDR'
            call compEpvcElgaHexa(elemProp, sshMatePara, epvcElga)
            epvcElgaAllCmp(:, 5) = epvcElga(:, 1)
            option = 'EPVC_ELGA_PTOT'
            call compEpvcElgaHexa(elemProp, sshMatePara, epvcElga)
            epvcElgaAllCmp(:, 6) = epvcElga(:, 1)
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- Save strains
        call jevech('PDEFOPG', 'E', jvEpsi)
        do iIntePoint = 1, elemProp%elemInte%nbIntePoint
            do iCmp = 1, nbCmp
                zr(jvEpsi-1+nbCmp*(iIntePoint-1)+iCmp) = epvcElgaAllCmp(iIntePoint, iCmp)
            end do
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module SolidShell_Elementary_module
