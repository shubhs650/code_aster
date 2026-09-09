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
! Module for DEBUG for solid-shells elements
!
! ==================================================================================================
!
module SolidShell_Debug_module
! ==================================================================================================
    use SolidShell_type
! ==================================================================================================
    implicit none
! ==================================================================================================
    public  :: dbgObjElemProp, dbgObjElemInte, dbgObjMatePara, &
               dbgObjCellGeom, dbgObjGeomHexa, &
               dbgObjKineHexa, dbgMatrRigiKpg, &
               dbgObjStabHexa, dbgObjBehaPara, &
               dbgObjEpsgHexa, dbgObjEpslHexa
! ==================================================================================================
    private
#include "asterf_types.h"
#include "MeshTypes_type.h"
#include "asterfort/SolidShell_type.h"
#include "asterfort/assert.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! dbgObjElemProp
!
! --------------------------------------------------------------------------------------------------
    subroutine dbgObjElemProp(elemProp)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        type(SSH_ELEM_PROP), intent(in) :: elemProp
!   ------------------------------------------------------------------------------------------------
!
        SSH_DBG_STRG('Objet elemProp')
        WRITE (SSH_DBG_UNIT, *) ' *V* nbDof     :', elemProp%nbDof
        WRITE (SSH_DBG_UNIT, *) ' *V* nbDofGeom :', elemProp%nbDofGeom
        WRITE (SSH_DBG_UNIT, *) ' *V* nbNode    :', elemProp%nbNode
        WRITE (SSH_DBG_UNIT, *) ' *V* nbNodeGeom:', elemProp%nbNodeGeom
        call dbgObjElemInte(elemProp%elemInte)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! dbgObjElemInte
!
! --------------------------------------------------------------------------------------------------
    subroutine dbgObjElemInte(elemInte)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        type(SSH_ELEM_INTE), intent(in) :: elemInte
!   ------------------------------------------------------------------------------------------------
!
        SSH_DBG_STRG('Objet elemInte')
        WRITE (SSH_DBG_UNIT, *) ' *V* inteFami   :', elemInte%inteFami
        WRITE (SSH_DBG_UNIT, *) ' *V* nbIntePoint:', elemInte%nbIntePoint
        ASSERT(elemInte%jvWeight .ne. 0)
        ASSERT(elemInte%jvCoor .ne. 0)
        ASSERT(elemInte%jvShape .ne. 0)
        ASSERT(elemInte%jvDShape .ne. 0)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! dbgObjMatePara
!
! --------------------------------------------------------------------------------------------------
    subroutine dbgObjMatePara(sshMatePara)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        type(SSH_MATE_PARA), intent(in) :: sshMatePara
!   ------------------------------------------------------------------------------------------------
!
        SSH_DBG_STRG('Objet sshMatePara')
        WRITE (SSH_DBG_UNIT, *) ' *S* elemHookeMatrix: ', sum(sshMatePara%elemHookeMatrix)
        ASSERT(sshMatePara%materPara%jvMaterCode .ne. 0)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! dbgObjCellGeom
!
! --------------------------------------------------------------------------------------------------
    subroutine dbgObjCellGeom(cellGeom)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        type(SSH_CELL_GEOM), intent(in) :: cellGeom
!   ------------------------------------------------------------------------------------------------
!
        SSH_DBG_STRG('Objet cellGeom')
        WRITE (SSH_DBG_UNIT, *) ' *S* cellCenterCova: ', sum(cellGeom%cellCenterCova)
        WRITE (SSH_DBG_UNIT, *) ' *S* geomInit      : ', sum(cellGeom%geomInit)
        WRITE (SSH_DBG_UNIT, *) ' *S* geomInitX     : ', sum(cellGeom%geomInitX)
        WRITE (SSH_DBG_UNIT, *) ' *S* geomInitY     : ', sum(cellGeom%geomInitY)
        WRITE (SSH_DBG_UNIT, *) ' *S* geomInitZ     : ', sum(cellGeom%geomInitZ)
        WRITE (SSH_DBG_UNIT, *) ' *S* Jac0          : ', sum(cellGeom%Jac0)
        WRITE (SSH_DBG_UNIT, *) ' *S* JacInv0       : ', sum(cellGeom%JacInv0)
        WRITE (SSH_DBG_UNIT, *) ' *V* detJac0       : ', cellGeom%detJac0
        ASSERT(cellGeom%jvGeom .ne. 0)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! dbgObjBehaPara
!
! --------------------------------------------------------------------------------------------------
    subroutine dbgObjBehaPara(behaPara)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        type(SSH_BEHA_PARA), intent(in) :: behaPara
!   ------------------------------------------------------------------------------------------------
!
        SSH_DBG_STRG('Objet behaPara')
        WRITE (SSH_DBG_UNIT, *) ' *V* lMatr         : ', behaPara%lMatr
        WRITE (SSH_DBG_UNIT, *) ' *V* lVect         : ', behaPara%lVect
        WRITE (SSH_DBG_UNIT, *) ' *V* lSigm         : ', behaPara%lSigm
        WRITE (SSH_DBG_UNIT, *) ' *V* lVari         : ', behaPara%lVari
        WRITE (SSH_DBG_UNIT, *) ' *V* lLarge        : ', behaPara%lLarge
        WRITE (SSH_DBG_UNIT, *) ' *V* lMatrSyme     : ', behaPara%lMatrSyme
        WRITE (SSH_DBG_UNIT, *) ' *V* relaComp      : ', behaPara%relaComp
        WRITE (SSH_DBG_UNIT, *) ' *V* defoComp      : ', behaPara%defoComp
        WRITE (SSH_DBG_UNIT, *) ' *V* typeComp      : ', behaPara%typeComp
        ASSERT(associated(behaPara%compor))
        ASSERT(associated(behaPara%carcri))
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! dbgObjGeomHexa
!
! --------------------------------------------------------------------------------------------------
    subroutine dbgObjGeomHexa(geomHexa)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        type(SSH_GEOM_HEXA), intent(in) :: geomHexa
!   ------------------------------------------------------------------------------------------------
!
        SSH_DBG_STRG('Objet geomHexa')
        WRITE (SSH_DBG_UNIT, *) ' *S* Geom : ', sum(geomHexa%geomCurr)
        WRITE (SSH_DBG_UNIT, *) ' *S* T0   : ', sum(geomHexa%T0)
        WRITE (SSH_DBG_UNIT, *) ' *S* TXI  : ', sum(geomHexa%TXI)
        WRITE (SSH_DBG_UNIT, *) ' *S* TETA : ', sum(geomHexa%TETA)
        WRITE (SSH_DBG_UNIT, *) ' *S* TZETA: ', sum(geomHexa%TZETA)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! dbgObjKineHexa
!
! --------------------------------------------------------------------------------------------------
    subroutine dbgObjKineHexa(kineHexa, smallCstPart_, smallVarPart_)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        type(SSH_KINE_HEXA), intent(in)     :: kineHexa
        aster_logical, optional, intent(in) :: smallCstPart_, smallVarPart_
!   ------------------------------------------------------------------------------------------------
!
        if (present(smallCstPart_)) then
            SSH_DBG_STRG('Objet kineHexa (Small strains - Constant part)')
            WRITE (SSH_DBG_UNIT, *) ' *S* BCova0        : ', sum(kineHexa%BCova0)
            WRITE (SSH_DBG_UNIT, *) ' *S* BCovaXI       : ', sum(kineHexa%BCovaXI)
            WRITE (SSH_DBG_UNIT, *) ' *S* BCovaETA      : ', sum(kineHexa%BCovaETA)
            WRITE (SSH_DBG_UNIT, *) ' *S* BCovaZETA     : ', sum(kineHexa%BCovaZETA)
            WRITE (SSH_DBG_UNIT, *) ' *S* BCovaXIZETA   : ', sum(kineHexa%BCovaXIZETA)
            WRITE (SSH_DBG_UNIT, *) ' *S* BCovaETAZETA  : ', sum(kineHexa%BCovaETAZETA)
            WRITE (SSH_DBG_UNIT, *) ' *S* BCovaZETAZETA : ', sum(kineHexa%BCovaZETAZETA)
            WRITE (SSH_DBG_UNIT, *) ' *S* BCart0        : ', sum(kineHexa%BCart0)
            WRITE (SSH_DBG_UNIT, *) ' *S* BCartZ        : ', sum(kineHexa%BCartZ)
            WRITE (SSH_DBG_UNIT, *) ' *S* BCartZZ       : ', sum(kineHexa%BCartZZ)
            WRITE (SSH_DBG_UNIT, *) ' *S* BCartX        : ', sum(kineHexa%BCartX)
            WRITE (SSH_DBG_UNIT, *) ' *S* BCartY        : ', sum(kineHexa%BCartY)
            WRITE (SSH_DBG_UNIT, *) ' *S* BCartYZ       : ', sum(kineHexa%BCartYZ)
            WRITE (SSH_DBG_UNIT, *) ' *S* BCartXZ       : ', sum(kineHexa%BCartXZ)
        end if
        if (present(smallVarPart_)) then
            SSH_DBG_STRG('Objet kineHexa (Small strains - Variable part)')
            WRITE (SSH_DBG_UNIT, *) ' *S* BCovaEAS      : ', sum(kineHexa%BCovaEAS)
            WRITE (SSH_DBG_UNIT, *) ' *S* BCartEAS      : ', sum(kineHexa%BCartEAS)
            WRITE (SSH_DBG_UNIT, *) ' *S* B             : ', sum(kineHexa%B)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! dbgObjEpsgHexa
!
! --------------------------------------------------------------------------------------------------
    subroutine dbgObjEpsgHexa(epsgHexa)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        type(SSH_EPSG_HEXA), intent(in)     :: epsgHexa
!   ------------------------------------------------------------------------------------------------
!
        SSH_DBG_STRG('Objet epsgHexa')
        WRITE (SSH_DBG_UNIT, *) ' *S* epsgHexa%ECova0        : ', sum(epsgHexa%ECova0)
        WRITE (SSH_DBG_UNIT, *) ' *S* epsgHexa%ECovaXI       : ', sum(epsgHexa%ECovaXI)
        WRITE (SSH_DBG_UNIT, *) ' *S* epsgHexa%ECovaETA      : ', sum(epsgHexa%ECovaETA)
        WRITE (SSH_DBG_UNIT, *) ' *S* epsgHexa%ECovaZETA     : ', sum(epsgHexa%ECovaZETA)
        WRITE (SSH_DBG_UNIT, *) ' *S* epsgHexa%ECovaETAZETA  : ', sum(epsgHexa%ECovaETAZETA)
        WRITE (SSH_DBG_UNIT, *) ' *S* epsgHexa%ECovaXIZETA   : ', sum(epsgHexa%ECovaXIZETA)
        WRITE (SSH_DBG_UNIT, *) ' *S* epsgHexa%ECovaZETAZETA : ', sum(epsgHexa%ECovaZETAZETA)
        WRITE (SSH_DBG_UNIT, *) ' *S* epsgHexa%vale          : ', sum(epsgHexa%vale)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! dbgObjEpslHexa
!
! --------------------------------------------------------------------------------------------------
    subroutine dbgObjEpslHexa(epslHexa)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        type(SSH_EPSL_HEXA), intent(in)     :: epslHexa
!   ------------------------------------------------------------------------------------------------
!
        SSH_DBG_STRG('Objet epslHexa')
        WRITE (SSH_DBG_UNIT, *) ' *S* epslHexa%eigenVale     : ', sum(epslHexa%eigenVale)
        WRITE (SSH_DBG_UNIT, *) ' *S* epslHexa%eigenVect     : ', sum(epslHexa%eigenVect)
        WRITE (SSH_DBG_UNIT, *) ' *S* epslHexa%vale          : ', sum(epslHexa%vale)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! dbgMatrRigiKpg
!
! --------------------------------------------------------------------------------------------------
    subroutine dbgMatrRigiKpg(kpg, zeta, poids, jacob, matrRigi)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        integer(kind=8), intent(in)      :: kpg
        real(kind=8), intent(in) :: zeta, poids, jacob
        real(kind=8), intent(in) :: matrRigi(SSH_NBDOF_MAX, SSH_NBDOF_MAX)
!   ------------------------------------------------------------------------------------------------
!
        WRITE (SSH_DBG_UNIT, *) ' Sur KPG: ', kpg, &
            ' - zeta : ', zeta, ' - poids : ', poids, ' - jacob : ', jacob
        WRITE (SSH_DBG_UNIT, *) '   matrRigi : ', sum(matrRigi)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! dbgObjStabHexa
!
! --------------------------------------------------------------------------------------------------
    subroutine dbgObjStabHexa(Ueff, stabHexa)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        real(kind=8), intent(in)        :: Ueff
        type(SSH_STAB_HEXA), intent(in) :: stabHexa
!   ------------------------------------------------------------------------------------------------
!
        SSH_DBG_STRG('Objet stabHexa')
        WRITE (SSH_DBG_UNIT, *) ' *V* Ueff         : ', Ueff
        WRITE (SSH_DBG_UNIT, *) ' *S* SXI          : ', sum(stabHexa%SXI)
        WRITE (SSH_DBG_UNIT, *) ' *S* SETA         : ', sum(stabHexa%SETA)
        WRITE (SSH_DBG_UNIT, *) ' *S* SETAZETA     : ', sum(stabHexa%SETAZETA)
        WRITE (SSH_DBG_UNIT, *) ' *S* SXIZETA      : ', sum(stabHexa%SXIZETA)
        WRITE (SSH_DBG_UNIT, *) ' *S* matrStabMate : ', sum(stabHexa%matrStabMate)
        WRITE (SSH_DBG_UNIT, *) ' *S* forcStab     : ', sum(stabHexa%forcStab)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module SolidShell_Debug_module
