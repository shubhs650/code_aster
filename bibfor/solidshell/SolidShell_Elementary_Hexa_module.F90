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
! Module for elementary computations of HEXA cells for solid-shells elements
!
! ==================================================================================================
!
module SolidShell_Elementary_Hexa_module
! ==================================================================================================
    use BehaviourStrain_type
    use MaterialPara_module
    use MaterialPara_type
    use SolidShell_Debug_module
    use SolidShell_Geometry_Hexa_module
    use SolidShell_Kinematic_Hexa_module
    use SolidShell_Stabilization_Hexa_module
    use SolidShell_type
    use SolidShell_Utilities_module
! ==================================================================================================
    implicit none
! ==================================================================================================
    public  :: compRigiMatrHexa, compSiefElgaHexa, compForcNodaHexa, &
               compRigiGeomHexaKpg, &
               compEpsgElgaHexa, compEpsiElgaHexa, compEpslElgaHexa, &
               compLoadHexa, compMassMatrHexa, compRigiGeomMatrHexa, &
               compRefeForcNodaHexa, &
               compLoadExteStatVariHexa, compEpvcElgaHexa
    private :: prodBTSigm, compSiefExteStatVariHexa
! ==================================================================================================
    private
#include "asterc/r8vide.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/btsig.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/epstmc.h"
#include "asterfort/fointe.h"
#include "asterfort/jevecd.h"
#include "asterfort/jevech.h"
#include "asterfort/rcvalb.h"
#include "asterfort/SolidShell_type.h"
#include "asterfort/tecach.h"
#include "asterfort/tefrep.h"
#include "jeveux.h"
#include "MeshTypes_type.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! compRigiMatrHexa
!
! Compute rigidity matrix for HEXA - RIGI_MECA
!
! In  elemProp         : general properties of element
! In  cellGeom         : general geometric properties of cell
! In  sshMatePara      : parameters of material
! Out matrRigi         : rigidity matrix
!
! --------------------------------------------------------------------------------------------------
    subroutine compRigiMatrHexa(elemProp, cellGeom, sshMatePara, matrRigi)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(SSH_ELEM_PROP), intent(in) :: elemProp
        type(SSH_CELL_GEOM), intent(in) :: cellGeom
        type(SSH_MATE_PARA), intent(in) :: sshMatePara
        real(kind=8), intent(out) :: matrRigi(SSH_NBDOF_MAX, SSH_NBDOF_MAX)
! ----- Local
        type(SSH_GEOM_HEXA) :: geomHexa
        type(SSH_KINE_HEXA) :: kineHexa
        type(SSH_STAB_HEXA) :: stabHexa
        real(kind=8) :: zeta, poids, jacob, Ueff
        integer(kind=8) :: nbIntePoint, kpg, jvCoor, jvWeight
        real(kind=8) :: tBDB(SSH_NBDOF_MAX, SSH_NBDOF_MAX)
!   ------------------------------------------------------------------------------------------------
!
        nbIntePoint = elemProp%elemInte%nbIntePoint
        jvCoor = elemProp%elemInte%jvCoor
        jvWeight = elemProp%elemInte%jvWeight

! ----- Prepare geometric quantities
        call initGeomCellHexa(cellGeom, geomHexa)
        if (SSH_DBG_GEOM) call dbgObjGeomHexa(geomHexa)

! ----- Compute gradient matrix in covariant basis
        call compBCovaMatrHexa(geomHexa, kineHexa)

! ----- Compute gradient matrix in cartesian frame
        call compBCartMatrHexa(geomHexa, kineHexa)

        if (SSH_DBG_KINE) call dbgObjKineHexa(kineHexa, smallCstPart_=ASTER_TRUE)

! ----- Loop on Gauss points
        do kpg = 1, nbIntePoint
            zeta = zr(jvCoor-1+3*kpg)
            poids = zr(jvWeight-1+kpg)
            jacob = poids*cellGeom%detJac0

! --------- Compute EAS B matrix in cartesian frame at current Gauss point
            call compBCartEASMatrHexa(zeta, geomHexa, kineHexa)

! --------- Compute B matrix
            call compBMatrHexa(zeta, kineHexa)
            if (SSH_DBG_KINE) call dbgObjKineHexa(kineHexa, smallVarPart_=ASTER_TRUE)

! --------- Compute product tBSB
            call prodBTDB(sshMatePara%elemHookeMatrix, SSH_SIZE_TENS, elemProp%nbDof, kineHexa%B, &
                          tBDB)

! --------- Update matrix
            matrRigi = matrRigi+jacob*tBDB
            if (SSH_DBG_ELEM) call dbgMatrRigiKpg(kpg, zeta, poids, jacob, matrRigi)

        end do

! ----- Effective shear modulus for stabilization (elasticity)
        Ueff = sshMatePara%elemHookeMatrix(5, 5)

! ----- Compute stabilization matrix (material part, elasticity)
        call compStabMatrMateHexa(geomHexa, kineHexa, Ueff, stabHexa)
        if (SSH_DBG_STAB) call dbgObjStabHexa(Ueff, stabHexa)

! ----- Compute matrix
        matrRigi(1:SSH_NBDOFG_HEXA, 1:SSH_NBDOFG_HEXA) = &
            matrRigi(1:SSH_NBDOFG_HEXA, 1:SSH_NBDOFG_HEXA)+stabHexa%matrStabMate
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compSiefElgaHexa
!
! Compute stresses for HEXA - SIEF_ELGA
!
! In  elemProp         : general properties of element
! In  cellGeom         : general geometric properties of cell
! IO  sshMatePara      : parameters of material
! In  disp             : current displacements
! Out siefElga         : stresses at Gauss points
!
! --------------------------------------------------------------------------------------------------
    subroutine compSiefElgaHexa(elemProp, cellGeom, sshMatePara, disp, &
                                siefElga)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(SSH_ELEM_PROP), intent(in) :: elemProp
        type(SSH_CELL_GEOM), intent(in) :: cellGeom
        type(SSH_MATE_PARA), intent(inout) :: sshMatePara
        real(kind=8), intent(in) :: disp(SSH_NBDOF_HEXA)
        real(kind=8), intent(out) :: siefElga(SSH_SIZE_TENS*SSH_NBPG_MAX)
! ----- Local
        type(SSH_GEOM_HEXA) :: geomHexa
        type(SSH_KINE_HEXA) :: kineHexa
        real(kind=8) :: zeta, epsi(SSH_SIZE_TENS)
        integer(kind=8) :: nbIntePoint, kpg, jvCoor
        real(kind=8) :: epvcElga(SSH_NBPG_MAX, SSH_SIZE_TENS)
!   ------------------------------------------------------------------------------------------------
!
        nbIntePoint = elemProp%elemInte%nbIntePoint
        jvCoor = elemProp%elemInte%jvCoor

! ----- Prepare geometric quantities
        call initGeomCellHexa(cellGeom, geomHexa)
        if (SSH_DBG_GEOM) call dbgObjGeomHexa(geomHexa)

! ----- Compute gradient matrix in covariant basis
        call compBCovaMatrHexa(geomHexa, kineHexa)

! ----- Compute gradient matrix in cartesian frame
        call compBCartMatrHexa(geomHexa, kineHexa)
        if (SSH_DBG_KINE) call dbgObjKineHexa(kineHexa, smallCstPart_=ASTER_TRUE)

! ----- Compute strains from external state variables
        call compEpvcElgaHexa(elemProp, sshMatePara, epvcElga)

! ----- Loop on Gauss points
        do kpg = 1, nbIntePoint
            zeta = zr(jvCoor-1+3*kpg)

! --------- Compute EAS B matrix in cartesian frame at current Gauss point
            call compBCartEASMatrHexa(zeta, geomHexa, kineHexa)

! --------- Compute B matrix
            call compBMatrHexa(zeta, kineHexa)
            if (SSH_DBG_KINE) call dbgObjKineHexa(kineHexa, smallVarPart_=ASTER_TRUE)

! --------- Compute small strains
            call compEpsiHexa(kineHexa, disp, epsi)
            epsi = epsi-epvcElga(kpg, :)

! --------- Compute stresses
            siefElga(1+(kpg-1)*SSH_SIZE_TENS:SSH_SIZE_TENS*kpg) = &
                matmul(sshMatePara%elemHookeMatrix, epsi)

        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compForcNodaHexa
!
! Compute nodal forces for HEXA - FORC_NODA
!
! In  elemProp         : general properties of element
! In  cellGeom         : general geometric properties of cell
! In  siefElga         : stresses at Gauss points
! Out forcNoda         : nodal forces
!
! --------------------------------------------------------------------------------------------------
    subroutine compForcNodaHexa(elemProp, cellGeom, &
                                siefElga, forcNoda)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(SSH_ELEM_PROP), intent(in) :: elemProp
        type(SSH_CELL_GEOM), intent(in) :: cellGeom
        real(kind=8), intent(in) :: siefElga(SSH_SIZE_TENS*SSH_NBPG_MAX)
        real(kind=8), intent(out) :: forcNoda(SSH_NBDOF_MAX)
! ----- Locals
        type(SSH_GEOM_HEXA) :: geomHexa
        type(SSH_KINE_HEXA) :: kineHexa
        real(kind=8) :: disp(SSH_NBDOF_HEXA)
        integer(kind=8) :: jvCompor, jvDisp, iretc, iDof
        character(len=16) :: defoComp
!   ------------------------------------------------------------------------------------------------
!

! ----- Select configuration
        call tecach('ONO', 'PCOMPOR', 'L', iretc, iad=jvCompor)
        defoComp = 'PETIT'
        if (iretc .eq. 0) then
            defoComp = zk16(jvCompor-1+DEFO)
        end if

! ----- Update configuration
        if (defoComp .eq. 'PETIT') then
            call initGeomCellHexa(cellGeom, geomHexa)
        else
            call jevech('PDEPLAR', 'L', jvDisp)
            do iDof = 1, SSH_NBDOF_HEXA
                disp(iDof) = zr(jvDisp-1+iDof)
            end do
            call initGeomCellHexa(cellGeom, geomHexa, disp)
        end if
        if (SSH_DBG_GEOM) call dbgObjGeomHexa(geomHexa)

! ----- Compute gradient matrix in covariant basis
        call compBCovaMatrHexa(geomHexa, kineHexa)

! ----- Compute gradient matrix in cartesian frame
        call compBCartMatrHexa(geomHexa, kineHexa)
        if (SSH_DBG_KINE) call dbgObjKineHexa(kineHexa, smallCstPart_=ASTER_TRUE)

! ----- Compute Bt.Sigma
        call prodBTSigm(elemProp, cellGeom, geomHexa, kineHexa, siefElga, forcNoda)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! prodBTSigm
!
! Compute B.Sigm for HEXA - Nodal forces
!
! In  elemProp         : general properties of element
! In  cellGeom         : general geometric properties of cell
! In  siefElga         : stresses at Gauss points
! Out forcNoda         : nodal forces
!
! --------------------------------------------------------------------------------------------------
    subroutine prodBTSigm(elemProp, cellGeom, geomHexa, kineHexa, siefElga, forcNoda)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(SSH_ELEM_PROP), intent(in) :: elemProp
        type(SSH_CELL_GEOM), intent(in) :: cellGeom
        type(SSH_GEOM_HEXA), intent(in) :: geomHexa
        type(SSH_KINE_HEXA), intent(inout) :: kineHexa
        real(kind=8), intent(in) :: siefElga(SSH_SIZE_TENS*SSH_NBPG_MAX)
        real(kind=8), intent(out) :: forcNoda(SSH_NBDOF_MAX)
! ----- Locals
        real(kind=8) :: zeta, poids, jacob
        integer(kind=8) :: nbIntePoint, kpg, jvCoor, jvWeight
!   ------------------------------------------------------------------------------------------------
!
        nbIntePoint = elemProp%elemInte%nbIntePoint
        jvCoor = elemProp%elemInte%jvCoor
        jvWeight = elemProp%elemInte%jvWeight
        forcNoda = 0.d0

! ----- Loop on Gauss points
        do kpg = 1, nbIntePoint
            zeta = zr(jvCoor-1+3*kpg)
            poids = zr(jvWeight-1+kpg)
            jacob = poids*cellGeom%detJac0

! --------- Compute EAS B matrix in cartesian frame at current Gauss point
            call compBCartEASMatrHexa(zeta, geomHexa, kineHexa)

! --------- Compute B matrix
            call compBMatrHexa(zeta, kineHexa)
            if (SSH_DBG_KINE) call dbgObjKineHexa(kineHexa, smallVarPart_=ASTER_TRUE)

! --------- Product BT . sigmEner CUMULATED (see btsig subroutine)
            call btsig(elemProp%nbDof, SSH_SIZE_TENS, jacob, &
                       kineHexa%B, siefElga(1+SSH_SIZE_TENS*(kpg-1)), forcNoda)

        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compRigiGeomHexaKpg
!
! Compute geometric matrix at current Gauss point for HEXA
!
! In  geomHexa         : geometric properties for HEXA cell
! In  cellGeom         : general geometric properties of cell
! In  sigm             : stress tensor at current Gauss point
! Out matrRigi         : rigidity matrix
!
! --------------------------------------------------------------------------------------------------
    subroutine compRigiGeomHexaKpg(geomHexa, zeta, sigm, matrGeom)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(SSH_GEOM_HEXA), intent(in) :: geomHexa
        real(kind=8), intent(in) :: zeta, sigm(SSH_SIZE_TENS)
        real(kind=8), intent(out) :: matrGeom(SSH_NBDOF_MAX, SSH_NBDOF_MAX)
! ----- Locals
        integer(kind=8), parameter :: nbNodeGeom = SSH_NBNODEG_HEXA
        integer(kind=8) :: iNodeGeom, jNodeGeom
        real(kind=8) :: const(SSH_SIZE_TENS)
        real(kind=8) :: GCova0(SSH_SIZE_TENS), GCovaZETA(SSH_SIZE_TENS)
        real(kind=8) :: GCovaZETAZETA(SSH_SIZE_TENS)
        real(kind=8) :: GPinchZETA(SSH_SIZE_TENS), GPinchZZETA(SSH_SIZE_TENS)
        real(kind=8) :: GPinchZZ(SSH_SIZE_TENS)
!   ------------------------------------------------------------------------------------------------
!
        matrGeom = 0.d0

! ----- For "standard" nodes
        do iNodeGeom = 1, nbNodeGeom
            do jNodeGeom = 1, nbNodeGeom

! ------------- Compute gradients for geometric matrix
                call compGCovaMatrHexa(iNodeGeom, jNodeGeom, GCova0, GCovaZETA, GCovaZETAZETA)

! ------------- Compute matrix
                const = matmul(geomHexa%T0, GCova0)+ &
                        zeta*(matmul(geomHexa%T0, GCovaZETA)+ &
                              matmul(geomHexa%TZETA, GCova0))+ &
                        zeta*zeta*(matmul(geomHexa%T0, GCovaZETAZETA)+ &
                                   matmul(geomHexa%TZETA, GCovaZETA))
                matrGeom(3*(iNodeGeom-1)+1:3*(iNodeGeom-1)+3, &
                         3*(jNodeGeom-1)+1:3*(jNodeGeom-1)+3) = &
                    sum(const*sigm)*matr3Iden
            end do
        end do

! ----- For "pinch" node
        GPinchZETA = 0.d0
        GPinchZZETA = 0.d0
        GPinchZZ = 0.d0
        do iNodeGeom = 1, nbNodeGeom
            GPinchZETA(3) = -2.d0*hexaVectG3(iNodeGeom)
            GPinchZETA(5) = -2.d0*hexaVectG1(iNodeGeom)
            GPinchZETA(6) = -2.d0*hexaVectG2(iNodeGeom)
            GPinchZZETA(5) = -2.d0*hexaVectH3(iNodeGeom)
            GPinchZZETA(6) = -2.d0*hexaVectH2(iNodeGeom)
            const = zeta*matmul(geomHexa%T0, GPinchZETA)+ &
                    zeta*zeta*(matmul(geomHexa%T0, GPinchZZETA)+matmul(geomHexa%TZETA, GPinchZETA))
            matrGeom(25, 3*(iNodeGeom-1)+3) = sum(const*sigm)
            matrGeom(3*(iNodeGeom-1)+3, 25) = sum(const*sigm)
        end do
        GPinchZZ(3) = 4.d0
        const = matmul(geomHexa%T0, GPinchZZ)*zeta*zeta
        matrGeom(25, 25) = sum(const*sigm)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compEpsiElgaHexa
!
! Compute small strains for HEXA - EPSI_ELGA
!
! In  elemProp         : general properties of element
! In  cellGeom         : general geometric properties of cell
! In  disp             : current displacements
! Out epsiElga         : small strains at Gauss points
!
! --------------------------------------------------------------------------------------------------
    subroutine compEpsiElgaHexa(elemProp, cellGeom, disp, epsiElga)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(SSH_ELEM_PROP), intent(in) :: elemProp
        type(SSH_CELL_GEOM), intent(in) :: cellGeom
        real(kind=8), intent(in) :: disp(SSH_NBDOF_HEXA)
        real(kind=8), intent(out) :: epsiElga(SSH_SIZE_TENS*SSH_NBPG_MAX)
! ----- Local
        type(SSH_GEOM_HEXA) :: geomHexa
        type(SSH_KINE_HEXA) :: kineHexa
        real(kind=8) :: zeta, epsi(SSH_SIZE_TENS)
        integer(kind=8) :: nbIntePoint, kpg, jvCoor
!   ------------------------------------------------------------------------------------------------
!
        nbIntePoint = elemProp%elemInte%nbIntePoint
        jvCoor = elemProp%elemInte%jvCoor

! ----- Prepare geometric quantities
        call initGeomCellHexa(cellGeom, geomHexa)
        if (SSH_DBG_GEOM) call dbgObjGeomHexa(geomHexa)

! ----- Compute gradient matrix in covariant basis
        call compBCovaMatrHexa(geomHexa, kineHexa)

! ----- Compute gradient matrix in cartesian frame
        call compBCartMatrHexa(geomHexa, kineHexa)
        if (SSH_DBG_KINE) call dbgObjKineHexa(kineHexa, smallCstPart_=ASTER_TRUE)

! ----- Loop on Gauss points
        do kpg = 1, nbIntePoint
            zeta = zr(jvCoor-1+3*kpg)

! --------- Compute EAS B matrix in cartesian frame at current Gauss point
            call compBCartEASMatrHexa(zeta, geomHexa, kineHexa)

! --------- Compute B matrix
            call compBMatrHexa(zeta, kineHexa)
            if (SSH_DBG_KINE) call dbgObjKineHexa(kineHexa, smallVarPart_=ASTER_TRUE)

! --------- Compute small strains
            call compEpsiHexa(kineHexa, disp, epsi)
            epsi(4:6) = 0.5*epsi(4:6)
            epsiElga(1+(kpg-1)*SSH_SIZE_TENS:SSH_SIZE_TENS*kpg) = epsi

        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compEpsgElgaHexa
!
! Compute stresses for HEXA - EPSG_ELGA
!
! In  elemProp         : general properties of element
! In  cellGeom         : general geometric properties of cell
! In  disp             : current displacements
! Out epsgElga         : Euler-Lagrange strains at Gauss points
!
! --------------------------------------------------------------------------------------------------
    subroutine compEpsgElgaHexa(elemProp, cellGeom, disp, epsgElga)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(SSH_ELEM_PROP), intent(in) :: elemProp
        type(SSH_CELL_GEOM), intent(in) :: cellGeom
        real(kind=8), intent(in)        :: disp(SSH_NBDOF_HEXA)
        real(kind=8), intent(out)       :: epsgElga(SSH_SIZE_TENS*SSH_NBPG_MAX)
! ----- Local
        type(SSH_GEOM_HEXA) :: geomHexa
        type(SSH_KINE_HEXA) :: kineHexa
        type(SSH_EPSG_HEXA) :: epsgHexa
        real(kind=8) :: zeta
        integer(kind=8) :: nbIntePoint, kpg, jvCoor
!   ------------------------------------------------------------------------------------------------
!
        nbIntePoint = elemProp%elemInte%nbIntePoint
        jvCoor = elemProp%elemInte%jvCoor

! ----- Prepare geometric quantities
        call initGeomCellHexa(cellGeom, geomHexa)
        if (SSH_DBG_GEOM) call dbgObjGeomHexa(geomHexa)

! ----- Flag for large strains
        kineHexa%lLarge = ASTER_TRUE

        do kpg = 1, nbIntePoint
            zeta = zr(jvCoor-1+3*kpg)

! --------- Compute EAS B matrix in cartesian frame at current Gauss point
            call compBCartEASMatrHexa(zeta, geomHexa, kineHexa)

! --------- Compute Green-Lagrange strains
            call compECovaMatrHexa(cellGeom, disp, epsgHexa)
            call compEpsgHexa(zeta, geomHexa, epsgHexa)
            epsgHexa%vale = epsgHexa%vale+kineHexa%BCartEAS*disp(25)
            if (SSH_DBG_KINE) call dbgObjEpsgHexa(epsgHexa)

            epsgElga(1+(kpg-1)*SSH_SIZE_TENS:SSH_SIZE_TENS*kpg) = epsgHexa%vale(1:SSH_SIZE_TENS)

        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compEpslElgaHexa
!
! Compute stresses for HEXA - EPSL_ELGA
!
! In  elemProp         : general properties of element
! In  cellGeom         : general geometric properties of cell
! In  disp             : current displacements
! Out epslElga         : logarithmic strains at Gauss points
!
! --------------------------------------------------------------------------------------------------
    subroutine compEpslElgaHexa(elemProp, cellGeom, disp, epslElga)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(SSH_ELEM_PROP), intent(in) :: elemProp
        type(SSH_CELL_GEOM), intent(in) :: cellGeom
        real(kind=8), intent(in)        :: disp(SSH_NBDOF_HEXA)
        real(kind=8), intent(out)       :: epslElga(SSH_SIZE_TENS*SSH_NBPG_MAX)
! ----- Local
        type(SSH_GEOM_HEXA) :: geomHexa
        type(SSH_KINE_HEXA) :: kineHexa
        type(SSH_EPSG_HEXA) :: epsgHexa
        type(SSH_EPSL_HEXA) :: epslHexa
        real(kind=8) :: zeta
        integer(kind=8) :: cod, nbIntePoint, kpg, jvCoor
!   ------------------------------------------------------------------------------------------------
!
        nbIntePoint = elemProp%elemInte%nbIntePoint
        jvCoor = elemProp%elemInte%jvCoor

! ----- Prepare geometric quantities
        call initGeomCellHexa(cellGeom, geomHexa)
        if (SSH_DBG_GEOM) call dbgObjGeomHexa(geomHexa)

! ----- Flag for large strains
        kineHexa%lLarge = ASTER_TRUE

        do kpg = 1, nbIntePoint
            zeta = zr(jvCoor-1+3*kpg)

! --------- Compute EAS B matrix in cartesian frame at current Gauss point
            call compBCartEASMatrHexa(zeta, geomHexa, kineHexa)

! --------- Compute Green-Lagrange strains
            call compECovaMatrHexa(cellGeom, disp, epsgHexa)
            call compEpsgHexa(zeta, geomHexa, epsgHexa)
            epsgHexa%vale = epsgHexa%vale+kineHexa%BCartEAS*disp(25)
            if (SSH_DBG_KINE) call dbgObjEpsgHexa(epsgHexa)

! --------- Compute Logarithmic strains
            call compEpslHexa(epsgHexa, epslHexa, cod)
            if (SSH_DBG_KINE) call dbgObjEpslHexa(epslHexa)

            epslElga(1+(kpg-1)*SSH_SIZE_TENS:SSH_SIZE_TENS*kpg) = epslHexa%vale

        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compLoadHexa
!
! Compute load for HEXA - CHAR_MECA_PRES_R / CHAR_MECA_PESA_R / CHAR_MECA_FF3D3D / CHAR_MECA_FR3D3D
!
! In  elemProp         : general properties of element
! In  cellGeom         : general geometric properties of cell
! IO  sshMatePara      : parameters of material
! In  option           : name of option to compute
! Out loadNoda         : nodal force from loads (Neumann)
!
! --------------------------------------------------------------------------------------------------
    subroutine compLoadHexa(elemProp, cellGeom, sshMatePara, option, loadNoda)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(SSH_ELEM_PROP), intent(in) :: elemProp
        type(SSH_CELL_GEOM), intent(in) :: cellGeom
        type(SSH_MATE_PARA), intent(inout) :: sshMatePara
        character(len=16), intent(in) :: option
        real(kind=8), intent(out) :: loadNoda(SSH_NBDOF_MAX)
! ----- Local
        integer(kind=8), parameter :: ksp = 1
        character(len=4) :: inteFami
        integer(kind=8), parameter :: nbNode = SSH_NBNODE_HEXA, nbNodeGeom = SSH_NBNODEG_HEXA
        integer(kind=8) :: iNodeGeom, iDim, nbIntePoint, kpg, kdec, ldec, iret
        integer(kind=8) :: jvPres, jvPesa, jvForc, jvTime
        integer(kind=8) :: jvShape, jvDShape, jvWeight
        real(kind=8) :: presSup, presInf, area, fx, fy, fz, xx, yy, zz
        real(kind=8) :: coefGrav, jacob
        real(kind=8) :: rho(1)
        integer(kind=8) :: valeIret(1)
        integer(kind=8), parameter :: nbPara = 4
        real(kind=8) :: paraVale(nbPara)
        character(len=8), parameter :: paraName(nbPara) = (/'X   ', 'Y   ', 'Z   ', 'INST'/)
!   ------------------------------------------------------------------------------------------------
!
        loadNoda = 0.d0
        inteFami = elemProp%elemInte%inteFami

        if (option .eq. 'CHAR_MECA_PRES_R') then
! --------- Get input fields: for pressure, no node affected -> 0
            call jevecd('PPRESSR', jvPres, 0.d0)

! --------- Compute quantities in Ahmad frame for pinch quantities
            call compAhmadFrame(cellGeom, area)

! --------- Get pressures
            presSup = zr(jvPres-1+1)
            presInf = zr(jvPres-1+2)

! --------- Compute
            loadNoda(25) = loadNoda(25)+ &
                           2.d0*(presInf-presSup)*area/3.d0

        elseif (option .eq. 'CHAR_MECA_PESA_R') then

! --------- Get references to integration scheme
            inteFami = elemProp%elemInte%inteFami
            nbIntePoint = elemProp%elemInte%nbIntePoint
            jvShape = elemProp%elemInte%jvShape
            jvDShape = elemProp%elemInte%jvDShape
            jvWeight = elemProp%elemInte%jvWeight

! --------- Get input field: for gravity
            call jevech('PPESANR', 'L', jvPesa)

! --------- Loop on Gauss points
            do kpg = 1, nbIntePoint
! ------------- Initializations of material parameters on current integration point
                call initParaPoin(kpg, ksp, sshMatePara%materPara)

! ------------- Get density
                call rcvalb(sshMatePara%materPara%schemePara%fami, &
                            sshMatePara%materPara%schemePara%kpg, &
                            sshMatePara%materPara%schemePara%ksp, &
                            '+', sshMatePara%materPara%jvMaterCode, &
                            ' ', 'ELAS', 0, ' ', [0.d0], &
                            1, 'RHO', rho, valeIret(1), 1)

! ------------- Compute gradient matrix
                call dfdm3d(nbNode, kpg, jvWeight, jvDShape, zr(cellGeom%jvGeom), jacob)

! ------------- Compute
                coefGrav = rho(1)*jacob*zr(jvPesa)
                do iNodeGeom = 1, nbNodeGeom
                    do iDim = 1, 3
                        loadNoda(3*(iNodeGeom-1)+iDim) = loadNoda(3*(iNodeGeom-1)+iDim)+ &
                                                         coefGrav* &
                                                         zr(jvShape+(kpg-1)*nbNode+iNodeGeom-1)* &
                                                         zr(jvPesa+iDim)
                    end do
                end do
            end do

        elseif (option .eq. 'CHAR_MECA_FR3D3D') then
! --------- Get references to integration scheme
            nbIntePoint = elemProp%elemInte%nbIntePoint
            jvShape = elemProp%elemInte%jvShape
            jvDShape = elemProp%elemInte%jvDShape
            jvWeight = elemProp%elemInte%jvWeight

! --------- Get input fields
            call tefrep(option, 'PFR3D3D', jvForc)
            call jevech('PFR3D3D', 'L', jvForc)

! --------- Loop on Gauss points
            do kpg = 1, nbIntePoint
                ldec = (kpg-1)*nbNode

! ------------- Compute gradient matrix
                call dfdm3d(nbNode, kpg, jvWeight, jvDShape, zr(cellGeom%jvGeom), jacob)

! ------------- Evaluate force at Gauss points from node values
                fx = 0.d0
                fy = 0.d0
                fz = 0.d0
                do iNodeGeom = 1, nbNodeGeom
                    kdec = SSH_NDIM*(iNodeGeom-1)
                    fx = fx+zr(jvShape-1+ldec+iNodeGeom)*zr(jvForc+kdec)
                    fy = fy+zr(jvShape-1+ldec+iNodeGeom)*zr(jvForc+kdec+1)
                    fz = fz+zr(jvShape-1+ldec+iNodeGeom)*zr(jvForc+kdec+2)
                end do

! ------------- Compute load
                do iNodeGeom = 1, nbNodeGeom
                    kdec = SSH_NDIM*(iNodeGeom-1)
                    loadNoda(kdec+1) = loadNoda(kdec+1)+ &
                                       jacob*fx*zr(jvShape-1+ldec+iNodeGeom)
                    loadNoda(kdec+2) = loadNoda(kdec+2)+ &
                                       jacob*fy*zr(jvShape-1+ldec+iNodeGeom)
                    loadNoda(kdec+3) = loadNoda(kdec+3)+ &
                                       jacob*fz*zr(jvShape-1+ldec+iNodeGeom)
                end do
            end do

        elseif (option .eq. 'CHAR_MECA_FF3D3D') then
! --------- Get references to integration scheme
            nbIntePoint = elemProp%elemInte%nbIntePoint
            jvShape = elemProp%elemInte%jvShape
            jvDShape = elemProp%elemInte%jvDShape
            jvWeight = elemProp%elemInte%jvWeight

! --------- Get input fields
            call jevech('PFF3D3D', 'L', jvForc)
            call jevech('PINSTR', 'L', jvTime)
            paraVale(4) = zr(jvTime)

! --------- Loop on Gauss points
            do kpg = 1, nbIntePoint
                ldec = (kpg-1)*nbNode

! ------------- Compute gradient matrix
                call dfdm3d(nbNode, kpg, jvWeight, jvDShape, zr(cellGeom%jvGeom), jacob)

! ------------- Compute coordinates of current Gauss point
                xx = 0.d0
                yy = 0.d0
                zz = 0.d0
                do iNodeGeom = 1, nbNodeGeom
                    xx = xx+zr(cellGeom%jvGeom+3*iNodeGeom-3)*zr(jvShape-1+ldec+iNodeGeom)
                    yy = yy+zr(cellGeom%jvGeom+3*iNodeGeom-2)*zr(jvShape-1+ldec+iNodeGeom)
                    zz = zz+zr(cellGeom%jvGeom+3*iNodeGeom-1)*zr(jvShape-1+ldec+iNodeGeom)
                end do
                paraVale(1) = xx
                paraVale(2) = yy
                paraVale(3) = zz

! ------------- Evaluate force at Gauss points from function
                call fointe('FM', zk8(jvForc), nbPara, paraName, paraVale, fx, iret)
                call fointe('FM', zk8(jvForc+1), nbPara, paraName, paraVale, fy, iret)
                call fointe('FM', zk8(jvForc+2), nbPara, paraName, paraVale, fz, iret)

! ------------- Compute load
                do iNodeGeom = 1, nbNodeGeom
                    kdec = SSH_NDIM*(iNodeGeom-1)
                    loadNoda(kdec+1) = loadNoda(kdec+1)+ &
                                       jacob*fx*zr(jvShape-1+ldec+iNodeGeom)
                    loadNoda(kdec+2) = loadNoda(kdec+2)+ &
                                       jacob*fy*zr(jvShape-1+ldec+iNodeGeom)
                    loadNoda(kdec+3) = loadNoda(kdec+3)+ &
                                       jacob*fz*zr(jvShape-1+ldec+iNodeGeom)
                end do
            end do

        else
            ASSERT(ASTER_FALSE)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compMassMatrHexa
!
! Compute mass matrix for HEXA - MASS_MECA
!
! In  elemProp         : general properties of element
! In  cellGeom         : general geometric properties of cell
! IO  sshMatePara      : parameters of material
! Out matrMass         : mass matrix
!
! --------------------------------------------------------------------------------------------------
    subroutine compMassMatrHexa(elemProp, cellGeom, sshMatePara, matrMass)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(SSH_ELEM_PROP), intent(in) :: elemProp
        type(SSH_CELL_GEOM), intent(in) :: cellGeom
        type(SSH_MATE_PARA), intent(inout) :: sshMatePara
        real(kind=8), intent(out) :: matrMass(SSH_NBDOF_MAX, SSH_NBDOF_MAX)
! ----- Local
        integer(kind=8), parameter :: nbNodeGeom = SSH_NBNODEG_HEXA, ksp = 1
        integer(kind=8) :: iNodeGeom, jNodeGeom
        real(kind=8) :: poids, jacob, XI(3)
        real(kind=8) :: rho(1), N(SSH_NBNODEG_HEXA), NPinch
        integer(kind=8) :: valeIret(1)
        character(len=4) :: inteFami
        integer(kind=8) :: nbIntePoint, kpg, jvCoor, jvWeight
        real(kind=8) :: matrMassPt(SSH_NBDOF_MAX, SSH_NBDOF_MAX)
!   ------------------------------------------------------------------------------------------------
!
        matrMass = 0.d0
        nbIntePoint = elemProp%elemInte%nbIntePoint
        inteFami = elemProp%elemInte%inteFami
        jvCoor = elemProp%elemInte%jvCoor
        jvWeight = elemProp%elemInte%jvWeight

! ----- Loop on Gauss points
        do kpg = 1, nbIntePoint
            XI(1) = zr(jvCoor+3*(kpg-1)-1+1)
            XI(2) = zr(jvCoor+3*(kpg-1)-1+2)
            XI(3) = zr(jvCoor+3*(kpg-1)-1+3)
            poids = zr(jvWeight-1+kpg)
            jacob = poids*cellGeom%detJac0

! --------- Initializations of material parameters on current integration point
            call initParaPoin(kpg, ksp, sshMatePara%materPara)

! --------- Get density
            call rcvalb(sshMatePara%materPara%schemePara%fami, &
                        sshMatePara%materPara%schemePara%kpg, &
                        sshMatePara%materPara%schemePara%ksp, &
                        '+', sshMatePara%materPara%jvMaterCode, &
                        ' ', 'ELAS', 0, ' ', [0.d0], &
                        1, 'RHO', rho, valeIret(1), 1)

! --------- Construction of shape functions
            N = hexaVectS1+ &
                XI(1)*hexaVectG1+XI(2)*hexaVectG2+XI(3)*hexaVectG3+ &
                XI(1)*XI(2)*hexaVectH1+ &
                XI(3)*XI(2)*hexaVectH2+ &
                XI(1)*XI(3)*hexaVectH3+ &
                XI(1)*XI(2)*XI(3)*hexaVectH4
            NPinch = 1.d0-XI(3)*XI(3)

! --------- Compute matrix on volumic nodes
            matrMassPt = 0.d0
            do iNodeGeom = 1, nbNodeGeom
                do jNodeGeom = 1, nbNodeGeom
                    matrMassPt(3*(iNodeGeom-1)+1:3*(iNodeGeom-1)+3, &
                               3*(jNodeGeom-1)+1:3*(jNodeGeom-1)+3) = &
                        rho(1)*jacob*N(iNodeGeom)*N(jNodeGeom)*matr3Iden
                end do
            end do

! --------- Compute matrix on pinch node
            do iNodeGeom = 1, nbNodeGeom
                matrMassPt(SSH_NBDOF_HEXA, 3*(iNodeGeom-1)+3) = rho(1)*jacob*N(iNodeGeom)*NPinch
                matrMassPt(3*(iNodeGeom-1)+3, SSH_NBDOF_HEXA) = rho(1)*jacob*N(iNodeGeom)*NPinch
            end do
            matrMassPt(SSH_NBDOF_HEXA, SSH_NBDOF_HEXA) = rho(1)*jacob*NPinch*NPinch

! --------- Update matrix
            matrMass = matrMass+matrMassPt

        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compRigiGeomMatrHexa
!
! Compute geometric rigidity matrix for HEXA - RIGI_GEOM
!
! In  elemProp         : general properties of element
! In  cellGeom         : general geometric properties of cell
! In  nbIntePoint      : number of integration points on cell
! In  sigm             : stress tensor at current Gauss point
! Out matrRigiGeom     : geometric rigidity matrix
!
! --------------------------------------------------------------------------------------------------
    subroutine compRigiGeomMatrHexa(elemProp, cellGeom, nbIntePoint, sigm, matrRigiGeom)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(SSH_ELEM_PROP), intent(in) :: elemProp
        type(SSH_CELL_GEOM), intent(in) :: cellGeom
        integer(kind=8), intent(in)             :: nbIntePoint
        real(kind=8), intent(in)        :: sigm(SSH_SIZE_TENS, nbIntePoint)
        real(kind=8), intent(out)       :: matrRigiGeom(SSH_NBDOF_MAX, SSH_NBDOF_MAX)
! ----- Local
        type(SSH_GEOM_HEXA) :: geomHexa
        real(kind=8) :: zeta, poids, jacob
        integer(kind=8) :: kpg, jvCoor, jvWeight
        real(kind=8) :: matrRigiGeomPt(SSH_NBDOF_MAX, SSH_NBDOF_MAX)
!   ------------------------------------------------------------------------------------------------
!
        matrRigiGeom = 0.d0
        jvCoor = elemProp%elemInte%jvCoor
        jvWeight = elemProp%elemInte%jvWeight

! ----- Prepare geometric quantities
        call initGeomCellHexa(cellGeom, geomHexa)
        if (SSH_DBG_GEOM) call dbgObjGeomHexa(geomHexa)

! ----- Loop on Gauss points
        do kpg = 1, nbIntePoint
            zeta = zr(jvCoor-1+3*kpg)
            poids = zr(jvWeight-1+kpg)
            jacob = poids*cellGeom%detJac0

! --------- Compute geometric matrix at current Gauss point for HEXA
            call compRigiGeomHexaKpg(geomHexa, zeta, sigm(:, kpg), matrRigiGeomPt)

! --------- Update matrix
            matrRigiGeom = matrRigiGeom+jacob*matrRigiGeomPt

        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compRefeForcNodaHexa
!
! Compute reference force for HEXA - REFE_FORC_NODA
!
! In  elemProp         : general properties of element
! In  cellGeom         : general geometric properties of cell
! In  sigeRefe         : reference value for stress
! Out refeForcNoda     : reference force
!
! --------------------------------------------------------------------------------------------------
    subroutine compRefeForcNodaHexa(elemProp, cellGeom, sigmRefe, refeForcNoda)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(SSH_ELEM_PROP), intent(in) :: elemProp
        type(SSH_CELL_GEOM), intent(in) :: cellGeom
        real(kind=8), intent(in)        :: sigmRefe
        real(kind=8), intent(out)       :: refeForcNoda(SSH_NBDOF_MAX)
! ----- Locals
        type(SSH_GEOM_HEXA) :: geomHexa
        type(SSH_KINE_HEXA) :: kineHexa
        real(kind=8) :: siefElga(SSH_SIZE_TENS*SSH_NBPG_MAX)
        real(kind=8) :: forcNoda(SSH_NBDOF_MAX)
        integer(kind=8) :: iDof, iCmp, nbIntePoint
!   ------------------------------------------------------------------------------------------------
!
        refeForcNoda = 0.d0
        nbIntePoint = elemProp%elemInte%nbIntePoint

! ----- Prepare geometric quantities
        call initGeomCellHexa(cellGeom, geomHexa)
        if (SSH_DBG_GEOM) call dbgObjGeomHexa(geomHexa)

! ----- Compute gradient matrix in covariant basis
        call compBCovaMatrHexa(geomHexa, kineHexa)

! ----- Compute gradient matrix in cartesian frame
        call compBCartMatrHexa(geomHexa, kineHexa)
        if (SSH_DBG_KINE) call dbgObjKineHexa(kineHexa, smallCstPart_=ASTER_TRUE)

! ----- Compute force
        do iCmp = 1, SSH_SIZE_TENS*nbIntePoint
            siefElga = 0.d0
            siefElga(iCmp) = sigmRefe
            call prodBTSigm(elemProp, cellGeom, geomHexa, kineHexa, siefElga, forcNoda)
            do iDof = 1, elemProp%nbDof
                refeForcNoda(iDof) = refeForcNoda(iDof)+abs(forcNoda(iDof))
            end do
        end do
        do iDof = 1, elemProp%nbDof
            refeForcNoda(iDof) = refeForcNoda(iDof)/nbIntePoint
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compLoadExteStatVariHexa
!
! Compute external state variable load for HEXA - CHAR_MECA_TEMP_R
!
! In  elemProp         : general properties of element
! In  cellGeom         : general geometric properties of cell
! IO  sshMatePara         : parameters of material
! Out loadNoda         : nodal force from loads (Neumann)
!
! --------------------------------------------------------------------------------------------------
    subroutine compLoadExteStatVariHexa(elemProp, cellGeom, sshMatePara, indxVarcStrain, loadNoda)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(SSH_ELEM_PROP), intent(in) :: elemProp
        type(SSH_CELL_GEOM), intent(in) :: cellGeom
        type(SSH_MATE_PARA), intent(inout) :: sshMatePara
        integer(kind=8), intent(in) :: indxVarcStrain
        real(kind=8), intent(out) :: loadNoda(SSH_NBDOF_MAX)
! ----- Local
        type(SSH_GEOM_HEXA) :: geomHexa
        type(SSH_KINE_HEXA) :: kineHexa
        real(kind=8) :: siefElga(SSH_SIZE_TENS*SSH_NBPG_MAX)
!   ------------------------------------------------------------------------------------------------
!
        loadNoda = 0.d0

! ----- Compute stresses from external state variables
        call compSiefExteStatVariHexa(elemProp, sshMatePara, indxVarcStrain, siefElga)

! ----- Prepare geometric quantities
        call initGeomCellHexa(cellGeom, geomHexa)
        if (SSH_DBG_GEOM) call dbgObjGeomHexa(geomHexa)

! ----- Compute gradient matrix in covariant basis
        call compBCovaMatrHexa(geomHexa, kineHexa)

! ----- Compute gradient matrix in cartesian frame
        call compBCartMatrHexa(geomHexa, kineHexa)
        if (SSH_DBG_KINE) call dbgObjKineHexa(kineHexa, smallCstPart_=ASTER_TRUE)

! ----- Compute Bt.Sigma
        call prodBTSigm(elemProp, cellGeom, geomHexa, kineHexa, siefElga, loadNoda)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compSiefExteStatVariHexa
!
! Compute stresses from external state variables for HEXA - SIEF_ELGA
!
! In  elemProp         : general properties of element
! IO  sshMatePara      : parameters of material
! Out siefElga         : stresses at Gauss points
!
! --------------------------------------------------------------------------------------------------
    subroutine compSiefExteStatVariHexa(elemProp, sshMatePara, indxVarcStrain, siefElga)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(SSH_ELEM_PROP), intent(in) :: elemProp
        type(SSH_MATE_PARA), intent(inout) :: sshMatePara
        integer(kind=8), intent(in) :: indxVarcStrain
        real(kind=8), intent(out) :: siefElga(SSH_SIZE_TENS*SSH_NBPG_MAX)
! ----- Local
        integer(kind=8), parameter :: ksp = 1, ndim = 3
        integer(kind=8) :: nbIntePoint, kpg
        real(kind=8) :: epsiVarc(SSH_SIZE_TENS), timeCurr
        type(All_Varc_Strain) :: allVarcStrain
!   ------------------------------------------------------------------------------------------------
!
        nbIntePoint = elemProp%elemInte%nbIntePoint
        siefElga = 0.d0

! ----- Non-sense ! To suppress (see issue30887)
        timeCurr = r8vide()

! ----- Loop on Gauss points
        do kpg = 1, nbIntePoint
! --------- Initializations of material parameters on current integration point
            call initParaPoin(kpg, ksp, sshMatePara%materPara)

! --------- Compute strains from external state variables
            epsiVarc = 0.d0
            call epstmc(sshMatePara%materPara, '+', timeCurr, ndim, &
                        indxVarcStrain, allVarcStrain, &
                        epsiVarc)

! --------- Shear components with sqrt(2)
            epsiVarc(4:6) = 2.d0*epsiVarc(4:6)

! --------- Compute stresses from external state variables
            siefElga(1+(kpg-1)*SSH_SIZE_TENS:SSH_SIZE_TENS*kpg) = &
                matmul(sshMatePara%elemHookeMatrix, epsiVarc)

        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compEpvcElgaHexa
!
! Compute strains from external state variables for HEXA - EPVC_ELGA
!
! In  elemProp         : general properties of element
! IO  sshMatePara      : parameters of material
! Out epvcElga         : strains from external state variables
!
! --------------------------------------------------------------------------------------------------
    subroutine compEpvcElgaHexa(elemProp, sshMatePara, epvcElga)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(SSH_ELEM_PROP), intent(in) :: elemProp
        type(SSH_MATE_PARA), intent(inout) :: sshMatePara
        real(kind=8), intent(out) :: epvcElga(SSH_NBPG_MAX, SSH_SIZE_TENS)
! ----- Local
        integer(kind=8), parameter :: ksp = 1, ndim = 3
        integer(kind=8) :: nbIntePoint, kpg
        real(kind=8) :: epsiVarc(SSH_SIZE_TENS), timeCurr
        type(All_Varc_Strain) :: allVarcStrain
!   ------------------------------------------------------------------------------------------------
!
        nbIntePoint = elemProp%elemInte%nbIntePoint
        epvcElga = 0.d0

! ----- Non-sense ! To suppress (see issue30887)
        timeCurr = r8vide()

! ----- Loop on Gauss points
        do kpg = 1, nbIntePoint
! --------- Initializations of material parameters on current integration point
            call initParaPoin(kpg, ksp, sshMatePara%materPara)

! --------- Compute strains from external state variables
            epsiVarc = 0.d0
            call epstmc(sshMatePara%materPara, '+', timeCurr, ndim, &
                        VARC_STRAIN_ALL, allVarcStrain, &
                        epsiVarc)

! --------- Shear components with sqrt(2)
            epsiVarc(4:6) = 2.d0*epsiVarc(4:6)

! --------- Copy strains
            epvcElga(kpg, 1:SSH_SIZE_TENS) = epsiVarc

        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module SolidShell_Elementary_Hexa_module
