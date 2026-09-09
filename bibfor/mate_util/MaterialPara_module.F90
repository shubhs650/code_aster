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
! Types for the management of material parameters in behaviour
!
! ==================================================================================================
!
module MaterialPara_module
! ==================================================================================================
    use MaterialPara_type
    use coorSyst_module, only: hasOrieField
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: initParaCell, initParaPoin, chckLCSValid, initLCSPg, chckLCSDefine
    public :: initParaCsteCell, setMaterPara
    public :: getUserLCS, getUserLCSWithBaryCenter, initLCSZero, initLCSNone
    public :: compCellBary, getUserLCSCommand
    public :: lMaterVisc, copyMaterPara
! ==================================================================================================
    private
#include "asterc/r8dgrd.h"
#include "asterc/r8nnem.h"
#include "asterf_types.h"
#include "asterfort/angvx.h"
#include "asterfort/angvxy.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/ElasticityMaterial_type.h"
#include "asterfort/eiangl.h"
#include "asterfort/eulnau.h"
#include "asterfort/get_elas_id.h"
#include "asterfort/getvr8.h"
#include "asterfort/MaterialPara_type.h"
#include "asterfort/rccoma.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/utrcyl.h"
#include "jeveux.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! initParaCell
!
! Initializations of parameters on current cell
!
! In  fami             : integration point type
! In  jvMaterCode      : adress for material parameters
! IO  materPara        : parameters of material
!
! --------------------------------------------------------------------------------------------------
    subroutine initParaCell(famiZ, jvMaterCode, materPara)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: famiZ
        integer(kind=8), intent(in) :: jvMaterCode
        type(Material_Para), intent(inout) :: materPara
! ----- Local
        integer(kind=8) :: elasID, propCode
        character(len=8) :: fami
        character(len=16) :: elasKeyword
!   ------------------------------------------------------------------------------------------------
!
        fami = famiZ
        materPara%schemePara%fami = fami
        materPara%jvMaterCode = jvMaterCode
        call rccoma(jvMaterCode, 'ELAS', 0, elasKeyword, propCode)
        if (propCode .eq. 0) then
            call get_elas_id(jvMaterCode, elasID, elasKeyword)
            materPara%lElasIsMeta = (elasKeyword == 'ELAS_META')
            materPara%elasID = elasID
            materPara%elasKeyword = elasKeyword
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! setMaterPara
!
! Initializations of material parameters on current cell
!
! In  fami             : integration point type
! IO  materPara        : parameters of material
!
! --------------------------------------------------------------------------------------------------
    subroutine setMaterPara(fami, materPara)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: fami
        type(Material_Para), intent(inout) :: materPara
! ----- Local
        integer(kind=8) :: jvMaterc, jvMaterCode, iret
!   ------------------------------------------------------------------------------------------------
!
        call tecach('NNO', 'PMATERC', 'L', iret, iad=jvMaterc)
        if (iret .eq. 0) then
            jvMaterCode = zi(jvMaterc)
            call initParaCell(fami, jvMaterCode, materPara)
        else
            materPara%schemePara%fami = fami
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! initParaPoin
!
! Initializations of material parameters on current integration point
!
! In  kpg              : index of quadrature point
! In  ksp              : index of "sub"-point (plates, pipes, beams, etc.)
! IO  materPara        : parameters of material
!
! --------------------------------------------------------------------------------------------------
    subroutine initParaPoin(kpg, ksp, materPara)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: kpg, ksp
        type(Material_Para), intent(inout) :: materPara
!   ------------------------------------------------------------------------------------------------
!
        materPara%schemePara%kpg = kpg
        materPara%schemePara%ksp = ksp
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! initParaCsteCell
!
! Initializations of material parameters as constants on current cell
!
! IO  materPara        : parameters of material
!
! --------------------------------------------------------------------------------------------------
    subroutine initParaCsteCell(materPara)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(Material_Para), intent(inout) :: materPara
! ----- Locals
        integer(kind=8), parameter :: kpg = 1, ksp = 1
!   ------------------------------------------------------------------------------------------------
!
        materPara%schemePara%kpg = kpg
        materPara%schemePara%ksp = ksp
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! initLCSNone
!
! Set local coordinate system to None
!
! IO  materPara        : parameters of material
!
! --------------------------------------------------------------------------------------------------
    subroutine initLCSNone(materPara)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(Material_Para), intent(inout) :: materPara
!   ------------------------------------------------------------------------------------------------
!
        materPara%lcsPara%lcsType = MATER_LCS_NONE
        materPara%lcsPara%lcsAngle = r8nnem()
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! initLCSZero
!
! Set local coordinate system to Zero
!
! IO  materPara        : parameters of material
!
! --------------------------------------------------------------------------------------------------
    subroutine initLCSZero(materPara)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(Material_Para), intent(inout) :: materPara
!   ------------------------------------------------------------------------------------------------
!
        materPara%lcsPara%lcsType = MATER_LCS_ZERO
        materPara%lcsPara%lcsAngle = 0.d0
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! initLCSPg
!
! Set local coordinate system to Zero
!
! IO  materPara        : parameters of material
!
! --------------------------------------------------------------------------------------------------
    subroutine initLCSPg(ndim, nbNode, materPara)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: ndim, nbNode
        type(Material_Para), intent(inout) :: materPara
! ----- Locals
        integer(kind=8) :: jvCamass
        real(kind=8) :: anglNautPg(3*nbNode)
!   ------------------------------------------------------------------------------------------------
!
        if (.not. hasOrieField(jvCamass)) then
            call utmess('F', 'JOINT1_3')
        end if
        if (zr(jvCamass) .lt. 0.d0) then
            call utmess('F', 'JOINT1_47')
        end if
        call eiangl(ndim, nbNode, zr(jvCamass+1), anglNautPg)
        materPara%lcsPara%lcsType = MATER_LCS_PG
        materPara%lcsPara%lcsAngle = r8nnem()
        materPara%lcsPara%lcsAnglePg = 0.d0
        materPara%lcsPara%lcsAnglePg(1:3*nbNode) = anglNautPg
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! chckLCSValid
!
! In  materPara        : parameters of material
!
! --------------------------------------------------------------------------------------------------
    subroutine chckLCSValid(materPara)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(Material_Para), intent(in) :: materPara
!   ------------------------------------------------------------------------------------------------
!
        if (materPara%lcsPara%lcsType .eq. MATER_LCS_NONE) then
            call utmess('F', 'ALGORITH8_20')
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getUserLCS
!
! Get local coordinate system from user
!
! In  ndim             : dimension of element (2 ou 3)
! In  nbNode           : number of nodes
! In  jvGeom           : JEVEUX adress to initial geometry (mesh)
! Out lcsPara          : parameters of local coordinate system
!
! --------------------------------------------------------------------------------------------------
    subroutine getUserLCS(ndim, nbNode, jvGeom, lcsPara)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: ndim, nbNode, jvGeom
        type(LCS_Para), intent(out) :: lcsPara
! ----- Locals
        integer(kind=8) :: jvCamass, iDim
        real(kind=8) :: coorBary(3), anglNaut(3)
        real(kind=8) :: p(3, 3), xg(3), yg(3), orig(3), dire(3)
        real(kind=8) :: alpha, beta, xu, yu, xnorm
!   ------------------------------------------------------------------------------------------------
!
        call compCellBary(ndim, nbNode, jvGeom, coorBary)
        anglNaut = 0.d0
        lcsPara%lcsType = MATER_LCS_ZERO
        if (hasOrieField(jvCamass)) then
            if (zr(jvCamass) .gt. 0.d0) then
                lcsPara%lcsType = MATER_LCS_NAUT
                anglNaut(1) = zr(jvCamass+1)*r8dgrd()
                if (ndim .eq. 3) then
                    anglNaut(2) = zr(jvCamass+2)*r8dgrd()
                    anglNaut(3) = zr(jvCamass+3)*r8dgrd()
                end if

            else if (abs(zr(jvCamass)+1.d0) .lt. 1.d-3) then
                lcsPara%lcsType = MATER_LCS_CYL
! ON TRANSFORME LA DONNEE DU REPERE CYLINDRIQUE EN ANGLE NAUTIQUE
                orig(1:ndim) = zr(jvCamass+3+1:jvCamass+3+ndim)
                if (ndim .eq. 3) then
                    alpha = zr(jvCamass+1)*r8dgrd()
                    beta = zr(jvCamass+2)*r8dgrd()
                    dire(1) = cos(alpha)*cos(beta)
                    dire(2) = sin(alpha)*cos(beta)
                    dire(3) = -sin(beta)
                    call utrcyl(coorBary, dire, orig, p)
                    do iDim = 1, 3
                        xg(iDim) = p(1, iDim)
                        yg(iDim) = p(2, iDim)
                    end do
                    call angvxy(xg, yg, anglNaut)
                else
                    xu = coorBary(1)-orig(1)
                    yu = coorBary(2)-orig(2)
                    xnorm = sqrt(xu**2+yu**2)
                    xu = xu/xnorm
                    yu = yu/xnorm
                    p(1, 1) = xu
                    p(2, 1) = yu
                    p(1, 2) = -yu
                    p(2, 2) = xu
                    xg(1) = xu
                    xg(2) = yu
                    xg(3) = 0.d0
                    call angvx(xg, alpha, beta)
                    anglNaut(1) = alpha
                end if
            end if
        end if

! ----- Set parameters
        lcsPara%lcsAngle = anglNaut
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compCellBary
!
! Compute barycenter of cell from coordinates of nodes
!
! In  ndim             : dimension of element (2 ou 3)
! In  nbNode           : number of nodes
! In  jvGeom           : JEVEUX adress to initial geometry (mesh)
! Out coorBary         : coordinates of barycenter of cell
!
! --------------------------------------------------------------------------------------------------
    subroutine compCellBary(ndim, nbNode, jvGeom, coorBary)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: ndim, nbNode, jvGeom
        real(kind=8), intent(out) :: coorBary(3)
! ----- Locals
        integer(kind=8) :: iNode, iDim
!   ------------------------------------------------------------------------------------------------
!
        coorBary = 0.d0
        do iNode = 1, nbNode
            do iDim = 1, ndim
                coorBary(iDim) = coorBary(iDim)+ &
                                 zr(jvGeom+iDim+ndim*(iNode-1)-1)/nbNode
            end do
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getUserLCSWithBaryCenter
!
! Get local coordinate system from user and given barycenter of cell
!
! In  ndim             : dimension of element (2 ou 3)
! In  coorBary         : coordinates of barycenter of cell
! Out lcsPara          : parameters of local coordinate system
!
! --------------------------------------------------------------------------------------------------
    subroutine getUserLCSWithBaryCenter(ndim, coorBary, lcsPara)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: ndim
        real(kind=8), intent(in) :: coorBary(3)
        type(LCS_Para), intent(out) :: lcsPara
! ----- Locals
        integer(kind=8) :: jvCamass, iDim
        real(kind=8) :: anglNaut(3)
        real(kind=8) :: p(3, 3), xg(3), yg(3), orig(3), dire(3)
        real(kind=8) :: alpha, beta, xu, yu, xnorm
!   ------------------------------------------------------------------------------------------------
!
        anglNaut = 0.d0
        if (hasOrieField(jvCamass)) then
            if (zr(jvCamass) .gt. 0.d0) then
                lcsPara%lcsType = MATER_LCS_NAUT
                anglNaut(1) = zr(jvCamass+1)*r8dgrd()
                if (ndim .eq. 3) then
                    anglNaut(2) = zr(jvCamass+2)*r8dgrd()
                    anglNaut(3) = zr(jvCamass+3)*r8dgrd()
                end if

            else if (abs(zr(jvCamass)+1.d0) .lt. 1.d-3) then
                lcsPara%lcsType = MATER_LCS_CYL
! ON TRANSFORME LA DONNEE DU REPERE CYLINDRIQUE EN ANGLE NAUTIQUE
                orig(1:ndim) = zr(jvCamass+3+1:jvCamass+3+ndim)
                if (ndim .eq. 3) then
                    alpha = zr(jvCamass+1)*r8dgrd()
                    beta = zr(jvCamass+2)*r8dgrd()
                    dire(1) = cos(alpha)*cos(beta)
                    dire(2) = sin(alpha)*cos(beta)
                    dire(3) = -sin(beta)
                    call utrcyl(coorBary, dire, orig, p)
                    do iDim = 1, 3
                        xg(iDim) = p(1, iDim)
                        yg(iDim) = p(2, iDim)
                    end do
                    call angvxy(xg, yg, anglNaut)
                else
                    xu = coorBary(1)-orig(1)
                    yu = coorBary(2)-orig(2)
                    xnorm = sqrt(xu**2+yu**2)
                    xu = xu/xnorm
                    yu = yu/xnorm
                    p(1, 1) = xu
                    p(2, 1) = yu
                    p(1, 2) = -yu
                    p(2, 2) = xu
                    xg(1) = xu
                    xg(2) = yu
                    xg(3) = 0.d0
                    call angvx(xg, alpha, beta)
                    anglNaut(1) = alpha
                end if
            end if
        end if

! ----- Set parameters
        lcsPara%lcsAngle = anglNaut
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! lMaterVisc
!
! Detect viscoelasticity
!
! In  materPara        : parameters of material
!
! --------------------------------------------------------------------------------------------------
    aster_logical function lMaterVisc(materPara)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(Material_Para), intent(in) :: materPara
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(materPara%elasID .ne. ELAS_UNDEF)
        lMaterVisc = materPara%elasID .eq. ELAS_VISC_ISOT .or. &
                     materPara%elasID .eq. ELAS_VISC_ISTR .or. &
                     materPara%elasID .eq. ELAS_VISC_ORTH
!
!   ------------------------------------------------------------------------------------------------
    end function
! --------------------------------------------------------------------------------------------------
!
! chckLCSDefine
!
! Detect if anisotropic local coordinate system has been defined
!
! In  materPara        : parameters of material
!
! --------------------------------------------------------------------------------------------------
    aster_logical function chckLCSDefine(lcsPara)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(LCS_Para), intent(in) :: lcsPara
!   ------------------------------------------------------------------------------------------------
!
        chckLCSDefine = lcsPara%lcsType .eq. MATER_LCS_NAUT .or. &
                        lcsPara%lcsType .eq. MATER_LCS_CYL .or. &
                        lcsPara%lcsType .eq. MATER_LCS_EULER .or. &
                        lcsPara%lcsType .eq. MATER_LCS_ZERO
!
!   ------------------------------------------------------------------------------------------------
    end function
! --------------------------------------------------------------------------------------------------
!
! getUserLCSCommand
!
! Get local coordinate system from user (command)
!
! IO  lcsPara          : parameters of local coordinate system
!
! --------------------------------------------------------------------------------------------------
    subroutine getUserLCSCommand(ndim, lcsPara)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: ndim
        type(LCS_Para), intent(inout) :: lcsPara
! ----- Locals
        real(kind=8) :: anglNaut(3), anglEuler(3), angd(3)
        integer(kind=8) :: n1, n2
!   ------------------------------------------------------------------------------------------------
!
        anglNaut = 0.d0
        anglEuler = 0.d0
        call getvr8('MASSIF', 'ANGL_REP', iocc=1, nbval=3, vect=anglNaut, nbret=n1)
        call getvr8('MASSIF', 'ANGL_EULER', iocc=1, nbval=3, vect=anglEuler, nbret=n2)
        if (n1 .gt. 0) then
            anglNaut(1) = anglNaut(1)*r8dgrd()
            if (ndim .eq. 3) then
                anglNaut(2) = anglNaut(2)*r8dgrd()
                anglNaut(3) = anglNaut(3)*r8dgrd()
            end if
            lcsPara%lcsType = MATER_LCS_NAUT
            lcsPara%lcsAngle = anglNaut
        else if (n2 .gt. 0) then
            call eulnau(anglEuler, angd)
            anglNaut(1) = angd(1)*r8dgrd()
            if (ndim .eq. 3) then
                anglNaut(2) = angd(2)*r8dgrd()
                anglNaut(3) = angd(3)*r8dgrd()
            end if
            lcsPara%lcsType = MATER_LCS_EULER
            lcsPara%lcsAngle = anglNaut
        else
            lcsPara%lcsType = MATER_LCS_ZERO
            lcsPara%lcsAngle = 0.d0

        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! copyMaterPara
!
! Copy material parameters
!
! --------------------------------------------------------------------------------------------------
    subroutine copyMaterPara(materParaIn, fami, kpg, ksp, &
                             materParaOut)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(Material_Para), intent(in) :: materParaIn
        character(len=*), intent(in) :: fami
        integer(kind=8), intent(in) :: kpg, ksp
        type(Material_Para), intent(out) :: materParaOut
!   ------------------------------------------------------------------------------------------------
!
        materParaOut%elasID = materParaIn%elasID
        materParaOut%elasKeyword = materParaIn%elasKeyword
        materParaOut%jvMaterCode = materParaIn%jvMaterCode
        materParaOut%lElasIsMeta = materParaIn%lElasIsMeta
        materParaOut%lMetaLemaAni = materParaIn%lMetaLemaAni
        materParaOut%lcsPara = materParaIn%lcsPara
        materParaOut%schemePara = materParaIn%schemePara
        materParaOut%schemePara%fami = fami
        materParaOut%schemePara%kpg = kpg
        materParaOut%schemePara%ksp = ksp
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module MaterialPara_module
