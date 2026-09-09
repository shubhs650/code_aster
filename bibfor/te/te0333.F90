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
subroutine te0333(option, nomte)
!
    use BehaviourStrain_module
    use BehaviourStrain_type
    use FE_basis_module
    use FE_eval_module
    use FE_quadrature_module
    use FE_topo_module
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterc/r8nnem.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/calcgr.h"
#include "asterfort/ElasticityMaterial_type.h"
#include "asterfort/get_elas_para.h"
#include "asterfort/granvi.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/nbsigm.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "FE_module.h"
#include "jeveux.h"
#include "MeshTypes_type.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: 3D
! Option: EPSP_ELGA
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = "RIGI"
    integer(kind=8), parameter :: ksp = 1, mxcmel = 162
    integer(kind=8), parameter :: nbsgm = 6
    real(kind=8) :: epsiPlas(mxcmel), epsiCreep(nbsgm)
    real(kind=8) :: epsiTota(6), epsiVarc(6), epsiMeca(6), sigmEner(6)
    integer(kind=8) :: nbVari, variIndxTemp, nbVariGranger
    real(kind=8) :: e, nu, c1, c2, trsig
    aster_logical :: l_creep, lTempInVari, lMetaLemaAni, lTHM
    integer(kind=8) :: jvMaterc, jvDisp, jvVari, jvSigm, jvEpsi, jvTime
    integer(kind=8) :: kpg, npg, ndim, nno, iSig, iEps, nbSig, nbEps
    integer(kind=8) :: jtab(7), iret
    type(FE_Cell) :: FECell
    type(FE_Quadrature) :: FEQuad
    type(FE_basis) :: FEBasis
    real(kind=8) :: BGSEval(3, MAX_BS_CG), coorpg(3)
    character(len=16), pointer :: compor(:) => null()
    character(len=16) :: relaComp, comporFlua, comporPlas
    type(All_Varc_Strain) :: allVarcStrain
    real(kind=8) :: tempkpg
    type(Material_Para) :: materPara
!
! --------------------------------------------------------------------------------------------------
!
    lTHM = lteatt('TYPMOD2', 'THM')

! - Initialize a FE cell
    call FECell%init()
    ndim = FECell%ndim
    nno = FECell%nbnodes
    ASSERT(nno .le. MT_NNOMAX3D)

! - Initialization of quadrature
    call FEQuad%initCell(FECell, fami)
    npg = FEQuad%nbQuadPoints
    ASSERT(npg .le. MT_NNOMAX3D)

! - Initialization of basis functions
    call FEBasis%initCell(FECell)

! - Current displacements (nodes)
    call jevech('PDEPLAR', 'L', jvDisp)

! - Internal variables
    call jevech('PVARIGR', 'L', jvVari)
    call tecach('OOO', 'PVARIGR', 'L', iret, nval=7, itab=jtab)
    nbVari = max(jtab(6), 1)*jtab(7)

! - Stresses
    call jevech('PCONTRR', 'L', jvSigm)
    nbSig = nbsigm()
    nbEps = nbSig
    ASSERT(nbSig .eq. nbsgm)

! - Current time
    call jevech('PINSTR', 'L', jvTime)
    allVarcStrain%time = zr(jvTime)
    allVarcStrain%hasTime = ASTER_TRUE

! - Get material parameters
    call jevech('PMATERC', 'L', jvMaterc)
    call initParaCell(fami, zi(jvMaterc), materPara)
    call initLCSNone(materPara)
    if (materPara%elasID .ne. ELAS_ISOT) then
        call utmess('F', 'ELEMENTS6_2')
    end if

! - Get fields for non-linear behaviour
    call jevech('PCOMPOR', 'L', vk16=compor)

! - Properties of behaviour
    relaComp = compor(RELA_NAME)
    comporFlua = compor(CREEP_NAME)
    comporPlas = compor(PLAS_NAME)
    lMetaLemaAni = (relaComp == "META_LEMA_ANI")

! - Detect Granger law (creep)
    if (relaComp(1:13) .ne. 'BETON_GRANGER' .and. &
        (relaComp .ne. 'KIT_DDI' .or. comporFlua(1:13) .ne. 'BETON_GRANGER')) then
        l_creep = ASTER_FALSE
    else
        call granvi("3D", nvi_=nbVariGranger)
        l_creep = ASTER_TRUE
    end if

! - Has temperature in internal state variable (== maximum)
    lTempInVari = ASTER_FALSE
    if (relaComp .eq. 'BETON_DOUBLE_DP') then
        variIndxTemp = 3
        lTempInVari = ASTER_TRUE
    else if (relaComp .eq. 'KIT_DDI') then
        if (comporPlas .eq. 'BETON_DOUBLE_DP') then
            if (comporFlua(1:13) .eq. 'BETON_GRANGER') then
                variIndxTemp = nbVariGranger+3
                lTempInVari = ASTER_TRUE
            else
                call utmess('F', 'COMPOR5_76')
            end if
        end if
    end if

! - Loop on Gauss points
    do kpg = 1, npg
! ----- Current coordinates of Gauss point
        coorpg = FEQuad%points_param(1:3, kpg)

! ----- Compute the gradient of the scalar basis
        BGSEval = FEBasis%grad(coorpg, FEQuad%jacob(1:3, 1:3, kpg))

! ----- Kinematic - Total strains
        epsiTota = FEEvalGradSymMat(FEBasis, zr(jvDisp), coorpg, BGSEval)
        epsiTota(4:6) = epsiTota(4:6)/sqrt(2.d0)

! ----- Initializations of material parameters on current integration point
        call initParaPoin(kpg, ksp, materPara)

! ----- Detect external state variable
        call strainDetectVarc('+', lTHM, materPara, allVarcStrain)

! ----- Get current temperature
        tempkpg = r8nnem()
        if (allVarcStrain%list(VARC_STRAIN_TEMP)%exist) then
            tempkpg = allVarcStrain%list(VARC_STRAIN_TEMP)%varcCurr(1)
        end if

! ----- Change temperature from internal variable (maximum) for BETON_DOUBLE_DP/BETON_GRANGER
        if (lTempInVari) then
            if (tempkpg .lt. zr(jvVari+(kpg-1)*nbVari+variIndxTemp-1)) then
                tempkpg = zr(jvVari+(kpg-1)*nbVari+variIndxTemp-1)
            end if
        end if

! ----- Set temperature
        allVarcStrain%hasTemp = ASTER_TRUE
        allVarcStrain%temp = tempkpg

! ----- Get elastic parameters (only isotropic elasticity)
        call get_elas_para(fami, zi(jvMaterc), '+', kpg, ksp, &
                           materPara%elasID, materPara%elasKeyword, &
                           time=allVarcStrain%time, temp=allVarcStrain%temp, e_=e, nu_=nu)

! ----- Compute non-mechanical strains (epsiVarc) for some external state variables
        call compVarcStrain('+', materPara, allVarcStrain)
        call getVarcStrain('+', VARC_STRAIN_ALL, allVarcStrain, epsiVarc)
        epsiVarc(4:6) = epsiVarc(4:6)/sqrt(2.d0)

! ----- Compute mechanical strains epsiMeca = epsiTota - epsiVarc
        epsiMeca = 0.d0
        epsiMeca = epsiTota-epsiVarc

! ----- Compute creep strains
        epsiCreep = 0.d0
        if (l_creep) then
            call calcgr(kpg, nbSig, nbVari, zr(jvVari), nu, epsiCreep)
        end if

! ----- Compute stresses
        do iSig = 1, nbSig
            sigmEner(iSig) = zr(jvSigm+(kpg-1)*nbSig+iSig-1)
        end do
        trsig = sigmEner(1)+sigmEner(2)+sigmEner(3)

! ----- Compute plastic strains epsiPlas = epsi_tota - epsi_elas - epsiCreep
        c1 = (1.d0+nu)/e
        c2 = nu/e
        epsiPlas(nbEps*(kpg-1)+1) = epsiMeca(1)-(c1*sigmEner(1)-c2*trsig)-epsiCreep(1)
        epsiPlas(nbEps*(kpg-1)+2) = epsiMeca(2)-(c1*sigmEner(2)-c2*trsig)-epsiCreep(2)
        epsiPlas(nbEps*(kpg-1)+3) = epsiMeca(3)-(c1*sigmEner(3)-c2*trsig)-epsiCreep(3)
        epsiPlas(nbEps*(kpg-1)+4) = epsiMeca(4)-c1*sigmEner(4)-epsiCreep(4)
        epsiPlas(nbEps*(kpg-1)+5) = epsiMeca(5)-c1*sigmEner(5)-epsiCreep(5)
        epsiPlas(nbEps*(kpg-1)+6) = epsiMeca(6)-c1*sigmEner(6)-epsiCreep(6)
    end do

! - Plastic strain output
    call jevech('PDEFOPG', 'E', jvEpsi)
    do kpg = 1, npg
        do iEps = 1, nbEps
            zr(jvEpsi+nbEps*(kpg-1)+iEps-1) = epsiPlas(nbEps*(kpg-1)+iEps)
        end do
    end do
!
end subroutine
