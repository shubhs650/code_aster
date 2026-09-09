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
! aslint: disable=W1306
!
! ==================================================================================================
!
! Module for the management of integration of behaviour
!
! ==================================================================================================
!
module Behaviour_module
! ==================================================================================================
    use Behaviour_type
    use BehaviourMGIS_type
    use BehaviourStrain_module
    use BehaviourStrain_type
    use HHO_quadrature_module
    use HHO_type
    use MaterialPara_module
    use MaterialPara_type
    use calcul_module, only: ca_jvcnom_, ca_nbcvrc_
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: behaviourOption, behaviourInit, behaviourInitPoint
    private :: initPara, initESVA
    public :: behaviourSetParaCell, behaviourSetParaPoin, setFromCompor
    public :: behaviourPrepESVAPoin, behaviourPrepStrain, behaviourPrepESVAExte
    public :: behaviourPrepESVAGeom, behaviourPrepModel, behaviourPredictionStress
    public :: behaviourCoorGauss, behaviourPrepESVAGeomHHO
    private :: computeStrainESVA, computeStrainMeca
    private :: varcIsGEOM, isSolverIsExte
    private :: prepEltSize1, prepGradVelo
    public :: getAsterVariableName, getMFrontVariableName
    private :: setFromOption, setFromCarcri, chckBounds
    public :: detectVarc, behaviourGetParameters
! ==================================================================================================
    private
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/indik8.h"
#include "asterc/mgis_get_esvs.h"
#include "asterc/mgis_get_number_of_esvs.h"
#include "asterc/r8nnem.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/BehaviourMGIS_type.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/ElasticityMaterial_type.h"
#include "asterfort/isdeco.h"
#include "asterfort/jevech.h"
#include "asterfort/leverettIsotMeca.h"
#include "asterfort/matinv.h"
#include "asterfort/nmgeom.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
#include "MeshTypes_type.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! detectVarc
!
! Detect external state variables
!
! IO  BEHInteg         : main object for managing the integration of behavior laws
!
! --------------------------------------------------------------------------------------------------
    subroutine detectVarc(BEHInteg)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(Behaviour_Integ), intent(inout) :: BEHInteg
! ----- Local
        character(len=1), parameter :: poum = "T"
!   ------------------------------------------------------------------------------------------------
!
        call strainDetectVarc(poum, BEHInteg%behavPara%lTHM, BEHInteg%materPara, &
                              BEHInteg%allVarcStrain)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! behaviourInit
!
! Initialisation of behaviour datastructure
!
! Out BEHInteg         : main object for managing the integration of behavior laws
!
! --------------------------------------------------------------------------------------------------
    subroutine behaviourInit(BEHInteg)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(Behaviour_Integ), intent(out) :: BEHInteg
!   ------------------------------------------------------------------------------------------------
!
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG> Initialization of datastructures'
        end if

! ----- Initialization of parameters
        call initPara(BEHInteg%behavPara)

! ----- Initialization of parameters for external state variables
        call initESVA(BEHInteg%behavESVA)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! initPara
!
! Initialisation of parameters
!
! Out BEHInteg         : main object for managing the integration of behavior laws
!
! --------------------------------------------------------------------------------------------------
    subroutine initPara(behavPara)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(Behaviour_Para), intent(out) :: behavPara
!   ------------------------------------------------------------------------------------------------
!
        behavPara%timeCurr = r8nnem()
        behavPara%timePrev = r8nnem()
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! initESVA
!
! Initialization of parameters for external state variables
!
! Out behavESVA        : parameters for external state variables
!
! --------------------------------------------------------------------------------------------------
    subroutine initESVA(behavESVA)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(BehaviourESVA), intent(out) :: behavESVA
! ----- Locam
        aster_logical :: lGeomInESVA
!   ------------------------------------------------------------------------------------------------
!

! ----- Detect if geometry is in external state variables
        call varcIsGEOM(lGeomInESVA)
        behavESVA%lGeomInESVA = lGeomInESVA

! ----- For geometric properties
        behavESVA%behavESVAGeom%elemSize1 = r8nnem()
        behavESVA%behavESVAGeom%gradVelo(9) = r8nnem()
        behavESVA%behavESVAGeom%coorElga(27, 3) = r8nnem()

! ----- For other properties
        behavESVA%behavESVAOther%rotpg(3*3) = r8nnem()
        behavESVA%behavESVAOther%tenscab = r8nnem()
        behavESVA%behavESVAOther%curvcab = r8nnem()
        behavESVA%behavESVAOther%nonloc(2) = r8nnem()
        behavESVA%behavESVAOther%r = r8nnem()
        behavESVA%behavESVAOther%time = r8nnem()
        behavESVA%behavESVAOther%hygrPrev = r8nnem()
        behavESVA%behavESVAOther%hygrIncr = r8nnem()
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! behaviourSetParaCell
!
! Set parameters on a cell
!
! In  typmod           : finite element model
! In  option           : option to compute
! In  compor           : map for behaviour
! In  carcri           : parameters for comportment
! In  timePrev         : time at beginning of time step
! In  timeCurr         : time at end of time step
! In  materPara        : parameters of material
! IO  BEHInteg         : main object for managing the integration of behavior laws
!
! --------------------------------------------------------------------------------------------------
    subroutine behaviourSetParaCell(typmod, option, &
                                    compor, carcri, &
                                    timePrev, timeCurr, &
                                    materPara, BEHInteg)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=8), intent(in) :: typmod(2)
        character(len=16), intent(in) :: option
        character(len=16), intent(in) :: compor(COMPOR_SIZE)
        real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
        real(kind=8), intent(in) :: timePrev, timeCurr
        type(Material_Para), intent(in) :: materPara
        type(Behaviour_Integ), intent(inout) :: BEHInteg
!   ------------------------------------------------------------------------------------------------
!
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG>  Paramètres constants sur la cellule'
        end if

! ----- Set general parameters on cell
        BEHInteg%materPara = materPara
        BEHInteg%behavPara%timePrev = timePrev
        BEHInteg%behavPara%timeCurr = timeCurr

! ----- Set parameters for material properties
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG>  Récupération de l élasticité'
        end if
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG>  Type d élasticité      : ', BEHInteg%materPara%elasID
            WRITE (6, *) '<DEBUG>  Mot-clef               : ', BEHInteg%materPara%elasKeyword
            WRITE (6, *) '<DEBUG>  Présence de métallurgie: ', BEHInteg%materPara%lElasIsMeta
        end if

! ----- Detect VERI_BORNE
        call chckBounds(BEHInteg)

! ----- Set parameters from option
        call setFromOption(option, BEHInteg)

! ----- Set parameters from COMPOR map
        call setFromCompor(compor, BEHInteg)

! ----- Set parameters from CARCRI map
        if (BEHInteg%behavPara%lNonLinear) then
            call setFromCarcri(carcri, BEHInteg)
        end if

! ----- Set finite element model
        call behaviourPrepModel(typmod, BEHInteg)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! setFromOption
!
! In  option           : option to compute
! IO  BEHInteg         : main object for managing the integration of behavior laws
!
! --------------------------------------------------------------------------------------------------
    subroutine setFromOption(option, BEHInteg)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=16), intent(in) :: option
        type(Behaviour_Integ), intent(inout) :: BEHInteg
!   ------------------------------------------------------------------------------------------------
!
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG>  From OPTION'
        end if

        BEHInteg%behavPara%lImplex = option .eq. "RIGI_MECA_IMPLEX" .or. &
                                     option .eq. "RAPH_MECA_IMPLEX"
        BEHInteg%behavPara%lVari = L_VARI(option)
        BEHInteg%behavPara%lSigm = L_SIGM(option)
        BEHInteg%behavPara%lMatr = L_MATR(option)
        BEHInteg%behavPara%lPred = L_PRED(option)
        BEHInteg%behavPara%lNonLinear = option .ne. "FORC_NODA"
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! setFromCompor
!
! In  compor           : map for behaviour
! IO  BEHInteg         : main object for managing the integration of behavior laws
!
! --------------------------------------------------------------------------------------------------
    subroutine setFromCompor(compor, BEHInteg)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=16), intent(in) :: compor(COMPOR_SIZE)
        type(Behaviour_Integ), intent(inout) :: BEHInteg
! ----- Locals
        character(len=16) :: relaComp, defoLDC, defoComp, reguVisc, postIncr, mgisAddr
        integer(kind=8) :: nvi, numlc
        aster_logical :: lDeborst
!   ------------------------------------------------------------------------------------------------
!
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG>  From COMPOR'
        end if

! ----- Get parameters from behaviour
        defoLDC = compor(DEFO_LDC)
        relaComp = compor(RELA_NAME)
        defoComp = compor(DEFO)
        reguVisc = compor(REGUVISC)
        postIncr = compor(POSTINCR)
        mgisAddr = compor(MGIS_ADDR)
        read (compor(NUME), '(I16)') numlc
        read (compor(NVAR), '(I16)') nvi

! ----- Set flags from behaviour
        BEHInteg%behavPara%lFiniteStrain = defoComp .eq. 'SIMO_MIEHE' .or. &
                                           defoComp .eq. 'GROT_GDEP' .or. &
                                           defoComp .eq. 'GREEN_LAGRANGE'
        BEHInteg%behavPara%lGdefLog = defoComp .eq. 'GDEF_LOG'
        BEHInteg%behavPara%lAnnealing = postIncr .eq. "REST_ECRO"
        BEHInteg%behavPara%lStrainMeca = defoLDC .eq. 'MECANIQUE'
        BEHInteg%behavPara%lStrainAll = defoLDC .eq. 'TOTALE'
        BEHInteg%behavPara%lStrainOld = defoLDC .eq. 'OLD'
        BEHInteg%behavPara%lReguVisc = reguVisc .eq. 'REGU_VISC_ELAS'
        BEHInteg%materPara%lMetaLemaAni = relaComp .eq. 'META_LEMA_ANI'
        BEHInteg%behavESVA%behavESVAExte%mgisAddr = mgisAddr
        BEHInteg%behavPara%nvi = nvi
        BEHInteg%behavPara%numlc = numlc
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG>  From COMPOR - lFiniteStrain: ', BEHInteg%behavPara%lFiniteStrain
            WRITE (6, *) '<DEBUG>  From COMPOR - lGdefLog: ', BEHInteg%behavPara%lGdefLog
            WRITE (6, *) '<DEBUG>  From COMPOR - lAnnealing: ', BEHInteg%behavPara%lAnnealing
            WRITE (6, *) '<DEBUG>  From COMPOR - lStrainMeca: ', BEHInteg%behavPara%lStrainMeca
            WRITE (6, *) '<DEBUG>  From COMPOR - lStrainAll: ', BEHInteg%behavPara%lStrainAll
            WRITE (6, *) '<DEBUG>  From COMPOR - lStrainOld: ', BEHInteg%behavPara%lStrainOld
            WRITE (6, *) '<DEBUG>  From COMPOR - lReguVisc: ', BEHInteg%behavPara%lReguVisc
            WRITE (6, *) '<DEBUG>  From COMPOR - mgisAddr: ', mgisAddr
            if (lDeborst) then
                WRITE (6, *) '<DEBUG>  From COMPOR - De Borst algorithm'
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! setFromCarcri
!
! In  carcri           : parameters for comportment
! IO  BEHInteg         : main object for managing the integration of behavior laws
!
! --------------------------------------------------------------------------------------------------
    subroutine setFromCarcri(carcri, BEHInteg)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
        type(Behaviour_Integ), intent(inout) :: BEHInteg
! ----- Locals
        aster_logical :: lMGIS, lUMAT
!   ------------------------------------------------------------------------------------------------
!
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG>  From CARCRI'
        end if

! ----- External solvers ?
        call isSolverIsExte(carcri, lMGIS, lUMAT)
        BEHInteg%behavPara%lMGIS = lMGIS
        BEHInteg%behavPara%lUMAT = lUMAT
        BEHInteg%behavPara%lExteSolver = lUMAT .or. lMGIS

! ----- DEBUG
        if (LDC_PREP_DEBUG .eq. 1) then
            if (lMGIS) then
                WRITE (6, *) '<DEBUG>  External solver: MFront'
            elseif (lUMAT) then
                WRITE (6, *) '<DEBUG>  External solver: UMAT'
            else
                WRITE (6, *) '<DEBUG>  Internal solver'
            end if
        end if

! ----- Get list of external state variables from user
        call getListUserESVA(carcri, BEHInteg%behavESVA%tabcod)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! behaviourPrepModel
!
! In  typmod           : finite element model
! IO  BEHInteg         : main object for managing the integration of behavior laws
!
! --------------------------------------------------------------------------------------------------
    subroutine behaviourPrepModel(typmod, BEHInteg)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=8), intent(in) :: typmod(2)
        type(Behaviour_Integ), intent(inout) :: BEHInteg
! ----- Locals
        character(len=8), parameter:: &
            allowedTypmod1(*) = [character(len=8):: '1D', 'PLAN', 'C_PLAN', 'D_PLAN', 'AXIS', '3D']
        character(len=8), parameter:: &
            allowedTypmod2(*) = [character(len=8):: ' ', 'GRADVARI', 'HHO', 'THM', 'GRADSIGM', &
                                                                         'EJ_HYME', 'ELEMJOIN', &
                                                                         'INTERFAC', 'INTSOLPI', &
                                                                         'CABLE_GA', 'JHMS']
        character(len=8), parameter:: featNdimsiEnabled(*) = &
                                      [character(len=8):: ' ', 'GRADVARI', 'HHO', 'THM']
        character(len=8), parameter:: featEpsComEnabled(*) = &
                                      [character(len=8):: ' ', 'GRADVARI', 'HHO', 'THM']
        character(len=8), parameter:: featDeborstEnabled(*) = [character(len=8):: ' ', 'HHO']
        character(len=8), parameter:: featCtrlEpsEnabled(*) = [character(len=8):: ' ', 'HHO', 'THM']
        character(len=8), parameter:: featReguViscEnabled(*) = &
                                      [character(len=8):: ' ', 'GRADVARI', 'HHO']
!   ------------------------------------------------------------------------------------------------
        integer(kind=8) :: lawOffset
!   ------------------------------------------------------------------------------------------------
!
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG>  Paramètres du modèle'
        end if

        ! Exclusions
        ! if (.not. any(typmod(1) == allowedTypmod1)) write(6,*) 'typmod1 =',typmod(1)
        ! if (.not. any(typmod(2) == allowedTypmod2)) write(6,*) 'typmod2 =',typmod(2)
        ASSERT(any(typmod(1) == allowedTypmod1))
        ASSERT(any(typmod(2) == allowedTypmod2))

        ! Type of element
        BEHInteg%behavPara%lTHM = typmod(2) .eq. 'THM'
        BEHInteg%behavPara%lGradVari = typmod(2) .eq. 'GRADVARI' .or. typmod(2) .eq. 'HHO_GRAD'

        BEHInteg%behavPara%lAxis = typmod(1) .eq. 'AXIS'
        BEHInteg%behavPara%lPlaneStrain = typmod(1) (1:6) .eq. 'D_PLAN'
        BEHInteg%behavPara%lPlaneStress = typmod(1) (1:6) .eq. 'C_PLAN'
        BEHInteg%BehavPara%lTwoDim = any(typmod(1) == ['PLAN  ', 'C_PLAN', 'D_PLAN', 'AXIS  '])
        BEHInteg%BehavPara%lOneDim = typmod(1) .eq. '1D'
        BEHInteg%behavPara%lThreeDim = typmod(1) (1:2) .eq. '3D'

        ! Allowed features
        BEHInteg%behavPara%lNdimsiEnabled = any(typmod(2) == featNdimsiEnabled)
        BEHInteg%behavPara%lEpsComEnabled = any(typmod(2) == featEpsComEnabled)
        BEHInteg%behavPara%lDeborstEnabled = any(typmod(2) == featDeborstEnabled)
        BEHInteg%behavPara%lCtrlEpsEnabled = any(typmod(2) == featCtrlEpsEnabled)
        BEHInteg%behavPara%lReguViscEnabled = any(typmod(2) == featReguViscEnabled)

        ASSERT(.not. BEHInteg%behavPara%lDeborstEnabled .or. BEHInteg%behavPara%lNdimsiEnabled)
        ASSERT(.not. BEHInteg%behavPara%lEpsComEnabled .or. BEHInteg%behavPara%lNdimsiEnabled)
        ASSERT(.not. BEHInteg%behavPara%lCtrlEpsEnabled .or. BEHInteg%behavPara%lNdimsiEnabled)
        ASSERT(.not. BEHInteg%behavPara%lReguViscEnabled .or. BEHInteg%behavPara%lNdimsiEnabled)

        ASSERT(.not. BEHInteg%behavPara%lStrainMeca .or. BEHInteg%behavPara%lEpsComEnabled)
        ASSERT(.not. BEHInteg%behavPara%lTHM .or. BEHInteg%behavPara%lEpsComEnabled)
        ASSERT(.not. BEHInteg%behavPara%lReguVisc .or. BEHInteg%behavPara%lReguViscEnabled)

        if (BEHInteg%behavPara%lFiniteStrain) then
            ASSERT(.not. BEHInteg%behavPara%lTHM)
            ASSERT(.not. BEHInteg%behavPara%lReguVisc)
            ASSERT(.not. BEHInteg%behavPara%lStrainMeca .or. ca_nbcvrc_ .eq. 0)
        end if

        ! Default size of strain and stress tensors
        if (BEHInteg%behavPara%lNdimsiEnabled) then
            if (BEHInteg%behavPara%lThreeDim) then
                BEHInteg%behavPara%ndimsi = 6
            else if (BEHInteg%behavPara%lTwoDim) then
                BEHInteg%behavPara%ndimsi = 4
            else if (BEHInteg%BehavPara%lOneDim) then
                BEHInteg%behavPara%ndimsi = 1
            else
                write (6, *) 'PROBLEME: typmod  -> ', typmod(1:2)
                ASSERT(ASTER_FALSE)
            end if
        else
            BEHInteg%behavPara%ndimsi = 0
        end if

        lawOffset = 0
        if (BEHInteg%behavPara%lImplex) then
            lawOffset = lawOffset+2000
        end if
        if (typmod(2) .eq. 'GRADSIGM') then
            lawOffset = lawOffset+4000
        end if
        if (typmod(1) .eq. '1D') then
            lawOffset = lawOffset+5000
        end if
        if (typmod(2) .eq. 'GRADVARI') then
            lawOffset = lawOffset+6000
        end if
        if (typmod(2) .eq. 'EJ_HYME' .or. typmod(2) .eq. 'ELEMJOIN') then
            lawOffset = lawOffset+7000
        end if
        if (typmod(2) .eq. 'INTERFAC') then
            lawOffset = lawOffset+9000
        end if
        if (typmod(2) .eq. 'INTSOLPI') then
            lawOffset = lawOffset+9500
        end if
        BEHInteg%behavPara%lawIndexOffset = lawOffset
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! behaviourPrepStrain
!
! Prepare input strains for the behaviour law
!    -> If defoLDC = 'MECANIQUE', prepare mechanical strain
!    -> If defoLDC = 'TOTALE' or 'OLD', keep total strain
!
! In  BEHInteg         : main object for managing the integration of behavior laws
! IO  epsm             : In : total strains at beginning of current step time
!                        Out : mechanical strains at beginning of current step time
! IO  deps             : In : increment of total strains during current step time
!                        Out : increment of mechanical strains during current step time
!
! --------------------------------------------------------------------------------------------------
    subroutine behaviourPrepStrain(BEHInteg, epsm, deps)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(Behaviour_Integ), intent(inout) :: BEHInteg
        real(kind=8), intent(inout) :: epsm(:), deps(:)
! ----- Local
        aster_logical :: lPtot, lStrainMeca
!   ------------------------------------------------------------------------------------------------
!
        if (ca_nbcvrc_ .ne. 0) then
            lStrainMeca = BEHInteg%behavPara%lStrainMeca
            lPtot = BEHInteg%allVarcStrain%list(VARC_STRAIN_PTOT)%exist
            if (lStrainMeca .or. lPtot) then
                if (LDC_PREP_DEBUG .eq. 1) then
                    WRITE (6, *) '<DEBUG>  Présence de VARC avec nouveau système ou PTOT'
                end if
! ------------- Compute non-mechanic strains for some external state variables
                call computeStrainESVA(BEHInteg%allVarcStrain, BEHInteg%behavESVA)

! ------------- Subtract to get mechanical strain epsm and deps become mechanical strains
                call computeStrainMeca(BEHInteg, epsm, deps)
            end if
        else
            if (LDC_PREP_DEBUG .eq. 1) then
                WRITE (6, *) '<DEBUG>  Absence de VARC'
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! behaviourSetParaPoin
!
! In  kpg              : index of quadrature point
! In  ksp              : index of "sub"-point (plates, pipes, beams, etc.)
! IO  BEHInteg         : main object for managing the integration of behavior laws
!
! --------------------------------------------------------------------------------------------------
    subroutine behaviourSetParaPoin(kpg, ksp, BEHInteg)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: kpg, ksp
        type(Behaviour_Integ), intent(inout) :: BEHInteg
!   ------------------------------------------------------------------------------------------------
!
        BEHInteg%materPara%schemePara%kpg = kpg
        BEHInteg%materPara%schemePara%ksp = ksp
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! behaviourInitPoint
!
! Initialisation of behaviour datastructure - Special for SIMU_POINT_MAT
!
! In  relaComp         : RELATION in COMPORTemENT
! IO  BEHInteg         : main object for managing the integration of behavior laws
!
! --------------------------------------------------------------------------------------------------
    subroutine behaviourInitPoint(relaComp, BEHInteg)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=16), intent(in) :: relaComp
        type(Behaviour_Integ), intent(inout) :: BEHInteg
!   ------------------------------------------------------------------------------------------------
!

! ----- Set "real" zero
        BEHInteg%behavESVA%behavESVAGeom%coorElga = 0.d0

! ----- Special for GRAD_VELO
        if (BEHInteg%behavESVA%tabcod(GRADVELO) .eq. 1) then
            call utmess('A', 'COMPOR2_39')
            BEHInteg%behavESVA%behavESVAGeom%gradVelo = 0.d0
        end if

! ----- Special for ELTSIZE1
        if (BEHInteg%behavESVA%tabcod(ELTSIZE1) .eq. 1) then
            if (relaComp .ne. 'BETON_DOUBLE_DP') then
                call utmess('F', 'COMPOR2_12')
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getListUserESVA
!
! Get list of external state variables from user (AFFE_VARC)
!
! In  carcri           : parameters for comportment
! Out tabcod           : list of integers to detect external state variables
!
!   ------------------------------------------------------------------------------------------------
    subroutine getListUserESVA(carcri, tabcod)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
        integer(kind=8), intent(out) :: tabcod(60)
! ----- Local
        integer(kind=8) :: jvariext1, jvariext2, variextecode(2)
!   ------------------------------------------------------------------------------------------------
!
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG>  Get list of external state variables from user'
        end if

        jvariext1 = nint(carcri(IVARIEXT1))
        jvariext2 = nint(carcri(IVARIEXT2))
        tabcod = 0
        variextecode(1) = jvariext1
        variextecode(2) = jvariext2
        call isdeco(variextecode, tabcod, 60)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! behaviourPrepESVAGeom
!
! Prepare external state variables - Geometry
!
! In  nno              : number of nodes
! In  npg              : number of Gauss points
! In  ndim             : dimension of problem (2 or 3)
! In  jv_poids         : JEVEUX adress for weight of Gauss points
! In  jv_func          : JEVEUX adress for shape functions
! In  jv_dfunc         : JEVEUX adress for derivative of shape functions
! In  geom             : initial coordinates of nodes
! IO  BEHInteg         : main object for managing the integration of behavior laws
! In  deplm            : displacements of nodes at beginning of time step
! In  ddepl            : displacements of nodes since beginning of time step
!
! --------------------------------------------------------------------------------------------------
    subroutine behaviourPrepESVAGeom(nno, npg, ndim, &
                                     jv_poids, jv_func, jv_dfunc, &
                                     geom, BEHInteg, &
                                     deplm_, ddepl_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: nno, npg, ndim
        integer(kind=8), intent(in) :: jv_poids, jv_func, jv_dfunc
        real(kind=8), intent(in) :: geom(ndim, nno)
        type(Behaviour_Integ), intent(inout) :: BEHInteg
        real(kind=8), optional, intent(in) :: deplm_(ndim, nno), ddepl_(ndim, nno)
!   ------------------------------------------------------------------------------------------------
!
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG> Preparation of external state variable for each element'
        end if

! ----- Compute element size 1
        if (BEHInteg%behavESVA%tabcod(ELTSIZE1) .eq. 1) then
            call prepEltSize1(nno, npg, ndim, &
                              jv_poids, jv_func, jv_dfunc, &
                              geom, BEHInteg%behavPara, &
                              BEHInteg%behavESVA%behavESVAGeom)
        end if

! ----- Compute gradient of velocity
        if (BEHInteg%behavESVA%tabcod(GRADVELO) .eq. 1) then
            if (.not. present(deplm_) .or. .not. present(ddepl_)) then
                call utmess('F', 'COMPOR2_26')
            end if
            call prepGradVelo(nno, npg, ndim, &
                              jv_poids, jv_func, jv_dfunc, &
                              geom, deplm_, ddepl_, &
                              BEHInteg%behavESVA%behavESVAGeom)
        end if

! ----- Coordinates of Gauss points (always)
        call behaviourCoorGauss(nno, npg, ndim, &
                                jv_func, geom, &
                                BEHInteg%behavESVA%behavESVAGeom)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! behaviourPrepESVAGeomHHO
!
! Prepare external state variables - Geometry
!
! In  nno              : number of nodes
! In  npg              : number of Gauss points
! In  ndim             : dimension of problem (2 or 3)
! In  jv_poids         : JEVEUX adress for weight of Gauss points
! In  jv_func          : JEVEUX adress for shape functions
! In  jv_dfunc         : JEVEUX adress for derivative of shape functions
! In  geom             : initial coordinates of nodes
! IO  BEHInteg         : main object for managing the integration of behavior laws
! In  deplm            : displacements of nodes at beginning of time step
! In  ddepl            : displacements of nodes since beginning of time step
!
! --------------------------------------------------------------------------------------------------
    subroutine behaviourPrepESVAGeomHHO(hhoCell, hhoQuad, BEHInteg)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Quadrature), intent(in) :: hhoQuad
        type(Behaviour_Integ), intent(inout) :: BEHInteg
!   ------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: ndim, npg
        real(kind=8) :: lc
!
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG> Preparation of external state variable for each HHO element'
        end if
!
        ndim = hhoCell%ndim
        npg = hhoQuad%nbQuadPoints

! ----- Compute element size 1
        if (BEHInteg%behavESVA%tabcod(ELTSIZE1) .eq. 1) then
            if (ndim == 3) then
                if (npg .ge. 9) then
                    lc = hhoCell%measure**0.33333333333333d0
                else
                    lc = sqrt(2.0)*hhoCell%measure**0.33333333333333d0
                end if
            else
                if (npg .ge. 5) then
                    lc = sqrt(hhoCell%measure)
                else
                    lc = sqrt(2.0*hhoCell%measure)
                end if
            end if
            ! ------ Save values
            BEHInteg%behavESVA%behavESVAGeom%lElemSize1 = ASTER_TRUE
            BEHInteg%behavESVA%behavESVAGeom%elemSize1 = lc
        end if

! ----- Compute gradient of velocity
        if (BEHInteg%behavESVA%tabcod(GRADVELO) .eq. 1) then
            ASSERT(ASTER_FALSE)
        end if

! ----- Coordinates of Gauss points (always)
        BEHInteg%behavESVA%behavESVAGeom%coorElga = 0.d0
        ASSERT(npg .le. VARC_GEOM_NBMAXI)
!
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG>  Compute coordinates of Gauss points'
        end if

        BEHInteg%behavESVA%behavESVAGeom%coorElga(1:npg, 1:ndim) = &
            hhoQuad%points(1:npg, 1:ndim)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! behaviourPrepESVAPoin
!
! Prepare external state variables at Gauss point
!
! IO  BEHInteg         : main object for managing the integration of behavior laws
!
! --------------------------------------------------------------------------------------------------
    subroutine behaviourPrepESVAPoin(BEHInteg)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(Behaviour_Integ), intent(inout) :: BEHInteg
! ----- Local
        aster_logical :: lExteSolver, lStrainMeca, lhasInelasticStrains
!   ------------------------------------------------------------------------------------------------
!
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG> Preparation of external state variable for current point'
        end if

! ----- Flags for external solvers
        lExteSolver = BEHInteg%behavPara%lExteSolver

! ----- Flag for strain decomposition
        lStrainMeca = BEHInteg%behavPara%lStrainMeca

! ----- Prepare hygrometry (not in AFFE_VARC)
        if (BEHInteg%behavESVA%tabcod(HYGR) .eq. 1) then
            call prepHygr(BEHInteg%materPara, BEHInteg%behavESVA)
        end if

! ----- Prepare external PTOT state variables (in AFFE_VARC)
        if (ca_nbcvrc_ .ne. 0) then
            if (BEHInteg%allVarcStrain%list(VARC_STRAIN_PTOT)%exist) then
                call compPtotStrainField('T', BEHInteg%materPara, &
                                         BEHInteg%allVarcStrain%list(VARC_STRAIN_PTOT))
            end if
        end if

! ----- Prepare external state variables from fields (in AFFE_VARC)
        if (ca_nbcvrc_ .ne. 0 .or. BEHInteg%behavPara%lTHM) then
            if (lExteSolver .or. lStrainMeca) then
                call compVarcStrain('T', BEHInteg%materPara, &
                                    BEHInteg%allVarcStrain)
            end if
        end if

! ----- Detect inelastic strains
        lhasInelasticStrains = ASTER_FALSE
        if (ca_nbcvrc_ .ne. 0) then
            lhasInelasticStrains = BEHInteg%allVarcStrain%hasInelasticStrains
        end if
        BEHInteg%behavESVA%lhasInelasticStrains = lhasInelasticStrains

! ----- Prepare other external state variables (For temperature: see preparation of fields)
        if (lExteSolver) then
            BEHInteg%behavESVA%behavESVAOther%time = BEHInteg%behavPara%timePrev
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getMFrontVariableName
!
! Get MFront glossary or entry name from code_aster external state variable name
!
! In  exteNameAster    : code_aster external state variable name
!
! --------------------------------------------------------------------------------------------------
    character(len=64) function getMFrontVariableName(exteNameAster)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=64), intent(in) :: exteNameAster
! ----- Local
        integer(kind=8) :: iEsva
!   ------------------------------------------------------------------------------------------------
!
        getMFrontVariableName = exteNameAster
        do iEsva = 1, ESVA_EXTE_MGIS_NBMAXI
            if (fromAsterToMFront(1, iEsva) == exteNameAster) then
                getMFrontVariableName = fromAsterToMFront(2, iEsva)
            end if
        end do
!
!   ------------------------------------------------------------------------------------------------
    end function
! --------------------------------------------------------------------------------------------------
!
! getAsterVariableName
!
! Get Aster name from MFront variable name
!
! In  exteNameMGIS     : MGIS external state variable name
!
! --------------------------------------------------------------------------------------------------
    character(len=8) function getAsterVariableName(exteNameMGIS)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=64), intent(in) :: exteNameMGIS
! ----- Local
        character(len=64) :: exteNameAster
        integer(kind=8) :: iEsva
!   ------------------------------------------------------------------------------------------------
!
        exteNameAster = exteNameMGIS
        do iEsva = 1, ESVA_EXTE_MGIS_NBMAXI
            if (fromAsterToMFront(2, iEsva) == exteNameMGIS) then
                exteNameAster = fromAsterToMFront(1, iEsva)
            end if
        end do
        getAsterVariableName = exteNameAster(1:8)
!
!   ------------------------------------------------------------------------------------------------
    end function
! --------------------------------------------------------------------------------------------------
!
! logUndefinedVariable
!
! Log the lack of definition of external state variable
!
! In  exteNameAster    : variable name
!
! --------------------------------------------------------------------------------------------------
    subroutine logUndefinedVariable(exteNameAster)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=8), intent(in) :: exteNameAster
! ----- Local
        character(len=64) :: exteNameAster64, MGISName
!   ------------------------------------------------------------------------------------------------
!
        exteNameAster64 = exteNameAster
        MGISName = getMFrontVariableName(exteNameAster64)
        if (MGISName .eq. exteNameAster64) then
            call utmess('F', 'COMPOR4_23', sk=exteNameAster)
        else
            call utmess('F', 'COMPOR4_75', nk=2, valk=[exteNameAster, MGISName])
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! behaviourPrepESVAExte
!
! Prepare external state variables for external solvers (UMAT/MFRONT)
!
! IO  BEHInteg         : main object for managing the integration of behavior laws
!
! --------------------------------------------------------------------------------------------------
    subroutine behaviourPrepESVAExte(BEHInteg)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(Behaviour_Integ), intent(inout) :: BEHInteg
! ----- Local
        character(len=8), parameter :: nameUMAT(VARC_EXTE_NBMAXI) = (/ &
                                       'SECH    ', 'HYDR    ', 'IRRA    ', &
                                       'NEUT1   ', 'NEUT2   ', 'CORR    ', &
                                       'ALPHPUR ', 'ALPHBETA'/)
        character(len=64) :: exteNameESVA(VARC_EXTE_NBMAXI)
        character(len=64) :: exteNameMGIS
        character(len=8) :: exteNameAster
        real(kind=8) :: valePrev, valeCurr
        integer(kind=8) :: iESVA, nbESVA, indxVarcStrain, iret
        aster_logical :: lMGIS, lUMAT, exist
        character(len=16) :: mgisAddr
!   ------------------------------------------------------------------------------------------------
!
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG> Preparation of external state variable for external solvers'
        end if

! ----- Get parameters
        lMGIS = BEHInteg%behavPara%lMGIS
        lUMAT = BEHInteg%behavPara%lUMAT

! ----- Get list of external state variables in external solvers
        nbESVA = 0
        if (lMGIS) then
            mgisAddr = BEHInteg%behavESVA%behavESVAExte%mgisAddr
            call mgis_get_number_of_esvs(mgisAddr, nbESVA)
            ASSERT(nbESVA .le. VARC_EXTE_NBMAXI)
            call mgis_get_esvs(mgisAddr, exteNameESVA)
        elseif (lUMAT) then
            exteNameESVA = nameUMAT
            nbESVA = VARC_EXTE_NBMAXI
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- Default: all external states variables are scalar (not strain )
        BEHInteg%behavESVA%behavESVAExte%nbESVAScal = nbESVA

        if (LDC_PREP_DEBUG .eq. 1) then
            if (nbESVA .eq. 0) then
                WRITE (6, *) '<DEBUG> No external state variables defined in MFront/UMAT'
            else
                WRITE (6, *) '<DEBUG> Number of external state variables defined in MFront/UMAT:', &
                    nbESVA
            end if
        end if

! ----- Set values of ExternalStateVariables
        BEHInteg%behavESVA%behavESVAExte%scalESVAPrev = 0.d0
        BEHInteg%behavESVA%behavESVAExte%scalESVAIncr = 0.d0

        do iESVA = 1, nbESVA
! --------- Translate names of external state variable
            exteNameAster = getAsterVariableName(exteNameESVA(iESVA))
            exteNameMGIS = getMFrontVariableName(exteNameESVA(iESVA))

            if (LDC_PREP_DEBUG .eq. 1) then
                WRITE (6, *) '<DEBUG> External state variable:', iESVA
                WRITE (6, *) '<DEBUG>  Aster name: ', exteNameAster
                if (lMGIS) then
                    WRITE (6, *) '<DEBUG>  MGIS name: ', exteNameMGIS
                end if
                if (lUMAT) then
                    WRITE (6, *) '<DEBUG>  UMAT name: ', exteNameESVA(iESVA)
                end if
            end if

            select case (exteNameMGIS)
            case ('GRADVELO')
                call utmess('F', 'COMPOR4_25', sk=exteNameAster)

            case ('ConcreteDrying')
                indxVarcStrain = VARC_STRAIN_SECH
                exist = BEHInteg%allVarcStrain%list(indxVarcStrain)%exist
                if (exist) then
                    BEHInteg%behavESVA%behavESVAExte%scalESVAPrev(iESVA) = &
                        BEHInteg%allVarcStrain%list(indxVarcStrain)%varcPrev(1)
                    BEHInteg%behavESVA%behavESVAExte%scalESVAIncr(iESVA) = &
                        BEHInteg%allVarcStrain%list(indxVarcStrain)%varcIncr(1)
                else
                    if (.not. lUMAT) then
                        call logUndefinedVariable(exteNameAster)
                    end if
                end if

            case ('ConcreteHydration')
                indxVarcStrain = VARC_STRAIN_HYDR
                exist = BEHInteg%allVarcStrain%list(indxVarcStrain)%exist
                if (exist) then
                    BEHInteg%behavESVA%behavESVAExte%scalESVAPrev(iESVA) = &
                        BEHInteg%allVarcStrain%list(indxVarcStrain)%varcPrev(1)
                    BEHInteg%behavESVA%behavESVAExte%scalESVAIncr(iESVA) = &
                        BEHInteg%allVarcStrain%list(indxVarcStrain)%varcIncr(1)
                else
                    if (BEHInteg%behavESVA%behavESVAOther%lHygr) then
                        BEHInteg%behavESVA%behavESVAExte%scalESVAPrev(iESVA) = 0.d0
                        BEHInteg%behavESVA%behavESVAExte%scalESVAIncr(iESVA) = 0.d0
                    else
                        if (.not. lUMAT) then
                            call logUndefinedVariable(exteNameAster)
                        end if
                    end if
                end if

            case ('Hygrometry')
                ASSERT(lMGIS)
                exist = BEHInteg%behavESVA%behavESVAOther%lHygr
                if (exist) then
                    BEHInteg%behavESVA%behavESVAExte%scalESVAPrev(iESVA) = &
                        BEHInteg%behavESVA%behavESVAOther%hygrPrev
                    BEHInteg%behavESVA%behavESVAExte%scalESVAIncr(iESVA) = &
                        BEHInteg%behavESVA%behavESVAOther%hygrIncr
                else
                    call utmess('F', 'COMPOR4_26', sk=exteNameAster)
                end if

            CASE ('Temperature')
                ASSERT(lMGIS)
                indxVarcStrain = VARC_STRAIN_TEMP
                exist = BEHInteg%allVarcStrain%list(indxVarcStrain)%exist
                if (exist) then
                    BEHInteg%behavESVA%behavESVAExte%scalESVAPrev(iESVA) = &
                        BEHInteg%allVarcStrain%list(indxVarcStrain)%varcPrev(1)
                    BEHInteg%behavESVA%behavESVAExte%scalESVAIncr(iESVA) = &
                        BEHInteg%allVarcStrain%list(indxVarcStrain)%varcIncr(1)
                else
                    call utmess('A', 'COMPOR4_26', sk=exteNameAster)
                end if

            case ('ElementSize')
                ASSERT(lMGIS)
                ASSERT(BEHInteg%behavESVA%behavESVAGeom%lElemSize1)
                BEHInteg%behavESVA%behavESVAExte%scalESVAPrev(iESVA) = &
                    BEHInteg%behavESVA%behavESVAGeom%elemSize1

            CASE ('ReferenceTemperature')
                ASSERT(lMGIS)
                indxVarcStrain = VARC_STRAIN_TEMP
                exist = BEHInteg%allVarcStrain%list(indxVarcStrain)%exist
                if (exist) then
                    BEHInteg%behavESVA%behavESVAExte%scalESVAPrev(iESVA) = &
                        BEHInteg%allVarcStrain%list(indxVarcStrain)%varcRefe
                else
                    call utmess('F', 'COMPOR4_26', sk=exteNameAster)
                end if

            CASE ('Time')
                ASSERT(lMGIS)
                BEHInteg%behavESVA%behavESVAExte%scalESVAPrev(iESVA) = &
                    BEHInteg%behavESVA%behavESVAOther%time

            CASE ('TEMP')
                ASSERT(lUMAT)
                indxVarcStrain = VARC_STRAIN_TEMP
                exist = BEHInteg%allVarcStrain%list(indxVarcStrain)%exist
                if (exist) then
                    BEHInteg%behavESVA%behavESVAExte%scalESVAPrev(iESVA) = &
                        BEHInteg%allVarcStrain%list(indxVarcStrain)%varcPrev(1)
                    BEHInteg%behavESVA%behavESVAExte%scalESVAIncr(iESVA) = &
                        BEHInteg%allVarcStrain%list(indxVarcStrain)%varcIncr(1)
                else
                    call utmess('F', "COMPOR4_76")
                end if

            case default
                call rcvarc(' ', exteNameAster, '-', &
                            BEHInteg%materPara%schemePara%fami, &
                            BEHInteg%materPara%schemePara%kpg, &
                            BEHInteg%materPara%schemePara%ksp, &
                            valePrev, iret)

                call rcvarc(' ', exteNameAster, '-', &
                            BEHInteg%materPara%schemePara%fami, &
                            BEHInteg%materPara%schemePara%kpg, &
                            BEHInteg%materPara%schemePara%ksp, &
                            valePrev, iret)
                if (iret .eq. 0) then
                    iret = 0
                    call rcvarc('F', exteNameAster, '+', &
                                BEHInteg%materPara%schemePara%fami, &
                                BEHInteg%materPara%schemePara%kpg, &
                                BEHInteg%materPara%schemePara%ksp, &
                                valeCurr, iret)
                    BEHInteg%behavESVA%behavESVAExte%scalESVAPrev(iESVA) = valePrev
                    BEHInteg%behavESVA%behavESVAExte%scalESVAIncr(iESVA) = valeCurr-valePrev
                else
                    if (.not. lUMAT) then
                        call logUndefinedVariable(exteNameAster)
                    end if
                end if
            end select
            if (LDC_PREP_DEBUG .eq. 1) then
                WRITE (6, *) '<DEBUG>  Valeur précédente   : ', &
                    BEHInteg%behavESVA%behavESVAExte%scalESVAPrev(iESVA)
                WRITE (6, *) '<DEBUG>  Valeur incrémentale : ', &
                    BEHInteg%behavESVA%behavESVAExte%scalESVAIncr(iESVA)
            end if
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! prepHygr
!
! Prepare external state variable HYGR (specific)
!
! In  materPara        : parameters of material
! IO  behavESVA        : parameters for External State Variables
!
! --------------------------------------------------------------------------------------------------
    subroutine prepHygr(materPara, behavESVA)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(Material_Para), intent(in) :: materPara
        type(BehaviourESVA), intent(inout) :: behavESVA
! ----- Local
        integer(kind=8), parameter :: nbPara = 1
        integer(kind=8)  :: codret(nbPara)
        real(kind=8) :: paraVale(nbPara)
        character(len=16), parameter :: paraName(nbPara) = (/'FONC_DESORP'/)
        character(len=8) :: fami
        integer(kind=8) :: kpg, ksp, jvMaterCode
        character(len=16) :: phenom
        real(kind=8) :: funcDesorpPrev, funcDesorpCurr
        aster_logical :: exist
!   ------------------------------------------------------------------------------------------------
!
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG>  Prepare HYGR'
        end if

        fami = materPara%schemePara%fami
        kpg = materPara%schemePara%kpg
        ksp = materPara%schemePara%ksp
        jvMaterCode = materPara%jvMaterCode

! ----- Get parameters
        funcDesorpPrev = 0.d0
        funcDesorpCurr = 0.d0
        exist = ASTER_FALSE
        call rccoma(jvMaterCode, 'BETON_DESORP', 0, phenom, codret(1))
        if (codret(1) .ne. 0) then
            call utmess('F', 'COMPOR2_94')
        end if
        exist = ASTER_TRUE
        call rcvalb(fami, kpg, ksp, '-', jvMaterCode, &
                    ' ', 'BETON_DESORP', 0, ' ', [0.d0], &
                    nbPara, paraName, paraVale, codret, 0)
        if (codret(1) .eq. 0) then
            funcDesorpPrev = paraVale(1)
            call rcvalb(fami, kpg, ksp, '+', jvMaterCode, &
                        ' ', 'BETON_DESORP', 0, ' ', [0.d0], &
                        nbPara, paraName, paraVale, codret, 0)
            ASSERT(codret(1) .eq. 0)
            funcDesorpCurr = paraVale(1)
        else
            !   leverett isotherm
            call leverettIsotMeca(fami, kpg, ksp, jvMaterCode, funcDesorpPrev, &
                                  funcDesorpCurr)
        end if

! ----- Save values
        if (exist) then
            behavESVA%behavESVAOther%lHygr = exist
            behavESVA%behavESVAOther%hygrPrev = funcDesorpPrev
            behavESVA%behavESVAOther%hygrIncr = funcDesorpCurr-funcDesorpPrev
        end if

! ----- Debug
        if (LDC_PREP_DEBUG .eq. 1) then
            if (behavESVA%behavESVAOther%lHygr) then
                WRITE (6, *) '<DEBUG>  Values of HYGR: ', &
                    behavESVA%behavESVAOther%hygrPrev, &
                    behavESVA%behavESVAOther%hygrIncr
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! behaviourPredictionStress
!
! Compute the contribution of the "thermal stress increment" in prediction
!
! in  behavESVA        : parameters for external state variables
! in  dsidep           : tangent operator in prediction
! IO  sig              : stress (Taylor term of order 0)
!
! --------------------------------------------------------------------------------------------------
    subroutine behaviourPredictionStress(behavESVA, dsidep, sig)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(BehaviourESVA), intent(in) :: behavESVA
        real(kind=8), intent(in) :: dsidep(:, :)
        real(kind=8), intent(inout) :: sig(:)
! ----- Local
        integer(kind=8) :: ndimsi
!   ------------------------------------------------------------------------------------------------
!
        if (ca_nbcvrc_ .ne. 0) then
            ndimsi = size(sig)
            sig = sig-matmul(dsidep, behavESVA%depsi_varc(1:ndimsi))
            if (LDC_PREP_DEBUG .eq. 1) then
                WRITE (6, *) '<DEBUG>  Prediction stress: ', sig
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! prepEltSize1
!
! Compute external state variables - Geometry: size of element (ELTSIZE1)
!
! In  nno              : number of nodes
! In  npg              : number of Gauss points
! In  ndim             : dimension of problem (2 or 3)
! In  jv_poids         : JEVEUX adress for weight of Gauss points
! In  jv_func          : JEVEUX adress for shape functions
! In  jv_dfunc         : JEVEUX adress for derivative of shape functions
! In  geom             : initial coordinates of nodes
! In  behavPara        : parameters for integration of behaviour
! IO  behavESVAGeom    : geometric properties
!
! --------------------------------------------------------------------------------------------------
    subroutine prepEltSize1(nno, npg, ndim, &
                            jv_poids, jv_func, jv_dfunc, &
                            geom, behavPara, &
                            behavESVAGeom)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: nno, npg, ndim
        integer(kind=8), intent(in) :: jv_poids, jv_func, jv_dfunc
        real(kind=8), intent(in) :: geom(ndim, nno)
        type(Behaviour_Para), intent(in) :: behavPara
        type(BehaviourESVA_Geom), intent(inout) :: behavESVAGeom
! ----- Local
        integer(kind=8) :: kpg, iComp, i
        real(kind=8) :: lc, dfdx(MT_NNOMAX3D), dfdy(MT_NNOMAX3D), dfdz(MT_NNOMAX3D), poids, r
        real(kind=8) :: volume, surfac
        real(kind=8), parameter :: rac2 = sqrt(2.d0)
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(nno .le. MT_NNOMAX3D)
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG>  Compute ELTSIZE1'
        end if
!
        if (behavPara%lThreeDim) then
            volume = 0.d0
            do kpg = 1, npg
                call dfdm3d(nno, kpg, jv_poids, jv_dfunc, geom, &
                            poids, dfdx, dfdy, dfdz)
                volume = volume+poids
            end do
            if (npg .ge. 9) then
                lc = volume**0.33333333333333d0
            else
                lc = rac2*volume**0.33333333333333d0
            end if
        elseif (behavPara%lAxis .or. behavPara%lPlaneStrain) then
            surfac = 0.d0
            do kpg = 1, npg
                iComp = (kpg-1)*nno
                call dfdm2d(nno, kpg, jv_poids, jv_dfunc, geom, &
                            poids, dfdx, dfdy)
                if (behavPara%lAxis) then
                    r = 0.d0
                    do i = 1, nno
                        r = r+geom(1, i)*zr(jv_func+i+iComp-1)
                    end do
                    poids = poids*r
                end if
                surfac = surfac+poids
            end do
            if (npg .ge. 5) then
                lc = surfac**0.5d0
            else
                lc = rac2*surfac**0.5d0
            end if
        elseif (behavPara%lPlaneStress) then
            lc = r8nnem()
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- Save values
        behavESVAGeom%lElemSize1 = ASTER_TRUE
        behavESVAGeom%elemSize1 = lc
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! prepGradVelo
!
! Compute external state variables - Geometry: gradient of velocity (GRADVELO)
!
! In  nno              : number of nodes
! In  npg              : number of Gauss points
! In  ndim             : dimension of problem (2 or 3)
! In  jv_poids         : JEVEUX adress for weight of Gauss points
! In  jv_func          : JEVEUX adress for shape functions
! In  jv_dfunc         : JEVEUX adress for derivative of shape functions
! In  geom             : initial coordinates of nodes
! In  deplm            : displacements of nodes at beginning of time step
! In  ddepl            : displacements of nodes since beginning of time step
! IO  behavESVAGeom    : geometric properties
!
! --------------------------------------------------------------------------------------------------
    subroutine prepGradVelo(nno, npg, ndim, &
                            jv_poids, jv_func, jv_dfunc, &
                            geom, deplm, ddepl, &
                            behavESVAGeom)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: nno, npg, ndim
        integer(kind=8), intent(in) :: jv_poids, jv_func, jv_dfunc
        real(kind=8), intent(in) :: geom(ndim, nno)
        real(kind=8), intent(in) :: deplm(ndim, nno), ddepl(ndim, nno)
        type(BehaviourESVA_Geom), intent(inout) :: behavESVAGeom
! ----- Local
        integer(kind=8) :: nddl, kpg, i, j
        real(kind=8) :: l(3, 3), fmm(3, 3), df(3, 3), f(3, 3), r8bid, r
        real(kind=8) :: deplp(3, MT_NNOMAX3D), geomm(3, MT_NNOMAX3D), epsbid(6)
        real(kind=8) :: dfdi(nno, 3)
        real(kind=8), parameter :: id(9) = (/1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0/)
        blas_int :: b_incx, b_incy, b_n
!   ------------------------------------------------------------------------------------------------
!
        nddl = ndim*nno
        behavESVAGeom%gradVelo = 0.d0
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG>  Compute GRADVELO'
        end if
!
        b_n = to_blas_int(nddl)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, geom, b_incx, geomm, b_incy)
        b_n = to_blas_int(nddl)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0, deplm, b_incx, geomm, &
                   b_incy)
        b_n = to_blas_int(nddl)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call dcopy(b_n, deplm, b_incx, deplp, b_incy)
        b_n = to_blas_int(nddl)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        call daxpy(b_n, 1.d0, ddepl, b_incx, deplp, &
                   b_incy)
        do kpg = 1, npg
            call nmgeom(ndim, nno, .false._1, .true._1, geom, &
                        kpg, jv_poids, jv_func, jv_dfunc, deplp, &
                        .true._1, r8bid, dfdi, f, epsbid, &
                        r)
            call nmgeom(ndim, nno, .false._1, .true._1, geomm, &
                        kpg, jv_poids, jv_func, jv_dfunc, ddepl, &
                        .true._1, r8bid, dfdi, df, epsbid, &
                        r)
            b_n = to_blas_int(9)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call daxpy(b_n, -1.d0, id, b_incx, df, &
                       b_incy)
            call matinv('S', 3, f, fmm, r8bid)
            l = matmul(df, fmm)
            do i = 1, 3
                do j = 1, 3
                    behavESVAGeom%gradVelo(3*(i-1)+j) = l(i, j)
                end do
            end do
        end do
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! behaviourCoorGauss
!
! Compute external state variables - Geometry: coordinates of Gauss points
!
! In  nno              : number of nodes
! In  npg              : number of Gauss points
! In  ndim             : dimension of problem (2 or 3)
! In  jv_func          : JEVEUX adress for shape functions
! In  geom             : initial coordinates of nodes
! IO  behavESVAGeom    : geometric properties
!
! --------------------------------------------------------------------------------------------------
    subroutine behaviourCoorGauss(nno, npg, ndim, &
                                  jv_func, geom, &
                                  behavESVAGeom)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: nno, npg, ndim
        integer(kind=8), intent(in) :: jv_func
        real(kind=8), intent(in) :: geom(ndim, nno)
        type(BehaviourESVA_Geom), intent(inout) :: behavESVAGeom
! ----- Local
        integer(kind=8) :: i, iComp, kpg
!   ------------------------------------------------------------------------------------------------
!
        behavESVAGeom%coorElga = 0.d0
        ASSERT(npg .le. VARC_GEOM_NBMAXI)
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG>  Compute coordinates of Gauss points'
        end if
!
        do kpg = 1, npg
            do i = 1, ndim
                do iComp = 1, nno
                    behavESVAGeom%coorElga(kpg, i) = &
                        behavESVAGeom%coorElga(kpg, i)+ &
                        geom(i, iComp)*zr(jv_func-1+nno*(kpg-1)+iComp)
                end do
            end do
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! computeStrainESVA
!
! Pre-compute non-mechanical strains from external state variables
!
! In  allVarcStrain    : all external state variables for anelastic strains
! IO  behavESVA        : parameters for External State Variables
! In  neps             : number of components of strains
!
! --------------------------------------------------------------------------------------------------
    subroutine computeStrainESVA(allVarcStrain, behavESVA)
!   -----------------------------------------------------------------------------------------------
! ----- Parameters
        type(All_Varc_Strain), intent(in) :: allVarcStrain
        type(BehaviourESVA), intent(inout) :: behavESVA
! ----- Local
        aster_logical :: hasAnelastiStrains
!   ------------------------------------------------------------------------------------------------
!
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG>  Pre-compute non-mechanical strains from external state variables'
        end if

        behavESVA%depsi_varc = 0.d0
        behavESVA%epsi_varc = 0.d0
        hasAnelastiStrains = ASTER_FALSE

        call getVarcStrain('-', VARC_STRAIN_ALL, allVarcStrain, behavESVA%epsi_varc)
        call getVarcStrain('T', VARC_STRAIN_ALL, allVarcStrain, behavESVA%depsi_varc)

! ----- DEBUG
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG>  Strains from external state variables - Prev: ', &
                behavESVA%epsi_varc
            WRITE (6, *) '<DEBUG>  Strains from external state variables - Incr: ', &
                behavESVA%depsi_varc
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! computeStrainMeca
!
! Prepare strains (substracting "thermic" strains to total strains to get mechanical part)
!
! In  BEHInteg         : main object for managing the integration of behavior laws
! IO  epsm             : In : total strains at beginning of current step time
!                        Out : mechanical strains at beginning of current step time
! IO  deps             : In : increment of total strains during current step time
!                        Out : increment of mechanical strains during current step time
!
! --------------------------------------------------------------------------------------------------
    subroutine computeStrainMeca(BEHInteg, epsm, deps)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(Behaviour_Integ), intent(in) :: BEHInteg
        real(kind=8), intent(inout) :: epsm(:), deps(:)
! ----- Local
        integer(kind=8):: ndimsi
!   ------------------------------------------------------------------------------------------------
!
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG> Calcul des déformations mécaniques'
        end if
        ASSERT(BEHInteg%behavPara%lEpsComEnabled)

        ndimsi = BEHInteg%behavPara%ndimsi
        epsm(1:ndimsi) = epsm(1:ndimsi)-BEHInteg%behavESVA%epsi_varc(1:ndimsi)
        deps(1:ndimsi) = deps(1:ndimsi)-BEHInteg%behavESVA%depsi_varc(1:ndimsi)

        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG>  Prepare strains for integration: ', epsm, deps
        end if
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! varcIsGEOM
!
! Detect 'GEOM' in external state variables
!
! Out lGeomInESVA      : flag for GEOM in external state variables (AFFE_VARC)
!
! --------------------------------------------------------------------------------------------------
    subroutine varcIsGEOM(lGeomInESVA)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        aster_logical, intent(inout) :: lGeomInESVA
! ----- Local
        character(len=8), parameter :: varc_geom = 'X'
        integer(kind=8) :: varc_indx
!   ------------------------------------------------------------------------------------------------
!
        lGeomInESVA = ASTER_FALSE

! ----- Detect 'GEOM' external state variables
        if (ca_nbcvrc_ .eq. 0) then
            lGeomInESVA = ASTER_FALSE
        else
            varc_indx = indik8(zk8(ca_jvcnom_), varc_geom, 1, ca_nbcvrc_)
            lGeomInESVA = varc_indx .ne. 0
        end if
!
        if (LDC_PREP_DEBUG .eq. 1) then
            WRITE (6, *) '<DEBUG> Detect GEOM in external state variables: ', lGeomInESVA
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! isSolverIsExte
!
! Detect if external solver (MFront, UMAT) is used
!
! In  carcri           : parameters for comportment
! Out lMGIS            : logical for MFront/MGIS behaviour
! Out lUMAT            : logical for UMAT behaviour
!
! --------------------------------------------------------------------------------------------------
    subroutine isSolverIsExte(carcri, lMGIS, lUMAT)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
        aster_logical, intent(out) :: lMGIS, lUMAT
!   ------------------------------------------------------------------------------------------------
        lMGIS = ASTER_FALSE
        lUMAT = ASTER_FALSE
        if (nint(carcri(EXTE_TYPE)) .eq. 1) then
!       MFront official
            lMGIS = ASTER_TRUE
        elseif (nint(carcri(EXTE_TYPE)) .eq. 2) then
!       MFront proto
            lMGIS = ASTER_TRUE
        elseif (nint(carcri(EXTE_TYPE)) .eq. 4) then
            lUMAT = ASTER_TRUE
        end if
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! behaviourOption
!
! Select objects to construct from option name
!
! In  option           : name of option to compute
! In  compor           : name of comportment definition (field)
! Out lMatr            : flag when tangent matrix
! Out lVect            : flag when internal forces vector
! Out lVari            : flag when internal state variables
! Out lSigm            : flag when stress and return code error
! Out codret           : return code when integrate behaviour
!
! --------------------------------------------------------------------------------------------------
    subroutine behaviourOption(option, compor, &
                               lMatr, lVect, &
                               lVari, lSigm, &
                               codret_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=16), intent(in) :: option
        character(len=16), intent(in) :: compor(*)
        aster_logical, intent(out) :: lMatr, lVect, lVari, lSigm
        integer(kind=8), optional, intent(out) :: codret_
! ----- Local
!   ------------------------------------------------------------------------------------------------
        aster_logical :: lPred
        integer(kind=8) :: copred, jv_copred, codret
!   ------------------------------------------------------------------------------------------------
        lVari = L_VARI(option)
        lSigm = L_SIGM(option)
        lVect = L_VECT(option)
        lMatr = L_MATR(option)
        lPred = L_PRED(option)
        if (lPred) then
            copred = 0
            if (compor(DEFO_LDC) .eq. 'MECANIQUE') then
                copred = 1
            end if
            call jevech('PCOPRED', 'E', jv_copred)
            zi(jv_copred) = copred
        end if
        codret = 0
        if (present(codret_)) then
            codret_ = codret
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! behaviourGetParameters
!
! Get parameters from behaviour (linear and non-linear cases)
!
! Out relaComp         : RELATION in behaviour
! Out defoComp         : DEFORMATION in behaviour
!
! --------------------------------------------------------------------------------------------------
    subroutine behaviourGetParameters(relaName, defoComp)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=16), intent(out) :: relaName, defoComp
! ----- Local
        integer(kind=8) :: iret, jtab(7)
!   ------------------------------------------------------------------------------------------------
!
        relaName = " "
        defoComp = " "
        call tecach('NNO', 'PCOMPOR', 'L', iret, nval=7, itab=jtab)
        relaName = 'ELAS'
        defoComp = 'PETIT'
        if (iret .eq. 0) then
            relaName = zk16(jtab(1)-1+RELA_NAME)
            defoComp = zk16(jtab(1)-1+DEFO)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! chckBounds
!
! Detect VERI_BORNE
!
! IO  BEHInteg         : main object for managing the integration of behavior laws
!
! --------------------------------------------------------------------------------------------------
    subroutine chckBounds(BEHInteg)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(Behaviour_Integ), intent(inout) :: BEHInteg
! ----- Localfami, jvMaterCode,
        integer(kind=8), parameter :: kpg = 1, ksp = 1
        character(len=8), parameter :: materPoin = ' '
        integer(kind=8), parameter :: nbProp = 4
        character(len=16), parameter :: propName(nbProp) = (/'EPSI_MAXI', 'VEPS_MAXI', &
                                                             'TEMP_MINI', 'TEMP_MAXI'/)
        real(kind=8) :: propVale(nbProp)
        integer(kind=8) :: propCode(nbProp)
!   ------------------------------------------------------------------------------------------------
!
        call rcvalb(BEHInteg%materPara%schemePara%fami, &
                    kpg, ksp, '+', &
                    BEHInteg%materPara%jvMaterCode, &
                    materPoin, 'VERI_BORNE', 0, ' ', [0.d0], &
                    nbProp, propName, propVale, propCode, 0)
        BEHInteg%behavPara%lChckBounds = propCode(1) .eq. 0 .or. &
                                         propCode(2) .eq. 0 .or. &
                                         propCode(3) .eq. 0 .or. &
                                         propCode(4) .eq. 0
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module Behaviour_module
