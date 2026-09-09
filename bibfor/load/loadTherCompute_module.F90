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
! Module for the management of computation of loads (thermic)
!
! ==================================================================================================
!
module loadTherCompute_module
! ==================================================================================================
    use loadTherCompute_type
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: hasFuncLoad, isTherLoadExist
    public :: prepGeneralFields
    public :: compLoadVect, compLoadEvolVect
    public :: compLoadResi, compLoadEvolResi
    public :: compLoadMatr, compLoadEvolMatr
    private :: compLoadVectType, compLoadResiType, compLoadMatrType, getTherNeumField
    private :: getNeumLoadType, getLigrelToUse, getFieldFromEvol
    private :: getRHSOption, getResiOption, getMatrOption
    private :: prepVectFields, prepResiFields, prepMatrFields
! ==================================================================================================
    private
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/corich.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/exixfe.h"
#include "asterfort/gcnco2.h"
#include "asterfort/gettco.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/megeom.h"
#include "asterfort/multResuElem.h"
#include "asterfort/reajre.h"
#include "asterfort/rsinch.h"
#include "asterfort/utmess.h"
#include "asterfort/xajcin.h"
#include "LoadTypes_type.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! hasFuncLoad
!
! Detect non-constant loads
!
! In  listLoad          : list of loads
! Out coefCste          : flag for non-constant loads
!
! --------------------------------------------------------------------------------------------------
    subroutine hasFuncLoad(listLoadZ, coefCste)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: listLoadZ
        aster_logical, intent(out) :: coefCste
! ----- Local
        integer(kind=8), pointer :: listLoadInfo(:) => null()
        character(len=24) :: loadInfoJv
        integer(kind=8) :: nbLoad, iLoad
!   ------------------------------------------------------------------------------------------------
!
        coefCste = ASTER_TRUE
        loadInfoJv = listLoadZ(1:19)//'.INFC'
        call jeveuo(loadInfoJv, "L", vi=listLoadInfo)
        nbLoad = listLoadInfo(1)
        do iLoad = 1, nbLoad
            if (listLoadInfo(nbLoad+iLoad+1) .eq. 3) then
                coefCste = ASTER_FALSE
            end if
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! prepGeneralFields
!
! Prepare input fields for computation of thermal loads
!
! In  model             : name of model
! In  mateco            : mane of coded material
! In  varcCurr          : command variable for current time
! In  tempPrev          : previous temperature
! In  tempIter          : temperature field at current Newton iteration
! Out nbFieldInGene     : number of input fields (generic for all loads)
! Out lpain             : list of input parameters
! Out lchin             : list of input fields
!
! --------------------------------------------------------------------------------------------------
    subroutine prepGeneralFields(modelZ, matecoZ, &
                                 varcCurrZ, tempPrevZ, tempIterZ, &
                                 nbFieldInGene, lpain, lchin)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: modelZ, matecoZ
        character(len=*), intent(in) :: varcCurrZ, tempPrevZ, tempIterZ
        integer(kind=8), intent(out) :: nbFieldInGene
        character(len=*), intent(out) :: lpain(LOAD_NEUT_NBMAXIN), lchin(LOAD_NEUT_NBMAXIN)
! ----- Local
        character(len=24) :: chgeom
!   ------------------------------------------------------------------------------------------------
!
        nbFieldInGene = 0
        lpain = " "
        lchin = " "

! ----- Prepare field for geometry
        call megeom(modelZ, chgeom)
        nbFieldInGene = 0

! ----- Standard fields
        lpain(1) = 'PGEOMER'
        lchin(1) = chgeom
        lpain(2) = 'PTEMPER'
        lchin(2) = tempPrevZ
        lpain(3) = 'PTEMPEI'
        lchin(3) = tempIterZ
        lpain(4) = 'PMATERC'
        lchin(4) = matecoZ
        lpain(5) = 'PVARCPR'
        lchin(5) = varcCurrZ
        nbFieldInGene = 5

        ASSERT(nbFieldInGene .le. LOAD_NEUT_NBMAXIN)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compLoadVect
!
! Computation of thermal loads - Vector
!
! In  typeTher          : type of thermics
!                         'MOVE' for moving sources
!                         'STAT' if not
! In  modelZ            ; model
! In  timeMap           : time (<CARTE>)
! In  timeMove          : modified time (<CARTE>) for THER_NON_LINE_MO
! In  iLoad             : index of current load
! In  loadNume          : identification of load type
! In  loadPreObject     : base JEVEUX name for object
! In  loadLigrel        : ligrel for load
! In  nbFieldInGene     : number of input fields (generic)
! IO  lpain             : list of input parameters
! IO  lchin             : list of input fields
! In  jvBase            : JEVEUX base to create vector
! IO  resuElem          : name of elementary results
! In  vectElem          : name of elementary vectors
!
! --------------------------------------------------------------------------------------------------
    subroutine compLoadVect(typeTher, &
                            modelZ, timeMapZ, timeMoveZ, &
                            iLoad, loadNume, &
                            loadPreObjectZ, loadLigrelZ, &
                            nbFieldInGene, lpain, lchin, &
                            jvBase, resuElemZ, vectElemZ)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=4), intent(in) :: typeTher
        integer(kind=8), intent(in) :: iLoad, loadNume
        character(len=*), intent(in) :: modelZ, timeMapZ, timeMoveZ
        character(len=*), intent(in) :: loadPreObjectZ, loadLigrelZ
        integer(kind=8), intent(in) :: nbFieldInGene
        character(len=*), intent(inout) :: lpain(LOAD_NEUT_NBMAXIN), lchin(LOAD_NEUT_NBMAXIN)
        character(len=1), intent(in) :: jvBase
        character(len=*), intent(inout) :: resuElemZ
        character(len=*), intent(in) :: vectElemZ
! ----- Local
        integer(kind=8), parameter :: nbInputField = 0
        character(len=24), parameter :: inputLoadField(2) = &
                                        (/'                        ', &
                                          '                        '/)
        integer(kind=8) :: indxNeutType
!   ------------------------------------------------------------------------------------------------
!
        do indxNeutType = 1, LOAD_NEUT_NBTYPE
            call compLoadVectType(typeTher, &
                                  modelZ, timeMapZ, timeMoveZ, &
                                  indxNeutType, iLoad, loadNume, &
                                  loadPreObjectZ, loadLigrelZ, &
                                  nbInputField, inputLoadField, &
                                  nbFieldInGene, lpain, lchin, &
                                  jvBase, resuElemZ, vectElemZ)
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compLoadVectType
!
! Computation of specific thermal load type - Vector
!
! In  indxNeutType      : index of the type
! In  typeTher          : type of thermics
!                         'MOVE' for moving sources
!                         'STAT' if not
! In  modelZ            ; model
! In  timeMap           : time (<CARTE>)
! In  timeMove          : modified time (<CARTE>) for THER_NON_LINE_MO
! In  indxNeutType      : index of the type
! In  iLoad             : index of current load
! In  loadNume          : identification of load type
! In  loadPreObject     : base JEVEUX name for object
! In  loadLigrel        : ligrel for load
! In  nbInputField      : number of input fields given by load (for EVOL_CHAR)
! In  inputLoadFieldZ   : name of input fields given by load (for EVOL_CHAR)
! In  nbFieldInGene     : number of input fields (generic)
! IO  lpain             : list of input parameters
! IO  lchin             : list of input fields
! In  jvBase            : JEVEUX base to create vectors
! IO  resuElem          : name of elementary results
! In  vectElem          : name of elementary vectors
!
! --------------------------------------------------------------------------------------------------
    subroutine compLoadVectType(typeTher, &
                                modelZ, timeMapZ, timeMoveZ, &
                                indxNeutType, iLoad, loadNume, &
                                loadPreObjectZ, loadLigrelZ, &
                                nbInputField, inputLoadFieldZ, &
                                nbFieldInGene, lpain, lchin, &
                                jvBase, resuElemZ, vectElemZ)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=4), intent(in) :: typeTher
        integer(kind=8), intent(in) :: indxNeutType, iLoad, loadNume
        character(len=*), intent(in) :: modelZ, timeMapZ, timeMoveZ
        character(len=*), intent(in) :: loadPreObjectZ, loadLigrelZ
        integer(kind=8), intent(in) :: nbInputField
        character(len=*), intent(in) :: inputLoadFieldZ(2)
        integer(kind=8), intent(in) :: nbFieldInGene
        character(len=*), intent(inout) :: lpain(LOAD_NEUT_NBMAXIN), lchin(LOAD_NEUT_NBMAXIN)
        character(len=1), intent(in) :: jvBase
        character(len=*), intent(inout) :: resuElemZ
        character(len=*), intent(in) :: vectElemZ
! ----- Local
        integer(kind=8), parameter ::  nbFieldOut = 1
        character(len=8), parameter :: lpaout(nbFieldOut) = 'PVECTTR'
        character(len=24) :: lchout(nbFieldOut)
        aster_logical :: loadExist, loadIsFunc, vectCalc
        character(len=8) :: newnom
        character(len=16) :: loadRHSOption
        character(len=24) :: loadField(2), ligrelToUse
        integer(kind=8) :: nbFieldIn
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(indxNeutType .ge. 1)
        ASSERT(indxNeutType .le. LOAD_NEUT_NBTYPE)

! ----- Detect thermal load
        call getNeumLoadType(indxNeutType, &
                             loadNume, loadPreObjectZ, &
                             nbInputField, inputLoadFieldZ, &
                             loadExist, loadIsFunc, &
                             loadField)

        vectCalc = ASTER_FALSE
        if (loadExist) then
! --------- Get option to compute
            call getRHSOption(indxNeutType, loadIsFunc, loadRHSOption)
            if (loadRHSOption .ne. "NoVector") then
                vectCalc = ASTER_TRUE
            end if
        end if

        if (vectCalc) then

! --------- Add specific input fields
            call prepVectFields(indxNeutType, typeTher, &
                                modelZ, timeMapZ, timeMoveZ, &
                                loadField, loadRHSOption, &
                                loadIsFunc, &
                                nbFieldInGene, nbFieldIn, lpain, lchin)

! --------- Get LIGREL to use
            call getLigrelToUse(indxNeutType, &
                                modelZ, loadLigrelZ, &
                                ligrelToUse)

! --------- Generate new RESU_ELEM name
            newnom = resuElemZ(10:16)
            call gcnco2(newnom)
            resuElemZ(10:16) = newnom(2:8)
            lchout(1) = resuElemZ
            call corich('E', resuElemZ, ichin_=iLoad)

! --------- Compute
            call calcul("S", loadRHSOption, ligrelToUse, &
                        nbFieldIn, lchin, lpain, &
                        nbFieldOut, lchout, lpaout, &
                        jvBase, 'OUI')

! --------- Add RESU_ELEM
            call reajre(vectElemZ, resuElemZ, jvBase)

        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getRHSOption
!
! Get name of option for vector (right-hand side)
!
! In  indxNeutType      : index of the type
! In  loadIsFunc        : flag if load is a function
! Out option            : option for RHS
!
! --------------------------------------------------------------------------------------------------
    subroutine getRHSOption(indxNeutType, loadIsFunc, option)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: indxNeutType
        aster_logical, intent(in) :: loadIsFunc
        character(len=16), intent(out) :: option
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(indxNeutType .ge. 1 .and. indxNeutType .le. LOAD_NEUT_NBTYPE)
        option = "NoVector"
        if (loadIsFunc) then
            option = therLoadVectF(indxNeutType)
        else
            option = therLoadVectR(indxNeutType)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! prepVectFields
!
! Prepare specific fields for vector
!
! In  indxNeutType      : index of the type
! In  typeTher          : type of thermics
!                         'MOVE' for moving sources
!                         'STAT' if not
! In  modelZ            ; model
! In  loadOption        : option to compute
! In  timeMap           : time (<CARTE>)
! In  timeMove          : modified time (<CARTE>) for THER_NON_LINE_MO
! In  loadField         : standard input fields
! In  loadIsFunc        : flag if load is a function
! In  nbFieldInGene     : number of input fields (generic)
! Out nbFieldIn         : number of input fields
! IO  lpain             : list of input parameters
! IO  lchin             : list of input fields
!
! --------------------------------------------------------------------------------------------------
    subroutine prepVectFields(indxNeutType, &
                              typeTher, &
                              modelZ, timeMapZ, timeMoveZ, &
                              loadFieldZ, loadOptionZ, &
                              loadIsFunc, &
                              nbFieldInGene, nbFieldIn, lpain, lchin)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: indxNeutType
        character(len=4), intent(in) :: typeTher
        character(len=*), intent(in) :: modelZ, loadOptionZ, timeMapZ, timeMoveZ
        character(len=*), intent(in) :: loadFieldZ(2)
        aster_logical, intent(in) :: loadIsFunc
        integer(kind=8), intent(in) :: nbFieldInGene
        integer(kind=8), intent(out) :: nbFieldIn
        character(len=*), intent(inout) :: lpain(LOAD_NEUT_NBMAXIN), lchin(LOAD_NEUT_NBMAXIN)
! ----- Local
        integer(kind=8) :: ier
        aster_logical :: lXfem
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(indxNeutType .ge. 1 .and. indxNeutType .le. LOAD_NEUT_NBTYPE)
        nbFieldIn = nbFieldInGene

! ----- Detect XFEM
        call exixfe(modelZ, ier)
        lXfem = ier .ne. 0

! ----- Name of input fields
        if (loadIsFunc) then
            ASSERT(therLoadParaF(indxNeutType) .ne. 'NoInput')
            nbFieldIn = nbFieldIn+1
            lpain(nbFieldIn) = therLoadParaF(indxNeutType)
            lchin(nbFieldIn) = loadFieldZ(1)
! --------- Some loads need two input fields
            if (indxNeutType .eq. LOAD_NEUT_ECHANGE) then
                nbFieldIn = nbFieldIn+1
                lpain(nbFieldIn) = 'PT_EXTF'
                lchin(nbFieldIn) = loadFieldZ(2)
            end if

        else
            ASSERT(therLoadParaR(indxNeutType) .ne. 'NoInput')
            nbFieldIn = nbFieldIn+1
            lpain(nbFieldIn) = therLoadParaR(indxNeutType)
            lchin(nbFieldIn) = loadFieldZ(1)

! --------- Some loads need two input fields
            if (indxNeutType .eq. LOAD_NEUT_ECHANGE) then
                nbFieldIn = nbFieldIn+1
                lpain(nbFieldIn) = 'PT_EXTR'
                lchin(nbFieldIn) = loadFieldZ(2)
            end if
        end if

! ----- Specific fields for XFEM
        if (lXfem) then
            call xajcin(modelZ, loadOptionZ, LOAD_NEUT_NBMAXIN, lchin, lpain, &
                        nbFieldIn)
        end if

! ----- Select time for ECHANGE_PAROI load
        nbFieldIn = nbFieldIn+1
        lpain(nbFieldIn) = 'PINSTR'
        lchin(nbFieldIn) = timeMapZ
        if (indxNeutType .eq. LOAD_NEUT_ECH_PAROI) then
            if (typeTher .eq. 'MOVE') then
                lchin(nbFieldIn) = timeMoveZ
            end if
        end if

        ASSERT(nbFieldIn .le. LOAD_NEUT_NBMAXIN)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getNeumLoadType
!
! Detect thermal load
!
! In  indxNeutType      : index of the type
! In  loadNume          : identification of load type
! In  loadPreObject     : base JEVEUX name for object
! In  nbInputField      : number of input fields given by load (for EVOL_CHAR)
! In  inputLoadFieldZ   : name of input fields given by load (for EVOL_CHAR)
! Out loadExist         : flag if load exists
! Out loadIsFunc        : flag if load is a function
! Out loadField         : standard input field
!
! --------------------------------------------------------------------------------------------------
    subroutine getNeumLoadType(indxNeutType, &
                               loadNume, loadPreObjectZ, &
                               nbInputField, inputLoadFieldZ, &
                               loadExist, loadIsFunc, &
                               loadField)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: indxNeutType, loadNume
        character(len=*), intent(in) :: loadPreObjectZ
        integer(kind=8), intent(in) :: nbInputField
        character(len=*), intent(in) :: inputLoadFieldZ(2)
        aster_logical, intent(out) :: loadExist, loadIsFunc
        character(len=24), intent(out) :: loadField(2)
! ----- Local
        integer(kind=8) :: iret
!   ------------------------------------------------------------------------------------------------
!
        loadExist = ASTER_FALSE
        loadIsFunc = ASTER_FALSE

! ----- Identify current load: get name of input fields for this load
        loadField = " "
        if (nbInputField .eq. 2) then
            ASSERT(indxNeutType .eq. LOAD_NEUT_ECHANGE)
            loadField(1) = inputLoadFieldZ(1)
            loadField(2) = inputLoadFieldZ(2)
        elseif (nbInputField .eq. 1) then
            loadField(1) = inputLoadFieldZ(1)
            loadField(2) = " "
        elseif (nbInputField .eq. 0) then
            if (indxNeutType .eq. LOAD_NEUT_ECHANGE) then
                loadField(1) = loadPreObjectZ(1:13)//therLoadField(indxNeutType)
                loadField(2) = loadPreObjectZ(1:13)//".T_EXT"
            else
                loadField(1) = loadPreObjectZ(1:13)//therLoadField(indxNeutType)
                loadField(2) = " "
            end if
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- Is this load exists ?
        iret = 0
        if (loadField(1) .ne. " ") then
            call exisd('CHAMP_GD', loadField(1), iret)
        end if
        loadExist = iret .ne. 0

! ----- For ECHANGE: two input fields are required
        if (loadField(2) .ne. " ") then
            ASSERT(indxNeutType .eq. LOAD_NEUT_ECHANGE)
            if (loadExist) then
                call exisd('CHAMP_GD', loadField(2), iret)
                loadExist = iret .ne. 0
                ASSERT(loadExist)
            end if
        end if

! ----- Detect main type
        if (loadExist) then
            if (loadNume .eq. 2) then
                ASSERT(nbInputField .eq. 0)
                loadIsFunc = ASTER_TRUE
            else if (loadNume .eq. 3) then
                ASSERT(nbInputField .eq. 0)
                loadIsFunc = ASTER_TRUE
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getLigrelToUse
!
! Get LIGREL to use
!
! In  indxNeutType      : index of the type
! In  model             : model
! In  loadLigrel        : ligrel for load
! Out ligrelToUse       : ligrel to use
!
! --------------------------------------------------------------------------------------------------
    subroutine getLigrelToUse(indxNeutType, &
                              modelZ, loadLigrelZ, &
                              ligrelToUse)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: indxNeutType
        character(len=*), intent(in) :: modelZ, loadLigrelZ
        character(len=24), intent(out) :: ligrelToUse
! ----- Locals
        character(len=24) :: modelLigrel
        integer(kind=8) :: iret
        aster_logical :: lXfem
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(indxNeutType .ge. 1 .and. indxNeutType .le. LOAD_NEUT_NBTYPE)
        ligrelToUse = " "

! ----- Detect XFEM
        call exixfe(modelZ, iret)
        lXfem = iret .ne. 0
        call dismoi('NOM_LIGREL', modelZ, 'MODELE', repk=modelLigrel)

! ----- Select LIGREL
        if (lXfem) then
            ligrelToUse = modelLigrel
        else
            if (indxNeutType .eq. LOAD_NEUT_ECH_PAROI) then
                ligrelToUse = loadLigrelZ
            else
                ligrelToUse = modelLigrel
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getFieldFromEvol
!
! Get field from EVOL_CHAR datastructure
!
! In  evolChar          : name of datastructure from AFFE_CHAR_*/EVOL_CHAR
! In  time              : time
! In  nameInEvol        : type of load in EVOL_CHAR
! In  inputLoadField    : name of datastructure for defining load
! Out exist             : flag if nameInEvol is in evol_char
!
! --------------------------------------------------------------------------------------------------
    subroutine getFieldFromEvol(evolChar, time, &
                                nameInEvolZ, inputLoadField, &
                                exist)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=8), intent(in) :: evolChar
        real(kind=8), intent(in) :: time
        character(len=*), intent(in) :: nameInEvolZ
        character(len=19), intent(in) :: inputLoadField
        aster_logical, intent(out):: exist
! ----- Local
        integer(kind=8) :: ier
        character(len=16) :: nameInEvol
        real(kind=8), parameter :: prec = 1.0d-10
        character(len=8), parameter :: crit = 'ABSOLU'
!   ------------------------------------------------------------------------------------------------
!
        exist = ASTER_FALSE
        nameInEvol = nameInEvolZ
        call rsinch(evolChar, nameInEvol, 'INST', time, inputLoadField, &
                    'EXCLU', 'EXCLU', 0, 'V', prec, crit, ier)
        if (ier .le. 2) then
            exist = ASTER_TRUE
        else if (ier .eq. 11 .or. ier .eq. 12 .or. ier .eq. 20) then
            call utmess('F', 'CHARGES8_2', sr=time, sk=nameInEvol)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compLoadEvolVect
!
! Compute loads from EVOL_CHAR keyword - Vector
!
! In  typeTher          : type of thermics
!                         'MOVE' for moving sources
!                         'STAT' if not
! In  time              : time
! In  modelZ            ; model
! In  timeMap           : time (<CARTE>)
! In  timeMove          : modified time (<CARTE>) for THER_NON_LINE_MO
! In  iLoad             : index of current load
! In  loadPreObject     : base JEVEUX name for object
! In  loadLigrel        : ligrel for load
! In  nbFieldInGene     : number of input fields (generic)
! IO  lpain             : list of input parameters
! IO  lchin             : list of input fields
! In  jvBase            : JEVEUX base to create vectors
! IO  resuElem          : name of elementary results
! In  vectElem          : name of elementary vectors
!
! --------------------------------------------------------------------------------------------------
    subroutine compLoadEvolVect(typeTher, &
                                time, modelZ, timeMapZ, timeMoveZ, &
                                iLoad, loadPreObjectZ, loadLigrelZ, &
                                nbFieldInGene, lpain, lchin, &
                                jvBase, resuElemZ, vectElemZ)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=4), intent(in) :: typeTher
        real(kind=8), intent(in) :: time
        character(len=*), intent(in) :: modelZ, timeMapZ, timeMoveZ
        integer(kind=8), intent(in) :: iLoad
        character(len=*), intent(in) :: loadPreObjectZ, loadLigrelZ
        integer(kind=8), intent(in) :: nbFieldInGene
        character(len=*), intent(inout) :: lpain(LOAD_NEUT_NBMAXIN)
        character(len=*), intent(inout) :: lchin(LOAD_NEUT_NBMAXIN)
        character(len=1), intent(in) :: jvBase
        character(len=*), intent(inout) :: resuElemZ
        character(len=*), intent(in) :: vectElemZ
! ----- Local
        character(len=24), parameter :: inputLoadField(2) = &
                                        (/'&&NTDEPR.OBJECT1        ', &
                                          '&&NTDEPR.OBJECT2        '/)
! ----- Loads in EVOl_CHAR: no function
        integer(kind=8), parameter :: loadNume = 1
        integer(kind=8) :: ier, nbField, nbInputField
        character(len=8) :: evol_char
        character(len=16) :: type_sd
        character(len=24) :: loadObjectJv
        character(len=8), pointer :: loadObject(:) => null()
        integer(kind=8) :: indxNeutType
        aster_logical :: existFluxNorm, existCoefH, existTempExt, hasLoad
!   ------------------------------------------------------------------------------------------------
!
        call jemarq()

! ----- Get object from AFFE_CHAR_THER
        loadObjectJv = loadPreObjectZ(1:13)//'.EVOL.CHAR'
        call jeexin(loadObjectJv, ier)
        if (ier .eq. 0) then
            goto 99
        end if
        call jeveuo(loadObjectJv, 'L', vk8=loadObject)
        evol_char = loadObject(1)

! ----- No load in EVOL_CHAR
        indxNeutType = LOAD_NEUT_UNKNOWN
        hasLoad = ASTER_FALSE
        nbInputField = 0

! ----- Some checks
        call dismoi('NB_CHAMP_UTI', evol_char, 'RESULTAT', repi=nbField)
        ASSERT(nbField .gt. 0)
        call gettco(evol_char, type_sd)
        ASSERT(type_sd .eq. 'EVOL_CHAR')

! ----- Get FLUX_REP
        call getFieldFromEvol(evol_char, time, &
                              "FLUN", inputLoadField(1), &
                              existFluxNorm)

        if (existFluxNorm) then
            hasLoad = ASTER_TRUE
            indxNeutType = LOAD_NEUT_FLUX_XYZ
            nbInputField = 1
        end if

! ----- Compute FLUX_REP
        if (existFluxNorm) then
            call compLoadVectType(typeTher, &
                                  modelZ, timeMapZ, timeMoveZ, &
                                  indxNeutType, iLoad, loadNume, &
                                  loadPreObjectZ, loadLigrelZ, &
                                  nbInputField, inputLoadField, &
                                  nbFieldInGene, lpain, lchin, &
                                  jvBase, resuElemZ, vectElemZ)
        end if

! ----- Get ECHANGE with COEF_H and T_EXT
        call getFieldFromEvol(evol_char, time, &
                              "COEF_H", inputLoadField(1), &
                              existCoefH)

! ----- Get ECHANGE with COEF_H and T_EXT
        call getFieldFromEvol(evol_char, time, &
                              "T_EXT", inputLoadField(2), &
                              existTempExt)

        if ((existCoefH .and. .not. existTempExt) .or. &
            (existTempExt .and. .not. existCoefH)) then
            call utmess('F', 'CHARGES8_12', sr=time)
        end if

        if (existTempExt .and. existCoefH) then
            hasLoad = ASTER_TRUE
            indxNeutType = LOAD_NEUT_ECHANGE
            nbInputField = 2
        end if

! ----- Compute ECHANGE with COEF_H and T_EXT
        if (existTempExt .and. existCoefH) then
            call compLoadVectType(typeTher, &
                                  modelZ, timeMapZ, timeMoveZ, &
                                  indxNeutType, iLoad, loadNume, &
                                  loadPreObjectZ, loadLigrelZ, &
                                  nbInputField, inputLoadField, &
                                  nbFieldInGene, lpain, lchin, &
                                  jvBase, resuElemZ, vectElemZ)
        end if
!
        if (.not. hasLoad) then
            call utmess('A', 'CHARGES8_1', sr=time)
        end if
!
99      continue
!
        call jedema()
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compLoadResi
!
! Computation of thermal loads - Residual
!
! In  l_stat            : .true. if stationnary
! In  theta             : value of coefficient for theta scheme
! In  modelZ            ; model
! In  timeMap           : time (<CARTE>)
! In  loadNume          : identification of load type
! In  loadPreObject     : base JEVEUX name for object
! In  loadLigrel        : ligrel for load
! In  nbFieldInGene     : number of input fields (generic)
! IO  lpain             : list of input parameters
! IO  lchin             : list of input fields
! In  jvBase            : JEVEUX base to create vector
! IO  resuElem          : name of elementary results
! In  vectElem          : name of elementary vectors
!
! --------------------------------------------------------------------------------------------------
    subroutine compLoadResi(l_stat, theta, &
                            modelZ, timeMapZ, &
                            loadNume, &
                            loadPreObjectZ, loadLigrelZ, &
                            nbFieldInGene, lpain, lchin, &
                            jvBase, resuElemZ, vectElemZ)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        aster_logical, intent(in) :: l_stat
        real(kind=8), intent(in) :: theta
        character(len=*), intent(in) :: modelZ, timeMapZ
        integer(kind=8), intent(in) :: loadNume
        character(len=*), intent(in) :: loadPreObjectZ, loadLigrelZ
        integer(kind=8), intent(in) :: nbFieldInGene
        character(len=*), intent(inout) :: lpain(LOAD_NEUT_NBMAXIN), lchin(LOAD_NEUT_NBMAXIN)
        character(len=1), intent(in) :: jvBase
        character(len=*), intent(inout) :: resuElemZ
        character(len=*), intent(in) :: vectElemZ
! ----- Local
        integer(kind=8), parameter :: nbInputField = 0
        character(len=24), parameter :: inputLoadField(2) = &
                                        (/'                        ', &
                                          '                        '/)
        integer(kind=8) :: indxNeutType
!   ------------------------------------------------------------------------------------------------
!
        do indxNeutType = 1, LOAD_NEUT_NBTYPE
            call compLoadResiType(l_stat, theta, &
                                  modelZ, timeMapZ, &
                                  indxNeutType, loadNume, &
                                  loadPreObjectZ, loadLigrelZ, &
                                  nbInputField, inputLoadField, &
                                  nbFieldInGene, lpain, lchin, &
                                  jvBase, resuElemZ, vectElemZ)
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compLoadResiType
!
! Computation of specific thermal load type - Residual
!
! In  l_stat            : .true. if stationnary
! In  theta             : value of coefficient for theta scheme
! In  modelZ            ; model
! In  timeMap           : time (<CARTE>)
! In  indxNeutType      : index of the type
! In  loadNume          : identification of load type
! In  loadPreObject     : base JEVEUX name for object
! In  loadLigrel        : ligrel for load
! In  nbInputField      : number of input fields given by load (for EVOL_CHAR)
! In  inputLoadFieldZ   : name of input fields given by load (for EVOL_CHAR)
! In  nbFieldInGene     : number of input fields (generic)
! IO  lpain             : list of input parameters
! IO  lchin             : list of input fields
! In  jvBase            : JEVEUX base to create vectors
! IO  resuElem          : name of elementary results
! In  vectElem          : name of elementary vectors
!
! --------------------------------------------------------------------------------------------------
    subroutine compLoadResiType(l_stat, theta, &
                                modelZ, timeMapZ, &
                                indxNeutType, loadNume, &
                                loadPreObjectZ, loadLigrelZ, &
                                nbInputField, inputLoadFieldZ, &
                                nbFieldInGene, lpain, lchin, &
                                jvBase, resuElemZ, vectElemZ)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        aster_logical, intent(in) :: l_stat
        real(kind=8), intent(in) :: theta
        integer(kind=8), intent(in) :: indxNeutType, loadNume
        character(len=*), intent(in) :: modelZ, timeMapZ
        character(len=*), intent(in) :: loadPreObjectZ, loadLigrelZ
        integer(kind=8), intent(in) :: nbInputField
        character(len=*), intent(in) :: inputLoadFieldZ(2)
        integer(kind=8), intent(in) :: nbFieldInGene
        character(len=*), intent(inout) :: lpain(LOAD_NEUT_NBMAXIN), lchin(LOAD_NEUT_NBMAXIN)
        character(len=1), intent(in) :: jvBase
        character(len=*), intent(inout) :: resuElemZ
        character(len=*), intent(in) :: vectElemZ
! ----- Local
        integer(kind=8), parameter ::  nbFieldOut = 1
        character(len=8), parameter :: lpaout(nbFieldOut) = 'PRESIDU'
        character(len=24) :: lchout(nbFieldOut)
        aster_logical :: loadExist, loadIsFunc, resiCalc
        character(len=8) :: newnom
        character(len=16) :: loadResiOption
        character(len=24) :: loadField(2), ligrelToUse
        integer(kind=8) :: nbFieldIn
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(indxNeutType .ge. 1)
        ASSERT(indxNeutType .le. LOAD_NEUT_NBTYPE)

! ----- Detect thermal load
        call getNeumLoadType(indxNeutType, &
                             loadNume, loadPreObjectZ, &
                             nbInputField, inputLoadFieldZ, &
                             loadExist, loadIsFunc, &
                             loadField)

! ----- Get option to compute
        resiCalc = ASTER_FALSE
        if (loadExist) then
            call getResiOption(indxNeutType, loadIsFunc, loadResiOption)
            if (loadResiOption .ne. "NoVector") then
                resiCalc = ASTER_TRUE
            end if
        end if

        if (resiCalc) then
! --------- Add specific input fields
            call prepResiFields(indxNeutType, &
                                timeMapZ, &
                                loadField(1), &
                                loadIsFunc, &
                                nbFieldInGene, nbFieldIn, lpain, lchin)

! --------- Get LIGREL to use
            call getLigrelToUse(indxNeutType, &
                                modelZ, loadLigrelZ, &
                                ligrelToUse)

! --------- Generate new RESU_ELEM name
            newnom = resuElemZ(10:16)
            call gcnco2(newnom)
            resuElemZ(10:16) = newnom(2:8)
            lchout(1) = resuElemZ
            call corich('E', resuElemZ, ichin_=-1)

! --------- Compute
            call calcul("S", loadResiOption, ligrelToUse, &
                        nbFieldIn, lchin, lpain, &
                        nbFieldOut, lchout, lpaout, &
                        jvBase, 'OUI')

            if (.not. l_stat) then
                call multResuElem(lchout(1), theta)
            end if

! --------- Add RESU_ELEM
            call reajre(vectElemZ, resuElemZ, jvBase)

        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getResiOption
!
! Get name of option for residual
!
! In  indxNeutType      : index of the type
! In  loadIsFunc        : flag if load is a function
! Out option            : option for residual
!
! --------------------------------------------------------------------------------------------------
    subroutine getResiOption(indxNeutType, loadIsFunc, option)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: indxNeutType
        aster_logical, intent(in) :: loadIsFunc
        character(len=16), intent(out) :: option
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(indxNeutType .ge. 1 .and. indxNeutType .le. LOAD_NEUT_NBTYPE)
        option = "NoVector"
        if (loadIsFunc) then
            option = therLoadResiF(indxNeutType)
        else
            option = therLoadResiR(indxNeutType)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! prepResiFields
!
! Prepare specific fields for residuals
!
! In  indxNeutType      : index of the type
! In  timeMap           : time (<CARTE>)
! In  loadField         : standard input field
! In  loadIsFunc        : flag if load is a function
! In  nbFieldInGene     : number of input fields (generic)
! Out nbFieldIn         : number of input fields
! IO  lpain             : list of input parameters
! IO  lchin             : list of input fields
!
! --------------------------------------------------------------------------------------------------
    subroutine prepResiFields(indxNeutType, &
                              timeMapZ, &
                              loadFieldZ, &
                              loadIsFunc, &
                              nbFieldInGene, nbFieldIn, lpain, lchin)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: indxNeutType
        character(len=*), intent(in) :: timeMapZ
        character(len=*), intent(in) :: loadFieldZ
        aster_logical, intent(in) :: loadIsFunc
        integer(kind=8), intent(in) :: nbFieldInGene
        integer(kind=8), intent(out) :: nbFieldIn
        character(len=*), intent(inout) :: lpain(LOAD_NEUT_NBMAXIN), lchin(LOAD_NEUT_NBMAXIN)
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(indxNeutType .ge. 1 .and. indxNeutType .le. LOAD_NEUT_NBTYPE)
        nbFieldIn = nbFieldInGene

! ----- Name of input fields
        if (loadIsFunc) then
            ASSERT(therResiParaF(indxNeutType) .ne. 'NoInput')
            nbFieldIn = nbFieldIn+1
            lpain(nbFieldIn) = therResiParaF(indxNeutType)
            lchin(nbFieldIn) = loadFieldZ
        else
            ASSERT(therResiParaR(indxNeutType) .ne. 'NoInput')
            nbFieldIn = nbFieldIn+1
            lpain(nbFieldIn) = therResiParaR(indxNeutType)
            lchin(nbFieldIn) = loadFieldZ
        end if

! ----- Select time for ECHANGE_PAROI load
        nbFieldIn = nbFieldIn+1
        lpain(nbFieldIn) = 'PINSTR'
        lchin(nbFieldIn) = timeMapZ

        ASSERT(nbFieldIn .le. LOAD_NEUT_NBMAXIN)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compLoadEvolResi
!
! Compute loads from EVOL_CHAR keyword - Residuals
!
! In  l_stat            : .true. if stationnary
! In  theta             : value of coefficient for theta scheme
! In  time              : time
! In  modelZ            ; model
! In  timeMap           : time (<CARTE>)
! In  loadPreObject     : base JEVEUX name for object
! In  loadLigrel        : ligrel for load
! In  nbFieldInGene     : number of input fields (generic)
! IO  lpain             : list of input parameters
! IO  lchin             : list of input fields
! In  jvBase            : JEVEUX base to create vectors
! IO  resuElem          : name of elementary results
! In  vectElem          : name of elementary vectors
!
! --------------------------------------------------------------------------------------------------
    subroutine compLoadEvolResi(l_Stat, theta, time, &
                                modelZ, timeMapZ, &
                                loadPreObjectZ, loadLigrelZ, &
                                nbFieldInGene, lpain, lchin, &
                                jvBase, resuElemZ, vectElemZ)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        aster_logical, intent(in) :: l_stat
        real(kind=8), intent(in) :: theta
        real(kind=8), intent(in) :: time
        character(len=*), intent(in) :: modelZ, timeMapZ
        character(len=*), intent(in) :: loadPreObjectZ, loadLigrelZ
        integer(kind=8), intent(in) :: nbFieldInGene
        character(len=*), intent(inout) :: lpain(LOAD_NEUT_NBMAXIN)
        character(len=*), intent(inout) :: lchin(LOAD_NEUT_NBMAXIN)
        character(len=1), intent(in) :: jvBase
        character(len=*), intent(inout) :: resuElemZ
        character(len=*), intent(in) :: vectElemZ
! ----- Local
        character(len=24), parameter :: inputLoadField(2) = &
                                        (/'&&NTDEPR.OBJECT1        ', &
                                          '&&NTDEPR.OBJECT2        '/)
! ----- Loads in EVOl_CHAR: no function
        integer(kind=8), parameter :: loadNume = 1
        integer(kind=8) :: ier, nbField, nbInputField
        character(len=8) :: evol_char
        character(len=16) :: type_sd
        character(len=24) :: loadObjectJv
        character(len=8), pointer :: loadObject(:) => null()
        integer(kind=8) :: indxNeutType
        aster_logical :: existFluxNorm, existCoefH, hasLoad
!   ------------------------------------------------------------------------------------------------
!
        call jemarq()

! ----- Get object from AFFE_CHAR_THER
        loadObjectJv = loadPreObjectZ(1:13)//'.EVOL.CHAR'
        call jeexin(loadObjectJv, ier)
        if (ier .eq. 0) then
            goto 99
        end if
        call jeveuo(loadObjectJv, 'L', vk8=loadObject)
        evol_char = loadObject(1)

! ----- No load in EVOL_CHAR
        indxNeutType = LOAD_NEUT_UNKNOWN
        hasLoad = ASTER_FALSE
        nbInputField = 0

! ----- Some checks
        call dismoi('NB_CHAMP_UTI', evol_char, 'RESULTAT', repi=nbField)
        ASSERT(nbField .gt. 0)
        call gettco(evol_char, type_sd)
        ASSERT(type_sd .eq. 'EVOL_CHAR')

! ----- Get FLUX_REP
        call getFieldFromEvol(evol_char, time, &
                              "FLUN", inputLoadField(1), &
                              existFluxNorm)

        if (existFluxNorm) then
            hasLoad = ASTER_TRUE
            indxNeutType = LOAD_NEUT_FLUX_XYZ
            nbInputField = 1
        end if

! ----- Compute FLUX_REP
        if (existFluxNorm) then
            call compLoadResiType(l_stat, theta, &
                                  modelZ, timeMapZ, &
                                  indxNeutType, loadNume, &
                                  loadPreObjectZ, loadLigrelZ, &
                                  nbInputField, inputLoadField, &
                                  nbFieldInGene, lpain, lchin, &
                                  jvBase, resuElemZ, vectElemZ)
        end if

! ----- Get ECHANGE with COEF_H
        call getFieldFromEvol(evol_char, time, &
                              "COEF_H", inputLoadField(1), &
                              existCoefH)

        if (existCoefH) then
            hasLoad = ASTER_TRUE
            indxNeutType = LOAD_NEUT_ECHANGE
            nbInputField = 1
        end if

! ----- Compute ECHANGE with COEF_H
        if (existCoefH) then
            call compLoadResiType(l_stat, theta, &
                                  modelZ, timeMapZ, &
                                  indxNeutType, loadNume, &
                                  loadPreObjectZ, loadLigrelZ, &
                                  nbInputField, inputLoadField, &
                                  nbFieldInGene, lpain, lchin, &
                                  jvBase, resuElemZ, vectElemZ)
        end if
!
        if (.not. hasLoad) then
            call utmess('A', 'CHARGES8_1', sr=time)
        end if
!
99      continue
!
        call jedema()
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compLoadMatr
!
! Computation of thermal loads - Matrix
!
! In  l_stat            : .true. if stationnary
! In  theta             : value of coefficient for theta scheme
! In  modelZ            ; model
! In  timeMap           : time (<CARTE>)
! In  loadNume          : identification of load type
! In  loadPreObject     : base JEVEUX name for object
! In  loadLigrel        : ligrel for load
! In  nbFieldInGene     : number of input fields (generic)
! IO  lpain             : list of input parameters
! IO  lchin             : list of input fields
! In  jvBase            : JEVEUX base to create vector
! IO  resuElem          : name of elementary results
! In  matrElem          : name of elementary matrices
!
! --------------------------------------------------------------------------------------------------
    subroutine compLoadMatr(l_stat, theta, &
                            modelZ, timeMapZ, &
                            loadNume, &
                            loadPreObjectZ, loadLigrelZ, &
                            nbFieldInGene, lpain, lchin, &
                            jvBase, resuElemZ, matrElemZ, &
                            lMove_, timeMapMoveZ_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        aster_logical, intent(in) :: l_stat
        real(kind=8), intent(in) :: theta
        character(len=*), intent(in) :: modelZ, timeMapZ
        integer(kind=8), intent(in) :: loadNume
        character(len=*), intent(in) :: loadPreObjectZ, loadLigrelZ
        integer(kind=8), intent(in) :: nbFieldInGene
        character(len=*), intent(inout) :: lpain(LOAD_NEUT_NBMAXIN), lchin(LOAD_NEUT_NBMAXIN)
        character(len=1), intent(in) :: jvBase
        character(len=*), intent(inout) :: resuElemZ
        character(len=*), intent(in) :: matrElemZ
        aster_logical, optional, intent(in) :: lMove_
        character(len=*), optional, intent(in) :: timeMapMoveZ_
! ----- Local
        integer(kind=8), parameter :: nbInputField = 0
        character(len=24), parameter :: inputLoadField(2) = &
                                        (/'                        ', &
                                          '                        '/)
        integer(kind=8) :: indxNeutType
        aster_logical :: lMove
        character(len=24) :: timeMapMove
!   ------------------------------------------------------------------------------------------------
!
        lMove = ASTER_FALSE
        timeMapMove = " "
        if (present(lMove_)) then
            lMove = lMove_
            timeMapMove = timeMapMoveZ_
        end if
        do indxNeutType = 1, LOAD_NEUT_NBTYPE
            call compLoadMatrType(l_stat, theta, &
                                  modelZ, timeMapZ, &
                                  indxNeutType, loadNume, &
                                  loadPreObjectZ, loadLigrelZ, &
                                  nbInputField, inputLoadField, &
                                  nbFieldInGene, lpain, lchin, &
                                  jvBase, resuElemZ, matrElemZ, &
                                  lMove, timeMapMove)
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compLoadMatrType
!
! Computation of specific thermal load type - Matrix
!
! In  l_stat            : .true. if stationnary
! In  theta             : value of coefficient for theta scheme
! In  modelZ            ; model
! In  timeMap           : time (<CARTE>)
! In  indxNeutType      : index of the type
! In  loadNume          : identification of load type
! In  loadPreObject     : base JEVEUX name for object
! In  loadLigrel        : ligrel for load
! In  nbInputField      : number of input fields given by load (for EVOL_CHAR)
! In  inputLoadFieldZ   : name of input fields given by load (for EVOL_CHAR)
! In  nbFieldInGene     : number of input fields (generic)
! IO  lpain             : list of input parameters
! IO  lchin             : list of input fields
! In  jvBase            : JEVEUX base to create matrices
! IO  resuElem          : name of elementary results
! In  matrElem          : name of elementary matrices
!
! --------------------------------------------------------------------------------------------------
    subroutine compLoadMatrType(l_stat, theta, &
                                modelZ, timeMapZ, &
                                indxNeutType, loadNume, &
                                loadPreObjectZ, loadLigrelZ, &
                                nbInputField, inputLoadFieldZ, &
                                nbFieldInGene, lpain, lchin, &
                                jvBase, resuElemZ, matrElemZ, &
                                lMove_, timeMapMoveZ_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        aster_logical, intent(in) :: l_stat
        real(kind=8), intent(in) :: theta
        integer(kind=8), intent(in) :: indxNeutType, loadNume
        character(len=*), intent(in) :: modelZ, timeMapZ
        character(len=*), intent(in) :: loadPreObjectZ, loadLigrelZ
        integer(kind=8), intent(in) :: nbInputField
        character(len=*), intent(in) :: inputLoadFieldZ(2)
        integer(kind=8), intent(in) :: nbFieldInGene
        character(len=*), intent(inout) :: lpain(LOAD_NEUT_NBMAXIN), lchin(LOAD_NEUT_NBMAXIN)
        character(len=1), intent(in) :: jvBase
        character(len=*), intent(inout) :: resuElemZ
        character(len=*), intent(in) :: matrElemZ
        aster_logical, optional, intent(in) :: lMove_
        character(len=*), optional, intent(in) :: timeMapMoveZ_
! ----- Local
        integer(kind=8), parameter ::  nbFieldOut = 1
        character(len=8), parameter :: lpaout(nbFieldOut) = 'PMATTTR'
        character(len=24) :: lchout(nbFieldOut)
        aster_logical :: loadExist, loadIsFunc, matrCalc
        character(len=8) :: newnom
        character(len=16) :: loadMatrOption
        character(len=24) :: loadField(2), ligrelToUse
        integer(kind=8) :: nbFieldIn
        aster_logical :: lMove
        character(len=24) :: timeMapMove
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(indxNeutType .ge. 1)
        ASSERT(indxNeutType .le. LOAD_NEUT_NBTYPE)
        lMove = ASTER_FALSE
        timeMapMove = " "
        if (present(lMove_)) then
            lMove = lMove_
            timeMapMove = timeMapMoveZ_
        end if

! ----- Detect thermal load
        call getNeumLoadType(indxNeutType, &
                             loadNume, loadPreObjectZ, &
                             nbInputField, inputLoadFieldZ, &
                             loadExist, loadIsFunc, &
                             loadField)
        ASSERT(nbInputField .le. 1)

! ----- Get option to compute
        matrCalc = ASTER_FALSE
        if (loadExist) then
            call getMatrOption(indxNeutType, loadIsFunc, loadMatrOption)
            if (loadMatrOption .ne. "NoMatrix") then
                matrCalc = ASTER_TRUE
            end if
        end if

        if (matrCalc) then
! --------- Add specific input fields
            call prepMatrFields(indxNeutType, &
                                timeMapZ, &
                                loadField(1), &
                                loadIsFunc, &
                                nbFieldInGene, nbFieldIn, lpain, lchin, &
                                lMove, timeMapMove)

! --------- Get LIGREL to use
            call getLigrelToUse(indxNeutType, &
                                modelZ, loadLigrelZ, &
                                ligrelToUse)

! --------- Generate new RESU_ELEM name
            newnom = resuElemZ(12:16)
            call gcnco2(newnom)
            resuElemZ(10:16) = newnom(2:8)
            lchout(1) = resuElemZ
            call corich('E', resuElemZ, ichin_=-1)

! --------- Compute
            call calcul("S", loadMatrOption, ligrelToUse, &
                        nbFieldIn, lchin, lpain, &
                        nbFieldOut, lchout, lpaout, &
                        jvBase, 'OUI')
            if (.not. l_stat) then
                call multResuElem(lchout(1), theta)
            end if

! --------- Add RESU_ELEM
            call reajre(matrElemZ, resuElemZ, jvBase)

        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getMatrOption
!
! Get name of option for matrix
!
! In  indxNeutType      : index of the type
! In  loadIsFunc        : flag if load is a function
! Out option            : option for matrix
!
! --------------------------------------------------------------------------------------------------
    subroutine getMatrOption(indxNeutType, loadIsFunc, option)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: indxNeutType
        aster_logical, intent(in) :: loadIsFunc
        character(len=16), intent(out) :: option
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(indxNeutType .ge. 1 .and. indxNeutType .le. LOAD_NEUT_NBTYPE)
        option = "NoMatrix"
        if (loadIsFunc) then
            option = therLoadMatrF(indxNeutType)
        else
            option = therLoadMatrR(indxNeutType)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! prepMatrFields
!
! Prepare specific fields for matrices
!
! In  indxNeutType      : index of the type
! In  timeMap           : time (<CARTE>)
! In  loadField         : standard input field
! In  loadIsFunc        : flag if load is a function
! In  nbFieldInGene     : number of input fields (generic)
! Out nbFieldIn         : number of input fields
! IO  lpain             : list of input parameters
! IO  lchin             : list of input fields
!
! --------------------------------------------------------------------------------------------------
    subroutine prepMatrFields(indxNeutType, &
                              timeMapZ, &
                              loadFieldZ, &
                              loadIsFunc, &
                              nbFieldInGene, nbFieldIn, lpain, lchin, &
                              lMove_, timeMapMoveZ_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: indxNeutType
        character(len=*), intent(in) :: timeMapZ
        character(len=*), intent(in) :: loadFieldZ
        aster_logical, intent(in) :: loadIsFunc
        integer(kind=8), intent(in) :: nbFieldInGene
        integer(kind=8), intent(out) :: nbFieldIn
        character(len=*), intent(inout) :: lpain(LOAD_NEUT_NBMAXIN), lchin(LOAD_NEUT_NBMAXIN)
        aster_logical, optional, intent(in) :: lMove_
        character(len=*), optional, intent(in) :: timeMapMoveZ_
! ----- Locals
        aster_logical :: lMove
        character(len=24) :: timeMapMove
!   ------------------------------------------------------------------------------------------------
!
        lMove = ASTER_FALSE
        timeMapMove = " "
        if (present(lMove_)) then
            lMove = lMove_
            timeMapMove = timeMapMoveZ_
        end if
        ASSERT(indxNeutType .ge. 1 .and. indxNeutType .le. LOAD_NEUT_NBTYPE)
        nbFieldIn = nbFieldInGene

! ----- Name of input fields
        if (loadIsFunc) then
            ASSERT(therMatrParaF(indxNeutType) .ne. 'NoInput')
            nbFieldIn = nbFieldIn+1
            lpain(nbFieldIn) = therMatrParaF(indxNeutType)
            lchin(nbFieldIn) = loadFieldZ

        else
            ASSERT(therMatrParaR(indxNeutType) .ne. 'NoInput')
            nbFieldIn = nbFieldIn+1
            lpain(nbFieldIn) = therMatrParaR(indxNeutType)
            lchin(nbFieldIn) = loadFieldZ

        end if

! ----- Select time for ECHANGE_PAROI load
        nbFieldIn = nbFieldIn+1
        lpain(nbFieldIn) = 'PINSTR'
        lchin(nbFieldIn) = timeMapZ
        if (indxNeutType .eq. LOAD_NEUT_ECH_PAROI) then
            if (lMove) then
                lchin(nbFieldIn) = timeMapMove
            end if
        end if

        ASSERT(nbFieldIn .le. LOAD_NEUT_NBMAXIN)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compLoadEvolMatr
!
! Compute loads from EVOL_CHAR keyword - Matrix
!
! In  l_stat            : .true. if stationnary
! In  theta             : value of coefficient for theta scheme
! In  time              : time
! In  modelZ            ; model
! In  timeMap           : time (<CARTE>)
! In  loadPreObject     : base JEVEUX name for object
! In  loadLigrel        : ligrel for load
! In  nbFieldInGene     : number of input fields (generic)
! IO  lpain             : list of input parameters
! IO  lchin             : list of input fields
! In  jvBase            : JEVEUX base to create matrices
! IO  resuElem          : name of elementary results
! In  vectElem          : name of elementary matrices
!
! --------------------------------------------------------------------------------------------------
    subroutine compLoadEvolMatr(l_Stat, theta, time, &
                                modelZ, timeMapZ, &
                                loadPreObjectZ, loadLigrelZ, &
                                nbFieldInGene, lpain, lchin, &
                                jvBase, resuElemZ, matrElemZ)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        aster_logical, intent(in) :: l_stat
        real(kind=8), intent(in) :: theta
        real(kind=8), intent(in) :: time
        character(len=*), intent(in) :: modelZ, timeMapZ
        character(len=*), intent(in) :: loadPreObjectZ, loadLigrelZ
        integer(kind=8), intent(in) :: nbFieldInGene
        character(len=*), intent(inout) :: lpain(LOAD_NEUT_NBMAXIN)
        character(len=*), intent(inout) :: lchin(LOAD_NEUT_NBMAXIN)
        character(len=1), intent(in) :: jvBase
        character(len=*), intent(inout) :: resuElemZ
        character(len=*), intent(in) :: matrElemZ
! ----- Local
        character(len=24), parameter :: inputLoadField = "&&NTDEPR.OBJECT1"
! ----- Loads in EVOl_CHAR: no function
        integer(kind=8), parameter :: loadNume = 1
        integer(kind=8) :: ier, nbField, nbInputField
        character(len=8) :: evol_char
        character(len=16) :: type_sd
        character(len=24) :: loadObjectJv
        character(len=8), pointer :: loadObject(:) => null()
        integer(kind=8) :: indxNeutType
        aster_logical :: existFluxNorm, existCoefH, hasLoad
!   ------------------------------------------------------------------------------------------------
!
        call jemarq()

! ----- Get object from AFFE_CHAR_THER
        loadObjectJv = loadPreObjectZ(1:13)//'.EVOL.CHAR'
        call jeexin(loadObjectJv, ier)
        if (ier .eq. 0) then
            goto 99
        end if
        call jeveuo(loadObjectJv, 'L', vk8=loadObject)
        evol_char = loadObject(1)

! ----- No load in EVOL_CHAR
        indxNeutType = LOAD_NEUT_UNKNOWN
        hasLoad = ASTER_FALSE
        nbInputField = 0

! ----- Some checks
        call dismoi('NB_CHAMP_UTI', evol_char, 'RESULTAT', repi=nbField)
        ASSERT(nbField .gt. 0)
        call gettco(evol_char, type_sd)
        ASSERT(type_sd .eq. 'EVOL_CHAR')

! ----- Get FLUX_REP
        call getFieldFromEvol(evol_char, time, &
                              "FLUN", inputLoadField, &
                              existFluxNorm)

        if (existFluxNorm) then
            hasLoad = ASTER_TRUE
            indxNeutType = LOAD_NEUT_FLUX_XYZ
            nbInputField = 1
        end if

! ----- Compute FLUX_REP
        if (existFluxNorm) then
            call compLoadMatrType(l_stat, theta, &
                                  modelZ, timeMapZ, &
                                  indxNeutType, loadNume, &
                                  loadPreObjectZ, loadLigrelZ, &
                                  nbInputField, inputLoadField, &
                                  nbFieldInGene, lpain, lchin, &
                                  jvBase, resuElemZ, matrElemZ)
        end if

! ----- Get ECHANGE with COEF_H
        call getFieldFromEvol(evol_char, time, &
                              "COEF_H", inputLoadField, &
                              existCoefH)

        if (existCoefH) then
            hasLoad = ASTER_TRUE
            indxNeutType = LOAD_NEUT_ECHANGE
            nbInputField = 1
        end if

! ----- Compute ECHANGE with COEF_H
        if (existCoefH) then
            call compLoadMatrType(l_stat, theta, &
                                  modelZ, timeMapZ, &
                                  indxNeutType, loadNume, &
                                  loadPreObjectZ, loadLigrelZ, &
                                  nbInputField, inputLoadField, &
                                  nbFieldInGene, lpain, lchin, &
                                  jvBase, resuElemZ, matrElemZ)
        end if
!
        if (.not. hasLoad) then
            call utmess('A', 'CHARGES8_1', sr=time)
        end if
!
99      continue
!
        call jedema()
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getTherNeumField
!
! Get name of input field to define thermal Neumann loads
!
! In  indxNeutType      : index of the type
! In  loadPreObject     : base JEVEUX name for object
! Out loadField         : name of input field to define load
!
! --------------------------------------------------------------------------------------------------
    subroutine getTherNeumField(indxNeutType, loadPreObjectZ, loadField)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: indxNeutType
        character(len=*), intent(in) :: loadPreObjectZ
        character(len=24), intent(out) :: loadField
!   ------------------------------------------------------------------------------------------------
!
        loadField = " "
        ASSERT(indxNeutType .ge. 1 .and. indxNeutType .le. LOAD_NEUT_NBTYPE)
        loadField = loadPreObjectZ(1:13)//therLoadField(indxNeutType)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! isTherLoadExist
!
! Detect existence of thermal load
!
! In  indxNeutType      : index of the type
! In  loadPreObject     : base JEVEUX name for object
! Out loadExist         : flag if load exists
! Out loadField         : standard input field
! In  hasInputField     : input field given by load (for EVOL_CHAR)
! In  inputLoadFieldZ   : name of input field given by load (for EVOL_CHAR)
!
! --------------------------------------------------------------------------------------------------
    subroutine isTherLoadExist(indxNeutType, loadPreObjectZ, &
                               loadExist, &
                               loadField_, hasInputField_, inputLoadFieldZ_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: indxNeutType
        character(len=*), intent(in) :: loadPreObjectZ
        aster_logical, intent(out) :: loadExist
        character(len=24), optional, intent(out) :: loadField_
        aster_logical, optional, intent(in) :: hasInputField_
        character(len=*), optional, intent(in) :: inputLoadFieldZ_
! ----- Local
        integer(kind=8) :: iret
        character(len=24) :: loadField, inputLoadField
        aster_logical :: hasInputField
!   ------------------------------------------------------------------------------------------------
!
        loadExist = ASTER_FALSE
        loadField = " "

! ----- Inputs
        hasInputField = ASTER_FALSE
        inputLoadField = " "
        if (present(hasInputField_)) then
            hasInputField = hasInputField_
        end if
        if (present(inputLoadFieldZ_)) then
            inputLoadField = inputLoadFieldZ_
        end if

! ----- Identify current load: get name of input field for this load
        if (hasInputField) then
            loadField = inputLoadField
        else
! --------- Field to detect
            if (indxNeutType .le. LOAD_NEUT_NBTYPE) then
                call getTherNeumField(indxNeutType, loadPreObjectZ, loadField)
            else
                ASSERT(ASTER_FALSE)
            end if
        end if

! ----- Is this load exists ?
        iret = 0
        if (loadField .ne. " ") then
            call exisd('CHAMP_GD', loadField, iret)
        end if

! ----- Set outputs
        loadExist = iret .ne. 0
        if (present(loadField_)) then
            loadField_ = loadField
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module loadTherCompute_module
