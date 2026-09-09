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
! aslint: disable=W1504
!
! ==================================================================================================
!
! Module for the management of computation of loads (mechanics)
!
! ==================================================================================================
!
module loadMecaCompute_module
! ==================================================================================================
    use HHO_precalc_module, only: hhoAddInputField
    use loadMecaCompute_type
    use coorSyst_module, only: setOrieFields
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: compLoadVect, compLoadEvolVect, compLoadMatr, compLoadEvolMatr
    public :: getApplyTypeForce, getRHSOption, getLHSOption
    public :: prepGeneralFields, prepSpecificFields, isMecaLoadExist
    public :: getMecaNeumField
    private :: getFieldFromEvol, compLoadWind, getNeumLoadType, getLigrelToUse
    private :: compLoadVectType, compLoadMatrType
! ==================================================================================================
    private
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/barych.h"
#include "asterfort/calcul.h"
#include "asterfort/chpnua.h"
#include "asterfort/cnocre.h"
#include "asterfort/codent.h"
#include "asterfort/copisd.h"
#include "asterfort/corich.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/exisd.h"
#include "asterfort/exixfe.h"
#include "asterfort/gcnco2.h"
#include "asterfort/gcncon.h"
#include "asterfort/gettco.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetc.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/mecact.h"
#include "asterfort/setStructFields.h"
#include "asterfort/megeom.h"
#include "asterfort/meharm.h"
#include "asterfort/nuachp.h"
#include "asterfort/pronua.h"
#include "asterfort/reajre.h"
#include "asterfort/rsinch.h"
#include "asterfort/utmess.h"
#include "asterfort/vtgpld.h"
#include "asterfort/xajcin.h"
#include "LoadTypes_type.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! getMecaNeumField
!
! Get name of input field to define mechanical Neumann loads
!
! In  indxNeumType      : index of the type
! In  loadPreObject     : base JEVEUX name for object
! Out loadField         : name of input field to define load
!
! --------------------------------------------------------------------------------------------------
    subroutine getMecaNeumField(indxNeumType, loadPreObjectZ, loadField)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: indxNeumType
        character(len=*), intent(in) :: loadPreObjectZ
        character(len=24), intent(out) :: loadField
!   ------------------------------------------------------------------------------------------------
!
        loadField = " "
        ASSERT(indxNeumType .ge. 1 .and. indxNeumType .le. LOAD_NEUM_NBTYPE)
        loadField = loadPreObjectZ(1:13)//mecaLoadField(indxNeumType)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compLoadEvolVect
!
! Compute loads from EVOL_CHAR keyword - Vector
!
! In  model             : name of model
! In  caraElem          : name of elementary characteristics (field)
! In  materField        : name of material characteristics (field)
! In  compor            : name of comportment definition (field)
! In  time              : time
! In  applySuiv         : flag for undead loads
! In  iLoad             : index of load in list of loads datastructure
! In  timePrev          : previous time
! In  timeCurr          : current time
! In  timeTheta         : parameter theta
! In  loadPreObject     : base JEVEUX name for object
! In  loadLigrel        : ligrel for load
! In  ligrelCalc        : LIGREL to compute
! In  nbFieldInGene     : number of input fields (before specific ones)
! IO  lpain             : list of input parameters
! IO  lchin             : list of input fields
! IO  resuElem          : name of elementary results
! In  vectElem          : name of elementary vectors
! In  dispPrev          : displacement at beginning of current time
! In  dispCumuInst      : displacement increment from beginning of current time
! In  strxPrev          : fibers information at beginning of current time
! In  viteCurr          : speed at current time
! In  acceCurr          : acceleration at current time
!
! --------------------------------------------------------------------------------------------------
    subroutine compLoadEvolVect(modelZ, caraElemZ, materFieldZ, compor, &
                                time, jvBase, &
                                applySuiv, iLoad, &
                                timePrev, timeCurr, timeTheta, &
                                loadPreObjectZ, loadLigrelZ, ligrelCalc, &
                                nbFieldInGene, lpain, lchin, &
                                resuElem, vectElem, &
                                dispPrev_, dispCumuInst_, strxPrev_, viteCurr_, acceCurr_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: modelZ, caraElemZ, materFieldZ
        character(len=24), intent(in) :: compor
        real(kind=8), intent(in) :: time
        character(len=1), intent(in) :: jvBase
        integer(kind=8), intent(in) :: iLoad
        real(kind=8), intent(in) :: timePrev, timeCurr, timeTheta
        character(len=*), intent(in) :: loadPreObjectZ, loadLigrelZ
        aster_logical, intent(in) :: applySuiv
        character(len=24), intent(in)  :: ligrelCalc
        integer(kind=8), intent(in) :: nbFieldInGene
        character(len=*), intent(inout) :: lpain(LOAD_NEUM_NBMAXIN)
        character(len=*), intent(inout) :: lchin(LOAD_NEUM_NBMAXIN)
        character(len=19), intent(inout) :: resuElem
        character(len=19), intent(in) :: vectElem
        character(len=19), optional, intent(in) :: dispPrev_, dispCumuInst_
        character(len=19), optional, intent(in) :: strxPrev_, viteCurr_, acceCurr_
! ----- Local
        character(len=19), parameter :: inputLoadField = '&&NMDEPR'
        character(len=1), parameter :: stop = 'S'
! ----- Loads in EVOl_CHAR: no function
        aster_logical, parameter :: hasInputField = ASTER_TRUE
        aster_logical, parameter :: applyPilo = ASTER_FALSE
        integer(kind=8), parameter :: loadNume = 1
        integer(kind=8) :: ier, nbField, iexist
        character(len=8) :: evol_char, newnom
        character(len=16) :: type_sd
        character(len=24) :: loadObjectJv
        character(len=8), pointer :: loadObject(:) => null()
        integer(kind=8) :: indxNeumType
        character(len=24) :: dispPrev, dispCumuInst, strxPrev, viteCurr, acceCurr
        aster_logical :: existBody2D, existBody3D, existSurf, existLine, existPressure
        aster_logical :: existNodalForce, existWind, hasLoad
        aster_logical :: lSuivLoad
!   ------------------------------------------------------------------------------------------------
!
        call jemarq()

! ----- Initializations
        dispPrev = " "
        if (present(dispPrev_)) then
            dispPrev = dispPrev_
        end if
        dispCumuInst = " "
        if (present(dispCumuInst_)) then
            dispCumuInst = dispCumuInst_
        end if
        strxPrev = " "
        if (present(strxPrev_)) then
            strxPrev = strxPrev_
        end if
        viteCurr = " "
        if (present(viteCurr_)) then
            viteCurr = viteCurr_
        end if
        acceCurr = " "
        if (present(acceCurr_)) then
            acceCurr = acceCurr_
        end if

! ----- Get object from AFFE_CHAR_MECA
        loadObjectJv = loadPreObjectZ(1:13)//'.EVOL.CHAR'
        call jeexin(loadObjectJv, ier)
        if (ier .eq. 0) then
            goto 99
        end if
        call jeveuo(loadObjectJv, 'L', vk8=loadObject)
        evol_char = loadObject(1)

! ----- No loads
        indxNeumType = LOAD_NEUM_UNKNOWN
        hasLoad = ASTER_FALSE

! ----- Some checks
        call dismoi('NB_CHAMP_UTI', evol_char, 'RESULTAT', repi=nbField)
        ASSERT(nbField .gt. 0)
        call gettco(evol_char, type_sd)
        ASSERT(type_sd .eq. 'EVOL_CHAR')

! ----- Get body force (3D)
        call getFieldFromEvol(evol_char, time, &
                              "FVOL_3D", inputLoadField, &
                              existBody3D)
        if (existBody3D) then
            indxNeumType = LOAD_NEUM_INTE_3D
            hasLoad = ASTER_TRUE
        end if

! ----- Get body force (2d)
        call getFieldFromEvol(evol_char, time, &
                              "FVOL_2D", inputLoadField, &
                              existBody2D)
        if (existBody2D) then
            indxNeumType = LOAD_NEUM_INTE_2D
            hasLoad = ASTER_TRUE
        end if

        if (existBody2D .and. existBody3D) then
            call utmess('F', 'CHARGES8_13')
        end if

! ----- Compute body forces (CHAR_MECA_FR2D2D / CHAR_MECA_FR3D3D)
        if (existBody3D .or. existBody2D) then
            call getApplyTypeForce(indxNeumType, lSuivLoad)
            if (applySuiv .and. .not. lSuivLoad) then
                call utmess('F', 'CHARGES8_14')
            end if
            call compLoadVectType(applyPilo, applySuiv, &
                                  indxNeumType, iLoad, loadNume, &
                                  modelZ, materFieldZ, compor, &
                                  strxPrev, viteCurr, acceCurr, &
                                  timePrev, timeCurr, timeTheta, &
                                  loadPreObjectZ, loadLigrelZ, &
                                  hasInputField, inputLoadField, &
                                  stop, nbFieldInGene, lpain, lchin, &
                                  ligrelCalc, jvBase, resuElem, vectElem)
        end if

! ----- Get surfacic force
        call getFieldFromEvol(evol_char, time, &
                              "FSUR_3D", inputLoadField, &
                              existSurf)
        if (existSurf) then
            indxNeumType = LOAD_NEUM_FORC_FACE
            hasLoad = ASTER_TRUE
        end if

! ----- Get lineic force
        call getFieldFromEvol(evol_char, time, &
                              "FSUR_2D", inputLoadField, &
                              existLine)
        if (existLine) then
            indxNeumType = LOAD_NEUM_FORC_WIRE
            hasLoad = ASTER_TRUE
        end if

        if (existSurf .and. existLine) then
            call utmess('F', 'CHARGES8_13')
        end if

! ----- Compute surfacic forces (CHAR_MECA_FR2D3D / CHAR_MECA_FR1D2D)
        if (existSurf .or. existLine) then
            call getApplyTypeForce(indxNeumType, lSuivLoad)
            if (applySuiv .and. .not. lSuivLoad) then
                call utmess('F', 'CHARGES8_14')
            end if
            call compLoadVectType(applyPilo, applySuiv, &
                                  indxNeumType, iLoad, loadNume, &
                                  modelZ, materFieldZ, compor, &
                                  strxPrev, viteCurr, acceCurr, &
                                  timePrev, timeCurr, timeTheta, &
                                  loadPreObjectZ, loadLigrelZ, &
                                  hasInputField, inputLoadField, &
                                  stop, nbFieldInGene, lpain, lchin, &
                                  ligrelCalc, jvBase, resuElem, vectElem)
        end if

! ----- Get pressure
        call getFieldFromEvol(evol_char, time, &
                              "PRES", inputLoadField, &
                              existPressure)
        if (existPressure) then
            indxNeumType = LOAD_NEUM_PRESSURE
            hasLoad = ASTER_TRUE
        end if

! ----- Compute pressure (CHAR_MECA_PRES_R)
        if (existPressure) then
            call getApplyTypeForce(indxNeumType, lSuivLoad)
            if (applySuiv .and. .not. lSuivLoad) then
                call utmess('F', 'CHARGES8_14')
            end if
            call compLoadVectType(applyPilo, applySuiv, &
                                  indxNeumType, iLoad, loadNume, &
                                  modelZ, materFieldZ, compor, &
                                  strxPrev, viteCurr, acceCurr, &
                                  timePrev, timeCurr, timeTheta, &
                                  loadPreObjectZ, loadLigrelZ, &
                                  hasInputField, inputLoadField, &
                                  stop, nbFieldInGene, lpain, lchin, &
                                  ligrelCalc, jvBase, resuElem, vectElem)
        end if

! ----- Get nodal force (VECT_ASSE)
        call getFieldFromEvol(evol_char, time, &
                              "FORC_NODA", inputLoadField, &
                              existNodalForce)
        if (existNodalForce) then
            indxNeumType = LOAD_NEUM_FORC_NODA
            hasLoad = ASTER_TRUE
        end if

! ----- Compute nodal force (VECT_ASSE)
        if (existNodalForce) then
            if (applySuiv) then
                call utmess('F', 'CHARGES8_14')
            end if
            newnom = resuElem(10:16)
            call gcnco2(newnom)
            resuElem(10:16) = newnom(2:8)
            call corich('E', resuElem, ichin_=iLoad)
            call copisd('CHAMP_GD', jvBase, inputLoadField, resuElem)
            call exisd('CHAMP_GD', resuElem, iexist)
            ASSERT((iexist .gt. 0) .or. (stop .eq. 'C'))
            call reajre(vectElem, resuElem, jvBase)
        end if

! ----- Get lineic forces for wind (CHAR_MECA_SR1D1D)
        call getFieldFromEvol(evol_char, time, &
                              "VITE_VENT", inputLoadField, &
                              existWind)

! ----- Compute pressure (CHAR_MECA_SR1D1D)
        if (applySuiv .and. existWind) then
            hasLoad = ASTER_TRUE
            call compLoadWind(modelZ, caraElemZ, &
                              dispPrev, dispCumuInst, strxPrev, viteCurr, &
                              iLoad, inputLoadField, ligrelCalc, &
                              resuElem, vectElem)
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
! getFieldFromEvol
!
! Get field from EVOL_CHAR datastructure
!
! In  ndim              : space dimension (2 or 3)
! In  evol_char         : name of datastructure from AFFE_CHAR_MECA for EVOL_CHAR
! In  time              : time
! In  nameInEvol        : type of load in EVOL_CHAR
! In  inputLoadField    : name of datastructure for defining load
! Out exist             : flag if nameInEvol is in evol_char
!
! --------------------------------------------------------------------------------------------------
    subroutine getFieldFromEvol(evol_char, time, &
                                nameInEvolZ, inputLoadField, &
                                exist)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=8), intent(in) :: evol_char
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
        call rsinch(evol_char, nameInEvol, 'INST', time, inputLoadField, &
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
! getApplyTypeForce
!
! Get status of application of this kind of force
!
! In  indxNeumType      : index of the type
! Out lSuivLoad         : Flag if load can been undead load type
! Out lPiloLoad         : Flag if load can been used for continuation methods
!
! --------------------------------------------------------------------------------------------------
    subroutine getApplyTypeForce(indxNeumType, lSuivLoad_, lPiloLoad_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: indxNeumType
        aster_logical, optional, intent(out) :: lSuivLoad_, lPiloLoad_
! ----- Local
        aster_logical :: lSuivLoad, lPiloLoad
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(indxNeumType .ge. 1 .and. indxNeumType .le. LOAD_NEUM_NBTYPE)
        lSuivLoad = mecaLoadSuiv(indxNeumType)
        lPiloLoad = mecaLoadPilo(indxNeumType)
        if (present(lSuivLoad_)) then
            lSuivLoad_ = lSuivLoad
        end if
        if (present(lPiloLoad_)) then
            lPiloLoad_ = lPiloLoad
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
! In  indxNeumType      : index of the type
! In  applySuiv         : flag for undead loads
! In  loadIsFunc        : flag if load is a function
! In  loadIsSigm        : flag if load is PRE_SIGM
! Out option            : option for RHS
!
! --------------------------------------------------------------------------------------------------
    subroutine getRHSOption(indxNeumType, applySuiv, loadIsFunc, loadIsSigm, option)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: indxNeumType
        aster_logical, intent(in) :: applySuiv, loadIsFunc, loadIsSigm
        character(len=16), intent(out) :: option
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(indxNeumType .ge. 1 .and. indxNeumType .le. LOAD_NEUM_NBTYPE)
        option = "NoVector"
        if (loadIsFunc) then
            if (applySuiv) then
                option = mecaLoadVectSF(indxNeumType)
            else
                option = mecaLoadVectF(indxNeumType)
            end if
        elseif (loadIsSigm) then
            option = mecaLoadVectR(indxNeumType)
            ASSERT(option .eq. 'FORC_NODA')
        else
            if (applySuiv) then
                option = mecaLoadVectSR(indxNeumType)
            else
                option = mecaLoadVectR(indxNeumType)
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compLoadWind
!
! Compute load for wind
!
! In  model             : name of model
! In  caraElem          : name of elementary characteristics (field)
! In  dispPrev          : displacement at beginning of current time
! In  dispCumuInst      : displacement increment from beginning of current time
! In  strxPrev          : fibers information at beginning of current time
! In  viteCurr          : speed at current time
! In  iLoad             : index of load in loads datastructure
! In  inputLoadField    : name of datastructure for defining load
! In  ligrelCalc        : LIGREL to compute load
! IO  resuElem          : number of elementary results
! In  vectElem          : name of elementary vectors
!
! --------------------------------------------------------------------------------------------------
    subroutine compLoadWind(model, caraElem, &
                            dispPrev, dispCumuInst, strxPrev, viteCurr, &
                            iLoad, inputLoadField, ligrelCalc, &
                            resuElem, vectElem)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=24), intent(in) :: model, caraElem
        character(len=24), intent(in) :: dispPrev, dispCumuInst, strxPrev, viteCurr
        integer(kind=8), intent(in) :: iLoad
        character(len=19), intent(in) :: inputLoadField
        character(len=24), intent(in) :: ligrelCalc
        character(len=19), intent(inout) :: resuElem
        character(len=19), intent(in) :: vectElem
! ----- Local
        integer(kind=8), parameter :: nbCmp = 3
        character(len=8), parameter :: cmpName(3) = (/'DX', 'DY', 'DZ'/)
        character(len=1), parameter :: jvBaseTemporary = "V"
        character(len=16), parameter :: option = 'CHAR_MECA_SR1D1D'
        character(len=19), parameter :: field_no_refe = '&&MNVGME.RESU_PROJE'
        character(len=19), parameter :: nuage1 = '&&NUAGE1', nuage2 = '&&NUAGE2'
        character(len=19), parameter :: method = 'NUAGE_DEG_1'
        integer(kind=8), parameter :: nbFieldInMax = 100, nbFieldOut = 1
        character(len=8) :: lpain(nbFieldInMax), lpaout(nbFieldOut)
        character(len=24) :: lchin(nbFieldInMax)
        integer(kind=8) :: nbEqua, ndim, nbno, dime, ibid, ier, nbFieldIn
        character(len=8) :: mesh_1, mesh_2, mesh_defo, answer, newnom
        integer(kind=8), pointer :: mesh1Dime(:) => null()
        character(len=19) :: nume_equa, field_no_refe1
        character(len=24) :: chgeom
        character(len=24), pointer :: fieldRefe(:) => null()
!   ------------------------------------------------------------------------------------------------
!
        mesh_defo = '.0000000'
        field_no_refe1 = '.0000000'
        newnom = '.0000000'

! ----- Mesh information
        call jelira(inputLoadField//'.VALE', 'LONMAX', ival=nbequa)
        call dismoi('NOM_MAILLA', inputLoadField, 'CHAMP', repk=mesh_1)
        call dismoi('Z_CST', mesh_1, 'MAILLAGE', repk=answer)
        ndim = 3
        if (answer .eq. 'OUI') then
            ndim = 2
        end if

! ----- Check
        call jeveuo(mesh_1//'.DIME', 'E', vi=mesh1Dime)
        nbno = mesh1Dime(1)
        dime = mesh1Dime(6)
        if (nbno*dime .ne. nbequa) then
            call utmess('F', 'CHARGES8_10')
        end if

! ----- Mesh deformation
        call gcncon('.', mesh_defo)
        call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh_2)
        call copisd('MAILLAGE', jvBaseTemporary, mesh_2, mesh_defo)
        call vtgpld('CUMU', 1.d0, mesh_2//'.COORDO', dispPrev, jvBaseTemporary, &
                    mesh_defo//'.COORDO1')
        call vtgpld('CUMU', 1.d0, mesh_defo//'.COORDO1', dispCumuInst, jvBaseTemporary, &
                    mesh_defo//'.COORDO')
        call detrsd('CHAMP_GD', mesh_defo//'.COORDO1')

! ----- Create reference field
        ibid = 0
        call cnocre(mesh_defo, 'DEPL_R', 0, [ibid], nbCmp, &
                    cmpName, [ibid], jvBaseTemporary, ' ', field_no_refe)

! ----- Create SD NUAGE
        call chpnua(ndim, inputLoadField, ' ', nuage1)
        call chpnua(ndim, field_no_refe, ' ', nuage2)

! ----- Projection on deformed mesh
        call pronua(method, nuage1, nuage2)
        call nuachp(nuage2, ' ', field_no_refe)

! ----- Set right mesh
        call dismoi("NUME_EQUA", field_no_refe, "CHAM_NO", repk=nume_equa)
        call jeveuo(nume_equa//'.REFN', 'E', vk24=fieldRefe)
        fieldRefe(1) = mesh_2

! ----- Generate relative speed field
        call jeexin(viteCurr(1:19)//'.VALE', ier)
        if (ier .gt. 0) then
            call gcncon('.', field_no_refe1)
            call copisd('CHAMP_GD', jvBaseTemporary, field_no_refe, field_no_refe1)
            call barych(field_no_refe1, viteCurr(1:19), 1.0d0, -1.0d0, &
                        field_no_refe, jvBaseTemporary)
        end if

! ----- Add input fields
        call megeom(model, chgeom)
        lpain(1) = 'PGEOMER'
        lchin(1) = chgeom
        lpain(2) = 'PVITER'
        lchin(2) = field_no_refe
        lpain(3) = 'PDEPLMR'
        lchin(3) = dispPrev
        lpain(4) = 'PDEPLPR'
        lchin(4) = dispCumuInst
        lpain(5) = 'PSTRXMR'
        lchin(5) = strxPrev
        nbFieldIn = 5

! ----- Add fields for structural elements
        call setStructFields(caraElem, nbFieldInMax, lchin, lpain, nbFieldIn)

! ----- Add fields for orientation
        call setOrieFields(nbFieldInMax, lpain, lchin, &
                           nbFieldIn, caraElem)

! ----- Output fields
        lpaout(1) = 'PVECTUR'

! ----- Generate new RESU_ELEM name
        newnom = resuElem(10:16)
        call gcnco2(newnom)
        resuElem(10:16) = newnom(2:8)
        call corich('E', resuElem, ichin_=iLoad)

! ----- Compute
        call calcul('S', option, ligrelCalc, &
                    nbFieldIn, lchin, lpain, &
                    nbFieldOut, resuElem, lpaout, &
                    jvBaseTemporary, 'OUI')
        call reajre(vectElem, resuElem, jvBaseTemporary)

! ----- Clean
        call detrsd('NUAGE', nuage1)
        call detrsd('NUAGE', nuage2)
        call jedetc(jvBaseTemporary, mesh_defo, 1)
        call detrsd('CHAMP_GD', field_no_refe1)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! prepGeneralFields
!
! Prepare input fields for computation of mechanical loads
!
! In  model             : name of model
! In  mateco            : mane of coded material
! In  caraElem          : name of elementary characteristics (field)
! In  numeHarm             : Fourier mode
! In  varcCurr          : command variable for current time
! In  dispPrev          : displacement at beginning of current time
! In  dispCumuInst      : displacement increment from beginning of current time
! Out nbFieldInGene     : number of input fields (generic for all loads)
! Out lpain             : list of input parameters
! Out lchin             : list of input fields
!
! --------------------------------------------------------------------------------------------------
    subroutine prepGeneralFields(modelZ, caraElemZ, materCodeZ, &
                                 numeHarm, &
                                 varcCurr, dispPrev, dispCumuInst, &
                                 nbFieldInGene, lpain, lchin)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=*), intent(in) :: modelZ, caraElemZ
        character(len=*), intent(in) ::  materCodeZ
        integer(kind=8), intent(in) :: numeHarm
        character(len=24), intent(in) :: varcCurr, dispPrev, dispCumuInst
        integer(kind=8), intent(out) :: nbFieldInGene
        character(len=*), intent(out) :: lpain(LOAD_NEUM_NBMAXIN), lchin(LOAD_NEUM_NBMAXIN)
! ----- Local
        integer(kind=8) :: ier
        aster_logical :: lXfem
        character(len=8) :: mesh
        character(len=24) :: chgeom, chharm
!   ------------------------------------------------------------------------------------------------
!
        nbFieldInGene = 0
        lpain = " "
        lchin = " "

! ----- Parameters from model
        call exixfe(modelZ, ier)
        lXfem = ier .ne. 0
        call dismoi('NOM_MAILLA', modelZ, 'MODELE', repk=mesh)

! ----- Get field for geometry
        call megeom(modelZ, chgeom)

! ----- Create field for Fourier
        call meharm(modelZ, numeHarm, chharm)

! ----- Standard fields
        lpain(1) = 'PGEOMER'
        lchin(1) = chgeom
        lpain(2) = 'PMATERC'
        lchin(2) = materCodeZ
        lpain(3) = 'PHARMON'
        lchin(3) = chharm
        lpain(4) = 'PVARCPR'
        lchin(4) = varcCurr
        lpain(5) = 'PDEPLMR'
        lchin(5) = dispPrev
        lpain(6) = 'PDEPLPR'
        lchin(6) = dispCumuInst
        lpain(7) = 'PABSCUR'
        lchin(7) = mesh(1:8)//'.ABSC_CURV'
        nbFieldInGene = 7

! ----- Add fields for structural elements
        call setStructFields(caraElemZ, LOAD_NEUM_NBMAXIN, lchin, lpain, nbFieldInGene)

! ----- Add fields for orientation
        call setOrieFields(LOAD_NEUM_NBMAXIN, lpain, lchin, &
                           nbFieldInGene, caraElemZ)

! ----- Specific fields for HHO
        call hhoAddInputField(modelZ, LOAD_NEUM_NBMAXIN, lchin, lpain, nbFieldInGene)

! ----- Specific fields for XFEM
        if (lXfem) then
            call xajcin(modelZ, 'CHAR_MECA_NEUM', LOAD_NEUM_NBMAXIN, lchin, lpain, nbFieldInGene)
        end if

        ASSERT(nbFieldInGene .le. LOAD_NEUM_NBMAXIN)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compLoadVect
!
! Computation of mechanical loads - Vector
!
! In  applyPilo         : flag for "PILOTAGE" loads (continuation methods)
! In  applySuiv         : flag for undead loads
! In  iLoad             : index of current load
! In  loadNume          : identification of load type
! In  model             : name of model
! In  materField        : name of material characteristics (field)
! In  compor            : name of comportment definition (field)
! In  strxPrev          : fibers information at beginning of current time
! In  viteCurr          : speed at current time
! In  acceCurr          : acceleration at current time
! In  timePrev          : previous time
! In  timeCurr          : current time
! In  timeTheta         : parameter theta
! In  loadPreObject     : base JEVEUX name for object
! In  loadLigrel        : ligrel for load
! In  stop              : CALCUL subroutine comportement
! In  nbFieldInGene     : number of input fields (generic)
! IO  lpain             : list of input parameters
! IO  lchin             : list of input fields
! In  ligrelCalc        : LIGREL to compute
! In  jvBase            : JEVEUX base to create vector
! IO  resuElem          : name of elementary results
! In  vectElem          : name of elementary vectors
!
! --------------------------------------------------------------------------------------------------
    subroutine compLoadVect(applyPilo, applySuiv, &
                            iLoad, loadNume, &
                            modelZ, materFieldZ, compor, &
                            strxPrev, viteCurr, acceCurr, &
                            timePrev, timeCurr, timeTheta, &
                            loadPreObjectZ, loadLigrelZ, &
                            stop, nbFieldInGene, lpain, lchin, &
                            ligrelCalcZ, jvBase, resuElemZ, vectElemZ)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        aster_logical, intent(in) :: applyPilo, applySuiv
        integer(kind=8), intent(in) :: iLoad, loadNume
        character(len=*), intent(in) :: modelZ, materFieldZ
        character(len=24), intent(in) :: compor, strxPrev, viteCurr, acceCurr
        real(kind=8), intent(in) :: timePrev, timeCurr, timeTheta
        character(len=*), intent(in) :: loadPreObjectZ, loadLigrelZ
        character(len=1), intent(in) :: stop
        integer(kind=8), intent(in) :: nbFieldInGene
        character(len=*), intent(inout) :: lpain(LOAD_NEUM_NBMAXIN), lchin(LOAD_NEUM_NBMAXIN)
        character(len=*), intent(in) :: ligrelCalcZ
        character(len=1), intent(in) :: jvBase
        character(len=*), intent(inout) :: resuElemZ
        character(len=*), intent(in) :: vectElemZ
! ----- Local
        aster_logical, parameter :: hasInputField = ASTER_FALSE
        character(len=24), parameter :: inputLoadField = " "
        integer(kind=8) :: indxNeumType
!   ------------------------------------------------------------------------------------------------
!
        do indxNeumType = 1, LOAD_NEUM_NBTYPE
            call compLoadVectType(applyPilo, applySuiv, &
                                  indxNeumType, iLoad, loadNume, &
                                  modelZ, materFieldZ, compor, &
                                  strxPrev, viteCurr, acceCurr, &
                                  timePrev, timeCurr, timeTheta, &
                                  loadPreObjectZ, loadLigrelZ, &
                                  hasInputField, inputLoadField, &
                                  stop, nbFieldInGene, lpain, lchin, &
                                  ligrelCalcZ, jvBase, resuElemZ, vectElemZ)
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getNeumLoadType
!
! Get type of Neumann loads
!
! In  indxNeumType      : index of the type
! In  loadNume          : identification of load type
! In  loadPreObject     : base JEVEUX name for object
! In  hasInputField     : input field given by load (for EVOL_CHAR)
! In  inputLoadFieldZ   : name of input field given by load (for EVOL_CHAR)
! Out loadExist         : flag if load exists
! Out loadIsFunc        : flag if load is a function
! Out loadIsSigm        : flag if load is PRE_SIGM
! Out loadIsVectAsse    : flag if load is VECT_ASSE
! Out loadField         : standard input field
!
! --------------------------------------------------------------------------------------------------
    subroutine getNeumLoadType(indxNeumType, &
                               loadNume, loadPreObjectZ, &
                               hasInputField, inputLoadFieldZ, &
                               loadExist, loadIsSigm, loadIsFunc, loadIsVectAsse, &
                               loadField)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: indxNeumType, loadNume
        character(len=*), intent(in) :: loadPreObjectZ
        aster_logical, intent(in) :: hasInputField
        character(len=*), intent(in) :: inputLoadFieldZ
        aster_logical, intent(out) :: loadExist, loadIsSigm, loadIsFunc, loadIsVectAsse
        character(len=24), intent(out) :: loadField
! ----- Local
        character(len=8), pointer :: vectAsse(:) => null()
        character(len=8) :: loadName
        character(len=16) :: loadCommand
!   ------------------------------------------------------------------------------------------------
!
        loadExist = ASTER_FALSE
        loadIsSigm = ASTER_FALSE
        loadIsFunc = ASTER_FALSE
        loadIsVectAsse = ASTER_FALSE

! ----- Is this load exists ?
        call isMecaLoadExist(indxNeumType, loadPreObjectZ, &
                             loadExist, loadField, &
                             hasInputField, inputLoadFieldZ)

! ----- Detect main type
        if (loadExist) then
            if (loadNume .eq. 8) then
                loadIsFunc = ASTER_TRUE
            else if (loadNume .eq. 2) then
                loadIsFunc = ASTER_TRUE
            else if (loadNume .eq. 3) then
                loadIsFunc = ASTER_TRUE
            end if

! --------- Specific
            if (loadNume .eq. 55) then
                loadIsSigm = ASTER_TRUE
            end if

! --------- Specific for VECT_ASSE
            if (indxNeumType .eq. LOAD_NEUM_VECT_ASSE) then
                loadIsVectAsse = ASTER_TRUE
                call jeveuo(loadField, 'L', vk8=vectAsse)
                loadField = vectAsse(1)
            end if

! --------- Specific for undead loads
            if (loadNume .eq. 4 .or. loadNume .eq. 9 .or. loadNume .eq. 11) then
                loadName = loadPreObjectZ(1:8)
                call dismoi('TYPE_CHARGE', loadName, 'CHARGE', repk=loadCommand)
                if (loadCommand(5:6) .eq. '_F') then
                    loadIsFunc = ASTER_TRUE
                end if
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! prepSpecificFields
!
! Prepare specific fields
!
! In  indxNeumType      : index of the type
! In  applySuiv         : flag for undead loads
! In  model             : name of model
! In  materField        : name of material characteristics (field)
! In  compor            : name of comportment definition (field)
! In  strxPrev          : fibers information at beginning of current time
! In  viteCurr          : speed at current time
! In  acceCurr          : acceleration at current time
! In  timePrev          : previous time
! In  timeCurr          : current time
! In  timeTheta         : parameter theta
! In  loadField         : standard input field
! In  loadPreObject     : base JEVEUX name for object
! In  loadIsFunc        : flag if load is a function
! In  loadIsSigm        : flag if load is PRE_SIGM
! In  loadIsVectAsse    : flag if load is VECT_ASSE
! In  nbFieldInGene     : number of input fields (generic)
! Out nbFieldIn         : number of input fields
! IO  lpain             : list of input parameters
! IO  lchin             : list of input fields
!
! --------------------------------------------------------------------------------------------------
    subroutine prepSpecificFields(indxNeumType, applySuiv, &
                                  modelZ, materFieldZ, compor, &
                                  strxPrev, viteCurr, acceCurr, &
                                  timePrev, timeCurr, timeTheta, &
                                  loadFieldZ, loadPreObjectZ, &
                                  loadIsSigm, loadIsFunc, loadIsVectAsse, &
                                  nbFieldInGene, nbFieldIn, lpain, lchin)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: indxNeumType
        aster_logical, intent(in) :: applySuiv
        character(len=*), intent(in) :: modelZ, materFieldZ
        character(len=24), intent(in) :: compor, strxPrev, viteCurr, acceCurr
        real(kind=8), intent(in) :: timePrev, timeCurr, timeTheta
        character(len=*), intent(in) :: loadFieldZ, loadPreObjectZ
        aster_logical, intent(in) :: loadIsSigm, loadIsFunc, loadIsVectAsse
        integer(kind=8), intent(in) :: nbFieldInGene
        integer(kind=8), intent(out) :: nbFieldIn
        character(len=*), intent(inout) :: lpain(LOAD_NEUM_NBMAXIN), lchin(LOAD_NEUM_NBMAXIN)
! ----- Local
        integer(kind=8) :: ig
        character(len=8) :: ng
        character(len=8), pointer :: preSigm(:) => null()
        integer(kind=8), pointer :: desc(:) => null()
        character(len=24) :: ligrel_model
        character(len=24), parameter :: chtime = "&&VECHME.CHTIME"
        integer(kind=8), parameter :: nbCmp = 3
        character(len=8), parameter :: cmpName(nbCmp) = (/ &
                                       "INST    ", "DELTAT  ", "THETA   "/)
        real(kind=8) :: cmpValue(nbCmp)
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(indxNeumType .ge. 1 .and. indxNeumType .le. LOAD_NEUM_NBTYPE)
        nbFieldIn = nbFieldInGene

! ----- Parameters from model
        ligrel_model = modelZ(1:8)//'.MODELE'

! ----- Prepare field for time
        if (applySuiv) then
            cmpValue(1) = timeCurr
        else
            cmpValue(1) = timePrev
        end if
        cmpValue(2) = timeCurr-timePrev
        cmpValue(3) = timeTheta
        call mecact('V', chtime, 'LIGREL', ligrel_model, 'INST_R  ', &
                    ncmp=nbCmp, lnomcmp=cmpName, vr=cmpValue)
        nbFieldIn = nbFieldIn+1
        lpain(nbFieldIn) = 'PINSTR'
        lchin(nbFieldIn) = chtime

! ----- Name of input fields
        if (loadIsFunc) then
            ASSERT(mecaLoadParaF(indxNeumType) .ne. 'NoInput')
            nbFieldIn = nbFieldIn+1
            lpain(nbFieldIn) = mecaLoadParaF(indxNeumType)
            lchin(nbFieldIn) = loadFieldZ(1:24)
! --------- Some loads need two input fields
            if (indxNeumType .eq. LOAD_NEUM_EFFE_FOND) then
                nbFieldIn = nbFieldIn+1
                lpain(nbFieldIn) = 'PPREFFF'
                lchin(nbFieldIn) = loadPreObjectZ(1:13)//'.PREFF'
            end if
            if (indxNeumType .eq. LOAD_NEUM_THM_ECHA) then
                nbFieldIn = nbFieldIn+1
                lpain(nbFieldIn) = 'PDEPLMR'
                lchin(nbFieldIn) = loadPreObjectZ(1:13)//'.DEPL_R'
            end if
            if (indxNeumType .eq. LOAD_NEUM_THM_ECHAH) then
                nbFieldIn = nbFieldIn+1
                lpain(nbFieldIn) = 'PDEPLMR'
                lchin(nbFieldIn) = loadPreObjectZ(1:13)//'.DEPL_R'
            end if

        elseif (loadIsVectAsse) then
! --------- No input field
            ASSERT(mecaLoadParaR(indxNeumType) .eq. 'NoInput')

        else if (loadIsSigm) then
            call jeveuo(loadPreObjectZ(1:13)//'.SIINT.VALE', 'L', vk8=preSigm)
            nbFieldIn = nbFieldIn+1
            lpain(nbFieldIn) = mecaLoadParaR(indxNeumType)
            lchin(nbFieldIn) = preSigm(1)

        else
            ASSERT(mecaLoadParaR(indxNeumType) .ne. 'NoInput')
            nbFieldIn = nbFieldIn+1
            lpain(nbFieldIn) = mecaLoadParaR(indxNeumType)
            lchin(nbFieldIn) = loadFieldZ(1:24)

! --------- Input field can been different
            if (indxNeumType .eq. LOAD_NEUM_PRE_EPSI) then
                call jeveuo(loadPreObjectZ(1:13)//'.EPSIN.DESC', 'L', vi=desc)
                ig = desc(1)
                call jenuno(jexnum('&CATA.GD.NOMGD', ig), ng)
                if (ng .eq. 'NEUT_K8') then
                    call jeveuo(loadPreObjectZ(1:13)//'.EPSIN.VALE', 'L', vk8=preSigm)
                    lchin(nbFieldIn) = preSigm(1)
                end if
            end if

! --------- Some loads need two input fields
            if (indxNeumType .eq. LOAD_NEUM_EFFE_FOND) then
                nbFieldIn = nbFieldIn+1
                lpain(nbFieldIn) = 'PPREFFR'
                lchin(nbFieldIn) = loadPreObjectZ(1:13)//'.PREFF'
            end if
            if (indxNeumType .eq. LOAD_NEUM_THM_ECHA) then
                nbFieldIn = nbFieldIn+1
                lpain(nbFieldIn) = 'PDEPLMR'
                lchin(nbFieldIn) = loadPreObjectZ(1:13)//'.DEPL_R'
            end if
            if (indxNeumType .eq. LOAD_NEUM_THM_ECHAH) then
                nbFieldIn = nbFieldIn+1
                lpain(nbFieldIn) = 'PDEPLMR'
                lchin(nbFieldIn) = loadPreObjectZ(1:13)//'.DEPL_R'
            end if
        end if

! ----- Fields for SUIV
        if (applySuiv) then
            nbFieldIn = nbFieldIn+1
            lpain(nbFieldIn) = 'PVITPLU'
            lchin(nbFieldIn) = viteCurr
            nbFieldIn = nbFieldIn+1
            lpain(nbFieldIn) = 'PACCPLU'
            lchin(nbFieldIn) = acceCurr
            nbFieldIn = nbFieldIn+1
            lpain(nbFieldIn) = 'PSTRXMR'
            lchin(nbFieldIn) = strxPrev
        end if

! ----- Behaviour: from DEFI_COMPOR or from non-linear
        if (applySuiv) then
            nbFieldIn = nbFieldIn+1
            lpain(nbFieldIn) = 'PCOMPOR'
            lchin(nbFieldIn) = compor
        else
            nbFieldIn = nbFieldIn+1
            lpain(nbFieldIn) = 'PCOMPOR'
            lchin(nbFieldIn) = materFieldZ(1:8)//'.COMPOR'
        end if

        ASSERT(nbFieldIn .le. LOAD_NEUM_NBMAXIN)
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getLigrelToUse
!
! Get LIGREL to use
!
! In  indxNeumType      : index of the type
! In  loadField         : field to describe load
! In  loadLigrel        : ligrel for load
! In  loadCalc          : ligrel to compute
! Out ligrelToUse       : ligrel to use
!
! --------------------------------------------------------------------------------------------------
    subroutine getLigrelToUse(indxNeumType, &
                              loadFieldZ, loadLigrelZ, ligrelCalcZ, &
                              ligrelToUse)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: indxNeumType
        character(len=*), intent(in) :: loadFieldZ, loadLigrelZ, ligrelCalcZ
        character(len=24), intent(out) :: ligrelToUse
! ----- Locals
        character(len=24) :: fieldType, fieldLigrel
        integer(kind=8) :: iret
        aster_logical :: l_find
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(indxNeumType .ge. 1 .and. indxNeumType .le. LOAD_NEUM_NBTYPE)
        ligrelToUse = " "

! ----- Select LIGREL
        if (indxNeumType .eq. LOAD_NEUM_FORC_NODA) then
            ligrelToUse = loadLigrelZ(1:24)
        elseif (indxNeumType .eq. LOAD_NEUM_VECT_ASSE) then
! --------- No ligrel
        else
            l_find = ASTER_FALSE
            call dismoi('TYPE_CHAMP', loadFieldZ, 'CHAMP', repk=fieldType)
            if (fieldType(1:2) == "EL") then
                call dismoi('NOM_LIGREL', loadFieldZ, 'CHAMP', repk=fieldLigrel)
                call exisd('LIGREL', fieldLigrel, iret)
                if (iret == 1) then
                    ligrelToUse = fieldLigrel
                    l_find = ASTER_TRUE
                end if
            end if
            if (.not. l_find) then
                ligrelToUse = ligrelCalcZ(1:24)
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compLoadVectType
!
! Computation of specific mechanical load type - Vector
!
! In  applyPilo         : flag for "PILOTAGE" loads (continuation methods)
! In  applySuiv         : flag for undead loads
! In  indxNeumType      : index of the type
! In  iLoad             : index of current load
! In  loadNume          : identification of load type
! In  model             : name of model
! In  materField        : name of material characteristics (field)
! In  compor            : name of comportment definition (field)
! In  strxPrev          : fibers information at beginning of current time
! In  viteCurr          : speed at current time
! In  acceCurr          : acceleration at current time
! In  timePrev          : previous time
! In  timeCurr          : current time
! In  timeTheta         : parameter theta
! In  loadPreObject     : base JEVEUX name for object
! In  loadLigrel        : ligrel for load
! In  hasInputField     : input field given by load (for EVOL_CHAR)
! In  inputLoadFieldZ   : name of input field given by load (for EVOL_CHAR)
! In  stop              : CALCUL subroutine comportement
! In  nbFieldInGene     : number of input fields (generic)
! IO  lpain             : list of input parameters
! IO  lchin             : list of input fields
! In  ligrelCalc        : LIGREL to compute
! In  jvBase            : JEVEUX base to create vectors
! IO  resuElem          : name of elementary results
! In  vectElem          : name of elementary vectors
!
! --------------------------------------------------------------------------------------------------
    subroutine compLoadVectType(applyPilo, applySuiv, &
                                indxNeumType, iLoad, loadNume, &
                                modelZ, materFieldZ, compor, &
                                strxPrev, viteCurr, acceCurr, &
                                timePrev, timeCurr, timeTheta, &
                                loadPreObjectZ, loadLigrelZ, &
                                hasInputField, inputLoadFieldZ, &
                                stop, nbFieldInGene, lpain, lchin, &
                                ligrelCalcZ, jvBase, resuElemZ, vectElemZ)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        aster_logical, intent(in) :: applyPilo, applySuiv
        integer(kind=8), intent(in) :: indxNeumType, iLoad, loadNume
        character(len=*), intent(in) :: modelZ, materFieldZ
        character(len=24), intent(in) :: compor, strxPrev, viteCurr, acceCurr
        real(kind=8), intent(in) :: timePrev, timeCurr, timeTheta
        character(len=*), intent(in) :: loadPreObjectZ, loadLigrelZ
        aster_logical, intent(in) :: hasInputField
        character(len=*), intent(in) :: inputLoadFieldZ
        character(len=1), intent(in) :: stop
        integer(kind=8), intent(in) :: nbFieldInGene
        character(len=*), intent(inout) :: lpain(LOAD_NEUM_NBMAXIN), lchin(LOAD_NEUM_NBMAXIN)
        character(len=*), intent(in) :: ligrelCalcZ
        character(len=1), intent(in) :: jvBase
        character(len=*), intent(inout) :: resuElemZ
        character(len=*), intent(in) :: vectElemZ
! ----- Local
        integer(kind=8), parameter ::  nbFieldOut = 1
        character(len=8), parameter :: lpaout = 'PVECTUR'
        aster_logical :: loadExist, loadIsSigm, loadIsFunc, loadIsVectAsse
        character(len=8) :: newnom, answer, mesh
        character(len=16) :: loadRHSOption
        character(len=24) :: loadField, ligrelToUse
        integer(kind=8) :: iexist, nbFieldIn
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(indxNeumType .ge. 1)
        ASSERT(indxNeumType .le. LOAD_NEUM_NBTYPE)

! ----- Detect mechanical load
        call getNeumLoadType(indxNeumType, &
                             loadNume, loadPreObjectZ, &
                             hasInputField, inputLoadFieldZ, &
                             loadExist, loadIsSigm, loadIsFunc, loadIsVectAsse, &
                             loadField)

        if (loadExist) then
! --------- Checks
            if (applySuiv) then
                if (.not. mecaLoadSuiv(indxNeumType)) then
                    call utmess('F', 'CHARGES8_3', sk=mecaLoadKeyword(indxNeumType))
                end if
            end if
            if (applyPilo) then
                if (.not. mecaLoadPilo(indxNeumType)) then
                    call utmess('F', 'CHARGES8_4', sk=mecaLoadKeyword(indxNeumType))
                end if
            end if

! --------- Get option to compute
            call getRHSOption(indxNeumType, applySuiv, loadIsFunc, loadIsSigm, &
                              loadRHSOption)

! --------- Add specific input fields
            call prepSpecificFields(indxNeumType, applySuiv, &
                                    modelZ, materFieldZ, compor, &
                                    strxPrev, viteCurr, acceCurr, &
                                    timePrev, timeCurr, timeTheta, &
                                    loadField, loadPreObjectZ, &
                                    loadIsSigm, loadIsFunc, loadIsVectAsse, &
                                    nbFieldInGene, nbFieldIn, lpain, lchin)

! --------- Get LIGREL to use
            call getLigrelToUse(indxNeumType, &
                                loadField, loadLigrelZ, ligrelCalcZ, &
                                ligrelToUse)

! --------- Generate new RESU_ELEM name
            newnom = resuElemZ(10:16)
            call gcnco2(newnom)
            resuElemZ(10:16) = newnom(2:8)
            call corich('E', resuElemZ, ichin_=iLoad)

            if (loadRHSOption .eq. 'Copy_Load') then
                ASSERT(loadIsVectAsse)
                call copisd('CHAMP_GD', jvBase, loadField, resuElemZ)
                call exisd('CHAMP_GD', resuElemZ, iexist)
                ASSERT(iexist .gt. 0)
                call reajre(vectElemZ, resuElemZ, jvBase)
            else
                call calcul(stop, loadRHSOption, ligrelToUse, &
                            nbFieldIn, lchin, lpain, &
                            nbFieldOut, resuElemZ, lpaout, jvBase, &
                            'OUI')

! ------------- Copying output field
                call exisd('CHAMP_GD', resuElemZ, iexist)
                call dismoi('NOM_MAILLA', ligrelToUse, 'LIGREL', repk=mesh)
                call dismoi('PARALLEL_MESH', mesh, 'MAILLAGE', repk=answer)
                if (answer .eq. 'NON') then
                    ASSERT((iexist .gt. 0) .or. (stop .eq. 'C'))
                end if
                call reajre(vectElemZ, resuElemZ, jvBase)
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compLoadMatr
!
! Computation of mechanical loads - Matrix
!
! In  applyPilo         : flag for "PILOTAGE" loads (continuation methods)
! In  applySuiv         : flag for undead loads
! In  iLoad             : index of current load
! IO  iMatr             : index of current matrix
! In  loadNume          : identification of load type
! In  model             : name of model
! In  materField        : name of material characteristics (field)
! In  compor            : name of comportment definition (field)
! In  strxPrev          : fibers information at beginning of current time
! In  viteCurr          : speed at current time
! In  acceCurr          : acceleration at current time
! In  timePrev          : previous time
! In  timeCurr          : current time
! In  timeTheta         : parameter theta
! In  loadLigrel        : ligrel for load
! In  loadPreObject     : base JEVEUX name for object
! In  stop              : CALCUL subroutine comportement
! IO  nbFieldIn         : number of input fields
! IO  lpain             : list of input parameters
! IO  lchin             : list of input fields
! In  ligrelCalc        : LIGREL to compute
! In  jvBase            : JEVEUX base to create matrix
! In  matrElemZ         : name of elementary matrix
!
! --------------------------------------------------------------------------------------------------
    subroutine compLoadMatr(applyPilo, applySuiv, &
                            iLoad, iMatr, loadNume, &
                            modelZ, materFieldZ, compor, &
                            strxPrev, viteCurr, acceCurr, &
                            timePrev, timeCurr, timeTheta, &
                            loadLigrelZ, loadPreObjectZ, &
                            stop, nbFieldIn, lpain, lchin, &
                            ligrelCalcZ, jvBase, matrElemZ)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        aster_logical, intent(in) :: applyPilo, applySuiv
        integer(kind=8), intent(in) :: iLoad, loadNume
        integer(kind=8), intent(inout) :: iMatr
        character(len=*), intent(in) :: modelZ, materFieldZ
        character(len=24), intent(in) :: compor, strxPrev, viteCurr, acceCurr
        real(kind=8), intent(in) :: timePrev, timeCurr, timeTheta
        character(len=*), intent(in) :: loadLigrelZ, loadPreObjectZ
        character(len=1), intent(in) :: stop
        integer(kind=8), intent(inout) :: nbFieldIn
        character(len=*), intent(inout) :: lpain(LOAD_NEUM_NBMAXIN), lchin(LOAD_NEUM_NBMAXIN)
        character(len=*), intent(in) :: ligrelCalcZ
        character(len=1), intent(in) :: jvBase
        character(len=*), intent(in) :: matrElemZ
! ----- Local
        aster_logical, parameter :: hasInputField = ASTER_FALSE
        character(len=24), parameter :: inputLoadField = " "
        integer(kind=8) :: indxNeumType
!   ------------------------------------------------------------------------------------------------
!
        do indxNeumType = 1, LOAD_NEUM_NBTYPE
            call compLoadMatrType(applyPilo, applySuiv, &
                                  indxNeumType, iLoad, iMatr, loadNume, &
                                  modelZ, materFieldZ, compor, &
                                  strxPrev, viteCurr, acceCurr, &
                                  timePrev, timeCurr, timeTheta, &
                                  loadLigrelZ, loadPreObjectZ, &
                                  hasInputField, inputLoadField, &
                                  stop, nbFieldIn, lpain, lchin, &
                                  ligrelCalcZ, jvBase, matrElemZ)
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compLoadMatrType
!
! Computation of specific mechanical load type - Matrix
!
! In  applyPilo         : flag for "PILOTAGE" loads (continuation methods)
! In  applySuiv         : flag for undead loads
! In  indxNeumType      : index of the type
! In  iLoad             : index of current load
! IO  iMatr             : index of current matrix
! In  loadNume          : identification of load type
! In  loadLigrel        : ligrel for load
! In  loadPreObject     : base JEVEUX name for object
! In  hasInputField     : input field given by load (for EVOL_CHAR)
! In  inputLoadFieldZ   : name of input field given by load (for EVOL_CHAR)
! In  stop              : CALCUL subroutine comportement
! IO  nbFieldInGene     : number of input fields (generic)
! IO  lpain             : list of input parameters
! IO  lchin             : list of input fields
! In  ligrelCalc        : LIGREL to compute
! In  jvBase            : JEVEUX base to create matrix
! In  matrElem          : name of elementary matrix
!
! --------------------------------------------------------------------------------------------------
    subroutine compLoadMatrType(applyPilo, applySuiv, &
                                indxNeumType, iLoad, iMatr, loadNume, &
                                modelZ, materFieldZ, compor, &
                                strxPrev, viteCurr, acceCurr, &
                                timePrev, timeCurr, timeTheta, &
                                loadLigrelZ, loadPreObjectZ, &
                                hasInputField, inputLoadFieldZ, &
                                stop, nbFieldInGene, lpain, lchin, &
                                ligrelCalcZ, jvBase, matrElemZ)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        aster_logical, intent(in) :: applyPilo, applySuiv
        integer(kind=8), intent(in) :: indxNeumType, iLoad, loadNume
        integer(kind=8), intent(inout) :: iMatr
        character(len=*), intent(in) :: modelZ, materFieldZ
        character(len=24), intent(in) :: compor, strxPrev, viteCurr, acceCurr
        real(kind=8), intent(in) :: timePrev, timeCurr, timeTheta
        character(len=*), intent(in) :: loadLigrelZ, loadPreObjectZ
        aster_logical, intent(in) :: hasInputField
        character(len=*), intent(in) :: inputLoadFieldZ
        character(len=1), intent(in) :: stop
        integer(kind=8), intent(in) :: nbFieldInGene
        character(len=*), intent(inout) :: lpain(LOAD_NEUM_NBMAXIN), lchin(LOAD_NEUM_NBMAXIN)
        character(len=*), intent(in) :: ligrelCalcZ
        character(len=1), intent(in) :: jvBase
        character(len=*), intent(in) :: matrElemZ
! ----- Local
        integer(kind=8), parameter :: nbFieldOut = 1
        character(len=8) :: lpaout
        character(len=24) :: lchout
        aster_logical :: loadExist, loadIsSigm, loadIsFunc, loadIsVectAsse
        character(len=16) :: loadLHSOption
        character(len=24) :: loadField, ligrelToUse
        character(len=8) :: loadMatrType
        character(len=24), pointer :: listResuElem(:) => null()
        integer(kind=8) :: nbFieldIn
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(indxNeumType .ge. 1)
        ASSERT(indxNeumType .le. LOAD_NEUM_NBTYPE)
        lchout = '&&MECGME.00000000'

! ----- Detect mechanical load
        call getNeumLoadType(indxNeumType, &
                             loadNume, loadPreObjectZ, &
                             hasInputField, inputLoadFieldZ, &
                             loadExist, loadIsSigm, loadIsFunc, loadIsVectAsse, &
                             loadField)

        if (loadExist) then
! --------- Checks
            if (.not. mecaLoadSuiv(indxNeumType)) then
                call utmess('F', 'CHARGES8_3', sk=mecaLoadKeyword(indxNeumType))
            end if
            if (applyPilo) then
                if (.not. mecaLoadPilo(indxNeumType)) then
                    call utmess('F', 'CHARGES8_4', sk=mecaLoadKeyword(indxNeumType))
                end if
            end if

! --------- Get option to compute
            call getLHSOption(indxNeumType, loadIsFunc, loadIsSigm, &
                              loadLHSOption, loadMatrType)

            if (loadLHSOption .ne. 'NoMatrix') then
                ASSERT(.not. loadIsVectAsse)
                ASSERT(.not. loadIsSigm)

! ------------- Add specific input fields
                call prepSpecificFields(indxNeumType, applySuiv, &
                                        modelZ, materFieldZ, compor, &
                                        strxPrev, viteCurr, acceCurr, &
                                        timePrev, timeCurr, timeTheta, &
                                        loadField, loadPreObjectZ, &
                                        loadIsSigm, loadIsFunc, loadIsVectAsse, &
                                        nbFieldInGene, nbFieldIn, lpain, lchin)

! ------------- Get LIGREL to use
                call getLigrelToUse(indxNeumType, &
                                    loadField, loadLigrelZ, ligrelCalcZ, &
                                    ligrelToUse)

! ------------- Generate new RESU_ELEM name
                if (iMatr .ge. 0) then
                    lpaout = loadMatrType
                    lchout(10:10) = 'G'
                    iMatr = iMatr+1
                    call codent(iload, 'D0', lchout(7:8))
                    call codent(iMatr, 'D0', lchout(12:14))
                end if

! ------------- Get old RESU_ELEM
                if (iMatr .lt. 0) then
                    lpaout = loadMatrType
                    call jeveuo(matrElemZ(1:19)//'.RELR', 'L', vk24=listResuElem)
                    lchout = listResuElem(abs(iMatr))
                end if

! ------------- Computation
                call calcul(stop, loadLHSOption, ligrelToUse, &
                            nbFieldIn, lchin, lpain, &
                            nbFieldOut, lchout, lpaout, jvBase, &
                            'OUI')

! ------------- Copying output resu Elem
                if (iMatr .ge. 0) then
                    call reajre(matrElemZ, lchout, 'V')
                end if

            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getLHSOption
!
! Get name of option for matrix (left-hand side)
!
! In  indxNeumType      : index of the type
! In  loadIsFunc        : flag if load is a function
! In  loadIsSigm        : flag if load is PRE_SIGM
! Out option            : option for LHS
! Our loadMatrType      : type of matrix (symmetric or not)
!
! --------------------------------------------------------------------------------------------------
    subroutine getLHSOption(indxNeumType, loadIsFunc, loadIsSigm, &
                            option, loadMatrType)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: indxNeumType
        aster_logical, intent(in) :: loadIsFunc, loadIsSigm
        character(len=16), intent(out) :: option
        character(len=8), intent(out) :: loadMatrType
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(indxNeumType .ge. 1 .and. indxNeumType .le. LOAD_NEUM_NBTYPE)
        option = "NoMatrix"
        loadMatrType = mecaLoadParaM(indxNeumType)
        if (loadIsFunc) then
            option = mecaLoadMatrF(indxNeumType)
        elseif (loadIsSigm) then
            ASSERT(indxNeumType .eq. LOAD_NEUM_PRE_SIGM)
            option = mecaLoadMatrF(LOAD_NEUM_PRE_SIGM)
        else
            option = mecaLoadMatrR(indxNeumType)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compLoadEvolMatr
!
! Compute loads from EVOL_CHAR keyword - Matrix
!
! In  time              : current time
! In  iLoad             : index of current load
! IO  iMatr             : index of current matrix
! In  applySuiv         : flag for undead loads
! In  model             : name of model
! In  materField        : name of material characteristics (field)
! In  compor            : name of comportment definition (field)
! In  strxPrev          : fibers information at beginning of current time
! In  viteCurr          : speed at current time
! In  acceCurr          : acceleration at current time
! In  timePrev          : previous time
! In  timeCurr          : current time
! In  timeTheta         : parameter theta
! In  loadPreObject     : base JEVEUX name for object
! In  loadLigrel        : ligrel for load
! In  ligrelCalc        : LIGREL to compute
! IO  nbFieldIn         : number of input fields (before specific ones)
! IO  lpain             : list of input parameters
! IO  lchin             : list of input fields
! In  matrElem          : name of elementary matrix
!
! --------------------------------------------------------------------------------------------------
    subroutine compLoadEvolMatr(time, jvBase, &
                                applySuiv, iLoad, iMatr, &
                                modelZ, materFieldZ, compor, &
                                strxPrev, viteCurr, acceCurr, &
                                timePrev, timeCurr, timeTheta, &
                                loadPreObjectZ, loadLigrelZ, ligrelCalc, &
                                nbFieldIn, lpain, lchin, &
                                matrElemZ)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        real(kind=8), intent(in) :: time
        character(len=1), intent(in) :: jvBase
        integer(kind=8), intent(in) :: iLoad
        integer(kind=8), intent(inout) :: iMatr
        character(len=*), intent(in) :: modelZ, materFieldZ
        character(len=24), intent(in) :: compor, strxPrev, viteCurr, acceCurr
        real(kind=8), intent(in) :: timePrev, timeCurr, timeTheta
        character(len=*), intent(in) :: loadPreObjectZ, loadLigrelZ
        aster_logical, intent(in) :: applySuiv
        character(len=24), intent(in)  :: ligrelCalc
        integer(kind=8), intent(inout) :: nbFieldIn
        character(len=*), intent(inout) :: lpain(LOAD_NEUM_NBMAXIN)
        character(len=*), intent(inout) :: lchin(LOAD_NEUM_NBMAXIN)
        character(len=*), intent(in) :: matrElemZ
! ----- Local
        character(len=19), parameter :: inputLoadField = '&&NMDEPR'
        character(len=1), parameter :: stop = 'S'
! ----- Loads in EVOl_CHAR: no function
        aster_logical, parameter :: hasInputField = ASTER_TRUE
        aster_logical, parameter :: applyPilo = ASTER_FALSE
        integer(kind=8), parameter :: loadNume = 1
        integer(kind=8) :: ier, nbField
        character(len=8) :: evol_char
        character(len=16) :: type_sd
        character(len=24) :: loadObjectJv
        character(len=8), pointer :: loadObject(:) => null()
        integer(kind=8) :: indxNeumType
        aster_logical :: existPressure, lSuivLoad
!   ------------------------------------------------------------------------------------------------
!
        call jemarq()

! ----- Get object from AFFE_CHAR_MECA
        loadObjectJv = loadPreObjectZ(1:13)//'.EVOL.CHAR'
        call jeexin(loadObjectJv, ier)
        if (ier .eq. 0) then
            goto 99
        end if
        call jeveuo(loadObjectJv, 'L', vk8=loadObject)
        evol_char = loadObject(1)
        indxNeumType = LOAD_NEUM_UNKNOWN

! ----- Some checks
        call dismoi('NB_CHAMP_UTI', evol_char, 'RESULTAT', repi=nbField)
        ASSERT(nbField .gt. 0)
        call gettco(evol_char, type_sd)
        ASSERT(type_sd .eq. 'EVOL_CHAR')

! ----- Get pressure
        call getFieldFromEvol(evol_char, time, &
                              "PRES", inputLoadField, &
                              existPressure)
        if (existPressure) then
            indxNeumType = LOAD_NEUM_PRESSURE
        end if

! ----- Compute pressure (CHAR_MECA_PRES_R)
        if (existPressure) then
            call getApplyTypeForce(indxNeumType, lSuivLoad)
            if (applySuiv .and. .not. lSuivLoad) then
                call utmess('F', 'CHARGES8_14')
            end if
            call compLoadMatrType(applyPilo, applySuiv, &
                                  indxNeumType, iLoad, iMatr, loadNume, &
                                  modelZ, materFieldZ, compor, &
                                  strxPrev, viteCurr, acceCurr, &
                                  timePrev, timeCurr, timeTheta, &
                                  loadLigrelZ, loadPreObjectZ, &
                                  hasInputField, inputLoadField, &
                                  stop, nbFieldIn, lpain, lchin, &
                                  ligrelCalc, jvBase, matrElemZ)
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
! isMecaLoadExist
!
! Detect existence of mechanical load
!
! In  indxNeumType      : index of the type
! In  loadPreObject     : base JEVEUX name for object
! Out loadExist         : flag if load exists
! Out loadField         : standard input field
! In  hasInputField     : input field given by load (for EVOL_CHAR)
! In  inputLoadFieldZ   : name of input field given by load (for EVOL_CHAR)
!
! --------------------------------------------------------------------------------------------------
    subroutine isMecaLoadExist(indxNeumType, loadPreObjectZ, &
                               loadExist, &
                               loadField_, hasInputField_, inputLoadFieldZ_)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        integer(kind=8), intent(in) :: indxNeumType
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
            if (indxNeumType .le. LOAD_NEUM_NBTYPE) then
                call getMecaNeumField(indxNeumType, loadPreObjectZ, loadField)
            elseif (indxNeumType .eq. LOAD_NEUM_PWAVE) then
                loadField = loadPreObjectZ(1:13)//'.ONDPL'
            elseif (indxNeumType .eq. LOAD_NEUM_VITE_FACE) then
                loadField = loadPreObjectZ(1:13)//'.VFACE'
            elseif (indxNeumType .eq. LOAD_NEUM_WAVE) then
                loadField = loadPreObjectZ(1:13)//'.ONDE'
            elseif (indxNeumType .eq. LOAD_CPL) then
                loadField = loadPreObjectZ(1:13)//'.PAIRS'
            else
                ASSERT(ASTER_FALSE)
            end if
        end if

! ----- Is this load exists ?
        iret = 0
        if (loadField .ne. " ") then
            if (indxNeumType .eq. LOAD_NEUM_VECT_ASSE) then
                call jeexin(loadField, iret)
            else
                call exisd('CHAMP_GD', loadField, iret)
            end if
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
end module loadMecaCompute_module
