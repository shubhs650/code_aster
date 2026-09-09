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
subroutine comp_meca_read(lInitialState, prepMapCompor, model)
!
    use BehaviourPrepare_type
    implicit none
!
#include "asterc/lccree.h"
#include "asterc/lcdiscard.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/comp_meca_deflc.h"
#include "asterfort/comp_meca_incr.h"
#include "asterfort/comp_meca_l.h"
#include "asterfort/comp_meca_rkit.h"
#include "asterfort/comp_read_mesh.h"
#include "asterfort/compGetMecaPart.h"
#include "asterfort/compGetRelation.h"
#include "asterfort/dismoi.h"
#include "asterfort/getExternalBehaviourPara.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/jeveuo.h"
!
    aster_logical, intent(in) :: lInitialState
    type(BehaviourPrep_MapCompor), intent(inout) :: prepMapCompor
    character(len=8), intent(in) :: model
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of constitutive laws (mechanics)
!
! Read from command file
!
! --------------------------------------------------------------------------------------------------
!
! In  lInitialState    : .true. if initial state is defined
! IO  prepMapCompor    : datastructure to construct COMPOR map
! In  model            : model
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter:: factorKeyword = 'COMPORTEMENT'
    character(len=24), parameter :: cellAffeJv = '&&CARCREAD.LIST'
    aster_logical :: lAllCellAffe
    integer(kind=8) :: nbCellAffe
    character(len=8) :: mesh
    integer(kind=8) :: iFactorKeyword, nbFactorKeyword, iret
    character(len=16) :: defoComp, relaComp, typeCpla, multComp, typeComp, relaMeca
    character(len=16) :: postIter, defoLdc, rigiGeom, reguVisc, postIncr
    character(len=16) :: kitComp(4), answer
    character(len=16) :: relaCompPY
    character(len=19) :: modelFED
    aster_logical :: l_cristal, l_kit, lTotalStrain
    integer(kind=8), pointer :: modelCell(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    nbFactorKeyword = prepMapCompor%nb_comp
    mesh = ' '
    lTotalStrain = ASTER_FALSE

! - Pointer to list of elements in model
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelFED)
    call jeveuo(modelFED//'.TYFE', 'L', vi=modelCell)
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)

! - Read informations
    do iFactorKeyword = 1, nbFactorKeyword
! ----- Get RELATION from command file
        relaComp = 'VIDE'
        call compGetRelation(factorKeyword, iFactorKeyword, relaComp)

! ----- Detection of specific cases
        call comp_meca_l(relaComp, 'KIT', l_kit)
        call comp_meca_l(relaComp, 'CRISTAL', l_cristal)

! ----- Get DEFORMATION from command file
        defoComp = 'VIDE'
        call getvtx(factorKeyword, 'DEFORMATION', iocc=iFactorKeyword, scal=defoComp)

! ----- Get RIGI_GEOM from command file
        rigiGeom = ' '
        call getvtx(factorKeyword, 'RIGI_GEOM', iocc=iFactorKeyword, &
                    scal=rigiGeom, nbret=iret)
        if (iret .eq. 0) then
            rigiGeom = 'VIDE'
        end if

! ----- Post-treatment at each Newton iteration
        postIter = 'VIDE'
        call getvtx(factorKeyword, 'POST_ITER', iocc=iFactorKeyword, &
                    scal=postIter, nbret=iret)
        if (iret .eq. 0) then
            postIter = 'VIDE'
        end if

! ----- Viscuous regularization
        reguVisc = 'VIDE'
        call getvtx(factorKeyword, 'REGU_VISC', iocc=iFactorKeyword, scal=answer, nbret=iret)
        if (iret .eq. 1) then
            if (answer .eq. 'OUI') then
                reguVisc = 'REGU_VISC_ELAS'
            elseif (answer .eq. 'NON') then
                reguVisc = 'VIDE'
            else
                ASSERT(ASTER_FALSE)
            end if
        end if

! ----- Post-treatment at each time step
        postIncr = "VIDE"
        call getvtx(factorKeyword, 'POST_INCR', iocc=iFactorKeyword, &
                    scal=postIncr, nbret=iret)
        if (iret .eq. 0) then
            postIncr = 'VIDE'
        end if

! ----- For KIT
        kitComp = 'VIDE'
        if (l_kit) then
            call comp_meca_rkit(factorKeyword, iFactorKeyword, relaComp, kitComp)
        end if

! ----- Get mechanical part of behaviour
        relaMeca = 'VIDE'
        call compGetMecaPart(relaComp, kitComp, relaMeca)

! ----- Coding comportment (Python)
        call lccree(1, relaComp, relaCompPY)

! ----- Get multi-material *CRISTAL
        multComp = 'VIDE'
        if (l_cristal) then
            call getvid(factorKeyword, 'COMPOR', iocc=iFactorKeyword, scal=multComp)
        end if

! ----- Get affectation
        call comp_read_mesh(mesh, factorKeyword, iFactorKeyword, &
                            cellAffeJv, lAllCellAffe, nbCellAffe)

! ----- Get parameters for external programs (MFRONT/UMAT)
        typeCpla = 'VIDE'
        call getExternalBehaviourPara(mesh, modelCell, &
                                      cellAffeJv, lAllCellAffe, nbCellAffe, &
                                      relaComp, relaCompPY, relaMeca, defoComp, &
                                      factorKeyword, iFactorKeyword, &
                                      prepMapCompor%prepExte(iFactorKeyword))
        typeCpla = prepMapCompor%prepExte(iFactorKeyword)%cplaMGIS

! ----- Select type of behaviour (incremental or total)
        typeComp = 'VIDE'
        call comp_meca_incr(lInitialState, relaComp, defoComp, typeComp)

! ----- Select type of strain (mechanical or total) from catalog
        defoLdc = 'VIDE'
        call comp_meca_deflc(relaComp, defoComp, defoLdc)
        lTotalStrain = defoLdc .eq. 'TOTALE'

! ----- Discard
        call lcdiscard(relaCompPY)

! ----- Save parameters
        prepMapCompor%prepPara(iFactorKeyword)%rela_Comp = relaComp
        prepMapCompor%prepPara(iFactorKeyword)%meca_comp = relaMeca
        prepMapCompor%prepPara(iFactorKeyword)%defo_Comp = defoComp
        prepMapCompor%prepPara(iFactorKeyword)%type_comp = typeComp
        prepMapCompor%prepPara(iFactorKeyword)%type_cpla = typeCpla
        prepMapCompor%prepPara(iFactorKeyword)%kit_Comp = kitComp
        prepMapCompor%prepPara(iFactorKeyword)%mult_comp = multComp
        prepMapCompor%prepPara(iFactorKeyword)%post_iter = postIter
        prepMapCompor%prepPara(iFactorKeyword)%defo_ldc = defoLdc
        prepMapCompor%prepPara(iFactorKeyword)%rigi_geom = rigiGeom
        prepMapCompor%prepPara(iFactorKeyword)%regu_visc = reguVisc
        prepMapCompor%prepPara(iFactorKeyword)%post_incr = postIncr
        prepMapCompor%prepPara(iFactorKeyword)%lTotalStrain = lTotalStrain
    end do

    if (prepMapCompor%lDebug) then
        WRITE (6, *) "Données lues: ", nbFactorKeyword, " occurrences."
        do iFactorKeyword = 1, nbFactorKeyword
            WRITE (6, *) "- Occurrence : ", iFactorKeyword
            WRITE (6, *) "--- relaComp : ", prepMapCompor%prepPara(iFactorKeyword)%rela_Comp
            WRITE (6, *) "--- relaMeca : ", prepMapCompor%prepPara(iFactorKeyword)%meca_comp
            WRITE (6, *) "--- defoComp : ", prepMapCompor%prepPara(iFactorKeyword)%defo_Comp
            WRITE (6, *) "--- type_comp : ", prepMapCompor%prepPara(iFactorKeyword)%type_comp
            WRITE (6, *) "--- type_cpla : ", prepMapCompor%prepPara(iFactorKeyword)%type_cpla
            WRITE (6, *) "--- kitComp  : ", prepMapCompor%prepPara(iFactorKeyword)%kit_Comp
            WRITE (6, *) "--- mult_comp : ", prepMapCompor%prepPara(iFactorKeyword)%mult_comp
            WRITE (6, *) "--- post_iter : ", prepMapCompor%prepPara(iFactorKeyword)%post_iter
            WRITE (6, *) "--- defo_ldc  : ", prepMapCompor%prepPara(iFactorKeyword)%defo_ldc
            WRITE (6, *) "--- rigi_geom : ", prepMapCompor%prepPara(iFactorKeyword)%rigi_geom
            WRITE (6, *) "--- regu_visc : ", prepMapCompor%prepPara(iFactorKeyword)%regu_visc
            WRITE (6, *) "--- post_incr : ", prepMapCompor%prepPara(iFactorKeyword)%post_incr
            WRITE (6, *) "--- total strain : ", prepMapCompor%prepPara(iFactorKeyword)%lTotalStrain
        end do
    end if
!
end subroutine
