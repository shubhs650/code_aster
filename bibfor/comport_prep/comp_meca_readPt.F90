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
subroutine comp_meca_readPt(lInitialState, prepMapCompor)
!
    use BehaviourPrepare_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/comp_meca_deflc.h"
#include "asterfort/comp_meca_incr.h"
#include "asterfort/comp_meca_l.h"
#include "asterfort/comp_meca_rkit.h"
#include "asterfort/compGetMecaPart.h"
#include "asterfort/compGetRelation.h"
#include "asterfort/getExternalBehaviourParaPt.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
!
    aster_logical, intent(in) :: lInitialState
    type(BehaviourPrep_MapCompor), intent(inout) :: prepMapCompor
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
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter:: factorKeyword = 'COMPORTEMENT'
    integer(kind=8) :: iFactorKeyword, nbFactorKeyword, iret
    character(len=16) :: defoComp, relaComp, typeCpla, multComp, typeComp, relaMeca
    character(len=16) :: postIter, defoLdc, rigiGeom, reguVisc, postIncr
    character(len=16) :: kitComp(4), answer
    aster_logical :: l_cristal, l_kit, lTotalStrain
!
! --------------------------------------------------------------------------------------------------
!
    nbFactorKeyword = prepMapCompor%nb_comp
    lTotalStrain = ASTER_FALSE

! - Read informations
    do iFactorKeyword = 1, nbFactorKeyword
! ----- Get RELATION from command file
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

! ----- Get multi-material *CRISTAL
        multComp = 'VIDE'
        if (l_cristal) then
            call getvid(factorKeyword, 'COMPOR', iocc=iFactorKeyword, scal=multComp)
        end if

! ----- Get parameters for external programs (MFRONT/UMAT)
        typeCpla = 'VIDE'
        call getExternalBehaviourParaPt(relaMeca, defoComp, &
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

! ----- Save parameters
        prepMapCompor%prepPara(iFactorKeyword)%rela_comp = relaComp
        prepMapCompor%prepPara(iFactorKeyword)%meca_comp = relaMeca
        prepMapCompor%prepPara(iFactorKeyword)%defo_comp = defocomp
        prepMapCompor%prepPara(iFactorKeyword)%type_comp = typeComp
        prepMapCompor%prepPara(iFactorKeyword)%type_cpla = typeCpla
        prepMapCompor%prepPara(iFactorKeyword)%kit_comp = kitComp
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
            WRITE (6, *) "--- rela_comp : ", prepMapCompor%prepPara(iFactorKeyword)%rela_comp
            WRITE (6, *) "--- meca_comp : ", prepMapCompor%prepPara(iFactorKeyword)%meca_comp
            WRITE (6, *) "--- defo_comp : ", prepMapCompor%prepPara(iFactorKeyword)%defo_comp
            WRITE (6, *) "--- type_comp : ", prepMapCompor%prepPara(iFactorKeyword)%type_comp
            WRITE (6, *) "--- type_cpla : ", prepMapCompor%prepPara(iFactorKeyword)%type_cpla
            WRITE (6, *) "--- kit_comp  : ", prepMapCompor%prepPara(iFactorKeyword)%kit_comp
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
