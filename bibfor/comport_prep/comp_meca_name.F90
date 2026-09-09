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
subroutine comp_meca_name(nbVari, nbVariMeca, &
                          l_excl, variExcl, l_kit_meta, &
                          relaComp, defoComp, kitComp, typeCpla, postIter, &
                          reguVisc, postIncr, &
                          adrsMGIS, solvBehavType, comporInfoVari)
!
    implicit none
!
#include "asterc/lccree.h"
#include "asterc/lcdiscard.h"
#include "asterc/lcinfo.h"
#include "asterc/lcvari.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/comp_meca_code.h"
#include "asterfort/comp_mfront_vname.h"
#include "asterfort/Metallurgy_type.h"
!
    integer(kind=8), intent(in) :: nbVari, nbVariMeca
    aster_logical, intent(in) :: l_excl
    character(len=16), intent(in) :: variExcl
    aster_logical, intent(in) :: l_kit_meta
    character(len=16), intent(in) :: relaComp, defoComp, kitComp(4)
    character(len=16), intent(in) :: typeCpla, postIter, reguVisc, postIncr
    character(len=16), intent(in) :: adrsMGIS
    integer(kind=8), intent(in) :: solvBehavType
    character(len=16), pointer :: comporInfoVari(:)
!
! --------------------------------------------------------------------------------------------------
!
! Preparation of comportment (mechanics)
!
! Names of internal state variables
!
! --------------------------------------------------------------------------------------------------
!
! In  nbVari           : number of internal variables
! In  nbVariMeca       : number of internal variables for mechanic
! In  l_excl           : .true. if exception case (no names for internal variables)
! In  variExcl         : name of internal variables if l_excl
! In  l_kit_meta       : .true. if metallurgy
! In  relaComp         : behaviour (RELATION keyword)
! In  defoComp         : model of strain (DEFORMATION keyword)
! In  typeCpla         : plane stress method (analytical or De Borst algorithm)
! In  kitComp          : KIT behaviour
! In  postIter         : type of post_treatment at each Newton iteration (POST_ITER keyword)
! In  reguVisc         : keyword for viscuous regularization (REGU_VISC keyword)
! In  postIncr         : type of post-treatment at end of time step (POST_INCR keyword)
! Ptr comporInfoVari   : pointer to names of internal state variables
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: metaNbVariMaxi = 30
    character(len=7) :: metaPhasName(META_MECA_NBPHASE_MAXI)
    character(len=8) :: metaRelaName(metaNbVariMaxi)
    character(len=16) :: metaGlobName(metaNbVariMaxi)
    integer(kind=8) :: idummy, idummy2, nbVariOther, iVariMeca, iVari
    character(len=16) :: compCodePy
    character(len=16) :: metaPhas, metaRela, metaGlob
    character(len=16) :: metaPhasPy, metaRelaPy, metaGlobPy
    integer(kind=8) :: nbMetaPhas, nbVariMetaRela, nbVariMetaGlob
    integer(kind=8) :: iMetaPhas, iVariMetaRela, iVariMetaGlob
!
! --------------------------------------------------------------------------------------------------
!
    if (l_excl) then
        comporInfoVari(1:nbVari) = variExcl
    else
! ----- Name of internal state variables
        if (solvBehavType .eq. SOLV_BEHAV_UMAT) then
            call comp_meca_code(relaComp, defoComp, typeCpla, kitComp, &
                                postIter, reguVisc, postIncr, &
                                compCodePy)
            nbVariOther = nbVari-nbVariMeca
            do iVariMeca = 1, nbVariMeca
                comporInfoVari(iVariMeca) = 'NoName'
            end do
            if (nbVariOther .ne. 0) then
                call lcvari(compCodePy, nbVariOther, comporInfoVari(nbVariMeca+1:nbVari))
            end if
            call lcdiscard(compCodePy)

        else if (solvBehavType .eq. SOLV_BEHAV_MGIS_OFFI .or. &
                 solvBehavType .eq. SOLV_BEHAV_MGIS_PROTO) then
            ASSERT(adrsMGIS .ne. ' ')
            call comp_meca_code(relaComp, defoComp, typeCpla, kitComp, &
                                postIter, reguVisc, postIncr, &
                                compCodePy)
            nbVariOther = nbVari-nbVariMeca
            call comp_mfront_vname(adrsMGIS, nbVariMeca, comporInfoVari)
            if (nbVariOther .ne. 0) then
                call lcvari(compCodePy, nbVariOther, comporInfoVari(nbVariMeca+1:nbVari))
            end if
            call lcdiscard(compCodePy)

        else
            if (l_kit_meta) then
! ------------- metaPhas: ACIER, ZIRC, ...
! ------------- metaRela: internal state variables (by phase)
! -------------           META_P_CINE_LINE, META_P_ISOT_LINE, META_P_ISOT_TRAC, META_V_CINE_LINE
! -------------           META_V_ISOT_LINE, META_V_ISOT_TRAC
! ------------- metaGlob: internal state variables (global)
! -------------           META_G_ISOT_*, META_G_CINE_*
                metaPhas = kitComp(1)
                metaRela = kitComp(2)
                metaGlob = kitComp(3)
                call lccree(1, metaPhas, metaPhasPy)
                call lccree(1, metaRela, metaRelaPy)
                call lccree(1, metaGlob, metaGlobPy)
                call lcinfo(metaPhasPy, idummy, nbMetaPhas, idummy2)
                call lcinfo(metaRelaPy, idummy, nbVariMetaRela, idummy2)
                call lcinfo(metaGlobPy, idummy, nbVariMetaGlob, idummy2)
                ASSERT(nbMetaPhas .le. META_MECA_NBPHASE_MAXI)
                ASSERT(nbVariMetaRela .le. metaNbVariMaxi)
                ASSERT(nbVariMetaGlob .le. metaNbVariMaxi)
                call lcvari(metaPhasPy, nbMetaPhas, metaPhasName)
                call lcvari(metaRelaPy, nbVariMetaRela, metaRelaName)
                call lcvari(metaGlobPy, nbVariMetaGlob, metaGlobName)
                iVari = 0

! ------------- Add internal state variables (by phase)
                do iMetaPhas = 1, nbMetaPhas
                    do iVariMetaRela = 1, nbVariMetaRela
                        iVari = iVari+1
                        comporInfoVari(iVari) = &
                            metaPhasName(iMetaPhas)//'#'//metaRelaName(iVariMetaRela)
                    end do
                end do

! ------------- Add internal state variables (global)
                do iVariMetaGlob = 1, nbVariMetaGlob
                    iVari = iVari+1
                    comporInfoVari(iVari) = metaGlobName(iVariMetaGlob)
                end do
                ASSERT(iVari .eq. nbVariMetaGlob+nbVariMetaRela*nbMetaPhas)
                call lcdiscard(metaPhasPy)
                call lcdiscard(metaRelaPy)
                call lcdiscard(metaGlobPy)

! ------------- Other internal state variables (GDEF_LOG, etc.)
                nbVariOther = nbVari-iVari
                if (nbVariOther .ne. 0) then
                    call comp_meca_code(relaComp, defoComp, typeCpla, kitComp, &
                                        postIter, reguVisc, postIncr, &
                                        compCodePy)
                    call lcvari(compCodePy, nbVariOther, comporInfoVari(iVari+1:nbVari))
                    call lcdiscard(compCodePy)
                end if

            else
                call comp_meca_code(relaComp, defoComp, typeCpla, kitComp, &
                                    postIter, reguVisc, postIncr, &
                                    compCodePy)
                call lcvari(compCodePy, nbVari, comporInfoVari)
                call lcdiscard(compCodePy)

            end if
        end if
    end if
!
end subroutine
