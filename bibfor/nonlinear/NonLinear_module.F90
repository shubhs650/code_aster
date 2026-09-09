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
module NonLinear_module
! ==================================================================================================
    use NonLin_Datastructure_type
    use NonLinearDyna_type
    use Rom_Datastructure_type
    use ldlt_xp_data_module
    use NonLinearElem_module, only: elemDiri, elemNeum, elemGeom
! ==================================================================================================
    implicit none
! ==================================================================================================
    private :: swapMatrToSecant, getMatrTypePred, getMatrTypeCorr, &
               getOptionPred, getOptionCorr, getOptionPost, &
               isMatrUpdatePred, isMatrUpdateCorr, &
               isElasMatr
    public  :: getMatrType, getOption, &
               isMatrUpdate, &
               isRigiMatrCompute, isInteVectCompute, &
               factorSystem, &
               inteForceGetOption, &
               updateLoadBCMatrix, compElemGeom
! ==================================================================================================
    private
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/diinst.h"
#include "asterfort/infdbg.h"
#include "asterfort/isfonc.h"
#include "asterfort/mtdscr.h"
#include "asterfort/ndynlo.h"
#include "asterfort/nmchex.h"
#include "asterfort/nmrinc.h"
#include "asterfort/nmtime.h"
#include "asterfort/NonLinear_type.h"
#include "asterfort/preres.h"
#include "asterfort/romAlgoNLCorrEFMatrixModify.h"
#include "asterfort/utmess.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! getOptionPred
!
! Get name of option for non-linear computation - Prediction
!
! In  matrType         : type of matrix
! In  lImplex          : flag for IMPLEX method
! Out nonLinearOption  : name of option for non-linear computation
!
! --------------------------------------------------------------------------------------------------
    subroutine getOptionPred(matrType, lImplex, nonLinearOption)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=16), intent(in) :: matrType
        aster_logical, intent(in) :: lImplex
        character(len=16), intent(out) :: nonLinearOption
!   ------------------------------------------------------------------------------------------------
!
        nonLinearOption = ' '
!
        if (matrType .eq. 'TANGENTE') then
            nonLinearOption = 'RIGI_MECA_TANG'
        else if (matrType .eq. 'SECANTE') then
            nonLinearOption = 'RIGI_MECA_ELAS'
        else if (matrType .eq. 'ELASTIQUE') then
            nonLinearOption = 'RIGI_MECA'
        else if (matrType .eq. 'DEPL_CALCULE') then
            nonLinearOption = 'RIGI_MECA_TANG'
        else if (matrType .eq. 'EXTRAPOLE') then
            nonLinearOption = 'RIGI_MECA_TANG'
        else
            ASSERT(ASTER_FALSE)
        end if
        if (lImplex) then
            nonLinearOption = 'RIGI_MECA_IMPLEX'
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getOptionCorr
!
! Get name of option for non-linear computation - Correction
!
! In  matrType         : type of matrix
! In  l_update_matr    : flag to update matrix
! In  lImplex          : flag for IMPLEX method
! Out nonLinearOption  : name of option for non-linear computation
!
! --------------------------------------------------------------------------------------------------
    subroutine getOptionCorr(matrType, l_update_matr, lImplex, nonLinearOption)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=16), intent(in) :: matrType
        aster_logical, intent(in) :: l_update_matr, lImplex
        character(len=16), intent(out) :: nonLinearOption
!   ------------------------------------------------------------------------------------------------
!
        nonLinearOption = ' '
!
        if (l_update_matr) then
            if (matrType .eq. 'TANGENTE') then
                nonLinearOption = 'FULL_MECA'
            elseif (matrType .eq. 'ELASTIQUE') then
                nonLinearOption = 'FULL_MECA_ELAS'
            elseif (matrType .eq. 'SECANTE') then
                nonLinearOption = 'FULL_MECA_ELAS'
            else
                ASSERT(ASTER_FALSE)
            end if
        else
            nonLinearOption = 'RAPH_MECA'
            if (lImplex) then
                nonLinearOption = 'RAPH_MECA_IMPLEX'
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getOptionPost
!
! Get name of option for non-linear computation - Buckling
!
! In  matrType           : type of matrix
! In  lImplex            : flag for IMPLEX method
! Out nonLinearOption    : name of option for non-linear computation
!
! --------------------------------------------------------------------------------------------------
    subroutine getOptionPost(matrType, nonLinearOption)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=16), intent(in) :: matrType
        character(len=16), intent(out) :: nonLinearOption
!   ------------------------------------------------------------------------------------------------
!
        nonLinearOption = ' '
!
        if (matrType .eq. 'TANGENTE') then
            nonLinearOption = 'RIGI_MECA_TANG'
        else if (matrType .eq. 'SECANTE') then
            nonLinearOption = 'RIGI_MECA_ELAS'
        else if (matrType .eq. 'ELASTIQUE') then
            nonLinearOption = 'RIGI_MECA'
        else
            ASSERT(ASTER_FALSE)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getOption
!
! Get name of option for non-linear computation
!
! In  phaseType        : name of current phase (prediction/correction/internal forces)
! In  listFuncActi     : list of active functionnalities
! In  matrType         : type of matrix
! Out nonLinearOption  : name of option for non-linear computation
! In  l_update_matr    : flag to update matrix
!
! --------------------------------------------------------------------------------------------------
    subroutine getOption(phaseType, listFuncActi, matrType, nonLinearOption, l_update_matr_)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        integer(kind=8), intent(in) :: phaseType, listFuncActi(*)
        character(len=16), intent(in) :: matrType
        character(len=16), intent(out) :: nonLinearOption
        aster_logical, optional, intent(in) :: l_update_matr_
! - Local
        aster_logical :: lImplex
!   ------------------------------------------------------------------------------------------------
        nonLinearOption = ' '
        lImplex = isfonc(listFuncActi, 'IMPLEX')
        if (phaseType .eq. PRED_EULER) then
            call getOptionPred(matrType, lImplex, nonLinearOption)
        else if (phaseType .eq. CORR_NEWTON) then
            call getOptionCorr(matrType, l_update_matr_, lImplex, nonLinearOption)
        else if (phaseType .eq. POST_BUCKLING) then
            call getOptionPost(matrType, nonLinearOption)
        else
            ASSERT(ASTER_FALSE)
        end if
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getMatrType
!
! Get type of matrix
!
! In  phaseType        : name of current phase (prediction/correction/internal forces)
! In  listFuncActi     : list of active functionnalities
! In  sddisc           : datastructure for discretization
! In  numeTime         : index of current time step
! In  ds_algopara      : datastructure for algorithm parameters
! Out matrType         : type of matrix
! Out reac_iter        : frequency to update matrix (Newton iteration)
! Out reac_incr        : frequency to update matrix (time step)
!
! --------------------------------------------------------------------------------------------------
    subroutine getMatrType(phaseType, listFuncActi, sddisc, numeTime, ds_algopara, &
                           matrType, reac_iter_, reac_incr_)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        integer(kind=8), intent(in) :: phaseType, listFuncActi(*)
        character(len=19), intent(in) :: sddisc
        integer(kind=8), intent(in) :: numeTime
        type(NL_DS_AlgoPara), intent(in) :: ds_algopara
        character(len=16), intent(out) :: matrType
        integer(kind=8), optional, intent(out) :: reac_iter_, reac_incr_
! - Local
        aster_logical :: l_dischoc
!   ------------------------------------------------------------------------------------------------
!
        matrType = ' '
        if (phaseType .eq. PRED_EULER) then
            call getMatrTypePred(sddisc, numeTime, ds_algopara, matrType, reac_incr_)
        else if (phaseType .eq. CORR_NEWTON) then
            l_dischoc = isfonc(listFuncActi, 'DIS_CHOC')
            call getMatrTypeCorr(l_dischoc, sddisc, numeTime, ds_algopara, matrType, reac_iter_)
        else
            ASSERT(ASTER_FALSE)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getMatrTypePred
!
! Get type of matrix for prediction
!
! In  sddisc           : datastructure for discretization
! In  numeTime         : index of current time step
! In  ds_algopara      : datastructure for algorithm parameters
! Out matrType         : type of matrix
! Out reac_incr        : frequency to update matrix (time step)
!
! --------------------------------------------------------------------------------------------------
    subroutine getMatrTypePred(sddisc, numeTime, ds_algopara, matrType, reac_incr)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=19), intent(in) :: sddisc
        integer(kind=8), intent(in) :: numeTime
        type(NL_DS_AlgoPara), intent(in) :: ds_algopara
        character(len=16), intent(out) :: matrType
        integer(kind=8), intent(out) :: reac_incr
! - Local
        aster_logical :: l_swap
!   ------------------------------------------------------------------------------------------------
!
        matrType = ds_algopara%matrix_pred
        reac_incr = ds_algopara%reac_incr

! - Swap matrix to secant ?
        call swapMatrToSecant(sddisc, numeTime, ds_algopara, l_swap)
        if (l_swap) then
            matrType = 'SECANTE'
            reac_incr = 1
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! getMatrTypeCorr
!
! Get type of matrix for correction
!
! In  l_dischoc        : flag if DIS_CHOC elements are present
! In  sddisc           : datastructure for discretization
! In  numeTime         : index of current time step
! In  ds_algopara      : datastructure for algorithm parameters
! Out matrType         : type of matrix
! Out reac_iter        : frequency to update matrix (Newton iteration)
!
! --------------------------------------------------------------------------------------------------
    subroutine getMatrTypeCorr(l_dischoc, sddisc, numeTime, ds_algopara, matrType, reac_iter)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        aster_logical, intent(in) :: l_dischoc
        character(len=19), intent(in) :: sddisc
        integer(kind=8), intent(in) :: numeTime
        type(NL_DS_AlgoPara), intent(in) :: ds_algopara
        character(len=16), intent(out) :: matrType
        integer(kind=8), intent(out) :: reac_iter
! - Local
        aster_logical :: l_swap, l_swapToElastic
!   ------------------------------------------------------------------------------------------------
!
        matrType = ds_algopara%matrix_corr
        reac_iter = ds_algopara%reac_iter
        l_swapToElastic = ds_algopara%l_swapToElastic

! - Swap matrix to secant ?
        call swapMatrToSecant(sddisc, numeTime, ds_algopara, l_swap)
        if (l_swap) then
            matrType = 'SECANTE'
            reac_iter = ds_algopara%reac_iter_elas
        end if

! - Change matrix if DIS_CHOC
        if (l_dischoc) then
            matrType = 'TANGENTE'
        end if

! - Swap matrix to elastic if contact is not stabilized and PRED_CONTACT = OUI
        if (l_swapToElastic) then
            matrType = 'ELASTIQUE'
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! swapMatrToSecant
!
! Swap matrix to secant
!
! In  sddisc           : datastructure for discretization
! In  numeTime         : index of current time step
! In  ds_algopara      : datastructure for algorithm parameters
! Out l_swap           : flag to swap matrix to secant
!
! --------------------------------------------------------------------------------------------------
    subroutine swapMatrToSecant(sddisc, numeTime, ds_algopara, l_swap)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=19), intent(in) :: sddisc
        integer(kind=8), intent(in) :: numeTime
        type(NL_DS_AlgoPara), intent(in) :: ds_algopara
        aster_logical, intent(out) :: l_swap
! - Local
        real(kind=8) :: pas_mini_elas, timeIncr, timePrev, timeCurr
!   ------------------------------------------------------------------------------------------------
!
        l_swap = ASTER_FALSE
!
        pas_mini_elas = ds_algopara%pas_mini_elas
        timePrev = diinst(sddisc, numeTime-1)
        timeCurr = diinst(sddisc, numeTime)
        timeIncr = timeCurr-timePrev

! - Elastic matrix for timeIncr < PAS_MINI_ELAS
        if (pas_mini_elas .ge. 0.d0) then
            if (abs(timeIncr) .lt. pas_mini_elas) then
                l_swap = ASTER_TRUE
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! isMatrUpdate
!
! Update global matrix ?
!
! In  phaseType        : name of current phase of algorithm
! In  matrType         : type of matrix
! In  listFuncActi     : list of active functionnalities
! In  nlDynaDamping    : damping parameters
! In  ds_system        : datastructure for non-linear system management
! Out l_update_matr    : flag to update matrix
! In  numeTime         : index of current time step
! In  iter_newt        : index of current Newton iteration
! In  reac_iter        : frequency to update matrix (Newton iteration)
! In  reac_incr        : frequency to update matrix (time step)
!
! --------------------------------------------------------------------------------------------------
    subroutine isMatrUpdate(phaseType, matrType, listFuncActi, &
                            nlDynaDamping, ds_system, &
                            l_update_matr, &
                            nume_inst_, iter_newt_, &
                            reac_iter_, reac_incr_)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        integer(kind=8), intent(in) :: phaseType
        character(len=16), intent(in) :: matrType
        integer(kind=8), intent(in) :: listFuncActi(*)
        type(NLDYNA_DAMPING), intent(in) :: nlDynaDamping
        type(NL_DS_System), intent(in) :: ds_system
        aster_logical, intent(out) :: l_update_matr
        integer(kind=8), optional, intent(in) :: nume_inst_, iter_newt_
        integer(kind=8), optional, intent(in) :: reac_iter_, reac_incr_
! - Local
        aster_logical :: l_matr_elas
        aster_logical :: l_cont_elem, lDyna
        aster_logical :: lDampMatrix, l_dischoc, l_varc, l_elas_fo
!   ------------------------------------------------------------------------------------------------
!
        l_update_matr = ASTER_FALSE

! - Active functionnalities
        lDyna = isfonc(listFuncActi, 'DYNAMIQUE')
        lDampMatrix = nlDynaDamping%hasMatrDamp
        l_cont_elem = isfonc(listFuncActi, 'ELT_CONTACT')
        l_dischoc = isfonc(listFuncActi, 'DIS_CHOC')
        l_varc = isfonc(listFuncActi, 'EXI_VARC')
        l_elas_fo = isfonc(listFuncActi, 'ELAS_FO')

! - Is matrix is elastic ?
        call isElasMatr(matrType, l_matr_elas)

! - Update matrix ?
        if (phaseType .eq. PRED_EULER) then
            call isMatrUpdatePred(lDyna, lDampMatrix, l_matr_elas, &
                                  l_varc, l_elas_fo, &
                                  reac_incr_, nume_inst_, &
                                  l_update_matr)
        else if (phaseType .eq. CORR_NEWTON) then
            call isMatrUpdateCorr(matrType, iter_newt_, reac_iter_, l_update_matr)
        else
            ASSERT(ASTER_FALSE)
        end if

! - Update matrix if elementary contact (CONTINUE/XFEM/LAC)
        if (l_cont_elem) then
            if (.not. l_update_matr) then
                if (matrType .ne. 'ELASTIQUE') then
                    call utmess('A', 'MECANONLINE5_4')
                end if
                l_update_matr = ASTER_TRUE
            end if
        end if

! - Update matrix if DIS_CHOC elements
        if (l_dischoc) then
            if (.not. l_update_matr) then
                if (matrType .ne. 'ELASTIQUE') then
                    call utmess('A', 'MECANONLINE5_5')
                end if
                l_update_matr = ASTER_TRUE
            end if
        end if

! - Update if contact matrix in global matrix
        if (ds_system%l_matr_cont) then
            if (.not. l_update_matr) then
                if (matrType .ne. 'ELASTIQUE') then
                    call utmess('A', 'MECANONLINE5_5')
                end if
                l_update_matr = ASTER_TRUE
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! isMatrUpdatePred
!
! Update global matrix at prediction ?
!
! In  l_dyna           : flag for dynamic
! In  lDampMatrix      : flag for damping matrix
! In  l_matr_elas      : flag if matrix is elastic
! In  l_varc           : flag for external state variables (AFFE_VARC)
! In  l_elas_fo        : flag if elasticity parameters are functions
! In  reac_incr        : frequency of matrix update for time step
! In  numeTime         : index of current time step
! Out l_update_matr    : flag to update matrix
!
! --------------------------------------------------------------------------------------------------
    subroutine isMatrUpdatePred(l_dyna, lDampMatrix, l_matr_elas, &
                                l_varc, l_elas_fo, &
                                reac_incr, numeTime, &
                                l_update_matr)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        aster_logical, intent(in) :: l_dyna, lDampMatrix, l_matr_elas, l_varc, l_elas_fo
        integer(kind=8), intent(in) :: reac_incr, numeTime
        aster_logical, intent(out) :: l_update_matr
! - Local
        aster_logical :: l_first_step
!   ------------------------------------------------------------------------------------------------
!
        l_update_matr = ASTER_FALSE
        l_first_step = numeTime .le. 1
!
        if ((reac_incr .eq. 0) .and. (numeTime .ne. 1)) then
            l_update_matr = ASTER_FALSE
        end if
        if (numeTime .eq. 1) then
            l_update_matr = ASTER_TRUE
        end if
        if ((reac_incr .ne. 0) .and. (numeTime .ne. 1)) then
            l_update_matr = mod(numeTime-1, reac_incr) .eq. 0
        end if

! - Update matrix if damping matrix
        if (l_dyna .and. lDampMatrix .and. l_first_step) then
            l_update_matr = ASTER_TRUE
        end if

! - Update matrix if command variables and elastic function
        if (l_matr_elas .and. l_varc .and. l_elas_fo) then
            l_update_matr = ASTER_TRUE
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! isMatrUpdateCorr
!
! Update global matrix at correction ?
!
! In  matrType         : type of matrix
! In  iter_newt        : index of current Newton iteration
! In  reac_iter        : frequency of matrix update for Newton
! Out l_update_matr    : flag to update matrix
!
! --------------------------------------------------------------------------------------------------
    subroutine isMatrUpdateCorr(matrType, iter_newt, reac_iter, l_update_matr)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=16), intent(in) :: matrType
        integer(kind=8), intent(in) :: iter_newt, reac_iter
        aster_logical, intent(out) :: l_update_matr
!   ------------------------------------------------------------------------------------------------
        l_update_matr = ASTER_FALSE
!
        if ((matrType .eq. 'TANGENTE') .or. (matrType .eq. 'SECANTE')) then
            l_update_matr = ASTER_FALSE
            if (reac_iter .ne. 0) then
                l_update_matr = mod(iter_newt+1, reac_iter) .eq. 0
            end if
        else
            l_update_matr = ASTER_FALSE
        end if
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! isElasMatr
!
! Is matrix is elastic ?
!
! In  matrType         : type of matrix
! Out l_matr_elas      : flag if matrix is elastic
!
! --------------------------------------------------------------------------------------------------
    subroutine isElasMatr(matrType, l_matr_elas)
! - Parameters
        character(len=16), intent(in) :: matrType
        aster_logical, intent(out) :: l_matr_elas
!   ------------------------------------------------------------------------------------------------
        l_matr_elas = ASTER_FALSE
        if (matrType .eq. 'ELASTIQUE') then
            l_matr_elas = ASTER_TRUE
        else
            l_matr_elas = ASTER_FALSE
        end if
!   -----------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! isRigiMatrCompute
!
! Do the rigidity matrices have to be calculated/assembled ?
!
! In  phaseType        : name of current phase (prediction/correction/internal forces)
! In  sddyna           : name of dynamic parameters datastructure
! In  numeTime         : index of current time step
! In  l_update_matr    : flag to update matrix
! In  l_comp_damp      : flag if damp elementary matrices have to be calculated
! Out l_comp_rigi      : flag if rigidity elementary matrices have to be calculated
! Out l_asse_rigi      : flag if rigidity elementary matrices have to be assembled
!
! --------------------------------------------------------------------------------------------------
    subroutine isRigiMatrCompute(phaseType, &
                                 sddyna, numeTime, &
                                 l_update_matr, lDampCompute, &
                                 lRigiCompute, lRigiAssemble)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        integer(kind=8), intent(in) :: phaseType
        character(len=19), intent(in) :: sddyna
        integer(kind=8), intent(in) :: numeTime
        aster_logical, intent(in) :: l_update_matr, lDampCompute
        aster_logical, intent(out) :: lRigiCompute, lRigiAssemble
! - Local
        aster_logical :: l_first_step, l_shift_mass
!   ------------------------------------------------------------------------------------------------
!
        lRigiCompute = ASTER_FALSE
        lRigiAssemble = ASTER_FALSE

! ----- First step ?
        l_first_step = numeTime .le. 1

! ----- Active functionnalities
        l_shift_mass = ndynlo(sddyna, 'COEF_MASS_SHIFT')

! ----- Rigidity matrices have to be calculated ?
        if (phaseType .eq. PRED_EULER) then
            lRigiCompute = l_update_matr
            lRigiAssemble = ASTER_TRUE
        end if

! ----- Rayleigh: need update of rigidity matrices
        if (lDampCompute) then
            lRigiCompute = ASTER_TRUE
        end if

! ----- For COEF_MASS_SHIFT
        if (l_shift_mass .and. l_first_step) then
            lRigiCompute = ASTER_TRUE
            lRigiAssemble = ASTER_TRUE
        end if

! ----- From status of global matrix
        if (l_update_matr) then
            lRigiAssemble = ASTER_TRUE
        end if
!
!   -----------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! isInteVectCompute
!
! Do the internal forces vectors have to be calculated ?
!
! In  phaseType        : name of current phase of algorithm
! In  listFuncActi     : list of active functionnalities
! In  nonLinearOption  : name of option for non-linear computation
! In  iter_newt        : index of current Newton iteration
! In  l_comp_rigi      : flag if rigidity elementary matrices have to be calculated
! Out l_comp_fint      : flag if internal forces elementary vectors have to be calculated
!
! --------------------------------------------------------------------------------------------------
    subroutine isInteVectCompute(phaseType, listFuncActi, &
                                 nonLinearOption, iter_newt, &
                                 l_comp_rigi, l_comp_fint)
! - Parameters
        integer(kind=8), intent(in) :: phaseType, listFuncActi(*)
        character(len=16), intent(in) :: nonLinearOption
        integer(kind=8), intent(in) :: iter_newt
        aster_logical, intent(in) :: l_comp_rigi
        aster_logical, intent(out) :: l_comp_fint
!   ------------------------------------------------------------------------------------------------
! - Local
        aster_logical :: l_unil, l_cont_disc, l_line_search
!   ------------------------------------------------------------------------------------------------
        l_comp_fint = ASTER_FALSE
!
! - Active functionnalities
!
        l_unil = isfonc(listFuncActi, 'LIAISON_UNILATER')
        l_cont_disc = isfonc(listFuncActi, 'CONT_DISCRET')
        l_line_search = isfonc(listFuncActi, 'RECH_LINE')
!
! - Is internal forces elementary vectors have to be calculated ?
!
        if (phaseType .eq. PRED_EULER) then
            if (nonLinearOption(1:9) .eq. 'FULL_MECA') then
                l_comp_fint = ASTER_TRUE
            else if (nonLinearOption(1:10) .eq. 'RIGI_MECA ') then
                l_comp_fint = ASTER_FALSE
            else if (nonLinearOption(1:10) .eq. 'RIGI_MECA_') then
                l_comp_fint = ASTER_FALSE
            else if (nonLinearOption(1:9) .eq. 'RAPH_MECA') then
                l_comp_fint = ASTER_TRUE
            else
                ASSERT(ASTER_FALSE)
            end if
        else if (phaseType .eq. CORR_NEWTON) then
            l_comp_fint = l_comp_rigi
        else if (phaseType .eq. INTE_FORCE) then
            if (.not. l_line_search .or. iter_newt .eq. 0) then
                l_comp_fint = ASTER_TRUE
            else
                if (nonLinearOption .eq. 'FULL_MECA') then
                    l_comp_fint = ASTER_TRUE
                else
                    l_comp_fint = ASTER_FALSE
                end if
            end if
            if (l_cont_disc .or. l_unil) then
                l_comp_fint = ASTER_TRUE
            end if
        else
            ASSERT(ASTER_FALSE)
        end if
!   -----------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! factorSystem
!
! Factorization of linear system
!
! In  listFuncActi     : list of active functionnalities
! IO  ds_measure       : datastructure for measure and statistics management
! In  ds_algorom       : datastructure for ROM parameters
! In  numeDof          : name of numeDof object (numbering equation)
! In  solveu           : name of datastructure for solver
! In  ds_system        : datastructure for non-linear system management
! In  maprec           : matrix for pre-conditionning
! In  matass           : matrix of linear system
! Out faccvg           : error from factorization
!
! --------------------------------------------------------------------------------------------------
    subroutine factorSystem(listFuncActi, ds_measure, ds_algorom, &
                            numeDof, solveu, maprec, matass, &
                            faccvg)
! - Parameters
        integer(kind=8), intent(in) :: listFuncActi(*)
        type(NL_DS_Measure), intent(inout) :: ds_measure
        type(ROM_DS_AlgoPara), intent(in) :: ds_algorom
        character(len=19), intent(in) :: maprec, matass, solveu
        character(len=24), intent(in) :: numeDof
        integer(kind=8), intent(out) :: faccvg
!   ------------------------------------------------------------------------------------------------
! - Local
        aster_logical :: l_rom
        integer(kind=8) :: npvneg, ifm, niv
!   ------------------------------------------------------------------------------------------------
!
        call infdbg('MECANONLINE', ifm, niv)
        l_rom = isfonc(listFuncActi, 'ROM')
        faccvg = 0
!
        call nmtime(ds_measure, 'Init', 'Factor')
        call nmtime(ds_measure, 'Launch', 'Factor')
        if (l_rom .and. ds_algorom%phase .eq. 'HROM') then
            call mtdscr(matass)
        elseif (l_rom .and. ds_algorom%phase .eq. 'CORR_EF') then
            call mtdscr(matass)
            call romAlgoNLCorrEFMatrixModify(numeDof, matass, ds_algorom)
            call preres(solveu, 'V', faccvg, maprec, matass, npvneg, -9999)
            if (niv .ge. 2) then
                call utmess('I', 'MECANONLINE13_42')
            end if
        else
            call preres(solveu, 'V', faccvg, maprec, matass, npvneg, -9999)
            if (niv .ge. 2) then
                call utmess('I', 'MECANONLINE13_42')
            end if
        end if
        call nmtime(ds_measure, 'Stop', 'Factor')
!   only increase nb of factorization if the matrix has indeed been factored
!   this can not be the case with the use of LDLT_SP/DP
        if (really_factored) call nmrinc(ds_measure, 'Factor')
!
!   -----------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! inteForceGetOption
!
! Get options to compute internal forces
!
! In  phaseType        : name of current phase of algorithm
! In  listFuncActi     : list of active functionnalities
! In  ds_algorom       : datastructure for ROM parameters
! Out lNodeComp        : compute internal forces without integration (FORC_NODA)
! Out lInteComp        : compute internal forces with integration (RAPH_NODA, RIGI_MECA_TANG)
! Out typeAsse         : type of assembly for internal forces
!
! --------------------------------------------------------------------------------------------------
    subroutine inteForceGetOption(phaseType, listFuncActi, ds_algorom, &
                                  lNodeComp, lInteComp, typeAsse)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        integer(kind=8), intent(in) :: phaseType, listFuncActi(*)
        aster_logical, intent(out) :: lNodeComp, lInteComp
        integer(kind=8), intent(out) :: typeAsse
        type(ROM_DS_AlgoPara), intent(in) :: ds_algorom
! - Local
        aster_logical :: lDyna, lRomCorr, lHHO
!   ------------------------------------------------------------------------------------------------
!
        lNodeComp = ASTER_FALSE
        lInteComp = ASTER_FALSE
        typeAsse = INTE_FORCE_NONE

! - Active functionnalites
        lDyna = isfonc(listFuncActi, 'DYNAMIQUE')
        lHHO = isfonc(listFuncActi, 'HHO')
        lRomCorr = ASTER_FALSE
        if (ds_algorom%l_rom) then
            lRomCorr = ds_algorom%phase .eq. 'CORR_EF'
        end if

! - Options to compute internal forces
        if (phaseType .eq. PRED_EULER) then
            lNodeComp = ASTER_TRUE
            lInteComp = (lDyna .or. lRomCorr .or. lHHO)
        elseif (phaseType .eq. CORR_NEWTON) then
            lNodeComp = ASTER_FALSE
            lInteComp = ASTER_TRUE
        else
            ASSERT(ASTER_FALSE)
        end if

! - Which second member ?
        if (phaseType .eq. PRED_EULER) then
            typeAsse = INTE_FORCE_COMB
            if (lDyna .or. lRomCorr .or. lHHO) then
                typeAsse = INTE_FORCE_INTE
            end if
        elseif (phaseType .eq. CORR_NEWTON) then
            typeAsse = INTE_FORCE_INTE
        else
            ASSERT(ASTER_FALSE)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! updateLoadBCMatrix
!
! Update elementary matrices for loads and boundary conditions (undead cases)
!
! In  listFuncActi     : list of active functionnalities
! In  listLoad         : name of datastructure for list of loads
! In  sddisc           : datastructure for discretization
! In  numeTime         : index of current time step
! In  model            : model
! In  caraElem         : field for elementary characteristics
! In  ds_material      : datastructure for material parameters
! In  ds_constitutive  : datastructure for constitutive laws management
! In  hval_incr        : hat-variable for incremental values fields
! In  hval_algo        : hat-variable for algorithms fields
! In  hval_meelem      : hat-variable for elementary matrix
!
! --------------------------------------------------------------------------------------------------
    subroutine updateLoadBCMatrix(listFuncActi, listLoad, &
                                  sddisc, numeTime, &
                                  modelZ, caraElemZ, &
                                  ds_material, ds_constitutive, &
                                  hval_incr, hval_algo, &
                                  hval_meelem)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        integer(kind=8), intent(in) :: listFuncActi(*)
        character(len=19), intent(in) :: listLoad, sddisc
        integer(kind=8), intent(in) :: numeTime
        character(len=24), intent(in) :: modelZ, caraElemZ
        type(NL_DS_Material), intent(in) :: ds_material
        type(NL_DS_Constitutive), intent(in) :: ds_constitutive
        character(len=19), intent(in) :: hval_incr(*), hval_algo(*)
        character(len=19), intent(in) :: hval_meelem(*)
! - Local
        character(len=24) :: model, caraElem
        aster_logical :: lDiriUndead, lNeumUndead
        character(len=24) :: diriElem, neumElem
        character(len=19) :: dispPrev, dispCumu
        real(kind=8) :: timePrev, timeCurr
!   ------------------------------------------------------------------------------------------------
!
        model = modelZ
        caraElem = caraElemZ

! - Active functionnalites
        lNeumUndead = isfonc(listFuncActi, 'NEUM_UNDEAD')
        lDiriUndead = isfonc(listFuncActi, 'DIRI_UNDEAD')

! - Upddate elementary matrices for BC
        if (lDiriUndead) then
            call nmchex(hval_meelem, 'MEELEM', 'MEDIRI', diriElem)
            call elemDiri(model, listLoad, diriElem)
        end if

! - Update elementary matrices for Neumann loads (undead)
        if (lNeumUndead) then
            call nmchex(hval_incr, 'VALINC', 'DEPMOI', dispPrev)
            call nmchex(hval_algo, 'SOLALG', 'DEPDEL', dispCumu)
            call nmchex(hval_meelem, 'MEELEM', 'MESUIV', neumElem)
            timePrev = diinst(sddisc, numeTime-1)
            timeCurr = diinst(sddisc, numeTime)
            call elemNeum(listLoad, model, caraElem, &
                          ds_material%mater, ds_material%mateco, &
                          ds_constitutive%compor, &
                          timePrev, timeCurr, &
                          dispPrev, dispCumu, &
                          neumElem)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compElemGeom
!
! Compute elementary matrices for geometry (buckling)
!
! In  model            : model
! In  caraElem         : field for elementary characteristics
! In  ds_material      : datastructure for material parameters
! In  hval_incr        : hat-variable for incremental values fields
! In  hval_meelem      : hat-variable for elementary matrix
!
! --------------------------------------------------------------------------------------------------
    subroutine compElemGeom(modelZ, caraElemZ, &
                            ds_material, &
                            hval_incr, hval_meelem)
!   ------------------------------------------------------------------------------------------------
! - Parameters
        character(len=24), intent(in) :: modelZ, caraElemZ
        type(NL_DS_Material), intent(in) :: ds_material
        character(len=19), intent(in) :: hval_incr(*), hval_meelem(*)
! - Local
        character(len=24) :: model, caraElem
        character(len=24) :: geomElem
        character(len=19) :: sigmCurr, strxCurr
!   ------------------------------------------------------------------------------------------------
!
        model = modelZ
        caraElem = caraElemZ

        call nmchex(hval_incr, 'VALINC', 'SIGPLU', sigmCurr)
        call nmchex(hval_incr, 'VALINC', 'STRPLU', strxCurr)
        call nmchex(hval_meelem, 'MEELEM', 'MEGEOM', geomElem)
        call elemGeom(model, caraElem, &
                      ds_material%mateco, &
                      strxCurr, sigmCurr, &
                      geomElem)
!
!   ----------------------------------------------------------------------------------------------
    end subroutine
!
end module NonLinear_module
