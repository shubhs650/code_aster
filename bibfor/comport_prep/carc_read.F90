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
subroutine carc_read(prepMapCarcri, model)
!
    use BehaviourPrepare_type
    implicit none
!
#include "asterc/lccree.h"
#include "asterc/lcdiscard.h"
#include "asterc/lcsymm.h"
#include "asterc/lctest.h"
#include "asterc/r8vide.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/comp_meca_l.h"
#include "asterfort/comp_meca_rkit.h"
#include "asterfort/comp_read_mesh.h"
#include "asterfort/compGetMecaPart.h"
#include "asterfort/compGetRelation.h"
#include "asterfort/dismoi.h"
#include "asterfort/exicp.h"
#include "asterfort/getBehaviourAlgo.h"
#include "asterfort/getBehaviourPara.h"
#include "asterfort/getExternalBehaviourPara.h"
#include "asterfort/getExternalStateVariable.h"
#include "asterfort/getTHMPara.h"
#include "asterfort/getvis.h"
#include "asterfort/getvr8.h"
#include "asterfort/getvtx.h"
#include "asterfort/jeveuo.h"
#include "asterfort/utmess.h"
!
    type(BehaviourPrep_MapCarcri), intent(inout) :: prepMapCarcri
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
! IO  prepMapCarcri    : datastructure to construct CARCRI map
! In  model            : name of model
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: factorKeyword = 'COMPORTEMENT'
    character(len=24), parameter :: cellAffeJv = '&&CARCREAD.LIST'
    aster_logical :: lAllCellAffe
    integer(kind=8) :: nbCellAffe
    aster_logical :: lPlaneStress
    character(len=16) :: answer
    integer(kind=8) :: iFactorKeyword, iret, nbFactorKeyword
    character(len=16) :: type_matr_tang, method, post_iter
    real(kind=8) :: parm_theta, vale_pert_rela
    real(kind=8) :: resi_deborst_max
    real(kind=8) :: resi_radi_rela
    real(kind=8) :: parm_theta_thm, parm_alpha_thm
    integer(kind=8) :: type_matr_t, iter_inte_pas, iter_deborst_max
    integer(kind=8) :: ipostiter, iveriborne
    character(len=8) :: mesh
    character(len=16) :: relaCompPY, defoCompPY, relaMecaPY
    character(len=16) :: veri_borne
    character(len=16) :: kitComp(4)
    character(len=16) :: defoComp, relaComp
    character(len=16) :: thmc_comp, hydr_comp, ther_comp, relaMeca
    character(len=19) :: modelFED
    aster_logical :: l_kit_thm, l_kit_ddi, l_exist_thm
    aster_logical :: l_kit
    aster_logical :: l_mfront_proto, l_mfront_offi
    aster_logical ::  l_matr_unsymm
    integer(kind=8) :: variExteCode(2)
    integer(kind=8), pointer :: modelCell(:) => null()
    character(len=16) :: algo_inte
    real(kind=8) :: algo_inte_r
    real(kind=8), pointer :: resi_inte_p => null()
    integer(kind=8), pointer :: iter_inte_maxi_p => null()
    type(BehaviourPrep_Exte) :: prepExte
!
! --------------------------------------------------------------------------------------------------
!
    iFactorKeyword = 0
    iret = 0
    type_matr_tang = ' '
    method = ' '
    post_iter = ' '
    parm_theta = 0.d0
    vale_pert_rela = 0.d0
    resi_deborst_max = 0.d0
    resi_radi_rela = 0.d0
    parm_theta_thm = 0.d0
    parm_alpha_thm = 0.d0
    type_matr_t = 0
    iter_inte_pas = 0
    iter_deborst_max = 0
    ipostiter = 0
    iveriborne = 0
    veri_borne = ' '
    kitComp(1:4) = (/'VIDE', 'VIDE', 'VIDE', 'VIDE'/)
    defoComp = ' '
    relaComp = ' '
    thmc_comp = ' '
    hydr_comp = ' '
    ther_comp = ' '
    relaMeca = ' '
    l_kit_thm = ASTER_FALSE
    l_kit_ddi = ASTER_FALSE
    l_exist_thm = ASTER_FALSE
    l_kit = ASTER_FALSE
    nbFactorKeyword = prepMapCarcri%nb_comp
    l_exist_thm = ASTER_FALSE

! - Pointer to the list of cells in model
    mesh = ' '
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=modelFED)
    call jeveuo(modelFED//'.TYFE', 'L', vi=modelCell)

! - Read informations
    do iFactorKeyword = 1, nbFactorKeyword
! ----- Get RELATION from command file
        call compGetRelation(factorKeyword, iFactorKeyword, relaComp)

! ----- Get DEFORMATION from command file
        call getvtx(factorKeyword, 'DEFORMATION', iocc=iFactorKeyword, scal=defoComp)

! ----- Detection of specific cases
        call comp_meca_l(relaComp, 'KIT', l_kit)
        call comp_meca_l(relaComp, 'KIT_THM', l_kit_thm)
        call comp_meca_l(relaComp, 'KIT_DDI', l_kit_ddi)
        call comp_meca_l(relaComp, 'MFRONT_OFFI', l_mfront_offi)
        call comp_meca_l(relaComp, 'MFRONT_PROTO', l_mfront_proto)
        if (l_kit_thm) then
            l_exist_thm = ASTER_TRUE
        end if

! ----- For KIT
        if (l_kit) then
            call comp_meca_rkit(factorKeyword, iFactorKeyword, relaComp, kitComp)
        end if

! ----- Get mechanical part of behaviour
        call compGetMecaPart(relaComp, kitComp, relaMeca)

! ----- Coding comportment (Python)
        call lccree(1, relaComp, relaCompPY)
        call lccree(1, defoComp, defoCompPY)
        call lccree(1, relaMeca, relaMecaPY)

! ----- Symmetric or not ?
        l_matr_unsymm = ASTER_FALSE
        call lcsymm(relaCompPY, answer)
        l_matr_unsymm = l_matr_unsymm .or. answer .eq. 'No'
        call lcsymm(relaMecaPY, answer)
        l_matr_unsymm = l_matr_unsymm .or. answer .eq. 'No'
        call lcsymm(defoCompPY, answer)
        l_matr_unsymm = l_matr_unsymm .or. answer .eq. 'No'
        if (l_mfront_proto) then
            call getvtx(factorKeyword, 'SYME_MATR_TANG', iocc=iFactorKeyword, &
                        scal=answer, nbret=iret)
            if (iret .ne. 0) then
                l_matr_unsymm = l_matr_unsymm .or. answer .eq. 'NON'
            end if
        end if

! ----- Get ITER_INTE_PAS
        call getvis(factorKeyword, 'ITER_INTE_PAS', iocc=iFactorKeyword, &
                    scal=iter_inte_pas, nbret=iret)
        if (iret .eq. 0) then
            iter_inte_pas = 0
        end if

! ----- Get ITER_CPLAN_MAXI/RESI_CPLAN_MAXI/RESI_CPLAN_RELA (Deborst method)
        iter_deborst_max = 1
        call getvis(factorKeyword, 'ITER_CPLAN_MAXI', iocc=iFactorKeyword, &
                    scal=iter_deborst_max)

        call getvr8(factorKeyword, 'RESI_CPLAN_MAXI', iocc=iFactorKeyword, &
                    scal=resi_deborst_max, nbret=iret)
        if (iret .eq. 0) then
            call getvr8(factorKeyword, 'RESI_CPLAN_RELA', iocc=iFactorKeyword, &
                        scal=resi_deborst_max, nbret=iret)
            if (iret .eq. 0) then
                resi_deborst_max = r8vide()
            else
                resi_deborst_max = -resi_deborst_max
            end if
        end if

! ----- Get TYPE_MATR_TANG/VALE_PERT_RELA
        vale_pert_rela = 0.d0
        type_matr_t = 0
        type_matr_tang = ' '
        call getvtx(factorKeyword, 'TYPE_MATR_TANG', iocc=iFactorKeyword, &
                    scal=type_matr_tang, nbret=iret)
        if (iret .eq. 0) then
            type_matr_t = 0
        else
            if (type_matr_tang .eq. 'PERTURBATION') then
                type_matr_t = 1
                call getvr8(factorKeyword, 'VALE_PERT_RELA', iocc=iFactorKeyword, &
                            scal=vale_pert_rela)
            else if (type_matr_tang .eq. 'VERIFICATION') then
                type_matr_t = 2
                call getvr8(factorKeyword, 'VALE_PERT_RELA', iocc=iFactorKeyword, &
                            scal=vale_pert_rela)
            else if (type_matr_tang .eq. 'MATR_ELAS') then
                type_matr_t = 3
            else if (type_matr_tang .eq. 'MATR_ENDO') then
                type_matr_t = 4
            else
                ASSERT(ASTER_FALSE)
            end if
            call lctest(relaCompPY, 'TYPE_MATR_TANG', type_matr_tang, iret)
            if (iret .eq. 0) then
                call utmess('F', 'COMPOR1_46', nk=2, valk=[type_matr_tang, relaComp])
            end if
        end if

! ----- Get PARM_THETA (for viscous laws)
        parm_theta = 1.d0
        call getvr8(factorKeyword, 'PARM_THETA', iocc=iFactorKeyword, scal=parm_theta)

! ----- Get RESI_RADI_RELA
        if (type_matr_t .eq. 0) then
            call getvr8(factorKeyword, 'RESI_RADI_RELA', iocc=iFactorKeyword, &
                        scal=resi_radi_rela, nbret=iret)
            if (iret .eq. 0) then
                resi_radi_rela = -10.d0
            end if
        end if

! ----- Get POST_ITER
        ipostiter = 0
        post_iter = ' '
        if (type_matr_t .eq. 0) then
            call getvtx(factorKeyword, 'POST_ITER', iocc=iFactorKeyword, &
                        scal=post_iter, nbret=iret)
            if (iret .eq. 1) then
                if (post_iter .eq. 'CRIT_RUPT') then
                    ipostiter = 1
                end if
            end if
        end if

! ----- Get VERI_BORNE
        iveriborne = 0
        call getvtx(factorKeyword, 'VERI_BORNE', iocc=iFactorKeyword, &
                    scal=veri_borne, nbret=iret)
        if (iret .eq. 0) then
            iveriborne = 2
        else
            if (veri_borne .eq. 'ARRET') then
                iveriborne = 2
            elseif (veri_borne .eq. 'MESSAGE') then
                iveriborne = 1
            else
                iveriborne = 0
            end if
        end if

! ----- Get affectation
        call comp_read_mesh(mesh, factorKeyword, iFactorKeyword, &
                            cellAffeJv, lAllCellAffe, nbCellAffe)

! ----- Get parameters for external programs (MFRONT/UMAT)
        call getExternalBehaviourPara(mesh, modelCell, &
                                      cellAffeJv, lAllCellAffe, nbCellAffe, &
                                      relaComp, relaCompPY, relaMeca, defoComp, &
                                      factorKeyword, iFactorKeyword, &
                                      prepExte)

! ----- Detect plane stress elements
        lPlaneStress = exicp(model, lAllCellAffe, cellAffeJv, nbCellAffe)

! ----- Get ALGO_INTE
        algo_inte_r = 0.d0
        call getBehaviourAlgo(lPlaneStress, relaComp, &
                              relaCompPY, relaMecaPY, &
                              factorKeyword, iFactorKeyword, &
                              algo_inte, algo_inte_r)

! ----- Get RESI_INTE/ITER_INTE_MAXI
        call getBehaviourPara(l_mfront_proto, l_kit_thm, &
                              factorKeyword, iFactorKeyword, algo_inte, &
                              iter_inte_maxi_p, resi_inte_p)

! ----- Get external state variables
        call getExternalStateVariable(relaComp, relaCompPY, &
                                      l_mfront_offi, l_mfront_proto, &
                                      prepExte%adrsMGIS, variExteCode)

! ----- Discard
        call lcdiscard(relaMecaPY)
        call lcdiscard(relaCompPY)
        call lcdiscard(defoCompPY)

! ----- Save parameters
        prepMapCarcri%prepCrit(iFactorKeyword)%rela_comp = relaComp
        prepMapCarcri%prepCrit(iFactorKeyword)%meca_comp = relaMeca
        prepMapCarcri%prepCrit(iFactorKeyword)%type_matr_t = type_matr_t
        prepMapCarcri%prepCrit(iFactorKeyword)%parm_theta = parm_theta
        prepMapCarcri%prepCrit(iFactorKeyword)%iter_inte_pas = iter_inte_pas
        prepMapCarcri%prepCrit(iFactorKeyword)%vale_pert_rela = vale_pert_rela
        prepMapCarcri%prepCrit(iFactorKeyword)%resi_deborst_max = resi_deborst_max
        prepMapCarcri%prepCrit(iFactorKeyword)%iter_deborst_max = iter_deborst_max
        prepMapCarcri%prepCrit(iFactorKeyword)%resi_radi_rela = resi_radi_rela
        prepMapCarcri%prepCrit(iFactorKeyword)%ipostiter = ipostiter
        prepMapCarcri%prepCrit(iFactorKeyword)%iveriborne = iveriborne
        prepMapCarcri%prepCrit(iFactorKeyword)%l_matr_unsymm = l_matr_unsymm
        prepMapCarcri%prepCrit(iFactorKeyword)%algo_inte_r = algo_inte_r
        if (associated(resi_inte_p)) then
            allocate (prepMapCarcri%prepCrit(iFactorKeyword)%resi_inte)
            prepMapCarcri%prepCrit(iFactorKeyword)%resi_inte = resi_inte_p
            deallocate (resi_inte_p)
        end if
        if (associated(iter_inte_maxi_p)) then
            allocate (prepMapCarcri%prepCrit(iFactorKeyword)%iter_inte_maxi)
            prepMapCarcri%prepCrit(iFactorKeyword)%iter_inte_maxi = iter_inte_maxi_p
            deallocate (iter_inte_maxi_p)
        end if
        prepMapCarcri%prepCrit(iFactorKeyword)%jvariext1 = variExteCode(1)
        prepMapCarcri%prepCrit(iFactorKeyword)%jvariext2 = variExteCode(2)
        prepMapCarcri%prepCrit(iFactorKeyword)%prepExte = prepExte
    end do

! - Get SCHEMA_THM parameters
    if (l_exist_thm) then
        call getTHMPara(prepMapCarcri)
    end if
!
end subroutine
