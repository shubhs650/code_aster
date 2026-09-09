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
subroutine nxacmv(model, materField, materCode, caraElem, listLoad, nume_dof, &
                  solver, l_stat, timeMap, timeParaIn, tempIter, &
                  vhydr, varcPrev, varcCurr, cn2mbr_stat, &
                  cn2mbr_tran, matass, maprec, cndiri, cncine, &
                  mediri, comporTher, ds_algorom_)
!
    use Rom_Datastructure_type
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/asasve.h"
#include "asterfort/ascavc.h"
#include "asterfort/ascova.h"
#include "asterfort/asmatr.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mecact.h"
#include "asterfort/merxth.h"
#include "asterfort/mtdscr.h"
#include "asterfort/preres.h"
#include "asterfort/vechnl.h"
#include "asterfort/vechth.h"
#include "asterfort/vedith.h"
#include "asterfort/vetnth_nonl.h"
#include "asterfort/vrcins.h"
#include "asterfort/vtaxpy.h"
#include "asterfort/vtzero.h"
!
    character(len=8), intent(in) :: model, materField, caraElem
    character(len=24), intent(in) :: materCode, listLoad
    character(len=24), intent(in) :: nume_dof
    character(len=19), intent(in) :: solver
    character(len=24), intent(in) :: timeMap
    character(len=19), intent(in) :: varcPrev
    character(len=19), intent(in) :: varcCurr
    aster_logical, intent(in) :: l_stat
    real(kind=8), intent(in) :: timeParaIn(6)
    character(len=24), intent(in) :: tempIter
    character(len=24), intent(in) :: vhydr
    character(len=24), intent(in) :: cn2mbr_stat
    character(len=24), intent(in) :: cn2mbr_tran
    character(len=24), intent(in) :: matass
    character(len=19), intent(in) :: maprec
    character(len=24), intent(in) :: cndiri
    character(len=24), intent(out) :: cncine
    character(len=24), intent(in) :: mediri
    character(len=24), intent(in) :: comporTher
    type(ROM_DS_AlgoPara), optional, intent(in) :: ds_algorom_
!
! --------------------------------------------------------------------------------------------------
!
! THER_NON_LINE - Algorithm
!
! Compute second members and tangent matrix
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : name of model
! In  materField       : name of material characteristics (field)
! In  caraElem         : name of elementary characteristics (field)
! In  listLoad         : name of datastructure for list of loads
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ibid, ierr, iret
    integer(kind=8) :: jtn, i_vect
    character(len=2) :: codret
    real(kind=8) :: timeCurr, timePrev
    character(len=8), parameter :: nomcmp(6) = (/'INST    ', 'DELTAT  ', &
                                                 'THETA   ', 'KHI     ', &
                                                 'R       ', 'RHO     '/)
    character(len=24) :: ligrmo
    character(len=24) :: vadiri, vachtp, vatntp, vatnti, vachtn
    character(len=24), parameter  :: merigi = '&&METRIG           .RELR'
    character(len=24), parameter  :: memass = '&&METMAS           .RELR'
    character(len=24) :: vediri
    character(len=24) :: vechtp
    character(len=24), parameter  :: vetntp = '&&VETNTP           .RELR'
    character(len=24), parameter  :: vetnti = '&&VETNTI           .RELR'
    character(len=24), parameter  :: vechtn = '&&VECHTN           .RELR'
    character(len=24) :: cntntp
    character(len=24) :: cnchtp
    character(len=24) :: cnchnl
    character(len=24) :: cntnti
    real(kind=8) :: timePara(6)
    character(len=24) :: loadNameJv, loadInfoJv, loadFuncJv
    character(len=24), pointer :: v_resu_elem(:) => null()
    integer(kind=8), parameter :: nb_max = 9
    integer(kind=8) :: nb_vect, nb_matr
    real(kind=8) :: vect_coef(nb_max)
    character(len=19) :: vect_name(nb_max), matr_name(nb_max)
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
!
! - Initializations
!
    vect_coef(:) = 0.d0
    vect_name(:) = ' '
    matr_name(:) = ' '
    cntntp = ' '
    cnchtp = ' '
    cnchnl = ' '
    cntnti = ' '
    call dismoi('NOM_LIGREL', model, 'MODELE', repk=ligrmo)
    vediri = '&&VEDIRI           .RELR'
    vechtp = '&&VECHTP           .RELR'
    vadiri = '&&NTACMV.VADIRI'
    vachtp = '&&NTACMV.VACHTP'
    vatntp = '&&NTACMV.VATNTP'
    vatnti = '&&NTACMV.VATNTI'
    vachtn = '&&NTACMV.VACHTN'
    timeCurr = timeParaIn(1)
    timePrev = timeCurr-timeParaIn(2)

! - Access to datastructure of list of loads
    loadNameJv = listLoad(1:19)//'.LCHA'
    loadInfoJv = listLoad(1:19)//'.INFC'
    loadFuncJv = listLoad(1:19)//'.FCHA'

! - Construct command variables fields
    call vrcins(model, materField, caraElem, timeCurr, varcCurr, codret)
    call vrcins(model, materField, caraElem, timePrev, varcPrev, codret)

! - Create <CARTE> for time
    timePara = timeParaIn
    if (l_stat) then
        timePara(3) = -1
    end if
    call mecact('V', timeMap, 'MODELE', ligrmo, 'INST_R', &
                ncmp=6, lnomcmp=nomcmp, vr=timePara)
!
! - Compute Dirichlet loads (AFFE_CHAR_THER)
!
    call vedith(model, loadNameJv, loadInfoJv, timeMap, vediri)
    call asasve(vediri, nume_dof, 'R', 'D', vadiri)
    call ascova('D', vadiri, loadFuncJv, 'INST', timeCurr, &
                'R', cndiri)
!
! - Compute Dirichlet loads (AFFE_CHAR_CINE)
!
    cncine = ' '
    call ascavc(loadNameJv, loadInfoJv, loadFuncJv, nume_dof, timeCurr, &
                cncine)
!
! - Compute CHAR_THER_EVOLNI
!
    if (.not. l_stat) then
        call vetnth_nonl(model, caraElem, materCode, timeMap, comporTher, &
                         tempIter, vhydr, &
                         varcPrev, varcCurr, &
                         'V', vetntp, vetnti)
        call asasve(vetnti, nume_dof, 'R', 'D', vatnti)
        call jeveuo(vatnti, 'L', jtn)
        cntnti = zk24(jtn)
        call asasve(vetntp, nume_dof, 'R', 'D', vatntp)
        call jeveuo(vatntp, 'L', jtn)
        cntntp = zk24(jtn)
    end if

! - Compute Neumann loads (second member) - Linear part
    call vechth('STAT', &
                model, materCode, &
                loadNameJv, loadInfoJv, &
                timeCurr, &
                vechtp, &
                varcCurrZ_=varcCurr, timeMapZ_=timeMap, tempPrevZ_=tempIter)

    call asasve(vechtp, nume_dof, 'R', 'D', vachtp)
    call ascova('D', vachtp, loadFuncJv, 'INST', timeCurr, &
                'R', cnchtp)
    if (l_stat) then
        call jedetr(vechtp)
    end if
!
! - Compute Neumann loads (second member) - Nonlinear part
!
    call vechnl(model, loadNameJv, loadInfoJv, timeMap, &
                tempIter, vechtn, 'V')
    call asasve(vechtn, nume_dof, 'R', 'D', vachtn)
    call ascova('D', vachtn, ' ', 'INST', timeCurr, &
                'R', cnchnl)
    if (l_stat) then
        call jedetr(vechtn)
    end if
!
! - Compute second members
!
    call vtzero(cn2mbr_stat)
    call vtzero(cn2mbr_tran)
    if (l_stat) then
        nb_vect = 2
        vect_coef(1) = 1.d0
        vect_coef(2) = 1.d0
        vect_name(1) = cnchtp(1:19)
        vect_name(2) = cnchnl(1:19)
        do i_vect = 1, nb_vect
            call vtaxpy(vect_coef(i_vect), vect_name(i_vect), cn2mbr_stat)
        end do
    else
        nb_vect = 3
        vect_coef(1) = 1.d0
        vect_coef(2) = 1.d0
        vect_coef(3) = 1.d0
        vect_name(1) = cnchtp(1:19)
        vect_name(2) = cnchnl(1:19)
        vect_name(3) = cntntp(1:19)
        do i_vect = 1, nb_vect
            call vtaxpy(vect_coef(i_vect), vect_name(i_vect), cn2mbr_stat)
        end do
        vect_name(1) = cnchtp(1:19)
        vect_name(2) = cnchnl(1:19)
        vect_name(3) = cntnti(1:19)
        do i_vect = 1, nb_vect
            call vtaxpy(vect_coef(i_vect), vect_name(i_vect), cn2mbr_tran)
        end do
    end if

! - New <CARTE> for time
    timePara = timeParaIn
    call detrsd("CARTE", timeMap)
    call mecact('V', timeMap, 'MODELE', ligrmo, 'INST_R', &
                ncmp=6, lnomcmp=nomcmp, vr=timePara)

! - Tangent matrix (non-linear) - Material and loads
    call merxth(l_stat, &
                model, caraElem, materCode, &
                loadNameJv, loadInfoJv, &
                timePara, timeMap, &
                tempIter, comporTher, varcCurr, &
                merigi, 'V')
    nb_matr = 0
    call jeexin(merigi(1:8)//'           .RELR', iret)
    if (iret .gt. 0) then
        call jeveuo(merigi(1:8)//'           .RELR', 'L', vk24=v_resu_elem)
        if (v_resu_elem(1) .ne. ' ') then
            nb_matr = nb_matr+1
            matr_name(nb_matr) = merigi(1:19)
        end if
    end if
    call jeexin(mediri(1:8)//'           .RELR', iret)
    if (iret .gt. 0) then
        call jeveuo(mediri(1:8)//'           .RELR', 'L', vk24=v_resu_elem)
        if (v_resu_elem(1) .ne. ' ') then
            nb_matr = nb_matr+1
            matr_name(nb_matr) = mediri(1:19)
        end if
    end if
    call jeexin(memass(1:8)//'           .RELR', iret)
    if (iret .gt. 0) then
        call jeveuo(memass(1:8)//'           .RELR', 'L', vk24=v_resu_elem)
        if (v_resu_elem(1) .ne. ' ') then
            nb_matr = nb_matr+1
            matr_name(nb_matr) = memass(1:19)
        end if
    end if
    call asmatr(nb_matr, matr_name, ' ', nume_dof, &
                loadInfoJv, 'ZERO', 'V', 1, matass)
!
! - Factorization of matrix
!
    if (present(ds_algorom_)) then
        if (ds_algorom_%l_rom) then
            call mtdscr(matass)
        else
            call preres(solver, 'V', ierr, maprec, matass, &
                        ibid, -9999)
        end if
    else
        call preres(solver, 'V', ierr, maprec, matass, &
                    ibid, -9999)
    end if
!
    call jedema()
end subroutine
