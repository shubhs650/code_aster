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

subroutine ccpara(option, &
                  modelZ, materFieldZ, caraElemZ, &
                  resultIn, resultOut, &
                  numeStore, numeStorePrev, isTransient)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8nnem.h"
#include "asterfort/dismoi.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/mearcc.h"
#include "asterfort/mecact.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rsexch.h"
#include "asterfort/vrcins.h"
#include "asterfort/vrcref.h"
#include "asterc/isnnem.h"
!
    character(len=16), intent(in) :: option
    character(len=*), intent(in) :: modelZ, caraElemZ, materFieldZ
    character(len=8), intent(in):: resultIn, resultOut
    integer(kind=8), intent(in) :: numeStore, numeStorePrev
    aster_logical, intent(in) :: isTransient
!
! --------------------------------------------------------------------------------------------------
!
! CALC_CHAMP
!
! Compute ELEM, ELNO and ELGA fields - Create generic input fields
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: timeCmpNb = 6
    character(len=8), parameter :: timeCmpName(timeCmpNb) = (/'INST    ', 'DELTAT  ', 'THETA   ', &
                                                              'KHI     ', 'R       ', 'RHO     '/)
    real(kind=8) :: timeCmmVale(timeCmpNb)
    aster_logical :: lPlate, lCplan, lCompVarc
    integer(kind=8) :: nbParaIn, iret, jvPara, iParaIn
    integer(kind=8) :: optionNume, nh
    integer(kind=8), parameter :: massDiag = 1
    real(kind=8) :: omega2, freq, time, timePrev
    character(len=2) :: chdret
    character(len=8) :: mesh, materField, model, caraElem, answer
    character(len=24) :: sigmElno
    character(len=24) :: paraType
    character(len=24), parameter :: chtime = '&&CCPARA.CH_INST_R'
    character(len=24), parameter :: chfreq = '&&CCPARA.FREQ'
    character(len=24), parameter :: chome2 = '&&CCPARA.OMEGA2'
    character(len=24), parameter :: chharm = '&&CCPARA.NUME_MODE'
    character(len=24), parameter :: chvarc = '&&CCPARA.VARI_INT_N'
    character(len=24), parameter :: chvref = '&&CCPARA.VARI_INT_REF'
    character(len=24), parameter :: chvac2 = '&&CCPARA.VARI_INT_NM1'
    character(len=24), parameter :: chmass = '&&CCPARA.MASS_MECA_D'
    character(len=24), parameter :: chsigf = '&&CCPARA.CHAM_SI2D'
    integer(kind=8), pointer :: cataDescopt(:) => null()
    character(len=24), pointer :: cataLocalis(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    model = modelZ
    materField = materFieldZ
    caraElem = caraElemZ

! - Access to catalog of options
    call jenonu(jexnom('&CATA.OP.NOMOPT', option), optionNume)
    call jeveuo(jexnum('&CATA.OP.DESCOPT', optionNume), 'L', vi=cataDescopt)
    call jeveuo(jexnum('&CATA.OP.LOCALIS', optionNume), 'L', vk24=cataLocalis)
    nbParaIn = cataDescopt(2)

! - Access to mesh
    call dismoi('NOM_MAILLA', model, 'MODELE', repk=mesh)

! - Some specific elements
    call dismoi('EXI_PLAQUE', model, 'MODELE', repk=answer)
    lPlate = answer .eq. "OUI"
    call dismoi('EXI_C_PLAN', model, 'MODELE', repk=answer)
    lCPlan = answer .eq. "OUI"

! - For issue35203
    lCompVarc = ASTER_TRUE
    if (option .eq. "EPSI_ELGA") then
        if (.not. lPlate .and. .not. lCPlan) then
            lCompVarc = ASTER_FALSE
        end if
    end if

! - Get current time (if exist)
    time = 0.d0
    if (isTransient) then
        call rsadpa(resultIn, 'L', 1, 'INST', numeStore, 0, sjv=jvPara)
        time = zr(jvPara)
    end if

! - Get external state variables: reference field only
    call vrcref(model, materField, caraElem, chvref(1:19))

! - Loop on input parameters for input fields
    do iParaIn = 1, nbParaIn
        paraType = cataLocalis(3*iParaIn-1)
        if (paraType .eq. chtime) then
! --------- Input field: time
            if (isTransient) then
                timeCmmVale(1) = time
                timeCmmVale(2) = r8nnem()
                timeCmmVale(3) = r8nnem()
                timeCmmVale(4) = r8nnem()
                timeCmmVale(5) = r8nnem()
                timeCmmVale(6) = r8nnem()
                call mecact('V', chtime, 'MAILLA', mesh, 'INST_R', &
                            ncmp=timeCmpNb, lnomcmp=timeCmpName, vr=timeCmmVale)
            end if
!
        else if (paraType .eq. chfreq) then
! --------- Input field: frequency
            call jenonu(jexnom(resultIn//'           .NOVA', 'FREQ'), iret)
            if (iret .ne. 0) then
                call rsadpa(resultIn, 'L', 1, 'FREQ', numeStore, 0, sjv=jvPara)
                freq = zr(jvPara)
            else
                freq = 1.d0
            end if
            call mecact('V', chfreq, 'MAILLA', mesh, 'FREQ_R', &
                        ncmp=1, nomcmp='FREQ', sr=freq)
!
        else if (paraType .eq. chome2) then
! --------- Input field: pulsation
            call jenonu(jexnom(resultIn//'           .NOVA', 'OMEGA2'), iret)
            if (iret .ne. 0) then
                call rsadpa(resultIn, 'L', 1, 'OMEGA2', numeStore, 0, sjv=jvPara)
                omega2 = zr(jvPara)
            else
                omega2 = 1.0d0
            end if
            call mecact('V', chome2, 'MAILLA', mesh, 'OME2_R', &
                        ncmp=1, nomcmp='OMEG2', sr=omega2)
!
        else if (paraType .eq. chharm) then
! --------- Input field: Fourier mode
            call jenonu(jexnom(resultIn//'           .NOVA', 'NUME_MODE'), iret)
            if (iret .ne. 0) then
                call rsadpa(resultIn, 'L', 1, 'NUME_MODE', numeStore, &
                            0, sjv=jvPara, istop=0)
                nh = zi(jvPara)
                if (nh .ne. isnnem()) then
                    call mecact('V', chharm, 'MAILLA', mesh, 'HARMON', &
                                ncmp=1, nomcmp='NH', si=nh)

                end if
            end if
!
        else if (paraType .eq. chmass) then
! --------- Input field: diagonal mass matrix
            call mecact('V', chmass, 'MAILLA', mesh, 'POSI', &
                        ncmp=1, nomcmp='POS', si=massDiag)
!
        else if (paraType .eq. chvac2) then
! --------- Input field: external state variable at previous time step
            if (isTransient) then
                call rsadpa(resultIn, 'L', 1, 'INST', numeStorePrev, 0, sjv=jvPara)
                timePrev = zr(jvPara)
            else
                timePrev = 0.d0
            end if
            if (lCompVarc) then
                call vrcins(model, materField, caraElem, timePrev, chvac2(1:19), chdret)
            end if
!
        else if (paraType .eq. chvarc) then
! --------- Input field: external state variable at current time step
            if (lCompVarc) then
                call vrcins(model, materField, caraElem, time, chvarc(1:19), chdret)
            end if

        else if (paraType .eq. chsigf) then
! --------- Input field: stress tensor for SIRO_ELEM
            call rsexch(' ', resultIn, 'SIGM_ELNO', numeStore, sigmElno, iret)
            if (iret .ne. 0) then
                call rsexch('F', resultOut, 'SIGM_ELNO', numeStore, sigmElno, iret)
            end if
            call mearcc(option, model, sigmElno, chsigf)
!
        end if
    end do
!
end subroutine
