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
subroutine utdidt(getset, sddisc, ques_type, question, index_, &
                  valr_, vali_, valk_)
!
    implicit none
!
#include "asterf_types.h"
#include "event_def.h"
#include "asterfort/assert.h"
#include "asterfort/jeveuo.h"
!
    character(len=1), intent(in) :: getset
    character(len=19), intent(in) :: sddisc
    character(len=4), intent(in) :: ques_type
    character(len=*), intent(in) :: question
    integer(kind=8), intent(in), optional :: index_
    integer(kind=8), intent(inout), optional :: vali_
    real(kind=8), intent(inout), optional :: valr_
    character(len=*), intent(inout), optional :: valk_
!
! --------------------------------------------------------------------------------------------------
!
! Utility for discretization datastructure
!
! --------------------------------------------------------------------------------------------------
!
! IN  GETSET : 'L' -> LECTURE
!              'E' -> ECRITURE
! IN  SDDISC : SDDISC LOCALE A OP00700
! IN  TYPQUE : TYPE DE DEMANDE (LIST, ECHE OU ADAP)
! IN  IOCC   : NUMERO OCCURRENCE (POUR ECHEC/ADAPT)
! IN  QUEST  : QUESTION
! I/O VALI   : VALEUR ENTIERE
! I/O VALR   : VALEUR REELLE
! I/O VALK   : VALEUR CHAINE
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: iechec, i_adap
    character(len=16) :: valk
    real(kind=8) :: valr
    integer(kind=8) :: vali
    character(len=24) :: sddiscLinfName
    real(kind=8), pointer :: sddiscLinf(:) => null()
    character(len=24) :: sddiscEevrName
    real(kind=8), pointer :: sddiscEevr(:) => null()
    character(len=24) :: sddiscEevkName
    character(len=16), pointer :: sddiscEevk(:) => null()
    character(len=24) :: sddiscEsurName
    real(kind=8), pointer :: sddiscEsur(:) => null()
    character(len=24) :: sddiscAevrName
    real(kind=8), pointer :: sddiscAevr(:) => null()
    character(len=24) :: sddiscAtprName
    real(kind=8), pointer :: sddiscAtpr(:) => null()
    character(len=24) :: sddiscAtpkName
    character(len=16), pointer :: sddiscAtpk(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(ques_type .eq. 'LIST' .or. ques_type .eq. 'ECHE' .or. ques_type .eq. 'ADAP')
    ASSERT(getset .eq. 'L' .or. getset .eq. 'E')

! - Initializations
    if (getset .eq. 'L') then
        valk = ' '
        vali = 0
        valr = 0.d0
    else
        if (present(valk_)) then
            valk = valk_
        end if
        if (present(vali_)) then
            vali = vali_
        end if
        if (present(valr_)) then
            valr = valr_
        end if
    end if
!
! - Questions about LIST
!
    if (ques_type .eq. 'LIST') then
        sddiscLinfName = sddisc(1:19)//'.LINF'
        call jeveuo(sddiscLinfName, getset, vr=sddiscLinf)
        if (question .eq. 'METHODE') then
            if (getset .eq. 'L') then
                vali = nint(sddiscLinf(1))
                if (vali .eq. 1) valk = 'MANUEL'
                if (vali .eq. 2) valk = 'AUTO'
            else if (getset .eq. 'E') then
                if (valk .eq. 'MANUEL') then
                    sddiscLinf(1) = 1
                else if (valk .eq. 'AUTO') then
                    sddiscLinf(1) = 2
                else
                    ASSERT(ASTER_FALSE)
                end if
            end if
        else if (question .eq. 'PAS_MINI') then
            if (getset .eq. 'L') then
                valr = sddiscLinf(2)
            else if (getset .eq. 'E') then
                sddiscLinf(2) = valr
            end if
        else if (question .eq. 'PAS_MAXI') then
            if (getset .eq. 'L') then
                valr = sddiscLinf(3)
            else if (getset .eq. 'E') then
                sddiscLinf(3) = valr
            end if
        else if (question .eq. 'NB_PAS_MAXI') then
            if (getset .eq. 'L') then
                vali = nint(sddiscLinf(4))
            else if (getset .eq. 'E') then
                sddiscLinf(4) = vali
            end if
        else if (question .eq. 'DTMIN') then
            if (getset .eq. 'L') then
                valr = sddiscLinf(5)
            else if (getset .eq. 'E') then
                sddiscLinf(5) = valr
            end if
        else if (question .eq. 'DT-') then
            if (getset .eq. 'L') then
                valr = sddiscLinf(6)
            else if (getset .eq. 'E') then
                sddiscLinf(6) = valr
            end if
        else if (question .eq. 'EXIS_DECOUPE') then
            if (getset .eq. 'L') then
                vali = nint(sddiscLinf(7))
                if (vali .eq. 0) valk = 'NON'
                if (vali .eq. 1) valk = 'OUI'
            else if (getset .eq. 'E') then
                if (valk .eq. 'NON') then
                    sddiscLinf(7) = 0
                else if (valk .eq. 'OUI') then
                    sddiscLinf(7) = 1
                else
                    write (6, *) 'VALK: ', valk
                    ASSERT(ASTER_FALSE)
                end if
            end if
        else if (question .eq. 'NBINST') then
            if (getset .eq. 'L') then
                vali = nint(sddiscLinf(8))
            else if (getset .eq. 'E') then
                sddiscLinf(8) = vali
            end if
        else if (question .eq. 'NECHEC') then
            if (getset .eq. 'L') then
                vali = nint(sddiscLinf(9))
            else if (getset .eq. 'E') then
                sddiscLinf(9) = vali
            end if
        else if (question .eq. 'NADAPT') then
            if (getset .eq. 'L') then
                vali = nint(sddiscLinf(10))
            else if (getset .eq. 'E') then
                sddiscLinf(10) = vali
            end if
        else
            ASSERT(ASTER_FALSE)
        end if

! - Questions about ECHEC
    else if (ques_type .eq. 'ECHE') then
        sddiscEevrName = sddisc(1:19)//'.EEVR'
        sddiscEevkName = sddisc(1:19)//'.EEVK'
        sddiscEsurName = sddisc(1:19)//'.ESUR'
        call jeveuo(sddiscEevrName, getset, vr=sddiscEevr)
        call jeveuo(sddiscEevkName, getset, vk16=sddiscEevk)
        call jeveuo(sddiscEsurName, getset, vr=sddiscEsur)
        if (present(index_)) then
            iechec = index_
        end if
        if (question .eq. 'VERIF_EVEN') then
            if (getset .eq. 'L') then
                vali = nint(sddiscEevr(SIZE_LEEVR*(iechec-1)+3))
                if (vali .eq. 0) valk = 'OUI'
                if (vali .eq. 1) valk = 'NON'
            else if (getset .eq. 'E') then
                if (valk .eq. 'OUI') then
                    sddiscEevr(SIZE_LEEVR*(iechec-1)+3) = 0
                else if (valk .eq. 'NON') then
                    sddiscEevr(SIZE_LEEVR*(iechec-1)+3) = 1
                else
                    ASSERT(ASTER_FALSE)
                end if
            end if

! ----- Parameters for DELTA_GRANDEUR
        else if (question .eq. 'NOM_CHAM') then
            if (getset .eq. 'L') then
                valk = sddiscEevk(SIZE_LEEVK*(iechec-1)+1)
            else if (getset .eq. 'E') then
                sddiscEevk(SIZE_LEEVK*(iechec-1)+1) = valk
            end if
        else if (question .eq. 'NOM_CMP') then
            if (getset .eq. 'L') then
                valk = sddiscEevk(SIZE_LEEVK*(iechec-1)+2)
            else if (getset .eq. 'E') then
                sddiscEevk(SIZE_LEEVK*(iechec-1)+2) = valk
            end if
        else if (question .eq. 'CRIT_COMP') then
            if (getset .eq. 'L') then
                valk = sddiscEevk(SIZE_LEEVK*(iechec-1)+3)
            else if (getset .eq. 'E') then
                sddiscEevk(SIZE_LEEVK*(iechec-1)+3) = valk
            end if
        else if (question .eq. 'VALE_REF') then
            if (getset .eq. 'L') then
                valr = sddiscEevr(SIZE_LEEVR*(iechec-1)+5)
            else if (getset .eq. 'E') then
                sddiscEevr(SIZE_LEEVR*(iechec-1)+5) = valr
            end if

! ----- Parameters for INTERPENETRATION
        else if (question .eq. 'PENE_MAXI') then
            if (getset .eq. 'L') then
                valr = sddiscEevr(SIZE_LEEVR*(iechec-1)+6)
            else if (getset .eq. 'E') then
                sddiscEevr(SIZE_LEEVR*(iechec-1)+6) = valr
            end if

! ----- Parameters for RESI_MAXI
        else if (question .eq. 'RESI_GLOB_MAXI') then
            if (getset .eq. 'L') then
                valr = sddiscEevr(SIZE_LEEVR*(iechec-1)+7)
            else if (getset .eq. 'E') then
                sddiscEevr(SIZE_LEEVR*(iechec-1)+7) = valr
            end if

! ----- Parameters for DECOUPE
        else if (question .eq. 'SUBD_METHODE') then
            if (getset .eq. 'L') then
                vali = nint(sddiscEsur(SIZE_LESUR*(iechec-1)+1))
                if (vali .eq. 0) valk = 'AUCUNE'
                if (vali .eq. 1) valk = 'MANUEL'
                if (vali .eq. 2) valk = 'AUTO'
            else if (getset .eq. 'E') then
                if (valk .eq. 'AUCUNE') then
                    sddiscEsur(SIZE_LESUR*(iechec-1)+1) = 0
                else if (valk .eq. 'MANUEL') then
                    sddiscEsur(SIZE_LESUR*(iechec-1)+1) = 1
                else if (valk .eq. 'AUTO') then
                    sddiscEsur(SIZE_LESUR*(iechec-1)+1) = 2
                else
                    ASSERT(ASTER_FALSE)
                end if
            end if
        else if (question .eq. 'SUBD_METHODE_AUTO') then
            if (getset .eq. 'L') then
                vali = nint(sddiscEsur(SIZE_LESUR*(iechec-1)+10))
                if (vali .eq. 2) valk = 'EXTRAPOLE'
            else if (getset .eq. 'E') then
                if (valk .eq. 'EXTRAPOLE') then
                    sddiscEsur(SIZE_LESUR*(iechec-1)+10) = 21
                else
                    ASSERT(ASTER_FALSE)
                end if
            end if
        else if (question .eq. 'SUBD_PAS') then
            if (getset .eq. 'L') then
                vali = nint(sddiscEsur(SIZE_LESUR*(iechec-1)+2))
            else if (getset .eq. 'E') then
                sddiscEsur(SIZE_LESUR*(iechec-1)+2) = vali
            end if
        else if (question .eq. 'SUBD_PAS_MINI') then
            if (getset .eq. 'L') then
                valr = sddiscEsur(SIZE_LESUR*(iechec-1)+3)
            else if (getset .eq. 'E') then
                sddiscEsur(SIZE_LESUR*(iechec-1)+3) = valr
            end if
        else if (question .eq. 'SUBD_NIVEAU') then
            if (getset .eq. 'L') then
                vali = nint(sddiscEsur(SIZE_LESUR*(iechec-1)+4))
            else if (getset .eq. 'E') then
                sddiscEsur(SIZE_LESUR*(iechec-1)+4) = vali
            end if
        else if (question .eq. 'SUBD_INST') then
            if (getset .eq. 'L') then
                valr = sddiscEsur(SIZE_LESUR*(iechec-1)+5)
            else if (getset .eq. 'E') then
                sddiscEsur(SIZE_LESUR*(iechec-1)+5) = valr
            end if
        else if (question .eq. 'SUBD_RATIO') then
            if (getset .eq. 'L') then
                vali = nint(sddiscEsur(SIZE_LESUR*(iechec-1)+9))
            else if (getset .eq. 'E') then
                sddiscEsur(SIZE_LESUR*(iechec-1)+9) = vali
            end if
        else
            ASSERT(ASTER_FALSE)
        end if

! - Questions about ADAPTATION
    else if (ques_type .eq. 'ADAP') then
        sddiscAevrName = sddisc(1:19)//'.AEVR'
        sddiscAtprName = sddisc(1:19)//'.ATPR'
        sddiscAtpkName = sddisc(1:19)//'.ATPK'
        call jeveuo(sddiscAevrName, getset, vr=sddiscAevr)
        call jeveuo(sddiscAtprName, getset, vr=sddiscAtpr)
        call jeveuo(sddiscAtpkName, getset, vk16=sddiscAtpk)
        i_adap = index_
        if (question .eq. 'NB_INCR_SEUIL') then
            if (getset .eq. 'L') then
                vali = nint(sddiscAevr(SIZE_LAEVR*(i_adap-1)+2))
            else if (getset .eq. 'E') then
                sddiscAevr(SIZE_LAEVR*(i_adap-1)+2) = 1
            end if
        else if (question .eq. 'NOM_PARA') then
            if (getset .eq. 'L') then
                vali = nint(sddiscAevr(SIZE_LAEVR*(i_adap-1)+3))
                if (vali .eq. 1) valk = 'NB_ITER_NEWT'
                if (vali .eq. 2) valk = 'DP'
            else if (getset .eq. 'E') then
                if (valk .eq. 'NB_ITER_NEWT') then
                    sddiscAevr(SIZE_LAEVR*(i_adap-1)+3) = 1
                else if (valk .eq. 'SEUIL_AVEC_FORMU') then
                    sddiscAevr(SIZE_LAEVR*(i_adap-1)+3) = 2
                else
                    ASSERT(ASTER_FALSE)
                end if
            end if
        else if (question .eq. 'CRIT_COMP') then
            if (getset .eq. 'L') then
                vali = nint(sddiscAevr(SIZE_LAEVR*(i_adap-1)+4))
                if (vali .eq. 1) valk = 'LT'
                if (vali .eq. 2) valk = 'GT'
                if (vali .eq. 3) valk = 'LE'
                if (vali .eq. 4) valk = 'GE'
            else if (getset .eq. 'E') then
                if (valk .eq. 'LT') then
                    sddiscAevr(SIZE_LAEVR*(i_adap-1)+4) = 1
                else if (valk .eq. 'GT') then
                    sddiscAevr(SIZE_LAEVR*(i_adap-1)+4) = 2
                else if (valk .eq. 'LE') then
                    sddiscAevr(SIZE_LAEVR*(i_adap-1)+4) = 3
                else if (valk .eq. 'GE') then
                    sddiscAevr(SIZE_LAEVR*(i_adap-1)+4) = 4
                else
                    ASSERT(ASTER_FALSE)
                end if
            end if
        else if (question .eq. 'VALE') then
            if (getset .eq. 'L') then
                valr = sddiscAevr(SIZE_LAEVR*(i_adap-1)+5)
                vali = nint(sddiscAevr(SIZE_LAEVR*(i_adap-1)+5))
            else if (getset .eq. 'E') then
                sddiscAevr(SIZE_LAEVR*(i_adap-1)+5) = valr
            end if
        else if (question .eq. 'NB_EVEN_OK') then
            if (getset .eq. 'L') then
                valr = sddiscAevr(SIZE_LAEVR*(i_adap-1)+6)
                vali = nint(valr)
            else if (getset .eq. 'E') then
                sddiscAevr(SIZE_LAEVR*(i_adap-1)+6) = vali
            end if
        else if (question .eq. 'PCENT_AUGM') then
            if (getset .eq. 'L') then
                valr = sddiscAtpr(SIZE_LATPR*(i_adap-1)+2)
            else if (getset .eq. 'E') then
                sddiscAtpr(SIZE_LATPR*(i_adap-1)+2) = valr
            end if
        else if (question .eq. 'VALE_REF') then
            if (getset .eq. 'L') then
                valr = sddiscAtpr(SIZE_LATPR*(i_adap-1)+3)
            else if (getset .eq. 'E') then
                sddiscAtpr(SIZE_LATPR*(i_adap-1)+3) = valr
            end if
        else if (question .eq. 'NU_CMP') then
            if (getset .eq. 'L') then
                valr = sddiscAtpr(SIZE_LATPR*(i_adap-1)+4)
                vali = nint(valr)
            else if (getset .eq. 'E') then
                sddiscAtpr(SIZE_LATPR*(i_adap-1)+4) = vali
            end if
        else if (question .eq. 'NB_ITER_NEWTON_REF') then
            if (getset .eq. 'L') then
                valr = sddiscAtpr(SIZE_LATPR*(i_adap-1)+5)
                vali = nint(valr)
            else if (getset .eq. 'E') then
                sddiscAtpr(SIZE_LATPR*(i_adap-1)+5) = vali
            end if
        else if (question .eq. 'NOM_CHAM') then
            if (getset .eq. 'L') then
                valk = sddiscAtpk(SIZE_LATPK*(i_adap-1)+2)
            else if (getset .eq. 'E') then
                sddiscAtpk(SIZE_LATPK*(i_adap-1)+2) = valk
            end if
        else if (question .eq. 'NOM_CMP') then
            if (getset .eq. 'L') then
                valk = sddiscAtpk(SIZE_LATPK*(i_adap-1)+3)
            else if (getset .eq. 'E') then
                sddiscAtpk(SIZE_LATPK*(i_adap-1)+3) = valk
            end if
        else
            ASSERT(ASTER_FALSE)
        end if
    end if
!
    if (present(vali_)) then
        vali_ = vali
    end if
    if (present(valr_)) then
        valr_ = valr
    end if
    if (present(valk_)) then
        valk_ = valk
    end if
!
end subroutine
