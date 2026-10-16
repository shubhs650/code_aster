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
subroutine rcvarc(arret, varc_name_, poum, &
                  famiZ, kpg, ksp, &
                  varc_vale, iret)
!
    use calcul_module, only: ca_decala_, ca_iactif_, ca_iredec_, &
                             ca_jfpgl_, ca_jvcnom_, ca_nbcvrc_, ca_nfpg_, &
                             ca_nomte_, ca_option_, &
                             ca_td1_, ca_tf1_, ca_timed1_, ca_timef1_, &
                             ca_ctempl_, ca_ctempr_, ca_ctempm_, ca_ctempp_, &
                             ca_cpcapl_, ca_cpcapm_, ca_cpcapp_
!
    implicit none
!
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterc/r8nnem.h"
#include "asterfort/assert.h"
#include "asterfort/rcvarp.h"
#include "asterfort/tecach.h"
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"
!
    character(len=1), intent(in) :: arret
    character(len=*), intent(in) :: varc_name_
    character(len=*), intent(in) :: poum
    character(len=*), intent(in) :: famiZ
    integer(kind=8), intent(in) :: kpg, ksp
    real(kind=8), intent(out) :: varc_vale
    integer(kind=8), intent(out) :: iret
!
! --------------------------------------------------------------------------------------------------
!
! Material - External state variables (VARC)
!
! Get value of external state variable
!
! --------------------------------------------------------------------------------------------------
!
! In  arret            : in case of problem
!                     ' '   no message and output return code (iret)
!                     'F'   fatal error
! In  varc_name        : name of external state variable
! In  poum             : parameters evaluation
!                     '-' at beginning of step time
!                     '+' at end of step time
!                   'REF' for reference value
! In  fami             : Gauss family for integration point rule
! In  kpg              : current point gauss
! In  ksp              : current "sous-point" gauss
! Out varc_vale        : value of external state variable
! Out iret             : code if error
!                      0    if OK
!                      1    if not found
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: varc_name, fami
    integer(kind=8) :: nb2vrc, k, i
    integer(kind=8) :: varc_indx, ibid, nbsp, kpgvrc
    integer(kind=8) :: iadzi, iazk24, kpgmat
    integer(kind=8) :: itabm(7), itabp(7), itabr(7)
    character(len=24) :: valk(4)
    real(kind=8) :: rundf, valvrm, valvrp
!
! --------------------------------------------------------------------------------------------------
!
    rundf = r8nnem()
    varc_name = varc_name_
    fami = famiZ
!
! - From SIMU_POINT_MAT
!
    if (ca_iactif_ .eq. 2) then
        ASSERT(fami .eq. 'PMAT')
        call rcvarp(arret, varc_name, poum, varc_vale, iret)
        goto 999
    end if

! - For coupled problems
    if (ca_ctempl_ .eq. 1) then
        if (varc_name .eq. 'TEMP') then
            if (poum .eq. '-') then
                varc_vale = ca_ctempm_
            elseif (poum .eq. '+') then
                varc_vale = ca_ctempp_
            elseif (poum .eq. 'REF') then
                varc_vale = ca_ctempr_
            else
                ASSERT(.false.)
            end if
            iret = 0
            goto 999
        end if
    end if
    if (ca_cpcapl_ .eq. 1) then
        if (varc_name .eq. 'PCAP') then
            if (poum .eq. '-') then
                varc_vale = ca_cpcapm_
            elseif (poum .eq. '+') then
                varc_vale = ca_cpcapp_
            else
                ASSERT(.false.)
            end if
            iret = 0
            goto 999
        end if
    end if

!
! - No external state variable
!
    if (ca_nbcvrc_ .eq. 0) then
        goto 998
    end if
!
! - Get index of external state variable
!
    varc_indx = indik8(zk8(ca_jvcnom_), varc_name, 1, ca_nbcvrc_)
    if (varc_indx .eq. 0) then
        iret = 1
        if (arret .eq. ' ') then
            varc_vale = rundf
            goto 999
        else
            call tecael(iadzi, iazk24)
            valk(1) = varc_name
            valk(2) = poum
            call utmess('F', 'CALCUL_26', nk=2, valk=valk, si=zi(iadzi))
        end if
    end if
!
! - Get index of gauss point in MATER integration rule
!
    k = indik8(zk8(ca_jfpgl_), fami, 1, ca_nfpg_)

    if (k .eq. 0) then
        do i = 1, ca_nfpg_
            WRITE (6, *) "RCVRAC <", fami(1:8), ">:", i, zk8(ca_jfpgl_-1+i)
        end do
        valk(1) = varc_name
        valk(2) = fami
        valk(3) = ca_option_
        valk(4) = ca_nomte_
        call utmess('F', 'CALCUL_31', nk=4, valk=valk)
    end if
    kpgmat = ca_decala_(k)+kpg
!
! - Get information about current element.
!   Keep cache metadata local to this call so rcvarc stays reentrant.
!
    if (poum .eq. '-' .or. (poum .eq. '+' .and. ca_iredec_ .eq. 1)) then
        if (arret .ne. ' ') then
            call tecach('OOO', 'PVARCMR', 'L', ibid, nval=7, itab=itabm)
        else
            call tecach('NNN', 'PVARCMR', 'L', iret, nval=7, itab=itabm)
            if (iret .ne. 0) goto 998
        end if
    end if

    if (poum .eq. '+' .or. (poum .eq. '-' .and. ca_iredec_ .eq. 1)) then
        if (arret .ne. ' ') then
            call tecach('OOO', 'PVARCPR', 'L', ibid, nval=7, itab=itabp)
        else
            call tecach('NNN', 'PVARCPR', 'L', iret, nval=7, itab=itabp)
            if (iret .ne. 0) goto 998
        end if
    end if

    if (poum .eq. 'REF') then
        if (arret .ne. ' ') then
            call tecach('OOO', 'PVARCRR', 'L', ibid, nval=7, itab=itabr)
        else
            call tecach('NNN', 'PVARCRR', 'L', iret, nval=7, itab=itabr)
            if (iret .ne. 0) goto 998
        end if
    end if
!
! - Get value
!
    if (poum .eq. 'REF') then
        nb2vrc = itabr(6)
        nbsp = itabr(7)
        kpgvrc = (kpgmat-1)*nbsp+ksp
        varc_vale = zr(itabr(1)-1+(kpgvrc-1)*nb2vrc+varc_indx)

    else if (poum .eq. '+' .and. ca_iredec_ .eq. 0) then
        nb2vrc = itabp(6)
        nbsp = itabp(7)
        kpgvrc = (kpgmat-1)*nbsp+ksp
        varc_vale = zr(itabp(1)-1+(kpgvrc-1)*nb2vrc+varc_indx)

    else if (poum .eq. '-' .and. ca_iredec_ .eq. 0) then
        nb2vrc = itabm(6)
        nbsp = itabm(7)
        kpgvrc = (kpgmat-1)*nbsp+ksp
        varc_vale = zr(itabm(1)-1+(kpgvrc-1)*nb2vrc+varc_indx)

    else if (ca_iredec_ .eq. 1) then
        nb2vrc = itabm(6)
        nbsp = itabm(7)
        kpgvrc = (kpgmat-1)*nbsp+ksp
        valvrm = zr(itabm(1)-1+(kpgvrc-1)*nb2vrc+varc_indx)

        nb2vrc = itabp(6)
        nbsp = itabp(7)
        kpgvrc = (kpgmat-1)*nbsp+ksp
        valvrp = zr(itabp(1)-1+(kpgvrc-1)*nb2vrc+varc_indx)

        if ((.not. isnan(valvrm)) .and. (.not. isnan(valvrp))) then
            if (poum .eq. '-') then
                varc_vale = valvrm+(ca_td1_-ca_timed1_)*(valvrp-valvrm)/(ca_timef1_- &
                                                                         ca_timed1_)
            else if (poum .eq. '+') then
                varc_vale = valvrm+(ca_tf1_-ca_timed1_)*(valvrp-valvrm)/(ca_timef1_- &
                                                                         ca_timed1_)
            else
                ASSERT(.false.)
            end if
        else
            varc_vale = rundf
        end if

    else
        ASSERT(.false.)
    end if
!
    iret = 0
    if (isnan(varc_vale)) then
        iret = 1
    end if
!
! - Manage error
!
    if (iret .eq. 1) then
        if (arret .eq. ' ') then
            varc_vale = rundf
        else
            call tecael(iadzi, iazk24)
            valk(1) = varc_name
            valk(2) = "+"
            call utmess('F', 'CALCUL_26', nk=2, valk=valk, si=zi(iadzi))
        end if
    end if
    goto 999
!
998 continue
    if (arret .eq. ' ') then
        varc_vale = rundf
        iret = 1
    else
        call tecael(iadzi, iazk24)
        valk(1) = zk24(iazk24-1+3)
        call utmess('F', 'CALCUL_32', sk=valk(1))
    end if
!
999 continue
!
end subroutine
