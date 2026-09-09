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
subroutine nm1dis(materPara, &
                  option, relaComp, materPoin, &
                  em, ep, sigm, deps, vim, &
                  sigp, vip, dsde)
!
    use MaterialPara_type
    use Behaviour_type
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/rcfonc.h"
#include "asterfort/rctrac.h"
#include "asterfort/rctype.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/utmess.h"
!
    type(Material_Para), intent(in) :: materPara
    character(len=16) :: option, relaComp
    character(len=*) :: materPoin
    real(kind=8) :: em, ep, sigm, deps, vim(*), sigy
    real(kind=8) :: vip(*), sigp, dsde
!
! --------------------------------------------------------------------------------------------------
!
!          PLASTICITE VON MISES ISOTROPE BILINEAIRE MONODIM
!
! --------------------------------------------------------------------------------------------------
!
!   em      : module d'Young à t-
!   ep      : module d'Young à t+
!   sigm    : contrainte à t-
!   vim     : variables internes à t-
!   deps    : déformation totale plus - déformation moins - incrément déformation thermique
!   sigp    : contraintes à t+
!   vip     : variables internes à t+
!   epsp    : deformation  plastique plus
!   dsde    : dsig/deps
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbProp = 2
    character(len=16), parameter :: propName(nbProp) = (/'D_SIGM_EPSI', &
                                                         'SY         '/)
    real(kind=8) :: propVale(nbProp)
    integer(kind=8) :: propCode(nbProp)
    integer(kind=8), parameter :: nbPara = 1
    character(len=16), parameter :: paraName(nbPara) = (/'TEMP'/)
    real(kind=8) :: paraVale
    integer(kind=8) :: jprolm, jvalem, nbvalm, nbvalp, jprolp, jvalep, iret
    real(kind=8) :: rprim, rm, sige, airerp, sieleq, rp, dp, nu, asige, pm, et
    character(len=8) :: tracParaType
    real(kind=8) :: tracParaVale
!
! --------------------------------------------------------------------------------------------------
!
    pm = vim(1)
    et = 0.0d0

!   caractéristiques écrouissage linéaire
    if ((relaComp .eq. 'VMIS_ISOT_LINE') .or. (relaComp .eq. 'GRILLE_ISOT_LINE')) then
        call rcvalb(materPara%schemePara%fami, &
                    materPara%schemePara%kpg, &
                    materPara%schemePara%ksp, &
                    '+', materPara%jvMaterCode, &
                    materPoin, 'ECRO_LINE', &
                    0, ' ', [0.d0], &
                    1, propName, propVale, propCode, 1)
        call rcvalb(materPara%schemePara%fami, &
                    materPara%schemePara%kpg, &
                    materPara%schemePara%ksp, &
                    '+', materPara%jvMaterCode, &
                    materPoin, 'ECRO_LINE', &
                    0, ' ', [0.d0], &
                    1, propName(2), propVale(2), propCode(2), 0)
        if (propCode(2) .ne. 0) then
            propVale(2) = 0.d0
        end if
        et = propVale(1)
        sigy = propVale(2)
        rprim = ep*et/(ep-et)
        rm = rprim*vim(1)+sigy

!   caractéristiques écrouissage donné par courbe de traction
    else if (relaComp .eq. 'VMIS_ISOT_TRAC') then
        call rcvarc(' ', 'TEMP', '-', &
                    materPara%schemePara%fami, &
                    materPara%schemePara%kpg, &
                    materPara%schemePara%ksp, &
                    paraVale, iret)
        call rctype(materPara%jvMaterCode, &
                    nbPara, paraName, [paraVale], &
                    tracParaVale, tracParaType, &
                    materi=materPoin)
        if ((tracParaType .eq. 'TEMP') .and. (iret .eq. 1)) then
            call utmess('F', 'COMPOR5_5', sk=tracParaType)
        end if
        call rctrac(materPara%jvMaterCode, &
                    1, 'SIGM', tracParaVale, &
                    jprolm, jvalem, nbvalm, em, &
                    materi=materPoin)
        call rcvarc(' ', 'TEMP', '+', &
                    materPara%schemePara%fami, &
                    materPara%schemePara%kpg, &
                    materPara%schemePara%ksp, &
                    paraVale, iret)
        call rctype(materPara%jvMaterCode, &
                    nbPara, paraName, [paraVale], &
                    tracParaVale, tracParaType, &
                    materi=materPoin)
        if ((tracParaType .eq. 'TEMP') .and. (iret .eq. 1)) then
            call utmess('F', 'COMPOR5_5', sk=tracParaType)
        end if

        call rctrac(materPara%jvMaterCode, 1, &
                    'SIGM', tracParaVale, &
                    jprolp, jvalep, nbvalp, ep, &
                    materi=materPoin)
        call rcfonc('S', 1, jprolp, jvalep, nbvalp, sigy=sigy)
        call rcfonc('V', 1, jprolp, jvalep, nbvalp, p=vim(1), rp=rm, rprim=rprim, airerp=airerp)
        et = rprim
    else
        ASSERT(ASTER_FALSE)
    end if

!   estimation élastique
    sige = ep*(sigm/em+deps)
    sieleq = abs(sige)

!   calcul epsp, p , sig
    if (option(1:9) .eq. 'FULL_MECA' .or. option(1:9) .eq. 'RAPH_MECA') then
        if (sieleq .le. rm) then
            dp = 0.d0
            sigp = sige
            dsde = ep
            vip(2) = 0.d0
            vip(1) = vim(1)
            sigp = sige
        else
            vip(2) = 1.d0
            if ((relaComp .eq. 'VMIS_ISOT_LINE') .or. (relaComp .eq. 'GRILLE_ISOT_LINE')) then
                dp = abs(sige)-rm
                dp = dp/(rprim+ep)
                rp = sigy+rprim*(pm+dp)
                if (option .eq. 'FULL_MECA_ELAS') then
                    dsde = ep
                else
                    dsde = et
                end if
            else
                nu = 0.5d0
                asige = abs(sige)
                call rcfonc('E', 1, jprolp, jvalep, nbvalp, e=ep, nu=nu, p=vim(1), rp=rp, &
                            rprim=rprim, airerp=airerp, sieleq=asige, dp=dp)
                if (option .eq. 'FULL_MECA_ELAS') then
                    dsde = ep
                else
                    dsde = ep*rprim/(ep+rprim)
                end if
            end if
            vip(1) = vim(1)+dp
            sigp = sige/(1.d0+ep*dp/rp)
        end if
    end if
    if (option(1:10) .eq. 'RIGI_MECA_') then
        if ((vim(2) .lt. 0.5d0) .or. (option .eq. 'RIGI_MECA_ELAS')) then
            dsde = ep
        else
            dsde = et
        end if
    end if
end subroutine
