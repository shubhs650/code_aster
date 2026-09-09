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
subroutine nm1dci(materPara, &
                  option, materPoin, &
                  em, ep, sigm, deps, vim, &
                  sigp, vip, dsde)
!
    use MaterialPara_type
    implicit none
!
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
!
    type(Material_Para), intent(in) :: materPara
    character(len=16) :: option
    character(len=*) :: materPoin
    real(kind=8) :: ep, em
    real(kind=8) :: sigm, deps, vim(2)
    real(kind=8) :: sigp, vip(2), dsde, sieleq
!
! --------------------------------------------------------------------------------------------------
!
!          PLASTICITE VON MISES CINEMATIQUE BILINEAIRE MONODIM
!          ON PEUT AVOIR T0 DIFF TREF
!
! --------------------------------------------------------------------------------------------------
!
! IN  EM        : MODULE D YOUNG MOINS
! IN  EP        : MODULE D YOUNG PLUS
!
! IN  SIGM    : CONTRAINTE AU TEMPS MOINS
! IN  DEPS    : DEFORMATION  TOTALE PLUS - DEFORMATION MOINS
!                       - INCREMENT DEFORMATION THERMIQUE
! IN  VIM     : VARIABLE INTERNES MOINS
! IN  OPTION     : OPTION DE CALCUL
!
! OUT SIG     : CONTRAINTES PLUS
! OUT VIP     : VARIABLE INTERNES PLUS
! OUT DSDE    : DSIG/DEPS
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbProp = 2
    character(len=16), parameter :: propName(nbProp) = (/'D_SIGM_EPSI', &
                                                         'SY         '/)
    real(kind=8) :: propVale(nbProp)
    integer(kind=8) :: propCode(nbProp)
    real(kind=8) :: sige, dp, etm, etp, xp, xm, hm, hp, sigy
!
! --------------------------------------------------------------------------------------------------
!
    call rcvalb(materPara%schemePara%fami, &
                materPara%schemePara%kpg, &
                materPara%schemePara%ksp, &
                '-', materPara%jvMaterCode, &
                materPoin, 'ECRO_LINE', &
                0, ' ', [0.d0], &
                1, propName, propVale, &
                propCode, 1)
    etm = propVale(1)
!
    if (etm .le. 0.) then
        call utmess('F', 'COMPOR1_53', nr=2, valr=[etm, em])
    end if
!
    hm = em*etm/(em-etm)
!
    call rcvalb(materPara%schemePara%fami, &
                materPara%schemePara%kpg, &
                materPara%schemePara%ksp, &
                '+', materPara%jvMaterCode, &
                materPoin, 'ECRO_LINE', &
                0, ' ', [0.d0], &
                nbProp, propName, propVale, &
                propCode, 1)
    etp = propVale(1)
    hp = ep*etp/(ep-etp)
    sigy = propVale(2)
    xm = vim(1)
!     ------------------------------------------------------------------
    sige = ep*(sigm/em+deps)-hp/hm*xm
    sieleq = abs(sige)
!     ------------------------------------------------------------------
!     CALCUL EPSP, P , SIG
!     ------------------------------------------------------------------
    if (option(1:9) .eq. 'FULL_MECA' .or. option(1:9) .eq. 'RAPH_MECA') then
        if (sieleq .le. sigy) then
            vip(2) = 0.d0
            dsde = ep
            dp = 0.d0
            xp = hp/hm*xm
            sigp = ep*(sigm/em+deps)
            vip(1) = xp
        else
            vip(2) = 1.d0
            dp = (sieleq-sigy)/(ep+hp)
            if (option .eq. 'FULL_MECA_ELAS') then
                dsde = ep
            else
                dsde = etp
            end if
            xp = hp/hm*xm+hp*dp*sige/sieleq
            sigp = xp+sigy*sige/sieleq
            vip(1) = xp
        end if
    end if
    if (option(1:10) .eq. 'RIGI_MECA_') then
        if ((vim(2) .lt. 0.5d0) .or. (option .eq. 'RIGI_MECA_ELAS')) then
            dsde = ep
        else
            dsde = etp
        end if
    end if
end subroutine
