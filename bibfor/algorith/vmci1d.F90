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
subroutine vmci1d(materPara, &
                  option, materPoin, &
                  em, ep, &
                  sigm, deps, vim, &
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
    real(kind=8) :: sigm, deps, vim(*)
    real(kind=8) :: sigp, vip(*), dsde
!
! --------------------------------------------------------------------------------------------------
!
!           PLASTICITE VON MISES CINEMATIQUE LINEAIRE EN 1D
!              FORTEMENT INSPIRE DE NM1DCI
!
! --------------------------------------------------------------------------------------------------
!
!        EM     : MODULE D YOUNG MOINS
!        EP     : MODULE D YOUNG PLUS
!        SIGM   : CONTRAINTE AU TEMPS MOINS
!        DEPS   : DEFORMATION TOTALE PLUS - DEFORMATION MOINS
!                       - INCREMENT DEFORMATION THERMIQUE
!        VIM    : VARIABLE INTERNES MOINS
!        OPTION : OPTION DE CALCUL
!  OUT
!        SIGP   : CONTRAINTES PLUS
!        VIP    : VARIABLE INTERNES PLUS
!        DSDE   : DSIG/DEPS
! --------------------------------------------------------------------------------------------------
!     Variables internes
!       icels : critère sigma
!       icelu : critère epsi
!       iepsq : déformation équivalente
!       iplas : indicateur plastique
!       idiss : dissipation plastique
!       iwthe : dissipation thermodynamique
!       i..m  : ecrouissage cinematique
!
! --------------------------------------------------------------------------------------------------
!
!   index des variables internes
!           'CRITSIG', 'CRITEPS', 'EPSPEQ', 'INDIPLAS', 'DISSIP', 'DISSTHER',
!           'XCINXX',  'XCINYY',  'XCINZZ', 'XCINXY', 'XCINXZ', 'XCINYZ',
    integer(kind=8), parameter :: icels = 1, icelu = 2, iepsq = 3, iplas = 4, idiss = 5, iwthe = 6
    integer(kind=8), parameter :: ixxm = 7
    integer(kind=8), parameter :: nbvari = 12
!
    integer(kind=8), parameter :: nbProp = 4
    character(len=16), parameter :: propName(nbProp) = (/'D_SIGM_EPSI', 'SY         ', &
                                                         'SIGM_LIM   ', 'EPSI_LIM   '/)
    real(kind=8) :: propVale(nbProp)
    integer(kind=8) :: propCode(nbProp)
    real(kind=8) :: sigy, sieleq, sige, dp, etm, etp, xp, xm, hm, hp, sgels, epelu
    character(len=16) :: valkm(3)
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
    hm = em*etm/(em-etm)

    call rcvalb(materPara%schemePara%fami, &
                materPara%schemePara%kpg, &
                materPara%schemePara%ksp, &
                '+', materPara%jvMaterCode, &
                materPoin, 'ECRO_LINE', &
                0, ' ', [0.d0], &
                nbProp, propName, propVale, &
                propCode, 1)

!   vérification que SIGM_LIM, EPSI_LIM sont présents
    if (propCode(3)+propCode(4) .ne. 0) then
        valkm(1) = 'VMIS_CINE_GC'
        valkm(2) = propName(3)
        valkm(3) = propName(4)
        call utmess('F', 'COMPOR1_76', nk=3, valk=valkm)
    end if
    etp = propVale(1)
    sigy = propVale(2)
    sgels = propVale(3)
    epelu = propVale(4)
!
    hp = ep*etp/(ep-etp)
    xm = vim(ixxm)
!
    sige = ep*(sigm/em+deps)-hp*xm/hm
    sieleq = abs(sige)
!
! --------------------------------------------------------------------------------------------------
!   calcul : EPSP, P , SIG
    if ((option(1:9) .eq. 'FULL_MECA') .or. (option(1:9) .eq. 'RAPH_MECA')) then
        vip(1:nbvari) = vim(1:nbvari)
        if (sieleq .le. sigy) then
            vip(iplas) = 0.d0
            dsde = ep
            dp = 0.d0
            xp = hp*xm/hm
            sigp = ep*(sigm/em+deps)
            vip(ixxm) = xp
            vip(icelu) = (sigm/em+deps)/epelu
            vip(iepsq) = (sigm/em+deps)
        else
            vip(iplas) = 1.d0
            dp = (sieleq-sigy)/(ep+hp)
            if (option .eq. 'FULL_MECA_ELAS') then
                dsde = ep
            else
                dsde = etp
            end if
            xp = hp*xm/hm+hp*dp*sige/sieleq
            sigp = xp+sigy*sige/sieleq
            vip(ixxm) = xp
            vip(icelu) = ((sigp-sigy)/etp+sigy/ep)/epelu
            vip(iepsq) = ((sigp-sigy)/etp+sigy/ep)
        end if
        vip(icels) = sigp/sgels
!       dissipation thermodynamique
        vip(iwthe) = vim(iwthe)+sigy*dp
!       dissipation irréversible
        vip(idiss) = vim(idiss)+(dsde*deps-(sigp-sigm))*deps/2.0d0
    end if
    if (option(1:10) .eq. 'RIGI_MECA_') then
        if ((vim(iplas) .lt. 0.5d0) .or. (option .eq. 'RIGI_MECA_ELAS')) then
            dsde = ep
        else
            dsde = etp
        end if
    end if
end subroutine
