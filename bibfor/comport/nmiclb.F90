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
subroutine nmiclb(materPara, &
                  option, relaComp, carcri, &
                  xlong0, aire, tmoins, tplus, &
                  dlong0, effnom, vim, effnop, vip, &
                  klv, fono, epsm, codret)
!
    use MaterialPara_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/lcimpl.h"
#include "asterfort/nm1dci.h"
#include "asterfort/nm1dco.h"
#include "asterfort/nm1dis.h"
#include "asterfort/rcvalb.h"
#include "asterfort/relax_acier_cable.h"
#include "asterfort/utmess.h"
#include "asterfort/verift.h"
!
    type(Material_Para), intent(in) :: materPara
    character(len=16), intent(in) :: option, relaComp
    real(kind=8) :: carcri(CARCRI_SIZE)
    real(kind=8) :: xlong0
    real(kind=8) :: aire
    real(kind=8) :: tmoins
    real(kind=8) :: tplus
    real(kind=8) :: dlong0
    real(kind=8) :: effnom
    real(kind=8) :: vim(*)
    real(kind=8) :: effnop
    real(kind=8) :: vip(*)
    real(kind=8) :: klv(21)
    real(kind=8) :: fono(6)
    real(kind=8) :: epsm
    integer(kind=8) :: codret
!
! --------------------------------------------------------------------------------------------------
!
!    TRAITEMENT DE LA RELATION DE COMPORTEMENT -ELASTOPLASTICITE-
!    ECROUISSAGE ISOTROPE ET CINEMATIQUE- LINEAIRE - VON MISES-
!    POUR UN MODELE BARRE ELEMENT MECA_BARRE
!
! --------------------------------------------------------------------------------------------------
!
!       XLONG0 : LONGUEUR DE L'ELEMENT DE BARRE AU REPOS
!       aire   : SECTION DE LA BARRE
!       TMOINS : INSTANT PRECEDENT
!       TPLUS  : INSTANT COURANT
!       DLONG0 : INCREMENT D'ALLONGEMENT DE L'ELEMENT
!       EFFNOM : EFFORT NORMAL PRECEDENT
!       TREF   : TEMPERATURE DE REFERENCE
!       TEMPM  : TEMPERATURE IMPOSEE A L'INSTANT PRECEDENT
!       TEMPP  : TEMPERATURE IMPOSEE A L'INSTANT COURANT
!       OPTION : OPTION DEMANDEE (R_M_T,FULL OU RAPH_MECA)
! OUT : EFFNOP : CONTRAINTE A L'INSTANT ACTUEL
!       VIP    : VARIABLE INTERNE A L'INSTANT ACTUEL
!       FONO   : FORCES NODALES COURANTES
!       KLV    : MATRICE TANGENTE
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: materPoin = " "
    integer(kind=8) :: propCode(1)
    real(kind=8)  :: sigm, deps, depsth, depsm, em, ep
    real(kind=8)  :: sigp, xrig, propVale(1), dsde
    aster_logical :: isot, cine, elas, corr, implex, isotli, relax
!
! --------------------------------------------------------------------------------------------------
!
    elas = ASTER_FALSE
    isot = ASTER_FALSE
    cine = ASTER_FALSE
    corr = ASTER_FALSE
    implex = option .eq. 'RIGI_MECA_IMPLEX' .or. option .eq. 'RAPH_MECA_IMPLEX'
    isotli = ASTER_FALSE
    relax = ASTER_FALSE
    if (relaComp .eq. 'ELAS') then
        elas = ASTER_TRUE
    else if ((relaComp .eq. 'VMIS_ISOT_LINE') .or. &
             (relaComp .eq. 'VMIS_ISOT_TRAC')) then
        isot = ASTER_TRUE
        if (relaComp .eq. 'VMIS_ISOT_LINE') then
            isotli = ASTER_TRUE
        end if
    else if (relaComp .eq. 'VMIS_CINE_LINE') then
        cine = ASTER_TRUE
    else if (relaComp .eq. 'CORR_ACIER') then
        corr = ASTER_TRUE
    else if (relaComp .eq. 'RELAX_ACIER') then
        relax = ASTER_TRUE
    end if
    if (implex) then
        if ((.not. elas) .and. (.not. isotli)) then
            call utmess('F', 'POUTRE0_49', sk=relaComp)
        end if
    end if
!
    klv = 0.d0
    fono = 0.d0
!
!   Récupération des caractéristiques
    deps = dlong0/xlong0
    sigm = effnom/aire
!
    if (isot .and. (.not. implex)) then
!       Caractéristiques élastiques a t-
        call rcvalb(materPara%schemePara%fami, &
                    materPara%schemePara%kpg, &
                    materPara%schemePara%ksp, &
                    '-', &
                    materPara%jvMaterCode, &
                    materPoin, 'ELAS', &
                    0, ' ', [0.d0], &
                    1, 'E', propVale, propCode, 1)
        em = propVale(1)

!       Caractéristiques élastiques a t+
        call rcvalb(materPara%schemePara%fami, &
                    materPara%schemePara%kpg, &
                    materPara%schemePara%ksp, &
                    '+', &
                    materPara%jvMaterCode, &
                    materPoin, 'ELAS', &
                    0, ' ', [0.d0], &
                    1, 'E', propVale, propCode, 1)
        ep = propVale(1)
!
        call verift(materPara%schemePara%fami, &
                    materPara%schemePara%kpg, &
                    materPara%schemePara%ksp, &
                    'T', &
                    materPara%jvMaterCode, &
                    epsth_=depsth)
        depsm = deps-depsth
        call nm1dis(materPara, &
                    option, relaComp, materPoin, &
                    em, ep, sigm, depsm, vim, &
                    sigp, vip, dsde)

    else if (cine) then
!       Caractéristiques élastiques a t-
        call rcvalb(materPara%schemePara%fami, &
                    materPara%schemePara%kpg, &
                    materPara%schemePara%ksp, &
                    '-', &
                    materPara%jvMaterCode, &
                    materPoin, 'ELAS', &
                    0, ' ', [0.d0], &
                    1, 'E', propVale, propCode, 1)
        em = propVale(1)
!       Caractéristiques élastiques a t+
        call rcvalb(materPara%schemePara%fami, &
                    materPara%schemePara%kpg, &
                    materPara%schemePara%ksp, &
                    '+', &
                    materPara%jvMaterCode, &
                    materPoin, 'ELAS', &
                    0, ' ', [0.d0], &
                    1, 'E', propVale, propCode, 1)
        ep = propVale(1)
!
        call verift(materPara%schemePara%fami, &
                    materPara%schemePara%kpg, &
                    materPara%schemePara%ksp, &
                    'T', &
                    materPara%jvMaterCode, &
                    epsth_=depsth)
        depsm = deps-depsth
        call nm1dci(materPara, &
                    option, materPoin, &
                    em, ep, sigm, depsm, vim, &
                    sigp, vip, dsde)

    else if (relax) then
        call verift(materPara%schemePara%fami, &
                    materPara%schemePara%kpg, &
                    materPara%schemePara%ksp, &
                    'T', &
                    materPara%jvMaterCode, &
                    epsth_=depsth)
        depsm = deps-depsth
        call relax_acier_cable(materPara%schemePara%fami, &
                               materPara%schemePara%kpg, &
                               materPara%schemePara%ksp, &
                               materPara%jvMaterCode, &
                               sigm, epsm, depsm, vim, &
                               sigp, vip, dsde)

    else if (elas) then
!       Caractéristiques élastiques a t-
        call rcvalb(materPara%schemePara%fami, &
                    materPara%schemePara%kpg, &
                    materPara%schemePara%ksp, &
                    '-', &
                    materPara%jvMaterCode, &
                    materPoin, 'ELAS', &
                    0, ' ', [0.d0], &
                    1, 'E', propVale, propCode, 1)
        em = propVale(1)

!       Caractéristiques élastiques a t+
        call rcvalb(materPara%schemePara%fami, &
                    materPara%schemePara%kpg, &
                    materPara%schemePara%ksp, &
                    '+', &
                    materPara%jvMaterCode, &
                    materPoin, 'ELAS', &
                    0, ' ', [0.d0], &
                    1, 'E', propVale, propCode, 1)
        ep = propVale(1)
!
        dsde = ep
        vip(1) = 0.d0
        call verift(materPara%schemePara%fami, &
                    materPara%schemePara%kpg, &
                    materPara%schemePara%ksp, &
                    'T', &
                    materPara%jvMaterCode, &
                    epsth_=depsth)
        sigp = ep*(sigm/em+deps-depsth)

    else if (corr) then
!       Caractéristiques élastiques a t-
        call rcvalb(materPara%schemePara%fami, &
                    materPara%schemePara%kpg, &
                    materPara%schemePara%ksp, &
                    '-', &
                    materPara%jvMaterCode, &
                    materPoin, 'ELAS', &
                    0, ' ', [0.d0], &
                    1, 'E', propVale, propCode, 1)
        em = propVale(1)
!       Caractéristiques élastiques a t+
        call rcvalb(materPara%schemePara%fami, &
                    materPara%schemePara%kpg, &
                    materPara%schemePara%ksp, &
                    '+', &
                    materPara%jvMaterCode, &
                    materPoin, 'ELAS', &
                    0, ' ', [0.d0], &
                    1, 'E', propVale, propCode, 1)
        ep = propVale(1)
!
        call nm1dco(materPara, option, carcri, &
                    materPoin, &
                    ep, sigm, epsm, deps, &
                    vim, sigp, vip, dsde, &
                    codret)

    else if (implex) then
!       Caractéristiques élastiques a t-
        call rcvalb(materPara%schemePara%fami, &
                    materPara%schemePara%kpg, &
                    materPara%schemePara%ksp, &
                    '-', &
                    materPara%jvMaterCode, &
                    materPoin, 'ELAS', &
                    0, ' ', [0.d0], &
                    1, 'E', propVale, propCode, 1)
        em = propVale(1)
!       Caractéristiques élastiques a t+
        call rcvalb(materPara%schemePara%fami, &
                    materPara%schemePara%kpg, &
                    materPara%schemePara%ksp, &
                    '+', &
                    materPara%jvMaterCode, &
                    materPoin, 'ELAS', &
                    0, ' ', [0.d0], &
                    1, 'E', propVale, propCode, 1)
        ep = propVale(1)
        call lcimpl(materPara%schemePara%fami, &
                    materPara%schemePara%kpg, &
                    materPara%schemePara%ksp, &
                    materPara%jvMaterCode, &
                    em, ep, sigm, tmoins, tplus, deps, &
                    vim, option, sigp, vip, dsde)

    else
        ASSERT(ASTER_FALSE)
    end if
!
!   Calcul du coefficient non nul de la matrice tangente
    if (option(1:10) .eq. 'RIGI_MECA_' .or. option(1:9) .eq. 'FULL_MECA') then
        xrig = dsde*aire/xlong0
        klv(1) = xrig
        klv(7) = -xrig
        klv(10) = xrig
    end if
!
!   Calcul des forces nodales
    if (option(1:9) .eq. 'RAPH_MECA' .or. option(1:9) .eq. 'FULL_MECA') then
        effnop = sigp*aire
        fono(1) = -effnop
        fono(4) = effnop
    end if
!
    if (implex) then
        effnop = sigp*aire
        fono(1) = -effnop
        fono(4) = effnop
    end if
!
end subroutine
