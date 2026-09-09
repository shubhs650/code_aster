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
subroutine nmiclg(materPara, &
                  option, relaComp, carcri, &
                  epsm, deps, sigm, vim, &
                  sigp, vip, dsde, &
                  codret)
!
    use MaterialPara_type
    use MaterialPara_module
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/lcimpl.h"
#include "asterfort/nm1das.h"
#include "asterfort/nm1dci.h"
#include "asterfort/nm1dco.h"
#include "asterfort/nm1dis.h"
#include "asterfort/nmmaba.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
#include "asterfort/verift.h"
#include "asterfort/Behaviour_type.h"
!
    type(Material_Para), intent(in) :: materPara
    character(len=16), intent(in) :: option, relaComp
    real(kind=8) :: carcri(CARCRI_SIZE)
    real(kind=8) :: vim(*)
    real(kind=8) :: vip(*)
    real(kind=8) :: sigy, sigm, deps, sigp
    real(kind=8) :: dsde, epsm
    integer(kind=8) :: codret
!
! --------------------------------------------------------------------------------------------------
!
!    TRAITEMENT DE LA RELATION DE COMPORTEMENT -ELASTOPLASTICITE-
!    ECROUISSAGE ISOTROPE ET CINEMATIQUE- LINEAIRE - VON MISES-
!    POUR UN MODELE CABLE_GAINE (EQUIVALENT A MECA_BARRE)
!
! --------------------------------------------------------------------------------------------------
!
!       OPTION : OPTION DEMANDEE (R_M_T,FULL OU RAPH_MECA)
!       IMATE : POINTEUR MATERIAU CODE
!       EPSM  : DEFORMATION A L'INSTANT MOINS
!       DEPS  : INCREMENT DE DEFORMATION
!       SIGM  : CONTRAINTE A L'INSTANT MOINS
!       VIM   : VARIABLE INTERNE A L'INSTANT MOINS
! OUT : SIGP  : CONTRAINTE A L'INSTANT ACTUEL
!       VIP    : VARIABLE INTERNE A L'INSTANT ACTUEL
!       DSDE   : MATRICE TANGENTE
!       CRILDC : CRITERE LOI DE COMPORTEMENT
!       CODRET : CODE RETOUR
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: materPoin = " "
    integer(kind=8), parameter :: kpgFPG1 = 1, kspFPG1 = 1
    character(len=8), parameter :: famiFPG1 = "FPG1"
    type(Material_Para) :: materParaFPG1
    integer(kind=8), parameter :: nbProp = 4
    integer(kind=8) :: propCode(nbProp)
    real(kind=8) :: propVale(nbProp)
    character(len=16), parameter :: propName(nbProp) = &
                                    (/'SY_C        ', 'DC_SIGM_EPSI', &
                                      'SY_T        ', 'DT_SIGM_EPSI'/)
    real(kind=8) :: depsth, depsm, tmoins, tplus
    real(kind=8) :: em, ep, dsdem, dsdep
    real(kind=8) :: syc, etc, syt, ett
    aster_logical :: isot, cine, elas, corr, implex, isotli, asyml, sans
!
! --------------------------------------------------------------------------------------------------
!
    elas = ASTER_FALSE
    isot = ASTER_FALSE
    cine = ASTER_FALSE
    corr = ASTER_FALSE
    implex = option .eq. 'RIGI_MECA_IMPLEX'
    isotli = ASTER_FALSE
    asyml = ASTER_FALSE
    if (relaComp .eq. 'ELAS') then
        elas = ASTER_TRUE
    else if ((relaComp .eq. 'VMIS_ISOT_LINE') .or. (relaComp .eq. 'VMIS_ISOT_TRAC')) then
        isot = ASTER_TRUE
        if (relaComp .eq. 'VMIS_ISOT_LINE') then
            isotli = ASTER_TRUE
        end if
    else if (relaComp .eq. 'VMIS_CINE_LINE') then
        cine = ASTER_TRUE
    else if (relaComp .eq. 'CORR_ACIER') then
        corr = ASTER_TRUE
    else if (relaComp .eq. 'VMIS_ASYM_LINE') then
        asyml = ASTER_TRUE
    else if (relaComp .eq. 'SANS') then
        sans = ASTER_TRUE
    end if
    if (implex) then
        if ((.not. elas) .and. (.not. isotli)) then
            call utmess('F', 'POUTRE0_49', sk=relaComp)
        end if
    end if

! - CARACTERISTIQUES ELASTIQUES A TMOINS
    call rcvalb(materPara%schemePara%fami, &
                materPara%schemePara%kpg, &
                materPara%schemePara%ksp, &
                '-', &
                materPara%jvMaterCode, &
                materPoin, 'ELAS', &
                0, ' ', [0.d0], &
                1, 'E', propVale, propCode, 1)
    em = propVale(1)

! - CARACTERISTIQUES ELASTIQUES A TPLUS
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
    if (isot .and. (.not. implex)) then
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
        call verift(materPara%schemePara%fami, &
                    materPara%schemePara%kpg, &
                    materPara%schemePara%ksp, &
                    'T', &
                    materPara%jvMaterCode, &
                    epsth_=depsth)
        depsm = deps-depsth
        call nm1dci(materPara, &
                    option, ' ', &
                    em, ep, sigm, depsm, vim, &
                    sigp, vip, dsde)

    else if (elas) then
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
        call nm1dco(materPara, option, carcri, &
                    ' ', &
                    ep, sigm, epsm, deps, &
                    vim, sigp, vip, dsde, &
                    codret)

    else if (implex) then
        call lcimpl(materPara%schemePara%fami, &
                    materPara%schemePara%kpg, &
                    materPara%schemePara%ksp, &
                    materPara%jvMaterCode, &
                    em, ep, sigm, tmoins, tplus, deps, &
                    vim, option, sigp, vip, dsde)

    else if (asyml) then
        call nmmaba(materPara%jvMaterCode, relaComp, &
                    ep, dsde, sigy)

! ----- Copy material parameters with other scheme parameters
        call copyMaterPara(materPara, famiFPG1, kpgFPG1, kspFPG1, &
                           materParaFPG1)

! ----- CARACTERISTIQUES ECROUISSAGE LINEAIRE ASYMETRIQUE
        call rcvalb(materParaFPG1%schemePara%fami, &
                    materParaFPG1%schemePara%kpg, &
                    materParaFPG1%schemePara%ksp, &
                    '+', &
                    materParaFPG1%jvMaterCode, &
                    ' ', 'ECRO_ASYM_LINE', &
                    0, ' ', [0.d0], &
                    nbProp, propName, propVale, &
                    propCode, 1)
        syc = propVale(1)
        etc = propVale(2)
        syt = propVale(3)
        ett = propVale(4)
        call nm1das(materPara, &
                    ep, syc, &
                    syt, etc, ett, &
                    sigm, deps, vim, &
                    sigp, vip, dsdem, dsdep)

        if (option(1:10) .eq. 'RIGI_MECA_' .or. option(1:9) .eq. 'FULL_MECA') then
            if (option(11:14) .eq. 'ELAS') then
                dsde = ep
            else
                if (option(1:14) .eq. 'RIGI_MECA_TANG') then
                    dsde = dsdem
                else
                    dsde = dsdep
                end if
            end if
        end if

    else if (sans) then
        sigp = 0.d0
        dsde = 0.d0

    else
        ASSERT(ASTER_FALSE)

    end if
!
end subroutine
