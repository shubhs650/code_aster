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
subroutine nmco1d(BEHInteg, &
                  relaComp, &
                  option, deps, sigm, &
                  vim, sigp, vip, dsidep, codret)
!
    use Behaviour_type
    use Behaviour_module
    implicit none
!
#include "asterf_types.h"
#include "asterfort/nm1dci.h"
#include "asterfort/nm1dis.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
#include "asterfort/verift.h"
#include "asterfort/vmci1d.h"
#include "jeveux.h"
!
    type(Behaviour_Integ), intent(inout) :: BEHInteg
    integer(kind=8) :: codret
    character(len=16) :: option, relaComp
    real(kind=8) :: deps, sigm, vim(*)
    real(kind=8) :: sigp, vip(*), dsidep
!
! --------------------------------------------------------------------------------------------------
!
!          REALISE LES LOIS 1D (DEBORST OU EXPLICITEMENT 1D)
!
! --------------------------------------------------------------------------------------------------
!
! IN  IMATE   : ADRESSE DU MATERIAU CODE
! IN  OPTION  : OPTION DEMANDEE : RIGI_MECA_TANG , FULL_MECA , RAPH_MECA
! IN  EPSM    : DEFORMATIONS A L'INSTANT DU CALCUL PRECEDENT
! IN  DEPS    : INCREMENT DE DEFORMATION (SCALAIRE DANS CE CAS)
! IN  SIGM    : CONTRAINTE A L'INSTANT DU CALCUL PRECEDENT
! IN  VIM     : VARIABLES INTERNES A L'INSTANT DU CALCUL PRECEDENT
! IN   TM     : TEMPERATURE L'INSTANT DU CALCUL PRECEDENT
! IN   TP     : TEMPERATURE A L'INSTANT DU
! IN  TREF    : TEMPERATURE DE REFERENCE
! OUT SIGP    : CONTRAINTE A L'INSTANT ACTUEL
!     VIP     : VARIABLES INTERNES A L'INSTANT ACTUEL
!     DSIDEP  : RIGIDITE (SCALAIRE DANS CE CAS)
!     CODRET  : CODE RETOUR NON NUL SI SIGYY OU SIGZZ NON NULS
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbProp = 1
    character(len=16), parameter :: propName(nbProp) = (/'E'/)
    real(kind=8) :: propVale(nbProp)
    integer(kind=8) :: propCode(nbProp)
    character(len=8), parameter :: materPoin = ' '
    aster_logical :: cine, isot, elas, cinegc
    real(kind=8) :: em, ep, depsth, depsm
! --------------------------------------------------------------------------------------------------
!
    elas = ASTER_FALSE
    isot = ASTER_FALSE
    cine = ASTER_FALSE
    cinegc = ASTER_FALSE
    codret = 0
    sigp = 0.d0
!
    if (relaComp(1:16) .eq. 'GRILLE_ISOT_LINE') then
        isot = ASTER_TRUE
    else if (relaComp(1:16) .eq. 'GRILLE_CINE_LINE') then
        cine = ASTER_TRUE
    else if (relaComp(1:12) .eq. 'VMIS_CINE_GC') then
        cinegc = ASTER_TRUE
    else if (relaComp(1:4) .eq. 'ELAS') then
        elas = ASTER_TRUE
    else
        call utmess('F', 'COMPOR4_32', sk=relaComp)
    end if
!
    call rcvalb(BEHInteg%materPara%schemePara%fami, &
                BEHInteg%materPara%schemePara%kpg, &
                BEHInteg%materPara%schemePara%ksp, &
                '-', &
                BEHInteg%materPara%jvMaterCode, &
                ' ', 'ELAS', &
                0, ' ', [0.d0], &
                nbProp, propName, propVale, &
                propCode, 1)
    em = propVale(1)
    call rcvalb(BEHInteg%materPara%schemePara%fami, &
                BEHInteg%materPara%schemePara%kpg, &
                BEHInteg%materPara%schemePara%ksp, &
                '+', &
                BEHInteg%materPara%jvMaterCode, &
                ' ', 'ELAS', &
                0, ' ', [0.d0], &
                nbProp, propName, propVale, &
                propCode, 1)
    ep = propVale(1)
!
    if (isot) then
        call verift(BEHInteg%materPara%schemePara%fami, &
                    BEHInteg%materPara%schemePara%kpg, &
                    BEHInteg%materPara%schemePara%ksp, &
                    'T', BEHInteg%materPara%jvMaterCode, epsth_=depsth)
        depsm = deps-depsth
        call nm1dis(BEHINteg%materPara, &
                    option, relaComp, materPoin, &
                    em, ep, &
                    sigm, depsm, vim, &
                    sigp, vip, dsidep)

    else if (cine) then
        call verift(BEHInteg%materPara%schemePara%fami, &
                    BEHInteg%materPara%schemePara%kpg, &
                    BEHInteg%materPara%schemePara%ksp, &
                    'T', BEHInteg%materPara%jvMaterCode, epsth_=depsth)
        depsm = deps-depsth
        call nm1dci(BEHInteg%materPara, &
                    option, materPoin, &
                    em, ep, &
                    sigm, depsm, vim, &
                    sigp, vip, dsidep)

    else if (cinegc) then
        call verift(BEHInteg%materPara%schemePara%fami, &
                    BEHInteg%materPara%schemePara%kpg, &
                    BEHInteg%materPara%schemePara%ksp, &
                    'T', BEHInteg%materPara%jvMaterCode, epsth_=depsth)
        depsm = deps-depsth
        call vmci1d(BEHINteg%materPara, &
                    option, materPoin, &
                    em, ep, &
                    sigm, depsm, vim, &
                    sigp, vip, dsidep)

    else if (elas) then
        if (option(1:9) .eq. 'FULL_MECA' .or. option(1:10) .eq. 'RIGI_MECA_') then
            dsidep = ep
        end if
        if (option(1:9) .eq. 'RAPH_MECA' .or. option(1:9) .eq. 'FULL_MECA') then
            vip(1) = 0.d0
            call verift(BEHInteg%materPara%schemePara%fami, &
                        BEHInteg%materPara%schemePara%kpg, &
                        BEHInteg%materPara%schemePara%ksp, &
                        'T', BEHInteg%materPara%jvMaterCode, &
                        epsth_=depsth)
            sigp = ep*(sigm/em+deps-depsth)
        end if

    end if
end subroutine
