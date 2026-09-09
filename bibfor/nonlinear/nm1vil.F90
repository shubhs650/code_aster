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
! aslint: disable=W0413
!
subroutine nm1vil(materPara, &
                  relaComp, carcri, &
                  materPoin, &
                  instam, instap, tm, tp, &
                  deps, sigm, vim, &
                  defam, defap, sigp, vip, &
                  dsidep, iret, nbvalc)
!
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterc/r8t0.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/granac.h"
#include "asterfort/MaterialPara_type.h"
#include "asterfort/nmasse.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    type(Material_Para), intent(inout) :: materPara
    character(len=16), intent(in) :: relaComp
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    character(len=*), intent(in) :: materPoin
    integer(kind=8) :: iret, nbvalc
    real(kind=8) :: instam, instap
    real(kind=8) :: tm, tp
    real(kind=8) :: irram, irrap
    real(kind=8) :: deps
    real(kind=8) :: sigm, vim(nbvalc)
    real(kind=8) :: defam, defap
    real(kind=8) :: sigp, vip(nbvalc), dsidep, alpha
!
! --------------------------------------------------------------------------------------------------
!
!      VISCO_PLASTICITE FLUAGE SOUS IRRADIATION AVEC GRANDISSEMENT
!      VISC_IRRA_LOG OU GRAN_IRRA_LOG
!      LOI 1D PURE. MODIF JMP POUR ECRIRE SIMPLEMENT :
! DEPSVP=SIGMA+.EXP(-Q/T)*(A.OMEGA/(1+OMEGA*FLUENCE)+B*FLUENCE)*DFLUENCE
!
! --------------------------------------------------------------------------------------------------
!
! IN  INSTAM  : INSTANT DU CALCUL PRECEDENT
! IN  INSTAP  : INSTANT DU CALCUL
! IN  TM      : TEMPERATURE A L'INSTANT PRECEDENT
! IN  TP      : TEMPERATURE A L'INSTANT DU CALCUL
! IN  TREF    : TEMPERATURE DE REFERENCE
! IN  DEPS    : INCREMENT DE DEFORMATION-INCREMENT DEFORMATION THERMIQUE
! IN  SIGM    : CONTRAINTES A L'INSTANT DU CALCUL PRECEDENT
! IN  VIM     : VARIABLES INTERNES A L'INSTANT DU CALCUL PRECEDENT
! IN  OPTION  : OPTION DEMANDEE : RIGI_MECA_TANG , FULL_MECA , RAPH_MECA
! IN  DEFAM   : DEFORMATIONS ANELASTIQUES A L'INSTANT PRECEDENT
! IN  DEFAP   : DEFORMATIONS ANELASTIQUES A L'INSTANT DU CALCUL
! IN  ANGMAS  : LES TROIS ANGLES DU MOT_CLEF MASSIF (AFFE_CARA_ELEM)
! OUT SIGP    : CONTRAINTES A L'INSTANT ACTUEL
! OUT VIP     : VARIABLES INTERNES A L'INSTANT ACTUEL
! OUT DSIDEP  : MODULE TANGENT
! OUT IRET    : CODE RETOUR DE LA RECHERCHE DE ZERO DE F(X)=0
!                   IRET=0 => PAS DE PROBLEME
!                   IRET=1 => ECHEC
!
! --------------------------------------------------------------------------------------------------
!
!     COMMON POUR LES PARAMETRES DES LOIS VISCOPLASTIQUES
    common/nmpavp/dpc, sieleq, deuxmu, deltat, tschem, prec, theta, niter
    real(kind=8) :: dpc, sieleq, deuxmu, deltat, tschem, prec, theta, niter
!     COMMON POUR LES PARAMETRES DES LOIS DE FLUAGE SOUS IRRADIATION
!     VISC_IRRA_LOG: A      B      CTPS    ENER
    common/nmpair/a, b, ctps, ener
    real(kind=8) :: a, b, c, ctps, ener
! PARAMETRES MATERIAUX
! ELASTIQUES
    real(kind=8) :: ep, nup, troikp, deumup
    real(kind=8) :: em, num, troikm, deumum
! AUTRES
    integer(kind=8) :: iret2
    integer(kind=8), parameter  :: nbProp = 5
    real(kind=8) :: propVale(nbProp)
    character(len=8), parameter :: propName(nbProp) = &
                                   (/'A       ', 'B       ', &
                                     'CSTE_TPS', 'ENER_ACT', &
                                     'C       '/)
    integer(kind=8) :: propCode(nbProp)
    real(kind=8) :: t1, t2
    real(kind=8) :: degran, depsan, depsim, depsgr
    real(kind=8) :: coef1, coefb, expqt
    real(kind=8) :: fluphi
!
! --------------------------------------------------------------------------------------------------
!
    iret = 0
!     PARAMETRE THETA D'INTEGRATION
!
    theta = carcri(4)
    t1 = abs(theta-0.5d0)
    t2 = abs(theta-1.d0)
    prec = 0.000001d0
    if ((t1 .gt. prec) .and. (t2 .gt. prec)) then
        call utmess('F', 'ALGORITH6_55')
    end if
!
! TEMPERATURE AU MILIEU DU PAS DE TEMPS  (DANS COMMON / NMPAVP /)
    tschem = tm*(1.d0-theta)+tp*theta
! DEFORMATION PLASTIQUE CUMULEE  (DANS COMMON / NMPAVP /)
    dpc = vim(1)
! INCREMENT DE TEMPS (DANS COMMON / NMPAVP /)
    deltat = instap-instam
! CARACTERISTIQUES ELASTIQUES VARIABLES

    call nmasse(materPara%schemePara%fami, &
                materPara%schemePara%kpg, &
                materPara%schemePara%ksp, &
                '-', &
                materPara%jvMaterCode, &
                materPoin, instam, em, num, deumum, &
                troikm)
!
    call nmasse(materPara%schemePara%fami, &
                materPara%schemePara%kpg, &
                materPara%schemePara%ksp, &
                '+', &
                materPara%jvMaterCode, &
                materPoin, instap, ep, nup, deumup, &
                troikp)
!
!     IRRADIATION AU POINT CONSIDERE
!     FLUX NEUTRONIQUE
    call rcvarc('F', 'IRRA', '-', &
                materPara%schemePara%fami, &
                materPara%schemePara%kpg, &
                materPara%schemePara%ksp, &
                irram, iret2)

    if (iret2 .gt. 0) irram = 0.d0
    call rcvarc('F', 'IRRA', '+', &
                materPara%schemePara%fami, &
                materPara%schemePara%kpg, &
                materPara%schemePara%ksp, &
                irrap, iret2)
    if (iret2 .gt. 0) irrap = 0.d0
    irrap = irrap-irram+vim(2)
    irram = vim(2)
!
    fluphi = (irrap-irram)/deltat

! - RECUPERATION DES CARACTERISTIQUES DES LOIS DE FLUAGE
    call rcvalb(materPara%schemePara%fami, &
                materPara%schemePara%kpg, &
                materPara%schemePara%ksp, &
                '+', &
                materPara%jvMaterCode, &
                materPoin, relaComp, &
                0, ' ', [0.d0], &
                nbProp, propName, propVale, &
                propCode, 0)

!     TRAITEMENT DES PARAMETRES DE LA LOI DE FLUAGE
    if (propCode(1) .eq. 0) then
        a = propVale(1)
        b = propVale(2)
        ctps = propVale(3)
        ener = propVale(4)
        if (relaComp(1:10) .eq. 'GRAN_IRRA_') then
            c = propVale(5)
        else
            c = 0.0d0
        end if
        if (fluphi .lt. -prec) then
            call utmess('F', 'ALGORITH6_57')
        end if
    else
        call utmess('F', 'ALGORITH6_58')
    end if
!
!     CALCUL DE LA DEFORMATION DE GRANDISSEMENT
    degran = 0.0d0
    call granac(materPara%schemePara%fami, &
                materPara%schemePara%kpg, &
                materPara%schemePara%ksp, &
                materPara%jvMaterCode, &
                materPoin, &
                relaComp, irrap, irram, tm, tp, &
                depsgr)
!
    if (relaComp(1:10) .eq. 'GRAN_IRRA_') then
        vip(3) = vim(3)+depsgr
        if (depsgr .ne. 0.d0) then
            ASSERT(materPara%lcsPara%lcsType .eq. MATER_LCS_ZERO)
            alpha = 0.d0

!        INCREMENT DEFORMATION DE GRANDISSEMENT DANS LE REPERE
            degran = depsgr*cos(alpha)*cos(alpha)
        end if
    end if
!     INCREMENT DEFORMATION ANELASTIQUE
    depsan = defap-defam
!     INCREMENT DEFORMATION IMPOSEE
    depsim = depsan+degran
!
    expqt = exp(-ener/(tp+r8t0()))
!
!     coefb=expqt*((a*ctps/(1.d0+ctps*irrap))+b+c*ctps*exp(-ctps*irrap))*(irrap-irram)
    coefb = expqt*(a*(log(1.d0+ctps*irrap)-log(1.d0+ctps*irram))+ &
                   b*(irrap-irram)+ &
                   c*(exp(-ctps*irram)-exp(-ctps*irrap)))
    coef1 = ep/(1.d0+ep*coefb)
!
! CONTRAINTE ACTUALISEE
!
!
    sigp = coef1*(sigm/em+deps-depsim)
!
! DEFORMATION PLASTIQUE CUMULEE ACTUALISEE
!
    vip(1) = vim(1)+(abs(sigp)*coefb)
! MODULE TANGENT POUR MATRICE TANGENTE
!
    dsidep = coef1
    vip(2) = irrap
!
end subroutine
