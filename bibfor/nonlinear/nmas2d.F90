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
! aslint: disable=W1504
!
subroutine nmas2d(BEHInteg, &
                  nno, npg, &
                  ipoids, ivf, idfde, &
                  geom, typmod, option, &
                  compor, multComp, lgpg, carcri, instam, instap, &
                  deplm, deplp, sigm, vim, &
                  dfdi, def, sigp, vip, matuu, &
                  vectu, codret)
!
    use Behaviour_type
    use Behaviour_module
    implicit none
!
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/calsta.h"
#include "asterfort/codere.h"
#include "asterfort/dfda2d.h"
#include "asterfort/iniqs4.h"
#include "asterfort/nmcomp.h"
#include "asterfort/nmgeom.h"
#include "jeveux.h"
!
    type(Behaviour_Integ), intent(inout) :: BEHInteg
    character(len=8), intent(in) :: typmod(2)
    character(len=16), intent(in) :: compor(COMPOR_SIZE), multComp
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    integer(kind=8) :: nno, npg, lgpg, codret
    integer(kind=8) :: ipoids, ivf, idfde
    character(len=16) :: option
    real(kind=8) :: instam, instap
    real(kind=8) :: geom(2, nno)
    real(kind=8) :: deplm(2, nno), deplp(2, nno), dfdi(nno, 2)
    real(kind=8) :: def(4, nno, 2)
    real(kind=8) :: sigm(10, npg), sigp(10, npg)
    real(kind=8) :: vim(lgpg, npg), vip(lgpg, npg)
    real(kind=8) :: matuu(*), vectu(2, nno)
!
! --------------------------------------------------------------------------------------------------
!
!     BUT:  CALCUL  DES OPTIONS RIGI_MECA_TANG, RAPH_MECA ET FULL_MECA
!           EN HYPO-ELASTICITE EN 2D POUR LE QUAD4 SOUS INTEGRE
!           STABILITE PAR ASSUMED STRAIN
!
! --------------------------------------------------------------------------------------------------
!
! IN  NNO     : NOMBRE DE NOEUDS DE L'ELEMENT
! IN  NPG     : NOMBRE DE POINTS DE GAUSS
! IN  POIDSG  : POIDS DES POINTS DE GAUSS
! IN  VFF     : VALEUR  DES FONCTIONS DE FORME
! IN  DFDE    : DERIVEE DES FONCTIONS DE FORME ELEMENT DE REFERENCE
! IN  DFDK    : DERIVEE DES FONCTIONS DE FORME ELEMENT DE REFERENCE
! IN  GEOM    : COORDONEES DES NOEUDS
! IN  TYPMOD  : TYPE DE MODELISATION
! IN  OPTION  : OPTION DE CALCUL
! IN  COMPOR  : COMPORTEMENT
! IN  LGPG    : "LONGUEUR" DES VARIABLES INTERNES POUR 1 POINT DE GAUSS
!               CETTE LONGUEUR EST UN MAJORANT DU NBRE REEL DE VAR. INT.
! IN  INSTAM  : INSTANT PRECEDENT
! IN  INSTAP  : INSTANT DE CALCUL
! IN  DEPLM   : DEPLACEMENT A L'INSTANT PRECEDENT
! IN  DEPLP   : INCREMENT DE DEPLACEMENT
! IN  SIGM    : CONTRAINTES A L'INSTANT PRECEDENT
! IN  VIM     : VARIABLES INTERNES A L'INSTANT PRECEDENT
! OUT DFDI    : DERIVEE DES FONCTIONS DE FORME  AU DERNIER PT DE GAUSS
! OUT DEF     : PRODUIT DER. FCT. FORME PAR F   AU DERNIER PT DE GAUSS
! OUT SIGP    : CONTRAINTES DE CAUCHY (RAPH_MECA ET FULL_MECA)
! OUT VIP     : VARIABLES INTERNES    (RAPH_MECA ET FULL_MECA)
! OUT MATUU   : MATRICE DE RIGIDITE PROFIL (RIGI_MECA_TANG ET FULL_MECA)
! OUT VECTU   : FORCES NODALES (RAPH_MECA ET FULL_MECA)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: ksp = 1, kpg = 1
    aster_logical, parameter :: grand = ASTER_FALSE
    real(kind=8), parameter :: rac2 = sqrt(2.d0)
    integer(kind=8), parameter :: ndimLdc = 2
    aster_logical :: axi
    integer(kind=8) :: kk, kkd, n, i, m, j, j1, kl, kpgs
    integer(kind=8) :: cod(9), npgs
    real(kind=8) :: dsidep(6, 6), f(3, 3), eps(6), deps(6), r, sigma(6), sign(6)
    real(kind=8) :: poids, tmp, sig(6)
    real(kind=8) :: jac, sigas(4, 4), pqx, pqy, defc(4, 4, 2)
    real(kind=8) :: dh(8), gamma(8), coopg(8)
    real(kind=8) :: sdkdx(4), sdkdy(4), sdedx(4), sdedy(4), poi2sg(4)
    real(kind=8) :: sdfdy(4, 4), sdfdx(4, 4), sdfde(4, 4), sdfdk(4, 4)
    real(kind=8) :: qplus(6), qmoins(6), dq(6), defn(4, 4, 2), kron(3, 3)
    character(len=16) :: optios
    data kron/1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0/
!
! --------------------------------------------------------------------------------------------------
!
    axi = typmod(1) .eq. 'AXIS'
    cod = 0
!
    do i = 1, 3
        do j = 1, 3
            f(i, j) = kron(i, j)
        end do
    end do

! - Prepare external state variables (geometry)
    call behaviourPrepESVAGeom(nno, npg, ndimLdc, &
                               ipoids, ivf, idfde, &
                               geom, BEHInteg, &
                               deplm, deplp)

! - INITIALISATION QUAS4
    call iniqs4(nno, sdfde, sdfdk, poi2sg, coopg)
!
! - CALCUL DU VECTEUR GAMMA
    gamma(1) = (geom(1, 4)*(geom(2, 2)-geom(2, 3))+ &
                geom(1, 2)*(geom(2, 3)-geom(2, 4))+ &
                geom(1, 3)*(geom(2, 4)-geom(2, 2)))/ &
               (2*(((geom(1, 4)-geom(1, 2))*(geom(2, 1)-geom(2, 3)))+ &
                   (geom(1, 1)-geom(1, 3))*(geom(2, 2)-geom(2, 4))))

    gamma(2) = (geom(1, 4)*(geom(2, 3)-geom(2, 1))+ &
                geom(1, 3)*(geom(2, 1)-geom(2, 4))+ &
                geom(1, 1)*(geom(2, 4)-geom(2, 3)))/ &
               (2*(((geom(1, 4)-geom(1, 2))*(geom(2, 1)-geom(2, 3)))+ &
                   (geom(1, 1)-geom(1, 3))*(geom(2, 2)-geom(2, 4))))

    gamma(3) = (geom(1, 4)*(geom(2, 1)-geom(2, 2))+ &
                geom(1, 1)*(geom(2, 2)-geom(2, 4))+ &
                geom(1, 2)*(geom(2, 4)-geom(2, 1)))/ &
               (2*(((geom(1, 4)-geom(1, 2))*(geom(2, 1)-geom(2, 3)))+ &
                   (geom(1, 1)-geom(1, 3))*(geom(2, 2)-geom(2, 4))))

    gamma(4) = (geom(1, 3)*(geom(2, 1)-geom(2, 2))+ &
                geom(1, 1)*(geom(2, 2)-geom(2, 3))+ &
                geom(1, 2)*(geom(2, 3)-geom(2, 1)))/ &
               (2*(((geom(1, 2)-geom(1, 4))*(geom(2, 1)-geom(2, 3)))- &
                   (geom(1, 1)-geom(1, 3))*(geom(2, 2)-geom(2, 4))))

!   CALCUL DE DFDI,F,EPS,DEPS,R(EN AXI) ET POIDS
    eps = 0.d0
    call nmgeom(2, nno, axi, grand, geom, &
                kpg, ipoids, ivf, idfde, deplm, &
                .true._1, poids, dfdi, f, eps, &
                r)

!   CALCUL DE DEPS
    deps = 0.d0
    call nmgeom(2, nno, axi, grand, geom, &
                kpg, ipoids, ivf, idfde, deplp, &
                .true._1, poids, dfdi, f, deps, &
                r)

!   CALCUL DES PRODUITS SYMETR. DE F PAR N,
    do n = 1, nno
        do i = 1, 2
            def(1, n, i) = f(i, 1)*dfdi(n, 1)
            def(2, n, i) = f(i, 2)*dfdi(n, 2)
            def(3, n, i) = 0.d0
            def(4, n, i) = (f(i, 1)*dfdi(n, 2)+f(i, 2)*dfdi(n, 1))/rac2
        end do
    end do
!
    do i = 1, 3
        sign(i) = sigm(i, kpg)
    end do
    sign(4) = sigm(4, kpg)*rac2

! - LOI DE COMPORTEMENT
    if (option(1:9) .eq. 'RAPH_MECA') then
        optios = 'FULL_MECA'
    else
        optios = option
    end if

! - Set main parameters for behaviour (on point) - (only one gauss point !)
    call behaviourSetParaPoin(kpg, ksp, BEHInteg)

! - Integrator
    sigma = 0.d0
    call nmcomp(BEHInteg, &
                ndimLdc, optios, typmod, &
                instam, instap, &
                compor, carcri, multComp, &
                6, eps, deps, &
                6, sign, &
                vim(1, kpg), &
                sigma, vip(1, kpg), &
                36, dsidep, cod(kpg))

! - ERREUR D'INTEGRATION
    if (cod(kpg) .eq. 1) then
        goto 999
    end if

    if (option(1:10) .eq. 'RIGI_MECA_' .or. option(1:9) .eq. 'FULL_MECA') then
!       CALCUL DE KC (MATRICE DE RIGIDITE AU CENTRE)
        do n = 1, nno
            do i = 1, 2
                do kl = 1, 4
                    sig(kl) = 0.d0
                    sig(kl) = sig(kl)+def(1, n, i)*dsidep(1, kl)
                    sig(kl) = sig(kl)+def(2, n, i)*dsidep(2, kl)
                    sig(kl) = sig(kl)+def(3, n, i)*dsidep(3, kl)
                    sig(kl) = sig(kl)+def(4, n, i)*dsidep(4, kl)
                end do
                do j = 1, 2
                    do m = 1, n
                        if (m .eq. n) then
                            j1 = i
                        else
                            j1 = 2
                        end if
!               RIGIDITE ELASTIQUE
                        tmp = 0.d0
                        tmp = tmp+sig(1)*def(1, m, j)
                        tmp = tmp+sig(2)*def(2, m, j)
                        tmp = tmp+sig(3)*def(3, m, j)
                        tmp = tmp+sig(4)*def(4, m, j)
!               STOCKAGE EN TENANT COMPTE DE LA SYMETRIE
                        if (j .le. j1) then
                            kkd = (2*(n-1)+i-1)*(2*(n-1)+i)/2
                            kk = kkd+2*(m-1)+j
                            matuu(kk) = matuu(kk)+tmp*poids
                        end if
                    end do
                end do
            end do
        end do
!
!           CORRECTION DE LA MATRICE DE RIGIDITE
!                 CALCUL DE KSTAB
!     --------------------------------------------
        npgs = 4

!        CALCUL DES TERMES EVALUES AUX 4 POINTS DE GAUSS
        do kpgs = 1, npgs
            call dfda2d(kpgs, nno, poi2sg(kpgs), sdfde, sdfdk, &
                        sdedx, sdedy, sdkdx, sdkdy, sdfdx, &
                        sdfdy, geom, jac)
            dh(2*kpgs-1) = coopg(2*kpgs-1)*sdkdx(kpgs)+coopg(2*kpgs)*sdedx(kpgs)
            dh(2*kpgs) = coopg(2*kpgs-1)*sdkdy(kpgs)+coopg(2*kpgs)*sdedy(kpgs)
            call calsta(gamma, dh, def, nno, &
                        kpgs, sig, tmp, kk, kkd, &
                        matuu, dsidep, jac)
        end do
    end if

! - CALCUL DES FORCES INTERNES ET DES CONTRAINTES DE CAUCHY
    if (option(1:9) .eq. 'FULL_MECA' .or. option(1:9) .eq. 'RAPH_MECA') then
!
!     INITIALISATION
        npgs = 4
        pqx = 0.d0
        pqy = 0.d0

!       DEPLACEMENTS GENERALISES
        do kl = 1, nno
            pqx = pqx+gamma(kl)*deplp(1, kl)
            pqy = pqy+gamma(kl)*deplp(2, kl)
        end do

!      INCREMENT DES CONTRAINTES GENERALISEES
        do i = 1, 6
            qmoins(i) = sigm(i+4, kpg)
            dq(1) = (dsidep(1, 1)-dsidep(2, 1))*pqx
            dq(2) = (dsidep(2, 1)-dsidep(1, 1))*pqy
            dq(3) = (dsidep(1, 2)-dsidep(2, 2))*pqx
            dq(4) = (dsidep(2, 2)-dsidep(1, 2))*pqy
            dq(5) = 0.d0
            dq(6) = 0.d0
            qplus(i) = qmoins(i)+dq(i)
        end do
!
!
!      OPERATEUR DE GRADIENT AU CENTRE
        do n = 1, nno
            do i = 1, 2
                defc(1, n, i) = def(1, n, i)
                defc(2, n, i) = def(2, n, i)
                defc(3, n, i) = def(3, n, i)
                defc(4, n, i) = def(4, n, i)
            end do
        end do
!
!
!      OPERATEUR DE STABILISATION DU GRADIENT AU 4 POINTS DE GAUSS
        do kpgs = 1, npgs
            call dfda2d(kpgs, nno, poi2sg(kpgs), sdfde, sdfdk, &
                        sdedx, sdedy, sdkdx, sdkdy, sdfdx, &
                        sdfdy, geom, jac)
!
            dh(2*kpgs-1) = coopg(2*kpgs-1)*sdkdx(kpgs)+coopg(2*kpgs)*sdedx(kpgs)
            dh(2*kpgs) = coopg(2*kpgs-1)*sdkdy(kpgs)+coopg(2*kpgs)*sdedy(kpgs)
!
!
            do n = 1, nno
                do i = 1, 2
                    defn(1, n, i) = f(i, 1)*gamma(n)*dh(2*kpgs-1)*(0.5d0)+ &
                                    f(i, 2)*gamma(n)*dh(2*kpgs)*(-0.5d0)
                    defn(2, n, i) = f(i, 2)*gamma(n)*dh(2*kpgs)*0.5d0+ &
                                    f(i, 1)*gamma(n)*dh(2*kpgs-1)*(-0.5d0)
                    defn(3, n, i) = 0.d0
                    defn(4, n, i) = 0.d0
                end do
            end do

!           CONTRAINTES DE HOURGLASS
            sigas(1, kpgs) = (qplus(1)*dh(2*kpgs-1)+qplus(2)*dh(2*kpgs))
            sigas(2, kpgs) = (qplus(3)*dh(2*kpgs-1)+qplus(4)*dh(2*kpgs))
            sigas(3, kpgs) = 0.d0
            sigas(4, kpgs) = 0.d0

!           CALCUL DES FORCES INTERNES
            do n = 1, nno
                do i = 1, 2
                    do kl = 1, 3
                        vectu(i, n) = vectu(i, n)+defc(kl, n, i)*sigas(kl, kpgs)*jac+ &
                                      defn(kl, n, i)*sigas(kl, kpgs)*jac
                    end do
                    vectu(i, n) = vectu(i, n)+defc(4, n, i)*sigas(4, kpgs)*jac*rac2+ &
                                  defn(4, n, i)*sigas(4, kpgs)*jac
                end do
            end do
!
            do n = 1, nno
                do i = 1, 2
                    do kl = 1, 3
                        vectu(i, n) = vectu(i, n)+defc(kl, n, i)*sigma(kl)*jac+ &
                                      defn(kl, n, i)*sigma(kl)*jac
                    end do
                    vectu(i, n) = vectu(i, n)+defc(4, n, i)*sigma(4)*jac+ &
                                  defn(4, n, i)*sigma(4)*jac/rac2
                end do
            end do
        end do
        do kl = 1, 3
            sigp(kl, kpg) = sigma(kl)
        end do
        sigp(4, kpg) = sigma(4)/rac2
        do i = 1, 6
            sigp(i+4, kpg) = qplus(i)
        end do
    end if
!
999 continue
! - SYNTHESE DES CODES RETOURS
    call codere(cod, npg, codret)
!
end subroutine
