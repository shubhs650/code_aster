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
! aslint: disable=W1504,W1306
!
subroutine cgfint(BEHInteg, &
                  typmod, option, &
                  comporKit, carcrikit, &
                  ndim, nno1, nno2, npg, &
                  wref, vff1, vff2, dffr1, &
                  geom, tang, &
                  instam, instap, &
                  ddlm, ddld, &
                  iu, iuc, im, a, &
                  lgpg, sigm, vim, &
                  sigp, vip, &
                  matr, vect, &
                  codret)
!
    use Behaviour_type
    use Behaviour_module
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/cgcine.h"
#include "asterfort/codere.h"
#include "asterfort/nmcomp.h"
#include "asterfort/nmiclg.h"
#include "asterfort/r8inir.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
!
    type(Behaviour_Integ), intent(inout) :: BEHInteg
    character(len=8), intent(in) :: typmod(2)
    character(len=16), intent(in) :: option, comporKit(COMPOR_SIZE)
    real(kind=8), intent(in) :: carcrikit(CARCRI_SIZE)
    integer(kind=8) :: ndim, nno1, nno2, npg, lgpg, iu(3, 3), iuc(3), im(3)
    integer(kind=8) :: codret
    real(kind=8) :: vff1(nno1, npg), vff2(nno2, npg), geom(ndim, nno1), wref(npg)
    real(kind=8) :: instam, instap, sigm(3, npg)
    real(kind=8) :: ddlm(nno1*(ndim+1)+nno2), ddld(nno1*(ndim+1)+nno2)
    real(kind=8) :: vim(lgpg, npg), vip(lgpg, npg), vect(nno1*(ndim+1)+nno2)
    real(kind=8) :: dffr1(nno1, npg), tang(*), sigp(3, npg), matr(*)
    real(kind=8) :: a
!
! --------------------------------------------------------------------------------------------------
!
!   RAPH_MECA, RIGI_MECA_* ET FULL_MECA_* POUR L'ELEMENT CABLE/GAINE
!
! --------------------------------------------------------------------------------------------------
!
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  NNO1    : NOMBRE DE NOEUDS (FAMILLE U)
! IN  NNO2    : NOMBRE DE NOEUDS (FAMILLE L)
! IN  NPG     : NOMBRE DE POINTS DE GAUSS
! IN  WREF    : POIDS DES POINTS DE GAUSS DE REFERENCE
! IN  VFF1    : VALEUR DES FONCTIONS DE FORME (FAMILLE U)
! IN  VFF2    : VALEUR DES FONCTIONS DE FORME (FAMILLE L)
! IN  DFFR1   : DERIVEES DES FONCTIONS DE FORME DE REFERENCE (FAMILLE U)
! IN  GEOM    : COORDONNEES DES NOEUDS
! IN  TANG    : TANGENTE AUX NOEUDS
! IN  TYPMOD  : TYPE DE MODELISATION
! IN  OPTION  : OPTION DE CALCUL
! IN  LGPG    : LONGUEUR DU TABLEAU DES VARIABLES INTERNES
! IN  CRIT    : CRITERES DE CONVERGENCE LOCAUX
! IN  INSTAM  : INSTANT PRECEDENT
! IN  INSTAP  : INSTANT DE CALCUL
! IN  DDLM    : DDL A L'INSTANT PRECEDENT
! IN  DDLD    : INCREMENT DES DDL
! IN  IU      : DECALAGE D'INDICE POUR ACCEDER AUX DDL DE DEPLACEMENT GA
! IN  IUC     : DECALAGE D'INDICE POUR ACCEDER AUX DDL DE DEPL RELATIF
! IN  IM      : DECALAGE D'INDICE POUR ACCEDER AUX DDL DE LAGRANGE
! IN  A       : SECTION DE LA BARRE
! IN  SIGM    : CONTRAINTES A L'INSTANT PRECEDENT
! IN  VIM     : VARIABLES INTERNES A L'INSTANT PRECEDENT
! OUT SIGP    : CONTRAINTES DE CAUCHY (RAPH_MECA   ET FULL_MECA_*)
! OUT VIP     : VARIABLES INTERNES    (RAPH_MECA   ET FULL_MECA_*)
! OUT MATR    : MATRICE DE RIGIDITE   (RIGI_MECA_* ET FULL_MECA_*)
! OUT VECT    : FORCES INTERIEURES    (RAPH_MECA   ET FULL_MECA_*)
! OUT CODRET  : CODE RETOUR
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: multComp = " "
    integer(kind=8), parameter:: nepsSheath = 2, nsigSheath = 1, ndsdeSheath = 2
    integer(kind=8), parameter :: ksp = 1
    integer(kind=8), parameter :: nbProp = 1
    integer(kind=8) :: propCode(nbProp)
    real(kind=8) :: propVale(nbProp)
    character(len=16), parameter :: propName(nbProp) = (/'PENA_LAGR'/)
    aster_logical  :: lMatr, lSigm, lVari
    character(len=16) :: relaSheath, relaCable
    integer(kind=8) :: numeSheath, nbviSheath, nbviCable
    aster_logical :: resi, rigi
    integer(kind=8) :: nddl, kpg, cod(npg), n, i, m, j, kk
    integer(kind=8) :: nume
    real(kind=8) :: r, mu, epsm, deps, wg, l(3), de, ddedt, t1
    real(kind=8) :: epsmSheath(nepsSheath), depsSheath(nepsSheath)
    real(kind=8) :: sigmSheath(nsigSheath), sigpSheath(nsigSheath)
    real(kind=8) :: dsdeSheath(nsigSheath, nepsSheath)
    real(kind=8) :: b(4, 3), gliss
    real(kind=8) :: sigcab, dsidep, dde(nepsSheath), ddedn, courb
    character(len=16) :: comporSheath(COMPOR_SIZE)
    character(len=1) :: poum
! ----------------------------------------------------------------------
!

! - INITIALISATION
    resi = option(1:4) .eq. 'FULL' .or. option(1:4) .eq. 'RAPH'
    rigi = option(1:4) .eq. 'FULL' .or. option(1:4) .eq. 'RIGI'
    nddl = nno1*(ndim+1)+nno2
    cod = 0

! - Prepare compor maps
    ASSERT(comporKit(RELA_NAME) .eq. 'KIT_CG')
    relaSheath = comporKit(SHEATH_NAME)
    relaCable = comporKit(CABLE_NAME)

! - Get parameters
    read (comporKit(NUME), '(I16)') nume
    read (comporKit(SHEATH_NUME), '(I16)') numeSheath
    read (comporKit(SHEATH_NVAR), '(I16)') nbviSheath
    read (comporKit(CABLE_NVAR), '(I16)') nbviCable
    comporSheath = comporKit
    comporSheath(RELA_NAME) = relaSheath
    write (comporSheath(NUME), '(I16)') numeSheath
    write (comporSheath(NVAR), '(I16)') nbviSheath
    BEHInteg%behavPara%nvi = nbviSheath
    BEHInteg%behavPara%numlc = numeSheath
!
    if (rigi) call r8inir(nddl*nddl, 0.d0, matr, 1)
    if (resi) call r8inir(nddl, 0.d0, vect, 1)

    lSigm = L_SIGM(option)
    lVari = L_VARI(option)
    lMatr = L_MATR(option)

    do kpg = 1, npg

! ----- Set main parameters for behaviour (on point)
        call behaviourSetParaPoin(kpg, ksp, BEHInteg)

! ----- CALCUL DES ELEMENTS GEOM DE L'EF AU POINT DE GAUSS CONSIDERE
        call cgcine(ndim, nno1, vff1(1, kpg), wref(kpg), dffr1(1, kpg), &
                    geom, tang, wg, l, b, &
                    courb)

!      CALCUL DU DEPLACEMENT TANGENTIEL DE LA GAINE (UPROJ)
!                DEPLACEMENT D
! ---------------U CABLE               (UCAB)
!             DE LA DEFORMATION DANS LE CABLE       (EPSCAB)
!             DU LAGRANGE                           (MU)
!         AU POINT DE GAUSS CONSIDERE
        gliss = 0.d0
        epsm = 0.d0
        deps = 0.d0
        mu = 0.d0
        do n = 1, nno1
            do j = 1, ndim
                epsm = epsm+b(j, n)*ddlm(iu(j, n))
                deps = deps+b(j, n)*ddld(iu(j, n))
            end do
            gliss = gliss+l(n)*(ddlm(iuc(n))+ddld(iuc(n)))
            epsm = epsm+b(ndim+1, n)*ddlm(iuc(n))
            deps = deps+b(ndim+1, n)*ddld(iuc(n))
        end do
!
        mu = 0.d0
        do n = 1, nno2
            mu = mu+vff2(n, kpg)*(ddlm(im(n))+ddld(im(n)))
        end do
!
!
!      LOI DE COMPORTEMENT DU CABLE
!         SORTIE DE CONTRAINTE DANS LE CABLE     (SIGCAB)
!                DES VARIABLES INTERNES DU CABLE (VIP)
!                DE LA TANGENT
! ---------------E D(SIGCAB)/D(EPSCAB)
!
        if (relaCable .eq. 'ELAS' .or. relaCable .eq. 'VMIS_ISOT_TRAC' .or. &
            relaCable .eq. 'VMIS_ISOT_LINE' .or. &
            relaCable .eq. 'CORR_ACIER' .or. relaCable .eq. 'VMIS_CINE_LINE' .or. &
            relaCable .eq. 'VMIS_ASYM_LINE' .or. relaCable .eq. 'SANS') then
            call nmiclg(BEHInteg%materPara, &
                        option, relaCable, carcriKit, &
                        epsm, deps, &
                        sigm(1, kpg)/a, vim(1, kpg), &
                        sigcab, vip(1, kpg), dsidep, &
                        codret)
        else
            call utmess('F', 'CABLE0_26')
        end if
!
!      LOI DE COMPORTEMENT: GLISSEMENT GAINE CABLE
!         ADHERENCE PARFAITE OU GLISSEMENT PARFAIT
!         SORTIE DE GLISSEMENT RELATIF AU POINT DE GAUSS (DE)
!                DE LA TANGENTE D(DE)/D(TAU) AVEC
!                   TAU=L+MU(GLISS)
!      CONVENTIONS :
!       1. MU EST RANGE DANS EPSM(1)
!       2. GLISS EST RANGE DANS EPSD(1)
!       3. DELTA EST RENVOYE DANS SIGP(1)             : DE
!       4. D(DE)/DTAU EST RENVOYE DANS DSIDEP(1,1) : DDEDT
!       5. TENSION ET COURBURE : DANS TABLEAU DE TRAVAIL
!
!    RECUPERATION DES PARAMETRES PHYSIQUES
        if (option .eq. 'RIGI_MECA_TANG') then
            poum = '-'
        else
            poum = '+'
        end if
!
        call rcvalb(BEHInteg%materPara%schemePara%fami, &
                    BEHInteg%materPara%schemePara%kpg, &
                    BEHInteg%materPara%schemePara%ksp, &
                    poum, &
                    BEHInteg%materPara%jvMaterCode, &
                    ' ', 'CABLE_GAINE_FROT', &
                    0, ' ', [0.d0], &
                    nbProp, propName, propVale, &
                    propCode, 2)
        r = propVale(1)

! ----- Set cable parameters (on point)
        BEHInteg%behavESVA%behavESVAOther%tenscab = a*sigcab
        BEHInteg%behavESVA%behavESVAOther%curvcab = courb
!
        de = 0.d0
        epsmSheath = [mu, gliss]
        depsSheath = [0.d0, 0.d0]
        sigmSheath = [0.d0]
        call nmcomp(BEHInteg, &
                    ndim, option, typmod, &
                    instam, instap, &
                    comporSheath, carcrikit, multComp, &
                    nepsSheath, epsmSheath, depsSheath, &
                    nsigSheath, sigmSheath, &
                    vim(nbviCable+1, kpg), &
                    sigpSheath, vip(nbviCable+1, kpg), &
                    ndsdeSheath, dsdeSheath, cod(kpg))

        if (cod(kpg) .ne. 0) goto 999
        if (lSigm) de = sigpSheath(1)
        if (lMatr) dde = dsdeSheath(1, :)
!
!      FORCE INTERIEURE ET CONTRAINTES DE CAUCHY
!
        if (resi) then
!
!        STOCKAGE DES CONTRAINTES
!        CONVENTION DE RANGEMENT SIGP(1,2,3) EXPLICITE CI-DESSOUS
            sigp(1, kpg) = sigcab*a
            sigp(2, kpg) = mu+r*(gliss-de)
            sigp(3, kpg) = gliss-de
!
!        VECTEUR FINT:U ET UC
            do n = 1, nno1
                do i = 1, ndim
                    kk = iu(i, n)
                    t1 = b(i, n)*sigp(1, kpg)
                    vect(kk) = vect(kk)+wg*t1
                end do
                kk = iuc(n)
                t1 = b(4, n)*sigp(1, kpg)+l(n)*sigp(2, kpg)
                vect(kk) = vect(kk)+wg*t1
            end do
!
!        VECTEUR FINT:M
            do n = 1, nno2
                kk = im(n)
                t1 = vff2(n, kpg)*sigp(3, kpg)
                vect(kk) = vect(kk)+wg*t1
            end do
!
!
        end if
!
!
! - CALCUL DE LA MATRICE DE RIGIDITE
!   STOCKAGE TRIANGLE INFERIEUR LIGNE DE DFI/DUJ
!
        if (rigi) then
!        MATRICE K:U(I,N),U(J,M)
            ddedt = dde(1)
            ddedn = a*dde(2)
            do n = 1, nno1
                do i = 1, ndim
                    do m = 1, nno1
                        do j = 1, ndim
                            kk = (iu(i, n)-1)*nddl+iu(j, m)
                            t1 = b(i, n)*b(j, m)*dsidep*a
                            matr(kk) = matr(kk)+wg*t1
                        end do
                    end do
                end do
            end do
!
!        MATRICE K:U(I,N),UC(M)
            do n = 1, nno1
                do i = 1, ndim
                    do m = 1, nno1
                        t1 = b(i, n)*b(4, m)*dsidep*a
                        kk = (iu(i, n)-1)*nddl+iuc(m)
                        matr(kk) = matr(kk)+wg*t1
                        t1 = t1-r*l(m)*ddedn*dsidep*b(i, n)
                        kk = (iuc(m)-1)*nddl+iu(i, n)
                        matr(kk) = matr(kk)+wg*t1
                    end do
                end do
            end do
!
!
!        MATRICES K:UC(N),UC(M)
            do n = 1, nno1
                do m = 1, nno1
                    t1 = b(4, n)*b(4, m)*dsidep*a+r*l(n)*l(m)*(1.d0-r* &
                                                               ddedt)-r*l(n)*ddedn*dsidep*b(4, m)
                    kk = (iuc(n)-1)*nddl+iuc(m)
                    matr(kk) = matr(kk)+wg*t1
                end do
            end do
!
!        MATRICE K:UC(N),MU(M)
            do n = 1, nno1
                do m = 1, nno2
                    t1 = l(n)*vff2(m, kpg)*(1.d0-r*ddedt)
                    kk = (iuc(n)-1)*nddl+im(m)
                    matr(kk) = matr(kk)+wg*t1
                    t1 = t1-vff2(m, kpg)*ddedn*dsidep*b(4, n)
                    kk = (im(m)-1)*nddl+iuc(n)
                    matr(kk) = matr(kk)+wg*t1
                end do
            end do
!
!        MATRICES K:MU(N),MU(M)
            do n = 1, nno2
                do m = 1, nno2
                    t1 = -vff2(n, kpg)*ddedt*vff2(m, kpg)
                    kk = (im(n)-1)*nddl+im(m)
                    matr(kk) = matr(kk)+wg*t1
                end do
            end do
!
!        MATRICE K:U(N),MU(M)
            do n = 1, nno1
                do i = 1, ndim
                    do m = 1, nno2
                        t1 = -vff2(m, kpg)*ddedn*dsidep*b(i, n)
                        kk = (im(m)-1)*nddl+iu(i, n)
                        matr(kk) = matr(kk)+wg*t1
                    end do
                end do
            end do
!
        end if
!
    end do
!
! - SYNTHESE DES CODES RETOUR
!
999 continue
    call codere(cod, npg, codret)
!
end subroutine
