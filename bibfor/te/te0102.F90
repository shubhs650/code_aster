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
subroutine te0102(option, nomte)
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystNone
    implicit none
!
#include "jeveux.h"
#include "asterfort/codent.h"
#include "asterfort/cq3d2d.h"
#include "asterfort/dfdm1d.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
!
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:  CALCUL DES MATRICES ELEMENTAIRES
!                          OPTION : 'MASS_THER      '
!                          CAS COQUE
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
! ......................................................................
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, ind, j, jgano, nbddl, nbres
    integer(kind=8) :: nbv, nbvar, ndimax
    real(kind=8) :: rocp, un, undemi
    parameter(ndimax=27)
    parameter(nbres=6)
    parameter(nbvar=2)
    integer(kind=8) :: icodre(nbres)
    character(len=2) :: num
    character(len=8) :: nompar(nbvar), fami, poum
    character(len=16) :: nomres(nbres)
    character(len=32) :: phenom
    real(kind=8) :: m(3, 3), h
    real(kind=8) :: valres(nbres)
    real(kind=8) :: coor2d(18)
    real(kind=8) :: dfdx(9), dfdy(9), poids, pm
    real(kind=8) :: mun, zero, deux, quatre
    real(kind=8) :: quinze, seize, cour, cosa, sina, r
    real(kind=8) :: valpar(nbvar), tempe, instan
    real(kind=8) :: masse(ndimax, ndimax), thickness
    integer(kind=8) :: nno, kp, npg2, gi, pi, gj, pj, k, imattt
    integer(kind=8) :: ipoids, ivf, idfde, igeom, kpg, spt
    integer(kind=8) :: imate, itemps, nnos, ndim
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!
    if (nomte .ne. 'THCPSE3 ' .and. nomte .ne. 'THCASE3 ') then
        call elrefe_info(fami='MASS', ndim=ndim, nno=nno, nnos=nnos, npg=npg2, &
                         jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
    else
        call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg2, &
                         jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
    end if
!
!
! --- INITIALISATIONS :
!     ---------------
    mun = -1.0d0
    zero = 0.0d0
    undemi = 0.5d0
    un = 1.0d0
    deux = 2.0d0
    quatre = 4.0d0
    quinze = 15.0d0
    seize = 16.0d0
    fami = 'FPG1'
    kpg = 1
    spt = 1
    poum = '+'
!
    tempe = zero
    instan = zero
    nompar(1) = 'INST'
    nompar(2) = 'TEMP'
    valpar(1) = instan
    valpar(2) = tempe
!
    do i = 1, 3
        do j = 1, 3
            m(i, j) = zero
        end do
    end do
!
    do i = 1, ndimax
        do j = 1, ndimax
            masse(i, j) = zero
        end do
    end do
!
!
! --- RECUPERATION DES COORDONNEES DES NOEUDS DE L'ELEMENT :
!     ----------------------------------------------------
    call jevech('PGEOMER', 'L', igeom)
!
! --- RECUPERATION DU MATERIAU :
!     ------------------------
    call jevech('PMATERC', 'L', imate)

! - Get plate parameters
    call getCara(plateCara, plateOrie)

! - No global<=>local transformation
    call compCoorSystNone(plateOrie)

    thickness = plateCara%thick

! --- RECUPERATION DE L'INSTANT DU CALCUL
    call jevech('PINSTR', 'L', itemps)
    valpar(1) = zr(itemps)
!
! --- RECUPERATION DE LA NATURE DU MATERIAU DANS PHENOM :
!     -------------------------------------------------
    call rccoma(zi(imate), 'THER', 1, phenom, icodre(1))
!
! --- DETERMINATION DU TENSEUR DE CAPACITE THERMIQUE :
!     ==============================================
!
! --- CAS DES COQUES THERMIQUES MULTI-COUCHES :
!     ---------------------------------------
    if (phenom .eq. 'THER_COQMU') then
!
! ---   NOM DES COMPOSANTES DU TENSEUR DE CAPACITE
! ---   THERMIQUE HOMOGENEISE :
!       ---------------------
        do i = 1, nbres
            call codent(i+24, 'G', num)
            nomres(i) = 'HOM_'//num
        end do
!
! ---   INTERPOLATION DES TERMES DU TENSEUR DE CAPACITE THERMIQUE
! ---   EN FONCTION DU TEMPS ET DE LA TEMPERATURE
! ---   (L'INTERPOLATION EN FONCTION DE LA TEMPERATURE EST
! ---    INACTIVE POUR LE MOMENT) :
!       -------------------------
        call rcvalb(fami, kpg, spt, poum, zi(imate), &
                    ' ', 'THER_COQMU', nbvar, nompar, valpar, &
                    nbres, nomres, valres, icodre, 1)
!
! ---   TENSEUR DE CAPACITE THERMIQUE :
!       -----------------------------
        m(1, 1) = valres(1)
        m(2, 1) = valres(2)
        m(3, 1) = valres(3)
        m(2, 2) = valres(4)
        m(3, 2) = valres(5)
        m(3, 3) = valres(6)
!
! --- CAS DES COQUES THERMIQUES ISOTROPES :
!     ===================================
    else if (phenom .eq. 'THER') then
!
! ---   INTERPOLATION DE LA CAPACITE THERMIQUE EN FONCTION DU TEMPS
! ---   ET DE LA TEMPERATURE
! ---   (L'INTERPOLATION EN FONCTION DE LA TEMPERATURE EST
! ---    INACTIVE POUR LE MOMENT) :
!       -------------------------
        nbv = 1
        nomres(1) = 'RHO_CP'
        call rcvalb(fami, kpg, spt, poum, zi(imate), &
                    ' ', 'THER', nbvar, nompar, valpar, &
                    nbv, nomres, valres, icodre, 1)
!
! ---   CAPACITE THERMIQUE :
!       ------------------
        rocp = valres(1)
!
! ---   DEMI-EPAISSEUR  :
!       --------------
        h = undemi*thickness
!
! ---   TENSEUR DE CAPACITE THERMIQUE :
!       -----------------------------
        m(1, 1) = seize*rocp*h/quinze
        m(2, 1) = deux*rocp*h/quinze
        m(3, 1) = deux*rocp*h/quinze
        m(2, 2) = quatre*rocp*h/quinze
        m(3, 2) = mun*rocp*h/quinze
        m(3, 3) = quatre*rocp*h/quinze
!
! --- CAS DES COQUES THERMIQUES HETEROGENES :
!     -------------------------------------
    else if (phenom .eq. 'THER_COQUE') then
!
! ---   NOM DES COMPOSANTES DU TENSEUR DE CAPACITE
! ---   THERMIQUE HOMOGENEISE
! ---   EN NOTANT RHOCP LA CAPACITE THERMIQUE EN CHAQUE POINT
! ---             P1(X3), P2(X3), P3(X3) LES POLYNOMES
! ---   DE LAGRANGE (OU AUTRES) DE DEGRE 2 RELATIFS A L'INTERPOLATION
! ---   DE LA TEMPERATURE DANS L'EPAISSEUR TELS QUE
! ---   P1 EST RELATIF A LA TEMPERATURE MOYENNE
! ---   P2 EST RELATIF A LA TEMPERATURE INFERIEURE
! ---   P3 EST RELATIF A LA TEMPERATURE SUPERIEURE
! ---   (I.E. T(X1,X2,X3) =    P1(X3)*TMOY(X1,X2)
! ---                        + P2(X3)*TINF(X1,X2)
! ---                        + P3(X3)*TSUP(X1,X2)
! ---   LES TERMES DU TENSEUR DE CAPACITE THERMIQUE HOMOGENEISE
! ---   SONT ALORS :
!       ----------
! ---   TERME SOMME_EPAISSEUR(RHOCP*P1(X3)*P1(X3).DX3) :
        nomres(1) = 'CMAS_MM'
! ---   TERME SOMME_EPAISSEUR(RHOCP*P1(X3)*P2(X3).DX3) :
        nomres(2) = 'CMAS_MP'
! ---   TERME SOMME_EPAISSEUR(RHOCP*P2(X3)*P2(X3).DX3) :
        nomres(3) = 'CMAS_PP'
! ---   TERME SOMME_EPAISSEUR(RHOCP*P2(X3)*P3(X3).DX3) :
        nomres(4) = 'CMAS_SI'
!
! ---   INTERPOLATION DES COMPOSANTES DU TENSEUR DE CAPACITE THERMIQUE
! ---   EN FONCTION DU TEMPS ET DE LA TEMPERATURE
! ---   (L'INTERPOLATION EN FONCTION DE LA TEMPERATURE EST
! ---    INACTIVE POUR LE MOMENT) :
!       -------------------------
        nbv = 4
        call rcvalb(fami, kpg, spt, poum, zi(imate), &
                    ' ', phenom, nbvar, nompar, valpar, &
                    nbv, nomres, valres, icodre, 1)
!
        m(1, 1) = valres(1)
        m(1, 2) = valres(2)
        m(1, 3) = valres(2)
        m(2, 2) = valres(3)
        m(2, 3) = valres(4)
        m(3, 3) = valres(3)
        m(2, 1) = m(1, 2)
        m(3, 1) = m(1, 3)
        m(3, 2) = m(2, 3)
!
    else
        call utmess('F', 'ELEMENTS3_17', sk=phenom)
    end if
!
!
!===================================
! --- CALCUL DE LA MASSE THERMIQUE =
!===================================
!
! --- CAS DES COQUES SURFACIQUES :
!     --------------------------
    if (nomte .ne. 'THCPSE3' .and. nomte .ne. 'THCASE3') then
!
! --- DETERMINATION DES COORDONNEES COOR2D DES NOEUDS DE L'ELEMENT
! --- DANS LE REPERE DE L'ELEMENT :
!     ---------------------------
        call cq3d2d(nno, zr(igeom), un, zero, coor2d)
!
! --- BOUCLE SUR LES POINTS D'INTEGRATION :
!     -----------------------------------
        do kp = 1, npg2
            k = (kp-1)*nno
!
! ---   DERIVEES DES FONCTIONS DE FORME ET PRODUIT JACOBIEN*POIDS
! ---   (DANS POIDS)  SUR L'ELEMENT :
!       ---------------------------
            call dfdm2d(nno, kp, ipoids, idfde, coor2d, &
                        poids, dfdx, dfdy)
            do gi = 1, nno
                do gj = 1, gi
                    do pi = 1, 3
                        do pj = 1, pi
                            pm = m(pi, pj)*zr(ivf+k+gi-1)*zr(ivf+k+gj-1)*poids
!
! ---     AFFECTATION DES TERMES HORS DIAGONAUX DE LA TRIANGULAIRE
! ---     INFERIEURE DE LA SOUS-MATRICE :
!         -----------------------------
                            if ((pi .ne. pj) .and. (gi .ne. gj)) then
                                i = 3*(gi-1)+pj
                                j = 3*(gj-1)+pi
                                masse(i, j) = masse(i, j)+pm
                            end if
!
! ---     AFFECTATION DES TERMES DE LA TRIANGULAIRE SUPERIEURE
! ---     DE LA SOUS-MATRICE :
!         ------------------
                            i = 3*(gi-1)+pi
                            j = 3*(gj-1)+pj
                            masse(i, j) = masse(i, j)+pm
                        end do
                    end do
                end do
            end do
        end do
!
    else
!
! --- CAS DES COQUES LINEIQUES :
!     ------------------------
!
! ---  BOUCLE SUR LES POINTS D'INTEGRATION :
!      -----------------------------------
        do kp = 1, npg2
            k = (kp-1)*nno
            call dfdm1d(nno, zr(ipoids+kp-1), zr(idfde+k), zr(igeom), dfdx, &
                        cour, poids, cosa, sina)
!
            if (nomte .eq. 'THCASE3') then
                r = zero
                do i = 1, nno
                    r = r+zr(igeom+2*i-2)*zr(ivf+k+i-1)
                end do
                poids = poids*r
            end if
!
            do gi = 1, nno
                do gj = 1, gi
                    do pi = 1, 3
                        do pj = 1, pi
                            pm = m(pi, pj)*zr(ivf+k+gi-1)*zr(ivf+k+gj-1)*poids
!
! ---     AFFECTATION DES TERMES HORS DIAGONAUX DE LA TRIANGULAIRE
! ---     INFERIEURE DE LA SOUS-MATRICE :
!         -----------------------------
                            if ((pi .ne. pj) .and. (gi .ne. gj)) then
                                i = 3*(gi-1)+pj
                                j = 3*(gj-1)+pi
                                masse(i, j) = masse(i, j)+pm
                            end if
!
! ---     AFFECTATION DES TERMES DE LA TRIANGULAIRE SUPERIEURE
! ---     DE LA SOUS-MATRICE :
!         ------------------
                            i = 3*(gi-1)+pi
                            j = 3*(gj-1)+pj
                            masse(i, j) = masse(i, j)+pm
                        end do
                    end do
                end do
            end do
        end do
!
    end if
!
! --- RECUPERATION DE LA MATRICE DE MASSE THERMIQUE EN SORTIE DU TE :
!     -------------------------------------------------------------
    call jevech('PMATTTR', 'E', imattt)
!
! --- AFFECTATION DE LA MATRICE DE MASSE THERMIQUE EN SORTIE DU TE :
!     ------------------------------------------------------------
    nbddl = 3*nno
    ind = 0
    do i = 1, nbddl
        do j = 1, i
            ind = ind+1
            zr(imattt+ind-1) = masse(i, j)
        end do
    end do
!
end subroutine
