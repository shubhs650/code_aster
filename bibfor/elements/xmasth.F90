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
subroutine xmasth(ndim, elrefp, nnop, imate, itemps, &
                  igeom, lonch, cnset, jpintt, lsn, &
                  lst, heavn, basloc, heavt, nfh, &
                  nfe, mattt)
!.......................................................................
    implicit none
!
!     BUT: MATRICE DE MASSE ELEMENTAIRE EN THERMIQUE LINEAIRE
!          ELEMENTS X-FEM LINEAIRES
!
!          OPTION : 'MASS_THER'
!
! IN :
! ---
! NDIM   --> DIMENSION DE L'ESPACE (2 OU 3)
! ELREFP --> NOM DE L'ELT PARENT DE REFERENCE
! NNOP   --> NBRE DE NOEUDS DE L'ELT PARENT DE REFERENCE
! IMATE  --> ADRESSE DU MATERIAU
! ITEMPS --> ADRESSE DES PARAMETRES DE LA DICRETISATION EN TEMPS
! IGEOM  --> ADRESSE DES COORDONEES DES NOEUDS DE L'ELT PARENT
! LONCH  --> LONGUEURS DES CHAMPS UTILISES
! CNSET  --> CONNECTIVITE DES SOUS-ELEMENTS
! PINTT  --> ADRESSE DES COORDONEES DES POINTS D'INTERSECTION
! LSN    --> VALEUR DE LA LEVEL SET NORMALE AUX NOEUDS PARENTS
! LST    --> VALEUR DE LA LEVEL SET TANGENTIELLE AUX NOEUDS PARENTS
! BASLOC --> BASE LOCALE AU FOND DE FISSURE
! HEAVT  --> VALEURS DE L'HEAVISIDE SUR LES SS-ELTS
! NFH    --> NBRE DE FONCTION D'ENRICHISSEMENT HEAVISIDE (0 OU 1)
! NFE    --> NBRE DE FONCTION D'ENRICHISSEMENT CRACKTIP  (0 OU 1)
!
! OUT :
! ----
! MATTT  --> MATRICE DE MASSE ELEMENTAIRE
!.......................................................................
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/lteatt.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/reeref.h"
#include "asterfort/utmess.h"
#include "asterfort/xcalc_code.h"
#include "asterfort/xcalc_heav.h"
#include "asterfort/xcalf2.h"
#include "MeshTypes_type.h"
#include "asterfort/xcalfe.h"
!-----------------------------------------------------------------------
!
    character(len=8) :: elrefp
    integer(kind=8) :: ndim, nnop, imate, itemps, igeom, nfh, nfe, jpintt
    integer(kind=8) :: lonch(10), cnset(4*32), heavt(36), heavn(MT_NNOMAX3D, 5)
    real(kind=8) :: lsn(nnop), lst(nnop), basloc(*), mattt(*)
!
!-----------------------------------------------------------------------
!
    character(len=8) :: elrese(3), fami(3), poum
    character(len=16) :: phenom
    aster_logical :: axi
    real(kind=8) :: baslog(3*ndim), lsng, lstg, coorse(81), xg(ndim)
    real(kind=8) :: xe(ndim)
    real(kind=8) :: femec(4), dgdmec(4, ndim), feth, ff(nnop)
    real(kind=8) :: he
    real(kind=8) :: ffenr(nnop, 1+nfh+nfe), valpar(1), valres(1), jac
    real(kind=8) :: rhocp
    real(kind=8) :: r
    integer(kind=8) :: ivf, kpg, nno, npg, j, iret, nse, ise, inp, in, ino, kddl
    integer(kind=8) :: nbddl, hea_se
    integer(kind=8) :: mxstac, icodre(1), spt, ipoids, idfde, nosema, ind1, lddl, jnp
    integer(kind=8) :: iddlma
    integer(kind=8) :: ind2
!
    parameter(mxstac=1000)
!     NBRE MAX DE NOEUDS D'UN SOUS-ELEMENT (TRIA3,TETRA4,TRIA6 -> 6)
    parameter(nosema=6)
!
    real(kind=8) :: r8bid1(nosema), r8bid2(nosema), r8bid3(nosema)
!
    data elrese/'SE2', 'TR3', 'TE4'/
    data fami/'BID', 'XINT', 'XINT'/
!
!-----------------------------------------------------------------------
!
!     VERIF QUE LES TABLEAUX LOCAUX DYNAMIQUES NE SONT PAS TROP GRANDS
!     (VOIR CRS 1404)
    ASSERT(nnop .le. mxstac)
!
!     S'AGIT-IL D'UNE MODELISATION AXIS
    axi = .false.
    if (lteatt('AXIS', 'OUI')) axi = .true.
!
!     NBRE DE DDLS PAR NOEUD
    nbddl = 1+nfh+nfe
!
!
!     POUR PREPARER L'APPEL A RCVALB
    call rccoma(zi(imate), 'THER', 1, phenom, icodre(1))
    if (icodre(1) .ne. 0) then
        call utmess('F', 'ELEMENTS2_63')
    end if
!     POUR L'INSTANT ON NE TRAITE PAS 'THER_ORTH'
    ASSERT(phenom .eq. 'THER')
    valpar(1) = zr(itemps-1+1)
    spt = 1
    poum = '+'
!
!     SOUS-ELEMENT DE REFERENCE : RECUP DE NNO,NPG,IPOIDS,IVF,IDFDE
    call elrefe_info(elrefe=elrese(ndim), fami=fami(ndim), nno=nno, npg=npg, jpoids=ipoids, &
                     jvf=ivf, jdfde=idfde)
!
!     RECUPERATION DE LA SUBDIVISION DE L'ELEMENT EN NSE SOUS ELEMENT
    nse = lonch(1)
!
! ----------------------------------------------------------------------
! --- BOUCLE SUR LES NSE SOUS-ELEMENTS
! ----------------------------------------------------------------------
!
    do ise = 1, nse
!
!       VALEUR (CSTE) DE LA FONCTION HEAVISIDE SUR LE SS-ELT
        he = 1.d0*heavt(ise)
        hea_se = xcalc_code(1, he_real=[he])
!
!       BOUCLE SUR LES SOMMETS DU SOUS-TETRA/TRIA -> COORDS NOEUDS
        do in = 1, nno
            ino = cnset(nno*(ise-1)+in)
            do j = 1, ndim
                if (ino .lt. 1000) then
                    coorse(ndim*(in-1)+j) = zr(igeom-1+ndim*(ino-1)+j)
                else if (ino .gt. 1000 .and. ino .lt. 2000) then
                    coorse(ndim*(in-1)+j) = zr(jpintt-1+ndim*(ino-1000- &
                                                              1)+j)
                else
                    ASSERT(.false.)
                end if
            end do
        end do
!
! ----------------------------------------------------------------------
! ----- BOUCLE SUR LES POINTS DE GAUSS
! ----------------------------------------------------------------------
!
        do kpg = 1, npg
!
!         COORDONNÉES DU PT DE GAUSS DANS LE REPÈRE RÉEL : XG
            xg(:) = 0.d0
            do j = 1, ndim
                do in = 1, nno
                    xg(j) = xg(j)+zr(ivf-1+nno*(kpg-1)+in)*coorse(ndim*( &
                                                                  in-1)+j)
                end do
            end do
!
!         XG -> XE (DANS LE REPERE DE l'ELREFP) ET VALEURS DES FF EN XE
            call reeref(elrefp, nnop, zr(igeom), xg, ndim, &
                        xe, ff)
!
! ------- SI ENRICHISSEMENT SINGULIER
            if (nfe .gt. 0) then
!           BASE LOCALE ET LEVEL SETS AU POINT DE GAUSS
                baslog(:) = 0.d0
                lsng = 0.d0
                lstg = 0.d0
                do inp = 1, nnop
                    lsng = lsng+lsn(inp)*ff(inp)
                    lstg = lstg+lst(inp)*ff(inp)
                    do j = 1, 3*ndim
                        baslog(j) = baslog(j)+basloc(3*ndim*(inp-1)+j)*ff(inp)
                    end do
                end do
!           FONCTION D'ENRICHISSEMENT (MECA) AU PG ET DÉRIVÉES
                if (ndim .eq. 2) then
                    call xcalf2(he, lsng, lstg, baslog, femec, &
                                dgdmec, iret)
                else if (ndim .eq. 3) then
                    call xcalfe(he, lsng, lstg, baslog, femec, &
                                dgdmec, iret)
                end if
!           PB DE CALCUL DES DERIVEES DES FONCTIONS SINGULIERES
!           CAR ON SE TROUVE SUR LE FOND DE FISSURE
                ASSERT(iret .ne. 0)
!           ON NE GARDE QUE LES ENRICHISSEMENTS UTILES EN THERMIQUE
                feth = femec(1)
            end if
! ------- FIN SI ENRICHISSEMENT SINGULIER
!
!         CALCULER LE JACOBIEN DE LA TRANSFO SSTET->SSTET REF
!         AVEC LES COORDONNEES DU SOUS-ELEMENT
            if (ndim .eq. 2) then
                call dfdm2d(nno, kpg, ipoids, idfde, coorse, &
                            jac, r8bid1, r8bid2)
            else if (ndim .eq. 3) then
                call dfdm3d(nno, kpg, ipoids, idfde, coorse, &
                            jac, r8bid1, r8bid2, r8bid3)
            end if
!
!         MODIFICATION DU JACOBIEN SI AXI
            if (axi) then
                r = 0.d0
                do inp = 1, nnop
                    r = r+ff(inp)*zr(igeom-1+2*(inp-1)+1)
                end do
                ASSERT(r .gt. 0d0)
                jac = jac*r
            end if
!
!         RECUPERER LES PARAMETRES MATERIAUX
            call rcvalb('XFEM', kpg, spt, poum, zi(imate), &
                        ' ', phenom, 1, 'INST', [valpar], &
                        1, 'RHO_CP', valres, icodre, 1)
            rhocp = valres(1)
!
!         FFENR : TABLEAU DES FF ENRICHIES
            do inp = 1, nnop
!           DDL CLASSIQUE (TEMP)
                ffenr(inp, 1) = ff(inp)
!           DDL HEAVISIDE (H1)
                if (nfh .eq. 1) then
                    ffenr(inp, 1+nfh) = xcalc_heav(heavn(inp, 1), hea_se, heavn(inp, 5))*ff(inp)
                end if
!           DDL CRACK-TIP (E1)
                if (nfe .eq. 1) then
                    ffenr(inp, 1+nfh+nfe) = feth*ff(inp)
                end if
            end do
!
! ------- REMPLISSAGE DE LA MATRICE DE MASSE
!
            do inp = 1, nnop
!
                do kddl = 1, nbddl
!
                    ind1 = (nbddl*(inp-1)+kddl-1)*(nbddl*(inp-1)+kddl)/2
!
                    do lddl = 1, nbddl
!
                        do jnp = 1, inp
!
!                 IDDLMA : NUMERO DE DDL MAX PR NE PAS DEPASSER LA
!                 DIGAONALE (ON STOCKE LA PARTIE TRIANGULAIRE INF)
                            if (jnp .eq. inp) then
                                iddlma = kddl
                            else
                                iddlma = nbddl
                            end if
!
!                 ON NE DEPASSE PAS PAS LA DIAGONALE
                            if (lddl .le. iddlma) then
                                ind2 = ind1+nbddl*(jnp-1)+lddl
                                mattt(ind2) = mattt(ind2)+rhocp*jac*ffenr(inp, kddl)*ffen&
                                              &r(jnp, lddl)
                            end if
!
                        end do
                    end do
                end do
            end do
!
        end do
!
! ----------------------------------------------------------------------
! ----- FIN BOUCLE SUR LES POINTS DE GAUSS
! ----------------------------------------------------------------------
!
    end do
!
! ----------------------------------------------------------------------
! --- FIN BOUCLE SUR LES SOUS-ELEMENTS
! ----------------------------------------------------------------------
!
end subroutine
