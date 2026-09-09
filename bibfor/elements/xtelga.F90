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
subroutine xtelga(ndim, elrefp, nnop, igeom, tempno, &
                  lonch, cnset, jpintt, lsn, lst, &
                  heavn, basloc, heavt, nfh, nfe, &
                  temppg)
! aslint: disable=W1306
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/reeref.h"
#include "asterfort/xcalc_code.h"
#include "asterfort/xcalc_heav.h"
#include "MeshTypes_type.h"
#include "asterfort/xcalf2.h"
#include "asterfort/xcalfe.h"
!
    character(len=8) :: elrefp
    integer(kind=8) :: ndim, nnop, igeom, nfh, nfe, jpintt
    integer(kind=8) :: lonch(10), cnset(4*32), heavt(36), heavn(MT_NNOMAX3D, 5)
    real(kind=8) :: tempno(nnop*(1+nfh+nfe)), lsn(nnop), lst(nnop)
    real(kind=8) :: basloc(*), temppg(*)
!
!-----------------------------------------------------------------------
!
!     BUT: THERMIQUE LINEAIRE / ELEMENTS X-FEM LINEAIRES
!          CALCUL DE L'OPTION : 'TEMP_ELGA'
!
! IN :
! ---
! NDIM   --> DIMENSION DE L'ESPACE (2 OU 3)
! ELREFP --> NOM DE L'ELT PARENT DE REFERENCE
! NNOP   --> NBRE DE NOEUDS DE L'ELT PARENT DE REFERENCE
! IGEOM  --> ADRESSE DES COORDONEES DES NOEUDS DE L'ELT PARENT
! TEMPNO --> TEMPERATURE AUX NOEUDS
! LONCH  --> LONGUEURS DES CHAMPS UTILISES
! CNSET  --> CONNECTIVITE DES SOUS-ELEMENTS
! JPINTT --> ADRESSE DES COORDONEES DES POINTS D'INTERSECTION
! LSN    --> VALEUR DE LA LEVEL SET NORMALE AUX NOEUDS PARENTS
! LST    --> VALEUR DE LA LEVEL SET TANGENTIELLE AUX NOEUDS PARENTS
! BASLOC --> BASE LOCALE AU FOND DE FISSURE
! HEAVT  --> VALEURS DE L'HEAVISIDE SUR LES SS-ELTS
! NFH    --> NBRE DE FONCTION D'ENRICHISSEMENT HEAVISIDE (0 OU 1)
! NFE    --> NBRE DE FONCTION D'ENRICHISSEMENT CRACKTIP  (0 OU 1)
!
! OUT :
! ----
! TEMPPG --> TEMPERATURE ET SES DERIVEES PARTIELLES AUX POINTS DE GAUSS
!            DANS L'ORDRE : 'TEMP', 'DTX', 'DTY'          EN 2D
!                           'TEMP', 'DTX', 'DTY', 'DTZ'   EN 3D
!
!-----------------------------------------------------------------------
!
    character(len=8) :: elrese(3), fami(3)
    real(kind=8) :: baslog(3*ndim), tem, dtem(ndim), lsng, lstg, coorse(81), xg(ndim)
    real(kind=8) :: xe(ndim)
    real(kind=8) :: femec(4), dgdmec(4, ndim), feth, ff(nnop)
    real(kind=8) :: dgdth(ndim), dffenr(nnop, 1+nfh+nfe, ndim)
    real(kind=8) :: he
    real(kind=8) :: ffenr(nnop, 1+nfh+nfe), dfdi(nnop, ndim)
    integer(kind=8) :: ivf, kpg, nno, npg, j, iret, nse, ise, inp, in, ino, kddl
    integer(kind=8) :: nbddl, hea_se
    integer(kind=8) :: mxstac
!
    parameter(mxstac=1000)
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
!     NBRE DE DDLS PAR NOEUD
    nbddl = 1+nfh+nfe
!
!     SOUS-ELEMENT DE REFERENCE : RECUP DE NNO,NPG,IVF
    call elrefe_info(elrefe=elrese(ndim), fami=fami(ndim), nno=nno, npg=npg, jvf=ivf)
!
!     RECUPERATION DE LA SUBDIVISION DE L'ELEMENT EN NSE SOUS ELEMENT
    nse = lonch(1)
!
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
                        xe, ff, dfdi=dfdi)
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
                do j = 1, ndim
                    dgdth(j) = dgdmec(1, j)
                end do
            end if
! ------- FIN SI ENRICHISSEMENT SINGULIER
!
!         FFENR : TABLEAU DES FF ENRICHIES
            do inp = 1, nnop
!           DDL CLASSIQUE (TEMP)
                ffenr(inp, 1) = ff(inp)
                do j = 1, ndim
                    dffenr(inp, 1, j) = dfdi(inp, j)
                end do
!           DDL HEAVISIDE (H1)
                if (nfh .eq. 1) then
                    ffenr(inp, 1+nfh) = xcalc_heav(heavn(inp, 1), hea_se, heavn(inp, 5))*ff(inp)
                    do j = 1, ndim
                        dffenr(inp, 1+nfh, j) = xcalc_heav( &
                                                heavn(inp, 1), hea_se, heavn(inp, 5))*dfdi(inp, j)
                    end do
                end if
!           DDL CRACK-TIP (E1)
                if (nfe .eq. 1) then
                    ffenr(inp, 1+nfh+nfe) = feth*ff(inp)
                    do j = 1, ndim
                        dffenr(inp, 1+nfh+nfe, j) = feth*dfdi(inp, j)+ff(inp)*dgdth(j)
                    end do
                end if
            end do
!
!         CALCUL DE TEMP AU PG ET DE SES DERIVEES
            tem = 0.d0
            dtem = 0.d0
            do inp = 1, nnop
                do kddl = 1, nbddl
                    tem = tem+tempno(nbddl*(inp-1)+kddl)*ffenr(inp, kddl)
                    do j = 1, ndim
                        dtem(j) = dtem(j)+tempno(nbddl*(inp-1)+kddl)*dffenr(inp, kddl, j)
                    end do
                end do
            end do
!
!         ECRITURE DE TEMP ET DE SES DERIVEES AU PG
            temppg(npg*(ndim+1)*(ise-1)+(kpg-1)*ndim+kpg) = tem
            do j = 1, ndim
                temppg(npg*(ndim+1)*(ise-1)+(kpg-1)*ndim+kpg+j) = dtem(j)
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
