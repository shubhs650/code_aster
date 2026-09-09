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
subroutine xsifle(ndim, ifa, jptint, cface, igeom, &
                  nfh, jheavn, singu, nfe, ddlc, &
                  ddlm, jlsn, jlst, jstno, ipres, &
                  ipref, itemps, idepl, nnop, valres, &
                  basloc, ithet, nompar, option, igthet, &
                  jbasec, contac)
!
!
! aslint: disable=W1306,W1504
    implicit none
#include "jeveux.h"
#include "asterc/r8pi.h"
#include "asterfort/assert.h"
#include "asterfort/elelin.h"
#include "asterfort/elref1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/lteatt.h"
#include "asterfort/tecach.h"
#include "asterfort/xjacf2.h"
#include "MeshTypes_type.h"
#include "asterfort/xjacff.h"
#include "asterfort/xsifl1.h"
#include "asterfort/xxmmvd.h"
!
    character(len=8) :: nompar(4)
    character(len=16) :: option
    integer(kind=8) :: ndim, ifa, cface(30, 6), igeom, nfh, singu, jlst, ipres
    integer(kind=8) :: nfe, ddlc, ipref, itemps, nnop, ithet, jptint, igthet, idepl
    integer(kind=8) :: ddlm, jbasec, contac, jheavn
    integer(kind=8) :: jstno, jlsn
    real(kind=8) :: valres(3)
    real(kind=8) :: basloc(9*nnop)
!
!
!    - FONCTION REALISEE:  CALCUL DES OPTIONS DE POST-TRAITEMENT
!                          EN MÉCANIQUE DE LA RUPTURE
!                          SUR LES LEVRES DES FISSURES X-FEM
!
! IN  NDIM  : DIMENSION DU PROBLEME
! IN  IFA   : NUMERO DE LA FACETTE DE CONTACT DS L ELEMENT PARENT
! IN JPTINT : ADRESSE TOPOFAC.PI - COORDONNEES PTS D'INTERSECTION
! IN JAINT  : ADRESSE TOPOFAC.AI - ARETES INTERSECTEES
! IN CFACE  : TABLEAU TOPOFAC.CF - CONNECTIVITE DES FACETTES
! IN IGEOM  : ADRESSE COORDONNEES NOEUDS MAILLE PARENT
! IN NFH    : NOMBRE ENRICHISSEMENT HEAVISIDE
! IN SINGU  : PRESENCE ENRICHISSEMENT CRACK TIP
! IN NFE    : NOMBRE ENRICHISSEMENTS CRACK TIP
! IN DDLC   : NOMBRE DDLS CONTACT SUR CHAQUE NOEUD QUI EN PORTE
! IN DDLM   : NOMBRE DDL SUR CHAQUE NOEUDS MILIEU
! IN JLST   : ADRESSE LEVEL-SET TANGENTE
! IN IPRES  : PRESSION SI CONSTANTE
! IN IPREF  : PRESSION SI FONCTION
! IN ITEMPS : CARTE (CONSTANTE) AVEC L INSTANT
!        UTILISE SI PRESSION FONCTION DU TEMPS
! IN IDEPL  : ADRESSE DEPLACEMENT
! IN NNOP   : NOMBRE DE NOEUDS MAILLE PARENT
! IN VALRES : CONSTANTES ELASTIQUES MATERIAU /E, NU, ALPHA/
! IN BASLOC : FISS.BASLOC -  BASE DU PROJETE AU NOEUD
! IN ITHET  : CHAMP THETA AUX NOEUDS DU MAILLAGE PARENT
! IN NOMPAR : NOM DES PARAMETRES SI PRESSION VARIABLE
! IN PRESN  : PRESSION AUX NOEUDS PARENT
! IN OPTION : OPTION DE CALCUL
! OUT IGTHET: G (OPTION CALC_G) ET K1, K2, K3 (SI OPTION CALC_K_G)
! IN JBASEC : ADRESSE TOPOFAC.BA
! IN CONTAC : TYPE DE CONTACT P1P1 ou P2P1
!
!
    integer(kind=8) :: nnof, npgf, ipoidf, ivff, idfdef
    integer(kind=8) :: ipgf, zxain, heavn(nnop, 5)
    integer(kind=8) :: ddld, ddls, ncompn, ino, ig, iret, jtab(7)
    real(kind=8) :: xg(3), jac, ff(MT_NNOMAX3D), nd(3)
    real(kind=8) :: angl(2)
    real(kind=8) :: e, nu, mu, ka, coeff, coeff3, rnbid(MT_NNOMAX3D)
    real(kind=8) :: dfdi(nnop, ndim)
    real(kind=8) :: he(2)
!
    integer(kind=8) :: nnos, nno
    character(len=8) :: elref, fpg, elc, elrefc
    real(kind=8) :: tau1(3), tau2(3)
    data he/-1.d0, 1.d0/
!
!
    zxain = xxmmvd('ZXAIN')
!     PAR CONVENTION :
!     LEVRE INFERIEURE (HE=-1) EST LA LEVRE 1, DE NORMALE SORTANTE  ND
!     LEVRE SUPERIEURE (HE=+1) EST LA LEVRE 2, DE NORMALE SORTANTE -ND
    angl(1) = -r8pi()
    angl(2) = r8pi()
!
!     RECUPERATION DE LA DEFINITION DES DDL HEAVISIDES
    if (nfh .gt. 0) then
        call tecach('OOO', 'PHEA_NO', 'L', iret, nval=7, &
                    itab=jtab)
        ncompn = jtab(2)/jtab(3)
        ASSERT(ncompn .eq. 5)
        do ino = 1, nnop
            do ig = 1, ncompn
                heavn(ino, ig) = zi(jheavn-1+ncompn*(ino-1)+ig)
            end do
        end do
    end if
!
!     NOMBRE DE DDL DE DEPLACEMENT À CHAQUE NOEUD SOMMET
    ddld = ndim*(1+nfh+nfe)
!
!     NOMBRE DE DDL TOTAL (DEPL+CONTACT) À CHAQUE NOEUD SOMMET
    ddls = ddld+ddlc
!
    call elref1(elref)
!
    if (ndim .eq. 3) then
        elc = 'TR3'
        fpg = 'XCON'
    else if (ndim .eq. 2) then
        elc = 'SE2'
        fpg = 'MASS'
    end if
!
    call elrefe_info(elrefe=elc, fami=fpg, nno=nnof, npg=npgf, jpoids=ipoidf, &
                     jvf=ivff, jdfde=idfdef)
!
!     MATÉRIAU HOMOGENE
!     ON PEUT PAS LE RECUPERER SUR LES POINTS DE GAUSS DES FACETTES
!     CAR LA FAMILLE CONCATENEE DES PG DES FACETTES N'EXISTE PAS
    e = valres(1)
    nu = valres(2)
    mu = e/(2.d0*(1.d0+nu))
!
!   DEFINITION DE KAPPA
    if (ndim .eq. 2 .and. lteatt('C_PLAN', 'OUI')) then
        ka = (3.d0-nu)/(1.d0+nu)
    else
        ka = 3.d0-4.d0*nu
    end if
    coeff = e/(1.d0-nu*nu)
    coeff3 = 2.d0*mu
!
!     ----------------------------------------------------------------
!     BOUCLE SUR LES POINTS DE GAUSS DES FACETTES
    do ipgf = 1, npgf
!
!       CALCUL DE JAC (PRODUIT DU JACOBIEN ET DU POIDS)
!       ET DES FF DE L'ÉLÉMENT PARENT AU POINT DE GAUSS
!       ET LA NORMALE ND ORIENTÉE DE ESCL -> MAIT
!       ET DE XG : COORDONNEES REELLES DU POINT DE GAUSS
!       ET DE DFDI : DERIVES DES FF PARENT
        call elelin(contac, elref, elrefc, nno, nnos)
        elrefc = 'NON'
        ff(:) = 0.d0
        if (ndim .eq. 3) then
            ASSERT(nno .eq. nnop)
            call xjacff(elref, elrefc, elc, ndim, fpg, &
                        jptint, ifa, cface, ipgf, nnop, &
                        nnos, igeom, jbasec, xg, jac, &
                        ff, rnbid, dfdi, nd, tau1, &
                        tau2)
        else if (ndim .eq. 2) then
            call xjacf2(elref, elrefc, elc, ndim, fpg, &
                        jptint, ifa, cface, ndim, ipgf, &
                        nnop, nnos, igeom, jbasec, xg, &
                        jac, ff, rnbid, dfdi, nd, &
                        tau1)
        end if
        call xsifl1(elref, angl, basloc, coeff, coeff3, &
                    ddlm, ddls, dfdi, ff, he, &
                    heavn, idepl, igthet, ipref, ipres, &
                    ithet, jac, jlsn, jlst, jstno, &
                    ka, mu, nd, ndim, nfh, &
                    nnop, nnos, itemps, nompar, option, &
                    singu, xg, igeom)
    end do
!     FIN DE BOUCLE SUR LES POINTS DE GAUSS DES FACETTES
!     ----------------------------------------------------------------
!
end subroutine
