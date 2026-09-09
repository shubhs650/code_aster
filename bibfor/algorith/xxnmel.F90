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
! aslint: disable=W1306,W1504
!
subroutine xxnmel(typmod, materPara, &
                  elrefp, elrese, ndim, coorse, &
                  jvGeom, he, nfh, ddlc, ddlm, &
                  nnops, nfe, basloc, nnop, npg, &
                  lsn, lst, idecpg, &
                  matuu, &
                  nfiss, heavn, jstno)
!
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!lgpg,
#include "asterc/r8vide.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/dmatmc.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/iimatu.h"
#include "asterfort/indent.h"
#include "asterfort/reeref.h"
#include "asterfort/xcalc_code.h"
#include "asterfort/xcalc_heav.h"
#include "asterfort/xcalfev_wrap.h"
#include "asterfort/xkamat.h"
#include "asterfort/xnbddl.h"
#include "jeveux.h"
#include "MeshTypes_type.h"
!
    character(len=8), intent(in) :: typmod(2)
    type(Material_Para), intent(inout) :: materPara
    integer(kind=8) :: nnop, nfiss, ddlc, ddlm
    integer(kind=8) :: idecpg, jvGeom, nnops
    integer(kind=8) :: ndim, nfe, nfh, npg, heavn(nnop, 5)
    integer(kind=8) :: jstno
    real(kind=8) :: basloc(3*ndim*nnop), coorse(*), he(nfiss)
    real(kind=8) :: lsn(nnop), lst(nnop)
    real(kind=8) :: matuu(*)
    character(len=8) :: elrefp, elrese, fami_se
!
! --------------------------------------------------------------------------------------------------
!
!     BUT:  CALCUL  DES OPTIONS RIGI_MECA_TANG, RAPH_MECA ET FULL_MECA
!           EN HYPER-ELASTICITE AVEC X-FEM EN 2D ET EN 3D
!
! --------------------------------------------------------------------------------------------------
!
! IN  ELREFP  : ÉLÉMENT DE RÉFÉRENCE PARENT
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  COORSE  : COORDONNÉES DES SOMMETS DU SOUS-ÉLÉMENT
! IN  IGEOM   : COORDONNÉES DES NOEUDS DE L'ÉLÉMENT PARENT
! IN  HE      : VALEUR DE LA FONCTION HEAVISIDE SUR LE SOUS-ÉLT
! IN  NFH     : NOMBRE DE DDL HEAVYSIDE (PAR NOEUD)
! IN  DDLC    : NOMBRE DE DDL DE CONTACT (PAR NOEUD)
! IN  DDLM    : NOMBRE DE DDL PAR NOEUD MILIEU
! IN  NNOPS   : NOMBRE DE NOEUDS SOMMET ELEMENTS PARENT
! IN  NFE     : NOMBRE DE FONCTIONS SINGULIÈRES D'ENRICHISSEMENT
! IN  BASLOC  : BASE LOCALE AU FOND DE FISSURE
! IN  NNOP    : NOMBRE DE NOEUDS DE L'ELEMENT PARENT
! IN  NPG     : NOMBRE DE POINTS DE GAUSS DU SOUS-ÉLÉMENT
! IN  LGPG    : "LONGUEUR" DES VARIABLES INTERNES POUR 1 POINT DE GAUSS
!               CETTE LONGUEUR EST UN MAJORANT DU NBRE REEL DE VAR. INT.
! IN  IDEPL   : ADRESSE DU DEPLACEMENT A PARTIR DE LA CONF DE REF
! IN  LSN     : VALEUR DE LA LEVEL SET NORMALE AUX NOEUDS PARENTS
! IN  LST     : VALEUR DE LA LEVEL SET TANGENTE AUX NOEUDS PARENTS
! IN  IDECPG  : POSITION DANS LA FAMILLE 'XFEM' DU 1ER POINT DE GAUSS
!               DU SOUS ELEMENT COURRANT (EN FAIT IDECPG+1)
!
! OUT SIG     : CONTRAINTES DE CAUCHY (RAPH_MECA ET FULL_MECA)
! OUT VI      : VARIABLES INTERNES    (RAPH_MECA ET FULL_MECA)
! OUT MATUU   : MATRICE DE RIGIDITE PROFIL (RIGI_MECA_TANG ET FULL_MECA)
! OUT IVECTU  : VECTEUR FORCES NODALES (RAPH_MECA ET FULL_MECA)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: ksp = 1
    integer(kind=8) :: kpg, i, ig, n, nn, m, mn, j, j1, kl, l, kkd, ipg
    integer(kind=8) :: ddld, ddls, nno, npgbis, cpt, ndimb, dec(nnop)
    integer(kind=8) :: idfde, ipoids, ivf, hea_se
    integer(kind=8) :: singu, alp, ii, jj
    real(kind=8) :: dsidep(6, 6), eps(6)
    real(kind=8) :: tmp2, sigp(6, 3*(1+nfh+nfe*ndim))
    real(kind=8) :: xg(ndim), xe(ndim), ff(nnop), jac
    real(kind=8) :: dfdi(nnop, ndim), f(3, 3)
    real(kind=8) :: def(6, ndim*(1+nfh+ndim), nnop)
    real(kind=8) :: r
    real(kind=8) :: fk(MT_NNOMAX, 3, 3), dkdgl(MT_NNOMAX, 3, 3, 3), ka, mu
    integer(kind=8) :: nbsig
    real(kind=8) :: d(36), instan
    aster_logical :: axi, cplan
    real(kind=8), parameter :: rac2 = 1.4142135623731d0
!
! --------------------------------------------------------------------------------------------------
!

!
!     ATTENTION, EN 3D, ZR(IDEPL) ET ZR(VECTU) SONT DIMENSIONNÉS DE
!     TELLE SORTE QU'ILS NE PRENNENT PAS EN COMPTE LES DDL SUR LES
!     NOEUDS MILIEU
!
!     NOMBRE DE DDL DE DEPLACEMENT À CHAQUE NOEUD
    call xnbddl(ndim, nfh, nfe, ddlc, ddld, ddls, singu)
    call xkamat(materPara%jvMaterCode, ndim, axi, ka, mu)

! - Type of finite element
    axi = typmod(1) .eq. 'AXIS'
    cplan = typmod(1) .eq. 'C_PLAN'

! - ADRESSE DES COORD DU SOUS ELT EN QUESTION
    fami_se = 'XINT'
    if (nfe .gt. 0) then
        if (ndim .eq. 3 .and. &
            (count(zi((jstno-1+1):(jstno-1+nnop)) .eq. 2)+ &
             count(zi((jstno-1+1):(jstno-1+nnop)) .eq. 0)) .eq. nnop) then
            fami_se = 'XGEO'
        end if
    end if

! - Get element parameters
    call elrefe_info(elrefe=elrese, fami=fami_se, ndim=ndimb, nno=nno, &
                     npg=npgbis, jpoids=ipoids, jvf=ivf, jdfde=idfde)
    ASSERT(npg .eq. npgbis .and. ndim .eq. ndimb)

! - DECALAGES CALCULES EN AMONT: PERF
    do n = 1, nnop
        call indent(n, ddls, ddlm, nnops, dec(n))
    end do

! - CALCUL DE L IDENTIFIANT DU SS ELEMENT
    hea_se = xcalc_code(nfiss, he_real=[he])

! - Loop on Gauss points
    do kpg = 1, npg
! ----- Initializations of material parameters on current integration point
        call initParaPoin(kpg, ksp, materPara)

!       COORDONNÉES DU PT DE GAUSS DANS LE REPÈRE RÉEL : XG
        xg = 0.d0
        do i = 1, ndim
            do n = 1, nno
                xg(i) = xg(i)+zr(ivf-1+nno*(kpg-1)+n)*coorse(ndim*(n-1)+i)
            end do
        end do
!
!       COORDONNÉES DU POINT DE GAUSS DANS L'ÉLÉMENT DE RÉF PARENT : XE
!       ET CALCUL DE FF ET DFDI
        call reeref(elrefp, nnop, zr(jvGeom), xg, ndim, &
                    xe, ff, dfdi=dfdi)
!
!       FONCTION D'ENRICHISSEMENT AU POINT DE GAUSS ET LEURS DÉRIVÉES
        if (singu .gt. 0) then
            call xcalfev_wrap(ndim, nnop, basloc, zi(jstno), he(1), &
                              lsn, lst, zr(jvGeom), ka, mu, ff, fk, dfdi, dkdgl, &
                              elref=elrefp, kstop='C')
        end if
! -     CALCUL DE LA DISTANCE A L'AXE (AXISYMETRIQUE)
        if (axi) then
            r = 0.d0
            do n = 1, nnop
                r = r+ff(n)*zr(jvGeom-1+2*(n-1)+1)
            end do
!
            ASSERT(r .gt. 0d0)
!           ATTENTION : LE POIDS N'EST PAS X R
!           CE SERA FAIT PLUS TARD AVEC JAC = JAC X R
        end if
!
!       COORDONNÉES DU POINT DE GAUSS DANS L'ÉLÉMENT DE RÉF PARENT : XE
!       ET CALCUL DE FF, DFDI, ET EPS

        f = 0.d0
        eps = 0.d0
        do i = 1, 3
            f(i, i) = 1.d0
        end do

!
! - CALCUL DES ELEMENTS GEOMETRIQUES
!
!       CALCUL DES PRODUITS SYMETR. DE F PAR N,
        def(:, :, :) = 0.d0
        do n = 1, nnop
            cpt = 0
!         FONCTIONS DE FORME CLASSIQUES
            do i = 1, ndim
                cpt = cpt+1
                def(1, i, n) = f(i, 1)*dfdi(n, 1)
                def(2, i, n) = f(i, 2)*dfdi(n, 2)
                def(3, i, n) = 0.d0
                def(4, i, n) = (f(i, 1)*dfdi(n, 2)+f(i, 2)*dfdi(n, 1))/rac2
                if (ndim .eq. 3) then
                    def(3, i, n) = f(i, 3)*dfdi(n, 3)
                    def(5, i, n) = (f(i, 1)*dfdi(n, 3)+f(i, 3)*dfdi(n, 1))/rac2
                    def(6, i, n) = (f(i, 2)*dfdi(n, 3)+f(i, 3)*dfdi(n, 2))/rac2
                end if
            end do
!         TERME DE CORRECTION (3,3) AXI QUI PORTE EN FAIT SUR LE DDL 1
            if (axi) then
                def(3, 1, n) = f(3, 3)*ff(n)/r
            end if
!         ENRICHISSEMENT PAR HEAVISIDE
            do ig = 1, nfh
                do i = 1, ndim
                    cpt = cpt+1
                    do m = 1, 2*ndim
                        def(m, cpt, n) = &
                            def(m, i, n)*xcalc_heav(heavn(n, ig), hea_se, heavn(n, 5))
                    end do
                    if (ndim .eq. 2) then
                        def(3, cpt, n) = 0.d0
                    end if
                end do
!               TERME DE CORRECTION (3,3) A PORTER SUR LE DDL 1+NDIM*IG
                if (axi) then
                    def(3, 1+ndim*ig, n) = &
                        f(3, 3)*ff(n)/r*xcalc_heav(heavn(n, ig), hea_se, heavn(n, 5))
                end if
            end do
!         ENRICHISSEMENT PAR LES NFE FONTIONS SINGULIÈRES
            do alp = 1, ndim*nfe
                do i = 1, ndim
                    cpt = cpt+1
                    def(1, cpt, n) = f(i, 1)*dkdgl(n, alp, i, 1)
                    def(2, cpt, n) = f(i, 2)*dkdgl(n, alp, i, 2)
                    def(3, cpt, n) = 0.d0
                    def(4, cpt, n) = (f(i, 1)*dkdgl(n, alp, i, 2)+ &
                                      f(i, 2)*dkdgl(n, alp, i, 1))/rac2
                    if (ndim .eq. 3) then
                        def(3, cpt, n) = f(i, 3)*dkdgl(n, alp, i, 3)
                        def(5, cpt, n) = (f(i, 1)*dkdgl(n, alp, i, 3)+ &
                                          f(i, 3)*dkdgl(n, alp, i, 1))/rac2
                        def(6, cpt, n) = (f(i, 3)*dkdgl(n, alp, i, 2)+ &
                                          f(i, 2)*dkdgl(n, alp, i, 3))/rac2
                    end if
                end do
            end do
!   TERME DE CORRECTION (3,3) AXI PORTE SUR LE DDL 1+NDIM*(NFH+ALP)
!      EN AXI: ON PROJETTE L ENRICHISSEMENT VECTORIEL SUIVANT X
            if (axi) then
                do alp = 1, ndim*nfe
                    def(3, 1+ndim*(nfh+alp), n) = f(3, 3)*fk(n, alp, 1)/r
                end do
            end if
            ASSERT(cpt .eq. ddld)
        end do
!       CALCULER LE JACOBIEN DE LA TRANSFO SSTET->SSTET REF
!       AVEC LES COORDONNEES DU SOUS-ELEMENT
        if (ndim .eq. 2) then
            call dfdm2d(nno, kpg, ipoids, idfde, coorse, jac)
        else if (ndim .eq. 3) then
            call dfdm3d(nno, kpg, ipoids, idfde, coorse, jac)
        end if
        if (axi) then
            jac = jac*r
        end if
!
! - CALCUL DE LA MATRICE DE RIGIDITE POUR L'OPTION RIGI_MECA
!

!           Calcul du tenseur de Hooke [D]
!              {sxx, syy, szz, sxy, sxz, syz}^T = [Ð]{exx, eyy, ezz, 2*exy, 2*exz, 2*eyz}^T
        ipg = idecpg+kpg
        instan = r8vide()
        nbsig = 2*ndim
!
        call dmatmc(materPara, "+", instan, &
                    nbsig, d)
!
!           Calcul du tenseur de comportement tangent [D']
!              {sxx, syy, szz, sqrt(2)*sxy, sqrt(2)*sxz, sqrt(2)*syz}^T =
!                  [D']{exx, eyy, ezz, sqrt(2)*exy, sqrt(2)*exz, sqrt(2)*eyz}^T
        do j = 1, nbsig
            do i = 1, nbsig
                dsidep(i, j) = d((j-1)*nbsig+i)
            end do
        end do
        do j = 4, nbsig
            do i = 4, nbsig
                dsidep(i, j) = 2.d0*d((j-1)*nbsig+i)
            end do
        end do
!
!           Calcul de la matrice de rigidité : [K] = [B]^T[D'][B]
!
        do n = 1, nnop
            nn = dec(n)
            do i = 1, ddld
                do kl = 1, 2*ndim
                    sigp(kl, i) = 0.d0
                    do l = 1, 2*ndim
                        sigp(kl, i) = sigp(kl, i)+def(l, i, n)*dsidep(l, kl)
                    end do
                end do
            end do
!
            do m = 1, n
                mn = dec(m)
                do i = 1, ddld
                    ii = iimatu(i, ndim, nfh, nfe)
                    kkd = (nn+ii-1)*(nn+ii)/2
                    if (m .eq. n) then
                        j1 = ii
                    else
                        j1 = ddld
                    end if
!
                    do j = 1, ddld
!
                        jj = iimatu(j, ndim, nfh, nfe)
!
!                 RIGIDITE ELASTIQUE
                        tmp2 = 0.d0
                        do l = 1, 2*ndim
                            tmp2 = tmp2+sigp(l, i)*def(l, j, m)
                        end do
!
!                 STOCKAGE EN TENANT COMPTE DE LA SYMETRIE
                        if (jj .le. j1) then
                            matuu(kkd+mn+jj) = matuu(kkd+mn+jj)+tmp2*jac
                        end if
!
                    end do
                end do
            end do
        end do

    end do
!
end subroutine
