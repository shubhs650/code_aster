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
subroutine xxnmpl(BEHInteg, &
                  option, typmod, &
                  compor, carcri, &
                  elrefp, elrese, ndim, coorse, jvGeom, &
                  he, nfh, ddlc, ddlm, nfe, &
                  instam, instap, ideplp, sigm, vip, &
                  basloc, nnop, npg, &
                  lgpg, idepl, &
                  lsn, lst, idecpg, sig, vi, &
                  matuu, ivectu, codret, nfiss, heavn, jstno, &
                  lMatr, lVect, lSigm)
!
    use Behaviour_type
    use Behaviour_module
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/dfdm2d.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/iimatu.h"
#include "asterfort/indent.h"
#include "asterfort/nmcomp.h"
#include "asterfort/reeref.h"
#include "asterfort/xcalc_code.h"
#include "asterfort/xcalc_heav.h"
#include "asterfort/xcalfev_wrap.h"
#include "asterfort/xcinem.h"
#include "asterfort/xkamat.h"
#include "asterfort/xnbddl.h"
#include "jeveux.h"
#include "MeshTypes_type.h"
!
    type(Behaviour_Integ), intent(inout) :: BEHInteg
    character(len=8), intent(in) :: typmod(2)
    character(len=16), intent(in) :: option, compor(COMPOR_SIZE)
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    integer(kind=8) :: ndim, jvGeom, lgpg, codret, nnop, npg
    integer(kind=8) :: nfh, ddlc, ddlm, nfe, idepl, ivectu, ideplp
    integer(kind=8) :: nfiss, heavn(nnop, 5), idecpg
    integer(kind=8) :: jstno
    character(len=8) :: elrefp, elrese
    real(kind=8) :: basloc(3*ndim*nnop), he(nfiss)
    real(kind=8) :: lsn(nnop), lst(nnop), coorse(*)
    real(kind=8) :: vi(lgpg, npg), vip(lgpg, npg), sig(2*ndim, npg), matuu(*)
    real(kind=8) :: instam, instap, sigm(2*ndim, npg), sign(6)
    aster_logical, intent(in) :: lMatr, lVect, lSigm
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
! IN  NFH     : NOMBRE DE FONCTIONS HEAVYSIDE
! IN  DDLC    : NOMBRE DE DDL DE CONTACT (PAR NOEUD)
! IN  DDLM    : NOMBRE DE DDL PAR NOEUD MILIEU (EN 2D)
! IN  NFE     : NOMBRE DE FONCTIONS SINGULIÈRES D'ENRICHISSEMENT
! IN  BASLOC  : BASE LOCALE AU FOND DE FISSURE AUX NOEUDS
! IN  NNOP    : NOMBRE DE NOEUDS DE L'ELEMENT PARENT
! IN  NPG     : NOMBRE DE POINTS DE GAUSS DU SOUS-ÉLÉMENT
! IN  LGPG    : "LONGUEUR" DES VARIABLES INTERNES POUR 1 POINT DE GAUSS
!               CETTE LONGUEUR EST UN MAJORANT DU NBRE REEL DE VAR. INT.
! IN  IDEPL   : ADRESSE DU DEPLACEMENT A PARTIR DE LA CONF DE REF
! IN  LSN     : VALEUR DE LA LEVEL SET NORMALE AUX NOEUDS PARENTS
! IN  LST     : VALEUR DE LA LEVEL SET TANGENTE AUX NOEUDS PARENTS
! OUT SIG     : CONTRAINTES DE CAUCHY (RAPH_MECA ET FULL_MECA)
! OUT VI      : VARIABLES INTERNES    (RAPH_MECA ET FULL_MECA)
! OUT MATUU   : MATRICE DE RIGIDITE PROFIL (RIGI_MECA_TANG ET FULL_MECA)
! OUT IVECTU  : VECTEUR FORCES NODALES (RAPH_MECA ET FULL_MECA)
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: multComp = " "
    integer(kind=8), parameter :: ksp = 1
    integer(kind=8) :: i, ig, j, j1, kkd, kl, kpg, l, m, n, nn, mn
    integer(kind=8) :: ddls, ddld, cpt, idfde, ipoids, ivf, dec(nnop)
    integer(kind=8) :: ndimb, nno, nnops, npgbis, hea_se
    integer(kind=8) :: singu, alp, ii, jj
    real(kind=8) :: dsidep(6, 6), f(3, 3), deps(6), sigma(6)
    real(kind=8) :: eps(6), sigp(6)
    real(kind=8) :: tmp2
    real(kind=8) :: xg(ndim), xe(ndim), ff(nnop), jac
    real(kind=8) :: rbid33(3, 3)
    real(kind=8) :: dfdi(nnop, ndim)
    real(kind=8) :: def(6, nnop, ndim*(1+nfh+nfe*ndim)), r
    real(kind=8) :: fk(MT_NNOMAX, 3, 3), dkdgl(MT_NNOMAX, 3, 3, 3), ka, mu
    aster_logical :: axi, cplan
    real(kind=8), parameter :: rac2 = 1.4142135623731d0
!
! --------------------------------------------------------------------------------------------------
!

!     ATTENTION, EN 3D, ZR(IDEPL) ET ZR(VECTU) SONT DIMENSIONNÉS DE
!     TELLE SORTE QU'ILS NE PRENNENT PAS EN COMPTE LES DDL SUR LES
!     NOEUDS MILIEU
!
!     NOMBRE DE DDL DE DEPLACEMENT À CHAQUE NOEUD
    call xnbddl(ndim, nfh, nfe, ddlc, ddld, ddls, singu)
!     RECUPERATION DU NOMBRE DE NOEUDS SOMMETS DE L'ELEMENT PARENT
    call elrefe_info(fami='RIGI', nnos=nnops)
!
! - Type of finite element
!
    axi = typmod(1) .eq. 'AXIS'
    cplan = typmod(1) .eq. 'C_PLAN'

! - Get element parameters
    call elrefe_info(elrefe=elrese, fami='XINT', ndim=ndimb, nno=nno, &
                     npg=npgbis, jpoids=ipoids, jvf=ivf, jdfde=idfde)
    ASSERT(npg .eq. npgbis .and. ndim .eq. ndimb)

! - Prepare external state variables
    call behaviourPrepESVAGeom(nno, npg, ndim, &
                               ipoids, ivf, idfde, &
                               zr(jvGeom), BEHInteg)

! - DECALAGES CALCULES EN AMONT: PERF
    do n = 1, nnop
        call indent(n, ddls, ddlm, nnops, dec(n))
    end do
!
! - CALCUL DE L IDENTIFIANT DU SS ELEMENT
!
    hea_se = xcalc_code(nfiss, he_real=[he])

! - Loop on Gauss points
    do kpg = 1, npg
!
!       COORDONNÉES DU PT DE GAUSS DANS LE REPÈRE RÉEL : XG
        xg = 0.d0
        do i = 1, ndim
            do n = 1, nno
                xg(i) = xg(i)+zr(ivf-1+nno*(kpg-1)+n)*coorse(ndim*(n-1)+i)
            end do
        end do
!
!       COORDONNÉES DU POINT DE GAUSS DANS L'ÉLÉMENT DE RÉF PARENT : XE
!       CALCUL DE FF ET DFDI
        call reeref(elrefp, nnop, zr(jvGeom), xg, ndim, &
                    xe, ff, dfdi=dfdi)
!
!       FONCTION D'ENRICHISSEMENT AU POINT DE GAUSS ET LEURS DÉRIVÉES
        if (singu .gt. 0) then
            call xkamat(BEHInteg%materPara%jvMaterCode, ndim, axi, ka, mu)
            call xcalfev_wrap(ndim, nnop, basloc, zi(jstno), he(1), &
                              lsn, lst, zr(jvGeom), ka, mu, ff, fk, dfdi, dkdgl)
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
! -     CALCUL DE DEPS
        call xcinem(axi, jvGeom, nnop, nnops, ideplp, &
                    ndim, he, &
                    nfiss, nfh, singu, ddls, ddlm, &
                    fk, dkdgl, ff, dfdi, f, &
                    deps, rbid33, heavn)
! -     CALCUL DE EPS
        call xcinem(axi, jvGeom, nnop, nnops, idepl, &
                    ndim, he, &
                    nfiss, nfh, singu, ddls, ddlm, &
                    fk, dkdgl, ff, dfdi, f, &
                    eps, rbid33, heavn)
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
                def(1, n, i) = f(i, 1)*dfdi(n, 1)
                def(2, n, i) = f(i, 2)*dfdi(n, 2)
                def(3, n, i) = 0.d0
                def(4, n, i) = (f(i, 1)*dfdi(n, 2)+f(i, 2)*dfdi(n, 1))/rac2
                if (ndim .eq. 3) then
                    def(3, n, i) = f(i, 3)*dfdi(n, 3)
                    def(5, n, i) = (f(i, 1)*dfdi(n, 3)+f(i, 3)*dfdi(n, 1))/rac2
                    def(6, n, i) = (f(i, 2)*dfdi(n, 3)+f(i, 3)*dfdi(n, 2))/rac2
                end if
            end do
!         TERME DE CORRECTION (3,3) AXI QUI PORTE EN FAIT SUR LE DDL 1
            if (axi) then
                def(3, n, 1) = f(3, 3)*ff(n)/r
            end if
!         ENRICHISSEMENT PAR HEAVISIDE
            do ig = 1, nfh
                do i = 1, ndim
                    cpt = cpt+1
                    do m = 1, 2*ndim
                        def(m, n, cpt) = def(m, n, i)*xcalc_heav(heavn(n, ig), hea_se, heavn(n, 5))
                    end do
                    if (ndim .eq. 2) then
                        def(3, n, cpt) = 0.d0
                    end if
                end do
!               TERME DE CORRECTION (3,3) A PORTER SUR LE DDL 1+NDIM*IG
                if (axi) then
                    def(3, n, (1+ndim*ig)) = &
                        f(3, 3)*ff(n)/r*xcalc_heav(heavn(n, ig), hea_se, heavn(n, 5))
                end if
            end do
!         ENRICHISSEMENT PAR LES NFE FONTIONS SINGULIÈRES
            do ig = 1, singu
                do alp = 1, ndim
                    do i = 1, ndim
                        cpt = cpt+1
                        def(1, n, cpt) = f(i, 1)*dkdgl(n, alp, i, 1)
                        def(2, n, cpt) = f(i, 2)*dkdgl(n, alp, i, 2)
                        def(3, n, cpt) = 0.d0
                        def(4, n, cpt) = (f(i, 1)*dkdgl(n, alp, i, 2)+ &
                                          f(i, 2)*dkdgl(n, alp, i, 1))/rac2
                        if (ndim .eq. 3) then
                            def(3, n, cpt) = f(i, 3)*dkdgl(n, alp, i, 3)
                            def(5, n, cpt) = (f(i, 1)*dkdgl(n, alp, i, 3)+ &
                                              f(i, 3)*dkdgl(n, alp, i, 1))/rac2
                            def(6, n, cpt) = (f(i, 3)*dkdgl(n, alp, i, 2)+ &
                                              f(i, 2)*dkdgl(n, alp, i, 3))/rac2
                        end if
                    end do
                end do
!   TERME DE CORRECTION (3,3) AXI PORTE SUR LE DDL 1+NDIM*(NFH+ALP)
!      EN AXI: ON PROJETTE L ENRICHISSEMENT VECTORIEL SUIVANT X
                if (axi) then
                    do alp = 1, ndim
                        def(3, n, (1+ndim*(nfh+alp))) = f(3, 3)*fk(n, alp, 1)/r
                    end do
                end if
            end do
            ASSERT(cpt .eq. ddld)
        end do
!       CALCULER LE JACOBIEN DE LA TRANSFO SSTET->SSTET REF
!       AVEC LES COORDONNEES DU SOUS-ELEMENT
        if (ndim .eq. 2) then
            call dfdm2d(nno, kpg, ipoids, idfde, coorse, &
                        jac)
        else if (ndim .eq. 3) then
            call dfdm3d(nno, kpg, ipoids, idfde, coorse, &
                        jac)
        end if
        if (axi) then
            jac = jac*r
        end if
! ----- Preparation for behaviour
        do m = 1, 3
            sign(m) = sigm(m, kpg)
        end do
        do m = 4, 2*ndim
            sign(m) = sigm(m, kpg)*rac2
        end do
        sigma = 0.d0

! ----- Set main parameters for behaviour (on point)
        call behaviourSetParaPoin(idecpg+kpg, ksp, BEHInteg)

! ----- Integrate
        call nmcomp(BEHInteg, &
                    ndim, option, typmod, &
                    instam, instap, &
                    compor, carcri, multComp, &
                    6, eps, deps, &
                    6, sign, &
                    vi(1, kpg), &
                    sigma, vip(1, kpg), &
                    36, dsidep, codret)

! ----- Rigidity matrix
        if (lMatr) then
            do n = 1, nnop
                nn = dec(n)
                do i = 1, ddld
                    ii = iimatu(i, ndim, nfh, nfe)
                    kkd = (nn+ii-1)*(nn+ii)/2
                    do kl = 1, 2*ndim
                        sigp(kl) = 0.d0
                        do l = 1, 2*ndim
                            sigp(kl) = sigp(kl)+def(l, n, i)*dsidep(l, kl)
                        end do
                    end do
                    do j = 1, ddld
                        jj = iimatu(j, ndim, nfh, nfe)
                        do m = 1, n
                            mn = dec(m)
!
                            if (m .eq. n) then
                                j1 = ii
                            else
                                j1 = ddld
                            end if
!
!                 RIGIDITE ELASTIQUE
                            tmp2 = 0.d0
                            do l = 1, 2*ndim
                                tmp2 = tmp2+sigp(l)*def(l, m, j)
                            end do
!
!                 STOCKAGE EN TENANT COMPTE DE LA SYMETRIE
                            if (jj .le. j1) then
                                matuu(kkd+mn+jj) = matuu(kkd+mn+jj)+(tmp2)*jac
                            end if
!
                        end do
                    end do
                end do
            end do
        end if
! ----- Internal forces
        if (lVect) then
            do n = 1, nnop
                nn = dec(n)
                do i = 1, ddld
                    ii = iimatu(i, ndim, nfh, nfe)
                    do l = 1, 2*ndim
                        zr(ivectu-1+nn+ii) = zr(ivectu-1+nn+ii)+def(l, n, i)*sigma(l)*jac
                    end do
                end do
            end do
        end if
! ----- Stress
        if (lSigm) then
            do l = 1, 3
                sig(l, kpg) = sigma(l)
            end do
            do l = 4, 2*ndim
                sig(l, kpg) = sigma(l)/rac2
            end do
        end if
!
    end do
!
end subroutine
