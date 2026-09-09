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
subroutine xsifl1(elrefp, angl, basloc, coeff, coeff3, &
                  ddlm, ddls, dfdi, ff, he, &
                  heavn, idepl, igthet, ipref, ipres, &
                  ithet, jac, jlsn, jlst, jstno, &
                  ka, mu, nd, ndim, nfh, &
                  nnop, nnops, itemps, nompar, option, &
                  singu, xg, igeom)
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/coor_cyl.h"
#include "asterfort/fointe.h"
#include "asterfort/indent.h"
#include "asterfort/lteatt.h"
#include "asterfort/utmess.h"
#include "asterfort/xcalc_code.h"
#include "asterfort/xcalc_heav.h"
#include "MeshTypes_type.h"
#include "asterfort/xcalfev_wrap.h"
#include "asterfort/xdeffk.h"
!
! Calcul de G avec forces de pression XFEM sur les levres
!   de la fissure
!
    character(len=8) :: elrefp
    integer(kind=8) :: nnop, ndim, heavn(nnop, 5)
    integer(kind=8) :: jstno, jlsn
    real(kind=8) :: angl(2), basloc(9*nnop), cisa, coeff, coeff3
    integer(kind=8) :: cpt, ddlm, ddls
    real(kind=8) :: depla(3), dfdi(nnop, ndim), dfor(3), divt
    real(kind=8) :: dtdm(3, 3), ff(MT_NNOMAX3D)
    real(kind=8) :: forrep(3, 2), g, he(2)
    integer(kind=8) :: i, idepl, ier, igthet, ilev, indi, ino, ig, hea_fa(2)
    integer(kind=8) :: ipref, ipres, ithet, j
    real(kind=8) :: jac
    integer(kind=8) :: jlst
    real(kind=8) :: jm, var(ndim+1)
    real(kind=8) :: k1, k2, k3, ka, mu, nd(3)
    integer(kind=8) :: nfh, nnops, itemps
    character(len=8) :: nompar(4)
    character(len=16) :: option
    real(kind=8) :: p(3, 3), pres, invp(3, 3)
    real(kind=8) :: pres_test, cisa_test, r8pre
    integer(kind=8) :: singu
    real(kind=8) :: theta(3), u1(3), u2(3), u3(3)
    real(kind=8) :: xg(3), r
    aster_logical :: axi, l_pres_var, l_cisa_var, l_not_zero
    real(kind=8) :: fk(MT_NNOMAX3D, 3, 3), fkpo(3, 3)
    integer(kind=8) :: alp, igeom
    real(kind=8) :: rg, tg
!
!
! CALCUL DE L IDENTIFIANT DE FACETTES MAITRE/ESCLAVE
    do ilev = 1, 2
        hea_fa(ilev) = xcalc_code(1, he_real=[he(ilev)])
    end do
!
!     -----------------------------------------------
!     1) CALCUL DES FORCES SUIVANT LES OPTIONS
!     -----------------------------------------------
!
    forrep(:, :) = 0.d0
!
    if ((option .eq. 'CALC_K_G_XFEM') .or. (option .eq. 'CALC_G_XFEM')) then
!
!         CALCUL DE LA PRESSION AUX POINTS DE GAUSS
        pres = 0.d0
        cisa = 0.d0
        do ino = 1, nnop
            if (ndim .eq. 3) pres = pres+zr(ipres-1+ino)*ff(ino)
            if (ndim .eq. 2) then
                pres = pres+zr(ipres-1+2*(ino-1)+1)*ff(ino)
                cisa = cisa+zr(ipres-1+2*(ino-1)+2)*ff(ino)
            end if
        end do
        do j = 1, ndim
            forrep(j, 1) = -pres*nd(j)
            forrep(j, 2) = -pres*(-nd(j))
        end do
        if (ndim .eq. 2) then
            forrep(1, 1) = forrep(1, 1)-cisa*nd(2)
            forrep(2, 1) = forrep(2, 1)+cisa*nd(1)
            forrep(1, 2) = forrep(1, 2)-cisa*(-nd(2))
            forrep(2, 2) = forrep(2, 2)+cisa*(-nd(1))
        end if
!
    else if ((option .eq. 'CALC_K_G_XFEM_F') .or. (option .eq. 'CALC_G_XFEM_F')) then
!
!         VALEUR DE LA PRESSION
        var(:) = 0.d0
        do j = 1, ndim
            var(j) = xg(j)
        end do
        var(ndim+1) = zr(itemps)
        call fointe('FM', zk8(ipref), ndim+1, nompar, var, &
                    pres, ier)
        if (ndim .eq. 2) call fointe('FM', zk8(ipref+1), ndim+1, nompar, var, &
                                     cisa, ier)
        do j = 1, ndim
            forrep(j, 1) = -pres*nd(j)
            forrep(j, 2) = -pres*(-nd(j))
        end do
        if (ndim .eq. 2) then
            forrep(1, 1) = forrep(1, 1)-cisa*nd(2)
            forrep(2, 1) = forrep(2, 1)+cisa*nd(1)
            forrep(1, 2) = forrep(1, 2)-cisa*(-nd(2))
            forrep(2, 2) = forrep(2, 2)+cisa*(-nd(1))
        end if
    else
        call utmess('F', 'XFEM_15')
    end if
!
!     -----------------------------------
!     2) CALCUL DE THETA ET DE DIV(THETA)
!     -----------------------------------
    divt = 0.d0
    dtdm(:, :) = 0.d0
!
    do i = 1, ndim
        theta(i) = 0.d0
        do ino = 1, nnop
            theta(i) = theta(i)+ff(ino)*zr(ithet-1+ndim*(ino-1)+i)
        end do
!
        do j = 1, ndim
            do ino = 1, nnop
                dtdm(i, j) = dtdm(i, j)+zr(ithet-1+ndim*(ino-1)+i)*dfdi(ino, j)
            end do
        end do
!
        divt = divt+dtdm(i, i)
!
    end do
!
    axi = lteatt('AXIS', 'OUI')
    if (axi) then
        r = 0.d0
        do ino = 1, nnop
            r = r+ff(ino)*zr(igeom-1+2*(ino-1)+1)
        end do
        ASSERT(r .gt. 0d0)
        divt = divt+theta(1)/r
        jac = jac*r
    end if
!
!
!     BOUCLE SUR LES DEUX LEVRES
    do ilev = 1, 2
!
!
!       FONCTION D'ENRICHISSEMENT AU POINT DE GAUSS
        if (singu .gt. 0) then
            if (he(ilev) .gt. 0) then
                call xcalfev_wrap(ndim, nnop, basloc, zi(jstno), he(ilev), &
                                  zr(jlsn), zr(jlst), zr(igeom), ka, mu, &
                                  ff, fk, face='MAIT', elref=elrefp, kstop='C')
            else
                call xcalfev_wrap(ndim, nnop, basloc, zi(jstno), he(ilev), &
                                  zr(jlsn), zr(jlst), zr(igeom), ka, mu, &
                                  ff, fk, face='ESCL', elref=elrefp, kstop='C')
            end if
        end if
!       CALCUL DES COORDONNEES CYLINDRIQUES
        call coor_cyl(ndim, nnop, basloc, zr(igeom), ff, &
                      p, invp, rg, tg, l_not_zero)
!       ---------------------------------------------
!       3) CALCUL DU DEPLACEMENT
!       ---------------------------------------------
        depla(:) = 0.d0
        do ino = 1, nnop
            call indent(ino, ddls, ddlm, nnops, indi)
            cpt = 0
!         DDLS CLASSIQUES
            do i = 1, ndim
                cpt = cpt+1
                depla(i) = depla(i)+ff(ino)*zr(idepl-1+indi+cpt)
            end do
!         DDLS HEAVISIDE
            do i = 1, ndim
                do ig = 1, nfh
                    cpt = cpt+1
                    depla(i) = depla(i)+xcalc_heav(heavn(ino, ig), hea_fa(ilev), heavn(ino, 5))*f&
                               &f(ino)*zr(idepl-1+indi+cpt)
                end do
            end do
!         DDL ENRICHIS EN FOND DE FISSURE
            do alp = 1, singu*ndim
                cpt = cpt+1
                do i = 1, ndim
                    depla(i) = depla(i)+fk(ino, alp, i)*zr(idepl-1+indi+cpt)
                end do
            end do
        end do
!
!       --------------------------------
!       4) CALCUL DES CHAMPS AUXILIAIRES
!       --------------------------------
!
        if (option(1:8) .eq. 'CALC_K_G') then
!
! --------- champs singuliers
            call xdeffk(ka, mu, rg, angl(ilev), ndim, &
                        fkpo(1:ndim, 1:ndim))
!
!         CHAMPS AUXILIARES DANS LA BASE GLOBALE : U1,U2,U3
            u1(:) = 0.d0
            u2(:) = 0.d0
            u3(:) = 0.d0
            do i = 1, ndim
                do j = 1, ndim
                    u1(i) = u1(i)+p(i, j)*fkpo(1, j)
                    u2(i) = u2(i)+p(i, j)*fkpo(2, j)
                    if (ndim .eq. 3) u3(i) = u3(i)+p(i, j)*fkpo(3, j)
                end do
            end do
        end if
!
!       -----------------------------------------
!       5) CALCUL DE 'DFOR' =  D(PRES)/DI . THETA
!       -----------------------------------------
        dfor(:) = 0.d0
!
!       D(PRES)/DI n'etait pas correctement calcule dans cette routine.
!       issue24174 supprime le calcul de cette quantite (DFOR reste nul)
!       et interdit toute autre chose qu'un chargement constant.
        if ((option .eq. 'CALC_K_G_XFEM') .or. (option .eq. 'CALC_G_XFEM')) then
!
!           Tester le nom de l'option (CALC_*G ou CALC_*G_F) ne suffit
!           pas pour detecter le caractere constant du chargement si on
!           a un evol_char. Le test ci-dessous est plus robuste.
            r8pre = r8prem()
!           en 2D on a les composantes PRES, CISA
            if (ndim .eq. 2) then
                pres_test = zr(ipres-1+2*(1-1)+1)
                cisa_test = zr(ipres-1+2*(1-1)+2)
                do ino = 1, nnop
                    l_pres_var = abs(zr(ipres-1+2*(ino-1)+1)-pres_test) .ge. r8pre
                    l_cisa_var = abs(zr(ipres-1+2*(ino-1)+2)-cisa_test) .ge. r8pre
                    if (l_pres_var .or. l_cisa_var) then
                        ASSERT(.false.)
                    end if
                end do
            end if
!           en 3D on a uniquement la composante PRES
            if (ndim .eq. 3) then
                pres_test = zr(ipres-1+1)
                do ino = 1, nnop
                    l_pres_var = abs(zr(ipres-1+ino)-pres_test) .ge. r8pre
                    if (l_pres_var) then
                        ASSERT(.false.)
                    end if
                end do
            end if
!
        elseif ((option .eq. 'CALC_K_G_XFEM_F') .or. &
                (option .eq. 'CALC_G_XFEM_F')) then
!
            call utmess('F', 'XFEM_99')
!
        else
            ASSERT(.false.)
        end if
!
!       -----------------------------------
!       6) CALCUL EFFECTIF DE G, K1, K2, K3
!       -----------------------------------
        g = 0.d0
        k1 = 0.d0
        k2 = 0.d0
        k3 = 0.d0
        do j = 1, ndim
            g = g+(forrep(j, ilev)*divt+dfor(j))*depla(j)
            if (option(1:8) .eq. 'CALC_K_G') then
                k1 = k1+(forrep(j, ilev)*divt+dfor(j))*u1(j)
                k2 = k2+(forrep(j, ilev)*divt+dfor(j))*u2(j)
                if (ndim .eq. 3) k3 = k3+(forrep(j, ilev)*divt+dfor(j))*u3(j)
            end if
        end do
!
        jm = jac*0.5d0
!
        if (ndim .eq. 3) then
            zr(igthet-1+1) = zr(igthet-1+1)+g*jac
            if (option(1:8) .eq. 'CALC_K_G') then
                zr(igthet-1+2) = zr(igthet-1+2)+k1*jm*sqrt(coeff)
                zr(igthet-1+3) = zr(igthet-1+3)+k2*jm*sqrt(coeff)
                zr(igthet-1+4) = zr(igthet-1+4)+k3*jm*sqrt(coeff3)
                zr(igthet-1+5) = zr(igthet-1+5)+k1*jm*coeff
                zr(igthet-1+6) = zr(igthet-1+6)+k2*jm*coeff
                zr(igthet-1+7) = zr(igthet-1+7)+k3*jm*coeff3
            end if
        else if (ndim .eq. 2) then
!
            zr(igthet-1+1) = zr(igthet-1+1)+g*jac
!
            if (option(1:8) .eq. 'CALC_K_G') then
                zr(igthet-1+2) = zr(igthet-1+2)+k1*jm*sqrt(coeff)
                zr(igthet-1+3) = zr(igthet-1+3)+k2*jm*sqrt(coeff)
                zr(igthet-1+4) = zr(igthet-1+4)+k1*jm*coeff
                zr(igthet-1+5) = zr(igthet-1+5)+k2*jm*coeff
            end if
!
        end if
!
    end do
!     FIN DE BOUCLE SUR LES DEUX LEVRES
end subroutine
