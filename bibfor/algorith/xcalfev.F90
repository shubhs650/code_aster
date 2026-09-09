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
! aslint: disable=W1306
!
subroutine xcalfev(elrefp, ndim, nnop, basloc, stano, &
                   he, geom, kappa, mu, ff, &
                   fk, dfdi, dkdgl, face, nnop_lin, &
                   ff_lin, dfdi_lin)
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8pi.h"
#include "asterfort/coor_cyl.h"
#include "asterfort/elrfvf.h"
#include "asterfort/is_enr_line.h"
#include "asterfort/utmess.h"
#include "asterfort/xdeffk_wrap.h"
#include "asterfort/xderfk_wrap.h"
#include "asterfort/xelrex.h"
#include "MeshTypes_type.h"
!
    character(len=8), intent(in) :: elrefp
    integer(kind=8) :: ndim, nnop, stano(*)
    real(kind=8) :: he, basloc(*), fk(MT_NNOMAX, 3, 3)
    real(kind=8) :: kappa, ff(*), geom(*), mu
    real(kind=8), optional :: dkdgl(MT_NNOMAX, 3, 3, 3)
    real(kind=8), optional :: dfdi(nnop, ndim)
    character(len=4), optional :: face
    integer(kind=8), optional :: nnop_lin
    real(kind=8), optional :: ff_lin(:)
    real(kind=8), optional :: dfdi_lin(:, :)
!
!
!
!     BUT:  CALCUL DES FONCTIONS D'ENRICHISSEMENT <VECTORIEL> EN UN POINT DE GAUSS
!            DANS LA BASE <GLOBALE>
!
! IN  HE      : VALEUR DE LA FONCTION HEAVYSIDE CSTE LE SS-ELT
! IN  BASLOC  : BASE LOCALE AU FOND DE FISSURE
! IN  KA, MU  : PARAMETRES MATERIAU / LES FONCTIONS ASYMPTIQUE POUR UN MATERIAU ELASTIQUE
! IN  FF      : FONCTIONS DE FORMES DE L ELEMENT PARENT
! IN  DFDI    : DERIVEES DES FONCTIONS DE FORMES DE L ELEMENT PARENT
!
! OUT FK      : VALEURS DES FONCTIONS D'ENRICHISSEMENT <VECTORIEL> DANS LA BASE <GLOBALE>
! OUT DKDGL   : DERIVEES DES FONCTIONS D'ENRICHISSEMENT <VECTORIEL> DANS LA BASE <GLOBALE>
!
!----------------------------------------------------------------
!
    integer(kind=8) :: i, j, k, ino, l, alp, nnops
    integer(kind=8) :: ndime, nno
    real(kind=8) :: p(MT_NNOMAX, 3, 3), invp(MT_NNOMAX, 3, 3), p_g(3, 3), invp_g(3, 3)
    real(kind=8) :: dkdpo(ndim, ndim, 2), dkdlo(3, 3, 2), fkpo(ndim, ndim), fk_gl(ndim, ndim)
    real(kind=8) :: rr, th, r_n(MT_NNOMAX), t_n(MT_NNOMAX), fkpo_n(MT_NNOMAX, 3, 3)
    real(kind=8) :: fkpo_g(3, 3), dkdgl_g(3, 3, 3), signe
    real(kind=8) :: ff1(MT_NNOMAX), dfdi1(MT_NNOMAX, 3)
    real(kind=8) :: xref(3*MT_NNOMAX3D), ff_n(MT_NNOMAX)
    aster_logical :: lderiv, l_not_zero, lshift, lctlin, lbid
    aster_logical :: lcourb
    character(len=8) :: method
!----------------------------------------------------------------
!
    lctlin = is_enr_line()
    lshift = .not. (count(stano(1:nnop) .eq. -2) .eq. nnop)
    method = 'DEFAULT'
    lcourb = .false.
    lbid = ASTER_FALSE
!
    xref(:) = 0.d0
    ff_n(:) = 0.d0
    p(:, :, :) = 0.d0
    invp(:, :, :) = 0.d0
    r_n(:) = 0.d0
    t_n(:) = 0.d0
    p_g(:, :) = 0.d0
    invp_g(:, :) = 0.d0
    ff1(:) = 0.d0
    dfdi1(:, :) = 0.d0
    th = 0.d0
    fk(:, :, :) = 0.d0
    if (.not. present(dkdgl)) then
        lderiv = .false.
    else
        lderiv = .true.
        if (.not. present(dfdi)) then
            call utmess('F', 'ELEMENTS6_6', sk='dfdi')
        end if
        dkdgl(:, :, :, :) = 0.d0
    end if
!
    if (present(nnop_lin) .and. lctlin) then
        if (.not. present(ff_lin)) then
            call utmess('F', 'ELEMENTS6_6', sk='ff_lin')
        end if
        nnops = nnop_lin
        ff1(1:nnops) = ff_lin(1:nnops)
        if (lderiv) then
            if (.not. present(dfdi_lin)) then
                call utmess('F', 'ELEMENTS6_6', sk='dfdi_lin')
            end if
            dfdi1(1:nnops, 1:ndim) = dfdi_lin(1:nnops, 1:ndim)
        end if
    else
        nnops = nnop
        ff1(1:nnops) = ff(1:nnops)
        if (lderiv) dfdi1(1:nnops, 1:ndim) = dfdi(1:nnops, 1:ndim)
    end if
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  PREPARATION DES COORDONNEES CYLINDRIQUES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    lcourb = lcourb .and. lderiv
    if (lderiv) then
        call coor_cyl(ndim, nnop, basloc, geom, ff, &
                      p_g, invp_g, rr, th, l_not_zero, lcourb)
    else
        call coor_cyl(ndim, nnop, basloc, geom, ff, &
                      p_g, invp_g, rr, th, l_not_zero)
    end if
!
    if (.not. l_not_zero) goto 999
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  CORRECTION POUR LE CONTACT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (present(face)) then
        if (.not. (face .eq. 'MAIT' .or. face .eq. 'ESCL' .or. face .eq. ' ')) then
            call utmess('F', 'ELEMENTS6_6', sk='face')
        end if
        if (face .eq. 'MAIT') th = +1.d0*r8pi()
        if (face .eq. 'ESCL') th = -1.d0*r8pi()
        if (face .eq. ' ') th = he*abs(th)
    else
        th = he*abs(th)
    end if
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  CALCUL DU SIGNE POUR L INTERPOLATION FANTOME
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    signe = sign(1.d0, th)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  PREPARATION DU SHIFT: INTERPOLATION AUX NOEUDS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    fkpo_g(:, :) = 0.d0
    fkpo_n(:, :, :) = 0.d0
    if (lshift) then
        call xelrex(elrefp, nno, xref, ndime=ndime)
        do ino = 1, nnop
            call elrfvf(elrefp, xref((ndime*(ino-1)+1):(ndime*(ino-1)+ndime)), ff_n)
            call coor_cyl(ndim, nnop, basloc, geom, ff_n, &
                          p(ino, 1:3, 1:3), invp(ino, 1:3, 1:3), r_n(ino), t_n(ino), lbid)
        end do
!
        do ino = 1, nnop
!     INTERPOLATION FANTOME DANS LES MAILLES XHT
            if (stano(ino) .eq. 1 .or. stano(ino) .eq. 3) then
                call xdeffk_wrap(kappa, mu, r_n(ino), signe*abs(t_n(ino)), ndim, &
                                 fkpo_n(ino, 1:ndim, 1:ndim), method, stano(ino))
!     INTERPOLATION STANDARD DANS LES MAILLES XT
            else if (stano(ino) .eq. 0 .or. stano(ino) .eq. 2) then
                call xdeffk_wrap(kappa, mu, r_n(ino), t_n(ino), ndim, &
                                 fkpo_n(ino, 1:ndim, 1:ndim), method, stano(ino))
            else if (stano(ino) .eq. -2) then
                goto 5
            else
                call utmess('F', 'ELEMENTS6_6', sk='stano')
            end if
            do alp = 1, ndim
                do i = 1, ndim
                    fkpo_g(alp, i) = fkpo_g(alp, i)+fkpo_n(ino, alp, i)*ff(ino)
                end do
            end do
5           continue
        end do
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     FONCTIONS D'ENRICHISSEMENT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  * AU POINT DE GAUSS
    call xdeffk_wrap(kappa, mu, rr, th, ndim, &
                     fkpo, method, 0)
!  * CONVERSION DANS LA BASE GLOBALE
    fk_gl(:, :) = 0.d0
    do alp = 1, ndim
        do i = 1, ndim
            do j = 1, ndim
                fk_gl(alp, i) = fk_gl(alp, i)+p_g(i, j)*(fkpo(alp, j)-fkpo_g(alp, j))
            end do
        end do
    end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  PREPARATION DE LA DERIVATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    lderiv = lderiv .and. l_not_zero
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  PREPARATION DE LA DERIVATION DE LA FONCTION SHIFT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    dkdgl_g(:, :, :) = 0.d0
    if (lderiv) then
        if (lshift) then
            do ino = 1, nnop
                if (stano(ino) .eq. -2) goto 7
                do alp = 1, ndim
                    do i = 1, ndim
                        do j = 1, ndim
                            do k = 1, ndim
                                dkdgl_g(alp, i, j) = dkdgl_g(alp, i, j)+ &
                                                     p_g(i, k)*fkpo_n(ino, alp, k)* &
                                                     dfdi(ino, j)
                            end do
                        end do
                    end do
                end do
7               continue
            end do
        end if
!  ON ROGNE SUR TOUT / CALCUL DE LA DERIVEE EN AMONT
!  *  DERIVEES DES FONCTIONS D'ENRICHISSEMENT DANS LA BASE POLAIRE
        call xderfk_wrap(kappa, mu, rr, th, ndim, &
                         dkdpo, method, 0)
!  *  DERIVEES DES FONCTIONS D'ENRICHISSEMENT DANS LA BASE LOCALE
        do alp = 1, ndim
            do i = 1, ndim
                dkdlo(alp, i, 1) = dkdpo(alp, i, 1)*cos(th)-dkdpo(alp, i, 2)*sin(th)/rr
                dkdlo(alp, i, 2) = dkdpo(alp, i, 1)*sin(th)+dkdpo(alp, i, 2)*cos(th)/rr
            end do
        end do
    end if
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  CALCUL DES FONCTIONS VECTORIELLES AU POINT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do ino = 1, nnop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  POUR LE QUADRATIQUE => BASCULEMENT EN LINEAIRE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (ino .gt. nnops) goto 10
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     DERIVEES DES FONCTIONS D'ENRICHISSEMENT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (.not. lderiv) goto 11
!  *  CONVERSION DANS LA BASE GLOBALE
        do alp = 1, ndim
            do i = 1, ndim
                do j = 1, ndim
                    do k = 1, ndim
                        do l = 1, 2
                            dkdgl(ino, alp, i, j) = dkdgl(ino, alp, i, j)+ &
                                                    p_g(i, k)*dkdlo(alp, k, l)* &
                                                    invp_g(l, j)
                        end do
                    end do
                    dkdgl(ino, alp, i, j) = (dkdgl(ino, alp, i, j)-dkdgl_g(alp, i, j))*ff1(ino)+ &
                                            fk_gl(alp, i)*dfdi1(ino, j)
                end do
            end do
        end do
!
11      continue
!  *  MULTIPLICATION DE FK PAR FF
        do alp = 1, ndim
            do i = 1, ndim
                fk(ino, alp, i) = fk_gl(alp, i)*ff1(ino)
            end do
        end do
!
10      continue
    end do
!
999 continue
!
end subroutine
