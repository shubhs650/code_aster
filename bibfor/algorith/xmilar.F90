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
subroutine xmilar(ndim, ndime, elrefp, geom, pinref, &
                  ia, ib, im, ip, ksia, &
                  ksib, milara, milarb, pintt, pmitt)
!
    implicit none
!
#include "asterfort/xelrex.h"
#include "asterfort/reeref.h"
#include "asterfort/reerel.h"
#include "asterfort/xnormv.h"
#include "blas/ddot.h"
#include "MeshTypes_type.h"
    integer(kind=8) :: ndim, ndime, ia, ib, im, ip
    real(kind=8) :: milara(3), milarb(3), pinref(*), geom(*)
    real(kind=8) :: ksia(ndime), ksib(ndime), pintt(*), pmitt(*)
    character(len=8) :: elrefp
!
!                      TROUVER LES PTS MILIEUX ENTRE LES EXTREMITES DE
!                      L'ARETE ET LE POINT D'INTERSECTION
!
!     ENTREE
!       NDIM    : DIMENSION TOPOLOGIQUE DU MAILLAGE
!       PINTER  : COORDONNÉES DES POINTS D'INTERSECTION
!       TABAR   : COORDONNEES DES 3 NOEUDS DE L'ARETE
!       AREINT  : POSITION DU PT INTER DE L'ARETE DANS LA LISTE
!       PINTT   : COORDONNEES REELES DES POINTS D'INTERSECTION
!       PINTT   : COORDONNEES REELES DES POINTS MILIEUX
!
!     SORTIE
!       MILARA  : COOR DU PT MILIEU ENTRE 1ER PT DE COORSG ET PT INTER
!       MILARB  : COOR DU PT MILIEU ENTRE 2EM PT DE COORSG ET PT INTER
!     ----------------------------------------------------------------
!
    integer(kind=8) :: nno, j
    real(kind=8) :: x(3*MT_NNOMAX3D), newpt(ndim), pta(ndim), ptb(ndim), ptm(ndim), ff(MT_NNOMAX)
    real(kind=8) :: ab(ndime), aip(ndime), normab, normaip, s
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------
!
    call xelrex(elrefp, nno, x)
    if (ia .lt. 1000) then
        do j = 1, ndime
            pta(j) = x(ndime*(ia-1)+j)
        end do
    else
        do j = 1, ndim
            newpt(j) = pintt(ndim*(ia-1001)+j)
        end do
        call reeref(elrefp, nno, geom, newpt, ndim, &
                    pta, ff)
    end if
    if (ib .lt. 1000) then
        do j = 1, ndime
            ptb(j) = x(ndime*(ib-1)+j)
        end do
    else
        do j = 1, ndim
            newpt(j) = pintt(ndim*(ib-1001)+j)
        end do
        call reeref(elrefp, nno, geom, newpt, ndim, &
                    ptb, ff)
    end if
    if (im .lt. 2000) then
        do j = 1, ndime
            ptm(j) = x(ndime*(im-1)+j)
        end do
    else
        do j = 1, ndim
            newpt(j) = pmitt(ndim*(im-2001)+j)
        end do
        call reeref(elrefp, nno, geom, newpt, ndim, &
                    ptm, ff)
    end if
!
    if (im .lt. 2000) then
        do j = 1, ndime
            ksia(j) = (pinref(ndime*(ip-1)+j)+pta(j))/2.d0
            ksib(j) = (pinref(ndime*(ip-1)+j)+ptb(j))/2.d0
        end do
    else
!   RECHERCHE DE L'ABSCISSE CURVILIGNE DE IP SUR LE SEGMENT AB
        do j = 1, ndime
            ab(j) = ptb(j)-pta(j)
            aip(j) = pinref(ndime*(ip-1)+j)-pta(j)
        end do
        call xnormv(ndime, ab, normab)
        call xnormv(ndime, aip, normaip)
        b_n = to_blas_int(ndime)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        s = normaip/normab*ddot(b_n, ab, b_incx, aip, b_incy)
        do j = 1, ndime
            ksia(j) = (1.d0-s)*(1.d0-s/2.d0)*pta(j)+s/2.d0*(s-1.d0)*ptb(j)+s*(2.d0-s)*ptm(j)
            ksib(j) = s/2.d0*(s-1.d0)*pta(j)+s/2.d0*(s+1.d0)*ptb(j)+(s+1.d0)*(1.d0-s)*ptm(j)
        end do
    end if
!
    call reerel(elrefp, nno, ndim, geom, ksia, &
                milara)
    call reerel(elrefp, nno, ndim, geom, ksib, &
                milarb)
!
end subroutine
