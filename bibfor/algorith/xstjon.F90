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
subroutine xstjon(elrefp, ndim, joncno, jlsn, igeom, &
                  nfiss, nfisc, fisco, nnops, txlsn, &
                  n, c)
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "blas/ddot.h"
#include "MeshTypes_type.h"
#include "asterfort/conare.h"
#include "asterfort/confac.h"
#include "asterfort/provec.h"
#include "asterfort/reeref.h"
#include "asterfort/tecael.h"
#include "asterfort/xelrex.h"
#include "asterfort/xnormv.h"
#include "MeshTypes_type.h"
    character(len=8) :: elrefp
    integer(kind=8) :: ndim, joncno, jlsn, nfiss, nfisc, nnops, fisco(*), igeom
    real(kind=8) :: txlsn(28)
    integer(kind=8), intent(in), optional :: n
    real(kind=8), intent(in), optional :: c(ndim)
!
! ======================================================================
!             DETERMINER LES NOEUDS SOMMETS DE L'ELEMENT PARENT
!             QUI VOIENT UNE JONCTION (POUR LE CONTACT MULTI)
!
!     ENTREE  N     : NUMERO DU NOEUD DANS L'ELEMENT PARENT SI LE POINT
!                     D'INTERSECTION EST CONFONDU AVEC UN NOEUD
!             C     : COORDONNES DE REFERENCE DU POINT D'INTERSECTION
!
!     SORTIE  PTJONC: LE POINT CONSIDERE EST-IL UN POINT DE JONCTION?
!
!---------------------------------------------------------------------
!
    integer(kind=8) :: i, iadzi, iazk24, ft(12, 3), nbft, f(6, 8), nbf
    integer(kind=8) :: ar(12, 3), nbar, iar, nno, ino
    real(kind=8) :: xref(3*MT_NNOMAX3D), u(3), v(3), w(3), norme, normal(3), cridist
    real(kind=8) :: cref(ndim), ff(MT_NNOMAX), val
    character(len=8) :: typma
    aster_logical :: jonc, arete, face
    blas_int :: b_incx, b_incy, b_n
    parameter(cridist=1.d-8)
!
!---------------------------------------------------------------------
!     DEBUT
!---------------------------------------------------------------------
!
    call tecael(iadzi, iazk24, noms=0)
    typma = zk24(iazk24-1+3+zi(iadzi-1+2)+3) (1:8)
    call xelrex(elrefp, nno, xref)
!
    jonc = .false.
    arete = .false.
    face = .false.
!
    if (present(n)) then
        if (n .le. nnops) then
            do i = 1, nfisc
                if (zr(jlsn-1+(n-1)*nfiss+fisco(2*i-1)) .eq. 0.d0) then
                    zi(joncno-1+n) = 1
                end if
            end do
        else if (n .le. nno) then
            call conare(typma, ar, nbar)
            do i = 1, nfisc
                if (zr(jlsn-1+(n-1)*nfiss+fisco(2*i-1)) .eq. 0.d0) then
                    do iar = 1, nbar
                        if (n .eq. ar(iar, 3)) then
                            zi(joncno-1+ar(iar, 1)) = 1
                            zi(joncno-1+ar(iar, 2)) = 1
                        end if
                    end do
                end if
            end do
        else if (n .gt. nno .and. n .lt. 1000) then
            do i = 1, nfisc
                if (txlsn((n-nno-1)*nfiss+fisco(2*i-1)) .eq. 0.d0) then
                    zi(joncno-1+n) = 1
                end if
            end do
        else if (n .lt. 2000) then
            call reeref(elrefp, nno, zr(igeom), c, ndim, &
                        cref, ff)
            do i = 1, nfisc
                val = 0.d0
                do ino = 1, nno
                    val = val+zr(jlsn-1+(ino-1)*nfiss+fisco(2*i-1))*ff(ino)
                end do
                if (abs(val) .lt. cridist) jonc = .true.
            end do
        else if (n .lt. 3000) then
            call reeref(elrefp, nno, zr(igeom), c, ndim, &
                        cref, ff)
            do i = 1, nfisc
                val = 0.d0
                do ino = 1, nno
                    val = val+zr(jlsn-1+(ino-1)*nfiss+fisco(2*i-1))*ff(ino)
                end do
                if (abs(val) .lt. cridist) jonc = .true.
            end do
        end if
    else
        jonc = .true.
    end if
!
    if (jonc) then
        call conare(typma, ar, nbar)
        do iar = 1, nbar
            u(:) = 0.d0
            v(:) = 0.d0
            do i = 1, ndim
                u(i) = xref((ar(iar, 2)-1)*ndim+i)-xref((ar(iar, 1)-1)*ndim+i)
                v(i) = c(i)-xref((ar(iar, 1)-1)*ndim+i)
            end do
            call xnormv(ndim, u, norme)
            call xnormv(ndim, v, norme)
            call provec(u, v, normal)
            call xnormv(ndim, normal, norme)
            if (norme .lt. cridist) then
                zi(joncno-1+ar(iar, 1)) = 1
                zi(joncno-1+ar(iar, 2)) = 1
                arete = .true.
            end if
        end do
        if (ndim .eq. 3) then
            call confac(typma, ft, nbft, f, nbf)
            do iar = 1, nbf
                u(:) = 0.d0
                v(:) = 0.d0
                do i = 1, ndim
                    u(i) = xref((f(iar, 2)-1)*ndim+i)-xref((f(iar, 1)-1)*ndim+i)
                    v(i) = xref((f(iar, 3)-1)*ndim+i)-xref((f(iar, 1)-1)*ndim+i)
                    w(i) = c(i)-xref((f(iar, 1)-1)*ndim+i)
                end do
                call xnormv(ndim, u, norme)
                call xnormv(ndim, v, norme)
                call xnormv(ndim, w, norme)
                call provec(u, v, normal)
                call xnormv(ndim, normal, norme)
                b_n = to_blas_int(ndim)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                if (abs(ddot(b_n, w, b_incx, normal, b_incy)) .lt. cridist) then
                    zi(joncno-1+f(iar, 1)) = 1
                    zi(joncno-1+f(iar, 2)) = 1
                    zi(joncno-1+f(iar, 3)) = 1
                    if (f(iar, 4) .ge. 1) then
                        zi(joncno-1+f(iar, 4)) = 1
                    end if
                    face = .true.
                end if
            end do
        end if
    end if
!
    if (jonc .and. .not. arete .and. .not. face) then
        do i = 1, nnops
            zi(joncno-1+i) = 1
        end do
    end if
!---------------------------------------------------------------------
!     FIN
!---------------------------------------------------------------------
end subroutine
