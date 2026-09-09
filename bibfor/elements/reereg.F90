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
subroutine reereg(stop, elrefp, nnop, coor, xg, &
                  ndim, xe, iret, toler, ndim_coor_)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/elrfdf.h"
#include "asterfort/elrfvf.h"
#include "asterfort/invjax.h"
#include "asterfort/utmess.h"
#include "MeshTypes_type.h"
#include "blas/ddot.h"
!
    character(len=1) :: stop
    character(len=8) :: elrefp
    integer(kind=8) :: nnop, ndim
    real(kind=8) :: coor(*)
    real(kind=8) :: xg(ndim)
    real(kind=8) :: xe(ndim)
    real(kind=8), optional, intent(in) :: toler
    integer(kind=8) :: iret
    integer(kind=8), optional, intent(in) :: ndim_coor_
!
! ----------------------------------------------------------------------
!
! TROUVER LES COORDONNEES DANS L'ELEMENT DE REFERENCE D'UN
! POINT DONNE DANS L'ESPACE REEL PAR LA METHODE DE NEWTON
!
! ----------------------------------------------------------------------
!
!
! IN  STOP   : /'S' : ON S'ARRETE EN ERREUR <F> EN CAS D'ECHEC
!              /'C' : ON CONTINUE EN CAS D'ECHEC (IRET=1)
! IN  ELREFP : TYPE DE L'ELEMENT
! IN  NNOP   : NOMBRE DE NOEUDS DE L'ELEMENT
! IN  COOR   : COORDONNEES DS ESPACE REEL DES NOEUDS DE L'ELEMENT
! IN  XG     : COORDONNEES DU POINT DANS L'ESPACE REEL
! IN  NDIM   : DIMENSION DE L'ESPACE
! OUT XE     : COORDONNEES DU POINT DANS L'ESPACE PARA DE L'ELEMENT
! OUT IRET   : 0 : ON A CONVERGE
!            : 1 : ON N'A PAS CONVERGE
!
! ----------------------------------------------------------------------
!
    integer(kind=8) :: nbnomx, itermx, ndim_coor
    parameter(nbnomx=MT_NNOMAX3D, itermx=50)
!
    real(kind=8) :: zero, tolerc
    integer(kind=8) :: iter, i, k, idim, ino, ipb
    integer(kind=8) :: nno, nderiv
    real(kind=8) :: etmp(3), err
    real(kind=8) :: point(ndim), xenew(ndim), invjac(3, 3)
    real(kind=8) :: dff(3, nbnomx)
    real(kind=8) :: ff(nnop)
    blas_int :: b_incx, b_incy, b_n
!
! ----------------------------------------------------------------------
!
! --- si tolerance non precisee, defaut a 10E-8
!
    if (present(toler)) then
        tolerc = toler
    else
        tolerc = 1.d-8
    end if
!
    if (present(ndim_coor_)) then
        ndim_coor = ndim_coor_
    else
        ndim_coor = ndim
    end if
!
! --- INITIALISATIONS
!
    zero = 0.d0
    iter = 0
    iret = 0
    xe(:) = zero
!
100 continue
    iter = iter+1
!
! --- VALEURS DES FONCTIONS DE FORME EN XE: FF
!
    call elrfvf(elrefp, xe, ff, nno)
    ASSERT(nno .eq. nnop)
!
! --- DERIVEES PREMIERES DES FONCTIONS DE FORME EN XE: DFF
!
    call elrfdf(elrefp, xe, dff, nno, nderiv)
!      ASSERT(NDERIV.EQ.NDIM)
!
! --- CALCUL DES COORDONNEES DU POINT: POINT
!
    point(:) = zero
    do idim = 1, ndim
        do ino = 1, nno
            point(idim) = point(idim)+ff(ino)*coor(ndim_coor*(ino-1)+idim)
        end do
    end do
!
! --- CALCUL DE L'INVERSE DE LA JACOBIENNE EN XE: INVJAC
!
    call invjax(stop, nno, ndim, nderiv, dff, &
                coor, invjac, ipb, ndim_coor)
    if (ipb .eq. 1) then
        if (stop .eq. 'S') then
            call utmess('F', 'ALGORITH5_19')
        else if (stop .eq. 'C') then
            iret = 1
            goto 999
        else
            ASSERT(.false.)
        end if
    end if
!
! --- UPDATE XE
!
    do i = 1, ndim
        xenew(i) = xe(i)
        do k = 1, ndim
            xenew(i) = xenew(i)-invjac(i, k)*(point(k)-xg(k))
        end do
    end do
!
! --- CALCUL DE L'ERREUR: ERR
!
    do i = 1, ndim
        etmp(i) = xenew(i)-xe(i)
    end do
    b_n = to_blas_int(nderiv)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    err = ddot(b_n, etmp, b_incx, etmp, b_incy)
!
! --- NOUVELLE VALEUR DE XE
!
    xe(1:ndim) = xenew(1:ndim)
!
! --- TEST DE FIN DE BOUCLE
!
    if (err .le. tolerc) then
        goto 999
    else if (iter .lt. itermx) then
        goto 100
    else
        if (stop .eq. 'S') then
            call utmess('F', 'ELEMENTS2_58')
        else
            ASSERT(stop .eq. 'C')
            iret = 1
        end if
    end if
!
999 continue
!
end subroutine
