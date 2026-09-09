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
subroutine fpesa(plateCara, &
                 nomte, xi, nb1, &
                 vecl)
!
    use plate_type
    implicit none
!
#include "jeveux.h"
#include "asterfort/dxroep.h"
#include "asterfort/forpes.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/r8inir.h"
#include "asterfort/vectci.h"
#include "asterfort/vexpan.h"
!
    type(plateCara_Para), intent(in) :: plateCara
    character(len=16), intent(in) :: nomte
    real(kind=8), intent(in) :: xi(3, *)
    integer(kind=8), intent(in) :: nb1
    real(kind=8), intent(out) :: vecl(51)
!
    integer(kind=8) :: npgsn
    real(kind=8) :: rho, epais, pesan, rnormc
    real(kind=8) :: vpesan(3), vecl1(42)

!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, intsn, jpesa, lzi, lzr
!-----------------------------------------------------------------------
    call jevech('PPESANR', 'L', jpesa)
    pesan = zr(jpesa)
    do i = 1, 3
        vpesan(i) = pesan*zr(jpesa+i)
    end do
!
    call jevete('&INEL.'//nomte(1:8)//'.DESI', ' ', lzi)
    npgsn = zi(lzi-1+4)
!
    call jevete('&INEL.'//nomte(1:8)//'.DESR', ' ', lzr)

! - Get plate parameters
    call dxroep(plateCara, rho, epais)
!
    call r8inir(42, 0.d0, vecl1, 1)
!
    do intsn = 1, npgsn
        call vectci(intsn, nb1, xi, zr(lzr), rnormc)
        call forpes(intsn, nb1, zr(lzr), rho, epais, &
                    vpesan, rnormc, vecl1)
    end do
!
    call vexpan(nb1, vecl1, vecl)
    do i = 1, 3
        vecl(6*nb1+i) = 0.d0
    end do
!
end subroutine
