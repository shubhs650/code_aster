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
subroutine vdxedg(plateCara, plateOrie, &
                  nomte, option, nodeCoor, &
                  degeElga, degeElno)
!
    use plate_type
    implicit none
!
#include "asterfort/btdfn.h"
#include "asterfort/btdmsn.h"
#include "asterfort/btdmsr.h"
#include "asterfort/hsj1f.h"
#include "asterfort/hsj1ms.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/mahsf.h"
#include "asterfort/mahsms.h"
#include "asterfort/trndgl.h"
#include "asterfort/vddege.h"
#include "asterfort/vectan.h"
#include "asterfort/vectgt.h"
#include "jeveux.h"
!
    type(plateCara_Para), intent(in) :: plateCara
    type(plateOrie_Para), intent(in) :: plateOrie
    character(len=16), intent(in) :: nomte, option
    real(kind=8), intent(in) :: nodeCoor(3, 9)
    real(kind=8), intent(out) :: degeElga(72), degeElno(8, 9)
!
! --------------------------------------------------------------------------------------------------
!
! COQUE_3D
!
! Compute DEGE_ELGA and DEGE_ELNO
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), parameter :: zero = 0.d0, un = 1.d0
    real(kind=8), parameter :: btdfMembZero(4, 3, 42) = 0.d0
    real(kind=8), parameter :: btdfShearZero(4, 2, 42) = 0.d0
    real(kind=8), parameter :: btdfBendZero(3, 42) = 0.d0
    integer(kind=8) :: nb1, nb2, npgsn, npgsr
    integer(kind=8) :: kpgsn, kpgsr, kwgt
    integer(kind=8) :: jvDisp
    integer(kind=8) :: i, k
    integer(kind=8) :: lzi, lzr
    real(kind=8) :: vectTangKpg(2, 3), vectBaseKpg(3, 3)
    real(kind=8) :: hsfm(3, 9), hss(2, 9), hsj1m(3, 9), hsj1s(2, 9)
    real(kind=8) :: hsf(3, 9), hsj1fx(3, 9), wgt
    real(kind=8) :: btdm(4, 3, 42), btdf(3, 42), btds(4, 2, 42)
    real(kind=8) :: btild(5, 42), btild1(5, 42)
    real(kind=8) :: disp(42), rotf(9)
    real(kind=8) :: epsiBend(5), epsiMemb(5)
    real(kind=8) :: epais
!
! --------------------------------------------------------------------------------------------------
!
    degeElga = 0.d0
    degeElno = 0.d0

! - Access to static objects of COQUE_3D
    call jevete('&INEL.'//nomte(1:8)//'.DESI', ' ', lzi)
    nb1 = zi(lzi-1+1)
    nb2 = zi(lzi-1+2)
    npgsr = zi(lzi-1+3)
    npgsn = zi(lzi-1+4)
    call jevete('&INEL.'//nomte(1:8)//'.DESR', ' ', lzr)

! - Get plate parameters
    epais = plateCara%thick

! - Get displacements
    call jevech('PDEPLAR', 'L', jvDisp)

! - Change coordinates of displacements/rotations
    call trndgl(nb2, plateOrie%vectNorm, plateOrie%vectTang, zr(jvDisp), &
                disp, rotf)
!
    kwgt = 0
! - MEMBRANE ET CISAILLEMENT
    do kpgsr = 1, npgsr
        call mahsms(plateOrie, &
                    0, nb1, &
                    nodeCoor, un, kpgsr, &
                    zr(lzr), epais, &
                    vectBaseKpg, vectTangKpg, &
                    hsfm, hss)
        call hsj1ms(epais, vectTangKpg, vectBaseKpg, hsfm, hss, &
                    hsj1m, hsj1s)
        call btdmsr(nb1, nb2, un, kpgsr, zr(lzr), &
                    epais, plateOrie%vectTang, hsj1m, hsj1s, btdm, &
                    btds)
    end do

! - FLEXION
    do kpgsn = 1, npgsn
        call mahsf(plateOrie, &
                   1, nb1, &
                   nodeCoor, un, kpgsn, &
                   zr(lzr), epais, &
                   vectBaseKpg, vectTangKpg, &
                   hsf)
        call hsj1f(kpgsn, zr(lzr), epais, vectTangKpg, vectBaseKpg, &
                   hsf, kwgt, hsj1fx, wgt)
        call btdfn(1, nb1, nb2, un, kpgsn, &
                   zr(lzr), epais, plateOrie%vectTang, hsj1fx, btdf)

! ----- Final btild/btild1
        call btdmsn(1, nb1, kpgsn, npgsr, zr(lzr), &
                    btdm, btdfBendZero, btdfShearZero, btild1)
        call btdmsn(1, nb1, kpgsn, npgsr, zr(lzr), &
                    btdfMembZero, btdf, btds, btild)

! ----- Strain (membrane)
        epsiMemb = zero
        do i = 1, 5
            do k = 1, 5*nb1+2
                epsiMemb(i) = epsiMemb(i)+btild1(i, k)*disp(k)
            end do
        end do

! ----- Strain (bending)
        epsiBend = zero
        do i = 1, 5
            do k = 1, 5*nb1+2
                epsiBend(i) = epsiBend(i)+btild(i, k)*disp(k)
            end do
        end do
        epsiBend(1) = epsiBend(1)/epais
        epsiBend(2) = epsiBend(2)/epais
        epsiBend(3) = epsiBend(3)/epais
        epsiBend(4) = epsiBend(4)
        epsiBend(5) = epsiBend(5)

! ----- Save DEGE_ELGA
        degeElga((kpgsn-1)*8+1) = epsiMemb(1)
        degeElga((kpgsn-1)*8+2) = epsiMemb(2)
        degeElga((kpgsn-1)*8+3) = epsiMemb(3)
        degeElga((kpgsn-1)*8+4) = epsiBend(1)
        degeElga((kpgsn-1)*8+5) = epsiBend(2)
        degeElga((kpgsn-1)*8+6) = epsiBend(3)
        degeElga((kpgsn-1)*8+7) = epsiBend(4)
        degeElga((kpgsn-1)*8+8) = epsiBend(5)

    end do

! - Compute DEGE_ELNO
    if (option .eq. 'DEGE_ELNO') then
        call vddege(nomte, nb2, npgsn, zr(lzr), degeElga, &
                    degeElno)
    end if
!
end subroutine
