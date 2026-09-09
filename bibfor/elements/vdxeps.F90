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
subroutine vdxeps(plateCara, plateOrie, &
                  nomte, nodeCoor, &
                  nbLayer, epsiElga)
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
#include "asterfort/vdesga.h"
#include "asterfort/vectgt.h"
#include "jeveux.h"
!
    type(plateCara_Para), intent(in) :: plateCara
    type(plateOrie_Para), intent(in) :: plateOrie
    character(len=16), intent(in) :: nomte
    real(kind=8), intent(in) :: nodeCoor(3, 9)
    integer(kind=8), intent(in) :: nbLayer
    real(kind=8), intent(out) :: epsiElga(6*27*nbLayer)
!
! --------------------------------------------------------------------------------------------------
!
! COQUE_3D
!
! Compute EPSI_ELGA
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), parameter :: un = 1.d0, deux = 2.d0
    integer(kind=8), parameter :: npgLayer = 3
    integer(kind=8) :: nb1, nb2, npgsn, npgsr
    integer(kind=8) :: i
    integer(kind=8) :: jvDisp
    integer(kind=8) :: lzi, lzr
    integer(kind=8) :: iLayer, kpgLayer, kpgsn, kpgsr, kpgs, kwgt
    real(kind=8) :: epsiKpg(6*27*nbLayer)
    real(kind=8) :: vectTangKpg(2, 3), vectBaseKpg(3, 3)
    real(kind=8) :: hsfm(3, 9), hss(2, 9), hsj1m(3, 9), hsj1s(2, 9)
    real(kind=8) :: btdm(4, 3, 42), btds(4, 2, 42)
    real(kind=8) :: hsf(3, 9), hsj1fx(3, 9), wgt
    real(kind=8) :: btdf(3, 42), btild(5, 42)
    real(kind=8) :: disp(42), rotf(9)
    real(kind=8) :: epais
    real(kind=8) :: ksi3s2, hic, zmin
!
! --------------------------------------------------------------------------------------------------
!
    epsiElga = 0.d0

! - Access to static objects of COQUE_3D
    call jevete('&INEL.'//nomte(1:8)//'.DESI', ' ', lzi)
    nb1 = zi(lzi-1+1)
    nb2 = zi(lzi-1+2)
    npgsr = zi(lzi-1+3)
    npgsn = zi(lzi-1+4)
    call jevete('&INEL.'//nomte(1:8)//'.DESR', ' ', lzr)

! - Get plate parameters
    epais = plateCara%thick

! - Get properties of shell
    hic = un/nbLayer
    zmin = -0.5d0

! - Get displacements
    call jevech('PDEPLAR', 'L', jvDisp)

! - Compute local basis
    call trndgl(nb2, plateOrie%vectNorm, plateOrie%vectTang, zr(jvDisp), disp, &
                rotf)
!
    kwgt = 0
    kpgs = 0
    do iLayer = 1, nbLayer
        do kpgLayer = 1, npgLayer
! --------- COTE DES POINTS D'INTEGRATION
            if (kpgLayer .eq. 1) then
                ksi3s2 = zmin+(iLayer-1)*hic
            else if (kpgLayer .eq. 2) then
                ksi3s2 = zmin+hic/deux+(iLayer-1)*hic
            else
                ksi3s2 = zmin+hic+(iLayer-1)*hic
            end if

! --------- MEMBRANE ET CISAILLEMENT
            do kpgsr = 1, npgsr
                call mahsms(plateOrie, &
                            0, nb1, &
                            nodeCoor, ksi3s2, kpgsr, &
                            zr(lzr), epais, &
                            vectBaseKpg, vectTangKpg, &
                            hsfm, hss)
                call hsj1ms(epais, vectTangKpg, vectBaseKpg, hsfm, hss, &
                            hsj1m, hsj1s)
                call btdmsr(nb1, nb2, ksi3s2, kpgsr, zr(lzr), &
                            epais, plateOrie%vectTang, hsj1m, hsj1s, btdm, &
                            btds)
            end do

            do kpgsn = 1, npgsn
! ------------- FLEXION, CISAILLEMENT, NORMAL
                call mahsf(plateOrie, &
                           1, nb1, &
                           nodeCoor, ksi3s2, kpgsn, &
                           zr(lzr), epais, &
                           vectBaseKpg, vectTangKpg, &
                           hsf)
                call hsj1f(kpgsn, zr(lzr), epais, vectTangKpg, vectBaseKpg, &
                           hsf, kwgt, hsj1fx, wgt)
                call btdfn(1, nb1, nb2, ksi3s2, kpgsn, &
                           zr(lzr), epais, plateOrie%vectTang, hsj1fx, btdf)

! ------------- Final btild
                call btdmsn(1, nb1, kpgsn, npgsr, zr(lzr), &
                            btdm, btdf, btds, btild)

! ------------- Compute strain
                call vdesga(plateCara, plateOrie, &
                            kwgt, nb1, &
                            disp, btild, vectBaseKpg, &
                            epsiKpg_=epsiKpg)
!
                kpgs = kpgs+1
                do i = 1, 6
                    epsiElga(6*((kpgsn-1)*npgLayer*nbLayer+npgLayer*(iLayer-1)+kpgLayer-1)+i) = &
                        epsiKpg(i+6*(kpgs-1))
                end do
            end do
        end do
    end do
!
end subroutine
