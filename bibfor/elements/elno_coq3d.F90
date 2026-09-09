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
subroutine elno_coq3d(plateCara, plateOrie, &
                      lgreen, option, nomte, &
                      nb2, &
                      npgsn, nso, nbLayer, &
                      matrGano, fieldElgaVale, fieldElnoVale)
!
    use plate_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/caurtg.h"
#include "asterfort/pk2cau.h"
#include "asterfort/vdsiro.h"
!
    type(plateCara_Para), intent(in) :: plateCara
    type(plateOrie_Para), intent(in) :: plateOrie
    aster_logical, intent(in):: lgreen
    character(len=16), intent(in) :: option, nomte
    integer(kind=8), intent(in) :: nb2, npgsn, nso, nbLayer
    real(kind=8), intent(in) :: matrGano(*), fieldElgaVale(*)
    real(kind=8), intent(out) :: fieldElnoVale(*)
!
! --------------------------------------------------------------------------------------------------
!
!     CALCUL DES OPTIONS DES ELEMENTS DE COQUE 3D
!     OPTIONS : EPSI_ELNO
!               SIEF_ELNO
!               SIGM_ELNO
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i, ic, icmp, ii, ino
    integer(kind=8) :: inte, intsn, isp, j
    integer(kind=8) :: jj, k1, kpgs, l
    real(kind=8) :: s, thick
    integer(kind=8), parameter :: ncmp = 6, npge = 3
    integer(kind=8) :: iLayer, nordo
    real(kind=8) :: fieldElga(6, 27*nbLayer), fieldElno(6, 12*nbLayer)
    real(kind=8) :: fieldElnoLoca(6, 12*nbLayer)
    real(kind=8) :: pk2(6, 27*nbLayer), matgnu(6, 12*nbLayer), fieldElnoGlob(6, 12*nbLayer)

!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nbLayer .ge. 1)
    thick = plateCara%thick

! - Get values at Gauss points of input field
    kpgs = 0
    do iLayer = 1, nbLayer
        do inte = 1, npge
            do intsn = 1, npgsn
                kpgs = kpgs+1
                k1 = 6*((intsn-1)*npge*nbLayer+ &
                        (iLayer-1)*npge+inte-1)
                do i = 1, 6
                    fieldElga(i, kpgs) = fieldElgaVale(k1+i)
                end do
            end do
        end do
    end do

! - If lGreen, input field is stress !
    if (lgreen) then
! ---- AFFECTATION DES CONTRAINTES DE PIOLA-KIRCHHOFF DE SECONDE ESPECE
        do i = 1, 6
            do j = 1, kpgs
                pk2(i, j) = fieldElga(i, j)
            end do
        end do

! ----- TRANSFORMATION DES CONTRAINTES DE PIOLA-KIRCHHOFF DE
! ----- SECONDE ESPECE PK2 EN CONTRAINTES DE CAUCHY
        call pk2cau(plateOrie, &
                    nbLayer, thick, &
                    nomte, ncmp, &
                    pk2, fieldElga)
    end if

! - Compute ELNO field in local frame
    do iLayer = 1, nbLayer
        do ic = 1, ncmp
            do i = 1, npge*nso
                l = npge*npgsn*(i-1)
                s = 0.d0
                do j = 1, npge*npgsn
                    jj = (iLayer-1)*npge*npgsn+j
                    s = s+matrGano(l+j)*fieldElga(ic, jj)
                end do
                ii = (iLayer-1)*npge*nso+i
                fieldElno(ic, ii) = s
            end do
        end do
    end do

! - PASSAGE DU VECTEUR DES CONTRAINTES DEFINI AUX NOEUDS
! - DE L'ELEMENT DU REPERE INTRINSEQUE AU REPERE UTILISATEUR
    do iLayer = 1, nbLayer
        do nordo = -1, 1
            isp = npge*(iLayer-1)+nordo+2
            do i = 1, ncmp
                do j = 1, nso
                    jj = nso*(nordo+1)+nso*npge*(iLayer-1)+j
                    fieldElnoLoca(i, j) = fieldElno(i, jj)
                end do
                if (nomte .eq. 'MEC3QU9H') then
                    fieldElnoLoca(i, 5) = (fieldElnoLoca(i, 1)+fieldElnoLoca(i, 2))/2.d0
                    fieldElnoLoca(i, 6) = (fieldElnoLoca(i, 2)+fieldElnoLoca(i, 3))/2.d0
                    fieldElnoLoca(i, 7) = (fieldElnoLoca(i, 3)+fieldElnoLoca(i, 4))/2.d0
                    fieldElnoLoca(i, 8) = (fieldElnoLoca(i, 4)+fieldElnoLoca(i, 1))/2.d0
                    fieldElnoLoca(i, 9) = (fieldElnoLoca(i, 1)+fieldElnoLoca(i, 2)+ &
                                           fieldElnoLoca(i, 3)+fieldElnoLoca(i, 4))/4.d0
                else if (nomte .eq. 'MEC3TR7H') then
                    fieldElnoLoca(i, 4) = (fieldElnoLoca(i, 1)+fieldElnoLoca(i, 2))/2.d0
                    fieldElnoLoca(i, 5) = (fieldElnoLoca(i, 2)+fieldElnoLoca(i, 3))/2.d0
                    fieldElnoLoca(i, 6) = (fieldElnoLoca(i, 3)+fieldElnoLoca(i, 1))/2.d0
                    fieldElnoLoca(i, 7) = (fieldElnoLoca(i, 1)+fieldElnoLoca(i, 2)+ &
                                           fieldElnoLoca(i, 3))/3.d0
                end if
            end do

            if (lgreen) then
                call vdsiro(nb2, 1, plateOrie%matevn, 'IU', 'N', &
                            fieldElnoLoca, matgnu)
                call caurtg(nomte, ncmp, matgnu, fieldElnoGlob)
            else
                call vdsiro(nb2, 1, plateOrie%matevn, 'IU', 'N', &
                            fieldElnoLoca, fieldElnoGlob)
            end if
!
            if (option .eq. 'EPSI_ELNO') then
                do icmp = 1, ncmp
                    do ino = 1, nb2
                        fieldElnoVale((ino-1)*ncmp*nbLayer*npge+(isp-1)*ncmp+icmp) = &
                            fieldElnoLoca(icmp, ino)
                    end do
                end do
            else if ((option .eq. 'SIEF_ELNO') .or. ( &
                     option .eq. 'SIGM_ELNO')) then
                do icmp = 1, ncmp
                    do ino = 1, nb2
                        fieldElnoVale((ino-1)*ncmp*nbLayer*npge+(isp-1)*ncmp+icmp) = &
                            fieldElnoGlob(icmp, ino)
                    end do
                end do
            else
                ASSERT(.false.)
            end if
        end do
    end do
end subroutine
