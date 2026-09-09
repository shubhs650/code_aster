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
subroutine efcoq3d(plateCara, plateOrie, &
                   nomte, nb1, nb2, &
                   npgsn, npgsr, npge, nso, &
                   nodeCoor, &
                   desr, siefElga, matrGano, &
                   efgeElno)
!
    use plate_type
    implicit none
!
#include "asterfort/utmess.h"
#include "asterfort/vdefgn.h"
#include "asterfort/vdefro.h"
#include "asterfort/vectgt.h"
!
    type(plateCara_Para), intent(in) :: plateCara
    type(plateOrie_Para), intent(in) :: plateOrie
    character(len=16), intent(in) :: nomte
    integer(kind=8), intent(in) :: nb1, nb2
    integer(kind=8), intent(in) :: npgsn, npgsr, npge, nso
    real(kind=8), intent(in) :: nodeCoor(*)
    real(kind=8), intent(inout) :: desr(*)
    real(kind=8), intent(in) :: siefElga(*), matrGano(*)
    real(kind=8), intent(out) :: efgeElno(*)
!
! --------------------------------------------------------------------------------------------------
!
!     CALCUL DE EFGE_ELNO
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i, ic, icomp, ii
    integer(kind=8) :: inte, intsn, intsr, isom
    integer(kind=8) :: j, jj
    integer(kind=8) :: k, k1, kpgs, l
    integer(kind=8), parameter :: ncmp = 6
    real(kind=8) :: hLayer, s, zic, zmin, epais
    integer(kind=8) :: iLayer, nblayer
    real(kind=8) :: vectBaseKpg(3, 3)
    real(kind=8) :: sigm(6, 270), sigma(6, 120), effgc(8, 9), efgeElnoLoca(8, 9)
    real(kind=8), parameter :: zero = 0.d0
!
! --------------------------------------------------------------------------------------------------
!
    nbLayer = plateCara%nbLayer
    epais = plateCara%thick
    if (nbLayer .le. 0) then
        call utmess('F', 'PLATE1_10')
    end if
    zmin = -epais/2.d0
    hLayer = epais/nbLayer
!
    kpgs = 0
    do iLayer = 1, nbLayer
        do inte = 1, npge
            if (inte .eq. 1) then
                zic = zmin+(iLayer-1)*hLayer
            else if (inte .eq. 2) then
                zic = zmin+hLayer/2.d0+(iLayer-1)*hLayer
            else
                zic = zmin+hLayer+(iLayer-1)*hLayer
            end if
!
            do intsn = 1, npgsn
                kpgs = kpgs+1
                k1 = 6*((intsn-1)*npge*nbLayer+ &
                        (iLayer-1)*npge+inte-1)
                do i = 1, 6
                    sigm(i, kpgs) = siefElga(k1+i)
                end do
            end do
        end do
    end do

! - Compute local basis at integration points
    k = 0
    do intsr = 1, npgsr
        call vectgt(plateOrie, 0, nb1, &
                    nodeCoor, zero, intsr, &
                    epais, desr, &
                    vectBaseKpg)
        do j = 1, 3
            do i = 1, 3
                k = k+1
                desr(1+2000+k-1) = vectBaseKpg(i, j)
            end do
        end do
    end do

! - EXTRAPOLATION VERS LES NOEUDS SOMMETS
    do iLayer = 1, nbLayer
        do ic = 1, ncmp
            do i = 1, npge*nso
                l = npge*npgsn*(i-1)
                s = 0.d0
                do j = 1, npge*npgsn
                    jj = (iLayer-1)*npge*npgsn+j
                    s = s+matrGano(l+j)*sigm(ic, jj)
                end do
                ii = (iLayer-1)*npge*nso+i
                sigma(ic, ii) = s
            end do
        end do
    end do

    do i = 1, nb2
        do j = 1, 8
            efgeElnoLoca(j, i) = 0.d0
        end do
    end do

    do ic = 1, nbLayer
        j = (ic-1)*npge*nso+1
        zic = zmin+(ic-1)*hLayer
        call vdefgn(nomte, nb2, hLayer, zic, sigma(1, j), &
                    effgc)
        do isom = 1, nb2
            do icomp = 1, 8
                efgeElnoLoca(icomp, isom) = efgeElnoLoca(icomp, isom)+effgc(icomp, isom)
            end do
        end do
    end do

! - PASSAGE DU VECTEUR DES EFFORTS GENERALISES DEFINI AUX NOEUDS
! - DE L'ELEMENT DU REPERE INTRINSEQUE AU REPERE UTILISATEUR :
    call vdefro(nb2, plateOrie%matevn, efgeElnoLoca, efgeElno)
!
end subroutine
