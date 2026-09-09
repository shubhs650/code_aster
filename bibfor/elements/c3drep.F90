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
subroutine c3drep(nomte, &
                  epais, alpha, beta, &
                  nodeCoor, nodeNume, &
                  pgl)
!
    use plate_type
    implicit none
!
#include "asterfort/jeveuo.h"
#include "asterfort/vdrep2.h"
#include "asterfort/vectan.h"
#include "asterfort/vectgt.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: nomte
    real(kind=8), intent(in) :: epais, alpha, beta
    real(kind=8), intent(in) :: nodeCoor(3, *)
    integer(kind=8), intent(in) :: nodeNume
    real(kind=8), intent(out) :: pgl(3, 3)
!
! --------------------------------------------------------------------------------------------------
!
!         CETTE ROUTINE REALISE LA MEME TACHE QUE COQREP MAIS POUR LES
!         COQUES 3D
!         CALCUL DE LA MATRICE DE PASSAGE DU REPERE DE L'ELEMENT A
!         LA VARIETE (LE REPERE DE LA VARIETE EST OBTENU PAR LA MATRICE
!         DE PASSAGE GLOBAL -> LOCAL) AINSI QUE SON INVERSE
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), parameter:: zero = 0.d0
    integer(kind=8), parameter :: ptType = 0
    integer(kind=8) :: nb1, nb2, npgsr, i, j, k, intsr
    real(kind=8) :: vectNorm(9, 3), vectBaseKpg(3, 3)
    real(kind=8) :: vectTang(9, 2, 3), vectmp(3, 3), pgltmp(3, 3)
    real(kind=8) :: matevn(2, 2, 10), v
    real(kind=8), pointer :: desr(:) => null()
    integer(kind=8), pointer :: desi(:) => null()
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!

! - Access to static objects of COQUE_3D
    call jeveuo('&INEL.'//nomte(1:8)//'.DESI', 'L', vi=desi)
    call jeveuo('&INEL.'//nomte(1:8)//'.DESR', 'L', vr=desr)
    nb1 = desi(1)
    nb2 = desi(2)
    npgsr = desi(3)

! - Compute local basis
    call vectan(nb1, nb2, &
                nodeCoor, desr, &
                vectNorm, vectTang)
    plateOrie%vectNorm = vectNorm
    plateOrie%vectTang = vectTang

! - Compute local basis at integration points
    k = 0
    do intsr = 1, npgsr
        call vectgt(plateOrie, ptType, nb1, &
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

! - Compute global<=>local matrices (at nodes and integration points)
    call vdrep2(alpha, beta, nb2, npgsr, desr, matevn)
!
    vectmp(1, 1) = matevn(1, 1, nodeNume)
    vectmp(1, 2) = matevn(1, 2, nodeNume)
    vectmp(2, 1) = matevn(2, 1, nodeNume)
    vectmp(2, 2) = matevn(2, 2, nodeNume)
    vectmp(1, 3) = 0.d0
    vectmp(2, 3) = 0.d0
    vectmp(3, 3) = 1.d0
    vectmp(3, 1) = 0.d0
    vectmp(3, 2) = 0.d0

! - Compute PGL
    k = 0
    do j = 1, 3
        do i = 1, 3
            k = k+1
            pgltmp(i, j) = desr(1+1090+(nodeNume-1)*9+k-1)
        end do
    end do
    do i = 1, 3
        do j = 1, 3
            v = 0.d0
            do k = 1, 3
                v = v+vectmp(i, k)*pgltmp(k, j)
            end do
            pgl(i, j) = v
        end do
    end do
!
end subroutine
