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
subroutine te0055(option, nomte)
!
    use plate_type
    use plateGeom_module, only: creaCaraMini
    use resi_refe_module, only: RESI_REFE
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/dxbsig.h"
#include "asterfort/dxqpgl.h"
#include "asterfort/dxtpgl.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/utmess.h"
#include "asterfort/utpvgl.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: DKT
!
! Options: REFE_FORC_NODA
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = 'RIGI'
    integer(kind=8), parameter :: nbEfgeNd = 8
    integer(kind=8) :: nno
    integer(kind=8) :: i, j, k
    integer(kind=8) :: jvGeom, jvVect
    real(kind=8) :: pgl(3, 3), xyzl(3, 4), bsigmEner(24), effgt(32)
    real(kind=8) :: effref, momref
    real(kind=8) :: foref, moref
    type(RESI_REFE) :: refe
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami=fami, nno=nno)
    ASSERT(nno .eq. 3 .or. nno .eq. 4)
    ASSERT(option .eq. 'REFE_FORC_NODA')

    call creaCaraMini(plateCara)

! - Get geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Compute matrix for local basis
    if (nno .eq. 3) then
        call dxtpgl(zr(jvGeom), pgl)
    else if (nno .eq. 4) then
        call dxqpgl(zr(jvGeom), pgl)
    end if

! - Change coordinates of geometry
    call utpvgl(nno, 3, pgl, zr(jvGeom), xyzl)

    call refe%Init(nomte)
    foref = refe%GetRef('EFFORT')
    moref = refe%GetRef('MOMENT')
    call refe%Check()
    do i = 1, nno
        do j = 1, 3
            effgt(nbEfgeNd*(i-1)+j) = foref
            effgt(nbEfgeNd*(i-1)+3+j) = moref
        end do
        effgt(nbEfgeNd*(i-1)+7) = 0.0d0
        effgt(nbEfgeNd*(i-1)+8) = 0.0d0
    end do

! - CALCUL DES EFFORTS INTERNES (I.E. SOMME_VOL(BT_SIG))
    call dxbsig(plateCara, plateOrie, &
                nomte, option, &
                xyzl, pgl, effgt, &
                bsigmEner)

! - AFFECTATION DES VALEURS DE BSIGMA AU VECTEUR EN SORTIE
    call jevech('PVECTUR', 'E', jvVect)
    k = 0
    do i = 1, nno
        effref = (abs(bsigmEner(k+1))+abs(bsigmEner(k+2))+abs(bsigmEner(k+3)))/3.d0
        momref = (abs(bsigmEner(k+4))+abs(bsigmEner(k+5))+abs(bsigmEner(k+6)))/3.d0
        do j = 1, 6
            k = k+1
            if (j .lt. 4) then
                zr(jvVect+k-1) = effref
            else
                zr(jvVect+k-1) = momref
            end if
        end do
    end do
!
end subroutine
