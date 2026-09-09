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
subroutine ef0436()
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystMemb
    implicit none
!
#include "asterc/r8dgrd.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/mbcine.h"
#include "asterfort/mbrigi.h"
#include "asterfort/ppgan2.h"
#include "asterfort/r8inir.h"
#include "asterfort/verift.h"
#include "jeveux.h"
!
! --------------------------------------------------------------------------------------------------
!
!  CALCUL DE EFGE_ELNO POUR LES MEMBRANES
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = 'RIGI'
    integer(kind=8), parameter :: nddl = 3, ncomp = 3
    integer(kind=8) :: nno, npg
    integer(kind=8) :: i, n, c, cc, kpg
    integer(kind=8) :: ivf, idfde, jgano, jvEfgeElno
    integer(kind=8) :: jvGeom, jvMaterc, jvDisp
    real(kind=8) :: dff(2, 8), vff(8), b(3, 3, 8), jac
    real(kind=8) :: epsm(3), epsthe, sig(3), sigg(3, 9), matrRigi(3, 3)
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!

! - Get plate parameters
    call getCara(plateCara, plateOrie)

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Compute global<=>local transformation
    call compCoorSystMemb(plateOrie)

    call elrefe_info(fami=fami, nno=nno, &
                     npg=npg, jvf=ivf, jdfde=idfde, jgano=jgano)

! - Input fields
    call jevech('PDEPLAR', 'L', jvDisp)
    call jevech('PMATERC', 'L', jvMaterc)

! - Output field
    call jevech('PEFFORR', 'E', jvEfgeElno)

    sigg = 0.d0
    do kpg = 1, npg
        do n = 1, nno
            vff(n) = zr(ivf+(kpg-1)*nno+n-1)
            dff(1, n) = zr(idfde+(kpg-1)*nno*2+(n-1)*2)
            dff(2, n) = zr(idfde+(kpg-1)*nno*2+(n-1)*2+1)
        end do

! ----- CALCUL DE LA MATRICE "B" :
        call mbcine(plateOrie, &
                    nno, zr(jvGeom), dff, &
                    b, jac)

!       -- CALCUL DE LA DEFORMATION MEMBRANAIRE DANS LE REPERE LOCAL
        epsm = 0.d0
        do n = 1, nno
            do i = 1, nddl
                do c = 1, ncomp
                    epsm(c) = epsm(c)+b(c, i, n)*zr(jvDisp+(n-1)*nddl+i-1)
                end do
            end do
        end do

!       -- RETRAIT DE LA DEFORMATION THERMIQUE
        call verift(fami, kpg, 1, '+', zi(jvMaterc), &
                    epsth_=epsthe)
        epsm(1) = epsm(1)-epsthe
        epsm(2) = epsm(2)-epsthe

!       --  CALCUL DE LA CONTRAINTE AU PG
        call mbrigi(fami, kpg, jvMaterc, matrRigi)
        sig = 0.d0
        do c = 1, ncomp
            do cc = 1, ncomp
                sig(c) = sig(c)+epsm(cc)*matrRigi(cc, c)
            end do
        end do
        do c = 1, ncomp
            sigg(c, kpg) = sig(c)
        end do

    end do

! - ELGA -> ELNO
    call ppgan2(jgano, 1, ncomp, sigg, zr(jvEfgeElno))

end subroutine
