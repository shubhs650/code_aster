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
subroutine matrc(plateOrie, vectBaseKpg, tempMoye, kcis, matrElas)
!
    use plate_type
    implicit none
!
#include "asterfort/coqrep.h"
#include "asterfort/jevech.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvala.h"
#include "asterfort/utbtab.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    type(plateOrie_Para), intent(in) :: plateOrie
    real(kind=8), intent(in) :: vectBaseKpg(3, 3)
    real(kind=8), intent(in) :: tempMoye, kcis
    real(kind=8), intent(out) :: matrElas(5, 5)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbPropMaxi = 5
    character(len=16) :: propName(nbPropMaxi)
    real(kind=8) :: propVale(nbPropMaxi)
    integer(kind=8) :: propCode(nbPropMaxi)
    integer(kind=8), parameter :: nbPara = 1
    character(len=8), parameter :: paraName(nbPara) = (/'TEMP'/)
    real(kind=8) :: paraVale(nbPara)
    character(len=32) :: elasKeyword
    real(kind=8) :: young, nu, nult, nutl
    real(kind=8) :: dorth(3, 3), work(3, 3), d(3, 3)
    real(kind=8) :: dcis(2, 2), d2(2, 2), el, et, glt, gtn, delta
    integer(kind=8) :: i, j, jvMaterc, nbProp
    real(kind=8) :: passag(3, 3), pas2(2, 2), c, s
!
! --------------------------------------------------------------------------------------------------
!
    matrElas = 0.d0
    paraVale(1) = tempMoye

! - Access to material
    call jevech('PMATERC', 'L', jvMaterc)
    call rccoma(zi(jvMaterc), 'ELAS', 1, elasKeyword)

    if (elasKeyword .eq. 'ELAS') then
        nbProp = 2
        propName(1) = 'E'
        propName(2) = 'NU'
        call rcvala(zi(jvMaterc), ' ', elasKeyword, &
                    nbPara, paraName, paraVale, &
                    nbProp, propName, &
                    propVale, propCode, 1)
        young = propVale(1)
        nu = propVale(2)

! ----- CONSTRUCTION DE LA MATRICE
        matrElas(1, 1) = young/(1.d0-nu*nu)
        matrElas(1, 2) = matrElas(1, 1)*nu
        matrElas(2, 1) = matrElas(1, 2)
        matrElas(2, 2) = matrElas(1, 1)
        matrElas(3, 3) = young/2.d0/(1.d0+nu)
        matrElas(4, 4) = matrElas(3, 3)*kcis
        matrElas(5, 5) = matrElas(4, 4)

    else if (elasKeyword .eq. 'ELAS_ORTH') then
        nbProp = 5
        propName(1) = 'E_L'
        propName(2) = 'E_T'
        propName(3) = 'NU_LT'
        propName(4) = 'G_LT'
        propName(5) = 'G_TN'
        call rcvala(zi(jvMaterc), ' ', elasKeyword, &
                    nbPara, paraName, paraVale, &
                    nbProp, propName, &
                    propVale, propCode, 1)
        el = propVale(1)
        et = propVale(2)
        nult = propVale(3)
        glt = propVale(4)
        gtn = propVale(5)
        nutl = et*nult/el
        delta = 1.d0-nult*nutl
        dorth(1, 1) = el/delta
        dorth(1, 2) = nult*et/delta
        dorth(1, 3) = 0.d0
        dorth(2, 2) = et/delta
        dorth(2, 1) = dorth(1, 2)
        dorth(2, 3) = 0.d0
        dorth(3, 1) = 0.d0
        dorth(3, 2) = 0.d0
        dorth(3, 3) = glt

        call coqrep(vectBaseKpg, plateOrie%alpha, plateOrie%beta, &
                    c_=c, s_=s)
        passag = 0.d0
        passag(1, 1) = c*c
        passag(2, 2) = c*c
        passag(1, 2) = s*s
        passag(2, 1) = s*s
        passag(1, 3) = c*s
        passag(3, 1) = -2.d0*c*s
        passag(2, 3) = -c*s
        passag(3, 2) = 2.d0*c*s
        passag(3, 3) = c*c-s*s
        call utbtab('ZERO', 3, 3, dorth, passag, work, d)
        do i = 1, 3
            do j = 1, 3
                matrElas(i, j) = d(i, j)
            end do
        end do
!
        pas2 = 0.d0
        pas2(1, 1) = c
        pas2(2, 2) = c
        pas2(1, 2) = s
        pas2(2, 1) = -s
        dcis(1, 1) = glt
        dcis(1, 2) = 0.d0
        dcis(2, 1) = 0.d0
        dcis(2, 2) = gtn
        call utbtab('ZERO', 2, 2, dcis, pas2, work, d2)
        do i = 1, 2
            do j = 1, 2
                matrElas(3+i, 3+j) = d2(i, j)
            end do
        end do
    else
        call utmess('F', 'ELEMENTS_45', sk=elasKeyword)
    end if
!
end subroutine
