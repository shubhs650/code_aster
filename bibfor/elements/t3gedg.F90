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
subroutine t3gedg(plateCara, plateOrie, &
                  xyzl, option, depl, edgl)
!
    use plate_type
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/dstbfb.h"
#include "asterfort/dxmate.h"
#include "asterfort/dxtbm.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/gtria3.h"
#include "asterfort/t3gbc.h"
!
    type(plateCara_Para), intent(in) :: plateCara
    type(plateOrie_Para), intent(in) :: plateOrie
    real(kind=8) :: xyzl(3, *), depl(*), edgl(*)
    character(len=16) :: option
!
! --------------------------------------------------------------------------------------------------
!
! T3G/T3GG
!
! DEGE_ELGA / DEGE_ELNO
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ndim, nno, nnos, npg, ipoids, icoopg, ivf, idfdx, idfd2, jgano
    integer(kind=8) :: multic, ne, k, j, i, ie
    real(kind=8) :: depf(9), depm(6)
    real(kind=8) :: df(3, 3), dm(3, 3), dmf(3, 3), dmc(3, 2), dfc(3, 2)
    real(kind=8) :: dci(2, 2), dc(2, 2)
    real(kind=8) :: bf(3, 9), bm(3, 6), bc(2, 9)
    real(kind=8) :: bdm(3), bdf(3), bcdf(2), dcis(2)
    real(kind=8) :: vf(3), vm(3), vt(2)
    real(kind=8) :: vfm(3), vmf(3), vmc(3), vfc(3), carat3(21)
    real(kind=8) :: qsi, eta
    aster_logical :: coupmf
    character(len=8) :: fami
!
! --------------------------------------------------------------------------------------------------
!
    fami = 'RIGI'
    if (option(6:9) .eq. 'ELGA') then
        call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                         jpoids=ipoids, jcoopg=icoopg, jvf=ivf, jdfde=idfdx, jdfd2=idfd2, &
                         jgano=jgano)
        ne = npg
        fami = 'RIGI'
    else if (option(6:9) .eq. 'ELNO') then
        call elrefe_info(fami='NOEU', ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                         jpoids=ipoids, jcoopg=icoopg, jvf=ivf, jdfde=idfdx, jdfd2=idfd2, &
                         jgano=jgano)
        ne = nno
        fami = 'NOEU'
    end if
!
!     ----- CALCUL DES MATRICES DE RIGIDITE DU MATERIAU EN FLEXION,
!           MEMBRANE ET CISAILLEMENT INVERSEES -------------------------
!
!     ----- CALCUL DU JACOBIEN SUR LE TRIANGLE -----------------
    call gtria3(xyzl, carat3)
!     ----- CARACTERISTIQUES DES MATERIAUX --------
    call dxmate(plateCara, plateOrie, &
                fami, df, dm, dmf, dc, &
                dci, dmc, dfc, &
                multic, coupmf)
!     ----- COMPOSANTES DEPLACEMENT MEMBRANE ET FLEXION ----------------
    do j = 1, nno
        do i = 1, 2
            depm(i+2*(j-1)) = depl(i+6*(j-1))
        end do
        depf(1+3*(j-1)) = depl(1+2+6*(j-1))
        depf(2+3*(j-1)) = depl(3+2+6*(j-1))
        depf(3+3*(j-1)) = -depl(2+2+6*(j-1))
    end do
!

!     ------ CALCUL DE LA MATRICE BM -----------------------------------
    call dxtbm(carat3(9), bm)
!
!     ------- CALCUL DE LA MATRICE BFB -------------------------------
    call dstbfb(carat3(9), bf)
!
!     ---- CALCUL DE LA MATRICE BC ----------------------------------
!
    if (option(1:4) .eq. 'DEGE') then
        qsi = 1.d0/3.d0
        eta = qsi
        call t3gbc(xyzl, qsi, eta, bc)
        do ie = 1, ne
            bcdf(1) = 0.d0
            bcdf(2) = 0.d0
            do j = 1, 9
                bcdf(1) = bcdf(1)+bc(1, j)*depf(j)
                bcdf(2) = bcdf(2)+bc(2, j)*depf(j)
            end do
            do k = 1, 3
                bdf(k) = 0.d0
                bdm(k) = 0.d0
            end do
            do i = 1, 3
                do j = 1, 9
                    bdf(i) = bdf(i)+bf(i, j)*depf(j)
                end do
                do j = 1, 6
                    bdm(i) = bdm(i)+bm(i, j)*depm(j)
                end do
            end do
            do i = 1, 3
                edgl(i+8*(ie-1)) = bdm(i)
                edgl(i+3+8*(ie-1)) = bdf(i)
            end do
            edgl(7+8*(ie-1)) = bcdf(1)
            edgl(8+8*(ie-1)) = bcdf(2)
!           --- PASSAGE DE LA DISTORSION A LA DEFORMATION DE CIS. ------
            edgl(3+8*(ie-1)) = edgl(3+8*(ie-1))/2.d0
            edgl(6+8*(ie-1)) = edgl(6+8*(ie-1))/2.d0
            edgl(7+8*(ie-1)) = bcdf(1)/2.d0
            edgl(8+8*(ie-1)) = bcdf(2)/2.d0
        end do
    else
        do ie = 1, ne
!           ------ VT = DC.BC.DEPF -------------------------------------
!     ---- CALCUL DE LA MATRICE BC ----------------------------------
            qsi = 1.d0/3.d0
            eta = qsi
!
            call t3gbc(xyzl, qsi, eta, bc)
!
            vt(1) = 0.d0
            vt(2) = 0.d0
            bcdf(1) = 0.d0
            bcdf(2) = 0.d0
            do j = 1, 9
                bcdf(1) = bcdf(1)+bc(1, j)*depf(j)
                bcdf(2) = bcdf(2)+bc(2, j)*depf(j)
            end do
            vt(1) = dc(1, 1)*bcdf(1)+dc(1, 2)*bcdf(2)
            vt(2) = dc(2, 1)*bcdf(1)+dc(2, 2)*bcdf(2)
            do k = 1, 3
                bdf(k) = 0.d0
                bdm(k) = 0.d0
                vf(k) = 0.d0
                vm(k) = 0.d0
                vfm(k) = 0.d0
                vmf(k) = 0.d0
                vmc(k) = 0.0d0
                vfc(k) = 0.0d0
            end do
!           ------ VF = DF.BF.DEPF , VFM = DMF.BM.DEPM ----------------
!           ------ VM = DM.BM.DEPM , VMF = DMF.BF.DEPF ----------------
            do i = 1, 3
                do j = 1, 9
                    bdf(i) = bdf(i)+bf(i, j)*depf(j)
                end do
                do j = 1, 6
                    bdm(i) = bdm(i)+bm(i, j)*depm(j)
                end do
            end do
            do i = 1, 3
                do j = 1, 3
                    vf(i) = vf(i)+df(i, j)*bdf(j)
                    vfm(i) = vfm(i)+dmf(i, j)*bdm(j)
                    vm(i) = vm(i)+dm(i, j)*bdm(j)
                    vmf(i) = vmf(i)+dmf(i, j)*bdf(j)
                end do
            end do
!
            dcis(1) = dci(1, 1)*vt(1)+dci(1, 2)*vt(2)
            dcis(2) = dci(2, 1)*vt(1)+dci(2, 2)*vt(2)
!
            vmc(1) = dmc(1, 1)*dcis(1)+dmc(1, 2)*dcis(2)
            vmc(2) = dmc(2, 1)*dcis(1)+dmc(2, 2)*dcis(2)
            vmc(3) = dmc(3, 1)*dcis(1)+dmc(3, 2)*dcis(2)
!
            vfc(1) = dfc(1, 1)*dcis(1)+dfc(1, 2)*dcis(2)
            vfc(2) = dfc(2, 1)*dcis(1)+dfc(2, 2)*dcis(2)
            vfc(3) = dfc(3, 1)*dcis(1)+dfc(3, 2)*dcis(2)
!
!
            do i = 1, 3
                edgl(i+8*(ie-1)) = vm(i)+vmf(i)+vmc(i)
                edgl(i+3+8*(ie-1)) = vf(i)+vfm(i)+vfc(i)
            end do
            edgl(7+8*(ie-1)) = vt(1)
            edgl(8+8*(ie-1)) = vt(2)
        end do
    end if
end subroutine
