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
! aslint: disable=W0413
!
subroutine t3grig(plateCara, plateOrie, &
                  xyzl, option, pgl, &
                  matrRigi_, ener_)
!
    use plate_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/bsthpl.h"
#include "asterfort/dstbfb.h"
#include "asterfort/dxmate.h"
#include "asterfort/dxtbm.h"
#include "asterfort/dxtloc.h"
#include "asterfort/dxtloe.h"
#include "asterfort/gtria3.h"
#include "asterfort/jevech.h"
#include "asterfort/t3gbc.h"
#include "asterfort/utbtab.h"
#include "asterfort/utmess.h"
#include "asterfort/utpvgl.h"
#include "jeveux.h"
!
    type(plateCara_Para), intent(in) :: plateCara
    type(plateOrie_Para), intent(in) :: plateOrie
    real(kind=8), intent(in) :: xyzl(3, *), pgl(*)
    character(len=16), intent(in) :: option
    real(kind=8), optional, intent(out) :: matrRigi_(300), ener_(3)
!
! --------------------------------------------------------------------------------------------------
!
!     MATRICE DE RIGIDITE DE L'ELEMENT T3GAMMA (AVEC CISAILLEMENT)
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), parameter :: zero = 0.d0
    integer(kind=8) :: multic
    real(kind=8) :: depl(18)
    real(kind=8) :: df(3, 3), dm(3, 3), dmf(3, 3), dc(2, 2), dci(2, 2)
    real(kind=8) :: dmc(3, 2), dfc(3, 2)
    real(kind=8) :: bfb(3, 9)
    real(kind=8) :: bc(2, 9)
    real(kind=8) :: bm(3, 6)
    real(kind=8) :: xab1(3, 6), xab2(3, 9), xab3(2, 9)
    real(kind=8) :: kc(81)
    real(kind=8) :: flexi(81), flex(81)
    real(kind=8) :: memb(36)
    real(kind=8) :: mefl(54)
    real(kind=8) :: enerTher, carat3(25)
    real(kind=8) :: qsi, eta
    aster_logical :: coupmf
    integer(kind=8) :: jvDisp, k
    real(kind=8) :: ctor, excent
    real(kind=8) :: aire
    real(kind=8) :: matrRigi(300), ener(3)
!
! --------------------------------------------------------------------------------------------------
!
    enerTher = zero

! - Get parameters
    ctor = plateCara%coefRigiDRZ
    excent = plateCara%offset
    if (excent .ne. zero) then
        call utmess('F', 'ELEMENTS2_57')
    end if

! - Geometric properties
    call gtria3(xyzl, carat3)

! - Get elementary matrix of rigidity
    call dxmate(plateCara, plateOrie, &
                'RIGI', df, dm, dmf, dc, &
                dci, dmc, dfc, &
                multic, coupmf)

!     ------------------------------------------------------------------
!     CALCUL DE LA MATRICE DE RIGIDITE DE L'ELEMENT EN MEMBRANE
!     ------------------------------------------------------------------
    memb = 0.d0
    flex = 0.d0
    mefl = 0.d0
    kc = 0.d0
!
!     ------ CALCUL DE LA MATRICE BM -----------------------------------
    call dxtbm(carat3(9), bm)
    aire = carat3(8)

!     ------ CALCUL DU PRODUIT BMT.DM.BM -------------------------------
    call utbtab('ZERO', 3, 6, dm, bm, &
                xab1, memb)
    do k = 1, 36
        memb(k) = memb(k)*aire
    end do
!
!     ------------------------------------------------------------------
!     CALCUL DES MATRICES DE RIGIDITE DE L'ELEMENT EN FLEXION ET
!     COUPLAGE MEMBRANE/FLEXION
!     ------------------------------------------------------------------
!     ------- CALCUL DE LA MATRICE BFB -------------------------------
    call dstbfb(carat3(9), bfb)
!
!     ------- CALCUL DU PRODUIT BFBT.DF.BFB --------------------------
    call utbtab('ZERO', 3, 9, df, bfb, &
                xab2, flex)
!
!        ---- CALCUL DE LA MATRICE BC ----------------------------------
    qsi = 1.d0/3.d0
    eta = qsi
    call t3gbc(xyzl, qsi, eta, bc)
!
!        ---- CALCUL DU PRODUIT BCT.DC.BC -----------------------------
    call utbtab('ZERO', 2, 9, dc, bc, &
                xab3, kc)
!
    do k = 1, 81
        flexi(k) = (flex(k)+kc(k))*aire
    end do
!
    if (option .eq. 'RIGI_MECA') then
        call dxtloc(flexi, memb, mefl, ctor, matrRigi)
    else if (option .eq. 'EPOT_ELEM') then
        call jevech('PDEPLAR', 'L', jvDisp)
        call utpvgl(3, 6, pgl, zr(jvDisp), depl)
        call dxtloe(flex, memb, mefl, ctor, coupmf, &
                    depl, ener)
        ! call bsthpl(plateCara, plateOrie, &
        !             jvGeom, nomte, xyzl, &
        !             bsigth)
        ! if (indith) then
        !     enerTher = 0.d0
        !     do i = 1, 18
        !         enerTher = enerTher+depl(i)*bsigth(i)
        !     end do
        !     ener(1) = ener(1)-enerTher
        ! end if
    end if
!
    if (present(matrRigi_)) then
        matrRigi_ = matrRigi
    end if
    if (present(ener_)) then
        ener_ = ener
    end if
!
end subroutine
