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
subroutine te0446(option, nomte)
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystPara, compCoorSystPlate, creaCaraMini
    use resi_refe_module, only: RESI_REFE
    implicit none
!
#include "asterc/r8dgrd.h"
#include "asterc/r8prem.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/dxbsig.h"
#include "asterfort/dxefro.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/utpvgl.h"
#include "blas/dcopy.h"
#include "jeveux.h"
!
    character(len=16) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: DKTG/Q4GG
!
! Options: FORC_NODA
!          REFE_FORC_NODA
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbEfgeNd = 8
    integer(kind=8) :: jvSief, jtab(7), jvDisp, jvGeom
    integer(kind=8) :: jvCompor, i, i1, i2, j, k, jvVect, ipg, iretc, iret
    integer(kind=8) :: nno, npg
    real(kind=8) :: pgl(3, 3), xyzl(3, 4), forcNoda(24)
    real(kind=8) :: effgt(32), effort(32)
    real(kind=8) :: effref, momref
    real(kind=8) :: foref, moref
    aster_logical :: reactu
    blas_int :: b_incx, b_incy, b_n
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
    type(RESI_REFE):: refe
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', nno=nno, npg=npg)

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

    if (option .eq. 'FORC_NODA') then
! ----- Get plate parameters
        call getCara(plateCara, plateOrie)

! ----- Calculate the transformation: global coordinate system/intrinsic coordinate system
        call compCoorSystPara(plateCara, zr(jvGeom), pgl)

! ----- Compute coordinate system for plate
        call compCoorSystPlate(pgl, plateCara, plateOrie)

! ----- Change coordinates of displacements
        call utpvgl(nno, 3, pgl, zr(jvGeom), xyzl)

! ----- VECTEUR DES EFFORTS GENERALISES
        call tecach('OOO', 'PSIEFR', 'L', iret, nval=7, itab=jtab)

! ----- PASSAGE DU VECTEUR DES EFFORTS GENERALISES DU REPERE LOCAL AU REPERE INTRINSEQUE
        do ipg = 1, npg
            jvSief = jtab(1)+8*(ipg-1)
            b_n = to_blas_int(8)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(jvSief), b_incx, effort(8*(ipg-1)+1), b_incy)
        end do
        call dxefro(npg, plateOrie%t2ui, effort, effgt)

        call tecach('ONO', 'PCOMPOR', 'L', iretc, iad=jvCompor)
        reactu = .false.
        if (iretc .eq. 0) then
            if (zk16(jvCompor+2) (6:10) .eq. '_REAC') call utmess('A', 'ELEMENTS2_72')
            reactu = (zk16(jvCompor+2) .eq. 'PETIT_REAC' .or. zk16(jvCompor+2) .eq. 'GROT_GDEP')
        end if
!
        if (reactu) then
            call jevech('PDEPLAR', 'L', jvDisp)
            do i = 1, nno
                i1 = 3*(i-1)
                i2 = 6*(i-1)
                zr(jvGeom+i1) = zr(jvGeom+i1)+zr(jvDisp+i2)
                zr(jvGeom+i1+1) = zr(jvGeom+i1+1)+zr(jvDisp+i2+1)
                zr(jvGeom+i1+2) = zr(jvGeom+i1+2)+zr(jvDisp+i2+2)
            end do

! --------- Calculate the transformation: global coordinate system/intrinsic coordinate system
            call compCoorSystPara(plateCara, zr(jvGeom), pgl)

! --------- Compute coordinate system for plate
            call compCoorSystPlate(pgl, plateCara, plateOrie)

! --------- Change coordinates of displacements
            call utpvgl(nno, 3, pgl, zr(jvGeom), xyzl)
        end if

! ----- CALCUL DES EFFORTS INTERNES (I.E. SOMME_VOL(BT_SIG))
        call dxbsig(plateCara, plateOrie, &
                    nomte, option, &
                    xyzl, pgl, effgt, &
                    forcNoda)

! ----- AFFECTATION DES VALEURS DE BSIGMA AU VECTEUR EN SORTIE
        call jevech('PVECTUR', 'E', jvVect)

        k = 0
        do i = 1, nno
            do j = 1, 6
                k = k+1
                zr(jvVect+k-1) = forcNoda(k)
            end do
        end do

    else if (option .eq. 'REFE_FORC_NODA') then
        call creaCaraMini(plateCara)
        call compCoorSystPara(plateCara, zr(jvGeom), pgl)
        call refe%Init(nomte)
        foref = refe%GetRef('EFFORT')
        moref = refe%GetRef('MOMENT')
        call refe%Check()
        do i = 1, nno
            do j = 1, 3
                effgt((i-1)*nbEfgeNd+j) = foref
                effgt((i-1)*nbEfgeNd+3+j) = moref
                effgt((i-1)*nbEfgeNd+7) = foref
                effgt((i-1)*nbEfgeNd+8) = foref
            end do
        end do

! ----- Change coordinates of displacements
        call utpvgl(nno, 3, pgl, zr(jvGeom), xyzl)

! ----- CALCUL DES EFFORTS INTERNES (I.E. SOMME_VOL(BT_SIG))
        call dxbsig(plateCara, plateOrie, &
                    nomte, option, &
                    xyzl, pgl, effgt, &
                    forcNoda)

! ----- AFFECTATION DES VALEURS DE BSIGMA AU VECTEUR EN SORTIE
        call jevech('PVECTUR', 'E', jvVect)
        k = 0
        do i = 1, nno
            effref = (abs(forcNoda(k+1))+abs(forcNoda(k+2))+abs(forcNoda(k+3)))/3.d0
            momref = (abs(forcNoda(k+4))+abs(forcNoda(k+5))+abs(forcNoda(k+6)))/3.d0
            ASSERT(abs(effref) .gt. r8prem())
            ASSERT(abs(momref) .gt. r8prem())
            do j = 1, 6
                k = k+1
                if (j .lt. 4) then
                    zr(jvVect+k-1) = effref
                else
                    zr(jvVect+k-1) = momref
                end if
            end do
        end do
    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
