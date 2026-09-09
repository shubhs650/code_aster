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
subroutine te0486(option, nomte)
!
    use plateGeom_module, only: getManifoldBase
    implicit none
!
#include "asterc/r8prem.h"
#include "asterf_types.h"
#include "asterfort/antisy.h"
#include "asterfort/assert.h"
#include "asterfort/b1tdb2.h"
#include "asterfort/btsig.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/matdn.h"
#include "asterfort/provec.h"
#include "asterfort/tecach.h"
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: COQUE_3D
!
! Options: CHAR_MECA_PRSU_R
!          RIGI_MECA_PRSU_R
!          CHAR_MECA_PRSU_F
!          RIGI_MECA_PRSU_F
!          CHAR_MECA_SRCO3D
!          RIGI_MECA_SRCO3D
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbPara = 4
    character(len=8), parameter :: paraName(nbPara) = (/"X   ", "Y   ", "Z   ", "INST"/)
    real(kind=8) :: paraVale(nbPara)
    integer(kind=8) :: nno
    integer(kind=8) :: jvGeom, jvDispM, jvDispP, jvPres, jpres
    integer(kind=8) :: lzi, lzr, iadzi, iazk24, iret
    integer(kind=8) :: i, j, in, kn, ii, komptn, nb1, nb2
    integer(kind=8) :: jvVect, jvMatr, intsn, npgsn
    integer(kind=8) :: irco3d, ifco3d, itemps, ierz
    real(kind=8) :: presKpg, presNode(9), madn(3, 51), nks1(3, 51), nks2(3, 51)
    real(kind=8) :: a1(3), a2(3), anta1(3, 3), anta2(3, 3), surf(3)
    real(kind=8) :: matrRigi(51*51), pr
    real(kind=8) :: vecta(9, 2, 3)
    real(kind=8) :: geom_reac(3*9)
    aster_logical :: lFrameLocal
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(nomte .eq. 'MEC3QU9H' .or. nomte .eq. 'MEC3TR7H')
    call elrefe_info(fami='RIGI', nno=nno)

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Get displacements
    call jevech('PDEPLMR', 'L', jvDispM)
    call jevech('PDEPLPR', 'L', jvDispP)

! - Get input fields
    call tecach('NNO', 'PPRESSR', 'L', iret, iad=jvPres)
    if (option .eq. 'CHAR_MECA_PRSU_R' .or. option .eq. 'RIGI_MECA_PRSU_R') then
        call jevech('PPRESSR', 'L', jpres)
    else if (option .eq. 'CHAR_MECA_PRSU_F' .or. option .eq. 'RIGI_MECA_PRSU_F') then
        call jevech('PPRESSF', 'L', jpres)
    end if

! - Access to static objects of COQUE_3D
    call jevete('&INEL.'//nomte(1:8)//'.DESI', ' ', lzi)
    nb1 = zi(lzi-1+1)
    nb2 = zi(lzi-1+2)
    npgsn = zi(lzi-1+4)
    call jevete('&INEL.'//nomte(1:8)//'.DESR', ' ', lzr)

! - Get output fields
    if (option(1:16) .eq. 'CHAR_MECA_SRCO3D' .or. option(1:16) .eq. 'CHAR_MECA_SFCO3D' .or. &
        option(1:16) .eq. 'CHAR_MECA_PRSU_R' .or. option(1:16) .eq. 'CHAR_MECA_PRSU_F') then
        call jevech('PVECTUR', 'E', jvVect)
    end if
    if (option(1:16) .eq. 'RIGI_MECA_SRCO3D' .or. option(1:16) .eq. 'RIGI_MECA_SFCO3D' .or. &
        option(1:16) .eq. 'RIGI_MECA_PRSU_R' .or. option(1:16) .eq. 'RIGI_MECA_PRSU_F') then
        call jevech('PMATUNS', 'E', jvMatr)
        matrRigi = 0.d0
    end if

! - Update geometry
    do in = 1, nb2-1
        do ii = 1, 3
            geom_reac(3*(in-1)+ii) = zr(jvGeom-1+3*(in-1)+ii)+ &
                                     zr(jvDispM-1+6*(in-1)+ii)+ &
                                     zr(jvDispP-1+6*(in-1)+ii)
        end do
    end do

! - Compute pressure at nodes
    presNode = 0.d0
    if (option(10:16) .eq. '_SRCO3D') then
        call jevech('PFRCO3D', 'L', irco3d)
        lFrameLocal = abs(zr(irco3d-1+7)-3.d0) .lt. 1.d-3
        if (lFrameLocal) then
            do kn = 1, nb2
                presNode(kn) = zr(irco3d-1+(kn-1)*7+3)
            end do
        else
            do kn = 1, nb2
                do i = 1, 6
                    if (abs(zr(irco3d-1+(kn-1)*7+i)) .gt. r8prem()) then
                        call utmess('F', 'CHARGES_11')
                    end if
                end do
            end do
            presNode(1:nb2) = 0.d0
        end if

    else if (option(10:16) .eq. '_SFCO3D') then
        call jevech('PFFCO3D', 'L', ifco3d)
        call jevech('PINSTR', 'L', itemps)
        paraVale(4) = zr(itemps)
        lFrameLocal = zk8(ifco3d-1+7) .eq. 'LOCAL_PR'
        do in = 0, nb2-1
            paraVale(1) = zr(jvGeom+3*in)
            paraVale(2) = zr(jvGeom+3*in+1)
            paraVale(3) = zr(jvGeom+3*in+2)
            if (lFrameLocal) then
                call fointe('FM', zk8(ifco3d+2), nbPara, paraName, paraVale, &
                            pr, ierz)
                presNode(in+1) = pr
                if (ierz .ne. 0) then
                    call utmess('F', 'ELEMENTS4_1')
                end if
            else
                do i = 1, 6
                    call fointe('FM', zk8(ifco3d-1+i), nbPara, paraName, paraVale, &
                                pr, ierz)
                    if (abs(pr) .gt. r8prem()) call utmess('F', 'CHARGES_11')
                end do
                presNode(in+1) = 0.d0
            end if
        end do

    else if (option(10:16) .eq. '_PRSU_R') then
        do kn = 1, nb2
            presNode(kn) = zr(jpres-1+(kn-1)*1+1)
        end do

    else if (option(10:16) .eq. '_PRSU_F') then
        call jevech('PINSTR', 'L', itemps)
        paraVale(4) = zr(itemps)
        do j = 0, nb1-1
            paraVale(1) = zr(jvGeom+3*j)
            paraVale(2) = zr(jvGeom+3*j+1)
            paraVale(3) = zr(jvGeom+3*j+2)
            call fointe('FM', zk8(jpres), nbPara, paraName, paraVale, &
                        pr, ierz)
            if (ierz .ne. 0) then
                call utmess('F', 'ELEMENTS4_1')
            end if
            if (pr .ne. 0.d0) then
                call tecael(iadzi, iazk24)
                call utmess('F', 'ELEMENTS4_92', si=zi(iadzi-1+1))
            end if
        end do
    end if

! - VECTEURS TANGENTS A1 ET A2 AUX NOEUDS NON NORMALISES
    call getManifoldBase(nb1, nb2, zr(lzr), &
                         geom_reac, &
                         vecta)

    do intsn = 1, npgsn
        a1 = 0.d0
        a2 = 0.d0
        presKpg = 0.d0
        do kn = 1, nb2
            presKpg = presKpg+zr(lzr-1+459+9*(intsn-1)+kn)*presNode(kn)
            do ii = 1, 3
                a1(ii) = a1(ii)+zr(lzr-1+459+9*(intsn-1)+kn)*vecta(kn, 1, ii)
                a2(ii) = a2(ii)+zr(lzr-1+459+9*(intsn-1)+kn)*vecta(kn, 2, ii)
            end do
        end do

! ----- A1 VECTORIEL A2
        call provec(a1, a2, surf)

! ----- MATRICE D INTERPOLATION POUR LES DEPLACEMENTS
        call matdn(nb1, zr(lzr), intsn, madn, nks1, &
                   nks2)

        if (option(1:16) .eq. 'CHAR_MECA_SRCO3D' .or. option(1:16) .eq. 'CHAR_MECA_SFCO3D' .or. &
            option(1:16) .eq. 'CHAR_MECA_PRSU_R') then
            call btsig(6*nb1+3, 3, -presKpg*zr(lzr-1+127+intsn-1), madn, surf, &
                       zr(jvVect))
        end if
        if (option(1:16) .eq. 'RIGI_MECA_SRCO3D' .or. option(1:16) .eq. 'RIGI_MECA_SFCO3D' .or. &
            option(1:16) .eq. 'RIGI_MECA_PRSU_R') then

! --------- MATRICE ANTISYM DE A1 ET DE A2
            call antisy(a1, 1.d0, anta1)
            call antisy(a2, 1.d0, anta2)

! --------- PREMIER TERME
            call b1tdb2(madn, nks2, anta1, presKpg*zr(lzr-1+127+intsn-1), 3, &
                        6*nb1+3, matrRigi)

! --------- DEUXIEME TERME
            call b1tdb2(madn, nks1, anta2, -presKpg*zr(lzr-1+127+intsn-1), 3, &
                        6*nb1+3, matrRigi)
        end if
    end do
!
    if (option(1:16) .eq. 'RIGI_MECA_SRCO3D' .or. option(1:16) .eq. 'RIGI_MECA_SFCO3D' .or. &
        option(1:16) .eq. 'RIGI_MECA_PRSU_R') then
        komptn = 0
        do j = 1, 6*nb1+3
            do i = 1, 6*nb1+3
                zr(jvMatr+komptn) = -matrRigi((6*nb1+3)*(i-1)+j)
                komptn = komptn+1
            end do
        end do
    end if

!
end subroutine
