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
! aslint: disable=C0110
!
subroutine te0443(option, nomte)
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystCO3D
    implicit none
!
#include "asterc/r8dgrd.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/vdefro.h"
#include "asterfort/vdrep2.h"
#include "asterfort/vdsiro.h"
#include "asterfort/vectgt.h"
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
! Options: REPE_TENS, REPE_GENE
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), parameter :: zero = 0.d0
    integer(kind=8), parameter :: nptmax = 9, ncpmax = 8, nspmax = 162
    integer(kind=8) :: nno, npg, iret(4)
    integer(kind=8) :: jvGeom, jvFieldIn, jvFieldOut, jvAngrep, np, i, itab(7), iret1, iret2, nbsp
    integer(kind=8) :: vali(2)
    integer(kind=8) :: intsn, j, k, lzi, lzr, nb1
    integer(kind=8) :: nb2, ncmp, npgsn, npgsr
    real(kind=8) :: epais, s
    real(kind=8) :: matvn1(2, 2, 10), matvg1(2, 2, 10)
    real(kind=8) :: matvn2(2, 2, 10), matvg2(2, 2, 10)
    real(kind=8) :: vectBaseKpg(3, 3), fieldInLoca(nptmax*ncpmax*nspmax)
    integer(kind=8) :: repType
    real(kind=8) :: repAlpha, repBeta
    character(len=8) :: paraInName, paraOutName
    character(len=24) :: messk(2)
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(option .eq. 'REPE_TENS' .or. option .eq. 'REPE_GENE')

! - Get plate parameters
    call getCara(plateCara, plateOrie)

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Compute global<=>local transformation
    call compCoorSystCO3D(nomte, jvGeom, &
                          plateCara, plateOrie)

! - Access to static objects of COQUE_3D
    call jevete('&INEL.'//nomte(1:8)//'.DESI', ' ', lzi)
    nb1 = zi(lzi-1+1)
    nb2 = zi(lzi-1+2)
    npgsr = zi(lzi-1+3)
    npgsn = zi(lzi-1+4)
    call jevete('&INEL.'//nomte(1:8)//'.DESR', ' ', lzr)

    if (option .eq. 'REPE_TENS') then
        ncmp = 6
        call tecach('ONO', 'PCOGAIN', 'L', iret(1), nval=7, itab=itab)
        call tecach('ONO', 'PCONOIN', 'L', iret(2), nval=7, itab=itab)
        call tecach('ONO', 'PDEGAIN', 'L', iret(3), nval=7, itab=itab)
        call tecach('ONO', 'PDENOIN', 'L', iret(4), nval=7, itab=itab)
        iret1 = iret(1)+iret(2)+iret(3)+iret(4)
        ASSERT(iret1 .eq. 6)
        if (iret(1) .eq. 0) then
            paraInName = 'PCOGAIN'
            paraOutName = 'PCOGAOUT'
        else if (iret(2) .eq. 0) then
            paraInName = 'PCONOIN'
            paraOutName = 'PCONOOUT'
        else if (iret(3) .eq. 0) then
            paraInName = 'PDEGAIN'
            paraOutName = 'PDEGAOUT'
        else if (iret(4) .eq. 0) then
            paraInName = 'PDENOIN'
            paraOutName = 'PDENOOUT'
        end if
    else if (option .eq. 'REPE_GENE') then
        ncmp = 8
        call tecach('ONO', 'PEFGAIN', 'L', iret(1), nval=7, itab=itab)
        call tecach('ONO', 'PEFNOIN', 'L', iret(2), nval=7, itab=itab)
        call tecach('ONO', 'PDGGAIN', 'L', iret(3), nval=7, itab=itab)
        call tecach('ONO', 'PDGNOIN', 'L', iret(4), nval=7, itab=itab)
        iret1 = iret(1)+iret(2)+iret(3)+iret(4)
        ASSERT(iret1 .eq. 6)
        if (iret(1) .eq. 0) then
            paraInName = 'PEFGAIN'
            paraOutName = 'PEFGAOUT'
        else if (iret(2) .eq. 0) then
            paraInName = 'PEFNOIN'
            paraOutName = 'PEFNOOUT'
        else if (iret(3) .eq. 0) then
            paraInName = 'PDGGAIN'
            paraOutName = 'PDGGAOUT'
        else if (iret(4) .eq. 0) then
            paraInName = 'PDGNOIN'
            paraOutName = 'PDGNOOUT'
        end if
    end if

    call elrefe_info(fami='MASS', nno=nno, npg=npg)
    if (paraInName(4:5) .eq. 'NO') then
        np = nno
    else if (paraInName(4:5) .eq. 'GA') then
        np = npg
    end if
    ASSERT(np .le. nptmax)

! - Get input and output fields
    call jevech(paraInName, 'L', jvFieldIn)
    call tecach('OOO', paraInName, 'L', iret2, nval=7, itab=itab)
    nbsp = itab(7)
    if ((nbsp .ne. 1) .and. (mod(nbsp, 3) .ne. 0)) then
        call utmess('F', 'ELEMENTS5_54', si=nbsp)
    end if
    ASSERT(ncmp .le. ncpmax)
    vali(1) = nspmax
    vali(2) = nbsp
    if (nbsp .gt. nspmax) then
        call utmess('F', 'ELEMENTS5_4', ni=2, vali=vali)
    end if
    call jevech(paraOutName, 'E', jvFieldOut)

! - Paramètres de la coque
    epais = plateCara%thick

! - Definition of new frame for output field
    call jevech('PANGREP', 'L', jvAngrep)
    repType = nint(zr(jvAngrep-1+3))
    repAlpha = zr(jvAngrep-1+1)
    repBeta = zr(jvAngrep-1+2)

! - Compute local basis at middle plane
    k = 0
    do intsn = 1, npgsn
        call vectgt(plateOrie, 1, nb1, &
                    zr(jvGeom), zero, intsn, &
                    epais, zr(lzr), &
                    vectBaseKpg)
        do j = 1, 3
            do i = 1, 3
                k = k+1
                zr(lzr+2000+k-1) = vectBaseKpg(i, j)
            end do
        end do
    end do

! - Get matrices for local => global from input field
    matvn1 = plateOrie%matevn
    matvg1 = plateOrie%matevg

! - Construct matrices for global => local from input field
    if ((repType .eq. 0) .or. (repType .eq. 2)) then
        if (paraInName(4:5) .eq. 'NO') then
            do i = 1, np
                s = matvn1(1, 2, i)
                matvn1(2, 1, i) = s
                matvn1(1, 2, i) = -s
            end do
        else if (paraInName(4:5) .eq. 'GA') then
            do i = 1, np
                s = matvg1(1, 2, i)
                matvg1(2, 1, i) = s
                matvg1(1, 2, i) = -s
            end do
        end if
    end if
!
    if (repType == 0) then
! ----- Local to global with inv(matvn1/matvg1): global=> local
        if (option .eq. 'REPE_TENS') then
            if (paraInName(4:5) .eq. 'NO') then
                call vdsiro(np, nbsp, matvn1, 'IU', 'N', zr(jvFieldIn), fieldInLoca)
            else if (paraInName(4:5) .eq. 'GA') then
                call vdsiro(np, nbsp, matvg1, 'IU', 'G', zr(jvFieldIn), fieldInLoca)
            end if
        else if (option .eq. 'REPE_GENE') then
            if (paraInName(4:5) .eq. 'NO') then
                call vdefro(np, matvn1, zr(jvFieldIn), fieldInLoca)
            else if (paraInName(4:5) .eq. 'GA') then
                call vdefro(np, matvg1, zr(jvFieldIn), fieldInLoca)
            end if
        end if

! ----- Compute matrix from old to new frame
        call vdrep2(repAlpha, repBeta, nb2, npgsr, zr(lzr), matvn2, matvg2)

! ----- Local Old => local New
        if (option .eq. 'REPE_TENS') then
            if (paraInName(4:5) .eq. 'NO') then
                call vdsiro(np, nbsp, matvn2, 'IU', 'N', fieldInLoca, zr(jvFieldOut))
            else if (paraInName(4:5) .eq. 'GA') then
                call vdsiro(np, nbsp, matvg2, 'IU', 'G', fieldInLoca, zr(jvFieldOut))
            end if
        else if (option .eq. 'REPE_GENE') then
            if (paraInName(4:5) .eq. 'NO') then
                call vdefro(np, matvn2, fieldInLoca, zr(jvFieldOut))
            else if (paraInName(4:5) .eq. 'GA') then
                call vdefro(np, matvg2, fieldInLoca, zr(jvFieldOut))
            end if
        end if

    else if (repType == 1) then
! ----- Local to global with matvn1/matvg1
        if (option .eq. 'REPE_TENS') then
            if (paraInName(4:5) .eq. 'NO') then
                call vdsiro(np, nbsp, matvn1, 'IU', 'N', zr(jvFieldIn), zr(jvFieldOut))
            else if (paraInName(4:5) .eq. 'GA') then
                call vdsiro(np, nbsp, matvg1, 'IU', 'G', zr(jvFieldIn), zr(jvFieldOut))
            end if
        else if (option .eq. 'REPE_GENE') then
            if (paraInName(4:5) .eq. 'NO') then
                call vdefro(np, matvn1, zr(jvFieldIn), zr(jvFieldOut))
            else if (paraInName(4:5) .eq. 'GA') then
                call vdefro(np, matvg1, zr(jvFieldIn), zr(jvFieldOut))
            end if
        end if

    else if (repType == 2) then
! ----- Local to global with inv(matvn1/matvg1): global=> local
        if (option .eq. 'REPE_TENS') then
            if (paraInName(4:5) .eq. 'NO') then
                call vdsiro(np, nbsp, matvn1, 'IU', 'N', zr(jvFieldIn), zr(jvFieldOut))
            else if (paraInName(4:5) .eq. 'GA') then
                call vdsiro(np, nbsp, matvg1, 'IU', 'G', zr(jvFieldIn), zr(jvFieldOut))
            end if
        else if (option .eq. 'REPE_GENE') then
            if (paraInName(4:5) .eq. 'NO') then
                call vdefro(np, matvn1, zr(jvFieldIn), zr(jvFieldOut))
            else if (paraInName(4:5) .eq. 'GA') then
                call vdefro(np, matvg1, zr(jvFieldIn), zr(jvFieldOut))
            end if
        end if
    else
        messk(1) = nomte
        messk(2) = 'COQUE_UTIL_CYL'
        call utmess('F', 'ALGORITH12_41', nk=2, valk=messk)
    end if
!
end subroutine
