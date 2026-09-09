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
subroutine te0012(option, nomte)
!
    use THM_type
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/thmGetGene.h"
#include "asterfort/thmGetElemModel.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/pmavec.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/vecma.h"
#include "blas/ddot.h"
!
    character(len=16) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
!     BUT: CALCUL DES MATRICES DE MASSE ELEMENTAIRES EN MECANIQUE
!          ELEMENTS ISOPARAMETRIQUES 3D
!
!          OPTION : 'MASS_MECA'
!          OPTION : 'MASS_MECA_DIAG'
!          OPTION : 'MASS_MECA_EXPLI'
!          OPTION : 'M_GAMMA'
!          OPTION : 'ECIN_ELEM'
!
!     ENTREES  ---> OPTION : OPTION DE CALCUL
!              ---> NOMTE  : NOM DU TYPE ELEMENT
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: icodre(1)
    character(len=16) :: phenom
    character(len=8) :: fami
    character(len=3) :: stopz
    real(kind=8) :: a(3, 3, 27, 27), matp(81, 81), matv(3321)
    real(kind=8) :: poids, rho(1)
    integer(kind=8) :: ipoids, ivf, idfde, igeom, imate
    integer(kind=8) :: nno, kp, i, j, k, imatuu, iacce, ivect
    integer(kind=8) :: ijkl, ik, l, npg, nddl, nvec
    integer(kind=8) :: n1, n2, i2, j2, k2
    integer(kind=8) :: idiag, nnos, iret
    integer(kind=8) :: idepl, ivite, iecin, ifreq
    real(kind=8) :: trace, alfa, wgt, masvit(81), masdep(81)
    real(kind=8) :: vect1(81), vect2(81)
    integer(kind=8) :: mecani(5), press1(7), press2(7), tempe(5), second(5), idec
    aster_logical, parameter :: l_vf = ASTER_FALSE
    type(THM_DS) :: ds_thm
    blas_int :: b_incx, b_incy, b_n
!
! --------------------------------------------------------------------------------------------------
!
    fami = 'MASS'
    call elrefe_info(fami=fami, nno=nno, nnos=nnos, npg=npg, jpoids=ipoids, &
                     jvf=ivf, jdfde=idfde)
    nddl = 3*nno
    nvec = nddl*(nddl+1)/2
!
! - Get model of finite element
!
    call thmGetElemModel(ds_thm)
!
! - Get generalized coordinates
!
    call thmGetGene(ds_thm, l_vf, 3, mecani, press1, &
                    press2, tempe, second)
!
    if (lteatt('TYPMOD2', 'THM')) then
        idec = press1(1)+press2(1)+tempe(1)
    else
        idec = 0
    end if
!
    call jevech('PGEOMER', 'L', igeom)
    call jevech('PMATERC', 'L', imate)
    call rccoma(zi(imate), 'ELAS', 1, phenom, icodre(1))
!
!
    a(:, :, :, :) = 0.0d0
!
!    BOUCLE SUR LES POINTS DE GAUSS
!
    do kp = 1, npg
        l = (kp-1)*nno
        call dfdm3d(nno, kp, ipoids, idfde, zr(igeom), &
                    poids)
        call rcvalb(fami, kp, 1, '+', zi(imate), &
                    ' ', phenom, 0, ' ', [0.d0], &
                    1, 'RHO', rho, icodre(1), 1)
        do i = 1, nno
            do j = 1, i
                a(1, 1, i, j) = a(1, 1, i, j)+rho(1)*poids*zr(ivf+l+i-1)*zr(ivf+l+j-1)
            end do
        end do
    end do
!
    do i = 1, nno
        do j = 1, i
            a(2, 2, i, j) = a(1, 1, i, j)
            a(3, 3, i, j) = a(1, 1, i, j)
        end do
    end do
!
    if (option .eq. 'MASS_MECA') then
!
        call jevech('PMATUUR', 'E', imatuu)
        do k = 1, nvec
            matv(k) = 0.0d0
        end do
!
! PASSAGE DU STOCKAGE RECTANGULAIRE (A) AU STOCKAGE TRIANGULAIRE (ZR)
!
        do k = 1, 3
            do l = 1, 3
                do i = 1, nno
                    ik = ((3*i+k-4)*(3*i+k-3))/2
                    do j = 1, i
                        ijkl = ik+3*(j-1)+l
                        matv(ijkl) = a(k, l, i, j)
                    end do
                end do
            end do
        end do
        if (idec .eq. 0) then
            do i = 1, nvec
                zr(imatuu+i-1) = matv(i)
            end do
        else
            do k = 1, nno
                do n1 = 1, 3
                    i = 3*k+n1-3
                    if (k .le. nnos) then
                        i2 = i+idec*(k-1)
                    else
                        i2 = i+idec*nnos
                    end if
                    do l = 1, nno
                        do n2 = 1, 3
                            j = 3*l+n2-3
                            if (j .gt. i) goto 405
                            if (l .le. nnos) then
                                j2 = j+idec*(l-1)
                            else
                                j2 = j+idec*nnos
                            end if
                            zr(imatuu+i2*(i2-1)/2+j2-1) = matv(i*(i-1)/2+j)
                        end do
                    end do
405                 continue
                end do
            end do
        end if
!
    else if (option .eq. 'MASS_MECA_DIAG' .or. option .eq. 'MASS_MECA_EXPLI') then
!
        call jevech('PMATUUR', 'E', imatuu)
!
!-- CALCUL DE LA MASSE DE L'ELEMENT
!
        wgt = a(1, 1, 1, 1)
        do i = 2, nno
            do j = 1, i-1
                wgt = wgt+2*a(1, 1, i, j)
            end do
            wgt = wgt+a(1, 1, i, i)
        end do
!
!-- CALCUL DE LA TRACE EN TRANSLATION SUIVANT X
!
        trace = 0.d0
        do i = 1, nno
            trace = trace+a(1, 1, i, i)
        end do
!
!-- CALCUL DU FACTEUR DE DIAGONALISATION
!
        alfa = wgt/trace
!
! PASSAGE DU STOCKAGE RECTANGULAIRE (A) AU STOCKAGE TRIANGULAIRE (ZR)
!
        k = 0
        do j = 1, nno
            do i = 1, 3
                k = k+1
                if (idec .eq. 0) then
                    idiag = k*(k+1)/2
                else
                    if (j .le. nnos) then
                        k2 = k+idec*(j-1)
                    else
                        k2 = k+idec*nnos
                    end if
                    idiag = k2*(k2+1)/2
                end if
                zr(imatuu+idiag-1) = a(i, i, j, j)*alfa
            end do
        end do
!
    else if (option .eq. 'M_GAMMA') then
!
        call jevech('PACCELR', 'L', iacce)
        call jevech('PVECTUR', 'E', ivect)
        do k = 1, nvec
            matv(k) = 0.0d0
        end do
        do k = 1, 3
            do l = 1, 3
                do i = 1, nno
                    ik = ((3*i+k-4)*(3*i+k-3))/2
                    do j = 1, i
                        ijkl = ik+3*(j-1)+l
                        matv(ijkl) = a(k, l, i, j)
                    end do
                end do
            end do
        end do
        call vecma(matv, nvec, matp, nddl)
        if (idec .eq. 0) then
            call pmavec('ZERO', nddl, matp, zr(iacce), zr(ivect))
        else
            do k = 1, nddl
                vect1(k) = 0.0d0
                vect2(k) = 0.0d0
            end do
            do k = 1, nno
                do n1 = 1, 3
                    i = 3*k+n1-3
                    if (k .le. nnos) then
                        i2 = i+idec*(k-1)
                    else
                        i2 = i+idec*nnos
                    end if
                    vect1(i) = zr(iacce+i2-1)
                end do
            end do
            call pmavec('ZERO', nddl, matp, vect1, vect2)
            do k = 1, nno
                do n1 = 1, 3
                    i = 3*k+n1-3
                    if (k .le. nnos) then
                        i2 = i+idec*(k-1)
                    else
                        i2 = i+idec*nnos
                    end if
                    zr(ivect+i2-1) = vect2(i)
                end do
            end do
        end if
!
! OPTION ECIN_ELEM : CALCUL DE L'ENERGIE CINETIQUE
!
    else if (option .eq. 'ECIN_ELEM') then
        stopz = 'ONO'
        call tecach(stopz, 'PVITESR', 'L', iret, iad=ivite)
! IRET NE PEUT VALOIR QUE 0 (TOUT EST OK) OU 2 (CHAMP NON FOURNI)
        if (iret .eq. 0) then
            call jevech('PENERCR', 'E', iecin)
            do k = 1, nvec
                matv(k) = 0.0d0
            end do
            do k = 1, 3
                do l = 1, 3
                    do i = 1, nno
                        ik = ((3*i+k-4)*(3*i+k-3))/2
                        do j = 1, i
                            ijkl = ik+3*(j-1)+l
                            matv(ijkl) = a(k, l, i, j)
                        end do
                    end do
                end do
            end do
            call vecma(matv, nvec, matp, nddl)
            call pmavec('ZERO', nddl, matp, zr(ivite), masvit)
            b_n = to_blas_int(nddl)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            zr(iecin) = .5d0*ddot(b_n, zr(ivite), b_incx, masvit, b_incy)
        else
            call tecach(stopz, 'PDEPLAR', 'L', iret, iad=idepl)
            if (iret .eq. 0) then
                call jevech('POMEGA2', 'L', ifreq)
                call jevech('PENERCR', 'E', iecin)
                do k = 1, nvec
                    matv(k) = 0.0d0
                end do
                do k = 1, 3
                    do l = 1, 3
                        do i = 1, nno
                            ik = ((3*i+k-4)*(3*i+k-3))/2
                            do j = 1, i
                                ijkl = ik+3*(j-1)+l
                                matv(ijkl) = a(k, l, i, j)
                            end do
                        end do
                    end do
                end do
                call vecma(matv, nvec, matp, nddl)
                call pmavec('ZERO', nddl, matp, zr(idepl), masdep)
                b_n = to_blas_int(nddl)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                zr(iecin) = .5d0*ddot(b_n, zr(idepl), b_incx, masdep, b_incy)*zr(ifreq)
            else
                call utmess('F', 'ELEMENTS2_1', sk=option)
            end if
        end if
!
    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
