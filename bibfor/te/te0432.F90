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
subroutine te0432(option, nomte)
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystPara, compCoorSystGrid
    implicit none
!
#include "asterf_types.h"
#include "asterfort/cargri.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/nmgrib.h"
#include "asterfort/pmavec.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/vecma.h"
#include "blas/ddot.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: GRILLE_MEMBRANE / GRILLE_EXCENTRE
!
! Elements: MASS_MECA* and M_GAMMA
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: propCode(2)
    character(len=8), parameter :: fami = 'MASS'
    character(len=3) :: stopz
    integer(kind=8) :: nno, npg, i, imatuu, ndim, nnos, jgano
    integer(kind=8) :: ipoids, ivf, idfde, jvGeom, jvMaterc
    integer(kind=8) :: iret, iretd, iretv
    integer(kind=8) :: kpg, n, j, kkd, m, k
    integer(kind=8) :: kk, nddl
    integer(kind=8) :: iacce, ivect, l, nvec, ivite, ifreq, iecin, idepl
    real(kind=8) :: dff(2, 8), p(3, 6), tempRefe
    real(kind=8) :: dir11(3), vff(8), b(6, 8), jac, rho(1)
    real(kind=8) :: densit
    real(kind=8) :: distn, pgl(3, 3), masdep(48)
    real(kind=8) :: aexc(3, 3, 8, 8), a(6, 6, 8, 8), coef, matv(1176)
    real(kind=8) :: matp(48, 48)
    real(kind=8) :: diag(3, 8), wgt, alfam(3), somme(3), masvit(48), ecin
    aster_logical :: lexc, ldiag
    blas_int :: b_incx, b_incy, b_n
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!
    lexc = (lteatt('MODELI', 'GRC'))
    ldiag = (option(1:10) .eq. 'MASS_MECA_')

! - FONCTIONS DE FORMES ET POINTS DE GAUSS
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)

! - Get reference temperature
    call rcvarc(' ', 'TEMP', 'REF', fami, 1, 1, tempRefe, iret)

! - Get plate parameters
    call getCara(plateCara, plateOrie)

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Calculate the transformation: global coordinate system/intrinsic coordinate system
    if (lexc) then
        call compCoorSystPara(plateCara, zr(jvGeom), pgl)
        nddl = 6
    else
        nddl = 3
    end if
    call compCoorSystGrid(pgl, plateCara, plateOrie)

! - Get input fields
    call jevech('PMATERC', 'L', jvMaterc)

    if (option .eq. 'MASS_MECA') then

    else if (option .eq. 'M_GAMMA') then
        call jevech('PACCELR', 'L', iacce)
    else if (option .eq. 'ECIN_ELEM') then
        stopz = 'ONO'
        call tecach(stopz, 'PVITESR', 'L', iretv, iad=ivite)
        if (iretv .ne. 0) then
            call tecach(stopz, 'PDEPLAR', 'L', iretd, iad=idepl)
            if (iretd .eq. 0) then
                call jevech('POMEGA2', 'L', ifreq)
            else
                call utmess('F', 'ELEMENTS2_1', sk=option)
            end if
        end if
    end if

! - Get output fields
    if (option(1:9) .eq. 'MASS_MECA') then
        call jevech('PMATUUR', 'E', imatuu)
    else if (option .eq. 'M_GAMMA') then
        call jevech('PVECTUR', 'E', ivect)
    else if (option .eq. 'ECIN_ELEM') then
        call jevech('PENERCR', 'E', iecin)
    end if

! - LECTURE DES CARACTERISTIQUES DE GRILLE ET CALCUL DE LA DIRECTION D'ARMATURE
    call cargri(plateCara, plateOrie, &
                densit, distn, dir11)

    wgt = 0.d0
    a = 0.D0
    aexc = 0.d0
    do kpg = 1, npg
        do n = 1, nno
            vff(n) = zr(ivf+(kpg-1)*nno+n-1)
            dff(1, n) = zr(idfde+(kpg-1)*nno*2+(n-1)*2)
            dff(2, n) = zr(idfde+(kpg-1)*nno*2+(n-1)*2+1)
        end do

! ----- MASS_MECA
        call rcvalb(fami, kpg, 1, '+', zi(jvMaterc), &
                    ' ', 'ELAS', 0, ' ', [0.d0], &
                    1, 'RHO', rho, propCode, 1)

! ----- CALCUL DE LA MATRICE "B" : DEPL NODAL --> EPS11 ET DU JACOBIEN
        call nmgrib(nno, zr(jvGeom), dff, dir11, lexc, &
                    plateOrie%gridNorm, b, jac, p)
        wgt = wgt+rho(1)*zr(ipoids+kpg-1)*jac*densit

!
        do n = 1, nno
            do i = 1, n
                coef = rho(1)*zr(ipoids+kpg-1)*jac*densit*vff(n)*vff(i)
                a(1, 1, n, i) = a(1, 1, n, i)+coef
                a(2, 2, n, i) = a(2, 2, n, i)+coef
                a(3, 3, n, i) = a(3, 3, n, i)+coef
            end do
        end do
!
        if (lexc) then
            do i = 1, 3
                do j = 1, 3
                    do n = 1, nno
                        do m = 1, n
                            aexc(i, j, n, m) = a(i, j, n, m)
                        end do
                    end do
                end do
            end do

            a = 0.d0
            do i = 1, 6
                do j = 1, 6
                    do n = 1, nno
                        do m = 1, n
                            do k = 1, 3
                                a(i, j, n, m) = a(i, j, n, m)+p(k, i)*p(k, j)*aexc(k, k, n, m)
                            end do
                        end do
                    end do
                end do
            end do
        end if
    end do

    if (ldiag) then
!-- CALCUL DE LA TRACE EN TRANSLATION SUIVANT X
        diag = 0.d0
        somme = 0.d0
        do i = 1, 3
            do j = 1, nno
                somme(i) = somme(i)+a(i, i, j, j)
            end do
            alfam(i) = wgt/somme(i)
        end do

! PASSAGE DU STOCKAGE RECTANGULAIRE (A) AU STOCKAGE TRIANGULAIRE (ZR)

        do j = 1, nno
            do i = 1, 3
                diag(i, j) = a(i, i, j, j)*alfam(i)
            end do
        end do
!
        do k = 1, nddl
            do l = 1, nddl
                do i = 1, nno
                    do j = 1, nno
                        a(k, l, i, j) = 0.d0
                    end do
                end do
            end do
        end do
        do k = 1, 3
            do i = 1, nno
                a(k, k, i, i) = diag(k, i)
            end do
        end do
        if (nddl .eq. 6) then
            do i = 1, nno
                a(4, 4, i, i) = a(4, 4, i, i)*alfam(1)
                a(5, 5, i, i) = a(4, 4, i, i)*alfam(2)
                a(6, 6, i, i) = a(4, 4, i, i)*alfam(3)
            end do
        end if
    end if
!
    if (option(1:9) .eq. 'MASS_MECA') then
        do k = 1, nddl
            do l = 1, nddl
                do i = 1, nno
                    kkd = ((nddl*(i-1)+k-1)*(nddl*(i-1)+k))/2
                    do j = 1, i
                        kk = kkd+nddl*(j-1)+l
                        zr(imatuu+kk-1) = a(k, l, i, j)
                    end do
                end do
            end do
        end do
    else if (option .eq. 'M_GAMMA' .or. option .eq. 'ECIN_ELEM') then
        nvec = nddl*nno*(nddl*nno+1)/2
        do k = 1, nvec
            matv(k) = 0.0d0
        end do
        do k = 1, nddl
            do l = 1, nddl
                do i = 1, nno
                    kkd = ((nddl*(i-1)+k-1)*(nddl*(i-1)+k))/2
                    do j = 1, i
                        kk = kkd+nddl*(j-1)+l
                        matv(kk) = a(k, l, i, j)
                    end do
                end do
            end do
        end do
        call vecma(matv, nvec, matp, nddl*nno)
        if (option .eq. 'M_GAMMA') then
            call pmavec('ZERO', nddl*nno, matp, zr(iacce), zr(ivect))
        else if (option .eq. 'ECIN_ELEM') then
            if (iretv .eq. 0) then
                call pmavec('ZERO', nddl*nno, matp, zr(ivite), masvit)
                b_n = to_blas_int(nddl*nno)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                ecin = .5d0*ddot(b_n, zr(ivite), b_incx, masvit, b_incy)
            else
                call pmavec('ZERO', nddl*nno, matp, zr(idepl), masdep)
                b_n = to_blas_int(nddl*nno)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                ecin = .5d0*ddot(b_n, zr(idepl), b_incx, masdep, b_incy)*zr(ifreq)
            end if
            zr(iecin) = ecin
        end if
    end if
!
end subroutine
