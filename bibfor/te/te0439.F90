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
subroutine te0439(option, nomte)
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystMemb
    implicit none
!
#include "asterc/r8dgrd.h"
#include "asterc/r8prem.h"
#include "asterf_types.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/mbcine.h"
#include "asterfort/rcvalb.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    character(len=16) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: MEMBRANE
!
! Options: MASS_MECA*
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = 'MASS'
    integer(kind=8), parameter :: nddl = 3
    integer(kind=8), parameter :: nbProp = 1
    character(len=8), parameter :: propName(nbProp) = (/'RHO'/)
    real(kind=8) :: propVale(nbProp)
    integer(kind=8) :: propCode(nbProp)
    integer(kind=8) :: nno, npg, i, jvMatr, iret_cmp
    integer(kind=8) :: ipoids, ivf, idfde, jvGeom, jvMaterc, jvCompor
    integer(kind=8) :: kpg, n, j, kkd, k
    integer(kind=8) :: kk, l
    real(kind=8) :: dff(2, 9)
    real(kind=8) :: vff(9), b(3, 3, 9), jac, rho
    real(kind=8) :: h, preten
    real(kind=8) :: a(3, 3, 9, 9), coef
    real(kind=8) :: diag(3, 9), wgt, alfam(3), somme(3)
    aster_logical :: ldiag, grdef
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

!
    call tecach('ONO', 'PCOMPOR', 'L', iret_cmp, iad=jvCompor)
    grdef = ASTER_FALSE
    if (iret_cmp == 0) then
        grdef = (zk16(jvCompor+2) (1:9) .eq. 'GROT_GDEP')
    end if
    ldiag = (option(1:10) .eq. 'MASS_MECA_')

! - FONCTIONS DE FORMES ET POINTS DE GAUSS
    call elrefe_info(fami=fami, nno=nno, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde)
    a = 0.d0

! - Input fields
    call jevech('PMATERC', 'L', jvMaterc)

! - Output field
    call jevech('PMATUUR', 'E', jvMatr)

! - EPAISSEUR ET PRETCONTRAINTES
    h = plateCara%thick
    if (h .lt. r8prem()) then
        call utmess('F', 'MEMBRANE_1')
    end if
    preten = plateCara%tension/h

    if (grdef) then
        call utmess('F', 'MEMBRANE_9')

    else
        if (iret_cmp .ne. 0) then
            if (abs(h-1.d0) .gt. r8prem()) then
                call utmess('F', 'MEMBRANE_11')
            end if
        end if
!
        h = 1.d0
    end if

    wgt = 0.d0
    do kpg = 1, npg
        do n = 1, nno
            vff(n) = zr(ivf+(kpg-1)*nno+n-1)
            dff(1, n) = zr(idfde+(kpg-1)*nno*2+(n-1)*2)
            dff(2, n) = zr(idfde+(kpg-1)*nno*2+(n-1)*2+1)
        end do
!
! - MASS_MECA
!
        if (grdef) then
            call rcvalb(fami, kpg, 1, '+', &
                        zi(jvMaterc), ' ', 'ELAS', &
                        0, ' ', [0.d0], &
                        nbProp, propName, propVale, &
                        propCode, 1)
        else
            call rcvalb(fami, kpg, 1, '+', &
                        zi(jvMaterc), ' ', 'ELAS_MEMBRANE', &
                        0, ' ', [0.d0], &
                        nbProp, propName, propVale, &
                        propCode, 1)
        end if
        rho = propVale(1)

! ----- CALCUL DE LA MATRICE "B" :
        call mbcine(plateOrie, &
                    nno, zr(jvGeom), dff, &
                    b, jac)
!
        wgt = wgt+rho*zr(ipoids+kpg-1)*jac*h
        do n = 1, nno
            do i = 1, n
                coef = rho*zr(ipoids+kpg-1)*jac*vff(n)*vff(i)*h
                a(1, 1, n, i) = a(1, 1, n, i)+coef
                a(2, 2, n, i) = a(2, 2, n, i)+coef
                a(3, 3, n, i) = a(3, 3, n, i)+coef
            end do
        end do
!
    end do

    if (ldiag) then
        diag = 0.d0
        somme = 0.d0
        do i = 1, 3
            do j = 1, nno
                somme(i) = somme(i)+a(i, i, j, j)
            end do
            alfam(i) = wgt/somme(i)
        end do
        do j = 1, nno
            do i = 1, 3
                diag(i, j) = a(i, i, j, j)*alfam(i)
            end do
        end do
        a = 0.d0
        do k = 1, 3
            do i = 1, nno
                a(k, k, i, i) = diag(k, i)
            end do
        end do
    end if
!
    do k = 1, nddl
        do l = 1, nddl
            do i = 1, nno
                kkd = ((nddl*(i-1)+k-1)*(nddl*(i-1)+k))/2
                do j = 1, i
                    kk = kkd+nddl*(j-1)+l
                    zr(jvMatr+kk-1) = a(k, l, i, j)
                end do
            end do
        end do
    end do
!
end subroutine
