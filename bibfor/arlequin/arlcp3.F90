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
subroutine arlcp3(nbma1, nbma2, numno1, numno2, m3dea, &
                  m1dea, numn1t, numn2t, len1, len2, &
                  lisrel, charge)
!
!
! ----------------------------------------------------------------------
!
! ROUTINE ARLEQUIN
!
! ASSEMBLAGE DANS LES MATRICES ELEMENTAIRES DE COUPLAGE ARLEQUIN
!
! ----------------------------------------------------------------------
!
    implicit none
!
#include "asterfort/afrela.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "MeshTypes_type.h"
!
!     ARGUMENTS:
!     ----------
!
    integer(kind=8) :: nbnomx
    parameter(nbnomx=MT_NNOMAX)
    integer(kind=8) :: nbma1, nbma2
    integer(kind=8) :: len1, len2
    real(kind=8) :: m3dea(12, 3*nbnomx, nbma1), m1dea(12, 12, nbma2)
    character(len=19) :: lisrel
    character(len=8) :: charge
    character(len=5), dimension(nbnomx+2, nbma1) :: numno1
    character(len=5), dimension(2, nbma2) :: numno2
    character(len=5), dimension(nbnomx*nbma1) :: numn1t
    character(len=5), dimension(2*nbma2) :: numn2t
!
!-----------------------------------------------------------------------
    integer(kind=8) :: nbterm
    complex(kind=8) :: betac, coefc(3*len1+6*len2)
    character(len=4) :: typcoe, typval
    integer(kind=8) :: jdime(3*len1+6*len2, 3)
    real(kind=8) :: beta, direct(3*len1+6*len2, 3)
    real(kind=8) :: coefri(3*len1+6*len2)
    integer(kind=8) :: i, j, k, m, n, p, q
    real(kind=8) :: m1dass(6*len2, 6*len2), m3dass(6*len2, 3*len1)
    character(len=8) :: ddl1(3*len1+6*len2)
    real(kind=8) :: mmixe1(6*len2, 3*len1+6*len2)
    character(len=8) :: noeud1(3*len1+6*len2)
!-----------------------------------------------------------------------
    call jemarq()
!
! --- INITIALISATIONS
!
    m1dass = 0.0
    m3dass = 0.0
    mmixe1 = 0.0
!
! --- ASSEMBLAGE 1D
!
    do i = 1, len2
        do j = 1, len2
            do k = 1, nbma2
                do m = 1, 2
                    if (numno2(m, k) == numn2t(i)) then
                        do n = 1, 2
                            if (numno2(n, k) == numn2t(j)) then
                                do p = 1, 6
                                    do q = 1, 6
                                        m1dass(6*(i-1)+p, 6*(j-1)+q) = &
                                            m1dass(6*(i-1)+p, 6*(j-1)+q) &
                                            +m1dea(6*(m-1)+p, 6*(n-1)+q, k)
                                    end do
                                end do
                            end if
                        end do
                    end if
                end do
            end do
        end do
    end do
!
! --- ASSEMBLAGE 3D
!
    do i = 1, len2
        do j = 1, len1
            do k = 1, nbma1
                do m = 1, 2
                    if (numno1(m, k) == numn2t(i)) then
                        do n = 1, nbnomx
                            if (numno1(2+n, k) == numn1t(j)) then
                                do p = 1, 6
                                    do q = 1, 3
                                        m3dass(6*(i-1)+p, 3*(j-1)+q) = &
                                            m3dass(6*(i-1)+p, 3*(j-1)+q) &
                                            +m3dea(6*(m-1)+p, 3*(n-1)+q, k)
                                    end do
                                end do
                            end if
                        end do
                    end if
                end do
            end do
        end do
    end do
!
! --- CONCATENATION DES MATRICES DE COUPLAGE ASSEMBLEES
!
    do i = 1, 6*len2
        do j = 1, 6*len2
            mmixe1(i, j) = m1dass(i, j)
        end do
        do j = 1, 3*len1
            mmixe1(i, 6*len2+j) = -1*m3dass(i, j)
        end do
    end do
!
!
! --- AFFECTATION DES RELATIONS ARLEQUIN
!
    nbterm = 6*len2+3*len1
!
    do i = 1, nbterm
        coefc(i) = 0
        jdime(i, :) = [0, 0, 0]
        direct(i, :) = [0, 0, 0]
        ddl1(i) = '000'
    end do
!
    do i = 1, len2
        do j = 1, 6
            noeud1(6*(i-1)+j) = 'N'//numn2t(i)
        end do
        ddl1(6*i-5) = 'DX '
        ddl1(6*i-4) = 'DY '
        ddl1(6*i-3) = 'DZ '
        ddl1(6*i-2) = 'DRX'
        ddl1(6*i-1) = 'DRY'
        ddl1(6*i) = 'DRZ'
    end do
!
    do i = 1, len1
        do j = 1, 3
            noeud1(6*len2+3*(i-1)+j) = 'N'//numn1t(i)
        end do
        ddl1(6*len2+3*i-2) = 'DX '
        ddl1(6*len2+3*i-1) = 'DY '
        ddl1(6*len2+3*i) = 'DZ '
    end do
!
    beta = 0.0d0
    betac = (1.0d0, 0.0d0)
    typcoe = 'REEL'
    typval = 'REEL'
!
    do i = 1, 6*len2
        coefri = mmixe1(i, :)
        call afrela(coefri(1:nbterm), coefc(1:nbterm), ddl1(1:nbterm), noeud1(1:nbterm), &
                    jdime(1:nbterm, :), direct(1:nbterm, :), nbterm, beta, betac, &
                    ' ', typcoe, typval, 0.d0, lisrel)
    end do
!
    call jedema()
!
end subroutine
