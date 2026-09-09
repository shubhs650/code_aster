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

module tenseur_dime_module

! =====================================================================
!  Utilitaires pour l'usage des tenseurs
! =====================================================================

    implicit none
    private
    public:: rs, kron, voigt, proten, identity, sph_norm, deviator, prod_vect, matrix_to_voigt
    public:: NDIM_TO_NDIMSI

#include "asterfort/assert.h"

    real(kind=8), parameter, dimension(6)::KRONECKER = [1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0]
    real(kind=8), parameter, dimension(6)::RACINE_2 = [1.d0, 1.d0, 1.d0, &
                                                       sqrt(2.d0), sqrt(2.d0), sqrt(2.d0)]
    real(kind=8), parameter             ::RAC2 = sqrt(2.d0)
    real(kind=8), parameter             ::RAC3 = sqrt(3.d0)
    integer(kind=8), parameter, dimension(3):: NDIM_TO_NDIMSI = [1, 4, 6]

contains

! =====================================================================
!  Extension d'un vecteur a la taille n complete par des zeros
! =====================================================================

    function rs(nout, vin) result(vout)
        implicit none
        integer(kind=8), intent(in)                  :: nout
        real(kind=8), dimension(:), intent(in) :: vin
        real(kind=8), dimension(nout)         :: vout
! ---------------------------------------------------------------------
        integer(kind=8) nin
! ---------------------------------------------------------------------
        nin = size(vin)
        if (nin .le. nout) then
            vout = 0
            vout(1:nin) = vin
        else
            vout = vin(1:nout)
        end if

    end function rs

! =====================================================================
!  Kronecker en representation vectorielle (1:ndimsi)
! =====================================================================

    function kron(ndimsi) result(kr)

        implicit none
        integer(kind=8), intent(in)            ::ndimsi
        real(kind=8), dimension(ndimsi):: kr
! ---------------------------------------------------------------------
        ASSERT(ndimsi .eq. 1 .or. ndimsi .eq. 4 .or. ndimsi .eq. 6)
        kr = KRONECKER(1:ndimsi)

    end function kron

! =====================================================================
!  Vecteur de transformation pour representation de Voigt (*rac2 sur cis)
! =====================================================================

    function voigt(ndimsi) result(rac2)

        implicit none
        integer(kind=8), intent(in)            ::ndimsi
        real(kind=8), dimension(ndimsi):: rac2
! ---------------------------------------------------------------------
        ASSERT(ndimsi .eq. 1 .or. ndimsi .eq. 4 .or. ndimsi .eq. 6)
        rac2 = RACINE_2(1:ndimsi)

    end function voigt

! =====================================================================
!  Passe la matrice symétrique mat en vecteur avec notation de voigt
! =====================================================================

    function matrix_to_voigt(mat, ndimsi) result(vect)

        implicit none

        real(kind=8), intent(in):: mat(:, :)
        integer(kind=8), intent(in)::ndimsi
        real(kind=8):: vect(ndimsi)
! ---------------------------------------------------------------------
        ASSERT(ndimsi .eq. 1 .or. ndimsi .eq. 4 .or. ndimsi .eq. 6)
        vect(1) = mat(1, 1)
        if (ndimsi .eq. 4 .or. ndimsi .eq. 6) then
            vect(2) = mat(2, 2)
            vect(3) = mat(3, 3)
            vect(4) = mat(1, 2)
        end if
        if (ndimsi .eq. 6) then
            vect(5) = mat(1, 3)
            vect(6) = mat(2, 3)
        end if
        vect = vect*voigt(ndimsi)

    end function matrix_to_voigt

! =====================================================================
!  matrice identite de taille n
! =====================================================================

    function identity(n) result(idm)
        implicit none
        integer(kind=8), intent(in) :: n
        real(kind=8), dimension(n, n) :: idm
! ---------------------------------------------------------------------
        integer(kind=8) :: i
! ---------------------------------------------------------------------
        idm = 0.d0
        do i = 1, n
            idm(i, i) = 1.d0
        end do

    end function identity

! =====================================================================
!  Deviateur d'un tenseur en representation de Voigt
! =====================================================================

    function deviator(u) result(w)
        implicit none
        real(kind=8), dimension(:), intent(in) :: u
        real(kind=8), dimension(size(u))      :: w
! ---------------------------------------------------------------------
        ASSERT(size(u) .ge. 3)
        w = u-sum(u(1:3))*kron(size(u))/3.d0

    end function deviator

! =====================================================================
!  Partie spherique normee d'un tenseur en representation de Voigt
! =====================================================================

    function sph_norm(u) result(y)
        implicit none
        real(kind=8), dimension(:), intent(in) :: u
        real(kind=8)                         :: y
! ---------------------------------------------------------------------
        ASSERT(size(u) .ge. 3)
        y = sum(u(1:3))/RAC3

    end function sph_norm

! =====================================================================
!  Produit tensoriel de deux vecteurs de dimension quelconque
! =====================================================================

    function proten(u, v) result(w)
        implicit none
        real(kind=8), dimension(:), intent(in) :: u, v
        real(kind=8), dimension(size(u), size(v)) :: w
! ---------------------------------------------------------------------
        integer(kind=8) :: i, j
! ---------------------------------------------------------------------
        do i = 1, size(u)
            do j = 1, size(v)
                w(i, j) = u(i)*v(j)
            end do
        end do

    end function proten

! =====================================================================
!  Produit vectoriel de deux vecteurs de dimension 3
! =====================================================================

    function prod_vect(u, v) result(w)
        implicit none
        real(kind=8), dimension(:), intent(in) :: u, v
        real(kind=8), dimension(size(u)) :: w
! ---------------------------------------------------------------------
        ASSERT(size(u) .eq. 3)
        ASSERT(size(v) .eq. 3)
        w(1) = u(2)*v(3)-u(3)*v(2)
        w(2) = u(3)*v(1)-u(1)*v(3)
        w(3) = u(1)*v(2)-u(2)*v(1)

    end function prod_vect

end module tenseur_dime_module
