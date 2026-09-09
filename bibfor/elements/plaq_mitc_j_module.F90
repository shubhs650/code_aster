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
!     RMQU9 - MODULE DE LA MATRICE J
! interface c++
! aslint: disable=W1403
module c_interface_plaq_mitc_j
    use iso_c_binding
    implicit none

    interface
        subroutine BP4_qu9_Fortran(w, nw, coordinate_dofs, ncd, &
            &entity_local_index, ne, cst, ncst, A) bind(C, name="BP4_qu9_Fortran")
            import :: c_double, c_int
            implicit none
!
            ! Longueurs des tableaux
            integer(c_int), value, intent(in) :: nw
            integer(c_int), value, intent(in) :: ncd
            integer(c_int), value, intent(in) :: ncst
            integer(c_int), value, intent(in) :: ne
!
            ! Tableaux d'entrée
            real(c_double), intent(in) :: w(*)
            real(c_double), intent(in) :: coordinate_dofs(*)
            real(c_double), intent(in) :: cst(*)
            integer(c_int), intent(in) :: entity_local_index(*)
!
            ! Tableau de sortie Matrice de rigidité (taille supposée connue)
            real(c_double), intent(out) :: A(*)
        end subroutine BP4_qu9_Fortran
        subroutine BP5_qu9_Fortran(w, nw, coordinate_dofs, ncd, &
            &entity_local_index, ne, cst, ncst, A) bind(C, name="BP5_qu9_Fortran")
            import :: c_double, c_int
            implicit none
!
            ! Longueurs des tableaux
            integer(c_int), value, intent(in) :: nw
            integer(c_int), value, intent(in) :: ncd
            integer(c_int), value, intent(in) :: ncst
            integer(c_int), value, intent(in) :: ne
!
            ! Tableaux d'entrée
            real(c_double), intent(in) :: w(*)
            real(c_double), intent(in) :: coordinate_dofs(*)
            real(c_double), intent(in) :: cst(*)
            integer(c_int), intent(in) :: entity_local_index(*)
!
            ! Tableau de sortie Matrice de rigidité (taille supposée connue)
            real(c_double), intent(out) :: A(*)
        end subroutine BP5_qu9_Fortran
    end interface

end module c_interface_plaq_mitc_j
