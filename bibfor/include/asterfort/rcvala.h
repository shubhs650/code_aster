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
!
          interface
              subroutine rcvala(jmat, nomat, phenom, nbpar, nompar, valpar,    &
       &nbres, nomres, valres, icodre, iarret, nan)
                  integer(kind=8), intent(in) :: nbres
                  integer(kind=8), intent(in) :: nbpar
                  integer(kind=8), intent(in) :: jmat
                  character(len=*), intent(in) :: nomat
                  character(len=*), intent(in) :: phenom
                  character(len=*), intent(in) :: nompar(nbpar)
                  real(kind=8), intent(in) :: valpar(nbpar)
                  character(len=*), intent(in) :: nomres(nbres)
                  real(kind=8), intent(out) :: valres(nbres)
                  integer(kind=8), intent(out) :: icodre(nbres)
                  integer(kind=8), intent(in) :: iarret
                  character(len=3), optional, intent(in) :: nan
              end subroutine rcvala
          end interface
