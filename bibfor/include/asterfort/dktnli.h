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
interface
    subroutine dktnli(plateCara, plateOrie, &
                      BEHInteg, option, typmod, &
                      instm, instp, &
                      xyzl, uml, dul, &
                      btsig, ktan, codret)
        use Behaviour_type
        use plate_type
        type(plateCara_Para), intent(in) :: plateCara
        type(plateOrie_Para), intent(in) :: plateOrie
        type(Behaviour_Integ), intent(inout) :: BEHInteg
        character(len=16), intent(in) :: option
        character(len=8), intent(in) :: typmod(2)
        real(kind=8), intent(in) :: instm, instp
        real(kind=8), intent(in) :: xyzl(3, 4), uml(6, 4), dul(6, 4)
        real(kind=8), intent(out) :: ktan(576), btsig(6, 4)
        integer(kind=8), intent(out) :: codret
    end subroutine dktnli
end interface
