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
    subroutine nmextr_comp(field, field_disc, field_type, &
                           meshZ, modelZ, caraElemZ, &
                           ds_material, ds_constitutive, &
                           dispCurrZ, strxCurrZ, varcCurrZ, time, ligrelZ_)
        use NonLin_Datastructure_type
        character(len=19), intent(in) :: field
        character(len=24), intent(in) :: field_type
        character(len=4), intent(in) :: field_disc
        character(len=*), intent(in) :: modelZ, meshZ, caraElemZ
        type(NL_DS_Material), intent(in) :: ds_material
        type(NL_DS_Constitutive), intent(in) :: ds_constitutive
        character(len=*), intent(in) :: dispCurrZ, strxCurrZ, varcCurrZ
        real(kind=8), intent(in) :: time
        character(len=*), optional, intent(in) :: ligrelZ_
    end subroutine nmextr_comp
end interface
