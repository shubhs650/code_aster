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
subroutine get_elas_id(jvMaterCode, elasID, elasKeyword_)
!
    implicit none
!
#include "asterfort/ElasticityMaterial_type.h"
#include "asterfort/rccoma.h"
#include "asterfort/utmess.h"
!
    integer(kind=8), intent(in) :: jvMaterCode
    integer(kind=8), intent(out) :: elasID
    character(len=*), optional, intent(out) :: elasKeyword_
!
! --------------------------------------------------------------------------------------------------
!
! Comportment utility
!
! Get elasticity type
!
! --------------------------------------------------------------------------------------------------
!
! In  jvMaterCode      : coded material address
! Out elasID           : type of elasticity
! Out elasKeyword      : factor keyword for type of elasticity parameters
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: elasKeyword
!
! --------------------------------------------------------------------------------------------------
!

! - Keyword for elasticity parameters in material
    call rccoma(jvMaterCode, 'ELAS', 1, elasKeyword)

! - Type of elasticity (Isotropic/Orthotropic/Transverse isotropic)
    if (elasKeyword .eq. 'ELAS' .or. &
        elasKeyword .eq. 'ELAS_HYPER' .or. &
        elasKeyword .eq. 'ELAS_HYPER_VISC' .or. &
        elasKeyword .eq. 'ELAS_META') then
        elasID = ELAS_ISOT

    elseif (elasKeyword .eq. 'ELAS_MEMBRANE') then
        elasID = ELAS_MEMBRANE

    elseif (elasKeyword .eq. 'ELAS_GLRC') then
        elasID = ELAS_GLRC

    elseif (elasKeyword .eq. 'ELAS_DHRC') then
        elasID = ELAS_DHRC

    elseif (elasKeyword .eq. 'ELAS_ORTH') then
        elasID = ELAS_ORTH

    elseif (elasKeyword .eq. 'ELAS_ISTR') then
        elasID = ELAS_ISTR

    elseif (elasKeyword .eq. 'ELAS_VISCO') then
        elasID = ELAS_VISC_ISOT

    elseif (elasKeyword .eq. 'ELAS_VISCO_ORTH') then
        elasID = ELAS_VISC_ORTH

    elseif (elasKeyword .eq. 'ELAS_VISCO_ISTR') then
        elasID = ELAS_VISC_ISTR

    elseif (elasKeyword .eq. 'ELAS_COQUE') then
        elasID = ELAS_SHELL

    elseif (elasKeyword .eq. 'ELAS_COQMU') then
        elasID = ELAS_COMPOSITE

    else
        call utmess('F', 'COMPOR5_15', sk=elasKeyword)
    end if
!
    if (present(elasKeyword_)) then
        elasKeyword_ = elasKeyword
    end if

end subroutine
