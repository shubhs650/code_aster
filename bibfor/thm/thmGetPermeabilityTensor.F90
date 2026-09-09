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
subroutine thmGetPermeabilityTensor(ds_thm, &
                                    ndim, phi, endo, &
                                    tperm)
!
    use Behaviour_type
    use MaterialPara_type
    use THM_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/rcvala.h"
#include "asterfort/tpermh.h"
#include "asterfort/utmess.h"
!
    type(THM_DS), intent(in) :: ds_thm
    integer(kind=8), intent(in) :: ndim
    real(kind=8), intent(in) :: phi, endo
    real(kind=8), intent(out) :: tperm(ndim, ndim)
!
! --------------------------------------------------------------------------------------------------
!
! THM
!
! Get permeability tensor
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_thm           : datastructure for THM
! In  ndim             : dimension of space (2 or 3)
! In  phi              : porosity
! In  endo             : damage
! Out tperm            : permeability tensor
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbProp = 4
    integer(kind=8) :: propCode(nbProp)
    real(kind=8) :: propVale(nbProp)
    character(len=16), parameter :: propName(nbProp) = (/'PERM_IN ', 'PERMIN_L', &
                                                         'PERMIN_N', 'PERMIN_T'/)
    integer(kind=8) :: aniso
!
! --------------------------------------------------------------------------------------------------
!
    tperm(1:ndim, 1:ndim) = 0.d0
    aniso = 0
    propVale = 0.d0

! - Read parameters (intrinsic permeability)
    propVale(1) = 1.d0
    if ((ds_thm%ds_behaviour%rela_hydr .eq. 'HYDR_UTIL') .or. &
        (ds_thm%ds_behaviour%rela_hydr .eq. 'HYDR_VGM') .or. &
        (ds_thm%ds_behaviour%rela_hydr .eq. 'HYDR_VGC') .or. &
        (ds_thm%ds_behaviour%rela_hydr .eq. 'HYDR_TABBAL')) then
        call rcvala(ds_thm%ds_behaviour%BEHInteg%materPara%jvMaterCode, &
                    ' ', 'THM_DIFFU', &
                    1, 'PORO', [phi], &
                    1, propName, propVale, &
                    propCode, 0, nan='NON')

        if (propCode(1) .eq. 1) then
! --------- Anisotropic
            call rcvala(ds_thm%ds_behaviour%BEHInteg%materPara%jvMaterCode, &
                        ' ', 'THM_DIFFU', &
                        1, 'PORO', [phi], &
                        1, propName(3), propVale(3), &
                        propCode, 0, nan='NON')
            if (propCode(1) .eq. 0) then
                aniso = 1
                call rcvala(ds_thm%ds_behaviour%BEHInteg%materPara%jvMaterCode, &
                            ' ', 'THM_DIFFU', &
                            1, 'PORO', [phi], &
                            1, propName(2), propVale(2), &
                            propCode, 0, nan='NON')
            else
                aniso = 2
                call rcvala(ds_thm%ds_behaviour%BEHInteg%materPara%jvMaterCode, &
                            ' ', 'THM_DIFFU', &
                            1, 'PORO', [phi], &
                            1, propName(2), propVale(2), &
                            propCode, 0, nan='NON')
                call rcvala(ds_thm%ds_behaviour%BEHInteg%materPara%jvMaterCode, &
                            ' ', 'THM_DIFFU', &
                            1, 'PORO', [phi], &
                            1, propName(4), propVale(4), &
                            propCode, 0, nan='NON')
            end if
        else if (propCode(1) .eq. 0) then
! --------- Isotropic
            aniso = 0
        end if
    else if (ds_thm%ds_behaviour%rela_hydr .eq. 'HYDR_ENDO') then
        if ((ds_thm%ds_behaviour%rela_meca .eq. 'MAZARS') .or. &
            (ds_thm%ds_behaviour%rela_meca .eq. 'ENDO_ISOT_BETON')) then
            aniso = 0
            call rcvala(ds_thm%ds_behaviour%BEHInteg%materPara%jvMaterCode, &
                        ' ', 'THM_DIFFU', &
                        1, 'ENDO', [endo], &
                        1, ['PERM_END'], propVale(1), &
                        propCode, 1)
            call rcvala(ds_thm%ds_behaviour%BEHInteg%materPara%jvMaterCode, &
                        ' ', 'THM_DIFFU', &
                        1, 'ENDO', [endo], &
                        1, ['PERM_END'], propVale(2), &
                        propCode, 1)
            call rcvala(ds_thm%ds_behaviour%BEHInteg%materPara%jvMaterCode, &
                        ' ', 'THM_DIFFU', &
                        1, 'ENDO', [endo], &
                        1, ['PERM_END'], propVale(3), &
                        propCode, 1)
        else
            call utmess('F', 'THM1_43', nk=2, &
                        valk=[ds_thm%ds_behaviour%rela_hydr, ds_thm%ds_behaviour%rela_meca])
        end if
    else
        ASSERT(ASTER_FALSE)
    end if

! - Compute permeability tensor
    call tpermh(ndim, ds_thm%ds_behaviour%BEHInteg%materPara%lcsPara%lcsAngle, aniso, propVale, &
                tperm)
!
end subroutine
