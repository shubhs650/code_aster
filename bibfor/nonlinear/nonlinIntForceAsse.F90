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
subroutine nonlinIntForceAsse(typeAsse, list_func_acti, sdnume, &
                              ds_material, ds_constitutive, ds_system)
!
    use NonLin_Datastructure_type
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/NonLinear_type.h"
#include "asterfort/vtaxpy.h"
#include "asterfort/vtzero.h"
#include "asterfort/assert.h"
#include "asterfort/assvec.h"
#include "asterfort/exisd.h"
#include "asterfort/isfonc.h"
#include "asterfort/infdbg.h"
#include "asterfort/utmess.h"
#include "asterfort/nmdebg.h"
#include "asterfort/copisd.h"
!
    integer(kind=8), intent(in) :: typeAsse, list_func_acti(*)
    character(len=19), intent(in) :: sdnume
    type(NL_DS_Material), intent(in) :: ds_material
    type(NL_DS_Constitutive), intent(in) :: ds_constitutive
    type(NL_DS_System), intent(in) :: ds_system
!
! --------------------------------------------------------------------------------------------------
!
! MECA_NON_LINE - Algorithm
!
! Assembly for internal forces vector
!
! --------------------------------------------------------------------------------------------------
!
! In  typeAsse         : type of assembly
! In  list_func_acti   : list of active functionnalities
! In  sdnume           : datastructure for dof positions
! In  ds_material      : datastructure for material parameters
! In  ds_constitutive  : datastructure for constitutive laws management
! In  ds_system        : datastructure for non-linear system management
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ifm, niv
    integer(kind=8) :: iRet
    real(kind=8) :: vectCoef(2)
    character(len=19) :: vectElem(2)
    character(len=8) :: vevcprCurr, vevcprPrev
    aster_logical :: l_resi_comp
    character(len=19) :: cnvcpr
!
! --------------------------------------------------------------------------------------------------
!
    call infdbg('MECANONLINE', ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'MECANONLINE13_39')
    end if

! - Active functionnalities
    l_resi_comp = isfonc(list_func_acti, 'RESI_COMP')

! - Parameters for external state variables vector
    cnvcpr = ds_material%fvarc_pred(1:19)
    vevcprPrev = ds_material%vevcprPrev
    vevcprCurr = ds_material%vevcprCurr
    vectCoef(1) = +1.d0
    vectCoef(2) = -1.d0
    vectElem(1) = vevcprCurr
    vectElem(2) = vevcprPrev

! - Assemble
    if (typeAsse .eq. INTE_FORCE_COMB) then
! ----- Assemblage special: In: cninte + cnfnod + ds_material%fvarc_pred
        call vtzero(ds_system%cnfint)
        call exisd('CHAM_ELEM', ds_constitutive%code_pred, iRet)
        if (iRet .eq. 0) then
! --------- Some options (RIGI_MECA_ELAS) has no CODE_PRED
            call assvec('V', ds_system%cnfnod, 1, ds_system%vefnod, [1.d0], ds_system%nume_dof)
            call vtaxpy(1.d0, ds_system%cnfnod, ds_system%cnfint)
            call assvec('V', cnvcpr, 2, vectElem, vectCoef, ds_system%nume_dof)
            call vtaxpy(-1.d0, cnvcpr, ds_system%cnfint)
        else
            call assvec('V', cnvcpr, 2, vectElem, vectCoef, ds_system%nume_dof, &
                        maskElem_=ds_constitutive%code_pred, &
                        maskInve_=ASTER_TRUE)
            call assvec('V', ds_system%cnfnod, 1, ds_system%vefnod, [1.d0], ds_system%nume_dof, &
                        maskElem_=ds_constitutive%code_pred, &
                        maskInve_=ASTER_TRUE)
            call assvec('V', ds_system%cninte, 1, ds_system%veinte, [1.d0], ds_system%nume_dof, &
                        maskElem_=ds_constitutive%code_pred, &
                        maskInve_=ASTER_FALSE)
            call vtaxpy(1.d0, ds_system%cnfnod, ds_system%cnfint)
            call vtaxpy(1.d0, ds_system%cninte, ds_system%cnfint)
            call vtaxpy(-1.d0, cnvcpr, ds_system%cnfint)
        end if

    elseif (typeAsse .eq. INTE_FORCE_INTE) then
!       Assemblage que de l'intégration (RAPH_MECA / FULL_MECA / RIGI_MECA_TANG)
        call assvec('V', ds_system%cnfint, 1, ds_system%veinte, [1.d0], ds_system%nume_dof)
    elseif (typeAsse .eq. INTE_FORCE_NONE) then
!       Nothing to do
        ASSERT(.not. l_resi_comp)
    else
        ASSERT(ASTER_FALSE)
    end if

! - If RESI_COMP_RELA
    if (l_resi_comp .and. typeAsse .ne. INTE_FORCE_INTE) then
        call copisd('CHAMP_GD', 'V', ds_system%cnfint, ds_system%cncomp)
    end if

! - Debug
    if (niv .ge. 2) then
        call nmdebg('VECT', ds_system%cnfint, 6)
    end if
!
end subroutine
