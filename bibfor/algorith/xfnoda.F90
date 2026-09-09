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
subroutine xfnoda(ds_thm, &
                  mecani, press1, enrmec, dimenr, &
                  dimcon, ndim, congem, &
                  r, enrhyd, nfh)
!
    use MaterialPara_module
    use MaterialPara_type
    use THM_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/rcvalb.h"
!
    type(THM_DS), intent(in) :: ds_thm
    integer(kind=8) :: mecani(5), press1(7), enrmec(3), dimenr, dimcon
    integer(kind=8) :: ndim, yaenrm, adenme
    integer(kind=8) :: enrhyd(3), yaenrh, adenhy, nfh
    real(kind=8) :: congem(dimcon), r(dimenr)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: kpgFPG1 = 1, kspFPG1 = 1
    character(len=8), parameter :: famiFPG1 = "FPG1"
    type(Material_Para) :: materParaFPG1
    character(len=8), parameter :: poum = "+"
    real(kind=8), parameter :: rac2 = sqrt(2.d0)
    integer(kind=8), parameter :: nbProp = 3
    real(kind=8) :: propVale(nbProp)
    integer(kind=8) :: propCode(nbProp)
    character(len=8), parameter :: propName(nbProp) = (/'PESA_X', 'PESA_Y', 'PESA_Z'/)
    integer(kind=8) :: addeme, adcome
    integer(kind=8) :: addep1, adcp11, i, ifh
    real(kind=8) :: gravity(3)
    type(Material_Para) :: materPara
!
! --------------------------------------------------------------------------------------------------
!
    materPara = ds_thm%ds_behaviour%BEHInteg%materPara

! - Copy material parameters with other scheme parameters
    call copyMaterPara(materPara, famiFPG1, kpgFPG1, kspFPG1, &
                       materParaFPG1)

! - Get parameters of gravity
    call rcvalb(materParaFPG1%schemePara%fami, &
                materParaFPG1%schemePara%kpg, &
                materParaFPG1%schemePara%ksp, &
                poum, &
                materParaFPG1%jvMaterCode, &
                ' ', 'THM_DIFFU', &
                0, ' ', [0.d0], &
                nbProp, propName, propVale, &
                propCode, 1)
    gravity(1) = propVale(1)
    gravity(2) = propVale(2)
    gravity(3) = propVale(3)

! ======================================================================
! --- DETERMINATION DES VARIABLES CARACTERISANT LE MILIEU --------------
! ======================================================================
    addeme = mecani(2)
    addep1 = press1(3)
    adcp11 = press1(4)
    adcome = mecani(3)
    yaenrm = enrmec(1)
    adenme = enrmec(2)
    yaenrh = enrhyd(1)
    adenhy = enrhyd(2)
! ======================================================================
! --- COMME CONGEM CONTIENT LES VRAIES CONTRAINTES ET ------------------
! --- COMME PAR LA SUITE ON TRAVAILLE AVEC SQRT(2)*SXY -----------------
! --- ON COMMENCE PAR MODIFIER LES CONGEM EN CONSEQUENCE ---------------
! ======================================================================
    if (ds_thm%ds_elem%l_dof_meca) then
        do i = 4, 6
            congem(adcome+6+i-1) = congem(adcome+6+i-1)*rac2
            congem(adcome+i-1) = congem(adcome+i-1)*rac2
        end do
    end if
! ======================================================================
! --- CALCUL DU RESIDU R (TERMES CLASSIQUES ) --------------------------
! ======================================================================
    if (ds_thm%ds_elem%l_dof_meca) then
        do i = 1, 6
            r(addeme+ndim+i-1) = r(addeme+ndim+i-1)+congem(adcome-1+i)
        end do
        do i = 1, 6
            r(addeme+ndim-1+i) = r(addeme+ndim-1+i)+congem(adcome+6+i-1)
        end do
        if (ds_thm%ds_elem%l_dof_pre1) then
            do i = 1, ndim
                r(addeme+i-1) = r(addeme+i-1)-gravity(i)*congem(adcp11)
            end do
        end if
    end if
! ======================================================================
! --- CALCUL DU RESIDU R (TERMES HEAVISIDE ) ---------------------------
! ======================================================================
    if (yaenrm .eq. 1) then
        if (ds_thm%ds_elem%l_dof_meca) then
            if (ds_thm%ds_elem%l_dof_pre1) then
                do ifh = 1, nfh
                    do i = 1, ndim
                        r(adenme+i-1+(ifh-1)*(ndim+1)) = &
                            r(adenme+i-1+(ifh-1)*(ndim+1))-gravity(i)*congem(adcp11)
                    end do
                end do
            end if
        end if
    end if
!
end subroutine
