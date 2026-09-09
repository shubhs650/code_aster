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
! aslint: disable=W1504,W1306
!
subroutine coeime(ds_thm, &
                  lSigm, lVari, lMatr, &
                  option, &
                  ndim, dimdef, dimcon, &
                  addeme, addep1, &
                  nbvari, npg, npi, &
                  defgep, defgem, sigm, sigp, varim, &
                  varip, ouvh, tlint, drde, kpi, &
                  retcom)
!
    use MaterialPara_module
    use MaterialPara_type
    use THM_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/lcejex.h"
#include "asterfort/lcejli.h"
#include "asterfort/lcjohm.h"
#include "asterfort/rcvalb.h"
!
    type(THM_DS), intent(in) :: ds_thm
    character(len=16), intent(in) :: option
    aster_logical, intent(in) :: lSigm, lVari, lMatr
    integer(kind=8), intent(in) :: ndim, dimcon, dimdef
    integer(kind=8), intent(in) :: addeme, addep1, npg, kpi, npi, nbvari
    real(kind=8), intent(in) :: defgem(dimdef), defgep(dimdef)
    real(kind=8), intent(in) :: sigm(dimcon)
    real(kind=8), intent(inout) :: sigp(dimcon)
    real(kind=8), intent(in) :: varim(nbvari)
    real(kind=8), intent(inout) :: varip(nbvari)
    real(kind=8), intent(out) :: ouvh, tlint
    real(kind=8), intent(inout) :: drde(dimdef, dimdef)
    integer(kind=8), intent(out) :: retcom
!
! --------------------------------------------------------------------------------------------------
!
!  INTEGRATION DE LA LOI DE COMPORTEMENT MECANIQUE ET RENVOI DE LA LOI CUBIQUE
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_thm           : datastructure for THM
! IN MECA   : COMPORTEMENT MECA
! IN IMATE  : CODE MATERIAU
! IN RESI   : FULL_MECA OU RAPH_MECA
! IN RIGI   : FULL_MECA OU RIGI_MECA
! IN NDIM   : DIMENSION ESPACE
! IN DIMDEF : DIMENSION DEFORMATION GENERALISEE
! IN DIMCON : DIMENSION VECTEUR CONTRAINTES GENERALISEES
! IN ADDEME : ADRESSE DES DEFORMATIONS MECANIQUES
! IN ADDEP1 : ADRESSE DES DEFORMATIONS PRESSION 1
! IN NBVARI : NOMBRE DE VARIABLES INTERNES
! IN ADVIME : ADRESSE DES VI MECANIQUES
! IN ADVICO : ADRESSE DES VI DE COUPLAGE
! IN NPG    : NOMBRE DE POINTS DE GAUSS
! IN DEFGEP : DEFORMATIONS AU TEMPS PLUS
! IN DERGEM : DEFORMATIONS AU TEMPS MOINS
! IN SIGM   : CONTRAINTES AU TEMPS MOINS
! IN VARIM  : VARIABLES INTERNES AU TEMPS MOINS
! IN KPI    : POINT D'INTEGRATION
! =====================================================================
! OUT SIGP  : CONTRAINTES AU TEMPS PLUS
! OUT VARIP : VARIABLES INTERNES AU TEMPS PLUS
! OUT OUVH  : OUVERTURE NORMALE DU JOINT
! OUT TLINT : PERMEABILITE LONGITUDINALE
! OUT DRDE  : MATRICE DE RIGIDITE
! OUT RETCOM : RETOUR LOI DE COMPORTEMENT
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: kpgFPG1 = 1, kspFPG1 = 1
    character(len=8), parameter :: famiFPG1 = "FPG1"
    type(Material_Para) :: materParaFPG1
    character(len=8), parameter :: poum = "+"
    integer(kind=8), parameter :: nbProp = 2
    character(len=16), parameter :: propName(nbProp) = (/'OUV_FICT', 'UN_SUR_N'/)
    real(kind=8) :: propVale(nbProp)
    integer(kind=8) :: propCode(nbProp)
    integer(kind=8) :: i, j
    real(kind=8) :: da(ndim), dsidep(6, 6), ouvfic, unsurn
    character(len=16) :: relaMeca
    integer(kind=8) :: advime, advico, vicphi
    type(Material_Para) :: materPara
!
! --------------------------------------------------------------------------------------------------
!
    ouvh = 0.d0
    tlint = 0.d0

! - Copy material parameters with other scheme parameters
    materPara = ds_thm%ds_behaviour%BEHInteg%materPara
    call copyMaterPara(materPara, famiFPG1, kpgFPG1, kspFPG1, &
                       materParaFPG1)

    relaMeca = ds_thm%ds_behaviour%rela_meca
    advime = ds_thm%ds_behaviour%advime
    advico = ds_thm%ds_behaviour%advico
    vicphi = ds_thm%ds_behaviour%vicphi

    if (relaMeca .eq. 'JOINT_BANDIS') then
        call lcjohm(materParaFPG1, &
                    lSigm, lMatr, lVari, &
                    kpi, npg, &
                    addeme, advico, ndim, dimdef, &
                    dimcon, nbvari, defgem, defgep, varim, &
                    varip, sigm, sigp, drde, ouvh, &
                    retcom)
        tlint = ouvh**2/12.d0
        if (lVari) then
            varip(advime) = tlint
        end if
        if (lSigm) then
            if (ds_thm%ds_elem%l_dof_pre1) then
                sigp(1+ndim) = -defgep(addep1)
            end if
        end if
        if ((lMatr) .and. (kpi .le. npg)) then
            if (ds_thm%ds_elem%l_dof_pre1) then
                drde(addeme, addep1) = -1.d0
            end if
        end if
    end if

    if (relaMeca .eq. 'CZM_LIN_REG') then
        do i = 1, ndim
            da(i) = defgep(i)-defgem(i)
        end do
        call lcejli(materPara%schemePara%fami, &
                    materPara%schemePara%kpg, &
                    materPara%schemePara%ksp, &
                    ndim, &
                    materPara%jvMaterCode, &
                    option, defgem, da, sigp, dsidep, &
                    varim(advime), varip(advime))

        if (nint(varip(advime)) .eq. 2) then
            unsurn = 0.d0
        else
            call rcvalb(materParaFPG1%schemePara%fami, &
                        materParaFPG1%schemePara%kpg, &
                        materParaFPG1%schemePara%ksp, &
                        poum, &
                        materParaFPG1%jvMaterCode, &
                        ' ', 'THM_RUPT', &
                        0, ' ', [0.d0], &
                        nbProp, propName, propVale, &
                        propCode, 1)
            ouvfic = propVale(1)
            unsurn = propVale(2)
        end if
        if (lMatr) then
            if (kpi .le. npg) then
                do i = 1, ndim
                    do j = 1, ndim
                        drde(i, j) = dsidep(i, j)
                    end do
                end do
                if ((ds_thm%ds_elem%l_dof_pre1) .and. &
                    ((nint(varip(advime+2)) .eq. 1) .or. (nint(varip(advime+2)) .eq. 2))) then
                    drde(1, addep1) = -1.d0
                end if
            end if
            if ((kpi .gt. npg) .or. (npi .eq. npg)) then
                drde(addep1, addep1) = drde(addep1, addep1)-unsurn
            end if

            ouvh = varim(advico+vicphi)
            if (nint(varim(3)) .eq. 0) then
                ouvh = ouvfic
            end if
            tlint = ouvh**2/12
        end if
        if (lSigm) then
            if ((ds_thm%ds_elem%l_dof_pre1) .and. &
                ((nint(varip(advime+2)) .eq. 1) .or. (nint(varip(advime+2)) .eq. 2))) then
                sigp(1+ndim) = -defgep(addep1)
            end if
        end if
        if (lVari) then
            varip(advico+vicphi) = defgep(1)
            ouvh = varip(advico+vicphi)
            if ((nint(varip(3)) .eq. 0)) then
                ouvh = ouvfic
            end if
            tlint = ouvh**2/12
            varip(advico+vicphi) = defgep(1)+defgep(addep1)*unsurn
        end if
    end if

    if (relaMeca .eq. 'CZM_EXP_REG') then
        do i = 1, ndim
            da(i) = defgep(i)-defgem(i)
        end do
        call lcejex(materPara%schemePara%fami, &
                    materPara%schemePara%kpg, &
                    materPara%schemePara%ksp, &
                    ndim, &
                    materPara%jvMaterCode, &
                    option, defgem, da, sigp, dsidep, &
                    varim(advime), varip(advime))

        call rcvalb(materParaFPG1%schemePara%fami, &
                    materParaFPG1%schemePara%kpg, &
                    materParaFPG1%schemePara%ksp, &
                    poum, &
                    materParaFPG1%jvMaterCode, &
                    ' ', 'THM_RUPT', &
                    0, ' ', [0.d0], &
                    nbProp, propName, propVale, &
                    propCode, 1)
        ouvfic = propVale(1)
        unsurn = propVale(2)
        if (lMatr) then
            if (kpi .le. npg) then
                do i = 1, ndim
                    do j = 1, ndim
                        drde(i, j) = dsidep(i, j)
                    end do
                end do
                if ((ds_thm%ds_elem%l_dof_pre1) .and. (nint(varip(advime+2)) .eq. 1)) then
                    drde(1, addep1) = -1.d0
                end if
            end if
            if ((kpi .gt. npg) .or. (npi .eq. npg)) then
                drde(addep1, addep1) = drde(addep1, addep1)-unsurn
            end if
            ouvh = varim(advico+vicphi)
            if (nint(varim(3)) .eq. 0) then
                ouvh = ouvfic
            end if
            tlint = ouvh**2/12
        end if
        if (lSigm) then
            if ((ds_thm%ds_elem%l_dof_pre1) .and. (nint(varip(advime+2)) .eq. 1)) then
                sigp(1+ndim) = -defgep(addep1)
            end if
        end if
        if (lVari) then
            varip(advico+vicphi) = defgep(1)
            ouvh = varip(advico+vicphi)
            if (nint(varip(3)) .eq. 0) then
                ouvh = ouvfic
            end if
            tlint = ouvh**2/12
            varip(advico+vicphi) = defgep(1)+defgep(addep1)*unsurn
        end if
    end if
!
end subroutine
