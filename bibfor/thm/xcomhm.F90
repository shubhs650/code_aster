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
subroutine xcomhm(ds_thm, &
                  option, time_curr, &
                  ndim, dimdef, dimcon, nbvari, &
                  addeme, adcome, addep1, adcp11, &
                  addep2, addete, defgem, &
                  defgep, congem, congep, vintm, &
                  vintp, dsde, gravity, retcom, &
                  kpg, npg, dimenr, &
                  yaenrh, adenhy, nfh)
!
    use MaterialPara_module
    use THM_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/calcva.h"
#include "asterfort/tebiot.h"
#include "asterfort/thmEvalGravity.h"
#include "asterfort/thmGetParaBiot.h"
#include "asterfort/thmGetParaElas.h"
#include "asterfort/thmGetParaHydr.h"
#include "asterfort/thmGetParaTher.h"
#include "asterfort/thmGetPermeabilityTensor.h"
#include "asterfort/thmMatrHooke.h"
#include "asterfort/xcalfh.h"
#include "asterfort/xcalme.h"
#include "asterfort/xhmsat.h"
!
    type(THM_DS), intent(inout) :: ds_thm
    integer(kind=8) :: retcom, kpg, npg, nfh
    integer(kind=8) :: ndim, dimdef, dimcon, nbvari
    integer(kind=8) :: addeme, addep1, addep2, addete
    integer(kind=8) :: adcome, adcp11
    real(kind=8) :: defgem(1:dimdef), defgep(1:dimdef), congep(1:dimcon)
    real(kind=8) :: congem(1:dimcon), vintm(1:nbvari), vintp(1:nbvari)
    real(kind=8) :: time_curr
    character(len=16) :: option
    integer(kind=8) :: dimenr
    integer(kind=8) :: yaenrh, adenhy
    real(kind=8) :: dsde(1:dimcon, 1:dimenr)
    real(kind=8) :: gravity(3)
!
! --------------------------------------------------------------------------------------------------
!
! CALCULE LES CONTRAINTES GENERALISEES ET LA MATRICE TANGENTE AU POINT
! DE GAUSS SUIVANT LES OPTIONS DEFINIES
!
! --------------------------------------------------------------------------------------------------
!
! IN OPTION : OPTION DE CALCUL
! IN COMPOR : COMPORTEMENT
! IN IMATE  : MATERIAU CODE
! IN NDIM   : DIMENSION DE L'ESPACE
! IN DIMDEF : DIMENSION DU TABLEAU DES DEFORMATIONS GENERALISEES
!             AU POINT DE GAUSS CONSIDERE
! IN DIMCON : DIMENSION DU TABLEAU DES CONTRAINTES GENERALISEES
!             AU POINT DE GAUSS CONSIDERE
! IN NBVARI : NOMBRE TOTAL DE VARIABLES INTERNES AU POINT DE GAUSS
! IN ADDEME : ADRESSE DES DEFORMATIONS MECANIQUES
! IN ADDEP1 : ADRESSE DES DEFORMATIONS CORRESPONDANT A LA PRESSION 1
! IN ADDEP2 : ADRESSE DES DEFORMATIONS CORRESPONDANT A LA PRESSION 2
! IN ADDETE : ADRESSE DES DEFORMATIONS THERMIQUES
! IN ADCOME : ADRESSE DES CONTRAINTES MECANIQUES
! IN ADCP11 : ADRESSE DES CONTRAINTES FLUIDE 1 PHASE 1
! IN ADCP11 : ADRESSE DES CONTRAINTES FLUIDE 1 PHASE 2
! IN ADCP11 : ADRESSE DES CONTRAINTES FLUIDE 2 PHASE 1
! IN ADCP11 : ADRESSE DES CONTRAINTES FLUIDE 2 PHASE 2
! IN ADCOTE : ADRESSE DES CONTRAINTES THERMIQUES
! IN DEFGEM : DEFORMATIONS GENERALISEES A L'INSTANT MOINS
! IN DEFGEP : DEFORMATIONS GENERALISEES A L'INSTANT PLUS
! IN CONGEM : CONTRAINTES GENERALISEES A L'INSTANT MOINS
! IN VINTM  : VARIABLES INTERNES A L'INSTANT MOINS
! IN TYPMOD : MODELISATION (D_PLAN, AXI, 3D ?)
!
! OUT CONGEP : CONTRAINTES GENERALISEES A L'INSTANT PLUS
! OUT VINTP  : VARIABLES INTERNES A L'INSTANT PLUS
! OUT DSDE   : MATRICE TANGENTE CONTRAINTES DEFORMATIONS
!
! OUT RETCOM : RETOUR LOI DE COMPORTEMENT
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: p1, dp1, gradP1(3), p2, dp2, gradP2(3), temp, dtemp, gradTemp(3)
    real(kind=8) :: phi, rho11, epsv, deps(6), depsv
    real(kind=8) :: satur, endo
    real(kind=8) :: tbiot(6)
    real(kind=8) :: tperm(ndim, ndim)
!
! --------------------------------------------------------------------------------------------------
!
    retcom = 0

! - Update unknowns
    call calcva(ds_thm, ndim, &
                defgem, defgep, &
                addeme, addep1, addep2, addete, &
                depsv, epsv, deps, &
                temp, dtemp, gradTemp, &
                p1, dp1, gradP1, &
                p2, dp2, gradP2, &
                retcom)
    if (retcom .ne. 0) then
        goto 99
    end if

! - Get hydraulic parameters
    call thmGetParaHydr(ds_thm)

! - Get Biot parameters (for porosity evolution), paraThetaCpl
    call thmGetParaBiot(ds_thm)

! - Compute Biot tensor
    call tebiot(ds_thm, tbiot)

! - Get elastic parameters
    call thmGetParaElas(temp, ndim, ds_thm)
    call thmMatrHooke(ds_thm)

! - Get thermic parameters
    call thmGetParaTher(temp, ds_thm)

! - Compute generalized stresses and matrix for coupled quantities
    call xhmsat(ds_thm, option, &
                ndim, dimenr, &
                dimcon, nbvari, addeme, &
                adcome, &
                addep1, adcp11, congem, congep, vintm, &
                vintp, dsde, epsv, depsv, &
                dp1, phi, rho11, &
                satur, retcom, tbiot, &
                yaenrh, adenhy, nfh)
    if (retcom .ne. 0) then
        goto 99
    end if

! - Main select subroutine to integrate mechanical behaviour
    if (ds_thm%ds_elem%l_dof_meca .and. kpg .le. npg) then
        call xcalme(ds_thm, &
                    option, ndim, dimenr, &
                    dimcon, addeme, adcome, congep, &
                    dsde, deps)
        if (retcom .ne. 0) then
            goto 99
        end if
    end if

! - Get permeability tensor
    if ((option(1:9) .eq. 'FULL_MECA') .or. (option(1:9) .eq. 'RAPH_MECA')) then
        endo = vintp(1)
    else
        endo = vintm(1)
    end if
    call thmGetPermeabilityTensor(ds_thm, &
                                  ndim, phi, endo, &
                                  tperm)

! - Compute gravity
    call thmEvalGravity(ds_thm, time_curr, gravity)

! - Compute flux and stress for hydraulic
    if ((ds_thm%ds_elem%l_dof_pre1) .and. (yaenrh .eq. 1)) then
        call xcalfh(ds_thm, &
                    option, ndim, dimcon, &
                    addep1, adcp11, addeme, congep, dsde, &
                    gradP1, rho11, gravity, tperm, &
                    dimenr, &
                    adenhy, nfh)
        if (retcom .ne. 0) then
            goto 99
        end if
    end if
!
99  continue
!
end subroutine
