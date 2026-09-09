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
subroutine te0529(option, nomte)
!
    use BehaviourStrain_module
    use BehaviourStrain_type
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/epstmc.h"
#include "asterfort/jevech.h"
#include "asterfort/tecach.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
!     BUT: CALCUL DES DEFORMATIONS LIEES AUX VARIABLES DE COMMANDE
!          AUX POINTS D'INTEGRATION DES ELEMENTS ISOPARAMETRIQUES 3D
!
!          OPTIONS : 'EPVC_ELGA'
!    CINQ COMPOSANTES :
!    EPTHER_L = DILATATION THERMIQUE (LONGI)   : ALPHA_L*(T-TREF)
!    EPTHER_T = DILATATION THERMIQUE (TRANSV)   : ALPHA_T*(T-TREF)
!    EPTHER_N = DILATATION THERMIQUE (NORMLALE)   : ALPHA_N*(T-TREF)
!    EPSECH = RETRAIT DE DESSICCATION : -K_DESSIC(SREF-SECH)
!    EPHYDR = RETRAIT ENDOGENE        : -B_ENDOGE*HYDR
!    EPPTOT = RETRAIT DU A LA PRESSION DE FLUIDE EN THM CHAINEE :
!             -BIOT*PTOT
!
!     ENTREES  ---> OPTION : OPTION DE CALCUL
!              ---> NOMTE  : NOM DU TYPE ELEMENT
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: ksp = 1, nbEpsi = 6
    character(len=8), parameter :: fami = 'RIGI'
    integer(kind=8) :: ndim, nno, npg, kpg, iEpsi, iret
    integer(kind=8) :: jvGeom, jvTime, jvEpsi, jvMaterc
    real(kind=8) :: time
    real(kind=8) :: epsiVarc(162)
    real(kind=8) :: epsiSech(nbEpsi), epsiTher(nbEpsi), epsiHydr(nbEpsi), epsiPtot(nbEpsi)
    type(All_Varc_Strain) :: allVarcStrain
    type(Material_Para) :: materPara
!
! --------------------------------------------------------------------------------------------------
!
    epsiVarc = 0.d0

! - Get element parameters
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, npg=npg)
    ASSERT(npg .le. 27)

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Get material parameters
    call jevech('PMATERC', 'L', jvMaterc)

! - Initializations of material parameters on current cell
    call initParaCell(fami, zi(jvMaterc), materPara)

! - Set local coordinate system from user
    call getUserLCS(ndim, nno, jvGeom, materPara%lcsPara)

! - Get current time
    call tecach('NNO', 'PINSTR', 'L', iret, iad=jvTime)
    if (jvTime .ne. 0) then
        time = zr(jvTime)
    else
        time = r8vide()
    end if

! - Compute
    do kpg = 1, npg
! ----- Initializations of material parameters on current integration point
        call initParaPoin(kpg, ksp, materPara)

! ----- Compute inelastic strains
        call epstmc(materPara, "+", time, ndim, &
                    VARC_STRAIN_ALL, allVarcStrain)

! ----- Get TEMP strains
        call getVarcStrain('+', VARC_STRAIN_TEMP, allVarcStrain, epsiTher)

! ----- Get SECH strains
        call getVarcStrain('+', VARC_STRAIN_SECH, allVarcStrain, epsiSech)

! ----- Get HYDR strains
        call getVarcStrain('+', VARC_STRAIN_HYDR, allVarcStrain, epsiHydr)

! ----- Get PTOT strains
        call getVarcStrain('+', VARC_STRAIN_PTOT, allVarcStrain, epsiPtot)

! ----- Get EPSA strains
        ! call getVarcStrain('+', VARC_STRAIN_EPSA, allVarcStrain, epsiEpsa)
        epsiVarc(1+nbEpsi*(kpg-1)) = epsiTher(1)
        epsiVarc(2+nbEpsi*(kpg-1)) = epsiTher(2)
        epsiVarc(3+nbEpsi*(kpg-1)) = epsiTher(3)
        epsiVarc(4+nbEpsi*(kpg-1)) = epsiSech(1)
        epsiVarc(5+nbEpsi*(kpg-1)) = epsiHydr(1)
        epsiVarc(6+nbEpsi*(kpg-1)) = epsiPtot(1)

    end do

! - Set output
    call jevech('PDEFOPG', 'E', jvEpsi)
    do kpg = 1, npg
        do iEpsi = 1, nbEpsi
            zr(jvEpsi+nbEpsi*(kpg-1)+iEpsi-1) = epsiVarc(nbEpsi*(kpg-1)+iEpsi)
        end do
    end do
!
end subroutine
