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
subroutine te0242(option, nomte)
!
    use FE_topo_module
    use FE_quadrature_module
    use FE_basis_module
    use FE_eval_module
!
    implicit none
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/jevech.h"
#include "asterfort/leverettIsotTher.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcdiff.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvarc.h"
#include "asterfort/utmess.h"
#include "FE_module.h"
!
    character(len=16) :: option, nomte
! ......................................................................
!    - FONCTION REALISEE:   OPTION : 'DIFF_ELGA' , 'HYGR_ELGA'
!
!      DIFF_ELGA : CALCUL DU OU DES COEFFICIENTS DE DIFFUSION (SECHAGE)
!      HYGR_ELGA : CALCUL DES GRANDEURS LIEES A L'HYGROMETRIE (SECHAGE)
!
!    - ARGUMENTS:
!        DONNEES:      OPTION       -->  OPTION DE CALCUL
!                      NOMTE        -->  NOM DU TYPE ELEMENT
!
! ......................................................................
!
    type(FE_Cell) :: FECell
    type(FE_Quadrature) :: FEQuadCell
    type(FE_basis) :: FEBasis
!
    integer(kind=8) :: nbres
    parameter(nbres=1)
    integer(kind=8) :: nbpar
    parameter(nbpar=1)
    integer(kind=8) :: icodre(nbres)
    real(kind=8) :: tpg, diff, difl, difv, hygr, pcap
    real(kind=8) :: sechpg, valres(nbres), valpar(nbpar)
    real(kind=8), pointer :: fieldOutGauss(:) => null()
    real(kind=8), pointer :: sechr(:) => null()

    integer(kind=8) ::  kp, nbcmp
    integer(kind=8) ::  imate, iret
    character(len=8) :: nompar(nbpar)
    character(len=16) :: rela_name, phenom, nomres(nbres)
    character(len=16), pointer :: compor(:) => null()
! ----------------------------------------------------------------------
    call FECell%init()
    call FEQuadCell%initCell(FECell, "RIGI")
    call FEBasis%initCell(FECell)
!
    call jevech('PMATERC', 'L', imate)
    call jevech('PSECHRR', 'L', vr=sechr)

    if (option == "DIFF_ELGA") then
        call jevech('PCOMPOR', 'L', vk16=compor)
        rela_name = compor(RELA_NAME)
        call jevech('PDIFFPG', 'E', vr=fieldOutGauss)
        nbcmp = 3
        do kp = 1, FEQuadCell%nbQuadPoints
            if (rela_name(1:5) .eq. 'SECH_') then
                sechpg = FEEvalFuncRScal(FEBasis, sechr, FEQuadCell%points_param(1:3, kp))
                call rcvarc(' ', 'TEMP', '+', 'RIGI', kp, 1, tpg, iret)
                if (iret .ne. 0) call utmess('F', 'THERMIQUE1_2')
                call rcdiff(zi(imate), rela_name, tpg, sechpg, diff, difl, difv)
                fieldOutGauss(nbcmp*(kp-1)+1) = diff
                fieldOutGauss(nbcmp*(kp-1)+2) = difl
                fieldOutGauss(nbcmp*(kp-1)+3) = difv
            else
                ASSERT(ASTER_FALSE)
            end if
        end do
    elseif (option == "HYGR_ELGA") then
        call rccoma(zi(imate), 'BETON_DESORP', 0, phenom, icodre(1))
        if (icodre(1) .ne. 0) then
            call utmess('A', 'CALCCHAMP1_3')
        else
            call jevech('PHYGRPG', 'E', vr=fieldOutGauss)
            nbcmp = 2
            do kp = 1, FEQuadCell%nbQuadPoints
                sechpg = FEEvalFuncRScal(FEBasis, sechr, FEQuadCell%points_param(1:3, kp))
                nomres(1) = 'FONC_DESORP'
                nompar(1) = 'SECH'
                valpar(1) = sechpg
                call rcvalb('RIGI', kp, 1, '+', zi(imate), &
                            ' ', 'BETON_DESORP', nbpar, nompar, valpar, &
                            1, nomres, valres, icodre, 0)

                if (icodre(1) .eq. 0) then
                    hygr = valres(1)
                    pcap = 0.d0
                else
                    call rcvarc(' ', 'TEMP', '+', 'RIGI', kp, 1, tpg, iret)
                    if (iret .ne. 0) call utmess('F', 'THERMIQUE1_2')
                    call leverettIsotTher(sechpg, tpg, zi(imate), hygr, pc_=pcap)
                end if
                fieldOutGauss(nbcmp*(kp-1)+1) = hygr
                fieldOutGauss(nbcmp*(kp-1)+2) = pcap
            end do
        end if
    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
