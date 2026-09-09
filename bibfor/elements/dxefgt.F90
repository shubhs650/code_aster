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
subroutine dxefgt(plateCara, plateOrie, &
                  sigmTher)
!
    use plate_type
    implicit none
!
#include "asterfort/dxmath.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/rcvarc.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    type(plateCara_Para), intent(in) :: plateCara
    type(plateOrie_Para), intent(in) :: plateOrie
    real(kind=8), intent(out) :: sigmTher(32)
!
! --------------------------------------------------------------------------------------------------
!
! --- EFFORTS GENERALISES N, M, V D'ORIGINE THERMIQUE AUX POINTS
! --- D'INTEGRATION POUR LES ELEMENTS COQUES A FACETTES PLANES :
! --- DST, DKT, DSQ, DKQ, Q4G DUS :
! ---  .A UN CHAMP DE TEMPERATURES SUR LE PLAN MOYEN DONNANT
! ---        DES EFFORTS DE MEMBRANE
! ---  .A UN GRADIENT DE TEMPERATURES DANS L'EPAISSEUR DE LA COQUE
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: npgh = 3
    integer(kind=8) :: npg, nbLayer
    integer(kind=8) :: kpgMid, kpgSup
    real(kind=8) :: df(3, 3), dm(3, 3), dmf(3, 3)
    real(kind=8) :: tempMidKpg, tempSupKpg, tempInfKpg
    character(len=8), parameter :: fami = "RIGI"
    integer(kind=8) :: kpg, indith, iretRefe, iretMid, iretInf, iretSup
    real(kind=8) :: coe1, coe2, epais, tempRefe
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami=fami, npg=npg)
    sigmTher = 0.d0

! - RECUPERATION DE LA TEMPERATURE DE REFERENCE
    call rcvarc(' ', 'TEMP', 'REF', fami, 1, 1, tempRefe, iretRefe)

! - CALCUL DES COEFFICIENTS THERMOELASTIQUES DE FLEXION,MEMBRANE, MEMBRANE-FLEXION
    call dxmath(plateCara, plateOrie, &
                'RIGI', npg, &
                df, dm, dmf, &
                indith)

! - Get plate properties
    epais = plateCara%thick
    nbLayer = plateCara%nbLayer

    if (indith .ne. -1) then
        kpgMid = (3*nbLayer+1)/2
        kpgSup = npgh*nbLayer

        do kpg = 1, npg
!  --      TEMPERATURES SUR LES FEUILLETS MOYEN, SUPERIEUR ET INFERIEUR
            call rcvarc(' ', 'TEMP', '+', fami, kpg, kpgMid, tempMidKpg, iretMid)
            call rcvarc(' ', 'TEMP', '+', fami, kpg, 1, tempInfKpg, iretInf)
            call rcvarc(' ', 'TEMP', '+', fami, kpg, kpgSup, tempSupKpg, iretSup)
            if (iretMid+iretInf+iretSup .eq. 0) then
                if (iretRefe .eq. 1) then
                    call utmess('F', 'COMPOR5_43')
                else

!  --      LES COEFFICIENTS SUIVANTS RESULTENT DE L'HYPOTHESE SELON
!  --      LAQUELLE LA TEMPERATURE EST PARABOLIQUE DANS L'EPAISSEUR.
!  --      ON NE PREJUGE EN RIEN DE LA NATURE DU MATERIAU.
!  --      CETTE INFORMATION EST CONTENUE DANS LES MATRICES QUI
!  --      SONT LES RESULTATS DE LA ROUTINE DXMATH
                    coe1 = (tempSupKpg+tempInfKpg+4.d0*tempMidKpg)/6.d0-tempRefe
                    coe2 = (tempSupKpg-tempInfKpg)/epais
                    sigmTher(1+8*(kpg-1)) = coe1*(dm(1, 1)+dm(1, 2))+coe2*(dmf(1, 1)+dmf(1, 2))
                    sigmTher(2+8*(kpg-1)) = coe1*(dm(2, 1)+dm(2, 2))+coe2*(dmf(2, 1)+dmf(2, 2))
                    sigmTher(3+8*(kpg-1)) = coe1*(dm(3, 1)+dm(3, 2))+coe2*(dmf(3, 1)+dmf(3, 2))
                    sigmTher(4+8*(kpg-1)) = coe2*(df(1, 1)+df(1, 2))+coe1*(dmf(1, 1)+dmf(1, 2))
                    sigmTher(5+8*(kpg-1)) = coe2*(df(2, 1)+df(2, 2))+coe1*(dmf(2, 1)+dmf(2, 2))
                    sigmTher(6+8*(kpg-1)) = coe2*(df(3, 1)+df(3, 2))+coe1*(dmf(3, 1)+dmf(3, 2))
                end if
            end if
        end do
    end if
end subroutine
