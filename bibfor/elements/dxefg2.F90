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
subroutine dxefg2(plateCara, plateOrie, &
                  pgl, sigt)
!
    use plate_type
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/dxmat1.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/rcvarc.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    type(plateCara_Para), intent(in) :: plateCara
    type(plateOrie_Para), intent(in) :: plateOrie
    real(kind=8), intent(in) :: pgl(3, 3)
    real(kind=8), intent(out) :: sigt(32)
!
! --------------------------------------------------------------------------------------------------
!
! --- EFFORTS GENERALISES N, M, V D'ORIGINE THERMIQUE AUX POINTS
! --- D'INTEGRATION POUR LES ELEMENTS COQUES A FACETTES PLANES :
! --- DKTG DUS :
! ---  .A UN CHAMP DE TEMPERATURES SUR LE PLAN MOYEN DONNANT
! ---        DES EFFORTS DE MEMBRANE
! ---  .A UN GRADIENT DE TEMPERATURES DANS L'EPAISSEUR DE LA COQUE
!     ------------------------------------------------------------------
!     IN  PGL(3,3)     : MATRICE DE PASSAGE DU REPERE GLOBAL AU REPERE
!                        LOCAL
!     OUT SIGT(1)      : EFFORTS  GENERALISES D'ORIGINE THERMIQUE
!                        AUX POINTS D'INTEGRATION
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: npg
    real(kind=8) :: df(3, 3), dm(3, 3), dmf(3, 3)
    real(kind=8) :: tmoypg, tsuppg, tinfpg
    character(len=8), parameter :: fami = 'RIGI'
    integer(kind=8) :: kpg, indith, iret, iret1, ireti, irets, iretm
    real(kind=8) :: coe1, coe2, epais, tempRefe, rbid
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami=fami, npg=npg)

    sigt = 0.d0

! - RECUPERATION DE LA TEMPERATURE DE REFERENCE ET DE L'EPAISSEUR DE LA COQUE
    epais = plateCara%thick
    call rcvarc(' ', 'TEMP_MIL', 'REF', fami, 1, 1, tempRefe, iret1)

! - CALCUL DES COEFFICIENTS THERMOELASTIQUES DE FLEXION, MEMBRANE, MEMBRANE-FLEXION
    call dxmat1(plateOrie, &
                'RIGI', epais, df, dm, dmf, pgl, indith, npg)

    if (indith .ne. -1) then
        do kpg = 1, npg
!  --      TEMPERATURES SUR LES FEUILLETS MOYEN, SUPERIEUR ET INFERIEUR
            call rcvarc(' ', 'TEMP_INF', '+', fami, kpg, 1, tinfpg, ireti)
            call rcvarc(' ', 'TEMP_SUP', '+', fami, kpg, 1, tsuppg, irets)
            call rcvarc(' ', 'TEMP_MIL', '+', fami, kpg, 1, tmoypg, iretm)
            ASSERT(ireti .eq. irets)

!           -- si il n'existe ni TEMP_INF, ni TEMP_SUP :
            if (ireti .ne. 0) then
!               -- si on trouve 'TEMP' : c'est probablement une erreur d'utilisation :
                call rcvarc(' ', 'TEMP', '+', fami, kpg, 1, rbid, iret)
                if (iret .eq. 0) call utmess('F', 'CALCULEL3_18')
!               -- sinon, il n'y a rien a calculer
                ASSERT(kpg .eq. 1)
                goto 999
            end if

!           -- si on ne trouve pas TEMP_MIL, on prend la moyenne de TEM_INF te TEMP_SUP :
            if (iretm .ne. 0) then
                tmoypg = (tinfpg+tsuppg)/2.d0
            end if
!
            if (iret1 .eq. 1) then
                call utmess('F', 'COMPOR5_43')
            else
!  --          LES COEFFICIENTS SUIVANTS RESULTENT DE L'HYPOTHESE SELON
!  --          LAQUELLE LA TEMPERATURE EST PARABOLIQUE DANS L'EPAISSEUR.
!  --          ON NE PREJUGE EN RIEN DE LA NATURE DU MATERIAU.
!  --          CETTE INFORMATION EST CONTENUE DANS LES MATRICES QUI
!  --          SONT LES RESULTATS DE LA ROUTINE DXMATH
                coe1 = (tsuppg+tinfpg+4.d0*tmoypg)/6.d0-tempRefe
                coe2 = (tsuppg-tinfpg)/epais
                sigt(1+8*(kpg-1)) = coe1*(dm(1, 1)+dm(1, 2))+coe2*(dmf(1, 1)+dmf(1, 2))
                sigt(2+8*(kpg-1)) = coe1*(dm(2, 1)+dm(2, 2))+coe2*(dmf(2, 1)+dmf(2, 2))
                sigt(3+8*(kpg-1)) = coe1*(dm(3, 1)+dm(3, 2))+coe2*(dmf(3, 1)+dmf(3, 2))
                sigt(4+8*(kpg-1)) = coe2*(df(1, 1)+df(1, 2))+coe1*(dmf(1, 1)+dmf(1, 2))
                sigt(5+8*(kpg-1)) = coe2*(df(2, 1)+df(2, 2))+coe1*(dmf(2, 1)+dmf(2, 2))
                sigt(6+8*(kpg-1)) = coe2*(df(3, 1)+df(3, 2))+coe1*(dmf(3, 1)+dmf(3, 2))
            end if
        end do
    end if

999 continue

end subroutine
