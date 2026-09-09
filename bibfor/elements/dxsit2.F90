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
subroutine dxsit2(plateCara, plateOrie, &
                  sigma)
!
    use plate_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/dxmat2.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/plate_type.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvarc.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    type(plateCara_Para), intent(in) :: plateCara
    type(plateOrie_Para), intent(in) :: plateOrie
    real(kind=8) :: sigma(*)
!
! --------------------------------------------------------------------------------------------------
!
!       CALCUL LES CONTRAINTES VRAIES AUX POINTS DE GAUSS
!       RETRANCHE LES CONTRAINTES PLANES D'ORIGINE THERMIQUE AUX POINTS
!       DE GAUSS, POUR LES ELEMENTS COQUES A FACETTES PLANES :
!       DST, DKT, DSQ, DKQ, Q4G DUS :
!       .A UN CHAMP DE TEMPERATURES MOYEN ET
!       .A UN GRADIENT DE TEMPERATURES DANS L'EPAISSEUR DE LA COQUE
!       DANS LE CAS ELASTIQUE ISOTROPE HOMOGENE
!       CAS ELAS_COQMU
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbcmp = 6
    character(len=8), parameter :: fami = "RIGI"
    integer(kind=8) :: npg
    integer(kind=8) :: iret1, iret2, iret3, iret4, iret5
    integer(kind=8) :: iLayer, nbLayer, ipg, igauh, npgh, icpg, imoy
    integer(kind=8) :: jvMaterc
    real(kind=8) :: dm(3, 3), tempRefe
    real(kind=8) :: tinf(4), tmoy(4), tsup(4)
    real(kind=8) :: ordi, epi, epais, coe1, coe2
    character(len=10) :: elasKeyword
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami=fami, npg=npg)
!
    iret1 = 0
    iret2 = 0
    iret3 = 0
    iret4 = 0
    iret5 = 0

! - Get layers of shell
    nbLayer = plateCara%nbLayer
    if (plateCara%type .eq. PLATE_DKTG) then
        ASSERT(nbLayer .eq. 1)
        npgh = 1
    else
        npgh = 3
    end if
    ASSERT(nbLayer .ge. 1)
    imoy = (3*nbLayer+1)/2

! - RECUPERATION DE LA TEMPERATURE DE REFERENCE
    call rcvarc(' ', 'TEMP', 'REF', 'RIGI', 1, &
                1, tempRefe, iret1)

!   S'IL N'Y A PAS DE TEMPERATURE DE REFERENCE, ON NE FAIT RIEN
    if (iret1 .eq. 1) goto 999

! - RECUPERATION DE LA TEMPERATURE SUR LES FEUILLETS
    do ipg = 1, npg
        call rcvarc(' ', 'TEMP', '+', 'RIGI', ipg, &
                    1, tinf(ipg), iret2)
        call rcvarc(' ', 'TEMP', '+', 'RIGI', ipg, &
                    imoy, tmoy(ipg), iret3)
        call rcvarc(' ', 'TEMP', '+', 'RIGI', ipg, &
                    3*nbLayer, tsup(ipg), iret4)
        iret5 = iret5+iret2+iret3+iret4
    end do
!
    call jevech('PMATERC', 'L', jvMaterc)
    call rccoma(zi(jvMaterc), 'ELAS', 1, elasKeyword)
    if ((elasKeyword .eq. 'ELAS') .or. (elasKeyword .eq. 'ELAS_ISTR') .or. &
        (elasKeyword .eq. 'ELAS_ORTH') .or. (elasKeyword .eq. 'ELAS_COQUE')) then
        call utmess('F', 'ELEMENTS_52', sk=elasKeyword(1:10))
    end if

! - CALCUL DES MATRICES DE HOOKE DE FLEXION, MEMBRANE,
! - MEMBRANE-FLEXION, CISAILLEMENT, CISAILLEMENT INVERSE
    do ipg = 1, npg
        do iLayer = 1, nbLayer
            do igauh = 1, npgh
                icpg = nbcmp*npgh*nbLayer*(ipg-1)+ &
                       nbcmp*npgh*(iLayer-1)+ &
                       nbcmp*(igauh-1)

                call dxmat2(plateCara, plateOrie, &
                            iLayer, npg, &
                            ordi, epi, &
                            epais, dm)

                if (iret5 .eq. 0) then
                    if (iret1 .eq. 1) then
                        call utmess('F', 'CALCULEL_15')
                    else
!  --      LES COEFFICIENTS SUIVANTS RESULTENT DE L'HYPOTHESE SELON
!  --      LAQUELLE LA TEMPERATURE EST PARABOLIQUE DANS L'EPAISSEUR.
!  --      LES COEFFICIENTS THERMOELASTIQUES PROVIENNENT DES
!  --      MATRICES QUI SONT LES RESULTATS DE LA ROUTINE DXMATL
                        coe1 = (tsup(ipg)+tinf(ipg)+4.d0*tmoy(ipg))/6.d0-tempRefe
                        coe2 = (tsup(ipg)-tinf(ipg))*(ordi+dble(igauh-2)*epi/2.d0)/epais
                        sigma(1+icpg) = sigma(1+icpg)-((dm(1, 1)+dm(1, 2))/epi)*(coe1+coe2)
                        sigma(2+icpg) = sigma(2+icpg)-((dm(2, 1)+dm(2, 2))/epi)*(coe1+coe2)
                        sigma(4+icpg) = sigma(4+icpg)-((dm(3, 1)+dm(3, 2))/epi)*(coe1+coe2)
                    end if
                end if
            end do
        end do
    end do
!
999 continue
!
end subroutine
