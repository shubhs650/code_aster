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
subroutine dxefnt(plateCara, plateOrie, &
                  sigmTher)
!
    use plateGeom_module, only: isPlateQuad, isPlateTria
    use plate_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/dxmath.h"
#include "asterfort/jevech.h"
#include "asterfort/r8inir.h"
#include "asterfort/rccoma.h"
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
! --- EFFORTS GENERALISES D'ORIGINE THERMIQUE AUX NOEUDS
! --- POUR LES ELEMENTS COQUES A FACETTES PLANES :
! --- DST, DKT, DSQ, DKQ, Q4G DUS :
! ---  .A UN CHAMP DE TEMPERATURES SUR LE PLAN MOYEN DONNANT
! ---        DES EFFORTS DE MEMBRANE
! ---  .A UN GRADIENT DE TEMPERATURES DANS L'EPAISSEUR DE LA COQUE
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: npgh = 3
    character(len=10) :: elasKeyword
    real(kind=8) :: df(3, 3), dm(3, 3), dmf(3, 3)
    real(kind=8) :: tempSup(4), tempInf(4), tempMid(4), rbid
    integer(kind=8) :: nbLayer, kpgMid, kpgSup
    integer(kind=8) :: indith, ino, iretRefe, iretInf, iretSup, iretMid
    integer(kind=8) :: jvMaterc, nno
    real(kind=8) :: coe1, coe2, epais, tempRefe
!
! --------------------------------------------------------------------------------------------------
!
    iretRefe = 0
    iretInf = 0
    iretSup = 0
    iretMid = 0
    sigmTher = 0.d0

! - S'IL N'Y A PAS DE TEMPERATURE, IL N'Y A RIEN A CALCULER :
    call rcvarc(' ', 'TEMP', '+', 'NOEU', 1, 1, rbid, iretMid)
    if (iretMid .ne. 0) goto 30

! - Get plate properties
    epais = plateCara%thick
    nbLayer = plateCara%nbLayer

    call rcvarc(' ', 'TEMP', 'REF', 'NOEU', 1, 1, tempRefe, iretRefe)
!
    if (isPlateQuad(plateCara)) then
        nno = 4
    else if (isPlateTria(plateCara)) then
        nno = 3
    else
        ASSERT(ASTER_FALSE)
    end if
!
    kpgMid = (3*nbLayer+1)/2
    kpgSup = npgh*nbLayer
    do ino = 1, nno
        call rcvarc(' ', 'TEMP', '+', 'NOEU', ino, kpgMid, tempMid(ino), iretMid)
        call rcvarc(' ', 'TEMP', '+', 'NOEU', ino, 1, tempInf(ino), iretInf)
        call rcvarc(' ', 'TEMP', '+', 'NOEU', ino, kpgSup, tempSup(ino), iretSup)
    end do
!
    call jevech('PMATERC', 'L', jvMaterc)
    call rccoma(zi(jvMaterc), 'ELAS', 1, elasKeyword)
!
    if ((elasKeyword .eq. 'ELAS') .or. (elasKeyword .eq. 'ELAS_COQUE') .or. &
        (elasKeyword .eq. 'ELAS_COQMU') .or. (elasKeyword .eq. 'ELAS_GLRC')) then

! - CALCUL DES MATRICES DE HOOKE DE FLEXION, MEMBRANE, MEMBRANE-FLEXION
! - CISAILLEMENT, CISAILLEMENT INVERSE
        call dxmath(plateCara, plateOrie, &
                    'NOEU', nno, &
                    df, dm, dmf, &
                    indith)
        if (indith .ne. -1) then
            if (iretMid+iretInf+iretSup .eq. 0) then
                if (iretRefe .eq. 1) then
                    call utmess('F', 'COMPOR5_43')
                else
                    do ino = 1, nno
!  --      LES COEFFICIENTS SUIVANTS RESULTENT DE L'HYPOTHESE SELON
!  --      LAQUELLE LA TEMPERATURE EST PARABOLIQUE DANS L'EPAISSEUR.
!  --      ON NE PREJUGE EN RIEN DE LA NATURE DU MATERIAU.
!  --      CETTE INFORMATION EST CONTENUE DANS LES MATRICES QUI
!  --      SONT LES RESULTATS DE LA ROUTINE DXMATH
                        coe1 = (tempSup(ino)+tempInf(ino)+4.d0*tempMid(ino))/6.d0-tempRefe
                        coe2 = (tempSup(ino)-tempInf(ino))/epais
                        sigmTher(1+8*(ino-1)) = coe1*(dm(1, 1)+dm(1, 2))+coe2*(dmf(1, 1)+dmf(1, 2))
                        sigmTher(2+8*(ino-1)) = coe1*(dm(2, 1)+dm(2, 2))+coe2*(dmf(2, 1)+dmf(2, 2))
                        sigmTher(3+8*(ino-1)) = coe1*(dm(3, 1)+dm(3, 2))+coe2*(dmf(3, 1)+dmf(3, 2))
                        sigmTher(4+8*(ino-1)) = coe2*(df(1, 1)+df(1, 2))+coe1*(dmf(1, 1)+dmf(1, 2))
                        sigmTher(5+8*(ino-1)) = coe2*(df(2, 1)+df(2, 2))+coe1*(dmf(2, 1)+dmf(2, 2))
                        sigmTher(6+8*(ino-1)) = coe2*(df(3, 1)+df(3, 2))+coe1*(dmf(3, 1)+dmf(3, 2))
                    end do
                end if
            end if
        end if
    end if
30  continue
end subroutine
