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
subroutine postcoq3d(plateCara, plateOrie, &
                     optionZ, nomteZ, nbLayer)
!
    use plateGeom_module, only: updateCoorSystCO3D
    use plate_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/tecach.h"
#include "asterfort/vdefro.h"
#include "asterfort/vdsiro.h"
#include "asterfort/vdxedg.h"
#include "asterfort/vdxeps.h"
#include "asterfort/vdxsig.h"
#include "jeveux.h"
!
    type(plateCara_Para), intent(in) :: plateCara
    type(plateOrie_Para), intent(inout) :: plateOrie
    character(len=*), intent(in) :: optionZ, nomteZ
    integer(kind=8), intent(in) :: nbLayer
!
! --------------------------------------------------------------------------------------------------
!
! COQUE_3D
!
! Compute DEGE_ELGA, DEGE_ELNO, EPSI_ELGA, SIEF_ELGA
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
! In  nbLayer          : number of layers
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: npgt = 10
    integer(kind=8) :: jvSigm, jvGeom, lzi, lzr, jvDege, jvEpsi
    integer(kind=8) :: nb1, npgsn, npgsr
    character(len=16) :: option, nomte
    integer(kind=8) :: itab(7), iret
    real(kind=8) :: degeElga(72), degeElno(8, 9)
    real(kind=8) :: siefElga(6*27*nbLayer), epsiElga(6*27*nbLayer)
    real(kind=8) :: matevn(2, 2, npgt), matevg(2, 2, npgt)
!
! --------------------------------------------------------------------------------------------------
!
    option = optionZ
    nomte = nomteZ

! - Acces to geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Access to static objects of COQUE_3D
    call jevete('&INEL.'//nomte(1:8)//'.DESI', ' ', lzi)
    nb1 = zi(lzi-1+1)
    npgsr = zi(lzi-1+3)
    npgsn = zi(lzi-1+4)
    call jevete('&INEL.'//nomte(1:8)//'.DESR', ' ', lzr)
    ASSERT(npgsr .le. 10)

! - Compute option
    if (option .eq. 'DEGE_ELGA' .or. option .eq. 'DEGE_ELNO') then
        call vdxedg(plateCara, plateOrie, &
                    nomte, option, zr(jvGeom), &
                    degeElga, degeElno)

    elseif (option .eq. 'SIEF_ELGA') then
        call vdxsig(plateCara, plateOrie, &
                    nomte, zr(jvGeom), &
                    nbLayer, siefElga)

    elseif (option .eq. 'EPSI_ELGA') then
        call vdxeps(plateCara, plateOrie, &
                    nomte, zr(jvGeom), &
                    nbLayer, epsiElga)

    else
        ASSERT(ASTER_FALSE)
    end if

! - Update basis
    call updateCoorSystCO3D(plateCara, plateOrie, &
                            nomte, zr(jvGeom), &
                            npgsr, nb1, lzr, &
                            matevn, matevg)

    if (option .eq. 'EPSI_ELGA') then
        call tecach('OOO', 'PDEFOPG', 'E', iret, nval=7, itab=itab)
        jvEpsi = itab(1)
        call vdsiro(itab(3), itab(7), matevg, 'IU', 'G', &
                    epsiElga, zr(jvEpsi))

    else if (option .eq. 'SIEF_ELGA') then
        call tecach('OOO', 'PCONTRR', 'E', iret, nval=7, itab=itab)
        jvSigm = itab(1)
        call vdsiro(itab(3), itab(7), matevg, 'IU', 'G', &
                    siefElga, zr(jvSigm))

    else if (option(1:9) .eq. 'DEGE_ELGA') then
        call tecach('OOO', 'PDEFOPG', 'E', iret, nval=7, itab=itab)
        jvDege = itab(1)
        call vdefro(npgsn, matevn, degeElga, zr(jvDege))

    else if (option(1:9) .eq. 'DEGE_ELNO') then
        call tecach('OOO', 'PDEFOGR', 'E', iret, nval=7, itab=itab)
        jvDege = itab(1)
        call vdefro((nb1+1), matevn, degeElno, zr(jvDege))

    end if
!
end subroutine
