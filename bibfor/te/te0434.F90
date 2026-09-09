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
subroutine te0434(option, nomte)
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystMemb
    implicit none
!
#include "asterc/r8dgrd.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/mbgchg.h"
#include "asterfort/mbxchg.h"
#include "asterfort/rccoma.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    character(len=16) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: MEMBRANE
!
! Options: CHAR_MECA_TEMP*, CHAR_MECA_EPSI*, FORC_NODA, REFE_FORC_NODA
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = 'RIGI'
    integer(kind=8), parameter :: nddl = 3, ncomp = 3
    character(len=32) :: elasKeyword
    integer(kind=8) :: nno, npg
    integer(kind=8) :: n, kpg
    integer(kind=8) :: ipoids, ivf, idfde, iret, jvCompor, itab(1), jvInst
    integer(kind=8) :: jvGeom, jvMaterc, jvSief, jvPesa, jvEpsi, jvVect
    integer(kind=8) :: icodre1, icodre2
    real(kind=8) :: dff(2, 9), vff(9)
    real(kind=8) :: h, preten
    aster_logical :: lGravity
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!
    lGravity = (option .eq. 'CHAR_MECA_PESA_R')

! - Get plate parameters
    call getCara(plateCara, plateOrie)

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Compute global<=>local transformation
    call compCoorSystMemb(plateOrie)

    call elrefe_info(fami=fami, nno=nno, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde)
!
! - Input fields
    call tecach('N', 'PCOMPOR', 'L', iret, 1, itab)
    jvCompor = itab(1)
    if (option .eq. 'FORC_NODA') then
        call jevech('PSIEFR', 'L', jvSief)
        call jevech('PMATERC', 'L', jvMaterc)
    else if (option .eq. 'REFE_FORC_NODA') then
        call jevech('PMATERC', 'L', jvMaterc)
    else if (option .eq. 'CHAR_MECA_EPSI_R') then
        call jevech('PMATERC', 'L', jvMaterc)
        call jevech('PEPSINR', 'L', jvEpsi)
    else if (option .eq. 'CHAR_MECA_EPSI_F') then
        call jevech('PMATERC', 'L', jvMaterc)
        call jevech('PEPSINF', 'L', jvEpsi)
        call jevech('PINSTR', 'L', jvInst)
    else if (lGravity) then
        call jevech('PMATERC', 'L', jvMaterc)
        call jevech('PPESANR', 'L', jvPesa)
    else if (option .eq. 'CHAR_MECA_TEMP_R') then
        call jevech('PMATERC', 'L', jvMaterc)
    else
        ASSERT(ASTER_FALSE)
    end if

! - Output field
    call jevech('PVECTUR', 'E', jvVect)

! - EPAISSEUR ET PRETCONTRAINTES
    h = plateCara%thick
    if (h .lt. r8prem()) then
        call utmess('F', 'MEMBRANE_1')
    end if
    preten = plateCara%tension/h

! -----------------------------------------------------------------
! ---  VERIFICATION DE LA CORRESPONDANCE MATERIAU / COMPORTMENT ---
! -----------------------------------------------------------------
    call rccoma(zi(jvMaterc), 'ELAS_MEMBRANE', 0, elasKeyword, icodre1)
    call rccoma(zi(jvMaterc), 'ELAS', 0, elasKeyword, icodre2)
    if (icodre1 .eq. 0) then
        if ((jvCompor .ne. 0) .and. (zk16(jvCompor+2) (1:5) .ne. 'PETIT')) then
            call utmess('F', 'MEMBRANE_10')
        end if
    elseif (icodre2 .eq. 0) then
        if (((jvCompor .eq. 0) .or. (zk16(jvCompor+2) (1:9) .ne. 'GROT_GDEP')) &
            .and. (.not. lGravity)) then
            call utmess('F', 'MEMBRANE_10')
        end if
    end if

    do kpg = 1, npg
        do n = 1, nno
            vff(n) = zr(ivf+(kpg-1)*nno+n-1)
            dff(1, n) = zr(idfde+(kpg-1)*nno*2+(n-1)*2)
            dff(2, n) = zr(idfde+(kpg-1)*nno*2+(n-1)*2+1)
        end do

        if (icodre1 .eq. 0) then
            call mbxchg(plateOrie, &
                        option, fami, &
                        nddl, nno, ncomp, kpg, npg, &
                        jvEpsi, jvInst, ipoids, jvGeom, &
                        jvMaterc, jvPesa, jvVect, jvSief, &
                        vff, dff)

        elseif (icodre2 .eq. 0) then
            if ((option .ne. 'FORC_NODA') .and. (option .ne. 'CHAR_MECA_PESA_R')) then
                call utmess('F', 'MEMBRANE_7')
            end if
            call mbgchg(plateOrie, &
                        option, fami, &
                        nddl, nno, ncomp, kpg, &
                        jvMaterc, jvSief, &
                        ipoids, jvPesa, jvGeom, jvVect, &
                        vff, dff, &
                        h, preten)
        end if
    end do
!
end subroutine
