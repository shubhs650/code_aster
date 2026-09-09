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
subroutine te0452(option, nomte)
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystNone
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/excent.h"
#include "asterfort/jevech.h"
#include "asterfort/tecach.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: DKT/DKTG/DST/Q4G/Q4GG
!
! Options: EFGE_EXCENT
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: itab1(7), itab2(7), iret, jvEfgeIn, jvEfgeOut, lgcata
    integer(kind=8) :: nbpoin, nbcmp, ibid
    aster_logical :: lreel
    real(kind=8) :: excen
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(option .eq. 'EFGE_EXCENT')
    call getCara(plateCara, plateOrie)
    call compCoorSystNone(plateOrie)

    call tecach('ONO', 'PEFFONR', 'L', iret, nval=7, itab=itab1)
    if (iret .eq. 0) then
        lreel = .true.
        call tecach('OOO', 'PEFFOENR', 'E', ibid, nval=7, itab=itab2)
    else
        call tecach('ONO', 'PEFFONC', 'L', iret, nval=7, itab=itab1)
        if (iret .eq. 0) then
            lreel = .false.
            call tecach('OOO', 'PEFFOENC', 'E', ibid, nval=7, itab=itab2)
        else
            call tecach('ONO', 'PEFFOGR', 'L', iret, nval=7, itab=itab1)
            if (iret .eq. 0) then
                lreel = .true.
                call tecach('OOO', 'PEFFOEGR', 'E', ibid, nval=7, itab=itab2)
            else
                lreel = .false.
                call tecach('OOO', 'PEFFOGC', 'L', ibid, nval=7, itab=itab1)
                call tecach('OOO', 'PEFFOEGC', 'E', ibid, nval=7, itab=itab2)
            end if
        end if
    end if

! - Input field
    jvEfgeIn = itab1(1)
    nbpoin = itab1(3)
    lgcata = itab1(2)
    nbcmp = lgcata/nbpoin
    ASSERT(lgcata .eq. nbpoin*nbcmp)
    ASSERT(nbcmp .eq. 6 .or. nbcmp .eq. 8)

! - Output field
    jvEfgeOut = itab2(1)
    ASSERT(itab2(2) .eq. lgcata)
!
    excen = plateCara%offset
    call excent('MOY', excen, nbpoin, nbcmp, lreel, &
                zr(jvEfgeIn), zr(jvEfgeOut), zc(jvEfgeIn), zc(jvEfgeOut))
!
end subroutine
