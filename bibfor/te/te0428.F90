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
subroutine te0428(option, nomte)
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystPara, compCoorSystPlate, &
                                isPlateTria, isPlateQuad
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/dkqrge.h"
#include "asterfort/dktrge.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/plate_type.h"
#include "asterfort/utpslg.h"
#include "asterfort/utpvgl.h"
#include "jeveux.h"
!
    character(len=16) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: DKT, DKTG
!
! Options: RIGI_GEOM
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: nno
    integer(kind=8) :: jvGeom, jvMatr
    real(kind=8) :: pgl(3, 3), xyzl(3, 4)
    real(kind=8) :: matrGeom(300)
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(option .eq. 'RIGI_GEOM')
    call elrefe_info(fami='RIGI', nno=nno)

! - Get plate parameters
    call getCara(plateCara, plateOrie)
    ASSERT(plateCara%type .eq. PLATE_DKT .or. plateCara%type .eq. PLATE_DKTG)

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Calculate the transformation: global coordinate system/intrinsic coordinate system
    call compCoorSystPara(plateCara, zr(jvGeom), pgl)

! - Compute coordinate system for plate
    call compCoorSystPlate(pgl, plateCara, plateOrie)

! - Change coordinates of displacements
    call utpvgl(nno, 3, pgl, zr(jvGeom), xyzl)
!
    if (option .eq. 'RIGI_GEOM') then
        if (isPlateTria(plateCara)) then
            call dktrge(plateCara, plateOrie, &
                        xyzl, matrGeom)
        elseif (isPlateQuad(plateCara)) then
            call dkqrge(plateCara, plateOrie, &
                        xyzl, matrGeom)
        else
            ASSERT(ASTER_FALSE)
        end if
        call jevech('PMATUUR', 'E', jvMatr)
        call utpslg(nno, 6, pgl, matrGeom, zr(jvMatr))
    else
        ASSERT(ASTER_FALSE)
    end if
!
end
