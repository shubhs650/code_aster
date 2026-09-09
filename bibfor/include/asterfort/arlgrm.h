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
!
interface
#include "MeshTypes_type.h"
    subroutine arlgrm(mail  ,nomgrp,dime  ,ima  ,connex,loncum, &
                      nummai,nommai,itypm ,nbno , cxno)
        integer(kind=8) :: ima
        integer(kind=8) :: dime
        integer(kind=8) :: connex(*)
        integer(kind=8) :: loncum(*)
        integer(kind=8) :: nummai
        integer(kind=8) :: nbno
        integer(kind=8) :: itypm
        integer(kind=8) :: cxno(MT_NNOMAX3D)
        character(len=8) :: mail
        character(len=8) :: nommai
        character(len=19) :: nomgrp
    end subroutine arlgrm
end interface
