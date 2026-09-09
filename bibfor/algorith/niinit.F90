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
subroutine niinit(typmod, &
                  ndim, nno1, nno2, nno3, nno4, &
                  vu, vg, vp, vpi)
!
    implicit none
!
#include "asterfort/assert.h"
#include "MeshTypes_type.h"
!
    character(len=8), intent(in) :: typmod(*)
    integer(kind=8), intent(in) :: ndim, nno1, nno2, nno3, nno4
    integer(kind=8), intent(out) :: vu(3, MT_NNOMAX), vg(MT_NNOMAX)
    integer(kind=8), intent(out) :: vp(MT_NNOMAX), vpi(3, MT_NNOMAX)
!
! --------------------------------------------------------------------------------------------------
!
! Elements - Mixed formulation
!
! Get index of dof
!
! --------------------------------------------------------------------------------------------------
!
! In  typmod           : type of modelization
! In  ndim             : dimension of element (2 or 3)
! In  nno1             : number of nodes for displacements
! In  nno2             : number of nodes for GONF
! In  nno3             : number of nodes for PRES
! In  nno4             : number of nodes for GRAD_PRES
! Out vu               : index of dof for displacements
! Out vg               : index of dof for GONF
! Out vp               : index of dof for PRES
! Out vpi              : index of dof for GRAD_PRES
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i_node, os, iefm
!
! --------------------------------------------------------------------------------------------------
!
    vu(:, :) = 0
    vg(:) = 0
    vp(:) = 0
    vpi(:, :) = 0
!
    if (nno2 .ne. 0) then
        if (nno1 .eq. nno3) then
! --------- 2 = P2-P1-P2
            iefm = 2
        else if (nno2 .eq. nno3) then
! --------- 1 = P2-P1-P1
            iefm = 1
        end if
    else
        if (nno4 .eq. 0) then
            if (nno1 .eq. nno3) then
! ------------- 7 = P1+-P1
                iefm = 7
            else
! ------------- 8 = P2-P1
                iefm = 8
            end if
        else
! --------- 9 = P1 OSGS-P1
            iefm = 9
        end if
    end if
!
    if (ndim .eq. 3) then
        if (iefm .eq. 1) then
            do i_node = 1, nno2
                vu(1, i_node) = 1+(i_node-1)*5
                vu(2, i_node) = 2+(i_node-1)*5
                vu(3, i_node) = 3+(i_node-1)*5
                vp(i_node) = 4+(i_node-1)*5
                vg(i_node) = 5+(i_node-1)*5
            end do
            os = 5*nno2
            do i_node = 1, nno1-nno2
                vu(1, i_node+nno2) = 1+(i_node-1)*3+os
                vu(2, i_node+nno2) = 2+(i_node-1)*3+os
                vu(3, i_node+nno2) = 3+(i_node-1)*3+os
            end do
        else if (iefm .eq. 2) then
            do i_node = 1, nno2
                vu(1, i_node) = 1+(i_node-1)*5
                vu(2, i_node) = 2+(i_node-1)*5
                vu(3, i_node) = 3+(i_node-1)*5
                vp(i_node) = 4+(i_node-1)*5
                vg(i_node) = 5+(i_node-1)*5
            end do
            os = 5*nno2
            do i_node = 1, nno1-nno2
                vu(1, i_node+nno2) = 1+(i_node-1)*4+os
                vu(2, i_node+nno2) = 2+(i_node-1)*4+os
                vu(3, i_node+nno2) = 3+(i_node-1)*4+os
                vp(i_node+nno2) = 4+(i_node-1)*4+os
            end do
        else if (iefm .eq. 7) then
            do i_node = 1, nno1
                vu(1, i_node) = 1+(i_node-1)*4
                vu(2, i_node) = 2+(i_node-1)*4
                vu(3, i_node) = 3+(i_node-1)*4
                vp(i_node) = 4+(i_node-1)*4
            end do
        else if (iefm .eq. 8) then
            do i_node = 1, nno3
                vu(1, i_node) = 1+(i_node-1)*4
                vu(2, i_node) = 2+(i_node-1)*4
                vu(3, i_node) = 3+(i_node-1)*4
                vp(i_node) = 4+(i_node-1)*4
            end do
            os = 4*nno3
            do i_node = 1, nno1-nno3
                vu(1, i_node+nno3) = 1+(i_node-1)*3+os
                vu(2, i_node+nno3) = 2+(i_node-1)*3+os
                vu(3, i_node+nno3) = 3+(i_node-1)*3+os
            end do
        else if (iefm .eq. 9) then
            do i_node = 1, nno1
                vu(1, i_node) = 1+(i_node-1)*7
                vu(2, i_node) = 2+(i_node-1)*7
                vu(3, i_node) = 3+(i_node-1)*7
                vp(i_node) = 4+(i_node-1)*7
                vpi(1, i_node) = 5+(i_node-1)*7
                vpi(2, i_node) = 6+(i_node-1)*7
                vpi(3, i_node) = 7+(i_node-1)*7
            end do
        else
            ASSERT(ASTER_FALSE)
        end if
    else if (ndim .eq. 2) then
        if (iefm .eq. 1) then
            do i_node = 1, nno2
                vu(1, i_node) = 1+(i_node-1)*4
                vu(2, i_node) = 2+(i_node-1)*4
                vu(3, i_node) = 0
                vp(i_node) = 3+(i_node-1)*4
                vg(i_node) = 4+(i_node-1)*4
            end do
            os = 4*nno2
            do i_node = 1, nno1-nno2
                vu(1, i_node+nno2) = 1+(i_node-1)*2+os
                vu(2, i_node+nno2) = 2+(i_node-1)*2+os
                vu(3, i_node) = 0
            end do
        else if (iefm .eq. 2) then
            do i_node = 1, nno2
                vu(1, i_node) = 1+(i_node-1)*4
                vu(2, i_node) = 2+(i_node-1)*4
                vu(3, i_node) = 0
                vp(i_node) = 3+(i_node-1)*4
                vg(i_node) = 4+(i_node-1)*4
            end do
            os = 4*nno2
            do i_node = 1, nno1-nno2
                vu(1, i_node+nno2) = 1+(i_node-1)*3+os
                vu(2, i_node+nno2) = 2+(i_node-1)*3+os
                vu(3, i_node) = 0
                vp(i_node+nno2) = 3+(i_node-1)*3+os
            end do
        else if (iefm .eq. 7) then
            do i_node = 1, nno1
                vu(1, i_node) = 1+(i_node-1)*3
                vu(2, i_node) = 2+(i_node-1)*3
                vu(3, i_node) = 0
                vp(i_node) = 3+(i_node-1)*3
            end do
        else if (iefm .eq. 8) then
            do i_node = 1, nno3
                vu(1, i_node) = 1+(i_node-1)*3
                vu(2, i_node) = 2+(i_node-1)*3
                vu(3, i_node) = 0
                vp(i_node) = 3+(i_node-1)*3
            end do
            os = 3*nno3
            do i_node = 1, nno1-nno3
                vu(1, i_node+nno3) = 1+(i_node-1)*2+os
                vu(2, i_node+nno3) = 2+(i_node-1)*2+os
                vu(3, i_node+nno3) = 0
            end do
        else if (iefm .eq. 9) then
            do i_node = 1, nno1
                vu(1, i_node) = 1+(i_node-1)*5
                vu(2, i_node) = 2+(i_node-1)*5
                vu(3, i_node) = 0
                vp(i_node) = 3+(i_node-1)*5
                vpi(1, i_node) = 4+(i_node-1)*5
                vpi(2, i_node) = 5+(i_node-1)*5
                vpi(3, i_node) = 0
            end do
        else
            ASSERT(ASTER_FALSE)
        end if
    else
        ASSERT(ASTER_FALSE)
    end if
!
    if (typmod(1) .eq. 'AXIS') then
        do i_node = 1, nno1
            vu(3, i_node) = vu(1, i_node)
        end do
    end if
!
end subroutine
