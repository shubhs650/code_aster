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
subroutine plateFollowerPressure(option, nomte)
!
    use plate_type, only: plateOrie_Para
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dxqfor.h"
#include "asterfort/dxqpgl.h"
#include "asterfort/dxtfor.h"
#include "asterfort/dxtpgl.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe.h"
#include "asterfort/jevech.h"
#include "asterfort/tecach.h"
#include "asterfort/utpvgl.h"
#include "asterfort/utpvlg.h"
!
    character(len=16), intent(in) :: option, nomte
!
    integer(kind=8), parameter :: mxnoeu = 4, ndofbynode = 6
    integer(kind=8), parameter :: maxdof = mxnoeu*ndofbynode
    real(kind=8), parameter :: eps_fd = 1.d-7
!
    aster_logical :: l_func, global
    integer(kind=8) :: jv_geom, jv_time, jv_pres, jv_depm, jv_depp
    integer(kind=8) :: jv_vect, jv_matr
    integer(kind=8) :: nno, npg, ndim, nnos, jgano, ipoids, ivf, idfde
    integer(kind=8) :: iret, idof
    real(kind=8) :: time
    real(kind=8) :: depl_tot(maxdof), delta(maxdof), vect0(maxdof), vectp(maxdof), vectm(maxdof)
!
    l_func = (option .eq. 'CHAR_MECA_PRSU_F') .or. (option .eq. 'RIGI_MECA_PRSU_F')
    global = ASTER_FALSE
!
    call jevech('PGEOMER', 'L', jv_geom)
    call jevech('PDEPLMR', 'L', jv_depm)
    call jevech('PDEPLPR', 'L', jv_depp)
    if (l_func) then
        call jevech('PPRESSF', 'L', jv_pres)
    else
        call jevech('PPRESSR', 'L', jv_pres)
    end if
!
    time = 0.d0
    if (l_func) then
        call tecach('NNO', 'PINSTR', 'L', iret, iad=jv_time)
        if (jv_time .ne. 0) time = zr(jv_time)
    end if
!
    call elrefe_info(fami='RIGI', nno=nno, nnos=nnos, npg=npg, ndim=ndim, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
    ASSERT(ndim .eq. 2)
    ASSERT(nno .ge. 3 .and. nno .le. mxnoeu)
!
    depl_tot = 0.d0
    do idof = 1, nno*ndofbynode
        depl_tot(idof) = zr(jv_depm+idof-1) + zr(jv_depp+idof-1)
    end do
!
    call compute_follower_vector(nomte, nno, l_func, time, jv_geom, jv_pres, depl_tot, vect0)
!
    if (option(1:9) .eq. 'CHAR_MECA') then
        call jevech('PVECTUR', 'E', jv_vect)
        do idof = 1, nno*ndofbynode
            zr(jv_vect+idof-1) = vect0(idof)
        end do
    else if (option(1:9) .eq. 'RIGI_MECA') then
        call jevech('PMATUNS', 'E', jv_matr)
        do idof = 1, nno*ndofbynode
            delta = 0.d0
            delta(idof) = eps_fd
            call compute_follower_vector(nomte, nno, l_func, time, jv_geom, jv_pres, depl_tot+delta, vectp)
            call compute_follower_vector(nomte, nno, l_func, time, jv_geom, jv_pres, depl_tot-delta, vectm)
            call store_fd_column(nno, idof, eps_fd, vectp, vectm, zr(jv_matr))
        end do
    else
        ASSERT(ASTER_FALSE)
    end if
!
contains
!
    subroutine compute_follower_vector(nomte, nno, l_func, time, jv_geom, jv_pres, depl, vect_out)
        character(len=16), intent(in) :: nomte
        integer(kind=8), intent(in) :: nno, jv_geom, jv_pres
        aster_logical, intent(in) :: l_func
        real(kind=8), intent(in) :: time, depl(maxdof)
        real(kind=8), intent(out) :: vect_out(maxdof)
!
        integer(kind=8) :: inode, ier
        real(kind=8) :: geom_curr(3, mxnoeu), xyzl(3, mxnoeu)
        real(kind=8) :: pgl_geom(3, 3)
        real(kind=8) :: force_nodal(6, mxnoeu), vecl(maxdof)
        real(kind=8) :: press_nodal(mxnoeu), valpar(4), pr
        character(len=8) :: nompar(4)
        type(plateOrie_Para) :: plateOrie
!
        vect_out = 0.d0
        geom_curr = 0.d0
        xyzl = 0.d0
        force_nodal = 0.d0
        vecl = 0.d0
        press_nodal = 0.d0
!
        do inode = 1, nno
            geom_curr(1, inode) = zr(jv_geom+3*(inode-1)) + depl(6*(inode-1)+1)
            geom_curr(2, inode) = zr(jv_geom+3*(inode-1)+1) + depl(6*(inode-1)+2)
            geom_curr(3, inode) = zr(jv_geom+3*(inode-1)+2) + depl(6*(inode-1)+3)
        end do
!
        if (nno .eq. 3) then
            call dxtpgl(geom_curr, pgl_geom)
        else
            call dxqpgl(geom_curr, pgl_geom)
        end if
!
        call utpvgl(nno, 3, pgl_geom, geom_curr, xyzl)
!
        if (l_func) then
            nompar(1) = 'X'
            nompar(2) = 'Y'
            nompar(3) = 'Z'
            nompar(4) = 'INST'
            do inode = 1, nno
                valpar(1) = geom_curr(1, inode)
                valpar(2) = geom_curr(2, inode)
                valpar(3) = geom_curr(3, inode)
                valpar(4) = time
                call fointe('FM', zk8(jv_pres), 4, nompar, valpar, pr, ier)
                ASSERT(ier .eq. 0)
                press_nodal(inode) = pr
            end do
        else
            do inode = 1, nno
                press_nodal(inode) = zr(jv_pres+inode-1)
            end do
        end if
!
        do inode = 1, nno
            force_nodal(3, inode) = -press_nodal(inode)
        end do
!
        plateOrie%t2iu = [1.d0, 0.d0, 0.d0, 1.d0]
        if (nno .eq. 3) then
            call dxtfor(plateOrie, global, xyzl, force_nodal, vecl)
        else
            call dxqfor(plateOrie, global, xyzl, force_nodal, vecl)
        end if
        call utpvlg(nno, ndofbynode, pgl_geom, vecl, vect_out)
    end subroutine compute_follower_vector
!
    subroutine store_fd_column(nno, jdof, eps_fd, vectp, vectm, matr_uns)
        integer(kind=8), intent(in) :: nno, jdof
        real(kind=8), intent(in) :: eps_fd, vectp(maxdof), vectm(maxdof)
        real(kind=8), intent(inout) :: matr_uns(*)
        integer(kind=8) :: idof, k
        k = jdof
        do idof = 1, nno*ndofbynode
            matr_uns(k) = -(vectp(idof)-vectm(idof))/(2.d0*eps_fd)
            k = k + nno*ndofbynode
        end do
    end subroutine store_fd_column
!
end subroutine
