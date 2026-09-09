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
subroutine nmecrr(sddisc, paraNameZ, paraValeR_, paraValeI_)
!
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/jeveuo.h"
#include "jeveux.h"
!
    character(len=19), intent(in) :: sddisc
    character(len=*), intent(in) :: paraNameZ
    integer(kind=8), optional, intent(in) :: paraValeI_
    real(kind=8), optional, intent(in) :: paraValeR_
!
! --------------------------------------------------------------------------------------------------
!
! Non-linear algorithm - Discretization management
!
! Management of parameters - Write
!
! --------------------------------------------------------------------------------------------------
!
! In  sddisc           : datastructure for time discretization
! In  paraName         : parameter to manage
!   1 MXITER               : MAX( ITER_GLOB_MAXI , ITER_GLOB_ELAS )
!   2 MNITER               : MIN( ITER_GLOB_MAXI , ITER_GLOB_ELAS )
!   3 NBITER               : NOMBRE MAX ITERATIONS (Y COMPRIS EXTRAPOL)
!   4 PAS_MINI_ELAS        : PAS_MINI_ELAS
!   5 RESI_GLOB_RELA       : RESI_GLOB_RELA DONNE
!   6 RESI_GLOB_MAXI       : RESI_GLOB_MAXI
!   7 TYPE_RESI            :  =1 ON A DONNE RESI_GLOB_RELA
!                             =2 ON A DONNE RESI_GLOB_MAXI
!                             =3 C'EST (1) ET (2)
!                             =0 ON A RIEN DONNE ==> =1 (DEFAUT)
!   8 INIT_NEWTON_KRYLOV   : RESIDU INITIAL POUR NEWTON KRYLOV
!   9 ITER_NEWTON_KRYLOV   : RESIDU COURANT POUR NEWTON KRYLOV
!  10 ITERSUP              : =3 ON AUTORISE DES ITERATIONS EN PLUS
! Out paraValeI        : parameter (integer)
! Out paraValeR        : parameter (real)
!
! --------------------------------------------------------------------------------------------------
!
    character(len=24) :: paraName
    integer(kind=8):: paraValeI
    real(kind=8) :: paraValeR
    character(len=24) :: sddiscIfcvName
    real(kind=8), pointer :: sddiscIfcv(:) => null()
!
! --------------------------------------------------------------------------------------------------
!
    paraName = paraNameZ

! - Access to datastructure
    sddiscIfcvName = sddisc(1:19)//'.IFCV'
    call jeveuo(sddiscIfcvName, 'E', vr=sddiscIfcv)

! - Set
    paraValeI = 0
    paraValeR = 0.d0
    if (present(paraValeI_)) then
        paraValeI = paraValeI_
    end if
    if (present(paraValeR_)) then
        paraValeR = paraValeR_
    end if

    if (paraName .eq. 'MXITER') then
        sddiscIfcv(1) = paraValeI

    else if (paraName .eq. 'MNITER') then
        sddiscIfcv(2) = paraValeI

    else if (paraName .eq. 'NBITER') then
        sddiscIfcv(3) = paraValeI

    else if (paraName .eq. 'PAS_MINI_ELAS') then
        sddiscIfcv(4) = paraValeR

    else if (paraName .eq. 'RESI_GLOB_RELA') then
        sddiscIfcv(5) = paraValeR

    else if (paraName .eq. 'RESI_GLOB_MAXI') then
        sddiscIfcv(6) = paraValeR

    else if (paraName .eq. 'TYPE_RESI') then
        sddiscIfcv(7) = paraValeI

    else if (paraName .eq. 'INIT_NEWTON_KRYLOV') then
        sddiscIfcv(8) = paraValeR

    else if (paraName .eq. 'ITER_NEWTON_KRYLOV') then
        sddiscIfcv(9) = paraValeR

    else if (paraName .eq. 'ITERSUP') then
        sddiscIfcv(10) = paraValeI

    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
