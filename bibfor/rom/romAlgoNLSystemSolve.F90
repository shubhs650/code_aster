! --------------------------------------------------------------------
! Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
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
! person_in_charge: mickael.abbas at edf.fr
!
subroutine romAlgoNLSystemSolve(matr_asse, vect_2mbr, vect_cine, ds_algorom, vect_solu, &
                                l_update_redu_)
!
    use Rom_Datastructure_type
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/infniv.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mgauss.h"
#include "asterfort/mrmult.h"
#include "asterfort/jelira.h"
#include "asterfort/rsexch.h"
#include "asterfort/vtaxpy.h"
#include "asterfort/vtzero.h"
#include "asterfort/csmbgg.h"
#include "asterfort/mrconl.h"
#include "asterfort/mtmchc.h"
#include "asterfort/utmess.h"
#include "blas/ddot.h"
!
    character(len=24), intent(in) :: matr_asse
    character(len=24), intent(in) :: vect_2mbr
    character(len=24), intent(in) :: vect_cine
    type(ROM_DS_AlgoPara), intent(in) :: ds_algorom
    character(len=19), intent(in) :: vect_solu
    aster_logical, optional, intent(in) :: l_update_redu_
!
! --------------------------------------------------------------------------------------------------
!
! Model reduction - Solving non-linear problem
!
! Solve reduced system
!
! --------------------------------------------------------------------------------------------------
!
! In  matr_asse        : matrix
! In  vect_2mbr        : second member
! In  vect_cine        : vector for AFFE_CHAR_CINE load
! In  ds_algorom       : datastructure for ROM parameters
! In  vect_solu        : solution
! In  l_update_redu    : flag for update reduced coordinates
!
! --------------------------------------------------------------------------------------------------
!
    integer :: ifm, niv
    character(len=24) :: gamma, fieldName
    real(kind=8), pointer :: v_gamma(:) => null()
    real(kind=8), pointer :: v_vect_2mbr(:) => null()
    integer :: nbEqua_2mbr, nbEqua_matr, nbEqua, nbMode
    integer :: iMode, jMode, iEqua
    integer :: jv_matr, iret
    aster_logical :: l_hrom, l_rom, l_update_redu
    character(len=8) :: resultName
    complex(kind=8) :: cbid
    character(len=19) :: mode, vcine19
    real(kind=8) :: term1, term2, det, term
    real(kind=8), pointer :: v_matr_rom(:) => null()
    real(kind=8), pointer :: v_vect_rom(:) => null()
    real(kind=8), pointer :: v_mrmult(:) => null()
    real(kind=8), pointer :: v_mode(:) => null()
    real(kind=8), pointer :: v_vect_cine(:) => null()
    real(kind=8), pointer :: v_vect_solu(:) => null()
    character(len=24), pointer :: refa(:) => null()
    character(len=16) :: syme
    blas_int :: b_incx, b_incy, b_n
    cbid = dcmplx(0.d0, 0.d0)
!
! --------------------------------------------------------------------------------------------------
!
    call infniv(ifm, niv)
    if (niv .ge. 2) then
        call utmess('I', 'ROM5_40')
    end if
!
! - Get parameters
!
    l_rom = ds_algorom%l_rom
    l_hrom = ds_algorom%l_hrom
    gamma = ds_algorom%gamma
    resultName = ds_algorom%ds_empi%resultName
    nbMode = ds_algorom%ds_empi%nbMode
    nbEqua = ds_algorom%ds_empi%mode%nbEqua
    fieldName = ds_algorom%ds_empi%mode%fieldName
    vcine19 = vect_cine(1:19)
    ASSERT(ds_algorom%ds_empi%mode%fieldSupp .eq. 'NOEU')
    ASSERT(l_rom)
    l_update_redu = ASTER_TRUE
    if (present(l_update_redu_)) then
        l_update_redu = l_update_redu_
    end if
!
! - Access to reduced coordinates
!
    call jeveuo(gamma, 'E', vr=v_gamma)
!
! - Access to second member
!
    call jeveuo(vect_2mbr(1:19)//'.VALE', 'E', vr=v_vect_2mbr)
    call jelira(vect_2mbr(1:19)//'.VALE', 'LONMAX', nbEqua_2mbr)
    ASSERT(nbEqua .eq. nbEqua_2mbr)
!
! - Access to matrix
    call jeveuo(matr_asse(1:19)//'.REFA', 'L', vk24=refa)
    if (refa(3) .eq. 'ELIML') then
        call mtmchc(matr_asse, 'ELIMF')
    end if
    call jeveuo(matr_asse(1:19)//'.&INT', 'L', jv_matr)
    call dismoi('NB_EQUA', matr_asse, 'MATR_ASSE', repi=nbEqua_matr)
    ASSERT(nbEqua .eq. zi(jv_matr+2))
! 
! - Matrix is symetric ?
    call dismoi('TYPE_MATRICE', matr_asse, 'MATR_ASSE', repk=syme)
!
! - Second member correction for AFFE_CHAR_CINE
!
    call jeveuo(vcine19//'.VALE', 'L', vr=v_vect_cine)
    call mrconl('MULT', jv_matr, 0, 'R', v_vect_2mbr, &
                1)
    call csmbgg(jv_matr, v_vect_2mbr, v_vect_cine, [cbid], [cbid], &
                'R')
!
! - Truncation of second member
!
    if (l_hrom) then
        do iEqua = 1, nbEqua
            if (ds_algorom%v_equa_int(iEqua) .eq. 1) then
                v_vect_2mbr(iEqua) = 0.d0
            end if
        end do
    end if
!
! - Allocate objects
!
    AS_ALLOCATE(vr=v_matr_rom, size=nbMode*nbMode)
    AS_ALLOCATE(vr=v_vect_rom, size=nbMode)
    AS_ALLOCATE(vr=v_mrmult, size=nbEqua)
!
! - Compute reduced objects
!
    if (syme .eq. 'NON_SYM' .or. l_hrom) then
        do iMode = 1, nbMode
            call rsexch(' ', resultName, fieldName, iMode, mode, &
                        iret)
            call jeveuo(mode(1:19)//'.VALE', 'L', vr=v_mode)
            b_n = to_blas_int(nbEqua)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            term1 = ddot(b_n, v_mode, b_incx, v_vect_2mbr, b_incy)
            v_vect_rom(iMode) = term1
            call mrmult('ZERO', jv_matr, v_mode, v_mrmult, 1, &
                        .false._1, l_rom)
            if (l_hrom) then
                do iEqua = 1, nbEqua
                    if (ds_algorom%v_equa_int(iEqua) .eq. 1) then
                        v_mrmult(iEqua) = 0.d0
                    end if
                end do
            end if
            do jMode = 1, nbMode
                call rsexch(' ', resultName, fieldName, jMode, mode, &
                            iret)
                call jeveuo(mode(1:19)//'.VALE', 'L', vr=v_mode)
                b_n = to_blas_int(nbEqua)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                term2 = ddot(b_n, v_mode, b_incx, v_mrmult, b_incy)
                v_matr_rom(nbMode*(iMode-1)+jMode) = term2
            end do
        end do
! - In case of symetric matrix, only half of coefficient of reduced matrix is computed 
    elseif (syme .eq. 'SYMETRI') then
        do iMode = 1, nbMode
            call rsexch(' ', resultName, fieldName, iMode, mode, &
                        iret)
            call jeveuo(mode(1:19)//'.VALE', 'L', vr=v_mode)
            b_n = to_blas_int(nbEqua)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            term1 = ddot(b_n, v_mode, b_incx, v_vect_2mbr, b_incy)
            v_vect_rom(iMode) = term1
            call mrmult('ZERO', jv_matr, v_mode, v_mrmult, 1, &
                        .false._1, l_rom)
            do jMode = iMode, nbMode
                call rsexch(' ', resultName, fieldName, jMode, mode, &
                            iret)
                call jeveuo(mode(1:19)//'.VALE', 'L', vr=v_mode)
                b_n = to_blas_int(nbEqua)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                term2 = ddot(b_n, v_mode, b_incx, v_mrmult, b_incy)
                v_matr_rom(nbMode*(iMode-1)+jMode) = term2
                v_matr_rom(nbMode*(jMode-1)+iMode) = term2
            end do
        end do
    else 
         ASSERT(ASTER_FALSE)
    end if 
!
! - Solve system
!
    call mgauss('NFSP', v_matr_rom, v_vect_rom, nbMode, nbMode, &
                1, det, iret)
    if (l_update_redu) then
        v_gamma = v_gamma+v_vect_rom
    end if
!
! - Project in physical space
!
    call vtzero(vect_solu)
    do iMode = 1, nbMode
        term = v_vect_rom(iMode)
        call rsexch(' ', resultName, fieldName, iMode, mode, &
                    iret)
        call vtaxpy(term, mode, vect_solu)
    end do
    call jeveuo(vect_solu(1:19)//'.VALE', 'E', vr=v_vect_solu)
    call mrconl('MULT', jv_matr, 0, 'R', v_vect_solu, &
                1)
!
! - Clean
!
    AS_DEALLOCATE(vr=v_matr_rom)
    AS_DEALLOCATE(vr=v_vect_rom)
    AS_DEALLOCATE(vr=v_mrmult)
!
end subroutine
