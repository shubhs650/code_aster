! --------------------------------------------------------------------
! Copyright (C) 2016 - 2026 - EDF - www.code-aster.org
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

module saddle_point_data_module
!
!
! aslint: disable=W1403
!
    use saddle_point_context_type
    use augmented_lagrangian_context_type
!
    implicit none
!
    private
#include "asterfort/assert.h"

    integer(kind=8), parameter :: max_sp_ctxt = 5

    character(len=19), save :: sp_context_name(max_sp_ctxt) = " "
    character(len=19), save :: sp_pc_context_name(max_sp_ctxt) = " "

    type(saddlepoint_ctxt), target, save :: sp_context(max_sp_ctxt)
    type(augm_lagr_ctxt), target, save :: sp_pc_context(max_sp_ctxt)

    public :: set_saddle_point_context, get_saddle_point_context, clear_saddle_point_context
    public :: set_augmented_lagrangian_context, get_augmented_lagrangian_context
    public :: clear_augmented_lagrangian_context

contains

    subroutine set_saddle_point_context(matas, ctxt)
        character(len=*), intent(in) :: matas
        type(saddlepoint_ctxt), intent(in) :: ctxt

        integer(kind=8) :: idx
        character(len=19) :: matas_key

        matas_key = normalize_name(matas)
        ASSERT(len_trim(matas_key) > 0)
        idx = find_slot(sp_context_name, matas_key)
        if (idx == 0) idx = find_free_slot(sp_context_name)
        ASSERT(idx > 0)

        if (sp_context_name(idx) /= " ") then
            call clear_augmented_lagrangian_context(matas_key)
            call free_saddle_point_context(sp_context(idx))
        end if

        sp_context_name(idx) = matas_key
        sp_context(idx) = ctxt
    end subroutine set_saddle_point_context

    function get_saddle_point_context(matas) result(ctxt)
        character(len=*), intent(in) :: matas
        type(saddlepoint_ctxt), pointer :: ctxt

        integer(kind=8) :: idx
        character(len=19) :: matas_key

        nullify (ctxt)
        matas_key = normalize_name(matas)
        if (len_trim(matas_key) == 0) return
        idx = find_slot(sp_context_name, matas_key)
        if (idx /= 0) ctxt => sp_context(idx)
    end function get_saddle_point_context

    subroutine clear_saddle_point_context(matas)
        character(len=*), intent(in) :: matas

        integer(kind=8) :: idx
        character(len=19) :: matas_key

        matas_key = normalize_name(matas)
        if (len_trim(matas_key) == 0) return
        idx = find_slot(sp_context_name, matas_key)
        if (idx == 0) return

        call clear_augmented_lagrangian_context(matas_key)
        call free_saddle_point_context(sp_context(idx))
        sp_context_name(idx) = " "
    end subroutine clear_saddle_point_context

    subroutine set_augmented_lagrangian_context(matas, ctxt)
        character(len=*), intent(in) :: matas
        type(augm_lagr_ctxt), intent(in) :: ctxt

        integer(kind=8) :: idx
        character(len=19) :: matas_key

        matas_key = normalize_name(matas)
        ASSERT(len_trim(matas_key) > 0)
        idx = find_slot(sp_pc_context_name, matas_key)
        if (idx == 0) idx = find_free_slot(sp_pc_context_name)
        ASSERT(idx > 0)

        if (sp_pc_context_name(idx) /= " ") then
            call free_augm_lagrangian_context(sp_pc_context(idx))
        end if

        sp_pc_context_name(idx) = matas_key
        sp_pc_context(idx) = ctxt
    end subroutine set_augmented_lagrangian_context

    function get_augmented_lagrangian_context(matas) result(ctxt)
        character(len=*), intent(in) :: matas
        type(augm_lagr_ctxt), pointer :: ctxt

        integer(kind=8) :: idx
        character(len=19) :: matas_key

        nullify (ctxt)
        matas_key = normalize_name(matas)
        if (len_trim(matas_key) == 0) return
        idx = find_slot(sp_pc_context_name, matas_key)
        if (idx /= 0) ctxt => sp_pc_context(idx)
    end function get_augmented_lagrangian_context

    subroutine clear_augmented_lagrangian_context(matas)
        character(len=*), intent(in) :: matas

        integer(kind=8) :: idx
        character(len=19) :: matas_key

        matas_key = normalize_name(matas)
        if (len_trim(matas_key) == 0) return
        idx = find_slot(sp_pc_context_name, matas_key)
        if (idx == 0) return

        call free_augm_lagrangian_context(sp_pc_context(idx))
        sp_pc_context_name(idx) = " "
    end subroutine clear_augmented_lagrangian_context

    function normalize_name(matas) result(matas_key)
        character(len=*), intent(in) :: matas
        character(len=19) :: matas_key

        matas_key = " "
        matas_key = trim(matas)
    end function normalize_name

    function find_slot(name_list, matas) result(idx)
        character(len=19), intent(in) :: name_list(max_sp_ctxt)
        character(len=19), intent(in) :: matas
        integer(kind=8) :: idx
        integer(kind=8) :: ii

        idx = 0
        do ii = 1, max_sp_ctxt
            if (name_list(ii) == matas) then
                idx = ii
                exit
            end if
        end do
    end function find_slot

    function find_free_slot(name_list) result(idx)
        character(len=19), intent(in) :: name_list(max_sp_ctxt)
        integer(kind=8) :: idx
        integer(kind=8) :: ii

        idx = 0
        do ii = 1, max_sp_ctxt
            if (name_list(ii) == " ") then
                idx = ii
                exit
            end if
        end do
    end function find_free_slot
!
end module saddle_point_data_module
