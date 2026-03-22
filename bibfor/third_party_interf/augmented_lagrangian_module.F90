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

module augmented_lagrangian_module
!
#include "asterf_types.h"
#include "asterf_petsc.h"
!
!
    use aster_petsc_module
    use petsc_data_module, only: nomat_courant
    use saddle_point_context_type
    use augmented_lagrangian_context_type
    use saddle_point_data_module, only: get_saddle_point_context, set_augmented_lagrangian_context, &
                                        get_augmented_lagrangian_context, clear_augmented_lagrangian_context
!
    implicit none
!
    private
#include "asterf.h"
#include "asterfort/assert.h"
!
    public :: augmented_lagrangian_setup, augmented_lagrangian_apply, augmented_lagrangian_destroy
!
#ifdef ASTER_HAVE_PETSC
contains
!
    subroutine augmented_lagrangian_setup(al_pc, ierr)
        PC, intent(inout)           :: al_pc
        PetscErrorCode, intent(out) :: ierr
        type(saddlepoint_ctxt), pointer :: sp_ctxt => null()
        !
        sp_ctxt => get_saddle_point_context(nomat_courant)
        ASSERT(associated(sp_ctxt))
        call set_augmented_lagrangian_context(nomat_courant, &
                                              new_augmented_lagrangian_context(sp_ctxt))
        !
        ierr = 0
        !
    end subroutine augmented_lagrangian_setup
!
! PCApply operation used for the PCShell based on the Augmented Lagrangian Method
!
    subroutine augmented_lagrangian_apply(al_pc, x, y, ierr)
        ! Dummy arguments
        !
        PC, intent(in)    :: al_pc
        Vec, intent(in)   :: x
        Vec, intent(inout):: y
        PetscErrorCode, intent(out) :: ierr
        !
        ! Local Variables
        PetscReal :: beta
        type(augm_lagr_ctxt), pointer :: sp_pc_ctxt => null()
        type(saddlepoint_ctxt), pointer :: sp_ctxt => null()
        !
        ! Use work space and scatter context defined in the underlying saddle_point_context
        sp_pc_ctxt => get_augmented_lagrangian_context(nomat_courant)
        ASSERT(associated(sp_pc_ctxt))
        sp_ctxt => sp_pc_ctxt%sp_ctxt
        ASSERT(associated(sp_ctxt))
        !
        ! Load values of  x1, x2 from x : use SCATTER_REVERSE
        !
        call VecScatterBegin(sp_ctxt%scatter_to_phys, x, sp_ctxt%x1, &
       &  INSERT_VALUES, SCATTER_REVERSE, ierr)
        ASSERT(ierr == 0)
        call VecScatterEnd(sp_ctxt%scatter_to_phys, x, sp_ctxt%x1, &
       & INSERT_VALUES, SCATTER_REVERSE, ierr)
        ASSERT(ierr == 0)
        !
        call VecScatterBegin(sp_ctxt%scatter_to_lag1, x, sp_ctxt%x2, &
       & INSERT_VALUES, SCATTER_REVERSE, ierr)
        ASSERT(ierr == 0)
        call VecScatterEnd(sp_ctxt%scatter_to_lag1, x, sp_ctxt%x2, &
       & INSERT_VALUES, SCATTER_REVERSE, ierr)
        ASSERT(ierr == 0)
        !
        ! Apply Augmented Lgrangian preconditioner :
        ! Thèse Sylvain p. 113
        ! (y1) = ( LL^T  2 C^T     )^{-1} (x1) = ( L^{-T} L^{-1}  0      ) ( I  2*gamma*C^T ) (x1)
        ! (y2)   ( 0    -1/gamma*I )      (x2)   ( 0             -gamma*I) ( 0  I           ) (x2)
        !
        ! c'est à dire : y1 = L^{-T} L^{-1} (x1 + 2*gamma*C^T*x2)
        !                y2 = -gamma*x2
        !
        ! y2 = -gamma*x2
        call VecCopy(sp_ctxt%x2, sp_ctxt%y2, ierr)
        ASSERT(ierr == 0)
        beta = -sp_pc_ctxt%gamma
        call VecScale(sp_ctxt%y2, beta, ierr)
        ASSERT(ierr == 0)
        !
        ! x_1 = x1 + 2*gamma*C^T*x2
        !
        ! xtmp = C^T x2
        call MatMultTranspose(sp_ctxt%c_mat, sp_ctxt%x2, sp_ctxt%xtmp, ierr)
        ASSERT(ierr == 0)
        ! x1 = x1 + (2*gamma)* xtmp
        beta = 2.0d0*sp_pc_ctxt%gamma
        call VecAXPY(sp_ctxt%x1, beta, sp_ctxt%xtmp, ierr)
        ASSERT(ierr == 0)
        !
        ! Apply pcphy
        call PCApply(sp_pc_ctxt%pcphy, sp_ctxt%x1, sp_ctxt%y1, ierr)
        ASSERT(ierr == 0)
        !
        ! Load result to y : this time, use SCATTER_FORWARD
        !
        call VecScatterBegin(sp_ctxt%scatter_to_phys, sp_ctxt%y1, y, &
          &  INSERT_VALUES, SCATTER_FORWARD, ierr)
        ASSERT(ierr == 0)
        call VecScatterEnd(sp_ctxt%scatter_to_phys, sp_ctxt%y1, y, &
          &  INSERT_VALUES, SCATTER_FORWARD, ierr)
        ASSERT(ierr == 0)
        !
        call VecScatterBegin(sp_ctxt%scatter_to_lag1, sp_ctxt%y2, y, &
          &  INSERT_VALUES, SCATTER_FORWARD, ierr)
        ASSERT(ierr == 0)
        call VecScatterEnd(sp_ctxt%scatter_to_lag1, sp_ctxt%y2, y, &
          &  INSERT_VALUES, SCATTER_FORWARD, ierr)
        ASSERT(ierr == 0)
        !
    end subroutine augmented_lagrangian_apply
!
    subroutine augmented_lagrangian_destroy(al_pc, ierr)
        ! Dummy arguments
        !
        PC, intent(inout)    :: al_pc
        PetscErrorCode, intent(out) :: ierr
        !
        call clear_augmented_lagrangian_context(nomat_courant)
        !
        ierr = 0
        !
    end subroutine augmented_lagrangian_destroy
!
#else
contains
!
    subroutine augmented_lagrangian_setup(al_pc, ierr)
        integer(kind=8), intent(inout):: al_pc
        integer(kind=8), intent(out)  :: ierr
        integer(kind=8) :: idummy
        idummy = al_pc
        ierr = 1
    end subroutine augmented_lagrangian_setup
!
    subroutine augmented_lagrangian_apply(al_pc, x, y, ierr)
        integer(kind=8), intent(in)    :: al_pc
        integer(kind=8), intent(in)   :: x
        integer(kind=8), intent(inout):: y
        integer(kind=8), intent(out) :: ierr
        integer(kind=8) :: idummy
        idummy = al_pc+x+y
        ierr = 1
    end subroutine augmented_lagrangian_apply
!
    subroutine augmented_lagrangian_destroy(al_pc, ierr)
        integer(kind=8), intent(inout)    :: al_pc
        integer(kind=8), intent(out) :: ierr
        integer(kind=8) :: idummy
        idummy = al_pc
        ierr = 1
    end subroutine augmented_lagrangian_destroy
#endif
end module augmented_lagrangian_module
