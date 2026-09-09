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

module HHO_utils_module
!
    use HHO_type
    use HHO_size_module
    use HHO_matrix_module
!
    implicit none
!
    private
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/HHO_size_module.h"
#include "asterfort/symt46.h"
#include "MeshTypes_type.h"
!
! --------------------------------------------------------------------------------------------------
!
! HHO - utilitaries
!
! Some common utilitaries
!
! --------------------------------------------------------------------------------------------------
    public :: hhoCopySymPartMat, hhoPrintMat, hhoNorm2Mat, hhoProdSmatVec
    public :: hhoPrintTensor4, hhoPrintTensor4Mangle
    public :: hhoGetTypeFromModel, MatScal2Vec, MatCellScal2Vec
    public :: CellNameL2S, CellNameS2L
    public :: hhoIsIdentityMat
!    private  ::
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoGetTypeFromModel(model, hhoData, ndim)
!
        implicit none
!
        character(len=8), intent(in) :: model
        type(HHO_Data), intent(out) :: hhoData
        integer(kind=8), intent(out) :: ndim
!
! --------------------------------------------------------------------------------------------------
!
! HHO - Utilities
!
! Get type from model
!
! --------------------------------------------------------------------------------------------------
!
! In  model            : name of model
! Out hhoData          : information about the HHO formulation
! Out ndim             : topogical dimension of the HHO formulation
! --------------------------------------------------------------------------------------------------
!
        character(len=8) :: answer
!
        call dismoi('EXI_HHO', model, 'MODELE', repk=answer)
        if (answer .eq. 'OUI') then
            call dismoi('EXI_HHO_QUAR', model, 'MODELE', repk=answer)
            if (answer .eq. 'OUI') then
                call hhoData%initialize(4, 4, 4, 0.d0, ASTER_FALSE, &
                                        ASTER_FALSE)
            else
                call dismoi('EXI_HHO_MQUAR', model, 'MODELE', repk=answer)
                if (answer .eq. 'OUI') then
                    call hhoData%initialize(4, 5, 4, 0.d0, ASTER_FALSE, &
                                            ASTER_FALSE)
                else
                    call dismoi('EXI_HHO_CUBI', model, 'MODELE', repk=answer)
                    if (answer .eq. 'OUI') then
                        call hhoData%initialize(3, 3, 3, 0.d0, ASTER_FALSE, &
                                                ASTER_FALSE)
                    else
                        call dismoi('EXI_HHO_MCUBI', model, 'MODELE', repk=answer)
                        if (answer .eq. 'OUI') then
                            call hhoData%initialize(3, 4, 3, 0.d0, ASTER_FALSE, &
                                                    ASTER_FALSE)
                        else
                            call dismoi('EXI_HHO_QUAD', model, 'MODELE', repk=answer)
                            if (answer .eq. 'OUI') then
                                call hhoData%initialize(2, 2, 2, 0.d0, ASTER_FALSE, &
                                                        ASTER_FALSE)
                            else
                                call dismoi('EXI_HHO_MQUAD', model, 'MODELE', repk=answer)
                                if (answer .eq. 'OUI') then
                                    call hhoData%initialize(2, 3, 2, 0.d0, ASTER_FALSE, &
                                                            ASTER_FALSE)
                                else
                                    call dismoi('EXI_HHO_LINE', model, 'MODELE', repk=answer)
                                    if (answer .eq. 'OUI') then
                                        call hhoData%initialize(1, 1, 1, 0.d0, ASTER_FALSE, &
                                                                ASTER_FALSE)
                                    else
                                        call dismoi('EXI_HHO_MLINE', model, 'MODELE', repk=answer)
                                        if (answer .eq. 'OUI') then
                                            call hhoData%initialize(1, 2, 1, 0.d0, ASTER_FALSE, &
                                                                    ASTER_FALSE)
                                        else
                                            call dismoi('EXI_HHO_CSTE', model, 'MODELE', &
                                                        repk=answer)
                                            if (answer .eq. 'OUI') then
                                                call hhoData%initialize(0, 0, 0, 0.d0, &
                                                                        ASTER_FALSE, ASTER_FALSE)
                                            else
                                                ASSERT(ASTER_FALSE)
                                            end if
                                        end if
                                    end if
                                end if
                            end if
                        end if
                    end if
                end if
            end if
        else
            ASSERT(ASTER_FALSE)
        end if
!
        ndim = 0
        call dismoi('DIM_GEOM', model, 'MODELE', repi=ndim)
        ASSERT(ndim == 2 .or. ndim == 3)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoCopySymPartMat(uplo, mat, size_mat)
!
        implicit none
!
        character(len=1), intent(in) :: uplo
        real(kind=8), dimension(:, :), intent(inout) :: mat
        integer(kind=8), intent(in), optional :: size_mat
!
! --------------------------------------------------------------------------------------------------
!
!   Copy the other part of a symetric matrix
!   In uplo     : 'L' the lower part is given or 'U' for the upper part
!   In size     : size of the matrix
!   Inout mat   : matrix to copy
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: ncols, nrows, i
!
        if (present(size_mat)) then
            nrows = size_mat
            ncols = size_mat
        else
            nrows = size(mat, 1)
            ncols = size(mat, 2)
        end if
!
        ASSERT(nrows == ncols)
!
        if (uplo == 'L') then
            do i = 1, ncols-1
                mat(i, i+1:ncols) = mat(i+1:ncols, i)
            end do
        else if (uplo == 'U') then
            do i = 1, ncols-1
                mat(i+1:ncols, i) = mat(i, i+1:ncols)
            end do
        else
            ASSERT(ASTER_FALSE)
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine CellNameL2S(long, short)
!
        implicit none
!
        integer(kind=8), intent(in) :: long
        character(len=8), intent(out) :: short
!
! --------------------------------------------------------------------------------------------------
!
!   Copy the other part of a symetric matrix
!   In uplo     : 'L' the lower part is given or 'U' for the upper part
!   In size     : size of the matrix
!   Inout mat   : matrix to copy
! --------------------------------------------------------------------------------------------------
!
!
        select case (long)
        case (MT_HEXA8)
            short = "HE8"
        case (MT_TETRA4)
            short = "TE4"
        case (MT_PYRAM5)
            short = "PY5"
        case (MT_PENTA6)
            short = "PE6"
        case (MT_QUAD4)
            short = "QU4"
        case (MT_TRIA3)
            short = "TR3"
        case (MT_SEG2)
            short = "SE2"
        case default
            ASSERT(ASTER_FALSE)
        end select
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine CellNameS2L(short, long)
!
        implicit none
!
        character(len=8), intent(out) :: long
        character(len=8), intent(in) :: short
!
! --------------------------------------------------------------------------------------------------
!
!   Copy the other part of a symetric matrix
!   In uplo     : 'L' the lower part is given or 'U' for the upper part
!   In size     : size of the matrix
!   Inout mat   : matrix to copy
! --------------------------------------------------------------------------------------------------
!
!
        select case (short)
        case ("SE2")
            long = "SEG2"
        case ("SE3")
            long = "SEG3"
        case ("SE4")
            long = "SEG4"
        case ("TR3")
            long = "TRIA3"
        case ("TR6")
            long = "TRIA6"
        case ("TR7")
            long = "TRIA7"
        case ("TR1")
            long = "TRIA10"
        case ("QU4")
            long = "QUAD4"
        case ("QU8")
            long = "QUAD8"
        case ("QU9")
            long = "QUAD9"
        case ("Q12")
            long = "QUAD12"
        case ("HE8")
            long = "HEXA8"
        case ("H20")
            long = "HEXA20"
        case ("H27")
            long = "HEXA27"
        case ("TE4")
            long = "TETRA4"
        case ("T10")
            long = "TETRA10"
        case ("T15")
            long = "TETRA15"
        case ("T20")
            long = "TETRA20"
        case ("PE6")
            long = "PENTA6"
        case ("P15")
            long = "PENTA15"
        case ("P18")
            long = "PENTA18"
        case ("P21")
            long = "PENTA21"
        case ("PY5")
            long = "PYRAM5"
        case ("P13")
            long = "PYRAM13"
        case ("P19")
            long = "PYRAM19"
        case default
            ASSERT(ASTER_FALSE)
        end select
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoPrintMat(mat)
!
        implicit none
!
        real(kind=8), dimension(:, :), intent(in) :: mat
!
! --------------------------------------------------------------------------------------------------
!
!   print matrix
!   In mat   : matrix to print
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: ncols, nrows, i
!
        nrows = size(mat, 1)
        ncols = size(mat, 2)
!
!
        write (6, *) "matrix of", nrows, "rows x ", ncols, "cols"
        do i = 1, nrows
            write (6, '(50ES13.6)') mat(i, 1:ncols)
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoPrintTensor4(tens)
!
        implicit none
!
        real(kind=8), dimension(3, 3, 3, 3), intent(in) :: tens
!
! --------------------------------------------------------------------------------------------------
!
!   print fourth order tensor
!   In tens   : tensor to print
! --------------------------------------------------------------------------------------------------
!
!
        integer(kind=8) :: i, k
!
!
        write (6, *) "tensor of 9 rows x 9 cols"
        do i = 1, 3
            do k = 1, 3
                write (6, '(50F14.7)') tens(i, 1, k, 1), tens(i, 1, k, 2), tens(i, 1, k, 3), &
                    tens(i, 2, k, 1), tens(i, 2, k, 2), tens(i, 2, k, 3), &
                    tens(i, 3, k, 1), tens(i, 3, k, 2), tens(i, 3, k, 3)
            end do
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoPrintTensor4Mangle(tens)
!
        implicit none
!
        real(kind=8), dimension(3, 3, 3, 3), intent(in) :: tens
!
! --------------------------------------------------------------------------------------------------
!
!   print fourth order tensor
!   In tens   : tensor to print
! --------------------------------------------------------------------------------------------------
!
!
        real(kind=8) :: tens6(6, 6)
!
        tens6 = 0.d0
        call symt46(tens, tens6)
        call hhoPrintMat(tens6)
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    function hhoNorm2Mat(mat) result(norm)
!
        implicit none
!
        real(kind=8), dimension(:, :), intent(in) :: mat
        real(kind=8) :: norm
!
! --------------------------------------------------------------------------------------------------
!
!   compute the frobenius-norm of a matrix
!   In mat   : matrix to evaluate
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: nrows, ncols, irow, jcol
!
        nrows = size(mat, 1)
        ncols = size(mat, 2)
!
        norm = 0.d0
        do jcol = 1, ncols
            do irow = 1, nrows
                norm = norm+mat(irow, jcol)**2
            end do
        end do
!
        norm = sqrt(norm)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine hhoProdSmatVec(SymMat, vect, ndim, resu)
!
        implicit none
!
        real(kind=8), dimension(6), intent(in) :: SymMat
        real(kind=8), dimension(3), intent(in) :: vect
        integer(kind=8), intent(in) :: ndim
        real(kind=8), dimension(3), intent(out) :: resu
!
! --------------------------------------------------------------------------------------------------
!
!   evaluate the product of a symetrix matrix with a vector
!   In SymMat       : symmetric matrix (M11, M22, M33, Rac2*M12, Rac2*M13, Rac2*M23)
!   In vector       : vector
!   In ndim         : dim of the matrix
!   Out resu        : result of the product
! --------------------------------------------------------------------------------------------------
!
        real(kind=8), parameter :: un_rac2 = 1.d0/sqrt(2.d0)
!
        resu = 0.d0
!
        if (ndim == 2) then
            resu(1) = SymMat(1)*vect(1)+SymMat(4)*un_rac2*vect(2)
            resu(2) = SymMat(2)*vect(2)+SymMat(4)*un_rac2*vect(1)
        else if (ndim == 3) then
            resu(1) = SymMat(1)*vect(1)+(SymMat(4)*vect(2)+SymMat(5)*vect(3))*un_rac2
            resu(2) = SymMat(2)*vect(2)+(SymMat(4)*vect(1)+SymMat(6)*vect(3))*un_rac2
            resu(3) = SymMat(3)*vect(3)+(SymMat(5)*vect(1)+SymMat(6)*vect(2))*un_rac2
        else
            ASSERT(ASTER_FALSE)
        end if
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine MatScal2Vec(hhoCell, hhoData, mat_scal, mat_vec)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Data), intent(in) :: hhoData
        type(HHO_matrix), intent(in) :: mat_scal
        type(HHO_matrix), intent(out) :: mat_vec
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Create the vectorial matrix by tensorization of the scalar matrix
!   In hhoCell      : the current HHO Cell
!   In hhoData      : information on HHO methods
!   In mat_scal     : matrix for the scalar problem
!   Out mat_vec     : matrix for the vectorial problem
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: cbs_comp, fbs_comp, cbs, fbs, jbeginVec, jendVec
        integer(kind=8) :: total_dofs, idim, jbeginCell, jendCell, jbeginFace, jendFace, iFace
        integer(kind=8) :: jFace, ibeginFace, iendFace, ibeginVec, iendVec
        integer(kind=8) :: faces_dofs, faces_dofs_comp, total_dofs_comp
! --------------------------------------------------------------------------------------------------
!
! -- number of dofs
        call hhoMecaDofs(hhoCell, hhoData, cbs, fbs, total_dofs)
        faces_dofs = total_dofs-cbs
!
        cbs_comp = cbs/hhoCell%ndim
        fbs_comp = fbs/hhoCell%ndim
        faces_dofs_comp = faces_dofs/hhoCell%ndim
        total_dofs_comp = total_dofs/hhoCell%ndim
!
! -- copy the scalar matrix in the vectorial matrix
        call mat_vec%initialize(hhoCell%ndim*mat_scal%nrows, hhoCell%ndim*mat_scal%ncols, 0.d0)
!
        do idim = 1, hhoCell%ndim
! --------- copy volumetric part
            jbeginCell = faces_dofs+(idim-1)*cbs_comp+1
            jendCell = jbeginCell+cbs_comp-1
!
            mat_vec%m(jbeginCell:jendCell, jbeginCell:jendCell) &
                = mat_scal%m(faces_dofs_comp+1:total_dofs_comp, faces_dofs_comp+1:total_dofs_comp)
!
! --------- copy faces part
            do iFace = 1, hhoCell%nbfaces
                ibeginFace = (iFace-1)*fbs+(idim-1)*fbs_comp+1
                iendFace = ibeginFace+fbs_comp-1
                ibeginVec = (iFace-1)*fbs_comp+1
                iendVec = ibeginVec+fbs_comp-1
!
                do jFace = 1, hhoCell%nbfaces
                    jbeginFace = (jFace-1)*fbs+(idim-1)*fbs_comp+1
                    jendFace = jbeginFace+fbs_comp-1
                    jbeginVec = (jFace-1)*fbs_comp+1
                    jendVec = jbeginVec+fbs_comp-1
!
                    mat_vec%m(ibeginFace:iendFace, jbeginFace:jendFace) &
                        = mat_scal%m(ibeginVec:iendVec, jbeginVec:jendVec)
                end do
!
! --------- copy coupled part
                mat_vec%m(jbeginCell:jendCell, ibeginFace:iendFace) &
                    = mat_scal%m(faces_dofs_comp+1:total_dofs_comp, ibeginVec:iendVec)
                mat_vec%m(ibeginFace:iendFace, jbeginCell:jendCell) &
                    = mat_scal%m(ibeginVec:iendVec, faces_dofs_comp+1:total_dofs_comp)
            end do
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    subroutine MatCellScal2Vec(hhoCell, hhoData, mat_scal, mat_vec)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: hhoCell
        type(HHO_Data), intent(in) :: hhoData
        real(kind=8), intent(in) :: mat_scal(MSIZE_CELL_SCAL, MSIZE_CELL_SCAL)
        type(HHO_matrix), intent(out) :: mat_vec
!
! --------------------------------------------------------------------------------------------------
!   HHO
!
!   Create the vectorial matrix by tensorization of the scalar matrix
!   In hhoCell      : the current HHO Cell
!   In hhoData      : information on HHO methods
!   In mat_scal     : matrix for the scalar problem
!   Out mat_vec     : matrix for the vectorial problem
!
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: cbs_comp, cbs, fbs
        integer(kind=8) :: total_dofs, idim, jbeginCell, jendCell
! --------------------------------------------------------------------------------------------------
!
! -- number of dofs
        call hhoMecaDofs(hhoCell, hhoData, cbs, fbs, total_dofs)
!
        cbs_comp = cbs/hhoCell%ndim
!
! -- copy the scalar matrix in the vectorial matrix
        call mat_vec%initialize(cbs, cbs, 0.d0)
!
        do idim = 1, hhoCell%ndim
! --------- copy volumetric part
            jbeginCell = (idim-1)*cbs_comp+1
            jendCell = jbeginCell+cbs_comp-1
!
            mat_vec%m(jbeginCell:jendCell, jbeginCell:jendCell) = mat_scal(1:cbs_comp, 1:cbs_comp)
        end do
!
    end subroutine
!
!===================================================================================================
!
!===================================================================================================
!
    function hhoIsIdentityMat(mat, size) result(id)
!
        implicit none
!
        real(kind=8), dimension(:, :), intent(in) :: mat
        integer(kind=8), intent(in) :: size
        aster_logical :: id
!
! --------------------------------------------------------------------------------------------------
!
!   Return True if the matrix is Identity
!   In mat   : matrix to evaluate
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: i, j
!
        id = ASTER_TRUE
        do j = 1, size
            if (abs(mat(j, j)-1.d0) .ge. 1.d-10) then
                id = ASTER_FALSE
                exit
            end if
            do i = 1, j-1
                if (abs(mat(i, j)) .ge. 1.d-10) then
                    id = ASTER_FALSE
                    exit
                end if
            end do
        end do
!
    end function
!
end module
