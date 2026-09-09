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
module HHO_measure_module
!
    use HHO_type
    use compensated_ops_module, only: sum, dot_product, matmul, norm2
!
    implicit none
!
    private
#include "asterc/r8maem.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "MeshTypes_type.h"
!
! --------------------------------------------------------------------------------------------------
!
! HHO
!
! Compute measure for elements:
!  - volume for 3D elements
!  - surface for 2D elements
!  - length for 1D element
!
! --------------------------------------------------------------------------------------------------
!
    public   :: hhoMeasureCell, hhoMeasureFace, hhoDiameterCell, hhoDiameterFace
    public   :: hhoLengthBoundingBoxCell, hhoLengthBoundingBoxFace
    public   :: hhoCenterBoundingBoxCell, hhoCenterBoundingBoxFace
    public   :: hho_vol_tetra, hho_surface_tri
    private  :: hho_vol_hexa, hho_surface_quad, hho_length_edge
    private  :: hho_vol_prism, hho_vol_pyram
    private  :: hhoDiameter, prod_vec
    private  :: hhoBoundingBoxCell, hhoBoundingBoxFace
!
contains
!
!===================================================================================================
!
!===================================================================================================
!
    function prod_vec(v0, v1) result(v2)
!
        implicit none
!
        real(kind=8), dimension(3), intent(in) :: v0, v1
        real(kind=8), dimension(3) :: v2
!
! --------------------------------------------------------------------------------------------------
!  In v0        :: vector 0
!  In v1        :: vector 1
!  Out v2       :: vector 2
! --------------------------------------------------------------------------------------------------
!
        v2 = 0.d0
!
        v2(1) = v0(2)*v1(3)-v0(3)*v1(2)
        v2(2) = -(v0(1)*v1(3)-v0(3)*v1(1))
        v2(3) = v0(1)*v1(2)-v0(2)*v1(1)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function hho_vol_hexa(nodes) result(vol)
!
        implicit none
!
        real(kind=8), dimension(3, 8), intent(in) :: nodes
        real(kind=8) :: vol
!
! --------------------------------------------------------------------------------------------------
!  In nodes        :: list of nodes
!  Out vol         :: volume
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8), dimension(4, 5) :: tets
        integer(kind=8) :: i, j
        real(kind=8) :: nodestet(3, 4), vol_tetra(5)
! --------------------------------------------------------------------------------------------------
!
! --- split the hexa in 5 tets - see hhoSplitSimplex
        tets(1:4, 1) = (/1, 2, 4, 5/)
        tets(1:4, 2) = (/2, 3, 4, 7/)
        tets(1:4, 3) = (/2, 4, 5, 7/)
        tets(1:4, 4) = (/2, 5, 6, 7/)
        tets(1:4, 5) = (/4, 5, 7, 8/)
!
        vol = 0.d0
        do i = 1, 5
            do j = 1, 4
                nodestet(1:3, j) = nodes(1:3, tets(j, i))
            end do
            vol_tetra(i) = hho_vol_tetra(nodestet)
        end do
        vol = sum(vol_tetra)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function hho_vol_prism(nodes) result(vol)
!
        implicit none
!
        real(kind=8), dimension(3, 6), intent(in) :: nodes
        real(kind=8) :: vol, vol_tetra(3)
!
! --------------------------------------------------------------------------------------------------
!  In nodes        :: list of nodes
!  Out vol         :: volume
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8), dimension(4, 3) :: tets
        integer(kind=8) :: i, j
        real(kind=8) :: nodestet(3, 4)
! --------------------------------------------------------------------------------------------------
!
! --- split the prsim in 3 tets - see hhoSplitSimplex
        tets(1:4, 1) = (/1, 2, 3, 4/)
        tets(1:4, 2) = (/2, 3, 4, 5/)
        tets(1:4, 3) = (/3, 4, 5, 6/)
!
        vol = 0.d0
        do i = 1, 3
            do j = 1, 4
                nodestet(1:3, j) = nodes(1:3, tets(j, i))
            end do
            vol_tetra(i) = hho_vol_tetra(nodestet)
        end do
        vol = sum(vol_tetra)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function hho_vol_pyram(nodes) result(vol)
!
        implicit none
!
        real(kind=8), dimension(3, 5), intent(in) :: nodes
        real(kind=8) :: vol
!
! --------------------------------------------------------------------------------------------------
!  In nodes        :: list of nodes
!  Out vol         :: volume
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8), dimension(4, 2) :: tets
        integer(kind=8) :: i, j
        real(kind=8) :: nodestet(3, 4), vol_tetra(2)
! --------------------------------------------------------------------------------------------------
!
! --- split the pyramid in 2 tets - see hhoSplitSimplex
        tets(1:4, 1) = (/1, 2, 3, 5/)
        tets(1:4, 2) = (/1, 3, 4, 5/)
!
        vol = 0.d0
        do i = 1, 2
            do j = 1, 4
                nodestet(1:3, j) = nodes(1:3, tets(j, i))
            end do
            vol_tetra(i) = hho_vol_tetra(nodestet)
        end do
        vol = sum(vol_tetra)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function hho_vol_tetra(nodes) result(vol)
!
        implicit none
!
        real(kind=8), dimension(3, 4), intent(in) :: nodes
        real(kind=8) :: vol
!
! --------------------------------------------------------------------------------------------------
!  In nodes        :: list of nodes
!  Out vol         :: volume
! --------------------------------------------------------------------------------------------------
!
        real(kind=8), dimension(3) :: v0, v1, v2, cross
! --------------------------------------------------------------------------------------------------
!
        vol = 0.d0
        v0(1:3) = nodes(1:3, 2)-nodes(1:3, 1)
        v1(1:3) = nodes(1:3, 3)-nodes(1:3, 1)
        v2(1:3) = nodes(1:3, 4)-nodes(1:3, 1)
!
        cross(1:3) = prod_vec(v1, v2)
        vol = abs(dot_product(v0, cross)/6.d0)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function hho_surface_quad(nodes) result(surface)
!
        implicit none
!
        real(kind=8), dimension(3, 4), intent(in) :: nodes
        real(kind=8) :: surface
!
! --------------------------------------------------------------------------------------------------
!  In nodes        :: list of nodes
!  Out vol         :: surface
! --------------------------------------------------------------------------------------------------
!
        real(kind=8), dimension(3) :: e1, e2, e3, e4, d1, d2
        real(kind=8) :: l1, l2, l3, l4, ld1, ld2
!
! ---- edge
        e1(1:3) = nodes(1:3, 2)-nodes(1:3, 1)
        e2(1:3) = nodes(1:3, 3)-nodes(1:3, 2)
        e3(1:3) = nodes(1:3, 4)-nodes(1:3, 3)
        e4(1:3) = nodes(1:3, 1)-nodes(1:3, 4)
! ---- diagonals
        d1(1:3) = nodes(1:3, 3)-nodes(1:3, 1)
        d2(1:3) = nodes(1:3, 4)-nodes(1:3, 2)
!
! ---- lengths
        l1 = dot_product(e1, e1)
        l2 = dot_product(e2, e2)
        l3 = dot_product(e3, e3)
        l4 = dot_product(e4, e4)
        ld1 = dot_product(d1, d1)
        ld2 = dot_product(d2, d2)
!
        surface = sqrt(4.d0*ld1*ld2-(l1-l2+l3-l4)**2)/4.d0
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function hho_surface_tri(nodes) result(surface)
!
        implicit none
!
        real(kind=8), dimension(3, 3), intent(in) :: nodes
        real(kind=8) :: surface
!
!---------------------------------------------------------------------------------------------------
!  In nodes        :: list of nodes
!  Out vol         :: surface
!---------------------------------------------------------------------------------------------------
!
        real(kind=8) :: a, b, c
        real(kind=8) :: x, y, z, t, s
!
        a = norm2(nodes(1:3, 2)-nodes(1:3, 1))
        b = norm2(nodes(1:3, 3)-nodes(1:3, 1))
        c = norm2(nodes(1:3, 3)-nodes(1:3, 2))
!
        if (a < b) then
            t = a
            a = b
            b = t
        end if
        if (b < c) then
            t = b
            b = c
            c = t
        end if
        if (a < b) then
            t = a
            a = b
            b = t
        end if
!
        x = a+(b+c)
        y = c-(a-b)
        z = c+(a-b)
        t = a+(b-c)
        s = x*y*z*t
!
        surface = sqrt(s)/4.d0
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function hho_length_edge(nodes) result(length)
!
        implicit none
!
        real(kind=8), dimension(3, 2), intent(in) :: nodes
        real(kind=8) :: length
!
! --------------------------------------------------------------------------------------------------
!  In nodes        :: list of nodes
!  Out vol         :: length
! --------------------------------------------------------------------------------------------------
!
        real(kind=8), dimension(3) :: v0
! --------------------------------------------------------------------------------------------------
!
        v0(1:3) = nodes(1:3, 2)-nodes(1:3, 1)
!
        length = norm2(v0)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function hhoMeasureCell(cell) result(measure)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: cell
        real(kind=8) :: measure
!
! --------------------------------------------------------------------------------------------------
!  In HHO_Cell           :: cell HHO
!  Out measure           :: measure of the cell
! --------------------------------------------------------------------------------------------------
!
!
        measure = 0.d0
!
        if (cell%typema == MT_HEXA8) then
            measure = hho_vol_hexa(cell%coorno(1:3, 1:8))
        else if (cell%typema == MT_TETRA4) then
            measure = hho_vol_tetra(cell%coorno(1:3, 1:4))
        else if (cell%typema == MT_PYRAM5) then
            measure = hho_vol_pyram(cell%coorno(1:3, 1:5))
        else if (cell%typema == MT_PENTA6) then
            measure = hho_vol_prism(cell%coorno(1:3, 1:6))
        else if (cell%typema == MT_QUAD4) then
            measure = hho_surface_quad(cell%coorno(1:3, 1:4))
        else if (cell%typema == MT_TRIA3) then
            measure = hho_surface_tri(cell%coorno(1:3, 1:3))
        else
            ASSERT(ASTER_FALSE)
        end if
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function hhoMeasureFace(face) result(measure)
!
        implicit none
!
        type(HHO_Face), intent(in) :: face
        real(kind=8) :: measure
!
! --------------------------------------------------------------------------------------------------
!  In HHO_Face           :: face HHO
!  Out measure           :: measure of the face
! --------------------------------------------------------------------------------------------------
!
        measure = 0.d0
!
        if (face%typema == MT_QUAD4) then
            measure = hho_surface_quad(face%coorno(1:3, 1:4))
        else if (face%typema == MT_TRIA3) then
            measure = hho_surface_tri(face%coorno(1:3, 1:3))
        else if (face%typema == MT_SEG2) then
            measure = hho_length_edge(face%coorno(1:3, 1:2))
        else
            ASSERT(ASTER_FALSE)
        end if
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function hhoDiameter(coorno, nbnodes) result(diam)
!
        implicit none
        integer(kind=8), intent(in) :: nbnodes
        real(kind=8), dimension(3, nbnodes), intent(in) :: coorno
        real(kind=8) :: diam
!
! --------------------------------------------------------------------------------------------------
!  In coorno            :: coordinates of the nodes
!  In nbnodes           :: number of nodes
!  Out diam              :: maximum length of two distinc nodes
! --------------------------------------------------------------------------------------------------
!
        integer(kind=8) :: inode, jnode
        real(kind=8), dimension(3) :: vector
!
        diam = 0.d0
!
        do inode = 1, nbnodes
            do jnode = inode+1, nbnodes
                vector(1:3) = coorno(1:3, inode)-coorno(1:3, jnode)
                diam = max(norm2(vector), diam)
            end do
        end do
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function hhoDiameterCell(cell) result(measure)
!
        implicit none
!
        type(HHO_Cell), intent(in) :: cell
        real(kind=8) :: measure
!
! --------------------------------------------------------------------------------------------------
!  In HHO_Cell              :: cell HHO
!  Out measure              :: maximum length of two distinc nodes
! --------------------------------------------------------------------------------------------------
!
!
        measure = 0.d0
!
        if (cell%typema == MT_HEXA8) then
            measure = hhoDiameter(cell%coorno(1:3, 1:8), 8)
        else if (cell%typema == MT_TETRA4) then
            measure = hhoDiameter(cell%coorno(1:3, 1:4), 4)
        else if (cell%typema == MT_PYRAM5) then
            measure = hhoDiameter(cell%coorno(1:3, 1:5), 5)
        else if (cell%typema == MT_PENTA6) then
            measure = hhoDiameter(cell%coorno(1:3, 1:6), 6)
        else if (cell%typema == MT_QUAD4) then
            measure = hhoDiameter(cell%coorno(1:3, 1:4), 4)
        else if (cell%typema == MT_TRIA3) then
            measure = hhoDiameter(cell%coorno(1:3, 1:3), 3)
        else
            ASSERT(ASTER_FALSE)
        end if
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function hhoDiameterFace(face) result(measure)
!
        implicit none
!
        type(HHO_Face), intent(in) :: face
        real(kind=8) :: measure
!
! --------------------------------------------------------------------------------------------------
!  In HHO_Face              :: face HHO
!  Out measure              :: maximum length of two distinc nodes
! --------------------------------------------------------------------------------------------------
!
        measure = 0.d0
        if (face%typema == MT_QUAD4) then
            measure = hhoDiameter(face%coorno(1:3, 1:4), 4)
        else if (face%typema == MT_TRIA3) then
            measure = hhoDiameter(face%coorno(1:3, 1:3), 3)
        else if (face%typema == MT_SEG2) then
            measure = hhoDiameter(face%coorno(1:3, 1:2), 2)
        else
            ASSERT(ASTER_FALSE)
        end if
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function hhoBoundingBoxCell(hhoCell, axes) result(box)
!
        implicit none
!
        type(HHO_Cell), intent(in)    :: hhoCell
        real(kind=8), intent(in)      :: axes(3, 3)
        real(kind=8), dimension(3, 2) :: box
!
! --------------------------------------------------------------------------------------------------
!  In HHO_Cell              :: cell HHO
!  In axes                  :: axes to use
!  Out box                  :: bounding box of the cell
! --------------------------------------------------------------------------------------------------
        real(kind=8), dimension(3) :: pt
        real(kind=8) :: rotmat(3, 3)
        integer(kind=8) :: inode, idim
! --------------------------------------------------------------------------------------------------
!
        box(:, 1) = r8maem()
        box(:, 2) = -r8maem()
        pt = 0.d0
        rotmat = transpose(axes)
!
        do inode = 1, hhoCell%nbnodes
            pt = matmul(rotmat, hhoCell%coorno(1:3, inode))
            do idim = 1, 3
                box(idim, 1) = min(box(idim, 1), pt(idim))
                box(idim, 2) = max(box(idim, 2), pt(idim))
            end do
        end do
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function hhoBoundingBoxFace(hhoFace, axes) result(box)
!
        implicit none
!
        type(HHO_Face), intent(in)    :: hhoFace
        real(kind=8), intent(in)      :: axes(3, 2)
        real(kind=8), dimension(2, 2)  :: box
!
! --------------------------------------------------------------------------------------------------
!  In HHO_Face              :: face HHO
!  In axes                  :: axes to use
!  Out box                  :: bounding box of the face
! --------------------------------------------------------------------------------------------------
        real(kind=8) :: pt(2)
        real(kind=8) :: rotmat(2, 3)
        integer(kind=8) :: inode, idim
! --------------------------------------------------------------------------------------------------
!
        rotmat = transpose(axes)
        box(:, 1) = r8maem()
        box(:, 2) = -r8maem()
        pt = 0.d0
!
        do inode = 1, hhoFace%nbnodes
            pt = matmul(rotmat, hhoFace%coorno(1:3, inode))
            do idim = 1, 2
                box(idim, 1) = min(box(idim, 1), pt(idim))
                box(idim, 2) = max(box(idim, 2), pt(idim))
            end do
        end do
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function hhoLengthBoundingBoxCell(hhoCell, axes) result(length)
!
        implicit none
!
        type(HHO_Cell), intent(in)    :: hhoCell
        real(kind=8), intent(in)      :: axes(3, 3)
        real(kind=8), dimension(3)    :: length
!
! --------------------------------------------------------------------------------------------------
!  In HHO_Cell              :: cell HHO
!  In axes                  :: axes to use
!  Out length               :: length of the bounding box of the cell
! --------------------------------------------------------------------------------------------------
        real(kind=8), dimension(3, 2) :: box
        integer(kind=8) :: idim, ndim
! --------------------------------------------------------------------------------------------------
        length = 1.d0
        ndim = hhoCell%ndim
!
        box = hhoBoundingBoxCell(hhoCell, axes)
!
        do idim = 1, ndim
            length(idim) = abs(box(idim, 2)-box(idim, 1))
        end do
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function hhoLengthBoundingBoxFace(hhoFace, axes) result(length)
!
        implicit none
!
        type(HHO_Face), intent(in)    :: hhoFace
        real(kind=8), intent(in)      :: axes(3, 2)
        real(kind=8), dimension(2)    :: length
!
! --------------------------------------------------------------------------------------------------
!  In HHO_Face              :: face HHO
!  In axes                  :: axes to use
!  Out length               :: length of the bounding box of the face
! --------------------------------------------------------------------------------------------------
        real(kind=8) :: box(2, 2)
        integer(kind=8) :: idim, ndim
! --------------------------------------------------------------------------------------------------
        length = 1.d0
        ndim = hhoFace%ndim
!
        if (ndim == 1) then
            length(1) = hhoDiameterFace(hhoFace)
        else
            box = hhoBoundingBoxFace(hhoFace, axes)
!
            do idim = 1, ndim
                length(idim) = abs(box(idim, 2)-box(idim, 1))
            end do
        end if
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function hhoCenterBoundingBoxCell(hhoCell, axes) result(center)
!
        implicit none
!
        type(HHO_Cell), intent(in)    :: hhoCell
        real(kind=8), intent(in)      :: axes(3, 3)
        real(kind=8), dimension(3)    :: center
!
! --------------------------------------------------------------------------------------------------
!  In HHO_Cell              :: cell HHO
!  In axes                  :: axes to use
!  Out center               :: center of the bounding box of the cell
! --------------------------------------------------------------------------------------------------
        real(kind=8), dimension(3, 2) :: box, c(3)
        integer(kind=8) :: idim, ndim
! --------------------------------------------------------------------------------------------------
        c = 0.d0
        ndim = hhoCell%ndim
!
        box = hhoBoundingBoxCell(hhoCell, axes)
!
        do idim = 1, ndim
            c(idim) = (box(idim, 1)+box(idim, 2))/2.d0
        end do

        center = matmul(axes, c)
!
    end function
!
!===================================================================================================
!
!===================================================================================================
!
    function hhoCenterBoundingBoxFace(hhoFace, axes) result(center)
!
        implicit none
!
        type(HHO_Face), intent(in)    :: hhoFace
        real(kind=8), intent(in)      :: axes(3, 2)
        real(kind=8), dimension(3)    :: center
!
! --------------------------------------------------------------------------------------------------
!  In HHO_Face              :: face HHO
!  In axes                  :: axes to use
!  Out center               :: center of the bounding box of the face
! --------------------------------------------------------------------------------------------------
        real(kind=8) :: box(2, 2), c(2)
        integer(kind=8) :: idim, ndim
! --------------------------------------------------------------------------------------------------
        ndim = hhoFace%ndim
!
        if (ndim == 1) then
            center = hhoFace%barycenter
        else
            box = hhoBoundingBoxFace(hhoFace, axes)
!
            do idim = 1, ndim
                c(idim) = (box(idim, 2)+box(idim, 1))/2.d0
            end do

            center = matmul(axes, c)
        end if
!
    end function
!
end module
