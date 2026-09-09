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
subroutine te0028(option, nomte)
!
    use Behaviour_module, only: behaviourOption
    use c_interface_plaq_mitc_j
    use FE_basis_module
    use FE_quadrature_module
    use FE_topo_module
    use iso_c_binding
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystNone
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/dxroep.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/get_elas_id.h"
#include "asterfort/get_elas_para.h"
#include "asterfort/jevech.h"
#include "asterfort/writeMatrix.h"
#include "FE_module.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: PLAQ_MITC
!
! Options: FULL_MECA_*, RIGI_MECA_*, RAPH_MECA
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ndim, nno, nnos, npg, jpoids, jvf, jdfde, jvGeom, jvMaterc
    integer(kind=8) :: i, j, k, elas_id, igau
    real(c_double) :: e, nu, epais, rho, e_, nu_
    real(c_double) :: tiny
    character(len=8), parameter :: fami = 'RIGI'
    character(len=16) :: elas_keyword
    real(c_double) :: cst(4), coor(27), cdofs_f(12)
    ! NOTICE: see the size of the arrays in the C file: c_interface_plaq_mitc_j
    integer(c_int) :: ncst, ncd, ne0, ne1, ne2, ne3
    integer(c_int) :: entities0(1), entities1(1), entities2(1), entities3(1)
    integer(kind=4), parameter :: size_init = 30, size_final = 22
    real(c_double) :: w_0(size_init), signs(size_final)
    integer(kind=8) :: reorder(size_final)
    real(c_double), dimension(size_init, size_init) :: M_elem
    real(c_double), dimension(size_final, size_final) :: M_cond, M_final
    real(c_double), dimension(size_init*size_init) :: M0, M1, M2, M3, M4
    real(c_double) :: AA(22, 22), BB(4, 4), CC(22, 4), DD(4, 4), CDinv(22, 4)
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!
    w_0 = 0.d0
    tiny = 1.d-8

! - Get plate parameters
    call getCara(plateCara, plateOrie)

! - No global<=>local transformation
    call compCoorSystNone(plateOrie)

! - Finite element informations
    call elrefe_info(fami=fami, ndim=ndim, nno=nno, nnos=nnos, &
                     npg=npg, jpoids=jpoids, jvf=jvf, jdfde=jdfde)

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Material parameters
    call jevech('PMATERC', 'L', jvMaterc)

! ----- Get elastic parameters (only isotropic elasticity)
! FIXME: for instance E and NU are supposed to be equal for all quadpoints
!
    call get_elas_id(zi(jvMaterc), elas_id, elas_keyword)
    call get_elas_para(fami, zi(jvMaterc), '+', 1, 1, &
                       elas_id, elas_keyword, &
                       e_=e, nu_=nu)

    do igau = 2, npg
        call get_elas_id(zi(jvMaterc), elas_id, elas_keyword)
        call get_elas_para(fami, zi(jvMaterc), '+', igau, 1, &
                           elas_id, elas_keyword, &
                           e_=e_, nu_=nu_)
        ASSERT(abs(e-e_) .le. tiny .and. abs(nu-nu_) .le. tiny)
    end do

    call dxroep(plateCara, rho, epais)

! Fill material parameters vector
    cst(1) = e
    cst(2) = nu
    cst(3) = 5.0/6.0
    cst(4) = epais
    ncst = size(cst)
!
! Remplissage du vecteur de coordonnées (3 coordonnées par nœud)
    do i = 0, 8
        coor(3*i+1) = zr(jvGeom+3*i)
        coor(3*i+2) = zr(jvGeom+3*i+1)
        coor(3*i+3) = zr(jvGeom+3*i+2)
    end do

! Remplissage des coordonées (N1, N4, N2, N3)
    cdofs_f(1) = coor(1)
    cdofs_f(2) = coor(2)
    cdofs_f(3) = coor(3)
!
    cdofs_f(4) = coor(10)
    cdofs_f(5) = coor(11)
    cdofs_f(6) = coor(12)
!
    cdofs_f(7) = coor(4)
    cdofs_f(8) = coor(5)
    cdofs_f(9) = coor(6)
!
    cdofs_f(10) = coor(7)
    cdofs_f(11) = coor(8)
    cdofs_f(12) = coor(9)
!
    ncd = size(cdofs_f)
!
! On définit 4 entities car 4 arêtes dans le quadrangle
    entities0(1) = 0
    entities1(1) = 1
    entities2(1) = 2
    entities3(1) = 3
    ne0 = 1
    ne1 = 1
    ne2 = 1
    ne3 = 1
!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!MATRICE DE RIGIDITÉ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Appel de la fonction C++ - Matrice de rigidité -> Part1
    call BP4_qu9_Fortran(w_0, size_init, cdofs_f, ncd, entities0, ne0, cst, ncst, M0)

! Appel de la fonction C++ - Matrice de rigidité -> Part2
    call BP5_qu9_Fortran(w_0, size_init, cdofs_f, ncd, entities0, ne0, cst, ncst, M1)
    call BP5_qu9_Fortran(w_0, size_init, cdofs_f, ncd, entities1, ne1, cst, ncst, M2)
    call BP5_qu9_Fortran(w_0, size_init, cdofs_f, ncd, entities2, ne2, cst, ncst, M3)
    call BP5_qu9_Fortran(w_0, size_init, cdofs_f, ncd, entities3, ne3, cst, ncst, M4)

! Remplissage de la matrice
    do i = 1, size_init
        do j = 1, size_init
            k = (j-1)*size_init+i
            M_elem(i, j) = M0(k)+M1(k)+M2(k)+M3(k)+M4(k)
        end do
    end do

! Recupérer les bloques
    ! (z, z) où z = (thetha_x theta_y w)
    AA = M_elem(1:22, 1:22)

    ! (gamma, gamma)
    BB = M_elem(23:26, 23:26)

    ! (z, p) où z = (thetha_x theta_y w)
    CC = M_elem(1:22, 27:30)

    ! (gamma, p)
    DD = M_elem(23:26, 27:30)

    do j = 1, 4
        CDinv(:, j) = CC(:, j)/DD(j, j)
    end do

    ! Condensation statique
    M_cond = AA+matmul(matmul(CDinv, BB), transpose(CDinv))

! Reorganisation de la matrice
    ![w1, θ_y1, -θ_x1, w2, θ_y2, -θ_x2, w3, θ_y3, -θ_x3, w4, θ_y4, -θ_x4,
    ! θ_y5, -θ_x5, θ_y6, -θ_x6, θ_y7, -θ_x7, θ_y8, -θ_x8, θ_y9, -θ_x9]
    !
    reorder = (/ &
              19, 2, 1, &
              21, 6, 5, &
              22, 8, 7, &
              20, 4, 3, &
              12, 11, &
              16, 15, &
              14, 13, &
              10, 9, &
              18, 17 &
              /)
    signs = (/ &
            1.d0, 1.d0, -1.d0, &
            1.d0, 1.d0, -1.d0, &
            1.d0, 1.d0, -1.d0, &
            1.d0, 1.d0, -1.d0, &
            1.d0, -1.d0, &
            1.d0, -1.d0, &
            1.d0, -1.d0, &
            1.d0, -1.d0, &
            1.d0, -1.d0 &
            /)
!
    do i = 1, size_final
        do j = 1, size_final
            M_final(i, j) = signs(i)*signs(j)*M_cond(reorder(i), reorder(j))
        end do
    end do
!
    call writeMatrix('PMATUUR', size(M_final, dim=1), size(M_final, dim=2), ASTER_TRUE, M_final)
!
end subroutine
