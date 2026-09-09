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
subroutine mpglcp(typecp, nbnolo, coordo, alpha, beta, &
                  gamma, pgl)
!
    implicit none
!
#include "asterfort/angvx.h"
#include "asterfort/assert.h"
#include "asterfort/coqrep.h"
#include "asterfort/dxqpgl.h"
#include "asterfort/dxtpgl.h"
#include "asterfort/matrot.h"
#include "jeveux.h"
!
    character(len=1) :: typecp
    integer(kind=8) :: nbnolo
    real(kind=8) :: coordo(*), alpha, beta, gamma, pgl(3, 3)
!
! --------------------------------------------------------------------------------------------------
!
!  CALCUL DE LA MATRICE DE PASSAGE GLOBAL -> LOCAL COQUES ET POUTRES
!
! --------------------------------------------------------------------------------------------------
!
!
!  ROUTINE CALCUL DE LA MATRICE DE PASSAGE DU REPERE GLOBAL AU REPERE
!    LOCAL DANS LE CAS DES COQUES ET DES POUTRES
!
! IN  :
!   TYPECP  K1   'P' POUR POUTRES
!                'D' pour discrets à 2 noeuds (idem poutres)
!                '1' pour les elements à 1 noeud (alpha,beta,gamma sont donnés)
!                'C' POUR COQUES
!   NBNOLO  I    NOMBRE DE NOEUDS.
!                   2 pour les poutres
!                   2 ou 1 pour les discrets
!                   3 ou 4 pour les coques
!   COORDO  R*   TABLEAU CONTENANT LES COORDOONEES DES NOEUDS
!                  DIMENSIONNE A NBNOLO*3
!   ALPHA   R    PREMIER ANGLE NAUTIQUE
!   BETA    R    DEUXIEME ANGLE NAUTIQUE
!   GAMMA   R    TROISIEME ANGLE NAUTIQUE
!
!  POUR LES POUTRES, SEUL GAMMA EST A FOURNIR (CAR ALPHA ET BETA SONT
!    RECALCULES A PARTIR DES COORDONNEES)
!  POUR LES COQUES, SEULS ALPHA ET BETA SONT A FOURNIR
!
! OUT :
!   PGL     R*   LA MATRICE DE PASSAGE DE DIMENSION 3*3
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: xd(3), angl(3), alphal, betal, t2ui(2, 2)
    real(kind=8) :: mat1(3, 3), mat2(3, 3)
!
! --------------------------------------------------------------------------------------------------
!
    if ((typecp .eq. 'P') .or. (typecp .eq. 'D')) then
        ASSERT((nbnolo .eq. 2) .or. (nbnolo .eq. 3))
!       CALCUL DE ALPHA ET BETA
        xd = coordo(1:3)-coordo(4:6)
        call angvx(xd, alphal, betal)
        angl(1) = alphal
        angl(2) = betal
        angl(3) = gamma
        call matrot(angl, pgl)
!
    else if (typecp .eq. '1') then
        ASSERT(nbnolo .eq. 1)
        angl(1) = alpha
        angl(2) = beta
        angl(3) = gamma
        call matrot(angl, pgl)
!
    else if (typecp .eq. 'C') then
! ----- Calculate the transformation: global coordinate system/intrinsic coordinate system
        pgl = 0.d0
        if (nbnolo .eq. 3) then
            call dxtpgl(coordo, pgl)
        elseif (nbnolo .eq. 4) then
            call dxqpgl(coordo, pgl)
        else
            ASSERT(ASTER_FALSE)
        end if

! ----- Compute operators for coordinate transformation
        call coqrep(pgl, alpha, beta, t2ui_=t2ui)

! ----- Inverse: from intrinsic to global
        mat1(1, 1) = pgl(1, 1)
        mat1(1, 2) = pgl(2, 1)
        mat1(1, 3) = pgl(3, 1)
        mat1(2, 1) = pgl(1, 2)
        mat1(2, 2) = pgl(2, 2)
        mat1(2, 3) = pgl(3, 2)
        mat1(3, 1) = pgl(1, 3)
        mat1(3, 2) = pgl(2, 3)
        mat1(3, 3) = pgl(3, 3)

        mat2(1, 1) = t2ui(1, 1)
        mat2(1, 2) = t2ui(2, 1)
        mat2(1, 3) = 0.d0
        mat2(2, 1) = t2ui(1, 2)
        mat2(2, 2) = t2ui(2, 2)
        mat2(2, 3) = 0.d0
        mat2(3, 1) = 0.d0
        mat2(3, 2) = 0.d0
        mat2(3, 3) = 1.d0
        pgl = matmul(mat1, mat2)
    end if
!
end subroutine
