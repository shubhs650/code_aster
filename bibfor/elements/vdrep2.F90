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
subroutine vdrep2(alphaIn, betaIn, nb2, npgsr, desr, &
                  matevn, matevg_)
!
    implicit none
!
#include "asterf_types.h"
#include "asterc/r8dgrd.h"
#include "asterfort/coqrep.h"
!
    real(kind=8), intent(in) :: alphaIn, betaIn
    integer(kind=8), intent(in) :: nb2, npgsr
    real(kind=8), intent(in) :: desr(*)
    real(kind=8), intent(out) :: matevn(2, 2, 1)
    real(kind=8), optional, intent(out) :: matevg_(2, 2, 1)
!
! --------------------------------------------------------------------------------------------------
!
!      VDREP2   -- DETERMINATION DES MATRICES DE PASSAGE
!                  DES REPERES INTRINSEQUES AUX NOEUDS  DE L'ELEMENT
!                  AU REPERE UTILISATEUR (MATRICE MATEVN)
!                  ET DES REPERES INTRINSEQUES AUX POINTS D'INTEGRATION
!                  DE L'ELEMENT AU REPERE UTILISATEUR (MATRICE MATEVG)
!                  POUR LES ELEMENTS DE COQUE EPAISSE 3D .
!
!   ARGUMENT        E/S   TYPE         ROLE
!    ALPHA, BETA    IN     R    ANGLES DETERMINANT LE REPERE UTILISATEUR
!    MATEVN(2,2,10) OUT    R        MATRICES DE PASSAGE DES REPERES
!                                   INTRINSEQUES AUX NOEUDS  DE
!                                   L'ELEMENT AU REPERE UTILISATEUR
!    MATEVG(2,2,10) OUT    R        MATRICES DE PASSAGE DES REPERES
!                                   INTRINSEQUES AUX POINTS
!                                   D'INTEGRATION DE L'ELEMENT AU
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i, idec, ino, j, k, kpgsr
    real(kind=8) :: c, s
    real(kind=8) :: pgl(3, 3), alpha, beta
!
! --------------------------------------------------------------------------------------------------
!
    alpha = alphaIn*r8dgrd()
    beta = betaIn*r8dgrd()
    matevn = 0.d0

! - DETERMINATION DES MATRICES DE PASSAGE DES REPERES INTRINSEQUES
! - AUX NOEUDS DE L'ELEMENT AU REPERE UTILISATEUR
    idec = 1090

    do ino = 1, nb2
! ---   RECUPERATION DE LA MATRICE DE PASSAGE AU NOEUD COURANT
        k = 0
        do j = 1, 3
            do i = 1, 3
                k = k+1
                pgl(i, j) = desr(idec+(ino-1)*9+k)
            end do
        end do

! ----- Compute operators for coordinate transformation
        call coqrep(pgl, alpha, beta, &
                    c_=c, s_=s)
        matevn(1, 1, ino) = c
        matevn(2, 1, ino) = s
        matevn(1, 2, ino) = -s
        matevn(2, 2, ino) = c
!
    end do

! - Matrix for point quantities (reduced integration) - Compute local => global
    if (present(matevg_)) then
        idec = 2000
        matevg_ = 0.d0
        do kpgsr = 1, npgsr
! ----- RECUPERATION DE LA MATRICE DE PASSAGE AU POINT D'INTEGRATION COURANT
            k = 0
            do j = 1, 3
                do i = 1, 3
                    k = k+1
                    pgl(i, j) = desr(idec+(kpgsr-1)*9+k)
                end do
            end do

! ----- Compute operators for coordinate transformation
            call coqrep(pgl, &
                        alpha, beta, &
                        c_=c, s_=s)
!
            matevg_(1, 1, kpgsr) = c
            matevg_(2, 1, kpgsr) = s
            matevg_(1, 2, kpgsr) = -s
            matevg_(2, 2, kpgsr) = c
!
        end do
    end if
!
end subroutine
