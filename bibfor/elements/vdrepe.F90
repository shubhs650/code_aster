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
subroutine vdrepe(plateOrie, nomtez, matevn, matevg)
!
    use plate_type
    implicit none
!
#include "asterc/r8dgrd.h"
#include "asterfort/assert.h"
#include "asterfort/coqrep.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "jeveux.h"
!
    type(plateOrie_Para), intent(in) :: plateOrie
    character(len=*), intent(in) :: nomtez
    real(kind=8), intent(out) :: matevn(2, 2, 1), matevg(2, 2, 1)
!
! --------------------------------------------------------------------------------------------------
!
!      VDREPE   -- DETERMINATION DES MATRICES DE PASSAGE
!                  DES REPERES INTRINSEQUES AUX NOEUDS  DE L'ELEMENT
!                  AU REPERE UTILISATEUR (MATRICE MATEVN)
!                  ET DES REPERES INTRINSEQUES AUX POINTS D'INTEGRATION
!                  DE L'ELEMENT AU REPERE UTILISATEUR (MATRICE MATEVG)
!                  POUR LES ELEMENTS DE COQUE EPAISSE 3D .
!
!   ARGUMENT        E/S   TYPE         ROLE
!    NOMTE          IN     K*       NOM DU TYPE D'ELEMENT
!    MATEVN(2,2,10) OUT    R        MATRICES DE PASSAGE DES REPERES
!                                   INTRINSEQUES AUX NOEUDS  DE
!                                   L'ELEMENT AU REPERE UTILISATEUR
!    MATEVG(2,2,10) OUT    R        MATRICES DE PASSAGE DES REPERES
!                                   INTRINSEQUES AUX POINTS
!                                   D'INTEGRATION DE L'ELEMENT AU
!                                   REPERE UTILISATEUR
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16) :: nomte
    real(kind=8) :: pgl(3, 3)
    integer(kind=8) :: i, idec, kpgsr, ino, j, k
    integer(kind=8) :: lzi, lzr, nb2, npgsr
    real(kind=8) :: alpha, beta, c, s
!
! --------------------------------------------------------------------------------------------------
!
    nomte = nomtez

! - Access objects
    call jevete('&INEL.'//nomte(1:8)//'.DESI', ' ', lzi)
    nb2 = zi(lzi-1+2)
    npgsr = zi(lzi-1+3)
    call jevete('&INEL.'//nomte(1:8)//'.DESR', ' ', lzr)

! - RECUPERATION DES ANGLES DETERMINANT LE REPERE UTILISATEUR
! - PAR RAPPORT AU REPERE GLOBAL :
    ASSERT(plateOrie%lRead)
    alpha = plateOrie%alpha
    beta = plateOrie%beta

! - Matrix for nodal quantities - Compute local => global
    idec = 1090
    do ino = 1, nb2
! ----- RECUPERATION DE LA MATRICE DE PASSAGE AU NOEUD COURANT
        k = 0
        do j = 1, 3
            do i = 1, 3
                k = k+1
                pgl(i, j) = zr(lzr+idec+(ino-1)*9+k-1)
            end do
        end do

! ----- Compute operators for coordinate transformation
        call coqrep(pgl, &
                    alpha, beta, &
                    c_=c, s_=s)

!       -- (C,S) N'EST PAS TOUJOURS EXACTEMENT DE NORME=1:
        c = c/sqrt(c*c+s*s)
        s = s/sqrt(c*c+s*s)
!
        matevn(1, 1, ino) = c
        matevn(2, 1, ino) = s
        matevn(1, 2, ino) = -s
        matevn(2, 2, ino) = c
    end do

! - Matrix for point quantities (reduced integration) - Compute local => global
    idec = 2000
    do kpgsr = 1, npgsr
! ----- RECUPERATION DE LA MATRICE DE PASSAGE AU POINT D'INTEGRATION COURANT
        k = 0
        do j = 1, 3
            do i = 1, 3
                k = k+1
                pgl(i, j) = zr(lzr+idec+(kpgsr-1)*9+k-1)
            end do
        end do

! ----- Compute operators for coordinate transformation
        call coqrep(pgl, &
                    alpha, beta, &
                    c_=c, s_=s)

!       -- (C,S) N'EST PAS TOUJOURS EXACTEMENT DE NORME=1:
        c = c/sqrt(c*c+s*s)
        s = s/sqrt(c*c+s*s)
!
        matevg(1, 1, kpgsr) = c
        matevg(2, 1, kpgsr) = s
        matevg(1, 2, kpgsr) = -s
        matevg(2, 2, kpgsr) = c
!
    end do
!
end subroutine
