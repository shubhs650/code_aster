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
subroutine pipetu(ndim, mate, sup, sud, vim, &
                  dtau, copilo)
!
!
    implicit none
#include "asterc/r8gaem.h"
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/utmess.h"
#include "asterfort/rcvalb.h"
#include "asterfort/zerop2.h"
#include "blas/ddot.h"
    integer(kind=8), intent(in) :: ndim
    integer(kind=8), intent(in) :: mate
    real(kind=8), intent(in) :: sup(:)
    real(kind=8), intent(in) :: sud(:)
    real(kind=8), intent(in) :: vim(:)
    real(kind=8), intent(in) :: dtau
    real(kind=8), intent(out) :: copilo(:)
!
!-----------------------------------------------------------------------
!
! PILOTAGE PRED_ELAS POUR LA LOI COHESIVE CZM_TURON
! DE L'ELEMENT DE JOINT (2D ET 3D)
! SAUT DE DEPLACMENT : SU = SUP + eta*SUD
!
!-----------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbpa = 6
    integer(kind=8) :: i, j, nrac, ok(4), nsol
    real(kind=8) :: p0, p1, p2, rac(2), eta(4), a0(4), a1(4), tmp
    real(kind=8) :: lc, k0, ka, kref, c, val(nbpa), etasol(4), xn
    real(kind=8) :: k, delta_N_0, delta_T_0, delta_N_f, delta_T_f, para_bk
    real(kind=8) :: delta_0, delta_f
    integer(kind=8) :: cod(nbpa), kpg, spt
    character(len=16) :: nom(nbpa)
    character(len=8) :: fami, poum
    blas_int :: b_incx, b_incy, b_n
!
    data nom/'K', 'SIGM_C_N', 'SIGM_C_T', 'GC_N', 'GC_T', 'ETA_BK'/
!
!-----------------------------------------------------------------------
!
! INITIALISATION
    ok(1) = 0
    ok(2) = 0
    ok(3) = 0
    ok(4) = 0
    fami = 'FPG1'
    kpg = 1
    spt = 1
    poum = '+'
!
! RECUPERATION DES PARAMETRES PHYSIQUES
    call rcvalb(fami, kpg, spt, poum, mate, &
                ' ', 'RUPT_TURON', 0, ' ', [0.d0], &
                5, nom, val, cod, 2)
    k = val(1)
    delta_N_0 = val(2)/k
    delta_T_0 = val(3)/k
    delta_N_f = 2*val(4)/(k*delta_N_0)
    delta_T_f = 2*val(5)/(k*delta_T_0)
    para_bk = val(6)
!!
!! SAUTS PILOTES : PARTIES NORMALE POSITIVE ET TANGENTIELLE
!    sud_T = sqrt(sud(2)**2 + sud(3)**2)
!    sud_N_pos = max(0.d0,sud(1))
!    l_sud0 = ((sud_N_pos + sud_T) .lt. r8prem())
!!
!! CALCUL DU TAUX DE MIXITE A T+ (PAR LES SAUTS DE LA PARTIE PILOTEE)
!    if (.not. l_sud0) then
!        beta = sud_T / (sud_N_pos + sud_T)
!!
!! CALCUL DU TAUX DE MIXITE A T+ (PAR LES TAUX DE RESTITUTION D'ENERGIE)
!        b = beta**2/(1-2*beta+2*beta**2)
!!
!! CALCUL DU SEUIL D'INITIATION DE L'ENDOMMAGEMENT EN MODE MIXTE A T+
!        delta_0 = sqrt(delta_N_0**2 + (delta_T_0**2-delta_N_0**2)*b**para_bk)
!        ASSERT(delta_0 .gt. r8prem())
!!   * version alternative (critère de Ye) --> à implémenter
!!
!! CALCUL DU SEUIL DE PROPAGATION DE LA FISSURE EN MODE MIXTE A T+
!        delta_f = 1/delta_0 * (delta_N_0*delta_N_f + &
!                (delta_T_0*delta_T_f - delta_N_0*delta_N_f )*b**para_bk)
!        ASSERT(delta_f .gt. r8prem())
!    else
!!   * CAS PARTICULIERS SI ON (RE)PASSE PAR UN ETAT DE SAUT NUL :
!!   * ON RESTE AVEC LES VALEURS PRECEDENTES
!        beta = vim(13)
!        b = vim(14)
!        delta_0 = vim(15)
!        delta_f = vim(16)
!!   * POUR LE PREMIER PAS DE TEMPS SANS ETAT INIT
!        if (delta_0 .lt. r8prem()) then
!            call utmess('A', 'COMPOR4_73')
!            delta_0 = min(delta_N_0, delta_T_0)
!            delta_f = min(delta_N_f, delta_T_f)
!        endif
!    endif
    delta_0 = vim(15)
    delta_f = vim(16)
!
!! PARAMETRES POUR L'EQUATION DE PILOTAGE (inspiré de pipeba.F90)
!    lc = delta_f/2
!    k0 = delta_0/2
    lc = delta_f
    k0 = delta_0
    ka = max(k0, vim(1))
    kref = max(lc, ka)
!
    c = dtau*kref+ka
!
!    RESOLUTION FEL(ETA) = DTAU
!    OU FEL(ETA) = ( SQRT(P0 + 2 P1 ETA + P2 ETA**2) - KA) / KREF
!    PORTION EN COMPRESSION : FEL = (ABS(SU(2)) - KA ) / KREF
!    ON INCLUT EGALEMENT UN SAFE-GUARD SU_N > -KREF CAR AU-DELA CE SONT
!    DES SOLUTIONS TRES FORTEMENT EN COMPRESSION QUI FONT EXPLOSER LA
!    PENALISATION
!
    p0 = 0.d0
    p1 = 0.d0
    p2 = 0.d0
    do i = 2, ndim
        p2 = p2+sud(i)*sud(i)
        p1 = p1+sud(i)*sup(i)
        p0 = p0+sup(i)*sup(i)
    end do
!
!    PAS DE SOLUTION
    if (p2 .lt. (1.d0/r8gaem()**0.5d0)) goto 100
!
!    RECHERCHE DES SOLUTIONS
    call zerop2(2*p1/p2, (p0-c**2)/p2, rac, nrac)
    if (nrac .le. 1) goto 100
!
    xn = sup(1)+rac(2)*sud(1)
    if (xn .le. 0 .and. xn .ge. -kref) then
        ok(1) = 1
        eta(1) = rac(2)
        a1(1) = (p1+p2*eta(1))/(kref*c)
        a0(1) = dtau-eta(1)*a1(1)
    end if
!
    xn = sup(1)+rac(1)*sud(1)
    if (xn .le. 0 .and. xn .ge. -kref) then
        ok(2) = 1
        eta(2) = rac(1)
        a1(2) = (p1+p2*eta(2))/(kref*c)
        a0(2) = dtau-eta(2)*a1(2)
    end if
!
100 continue
!
!
!    PORTION EN TRACTION : FEL = (SQR(SU(1)**2 + SU(2)**2) - KA) / KREF
!
    b_n = to_blas_int(ndim)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    p2 = ddot(b_n, sud, b_incx, sud, b_incy)
    b_n = to_blas_int(ndim)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    p1 = ddot(b_n, sud, b_incx, sup, b_incy)
    b_n = to_blas_int(ndim)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    p0 = ddot(b_n, sup, b_incx, sup, b_incy)
!
!    PAS DE SOLUTION
    if (p2 .lt. (1.d0/r8gaem()**0.5d0)) goto 200
!
!    RECHERCHE DES SOLUTIONS
    call zerop2(2*p1/p2, (p0-c**2)/p2, rac, nrac)
    if (nrac .le. 1) goto 200
!
    if (sup(1)+rac(2)*sud(1) .gt. 0) then
        ok(3) = 1
        eta(3) = rac(2)
        a1(3) = (p1+p2*eta(3))/(kref*c)
        a0(3) = dtau-eta(3)*a1(3)
    end if
!
    if (sup(1)+rac(1)*sud(1) .gt. 0) then
        ok(4) = 1
        eta(4) = rac(1)
        a1(4) = (p1+p2*eta(4))/(kref*c)
        a0(4) = dtau-eta(4)*a1(4)
    end if
!
200 continue
!
!
! -- CLASSEMENT DES SOLUTIONS
!
    nsol = ok(1)+ok(2)+ok(3)+ok(4)
    ASSERT(nsol .le. 2)
!
    j = 0
    do i = 1, 4
        if (ok(i) .eq. 1) then
            j = j+1
            etasol(j) = eta(i)
            copilo(1+2*(j-1)) = a0(i)
            copilo(2+2*(j-1)) = a1(i)
        end if
    end do
!
!    ON RANGE LES SOLUTIONS DANS L'ORDRE CROISSANT (SI NECESSAIRE)
    if (nsol .eq. 2) then
        if (etasol(2) .lt. etasol(1)) then
            tmp = etasol(2)
            etasol(2) = etasol(1)
            etasol(1) = tmp
!
            tmp = copilo(1)
            copilo(1) = copilo(3)
            copilo(3) = tmp
!
            tmp = copilo(2)
            copilo(2) = copilo(4)
            copilo(4) = tmp
        end if
    end if
!
!
!    TRAITEMENT EN L'ABSENCE DE SOLUTION
    if (nsol .eq. 0) then
!    SI DEPLACEMENT PILOTE NUL ET SAUT EQ INFERIEUR A DTAU
!    ON IGNORE LE POINT POUR LA RESOLUTION GLOBALE
        if (p2 .le. (1.d0/r8gaem()**0.5d0) .and. (sqrt(p0)) .le. dtau) then
            copilo(1) = 0.d0
            copilo(2) = 0.d0
            copilo(3) = 0.d0
            copilo(4) = 0.d0
! DANS LES AUTRE CAS COURBE TJS SUPERIEURE A DTAU : ON PLANTE
        else
            copilo(1) = 1.d0
            copilo(5) = 1.d0
        end if
    end if
!
end subroutine
