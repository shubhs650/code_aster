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
! aslint: disable=W0413
!
subroutine nmcham(fami, kpg, ksp, imate, relaComp, &
                  matel, mat, nbvar, memo, visc, &
                  idelta, coef)
!
    implicit none
!
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
#include "asterfort/verift.h"
!
    integer(kind=8) :: imate, nbvar, kpg, ksp, memo, visc, idelta, nrad
    character(len=16), intent(in) :: relaComp
    real(kind=8) :: mat(18), matel(4)
    character(len=*) :: fami
!
!      NMCHAM   -- COEFFICIENTS MATERIAU DES LOIS DE COMPORTEMENT
!                  'VMIS_CINx_CHAB'  'VISC_CINx_CHAB'
!                  'VMIS_CIN2_MEMO'  'VISC_CIN2_MEMO'
!
!     ARGUMENT  E/S  TYPE         ROLE
!     IMATE     IN    I    ADRESSE DU MATERIAU CODE
!     COMPOR    IN    K16  COMPOR(1) NOM DU COMPORTEMENT
!     MAT       OUT   R    COEF MATERIAU
!

    real(kind=8) :: coef, valres(12), c2inf, gamm20, delta1, delta2
    real(kind=8) :: r0, rinf, b, cinf, k, w, gamma0, epsi
    real(kind=8) :: un, ainf, kvi, valden, unskvi
    integer(kind=8) :: icodre(12)
    character(len=16) :: valk(2), texte(2), nomres(12)
    character(len=8) :: nomemo(4)
!.========================= DEBUT DU CODE EXECUTABLE ==================
!
    nbvar = 0
    if (relaComp(6:9) .eq. 'CIN1') then
        nbvar = 1
    else if (relaComp(6:9) .eq. 'CIN2') then
        nbvar = 2
    else if (relaComp(6:9) .eq. 'MEMO') then
        nbvar = 2
    else
        ASSERT(ASTER_FALSE)
    end if
!
    if (relaComp(1:4) .eq. 'VMIS') then
        visc = 0
    else if (relaComp(1:4) .eq. 'VISC') then
        visc = 1
    else
        ASSERT(ASTER_FALSE)
    end if
!
    memo = 0
    if (relaComp(11:14) .eq. 'MEMO') then
        memo = 1
    else if (relaComp(6:9) .eq. 'MEMO') then
        memo = 1
    end if
!
    if (relaComp(11:14) .eq. 'NRAD') then
        nrad = 1
    else
        nrad = 0
    end if
!
! --- INITIALISATIONS :
!     ===============
    un = 1.0d0
!
    call verift(fami, kpg, ksp, 'T', imate, &
                epsth_=coef)
!
! --- RECUPERATION DES CARACTERISTIQUES ELASTIQUES :
!     ============================================
    nomres(1) = 'E'
    nomres(2) = 'NU'
!
! ---  CARACTERISTIQUES A L'INSTANT PRECEDENT :
!      --------------------------------------
    call rcvalb(fami, kpg, ksp, '-', imate, &
                ' ', 'ELAS', 0, ' ', [0.d0], &
                2, nomres(1), matel(1), icodre(1), 1)
!
! ---  CARACTERISTIQUES A L'INSTANT ACTUEL :
!      -----------------------------------
!
    call rcvalb(fami, kpg, ksp, '+', imate, &
                ' ', 'ELAS', 0, ' ', [0.d0], &
                2, nomres(1), matel(3), icodre(1), 1)
!
! --- RECUPERATION DES CARACTERISTIQUES D'ECROUISSAGE :
!     ===============================================
    nomres(1) = 'R_0'
    nomres(2) = 'R_I'
    nomres(3) = 'B'
    nomres(5) = 'K'
    nomres(6) = 'W'
    nomres(8) = 'A_I'
    if (nbvar .eq. 1) then
        nomres(4) = 'C_I'
        nomres(7) = 'G_0'
        call rcvalb(fami, kpg, ksp, '+', imate, &
                    ' ', 'CIN1_CHAB', 0, ' ', [0.d0], &
                    8, nomres, valres, icodre, 1)
    else if (nbvar .eq. 2) then
        nomres(4) = 'C1_I'
        nomres(7) = 'G1_0'
        nomres(9) = 'C2_I'
        nomres(10) = 'G2_0'
        call rcvalb(fami, kpg, ksp, '+', imate, &
                    ' ', 'CIN2_CHAB', 0, ' ', [0.d0], &
                    10, nomres, valres, icodre, 1)
    end if
!
    r0 = valres(1)
    rinf = valres(2)
    b = valres(3)
    cinf = valres(4)
    k = valres(5)
    w = valres(6)
    gamma0 = valres(7)
    ainf = valres(8)
    c2inf = 0.d0
    gamm20 = 0.d0
    if (nbvar .eq. 2) then
        c2inf = valres(9)
        gamm20 = valres(10)
    end if
!
    mat(1) = r0
    mat(2) = rinf
    mat(3) = b
    mat(4) = cinf
    mat(5) = k
    mat(6) = w
    mat(7) = gamma0
    mat(8) = ainf
!
    if (b < 0.d0) then
        valk(1) = relaComp
        valk(2) = 'b'
        call utmess('A', 'COMPOR1_84', nk=2, valk=valk, sr=b)
    end if
!
    if (w < 0.d0) then
        valk(1) = relaComp
        valk(2) = 'w'
        call utmess('A', 'COMPOR1_84', nk=2, valk=valk, sr=w)
    end if
!
!     IDELTA : TYPE DE NON PROPORTIONNALITE
    if (nrad .eq. 1) then
        nomres(1) = 'DELTA1'
        nomres(2) = 'DELTA2'
        call rcvalb(fami, kpg, ksp, '+', imate, &
                    ' ', 'CIN2_NRAD', 0, ' ', [0.d0], &
                    2, nomres, valres, icodre, 1)
        delta1 = valres(1)
        delta2 = valres(2)
        mat(17) = delta1
        mat(18) = delta2
        idelta = 0
        epsi = r8prem()
        if (abs(delta1-1.d0) .gt. epsi) then
            if (abs(delta2-1.d0) .gt. epsi) then
                idelta = 3
            else
                idelta = 1
            end if
        else
            if (abs(delta2-1.d0) .gt. epsi) then
                idelta = 2
            end if
        end if
!        UTILE POUR LE CAS DES FONCTIONS
        if ((delta1 .gt. (1.d0+epsi)) .or. (delta1 .lt. -epsi)) then
            call utmess('F', 'COMPOR1_80', sr=delta1)
        end if
        if ((delta2 .gt. (1.d0+epsi)) .or. (delta2 .lt. -epsi)) then
            call utmess('F', 'COMPOR1_81', sr=delta2)
        end if
    else
        delta1 = 1.d0
        delta2 = 1.d0
        mat(17) = delta1
        mat(18) = delta2
    end if
!
!
    if (nbvar .eq. 2) then
        mat(9) = c2inf
        mat(10) = gamm20
    else
        mat(9) = 0.d0
        mat(10) = 0.d0
    end if
!
    if (visc .eq. 1) then
!
! ---    RECUPERATION DES CARACTERISTIQUES VISQUEUSES :
!        ============================================
        nomres(1) = 'N'
        nomres(2) = 'UN_SUR_K'
        nomres(3) = 'UN_SUR_M'
        call rcvalb(fami, kpg, ksp, '+', imate, &
                    ' ', 'LEMAITRE', 0, ' ', [0.d0], &
                    3, nomres, valres, icodre, 0)
!
        if (icodre(1) .eq. 0) then
            valden = valres(1)
            unskvi = valres(2)
            if (valden .le. 0.d0) then
                call utmess('F', 'ALGORITH6_67')
            end if
            if (unskvi .eq. 0.d0) then
                call utmess('F', 'ALGORITH6_68')
            end if
            kvi = un/unskvi
            if (valres(3) .ne. 0.d0) then
                call utmess('F', 'ALGORITH6_69')
            end if
        else
            texte(1) = relaComp
            texte(2) = "VMIS"//relaComp(5:16)
            call utmess('F', 'COMPOR1_32', nk=2, valk=texte)
        end if
    else
        valden = un
        kvi = 0.d0
    end if
    mat(11) = valden
    mat(12) = kvi
    mat(13) = 0.d0
    mat(14) = 0.d0
    mat(15) = 0.d0
    mat(16) = 0.d0
!
    if (memo .eq. 1) then
! ---    RECUPERATION DES CARACTERISTIQUES MEMOIRE :
!        ============================================
        nomemo(1) = 'ETA     '
        nomemo(2) = 'Q_0     '
        nomemo(3) = 'Q_M     '
        nomemo(4) = 'MU      '
        call rcvalb(fami, kpg, ksp, '+', imate, &
                    ' ', 'MEMO_ECRO', 0, ' ', [0.d0], &
                    4, nomemo, mat(13), icodre, 1)
!
    end if
end subroutine
