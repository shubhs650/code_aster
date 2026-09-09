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
subroutine caurtg(nomte, ncmp, sigmau, sigrtg)
!
    use plate_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/marota.h"
#include "asterfort/tecach.h"
#include "asterfort/utbtab.h"
#include "jeveux.h"
!
    integer(kind=8), intent(in) :: ncmp
    character(len=16), intent(in) :: nomte
    real(kind=8), intent(in) :: sigmau(ncmp, 1)
    real(kind=8), intent(out) :: sigrtg(ncmp, 1)
!
! --------------------------------------------------------------------------------------------------
!
!      CAURTG  -- PASSAGE DES CONTRAINTES DE CAUCHY SIGMAU
!                 CALCULEES DANS LE REPERE UTILISATEUR
!                 VERS LE REPERE UTILISATEUR TOURNE DE
!                 LA ROTATION FAISANT PASSER DE L'ETAT
!                 INITIAL A L'ETAT DEFORME DANS LE CAS GROT_GDEP .
!                 SIGRTG DESIGNE LES CONTRAINTES DE CAUCHY DANS
!                 CE DERNIER REPERE .
!
!   ARGUMENT        E/S  TYPE         ROLE
!    NOMTE          IN     K16      NOM DU TYPE D'ELEMENT
!    NCMP           IN     I        NOMBRE DE COMPOSANTES DU TENSEUR
!                                   DES CONTRAINTES
!    SIGMAU(NCMP,1) IN     R        VECTEUR DES CONTRAINTES
!                                   DE CAUCHY DANS LE REPERE UTILISATEUR
!    SIGRTG(NCMP,1) VAR    R        VECTEUR DES CONTRAINTES DE CAUCHY
!                                   TOURNE DU REPERE UTILISATEUR VERS
!                                   LE REPERE TRANSFORME DU REPERE
!                                   UTILISATEUR PAR LA GRANDE ROTATION
!                                   CALCULEE POUR COMP_ELAS - GROT_GDEP
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8) :: vectRota(9, 3)
    real(kind=8) :: xab(3, 3), sigmad(3, 3), sigmat(3, 3)
    real(kind=8) :: drot(3, 3), tetag(3)
    integer(kind=8) :: i, jvDisp, jvGeom, ii, in, iret
    integer(kind=8) :: lzi, lzr, nb1, nb2
!
! --------------------------------------------------------------------------------------------------
!

! - Access to static objects of COQUE_3D
    call jevete('&INEL.'//nomte(1:8)//'.DESI', ' ', lzi)
    nb1 = zi(lzi+1-1)
    nb2 = zi(lzi+2-1)
    call jevete('&INEL.'//nomte(1:8)//'.DESR', ' ', lzr)

! - Get displacements
    call tecach('OOO', 'PDEPLAR', 'L', iret, iad=jvDisp)

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Separation of displacements and rotations
    vectRota = 0.d0
    do in = 1, nb1
        do ii = 1, 3
            vectRota(in, ii) = zr(jvDisp+6*(in-1)+ii+3-1)
        end do
    end do
    do ii = 1, 3
        vectRota(nb2, ii) = zr(jvDisp+6*nb1+ii-1)
    end do

    do i = 1, nb2
        tetag(1) = vectRota(i, 1)
        tetag(2) = vectRota(i, 2)
        tetag(3) = vectRota(i, 3)
        call marota(tetag, drot)
        sigmat(1, 1) = sigmau(1, i)
        sigmat(2, 2) = sigmau(2, i)
        sigmat(3, 3) = sigmau(3, i)
        sigmat(1, 2) = sigmau(4, i)
        sigmat(2, 1) = sigmat(1, 2)
        if (ncmp .eq. 6) then
            sigmat(1, 3) = sigmau(5, i)
            sigmat(2, 3) = sigmau(6, i)
            sigmat(3, 1) = sigmat(1, 3)
            sigmat(3, 2) = sigmat(2, 3)
        end if

! ---   ROTATION DU TENSEUR DES CONTRAINTES DE CAUCHY DE LA
! ---   ROTATION FAISANT PASSER DE L'ETAT INITAL A L'ETAT DEFORME
        call utbtab('ZERO', 3, 3, sigmat, drot, xab, sigmad)

! ---   AFFECTATION DU VECTEUR EN SORTIE DES CONTRAINTES
! ---   DE CAUCHY DANS LE REPERE UTILISATEUR TOURNE
        sigrtg(1, i) = sigmad(1, 1)
        sigrtg(2, i) = sigmad(2, 2)
        sigrtg(3, i) = sigmad(3, 3)
        sigrtg(4, i) = sigmad(1, 2)
        if (ncmp .eq. 6) then
            sigrtg(5, i) = sigmad(1, 3)
            sigrtg(6, i) = sigmad(2, 3)
        end if

    end do
!
end subroutine
