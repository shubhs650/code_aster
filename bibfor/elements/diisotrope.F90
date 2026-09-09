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
subroutine diisotrope(for_discret, iret)
!
! --------------------------------------------------------------------------------------------------
!
!        COMPORTEMENT ISOTROPE
!
! --------------------------------------------------------------------------------------------------
!
! IN    for_discret : voir l'appel
! OUT   iret        : code retour
!
! --------------------------------------------------------------------------------------------------
!
    use te0047_type
    implicit none
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/diraidklv.h"
#include "asterfort/infdis.h"
#include "asterfort/jevech.h"
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"
#include "asterfort/diisotrope_uni.h"
#include "asterfort/diisotrope_bid.h"
#include "blas/dcopy.h"
!
    type(te0047_dscr), intent(in) :: for_discret
    integer(kind=8), intent(out) :: iret
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) ::  ivarim, jdc, irep, iadzi, iazk24
    integer(kind=8) ::   ipi, imate, imater, jmat, nbmat
    integer(kind=8) :: icontm, neq, kk, nb_fonc
    real(kind=8) :: r8bid, klv(78), raide(6)
    character(len=8) :: k8bid
    character(len=24) :: messak(6)
!   pour le matériau
    character(len=16) :: materiau
!   pour la loi de comportement
!!!   si on change nbfct, penser à changer la dimension de ldcfct dans les .h
!   Équations du système
!
!
! --------------------------------------------------------------------------------------------------
!   Paramètres associés au matériau codé
    blas_int :: b_incx, b_incy, b_n
! --------------------------------------------------------------------------------------------------
!
    iret = 0
    neq = for_discret%nno*for_discret%nc
! Récupération du matériau
    call jevech('PMATERC', 'L', imater)
! Variables internes a t- : Force  Up  Puiss  tangente
    call jevech('PVARIMR', 'L', ivarim)
! Effort à t-
    call jevech('PCONTMR', 'L', icontm)
! récupération des caractéristiques élastique
    call jevech('PCADISK', 'L', jdc)
    call infdis('REPK', irep, r8bid, k8bid)
! Seulement en 3D
    if (for_discret%nomte(1:10) .ne. 'MECA_DIS_T') then
        messak(1) = for_discret%nomte
        messak(2) = 'NON_LINEAR'
        messak(3) = for_discret%type_comp
        messak(4) = for_discret%rela_comp
        call tecael(iadzi, iazk24)
        messak(5) = zk24(iazk24-1+3)
        call utmess('F', 'DISCRETS_22', nk=5, valk=messak)
    end if
! Seulement en repère local : irep = 2
    if (irep .ne. 2) then
        messak(1) = for_discret%nomte
        messak(2) = 'NON_LINEAR'
        messak(3) = for_discret%type_comp
        messak(4) = for_discret%rela_comp
        call tecael(iadzi, iazk24)
        messak(5) = zk24(iazk24-1+3)
        call utmess('F', 'DISCRETS_5', nk=5, valk=messak)
    end if
! les caractéristiques sont toujours dans le repère local. on fait seulement une copie
    b_n = to_blas_int(for_discret%nbt)
    b_incx = to_blas_int(1)
    b_incy = to_blas_int(1)
    call dcopy(b_n, zr(jdc), b_incx, klv, b_incy)

! Récupère les termes diagonaux de la matrice de raideur
    call diraidklv(for_discret%nomte, raide, klv)
    !
! Adresse de la SD mater
    jmat = zi(imater)
! Nombre de matériau sur la maille : 1 seul autorisé
    nbmat = zi(jmat)
    ASSERT(nbmat .eq. 1)
! Adresse du matériau codé
    imate = jmat+zi(jmat+nbmat+1)

! Recherche du matériau dans la SD compor
    materiau = 'DIS_ECRO_TRAC'
    ipi = 0
    do kk = 1, zi(imate+1)
        if (zk32(zi(imate)+kk-1) (1:16) .eq. materiau) then
            ipi = zi(imate+2+kk-1)
            goto 10
        end if
    end do
    messak(1) = for_discret%nomte
    messak(2) = 'NON_LINEAR'
    messak(3) = for_discret%type_comp
    messak(4) = materiau
    call tecael(iadzi, iazk24)
    messak(5) = zk24(iazk24-1+3)
    call utmess('F', 'DISCRETS_7', nk=5, valk=messak)
10  continue
    !

    nb_fonc = zi(ipi+2)
    ASSERT(nb_fonc == 1 .or. nb_fonc == 2)
    ! --- Appel du bloc extrait vers diisotrope_uni ---
    if (nb_fonc == 1) then
        ! si une seule direction (TAN ou X) est renseignée

        call diisotrope_uni(for_discret, iret, ipi, jmat, ivarim, icontm, klv, raide)

    elseif (nb_fonc == 2) then
        ! si deux directions (TAN et X) sont renseignées
        call diisotrope_bid(for_discret, iret, ipi, jmat, ivarim, icontm, klv, raide)

    end if
end subroutine
