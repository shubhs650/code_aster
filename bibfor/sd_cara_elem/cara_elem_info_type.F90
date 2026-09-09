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

module cara_elem_info_type
!
!
!   Définition des "user_type" utilisés que par AFFE_CARA_ELEM.
!
! --------------------------------------------------------------------------------------------------
!
!   cara_elem_info  : Variable globale définit pour tout AFFE_CARA_ELEM
!       nomu            : nom du concept en sortie de la commande.  getres(nomu , x , x )
!       concept         : nom du concept résultat.                  getres( x , concept , x )
!       commande        : nom de la commande.                       getres( x , x , commande )
!       modele          : nom du modèle.
!       maillage        : nom du maillage.
!       modmail         : nom de la sd des mailles du maillages  : modele//'.MAILLE'
!       jmodmail        : pointeur sur la SD modmail !!! Doit être fait avec jeveut !!!
!       nbnoeu          : nombre de noeud du maillage.
!       nbmail          : nombre de maille du maillage.
!       GroupeMaxOccur   : nombre d'entité sous les mots clefs simple GROUP_MA
!       dimmod          : dimension topologique du modèle.      dismoi('DIM_GEOM' sur 'MODELE')
!
!       Pour certains mot clef : RIGI_PARASOL, ...
!           MailleMaxOccur  : Le nombre de maille traitée dans les occurrences
!           NoeudMaxMaille  : Le nombre max de noeud par maille
!
!       ivr         : Pour faire des vérifications de cohérences et des impressions
!           ivr(1)=1    : vérification MAILLES
!           ivr(2)      : libre
!           ivr(3)=niv  : niveau d'impression
!           ivr(4)=ifm  : unité d'impression
!
! --------------------------------------------------------------------------------------------------
!
#include "asterf_types.h"
!
    implicit none
!
    type cara_elem_info
        character(len=8)  :: nomu
        character(len=16) :: concept
        character(len=8)  :: modele
        character(len=8)  :: maillage
        character(len=24) :: modmail
        aster_logical     :: IsParaMesh
        integer(kind=8)   :: jmodmail
        integer(kind=8)   :: nbmail
        integer(kind=8)   :: nbnoeu
        integer(kind=8)   :: dimmod
        integer(kind=8)   :: GroupeMaxOccur
        integer(kind=8)   :: NoeudMaxMaille, MailleMaxOccur
        integer(kind=8)   :: ivr(4)
        logical           :: VerifMaille
    end type cara_elem_info
!
end module
