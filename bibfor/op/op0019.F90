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

subroutine op0019()
!
! --------------------------------------------------------------------------------------------------
!
!                O P E R A T E U R    AFFE_CARA_ELEM
!
! --------------------------------------------------------------------------------------------------
!
    use cara_elem_parameter_module
    use cara_elem_info_type
    use cara_elem_carte_type
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterfort/aceaco.h"
#include "asterfort/aceadi.h"
#include "asterfort/aceagb.h"
#include "asterfort/aceama.h"
#include "asterfort/aceamb.h"
#include "asterfort/aceamr.h"
#include "asterfort/aceaor.h"
#include "asterfort/aceapc.h"
#include "asterfort/aceapf.h"
#include "asterfort/aceapo.h"
#include "asterfort/acearg.h"
#include "asterfort/acearm.h"
#include "asterfort/acearp.h"
#include "asterfort/ace_crea_carte.h"
#include "asterfort/ace_affe_barre.h"
#include "asterfort/ace_affe_cable.h"
#include "asterfort/ace_masse_repartie.h"
#include "asterfort/ace_verif_affe.h"
#include "asterfort/ace_get_node_reparti.h"
#include "asterfort/as_allocate.h"
#include "asterfort/as_deallocate.h"
#include "asterfort/checkCaraElem.h"
#include "asterfort/coqucf.h"
#include "asterfort/detrsd.h"
#include "asterfort/detrsd_vide.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvid.h"
#include "asterfort/getvtx.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/isParallelMesh.h"
#include "asterfort/jedema.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jeveut.h"
#include "asterfort/jeexin.h"
#include "asterfort/pmfd00.h"
#include "asterfort/tecart.h"
#include "asterfort/utmess.h"
#include "asterfort/verif_affe.h"
#include "asterfort/verima.h"
#include "asterfort/wkvect.h"
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8)   :: nb_type_elem(ACE_NB_ELEMENT)
    integer(kind=8)   :: nbocc(ACE_NB_MCLEF)
! --------------------------------------------------------------------------------------------------
!   Pour les cartes :
    type(cara_elem_carte)   :: info_carte(ACE_NB_CARTE)
!   Infomation sur le concept : maillage, modele, nb noeuds, nb mailles, ...
    type(cara_elem_info) :: info_concept
!
! --------------------------------------------------------------------------------------------------
    character(len=8) :: mclef_type
!
    integer(kind=8) :: ivr(4), iret, jadr, ii, ireponse
    integer(kind=8) :: nbelemdisc, nbelemstrx, nbelemori
    integer(kind=8) :: ixma
    integer(kind=8) :: GroupeMaxOccur, ifm, niv, nbvm
    integer(kind=8) :: nbmail, nbnoeu
    integer(kind=8) :: npoutr, ncable, nbarre, nbdisc
    integer(kind=8) :: iclf, ioc, icle, ng
    integer(kind=8) :: jdnm
    aster_logical :: locagb, locamb
    character(len=3)  :: verif
    character(len=8)  :: nomu, nomo, noma
    character(len=16) :: concep, cmd, mclef
    character(len=19) :: ligrel
    character(len=24) :: modnom, mlgnma, mlgnno
!
! --------------------------------------------------------------------------------------------------
    integer(kind=8), pointer    :: affe_mail(:) => null()
    character(len=24), pointer  :: grp_lmax(:) => null()
    integer(kind=8), pointer    :: grp_nbma(:) => null()
! --------------------------------------------------------------------------------------------------
    call jemarq()
!   CALL ONERRF('ABORT', K16BID, IRET)
    iret = 0
! --------------------------------------------------------------------------------------------------
!   Récupération des arguments de la commande
    call getres(nomu, concep, cmd)
! --------------------------------------------------------------------------------------------------
!   Modèle
    call getvid(' ', 'MODELE', scal=nomo, nbret=nbvm)
!   Enregistre le nom du modèle dans la SD de AFFE_CARA_ELEM
    call wkvect(nomu//'.MODELE', 'G V K8', 1, jadr)
    zk8(jadr) = nomo
!   Construction des noms jeveux du concept modèle
    modnom = nomo//'.MODELE    .LGRF'
!   Récupération du nom du maillage associé
    call jeveuo(modnom, 'L', jdnm)
    noma = zk8(jdnm)
!   Construction des noms jeveux du concept maillage associé
    mlgnma = noma//'.TYPMAIL'
    mlgnno = noma//'.COORDO    .VALE'
!   Nombre de mailles du maillage
    call jelira(mlgnma, 'LONMAX', nbmail)
!   Nombre de noeuds du maillage
    call jelira(mlgnno, 'LONMAX', nbnoeu)
    nbnoeu = nbnoeu/3
    info_concept%nbmail = nbmail
    info_concept%nbnoeu = nbnoeu
!   Récupération de la dimension géométrique du modèle
    call dismoi('DIM_GEOM', nomo, 'MODELE', repi=ireponse)
    info_concept%dimmod = ireponse
    if (ireponse .ge. 100) then
        ireponse = ireponse-100
        info_concept%dimmod = 1
    end if
    if (ireponse .ge. 20) then
        ireponse = ireponse-20
        info_concept%dimmod = 2
    end if
    if (ireponse .eq. 3) info_concept%dimmod = 3
!
!   Mémorisation des informations pour ne plus le refaire
    info_concept%nomu = nomu
    info_concept%concept = concep
    info_concept%modele = nomo
    info_concept%maillage = noma
    info_concept%IsParaMesh = isParallelMesh(noma)
!
    call dismoi('NOM_LIGREL', nomo, 'MODELE', repk=ligrel)
    info_concept%modmail = ligrel//'.TYFE'
    call jeexin(info_concept%modmail, ixma)
    if (ixma .ne. 0) then
!       C'est un jeveut : c'est volontaire
        call jeveut(info_concept%modmail, 'L', info_concept%jmodmail)
    end if
!
    info_concept%NoeudMaxMaille = 0
    info_concept%MailleMaxOccur = 0
!
! --------------------------------------------------------------------------------------------------
    ivr(:) = 0
!   Pour faire des vérifications de cohérence d'affectation
    info_concept%VerifMaille = .true.
    verif = 'OUI'
    call getvtx(' ', 'VERIF', scal=verif, nbret=ng)
!   ivr(1)      : vérification MAILLES
!                   0 : On fait la vérification, par défaut
!                   1 : On ne fait pas la vérifications
!   ivr(2)      : libre
!   ivr(3)=niv  : niveau d'impression
!   ivr(4)=ifm  : unité d'impression
    if ((ng .gt. 0) .and. (verif .eq. 'NON')) then
        ivr(1) = 1
        info_concept%VerifMaille = .false.
    end if
!   Récupération du niveau d'impression
    call infmaj()
    call infniv(ifm, niv)
    ivr(3) = niv
    ivr(4) = ifm
    info_concept%ivr(:) = ivr(:)
!
! --------------------------------------------------------------------------------------------------
!   Occurence des mots clefs facteur
    do ii = 1, ACE_NB_MCLEF
        call getfac(ACE_MCLEF(ii), nbocc(ii))
    end do
!
! --------------------------------------------------------------------------------------------------
!   Initialisation des éléments pouvant être affetés :
!       elem_supp%[ catanom | acenum | catanum | aceind]
    call ACE_Init_elem_affe()
! --------------------------------------------------------------------------------------------------
!   Vérification de l'existence des GROUP_MA*
!       Après cette vérification il n'est plus nécessaire d'utiliser les routines
!           En sequentiel plus besoin de : verima , getvem(getvtx+verima)
!           En //Mesh les GROUP_MA peuvent ne pas exister sur le processeur
!               Il faut donc faire un verima(group_ma,nbgrp)
!                   En sortie group_ma,nbgrp sont modifiés pour ne contenir que les groupes
!                   présents sur le processeur actuel
!   Comptage des GROUP_MA* ... pour ne pas faire des ALLOCATE dans les boucles.
    GroupeMaxOccur = 10
    do iclf = 1, ACE_NB_MCLEF
        do ioc = 1, nbocc(iclf)
            doicle0: do icle = 1, ACE_NB_GRMA_MA
                ii = MCLEF_GRP_MA(icle+(iclf-1)*ACE_NB_GRMA_MA)
                if (ii .eq. ACE_NOTHING) cycle doicle0
                mclef = ACE_GRMA_MA(ii)
                call getvtx(ACE_MCLEF(iclf), mclef, iocc=ioc, nbval=0, nbret=ng)
                GroupeMaxOccur = max(GroupeMaxOccur, -ng)
            end do doicle0
        end do
    end do
    if (GroupeMaxOccur .lt. 10) then
        call utmess('F', 'AFFECARAELEM_2')
    end if
!   Vérifications
    AS_ALLOCATE(vk24=grp_lmax, size=GroupeMaxOccur)
    AS_ALLOCATE(vi=grp_nbma, size=GroupeMaxOccur)
!
!   Si vérification
    if (info_concept%VerifMaille) then
        do iclf = 1, ACE_NB_MCLEF
            do ioc = 1, nbocc(iclf)
                doicle1: do icle = 1, ACE_NB_GRMA_MA
                    ii = MCLEF_GRP_MA(icle+(iclf-1)*ACE_NB_GRMA_MA)
                    if (ii .eq. ACE_NOTHING) cycle doicle1
                    mclef = ACE_GRMA_MA(ii)
                    mclef_type = ACE_GRMA_TY(ii)
                    call getvtx(ACE_MCLEF(iclf), mclef, iocc=ioc, nbval=GroupeMaxOccur, &
                                vect=grp_lmax, nbret=ng)
                    call verima(noma, grp_lmax, ng, mclef_type)
                end do doicle1
            end do
        end do
    end if
!
    info_concept%GroupeMaxOccur = GroupeMaxOccur
!
! --------------------------------------------------------------------------------------------------
!   Pour mémoriser les mailles affectées.
!       Si affe_mail(i) =  elem_supp%catanum    la maille est dans le modèle mais pas affectée
!                       < 0                     la maille est dans le modèle et est affectée
!                       =  0                    la maille n'est pas affectée
!                                               ou n'est pas connue de affe_cara_elem
    AS_ALLOCATE(vi=affe_mail, size=nbmail)
! --------------------------------------------------------------------------------------------------
!   Compteur d'éléments et vérification cohérence des affectations
!       Prise en compte du //
    affe_mail(:) = 0
    call ace_verif_affe(info_concept, nbocc, nb_type_elem, affe_mail)
!   Nombre d'éléments sur le processeur actuel (cas du //)
    npoutr = nb_type_elem(ACE_NU_POUTRE)
    ncable = nb_type_elem(ACE_NU_CABLE)
    nbarre = nb_type_elem(ACE_NU_BARRE)
    nbdisc = nb_type_elem(ACE_NU_DISCRET)
!
    locagb = nb_type_elem(ACE_NU_GRILLE) .ne. 0
    locamb = nb_type_elem(ACE_NU_MEMBRANE) .ne. 0
!
! --------------------------------------------------------------------------------------------------
!   Création des cartes utilisées par la suite
    call ace_crea_carte(info_concept, info_carte)
!
! --------------------------------------------------------------------------------------------------
!   Affectation des éléments CABLE
    call ace_affe_cable(nbocc(ACE_CABLE), info_concept, info_carte, grp_lmax, grp_nbma, affe_mail)
!
! --------------------------------------------------------------------------------------------------
!   Affectation des éléments BARRE
    call ace_affe_barre(nbocc(ACE_BARRE), info_concept, info_carte, grp_lmax, grp_nbma, affe_mail)
!
! --------------------------------------------------------------------------------------------------
!   Traitement des masses réparties
    call ace_masse_repartie(nbocc(ACE_MASS_REP), info_concept, info_carte, grp_lmax, grp_nbma, &
                            affe_mail, nbdisc)
!
! --------------------------------------------------------------------------------------------------
!   Prédétermination du nombre de noeuds des discrets affectés par des raideurs réparties
    call ace_get_node_reparti(nbocc(ACE_RIGI_PARASOL), ACE_MCLEF(ACE_RIGI_PARASOL), &
                              info_concept, grp_lmax)
    call ace_get_node_reparti(nbocc(ACE_MASS_AJOU), ACE_MCLEF(ACE_MASS_AJOU), &
                              info_concept, grp_lmax)
    call ace_get_node_reparti(nbocc(ACE_RIGI_MISS_3D), ACE_MCLEF(ACE_RIGI_MISS_3D), &
                              info_concept, grp_lmax)
!
! --------------------------------------------------------------------------------------------------
!   Affectation des orientations aux éléments 'orientable'
!       Si il existe un élément orientable mais pas orienté avec le mot clef ORIENTATION
!       la carte doit être créée et être nulle
    nbelemori = nbocc(ACE_POUTRE)+nbocc(ACE_DISCRET)+nbocc(ACE_DISCRET_2D)+ &
                nbocc(ACE_BARRE)+nbocc(ACE_RIGI_PARASOL)+nbocc(ACE_CABLE)
    if (nbelemori .ne. 0) then
        call aceaor(nbocc(ACE_ORIENTATION), info_concept)
    end if
! --------------------------------------------------------------------------------------------------
!   Affectation des caractéristiques aux éléments poutres
    if (nbocc(ACE_POUTRE) .ne. 0) then
        if (npoutr .gt. 0) then
            call aceapo(nbocc(ACE_POUTRE), npoutr, info_concept, info_carte, affe_mail)
        end if
    end if
! --------------------------------------------------------------------------------------------------
!   AFFECTATION DES EPAISSEURS/COURBURES/ANGLES AUX ELEMENTS COQUES
    if (nbocc(ACE_COQUE) .ne. 0) then
        call aceaco(nomu, noma, GroupeMaxOccur, locagb, locamb, nbocc(ACE_COQUE))
    end if
! --------------------------------------------------------------------------------------------------
!   Affectation des coefficients de correction pour les COUDES
    if (nbocc(ACE_POUTRE) .ne. 0) then
        call aceapc(nomu, noma, GroupeMaxOccur, nbocc(ACE_POUTRE))
    end if
!
! --------------------------------------------------------------------------------------------------
!   AFFECTATION DES REPERES AUX ELEMENTS THERMIQUES ET MECANIQUES
    if (nbocc(ACE_MASSIF) .ne. 0) then
        call aceama(nomu, noma, GroupeMaxOccur, nbocc(ACE_MASSIF))
    end if
! --------------------------------------------------------------------------------------------------
!   AFFECTATION DES REPERES AUX ELEMENTS POUTRE_FLUI
    if (nbocc(ACE_POUTRE_FLUI) .ne. 0) then
        call aceapf(nomu, noma, GroupeMaxOccur, nbocc(ACE_POUTRE_FLUI))
    end if
! --------------------------------------------------------------------------------------------------
!   AFFECTATION DES MATRICES AUX RAIDEURS REPARTIES
    if (nbocc(ACE_RIGI_PARASOL) .ne. 0) then
        call acearp(nbocc(ACE_RIGI_PARASOL), info_concept, info_carte, grp_lmax, affe_mail)
    end if
    if (nbocc(ACE_RIGI_GRILLE) .ne. 0) then
        call acearg(nbocc(ACE_RIGI_GRILLE), info_concept, info_carte, affe_mail)
    end if
! --------------------------------------------------------------------------------------------------
!   AFFECTATION DES MATRICES AUX ELEMENTS DISCRETS
    if (nbocc(ACE_DISCRET) .ne. 0) then
        call aceadi(nbocc(ACE_DISCRET), info_concept, info_carte, ACE_MCLEF(ACE_DISCRET))
    else if (nbocc(ACE_DISCRET_2D) .ne. 0) then
        call aceadi(nbocc(ACE_DISCRET_2D), info_concept, info_carte, ACE_MCLEF(ACE_DISCRET_2D))
    end if
! --------------------------------------------------------------------------------------------------
!   AFFECTATION DES CARACTERISTIQUES POUR L'ELEMENT "GRILLE"
    if (nbocc(ACE_GRILLE) .ne. 0) then
        call aceagb(nomu, noma, locamb, nbocc(ACE_GRILLE))
    end if
! --------------------------------------------------------------------------------------------------
!   AFFECTATION DES MATRICES AUX RAIDEURS MISS
    if (nbocc(ACE_RIGI_MISS_3D) .ne. 0) then
        call acearm(nbocc(ACE_RIGI_MISS_3D), info_concept, info_carte)
    end if
! --------------------------------------------------------------------------------------------------
!   AFFECTATION DES CARACTERISTIQUES POUR L'ELEMENT "MEMBRANE"
    if (nbocc(ACE_MEMBRANE) .ne. 0) then
        call aceamb(nomu, noma, GroupeMaxOccur, nbocc(ACE_MEMBRANE))
    end if
! --------------------------------------------------------------------------------------------------
!   AFFECTATION DES MATRICES AUX MASSES REPARTIES
    if (nbocc(ACE_MASS_AJOU) .ne. 0) then
        call aceamr(nbocc(ACE_MASS_AJOU), info_concept, info_carte, grp_lmax)
    end if
!
!   POUR LES COQUES, GRILLES IL PEUT EXISTER UNE CARTE FONCTION
!   IL FAUT L'EVALUER ET METTRE LE RESULTAT DANS LA CARTE DES REELS
    if ((nbocc(ACE_COQUE) .ne. 0) .or. (nbocc(ACE_GRILLE) .ne. 0)) then
        call coqucf(nomu)
    end if
!
! --------------------------------------------------------------------------------------------------
!   TRAITEMENT DES MOTS CLES
!           MULTIFIBRE  /  GEOM_FIBRE
!           COQUE       /  COQUE_NCOU
!           GRILLE      /  COQUE_NCOU
!           MEMBRANE    /  COQUE_NCOU
!           POUTRE      /  TUYAU_NCOU  TUYAU_NSEC
    call pmfd00()
! --------------------------------------------------------------------------------------------------
!   Appel de l'option de verification VERI_CARA_ELEM
    if (info_concept%VerifMaille) then
        call checkCaraElem(nomo, nomu)
    end if
!
! --------------------------------------------------------------------------------------------------
!   Destruction sélective des CARTES  !! attention au cas //Mesh
!       Si elles n'ont pas lieu d'être
!       Si elles sont inutilisées
    nbelemdisc = nbdisc
    nbelemstrx = npoutr+nbarre+ncable
!
    do ii = 1, ACE_NB_CARTE
        if (ACE_CARTE(3+(ii-1)*ACE_NB_CARTE_CMP) .eq. 'DISCRET') then
            if (nbelemdisc .eq. 0) then
                call detrsd('CHAMP', info_carte(ii)%nom_carte)
            end if
        else if (ACE_CARTE(3+(ii-1)*ACE_NB_CARTE_CMP) .eq. 'STRX_LINEI') then
            if (nbelemstrx .eq. 0) then
                call detrsd('CHAMP', info_carte(ii)%nom_carte)
            else
!               Compactage de la carte, si elle est utilisée sinon destruction
                if (info_carte(ii)%utilise) then
                    call tecart(info_carte(ii)%nom_carte)
                else
                    call detrsd('CHAMP', info_carte(ii)%nom_carte)
                end if
            end if
        else
!           Si la carte n'a pas été utilisée
            if (.not. info_carte(ii)%utilise) then
                call detrsd('CHAMP', info_carte(ii)%nom_carte)
            end if
        end if
    end do
!   Certaines cartes peuvent être vides : il faut les detruire.
    do ii = 1, ACE_NB_CARTE
        call detrsd_vide('CARTE', info_carte(ii)%nom_carte)
    end do
!
! --------------------------------------------------------------------------------------------------
!   Audit assignments
    if (info_concept%VerifMaille) then
        call verif_affe(modele=nomo, sd=nomu)
    end if
!
    AS_DEALLOCATE(vi=affe_mail)
    AS_DEALLOCATE(vk24=grp_lmax)
    AS_DEALLOCATE(vi=grp_nbma)
!
    call jedema()
end subroutine
