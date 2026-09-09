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
subroutine pjrisp(model2, masp, corres, caraElem)
!
    use coorSyst_module, only: setOrieFields
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/celces.h"
#include "asterfort/cesexi.h"
#include "asterfort/cesvar.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jecrec.h"
#include "asterfort/jecreo.h"
#include "asterfort/jecroc.h"
#include "asterfort/jedema.h"
#include "asterfort/jeecra.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/setStructFields.h"
#include "asterfort/wkvect.h"
#include "jeveux.h"
!
    character(len=8) :: masp, model2, caraElem
    character(len=16) :: corres
!
! --------------------------------------------------------------------------------------------------
!
!           COMMANDE PROJ_CHAMP / METHODE='SOUS_POINT_RIGI'
!
!   créer un maillage dont les noeuds sont positionnés sur les sous-points de gauss
!   d'un modele (model2) pour la famille de points RIGI
!
! --------------------------------------------------------------------------------------------------
!
!   in :
!       model2    : modele "2".
!       masp    : nom du maillage des points de Gauss. La SD correspondante est vide.
!       corres  : nom de l'objet qui contient les données de correspondance des mailles.
!       caraElem    : nom du CARA_ELEM.
!
!   out :
!       masp    : SD contenant le maillage de POI1 correspondant aux points de Gauss
!       corres  : création de l'objet corres.PJEF_SP
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: option = 'COOR_ELGA'
    integer(kind=8), parameter :: nbFieldInMax = 100, nbFieldOut = 1
    character(len=8) :: lpain(nbFieldInMax), lpaout(nbFieldOut)
    character(len=19) :: lchin(nbFieldInMax), lchout(nbFieldOut)
!
    integer(kind=8) :: nbFieldIn
    integer(kind=8) :: ntgeo, ipo, ipg, nuno2
    integer(kind=8) :: ibid, nbnosp, nno2, ino2p
    integer(kind=8) :: j1, ipoi1
    integer(kind=8) :: nbCell, nbpt, nbsp, nbcmp
    integer(kind=8) :: ima, ipt, isp, icmp, iad, iadime
    integer(kind=8) :: jpo2, jcesd, jcesl, jcesv, iatypm
    integer(kind=8) :: nbpgmx, nbspmx
    character(len=8) :: mesh2
    character(len=19) :: chgeom, modelLigrel
    character(len=19), parameter :: chamg = '&&PJRISP.PGCOOR', ces = '&&PJRISP.PGCORS'
    character(len=24) :: coodsc
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()

! - Initializations
    lpain = ' '
    lchin = ' '
    lpaout = ' '
    lchout = ' '
    call dismoi('NOM_LIGREL', model2, 'MODELE', repk=modelLigrel)

!   récuperation du nom du maillage 2
    call dismoi('NOM_MAILLA', model2, 'MODELE', repk=mesh2)
    chgeom = mesh2//'.COORDO'

! - Add input fields
    lchin(1) = chgeom(1:19)
    lpain(1) = 'PGEOMER'
    nbFieldIn = 1

! - Add fields for structural elements
    call setStructFields(caraElem, nbFieldInMax, lchin, lpain, nbFieldIn)

! - Add fields for orientation
    call setOrieFields(nbFieldInMax, lpain, lchin, &
                       nbFieldIn, caraElem)

    call cesvar(caraElem, ' ', modelLigrel, chamg)
    call calcul('S', option, modelLigrel, &
                nbFieldIn, lchin, lpain, &
                nbFieldOut, chamg, 'PCOORPG', &
                'V', 'OUI')

!   chamg : 4 composantes X,Y,Z,W
!   Transformation en CHAM_ELEM_S
    call celces(chamg, 'V', ces)
!
    call jeveuo(ces//'.CESD', 'L', jcesd)
    call jeveuo(ces//'.CESL', 'L', jcesl)
    call jeveuo(ces//'.CESV', 'E', jcesv)
    nbCell = zi(jcesd-1+1)
!
! --------------------------------------------------------------------------------------------------
!   Calcul de nbnosp : nombre de noeuds (et de mailles) de masp
    nbnosp = 0
    nbpgmx = zi(jcesd-1+3)
    nbspmx = zi(jcesd-1+4)
!
!   On crée un tableau :
!      * la premiere valeur est le numero de la maille
!      * la deuxieme valeur est le numero du pg dans cette maille
!      * la troisieme valeur est le numero du sous-point
!   Dimension : (NBMA*NBPGMX*NBSPMX)*3 = (NB DE MAILLES * NB DE PG MAX  * NB DE SP MAX) * 3
!
    call wkvect(corres//'.PJEF_SP', 'V V I', nbCell*nbpgmx*nbspmx*3, jpo2)
!
    ipo = 1
    do ima = 1, nbCell
        nbpt = zi(jcesd-1+5+4*(ima-1)+1)
        nbsp = zi(jcesd-1+5+4*(ima-1)+2)
        if (nbsp .lt. 1) goto 100
        do ipg = 1, nbpt
            do isp = 1, nbsp
                zi(jpo2-1+ipo) = ima
                zi(jpo2-1+ipo+1) = ipg
                zi(jpo2-1+ipo+2) = isp
                ipo = ipo+3
            end do
        end do
        nbnosp = nbnosp+nbpt*nbsp
100     continue
    end do
!
! --------------------------------------------------------------------------------------------------
!   Création du .DIME du nouveau maillage
!   Il y a autant de mailles que de noeuds car toutes les mailles sont des poi1
!
    call wkvect(masp//'.DIME', 'V V I', 6, iadime)
    zi(iadime-1+1) = nbnosp
    zi(iadime-1+3) = nbnosp
    zi(iadime-1+6) = 3
!
! --------------------------------------------------------------------------------------------------
!   Création du .CONNEX et du .TYPMAIL du nouveau maillage
!
    call jecrec(masp//'.CONNEX', 'V V I', 'NU', 'CONTIG', 'VARIABLE', nbnosp)
    call jeecra(masp//'.CONNEX', 'LONT', nbnosp, ' ')
    call jeveuo(masp//'.CONNEX', 'E', ibid)
!
    call wkvect(masp//'.TYPMAIL', 'V V I', nbnosp, iatypm)
    call jenonu(jexnom('&CATA.TM.NOMTM', 'POI1'), ipoi1)
!
    nuno2 = 0
    do ima = 1, nbnosp
        zi(iatypm-1+ima) = ipoi1
        nno2 = 1
        call jecroc(jexnum(masp//'.CONNEX', ima))
        call jeecra(jexnum(masp//'.CONNEX', ima), 'LONMAX', nno2)
        nuno2 = nuno2+1
        zi(ibid-1+nuno2) = nuno2
    end do
!
! --------------------------------------------------------------------------------------------------
!
!   COORDO.VALE du nouveau maillage
    call wkvect(masp//'.COORDO    .VALE', 'V V R', 3*nbnosp, j1)
!
    ino2p = 0
    do ima = 1, nbCell
        nbpt = zi(jcesd-1+5+4*(ima-1)+1)
        nbsp = zi(jcesd-1+5+4*(ima-1)+2)
        nbcmp = zi(jcesd-1+5+4*(ima-1)+3)
        if (nbpt .eq. 0) goto 160
        ASSERT(nbcmp .ge. 3)
        do ipt = 1, nbpt
            do isp = 1, nbsp
                ino2p = ino2p+1
                do icmp = 1, 3
                    call cesexi('C', jcesd, jcesl, ima, ipt, &
                                isp, icmp, iad)
                    if (iad .gt. 0) then
                        zr(j1-1+3*(ino2p-1)+icmp) = zr(jcesv-1+iad)
                    end if
                end do
            end do
        end do
160     continue
    end do
    ASSERT(ino2p .eq. nbnosp)
!
! --------------------------------------------------------------------------------------------------
!   Création du .DESC du nouveau maillage
!
    coodsc = masp//'.COORDO    .DESC'
    call jenonu(jexnom('&CATA.GD.NOMGD', 'GEOM_R'), ntgeo)
    call jecreo(coodsc, 'V V I')
    call jeecra(coodsc, 'LONMAX', 3)
    call jeecra(coodsc, 'DOCU', cval='CHGO')
    call jeveuo(coodsc, 'E', iad)
    zi(iad) = ntgeo
    zi(iad+1) = -3
    zi(iad+2) = 14
!
    call detrsd('CHAM_ELEM', chamg)
    call detrsd('CHAM_ELEM_S', ces)
    call jedema()
end subroutine
