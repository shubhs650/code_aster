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

subroutine chpchd(chin, type, celmod, prol0, base, &
                  chou, model_)
    implicit none
#include "jeveux.h"
#include "asterf_types.h"
#include "asterc/getres.h"
#include "asterfort/assert.h"
#include "asterfort/carces.h"
#include "asterfort/celces.h"
#include "asterfort/celfpg.h"
#include "asterfort/cescar.h"
#include "asterfort/cescel.h"
#include "asterfort/cesces.h"
#include "asterfort/cescns.h"
#include "asterfort/cgocns.h"
#include "asterfort/cnocns.h"
#include "asterfort/cnsces.h"
#include "asterfort/cnscno.h"
#include "asterfort/global_numbering_compute.h"
#include "asterfort/global_numbering_communicate.h"
#include "asterfort/detrsd.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedetr.h"
#include "asterfort/manopg.h"
#include "asterfort/utmess.h"
#include "asterfort/xcesrd.h"
#include "asterfort/xnpgxx.h"

    character(len=*) :: chin, chou, base, celmod, type
    character(len=8), optional, intent(in) :: model_

! -----------------------------------------------------------------
!  BUT : CHANGER LE SUPPORT GEOMETRIQUE D'UN CHAMP
! -----------------------------------------------------------------
!
! CHIN IN/JXIN  K19 : NOM DU CHAMP A CHANGER
!        TYPE AUTORISE POUR CHIN :
!           NOEU, CART, ELGA, ELNO, ELEM, CESE
!
! CHOU IN/JXOUT K19 : NOM DU CHAMP RESULTAT
! BASE IN       K1  : BASE DE CREATION DE CHOU : /'G' / 'V'
! TYPE IN       K19 : TYPE DE SUPPORT GEOMETRIQUE VOULU POUR CHOU
!                     /'NOEU' /'CART' /'ELNO' /ELGA' /'ELEM'
!
! ARGUMENTS UTILISES SI TYPE=ELNO/ELGA/ELEM :
!   PROL0 IN   K3  :
!        /'OUI' : LE CHAM_ELEM CHOU EST PROLONGE
!         PAR DES VALEURS NULLES LA OU IL N'EST PAS DEFINI.
!        /'NON' : ERREUR <F> SI IL EXISTE DES
!         DES VALEURS DE CHOU QUI NE SONT PAS AFFECTEES DANS CHIN
!   CELMOD IN/JXIN  K19 : NOM D'UN CHAM_ELEM "MODELE" SI TYPE='EL..'
!                         POUR XFEM UNIQUEMENT
!                       : NOM DU MAILLAGE SI TYPE='GEOM.'
!
!  LES CAS TRAITES AUJOURD'HUI SONT :
!
!         /'GEOM->NOEU'   : CHAM_GEOM -> CHAM_NO
!         /'GEOM->ELNO'   : CHAM_GEOM -> CHAM_ELEM/ELNO

!         /'NOEU->ELNO'   : CHAM_NO -> CHAM_ELEM/ELNO
!         /'NOEU->ELGA'   : CHAM_NO -> CHAM_ELEM/ELGA
!
!         /'CART->ELEM'   : CARTE   -> CHAM_ELEM/ELEM
!         /'CART->ELGA'   : CARTE   -> CHAM_ELEM/ELGA
!         /'CART->ELNO'   : CARTE   -> CHAM_ELEM/ELNO
!         /'CART->NOEU'   : CARTE   -> CHAM_NO
!
!         /'ELGA->NOEU'   : CHAM_ELEM/ELGA    -> CHAM_NO
!         /'ELGA->ELNO'   : CHAM_ELEM/ELGA    -> CHAM_ELEM/ELNO

!         /'ELNO->NOEU'   : CHAM_ELEM/ELNO    -> CHAM_NO
!         /'ELNO->ELGA'   : CHAM_ELEM/ELNO    -> CHAM_ELEM/ELGA
!
!         /'ELEM->ELNO'   : CHAM_ELEM/ELEM    -> CHAM_ELEM/ELNO
!         /'ELEM->ELGA'   : CHAM_ELEM/ELEM    -> CHAM_ELEM/ELGA
!
!         /'CESE->ELNO'   : CHAM_ELEM_S/ELEM  -> CHAM_ELEM/ELNO
!         /'CESE->ELGA'   : CHAM_ELEM_S/ELEM  -> CHAM_ELEM/ELGA
!         /'CESE->ELEM'   : CHAM_ELEM_S/ELEM  -> CHAM_ELEM/ELEM
! -----------------------------------------------------------------
!
    integer(kind=8) :: ib, iret, nncp, ibid
    character(len=3) :: prol0, exixfm
    character(len=8) :: ma, ma2, tychi, nomgd, param, model
    character(len=16) :: cas, option, nomcmd, kbid
    character(len=19) :: cesmod, ces1, cns1, mnoga, ligrel, ces2, chsnpg, nume_equa
    character(len=24) :: valk(4)
    aster_logical :: bool
!
!     ------------------------------------------------------------------
    mnoga = '&&CHPCHD.MANOGA'
    cesmod = '&&CHPCHD.CESMOD'
    ces1 = '&&CHPCHD.CES1'
    ces2 = '&&CHPCHD.CES2'
    cns1 = '&&CHPCHD.CNS1'

    model = ' '
    if (present(model_)) then
        model = model_
    end if
!
!
! 1- CALCUL DE:
!      MA    : MAILLAGE ASSOCIE A CHIN
!      TYCHI : TYPE DU CHAMP CHIN (CART/NOEU/ELNO/ELGA/CESE)
!      NOMGD : NOM DE LA GRANDEUR ASSOCIEE A CHIN
! ------------------------------------------------------------------
!
    call dismoi('TYPE_CHAMP', chin, 'CHAMP', repk=tychi)
    call dismoi('NOM_GD', chin, 'CHAMP', repk=nomgd)
    bool = tychi .eq. 'NOEU' .or. tychi .eq. 'CART' .or. tychi .eq. 'ELNO' .or. tychi .eq. 'ELGA' &
           .or. tychi .eq. 'ELEM' .or. tychi .eq. 'CESE' .or. tychi .eq. 'GEOM'
    ASSERT(bool)
!
!
! 2.  -- SI TYPE = 'EL..' : ON CREE UN CHAM_ELEM_S "MODELE" : CESMOD
!         LIGREL: NOM DU LIGREL ASSOCIE A CHOU
! ---------------------------------------------------------------
    if (type(1:2) .eq. 'EL') then
        ASSERT(celmod .ne. ' ')
        call dismoi('NOM_LIGREL', celmod, 'CHAM_ELEM', repk=ligrel)
        call dismoi('NOM_OPTION', celmod, 'CHAM_ELEM', repk=option)
        call dismoi('NOM_PARAM', celmod, 'CHAM_ELEM', repk=param)
        call dismoi('NOM_MAILLA', ligrel, 'LIGREL', repk=ma2)
        call dismoi('NOM_MAILLA', chin, 'CHAMP', repk=ma)
        if (ma .ne. ma2) then
            call utmess('F', 'CALCULEL4_59')
        end if
        call celces(celmod, 'V', cesmod)
    end if
!
!
! 3.  -- CALCUL DE CAS :
! ---------------------------------------
!     SONT TRAITES AUJOURD'HUI :
!
!         /'GEOM->NOEU'   : CHAM_GEOM -> CHAM_NO
!         /'GEOM->ELNO'   : CHAM_GEOM -> CHAM_ELEM/ELNO
!
!         /'NOEU->ELNO'   : CHAM_NO -> CHAM_ELEM/ELNO
!         /'NOEU->ELGA'   : CHAM_NO -> CHAM_ELEM/ELGA
!
!         /'CART->ELEM'   : CARTE   -> CHAM_ELEM/ELEM
!         /'CART->ELGA'   : CARTE   -> CHAM_ELEM/ELGA
!         /'CART->ELNO'   : CARTE   -> CHAM_ELEM/ELNO
!         /'CART->NOEU'   : CARTE   -> CHAM_NO
!
!         /'ELGA->NOEU'   : CHAM_ELEM/ELGA    -> CHAM_NO
!         /'ELGA->ELNO'   : CHAM_ELEM/ELGA    -> CHAM_ELEM/ELNO

!         /'ELNO->NOEU'   : CHAM_ELEM/ELNO    -> CHAM_NO
!         /'ELNO->ELGA'   : CHAM_ELEM/ELNO    -> CHAM_ELEM/ELGA
!
!         /'ELEM->ELNO'   : CHAM_ELEM/ELEM    -> CHAM_ELEM/ELNO
!         /'ELEM->ELGA'   : CHAM_ELEM/ELEM    -> CHAM_ELEM/ELGA
!
!         /'CESE->ELNO'   : CHAM_ELEM_S/ELEM  -> CHAM_ELEM/ELNO
!         /'CESE->ELGA'   : CHAM_ELEM_S/ELEM  -> CHAM_ELEM/ELGA
!         /'CESE->ELEM'   : CHAM_ELEM_S/ELEM  -> CHAM_ELEM/ELEM
!
    cas = ' '
    cas(1:4) = tychi
    cas(5:6) = '->'
    cas(7:10) = type
!
!
! 4.  TRAITEMENT DES DIFFERENTS CAS DE FIGURE :
! ----------------------------------------------
    nncp = 0
    if (cas .eq. 'NOEU->ELGA') then
!     ----------------------------------
        call manopg(model, ligrel, option, param, mnoga)
!
        call cnocns(chin, 'V', cns1)
        call cnsces(cns1, 'ELGA', cesmod, mnoga, 'V', &
                    ces1)
        call detrsd('CHAM_NO_S', cns1)
        call detrsd('CHAM_ELEM_S', mnoga)
!
!       construction du CHAM_ELEM_S conteant le nombre de points
!       de Gauss reellement utilises par chaque element, dans le
!       cas d'un champs ELGA, base sur la famille "XFEM"
        chsnpg = '&&CHPCHD.CHSNPG'
        call xnpgxx(model, ligrel, option, param, chsnpg, exixfm)
!
        if (exixfm .eq. 'OUI') then
!            si le champ ELGA s'appuie sur la famille "XFEM", on
!            desaffecte toutes les composantes associes aux points
!            de Gauss inutilises
            call xcesrd(ces1, chsnpg)
        end if
!
        call detrsd('CHAM_ELEM_S', chsnpg)
!
        call cescel(ces1, ligrel, option, param, prol0, &
                    nncp, base, chou, 'F', ibid)
        call detrsd('CHAM_ELEM_S', ces1)
!
!
    else if (cas .eq. 'ELNO->ELGA') then
!     ----------------------------------
        call manopg(model, ligrel, option, param, mnoga)
!
        call celces(chin, 'V', ces1)
        call cesces(ces1, 'ELGA', cesmod, mnoga, ' ', &
                    'V', ces2)
        call detrsd('CHAM_ELEM_S', ces1)
        call detrsd('CHAM_ELEM_S', mnoga)
!
        call cescel(ces2, ligrel, option, param, prol0, &
                    nncp, base, chou, 'F', ibid)
        call detrsd('CHAM_ELEM_S', ces2)
!
!
    else if (cas .eq. 'NOEU->ELNO') then
!     ----------------------------------------------------------------
!
        call cnocns(chin, 'V', cns1)
        call cnsces(cns1, 'ELNO', cesmod, ' ', 'V', &
                    ces1)
        call detrsd('CHAM_NO_S', cns1)
!
        call cescel(ces1, ligrel, option, param, prol0, &
                    nncp, base, chou, 'F', ibid)
        call detrsd('CHAM_ELEM_S', ces1)
!
    else if (cas(1:4) .eq. 'GEOM') then
!     ----------------------------------------------------------------
!
        call cgocns(chin, 'V', cns1, celmod)
        if (cas(7:10) == "NOEU") then
            call cnscno(cns1, ' ', 'NON', base, chou, 'F', ibid, lprofconst=ASTER_FALSE)
            call dismoi('NUME_EQUA', chou, 'CHAM_NO', repk=nume_equa)
            ! create numbering
            call global_numbering_compute(nume_equa)
            ! communicate numbering
            call global_numbering_communicate(nume_equa)
        elseif (cas(7:10) == "ELNO") then
            call cnsces(cns1, 'ELNO', cesmod, ' ', 'V', ces1)
!
            call cescel(ces1, ligrel, option, param, prol0, &
                        nncp, base, chou, 'F', ibid)
            call detrsd('CHAM_ELEM_S', ces1)
        else
            ASSERT(ASTER_FALSE)
        end if
!
        call detrsd('CHAM_NO_S', cns1)
!
!
    elseif ((cas .eq. 'ELNO->NOEU') .or. (cas .eq. 'ELGA->NOEU') .or. &
            (cas .eq. 'CART->NOEU')) then
!     ----------------------------------------------------------------
!
        if (cas(1:4) .eq. 'ELNO') then
            call celces(chin, 'V', ces1)
!
        else if (cas(1:4) .eq. 'ELGA') then
            call celces(chin, 'V', ces1)
            call celfpg(chin, '&&CHPCHD.CELFPG', iret)
!
        else if (cas(1:4) .eq. 'CART') then
            call carces(chin, 'ELNO', ' ', 'V', ces1, &
                        'A', iret)
!
        else
            ASSERT(.false.)
        end if
!
        call cescns(ces1, '&&CHPCHD.CELFPG', 'V', cns1, ' ', &
                    iret)
        call cnscno(cns1, ' ', 'NON', base, chou, &
                    'F', ibid)
!
        call detrsd('CHAM_NO_S', cns1)
        call detrsd('CHAM_ELEM_S', ces1)
        call jedetr('&&CHPCHD.CELFPG')
!
!
    else if (cas(1:8) .eq. 'CART->EL' .or. cas(1:8) .eq. 'ELEM->EL') then
!     ----------------------------------------------------------------
        ASSERT(ligrel .ne. ' ')
        if (cas(1:4) .eq. 'CART') then
            call carces(chin, cas(7:10), cesmod, 'V', ces1, &
                        'A', ib)
        else
            ASSERT(cas(1:4) .eq. 'ELEM')
            call celces(chin, 'V', ces2)
            call cesces(ces2, cas(7:10), cesmod, ' ', ' ', 'V', ces1)
            call detrsd('CHAM_ELEM_S', ces2)
        end if
!
        if (type .eq. 'ELGA') then
!           construction du CHAM_ELEM_S conteant le nombre de points
!           de Gauss reellement utilises par chaque element, dans le
!           cas d'un champs ELGA, base sur la famille "XFEM"
            chsnpg = '&&CHPCHD.CHSNPG'
            call xnpgxx(model, ligrel, option, param, chsnpg, exixfm)
!
            if (exixfm .eq. 'OUI') then
!                si le champ ELGA s'appuie sur la famille "XFEM", on
!                desaffecte toutes les composantes associes aux points
!                de Gauss inutilises
                call xcesrd(ces1, chsnpg)
            end if
!
            call detrsd('CHAM_ELEM_S', chsnpg)
        end if
!
        call cescel(ces1, ligrel, option, param, prol0, &
                    nncp, base, chou, 'F', ibid)
        call detrsd('CHAM_ELEM_S', ces1)
!
        if (nncp .ne. 0) then
            call getres(kbid, kbid, nomcmd)
            if (nomcmd .eq. 'CREA_CHAMP') then
                valk(1) = chou(1:8)
                valk(2) = option
                valk(3) = param
                call utmess('A', 'CALCULEL6_77', nk=3, valk=valk)
            end if
        end if
!
!
    else if (cas(1:8) .eq. 'CESE->EL') then
!     ----------------------------------------------------------------
        ASSERT(ligrel .ne. ' ')
        call cesces(chin, cas(7:10), cesmod, ' ', ' ', 'V', ces1)
        call cescel(ces1, ligrel, option, param, prol0, &
                    nncp, base, chou, 'F', ibid)
        call detrsd('CHAM_ELEM_S', ces1)
!
        if (nncp .ne. 0) then
            call getres(kbid, kbid, nomcmd)
            if (nomcmd .eq. 'CREA_CHAMP') then
                valk(1) = chou(1:8)
                valk(2) = option
                valk(3) = param
                call utmess('A', 'CALCULEL6_77', nk=3, valk=valk)
            end if
        end if
!
!
    else if (cas .eq. 'ELGA->ELNO') then
!     ----------------------------------------------------------------
!
        call celces(chin, 'V', ces1)
        call celfpg(chin, '&&CHPCHD.CELFPG', iret)
        ASSERT(iret .eq. 0)
        call cesces(ces1, 'ELNO', cesmod, ' ', '&&CHPCHD.CELFPG', &
                    'V', ces2)
        call cescel(ces2, ligrel, option, param, prol0, &
                    nncp, base, chou, 'F', ibid)
!
        call detrsd('CHAM_ELEM_S', ces1)
        call detrsd('CHAM_ELEM_S', ces2)
        call jedetr('&&CHPCHD.CELFPG')
!
!
    else if (cas .eq. 'NOEU->ELEM') then
!     ----------------------------------------------------------------
!
        call cnocns(chin, 'V', cns1)
        call cnsces(cns1, 'ELEM', cesmod, ' ', 'V', &
                    ces1)
        call detrsd('CHAM_NO_S', cns1)
!
        call cescel(ces1, ligrel, option, param, prol0, &
                    nncp, base, chou, 'F', ibid)
        call detrsd('CHAM_ELEM_S', ces1)
!
!
    else if (cas .eq. 'ELEM->NOEU') then
        call celces(chin, 'V', ces1)
        call cescns(ces1, ' ', 'V', cns1, 'F', ibid)
        call detrsd('CHAM_ELEM_S', ces1)
        call cnscno(cns1, ' ', 'NON', base, chou, &
                    'F', ibid)
        call detrsd('CHAM_NO_S', cns1)

    else if (cas .eq. 'ELEM->CART') then
        call celces(chin, 'V', ces1)
        call cescar(ces1, chou, "G")
        call detrsd('CHAM_ELEM_S', ces1)

    else
!       CAS NON ENCORE PROGRAMME
        call utmess('F', 'CALCULEL_5', sk=cas)
    end if
!
!
!
!
!     -- MENAGE :
!     ------------
    if (type(1:2) .eq. 'EL') call detrsd('CHAM_ELEM_S', cesmod)
!
end subroutine
