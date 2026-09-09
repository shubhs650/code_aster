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

subroutine exlim5(motfaz, motcleZ, toutZ, nomsd, modelz, ligrel)
    implicit none
#include "jeveux.h"
#include "asterc/getexm.h"
#include "asterc/getfac.h"
#include "asterc/getres.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/exlim1.h"
#include "asterfort/getvtx.h"
#include "asterfort/gnomsd.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeexin.h"
#include "asterfort/jeveuo.h"
#include "asterfort/juveca.h"
#include "asterfort/reliem.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    character(len=*) :: motfaz, motcleZ, toutZ, nomsd, modelz, ligrel
! but  :  scruter les mots cle tout/group_ma/maille pour creer
!         un ligrel "reduit" a partir du ligrel du modele modelz
!
! in  : modelz : nom du modele
!
! out/jxout   : ligrel  : ligrel reduit
!     attention :
!          - le nom de ligrel est toujours "out"
!          - parfois on rend ligrel=ligrel(modele) :
!             - alors on ne tient donc pas compte de 'base'
!             - il ne faut pas le detruire !
!          - parfois on en cree un nouveau sur la base 'base'
!             - le nom du ligrel est obtenu par gnomsd
!  -----------------------------------------------------------------
!
    integer(kind=8) :: n1, jma, nbma, iocc, nbocc, nbmatot, iexi, ima
    character(len=1) :: base
    character(len=8) :: modele, noma
    character(len=16) :: motfac, motcle(2), typmcl(2), oper, k16b, tout
    character(len=19) :: ligrmo
    character(len=24) :: lismai, noojb, lismaiT
    integer(kind=8), pointer :: v_lismai(:) => null()
    integer(kind=8), pointer :: v_lismai_tot(:) => null()
    parameter(base='G')
!  -----------------------------------------------------------------
!
    motfac = motfaz
    modele = modelz
    tout = toutZ

    if (modele .eq. ' ') then
        call utmess('F', 'UTILITAI8_10')
    end if
!
    call dismoi('NOM_LIGREL', modele, 'MODELE', repk=ligrmo)
    call dismoi('NOM_MAILLA', modele, 'MODELE', repk=noma)
    lismai = '&&EXLIM5.LISTE_MAILLES'
    lismaiT = '&&EXLIM5.LISTE_MAILLES_T'
!
!     --  SI ON DOIT TOUT PRENDRE , LIGREL = LIGRMO
!     ------------------------------------------------------
    nbmatot = 0
    call getfac(motfac, nbocc)
    do iocc = 1, nbocc
        if (motfac .ne. ' ') then
            if (getexm(motfac, tout) .eq. 1) then
                call getvtx(motfac, tout, iocc=iocc, nbval=0, nbret=n1)
                if (n1 .ne. 0) then
                    call jedetr(lismai)
                    call jedetr(lismaiT)
                    goto 9998
                end if
            end if
        else
            if (getexm(' ', tout) .eq. 1) then
                call getvtx(' ', tout, nbval=0, nbret=n1)
                if (n1 .ne. 0) then
                    call jedetr(lismai)
                    call jedetr(lismaiT)
                    goto 9998
                end if
            end if
        end if
!
        motcle(1) = motcleZ
        motcle(2) = 'MAILLE'
        typmcl(1) = 'GROUP_MA'
        typmcl(2) = 'MAILLE'
!
!     --- CREATION ET AFFECTATION DU VECTEUR DE K8 DE NOM LISMAI
!         CONTENANT LES NOMS DES MAILLES FORMANT LE LIGREL A CREER
!         --------------------------------------------------------
        call reliem(modele, noma, 'NU_MAILLE', motfac, iocc, &
                    2, motcle(1), typmcl(1), lismai, nbma)
!
!         -- SI LES MOTS CLES GROUP_MA ET MAILLE N'ONT PAS ETE UTILISES:
        if (nbma .eq. 0) then
            call jedetr(lismai)
            cycle
        end if
        call jeexin(lismaiT, iexi)
        if (iexi .eq. 0) then
            call wkvect(lismaiT, 'V V I', nbma, vi=v_lismai_tot)
        else
            call juveca(lismaiT, nbmatot+nbma)
            call jeveuo(lismaiT, 'L', vi=v_lismai_tot)
        end if
        call jeveuo(lismai, 'L', vi=v_lismai)
        do ima = 1, nbma
            v_lismai_tot(nbmatot+ima) = v_lismai(ima)
        end do
        nbmatot = nbmatot+nbma
        call jedetr(lismai)
    end do
!
!
!
! --- CREATION DU LIGREL
!     ---------------------------------
    if (nbmatot .ne. 0) then
        call getres(k16b, k16b, oper)
        if (oper .ne. 'IMPR_RESU') then
            noojb = nomsd//'.LIGR000000.LIEL'
            call gnomsd(' ', noojb, 14, 19)
        else
            ASSERT(.false.)
        end if
        ligrel = noojb(1:19)
        ASSERT(ligrel(1:8) .ne. ' ')
        call jeveuo(lismaiT, 'L', jma)
        call exlim1(zi(jma), nbma, modele, base, ligrel)
        call jedetr(lismaiT)
    end if
    goto 999
!
!
9998 continue
    ligrel = ligrmo
!
999 continue
!
!
end subroutine
