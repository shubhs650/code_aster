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
subroutine nttcmv(model, materCode, caraElem, listLoad, nume_dof, &
                  solver, timeMap, tpsthe, tpsnp1, reasvt, &
                  reasmt, creas, vtemp, vtempm, vec2nd, &
                  matass, maprec, cndirp, cnchci, cnchtp)
!
    implicit none
!
#include "jeveux.h"
#include "asterf_types.h"
#include "asterfort/asasve.h"
#include "asterfort/ascavc.h"
#include "asterfort/ascova.h"
#include "asterfort/asmatr.h"
#include "asterfort/dismoi.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/mecact.h"
#include "asterfort/medith.h"
#include "asterfort/mertth.h"
#include "asterfort/metnth.h"
#include "asterfort/preres.h"
#include "asterfort/vechth.h"
#include "asterfort/vedith.h"
!
    character(len=8), intent(in) :: model, caraElem
    character(len=24), intent(in) :: materCode, listLoad
    character(len=24), intent(in) :: nume_dof
    character(len=19), intent(in) :: solver
    character(len=24), intent(in) :: timeMap
    aster_logical :: reasvt, reasmt
    real(kind=8) :: tpsthe(6), tpsnp1
    character(len=1) :: creas
    character(len=19) :: maprec
    character(len=24) :: timeMapMove
    character(len=24) :: vtemp, vtempm, vec2nd
    character(len=24) :: matass, cndirp, cnchci, cnchtp
!
! --------------------------------------------------------------------------------------------------
!
! COMMANDE THER_MOBI_NLINE : ACTUALISATION
!   - DES VECTEURS CONTRIBUANT AU SECOND MEMBRE
!   - DE LA MATRICE ASSEMBLEE (EVENTUELLEMENT)
!
! --------------------------------------------------------------------------------------------------
!
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: ibid, k, iret, ierr, nbmat
    integer(kind=8) :: j2nd, lonch
    character(len=1) :: typres
    character(len=8) :: nomcmp(6)
    character(len=19) :: merigi
    character(len=24) :: ligrmo, mediri
    character(len=19) ::  tlimat(3)
    character(len=24) :: vediri, vechtp, vadirp, vachtp, timeMapMatr
    character(len=24), parameter :: metrnl = "&&METRNL"
    character(len=19) :: resuElem
    real(kind=8) :: time_curr
    character(len=24), pointer :: resuElemRelr(:) => null()
    real(kind=8), pointer :: chtp(:) => null()
    real(kind=8), pointer :: dirp(:) => null()
    character(len=24) :: loadNameJv, loadInfoJv, loadFuncJv
!
    data typres/'R'/
    data nomcmp/'INST    ', 'DELTAT  ', 'THETA   ', 'KHI     ', &
        'R       ', 'RHO     '/
    data mediri/'&&MEDIRI           .RELR'/
    data vediri/'&&VETDIR           .RELR'/
    data vechtp/'&&VETCHA           .RELR'/
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    vadirp = '&&VATDIR'
    vachtp = '&&VATCHA'
    timeMapMove = '&&NTTCMV.TIMEMO'
    timeMapMatr = '&&NTTCMV.TIMEMA'
    merigi = '&&METRIG'
    creas = ' '
    time_curr = tpsthe(1)

! - Access to datastructure of list of loads
    loadNameJv = listLoad(1:19)//'.LCHA'
    loadInfoJv = listLoad(1:19)//'.INFC'
    loadFuncJv = listLoad(1:19)//'.FCHA'
!
! ======================================================================
!         VECTEURS (CHARGEMENTS) CONTRIBUANT AU SECOND MEMBRE
! ======================================================================
!
    if (reasvt) then
!
! ----- Field for timeMap
        call dismoi('NOM_LIGREL', model, 'MODELE', repk=ligrmo)
        call mecact('V', timeMap, 'MODELE', ligrmo, 'INST_R', &
                    ncmp=6, lnomcmp=nomcmp, vr=tpsthe)
!
! ----- Field for shifted timeMap with 1-THETA
!
        tpsthe(3) = 1.d0
        call mecact('V', timeMapMove, 'MODELE', ligrmo, 'INST_R', &
                    ncmp=6, lnomcmp=nomcmp, vr=tpsthe)
        tpsthe(3) = -1.d0
        call mecact('V', timeMapMatr, 'MODELE', ligrmo, 'INST_R', &
                    ncmp=6, lnomcmp=nomcmp, vr=tpsthe)
        tpsthe(3) = 0.d0
!
! ----- TEMPERATURES IMPOSEES                                  ---> CNDIRP
!
        call vedith(model, loadNameJv, loadInfoJv, timeMap, vediri)
        call asasve(vediri, nume_dof, typres, 'D', vadirp)
        call ascova('D', vadirp, loadFuncJv, 'INST', time_curr, &
                    typres, cndirp)
        call jeveuo(cndirp(1:19)//'.VALE', 'E', vr=dirp)
!
! --- CHARGES CINEMATIQUES                                   ---> CNCHCI
!
        cnchci = ' '
        call ascavc(loadNameJv, loadInfoJv, loadFuncJv, nume_dof, tpsnp1, &
                    cnchci, l_hho_=ASTER_FALSE)
!
! --- CHARGEMENTS THERMIQUES                                 ---> CNCHTP
!            RQ : POUR LE CALCUL THERMIQUE, LES ARGUMENTS VTEMPP,
!                 VTEMPD ET THETA SONT INUTILISES.
!
        call vechth('MOVE', &
                    model, materCode, &
                    loadNameJv, loadInfoJv, &
                    time_curr, &
                    vechtp, &
                    timeMapZ_=timeMap, tempPrevZ_=vtemp, timeMoveZ_=timeMapMove)
        call asasve(vechtp, nume_dof, typres, 'D', vachtp)
        call ascova('D', vachtp, loadFuncJv, 'INST', time_curr, &
                    typres, cnchtp)
        call jeveuo(cnchtp(1:19)//'.VALE', 'E', vr=chtp)
        call jelira(cnchtp(1:19)//'.VALE', 'LONMAX', lonch)
!
! --- SECOND MEMBRE COMPLET                                  ---> VEC2ND
!
        call jeveuo(vec2nd(1:19)//'.VALE', 'E', j2nd)
        do k = 1, lonch
            zr(j2nd+k-1) = chtp(k)+dirp(k)
        end do
!
    end if
!
! ======================================================================
!              MATRICE ASSEMBLEE
! ======================================================================
!
    if (reasmt) then

! --- (RE)CALCUL DE LA MATRICE DES DIRICHLET POUR L'ASSEMBLER
        call medith('V', 'ZERO', model, listLoad, mediri)

! ----- Elementary matrix for transport (volumic and surfacic terms)
        creas = 'M'
        call mertth(model, loadNameJv, loadInfoJv, &
                    caraElem, materCode, &
                    timeMapMatr, timeMapMove, vtemp, vtempm, merigi)

! ----- Elementary matrix for boundary conditions
        call metnth(model, loadNameJv, loadInfoJv, &
                    caraElem, materCode, &
                    timeMap, vtempm, metrnl)
!
        nbmat = 0
        call jeveuo(merigi(1:19)//'.RELR', 'L', vk24=resuElemRelr)
        resuElem = resuElemRelr(1) (1:19)
        if (resuElem .ne. ' ') then
            nbmat = nbmat+1
            tlimat(nbmat) = merigi(1:19)
        end if
!
        call jeexin(metrnl(1:19)//'.RELR', iret)
        if (iret .gt. 0) then
            call jeveuo(metrnl(1:19)//'.RELR', 'L', vk24=resuElemRelr)
            resuElem = resuElemRelr(1) (1:19)
            if (resuElem .ne. ' ') then
                nbmat = nbmat+1
                tlimat(nbmat) = metrnl(1:19)
            end if
        end if
!
        call jeexin(mediri(1:8)//'           .RELR', iret)
        if (iret .gt. 0) then
            call jeveuo(mediri(1:8)//'           .RELR', 'L', vk24=resuElemRelr)
            resuElem = resuElemRelr(1) (1:19)
            if (resuElem .ne. ' ') then
                nbmat = nbmat+1
                tlimat(nbmat) = mediri(1:19)
            end if
        end if
!
! --- ASSEMBLAGE DE LA MATRICE
!
        call asmatr(nbmat, tlimat, ' ', nume_dof, &
                    listLoad, 'ZERO', 'V', 1, matass)
!
! --- DECOMPOSITION OU CALCUL DE LA MATRICE DE PRECONDITIONNEMENT
!
        call preres(solver, 'V', ierr, maprec, matass, &
                    ibid, -9999)
!
    end if
!
    call jedema()
end subroutine
