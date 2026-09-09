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

subroutine mdchdl(lnoue2, iliai, ddlcho, ier, l_rotaz)
    implicit none
#include "asterf_types.h"
#include "asterfort/posddl.h"
#include "asterfort/utmess.h"
#include "asterfort/nlget.h"

    aster_logical, intent(in) :: lnoue2
    integer(kind=8), intent(in) :: iliai
    integer(kind=8), intent(out) :: ddlcho(*), ier
    aster_logical, intent(in), optional :: l_rotaz
!
!     ROUTINE APPELEE PAR MDCHOC,
!     TRAITEMENT DES DDL
!
! IN  : NBNLI  : DIMENSION DES TABLEAUX (NBCHOC+NBSISM+NBFLAM)
! IN  : NOECHO : DEFINITION DES NOEUDS DE CHOC
! IN  : LNOUE2 : CHOC DEFINIT ENTRE 2 NOEUDS
! IN  : ILIAI  : NUMERO DE LA LIAISON TRAITEE
! OUT : DDLCHO : TABLEAU DES NUMEROTATIONS DES NOEUDS DE CHOC
! OUT : IER    : NIVEAU D'ERREUR
!     ------------------------------------------------------------------
    integer(kind=8) :: nunoe, nuddl, nbcomp
    character(len=8) :: nume1, noeu1, nume2, noeu2, sd_nl
    character(len=24) :: valk(2)
    aster_logical :: l_rota
!     ------------------------------------------------------------------
!
    if (present(l_rotaz)) then
        l_rota = l_rotaz
    else
        l_rota = ASTER_FALSE
    end if

    sd_nl = '&&OP29NL'
    call nlget(sd_nl, _NUMDDL_1, iocc=iliai, kscal=nume1)
    call nlget(sd_nl, _NO1_NAME, iocc=iliai, kscal=noeu1)
!
    call posddl('NUME_DDL', nume1, noeu1, 'DX', nunoe, &
                nuddl)
    if (nunoe .eq. 0) then
        ier = ier+1
        valk(1) = noeu1
        call nlget(sd_nl, _MESH_1, iocc=iliai, kscal=valk(2))
        call utmess('E', 'ALGORITH5_27', nk=2, valk=valk)
    end if
    if (nuddl .eq. 0) then
        ier = ier+1
        call utmess('E', 'ALGORITH5_28', sk=noeu1)
    end if
    ddlcho(1) = nuddl
!
    call posddl('NUME_DDL', nume1, noeu1, 'DY', nunoe, &
                nuddl)
    if (nuddl .eq. 0) then
        ier = ier+1
        call utmess('E', 'ALGORITH5_29', sk=noeu1)
    end if
    ddlcho(2) = nuddl
!
    call posddl('NUME_DDL', nume1, noeu1, 'DZ', nunoe, &
                nuddl)
    if (nuddl .eq. 0) then
        ier = ier+1
        call utmess('E', 'ALGORITH5_30', sk=noeu1)
    end if
    ddlcho(3) = nuddl
    nbcomp = 3

    if (l_rota) then
        call posddl('NUME_DDL', nume1, noeu1, 'DRX', nunoe, &
                    nuddl)
        if (nuddl .eq. 0) then
            ier = ier+1
            call utmess('E', 'ALGORITH5_28', sk=noeu1)
        end if
        ddlcho(4) = nuddl
        !
        call posddl('NUME_DDL', nume1, noeu1, 'DRY', nunoe, &
                    nuddl)
        if (nuddl .eq. 0) then
            ier = ier+1
            call utmess('E', 'ALGORITH5_29', sk=noeu1)
        end if
        ddlcho(5) = nuddl
        !
        call posddl('NUME_DDL', nume1, noeu1, 'DRZ', nunoe, &
                    nuddl)
        if (nuddl .eq. 0) then
            ier = ier+1
            call utmess('E', 'ALGORITH5_30', sk=noeu1)
        end if
        ddlcho(6) = nuddl
        nbcomp = 6
    end if
!
    if (lnoue2) then

        call nlget(sd_nl, _NUMDDL_2, iocc=iliai, kscal=nume2)
        call nlget(sd_nl, _NO2_NAME, iocc=iliai, kscal=noeu2)

        call posddl('NUME_DDL', nume2, noeu2, 'DX', nunoe, &
                    nuddl)
        if (nunoe .eq. 0) then
            ier = ier+1
            valk(1) = noeu2
            call nlget(sd_nl, _MESH_2, iocc=iliai, kscal=valk(2))
            call utmess('E', 'ALGORITH5_27', nk=2, valk=valk)
        end if
        if (nuddl .eq. 0) then
            ier = ier+1
            call utmess('E', 'ALGORITH5_28', sk=noeu2)
        end if
        ddlcho(nbcomp+1) = nuddl
!
        call posddl('NUME_DDL', nume2, noeu2, 'DY', nunoe, &
                    nuddl)
        if (nuddl .eq. 0) then
            ier = ier+1
            call utmess('E', 'ALGORITH5_29', sk=noeu2)
        end if
        ddlcho(nbcomp+2) = nuddl
!
        call posddl('NUME_DDL', nume2, noeu2, 'DZ', nunoe, &
                    nuddl)
        if (nuddl .eq. 0) then
            ier = ier+1
            call utmess('E', 'ALGORITH5_30', sk=noeu2)
        end if
        ddlcho(nbcomp+3) = nuddl

        if (l_rota) then
            call posddl('NUME_DDL', nume2, noeu2, 'DRX', nunoe, &
                        nuddl)
            if (nuddl .eq. 0) then
                ier = ier+1
                call utmess('E', 'ALGORITH5_28', sk=noeu2)
            end if
            ddlcho(nbcomp+4) = nuddl
            !
            call posddl('NUME_DDL', nume2, noeu2, 'DRY', nunoe, &
                        nuddl)
            if (nuddl .eq. 0) then
                ier = ier+1
                call utmess('E', 'ALGORITH5_29', sk=noeu2)
            end if
            ddlcho(nbcomp+5) = nuddl
            !
            call posddl('NUME_DDL', nume2, noeu2, 'DRZ', nunoe, &
                        nuddl)
            if (nuddl .eq. 0) then
                ier = ier+1
                call utmess('E', 'ALGORITH5_30', sk=noeu2)
            end if
            ddlcho(nbcomp+6) = nuddl
        end if
    else
        ddlcho(nbcomp+1:2*nbcomp) = ddlcho(1:nbcomp)
    end if
!
end subroutine
