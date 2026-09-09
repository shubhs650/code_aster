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
subroutine te0415(optioz, nomtz)
    implicit none
#include "jeveux.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/tecach.h"
#include "asterfort/utmess.h"
#include "asterfort/Behaviour_type.h"
!
    character(len=*) :: optioz, nomtz
    character(len=16) :: option, nomte
!     ----------------------------------------------------------------
!     CALCUL DES OPTIONS DES ELEMENTS DE COQUE 3D
!     OPTIONS : VARI_ELNO
!          -----------------------------------------------------------
!
!-----------------------------------------------------------------------
    integer(kind=8) :: i, i1, ic, ichg
    integer(kind=8) :: ino, inp, iret
    integer(kind=8) :: j, j1
    integer(kind=8) :: jvari, k1, k2, l, lgpg
    integer(kind=8) :: lzi, lzr, nbcou, nbvari, nep, np1
    integer(kind=8) :: np2, np3, np4, npge, npo, npp
    integer(kind=8) :: nso
    real(kind=8) :: s
!-----------------------------------------------------------------------
    character(len=16), pointer :: compor(:) => null()
    parameter(npge=3)
    integer(kind=8) :: icou, jvGano, jnbspi
    integer(kind=8) :: nb2, npgsn, jtab(7)
!
    option = optioz
    nomte = nomtz
    call jevete('&INEL.'//nomte(1:8)//'.DESI', ' ', lzi)
    nb2 = zi(lzi-1+2)
    npgsn = zi(lzi-1+4)
    call jevete('&INEL.'//nomte(1:8)//'.DESR', ' ', lzr)
    if (nomte .eq. 'MEC3QU9H') then
        nso = 4
    else if (nomte .eq. 'MEC3TR7H') then
        nso = 3
    end if
!
    if (option .eq. 'VARI_ELNO') then
!
        call jevech('PVARIGR', 'L', ichg)
        call jevech('PCOMPOR', 'L', vk16=compor)
        read (compor(NVAR), '(I16)') nbvari
        call tecach('OOO', 'PVARIGR', 'L', iret, nval=7, &
                    itab=jtab)
        lgpg = max(jtab(6), 1)*jtab(7)
        call jevech('PNBSP_I', 'L', jnbspi)
        nbcou = zi(jnbspi-1+1)
        if (nbcou .le. 0) then
            call utmess('F', 'ELEMENTS_12')
        end if
!
!
!
! -- RECUPERATION DES VARIABLES INTERNES
! -- NBVARI = NOMBRES DE VARIABLES INTERNES
! -- STOCKAGE DANS PVARIGR : PAR POINT DE GAUSS DU PREMIER AU DERNIER
!
        call jevete('&INEL.'//nomte//'.B', ' ', jvGano)
!
!-- EXTRAPOLATION AUX NOEUDS SOMMETS (3 OU 4)
!
        call jevech('PVARINR', 'E', jvari)
!
        do icou = 1, nbcou
            do ic = 1, nbvari
                do i = 1, npge*nso
                    l = npge*npgsn*(i-1)
                    s = 0.d0
                    do j = 1, npge*npgsn
! -- DETERMINATION DU PT DE GAUSS A PARTIR DE LA POSITION JJ
                        do k1 = 1, npgsn
                            do k2 = 1, npge
                                j1 = (k1-1)*npge+k2
                                if (j1 .eq. j) then
                                    inp = k1
                                    nep = k2-1
                                end if
                            end do
                        end do
                        npp = (inp-1)*lgpg
                        npp = npp+ic+nbvari*((icou-1)*npge+nep)
! -- ZR(ICHG-1+NPP) = VARI(IC,JJ)
!                JJ = (ICOU-1)*NPGE*NPGSN + J
                        s = s+zr(jvGano-1+l+j)*zr(ichg-1+npp)
                    end do
! -- DETERMINATION DU NOEUD SOMMET A PARTIR DE LA POSITION II
                    do k1 = 1, nso
                        do k2 = 1, npge
                            i1 = (k1-1)*npge+k2
                            if (i1 .eq. i) then
                                ino = k1
                                nep = k2-1
                            end if
                        end do
                    end do
                    npo = (ino-1)*lgpg
                    npo = npo+ic+nbvari*((icou-1)*npge+nep)
                    zr(jvari-1+npo) = s
                end do
            end do
        end do
!
! -- CREATION DU CHAMP DE VARIABLES INTERNES POUR LES POINTS
! -- MILIEUX ET LE CENTRE
! -- STOCKAGE DANS PVARINR : PAR NOEUD DU PREMIER AU DERNIER
!
        if (nomte .eq. 'MEC3QU9H') then
            do ic = 5, nb2
                npo = (ic-1)*lgpg
                if (ic .eq. 5) then
                    np1 = 0
                    np2 = lgpg
                else if (ic .eq. 6) then
                    np1 = lgpg
                    np2 = 2*lgpg
                else if (ic .eq. 7) then
                    np1 = 2*lgpg
                    np2 = 3*lgpg
                else if (ic .eq. 8) then
                    np1 = 3*lgpg
                    np2 = 0
                else if (ic .eq. 9) then
                    np1 = 0
                    np2 = lgpg
                    np3 = 2*lgpg
                    np4 = 3*lgpg
                end if
                if (ic .ne. 9) then
                    do i = 1, lgpg
                        zr(jvari-1+npo+i) = zr(jvari-1+np1+i)+zr(jvari-1+np2+i)
                        zr(jvari-1+npo+i) = zr(jvari-1+npo+i)/2.d0
                    end do
                else
                    do i = 1, lgpg
                        zr(jvari-1+npo+i) = zr(jvari-1+np1+i)+zr(jvari-1+np2+i)
                        zr(jvari-1+npo+i) = zr(jvari-1+npo+i)+zr(jvari-1+np3+i)
                        zr(jvari-1+npo+i) = zr(jvari-1+npo+i)+zr(jvari-1+np4+i)
                        zr(jvari-1+npo+i) = zr(jvari-1+npo+i)/4.d0
                    end do
                end if
            end do
        else if (nomte .eq. 'MEC3TR7H') then
            do ic = 4, nb2
                npo = (ic-1)*lgpg
                if (ic .eq. 4) then
                    np1 = 0
                    np2 = lgpg
                else if (ic .eq. 5) then
                    np1 = lgpg
                    np2 = 2*lgpg
                else if (ic .eq. 6) then
                    np1 = 2*lgpg
                    np2 = 0
                else if (ic .eq. 7) then
                    np1 = 0
                    np2 = lgpg
                    np3 = 2*lgpg
                end if
                if (ic .ne. 7) then
                    do i = 1, lgpg
                        zr(jvari-1+npo+i) = zr(jvari-1+np1+i)+zr(jvari-1+np2+i)
                        zr(jvari-1+npo+i) = zr(jvari-1+npo+i)/2.d0
                    end do
                else
                    do i = 1, lgpg
                        zr(jvari-1+npo+i) = zr(jvari-1+np1+i)+zr(jvari-1+np2+i)
                        zr(jvari-1+npo+i) = zr(jvari-1+npo+i)+zr(jvari-1+np3+i)
                        zr(jvari-1+npo+i) = zr(jvari-1+npo+i)/3.d0
                    end do
                end if
            end do
        end if
!
! ------------------------------------------------------------
!
    end if
end subroutine
