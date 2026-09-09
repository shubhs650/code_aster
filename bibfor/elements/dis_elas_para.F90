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
subroutine dis_elas_para(for_discret, nomphe)
!
    use te0047_type
    implicit none
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterfort/dikpkt.h"
#include "asterfort/dis_elas_para_klfl.h"
#include "asterfort/jevech.h"
#include "asterfort/tecach.h"
#include "asterfort/ut2mlg.h"
#include "asterfort/ut2vgl.h"
#include "asterfort/ut2vlg.h"
#include "asterfort/utpslg.h"
#include "asterfort/utpvgl.h"
#include "asterfort/utpvlg.h"
!
    type(te0047_dscr), intent(in) :: for_discret
    character(len=*), intent(in) :: nomphe
!
! --------------------------------------------------------------------------------------------------
!
! Ajout de la contribution d'un élément élastique parallèle à un élément discret
!
! --------------------------------------------------------------------------------------------------
! in :
!       for_discret : type dérivé de l'élément discret
!       nomphe      : nom du phénomène dans DEFI_MATERIAU (e.g. "DIS_CONTACT") pour lire (KP,KT)
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: imater, ii, jj, k, i, j
    integer(kind=8) :: iretlc, idepen
    integer(kind=8) :: imatr, icont, ifono
    integer(kind=8) :: nc, nn
    real(kind=8) :: kp, kt1, kt2
    real(kind=8) :: utotxyz(3), flp(3), ulp(12), dpe(12), fl(12), fg(12)
!
    real(kind=8) :: klvp(for_discret%nsym), kgvp(for_discret%nsym)
!
    aster_logical :: IsSymetrique
!
! --------------------------------------------------------------------------------------------------
!   Extraction des raideurs (KP,KT) à partir des caractéristiques du matériau
    call jevech('PMATERC', 'L', imater)
    call dikpkt(zi(imater), nomphe, kp, kt1, kt2)
!
    if ((abs(kp) .ge. r8prem()) .or. (abs(kt1) .ge. r8prem()) .or. (abs(kt2) .ge. r8prem())) then
        nn = for_discret%nno
        nc = for_discret%nc
        !   Déplacements totaux dans le repère local
        ! --- Déplacement
        ulp(:) = for_discret%ulm(:)+for_discret%dul(:)
        ! --- Déplacement d'entraînement
        dpe(:) = 0.0d0
        call tecach('ONO', 'PDEPENT', 'L', iretlc, iad=idepen)
        if (iretlc .eq. 0) then
            if (for_discret%ndim .eq. 3) then
                call utpvgl(nn, nc, for_discret%pgl, zr(idepen), dpe)
            else
                call ut2vgl(nn, nc, for_discret%pgl, zr(idepen), dpe)
            end if
        end if
        ! --- Déplacement total
        utotxyz(:) = 0.d0
        if (nn .eq. 2) then
            utotxyz(1) = ulp(1+nc)-ulp(1)+dpe(1+nc)-dpe(1)
            utotxyz(2) = ulp(2+nc)-ulp(2)+dpe(2+nc)-dpe(2)
            if (for_discret%ndim .eq. 3) then
                utotxyz(3) = ulp(3+nc)-ulp(3)+dpe(3+nc)-dpe(3)
            end if
        else if (nn .eq. 1) then
            utotxyz(1) = ulp(1)+dpe(1)
            utotxyz(2) = ulp(2)+dpe(2)
            if (for_discret%ndim .eq. 3) then
                utotxyz(3) = ulp(3)+dpe(3)
            end if
        end if

        kgvp = 0.d0
        klvp = 0.d0
        !   Matrice de raideur et vecteur force de l'élément élastique dans le repère local
        call dis_elas_para_klfl(for_discret, kp, kt1, kt2, utotxyz, klvp, flp)

        !   Mise à jour de la matrice tangente de l'élément élastique (dans le repère global)
        call hasSymmetricTangentMatrix(for_discret, IsSymetrique)
        if (for_discret%lMatr) then
            if (for_discret%ndim .eq. 3) then
                call utpslg(for_discret%nno, for_discret%nc, for_discret%pgl, klvp, kgvp)
            else
                call ut2mlg(for_discret%nno, for_discret%nc, for_discret%pgl, klvp, kgvp)
            end if
            if (IsSymetrique) then
                call jevech('PMATUUR', 'E', imatr)
                do ii = 1, for_discret%nsym
                    zr(imatr-1+ii) = zr(imatr-1+ii)+kgvp(ii)
                end do
            else
                call jevech('PMATUNS', 'E', imatr)
                k = 1
                do i = 1, for_discret%neq
                    do j = 1, i
                        ii = for_discret%neq*(i-1)
                        jj = for_discret%neq*(j-1)
                        zr(imatr-1+ii+j) = zr(imatr-1+ii+j)+kgvp(k)
                        if (i /= j) then
                            zr(imatr-1+jj+i) = zr(imatr-1+jj+i)+kgvp(k)
                        end if
                        k = k+1
                    end do
                end do
            end if

        end if

        !   Mise à jour des efforts généralisés
        if (for_discret%lSigm) then
            call jevech('PCONTPR', 'E', icont)
            if (nn .eq. 1) then
                zr(icont-1+1) = zr(icont-1+1)+flp(1)
                zr(icont-1+2) = zr(icont-1+2)+flp(2)
                if (for_discret%ndim .eq. 3) then
                    zr(icont-1+3) = zr(icont-1+3)+flp(3)
                end if
            else if (nn .eq. 2) then
                zr(icont-1+1) = zr(icont-1+1)+flp(1)
                zr(icont-1+1+nc) = zr(icont-1+1+nc)+flp(1)
                zr(icont-1+2) = zr(icont-1+2)+flp(2)
                zr(icont-1+2+nc) = zr(icont-1+2+nc)+flp(2)
                if (for_discret%ndim .eq. 3) then
                    zr(icont-1+3) = zr(icont-1+3)+flp(3)
                    zr(icont-1+3+nc) = zr(icont-1+3+nc)+flp(3)
                end if
            end if
        end if

        !   Mise à jour des forces nodales
        if (for_discret%lVect) then
            call jevech('PVECTUR', 'E', ifono)
            fl(1:12) = 0.d0
            if (nn .eq. 1) then
                fl(1) = flp(1)
                fl(2) = flp(2)
                if (for_discret%ndim .eq. 3) then
                    fl(3) = flp(3)
                end if
                ! Attention aux signes des efforts sur le premier noeud pour MECA_DIS_*_L
            else if (nn .eq. 2) then
                fl(1) = -flp(1)
                fl(1+nc) = flp(1)
                fl(2) = -flp(2)
                fl(2+nc) = flp(2)
                if (for_discret%ndim .eq. 3) then
                    fl(3) = -flp(3)
                    fl(3+nc) = flp(3)
                end if
            end if
            ! Passage dans le repère global
            fg(1:12) = 0.d0
            if (for_discret%ndim .eq. 3) then
                call utpvlg(nn, nc, for_discret%pgl, fl, fg)
            else
                call ut2vlg(nn, nc, for_discret%pgl, fl, fg)
            end if
            ! Mise à jour des forces nodales dans le repère global
            if (nn .eq. 1) then
                zr(ifono-1+1) = zr(ifono-1+1)+fg(1)
                zr(ifono-1+2) = zr(ifono-1+2)+fg(2)
                if (for_discret%ndim .eq. 3) then
                    zr(ifono-1+3) = zr(ifono-1+3)+fg(3)
                end if
            else if (nn .eq. 2) then
                zr(ifono-1+1) = zr(ifono-1+1)+fg(1)
                zr(ifono-1+1+nc) = zr(ifono-1+1+nc)+fg(1+nc)
                zr(ifono-1+2) = zr(ifono-1+2)+fg(2)
                zr(ifono-1+2+nc) = zr(ifono-1+2+nc)+fg(2+nc)
                if (for_discret%ndim .eq. 3) then
                    zr(ifono-1+3) = zr(ifono-1+3)+fg(3)
                    zr(ifono-1+3+nc) = zr(ifono-1+3+nc)+fg(3+nc)
                end if
            end if
        end if
    end if
!
end subroutine
