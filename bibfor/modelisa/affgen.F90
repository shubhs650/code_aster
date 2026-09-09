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

subroutine affgen(tmp, nom, nel, napcis, foncis)
!
    use cara_elem_parameter_module
!
    implicit none
#include "jeveux.h"
#include "asterc/r8pi.h"
#include "asterfort/assert.h"
#include "asterfort/fointe.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/utmess.h"
!
    character(len=19) :: napcis, foncis
    character(len=24) :: tmp, nom
!       AFFECTATION DES CARACTÉRISTIQUES GÉNÉRALES CALCULÉES
!       A PARTIR DES DONNÉES GÉOMÉTRIQUES (RECTANGLE ET CERCLE)
!       RQ : NTEL(1) = NUMÉRO DU TYPE ÉLÉMENT MECA_POU_D_T
!            NTEL(2) = NUMÉRO DU TYPE ÉLÉMENT MECA_POU_D_E
!            NTEL(4) = NUMÉRO DU TYPE ÉLÉMENT MEFS_POU_D_T
!            NTEL(5) = NUMÉRO DU TYPE ÉLÉMENT MECA_POU_D_TG
!            NTEL(14) = NUMÉRO DU TYPE ÉLÉMENT MECA_S8P2
!     ------------------------------------------------------------------
    real(kind=8) :: eps, pi, alpha, beta, ccis
    real(kind=8) :: hy, hz, epy, epz, hyi, hzi
    real(kind=8) :: a, b, a4, b4, b3
    real(kind=8) :: ct, cd, jx
    real(kind=8) :: re, ri, e
    real(kind=8) :: valpay(2), valpaz(2), valpaf
    character(len=24) :: nompa(2), nompaf
!-----------------------------------------------------------------------
    integer(kind=8) :: ier, igen, igen2, igeoc, igeor, isec, ii
    integer(kind=8) :: jdge, nel, idx1poutre, idx2poutre
    real(kind=8) :: aireint, ay, az
!-----------------------------------------------------------------------
    data eps/1.d-3/
!     ------------------------------------------------------------------
!
    call jemarq()
    pi = r8pi()
!
!   Index de départ des poutres dans la SD cara_elem
    idx1poutre = elem_supp%aceind(ACE_NU_POUTRE, 1)
    idx2poutre = elem_supp%aceind(ACE_NU_POUTRE, 2)
    do ii = idx1poutre, idx2poutre
!       Si on trouve une poutre c'est OK sinon PLOUF
        if (nel .eq. elem_supp%catanum(ii)) goto 100
    end do
    call utmess('F', 'MODELISA_86')
100 continue
!   On décale de -1 pour éviter les +1 après
    idx1poutre = idx1poutre-1
!
    call jeveuo(jexnom(tmp, nom), 'E', jdge)
    isec = nint(zr(jdge+35))
!
!   calcul des caractéristiques générales section rectangulaire
!       -  erreur   si  hy  <= 0  ou  hz  <= 0  ou
!                       epy <= 0  ou  epz <= 0  (test dans affdef)
    if (isec .eq. 1) then
        do ii = 1, 2
            if (ii .eq. 1) then
                igeor = 24
                igeoc = 32
                igen = 1
                igen2 = 37
            else
                igeor = 28
                igeoc = 34
                igen = 12
                igen2 = 38
            end if
            zr(jdge+igeoc-1) = 0.d0
            zr(jdge+igeoc) = 0.d0
            hy = zr(jdge+igeor-1)
            hz = zr(jdge+igeor)
            epy = zr(jdge+igeor+1)
            epz = zr(jdge+igeor+2)
            hyi = hy-2.d0*epy
            hzi = hz-2.d0*epz
!           A
            zr(jdge+igen-1) = hy*hz-hyi*hzi
!           IY
            zr(jdge+igen) = (hy*(hz*hz*hz)-hyi*(hzi*hzi*hzi))/12.d0
!           IZ
            zr(jdge+igen+1) = (hz*(hy*hy*hy)-hzi*(hyi*hyi*hyi))/12.d0
!           EY
            zr(jdge+igen+4) = 0.d0
!           EZ
            zr(jdge+igen+5) = 0.d0
!           RY
            zr(jdge+igen+7) = hy/2.d0
!           RZ
            zr(jdge+igen+8) = hz/2.d0
!
!           cas de la section rectangulaire pleine
            if (abs(hyi/hy) .le. eps .or. abs(hzi/hz) .le. eps) then
                a = hy/2.d0
                b = hz/2.d0
                if (a/b .lt. 1.d0) then
                    a = hz/2.d0
                    b = hy/2.d0
                end if
                a4 = a**4*12.d0
                b4 = b**4
                b3 = b**3
                jx = a*b3*(16.d0/3.d0-3.36d0*b*(1.d0-b4/a4)/a)
!               AY & AZ
                if (nel .eq. elem_supp%catanum(idx1poutre+ACE_MECA_POU_D_T)) then
                    zr(jdge+igen+2) = 1.2d0
                    zr(jdge+igen+3) = 1.2d0
                else if (nel .eq. elem_supp%catanum(idx1poutre+ACE_MECA_POU_D_E)) then
                    zr(jdge+igen+2) = 0.d0
                    zr(jdge+igen+3) = 0.d0
                else if (nel .eq. elem_supp%catanum(idx1poutre+ACE_MECA_POU_D_T_GD)) then
                    zr(jdge+igen+2) = 1.2d0
                    zr(jdge+igen+3) = 1.2d0
                else if (nel .eq. elem_supp%catanum(idx1poutre+ACE_MEFS_POU_D_T)) then
                    zr(jdge+igen+2) = 1.2d0
                    zr(jdge+igen+3) = 1.2d0
                else if (nel .eq. elem_supp%catanum(idx1poutre+ACE_MECA_POU_D_TG)) then
                    zr(jdge+igen+2) = 1.2d0
                    zr(jdge+igen+3) = 1.2d0
                else if (nel .eq. elem_supp%catanum(idx1poutre+ACE_MECA_POU_D_EM)) then
                    zr(jdge+igen+2) = 0.d0
                    zr(jdge+igen+3) = 0.d0
                else if (nel .eq. elem_supp%catanum(idx1poutre+ACE_MECA_POU_D_TGM)) then
                    zr(jdge+igen+2) = 1.2d0
                    zr(jdge+igen+3) = 1.2d0
                else if (nel .eq. elem_supp%catanum(idx1poutre+ACE_MECA_S8P2)) then
                    zr(jdge+igen+2) = 1.2d0
                    zr(jdge+igen+3) = 1.2d0
                end if
!               JX
                zr(jdge+igen+6) = jx
!               RT
                zr(jdge+igen+9) = jx*(3.d0*a+1.8d0*b)/(8.d0*a*a*b*b)
!               AI
                zr(jdge+igen2-1) = 0.d0
            else
!               cas du tube rectangulaire ---
!               JX
                ct = 2.d0*epy*epz*(hy-epy)*(hy-epy)*(hz-epz)*(hz-epz)
                cd = hy*epy+hz*epz-epy*epy-epz*epz
                jx = ct/cd
                zr(jdge+igen+6) = jx
!               AY  : jdge+igen+2
!               AZ  : jdge+igen+3
                if (nel .eq. elem_supp%catanum(idx1poutre+ACE_MECA_POU_D_E) .or. &
                    nel .eq. elem_supp%catanum(idx1poutre+ACE_MECA_POU_D_EM)) then
                    zr(jdge+igen+2) = 0.d0
                    zr(jdge+igen+3) = 0.d0
                else
!                   interpolation des coefficients de cisaillement
                    alpha = (hy-2.d0*epy)/hy
                    beta = (hz-2.d0*epz)/hz
                    ASSERT((alpha .ge. 0.d0) .or. (beta .ge. 0.d0))
                    if (alpha .gt. 0.95d0 .or. beta .gt. 0.95d0) then
                        call utmess('F', 'MODELISA10_15')
                    end if
                    nompa(1) = 'ALPHA'
                    nompa(2) = 'BETA'
                    valpay(1) = alpha
                    valpay(2) = beta
                    valpaz(1) = beta
                    valpaz(2) = alpha
                    call fointe('A', napcis, 2, nompa, valpay, ay, ier)
                    call fointe('A', napcis, 2, nompa, valpaz, az, ier)
                    zr(jdge+igen+2) = ay
                    zr(jdge+igen+3) = az
                end if
!               RT. TUBE RECTANGULAIRE MINCE D’ÉPAISSEUR CONSTANTE. RT=JX/2.E.AINT
                aireint = (hy-2.d0*epy)*(hz-2.d0*epz)
                zr(jdge+igen+9) = jx/(2.d0*epz*aireint)
!               AI
                zr(jdge+igen2-1) = hyi*hzi
            end if
!           AY/AZ POUR TUYAUX ET 3D_FAISCEAU
            if (nel .eq. elem_supp%catanum(idx1poutre+ACE_MET3SEG3) .or. &
                nel .eq. elem_supp%catanum(idx1poutre+ACE_MET6SEG3) .or. &
                nel .eq. elem_supp%catanum(idx1poutre+ACE_MET3SEG4) .or. &
                nel .eq. elem_supp%catanum(idx1poutre+ACE_MECA_POHO_HEXA8) .or. &
                nel .eq. elem_supp%catanum(idx1poutre+ACE_MECA_POHO_HEXA20)) then
                zr(jdge+igen+2) = 0.d0
                zr(jdge+igen+3) = 0.d0
            end if
        end do
!       JG1,JG2,IYR21,IYR22,IZR21,IZR22 :
        do ii = 1, 6
            zr(jdge-1+38+ii) = 0.d0
        end do
    end if
!
!   calcul des caractéristiques générales section circulaire
!       erreur   si re <= 0  ou  e > re ou e <= 0  (test dans affdef)
    if (isec .eq. 2) then
        do ii = 1, 2
            if (ii .eq. 1) then
                igeor = 24
                igeoc = 32
                igen = 1
                igen2 = 37
            else
                igeor = 28
                igeoc = 34
                igen = 12
                igen2 = 38
            end if
            zr(jdge+igeor-1) = 0.d0
            zr(jdge+igeor) = 0.d0
            zr(jdge+igeor+1) = 0.d0
            zr(jdge+igeor+2) = 0.d0
            re = zr(jdge+igeoc-1)
            e = zr(jdge+igeoc)
            ri = re-e
!           A
            zr(jdge+igen-1) = pi*(re*re-ri*ri)
!           IY
            zr(jdge+igen) = pi*(re**4-ri**4)/4.d0
!           IZ
            zr(jdge+igen+1) = zr(jdge+igen)
!           EY
            zr(jdge+igen+4) = 0.d0
!           EZ
            zr(jdge+igen+5) = 0.d0
!           JX
            zr(jdge+igen+6) = zr(jdge+igen)*2.d0
!           RY
            zr(jdge+igen+7) = re
!           RZ
            zr(jdge+igen+8) = re
!           RT
            zr(jdge+igen+9) = re
!           AY et AZ
            if (nel .eq. elem_supp%catanum(idx1poutre+ACE_MECA_POU_D_E) .or. &
                nel .eq. elem_supp%catanum(idx1poutre+ACE_MECA_POU_D_EM)) then
                zr(jdge+igen+2) = 0.0d0
                zr(jdge+igen+3) = 0.0d0
            else
!               interpolation des coefficients de cisaillement
                alpha = ri/re
                ASSERT((alpha .ge. 0.d0) .or. (alpha .le. 1.d0))
                nompaf = 'ALPHA'
                valpaf = alpha
                call fointe('A', foncis, 1, [nompaf], [valpaf], ccis, ier)
                zr(jdge+igen+2) = ccis
                zr(jdge+igen+3) = ccis
            end if
!           AI
            zr(jdge+igen2-1) = pi*ri*ri
!           AY/AZ POUR TUYAUX ET 3D_FAISCEAU
            if (nel .eq. elem_supp%catanum(idx1poutre+ACE_MET3SEG3) .or. &
                nel .eq. elem_supp%catanum(idx1poutre+ACE_MET6SEG3) .or. &
                nel .eq. elem_supp%catanum(idx1poutre+ACE_MET3SEG4) .or. &
                nel .eq. elem_supp%catanum(idx1poutre+ACE_MECA_POHO_HEXA8) .or. &
                nel .eq. elem_supp%catanum(idx1poutre+ACE_MECA_POHO_HEXA20)) then
                zr(jdge+igen+2) = 0.d0
                zr(jdge+igen+3) = 0.d0
            end if
!
        end do
!       JG1,JG2,IYR21,IYR22,IZR21,IZR22 :
        do ii = 1, 6
            zr(jdge-1+38+ii) = 0.d0
        end do
    end if
!
    call jedema()
end subroutine
