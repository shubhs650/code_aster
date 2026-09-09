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

subroutine affdef(tmp, nom, nel, tab, ier)
!
    use cara_elem_parameter_module
!
    implicit none
#include "jeveux.h"
#include "asterc/r8maem.h"
#include "asterfort/jedema.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnom.h"
#include "asterfort/utmess.h"
!
    character(len=8) :: tab(*)
    character(len=24) :: nom, tmp
!
!  - VERIFICATION DE LA COMPLETUDE DES DONNES OBLIGATOIRES A
!    ENTRER , ET DE LA POSITIVITE DE CERTAINES VALEURS
!  - AFFECTATION DES VALEURS PAR DEFAUT AUX CARACTERISTIQUES
!    GENERALES ET GEOMETRIQUES ABSENTES
!  - NTEL(I)  = NUMERO DU TYPE ELEMENT
!       I = 1 : MECA_POU_D_T
!           2 : MECA_POU_D_E
!           4 : MEFS_POU_D_T
!           5 : MECA_POU_D_TG
!          11 : MECA_POU_D_EM
!          12 : MECA_POU_D_TGM
!          14 : MECA_S8P2
!
! --- ------------------------------------------------------------------
!     TAB  1      2      3      4      5     6     7     8    9    10
!     0    A1     IY1    IZ1    AY1    AZ1   EY1   EZ1   JX1  RY1  RZ1
!     1    RT1    A2     IY2    IZ2    AY2   AZ2   EY2   EZ2  JX2  RY2
!     2    RZ2    RT2    TVAR   HY1    HZ1   EPY1  EPZ1  HY2  HZ2  EPY2
!     3    EPZ2   R1     E1     R2     E2    TSEC  AI1   AI2  JG1  JG2
!     4    IYR21  IYR22  IZR21  IZR22
!--- -------------------------------------------------------------------
    integer(kind=8) :: ier, isec, j, jdge, nc, nd, ne, idx1poutre
    integer(kind=8) :: nel, ng, nr, nt, nx, ny, nz
!--- -------------------------------------------------------------------
    parameter(nr=4, nc=2, ng=8)
    parameter(nt=4, ne=12, nd=6)
    parameter(nx=10, ny=8, nz=4)
    integer(kind=8) :: ogen(ng), orec(nr), ocer(nc), otpe(nt)
    integer(kind=8) :: dexc(ne), ddfx(nd), drec(nr), dcer(nc)
    integer(kind=8) :: pgen(nx), prec(ny), pcer(nz)
    character(len=24) :: valk(3)
    real(kind=8) :: tst
!--- -------------------------------------------------------------------
    data ogen/1, 2, 3, 8, 12, 13, 14, 19/
    data orec/24, 25, 28, 29/
    data ocer/32, 34/
    data otpe/4, 5, 15, 16/
    data dexc/6, 7, 17, 18, 37, 38, 39, 40, 41, 42, 43, 44/
    data ddfx/9, 10, 11, 20, 21, 22/
    data drec/26, 27, 30, 31/
    data dcer/33, 35/
    data pgen/1, 9, 10, 11, 12, 20, 21, 22, 37, 38/
    data prec/24, 25, 26, 27, 28, 29, 30, 31/
    data pcer/32, 33, 34, 35/
!--- -------------------------------------------------------------------
!
    call jemarq()
!   On enleve 1% au Real maxi : il reste maxi
!       Modifs des tests    if (real .eq. tst)  ==> if (real .ge. tst)
!                           if (real .ne. tst)  ==> if (real .le. tst)
    tst = r8maem()*0.99
!
    call jeveuo(jexnom(tmp, nom), 'E', jdge)
    isec = nint(zr(jdge+35))
!
!   Index de départ des poutres dans la SD cara_elem
!       On décale de -1, pour éviter les +1 après
    idx1poutre = elem_supp%aceind(ACE_NU_POUTRE, 1)-1
!
!    1:MECA_POU_D_T
!    4:MEFS_POU_D_T
!    5:MECA_POU_D_TG
!   12:MECA_POU_D_TGM
!    2:MECA_POU_D_E
!   11:MECA_POU_D_EM

!   COMPLETUDE DES DONNES GENERALES
    if (isec .eq. 0) then
        do j = 1, ng
            if (zr(jdge+ogen(j)-1) .ge. tst) then
                valk(1) = nom
                valk(2) = tab(ogen(j))
                call utmess('A', 'MODELISA_77', nk=2, valk=valk)
                ier = ier+1
            end if
        end do
        if ((nel .eq. elem_supp%catanum(idx1poutre+ACE_MECA_POU_D_T)) .or. &
            (nel .eq. elem_supp%catanum(idx1poutre+ACE_MEFS_POU_D_T)) .or. &
            (nel .eq. elem_supp%catanum(idx1poutre+ACE_MECA_POU_D_TG)) .or. &
            (nel .eq. elem_supp%catanum(idx1poutre+ACE_MECA_POU_D_TGM)) .or. &
            (nel .eq. elem_supp%catanum(idx1poutre+ACE_MECA_S8P2))) then
            do j = 1, nt
                if (zr(jdge+otpe(j)-1) .ge. tst) then
                    valk(1) = nom
                    valk(2) = tab(otpe(j))
                    call utmess('A', 'MODELISA_78', nk=2, valk=valk)
                end if
            end do
        end if
    end if
!
!   COMPLETUDE DES DONNES GEOMETRIQUES RECTANGLE
    if (isec .eq. 1) then
        do j = 1, nr
            if (zr(jdge+orec(j)-1) .ge. tst) then
                valk(1) = nom
                valk(2) = tab(orec(j))
                call utmess('A', 'MODELISA_79', nk=2, valk=valk)
                ier = ier+1
            end if
        end do
    end if
!
!   COMPLETUDE DES DONNES GEOMETRIQUES CERCLE
    if (isec .eq. 2) then
        do j = 1, nc
            if (zr(jdge+ocer(j)-1) .ge. tst) then
                valk(1) = nom
                valk(2) = tab(ocer(j))
                call utmess('A', 'MODELISA_80', nk=2, valk=valk)
                ier = ier+1
            end if
        end do
    end if
!
!   VERIFICATION DE LA STRICTE POSITIVITE DE  VALEURS GENERALE
    if (isec .eq. 0) then
        do j = 1, nx
            if (zr(jdge+pgen(j)-1) .le. tst) then
                if (zr(jdge+pgen(j)-1) .le. 0.d0) then
                    valk(1) = nom
                    valk(2) = tab(pgen(j))
                    call utmess('A', 'MODELISA_81', nk=2, valk=valk)
                    ier = ier+1
                end if
            end if
        end do
    end if
!
!     VERIFICATION DE LA STRICTE POSITIVITE DE VALEURS RECTANGLE
    if (isec .eq. 1) then
        do j = 1, ny
            if (zr(jdge+prec(j)-1) .le. tst) then
                if (zr(jdge+prec(j)-1) .le. 0.d0) then
                    valk(1) = nom
                    valk(2) = tab(prec(j))
                    call utmess('A', 'MODELISA_82', nk=2, valk=valk)
                    ier = ier+1
                end if
            end if
        end do
    end if
!
!     VERIFICATION DE LA STRICTE POSITIVITE DE VALEURS CERCLE
    if (isec .eq. 2) then
        do j = 1, nz
            if (zr(jdge+pcer(j)-1) .le. tst) then
                if (zr(jdge+pcer(j)-1) .le. 0.d0) then
                    valk(1) = nom
                    valk(2) = tab(pcer(j))
                    call utmess('A', 'MODELISA_83', nk=2, valk=valk)
                    ier = ier+1
                end if
            end if
        end do
    end if
!
    if (ier .ne. 0) goto 999
!
!     AFFECTATION DES VALEURS PAR DEFAUT POUR LES DONNEES GENERALES
    if (isec .eq. 0) then
!        EXCENTREMENTS, AIRES INTERIEURES, CONSTANTES DE GAUCHISSEMENT
        do j = 1, ne
            if (zr(jdge+dexc(j)-1) .ge. tst) zr(jdge+dexc(j)-1) = 0.d0
        end do
!        DIST. FIBRE EXT.+ RAYON TORSION
        do j = 1, nd
            if (zr(jdge+ddfx(j)-1) .ge. tst) zr(jdge+ddfx(j)-1) = 1.d0
        end do
!        EULER
        if ((nel .eq. elem_supp%catanum(idx1poutre+ACE_MECA_POU_D_E)) .or. &
            (nel .eq. elem_supp%catanum(idx1poutre+ACE_MECA_POU_D_EM))) then
            do j = 1, nt
                zr(jdge+otpe(j)-1) = 0.d0
            end do
        end if
    end if
!
!     AFFECTATION DES VALEURS PAR DEFAUT POUR LES DONNEES RECTANGLE
    if (isec .eq. 1) then
        do j = 1, nr
            if (zr(jdge+drec(j)-1) .ge. tst) then
                zr(jdge+drec(j)-1) = zr(jdge+orec(j)-1)/2.d0
            else
                if (zr(jdge+drec(j)-1) .gt. (zr(jdge+orec(j)-1)/2.d0)) then
                    valk(1) = nom
                    valk(2) = tab(drec(j))
                    valk(3) = tab(orec(j))
                    call utmess('A', 'MODELISA_84', nk=3, valk=valk)
                    ier = ier+1
                end if
            end if
        end do
    end if
!
!     AFFECTATION DES VALEURS PAR DEFAUT POUR LES DONNEES CERCLE
    if (isec .eq. 2) then
        do j = 1, nc
            if (zr(jdge+dcer(j)-1) .ge. tst) then
                zr(jdge+dcer(j)-1) = zr(jdge+ocer(j)-1)
            else
                if (zr(jdge+dcer(j)-1) .gt. zr(jdge+ocer(j)-1)) then
                    valk(1) = nom
                    valk(2) = tab(dcer(j))
                    valk(3) = tab(ocer(j))
                    call utmess('A', 'MODELISA_85', nk=3, valk=valk)
                    ier = ier+1
                end if
            end if
        end do
    end if
!
999 continue
    call jedema()
end subroutine
