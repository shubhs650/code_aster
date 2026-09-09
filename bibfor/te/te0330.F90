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
subroutine te0330(option, nomte)
    implicit none
!
!     CALCUL DES TERMES ELEMENTAIRES DE LAGRANGE
!     ELEMENTS  DDEPL_R_.... / DTEMP_R_.... / DPRES_C_....
!
!     EN ENTREE :
!        OPTION : NOM DE L'OPTION A CALCULER
!        NOMTE  : NOM DU TYPE_ELEMENT
!
!
!
!
!
!
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/fointe.h"
#include "asterfort/jevech.h"
!
    real(kind=8) :: k(10, 10), valpar(4), result
    integer(kind=8) :: i, ier, j, ndl, ip, jmat, jdimp, jdmul, jgeom, jlagr
    integer(kind=8) :: jvec, jtime, nbpar
    complex(kind=8) :: kc(10, 10)
!
    character(len=8) :: nompar(4)
    character(len=19) :: nomfon
    character(len=16), intent(in) :: option
    character(len=16), intent(in) :: nomte
! ----------------------------------------------------------------------
!
    ndl = 1
!
    if (option(6:12) .eq. 'DDLM_R') then
        call jevech('PDDLMUR', 'L', jdmul)
        k = 0.d0
!
        do i = 1, ndl
            k(i, ndl+1) = zr(jdmul-1+i)
        end do
        if (option(1:4) .eq. 'MECA') call jevech('PMATUUR', 'E', jmat)
        if (option(1:4) .eq. 'THER') call jevech('PMATTTR', 'E', jmat)
        ip = 0
        do i = 1, ndl+1
            do j = 1, i
                zr(jmat-1+ip+j) = k(j, i)
            end do
            ip = ip+i
        end do
!
    else if (option(6:12) .eq. 'DDLM_C') then
        ASSERT(.false.)
        call jevech('PDDLMUC', 'L', jdmul)
        kc = (0.d0, 0.d0)
!
        do i = 1, ndl
            kc(i, ndl+1) = zc(jdmul-1+i)
        end do
        if (option(1:4) .eq. 'ACOU') call jevech('PMATTTC', 'E', jmat)
        ip = 0
        do i = 1, ndl+1
            do j = 1, i
                zc(jmat-1+ip+j) = kc(j, i)
            end do
            ip = ip+i
        end do
!
    else if (option(6:11) .eq. 'BTLA_R') then
        call jevech('PDDLMUR', 'L', jdmul)
        call jevech('PLAGRAR', 'L', jlagr)
        if (option(1:4) .eq. 'MECA') call jevech('PVECTUR', 'E', jvec)
        if (option(1:4) .eq. 'THER') call jevech('PVECTTR', 'E', jvec)
        zr(jvec) = zr(jdmul)*zr(jlagr+1)
        zr(jvec+1) = 0.d0
!
    else if (option(6:9) .eq. 'BU_R') then
        call jevech('PDDLMUR', 'L', jdmul)
        call jevech('PDDLIMR', 'L', jdimp)
        if (option(1:4) .eq. 'MECA') call jevech('PVECTUR', 'E', jvec)
        if (option(1:4) .eq. 'THER') call jevech('PVECTTR', 'E', jvec)
        zr(jvec) = 0.d0
        zr(jvec+1) = zr(jdmul)*zr(jdimp)
!
    else if (option(6:11) .eq. 'DDLI_R') then
        call jevech('PDDLIMR', 'L', jdimp)
        if (option(1:4) .eq. 'MECA') call jevech('PVECTUR', 'E', jvec)
        if (option(1:4) .eq. 'THER') call jevech('PVECTTR', 'E', jvec)
        zr(jvec-1+ndl+1) = zr(jdimp-1+1)
!
    else if (option(6:11) .eq. 'DDLI_C') then
        call jevech('PDDLIMC', 'L', jdimp)
        if (option(1:4) .eq. 'MECA') call jevech('PVECTUC', 'E', jvec)
        if (option(1:4) .eq. 'THER') call jevech('PVECTTC', 'E', jvec)
        if (option(1:4) .eq. 'ACOU') call jevech('PVECTTC', 'E', jvec)
        zc(jvec-1+ndl+1) = zc(jdimp-1+1)
!
    else if (option(6:11) .eq. 'DDLI_F') then
        call jevech('PDDLIMF', 'L', jdimp)
        call jevech('PGEOMER', 'L', jgeom)
        call jevech('PINSTR', 'L', jtime)
        if (option(1:4) .eq. 'MECA') call jevech('PVECTUR', 'E', jvec)
        if (option(1:4) .eq. 'THER') call jevech('PVECTTR', 'E', jvec)
        if (option(1:4) .eq. 'ACOU') call jevech('PVECTTC', 'E', jvec)
        nomfon = zk24(jdimp-1+1) (1:19)
        nbpar = 4
        nompar(1) = 'X'
        nompar(2) = 'Y'
        nompar(3) = 'Z'
        nompar(4) = 'INST'
        valpar(1) = zr(jgeom-1+1)
        valpar(2) = zr(jgeom-1+2)
        valpar(3) = zr(jgeom-1+3)
        valpar(4) = zr(jtime-1+1)
        call fointe('FM', nomfon, nbpar, nompar, valpar, &
                    result, ier)
        if (option(1:4) .eq. 'ACOU') then
            zc(jvec-1+ndl+1) = result
        else
            zr(jvec-1+ndl+1) = result
        end if
    end if
! ----------------------------------------------------------------------
end subroutine
