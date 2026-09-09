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
subroutine mbxnlr(plateOrie, &
                  option, fami, &
                  nddl, nno, ncomp, kpg, &
                  ipoids, jvGeom, &
                  jvMaterc, jvDispM, jvDispIncr, jvVect, jvSigm, &
                  jvMatrSyme, dff, &
                  lVect, lMatr)
!
    use plate_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/mbcine.h"
#include "asterfort/mbrigi.h"
#include "asterfort/verift.h"
#include "jeveux.h"
!
    type(plateOrie_Para), intent(in) :: plateOrie
    character(len=16), intent(in) :: option
    character(len=8), intent(in) :: fami
    integer(kind=8), intent(in) :: nddl, nno, ncomp, kpg
    integer(kind=8), intent(in) :: ipoids, jvGeom, jvMaterc, jvDispM, jvDispIncr
    integer(kind=8), intent(in) :: jvVect, jvSigm, jvMatrSyme
    real(kind=8), intent(in) :: dff(2, nno)
    aster_logical, intent(in) :: lVect, lMatr
!
! --------------------------------------------------------------------------------------------------
!
!    - FONCTION REALISEE:  CALCUL DES OPTIONS DE COMPORTEMENT
!                          POUR LES MEMBRANES EN PETITES DEFORMATIONS
!
! --------------------------------------------------------------------------------------------------
!
! IN  OPTION       OPTION DE CALCUL
! IN  FAMI         NOM DE LA FAMILLE DE POINTS DE GAUSS :
!                  'RIGI','MASS',..
! IN  NDDL         NOMBRE DE DERGES DE LIBERTE AUX NOEUDS
! IN  NNO          NOMBRE DE NOEUDS
! IN  NCOMP        NOMBRE DE COMPOSANTS DANS LES VECTEURS COLONNES
!                  DE CONTRAINTE ET DEFORMATION
! IN  KPG          INCREMENT SUR LA BOUCLE DES PTS DE GAUSS
! IN  IPOIDS       ADRESSE DANS ZR DU TABLEAU POIDS
! IN  IGEOM        ADRESSE DANS ZR DU TABLEAU PGEOMER
! IN  IMATE        ADRESSE DANS ZI DU TABLEAU PMATERC
! IN  IDEPLM       ADRESSE DANS ZR DU TABLEAU PDEPLMR
! IN  IDEPLP       ADRESSE DANS ZR DU TABLEAU PDEPLPR
! IN  IVECTU       ADRESSE DANS ZR DU TABLEAU PVECTUR
! IN  ICONTP       ADRESSE DANS ZR DU TABLEAU PCONTPR
! IN  IMATUU       ADRESSE DANS ZR DU TABLEAU PMATUUR
! IN  DFF          DERIVEE DES F. DE FORME
! IN  ALPHA, BETA  ANGLES NAUTIQUES ORIENTANT LE COMPORTEMENT
!                        ORTHOTROPE DE LA MEMBRANE (EN RADIAN)
! IN  VECTEU       BOOL: 1 SI FULL_MECA OU RAPH_MECA
! IN  MATRIC       BOOL: 1 SI FULL_MECA OU RIGI_MECA
!
! OUT ***          ***
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i, j, j1, n, m, c, cc, kk, kkd
    real(kind=8) :: b(3, 3, 9), jac
    real(kind=8) :: epsm(3), deps(3), epsth(3), epsthe, sigp(3), tmp, matrRigi(3, 3)
!
! --------------------------------------------------------------------------------------------------
!

! - CALCUL DE LA MATRICE "B"
    call mbcine(plateOrie, &
                nno, zr(jvGeom), dff, &
                b, jac)

    if (option .eq. 'RIGI_MECA') then
        call mbrigi(fami, kpg, jvMaterc, matrRigi)

    elseif ((option(1:9) .eq. 'RAPH_MECA') .or. (option(1:9) .eq. 'FULL_MECA') .or. &
            (option(1:10) .eq. 'RIGI_MECA_')) then

        epsm = 0.d0
        deps = 0.d0
        do n = 1, nno
            do i = 1, nddl
                do c = 1, ncomp
                    epsm(c) = epsm(c)+b(c, i, n)*zr(jvDispM+(n-1)*nddl+i-1)
                    deps(c) = deps(c)+b(c, i, n)*zr(jvDispIncr+(n-1)*nddl+i-1)
                end do
            end do
        end do
!
        call verift(fami, kpg, 1, '+', zi(jvMaterc), &
                    epsth_=epsthe)
        epsth = 0.d0
        epsth(1) = epsthe
        epsth(2) = epsthe
!
        call mbrigi(fami, kpg, jvMaterc, matrRigi)
        sigp = 0.d0
        do c = 1, ncomp
            do cc = 1, ncomp
                sigp(c) = sigp(c)+ &
                          (epsm(cc)+deps(cc)-epsth(cc))*matrRigi(cc, c)
            end do
        end do
        if ((option(1:9) .eq. 'RAPH_MECA') .or. (option(1:9) .eq. 'FULL_MECA')) then
            do c = 1, ncomp
                zr(jvSigm+(kpg-1)*ncomp+c-1) = sigp(c)
            end do
        end if
    end if

    if (lVect) then
        do n = 1, nno
            do i = 1, nddl
                do c = 1, ncomp
                    zr(jvVect+(n-1)*nddl+i-1) = zr(jvVect+(n-1)*nddl+i-1)+ &
                                                b(c, i, n)*sigp(c)*zr(ipoids+kpg-1)*jac
                end do
            end do
        end do
    end if
!
    if (lMatr) then
        do n = 1, nno
            do i = 1, nddl
                kkd = (nddl*(n-1)+i-1)*(nddl*(n-1)+i)/2
                do j = 1, nddl
                    do m = 1, n
                        if (m .eq. n) then
                            j1 = i
                        else
                            j1 = nddl
                        end if
                        tmp = 0.d0
                        do c = 1, ncomp
                            do cc = 1, ncomp
                                tmp = tmp+ &
                                      b(cc, i, n)*matrRigi(cc, c)*b(c, j, m)*zr(ipoids+kpg-1)*jac
                            end do
                        end do
                        if (j .le. j1) then
                            kk = kkd+nddl*(m-1)+j
                            zr(jvMatrSyme+kk-1) = zr(jvMatrSyme+kk-1)+tmp
                        end if
                    end do
                end do
            end do
        end do
    end if

end subroutine
