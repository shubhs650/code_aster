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

subroutine nmspse(ndim, nno, nddl, &
                  nno_p, nno_s, nddl_s, npg, &
                  vff_s, vf_p, pgl, geom, neps, mate, matpou, &
                  deplm, urpg)
!
! aslint: disable=W1306
!
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterfort/nmspci.h"
#include "asterfort/ffpoutimo.h"
#include "asterfort/utpvgl.h"
#include "asterfort/utpvlg.h"
#include "asterfort/lonelesp.h"
!
    integer(kind=8) :: ndim, nno, nddl, nno_p, nno_s, nddl_s, npg
    integer(kind=8) :: neps, mate
    real(kind=8) :: geom(3*nno), vf_p(ndim, npg), vff_s(nno, npg)
    real(kind=8) :: pgl(3, 3), deplm(nddl)
    real(kind=8) :: urpg(neps, npg)
    character(len=8), intent(in)   :: matpou
!
!-----------------------------------------------------------------------
!  CALCUL DE SAUT_ELGA POUR LES ELEMENTS 3D_INTERF_POU (TE0383)
!-----------------------------------------------------------------------
! IN  NDIM    DIMENSION DU PROBLEME (=3)
! IN  NNO     NOMBRE DE NOEUDS TOTAL DE L'ELEMENT
! IN  NDDL    NOMBRE DE DEGRES DE LIBERTE EN DEPL TOTAL (3 PAR NOEUD)
! IN  NNO_P   NOMBRE DE NOEUDS DE LA PARTIE POUTRE
! IN  NNO_S   NOMBRE DE NOEUDS DE LA PARTIE SOL
! IN  NDDL_S  NOMBRE DE DEGRES DE LIBERTE EN DEPL DU SOL (3 PAR NOEUD)
! IN  NPG     NOMBRE DE POINTS DE GAUSS
! IN  VF_P    COORDONNES DES POINTS DE GAUSS
! IN  VFF_S   VALEUR DES FONCTIONS DE FORME DE L'ELEMENT SOL
! IN  PGL     MATRICE DE PROJECTION DU REPERE GLOBAL --> LOCAL
! IN  GEOM    COORDONNEES GEOMETRIQUES DES NOEUDS
! IN  NEPS    DIMENSION DES DEFORMATIONS/DEPLACEMENTS RELATIFS
! IN  MATE    MATERIAU CODE
! IN  MATPOU  NOM DU MATERIAU DE LA POUTRE
! IN  DEPLM   DEPLACEMENTS AUX NOEUDS
! OUT URPG    DEPLACEMENTS RELATIFS AUX POINTS DE GAUSS
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: kpg, n, iddl
    integer(kind=8) :: ddl_sp(3*nno_s+6*nno_p), ddl_s(3*nno_s)
    integer(kind=8) :: ddl_p(6*nno_p)
    integer(kind=8) :: no1, no2
    real(kind=8) :: xl
    real(kind=8) :: vff_p(18)
    real(kind=8) :: ur(3)
    real(kind=8) :: b(3, nddl)
    real(kind=8) :: deplmrot(nddl), urpgrot(neps, npg)

    ur = 0.d0
    ddl_s = [(iddl, iddl=1, 3*nno_s, 1)]
    ddl_p = (/7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6/)
    ddl_sp = [ddl_s, ddl_p+3*nno_s]

! - Get pile length
    no1 = 23
    no2 = 25
    xl = lonelesp(geom, 27, no1, no2)

! - Rotate displacements
    call utpvgl(nno_s, 3, pgl, deplm, deplmrot)
    call utpvgl(nno_p, 6, pgl, deplm(nddl_s+1), deplmrot(nddl_s+1))

!
! - Loop on Gauss points
!
    do kpg = 1, npg

! ----- Get matrix B giving relative displacement from node displacements
        call ffpoutimo(vf_p(1, kpg), xl, mate, matpou, vff_p)
        call nmspci(nno_p, nno_s, vff_p, vff_s(1, kpg), b)

! ----- Calculation of relative displacement
        ur = matmul(b, deplmrot(ddl_sp))

! ----- Out
        do n = 1, 3
            urpgrot(n, kpg) = ur(n)
        end do

    end do

! - Rotate relative displacement to global
    call utpvlg(npg, neps, pgl, urpgrot, urpg)

end subroutine
