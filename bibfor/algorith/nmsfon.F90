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

subroutine nmsfon(refe, ndim, nno, npg, nddl, &
                  geomi, vff, idff, iw, sief, fint)
!
    use bloc_fe_module, only: prod_bd, prod_sb, prod_bkb, add_fint, add_matr
!
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/dfdmip.h"
#include "asterfort/nmbeps.h"
#include "asterfort/teattr.h"
#include "jeveux.h"

    integer(kind=8), intent(in)             :: ndim, nno, npg, nddl
    integer(kind=8), intent(in)             :: iw, idff
    real(kind=8), intent(in)        :: geomi(ndim, nno), vff(nno, npg)
    real(kind=8), intent(in)        :: sief(4*ndim, npg)
    aster_logical, intent(in)       :: refe
    real(kind=8), intent(out)       :: fint(nddl)
!
! --------------------------------------------------------------------------------------------------
!
!     RAPH_MECA, RIGI_MECA_* ET FULL_MECA_* , ELEMENTS MIX_STA
!
! --------------------------------------------------------------------------------------------------
!
! IN  FAMI    : FAMILLE DE POINTS DE GAUSS
! IN  OPTION  : OPTION DE CALCUL
! IN  TYPMOD  : TYPE DE MODELISATION
! IN  NDIM    : DIMENSION DE L'ESPACE
! IN  NNO     : NOMBRE DE NOEUDS STANDARDS D'UN ELEMENT
! IN  NPG     : NOMBRE DE POINTS DE GAUSS
! IN  NDDL    : DEGRES DE LIBERTE D'UN ELEMENT ENRICHI
! IN  IW      : PTR. POIDS DES POINTS DE GAUSS
! IN  VFF     : VALEUR  DES FONCTIONS DE FORME DE DEPLACEMENT
! IN  IDFF    : PTR. DERIVEE DES FONCTIONS DE FORME DE DEPLACEMENT ELEMENT DE REF.
! IN  GEOMI   : COORDONNEES DES NOEUDS (CONFIGURATION INITIALE)
! IN  COMPOR  : COMPORTEMENT
! IN  MATE    : MATERIAU CODE
! IN  LGPG    : DIMENSION DU VECTEUR DES VAR. INTERNES POUR 1 PT GAUSS
! IN  CRIT    : CRITERES DE CONVERGENCE LOCAUX
! IN  ANGMAS  : LES TROIS ANGLES DU MOT_CLEF MASSIF (AFFE_CARA_ELEM)
! IN  INSTM   : VALEUR DE L'INSTANT T-
! IN  INSTP   : VALEUR DE L'INSTANT T+
! IN  MATSYM  : .TRUE. SI MATRICE SYMETRIQUE
! IN  DDLM    : DDL AU PAS T-
! IN  DDLD    : INCREMENT DE DDL ENTRE T- ET T+
! IN  SIGMG    : CONTRAINTES GENERALISEES EN T-
!                SIGMG(1:2*NDIM) CAUCHY
!                SIGMG(2*NDIM,NPES) : SIG_A, SIG_LAM
! IN  VIM     : VARIABLES INTERNES EN T-
! OUT SIGPG    : CONTRAINTES GENERALIEES (RAPH_MECA ET FULL_MECA_*)
! OUT VIP     : VARIABLES INTERNES    (RAPH_MECA ET FULL_MECA_*)
! OUT FINT    : FORCES INTERIEURES (RAPH_MECA ET FULL_MECA_*)
! OUT MATR   : MATR. DE RIGIDITE NON SYM. (RIGI_MECA_* ET FULL_MECA_*)
! OUT CODRET  : CODE RETOUR DE L'INTEGRATION DE LA LDC
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), parameter:: rac2 = sqrt(2.d0)
    real(kind=8), dimension(6), parameter  :: vrac2 = (/1.d0, 1.d0, 1.d0, &
                                                        sqrt(2.d0), sqrt(2.d0), sqrt(2.d0)/)
! ----------------------------------------------------------------------
    aster_logical :: axi
    character(len=16) :: formulation
    integer(kind=8)       :: g, n, i, j, cod(npg)
    integer(kind=8)       :: xu(ndim, nno), xe(2*ndim, nno)
    integer(kind=8)       :: ndu, nde
    real(kind=8)  :: r, dff(nno, ndim), poids, coeff_trace
    real(kind=8)  :: bu(2*ndim, ndim, nno), be(2*ndim, 2*ndim, nno)
    real(kind=8)  :: b_dev(2*ndim, ndim, nno), b_vol(2*ndim, ndim, nno)
    real(kind=8)  :: b(2*ndim, ndim, nno), bbar_vol(2*ndim, ndim, nno)
    real(kind=8)  :: siefu(2*ndim), siefe(2*ndim)
! --------------------------------------------------------------------------------------------------
!

! --- INITIALISATION ---
    axi = ASTER_FALSE
    call teattr('S', 'FORMULATION', formulation)
    cod = 0.d0
    coeff_trace = 1.d0/ndim

    ! Nombre de ddls
    ndu = ndim
    nde = 2*ndim

    fint = 0

    ! tableaux de reference bloc (depl,epsi) -> numero du ddl
    forall (i=1:ndu, n=1:nno) xu(i, n) = (n-1)*(ndu+nde)+i
    forall (i=1:nde, n=1:nno) xe(i, n) = (n-1)*(ndu+nde)+ndu+i

    ! Formulation STA_INCO : calcul de Bbar_vol
    if (formulation .eq. "STA_INCO") then
        bbar_vol = 0.d0
        do g = 1, npg
            b_vol = 0.0d0
            b = 0.d0
            ! Calcul des derivees des fonctions de forme P1, du rayon r et des poids
            call dfdmip(ndim, nno, axi, geomi, g, iw, vff(1, g), idff, r, poids, dff)
            ! Calcul de la partie volumique de B
            forall (i=1:ndim, j=1:ndim) b_vol(i, j, :) = coeff_trace*dff(:, j)
            ! Calcul de Bbar_vol
            bbar_vol = bbar_vol+1.d0/npg*b_vol
        end do
    end if

    gauss: do g = 1, npg
        bu = 0.d0
        ! Calcul des derivees des fonctions de forme P1, du rayon r et des poids
        call dfdmip(ndim, nno, axi, geomi, g, iw, vff(1, g), idff, r, poids, dff)

        ! Formulation STA : calcul de la matrice Bu = B
        if (formulation .eq. "STA") then
            call nmbeps(axi, r, vff(:, g), dff, bu)

            ! Formulation STA_INCO : calcul de la matrice Bu = B_dev + Bbar_vol
        else if (formulation .eq. "STA_INCO") then
            b = 0.d0
            b_vol = 0.d0
            b_dev = 0.d0
            ! Calcul de la matrice B
            call nmbeps(axi, r, vff(:, g), dff, b)
            ! Calcul de la partie volumique de B
            forall (i=1:ndim, j=1:ndim) b_vol(i, j, :) = coeff_trace*dff(:, j)
            ! Calcul de la partie déviatorique de B
            b_dev = b-b_vol
            ! Calcul de Bu = B_dev + Bbar_vol
            bu = b_dev+bbar_vol
        else
            ASSERT(ASTER_FALSE)
        end if

        ! Calcul de la matrice N
        be = 0.d0
        if (ndim .eq. 2) then
            be(1, 1, :) = vff(:, g)
            be(2, 2, :) = vff(:, g)
            be(3, 3, :) = vff(:, g)
            be(4, 4, :) = rac2*vff(:, g)
        else if (ndim .eq. 3) then
            be(1, 1, :) = vff(:, g)
            be(2, 2, :) = vff(:, g)
            be(3, 3, :) = vff(:, g)
            be(4, 4, :) = rac2*vff(:, g)
            be(5, 5, :) = rac2*vff(:, g)
            be(6, 6, :) = rac2*vff(:, g)
        end if

        ! Extraction des blocs de contraintes generalisees
        siefu = sief(1:2*ndim, g)*vrac2(1:2*ndim)
        siefe = sief(2*ndim+1:4*ndim, g)

        ! Correction matrices B si REFE_FORC_NODA
        if (refe) then
            bu = abs(bu)
            be = abs(be)
        end if

        ! Calcul des contributions aux forces interieures
        call add_fint(fint, xu, poids*prod_sb(siefu, bu))
        call add_fint(fint, xe, poids*prod_sb(siefe, be))

    end do gauss

end subroutine
