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
! aslint: disable=W1306,W1504
!
subroutine fneihm(ds_thm, &
                  nno1, nno2, &
                  npi, npg, wref, iu, ip, &
                  ipf, iq, vff1, vff2, dffr2, &
                  geom, congem, r, vectu, &
                  mecani, press1, press2, dimdef, &
                  dimcon, dimuel, ndim, axi)
!
    use THM_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/fonoei.h"
#include "asterfort/matthm.h"
!
    type(THM_DS), intent(in) :: ds_thm
    aster_logical :: axi
    integer(kind=8) :: dimdef, dimcon, nno1, nno2
    integer(kind=8) :: dimuel, ndim
    integer(kind=8) :: npi, npg, mecani(8), press1(9), press2(9)
    integer(kind=8) :: addeme, addep1, addep2
    integer(kind=8) :: iu(3, 18), ip(2, 9), ipf(2, 2, 9), iq(2, 2, 9)
    real(kind=8) :: geom(ndim, nno2), dffr2(ndim-1, nno2, npi)
    real(kind=8) :: congem(dimcon, npi), vff1(nno1, npi), vff2(nno2, npi)
    real(kind=8) :: vectu(dimuel), r(dimdef), wref(npg)
!
! --------------------------------------------------------------------------------------------------
!
! BUT:  CALCUL  DE L'OPTION FORC_NODA POUR JOINT AVEC COUPLAGE HM
!
! --------------------------------------------------------------------------------------------------
!
! In  ds_thm           : datastructure for THM
! AXI       AXISYMETRIQUE ?
! TYPMOD    MODELISATION (D_PLAN, AXI, 3D ?)
! MODINT    METHODE D'INTEGRATION (CLASSIQUE,LUMPEE(D),REDUITE(R) ?)
! NNO       NB DE NOEUDS DE L'ELEMENT
! NNOS      NB DE NOEUDS SOMMETS DE L'ELEMENT
! NNOM      NB DE NOEUDS MILIEUX DE L'ELEMENT
! NDDLS     NB DE DDL SUR LES SOMMETS
! NDDLM     NB DE DDL SUR LES MILIEUX
! NPI       NB DE POINTS D'INTEGRATION DE L'ELEMENT
!
! NDIM      DIMENSION DE L'ESPACE
! DIMUEL    NB DE DDL TOTAL DE L'ELEMENT
! DIMCON    DIMENSION DES CONTRAINTES GENERALISEES ELEMENTAIRES
! DIMDEF    DIMENSION DES DEFORMATIONS GENERALISEES ELEMENTAIRES
! IVF       FONCTIONS DE FORMES QUADRATIQUES
! IVF2      FONCTIONS DE FORMES LINEAIRES
! OUT R       : TABLEAU DES RESIDUS
! OUT VECTU   : FORCES NODALES
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: adcome, addlh1, adcop1, adcop2
    integer(kind=8) :: kpi, i, n
    real(kind=8) :: wi, q(dimdef, dimuel)
!
! --------------------------------------------------------------------------------------------------
!

! - DETERMINATION DES VARIABLES CARACTERISANT LE MILIEU
    addeme = mecani(2)
    adcome = mecani(3)
    addep1 = press1(3)
    addlh1 = press1(4)
    adcop1 = press1(7)
    addep2 = press2(3)
    adcop2 = press2(6)
!
! ======================================================================
! --- INITIALISATION DE VECTU ------------------------------------------
! ======================================================================
    do i = 1, dimuel
        vectu(i) = 0.d0
    end do
! ======================================================================
! --- CALCUL POUR CHAQUE POINT DE GAUSS : BOUCLE SUR KPG ---------------
! ======================================================================
    do kpi = 1, npg
!
! ======================================================================
! --- INITIALISATION DE R ----------------------------------------------
! ======================================================================
        do i = 1, dimdef
            r(i) = 0.d0
        end do

! ----- Compute [B] matrix for generalized strains
        call matthm(ds_thm, &
                    ndim, axi, nno1, nno2, dimuel, &
                    dimdef, iu, ip, ipf, iq, &
                    addep1, &
                    addlh1, vff1(1, kpi), vff2(1, kpi), dffr2(1, 1, kpi), wref(kpi), &
                    geom, wi, q)

        call fonoei(ds_thm, &
                    ndim, dimdef, dimcon, &
                    addeme, &
                    addep1, addep2, addlh1, adcome, &
                    adcop1, adcop2, congem(1, kpi), &
                    r)
!
! ======================================================================
! --- CONTRIBUTION DU POINT D'INTEGRATION KPI AU RESIDU ----------------
! ======================================================================
!
        do i = 1, dimuel
            do n = 1, dimdef
                vectu(i) = vectu(i)+q(n, i)*r(n)*wi
            end do
        end do
    end do
!
end subroutine
