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
subroutine aseihm(ds_thm, option, &
                  lSigm, lVari, lMatr, lVect, &
                  l_axi, ndim, nno1, nno2, &
                  npi, npg, dimuel, dimdef, dimcon, &
                  nbvari, iu, ip, ipf, &
                  iq, mecani, press1, press2, tempe, &
                  vff1, vff2, dffr2, time_prev, time_curr, &
                  deplm, deplp, sigm, sigp, varim, &
                  varip, wref, geom, &
                  compor, vectu, matuu, &
                  retcom)
!
    use Behaviour_module
    use Behaviour_type
    use THM_type
    implicit none
!
#include "asterf_types.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/coeihm.h"
#include "asterfort/matthm.h"
#include "asterfort/thmGetBehaviour.h"
#include "asterfort/thmGetBehaviourChck.h"
#include "asterfort/thmGetBehaviourVari.h"
#include "asterfort/thmGetParaInit.h"
!
    type(THM_DS), intent(inout) :: ds_thm
    aster_logical, intent(in) :: lSigm, lVari, lMatr, lVect
    integer(kind=8) :: ndim, nno1, nno2, npi, npg, dimuel, dimdef, dimcon, nbvari
    integer(kind=8) :: mecani(8), press1(9), press2(9), tempe(5)
    integer(kind=8) :: iu(3, 18), ip(2, 9), ipf(2, 2, 9), iq(2, 2, 9)
    real(kind=8) :: vff1(nno1, npi), vff2(nno2, npi), dffr2(ndim-1, nno2, npi)
    real(kind=8) :: wref(npi)
    real(kind=8) :: time_prev, time_curr, deplm(dimuel), deplp(dimuel)
    real(kind=8) :: geom(ndim, nno2)
    real(kind=8) :: sigm(dimcon, npi), varim(nbvari, npi)
    character(len=16) :: option
    aster_logical :: l_axi
    character(len=16), intent(in) :: compor(COMPOR_SIZE)
    integer(kind=8) :: retcom
    real(kind=8) :: vectu(dimuel), varip(nbvari, npi), sigp(dimcon, npi)
    real(kind=8) :: matuu(dimuel*dimuel)
!
! --------------------------------------------------------------------------------------------------
!
!     BUT:  CALCUL DU VECTEUR FORCES INTERNES ELEMENTAIRE, DES
!           CONTRAINTES GENERALISEES, DES VARIABLES INTERNES
!           ET/OU DE L'OPERATEUR TANGENT ELEMENTAIRE
!
! --------------------------------------------------------------------------------------------------
!
! IO  ds_thm           : datastructure for THM
! IN OPTION  : OPTION DE CALCUL
! IN AXI     : AXISYMETRIQUE ?
! IN NDIM    : DIMENSION DE L'ESPACE
! IN NNO1    : NOMBRE DE NOEUDS DE LA FAMILLE 1
! IN NNO2    : NOMBRE DE NOEUDS DE LA FAMILLE 2
! IN NPI     : NOMBRE DE POINTS D'INTEGRATION
! IN NPG     : NOMBRE DE POINTS DE GAUSS
! IN DIMUEL  : NOMBRE DE DDL
! IN DIMDEF  : DIMENSION DU VECTEUR DEFORMATIONS GENERALISEES
! IN DIMCON  : DIMENSION DU VECTEUR CONTRAINTES GENERALISEES
! IN IMATE   : MATERIAU CODE
! IN IU      : DECALAGE D'INDICE POUR ACCEDER AUX DDL DE DEPLACEMENT
! IN IP      : DECALAGE D'INDICE POUR ACCEDER AUX DDL DE PRESSION MILIEU
! IN IPF     : DECALAGE D'INDICE POUR ACCEDER AUX DDL DE PRESSION FACES
! IN IQ      : DECALAGE D'INDICE POUR ACCEDER AUX DDL DE LAGRANGE HYDRO
! IN MECANI  : INFOS MECANIQUE
! IN PRESS1  : INFOS CONSTITUANT 1
! IN PRESS2  : INFOS CONSTITUANT 2
! IN TEMPE   : INFOS TEMPERATURE
! IN VFF1    : VALEUR DES FONCTIONS DE FORME (FAMILLE 1)
! IN VFF2    : VALEUR DES FONCTIONS DE FORME (FAMILLE 2)
! IN DFFR2   : DERIVEES DES FONCTIONS DE FORME (FAMILLE 2)
! IN INSTAM  : INSTANT PRECEDENT
! IN INSTAP  : INSTANT ACTUEL
! IN DEPLM   : DDL A L'INSTANT PRECEDENT
! IN DEPLP   : DDL A L'INSTANT ACTUEL
! IN SIGM    : CONTRAINTES GENERALISEES AUX POINTS D'INTEGRATION
! IN VARIM   : VARIABLES INTERNES AUX POINTS D'INTEGRATION
! IN NOMAIL  : NUMERO DE LA MAILLE
! IN WREF    : POIDS DE REFERENCE DES POINTS D'INTEGRATIONS
! IN GEOM    : COORDONNEES DES NOEUDS (FAMILLE 1)
! IN ANG     : ANGLES D'EULER NODAUX (FAMILLE 1)
! IN COMPOR  : COMPORTEMENT
! IN PERMAN  : REGIME PERMANENT ?
! =====================================================================
! OUT RETCOM : CODE RETOUR LOI DE COMPORTEMENT
! OUT VECTU  : VECTEUR FORCE INTERNE ELEMENTAIRE
! OUT MATUU  : OPERATEUR TANGENT ELEMENTAIRE
! OUT SIGP   : CONTRAINTES GENERALISEES AU TEMPS PLUS AUX POINTS
!              D'INTEGRATION
! OUT VARIP  : VARIABLES INTERNES AU TEMPS PLUS AUX POINTS D'INTEGRATION
! --- POUR L'HYDRAULIQUE : VAR. INT. 1 : RHO_LIQUIDE - RHO_0
! --- POUR LE COUPLAGE   : VAR. INT. 1 : PHI - PHI_0
! ---                    : VAR. INT. 2 : PVP - PVP_0 SI VAPEUR
! ---                    : VAR. INT. 3 : SATURATION SI LOI NON SATUREE
!                        : VAR. INT. 4 : OUVH
! --- POUR LA MECANIQUE  : VAR. INT. 1 : TLINT
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: ksp = 1
    integer(kind=8) :: addeme, adcome, addep1, addep2, addete
    integer(kind=8) :: adcp11, adcp12, adcp21, adcp22, adcote, adcop1
    integer(kind=8) :: i, j, m, k, km, kpi, addlh1
    real(kind=8) :: q(dimdef, dimuel), res(dimdef), drde(dimdef, dimdef), wi
    real(kind=8) :: defgem(dimdef), defgep(dimdef), matri
!
! --------------------------------------------------------------------------------------------------
!
    if (lSigm) then
        sigp = sigm
    end if

! - Get parameters for behaviour
    call thmGetBehaviour(compor, ds_thm)

! - Get parameters for internal variables
    call thmGetBehaviourVari(ds_thm)

! - Some checks between behaviour and model
    call thmGetBehaviourChck(ds_thm)

! - Get initial parameters (THM_INIT)
    call thmGetParaInit(ds_thm, l_check_=ASTER_TRUE)
!
! =====================================================================
! --- DETERMINATION DES VARIABLES CARACTERISANT LE MILIEU ET OPTION ---
! =====================================================================
    addeme = mecani(2)
    adcome = mecani(3)
    addep1 = press1(3)
    addlh1 = press1(4)
    adcp11 = press1(5)
    adcp12 = press1(6)
    adcop1 = press1(7)
    addep2 = press2(3)
    adcp21 = press2(4)
    adcp22 = press2(5)
    addete = tempe(2)
    adcote = tempe(3)
!
! ======================================================================
! --- INITIALISATION DE VECTU, MATUU A 0 SUIVANT OPTION ----------------
! ======================================================================
    if (lVect) then
        vectu(:) = 0.d0
    end if
    if (lMatr) then
        matuu(:) = 0.d0
    end if

! - Loop on integration points
    do kpi = 1, npi

! ----- Compute [B] matrix for generalized strains
        call matthm(ds_thm, &
                    ndim, l_axi, nno1, nno2, dimuel, &
                    dimdef, iu, ip, ipf, iq, &
                    addep1, &
                    addlh1, vff1(1, kpi), vff2(1, kpi), dffr2(1, 1, kpi), wref(kpi), &
                    geom, wi, q)

! ----- Compute generalized strains
        do i = 1, dimdef
            defgem(i) = 0.d0
            defgep(i) = 0.d0
            do j = 1, dimuel
                defgem(i) = defgem(i)+q(i, j)*deplm(j)
                defgep(i) = defgep(i)+q(i, j)*deplp(j)
            end do
        end do

! ----- Set main parameters for behaviour (on point)
        call behaviourSetParaPoin(kpi, ksp, ds_thm%ds_behaviour%BEHInteg)

! ----- Integration of behaviour
        call coeihm(ds_thm, &
                    lSigm, lVari, lMatr, lVect, &
                    option, &
                    time_prev, time_curr, &
                    ndim, dimdef, dimcon, nbvari, &
                    addeme, adcome, &
                    addep1, adcp11, adcp12, addlh1, adcop1, &
                    addep2, adcp21, adcp22, addete, adcote, &
                    defgem, defgep, kpi, npg, npi, &
                    sigm(1, kpi), sigp(1, kpi), varim(1, kpi), varip(1, kpi), res, &
                    drde, retcom)
!
! =====================================================================
! --- CALCUL DES FORCES INTERIEURES ET DE L'OPERATEUR TANGENT ---------
! =====================================================================
!
        if (lVect) then
            do k = 1, dimuel
                do i = 1, dimdef
                    vectu(k) = vectu(k)+wi*q(i, k)*res(i)
                end do
            end do
        end if
!
        if (lMatr) then
            km = 1
            do k = 1, dimuel
                do m = 1, dimuel
                    matri = 0.d0
                    do i = 1, dimdef
                        do j = 1, dimdef
                            matri = matri+wi*q(i, k)*drde(i, j)*q(j, m)
                        end do
                    end do
                    matuu(km) = matuu(km)+matri
                    km = km+1
                end do
            end do
        end if
    end do
!
end subroutine
