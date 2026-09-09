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
! aslint: disable=W1501
!
subroutine dxmate(plateCara, plateOrie, &
                  famiZ, df, dm, dmf, dc, &
                  dci, dmc, dfc, &
                  multic, coupmf)
!
    use plate_type
    implicit none
!
#include "asterc/r8prem.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/codent.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/moyte2.h"
#include "asterfort/moytem.h"
#include "asterfort/plate_type.h"
#include "asterfort/rcadlv.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/rcvalt.h"
#include "asterfort/tecach.h"
#include "asterfort/utbtab.h"
#include "asterfort/utdtab.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    type(plateCara_Para), intent(in) :: plateCara
    type(plateOrie_Para), intent(in) :: plateOrie
    character(len=*), intent(in) :: famiZ
    real(kind=8), intent(out) :: df(3, 3), dm(3, 3), dmf(3, 3)
    real(kind=8), intent(out) :: dc(2, 2), dci(2, 2)
    real(kind=8), intent(out) :: dmc(3, 2), dfc(3, 2)
    integer(kind=8), intent(out) :: multic
    aster_logical, intent(out) :: coupmf
!
! --------------------------------------------------------------------------------------------------
!
!     CALCUL DES MATRICES DE RIGIDITE DE FLEXION, MEMBRANE , COUPLAGE
!     MEMBRANE-FLEXION ET CISAILLEMENT POUR UN MATERIAU ISOTROPE OU
!     MULTICOUCHE
!     OUT MULTIC :
!        1 POUR UN MATERIAU MULTICOUCHE SANS COUPLAGE MEMBRANE-FLEXION
!        0 DANS LES AUTRES CAS
!     OUT COUPMF :
!        .TRUE. POUR UN MATERIAU AVEC COUPLAGE MEMBRANE-FLEXION
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: npgh = 3
    real(kind=8), parameter :: deux = 2.d0
    integer(kind=8), parameter :: nbPara = 1
    character(len=8), parameter :: paraName(nbPara) = 'TEMP'
    integer(kind=8), parameter :: nbPropMaxi = 33
    character(len=16) :: propName(nbPropMaxi)
    integer(kind=8) :: propCode(nbPropMaxi)
    real(kind=8) :: propVale(nbPropMaxi)
    integer(kind=8) :: nbProp
    integer(kind=8), parameter :: nbElasCoque = 45
    integer(kind=8) :: elasCoqueCode(nbElasCoque)
    real(kind=8) :: elasCoqueVale(nbElasCoque)
    integer(kind=8) :: i, j, k, elasco
    integer(kind=8) :: npg, nbLayer, iret
    integer(kind=8) :: jvMaterc
    integer(kind=8) ::  propCodeDHRC, jadr, n1
    real(kind=8) :: kcis, cdf, cdm, cdc, gcis
    real(kind=8) :: young, nu, epais, temp, excent
    real(kind=8) :: xab1(3, 3), xab2(2, 2), xab3(3, 2)
    real(kind=8) :: det
    real(kind=8) :: em, ef, num, nuf
    character(len=3) :: nume
    character(len=32) :: elasKeyword
    character(len=8) :: fami
    aster_logical :: lDKTG
!
! --------------------------------------------------------------------------------------------------
!
    fami = famiZ
    elasco = 0
    multic = 0
    coupmf = ASTER_FALSE
    dmf = 0.d0
    dmc = 0.d0
    dfc = 0.d0
    dc = 0.d0
    dci = 0.d0
    lDKTG = plateCara%type .eq. PLATE_DKTG
!
    call elrefe_info(fami=fami, npg=npg)

! - Get properties of shell
    epais = plateCara%thick
    excent = plateCara%offset
    nbLayer = plateCara%nbLayer
    if (plateCara%type .eq. PLATE_DKTG) then
        ASSERT(nbLayer .eq. 1)
    end if

! - Get access to material properties
    call tecach('NNO', 'PMATERC', 'L', iret, iad=jvMaterc)
    if (iret .ne. 0) then
        goto 999
    end if
    call rccoma(zi(jvMaterc), 'ELAS', 1, elasKeyword)

!   -- calcul de nbProp et propName :
!   ----------------------------
    if (elasKeyword .eq. 'ELAS_COQMU') then
        nbProp = 26
        do i = 1, nbProp
            call codent(i, 'G', nume)
            propName(i) = 'HOM_'//nume
        end do

    else if (elasKeyword .eq. 'ELAS') then
        nbProp = 2
        propName(1) = 'E'
        propName(2) = 'NU'

    else if (elasKeyword .eq. 'ELAS_GLRC') then
        nbProp = 6
        propName(1) = 'E_M'
        propName(2) = 'NU_M'
        propName(3) = 'E_F'
        propName(4) = 'NU_F'
        propName(5) = 'BT1'
        propName(6) = 'BT2'

    else if (elasKeyword .eq. 'ELAS_COQUE') then
!        call utmess('A', 'ELEMENTS_93', sk=elasKeyword)
!       -- on remplit propName plus tard ...

    else if (elasKeyword .eq. 'ELAS_DHRC') then
!      -- pour ELAS_DHRC, on n'utilise pas propName

    else if (elasKeyword .eq. 'ELAS_ORTH') then
        call utmess('F', 'ELEMENTS_91', sk=elasKeyword)

    else if (elasKeyword .eq. 'ELAS_ISTR') then
        call utmess('F', 'ELEMENTS_92', sk=elasKeyword)

    else
        call utmess('F', 'ELEMENTS_42', sk=elasKeyword)
    end if

! - Compute mean temperature
    if (lDKTG) then
        call moyte2(fami, npg, '+', temp, iret)
    else
        call moytem(fami, npg, npgh*nbLayer, '+', temp, iret)
    end if

    if (elasKeyword .eq. 'ELAS') then
!
!        ------ MATERIAU ISOTROPE AVEC DECOUPLAGE MEMBRANE FLEXION----
!
        call rcvalb(fami, 1, 1, '+', zi(jvMaterc), &
                    ' ', elasKeyword, nbPara, paraName, [temp], &
                    nbProp, propName, propVale, propCode, 1)
!
        young = propVale(1)
        nu = propVale(2)
!
        multic = 0
        kcis = 5.d0/6.d0
!
!      ---- CALCUL DE LA MATRICE DE RIGIDITE EN FLEXION --------------
        cdf = young*epais*epais*epais/12.d0/(1.d0-nu*nu)
        df(:, :) = 0.d0
        dmf(:, :) = 0.d0
        df(1, 1) = cdf
        df(1, 2) = cdf*nu
        df(2, 1) = df(1, 2)
        df(2, 2) = df(1, 1)
        df(3, 3) = cdf*(1.d0-nu)/2.d0
!      ---- CALCUL DE LA MATRICE DE RIGIDITE EN MEMBRANE -------------
        cdm = epais*young/(1.d0-nu*nu)
        dm(:, :) = 0.d0
        dm(1, 1) = cdm
        dm(1, 2) = cdm*nu
        dm(2, 1) = dm(1, 2)
        dm(2, 2) = dm(1, 1)
        dm(3, 3) = cdm*(1.d0-nu)/2.d0
!      --- CALCUL DE LA MATRICE DE RIGIDITE EN CISAILLEMENT ----------
        gcis = young/2.d0/(1.d0+nu)
        cdc = gcis*kcis*epais
        dc(1, 1) = cdc
        dc(2, 2) = dc(1, 1)
        dc(1, 2) = 0.d0
        dc(2, 1) = 0.d0
!      --- CALCUL DE SON INVERSE ------------------------------------
        dci(1, 1) = 1.d0/dc(1, 1)
        dci(2, 2) = dci(1, 1)
        dci(1, 2) = 0.d0
        dci(2, 1) = 0.d0
!      --- CALCUL DE LA MATRICE DE COUPLAGE MEMBRANE-FLEXION --------
!      --- ET REACTUALISATION DE LA MATRICE DE FLEXION       --------
!      --- DANS LE CAS D'UN EXCENTREMENT                     --------
        do i = 1, 3
            do j = 1, 3
                dmf(i, j) = excent*dm(i, j)
                df(i, j) = df(i, j)+excent*excent*dm(i, j)
            end do
        end do
!
    else if (elasKeyword .eq. 'ELAS_GLRC') then
!
!        ------ MATERIAU ISOTROPE --------------------------------------
!
        call rcvalb(fami, 1, 1, '+', zi(jvMaterc), &
                    ' ', elasKeyword, nbPara, paraName, [temp], &
                    2, propName, propVale, propCode, 1)
!
        em = propVale(1)
        num = propVale(2)
!
        call rcvalb(fami, 1, 1, '+', zi(jvMaterc), &
                    ' ', elasKeyword, nbPara, paraName, [temp], &
                    2, propName(3), propVale, propCode, 0)
!
        if (propCode(1) .eq. 0) then
            ef = propVale(1)
        else
            ef = em
        end if
!
        if (propCode(2) .eq. 0) then
            nuf = propVale(2)
        else
            nuf = num
        end if
!
        multic = 0
        kcis = 5.d0/6.d0
!
!      ---- CALCUL DE LA MATRICE DE RIGIDITE EN FLEXION --------------
        cdf = ef*epais*epais*epais/12.d0/(1.d0-nuf*nuf)
        df(:, :) = 0.d0
        dmf(:, :) = 0.d0
        df(1, 1) = cdf
        df(1, 2) = cdf*nuf
        df(2, 1) = df(1, 2)
        df(2, 2) = df(1, 1)
        df(3, 3) = cdf*(1.d0-nuf)/2.d0
!      ---- CALCUL DE LA MATRICE DE RIGIDITE EN MEMBRANE -------------
        cdm = epais*em/(1.d0-num*num)
        dm(:, :) = 0.d0
        dm(1, 1) = cdm
        dm(1, 2) = cdm*num
        dm(2, 1) = dm(1, 2)
        dm(2, 2) = dm(1, 1)
        dm(3, 3) = cdm*(1.d0-num)/2.d0
!      --- CALCUL DE LA MATRICE DE RIGIDITE EN CISAILLEMENT ----------
        call rcvalb(fami, 1, 1, '+', zi(jvMaterc), &
                    ' ', elasKeyword, nbPara, paraName, [temp], &
                    2, propName(5), propVale, propCode, 0)
!
        if (propCode(1) .eq. 0) then
            dc(1, 1) = propVale(1)
            dc(2, 2) = propVale(2)
        else
            gcis = em/2.d0/(1.d0+num)
            cdc = gcis*kcis*epais
            dc(1, 1) = cdc
            dc(2, 2) = dc(1, 1)
        end if
!
        dc(1, 2) = 0.d0
        dc(2, 1) = 0.d0
!      --- CALCUL DE SON INVERSE ------------------------------------
        dci(1, 1) = 1.d0/dc(1, 1)
        dci(2, 2) = dci(1, 1)
        dci(1, 2) = 0.d0
        dci(2, 1) = 0.d0
!      --- CALCUL DE LA MATRICE DE COUPLAGE MEMBRANE-FLEXION --------
!      --- ET REACTUALISATION DE LA MATRICE DE FLEXION       --------
!      --- DANS LE CAS D'UN EXCENTREMENT                     --------
        do i = 1, 3
            do j = 1, 3
                dmf(i, j) = excent*dm(i, j)
                df(i, j) = df(i, j)+excent*excent*dm(i, j)
            end do
        end do
!
    else if (elasKeyword .eq. 'ELAS_COQUE') then

!       -- on recupere TOUS les parametres de ELAS_COQUE :
        call rcvalt(fami, 1, 1, '+', zi(jvMaterc), ' ', &
                    'ELAS_COQUE', nbPara, paraName, [temp], &
                    nbElasCoque, elasCoqueVale, elasCoqueCode, 1)

!       -- selon le type d'elasticite :
        if (elasCoqueCode(1) .eq. 0) then
            elasco = 1
            nbProp = 10
            propName(1) = 'MEMB_L  '
            propName(2) = 'MEMB_LT '
            propName(3) = 'MEMB_T  '
            propName(4) = 'MEMB_G_LT'
            propName(5) = 'FLEX_L  '
            propName(6) = 'FLEX_LT '
            propName(7) = 'FLEX_T  '
            propName(8) = 'FLEX_G_LT'
            propName(9) = 'CISA_L  '
            propName(10) = 'CISA_T  '
            do k = 1, nbProp
                ASSERT(elasCoqueCode(k) .eq. 0)
                propVale(k) = elasCoqueVale(k)
            end do
        else
            ASSERT(elasCoqueCode(11) .eq. 0)
            elasco = 2
            nbProp = 33
            coupmf = .true.
            propName(1) = 'M_LLLL  '
            propName(2) = 'M_LLTT  '
            propName(3) = 'M_LLLT  '
            propName(4) = 'M_TTTT  '
            propName(5) = 'M_TTLT  '
            propName(6) = 'M_LTLT  '
            propName(7) = 'F_LLLL  '
            propName(8) = 'F_LLTT  '
            propName(9) = 'F_LLLT  '
            propName(10) = 'F_TTTT  '
            propName(11) = 'F_TTLT  '
            propName(12) = 'F_LTLT  '
            propName(13) = 'MF_LLLL '
            propName(14) = 'MF_LLTT '
            propName(15) = 'MF_LLLT '
            propName(16) = 'MF_TTTT '
            propName(17) = 'MF_TTLT '
            propName(18) = 'MF_LTLT '
            propName(19) = 'MC_LLLZ '
            propName(20) = 'MC_LLTZ '
            propName(21) = 'MC_TTLZ '
            propName(22) = 'MC_TTTZ '
            propName(23) = 'MC_LTLZ '
            propName(24) = 'MC_LTTZ '
            propName(25) = 'FC_LLLZ '
            propName(26) = 'FC_LLTZ '
            propName(27) = 'FC_TTLZ '
            propName(28) = 'FC_TTTZ '
            propName(29) = 'FC_LTLZ '
            propName(30) = 'FC_LTTZ '
            propName(31) = 'C_LZLZ  '
            propName(32) = 'C_LZTZ  '
            propName(33) = 'C_TZTZ  '
            do k = 1, nbProp
                ASSERT(elasCoqueCode(10+k) .eq. 0)
                propVale(k) = elasCoqueVale(10+k)
            end do
        end if

        if (elasco .eq. 1) then
            multic = 0
!
!        ---- CALCUL DE LA MATRICE DE RIGIDITE EN MEMBRANE -------------
            dm(1, 1) = propVale(1)
            dm(1, 2) = propVale(2)
            dm(1, 3) = 0.d0
            dm(2, 1) = dm(1, 2)
            dm(2, 2) = propVale(3)
            dm(2, 3) = 0.d0
            dm(3, 1) = 0.d0
            dm(3, 2) = 0.d0
            dm(3, 3) = propVale(4)
!        ---- CALCUL DE LA MATRICE DE RIGIDITE EN FLEXION --------------
            df(1, 1) = propVale(5)
            df(1, 2) = propVale(6)
            df(1, 3) = 0.d0
            df(2, 1) = df(1, 2)
            df(2, 2) = propVale(7)
            df(2, 3) = 0.d0
            df(3, 1) = 0.d0
            df(3, 2) = 0.d0
            df(3, 3) = propVale(8)
!        --- COUPLAGE  MEMBRANE FLEXION --------------------------------
            dmf(1, 1) = 0.d0
            dmf(1, 2) = 0.d0
            dmf(1, 3) = 0.d0
            dmf(2, 1) = 0.d0
            dmf(2, 2) = 0.d0
            dmf(2, 3) = 0.d0
            dmf(3, 1) = 0.d0
            dmf(3, 2) = 0.d0
            dmf(3, 3) = 0.d0
!        --- CALCUL DE LA MATRICE DE RIGIDITE EN CISAILLEMENT ----------
            dc(1, 1) = propVale(9)
            dc(1, 2) = 0.d0
            dc(2, 1) = 0.d0
            dc(2, 2) = propVale(10)
!        --- CALCUL DE SON INVERSE -------------------------------------
            dci(1, 1) = 1/propVale(9)
            dci(1, 2) = 0.d0
            dci(2, 1) = 0.d0
            dci(2, 2) = 1/propVale(10)
!
        else if (elasco .eq. 2) then
            multic = 0
!        ---- CALCUL DE LA MATRICE DE RIGIDITE EN MEMBRANE -------------
            dm(1, 1) = propVale(1)
            dm(1, 2) = propVale(2)
            dm(1, 3) = propVale(3)
            dm(2, 1) = dm(1, 2)
            dm(2, 2) = propVale(4)
            dm(2, 3) = propVale(5)
            dm(3, 1) = dm(1, 3)
            dm(3, 2) = dm(2, 3)
            dm(3, 3) = propVale(6)
!        ---- CALCUL DE LA MATRICE DE RIGIDITE EN FLEXION --------------
            df(1, 1) = propVale(7)
            df(1, 2) = propVale(8)
            df(1, 3) = propVale(9)
            df(2, 1) = df(1, 2)
            df(2, 2) = propVale(10)
            df(2, 3) = propVale(11)
            df(3, 1) = df(1, 3)
            df(3, 2) = df(2, 3)
            df(3, 3) = propVale(12)
!        --- COUPLAGE  MEMBRANE FLEXION --------------------------------
            dmf(1, 1) = propVale(13)
            dmf(1, 2) = propVale(14)
            dmf(1, 3) = propVale(15)
            dmf(2, 1) = dmf(1, 2)
            dmf(2, 2) = propVale(16)
            dmf(2, 3) = propVale(17)
            dmf(3, 1) = dmf(1, 3)
            dmf(3, 2) = dmf(2, 3)
            dmf(3, 3) = propVale(18)
!        --- COUPLAGE  MEMBRANE CISAILLEMENT ---------------------------
            dmc(1, 1) = propVale(19)
            dmc(1, 2) = propVale(20)
            dmc(2, 1) = propVale(21)
            dmc(2, 2) = propVale(22)
            dmc(3, 1) = propVale(23)
            dmc(3, 2) = propVale(24)
!        --- COUPLAGE  FLEXION CISAILLEMENT ---------------------------
            dfc(1, 1) = propVale(25)
            dfc(1, 2) = propVale(26)
            dfc(2, 1) = propVale(27)
            dfc(2, 2) = propVale(28)
            dfc(3, 1) = propVale(29)
            dfc(3, 2) = propVale(30)
!        --- CALCUL DE LA MATRICE DE RIGIDITE EN CISAILLEMENT ----------
            dc(1, 1) = propVale(31)
            dc(1, 2) = propVale(32)
            dc(2, 1) = dc(1, 2)
            dc(2, 2) = propVale(33)
!        --- CALCUL DE SON INVERSE -------------------------------------
            det = dc(1, 1)*dc(2, 2)-dc(1, 2)*dc(2, 1)
            if (det .gt. r8prem()) then
                dci(1, 1) = dc(2, 2)/det
                dci(1, 2) = -dc(1, 2)/det
                dci(2, 1) = -dc(2, 1)/det
                dci(2, 2) = dc(1, 1)/det
            else
                call utmess('F', 'ELEMENTS_43')
            end if
        end if
!        --- CALCUL DE LA MATRICE DE COUPLAGE MEMBRANE-FLEXION --------
!        --- REACTUALISATION DE LA MATRICE DE FLEXION ET DE LA --------
!        --- MATRICE DE COUPLAGE FLEXION-CISAILLEMENT          --------
!        --- DANS LE CAS D'UN EXCENTREMENT                     --------
        do i = 1, 3
            do j = 1, 3
                df(i, j) = df(i, j)+deux*excent*dmf(i, j)+excent*excent*dm(i, j)
                dmf(i, j) = dmf(i, j)+excent*dm(i, j)
            end do
        end do
!
        do i = 1, 3
            do j = 1, 2
                dfc(i, j) = dfc(i, j)+excent*dmc(i, j)
            end do
        end do
!
        if (.not. lDKTG) then
            call utbtab('ZERO', 3, 3, dm, plateOrie%t1ve, xab1, dm)
            call utbtab('ZERO', 3, 3, df, plateOrie%t1ve, xab1, df)
            call utbtab('ZERO', 3, 3, dmf, plateOrie%t1ve, xab1, dmf)
            call utbtab('ZERO', 2, 2, dc, plateOrie%t2ui, xab2, dc)
            call utbtab('ZERO', 2, 2, dci, plateOrie%t2ui, xab2, dci)
            if (elasco .eq. 2) then
                call utdtab('ZERO', 3, 2, 2, 3, &
                            dmc, plateOrie%t2ui, plateOrie%t1ve, xab3, dmc)
                call utdtab('ZERO', 3, 2, 2, 3, &
                            dfc, plateOrie%t2ui, plateOrie%t1ve, xab3, dfc)
            end if
        end if
!
    else if (elasKeyword .eq. 'ELAS_COQMU') then
!        ------ MATERIAU MULTICOUCHE -----------------------------------
        call rcvalb(fami, 1, 1, '+', zi(jvMaterc), &
                    ' ', elasKeyword, nbPara, paraName, [temp], &
                    18, propName, propVale, propCode, 1)
        dm(1, 1) = propVale(1)
        dm(1, 2) = propVale(2)
        dm(1, 3) = propVale(3)
        dm(2, 2) = propVale(4)
        dm(2, 3) = propVale(5)
        dm(3, 3) = propVale(6)
        dm(2, 1) = dm(1, 2)
        dm(3, 1) = dm(1, 3)
        dm(3, 2) = dm(2, 3)
        dmf(1, 1) = propVale(7)
        dmf(1, 2) = propVale(8)
        dmf(1, 3) = propVale(9)
        dmf(2, 2) = propVale(10)
        dmf(2, 3) = propVale(11)
        dmf(3, 3) = propVale(12)
        dmf(2, 1) = dmf(1, 2)
        dmf(3, 1) = dmf(1, 3)
        dmf(3, 2) = dmf(2, 3)
        df(1, 1) = propVale(13)
        df(1, 2) = propVale(14)
        df(1, 3) = propVale(15)
        df(2, 2) = propVale(16)
        df(2, 3) = propVale(17)
        df(3, 3) = propVale(18)
        df(2, 1) = df(1, 2)
        df(3, 1) = df(1, 3)
        df(3, 2) = df(2, 3)
        call rcvalb(fami, 1, 1, '+', zi(jvMaterc), &
                    ' ', elasKeyword, nbPara, paraName, [temp], &
                    6, propName(21), propVale(21), propCode(21), 1)
        dci(1, 1) = propVale(21)
        dci(2, 2) = propVale(22)
        dci(1, 2) = propVale(23)
        dci(2, 1) = dci(1, 2)
        dc(1, 1) = propVale(24)
        dc(2, 2) = propVale(25)
        dc(1, 2) = propVale(26)
        dc(2, 1) = dc(1, 2)
!
!        --- CALCUL DE LA MATRICE DE COUPLAGE MEMBRANE-FLEXION --------
!        --- REACTUALISATION DE LA MATRICE DE FLEXION          --------
!        --- DANS LE CAS D'UN EXCENTREMENT                     --------
        do i = 1, 3
            do j = 1, 3
                df(i, j) = df(i, j)+deux*excent*dmf(i, j)+excent*excent*dm(i, j)
                dmf(i, j) = dmf(i, j)+excent*dm(i, j)
            end do
        end do
!
!        ----------- MATRICES DANS LE REPERE INTRINSEQUE DE L'ELEMENT --
        call utbtab('ZERO', 3, 3, dm, plateOrie%t1ve, &
                    xab1, dm)
        call utbtab('ZERO', 3, 3, df, plateOrie%t1ve, &
                    xab1, df)
        call utbtab('ZERO', 3, 3, dmf, plateOrie%t1ve, &
                    xab1, dmf)
        call utbtab('ZERO', 2, 2, dc, plateOrie%t2ui, &
                    xab2, dc)
        call utbtab('ZERO', 2, 2, dci, plateOrie%t2ui, &
                    xab2, dci)
!
        multic = 1
!
    else if (elasKeyword .eq. 'ELAS_DHRC') then
        multic = 0
        coupmf = .true.
!        ------ MATERIAU ELAS_DHRC -----------------------------------
        call rcadlv(fami, 1, 1, '+', zi(jvMaterc), ' ', &
                    'ELAS_DHRC', 'A0', nbPara, paraName, [temp], jadr, n1, propCodeDHRC, 1)
        ASSERT(propCodeDHRC .eq. 0 .and. n1 .eq. 21)

        dm(1, 1) = zr(jadr-1+1)
        dm(1, 2) = zr(jadr-1+2)
        dm(1, 3) = zr(jadr-1+3)
        dm(2, 2) = zr(jadr-1+7)
        dm(2, 3) = zr(jadr-1+8)
        dm(3, 3) = zr(jadr-1+12)
        dm(2, 1) = dm(1, 2)
        dm(3, 1) = dm(1, 3)
        dm(3, 2) = dm(2, 3)
        dmf(1, 1) = zr(jadr-1+4)
        dmf(1, 2) = zr(jadr-1+5)
        dmf(1, 3) = zr(jadr-1+6)
        dmf(2, 1) = zr(jadr-1+9)
        dmf(2, 2) = zr(jadr-1+10)
        dmf(2, 3) = zr(jadr-1+11)
        dmf(3, 1) = zr(jadr-1+13)
        dmf(3, 2) = zr(jadr-1+14)
        dmf(3, 3) = zr(jadr-1+15)
        df(1, 1) = zr(jadr-1+16)
        df(1, 2) = zr(jadr-1+17)
        df(1, 3) = zr(jadr-1+18)
        df(2, 2) = zr(jadr-1+19)
        df(2, 3) = zr(jadr-1+20)
        df(3, 3) = zr(jadr-1+21)
        df(2, 1) = df(1, 2)
        df(3, 1) = df(1, 3)
        df(3, 2) = df(2, 3)
!
!        --- CALCUL DE LA MATRICE DE COUPLAGE MEMBRANE-FLEXION --------
!        --- REACTUALISATION DE LA MATRICE DE FLEXION          --------
!        --- DANS LE CAS D'UN EXCENTREMENT                     --------
        do i = 1, 3
            do j = 1, 3
                df(i, j) = df(i, j)+deux*excent*dmf(i, j)+excent*excent*dm(i, j)
                dmf(i, j) = dmf(i, j)+excent*dm(i, j)
            end do
        end do
!
!        ----------- MATRICES DANS LE REPERE INTRINSEQUE DE L'ELEMENT --
        call utbtab('ZERO', 3, 3, dm, plateOrie%t1ve, &
                    xab1, dm)
        call utbtab('ZERO', 3, 3, df, plateOrie%t1ve, &
                    xab1, df)
        call utbtab('ZERO', 3, 3, dmf, plateOrie%t1ve, &
                    xab1, dmf)
    end if
!
    do i = 1, 3
        do j = 1, 3
            if (abs(dmf(i, j)) .gt. 1.d-10) coupmf = .true.
        end do
    end do
!
999 continue
!
end subroutine
