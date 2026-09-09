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
! aslint: disable=W0413
!
subroutine dxmath(plateCara, plateOrie, &
                  famiZ, npg, &
                  df, dm, dmf, &
                  indith)
!
    use plate_type
    implicit none
!
#include "asterc/r8dgrd.h"
#include "asterfort/codent.h"
#include "asterfort/jevech.h"
#include "asterfort/moytem.h"
#include "asterfort/r8inir.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utbtab.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    type(plateCara_Para), intent(in) :: plateCara
    type(plateOrie_Para), intent(in) :: plateOrie
    character(len=*), intent(in) :: famiZ
    integer(kind=8), intent(in) :: npg
    real(kind=8), intent(out) :: df(3, 3), dm(3, 3), dmf(3, 3)
    integer(kind=8), intent(out) :: indith
!
! --------------------------------------------------------------------------------------------------
!
!     CALCUL DES MATRICES DE COEFFCIENTS THERMOELASTIQUES DE FLEXION,
!     MEMBRANE, COUPLAGE MEMBRANE-FLEXION POUR UN MATERIAU ISOTROPE OU
!     MULTICOUCHE
!     OUT MULTIC :
!        1 POUR UN MATERIAU MULTICOUCHE SANS COUPLAGE MEMBRANE-FLEXION
!        2 POUR UN MATERIAU MULTICOUCHE AVEC COUPLAGE MEMBRANE-FLEXION
!        0 DANS LES AUTRES CAS
!     LA VARIABLE INDITH EST INITIALISEE A 0
!     DANS LE CAS OU LE COEFFICIENT DE DILATATION ALPHA N'A
!     PAS ETE DONNE, INDITH VAUT -1 ET ON  NE CALCULE PAS LES
!     CONTRAINTES THERMIQUES
!
! --------------------------------------------------------------------------------------------------
!
    real(kind=8), parameter :: deux = 2.d0
    integer(kind=8), parameter :: npgh = 3
    integer(kind=8), parameter :: nbPara = 1
    character(len=8), parameter:: paraName = 'TEMP'
    real(kind=8) :: paraVale
    integer(kind=8), parameter :: nbPropMaxi = 56
    real(kind=8) :: propVale(nbPropMaxi)
    integer(kind=8) :: propCode(nbPropMaxi)
    character(len=16) :: propName(nbPropMaxi)
    integer(kind=8) :: nbProp
    integer(kind=8) :: jvMaterc, iret, multic
    integer(kind=8) :: i, j, elasco, indalf
    real(kind=8) :: cdf, cdm
    real(kind=8) :: young, nu, epais, excent
    real(kind=8) :: xab1(3, 3)
    real(kind=8) :: alphat
    real(kind=8) :: em, ef, num, nuf
    character(len=3) :: iValStr
    character(len=8) :: fami
    character(len=32) :: elasKeyword
    real(kind=8) :: dmc(3, 2), dfc(3, 2)
!
! --------------------------------------------------------------------------------------------------
!
    fami = famiZ
    dm = 0.d0
    df = 0.d0
    dmf = 0.d0
    indith = 0

! - Get properties of shell
    epais = plateCara%thick
    excent = plateCara%offset

    call jevech('PMATERC', 'L', jvMaterc)
    call rccoma(zi(jvMaterc), 'ELAS', 1, elasKeyword)
    if (elasKeyword .eq. 'ELAS_COQMU') then
        nbProp = 56
        do i = 1, nbProp
            call codent(i, 'G', iValStr)
            propName(i) = 'HOM_'//iValStr
        end do

    else if (elasKeyword .eq. 'ELAS') then
        nbProp = 3
        propName(1) = 'E'
        propName(2) = 'NU'
        propName(3) = 'ALPHA'

    else if (elasKeyword .eq. 'ELAS_GLRC') then
        nbProp = 5
        propName(1) = 'E_M'
        propName(2) = 'NU_M'
        propName(3) = 'E_F'
        propName(4) = 'NU_F'
        propName(5) = 'ALPHA'

    else if (elasKeyword .eq. 'ELAS_COQUE') then
        call rcvalb(fami, 1, 1, '+', &
                    zi(jvMaterc), ' ', elasKeyword, &
                    0, ' ', [0.0d0], &
                    1, 'MEMB_L  ', propVale(1), &
                    propCode, 0)
        if (propCode(1) .eq. 1) then
            call rcvalb(fami, 1, 1, '+', &
                        zi(jvMaterc), ' ', elasKeyword, &
                        0, ' ', [0.0d0], &
                        1, 'M_LLLL  ', propVale(1), &
                        propCode, 0)
            if (propCode(1) .eq. 1) then
                call utmess('F', 'ELEMENTS_41')
            else
                elasco = 2
            end if
        else
            elasco = 1
        end if
        if (elasco .eq. 1) then
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
            propName(11) = 'ALPHA   '
        else if (elasco .eq. 2) then
            nbProp = 33
            multic = 2
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
            propName(34) = 'ALPHA   '
        end if
    else if (elasKeyword .eq. 'ELAS_ORTH') then
        call utmess('F', 'ELEMENTS_91', sk=elasKeyword)
    else if (elasKeyword .eq. 'ELAS_ISTR') then
        call utmess('F', 'ELEMENTS_92', sk=elasKeyword)
    else
        call utmess('F', 'ELEMENTS_42', sk=elasKeyword)
    end if

! - Compute mean temperature (on all point and "sous-point" gauss)
    call moytem(fami, npg, npgh, '+', paraVale, iret)

! - COmpute elasticity matrix
    dmc = 0.d0
    dfc = 0.d0
    if (elasKeyword .eq. 'ELAS') then
        multic = 0
        call rcvalb(fami, 1, 1, '+', &
                    zi(jvMaterc), ' ', elasKeyword, &
                    nbPara, paraName, [paraVale], &
                    2, propName, propVale, &
                    propCode, 1)
        call rcvalb(fami, 1, 1, '+', &
                    zi(jvMaterc), ' ', elasKeyword, &
                    nbPara, paraName, [paraVale], &
                    1, propName(3), propVale(3), &
                    propCode(3), 0)
        if ((propCode(3) .ne. 0) .or. (propVale(3) .eq. 0.d0)) then
            indith = -1
            goto 90
        end if
        young = propVale(1)
        nu = propVale(2)
        alphat = propVale(3)
        young = young*alphat
!
!      ---- CALCUL DE LA MATRICE DE RIGIDITE EN FLEXION --------------
        cdf = young*epais*epais*epais/12.d0/(1.d0-nu*nu)
        df(1, 1) = cdf
        df(1, 2) = cdf*nu
        df(2, 1) = df(1, 2)
        df(2, 2) = df(1, 1)
!      ---- CALCUL DE LA MATRICE DE RIGIDITE EN MEMBRANE -------------
        cdm = epais*young/(1.d0-nu*nu)
        dm(1, 1) = cdm
        dm(1, 2) = cdm*nu
        dm(2, 1) = dm(1, 2)
        dm(2, 2) = dm(1, 1)
!      --- CALCUL DE LA MATRICE DE COUPLAGE MEMBRANE-FLEXION --------
!      --- ET REACTUALISATION DE LA MATRICE DE FLEXION       --------
!      --- DANS LE CAS D'UN EXCENTREMENT                     --------
        do i = 1, 3
            do j = 1, 3
                dmf(i, j) = excent*dm(i, j)
                df(i, j) = df(i, j)+excent*excent*dm(i, j)
            end do
        end do
    else if (elasKeyword .eq. 'ELAS_GLRC') then
!        ------ MATERIAU ISOTROPE ------------------------------------
!
        multic = 0
!
        call rcvalb(fami, 1, 1, '+', &
                    zi(jvMaterc), ' ', elasKeyword, &
                    nbPara, paraName, [paraVale], &
                    2, propName, propVale, &
                    propCode, 1)
!
        em = propVale(1)
        num = propVale(2)
!
        call rcvalb(fami, 1, 1, '+', &
                    zi(jvMaterc), ' ', elasKeyword, &
                    nbPara, paraName, [paraVale], &
                    3, propName(3), propVale(3), &
                    propCode(3), 0)
        if ((propCode(5) .ne. 0) .or. (propVale(5) .eq. 0.d0)) then
            indith = -1
            goto 90
        end if
!
        if (propCode(3) .eq. 0) then
            ef = propVale(3)
        else
            ef = em
        end if
!
        if (propCode(4) .eq. 0) then
            nuf = propVale(4)
        else
            nuf = num
        end if
!
        alphat = propVale(5)
        em = em*alphat
        ef = ef*alphat
!
!      ---- CALCUL DE LA MATRICE DE RIGIDITE EN FLEXION --------------
        cdf = ef*epais*epais*epais/12.d0/(1.d0-nuf*nuf)
        df(1, 1) = cdf
        df(1, 2) = cdf*nuf
        df(2, 1) = df(1, 2)
        df(2, 2) = df(1, 1)
!      ---- CALCUL DE LA MATRICE DE RIGIDITE EN MEMBRANE -------------
        cdm = em*em/(1.d0-num*num)
        dm(1, 1) = cdm
        dm(1, 2) = cdm*num
        dm(2, 1) = dm(1, 2)
        dm(2, 2) = dm(1, 1)
!      --- CALCUL DE LA MATRICE DE COUPLAGE MEMBRANE-FLEXION --------
!      --- ET REACTUALISATION DE LA MATRICE DE FLEXION       --------
!      --- DANS LE CAS D'UN EXCENTREMENT                     --------
        do i = 1, 3
            do j = 1, 3
                dmf(i, j) = excent*dm(i, j)
                df(i, j) = df(i, j)+excent*excent*dm(i, j)
            end do
        end do

    else if (elasKeyword .eq. 'ELAS_COQUE') then
        multic = 0
        call rcvalb(fami, 1, 1, '+', &
                    zi(jvMaterc), ' ', elasKeyword, &
                    nbPara, paraName, [paraVale], &
                    nbProp, propName, propVale, &
                    propCode, 1)
        if (elasco .eq. 1) then
            indalf = 11
        else if (elasco .eq. 2) then
            indalf = 34
        end if
        call rcvalb(fami, 1, 1, '+', &
                    zi(jvMaterc), ' ', elasKeyword, &
                    nbPara, paraName, [paraVale], &
                    1, propName(indalf), propVale(indalf), &
                    propCode(indalf), 0)
        if ((propCode(indalf) .ne. 0) .or. (propVale(indalf) .eq. 0.d0)) then
            indith = -1
            goto 90
        end if
        alphat = propVale(indalf)
!
        if (elasco .eq. 1) then
!        ---- CALCUL DE LA MATRICE DE RIGIDITE EN MEMBRANE -------------
            dm(1, 1) = propVale(1)*alphat
            dm(1, 2) = propVale(2)*alphat
            dm(2, 1) = dm(1, 2)
            dm(2, 2) = propVale(3)*alphat
!        ---- CALCUL DE LA MATRICE DE RIGIDITE EN FLEXION --------------
            df(1, 1) = propVale(5)*alphat
            df(1, 2) = propVale(6)*alphat
            df(2, 1) = df(1, 2)
            df(2, 2) = propVale(7)*alphat
!
        else if (elasco .eq. 2) then
!
            multic = 2
!        ---- CALCUL DE LA MATRICE DE RIGIDITE EN MEMBRANE -------------
            dm(1, 1) = propVale(1)*alphat
            dm(1, 2) = propVale(2)*alphat
            dm(1, 3) = propVale(3)*alphat
            dm(2, 1) = dm(1, 2)
            dm(3, 1) = dm(1, 3)
            dm(2, 2) = propVale(4)*alphat
            dm(2, 3) = propVale(5)*alphat
            dm(3, 3) = propVale(6)*alphat
!        ---- CALCUL DE LA MATRICE DE RIGIDITE EN FLEXION --------------
            df(1, 1) = propVale(7)*alphat
            df(1, 2) = propVale(8)*alphat
            df(1, 3) = propVale(9)*alphat
            df(2, 1) = df(1, 2)
            df(3, 1) = df(1, 3)
            df(2, 2) = propVale(10)*alphat
            df(2, 3) = propVale(11)*alphat
            df(3, 2) = df(2, 3)
            df(3, 3) = propVale(12)*alphat
!        --- COUPLAGE  MEMBRANE FLEXION --------------------------------
            dmf(1, 1) = propVale(13)*alphat
            dmf(1, 2) = propVale(14)*alphat
            dmf(1, 3) = propVale(15)*alphat
            dmf(2, 1) = dmf(1, 2)
            dmf(3, 1) = dmf(1, 3)
            dmf(2, 2) = propVale(16)*alphat
            dmf(2, 3) = propVale(17)*alphat
            dmf(3, 2) = dmf(2, 3)
            dmf(3, 3) = propVale(18)*alphat
!
        end if
!        --- CALCUL DE LA MATRICE DE COUPLAGE MEMBRANE-FLEXION --------
!        --- REACTUALISATION DE LA MATRICE DE FLEXION          --------
!        --- DANS LE CAS D'UN EXCENTREMENT                     --------
        do i = 1, 3
            do j = 1, 3
                df(i, j) = df(i, j)+deux*excent*dmf(i, j)+excent*excent*dm(i, j)
                dmf(i, j) = dmf(i, j)+excent*dm(i, j)
            end do
        end do
!        ----------- MATRICES DANS LE REPERE INTRINSEQUE DE L'ELEMENT --
!
        call utbtab('ZERO', 3, 3, dm, plateOrie%t1ve, xab1, dm)
        call utbtab('ZERO', 3, 3, df, plateOrie%t1ve, xab1, df)
        call utbtab('ZERO', 3, 3, dmf, plateOrie%t1ve, xab1, dmf)
!
    else if (elasKeyword .eq. 'ELAS_COQMU') then
!        ------ MATERIAU MULTICOUCHE -----------------------------------
        call rcvalb(fami, 1, 1, '+', &
                    zi(jvMaterc), ' ', elasKeyword, &
                    nbPara, paraName, [paraVale], &
                    1, propName(19), propVale(19), &
                    propCode(19), 0)
        epais = propVale(19)
        call rcvalb(fami, 1, 1, '+', &
                    zi(jvMaterc), ' ', elasKeyword, &
                    nbPara, paraName, [paraVale], &
                    27, propName(30), propVale(30), &
                    propCode(30), 0)
        dm(1, 1) = propVale(30)
        dm(1, 2) = propVale(31)
        dm(1, 3) = propVale(32)
        dm(2, 1) = propVale(33)
        dm(2, 2) = propVale(34)
        dm(2, 3) = propVale(35)
        dm(3, 1) = propVale(36)
        dm(3, 2) = propVale(37)
        dm(3, 3) = propVale(38)
        dmf(1, 1) = propVale(39)
        dmf(1, 2) = propVale(40)
        dmf(1, 3) = propVale(41)
        dmf(2, 1) = propVale(42)
        dmf(2, 2) = propVale(43)
        dmf(2, 3) = propVale(44)
        dmf(3, 1) = propVale(45)
        dmf(3, 2) = propVale(46)
        dmf(3, 3) = propVale(47)
        df(1, 1) = propVale(48)
        df(1, 2) = propVale(49)
        df(1, 3) = propVale(50)
        df(2, 1) = propVale(51)
        df(2, 2) = propVale(52)
        df(2, 3) = propVale(53)
        df(3, 1) = propVale(54)
        df(3, 2) = propVale(55)
        df(3, 3) = propVale(56)
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
!        ----------- MATRICES DANS LE REPERE INTRINSEQUE DE L'ELEMENT --
!
        call utbtab('ZERO', 3, 3, dm, plateOrie%t1ve, xab1, dm)
        call utbtab('ZERO', 3, 3, df, plateOrie%t1ve, xab1, df)
        call utbtab('ZERO', 3, 3, dmf, plateOrie%t1ve, xab1, dmf)
!
        multic = 1
!
    end if
!
    do i = 1, 3
        do j = 1, 3
            if (abs(dmf(i, j)) .gt. 1.d-10) multic = 2
        end do
    end do
!
90  continue
end subroutine
