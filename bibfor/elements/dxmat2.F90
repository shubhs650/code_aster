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
subroutine dxmat2(plateCara, plateOrie, &
                  iLayer, npg, &
                  ordi, epi, &
                  epais, dm)
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
    integer(kind=8), intent(in) :: iLayer, npg
    real(kind=8), intent(out) :: ordi, epi
    real(kind=8), intent(out) :: epais, dm(3, 3)
!
! --------------------------------------------------------------------------------------------------
!
!       CALCUL DES MATRICES DE COEFFCIENTS THERMOELASTIQUES DE MEMBRANE,
!       POUR UN MATERIAU ISOTROPE OU MULTICOUCHE
!       LA VARIABLE INDITH EST INITIALISEE A 0
!       DANS LE CAS OU LE COEFFICIENT DE DILATATION ALPHA N'A
!       PAS ETE DONNE, INDITH VAUT -1 ET ON  NE CALCULE PAS LES
!       CONTRAINTES THERMIQUES
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbPara = 1
    character(len=8), parameter :: paraName(nbPara) = (/'TEMP'/)
    real(kind=8) :: paraVale
    integer(kind=8), parameter :: nbPropMaxi = 134
    integer(kind=8) :: propCode(nbPropMaxi)
    real(kind=8) :: propVale(nbPropMaxi)
    character(len=16) :: propName(nbPropMaxi)
    integer(kind=8) :: jvMaterc, iret
    integer(kind=8) :: nbProp, i, nbLayer
    real(kind=8) :: cdf, cdm, df(3, 3), dmf(3, 3)
    real(kind=8) :: young, nu, alpha
    real(kind=8) :: xab1(3, 3), dh(3, 3)
    character(len=2) :: iValeStre
    character(len=3) :: iPropStr, iLayerStr
    character(len=32) :: elasKeyword
!
! --------------------------------------------------------------------------------------------------
!
    call r8inir(9, 0.d0, dm, 1)
    call r8inir(9, 0.d0, df, 1)
    call r8inir(9, 0.d0, dh, 1)
    call r8inir(9, 0.d0, dmf, 1)

! - Get parameters of shells
    epais = plateCara%thick
    nbLayer = plateCara%nbLayer

! - We are on current layer
    epi = epais
    ordi = 0.d0

    call jevech('PMATERC', 'L', jvMaterc)
    call rccoma(zi(jvMaterc), 'ELAS', 1, elasKeyword)
!
    if (elasKeyword .eq. 'ELAS_COQMU') then
        nbProp = 56
        do i = 1, nbProp
            call codent(i, 'G', iPropStr)
            propName(i) = 'HOM_'//iPropStr
        end do
        call codent(iLayer, 'G', iLayerStr)
        do i = 1, 78
            call codent(i, 'G', iValeStre)
            propName(56+i) = 'C'//iLayerStr//'_V'//iValeStre
        end do

    else if (elasKeyword .eq. 'ELAS') then
        nbProp = 3
        propName(1) = 'E'
        propName(2) = 'NU'
        propName(3) = 'ALPHA'

    else if (elasKeyword .eq. 'ELAS_COQUE') then
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

    else
        call utmess('F', 'ELEMENTS_44', sk=elasKeyword)
    end if

! - RECUPERATION DE LA TEMPERATURE POUR LE MATERIAU
    call moytem('RIGI', npg, 3*nbLayer, '+', paraVale, iret)

    if (elasKeyword .eq. 'ELAS') then
        call rcvalb('RIGI', 1, 1, '+', &
                    zi(jvMaterc), ' ', elasKeyword, &
                    nbPara, paraName, [paraVale], &
                    2, propName, propVale, &
                    propCode, 1)
        call rcvalb('RIGI', 1, 1, '+', &
                    zi(jvMaterc), ' ', elasKeyword, &
                    nbPara, paraName, [paraVale], &
                    1, propName(3), propVale(3), &
                    propCode(3), 0)
        if ((propCode(3) .ne. 0) .or. (propVale(3) .eq. 0.d0)) then
            goto 70
        else if ((iret .eq. 1) .and. (propCode(3) .ne. 0)) then
            call utmess('F', 'CALCULEL_15')
        end if
        young = propVale(1)
        nu = propVale(2)
        alpha = propVale(3)
        young = young*alpha

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

    else if (elasKeyword .eq. 'ELAS_COQUE') then
        call rcvalb('RIGI', 1, 1, '+', &
                    zi(jvMaterc), ' ', elasKeyword, &
                    nbPara, paraName, [paraVale], &
                    nbProp, propName, propVale, &
                    propCode, 1)
        call rcvalb('RIGI', 1, 1, '+', &
                    zi(jvMaterc), ' ', elasKeyword, &
                    nbPara, paraName, [paraVale], &
                    1, propName(11), propVale(11), &
                    propCode(11), 0)
        if ((propCode(11) .ne. 0) .or. (propVale(11) .eq. 0.d0)) then
            goto 70
        else if ((iret .eq. 1) .and. (propCode(11) .ne. 0)) then
            call utmess('F', 'CALCULEL_15')
        end if
        alpha = propVale(11)

!        ---- CALCUL DE LA MATRICE DE RIGIDITE EN MEMBRANE -------------
        dm(1, 1) = propVale(1)*alpha
        dm(1, 2) = propVale(2)*alpha
        dm(2, 1) = dm(1, 2)
        dm(2, 2) = propVale(3)*alpha

!        ---- CALCUL DE LA MATRICE DE RIGIDITE EN FLEXION --------------
        df(1, 1) = propVale(5)*alpha
        df(1, 2) = propVale(6)*alpha
        df(2, 1) = df(1, 2)
        df(2, 2) = propVale(7)*alpha

!        ----------- MATRICES DANS LE REPERE INTRINSEQUE DE L'ELEMENT --
        call utbtab('ZERO', 3, 3, dm, plateOrie%t1ve, xab1, dm)
        call utbtab('ZERO', 3, 3, df, plateOrie%t1ve, xab1, df)
        call utbtab('ZERO', 3, 3, dmf, plateOrie%t1ve, xab1, dmf)

    else if (elasKeyword .eq. 'ELAS_COQMU') then
        call rcvalb('RIGI', 1, 1, '+', zi(jvMaterc), &
                    ' ', elasKeyword, nbPara, paraName, [paraVale], &
                    1, propName(19), propVale(19), propCode(19), 1)
        epais = propVale(19)
        call rcvalb('RIGI', 1, 1, '+', zi(jvMaterc), &
                    ' ', elasKeyword, nbPara, paraName, [paraVale], &
                    1, propName(57), propVale(57), propCode(57), 1)
        epi = propVale(57)
        call rcvalb('RIGI', 1, 1, '+', zi(jvMaterc), &
                    ' ', elasKeyword, nbPara, paraName, [paraVale], &
                    1, propName(59), propVale(59), propCode(59), 1)
        ordi = propVale(59)
        call rcvalb('RIGI', 1, 1, '+', zi(jvMaterc), &
                    ' ', elasKeyword, nbPara, paraName, [paraVale], &
                    27, propName(102), propVale(102), propCode(102), 1)
        dm(1, 1) = propVale(102)
        dm(1, 2) = propVale(103)
        dm(1, 3) = propVale(104)
        dm(2, 1) = propVale(105)
        dm(2, 2) = propVale(106)
        dm(2, 3) = propVale(107)
        dm(3, 1) = propVale(108)
        dm(3, 2) = propVale(109)
        dm(3, 3) = propVale(110)
        dmf(1, 1) = propVale(111)
        dmf(1, 2) = propVale(112)
        dmf(1, 3) = propVale(113)
        dmf(2, 1) = propVale(114)
        dmf(2, 2) = propVale(115)
        dmf(2, 3) = propVale(116)
        dmf(3, 1) = propVale(117)
        dmf(3, 2) = propVale(118)
        dmf(3, 3) = propVale(119)
        df(1, 1) = propVale(120)
        df(1, 2) = propVale(121)
        df(1, 3) = propVale(122)
        df(2, 1) = propVale(123)
        df(2, 2) = propVale(124)
        df(2, 3) = propVale(125)
        df(3, 1) = propVale(126)
        df(3, 2) = propVale(127)
        df(3, 3) = propVale(128)

! ----- MATRICES DANS LE REPERE INTRINSEQUE DE L'ELEMENT
        call utbtab('ZERO', 3, 3, dm, plateOrie%t1ve, xab1, dm)
        call utbtab('ZERO', 3, 3, df, plateOrie%t1ve, xab1, df)
        call utbtab('ZERO', 3, 3, dmf, plateOrie%t1ve, xab1, dmf)
!
    end if
70  continue
end subroutine
