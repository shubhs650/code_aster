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
subroutine nmfi2d(BEHInteg, &
                  npg, lgpg, option, geom, &
                  deplm, ddepl, sigmo, sigma, fint, &
                  kTang, vim, vip, &
                  tm, tp, &
                  carcri, compor, typmod, &
                  lMatr, lVect, lSigm, &
                  codret)
!
    use Behaviour_type
    use Behaviour_module
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/codere.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/gedisc.h"
#include "asterfort/nmcomp.h"
#include "asterfort/nmfisa.h"
#include "jeveux.h"
!
    type(Behaviour_Integ), intent(inout) :: BEHInteg
    integer(kind=8) :: npg, lgpg, codret
    real(kind=8) :: geom(2, 4), deplm(8), ddepl(8), tm, tp
    real(kind=8) :: fint(8), kTang(8, 8), sigmo(6, npg), sigma(6, npg)
    real(kind=8) :: vim(lgpg, npg), vip(lgpg, npg)
    character(len=8), intent(in) :: typmod(2)
    character(len=16), intent(in) :: option, compor(COMPOR_SIZE)
    real(kind=8), intent(in) :: carcri(CARCRI_SIZE)
    aster_logical, intent(in) :: lMatr, lVect, lSigm
!
! --------------------------------------------------------------------------------------------------
!
! ELEMENT DE JOINT
!      CALCUL DU SAUT DANS L'ELEMENT
!             DE LA CONTRAINTE A PARTIR D'UNE LDC
!             DE FINT ET KTAN : EFFORTS INTERIEURS ET MATRICE TANGENTE.
!
!      OPTION : OPTIONS DE CALCUL EN FONCTION DE LA SUBROUTINE LANCEE
!       * RAPH_MECA      : U = U- + DU  ->   SIGMA , FINT
!       * FULL_MECA      : U = U- + DU  ->   SIGMA , FINT , KTAN
!       * RIGI_MECA_TANG : U = U-       ->                  KTAN
!       * FORC_NODA      : TRAITE DANS NMFIFI.F
!
! SUBROUTINE APPELEE DANS LE TE0201
!
! IN  : OPTION,COMPOR,GEOM,DEPLM,DDEPL,VIM,NPG,TYPMOD,MATE
! IN  : TM INSTANT MOINS
! IN  : TP INSTANT PLUS
! OUT : SIGMA,FINT,KTAN,VIP,CODRET
! I/O :
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: multComp = " "
    integer(kind=8), parameter :: ksp = 1
    aster_logical :: axi
    integer(kind=8) :: cod(9), i, j, q, s, kpg
    integer(kind=8) :: ndim, nno, ipoids, ivf, idfde
    real(kind=8) :: coopg(3, npg)
    real(kind=8) :: dsidep(6, 6), b(2, 8), sigmPost(6)
    real(kind=8) :: sum(2), dsu(2), poids
!
! --------------------------------------------------------------------------------------------------
!
    cod = 0
    axi = typmod(1) .eq. 'AXIS'
    if (lVect) then
        fint = 0.d0
    end if
    if (lMatr) then
        kTang = 0.d0
    end if

! - Get element parameters
    call elrefe_info(fami=BEHInteg%materPara%schemePara%fami, &
                     ndim=ndim, nno=nno, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde)

!   CALCUL DES COORDONNEES DES POINTS DE GAUSS
    call gedisc(2, nno, npg, zr(ivf), geom, coopg)

! - Loop on Gauss points
    do kpg = 1, npg
!
! CALCUL DE LA MATRICE B DONNANT LES SAUT PAR ELEMENTS A PARTIR DES
! DEPLACEMENTS AUX NOEUDS , AINSI QUE LE POIDS DES PG :
! LE CHANGEMENT DE REPERE EST INTEGRE DANS LA MATRICE B (VOIR NMFISA)
!
        call nmfisa(axi, geom, kpg, poids, b)
!
! CALCUL DU SAUT DE DEPLACEMENT DANS L'ELEMENT (SU_N,SU_T) = B U
!
        sum = 0.d0
        dsu = 0.d0
        do j = 1, 8
            sum(1) = sum(1)+b(1, j)*deplm(j)
            sum(2) = sum(2)+b(2, j)*deplm(j)
        end do
        if (lVect) then
            do j = 1, 8
                dsu(1) = dsu(1)+b(1, j)*ddepl(j)
                dsu(2) = dsu(2)+b(2, j)*ddepl(j)
            end do
        end if

! ----- Set main parameters for behaviour (on point)
        call behaviourSetParaPoin(kpg, ksp, BEHInteg)
        BEHInteg%behavESVA%behavESVAGeom%coorElga(kpg, 1:2) = coopg(1:2, kpg)

! ----- Integrator
        sigmPost = 0.d0
        call nmcomp(BEHInteg, &
                    ndim, option, typmod, &
                    tm, tp, &
                    compor, carcri, multComp, &
                    2, sum, dsu, &
                    2, sigmo(:, kpg), &
                    vim(:, kpg), &
                    sigmPost, vip(:, kpg), &
                    36, dsidep, cod(kpg))
        if (cod(kpg) .eq. 1) goto 900

! ----- Stresses
        if (lSigm) then
            sigma(1, kpg) = sigmPost(1)
            sigma(2, kpg) = sigmPost(2)
        end if

! ----- Internal forces
        if (lVect) then
!       Il faudrait séparer les deux => petit travail de réflexion
            ASSERT(lSigm)
            do i = 1, 8
                do q = 1, 2
                    fint(i) = fint(i)+poids*b(q, i)*sigma(q, kpg)
                end do
            end do
        end if

! ----- Rigidity matrix
        if (lMatr) then
            do i = 1, 8
                do j = 1, 8
                    do q = 1, 2
                        do s = 1, 2
                            kTang(i, j) = kTang(i, j)+poids*b(q, i)*dsidep(q, s)*b(s, j)
                        end do
                    end do
                end do
            end do
        end if
    end do
!
900 continue
    call codere(cod, npg, codret)
!
end subroutine
