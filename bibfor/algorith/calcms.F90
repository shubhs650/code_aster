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
subroutine calcms(materPara, &
                  nbphas, nbcomm, cpmono, nmat, pgl2, &
                  coeft, nfs, nsg, toutms)
!
    use MaterialPara_type
    implicit none
!
#include "asterc/r8dgrd.h"
#include "asterfort/assert.h"
#include "asterfort/lcmmsg.h"
#include "asterfort/matrot.h"
#include "asterfort/promat.h"
#include "asterfort/utmess.h"
!
    type(Material_Para), intent(in) :: materPara
    integer(kind=8) :: nmat, nbcomm(nmat, 3), nfs, nbphas, nsg
    real(kind=8) :: pgl(3, 3), toutms(nbphas, nfs, nsg, 7), coeft(nmat)
    real(kind=8) :: q(3, 3)
!
!       IN
!         NBCOMM :  NOMBRE DE COEF COEFTIAU PAR FAMILLE
!         CPMONO :  NOMS DES LOIS COEFTIAU PAR FAMILLE
!   PGL2   : MATRICE DE PASSAGE GLOBAL LOCAL DU POLYCRSTAL (a faire)
!     OUT:
!           TOUTMS  :  TOUS LES TENSEURS MS
! INTEGRATION DES LOIS POLYCRISTALLINES PAR UNE METHODE DE RUNGE KUTTA
!
!     CETTE ROUTINE FOURNIT LA DERIVEE DE L ENSEMBLE DES VARIABLES
!     INTERNES DU MODELE
!
!     ----------------------------------------------------------------
    character(len=16) :: nomfam
    character(len=24) :: cpmono(5*nmat+1)
    real(kind=8) :: anglNautPhase(3), pgl1(3, 3), pgl2(3, 3)
    real(kind=8) :: ms(6), ng(3), lg(3)
    integer(kind=8) :: nbfsys, i, ifa, nbsys, is, indori, indcp, ir
    integer(kind=8) :: indpha, iphas
!     ----------------------------------------------------------------
    ir = 0
!         CALCUl DES TENSEURS MS POUR GAGNER DU TEMPS
    do iphas = 1, nbphas
!        INDPHA indice debut phase IPHAS dans NBCOMM
        indpha = nbcomm(1+iphas, 1)
!         recuperer l'orientation de la phase et la proportion
        indori = nbcomm(1+iphas, 3)+1
        anglNautPhase(1) = coeft(indori)*r8dgrd()
        anglNautPhase(2) = coeft(indori+1)*r8dgrd()
        anglNautPhase(3) = coeft(indori+2)*r8dgrd()
        call matrot(anglNautPhase, pgl1)
        call matrot(materPara%lcsPara%lcsAngle, pgl2)
        call promat(pgl1, 3, 3, 3, pgl2, &
                    3, 3, 3, pgl)
        nbfsys = nbcomm(indpha, 1)
        indcp = nbcomm(1+iphas, 2)
        ASSERT(nbfsys .le. nfs)
!        Nombre de variables internes de la phase (=monocristal)
        do ifa = 1, nbfsys
            nomfam = cpmono(indcp+5*(ifa-1)+1) (1:16)
            call lcmmsg(nomfam, nbsys, 0, pgl, ms, &
                        ng, lg, ir, q)
            if (nbsys .eq. 0) then
                call utmess('F', 'ALGORITH_70')
            end if
!           indice de la famille IFA
!            INDFA=INDPHA+IFA
!
            do is = 1, nbsys
!              CALCUL DE LA SCISSION REDUITE =
!              PROJECTION DE SIG SUR LE SYSTEME DE GLISSEMENT
!              TAU      : SCISSION REDUITE TAU=SIG:MS
                call lcmmsg(nomfam, nbsys, is, pgl, ms, &
                            ng, lg, ir, q)
                do i = 1, 6
                    toutms(iphas, ifa, is, i) = ms(i)
                end do
            end do
            toutms(iphas, ifa, 1, 7) = nbsys
        end do
!
    end do
end subroutine
