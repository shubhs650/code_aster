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
subroutine mbxchg(plateOrie, &
                  option, fami, &
                  nddl, nno, ncomp, kpg, npg, &
                  jvEpsi, jvInst, ipoids, jvGeom, &
                  jvMaterc, jvPesa, jvVect, jvSief, &
                  vff, dff)
!
    use plate_type
    use resi_refe_module, only: RESI_REFE
    implicit none
!
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/fointe.h"
#include "asterfort/mbcine.h"
#include "asterfort/mbrigi.h"
#include "asterfort/rcvalb.h"
#include "asterfort/verift.h"
#include "jeveux.h"
!
    type(plateOrie_Para), intent(in) :: plateOrie
    character(len=16), intent(in) :: option
    character(len=8), intent(in) :: fami
    integer(kind=8), intent(in) :: nddl, nno, ncomp, npg
    integer(kind=8), intent(in) :: kpg
    integer(kind=8), intent(in) :: ipoids, jvGeom, jvMaterc, jvPesa, jvEpsi, jvInst
    integer(kind=8), intent(in) :: jvVect, jvSief
    real(kind=8), intent(in) :: dff(2, nno), vff(nno)
!
! --------------------------------------------------------------------------------------------------
!
!    - FONCTION REALISEE:  CALCUL DES OPTIONS DE DE CHARGEMENT :
!                                  - CHAR_MECA_EPSI_R
!                                  - CHAR_MECA_EPSI_F
!                                  - CHAR_MECA_PESA_R
!                                  - CHAR_MECA_TEMP_R
!                                  - FORC_NODA
!                                  - REFE_FORC_NODA
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
! IN  NPG          NOMBRE DE POINT DE GAUSS
! IN  IEPSIN       ADRESSE DANS ZR DU TABLEAU PEPSINR
! IN  ITEMPS       ADRESSE DANS ZR DU TABLEAU PINSTR
! IN  IPOIDS       ADRESSE DANS ZR DU TABLEAU POIDS
! IN  IGEOM        ADRESSE DANS ZR DU TABLEAU PGEOMER
! IN  IMATE        ADRESSE DANS ZI DU TABLEAU PMATERC
! IN  IPESA        ADRESSE DANS ZR DU TABLEAU PPESANR
! IN  IVECTU       ADRESSE DANS ZR DU TABLEAU PVECTUR
! IN  ICONTM       ADRESSE DANS ZR DU TABLEAU PCONMR
! IN  VFF          VALEURS DES FONCTIONS DE FORME
! IN  DFF          DERIVEE DES F. DE FORME
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbPara = 4
    character(len=8), parameter :: paraName(nbPara) = (/'X   ', 'Y   ', 'Z   ', &
                                                        'INST'/)
    real(kind=8) :: paraVale(nbPara)
    integer(kind=8), parameter :: nbProp = 1
    character(len=8), parameter :: propName(nbProp) = (/'RHO'/)
    real(kind=8) :: propVale(nbProp)
    integer(kind=8) :: propCode(nbProp)
    integer(kind=8) :: i, n, c, cc, ier
    real(kind=8) :: b(3, 3, 9), jac
    real(kind=8) :: matrRigi(3, 3), rho
    real(kind=8) :: epsthe, sgmref, sig(3), alpha, beta
    real(kind=8) :: xgau, ygau, zgau, epsinif(3)
    type(RESI_REFE):: refe
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(plateOrie%lUpdate)
    alpha = plateOrie%alpha
    beta = plateOrie%beta

! - CALCUL DE LA MATRICE "B"
    call mbcine(plateOrie, &
                nno, zr(jvGeom), dff, &
                b, jac)

    if ((option .eq. 'FORC_NODA') .or. (option .eq. 'CHAR_MECA_TEMP_R') .or. &
        (option(1:15) .eq. 'CHAR_MECA_EPSI_')) then
        if (option .eq. 'FORC_NODA') then
            do c = 1, ncomp
                sig(c) = zr(jvSief+(kpg-1)*ncomp+c-1)
            end do

        else if (option .eq. 'CHAR_MECA_EPSI_R') then
            call mbrigi(fami, kpg, jvMaterc, matrRigi)
            sig = 0.d0
            do c = 1, ncomp
                do cc = 1, ncomp
                    sig(c) = sig(c)+zr(jvEpsi+ncomp*(kpg-1)+cc-1)*matrRigi(cc, c)
                end do
            end do

        else if (option .eq. 'CHAR_MECA_EPSI_F') then
            call mbrigi(fami, kpg, jvMaterc, matrRigi)
            sig = 0.d0
            paraVale(4) = zr(jvInst)
            xgau = 0.d0
            ygau = 0.d0
            zgau = 0.d0
            do i = 1, nno
                xgau = xgau+vff(i)*zr(jvGeom-1+1+3*(i-1))
                ygau = ygau+vff(i)*zr(jvGeom-1+2+3*(i-1))
                zgau = zgau+vff(i)*zr(jvGeom-1+3+3*(i-1))
            end do
            paraVale(1) = xgau
            paraVale(2) = ygau
            paraVale(3) = zgau
            call fointe('FM', zk8(jvEpsi), 4, paraName, paraVale, epsinif(1), ier)
            call fointe('FM', zk8(jvEpsi+1), 4, paraName, paraVale, epsinif(2), ier)
            call fointe('FM', zk8(jvEpsi+2), 4, paraName, paraVale, epsinif(3), ier)
            do c = 1, ncomp
                do cc = 1, ncomp
                    sig(c) = sig(c)+epsinif(cc)*matrRigi(cc, c)
                end do
            end do

        else if (option .eq. 'CHAR_MECA_TEMP_R') then
            call verift(fami, kpg, 1, '+', zi(jvMaterc), &
                        epsth_=epsthe)
            call mbrigi(fami, kpg, jvMaterc, matrRigi)
            sig = 0.d0
            do c = 1, ncomp
                sig(c) = epsthe*(matrRigi(1, c)+matrRigi(2, c))
            end do
        end if
        do n = 1, nno
            do i = 1, nddl
                do c = 1, ncomp
                    zr(jvVect+(n-1)*nddl+i-1) = zr(jvVect+(n-1)* &
                                                   nddl+i-1)+b(c, i, n)*sig(c)*zr(ipoids+kpg-1)*jac
                end do
            end do
        end do

    else if (option .eq. 'REFE_FORC_NODA') then
        call refe%Init('MEMBRANE')
        sgmref = refe%GetRef('SIGM')
        call refe%Check()
        ASSERT(sgmref .gt. 0.d0)
        do n = 1, nno
            do i = 1, nddl
                zr(jvVect+(n-1)*nddl+i-1) = zr(jvVect+(n-1)*nddl+i-1)+ &
                                            sgmref*sqrt(abs(jac))/npg
            end do
        end do

    else if (option .eq. 'CHAR_MECA_PESA_R') then
        call rcvalb(fami, kpg, 1, '+', &
                    zi(jvMaterc), ' ', 'ELAS_MEMBRANE', &
                    0, ' ', [0.d0], &
                    nbProp, propName, propVale, &
                    propCode, 1)
        rho = propVale(1)
        do n = 1, nno
            do i = 1, nddl
                zr(jvVect+(n-1)*nddl+i-1) = zr(jvVect+(n-1)*nddl+i-1)+ &
                                            rho*zr(jvPesa)*zr(jvPesa+i)*vff(n)*zr(ipoids+kpg-1)*jac
            end do
        end do
    end if

end subroutine
