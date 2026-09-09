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
subroutine te0404(option, nomte)
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystPara, compCoorSystPlate
    implicit none
!
#include "asterf_types.h"
#include "asterfort/dxmate.h"
#include "asterfort/ElasticityMaterial_type.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/get_elas_id.h"
#include "asterfort/jevech.h"
#include "asterfort/rcvalb.h"
#include "asterfort/teattr.h"
#include "asterfort/utmess.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: TOUS
!
! Options: PAS_COURANT
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8) :: fami
    integer(kind=8) :: propCode(2)
    character(len=16) :: propName(2)
    character(len=8) :: cnd
    integer(kind=8) :: icour, jvMaterc, jvGeom, nd, ndim, nno, nnos, npg
    integer(kind=8) :: i, j, ipoids, ivf, idfde, jgano, ier
    integer(kind=8) :: multic
    real(kind=8) :: dmin, distij, xi, yi, zii, xj, yj, zj
    real(kind=8) :: e, nu, vitmat, epais
    real(kind=8) :: df(3, 3), dm(3, 3), dmf(3, 3), dc(2, 2), dci(2, 2)
    real(kind=8) :: dmc(3, 2), dfc(3, 2)
    real(kind=8) :: pgl(3, 3), propVale(2)
    aster_logical :: coupmf
    integer(kind=8) :: elasID
    character(len=16) :: elasKeyword
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!
    call jevech('PCOURAN', 'E', icour)
!
!     RECUPERATION DES COORDONNEES DES NOEUDS
    call teattr('S', 'DIM_COOR_MODELI', cnd, ier)
    read (cnd, '(I8)') nd
    call jevech('PGEOMER', 'L', jvGeom)
    fami = 'RIGI'

    call elrefe_info(fami=fami, ndim=ndim, nno=nno, nnos=nnos, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde, jgano=jgano)
!
!   CALCUL DE LA PLUS PETITE DISTANCE ENTRE LES NOEUDS SOMMETS
    dmin = sqrt((zr(jvGeom-1+nd*(2-1)+1)-zr(jvGeom-1+1))**2+ &
                (zr(jvGeom-1+nd*(2-1)+2)-zr(jvGeom-1+2))**2+ &
                (zr(jvGeom-1+nd*(2-1)+3)-zr(jvGeom-1+3))**2)
    do i = 1, nnos-1
        do j = i+1, nnos
            xi = zr(jvGeom-1+nd*(i-1)+1)
            yi = zr(jvGeom-1+nd*(i-1)+2)
            xj = zr(jvGeom-1+nd*(j-1)+1)
            yj = zr(jvGeom-1+nd*(j-1)+2)
            if (nd .eq. 3) then
                zii = zr(jvGeom-1+nd*(i-1)+3)
                zj = zr(jvGeom-1+nd*(j-1)+3)
            else
                zii = 0.d0
                zj = 0.d0
            end if
            distij = sqrt((xj-xi)**2+(yj-yi)**2+(zj-zii)**2)
            if ((distij .le. dmin) .and. (distij .ne. 0)) then
                dmin = distij
            end if
        end do
    end do

! - Get material parameters
    call jevech('PMATERC', 'L', jvMaterc)
    call get_elas_id(zi(jvMaterc), elasID, elasKeyword)
!
    if (elasID .eq. ELAS_ISOT) then
        propName(1) = 'E'
        propName(2) = 'NU'
        fami = 'FPG1'
        call rcvalb(fami, 1, 1, '+', zi(jvMaterc), &
                    ' ', elasKeyword, 0, ' ', [0.d0], &
                    2, propName, propVale, propCode, 1)
        e = propVale(1)
        nu = propVale(2)

    elseif (elasID .eq. ELAS_GLRC) then
        propName(1) = 'E_M'
        propName(2) = 'NU_M'
        call rcvalb(fami, 1, 1, '+', zi(jvMaterc), &
                    ' ', elasKeyword, 0, ' ', [0.d0], &
                    2, propName, propVale, propCode, 1)
        e = propVale(1)
        nu = propVale(2)

    elseif (elasID .eq. ELAS_SHELL .or. elasID .eq. ELAS_DHRC) then
        call getCara(plateCara, plateOrie)
        epais = plateCara%thick
        call compCoorSystPara(plateCara, zr(jvGeom), pgl)
        call compCoorSystPlate(pgl, plateCara, plateOrie)
        call dxmate(plateCara, plateOrie, &
                    fami, df, dm, dmf, dc, &
                    dci, dmc, dfc, &
                    multic, coupmf)
        nu = dm(1, 2)/dm(1, 1)
        e = (1.d0-nu**2)*dm(1, 1)/epais

    else if (elasID .eq. ELAS_MEMBRANE) then
        propName(1) = 'M_LLLL'
        call rcvalb(fami, 1, 1, '+', zi(jvMaterc), &
                    ' ', elasKeyword, 0, ' ', [0.d0], &
                    1, propName, propVale, propCode, 1)
        e = propVale(1)

    else
        call utmess('F', 'DYNAMIQUE_32')
    end if

! - Get density
    call rcvalb(fami, 1, 1, '+', zi(jvMaterc), &
                ' ', elasKeyword, 0, ' ', [0.d0], &
                1, 'RHO', propVale, propCode, 1)

!  CALCUL DE LA CELERITE DES ONDES DANS LE MATERIAU
    vitmat = sqrt(e/propVale(1))

! - CALCUL DU PAS DE TEMPS DE LA CONDITION DE COURANT
    zr(icour) = dmin/vitmat
!
end subroutine
