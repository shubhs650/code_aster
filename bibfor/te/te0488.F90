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
subroutine te0488(option, nomte)
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystPara, compCoorSystNone, compCoorSystCO3D
    implicit none
!
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/dfdm3d.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/lteatt.h"
#include "asterfort/plate_type.h"
#include "asterfort/subaco.h"
#include "asterfort/sumetr.h"
#include "asterfort/tecach.h"
#include "asterfort/utpvlg.h"
#include "asterfort/vectgt.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: all
!
! Options: COOR_ELGA
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: jv_poids, jv_vf, jvGeom, jv_coopg, jv_dfde
    integer(kind=8) :: nno, kpg, npg, ino, ndim
    real(kind=8) :: xx, yy, zz, poids, cova(3, 3), metr(2, 2), jac
    integer(kind=8) ::  nbLayer, decpo, iLayer, ispc, jtab(7), nbsp, iret
    real(kind=8) :: epais, excen, bas, epc, pgl(3, 3), gm2(3)
    integer(kind=8) :: lzi, lzr, nb1
    real(kind=8) :: vectBaseKpg(3, 3)
    real(kind=8) :: hh
    real(kind=8), parameter :: zero = 0.d0
    aster_logical :: l_coq3d, l_grille, l_solid_shell, lPlate
    real(kind=8), parameter :: gm1(3) = (/0.d0, 0.d0, 1.d0/)
    real(kind=8), parameter :: poidc(3) = (/0.16666666666666666d0, 0.66666666666666663d0, &
                                            0.16666666666666666d0/)
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!
    l_coq3d = lteatt('MODELI', 'CQ3')
    l_grille = lteatt('MODELI', 'GRC')
    l_solid_shell = lteatt('MODELI', 'SSH')
    lPlate = lteatt('COQUE', 'OUI') .or. lteatt('PLAQUE', 'OUI')

    if (l_coq3d) then
        call elrefe_info(fami='MASS', ndim=ndim, nno=nno, npg=npg, &
                         jpoids=jv_poids, jvf=jv_vf, jdfde=jv_dfde)
    else
        call elrefe_info(fami='RIGI', ndim=ndim, nno=nno, npg=npg, &
                         jpoids=jv_poids, jvf=jv_vf, jdfde=jv_dfde)
    end if

! - Access to input fields
    call jevech('PGEOMER', 'L', jvGeom)

! - Access to output fields
    call tecach('OOO', 'PCOORPG', 'E', iret, nval=7, itab=jtab)
    jv_coopg = jtab(1)
    nbsp = jtab(7)

! - Get parameters for structural elements
    if (lPlate) then
        call getCara(plateCara, plateOrie)
        call compCoorSystPara(plateCara, zr(jvGeom), pgl)
        call compCoorSystNone(plateOrie)
    end if

    if (l_grille) then
        nbsp = 1
        call utpvlg(1, 3, pgl, gm1, gm2)
        excen = plateCara%offset
    end if
    if (nbsp .ne. 1) then
        ASSERT(plateCara%type .ne. PLATE_UNKW)
        nbLayer = plateCara%nbLayer
        epais = plateCara%thick
        excen = plateCara%offset
        if (l_coq3d) then
            call jevete('&INEL.'//nomte(1:8)//'.DESI', ' ', lzi)
            nb1 = zi(lzi-1+1)
            npg = zi(lzi-1+4)
            call jevete('&INEL.'//nomte(1:8)//'.DESR', ' ', lzr)
            call compCoorSystCO3D(nomte, jvGeom, &
                                  plateCara, plateOrie)
        else
            call compCoorSystPara(plateCara, zr(jvGeom), pgl)
            call compCoorSystNone(plateOrie)
            call utpvlg(1, 3, pgl, gm1, gm2)
        end if
        bas = -epais/2.d0+excen
        epc = epais/nbLayer
    end if
!
    do kpg = 1, npg
        xx = zero
        yy = zero
        zz = zero
        if (l_solid_shell) then
            do ino = 1, nno-1
                xx = xx+zr(jvGeom+3*(ino-1)+0)*zr(jv_vf+(kpg-1)*nno+ino-1)
                yy = yy+zr(jvGeom+3*(ino-1)+1)*zr(jv_vf+(kpg-1)*nno+ino-1)
                zz = zz+zr(jvGeom+3*(ino-1)+2)*zr(jv_vf+(kpg-1)*nno+ino-1)
            end do
        else
            do ino = 1, nno
                xx = xx+zr(jvGeom+3*(ino-1)+0)*zr(jv_vf+(kpg-1)*nno+ino-1)
                yy = yy+zr(jvGeom+3*(ino-1)+1)*zr(jv_vf+(kpg-1)*nno+ino-1)
                zz = zz+zr(jvGeom+3*(ino-1)+2)*zr(jv_vf+(kpg-1)*nno+ino-1)
            end do
        end if
        if (ndim .eq. 3) then
            call dfdm3d(nno, kpg, jv_poids, jv_dfde, zr(jvGeom), poids)
        else if (ndim .eq. 2) then
            call subaco(nno, zr(jv_dfde+(kpg-1)*ndim*nno), zr(jvGeom), cova)
            call sumetr(cova, metr, jac)
            poids = jac*zr(jv_poids-1+kpg)
        else
            ASSERT(ASTER_FALSE)
        end if
!
        if (nbsp .ne. 1) then
            decpo = 4*3*nbLayer*(kpg-1)
            if (l_coq3d) then
                call vectgt(plateOrie, 1, nb1, &
                            zr(jvGeom), 0.d0, kpg, &
                            epais, zr(lzr), &
                            vectBaseKpg)
                gm2(1) = vectBaseKpg(3, 1)
                gm2(2) = vectBaseKpg(3, 2)
                gm2(3) = vectBaseKpg(3, 3)
            end if
            do iLayer = 1, nbLayer
                do ispc = 1, 3
                    hh = bas+dble(iLayer-1)*epc+dble(ispc-1)*epc/2.d0
                    zr(jv_coopg+decpo+(iLayer-1)*12+(ispc-1)*4+0) = xx+hh*gm2(1)
                    zr(jv_coopg+decpo+(iLayer-1)*12+(ispc-1)*4+1) = yy+hh*gm2(2)
                    zr(jv_coopg+decpo+(iLayer-1)*12+(ispc-1)*4+2) = zz+hh*gm2(3)
                    zr(jv_coopg+decpo+(iLayer-1)*12+(ispc-1)*4+3) = poids*epc*poidc(ispc)
                end do
            end do
        else if (l_grille) then
            zr(jv_coopg+4*(kpg-1)+0) = xx+excen*gm2(1)
            zr(jv_coopg+4*(kpg-1)+1) = yy+excen*gm2(2)
            zr(jv_coopg+4*(kpg-1)+2) = zz+excen*gm2(3)
            zr(jv_coopg+4*(kpg-1)+3) = poids
        else
            zr(jv_coopg+4*(kpg-1)+0) = xx
            zr(jv_coopg+4*(kpg-1)+1) = yy
            zr(jv_coopg+4*(kpg-1)+2) = zz
            zr(jv_coopg+4*(kpg-1)+3) = poids
        end if
    end do
!
end subroutine
