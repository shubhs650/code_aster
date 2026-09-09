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
subroutine te0403(option, nomte)
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystCO3D
    implicit none
!
#include "asterfort/fcent.h"
#include "asterfort/fointe.h"
#include "asterfort/fpesa.h"
#include "asterfort/fpres.h"
#include "asterfort/fsurf.h"
#include "asterfort/jevech.h"
#include "asterfort/jevete.h"
#include "asterfort/tecael.h"
#include "asterfort/trnflg.h"
#include "asterfort/utmess.h"
#include "asterfort/vectan.h"
#include "jeveux.h"
!
    character(len=16) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: COQUE_3D
!
! Options: CHAR_MECA_FRCO3D
!          CHAR_MECA_FFCO3D
!          CHAR_MECA_PRES_R
!          CHAR_MECA_PRES_F
!          CHAR_MECA_PESA_R
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbPara = 4
    character(len=8), parameter :: paraName(nbPara) = (/"X   ", "Y   ", "Z   ", "INST"/)
    real(kind=8) :: paraVale(nbPara)
    integer(kind=8) :: nb1
    real(kind=8) :: vectNorm(9, 3), vectTang(9, 2, 3), vectBase(9, 3, 3)
    real(kind=8) :: vecl(51)
    real(kind=8) :: pr
    integer(kind=8) :: iadzi, iazk24, iNode, ier, itemps, j
    integer(kind=8) :: jvGeom, jpres, jvecg, lzi, lzr, nb2
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!

! - Access objects
    call jevete('&INEL.'//nomte(1:8)//'.DESR', ' ', lzr)
    call jevete('&INEL.'//nomte(1:8)//'.DESI', ' ', lzi)
    nb1 = zi(lzi-1+1)
    nb2 = zi(lzi-1+2)

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

    if ((option .eq. 'CHAR_MECA_PESA_R') .or. (option .eq. 'CHAR_MECA_ROTA_R')) then
! ----- Get plate parameters
        call getCara(plateCara, plateOrie)

! ----- Compute global<=>local transformation
        call compCoorSystCO3D(nomte, jvGeom, &
                              plateCara, plateOrie)

        vectNorm = plateOrie%vectNorm
        vectTang = plateOrie%vectTang
    else
! ---- Compute local basis at nodes
        call vectan(nb1, nb2, &
                    zr(jvGeom), zr(lzr), &
                    vectNorm, vectTang)

    end if

! - Fuse tangents and normal in same object
    do iNode = 1, nb2
        vectBase(iNode, 1:2, 1:3) = vectTang(iNode, 1:2, 1:3)
        vectBase(iNode, 3, 1:3) = vectNorm(iNode, 1:3)
    end do

    if (option .eq. 'CHAR_MECA_FRCO3D' .or. option .eq. 'CHAR_MECA_FFCO3D') then
!------------------------------------------------------
!      PAS DE CHANGEMENT DE SIGNE POUR LES FORCES REPARTIES
!------------------------------------------------------
        call fsurf(option, nomte, zr(jvGeom), nb1, vecl, &
                   vectBase)

    else if (option .eq. 'CHAR_MECA_PESA_R') then
        call fpesa(plateCara, &
                   nomte, zr(jvGeom), nb1, &
                   vecl)

    else if (option .eq. 'CHAR_MECA_ROTA_R') then
        call fcent(plateCara, &
                   nomte, zr(jvGeom), nb1, &
                   vecl)

    else if (option .eq. 'CHAR_MECA_PRES_R') then
!------------------------------------------------------
!      CHANGEMENT DE SIGNE POUR LES PRESSIONS DANS FPRES
!------------------------------------------------------
        call fpres(nomte, zr(jvGeom), nb1, vecl, vectBase)

    else if (option .eq. 'CHAR_MECA_PRES_F') then
! ----- Only zero is allowed
        call jevech('PPRESSF', 'L', jpres)
        if (zk8(jpres) .eq. '&FOZERO') goto 999
        call jevech('PINSTR', 'L', itemps)
        paraVale(4) = zr(itemps)
        do j = 0, nb1-1
            paraVale(1) = zr(jvGeom+3*j)
            paraVale(2) = zr(jvGeom+3*j+1)
            paraVale(3) = zr(jvGeom+3*j+2)
            call fointe('FM', zk8(jpres), nbPara, paraName, paraVale, pr, ier)
            if (pr .ne. 0.d0) then
                call tecael(iadzi, iazk24)
                call utmess('F', 'ELEMENTS4_92', si=zi(iadzi-1+1))
            end if
        end do
        goto 999
    end if
!
    call jevech('PVECTUR', 'E', jvecg)
    call trnflg(nb2, vectBase, vecl, zr(jvecg))
!
999 continue
end subroutine
