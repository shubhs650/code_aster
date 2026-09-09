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
subroutine te0430(option, nomte)
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystPara, compCoorSystGrid
    use resi_refe_module, only: RESI_REFE
    implicit none
!
#include "asterc/r8vide.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/cargri.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/fointe.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/nmgrib.h"
#include "asterfort/rcvalb.h"
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"
#include "asterfort/verift.h"
#include "jeveux.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: GRILLE_MEMBRANE / GRILLE_EXCENTRE
!
! Elements: CHAR_MECA_EPSI_R
!           CHAR_MECA_PESA_R
!           CHAR_MECA_TEMP_R
!           FORC_NODA
!           REFE_FORC_NODA
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8), parameter :: nbProp = 1
    integer(kind=8) :: propCode(nbProp)
    character(len=16), parameter :: propName(nbProp) = (/'E'/)
    real(kind=8) :: propVale(nbProp)
    integer(kind=8), parameter :: nbPara = 4
    character(len=8), parameter :: paraName(nbPara) = (/'X   ', 'Y   ', 'Z   ', 'INST'/)
    real(kind=8) :: paraVale(nbPara)
    integer(kind=8) :: nddl, nno, npg, i, kpg, n, iret, ier
    integer(kind=8) :: ipoids, ivf, idfde, jvGeom, jvMaterc, jvSigm, jvVect
    integer(kind=8) :: jvPesa, jvEpsi, jvInstr
    integer(kind=8) :: iadzi, iazk24
    real(kind=8) :: dff(2, 8), vff(8), b(6, 8), p(3, 6), jac, epsthe
    real(kind=8) :: dir11(3), densit, pgl(3, 3), distn
    character(len=8), parameter :: fami = 'RIGI'
    real(kind=8) :: sig, rho, b_max_rot
    aster_logical :: lexc
    real(kind=8) :: xgau, ygau, zgau, exx
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
    type(RESI_REFE):: refe
!
! --------------------------------------------------------------------------------------------------
!
    lexc = (lteatt('MODELI', 'GRC'))

! - Get plate parameters
    call getCara(plateCara, plateOrie)

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Calculate the transformation: global coordinate system/intrinsic coordinate system
    if (lexc) then
        call compCoorSystPara(plateCara, zr(jvGeom), pgl)
        nddl = 6
    else
        nddl = 3
    end if
    call compCoorSystGrid(pgl, plateCara, plateOrie)

! - FONCTIONS DE FORMES ET POINTS DE GAUSS
    call elrefe_info(fami='RIGI', nno=nno, npg=npg, &
                     jpoids=ipoids, jvf=ivf, jdfde=idfde)

! - Input fields
    if (option .eq. 'FORC_NODA') then
        call jevech('PSIEFR', 'L', jvSigm)
    else if (option .eq. 'REFE_FORC_NODA') then
        call jevech('PMATERC', 'L', jvMaterc)
    else if (option .eq. 'CHAR_MECA_EPSI_R') then
        call jevech('PMATERC', 'L', jvMaterc)
        call jevech('PEPSINR', 'L', jvEpsi)
    else if (option .eq. 'CHAR_MECA_EPSI_F') then
        call jevech('PMATERC', 'L', jvMaterc)
        call jevech('PEPSINF', 'L', jvEpsi)
        call jevech('PINSTR', 'L', jvInstr)
    else if (option .eq. 'CHAR_MECA_PESA_R') then
        call jevech('PMATERC', 'L', jvMaterc)
        call jevech('PPESANR', 'L', jvPesa)
    else if (option .eq. 'CHAR_MECA_TEMP_R') then
        call jevech('PMATERC', 'L', jvMaterc)
    end if

! - Output fields
    call jevech('PVECTUR', 'E', jvVect)

! - LECTURE DES CARACTERISTIQUES DE GRILLE ET CALCUL DE LA DIRECTION D'ARMATURE
    call cargri(plateCara, plateOrie, &
                densit, distn, dir11)

    do kpg = 1, npg
        do n = 1, nno
            vff(n) = zr(ivf+(kpg-1)*nno+n-1)
            dff(1, n) = zr(idfde+(kpg-1)*nno*2+(n-1)*2)
            dff(2, n) = zr(idfde+(kpg-1)*nno*2+(n-1)*2+1)
        end do

! ----- CALCUL DE LA MATRICE "B" : DEPL NODAL --> EPS11 ET DU JACOBIEN
        call nmgrib(nno, zr(jvGeom), dff, dir11, lexc, &
                    plateOrie%gridNorm, b, jac, p)

        if ((option .eq. 'FORC_NODA') .or. (option .eq. 'CHAR_MECA_TEMP_R') .or. &
            (option(1:15) .eq. 'CHAR_MECA_EPSI_')) then

            if (option .eq. 'FORC_NODA') then
                sig = zr(jvSigm+kpg-1)

            else if (option .eq. 'CHAR_MECA_EPSI_R') then
                call rcvalb(fami, kpg, 1, '+', &
                            zi(jvMaterc), ' ', 'ELAS', &
                            0, ' ', [0.d0], &
                            nbProp, propName, propVale, &
                            propCode, 1)
                sig = propVale(1)*zr(jvEpsi+kpg-1)

            else if (option .eq. 'CHAR_MECA_EPSI_F') then
                call rcvalb(fami, kpg, 1, '+', &
                            zi(jvMaterc), ' ', 'ELAS', &
                            0, ' ', [0.d0], &
                            nbProp, propName, propVale, &
                            propCode, 1)
                xgau = 0.d0
                ygau = 0.d0
                zgau = 0.d0
                do i = 1, nno
                    xgau = xgau+zr(ivf-1+i+nno*(kpg-1))*zr(jvGeom-1+1+3*(i-1))
                    ygau = ygau+zr(ivf-1+i+nno*(kpg-1))*zr(jvGeom-1+2+3*(i-1))
                    zgau = zgau+zr(ivf-1+i+nno*(kpg-1))*zr(jvGeom-1+3+3*(i-1))
                end do
                paraVale(1) = xgau
                paraVale(2) = ygau
                paraVale(3) = zgau
                paraVale(4) = zr(jvInstr)

                call fointe('FM', zk8(jvEpsi), nbPara, paraName, paraVale, exx, ier)
                sig = propVale(1)*exx

            else if (option .eq. 'CHAR_MECA_TEMP_R') then
                call verift(fami, kpg, 1, '+', zi(jvMaterc), &
                            iret_=iret, epsth_=epsthe)
                if (iret .ne. 0) then
                    call tecael(iadzi, iazk24)
                    call utmess('S', 'CALCULEL2_81', si=zi(iadzi-1+1))
                end if
                call rcvalb(fami, kpg, 1, '+', &
                            zi(jvMaterc), ' ', 'ELAS', &
                            0, ' ', [0.d0], &
                            nbProp, propName, &
                            propVale, propCode, 1)
                sig = propVale(1)*epsthe
            end if
            do n = 1, nno
                do i = 1, nddl
                    zr(jvVect+(n-1)*nddl+i-1) = zr(jvVect+(n-1)*nddl+i-1)+ &
                                                b(i, n)*sig*zr(ipoids+kpg-1)*jac*densit
                end do
            end do
!
! - REFE_FORC_NODA : ON CALCULE DES FORCES DE REFERENCE
!      (N'EST VALABLE QUE POUR LES GRILLES MEMBRANES)
!
        else if (option .eq. 'REFE_FORC_NODA') then
            call refe%Init(nomte)
            sig = refe%GetRef('SIGM')
            call refe%Check()
!
            do n = 1, nno
                do i = 1, 3
                    zr(jvVect+(n-1)*nddl+i-1) = zr(jvVect+(n-1)*nddl+i-1)+ &
                                                sig*sqrt(abs(jac))*densit/npg
                end do
                b_max_rot = 0.d0
                do i = 4, nddl
                    if (abs(b(i, n)) .gt. b_max_rot) then
                        b_max_rot = abs(b(i, n))
                    end if
                end do
                do i = 4, nddl
                    zr(jvVect+(n-1)*nddl+i-1) = zr(jvVect+(n-1)*nddl+i-1)+ &
                                                b_max_rot*sig*sqrt(abs(jac))*densit/npg
                end do
            end do

        else if (option .eq. 'CHAR_MECA_PESA_R') then
            call rcvalb(fami, kpg, 1, '+', &
                        zi(jvMaterc), ' ', 'ELAS', &
                        0, ' ', [0.d0], &
                        1, 'RHO', &
                        propVale, propCode, 1)
            rho = propVale(1)
            do n = 1, nno
                do i = 1, 3
                    zr(jvVect+(n-1)*nddl+i-1) = zr(jvVect+(n-1)*nddl+i-1)+ &
                                                rho*zr(ipoids+kpg-1)*zr(jvPesa)*zr(jvPesa+i)* &
                                                vff(n)*densit*jac
                end do
            end do
        end if
    end do
!
end subroutine
