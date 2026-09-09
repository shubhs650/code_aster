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
subroutine te0444(option, nomte)
!
    use plate_type
    use plateGeom_module, only: getCara, compCoorSystPara, compCoorSystPlate
    implicit none
!
#include "asterfort/assert.h"
#include "asterfort/dkqmas.h"
#include "asterfort/dkqrig.h"
#include "asterfort/dktmas.h"
#include "asterfort/dktrig.h"
#include "asterfort/dxiner.h"
#include "asterfort/dxroep.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/pmavec.h"
#include "asterfort/q4grig.h"
#include "asterfort/t3grig.h"
#include "asterfort/utpslg.h"
#include "asterfort/utpvgl.h"
#include "asterfort/vecma.h"
#include "jeveux.h"
!
    character(len=16) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: DKTG/Q4GG
!
! Options: ECIN_ELEM
!          EPOT_ELEM
!          MASS_INER
!          MASS_MECA*
!          M_GAMMA
!
! --------------------------------------------------------------------------------------------------
!
! In  option           : name of option to compute
! In  nomte            : type of finite element
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: i, j
    integer(kind=8) :: jvGeom, jvMatr, jvEner, jvOmega, jvAcce, jvVect, jvMassIner
    integer(kind=8) :: nno, nddl, n1, ni, n2, nbTermSyme
    real(kind=8) :: rho, epais
    real(kind=8) :: pgl(3, 3), xyzl(3, 4)
    real(kind=8) :: ener(3), matrFull(24, 24), matrMassGlob(300)
!     ---> POUR DKT MATELEM = 3 * 6 DDL = 171 TERMES STOCKAGE SYME
!     ---> POUR DKQ MATELEM = 4 * 6 DDL = 300 TERMES STOCKAGE SYME
    real(kind=8) :: matrMass(300)
    type(plateCara_Para) :: plateCara
    type(plateOrie_Para) :: plateOrie
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami='RIGI', nno=nno)

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Get plate parameters
    call getCara(plateCara, plateOrie)

! - Calculate the transformation: global coordinate system/intrinsic coordinate system
    call compCoorSystPara(plateCara, zr(jvGeom), pgl)

! - Compute coordinate system for plate
    call compCoorSystPlate(pgl, plateCara, plateOrie)

! - Change coordinates of displacements
    call utpvgl(nno, 3, pgl, zr(jvGeom), xyzl)

    if (option .eq. 'EPOT_ELEM') then
        if (nomte .eq. 'MEDKTG3') then
            call dktrig(plateCara, plateOrie, &
                        xyzl, option, pgl, &
                        ener_=ener)
        else if (nomte .eq. 'MEDKQG4') then
            call dkqrig(plateCara, plateOrie, &
                        xyzl, option, pgl, &
                        ener_=ener)
        else if (nomte .eq. 'MET3GG3') then
            call t3grig(plateCara, plateOrie, &
                        xyzl, option, pgl, &
                        ener_=ener)
        else if (nomte .eq. 'MEQ4GG4') then
            call q4grig(plateCara, plateOrie, &
                        xyzl, option, pgl, &
                        ener_=ener)
        end if
        call jevech('PENERDR', 'E', jvEner)
        do i = 1, 3
            zr(jvEner-1+i) = ener(i)
        end do

    else if (option .eq. 'MASS_MECA' .or. option .eq. 'MASS_MECA_DIAG' &
             .or. option .eq. 'MASS_MECA_EXPLI' .or. option .eq. 'M_GAMMA' &
             .or. option .eq. 'ECIN_ELEM') then
        if (nomte .eq. 'MEDKTG3' .or. nomte .eq. 'MET3GG3') then
            call dktmas(plateCara, &
                        xyzl, option, pgl, &
                        matrMass, ener)
        else if (nomte .eq. 'MEDKQG4' .or. nomte .eq. 'MEQ4GG4') then
            call dkqmas(plateCara, plateOrie, &
                        xyzl, option, pgl, &
                        matrMass, ener)
        end if

        if (option .eq. 'MASS_MECA') then
            call jevech('PMATUUR', 'E', jvMatr)
            call utpslg(nno, 6, pgl, matrMass, zr(jvMatr))
        else if (option .eq. 'ECIN_ELEM') then
            call jevech('PENERCR', 'E', jvEner)
            call jevech('POMEGA2', 'L', jvOmega)
            do i = 1, 3
                zr(jvEner-1+i) = zr(jvOmega)*ener(i)
            end do
        else if (option .eq. 'M_GAMMA') then
            call jevech('PACCELR', 'L', jvAcce)
            call jevech('PVECTUR', 'E', jvVect)
            nddl = 6*nno
            nbTermSyme = nddl*(nddl+1)/2
            call utpslg(nno, 6, pgl, matrMass, matrMassGlob)
            call vecma(matrMassGlob, nbTermSyme, matrFull, nddl)
            call pmavec('ZERO', nddl, matrFull, zr(jvAcce), zr(jvVect))
        else if (option .eq. 'MASS_MECA_DIAG' .or. option .eq. 'MASS_MECA_EXPLI') then
            call jevech('PMATUUR', 'E', jvMatr)
            nddl = 6*nno
            nbTermSyme = nddl*(nddl+1)/2
            do i = 1, nbTermSyme
                zr(jvMatr-1+i) = matrMass(i)
            end do
            if (option .eq. 'MASS_MECA_EXPLI') then
!     CORRECTION DES TERMES CORRESPONDANT AU DDL 6
!     NON PREVU PAR LA THEORIE DKT. ON RAJOUTE
!     UN TERME DIAGONAL NON ZERO EGAL A CELUI DU DDL 5.
!     CETTE CORRECTION A ETE INSPIRE PAR LA DEMARCHE DANS EUROPLEXUS
                do j = 1, nno
                    n1 = 6*(j-1)+5
                    n2 = 6*(j-1)+4
                    ni = 6*j
                    nbTermSyme = (ni+1)*ni/2
                    n1 = (n1+1)*n1/2
                    n2 = (n2+1)*n2/2
                    zr(jvMatr-1+nbTermSyme) = (zr(jvMatr-1+n1)+zr(jvMatr-1+n2))*0.5d0
                end do
            end if
        end if
    else if (option .eq. 'MASS_INER') then
        call jevech('PMASSINE', 'E', jvMassIner)
        call dxroep(plateCara, rho, epais)
        call dxiner(plateCara, &
                    zr(jvGeom), rho, epais, &
                    zr(jvMassIner), zr(jvMassIner+1), zr(jvMassIner+4))
    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
