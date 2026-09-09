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
subroutine te0038(option, nomte)
!
    use pipeElem_module
    implicit none
!
#include "asterc/r8prem.h"
#include "asterfort/assert.h"
#include "asterfort/jevech.h"
#include "asterfort/lonele.h"
#include "asterfort/lteatt.h"
#include "asterfort/matrot.h"
#include "asterfort/pmfitx.h"
#include "asterfort/poutre_modloc.h"
#include "asterfort/rccoma.h"
#include "asterfort/rcvalb.h"
#include "asterfort/utmess.h"
#include "asterfort/utpslg.h"
#include "asterfort/utpvlg.h"
#include "jeveux.h"
!
    character(len=*), intent(in):: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: TUYAU_*
!
! Option: MASS_INER
!
! --------------------------------------------------------------------------------------------------
!
    integer(kind=8) :: codres(1), nbNode, nbFourier
    character(len=16) :: elasKeyword
    real(kind=8) :: rho, a1, iy1, iz1, a2, ab2, ab3, ab4, amb, apb, valres(1)
    real(kind=8) :: xl, xl2, massInerLoca(6)
    real(kind=8) :: massInerGlob(6), pgl(3, 3)
    real(kind=8) :: po, poxi2, rr
    real(kind=8) :: ry1, ry2, rz1, rz2, unpr2, unpr4, unprr, xa, xi
    real(kind=8) :: xig, xisl, xixx, xixz, xizz, xzig, yig, zig
    real(kind=8) :: prec
    real(kind=8) :: casect(6), yg, zg, p1gl(3), p1gg(3), rbid
    integer(kind=8), parameter :: nno = 1, nc = 3
    integer(kind=8) :: jvMaterCode, igeom, jvOrien, jvMassIner, beamType
    integer(kind=8), parameter :: nbCara = 9
    real(kind=8) :: caraVale(nbCara)
    character(len=8), parameter :: caraName(nbCara) = (/'A1  ', 'IY1 ', 'IZ1 ', &
                                                        'RY1 ', 'RZ1 ', 'A2  ', &
                                                        'RY2 ', 'RZ2 ', 'TVAR'/)

    aster_logical :: lPipe, lBeamPMF, lBeamPoux, lBeam
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(option .eq. 'MASS_INER')
    prec = r8prem()

    lPipe = lteatt('TUYAU', 'OUI')
    lBeamPoux = lteatt('POUX', 'OUI')
    lBeamPMF = lteatt('TYPMOD2', 'PMF')
    lBeam = lBeamPoux .or. lBeamPMF
    ASSERT(lPipe .or. lBeamPoux .or. lBeamPMF)

! - Get density
    call jevech('PMATERC', 'L', jvMaterCode)
    if (.not. lBeamPMF) then
        call rccoma(zi(jvMaterCode), 'ELAS', 1, elasKeyword, codres(1))
        if ((elasKeyword .eq. 'ELAS') .or. (elasKeyword .eq. 'ELAS_ISTR') .or. &
            (elasKeyword .eq. 'ELAS_FLUI') .or. (elasKeyword .eq. 'ELAS_ORTH')) then
            call rcvalb('FPG1', 1, 1, '+', zi(jvMaterCode), ' ', elasKeyword, 0, ' ', [0.d0], &
                        1, 'RHO', valres, codres, 1)
            rho = valres(1)
        else
            if (lBeamPoux) then
                call utmess('F', 'POUTRE0_50')
            else if (lPipe) then
                call utmess('F', 'PIPE1_50')
            else
                ASSERT(ASTER_FALSE)
            end if
        end if
    end if

! - Get access to output field
    call jevech('PMASSINE', 'E', jvMassIner)

    if (lPipe) then
        call pipeGetDime(nomte, 'RIGI', &
                         nbNode, nbFourier)
        call pipeMassIner(nbNode, nbFourier, rho, jvMassIner)
    elseif (lBeam) then
! ----- Get geometry of beam + length
        xl = lonele(igeom=igeom)
        xl2 = xl*xl

! ----- Get beam orientation
        call jevech('PCAORIE', 'L', jvOrien)
        call matrot(zr(jvOrien), pgl)

! ----- Main properties of beam section
        call poutre_modloc('CAGNPO', caraName, nbCara, lvaleur=caraVale)
        a1 = caraVale(1)
        iy1 = caraVale(2)
        iz1 = caraVale(3)
        ry1 = caraVale(4)
        rz1 = caraVale(5)
        a2 = caraVale(6)
        ry2 = caraVale(7)
        rz2 = caraVale(8)
        beamType = nint(caraVale(9))

!       calcul des caracteristiques elementaires 'MASS_INER'
        massInerGlob = 0.d0
        massInerLoca = 0.d0
        massInerLoca(3) = iy1
        massInerLoca(6) = iz1

        if (beamType .eq. 0) then
!           poutre a section constante
!           masse
            if (lBeamPMF) then
                call pmfitx(zi(jvMaterCode), 2, casect, rbid)
                zr(jvMassIner) = casect(1)*xl
!               correction excentrement
                if (casect(1) .gt. prec) then
                    yg = casect(2)/casect(1)
                    zg = casect(3)/casect(1)
                    p1gl(1) = xl/2.d0
                    p1gl(2) = yg
                    p1gl(3) = zg
                    call utpvlg(1, 3, pgl, p1gl, p1gg)
                    iy1 = casect(5)-casect(1)*zg*zg
                    iz1 = casect(4)-casect(1)*yg*yg
!                   cdg
                    zr(jvMassIner+1) = zr(igeom+1)+p1gg(1)
                    zr(jvMassIner+2) = zr(igeom+2)+p1gg(2)
                    zr(jvMassIner+3) = zr(igeom+3)+p1gg(3)
                else
!                   cdg
                    zr(jvMassIner+1) = (zr(igeom+4)+zr(igeom+1))/2.d0
                    zr(jvMassIner+2) = (zr(igeom+5)+zr(igeom+2))/2.d0
                    zr(jvMassIner+3) = (zr(igeom+6)+zr(igeom+3))/2.d0
                end if
!               inertie
                massInerLoca(1) = (iy1+iz1)*xl
                massInerLoca(2) = 0.d0
                massInerLoca(3) = xl*iy1+casect(1)*xl*xl2/12.d0
                massInerLoca(4) = 0.d0
                massInerLoca(5) = 0.d0
                massInerLoca(6) = xl*iz1+casect(1)*xl*xl2/12.d0
            else
                zr(jvMassIner) = rho*a1*xl
!               cdg
                zr(jvMassIner+1) = (zr(igeom+4)+zr(igeom+1))/2.d0
                zr(jvMassIner+2) = (zr(igeom+5)+zr(igeom+2))/2.d0
                zr(jvMassIner+3) = (zr(igeom+6)+zr(igeom+3))/2.d0
!               inertie
                massInerLoca(1) = rho*(iy1+iz1)*xl
                massInerLoca(2) = 0.d0
                massInerLoca(3) = rho*xl*(iy1+a1*xl2/12.d0)
                massInerLoca(4) = 0.d0
                massInerLoca(5) = 0.d0
                massInerLoca(6) = rho*xl*(iz1+a1*xl2/12.d0)
            end if
            call utpslg(nno, nc, pgl, massInerLoca, massInerGlob)
!
        else if (beamType .eq. 1) then
!           poutre a section variable affine
            if ((abs(a1-(4.d0*ry1*rz1)) .gt. (a1*prec)) .or. &
                (abs(a2-(4.d0*ry2*rz2)) .gt. (a2*prec))) then
                call utmess('F', 'POUTRE0_81')
            end if
!           masse
            zr(jvMassIner) = rho*xl*(a1+a2)/2.d0
!           cdg
            xisl = (rz1+2.d0*rz2)/(3.d0*(rz1+rz2))
            zr(jvMassIner+1) = zr(igeom+1)+(zr(igeom+4)-zr(igeom+1))*xisl
            zr(jvMassIner+2) = zr(igeom+2)+(zr(igeom+5)-zr(igeom+2))*xisl
            zr(jvMassIner+3) = zr(igeom+3)+(zr(igeom+6)-zr(igeom+3))*xisl
!           inertie
            xa = xl*(rz1+rz2)
            amb = rz1-rz2
            apb = rz1+rz2
            ab2 = rz1**2+rz2**2+4.d0*rz1*rz2
            ab3 = rz1**3+3.d0*rz1**2*rz2-3.d0*rz1*rz2**2-rz2**3
            ab4 = rz1**4+rz2**4+2.d0*rz1*rz2*(rz1**2+rz2**2)
!
            xixx = xl*(4.d0*ab4-2.d0*amb*ab3+amb**2*ab2)/(18.d0*apb)
            xizz = (xl**3)*ab2/(18.d0*apb)
            xixz = (xl**2)*(ab3-amb*ab2)/(18.d0*apb)
            xig = rho*((xa*2.d0*ry1**3/3.d0)+2.d0*ry1*xixx)
            yig = rho*(2.d0*ry1*(xixx+xizz))
            zig = rho*((xa*2.d0*ry1**3/3.d0)+2.d0*ry1*xizz)
            xzig = rho*(2.d0*ry1*xixz)
            massInerLoca(1) = xig
            massInerLoca(2) = 0.d0
            massInerLoca(3) = yig
            massInerLoca(4) = xzig
            massInerLoca(5) = 0.d0
            massInerLoca(6) = zig
            call utpslg(nno, nc, pgl, massInerLoca, massInerGlob)
!
        else if (beamType .eq. 2) then
!           poutre a section variable homothetique
            if (a1 .eq. 0.d0) then
                call utmess('F', 'POUTRE0_82')
            end if
!           masse
            zr(jvMassIner) = rho*(a1+a2+sqrt(a1*a2))*xl/3.d0
!           CDG
            rr = sqrt(a2/a1)
            unprr = 1.d0+rr+rr**2
            xi = (1.d0+2.d0*rr+3.d0*(rr**2))/(4.d0*unprr)
            zr(jvMassIner+1) = zr(igeom+1)*(1.d0-xi)+zr(igeom+4)*xi
            zr(jvMassIner+2) = zr(igeom+2)*(1.d0-xi)+zr(igeom+5)*xi
            zr(jvMassIner+3) = zr(igeom+3)*(1.d0-xi)+zr(igeom+6)*xi
!           inertie
            unpr4 = unprr+rr**3+rr**4
            unpr2 = 1.d0+3.d0*rr+6.d0*rr**2
            po = rho*xl*a1*unprr/3.d0
            xig = rho*xl*(iy1+iz1)*unpr4/5.d0
            poxi2 = rho*(xl**3)*a1*unpr2/30.d0-po*((xi*xl)**2)
            yig = rho*xl*iy1*unpr4/5.d0+poxi2
            zig = rho*xl*iz1*unpr4/5.d0+poxi2
            massInerLoca(1) = xig
            massInerLoca(2) = 0.d0
            massInerLoca(3) = yig
            massInerLoca(4) = 0.d0
            massInerLoca(5) = 0.d0
            massInerLoca(6) = zig
            call utpslg(nno, nc, pgl, massInerLoca, massInerGlob)
        end if

        zr(jvMassIner+3+1) = massInerGlob(1)
        zr(jvMassIner+3+2) = massInerGlob(3)
        zr(jvMassIner+3+3) = massInerGlob(6)
        zr(jvMassIner+3+4) = massInerGlob(2)
        zr(jvMassIner+3+5) = massInerGlob(4)
        zr(jvMassIner+3+6) = massInerGlob(5)
!
    else
        ASSERT(ASTER_FALSE)
    end if
!
end subroutine
