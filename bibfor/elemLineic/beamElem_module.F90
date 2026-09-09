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
! ==================================================================================================
!
! Module for beam elements (EF)
!
! ==================================================================================================
!
module beamElem_module
! ==================================================================================================
    use beamElem_type
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: beamBMatrStraight, beamBMatrCurved
    public :: beamNMatr
    private :: beamNMatrStraight, beamNMatrCurved
! ==================================================================================================
    private
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/beamElem_type.h"
#include "jeveux.h"
! ==================================================================================================
contains
! --------------------------------------------------------------------------------------------------
!
! beamBMatrStraight
!
! Compute B matrix (u => epsi) for straight beam
!
! In  beamElemProp     : properties of beam
! In  nbDof            : number of DOF in element
! In  kpg              : current Gauss point
! In  phi              : coordinate of sub-point in section
! In  radiusSection    : radius of section (only circle)
! In  shiftDOF         : shift for B matrix
! In  ff               : shape functions on beam
! In  df1              : first derivative of shape functions on beam
! IO  b                : B matrix (u => epsi)
!
! --------------------------------------------------------------------------------------------------
    subroutine beamBMatrStraight(beamElemProp, nbDof, kpg, &
                                 phi, radiusSection, shiftDOF, &
                                 ff, df1, b)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(beamElem_Prop), intent(in) :: beamElemProp
        integer(kind=8), intent(in) :: nbDof, kpg
        real(kind=8), intent(in) :: phi, radiusSection
        integer(kind=8), intent(in) :: shiftDOF
        real(kind=8), intent(in) :: ff(*), df1(*)
        real(kind=8), intent(inout) :: b(4, nbDof)
! ----- Local
        integer(kind=8) :: iNode, iBloc, nbNode
        real(kind=8) :: hk, dhk, cosfi, sinfi
!   ------------------------------------------------------------------------------------------------
!
        ASSERT(beamElemProp%sectType .eq. BEAM_SECT_PIPE)
        cosfi = cos(phi)
        sinfi = sin(phi)
        nbNode = beamElemProp%nbNode

        do iNode = 1, nbNode
! --------- Get shape functions
            hk = ff(nbNode*(kpg-1)+iNode)
            dhk = df1(nbNode*(kpg-1)+iNode)*(2.d0/beamElemProp%elemLength)

! --------- Select column
            iBloc = (shiftDOF)*(iNode-1)
            ASSERT(iBloc+6 .le. nbDof)

! --------- Set values
            b(1, iBloc+1) = dhk
            b(1, iBloc+2) = 0.d0
            b(1, iBloc+3) = 0.d0
            b(1, iBloc+4) = 0.d0
            b(1, iBloc+5) = -radiusSection*cosfi*dhk
            b(1, iBloc+6) = radiusSection*sinfi*dhk

            b(2, iBloc+1) = 0.d0
            b(2, iBloc+2) = 0.d0
            b(2, iBloc+3) = 0.d0
            b(2, iBloc+4) = 0.d0
            b(2, iBloc+5) = 0.d0
            b(2, iBloc+6) = 0.d0

            b(3, iBloc+1) = 0.d0
            b(3, iBloc+2) = -cosfi*dhk
            b(3, iBloc+3) = sinfi*dhk
            b(3, iBloc+4) = -radiusSection*dhk
            b(3, iBloc+5) = sinfi*hk
            b(3, iBloc+6) = cosfi*hk

            b(4, iBloc+1) = 0.d0
            b(4, iBloc+2) = -sinfi*dhk
            b(4, iBloc+3) = -cosfi*dhk
            b(4, iBloc+4) = 0.d0
            b(4, iBloc+5) = -cosfi*hk
            b(4, iBloc+6) = sinfi*hk

        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! beamBMatrCurved
!
! Compute B matrix (u => epsi) for curved beam
!
! In  beamElemProp     : properties of beam
! In  nbDof            : number of DOF in element
! In  npg              : number of Gauss point
! In  kpg              : index of Gauss point
! In  xpg              : coordinates of Gauss points
! In  phi              : coordinate of sub-point in section
! In  radiusSection    : radius of section (only circle)
! In  shiftDOF         : shift for B matrix
! In  ff               : shape functions on beam
! In  df1              : first derivative of shape functions on beam
! IO  b                : B matrix (u => epsi)
!
! --------------------------------------------------------------------------------------------------
    subroutine beamBMatrCurved(beamElemProp, nbDof, &
                               npg, kpg, xpg, &
                               phi, radiusSection, shiftDOF, &
                               ff, df1, b)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(beamElem_Prop), intent(in) :: beamElemProp
        integer(kind=8), intent(in) :: nbDof, kpg, npg
        real(kind=8), intent(in) :: xpg(npg)
        real(kind=8), intent(in) :: phi, radiusSection
        integer(kind=8), intent(in) :: shiftDOF
        real(kind=8), intent(in) :: ff(*), df1(*)
        real(kind=8), intent(inout) :: b(4, nbDof)
! ----- Local
        integer(kind=8) :: iNode, iBloc, nbNode
        real(kind=8) :: hk, dhk, cosfi, sinfi, ck, sk
        real(kind=8) :: denr, tk(4), radiusElbow, thetaElbow
!   ------------------------------------------------------------------------------------------------
!
        nbNode = beamElemProp%nbNode
        tk = beamElemProp%tk
        radiusElbow = beamElemProp%radiusElbow
        thetaElbow = beamElemProp%thetaElbow
        cosfi = cos(phi)
        sinfi = sin(phi)
        denr = radiusElbow+radiusSection*sinfi

        do iNode = 1, nbNode
            hk = ff(nbNode*(kpg-1)+iNode)
            dhk = df1(nbNode*(kpg-1)+iNode)*(2.d0/thetaElbow)
            ck = cos((1.d0+xpg(kpg))*thetaElbow/2.d0-tk(iNode))
            sk = sin((1.d0+xpg(kpg))*thetaElbow/2.d0-tk(iNode))

! --------- Select column
            iBloc = (shiftDOF)*(iNode-1)
            ASSERT(iBloc+6 .le. nbDof)

! --------- Set values
            b(1, iBloc+1) = dhk*ck/denr
            b(1, iBloc+2) = dhk*sk/denr
            b(1, iBloc+3) = 0.d0
            b(1, iBloc+4) = radiusSection*cosfi*dhk*sk/denr
            b(1, iBloc+5) = -radiusSection*cosfi*dhk*ck/denr
            b(1, iBloc+6) = radiusSection*sinfi*dhk/denr

            b(2, iBloc+1) = 0.d0
            b(2, iBloc+2) = 0.d0
            b(2, iBloc+3) = 0.d0
            b(2, iBloc+4) = 0.d0
            b(2, iBloc+5) = 0.d0
            b(2, iBloc+6) = 0.d0

            b(3, iBloc+1) = cosfi*dhk*sk/denr
            b(3, iBloc+2) = -cosfi*dhk*ck/denr
            b(3, iBloc+3) = sinfi*dhk/denr
            b(3, iBloc+4) = -(hk*sk*radiusElbow*sinfi+radiusSection*dhk*ck)/denr
            b(3, iBloc+5) = (hk*ck*radiusElbow*sinfi-radiusSection*dhk*sk)/denr
            b(3, iBloc+6) = radiusElbow*cosfi*hk/denr

            b(4, iBloc+1) = sinfi*dhk*sk/denr
            b(4, iBloc+2) = -sinfi*dhk*ck/denr
            b(4, iBloc+3) = -cosfi*dhk/denr
            b(4, iBloc+4) = radiusElbow*cosfi*hk*sk/denr
            b(4, iBloc+5) = -radiusElbow*cosfi*hk*ck/denr
            b(4, iBloc+6) = radiusElbow*sinfi*hk/denr
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! beamNMatrStraight
!
! Get N matrix of shape functions and its transpose - Straight case
!
! In  beamElemProp     : properties of beam
! In  nbDof            : number of DOF in element
! In  kpg              : current Gauss point
! In  phi              : coordinate of sub-point in section
! In  radiusSection    : radius of section (only circle)
! In  shiftDOF         : shift for B matrix
! In  ff               : shape functions on beam
! IO  nvec             : vector of discrete disp/rota
! IO  tnvec              transposed vector of discrete disp/rota
!
! --------------------------------------------------------------------------------------------------
    subroutine beamNMatrStraight(beamElemProp, nbDof, kpg, &
                                 phi, radiusSection, shiftDOF, &
                                 ff, nvec, tnvec)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(beamElem_Prop), intent(in) :: beamElemProp
        integer(kind=8), intent(in) :: nbDof, kpg
        real(kind=8), intent(in) :: phi, radiusSection
        integer(kind=8), intent(in) :: shiftDOF
        real(kind=8), intent(in) :: ff(*)
        real(kind=8), intent(inout) :: nvec(BEAM_NBDOF, nbDof), tnvec(nbDof, BEAM_NBDOF)
! ----- Local
        integer(kind=8) :: iNode, iBloc, iDisp, nbNode
        real(kind=8) :: hk, cosfi, sinfi
!   ------------------------------------------------------------------------------------------------
!
        nbNode = beamElemProp%nbNode
        cosfi = cos(phi)
        sinfi = sin(phi)

        do iNode = 1, nbNode
            hk = ff(nbNode*(kpg-1)+iNode)

! --------- Select column
            iBloc = (shiftDOF)*(iNode-1)
            ASSERT(iBloc+6 .le. nbDof)

! --------- Set values
            do iDisp = 1, 3
                nvec(iDisp, iBloc+iDisp) = hk
                tnvec(iBloc+iDisp, iDisp) = hk
            end do

            nvec(1, iBloc+2) = 0.d0
            nvec(1, iBloc+4) = 0.d0
            nvec(1, iBloc+5) = -hk*radiusSection*cosfi
            nvec(1, iBloc+6) = hk*radiusSection*sinfi

            nvec(2, iBloc+1) = 0.d0
            nvec(2, iBloc+4) = hk*radiusSection*cosfi
            nvec(2, iBloc+5) = 0.d0

            nvec(3, iBloc+3) = hk
            nvec(3, iBloc+4) = -hk*radiusSection*sinfi
            nvec(3, iBloc+5) = 0.d0

            tnvec(iBloc+2, 1) = 0.d0
            tnvec(iBloc+4, 1) = 0.d0
            tnvec(iBloc+5, 1) = -hk*radiusSection*cosfi
            tnvec(iBloc+6, 1) = hk*radiusSection*sinfi

            tnvec(iBloc+1, 2) = 0.d0
            tnvec(iBloc+4, 2) = hk*radiusSection*cosfi
            tnvec(iBloc+5, 2) = 0.d0

            tnvec(iBloc+3, 3) = hk
            tnvec(iBloc+4, 3) = -hk*radiusSection*sinfi
            tnvec(iBloc+5, 3) = 0.d0
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! beamNMatrCurved
!
! Get N matrix of shape functions and its transpose - Curved case
!
! In  beamElemProp     : properties of beam
! In  nbDof            : number of DOF in element
! In  npg              : number of Gauss point
! In  kpg              : index of Gauss point
! In  xpg              : coordinates of Gauss points
! In  phi              : coordinate of sub-point in section
! In  radiusSection    : radius of section (only circle)
! In  shiftDOF         : shift for B matrix
! In  ff               : shape functions on beam
! IO  nvec             : vector of discrete disp/rota
! IO  tnvec           :  transposed vector of discrete disp/rota
!
! --------------------------------------------------------------------------------------------------
    subroutine beamNMatrCurved(beamElemProp, nbDof, npg, &
                               kpg, xpg, &
                               phi, radiusSection, shiftDOF, &
                               ff, nvec, tnvec)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(beamElem_Prop), intent(in) :: beamElemProp
        integer(kind=8), intent(in) :: nbDof, npg, kpg
        real(kind=8), intent(in) :: xpg(npg)
        real(kind=8), intent(in) :: phi, radiusSection
        integer(kind=8), intent(in) :: shiftDOF
        real(kind=8), intent(in) :: ff(*)
        real(kind=8), intent(inout) :: nvec(BEAM_NBDOF, nbDof), tnvec(nbDof, BEAM_NBDOF)
! ----- Local
        integer(kind=8) :: iNode, iBloc, iDisp, nbNode
        real(kind=8) :: hk, cosfi, sinfi, ck, sk, tk(4)
        real(kind=8) :: thetaElbow
!   ------------------------------------------------------------------------------------------------
!
        nbNode = beamElemProp%nbNode
        tk = beamElemProp%tk
        thetaElbow = beamElemProp%thetaElbow
        cosfi = cos(phi)
        sinfi = sin(phi)

        do iNode = 1, nbNode
            hk = ff(nbNode*(kpg-1)+iNode)
            ck = cos((1.d0+xpg(kpg))*thetaElbow/2.d0-tk(iNode))
            sk = sin((1.d0+xpg(kpg))*thetaElbow/2.d0-tk(iNode))

! --------- Select column
            iBloc = (shiftDOF)*(iNode-1)
            ASSERT(iBloc+6 .le. nbDof)

! --------- Set values
            do iDisp = 1, 3
                nvec(iDisp, iBloc+iDisp) = hk*ck
                tnvec(iBloc+iDisp, iDisp) = hk*ck
            end do

            nvec(1, ibloc+2) = hk*sk
            nvec(1, ibloc+4) = hk*radiusSection*cosfi*sk
            nvec(1, ibloc+5) = -hk*radiusSection*cosfi*ck
            nvec(1, ibloc+6) = hk*radiusSection*sinfi

            nvec(2, ibloc+1) = -hk*sk
            nvec(2, ibloc+4) = hk*radiusSection*cosfi*ck
            nvec(2, ibloc+5) = hk*radiusSection*cosfi*sk

            nvec(3, ibloc+3) = hk
            nvec(3, ibloc+4) = -hk*radiusSection*sinfi*ck
            nvec(3, ibloc+5) = -hk*radiusSection*sinfi*sk

            tnvec(ibloc+2, 1) = hk*sk
            tnvec(ibloc+4, 1) = hk*radiusSection*cosfi*sk
            tnvec(ibloc+5, 1) = -hk*radiusSection*cosfi*ck
            tnvec(ibloc+6, 1) = hk*radiusSection*sinfi

            tnvec(ibloc+1, 2) = -hk*sk
            tnvec(ibloc+4, 2) = hk*radiusSection*cosfi*ck
            tnvec(ibloc+5, 2) = hk*radiusSection*cosfi*sk

            tnvec(ibloc+3, 3) = hk
            tnvec(ibloc+4, 3) = -hk*radiusSection*sinfi*ck
            tnvec(ibloc+5, 3) = -hk*radiusSection*sinfi*sk
        end do
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! beamNMatr
!
! Get N matrix of shape functions and its transpose
!
! In  beamElemProp     : properties of beam
! In  nbDof            : number of DOF in element
! In  npg              : number of Gauss point
! In  kpg              : index of Gauss point
! In  xpg              : coordinates of Gauss points
! In  phi              : coordinate of sub-point in section
! In  radiusSection    : radius of section (only circle)
! In  shiftDOF         : shift for B matrix
! In  ff               : shape functions on beam
! IO  nvec             : vector of discrete disp/rota
! IO  tnvec           :  transposed vector of discrete disp/rota
!
! --------------------------------------------------------------------------------------------------
    subroutine beamNMatr(beamElemProp, nbDof, npg, &
                         kpg, xpg, &
                         phi, radiusSection, shiftDOF, &
                         ff, nvec, tnvec)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        type(beamElem_Prop), intent(in) :: beamElemProp
        integer(kind=8), intent(in) :: nbDof, npg, kpg
        real(kind=8), intent(in) :: xpg(npg)
        real(kind=8), intent(in) :: phi, radiusSection
        integer(kind=8), intent(in) :: shiftDOF
        real(kind=8), intent(in) :: ff(*)
        real(kind=8), intent(inout) :: nvec(BEAM_NBDOF, nbDof), tnvec(nbDof, BEAM_NBDOF)
!   ------------------------------------------------------------------------------------------------
!
        if (beamElemProp%beamType .eq. BEAM_TYPE_ELBOW) then
            call beamNMatrCurved(beamElemProp, nbDof, npg, &
                                 kpg, xpg, &
                                 phi, radiusSection, shiftDOF, &
                                 ff, nvec, tnvec)

        else if (beamElemProp%beamType .eq. BEAM_TYPE_STRAIGHT) then
            call beamNMatrStraight(beamElemProp, nbDof, kpg, &
                                   phi, radiusSection, shiftDOF, &
                                   ff, nvec, tnvec)
        else
            ASSERT(ASTER_FALSE)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module beamElem_module
