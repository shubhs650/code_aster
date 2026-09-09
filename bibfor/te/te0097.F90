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
subroutine te0097(option, nomte)
!
    use MaterialPara_module
    use MaterialPara_type
    implicit none
!
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/elref2.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/lteatt.h"
#include "asterfort/nbsigm.h"
#include "asterfort/niinit.h"
#include "asterfort/sigvmc.h"
#include "jeveux.h"
#include "MeshTypes_type.h"
!
    character(len=16), intent(in) :: option, nomte
!
! --------------------------------------------------------------------------------------------------
!
! Elementary computation
!
! Elements: *INCO*
!
! Options: SIEF_ELGA
!
! --------------------------------------------------------------------------------------------------
!
    character(len=8), parameter :: fami = 'RIGI'
    real(kind=8), parameter :: nharm = 0.d0
    integer(kind=8), parameter :: npgMax = 27
    integer(kind=8) :: ndim, npg, jvWeightDisp, jvShapeDisp, jvDShapeDisp, jvShapePres
    integer(kind=8) :: idim
    integer(kind=8) :: jvSigm, jvDisp, jvGeom, jvMaterc, nbsig
    integer(kind=8) :: kpg, isig, iNodeDisp, iNodePres, iOSGS
    integer(kind=8) :: nbNodeDisp, nbNodePres, nbNodeGonf
    real(kind=8) :: sigmDisp(npgMax*6), time, sigmTrac(npgMax)
    real(kind=8) :: presGaus(npgMax), dispU(3*MT_NNOMAX), dispP(MT_NNOMAX)
    integer(kind=8) :: vu(3, MT_NNOMAX), vg(MT_NNOMAX), vp(MT_NNOMAX), vpi(3, MT_NNOMAX)
    integer(kind=8) :: nbElrefe, iRefePres, iRefeGonf
    character(len=8) :: listElrefe(10), typmod(2)
    aster_logical :: lGonf
    type(Material_Para) :: materPara
!
! --------------------------------------------------------------------------------------------------
!
    nbsig = nbsigm()
    ASSERT(nbsig .le. 6)

! - Get current time
    time = r8vide()

! - Some unknonws
    if (lteatt('INCO', 'C2 ')) then
        lGonf = ASTER_FALSE
        iOSGS = 0
    elseif (lteatt('INCO', 'C2O')) then
        lGonf = ASTER_FALSE
        iOSGS = 1
    elseif (lteatt('INCO', 'C3 ')) then
        lGonf = ASTER_TRUE
        iOSGS = 0
    else
        ASSERT(ASTER_FALSE)
    end if

! - List of ELREFE
    call elref2(nomte, 10, listElrefe, nbElrefe)
    ASSERT(nbElrefe .ge. 2)
    if (lGonf) then
        iRefePres = 3
        iRefeGonf = 2
    else
        iRefePres = 2
        iRefeGonf = 0
    end if

! - Get paramers of finite element for displacements
    call elrefe_info(elrefe=listElrefe(1), fami=fami, &
                     ndim=ndim, nno=nbNodeDisp, npg=npg, &
                     jpoids=jvWeightDisp, jvf=jvShapeDisp, jdfde=jvDShapeDisp)
    ASSERT(npg .le. npgMax)
    ASSERT(nbNodeDisp .le. MT_NNOMAX)

! - Get paramers of finite element for pres
    call elrefe_info(elrefe=listElrefe(iRefePres), fami=fami, &
                     nno=nbNodePres, &
                     jvf=jvShapePres)
    ASSERT(nbNodePres .le. MT_NNOMAX)

! - Get paramers of finite element for gonf
    nbNodeGonf = 0
    if (iRefeGonf .ne. 0) then
        call elrefe_info(elrefe=listElrefe(iRefeGonf), fami=fami, &
                         nno=nbNodeGonf)
    end if
    ASSERT(nbNodeGonf .le. MT_NNOMAX)

! - Modelling
    typmod = ' '
    if (ndim .eq. 2 .and. lteatt('AXIS', 'OUI')) then
        typmod(1) = 'AXIS'
    else if (ndim .eq. 2 .and. lteatt('D_PLAN', 'OUI')) then
        typmod(1) = 'D_PLAN'
    else if (ndim .eq. 3) then
        typmod(1) = '3D'
    else
        ASSERT(ASTER_FALSE)
    end if

! - Get index of dof
    call niinit(typmod, &
                ndim, nbNodeDisp, nbNodeGonf, nbNodePres, iOSGS, &
                vu, vg, vp, vpi)

! - Geometry
    call jevech('PGEOMER', 'L', jvGeom)

! - Current displacements (nodes)
    call jevech('PDEPLAR', 'L', jvDisp)

! - Get material parameters
    call jevech('PMATERC', 'L', jvMaterc)

! - Initializations of material parameters on current cell
    call initParaCell(fami, zi(jvMaterc), materPara)

! - Set local coordinate system from user
    call getUserLCS(ndim, nbNodeDisp, jvGeom, materPara%lcsPara)

! - Get displacements for u, v, w
    dispU = 0.d0
    do iNodeDisp = 1, nbNodeDisp
        do idim = 1, ndim
            dispU(idim+ndim*(iNodeDisp-1)) = zr(jvDisp-1+vu(iDim, iNodeDisp))
        end do
    end do

! - Get displacements for p
    dispP = 0.d0
    do iNodePres = 1, nbNodePres
        dispP(iNodePres) = zr(jvDisp-1+vp(iNodePres))
    end do

! - Compute mechanical stress (without effect of external state variables)
    sigmDisp = 0.d0
    call sigvmc(materPara, &
                nbNodeDisp, ndim, nbsig, npg, &
                jvWeightDisp, jvShapeDisp, jvDShapeDisp, &
                zr(jvGeom), dispU, &
                time, nharm, &
                sigmDisp)

    do kpg = 1, npg
        presGaus(kpg) = 0.d0
        do iNodePres = 1, nbNodePres
            presGaus(kpg) = presGaus(kpg)+ &
                            zr(jvShapePres-1+nbNodePres*(kpg-1)+iNodePres)*dispP(iNodePres)
        end do
        sigmTrac(kpg) = 0.d0
        do isig = 1, 3
            sigmTrac(kpg) = sigmTrac(kpg)+sigmDisp(nbsig*(kpg-1)+isig)
        end do
    end do

! - Output field
    call jevech('PCONTRR', 'E', jvSigm)

    do kpg = 1, npg
        do isig = 1, nbsig+1
            if (isig .le. 3) then
                zr(jvSigm+(nbsig+1)*(kpg-1)+isig-1) = sigmDisp(nbsig*(kpg-1)+isig)- &
                                                      sigmTrac(kpg)/3.d0+presGaus(kpg)
            elseif (isig .le. nbsig) then
                zr(jvSigm+(nbsig+1)*(kpg-1)+isig-1) = sigmDisp(nbsig*(kpg-1)+isig)
            else
                zr(jvSigm+(nbsig+1)*(kpg-1)+isig-1) = sigmTrac(kpg)/3.d0-presGaus(kpg)
            end if
        end do
    end do
!
end subroutine
