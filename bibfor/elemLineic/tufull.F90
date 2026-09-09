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
! aslint: disable=W1306
!
subroutine tufull(materPara, option, nbFourier, nbDof)
!
    use beamElem_type
    use Behaviour_module
    use Behaviour_type
    use MaterialPara_module
    use MaterialPara_type
    use pipeElem_module
    use pipeElem_type
    implicit none
!
#include "asterc/r8pi.h"
#include "asterf_types.h"
#include "asterfort/assert.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/elrefe_info.h"
#include "asterfort/jevech.h"
#include "asterfort/kcoude.h"
#include "asterfort/mavec.h"
#include "asterfort/nmcomp.h"
#include "asterfort/pipeElem_type.h"
#include "asterfort/tecach.h"
#include "jeveux.h"
!
    type(Material_Para), intent(inout) :: materPara
    character(len=16), intent(in) :: option
    integer(kind=8), intent(in) :: nbDof, nbFourier
!
! --------------------------------------------------------------------------------------------------
!
! Compute non-linear options for pipe element
!
! In  option           : option to compute
! In  nbDof            : number of DOF in element
! In  nbFourier        : number of Fourier modes
!
! --------------------------------------------------------------------------------------------------
!
    character(len=16), parameter :: multComp = " "
    integer(kind=8), parameter :: ndimLdc = 2
    character(len=8), parameter :: fami = "RIGI"
    real(kind=8), parameter :: rac2 = sqrt(2.d0)
    real(kind=8) :: dispPrev(nbDof), dispIncr(nbDof)
    integer(kind=8) :: jvf, jdfde, jdfd2, jcoopg, jpoids
    real(kind=8) :: radiusLayer
    real(kind=8) :: poids, weightLayer(2*PIPE_MAX_LAYERS+1), weightSect(2*PIPE_MAX_SECTORS+1)
    real(kind=8) :: cisail, gxz
    real(kind=8) :: jacobi, xpg(PIPE_MAX_NPG)
    real(kind=8) :: phi, zeta
    integer(kind=8) :: nbLayer, nbSect
    integer(kind=8) :: nspgLayer, nspgSect, npg, nspg
    integer(kind=8) :: kspgLayer, kspgSect, kpg, ksp, kspTot
    integer(kind=8) :: nbNode
    integer(kind=8) :: iDof, iSief, iTens, nc, iret, k2
    real(kind=8) :: b(PIPE_TENS_SIZE, nbDof)
    real(kind=8) :: ktild(nbDof, nbDof), effint(nbDof)
    real(kind=8) :: effinb, epsi(4), depsi(4), eps2d(6), deps2d(6)
    real(kind=8) :: sigmPrepInte(6), sigmPostInte(6), sgmtd(4)
    real(kind=8) :: dsidep(6, 6), dtild(PIPE_TENS_SIZE, PIPE_TENS_SIZE)
    real(kind=8) :: instm, instp
    integer(kind=8) :: ivarix
    integer(kind=8) :: nbvari, lgpg, jtab(7)
    integer(kind=8) :: imatuu
    integer(kind=8) :: ivarip, ivarim, icontm, icontp, ivectu, jcret
    integer(kind=8) :: iinstm, iinstp, jvDispPrev, jvDispIncr, jvCarcri
    integer(kind=8) :: codret, cod
    character(len=16), pointer :: compor(:) => null()
    character(len=16) :: defoComp
    aster_logical :: lVect, lMatr, lVari, lSigm
    integer(kind=8) :: variLen
    real(kind=8), allocatable :: varip(:)
    character(len=8), parameter :: typmod(2) = (/'C_PLAN  ', '        '/)
    type(pipeElem_Prop) :: pipeElem
    type(sectPipe_Prop) :: sectPipe
    type(beamElem_Prop) :: beamElem
    type(Behaviour_Integ) :: BEHInteg
!
! --------------------------------------------------------------------------------------------------
!
    call elrefe_info(fami=fami, &
                     nno=nbNode, npg=npg, &
                     jpoids=jpoids, jcoopg=jcoopg, jvf=jvf, jdfde=jdfde, jdfd2=jdfd2)
    ASSERT(npg .le. PIPE_MAX_NPG)
    nc = nbDof*(nbDof+1)/2
    codret = 0

! - Get input fields
    call jevech('PVARIMR', 'L', ivarim)
    call jevech('PINSTMR', 'L', iinstm)
    call jevech('PINSTPR', 'L', iinstp)
    call jevech('PDEPLMR', 'L', jvDispPrev)
    call jevech('PDEPLPR', 'L', jvDispIncr)
    call jevech('PCONTMR', 'L', icontm)
    call tecach('OOO', 'PVARIMR', 'L', iret, nval=7, itab=jtab)
    lgpg = max(jtab(6), 1)*jtab(7)

! - Get time
    instm = zr(iinstm)
    instp = zr(iinstp)

! - Initialisation of behaviour datastructure
    call behaviourInit(BEHInteg)

! - Get fields for non-linear behaviour
    call jevech('PCARCRI', 'L', jvCarcri)
    call jevech('PCOMPOR', 'L', vk16=compor)

! - Properties of behaviour
    read (compor(NVAR), '(I16)') nbvari
    defoComp = compor(DEFO)
    ASSERT(defoComp .eq. 'PETIT')

! - Set main parameters for behaviour (on cell)
    call behaviourSetParaCell(typmod, option, &
                              compor, zr(jvCarcri), &
                              instm, instp, &
                              materPara, BEHInteg)

! - Select objects to construct from option name
    call behaviourOption(option, compor, &
                         lMatr, lVect, &
                         lVari, lSigm, &
                         codret)

! - Get output fields
    if (lMatr) then
        call jevech('PMATUUR', 'E', imatuu)
    end if
    if (lVect) then
        call jevech('PVECTUR', 'E', ivectu)
    end if
    if (lSigm) then
        call jevech('PCONTPR', 'E', icontp)
        call jevech('PCODRET', 'E', jcret)
    end if
!
    variLen = npg*lgpg
    allocate (varip(variLen))
    if (lVari) then
        call jevech('PVARIMP', 'L', ivarix)
        varip = zr(ivarix:ivarix+variLen-1)
    else
        varip = zr(ivarim:ivarim+variLen-1)
    end if

! - Get parameters about layers and sections
    call pipeGetSubPoints(nbLayer, nbSect, &
                          nspgLayer, nspgSect, nspg, &
                          weightLayer, weightSect)

! - Get properties of pipe
    call pipeGetProperties(nbNode, nbFourier, pipeElem)
    beamElem = pipeElem%beamElem
    sectPipe = beamElem%sectPipe

! - For DEBUG
    ! if (PIPE_DEBUG) then
    !     call pipeElemDebug(pipeElem)
    ! end if

! - Get coordinates of Gauss points (on segment)
    do kpg = 1, npg
        xpg(kpg) = zr(jcoopg-1+kpg)
    end do

! - Get displacements (in local base)
    call pipeGetDisp(pipeElem, nbDof, zr(jvDispPrev), dispPrev)
    call pipeGetDisp(pipeElem, nbDof, zr(jvDispIncr), dispIncr)

! - Loop on Gauss points (on segment)
    kspTot = 0
    ktild = 0.d0
    effint = 0.d0
    do kpg = 1, npg
! ----- Loop on sub-points in thickness
        b = 0.d0
        do kspgLayer = 1, nspgLayer
! --------- Radius of layer
            call pipeGetRadiusLayer(pipeElem, kspgLayer, nbLayer, radiusLayer)

! --------- Coordinate of sub-point in thickness
            zeta = (kspgLayer-1)*pipeElem%beamElem%sectPipe%thickness/(2.d0*nbLayer)- &
                   pipeElem%beamElem%sectPipe%thickness/2.d0

! --------- Loop on sub-points in section
            do kspgSect = 1, nspgSect
                ksp = (kspgLayer-1)*(2*nbSect+1)+kspgSect
                kspTot = kspTot+1
                k2 = lgpg*(kpg-1)+((2*nbSect+1)*(kspgLayer-1)+(kspgSect-1))*nbvari

! ------------- Set main parameters for behaviour (on point)
                call behaviourSetParaPoin(kpg, ksp, BEHInteg)

! ------------- Coordinate of sub-point in section
                phi = (kspgSect-1)*2.d0*r8pi()/(2.d0*nbSect)

! ------------- Compute B matrix (u => epsi)
                call pipeBMatr(pipeElem, nbDof, &
                               npg, xpg, &
                               phi, zeta, radiusLayer, &
                               kpg, zr(jvf), zr(jdfde), zr(jdfd2), &
                               b)

! ------------- Jacobian
                call pipeGetJacobian(pipeElem, nbSect, nbLayer, phi, radiusLayer, &
                                     jacobi)

! ------------- Weighted jacobian of current sub-point
                poids = zr(jpoids-1+kpg)* &
                        weightLayer(kspgLayer)*weightSect(kspgSect)* &
                        jacobi

! ------------- Compute strains
                epsi = 0.d0
                do iTens = 1, 4
                    do iDof = 1, nbDof
                        epsi(iTens) = epsi(iTens)+b(iTens, iDof)*dispPrev(iDof)
                    end do
                end do
                eps2d = 0
                eps2d(1) = epsi(1)
                eps2d(2) = epsi(2)
                eps2d(3) = 0.d0
                eps2d(4) = epsi(3)/rac2

                depsi = 0.d0
                do iTens = 1, 4
                    do iDof = 1, nbDof
                        depsi(iTens) = depsi(iTens)+b(iTens, iDof)*dispIncr(iDof)
                    end do
                end do
                deps2d = 0
                deps2d(1) = depsi(1)
                deps2d(2) = depsi(2)
                deps2d(3) = 0.d0
                deps2d(4) = depsi(3)/rac2
                gxz = epsi(4)+depsi(4)

! ------------- Prepare stress
                sigmPrepInte = 0.d0
                do iSief = 1, 3
                    sigmPrepInte(iSief) = zr(icontm-1+6*(kspTot-1)+iSief)
                end do
                sigmPrepInte(4) = zr(icontm-1+6*(kspTot-1)+4)*rac2

! ------------- Integrate
                sigmPostInte = 0.d0
                call nmcomp(BEHInteg, &
                            ndimLdc, option, typmod, &
                            instm, instp, &
                            compor, zr(jvCarcri), multComp, &
                            6, eps2d, deps2d, &
                            6, sigmPrepInte, zr(ivarim+k2), &
                            sigmPostInte, varip(1+k2), &
                            36, dsidep, cod)

! ------------- Get shear parameter from elasticity
                call pipeGetElasProp(BEHINteg%materPara%jvMaterCode, &
                                     kpg_=kpg, kspg_=ksp, &
                                     cisail_=cisail)

                if (cod .ne. 0) then
                    if (codret .ne. 1) then
                        codret = cod
                    end if
                end if

! ------------- Compute tangent matrix
                if (lMatr) then
                    dtild(1, 1) = dsidep(1, 1)
                    dtild(1, 2) = dsidep(1, 2)
                    dtild(1, 3) = dsidep(1, 4)/rac2
                    dtild(1, 4) = 0.d0
!
                    dtild(2, 1) = dsidep(2, 1)
                    dtild(2, 2) = dsidep(2, 2)
                    dtild(2, 3) = dsidep(2, 4)/rac2
                    dtild(2, 4) = 0.d0
!
                    dtild(3, 1) = dsidep(4, 1)/rac2
                    dtild(3, 2) = dsidep(4, 2)/rac2
                    dtild(3, 3) = dsidep(4, 4)/2.d0
                    dtild(3, 4) = 0.d0
!
                    dtild(4, 1) = 0.d0
                    dtild(4, 2) = 0.d0
                    dtild(4, 3) = 0.d0
                    dtild(4, 4) = cisail/2.d0
!
                    call kcoude(nbDof, poids, b, dtild, ktild)
                end if
!
                if (lSigm) then
                    do iSief = 1, 3
                        zr(icontp-1+6*(kspTot-1)+iSief) = sigmPostInte(iSief)
                    end do
                    zr(icontp-1+6*(kspTot-1)+4) = sigmPostInte(4)/rac2
                    zr(icontp-1+6*(kspTot-1)+5) = cisail*gxz/2.d0
                    zr(icontp-1+6*(kspTot-1)+6) = 0.d0
                end if
                if (lVect) then
                    ASSERT(lSigm)
                    sgmtd(1) = zr(icontp-1+6*(kspTot-1)+1)
                    sgmtd(2) = zr(icontp-1+6*(kspTot-1)+2)
                    sgmtd(3) = zr(icontp-1+6*(kspTot-1)+4)
                    sgmtd(4) = cisail*gxz/2.d0
                    do iDof = 1, nbDof
                        effinb = 0.d0
                        do iTens = 1, PIPE_TENS_SIZE
                            effinb = effinb+b(iTens, iDof)*sgmtd(iTens)
                        end do
                        effint(iDof) = effint(iDof)+poids*effinb
                    end do
                end if
            end do
        end do
    end do

! - Save matrix
    if (lMatr) then
        call pipeBaseForMatr(pipeElem, nbDof, ktild)
        call mavec(ktild, nbDof, zr(imatuu), nc)
    end if

! - Save internal forces
    if (lVect) then
        call pipeBaseForVect("LG", pipeElem, nbDof, effint)
        do iDof = 1, nbDof
            zr(ivectu-1+iDof) = effint(iDof)
        end do
    end if

! - Save internal state variables
    if (lVari) then
        call jevech('PVARIPR', 'E', ivarip)
        zr(ivarip:ivarip+variLen-1) = varip(1:variLen)
    end if

! - Save behaviour integration code
    if (lSigm) then
        zi(jcret) = codret
    end if
!
    deallocate (varip)
!
end subroutine
