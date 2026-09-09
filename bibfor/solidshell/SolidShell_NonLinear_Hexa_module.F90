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
module SolidShell_NonLinear_Hexa_module
! ==================================================================================================
    use SolidShell_type
    use SolidShell_Debug_module
    use SolidShell_Geometry_Hexa_module
    use SolidShell_Elementary_Hexa_module
    use SolidShell_Kinematic_Hexa_module
    use SolidShell_Stabilization_Hexa_module
    use SolidShell_Utilities_module
    use Behaviour_module
! ==================================================================================================
    implicit none
! ==================================================================================================
    public :: compNonLinearHexa
    private :: compSmallStrainHexa, compGdefLogHexa, postLog
! ==================================================================================================
    private
#include "jeveux.h"
#include "asterf_types.h"
#include "MeshTypes_type.h"
#include "asterfort/Behaviour_type.h"
#include "asterfort/SolidShell_type.h"
#include "asterc/isnnem.h"
#include "asterfort/assert.h"
#include "asterfort/btsig.h"
#include "asterfort/codere.h"
#include "asterfort/deflg2.h"
#include "asterfort/deflg3.h"
#include "asterfort/jevech.h"
#include "asterfort/nmcomp.h"
#include "asterfort/symt46.h"
#include "asterfort/tecach.h"
#include "blas/daxpy.h"
#include "blas/dcopy.h"
! ==================================================================================================
contains
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!
! compNonLinearHexa
!
! Compute non-linear options for HEXA
!
! In  option           : name of option to compute
! In  elemProp         : general properties of element
! In  cellGeom         : general geometric properties of cell
! IO  behaPara         : parameters of behaviour
!
! --------------------------------------------------------------------------------------------------
    subroutine compNonLinearHexa(option, elemProp, cellGeom, behaPara)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=16), intent(in) :: option
        type(SSH_ELEM_PROP), intent(in) :: elemProp
        type(SSH_CELL_GEOM), intent(in) :: cellGeom
        type(SSH_BEHA_PARA), intent(inout) :: behaPara
! ----- Locals
        integer(kind=8) :: nbIntePoint, nbVari
        integer(kind=8) :: jtab(7), iret
        integer(kind=8) :: jvGeom
        integer(kind=8) :: jvMatr, jvVect, jvSigmP, jvVariP, jvVariX, jvCodret
        integer(kind=8) :: jvTimeM, jvTimeP, jvSigmM, jvVariM, jvDispM, jvDispIncr
        blas_int :: b_incx, b_incy, b_n
!   ------------------------------------------------------------------------------------------------
!
!
! - Properties of finite element
        nbIntePoint = elemProp%elemInte%nbIntePoint
!
! - Get input fields
        jvGeom = cellGeom%jvGeom
        call jevech('PINSTMR', 'L', jvTimeM)
        call jevech('PINSTPR', 'L', jvTimeP)
        call jevech('PCONTMR', 'L', jvSigmM)
        call jevech('PDEPLMR', 'L', jvDispM)
        call jevech('PDEPLPR', 'L', jvDispIncr)
        call tecach('OOO', 'PVARIMR', 'L', iret, nval=7, &
                    itab=jtab)
        jvVariM = jtab(1)
        nbVari = max(jtab(6), 1)*jtab(7)
!
! - Get output fields
        jvMatr = isnnem()
        if (behaPara%lMatr) then
            if (behaPara%lMatrSyme) then
                call jevech('PMATUUR', 'E', jvMatr)
            else
                call jevech('PMATUNS', 'E', jvMatr)
            end if
        end if
        jvVect = isnnem()
        if (behaPara%lVect) then
            call jevech('PVECTUR', 'E', jvVect)
        end if
        jvSigmP = isnnem()
        if (behaPara%lSigm) then
            call jevech('PCONTPR', 'E', jvSigmP)
            call jevech('PCODRET', 'E', jvCodret)
        end if
        jvVariP = isnnem()
        if (behaPara%lVari) then
            call jevech('PVARIPR', 'E', jvVariP)
            call jevech('PVARIMP', 'L', jvVariX)
            b_n = to_blas_int(nbIntePoint*nbVari)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, zr(jvVariX), b_incx, zr(jvVariP), b_incy)
        end if
!
! - Compute
        if (behaPara%defoComp .eq. 'PETIT') then
            call compSmallStrainHexa(option, elemProp, cellGeom, behaPara, &
                                     nbIntePoint, nbVari, zr(jvTimeM), zr(jvTimeP), zr(jvDispM), &
                                     zr(jvDispIncr), zr(jvSigmM), zr(jvVariM), zr(jvSigmP), &
                                     zr(jvVariP), zr(jvMatr), zr(jvVect), zi(jvCodret))
!
        else if (behaPara%defoComp .eq. 'GDEF_LOG') then
            call compGdefLogHexa(option, elemProp, cellGeom, behaPara, &
                                 nbIntePoint, nbVari, zr(jvTimeM), zr(jvTimeP), zr(jvDispM), &
                                 zr(jvDispIncr), zr(jvSigmM), zr(jvVariM), zr(jvSigmP), &
                                 zr(jvVariP), zr(jvMatr), zr(jvVect), zi(jvCodret))
!
        else
            ASSERT(ASTER_FALSE)
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compSmallStrainHexa
!
! Compute non-linear options for HEXA - Small strains
!
! In  option           : name of option to compute
! In  elemProp         : general properties of element
! In  cellGeom         : general geometric properties of cell
! IO  behaPara         : parameters of behaviour
! In  nbIntePoint      : number of integration points on cell
! In  nbVari           : number of internal state variables
! In  timePrev         : time of previous time step
! In  timeCurr         : time of current time step
! In  dispPrev         : displacement at beginning of current time step
! In  dispIncr         : displacement since beginning of current time step
! In  sigm             : stress at beginning of current time step (before integration)
! In  vim              : internal state vari. at beginning of current time step (before integration)
! Out sigp             : stress at end of current time step (after integration)
! Out vip              : internal state vari. at end of current time step (after integration)
! Out matr             : tangent matrix
! Out vect             : vector of internal forces
! Out codret           : error code from integration of behaviour
!
! --------------------------------------------------------------------------------------------------
    subroutine compSmallStrainHexa(option, elemProp, cellGeom, behaPara, &
                                   nbIntePoint, nbVari, timePrev, timeCurr, dispPrev, &
                                   dispIncr, sigm, vim, sigp, vip, &
                                   matr, vect, codret)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=16), intent(in) :: option
        type(SSH_ELEM_PROP), intent(in) :: elemProp
        type(SSH_CELL_GEOM), intent(in) :: cellGeom
        type(SSH_BEHA_PARA), intent(inout) :: behaPara
        integer(kind=8), intent(in) :: nbIntePoint, nbVari
        real(kind=8), intent(in) :: timePrev, timeCurr
        real(kind=8), intent(in) :: dispPrev(SSH_NBDOF_HEXA), dispIncr(SSH_NBDOF_HEXA)
        real(kind=8), intent(in) :: sigm(SSH_SIZE_TENS, nbIntePoint), vim(nbVari, nbIntePoint)
        real(kind=8), intent(out) :: sigp(SSH_SIZE_TENS, nbIntePoint), vip(nbVari, nbIntePoint)
        real(kind=8), intent(out) :: matr(*), vect(SSH_NBDOF_HEXA)
        integer(kind=8), intent(out) :: codret
! ----- Local
        character(len=16), parameter :: multComp = " "
        type(SSH_GEOM_HEXA) :: geomHexa
        integer(kind=8), parameter :: ksp = 1
        real(kind=8), parameter :: rac2 = sqrt(2.d0)
        integer(kind=8) :: cod(nbIntePoint), kpg, iTens, nbDof, nbDofGeom
        integer(kind=8) :: jvCoor, jvWeight
        type(SSH_KINE_HEXA) :: kineHexa
        type(SSH_STAB_HEXA) :: stabHexa
        integer(kind=8) :: i, j, ij
        real(kind=8) :: zeta, poids, jacob
        real(kind=8) :: UeffKpg, Ueff
        real(kind=8) :: dispCurr(SSH_NBDOF_HEXA)
        real(kind=8) :: epsiPrev(SSH_SIZE_TENS), epsiIncr(SSH_SIZE_TENS)
        real(kind=8) :: sigmPost(SSH_SIZE_TENS), sigmPrep(SSH_SIZE_TENS)
        real(kind=8), dimension(SSH_NBDOF_HEXA, SSH_NBDOF_HEXA) :: matrMate, matrTang
        real(kind=8), dimension(SSH_SIZE_TENS, SSH_SIZE_TENS) :: dsidep
!   ------------------------------------------------------------------------------------------------
!
        cod = 0
        if (behaPara%lSigm) then
            sigp = 0.d0
        end if
        if (behaPara%lVect) then
            vect = 0.d0
        end if

! ----- Properties of finite element
        nbDof = elemProp%nbDof
        nbDofGeom = elemProp%nbDofGeom
        jvCoor = elemProp%elemInte%jvCoor
        jvWeight = elemProp%elemInte%jvWeight

! ----- Properties of finite element
        ASSERT(.not. behaPara%lLarge)

! ----- Total displacement from initial configuration
        dispCurr = dispIncr+dispPrev

! ----- Prepare geometric quantities
        call initGeomCellHexa(cellGeom, geomHexa)
        if (SSH_DBG_GEOM) call dbgObjGeomHexa(geomHexa)

! ----- Compute gradient matrix in covariant basis
        kineHexa%lLarge = ASTER_FALSE
        call compBCovaMatrHexa(geomHexa, kineHexa)

! ----- Compute gradient matrix in cartesian frame
        call compBCartMatrHexa(geomHexa, kineHexa)
        if (SSH_DBG_KINE) call dbgObjKineHexa(kineHexa, smallCstPart_=ASTER_TRUE)

! ----- Loop on Gauss points
        Ueff = 0.d0
        matrTang = 0.d0
        do kpg = 1, nbIntePoint
            zeta = zr(jvCoor-1+3*kpg)
            poids = zr(jvWeight-1+kpg)
            jacob = poids*cellGeom%detJac0

! --------- Compute EAS B matrix in cartesian frame at current Gauss point
            call compBCartEASMatrHexa(zeta, geomHexa, kineHexa)

! --------- Compute B matrix
            call compBMatrHexa(zeta, kineHexa)
            if (SSH_DBG_KINE) call dbgObjKineHexa(kineHexa, smallVarPart_=ASTER_TRUE)

! --------- Compute small strains at beginning of time step
            call compEpsiHexa(kineHexa, dispPrev, epsiPrev)

! --------- Compute increment of small strains
            call compEpsiHexa(kineHexa, dispIncr, epsiIncr)

! --------- Pre-treatment of stresses and strains
            epsiPrev(4) = epsiPrev(4)/rac2
            epsiPrev(5) = epsiPrev(5)/rac2
            epsiPrev(6) = epsiPrev(6)/rac2
            epsiIncr(4) = epsiIncr(4)/rac2
            epsiIncr(5) = epsiIncr(5)/rac2
            epsiIncr(6) = epsiIncr(6)/rac2
            do iTens = 1, 3
                sigmPrep(iTens) = sigm(iTens, kpg)
                sigmPrep(iTens+3) = sigm(iTens+3, kpg)*rac2
            end do

! --------- Set main parameters for behaviour (on point)
            call behaviourSetParaPoin(kpg, ksp, behaPara%BEHInteg)

! --------- Integrator
            sigmPost = 0.d0
            dsidep = 0.d0
            cod(kpg) = 0
            call nmcomp(behaPara%BEHInteg, &
                        SSH_NDIM, option, typmod, &
                        timePrev, timeCurr, &
                        behaPara%compor, behaPara%carcri, multComp, &
                        SSH_SIZE_TENS, epsiPrev, epsiIncr, &
                        SSH_SIZE_TENS, sigmPrep, vim(1, kpg), &
                        sigmPost, vip(1, kpg), SSH_SIZE_TENS*SSH_SIZE_TENS, dsidep, cod(kpg))
            if (cod(kpg) .eq. 1) then
                goto 99
            end if

! --------- Post-treatment of stresses and matrix
            if (behaPara%lSigm) then
                sigmPost(4) = sigmPost(4)/rac2
                sigmPost(5) = sigmPost(5)/rac2
                sigmPost(6) = sigmPost(6)/rac2
            end if
            if (behaPara%lMatr) then
                dsidep(4:6, 4:6) = dsidep(4:6, 4:6)/2.d0
                dsidep(4:6, 1:3) = dsidep(4:6, 1:3)/rac2
                dsidep(1:3, 4:6) = dsidep(1:3, 4:6)/rac2
            end if

! --------- Compute effective shear modulus for stabilization
            call compStabModulusHexa(sigmPost, epsiIncr, dsidep, UeffKpg)
            Ueff = Ueff+UeffKpg*poids/8.d0

! --------- Compute material part  at current Gauss point
            if (behaPara%lMatr) then
                matrMate = 0.d0
                call prodBTDB(dsidep, SSH_SIZE_TENS, elemProp%nbDof, kineHexa%B, matrMate)
            end if

! --------- Update tangent matrix
            if (behaPara%lMatr) then
                matrTang = matrTang+jacob*matrMate
            end if

! --------- Update internal force at current Gauss point
            if (behaPara%lVect) then
                call btsig(elemProp%nbDof, SSH_SIZE_TENS, jacob, kineHexa%B, sigmPost, &
                           vect)
            end if

! --------- Save stresses
            if (behaPara%lSigm) then
                do iTens = 1, SSH_SIZE_TENS
                    sigp(iTens, kpg) = sigmPost(iTens)
                end do
            end if
!
        end do

! ----- Stabilization
        if (behaPara%lVect) then
            call compStabSigmHexa(geomHexa, kineHexa, Ueff, dispCurr, stabHexa)
        end if
        if (behaPara%lMatr) then
            call compStabMatrMateHexa(geomHexa, kineHexa, Ueff, stabHexa)
        end if
        if (behaPara%lVect) then
            call compStabForcHexa(kineHexa, stabHexa)
        end if

! ----- Save matrix and vector
        if (behaPara%lMatr) then
            matrTang(1:nbDofGeom, 1:nbDofGeom) = matrTang(1:nbDofGeom, 1:nbDofGeom &
                                                          )+stabHexa%matrStabMate
        end if
        if (behaPara%lVect) then
            vect(1:SSH_NBDOFG_HEXA) = vect(1:SSH_NBDOFG_HEXA)+cellGeom%detJac0*stabHexa%forcStab
        end if
!
99      continue

! ----- Return code summary
        if (behaPara%lSigm) then
            call codere(cod, nbIntePoint, codret)
        end if

! ----- Write matrix
        if (behaPara%lMatr) then
            if (behaPara%lMatrSyme) then
                do i = 1, SSH_NBDOF_HEXA
                    do j = 1, i
                        ij = (i-1)*i/2+j
                        matr(ij) = matrTang(i, j)
                    end do
                end do
            else
                do j = 1, SSH_NBDOF_HEXA
                    do i = 1, SSH_NBDOF_HEXA
                        ij = j+(i-1)*SSH_NBDOF_HEXA
                        matr(ij) = matrTang(i, j)
                    end do
                end do
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! compGdefLogHexa
!
! Compute non-linear options for HEXA - GDEF_LOG
!
! In  option           : name of option to compute
! In  elemProp         : general properties of element
! In  cellGeom         : general geometric properties of cell
! IO  behaPara         : parameters of behaviour
! In  nbIntePoint      : number of integration points on cell
! In  nbVari           : number of internal state variables
! In  timePrev         : time of previous time step
! In  timeCurr         : time of current time step
! In  dispPrev         : displacement at beginning of current time step
! In  dispIncr         : displacement since beginning of current time step
! In  sigm             : stress at beginning of current time step (before integration)
! In  vim              : internal state vari. at beginning of current time step (before integration)
! Out sigp             : stress at end of current time step (after integration)
! Out vip              : internal state vari. at end of current time step (after integration)
! Out matr             : tangent matrix
! Out vect             : vector of internal forces
! Out codret           : error code from integration of behaviour
!
! --------------------------------------------------------------------------------------------------
    subroutine compGdefLogHexa(option, elemProp, cellGeom, behaPara, &
                               nbIntePoint, nbVari, timePrev, timeCurr, dispPrev, &
                               dispIncr, sigm, vim, sigp, vip, &
                               matr, vect, codret)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        character(len=16), intent(in)   :: option
        type(SSH_ELEM_PROP), intent(in) :: elemProp
        type(SSH_CELL_GEOM), intent(in) :: cellGeom
        type(SSH_BEHA_PARA), intent(inout) :: behaPara
        integer(kind=8), intent(in) :: nbIntePoint, nbVari
        real(kind=8), intent(in) :: timePrev, timeCurr
        real(kind=8), intent(in) :: dispPrev(SSH_NBDOF_HEXA), dispIncr(SSH_NBDOF_HEXA)
        real(kind=8), intent(in) :: sigm(SSH_SIZE_TENS, nbIntePoint), vim(nbVari, nbIntePoint)
        real(kind=8), intent(out) :: sigp(SSH_SIZE_TENS, nbIntePoint), vip(nbVari, nbIntePoint)
        real(kind=8), intent(out) :: matr(*), vect(SSH_NBDOF_HEXA)
        integer(kind=8), intent(out) :: codret
! ----- Local
        character(len=16), parameter :: multComp = " "
        type(SSH_GEOM_HEXA) :: geomHexa
        integer(kind=8), parameter :: ksp = 1
        integer(kind=8) :: cod(nbIntePoint), kpg, iTens, nbDof, nbDofGeom
        integer(kind=8) :: jvCoor, jvWeight
        aster_logical :: lVect, lMatr, lSigm, lVari, lMatrPred
        type(SSH_KINE_HEXA) :: kineHexa
        type(SSH_STAB_HEXA) :: stabHexa
        integer(kind=8) :: i, j, ij
        real(kind=8) :: zeta, poids, jacob
        real(kind=8) :: UeffKpg, Ueff
        real(kind=8) :: dispCurr(SSH_NBDOF_HEXA)
        real(kind=8) :: epslIncr(SSH_SIZE_TENS)
        real(kind=8) :: sigmPrep(SSH_SIZE_TENS), pk2(SSH_SIZE_TENS)
        real(kind=8) :: tPrev(SSH_SIZE_TENS), tCurr(SSH_SIZE_TENS)
        real(kind=8), dimension(SSH_NBDOF_HEXA, SSH_NBDOF_HEXA) :: matrMate, matrGeom, matrTang
        real(kind=8), dimension(SSH_SIZE_TENS, SSH_SIZE_TENS) :: dtde, dsidep
        blas_int :: b_incx, b_incy, b_n
!   ------------------------------------------------------------------------------------------------
!
        cod = 0
        if (behaPara%lSigm) then
            sigp = 0.d0
        end if
        if (behaPara%lVect) then
            vect = 0.d0
        end if
        lSigm = L_SIGM(option)
        lVect = L_VECT(option)
        lMatr = L_MATR(option)
        lVari = L_VARI(option)
        lMatrPred = L_MATR_PRED(option)

! ----- Properties of finite element
        nbDof = elemProp%nbDof
        nbDofGeom = elemProp%nbDofGeom
        jvCoor = elemProp%elemInte%jvCoor
        jvWeight = elemProp%elemInte%jvWeight

! ----- Properties of finite element
        ASSERT(behaPara%lLarge)

! ----- Total displacement from initial configuration
        dispCurr = dispIncr+dispPrev

! ----- Prepare geometric quantities
        call initGeomCellHexa(cellGeom, geomHexa, dispCurr)
        if (SSH_DBG_GEOM) call dbgObjGeomHexa(geomHexa)

! ----- Compute gradient matrix in covariant basis
        kineHexa%lLarge = ASTER_TRUE
        call compBCovaMatrHexa(geomHexa, kineHexa)

! ----- Compute gradient matrix in cartesian frame
        call compBCartMatrHexa(geomHexa, kineHexa)
        if (SSH_DBG_KINE) call dbgObjKineHexa(kineHexa, smallCstPart_=ASTER_TRUE)

! ----- Loop on Gauss points
        Ueff = 0.d0
        matrTang = 0.d0
        do kpg = 1, nbIntePoint
            zeta = zr(jvCoor-1+3*kpg)
            poids = zr(jvWeight-1+kpg)
            jacob = poids*cellGeom%detJac0

! --------- Compute EAS B matrix in cartesian frame at current Gauss point
            call compBCartEASMatrHexa(zeta, geomHexa, kineHexa)

! --------- Compute B matrix
            call compBMatrHexa(zeta, kineHexa)

! --------- Compute Green-Lagrange strains at beginning of time step
            call compECovaMatrHexa(cellGeom, dispPrev, kineHexa%epsgPrev)
            call compEpsgHexa(zeta, geomHexa, kineHexa%epsgPrev)
            kineHexa%epsgPrev%vale = kineHexa%epsgPrev%vale+kineHexa%BCartEAS*dispPrev(25)
            if (SSH_DBG_KINE) call dbgObjEpsgHexa(kineHexa%epsgPrev)

! --------- Compute logarithmic strains at beginning of time step
            call compEpslHexa(kineHexa%epsgPrev, kineHexa%epslPrev, cod(kpg))
            if (cod(kpg) .ne. 0) then
                goto 99
            end if
            if (SSH_DBG_KINE) call dbgObjEpslHexa(kineHexa%epslPrev)

! --------- Compute Green-Lagrange strains at end of time step
            call compECovaMatrHexa(cellGeom, dispCurr, kineHexa%epsgCurr)
            call compEpsgHexa(zeta, geomHexa, kineHexa%epsgCurr)
            kineHexa%epsgCurr%vale = kineHexa%epsgCurr%vale+kineHexa%BCartEAS*dispCurr(25)
            if (SSH_DBG_KINE) call dbgObjEpsgHexa(kineHexa%epsgCurr)

! --------- Compute logarithmic strains at end of time step
            call compEpslHexa(kineHexa%epsgCurr, kineHexa%epslCurr, cod(kpg))
            if (cod(kpg) .ne. 0) then
                goto 99
            end if
            if (SSH_DBG_KINE) call dbgObjEpslHexa(kineHexa%epslCurr)

! --------- Compute increment of strains
            epslIncr = kineHexa%epslCurr%vale-kineHexa%epslPrev%vale

! --------- Get "logarithmic" stresses from internal state variables at previous time step
            b_n = to_blas_int(2*SSH_NDIM)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call dcopy(b_n, vim(nbVari-6+1, kpg), b_incx, tPrev, b_incy)

! --------- Get Cauchy stresses
            do iTens = 1, 6
                sigmPrep(iTens) = sigm(iTens, kpg)
            end do

! --------- Set main parameters for behaviour (on point)
            call behaviourSetParaPoin(kpg, ksp, behaPara%BEHInteg)

! --------- Integrator
            tCurr = 0.d0
            dtde = 0.d0
            cod(kpg) = 0
            call nmcomp(behaPara%BEHInteg, &
                        SSH_NDIM, option, typmod, &
                        timePrev, timeCurr, &
                        behaPara%compor, behaPara%carcri, multComp, &
                        SSH_SIZE_TENS, kineHexa%epslPrev%vale, epslIncr, &
                        SSH_SIZE_TENS, tPrev, vim(1, kpg), &
                        tCurr, vip(1, kpg), SSH_SIZE_TENS*SSH_SIZE_TENS, dtde, cod(kpg))
            if (cod(kpg) .eq. 1) then
                goto 99
            end if

! --------- Post-treatment for logarithmic quantities
            call postLog(lMatrPred, lMatr, lSigm, &
                         kineHexa, tPrev, tCurr, dtde, &
                         dsidep, pk2)

! --------- Save "logarithmic" stresses in internal state variables
            if (lVari) then
                b_n = to_blas_int(2*SSH_NDIM)
                b_incx = to_blas_int(1)
                b_incy = to_blas_int(1)
                call dcopy(b_n, tCurr, b_incx, vip(nbVari-6+1, kpg), b_incy)
            end if

! --------- Compute effective shear modulus for stabilization
            call compStabModulusHexa(sigmPrep, kineHexa%epsgPrev%vale, dsidep, UeffKpg)
            Ueff = Ueff+UeffKpg*poids/8.d0

! --------- Compute material part of matrix at current Gauss point
            if (behaPara%lMatr) then
                matrMate = 0.d0
                call prodBTDB(dsidep, SSH_SIZE_TENS, elemProp%nbDof, kineHexa%B, matrMate)
            end if

! --------- Compute geometric part of matrix at current Gauss point
            if (behaPara%lMatr) then
                matrGeom = 0.d0
                if (lMatrPred) then
                    call compRigiGeomHexaKpg(geomHexa, zeta, sigmPrep, matrGeom)
                else
                    call compRigiGeomHexaKpg(geomHexa, zeta, pk2, matrGeom)
                end if
            end if

! --------- Update tangent matrix
            if (behaPara%lMatr) then
                matrTang = matrTang+jacob*matrMate
                matrTang = matrTang+jacob*matrGeom
            end if

! --------- Update internal force at current Gauss point
            if (behaPara%lVect) then
                call btsig(elemProp%nbDof, SSH_SIZE_TENS, jacob, kineHexa%B, pk2, &
                           vect)
            end if

! --------- Save stresses
            if (behaPara%lSigm) then
                do iTens = 1, SSH_SIZE_TENS
                    sigp(iTens, kpg) = pk2(iTens)
                end do
            end if
!
        end do

! ----- Stabilization
        if (behaPara%lVect .or. behaPara%lMatr) then
            call compStabSigmHexa(geomHexa, kineHexa, Ueff, dispCurr, stabHexa, &
                                  kineHexa%epsgCurr)
        end if
        if (behaPara%lMatr) then
            call compStabMatrMateHexa(geomHexa, kineHexa, Ueff, stabHexa)
            call compStabMatrGeomHexa(geomHexa, kineHexa, stabHexa)
        end if
        if (behaPara%lVect) then
            call compStabForcHexa(kineHexa, stabHexa)
        end if

! ----- Save matrix and vector
        if (behaPara%lMatr) then
            matrTang(1:nbDofGeom, 1:nbDofGeom) = matrTang(1:nbDofGeom, 1:nbDofGeom &
                                                          )+stabHexa%matrStabMate
            matrTang(1:nbDofGeom, 1:nbDofGeom) = matrTang(1:nbDofGeom, 1:nbDofGeom &
                                                          )+stabHexa%matrStabGeom
        end if
!
        if (behaPara%lVect) then
            vect(1:SSH_NBDOFG_HEXA) = vect(1:SSH_NBDOFG_HEXA)+cellGeom%detJac0*stabHexa%forcStab
        end if
!
99      continue

! ----- Return code summary
        call codere(cod, nbIntePoint, codret)

! ----- Write matrix
        if (behaPara%lMatr) then
            if (behaPara%lMatrSyme) then
                do i = 1, SSH_NBDOF_HEXA
                    do j = 1, i
                        ij = (i-1)*i/2+j
                        matr(ij) = matrTang(i, j)
                    end do
                end do
            else
                do j = 1, SSH_NBDOF_HEXA
                    do i = 1, SSH_NBDOF_HEXA
                        ij = j+(i-1)*SSH_NBDOF_HEXA
                        matr(ij) = matrTang(i, j)
                    end do
                end do
                ASSERT(ij .le. SSH_NBDOF_HEXA*SSH_NBDOF_HEXA)
            end if
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
! --------------------------------------------------------------------------------------------------
!
! postLog
!
! Post-treatment for logarithmic quantities
!
! In  lMatrPred        : flag for prediction
! In  lMatr            : flag for compute matrix
! In  lSigm            : flag for compute stress
! In  kineHexa         : kinematic quantities (HEXA cell)
! In  tPrev            : logarithmic stress before integration of behaviour
! In  tCurr            : logarithmic stress after integration of behaviour
! In  dtde             : jacobian matrix in logarithmic space
! Out dsidep           : jacobian matrix in PK2 space
! Out pk2              : PK2 stress after integration of behaviour
!
! --------------------------------------------------------------------------------------------------
    subroutine postLog(lMatrPred, lMatr, lSigm, kineHexa, tPrev, &
                       tCurr, dtde, dsidep, pk2)
!   ------------------------------------------------------------------------------------------------
! ----- Parameters
        aster_logical, intent(in) :: lMatrPred, lMatr, lSigm
        type(SSH_KINE_HEXA), intent(in) :: kineHexa
        real(kind=8), intent(in) :: tPrev(SSH_SIZE_TENS), tCurr(SSH_SIZE_TENS)
        real(kind=8), intent(in) :: dtde(SSH_SIZE_TENS, SSH_SIZE_TENS)
        real(kind=8), intent(out) :: dsidep(SSH_SIZE_TENS, SSH_SIZE_TENS), pk2(SSH_SIZE_TENS)
! ----- Locals
        real(kind=8), parameter :: rac2 = sqrt(2.d0)
        integer(kind=8) :: i, j
        type(SSH_EPSL_HEXA) :: epsl
        real(kind=8) :: pSyme(6, 6), pSymeT(6, 6)
        real(kind=8) :: tWork(6)
        real(kind=8) :: me(3, 3, 3, 3), xi(3, 3), feta(4)
        real(kind=8) :: tl(3, 3, 3, 3), tlSyme(6, 6), trav2(6, 6)
        blas_int :: b_incx, b_incy, b_n
!   ------------------------------------------------------------------------------------------------
!
        pk2 = 0.d0
        dsidep = 0.d0
!
! - Current quantities
        if (lMatrPred) then
            epsl = kineHexa%epslPrev
            tWork = tPrev
        else
            epsl = kineHexa%epslCurr
            tWork = tCurr
        end if
!
! - Compute tensor P (symmetric)
        call deflg2(epsl%eigenVect, epsl%eigenVale, epsl%logl, pSyme, feta, &
                    xi, me)
        pSymeT = transpose(pSyme)
!
        if (lMatr) then
! ----- Compute tensor T:L
            call deflg3(epsl%eigenVect, feta, xi, me, tWork, &
                        tl)
!
! ----- Symmetric version of T:L
            call symt46(tl, tlSyme)
!
! ----- Compute dsidep
            trav2 = matmul(pSymeT, dtde)
            dsidep = matmul(trav2, pSyme)
            b_n = to_blas_int(36)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            call daxpy(b_n, 1.d0, tlSyme, b_incx, dsidep, &
                       b_incy)
!
            dsidep(4:6, 4:6) = dsidep(4:6, 4:6)/2.0d0
            dsidep(4:6, 1:3) = dsidep(4:6, 1:3)/rac2
            dsidep(1:3, 4:6) = dsidep(1:3, 4:6)/rac2
        end if
!
! - Compute PK2 stresses
        if (lSigm) then
            do i = 1, 6
                do j = 1, 6
                    pk2(i) = pk2(i)+tCurr(j)*pSyme(j, i)
                end do
            end do
            pk2(4) = pk2(4)/rac2
            pk2(5) = pk2(5)/rac2
            pk2(6) = pk2(6)/rac2
        end if
!
!   ------------------------------------------------------------------------------------------------
    end subroutine
!
end module SolidShell_NonLinear_Hexa_module
